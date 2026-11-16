#include "BareFusedHydrocarbonNamer.h"
#include "FusedSystemTopology.h"
#include "indigo.h"

#include <set>
#include <map>
#include <vector>
#include <algorithm>

BareFusedHydrocarbonResult tryNameBareFusedHydrocarbon(
    int mol, const Graph &g, const FusedRingSystemInput &input) {
    // Gate 1: Exactly 3 or 4 rings, no detector error
    if ((input.ringSizes.size() != 3 && input.ringSizes.size() != 4) || !input.error.empty()) {
        return {BareFusedDisposition::NotMatched, {}, ""};
    }

    // Gate 2: Cyclic ring graph (peri-fused, not tree)
    if (input.closingFusions.empty()) {
        return {BareFusedDisposition::NotMatched, {}, ""};
    }

    // Gate 3: Exactly 1 or 2 interior atoms
    if (input.interiorAtoms.size() != 1 && input.interiorAtoms.size() != 2) {
        return {BareFusedDisposition::NotMatched, {}, ""};
    }

    // Gate 4: All atoms carbon, net charge 0, no radicals, no isotope labels
    int atomIter = indigoIterateAtoms(mol);
    if (atomIter < 0) return {BareFusedDisposition::NotMatched, {}, ""};
    while (int a = indigoNext(atomIter)) {
        int z = indigoAtomicNumber(a);
        int charge = 0;
        indigoGetCharge(a, &charge);
        int rad = 0;
        indigoGetRadicalElectrons(a, &rad);
        int iso = indigoIsotope(a);
        indigoFree(a);
        if (charge != 0 || rad > 0 || iso != 0 || (z != 6 && z != 1)) {
            indigoFree(atomIter);
            return {BareFusedDisposition::NotMatched, {}, ""};
        }
    }
    indigoFree(atomIter);

    for (const auto &node : g.nodes) {
        if (node.atomicNumber != 6) {
            return {BareFusedDisposition::NotMatched, {}, ""};
        }
    }

    // Gate 5: g node set == union(ringAtoms); every bond internal (bare)
    std::map<int, int> indigoToGraph;
    for (const auto &node : g.nodes) {
        indigoToGraph[node.indigoIdx] = node.id;
    }

    std::set<int> ringGraphNodes;
    for (const auto &ring : input.ringAtoms) {
        for (int a : ring) {
            auto it = indigoToGraph.find(a);
            if (it == indigoToGraph.end()) {
                return {BareFusedDisposition::NotMatched, {}, ""};
            }
            ringGraphNodes.insert(it->second);
        }
    }

    if (ringGraphNodes.size() != g.nodes.size()) {
        return {BareFusedDisposition::NotMatched, {}, ""};
    }

    // Convert input to graph-index coordinate space for peripheral numbering
    FusedRingSystemInput inputGraph = input;
    for (auto &ring : inputGraph.ringAtoms) {
        std::set<int> converted;
        for (int idx : ring) {
            converted.insert(indigoToGraph.at(idx));
        }
        ring = converted;
    }
    for (auto &cyc : inputGraph.ringDirectedCycles) {
        for (int &idx : cyc) {
            idx = indigoToGraph.at(idx);
        }
    }
    std::set<int> interiorGraphNodes;
    for (int idx : input.interiorAtoms) {
        auto it = indigoToGraph.find(idx);
        if (it != indigoToGraph.end()) {
            interiorGraphNodes.insert(it->second);
        }
    }
    inputGraph.interiorAtoms = interiorGraphNodes;
    if (interiorGraphNodes.size() != 1 && interiorGraphNodes.size() != 2) {
        return {BareFusedDisposition::NotMatched, {}, ""};
    }

    // Identify sp3 atoms in g (atoms whose bonds are all single)
    std::set<int> sp3Atoms;
    for (const auto &node : g.nodes) {
        bool allSingle = true;
        for (int order : node.bondOrders) {
            if (order != 1) {
                allSingle = false;
                break;
            }
        }
        if (allSingle) {
            sp3Atoms.insert(node.id);
        }
    }

    // Canonical descriptor: perimeter length, interior count, total atoms (SSSR-independent)
    FusedSystemTopology topo(inputGraph.ringAtoms, g.bonds);
    auto perim = topo.perimeterCycle();
    int perimLen = static_cast<int>(perim.size());
    int interiorCount = static_cast<int>(interiorGraphNodes.size());
    int totalAtoms = static_cast<int>(g.nodes.size());

    // Acenaphthylene: 11-atom perimeter, 1 interior atom, 12 total carbons
    if (perimLen == 11 && interiorCount == 1 && totalAtoms == 12) {
        if (sp3Atoms.empty()) {
            return {BareFusedDisposition::Named, {true, "acenaphthylene", "", NamingStatus::SystematicIupacName}, ""};
        } else if (sp3Atoms.size() == 2) {
            int u = *sp3Atoms.begin();
            int v = *std::next(sp3Atoms.begin());
            bool adjacent = false;
            for (int nei : g.nodes[u].neighbors) {
                if (nei == v) {
                    adjacent = true;
                    break;
                }
            }
            if (adjacent) {
                return {BareFusedDisposition::Rejected, {}, "1,2-dihydroacenaphthylene is out of Stage C2 scope (bare parents only)."};
            } else {
                return {BareFusedDisposition::Rejected, {}, "unrecognized acenaphthylene hydrogenation state."};
            }
        } else {
            return {BareFusedDisposition::Rejected, {}, "unrecognized acenaphthylene hydrogenation state."};
        }
    }

    // Phenalene: 12-atom perimeter, 1 interior atom, 13 total carbons
    if (perimLen == 12 && interiorCount == 1 && totalAtoms == 13) {
        if (sp3Atoms.empty()) {
            return {BareFusedDisposition::Rejected, {}, "requires indicated-hydrogen tautomer"};
        }
        if (sp3Atoms.size() > 1) {
            return {BareFusedDisposition::Rejected, {}, "over-reduced phenalene skeleton"};
        }

        // Exactly one sp3 atom
        int sp3 = *sp3Atoms.begin();
        int interiorNode = *interiorGraphNodes.begin();
        if (sp3 == interiorNode) {
            return {BareFusedDisposition::Named, {true, "9bH-phenalene", "", NamingStatus::SystematicIupacName}, ""};
        }

        std::map<int, int> atomRingCount;
        for (const auto &ring : inputGraph.ringAtoms) {
            for (int a : ring) atomRingCount[a]++;
        }

        auto isCornerAtom = [&](int a) -> bool {
            if (atomRingCount[a] != 1) return false;
            for (int nei : g.nodes[a].neighbors) {
                if (atomRingCount[nei] >= 2) return true;
            }
            return false;
        };

        if (isCornerAtom(sp3)) {
            return {BareFusedDisposition::Named, {true, "1H-phenalene", "", NamingStatus::SystematicIupacName}, ""};
        }

        return {BareFusedDisposition::Rejected, {},
            "non-Kekule diradical phenalene isomer (sp3 at an unstarred position) is not a closed-shell parent hydride, not supported."};
    }

    // Pyrene: 14-atom perimeter, 2 interior atoms, 16 total carbons
    if (perimLen == 14 && interiorCount == 2 && totalAtoms == 16) {
        if (sp3Atoms.empty()) {
            return {BareFusedDisposition::Named, {true, "pyrene", "", NamingStatus::SystematicIupacName}, ""};
        }
        return {BareFusedDisposition::Rejected, {}, "hydro-pyrene derivatives are out of Stage C2 scope (bare parents only)."};
    }

    return {BareFusedDisposition::NotMatched, {}, ""};
}
