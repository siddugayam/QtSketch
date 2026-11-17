// tools/iupac/AtomOrderRandomizer.cpp
#include "AtomOrderRandomizer.h"
#include "indigo.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace {

struct BondSpec {
    int a;
    int b;
    int order;
};

} // namespace

bool randomizeAtomOrder(const std::string& smiles, int seed, std::string& outSmiles) {
    int mol = indigoLoadMoleculeFromString(smiles.c_str());
    if (mol < 0) return false;

    // Skip stereo-sensitive molecules: this reordering does not re-derive
    // stereo descriptors (tetrahedral parity, cis/trans) for the new atom
    // order, so a naive permutation could silently produce a molecule
    // with wrong stereochemistry. Bail out rather than risk that.
    if (indigoCountStereocenters(mol) > 0) {
        indigoFree(mol);
        return false;
    }
    bool hasBondStereo = false;
    {
        int it = indigoIterateBonds(mol);
        int b;
        while ((b = indigoNext(it)) != 0) {
            if (indigoBondStereo(b) != 0) hasBondStereo = true;
            indigoFree(b);
        }
        indigoFree(it);
    }
    if (hasBondStereo) {
        indigoFree(mol);
        return false;
    }

    int n = indigoCountAtoms(mol);
    if (n <= 0) {
        indigoFree(mol);
        return false;
    }

    std::vector<std::string> symbols(n);
    std::vector<int> charges(n, 0);
    std::vector<int> isotopes(n, 0);
    {
        int it = indigoIterateAtoms(mol);
        int a;
        while ((a = indigoNext(it)) != 0) {
            int idx = indigoIndex(a);
            if (idx >= 0 && idx < n) {
                const char* sym = indigoSymbol(a);
                symbols[idx] = sym ? sym : "C";
                int charge = 0;
                indigoGetCharge(a, &charge);
                charges[idx] = charge;
                isotopes[idx] = indigoIsotope(a);
            }
            indigoFree(a);
        }
        indigoFree(it);
    }

    // Capture the canonical form before any mutation -- this is what the
    // reordered molecule must match regardless of aromatic/Kekule
    // representation.
    const char* origCanonC = indigoCanonicalSmiles(mol);
    std::string origCanon = origCanonC ? origCanonC : "";
    if (origCanon.empty()) {
        indigoFree(mol);
        return false;
    }

    // Dearomatize so indigoBondOrder() returns concrete single/double
    // bond orders (1/2/3) that indigoAddBond() can reconstruct -- an
    // aromatic bond's order is not directly re-addable when building a
    // molecule from scratch atom-by-atom.
    indigoDearomatize(mol);

    std::vector<BondSpec> bonds;
    {
        int it = indigoIterateBonds(mol);
        int b;
        while ((b = indigoNext(it)) != 0) {
            int src = indigoSource(b);
            int dst = indigoDestination(b);
            bonds.push_back({indigoIndex(src), indigoIndex(dst), indigoBondOrder(b)});
            indigoFree(src);
            indigoFree(dst);
            indigoFree(b);
        }
        indigoFree(it);
    }

    indigoFree(mol);

    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937 rng(static_cast<unsigned int>(seed));
    std::shuffle(perm.begin(), perm.end(), rng);

    int newMol = indigoCreateMolecule();
    if (newMol < 0) return false;

    std::vector<int> oldIndexToNewHandle(n, -1);
    bool atomFailure = false;
    for (int newPos = 0; newPos < n && !atomFailure; newPos++) {
        int oldIdx = perm[newPos];
        int newAtom = indigoAddAtom(newMol, symbols[oldIdx].c_str());
        if (newAtom < 0) {
            atomFailure = true;
            break;
        }
        if (charges[oldIdx] != 0) indigoSetCharge(newAtom, charges[oldIdx]);
        if (isotopes[oldIdx] > 0) indigoSetIsotope(newAtom, isotopes[oldIdx]);
        oldIndexToNewHandle[oldIdx] = newAtom;
    }
    if (atomFailure) {
        indigoFree(newMol);
        return false;
    }

    for (const auto& bs : bonds) {
        if (bs.a < 0 || bs.a >= n || bs.b < 0 || bs.b >= n) continue;
        if (oldIndexToNewHandle[bs.a] < 0 || oldIndexToNewHandle[bs.b] < 0) continue;
        indigoAddBond(oldIndexToNewHandle[bs.a], oldIndexToNewHandle[bs.b], bs.order);
    }

    // Restore aromatic perception (we dearomatized the source molecule
    // to get re-addable concrete bond orders) so the emitted SMILES is
    // in the same representation style as the input.
    indigoAromatize(newMol);

    const char* reorderedC = indigoSmiles(newMol);
    std::string reordered = reorderedC ? reorderedC : "";
    indigoFree(newMol);

    // Indigo's SMILES writer emits "[Te]" (capitalized) even for aromatic
    // tellurium in rings, but Indigo's SMILES parser treats uppercase "[Te]"
    // as aliphatic, failing aromaticity perception and adding spurious
    // implicit hydrogens (turning tellurophene into dihydrotellurophene).
    // If [Te] was emitted in an aromatic structure, lowercase to [te].
    if (reordered.find("[Te]") != std::string::npos) {
        std::string candidate = reordered;
        size_t pos = 0;
        while ((pos = candidate.find("[Te]", pos)) != std::string::npos) {
            candidate.replace(pos, 4, "[te]");
            pos += 4;
        }
        int testMol = indigoLoadMoleculeFromString(candidate.c_str());
        if (testMol >= 0) {
            const char* testCanonC = indigoCanonicalSmiles(testMol);
            if (testCanonC && std::string(testCanonC) == origCanon) {
                reordered = candidate;
            }
            indigoFree(testMol);
        }
    }

    // Self-check: must still be the same molecule when re-parsed from the emitted SMILES.
    if (reordered.empty()) {
        return false;
    }
    int reparsedCheck = indigoLoadMoleculeFromString(reordered.c_str());
    if (reparsedCheck < 0) {
        return false;
    }
    const char* reparsedCanonC = indigoCanonicalSmiles(reparsedCheck);
    std::string reparsedCanon = reparsedCanonC ? reparsedCanonC : "";
    indigoFree(reparsedCheck);
    if (reparsedCanon != origCanon) {
        return false;
    }

    // Self-check: for molecules with more than 2 atoms, the reordered raw
    // SMILES must actually differ from a straight re-parse of the input --
    // otherwise this has silently degraded back into a no-op and must
    // report as "could not verify a real reordering", not a false pass.
    if (n > 2) {
        int reparsed = indigoLoadMoleculeFromString(smiles.c_str());
        std::string reparsedRaw;
        if (reparsed >= 0) {
            const char* rs = indigoSmiles(reparsed);
            reparsedRaw = rs ? rs : "";
            indigoFree(reparsed);
        }
        if (reordered == reparsedRaw) {
            // Symmetric molecules with no distinguishing substituents (e.g.
            // unsubstituted benzene) can legitimately produce the same
            // SMILES string under every permutation -- the string encodes
            // graph structure, not atom identity. Treat as "cannot verify
            // a reordering" rather than a false pass.
            return false;
        }
    }

    outSmiles = reordered;
    return true;
}
