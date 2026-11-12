#include "FusedRingDirectionDetector.h"
#include "FusedSystemTopology.h"
#include "indigo.h"
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <queue>
#include <cmath>

namespace {
constexpr double kPI = 3.14159265358979323846;
}

static std::vector<int> buildDirectedCycle(const std::set<int>& rNodes, int u, int v, const std::map<int, std::vector<int>>& neighbors) {
    std::vector<int> cycle;
    cycle.push_back(u);
    cycle.push_back(v);
    int prev = u;
    int curr = v;
    while (cycle.size() < rNodes.size()) {
        int next = -1;
        auto it = neighbors.find(curr);
        if (it != neighbors.end()) {
            for (int nei : it->second) {
                if (rNodes.count(nei) && nei != prev) {
                    next = nei;
                    break;
                }
            }
        }
        if (next == -1) break; // Should not happen for a valid simple cycle
        cycle.push_back(next);
        prev = curr;
        curr = next;
    }
    return cycle;
}

FusedRingSystemInput detectFusedRingDirections(int mol) {
    FusedRingSystemInput input;
    
    // 1. Build atom adjacency list
    std::map<int, std::vector<int>> adj;
    int atomIter = indigoIterateAtoms(mol);
    if (atomIter > 0) {
        while (int atom = indigoNext(atomIter)) {
            int idx = indigoIndex(atom);
            int neiIter = indigoIterateNeighbors(atom);
            if (neiIter > 0) {
                while (int nei = indigoNext(neiIter)) {
                    adj[idx].push_back(indigoIndex(nei));
                    indigoFree(nei);
                }
                indigoFree(neiIter);
            }
            indigoFree(atom);
        }
        indigoFree(atomIter);
    }
    
    // 2. Extract SSSR rings
    std::vector<std::set<int>> rings;
    int sssrIter = indigoIterateSSSR(mol);
    if (sssrIter > 0) {
        while (int subMol = indigoNext(sssrIter)) {
            std::set<int> rNodes;
            int rAtomIter = indigoIterateAtoms(subMol);
            if (rAtomIter > 0) {
                while (int a = indigoNext(rAtomIter)) {
                    rNodes.insert(indigoIndex(a));
                    indigoFree(a);
                }
                indigoFree(rAtomIter);
            }
            rings.push_back(rNodes);
            indigoFree(subMol);
        }
        indigoFree(sssrIter);
    }
    
    int numRings = static_cast<int>(rings.size());
    if (numRings < 1) return input;
    
    // Ring sizes 3-10 supported (matches this file's own monocyclic
    // Hantzsch-Widman range, RingType::GENERAL_HETEROCYCLE); outside that,
    // this detector's own downstream geometry (FusedRingOrientation) has no
    // defined regular-polygon behavior, so reject rather than guess.
    for (int i = 0; i < numRings; ++i) {
        if (rings[i].size() < 3 || rings[i].size() > 10) return input;
    }

    // 3. Topology perception via FusedSystemTopology
    std::vector<GraphBond> bonds;
    for (const auto &p : adj) {
        int u = p.first;
        for (int v : p.second) {
            if (u < v) {
                bonds.push_back({u, v, 1});
            }
        }
    }
    FusedSystemTopology topo(rings, bonds);

    if (!topo.bridgedJunctions().empty()) {
        input.error = "Bridged ring junctions are not supported.";
        return input;
    }
    if (!topo.spiroJunctions().empty()) {
        input.error = "Spiro ring junctions are not supported.";
        return input;
    }

    // 4. Find shared bonds and ring adjacency
    std::vector<std::vector<int>> ringAdj(numRings);
    std::map<std::pair<int, int>, std::pair<int, int>> sharedAtoms;
    
    for (int i = 0; i < numRings; ++i) {
        for (int j = i + 1; j < numRings; ++j) {
            std::vector<int> common;
            std::set_intersection(rings[i].begin(), rings[i].end(),
                                  rings[j].begin(), rings[j].end(),
                                  std::back_inserter(common));
            if (common.size() == 2) {
                ringAdj[i].push_back(j);
                ringAdj[j].push_back(i);
                sharedAtoms[{i, j}] = {common[0], common[1]};
                sharedAtoms[{j, i}] = {common[0], common[1]};
            } else if (common.size() > 2) {
                // Bridged ring pair sharing > 2 atoms
                input.error = "Bridged ring junctions sharing > 2 atoms are not supported.";
                return input; 
            }
        }
    }
    
    input.ringAtoms = rings;
    input.ringSizes.resize(numRings);
    for (int i = 0; i < numRings; ++i) {
        input.ringSizes[i] = static_cast<int>(rings[i].size());
    }
    input.ringDirectedCycles.resize(numRings);
    if (numRings == 1) {
        // Single ring, no fusions: still expose a directed cycle (starting
        // atom/direction are both arbitrary here -- no fusion bond to anchor
        // against -- but must be non-empty and real so callers like
        // walkPeripheryClockwise can place this ring's atoms consistently).
        int anyAtom = *rings[0].begin();
        int anyNeighbor = -1;
        auto it = adj.find(anyAtom);
        if (it != adj.end()) {
            for (int n : it->second) {
                if (rings[0].count(n)) { anyNeighbor = n; break; }
            }
        }
        if (anyNeighbor != -1) {
            input.ringDirectedCycles[0] = buildDirectedCycle(rings[0], anyAtom, anyNeighbor, adj);
        }
        return input;
    }

    // Detect when ring-adjacency graph contains a cycle (more fusion edges than spanning tree
    // or presence of peri-fusion junctions).
    bool hasRingCycle = (topo.orthoJunctions().size() >= static_cast<size_t>(numRings)) || (!topo.periJunctions().empty());

    if (!hasRingCycle) {
        std::vector<bool> visited(numRings, false);
        
        // Pick any ring with degree <= 1 as a traversal root if one exists (a leaf of the
        // tree), otherwise fall back to ring 0.
        int rootRing = 0;
        for (int i = 0; i < numRings; ++i) {
            if (ringAdj[i].size() <= 1) { rootRing = i; break; }
        }
        
        visited[rootRing] = true;
        if (ringAdj[rootRing].empty()) return input;
        int firstNeighbor = ringAdj[rootRing][0];
        std::pair<int,int> rootBond = sharedAtoms[{rootRing, firstNeighbor}];

        struct QueueEntry { int ring; int u; int v; };
        std::vector<QueueEntry> queue;
        queue.push_back({rootRing, rootBond.first, rootBond.second});
        size_t qi = 0;

        while (qi < queue.size()) {
            QueueEntry entry = queue[qi++];
            int r = entry.ring;
            int u = entry.u;
            int v = entry.v;

            std::vector<int> cycle = buildDirectedCycle(rings[r], u, v, adj);
            input.ringDirectedCycles[r] = cycle;

            for (int nextR : ringAdj[r]) {
                if (visited[nextR]) continue;
                visited[nextR] = true;

                std::pair<int,int> b_out = sharedAtoms[{r, nextR}];
                int out_u = -1, out_v = -1;
                int k = -1;
                for (size_t c = 0; c < cycle.size(); ++c) {
                    int c1 = cycle[c];
                    int c2 = cycle[(c + 1) % cycle.size()];
                    if ((c1 == b_out.first && c2 == b_out.second) || (c1 == b_out.second && c2 == b_out.first)) {
                        k = static_cast<int>(c);
                        out_u = c1;
                        out_v = c2;
                        break;
                    }
                }
                if (k == -1) {
                    input.ringSizes.clear();
                    return input;
                }

                int outgoingDir = k;
                input.fusions.push_back({r, nextR, outgoingDir});
                queue.push_back({nextR, out_v, out_u});
            }
        }
        
        return input;
    }

    // ---------------- Cyclic Ring-Graph Placement Branch ----------------
    std::vector<bool> visited(numRings, false);
    std::vector<int> parent(numRings, -1);
    int rootRing = 0;
    size_t maxDeg = 0;
    for (int i = 0; i < numRings; ++i) {
        if (ringAdj[i].size() > maxDeg ||
            (ringAdj[i].size() == maxDeg && rings[i].size() > rings[rootRing].size())) {
            maxDeg = ringAdj[i].size();
            rootRing = i;
        }
    }
    visited[rootRing] = true;
    if (ringAdj[rootRing].empty()) return input;
    int firstNeighbor = ringAdj[rootRing][0];
    std::pair<int,int> rootBond = sharedAtoms[{rootRing, firstNeighbor}];

    struct ClosingEdge {
        int r1, r2;
        int u, v;
    };
    std::vector<ClosingEdge> closingEdges;

    struct QueueEntry { int ring; int u; int v; };
    std::vector<QueueEntry> queue;
    queue.push_back({rootRing, rootBond.first, rootBond.second});
    size_t qi = 0;

    while (qi < queue.size()) {
        QueueEntry entry = queue[qi++];
        int r = entry.ring;
        int u = entry.u;
        int v = entry.v;

        std::vector<int> cycle = buildDirectedCycle(rings[r], u, v, adj);
        input.ringDirectedCycles[r] = cycle;

        for (int nextR : ringAdj[r]) {
            if (visited[nextR]) {
                if (nextR != parent[r] && r < nextR) {
                    std::pair<int,int> b_out = sharedAtoms[{r, nextR}];
                    closingEdges.push_back({r, nextR, b_out.first, b_out.second});
                }
                continue;
            }
            visited[nextR] = true;
            parent[nextR] = r;

            std::pair<int,int> b_out = sharedAtoms[{r, nextR}];
            int out_u = -1, out_v = -1;
            int k = -1;
            for (size_t c = 0; c < cycle.size(); ++c) {
                int c1 = cycle[c];
                int c2 = cycle[(c + 1) % cycle.size()];
                if ((c1 == b_out.first && c2 == b_out.second) || (c1 == b_out.second && c2 == b_out.first)) {
                    k = static_cast<int>(c);
                    out_u = c1;
                    out_v = c2;
                    break;
                }
            }
            if (k == -1) {
                input.ringSizes.clear();
                input.error = "Failed to locate fusion bond in ring cycle.";
                return input;
            }

            int outgoingDir = k;
            input.fusions.push_back({r, nextR, outgoingDir});
            queue.push_back({nextR, out_v, out_u});
        }
    }

    if (queue.size() != static_cast<size_t>(numRings)) {
        input.ringSizes.clear();
        input.error = "Disconnected cyclic ring system.";
        return input;
    }

    // Explicit geometric consistency verification of closing edge(s):
    // Place all rings in 2D using the spanning tree fusions in input.fusions.
    auto apothemFn = [](int ringSize) -> double {
        return 1.0 / (2.0 * std::tan(kPI / ringSize));
    };
    auto edgeMidpointAngleFn = [](double ringRot, int ringSize, int bondIdx) -> double {
        double vertexSpacing = 2.0 * kPI / ringSize;
        return ringRot + (bondIdx + 0.5) * vertexSpacing;
    };

    struct PlacedGeom {
        double x = 0.0, y = 0.0;
        double rotation = 0.0;
        bool placed = false;
    };
    std::vector<PlacedGeom> placed(numRings);
    placed[rootRing] = {0.0, 0.0, 0.0, true};

    struct FusionAdj { int neighbor; int bondIndex; bool isOutgoing; };
    std::vector<std::vector<FusionAdj>> treeAdj(numRings);
    for (const auto &f : input.fusions) {
        treeAdj[f.ring1].push_back({f.ring2, f.bondIndex.value, true});
        treeAdj[f.ring2].push_back({f.ring1, f.bondIndex.value, false});
    }

    std::queue<int> gq;
    gq.push(rootRing);
    while (!gq.empty()) {
        int u = gq.front();
        gq.pop();
        for (const auto &nbr : treeAdj[u]) {
            int v = nbr.neighbor;
            if (placed[v].placed) continue;
            double edgeAngle;
            if (nbr.isOutgoing) {
                edgeAngle = edgeMidpointAngleFn(placed[u].rotation, input.ringSizes[u], nbr.bondIndex);
            } else {
                edgeAngle = edgeMidpointAngleFn(placed[u].rotation, input.ringSizes[u], nbr.bondIndex) + kPI;
            }
            double dist = apothemFn(input.ringSizes[u]) + apothemFn(input.ringSizes[v]);
            double vx = placed[u].x + dist * std::cos(edgeAngle);
            double vy = placed[u].y + dist * std::sin(edgeAngle);
            double vRot = edgeAngle + kPI - 0.5 * (2.0 * kPI / input.ringSizes[v]);
            placed[v] = {vx, vy, vRot, true};
            gq.push(v);
        }
    }

    auto getAtomCoord = [&](int r, int atom) -> std::pair<double, double> {
        const auto &cyc = input.ringDirectedCycles[r];
        int idx = -1;
        for (size_t i = 0; i < cyc.size(); ++i) {
            if (cyc[i] == atom) { idx = static_cast<int>(i); break; }
        }
        if (idx == -1) return {1e9, 1e9};
        double radius = 1.0 / (2.0 * std::sin(kPI / input.ringSizes[r]));
        double angle = placed[r].rotation + idx * (2.0 * kPI / input.ringSizes[r]);
        return {placed[r].x + radius * std::cos(angle), placed[r].y + radius * std::sin(angle)};
    };

    const double CLOSING_TOLERANCE = 0.25;
    for (const auto &ce : closingEdges) {
        auto p1_u = getAtomCoord(ce.r1, ce.u);
        auto p2_u = getAtomCoord(ce.r2, ce.u);
        auto p1_v = getAtomCoord(ce.r1, ce.v);
        auto p2_v = getAtomCoord(ce.r2, ce.v);
        double d_u = std::hypot(p1_u.first - p2_u.first, p1_u.second - p2_u.second);
        double d_v = std::hypot(p1_v.first - p2_v.first, p1_v.second - p2_v.second);
        if (d_u > CLOSING_TOLERANCE || d_v > CLOSING_TOLERANCE) {
            input.ringSizes.clear();
            input.fusions.clear();
            input.closingFusions.clear();
            input.ringDirectedCycles.clear();
            input.error = "Geometric contradiction: Impossible cyclic fusion closure (mismatched ring sizes or distorted geometry)";
            return input;
        }
    }

    for (const auto &ce : closingEdges) {
        const auto &cyc = input.ringDirectedCycles[ce.r1];
        int k = -1;
        for (size_t c = 0; c < cyc.size(); ++c) {
            int c1 = cyc[c];
            int c2 = cyc[(c + 1) % cyc.size()];
            if ((c1 == ce.u && c2 == ce.v) || (c1 == ce.v && c2 == ce.u)) {
                k = static_cast<int>(c);
                break;
            }
        }
        if (k != -1) {
            input.closingFusions.push_back({ce.r1, ce.r2, k});
        }
    }

    return input;
}
