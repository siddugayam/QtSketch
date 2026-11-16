#include "FusedSystemTopology.h"
#include <algorithm>

FusedSystemTopology::FusedSystemTopology(const std::vector<std::set<int>> &rings,
                                         const std::vector<GraphBond> &bonds)
{
    compute(rings, bonds);
}

FusedSystemTopology::FusedSystemTopology(const std::vector<std::set<int>> &rings,
                                         const Graph &g)
{
    compute(rings, g.bonds);
}

int FusedSystemTopology::ringCountForAtom(int atom) const {
    auto it = m_atomToRings.find(atom);
    if (it != m_atomToRings.end()) {
        return static_cast<int>(it->second.size());
    }
    return 0;
}

std::vector<int> FusedSystemTopology::ringsForAtom(int atom) const {
    auto it = m_atomToRings.find(atom);
    if (it != m_atomToRings.end()) {
        return it->second;
    }
    return {};
}

bool FusedSystemTopology::containsAtom(int atom) const {
    return m_atomToRings.find(atom) != m_atomToRings.end();
}

std::vector<RingAdjacency> FusedSystemTopology::adjacentRings(int ringIndex) const {
    if (ringIndex >= 0 && ringIndex < static_cast<int>(m_adjacencyList.size())) {
        return m_adjacencyList[ringIndex];
    }
    return {};
}

FusedBondType FusedSystemTopology::junctionType(int ring1, int ring2) const {
    auto it = m_pairJunctionTypes.find({ring1, ring2});
    if (it != m_pairJunctionTypes.end()) {
        return it->second;
    }
    return FusedBondType::None;
}

std::vector<int> FusedSystemTopology::sharedAtoms(int ring1, int ring2) const {
    auto it = m_pairSharedAtoms.find({ring1, ring2});
    if (it != m_pairSharedAtoms.end()) {
        return it->second;
    }
    return {};
}

void FusedSystemTopology::compute(const std::vector<std::set<int>> &rings,
                                  const std::vector<GraphBond> &bonds)
{
    m_rings = rings;
    m_atomToRings.clear();
    m_interiorAtoms.clear();
    m_peripheralAtoms.clear();
    m_perimeterCycle.clear();
    m_hasSinglePerimeterComponent = false;
    m_orthoJunctions.clear();
    m_spiroJunctions.clear();
    m_periJunctions.clear();
    m_bridgedJunctions.clear();
    m_pairJunctionTypes.clear();
    m_pairSharedAtoms.clear();
    m_adjacencyList.clear();

    int numRings = static_cast<int>(rings.size());
    if (numRings == 0) return;

    // 1. Per-atom ring membership
    for (int r = 0; r < numRings; ++r) {
        for (int a : rings[r]) {
            m_atomToRings[a].push_back(r);
        }
    }

    // 2. Partition into interior (>=3 rings) and peripheral (1 or 2 rings)
    for (const auto &pair : m_atomToRings) {
        int a = pair.first;
        size_t count = pair.second.size();
        if (count >= 3) {
            m_interiorAtoms.insert(a);
        } else if (count >= 1) {
            m_peripheralAtoms.insert(a);
        }
    }

    // 3. Peri junctions: represented by peri center atom + ring set meeting there
    for (int a : m_interiorAtoms) {
        PeriJunction pj;
        pj.atom = a;
        pj.rings = m_atomToRings[a];
        m_periJunctions.push_back(pj);
    }
    std::sort(m_periJunctions.begin(), m_periJunctions.end(),
              [](const PeriJunction &x, const PeriJunction &y) {
                  return x.atom < y.atom;
              });

    // 4. Bond lookup table
    std::set<std::pair<int, int>> bondSet;
    for (const auto &b : bonds) {
        bondSet.insert({std::min(b.u, b.v), std::max(b.u, b.v)});
    }
    auto areBonded = [&](int u, int v) -> bool {
        return bondSet.count({std::min(u, v), std::max(u, v)}) > 0;
    };

    // 5. Classify ring-pair junctions
    for (int r1 = 0; r1 < numRings; ++r1) {
        for (int r2 = r1 + 1; r2 < numRings; ++r2) {
            std::vector<int> shared;
            for (int a : rings[r1]) {
                if (rings[r2].count(a)) {
                    shared.push_back(a);
                }
            }
            if (shared.empty()) continue;

            std::sort(shared.begin(), shared.end());
            m_pairSharedAtoms[{r1, r2}] = shared;
            m_pairSharedAtoms[{r2, r1}] = shared;

            if (shared.size() == 2 && areBonded(shared[0], shared[1])) {
                // Exactly 2 atoms shared and bonded -> Ortho-fusion
                OrthoJunction oj;
                oj.ring1 = r1;
                oj.ring2 = r2;
                oj.atom1 = std::min(shared[0], shared[1]);
                oj.atom2 = std::max(shared[0], shared[1]);
                m_orthoJunctions.push_back(oj);
                m_pairJunctionTypes[{r1, r2}] = FusedBondType::Ortho;
                m_pairJunctionTypes[{r2, r1}] = FusedBondType::Ortho;
            } else if (shared.size() == 1) {
                int a = shared[0];
                if (m_atomToRings[a].size() == 2) {
                    // Shared by ONLY these 2 rings -> Spiro
                    SpiroJunction sj;
                    sj.ring1 = r1;
                    sj.ring2 = r2;
                    sj.atom = a;
                    m_spiroJunctions.push_back(sj);
                    m_pairJunctionTypes[{r1, r2}] = FusedBondType::Spiro;
                    m_pairJunctionTypes[{r2, r1}] = FusedBondType::Spiro;
                } else {
                    // Shared by 3+ rings -> Peri center takes precedence over spiro
                    m_pairJunctionTypes[{r1, r2}] = FusedBondType::Peri;
                    m_pairJunctionTypes[{r2, r1}] = FusedBondType::Peri;
                }
            } else {
                // >2 shared atoms or 2 non-adjacent atoms -> Bridged
                BridgedJunction bj;
                bj.ring1 = r1;
                bj.ring2 = r2;
                bj.sharedAtoms = shared;
                m_bridgedJunctions.push_back(bj);
                m_pairJunctionTypes[{r1, r2}] = FusedBondType::Bridged;
                m_pairJunctionTypes[{r2, r1}] = FusedBondType::Bridged;
            }
        }
    }

    // Sort pairwise junction lists canonically
    std::sort(m_orthoJunctions.begin(), m_orthoJunctions.end(),
              [](const OrthoJunction &x, const OrthoJunction &y) {
                  if (x.ring1 != y.ring1) return x.ring1 < y.ring1;
                  return x.ring2 < y.ring2;
              });
    std::sort(m_spiroJunctions.begin(), m_spiroJunctions.end(),
              [](const SpiroJunction &x, const SpiroJunction &y) {
                  if (x.ring1 != y.ring1) return x.ring1 < y.ring1;
                  return x.ring2 < y.ring2;
              });
    std::sort(m_bridgedJunctions.begin(), m_bridgedJunctions.end(),
              [](const BridgedJunction &x, const BridgedJunction &y) {
                  if (x.ring1 != y.ring1) return x.ring1 < y.ring1;
                  return x.ring2 < y.ring2;
              });

    // 6. Build ring adjacency list
    m_adjacencyList.resize(numRings);
    for (int r = 0; r < numRings; ++r) {
        for (int other = 0; other < numRings; ++other) {
            if (other == r) continue;
            auto itType = m_pairJunctionTypes.find({r, other});
            if (itType != m_pairJunctionTypes.end()) {
                RingAdjacency adj;
                adj.neighborRing = other;
                adj.type = itType->second;
                adj.sharedAtoms = m_pairSharedAtoms.at({r, other});
                m_adjacencyList[r].push_back(adj);
            }
        }
    }

    // 7. Perimeter cycle extraction (FR-5.5 / P-25.3.3)
    if (m_peripheralAtoms.empty() || !m_bridgedJunctions.empty() || !m_spiroJunctions.empty()) {
        return;
    }

    // Build perimeter boundary edge graph:
    // Bonds between peripheral atoms that belong to EXACTLY 1 ring
    std::map<int, std::vector<int>> perimAdj;
    for (const auto &b : bonds) {
        int u = b.u;
        int v = b.v;
        if (m_peripheralAtoms.count(u) && m_peripheralAtoms.count(v)) {
            int ringsContainingBoth = 0;
            for (const auto &ring : m_rings) {
                if (ring.count(u) && ring.count(v)) {
                    ringsContainingBoth++;
                }
            }
            if (ringsContainingBoth == 1) {
                perimAdj[u].push_back(v);
                perimAdj[v].push_back(u);
            }
        }
    }

    // Every peripheral atom must have degree exactly 2 in the perimeter boundary graph
    for (int a : m_peripheralAtoms) {
        if (perimAdj[a].size() != 2) {
            return;
        }
    }

    // Walk the boundary cycle starting from the peripheral atom with lowest ID
    int startAtom = *m_peripheralAtoms.begin();
    std::sort(perimAdj[startAtom].begin(), perimAdj[startAtom].end());
    int nextAtom = perimAdj[startAtom][0];

    m_perimeterCycle.push_back(startAtom);
    int curr = nextAtom;
    int prev = startAtom;

    while (curr != startAtom && m_perimeterCycle.size() <= m_peripheralAtoms.size()) {
        m_perimeterCycle.push_back(curr);
        const auto &nbrs = perimAdj[curr];
        if (nbrs.size() != 2) {
            m_perimeterCycle.clear();
            return;
        }
        int next = (nbrs[0] == prev) ? nbrs[1] : nbrs[0];
        prev = curr;
        curr = next;
    }

    if (curr == startAtom && m_perimeterCycle.size() == m_peripheralAtoms.size()) {
        m_hasSinglePerimeterComponent = true;
    } else {
        m_perimeterCycle.clear();
        m_hasSinglePerimeterComponent = false;
    }
}
