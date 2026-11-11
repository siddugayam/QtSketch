#pragma once

#include "GraphTypes.h"
#include <map>
#include <set>
#include <utility>
#include <vector>

// Classification of fusion bonds/junctions between rings (P-23.2.5, P-24, P-25.3, P-25.4)
enum class FusedBondType {
    None,
    Ortho,   // Exactly 2 atoms shared between exactly 2 rings, forming one bond
    Peri,    // One atom shared among 3 or more rings
    Spiro,   // Exactly 1 atom shared between exactly 2 rings, no bond
    Bridged  // Rings sharing >2 atoms or 2 non-adjacent atoms
};

// Represents an ortho-fusion junction between two rings sharing a single bond
struct OrthoJunction {
    int ring1 = -1;
    int ring2 = -1;
    int atom1 = -1;
    int atom2 = -1;

    bool operator==(const OrthoJunction &o) const {
        return ring1 == o.ring1 && ring2 == o.ring2 &&
               atom1 == o.atom1 && atom2 == o.atom2;
    }
};

// Represents a spiro junction between two rings sharing exactly one atom (ring count == 2)
struct SpiroJunction {
    int ring1 = -1;
    int ring2 = -1;
    int atom = -1;

    bool operator==(const SpiroJunction &o) const {
        return ring1 == o.ring1 && ring2 == o.ring2 && atom == o.atom;
    }
};

// Represents a peri junction at an atom shared by 3 or more rings
struct PeriJunction {
    int atom = -1;
    std::vector<int> rings;

    bool operator==(const PeriJunction &o) const {
        return atom == o.atom && rings == o.rings;
    }
};

// Represents a bridged junction between two rings sharing >2 atoms or non-adjacent atoms
struct BridgedJunction {
    int ring1 = -1;
    int ring2 = -1;
    std::vector<int> sharedAtoms;

    bool operator==(const BridgedJunction &o) const {
        return ring1 == o.ring1 && ring2 == o.ring2 && sharedAtoms == o.sharedAtoms;
    }
};

// Ring neighbor record for adjacency queries
struct RingAdjacency {
    int neighborRing = -1;
    FusedBondType type = FusedBondType::None;
    std::vector<int> sharedAtoms;

    bool operator==(const RingAdjacency &o) const {
        return neighborRing == o.neighborRing && type == o.type && sharedAtoms == o.sharedAtoms;
    }
};

class FusedSystemTopology {
public:
    FusedSystemTopology() = default;
    FusedSystemTopology(const std::vector<std::set<int>> &rings, const std::vector<GraphBond> &bonds);
    FusedSystemTopology(const std::vector<std::set<int>> &rings, const Graph &g);

    // Ring accessors
    const std::vector<std::set<int>> &rings() const { return m_rings; }
    int ringCount() const { return static_cast<int>(m_rings.size()); }

    // Per-atom ring membership
    const std::map<int, std::vector<int>> &atomToRings() const { return m_atomToRings; }
    int ringCountForAtom(int atom) const;
    std::vector<int> ringsForAtom(int atom) const;
    bool containsAtom(int atom) const;

    // Interior vs Peripheral atom partition
    // Interior: atoms belonging to 3 or more rings
    const std::set<int> &interiorAtoms() const { return m_interiorAtoms; }
    bool isInterior(int atom) const { return m_interiorAtoms.count(atom) > 0; }

    // Peripheral: atoms belonging to 1 or 2 rings (atoms in 0 rings are excluded)
    const std::set<int> &peripheralAtoms() const { return m_peripheralAtoms; }
    bool isPeripheral(int atom) const { return m_peripheralAtoms.count(atom) > 0; }

    // Classified junctions
    const std::vector<OrthoJunction> &orthoJunctions() const { return m_orthoJunctions; }
    const std::vector<SpiroJunction> &spiroJunctions() const { return m_spiroJunctions; }
    const std::vector<PeriJunction> &periJunctions() const { return m_periJunctions; }
    const std::vector<BridgedJunction> &bridgedJunctions() const { return m_bridgedJunctions; }

    // Convenience accessors
    const std::set<int> &periCenters() const { return m_interiorAtoms; }

    // Ring adjacency
    std::vector<RingAdjacency> adjacentRings(int ringIndex) const;
    FusedBondType junctionType(int ring1, int ring2) const;
    std::vector<int> sharedAtoms(int ring1, int ring2) const;

private:
    void compute(const std::vector<std::set<int>> &rings, const std::vector<GraphBond> &bonds);

    std::vector<std::set<int>> m_rings;
    std::map<int, std::vector<int>> m_atomToRings;
    std::set<int> m_interiorAtoms;
    std::set<int> m_peripheralAtoms;

    std::vector<OrthoJunction> m_orthoJunctions;
    std::vector<SpiroJunction> m_spiroJunctions;
    std::vector<PeriJunction> m_periJunctions;
    std::vector<BridgedJunction> m_bridgedJunctions;

    std::map<std::pair<int, int>, FusedBondType> m_pairJunctionTypes;
    std::map<std::pair<int, int>, std::vector<int>> m_pairSharedAtoms;
    std::vector<std::vector<RingAdjacency>> m_adjacencyList;
};
