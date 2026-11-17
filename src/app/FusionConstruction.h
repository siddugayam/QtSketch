#pragma once

#include "GraphTypes.h"
#include "FusedRingOrientation.h"
#include <vector>
#include <string>
#include <set>
#include <map>
#include <QString>
#include <QStringList>

struct FusionComponent {
    int id{0};
    QString name;               // e.g. "phenanthrene", "benzo", "pyridine"
    QString prefix;             // e.g. "benzo", "furo", "thieno", "naphtho"
    int seniorityRank{0};       // Seniority rank (lower number = higher seniority)
    int ringCount{0};
    std::set<int> ringIndices;  // indices into the system's SSSR rings
    std::set<int> atomIndices;  // atom IDs in the graph
    int order{0};               // 0 = parent, 1 = first-order attached, 2 = second-order, etc.
    bool isMergedBenzoheterocycle{false};
};

struct FusionJunction {
    int comp1{-1};
    int comp2{-1};
    std::vector<int> sharedAtoms; // atoms in g shared between the two components
    std::vector<std::pair<int, int>> sharedBonds;
};

struct ComponentTree {
    std::vector<FusionComponent> components;
    std::vector<FusionJunction> junctions;
    int parentIndex{-1};
    bool valid{false};
    QString error;
};

struct ParentScore {
    bool heteroParent{false};
    int parentRingCount{0};
    int parentSeniority{0}; // lower number = more senior (FR-2.3)
    bool noFr7Fallback{true}; // FR-3.3.1 (a)
    int higherOrderCount{0}; // fewer is better (FR-3.3.1 b)
    int firstOrderCount{0};  // more is better (FR-3.3.1 c)
    int maxIdenticalAttached{0}; // more is better (FR-3.3.1 d)
    int interparentSeniority{9999}; // lower number = more senior interparent (FR-3.3.1 e)
    QString parentName;
    std::vector<int> sortedParentRings;

    bool operator<(const ParentScore &other) const {
        if (heteroParent != other.heteroParent) return heteroParent < other.heteroParent;
        if (parentRingCount != other.parentRingCount) return parentRingCount < other.parentRingCount;
        if (parentSeniority != other.parentSeniority) return parentSeniority > other.parentSeniority; // lower rank is better (FR-2.3)
        if (noFr7Fallback != other.noFr7Fallback) return noFr7Fallback < other.noFr7Fallback;
        if (higherOrderCount != other.higherOrderCount) return higherOrderCount > other.higherOrderCount; // fewer is better
        if (firstOrderCount != other.firstOrderCount) return firstOrderCount < other.firstOrderCount; // more is better
        if (maxIdenticalAttached != other.maxIdenticalAttached) return maxIdenticalAttached < other.maxIdenticalAttached;
        if (interparentSeniority != other.interparentSeniority) return interparentSeniority > other.interparentSeniority; // lower rank is better (FR-3.3.1 e)
        if (parentName != other.parentName) return parentName > other.parentName;
        return sortedParentRings > other.sortedParentRings;
    }
};

// Score a component tree for parent-selection
ParentScore scoreTree(const ComponentTree &tree, const Graph &g);

// Check if a component is fused to two mutually-fused components (FR-3.2 limitation)
// Under FR-3.2, standard fusion nomenclature provides no means of describing the fusion
// of one component to two other components which are themselves fused together.
bool hitsFr32Limitation(const ComponentTree &tree);

// Build component tree from ring system and molecular graph
ComponentTree decomposeComponentTree(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol = -1,
    bool allowFr32 = true);

// Order component citation according to FR-3.6:
// Higher-order components cited immediately before the host to which they are attached.
// Components at the same level cited in alphabetical order.
// Parent component cited last.
QStringList orderCitation(const ComponentTree &tree);

// Apply FR-3.5 benzoheterocycle merging when it doesn't break multiplicative prefixes
void applyBenzoheterocycleMerging(ComponentTree &tree, const Graph &g, int mol = -1);

// Encodes the fusion descriptor (e.g. "[1,2-a:2',1'-c]", "[3,2-f]", "[1,2-c:7,8-c']")
// applying parent side lettering (FR-4.1), attached numeric locants (FR-4.2),
// order primes (FR-4.3), and the 7-step locant preference ordering (FR-4.5).
QString encodeFusionDescriptor(const ComponentTree &tree, const Graph &g, int mol = -1);

// Constructs the complete IUPAC fusion name from the component tree,
// applying citation order (FR-3.6), descriptor encoding (FR-4.1-4.5),
// locant omission (FR-4.6), and identical component handling / punctuation hierarchy (FR-4.9).
QString constructFusionName(const ComponentTree &tree, const Graph &g, int mol = -1);

// Phase 5 Task 12: Alternative 1 (P-25.5.1 / FR-7.1)
// Skeletal replacement ('a' nomenclature) on an eligible fused hydrocarbon parent.
// Applies when hitsFr32Limitation(tree) is true and the underlying ring skeleton matches
// an eligible fused hydrocarbon (e.g. cyclopenta[cd]pentalene, phenalene).
QString nameBySkeletalReplacementFallback(
    const ComponentTree &tree,
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    int mol = -1);

// Phase 5 Task 13: Alternative 2 (P-25.5.2 / FR-7.1.1)
// Less preferred parent component fallback when standard senior parent hits FR-3.2 limitation.
QString nameByLessPreferredParentFallback(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol = -1);

// Phase 5 Task 13: Alternative 3 (P-25.5.3 / FR-7.2 / FR-8)
// Bridging nomenclature on eligible fused parent.
QString nameByBridgedFusedFallback(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol = -1);



