#pragma once
#include <optional>
#include <vector>
#include <set>
#include <map>
#include <QString>
#include <QStringList>
#include "NamingTypes.h"
#include "RetainedRingTypes.h"

// Result of a successful skeleton-only (bond-order-agnostic) match.
struct SkeletonMatchResult {
    std::map<int,int> atomMap;
    QStringList hydroLocants;
    QString indicatedHydrogenLocant;
};

// Returns full retained ring registry
const std::vector<RetainedRingSystem>& fullRegistry();

// Returns true and fills outAtomMap if sys matches the induced subgraph of ringAtoms
bool matchesRetainedSystem(int mol, const std::set<int> &ringAtoms,
                           const RetainedRingSystem &sys,
                           std::map<int,int> &outAtomMap);

// Fallback for matchesRetainedSystem: same topology + heteroatom identity, bond order ignored
bool skeletonMatchRetainedSystem(int mol, const std::set<int> &ringAtoms,
                                 const RetainedRingSystem &sys,
                                 SkeletonMatchResult &outResult);

// Computes all symmetry automorphisms of a reference molecule
std::vector<std::vector<int>> computeAutomorphisms(int refMol, bool preserveBondOrders = false);

// Tries to decompose the molecule's fused-ring system into retained
// polycyclic components (base + attached prefixes) per Blue Book
// FR-2.2/FR-2.3. Returns a result on success, std::nullopt if no
// retained-system match is found (caller should fall back to the
// existing monocyclic fusion logic).
std::optional<RingSystemNameResult> tryRetainedBaseDecomposition(
    int mol,
    const std::vector<std::set<int>> *subsetRings = nullptr,
    std::vector<RingSystemNameResult> *outSymmetryCandidates = nullptr);
