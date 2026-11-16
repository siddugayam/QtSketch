#pragma once

#include "GraphTypes.h"
#include "FusedRingOrientation.h"
#include "NamingTrace.h"
#include "SystemLocant.h"

#include <map>
#include <set>
#include <vector>
#include <QString>

std::map<int, QString> computePeripheralNumberingGeneral(
    const Graph &g,
    const FusedRingSystemInput &input,
    const std::set<int> &substituentBearingNodes = {},
    const std::set<int> &hydroLocantBearingNodes = {},
    std::vector<std::map<int, QString>> *tiedMapsOut = nullptr,
    bool restrictStartToPredecessorBridgehead = false,
    const std::set<int> &principalGroupBearingNodes = {});

// Number interior atoms per FR-5.5.2 / P-25.3.3.1.2-.3.
// Ordering contract: returns one SystemLocant per interior atom in the exact iteration order of interiorAtoms.
std::vector<SystemLocant> numberInteriorAtoms(
    const std::vector<int> &perimeterOrdering,
    const Graph &g,
    const std::set<int> &interiorAtoms);
