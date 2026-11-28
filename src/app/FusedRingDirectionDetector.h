#pragma once

#include "FusedRingOrientation.h"

// mol is an Indigo molecule handle.
// Returns a FusedRingSystemInput if the molecule contains a valid 
// linear chain of 6-membered ortho-fused rings. 
// Otherwise returns an empty/invalid FusedRingSystemInput (with ringSizes empty).
FusedRingSystemInput detectFusedRingDirections(int mol, const std::vector<std::set<int>> *ringsOverride = nullptr);
