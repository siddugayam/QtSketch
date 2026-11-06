// src/app/PlacementUtils.h
#pragma once
#include <QList>
#include <cmath>

// Shared placement utilities for snapping and steric avoidance
// Used by AtomPlacementEngine and FragmentPlacementEngine

class PlacementUtils {
public:
    // Constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double SNAP_INTERVAL = PI / 6.0;      // 30 degrees
    static constexpr double STERIC_THRESHOLD = PI / 12.0;  // 15 degrees
    static constexpr int MAX_NUDGES = 6;                 // Maximum nudge attempts

    // Snap an angle to 30-degree increments with steric avoidance
    // Returns the best snapped angle that avoids existing bond angles
    // rawAngle: the raw drag angle in radians
    // existingAngles: list of existing bond angles (in radians) from the start atom
    static double snapAngle(double rawAngle, const QList<double>& existingAngles);

    // Normalize angle to range [-PI, PI]
    static double normalizeAngle(double angle);

    // Compute angle difference (smallest angle between two angles)
    static double angleDifference(double a1, double a2);

    // Check if an angle is too close to any existing angle
    static bool isAngleTooClose(double angle, const QList<double>& existingAngles, double threshold = STERIC_THRESHOLD);
};
