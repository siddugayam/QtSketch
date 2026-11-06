// src/app/PlacementUtils.cpp
#include "PlacementUtils.h"
#include <algorithm>

// Normalize angle to range [-PI, PI]
double PlacementUtils::normalizeAngle(double angle) {
    while (angle > PI) angle -= 2 * PI;
    while (angle < -PI) angle += 2 * PI;
    return angle;
}

// Compute angle difference (smallest angle between two angles)
double PlacementUtils::angleDifference(double a1, double a2) {
    double diff = std::abs(normalizeAngle(a1) - normalizeAngle(a2));
    return std::min(diff, 2 * PI - diff);
}

// Check if an angle is too close to any existing angle
bool PlacementUtils::isAngleTooClose(double angle, const QList<double>& existingAngles, double threshold) {
    double normalized = normalizeAngle(angle);
    for (double ang : existingAngles) {
        if (angleDifference(normalized, ang) < threshold) {
            return true;
        }
    }
    return false;
}

// Snap an angle to 30-degree increments with steric avoidance
// Uses iterative nudging: keeps trying angles until one clears all existing angles
// or reaches MAX_NUDGES, then returns the best candidate found
double PlacementUtils::snapAngle(double rawAngle, const QList<double>& existingAngles) {
    // Basic snapping to 30° increments
    double snapInterval = SNAP_INTERVAL;
    double normalizedRaw = normalizeAngle(rawAngle);
    double snapped = std::round(normalizedRaw / snapInterval) * snapInterval;

    // If no existing angles, just return the snapped angle
    if (existingAngles.isEmpty()) {
        return snapped;
    }

    // Steric avoidance: find an angle that clears all existing angles
    // Try the snapped angle first, then nudge in both directions
    double bestCandidate = snapped;
    bool foundClear = !isAngleTooClose(snapped, existingAngles);

    if (!foundClear) {
        // Try nudging in positive direction first
        for (int i = 1; i <= MAX_NUDGES; ++i) {
            double candidate = snapped + i * snapInterval;
            if (!isAngleTooClose(candidate, existingAngles)) {
                bestCandidate = candidate;
                foundClear = true;
                break;
            }
        }

        // If still not clear, try negative direction
        if (!foundClear) {
            for (int i = 1; i <= MAX_NUDGES; ++i) {
                double candidate = snapped - i * snapInterval;
                if (!isAngleTooClose(candidate, existingAngles)) {
                    bestCandidate = candidate;
                    foundClear = true;
                    break;
                }
            }
        }

        // If we still haven't found a clear angle, use the snapped angle anyway
        // (better than nothing)
    }

    return normalizeAngle(bestCandidate);
}
