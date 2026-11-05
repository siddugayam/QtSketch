// src/app/molecule/PageSettings.h
#ifndef PAGESETTINGS_H
#define PAGESETTINGS_H

#include <QString>

// Document-level page settings that control output dimensions and margins.
// These settings are serialized as part of the .qtsketch project format.
//
// Design decisions:
// - Page size is in millimeters
// - Margins are in millimeters
// - Default to A4 page size (210 x 297 mm) per the plan's context

class PageSettings {
public:
    PageSettings() = default;
    explicit PageSettings(double widthMm, double heightMm,
                         double marginLeftMm = 0, double marginRightMm = 0,
                         double marginTopMm = 0, double marginBottomMm = 0);

    // Page dimensions in millimeters (default: A4 = 210 x 297 mm)
    double widthMm = 210.0;
    double heightMm = 297.0;

    // Margins in millimeters (default: 0)
    double marginLeftMm = 0.0;
    double marginRightMm = 0.0;
    double marginTopMm = 0.0;
    double marginBottomMm = 0.0;

    // Equality comparison
    bool operator==(const PageSettings& other) const;
    bool operator!=(const PageSettings& other) const { return !(*this == other); }

    // Create with defaults
    static PageSettings defaults();

    // Check if a point (in model coordinates) is within the page bounds
    // considering margins. Uses the provided scale to convert model to mm.
    bool isWithinBounds(double x, double y, double modelToMmScale) const;

    // Get the usable page area in model coordinates
    void getUsableArea(double* minX, double* minY, double* maxX, double* maxY, 
                       double modelToMmScale) const;
};

#endif // PAGESETTINGS_H
