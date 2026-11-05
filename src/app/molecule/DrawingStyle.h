// src/app/molecule/DrawingStyle.h
#ifndef DRAWINGSTYLE_H
#define DRAWINGSTYLE_H

#include <QString>

// Document-level drawing style settings that control artwork generation.
// These settings are serialized as part of the .qtsketch project format.
//
// Design decisions:
// - Bond length is in model coordinates (not pixels)
// - Line width is in points (1/72 inch)
// - Font family and size control text rendering
// - Colors are for publication output (monochrome by default per plan defaults)
// - Fixed physical scale: model coordinates map to physical dimensions via this style

class DrawingStyle {
public:
    DrawingStyle() = default;
    explicit DrawingStyle(double bondLength, double lineWidth, 
                         const QString& fontFamily, int fontSize,
                         bool monochrome);

    // Bond length in model coordinates (default: 1.5 model units = 5.08 mm at 1:1 scale)
    double bondLength = 1.5;

    // Line width in points (default: 0.6 pt)
    double lineWidth = 0.6;

    // Font settings (default: Arial 10 pt)
    QString fontFamily = QStringLiteral("Arial");
    int fontSize = 10;

    // Monochrome output for publication (default: true)
    bool monochrome = true;

    // Physical scale: model coordinates to millimeters
    // At default settings: 1.5 model units = 5.08 mm
    // So scale = 5.08 / 1.5 = 3.386666... mm per model unit
    // This is computed once during document creation from the legacy scale
    // and then fixed for the document's lifetime.
    double modelToMmScale = 5.08 / 1.5;

    // Equality comparison
    bool operator==(const DrawingStyle& other) const;
    bool operator!=(const DrawingStyle& other) const { return !(*this == other); }

    // Create with defaults per plan
    static DrawingStyle defaults();

    // Legacy scale capture: compute the fixed scale from an existing bond length
    // This is called once when importing existing work
    void captureLegacyScale(double legacyBondLengthMm);
};

#endif // DRAWINGSTYLE_H
