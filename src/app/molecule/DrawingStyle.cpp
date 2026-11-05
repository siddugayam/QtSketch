// src/app/molecule/DrawingStyle.cpp
#include "DrawingStyle.h"

DrawingStyle::DrawingStyle(double bondLength, double lineWidth,
                           const QString& fontFamily, int fontSize,
                           bool monochrome)
    : bondLength(bondLength)
    , lineWidth(lineWidth)
    , fontFamily(fontFamily)
    , fontSize(fontSize)
    , monochrome(monochrome)
    , modelToMmScale(5.08 / bondLength)  // Default physical scale
{
}

bool DrawingStyle::operator==(const DrawingStyle& other) const {
    return bondLength == other.bondLength &&
           lineWidth == other.lineWidth &&
           fontFamily == other.fontFamily &&
           fontSize == other.fontSize &&
           monochrome == other.monochrome &&
           modelToMmScale == other.modelToMmScale;
}

DrawingStyle DrawingStyle::defaults() {
    return DrawingStyle(1.5, 0.6, QStringLiteral("Arial"), 10, true);
}

void DrawingStyle::captureLegacyScale(double legacyBondLengthMm) {
    // Compute the fixed physical scale from the legacy bond length
    // This ensures existing coordinates are not reinterpreted when style changes
    if (legacyBondLengthMm > 0) {
        modelToMmScale = legacyBondLengthMm / bondLength;
    }
}
