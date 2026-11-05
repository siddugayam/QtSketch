// src/app/molecule/PageSettings.cpp
#include "PageSettings.h"

PageSettings::PageSettings(double widthMm, double heightMm,
                           double marginLeftMm, double marginRightMm,
                           double marginTopMm, double marginBottomMm)
    : widthMm(widthMm)
    , heightMm(heightMm)
    , marginLeftMm(marginLeftMm)
    , marginRightMm(marginRightMm)
    , marginTopMm(marginTopMm)
    , marginBottomMm(marginBottomMm)
{
}

bool PageSettings::operator==(const PageSettings& other) const {
    return widthMm == other.widthMm &&
           heightMm == other.heightMm &&
           marginLeftMm == other.marginLeftMm &&
           marginRightMm == other.marginRightMm &&
           marginTopMm == other.marginTopMm &&
           marginBottomMm == other.marginBottomMm;
}

PageSettings PageSettings::defaults() {
    return PageSettings(210.0, 297.0);  // A4
}

bool PageSettings::isWithinBounds(double x, double y, double modelToMmScale) const {
    // Convert model coordinates to mm
    double xMm = x * modelToMmScale;
    double yMm = y * modelToMmScale;
    
    // Check against page dimensions with margins
    return xMm >= marginLeftMm && xMm <= (widthMm - marginRightMm) &&
           yMm >= marginBottomMm && yMm <= (heightMm - marginTopMm);
}

void PageSettings::getUsableArea(double* minX, double* minY, double* maxX, double* maxY, 
                                  double modelToMmScale) const {
    // Convert page bounds from mm to model coordinates
    *minX = marginLeftMm / modelToMmScale;
    *minY = marginBottomMm / modelToMmScale;
    *maxX = (widthMm - marginRightMm) / modelToMmScale;
    *maxY = (heightMm - marginTopMm) / modelToMmScale;
}
