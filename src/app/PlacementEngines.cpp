#include "PlacementEngines.h"
#include "PlacementUtils.h"
#include <cmath>
#include <QDebug>

static const double PI = 3.14159265358979323846;

PlacementResult AtomPlacementEngine::compute(const QPointF& startPos, const QPointF& currentMouse, 
                                             const QString& toolAtom, const QList<double>& existingAngles, 
                                             double bondLength) {
    PlacementResult res;
    double dx = currentMouse.x() - startPos.x();
    double dy = currentMouse.y() - startPos.y();
    
    if (std::abs(dx) < 0.01 && std::abs(dy) < 0.01) {
        // Not enough drag, just default right
        dx = 1.0; dy = 0.0;
    }
    
    double angle = std::atan2(dy, dx);
    angle = PlacementUtils::snapAngle(angle, existingAngles);
    
    QPointF newPos(startPos.x() + std::cos(angle) * bondLength, 
                   startPos.y() + std::sin(angle) * bondLength);
                   
    res.atoms.append({newPos, toolAtom});
    res.bonds.append({startPos, newPos});
    res.valid = true;
    return res;
}

PlacementResult FragmentPlacementEngine::compute(const QPointF& startPos, const QPointF& currentMouse, 
                                                 const QString& fragmentId, const QList<double>& existingAngles, 
                                                 double bondLength) {
    PlacementResult res;
    double dx = currentMouse.x() - startPos.x();
    double dy = currentMouse.y() - startPos.y();
    bool noDrag = false;
    if (std::abs(dx) < 0.01 && std::abs(dy) < 0.01) { 
        dx = 1.0; dy = 0.0; 
        noDrag = true;
    }
    
    double rawAngle = std::atan2(dy, dx);
    
    // Use the same snapping logic as AtomPlacementEngine (with steric avoidance)
    double snappedAngle = PlacementUtils::snapAngle(rawAngle, existingAngles);
    
    // For no-drag placement (plain click), rotate by 90 degrees for vertex-up orientation
    // (standard chemistry drawing convention: vertex at top/bottom, not flat-top)
    if (noDrag) {
        snappedAngle += PI / 2.0;
    }

    // fragmentId is always a TEMPLATE_ id now -- PlacementPreviewManager::updatePreview routes
    // FG_/SS_/LIB_ tool ids to V8Process::getFunctionalGroupPlacementResult instead (real
    // template geometry, graft-aware), not to this function. The single-fake-atom fallback that
    // used to live here for FG_/SS_/LIB_ ids has been removed.
    int n = 6; // Default
    if (fragmentId == "TEMPLATE_BENZENE") {
        n = 6;
    } else {
        QString sizeStr = fragmentId.mid(9);
        bool ok;
        int parsed = sizeStr.toInt(&ok);
        if (ok && parsed >= 3 && parsed <= 24) n = parsed;
    }

    // Circumradius R = s / (2 * sin(PI / n))
    double R = bondLength / (2.0 * std::sin(PI / n));

    // We want the attachment atom at startPos, and the center of the ring at startPos + direction * R
    QPointF center(startPos.x() + std::cos(snappedAngle) * R,
                   startPos.y() + std::sin(snappedAngle) * R);

    // Generate the vertices of the regular n-gon
    for (int i = 0; i < n; ++i) {
        double a = snappedAngle + PI + (i * 2 * PI / n);
        QPointF p(center.x() + std::cos(a) * R, center.y() + std::sin(a) * R);
        res.atoms.append({p, "C"});
    }
    for (int i = 0; i < n; ++i) {
        res.bonds.append({res.atoms[i].pos, res.atoms[(i+1)%n].pos});
    }
    res.valid = true;

    return res;
}

PlacementResult ChainPlacementEngine::compute(const QPointF& startChem, const QPointF& currentChem, double bondLength) {
    PlacementResult res;
    if (bondLength <= 0) { res.valid = false; return res; }

    double dx = currentChem.x() - startChem.x();
    double dy = currentChem.y() - startChem.y();
    double dist = std::sqrt(dx * dx + dy * dy);
    int nBonds = std::min(200, std::max(1, (int)std::round(dist / bondLength)));
    double theta = std::atan2(dy, dx);
    double half = PI / 6.0;

    double px = startChem.x(), py = startChem.y();
    bool havePrev = false;
    QPointF prev;
    for (int i = 0; i <= nBonds; ++i) {
        QPointF atomPos(px, py);
        res.atoms.append({atomPos, QString()});  // unlabeled - carbon vertices render without letters
        if (havePrev) {
            res.bonds.append({prev, atomPos});
        }
        prev = atomPos;
        havePrev = true;
        double ang = theta + ((i % 2 == 0) ? half : -half);
        px += bondLength * std::cos(ang);
        py += bondLength * std::sin(ang);
    }
    res.valid = true;
    return res;
}
