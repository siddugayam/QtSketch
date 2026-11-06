#pragma once
#include <QObject>
#include <QPointF>
#include <QString>
#include <QList>
#include <QVariantMap>

struct PreviewAtom { 
    QPointF pos; 
    QString label; 
};

struct PreviewBond { 
    QPointF start; 
    QPointF end; 
};

struct PlacementResult {
    bool valid = false;
    QList<PreviewAtom> atoms;
    QList<PreviewBond> bonds;
};

class PlacementUtils;

class AtomPlacementEngine {
public:
    static PlacementResult compute(const QPointF& startPos, const QPointF& currentMouse, 
                                   const QString& toolAtom, const QList<double>& existingAngles, 
                                   double bondLength);
};

class FragmentPlacementEngine {
public:
    static PlacementResult compute(const QPointF& startPos, const QPointF& currentMouse,
                                   const QString& fragmentId, const QList<double>& existingAngles,
                                   double bondLength);
};

// Mirrors v8_worker.js's addChain(x1,y1,x2,y2) geometry exactly (same nBonds/theta/half-angle
// formula) so the live drag preview matches what actually gets committed. Preview-only; the
// commit path still calls V8Process::addChain directly rather than consuming this result, since
// addChain already recomputes the identical geometry worker-side and additionally fuses
// overlapping atoms - duplicating that here would risk the preview and the committed structure
// drifting apart.
class ChainPlacementEngine {
public:
    static PlacementResult compute(const QPointF& startChem, const QPointF& currentChem, double bondLength);
};
