// src/app/molecule/DocumentState.cpp
#include "DocumentState.h"
#include <memory>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <QPointF>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include "indigo.h"
#include <optional>
#include "BondAngleSuggester.h"
#include "PlacementCleanup.h"

// RAII guard for m_inCommand: execute/undo/redo previously set the flag true, ran an
// arbitrary EditCommand closure, then set it false -- if that closure ever threw (Indigo
// call, bad_alloc, ...), the flag stuck true forever and every future
// executeCommand/undo/redo silently no-op'd (their own `if (m_inCommand) return;` guards),
// bricking undo/redo for the rest of the session. The destructor resets the flag during
// stack unwinding too, so a throw still propagates to the caller but no longer leaves the
// document wedged.
namespace {
class InCommandGuard {
public:
    explicit InCommandGuard(bool& flag) : m_flag(flag) { m_flag = true; }
    ~InCommandGuard() { m_flag = false; }
private:
    bool& m_flag;
};
} // namespace

// ---- Pure geometry helpers for selectByRect/addSelectionByRect/selectByLasso ------
// Plain file-scope functions, not DocumentState members -- matches RenderPrimitives.cpp's
// own convention for its pure helpers (isCorrectStereoCenter, stereoLabelFor). Nothing
// outside this file needs them.

// Direct port of chem-core.js:9059-9065's isPointOnSegment. Only ever called from
// segmentsIntersect's own degenerate (all-four-points-collinear) branch below -- despite
// the name, this is a bounding-box containment check, not a true collinearity test (the
// real implementation's own behavior, kept exactly, not "fixed").
static bool isPointOnSegment(QPointF segA, QPointF segB, QPointF point) {
    double minX = std::min(segA.x(), segB.x()), maxX = std::max(segA.x(), segB.x());
    double minY = std::min(segA.y(), segB.y()), maxY = std::max(segA.y(), segB.y());
    return point.x() >= minX && point.x() <= maxX && point.y() >= minY && point.y() <= maxY;
}

// Direct port of chem-core.js:9047-9056's Box2Abs.segmentIntersection. Cross-product
// orientation test; falls back to isPointOnSegment above only in the exact-collinear
// degenerate case (all four cross products are exactly 0.0).
static bool segmentsIntersect(QPointF a, QPointF b, QPointF c, QPointF d) {
    double dc = (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
    double dd = (a.x() - d.x()) * (b.y() - d.y()) - (a.y() - d.y()) * (b.x() - d.x());
    double da = (c.x() - a.x()) * (d.y() - a.y()) - (c.y() - a.y()) * (d.x() - a.x());
    double db = (c.x() - b.x()) * (d.y() - b.y()) - (c.y() - b.y()) * (d.x() - b.x());
    if (dc == 0.0 && dd == 0.0 && da == 0.0 && db == 0.0) {
        return isPointOnSegment(a, b, c) || isPointOnSegment(a, b, d)
            || isPointOnSegment(c, d, a) || isPointOnSegment(c, d, b);
    }
    return dc * dd < 0 && da * db < 0;
}

// Direct port of selectByRect's own lineIntersectsRect (10-state.js:690-697): true if
// either endpoint lies inside the rect, or the segment crosses any of its 4 edges.
static bool segmentIntersectsRect(QPointF p1, QPointF p2, double minX, double minY, double maxX, double maxY) {
    if (p1.x() >= minX && p1.x() <= maxX && p1.y() >= minY && p1.y() <= maxY) return true;
    if (p2.x() >= minX && p2.x() <= maxX && p2.y() >= minY && p2.y() <= maxY) return true;
    QPointF c1(minX, minY), c2(maxX, minY), c3(maxX, maxY), c4(minX, maxY);
    return segmentsIntersect(p1, p2, c1, c2) || segmentsIntersect(p1, p2, c2, c3)
        || segmentsIntersect(p1, p2, c3, c4) || segmentsIntersect(p1, p2, c4, c1);
}

// Point-in-polygon via ray casting (even-odd rule), direct port of the real
// _pointInPolygon (10-state.js:735-745). poly is a list of (x,y) vertices.
static bool pointInPolygon(double px, double py, const QList<QPointF>& poly) {
    bool inside = false;
    for (int i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        double xi = poly[i].x(), yi = poly[i].y();
        double xj = poly[j].x(), yj = poly[j].y();
        bool intersect = ((yi > py) != (yj > py))
            && (px < (xj - xi) * (py - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
}

// Direct port of 10-state.js's shortestRingThroughBond (129-156): BFS over the whole bond
// graph WITH bondId itself excluded, from that bond's begin atom to its end atom. Returns
// an empty list if bondId doesn't exist or isn't part of any ring. -1 is used as the BFS
// "root" sentinel in the prev-map (matches this port's own null-id convention -- no real
// AtomId is ever -1), mirroring the real JS's `prev[start] = null`.
//
// Named bfsShortestRing (not shortestRingThroughBond) to avoid colliding with the
// PRE-EXISTING private DocumentState::shortestRingThroughBond(int moleculeHandle, int
// bondIdx) member function declared near the bottom of DocumentState.h (used by
// insertLibraryTemplateFused) -- that one operates on raw Indigo handles/indices, this one
// on this port's own EditableMolecule/AtomId/BondId wrapper types. Genuinely different
// functions serving different callers, not a duplicate to merge.
static QList<AtomId> bfsShortestRing(const EditableMolecule& mol, BondId bondId) {
    AtomId start = -1, goal = -1;
    if (!mol.bondEndpoints(bondId, start, goal)) return {};

    QHash<AtomId, QList<AtomId>> adj;
    for (BondId bid : mol.bondIds()) {
        if (bid == bondId) continue;
        AtomId a = -1, b = -1;
        if (!mol.bondEndpoints(bid, a, b)) continue;
        adj[a].append(b);
        adj[b].append(a);
    }

    QHash<AtomId, AtomId> prev;
    prev.insert(start, -1);
    QList<AtomId> queue;
    queue.append(start);
    int head = 0;
    while (head < queue.size()) {
        AtomId cur = queue[head++];
        if (cur == goal) break;
        for (AtomId nbr : adj.value(cur)) {
            if (!prev.contains(nbr)) {
                prev.insert(nbr, cur);
                queue.append(nbr);
            }
        }
    }
    if (!prev.contains(goal)) return {};

    QList<AtomId> cycle;
    for (AtomId a = goal; a != -1; a = prev.value(a, -1)) cycle.append(a);
    return cycle;
}

DocumentState::DocumentState(const QString& initialStructure)
    : m_molecule(initialStructure)
    , m_drawingStyle(DrawingStyle::defaults())
    , m_pageSettings(PageSettings::defaults()) {
}

EditableMolecule& DocumentState::molecule() {
    return m_molecule;
}

const EditableMolecule& DocumentState::molecule() const {
    return m_molecule;
}

void DocumentState::executeCommand(EditCommand cmd) {
    if (m_inCommand) return;
    // Remove any future redo states -- matches executeCommand's
    // `_history.splice(_historyPointer + 1)` exactly.
    if (m_historyPointer + 1 < static_cast<int>(m_history.size())) {
        m_history.resize(m_historyPointer + 1);
    }
    {
        InCommandGuard guard(m_inCommand);
        cmd.execute();
    }
    reconcileSelectionAfterCommand();
    m_history.push_back(std::move(cmd));
    ++m_historyPointer;
    if (static_cast<int>(m_history.size()) > kHistorySize) {
        m_history.erase(m_history.begin());
        --m_historyPointer;
    }
    m_dirty = true;
}

void DocumentState::undo() {
    if (m_inCommand) return;
    if (m_historyPointer < 0) return;
    {
        InCommandGuard guard(m_inCommand);
        m_history[m_historyPointer].invert();
    }
    --m_historyPointer;
    m_selection.clear();
    resetAllDragState();
    m_dirty = true;
}

void DocumentState::redo() {
    if (m_inCommand) return;
    if (m_historyPointer >= static_cast<int>(m_history.size()) - 1) return;
    ++m_historyPointer;
    {
        InCommandGuard guard(m_inCommand);
        m_history[m_historyPointer].execute();
    }
    m_selection.clear();
    resetAllDragState();
    m_dirty = true;
}

void DocumentState::reconcileSelectionAfterCommand() {
    if (!m_selection.atoms.isEmpty()) {
        QList<AtomId> allAtoms = m_molecule.atomIds();
        QList<SGroupId> allSgroups = m_molecule.sgroupIds();
        QSet<AtomId> validAtoms(allAtoms.begin(), allAtoms.end());
        QSet<SGroupId> validSgroups(allSgroups.begin(), allSgroups.end());
        QSet<AtomId> kept;
        for (AtomId id : m_selection.atoms) {
            if (validAtoms.contains(id) || validSgroups.contains(id)) kept.insert(id);
        }
        m_selection.atoms = kept;
    }
    if (!m_selection.bonds.isEmpty()) {
        QList<BondId> allBonds = m_molecule.bondIds();
        QSet<BondId> validBonds(allBonds.begin(), allBonds.end());
        QSet<BondId> kept;
        for (BondId id : m_selection.bonds) {
            if (validBonds.contains(id)) kept.insert(id);
        }
        m_selection.bonds = kept;
    }
    if (!m_selection.rxnArrows.isEmpty()) {
        QList<RxnArrowId> allArrows = m_molecule.rxnArrowIds();
        QSet<RxnArrowId> validArrows(allArrows.begin(), allArrows.end());
        QSet<RxnArrowId> kept;
        for (RxnArrowId id : m_selection.rxnArrows) {
            if (validArrows.contains(id)) kept.insert(id);
        }
        m_selection.rxnArrows = kept;
    }
    if (!m_selection.rxnPluses.isEmpty()) {
        QList<RxnPlusId> allPluses = m_molecule.rxnPlusIds();
        QSet<RxnPlusId> validPluses(allPluses.begin(), allPluses.end());
        QSet<RxnPlusId> kept;
        for (RxnPlusId id : m_selection.rxnPluses) {
            if (validPluses.contains(id)) kept.insert(id);
        }
        m_selection.rxnPluses = kept;
    }
    if (!m_selection.multitailArrows.isEmpty()) {
        QList<MultitailArrowId> allMtas = m_molecule.multitailArrowIds();
        QSet<MultitailArrowId> validMtas(allMtas.begin(), allMtas.end());
        QSet<MultitailArrowId> kept;
        for (MultitailArrowId id : m_selection.multitailArrows) {
            if (validMtas.contains(id)) kept.insert(id);
        }
        m_selection.multitailArrows = kept;
    }
}

bool DocumentState::isDirty() const {
    return m_dirty;
}

void DocumentState::markClean() {
    m_dirty = false;
}

void DocumentState::clearHistory() {
    m_history.clear();
    m_historyPointer = -1;
}

void DocumentState::setDrawingStyle(const DrawingStyle& style) {
    DrawingStyle newStyle = style;
    // Maintain fixed physical scale invariant: preserve current modelToMmScale
    // if incoming style doesn't provide a valid positive scale
    if (newStyle.modelToMmScale <= 0) {
        newStyle.modelToMmScale = m_drawingStyle.modelToMmScale;
    }
    if (m_drawingStyle == newStyle) {
        return;
    }

    DrawingStyle oldStyle = m_drawingStyle;
    EditCommand cmd;
    cmd.execute = [this, newStyle]() {
        m_drawingStyle = newStyle;
    };
    cmd.invert = [this, oldStyle]() {
        m_drawingStyle = oldStyle;
    };
    executeCommand(std::move(cmd));
}

void DocumentState::applyDrawingStyleWithRescale(const DrawingStyle& style) {
    // The rescale factor is new/old bond length and its inverse is what undo applies: a
    // zero or negative length would collapse every atom onto the centroid and leave a
    // zero divisor behind for the next call. Reject it outright.
    if (!(style.bondLength > 0.0) || !(m_drawingStyle.bondLength > 0.0)) return;

    // Check if bond length is changing
    if (m_drawingStyle.bondLength != style.bondLength) {
        // Compute scale factor: new/old bond length ratio
        double scaleFactor = style.bondLength / m_drawingStyle.bondLength;
        
        // Compute centroid as the scaling center
        QPointF centroid = m_molecule.computeCentroid();
        
        // Create a command that captures both the style change and coordinate rescale
        // We need to capture old state for undo
        DrawingStyle oldStyle = m_drawingStyle;
        
        // Save old atom positions
        QHash<AtomId, QPointF> oldPositions;
        for (AtomId id : m_molecule.atomIds()) {
            double x, y;
            if (m_molecule.atomPos(id, x, y)) {
                oldPositions[id] = QPointF(x, y);
            }
        }
        
        // Create the command
        EditCommand cmd;
        cmd.execute = [this, style, scaleFactor, centroid]() {
            // Rescale all coordinates
            m_molecule.scaleAllCoordinates(scaleFactor, centroid.x(), centroid.y());

            // modelToMmScale is fixed for the document's lifetime (style already
            // carries the document's true scale, seeded via SettingsService's
            // draft from the real DrawingStyle). The coordinate rescale above,
            // combined with the unchanged scale, is what makes a bond-length
            // change actually affect the final physical/export size -- do NOT
            // re-capture a legacy scale here, that would cancel the rescale out.
            m_drawingStyle = style;
        };
        cmd.invert = [this, oldStyle, oldPositions]() {
            // Restore old positions
            for (auto it = oldPositions.constBegin(); it != oldPositions.constEnd(); ++it) {
                m_molecule.setAtomPos(it.key(), it.value().x(), it.value().y());
            }
            
            // Restore old style
            m_drawingStyle = oldStyle;
        };
        
        executeCommand(std::move(cmd));
    } else {
        // Bond length not changing, apply the style (non-bond-length changes like font,
        // lineWidth, monochrome are wrapped in undo/dirty tracking inside setDrawingStyle)
        setDrawingStyle(style);
    }
}

void DocumentState::setPageSettings(const PageSettings& settings) {
    if (m_pageSettings == settings) {
        return;
    }

    PageSettings oldSettings = m_pageSettings;
    PageSettings newSettings = settings;
    EditCommand cmd;
    cmd.execute = [this, newSettings]() {
        m_pageSettings = newSettings;
    };
    cmd.invert = [this, oldSettings]() {
        m_pageSettings = oldSettings;
    };
    executeCommand(std::move(cmd));
}

bool DocumentState::canUndo() const {
    return m_historyPointer >= 0;
}

bool DocumentState::canRedo() const {
    return m_historyPointer < static_cast<int>(m_history.size()) - 1;
}

void DocumentState::setShowExplicitH(bool show) {
    m_showExplicitH = show;
}

bool DocumentState::showExplicitH() const {
    return m_showExplicitH;
}

void DocumentState::selectAtom(AtomId id) {
    m_selection.clear();
    m_selection.atoms.insert(id);
}
void DocumentState::selectBond(BondId id) {
    m_selection.clear();
    m_selection.bonds.insert(id);
}
void DocumentState::selectRxnArrow(RxnArrowId id) {
    m_selection.clear();
    m_selection.rxnArrows.insert(id);
}
void DocumentState::selectRxnPlus(RxnPlusId id) {
    m_selection.clear();
    m_selection.rxnPluses.insert(id);
}
void DocumentState::selectMultitailArrow(MultitailArrowId id) {
    m_selection.clear();
    m_selection.multitailArrows.insert(id);
}
void DocumentState::clearSelection() {
    m_selection.clear();
}

void DocumentState::addAtomToSelection(AtomId id) { m_selection.atoms.insert(id); }
void DocumentState::addBondToSelection(BondId id) { m_selection.bonds.insert(id); }
void DocumentState::addRxnArrowToSelection(RxnArrowId id) { m_selection.rxnArrows.insert(id); }
void DocumentState::addRxnPlusToSelection(RxnPlusId id) { m_selection.rxnPluses.insert(id); }
void DocumentState::addMultitailArrowToSelection(MultitailArrowId id) { m_selection.multitailArrows.insert(id); }

void DocumentState::removeAtomFromSelection(AtomId id) { m_selection.atoms.remove(id); }
void DocumentState::removeBondFromSelection(BondId id) { m_selection.bonds.remove(id); }
void DocumentState::removeRxnArrowFromSelection(RxnArrowId id) { m_selection.rxnArrows.remove(id); }
void DocumentState::removeRxnPlusFromSelection(RxnPlusId id) { m_selection.rxnPluses.remove(id); }
void DocumentState::removeMultitailArrowFromSelection(MultitailArrowId id) { m_selection.multitailArrows.remove(id); }

void DocumentState::selectAll() {
    m_selection.clear();
    for (AtomId id : m_molecule.atomIds()) m_selection.atoms.insert(id);
    for (BondId id : m_molecule.bondIds()) m_selection.bonds.insert(id);
    for (RxnArrowId id : m_molecule.rxnArrowIds()) m_selection.rxnArrows.insert(id);
    for (RxnPlusId id : m_molecule.rxnPlusIds()) m_selection.rxnPluses.insert(id);
    for (MultitailArrowId id : m_molecule.multitailArrowIds()) m_selection.multitailArrows.insert(id);
}

SelectionState& DocumentState::selection() {
    return m_selection;
}

void DocumentState::selectByRect(double x1, double y1, double x2, double y2) {
    double minX = std::min(x1, x2), maxX = std::max(x1, x2);
    double minY = std::min(y1, y2), maxY = std::max(y1, y2);
    m_selection.clear();

    for (AtomId id : m_molecule.atomIds()) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        if (x >= minX && x <= maxX && y >= minY && y <= maxY) m_selection.atoms.insert(id);
    }
    for (BondId id : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(id, ea, eb)) continue;
        if (m_selection.atoms.contains(ea) && m_selection.atoms.contains(eb)) {
            m_selection.bonds.insert(id);
            continue;
        }
        double x1p = 0, y1p = 0, x2p = 0, y2p = 0;
        if (m_molecule.atomPos(ea, x1p, y1p) && m_molecule.atomPos(eb, x2p, y2p)
            && segmentIntersectsRect(QPointF(x1p, y1p), QPointF(x2p, y2p), minX, minY, maxX, maxY)) {
            m_selection.bonds.insert(id);
        }
    }
    for (RxnArrowId id : m_molecule.rxnArrowIds()) {
        double x1p = 0, y1p = 0, x2p = 0, y2p = 0;
        if (!m_molecule.rxnArrowEndpoints(id, x1p, y1p, x2p, y2p)) continue;
        if (segmentIntersectsRect(QPointF(x1p, y1p), QPointF(x2p, y2p), minX, minY, maxX, maxY)) {
            m_selection.rxnArrows.insert(id);
        }
    }
    for (RxnPlusId id : m_molecule.rxnPlusIds()) {
        double x = 0, y = 0;
        if (!m_molecule.rxnPlusPos(id, x, y)) continue;
        if (x >= minX && x <= maxX && y >= minY && y <= maxY) m_selection.rxnPluses.insert(id);
    }
    // Multitail arrows deliberately excluded -- matches the real code's own comment
    // ("no rect-intersection hit-test exists for their spine+tails geometry").
}

void DocumentState::addSelectionByRect(double x1, double y1, double x2, double y2) {
    double minX = std::min(x1, x2), maxX = std::max(x1, x2);
    double minY = std::min(y1, y2), maxY = std::max(y1, y2);

    for (AtomId id : m_molecule.atomIds()) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        if (x >= minX && x <= maxX && y >= minY && y <= maxY) m_selection.atoms.insert(id);
    }
    for (BondId id : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(id, ea, eb)) continue;
        double x1p = 0, y1p = 0, x2p = 0, y2p = 0;
        if (!m_molecule.atomPos(ea, x1p, y1p) || !m_molecule.atomPos(eb, x2p, y2p)) continue;
        if (x1p >= minX && x1p <= maxX && y1p >= minY && y1p <= maxY
            && x2p >= minX && x2p <= maxX && y2p >= minY && y2p <= maxY) {
            m_selection.bonds.insert(id);
        }
    }
    // rxnArrows/rxnPluses/multitailArrows deliberately untouched -- matches the real
    // addSelectionByRect exactly (it never references them at all).
}

void DocumentState::selectByLasso(const QList<QPointF>& points) {
    m_selection.clear();
    if (points.size() < 3) return;

    for (AtomId id : m_molecule.atomIds()) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        if (pointInPolygon(x, y, points)) m_selection.atoms.insert(id);
    }
    for (BondId id : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(id, ea, eb)) continue;
        if (m_selection.atoms.contains(ea) && m_selection.atoms.contains(eb)) {
            m_selection.bonds.insert(id);
        }
    }
    for (RxnArrowId id : m_molecule.rxnArrowIds()) {
        double x1p = 0, y1p = 0, x2p = 0, y2p = 0;
        if (!m_molecule.rxnArrowEndpoints(id, x1p, y1p, x2p, y2p)) continue;
        if (pointInPolygon(x1p, y1p, points) && pointInPolygon(x2p, y2p, points)) {
            m_selection.rxnArrows.insert(id);
        }
    }
    for (RxnPlusId id : m_molecule.rxnPlusIds()) {
        double x = 0, y = 0;
        if (!m_molecule.rxnPlusPos(id, x, y)) continue;
        if (pointInPolygon(x, y, points)) m_selection.rxnPluses.insert(id);
    }
    // Multitail arrows deliberately excluded, same reasoning as selectByRect.
}

void DocumentState::selectFragment(AtomId atomId, BondId bondId) {
    // If bondId doesn't resolve to a real bond, atomId stays -1 and falls straight into
    // the clear-branch below -- same outcome as both args being null, matching the real
    // code's `if (atomId === null) { clear; return }` running AFTER this resolution
    // attempt (verified during spec review round 2: this is NOT the same as the
    // non-null-but-invalid no-op case further down).
    if (atomId == -1 && bondId != -1) {
        AtomId a = -1, b = -1;
        if (m_molecule.bondEndpoints(bondId, a, b)) atomId = a;
    }
    if (atomId == -1) {
        clearSelection();
        return;
    }

    QList<AtomId> allAtoms = m_molecule.atomIds();
    QSet<AtomId> validAtoms(allAtoms.begin(), allAtoms.end());
    if (!validAtoms.contains(atomId)) return; // no-op: not a real atom (sgroup-prim resolution out of scope)

    int fragIdx = m_molecule.atomFragmentIndex(atomId);
    QList<AtomId> fragAtoms = m_molecule.atomIdsInFragment(fragIdx);

    m_selection.clear();
    for (AtomId id : fragAtoms) m_selection.atoms.insert(id);

    QSet<AtomId> atomSet(fragAtoms.begin(), fragAtoms.end());
    for (BondId bid : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(bid, ea, eb)) continue;
        if (atomSet.contains(ea) && atomSet.contains(eb)) m_selection.bonds.insert(bid);
    }
}

void DocumentState::selectRing(AtomId atomId, BondId bondId) {
    if (atomId == -1 && bondId == -1) {
        clearSelection();
        return;
    }

    BondId startBondId = -1;
    if (bondId != -1) {
        startBondId = bondId; // bondId takes precedence when both given -- matches the real code
    } else {
        int minLen = -1;
        BondId bestBond = -1;
        for (BondId bid : m_molecule.bondIds()) {
            AtomId a = -1, b = -1;
            if (!m_molecule.bondEndpoints(bid, a, b)) continue;
            if (a != atomId && b != atomId) continue;
            QList<AtomId> ring = bfsShortestRing(m_molecule, bid);
            if (!ring.isEmpty() && (minLen == -1 || ring.size() < minLen)) {
                minLen = ring.size();
                bestBond = bid;
            }
        }
        startBondId = bestBond;
    }

    if (startBondId == -1) return; // no start bond found -- leave selection untouched

    QList<AtomId> ringAtoms = bfsShortestRing(m_molecule, startBondId);
    if (ringAtoms.isEmpty()) return; // not part of any ring -- leave selection untouched

    QSet<AtomId> atomSet(ringAtoms.begin(), ringAtoms.end());
    m_selection.clear();
    for (AtomId id : ringAtoms) m_selection.atoms.insert(id);
    for (BondId bid : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(bid, ea, eb)) continue;
        if (atomSet.contains(ea) && atomSet.contains(eb)) m_selection.bonds.insert(bid);
    }
}

void DocumentState::selectChain(AtomId atomId, BondId bondId) {
    if (atomId == -1 && bondId == -1) {
        clearSelection();
        return;
    }

    QList<AtomId> queue;
    QSet<AtomId> visited;
    AtomId seedA = -1, seedB = -1;

    if (atomId != -1) {
        // atomId takes precedence when both given -- matches the real code. No existence
        // check here: a bogus id is still pushed and selected as a lone 1-atom result (see
        // this method's own header comment for why that's the ported behavior, not a bug).
        queue.append(atomId);
        visited.insert(atomId);
    } else {
        AtomId a = -1, b = -1;
        if (!m_molecule.bondEndpoints(bondId, a, b)) return; // invalid bond -- leave selection untouched
        seedA = a;
        seedB = b;
        queue.append(a);
        queue.append(b);
        visited.insert(a);
        visited.insert(b);
    }

    // Ring-atom membership, computed once up-front over the WHOLE molecule -- matches the
    // real JS's own precompute exactly (same O(bonds^2) cost, not optimized away: fidelity
    // over performance for this port).
    QSet<AtomId> ringAtoms;
    for (BondId bid : m_molecule.bondIds()) {
        QList<AtomId> cycle = bfsShortestRing(m_molecule, bid);
        for (AtomId a : cycle) ringAtoms.insert(a);
    }

    QList<AtomId> resultAtoms;
    int head = 0;
    while (head < queue.size()) {
        AtomId cur = queue[head++];
        resultAtoms.append(cur);

        QString symbol = m_molecule.atomSymbol(cur);
        bool isH = (symbol == QStringLiteral("H"));
        bool isStart = (atomId != -1 && cur == atomId) || (atomId == -1 && (cur == seedA || cur == seedB));

        if (ringAtoms.contains(cur) && !isH && !isStart) {
            continue; // ring atom, not the seed: included but don't expand past it
        }

        QList<AtomId> neighbors = m_molecule.neighborAtomIds(cur);
        if (!isH) {
            int heavyCount = 0;
            for (AtomId n : neighbors) {
                if (m_molecule.atomSymbol(n) != QStringLiteral("H")) ++heavyCount;
            }
            if (heavyCount > 2) continue;             // branch point: unconditional, even at the seed
            if (heavyCount == 1 && !isStart) continue; // terminal, not the seed
            if (heavyCount == 0) continue;             // isolated
        }

        for (AtomId n : neighbors) {
            if (!visited.contains(n)) {
                visited.insert(n);
                queue.append(n);
            }
        }
    }

    QSet<AtomId> atomSet(resultAtoms.begin(), resultAtoms.end());
    m_selection.clear();
    for (AtomId id : resultAtoms) m_selection.atoms.insert(id);
    for (BondId bid : m_molecule.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!m_molecule.bondEndpoints(bid, ea, eb)) continue;
        if (atomSet.contains(ea) && atomSet.contains(eb)) m_selection.bonds.insert(bid);
    }
}

QString DocumentState::copySelection() const {
    if (m_selection.atoms.isEmpty()) return QString();
    QList<AtomId> atomIds(m_selection.atoms.begin(), m_selection.atoms.end());
    QList<BondId> bondIds(m_selection.bonds.begin(), m_selection.bonds.end());
    StringResult r = m_molecule.submoleculeMolfile(atomIds, bondIds);
    return r.success ? r.value : QString();
}

AtomId DocumentState::addAtom(const QString& symbol, double x, double y) {
    auto idBox = std::make_shared<AtomId>(-1);
    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, idBox, symbol, x, y]() { *idBox = mol.addAtom(symbol, x, y); };
    cmd.invert = [&mol, idBox]() { mol.removeAtom(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

BondId DocumentState::addBond(AtomId a, AtomId b, int order, int stereo) {
    auto idBox = std::make_shared<BondId>(-1);
    EditableMolecule& mol = m_molecule;
    EditableMolecule::Direction dir;
    switch (stereo) {
        case 1: dir = EditableMolecule::Direction::Up; break;
        case 6: dir = EditableMolecule::Direction::Down; break;
        case 4: dir = EditableMolecule::Direction::Either; break;
        default: dir = EditableMolecule::Direction::None; break;
    }
    EditCommand cmd;
    cmd.execute = [&mol, idBox, a, b, order, dir]() {
        *idBox = mol.addBond(a, b, order);
        if (dir != EditableMolecule::Direction::None) mol.setBondStereo(*idBox, dir);
    };
    cmd.invert = [&mol, idBox]() { mol.removeBond(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::deleteAtom(AtomId id) {
    // Plain-atom-and-incident-bonds case only in this sub-project (see the
    // sgroup-membership spike, Task 4). EditableMolecule::removeAtom already
    // cascades incident-bond removal internally (sub-project 1) for the
    // forward direction; the inverse must recreate the atom AND every bond
    // that existed on it, so incident-bond data is captured BEFORE removal
    // via bondEndpoints (Step 4) -- mirroring 20-edit.js's deleteAtomById
    // collecting bondsData first.
    EditableMolecule& mol = m_molecule;
    if (!mol.atomIds().contains(id)) return;

    QString label = mol.atomSymbol(id);
    double x = 0, y = 0;
    mol.atomPos(id, x, y);

    struct SavedBond { AtomId other; int order; bool idIsSource; };
    QList<SavedBond> savedBonds;
    for (BondId bid : mol.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!mol.bondEndpoints(bid, ea, eb)) continue;
        if (ea == id) savedBonds.append({eb, mol.bondOrder(bid), true});
        else if (eb == id) savedBonds.append({ea, mol.bondOrder(bid), false});
    }

    // invert() recreates the entity under a fresh id, so the id execute() removes is held
    // in a shared box that invert() refreshes -- otherwise redo after undo targets an id
    // that no longer exists and silently does nothing (same reason as deleteSelectionEntities).
    auto idBox = std::make_shared<AtomId>(id);
    EditCommand cmd;
    cmd.execute = [&mol, idBox]() { mol.removeAtom(*idBox); };
    cmd.invert = [&mol, idBox, label, x, y, savedBonds]() {
        AtomId newId = mol.addAtom(label, x, y);
        *idBox = newId;
        for (const SavedBond& sb : savedBonds) {
            if (sb.idIsSource) mol.addBond(newId, sb.other, sb.order);
            else mol.addBond(sb.other, newId, sb.order);
        }
    };
    executeCommand(std::move(cmd));
}

void DocumentState::deleteBond(BondId id) {
    EditableMolecule& mol = m_molecule;
    AtomId a = -1, b = -1;
    if (!mol.bondEndpoints(id, a, b)) return;
    int order = mol.bondOrder(id);

    EditCommand cmd;
    auto idBox = std::make_shared<BondId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeBond(*idBox); };
    cmd.invert = [&mol, idBox, a, b, order]() { *idBox = mol.addBond(a, b, order); };
    executeCommand(std::move(cmd));
}

void DocumentState::deleteSelectionEntities() {
    EditableMolecule& mol = m_molecule;

    // Step 1: disambiguate m_selection.atoms into real atoms vs. whole-pill sgroup ids.
    // Atom-existence is the gate, sgroup membership is the fallback -- AtomId/SGroupId are
    // independent, colliding counters, so checking sgroupIds() first would be unsafe.
    QSet<AtomId> realAtoms;
    QSet<SGroupId> wholePillSgroups;
    // Each call to atomIds()/sgroupIds() returns a FRESH QList -- capture into a named local
    // first rather than chaining .begin()/.end() off two separate temporaries (that mixes
    // iterators from two different QList instances, which is undefined behavior and crashed
    // reproducibly under gdb: SIGSEGV inside QSet's range constructor).
    QList<AtomId> allAtomIds = mol.atomIds();
    QList<SGroupId> allSgroupIds = mol.sgroupIds();
    QSet<AtomId> atomIdSet(allAtomIds.begin(), allAtomIds.end());
    QSet<SGroupId> sgroupIdSet(allSgroupIds.begin(), allSgroupIds.end());
    for (AtomId id : m_selection.atoms) {
        if (atomIdSet.contains(id)) {
            realAtoms.insert(id);
        } else if (sgroupIdSet.contains(id)) {
            wholePillSgroups.insert(static_cast<SGroupId>(id));
        }
        // else: stale/invalid id, silently skipped.
    }
    // Whole-pill selections expand to their member atoms too (deletion cascades to the
    // atoms in the general case, but whole-pill sgroups are handled separately below and
    // must NOT also appear in the resolved-atom set).
    QSet<AtomId> resolvedAtoms = realAtoms;

    // Step 2 covered by wholePillSgroups above.

    // Step 4: cascading incident-bond deletion.
    QSet<BondId> resolvedBonds = m_selection.bonds;
    for (BondId bid : mol.bondIds()) {
        if (resolvedBonds.contains(bid)) continue;
        AtomId ea = -1, eb = -1;
        if (!mol.bondEndpoints(bid, ea, eb)) continue;
        if (resolvedAtoms.contains(ea) || resolvedAtoms.contains(eb)) resolvedBonds.insert(bid);
    }

    // Step 5: for every resolved atom that is an individual member of some sgroup NOT
    // already covered by a whole-pill deletion, capture that sgroup's CURRENT full member
    // list (for undo's createSuperatomFromAtoms call).
    QHash<SGroupId, QList<AtomId>> affectedSgroupMembers;
    for (SGroupId sid : mol.sgroupIds()) {
        if (wholePillSgroups.contains(sid)) continue;
        QList<AtomId> members = mol.sgroupMemberAtomIds(sid);
        bool anyMemberDeleted = false;
        for (AtomId m : members) if (resolvedAtoms.contains(m)) { anyMemberDeleted = true; break; }
        if (anyMemberDeleted) affectedSgroupMembers.insert(sid, members);
    }

    // Capture bond data BEFORE deletion (endpoints + order), for undo.
    struct SavedBond { AtomId a, b; int order; };
    QList<SavedBond> savedBonds;
    for (BondId bid : resolvedBonds) {
        AtomId a = -1, b = -1;
        if (!mol.bondEndpoints(bid, a, b)) continue;
        savedBonds.append({a, b, mol.bondOrder(bid)});
    }

    // Capture atom data BEFORE deletion (label + position), for undo.
    struct SavedAtom { AtomId id; QString label; double x, y; };
    QList<SavedAtom> savedAtoms;
    for (AtomId aid : resolvedAtoms) {
        double x = 0, y = 0;
        mol.atomPos(aid, x, y);
        savedAtoms.append({aid, mol.atomSymbol(aid), x, y});
    }

    // Whole-pill sgroups: capture member atoms (all survive; no id translation needed).
    QHash<SGroupId, QList<AtomId>> wholePillMembers;
    for (SGroupId sid : wholePillSgroups) wholePillMembers.insert(sid, mol.sgroupMemberAtomIds(sid));

    // Step 6: capture selected rxnArrows/rxnPluses/multitailArrows (same accessors
    // deleteRxnArrow/deleteRxnPlus/deleteMultitailArrow already use).
    struct SavedArrow { RxnArrowId id; double x1, y1, x2, y2; QString mode, above, below; bool hasCurvature; double cx, cy; };
    QList<SavedArrow> savedArrows;
    for (RxnArrowId id : m_selection.rxnArrows) {
        SavedArrow sa; sa.id = id;
        if (!mol.rxnArrowEndpoints(id, sa.x1, sa.y1, sa.x2, sa.y2)) continue;
        sa.mode = mol.rxnArrowMode(id);
        sa.above = mol.rxnArrowConditionsAbove(id);
        sa.below = mol.rxnArrowConditionsBelow(id);
        sa.hasCurvature = mol.rxnArrowCurvature(id, sa.cx, sa.cy);
        savedArrows.append(sa);
    }
    struct SavedPlus { double x, y; };
    QList<SavedPlus> savedPluses;
    for (RxnPlusId id : m_selection.rxnPluses) {
        SavedPlus sp;
        if (!mol.rxnPlusPos(id, sp.x, sp.y)) continue;
        savedPluses.append(sp);
    }
    struct SavedMta { QList<double> pts; };
    QList<SavedMta> savedMtas;
    for (MultitailArrowId id : m_selection.multitailArrows) {
        QList<double> pts = mol.multitailArrowPoints(id);
        if (pts.isEmpty()) continue;
        savedMtas.append({pts});
    }

    // Mutable shared state (not plain captured-by-value lists): this port's own established
    // fresh-id-on-redo convention means every entity invert() recreates gets a BRAND NEW id,
    // never the original. Without updating these after each invert(), a SECOND execute()
    // (redo after undo) would try to remove ids that no longer exist -- EditableMolecule's
    // remove* methods silently return false on an unknown id, so this bug is completely
    // invisible except as "redo after undo does nothing" (found via live UI testing during
    // sub-project 7b; no existing test exercised delete -> undo -> redo specifically).
    auto bondsToRemove = std::make_shared<QList<BondId>>(resolvedBonds.values());
    auto atomsToRemove = std::make_shared<QList<AtomId>>(resolvedAtoms.values());
    auto pillsToRemove = std::make_shared<QList<SGroupId>>(wholePillSgroups.values());
    auto arrowsToRemove = std::make_shared<QList<RxnArrowId>>(m_selection.rxnArrows.begin(), m_selection.rxnArrows.end());
    auto plusesToRemove = std::make_shared<QList<RxnPlusId>>(m_selection.rxnPluses.begin(), m_selection.rxnPluses.end());
    auto mtasToRemove = std::make_shared<QList<MultitailArrowId>>(m_selection.multitailArrows.begin(), m_selection.multitailArrows.end());

    EditCommand cmd;
    cmd.execute = [&mol, bondsToRemove, atomsToRemove, pillsToRemove,
                   arrowsToRemove, plusesToRemove, mtasToRemove]() {
        // Order-independent per direct probe finding (Indigo auto-trims/auto-removes
        // affected sgroups correctly regardless of removal order) -- kept fixed for
        // readability: bonds, then atoms, then whole-pill sgroups, then non-molecular.
        for (BondId b : *bondsToRemove) mol.removeBond(b);
        for (AtomId a : *atomsToRemove) mol.removeAtom(a);
        for (SGroupId s : *pillsToRemove) mol.removeSuperatomOnly(s);
        for (RxnArrowId a : *arrowsToRemove) mol.removeRxnArrow(a);
        for (RxnPlusId p : *plusesToRemove) mol.removeRxnPlus(p);
        for (MultitailArrowId m : *mtasToRemove) mol.removeMultitailArrow(m);
    };
    cmd.invert = [&mol, savedAtoms, savedBonds, affectedSgroupMembers, wholePillMembers,
                  savedArrows, savedPluses, savedMtas,
                  bondsToRemove, atomsToRemove, pillsToRemove,
                  arrowsToRemove, plusesToRemove, mtasToRemove]() {
        // Step 1: recreate atoms fresh, building old-id -> new-id map.
        QHash<AtomId, AtomId> idMap;
        for (const SavedAtom& sa : savedAtoms) idMap.insert(sa.id, mol.addAtom(sa.label, sa.x, sa.y));

        // Refresh atomsToRemove with the ids just created, in the same order as savedAtoms,
        // so a subsequent execute() (redo) removes the CURRENT atoms, not the stale originals.
        QList<AtomId> newAtomsToRemove;
        for (const SavedAtom& sa : savedAtoms) newAtomsToRemove.append(idMap.value(sa.id));
        *atomsToRemove = newAtomsToRemove;

        // Step 2: recreate bonds via the mapped ids, capturing each fresh BondId for the
        // same reason as atomsToRemove above.
        QList<BondId> newBondsToRemove;
        for (const SavedBond& sb : savedBonds) {
            AtomId newA = idMap.value(sb.a, sb.a);
            AtomId newB = idMap.value(sb.b, sb.b);
            newBondsToRemove.append(mol.addBond(newA, newB, sb.order));
        }
        *bondsToRemove = newBondsToRemove;

        // Step 3: an individually-trimmed sgroup was NEVER removed by forward execute --
        // Indigo only auto-trimmed its member list (confirmed by probe), so the sgroup
        // (now holding fewer members) is still alive under its ORIGINAL captured id. Remove
        // that now-outdated survivor first, THEN recreate fresh with the full restored
        // member list -- creating the replacement without removing the survivor first would
        // leave TWO sgroups covering overlapping atoms (found via gdb: this exact omission
        // corrupted Indigo's sgroup state badly enough to crash on a later operation).
        // Remap each captured member id through idMap if it was deleted+recreated, or use it
        // directly if it survived untouched.
        for (auto it = affectedSgroupMembers.constBegin(); it != affectedSgroupMembers.constEnd(); ++it) {
            mol.removeSuperatomOnly(it.key());
            QList<AtomId> rebuiltMembers;
            for (AtomId originalId : it.value()) {
                rebuiltMembers.append(idMap.contains(originalId) ? idMap.value(originalId) : originalId);
            }
            mol.createSuperatomFromAtoms(rebuiltMembers);
        }

        // Step 4: recreate whole-pill sgroups -- member atoms all survived untouched, no
        // id translation needed for the MEMBERS, but the sgroup itself gets a fresh SGroupId
        // that pillsToRemove must be refreshed with (same reason as atomsToRemove above).
        QList<SGroupId> newPillsToRemove;
        for (auto it = wholePillMembers.constBegin(); it != wholePillMembers.constEnd(); ++it) {
            newPillsToRemove.append(mol.createSuperatomFromAtoms(it.value()));
        }
        *pillsToRemove = newPillsToRemove;

        // Step 5: recreate rxnArrows/rxnPluses/multitailArrows, refreshing their *ToRemove
        // lists with the fresh ids for the same reason as atomsToRemove above.
        QList<RxnArrowId> newArrowsToRemove;
        for (const SavedArrow& sa : savedArrows) {
            int newId = mol.addRxnArrow(sa.x1, sa.y1, sa.x2, sa.y2);
            mol.setRxnArrowMode(newId, sa.mode);
            mol.setRxnArrowConditions(newId, sa.above, sa.below);
            if (sa.hasCurvature) mol.setRxnArrowCurvature(newId, sa.cx, sa.cy);
            newArrowsToRemove.append(newId);
        }
        *arrowsToRemove = newArrowsToRemove;

        QList<RxnPlusId> newPlusesToRemove;
        for (const SavedPlus& sp : savedPluses) newPlusesToRemove.append(mol.addRxnPlus(sp.x, sp.y));
        *plusesToRemove = newPlusesToRemove;

        QList<MultitailArrowId> newMtasToRemove;
        for (const SavedMta& sm : savedMtas) newMtasToRemove.append(mol.addMultitailArrow(sm.pts));
        *mtasToRemove = newMtasToRemove;
    };
    executeCommand(std::move(cmd));
    clearSelection(); // matches 10-state.js's deleteSelection resetting _selection after executeCommand
}

QString DocumentState::cutSelection() {
    QString copied = copySelection();
    if (!copied.isEmpty()) deleteSelectionEntities();
    return copied;
}

void DocumentState::setStereoDescriptors(const QString& jsonMap) {
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(jsonMap.toUtf8(), &err);
    if (err.error != QJsonParseError::NoError || !doc.isObject()) return;
    QJsonObject root = doc.object();

    EditableMolecule& mol = m_molecule;
    QList<AtomId> allAtomIds = mol.atomIdsInIndigoOrder();

    for (AtomId id : mol.atomIds()) mol.setAtomStereoDescriptor(id, QString(), 0, 0);
    for (BondId id : mol.bondIds()) mol.setBondStereoCipLabel(id, QString());

    QJsonObject atomsObj = root.value(QStringLiteral("atoms")).toObject();
    for (auto it = atomsObj.constBegin(); it != atomsObj.constEnd(); ++it) {
        bool ok = false;
        int idx = it.key().toInt(&ok);
        if (!ok || idx < 0 || idx >= allAtomIds.size()) continue;
        AtomId atomId = allAtomIds[idx];
        QJsonObject entry = it.value().toObject();
        QString cipLabel = entry.value(QStringLiteral("cipLabel")).toString();
        int type = entry.value(QStringLiteral("type")).toInt(0);
        int group = entry.value(QStringLiteral("group")).toInt(0);
        mol.setAtomStereoDescriptor(atomId, cipLabel, type, group);
    }

    QJsonObject bondsObj = root.value(QStringLiteral("bonds")).toObject();
    for (auto it = bondsObj.constBegin(); it != bondsObj.constEnd(); ++it) {
        QStringList parts = it.key().split(QLatin1Char('-'));
        if (parts.size() != 2) continue;
        bool ok1 = false, ok2 = false;
        int idx1 = parts[0].toInt(&ok1);
        int idx2 = parts[1].toInt(&ok2);
        if (!ok1 || !ok2 || idx1 < 0 || idx1 >= allAtomIds.size() || idx2 < 0 || idx2 >= allAtomIds.size()) continue;
        AtomId atomId1 = allAtomIds[idx1];
        AtomId atomId2 = allAtomIds[idx2];
        BondId bondId = mol.findBond(atomId1, atomId2);
        if (bondId == -1) continue;
        QJsonObject entry = it.value().toObject();
        mol.setBondStereoCipLabel(bondId, entry.value(QStringLiteral("cipLabel")).toString());
    }
}

void DocumentState::selectSubstructureMatches(const QString& matchesJson) {
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(matchesJson.toUtf8(), &err);
    QJsonArray matches;
    if (err.error == QJsonParseError::NoError && doc.isObject()) {
        matches = doc.object().value(QStringLiteral("matches")).toArray();
    }

    EditableMolecule& mol = m_molecule;
    QList<AtomId> allAtomIds = mol.atomIdsInIndigoOrder();

    QSet<AtomId> atomIdSet;
    for (const QJsonValue& matchVal : matches) {
        QJsonArray match = matchVal.toArray();
        for (const QJsonValue& idxVal : match) {
            // QJsonValue::toInt() has no &ok overload (unlike QString/QVariant);
            // isDouble() is how a JSON value's numeric-ness is actually validated
            // here, so a null/missing/non-numeric entry doesn't silently pass as 0.
            bool ok = idxVal.isDouble();
            int idx = idxVal.toInt();
            if (ok && idx >= 1 && idx <= allAtomIds.size()) atomIdSet.insert(allAtomIds[idx - 1]);
        }
    }

    clearSelection();
    for (AtomId id : atomIdSet) addAtomToSelection(id);

    for (BondId bid : mol.bondIds()) {
        AtomId a = -1, b = -1;
        if (!mol.bondEndpoints(bid, a, b)) continue;
        if (atomIdSet.contains(a) && atomIdSet.contains(b)) addBondToSelection(bid);
    }
}

void DocumentState::setCheckIssues(const QString& jsonMap) {
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(jsonMap.toUtf8(), &err);
    if (err.error != QJsonParseError::NoError || !doc.isObject()) return;
    QJsonObject root = doc.object();

    EditableMolecule& mol = m_molecule;

    for (AtomId id : mol.atomIds()) mol.setAtomCheckWarningText(id, QString());
    for (BondId id : mol.bondIds()) mol.setBondCheckWarningText(id, QString());

    QJsonArray issues = root.value(QStringLiteral("issues")).toArray();
    for (const QJsonValue& issueVal : issues) {
        QJsonObject issue = issueVal.toObject();
        QString target = issue.value(QStringLiteral("target")).toString();
        if (target.isEmpty()) target = QStringLiteral("atom");
        QString type = issue.value(QStringLiteral("type")).toString();
        QJsonArray idsArr = issue.value(QStringLiteral("ids")).toArray();

        if (target == QStringLiteral("atom")) {
            QList<AtomId> atomOrder = mol.atomIdsInIndigoOrder();
            for (const QJsonValue& idVal : idsArr) {
                int idx = idVal.toInt(-1);
                if (idx < 0 || idx >= atomOrder.size()) continue;
                mol.setAtomCheckWarningText(atomOrder[idx], type);
            }
        } else {
            QList<BondId> bondOrder = mol.bondIdsInIndigoOrder();
            for (const QJsonValue& idVal : idsArr) {
                int idx = idVal.toInt(-1);
                if (idx < 0 || idx >= bondOrder.size()) continue;
                mol.setBondCheckWarningText(bondOrder[idx], type);
            }
        }
    }
}

void DocumentState::clearCanvas() {
    deserializeMol(QString());
}

void DocumentState::loadBenzene() {
    m_molecule.loadFrom(QString());
    m_selection.clear();
    m_history.clear();
    m_historyPointer = -1;
    resetAllDragState();
    m_molecule.addBenzeneRing(4.0, 4.0);
    markClean();
}

void DocumentState::buildBioSequenceView(const QString& sequenceText, const QString& seqType) {
    m_bioView.buildSequenceView(sequenceText, seqType);
}

void DocumentState::addBioMonomer(const QString& symbol, const QString& seqType) {
    m_bioView.addMonomer(symbol, seqType);
}

void DocumentState::deleteBioMonomer(int id) {
    m_bioView.deleteMonomer(id);
}

QList<BioMonomer> DocumentState::bioMonomers() const { return m_bioView.monomers(); }
QList<BioBond> DocumentState::bioBonds() const { return m_bioView.bonds(); }
QString DocumentState::bioSeqType() const { return m_bioView.seqType(); }

void DocumentState::changeAtomLabel(AtomId id, const QString& newLabel) {
    EditableMolecule& mol = m_molecule;
    QString oldLabel = mol.atomSymbol(id);
    if (oldLabel.isEmpty() || oldLabel == newLabel) return; // guarded, matches real JS
    EditCommand cmd;
    cmd.execute = [&mol, id, newLabel]() { mol.setAtomLabel(id, newLabel); };
    cmd.invert = [&mol, id, oldLabel]() { mol.setAtomLabel(id, oldLabel); };
    executeCommand(std::move(cmd));
}

void DocumentState::setAtomMapping(AtomId id, int mappingNumber) {
    EditableMolecule& mol = m_molecule;
    int oldAam = mol.atomAAM(id);
    if (oldAam == mappingNumber) return; // guarded, matches real JS
    EditCommand cmd;
    cmd.execute = [&mol, id, mappingNumber]() { mol.setAtomAAM(id, mappingNumber); };
    cmd.invert = [&mol, id, oldAam]() { mol.setAtomAAM(id, oldAam); };
    executeCommand(std::move(cmd));
}

void DocumentState::changeAtomCharge(AtomId id, int newCharge) {
    newCharge = std::max(-3, std::min(3, newCharge));
    EditableMolecule& mol = m_molecule;
    int oldCharge = mol.atomCharge(id);
    if (oldCharge == newCharge) return; // guarded, matches real JS
    EditCommand cmd;
    cmd.execute = [&mol, id, newCharge]() { mol.setAtomCharge(id, newCharge); };
    cmd.invert = [&mol, id, oldCharge]() { mol.setAtomCharge(id, oldCharge); };
    executeCommand(std::move(cmd));
}

void DocumentState::setAttachmentPoint(AtomId id, int order) {
    // NOT guarded -- real 20-edit.js's setAttachmentPoint always executes.
    EditableMolecule& mol = m_molecule;
    int oldOrder = mol.atomAttachmentOrder(id);
    EditCommand cmd;
    cmd.execute = [&mol, id, order]() { mol.setAtomAttachmentOrder(id, order); };
    cmd.invert = [&mol, id, oldOrder]() { mol.setAtomAttachmentOrder(id, oldOrder); };
    executeCommand(std::move(cmd));
}

void DocumentState::changeAtomIsotope(AtomId id, int isotope) {
    // NOT guarded -- real 20-edit.js's changeAtomIsotope always executes.
    isotope = std::max(0, isotope);
    EditableMolecule& mol = m_molecule;
    int oldIsotope = mol.atomIsotope(id);
    EditCommand cmd;
    cmd.execute = [&mol, id, isotope]() { mol.setAtomIsotope(id, isotope); };
    cmd.invert = [&mol, id, oldIsotope]() { mol.setAtomIsotope(id, oldIsotope); };
    executeCommand(std::move(cmd));
}

void DocumentState::changeAtomRadical(AtomId id, int radical) {
    // NOT guarded -- real 20-edit.js's changeAtomRadical always executes.
    // No clamp: unlike charge (-3..3, a real formal-charge range) or isotope
    // (>=0), radical is Indigo's own enum (INDIGO_SINGLET=101/DOUBLET=102/
    // TRIPLET=103, plus 0 for none) -- an [0,3] clamp here (MDL molfile RAD
    // code range, a different, incompatible encoding) would silently
    // truncate every real radical value before it ever reaches
    // EditableMolecule::setAtomRadical, which passes the value straight to
    // indigoSetRadical unmodified (Task 1).
    EditableMolecule& mol = m_molecule;
    int oldRadical = mol.atomRadical(id);
    EditCommand cmd;
    cmd.execute = [&mol, id, radical]() { mol.setAtomRadical(id, radical); };
    cmd.invert = [&mol, id, oldRadical]() { mol.setAtomRadical(id, oldRadical); };
    executeCommand(std::move(cmd));
}

void DocumentState::changeAtomValence(AtomId id, int valence) {
    // NOT guarded -- real 20-edit.js's changeAtomValence always executes.
    EditableMolecule& mol = m_molecule;
    int oldValence = mol.atomExplicitValence(id);
    EditCommand cmd;
    cmd.execute = [&mol, id, valence]() { mol.setAtomExplicitValence(id, valence); };
    cmd.invert = [&mol, id, oldValence]() { mol.setAtomExplicitValence(id, oldValence); };
    executeCommand(std::move(cmd));
}

DocumentState::AtomProperties DocumentState::atomProperties(AtomId id) const {
    AtomProperties props;
    props.label = m_molecule.atomSymbol(id);
    props.charge = m_molecule.atomCharge(id);
    props.isotope = m_molecule.atomIsotope(id);
    props.radical = m_molecule.atomRadical(id);
    props.explicitValence = m_molecule.atomExplicitValence(id);
    return props;
}

void DocumentState::changeBondOrder(BondId id, int newOrder) {
    EditableMolecule& mol = m_molecule;
    int oldOrder = mol.bondOrder(id);
    if (oldOrder == newOrder) return; // guarded, matches real JS's changeBondType
    EditCommand cmd;
    cmd.execute = [&mol, id, newOrder]() { mol.setBondOrderValue(id, newOrder); };
    cmd.invert = [&mol, id, oldOrder]() { mol.setBondOrderValue(id, oldOrder); };
    executeCommand(std::move(cmd));
}

void DocumentState::changeBondTypeAndStereo(BondId id, int newOrder, int newStereo) {
    EditableMolecule& mol = m_molecule;
    int oldOrder = mol.bondOrder(id);
    EditableMolecule::Direction oldDirection = mol.bondStereoDirectionEnum(id);

    EditableMolecule::Direction newDirection;
    switch (newStereo) {
        case 1: newDirection = EditableMolecule::Direction::Up; break;
        case 6: newDirection = EditableMolecule::Direction::Down; break;
        case 4: newDirection = EditableMolecule::Direction::Either; break;
        default: newDirection = EditableMolecule::Direction::None; break;
    }

    if (oldOrder == newOrder && oldDirection == newDirection) return; // matches real changeBondType's combined guard

    EditCommand cmd;
    cmd.execute = [&mol, id, newOrder, newDirection]() {
        mol.setBondOrderValue(id, newOrder);
        mol.setBondStereo(id, newDirection); // no-op if newOrder != 1, matching setBondStereo's own guard
    };
    cmd.invert = [&mol, id, oldOrder, oldDirection]() {
        mol.setBondOrderValue(id, oldOrder);
        mol.setBondStereo(id, oldDirection); // safe: order is set FIRST above, so this succeeds whenever oldOrder == 1
    };
    executeCommand(std::move(cmd));
}

void DocumentState::setBondStereo(BondId id, EditableMolecule::Direction direction) {
    EditableMolecule& mol = m_molecule;
    if (mol.bondOrder(id) != 1) return; // invalid id or non-single bond: no history entry
    EditableMolecule::Direction oldDirection = mol.bondStereoDirectionEnum(id);
    if (oldDirection == direction) return; // guarded, matches changeAtomCharge/changeBondOrder
    EditCommand cmd;
    cmd.execute = [&mol, id, direction]() { mol.setBondStereo(id, direction); };
    cmd.invert = [&mol, id, oldDirection]() { mol.setBondStereo(id, oldDirection); };
    executeCommand(std::move(cmd));
}

void DocumentState::setAtomQueryList(AtomId id, const QString& labelsCsv, bool notList) {
    EditableMolecule& mol = m_molecule;
    QString oldLabel = mol.atomSymbol(id);
    EditCommand cmd;
    cmd.execute = [&mol, id, labelsCsv, notList]() { mol.setAtomQueryList(id, labelsCsv, notList); };
    cmd.invert = [&mol, id, oldLabel]() { mol.clearAtomQueryList(id, oldLabel); };
    executeCommand(std::move(cmd));
}

// KNOWN LIMITATION: undo after clearAtomQueryList does NOT restore the
// original query list -- EditableMolecule has no "set query list directly
// from a saved numbers list" entry point bypassing the CSV/symbol parse
// (setAtomQueryList's only entry point takes a CSV string), so a correct
// invert would need that entry point added first. Out of scope for 3a;
// clearAtomQueryList's invert is documented as a no-op below rather than a
// silently-wrong "looks like it restores something" implementation. Revisit
// if a later sub-project needs this restored exactly.
void DocumentState::clearAtomQueryList(AtomId id, const QString& fallbackLabel) {
    EditableMolecule& mol = m_molecule;
    if (!mol.hasAtomQueryList(id)) return;
    EditCommand cmd;
    cmd.execute = [&mol, id, fallbackLabel]() { mol.clearAtomQueryList(id, fallbackLabel); };
    cmd.invert = []() { /* documented no-op -- see KNOWN LIMITATION comment above */ };
    executeCommand(std::move(cmd));
}

// KNOWN LIMITATION: chem-core.js's transformSelection also swaps wedge-bond
// stereo (codes 1<->6) for bonds fully inside the selection on a flip, so the
// depiction stays chemically consistent. Indigo exposes no way to do that --
// indigoBondStereo is read-only, the stereocenter-inference path was ruled out
// by sub-project 3a's Test 12, and indigoInvertStereo was ruled out by Test 14
// ("not a stereobond"; see EditableMolecule.h). Flips here apply the coordinate
// mirror ONLY; the stereo swap is an explicit, documented gap for a later
// sub-project.
void DocumentState::applyDiscreteTransform(DiscreteTransform mode) {
    EditableMolecule& mol = m_molecule;

    // Atoms only, exactly like the real transformSelection.
    struct SavedPos { AtomId id; double x, y; };
    QList<SavedPos> oldPos;
    double cx = 0, cy = 0;
    for (AtomId id : m_selection.atoms) {
        double x = 0, y = 0;
        if (!mol.atomPos(id, x, y)) continue;
        oldPos.append({id, x, y});
        cx += x; cy += y;
    }
    if (oldPos.size() < 2) return;   // matches the real ids.length < 2 guard
    cx /= oldPos.size();
    cy /= oldPos.size();

    EditCommand cmd;
    cmd.execute = [&mol, mode, oldPos, cx, cy]() {
        for (const SavedPos& p : oldPos) {
            double nx = p.x, ny = p.y;
            switch (mode) {
                case DiscreteTransform::RotateCW:  nx = cx + (p.y - cy); ny = cy - (p.x - cx); break;
                case DiscreteTransform::RotateCCW: nx = cx - (p.y - cy); ny = cy + (p.x - cx); break;
                case DiscreteTransform::FlipH:     nx = 2 * cx - p.x;    ny = p.y;             break;
                case DiscreteTransform::FlipV:     nx = p.x;             ny = 2 * cy - p.y;    break;
            }
            mol.setAtomPos(p.id, nx, ny);
        }
    };
    cmd.invert = [&mol, oldPos]() {
        for (const SavedPos& p : oldPos) mol.setAtomPos(p.id, p.x, p.y);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::rotateSelection90CW()      { applyDiscreteTransform(DiscreteTransform::RotateCW); }
void DocumentState::rotateSelection90CCW()     { applyDiscreteTransform(DiscreteTransform::RotateCCW); }
void DocumentState::flipSelectionHorizontal()  { applyDiscreteTransform(DiscreteTransform::FlipH); }
void DocumentState::flipSelectionVertical()    { applyDiscreteTransform(DiscreteTransform::FlipV); }

void DocumentState::moveImage(ImageId id, double dx, double dy) {
    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0, w = 0, h = 0;
    if (!mol.imageRect(id, x, y, w, h)) return;   // matches the real `if (!img) return`

    EditCommand cmd;
    cmd.execute = [&mol, id, dx, dy]() {
        double cx = 0, cy = 0, cw = 0, ch = 0;
        if (mol.imageRect(id, cx, cy, cw, ch)) mol.setImageRect(id, cx + dx, cy + dy, cw, ch);
    };
    cmd.invert = [&mol, id, dx, dy]() {
        double cx = 0, cy = 0, cw = 0, ch = 0;
        if (mol.imageRect(id, cx, cy, cw, ch)) mol.setImageRect(id, cx - dx, cy - dy, cw, ch);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::resizeImage(ImageId id, double scaleFactor) {
    // DELIBERATE ADDITION, not in the real JS: a zero factor would make
    // invert's 1/scaleFactor non-finite and permanently corrupt the image's
    // stored size. ImageRef holds plain doubles with no validation of its own.
    if (scaleFactor == 0.0) return;

    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0, w = 0, h = 0;
    if (!mol.imageRect(id, x, y, w, h)) return;

    EditCommand cmd;
    cmd.execute = [&mol, id, scaleFactor]() {
        double cx = 0, cy = 0, cw = 0, ch = 0;
        if (mol.imageRect(id, cx, cy, cw, ch)) mol.setImageRect(id, cx, cy, cw * scaleFactor, ch * scaleFactor);
    };
    cmd.invert = [&mol, id, scaleFactor]() {
        double cx = 0, cy = 0, cw = 0, ch = 0;
        if (mol.imageRect(id, cx, cy, cw, ch)) mol.setImageRect(id, cx, cy, cw / scaleFactor, ch / scaleFactor);
    };
    executeCommand(std::move(cmd));
}

ImageId DocumentState::addImage(const QByteArray& pngData, double cx, double cy, double halfW, double halfH) {
    if (pngData.isEmpty()) return -1;

    auto idBox = std::make_shared<ImageId>(-1);
    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, idBox, pngData, cx, cy, halfW, halfH]() {
        *idBox = mol.addImage(cx, cy, halfW, halfH, pngData);
    };
    cmd.invert = [&mol, idBox]() { mol.removeImage(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::deleteImage(ImageId id) {
    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0, w = 0, h = 0; QByteArray png;
    if (!mol.imageData(id, x, y, w, h, png)) return;

    EditCommand cmd;
    auto idBox = std::make_shared<ImageId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeImage(*idBox); };
    cmd.invert = [&mol, idBox, x, y, w, h, png]() { *idBox = mol.addImage(x, y, w, h, png); };
    executeCommand(std::move(cmd));
}

RxnArrowId DocumentState::addRxnArrow(double cx, double cy, const QString& mode) {
    auto idBox = std::make_shared<RxnArrowId>(-1);
    EditableMolecule& mol = m_molecule;
    double x2 = cx + m_drawingStyle.bondLength * 2.5;
    double y2 = cy;

    EditCommand cmd;
    cmd.execute = [&mol, idBox, cx, cy, x2, y2, mode]() {
        *idBox = mol.addRxnArrow(cx, cy, x2, y2);
        if (mode != QStringLiteral("filled-triangle")) mol.setRxnArrowMode(*idBox, mode);
    };
    cmd.invert = [&mol, idBox]() { mol.removeRxnArrow(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

RxnArrowId DocumentState::addCurvedArrow(double x1, double y1, double ctrlX, double ctrlY, double x2, double y2) {
    auto idBox = std::make_shared<RxnArrowId>(-1);
    EditableMolecule& mol = m_molecule;

    EditCommand cmd;
    cmd.execute = [&mol, idBox, x1, y1, ctrlX, ctrlY, x2, y2]() {
        *idBox = mol.addRxnArrow(x1, y1, x2, y2);
        mol.setRxnArrowMode(*idBox, QStringLiteral("curved-mechanism"));
        mol.setRxnArrowCurvature(*idBox, ctrlX, ctrlY);
    };
    cmd.invert = [&mol, idBox]() { mol.removeRxnArrow(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::setRxnArrowMode(RxnArrowId id, const QString& newMode) {
    EditableMolecule& mol = m_molecule;
    QString oldMode = mol.rxnArrowMode(id);
    if (!mol.rxnArrowIds().contains(id) || oldMode == newMode) return;

    EditCommand cmd;
    cmd.execute = [&mol, id, newMode]() { mol.setRxnArrowMode(id, newMode); };
    cmd.invert = [&mol, id, oldMode]() { mol.setRxnArrowMode(id, oldMode); };
    executeCommand(std::move(cmd));
}

void DocumentState::setRxnArrowConditions(RxnArrowId id, const QString& above, const QString& below) {
    EditableMolecule& mol = m_molecule;
    if (!mol.rxnArrowIds().contains(id)) return;
    QString oldAbove = mol.rxnArrowConditionsAbove(id);
    QString oldBelow = mol.rxnArrowConditionsBelow(id);
    if (oldAbove == above && oldBelow == below) return;

    EditCommand cmd;
    cmd.execute = [&mol, id, above, below]() { mol.setRxnArrowConditions(id, above, below); };
    cmd.invert = [&mol, id, oldAbove, oldBelow]() { mol.setRxnArrowConditions(id, oldAbove, oldBelow); };
    executeCommand(std::move(cmd));
}

void DocumentState::deleteRxnArrow(RxnArrowId id) {
    EditableMolecule& mol = m_molecule;
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    if (!mol.rxnArrowEndpoints(id, x1, y1, x2, y2)) return;
    QString mode = mol.rxnArrowMode(id);
    QString above = mol.rxnArrowConditionsAbove(id);
    QString below = mol.rxnArrowConditionsBelow(id);
    double cx = 0, cy = 0;
    bool hasCurvature = mol.rxnArrowCurvature(id, cx, cy);

    EditCommand cmd;
    auto idBox = std::make_shared<RxnArrowId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeRxnArrow(*idBox); };
    cmd.invert = [&mol, idBox, x1, y1, x2, y2, mode, above, below, hasCurvature, cx, cy]() {
        int newId = mol.addRxnArrow(x1, y1, x2, y2);
        *idBox = newId;
        mol.setRxnArrowMode(newId, mode);
        mol.setRxnArrowConditions(newId, above, below);
        if (hasCurvature) mol.setRxnArrowCurvature(newId, cx, cy);
    };
    executeCommand(std::move(cmd));
}

RxnPlusId DocumentState::addRxnPlus(double cx, double cy) {
    auto idBox = std::make_shared<RxnPlusId>(-1);
    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, idBox, cx, cy]() { *idBox = mol.addRxnPlus(cx, cy); };
    cmd.invert = [&mol, idBox]() { mol.removeRxnPlus(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::deleteRxnPlus(RxnPlusId id) {
    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0;
    if (!mol.rxnPlusPos(id, x, y)) return;

    EditCommand cmd;
    auto idBox = std::make_shared<RxnPlusId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeRxnPlus(*idBox); };
    cmd.invert = [&mol, idBox, x, y]() { *idBox = mol.addRxnPlus(x, y); };
    executeCommand(std::move(cmd));
}

MultitailArrowId DocumentState::addMultitailArrow(double cx, double cy) {
    double headX = cx + getMultitailHeadOffsetX();
    double headY = cy;

    auto idBox = std::make_shared<MultitailArrowId>(-1);
    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, idBox, headX, headY]() {
        *idBox = mol.addMultitailArrow(QList<double>{headX, headY});
    };
    cmd.invert = [&mol, idBox]() { mol.removeMultitailArrow(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::deleteMultitailArrow(MultitailArrowId id) {
    EditableMolecule& mol = m_molecule;
    QList<double> pts = mol.multitailArrowPoints(id);
    if (pts.isEmpty()) return;

    EditCommand cmd;
    auto idBox = std::make_shared<MultitailArrowId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeMultitailArrow(*idBox); };
    cmd.invert = [&mol, idBox, pts]() { *idBox = mol.addMultitailArrow(pts); };
    executeCommand(std::move(cmd));
}

void DocumentState::addMultitailArrowTail(MultitailArrowId id) {
    EditableMolecule& mol = m_molecule;
    QList<double> pts = mol.multitailArrowPoints(id);
    if (pts.size() < 2) return;

    double headX = pts[0], headY = pts[1];
    QList<double> allY;
    allY.append(headY);
    allY.append(headY + getMultitailHeight());
    for (int i = 3; i < pts.size(); i += 2) allY.append(pts[i]);
    std::sort(allY.begin(), allY.end());

    double maxGap = 0.0, gapY = allY.first();
    for (int i = 1; i < allY.size(); ++i) {
        double gap = allY[i] - allY[i - 1];
        if (gap > maxGap) { maxGap = gap; gapY = (allY[i] + allY[i - 1]) / 2.0; }
    }

    double tailX = headX - getMultitailTailInset();
    QList<double> newPts = pts;
    newPts.append(tailX);
    newPts.append(gapY);

    EditCommand cmd;
    cmd.execute = [&mol, id, newPts]() { mol.setMultitailArrowPoints(id, newPts); };
    cmd.invert = [&mol, id, pts]() { mol.setMultitailArrowPoints(id, pts); };
    executeCommand(std::move(cmd));
}

void DocumentState::setStereoFlags(const QString& type, int groupId) {
    EditableMolecule& mol = m_molecule;
    QString newType = type.isEmpty() ? QStringLiteral("abs") : type;
    int newGroupId = groupId;
    QString oldType = mol.stereoFlagsType();
    int oldGroupId = mol.stereoFlagsGroupId();
    if (oldType == newType && oldGroupId == newGroupId) return;

    EditCommand cmd;
    cmd.execute = [&mol, newType, newGroupId]() { mol.setStereoFlagsDocument(newType, newGroupId); };
    cmd.invert = [&mol, oldType, oldGroupId]() { mol.setStereoFlagsDocument(oldType, oldGroupId); };
    executeCommand(std::move(cmd));
}

bool DocumentState::addRGroup(int rgroupNumber) {
    EditableMolecule& mol = m_molecule;
    if (mol.rgroupNumbers().contains(rgroupNumber)) return false;

    EditCommand cmd;
    cmd.execute = [&mol, rgroupNumber]() { mol.addRGroupEntry(rgroupNumber); };
    cmd.invert = [&mol, rgroupNumber]() { mol.removeRGroupEntry(rgroupNumber); };
    executeCommand(std::move(cmd));
    return true;
}

bool DocumentState::deleteRGroup(int rgroupNumber) {
    EditableMolecule& mol = m_molecule;
    if (!mol.rgroupNumbers().contains(rgroupNumber)) return false;

    QString range; bool resth = false; int ifthen = 0;
    mol.rgroupLogic(rgroupNumber, range, resth, ifthen);
    QList<int> fragIds = mol.rgroupFragmentIds(rgroupNumber);

    EditCommand cmd;
    cmd.execute = [&mol, rgroupNumber]() { mol.removeRGroupEntry(rgroupNumber); };
    cmd.invert = [&mol, rgroupNumber, range, resth, ifthen, fragIds]() {
        mol.addRGroupEntry(rgroupNumber);
        mol.setRGroupLogic(rgroupNumber, range, resth, ifthen);
        for (int fid : fragIds) mol.addRGroupFragment(rgroupNumber, fid);
    };
    executeCommand(std::move(cmd));
    return true;
}

void DocumentState::setRGroupLogic(int rgroupNumber, const QString& range, bool resth, int ifthen) {
    EditableMolecule& mol = m_molecule;
    QString oldRange; bool oldResth = false; int oldIfthen = 0;
    if (!mol.rgroupLogic(rgroupNumber, oldRange, oldResth, oldIfthen)) return;
    if (oldRange == range && oldResth == resth && oldIfthen == ifthen) return;

    EditCommand cmd;
    cmd.execute = [&mol, rgroupNumber, range, resth, ifthen]() { mol.setRGroupLogic(rgroupNumber, range, resth, ifthen); };
    cmd.invert = [&mol, rgroupNumber, oldRange, oldResth, oldIfthen]() {
        mol.setRGroupLogic(rgroupNumber, oldRange, oldResth, oldIfthen);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::addRGroupMember(int rgroupNumber) {
    EditableMolecule& mol = m_molecule;
    if (!mol.rgroupNumbers().contains(rgroupNumber)) return;

    QList<int> existing = mol.rgroupFragmentIds(rgroupNumber);
    QList<int> newFragIds;
    for (AtomId aid : m_selection.atoms) {
        int fragId = mol.atomFragmentIndex(aid);
        if (fragId >= 0 && !existing.contains(fragId) && !newFragIds.contains(fragId)) {
            newFragIds.append(fragId);
        }
    }
    if (newFragIds.isEmpty()) return;

    EditCommand cmd;
    cmd.execute = [&mol, rgroupNumber, newFragIds]() {
        for (int fid : newFragIds) mol.addRGroupFragment(rgroupNumber, fid);
    };
    cmd.invert = [&mol, rgroupNumber, newFragIds]() {
        for (int fid : newFragIds) mol.removeRGroupFragment(rgroupNumber, fid);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::removeRGroupMember(int rgroupNumber, int fragId) {
    EditableMolecule& mol = m_molecule;
    if (!mol.rgroupFragmentIds(rgroupNumber).contains(fragId)) return;

    EditCommand cmd;
    cmd.execute = [&mol, rgroupNumber, fragId]() { mol.removeRGroupFragment(rgroupNumber, fragId); };
    cmd.invert = [&mol, rgroupNumber, fragId]() { mol.addRGroupFragment(rgroupNumber, fragId); };
    executeCommand(std::move(cmd));
}

TextId DocumentState::addText(const QString& plainStr, double x, double y, bool bold, bool italic) {
    if (plainStr.trimmed().isEmpty()) return -1;

    auto idBox = std::make_shared<TextId>(-1);
    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, idBox, plainStr, x, y, bold, italic]() {
        *idBox = mol.addTextAnnotation(x, y, plainStr);
        mol.setTextAnnotation(*idBox, plainStr, bold, italic);
    };
    cmd.invert = [&mol, idBox]() { mol.removeTextAnnotation(*idBox); };
    executeCommand(std::move(cmd));
    return *idBox;
}

void DocumentState::updateText(TextId id, const QString& plainStr, bool bold, bool italic) {
    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0; QString oldContent; bool oldBold = false, oldItalic = false;
    if (!mol.textAnnotationContent(id, x, y, oldContent, oldBold, oldItalic)) return;
    if (oldContent == plainStr && oldBold == bold && oldItalic == italic) return;

    EditCommand cmd;
    cmd.execute = [&mol, id, plainStr, bold, italic]() { mol.setTextAnnotation(id, plainStr, bold, italic); };
    cmd.invert = [&mol, id, oldContent, oldBold, oldItalic]() { mol.setTextAnnotation(id, oldContent, oldBold, oldItalic); };
    executeCommand(std::move(cmd));
}

void DocumentState::deleteText(TextId id) {
    EditableMolecule& mol = m_molecule;
    double x = 0, y = 0; QString content; bool bold = false, italic = false;
    if (!mol.textAnnotationContent(id, x, y, content, bold, italic)) return;

    EditCommand cmd;
    auto idBox = std::make_shared<TextId>(id); // see deleteAtom
    cmd.execute = [&mol, idBox]() { mol.removeTextAnnotation(*idBox); };
    cmd.invert = [&mol, idBox, x, y, content, bold, italic]() {
        int newId = mol.addTextAnnotation(x, y, content);
        *idBox = newId;
        mol.setTextAnnotation(newId, content, bold, italic);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::addBracketSelection() {
    if (m_selection.atoms.isEmpty() && m_selection.bonds.isEmpty()) return;

    EditableMolecule& mol = m_molecule;
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    bool has = false;

    auto accumulate = [&](AtomId aid) {
        double x = 0, y = 0;
        if (!mol.atomPos(aid, x, y)) return;
        minX = std::min(minX, x); maxX = std::max(maxX, x);
        minY = std::min(minY, y); maxY = std::max(maxY, y);
        has = true;
    };
    for (AtomId aid : m_selection.atoms) accumulate(aid);
    for (BondId bid : m_selection.bonds) {
        AtomId a = -1, b = -1;
        if (!mol.bondEndpoints(bid, a, b)) continue;
        accumulate(a);
        accumulate(b);
    }
    if (!has) return;

    const double pad = 0.8;
    double bMinX = minX - pad, bMinY = minY - pad, bMaxX = maxX + pad, bMaxY = maxY + pad;

    EditCommand cmd;
    cmd.execute = [&mol, bMinX, bMinY, bMaxX, bMaxY]() { mol.pushBracket(bMinX, bMinY, bMaxX, bMaxY); };
    cmd.invert = [&mol]() { mol.popBracket(); };
    executeCommand(std::move(cmd));
}

void DocumentState::applyMoveDelta(const QList<AtomId>& atomIds, const QList<RxnArrowId>& arrowIds,
                                   const QList<RxnPlusId>& plusIds, const QList<MultitailArrowId>& mtaIds,
                                   double dx, double dy) {
    for (AtomId id : atomIds) {
        double x = 0, y = 0;
        if (m_molecule.atomPos(id, x, y)) m_molecule.setAtomPos(id, x + dx, y + dy);
    }
    for (RxnArrowId id : arrowIds) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (m_molecule.rxnArrowEndpoints(id, x1, y1, x2, y2))
            m_molecule.setRxnArrowEndpoints(id, x1 + dx, y1 + dy, x2 + dx, y2 + dy);
    }
    for (RxnPlusId id : plusIds) {
        double x = 0, y = 0;
        if (m_molecule.rxnPlusPos(id, x, y)) m_molecule.setRxnPlusPos(id, x + dx, y + dy);
    }
    for (MultitailArrowId id : mtaIds) {
        // The real code moves only chem-core's spineTopX/spineTopY anchor; our
        // MultitailArrow stores absolute points with no anchor, so translating
        // every point is the representation-equivalent result.
        QList<double> pts = m_molecule.multitailArrowPoints(id);
        for (int i = 0; i + 1 < pts.size(); i += 2) { pts[i] += dx; pts[i + 1] += dy; }
        m_molecule.setMultitailArrowPoints(id, pts);
    }
}

void DocumentState::resetMoveDragState() {
    m_dragDeltaX = 0.0;
    m_dragDeltaY = 0.0;
    m_dragAtomIds.clear();
    m_dragArrowIds.clear();
    m_dragPlusIds.clear();
    m_dragMultitailIds.clear();
    m_dragHasOrigBBox = false;
}

void DocumentState::moveSelectionLive(double dx, double dy) {
    if (m_selection.atoms.isEmpty() && m_selection.rxnArrows.isEmpty()
        && m_selection.rxnPluses.isEmpty() && m_selection.multitailArrows.isEmpty()) return;

    // First call of the gesture: snapshot the id lists and the ATOM-ONLY bbox.
    if (m_dragDeltaX == 0.0 && m_dragDeltaY == 0.0) {
        m_dragAtomIds.clear();
        for (AtomId id : m_selection.atoms) m_dragAtomIds.append(id);
        m_dragArrowIds.clear();
        for (RxnArrowId id : m_selection.rxnArrows) m_dragArrowIds.append(id);
        m_dragPlusIds.clear();
        for (RxnPlusId id : m_selection.rxnPluses) m_dragPlusIds.append(id);
        m_dragMultitailIds.clear();
        for (MultitailArrowId id : m_selection.multitailArrows) m_dragMultitailIds.append(id);

        m_dragHasOrigBBox = false;
        double minX = 0, minY = 0, maxX = 0, maxY = 0;
        bool any = false;
        for (AtomId id : m_dragAtomIds) {
            double x = 0, y = 0;
            if (!m_molecule.atomPos(id, x, y)) continue;
            if (!any) { minX = maxX = x; minY = maxY = y; any = true; }
            else {
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            }
        }
        if (any) {
            m_dragHasOrigBBox = true;
            m_dragBBoxMinX = minX; m_dragBBoxMinY = minY;
            m_dragBBoxMaxX = maxX; m_dragBBoxMaxY = maxY;
        }
    }

    // Page-boundary hard clamp: cap the CUMULATIVE delta against the ORIGINAL
    // bbox so the whole selection stops together at the edge, like dragging a
    // window against a screen edge, rather than each item clamping
    // independently and flattening the shape. A selection wider/taller than the
    // page leaves that axis unclamped rather than fighting it.
    double newDeltaX = m_dragDeltaX + dx;
    double newDeltaY = m_dragDeltaY + dy;
    if (m_dragHasOrigBBox) {
        double minDx = kPageMinX - m_dragBBoxMinX, maxDx = kPageMaxX - m_dragBBoxMaxX;
        double minDy = kPageMinY - m_dragBBoxMinY, maxDy = kPageMaxY - m_dragBBoxMaxY;
        if (minDx <= maxDx) newDeltaX = std::max(minDx, std::min(maxDx, newDeltaX));
        if (minDy <= maxDy) newDeltaY = std::max(minDy, std::min(maxDy, newDeltaY));
    }
    double appliedDx = newDeltaX - m_dragDeltaX;
    double appliedDy = newDeltaY - m_dragDeltaY;
    m_dragDeltaX = newDeltaX;
    m_dragDeltaY = newDeltaY;

    applyMoveDelta(m_dragAtomIds, m_dragArrowIds, m_dragPlusIds, m_dragMultitailIds, appliedDx, appliedDy);
    m_dirty = true;
}

void DocumentState::commitMove() {
    if (m_dragDeltaX != 0.0 || m_dragDeltaY != 0.0) {
        double dx = m_dragDeltaX, dy = m_dragDeltaY;
        QList<AtomId> atomIds = m_dragAtomIds;
        QList<RxnArrowId> arrowIds = m_dragArrowIds;
        QList<RxnPlusId> plusIds = m_dragPlusIds;
        QList<MultitailArrowId> mtaIds = m_dragMultitailIds;
        auto isFirst = std::make_shared<bool>(true);
        DocumentState* self = this;   // safe: m_history is a member, so `this` outlives every command in it

        EditCommand cmd;
        cmd.execute = [self, isFirst, atomIds, arrowIds, plusIds, mtaIds, dx, dy]() {
            // The live calls already applied the delta, so executeCommand's own
            // synchronous first execute() must be a no-op -- exactly the real
            // commitMove's isFirst guard.
            if (*isFirst) { *isFirst = false; return; }
            self->applyMoveDelta(atomIds, arrowIds, plusIds, mtaIds, dx, dy);
        };
        cmd.invert = [self, isFirst, atomIds, arrowIds, plusIds, mtaIds, dx, dy]() {
            *isFirst = false;
            self->applyMoveDelta(atomIds, arrowIds, plusIds, mtaIds, -dx, -dy);
        };
        executeCommand(std::move(cmd));
    }
    resetMoveDragState();
}

QList<DocumentState::TransformPoint> DocumentState::snapshotSelectionPoints() const {
    QList<TransformPoint> pts;
    for (AtomId id : m_selection.atoms) {
        double x = 0, y = 0;
        if (m_molecule.atomPos(id, x, y)) pts.append({TransformPoint::Kind::Atom, id, x, y});
    }
    for (RxnArrowId id : m_selection.rxnArrows) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (!m_molecule.rxnArrowEndpoints(id, x1, y1, x2, y2)) continue;
        pts.append({TransformPoint::Kind::RxnArrowP1, id, x1, y1});
        pts.append({TransformPoint::Kind::RxnArrowP2, id, x2, y2});
    }
    for (RxnPlusId id : m_selection.rxnPluses) {
        double x = 0, y = 0;
        if (m_molecule.rxnPlusPos(id, x, y)) pts.append({TransformPoint::Kind::RxnPlus, id, x, y});
    }
    // Multitail arrows deliberately omitted -- see the header's limitation note.
    return pts;
}

void DocumentState::writeTransformPoint(EditableMolecule& mol, const TransformPoint& pt, double nx, double ny) {
    switch (pt.kind) {
        case TransformPoint::Kind::Atom:
            mol.setAtomPos(pt.id, nx, ny);
            break;
        case TransformPoint::Kind::RxnArrowP1: {
            double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
            if (mol.rxnArrowEndpoints(pt.id, x1, y1, x2, y2)) mol.setRxnArrowEndpoints(pt.id, nx, ny, x2, y2);
            break;
        }
        case TransformPoint::Kind::RxnArrowP2: {
            double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
            if (mol.rxnArrowEndpoints(pt.id, x1, y1, x2, y2)) mol.setRxnArrowEndpoints(pt.id, x1, y1, nx, ny);
            break;
        }
        case TransformPoint::Kind::RxnPlus:
            mol.setRxnPlusPos(pt.id, nx, ny);
            break;
    }
}

void DocumentState::rotateSelectionLive(double angleDelta) {
    if (m_rotateOrigPos.isEmpty()) {
        QList<TransformPoint> pts = snapshotSelectionPoints();
        if (pts.size() < 2) return;   // matches the real pts.length < 2 guard
        double cx = 0, cy = 0;
        for (const TransformPoint& p : pts) { cx += p.x; cy += p.y; }
        m_rotateOrigPos = pts;
        m_rotateCenterX = cx / pts.size();
        m_rotateCenterY = cy / pts.size();
        m_rotateTotalAngle = 0.0;
    }
    const double cx = m_rotateCenterX, cy = m_rotateCenterY;
    const double candidateTotal = m_rotateTotalAngle + angleDelta;
    const double cosA = std::cos(candidateTotal), sinA = std::sin(candidateTotal);

    // Page-boundary hard clamp: a rotation's effect on the bbox is not a simple
    // linear delta, so rather than solving for an exact boundary angle, refuse
    // to advance past the point where any point would leave the page. The drag
    // freezes there; the user can still rotate back the other way.
    for (const TransformPoint& op : m_rotateOrigPos) {
        const double odx = op.x - cx, ody = op.y - cy;
        const double nx = cx + odx * cosA - ody * sinA;
        const double ny = cy + odx * sinA + ody * cosA;
        if (nx < kPageMinX || nx > kPageMaxX || ny < kPageMinY || ny > kPageMaxY) return;
    }

    m_rotateTotalAngle = candidateTotal;
    for (const TransformPoint& p : m_rotateOrigPos) {
        const double dx = p.x - cx, dy = p.y - cy;
        writeTransformPoint(m_molecule, p, cx + dx * cosA - dy * sinA, cy + dx * sinA + dy * cosA);
    }
    m_dirty = true;
}

void DocumentState::commitRotate() {
    if (!m_rotateOrigPos.isEmpty() && m_rotateTotalAngle != 0.0) {
        QList<TransformPoint> origPos = m_rotateOrigPos;
        const double totalAngle = m_rotateTotalAngle;
        const double cx = m_rotateCenterX, cy = m_rotateCenterY;
        EditableMolecule& mol = m_molecule;
        auto isFirst = std::make_shared<bool>(true);

        // applyAngle(0) writes every point back to its snapshot position
        // (cos0=1, sin0=0 is the identity about the centre) -- that is exactly
        // how the real commitRotate implements its undo.
        auto applyAngle = [&mol, origPos, cx, cy](double angle) {
            const double c = std::cos(angle), s = std::sin(angle);
            for (const TransformPoint& p : origPos) {
                const double dx = p.x - cx, dy = p.y - cy;
                writeTransformPoint(mol, p, cx + dx * c - dy * s, cy + dx * s + dy * c);
            }
        };

        EditCommand cmd;
        cmd.execute = [isFirst, applyAngle, totalAngle]() {
            if (*isFirst) { *isFirst = false; return; }
            applyAngle(totalAngle);
        };
        cmd.invert = [isFirst, applyAngle]() {
            *isFirst = false;
            applyAngle(0.0);
        };
        executeCommand(std::move(cmd));
    }
    m_rotateOrigPos.clear();
    m_rotateTotalAngle = 0.0;
    m_rotateCenterX = 0.0;
    m_rotateCenterY = 0.0;
}

void DocumentState::scaleSelectionLive(double factor, double anchorX, double anchorY) {
    if (m_scaleOrigPos.isEmpty()) {
        QList<TransformPoint> pts = snapshotSelectionPoints();
        if (pts.size() < 2) return;   // matches the real pts.length < 2 guard
        m_scaleOrigPos = pts;
        m_scaleAnchorX = anchorX;
        m_scaleAnchorY = anchorY;
        m_scaleTotalFactor = 1.0;
    }
    // Guard against a handle dragged through its own anchor (collapse/invert).
    const double candidateFactor = std::max(0.05, factor);
    const double ax = m_scaleAnchorX, ay = m_scaleAnchorY;

    // Page-boundary hard clamp: freeze before any point would leave the page,
    // the same "check the candidate, refuse if out of bounds" approach rotate
    // uses.
    for (const TransformPoint& op : m_scaleOrigPos) {
        const double nx = ax + (op.x - ax) * candidateFactor;
        const double ny = ay + (op.y - ay) * candidateFactor;
        if (nx < kPageMinX || nx > kPageMaxX || ny < kPageMinY || ny > kPageMaxY) return;
    }

    m_scaleTotalFactor = candidateFactor;
    for (const TransformPoint& p : m_scaleOrigPos)
        writeTransformPoint(m_molecule, p, ax + (p.x - ax) * candidateFactor, ay + (p.y - ay) * candidateFactor);
    m_dirty = true;
}

void DocumentState::commitScale() {
    if (!m_scaleOrigPos.isEmpty() && m_scaleTotalFactor != 1.0) {
        QList<TransformPoint> origPos = m_scaleOrigPos;
        const double totalFactor = m_scaleTotalFactor;
        const double ax = m_scaleAnchorX, ay = m_scaleAnchorY;
        EditableMolecule& mol = m_molecule;
        auto isFirst = std::make_shared<bool>(true);

        // applyFactor(1) writes every point back to its snapshot position --
        // exactly how the real commitScale implements its undo.
        auto applyFactor = [&mol, origPos, ax, ay](double factor) {
            for (const TransformPoint& p : origPos)
                writeTransformPoint(mol, p, ax + (p.x - ax) * factor, ay + (p.y - ay) * factor);
        };

        EditCommand cmd;
        cmd.execute = [isFirst, applyFactor, totalFactor]() {
            if (*isFirst) { *isFirst = false; return; }
            applyFactor(totalFactor);
        };
        cmd.invert = [isFirst, applyFactor]() {
            *isFirst = false;
            applyFactor(1.0);
        };
        executeCommand(std::move(cmd));
    }
    m_scaleOrigPos.clear();
    m_scaleAnchorX = 0.0;
    m_scaleAnchorY = 0.0;
    m_scaleTotalFactor = 1.0;
}

void DocumentState::alignAtoms(const QString& direction) {
    if (m_selection.atoms.size() < 2) return;

    const bool horizontal = (direction == QStringLiteral("left") || direction == QStringLiteral("right")
                              || direction == QStringLiteral("centerH"));
    const bool vertical = (direction == QStringLiteral("top") || direction == QStringLiteral("bottom")
                            || direction == QStringLiteral("centerV"));
    if (!horizontal && !vertical) return;

    double minV = 0, maxV = 0, sumV = 0;
    bool any = false;
    for (AtomId id : m_selection.atoms) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        double v = horizontal ? x : y;
        if (!any) { minV = maxV = v; any = true; }
        else {
            if (v < minV) minV = v;
            if (v > maxV) maxV = v;
        }
        sumV += v;
    }
    if (!any) return;

    double target;
    if (direction == QStringLiteral("left") || direction == QStringLiteral("top")) target = minV;
    else if (direction == QStringLiteral("right") || direction == QStringLiteral("bottom")) target = maxV;
    else target = sumV / m_selection.atoms.size();

    QList<AtomId> ids;
    QList<QPointF> oldPos;
    for (AtomId id : m_selection.atoms) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        ids.append(id);
        oldPos.append(QPointF(x, y));
    }

    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, ids, oldPos, horizontal, target]() {
        for (int i = 0; i < ids.size(); ++i) {
            double x = oldPos[i].x(), y = oldPos[i].y();
            if (horizontal) x = target; else y = target;
            mol.setAtomPos(ids[i], x, y);
        }
    };
    cmd.invert = [&mol, ids, oldPos]() {
        for (int i = 0; i < ids.size(); ++i) mol.setAtomPos(ids[i], oldPos[i].x(), oldPos[i].y());
    };
    executeCommand(std::move(cmd));
}

void DocumentState::distributeAtoms(const QString& direction) {
    if (m_selection.atoms.size() < 3) return;

    const bool horizontal = (direction == QStringLiteral("horizontal"));

    QList<AtomId> ids;
    for (AtomId id : m_selection.atoms) ids.append(id);

    QList<QPointF> oldPos;
    QList<double> coordFor;
    for (AtomId id : ids) {
        double x = 0, y = 0;
        m_molecule.atomPos(id, x, y);
        oldPos.append(QPointF(x, y));
        coordFor.append(horizontal ? x : y);
    }

    QList<int> order;
    for (int i = 0; i < ids.size(); ++i) order.append(i);
    std::sort(order.begin(), order.end(), [&coordFor](int a, int b) { return coordFor[a] < coordFor[b]; });

    const double minC = coordFor[order.first()];
    const double maxC = coordFor[order.last()];
    const double step = (maxC - minC) / (order.size() - 1);

    QList<AtomId> sortedIds;
    QList<QPointF> sortedOldPos;
    QList<double> newCoord;
    for (int i = 0; i < order.size(); ++i) {
        sortedIds.append(ids[order[i]]);
        sortedOldPos.append(oldPos[order[i]]);
        newCoord.append(minC + step * i);
    }

    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, sortedIds, sortedOldPos, newCoord, horizontal]() {
        for (int i = 0; i < sortedIds.size(); ++i) {
            double x = sortedOldPos[i].x(), y = sortedOldPos[i].y();
            if (horizontal) x = newCoord[i]; else y = newCoord[i];
            mol.setAtomPos(sortedIds[i], x, y);
        }
    };
    cmd.invert = [&mol, sortedIds, sortedOldPos]() {
        for (int i = 0; i < sortedIds.size(); ++i) mol.setAtomPos(sortedIds[i], sortedOldPos[i].x(), sortedOldPos[i].y());
    };
    executeCommand(std::move(cmd));
}

double DocumentState::largestEmptyAngleAt(AtomId id) const {
    const double kPi = 3.14159265358979323846;
    double ax = 0, ay = 0;
    if (!m_molecule.atomPos(id, ax, ay)) return 0.0;

    QList<double> neighborAngles;
    for (AtomId nid : m_molecule.neighborAtomIds(id)) {
        double nx = 0, ny = 0;
        if (!m_molecule.atomPos(nid, nx, ny)) continue;
        neighborAngles.append(std::atan2(ny - ay, nx - ax));
    }

    if (neighborAngles.isEmpty()) return 0.0;
    if (neighborAngles.size() == 1) return neighborAngles.first() + 2.61799;

    std::sort(neighborAngles.begin(), neighborAngles.end());
    double maxGap = 0.0, bestAngle = 0.0;
    for (int i = 0; i < neighborAngles.size(); ++i) {
        double a1 = neighborAngles[i];
        double a2 = neighborAngles[(i + 1) % neighborAngles.size()];
        double gap = a2 - a1;
        while (gap <= 0) gap += 2.0 * kPi;
        if (gap > maxGap) {
            maxGap = gap;
            bestAngle = a1 + gap / 2.0;
        }
    }
    return bestAngle;
}

void DocumentState::layoutSelectedChain() {
    QList<AtomId> sel;
    for (AtomId id : m_selection.atoms) sel.append(id);
    if (sel.isEmpty()) return;

    // Gate 1: classify every bond touching the selection. Both ends selected -> internal
    // (builds the chain graph). Exactly one end selected -> a candidate attachment point to
    // the rest of the structure. Exactly one such attachment bond must exist.
    QHash<AtomId, QList<AtomId>> internalAdj;
    for (AtomId id : sel) internalAdj[id] = QList<AtomId>();

    struct Anchor { AtomId insideId; AtomId outsideId; };
    QList<Anchor> anchorBonds;

    for (BondId bid : m_molecule.bondIds()) {
        AtomId a = -1, b = -1;
        if (!m_molecule.bondEndpoints(bid, a, b)) continue;
        bool aSel = m_selection.atoms.contains(a);
        bool bSel = m_selection.atoms.contains(b);
        if (aSel && bSel) {
            internalAdj[a].append(b);
            internalAdj[b].append(a);
        } else if (aSel != bSel) {
            anchorBonds.append(aSel ? Anchor{a, b} : Anchor{b, a});
        }
    }

    if (anchorBonds.size() != 1) return;

    // Gate 2: no selected atom may have more than 2 internal neighbors (rules out branching).
    for (AtomId id : sel) {
        if (internalAdj.value(id).size() > 2) return;
    }

    // Gate 3: walk the chain from the atom adjacent to the anchor, refusing to revisit any
    // atom. Must visit every selected atom exactly once. NOTE: this does not reject a clean
    // cycle (every atom at or under the degree-2 cap already enforced by gate 2) -- the walk
    // simply traces once around the ring and stops once it has visited sel.size() atoms,
    // without ever needing the ring-closing bond. Only a genuinely disconnected sub-group
    // (unreachable atoms) makes chain.size() end up short of sel.size(). This matches the real
    // code's actual behavior exactly -- see spec precision note.
    Anchor anchor = anchorBonds.first();
    QList<AtomId> chain;
    chain.append(anchor.insideId);
    QSet<AtomId> visited;
    visited.insert(anchor.insideId);
    AtomId prevId = -1;
    AtomId curId = anchor.insideId;
    while (chain.size() < sel.size()) {
        AtomId nextId = -1;
        for (AtomId n : internalAdj.value(curId)) {
            if (n != prevId && !visited.contains(n)) { nextId = n; break; }
        }
        if (nextId < 0) break;
        chain.append(nextId);
        visited.insert(nextId);
        prevId = curId;
        curId = nextId;
    }
    if (chain.size() != sel.size()) return;

    // Gate 4: re-fetch the anchor's outside atom (real code re-fetches anchorAtom after the
    // walk and bails if missing). In practice this atom always exists since it came from a
    // real bond endpoint found in gate 1, but the check doubles as fetching the coordinates
    // needed for placement, so it's not dead weight to skip.
    double ax = 0, ay = 0;
    if (!m_molecule.atomPos(anchor.outsideId, ax, ay)) return;

    const double kPi = 3.14159265358979323846;
    double theta = largestEmptyAngleAt(anchor.outsideId);
    double half = kPi / 6.0;

    QList<QPointF> oldPositions;
    for (AtomId id : chain) {
        double x = 0, y = 0;
        m_molecule.atomPos(id, x, y);
        oldPositions.append(QPointF(x, y));
    }

    QList<QPointF> newPositions;
    double px = ax, py = ay;
    for (int k = 0; k < chain.size(); ++k) {
        double ang = theta + ((k % 2 == 0) ? half : -half);
        px += m_drawingStyle.bondLength * std::cos(ang);
        py += m_drawingStyle.bondLength * std::sin(ang);
        newPositions.append(QPointF(px, py));
    }

    EditableMolecule& mol = m_molecule;
    EditCommand cmd;
    cmd.execute = [&mol, chain, newPositions]() {
        for (int i = 0; i < chain.size(); ++i) mol.setAtomPos(chain[i], newPositions[i].x(), newPositions[i].y());
    };
    cmd.invert = [&mol, chain, oldPositions]() {
        for (int i = 0; i < chain.size(); ++i) mol.setAtomPos(chain[i], oldPositions[i].x(), oldPositions[i].y());
    };
    executeCommand(std::move(cmd));
}

void DocumentState::resetAllDragState() {
    resetMoveDragState();
    m_rotateOrigPos.clear();
    m_rotateTotalAngle = 0.0;
    m_rotateCenterX = 0.0;
    m_rotateCenterY = 0.0;
    m_scaleOrigPos.clear();
    m_scaleAnchorX = 0.0;
    m_scaleAnchorY = 0.0;
    m_scaleTotalFactor = 1.0;
}

void DocumentState::applyRingAlternation(const QList<AtomId>& ringAtoms, const QList<OldBondType>& oldBonds, bool aromatic) {
    int n = ringAtoms.size();
    if (n < 3) return;
    EditableMolecule& mol = m_molecule;

    auto findOldOrder = [&oldBonds](AtomId a, AtomId b) -> int {
        for (const OldBondType& ob : oldBonds) {
            if ((ob.a == a && ob.b == b) || (ob.a == b && ob.b == a)) return ob.order;
        }
        return -1;
    };

    QList<BondId> edgeBondIds;
    QList<int> edgeAnchorType;
    int anchorIdx = -1;
    for (int k = 0; k < n; ++k) {
        AtomId a1 = ringAtoms[k], a2 = ringAtoms[(k + 1) % n];
        edgeBondIds.append(mol.findBond(a1, a2));
        int anchorType = findOldOrder(a1, a2);
        edgeAnchorType.append(anchorType);
        if (anchorType != -1 && anchorIdx == -1) anchorIdx = k;
    }

    if (anchorIdx == -1) {
        // Freestanding ring: no fusion seam.
        if (!aromatic) return;
        if (n != 6) return;
        for (int k = 0; k < n; ++k) {
            if (edgeBondIds[k] >= 0) mol.setBondOrderValue(edgeBondIds[k], (k % 2 == 0) ? 2 : 1);
        }
        return;
    }

    // Fused, but explicitly non-aromatic: every new edge was already created
    // single, and the anchor edge keeps its own pre-existing type untouched
    // below -- nothing to alternate.
    if (!aromatic) return;

    // Only a clean single-seam alternation is handled; odd-sized fused rings
    // keep whichever pre-existing bonds they already had.
    if (n % 2 != 0) return;

    for (int step = 1; step < n; ++step) {
        int k = (anchorIdx + step) % n;
        if (edgeAnchorType[k] != -1) continue;   // another pre-existing edge: leave untouched
        if (edgeBondIds[k] < 0) continue;
        mol.setBondOrderValue(edgeBondIds[k], (step % 2 == 1) ? 1 : 2);
    }
}

void DocumentState::addRing(const QList<double>& coords, bool aromatic) {
    if (coords.size() < 6 || coords.size() % 2 != 0) return;   // need >= 3 points, even length

    EditableMolecule& mol = m_molecule;
    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();

    EditCommand cmd;
    cmd.execute = [this, &mol, coords, aromatic, createdAtoms, createdBonds]() {
        createdAtoms->clear();
        createdBonds->clear();

        // Bond types the ring might fuse onto, captured from the structure
        // as it stood BEFORE this ring is added.
        QList<OldBondType> oldBonds;
        for (BondId bid : mol.bondIds()) {
            AtomId a = -1, b = -1;
            if (mol.bondEndpoints(bid, a, b)) oldBonds.append({a, b, mol.bondOrder(bid)});
        }

        int n = coords.size() / 2;
        QList<AtomId> ringAtoms;
        for (int i = 0; i < n; ++i) {
            AtomId aid = mol.addAtom(QStringLiteral("C"), coords[i * 2], coords[i * 2 + 1]);
            ringAtoms.append(aid);
            createdAtoms->append(aid);
        }
        // All ring bonds start single; correct alternation is derived after
        // fusion, once we know whether any edge coincides with a pre-existing
        // bond.
        for (int i = 0; i < n; ++i) {
            BondId bid = mol.addBond(ringAtoms[i], ringAtoms[(i + 1) % n], 1);
            if (bid >= 0) createdBonds->append(bid);
        }

        EditableMolecule::MergeResult mergeResult = mol.mergeOverlappingAtoms(0.1, *createdAtoms);
        QList<AtomId> survivors;
        for (AtomId a : *createdAtoms) {
            if (!mergeResult.mergedAway.contains(a)) survivors.append(a);
        }
        *createdAtoms = survivors;
        for (BondId b : mergeResult.createdBonds) createdBonds->append(b);

        QList<AtomId> finalRingAtoms;
        for (AtomId a : ringAtoms) finalRingAtoms.append(mergeResult.mergedAway.value(a, a));

        applyRingAlternation(finalRingAtoms, oldBonds, aromatic);

        // Post-placement cleanup: normalize bond lengths for newly added atoms
        double targetBondLength = m_drawingStyle.bondLength;
        PlacementCleanup::cleanupPlacement(mol, *createdAtoms, targetBondLength);
    };
    cmd.invert = [&mol, createdAtoms, createdBonds]() {
        // Bonds first, then atoms -- avoids depending on removal cascade
        // ordering. Stale ids here (already gone via a cascade during
        // execute) fail cleanly and harmlessly, per every EditableMolecule
        // accessor's established contract.
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::addChain(double x1, double y1, double x2, double y2) {
    // Page-boundary clamp of the two INPUT endpoints, before any geometry is
    // derived -- matches the real addChain's _clampToPage(x1,y1)/_clampToPage(x2,y2).
    x1 = std::max(kPageMinX, std::min(kPageMaxX, x1));
    y1 = std::max(kPageMinY, std::min(kPageMaxY, y1));
    x2 = std::max(kPageMinX, std::min(kPageMaxX, x2));
    y2 = std::max(kPageMinY, std::min(kPageMaxY, y2));

    const double kPi = 3.14159265358979323846;
    double dx = x2 - x1, dy = y2 - y1;
    double dist = std::sqrt(dx * dx + dy * dy);
    const int kMaxChainBonds = 200;
    int nBonds = std::min(kMaxChainBonds, std::max(1, static_cast<int>(std::round(dist / m_drawingStyle.bondLength))));
    double theta = std::atan2(dy, dx);
    double half = kPi / 6.0;

    EditableMolecule& mol = m_molecule;
    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();

    double bondLength = m_drawingStyle.bondLength;
    EditCommand cmd;
    cmd.execute = [&mol, x1, y1, nBonds, theta, half, bondLength, createdAtoms, createdBonds]() {
        createdAtoms->clear();
        createdBonds->clear();

        QList<AtomId> chainAtoms;
        double px = x1, py = y1;
        for (int i = 0; i <= nBonds; ++i) {
            AtomId aid = mol.addAtom(QStringLiteral("C"), px, py);
            chainAtoms.append(aid);
            createdAtoms->append(aid);
            if (i > 0) {
                BondId bid = mol.addBond(chainAtoms[i - 1], chainAtoms[i], 1);
                if (bid >= 0) createdBonds->append(bid);
            }
            double ang = theta + ((i % 2 == 0) ? half : -half);
            px += bondLength * std::cos(ang);
            py += bondLength * std::sin(ang);
        }

        EditableMolecule::MergeResult mergeResult = mol.mergeOverlappingAtoms(0.1, *createdAtoms);
        QList<AtomId> survivors;
        for (AtomId a : *createdAtoms) {
            if (!mergeResult.mergedAway.contains(a)) survivors.append(a);
        }
        *createdAtoms = survivors;
        for (BondId b : mergeResult.createdBonds) createdBonds->append(b);

        // Post-placement cleanup: normalize bond lengths for newly added atoms
        PlacementCleanup::cleanupPlacement(mol, *createdAtoms, bondLength);
    };
    cmd.invert = [&mol, createdAtoms, createdBonds]() {
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
}

DocumentState::FunctionalGroupPlacement DocumentState::computeFunctionalGroupPlacement(
    const TemplateLibrary& lib, const QString& fgName, double cx, double cy, AtomId targetAtomId) {
    cx = std::max(kPageMinX, std::min(kPageMaxX, cx));
    cy = std::max(kPageMinY, std::min(kPageMaxY, cy));

    FunctionalGroupPlacement result;

    int fgHandle = lib.functionalGroup(fgName);
    if (fgHandle < 0) fgHandle = lib.saltOrSolvent(fgName);
    if (fgHandle < 0) fgHandle = lib.libraryTemplate(fgName);

    if (fgHandle < 0 || indigoCountAtoms(fgHandle) == 0) return result;   // valid stays false

    EditableMolecule& mol = m_molecule;

    // IMPORTANT (see Global Constraints, "sequencing consequence for Tasks
    // 4/5"): TemplateLibrary owns its own Indigo session, separate from
    // m_molecule's. Every direct indigo* call against fgHandle below MUST
    // happen before the first call into `mol` (mol.atomPos, below) --
    // EditableMolecule's own methods activate THEIR session on entry, which
    // would silently invalidate fgHandle for any later indigo* call in this
    // scope. So: find the attach point, compute the template's own bounding
    // box, AND capture the Molfile text (the only thing safe to carry across
    // the session boundary and into the command lambda, which may run again
    // much later via redo()) all first -- only then touch `mol`.

    int templateAttachIdx = -1;
    if (indigoCountSuperatoms(fgHandle) > 0) {
        int sup = indigoGetSuperatom(fgHandle, 0);
        if (sup >= 0) {
            int apIter = indigoIterateSGroupAttachmentPoints(sup);
            if (apIter >= 0) {
                int ap = indigoNext(apIter);
                if (ap > 0) { templateAttachIdx = indigoGetSGroupAttachmentPointAtomIdx(ap); indigoFree(ap); }
                indigoFree(apIter);
            }
            indigoFree(sup);
        }
    }

    double minX = 0, maxX = 0, minY = 0, maxY = 0, attachX = 0, attachY = 0;
    bool haveAttachPos = false;
    {
        bool any = false;
        int aIter = indigoIterateAtoms(fgHandle);
        if (aIter >= 0) {
            int a;
            while ((a = indigoNext(aIter)) > 0) {
                float* xyz = indigoXYZ(a);
                if (xyz) {
                    if (!any) { minX = maxX = xyz[0]; minY = maxY = xyz[1]; any = true; }
                    else {
                        if (xyz[0] < minX) minX = xyz[0];
                        if (xyz[0] > maxX) maxX = xyz[0];
                        if (xyz[1] < minY) minY = xyz[1];
                        if (xyz[1] > maxY) maxY = xyz[1];
                    }
                    if (indigoIndex(a) == templateAttachIdx) { attachX = xyz[0]; attachY = xyz[1]; haveAttachPos = true; }
                }
                indigoFree(a);
            }
            indigoFree(aIter);
        }
    }

    // Find the attach atom's overall bonded-neighbor direction within the template itself, to
    // know which direction the template's own substituent "points" by default. "First bonded
    // neighbor" is only well-defined when the attach atom has exactly one template-internal bond;
    // when it has two or more (common: 39/61 real templates), picking just the first one is
    // arbitrary (whichever bond happens to come first in molfile atom order) and was found to
    // systematically place a grafted branch on top of an already-existing bond at the target atom
    // in many real cases -- defeating the point of this feature. Instead, accumulate ALL of the
    // attach atom's internal-bond unit-vector directions and use their mean direction (average
    // unit vector, re-normalized via atan2). With exactly one bond this reduces to the original
    // single-bond direction; with 2+ bonds it points along the group's overall "bulk" direction
    // instead of an arbitrary single bond.
    // std::nullopt (not 0.0) when the attach atom has no template-internal bonds at all -- a
    // genuinely atom-only template, no known real case, but this must not silently rotate by a
    // meaningless 0-degree angle in that case. Also std::nullopt if every neighbor direction
    // happens to cancel out exactly (a symmetric case, unlikely but possible) -- translate-only
    // is the safe fallback rather than dividing by ~zero into a meaningless angle.
    std::optional<double> templateAngle;
    if (templateAttachIdx >= 0) {
        double sumX = 0, sumY = 0;
        bool anyNeighbor = false;
        int bIter = indigoIterateBonds(fgHandle);
        if (bIter >= 0) {
            int b;
            while ((b = indigoNext(bIter)) > 0) {
                int src = indigoSource(b), dst = indigoDestination(b);
                int srcIdx = indigoIndex(src), dstIdx = indigoIndex(dst);
                int neighborAtomHandle = -1;
                if (srcIdx == templateAttachIdx) neighborAtomHandle = dst;
                else if (dstIdx == templateAttachIdx) neighborAtomHandle = src;
                if (neighborAtomHandle >= 0) {
                    float* nxyz = indigoXYZ(neighborAtomHandle);
                    if (nxyz) {
                        double a = std::atan2(nxyz[1] - attachY, nxyz[0] - attachX);
                        sumX += std::cos(a);
                        sumY += std::sin(a);
                        anyNeighbor = true;
                    }
                }
                indigoFree(src); indigoFree(dst); indigoFree(b);
                // do NOT break -- must visit every bond touching the attach atom, not just the first
            }
            indigoFree(bIter);
        }
        if (anyNeighbor && std::hypot(sumX, sumY) > 1e-6) {
            templateAngle = std::atan2(sumY, sumX);
        }
    }

    // Last thing done against TemplateLibrary's session: capture the
    // structure as Molfile text. Everything from here on touches only `mol`
    // (a different session) and plain captured values.
    QString fgMolfile = lib.molfileText(fgHandle);

    bool graft = (targetAtomId >= 0) && templateAttachIdx >= 0;
    double targetX = 0, targetY = 0;
    if (graft) graft = mol.atomPos(targetAtomId, targetX, targetY);

    // Reuse BondAngleSuggester::suggestAngle (1-or-2-neighbor targets) or, when it returns
    // nullopt (0-or-3-plus-neighbor targets), suggestFallbackAngle -- build both functions'
    // atomsById/bondsList inputs from EditableMolecule directly (mol already has full graph
    // access) rather than threading V8Process's m_primitives down into DocumentState, which has
    // no existing dependency on it. This is needed whenever `graft` is true now, regardless of
    // whether templateAngle has a value -- the placement-offset angle below no longer depends on
    // templateAngle at all (only whether the TEMPLATE'S OWN SHAPE gets rotated, via useRotation,
    // still does).
    //
    // NOTE: this graft-path angle computation reads directly from EditableMolecule (the real,
    // uncollapsed molecule graph), while the separate live-preview computation in
    // AppController.cpp / V8Process::suggestFragmentAttachPoint reads from m_primitives (the
    // rendered view, where a collapsed superatom group appears as a single synthetic atom at the
    // group's centroid with no internal bonds/grandparents visible). These two views can disagree
    // for a target atom adjacent to a collapsed group, so the live drag/hover preview angle and
    // the actual committed graft angle can differ in that specific case. This is a known,
    // currently-undocumented residual gap, not something this fix solves -- just something the
    // next person touching this code needs to know about.
    std::optional<double> suggestedAngle;
    if (graft) {
        QVariantMap atomsById;
        for (AtomId id : mol.atomIds()) {
            double x = 0, y = 0;
            if (!mol.atomPos(id, x, y)) continue;
            QVariantMap a;
            a[QStringLiteral("x")] = x;
            a[QStringLiteral("y")] = y;
            atomsById[QString::number(id)] = a;
        }
        QVariantList bondsList;
        for (BondId id : mol.bondIds()) {
            AtomId a1 = -1, a2 = -1;
            if (!mol.bondEndpoints(id, a1, a2)) continue;
            QVariantMap b;
            b[QStringLiteral("begin")] = a1;
            b[QStringLiteral("end")] = a2;
            bondsList.append(b);
        }
        suggestedAngle = BondAngleSuggester::suggestAngle(targetAtomId, atomsById, bondsList);
        if (!suggestedAngle.has_value()) {
            suggestedAngle = BondAngleSuggester::suggestFallbackAngle(targetAtomId, atomsById, bondsList);
        }
    }

    // The point the template's attach atom must actually land at: one bond length (m_drawingStyle.bondLength)
    // away from the target along the placement angle above, NOT exactly on top of the target --
    // sprout semantics, not fusion. Falls back to the target's own position (attachTargetX/Y ==
    // targetX/Y) only if graft is false or somehow no angle at all is available (defensive; in
    // practice suggestedAngle always has a value whenever graft is true, since
    // suggestFallbackAngle never returns nullopt for a known atom).
    double attachTargetX = targetX, attachTargetY = targetY;
    if (graft && suggestedAngle.has_value()) {
        attachTargetX = targetX + m_drawingStyle.bondLength * std::cos(*suggestedAngle);
        attachTargetY = targetY + m_drawingStyle.bondLength * std::sin(*suggestedAngle);
    }

    // Compute the placement offset: align attach-atom-to-attachTarget if grafting,
    // else center the template's bounding box on (cx, cy).
    double dx = 0, dy = 0;
    if (graft) { dx = attachTargetX - attachX; dy = attachTargetY - attachY; }
    else { dx = cx - (minX + maxX) / 2.0; dy = cy - (minY + maxY) / 2.0; }

    bool useRotation = graft && templateAngle.has_value() && suggestedAngle.has_value() && haveAttachPos;
    double rotateBy = useRotation ? (*suggestedAngle - *templateAngle) : 0.0;

    result.valid = true;
    result.fgMolfile = fgMolfile;
    result.templateAttachIdx = templateAttachIdx;
    result.graft = graft;
    result.useRotation = useRotation;
    result.cosR = std::cos(rotateBy);
    result.sinR = std::sin(rotateBy);
    result.attachX = attachX;
    result.attachY = attachY;
    result.targetX = targetX;
    result.targetY = targetY;
    result.attachTargetX = attachTargetX;
    result.attachTargetY = attachTargetY;
    result.dx = dx;
    result.dy = dy;
    return result;
}

void DocumentState::insertFunctionalGroup(const TemplateLibrary& lib, const QString& fgName,
                                           double cx, double cy, AtomId targetAtomId, bool fullStructure) {
    EditableMolecule& mol = m_molecule;
    FunctionalGroupPlacement p = computeFunctionalGroupPlacement(lib, fgName, cx, cy, targetAtomId);

    if (!p.valid) {
        // Fallback: a single placeholder atom labeled fgName. Matches
        // computeFunctionalGroupPlacement's own condition for returning valid=false (template not
        // found or structurally empty). cx/cy re-clamped here since
        // computeFunctionalGroupPlacement clamps its own local copy but does not return it.
        double ccx = std::max(kPageMinX, std::min(kPageMaxX, cx));
        double ccy = std::max(kPageMinY, std::min(kPageMaxY, cy));
        auto idBox = std::make_shared<AtomId>(-1);
        EditCommand cmd;
        cmd.execute = [&mol, idBox, fgName, ccx, ccy]() { *idBox = mol.addAtom(fgName, ccx, ccy); };
        cmd.invert = [&mol, idBox]() { mol.removeAtom(*idBox); };
        executeCommand(std::move(cmd));
        return;
    }

    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();
    auto createdSGroups = std::make_shared<QList<SGroupId>>();

    EditCommand cmd;
    cmd.execute = [this, &mol, fgMolfile = p.fgMolfile, fgName, dx = p.dx, dy = p.dy, useRotation = p.useRotation,
                   cosR = p.cosR, sinR = p.sinR, attachX = p.attachX, attachY = p.attachY,
                   attachTargetX = p.attachTargetX, attachTargetY = p.attachTargetY,
                   templateAttachIdx = p.templateAttachIdx,
                   graft = p.graft, targetAtomId, fullStructure,
                   createdAtoms, createdBonds, createdSGroups]() {
        createdAtoms->clear();
        createdBonds->clear();
        createdSGroups->clear();

        EditableMolecule::InsertResult result = useRotation
            ? mol.insertStructure(fgMolfile, [attachX, attachY, cosR, sinR, attachTargetX, attachTargetY](double x, double y) {
                  double rx = x - attachX, ry = y - attachY;
                  double rotX = rx * cosR - ry * sinR;
                  double rotY = rx * sinR + ry * cosR;
                  return QPointF(rotX + attachTargetX, rotY + attachTargetY);
              })
            : mol.insertStructure(fgMolfile, [dx, dy](double x, double y) {
                  return QPointF(x + dx, y + dy);
              });
        *createdAtoms = result.createdAtoms;
        *createdBonds = result.createdBonds;
        *createdSGroups = result.createdSGroups;

        if (graft && templateAttachIdx >= 0) {
            AtomId newAttachAtom = result.sourceIndexToNewAtomId.value(templateAttachIdx, -1);
            if (newAttachAtom >= 0) {
                BondId newBond = mol.addBond(targetAtomId, newAttachAtom, 1);
                if (newBond >= 0) createdBonds->append(newBond);
                // newAttachAtom is NOT removed from createdAtoms -- it survives as its own real
                // atom, connected to the target by one new single bond (sprout semantics), not
                // fused away. For templates whose attach atom itself has 2+ internal bonds, this
                // is a real, visible behavior change from the old fusion model: those used to
                // fuse the target into the template's own ring/branch structure; now they always
                // become a clean single-bond substituent, with the template's own internal
                // multi-bond structure completely untouched.
            }
        }

        for (SGroupId sg : *createdSGroups) {
            mol.setSGroupExpanded(sg, fullStructure);
            mol.setSGroupLabel(sg, fgName);
        }

        // Post-placement cleanup: normalize bond lengths for newly added atoms
        // Use a local variable for bondLength to avoid capturing this
        double targetBondLength = m_drawingStyle.bondLength;
        PlacementCleanup::cleanupPlacement(mol, *createdAtoms, targetBondLength);
    };
    // createdSGroups is intentionally never swept here: EditableMolecule::removeAtom's
    // rebuildIndexTables() cascade already prunes any m_sgroupIdx entry whose underlying Indigo
    // superatom no longer exists, whenever every atom of that superatom is among createdAtoms.
    // Since the attach atom is never removed from createdAtoms at all now (it's never fused
    // away), createdAtoms already covers literally every atom the template inserted -- the same
    // "every member atom of the sgroup is in createdAtoms" precondition holds, even more simply
    // than before this change.
    cmd.invert = [&mol, createdAtoms, createdBonds]() {
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
}

PlacementResult DocumentState::getFunctionalGroupPreview(
    const TemplateLibrary& lib, const QString& fgName, double cx, double cy, AtomId targetAtomId) {
    FunctionalGroupPlacement p = computeFunctionalGroupPlacement(lib, fgName, cx, cy, targetAtomId);
    PlacementResult result;
    if (!p.valid) return result;

    // Build a throwaway EditableMolecule from the SAME fgMolfile text and the SAME transform
    // insertFunctionalGroup would commit -- mirrors insertStructureAt's own "EditableMolecule
    // tmp(sourceMolfile)" idiom rather than hand-rolling a second, parallel Indigo atom/bond
    // walk. tmp is never attached to m_molecule and is discarded when this function returns --
    // the real document is untouched.
    EditableMolecule tmp;
    auto transform = [&p](double x, double y) {
        if (p.useRotation) {
            double rx = x - p.attachX, ry = y - p.attachY;
            double rotX = rx * p.cosR - ry * p.sinR;
            double rotY = rx * p.sinR + ry * p.cosR;
            return QPointF(rotX + p.attachTargetX, rotY + p.attachTargetY);
        }
        return QPointF(x + p.dx, y + p.dy);
    };
    tmp.insertStructure(p.fgMolfile, transform);

    for (AtomId id : tmp.atomIds()) {
        double x = 0, y = 0;
        if (!tmp.atomPos(id, x, y)) continue;
        result.atoms.append({QPointF(x, y), tmp.atomSymbol(id)});
    }
    for (BondId id : tmp.bondIds()) {
        AtomId a = -1, b = -1;
        if (!tmp.bondEndpoints(id, a, b)) continue;
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (!tmp.atomPos(a, x1, y1) || !tmp.atomPos(b, x2, y2)) continue;
        result.bonds.append({QPointF(x1, y1), QPointF(x2, y2)});
    }
    // Under sprout semantics the attach atom no longer sits on top of the real target atom --
    // it's offset by one bond length. The ghost's own bonds above only ever cover the template's
    // OWN internal bonds (from tmp), so without this, the hover/drag ghost would show a floating
    // gap between the target and the group instead of the real future connection. Add the
    // not-yet-committed graft bond explicitly, using the same targetX/Y and attachTargetX/Y
    // computeFunctionalGroupPlacement already computed -- both endpoints are guaranteed
    // consistent with where tmp's own transformed attach atom actually landed, by construction
    // (same transform, same inputs).
    if (p.graft) {
        result.bonds.append({QPointF(p.targetX, p.targetY), QPointF(p.attachTargetX, p.attachTargetY)});
    }
    result.valid = true;
    return result;
}

void DocumentState::insertStructureAt(const QString& sourceMolfile, double cx, double cy) {
    cx = std::max(kPageMinX, std::min(kPageMaxX, cx));
    cy = std::max(kPageMinY, std::min(kPageMaxY, cy));

    EditableMolecule tmp(sourceMolfile);
    if (!tmp.isValid() || tmp.atomCount() == 0) return;

    double minX = 0, maxX = 0, minY = 0, maxY = 0;
    bool any = false;
    for (AtomId id : tmp.atomIds()) {
        double x = 0, y = 0;
        if (!tmp.atomPos(id, x, y)) continue;
        if (!any) { minX = maxX = x; minY = maxY = y; any = true; }
        else {
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
    }
    double dx = cx - (minX + maxX) / 2.0;
    double dy = cy - (minY + maxY) / 2.0;

    EditableMolecule& mol = m_molecule;
    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();
    auto createdSGroups = std::make_shared<QList<SGroupId>>();

    EditCommand cmd;
    cmd.execute = [&mol, sourceMolfile, dx, dy, createdAtoms, createdBonds, createdSGroups]() {
        EditableMolecule::InsertResult result = mol.insertStructure(sourceMolfile, [dx, dy](double x, double y) {
            return QPointF(x + dx, y + dy);
        });
        *createdAtoms = result.createdAtoms;
        *createdBonds = result.createdBonds;
        *createdSGroups = result.createdSGroups;
    };
    // createdSGroups is intentionally never swept here, for the same reason as
    // insertFunctionalGroup's cmd.invert (see that function's comment): removeAtom's
    // rebuildIndexTables() cascade prunes stale SGroup bookkeeping automatically once every
    // member atom is removed, and this function's plain-paste path never merges an inserted
    // atom onto a pre-existing one, so createdAtoms always covers any SGroup's full membership.
    cmd.invert = [&mol, createdAtoms, createdBonds]() {
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
}

namespace {
struct Fragment {
    QString molfile;
    double width = 0;
    double height = 0;
    double minX = 0;
    double maxX = 0;
    double minY = 0;
    double maxY = 0;
};
}

bool DocumentState::importReaction(const QString& text) {
    EditableMolecule& mol = m_molecule;
    // Every raw indigo* call below runs against mol's own session; activate it explicitly first
    // since (unlike the rest of this class) this method calls indigo* directly instead of going
    // through an EditableMolecule method that would activate its session for us. Without this,
    // the thread-local "current session" is whatever a previous, possibly-since-destroyed
    // EditableMolecule left behind (e.g. after DocumentManager::closeDocument releases one), and
    // indigoLoadReactionFromString and everything after it would run under a stale session.
    mol.ensureSession();

    int rxn = indigoLoadReactionFromString(text.toUtf8().constData());
    if (rxn < 0) return false;

    QList<Fragment> reactants, products, catalysts;
    auto collect = [&](int iterHandle, QList<Fragment>& out) {
        int comp;
        while ((comp = indigoNext(iterHandle)) > 0) {
            // Each component is laid out individually rather than calling indigoLayout(rxn) on the
            // whole reaction so that this method can compute its own custom left-to-right packing,
            // with this app's own RxnPlus/RxnArrow spacing, instead of accepting Indigo's built-in
            // reaction layout spacing.
            indigoLayout(comp);
            const char* mf = indigoMolfile(comp);
            if (!mf) { indigoFree(comp); continue; }
            Fragment f;
            f.molfile = QString::fromUtf8(mf);
            double minX = 0, maxX = 0, minY = 0, maxY = 0;
            bool any = false;
            int aIter = indigoIterateAtoms(comp);
            if (aIter >= 0) {
                int a;
                while ((a = indigoNext(aIter)) > 0) {
                    float* xyz = indigoXYZ(a);
                    if (xyz) {
                        if (!any) { minX = maxX = xyz[0]; minY = maxY = xyz[1]; any = true; }
                        else {
                            if (xyz[0] < minX) minX = xyz[0];
                            if (xyz[0] > maxX) maxX = xyz[0];
                            if (xyz[1] < minY) minY = xyz[1];
                            if (xyz[1] > maxY) maxY = xyz[1];
                        }
                    }
                    indigoFree(a);
                }
                indigoFree(aIter);
            }
            if (any) {
                // A single-atom fragment (e.g. water, "O") has a zero-width atom-coordinate
                // bounding box, but its rendered text label ("OH2") has real visual extent this
                // class can't measure directly (DocumentState is headless/UI-independent, no
                // QFontMetrics dependency) -- a defensible minimum-width floor keeps the RxnPlus
                // gap-centering fix (sub-project 7) from placing a "+" sign right on top of a
                // wide single-atom label. m_drawingStyle.bondLength * 0.6 is roughly the width a short 2-4
                // character label occupies at default zoom, matching a real single-atom-fragment
                // screenshot inspected in sub-project 7.
                f.width = std::max(maxX - minX, m_drawingStyle.bondLength * 0.6);
                f.height = maxY - minY;
                f.minX = minX;
                f.maxX = maxX;
                f.minY = minY;
                f.maxY = maxY;
            }
            out.append(f);
            indigoFree(comp);
        }
    };
    int reactantsIter = indigoIterateReactants(rxn);
    if (reactantsIter >= 0) { collect(reactantsIter, reactants); indigoFree(reactantsIter); }
    int productsIter = indigoIterateProducts(rxn);
    if (productsIter >= 0) { collect(productsIter, products); indigoFree(productsIter); }
    int catalystsIter = indigoIterateCatalysts(rxn);
    if (catalystsIter >= 0) { collect(catalystsIter, catalysts); indigoFree(catalystsIter); }
    indigoFree(rxn);

    if (reactants.isEmpty() && products.isEmpty()) return false;


    // Compute a left-to-right layout: reactants, a gap for the arrow, products. All fragments
    // share one vertical center (y=0 here; insertStructureAt's own translate-to-center pattern
    // means each fragment's own internal y-spread is what matters, not an absolute page
    // position -- the whole assembly gets clamped/placed like any other paste).
    const double gap = m_drawingStyle.bondLength;
    const double arrowGap = m_drawingStyle.bondLength * 2.5;   // matches DocumentState::addRxnArrow's own default span
    QList<double> reactantCenters, productCenters;
    double x = 0;
    for (int i = 0; i < reactants.size(); ++i) {
        if (i > 0) x += gap;
        x += reactants[i].width / 2.0;
        reactantCenters.append(x);
        x += reactants[i].width / 2.0;
    }
    double reactantBlockEnd = x;
    double arrowX1 = reactantBlockEnd + gap / 2.0;
    double arrowX2 = arrowX1 + arrowGap;
    x = arrowX2 + gap / 2.0;
    for (int i = 0; i < products.size(); ++i) {
        if (i > 0) x += gap;
        x += products[i].width / 2.0;
        productCenters.append(x);
        x += products[i].width / 2.0;
    }

    // Catalysts (Indigo's own name for what reaction-SMILES syntax calls "agents" -- the middle
    // A>agent>B component; confirmed this session that Indigo fully parses and retains this via
    // indigoIterateCatalysts, contrary to an earlier assumption that Indigo had no such API) are
    // laid out in their own left-to-right row, centered horizontally over the arrow's midpoint,
    // offset above it vertically. Negative Y renders "up" on screen in this app: ChemCanvas.qml's
    // offsetY = cy - chem.y * fNew has no sign flip, so increasing chem-space Y maps to increasing
    // (downward) screen Y -- confirmed this session, not assumed. The offset also grows by the
    // tallest catalyst fragment's own half-height (Fragment::height, already computed by collect()
    // but previously unused) -- final review found the original fixed offset only worked by
    // coincidence for the single-atom catalyst ([Pd], zero height) it was tested with; a
    // realistic multi-atom catalyst/ligand extends far enough vertically to overlap the arrow and
    // reactant/product row without this term. A later review pass found that fix was still
    // one-sided: the reactant/product row is centered on y=0 but can itself extend upward (toward
    // negative Y, toward the catalyst row) by up to its own tallest fragment's half-height, so a
    // tall reactant/product fragment paired with an ordinary-height catalyst could still overlap
    // even with the term above in place. The clearance below therefore also folds in the
    // reactant/product row's own tallest half-height, so it accounts for BOTH rows' own vertical
    // extent, not just the catalyst row's.
    double maxCatalystHalfHeight = 0.0;
    for (const Fragment& f : catalysts) maxCatalystHalfHeight = std::max(maxCatalystHalfHeight, f.height / 2.0);
    double maxRowHalfHeight = 0.0;
    for (const Fragment& f : reactants) maxRowHalfHeight = std::max(maxRowHalfHeight, f.height / 2.0);
    for (const Fragment& f : products) maxRowHalfHeight = std::max(maxRowHalfHeight, f.height / 2.0);
    const double catalystYOffset = -(m_drawingStyle.bondLength * 1.5) - maxCatalystHalfHeight - maxRowHalfHeight;
    QList<double> catalystCenters;
    {
        double catalystBlockWidth = 0;
        for (int i = 0; i < catalysts.size(); ++i) {
            if (i > 0) catalystBlockWidth += gap;
            catalystBlockWidth += catalysts[i].width;
        }
        double arrowMidX = (arrowX1 + arrowX2) / 2.0;
        double cx = arrowMidX - catalystBlockWidth / 2.0;
        for (int i = 0; i < catalysts.size(); ++i) {
            if (i > 0) cx += gap;
            cx += catalysts[i].width / 2.0;
            catalystCenters.append(cx);
            cx += catalysts[i].width / 2.0;
        }
    }

    // Recenter the whole assembly's overall bounding box onto the page origin (0,0), using the
    // same clamp-a-target-point pattern insertStructureAt uses for its own paste target. No
    // caller of importReaction supplies a target position today (V8Process::importReaction is
    // the only call site, and it always just passes the reaction text), so targetCx below always
    // evaluates to a literal 0.0 -- the std::max(kPageMinX, std::min(kPageMaxX, ...)) clamp has
    // no effect as currently called and exists only so a future caller-supplied target position
    // would be honored without restructuring this block. This genuinely RECENTERS the assembly;
    // it does NOT shrink or clamp an oversized assembly to fit within the page -- a reaction
    // assembly wider than the page will still extend past both edges after recentering, just
    // symmetrically around the origin instead of running off to one side starting at an
    // arbitrary offset. Real shrink-to-fit clamping is a separate, out-of-scope feature.
    double assemblyMinX = reactants.isEmpty() ? arrowX1 : (reactantCenters.first() - reactants.first().width / 2.0);
    double assemblyMaxX = products.isEmpty() ? arrowX2 : (productCenters.last() + products.last().width / 2.0);
    if (!catalysts.isEmpty()) {
        double catalystMinX = catalystCenters.first() - catalysts.first().width / 2.0;
        double catalystMaxX = catalystCenters.last() + catalysts.last().width / 2.0;
        assemblyMinX = std::min(assemblyMinX, catalystMinX);
        assemblyMaxX = std::max(assemblyMaxX, catalystMaxX);
    }
    // No vertical page-bounds handling is done here: the reaction layout is always vertically
    // centered on y=0 by construction, so there is nothing to recenter or clamp on that axis.
    double targetCx = std::max(kPageMinX, std::min(kPageMaxX, 0.0));
    double asmDx = targetCx - (assemblyMinX + assemblyMaxX) / 2.0;

    for (double& c : reactantCenters) c += asmDx;
    for (double& c : productCenters) c += asmDx;
    for (double& c : catalystCenters) c += asmDx;
    arrowX1 += asmDx; arrowX2 += asmDx;

    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();
    auto createdSGroups = std::make_shared<QList<SGroupId>>();
    auto createdArrows = std::make_shared<QList<RxnArrowId>>();
    auto createdPluses = std::make_shared<QList<RxnPlusId>>();

    EditCommand cmd;
    cmd.execute = [&mol, reactants, products, catalysts, reactantCenters, productCenters, catalystCenters,
                   catalystYOffset, arrowX1, arrowX2, gap,
                   createdAtoms, createdBonds, createdSGroups, createdArrows, createdPluses]() {
        auto placeFragment = [&](const Fragment& f, double cx, double cy = 0.0) {
            double centroidX = (f.minX + f.maxX) / 2.0;
            double centroidY = (f.minY + f.maxY) / 2.0;
            double dx = cx - centroidX;
            double dy = cy - centroidY;
            EditableMolecule::InsertResult r = mol.insertStructure(f.molfile, [dx, dy](double ptX, double ptY) {
                return QPointF(ptX + dx, ptY + dy);
            });
            *createdAtoms += r.createdAtoms;
            *createdBonds += r.createdBonds;
            *createdSGroups += r.createdSGroups;
        };
        for (int i = 0; i < reactants.size(); ++i) placeFragment(reactants[i], reactantCenters[i]);
        for (int i = 0; i < products.size(); ++i) placeFragment(products[i], productCenters[i]);
        for (int i = 0; i < catalysts.size(); ++i) placeFragment(catalysts[i], catalystCenters[i], catalystYOffset);
        for (int i = 1; i < reactantCenters.size(); ++i) {
            double leftEdge = reactantCenters[i-1] + reactants[i-1].width / 2.0;
            double rightEdge = reactantCenters[i] - reactants[i].width / 2.0;
            createdPluses->append(mol.addRxnPlus((leftEdge + rightEdge) / 2.0, 0));
        }
        for (int i = 1; i < productCenters.size(); ++i) {
            double leftEdge = productCenters[i-1] + products[i-1].width / 2.0;
            double rightEdge = productCenters[i] - products[i].width / 2.0;
            createdPluses->append(mol.addRxnPlus((leftEdge + rightEdge) / 2.0, 0));
        }
        createdArrows->append(mol.addRxnArrow(arrowX1, 0, arrowX2, 0));
    };
    // createdSGroups is intentionally never swept here either, for the same reason as
    // insertFunctionalGroup's and insertStructureAt's cmd.invert (see those functions'
    // comments) -- but note this function's own precondition is UNTESTED (unlike the other
    // two): a reaction-SMILES import can never produce a superatom (SMILES has no such
    // construct) so this path is moot in practice today, but a future .rxn-file import
    // containing an SUP block would need the same "every member atom is in createdAtoms"
    // property verified before relying on this comment's reasoning here specifically.
    cmd.invert = [&mol, createdArrows, createdPluses, createdBonds, createdAtoms]() {
        for (RxnArrowId a : *createdArrows) mol.removeRxnArrow(a);
        for (RxnPlusId p : *createdPluses) mol.removeRxnPlus(p);
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
    return true;
}

void DocumentState::addBondAndAtom(AtomId startId, const QString& label, double x, double y, int type, int stereo) {
    EditableMolecule& mol = m_molecule;
    if (!mol.atomIds().contains(startId)) return;

    AtomId collisionId = -1;
    for (AtomId otherId : mol.atomIds()) {
        if (otherId == startId) continue;
        double ox = 0, oy = 0;
        mol.atomPos(otherId, ox, oy);
        double d = std::sqrt((x - ox) * (x - ox) + (y - oy) * (y - oy));
        if (d < 0.3) { collisionId = otherId; break; }
    }
    if (collisionId != -1) {
        addBond(startId, collisionId, type, stereo);
        return;
    }

    double cx = std::max(kPageMinX, std::min(kPageMaxX, x));
    double cy = std::max(kPageMinY, std::min(kPageMaxY, y));

    EditableMolecule::Direction dir;
    switch (stereo) {
        case 1: dir = EditableMolecule::Direction::Up; break;
        case 6: dir = EditableMolecule::Direction::Down; break;
        case 4: dir = EditableMolecule::Direction::Either; break;
        default: dir = EditableMolecule::Direction::None; break;
    }

    auto atomIdBox = std::make_shared<AtomId>(-1);
    auto bondIdBox = std::make_shared<BondId>(-1);
    EditCommand cmd;
    cmd.execute = [&mol, atomIdBox, bondIdBox, label, cx, cy, startId, type, dir]() {
        *atomIdBox = mol.addAtom(label, cx, cy);
        *bondIdBox = mol.addBond(startId, *atomIdBox, type);
        if (dir != EditableMolecule::Direction::None) mol.setBondStereo(*bondIdBox, dir);
    };
    cmd.invert = [&mol, atomIdBox, bondIdBox]() {
        mol.removeBond(*bondIdBox);
        mol.removeAtom(*atomIdBox);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::addBondBetweenCoords(double x1, double y1, double x2, double y2, int type, int stereo) {
    EditableMolecule& mol = m_molecule;
    double cx1 = std::max(kPageMinX, std::min(kPageMaxX, x1));
    double cy1 = std::max(kPageMinY, std::min(kPageMaxY, y1));
    double cx2 = std::max(kPageMinX, std::min(kPageMaxX, x2));
    double cy2 = std::max(kPageMinY, std::min(kPageMaxY, y2));

    EditableMolecule::Direction dir;
    switch (stereo) {
        case 1: dir = EditableMolecule::Direction::Up; break;
        case 6: dir = EditableMolecule::Direction::Down; break;
        case 4: dir = EditableMolecule::Direction::Either; break;
        default: dir = EditableMolecule::Direction::None; break;
    }

    auto atom1IdBox = std::make_shared<AtomId>(-1);
    auto atom2IdBox = std::make_shared<AtomId>(-1);
    auto bondIdBox = std::make_shared<BondId>(-1);
    EditCommand cmd;
    cmd.execute = [&mol, atom1IdBox, atom2IdBox, bondIdBox, cx1, cy1, cx2, cy2, type, dir]() {
        *atom1IdBox = mol.addAtom(QStringLiteral("C"), cx1, cy1);
        *atom2IdBox = mol.addAtom(QStringLiteral("C"), cx2, cy2);
        *bondIdBox = mol.addBond(*atom1IdBox, *atom2IdBox, type);
        if (dir != EditableMolecule::Direction::None) mol.setBondStereo(*bondIdBox, dir);
    };
    cmd.invert = [&mol, atom1IdBox, atom2IdBox, bondIdBox]() {
        mol.removeBond(*bondIdBox);
        mol.removeAtom(*atom2IdBox);
        mol.removeAtom(*atom1IdBox);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::selectSingleItem(AtomId atomId, BondId bondId, RxnArrowId rxnArrowId,
                                      RxnPlusId rxnPlusId, MultitailArrowId multitailArrowId) {
    clearSelection();
    if (atomId != -1) { selectAtom(atomId); return; }
    if (bondId != -1) { selectBond(bondId); return; }
    if (rxnArrowId != -1) { selectRxnArrow(rxnArrowId); return; }
    if (rxnPlusId != -1) { selectRxnPlus(rxnPlusId); return; }
    if (multitailArrowId != -1) { selectMultitailArrow(multitailArrowId); return; }
}

std::function<QPointF(double, double)> DocumentState::makeSimilarityTransform(
        double p1x, double p1y, double p2x, double p2y, double q1x, double q1y, double q2x, double q2y) {
    double dp = std::sqrt((p2x - p1x) * (p2x - p1x) + (p2y - p1y) * (p2y - p1y));
    if (dp < 1e-9) return nullptr;
    double dq = std::sqrt((q2x - q1x) * (q2x - q1x) + (q2y - q1y) * (q2y - q1y));
    double s = dq / dp;
    double rot = std::atan2(q2y - q1y, q2x - q1x) - std::atan2(p2y - p1y, p2x - p1x);
    double cosR = std::cos(rot), sinR = std::sin(rot);
    return [=](double x, double y) {
        double rx = x - p1x, ry = y - p1y;
        return QPointF(q1x + s * (rx * cosR - ry * sinR), q1y + s * (rx * sinR + ry * cosR));
    };
}

QList<int> DocumentState::shortestRingThroughBond(int moleculeHandle, int bondIdx) {
    int bond = indigoGetBond(moleculeHandle, bondIdx);
    if (bond < 0) return {};
    int src = indigoSource(bond), dst = indigoDestination(bond);
    int start = indigoIndex(src), goal = indigoIndex(dst);
    indigoFree(src); indigoFree(dst); indigoFree(bond);

    QHash<int, QList<int>> adj;
    int bIter = indigoIterateBonds(moleculeHandle);
    if (bIter >= 0) {
        int b;
        while ((b = indigoNext(bIter)) > 0) {
            if (indigoIndex(b) != bondIdx) {
                int s = indigoSource(b), d = indigoDestination(b);
                int si = indigoIndex(s), di = indigoIndex(d);
                indigoFree(s); indigoFree(d);
                adj[si].append(di);
                adj[di].append(si);
            }
            indigoFree(b);
        }
        indigoFree(bIter);
    }

    QHash<int, int> prev;
    prev.insert(start, -2);   // -2: sentinel meaning "no predecessor" (distinct from "unvisited")
    QList<int> queue = { start };
    int qi = 0;
    while (qi < queue.size()) {
        int cur = queue[qi++];
        if (cur == goal) break;
        for (int nb : adj.value(cur)) {
            if (!prev.contains(nb)) { prev.insert(nb, cur); queue.append(nb); }
        }
    }
    if (!prev.contains(goal)) return {};

    QList<int> cycle;
    for (int a = goal; a != -2; a = prev.value(a)) cycle.append(a);
    return cycle;
}

double DocumentState::chooseEmptySide(double qax, double qay, double qbx, double qby, AtomId excludeA, AtomId excludeB,
                                       double cursorX, double cursorY) const {
    double vx = qbx - qax, vy = qby - qay;
    double scorePos = 0, scoreNeg = 0;
    for (AtomId id : m_molecule.atomIds()) {
        if (id == excludeA || id == excludeB) continue;
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        double cross = vx * (y - qay) - vy * (x - qax);
        if (cross > 0) scorePos += cross;
        else if (cross < 0) scoreNeg -= cross;
    }
    if (scorePos == 0 && scoreNeg == 0) {
        // Both sides empty: the real code breaks the tie using the live
        // cursor position at drop time. This port has no live cursor input,
        // but insertLibraryTemplateFused's own (cx, cy) parameters ARE the
        // drop-position equivalent, so they are threaded through here and
        // used exactly as the real code's cursorCross tiebreak does.
        double cursorCross = vx * (cursorY - qay) - vy * (cursorX - qax);
        return cursorCross >= 0 ? 1 : -1;
    }
    return scorePos <= scoreNeg ? 1 : -1;
}

void DocumentState::insertLibraryTemplateFused(const TemplateLibrary& lib, const QString& fgName,
                                                double cx, double cy, BondId targetBondId) {
    int fgHandle = lib.libraryTemplate(fgName);
    int bondIdx = lib.libraryTemplateFusionBondIdx(fgName);
    // Capture as Molfile text immediately -- nothing below ever touches
    // fgHandle (TemplateLibrary's own session) again. See Global Constraints,
    // "sequencing consequence for Tasks 4/5": every subsequent indigo* call
    // on the template's structure uses a LOCAL reparse (localFg, below) made
    // in m_molecule's own session instead, since m_molecule.bondEndpoints
    // (next line) already switches the active session away from
    // TemplateLibrary's.
    QString fgMolfile = (fgHandle >= 0) ? lib.molfileText(fgHandle) : QString();

    AtomId ta = -1, tb = -1;
    bool haveTargetBond = m_molecule.bondEndpoints(targetBondId, ta, tb);   // activates m_molecule's session

    int localFg = -1;
    if (!fgMolfile.isEmpty()) localFg = indigoLoadMoleculeFromString(fgMolfile.toUtf8().constData());

    QList<int> ringCycle;
    if (localFg >= 0 && bondIdx >= 0) ringCycle = shortestRingThroughBond(localFg, bondIdx);

    if (fgHandle < 0 || bondIdx < 0 || !haveTargetBond || ringCycle.isEmpty()) {
        if (localFg >= 0) indigoFree(localFg);
        insertFunctionalGroup(lib, fgName, cx, cy);
        return;
    }

    int fusionBond = indigoGetBond(localFg, bondIdx);
    int paH = indigoSource(fusionBond), pbH = indigoDestination(fusionBond);
    // indigoXYZ() returns a pointer into a buffer Indigo reuses across
    // calls (confirmed by this exact bug: reading both paXY and pbXY AFTER
    // making both calls returned pbH's coordinates for both, since the
    // second call overwrote the first's buffer) -- each result must be
    // copied out into plain doubles immediately, before the next indigoXYZ
    // call, matching the read-immediately pattern every other indigoXYZ use
    // in this codebase already follows (e.g. the bounding-box loops above).
    float* paXY = indigoXYZ(paH);
    double pax = paXY[0], pay = paXY[1];
    float* pbXY = indigoXYZ(pbH);
    double pbx = pbXY[0], pby = pbXY[1];
    indigoFree(paH); indigoFree(pbH); indigoFree(fusionBond);
    // Note: the fusion bond's own two atom indices are NOT captured
    // separately here -- shortestRingThroughBond's returned ringCycle
    // already includes both of the fusion bond's endpoints as part of the
    // cycle (it walks from one endpoint to the other through the rest of the
    // ring), so the later sourceIndexToNewAtomId lookup over ringCycle covers
    // them without a separate paIdx/pbIdx variable.

    double qax = 0, qay = 0, qbx = 0, qby = 0;
    m_molecule.atomPos(ta, qax, qay);
    m_molecule.atomPos(tb, qbx, qby);

    // Which side of the fusion bond does the template's own ring mass sit on?
    double ccx = 0, ccy = 0;
    for (int aidx : ringCycle) {
        int a = indigoGetAtom(localFg, aidx);
        float* xyz = indigoXYZ(a);
        ccx += xyz[0]; ccy += xyz[1];
        indigoFree(a);
    }
    ccx /= ringCycle.size(); ccy /= ringCycle.size();
    double templSide = ((pbx - pax) * (ccy - pay) - (pby - pay) * (ccx - pax)) >= 0 ? 1 : -1;
    double targetSide = chooseEmptySide(qax, qay, qbx, qby, ta, tb, cx, cy);

    std::function<QPointF(double, double)> transform;
    if (templSide == targetSide) transform = makeSimilarityTransform(pax, pay, pbx, pby, qax, qay, qbx, qby);
    else transform = makeSimilarityTransform(pax, pay, pbx, pby, qbx, qby, qax, qay);
    if (!transform) {
        indigoFree(localFg);
        insertFunctionalGroup(lib, fgName, cx, cy);
        return;
    }

    // Derive aromatic from the template's own authored bond types around the
    // ring cycle: any non-single bond -> true; all single -> false. Collect
    // every ring-cycle edge's bond ONCE up front (single pass over the
    // template's bonds), rather than re-scanning all bonds per ring edge.
    bool aromatic = false;
    {
        int n = ringCycle.size();
        int bIter = indigoIterateBonds(localFg);
        if (bIter >= 0) {
            int b;
            while ((b = indigoNext(bIter)) > 0) {
                int s = indigoSource(b), d = indigoDestination(b);
                int si = indigoIndex(s), di = indigoIndex(d);
                indigoFree(s); indigoFree(d);
                for (int k = 0; k < n; ++k) {
                    int a1 = ringCycle[k], a2 = ringCycle[(k + 1) % n];
                    if ((si == a1 && di == a2) || (si == a2 && di == a1)) {
                        if (indigoBondOrder(b) != 1) aromatic = true;
                        break;
                    }
                }
                indigoFree(b);
            }
            indigoFree(bIter);
        }
    }

    // Done with the local working copy -- insertStructure (inside
    // cmd.execute, below) reparses fgMolfile independently, in whatever
    // session is active when the command actually runs (its own first line
    // activates m_molecule's session regardless).
    indigoFree(localFg);

    EditableMolecule& mol = m_molecule;
    auto createdAtoms = std::make_shared<QList<AtomId>>();
    auto createdBonds = std::make_shared<QList<BondId>>();

    EditCommand cmd;
    cmd.execute = [this, &mol, fgMolfile, transform, ringCycle, aromatic, createdAtoms, createdBonds]() {
        createdAtoms->clear();
        createdBonds->clear();

        QList<OldBondType> oldBonds;
        for (BondId bid : mol.bondIds()) {
            AtomId a = -1, b = -1;
            if (mol.bondEndpoints(bid, a, b)) oldBonds.append({a, b, mol.bondOrder(bid)});
        }

        EditableMolecule::InsertResult result = mol.insertStructure(fgMolfile, transform);
        *createdAtoms = result.createdAtoms;
        *createdBonds = result.createdBonds;

        EditableMolecule::MergeResult mergeResult = mol.mergeOverlappingAtoms(0.1, *createdAtoms);
        QList<AtomId> survivors;
        for (AtomId a : *createdAtoms) {
            if (!mergeResult.mergedAway.contains(a)) survivors.append(a);
        }
        *createdAtoms = survivors;
        for (BondId b : mergeResult.createdBonds) createdBonds->append(b);

        QList<AtomId> finalRingAtoms;
        for (int idx : ringCycle) {
            AtomId aid = result.sourceIndexToNewAtomId.value(idx, -1);
            if (aid < 0) continue;
            finalRingAtoms.append(mergeResult.mergedAway.value(aid, aid));
        }
        applyRingAlternation(finalRingAtoms, oldBonds, aromatic);
    };
    cmd.invert = [&mol, createdAtoms, createdBonds]() {
        for (BondId b : *createdBonds) mol.removeBond(b);
        for (AtomId a : *createdAtoms) mol.removeAtom(a);
    };
    executeCommand(std::move(cmd));
}

void DocumentState::toggleSgroupExpanded(SGroupId id) {
    EditableMolecule& mol = m_molecule;
    AtomId probe = -1;
    if (!mol.superatomAttachAtom(id, probe)) return;   // no such sgroup: no-op
    bool wasExpanded = mol.sgroupExpanded(id);

    EditCommand cmd;
    cmd.execute = [&mol, id, wasExpanded]() { mol.setSGroupExpanded(id, !wasExpanded); };
    cmd.invert = [&mol, id, wasExpanded]() { mol.setSGroupExpanded(id, wasExpanded); };
    executeCommand(std::move(cmd));
}

void DocumentState::renameSGroup(SGroupId id, const QString& newLabel) {
    EditableMolecule& mol = m_molecule;
    AtomId probe = -1;
    if (!mol.superatomAttachAtom(id, probe)) return;   // no such sgroup: no-op
    QString oldLabel = mol.sgroupLabel(id);
    if (oldLabel == newLabel) return;

    EditCommand cmd;
    cmd.execute = [&mol, id, newLabel]() { mol.setSGroupLabel(id, newLabel); };
    cmd.invert = [&mol, id, oldLabel]() { mol.setSGroupLabel(id, oldLabel); };
    executeCommand(std::move(cmd));
}

void DocumentState::deserializeMol(const QString& data, bool centerOnPage) {
    EditableMolecule& mol = m_molecule;
    auto before = std::make_shared<MoleculeSnapshot>(mol.snapshot());
    // centerOnPage doubles as "is this a same-document write-back" (layout/aromatize/
    // normalize results, per this function's own established convention) rather than a
    // genuinely new document being opened -- see loadFrom's preserveEzIntent parameter.
    if (!mol.loadFrom(data, centerOnPage)) return;

    EditCommand cmd;
    cmd.execute = [this, &mol, data, centerOnPage]() {
        mol.loadFrom(data, centerOnPage);
        m_selection.clear();
        if (centerOnPage) centerMoleculeOnOrigin();
    };
    cmd.invert = [&mol, before]() { mol.restore(*before); };
    executeCommand(std::move(cmd));
}

void DocumentState::centerMoleculeOnOrigin() {
    double minX = 0, minY = 0, maxX = 0, maxY = 0;
    bool first = true;
    for (AtomId id : m_molecule.atomIds()) {
        double x = 0, y = 0;
        if (!m_molecule.atomPos(id, x, y)) continue;
        if (first) { minX = maxX = x; minY = maxY = y; first = false; }
        else {
            minX = std::min(minX, x); maxX = std::max(maxX, x);
            minY = std::min(minY, y); maxY = std::max(maxY, y);
        }
    }
    if (first) return;
    double cdx = -(minX + maxX) / 2.0;
    double cdy = -(minY + maxY) / 2.0;
    for (AtomId id : m_molecule.atomIds()) {
        double x = 0, y = 0;
        if (m_molecule.atomPos(id, x, y)) m_molecule.setAtomPos(id, x + cdx, y + cdy);
    }
    for (RxnArrowId id : m_molecule.rxnArrowIds()) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (m_molecule.rxnArrowEndpoints(id, x1, y1, x2, y2))
            m_molecule.setRxnArrowEndpoints(id, x1 + cdx, y1 + cdy, x2 + cdx, y2 + cdy);
    }
    for (RxnPlusId id : m_molecule.rxnPlusIds()) {
        double x = 0, y = 0;
        if (m_molecule.rxnPlusPos(id, x, y)) m_molecule.setRxnPlusPos(id, x + cdx, y + cdy);
    }
}

void DocumentState::setMoleculeName(const QString& name) {
    EditableMolecule& mol = m_molecule;
    QString oldName = mol.name();
    if (oldName == name) return;

    EditCommand cmd;
    cmd.execute = [&mol, name]() { mol.setName(name); };
    cmd.invert = [&mol, oldName]() { mol.setName(oldName); };
    executeCommand(std::move(cmd));
}


