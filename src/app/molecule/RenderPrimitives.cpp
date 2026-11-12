// src/app/molecule/RenderPrimitives.cpp
#include "RenderPrimitives.h"
#include "ElementData.h"
#include <QRegularExpression>
#include <QStringList>
#include <QSet>
#include <QHash>
#include <cmath>
#include <algorithm>

namespace {

// Ported from chem-core.js:12-18's isRGroupLabel. A cosmetic per-atom label check, NOT the
// deferred r-group data model (fragment membership, range/resth/ifthen policy).
bool isRGroupLabel(const QString& label) {
    static const QRegularExpression kPattern(QStringLiteral("^R[1-8]$"));
    return kPattern.match(label).hasMatch();
}

// Ported from chem-core.js:9082-9103's isCorrectStereoCenter -- ~20 lines of pure
// neighbor-count/implicit-H-parity logic with no Indigo dependency. Takes plain values already
// available via EditableMolecule's accessors, rather than re-deriving anything from a raw handle.
//
// endOtherNeighborCount is the neighbor-count of end's OTHER neighbor (the one that isn't begin),
// only meaningful (and only read) when endNeighborCount == 2; pass 0 when not applicable.
bool isCorrectStereoCenter(int bondStereo, int beginNeighborCount, int endNeighborCount,
                           int beginImplicitH, int endOtherNeighborCount) {
    if (bondStereo <= 0) return false;
    if (endNeighborCount == 1 && beginNeighborCount == 2 && beginImplicitH % 2 == 0) return false;
    if (endNeighborCount == 2 && beginNeighborCount == 2 && beginImplicitH % 2 == 0
        && endOtherNeighborCount == 1) return false;
    if (beginNeighborCount == 1) return false;
    return true;
}

// Display-string mapping over (stereocenterType, stereocenterGroup), matching chem-core's own
// "abs"/"&N"/"orN" convention. INDIGO_ABS=1, INDIGO_OR=2, INDIGO_AND=3, INDIGO_EITHER=4
// (indigo.h:435-438).
QString stereoLabelFor(int type, int group) {
    switch (type) {
        case 1: return QStringLiteral("abs");
        case 3: return QStringLiteral("&") + QString::number(group);
        case 2: return QStringLiteral("or") + QString::number(group);
        default: return QString();   // INDIGO_EITHER or not a stereocenter
    }
}

// Maps Indigo's CIPDesc enum (NONE=0, UNKNOWN=1, s=2, r=3, S=4, R=5, E=6, Z=7 --
// molecule_cip_calculator.h:38-48) to its display string.
QString cipDescriptorLabel(int cip) {
    switch (cip) {
        case 2: return QStringLiteral("s");
        case 3: return QStringLiteral("r");
        case 4: return QStringLiteral("S");
        case 5: return QStringLiteral("R");
        case 6: return QStringLiteral("E");
        case 7: return QStringLiteral("Z");
        default: return QString();   // NONE or UNKNOWN
    }
}

} // namespace

RenderPrimitives RenderPrimitiveBuilder::build(const EditableMolecule& mol, bool showExplicitH) {
    RenderPrimitives result;

    // ---- Sgroup contraction: build the atom->sgroup map and synthetic SgroupPrim list first ----
    // sgroupIds() enumerates every currently-existing sgroup; sgroupMemberAtomIds() gives each
    // one's FULL member-atom set (not just its attach atom), so contraction hides every member
    // atom exactly like chem-core's real sg.atoms-based contraction. A HASH (atom -> owning
    // sgroup id), not just a boolean set, is needed so Task 6's bond-primitive step can
    // distinguish "both endpoints hidden by the SAME sgroup" (skip the bond entirely) from
    // "both endpoints hidden by DIFFERENT sgroups" or "exactly one endpoint hidden" (substitute
    // the hidden endpoint(s) with their sgroup's synthetic id, matching buildRenderPrimitives'
    // own beginSgId/endSgId substitution exactly).
    QHash<AtomId, SGroupId> atomToSgroup;
    for (SGroupId sgId : mol.sgroupIds()) {
        if (mol.sgroupExpanded(sgId)) continue;   // expanded: render its atoms normally, no contraction

        AtomId attach = -1;
        mol.superatomAttachAtom(sgId, attach);   // -1 if none resolves; position falls back below

        SgroupPrim prim;
        prim.id = sgId;
        prim.attachAtomId = attach;
        prim.label = mol.sgroupLabel(sgId);   // real label for insertFunctionalGroup-created
                                               // groups (sub-project 25); still empty for a
                                               // sgroup discovered by loading/parsing an
                                               // arbitrary molfile Indigo can't name for us --
                                               // pre-existing, unrelated limitation.

        QList<AtomId> members = mol.sgroupMemberAtomIds(sgId);
        double px = 0, py = 0;
        if (attach >= 0 && mol.atomPos(attach, px, py)) {
            prim.x = px; prim.y = py;
        } else if (!members.isEmpty()) {
            // No attachment point resolved: fall back to the member centroid, matching
            // buildRenderPrimitives' own attachAtomId-else-centroid fallback.
            double sx = 0, sy = 0; int n = 0;
            for (AtomId aid : members) {
                double ax = 0, ay = 0;
                if (mol.atomPos(aid, ax, ay)) { sx += ax; sy += ay; ++n; }
            }
            if (n > 0) { prim.x = sx / n; prim.y = sy / n; }
        }
        result.sgroups.append(prim);

        AtomPrim synthetic;
        synthetic.id = sgId;
        synthetic.x = prim.x; synthetic.y = prim.y;
        synthetic.label = prim.label;
        synthetic.isSgroup = true;
        result.atoms.append(synthetic);

        for (AtomId aid : members) atomToSgroup.insert(aid, sgId);
    }

    // ---- Atom primitives ----
    double minX = 0, maxX = 0, minY = 0, maxY = 0;
    bool haveBBox = false;
    for (const AtomPrim& sg : result.atoms) {
        if (!haveBBox) { minX = maxX = sg.x; minY = maxY = sg.y; haveBBox = true; }
        else {
            if (sg.x < minX) minX = sg.x; if (sg.x > maxX) maxX = sg.x;
            if (sg.y < minY) minY = sg.y; if (sg.y > maxY) maxY = sg.y;
        }
    }

    // One CIP computation for the whole molecule; asking per atom re-ran it for every
    // atom on every refresh (minutes for a large, stereo-rich structure).
    const QHash<AtomId, int> atomCip = mol.allAtomCipDescriptors();

    for (AtomId id : mol.atomIds()) {
        if (atomToSgroup.contains(id)) continue;

        AtomPrim prim;
        prim.id = id;
        double x = 0, y = 0;
        mol.atomPos(id, x, y);
        prim.x = x; prim.y = y;

        QString symbol = mol.atomSymbol(id);
        prim.element = symbol;
        prim.isAtomList = mol.hasAtomQueryList(id);

        bool rgroupFlag = isRGroupLabel(symbol);
        prim.isRGroup = rgroupFlag;

        QList<AtomId> neighbors = mol.neighborAtomIds(id);
        bool isHetero = (symbol != QStringLiteral("C") && symbol != QStringLiteral("H"));
        // A bonded chain-end atom (exactly 1 neighbor) is a bare skeletal vertex like any other
        // carbon -- its position is already marked by the bond line, so it does not need its
        // own implicit-H label. Only a truly isolated atom (0 neighbors, nothing else marks
        // where it is) must always show its label, or it would render completely invisible.
        bool isIsolated = neighbors.isEmpty();
        if (rgroupFlag || prim.isAtomList) { isHetero = false; isIsolated = false; }

        QString renderLabel = symbol;
        int implicitH = mol.implicitHydrogenCount(id);
        if (prim.isAtomList) {
            QList<int> nums = mol.atomQueryListNumbers(id);
            QStringList labels;
            for (int n : nums) {
                QString sym = ElementData::symbolForNumber(n);
                labels.append(sym.isEmpty() ? QStringLiteral("?") : sym);
            }
            prim.atomListElements = labels.join(QStringLiteral(","));
            prim.atomListNot = mol.atomQueryListIsNotList(id);
            renderLabel = (prim.atomListNot ? QStringLiteral("!") : QString())
                          + QStringLiteral("[") + prim.atomListElements + QStringLiteral("]");
        } else if ((showExplicitH || isHetero || isIsolated) && implicitH > 0) {
            renderLabel += (implicitH == 1) ? QStringLiteral("H")
                                             : (QStringLiteral("H") + QString::number(implicitH));
        }
        prim.label = renderLabel;
        prim.implicitHCount = implicitH > 0 ? implicitH : 0;

        const ElementData::ElementInfo* info = ElementData::infoFor(symbol);
        prim.atomicNum = info ? info->number : 0;
        prim.atomicTitle = info ? info->title : QString();
        prim.atomicMass = info ? (std::round(info->mass * 1000.0) / 1000.0) : 0.0;
        prim.color = rgroupFlag ? QStringLiteral("#7B68EE") : ElementData::colorFor(symbol);

        prim.charge = mol.atomCharge(id);
        prim.isotope = mol.atomIsotope(id);
        prim.radical = mol.atomRadical(id);
        prim.explicitValence = mol.atomExplicitValence(id);
        prim.aam = mol.atomAAM(id);
        prim.checkWarning = mol.atomCheckWarningText(id);
        int attachOrder = mol.atomAttachmentOrder(id);
        prim.attachmentPoints = attachOrder > 0 ? (1 << (attachOrder - 1)) : 0;

        int stType = mol.stereocenterType(id);
        int stGroup = mol.stereocenterGroup(id);
        prim.stereoType = stType;
        prim.stereoGroup = stGroup;
        prim.stereoLabel = stereoLabelFor(stType, stGroup);
        prim.cipLabel = cipDescriptorLabel(atomCip.value(id, 0));

        double nbXSum = 0; int nbCount = 0;
        for (AtomId nb : neighbors) {
            double nx = 0, ny = 0;
            if (mol.atomPos(nb, nx, ny)) { nbXSum += nx; ++nbCount; }
        }
        prim.hOnLeft = nbCount > 0 && (nbXSum / nbCount) > x;

        result.atoms.append(prim);

        if (!haveBBox) { minX = maxX = x; minY = maxY = y; haveBBox = true; }
        else {
            if (x < minX) minX = x; if (x > maxX) maxX = x;
            if (y < minY) minY = y; if (y > maxY) maxY = y;
        }
    }

    // ---- Ring primitives ----
    QList<EditableMolecule::RingMembership> rings = mol.ringMembership();
    QHash<BondId, int> bondToFirstRing;   // "first ring wins" tie-break, matching the real code
    QHash<int, QPointF> ringCenterById;
    int ringPrimId = 0;
    for (const EditableMolecule::RingMembership& rm : rings) {
        double cx = 0, cy = 0;
        for (AtomId aid : rm.atoms) {
            double ax = 0, ay = 0;
            if (mol.atomPos(aid, ax, ay)) { cx += ax; cy += ay; }
        }
        if (!rm.atoms.isEmpty()) { cx /= rm.atoms.size(); cy /= rm.atoms.size(); }

        bool hasBondType4 = false;
        int doubleCount = 0;
        for (BondId bid : rm.bonds) {
            int order = mol.bondOrder(bid);
            if (order == 4) hasBondType4 = true;
            if (order == 2) ++doubleCount;
        }
        bool isAromatic = hasBondType4;
        if (!isAromatic && rm.atoms.size() == 6 && doubleCount == 3) isAromatic = true;

        double radius = 0;
        if (!rm.atoms.isEmpty()) {
            double ax = 0, ay = 0;
            if (mol.atomPos(rm.atoms[0], ax, ay)) {
                double dx = ax - cx, dy = ay - cy;
                radius = std::sqrt(dx * dx + dy * dy);
            }
        }

        RingPrim ringPrim;
        ringPrim.id = ringPrimId;
        ringPrim.atoms = rm.atoms;
        ringPrim.x = cx; ringPrim.y = cy; ringPrim.radius = radius;
        ringPrim.isAromatic = isAromatic;
        ringPrim.hasBondType4 = hasBondType4;
        result.rings.append(ringPrim);

        if (isAromatic) {
            for (BondId bid : rm.bonds) {
                if (!bondToFirstRing.contains(bid)) bondToFirstRing.insert(bid, ringPrimId);
            }
        }
        // Ring-center offset uses the FIRST ring found containing the bond, regardless of
        // aromaticity (matches buildRenderPrimitives: bondRingCenter is populated for every ring,
        // not only aromatic ones), so track it separately from the aromatic-only bondToFirstRing.
        for (BondId bid : rm.bonds) {
            if (!ringCenterById.contains(bid)) ringCenterById.insert(bid, QPointF(cx, cy));
        }
        ++ringPrimId;
    }
    // aromaticBondIds (distinct from ringCenterById): only bonds in an AROMATIC ring get
    // inAromaticRing=true, matching the real code's separate aromaticBondIds map.
    QSet<BondId> aromaticBondIds;
    for (auto it = bondToFirstRing.constBegin(); it != bondToFirstRing.constEnd(); ++it) {
        aromaticBondIds.insert(it.key());
    }

    // ---- Bond primitives ----
    const QMap<BondId, QString> bondCipLabels = mol.allBondCipLabels();
    for (BondId bid : mol.bondIds()) {
        AtomId ea = -1, eb = -1;
        if (!mol.bondEndpoints(bid, ea, eb)) continue;

        bool beginSg = atomToSgroup.contains(ea);
        bool endSg = atomToSgroup.contains(eb);
        // Fully-internal bond: both endpoints hidden by the SAME sgroup -- skip entirely,
        // matching buildRenderPrimitives' "if (beginSgId !== undefined && beginSgId === endSgId)
        // return". A cross-bond (exactly one endpoint hidden, OR both hidden by DIFFERENT
        // sgroups) is NOT skipped -- its hidden endpoint(s) are substituted with the owning
        // sgroup's synthetic id below, exactly matching the real code's own substitution.
        if (beginSg && endSg && atomToSgroup.value(ea) == atomToSgroup.value(eb)) continue;

        AtomId begin = beginSg ? static_cast<AtomId>(atomToSgroup.value(ea)) : ea;
        AtomId end = endSg ? static_cast<AtomId>(atomToSgroup.value(eb)) : eb;

        BondPrim prim;
        prim.id = bid;
        prim.begin = begin;
        prim.end = end;
        prim.beginIsSgroup = beginSg;
        prim.endIsSgroup = endSg;
        prim.type = mol.bondOrder(bid);
        prim.stereo = (beginSg || endSg) ? 0 : mol.bondStereoDirectionV2000(bid);
        prim.checkWarning = mol.bondCheckWarningText(bid);
        prim.cipLabel = bondCipLabels.value(bid);
        prim.reactingCenterStatus = 0;       // always 0 -- no reaction objects in this port

        prim.inAromaticRing = aromaticBondIds.contains(bid);

        if (ringCenterById.contains(bid)) {
            QPointF c = ringCenterById.value(bid);
            prim.hasRingCenter = true;
            prim.ringCenterX = c.x();
            prim.ringCenterY = c.y();
        }

        // Matches buildRenderPrimitives exactly: invalidStereo is computed (and can only be
        // non-false) when the bond is a real stereobond (prim.stereo > 0) with neither endpoint
        // sgroup-contracted, and it is the NEGATION of isCorrectStereoCenter (the real code:
        // "invalidStereo = !StereoValidator.isCorrectStereoCenter(...)"). Missing this negation
        // would invert the flag's meaning entirely.
        prim.invalidStereo = false;
        if (!beginSg && !endSg && prim.stereo > 0) {
            QList<AtomId> beginNeighbors = mol.neighborAtomIds(ea);
            QList<AtomId> endNeighbors = mol.neighborAtomIds(eb);
            int endOtherNeighborCount = 0;
            if (endNeighbors.size() == 2) {
                AtomId endOther = (endNeighbors[0] == ea) ? endNeighbors[1] : endNeighbors[0];
                endOtherNeighborCount = mol.neighborAtomIds(endOther).size();
            }
            prim.invalidStereo = !isCorrectStereoCenter(
                prim.stereo, beginNeighbors.size(), endNeighbors.size(),
                mol.implicitHydrogenCount(ea), endOtherNeighborCount);
        }

        result.bonds.append(prim);
    }

    // ---- Auxiliary primitives: texts, images, rxn arrows/pluses, multitail arrows ----
    // These need id lists, which EditableMolecule already exposes; each field is a direct
    // accessor call, no derived logic.
    for (int id : mol.rxnPlusIds()) {
        double x = 0, y = 0;
        if (!mol.rxnPlusPos(id, x, y)) continue;
        result.rxnPluses.append(RxnPlusPrim{id, x, y});
    }
    for (int id : mol.rxnArrowIds()) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        if (!mol.rxnArrowEndpoints(id, x1, y1, x2, y2)) continue;
        RxnArrowPrim prim;
        prim.id = id; prim.p1x = x1; prim.p1y = y1; prim.p2x = x2; prim.p2y = y2;
        prim.mode = mol.rxnArrowMode(id);
        prim.conditionsAbove = mol.rxnArrowConditionsAbove(id);
        prim.conditionsBelow = mol.rxnArrowConditionsBelow(id);
        double cx = 0, cy = 0;
        prim.hasCurvature = mol.rxnArrowCurvature(id, cx, cy);
        prim.curvatureX = cx; prim.curvatureY = cy;
        result.rxnArrows.append(prim);
    }
    for (int id : mol.multitailArrowIds()) {
        QList<double> pts = mol.multitailArrowPoints(id);
        if (pts.size() < 2) continue;
        MultitailArrowPrim prim;
        prim.id = id;
        prim.headX = pts[0];
        prim.headY = pts[1];
        for (int i = 2; i + 1 < pts.size(); i += 2) {
            prim.tails.append(QPointF(pts[i], pts[i + 1]));
        }
        result.multitailArrows.append(prim);
    }
    for (int id : mol.textAnnotationIds()) {
        double x = 0, y = 0; QString content; bool bold = false, italic = false;
        if (!mol.textAnnotationContent(id, x, y, content, bold, italic)) continue;
        result.texts.append(TextPrim{id, x, y, content, bold, italic});
    }
    for (int id : mol.imageIds()) {
        double x = 0, y = 0, w = 0, h = 0; QByteArray png;
        if (!mol.imageData(id, x, y, w, h, png)) continue;
        result.images.append(ImagePrim{id, x, y, w, h, png});
    }

    QList<int> rgNumbers = mol.rgroupNumbers();
    std::sort(rgNumbers.begin(), rgNumbers.end());
    for (int number : rgNumbers) {
        QString range; bool resth = false; int ifthen = 0;
        if (!mol.rgroupLogic(number, range, resth, ifthen)) continue;
        RGroupPrim rg;
        rg.number = number; rg.range = range; rg.resth = resth; rg.ifthen = ifthen;
        for (int fragId : mol.rgroupFragmentIds(number)) {
            rg.members.append(RGroupMemberPrim{fragId, mol.atomIdsInFragment(fragId)});
        }
        result.rgroups.append(rg);
    }

    for (int i = 0; i < mol.bracketCount(); ++i) {
        double bMinX = 0, bMinY = 0, bMaxX = 0, bMaxY = 0;
        if (!mol.bracketAt(i, bMinX, bMinY, bMaxX, bMaxY)) continue;
        result.brackets.append(BracketPrim{bMinX, bMinY, bMaxX, bMaxY});
    }

    result.bbox.valid = haveBBox;
    if (haveBBox) { result.bbox.minX = minX; result.bbox.minY = minY; result.bbox.maxX = maxX; result.bbox.maxY = maxY; }

    result.stereoFlagsType = mol.stereoFlagsType();
    result.stereoFlagsGroupId = mol.stereoFlagsGroupId();

    return result;
}
