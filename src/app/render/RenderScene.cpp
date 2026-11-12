// src/app/render/RenderScene.cpp
#include "RenderScene.h"
#include <QPainter>
#include <QPolygonF>
#include <QFontDatabase>
#include <QDebug>
#include <QFile>
#include <QImage>
#include <cmath>
#include "molecule/RenderPrimitives.h"

namespace {
    // Constants for rendering
    const double kMmPerInch = 25.4;
    const double kPtsPerInch = 72.0;
    const double kDefaultDPI = 72.0;
    
    // Bond length constants (from chem-core.js)
    const double kStandardBondLength = 1.5;  // model units
    const double kBondSpacing = 0.15;  // spacing for double bonds
    const double kBondLineWidth = 0.6;  // points
    
    // Atom rendering
    const double kSgroupPillRadius = 0.8;
    const double kRGroupCircleRadius = 0.7;
    
    // Arrow rendering
    const double kArrowheadSize = 0.5;
    const double kArrowLineWidth = 0.6;
    const double kPlusSignSize = 0.6;
    
    // Bracket rendering
    const double kBracketPadding = 0.3;
    const double kBracketLineWidth = 0.6;

    // Wedge/hash stereo bond rendering. wedgeTaperRatio mirrors the live
    // canvas's own constant (MoleculeLayer.qml / StyleSheets.qml -- 0.139 is
    // the value shared by every built-in style sheet except RSC's 0.128) and
    // is a ratio of the bond's own length, so it is scale-invariant and
    // carries over to model-unit space unchanged. The live canvas instead
    // picks its hash count from a style-sheet pixel-spacing constant
    // calibrated for screen space at a given zoom; that constant has no
    // well-defined equivalent in model-unit export space, so a fixed count
    // is used here instead -- still a genuine tapered fan, just not
    // length-adaptive the way the live canvas's is.
    const double kWedgeTaperRatio = 0.139;
    const int kHashCount = 6;
}

RenderScene::RenderScene(QObject* parent) 
    : QObject(parent)
    , m_modelToMmScale(5.08 / 1.5)  // Default: 1.5 model units = 5.08 mm
    , m_regularFontMetrics(QFont())
    , m_boldFontMetrics(QFont())
{
}

RenderScene::~RenderScene() {
    clear();
}

bool RenderScene::build(const RenderPrimitives& prims,
                        const DrawingStyle& drawingStyle,
                        const PageSettings& pageSettings) {
    clear();
    
    m_prims = prims;
    m_drawingStyle = drawingStyle;
    m_pageSettings = pageSettings;
    m_modelToMmScale = drawingStyle.modelToMmScale;
    
    // Build model bounds from primitives
    if (prims.bbox.valid) {
        m_modelBounds = QRectF(
            prims.bbox.minX, prims.bbox.minY,
            prims.bbox.maxX - prims.bbox.minX,
            prims.bbox.maxY - prims.bbox.minY
        );
    } else {
        // Compute from atoms
        double minX = 0, minY = 0, maxX = 0, maxY = 0;
        bool first = true;
        for (const AtomPrim& atom : prims.atoms) {
            if (first) {
                minX = maxX = atom.x;
                minY = maxY = atom.y;
                first = false;
            } else {
                if (atom.x < minX) minX = atom.x;
                if (atom.x > maxX) maxX = atom.x;
                if (atom.y < minY) minY = atom.y;
                if (atom.y > maxY) maxY = atom.y;
            }
        }
        if (!first) {
            m_modelBounds = QRectF(minX, minY, maxX - minX, maxY - minY);
        }
    }
    
    // Set up fonts
    m_requestedFontFamily = drawingStyle.fontFamily;
    
    // Try to create the font
    m_regularFont = QFont(m_requestedFontFamily, drawingStyle.fontSize);
    m_boldFont = QFont(m_requestedFontFamily, drawingStyle.fontSize, QFont::Bold);
    
    // Check if font family was substituted
    QFontInfo regularInfo(m_regularFont);
    m_resolvedFontFamily = regularInfo.family();
    m_fontSubstituted = (m_requestedFontFamily != m_resolvedFontFamily);
    
    if (m_fontSubstituted) {
        m_fontWarnings.append(
            QStringLiteral("Font '%1' not available, using '%2'").arg(m_requestedFontFamily, m_resolvedFontFamily)
        );
    }
    
    // Measure on a 72 dpi device so metrics are in points.
    QImage pointDevice(1, 1, QImage::Format_ARGB32);
    pointDevice.setDotsPerMeterX(qRound(kPtsPerInch / 0.0254));
    pointDevice.setDotsPerMeterY(qRound(kPtsPerInch / 0.0254));
    m_regularFontMetrics = QFontMetricsF(m_regularFont, &pointDevice);
    m_boldFontMetrics = QFontMetricsF(m_boldFont, &pointDevice);

    // Expand bounds to cover every painted primitive, not just atom
    // positions -- an annotation, image, or bracket placed away from the
    // molecule (or an atom's own label, wider than a zero-size point) used
    // to be silently excluded and clipped in export. Requires fonts to
    // already be set up (measureText), so this runs after the block above.
    bool haveBounds = !m_prims.atoms.isEmpty() || prims.bbox.valid;
    auto includePoint = [&](double x, double y) {
        if (!haveBounds) {
            m_modelBounds = QRectF(x, y, 0, 0);
            haveBounds = true;
        } else {
            if (x < m_modelBounds.left()) m_modelBounds.setLeft(x);
            if (x > m_modelBounds.right()) m_modelBounds.setRight(x);
            if (y < m_modelBounds.top()) m_modelBounds.setTop(y);
            if (y > m_modelBounds.bottom()) m_modelBounds.setBottom(y);
        }
    };
    auto includeRect = [&](double x0, double y0, double x1, double y1) {
        includePoint(x0, y0);
        includePoint(x1, y1);
    };

    // Label boxes: points -> mm -> model units.
    const double ptToModelUnits = (m_modelToMmScale > 0.0)
        ? (kMmPerInch / kPtsPerInch) / m_modelToMmScale : 0.0;
    auto includeBox = [&](double cx, double cy, const QSizeF& sizePt) {
        const double hw = sizePt.width() * ptToModelUnits / 2.0;
        const double hh = sizePt.height() * ptToModelUnits / 2.0;
        includeRect(cx - hw, cy - hh, cx + hw, cy + hh);
    };
    for (const AtomPrim& atom : m_prims.atoms) {
        if (atom.isSgroup) continue;
        const double halfLabelPt = atom.label.isEmpty() ? 0.0 : measureText(atom.label).width() / 2;
        if (!atom.label.isEmpty())
            includeBox(atom.x, atom.y, measureText(atom.label));
        // Charge and isotope, placed as in paintAtoms.
        const double superY = atom.y + m_regularFontMetrics.ascent() / 2 * ptToModelUnits;
        if (atom.charge != 0) {
            const QSizeF sz = measureText(QString(std::abs(atom.charge), atom.charge > 0 ? '+' : '-'));
            includeBox(atom.x + (halfLabelPt + sz.width() / 2) * ptToModelUnits, superY, sz);
        }
        if (atom.isotope > 0) {
            const QSizeF sz = measureText(QString::number(atom.isotope));
            includeBox(atom.x - (halfLabelPt + sz.width() / 2) * ptToModelUnits, superY, sz);
        }
    }
    for (const SgroupPrim& sg : m_prims.sgroups) {
        includeBox(sg.x, sg.y, measureText(sg.label));
    }
    for (const TextPrim& text : m_prims.texts) {
        includeBox(text.x, text.y, measureText(text.content, text.bold) * 2.0);
    }
    for (const ImagePrim& img : m_prims.images) {
        includeRect(img.x, img.y, img.x + img.w, img.y + img.h);
    }
    for (const BracketPrim& br : m_prims.brackets) {
        includeRect(br.minX, br.minY, br.maxX, br.maxY);
    }
    for (const RxnArrowPrim& arrow : m_prims.rxnArrows) {
        includePoint(arrow.p1x, arrow.p1y);
        includePoint(arrow.p2x, arrow.p2y);
        // Arrowhead tip, drawn past p2, plus its outline.
        const QPointF d(arrow.p2x - arrow.p1x, arrow.p2y - arrow.p1y);
        const double len = std::hypot(d.x(), d.y());
        const double tip = kArrowheadSize + kArrowLineWidth * ptToModelUnits;
        if (len > 1e-9)
            includePoint(arrow.p2x + d.x() / len * tip, arrow.p2y + d.y() / len * tip);
        // Conditions text above and below the midpoint (2 pt gap each side).
        if (arrow.conditionsAbove.isEmpty() && arrow.conditionsBelow.isEmpty()) continue;
        const QSizeF above = measureText(arrow.conditionsAbove);
        const QSizeF below = measureText(arrow.conditionsBelow);
        const double textH = std::max(above.height(), below.height()) + 2.0;
        includeBox((arrow.p1x + arrow.p2x) / 2, (arrow.p1y + arrow.p2y) / 2,
                   QSizeF(std::max(above.width(), below.width()), textH * 2));
    }
    for (const RxnPlusPrim& plus : m_prims.rxnPluses) {
        includePoint(plus.x, plus.y);
    }
    for (const MultitailArrowPrim& ma : m_prims.multitailArrows) {
        includePoint(ma.headX, ma.headY);
        for (const QPointF& tail : ma.tails) {
            includePoint(tail.x(), tail.y());
        }
    }

    // A single atom, or atoms that are exactly collinear (e.g. a flat
    // horizontal 2-atom molecule), legitimately produce a zero-width and/or
    // zero-height bounds -- real geometry, not an error. But a degenerate
    // (zero-area) rect still allocates a degenerate (often zero-pixel)
    // export canvas despite there being real, visible content (the bond
    // line itself, line width, atom labels already accounted for above).
    // Pad by half the document's own bond length so every export has a
    // sane minimum extent.
    if (haveBounds) {
        const double minPad = m_drawingStyle.bondLength / 2.0;
        if (m_modelBounds.width() < minPad) {
            double cx = m_modelBounds.center().x();
            m_modelBounds.setLeft(cx - minPad / 2.0);
            m_modelBounds.setRight(cx + minPad / 2.0);
        }
        if (m_modelBounds.height() < minPad) {
            double cy = m_modelBounds.center().y();
            m_modelBounds.setTop(cy - minPad / 2.0);
            m_modelBounds.setBottom(cy + minPad / 2.0);
        }
    }

    return true;
}

QRectF RenderScene::mmBounds() const {
    return QRectF(
        m_modelBounds.left() * m_modelToMmScale,
        m_modelBounds.top() * m_modelToMmScale,
        m_modelBounds.width() * m_modelToMmScale,
        m_modelBounds.height() * m_modelToMmScale
    );
}

QRectF RenderScene::pxBounds(double dpi) const {
    const double mmToPx = dpi / kMmPerInch;
    QRectF mm = mmBounds();
    return QRectF(
        mm.left() * mmToPx,
        mm.top() * mmToPx,
        mm.width() * mmToPx,
        mm.height() * mmToPx
    );
}

void RenderScene::clear() {
    m_prims = RenderPrimitives();
    m_modelBounds = QRectF();
    m_requestedFontFamily.clear();
    m_resolvedFontFamily.clear();
    m_fontSubstituted = false;
    m_fontWarnings.clear();
}

void RenderScene::paint(QPainter* painter) const {
    if (!painter || !painter->isActive()) {
        return;
    }
    
    // Save painter state
    painter->save();
    
    // Set up coordinate system
    // Model coordinates are directly used, with scale applied
    // The caller is responsible for any additional transformations
    
    // Point-valued widths (DrawingStyle::lineWidth, kArrowLineWidth,
    // kBracketLineWidth, and the few hardcoded pen widths below) are stroke
    // widths in real points, NOT model units. The caller (Exporter) has
    // already scaled this painter's CTM by model-unit-to-device-unit before
    // calling paint(), so a raw point value handed straight to QPen would be
    // multiplied by that same scale again, inflating a 0.6pt stroke to
    // several points. ptToModel converts a point value into the painter's
    // current logical (model-unit) space so the CTM's scaling brings it back
    // to the correct absolute device width: points -> device units (via the
    // device's own DPI) -> model units (dividing by the CTM's own scale).
    const double ctmScale = qAbs(painter->worldTransform().m11());
    const double ptToModel = (ctmScale > 1e-9)
        ? (painter->device()->logicalDpiX() / kPtsPerInch) / ctmScale
        : 1.0;

    // Paint in the correct order (matching QML layering):
    // 1. Bonds (behind atoms)
    // 2. Rings
    // 3. Atoms
    // 4. Brackets
    // 5. R-group labels
    // 6. S-groups
    // 7. Text annotations
    // 8. Images
    // 9. Reaction arrows
    // 10. Reaction pluses
    // 11. Multitail arrows
    
    paintBonds(painter, ptToModel);
    paintRings(painter, ptToModel);
    paintAtoms(painter, ptToModel);
    paintBrackets(painter, ptToModel);
    paintRGroups(painter);
    paintSGroups(painter, ptToModel);
    paintTexts(painter);
    paintImages(painter, ptToModel);
    paintRxnArrows(painter, ptToModel);
    paintRxnPluses(painter, ptToModel);
    paintMultitailArrows(painter, ptToModel);
    
    // Restore painter state
    painter->restore();
}

void RenderScene::paintAtoms(QPainter* painter, double ptToModel) const {
    if (m_prims.atoms.isEmpty()) return;

    const double lineWidth = m_drawingStyle.lineWidth * ptToModel;
    const bool monochrome = m_drawingStyle.monochrome;
    
    for (const AtomPrim& atom : m_prims.atoms) {
        if (atom.isSgroup) {
            // S-group contracted atoms are handled in paintSGroups
            continue;
        }
        
        const QPointF pos(atom.x, atom.y);
        
        // Draw atom label
        if (!atom.label.isEmpty()) {
            painter->save();
            
            // Set color
            QColor color;
            if (monochrome) {
                color = Qt::black;
            } else if (!atom.color.isEmpty()) {
                color = QColor(atom.color);
            } else {
                color = Qt::black;  // default
            }
            painter->setPen(QPen(color, lineWidth));
            
            // Set font
            QFont font = m_regularFont;
            if (atom.element == QStringLiteral("C") && atom.label != QStringLiteral("C")) {
                // Non-carbon labels (H, O, N, etc.) or special labels
                // Keep regular font
            }
            // For chiral centers with implicit H, the label might be "CH" etc.
            drawLabel(painter, font, pos, atom.label, Qt::AlignHCenter | Qt::AlignVCenter);

            painter->restore();
        }

        // Charge and isotope sit at the label's upper corners.
        const double halfLabelPt = atom.label.isEmpty() ? 0.0 : measureText(atom.label).width() / 2;
        const double superPt = -m_regularFontMetrics.ascent() / 2;
        
        // Draw R-group circle for R-group atoms
        if (atom.isRGroup) {
            painter->save();
            QColor bgColor = QColor(QStringLiteral("#7B68EE"));  // Medium purple
            if (monochrome) bgColor = Qt::black;
            
            painter->setBrush(bgColor);
            painter->setPen(Qt::NoPen);
            
            const double radius = kRGroupCircleRadius * 0.5;
            painter->drawEllipse(pos, radius, radius);
            painter->restore();
        }
        
        // Draw charge, isotope, radical markers
        if (atom.charge != 0) {
            // Draw charge marker
            painter->save();
            QString chargeStr;
            if (atom.charge > 0) {
                chargeStr = QString(atom.charge, '+');
            } else {
                chargeStr = QString(-atom.charge, '-');
            }
            drawLabel(painter, m_regularFont, pos, chargeStr, Qt::AlignLeft | Qt::AlignVCenter,
                      QPointF(halfLabelPt, superPt));
            painter->restore();
        }
        
        if (atom.isotope > 0) {
            // Draw isotope number
            painter->save();
            QString isoStr = QString::number(atom.isotope);
            drawLabel(painter, m_regularFont, pos, isoStr, Qt::AlignRight | Qt::AlignVCenter,
                      QPointF(-halfLabelPt, superPt));
            painter->restore();
        }
    }
}

void RenderScene::paintBonds(QPainter* painter, double ptToModel) const {
    if (m_prims.bonds.isEmpty()) return;

    const double lineWidth = m_drawingStyle.lineWidth * ptToModel;
    const bool monochrome = m_drawingStyle.monochrome;
    const double bondLength = m_drawingStyle.bondLength;
    
    for (const BondPrim& bond : m_prims.bonds) {
        // Get atom positions
        QPointF startPos, endPos;
        QString startLabel, endLabel;
        bool startFound = false, endFound = false;

        for (const AtomPrim& atom : m_prims.atoms) {
            if (atom.id == bond.begin) {
                startPos = QPointF(atom.x, atom.y);
                startLabel = atom.isSgroup ? QString() : atom.label;
                startFound = true;
            }
            if (atom.id == bond.end) {
                endPos = QPointF(atom.x, atom.y);
                endLabel = atom.isSgroup ? QString() : atom.label;
                endFound = true;
            }
        }
        
        // Also check synthetic sgroup atoms
        if (!startFound || !endFound) {
            for (const SgroupPrim& sg : m_prims.sgroups) {
                if (sg.id == bond.begin) {
                    startPos = QPointF(sg.x, sg.y);
                    startFound = true;
                }
                if (sg.id == bond.end) {
                    endPos = QPointF(sg.x, sg.y);
                    endFound = true;
                }
            }
        }
        
        if (!startFound || !endFound) continue;

        // Stop the bond at the edge of a label box so it does not cross the text.
        // Returns the fraction of the bond hidden under the label.
        const QPointF d = endPos - startPos;
        auto hidden = [&](const QString& label) {
            if (label.isEmpty()) return 0.0;
            const double padPt = 1.0;
            const double hw = (measureText(label).width() / 2 + padPt) * ptToModel;
            const double hh = (m_regularFontMetrics.ascent() / 2 + padPt) * ptToModel;
            const double tx = std::abs(d.x()) > 1e-9 ? hw / std::abs(d.x()) : 1e9;
            const double ty = std::abs(d.y()) > 1e-9 ? hh / std::abs(d.y()) : 1e9;
            return std::min(tx, ty);
        };
        const double t0 = hidden(startLabel);
        const double t1 = hidden(endLabel);
        if (t0 + t1 >= 1.0) continue;
        endPos = startPos + d * (1.0 - t1);
        startPos = startPos + d * t0;

        painter->save();

        // Determine bond color
        QColor color = monochrome ? Qt::black : Qt::black;

        if (bond.stereo == 1) {
            // Wedge bond (up): a filled triangle, narrow at the begin atom and
            // widening toward the end atom -- matches the live canvas's own
            // wedge geometry (MoleculeLayer.qml), not a thick solid line.
            QPointF d = endPos - startPos;
            double len = std::sqrt(d.x() * d.x() + d.y() * d.y());
            double w = kWedgeTaperRatio * len;
            QPointF n = (len > 1e-9) ? QPointF(d.y() / len * w, -d.x() / len * w) : QPointF();
            QPolygonF wedge;
            wedge << startPos << (endPos + n) << (endPos - n);
            painter->setPen(Qt::NoPen);
            painter->setBrush(color);
            painter->drawPolygon(wedge);
        } else if (bond.stereo == 6) {
            // Hash bond (down): a fan of perpendicular marks widening toward
            // the end atom -- matches the live canvas's own hash geometry,
            // not a thick dashed line.
            QPointF d = endPos - startPos;
            double len = std::sqrt(d.x() * d.x() + d.y() * d.y());
            double maxW = kWedgeTaperRatio * len;
            painter->setPen(QPen(color, lineWidth));
            for (int j = 1; j <= kHashCount; ++j) {
                double t = double(j) / kHashCount;
                QPointF c = startPos + d * t;
                double hw = t * maxW;
                QPointF hn = (len > 1e-9) ? QPointF(d.y() / len * hw, -d.x() / len * hw) : QPointF();
                painter->drawLine(c + hn, c - hn);
            }
        } else {
            // Plain, double, or triple bond: unchanged QPen-based line(s).
            QPen pen(color, lineWidth);
            if (bond.type == 2) {
                pen.setWidthF(lineWidth * 1.5);
            } else if (bond.type == 3) {
                pen.setWidthF(lineWidth * 2);
            }
            painter->setPen(pen);
            painter->drawLine(startPos, endPos);

            if (bond.type == 2) {
                QPointF offset = computeDoubleBondOffset(startPos, endPos, bondLength);
                painter->drawLine(startPos + offset, endPos + offset);
            }
            if (bond.type == 3) {
                QPointF offset1 = computeDoubleBondOffset(startPos, endPos, bondLength);
                QPointF offset2 = -offset1;
                painter->drawLine(startPos + offset1, endPos + offset1);
                painter->drawLine(startPos + offset2, endPos + offset2);
            }
        }

        painter->restore();
    }
}

QPointF RenderScene::computeDoubleBondOffset(const QPointF& start, const QPointF& end, double bondLength) const {
    QPointF dir = end - start;
    double length = std::sqrt(dir.x() * dir.x() + dir.y() * dir.y());
    if (length < 0.001) return QPointF();
    
    // Normalize
    dir /= length;
    
    // Perpendicular direction (rotated 90 degrees)
    QPointF perp(-dir.y(), dir.x());
    
    // Scale by spacing
    return perp * kBondSpacing;
}

void RenderScene::paintRings(QPainter* painter, double ptToModel) const {
    if (m_prims.rings.isEmpty()) return;

    const double lineWidth = m_drawingStyle.lineWidth * ptToModel;
    const bool monochrome = m_drawingStyle.monochrome;
    
    for (const RingPrim& ring : m_prims.rings) {
        if (ring.isAromatic) {
            // Draw aromatic ring circle
            painter->save();
            QColor color = monochrome ? Qt::black : Qt::black;
            painter->setPen(QPen(color, lineWidth));
            painter->setBrush(Qt::NoBrush);
            painter->drawEllipse(QPointF(ring.x, ring.y), ring.radius, ring.radius);
            painter->restore();
        }
    }
}

void RenderScene::paintSGroups(QPainter* painter, double ptToModel) const {
    if (m_prims.sgroups.isEmpty()) return;
    
    const bool monochrome = m_drawingStyle.monochrome;
    
    for (const SgroupPrim& sg : m_prims.sgroups) {
        painter->save();
        
        // Draw sgroup pill
        QColor bgColor = QColor(200, 200, 200);  // Light gray
        QColor textColor = Qt::black;
        
        if (monochrome) {
            bgColor = Qt::black;
            textColor = Qt::white;
        }
        
        painter->setBrush(bgColor);
        painter->setPen(Qt::NoPen);
        
        // Calculate pill size based on label
        QSizeF labelSize = measureText(sg.label) * ptToModel;
        const double pillWidth = labelSize.width() + kSgroupPillRadius * 2;
        const double pillHeight = std::max(labelSize.height(), kSgroupPillRadius * 2);
        
        // Draw rounded rectangle (pill)
        QRectF pillRect(
            sg.x - pillWidth / 2,
            sg.y - pillHeight / 2,
            pillWidth,
            pillHeight
        );
        
        painter->drawRoundedRect(pillRect, kSgroupPillRadius, kSgroupPillRadius);
        
        // Draw label
        painter->setPen(QPen(textColor, 0.5 * ptToModel));
        drawLabel(painter, m_regularFont, QPointF(sg.x, sg.y), sg.label,
                  Qt::AlignHCenter | Qt::AlignVCenter);
        
        painter->restore();
    }
}

void RenderScene::paintTexts(QPainter* painter) const {
    if (m_prims.texts.isEmpty()) return;
    
    for (const TextPrim& text : m_prims.texts) {
        painter->save();
        
        QFont font = text.bold ? m_boldFont : m_regularFont;
        if (text.italic) {
            font.setItalic(true);
        }
        painter->setPen(Qt::black);
        drawLabel(painter, font, QPointF(text.x, text.y), text.content,
                  Qt::AlignHCenter | Qt::AlignBaseline);
        
        painter->restore();
    }
}

void RenderScene::paintImages(QPainter* painter, double ptToModel) const {
    if (m_prims.images.isEmpty()) return;

    for (const ImagePrim& image : m_prims.images) {
        QPointF pos(image.x, image.y);
        QSizeF size(image.w, image.h);

        // For now, just draw a placeholder rectangle
        // In a full implementation, we would decode the PNG data
        painter->save();
        painter->setBrush(QColor(220, 220, 220));
        painter->setPen(QPen(Qt::gray, 0.5 * ptToModel));
        painter->drawRect(QRectF(pos - QPointF(size.width() / 2, size.height() / 2), size));
        painter->restore();
    }
}

void RenderScene::paintRxnArrows(QPainter* painter, double ptToModel) const {
    if (m_prims.rxnArrows.isEmpty()) return;

    const double arrowheadSize = kArrowheadSize;
    const double lineWidth = kArrowLineWidth * ptToModel;
    
    for (const RxnArrowPrim& arrow : m_prims.rxnArrows) {
        QPointF p1(arrow.p1x, arrow.p1y);
        QPointF p2(arrow.p2x, arrow.p2y);
        
        painter->save();
        painter->setPen(QPen(Qt::black, lineWidth));
        painter->setBrush(Qt::black);
        
        // Draw line
        painter->drawLine(p1, p2);
        
        // Draw arrowhead
        if (arrow.mode == QStringLiteral("filled-triangle")) {
            drawArrowhead(painter, p2, p1, arrowheadSize);
        }
        
        // Conditions sit above and below the middle of the arrow.
        const QPointF mid = (p1 + p2) / 2;
        const double gapPt = 2.0;
        if (!arrow.conditionsAbove.isEmpty()) {
            drawLabel(painter, m_regularFont, mid, arrow.conditionsAbove,
                      Qt::AlignHCenter | Qt::AlignBottom, QPointF(0, -gapPt));
        }
        if (!arrow.conditionsBelow.isEmpty()) {
            drawLabel(painter, m_regularFont, mid, arrow.conditionsBelow,
                      Qt::AlignHCenter | Qt::AlignTop, QPointF(0, gapPt));
        }
        
        painter->restore();
    }
}

void RenderScene::paintRxnPluses(QPainter* painter, double ptToModel) const {
    if (m_prims.rxnPluses.isEmpty()) return;

    const double plusSize = kPlusSignSize;
    const double lineWidth = kArrowLineWidth * ptToModel;
    
    for (const RxnPlusPrim& plus : m_prims.rxnPluses) {
        QPointF center(plus.x, plus.y);
        double halfSize = plusSize / 2;
        
        painter->save();
        painter->setPen(QPen(Qt::black, lineWidth));
        
        // Draw plus sign
        painter->drawLine(
            center + QPointF(-halfSize, 0),
            center + QPointF(halfSize, 0)
        );
        painter->drawLine(
            center + QPointF(0, -halfSize),
            center + QPointF(0, halfSize)
        );
        
        painter->restore();
    }
}

void RenderScene::paintMultitailArrows(QPainter* painter, double ptToModel) const {
    if (m_prims.multitailArrows.isEmpty()) return;

    const double lineWidth = kArrowLineWidth * ptToModel;
    
    for (const MultitailArrowPrim& arrow : m_prims.multitailArrows) {
        QPointF head(arrow.headX, arrow.headY);
        
        painter->save();
        painter->setPen(QPen(Qt::black, lineWidth));
        
        // Draw head point
        painter->setBrush(Qt::black);
        painter->drawEllipse(head, 0.15, 0.15);
        
        // Draw tails
        for (const QPointF& tail : arrow.tails) {
            painter->drawLine(head, tail);
        }
        
        painter->restore();
    }
}

void RenderScene::paintBrackets(QPainter* painter, double ptToModel) const {
    if (m_prims.brackets.isEmpty()) return;

    const double lineWidth = kBracketLineWidth * ptToModel;
    const double padding = kBracketPadding;
    
    for (const BracketPrim& bracket : m_prims.brackets) {
        QRectF rect(bracket.minX, bracket.minY, bracket.maxX - bracket.minX, bracket.maxY - bracket.minY);
        rect = rect.marginsAdded(QMarginsF(padding, padding, padding, padding));
        
        painter->save();
        painter->setPen(QPen(Qt::black, lineWidth));
        
        // Draw rounded rectangle for bracket
        painter->setBrush(Qt::NoBrush);
        painter->drawRoundedRect(rect, 0.3, 0.3);
        
        painter->restore();
    }
}

void RenderScene::paintRGroups(QPainter* painter) const {
    if (m_prims.rgroups.isEmpty()) return;
    
    for (const RGroupPrim& rg : m_prims.rgroups) {
        // Brackets are drawn in paintBrackets; here only the label's anchor is needed.
        double minX = std::numeric_limits<double>::max();
        double maxY = -std::numeric_limits<double>::max();
        for (const RGroupMemberPrim& member : rg.members) {
            for (AtomId atomId : member.atomIds) {
                for (const AtomPrim& atom : m_prims.atoms) {
                    if (atom.id == atomId) {
                        minX = std::min(minX, atom.x);
                        maxY = std::max(maxY, atom.y);
                    }
                }
            }
        }

        // Label goes above the top-left of the member atoms.
        if (!rg.range.isEmpty() && minX != std::numeric_limits<double>::max()) {
            QString label = QStringLiteral("R") + QString::number(rg.number)
                          + QStringLiteral(" (") + rg.range + QStringLiteral(")");
            painter->save();
            painter->setPen(Qt::black);
            drawLabel(painter, m_regularFont, QPointF(minX, maxY), label,
                      Qt::AlignLeft | Qt::AlignBottom, QPointF(0, -m_regularFontMetrics.height()));
            painter->restore();
        }
    }
}

void RenderScene::drawArrowhead(QPainter* painter, const QPointF& tip, const QPointF& from, double size) const {
    // Calculate arrowhead direction
    QPointF dir = tip - from;
    double length = std::sqrt(dir.x() * dir.x() + dir.y() * dir.y());
    if (length < 0.001) return;
    
    dir /= length;
    
    // Perpendicular directions
    QPointF perp1(-dir.y(), dir.x());
    QPointF perp2(dir.y(), -dir.x());
    
    // Calculate arrowhead points
    QPointF p1 = tip + dir * size;
    QPointF p2 = tip + perp1 * size * 0.6;
    QPointF p3 = tip + perp2 * size * 0.6;
    
    // Draw filled triangle
    QPolygonF triangle;
    triangle << p1 << p2 << p3;
    painter->drawPolygon(triangle);
}

void RenderScene::drawLabel(QPainter* painter, const QFont& font, const QPointF& modelPos,
                            const QString& text, Qt::Alignment align,
                            const QPointF& offsetPt) const {
    // Fonts are in points, so draw in device space instead of the scaled, y-flipped model space.
    const double ptToDevice = painter->device()->logicalDpiY() / kPtsPerInch;
    const QPointF anchor = painter->worldTransform().map(modelPos) + offsetPt * ptToDevice;
    const QFontMetricsF fm(font, painter->device());

    double x = anchor.x();
    if (align & Qt::AlignHCenter) x -= fm.horizontalAdvance(text) / 2;
    else if (align & Qt::AlignRight) x -= fm.horizontalAdvance(text);

    double y = anchor.y();
    if (align & Qt::AlignVCenter) y += (fm.ascent() - fm.descent()) / 2;
    else if (align & Qt::AlignTop) y += fm.ascent();
    else if (align & Qt::AlignBottom) y -= fm.descent();

    painter->save();
    painter->resetTransform();
    painter->setFont(font);
    painter->drawText(QPointF(x, y), text);
    painter->restore();
}

QSizeF RenderScene::measureText(const QString& text, bool bold) const {
    const QFontMetricsF& metrics = bold ? m_boldFontMetrics : m_regularFontMetrics;
    return metrics.size(Qt::TextSingleLine, text);
}
