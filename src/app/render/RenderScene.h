// src/app/render/RenderScene.h
#ifndef RENDERSCENE_H
#define RENDERSCENE_H

#include <QObject>
#include <QRectF>
#include <QString>
#include <QList>
#include <QPointF>
#include <QFont>
#include <QFontMetricsF>
#include "molecule/RenderPrimitives.h"
#include "molecule/DrawingStyle.h"
#include "molecule/PageSettings.h"

// Shared rendering scene for canvas, preview, PDF, SVG, and PNG export.
// Built from RenderPrimitives and provides painting operations with a QPainter.
//
// Design:
// - Read-only: scene is built once from primitives and style, then painted
// - Physical dimensions: uses DrawingStyle's fixed coordinate scale
// - Font tracking: detects substitution between requested and resolved fonts

class QPainter;

class RenderScene : public QObject {
    Q_OBJECT

public:
    explicit RenderScene(QObject* parent = nullptr);
    ~RenderScene();

    // Build scene from primitives and style
    // Returns false if font cannot be loaded
    bool build(const RenderPrimitives& prims,
               const DrawingStyle& drawingStyle,
               const PageSettings& pageSettings);

    // Paint the scene with a QPainter
    // The painter should be active on a QPaintDevice (QImage, QPdfWriter, QSvgGenerator, etc.)
    void paint(QPainter* painter) const;

    // Scene bounds in model coordinates
    const QRectF& modelBounds() const { return m_modelBounds; }

    // Scene bounds in millimeters
    QRectF mmBounds() const;

    // Scene bounds in pixels at given DPI
    QRectF pxBounds(double dpi = 72.0) const;

    // Font information
    const QString& requestedFontFamily() const { return m_requestedFontFamily; }
    const QString& resolvedFontFamily() const { return m_resolvedFontFamily; }
    bool fontSubstituted() const { return m_fontSubstituted; }
    const QList<QString>& fontWarnings() const { return m_fontWarnings; }

    // Get the current font and metrics
    const QFont& regularFont() const { return m_regularFont; }
    const QFont& boldFont() const { return m_boldFont; }
    const QFontMetricsF& regularFontMetrics() const { return m_regularFontMetrics; }
    const QFontMetricsF& boldFontMetrics() const { return m_boldFontMetrics; }

    // Fixed scale from model to millimeters
    double modelToMmScale() const { return m_modelToMmScale; }

    // Clear the scene
    void clear();

private:
    // Helper to create a QFont with the requested family
    QFont createFont(const QString& family, int pointSize, bool bold = false) const;

    // Text size in points.
    QSizeF measureText(const QString& text, bool bold = false) const;

    // Draws upright text at a model point. offsetPt is in points, y down.
    void drawLabel(QPainter* painter, const QFont& font, const QPointF& modelPos,
                   const QString& text, Qt::Alignment align,
                   const QPointF& offsetPt = QPointF()) const;

    // Helper for bond rendering
    QPointF computeDoubleBondOffset(const QPointF& start, const QPointF& end, double bondLength) const;

    // Helper for arrowhead rendering
    void drawArrowhead(QPainter* painter, const QPointF& tip, const QPointF& from, double size) const;

    // Paint individual element types. `ptToModel` converts a point-valued stroke
    // width into the painter's current logical (model-unit) space, so that after
    // the caller's CTM scales back up to device units, the on-device width is the
    // intended point value -- see paint()'s own computation of it.
    void paintAtoms(QPainter* painter, double ptToModel) const;
    void paintBonds(QPainter* painter, double ptToModel) const;
    void paintRings(QPainter* painter, double ptToModel) const;
    void paintSGroups(QPainter* painter, double ptToModel) const;
    void paintTexts(QPainter* painter) const;
    void paintImages(QPainter* painter, double ptToModel) const;
    void paintRxnArrows(QPainter* painter, double ptToModel) const;
    void paintRxnPluses(QPainter* painter, double ptToModel) const;
    void paintMultitailArrows(QPainter* painter, double ptToModel) const;
    void paintBrackets(QPainter* painter, double ptToModel) const;
    void paintRGroups(QPainter* painter) const;

    // Member variables
    RenderPrimitives m_prims;
    DrawingStyle m_drawingStyle;
    PageSettings m_pageSettings;
    
    QRectF m_modelBounds;
    double m_modelToMmScale;
    
    QString m_requestedFontFamily;
    QString m_resolvedFontFamily;
    bool m_fontSubstituted;
    QList<QString> m_fontWarnings;
    
    QFont m_regularFont;
    QFont m_boldFont;
    QFontMetricsF m_regularFontMetrics;
    QFontMetricsF m_boldFontMetrics;
};

#endif // RENDERSCENE_H
