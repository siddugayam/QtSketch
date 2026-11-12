// src/app/render/Exporter.cpp
#include "Exporter.h"
#include <QPdfWriter>
#include <QSvgGenerator>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QDebug>

// ============================================================================
// PdfExporter
// ============================================================================

PdfExporter::PdfExporter(QObject* parent)
    : Exporter(parent)
    , m_pageSizeMm(210, 297)  // A4
    , m_landscape(false)
{
}

PdfExporter::~PdfExporter() {
}

bool PdfExporter::exportToFile(const RenderScene& scene, const QString& filePath) {
    clearError();
    
    // Get scene bounds in millimeters
    QRectF sceneBoundsMm = scene.mmBounds();
    
    // Create PDF writer
    QPdfWriter writer(filePath);
    
    // Set page size based on scene bounds or use configured page size
    // Add some margin
    const double marginMm = 10.0;
    
    QSizeF pageSizeMm = m_pageSizeMm;
    // m_pageSizeMm may arrive in either orientation; turn it only when it disagrees with
    // the requested one (the old `if (!m_landscape) transpose()` turned a portrait A4
    // request into a landscape page).
    if (m_landscape != (pageSizeMm.width() > pageSizeMm.height())) {
        pageSizeMm.transpose();
    }
    
    // If scene is larger than page, expand the page
    if (sceneBoundsMm.width() > pageSizeMm.width() - marginMm * 2 ||
        sceneBoundsMm.height() > pageSizeMm.height() - marginMm * 2) {
        // Use scene size + margins
        pageSizeMm = QSizeF(
            sceneBoundsMm.width() + marginMm * 2,
            sceneBoundsMm.height() + marginMm * 2
        );
    }
    
    // Set page size directly in mm
    writer.setPageSize(QPageSize(QSizeF(pageSizeMm.width(), pageSizeMm.height()), QPageSize::Millimeter));

    // Set resolution (DPI)
    writer.setResolution(300);

    // Create painter
    QPainter painter(&writer);
    if (!painter.isActive()) {
        setError(QStringLiteral("Failed to create PDF painter"));
        return false;
    }

    // Query the painter's REAL device DPI rather than assuming a fixed 72 --
    // QPdfWriter's painter coordinate space tracks setResolution() above, so
    // hardcoding 72/25.4 here would silently scale everything wrong whenever
    // the resolution is anything other than 72 (same pattern RenderScene::
    // paint() itself already uses for its own line-width calculation).
    const double mmToPoints = painter.device()->logicalDpiX() / 25.4;
    
    // Set up coordinate system
    // PDF coordinates: origin at top-left, y increases downward
    // Model coordinates: origin at center, y increases upward
    // We need to transform
    
    // Calculate scale to fit scene in page
    QRectF pageRectMm(0, 0, pageSizeMm.width(), pageSizeMm.height());
    QRectF centeredScene = sceneBoundsMm;
    
    // Center the scene
    double offsetX = (pageSizeMm.width() - sceneBoundsMm.width()) / 2 - sceneBoundsMm.left();
    double offsetY = (pageSizeMm.height() - sceneBoundsMm.height()) / 2 - sceneBoundsMm.top();
    
    // Convert mm to device coordinates
    // In PDF, the coordinate system is already in points
    // But we need to scale from mm to the PDF page
    // For now, use a simple approach: scale to fit
    
    painter.save();
    
    // Translate to center. The scale below flips Y, so the scene's largest model y (the
    // rect's bottom()) is the edge that lands at the top of the centred box: device y =
    // ty - y, hence ty = topMargin + bottom(). Anchoring on top() instead is only right
    // for a scene whose y range is symmetric about 0.
    painter.translate(
        (writer.width() - sceneBoundsMm.width() * mmToPoints) / 2 - sceneBoundsMm.left() * mmToPoints,
        (writer.height() - sceneBoundsMm.height() * mmToPoints) / 2 + sceneBoundsMm.bottom() * mmToPoints
    );
    
    // scene.paint() works in raw model units; the full chain to device
    // points is modelUnit -> mm (modelToMmScale) -> points (mmToPoints).
    const double modelToPoints = scene.modelToMmScale() * mmToPoints;
    painter.scale(modelToPoints, -modelToPoints);  // Negative Y to flip coordinate system

    // Paint the scene
    scene.paint(&painter);

    painter.restore();
    painter.end();

    return true;
}

QVariantMap PdfExporter::exportOptions() const {
    QVariantMap options;
    options["pageWidthMm"] = m_pageSizeMm.width();
    options["pageHeightMm"] = m_pageSizeMm.height();
    options["landscape"] = m_landscape;
    return options;
}

void PdfExporter::setExportOptions(const QVariantMap& options) {
    if (options.contains("pageWidthMm")) {
        m_pageSizeMm.setWidth(options["pageWidthMm"].toDouble());
    }
    if (options.contains("pageHeightMm")) {
        m_pageSizeMm.setHeight(options["pageHeightMm"].toDouble());
    }
    if (options.contains("landscape")) {
        m_landscape = options["landscape"].toBool();
    }
}

void PdfExporter::setPageSize(const QSizeF& sizeMm) {
    m_pageSizeMm = sizeMm;
}

void PdfExporter::setMargins(const QMarginsF& marginsMm) {
    m_marginsMm = marginsMm;
}

void PdfExporter::setOrientation(bool landscape) {
    m_landscape = landscape;
}

// ============================================================================
// SvgExporter
// ============================================================================

SvgExporter::SvgExporter(QObject* parent)
    : Exporter(parent)
    , m_addMetadata(true)
    , m_textSelectable(true)
{
}

SvgExporter::~SvgExporter() {
}

bool SvgExporter::exportToFile(const RenderScene& scene, const QString& filePath) {
    clearError();
    
    // Get scene bounds in millimeters
    QRectF sceneBoundsMm = scene.mmBounds();
    
    // Create SVG generator
    QSvgGenerator generator;
    generator.setFileName(filePath);

    // Query the generator's own logical DPI rather than assuming a fixed 72
    // -- QSvgGenerator does not necessarily default to 72 (same reasoning as
    // PdfExporter's identical fix: never hardcode a DPI a device didn't
    // confirm it's actually using).
    const double mmToPoints = generator.logicalDpiX() / 25.4;
    QSize sizePt(
        std::ceil(sceneBoundsMm.width() * mmToPoints),
        std::ceil(sceneBoundsMm.height() * mmToPoints)
    );
    generator.setSize(sizePt);
    
    // The painter transform below already maps the scene's bounds onto
    // (0,0)-(width,height) in device units, so the viewBox must start at the origin.
    // Offsetting it by the scene's own position shifted the visible window away from
    // the drawing for any scene not touching the model origin.
    generator.setViewBox(QRectF(
        0, 0,
        sceneBoundsMm.width() * mmToPoints,
        sceneBoundsMm.height() * mmToPoints
    ));
    
    // Add metadata if requested
    if (m_addMetadata) {
        // Add description with font information
        QString description = QStringLiteral("QtSketch Chemical Structure");
        if (scene.fontSubstituted()) {
            description += QStringLiteral(" - Font substituted: ") + 
                          scene.requestedFontFamily() + QStringLiteral(" -> ") + 
                          scene.resolvedFontFamily();
        }
        
        // SVG metadata would be added here
        // For now, we'll just note the font substitution in the console
        if (scene.fontSubstituted()) {
            qWarning() << "Font substitution:" << scene.requestedFontFamily()
                       << "->" << scene.resolvedFontFamily();
        }
    }
    
    // Create painter
    QPainter painter;
    if (!painter.begin(&generator)) {
        setError(QStringLiteral("Failed to create SVG painter"));
        return false;
    }
    
    // Set up coordinate system
    // SVG: origin at top-left, y increases downward
    // Model: origin can be anywhere, y increases upward
    // We need to flip Y
    
    painter.save();

    // scene.paint() works in raw model units -- the full chain to device
    // points is modelUnit -> mm (modelToMmScale) -> points (mmToPoints).
    const double modelToPoints = scene.modelToMmScale() * mmToPoints;
    painter.scale(modelToPoints, -modelToPoints);

    // Translate to position the scene. sceneBoundsMm is in mm; dividing by
    // modelToMmScale converts it back to model units so this translate (now
    // running in the model-to-points-scaled coordinate system set up above)
    // shifts by the correct device amount, not modelToMmScale times too much.
    painter.translate(-sceneBoundsMm.left() / scene.modelToMmScale(),
                       -sceneBoundsMm.bottom() / scene.modelToMmScale());

    // Paint the scene
    scene.paint(&painter);
    
    painter.restore();
    painter.end();
    
    return true;
}

QVariantMap SvgExporter::exportOptions() const {
    QVariantMap options;
    options["addMetadata"] = m_addMetadata;
    options["textSelectable"] = m_textSelectable;
    return options;
}

void SvgExporter::setExportOptions(const QVariantMap& options) {
    if (options.contains("addMetadata")) {
        m_addMetadata = options["addMetadata"].toBool();
    }
    if (options.contains("textSelectable")) {
        m_textSelectable = options["textSelectable"].toBool();
    }
}

void SvgExporter::setAddMetadata(bool add) {
    m_addMetadata = add;
}

void SvgExporter::setTextSelectable(bool selectable) {
    m_textSelectable = selectable;
}

// ============================================================================
// PngExporter
// ============================================================================

PngExporter::PngExporter(QObject* parent)
    : Exporter(parent)
    , m_dpi(600)
    , m_backgroundColor(Qt::white)
{
}

PngExporter::~PngExporter() {
}

bool PngExporter::exportToFile(const RenderScene& scene, const QString& filePath) {
    clearError();
    
    // Get scene bounds in pixels at the specified DPI
    QRectF sceneBoundsPx = scene.pxBounds(m_dpi);
    
    // Round to integer dimensions
    QSize imageSize(
        std::ceil(sceneBoundsPx.width()),
        std::ceil(sceneBoundsPx.height())
    );
    
    // Create image
    QImage image(imageSize, QImage::Format_ARGB32);
    image.fill(m_backgroundColor);
    
    // Create painter
    QPainter painter(&image);
    if (!painter.isActive()) {
        setError(QStringLiteral("Failed to create PNG painter"));
        return false;
    }
    
    // Set up coordinate system
    // PNG: origin at top-left, y increases downward
    // Model: y increases upward
    // We need to flip Y
    
    painter.save();
    
    // Set DPI for the painter
    painter.setRenderHint(QPainter::TextAntialiasing);
    painter.setRenderHint(QPainter::Antialiasing);
    
    // Scale (with the Y flip) first, then shift in model units so the scene's left edge
    // lands on pixel column 0 and its largest model y on pixel row 0 -- the same order
    // SvgExporter uses. The previous device-space translate (-left, -top + height) only
    // lined up for a scene whose smallest y is exactly 0; anything else was drawn partly
    // or wholly outside the image.
    // scene.paint() works in raw model units -- the full chain to pixels is
    // modelUnit -> mm (modelToMmScale) -> pixels (dpi / 25.4). A bare
    // scale(1, -1) (Y-flip only, no actual magnification) was the most
    // severe of the three exporters' unit bugs: every export came out at
    // 1 pixel per model unit regardless of the requested DPI.
    const double pxPerModelUnit = scene.modelToMmScale() * (m_dpi / 25.4);
    painter.scale(pxPerModelUnit, -pxPerModelUnit);
    painter.translate(-scene.modelBounds().left(), -scene.modelBounds().bottom());

    // Paint the scene
    scene.paint(&painter);
    
    painter.restore();
    painter.end();
    
    // Save image
    if (!image.save(filePath, "PNG")) {
        setError(QStringLiteral("Failed to save PNG file"));
        return false;
    }
    
    return true;
}

QVariantMap PngExporter::exportOptions() const {
    QVariantMap options;
    options["dpi"] = m_dpi;
    options["backgroundColor"] = m_backgroundColor.name();
    return options;
}

void PngExporter::setExportOptions(const QVariantMap& options) {
    if (options.contains("dpi")) {
        m_dpi = options["dpi"].toInt();
    }
    if (options.contains("backgroundColor")) {
        m_backgroundColor = QColor(options["backgroundColor"].toString());
    }
}

void PngExporter::setResolution(int dpi) {
    m_dpi = dpi;
}

void PngExporter::setBackgroundColor(const QColor& color) {
    m_backgroundColor = color;
}
