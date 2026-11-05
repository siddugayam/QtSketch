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
    if (!m_landscape) {
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
    
    // Convert mm to points (1 point = 1/72 inch, 25.4 mm = 1 inch)
    const double mmToPoints = 72.0 / 25.4;
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
    
    // Translate to center
    painter.translate(
        (writer.width() - sceneBoundsMm.width() * mmToPoints) / 2 - sceneBoundsMm.left() * mmToPoints,
        (writer.height() - sceneBoundsMm.height() * mmToPoints) / 2 - sceneBoundsMm.top() * mmToPoints
    );
    
    // Scale from mm to points
    painter.scale(mmToPoints, -mmToPoints);  // Negative Y to flip coordinate system
    
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
    
    // Set size (in points)
    const double mmToPoints = 72.0 / 25.4;
    QSize sizePt(
        std::ceil(sceneBoundsMm.width() * mmToPoints),
        std::ceil(sceneBoundsMm.height() * mmToPoints)
    );
    generator.setSize(sizePt);
    
    // Set viewbox to match the scene bounds
    generator.setViewBox(QRectF(
        sceneBoundsMm.left() * mmToPoints,
        -sceneBoundsMm.top() * mmToPoints,  // Flip Y for SVG
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
    
    // Scale from model mm to points
    painter.scale(mmToPoints, -mmToPoints);
    
    // Translate to position the scene
    painter.translate(-sceneBoundsMm.left(), -sceneBoundsMm.bottom());
    
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
    
    // Translate and scale
    // First, translate to position the scene
    painter.translate(-sceneBoundsPx.left(), -sceneBoundsPx.top());
    
    // Then flip Y to match the model coordinate system
    // The scene's paint method expects y-up, but the image has y-down
    // So we need to translate and scale
    painter.translate(0, sceneBoundsPx.height());
    painter.scale(1, -1);
    
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
