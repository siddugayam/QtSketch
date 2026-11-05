// src/app/render/Exporter.h
#ifndef EXPORTER_H
#define EXPORTER_H

#include <QObject>
#include <QString>
#include <QSizeF>
#include <QColor>
#include <QMarginsF>
#include "RenderScene.h"

// Base class for all exporters (PDF, SVG, PNG)
class Exporter : public QObject {
    Q_OBJECT

public:
    explicit Exporter(QObject* parent = nullptr) : QObject(parent) {}
    virtual ~Exporter() {}

    // Export the scene to a file
    // Returns true on success, false on failure
    virtual bool exportToFile(const RenderScene& scene, const QString& filePath) = 0;

    // Get export options as a JSON-serializable structure
    // Subclasses can override to provide custom options
    virtual QVariantMap exportOptions() const {
        return QVariantMap();
    }

    // Set export options from a JSON-serializable structure
    virtual void setExportOptions(const QVariantMap& options) {
        Q_UNUSED(options);
    }

    // Get the last error message
    QString lastError() const { return m_lastError; }

protected:
    void setError(const QString& error) { m_lastError = error; }
    void clearError() { m_lastError.clear(); }

private:
    QString m_lastError;
};

// PDF Exporter
class PdfExporter : public Exporter {
    Q_OBJECT

public:
    explicit PdfExporter(QObject* parent = nullptr);
    ~PdfExporter();

    bool exportToFile(const RenderScene& scene, const QString& filePath) override;

    // Options
    QVariantMap exportOptions() const override;
    void setExportOptions(const QVariantMap& options) override;

    // Page settings
    void setPageSize(const QSizeF& sizeMm);
    void setMargins(const QMarginsF& marginsMm);
    void setOrientation(bool landscape);

private:
    QSizeF m_pageSizeMm;
    QMarginsF m_marginsMm;
    bool m_landscape;
};

// SVG Exporter
class SvgExporter : public Exporter {
    Q_OBJECT

public:
    explicit SvgExporter(QObject* parent = nullptr);
    ~SvgExporter();

    bool exportToFile(const RenderScene& scene, const QString& filePath) override;

    // Options
    QVariantMap exportOptions() const override;
    void setExportOptions(const QVariantMap& options) override;

    // SVG-specific settings
    void setAddMetadata(bool add);
    void setTextSelectable(bool selectable);

private:
    bool m_addMetadata;
    bool m_textSelectable;
};

// PNG Exporter
class PngExporter : public Exporter {
    Q_OBJECT

public:
    explicit PngExporter(QObject* parent = nullptr);
    ~PngExporter();

    bool exportToFile(const RenderScene& scene, const QString& filePath) override;

    // Options
    QVariantMap exportOptions() const override;
    void setExportOptions(const QVariantMap& options) override;

    // PNG-specific settings
    void setResolution(int dpi);
    void setBackgroundColor(const QColor& color);

private:
    int m_dpi;
    QColor m_backgroundColor;
};

#endif // EXPORTER_H
