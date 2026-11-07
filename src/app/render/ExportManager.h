// src/app/render/ExportManager.h
#ifndef EXPORTMANAGER_H
#define EXPORTMANAGER_H

#include <QObject>
#include <QString>
#include <QUrl>
#include <QVariantMap>
#include <QtQml/qqml.h>
#include "render/Exporter.h"

// QML-accessible export manager singleton
class ExportManager : public QObject {
    Q_OBJECT
    QML_ELEMENT
    QML_SINGLETON

public:
    explicit ExportManager(QObject* parent = nullptr);
    ~ExportManager();

    // Export a molecule (via molfile) to a file. Chemistry only -- any text
    // annotations, images, reaction arrows/pluses, multitail arrows, or
    // brackets the live document holds are NOT representable in a bare
    // molfile and so never reach the renderer through this entry point.
    Q_INVOKABLE void exportToFile(const QString& molfile, const QUrl& fileUrl,
                                   const QString& format, const QVariantMap& options = QVariantMap());

    // Export a full document snapshot (the same JSON ProjectFile::toJson/
    // V8Process::requestProjectSerialization produce) to a file. Unlike
    // exportToFile, this reconstructs the document's extension data (text,
    // images, reaction arrows/pluses, multitail arrows, brackets) as well as
    // its chemistry, drawing style and page settings, so the export matches
    // what the live canvas actually shows.
    Q_INVOKABLE void exportProjectToFile(const QString& projectJson, const QUrl& fileUrl,
                                          const QString& format, const QVariantMap& options = QVariantMap());

    // Export reaction grid to file
    Q_INVOKABLE void exportReactionGridToFile(const QString& molfile, const QUrl& fileUrl,
                                              const QString& format, const QVariantMap& options = QVariantMap());

    // Export batch of molfiles to file (as grid)
    Q_INVOKABLE void exportBatchGridToFile(const QStringList& molfiles, const QUrl& fileUrl,
                                           const QString& format, const QVariantMap& options = QVariantMap());

    // Export batch of molfiles to file (individual files or combined)
    Q_INVOKABLE void exportBatchToFile(const QStringList& molfiles, const QUrl& fileUrl,
                                       const QString& format, const QVariantMap& options = QVariantMap());

    // Get available export formats
    Q_INVOKABLE QStringList availableFormats() const;

    // Get default options for a format
    Q_INVOKABLE QVariantMap defaultOptions(const QString& format) const;

Q_SIGNALS:
    void exportFinished(bool success, const QString& error);
    void exportProgress(int percent);

private:
    // Create the appropriate exporter based on format
    Exporter* createExporter(const QString& format);

    // Helper to get format from file extension
    QString detectFormatFromPath(const QString& filePath) const;

    // Build RenderScene from molfile. options may contain "drawingStyle"/
    // "pageSettings" keys matching V8Process::getDrawingStyle()/
    // getPageSettings()'s JSON shape -- falls back to DrawingStyle::
    // defaults()/PageSettings::defaults() when absent.
    bool buildSceneFromMolfile(const QString& molfile, RenderScene& scene,
                               const QVariantMap& options = QVariantMap());

    // Build RenderScene from a full project JSON snapshot (ProjectFile schema),
    // reconstructing the document's extension data (text/images/arrows/etc.)
    // along with its chemistry, drawing style and page settings -- all taken
    // from the snapshot itself, not from `options`.
    bool buildSceneFromProjectJson(const QString& projectJson, RenderScene& scene);
};

#endif // EXPORTMANAGER_H
