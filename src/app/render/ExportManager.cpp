// src/app/render/ExportManager.cpp
#include "ExportManager.h"
#include "Exporter.h"
#include "RenderScene.h"
#include "molecule/RenderPrimitives.h"
#include "molecule/EditableMolecule.h"
#include "molecule/DrawingStyle.h"
#include "molecule/PageSettings.h"
#include "indigo.h"
#include "indigo-renderer.h"
#include <QFile>
#include <QDebug>
#include <QtConcurrent>
#include <QPointer>
#include <QCoreApplication>
#include <QJsonObject>

// Forward declaration
static void emitOnGuiThread(QPointer<ExportManager> self, std::function<void(ExportManager*)> fn);

ExportManager::ExportManager(QObject* parent) : QObject(parent) {
}

ExportManager::~ExportManager() {
}

QStringList ExportManager::availableFormats() const {
    return {"pdf", "svg", "png"};
}

QVariantMap ExportManager::defaultOptions(const QString& format) const {
    QVariantMap options;
    if (format == "pdf") {
        options["pageWidthMm"] = 210;  // A4
        options["pageHeightMm"] = 297;
        options["landscape"] = false;
    } else if (format == "svg") {
        options["addMetadata"] = true;
        options["textSelectable"] = true;
    } else if (format == "png") {
        options["dpi"] = 600;
        options["backgroundColor"] = "white";
    }
    return options;
}

Exporter* ExportManager::createExporter(const QString& format) {
    if (format == "pdf") {
        return new PdfExporter(this);
    } else if (format == "svg") {
        return new SvgExporter(this);
    } else if (format == "png") {
        return new PngExporter(this);
    }
    return nullptr;
}

QString ExportManager::detectFormatFromPath(const QString& filePath) const {
    QString path = filePath.toLower();
    if (path.endsWith(".pdf")) return "pdf";
    if (path.endsWith(".svg")) return "svg";
    if (path.endsWith(".png")) return "png";
    return "png";  // default
}

namespace {
// Shared with buildSceneFromMolfile below: converts the document-style/page
// JSON the QML caller attaches to `options` (matching V8Process::
// getDrawingStyle()/getPageSettings()'s own JSON shape exactly) into real
// DrawingStyle/PageSettings values. Falls back to each field's own default
// when a key is absent, so a caller that doesn't supply style/page JSON at
// all (e.g. an existing test) still gets sane, documented defaults rather
// than failing.
DrawingStyle drawingStyleFromOptions(const QVariantMap& options) {
    DrawingStyle style = DrawingStyle::defaults();
    if (!options.contains("drawingStyle")) return style;
    // QML marshals a plain JS object into a QVariantMap (not QJsonObject)
    // when crossing into a QVariantMap-typed parameter -- .toMap() handles
    // that case, and also a QJsonObject passed explicitly from C++ test
    // code, unlike QVariant::toJsonObject() which only accepts the latter.
    QJsonObject obj = QJsonObject::fromVariantMap(options["drawingStyle"].toMap());
    if (obj.contains("bondLength")) style.bondLength = obj["bondLength"].toDouble();
    if (obj.contains("lineWidth")) style.lineWidth = obj["lineWidth"].toDouble();
    if (obj.contains("fontFamily")) style.fontFamily = obj["fontFamily"].toString();
    if (obj.contains("fontSize")) style.fontSize = obj["fontSize"].toInt();
    if (obj.contains("monochrome")) style.monochrome = obj["monochrome"].toBool();
    // modelToMmScale is fixed per-document (not derived from bondLength here)
    // -- same invariant DocumentState itself preserves.
    if (obj.contains("modelToMmScale")) style.modelToMmScale = obj["modelToMmScale"].toDouble();
    return style;
}

PageSettings pageSettingsFromOptions(const QVariantMap& options) {
    PageSettings page = PageSettings::defaults();
    if (!options.contains("pageSettings")) return page;
    QJsonObject obj = QJsonObject::fromVariantMap(options["pageSettings"].toMap());
    if (obj.contains("widthMm")) page.widthMm = obj["widthMm"].toDouble();
    if (obj.contains("heightMm")) page.heightMm = obj["heightMm"].toDouble();
    if (obj.contains("marginLeftMm")) page.marginLeftMm = obj["marginLeftMm"].toDouble();
    if (obj.contains("marginRightMm")) page.marginRightMm = obj["marginRightMm"].toDouble();
    if (obj.contains("marginTopMm")) page.marginTopMm = obj["marginTopMm"].toDouble();
    if (obj.contains("marginBottomMm")) page.marginBottomMm = obj["marginBottomMm"].toDouble();
    return page;
}
} // namespace

// Build RenderScene from a molfile, using the SAME EditableMolecule +
// RenderPrimitiveBuilder pipeline the live canvas already relies on --
// confirmed real, already-trusted conversion -- instead of a second,
// independent, lossy Indigo-to-primitives reparse (the old version here
// hardcoded every bond to single, dropped charge/isotope/stereo entirely,
// and used a hardcoded default style/page regardless of the real document).
// Reactions are explicitly out of scope for this single-molecule export path
// (exportReactionGridToFile/exportBatchGridToFile keep the existing
// Indigo-grid-based route for those, per this phase's own hybrid design).
bool ExportManager::buildSceneFromMolfile(const QString& molfile, RenderScene& scene,
                                           const QVariantMap& options) {
    if (molfile.isEmpty()) return false;

    EditableMolecule mol(molfile);
    if (!mol.isValid()) return false;

    RenderPrimitives prims = RenderPrimitiveBuilder::build(mol, /*showExplicitH=*/false);
    DrawingStyle drawingStyle = drawingStyleFromOptions(options);
    PageSettings pageSettings = pageSettingsFromOptions(options);

    return scene.build(prims, drawingStyle, pageSettings);
}

void ExportManager::exportToFile(const QString& molfile, const QUrl& fileUrl, 
                                  const QString& format, const QVariantMap& options) {
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) {
        emit exportFinished(false, "Invalid file path");
        return;
    }

    QPointer<ExportManager> self = this;
    QtConcurrent::run([self, molfile, path, format, options]() {
        bool success = false;
        QString error;

        // Detect format from path if not specified
        QString actualFormat = format;
        if (actualFormat.isEmpty()) {
            actualFormat = self ? self->detectFormatFromPath(path) : "png";
        }

        // Create exporter
        Exporter* exporter = self ? self->createExporter(actualFormat) : nullptr;
        if (!exporter) {
            error = "Unsupported format: " + actualFormat;
            emitOnGuiThread(self, [success, error](ExportManager* m) {
                emit m->exportFinished(success, error);
            });
            return;
        }

        // Set options
        exporter->setExportOptions(options);

        // Build scene from molfile
        RenderScene scene;
        if (!self || !self->buildSceneFromMolfile(molfile, scene, options)) {
            error = "Failed to build render scene from molfile";
            delete exporter;
            emitOnGuiThread(self, [success, error](ExportManager* m) {
                emit m->exportFinished(success, error);
            });
            return;
        }

        // Check for font substitution
        if (scene.fontSubstituted()) {
            qWarning() << "Font substitution detected:" 
                       << scene.requestedFontFamily() << "->" << scene.resolvedFontFamily();
        }

        // Export
        success = exporter->exportToFile(scene, path);
        if (!success) {
            error = exporter->lastError();
        }

        delete exporter;

        emitOnGuiThread(self, [success, error](ExportManager* m) {
            emit m->exportFinished(success, error);
        });
    });
}

void ExportManager::exportReactionGridToFile(const QString& molfile, const QUrl& fileUrl,
                                             const QString& format, const QVariantMap& options) {
    // For now, use the existing Indigo-based approach
    // In a full implementation, this would use the new renderer
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) {
        emit exportFinished(false, "Invalid file path");
        return;
    }

    QPointer<ExportManager> self = this;
    QtConcurrent::run([self, molfile, path, format, options]() {
        bool success = false;
        QString error;

        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoRendererInit(sid);

        try {
            int rx = indigoLoadReactionFromString(molfile.toUtf8().constData());
            if (rx >= 0) {
                int arr = indigoCreateArray();
                int count = 0;
                int comp;

                int reactants = indigoIterateReactants(rx);
                while ((comp = indigoNext(reactants)) > 0) { 
                    indigoArrayAdd(arr, comp); 
                    count++; 
                }
                int products = indigoIterateProducts(rx);
                while ((comp = indigoNext(products)) > 0) { 
                    indigoArrayAdd(arr, comp); 
                    count++; 
                }

                if (count > 0) {
                    indigoSetOption("render-output-format", format.toUtf8().constData());
                    int res = indigoRenderGridToFile(arr, nullptr, count, path.toUtf8().constData());
                    success = (res >= 0);
                    if (!success) {
                        error = QString::fromUtf8(indigoGetLastError());
                    }
                } else {
                    error = "Reaction has no reactants or products to render";
                }

                indigoFree(arr);
                indigoFree(rx);
            } else {
                error = QString("Failed to load reaction: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception& e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during grid render";
        }

        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);

        emitOnGuiThread(self, [success, error](ExportManager* m) {
            emit m->exportFinished(success, error);
        });
    });
}

void ExportManager::exportBatchGridToFile(const QStringList& molfiles, const QUrl& fileUrl,
                                          const QString& format, const QVariantMap& options) {
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty() || molfiles.isEmpty()) {
        emit exportFinished(false, molfiles.isEmpty() ? "No structures to export" : "Invalid file path");
        return;
    }

    QPointer<ExportManager> self = this;
    QtConcurrent::run([self, molfiles, path, format, options]() {
        bool success = false;
        QString error;

        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoRendererInit(sid);

        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;

            for (const QString& molfile : molfiles) {
                if (molfile.isEmpty()) continue;
                bool isRxn = molfile.trimmed().startsWith("$RXN") || molfile.contains(">>");
                int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                                : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    indigoArrayAdd(arr, mol);
                    loadedMols.append(mol);
                }
            }

            if (!loadedMols.isEmpty()) {
                indigoSetOption("render-output-format", format.toUtf8().constData());
                int res = indigoRenderGridToFile(arr, nullptr, loadedMols.size(), path.toUtf8().constData());
                success = (res >= 0);
                if (!success) {
                    error = QString::fromUtf8(indigoGetLastError());
                }
            } else {
                error = "Failed to load any molecules";
            }

            indigoFree(arr);
            for (int mol : loadedMols) {
                indigoFree(mol);
            }
        } catch (const std::exception& e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during batch grid render";
        }

        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);

        emitOnGuiThread(self, [success, error](ExportManager* m) {
            emit m->exportFinished(success, error);
        });
    });
}

void ExportManager::exportBatchToFile(const QStringList& molfiles, const QUrl& fileUrl,
                                     const QString& format, const QVariantMap& options) {
    // For batch export, use the new renderer for each molecule
    // This is different from grid export which puts everything on one page
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty() || molfiles.isEmpty()) {
        emit exportFinished(false, molfiles.isEmpty() ? "No structures to export" : "Invalid file path");
        return;
    }

    QPointer<ExportManager> self = this;
    
    // For batch export to a single file, we need to decide what to do
    // For now, just export the first molecule or use grid rendering
    // In a full implementation, this could create a multi-page PDF or multiple files
    
    if (molfiles.size() > 1) {
        // For multiple molecules, use Indigo's grid rendering
        exportBatchGridToFile(molfiles, fileUrl, format, options);
        return;
    }

    // Single molecule - use new renderer
    exportToFile(molfiles[0], fileUrl, format, options);
}

// Helper to emit signals on GUI thread
static void emitOnGuiThread(QPointer<ExportManager> self, std::function<void(ExportManager*)> fn) {
    QMetaObject::invokeMethod(qApp, [self, fn]() {
        if (self) fn(self.data());
    }, Qt::QueuedConnection);
}
