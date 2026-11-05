// src/app/render/ExportManager.cpp
#include "ExportManager.h"
#include "Exporter.h"
#include "RenderScene.h"
#include "molecule/RenderPrimitives.h"
#include "molecule/DrawingStyle.h"
#include "molecule/PageSettings.h"
#include "indigo.h"
#include "indigo-renderer.h"
#include <QFile>
#include <QDebug>
#include <QtConcurrent>
#include <QPointer>
#include <QCoreApplication>

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

// Build RenderScene from molfile using Indigo to parse and layout
bool ExportManager::buildSceneFromMolfile(const QString& molfile, RenderScene& scene) {
    if (molfile.isEmpty()) return false;

    // Check if it's a reaction
    bool isReaction = molfile.trimmed().startsWith("$RXN") || molfile.contains(">>");

    // Use Indigo to parse the molfile
    unsigned long long sid = indigoAllocSessionId();
    indigoSetSessionId(sid);
    indigoRendererInit(sid);

    int mol = -1;
    if (isReaction) {
        mol = indigoLoadReactionFromString(molfile.toUtf8().constData());
    } else {
        mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
    }

    if (mol < 0) {
        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);
        return false;
    }

    // Layout if needed
    if (isReaction) {
        indigoLayout(mol);
    } else {
        indigoLayout(mol);
    }

    // Get molfile back (now with coordinates)
    const char* resultMol = isReaction ? indigoRxnfile(mol) : indigoMolfile(mol);
    QString layoutedMol = resultMol ? QString::fromUtf8(resultMol) : molfile;

    // TODO: Convert molfile to RenderPrimitives
    // For now, we'll use a simplified approach
    // In a full implementation, we would parse the molfile and build RenderPrimitives

    // Set up drawing style
    DrawingStyle drawingStyle;
    drawingStyle.modelToMmScale = 5.08 / 1.5;  // 1.5 model units = 5.08 mm
    drawingStyle.fontFamily = "Segoe UI Variable Text";
    drawingStyle.fontSize = 10;
    drawingStyle.lineWidth = 0.6;
    drawingStyle.monochrome = false;
    drawingStyle.bondLength = 1.5;

    // Set up page settings
    PageSettings pageSettings;
    pageSettings.marginLeftMm = 10;
    pageSettings.marginRightMm = 10;
    pageSettings.marginTopMm = 10;
    pageSettings.marginBottomMm = 10;
    pageSettings.widthMm = 210;
    pageSettings.heightMm = 297;

    // Build RenderPrimitives from the molfile
    // For now, create empty primitives - in a full implementation,
    // this would parse the molfile and build all the primitives
    RenderPrimitives prims;

    // Try to extract atom and bond information from Indigo
    if (!isReaction) {
        int it = indigoIterateAtoms(mol);
        if (it >= 0) {
            int atom;
            while ((atom = indigoNext(it)) > 0) {
                const char* symbol = indigoSymbol(atom);
                float* xyz = indigoXYZ(atom);
                
                if (xyz) {
                    AtomPrim aprim;
                    aprim.id = atom;
                    aprim.x = xyz[0];
                    aprim.y = xyz[1];
                    aprim.label = symbol ? QString::fromUtf8(symbol) : "C";
                    aprim.element = symbol ? QString::fromUtf8(symbol) : "C";
                    prims.atoms.append(aprim);
                }
            }
        }

        it = indigoIterateBonds(mol);
        if (it >= 0) {
            int bond;
            while ((bond = indigoNext(it)) > 0) {
                int begin = indigoSource(bond);
                int end = indigoDestination(bond);

                BondPrim bprim;
                bprim.id = bond;
                bprim.begin = begin;
                bprim.end = end;
                bprim.type = 1; // default to single bond
                prims.bonds.append(bprim);
            }
        }
    }

    // Calculate bounds
    if (!prims.atoms.isEmpty()) {
        double minX = prims.atoms[0].x;
        double maxX = prims.atoms[0].x;
        double minY = prims.atoms[0].y;
        double maxY = prims.atoms[0].y;
        
        for (const AtomPrim& atom : prims.atoms) {
            if (atom.x < minX) minX = atom.x;
            if (atom.x > maxX) maxX = atom.x;
            if (atom.y < minY) minY = atom.y;
            if (atom.y > maxY) maxY = atom.y;
        }
        
        prims.bbox.valid = true;
        prims.bbox.minX = minX;
        prims.bbox.maxX = maxX;
        prims.bbox.minY = minY;
        prims.bbox.maxY = maxY;
    }

    indigoFree(mol);
    indigoRendererDispose(sid);
    indigoReleaseSessionId(sid);

    // Build the scene
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
        if (!self || !self->buildSceneFromMolfile(molfile, scene)) {
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
