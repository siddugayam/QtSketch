#include "IndigoService.h"
#include "IupacNamer.h"
#include "BondStereoPerception.h"
#include <cmath>
#include <QDebug>
#include <QGuiApplication>
#include <QClipboard>
#include <QCoreApplication>
#include <QFile>
#include <QSaveFile>
#include <QFileInfo>
#include <QDir>
#include <mutex>
#include <functional>
#include "indigo.h"
#include "indigo-inchi.h"
#include "indigo-renderer.h"
#include <QtConcurrent>
#include <QPointer>
#include <QPointF>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QRegularExpression>
#include <QSet>
#include <QJsonArray>

// Read the monomer_library.ket once at first use; cache the bytes.
// Must be called with an active Indigo session — returns a library handle or -1.
static const QByteArray& monomerLibraryContent() {
    static QByteArray content;
    static std::once_flag flag;
    std::call_once(flag, []() {
        const QString rel = "/indigo/data/molecules/basic/monomer_library.ket";
        const QString appDir = QCoreApplication::applicationDirPath();
        const QString cwd    = QDir::currentPath();
        QStringList candidates = {
            QDir::cleanPath(appDir + "/../.." + rel),   // build/<config>/ → project root
            QDir::cleanPath(appDir + "/.."    + rel),   // build/ → project root
            QDir::cleanPath(cwd              + rel),    // CWD is project root
            QDir::cleanPath(cwd   + "/.."    + rel),    // CWD is build/
        };
        for (const QString& path : candidates) {
            QFile f(path);
            if (f.open(QIODevice::ReadOnly)) {
                content = f.readAll();
                return;
            }
        }
        qWarning() << "Indigo: monomer library not found. Tried:" << candidates;
    });
    return content;
}

// Strips cosmetic "M  SDS" (Sgroup display/expanded-state) lines before re-parsing a molfile
// for analysis. These lines are display-only (contracted vs. expanded Sgroup brackets) and
// never affect atoms/bonds/properties -- but a batched multi-index "M  SDS EXP" line (more than
// one Sgroup per line) breaks Indigo's strict fixed-column V2000 reader if the count field spills
// past its 3-character width (happens once 10+ Sgroups share a line), which the app's own
// re-serialization of an expanded biopolymer document can produce. Safe to drop outright for
// this read-only analysis pass since only the original document (unaffected) is ever displayed
// or saved -- this sanitized copy is used solely for indigoLoadMoleculeFromString here.
static QString sanitizeMolfileForAnalysis(const QString &molfile) {
    static const QRegularExpression sdsLine(
        QStringLiteral("^M  SDS.*$\\n?"),
        QRegularExpression::MultilineOption);
    QString out = molfile;
    out.remove(sdsLine);
    return out;
}

// The worker's MOL serializer switches to $RXN (reaction) format whenever the
// structure has a reaction arrow — indigoLoadMoleculeFromString() cannot parse
// that (it expects a single-molecule counts line, not "$RXN"'s multi-$MOL
// header) and fails with a scanner "end of stream" error. Molecule-level
// properties/stereo-descriptors/validation don't have a single well-defined
// value for a whole multi-molecule reaction scheme, so these are skipped
// gracefully here rather than attempted — same treatment as an empty molfile.
static bool isReactionFormat(const QString &data) {
    QString trimmed = data.trimmed();
    // ">>" is reaction-SMILES's reactant/product separator. A real molfile's first four
    // lines are a fixed structural header (never containing ">>"), and real molecule
    // SMILES/InChI text essentially never contains a literal ">>" substring, so this check
    // is safe as a pure widening -- every caller of isReactionFormat (checkStructure,
    // calcProperties, renderReactionGridToFile, layout, and this task's new importReaction
    // path) benefits identically.
    return trimmed.startsWith(QLatin1String("$RXN")) || trimmed.contains(QLatin1String(">>"));
}

// The vendored indigo.dll silently crashes the whole process (not a catchable C++ exception --
// confirmed empirically: an unrecognized character like a digit produces a clean "Invalid
// symbols" error, but a real IUPAC ambiguity code like 'X' or 'B' crashes with zero exception/
// stderr output) when a sequence/FASTA load contains an IUPAC ambiguity code. Root cause is
// inside the prebuilt binary (built from a fork that doesn't match the vendored v1.45.0
// reference source available for reading), so it can't be fixed there -- reject unsafe input
// before it ever reaches Indigo instead.
static bool validateBioSequenceAlphabet(const QString &text, const QString &seqType, bool isFasta, QString &badChars) {
    static const QSet<QChar> peptideSafe = {
        'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'
    };
    static const QSet<QChar> rnaSafe = { 'A','C','G','U' };
    static const QSet<QChar> dnaSafe = { 'A','C','G','T' };

    const QSet<QChar>* safe = &peptideSafe;
    if (seqType == QLatin1String("RNA")) safe = &rnaSafe;
    else if (seqType == QLatin1String("DNA")) safe = &dnaSafe;

    QString body = text;
    if (isFasta) {
        // Strip FASTA header lines (start with '>') before validating -- their free-text
        // content (e.g. ">seq1") is not sequence data and must not be checked against the
        // residue alphabet. Same convention as src/worker/70-biopolymer.js's own FASTA handling.
        QStringList kept;
        for (const QString &line : text.split(QLatin1Char('\n'))) {
            if (line.trimmed().startsWith(QLatin1Char('>'))) continue;
            kept << line;
        }
        body = kept.join(QLatin1Char('\n'));
    }

    QSet<QChar> bad;
    for (QChar rawCh : body) {
        QChar ch = rawCh.toUpper();
        if (!ch.isLetter()) continue; // non-letters (digits, '-', etc.) already produce a clean,
                                       // catchable Indigo error -- confirmed safe, don't touch them
        if (!safe->contains(ch)) bad.insert(ch);
    }
    if (bad.isEmpty()) return true;

    QStringList list;
    for (QChar c : bad) list << QString(c);
    badChars = list.join(QLatin1Char(','));
    return false;
}

// QPointer's guard is documented as safe only when checked from the pointed-to
// object's own thread (or with extra synchronization) -- checking it directly from
// a QtConcurrent background thread and then emitting on it is a real, if narrow,
// race against ~IndigoService() running concurrently on the GUI thread (mainly at
// app shutdown, if a slow Indigo call is still in flight). This posts the
// check-and-emit onto the GUI thread instead, where both the check and any
// concurrent destruction are serialized by the same single-threaded event loop.
static void emitOnGuiThread(QPointer<IndigoService> self, std::function<void(IndigoService *)> fn) {
    QMetaObject::invokeMethod(qApp, [self, fn]() {
        if (self) fn(self.data());
    }, Qt::QueuedConnection);
}

static int loadMonomerLibrary() {
    const QByteArray& bytes = monomerLibraryContent();
    if (bytes.isEmpty()) return -1;
    return indigoLoadMonomerLibraryFromString(bytes.constData());
}

// Shared helper: expand monomer handles to atoms, run layout, export molfile.
// Must be called with the correct Indigo session already set.
static QString indigoBioExpand(int mol) {
    QString result;
    if (mol < 0) return result;
    indigoExpandMonomers(mol);
    indigoExpandedMonomersToAtoms(mol);
    indigoLayout(mol);
    const char* mf = indigoMolfile(mol);
    if (mf) result = QString::fromUtf8(mf);
    indigoFree(mol);
    return result;
}

// indigoLayout/indigoClean2d relayout atoms around whatever origin the algorithm
// picks internally (typically near 0,0), which silently relocates the whole
// structure on the page even though only its geometry, not its position, was
// asked to change. Averages atom xyz (molecule) or every reaction molecule's
// atoms (reaction) into one centroid so callers can shift the result back.
static QPointF indigoCentroid(int mol, bool isRxn) {
    double sumX = 0, sumY = 0;
    int count = 0;
    auto accumulate = [&](int m) {
        int it = indigoIterateAtoms(m);
        if (it < 0) return;
        int a;
        while ((a = indigoNext(it)) > 0) {
            float* xyz = indigoXYZ(a);
            if (xyz) {
                sumX += xyz[0];
                sumY += xyz[1];
                count++;
            }
            indigoFree(a);
        }
        indigoFree(it);
    };
    if (isRxn) {
        int mit = indigoIterateMolecules(mol);
        if (mit >= 0) {
            int m;
            while ((m = indigoNext(mit)) > 0) { accumulate(m); indigoFree(m); }
            indigoFree(mit);
        }
    } else {
        accumulate(mol);
    }
    if (count == 0) return QPointF(0, 0);
    return QPointF(sumX / count, sumY / count);
}

// Shifts every atom (and, for reactions, every molecule's atoms) by (dx, dy).
static void indigoShiftAll(int mol, bool isRxn, double dx, double dy) {
    if (dx == 0.0 && dy == 0.0) return;
    auto shift = [&](int m) {
        int it = indigoIterateAtoms(m);
        if (it < 0) return;
        int a;
        while ((a = indigoNext(it)) > 0) {
            float* xyz = indigoXYZ(a);
            if (xyz) indigoSetXYZ(a, xyz[0] + static_cast<float>(dx), xyz[1] + static_cast<float>(dy), xyz[2]);
            indigoFree(a);
        }
        indigoFree(it);
    };
    if (isRxn) {
        int mit = indigoIterateMolecules(mol);
        if (mit >= 0) {
            int m;
            while ((m = indigoNext(mit)) > 0) { shift(m); indigoFree(m); }
            indigoFree(mit);
        }
    } else {
        shift(mol);
    }
}

IndigoService::IndigoService(QObject *parent) : QObject(parent) {
    m_sessionId = indigoAllocSessionId();
}

IndigoService::~IndigoService() {
    indigoSetSessionId(m_sessionId);
    indigoFreeAllObjects();
    indigoReleaseSessionId(m_sessionId);
}

void IndigoService::layout(const QString &molfile) {
    if (molfile.isEmpty()) {
        emit layoutFinished("", "");
        return;
    }
    
    QPointer<IndigoService> self = this;
    
    (void)QtConcurrent::run([self, molfile]() {
        QString resultMol = molfile;
        QString errorMsg;
        unsigned long long threadSessionId = indigoAllocSessionId();
        
        try {
            indigoSetSessionId(threadSessionId);
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                const QPointF oldCenter = indigoCentroid(mol, isRxn);
                if (indigoLayout(mol) >= 0) {
                    const QPointF newCenter = indigoCentroid(mol, isRxn);
                    indigoShiftAll(mol, isRxn, oldCenter.x() - newCenter.x(), oldCenter.y() - newCenter.y());
                    const char* result = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                    if (result) {
                        resultMol = QString::fromUtf8(result);
                    } else {
                        resultMol = "";
                        errorMsg = QString::fromUtf8(indigoGetLastError());
                    }
                } else {
                    qWarning() << "Indigo: Layout failed:" << indigoGetLastError();
                    resultMol = "";
                    errorMsg = QString("Layout failed: %1").arg(indigoGetLastError());
                }
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for layout:" << indigoGetLastError();
                resultMol = "";
                errorMsg = QString("Failed to load molecule for layout: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ Exception in layout:" << e.what();
            resultMol = "";
            errorMsg = QString::fromUtf8(e.what());
        } catch (...) {
            qWarning() << "Indigo C++ Unknown Exception in layout";
            resultMol = "";
            errorMsg = "Unknown exception in layout";
        }
        
        indigoReleaseSessionId(threadSessionId);
        
        emitOnGuiThread(self, [resultMol, errorMsg](IndigoService *s) { emit s->layoutFinished(resultMol, errorMsg); });
    });
}

void IndigoService::clean2d(const QString &molfile) {
    if (molfile.isEmpty()) {
        emit clean2dFinished("");
        return;
    }
    
    QPointer<IndigoService> self = this;
    
    (void)QtConcurrent::run([self, molfile]() {
        QString resultMol = molfile;
        unsigned long long threadSessionId = indigoAllocSessionId();
        
        try {
            indigoSetSessionId(threadSessionId);
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                const QPointF oldCenter = indigoCentroid(mol, isRxn);
                if (indigoClean2d(mol) >= 0) {
                    const QPointF newCenter = indigoCentroid(mol, isRxn);
                    indigoShiftAll(mol, isRxn, oldCenter.x() - newCenter.x(), oldCenter.y() - newCenter.y());
                    const char* result = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                    if (result) {
                        resultMol = QString::fromUtf8(result);
                    }
                } else {
                    qWarning() << "Indigo: Clean2D failed:" << indigoGetLastError();
                }
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for clean2d:" << indigoGetLastError();
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ Exception in clean2d:" << e.what();
        } catch (...) {
            qWarning() << "Indigo C++ Unknown Exception in clean2d";
        }
        
        indigoReleaseSessionId(threadSessionId);
        
        emitOnGuiThread(self, [resultMol](IndigoService *s) { emit s->clean2dFinished(resultMol); });
    });
}

void IndigoService::aromatize(const QString &molfile) {
    if (molfile.isEmpty()) {
        emit aromatizeFinished("");
        return;
    }

    QPointer<IndigoService> self = this;

    (void)QtConcurrent::run([self, molfile]() {
        QString resultMol = molfile;
        unsigned long long threadSessionId = indigoAllocSessionId();

        try {
            indigoSetSessionId(threadSessionId);
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                int ret = indigoAromatize(mol);
                if (ret < 0) {
                    qWarning() << "Indigo: Aromatize failed:" << indigoGetLastError();
                } else {
                    const char* result = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                    if (result) resultMol = QString::fromUtf8(result);
                }
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for aromatize:" << indigoGetLastError();
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ exception in aromatize:" << e.what();
        } catch (...) {
            qWarning() << "Indigo unknown C++ exception in aromatize";
        }

        indigoReleaseSessionId(threadSessionId);

        emitOnGuiThread(self, [resultMol](IndigoService *s) { emit s->aromatizeFinished(resultMol); });
    });
}

void IndigoService::smiles(const QString &molfile) {
    if (molfile.isEmpty()) {
        emit smilesFinished("");
        return;
    }

    QPointer<IndigoService> self = this;

    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long threadSessionId = indigoAllocSessionId();

        try {
            indigoSetSessionId(threadSessionId);
            int mol = isReactionFormat(molfile) ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                                                 : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                // indigoSmiles produces reaction SMILES ("reactants>>products") natively
                // when given a reaction handle — no output-format switch needed here.
                const char* s = indigoSmiles(mol);
                result = s ? QString::fromUtf8(s) : "";
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for SMILES:" << indigoGetLastError();
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ exception in smiles:" << e.what();
        } catch (...) {
            qWarning() << "Indigo unknown C++ exception in smiles";
        }

        indigoReleaseSessionId(threadSessionId);

        emitOnGuiThread(self, [result](IndigoService *s) { emit s->smilesFinished(result); });
    });
}

void IndigoService::dearomatize(const QString &molfile) {
    if (molfile.isEmpty()) { emit dearomatizeFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = molfile;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoDearomatize(mol);
                const char* mf = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                if (mf) result = QString::fromUtf8(mf);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: dearomatize load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->dearomatizeFinished(result); });
    });
}

void IndigoService::unfoldHydrogens(const QString &molfile) {
    if (molfile.isEmpty()) { emit unfoldHydrogensFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = molfile;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoUnfoldHydrogens(mol);
                const char* mf = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                if (mf) result = QString::fromUtf8(mf);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: unfoldHydrogens load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->unfoldHydrogensFinished(result); });
    });
}

void IndigoService::foldHydrogens(const QString &molfile) {
    if (molfile.isEmpty()) { emit foldHydrogensFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = molfile;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoFoldHydrogens(mol);
                const char* mf = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                if (mf) result = QString::fromUtf8(mf);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: foldHydrogens load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->foldHydrogensFinished(result); });
    });
}

void IndigoService::canonicalSmiles(const QString &molfile) {
    if (molfile.isEmpty()) { emit canonicalSmilesFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = isReactionFormat(molfile) ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                                                 : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                // Reaction-aware, same as indigoSmiles above.
                const char* smi = indigoCanonicalSmiles(mol);
                if (smi) result = QString::fromUtf8(smi);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: canonicalSmiles load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->canonicalSmilesFinished(result); });
    });
}

void IndigoService::inchi(const QString &molfile) {
    // InChI has no reaction form (unlike SMILES's native "reactants>>products"),
    // so reactions are guarded out the same way checkStructure guards them below.
    if (molfile.isEmpty() || isReactionFormat(molfile)) { emit inchiFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoInchiInit(sid);
        try {
            int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                const char* inchiStr = indigoInchiGetInchi(mol);
                if (inchiStr) result = QString::fromUtf8(inchiStr);
                else qWarning() << "Indigo: InChI generation failed:" << indigoInchiGetWarning();
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for InChI:" << indigoGetLastError();
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ exception in inchi:" << e.what();
        } catch (...) {
            qWarning() << "Indigo unknown C++ exception in inchi";
        }
        indigoInchiDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->inchiFinished(result); });
    });
}

void IndigoService::inchiKey(const QString &molfile) {
    // InChIKey is a hash of the InChI string, not of the molecule directly --
    // computes the InChI first, then keys it, in one round trip.
    if (molfile.isEmpty() || isReactionFormat(molfile)) { emit inchiKeyFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoInchiInit(sid);
        try {
            int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                const char* inchiStr = indigoInchiGetInchi(mol);
                // Must copy before the next Indigo call: indigoInchiGetInchi's
                // return value points into a buffer indigoInchiGetInchiKey
                // itself reuses/overwrites, so passing the raw pointer straight
                // through as its argument reads already-clobbered memory and
                // silently returns null (confirmed via a standalone repro).
                QString inchiCopy = inchiStr ? QString::fromUtf8(inchiStr) : QString();
                if (!inchiCopy.isEmpty()) {
                    const char* key = indigoInchiGetInchiKey(inchiCopy.toUtf8().constData());
                    if (key) result = QString::fromUtf8(key);
                    else qWarning() << "Indigo: InChIKey generation failed:" << indigoInchiGetWarning();
                } else {
                    qWarning() << "Indigo: InChI generation (for key) failed:" << indigoInchiGetWarning();
                }
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: Failed to load molecule for InChIKey:" << indigoGetLastError();
            }
        } catch (const std::exception& e) {
            qWarning() << "Indigo C++ exception in inchiKey:" << e.what();
        } catch (...) {
            qWarning() << "Indigo unknown C++ exception in inchiKey";
        }
        indigoInchiDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->inchiKeyFinished(result); });
    });
}

void IndigoService::hash(const QString &molfile) {
    if (molfile.isEmpty()) { emit hashFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            bool isRxn = isReactionFormat(molfile);
            int obj = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                             : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (obj >= 0) {
                int64_t h = indigoHash(obj);
                result = QString::number(h);
                indigoFree(obj);
            } else {
                qWarning() << "Indigo: hash load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->hashFinished(result); });
    });
}

void IndigoService::similarity(const QString &molfile, const QString &refSmiles) {
    if (molfile.isEmpty() || refSmiles.isEmpty() || isReactionFormat(molfile)) {
        emit similarityFinished("");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile, refSmiles]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            int ref = indigoLoadMoleculeFromString(refSmiles.toUtf8().constData());
            if (mol >= 0 && ref >= 0) {
                // Both structures must be aromatized before comparing - the similarity
                // fingerprint is sensitive to Kekule-vs-aromatic bond representation, and the
                // two sides here almost never agree on that by default: the canvas structure
                // arrives via a molfile (drawn bonds, typically Kekule), while a user-pasted
                // reference SMILES is very commonly aromatic notation (e.g. "c1ccccc1" from
                // PubChem/ChemDraw). Verified empirically: without this, comparing benzene
                // against itself (aromatic reference) scored 0.0769, not ~1.0; with it, 1.0000.
                if (indigoAromatize(mol) < 0) qWarning() << "Indigo: similarity aromatize(mol) failed:" << indigoGetLastError();
                if (indigoAromatize(ref) < 0) qWarning() << "Indigo: similarity aromatize(ref) failed:" << indigoGetLastError();
                float sim = indigoSimilarity(mol, ref, "tanimoto");
                result = QString::number(sim, 'f', 4);
            } else {
                qWarning() << "Indigo: similarity load failed:" << indigoGetLastError();
            }
            if (mol >= 0) indigoFree(mol);
            if (ref >= 0) indigoFree(ref);
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->similarityFinished(result); });
    });
}

void IndigoService::massComposition(const QString &molfile) {
    if (molfile.isEmpty() || isReactionFormat(molfile)) { emit massCompositionFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(sanitizeMolfileForAnalysis(molfile).toUtf8().constData());
            if (mol >= 0) {
                const char* mc = indigoMassComposition(mol);
                if (mc) result = QString::fromUtf8(mc);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: massComposition load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->massCompositionFinished(result); });
    });
}

void IndigoService::pkaValues(const QString &molfile) {
    if (molfile.isEmpty() || isReactionFormat(molfile)) { emit pkaValuesFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(sanitizeMolfileForAnalysis(molfile).toUtf8().constData());
            if (mol >= 0) {
                const char* pv = indigoPkaValues(mol);
                if (pv) result = QString::fromUtf8(pv);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: pkaValues load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->pkaValuesFinished(result); });
    });
}

void IndigoService::renderToFile(const QString &molfile, const QUrl &fileUrl, const QString &format) {
    if (molfile.isEmpty()) { emit renderFinished(false, "No structure to render."); return; }
    QPointer<IndigoService> self = this;
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) { emit renderFinished(false, "Invalid file path."); return; }
    (void)QtConcurrent::run([self, molfile, path, format]() {
        bool ok = false;
        QString error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoRendererInit(sid);
        try {
            int mol = isReactionFormat(molfile) ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                                                 : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoSetOption("render-output-format", format.toUtf8().constData());
                int res = indigoRenderToFile(mol, path.toUtf8().constData());
                if (res >= 0) ok = true;
                else error = QString::fromUtf8(indigoGetLastError());
                indigoFree(mol);
            } else {
                error = QString("Failed to load structure: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during render.";
        }
        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [ok, error](IndigoService *s) { emit s->renderFinished(ok, error); });
    });
}

void IndigoService::renderReactionGridToFile(const QString &molfile, const QUrl &fileUrl, const QString &format) {
    if (molfile.isEmpty()) { emit renderFinished(false, "No structure to render."); return; }
    if (!isReactionFormat(molfile)) {
        emit renderFinished(false, "Grid export is for reactions only - the active document is a plain molecule.");
        return;
    }
    QPointer<IndigoService> self = this;
    QString path = fileUrl.toLocalFile();
    QString outFormat = QFileInfo(path).suffix().toLower();
    if (outFormat.isEmpty()) outFormat = format;
    if (outFormat != "pdf" && outFormat != "png" && outFormat != "svg") {
        emit renderFinished(false, QString("Unsupported reaction scheme file type \".%1\" - use .pdf, .png or .svg.").arg(outFormat));
        return;
    }
    (void)QtConcurrent::run([self, molfile, path, outFormat]() {
        bool ok = false;
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
                while ((comp = indigoNext(reactants)) > 0) { indigoArrayAdd(arr, comp); count++; }
                int products = indigoIterateProducts(rx);
                while ((comp = indigoNext(products)) > 0) { indigoArrayAdd(arr, comp); count++; }
                if (count > 0) {
                    indigoSetOption("render-output-format", outFormat.toUtf8().constData());
                    int res = indigoRenderGridToFile(arr, nullptr, count, path.toUtf8().constData());
                    if (res >= 0) ok = true;
                    else error = QString::fromUtf8(indigoGetLastError());
                } else {
                    error = "Reaction has no reactants or products to render.";
                }
                indigoFree(arr);
                indigoFree(rx);
            } else {
                error = QString("Failed to load reaction: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during grid render.";
        }
        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [ok, error](IndigoService *s) { emit s->renderFinished(ok, error); });
    });
}

void IndigoService::exportBatchGridToFile(const QStringList &molfiles, const QUrl &fileUrl, const QString &format) {
    if (molfiles.isEmpty()) { emit renderFinished(false, "No structures to export."); return; }
    QPointer<IndigoService> self = this;
    QString path = fileUrl.toLocalFile();
    (void)QtConcurrent::run([self, molfiles, path, format]() {
        bool ok = false;
        QString error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoRendererInit(sid);
        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;
            try {
                for (const QString &mf : molfiles) {
                    int mol = indigoLoadMoleculeFromString(mf.toUtf8().constData());
                    if (mol >= 0) {
                        indigoLayout(mol);
                        indigoArrayAdd(arr, mol);
                        loadedMols.append(mol);
                    }
                }
            } catch (...) {
                for (int m : loadedMols) indigoFree(m);
                throw;
            }
            if (!loadedMols.isEmpty()) {
                int nColumns = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(loadedMols.size()))));
                indigoSetOption("render-output-format", format.toUtf8().constData());
                int res = indigoRenderGridToFile(arr, nullptr, nColumns, path.toUtf8().constData());
                if (res >= 0) ok = true;
                else error = QString::fromUtf8(indigoGetLastError());
            } else {
                error = "No valid structures could be parsed from this batch.";
            }
            for (int m : loadedMols) indigoFree(m);
            indigoFree(arr);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during batch grid export.";
        }
        indigoRendererDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [ok, error](IndigoService *s) { emit s->renderFinished(ok, error); });
    });
}

void IndigoService::parseRdfBatch(const QUrl &fileUrl) {
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) { emit rdfBatchParsed("", "Invalid file path."); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, path]() {
        QString error;
        QJsonArray records;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int reader = indigoIterateRDFile(path.toUtf8().constData());
            if (reader >= 0) {
                int item;
                int count = 0;
                const int limit = 500;  // same cap the SDF batch path already uses
                while (count < limit && (item = indigoNext(reader)) > 0) {
                    bool isRxn = indigoCountReactants(item) >= 0;
                    const char* text = isRxn ? indigoRxnfile(item) : indigoMolfile(item);
                    if (text) {
                        QJsonObject rec;
                        // Copy immediately via QString::fromUtf8 - text is a pointer into an
                        // internal buffer invalidated by the next Indigo call (verified: this
                        // exact mistake produced "scanner: BufferScanner::read() error" when
                        // the copy was deferred past a subsequent Indigo call).
                        rec["molfile"] = QString::fromUtf8(text);
                        rec["label"] = QString("Record %1").arg(count + 1);
                        records.append(rec);
                    }
                    indigoFree(item);
                    count++;
                }
                indigoFree(reader);
            } else {
                error = QString("Failed to open RDF file: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception while parsing RDF file.";
        }
        indigoReleaseSessionId(sid);
        QString json;
        if (error.isEmpty()) json = QString::fromUtf8(QJsonDocument(records).toJson(QJsonDocument::Compact));
        emitOnGuiThread(self, [json, error](IndigoService *s) { emit s->rdfBatchParsed(json, error); });
    });
}

void IndigoService::parseIndigoBatchFile(const QUrl &fileUrl, const QString &format) {
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) { emit indigoBatchParsed("", "Invalid file path."); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, path, format]() {
        QString error;
        QJsonArray records;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int reader = -1;
            if (format == "smiles") reader = indigoIterateSmilesFile(path.toUtf8().constData());
            else if (format == "cml") reader = indigoIterateCMLFile(path.toUtf8().constData());
            else if (format == "cdx") reader = indigoIterateCDXFile(path.toUtf8().constData());
            if (reader >= 0) {
                int item;
                int count = 0;
                const int limit = 500;
                while (count < limit && (item = indigoNext(reader)) > 0) {
                    // SMILES (and some CDX/CML sources) carry no 2D coordinates - indigoLayout()
                    // must run before serializing, or every atom collapses onto (0,0) and the
                    // batch picker's thumbnail renders as a blank point (same bug class as the
                    // indigoExtractCommonScaffold collapse fixed earlier for scaffold detection).
                    indigoLayout(item);
                    bool isRxn = indigoCountReactants(item) >= 0;
                    const char* text = isRxn ? indigoRxnfile(item) : indigoMolfile(item);
                    if (text) {
                        QJsonObject rec;
                        rec["molfile"] = QString::fromUtf8(text);
                        rec["label"] = QString("Record %1").arg(count + 1);
                        records.append(rec);
                    }
                    indigoFree(item);
                    count++;
                }
                indigoFree(reader);
            } else {
                error = QString("Failed to open %1 file: %2").arg(format, QString::fromUtf8(indigoGetLastError()));
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception while parsing batch file.";
        }
        indigoReleaseSessionId(sid);
        QString json;
        if (error.isEmpty()) json = QString::fromUtf8(QJsonDocument(records).toJson(QJsonDocument::Compact));
        emitOnGuiThread(self, [json, error](IndigoService *s) { emit s->indigoBatchParsed(json, error); });
    });
}

void IndigoService::exportBatchToFile(const QStringList &molfiles, const QUrl &fileUrl, const QString &format) {
    if (molfiles.isEmpty()) { emit renderFinished(false, "No structures to export."); return; }
    QPointer<IndigoService> self = this;
    QString path = fileUrl.toLocalFile();
    (void)QtConcurrent::run([self, molfiles, path, format]() {
        bool ok = false;
        QString error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int buffer = indigoWriteBuffer();
            int saver = buffer >= 0 ? indigoCreateSaver(buffer, format.toUtf8().constData()) : -1;
            if (saver >= 0) {
                int written = 0;
                int recordNum = 0;
                for (const QString &mf : molfiles) {
                    recordNum++;
                    bool isRxn = isReactionFormat(mf);
                    int obj = isRxn ? indigoLoadReactionFromString(mf.toUtf8().constData())
                                    : indigoLoadMoleculeFromString(mf.toUtf8().constData());
                    if (obj >= 0) {
                        // A blank molecule-name header line makes chem-core.js's own
                        // MolSerializer.deserialize() reject the record on read-back
                        // (badHeaderRecover defaults to false) - always give every
                        // record a real name before writing so exported files stay
                        // fully round-trippable by this app's own SDF/RDF reader.
                        const char* existingName = indigoName(obj);
                        if (!existingName || existingName[0] == '\0') {
                            indigoSetName(obj, QString("Record %1").arg(recordNum).toUtf8().constData());
                        }
                        if (indigoAppend(saver, obj) >= 0) written++;
                        indigoFree(obj);
                    }
                }
                indigoClose(saver);
                indigoFree(saver);
                char *raw = nullptr;
                int size = 0;
                if (written == 0) {
                    error = "No valid structures could be written.";
                } else if (indigoToBuffer(buffer, &raw, &size) < 0) {
                    error = QString::fromUtf8(indigoGetLastError());
                } else {
                    QSaveFile out(path);
                    if (out.open(QIODevice::WriteOnly) && out.write(raw, size) == size && out.commit()) ok = true;
                    else error = QString("Failed to write %1: %2").arg(path, out.errorString());
                }
            } else {
                error = QString("Failed to create %1 saver: %2").arg(format, QString::fromUtf8(indigoGetLastError()));
            }
            if (buffer >= 0) indigoFree(buffer);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during batch export.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [ok, error](IndigoService *s) { emit s->renderFinished(ok, error); });
    });
}


void IndigoService::autoMapReaction(const QString &molfile) {
    if (molfile.isEmpty()) { emit reactionMappingFinished("", "No structure to map."); return; }
    if (!isReactionFormat(molfile)) {
        emit reactionMappingFinished("", "Atom mapping is for reactions only - the active document is a plain molecule.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int rx = indigoLoadReactionFromString(molfile.toUtf8().constData());
            if (rx >= 0) {
                if (indigoAutomap(rx, "discard") >= 0) {
                    const char* rf = indigoRxnfile(rx);
                    if (rf) result = QString::fromUtf8(rf);
                    else error = QString::fromUtf8(indigoGetLastError());
                } else {
                    error = QString("Auto-mapping failed: %1").arg(indigoGetLastError());
                }
                indigoFree(rx);
            } else {
                error = QString("Failed to load reaction: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during atom mapping.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->reactionMappingFinished(result, error); });
    });
}

void IndigoService::clearReactionMapping(const QString &molfile) {
    if (molfile.isEmpty()) { emit reactionMappingFinished("", "No structure to clear mapping from."); return; }
    if (!isReactionFormat(molfile)) {
        emit reactionMappingFinished("", "Atom mapping is for reactions only - the active document is a plain molecule.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int rx = indigoLoadReactionFromString(molfile.toUtf8().constData());
            if (rx >= 0) {
                if (indigoClearAAM(rx) >= 0) {
                    const char* rf = indigoRxnfile(rx);
                    if (rf) result = QString::fromUtf8(rf);
                    else error = QString::fromUtf8(indigoGetLastError());
                } else {
                    error = QString("Clear mapping failed: %1").arg(indigoGetLastError());
                }
                indigoFree(rx);
            } else {
                error = QString("Failed to load reaction: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception clearing atom mapping.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->reactionMappingFinished(result, error); });
    });
}

void IndigoService::correctReactingCenters(const QString &molfile) {
    if (molfile.isEmpty()) { emit reactionMappingFinished("", "No structure to analyze."); return; }
    if (!isReactionFormat(molfile)) {
        emit reactionMappingFinished("", "Reacting centers are for reactions only - the active document is a plain molecule.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int rx = indigoLoadReactionFromString(molfile.toUtf8().constData());
            if (rx >= 0) {
                if (indigoCorrectReactingCenters(rx) >= 0) {
                    const char* rf = indigoRxnfile(rx);
                    if (rf) result = QString::fromUtf8(rf);
                    else error = QString::fromUtf8(indigoGetLastError());
                } else {
                    error = QString("Reacting-center correction failed: %1").arg(indigoGetLastError());
                }
                indigoFree(rx);
            } else {
                error = QString("Failed to load reaction: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during reacting-center correction.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->reactionMappingFinished(result, error); });
    });
}

void IndigoService::normalize(const QString &molfile) {
    if (molfile.isEmpty()) { emit normalizeFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = molfile;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoNormalize(mol, "");
                const char* mf = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                if (mf) result = QString::fromUtf8(mf);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: normalize load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->normalizeFinished(result); });
    });
}

void IndigoService::standardize(const QString &molfile) {
    if (molfile.isEmpty()) { emit standardizeFinished(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = molfile;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int mol = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                indigoStandardize(mol);
                const char* mf = isRxn ? indigoRxnfile(mol) : indigoMolfile(mol);
                if (mf) result = QString::fromUtf8(mf);
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: standardize load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->standardizeFinished(result); });
    });
}

void IndigoService::ionizeAtPh(const QString &molfile, double pH) {
    if (molfile.isEmpty()) { emit ionizeFinished("", "No structure to ionize."); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile, pH]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            const bool isRxn = isReactionFormat(molfile);
            int obj = isRxn ? indigoLoadReactionFromString(molfile.toUtf8().constData())
                            : indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (obj >= 0) {
                if (indigoIonize(obj, static_cast<float>(pH), 1.0f) >= 0) {
                    const char* out = isRxn ? indigoRxnfile(obj) : indigoMolfile(obj);
                    if (out) result = QString::fromUtf8(out);
                    else error = QString::fromUtf8(indigoGetLastError());
                } else {
                    error = QString("Ionization failed: %1").arg(indigoGetLastError());
                }
                indigoFree(obj);
            } else {
                error = QString("Failed to load structure: %1").arg(indigoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during ionization.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->ionizeFinished(result, error); });
    });
}

void IndigoService::calcProperties(const QString &molfile) {
    if (molfile.isEmpty() || isReactionFormat(molfile)) {
        emit propertiesReady(0, 0, "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, 0, 0, 0);
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        double mw = 0, mono = 0, tpsa = 0, logp = 0, molarRefractivity = 0, pka = 0, mostAbundantMass = 0;
        QString mf;
        int atoms = 0, bonds = 0, hba = 0, hbd = 0, rotBonds = 0;
        int heavyAtoms = 0;
        bool isChiral = false;
        int fragmentCount = 0;
        int ringCount = 0;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(sanitizeMolfileForAnalysis(molfile).toUtf8().constData());
            if (mol >= 0) {
                double mwCalc = indigoMolecularWeight(mol);
                if (!std::isnan(mwCalc)) mw = mwCalc;
                double monoCalc = indigoMonoisotopicMass(mol);
                if (!std::isnan(monoCalc)) mono = monoCalc;
                int atomsCalc = indigoCountAtoms(mol);
                if (atomsCalc >= 0) atoms = atomsCalc;
                int bondsCalc = indigoCountBonds(mol);
                if (bondsCalc >= 0) bonds = bondsCalc;
                double tpsaCalc = indigoTPSA(mol, 1);
                if (!std::isnan(tpsaCalc)) tpsa = tpsaCalc;
                double logpCalc = indigoLogP(mol);
                if (!std::isnan(logpCalc)) logp = logpCalc;
                int hbaCalc = indigoNumHydrogenBondAcceptors(mol);
                if (hbaCalc >= 0) hba = hbaCalc;
                int hbdCalc = indigoNumHydrogenBondDonors(mol);
                if (hbdCalc >= 0) hbd = hbdCalc;
                int rotBondsCalc = indigoNumRotatableBonds(mol);
                if (rotBondsCalc >= 0) rotBonds = rotBondsCalc;
                double molarRefractivityCalc = indigoMolarRefractivity(mol);
                if (!std::isnan(molarRefractivityCalc)) molarRefractivity = molarRefractivityCalc;
                double pkaCalc = indigoPka(mol);
                if (!std::isnan(pkaCalc)) pka = pkaCalc;
                int heavyAtomsCalc = indigoCountHeavyAtoms(mol);
                if (heavyAtomsCalc >= 0) heavyAtoms = heavyAtomsCalc;
                isChiral = indigoIsChiral(mol) != 0;
                double mostAbundantMassCalc = indigoMostAbundantMass(mol);
                if (!std::isnan(mostAbundantMassCalc)) mostAbundantMass = mostAbundantMassCalc;
                int fragmentCountCalc = indigoCountComponents(mol);
                if (fragmentCountCalc >= 0) fragmentCount = fragmentCountCalc;
                int ringCountCalc = indigoCountSSSR(mol);
                if (ringCountCalc >= 0) ringCount = ringCountCalc;
                int fmHandle = indigoMolecularFormula(mol);
                if (fmHandle >= 0) {
                    const char* fmStr = indigoToString(fmHandle);
                    if (fmStr) mf = QString::fromUtf8(fmStr);
                    indigoFree(fmHandle);
                }
                indigoFree(mol);
            } else {
                qWarning() << "Indigo: calcProperties load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [=](IndigoService *s) { emit s->propertiesReady(mw, mono, mf, atoms, bonds, tpsa, logp, hba, hbd, rotBonds, molarRefractivity, pka, heavyAtoms, isChiral, mostAbundantMass, fragmentCount, ringCount); });
    });
}

void IndigoService::calcStereoDescriptors(const QString &molfile) {
    if (molfile.isEmpty() || isReactionFormat(molfile)) {
        emit stereoDescriptorsReady("{\"atoms\":{},\"bonds\":{}}");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result = "{\"atoms\":{},\"bonds\":{}}";
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(sanitizeMolfileForAnalysis(molfile).toUtf8().constData());
            if (mol >= 0) {
                indigoAddCIPStereoDescriptors(mol);
                QJsonObject atomMap;
                // R/S (and rare r/s pseudo-asymmetric) atom descriptors: these are
                // reachable directly via the atom-stereocenter iterator.
                int iter = indigoIterateStereocenters(mol);
                if (iter >= 0) {
                    int atom;
                    while ((atom = indigoNext(iter)) != 0) {
                        if (atom == -1) break;
                        int idx = indigoIndex(atom);
                        int cip = indigoStereocenterCIPDescriptor(atom);
                        int stereoType = indigoStereocenterType(atom);
                        int group = indigoStereocenterGroup(atom);
                        QString label;
                        switch (cip) {
                            case 2: label = "s"; break;
                            case 3: label = "r"; break;
                            case 4: label = "S"; break;
                            case 5: label = "R"; break;
                            // Note: cip values 6/7 (E/Z) never occur here - Indigo
                            // stores double-bond CIP descriptors in a separate,
                            // bond-indexed map that indigoStereocenterCIPDescriptor
                            // (atom-only) cannot read. See the sgroup-based pass below.
                            default: label = ""; break;
                        }
                        if (!label.isEmpty()) {
                            QJsonObject entry;
                            entry["cipLabel"] = label;
                            entry["type"] = stereoType;
                            entry["group"] = group;
                            atomMap[QString::number(idx)] = entry;
                        }
                        indigoFree(atom);
                    }
                    indigoFree(iter);
                }

                // E/Z bond descriptors: read via the same direct KET JSON "cip" bond-field
                // helper IupacNamer.cpp's name generation already uses successfully. The
                // previous DAT-sgroup search here was dead code -- indigoAddCIPSgroups
                // (which emits "INDIGO_CIP_DESC" DAT sgroups) is only ever invoked from
                // Indigo's molfile-save path, never from indigoJson/KET, so it never
                // populated anything for this app's actual call pattern (confirmed
                // 2026-08-08 during the IUPAC namer's P-93.4 CIP-reuse fix; DAT-sgroup
                // hunt removed 2026-08-18 -- see IUPAC Blue Book Coverage.md item 8).
                QJsonObject bondMap;
                for (const auto &entry : computeIndigoBondCIP(mol)) {
                    QString key = QString("%1-%2").arg(entry.first.first).arg(entry.first.second);
                    QJsonObject bondEntry;
                    bondEntry["cipLabel"] = QString(entry.second);
                    bondMap[key] = bondEntry;
                }

                indigoFree(mol);
                QJsonObject combined;
                combined["atoms"] = atomMap;
                combined["bonds"] = bondMap;
                result = QString::fromUtf8(QJsonDocument(combined).toJson(QJsonDocument::Compact));
            } else {
                qWarning() << "Indigo: calcStereoDescriptors load failed:" << indigoGetLastError();
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->stereoDescriptorsReady(result); });
    });
}

// Parses indigoCheckObj's report shape ({"checkType": "...(id1,id2,...)"}) into
// structured issues. Shared by checkStructure's molecule and per-reaction-
// component branches, since indigoCheckObj gives back this exact same shape
// for a lone molecule handle or a single reaction-component handle alike.
static QJsonArray parseCheckReport(const QString &rawReport) {
    QJsonArray issuesArray;
    QJsonParseError parseErr;
    QJsonDocument doc = QJsonDocument::fromJson(rawReport.toUtf8(), &parseErr);
    if (parseErr.error != QJsonParseError::NoError || !doc.isObject()) return issuesArray;
    static const QSet<QString> atomChecks = {
        "valence", "radical", "pseudoatom", "stereo",
        "ambiguous_h", "3d_coord", "overlap_atom"
    };
    static const QRegularExpression idRe("\\((\\d+(?:,\\d+)*)\\)\\s*$");
    // doc.object() must be held in a named variable — QJsonObject::const_iterator
    // stores a raw pointer back to the QJsonObject it came from, so calling
    // doc.object() separately for constBegin()/constEnd() (as this used to) creates
    // two temporaries; the one behind the begin iterator is destroyed at the end of
    // the for-loop's init-statement, leaving a dangling iterator used for the rest
    // of the loop.
    const QJsonObject checkObj = doc.object();
    for (auto it = checkObj.constBegin(); it != checkObj.constEnd(); ++it) {
        QString checkType = it.key();
        QString value = it.value().toString();
        auto match = idRe.match(value);
        if (match.hasMatch()) {
            QStringList idStrs = match.captured(1).split(',');
            QJsonArray ids;
            for (const QString& s : idStrs) {
                bool ok = false;
                int idx = s.trimmed().toInt(&ok);
                if (ok && idx >= 1) ids.append(idx - 1); // normalize 1-based to 0-based
            }
            if (!ids.isEmpty()) {
                QJsonObject issue;
                issue["type"] = checkType;
                issue["target"] = atomChecks.contains(checkType) ? "atom" : "bond";
                issue["ids"] = ids;
                issuesArray.append(issue);
            }
        }
    }
    return issuesArray;
}

void IndigoService::checkStructure(const QString &molfile) {
    if (molfile.isEmpty()) { emit checkFinished("No structure to check."); emit checkIssuesReady("{\"issues\":[]}"); return; }
    bool isRxn = isReactionFormat(molfile);
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile, isRxn]() {
        QString report;
        QString structured;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            if (isRxn) {
                // indigoCheckObj on the whole reaction handle gives back a report shape
                // (empty-string key, no component association) parseCheckReport can't
                // use -- checking each reactant/product component individually gives the
                // same clean per-type shape the molecule branch already parses, just with
                // ids local to that one component. Inline atom/bond highlighting
                // (checkIssuesReady below) deliberately stays off for reactions: mapping
                // those per-component-local ids back to global canvas ids would need
                // matching chem-core.js's own reactant/product fragment ordering when it
                // serializes a $RXN string, which isn't a safely verifiable assumption
                // here -- a wrong mapping would silently highlight the wrong atom.
                int rx = indigoLoadReactionFromString(molfile.toUtf8().constData());
                if (rx >= 0) {
                    // Same "show the raw per-component check JSON" style the molecule
                    // branch below already uses (it doesn't reformat indigoCheckObj's
                    // text into prose either) -- just labeled by component here.
                    QStringList lines;
                    auto checkComponents = [&](int iterHandle, const char* label) {
                        int idx = 0;
                        int comp;
                        while ((comp = indigoNext(iterHandle)) > 0) {
                            idx++;
                            const char* res = indigoCheckObj(comp, "");
                            QString rawIssue = res ? QString::fromUtf8(res).trimmed() : QString();
                            if (!rawIssue.isEmpty() && rawIssue != "{}") {
                                lines << QString("%1 %2: %3").arg(label).arg(idx).arg(rawIssue);
                            }
                        }
                    };
                    checkComponents(indigoIterateReactants(rx), "Reactant");
                    checkComponents(indigoIterateProducts(rx), "Product");
                    report = lines.isEmpty() ? "No problems found." : lines.join("\n");
                    indigoFree(rx);
                } else {
                    report = QString("Parse error: %1").arg(indigoGetLastError());
                }
                structured = "{\"issues\":[]}";
            } else {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* res = indigoCheckObj(mol, "");
                    report = res ? QString::fromUtf8(res).trimmed() : QString();
                    QJsonParseError perr;
                    QJsonDocument reportDoc = QJsonDocument::fromJson(report.toUtf8(), &perr);
                    QJsonObject reportObj = (perr.error == QJsonParseError::NoError && reportDoc.isObject())
                        ? reportDoc.object() : QJsonObject();
                    // indigoCheckChirality returns 0 when inconsistent (on error), 1 when consistent (or unset)
                    if (indigoCheckChirality(mol) == 0) {
                        reportObj["chirality"] = "Chiral flag is inconsistent with the drawn stereocenters -- "
                            "either the structure is marked chiral with no stereocenters defined, or the chiral "
                            "flag itself is unset/ambiguous. Check the CHIRAL setting and any wedge/hash bonds.";
                    }
                    if (indigoCheckStereo(mol) > 0) {
                        reportObj["stereocenters"] = "One or more marked stereocenters are inconsistent with the "
                            "molecule's own symmetry (invalid or redundant stereo descriptor) -- re-check the "
                            "wedge/hash bonds around them.";
                    }
                    report = QString::fromUtf8(QJsonDocument(reportObj).toJson(QJsonDocument::Compact));
                    QJsonObject structuredObj;
                    structuredObj["issues"] = parseCheckReport(report);
                    structured = QString::fromUtf8(
                        QJsonDocument(structuredObj).toJson(QJsonDocument::Compact));
                    indigoFree(mol);
                } else {
                    report = QString("Parse error: %1").arg(indigoGetLastError());
                }
            }
        } catch (const std::exception &e) {
            report = QString::fromUtf8(e.what());
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [report, structured](IndigoService *s) {
            emit s->checkFinished(report.isEmpty() ? "{}" : report);
            emit s->checkIssuesReady(structured.isEmpty() ? "{\"issues\":[]}" : structured);
        });
    });
}

void IndigoService::generateIupacName(const QString &molfile) {
    if (molfile.isEmpty()) {
        emit iupacNameReady("", "No structure provided.");
        return;
    }
    if (isReactionFormat(molfile)) {
        emit iupacNameReady("", "IUPAC name generation is not supported for reactions.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString name, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        indigoInchiInit(sid);
        try {
            int mol = indigoLoadMoleculeFromString(sanitizeMolfileForAnalysis(molfile).toUtf8().constData());
            if (mol >= 0) {
                IupacResult res = IupacNamer::generateName(mol);
                if (res.success) {
                    name = res.name;
                } else {
                    error = res.error;
                }
                indigoFree(mol);
            } else {
                error = QString("Failed to load structure: %1").arg(QString::fromUtf8(indigoGetLastError()));
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during IUPAC name generation.";
        }
        indigoInchiDispose(sid);
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [name, error](IndigoService *s) {
            emit s->iupacNameReady(name, error);
        });
    });
}

void IndigoService::copyToClipboard(const QString &text) {
    QGuiApplication::clipboard()->setText(text);
}

void IndigoService::substructureSearch(const QString &molfile, const QString &smarts) {
    if (molfile.isEmpty() || smarts.isEmpty() || isReactionFormat(molfile)) {
        emit substructureSearchFinished("{\"error\":\"No structure or empty pattern\"}");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile, smarts]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
            if (mol >= 0) {
                int query = indigoLoadSmartsFromString(smarts.toUtf8().constData());
                if (query >= 0) {
                    int matcher = indigoSubstructureMatcher(mol, "");
                    if (matcher >= 0) {
                        const int kMaxMatches = 200;
                        QJsonArray matchesArr;
                        int matchIter = indigoIterateMatches(matcher, query);
                        if (matchIter >= 0) {
                            int mapping;
                            while (matchesArr.size() < kMaxMatches && (mapping = indigoNext(matchIter)) != 0) {
                                if (mapping == -1) break;
                                QJsonArray oneMatch;
                                int qAtomIter = indigoIterateAtoms(query);
                                if (qAtomIter >= 0) {
                                    int qAtom;
                                    while ((qAtom = indigoNext(qAtomIter)) != 0) {
                                        if (qAtom == -1) break;
                                        int tAtom = indigoMapAtom(mapping, qAtom);
                                        if (tAtom >= 0) { oneMatch.append(indigoIndex(tAtom) + 1); indigoFree(tAtom); }
                                        indigoFree(qAtom);
                                    }
                                    indigoFree(qAtomIter);
                                }
                                matchesArr.append(oneMatch);
                                indigoFree(mapping);
                            }
                            indigoFree(matchIter);
                        }
                        QJsonObject obj;
                        obj["matchCount"] = matchesArr.size();
                        obj["truncated"] = (matchesArr.size() >= kMaxMatches);
                        obj["matches"] = matchesArr;
                        result = QString::fromUtf8(QJsonDocument(obj).toJson(QJsonDocument::Compact));
                        indigoFree(matcher);
                    } else {
                        result = "{\"error\":\"Failed to build substructure matcher\"}";
                    }
                    indigoFree(query);
                } else {
                    QString err = QString::fromUtf8(indigoGetLastError()).replace("\"", "'");
                    result = QString("{\"error\":\"Invalid SMARTS: %1\"}").arg(err);
                }
                indigoFree(mol);
            } else {
                result = "{\"error\":\"Failed to load structure\"}";
            }
        } catch (const std::exception& e) {
            result = QString("{\"error\":\"%1\"}").arg(QString::fromUtf8(e.what()).replace("\"", "'"));
        } catch (...) {
            result = "{\"error\":\"Unknown exception\"}";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->substructureSearchFinished(result); });
    });
}

// ── Biopolymer loading ────────────────────────────────────────────────────────
// Each function: parse notation → expand monomers to atoms → 2D layout → molfile
void IndigoService::loadBioSequence(const QString &text, const QString &seqType) {
    if (text.isEmpty()) return;
    QString badChars;
    if (!validateBioSequenceAlphabet(text, seqType, /*isFasta=*/false, badChars)) {
        emit biopolymerLoadError(QString("SEQUENCE loader: Invalid symbols in the sequence: %1").arg(badChars));
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, text, seqType]() {
        QString result; QString err;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib < 0) { err = "Monomer library could not be loaded."; }
            else {
                int mol = indigoLoadSequenceFromString(
                    text.toUtf8().constData(), seqType.toUtf8().constData(), lib);
                if (mol >= 0) result = indigoBioExpand(mol);
                else err = QString::fromUtf8(indigoGetLastError());
                indigoFree(lib);
            }
        } catch (const std::exception &e) { err = QString::fromUtf8(e.what()); }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, err](IndigoService *s) {
            if (!result.isEmpty()) emit s->biopolymerLoaded(result);
            else emit s->biopolymerLoadError(err.isEmpty() ? "Failed to parse sequence." : err);
        });
    });
}

void IndigoService::loadBioFasta(const QString &text, const QString &seqType) {
    if (text.isEmpty()) return;
    QString badChars;
    if (!validateBioSequenceAlphabet(text, seqType, /*isFasta=*/true, badChars)) {
        emit biopolymerLoadError(QString("FASTA loader: Invalid symbols in the sequence: %1").arg(badChars));
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, text, seqType]() {
        QString result; QString err;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib < 0) { err = "Monomer library could not be loaded."; }
            else {
                int mol = indigoLoadFastaFromString(
                    text.toUtf8().constData(), seqType.toUtf8().constData(), lib);
                if (mol >= 0) result = indigoBioExpand(mol);
                else err = QString::fromUtf8(indigoGetLastError());
                indigoFree(lib);
            }
        } catch (const std::exception &e) { err = QString::fromUtf8(e.what()); }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, err](IndigoService *s) {
            if (!result.isEmpty()) emit s->biopolymerLoaded(result);
            else emit s->biopolymerLoadError(err.isEmpty() ? "Failed to parse FASTA." : err);
        });
    });
}

void IndigoService::loadBioHelm(const QString &text) {
    if (text.isEmpty()) return;
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, text]() {
        QString result; QString err;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib < 0) { err = "Monomer library could not be loaded."; }
            else {
                int mol = indigoLoadHelmFromString(text.toUtf8().constData(), lib);
                if (mol >= 0) result = indigoBioExpand(mol);
                else err = QString::fromUtf8(indigoGetLastError());
                indigoFree(lib);
            }
        } catch (const std::exception &e) { err = QString::fromUtf8(e.what()); }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, err](IndigoService *s) {
            if (!result.isEmpty()) emit s->biopolymerLoaded(result);
            else emit s->biopolymerLoadError(err.isEmpty() ? "Failed to parse HELM." : err);
        });
    });
}

void IndigoService::loadBioIdt(const QString &text) {
    if (text.isEmpty()) return;
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, text]() {
        QString result; QString err;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib < 0) { err = "Monomer library could not be loaded."; }
            else {
                int mol = indigoLoadIdtFromString(text.toUtf8().constData(), lib);
                if (mol >= 0) result = indigoBioExpand(mol);
                else err = QString::fromUtf8(indigoGetLastError());
                indigoFree(lib);
            }
        } catch (const std::exception &e) { err = QString::fromUtf8(e.what()); }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, err](IndigoService *s) {
            if (!result.isEmpty()) emit s->biopolymerLoaded(result);
            else emit s->biopolymerLoadError(err.isEmpty() ? "Failed to parse IDT." : err);
        });
    });
}

void IndigoService::loadBioAxoLabs(const QString &text) {
    if (text.isEmpty()) return;
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, text]() {
        QString result; QString err;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib < 0) { err = "Monomer library could not be loaded."; }
            else {
                int mol = indigoLoadAxoLabsFromString(text.toUtf8().constData(), lib);
                if (mol >= 0) result = indigoBioExpand(mol);
                else err = QString::fromUtf8(indigoGetLastError());
                indigoFree(lib);
            }
        } catch (const std::exception &e) { err = QString::fromUtf8(e.what()); }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, err](IndigoService *s) {
            if (!result.isEmpty()) emit s->biopolymerLoaded(result);
            else emit s->biopolymerLoadError(err.isEmpty() ? "Failed to parse AxoLabs." : err);
        });
    });
}

// ── Biopolymer export ─────────────────────────────────────────────────────────
// Each function: load molfile → call Indigo serialiser → emit ready signal

void IndigoService::exportBioSequence(const QString &molfile) {
    if (molfile.isEmpty()) { emit bioSequenceReady(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib >= 0) {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* s = indigoSequence(mol, lib);
                    if (s) result = QString::fromUtf8(s);
                    indigoFree(mol);
                }
                indigoFree(lib);
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->bioSequenceReady(result); });
    });
}

void IndigoService::exportBioFasta(const QString &molfile) {
    if (molfile.isEmpty()) { emit bioFastaReady(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib >= 0) {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* s = indigoFasta(mol, lib);
                    if (s) result = QString::fromUtf8(s);
                    indigoFree(mol);
                }
                indigoFree(lib);
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->bioFastaReady(result); });
    });
}

void IndigoService::exportBioHelm(const QString &molfile) {
    if (molfile.isEmpty()) { emit bioHelmReady(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib >= 0) {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* s = indigoHelm(mol, lib);
                    if (s) result = QString::fromUtf8(s);
                    indigoFree(mol);
                }
                indigoFree(lib);
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->bioHelmReady(result); });
    });
}

void IndigoService::exportBioIdt(const QString &molfile) {
    if (molfile.isEmpty()) { emit bioIdtReady(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib >= 0) {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* s = indigoIdt(mol, lib);
                    if (s) result = QString::fromUtf8(s);
                    indigoFree(mol);
                }
                indigoFree(lib);
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->bioIdtReady(result); });
    });
}

void IndigoService::exportBioAxoLabs(const QString &molfile) {
    if (molfile.isEmpty()) { emit bioAxoLabsReady(""); return; }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfile]() {
        QString result;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int lib = loadMonomerLibrary();
            if (lib >= 0) {
                int mol = indigoLoadMoleculeFromString(molfile.toUtf8().constData());
                if (mol >= 0) {
                    const char* s = indigoAxoLabs(mol, lib);
                    if (s) result = QString::fromUtf8(s);
                    indigoFree(mol);
                }
                indigoFree(lib);
            }
        } catch (...) {}
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result](IndigoService *s) { emit s->bioAxoLabsReady(result); });
    });
}

void IndigoService::findCommonScaffold(const QStringList &molfiles) {
    if (molfiles.size() < 2) {
        emit commonScaffoldFinished("", "Need at least 2 structures to find a common scaffold.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfiles]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;
            for (const QString &mf : molfiles) {
                int mol = indigoLoadMoleculeFromString(mf.toUtf8().constData());
                if (mol >= 0) {
                    if (indigoAromatize(mol) < 0) qWarning() << "Indigo: findCommonScaffold aromatize failed:" << indigoGetLastError();
                    indigoArrayAdd(arr, mol);
                    loadedMols.append(mol);
                }
            }
            if (loadedMols.size() >= 2) {
                int scaffold = indigoExtractCommonScaffold(arr, "");
                if (scaffold > 0) {
                    // indigoExtractCommonScaffold never assigns 2D coordinates on its own -
                    // without an explicit layout, every atom lands at (0,0,0) and the result
                    // renders as a single collapsed point. Verified empirically.
                    if (indigoLayout(scaffold) < 0) qWarning() << "Indigo: findCommonScaffold layout failed:" << indigoGetLastError();
                    const char* mf = indigoMolfile(scaffold);
                    if (mf) result = QString::fromUtf8(mf);
                    else error = QString::fromUtf8(indigoGetLastError());
                    indigoFree(scaffold);
                } else {
                    error = "No common scaffold found across these structures.";
                }
            } else {
                error = "Fewer than 2 structures could be parsed.";
            }
            for (int m : loadedMols) indigoFree(m);
            indigoFree(arr);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during scaffold extraction.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->commonScaffoldFinished(result, error); });
    });
}

void IndigoService::decomposeToRGroups(const QStringList &molfiles) {
    if (molfiles.size() < 2) {
        emit rgroupDecompositionFinished("", "Need at least 2 structures to decompose.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfiles]() {
        QString result, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;
            for (const QString &mf : molfiles) {
                int mol = indigoLoadMoleculeFromString(mf.toUtf8().constData());
                if (mol >= 0) {
                    // Same aromatize-before-matching fix as findCommonScaffold/similarity().
                    if (indigoAromatize(mol) < 0) qWarning() << "Indigo: decomposeToRGroups aromatize failed:" << indigoGetLastError();
                    indigoArrayAdd(arr, mol);
                    loadedMols.append(mol);
                }
            }
            if (loadedMols.size() >= 2) {
                int scaffold = indigoExtractCommonScaffold(arr, "");
                if (scaffold > 0) {
                    int decomp = indigoDecomposeMolecules(scaffold, arr);
                    if (decomp >= 0) {
                        int scaffoldWithRSites = indigoDecomposedMoleculeScaffold(decomp);
                        if (scaffoldWithRSites >= 0) {
                            const char* mf = indigoMolfile(scaffoldWithRSites);
                            if (mf) result = QString::fromUtf8(mf);
                            else error = QString::fromUtf8(indigoGetLastError());
                        } else {
                            error = QString::fromUtf8(indigoGetLastError());
                        }
                        indigoFree(decomp);
                    } else {
                        error = QString("Decomposition failed: %1").arg(indigoGetLastError());
                    }
                    indigoFree(scaffold);
                } else {
                    error = "No common scaffold found - cannot decompose without one.";
                }
            } else {
                error = "Fewer than 2 structures could be parsed.";
            }
            for (int m : loadedMols) indigoFree(m);
            indigoFree(arr);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during R-group decomposition.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, error](IndigoService *s) { emit s->rgroupDecompositionFinished(result, error); });
    });
}

void IndigoService::decomposeToRGroupsPerMolecule(const QStringList &molfiles, const QStringList &labels) {
    if (molfiles.size() < 2) {
        emit rgroupPerMoleculeDecompositionFinished("", "Need at least 2 structures to decompose.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfiles, labels]() {
        QString resultsJson, error;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;
            for (const QString &mf : molfiles) {
                int mol = indigoLoadMoleculeFromString(mf.toUtf8().constData());
                if (mol >= 0) {
                    if (indigoAromatize(mol) < 0) qWarning() << "Indigo: decomposeToRGroupsPerMolecule aromatize failed:" << indigoGetLastError();
                    indigoArrayAdd(arr, mol);
                    loadedMols.append(mol);
                }
            }
            if (loadedMols.size() >= 2) {
                int scaffold = indigoExtractCommonScaffold(arr, "");
                if (scaffold > 0) {
                    int decomp = indigoDecomposeMolecules(scaffold, arr);
                    if (decomp >= 0) {
                        indigoSetOption("molfile-saving-mode", "3000");
                        int iter = indigoIterateDecomposedMolecules(decomp);
                        if (iter >= 0) {
                            QJsonArray resultsArray;
                            int idx = 0;
                            while (indigoHasNext(iter)) {
                                int item = indigoNext(iter);
                                if (item < 0) break;
                                int withR = indigoDecomposedMoleculeWithRGroups(item);
                                if (withR > 0) {
                                    indigoLayout(withR);
                                    const char* mf = indigoMolfile(withR);
                                    if (mf) {
                                        QJsonObject obj;
                                        obj["index"] = idx;
                                        obj["label"] = (idx < labels.size() && !labels[idx].isEmpty()) ? labels[idx] : QString("Record %1").arg(idx + 1);
                                        obj["molfile"] = QString::fromUtf8(mf);
                                        resultsArray.append(obj);
                                    } else {
                                        qWarning() << "Indigo: decomposeToRGroupsPerMolecule molfile null for item" << idx << ":" << indigoGetLastError();
                                    }
                                    indigoFree(withR);
                                } else {
                                    qWarning() << "Indigo: decomposeToRGroupsPerMolecule withRGroups failed for item" << idx << ":" << indigoGetLastError();
                                }
                                indigoFree(item);
                                idx++;
                            }
                            indigoFree(iter);
                            resultsJson = QString::fromUtf8(QJsonDocument(resultsArray).toJson(QJsonDocument::Compact));
                        } else {
                            error = QString("Failed to iterate decomposed molecules: %1").arg(indigoGetLastError());
                        }
                        indigoFree(decomp);
                    } else {
                        error = QString("Decomposition failed: %1").arg(indigoGetLastError());
                    }
                    indigoFree(scaffold);
                } else {
                    error = "No common scaffold found - cannot decompose without one.";
                }
            } else {
                error = "Fewer than 2 structures could be parsed.";
            }
            for (int m : loadedMols) indigoFree(m);
            indigoFree(arr);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during R-group decomposition.";
        }
        indigoReleaseSessionId(sid);
        emitOnGuiThread(self, [resultsJson, error](IndigoService *s) { emit s->rgroupPerMoleculeDecompositionFinished(resultsJson, error); });
    });
}

void IndigoService::rankBySimilarity(const QString &refMolfile, const QStringList &molfiles) {
    if (refMolfile.isEmpty() || molfiles.isEmpty()) {
        emit similarityRankFinished("", "Need an active structure and at least one candidate to rank.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, refMolfile, molfiles]() {
        QString error;
        QJsonArray results;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int ref = indigoLoadMoleculeFromString(refMolfile.toUtf8().constData());
            if (ref >= 0) {
                // Same aromatize-before-compare fix as similarity() - verified there that
                // comparing benzene to itself scores 0.0769 without it, 1.0000 with it.
                if (indigoAromatize(ref) < 0) qWarning() << "Indigo: rankBySimilarity aromatize(ref) failed:" << indigoGetLastError();
                for (int i = 0; i < molfiles.size(); ++i) {
                    int mol = indigoLoadMoleculeFromString(molfiles[i].toUtf8().constData());
                    if (mol >= 0) {
                        if (indigoAromatize(mol) < 0) qWarning() << "Indigo: rankBySimilarity aromatize(mol) failed:" << indigoGetLastError();
                        float sim = indigoSimilarity(ref, mol, "tanimoto");
                        QJsonObject o;
                        o["index"] = i;
                        o["score"] = sim;
                        results.append(o);
                        indigoFree(mol);
                    }
                }
                indigoFree(ref);
            } else {
                error = "Could not parse the active structure.";
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during similarity ranking.";
        }
        indigoReleaseSessionId(sid);
        QString json;
        if (error.isEmpty()) json = QString::fromUtf8(QJsonDocument(results).toJson(QJsonDocument::Compact));
        emitOnGuiThread(self, [json, error](IndigoService *s) { emit s->similarityRankFinished(json, error); });
    });
}

void IndigoService::alignBatchToScaffold(const QStringList &molfiles) {
    if (molfiles.size() < 2) {
        emit batchAlignFinished("", "Need at least 2 structures to align to a common scaffold.");
        return;
    }
    QPointer<IndigoService> self = this;
    (void)QtConcurrent::run([self, molfiles]() {
        QString error;
        QJsonArray results;
        unsigned long long sid = indigoAllocSessionId();
        indigoSetSessionId(sid);
        try {
            int arr = indigoCreateArray();
            QList<int> loadedMols;
            for (const QString &mf : molfiles) {
                int mol = indigoLoadMoleculeFromString(mf.toUtf8().constData());
                if (mol >= 0) {
                    if (indigoAromatize(mol) < 0) qWarning() << "Indigo: alignBatchToScaffold aromatize failed:" << indigoGetLastError();
                    indigoArrayAdd(arr, mol);
                    loadedMols.append(mol);
                } else {
                    loadedMols.append(-1);
                }
            }
            int validCount = 0;
            for (int m : loadedMols) {
                if (m >= 0) validCount++;
            }
            if (validCount >= 2) {
                int scaffold = indigoExtractCommonScaffold(arr, "");
                if (scaffold > 0) {
                    if (indigoLayout(scaffold) < 0) qWarning() << "Indigo: alignBatchToScaffold scaffold layout failed:" << indigoGetLastError();
                    QVector<float> scaffXYZ;
                    QVector<int> scaffAtoms;
                    int satoms = indigoIterateAtoms(scaffold);
                    int sa;
                    while ((sa = indigoNext(satoms)) > 0) {
                        float* xyz = indigoXYZ(sa);
                        if (!xyz) continue; // indigoXYZ returns null if the scaffold atom has no
                                             // coordinates (e.g. the layout above failed) -- every
                                             // other Indigo call in this function already checks
                                             // its return before use, this one didn't.
                        scaffXYZ << xyz[0] << xyz[1] << xyz[2];
                        scaffAtoms << sa;
                    }
                    if (satoms > 0) indigoFree(satoms);

                    for (int i = 0; i < molfiles.size(); ++i) {
                        int mol = loadedMols[i];
                        if (mol >= 0) {
                            int matcher = indigoSubstructureMatcher(mol, "");
                            int match = matcher >= 0 ? indigoMatch(matcher, scaffold) : -1;
                            QVector<int> atomIds;
                            QVector<float> desired;
                            if (match > 0) {
                                for (int j = 0; j < scaffAtoms.size(); ++j) {
                                    int mapped = indigoMapAtom(match, scaffAtoms[j]);
                                    if (mapped > 0) {
                                        atomIds << indigoIndex(mapped);
                                        desired << scaffXYZ[j*3] << scaffXYZ[j*3+1] << scaffXYZ[j*3+2];
                                        indigoFree(mapped);
                                    }
                                }
                                indigoFree(match);
                            }
                            if (matcher >= 0) indigoFree(matcher);

                            if (!atomIds.isEmpty()) {
                                indigoAlignAtoms(mol, atomIds.size(), atomIds.data(), desired.data());
                            }
                            const char* mf = indigoMolfile(mol);
                            results.append(mf ? QString::fromUtf8(mf) : molfiles[i]);
                        } else {
                            results.append(molfiles[i]);
                        }
                    }
                    indigoFree(scaffold);
                } else {
                    error = "No common scaffold found - cannot align without one.";
                }
            } else {
                error = "Fewer than 2 structures could be parsed.";
            }
            for (int m : loadedMols) {
                if (m >= 0) indigoFree(m);
            }
            indigoFree(arr);
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during batch alignment.";
        }
        indigoReleaseSessionId(sid);
        QString json;
        if (error.isEmpty()) json = QString::fromUtf8(QJsonDocument(results).toJson(QJsonDocument::Compact));
        emitOnGuiThread(self, [json, error](IndigoService *s) { emit s->batchAlignFinished(json, error); });
    });
}



