#include "v8_process.h"
#include <QJsonDocument>
#include <QJsonObject>
#include <QEventLoop>
#include <QTimer>
#include <QJsonArray>
#include <QDebug>
#include <QCoreApplication>
#include <QDir>
#include <QGuiApplication>
#include <QClipboard>
#include "app/molecule/RenderPrimitives.h"
#include "app/molecule/RenderPrimitivesToVariant.h"
#include "app/molecule/ClipboardPreview.h"
#include "app/molecule/EditableMolecule.h"
#include "app/PlacementEngines.h"
#include "app/BondAngleSuggester.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <QtConcurrent>
#include <QPointer>
#include "indigo.h"

V8Process::V8Process(QObject *parent, TemplateLibrary* templateLibrary)
    : QObject(parent), m_templateLibrary(templateLibrary) {
    m_docState = std::make_unique<DocumentState>();
}

V8Process::~V8Process() = default;

void V8Process::applyLocalState(const QString &resultTag) {
    RenderPrimitives prims = RenderPrimitiveBuilder::build(m_docState->molecule(), m_docState->showExplicitH());
    m_primitives = renderPrimitivesToVariant(prims);
    m_selection = selectionStateToVariant(m_docState->selection());
    emit primitivesChanged();
    emit selectionChanged();
    // resultTag lets a caller mark this update as not a structural edit (e.g.
    // "selectionOnly" for pure selection changes) -- MainWindow.qml's onStateUpdated
    // uses this to avoid invalidating a just-generated IUPAC name on a selection-only
    // update, the same mechanism already used for "stereoUpdated"/"checkUpdated".
    emit stateUpdated(m_primitives, m_selection, m_docState->isDirty(),
                       m_docState->canUndo(), m_docState->canRedo(),
                       resultTag.isEmpty() ? QVariant() : QVariant(resultTag));
}

QString V8Process::buildSdfBatchListJson() const {
    int count = m_sdfBatch.recordCount();
    int limit = std::min(count, 500);
    QJsonArray records;
    for (int i = 0; i < limit; ++i) {
        SdfBatch::Thumbnail thumb = m_sdfBatch.thumbnailAt(i);
        QJsonArray atomsArr;
        for (const SdfBatch::ThumbnailAtom& a : thumb.atoms) {
            QJsonObject o;
            o["x"] = a.x;
            o["y"] = a.y;
            o["label"] = a.label;
            atomsArr.append(o);
        }
        QJsonArray bondsArr;
        for (const SdfBatch::ThumbnailBond& b : thumb.bonds) {
            QJsonObject o;
            o["x1"] = b.x1;
            o["y1"] = b.y1;
            o["x2"] = b.x2;
            o["y2"] = b.y2;
            o["type"] = b.type;
            bondsArr.append(o);
        }
        QJsonObject thumbObj;
        thumbObj["atoms"] = atomsArr;
        thumbObj["bonds"] = bondsArr;
        QJsonObject rec;
        rec["index"] = i;
        rec["label"] = m_sdfBatch.labelAt(i);
        rec["thumb"] = thumbObj;
        records.append(rec);
    }
    QJsonObject root;
    root["count"] = count;
    root["records"] = records;
    return QString::fromUtf8(QJsonDocument(root).toJson(QJsonDocument::Compact));
}

void V8Process::handleDeserializeBatchFromMolfiles(const QString& recordsJson) {
    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(recordsJson.toUtf8(), &err);
    QStringList molfiles;
    if (err.error == QJsonParseError::NoError && doc.isArray()) {
        for (const QJsonValue& v : doc.array()) {
            molfiles.append(v.toObject().value(QStringLiteral("molfile")).toString());
        }
    }
    if (!m_sdfBatch.loadFromMolfileList(molfiles)) {
        emit structureReady(QStringLiteral("sdf_batch_list"), QStringLiteral("{\"count\":0,\"records\":[]}"));
        return;
    }
    emit structureReady(QStringLiteral("sdf_batch_list"), buildSdfBatchListJson());
}

// Mirrors 70-biopolymer.js's bioGetSequenceViewSnapshot(): the JS worker calls this
// after every bio mutation to push {monomers,bonds,seqType} back to
// BiopolymerSequenceView.qml's openSnapshot(). The C++ port's bio command handlers below
// were missing this emit entirely, leaving "Preview as Sequence View" a silent no-op.
static QString buildBioSequenceViewJson(const DocumentState& docState) {
    QJsonArray monomersArr;
    for (const BioMonomer& m : docState.bioMonomers()) {
        QJsonObject o;
        o["id"] = m.id;
        o["label"] = m.label;
        o["alias"] = m.alias;
        o["x"] = m.x;
        o["y"] = m.y;
        o["monomerClass"] = m.monomerClass;
        o["ambiguous"] = m.ambiguous;
        monomersArr.append(o);
    }
    QJsonArray bondsArr;
    for (const BioBond& b : docState.bioBonds()) {
        QJsonObject o;
        o["fromId"] = b.fromId;
        o["toId"] = b.toId;
        bondsArr.append(o);
    }
    QJsonObject root;
    root["monomers"] = monomersArr;
    root["bonds"] = bondsArr;
    root["seqType"] = docState.bioSeqType();
    return QString::fromUtf8(QJsonDocument(root).toJson(QJsonDocument::Compact));
}

void V8Process::sendCommand(const QString& cmd, const QVariantList& args) {
        if (cmd == "getMoleculeName") {
            // Read-only: emits directly, no applyLocalState() -- nothing mutated. reqId "mol_name"
            // and raw-string (not JSON) data confirmed against 40-serialize.js:104-106.
            emit structureReady(QStringLiteral("mol_name"), m_docState->molecule().name());
            return;
        }
        if (cmd == "getSaltsAndSolventsList") {
            if (!m_templateLibrary) return;
            QStringList names = m_templateLibrary->saltOrSolventNames();
            QJsonArray arr;
            for (const QString& n : names) {
                QJsonObject o;
                o["label"] = n;
                arr.append(o);
            }
            emit structureReady(QStringLiteral("salts"),
                                 QString::fromUtf8(QJsonDocument(arr).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "getFunctionalGroupsList") {
            if (!m_templateLibrary) return;
            QStringList names = m_templateLibrary->functionalGroupNames();
            std::sort(names.begin(), names.end());
            QJsonArray arr;
            for (const QString& n : names) {
                QJsonObject o;
                o["label"] = n;
                o["group"] = m_templateLibrary->functionalGroupGroup(n);
                arr.append(o);
            }
            emit structureReady(QStringLiteral("fg_list"),
                                 QString::fromUtf8(QJsonDocument(arr).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "getTemplateLibraryList") {
            if (!m_templateLibrary) return;
            QStringList names = m_templateLibrary->libraryTemplateNames();
            QJsonArray arr;
            for (const QString& n : names) {
                QJsonObject o;
                o["label"] = n;
                o["group"] = m_templateLibrary->libraryTemplateGroup(n);
                arr.append(o);
            }
            emit structureReady(QStringLiteral("library_list"),
                                 QString::fromUtf8(QJsonDocument(arr).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "getClipboardPreview") {
            if (m_docClipboardMol.isEmpty()) {
                emit structureReady(QStringLiteral("clipboard_preview"), QStringLiteral("null"));
                return;
            }
            ClipboardPreview preview = buildClipboardPreview(m_docClipboardMol);
            QJsonArray atomsArr;
            for (const ClipboardPreviewAtom& a : preview.atoms) {
                QJsonObject o;
                o["x"] = a.x;
                o["y"] = a.y;
                o["label"] = a.label;
                atomsArr.append(o);
            }
            QJsonArray bondsArr;
            for (const ClipboardPreviewBond& b : preview.bonds) {
                QJsonObject o;
                o["x1"] = b.x1;
                o["y1"] = b.y1;
                o["x2"] = b.x2;
                o["y2"] = b.y2;
                bondsArr.append(o);
            }
            QJsonObject root;
            root["atoms"] = atomsArr;
            root["bonds"] = bondsArr;
            root["cx"] = preview.cx;
            root["cy"] = preview.cy;
            emit structureReady(QStringLiteral("clipboard_preview"),
                                 QString::fromUtf8(QJsonDocument(root).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "deserializeRdfBatch" && !args.isEmpty()) {
            handleDeserializeBatchFromMolfiles(args[0].toString());
            return;
        }
        if (cmd == "deserializeIndigoBatch" && !args.isEmpty()) {
            handleDeserializeBatchFromMolfiles(args[0].toString());
            return;
        }
        if (cmd == "deserializeSdfBatch" && !args.isEmpty()) {
            if (!m_sdfBatch.loadFromSdfText(args[0].toString())) {
                emit structureReady(QStringLiteral("sdf_batch_list"), QStringLiteral("{\"count\":0,\"records\":[]}"));
                return;
            }
            emit structureReady(QStringLiteral("sdf_batch_list"), buildSdfBatchListJson());
            return;
        }
        if (cmd == "loadSdfBatchRecord" && !args.isEmpty()) {
            int index = args[0].toInt();
            QString mf = m_sdfBatch.molfileAt(index);
            if (mf.isEmpty()) return; // matches the real code's silent no-response early return
            m_docState->deserializeMol(mf);
            m_sdfProps = m_sdfBatch.propsAt(index);
            applyLocalState();
            emit structureReady(QStringLiteral("sdf_batch_load"), QString());
            return;
        }
        if (cmd == "getSdfBatchMolfiles") {
            int limit = std::min(m_sdfBatch.recordCount(), 500);
            QJsonArray molfilesArr, labelsArr;
            for (int i = 0; i < limit; ++i) {
                QString mf = m_sdfBatch.molfileAt(i);
                if (mf.isEmpty()) continue; // matches the real per-record skip-on-failure
                molfilesArr.append(mf);
                labelsArr.append(m_sdfBatch.labelAt(i));
            }
            QJsonObject root;
            root["molfiles"] = molfilesArr;
            root["labels"] = labelsArr;
            emit structureReady(QStringLiteral("sdf_batch_molfiles"),
                                 QString::fromUtf8(QJsonDocument(root).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "realignSdfBatch" && !args.isEmpty()) {
            QJsonParseError err;
            QJsonDocument doc = QJsonDocument::fromJson(args[0].toString().toUtf8(), &err);
            if (err.error != QJsonParseError::NoError || !doc.isArray()) return; // matches the
                // real code's try/catch: a JSON.parse throw on malformed input silently returns
                // with no response, regardless of the current batch's state -- must be checked
                // BEFORE calling realign(), since a malformed-input empty list could otherwise
                // coincidentally match an empty batch's length (0 == 0) and wrongly succeed
            QStringList aligned;
            for (const QJsonValue& v : doc.array()) aligned.append(v.toString());
            if (!m_sdfBatch.realign(aligned)) return; // matches the real code's silent no-response early return
            emit structureReady(QStringLiteral("sdf_batch_realigned"), buildSdfBatchListJson());
            return;
        }
        if (cmd == "getSdfProps") {
            QJsonObject obj;
            for (auto it = m_sdfProps.constBegin(); it != m_sdfProps.constEnd(); ++it) obj[it.key()] = it.value();
            emit structureReady(QStringLiteral("sdf_props"),
                                 QString::fromUtf8(QJsonDocument(obj).toJson(QJsonDocument::Compact)));
            return;
        }
        if (cmd == "deleteRxnArrow" && !args.isEmpty()) {
            m_docState->deleteRxnArrow(args[0].toInt());
        } else if (cmd == "deleteRxnPlus" && !args.isEmpty()) {
            m_docState->deleteRxnPlus(args[0].toInt());
        } else if (cmd == "moveImage" && args.size() >= 3) {
            m_docState->moveImage(args[0].toInt(), args[1].toDouble(), args[2].toDouble());
        } else if (cmd == "resizeImage" && args.size() >= 2) {
            m_docState->resizeImage(args[0].toInt(), args[1].toDouble());
        } else if (cmd == "toggleSgroupExpanded" && !args.isEmpty()) {
            m_docState->toggleSgroupExpanded(args[0].toInt());
        } else if (cmd == "renameSgroup" && args.size() >= 2) {
            m_docState->renameSGroup(args[0].toInt(), args[1].toString());
        } else if (cmd == "addBracketSelection") {
            m_docState->addBracketSelection();
        } else if (cmd == "insertRecognizedStructure" && args.size() >= 3) {
            m_docState->insertStructureAt(args[0].toString(), args[1].toDouble(), args[2].toDouble());
        } else if (cmd == "setMoleculeName" && !args.isEmpty()) {
            m_docState->setMoleculeName(args[0].toString());
        } else if (cmd == "bioBuildSequenceView" && args.size() >= 2) {
            m_docState->buildBioSequenceView(args[0].toString(), args[1].toString());
            emit structureReady(QStringLiteral("biopolymer_seq_view"), buildBioSequenceViewJson(*m_docState));
        } else if (cmd == "bioAddMonomer" && args.size() >= 2) {
            m_docState->addBioMonomer(args[0].toString(), args[1].toString());
            emit structureReady(QStringLiteral("biopolymer_seq_view"), buildBioSequenceViewJson(*m_docState));
        } else if (cmd == "bioDeleteMonomer" && !args.isEmpty()) {
            m_docState->deleteBioMonomer(args[0].toInt());
            emit structureReady(QStringLiteral("biopolymer_seq_view"), buildBioSequenceViewJson(*m_docState));
        } else if (cmd == "selectRing") {
            m_docState->selectRing(
                (!args.isEmpty() && args[0].isValid()) ? args[0].toInt() : -1,
                (args.size() > 1 && args[1].isValid()) ? args[1].toInt() : -1);
        } else if (cmd == "selectChain") {
            m_docState->selectChain(
                (!args.isEmpty() && args[0].isValid()) ? args[0].toInt() : -1,
                (args.size() > 1 && args[1].isValid()) ? args[1].toInt() : -1);
        } else if (cmd == "selectSubstructureMatches" && !args.isEmpty()) {
            m_docState->selectSubstructureMatches(args[0].toString());
        } else if (cmd == "setShowExplicitH" && !args.isEmpty()) {
            m_docState->setShowExplicitH(args[0].toBool());
        } else if (cmd == "layoutSelectedChain") {
            m_docState->layoutSelectedChain();
        } else {
            return; // unhandled command name on a C++-engine document: silent no-op,
                     // matching every other unsupported gesture's existing behavior
        }
        applyLocalState();
        return;
}

void V8Process::init() {
    applyLocalState(); return; // m_docState is always a freshly-constructed empty document (see constructor)
}
void V8Process::loadMol(const QString& molfile) {
    m_docState->deserializeMol(molfile); applyLocalState(); return;
}
void V8Process::addAtom(const QString& label, double x, double y, int charge) {
    AtomId id = m_docState->addAtom(label, x, y);
    if (charge != 0) m_docState->changeAtomCharge(id, charge);
    applyLocalState();
}
void V8Process::addBondAndAtom(int beginAtomId, const QString& endAtomLabel, double x, double y, int bondType, int stereo) {
    m_docState->addBondAndAtom(beginAtomId, endAtomLabel, x, y, bondType, stereo);
    applyLocalState();
}
void V8Process::addBondBetweenCoords(double x1, double y1, double x2, double y2, int type, int stereo) {
    m_docState->addBondBetweenCoords(x1, y1, x2, y2, type, stereo);
    applyLocalState();
}
void V8Process::addBond(int beginAtomId, int endAtomId, int bondType, int stereoDir) {
    m_docState->addBond(beginAtomId, endAtomId, bondType, stereoDir); applyLocalState(); return;
}
QVariantMap V8Process::suggestBondEndpoint(int fromAtomId, double bondLength) {
    QVariantMap atomsById = m_primitives.value(QStringLiteral("atomsById")).toMap();
    QVariantMap fromAtom = atomsById.value(QString::number(fromAtomId)).toMap();
    QVariantMap result;
    if (fromAtom.isEmpty()) return result;
    double cx = fromAtom.value(QStringLiteral("x")).toDouble();
    double cy = fromAtom.value(QStringLiteral("y")).toDouble();
    QVariantList bondsList = m_primitives.value(QStringLiteral("bonds")).toList();

    std::optional<double> angle = BondAngleSuggester::suggestAngle(fromAtomId, atomsById, bondsList);
    if (angle.has_value()) {
        result[QStringLiteral("x")] = cx + std::cos(*angle) * bondLength;
        result[QStringLiteral("y")] = cy + std::sin(*angle) * bondLength;
        return result;
    }

    // 0 or 3+ existing neighbors: unchanged fallback.
    QList<double> existingAngles;
    for (const QVariant& bVar : bondsList) {
        QVariantMap b = bVar.toMap();
        int beginId = b.value(QStringLiteral("begin")).toInt();
        int endId = b.value(QStringLiteral("end")).toInt();
        int otherId = -1;
        if (beginId == fromAtomId) otherId = endId;
        else if (endId == fromAtomId) otherId = beginId;
        if (otherId == -1) continue;
        QVariantMap other = atomsById.value(QString::number(otherId)).toMap();
        if (other.isEmpty()) continue;
        double ox = other.value(QStringLiteral("x")).toDouble();
        double oy = other.value(QStringLiteral("y")).toDouble();
        existingAngles.append(std::atan2(oy - cy, ox - cx));
    }
    QPointF start(cx, cy);
    PlacementResult placed = AtomPlacementEngine::compute(start, start, QStringLiteral("C"), existingAngles, bondLength);
    if (placed.valid && !placed.atoms.isEmpty()) {
        result[QStringLiteral("x")] = placed.atoms[0].pos.x();
        result[QStringLiteral("y")] = placed.atoms[0].pos.y();
    }
    return result;
}
// Same {x,y}-or-empty-map convention as suggestBondEndpoint, for the FG_/SS_/LIB_/TEMPLATE_
// click-to-place path (AppController::handleDragEnd). Returns an empty map for 0 or 3+
// existing neighbors -- the caller's own existing fallback (FragmentPlacementEngine::compute)
// covers those, unchanged.
QVariantMap V8Process::suggestFragmentAttachPoint(int fromAtomId, double bondLength) {
    QVariantMap atomsById = m_primitives.value(QStringLiteral("atomsById")).toMap();
    QVariantMap fromAtom = atomsById.value(QString::number(fromAtomId)).toMap();
    QVariantMap result;
    if (fromAtom.isEmpty()) return result;
    double cx = fromAtom.value(QStringLiteral("x")).toDouble();
    double cy = fromAtom.value(QStringLiteral("y")).toDouble();
    QVariantList bondsList = m_primitives.value(QStringLiteral("bonds")).toList();

    std::optional<double> angle = BondAngleSuggester::suggestAngle(fromAtomId, atomsById, bondsList);
    if (!angle.has_value()) return result;
    result[QStringLiteral("x")] = cx + std::cos(*angle) * bondLength;
    result[QStringLiteral("y")] = cy + std::sin(*angle) * bondLength;
    return result;
}
void V8Process::addRing(const QVariantList& coords, bool aromatic) {
    QList<double> pts;
    for (const QVariant& v : coords) pts.append(v.toDouble());
    m_docState->addRing(pts, aromatic);
    applyLocalState();
}
void V8Process::deleteAtomById(int id) {
    m_docState->deleteAtom(id); applyLocalState(); return;
}
void V8Process::deleteBondById(int id) {
    m_docState->deleteBond(id); applyLocalState(); return;
}
void V8Process::deleteSelection() {
    m_docState->deleteSelectionEntities();
    applyLocalState();
}
void V8Process::undo() {
    m_docState->undo();
    applyLocalState();
}
void V8Process::redo() {
    m_docState->redo();
    applyLocalState();
}
void V8Process::copySelection() {
    m_docClipboardMol = m_docState->copySelection(); return; // pure read, no applyLocalState() -- nothing mutated
}
void V8Process::cutSelection() {
    m_docClipboardMol = m_docState->cutSelection(); applyLocalState(); return;
}
void V8Process::pasteSelection(double cx, double cy) {
    if (m_docClipboardMol.isEmpty()) return; // matches the real "if (!_clipboard) return"
    m_docState->insertStructureAt(m_docClipboardMol, cx, cy);
    applyLocalState();
}
void V8Process::selectAll() {
    m_docState->selectAll(); applyLocalState("selectionOnly"); return;
}
void V8Process::clearCanvas() {
    m_docState->clearCanvas(); applyLocalState(); return;
}
void V8Process::loadBenzene() {
    m_docState->loadBenzene(); applyLocalState(); return;
}
void V8Process::deserializeMol(const QString& data) {
    m_docState->deserializeMol(data); applyLocalState(); return;
}
void V8Process::requestStructure(const QString& fmt, const QString& reqId) {
    // Async, signal-based structure export: the actual Save UI flow (MainWindow.qml's
    // saveActive()) calls THIS method -- not the synchronous getStructure() below -- and
    // listens for structureReady to actually write the file to disk. Emitting
    // structureReady synchronously here is what makes save-to-file work; this exactly
    // matches how getStructure() itself already works.
    if (fmt == "mol") {
        StringResult r = m_docState->molecule().toMolfile();
        emit structureReady(reqId, r.success ? r.value : QString());
    } else {
        emit structureReady(reqId, QString()); // only molfile export is in scope for this pilot
    }
}
void V8Process::requestSelectionStructure(const QString& reqId) {
    m_docClipboardMol = m_docState->copySelection();
    if (m_docClipboardMol.isEmpty()) {
        emit structureReady(reqId, QString());
        return;
    }
    EditableMolecule mol(m_docClipboardMol);
    emit structureReady(reqId, mol.isValid() ? mol.toKetJson() : QString());
}
void V8Process::requestSerialize(const QString& reqId) {
    requestStructure("mol", reqId);
}
void V8Process::requestMolfileForNaming(const QString& reqId) {
    QString mf = getMolfileForNaming();
    emit structureReady(reqId, mf);
}
void V8Process::requestMolfileForNamingSelection(const QString& reqId) {
    emit structureReady(reqId, m_docState->copySelection());
}
QString V8Process::getStructure(const QString& fmt) {
    if (fmt != "mol") return QString(); // only molfile export is in scope for this pilot
    StringResult r = m_docState->molecule().toMolfile();
    return r.success ? r.value : QString();
}
QString V8Process::serializeMol() {
    return getStructure("mol");
}
QString V8Process::getMolfileForNaming() {
    return m_docState->molecule().molfileForNaming();
}
bool V8Process::importReaction(const QString& text) {
    if (!m_docState->importReaction(text)) return false;
    applyLocalState();
    return true;
}
void V8Process::loadStructure(const QString& format, const QString& data, bool centerOnPage) {
    if (format == "mol") {
        m_docState->deserializeMol(data, centerOnPage);
        applyLocalState();
        return;
    }
    if (format == "ket") {
        EditableMolecule mol(data);
        if (!mol.isValid()) return;
        StringResult mf = mol.toMolfile();
        if (!mf.success) return;
        m_docState->deserializeMol(mf.value, centerOnPage);
        applyLocalState();
        return;
    }
    // sdf: out of scope -- deserializeSdf is confirmed dead (zero QML callers;
    // MainWindow.qml's own File > Open handler routes .sdf through
    // deserializeSdfBatch directly, never through here)
}
void V8Process::insertFunctionalGroup(const QString& fgName, double cx, double cy, int targetAtomId, bool fullStructure) {
    if (!m_templateLibrary) return;
    m_docState->insertFunctionalGroup(*m_templateLibrary, fgName, cx, cy, targetAtomId, fullStructure);
    applyLocalState();
}
void V8Process::insertLibraryTemplateFused(const QString& fgName, double cx, double cy, int targetBondId) {
    if (!m_templateLibrary) return;
    m_docState->insertLibraryTemplateFused(*m_templateLibrary, fgName, cx, cy, targetBondId);
    applyLocalState();
}
void V8Process::requestSaltsAndSolventsList() { sendCommand("getSaltsAndSolventsList"); }
void V8Process::requestFunctionalGroupsList() { sendCommand("getFunctionalGroupsList"); }
void V8Process::requestTemplateLibraryList() { sendCommand("getTemplateLibraryList"); }
// Render + normalize a single thumbnail's atoms/bonds. Runs entirely on whatever thread calls
// it, using whatever Indigo session is already active there -- caller's responsibility.
static QString renderThumbnailJson(const QString& molfile) {
    EditableMolecule tmp(molfile);
    RenderPrimitives prims = RenderPrimitiveBuilder::build(tmp, /*showExplicitH=*/false);

    QJsonArray atomsArr;
    QJsonArray bondsArr;
    if (prims.bbox.valid) {
        double w = prims.bbox.maxX - prims.bbox.minX;
        double h = prims.bbox.maxY - prims.bbox.minY;
        if (w <= 0) w = 1;
        if (h <= 0) h = 1;
        const double pad = 0.1;
        double scale = (1.0 - 2.0 * pad) / std::max(w, h);
        double offX = pad + (1.0 - 2.0 * pad - w * scale) / 2.0;
        double offY = pad + (1.0 - 2.0 * pad - h * scale) / 2.0;

        QHash<AtomId, QPointF> positions;
        for (const AtomPrim& a : prims.atoms) {
            double nx = offX + (a.x - prims.bbox.minX) * scale;
            double ny = offY + (a.y - prims.bbox.minY) * scale;
            positions.insert(a.id, QPointF(nx, ny));
            QJsonObject o;
            o["x"] = nx;
            o["y"] = ny;
            o["label"] = a.label;
            o["color"] = a.color;
            atomsArr.append(o);
        }
        for (const BondPrim& b : prims.bonds) {
            if (!positions.contains(b.begin) || !positions.contains(b.end)) continue;
            QPointF p1 = positions.value(b.begin);
            QPointF p2 = positions.value(b.end);
            QJsonObject o;
            o["x1"] = p1.x();
            o["y1"] = p1.y();
            o["x2"] = p2.x();
            o["y2"] = p2.y();
            o["type"] = b.type > 0 ? b.type : 1;
            o["stereo"] = b.stereo;
            bondsArr.append(o);
        }
    }

    QJsonObject result;
    result["atoms"] = atomsArr;
    result["bonds"] = bondsArr;
    return QString::fromUtf8(QJsonDocument(result).toJson(QJsonDocument::Compact));
}

void V8Process::requestTemplateThumbnail(const QString& name, const QString& reqId) {
    if (m_templateLibrary) {
        int handle = m_templateLibrary->functionalGroup(name);
        if (handle < 0) handle = m_templateLibrary->saltOrSolvent(name);
        if (handle < 0) handle = m_templateLibrary->libraryTemplate(name);
        if (handle < 0) {
            emit structureReady(reqId, QStringLiteral("{}"));
            return;
        }

        // Rendering happens on a background thread, not here: a handful of real library
        // templates (confirmed empirically -- e.g. "C20H20", a dodecahedrane cage) make
        // Indigo's own SSSR ring-perception (EditableMolecule::ringMembership, called from
        // RenderPrimitiveBuilder::build) take an unbounded amount of time for certain highly
        // symmetric polycyclic topologies -- a known hard case for ring-perception algorithms
        // in general, not a bug in this port. On the main thread this manifested as a real,
        // reproduced Windows AppHang (confirmed via Event Viewer) the moment the Library
        // picker opened and synchronously requested all 249 thumbnails. Backgrounding it, the
        // exact same pattern already proven throughout IndigoService.cpp (fresh session
        // per task, QPointer-guarded emit marshaled back via QMetaObject::invokeMethod on
        // qApp), means a single pathological template can only ever delay its OWN thumbnail
        // -- never block the UI thread or the rest of the app.
        QString molfile = m_templateLibrary->molfileText(handle);
        QPointer<V8Process> self = this;

        (void)QtConcurrent::run([self, reqId, molfile]() {
            unsigned long long threadSessionId = indigoAllocSessionId();
            indigoSetSessionId(threadSessionId);

            QString json = QStringLiteral("{\"atoms\":[],\"bonds\":[]}");
            try {
                json = renderThumbnailJson(molfile);
            } catch (...) {
            }

            indigoReleaseSessionId(threadSessionId);

            QMetaObject::invokeMethod(qApp, [self, reqId, json]() {
                if (self) emit self->structureReady(reqId, json);
            }, Qt::QueuedConnection);
        });
        return;
    }
}
void V8Process::addRxnArrow(double x, double y, const QString& mode) {
    m_docState->addRxnArrow(x, y, mode); applyLocalState(); return;
}
void V8Process::addRxnPlus(double x, double y) {
    m_docState->addRxnPlus(x, y); applyLocalState(); return;
}
void V8Process::addCurvedArrow(double x1, double y1, double ctrlX, double ctrlY, double x2, double y2) {
    m_docState->addCurvedArrow(x1, y1, ctrlX, ctrlY, x2, y2); applyLocalState(); return;
}
void V8Process::setRxnArrowMode(int id, const QString& mode) {
    m_docState->setRxnArrowMode(id, mode); applyLocalState(); return;
}
void V8Process::setRxnArrowConditions(int id, const QString& above, const QString& below) {
    m_docState->setRxnArrowConditions(id, above, below); applyLocalState(); return;
}
void V8Process::setStereoFlags(const QString& type, int groupId) {
    m_docState->setStereoFlags(type, groupId); applyLocalState(); return;
}
void V8Process::transformSelection(const QString& mode) {
    // Mode strings confirmed against 20-edit.js:982-990.
    if (mode == "rotate_cw") m_docState->rotateSelection90CW();
    else if (mode == "rotate_ccw") m_docState->rotateSelection90CCW();
    else if (mode == "flip_h") m_docState->flipSelectionHorizontal();
    else if (mode == "flip_v") m_docState->flipSelectionVertical();
    else { emit errorOccurred("Unknown transform mode: " + mode); return; }
    applyLocalState();
}
void V8Process::addChain(double x1, double y1, double x2, double y2) {
    m_docState->addChain(x1, y1, x2, y2); applyLocalState(); return;
}
void V8Process::addText(const QString& content, double x, double y, bool bold, bool italic) {
    m_docState->addText(content, x, y, bold, italic); applyLocalState(); return;
}
void V8Process::updateText(int id, const QString& content, bool bold, bool italic) {
    m_docState->updateText(id, content, bold, italic); applyLocalState(); return;
}
void V8Process::deleteText(int id) {
    m_docState->deleteText(id); applyLocalState(); return;
}
void V8Process::addImage(const QString& base64DataUri, double cx, double cy, double halfW, double halfH) {
    // Not a real base64 decode -- ImagePrim::bitmap is round-tripped via
    // QString::fromUtf8 (RenderPrimitivesToVariant.cpp:122), meaning this port stores the
    // raw data-URI TEXT as the "image data" everywhere, same as the JS side (50-reactions.js
    // stores base64DataUri verbatim too). toUtf8() is the correct, exact-match conversion.
    m_docState->addImage(base64DataUri.toUtf8(), cx, cy, halfW, halfH);
    applyLocalState();
}
void V8Process::deleteImage(int id) {
    m_docState->deleteImage(id); applyLocalState(); return;
}
void V8Process::addRGroup(int rgroupNumber) {
    m_docState->addRGroup(rgroupNumber); applyLocalState(); return;
}
void V8Process::deleteRGroup(int rgroupNumber) {
    m_docState->deleteRGroup(rgroupNumber); applyLocalState(); return;
}
void V8Process::setRGroupLogic(int rgroupNumber, const QString& range, bool resth, int ifthen) {
    m_docState->setRGroupLogic(rgroupNumber, range, resth, ifthen); applyLocalState(); return;
}
void V8Process::addRGroupMember(int rgroupNumber) {
    m_docState->addRGroupMember(rgroupNumber); applyLocalState(); return;
}
void V8Process::removeRGroupMember(int rgroupNumber, int fragId) {
    m_docState->removeRGroupMember(rgroupNumber, fragId); applyLocalState(); return;
}
void V8Process::setAtomQueryList(int atomId, const QString& elementsCsv, bool notList) {
    m_docState->setAtomQueryList(atomId, elementsCsv, notList); applyLocalState(); return;
}
void V8Process::clearAtomQueryList(int atomId, const QString& fallbackLabel) {
    m_docState->clearAtomQueryList(atomId, fallbackLabel); applyLocalState(); return;
}
void V8Process::addMultitailArrow(double x, double y) {
    m_docState->addMultitailArrow(x, y); applyLocalState(); return;
}
void V8Process::deleteMultitailArrow(int id) {
    m_docState->deleteMultitailArrow(id); applyLocalState(); return;
}
void V8Process::addMultitailArrowTail(int id) {
    m_docState->addMultitailArrowTail(id); applyLocalState(); return;
}
void V8Process::changeAtomLabel(int id, const QString& label) {
    m_docState->changeAtomLabel(id, label); applyLocalState(); return;
}
void V8Process::setAtomMapping(int id, int mapping) {
    m_docState->setAtomMapping(id, mapping); applyLocalState(); return;
}
void V8Process::changeBondType(int id, int type, int stereo) {
    m_docState->changeBondTypeAndStereo(id, type, stereo); applyLocalState(); return;
}
void V8Process::changeAtomCharge(int id, int charge) {
    m_docState->changeAtomCharge(id, charge); applyLocalState(); return;
}
void V8Process::setAttachmentPoint(int id, int order) {
    m_docState->setAttachmentPoint(id, order); applyLocalState(); return;
}
void V8Process::changeAtomIsotope(int id, int isotope) {
    m_docState->changeAtomIsotope(id, isotope); applyLocalState(); return;
}
void V8Process::changeAtomRadical(int id, int radical) {
    m_docState->changeAtomRadical(id, radical); applyLocalState(); return;
}
void V8Process::changeAtomValence(int id, int valence) {
    m_docState->changeAtomValence(id, valence); applyLocalState(); return;
}
void V8Process::requestAtomProperties(int id) {
    // Read-only: emits directly via the same structureReady/reqId convention
    // requestStructure already established (sub-project 7b), no applyLocalState() call.
    if (!m_docState->molecule().atomIds().contains(id)) {
        emit structureReady(QStringLiteral("atom_props"), QStringLiteral("{}")); // matches real "JSON.stringify(props || {})"
        return;
    }
    DocumentState::AtomProperties props = m_docState->atomProperties(id);
    QJsonObject obj;
    obj["id"] = id;
    obj["label"] = props.label;
    obj["charge"] = props.charge;
    obj["isotope"] = props.isotope;
    obj["radical"] = props.radical;
    obj["explicitValence"] = props.explicitValence;
    emit structureReady(QStringLiteral("atom_props"), QString::fromUtf8(QJsonDocument(obj).toJson(QJsonDocument::Compact)));
}
void V8Process::selectByRect(double x1, double y1, double x2, double y2) {
    m_docState->selectByRect(x1, y1, x2, y2); applyLocalState("selectionOnly"); return;
}
void V8Process::addSelectionByRect(double x1, double y1, double x2, double y2) {
    m_docState->addSelectionByRect(x1, y1, x2, y2); applyLocalState("selectionOnly"); return;
}
void V8Process::selectByLasso(const QVariantList& pointsFlat) {
    QList<QPointF> pts;
    for (int i = 0; i + 1 < pointsFlat.size(); i += 2) {
        pts.append(QPointF(pointsFlat[i].toDouble(), pointsFlat[i + 1].toDouble()));
    }
    m_docState->selectByLasso(pts);
    applyLocalState("selectionOnly");
}
void V8Process::selectItem(const QVariant& atomId, const QVariant& bondId, const QVariant& rxnArrowId, const QVariant& rxnPlusId, const QVariant& multitailArrowId) {
    m_docState->selectSingleItem(
        atomId.isValid() ? atomId.toInt() : -1,
        bondId.isValid() ? bondId.toInt() : -1,
        rxnArrowId.isValid() ? rxnArrowId.toInt() : -1,
        rxnPlusId.isValid() ? rxnPlusId.toInt() : -1,
        multitailArrowId.isValid() ? multitailArrowId.toInt() : -1);
    applyLocalState("selectionOnly");
}
void V8Process::addItemToSelection(const QVariant& atomId, const QVariant& bondId) {
    // Atom-priority-over-bond dispatch, confirmed against 10-state.js:406-414.
    if (atomId.isValid()) m_docState->addAtomToSelection(atomId.toInt());
    else if (bondId.isValid()) m_docState->addBondToSelection(bondId.toInt());
    applyLocalState("selectionOnly");
}
void V8Process::removeItemFromSelection(const QVariant& atomId, const QVariant& bondId) {
    // Confirmed against 10-state.js:416-422.
    if (atomId.isValid()) m_docState->removeAtomFromSelection(atomId.toInt());
    else if (bondId.isValid()) m_docState->removeBondFromSelection(bondId.toInt());
    applyLocalState("selectionOnly");
}
void V8Process::selectFragment(const QVariant& atomId, const QVariant& bondId) {
    m_docState->selectFragment(
        atomId.isValid() ? atomId.toInt() : -1,
        bondId.isValid() ? bondId.toInt() : -1);
    applyLocalState("selectionOnly");
}
void V8Process::moveSelection(double dx, double dy) {
    m_docState->moveSelectionLive(dx, dy);
    applyLocalState();
}
void V8Process::commitMove() {
    m_docState->commitMove();
    applyLocalState();
}
void V8Process::rotateSelectionLive(double angleDelta) {
    m_docState->rotateSelectionLive(angleDelta); applyLocalState();
}
void V8Process::commitRotate() {
    m_docState->commitRotate(); applyLocalState();
}
void V8Process::scaleSelectionLive(double factor, double anchorX, double anchorY) {
    m_docState->scaleSelectionLive(factor, anchorX, anchorY); applyLocalState();
}
void V8Process::commitScale() {
    m_docState->commitScale(); applyLocalState();
}
void V8Process::centerStructure() { sendCommand("centerStructure"); }
void V8Process::normalizeStructure() { sendCommand("normalizeStructure"); }

void V8Process::alignAtoms(const QString& direction) {
    m_docState->alignAtoms(direction); applyLocalState();
}
void V8Process::distributeAtoms(const QString& direction) {
    m_docState->distributeAtoms(direction); applyLocalState();
}

void V8Process::setStereoDescriptors(const QString& jsonMap) {
    // Tagged "stereoUpdated" (matching the original JS engine's 60-analysis.js
    // convention, see archive/js-engine/worker/60-analysis.js) so MainWindow.qml's
    // onStateUpdated doesn't treat this routine, automatic recalculation (fired by
    // the calc_props debounce cycle after every edit) as a real structural edit --
    // without this tag it was wiping a just-generated IUPAC name every ~600ms.
    m_docState->setStereoDescriptors(jsonMap); applyLocalState("stereoUpdated");
}
void V8Process::setCheckIssues(const QString& jsonMap) {
    // Tagged "checkUpdated" for the same reason as setStereoDescriptors above.
    m_docState->setCheckIssues(jsonMap); applyLocalState("checkUpdated");
}

QString V8Process::getOsClipboardText() const {
    return QGuiApplication::clipboard()->text();
}

void V8Process::setOsClipboardText(const QString& text) {
    QGuiApplication::clipboard()->setText(text);
}

void V8Process::requestClipboardKet() {
    if (m_docClipboardMol.isEmpty()) {
        emit structureReady(QStringLiteral("clipboard_ket"), QString());
        return;
    }
    EditableMolecule mol(m_docClipboardMol);
    emit structureReady(QStringLiteral("clipboard_ket"), mol.isValid() ? mol.toKetJson() : QString());
}

void V8Process::importKetAtPosition(const QString& ket, double cx, double cy) {
    EditableMolecule mol(ket);
    if (!mol.isValid()) return;
    StringResult mf = mol.toMolfile();
    if (!mf.success) return;
    m_docState->insertStructureAt(mf.value, cx, cy);
    applyLocalState();
}

void V8Process::setOverlayState(const QVariantMap& state) {
    if (m_overlayState != state) {
        m_overlayState = state;
        emit overlayStateChanged();
    }
}

static double getAngle(double x1, double y1, double x2, double y2) {
    return std::atan2(y2 - y1, x2 - x1);
}

double V8Process::getLargestEmptyAngle(const QVariant& atomId, const QVariantMap& primitives) {
    QVariantMap atomsById = primitives.value("atomsById").toMap();
    if (!atomsById.contains(atomId.toString())) return 0;
    
    QVariantMap atom = atomsById.value(atomId.toString()).toMap();
    double ax = atom.value("x").toDouble();
    double ay = atom.value("y").toDouble();
    
    QList<double> neighbors;
    QVariantList bonds = primitives.value("bonds").toList();
    for (const QVariant& bVar : bonds) {
        QVariantMap b = bVar.toMap();
        QVariant beginId = b.value("begin");
        QVariant endId = b.value("end");
        if (beginId == atomId) {
            QVariantMap n1 = atomsById.value(endId.toString()).toMap();
            if (!n1.isEmpty()) neighbors.append(getAngle(ax, ay, n1.value("x").toDouble(), n1.value("y").toDouble()));
        } else if (endId == atomId) {
            QVariantMap n2 = atomsById.value(beginId.toString()).toMap();
            if (!n2.isEmpty()) neighbors.append(getAngle(ax, ay, n2.value("x").toDouble(), n2.value("y").toDouble()));
        }
    }
    
    if (neighbors.isEmpty()) return 0;
    if (neighbors.size() == 1) return neighbors.first() + 2.61799;
    
    std::sort(neighbors.begin(), neighbors.end());
    double maxGap = 0;
    double bestAngle = 0;
    for (int i = 0; i < neighbors.size(); i++) {
        double a1 = neighbors[i];
        double a2 = neighbors[(i + 1) % neighbors.size()];
        double gap = a2 - a1;
        while (gap <= 0) gap += M_PI * 2;
        if (gap > maxGap) {
            maxGap = gap;
            bestAngle = a1 + gap / 2;
        }
    }
    return bestAngle;
}

QVariantList V8Process::getRingPreviewCoords(int n, double cx, double cy, const QVariant& hoverAtomId, const QVariant& hoverBondId) {
    // Page/canvas boundary hard clamp, mirroring src/v8_worker.js's PAGE_MIN_X/
    // MAX_X/MIN_Y/MAX_Y (kept in sync manually -- must match). Only affects the
    // free-floating ring case below: when fusing onto an existing hover atom/bond
    // (checked further down), the ring's geometry is derived from that atom's own
    // already-in-bounds position instead, not from cx/cy, so clamping here can't
    // misalign a fused ring's seam.
    cx = std::max(-30.0, std::min(30.0, cx));
    cy = std::max(-21.0, std::min(21.0, cy));
    QVariantList coords;
    // Must match chem-core.js's StandardBondLength (MonomerSize * 2 = 0.75 * 2 = 1.5),
    // not an arbitrary 1.0 -- this function's output is used directly as the final
    // ring geometry passed to addRing() (see ChemCanvas.qml's call sites), not just a
    // cosmetic hover preview, so every ring placed via a TEMPLATE_ tool was rendering
    // at 1.0/1.5 (~67%) of the standard bond length used everywhere else in the app.
    double L = 1.5;
    double R = L / (2 * std::sin(M_PI / n));
    double angleStep = (2 * M_PI) / n;
    
    if (!hoverBondId.isNull() && hoverBondId.isValid()) {
        QVariantMap bond;
        QVariantList bonds = m_primitives.value("bonds").toList();
        for (const QVariant& bVar : bonds) {
            if (bVar.toMap().value("id") == hoverBondId) {
                bond = bVar.toMap();
                break;
            }
        }
        if (!bond.isEmpty()) {
            QVariantMap atomsById = m_primitives.value("atomsById").toMap();
            QVariantMap a1 = atomsById.value(bond.value("begin").toString()).toMap();
            QVariantMap a2 = atomsById.value(bond.value("end").toString()).toMap();
            if (!a1.isEmpty() && !a2.isEmpty()) {
                double a1x = a1.value("x").toDouble();
                double a1y = a1.value("y").toDouble();
                double a2x = a2.value("x").toDouble();
                double a2y = a2.value("y").toDouble();
                double midX = (a1x + a2x) / 2;
                double midY = (a1y + a2y) / 2;
                double vx = a2x - a1x;
                double vy = a2y - a1y;
                double len = std::sqrt(vx*vx + vy*vy);
                if (len != 0) {
                    double ux = vx / len;
                    double uy = vy / len;
                    double nx1 = -uy, ny1 = ux;
                    double nx2 = uy, ny2 = -ux;
                    double score1 = 0, score2 = 0;
                    for (auto it = atomsById.constBegin(); it != atomsById.constEnd(); ++it) {
                        if (it.key() == bond.value("begin").toString() || it.key() == bond.value("end").toString()) continue;
                        QVariantMap a = it.value().toMap();
                        double dx = a.value("x").toDouble() - midX;
                        double dy = a.value("y").toDouble() - midY;
                        double dot1 = dx * nx1 + dy * ny1;
                        double dot2 = dx * nx2 + dy * ny2;
                        if (dot1 > 0) score1 += dot1;
                        if (dot2 > 0) score2 += dot2;
                    }
                    double nx = (score1 < score2) ? nx1 : nx2;
                    double ny = (score1 < score2) ? ny1 : ny2;
                    double h = L / (2 * std::tan(M_PI / n));
                    double center_x = midX + nx * h;
                    double center_y = midY + ny * h;
                    double startAngle = getAngle(center_x, center_y, a1x, a1y);
                    double angle2 = getAngle(center_x, center_y, a2x, a2y);
                    double angleDiff = angle2 - startAngle;
                    while (angleDiff <= -M_PI) angleDiff += 2 * M_PI;
                    while (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                    double step = (angleDiff > 0) ? angleStep : -angleStep;
                    
                    QVariantMap p1; p1["x"] = a2x; p1["y"] = a2y; coords.append(p1);
                    QVariantMap p2; p2["x"] = a1x; p2["y"] = a1y; coords.append(p2);
                    for (int i = 2; i < n; i++) {
                        double a = startAngle - (i - 1) * step;
                        QVariantMap p;
                        p["x"] = center_x + R * std::cos(a);
                        p["y"] = center_y + R * std::sin(a);
                        coords.append(p);
                    }
                    return coords;
                }
            }
        }
    }
    
    if (!hoverAtomId.isNull() && hoverAtomId.isValid()) {
        QVariantMap atomsById = m_primitives.value("atomsById").toMap();
        QVariantMap atom = atomsById.value(hoverAtomId.toString()).toMap();
        if (!atom.isEmpty()) {
            double ax = atom.value("x").toDouble();
            double ay = atom.value("y").toDouble();
            double emptyAngle = getLargestEmptyAngle(hoverAtomId, m_primitives);
            double center_x = ax + R * std::cos(emptyAngle);
            double center_y = ay + R * std::sin(emptyAngle);
            double startAngle = getAngle(center_x, center_y, ax, ay);
            for (int i = 0; i < n; i++) {
                double a = startAngle + i * angleStep;
                QVariantMap p;
                p["x"] = center_x + R * std::cos(a);
                p["y"] = center_y + R * std::sin(a);
                coords.append(p);
            }
            return coords;
        }
    }
    
    double startAngle = -M_PI / 2;
    if (n % 2 == 0) startAngle += angleStep / 2;
    for (int i = 0; i < n; i++) {
        double a = startAngle + i * angleStep;
        QVariantMap p;
        p["x"] = cx + R * std::cos(a);
        p["y"] = cy + R * std::sin(a);
        coords.append(p);
    }
    return coords;
}

PlacementResult V8Process::getFunctionalGroupPlacementResult(const QString& toolId, double cx, double cy, int targetAtomId) {
    QString fgName;
    if (toolId.startsWith("LIB_")) fgName = toolId.mid(4);
    else if (toolId.startsWith("FG_") || toolId.startsWith("SS_")) fgName = toolId.mid(3);
    else return PlacementResult{};
    if (!m_templateLibrary) return PlacementResult{};
    return m_docState->getFunctionalGroupPreview(*m_templateLibrary, fgName, cx, cy, targetAtomId);
}

QVariantList V8Process::getFunctionalGroupPreview(const QString& toolId, double cx, double cy, const QVariant& hoverAtomId) {
    int targetAtomId = -1;
    if (!hoverAtomId.isNull() && hoverAtomId.isValid()) targetAtomId = hoverAtomId.toInt();
    PlacementResult result = getFunctionalGroupPlacementResult(toolId, cx, cy, targetAtomId);

    QVariantList coords;
    for (const auto& a : result.atoms) {
        QVariantMap amap;
        amap["x"] = a.pos.x();
        amap["y"] = a.pos.y();
        amap["label"] = a.label;
        coords.append(amap);
    }
    for (const auto& b : result.bonds) {
        QVariantMap bmap;
        bmap["startX"] = b.start.x();
        bmap["startY"] = b.start.y();
        bmap["endX"] = b.end.x();
        bmap["endY"] = b.end.y();
        bmap["isBond"] = true;   // distinguishes bond entries from atom entries in the flat list
        coords.append(bmap);
    }
    return coords;
}

// Project file serialization methods
void V8Process::requestProjectSerialization(const QString& reqId) {
    QJsonObject json = ProjectFile::toJson(m_docState->molecule(), 
                                          m_docState->drawingStyle(),
                                          m_docState->pageSettings());
    QJsonDocument doc(json);
    emit structureReady(reqId, doc.toJson(QJsonDocument::Compact));
}

void V8Process::loadProject(const QUrl& fileUrl) {
    // fileUrl arrives as a real file:// URL from QML's FileDialog (selectedFile)
    // or as a plain local path wrapped in a QUrl -- isLocalFile() distinguishes
    // the two so either form resolves to a real filesystem path, not the
    // literal "file:///C:/..." string QFile cannot open.
    const QString filePath = fileUrl.isLocalFile() ? fileUrl.toLocalFile() : fileUrl.toString();
    DrawingStyle style = m_docState->drawingStyle();
    PageSettings page = m_docState->pageSettings();
    QString error = ProjectFile::loadProject(filePath, m_docState->molecule(), style, page);
    if (!error.isEmpty()) {
        emit errorOccurred(QStringLiteral("Failed to load project: ") + error);
        return;
    }
    m_docState->setDrawingStyle(style);
    m_docState->setPageSettings(page);
    m_docState->markClean();
    applyLocalState();
}

bool V8Process::saveProject(const QUrl& fileUrl) {
    const QString filePath = fileUrl.isLocalFile() ? fileUrl.toLocalFile() : fileUrl.toString();
    return ProjectFile::saveProject(filePath, m_docState->molecule(),
                                   m_docState->drawingStyle(),
                                   m_docState->pageSettings());
}

QJsonObject V8Process::getDrawingStyle() const {
    const DrawingStyle& style = m_docState->drawingStyle();
    QJsonObject obj;
    obj["bondLength"] = style.bondLength;
    obj["lineWidth"] = style.lineWidth;
    obj["fontFamily"] = style.fontFamily;
    obj["fontSize"] = style.fontSize;
    obj["monochrome"] = style.monochrome;
    obj["modelToMmScale"] = style.modelToMmScale;
    return obj;
}

void V8Process::setDrawingStyle(const QJsonObject& style) {
    DrawingStyle ds;
    if (style.contains("bondLength")) ds.bondLength = style["bondLength"].toDouble();
    if (style.contains("lineWidth")) ds.lineWidth = style["lineWidth"].toDouble();
    if (style.contains("fontFamily")) ds.fontFamily = style["fontFamily"].toString();
    if (style.contains("fontSize")) ds.fontSize = style["fontSize"].toInt();
    if (style.contains("monochrome")) ds.monochrome = style["monochrome"].toBool();
    m_docState->setDrawingStyle(ds);
    applyLocalState();
}

void V8Process::applyDrawingStyleWithRescale(const QJsonObject& style) {
    DrawingStyle ds;
    if (style.contains("bondLength")) ds.bondLength = style["bondLength"].toDouble();
    else ds.bondLength = m_docState->drawingStyle().bondLength;
    if (style.contains("lineWidth")) ds.lineWidth = style["lineWidth"].toDouble();
    else ds.lineWidth = m_docState->drawingStyle().lineWidth;
    if (style.contains("fontFamily")) ds.fontFamily = style["fontFamily"].toString();
    else ds.fontFamily = m_docState->drawingStyle().fontFamily;
    if (style.contains("fontSize")) ds.fontSize = style["fontSize"].toInt();
    else ds.fontSize = m_docState->drawingStyle().fontSize;
    if (style.contains("monochrome")) ds.monochrome = style["monochrome"].toBool();
    else ds.monochrome = m_docState->drawingStyle().monochrome;
    
    // Copy over the current modelToMmScale to preserve physical scale
    ds.modelToMmScale = m_docState->drawingStyle().modelToMmScale;
    
    m_docState->applyDrawingStyleWithRescale(ds);
    applyLocalState();
}

QJsonObject V8Process::getPageSettings() const {
    const PageSettings& settings = m_docState->pageSettings();
    QJsonObject obj;
    obj["widthMm"] = settings.widthMm;
    obj["heightMm"] = settings.heightMm;
    obj["marginLeftMm"] = settings.marginLeftMm;
    obj["marginRightMm"] = settings.marginRightMm;
    obj["marginTopMm"] = settings.marginTopMm;
    obj["marginBottomMm"] = settings.marginBottomMm;
    return obj;
}

void V8Process::setPageSettings(const QJsonObject& settings) {
    PageSettings ps;
    if (settings.contains("widthMm")) ps.widthMm = settings["widthMm"].toDouble();
    if (settings.contains("heightMm")) ps.heightMm = settings["heightMm"].toDouble();
    if (settings.contains("marginLeftMm")) ps.marginLeftMm = settings["marginLeftMm"].toDouble();
    if (settings.contains("marginRightMm")) ps.marginRightMm = settings["marginRightMm"].toDouble();
    if (settings.contains("marginTopMm")) ps.marginTopMm = settings["marginTopMm"].toDouble();
    if (settings.contains("marginBottomMm")) ps.marginBottomMm = settings["marginBottomMm"].toDouble();
    m_docState->setPageSettings(ps);
    applyLocalState();
}

