#pragma once
#include <QObject>
#include <QVariantMap>
#include <QVariantList>
#include <QList>
#include <QByteArray>
#include <QUrl>
#include <QtQml/qqml.h>
#include <memory>
#include "app/molecule/DocumentState.h"
#include "app/molecule/TemplateLibrary.h"
#include "app/molecule/SdfBatch.h"
#include "app/PlacementEngines.h"
#include "app/molecule/ProjectFile.h"

class V8Process : public QObject {
    Q_OBJECT
    QML_NAMED_ELEMENT(Sketch)

    Q_PROPERTY(QVariantMap primitives READ primitives NOTIFY primitivesChanged)
    Q_PROPERTY(QVariantMap selection READ selection NOTIFY selectionChanged)
    Q_PROPERTY(QVariantMap overlayState READ overlayState WRITE setOverlayState NOTIFY overlayStateChanged)

public:
    // Every document runs on DocumentState; there is no other engine.
    // templateLibrary is non-owning: DocumentManager owns the real instance and outlives
    // every V8Process it creates. Defaulted to nullptr so any hypothetical future
    // construction site that doesn't need template insertion (e.g. a test) still
    // compiles unchanged.
    explicit V8Process(QObject *parent = nullptr, TemplateLibrary* templateLibrary = nullptr);
    ~V8Process();

    QVariantMap primitives() const { return m_primitives; }
    QVariantMap selection() const { return m_selection; }
    QVariantMap overlayState() const { return m_overlayState; }
    void setOverlayState(const QVariantMap& state);

    Q_INVOKABLE void init();
    Q_INVOKABLE void loadMol(const QString& molfile);
    Q_INVOKABLE void addAtom(const QString& label, double x, double y, int charge);
    Q_INVOKABLE void addBondAndAtom(int beginAtomId, const QString& endAtomLabel, double x, double y, int bondType, int stereo);
    Q_INVOKABLE void addBondBetweenCoords(double x1, double y1, double x2, double y2, int type, int stereo);
    Q_INVOKABLE void addBond(int beginAtomId, int endAtomId, int bondType, int stereoDir = 0);
    // Read-only: for a plain click (no drag) with a bond tool starting on an existing atom,
    // suggests where the new atom should land -- reusing AtomPlacementEngine's own
    // existing-neighbor-aware angle snap (the same engine the ATOM_/FG_ tools' click-to-extend
    // gesture already uses), so the new bond doesn't land directly on top of one that's
    // already there. Returns {} if fromAtomId doesn't exist.
    Q_INVOKABLE QVariantMap suggestBondEndpoint(int fromAtomId, double bondLength);
    Q_INVOKABLE QVariantMap suggestFragmentAttachPoint(int fromAtomId, double bondLength);
    Q_INVOKABLE void addRing(const QVariantList& coords, bool aromatic = true);
    Q_INVOKABLE QVariantList getRingPreviewCoords(int n, double cx, double cy, const QVariant& hoverAtomId, const QVariant& hoverBondId);
    PlacementResult getFunctionalGroupPlacementResult(const QString& toolId, double cx, double cy, int targetAtomId);
    Q_INVOKABLE QVariantList getFunctionalGroupPreview(const QString& toolId, double cx, double cy, const QVariant& hoverAtomId);
    Q_INVOKABLE void deleteAtomById(int id);
    Q_INVOKABLE void deleteBondById(int id);
    Q_INVOKABLE void deleteSelection();
    Q_INVOKABLE void undo();
    Q_INVOKABLE void redo();
    Q_INVOKABLE void copySelection();
    Q_INVOKABLE void cutSelection();
    Q_INVOKABLE void pasteSelection(double cx, double cy);
    Q_INVOKABLE void selectAll();
    Q_INVOKABLE void clearCanvas();
    Q_INVOKABLE void loadBenzene();
    Q_INVOKABLE void deserializeMol(const QString& data);
    Q_INVOKABLE void requestStructure(const QString& fmt, const QString& reqId);
    Q_INVOKABLE void requestSelectionStructure(const QString& reqId);
    Q_INVOKABLE void requestSerialize(const QString& reqId);
    Q_INVOKABLE void requestMolfileForNaming(const QString& reqId);
    // Selection-scoped variant: names only the currently selected fragment
    // (ChemDraw's "Structure to Name" behavior -- selection required, only
    // the selection is named, not the whole document). Emits structureReady
    // with an empty string if the selection is empty.
    Q_INVOKABLE void requestMolfileForNamingSelection(const QString& reqId);
    Q_INVOKABLE QString getStructure(const QString& fmt);
    Q_INVOKABLE QString serializeMol();
    Q_INVOKABLE QString getMolfileForNaming();
    Q_INVOKABLE void loadStructure(const QString& format, const QString& data, bool centerOnPage = false);
    Q_INVOKABLE bool importReaction(const QString& text);
    Q_INVOKABLE void insertFunctionalGroup(const QString& fgName, double cx, double cy, int targetAtomId = -1, bool fullStructure = true);
    Q_INVOKABLE void insertLibraryTemplateFused(const QString& fgName, double cx, double cy, int targetBondId);
    Q_INVOKABLE void requestSaltsAndSolventsList();
    Q_INVOKABLE void requestFunctionalGroupsList();
    Q_INVOKABLE void requestTemplateLibraryList();
    Q_INVOKABLE void requestTemplateThumbnail(const QString& name, const QString& reqId);
    Q_INVOKABLE void addRxnArrow(double x, double y, const QString& mode = QStringLiteral("filled-triangle"));
    Q_INVOKABLE void addRxnPlus(double x, double y);
    Q_INVOKABLE void addCurvedArrow(double x1, double y1, double ctrlX, double ctrlY, double x2, double y2);
    Q_INVOKABLE void setRxnArrowMode(int id, const QString& mode);
    Q_INVOKABLE void setRxnArrowConditions(int id, const QString& above, const QString& below);
    Q_INVOKABLE void setStereoFlags(const QString& type, int groupId);
    Q_INVOKABLE void addRGroup(int rgroupNumber);
    Q_INVOKABLE void deleteRGroup(int rgroupNumber);
    Q_INVOKABLE void setRGroupLogic(int rgroupNumber, const QString& range, bool resth, int ifthen);
    Q_INVOKABLE void addRGroupMember(int rgroupNumber);
    Q_INVOKABLE void removeRGroupMember(int rgroupNumber, int fragId);
    Q_INVOKABLE void setAtomQueryList(int atomId, const QString& elementsCsv, bool notList);
    Q_INVOKABLE void clearAtomQueryList(int atomId, const QString& fallbackLabel = QStringLiteral("C"));
    Q_INVOKABLE void addMultitailArrow(double x, double y);
    Q_INVOKABLE void deleteMultitailArrow(int id);
    Q_INVOKABLE void addMultitailArrowTail(int id);
    Q_INVOKABLE void changeAtomLabel(int id, const QString& label);
    Q_INVOKABLE void setAtomMapping(int id, int mapping);
    Q_INVOKABLE void changeBondType(int id, int type, int stereo = 0);
    Q_INVOKABLE void changeAtomCharge(int id, int charge);
    Q_INVOKABLE void setAttachmentPoint(int id, int order);
    Q_INVOKABLE void changeAtomIsotope(int id, int isotope);
    Q_INVOKABLE void changeAtomRadical(int id, int radical);
    Q_INVOKABLE void changeAtomValence(int id, int valence);
    Q_INVOKABLE void requestAtomProperties(int id);
    Q_INVOKABLE void selectByRect(double x1, double y1, double x2, double y2);
    Q_INVOKABLE void addSelectionByRect(double x1, double y1, double x2, double y2);
    Q_INVOKABLE void selectByLasso(const QVariantList& pointsFlat);
    Q_INVOKABLE void selectItem(const QVariant& atomId, const QVariant& bondId, const QVariant& rxnArrowId = QVariant(), const QVariant& rxnPlusId = QVariant(), const QVariant& multitailArrowId = QVariant());
    Q_INVOKABLE void addItemToSelection(const QVariant& atomId, const QVariant& bondId);
    Q_INVOKABLE void removeItemFromSelection(const QVariant& atomId, const QVariant& bondId);
    Q_INVOKABLE void selectFragment(const QVariant& atomId, const QVariant& bondId);
    Q_INVOKABLE void moveSelection(double dx, double dy);
    Q_INVOKABLE void commitMove();
    Q_INVOKABLE void rotateSelectionLive(double angleDelta);
    Q_INVOKABLE void commitRotate();
    Q_INVOKABLE void scaleSelectionLive(double factor, double anchorX, double anchorY);
    Q_INVOKABLE void commitScale();
    Q_INVOKABLE void centerStructure();
    Q_INVOKABLE void normalizeStructure();

    Q_INVOKABLE void alignAtoms(const QString& direction);
    Q_INVOKABLE void distributeAtoms(const QString& direction);
    Q_INVOKABLE void transformSelection(const QString& mode);
    Q_INVOKABLE void addChain(double x1, double y1, double x2, double y2);
    Q_INVOKABLE void addText(const QString& content, double x, double y, bool bold = false, bool italic = false);
    Q_INVOKABLE void updateText(int id, const QString& content, bool bold = false, bool italic = false);
    Q_INVOKABLE void deleteText(int id);
    Q_INVOKABLE void addImage(const QString& base64DataUri, double cx, double cy, double halfW, double halfH);
    Q_INVOKABLE void deleteImage(int id);

    Q_INVOKABLE void setStereoDescriptors(const QString& jsonMap);

    Q_INVOKABLE void setCheckIssues(const QString& jsonMap);

    Q_INVOKABLE void sendCommand(const QString& cmd, const QVariantList& args = QVariantList());

    // OS clipboard helpers (synchronous, no IPC)
    Q_INVOKABLE QString getOsClipboardText() const;
    Q_INVOKABLE void setOsClipboardText(const QString& text);

    // Worker clipboard → OS clipboard
    Q_INVOKABLE void requestClipboardKet();
    // Import KET text as a paste at (cx, cy)
    Q_INVOKABLE void importKetAtPosition(const QString& ket, double cx, double cy);
    
    // Project file serialization
    Q_INVOKABLE void requestProjectSerialization(const QString& reqId);
    Q_INVOKABLE void loadProject(const QUrl& fileUrl);
    Q_INVOKABLE bool saveProject(const QUrl& fileUrl);
    
    // Get current drawing style as JSON
    Q_INVOKABLE QJsonObject getDrawingStyle() const;
    Q_INVOKABLE void setDrawingStyle(const QJsonObject& style);
    
    // Apply drawing style with whole-drawing rescale (for bond length changes)
    // When bondLength changes, all atom coordinates are scaled to maintain visual proportions
    Q_INVOKABLE void applyDrawingStyleWithRescale(const QJsonObject& style);
    
    // Get current page settings as JSON
    Q_INVOKABLE QJsonObject getPageSettings() const;
    Q_INVOKABLE void setPageSettings(const QJsonObject& settings);

signals:
    void stateUpdated(const QVariantMap& state, const QVariantMap& selection, bool isDirty, bool canUndo, bool canRedo, const QVariant& result);
    void primitivesChanged();
    void selectionChanged();
    void overlayStateChanged();
    void structureReady(const QString& reqId, const QString& data);
    void errorOccurred(const QString& error);

private:
    // Turns m_docState's current molecule/selection into the exact QVariantMap shape and
    // signal set QML already consumes from the JS path -- the single place that makes a
    // C++-mode document indistinguishable from a JS-mode one at the QML boundary. Always
    // renders with showExplicitH=false; the setShowExplicitH toggle is out of scope for this
    // pilot (see the approved spec's Scope section).
    void applyLocalState(const QString &resultTag = QString());

    // Ports 40-serialize.js's _buildBatchRecordsFromStructs: builds the {count, records:
    // [{index, label, thumb: {atoms, bonds}}]} JSON shape shared by deserializeRdfBatch/
    // deserializeIndigoBatch/deserializeSdfBatch/realignSdfBatch. `count` is m_sdfBatch's true,
    // uncapped record count; the `records` array itself is capped at the first 500 entries,
    // matching the real code's own Math.min(count, 500) exactly.
    QString buildSdfBatchListJson() const;

    // Shared by deserializeRdfBatch/deserializeIndigoBatch, whose real JS bodies are byte-
    // identical (confirmed by direct source read) -- avoids literally duplicating the branch
    // body twice the way the real JS does. recordsJson is JSON.stringify([{molfile: "..."}, ...]).
    void handleDeserializeBatchFromMolfiles(const QString& recordsJson);

    std::unique_ptr<DocumentState> m_docState;   // always non-null; every document is C++-engine now
    TemplateLibrary* m_templateLibrary = nullptr;   // non-owning; null unless passed in at construction
    QString m_docClipboardMol;   // C++-engine-only clipboard cache: set by copySelection/cutSelection,
                                  // read by pasteSelection. Mirrors the JS worker's module-level
                                  // _clipboard (90-dispatch.js:133-141, 40-serialize.js:120-123).
    SdfBatch m_sdfBatch;          // C++-engine-only staged batch, mirrors the JS worker's module-level
                                  // _sdfBatchRecords (40-serialize.js). Chemistry/Indigo work is
                                  // already fully implemented and tested in SdfBatch itself.
    QHash<QString, QString> m_sdfProps;   // mirrors _sdfProps (40-serialize.js:2, starts as {}),
                                            // populated by loadSdfBatchRecord, read by getSdfProps.

    QVariantMap m_primitives;
    QVariantMap m_selection;
    QVariantMap m_overlayState;
    double getLargestEmptyAngle(const QVariant& atomId, const QVariantMap& primitives);
};
