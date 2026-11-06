// src/app/molecule/DocumentState.h
#ifndef DOCUMENTSTATE_H
#define DOCUMENTSTATE_H

// One per open document: owns the EditableMolecule, the undo/redo history,
// and the current selection. C++ replacement for src/worker/10-state.js's
// module-level _struct/_history/_historyPointer/_selection/_dirty state
// (sub-project 2 of the chem-core.js migration; see
// docs/superpowers/specs/2026-08-01-document-state-cpp-design.md).
//
// executeCommand/undo/redo are a direct line-for-line port of the three real
// JS functions of the same name, confirmed by direct source read: truncate
// any redo tail, run/invert the command, cap history at kHistorySize
// (dropping the oldest and compensating the pointer), clear selection on
// undo/redo.

#include "EditableMolecule.h"
#include "BiopolymerSequenceView.h"
#include "TemplateLibrary.h"
#include "EditCommand.h"
#include "SelectionState.h"
#include "PlacementEngines.h"
#include "DrawingStyle.h"
#include "PageSettings.h"
#include <vector>

class DocumentState {
public:
    explicit DocumentState(const QString& initialStructure = {});

    EditableMolecule& molecule();
    const EditableMolecule& molecule() const;

    void executeCommand(EditCommand cmd);
    void undo();
    void redo();
    bool isDirty() const;
    // Grounded requirement (not speculative): src/worker/30-templates.js's
    // loadBenzene() explicitly does `init(); addBenzeneRing(); _dirty = false`
    // -- a new document pre-populated with the starter benzene is NOT dirty.
    void markClean();
    bool canUndo() const;
    bool canRedo() const;

    // Ports src/worker/90-dispatch.js:147's setShowExplicitH -- a bare module-level variable
    // assignment in the real JS, deliberately NOT wrapped in makeCmd/executeCommand: it's a
    // display preference, not an edit, so it is not undoable and does not mark the document
    // dirty. Read by V8Process::applyLocalState() on every render.
    void setShowExplicitH(bool show);
    bool showExplicitH() const;

    void selectAtom(AtomId id);
    void selectBond(BondId id);
    void selectRxnArrow(RxnArrowId id);
    void selectRxnPlus(RxnPlusId id);
    void selectMultitailArrow(MultitailArrowId id);
    void clearSelection();

    void addAtomToSelection(AtomId id);
    void addBondToSelection(BondId id);
    void addRxnArrowToSelection(RxnArrowId id);
    void addRxnPlusToSelection(RxnPlusId id);
    void addMultitailArrowToSelection(MultitailArrowId id);

    void removeAtomFromSelection(AtomId id);
    void removeBondFromSelection(BondId id);
    void removeRxnArrowFromSelection(RxnArrowId id);
    void removeRxnPlusFromSelection(RxnPlusId id);
    void removeMultitailArrowFromSelection(MultitailArrowId id);

    // Selects every current atom/bond/rxnArrow/rxnPlus/multitailArrow. KNOWN
    // SIMPLIFICATION: the real 10-state.js selectAll() excludes atoms hidden
    // inside a contracted SUP functional-group sgroup (selecting the sgroup's
    // own pill id instead) -- that needs sgroup-contraction awareness
    // (type/isExpanded) EditableMolecule doesn't have yet. This selectAll()
    // selects raw atoms/bonds with no sgroup-contraction logic; revisit once
    // a later sub-project models contracted sgroups.
    void selectAll();
    SelectionState& selection();

    // Ports 10-state.js's selectByRect (646-733). REPLACES the current selection.
    // Sgroup contraction not modeled -- see this file's own selectAll() comment above for
    // why (same already-accepted port-wide simplification).
    void selectByRect(double x1, double y1, double x2, double y2);

    // Ports 10-state.js's addSelectionByRect (817-860). ADDS to the current selection.
    // Deliberately asymmetric with selectByRect: bonds need BOTH endpoints inside (no
    // edge-crossing test), rxnArrows/rxnPluses are not touched at all -- this matches the
    // real function exactly, not a simplification.
    void addSelectionByRect(double x1, double y1, double x2, double y2);

    // Ports 10-state.js's selectByLasso (754-815) / the shared _pointInPolygon helper
    // (735-745). REPLACES the current selection. Clears the selection (no-op) on fewer
    // than 3 points, matching the real "pointsFlat.length < 6" guard (6 = 3 points x,y each).
    void selectByLasso(const QList<QPointF>& points);

    // Ports 10-state.js's selectFragment (424-488): connected-component select via
    // EditableMolecule::atomFragmentIndex/atomIdsInFragment. REPLACES the current
    // selection when atomId resolves to a real atom; CLEARS it when atomId is null
    // (including the "bondId given but doesn't resolve" sub-case -- these produce the
    // SAME outcome in the real code, see this method's own body comment in the .cpp);
    // NO-OPs (leaves selection untouched) when a non-null atomId doesn't exist. Pass -1
    // for "null" (this port's established sentinel). No sgroup-contraction awareness --
    // same accepted simplification as selectAll()/selectByRect() above.
    void selectFragment(AtomId atomId, BondId bondId);

    // Ports 10-state.js's selectRing (490-530) / its shortestRingThroughBond helper
    // (129-156). REPLACES the current selection with the shortest ring through the given
    // bond (or the shortest ring incident to the given atom, ties won by first-found in
    // bondIds() order). CLEARS the selection only when BOTH args are null; every other
    // failure path (no start bond found, or the start bond isn't part of any ring) LEAVES
    // the selection untouched -- deliberately asymmetric with the both-null case, matches
    // the real code exactly. Precedence when both args are given: bondId wins (opposite of
    // selectFragment/selectChain, a genuine inconsistency in the original app, preserved).
    void selectRing(AtomId atomId, BondId bondId);

    // Ports 10-state.js's selectChain (532-622): heavy-atom BFS chain walk, stopping at
    // ring atoms (included but not expanded past, unless the seed itself), branch points
    // (3+ heavy substituents, unconditionally, even at the seed), and terminal atoms
    // (unless the seed). REPLACES the current selection; CLEARS it only when both args are
    // null; LEAVES it untouched when a given bondId doesn't resolve to a real bond. No
    // existence check on a directly-given atomId -- matches the real code exactly (a
    // bogus id selects just itself and stops, since neighborAtomIds gracefully returns
    // empty for an unknown atom).
    void selectChain(AtomId atomId, BondId bondId);

    // Ports 40-serialize.js's selectSubstructureMatches (627-646). matchesJson is
    // {matches: [[1-indexed atom position, ...], ...]}, the exact JSON shape
    // IndigoService::substructureSearch already produces (indigoIndex(tAtom) + 1). REPLACES
    // the current selection with the union of every matched atom across all match arrays,
    // plus any bond whose both endpoints ended up selected. Malformed JSON, a missing
    // "matches" key, an empty matches array, or an out-of-range index all degrade
    // gracefully to selecting nothing / ignoring that one index -- never a crash, matching
    // the real code's own graceful degradation exactly (this also covers the real
    // "clear matches" UI case, an empty matches array).
    void selectSubstructureMatches(const QString& matchesJson);

    // 40-serialize.js's copySelection: extracts the current selection's atoms+bonds as
    // MOL-format text. Pure read -- pushes NO history entry (matching the real function, which
    // performs no mutation). Returns an empty string if the selection has no atoms (matching the
    // real `if (_selection.atom_ids.length === 0) return null` guard) or if extraction fails.
    // Deliberately does NOT include any selected rxnArrows/rxnPluses/multitailArrows -- the real
    // _struct.clone(atomSet, bondSet) never touches those entity types either.
    QString copySelection() const;

    // 10-state.js's deleteSelection (862-961) / this migration's own sgroup-aware-deletion
    // gap, deferred since sub-project 3a and explicitly assigned here. Deliberately has NO
    // "empty selection = no-op" guard -- matches the real function's own actual behavior
    // (an empty selection is a harmless no-op command, not skipped entirely). See
    // docs/superpowers/specs/2026-08-03-delete-selection-cpp-design.md for the full design.
    void deleteSelectionEntities();

    // 40-serialize.js's cutSelection (96-102): copySelection() then delete the selection if
    // the copy produced real content. The delete itself is what pushes the undo entry --
    // copySelection is a pure read with no history entry, per its own established contract.
    QString cutSelection();

    // 60-analysis.js's setStereoDescriptors: applies an external CIP-perception result onto
    // stored fields. NOT undoable (no makeCmd in the real function either) -- this applies an
    // analysis result, not a user edit. See
    // docs/superpowers/specs/2026-08-03-analysis-results-cpp-design.md for the full design,
    // especially why positional indices resolve through atomIdsInIndigoOrder(), not atomIds().
    void setStereoDescriptors(const QString& jsonMap);

    // 60-analysis.js's setCheckIssues: applies external structure-check results onto stored
    // checkWarning fields. NOT undoable, same reasoning as setStereoDescriptors above.
    void setCheckIssues(const QString& jsonMap);

    // 30-templates.js's clearCanvas: replaces the document with a blank one, undoable.
    // Delegates to deserializeMol -- verified equivalent (indigoLoadMoleculeFromString("")
    // succeeds; deserializeMol has no zero-atom guard). See
    // docs/superpowers/specs/2026-08-04-clear-canvas-load-benzene-cpp-design.md.
    void clearCanvas();

    // 30-templates.js's loadBenzene: init() + addBenzeneRing(4,4) + _dirty=false. NOT undoable
    // -- wipes the entire undo/redo history (a first for this port; every other command only
    // appends to or truncates history, never clears it wholesale). See
    // docs/superpowers/specs/2026-08-04-clear-canvas-load-benzene-cpp-design.md.
    void loadBenzene();

    // src/worker/70-biopolymer.js's 4 functions, ported via BiopolymerSequenceView (sub-project
    // 6d). None are undoable -- the real functions have no makeCmd/executeCommand either. See
    // docs/superpowers/specs/2026-08-04-biopolymer-sequence-view-cpp-design.md.
    void buildBioSequenceView(const QString& sequenceText, const QString& seqType);
    void addBioMonomer(const QString& symbol, const QString& seqType = QString());
    void deleteBioMonomer(int id);
    QList<BioMonomer> bioMonomers() const;
    QList<BioBond> bioBonds() const;
    QString bioSeqType() const;

    AtomId addAtom(const QString& symbol, double x, double y);
    BondId addBond(AtomId a, AtomId b, int order, int stereo = 0);
    void deleteAtom(AtomId id);   // plain-atom-and-incident-bonds case only; see the
                                  // sgroup-membership spike documented in EditableMolecule.h
                                  // for why the sgroup-aware branch is deferred
    void deleteBond(BondId id);

    void changeAtomLabel(AtomId id, const QString& newLabel);
    void setAtomMapping(AtomId id, int mappingNumber);
    void changeAtomCharge(AtomId id, int newCharge);
    void setAttachmentPoint(AtomId id, int order);
    void changeAtomIsotope(AtomId id, int isotope);
    void changeAtomRadical(AtomId id, int radical);
    void changeAtomValence(AtomId id, int valence);

    struct AtomProperties { QString label; int charge; int isotope; int radical; int explicitValence; };
    AtomProperties atomProperties(AtomId id) const;

    void changeBondOrder(BondId id, int newOrder);

    // Ports 20-edit.js's changeBondType (586-601): changes bond order AND stereo as ONE
    // atomic undo step (real JS: a single makeCmd covering both fields). newStereo is the RAW
    // V2000 code (0/1/4/6, matching the real JS's own b.stereo field convention) -- NOT an
    // EditableMolecule::Direction -- mapped internally to avoid exposing EditableMolecule's
    // private V2000<->Indigo conversion helpers to this unrelated class. execute()/invert()
    // both set order FIRST, then attempt the stereo call: EditableMolecule::setBondStereo's own
    // guard (bondOrder(id) != 1) makes that a safe no-op whenever the order just became
    // non-single, so an original wedge on a bond later changed away from order 1 is still
    // exactly restored on undo (order goes back to 1 first, THEN the stereo call succeeds).
    void changeBondTypeAndStereo(BondId id, int newOrder, int newStereo);

    // Sets or clears a single bond's wedge/hash/either stereo direction (sub-project 8). See
    // EditableMolecule::setBondStereo's own header comment for why this works via a molfile
    // round-trip rather than a live per-bond setter, and for the pre-existing double-bond and
    // caller-responsibility (at-most-one-wedge-per-stereocenter) limitations that apply here too.
    // No history entry is pushed if the underlying call would be rejected (unknown bond id, a
    // non-single bond, or an unchanged value) -- same guard-before-executeCommand convention as
    // changeAtomCharge/changeBondOrder above.
    void setBondStereo(BondId id, EditableMolecule::Direction direction);

    void setAtomQueryList(AtomId id, const QString& labelsCsv, bool notList);
    // KNOWN LIMITATION: undo after clearAtomQueryList does NOT restore the
    // original query list (EditableMolecule has no set-list-from-numbers
    // entry point bypassing the CSV parse) -- it only re-clears on redo.
    // Revisit if a later sub-project needs this restored exactly.
    void clearAtomQueryList(AtomId id, const QString& fallbackLabel);

    // Discrete transforms of the selected ATOMS about their unweighted
    // centroid -- a direct port of 20-edit.js's transformSelection(mode),
    // split into four named methods instead of its JS mode-string dispatch.
    // Each no-ops below 2 selected atoms, matching the real ids.length < 2
    // guard, and pushes exactly one undoable command.
    void rotateSelection90CW();
    void rotateSelection90CCW();
    void flipSelectionHorizontal();
    void flipSelectionVertical();

    // Per-image transforms (20-edit.js has none; these port 50-reactions.js's
    // moveImage/resizeImage). Take an explicit id rather than reading the
    // selection -- images are deliberately not a SelectionState entity type.
    // Both commit immediately as one undoable command; there is no live-drag
    // protocol for images in the real code either.
    void moveImage(ImageId id, double dx, double dy);
    void resizeImage(ImageId id, double scaleFactor);

    // Image lifecycle (50-reactions.js: addImage/deleteImage). addImage no-ops (returns -1) on
    // empty image data, matching the real `if (!base64DataUri) return`.
    ImageId addImage(const QByteArray& pngData, double cx, double cy, double halfW, double halfH);
    void deleteImage(ImageId id);

    // Rxn-arrow lifecycle (50-reactions.js: addRxnArrow/addCurvedArrow/setRxnArrowMode/
    // setRxnArrowConditions/deleteRxnArrow). addRxnArrow computes its default second endpoint
    // the same way the real JS does (2.5 bond-lengths to the right). deleteRxnArrow's undo
    // recreates a fresh-id arrow with the same mode/conditions/curvature, matching this class's
    // established deleteAtom precedent (fresh id, not the original).
    RxnArrowId addRxnArrow(double cx, double cy, const QString& mode = QStringLiteral("filled-triangle"));
    RxnArrowId addCurvedArrow(double x1, double y1, double ctrlX, double ctrlY, double x2, double y2);
    void setRxnArrowMode(RxnArrowId id, const QString& newMode);
    void setRxnArrowConditions(RxnArrowId id, const QString& above, const QString& below);
    void deleteRxnArrow(RxnArrowId id);

    // Rxn-plus lifecycle (50-reactions.js: addRxnPlus/deleteRxnPlus).
    RxnPlusId addRxnPlus(double cx, double cy);
    void deleteRxnPlus(RxnPlusId id);

    // Multitail-arrow lifecycle (50-reactions.js: addMultitailArrow/deleteMultitailArrow/
    // addMultitailArrowTail). addMultitailArrow starts with ZERO tails, matching the real
    // tailsYOffset=[] -- storing just the head point is already a valid state per
    // RenderPrimitives.cpp's own `if (pts.size() < 2) continue` guard. addMultitailArrowTail
    // ports the real widest-gap insertion algorithm, adapted to this port's flat, absolute-
    // coordinate storage (no separate height/tailsYOffset fields exist here).
    MultitailArrowId addMultitailArrow(double cx, double cy);
    void deleteMultitailArrow(MultitailArrowId id);
    void addMultitailArrowTail(MultitailArrowId id);

    // Document-level stereo-display toggle (50-reactions.js's setStereoFlags). Defaults type to
    // "abs" and groupId to 0 exactly like the real function's `type || 'abs', groupId || 0`.
    void setStereoFlags(const QString& type, int groupId);

    // R-groups (50-reactions.js: addRGroup/deleteRGroup/setRGroupLogic/addRGroupMember/
    // removeRGroupMember). addRGroupMember reads the current selection's atoms, maps each to its
    // live fragment index (EditableMolecule::atomFragmentIndex), and registers only the fragments
    // not already members -- matching the real function's diff-then-guard exactly.
    bool addRGroup(int rgroupNumber);
    bool deleteRGroup(int rgroupNumber);
    void setRGroupLogic(int rgroupNumber, const QString& range, bool resth, int ifthen);
    void addRGroupMember(int rgroupNumber);
    void removeRGroupMember(int rgroupNumber, int fragId);

    // Text annotations (50-reactions.js: addText/updateText/deleteText). addText no-ops
    // (returns -1) on an empty or whitespace-only string, matching the real
    // `if (!plainStr || !String(plainStr).trim()) return`.
    TextId addText(const QString& plainStr, double x, double y, bool bold, bool italic);
    void updateText(TextId id, const QString& plainStr, bool bold, bool italic);
    void deleteText(TextId id);

    // Bracket selection (50-reactions.js's addBracketSelection): bbox over the current
    // selection's atoms (and bond endpoints), padded by 0.8, pushed onto EditableMolecule's
    // bracket stack. No-ops if the selection is empty. Undo pops the bracket (stack semantics,
    // no id needed).
    void addBracketSelection();

    // Whole-document replace from a MOL-format string (50-reactions.js's sibling in
    // 40-serialize.js: deserializeMol via _applyLoadedStruct). No-ops (no history entry) on a
    // parse failure. Extension data (name, texts, rxnArrows, etc.) is wiped on success, matching
    // the real JS exactly -- it never carries anything across a document replace.
    void deserializeMol(const QString& data, bool centerOnPage = false);

    // Document name (40-serialize.js's setMoleculeName). No-op if unchanged, matching the real
    // JS's `if (oldName === newName) return`. getMoleculeName needs no wrapper -- pure read,
    // already covered by EditableMolecule::name().
    void setMoleculeName(const QString& name);

    // Document-level drawing style
    const DrawingStyle& drawingStyle() const { return m_drawingStyle; }
    void setDrawingStyle(const DrawingStyle& style);
    
    // Apply drawing style with whole-drawing rescale (for bond length changes)
    // When bondLength changes, all atom coordinates are scaled by (newBondLength/oldBondLength)
    // around the drawing's centroid. This is a single undoable operation.
    // If bondLength hasn't changed, just applies the style without rescaling.
    void applyDrawingStyleWithRescale(const DrawingStyle& style);
    
    // Document-level page settings
    const PageSettings& pageSettings() const { return m_pageSettings; }
    void setPageSettings(const PageSettings& settings);
    
    // Access the fixed physical scale for model-to-mm conversion
    double modelToMmScale() const { return m_drawingStyle.modelToMmScale; }

    // ---- Live-drag gestures ------------------------------------------------
    // Repeated *Live calls during a UI drag mutate positions directly and push
    // NOTHING onto the undo history; the matching commit* call at gesture end
    // pushes exactly one undoable command covering the whole gesture. Direct
    // port of 20-edit.js's moveSelection/commitMove pair.
    //
    // NOTE the deliberate per-gesture asymmetry, ported exactly and NOT
    // normalized: move takes an INCREMENTAL (dx,dy) applied to current
    // positions; rotate takes an INCREMENTAL angle accumulated into a running
    // total re-applied from the original snapshot; scale takes an ABSOLUTE
    // factor measured from gesture start. See each method's port note.
    void moveSelectionLive(double dx, double dy);
    void commitMove();

    // Rotate about the centroid of the snapshotted selection points.
    // angleDelta is INCREMENTAL (radians); the running total is re-applied to
    // the ORIGINAL snapshot every call, which is what keeps a long drag free
    // of floating-point drift. Multitail arrows are excluded -- see the
    // limitation note on snapshotSelectionPoints below.
    void rotateSelectionLive(double angleDelta);
    void commitRotate();

    // Uniform scale about a caller-supplied FIXED anchor (the dragged handle's
    // geometric opposite on the bbox), NOT the centroid -- matching real
    // corner-drag resize behaviour and deliberately unlike rotate.
    //
    // NOTE: `factor` is ABSOLUTE (total scale measured from gesture start), not
    // incremental like rotate's angleDelta. Calling this with 2.0 then 3.0
    // yields 3x from the original, not 6x. Ported exactly from the real code.
    void scaleSelectionLive(double factor, double anchorX, double anchorY);
    void commitScale();

    // Structure-cleanup commands ported from src/worker/20-edit.js:753-834. Both operate on
    // the current selection's atoms only (never rxn-arrows/pluses/etc, unlike commitMove/
    // commitScale) -- confirmed by direct read of the real JS, which only ever touches
    // _selection.atom_ids for these two commands.
    void alignAtoms(const QString& direction);
    void distributeAtoms(const QString& direction);

    // Ports src/worker/20-edit.js:836-956's layoutSelectedChain. Re-positions only the
    // selected atoms, leaving everything else fixed. 4 validation gates, each a silent no-op
    // on failure (never a partial/best-effort layout) -- see largestEmptyAngleAt below and the
    // spec's precision note on gate 3 (a clean cycle is NOT rejected, only a genuinely
    // disconnected sub-group is).
    void layoutSelectedChain();

    // Generic ring construction, ported from 30-templates.js's addRing +
    // perceiveRingAlternation. `coords` is a flat [x0,y0, x1,y1, ...] list in
    // exactly the caller's winding order -- NEVER canonicalized, since the
    // rendering layer already handles either winding direction
    // (10-state.js:1171-1177). No-ops on malformed input (odd length or
    // fewer than 3 points). `aromatic` is a plain bool, not a tri-state type:
    // the only reachable path into this method (V8Process::addRing) always
    // passes a real bool, so the JS's unreachable `undefined` case needs no
    // representation here.
    void addRing(const QList<double>& coords, bool aromatic = true);

    // Zig-zag chain tool, ported from 20-edit.js's addChain. Both endpoints
    // are clamped to the page bounds FIRST, before any geometry is derived.
    // Geometry matches PlacementEngines.cpp's ChainPlacementEngine::compute
    // exactly (same nBonds/theta/half-angle formula) -- that class is the
    // live-drag preview and is deliberately NOT reused here (see its own
    // header comment: duplicating the merge step there risks preview and
    // committed structure drifting apart). Reuses mergeOverlappingAtoms so
    // endpoints graft onto whatever they overlap, exactly like addRing.
    void addChain(double x1, double y1, double x2, double y2);

    // Ports 20-edit.js's addBondAndAtom (605-669): drag from an existing atom to a new point. If
    // the target point collides with a DIFFERENT existing atom (within 0.3 units), no atom is
    // created -- just a single addBond call to that atom. Otherwise creates ONE new atom (clamped
    // to page bounds) and ONE bond to it, as a single atomic undoable command (following the
    // insertStructureAt/addAtom idBox pattern already established in this class). `stereo` is
    // accepted for call-site parity with the real command but is a documented no-op -- this port
    // has no way to SET bond stereo anywhere except via EditableMolecule::setBondStereo's own
    // molfile-round-trip mechanism (sub-project 8), which is a different, separate entry point
    // not reused here. The real JS's short-drag ring-angle/largest-empty-angle auto-placement
    // (20-edit.js:616-633, triggered when the drag distance is under 0.5 units) is out of scope
    // -- no equivalent geometry helper exists in this port yet; the caller's raw endpoint is
    // always used, even for a very short drag.
    void addBondAndAtom(AtomId startId, const QString& label, double x, double y, int type, int stereo);

    // Ports 20-edit.js's addBondBetweenCoords (671-702): drag between two empty points. No
    // collision/merge check at all (the real function's own comment: "Creates a bond between two
    // empty coordinate points"). Always creates two new "C" atoms (each clamped to page bounds
    // independently) and one bond between them, as a single atomic undoable command. Same
    // documented bond-stereo no-op as addBondAndAtom above.
    void addBondBetweenCoords(double x1, double y1, double x2, double y2, int type, int stereo);

    // Ports 10-state.js's selectItem (390-404): replaces the WHOLE selection with at most one
    // item, priority order atom > bond > rxnArrow > rxnPlus > multitailArrow > none (all unset
    // clears the selection). Not undoable -- matches the real function (no makeCmd) and this
    // class's own selectAtom/selectBond etc. -1 is this class's established "unset" sentinel
    // (matching insertFunctionalGroup's targetAtomId = -1).
    void selectSingleItem(AtomId atomId = -1, BondId bondId = -1, RxnArrowId rxnArrowId = -1,
                           RxnPlusId rxnPlusId = -1, MultitailArrowId multitailArrowId = -1);

    // Ports 30-templates.js's insertFunctionalGroup. TemplateLibrary is
    // passed explicitly by the caller (this class has no module-level
    // template-registry globals, unlike the real worker). graft = targetAtomId
    // resolves to a real atom AND the template has a SUP attachment point;
    // grafting fuses the template's attach atom onto targetAtomId via
    // graftAtomOnto rather than a position-coincidence merge. Falls back to
    // a single placeholder atom labeled fgName if the template is missing or
    // structurally empty. One undoable command per call.
    void insertFunctionalGroup(const TemplateLibrary& lib, const QString& fgName,
                                double cx, double cy, AtomId targetAtomId = -1, bool fullStructure = true);

    // Non-mutating sibling of insertFunctionalGroup: computes the exact same placement geometry
    // a commit would produce (translate-only when free-floating, or the same rotate-then-
    // translate graft transform when targetAtomId is a real atom) without touching m_molecule or
    // pushing any undo history. Used by both the hover-ghost preview and the drag-ghost preview
    // for FG_/SS_/LIB_ tools, so neither can ever drift out of sync with what actually commits.
    // Returns PlacementResult{valid=false} if the template can't be found or is structurally
    // empty (mirrors insertFunctionalGroup's own single-placeholder-atom fallback condition, but
    // this preview method does not attempt to synthesize a placeholder ghost for that case).
    PlacementResult getFunctionalGroupPreview(const TemplateLibrary& lib, const QString& fgName,
                                               double cx, double cy, AtomId targetAtomId);

    // Shared core for BOTH pasteSelection (source = copySelection()'s output) and
    // insertRecognizedStructure (source = an externally-recognized molfile string) -- the real
    // 30-templates.js's _insertStructAt is exactly this shared helper already. Clamps (cx, cy) to
    // the page bounds first. No-ops (no history entry) on an empty or unparseable sourceMolfile.
    void insertStructureAt(const QString& sourceMolfile, double cx, double cy);

    // Imports a full reaction ($RXN-format or reaction-SMILES ">>" text, both auto-detected
    // by indigoLoadReactionFromString) as reactant fragments + a RxnArrow + RxnPlus entities
    // between multiples, laid out left to right. Built as ONE hand-constructed EditCommand
    // (mirroring insertStructureAt's own internal pattern -- calling EditableMolecule methods
    // directly, not the DocumentState-level insertStructureAt/addRxnArrow/addRxnPlus wrappers,
    // each of which would otherwise push its own separate undo entry) so the whole import
    // undoes as one step. Returns false (no history entry, no partial state) if text can't be
    // parsed as a reaction at all.
    bool importReaction(const QString& text);

    // Ports 30-templates.js's insertLibraryTemplateFused. Maps the template's
    // designated fusion bond onto targetBondId via a 2-point similarity
    // transform, then reuses the EXISTING mergeOverlappingAtoms +
    // applyRingAlternation pipeline (sub-project 3c) unchanged -- library.sdf
    // templates have zero SUP groups (confirmed), so this path never touches
    // sgroup logic. Falls back to insertFunctionalGroup(lib, fgName, cx, cy)
    // with NO target atom on any precondition failure (missing fusion
    // metadata, missing target bond, template bond not in a ring).
    void insertLibraryTemplateFused(const TemplateLibrary& lib, const QString& fgName,
                                     double cx, double cy, BondId targetBondId);

    // Ports 30-templates.js's toggleSgroupExpanded. No-ops if id doesn't
    // resolve to a real sgroup (matches the real function's early return).
    void toggleSgroupExpanded(SGroupId id);

    // New capability (sub-project 25, not a port -- neither engine previously supported
    // renaming a contracted group's abbreviation). No-ops if id doesn't resolve to a real
    // sgroup, or if newLabel matches the current label (mirrors toggleSgroupExpanded's own
    // no-op-on-invalid-id convention and setMoleculeName's own unchanged-value guard).
    void renameSGroup(SGroupId id, const QString& newLabel);

private:
    // Shared placement-computation core for insertFunctionalGroup (commits) and
    // getFunctionalGroupPreview (read-only preview) -- everything insertFunctionalGroup used to
    // compute inline before constructing its EditCommand, now returned as plain data so both
    // callers can use it without duplicating ~140 lines of delicate Indigo-session-boundary-
    // sensitive math. See insertFunctionalGroup's own body for why the ordering of raw indigo*
    // calls against fgHandle (TemplateLibrary's own session) relative to calls into mol
    // (EditableMolecule's session) matters -- this helper preserves that ordering exactly as it
    // existed before extraction.
    struct FunctionalGroupPlacement {
        bool valid = false;
        QString fgMolfile;
        int templateAttachIdx = -1;
        bool graft = false;
        bool useRotation = false;
        double cosR = 1.0, sinR = 0.0;
        double attachX = 0, attachY = 0, targetX = 0, targetY = 0;
        double attachTargetX = 0, attachTargetY = 0;
        double dx = 0, dy = 0;
    };
    FunctionalGroupPlacement computeFunctionalGroupPlacement(
        const TemplateLibrary& lib, const QString& fgName, double cx, double cy, AtomId targetAtomId);

    // Recenters the whole document's atom bounding box (plus any reaction arrows/pluses)
    // on the origin. Ports the real JS's loadMolfile(molStr, centerOnPage) recentering
    // step (src/worker/40-serialize.js:313-325) -- SMILES/InChI layout results and
    // Ketcher-derived write-back ops (Aromatize, Normalize, etc.) arrive with arbitrary,
    // often off-page coordinates; genuine file-opens skip this (centerOnPage=false).
    // No-op if the document has zero positioned atoms (matches real JS's
    // `loaded.atoms.size > 0` guard).
    void centerMoleculeOnOrigin();

    enum class DiscreteTransform { RotateCW, RotateCCW, FlipH, FlipV };
    void applyDiscreteTransform(DiscreteTransform mode);

    // DocumentState's own equivalent of V8Process::getLargestEmptyAngle (v8_process.cpp:768-808),
    // built directly on EditableMolecule instead of rendered QVariantMap primitives. That
    // existing method stays untouched -- this is an independent port for layoutSelectedChain's
    // placement step, not a refactor of it. Same algorithm: 0 neighbors -> 0; 1 neighbor ->
    // that neighbor's angle + 2.61799 (fixed ~150 degree offset); 2+ neighbors -> largest
    // circular gap between consecutive sorted neighbor angles, midpoint of that gap.
    double largestEmptyAngleAt(AtomId id) const;

    static constexpr double kBondLength = 1.5;   // chem-core's StandardBondLength

    // Multitail-arrow creation-time constants, computed from DrawingStyle.
    // These are now instance methods that use m_drawingStyle.bondLength
    double getMultitailHeadOffsetX() const { return 2.0 * m_drawingStyle.bondLength; }
    double getMultitailHeight() const { return 3.0 * m_drawingStyle.bondLength; }
    double getMultitailTailInset() const { return 4.0 * m_drawingStyle.bondLength; }

    // Shared by addRing's alternation step: one bond as it stood BEFORE this
    // operation began, keyed by its endpoint pair. Presence of a pair here
    // identifies a fusion seam.
    struct OldBondType { AtomId a, b; int order; };
    void applyRingAlternation(const QList<AtomId>& ringAtoms, const QList<OldBondType>& oldBonds, bool aromatic);

    // Page bounds, copied verbatim from 10-state.js:238-239.
    static constexpr double kPageMinX = -30.0, kPageMaxX = 30.0;
    static constexpr double kPageMinY = -21.0, kPageMaxY = 21.0;

    // Move-gesture drag state (mirrors the real _drag* module-level vars).
    double m_dragDeltaX = 0.0, m_dragDeltaY = 0.0;
    QList<AtomId> m_dragAtomIds;
    QList<RxnArrowId> m_dragArrowIds;
    QList<RxnPlusId> m_dragPlusIds;
    QList<MultitailArrowId> m_dragMultitailIds;
    bool m_dragHasOrigBBox = false;
    double m_dragBBoxMinX = 0.0, m_dragBBoxMinY = 0.0, m_dragBBoxMaxX = 0.0, m_dragBBoxMaxY = 0.0;

    void applyMoveDelta(const QList<AtomId>& atomIds, const QList<RxnArrowId>& arrowIds,
                        const QList<RxnPlusId>& plusIds, const QList<MultitailArrowId>& mtaIds,
                        double dx, double dy);
    void resetMoveDragState();

    // Clears move, rotate AND scale gesture state. DELIBERATE DEVIATION from
    // 10-state.js, which resets only four move vars in undo() and nothing at
    // all in redo(): stale rotate/scale snapshot state surviving an undo lets a
    // later live call re-apply a transform computed from positions the undo has
    // already invalidated, corrupting coordinates. Called by undo() and redo().
    void resetAllDragState();
    void reconcileSelectionAfterCommand();

    // One snapshotted, transformable point. Shared by the rotate and scale
    // gestures (move does not use it -- it works from raw id lists, exactly
    // like the real moveSelection). Mirrors the real _snapshotSelectionPoints
    // tagged-point records; an RxnArrow contributes two independent points.
    //
    // KNOWN LIMITATION (deliberate): multitail arrows are NOT represented here.
    // chem-core's rotate/scale transform only a multitail arrow's spineTop
    // anchor point, and our MultitailArrow has no anchor -- only absolute
    // points -- so there is no faithful mapping. They participate in move
    // (translation is representation-independent) and are excluded from
    // rotate/scale rather than given invented rigid-body behaviour.
    struct TransformPoint {
        enum class Kind { Atom, RxnArrowP1, RxnArrowP2, RxnPlus };
        Kind kind;
        int id;
        double x, y;
    };
    QList<TransformPoint> snapshotSelectionPoints() const;
    static void writeTransformPoint(EditableMolecule& mol, const TransformPoint& pt, double nx, double ny);

    // Rotate-gesture drag state (mirrors the real _rotateDrag* vars).
    QList<TransformPoint> m_rotateOrigPos;
    double m_rotateCenterX = 0.0, m_rotateCenterY = 0.0;
    double m_rotateTotalAngle = 0.0;

    // Scale-gesture drag state (mirrors the real _scaleDrag* vars).
    QList<TransformPoint> m_scaleOrigPos;
    double m_scaleAnchorX = 0.0, m_scaleAnchorY = 0.0;
    double m_scaleTotalFactor = 1.0;

    // Pure geometry/graph helpers ported from 10-state.js, used only by
    // insertLibraryTemplateFused.
    static std::function<QPointF(double, double)> makeSimilarityTransform(
        double p1x, double p1y, double p2x, double p2y, double q1x, double q1y, double q2x, double q2y);
    static QList<int> shortestRingThroughBond(int moleculeHandle, int bondIdx);   // atom indices, or empty if none
    double chooseEmptySide(double qax, double qay, double qbx, double qby, AtomId excludeA, AtomId excludeB,
                           double cursorX, double cursorY) const;

    EditableMolecule m_molecule;
    BiopolymerSequenceView m_bioView;
    std::vector<EditCommand> m_history;
    int m_historyPointer = -1;
    bool m_inCommand = false;
    bool m_dirty = false;
    bool m_showExplicitH = false;
    SelectionState m_selection;
    DrawingStyle m_drawingStyle;
    PageSettings m_pageSettings;
    static constexpr int kHistorySize = 50;
};

#endif // DOCUMENTSTATE_H
