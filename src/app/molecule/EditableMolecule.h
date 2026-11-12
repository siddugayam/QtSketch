// src/app/molecule/EditableMolecule.h
#ifndef EDITABLEMOLECULE_H
#define EDITABLEMOLECULE_H

// Long-lived editable wrapper around one Indigo molecule handle -- the C++
// replacement for chem-core.js's Struct (sub-project 1 of the migration; see
// docs/superpowers/specs/2026-08-01-editable-molecule-cpp-design.md).
//
// ID design: public AtomId/BondId are monotonic, NEVER-reused ints mapped to
// Indigo's internal indices via QHash. Chosen regardless of Indigo's own
// index behavior so that IDs held in undo snapshots can never alias a
// later-created atom.
// Spike findings (tests/editable_molecule_test.cpp, Test 0):
//   "Spike finding (tests/editable_molecule_test.cpp Test 0): surviving atoms keep their
//   pre-removal indices after another atom is removed, BUT a newly-added atom CAN reuse a
//   just-freed index (observed: freed index 2 was reused by the next add). This is exactly
//   why AtomId/BondId are a separate never-reused counter, not raw Indigo indices -- reuse
//   would let an old undo-snapshot's ID silently alias a different, newer atom."
// Verified hard requirement: indigoClone preserves atom/bond indices.
//
// Query atom lists (spec risk 3): plain molecule handles REJECT V2000 atom
// lists outright -- indigoLoadMoleculeFromString returns an error ("atom
// lists are allowed only for queries"); only indigoLoadQueryMoleculeFromString
// accepts them (spike: editable_molecule_test Test 6). Consequence: query
// atom lists cannot live on the plain Indigo handle this class wraps -- they
// must join ExtensionData (or a separate query-molecule-backed type) in a
// later sub-project when query-feature editing is designed.
//
// SGroup-atom-membership (sub-project 3a design question): spike findings
// (tests/editable_molecule_test.cpp Test 11, alanine-free CCO + a 2-atom data
// sgroup): "indigoIterateAtoms(sgroupHandle) iterator=valid, atom count=2" --
// calling indigoIterateAtoms on the SGROUP handle itself (not the molecule
// handle) DOES work and enumerates exactly that sgroup's member atoms, so
// sgroup-atom membership IS queryable through this path, contrary to this
// header's original open question. "indigoRemove(sgroupHandle) returned 1;
// atom count before=3 after=3" -- removing a sgroup handle succeeds and only
// detaches the sgroup record; the molecule's atom count is unchanged, i.e.
// indigoRemove does NOT cascade-delete a sgroup's member atoms. Consequence
// for Task 7's deleteAtom: deleting a sgroup-member atom does not need to
// special-case sgroup cleanup for atom survival (removing the atom directly
// via the existing removeAtom path is safe and does not corrupt the sgroup
// record's atom list on its own), but this class does NOT yet update a
// sgroup's atom-index list when one of its member atoms is removed via
// removeAtom -- that bookkeeping gap is a documented, explicit limitation
// carried into a later sub-project, not attempted here.
//
// Bond stereo (sub-project 3a design question): spike findings
// (tests/editable_molecule_test.cpp Test 12, alanine CC(N)C(=O)O loaded from
// a plain SMILES with no 2D coordinates): "stereocenter type before
// indigoAddStereocenter: 0" (none detected on load, as expected -- SMILES
// carries no wedge/dash). "indigoAddStereocenter(ABS) returned 1" (the call
// to force a specific tetrahedral arrangement at the stereocenter succeeds).
// "indigoMarkStereobonds returned 0" and "bonds with nonzero
// indigoBondStereo() after marking: 0" -- despite the forced stereocenter,
// NOT ONE incident bond ended up with a nonzero (wedge/dash) value from the
// read-only indigoBondStereo() query. Consequence: indigoAddStereocenter +
// indigoMarkStereobonds does NOT translate into a settable per-bond
// wedge/dash through this path on a molecule with no 2D layout to derive a
// direction from (plausible root cause: wedge/dash assignment needs existing
// atom coordinates to pick which bond to mark, which this spike's molecule
// never had). Bond stereo (wedge/dash) is therefore NOT supported by
// changeBondOrder in this sub-project; it stays an explicit, documented gap
// until serialization (sub-project 5) or a dedicated investigation with a
// coordinate-bearing molecule finds a working mechanism.
//
// Bond stereo INVERT (sub-project 3b design question): spike findings
// (tests/editable_molecule_test.cpp Test 14, a V2000 molfile WITH real 2D
// coordinates AND an explicit wedge flag in bond-block column 4 -- the setup
// 3a's Test 12 lacked):
//   "indigoBondStereo(bond0) as loaded = 5"
//   "indigoInvertStereo(bond0) returned -1 (REJECTED)"
//   "invert error: core: indigoInvertStereo: not a stereobond"
//   "indigoBondStereo(bond0) after 1 invert = 5"
//   "indigoBondStereo(bond0) after 2 inverts = 5 (round-trip OK)"
// Two consequences, both load-bearing:
// (1) A wedge bond IS readable -- unlike 3a's coordinate-less SMILES case,
//     indigoBondStereo reports a nonzero value once the molecule has real 2D
//     coordinates and an explicit molfile wedge flag. But the value is
//     Indigo's OWN enum, not the MDL molfile code: INDIGO_UP = 5 and
//     INDIGO_DOWN = 6 (indigo.h:439-440), whereas MDL/chem-core.js use 1 for
//     wedge-up and 6 for hash-down. The two encodings collide on 6 and
//     disagree on up, so any future code bridging them MUST translate rather
//     than pass the number through.
// (2) indigoInvertStereo REJECTS a wedge bond outright ("not a stereobond" --
//     it covers cis/trans stereobonds and stereocenters, not wedge/dash
//     direction). Combined with indigoBondStereo being read-only and 3a's
//     Test 12 ruling out the stereocenter-inference path, there is NO way to
//     change a bond's wedge/dash through the Indigo C API as vendored here.
// Therefore: DocumentState's flip transforms apply the coordinate mirror ONLY.
// chem-core.js's transformSelection also swaps wedge stereo (1 <-> 6) for
// bonds fully inside the selection so the depiction stays chemically
// consistent; that swap is an explicit, documented gap, not an oversight.
//
// Atom-merge mechanics (sub-project 3c design question): spike findings
// (tests/editable_molecule_test.cpp Test 16): (1) a raw indigoRemove(atom)
// DOES auto-remove that atom's incident bonds -- bond count dropped from 1
// to 0 after removing one of the bonded pair, with no separate bond-removal
// call. (2) In the ordinary (non-seam) case, creating the replacement bond
// onto the kept atom BEFORE removing the doomed atom succeeds, and the
// order passed to indigoAddBond is read back unchanged (no surprise there).
// (3) In the seam case -- kept already bonded to the same target -- creating
// the would-be duplicate bond is REJECTED ("already have edge between
// vertices"), exactly 3b's GT-T5 finding, and the pre-existing seam bond's
// order is left completely unchanged by the attempt or by the subsequent
// removal of the doomed atom. Consequence: mergeOverlappingAtoms can
// unconditionally try to create each replacement bond BEFORE removing the
// doomed atom, treat a negative return as the expected seam-dedupe outcome
// (not an error), and rely on removeAtom's cascade to clean up whatever of
// doomed's original bonds remain (either the ones that were never
// successfully rewired, or the ones that were, since the original bond
// object is a separate one from the replacement).

#include <QString>
#include <QHash>
#include <QList>
#include <QMap>
#include <QPointF>
#include <functional>
#include "ExtensionData.h"
#include "MoleculeSnapshot.h"

using AtomId = int;
using BondId = int;
using SGroupId = int;

struct StringResult {
    bool success = false;
    QString value;
    QString error;
};

class EditableMolecule {
public:
    // initialStructure: molfile or SMILES (indigoLoadMoleculeFromString
    // auto-detects); empty string = new blank molecule.
    explicit EditableMolecule(const QString& initialStructure = {});
    ~EditableMolecule();

    EditableMolecule(const EditableMolecule&) = delete;
    EditableMolecule& operator=(const EditableMolecule&) = delete;

    bool isValid() const { return m_mol >= 0; }
    QString lastError() const { return m_lastError; }

    // Public wrapper around the private activateSession() (indigoSetSessionId(m_session)), for
    // callers outside this class that need to make their own raw indigo* calls against this
    // instance's session -- e.g. DocumentState::importReaction, which iterates reaction
    // components via indigoNext/indigoIterateReactants directly rather than through a dedicated
    // EditableMolecule method. Every EditableMolecule method already activates its own session
    // internally before touching m_mol; this exists only for that one kind of external caller.
    void ensureSession() const { activateSession(); }

    int atomCount() const;
    int bondCount() const;
    StringResult toMolfile() const;

    // Real KET-format JSON export via Indigo's own indigoJson() (already proven correct via its
    // use in IndigoService.cpp) -- C++ equivalent of chem-core.js's getClipboardAsKet. Returns an
    // empty string on any failure (invalid molecule or null Indigo result), matching the real
    // JS's own bare `ketStr = ""` failure mode exactly -- deliberately not a StringResult, since
    // no distinct error is ever surfaced to the real caller either.
    QString toKetJson() const;

    // Extracts an explicit atom+bond selection as standalone MOL-format text -- the C++
    // equivalent of chem-core.js's Struct.clone(atomSet, bondSet), used by copySelection.
    // Preserves any sgroup whose member atoms fall (fully or partially) within the selection,
    // with a correctly trimmed member list under partial selection (confirmed by direct probe
    // against the real indigo.dll). Defensively expands the effective atom set to include both
    // endpoints of every bond in bondIds before extraction -- confirmed by direct probe that
    // indigoCreateEdgeSubmolecule fails outright on an inconsistent selection (a bond without
    // both its endpoint atoms present) otherwise; this expansion eliminates that failure case
    // entirely rather than surfacing it to the caller. Error StringResult if any id is unknown
    // or the resulting submolecule has zero atoms.
    StringResult submoleculeMolfile(const QList<AtomId>& atomIds, const QList<BondId>& bondIds) const;

    // Replaces this object's content IN PLACE, in its OWN existing Indigo session (never
    // touches m_session) -- this is what keeps snapshot()/restore() undo-compatible across a
    // document replace. Clears m_atomIdx/m_bondIdx/m_sgroupIdx/m_sgroupExpanded and resets m_ext
    // to a fresh ExtensionData{}, matching the real JS's _applyLoadedStruct exactly (it never
    // carries extension data across a document replace). Deliberately does NOT reset
    // m_nextAtomId/m_nextBondId/m_nextSGroupId -- ids stay globally unique for the object's
    // whole lifetime regardless of how many documents are loaded into it. On a parse failure:
    // sets lastError(), returns false, and leaves the object COMPLETELY untouched.
    // preserveEzIntent: true for a same-document write-back round-trip (e.g. an auto-layout
    // result, matching DocumentState::deserializeMol's centerOnPage=true case) -- this keeps
    // each double bond's already-recorded "was E/Z ever deliberately specified" value instead
    // of re-deriving it from the incoming data's own (possibly layout-fabricated) coordinates.
    // False (the default) for a genuinely new document, where fresh capture is correct.
    bool loadFrom(const QString& molfileOrSmiles, bool preserveEzIntent = false);

    // Returns a molfile with double bonds whose stereo was never specified
    // (m_ezWasSpecified == 0) made collinear to prevent coordinate-based inference
    // from fabricating E/Z labels. Only used for IUPAC naming, never mutates
    // the actual document.
    QString molfileForNaming() const;

    // Atom/bond CRUD with stable never-reused external IDs
    AtomId addAtom(const QString& symbol, double x, double y);
    bool removeAtom(AtomId id);
    BondId addBond(AtomId a, AtomId b, int order);
    bool removeBond(BondId id);
    QString atomSymbol(AtomId id) const;
    int atomCharge(AtomId id) const;
    bool setAtomCharge(AtomId id, int charge);
    int atomIsotope(AtomId id) const;
    bool setAtomIsotope(AtomId id, int isotope);
    int atomRadical(AtomId id) const;
    bool setAtomRadical(AtomId id, int radical);
    int atomExplicitValence(AtomId id) const;
    bool setAtomExplicitValence(AtomId id, int valence);
    bool setAtomLabel(AtomId id, const QString& newSymbol);
    bool setBondOrderValue(BondId id, int order);
    int atomAttachmentOrder(AtomId id) const;
    bool setAtomAttachmentOrder(AtomId id, int order);
    bool setAtomQueryList(AtomId id, const QString& labelsCsv, bool notList);
    bool clearAtomQueryList(AtomId id, const QString& fallbackLabel);
    bool hasAtomQueryList(AtomId id) const;
    QList<int> atomQueryListNumbers(AtomId id) const;
    bool atomQueryListIsNotList(AtomId id) const;
    bool atomPos(AtomId id, double& x, double& y) const;
    bool setAtomPos(AtomId id, double x, double y);
    
    // Bulk coordinate scaling: scale all atom positions by a factor around a fixed point
    // Returns the number of atoms scaled. Used for rescaling when bond length changes.
    // The fixed point (centerX, centerY) is the center of scaling (e.g., drawing centroid).
    int scaleAllCoordinates(double factor, double centerX = 0, double centerY = 0);
    
    // Compute the centroid of all atoms (for use as scaling center)
    // Returns (0,0) if no atoms exist
    QPointF computeCentroid() const;
    
    int bondOrder(BondId id) const;
    bool bondEndpoints(BondId id, AtomId& a, AtomId& b) const;
    BondId findBond(AtomId a, AtomId b) const;   // -1 if no such bond
    QList<AtomId> atomIds() const;
    QList<BondId> bondIds() const;

    // Extension data storage accessors
    void setName(const QString& n);
    QString name() const;
    int addTextAnnotation(double x, double y, const QString& content);
    int textAnnotationCount() const;
    int addRxnArrow(double x1, double y1, double x2, double y2);
    int rxnArrowCount() const;
    int addRxnPlus(double x, double y);
    int rxnPlusCount() const;
    int addMultitailArrow(const QList<double>& headAndTails);
    int multitailArrowCount() const;
    int addImage(double x, double y, double w, double h, const QByteArray& pngData);
    int imageCount() const;
    bool removeTextAnnotation(TextId id);
    bool removeRxnArrow(RxnArrowId id);
    bool removeRxnPlus(RxnPlusId id);
    bool removeMultitailArrow(MultitailArrowId id);
    bool removeImage(ImageId id);
    QList<RxnArrowId> rxnArrowIds() const;
    QList<RxnPlusId> rxnPlusIds() const;
    QList<MultitailArrowId> multitailArrowIds() const;
    QList<TextId> textAnnotationIds() const;
    QList<ImageId> imageIds() const;
    bool textAnnotationContent(TextId id, double& x, double& y, QString& content, bool& bold, bool& italic) const;
    bool setTextAnnotation(TextId id, const QString& content, bool bold, bool italic);
    bool imageData(ImageId id, double& x, double& y, double& w, double& h, QByteArray& pngData) const;
    bool rxnArrowEndpoints(RxnArrowId id, double& x1, double& y1, double& x2, double& y2) const;
    bool setRxnArrowEndpoints(RxnArrowId id, double x1, double y1, double x2, double y2);
    bool rxnPlusPos(RxnPlusId id, double& x, double& y) const;
    bool setRxnPlusPos(RxnPlusId id, double x, double y);
    QList<double> multitailArrowPoints(MultitailArrowId id) const;
    bool setMultitailArrowPoints(MultitailArrowId id, const QList<double>& points);
    bool imageRect(ImageId id, double& x, double& y, double& w, double& h) const;
    bool setImageRect(ImageId id, double x, double y, double w, double h);
    void setStereoFlag(int fragmentIndex, int flag);
    int stereoFlag(int fragmentIndex) const;
    void setAtomAAM(AtomId id, int mapNumber);
    int atomAAM(AtomId id) const;
    void setAtomCheckWarning(AtomId id, bool warn);
    bool atomCheckWarning(AtomId id) const;

    // Document-level stereo-display toggle (_struct.stereoFlags = {type, groupId}), genuinely
    // consumed by MainToolbar.qml's ABS/REL combo box. Deliberately SEPARATE from the
    // per-fragment stereoFlag/setStereoFlag above, which models an unrelated, real chem-core
    // concept (a per-fragment ABS/AND/OR/MIXED value that is actually COMPUTED from atom stereo
    // labels, never stored) -- reusing that field's QHash<int,int> shape here would be wrong: it
    // cannot represent one document-wide {type, groupId} pair.
    QString stereoFlagsType() const;
    int stereoFlagsGroupId() const;
    // Mutator (sub-project 3e: setStereoFlags, src/worker/50-reactions.js). Named
    // setStereoFlagsDocument, NOT setStereoFlags, to avoid confusion with the unrelated
    // per-fragment setStereoFlag(fragmentIndex, flag) above -- DocumentState::setStereoFlags is
    // the JS-name-matching public entry point that calls this.
    void setStereoFlagsDocument(const QString& type, int groupId);

    // Real string check-warning TYPE (e.g. "valence", set by external structure-check results --
    // src/worker/60-analysis.js's setCheckIssues), NOT the same as the boolean
    // atomCheckWarning/setAtomCheckWarning above (sub-project 1 mismodeled the real field as a
    // yes/no flag). Bonds have no check-warning storage at all before this. Read-only for the
    // same reason as stereoFlagsType.
    QString atomCheckWarningText(AtomId id) const;
    QString bondCheckWarningText(BondId id) const;

    // External CIP-perception / structure-check RESULT APPLICATION (60-analysis.js's
    // setStereoDescriptors/setCheckIssues) -- distinct from the live atomCipDescriptor/
    // stereocenterType/stereocenterGroup accessors above. No-op if id is unknown; getters
    // default to QString()/0 for an unknown id or one with no entry yet.
    void setAtomStereoDescriptor(AtomId id, const QString& cipLabel, int type, int group);
    QString atomStoredCipLabel(AtomId id) const;
    int atomStoredStereoType(AtomId id) const;
    int atomStoredStereoGroup(AtomId id) const;
    void setBondStereoCipLabel(BondId id, const QString& cipLabel);
    QString bondStoredCipLabel(BondId id) const;
    void setAtomCheckWarningText(AtomId id, const QString& text);
    void setBondCheckWarningText(BondId id, const QString& text);

    // Reaction-arrow display fields buildRenderPrimitives reads: display mode, above/below
    // condition text, and optional curvature control point.
    QString rxnArrowMode(RxnArrowId id) const;
    QString rxnArrowConditionsAbove(RxnArrowId id) const;
    QString rxnArrowConditionsBelow(RxnArrowId id) const;
    bool rxnArrowCurvature(RxnArrowId id, double& x, double& y) const;

    // Mutators for the fields above (sub-project 3e: setRxnArrowMode/setRxnArrowConditions from
    // src/worker/50-reactions.js). false if id unknown; setRxnArrowCurvature(has=false) clears
    // curvature back to none.
    bool setRxnArrowMode(RxnArrowId id, const QString& mode);
    bool setRxnArrowConditions(RxnArrowId id, const QString& above, const QString& below);
    bool setRxnArrowCurvature(RxnArrowId id, double x, double y, bool has = true);

    // Pure read-only Indigo queries -- no mutation. implicitHydrogenCount wraps
    // indigoCountHydrogens; neighborAtomIds wraps this class's own bond-endpoint
    // bookkeeping (not a raw Indigo neighbor-iteration call, since AtomId is our own stable id,
    // not Indigo's raw index).
    int implicitHydrogenCount(AtomId id) const;
    QList<AtomId> neighborAtomIds(AtomId id) const;

    // One SSSR ring per entry, atom/bond ids in the ring object's own iteration order.
    // Confirmed by direct probe: SSSR ring objects support indigoIterateAtoms/indigoIterateBonds
    // directly, exactly like a real molecule handle -- no atom-pair-to-bond re-derivation needed.
    struct RingMembership {
        QList<AtomId> atoms;
        QList<BondId> bonds;
    };
    QList<RingMembership> ringMembership() const;

    // Live fragment (connected-component) index for an atom -- wraps indigoComponentIndex.
    // Confirmed by direct probe against the real indigo.dll: pure read query, 0-based, no
    // indigoCountComponents prerequisite, no staleness (updates immediately on graph edits).
    // -1 if id unknown.
    int atomFragmentIndex(AtomId id) const;
    // Reverse of atomFragmentIndex: every atom whose fragment index matches. Needed by
    // RenderPrimitiveBuilder to resolve an R-group member fragment id into the atom set a
    // renderer draws a bracket around (mirrors 10-state.js's _struct.getFragmentIds(fid)).
    QList<AtomId> atomIdsInFragment(int fragIndex) const;

    // Read-only queries on whatever stereo info the molecule already carries from load/perception
    // -- no mutation. 0 (not a stereocenter/stereobond) is the invalid-id sentinel too, matching
    // Indigo's own "0 if not a stereocenter" convention for these two calls.
    int stereocenterType(AtomId id) const;      // 0, or INDIGO_ABS/OR/AND/EITHER (indigo.h)
    int stereocenterGroup(AtomId id) const;     // AND/OR group number
    int bondStereoDirection(BondId id) const;   // 0, or INDIGO_UP/DOWN/EITHER/CIS/TRANS

    // Semantic wedge/hash/either direction, independent of BOTH numeric encodings this feature
    // touches: Indigo's raw indigoBondStereo() value (INDIGO_UP=5, INDIGO_DOWN=6, INDIGO_EITHER=4,
    // indigo.h:438-441) and the V2000 molfile bond-block stereo flag (Up=1, Down=6, Either=4).
    // These disagree on Up (5 vs 1) and must never be passed across this boundary as a raw int --
    // confirmed by a direct probe run for this sub-project: writing V2000 stereo flag 1 and
    // reading back indigoBondStereo() gave 5 for the same bond.
    enum class Direction { None, Up, Down, Either };

    // Translates bondStereoDirection()'s raw Indigo-side int into the semantic enum above. Safe
    // on any input (maps INDIGO_CIS/INDIGO_TRANS or any other unexpected value to Direction::None
    // rather than asserting) -- this is a pure read with no order==1 guard of its own.
    Direction bondStereoDirectionEnum(BondId id) const;

    // Same semantic direction as bondStereoDirectionEnum(), re-encoded as the V2000 molfile
    // stereo flag (0/1/4/6) instead of the Indigo-raw int bondStereoDirection() returns --
    // what any RENDERER wants (matches what gets written to a saved .mol file and what
    // MoleculeLayer.qml's stereo===1/4/6 checks expect), as opposed to bondStereoDirection()'s
    // raw-Indigo contract, which nothing outside bondStereoDirectionEnum() should consume
    // directly. See the Direction-enum comment above for why these two numberings must never
    // be conflated.
    int bondStereoDirectionV2000(BondId id) const;

    // Sets or clears a SINGLE bond's wedge/hash/either stereo direction. Indigo exposes no live
    // per-bond wedge setter (confirmed absent from the vendored indigo.h). This works by patching
    // the bond's V2000 stereo flag directly in a fresh toMolfile() serialization and reloading
    // through indigoLoadMoleculeFromString, which lets Indigo re-perceive real chirality from the
    // wedge plus the molecule's existing 2D coordinates, exactly like any MDL-format toolkit.
    // Returns false, with no mutation and no history entry expected by callers, if `id` isn't a
    // real bond, if the bond's order isn't exactly 1 (wedge/hash is only meaningful on single
    // bonds -- double-bond cis/trans is a separate, unrelated mechanism, out of scope here), or
    // if this molecule's toMolfile() ever emits V3000 instead of V2000 (unconfirmed whether any
    // molecule this port can construct triggers that -- guarded defensively, not assumed
    // impossible). Does NOT check whether the bond's atom already has a different bond wedged;
    // calling this on a second bond of an already-wedged stereocenter is allowed and will not
    // fail, but is the caller's responsibility to avoid if it matters (a future UI wedge tool
    // would enforce at-most-one-wedge-per-stereocenter itself, matching how Sketcher does it by
    // flipping rather than adding a second wedge -- not attempted at this primitive level).
    bool setBondStereo(BondId id, Direction dir);

    // CIP descriptor (R/S/r/s/E/Z, or NONE/UNKNOWN) as Indigo's own CIPDesc enum value, cast to
    // int (molecule_cip_calculator.h: NONE=0, UNKNOWN=1, s=2, r=3, S=4, R=5, E=6, Z=7). Computed
    // via a throwaway clone since indigoAddCIPStereoDescriptors mutates -- this molecule (m_mol)
    // is never touched. Confirmed by direct probe: no bond-level (E/Z) path exists through this
    // API applied to a double bond's own atoms -- this accessor is atom-only, matching Indigo's
    // own "does not represent an atom" rejection of a bond handle.
    int atomCipDescriptor(AtomId id) const;
    // Same descriptor for every atom from ONE CIP computation. Callers that need more than
    // one atom's value (the render-primitive builder asks for all of them on every refresh)
    // must use this: the per-atom accessor above re-runs the whole-molecule CIP calculation
    // on each call.
    QHash<AtomId, int> allAtomCipDescriptors() const;
    QMap<BondId, QString> allBondCipLabels() const;

    int addDataSGroup(const QList<AtomId>& atoms, const QString& description, const QString& data);
    int dataSGroupCount() const;

    MoleculeSnapshot snapshot() const;
    bool restore(const MoleculeSnapshot&);

    // The "next id to hand out" value of every id family. DocumentState rewinds these
    // before re-running a command on redo, so the command hands out the same ids it did
    // the first time and the history steps recorded after it still address its entities.
    // Not for general use: lowering a counter is only safe when nothing alive holds an id
    // at or above it, which is exactly the state undo leaves behind.
    struct IdCounters {
        AtomId atom; BondId bond; SGroupId sgroup; TextId text; RxnArrowId rxnArrow;
        RxnPlusId rxnPlus; MultitailArrowId multitailArrow; ImageId image;
    };
    IdCounters idCounters() const;
    void setIdCounters(const IdCounters& counters);

    // Re-expresses fuseOverlappingAtoms (10-state.js:179-228) in terms Indigo
    // actually supports: Indigo has no "reassign a bond's endpoints" API, so a
    // coincident-atom merge is done as capture-incident-bonds, try to create
    // each as a replacement bond onto the kept atom, then remove the doomed
    // atom (removeAtom's cascade cleans up whatever of doomed's original
    // bonds remain). See the header's "Atom-merge mechanics" spike paragraph
    // for the three Indigo behaviors this depends on.
    struct MergeResult {
        QHash<AtomId, AtomId> mergedAway;   // doomed id -> kept id
        QList<BondId> createdBonds;         // NEW replacement bonds created while rewiring
    };
    MergeResult mergeOverlappingAtoms(double tolerance = 0.1, const QList<AtomId>& candidateAtoms = {});

    // Rewires `doomed`'s incident bonds onto `kept` (creating each as a new
    // bond, or dropping it as a rejected-seam parallel edge), then removes
    // `doomed`. Returns the newly created replacement bonds. Extracted from
    // mergeOverlappingAtoms's own per-pair body (sub-project 3c) so
    // insertFunctionalGroup's explicit-substitution graft (sub-project 3d)
    // can reuse the identical mechanism instead of duplicating it.
    QList<BondId> graftAtomOnto(AtomId doomed, AtomId kept);

    // Non-undoable document-seeding helper (30-templates.js's addBenzeneRing)
    // -- NOT a DocumentState command; mutates directly and has no history
    // entry, matching the real function exactly.
    QList<AtomId> addBenzeneRing(double cx, double cy);

    // Re-expresses chem-core's tempStruct.mergeInto(_struct, ..., aidMap) in
    // terms of indigoMerge, whose own return value is a plain status code,
    // not an id mapping. Confirmed by direct probe against the real
    // indigo.dll (not assumed): indigoMerge appends the source's atoms/bonds
    // to this molecule in the source's OWN iteration order, and automatically
    // propagates any SUP sgroups (with attachment points) the source carries.
    // Before/after diffing therefore correctly discovers every created id,
    // including sgroups -- no manual sgroup reconstruction is needed.
    //
    // Takes Molfile TEXT, not a raw Indigo handle: TemplateLibrary owns its
    // own Indigo session (separate from this EditableMolecule's), and a raw
    // handle is a small integer scoped to whichever session created it --
    // confirmed by direct probe that reusing one across sessions either fails
    // or silently aliases an unrelated object. TemplateLibrary::molfileText()
    // is the confirmed-safe way to cross that boundary.
    struct InsertResult {
        QList<AtomId> createdAtoms;
        QList<BondId> createdBonds;
        QList<SGroupId> createdSGroups;
        QHash<int, AtomId> sourceIndexToNewAtomId;   // source molecule's internal atom index -> new AtomId
    };
    InsertResult insertStructure(const QString& sourceMolfile,
                                  const std::function<QPointF(double, double)>& transform);

    // Superatom (SUP sgroup) support. A new stable, never-reused SGroupId
    // (like AtomId/BondId) rather than the existing addDataSGroup precedent's
    // raw, reusable indigoIndex() -- these ARE looked up again later
    // (toggleSgroupExpanded), so a reused index could alias a different,
    // newer sgroup. There is no public constructor: insertStructure discovers
    // and assigns ids for merge-propagated sgroups internally, and nothing in
    // this sub-project's scope needs to build one from bare atoms.
    bool superatomAttachAtom(SGroupId id, AtomId& out) const;   // first attachment point's atom, if any
    void setSGroupExpanded(SGroupId id, bool expanded);
    bool sgroupExpanded(SGroupId id) const;
    void setSGroupLabel(SGroupId id, const QString& label);
    QString sgroupLabel(SGroupId id) const;

    // Sgroup introspection beyond superatomAttachAtom's single-attachment-point view: every
    // currently-existing sgroup id, and a given sgroup's FULL member-atom set (not just its
    // attach atom). Needed for correct render-primitive contraction (sub-project 4), which must
    // hide every member atom of a collapsed sgroup, not only its attach atom.
    QList<SGroupId> sgroupIds() const;
    QList<AtomId> sgroupMemberAtomIds(SGroupId id) const;

    // Indigo's current internal 0..count-1 atom/bond order, reverse-mapped to external
    // AtomId/BondId via m_atomIdx/m_bondIdx -- NOT the same as atomIds()/bondIds()'s
    // sorted-by-id order, which does not survive a remove-then-add sequence (confirmed by
    // direct probe: a newly-added atom/bond can reuse a just-freed internal index instead of
    // being appended at the end). Needed to resolve a positional index from an externally
    // reparsed molfile's own numbering (e.g. an analysis-tool JSON payload) back to the
    // correct external id.
    QList<AtomId> atomIdsInIndigoOrder() const;
    QList<BondId> bondIdsInIndigoOrder() const;

    // Removes ONLY the sgroup grouping (indigoRemove on the superatom handle) -- member
    // atoms are left completely untouched. Confirmed by direct probe against the real
    // indigo.dll. Used for "the whole collapsed pill was selected directly" deletion.
    bool removeSuperatomOnly(SGroupId id);

    // Builds a new sgroup from a bare atom list -- needed purely so undo can recreate a
    // deleted sgroup. No `name`/label parameter: confirmed by direct probe that neither
    // indigoDescription nor indigoName recovers a superatom's original display label, so
    // there is nothing to pass through even if this method took one (see
    // RenderPrimitives.cpp:86-89's pre-existing sub-project-4 gap). Assigns a fresh
    // SGroupId via m_nextSGroupId, same counter insertStructure's own sgroup discovery
    // uses. Returns -1 if memberAtoms is empty or any id is unknown.
    SGroupId createSuperatomFromAtoms(const QList<AtomId>& memberAtoms);

    // R-group (Markush) storage: fragIds/range/resth/ifthen per R-group NUMBER. See
    // ExtensionData.h's RGroupEntry for field meanings.
    bool addRGroupEntry(int rgroupNumber);
    bool removeRGroupEntry(int rgroupNumber);
    bool setRGroupLogic(int rgroupNumber, const QString& range, bool resth, int ifthen);
    bool addRGroupFragment(int rgroupNumber, int fragId);
    bool removeRGroupFragment(int rgroupNumber, int fragId);
    QList<int> rgroupNumbers() const;
    bool rgroupLogic(int rgroupNumber, QString& range, bool& resth, int& ifthen) const;
    QList<int> rgroupFragmentIds(int rgroupNumber) const;

    // Bracket-selection annotations (50-reactions.js's addBracketSelection): a plain push/pop
    // stack, no id -- matches the real `_struct.brackets.push(bracket)` / `.pop()` exactly.
    void pushBracket(double minX, double minY, double maxX, double maxY);
    bool popBracket();
    int bracketCount() const;
    bool bracketAt(int index, double& minX, double& minY, double& maxX, double& maxY) const;

private:
    // One Indigo session per instance, held for the object's lifetime
    // (unlike IndigoService's alloc/release-per-call one-shot pattern).
    unsigned long long m_session = 0;
    int m_mol = -1;               // Indigo molecule handle, -1 = invalid
    QString m_lastError;

    QHash<AtomId, int> m_atomIdx; // external stable ID -> indigo atom index
    QHash<int, AtomId> m_idxToAtom; // indigo atom index -> external stable ID (reverse of m_atomIdx)
    QHash<BondId, int> m_bondIdx; // external stable ID -> indigo bond index
    QHash<int, BondId> m_idxToBond; // indigo bond index -> external stable ID (reverse of m_bondIdx)
    AtomId m_nextAtomId = 1;      // monotonic, never reused
    BondId m_nextBondId = 1;

    ExtensionData m_ext;

    // Track E/Z stereo intent at load time: BondId -> indigoBondStereo() value at initial load
    // (0 = unspecified, 7 = CIS, 8 = TRANS, etc.). Only set once at load time, never updated
    // through manual edits. Used to prevent fabricated E/Z labels in IUPAC naming.
    QMap<BondId, int> m_ezWasSpecified;

    QHash<SGroupId, int> m_sgroupIdx;   // external stable ID -> indigo superatom index
    SGroupId m_nextSGroupId = 1;        // monotonic, never reused
    QHash<int, bool> m_sgroupExpanded;  // SGroupId -> expanded flag (Indigo has no such concept)
    QHash<int, QString> m_sgroupLabels; // SGroupId -> display label (Indigo has no such concept either)

    void activateSession() const; // indigoSetSessionId(m_session)
    void rebuildIndexTables();
    // Records each double bond's E/Z intent. For a bond whose Indigo-internal index has an
    // entry in preservedIndexToStereo (see loadFrom), that value is kept instead of a fresh
    // indigoBondStereo() read; pass {} for a genuinely new load.
    void captureEzIntent(const QMap<int, int>& preservedIndexToStereo);

    // Same "Indigo's current 0..count-1 order, reverse-mapped to external SGroupId via
    // m_sgroupIdx" contract as atomIdsInIndigoOrder()/bondIdsInIndigoOrder() above, but for
    // sgroups. Private (unlike its atom/bond siblings, which are public for an existing
    // documented external use case) because nothing outside setBondStereo needs it yet. Needed
    // to remap m_sgroupIdx after setBondStereo replaces the WHOLE Indigo molecule handle --
    // rebuildIndexTables()'s own sgroup handling only PRUNES vanished sgroups (its own comment,
    // EditableMolecule.cpp:626-631, is explicitly about in-place removal on the SAME handle, not
    // about a full handle replacement) and does not reindex survivors, which setBondStereo needs.
    QList<SGroupId> sgroupIdsInIndigoOrder() const;

    // The two numeric-encoding conversions Direction exists to keep apart. Both private: nothing
    // outside setBondStereo/bondStereoDirectionEnum needs raw access to either numbering.
    static int directionToV2000Code(Direction dir);
    static Direction directionFromIndigoBondStereo(int indigoValue);
    void assignFreshAtomBondIds();   // extracted from the constructor's own loop, shared with loadFrom
};

#endif // EDITABLEMOLECULE_H
