// tests/document_state_test.cpp
// Standalone tests for DocumentState (chem-core.js migration, sub-project 2).
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "app/molecule/SelectionState.h"
#include "app/molecule/EditCommand.h"
#include "app/molecule/DocumentState.h"
#include "app/molecule/TemplateLibrary.h"
#include "indigo.h"
#include <QFile>
#include <QDir>

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

static QString writeTempSdf(const QString& content, const QString& baseName) {
    QString path = QDir::tempPath() + QStringLiteral("/") + baseName + QStringLiteral(".sdf");
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        std::printf("[ERROR] writeTempSdf: failed to open '%s' for writing (error: %s)\n",
                    path.toUtf8().constData(), f.errorString().toUtf8().constData());
        return path;
    }
    f.write(content.toUtf8());
    f.close();
    return path;
}

static void test_selectionStateBasics() {
    std::printf("--- Test 0: SelectionState basics ---\n");
    SelectionState s;
    CHECK(s.isEmpty(), "fresh SelectionState is empty");

    s.atoms.insert(1);
    s.atoms.insert(2);
    s.bonds.insert(5);
    s.rxnArrows.insert(9);
    s.rxnPluses.insert(3);
    s.multitailArrows.insert(7);
    CHECK(!s.isEmpty(), "non-empty after inserts");
    CHECK(s.atoms.size() == 2, "2 atoms in the set");
    CHECK(s.atoms.contains(1) && s.atoms.contains(2), "both atom ids present");

    s.clear();
    CHECK(s.isEmpty(), "empty after clear()");
    CHECK(s.atoms.isEmpty() && s.bonds.isEmpty() && s.rxnArrows.isEmpty()
          && s.rxnPluses.isEmpty() && s.multitailArrows.isEmpty(),
          "every individual set is empty after clear()");
}

static void test_editCommandBasics() {
    std::printf("--- Test 1: EditCommand basics ---\n");
    int counter = 0;
    EditCommand cmd;
    cmd.execute = [&counter]() { counter += 5; };
    cmd.invert = [&counter]() { counter -= 5; };

    cmd.execute();
    CHECK(counter == 5, "execute() runs the forward closure");
    cmd.invert();
    CHECK(counter == 0, "invert() runs the inverse closure");
}

static void test_documentStateUndoRedo() {
    std::printf("--- Test 2: DocumentState construction, undo/redo, dirty ---\n");
    DocumentState doc;
    CHECK(!doc.isDirty(), "fresh DocumentState is not dirty");
    CHECK(!doc.canUndo(), "fresh DocumentState cannot undo");
    CHECK(!doc.canRedo(), "fresh DocumentState cannot redo");

    EditableMolecule& mol = doc.molecule();
    AtomId addedId = -1;
    EditCommand addCarbon;
    addCarbon.execute = [&mol, &addedId]() {
        addedId = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    };
    addCarbon.invert = [&mol, &addedId]() {
        mol.removeAtom(addedId);
        addedId = -1;
    };

    doc.executeCommand(addCarbon);
    CHECK(doc.molecule().atomCount() == 1, "atom added by the command");
    CHECK(doc.isDirty(), "dirty after executeCommand");
    CHECK(doc.canUndo(), "canUndo true after one command");
    CHECK(!doc.canRedo(), "canRedo false with nothing undone yet");

    doc.undo();
    CHECK(doc.molecule().atomCount() == 0, "atom removed by undo (invert ran)");
    CHECK(!doc.canUndo(), "canUndo false after undoing the only command");
    CHECK(doc.canRedo(), "canRedo true after an undo");

    doc.redo();
    CHECK(doc.molecule().atomCount() == 1, "atom re-added by redo (execute reran)");
    CHECK(doc.canUndo(), "canUndo true again after redo");
    CHECK(!doc.canRedo(), "canRedo false again after redoing everything");

    doc.markClean();
    CHECK(!doc.isDirty(), "markClean() clears the dirty flag");

    EditCommand addOxygen;
    addOxygen.execute = [&mol]() { mol.addAtom(QStringLiteral("O"), 1.0, 0.0); };
    addOxygen.invert = [&mol]() { /* not exercised in this check */ };
    doc.executeCommand(addOxygen);
    CHECK(doc.isDirty(), "dirty again after a command following markClean");

    doc.undo();
    CHECK(doc.canRedo(), "redo available right after an undo");
    EditCommand addNitrogen;
    addNitrogen.execute = [&mol]() { mol.addAtom(QStringLiteral("N"), 2.0, 0.0); };
    addNitrogen.invert = [&mol]() { /* not exercised in this check */ };
    doc.executeCommand(addNitrogen);
    CHECK(!doc.canRedo(), "old redo branch truncated by executing a new command");
}

static void test_executeCommandReentrancyGuard() {
    std::printf("--- Test: executeCommand ignores a reentrant call from within its own execute ---\n");
    DocumentState doc;
    size_t historyBefore = 0; // tracked via canUndo()/redo() round trips below instead of a direct accessor

    int innerRanCount = 0;
    EditCommand inner;
    inner.execute = [&innerRanCount]() { ++innerRanCount; };
    inner.invert = [&innerRanCount]() { --innerRanCount; };

    int outerRanCount = 0;
    EditCommand outer;
    outer.execute = [&doc, &outerRanCount, inner]() mutable {
        ++outerRanCount;
        doc.executeCommand(inner); // reentrant -- must be a silent no-op
    };
    outer.invert = [&outerRanCount]() { --outerRanCount; };

    doc.executeCommand(outer);

    CHECK(outerRanCount == 1, "outer command's own execute ran exactly once");
    CHECK(innerRanCount == 0, "inner (reentrant) command's execute never ran -- guard blocked it before cmd.execute()");
    CHECK(doc.canUndo(), "outer command was still recorded in history");

    doc.undo();
    CHECK(outerRanCount == 0, "undoing the outer command runs its own invert normally");
    CHECK(!doc.canUndo(), "history has exactly one entry -- the inner reentrant call was never pushed");
    (void)historyBefore;
}

static void test_undoReentrancyGuard() {
    std::printf("--- Test: undo() ignores a reentrant call from within a command's own invert ---\n");
    DocumentState doc;
    AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
    (void)a;
    CHECK(doc.canUndo(), "setup: one command in history");

    int reentrantInvertCalls = 0;
    EditCommand reentrantCmd;
    reentrantCmd.execute = []() {};
    reentrantCmd.invert = [&doc, &reentrantInvertCalls]() {
        ++reentrantInvertCalls;
        doc.undo(); // reentrant -- must be a silent no-op
    };
    doc.executeCommand(reentrantCmd);
    CHECK(doc.canUndo(), "setup: two commands in history now");

    doc.undo(); // undoes reentrantCmd; its own invert tries to call undo() again
    CHECK(reentrantInvertCalls == 1, "reentrantCmd's invert ran exactly once");
    CHECK(doc.canUndo(), "history pointer moved back by exactly ONE step, not two -- "
                          "the reentrant undo() call inside invert() was a no-op");
}

static void test_selectionReconciliationRemovesFullyStaleId() {
    std::printf("--- Test: stale sgroup id (no numeric collision) is removed from selection ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 1, 0);
    AtomId a3 = doc.addAtom(QStringLiteral("C"), 2, 0);
    // Whether this sgroup's id numerically collides with a1/a2/a3 doesn't matter for this test:
    // once ALL of its members (a1, a2, a3) are deleted, no atom exists with that id either way,
    // so the assertion below (the id is gone from selection) holds regardless.
    SGroupId sid = mol.createSuperatomFromAtoms({a1, a2, a3});
    CHECK(sid != -1, "setup: sgroup created");

    doc.selection().atoms.insert(sid);
    CHECK(doc.selection().atoms.contains(sid), "setup: sgroup id is selected");

    doc.deleteAtom(a1);
    doc.deleteAtom(a2);
    doc.deleteAtom(a3); // last member removed -- Indigo auto-destroys the sgroup here

    CHECK(!doc.selection().atoms.contains(sid),
          "the now-fully-stale sgroup id was removed from selection after the command that "
          "destroyed it, even though deleteAtom itself never touches selection directly");
}

static void test_selectionReconciliationKeepsCollidingRealAtom() {
    std::printf("--- Test: a real atom's id survives even if a destroyed sgroup once shared it ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0); // AtomId 1
    AtomId b = doc.addAtom(QStringLiteral("C"), 1, 0);
    AtomId c = doc.addAtom(QStringLiteral("C"), 2, 0);
    AtomId d = doc.addAtom(QStringLiteral("C"), 3, 0);
    SGroupId sid = mol.createSuperatomFromAtoms({b, c, d}); // first sgroup -> SGroupId 1, colliding with AtomId 1 (atom `a`)
    CHECK(sid == a, "setup: the sgroup's id naturally collides with atom a's id (both are 1)");

    doc.selection().atoms.insert(a); // representing "atom a is selected"; a's id happens to equal sid

    doc.deleteAtom(b);
    doc.deleteAtom(c);
    doc.deleteAtom(d); // destroys the sgroup as a side effect

    CHECK(doc.selection().atoms.contains(a),
          "atom a's id survives reconciliation -- it's a real, currently-valid atom, "
          "unaffected by the unrelated sgroup that happened to share its number");
    CHECK(mol.atomIds().contains(a) && mol.atomSymbol(a) == QStringLiteral("C"),
          "atom a itself was never touched by deleting the sgroup's own members");
}

static void test_selectionReconciliationCleansUpAfterDeleteSelectionEntities() {
    std::printf("--- Test: deleteSelectionEntities's own selected ids are cleaned up afterward ---\n");
    DocumentState doc;
    AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
    doc.selectAtom(a);
    CHECK(!doc.selection().atoms.isEmpty(), "setup: atom is selected");

    doc.deleteSelectionEntities();

    CHECK(doc.selection().atoms.isEmpty(),
          "selection is empty after deleting the very atom it referenced -- a bonus side effect "
          "of the same reconciliation, not just the cross-command aliasing scenario above");
}

static void test_selectionReconciliationRemovesStaleBondId() {
    std::printf("--- Test: stale bond id is removed from selection after a direct deleteBond ---\n");
    DocumentState doc;
    AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b = doc.addAtom(QStringLiteral("C"), 1, 0);
    BondId bond = doc.addBond(a, b, 1);
    doc.selection().bonds.insert(bond);
    CHECK(doc.selection().bonds.contains(bond), "setup: bond is selected");

    doc.deleteBond(bond);

    CHECK(!doc.selection().bonds.contains(bond),
          "the now-stale bond id was removed from selection after the command that deleted "
          "it, even though deleteBond itself never touches selection directly -- this is "
          "the exact gap reconcileSelectionAfterCommand had for every entity type except atoms");
}

static void test_selectionReconciliationRemovesStaleRxnArrowId() {
    std::printf("--- Test: stale rxnArrow id is removed from selection after a direct deleteRxnArrow ---\n");
    DocumentState doc;
    RxnArrowId arrow = doc.addRxnArrow(0, 0);
    doc.selection().rxnArrows.insert(arrow);
    CHECK(doc.selection().rxnArrows.contains(arrow), "setup: rxnArrow is selected");

    doc.deleteRxnArrow(arrow);

    CHECK(!doc.selection().rxnArrows.contains(arrow),
          "the now-stale rxnArrow id was removed from selection after deletion");
}

static void test_selectionReconciliationRemovesStaleRxnPlusId() {
    std::printf("--- Test: stale rxnPlus id is removed from selection after a direct deleteRxnPlus ---\n");
    DocumentState doc;
    RxnPlusId plus = doc.addRxnPlus(0, 0);
    doc.selection().rxnPluses.insert(plus);
    CHECK(doc.selection().rxnPluses.contains(plus), "setup: rxnPlus is selected");

    doc.deleteRxnPlus(plus);

    CHECK(!doc.selection().rxnPluses.contains(plus),
          "the now-stale rxnPlus id was removed from selection after deletion");
}

static void test_selectionReconciliationRemovesStaleMultitailArrowId() {
    std::printf("--- Test: stale multitailArrow id is removed from selection after a direct deleteMultitailArrow ---\n");
    DocumentState doc;
    MultitailArrowId mta = doc.addMultitailArrow(0, 0);
    doc.selection().multitailArrows.insert(mta);
    CHECK(doc.selection().multitailArrows.contains(mta), "setup: multitailArrow is selected");

    doc.deleteMultitailArrow(mta);

    CHECK(!doc.selection().multitailArrows.contains(mta),
          "the now-stale multitailArrow id was removed from selection after deletion");
}

static void test_documentStateHistoryCap() {
    std::printf("--- Test 3: DocumentState history cap at 50 ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    for (int i = 0; i < 51; ++i) {
        EditCommand cmd;
        cmd.execute = [&mol]() { mol.addAtom(QStringLiteral("C"), 0.0, 0.0); };
        cmd.invert = [&mol]() {
            QList<AtomId> ids = mol.atomIds();
            if (!ids.isEmpty()) mol.removeAtom(ids.last());
        };
        doc.executeCommand(cmd);
    }
    CHECK(doc.molecule().atomCount() == 51, "51 atoms added");
    int undoCount = 0;
    while (doc.canUndo() && undoCount <= 51) {
        doc.undo();
        ++undoCount;
    }
    CHECK(undoCount == 50, "at most 50 undos possible after 51 commands (history capped)");
}

static void test_documentStateEditingOperations() {
    std::printf("--- Test 5: DocumentState editing operations ---\n");
    DocumentState doc;

    // addAtom / addBond / deleteBond / deleteAtom, with undo/redo.
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("O"), 1, 0);
    CHECK(doc.molecule().atomCount() == 2, "2 atoms added via DocumentState::addAtom");
    BondId b1 = doc.addBond(a1, a2, 1);
    CHECK(doc.molecule().bondCount() == 1, "bond added via DocumentState::addBond");

    doc.deleteBond(b1);
    CHECK(doc.molecule().bondCount() == 0, "deleteBond removes the bond");
    doc.undo();
    CHECK(doc.molecule().bondCount() == 1, "undo restores the deleted bond");

    doc.deleteAtom(a2);
    CHECK(doc.molecule().atomCount() == 1, "deleteAtom removes the atom");
    CHECK(doc.molecule().bondCount() == 0, "deleteAtom also removes its incident bond");
    doc.undo();
    CHECK(doc.molecule().atomCount() == 2, "undo restores the deleted atom");
    CHECK(doc.molecule().bondCount() == 1, "undo restores the incident bond too");

    // deleteAtom's undo recreates the atom under a NEW id (documented
    // limitation -- EditableMolecule has no restore-at-exact-id primitive),
    // so `a2`'s original value is now stale. Re-derive the current id of
    // the recreated oxygen atom for the rest of this test.
    for (AtomId id : doc.molecule().atomIds()) {
        if (id != a1) a2 = id;
    }

    // changeAtomLabel: guarded (no-op when unchanged).
    doc.changeAtomLabel(a1, QStringLiteral("N"));
    CHECK(doc.molecule().atomSymbol(a1) == QStringLiteral("N"), "changeAtomLabel changes the symbol");
    bool canUndoBeforeNoOpLabel = doc.canUndo();
    doc.changeAtomLabel(a1, QStringLiteral("N")); // same label: must no-op, no new history entry
    CHECK(doc.canUndo() == canUndoBeforeNoOpLabel, "changeAtomLabel to the SAME label pushes no history entry");
    doc.undo();
    CHECK(doc.molecule().atomSymbol(a1) == QStringLiteral("C"), "undo restores the original label");

    // changeAtomCharge: guarded.
    doc.changeAtomCharge(a1, 1);
    CHECK(doc.molecule().atomCharge(a1) == 1, "changeAtomCharge sets the charge");
    doc.undo();
    CHECK(doc.molecule().atomCharge(a1) == 0, "undo restores the original charge (0)");

    // changeAtomIsotope: NOT guarded (real JS always executes, even for an
    // unchanged value). Proof: call it twice with the SAME value, then undo
    // twice -- if it were guarded, the second call would push no history
    // entry and a single undo would already reach isotope 0.
    doc.changeAtomIsotope(a1, 5);
    CHECK(doc.molecule().atomIsotope(a1) == 5, "changeAtomIsotope sets isotope to 5");
    doc.changeAtomIsotope(a1, 5); // same value again -- must still push a history entry
    doc.undo();
    CHECK(doc.molecule().atomIsotope(a1) == 5,
          "first undo only undoes the second (no-op-value) call -- isotope still 5, proving it WAS pushed");
    doc.undo();
    CHECK(doc.molecule().atomIsotope(a1) == 0, "second undo undoes the real change, back to 0");

    // changeAtomRadical / changeAtomValence: real end-to-end undo of an actual change.
    doc.changeAtomRadical(a1, INDIGO_SINGLET);
    CHECK(doc.molecule().atomRadical(a1) == INDIGO_SINGLET, "changeAtomRadical sets the radical");
    doc.undo();
    CHECK(doc.molecule().atomRadical(a1) == 0, "undo restores radical to 0");

    doc.changeAtomValence(a1, 4);
    CHECK(doc.molecule().atomExplicitValence(a1) == 4, "changeAtomValence sets explicit valence");
    doc.undo();
    CHECK(doc.molecule().atomExplicitValence(a1) == -1, "undo restores valence to unset (-1)");

    // setAttachmentPoint: NOT guarded.
    doc.setAttachmentPoint(a1, 1);
    CHECK(doc.molecule().atomAttachmentOrder(a1) == 1, "setAttachmentPoint sets order 1");
    doc.undo();
    CHECK(doc.molecule().atomAttachmentOrder(a1) == 0, "undo restores attachment order to 0 (none)");

    // setAtomMapping: guarded.
    doc.setAtomMapping(a1, 5);
    CHECK(doc.molecule().atomAAM(a1) == 5, "setAtomMapping sets the AAM number");
    doc.undo();
    CHECK(doc.molecule().atomAAM(a1) == 0, "undo restores AAM to 0");

    // atomProperties: read-only, not a command.
    doc.changeAtomCharge(a1, 2);
    DocumentState::AtomProperties props = doc.atomProperties(a1);
    CHECK(props.label == QStringLiteral("C"), "atomProperties reports the current label");
    CHECK(props.charge == 2, "atomProperties reports the current charge");

    // changeBondOrder: guarded. Reuses the bond deleteAtom's undo already
    // restored between a1 and (the recreated) a2 -- adding another bond on
    // top of it would fail, since Indigo rejects a duplicate parallel edge
    // between the same atom pair (probe-confirmed: "already have edge
    // between vertices").
    BondId liveBond = doc.molecule().bondIds().first();
    doc.changeBondOrder(liveBond, 3);
    CHECK(doc.molecule().bondOrder(liveBond) == 3, "changeBondOrder sets the new order");
    doc.undo();
    CHECK(doc.molecule().bondOrder(liveBond) == 1, "undo restores the original order");

    // setBondStereo: guarded (unchanged value no-ops), one EditCommand, undo/redo round-trips
    // through the semantic Direction enum, never a raw int (sub-project 8's own found-and-fixed
    // encoding-mismatch bug: Direction::Up is 1 on the V2000 write side, 5 (INDIGO_UP) on the
    // Indigo read side -- a naive undo that skipped the translation would write an invalid V2000
    // stereo code).
    //
    // NOTE (found by direct testing, not anticipated in the approved spec): liveBond is the
    // single bond of a bare 2-atom molecule (built earlier in this test for changeBondOrder),
    // whose begin-atom has only one real substituent -- not a chemically valid wedge origin.
    // Indigo correctly rejects the reload ("direction of bond #0 makes no sense", the same real
    // validation found during Task 1's own execution), so setBondStereo silently no-ops. This
    // test instead builds its own dedicated 4-substituent stereocenter, the same proven shape
    // used throughout Task 1's tests.
    AtomId stereoCenter = doc.addAtom(QStringLiteral("C"), 10, 10);
    AtomId stF = doc.addAtom(QStringLiteral("F"), 10, 11);
    AtomId stCl = doc.addAtom(QStringLiteral("Cl"), 9, 9.5);
    AtomId stBr = doc.addAtom(QStringLiteral("Br"), 11, 9.5);
    doc.addBond(stereoCenter, stF, 1);
    doc.addBond(stereoCenter, stCl, 1);
    BondId stereoBond = doc.addBond(stereoCenter, stBr, 1);

    CHECK(doc.molecule().bondStereoDirectionEnum(stereoBond) == EditableMolecule::Direction::None,
          "stereoBond starts with no stereo");

    doc.setBondStereo(stereoBond, EditableMolecule::Direction::Up);
    CHECK(doc.molecule().bondStereoDirectionEnum(stereoBond) == EditableMolecule::Direction::Up,
          "setBondStereo sets Up");
    CHECK(doc.molecule().bondStereoDirection(stereoBond) == INDIGO_UP,
          "raw indigoBondStereo reports INDIGO_UP specifically");

    doc.undo();
    CHECK(doc.molecule().bondStereoDirectionEnum(stereoBond) == EditableMolecule::Direction::None,
          "undo restores None");

    doc.redo();
    CHECK(doc.molecule().bondStereoDirection(stereoBond) == INDIGO_UP,
          "redo restores INDIGO_UP exactly");

    // REQUIRED: the encoding-mismatch regression check. Up -> None -> undo must restore the
    // EXACT raw value INDIGO_UP (5), not merely a truthy/nonzero value and not the raw V2000
    // write code (1) -- a test that only checked truthiness would NOT catch this bug class.
    doc.setBondStereo(stereoBond, EditableMolecule::Direction::None);
    CHECK(doc.molecule().bondStereoDirectionEnum(stereoBond) == EditableMolecule::Direction::None,
          "setBondStereo back to None");
    doc.undo();
    CHECK(doc.molecule().bondStereoDirection(stereoBond) == INDIGO_UP,
          "undo of Up->None restores EXACTLY INDIGO_UP (5), not merely truthy/nonzero");
    StringResult stereoMf = doc.molecule().toMolfile();
    CHECK(stereoMf.success, "molfile serialization still succeeds after the undo round-trip");

    // Guarded: setting to the SAME current direction (Up, from the undo above) pushes no
    // history entry.
    doc.setBondStereo(stereoBond, EditableMolecule::Direction::Up); // same value again -- no-op
    doc.undo();
    CHECK(doc.molecule().bondStereoDirectionEnum(stereoBond) == EditableMolecule::Direction::None,
          "the guarded same-value call pushed nothing -- this undo is the earlier Up->None entry");
    doc.redo();

    // Attempt on a double bond -> no history entry, no mutation.
    AtomId dblA1 = doc.addAtom(QStringLiteral("C"), 5, 5);
    AtomId dblA2 = doc.addAtom(QStringLiteral("O"), 5, 6);
    BondId dblBond = doc.addBond(dblA1, dblA2, 2);
    bool canUndoBeforeDoubleAttempt = doc.canUndo();
    doc.setBondStereo(dblBond, EditableMolecule::Direction::Up);
    CHECK(doc.canUndo() == canUndoBeforeDoubleAttempt,
          "setBondStereo on a double bond pushes no history entry");
    CHECK(doc.molecule().bondStereoDirectionEnum(dblBond) == EditableMolecule::Direction::None,
          "double bond's stereo remains None -- the attempt was a true no-op");

    // addBondAndAtom: drag from an existing atom to a new point.
    AtomId abStart = doc.addAtom(QStringLiteral("C"), 20, 20);
    doc.addBondAndAtom(abStart, QStringLiteral("N"), 21.5, 20, 1, 0);
    QList<BondId> bondsAfterAb = doc.molecule().bondIds();
    BondId newBondAb = -1;
    AtomId newAtomAb = -1;
    for (BondId bid : bondsAfterAb) {
        AtomId ba = -1, bb = -1;
        if (doc.molecule().bondEndpoints(bid, ba, bb) && (ba == abStart || bb == abStart)) {
            newBondAb = bid;
            newAtomAb = (ba == abStart) ? bb : ba;
        }
    }
    CHECK(newBondAb != -1 && newAtomAb != -1, "addBondAndAtom created a new bonded atom");
    CHECK(doc.molecule().atomSymbol(newAtomAb) == QStringLiteral("N"), "the new atom has the requested label");
    double abX = 0, abY = 0;
    doc.molecule().atomPos(newAtomAb, abX, abY);
    CHECK(abX == 21.5 && abY == 20, "the new atom is placed at the requested (clamped) coordinates");
    doc.undo();
    CHECK(!doc.molecule().atomIds().contains(newAtomAb), "undo removes the created atom");
    CHECK(doc.molecule().findBond(abStart, newAtomAb) == -1, "undo removes the created bond");
    doc.redo();
    bool foundRecreatedN = false;
    for (BondId bid : doc.molecule().bondIds()) {
        AtomId ba = -1, bb = -1;
        if (doc.molecule().bondEndpoints(bid, ba, bb) && (ba == abStart || bb == abStart)) {
            AtomId other = (ba == abStart) ? bb : ba;
            if (doc.molecule().atomSymbol(other) == QStringLiteral("N")) { foundRecreatedN = true; break; }
        }
    }
    CHECK(foundRecreatedN, "redo re-creates the bonded N atom (with a fresh id, per this port's own established convention)");

    // addBondAndAtom: collision-merge case -- dragging onto an EXISTING atom creates only a bond.
    AtomId mergeStart = doc.addAtom(QStringLiteral("C"), 25, 25);
    AtomId mergeTarget = doc.addAtom(QStringLiteral("O"), 26, 25);
    int bondCountBeforeMerge = doc.molecule().bondIds().size();
    doc.addBondAndAtom(mergeStart, QStringLiteral("C"), 26, 25, 1, 0); // lands within 0.3 of mergeTarget
    CHECK(doc.molecule().bondIds().size() == bondCountBeforeMerge + 1,
          "collision case adds exactly one new bond, not a new atom");
    CHECK(doc.molecule().findBond(mergeStart, mergeTarget) != -1,
          "the new bond connects to the EXISTING colliding atom");

    // addBondAndAtom: page-bounds clamping.
    AtomId clampStart = doc.addAtom(QStringLiteral("C"), 29, 20);
    doc.addBondAndAtom(clampStart, QStringLiteral("C"), 999, 999, 1, 0);
    QList<AtomId> allAtomsAfterClamp = doc.molecule().atomIds();
    AtomId clampedAtom = -1;
    for (AtomId aid : allAtomsAfterClamp) {
        double px = 0, py = 0;
        doc.molecule().atomPos(aid, px, py);
        if (px == 30.0) { clampedAtom = aid; break; } // kPageMaxX
    }
    CHECK(clampedAtom != -1, "an out-of-bounds endpoint is clamped to the page's max X (30.0)");

    // addBondAndAtom: invalid startId is a true no-op.
    bool canUndoBeforeInvalidStart = doc.canUndo();
    int atomCountBeforeInvalidStart = doc.molecule().atomIds().size();
    doc.addBondAndAtom(999999, QStringLiteral("C"), 0, 0, 1, 0);
    CHECK(doc.canUndo() == canUndoBeforeInvalidStart, "invalid startId pushes no history entry");
    CHECK(doc.molecule().atomIds().size() == atomCountBeforeInvalidStart, "invalid startId creates no atom");

    // addBondBetweenCoords: always creates two new atoms + one bond, as one undoable command.
    int atomCountBeforeBc = doc.molecule().atomIds().size();
    int bondCountBeforeBc = doc.molecule().bondIds().size();
    doc.addBondBetweenCoords(40, 40, 41.5, 40, 1, 0);
    CHECK(doc.molecule().atomIds().size() == atomCountBeforeBc + 2, "addBondBetweenCoords creates exactly 2 atoms");
    CHECK(doc.molecule().bondIds().size() == bondCountBeforeBc + 1, "addBondBetweenCoords creates exactly 1 bond");
    doc.undo();
    CHECK(doc.molecule().atomIds().size() == atomCountBeforeBc, "undo removes both created atoms");
    CHECK(doc.molecule().bondIds().size() == bondCountBeforeBc, "undo removes the created bond");
    doc.redo();
    CHECK(doc.molecule().atomIds().size() == atomCountBeforeBc + 2, "redo re-creates both atoms");
    CHECK(doc.molecule().bondIds().size() == bondCountBeforeBc + 1, "redo re-creates the bond");

    // selectSingleItem: priority order and the all-unset-clears case.
    AtomId selA = doc.addAtom(QStringLiteral("C"), 50, 50);
    AtomId selB = doc.addAtom(QStringLiteral("C"), 51, 50);
    BondId selBond = doc.addBond(selA, selB, 1);
    doc.selectSingleItem(selA, selBond); // atom set AND bond set -> atom wins
    CHECK(doc.selection().atoms.contains(selA) && doc.selection().atoms.size() == 1,
          "selectSingleItem: atom wins when both atom and bond are set");
    CHECK(doc.selection().bonds.isEmpty(), "selectSingleItem: bond NOT selected when atom also set");
    doc.selectSingleItem(-1, selBond); // bond set, atom unset -> bond wins
    CHECK(doc.selection().bonds.contains(selBond) && doc.selection().bonds.size() == 1,
          "selectSingleItem: bond wins when atom is unset");
    CHECK(doc.selection().atoms.isEmpty(), "selectSingleItem: no atom selected in the bond-only case");
    doc.selectSingleItem(); // everything unset -> clears
    CHECK(doc.selection().atoms.isEmpty() && doc.selection().bonds.isEmpty(),
          "selectSingleItem: all-unset clears the selection");

    // setAtomQueryList / clearAtomQueryList.
    doc.setAtomQueryList(a2, QStringLiteral("C,N"), false);
    CHECK(doc.molecule().hasAtomQueryList(a2), "setAtomQueryList sets the query list");
    doc.undo();
    CHECK(!doc.molecule().hasAtomQueryList(a2), "undo removes the query list");
    CHECK(doc.molecule().atomSymbol(a2) == QStringLiteral("O"), "undo restores the original label too");
}

static void test_discreteTransforms() {
    std::printf("--- Test 6: discrete 90-degree rotate / flip ---\n");

    // NOTE on seeding: these tests create entities via doc.molecule().addAtom(...)
    // rather than doc.addAtom(...). The latter is itself undoable, which would
    // leave history non-empty and make canUndo() -- a plain bool -- useless for
    // asserting "this operation pushed exactly one entry". Seeding through the
    // molecule directly keeps history empty, so canUndo() is an exact signal.

    // Two atoms at (0,0) and (2,0) -> centroid (1,0).
    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        CHECK(!doc.canUndo(), "history starts empty (atoms seeded via the molecule)");

        doc.rotateSelection90CW();
        CHECK(doc.canUndo(), "rotate CW pushes exactly one history entry");
        double x = 0, y = 0;
        // rotate_cw: (cx + (y-cy), cy - (x-cx)); for a1 (0,0): (1+0, 0-(-1)) = (1,1)
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 1 && y == 1, "rotate CW moves a1 to (1,1)");
        // for a2 (2,0): (1+0, 0-(1)) = (1,-1)
        CHECK(doc.molecule().atomPos(a2, x, y) && x == 1 && y == -1, "rotate CW moves a2 to (1,-1)");
        doc.undo();
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 0 && y == 0, "undo restores a1 exactly");
        CHECK(doc.molecule().atomPos(a2, x, y) && x == 2 && y == 0, "undo restores a2 exactly");
        CHECK(!doc.canUndo(), "history is empty again after undoing the only entry");
    }

    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);

        doc.rotateSelection90CCW();
        double x = 0, y = 0;
        // rotate_ccw: (cx - (y-cy), cy + (x-cx)); for a1 (0,0): (1-0, 0+(-1)) = (1,-1)
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 1 && y == -1, "rotate CCW moves a1 to (1,-1)");
        CHECK(doc.molecule().atomPos(a2, x, y) && x == 1 && y == 1, "rotate CCW moves a2 to (1,1)");
    }

    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 3);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 5);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        // centroid (1,4). flip_h: (2*1 - x, y). flip_v: (x, 2*4 - y).
        doc.flipSelectionHorizontal();
        double x = 0, y = 0;
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 2 && y == 3, "flip H mirrors a1 across cx");
        CHECK(doc.molecule().atomPos(a2, x, y) && x == 0 && y == 5, "flip H mirrors a2 across cx");
        doc.undo();
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 0 && y == 3, "undo restores a1 after flip H");

        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        doc.flipSelectionVertical();
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 0 && y == 5, "flip V mirrors a1 across cy");
        CHECK(doc.molecule().atomPos(a2, x, y) && x == 2 && y == 3, "flip V mirrors a2 across cy");
    }

    // Fewer than 2 selected atoms: no-op, no history entry (real code's ids.length < 2 guard).
    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 7, 7);
        doc.selectAtom(a1);
        doc.rotateSelection90CW();
        CHECK(!doc.canUndo(), "single-atom selection pushes no history entry");
        double x = 0, y = 0;
        CHECK(doc.molecule().atomPos(a1, x, y) && x == 7 && y == 7, "single-atom selection is unmoved");

        doc.clearSelection();
        doc.flipSelectionHorizontal();
        CHECK(!doc.canUndo(), "empty selection pushes no history entry");
    }
}

static void test_imageTransforms() {
    std::printf("--- Test 7: image move / resize ---\n");
    DocumentState doc;
    ImageId img = doc.molecule().addImage(1, 2, 4, 3, QByteArray("fakepng"));

    double x = 0, y = 0, w = 0, h = 0;
    doc.moveImage(img, 2.5, -1.5);
    CHECK(doc.molecule().imageRect(img, x, y, w, h) && x == 3.5 && y == 0.5,
          "moveImage offsets the image position");
    CHECK(w == 4 && h == 3, "moveImage leaves the size untouched");
    doc.undo();
    CHECK(doc.molecule().imageRect(img, x, y, w, h) && x == 1 && y == 2,
          "undo restores the original image position exactly");

    doc.resizeImage(img, 2.0);
    CHECK(doc.molecule().imageRect(img, x, y, w, h) && w == 8 && h == 6,
          "resizeImage scales the image size");
    CHECK(x == 1 && y == 2, "resizeImage leaves the position untouched");
    doc.undo();
    CHECK(doc.molecule().imageRect(img, x, y, w, h) && w == 4 && h == 3,
          "undo restores the original image size exactly");

    // Both undos above emptied the history, so !canUndo() is an exact signal here.
    CHECK(!doc.canUndo(), "history is empty after undoing both image commands");

    // scaleFactor 0 would make invert's 1/0 non-finite: guarded no-op.
    doc.resizeImage(img, 0.0);
    CHECK(!doc.canUndo(), "resizeImage(0) pushes no history entry");
    CHECK(doc.molecule().imageRect(img, x, y, w, h) && w == 4 && h == 3,
          "resizeImage(0) leaves the size unchanged");

    // Invalid ids no-op cleanly.
    doc.moveImage(9999, 1, 1);
    doc.resizeImage(9999, 2.0);
    CHECK(!doc.canUndo(), "invalid ImageId pushes no history entry");
}

static void test_liveDragMove() {
    std::printf("--- Test 8: live-drag move gesture ---\n");
    // Seeded via doc.molecule() so history starts empty -- see the seeding note
    // in test_discreteTransforms.
    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);

    doc.moveSelectionLive(1, 1);
    doc.moveSelectionLive(0.5, -0.25);
    double x = 0, y = 0;
    CHECK(doc.molecule().atomPos(a1, x, y) && x == 1.5 && y == 0.75,
          "live move accumulates increments onto current positions");
    CHECK(!doc.canUndo(), "live move calls push NO history entry");

    doc.commitMove();
    CHECK(doc.canUndo(), "commitMove pushes exactly one history entry");
    CHECK(doc.molecule().atomPos(a1, x, y) && x == 1.5 && y == 0.75,
          "commitMove leaves the live-applied position in place (isFirst guard)");

    doc.undo();
    CHECK(doc.molecule().atomPos(a1, x, y) && x == 0 && y == 0, "undo restores a1 exactly");
    CHECK(doc.molecule().atomPos(a2, x, y) && x == 2 && y == 0, "undo restores a2 exactly");
    doc.redo();
    CHECK(doc.molecule().atomPos(a1, x, y) && x == 1.5 && y == 0.75, "redo re-applies the move");

    // commitMove with no movement pushes nothing.
    {
        DocumentState d2;
        AtomId b1 = d2.molecule().addAtom(QStringLiteral("C"), 0, 0);
        d2.selectAtom(b1);
        d2.commitMove();
        CHECK(!d2.canUndo(), "commitMove with zero delta pushes no history entry");
    }

    // Page-boundary clamp: the whole selection stops together at the edge.
    {
        DocumentState d3;
        AtomId c1 = d3.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId c2 = d3.molecule().addAtom(QStringLiteral("C"), 10, 0);
        d3.selectAtom(c1);
        d3.addAtomToSelection(c2);
        d3.moveSelectionLive(1000, 0);   // far past kPageMaxX (30)
        double cx1 = 0, cy1 = 0, cx2 = 0, cy2 = 0;
        d3.molecule().atomPos(c1, cx1, cy1);
        d3.molecule().atomPos(c2, cx2, cy2);
        CHECK(cx2 == 30, "clamped: the rightmost atom lands exactly on the page edge");
        CHECK(cx1 == 20, "clamped: the selection keeps its shape (10 units apart)");
    }

    // rxnArrow / rxnPlus / multitailArrow all translate with the selection.
    {
        DocumentState d4;
        RxnArrowId ar = d4.molecule().addRxnArrow(0, 0, 5, 0);
        RxnPlusId pl = d4.molecule().addRxnPlus(1, 1);
        MultitailArrowId mta = d4.molecule().addMultitailArrow({0, 0, 2, 2});
        d4.selectRxnArrow(ar);
        d4.addRxnPlusToSelection(pl);
        d4.addMultitailArrowToSelection(mta);

        d4.moveSelectionLive(3, 4);
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        CHECK(d4.molecule().rxnArrowEndpoints(ar, x1, y1, x2, y2) && x1 == 3 && y1 == 4 && x2 == 8 && y2 == 4,
              "live move translates BOTH rxnArrow endpoints");
        double px = 0, py = 0;
        CHECK(d4.molecule().rxnPlusPos(pl, px, py) && px == 4 && py == 5, "live move translates the rxnPlus");
        QList<double> pts = d4.molecule().multitailArrowPoints(mta);
        CHECK(pts.size() == 4 && pts[0] == 3 && pts[1] == 4 && pts[2] == 5 && pts[3] == 6,
              "live move translates EVERY multitailArrow point");

        d4.commitMove();
        d4.undo();
        CHECK(d4.molecule().rxnPlusPos(pl, px, py) && px == 1 && py == 1, "undo restores the rxnPlus");
        pts = d4.molecule().multitailArrowPoints(mta);
        CHECK(pts.size() == 4 && pts[0] == 0 && pts[3] == 2, "undo restores the multitailArrow points");
    }
}

static void test_liveDragRotate() {
    std::printf("--- Test 9: live-drag rotate gesture ---\n");
    DocumentState doc;
    // Symmetric pair about the origin -> centroid (0,0), easy exact math.
    // Seeded via doc.molecule() so history starts empty.
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1, 0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), -1, 0);
    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);

    const double kPi = 3.14159265358979323846;
    doc.rotateSelectionLive(kPi / 2);   // +90 degrees about (0,0): (1,0) -> (0,1)
    double x = 0, y = 0;
    CHECK(doc.molecule().atomPos(a1, x, y) && std::fabs(x - 0.0) < 1e-9 && std::fabs(y - 1.0) < 1e-9,
          "live rotate by +90 moves (1,0) to (0,1)");
    CHECK(!doc.canUndo(), "live rotate calls push NO history entry");

    // A second incremental call accumulates: another +90 => 180 total.
    doc.rotateSelectionLive(kPi / 2);
    CHECK(doc.molecule().atomPos(a1, x, y) && std::fabs(x + 1.0) < 1e-9 && std::fabs(y - 0.0) < 1e-9,
          "a second increment accumulates to 180 degrees total");

    doc.commitRotate();
    CHECK(doc.canUndo(), "commitRotate pushes exactly one history entry");
    doc.undo();
    CHECK(doc.molecule().atomPos(a1, x, y) && std::fabs(x - 1.0) < 1e-9 && std::fabs(y - 0.0) < 1e-9,
          "undo restores a1 to its exact original position");
    CHECK(doc.molecule().atomPos(a2, x, y) && std::fabs(x + 1.0) < 1e-9 && std::fabs(y - 0.0) < 1e-9,
          "undo restores a2 to its exact original position");
    doc.redo();
    CHECK(doc.molecule().atomPos(a1, x, y) && std::fabs(x + 1.0) < 1e-9, "redo re-applies the rotation");

    // Fewer than 2 snapshot points: no-op (real code's pts.length < 2 guard).
    {
        DocumentState d2;
        AtomId b1 = d2.molecule().addAtom(QStringLiteral("C"), 5, 5);
        d2.selectAtom(b1);
        d2.rotateSelectionLive(kPi / 2);
        double bx = 0, by = 0;
        CHECK(d2.molecule().atomPos(b1, bx, by) && bx == 5 && by == 5, "single point: rotate is a no-op");
        d2.commitRotate();
        CHECK(!d2.canUndo(), "commitRotate after a no-op gesture pushes nothing");
    }

    // Page-boundary clamp: a rotation that would push a point off-page freezes.
    {
        DocumentState d3;
        AtomId c1 = d3.molecule().addAtom(QStringLiteral("C"), 29, 0);
        AtomId c2 = d3.molecule().addAtom(QStringLiteral("C"), -29, 0);
        d3.selectAtom(c1);
        d3.addAtomToSelection(c2);
        // 90 degrees would put both at y = +-29, outside kPageMaxY (21).
        d3.rotateSelectionLive(kPi / 2);
        double cx = 0, cy = 0;
        CHECK(d3.molecule().atomPos(c1, cx, cy) && cx == 29 && cy == 0,
              "out-of-page rotation is refused: positions unchanged");
    }

    // Multitail arrows are excluded from rotate but do not block it.
    {
        DocumentState d4;
        AtomId e1 = d4.molecule().addAtom(QStringLiteral("C"), 1, 0);
        AtomId e2 = d4.molecule().addAtom(QStringLiteral("C"), -1, 0);
        MultitailArrowId mta = d4.molecule().addMultitailArrow({4, 4, 6, 6});
        d4.selectAtom(e1);
        d4.addAtomToSelection(e2);
        d4.addMultitailArrowToSelection(mta);
        d4.rotateSelectionLive(kPi / 2);
        double ex = 0, ey = 0;
        CHECK(d4.molecule().atomPos(e1, ex, ey) && std::fabs(ey - 1.0) < 1e-9,
              "atoms still rotate when a multitail arrow is also selected");
        QList<double> pts = d4.molecule().multitailArrowPoints(mta);
        CHECK(pts.size() == 4 && pts[0] == 4 && pts[3] == 6,
              "multitail arrow points are NOT rotated (documented limitation)");
    }
}

static void test_liveDragScale() {
    std::printf("--- Test 10: live-drag scale gesture ---\n");
    // Seeded via doc.molecule() so history starts empty.
    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);

    doc.scaleSelectionLive(2.0, 0.0, 0.0);   // anchor at the origin
    double x = 0, y = 0;
    CHECK(doc.molecule().atomPos(a1, x, y) && x == 0 && y == 0, "the anchor point itself stays fixed");
    CHECK(doc.molecule().atomPos(a2, x, y) && x == 4 && y == 0, "scale 2x doubles the distance from the anchor");
    CHECK(!doc.canUndo(), "live scale calls push NO history entry");

    // factor is ABSOLUTE, not incremental: 3.0 means 3x from gesture start,
    // NOT 2x then 3x (= 6x). This is the real code's semantics.
    doc.scaleSelectionLive(3.0, 0.0, 0.0);
    CHECK(doc.molecule().atomPos(a2, x, y) && x == 6 && y == 0,
          "a second call with factor 3 gives 3x from ORIGINAL, not 6x");

    doc.commitScale();
    CHECK(doc.canUndo(), "commitScale pushes exactly one history entry");
    doc.undo();
    CHECK(doc.molecule().atomPos(a2, x, y) && x == 2 && y == 0, "undo restores the exact original position");
    doc.redo();
    CHECK(doc.molecule().atomPos(a2, x, y) && x == 6 && y == 0, "redo re-applies the scale");

    // Factor is floored at 0.05 to stop a handle dragged through its anchor.
    {
        DocumentState d2;
        AtomId b1 = d2.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId b2 = d2.molecule().addAtom(QStringLiteral("C"), 10, 0);
        d2.selectAtom(b1);
        d2.addAtomToSelection(b2);
        d2.scaleSelectionLive(-5.0, 0.0, 0.0);   // negative would invert/collapse
        double bx = 0, by = 0;
        CHECK(d2.molecule().atomPos(b2, bx, by) && std::fabs(bx - 0.5) < 1e-9,
              "negative factor is floored to 0.05 (10 * 0.05 = 0.5)");
    }

    // Fewer than 2 snapshot points: no-op.
    {
        DocumentState d3;
        AtomId c1 = d3.molecule().addAtom(QStringLiteral("C"), 5, 5);
        d3.selectAtom(c1);
        d3.scaleSelectionLive(2.0, 0.0, 0.0);
        double cx = 0, cy = 0;
        CHECK(d3.molecule().atomPos(c1, cx, cy) && cx == 5 && cy == 5, "single point: scale is a no-op");
        d3.commitScale();
        CHECK(!d3.canUndo(), "commitScale after a no-op gesture pushes nothing");
    }

    // Page-boundary clamp: a scale that would push a point off-page freezes.
    {
        DocumentState d4;
        AtomId e1 = d4.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId e2 = d4.molecule().addAtom(QStringLiteral("C"), 20, 0);
        d4.selectAtom(e1);
        d4.addAtomToSelection(e2);
        d4.scaleSelectionLive(10.0, 0.0, 0.0);   // would put e2 at x=200, past kPageMaxX
        double ex = 0, ey = 0;
        CHECK(d4.molecule().atomPos(e2, ex, ey) && ex == 20 && ey == 0,
              "out-of-page scale is refused: positions unchanged");
    }
}

static void test_dragStateResetOnUndoRedo() {
    std::printf("--- Test 11: undo/redo reset ALL drag state (deliberate deviation) ---\n");
    const double kPi = 3.14159265358979323846;

    // Each case seeds atoms via doc.molecule() (no history), then creates
    // EXACTLY ONE real history entry with changeAtomCharge. After the undo,
    // history is empty, so `!canUndo()` after the commit is an exact signal:
    // if the drag state had survived, the commit WOULD push an entry and
    // canUndo() would flip to true. An equality-against-a-captured-bool check
    // cannot detect that, because both states are "true".

    // Move: an uncommitted gesture must not survive an undo.
    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        doc.changeAtomCharge(a1, 1);           // the one and only history entry
        doc.moveSelectionLive(1, 0);           // gesture in flight, uncommitted
        doc.undo();                            // undoes the charge AND must clear drag state
        CHECK(!doc.canUndo(), "history is empty after undoing the only entry");
        doc.commitMove();
        CHECK(!doc.canUndo(),
              "commitMove after undo pushes nothing -- move drag state was cleared");
    }

    // Rotate: the snapshot must be discarded too (the real JS never clears it).
    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1, 0);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), -1, 0);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        doc.changeAtomCharge(a1, 1);
        doc.rotateSelectionLive(kPi / 2);
        doc.undo();
        CHECK(!doc.canUndo(), "history is empty after undoing the only entry");
        doc.commitRotate();
        CHECK(!doc.canUndo(),
              "commitRotate after undo pushes nothing -- rotate drag state was cleared");
    }

    // Scale, via redo() rather than undo(). Here history is NOT empty after the
    // redo, so canUndo() cannot distinguish the two outcomes -- assert on WHAT
    // the next undo consumes instead.
    {
        DocumentState doc;
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0, 0);
        AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 2, 0);
        doc.selectAtom(a1);
        doc.addAtomToSelection(a2);
        doc.changeAtomCharge(a1, 1);
        doc.undo();                            // give redo something to do
        doc.selectAtom(a1);                    // undo cleared the selection
        doc.addAtomToSelection(a2);
        doc.scaleSelectionLive(2.0, 0.0, 0.0); // gesture in flight
        doc.redo();                            // re-applies the charge; must clear drag state
        CHECK(doc.molecule().atomCharge(a1) == 1, "redo re-applied the charge change");
        doc.commitScale();
        doc.undo();
        // If the scale state had survived, commitScale would have pushed a
        // command and this undo would have consumed IT, leaving charge at 1.
        CHECK(doc.molecule().atomCharge(a1) == 0,
              "undo consumed the charge entry, not a stale scale command -- scale drag state was cleared");
    }
}

static void test_addRing() {
    std::printf("--- Test 12: addRing (fusion-aware alternation) ---\n");
    const double kPi = 3.14159265358979323846;

    // Freestanding aromatic 6-ring: alternates 2/1/2/1/2/1.
    {
        DocumentState doc;
        QList<double> coords;
        for (int i = 0; i < 6; ++i) {
            double angle = (kPi / 3.0) * i;
            coords.append(std::cos(angle));
            coords.append(std::sin(angle));
        }
        CHECK(!doc.canUndo(), "history starts empty");
        doc.addRing(coords, true);
        CHECK(doc.canUndo(), "addRing pushes exactly one history entry");
        CHECK(doc.molecule().atomCount() == 6, "6 atoms created");
        CHECK(doc.molecule().bondCount() == 6, "6 bonds created");
        QList<AtomId> ids = doc.molecule().atomIds();
        for (int i = 0; i < 6; ++i) {
            BondId bid = doc.molecule().findBond(ids[i], ids[(i + 1) % 6]);
            CHECK(bid >= 0, "consecutive ring atoms are bonded");
        }
        doc.undo();
        CHECK(doc.molecule().atomCount() == 0, "undo removes every ring atom");
        CHECK(doc.molecule().bondCount() == 0, "undo removes every ring bond");
        CHECK(!doc.canUndo(), "history is empty again");
    }

    // Freestanding NonAromatic 6-ring stays all single.
    {
        DocumentState doc;
        QList<double> coords;
        for (int i = 0; i < 6; ++i) {
            double angle = (kPi / 3.0) * i;
            coords.append(std::cos(angle));
            coords.append(std::sin(angle));
        }
        doc.addRing(coords, false);
        QList<AtomId> ids = doc.molecule().atomIds();
        for (int i = 0; i < 6; ++i) {
            BondId bid = doc.molecule().findBond(ids[i], ids[(i + 1) % 6]);
            CHECK(doc.molecule().bondOrder(bid) == 1, "every bond stays single when aromatic=false");
        }
    }

    // A freestanding 5-ring stays all single even with aromatic=true (the
    // real code's `if (n !== 6) return` branch).
    {
        DocumentState doc;
        QList<double> coords;
        for (int i = 0; i < 5; ++i) {
            double angle = (2.0 * kPi / 5.0) * i;
            coords.append(std::cos(angle));
            coords.append(std::sin(angle));
        }
        doc.addRing(coords, true);
        QList<AtomId> ids = doc.molecule().atomIds();
        for (int i = 0; i < 5; ++i) {
            BondId bid = doc.molecule().findBond(ids[i], ids[(i + 1) % 5]);
            CHECK(doc.molecule().bondOrder(bid) == 1, "5-ring stays all single (n != 6 branch)");
        }
    }

    // FUSION SEAM REGRESSION (the exact naphthalene bug the real comment
    // describes): a benzene ring, then a second aromatic 6-ring sharing one
    // edge with it. The seam bond must keep its committed order, and the
    // new ring must alternate outward WITHOUT two adjacent double bonds at
    // the seam.
    {
        DocumentState doc;
        AtomId a0 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
        BondId seam = doc.molecule().addBond(a0, a1, 2);   // pre-existing, committed order 2
        CHECK(seam > 0, "setup: seam bond created");

        // Second ring shares the (a0,a1) edge: reuse their exact coordinates
        // for the first two ring vertices so mergeOverlappingAtoms fuses onto
        // them, then continue the hexagon outward.
        QList<double> coords;
        coords << 0.0 << 0.0 << 1.5 << 0.0;             // coincide with a0, a1
        coords << 2.25 << 1.3 << 1.5 << 2.6 << 0.0 << 2.6 << -0.75 << 1.3;
        doc.addRing(coords, true);

        // 2 pre-existing (a0,a1) + 6 new ring atoms - 2 fused away (the two
        // new atoms coinciding with a0/a1) = 6 total.
        CHECK(doc.molecule().atomCount() == 6, "2 pre-existing + 6 new - 2 fused onto a0/a1 = 6 total");
        // seam(a0-a1) + 2 replacement bonds (a0-to-new-neighbor,
        // a1-to-new-neighbor) + 3 untouched inner ring bonds = 6 total.
        CHECK(doc.molecule().bondCount() == 6, "6 bonds total in the fused bicyclic system");
        CHECK(doc.molecule().bondOrder(seam) == 2, "the seam bond's committed order is UNCHANGED");

        // Walk the new ring's atoms (everything except a0/a1) and confirm no
        // two bonds adjacent to the seam are both double (the naphthalene bug).
        QList<AtomId> allIds = doc.molecule().atomIds();
        QList<AtomId> newRingAtoms;
        for (AtomId id : allIds) if (id != a0 && id != a1) newRingAtoms.append(id);
        // Only 4, not 6: two of the new ring's 6 vertices coincided with
        // a0/a1 and were fused away, never surviving to this point.
        CHECK(newRingAtoms.size() == 4, "4 surviving non-seam atoms belong to the new ring");
        // The two edges leaving the seam (a0-to-new-neighbor and
        // a1-to-new-neighbor) must both be SINGLE -- alternating outward from
        // an already-double seam.
        int singleCount = 0, doubleCount = 0;
        for (AtomId id : newRingAtoms) {
            BondId toA0 = doc.molecule().findBond(a0, id);
            BondId toA1 = doc.molecule().findBond(a1, id);
            if (toA0 >= 0) { if (doc.molecule().bondOrder(toA0) == 1) ++singleCount; else ++doubleCount; }
            if (toA1 >= 0) { if (doc.molecule().bondOrder(toA1) == 1) ++singleCount; else ++doubleCount; }
        }
        CHECK(singleCount == 2 && doubleCount == 0,
              "both edges leaving the seam are single -- no adjacent double bonds at the fusion point");

        doc.undo();
        CHECK(doc.molecule().atomCount() == 2, "undo removes the 4 surviving new-ring atoms, leaving a0/a1");
        CHECK(doc.molecule().bondCount() == 1, "undo removes exactly the new-ring bonds, leaving the seam");
        CHECK(doc.molecule().bondOrder(seam) == 2, "the seam bond survives undo with its original order");
    }

    // Malformed coords no-op cleanly. Each case isolates ONE guard condition
    // so a bug in either check can't hide behind the other also tripping.
    {
        DocumentState doc;
        // Odd length but >= 6 doubles (7 total: 3.5 "points") -- exercises
        // the `% 2 != 0` check specifically, independent of the size check.
        doc.addRing({0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0});
        CHECK(!doc.canUndo(), "odd-length coords (>=6 doubles): no history entry");
        // Even length but only 2 points (4 doubles) -- exercises the `< 6`
        // check specifically, independent of the oddness check.
        doc.addRing({0.0, 0.0, 1.0, 0.0});
        CHECK(!doc.canUndo(), "fewer than 3 points (even length): no history entry");
    }
}

static void test_addChain() {
    std::printf("--- Test 13: addChain ---\n");

    // Horizontal chain, dist=6, bondLen=1.5 -> nBonds = round(6/1.5) = 4 ->
    // 5 atoms, 4 bonds. theta=0, half=PI/6.
    {
        DocumentState doc;
        CHECK(!doc.canUndo(), "history starts empty");
        doc.addChain(0.0, 0.0, 6.0, 0.0);
        CHECK(doc.canUndo(), "addChain pushes exactly one history entry");
        CHECK(doc.molecule().atomCount() == 5, "5 atoms for a 4-bond chain");
        CHECK(doc.molecule().bondCount() == 4, "4 bonds in the chain");

        // Hand-computed expected vertices (theta=0, half=PI/6, bondLen=1.5),
        // matching ChainPlacementEngine::compute's exact formula.
        const double kPi = 3.14159265358979323846;
        const double half = kPi / 6.0;
        const double bondLen = 1.5;
        double px = 0.0, py = 0.0;
        QList<AtomId> ids = doc.molecule().atomIds();
        CHECK(ids.size() == 5, "setup: 5 atom ids to check");
        for (int i = 0; i <= 4; ++i) {
            double x = 0, y = 0;
            CHECK(doc.molecule().atomPos(ids[i], x, y), "atom position readable");
            // 1e-5, not 1e-9: EditableMolecule::addAtom stores coordinates via
            // indigoSetXYZ, which takes float (32-bit), not double -- exact
            // round numbers (0.0, 1.0, -0.5, etc.) survive that round-trip
            // losslessly, but these expected values are irrational trig
            // results, which lose precision below float32's ~7 significant
            // digits. 1e-5 is comfortably above that rounding error for
            // values in this test's range (0 to ~5.2).
            CHECK(std::fabs(x - px) < 1e-5 && std::fabs(y - py) < 1e-5, "vertex position matches the formula");
            double ang = 0.0 + ((i % 2 == 0) ? half : -half);
            px += bondLen * std::cos(ang);
            py += bondLen * std::sin(ang);
        }

        doc.undo();
        CHECK(doc.molecule().atomCount() == 0, "undo removes every chain atom");
        CHECK(doc.molecule().bondCount() == 0, "undo removes every chain bond");
    }

    // A zero-length drag still lays at least one bond (max(1, ...)).
    {
        DocumentState doc;
        doc.addChain(2.0, 2.0, 2.0, 2.0);
        CHECK(doc.molecule().atomCount() == 2, "zero-length drag still produces a 2-atom chain");
        CHECK(doc.molecule().bondCount() == 1, "zero-length drag still produces 1 bond");
    }

    // Chain endpoint placed exactly on an existing atom grafts onto it
    // (total atom count reflects the merge, not a duplicate).
    {
        DocumentState doc;
        AtomId existing = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        doc.addChain(0.0, 0.0, 3.0, 0.0);   // dist=3 -> nBonds=2 -> 3 atoms, first one coincides
        // 1 pre-existing + 3 new - 1 fused = 3 total.
        CHECK(doc.molecule().atomCount() == 3, "chain endpoint fuses onto the existing atom");
        CHECK(doc.molecule().atomSymbol(existing) == QStringLiteral("C"), "the existing atom survives");
        int chainBondsFromExisting = 0;
        for (BondId bid : doc.molecule().bondIds()) {
            AtomId a = -1, b = -1;
            doc.molecule().bondEndpoints(bid, a, b);
            if (a == existing || b == existing) ++chainBondsFromExisting;
        }
        CHECK(chainBondsFromExisting == 1, "the existing atom gained exactly one chain bond");
    }
}

static void test_importReaction() {
    std::printf("--- Test: importReaction ---\n");
    DocumentState doc("");

    bool ok = doc.importReaction(QStringLiteral("CCO>>CC=O"));
    CHECK(ok, "importReaction succeeds on a valid reaction SMILES");
    CHECK(doc.molecule().atomIds().size() >= 5, "importReaction created atoms for both sides");
    CHECK(doc.molecule().rxnArrowIds().size() == 1, "importReaction added exactly one RxnArrow");
    CHECK(doc.canUndo(), "importReaction pushes at least one history entry");

    int atomsBefore = doc.molecule().atomIds().size();
    int arrowsBefore = doc.molecule().rxnArrowIds().size();
    doc.undo();
    CHECK(doc.molecule().atomIds().size() == 0, "undo removes every atom importReaction created");
    CHECK(doc.molecule().rxnArrowIds().size() == 0, "undo removes the RxnArrow too");
    (void)atomsBefore; (void)arrowsBefore;

    DocumentState doc2("");
    bool ok2 = doc2.importReaction(QStringLiteral("not a reaction"));
    CHECK(!ok2, "importReaction rejects unparseable text");
    CHECK(doc2.molecule().atomIds().size() == 0, "failed importReaction leaves no partial atoms");
    CHECK(!doc2.canUndo(), "failed importReaction pushes no history entry");

    DocumentState doc3("");
    bool ok3 = doc3.importReaction(QStringLiteral("CCO.CC(=O)O>>CC=O.O"));
    
    QList<QList<AtomId>> components;
    QList<AtomId> unvisited = doc3.molecule().atomIds();
    while (!unvisited.isEmpty()) {
        QList<AtomId> comp;
        QList<AtomId> queue;
        queue.append(unvisited.takeFirst());
        while (!queue.isEmpty()) {
            AtomId curr = queue.takeFirst();
            comp.append(curr);
            for (BondId bid : doc3.molecule().bondIds()) {
                AtomId a = -1, b = -1;
                doc3.molecule().bondEndpoints(bid, a, b);
                AtomId neighbor = -1;
                if (a == curr) neighbor = b;
                else if (b == curr) neighbor = a;
                if (neighbor != -1 && unvisited.contains(neighbor)) {
                    unvisited.removeAll(neighbor);
                    queue.append(neighbor);
                }
            }
        }
        components.append(comp);
    }
    
    bool valid = ok3 && components.size() == 4;
    bool interFragmentOverlap = false;
    if (valid) {
        for (int i = 0; i < components.size(); ++i) {
            for (int j = i + 1; j < components.size(); ++j) {
                for (AtomId a1 : components[i]) {
                    for (AtomId a2 : components[j]) {
                        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
                        doc3.molecule().atomPos(a1, x1, y1);
                        doc3.molecule().atomPos(a2, x2, y2);
                        double dx = x1 - x2;
                        double dy = y1 - y2;
                        if (dx * dx + dy * dy < 0.1) {
                            interFragmentOverlap = true;
                        }
                    }
                }
            }
        }
    }
    CHECK(valid && !interFragmentOverlap, "fragments in 2-reactant reaction do not overlap each other");
}

static void test_importReactionRxnPlusGapCentering() {
    std::printf("--- Test: importReaction centers RxnPlus on the visual gap, not fragment centers ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("O.c1ccccc1>>CC=O"));
    CHECK(ok, "importReaction succeeds on an asymmetric-width 2-reactant reaction");

    QList<RxnPlusId> pluses = doc.molecule().rxnPlusIds();
    CHECK(pluses.size() == 1, "exactly one RxnPlus between the two reactants (product side has only one fragment)");

    // Group atoms into connected components (same BFS approach as test_importReaction's own
    // 2-reactant scenario above) to find each reactant fragment's own bounding box directly
    // from its actual placed atom positions, without depending on Indigo's exact layout width
    // for benzene.
    QList<QList<AtomId>> components;
    QList<AtomId> unvisited = doc.molecule().atomIds();
    while (!unvisited.isEmpty()) {
        QList<AtomId> comp;
        QList<AtomId> queue;
        queue.append(unvisited.takeFirst());
        while (!queue.isEmpty()) {
            AtomId curr = queue.takeFirst();
            comp.append(curr);
            for (BondId bid : doc.molecule().bondIds()) {
                AtomId a = -1, b = -1;
                doc.molecule().bondEndpoints(bid, a, b);
                AtomId neighbor = -1;
                if (a == curr) neighbor = b;
                else if (b == curr) neighbor = a;
                if (neighbor != -1 && unvisited.contains(neighbor)) {
                    unvisited.removeAll(neighbor);
                    queue.append(neighbor);
                }
            }
        }
        components.append(comp);
    }
    CHECK(components.size() == 3, "three fragments total: water, benzene, acetaldehyde");

    // Identify the two REACTANT fragments (the two leftmost by center-x) and compute each
    // one's own bounding box directly from its atoms' actual placed positions.
    struct Frag { double minX, maxX, centerX, effectiveHalfWidth; };
    QList<Frag> frags;
    for (const QList<AtomId>& comp : components) {
        double minX = 0, maxX = 0;
        bool any = false;
        for (AtomId id : comp) {
            double x = 0, y = 0;
            doc.molecule().atomPos(id, x, y);
            if (!any) { minX = maxX = x; any = true; }
            else { minX = std::min(minX, x); maxX = std::max(maxX, x); }
        }
        // DocumentState::kBondLength is a private member and not reachable from this free
        // function, so its value (1.5, "chem-core's StandardBondLength", DocumentState.h) is
        // hardcoded here directly, matching this test file's existing convention elsewhere of
        // hardcoding the bond length literal rather than referencing the private constant.
        double effectiveHalfWidth = std::max(maxX - minX, 1.5 * 0.6) / 2.0;
        frags.append({minX, maxX, (minX + maxX) / 2.0, effectiveHalfWidth});
    }
    std::sort(frags.begin(), frags.end(), [](const Frag& a, const Frag& b) { return a.centerX < b.centerX; });
    const Frag& left = frags[0];
    const Frag& right = frags[1];

    double leftEffectiveEdge = left.centerX + left.effectiveHalfWidth;
    double rightEffectiveEdge = right.centerX - right.effectiveHalfWidth;
    double correctGapMidpoint = (leftEffectiveEdge + rightEffectiveEdge) / 2.0;
    double naiveCenterMidpoint = (left.centerX + right.centerX) / 2.0;
    // This sanity check's threshold was originally calibrated against the raw atom-bounding-box
    // edges (before the single-atom-fragment width floor existed), where the only source of
    // divergence from naiveCenterMidpoint was benzene's own real half-width (~1.0 here). Now that
    // correctGapMidpoint is computed from EFFECTIVE (floor-aware) edges to match production's own
    // formula, water's own effective half-width (0.45, from the 0.9 floor) also contributes to
    // leftEffectiveEdge, which legitimately narrows the gap between the two candidate formulas
    // (measured here at ~0.275 with the real benzene geometry this test produces) -- this is the
    // same width-floor change, not a new regression, so the threshold is lowered rather than the
    // check being removed: 0.1 is still comfortably above zero, so this remains a genuine
    // guard against the degenerate case where the two formulas would coincide by coincidence.
    CHECK(std::abs(correctGapMidpoint - naiveCenterMidpoint) > 0.1,
          "setup sanity check: water (1 atom, narrow) and benzene (6 atoms, wide) are asymmetric enough that the two candidate formulas genuinely differ");

    double plusX = 0, plusY = 0;
    doc.molecule().rxnPlusPos(pluses[0], plusX, plusY);
    CHECK(std::abs(plusX - correctGapMidpoint) < 0.01,
          "RxnPlus sits at the true edge-to-edge gap midpoint, not the naive center-to-center midpoint");
}

static void test_importReactionRxnPlusGapCenteringSingleAtomFloor() {
    std::printf("--- Test: importReaction's RxnPlus gap-centering respects the single-atom-fragment width floor ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("O.c1ccccc1>>CC=O"));
    CHECK(ok, "importReaction succeeds on an asymmetric-width 2-reactant reaction");

    QList<RxnPlusId> pluses = doc.molecule().rxnPlusIds();
    CHECK(pluses.size() == 1, "exactly one RxnPlus between the two reactants");

    // Group atoms into connected components (same BFS approach as
    // test_importReactionRxnPlusGapCentering) to find benzene's own real placed bounding box.
    QList<QList<AtomId>> components;
    QList<AtomId> unvisited = doc.molecule().atomIds();
    while (!unvisited.isEmpty()) {
        QList<AtomId> comp;
        QList<AtomId> queue;
        queue.append(unvisited.takeFirst());
        while (!queue.isEmpty()) {
            AtomId curr = queue.takeFirst();
            comp.append(curr);
            for (BondId bid : doc.molecule().bondIds()) {
                AtomId a = -1, b = -1;
                doc.molecule().bondEndpoints(bid, a, b);
                AtomId neighbor = -1;
                if (a == curr) neighbor = b;
                else if (b == curr) neighbor = a;
                if (neighbor != -1 && unvisited.contains(neighbor)) {
                    unvisited.removeAll(neighbor);
                    queue.append(neighbor);
                }
            }
        }
        components.append(comp);
    }
    CHECK(components.size() == 3, "three fragments total: water, benzene, acetaldehyde");

    // Water is a single atom: find the 1-atom component, and its own real placed X (a
    // single-atom fragment's placement translation always lands the atom exactly at its own
    // fragment center, regardless of what width value was used for spacing -- so this position
    // alone doesn't reveal the floor; it's the input the algebra below needs).
    double waterX = 0;
    bool foundWater = false;
    QList<AtomId> benzeneComp;
    for (const QList<AtomId>& comp : components) {
        if (comp.size() == 1) {
            double x = 0, y = 0;
            doc.molecule().atomPos(comp[0], x, y);
            waterX = x;
            foundWater = true;
        } else if (comp.size() == 6) {
            benzeneComp = comp;
        }
    }
    CHECK(foundWater, "found the single-atom water fragment");
    CHECK(!benzeneComp.isEmpty(), "found the 6-atom benzene fragment");

    double benzeneMinX = 0, benzeneMaxX = 0;
    bool any = false;
    for (AtomId id : benzeneComp) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        if (!any) { benzeneMinX = benzeneMaxX = x; any = true; }
        else { benzeneMinX = std::min(benzeneMinX, x); benzeneMaxX = std::max(benzeneMaxX, x); }
    }
    double benzeneRightEdge = benzeneMinX;   // benzene sits to the RIGHT of water, so its LEFT
                                              // edge (benzeneMinX) is the gap-formula's rightEdge

    double plusX = 0, plusY = 0;
    doc.molecule().rxnPlusPos(pluses[0], plusX, plusY);

    // Gap formula: plusX = (leftEdge_water + rightEdge_benzene) / 2, where
    // leftEdge_water = waterX + waterEffectiveWidth/2 (water's own placement always lands it
    // exactly at its fragment center, so waterX IS reactantCenters[0] from the production code's
    // own math). Solve for waterEffectiveWidth:
    double leftEdgeWater = 2.0 * plusX - benzeneRightEdge;
    double waterEffectiveWidth = 2.0 * (leftEdgeWater - waterX);

    CHECK(waterEffectiveWidth > 0.5,
          "water's effective width (implied by the real RxnPlus position and benzene's real bounding box) "
          "reflects the minimum-width floor (~0.9), not the raw zero-width atom-coordinate bounding box "
          "that would be implied without the fix");
}

static void test_importReactionRecentersOnPageOrigin() {
    std::printf("--- Test: importReaction recenters a wide assembly onto the page origin (does not clamp/shrink it to fit) ---\n");
    QStringList reactantAtoms;
    for (int i = 0; i < 50; ++i) reactantAtoms.append(QStringLiteral("C"));
    QString reactionText = reactantAtoms.join(QStringLiteral(".")) + QStringLiteral(">>O");

    DocumentState doc("");
    bool ok = doc.importReaction(reactionText);
    CHECK(ok, "importReaction succeeds on a wide 50-reactant reaction");

    double minX = 0, maxX = 0;
    bool any = false;
    for (AtomId id : doc.molecule().atomIds()) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        if (!any) { minX = maxX = x; any = true; }
        else { minX = std::min(minX, x); maxX = std::max(maxX, x); }
    }
    for (RxnArrowId aid : doc.molecule().rxnArrowIds()) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        doc.molecule().rxnArrowEndpoints(aid, x1, y1, x2, y2);
        minX = std::min({minX, x1, x2});
        maxX = std::max({maxX, x1, x2});
    }
    CHECK(any, "reaction actually placed atoms");
    CHECK(maxX - minX > 60.0, "setup sanity check: this reaction's natural un-clamped width genuinely exceeds the page's own 60-unit span");

    double assemblyCenter = (minX + maxX) / 2.0;
    CHECK(std::abs(assemblyCenter - 0.0) < 0.01,
          "the wide assembly's overall center is recentered onto the page origin (0,0), not left wherever the raw left-to-right layout happened to start -- note this assembly is deliberately wider than the page and is NOT clamped/shrunk to fit; only its center is repositioned");
}

static void test_importReactionSingleReactantSingleProductRecentering() {
    std::printf("--- Test: importReaction recenters even the trivial 1-reactant/1-product case ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("CCO>>CC=O"));
    CHECK(ok, "importReaction succeeds on a simple 1-reactant/1-product reaction");
    CHECK(doc.molecule().rxnPlusIds().isEmpty(), "no RxnPlus entities needed when there is only one fragment per side");

    double minX = 0, maxX = 0;
    bool any = false;
    for (AtomId id : doc.molecule().atomIds()) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        if (!any) { minX = maxX = x; any = true; }
        else { minX = std::min(minX, x); maxX = std::max(maxX, x); }
    }
    for (RxnArrowId aid : doc.molecule().rxnArrowIds()) {
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        doc.molecule().rxnArrowEndpoints(aid, x1, y1, x2, y2);
        minX = std::min({minX, x1, x2});
        maxX = std::max({maxX, x1, x2});
    }
    CHECK(any, "reaction actually placed atoms");

    double assemblyCenter = (minX + maxX) / 2.0;
    CHECK(std::abs(assemblyCenter - 0.0) < 0.01,
          "even the trivial 1-reactant/1-product case (no multi-fragment spacing math to exercise) is correctly recentered onto the page origin (0,0)");
}

static void test_importReactionWithCatalyst() {
    std::printf("--- Test: importReaction places a real catalyst fragment above the arrow ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("CCO>[Pd]>CC=O"));
    CHECK(ok, "importReaction succeeds on a reaction with a catalyst component");

    // "CCO" (ethanol, 3 heavy atoms) + "[Pd]" (1 atom) + "CC=O" (acetaldehyde, 3 heavy atoms) --
    // if the catalyst were silently dropped, only 6 atoms would exist; if it's genuinely placed,
    // there are 7.
    CHECK(doc.molecule().atomIds().size() == 7,
          "catalyst atom is genuinely placed, not silently dropped (3 reactant + 1 catalyst + 3 product)");

    // Find the catalyst atom: the one whose element symbol is "Pd" (unambiguous -- no other atom
    // in this reaction is palladium).
    AtomId catalystAtom = -1;
    for (AtomId id : doc.molecule().atomIds()) {
        if (doc.molecule().atomSymbol(id) == QStringLiteral("Pd")) catalystAtom = id;
    }
    CHECK(catalystAtom != -1, "found the catalyst atom by its element symbol");

    double catX = 0, catY = 0;
    doc.molecule().atomPos(catalystAtom, catX, catY);
    CHECK(catY < -0.01, "catalyst is placed above the arrow (negative Y, per this app's Y-down screen convention)");

    // Reactant/product atoms stay on the y=0 centerline (unaffected by the new catalyst row).
    bool anyReactantOrProductBelowThreshold = false;
    for (AtomId id : doc.molecule().atomIds()) {
        if (id == catalystAtom) continue;
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        if (std::abs(y) > 5.0) anyReactantOrProductBelowThreshold = true;   // sanity bound, not exact
    }
    CHECK(!anyReactantOrProductBelowThreshold, "reactant/product atoms remain near the y=0 centerline, unaffected by the catalyst row");
}

static void test_importReactionWithMultiAtomCatalyst() {
    std::printf("--- Test: importReaction keeps a tall multi-atom catalyst clear of the arrow/reactant row ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("CCO>P(c1ccccc1)(c1ccccc1)c1ccccc1>CC=O"));
    CHECK(ok, "importReaction succeeds on a reaction with a multi-atom catalyst");

    QList<QList<AtomId>> components;
    QList<AtomId> unvisited = doc.molecule().atomIds();
    while (!unvisited.isEmpty()) {
        QList<AtomId> comp;
        QList<AtomId> queue;
        queue.append(unvisited.takeFirst());
        while (!queue.isEmpty()) {
            AtomId curr = queue.takeFirst();
            comp.append(curr);
            for (BondId bid : doc.molecule().bondIds()) {
                AtomId a = -1, b = -1;
                doc.molecule().bondEndpoints(bid, a, b);
                AtomId neighbor = -1;
                if (a == curr) neighbor = b;
                else if (b == curr) neighbor = a;
                if (neighbor != -1 && unvisited.contains(neighbor)) {
                    unvisited.removeAll(neighbor);
                    queue.append(neighbor);
                }
            }
        }
        components.append(comp);
    }
    CHECK(components.size() == 3, "three fragments total: ethanol, catalyst, acetaldehyde");

    QList<AtomId> catalystComp;
    QList<AtomId> reactantProductAtoms;
    for (const QList<AtomId>& comp : components) {
        bool isCatalyst = false;
        for (AtomId id : comp) {
            if (doc.molecule().atomSymbol(id) == QStringLiteral("P")) isCatalyst = true;
        }
        if (isCatalyst) catalystComp = comp;
        else reactantProductAtoms += comp;
    }
    CHECK(!catalystComp.isEmpty(), "found the catalyst fragment (contains the P atom)");
    CHECK(!reactantProductAtoms.isEmpty(), "found the reactant/product atoms");

    double catalystMaxY = -1e9;
    for (AtomId id : catalystComp) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        catalystMaxY = std::max(catalystMaxY, y);
    }
    double reactantProductMinY = 1e9;
    for (AtomId id : reactantProductAtoms) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        reactantProductMinY = std::min(reactantProductMinY, y);
    }
    CHECK(catalystMaxY < reactantProductMinY,
          "the catalyst row's lowest (closest-to-zero) atom stays strictly above the reactant/product row's highest atom -- no vertical overlap even for a tall multi-atom catalyst");
}

static void test_importReactionWithTallReactant() {
    std::printf("--- Test: importReaction keeps a tall reactant/product fragment clear of the catalyst row ---\n");
    DocumentState doc("");
    bool ok = doc.importReaction(QStringLiteral("P(c1ccccc1)(c1ccccc1)c1ccccc1>[Pd]>CC=O"));
    CHECK(ok, "importReaction succeeds on a reaction with a tall reactant and an ordinary catalyst");

    QList<QList<AtomId>> components;
    QList<AtomId> unvisited = doc.molecule().atomIds();
    while (!unvisited.isEmpty()) {
        QList<AtomId> comp;
        QList<AtomId> queue;
        queue.append(unvisited.takeFirst());
        while (!queue.isEmpty()) {
            AtomId curr = queue.takeFirst();
            comp.append(curr);
            for (BondId bid : doc.molecule().bondIds()) {
                AtomId a = -1, b = -1;
                doc.molecule().bondEndpoints(bid, a, b);
                AtomId neighbor = -1;
                if (a == curr) neighbor = b;
                else if (b == curr) neighbor = a;
                if (neighbor != -1 && unvisited.contains(neighbor)) {
                    unvisited.removeAll(neighbor);
                    queue.append(neighbor);
                }
            }
        }
        components.append(comp);
    }
    CHECK(components.size() == 3, "three fragments total: PPh3 reactant, Pd catalyst, acetaldehyde product");

    QList<AtomId> catalystComp;
    QList<AtomId> otherAtoms;
    for (const QList<AtomId>& comp : components) {
        bool isCatalyst = (comp.size() == 1);
        for (AtomId id : comp) {
            if (doc.molecule().atomSymbol(id) == QStringLiteral("Pd")) isCatalyst = true;
        }
        if (isCatalyst) catalystComp = comp;
        else otherAtoms += comp;
    }
    CHECK(!catalystComp.isEmpty(), "found the catalyst fragment (the single Pd atom)");
    CHECK(!otherAtoms.isEmpty(), "found the reactant/product atoms");

    double catalystMaxY = -1e9;
    for (AtomId id : catalystComp) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        catalystMaxY = std::max(catalystMaxY, y);
    }
    double otherMinY = 1e9;
    for (AtomId id : otherAtoms) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        otherMinY = std::min(otherMinY, y);
    }
    CHECK(catalystMaxY < otherMinY,
          "the catalyst row stays strictly above even a tall reactant fragment -- no vertical overlap when the TALL fragment is in the reactant/product row instead of the catalyst row");
}

static void test_insertFunctionalGroup() {
    std::printf("--- Test 14: insertFunctionalGroup ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    // Plain paste: no target atom, centers on (cx, cy).
    {
        DocumentState doc;
        CHECK(!doc.canUndo(), "history starts empty");
        doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 10.0, 10.0);
        CHECK(doc.canUndo(), "insertFunctionalGroup pushes exactly one history entry");
        CHECK(doc.molecule().atomCount() == 3, "\"Ac\" plain paste creates 3 atoms");
        CHECK(doc.molecule().bondCount() == 2, "\"Ac\" plain paste creates 2 bonds");
        doc.undo();
        CHECK(doc.molecule().atomCount() == 0, "undo removes every pasted atom");
        CHECK(doc.molecule().bondCount() == 0, "undo removes every pasted bond");
    }

    // Graft onto an existing atom: "Ac"'s attach atom (SUP attachment point
    // at atom index 0, confirmed by direct file read) sprouts a new bond to
    // the target instead of fusing into it -- the target survives with
    // exactly one new bond, and the attach atom survives as its own atom,
    // still carrying both of its own original internal bonds.
    {
        DocumentState doc;
        AtomId target = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 0.0, 0.0, target, true);

        // "Ac" has 3 atoms; none are fused away anymore, so 1 pre-existing + 3 new = 4 total.
        CHECK(doc.molecule().atomCount() == 4, "graft: 1 pre-existing + 3 new (none fused) = 4 total");
        CHECK(doc.molecule().atomSymbol(target) == QStringLiteral("C"), "target atom survives");
        int bondsFromTarget = 0;
        AtomId attachAtom = -1;
        for (BondId bid : doc.molecule().bondIds()) {
            AtomId a = -1, b = -1;
            doc.molecule().bondEndpoints(bid, a, b);
            if (a == target || b == target) { ++bondsFromTarget; attachAtom = (a == target) ? b : a; }
        }
        // "Ac"'s attach atom (atom index 0, the carbonyl carbon) is bonded to BOTH other atoms
        // within the fragment (bond block: "2 1 1" and "1 3 2") -- but under sprout semantics
        // those two bonds stay the attach atom's OWN, they are never rewired onto the target.
        // The target gains exactly one new bond, to the attach atom itself.
        CHECK(bondsFromTarget == 1, "target atom gained exactly one new bond, to the attach atom (no rewiring)");
        CHECK(attachAtom != -1 && doc.molecule().atomSymbol(attachAtom) == QStringLiteral("C"),
              "the attach atom (Ac's carbonyl carbon) is the one bonded to the target");
        int attachBondCount = 0;
        for (BondId bid : doc.molecule().bondIds()) {
            AtomId a = -1, b = -1;
            doc.molecule().bondEndpoints(bid, a, b);
            if (a == attachAtom || b == attachAtom) ++attachBondCount;
        }
        CHECK(attachBondCount == 3, "attach atom keeps both of its own original template-internal bonds, plus the one new bond to the target");
        double tx = 0, ty = 0, ax = 0, ay = 0;
        doc.molecule().atomPos(target, tx, ty);
        doc.molecule().atomPos(attachAtom, ax, ay);
        CHECK(std::abs(std::hypot(ax - tx, ay - ty) - 1.5) < 0.01, "attach atom sits exactly one bond length (kBondLength) from the target");
        doc.undo();
        CHECK(doc.molecule().atomCount() == 1, "undo restores to just the pre-existing target atom");
    }

    // Missing template: placeholder-atom fallback.
    {
        DocumentState doc;
        doc.insertFunctionalGroup(lib, QStringLiteral("NotARealTemplateName"), 2.0, 3.0);
        CHECK(doc.molecule().atomCount() == 1, "missing template falls back to one placeholder atom");
        QList<AtomId> ids = doc.molecule().atomIds();
        CHECK(doc.molecule().atomSymbol(ids[0]) == QStringLiteral("NotARealTemplateName"),
              "placeholder atom is labeled with the requested (missing) template name");
        doc.undo();
        CHECK(doc.molecule().atomCount() == 0, "undo removes the placeholder atom");
    }
}

static void test_insertFunctionalGroupSGroupCleanupOnUndo() {
    std::printf("--- Test: insertFunctionalGroup's SGroup is genuinely gone after undo, not left dangling ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    AtomId target = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 0.0, 0.0, target, true);

    CHECK(doc.molecule().sgroupIds().size() == 1, "graft created exactly one sgroup");

    doc.undo();

    // The pre-existing target atom survives undo (only the grafted template's own atoms/bonds
    // are removed) -- matching test_insertFunctionalGroup's own established graft-undo pattern.
    CHECK(doc.molecule().atomIds().size() == 1, "pre-existing target atom survives undo");
    CHECK(doc.molecule().sgroupIds().isEmpty(),
          "the sgroup created by the graft is genuinely gone after undo, not left dangling -- "
          "confirms EditableMolecule::removeAtom's rebuildIndexTables() cascade already prunes "
          "stale m_sgroupIdx entries correctly, with no need for cmd.invert to call an explicit "
          "SGroup-removal step");
}

static void test_insertFunctionalGroupAllTemplatesConnectToTarget() {
    std::printf("--- Test: EVERY fg.sdf attachment-point template actually connects to its graft target ---\n");
    // Regression test for a real bug found this session: the verified root cause is that
    // several templates' "M  SAL" line had a fixed-column atom-count field one character
    // short (MDL SGroup format requires it to be exactly 3 characters wide). Indigo's
    // Scanner::readIntFix tolerates the resulting trailing whitespace instead of erroring on
    // it -- instead every subsequent 3-wide atom-index field in that line gets read one
    // column left of where it should be,
    // and the loop's last readIntFix(3) over-runs past the line's newline. The SGroup
    // handler's closing skipLine() then consumes the entire FOLLOWING line -- which is
    // "M  SAP" -- instead of parsing it, so indigoIterateSGroupAttachmentPoints yields zero
    // items, templateAttachIdx never resolves, graft never happens, and the template silently
    // falls back to free-floating placement instead of connecting to the target atom. (The 15
    // brief-listed templates were also given a M SMT/M SAL/M SAP reorder for consistency with
    // the other 47 working templates; that reorder is real and worth keeping, but two of the
    // 15 -- Cbz and TBDPS -- already had 3-digit SAL counts and were never actually broken by
    // this bug, confirming the SAL field width, not M-block order, is the necessary-and-
    // sufficient cause.) This iterates every template fg.sdf claims has an attachment point
    // and proves each one actually grafts.
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    const QStringList allAttachPointTemplates = {
        QStringLiteral("Ac"), QStringLiteral("Bn"), QStringLiteral("Boc"), QStringLiteral("Bu"),
        QStringLiteral("Bz"), QStringLiteral("Cbz"), QStringLiteral("C2H5"), QStringLiteral("CCl3"),
        QStringLiteral("CF3"), QStringLiteral("CN"), QStringLiteral("CO2Et"), QStringLiteral("CO2H"),
        QStringLiteral("CO2Me"), QStringLiteral("CONH2"), QStringLiteral("CO2Pr"), QStringLiteral("CO2tBu"),
        QStringLiteral("Cp"), QStringLiteral("CPh3"), QStringLiteral("Cy"), QStringLiteral("Et"),
        QStringLiteral("FMOC"), QStringLiteral("iBu"), QStringLiteral("Indole"), QStringLiteral("iPr"),
        QStringLiteral("Me"), QStringLiteral("Mes"), QStringLiteral("Ms"), QStringLiteral("NCO"),
        QStringLiteral("NCS"), QStringLiteral("NHPh"), QStringLiteral("NO2"), QStringLiteral("OAc"),
        QStringLiteral("OCF3"), QStringLiteral("OCN"), QStringLiteral("OEt"), QStringLiteral("OMe"),
        QStringLiteral("Ph"), QStringLiteral("PhCOOH"), QStringLiteral("Piv"), QStringLiteral("PO2"),
        QStringLiteral("PO3"), QStringLiteral("PO3H2"), QStringLiteral("PO4"), QStringLiteral("PO4H2"),
        QStringLiteral("Pr"), QStringLiteral("sBu"), QStringLiteral("SCN"), QStringLiteral("SO2"),
        QStringLiteral("SO2Cl"), QStringLiteral("SO2H"), QStringLiteral("SO3"), QStringLiteral("SO3H"),
        QStringLiteral("SO4"), QStringLiteral("SO4H"), QStringLiteral("ster"), QStringLiteral("TBDMS"),
        QStringLiteral("TBDPS"), QStringLiteral("tBu"), QStringLiteral("Tf"), QStringLiteral("TMS"),
        QStringLiteral("Tos"), QStringLiteral("Ts")
    };
    CHECK(allAttachPointTemplates.size() == 62, "test fixture lists all 62 known attachment-point templates");

    int failCount = 0;
    QStringList failedNames;
    for (const QString& name : allAttachPointTemplates) {
        DocumentState doc;
        AtomId target = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        doc.insertFunctionalGroup(lib, name, 0.0, 0.0, target, true);

        QList<AtomId> allAtoms = doc.molecule().atomIds();
        bool connected = allAtoms.size() > 1; // at minimum, something was inserted
        if (connected) {
            int targetFragment = doc.molecule().atomFragmentIndex(target);
            for (AtomId id : allAtoms) {
                if (id == target) continue;
                if (doc.molecule().atomFragmentIndex(id) != targetFragment) { connected = false; break; }
            }
        }
        if (!connected) { ++failCount; failedNames.append(name); }
    }
    if (failCount > 0) {
        std::printf("[FAIL DETAIL] templates that did not connect to their graft target: %s\n",
                    failedNames.join(QStringLiteral(", ")).toUtf8().constData());
    }
    CHECK(failCount == 0, "every fg.sdf attachment-point template grafts into the SAME fragment as its target atom");
}

static void test_insertFunctionalGroupConh2OneNeighborConnects() {
    std::printf("--- Test: CONH2 graft onto a 1-neighbor target connects (regression pin for the M SAL field-width bug) ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
    doc.molecule().addBond(a1, a2, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("CONH2"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 5, "graft added exactly three new atoms (CONH2's C, O, N)");

    bool a1GainedNewBond = false;
    for (BondId id : doc.molecule().bondIds()) {
        AtomId e1 = -1, e2 = -1;
        doc.molecule().bondEndpoints(id, e1, e2);
        if ((e1 == a1 || e2 == a1) && e1 != a2 && e2 != a2) a1GainedNewBond = true;
    }
    CHECK(a1GainedNewBond, "target atom a1 gained a new bond to the grafted CONH2 group (not just a2's pre-existing bond)");

    int targetFragment = doc.molecule().atomFragmentIndex(a1);
    bool allConnected = true;
    for (AtomId id : allAtoms) {
        if (doc.molecule().atomFragmentIndex(id) != targetFragment) { allConnected = false; break; }
    }
    CHECK(allConnected, "all 5 atoms (original chain + grafted CONH2) belong to the same fragment -- no disconnected group");
}

static void test_getFunctionalGroupPreviewFreeFloating() {
    std::printf("--- Test: getFunctionalGroupPreview matches insertFunctionalGroup's free-floating placement, without mutating the document ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    PlacementResult preview = doc.getFunctionalGroupPreview(lib, QStringLiteral("Ac"), 10.0, 10.0, -1);
    CHECK(preview.valid, "preview is valid for a real template ('Ac')");
    CHECK(preview.atoms.size() == 3, "preview has 3 atoms, matching 'Ac's real atom count");
    CHECK(preview.bonds.size() == 2, "preview has 2 bonds, matching 'Ac's real bond count");

    CHECK(doc.molecule().atomIds().isEmpty(), "the real document has zero atoms after computing a preview -- nothing was mutated");
    CHECK(!doc.canUndo(), "no undo history was pushed by computing a preview");

    // Commit the same call and confirm the committed atoms land at the same positions the
    // preview predicted (both centered on (10,10), same relative shape).
    doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 10.0, 10.0);
    CHECK(doc.molecule().atomIds().size() == 3, "commit creates 3 atoms, matching the preview's own atom count");

    double previewMinX = 1e9, previewMaxX = -1e9;
    for (const PreviewAtom& a : preview.atoms) {
        previewMinX = std::min(previewMinX, a.pos.x());
        previewMaxX = std::max(previewMaxX, a.pos.x());
    }
    double committedMinX = 1e9, committedMaxX = -1e9;
    for (AtomId id : doc.molecule().atomIds()) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        committedMinX = std::min(committedMinX, x);
        committedMaxX = std::max(committedMaxX, x);
    }
    CHECK(std::abs(previewMinX - committedMinX) < 0.01 && std::abs(previewMaxX - committedMaxX) < 0.01,
          "preview's predicted X-extent matches the actually-committed structure's X-extent exactly");
}

static void test_getFunctionalGroupPreviewGraft() {
    std::printf("--- Test: getFunctionalGroupPreview matches insertFunctionalGroup's graft-sprout placement, without mutating the document ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
    doc.molecule().addBond(a1, a2, 1);

    PlacementResult preview = doc.getFunctionalGroupPreview(lib, QStringLiteral("Ac"), 0.0, 0.0, a1);
    CHECK(preview.valid, "preview is valid for a graft onto a real target atom");
    CHECK(preview.atoms.size() == 3, "preview shows all 3 of 'Ac's atoms");

    CHECK(doc.molecule().atomIds().size() == 2, "the real document still has only its original 2 atoms -- preview did not mutate it");
    CHECK(!doc.canUndo(), "no undo history was pushed by computing a preview");

    // Under sprout semantics no preview atom sits on top of a1 anymore -- the attach atom is
    // offset by one bond length. Instead, the preview's own new graft-bond ghost segment should
    // connect a1's real position to the attach atom's predicted position.
    bool foundGraftBond = false;
    for (const PreviewBond& pb : preview.bonds) {
        bool startsAtTarget = std::abs(pb.start.x()) < 0.01 && std::abs(pb.start.y()) < 0.01;
        bool endsAtTarget = std::abs(pb.end.x()) < 0.01 && std::abs(pb.end.y()) < 0.01;
        if (startsAtTarget || endsAtTarget) foundGraftBond = true;
    }
    CHECK(foundGraftBond, "preview includes the not-yet-committed graft bond connecting a1's real position to the attach atom");

    doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 0.0, 0.0, a1, true);
    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 5, "commit: 2 pre-existing + 3 new (none fused) = 5 total");

    // Every committed new atom's position should appear among the preview's 3 atom positions --
    // this now holds for ALL 3 (including the attach atom, which no longer gets fused away).
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2) continue;
        double cx = 0, cy = 0;
        doc.molecule().atomPos(id, cx, cy);
        bool foundMatch = false;
        for (const PreviewAtom& pa : preview.atoms) {
            if (std::abs(pa.pos.x() - cx) < 0.01 && std::abs(pa.pos.y() - cy) < 0.01) foundMatch = true;
        }
        CHECK(foundMatch, "every committed new atom's position exactly matches one of the preview's predicted positions");
    }
}


static void test_documentStateSelection() {
    std::printf("--- Test 4: DocumentState selection ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = mol.addAtom(QStringLiteral("O"), 1, 0);
    BondId b1 = mol.addBond(a1, a2, 1);
    int ra1 = mol.addRxnArrow(0, 0, 5, 0);
    int rp1 = mol.addRxnPlus(2, 0);
    int mta1 = mol.addMultitailArrow({5,5, 0,4, 0,6});

    CHECK(doc.selection().isEmpty(), "fresh DocumentState has empty selection");

    // selectX REPLACES the whole selection (matches every real _selection
    // literal in 10-state.js -- each clears the other 4 sets).
    doc.selectAtom(a1);
    CHECK(doc.selection().atoms.contains(a1) && doc.selection().atoms.size() == 1,
          "selectAtom selects exactly that atom");
    CHECK(doc.selection().bonds.isEmpty(), "selectAtom clears bonds");

    doc.selectBond(b1);
    CHECK(doc.selection().bonds.contains(b1), "selectBond selects that bond");
    CHECK(doc.selection().atoms.isEmpty(), "selectBond clears atoms (replace, not add)");

    doc.selectRxnArrow(ra1);
    CHECK(doc.selection().rxnArrows.contains(ra1), "selectRxnArrow works");
    doc.selectRxnPlus(rp1);
    CHECK(doc.selection().rxnPluses.contains(rp1), "selectRxnPlus works");
    doc.selectMultitailArrow(mta1);
    CHECK(doc.selection().multitailArrows.contains(mta1), "selectMultitailArrow works");

    doc.clearSelection();
    CHECK(doc.selection().isEmpty(), "clearSelection empties everything");

    // add*ToSelection: additive, does not clear other types.
    doc.addAtomToSelection(a1);
    doc.addAtomToSelection(a2);
    doc.addBondToSelection(b1);
    CHECK(doc.selection().atoms.size() == 2, "addAtomToSelection is additive");
    CHECK(doc.selection().bonds.contains(b1), "addBondToSelection didn't clear atoms");

    doc.removeAtomFromSelection(a1);
    CHECK(!doc.selection().atoms.contains(a1) && doc.selection().atoms.contains(a2),
          "removeAtomFromSelection removes only the specified atom");
    CHECK(doc.selection().bonds.contains(b1), "removeAtomFromSelection didn't touch bonds");

    doc.clearSelection();
    doc.selectAll();
    CHECK(doc.selection().atoms.size() == 2, "selectAll selects both atoms");
    CHECK(doc.selection().bonds.size() == 1, "selectAll selects the bond");
    CHECK(doc.selection().rxnArrows.size() == 1, "selectAll selects the rxn arrow");
    CHECK(doc.selection().rxnPluses.size() == 1, "selectAll selects the rxn plus");
    CHECK(doc.selection().multitailArrows.size() == 1, "selectAll selects the multitail arrow");

    // Spec requirement: selection clears after both undo() and redo() (the
    // real 10-state.js undo()/redo() both reset _selection to all-empty).
    // undo()/redo() were implemented in Task 5 before selection() existed to
    // assert against -- verified here, now that it does.
    EditCommand addSulfur;
    addSulfur.execute = [&mol]() { mol.addAtom(QStringLiteral("S"), 3, 0); };
    addSulfur.invert = [&mol]() {
        QList<AtomId> ids = mol.atomIds();
        if (!ids.isEmpty()) mol.removeAtom(ids.last());
    };
    doc.executeCommand(addSulfur);
    doc.selectAll();
    CHECK(!doc.selection().isEmpty(), "selection is non-empty right before undo");
    doc.undo();
    CHECK(doc.selection().isEmpty(), "undo() clears the selection");

    doc.selectAll();
    CHECK(!doc.selection().isEmpty(), "selection is non-empty right before redo");
    doc.redo();
    CHECK(doc.selection().isEmpty(), "redo() clears the selection");
}

static void test_insertLibraryTemplateFused() {
    std::printf("--- Test 15: insertLibraryTemplateFused ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    // Real fusion success: "Indole" (9 atoms, 10 bonds, bondIdx=6) fused onto
    // a document bond. The seam bond's committed order must survive.
    {
        DocumentState doc;
        AtomId a0 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
        BondId seam = doc.molecule().addBond(a0, a1, 2);
        CHECK(!doc.canUndo(), "history starts empty (seeded via the molecule directly)");

        doc.insertLibraryTemplateFused(lib, QStringLiteral("Indole"), 0.75, 3.0, seam);
        CHECK(doc.canUndo(), "insertLibraryTemplateFused pushes exactly one history entry");
        CHECK(doc.molecule().bondOrder(seam) == 2, "the seam bond's committed order is UNCHANGED");
        // "Indole" has 9 atoms; a successful bond-fusion merges exactly 2 of
        // them onto the pre-existing a0/a1 (the fusion bond's own two atoms),
        // so 2 pre-existing + 9 new - 2 fused = 9 total.
        CHECK(doc.molecule().atomCount() == 9, "2 pre-existing + 9 new - 2 fused onto the seam = 9 total");
        doc.undo();
        CHECK(doc.molecule().atomCount() == 2, "undo restores to just a0/a1");
        CHECK(doc.molecule().bondCount() == 1, "undo restores to just the seam bond");
        CHECK(doc.molecule().bondOrder(seam) == 2, "the seam bond survives undo with its original order");
    }

    // Fallback: template with no fusion metadata (any fg.sdf entry).
    {
        DocumentState doc;
        AtomId a0 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
        BondId b = doc.molecule().addBond(a0, a1, 1);
        doc.insertLibraryTemplateFused(lib, QStringLiteral("Ac"), 0.0, 0.0, b);
        // Falls back to plain insertFunctionalGroup(lib, "Ac", cx, cy) --
        // "Ac" has 3 atoms, no graft (fallback never passes a target atom).
        CHECK(doc.molecule().atomCount() == 2 + 3, "no-fusion-metadata fallback pastes \"Ac\" plainly");
    }

    // Fallback: targetBondId that doesn't resolve.
    {
        DocumentState doc;
        doc.insertLibraryTemplateFused(lib, QStringLiteral("Indole"), 0.0, 0.0, 9999);
        CHECK(doc.molecule().atomCount() == 9, "missing target bond: falls back to plain \"Indole\" paste");
    }
}

static void test_toggleSgroupExpanded() {
    std::printf("--- Test 16: toggleSgroupExpanded ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    // Seeded via doc.molecule() directly (not doc.insertFunctionalGroup,
    // which is itself undoable) so history starts empty and canUndo() is an
    // exact, unambiguous signal for the checks below -- the same seeding
    // discipline used throughout this migration whenever a test needs to
    // prove "this call pushed exactly one/zero history entries."
    DocumentState doc;
    QString acMolfile = lib.molfileText(lib.functionalGroup(QStringLiteral("Ac")));
    EditableMolecule::InsertResult r = doc.molecule().insertStructure(
        acMolfile, [](double x, double y) { return QPointF(x, y); });
    CHECK(r.createdSGroups.size() == 1, "setup: \"Ac\" produces exactly one sgroup");
    SGroupId sgId = r.createdSGroups[0];
    CHECK(doc.molecule().sgroupExpanded(sgId), "newly-inserted sgroup starts expanded=true");
    CHECK(!doc.canUndo(), "history is empty (sgroup seeded via the molecule directly)");

    doc.toggleSgroupExpanded(sgId);
    CHECK(doc.canUndo(), "toggleSgroupExpanded pushes exactly one history entry");
    CHECK(!doc.molecule().sgroupExpanded(sgId), "toggle flips expanded true -> false");
    doc.undo();
    CHECK(doc.molecule().sgroupExpanded(sgId), "undo restores expanded=true");
    CHECK(!doc.canUndo(), "history is empty again after undoing the only entry");

    doc.toggleSgroupExpanded(9999);
    CHECK(!doc.canUndo(), "invalid SGroupId pushes no history entry");
    CHECK(doc.molecule().sgroupExpanded(sgId), "invalid-id toggle call left the real sgroup untouched");
}

static void test_rxnArrowLifecycle() {
    std::printf("--- Test 17: rxn-arrow lifecycle ---\n");
    DocumentState doc;

    RxnArrowId a1 = doc.addRxnArrow(0, 0);
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    CHECK(doc.molecule().rxnArrowEndpoints(a1, x1, y1, x2, y2), "addRxnArrow creates a real arrow");
    CHECK(x1 == 0 && y1 == 0 && x2 == 1.5 * 2.5 && y2 == 0, "default second endpoint is 2.5 bond-lengths to the right");
    CHECK(doc.molecule().rxnArrowMode(a1) == QStringLiteral("filled-triangle"), "default mode is filled-triangle");
    doc.undo();
    CHECK(doc.molecule().rxnArrowIds().isEmpty(), "undo removes the created arrow");

    // Recreate a1 for real -- the throwaway above only proved undo works; everything below
    // needs a persistent first arrow to compare counts against.
    a1 = doc.addRxnArrow(0, 0);

    RxnArrowId a2 = doc.addRxnArrow(1, 1, QStringLiteral("resonance"));
    CHECK(doc.molecule().rxnArrowMode(a2) == QStringLiteral("resonance"), "explicit mode is applied at creation");

    RxnArrowId curved = doc.addCurvedArrow(0, 0, 1, 2, 3, 0);
    CHECK(doc.molecule().rxnArrowMode(curved) == QStringLiteral("curved-mechanism"), "addCurvedArrow sets curved-mechanism mode");
    double cx = 0, cy = 0;
    CHECK(doc.molecule().rxnArrowCurvature(curved, cx, cy) && cx == 1 && cy == 2, "addCurvedArrow stores the control point");
    doc.undo();
    CHECK(!doc.molecule().rxnArrowIds().contains(curved), "undo removes the curved arrow too");

    doc.setRxnArrowMode(a2, QStringLiteral("hollow-triangle"));
    CHECK(doc.molecule().rxnArrowMode(a2) == QStringLiteral("hollow-triangle"), "setRxnArrowMode changes the mode");
    bool canUndoBefore = doc.canUndo();
    doc.setRxnArrowMode(a2, QStringLiteral("hollow-triangle")); // unchanged -- must no-op
    CHECK(doc.canUndo() == canUndoBefore, "setRxnArrowMode to the SAME mode pushes no history entry");
    doc.undo();
    CHECK(doc.molecule().rxnArrowMode(a2) == QStringLiteral("resonance"), "undo restores the original mode");

    doc.setRxnArrowConditions(a2, QStringLiteral("Pd/C"), QStringLiteral("H2"));
    CHECK(doc.molecule().rxnArrowConditionsAbove(a2) == QStringLiteral("Pd/C"), "conditions above set");
    CHECK(doc.molecule().rxnArrowConditionsBelow(a2) == QStringLiteral("H2"), "conditions below set");
    doc.undo();
    CHECK(doc.molecule().rxnArrowConditionsAbove(a2).isEmpty(), "undo restores conditions to empty");

    doc.setRxnArrowConditions(a2, QStringLiteral("cond"), QStringLiteral("cond2"));
    doc.molecule().setRxnArrowCurvature(a2, 9, 9, true); // no DocumentState-level command for this exists (matches the real JS: curvature is only ever set via addCurvedArrow) -- mutate the molecule directly so the deleteRxnArrow capture-before-delete below has real curvature to preserve
    doc.deleteRxnArrow(a2);
    CHECK(doc.molecule().rxnArrowIds().size() == 1, "deleteRxnArrow removes the arrow (only a1 remains)");
    doc.undo();
    CHECK(doc.molecule().rxnArrowIds().size() == 2, "undo recreates the deleted arrow (fresh id)");
    bool foundRestored = false;
    for (RxnArrowId id : doc.molecule().rxnArrowIds()) {
        if (id != a1 && doc.molecule().rxnArrowConditionsAbove(id) == QStringLiteral("cond")) foundRestored = true;
    }
    CHECK(foundRestored, "restored arrow keeps its mode/conditions/curvature from before deletion");

    doc.deleteRxnArrow(9999);
    CHECK(!doc.canUndo() || doc.molecule().rxnArrowIds().size() == 2, "deleteRxnArrow on an unknown id is a no-op");
}

static void test_rxnPlusCommands() {
    std::printf("--- Test 23: rxn-plus add/delete ---\n");
    DocumentState doc;

    RxnPlusId p1 = doc.addRxnPlus(2, 3);
    double x = 0, y = 0;
    CHECK(doc.molecule().rxnPlusPos(p1, x, y) && x == 2 && y == 3, "addRxnPlus creates a real rxn-plus at the given position");
    doc.undo();
    CHECK(doc.molecule().rxnPlusIds().isEmpty(), "undo removes the created rxn-plus");

    RxnPlusId p2 = doc.addRxnPlus(5, 5);
    doc.deleteRxnPlus(p2);
    CHECK(doc.molecule().rxnPlusIds().isEmpty(), "deleteRxnPlus removes it");
    doc.undo();
    CHECK(doc.molecule().rxnPlusIds().size() == 1, "undo recreates it (fresh id)");
    RxnPlusId restored = doc.molecule().rxnPlusIds().first();
    CHECK(doc.molecule().rxnPlusPos(restored, x, y) && x == 5 && y == 5, "restored rxn-plus keeps its original position");

    doc.deleteRxnPlus(9999);
    CHECK(!doc.canUndo() || doc.molecule().rxnPlusIds().size() == 1, "deleteRxnPlus on an unknown id is a no-op");
}

static void test_multitailArrowCommands() {
    std::printf("--- Test 24: multitail arrow add/delete/add-tail ---\n");
    DocumentState doc;

    MultitailArrowId id = doc.addMultitailArrow(0, 0);
    QList<double> pts = doc.molecule().multitailArrowPoints(id);
    CHECK(pts.size() == 2, "addMultitailArrow starts with head only, zero tails -- matches the real tailsYOffset=[]");
    double headX = pts[0], headY = pts[1];
    CHECK(headX == 2.0 * 1.5 && headY == 0.0, "head is offset (2*bondLength, 0) from the creation point");
    doc.undo();
    CHECK(doc.molecule().multitailArrowIds().isEmpty(), "undo removes the created arrow");

    MultitailArrowId id2 = doc.addMultitailArrow(0, 0);
    doc.addMultitailArrowTail(id2);
    QList<double> pts1 = doc.molecule().multitailArrowPoints(id2);
    CHECK(pts1.size() == 4, "one addMultitailArrowTail call adds exactly one tail (2 more numbers)");
    double firstTailY = pts1[3];
    CHECK(firstTailY == (0.0 + 3.0 * 1.5) / 2.0, "first tail is centered at the midpoint of the full spine (widest gap)");
    doc.undo();
    CHECK(doc.molecule().multitailArrowPoints(id2).size() == 2, "undo removes the added tail");

    doc.addMultitailArrowTail(id2);
    doc.addMultitailArrowTail(id2);
    QList<double> pts2 = doc.molecule().multitailArrowPoints(id2);
    CHECK(pts2.size() == 6, "two addMultitailArrowTail calls produce two tails");
    CHECK(pts2[3] != pts2[5], "the two tails land at different Y positions (widest-gap search, not a fixed spot)");
    doc.undo();
    CHECK(doc.molecule().multitailArrowPoints(id2).size() == 4, "undo removes only the SECOND tail");

    doc.addMultitailArrowTail(9999);
    CHECK(!doc.canUndo() || doc.molecule().multitailArrowPoints(id2).size() == 4,
          "addMultitailArrowTail on an unknown id is a no-op");

    doc.deleteMultitailArrow(id2);
    CHECK(doc.molecule().multitailArrowIds().isEmpty(), "deleteMultitailArrow removes it");
    doc.undo();
    CHECK(doc.molecule().multitailArrowIds().size() == 1, "undo recreates it (fresh id) with the same points");
    MultitailArrowId restored = doc.molecule().multitailArrowIds().first();
    CHECK(doc.molecule().multitailArrowPoints(restored).size() == 4, "restored arrow keeps its one tail");

    doc.deleteMultitailArrow(9999);
    CHECK(!doc.canUndo() || doc.molecule().multitailArrowIds().size() == 1, "deleteMultitailArrow on an unknown id is a no-op");
}

static void test_setStereoFlags() {
    std::printf("--- Test 18: setStereoFlags ---\n");
    DocumentState doc;
    CHECK(doc.molecule().stereoFlagsType() == QStringLiteral("abs"), "default type is abs");
    CHECK(doc.molecule().stereoFlagsGroupId() == 0, "default groupId is 0");

    doc.setStereoFlags(QStringLiteral("rel"), 2);
    CHECK(doc.molecule().stereoFlagsType() == QStringLiteral("rel"), "setStereoFlags sets type");
    CHECK(doc.molecule().stereoFlagsGroupId() == 2, "setStereoFlags sets groupId");
    bool canUndoBefore = doc.canUndo();
    doc.setStereoFlags(QStringLiteral("rel"), 2); // unchanged -- must no-op
    CHECK(doc.canUndo() == canUndoBefore, "unchanged setStereoFlags pushes no history entry");
    doc.undo();
    CHECK(doc.molecule().stereoFlagsType() == QStringLiteral("abs"), "undo restores type to abs");
    CHECK(doc.molecule().stereoFlagsGroupId() == 0, "undo restores groupId to 0");

    // Falsy args default like the real JS (type || 'abs', groupId || 0).
    doc.setStereoFlags(QString(), 0);
    canUndoBefore = doc.canUndo();
    doc.setStereoFlags(QString(), 0);
    CHECK(doc.canUndo() == canUndoBefore, "empty-string type still normalizes to the current abs/0 state -- no-op");
}

static void test_rgroupCommands() {
    std::printf("--- Test 19: R-group commands ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 5, 5);   // separate fragment (no bond to a1)

    CHECK(doc.addRGroup(1), "addRGroup(1) succeeds");
    CHECK(!doc.addRGroup(1), "addRGroup(1) again fails");
    doc.undo();
    CHECK(doc.molecule().rgroupNumbers().isEmpty(), "undo removes the R-group");
    doc.redo();
    CHECK(doc.molecule().rgroupNumbers() == QList<int>{1}, "redo recreates it");

    doc.setRGroupLogic(1, QStringLiteral("1,2"), true, 1);
    QString range; bool resth = false; int ifthen = 0;
    doc.molecule().rgroupLogic(1, range, resth, ifthen);
    CHECK(range == QStringLiteral("1,2") && resth && ifthen == 1, "setRGroupLogic applied");
    doc.undo();
    doc.molecule().rgroupLogic(1, range, resth, ifthen);
    CHECK(range.isEmpty() && !resth && ifthen == 0, "undo restores default logic fields");
    doc.setRGroupLogic(1, QStringLiteral("1,2"), true, 1); // reapply so the rest of this test has it set

    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);
    doc.addRGroupMember(1);
    QList<int> members = doc.molecule().rgroupFragmentIds(1);
    CHECK(members.size() == 2, "addRGroupMember registers both selected atoms' (distinct) fragments");
    doc.undo();
    CHECK(doc.molecule().rgroupFragmentIds(1).isEmpty(), "undo removes both member fragments");
    doc.redo();
    CHECK(doc.molecule().rgroupFragmentIds(1).size() == 2, "redo restores both");

    bool canUndoBefore = doc.canUndo();
    doc.addRGroupMember(1); // same selection again -- no NEW fragments -- must no-op
    CHECK(doc.canUndo() == canUndoBefore, "addRGroupMember with nothing new to add pushes no history entry");

    int fragToRemove = members.first();
    doc.removeRGroupMember(1, fragToRemove);
    CHECK(doc.molecule().rgroupFragmentIds(1).size() == 1, "removeRGroupMember removes one fragment");
    doc.undo();
    CHECK(doc.molecule().rgroupFragmentIds(1).size() == 2, "undo restores it");

    doc.removeRGroupMember(1, 9999);
    CHECK(!doc.canUndo() || doc.molecule().rgroupFragmentIds(1).size() == 2,
          "removeRGroupMember with a non-member fragId is a no-op");

    CHECK(doc.deleteRGroup(1), "deleteRGroup succeeds");
    CHECK(doc.molecule().rgroupNumbers().isEmpty(), "R-group gone");
    doc.undo();
    CHECK(doc.molecule().rgroupNumbers() == QList<int>{1}, "undo recreates the R-group");
    doc.molecule().rgroupLogic(1, range, resth, ifthen);
    CHECK(range == QStringLiteral("1,2") && resth && ifthen == 1, "undo restores its logic fields too");
    CHECK(doc.molecule().rgroupFragmentIds(1).size() == 2, "undo restores its member fragments too");

    CHECK(!doc.deleteRGroup(9999), "deleteRGroup fails for an unknown R-group number");
}

static void test_textCommands() {
    std::printf("--- Test 20: text commands ---\n");
    DocumentState doc;

    TextId t1 = doc.addText(QStringLiteral("hello"), 1, 2, false, false);
    double x = 0, y = 0; QString content; bool bold = false, italic = false;
    CHECK(doc.molecule().textAnnotationContent(t1, x, y, content, bold, italic), "addText creates a real annotation");
    CHECK(content == QStringLiteral("hello") && x == 1 && y == 2, "content and position match");
    doc.undo();
    CHECK(doc.molecule().textAnnotationIds().isEmpty(), "undo removes the created text");

    // Recreate t1 for real -- the throwaway above only proved undo works; everything below needs
    // a persistent first annotation to compare counts against.
    t1 = doc.addText(QStringLiteral("hello"), 1, 2, false, false);
    bool canUndoBeforeNoOps = doc.canUndo();

    CHECK(doc.addText(QString(), 0, 0, false, false) == -1, "addText on an empty string returns -1 (no-op)");
    CHECK(doc.addText(QStringLiteral("   "), 0, 0, false, false) == -1, "addText on a whitespace-only string returns -1 (no-op)");
    CHECK(doc.canUndo() == canUndoBeforeNoOps, "neither no-op addText call pushed a history entry");

    TextId t2 = doc.addText(QStringLiteral("world"), 0, 0, false, false);
    doc.updateText(t2, QStringLiteral("WORLD"), true, true);
    CHECK(doc.molecule().textAnnotationContent(t2, x, y, content, bold, italic), "re-read after updateText");
    CHECK(content == QStringLiteral("WORLD") && bold && italic, "updateText changed content/bold/italic");
    doc.undo();
    CHECK(doc.molecule().textAnnotationContent(t2, x, y, content, bold, italic), "re-read after undo");
    CHECK(content == QStringLiteral("world") && !bold && !italic, "undo restores the original content/bold/italic");

    doc.deleteText(t2);
    CHECK(doc.molecule().textAnnotationIds().size() == 1, "deleteText removes the annotation");
    doc.undo();
    CHECK(doc.molecule().textAnnotationIds().size() == 2, "undo recreates it (fresh id)");

    doc.deleteText(9999);
    CHECK(!doc.canUndo() || doc.molecule().textAnnotationIds().size() == 2, "deleteText on an unknown id is a no-op");
}

static void test_addBracketSelection() {
    std::printf("--- Test 21: addBracketSelection ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 4, 3);
    doc.addBond(a1, a2, 1);

    bool canUndoBeforeEmpty = doc.canUndo();
    doc.addBracketSelection(); // nothing selected -- no-op
    CHECK(doc.canUndo() == canUndoBeforeEmpty, "addBracketSelection with an empty selection is a no-op");

    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);
    doc.addBracketSelection();
    CHECK(doc.molecule().bracketCount() == 1, "addBracketSelection pushes one bracket");
    double minX = 0, minY = 0, maxX = 0, maxY = 0;
    doc.molecule().bracketAt(0, minX, minY, maxX, maxY);
    CHECK(minX == 0 - 0.8 && minY == 0 - 0.8 && maxX == 4 + 0.8 && maxY == 3 + 0.8,
          "bracket bbox covers both atoms padded by 0.8, matching the real function");
    doc.undo();
    CHECK(doc.molecule().bracketCount() == 0, "undo pops the bracket");

    // Bond-only selection also contributes both endpoints' positions.
    doc.clearSelection();
    BondId b1 = doc.molecule().bondIds().first();
    doc.addBondToSelection(b1);
    doc.addBracketSelection();
    CHECK(doc.molecule().bracketCount() == 1, "a bond-only selection also produces a bracket");
    doc.molecule().bracketAt(0, minX, minY, maxX, maxY);
    CHECK(maxX == 4 + 0.8, "bond's far endpoint contributes to the bbox");
}

static void test_imageCommands() {
    std::printf("--- Test 22: image add/delete ---\n");
    DocumentState doc;

    ImageId img1 = doc.addImage(QByteArray("fakepng"), 1, 2, 4, 3);
    double x = 0, y = 0, w = 0, h = 0; QByteArray png;
    CHECK(doc.molecule().imageData(img1, x, y, w, h, png), "addImage creates a real image");
    CHECK(x == 1 && y == 2 && w == 4 && h == 3, "center/half-extents match the arguments");
    CHECK(png == QByteArray("fakepng"), "image data round-trips");
    doc.undo();
    CHECK(doc.molecule().imageIds().isEmpty(), "undo removes the created image");

    CHECK(doc.addImage(QByteArray(), 0, 0, 1, 1) == -1, "addImage with empty data returns -1 (no-op)");
    CHECK(!doc.canUndo(), "the no-op addImage call pushed no history entry");

    ImageId img2 = doc.addImage(QByteArray("otherpng"), 0, 0, 1, 1);
    doc.deleteImage(img2);
    CHECK(doc.molecule().imageIds().isEmpty(), "deleteImage removes the image");
    doc.undo();
    CHECK(doc.molecule().imageIds().size() == 1, "undo recreates it (fresh id)");
    ImageId restored = doc.molecule().imageIds().first();
    CHECK(doc.molecule().imageData(restored, x, y, w, h, png) && png == QByteArray("otherpng"),
          "restored image keeps its original data/geometry");

    doc.deleteImage(9999);
    CHECK(!doc.canUndo() || doc.molecule().imageIds().size() == 1, "deleteImage on an unknown id is a no-op");
}

static void test_deserializeMol() {
    std::printf("--- Test 25: deserializeMol ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    doc.molecule().setName(QStringLiteral("original"));
    doc.selectAtom(a1);
    CHECK(!doc.selection().atoms.isEmpty(), "setup: something is selected before deserializeMol");

    doc.deserializeMol(QStringLiteral("c1ccccc1"));
    CHECK(doc.molecule().atomCount() == 6, "deserializeMol replaces the document with benzene's 6 atoms");
    CHECK(doc.molecule().name().isEmpty(), "the old name is gone (extension data wiped, matching the real JS)");
    CHECK(doc.selection().atoms.isEmpty(), "selection is cleared by deserializeMol itself, not just by a later undo/redo");

    doc.undo();
    CHECK(doc.molecule().atomCount() == 1, "undo restores the original 1-atom document");
    CHECK(doc.molecule().name() == QStringLiteral("original"), "undo restores the original name too");

    doc.redo();
    CHECK(doc.molecule().atomCount() == 6, "redo re-applies the load");

    bool canUndoBefore = doc.canUndo();
    doc.deserializeMol(QStringLiteral("not a valid molecule at all $$$"));
    CHECK(doc.canUndo() == canUndoBefore, "a parse failure pushes no history entry (no-op)");
    CHECK(doc.molecule().atomCount() == 6, "a parse failure leaves the current document completely unchanged");
}

static void test_deserializeMolCenterOnPage() {
    std::printf("--- Test: deserializeMol centerOnPage ---\n");

    QString offOriginMolfile =
        "\n"
        "  -INDIGO-\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "   98.0000   98.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  102.0000  102.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n";

    DocumentState doc;
    doc.deserializeMol(offOriginMolfile);
    CHECK(doc.molecule().atomIds().size() == 2, "setup: 2 atoms loaded");
    {
        QList<AtomId> ids = doc.molecule().atomIds();
        double x = 0, y = 0;
        CHECK(doc.molecule().atomPos(ids[0], x, y) && std::abs(x - 98.0) < 0.001,
              "centerOnPage=false (default) leaves coordinates untouched");
    }

    doc.deserializeMol(offOriginMolfile, true);
    {
        QList<AtomId> ids = doc.molecule().atomIds();
        CHECK(ids.size() == 2, "2 atoms after centered reload");
        double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
        doc.molecule().atomPos(ids[0], x1, y1);
        doc.molecule().atomPos(ids[1], x2, y2);
        double minX = std::min(x1, x2), maxX = std::max(x1, x2);
        double minY = std::min(y1, y2), maxY = std::max(y1, y2);
        CHECK(std::abs(minX + maxX) < 0.001, "centerOnPage=true: bbox center x is 0");
        CHECK(std::abs(minY + maxY) < 0.001, "centerOnPage=true: bbox center y is 0");
        CHECK(std::abs((maxX - minX) - 4.0) < 0.001, "shape preserved (atoms still 4 apart)");
    }

    DocumentState doc2;
    doc2.loadBenzene();
    int atomsBeforeLoad = doc2.molecule().atomIds().size();
    doc2.deserializeMol(offOriginMolfile, true);
    CHECK(doc2.molecule().atomIds().size() == 2, "centered load replaced benzene with 2 atoms");
    CHECK(doc2.canUndo(), "load pushed a history entry");
    doc2.undo();
    CHECK(doc2.molecule().atomIds().size() == atomsBeforeLoad,
          "single undo fully restores the pre-load document (one atomic step, not two)");

    DocumentState doc3;
    doc3.deserializeMol(QString(), true);
    CHECK(doc3.molecule().atomIds().size() == 0, "empty input with centerOnPage=true does not crash");
}

static const char* kTestTemplateOneBond =
    "TestFG\n"
    "Ketcher 11161713142D 1   1.00000     0.00000     0\n"
    "\n"
    "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "M  STY  1   1 SUP\n"
    "M  SLB  1   1   1\n"
    "M  SAL   1  2   1   2\n"
    "M  SAP   1  1   1   0\n"
    "M  SMT   1 TestFG\n"
    "M  END\n";

// Attach atom (F) has ZERO template-internal bonds -- disconnected from the template's other
// atom (C). This still exercises the "no rotation" fallback (templateAngle stays nullopt because
// the attach atom has no bonds), but unlike a genuinely single-atom template, the disconnected C
// atom is NOT bonded to the attach atom and so is NOT merged away by graftAtomOnto -- it survives
// as a real, checkable atom, so the test can assert an actual position instead of only an atom
// count (see test_graftAngleOrientation_noTemplateInternalBonds).
static const char* kTestTemplateDisconnected =
    "TestFG\n"
    "Ketcher 11161713142D 1   1.00000     0.00000     0\n"
    "\n"
    "  2  0  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "M  STY  1   1 SUP\n"
    "M  SLB  1   1   1\n"
    "M  SAL   1  2   1   2\n"
    "M  SAP   1  1   1   0\n"
    "M  SMT   1 TestFG\n"
    "M  END\n";

// Attach atom (N) has TWO template-internal bonds, to atoms at 0 degrees and 120 degrees from it
// -- the multi-bond-attach-atom case the "mean neighbor direction" fix addresses (see
// test_graftAngleOrientation_multiBondAttachAtom).
static const char* kTestTemplateTwoBondsAttach =
    "TestFG\n"
    "Ketcher 11161713142D 1   1.00000     0.00000     0\n"
    "\n"
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "M  STY  1   1 SUP\n"
    "M  SLB  1   1   1\n"
    "M  SAL   1  3   1   2   3\n"
    "M  SAP   1  1   1   0\n"
    "M  SMT   1 TestFG\n"
    "M  END\n";

static bool hasBond(DocumentState& doc, AtomId a, AtomId b) {
    for (BondId id : doc.molecule().bondIds()) {
        AtomId e1 = -1, e2 = -1;
        if (!doc.molecule().bondEndpoints(id, e1, e2)) continue;
        if ((e1 == a && e2 == b) || (e1 == b && e2 == a)) return true;
    }
    return false;
}

static void test_graftAngleOrientation_oneNeighbor() {
    std::printf("--- Test: insertFunctionalGroup graft rotates for 1-neighbor target ---\n");
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_1n"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
    doc.molecule().addBond(a1, a2, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    // Expect the same 60-degree kink BondAngleSuggester itself produces for a 1-neighbor atom:
    // continuation of the a2->a1 bond (a2 is at angle 0 from a1, so continuation past a1 points
    // along -x = PI) kinked by 60 degrees (PI/3). Under sprout semantics the attach atom (O)
    // survives as its own atom one bond length from a1, and the far atom (C) continues in the
    // same rotated direction, one more template-internal bond length out.
    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 4, "graft added exactly two new atoms (attach atom survives, no longer merged into a1)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("O")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    CHECK(hasBond(doc, a1, attachAtom), "target atom a1 is bonded directly to the attach atom");
    CHECK(hasBond(doc, attachAtom, farAtom), "attach atom keeps its own original template-internal bond to the far atom");
    CHECK(!hasBond(doc, a1, farAtom), "target atom a1 is NOT directly bonded to the far atom (no fusion, just one new bond)");

    double ax = 0, ay = 0, fx = 0, fy = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    double expected = M_PI + M_PI / 3.0; // continuation (PI) + 60-degree kink, matching
                                          // BondAngleSuggester's own default turn direction
    double attachAngle = std::atan2(ay, ax);
    double farAngle = std::atan2(fy, fx);
    while (attachAngle < 0) attachAngle += 2 * M_PI;
    while (farAngle < 0) farAngle += 2 * M_PI;
    while (expected >= 2 * M_PI) expected -= 2 * M_PI;
    CHECK(std::abs(attachAngle - expected) < 0.01, "attach atom lands at the 60-degree-kink angle, not the template's raw 0-degree direction");
    CHECK(std::abs(farAngle - expected) < 0.01, "far atom lands along the same rotated direction as the attach atom");
    CHECK(std::abs(std::hypot(ax, ay) - 1.5) < 0.01, "attach atom sits exactly one bond length (kBondLength=1.5) from the target atom");
    CHECK(std::abs(std::hypot(fx, fy) - 2.5) < 0.01, "far atom sits one template-internal bond length further out from the attach atom");
}

static void test_graftAngleOrientation_twoNeighbors() {
    std::printf("--- Test: insertFunctionalGroup graft rotates for 2-neighbor target (120-degree bisector) ---\n");
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_2n"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.0, 0.0);   // angle 0 from a1
    AtomId a3 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 1.0);   // angle PI/2 from a1
    doc.molecule().addBond(a1, a2, 1);
    doc.molecule().addBond(a1, a3, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 5, "graft added exactly two new atoms (attach atom survives)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2 || id == a3) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("O")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    CHECK(hasBond(doc, a1, attachAtom), "target atom a1 is bonded directly to the attach atom");
    CHECK(!hasBond(doc, a1, farAtom), "target atom a1 is NOT directly bonded to the far atom");

    double ax = 0, ay = 0, fx = 0, fy = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    double attachAngle = std::atan2(ay, ax);
    double farAngle = std::atan2(fy, fx);
    // Existing neighbors at 0 and PI/2: the gap going counterclockwise from 0 to PI/2 is
    // width PI/2; the gap going counterclockwise from PI/2 back around to 0 is width 3PI/2 --
    // the larger gap, so the bisector sits at PI/2 + 3PI/2/2 = PI/2 + 3PI/4 = 5PI/4.
    double expected = 5.0 * M_PI / 4.0;
    while (attachAngle < 0) attachAngle += 2 * M_PI;
    while (farAngle < 0) farAngle += 2 * M_PI;
    while (expected >= 2 * M_PI) expected -= 2 * M_PI;
    CHECK(std::abs(attachAngle - expected) < 0.01, "attach atom bisects the larger gap between the two existing neighbors");
    CHECK(std::abs(farAngle - expected) < 0.01, "far atom lands along the same rotated direction as the attach atom");
    CHECK(std::abs(std::hypot(ax, ay) - 1.5) < 0.01, "attach atom sits one bond length from the target atom");
    CHECK(std::abs(std::hypot(fx, fy) - 2.5) < 0.01, "far atom sits one more template-internal bond length further out");
}

static void test_graftAngleOrientation_ringAtomTarget() {
    std::printf("--- Test: insertFunctionalGroup graft on a RING atom points outward, away from ring center ---\n");
    // BondAngleSuggester::suggestAngle has no ring/chain distinction anywhere in its code (pure
    // neighbor-count + local-geometry bisection) -- every existing graft-angle test only ever
    // builds an open-chain fixture, so this specific topology (target atom is a ring member with
    // two ring-bond neighbors at a real ring's interior angle, not an arbitrary 90-degree chain
    // fixture) has never actually been exercised end-to-end through insertFunctionalGroup.
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_ring"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    // A regular hexagon (real ring geometry, 120-degree interior angle at every vertex),
    // vertices at 60-degree spacing around the origin, radius 1. Target the vertex at angle 0
    // (position (1,0)) -- its two ring neighbors are the adjacent vertices at +/-60 degrees.
    AtomId ring[6];
    for (int i = 0; i < 6; ++i) {
        double a = i * M_PI / 3.0;
        ring[i] = doc.molecule().addAtom(QStringLiteral("C"), std::cos(a), std::sin(a));
    }
    for (int i = 0; i < 6; ++i) doc.molecule().addBond(ring[i], ring[(i + 1) % 6], 1);
    AtomId target = ring[0]; // (1, 0) -- ring neighbors ring[1] (60 deg) and ring[5] (-60 deg)

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 2.0, 0.0, target, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 8, "graft added exactly two new atoms (attach atom survives, no longer fused into the ring vertex)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        bool isRing = false;
        for (int i = 0; i < 6; ++i) if (id == ring[i]) isRing = true;
        if (isRing) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("O")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    CHECK(hasBond(doc, target, attachAtom), "ring target atom is bonded directly to the attach atom");
    CHECK(hasBond(doc, attachAtom, farAtom), "attach atom keeps its own template-internal bond to the far atom");
    CHECK(!hasBond(doc, target, farAtom), "ring target atom is NOT directly bonded to the far atom");
    CHECK(!hasBond(doc, ring[1], attachAtom) && !hasBond(doc, ring[5], attachAtom),
          "attach atom is not fused into the ring -- it has no bond to the target's OWN ring neighbors");

    double ax = 0, ay = 0, fx = 0, fy = 0, tx = 0, ty = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    doc.molecule().atomPos(target, tx, ty);
    // Ring neighbors sit at +60 and -60 degrees from the target; the larger gap (240 degrees,
    // the arc swinging AWAY from the ring interior) is the one bisected, landing at exactly 0
    // degrees -- i.e. straight out along the same radial direction as the target vertex itself.
    double attachAngle = std::atan2(ay - ty, ax - tx);
    double farAngle = std::atan2(fy - ty, fx - tx);
    while (attachAngle < 0) attachAngle += 2 * M_PI;
    while (farAngle < 0) farAngle += 2 * M_PI;
    CHECK(std::abs(attachAngle) < 0.01 || std::abs(attachAngle - 2 * M_PI) < 0.01,
          "grafted attach atom points radially outward from the ring center, not inward");
    CHECK(std::abs(farAngle) < 0.01 || std::abs(farAngle - 2 * M_PI) < 0.01,
          "grafted far atom continues in the same outward direction");
    CHECK(std::abs(std::hypot(ax - tx, ay - ty) - 1.5) < 0.01, "attach atom sits one bond length from the ring target atom");
    CHECK(std::abs(std::hypot(fx - tx, fy - ty) - 2.5) < 0.01, "far atom sits one more template-internal bond length further out");
    double outwardDot = (ax - tx) * 1.0 + (ay - ty) * 0.0;
    CHECK(outwardDot > 0.5, "attach atom is on the outward side of the ring vertex, not the inward side");
}

static void test_graftAngleOrientation_noRotationCase() {
    std::printf("--- Test: insertFunctionalGroup graft on an isolated (0-neighbor) target uses the largest-empty-angle fallback ---\n");
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_0n"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);  // zero existing neighbors

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 3, "graft added exactly two new atoms (attach atom survives)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        if (id == a1) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("O")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    CHECK(hasBond(doc, a1, attachAtom), "target atom is bonded directly to the attach atom");

    double ax = 0, ay = 0, fx = 0, fy = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    // BondAngleSuggester::suggestAngle returns nullopt for a 0-neighbor target; the new
    // suggestFallbackAngle fallback returns exactly angle 0 for that case (matching
    // V8Process::getLargestEmptyAngle's own 0-neighbor branch), and the template's own
    // templateAngle is also 0 for this fixture, so rotateBy is 0 -- the template's raw
    // orientation happens to be preserved here, but now offset by a real bond length instead of
    // exact overlap.
    CHECK(std::abs(ax - 1.5) < 0.01 && std::abs(ay - 0.0) < 0.01, "attach atom lands one bond length (1.5) along the fallback angle (0 degrees)");
    CHECK(std::abs(fx - 2.5) < 0.01 && std::abs(fy - 0.0) < 0.01, "far atom continues in the same direction, one more template-internal bond length out");
}

static void test_graftAngleOrientation_noTemplateInternalBonds() {
    std::printf("--- Test: insertFunctionalGroup graft with a bond-less attach atom skips rotation but still offsets by one bond length ---\n");
    // Attach atom (F) has zero template-internal bonds -- templateAngle stays nullopt, so
    // rotation must not apply even though a suggested angle exists for the target (this
    // invariant is unchanged by the sprout-not-merge change). The disconnected second atom (C)
    // is translated along with the attach atom (same dx,dy), preserving their raw relative
    // template-space offset unrotated.
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateDisconnected), QStringLiteral("test_graft_nobonds"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
    doc.molecule().addBond(a1, a2, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 4, "graft added exactly two new atoms (attach atom F now survives instead of merging into a1)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("F")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom (F) and the disconnected template atom (C)");
    CHECK(hasBond(doc, a1, attachAtom), "target atom a1 is bonded directly to the attach atom");
    CHECK(!hasBond(doc, attachAtom, farAtom), "attach atom and the disconnected atom remain unbonded, matching the template's own topology");

    double ax = 0, ay = 0, fx = 0, fy = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    // 1-neighbor suggested angle: continuation (PI) + 60-degree kink (PI/3) = 4*PI/3, same
    // formula as test_graftAngleOrientation_oneNeighbor. templateAngle is nullopt here (no
    // internal bonds), so rotation is skipped -- but the attach atom still translates to land
    // exactly one bond length (1.5) from a1 along that angle, and the disconnected atom moves
    // by the same (dx, dy) translation, preserving its raw (2, 0) template-space offset from the
    // attach atom's own new position.
    double angle = 4.0 * M_PI / 3.0;
    double expectedAttachX = 1.5 * std::cos(angle);
    double expectedAttachY = 1.5 * std::sin(angle);
    CHECK(std::abs(ax - expectedAttachX) < 0.01 && std::abs(ay - expectedAttachY) < 0.01,
          "attach atom lands one bond length from a1 along the 1-neighbor suggested angle");
    CHECK(std::abs(fx - (expectedAttachX + 2.0)) < 0.01 && std::abs(fy - expectedAttachY) < 0.01,
          "disconnected atom keeps its raw (2,0) template-space offset from the attach atom's new position, unrotated");
}

static void test_graftAngleOrientation_threeNeighbors() {
    std::printf("--- Test: insertFunctionalGroup graft on a 3-neighbor target uses the largest-empty-angle fallback ---\n");
    // BondAngleSuggester::suggestAngle returns nullopt for 0-OR-3+-neighbor targets; the new
    // suggestFallbackAngle takes over for exactly this case.
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_3n"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.0, 0.0);   // angle 0
    AtomId a3 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 1.0);   // angle PI/2
    AtomId a4 = doc.molecule().addAtom(QStringLiteral("C"), -1.0, 0.0);  // angle PI
    doc.molecule().addBond(a1, a2, 1);
    doc.molecule().addBond(a1, a3, 1);
    doc.molecule().addBond(a1, a4, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 6, "graft added exactly two new atoms (attach atom survives)");
    AtomId attachAtom = -1, farAtom = -1;
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2 || id == a3 || id == a4) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("O")) attachAtom = id;
        else farAtom = id;
    }
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    CHECK(hasBond(doc, a1, attachAtom), "target atom is bonded directly to the attach atom");

    double ax = 0, ay = 0, fx = 0, fy = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    doc.molecule().atomPos(farAtom, fx, fy);
    // Existing neighbors at 0, PI/2, PI: gaps are PI/2 (0 to PI/2), PI/2 (PI/2 to PI), and PI
    // (PI wrapping around back to 0) -- the largest gap is PI, bisected at PI + PI/2 = 3*PI/2.
    double expected = 3.0 * M_PI / 2.0;
    double attachAngle = std::atan2(ay, ax);
    double farAngle = std::atan2(fy, fx);
    while (attachAngle < 0) attachAngle += 2 * M_PI;
    while (farAngle < 0) farAngle += 2 * M_PI;
    while (expected >= 2 * M_PI) expected -= 2 * M_PI;
    CHECK(std::abs(attachAngle - expected) < 0.01, "attach atom bisects the largest empty angle among the three existing neighbors");
    CHECK(std::abs(farAngle - expected) < 0.01, "far atom lands along the same rotated direction as the attach atom");
    CHECK(std::abs(std::hypot(ax, ay) - 1.5) < 0.01, "attach atom sits one bond length from the target atom");
    CHECK(std::abs(std::hypot(fx, fy) - 2.5) < 0.01, "far atom sits one more template-internal bond length further out");
}

static void test_graftAngleOrientation_multiBondAttachAtom() {
    std::printf("--- Test: insertFunctionalGroup graft uses mean neighbor direction for a multi-bond attach atom ---\n");
    // Regression test for the "first bonded neighbor" heuristic: the attach atom here has TWO
    // template-internal bonds (to atoms at 0 degrees and 120 degrees from it), not one. Picking
    // just the first of the two arbitrarily would rotate the template by the wrong amount. This
    // is checked directly: the template's own two internal branches, as seen FROM the surviving
    // attach atom, must point along rotateBy + {0, 120} degrees, where rotateBy uses the MEAN of
    // the attach atom's two internal-bond directions (60 degrees), not just the first (0
    // degrees). This also exercises the core sprout-not-merge behavior for a multi-bond attach
    // atom specifically: under the old fuse model, this attach atom's 2 internal bonds would
    // have been rewired onto a1 (fusing a1 into the template's own branch structure); under
    // sprout semantics a1 gains exactly one new bond (to the surviving attach atom), and the
    // attach atom's own 2 original bonds stay entirely its own.
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateTwoBondsAttach), QStringLiteral("test_graft_multibond"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.0, 0.0);     // angle 0 from a1
    AtomId a3 = doc.molecule().addAtom(QStringLiteral("C"), -0.5, 0.866); // angle 120 deg from a1
    doc.molecule().addBond(a1, a2, 1);
    doc.molecule().addBond(a1, a3, 1);

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 6, "graft added exactly three new atoms (attach atom N survives, plus its two original branch atoms)");

    AtomId attachAtom = -1;
    QList<AtomId> branchAtoms;
    for (AtomId id : allAtoms) {
        if (id == a1 || id == a2 || id == a3) continue;
        if (doc.molecule().atomSymbol(id) == QStringLiteral("N")) attachAtom = id;
        else branchAtoms.append(id);
    }
    CHECK(attachAtom != -1, "found the surviving attach atom (N)");
    CHECK(branchAtoms.size() == 2, "found both of the attach atom's original branch atoms");
    CHECK(hasBond(doc, a1, attachAtom), "target atom a1 is bonded directly to the attach atom");
    CHECK(hasBond(doc, attachAtom, branchAtoms[0]) && hasBond(doc, attachAtom, branchAtoms[1]),
          "attach atom keeps both of its own original template-internal bonds");
    CHECK(!hasBond(doc, a1, branchAtoms[0]) && !hasBond(doc, a1, branchAtoms[1]),
          "target atom a1 is NOT directly bonded to either branch atom (no rewiring, just one new bond)");

    double ax = 0, ay = 0;
    doc.molecule().atomPos(attachAtom, ax, ay);
    // suggestedAngle for a1's own 2 neighbors (at 0 and 120 degrees): bisects the larger gap,
    // landing at 4*PI/3 (240 degrees) -- same formula as test_graftAngleOrientation_twoNeighbors.
    // templateAngle (mean of the attach atom's own two 0-and-120-degree internal bonds) is
    // PI/3 (60 degrees). rotateBy = suggestedAngle - templateAngle = 4*PI/3 - PI/3 = PI.
    double rotateBy = M_PI;
    QList<double> expectedBranchAngles = { rotateBy, rotateBy + 2.0 * M_PI / 3.0 }; // rotated 0 deg and 120 deg
    for (double& e : expectedBranchAngles) { while (e < 0) e += 2 * M_PI; while (e >= 2 * M_PI) e -= 2 * M_PI; }

    for (AtomId id : branchAtoms) {
        double bx = 0, by = 0;
        doc.molecule().atomPos(id, bx, by);
        double angleFromAttach = std::atan2(by - ay, bx - ax);
        while (angleFromAttach < 0) angleFromAttach += 2 * M_PI;
        bool matchesEither = std::abs(angleFromAttach - expectedBranchAngles[0]) < 0.01
                           || std::abs(angleFromAttach - expectedBranchAngles[1]) < 0.01;
        CHECK(matchesEither, "branch atom's direction from the attach atom matches a rotation using the MEAN of both internal bond directions (60 degrees), not just the first (0 degrees)");
    }
}

static void test_graftAngleOrientation_undoRedoRoundTrip() {
    std::printf("--- Test: insertFunctionalGroup graft rotation survives undo/redo ---\n");
    QString path = writeTempSdf(QString::fromUtf8(kTestTemplateOneBond), QStringLiteral("test_graft_undoredo"));
    TemplateLibrary lib(path, path, path);

    DocumentState doc;
    AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
    doc.molecule().addBond(a1, a2, 1);
    int atomCountBeforeGraft = doc.molecule().atomIds().size();

    doc.insertFunctionalGroup(lib, QStringLiteral("TestFG"), 0.0, 0.0, a1, true);

    // AtomIds may be reassigned differently across an undo/redo cycle (freed index slots can
    // refill in a different order), so re-identify both new atoms by symbol each time rather
    // than assuming stable ids.
    auto findBySymbol = [&doc, a1, a2](const QString& symbol) -> AtomId {
        for (AtomId id : doc.molecule().atomIds()) {
            if (id == a1 || id == a2) continue;
            if (doc.molecule().atomSymbol(id) == symbol) return id;
        }
        return -1;
    };

    QList<AtomId> allAtoms = doc.molecule().atomIds();
    CHECK(allAtoms.size() == 4, "graft added exactly two new atoms (attach atom survives)");
    AtomId attachAtom = findBySymbol(QStringLiteral("O"));
    AtomId farAtom = findBySymbol(QStringLiteral("C"));
    CHECK(attachAtom != -1 && farAtom != -1, "found both the surviving attach atom and the far template atom");
    double origAX = 0, origAY = 0, origFX = 0, origFY = 0;
    doc.molecule().atomPos(attachAtom, origAX, origAY);
    doc.molecule().atomPos(farAtom, origFX, origFY);

    CHECK(doc.canUndo(), "graft pushed a history entry");
    doc.undo();
    CHECK(doc.molecule().atomIds().size() == atomCountBeforeGraft, "undo fully removes the grafted atoms");

    CHECK(doc.canRedo(), "undo enables redo");
    doc.redo();
    QList<AtomId> allAtomsAfterRedo = doc.molecule().atomIds();
    CHECK(allAtomsAfterRedo.size() == 4, "redo restores both grafted atoms");
    AtomId redoAttachAtom = findBySymbol(QStringLiteral("O"));
    AtomId redoFarAtom = findBySymbol(QStringLiteral("C"));
    CHECK(redoAttachAtom != -1 && redoFarAtom != -1, "found both re-grafted atoms after redo");
    double redoAX = 0, redoAY = 0, redoFX = 0, redoFY = 0;
    doc.molecule().atomPos(redoAttachAtom, redoAX, redoAY);
    doc.molecule().atomPos(redoFarAtom, redoFX, redoFY);
    CHECK(std::abs(redoAX - origAX) < 0.01 && std::abs(redoAY - origAY) < 0.01,
          "redo reproduces the exact same attach-atom position as the original graft (command replay is correctly snapshotted)");
    CHECK(std::abs(redoFX - origFX) < 0.01 && std::abs(redoFY - origFY) < 0.01,
          "redo reproduces the exact same far-atom position as the original graft");
}

static void test_insertFunctionalGroupLabel() {
    std::printf("--- Test: insertFunctionalGroup sets sgroup label ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 10.0, 10.0);
    QList<SGroupId> sgIds = doc.molecule().sgroupIds();
    CHECK(sgIds.size() == 1, "insertFunctionalGroup created exactly one sgroup");
    CHECK(doc.molecule().sgroupLabel(sgIds[0]) == QStringLiteral("Ac"),
          "the sgroup's label matches the inserted group's name");
}

static void test_renameSGroup() {
    std::printf("--- Test: renameSGroup ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    DocumentState doc;
    doc.insertFunctionalGroup(lib, QStringLiteral("Ac"), 10.0, 10.0);
    QList<SGroupId> sgIds = doc.molecule().sgroupIds();
    CHECK(sgIds.size() == 1, "setup: one sgroup exists");
    SGroupId sg = sgIds[0];
    CHECK(doc.molecule().sgroupLabel(sg) == QStringLiteral("Ac"), "setup: label starts as \"Ac\"");

    doc.renameSGroup(sg, QStringLiteral("Acetyl"));
    CHECK(doc.molecule().sgroupLabel(sg) == QStringLiteral("Acetyl"), "renameSGroup changed the label");
    CHECK(doc.canUndo(), "rename pushed a history entry");

    doc.undo();
    CHECK(doc.molecule().sgroupLabel(sg) == QStringLiteral("Ac"), "undo restores the old label");

    int historySizeBefore = doc.canUndo() ? 1 : 0;
    doc.renameSGroup(sg, QStringLiteral("Ac")); // same as current: no-op
    CHECK((doc.canUndo() ? 1 : 0) == historySizeBefore, "renaming to the same value pushes no history entry");

    doc.renameSGroup(9999, QStringLiteral("bogus")); // unknown id: no-op, no crash
    CHECK(doc.molecule().sgroupLabel(sg) == QStringLiteral("Ac"), "unknown id is a no-op");
}

static void test_clearCanvas() {
    std::printf("--- Test: clearCanvas ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("N"), 1, 0);
    doc.addBond(a1, a2, 1);
    doc.selectAtom(a1);
    EditableMolecule& mol = doc.molecule();
    CHECK(mol.atomCount() == 2, "setup: 2 atoms before clearCanvas");
    CHECK(!doc.selection().atoms.isEmpty(), "setup: something selected before clearCanvas");

    doc.clearCanvas();

    CHECK(mol.atomCount() == 0, "document is empty after clearCanvas");
    CHECK(doc.selection().isEmpty(), "selection cleared after clearCanvas");
    CHECK(doc.isDirty(), "document is dirty after clearCanvas (an executed command, not markClean)");
    CHECK(doc.canUndo(), "undo available");

    doc.undo();
    CHECK(mol.atomCount() == 2, "undo restores the original 2 atoms");
    QString sym1, sym2;
    for (AtomId id : mol.atomIds()) {
        if (sym1.isEmpty()) sym1 = mol.atomSymbol(id); else sym2 = mol.atomSymbol(id);
    }
    CHECK((sym1 == QStringLiteral("C") && sym2 == QStringLiteral("N")) ||
          (sym1 == QStringLiteral("N") && sym2 == QStringLiteral("C")),
          "restored atoms have their original symbols");
    CHECK(mol.bondCount() == 1, "restored bond too");
    CHECK(doc.selection().isEmpty(),
          "selection is EMPTY after undo, NOT restored -- undo()'s own generic post-invert "
          "clear always wins, matching the real JS's clearCanvas (its own selection-restore "
          "is dead code there for the identical reason)");
}

static void test_loadBenzene() {
    std::printf("--- Test: loadBenzene ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    (void)a1;
    EditableMolecule& mol = doc.molecule();
    CHECK(doc.canUndo(), "setup: history has an entry before loadBenzene");
    CHECK(doc.isDirty(), "setup: document is dirty before loadBenzene");

    doc.loadBenzene();

    CHECK(mol.atomCount() == 6, "6 atoms after loadBenzene");
    int carbonCount = 0;
    for (AtomId id : mol.atomIds()) if (mol.atomSymbol(id) == QStringLiteral("C")) ++carbonCount;
    CHECK(carbonCount == 6, "all 6 atoms are carbon");
    CHECK(mol.bondCount() == 6, "6 bonds (a ring)");
    CHECK(doc.selection().isEmpty(), "selection is empty after loadBenzene");

    CHECK(!doc.canUndo(), "history wiped -- cannot undo past loadBenzene");
    CHECK(!doc.canRedo(), "history wiped -- nothing to redo either");
    CHECK(!doc.isDirty(), "a freshly-loaded benzene document is NOT dirty");
}

static void test_documentStateBioSequenceViewForwarding() {
    std::printf("--- Test: DocumentState biopolymer sequence view forwarding ---\n");
    DocumentState doc;
    bool undoBefore = doc.canUndo();
    bool dirtyBefore = doc.isDirty();

    doc.buildBioSequenceView(QStringLiteral("ACDE"), QStringLiteral("PEPTIDE"));
    CHECK(doc.bioMonomers().size() == 4, "buildBioSequenceView forwards correctly");
    CHECK(doc.bioSeqType() == QStringLiteral("PEPTIDE"), "bioSeqType forwards correctly");

    doc.addBioMonomer(QStringLiteral("F"));
    CHECK(doc.bioMonomers().size() == 5, "addBioMonomer forwards correctly");
    CHECK(doc.bioBonds().size() == 4, "bioBonds forwards correctly");

    int lastId = doc.bioMonomers().last().id;
    doc.deleteBioMonomer(lastId);
    CHECK(doc.bioMonomers().size() == 4, "deleteBioMonomer forwards correctly");

    CHECK(doc.canUndo() == undoBefore, "none of these calls pushed a history entry");
    CHECK(doc.isDirty() == dirtyBefore, "none of these calls changed dirty state");
}

static void test_setMoleculeName() {
    std::printf("--- Test 26: setMoleculeName ---\n");
    DocumentState doc;
    CHECK(doc.molecule().name().isEmpty(), "default name is empty");

    doc.setMoleculeName(QStringLiteral("Aspirin"));
    CHECK(doc.molecule().name() == QStringLiteral("Aspirin"), "setMoleculeName sets the name");
    bool canUndoBefore = doc.canUndo();
    doc.setMoleculeName(QStringLiteral("Aspirin")); // unchanged -- must no-op
    CHECK(doc.canUndo() == canUndoBefore, "setMoleculeName to the SAME name pushes no history entry");
    doc.undo();
    CHECK(doc.molecule().name().isEmpty(), "undo restores the original (empty) name");

    doc.redo();
    CHECK(doc.molecule().name() == QStringLiteral("Aspirin"), "redo re-applies the name");
}

static void test_copySelection() {
    std::printf("--- Test 27: copySelection ---\n");
    DocumentState doc;
    CHECK(doc.copySelection().isEmpty(), "empty selection returns an empty string");

    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("O"), 1, 0);
    BondId b1 = doc.addBond(a1, a2, 1);
    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);
    doc.addBondToSelection(b1);

    bool canUndoBefore = doc.canUndo();
    QString copied = doc.copySelection();
    CHECK(!copied.isEmpty(), "copySelection returns real MOL text for a non-empty selection");
    int reparsed = indigoLoadMoleculeFromString(copied.toUtf8().constData());
    CHECK(reparsed >= 0, "copied text reparses cleanly");
    if (reparsed >= 0) {
        CHECK(indigoCountAtoms(reparsed) == 2 && indigoCountBonds(reparsed) == 1, "copied structure has exactly the selected atoms/bond");
        indigoFree(reparsed);
    }
    CHECK(doc.canUndo() == canUndoBefore, "copySelection is a pure read -- it pushes no history entry");
}

static void test_insertStructureAt() {
    std::printf("--- Test 28: insertStructureAt ---\n");
    DocumentState doc;

    doc.insertStructureAt(QStringLiteral("CCO"), 10, 5);
    CHECK(doc.molecule().atomCount() == 3, "insertStructureAt merges the source's atoms into the document");

    // Confirm the pasted structure's bbox center landed at (10, 5).
    double minX = 1e9, minY = 1e9, maxX = -1e9, maxY = -1e9;
    for (AtomId id : doc.molecule().atomIds()) {
        double x = 0, y = 0;
        doc.molecule().atomPos(id, x, y);
        minX = std::min(minX, x); maxX = std::max(maxX, x);
        minY = std::min(minY, y); maxY = std::max(maxY, y);
    }
    CHECK(std::abs((minX + maxX) / 2.0 - 10.0) < 1e-6, "pasted structure's bbox center X is at cx");
    CHECK(std::abs((minY + maxY) / 2.0 - 5.0) < 1e-6, "pasted structure's bbox center Y is at cy");

    doc.undo();
    CHECK(doc.molecule().atomCount() == 0, "undo removes every atom/bond the paste created");
    doc.redo();
    CHECK(doc.molecule().atomCount() == 3, "redo re-applies the paste");

    bool canUndoBefore = doc.canUndo();
    doc.insertStructureAt(QString(), 0, 0);
    CHECK(doc.canUndo() == canUndoBefore, "empty source is a no-op -- no history entry");
    doc.insertStructureAt(QStringLiteral("not a valid molecule $$$"), 0, 0);
    CHECK(doc.canUndo() == canUndoBefore, "unparseable source is a no-op -- no history entry");
    CHECK(doc.molecule().atomCount() == 3, "document is unchanged by either no-op call");
}

static void test_insertStructureAtClampsToPageBounds() {
    std::printf("--- Test 29: insertStructureAt clamps to page bounds ---\n");
    DocumentState doc;
    doc.insertStructureAt(QStringLiteral("C"), 9999.0, -9999.0);
    CHECK(doc.molecule().atomCount() == 1, "structure is still inserted even with wildly out-of-range coordinates");
    double x = 0, y = 0;
    doc.molecule().atomPos(doc.molecule().atomIds().first(), x, y);
    CHECK(x <= 30.0 && y >= -21.0, "resulting position is clamped within the page bounds, not left at the raw input");
}

static void test_copyPasteRoundTrip() {
    std::printf("--- Test 30: copySelection -> insertStructureAt round-trip ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("O"), 1, 0);
    doc.addBond(a1, a2, 1);
    doc.selectAtom(a1);
    doc.addAtomToSelection(a2);
    doc.addBondToSelection(doc.molecule().bondIds().first());

    QString copied = doc.copySelection();
    CHECK(!copied.isEmpty(), "setup: copySelection produced real text");

    doc.insertStructureAt(copied, 20, 20);
    CHECK(doc.molecule().atomCount() == 4, "the document now has both the original and the pasted atoms");
}

static void test_deleteSelectionEntitiesBondedPair() {
    std::printf("--- Test: deleteSelectionEntities on two bonded, individually-selected atoms ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = doc.addAtom(QStringLiteral("N"), 1.0, 0.0);
    BondId bond = doc.addBond(a, b, 1);
    CHECK(mol.atomCount() == 2 && mol.bondCount() == 1, "setup: 2 atoms, 1 bond");

    doc.selectAtom(a);
    doc.addAtomToSelection(b);
    doc.deleteSelectionEntities();

    CHECK(mol.atomCount() == 0, "both atoms gone");
    CHECK(mol.bondCount() == 0, "bond gone (cascade)");
    CHECK(doc.canUndo(), "undo available");

    doc.undo();
    CHECK(mol.atomCount() == 2, "both atoms restored (fresh ids)");
    CHECK(mol.bondCount() == 1, "bond restored between the new ids");
    QList<AtomId> restored = mol.atomIds();
    CHECK(restored.size() == 2, "exactly 2 atoms after undo");
    BondId restoredBond = mol.findBond(restored[0], restored[1]);
    CHECK(restoredBond != -1, "restored bond connects the two restored atoms");
    (void)bond;
}

static void test_deleteSelectionEntitiesWholePill() {
    std::printf("--- Test: deleteSelectionEntities on a whole sgroup pill selected directly ---\n");
    // Built programmatically, NOT via a molfile constructor: AtomId/SGroupId are independent
    // counters both starting at 1, so a molfile-loaded sgroup (the FIRST sgroup ever created
    // for a fresh molecule) always gets SGroupId=1 -- which collides with AtomId=1 (the first
    // atom) whenever one exists, exactly the numeric-collision hazard this spec's own
    // deleteSelectionEntities design explicitly guards against. A decoy atom is added and
    // removed first so SGroupId=1 does NOT alias any real, currently-existing AtomId, letting
    // this test actually exercise the whole-pill-selection path instead of silently
    // misclassifying it as "atom 1 selected" (which corrupts sgroup state and crashes on undo).
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId decoy = mol.addAtom(QStringLiteral("Xe"), 100.0, 100.0);
    mol.removeAtom(decoy);   // frees AtomId 1 so it no longer collides with the sgroup id below

    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = mol.addAtom(QStringLiteral("C"), 1.0, 0.0);
    AtomId a3 = mol.addAtom(QStringLiteral("O"), 2.0, 0.0);
    mol.addBond(a1, a2, 1);
    mol.addBond(a2, a3, 1);
    SGroupId sid = mol.createSuperatomFromAtoms({a1, a2, a3});
    CHECK(sid != -1, "setup: sgroup created");
    CHECK(mol.sgroupIds().size() == 1, "1 sgroup on setup");
    CHECK(!mol.atomIds().contains(sid), "setup: sgroup id does not collide with any real atom id");

    doc.selection().atoms.insert(sid);   // whole-pill selection: the SGroupId stands in
                                          // for its atoms in m_selection.atoms, per spec.
    doc.deleteSelectionEntities();

    CHECK(mol.sgroupIds().isEmpty(), "sgroup grouping gone");
    CHECK(mol.atomCount() == 3, "member atoms untouched");

    doc.undo();
    CHECK(mol.sgroupIds().size() == 1, "sgroup restored on undo");
    CHECK(mol.atomCount() == 3, "still exactly the original 3 atoms (never removed)");
}

static void test_deleteSelectionEntitiesPartialMember() {
    std::printf("--- Test: deleteSelectionEntities on ONE of a sgroup's several members ---\n");
    const char* molfile =
        "acFixture\n"
        "  Ketcher\n\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "M  STY  1   1 SUP\n"
        "M  SAL   1  3   1   2   3\n"
        "M  SMT   1 Ac\n"
        "M  END\n";
    DocumentState doc(molfile);
    EditableMolecule& mol = doc.molecule();
    QList<AtomId> allAtoms = mol.atomIds();
    CHECK(allAtoms.size() == 3, "3 atoms on load");
    AtomId victim = allAtoms.last();   // the terminal O, member of the sgroup

    doc.selectAtom(victim);
    doc.deleteSelectionEntities();

    CHECK(mol.atomCount() == 2, "1 atom removed");
    CHECK(mol.sgroupIds().size() == 1, "sgroup survives (trimmed, not fully removed)");
    SGroupId sid = mol.sgroupIds().first();
    CHECK(mol.sgroupMemberAtomIds(sid).size() == 2, "sgroup member list trimmed to 2");

    doc.undo();
    CHECK(mol.atomCount() == 3, "atom restored");
    CHECK(mol.sgroupIds().size() == 1, "still 1 sgroup after undo");
    CHECK(mol.sgroupMemberAtomIds(mol.sgroupIds().first()).size() == 3,
          "sgroup's full original 3-atom member set restored");
}

static void test_deleteSelectionEntitiesMixedTypes() {
    std::printf("--- Test: deleteSelectionEntities on a mix of atoms/bonds/rxnArrow/rxnPlus ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    RxnArrowId arrow = doc.addRxnArrow(5.0, 0.0);
    RxnPlusId plus = doc.addRxnPlus(-5.0, 0.0);
    CHECK(mol.atomCount() == 1 && mol.rxnArrowCount() == 1 && mol.rxnPlusCount() == 1,
          "setup: 1 atom, 1 arrow, 1 plus");

    doc.selectAtom(a);
    doc.addRxnArrowToSelection(arrow);
    doc.addRxnPlusToSelection(plus);
    doc.deleteSelectionEntities();

    CHECK(mol.atomCount() == 0, "atom gone");
    CHECK(mol.rxnArrowCount() == 0, "arrow gone");
    CHECK(mol.rxnPlusCount() == 0, "plus gone");

    doc.undo();
    CHECK(mol.atomCount() == 1, "atom restored");
    CHECK(mol.rxnArrowCount() == 1, "arrow restored");
    CHECK(mol.rxnPlusCount() == 1, "plus restored");
}

static void test_deleteSelectionEntitiesRedoAfterUndo() {
    std::printf("--- Test: deleteSelectionEntities redo-after-undo (regression) ---\n");
    // Real bug found via live UI testing during sub-project 7b: execute() captured the
    // ORIGINAL atom/bond/etc ids by value, but invert() recreates everything with BRAND NEW
    // ids (this port's own established fresh-id-on-redo convention) -- so a second execute()
    // (redo after undo) tried to remove ids that no longer existed, and EditableMolecule's
    // remove* methods silently return false on an unknown id, making the bug invisible except
    // as "redo after undo does nothing." No prior test exercised delete -> undo -> redo.
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = doc.addAtom(QStringLiteral("O"), 1.0, 0.0);
    doc.addBond(a, b, 1);
    RxnArrowId arrow = doc.addRxnArrow(5.0, 0.0);
    RxnPlusId plus = doc.addRxnPlus(-5.0, 0.0);
    MultitailArrowId mta = doc.addMultitailArrow(3.0, 3.0);
    CHECK(mol.atomCount() == 2 && mol.bondIds().size() == 1 && mol.rxnArrowCount() == 1 &&
          mol.rxnPlusCount() == 1 && mol.multitailArrowCount() == 1,
          "setup: 2 atoms, 1 bond, 1 arrow, 1 plus, 1 multitail arrow");

    doc.selectAtom(a);
    doc.addAtomToSelection(b);
    doc.addRxnArrowToSelection(arrow);
    doc.addRxnPlusToSelection(plus);
    doc.addMultitailArrowToSelection(mta);
    doc.deleteSelectionEntities();
    CHECK(mol.atomCount() == 0 && mol.rxnArrowCount() == 0 && mol.rxnPlusCount() == 0 &&
          mol.multitailArrowCount() == 0, "first delete: everything gone");

    doc.undo();
    CHECK(mol.atomCount() == 2 && mol.bondIds().size() == 1 && mol.rxnArrowCount() == 1 &&
          mol.rxnPlusCount() == 1 && mol.multitailArrowCount() == 1,
          "undo: everything restored (with fresh ids, per this port's own convention)");

    doc.redo();
    CHECK(mol.atomCount() == 0 && mol.rxnArrowCount() == 0 && mol.rxnPlusCount() == 0 &&
          mol.multitailArrowCount() == 0,
          "REQUIRED regression check: redo after undo actually re-deletes everything -- "
          "previously silently did nothing, since it tried to remove the stale original ids");

    // A second undo/redo cycle must also work, proving the refreshed ids produced by the
    // FIRST redo's own invert() call are themselves correctly tracked, not just the first.
    doc.undo();
    CHECK(mol.atomCount() == 2 && mol.bondIds().size() == 1 && mol.rxnArrowCount() == 1 &&
          mol.rxnPlusCount() == 1 && mol.multitailArrowCount() == 1,
          "second undo: everything restored again");
    doc.redo();
    CHECK(mol.atomCount() == 0 && mol.rxnArrowCount() == 0 && mol.rxnPlusCount() == 0 &&
          mol.multitailArrowCount() == 0, "second redo: everything gone again");
}

static void test_deleteSelectionEntitiesWholePillRedoAfterUndo() {
    std::printf("--- Test: deleteSelectionEntities whole-pill-sgroup redo-after-undo (regression) ---\n");
    // Built programmatically, not via a molfile constructor -- same decoy-atom pattern as
    // test_deleteSelectionEntitiesWholePill above, needed so SGroupId=1 does not collide with
    // AtomId=1 (independent counters, both start at 1) and get misclassified as "atom 1
    // selected" instead of a whole-pill selection.
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId decoy = mol.addAtom(QStringLiteral("Xe"), 100.0, 100.0);
    mol.removeAtom(decoy);

    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = mol.addAtom(QStringLiteral("C"), 1.0, 0.0);
    AtomId a3 = mol.addAtom(QStringLiteral("O"), 2.0, 0.0);
    mol.addBond(a1, a2, 1);
    mol.addBond(a2, a3, 1);
    SGroupId sid = mol.createSuperatomFromAtoms({a1, a2, a3});
    CHECK(sid != -1 && !mol.atomIds().contains(sid), "setup: sgroup created, id does not collide");
    CHECK(mol.sgroupIds().size() == 1, "1 sgroup on setup");

    doc.selection().atoms.insert(sid);   // whole-pill selection: the SGroupId stands in
                                          // for its atoms in m_selection.atoms, per spec.
    doc.deleteSelectionEntities();
    CHECK(mol.sgroupIds().isEmpty(), "first delete: sgroup gone");

    doc.undo();
    CHECK(mol.sgroupIds().size() == 1, "undo: sgroup restored (fresh SGroupId)");

    doc.redo();
    CHECK(mol.sgroupIds().isEmpty(),
          "REQUIRED regression check: redo after undo removes the RESTORED sgroup's fresh id, "
          "not the stale original one");
}

static void test_cutSelection() {
    std::printf("--- Test: cutSelection (copy-then-delete round trip) ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = doc.addAtom(QStringLiteral("O"), 1.0, 0.0);
    doc.addBond(a, b, 1);

    doc.selectAtom(a);
    doc.addAtomToSelection(b);
    QString expectedCopy = doc.copySelection();
    CHECK(!expectedCopy.isEmpty(), "setup: copySelection produces content before cut");

    QString cutResult = doc.cutSelection();
    CHECK(cutResult == expectedCopy, "cutSelection returns the same text copySelection would have");
    CHECK(mol.atomCount() == 0, "selected content is gone from the live document");

    doc.undo();
    CHECK(mol.atomCount() == 2, "single undo restores everything");
}

static void test_setStereoDescriptorsValid() {
    std::printf("--- Test: setStereoDescriptors applies valid JSON ---\n");
    DocumentState doc;
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 1, 0);
    doc.addBond(a1, a2, 1);
    EditableMolecule& mol = doc.molecule();

    QList<AtomId> order = mol.atomIdsInIndigoOrder();
    CHECK(order.size() == 2, "setup: 2 atoms");
    int idx0 = order.indexOf(a1), idx1 = order.indexOf(a2);

    QString json = QString(
        "{\"atoms\":{\"%1\":{\"cipLabel\":\"R\",\"type\":1,\"group\":2}},"
        "\"bonds\":{\"%2-%3\":{\"cipLabel\":\"E\"}}}"
    ).arg(idx0).arg(qMin(idx0, idx1)).arg(qMax(idx0, idx1));

    bool historyBefore = doc.canUndo();
    doc.setStereoDescriptors(json);

    CHECK(mol.atomStoredCipLabel(a1) == QStringLiteral("R"), "atom cipLabel applied");
    CHECK(mol.atomStoredStereoType(a1) == 1, "atom stereoType applied");
    CHECK(mol.atomStoredStereoGroup(a1) == 2, "atom stereoGroup applied");
    CHECK(mol.atomStoredCipLabel(a2).isEmpty(), "the OTHER atom is untouched");
    BondId bond = mol.findBond(a1, a2);
    CHECK(mol.bondStoredCipLabel(bond) == QStringLiteral("E"), "bond cipLabel applied");
    CHECK(doc.canUndo() == historyBefore, "no history entry pushed");
}

static void test_setStereoDescriptorsReuseOrderRegression() {
    std::printf("--- Test: setStereoDescriptors resolves correctly after an index-reuse edit ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("N"), 1, 0);
    AtomId a3 = doc.addAtom(QStringLiteral("O"), 2, 0);
    mol.removeAtom(a2);
    AtomId a4 = doc.addAtom(QStringLiteral("Cl"), 3, 0);

    QList<AtomId> order = mol.atomIdsInIndigoOrder();
    int idxOfA4 = order.indexOf(a4);
    QString json = QString("{\"atoms\":{\"%1\":{\"cipLabel\":\"S\",\"type\":1,\"group\":0}},\"bonds\":{}}")
                       .arg(idxOfA4);

    doc.setStereoDescriptors(json);

    CHECK(mol.atomStoredCipLabel(a4) == QStringLiteral("S"),
          "the atom that reused the freed index gets the label -- NOT a1/a3");
    CHECK(mol.atomStoredCipLabel(a1).isEmpty(), "a1 untouched");
    CHECK(mol.atomStoredCipLabel(a3).isEmpty(), "a3 untouched");
}

static void test_setStereoDescriptorsResetsOnEachCall() {
    std::printf("--- Test: setStereoDescriptors replaces, does not accumulate ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx0 = mol.atomIdsInIndigoOrder().indexOf(a1);

    doc.setStereoDescriptors(QString("{\"atoms\":{\"%1\":{\"cipLabel\":\"R\",\"type\":1,\"group\":0}},\"bonds\":{}}").arg(idx0));
    CHECK(mol.atomStoredCipLabel(a1) == QStringLiteral("R"), "first call applies R");

    doc.setStereoDescriptors(QString("{\"atoms\":{},\"bonds\":{}}"));
    CHECK(mol.atomStoredCipLabel(a1).isEmpty(), "second call with empty payload clears it (reset ran)");
}

static void test_setStereoDescriptorsMalformedJsonIsNoOp() {
    std::printf("--- Test: setStereoDescriptors ignores malformed JSON ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx0 = mol.atomIdsInIndigoOrder().indexOf(a1);
    doc.setStereoDescriptors(QString("{\"atoms\":{\"%1\":{\"cipLabel\":\"R\",\"type\":1,\"group\":0}},\"bonds\":{}}").arg(idx0));
    CHECK(mol.atomStoredCipLabel(a1) == QStringLiteral("R"), "setup: R applied");

    doc.setStereoDescriptors(QStringLiteral("not valid json{{{"));
    CHECK(mol.atomStoredCipLabel(a1) == QStringLiteral("R"), "malformed JSON leaves prior data unchanged");
}

static void test_setStereoDescriptorsOutOfRangeIndexSkipped() {
    std::printf("--- Test: setStereoDescriptors skips an out-of-range index safely ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx0 = mol.atomIdsInIndigoOrder().indexOf(a1);

    QString json = QString(
        "{\"atoms\":{\"%1\":{\"cipLabel\":\"R\",\"type\":1,\"group\":0},"
        "\"99\":{\"cipLabel\":\"S\",\"type\":1,\"group\":0}},\"bonds\":{}}"
    ).arg(idx0);
    doc.setStereoDescriptors(json);

    CHECK(mol.atomStoredCipLabel(a1) == QStringLiteral("R"),
          "the valid entry still applies even though the payload also has an out-of-range one");
}

static void test_setCheckIssuesValid() {
    std::printf("--- Test: setCheckIssues applies atom- and bond-target issues ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 1, 0);
    BondId bond = doc.addBond(a1, a2, 1);

    QList<AtomId> atomOrder = mol.atomIdsInIndigoOrder();
    QList<BondId> bondOrder = mol.bondIdsInIndigoOrder();
    int atomIdx = atomOrder.indexOf(a1);
    int bondIdx = bondOrder.indexOf(bond);

    bool historyBefore = doc.canUndo();
    QString json = QString(
        "{\"issues\":["
        "{\"target\":\"atom\",\"type\":\"valence\",\"ids\":[%1]},"
        "{\"target\":\"bond\",\"type\":\"overlap_atom\",\"ids\":[%2]}"
        "]}"
    ).arg(atomIdx).arg(bondIdx);
    doc.setCheckIssues(json);

    CHECK(mol.atomCheckWarningText(a1) == QStringLiteral("valence"), "atom checkWarning applied");
    CHECK(mol.atomCheckWarningText(a2).isEmpty(), "the OTHER atom is untouched");
    CHECK(mol.bondCheckWarningText(bond) == QStringLiteral("overlap_atom"), "bond checkWarning applied");
    CHECK(doc.canUndo() == historyBefore, "no history entry pushed");
}

static void test_setCheckIssuesDefaultTargetIsAtom() {
    std::printf("--- Test: setCheckIssues defaults missing/unrecognized target to atom ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx = mol.atomIdsInIndigoOrder().indexOf(a1);

    doc.setCheckIssues(QString("{\"issues\":[{\"type\":\"radical\",\"ids\":[%1]}]}").arg(idx));
    CHECK(mol.atomCheckWarningText(a1) == QStringLiteral("radical"),
          "an issue with no target field applies to the atom, matching the real JS's || \"atom\" default");
}

static void test_setCheckIssuesResetsOnEachCall() {
    std::printf("--- Test: setCheckIssues replaces, does not accumulate ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx = mol.atomIdsInIndigoOrder().indexOf(a1);

    doc.setCheckIssues(QString("{\"issues\":[{\"target\":\"atom\",\"type\":\"valence\",\"ids\":[%1]}]}").arg(idx));
    CHECK(mol.atomCheckWarningText(a1) == QStringLiteral("valence"), "first call applies");

    doc.setCheckIssues(QStringLiteral("{\"issues\":[]}"));
    CHECK(mol.atomCheckWarningText(a1).isEmpty(), "second call with no issues clears it (reset ran)");
}

static void test_setCheckIssuesMalformedJsonIsNoOp() {
    std::printf("--- Test: setCheckIssues ignores malformed JSON ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    int idx = mol.atomIdsInIndigoOrder().indexOf(a1);
    doc.setCheckIssues(QString("{\"issues\":[{\"target\":\"atom\",\"type\":\"valence\",\"ids\":[%1]}]}").arg(idx));
    CHECK(mol.atomCheckWarningText(a1) == QStringLiteral("valence"), "setup: applied");

    doc.setCheckIssues(QStringLiteral("{not json"));
    CHECK(mol.atomCheckWarningText(a1) == QStringLiteral("valence"), "malformed JSON leaves prior data unchanged");
}

static void test_setCheckIssuesReuseOrderRegression() {
    std::printf("--- Test: setCheckIssues resolves correctly after a bond index-reuse edit ---\n");
    DocumentState doc;
    EditableMolecule& mol = doc.molecule();
    AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = doc.addAtom(QStringLiteral("C"), 1, 0);
    AtomId a3 = doc.addAtom(QStringLiteral("C"), 2, 0);
    BondId b1 = doc.addBond(a1, a2, 1);
    doc.addBond(a2, a3, 1);
    mol.removeBond(b1);
    BondId b3 = doc.addBond(a1, a3, 2);

    QList<BondId> order = mol.bondIdsInIndigoOrder();
    int idx = order.indexOf(b3);
    doc.setCheckIssues(QString("{\"issues\":[{\"target\":\"bond\",\"type\":\"overlap_atom\",\"ids\":[%1]}]}").arg(idx));

    CHECK(mol.bondCheckWarningText(b3) == QStringLiteral("overlap_atom"),
          "the bond that reused the freed index gets the warning");
}

static void test_changeBondTypeAndStereo() {
    std::printf("--- Test: changeBondTypeAndStereo (atomic order+stereo undo) ---\n");
    DocumentState doc;

    // Same proven 4-substituent stereocenter shape as test_documentStateEditingOperations'
    // own setBondStereo test -- a 2-substituent origin is chemically invalid for a wedge and
    // Indigo silently rejects it (found during sub-project 8's own execution).
    AtomId center = doc.addAtom(QStringLiteral("C"), 10, 10);
    AtomId f = doc.addAtom(QStringLiteral("F"), 10, 11);
    AtomId cl = doc.addAtom(QStringLiteral("Cl"), 9, 9.5);
    AtomId br = doc.addAtom(QStringLiteral("Br"), 11, 9.5);
    doc.addBond(center, f, 1);
    doc.addBond(center, cl, 1);
    BondId bond = doc.addBond(center, br, 1);

    CHECK(doc.molecule().bondOrder(bond) == 1, "bond starts as a single bond");
    CHECK(doc.molecule().bondStereoDirectionEnum(bond) == EditableMolecule::Direction::None,
          "bond starts with no stereo");

    // 1) Setting stereo (order unchanged) sets Direction::Up from the raw V2000 code 1.
    doc.changeBondTypeAndStereo(bond, 1, 1);
    CHECK(doc.molecule().bondOrder(bond) == 1, "order unchanged (still 1)");
    CHECK(doc.molecule().bondStereoDirectionEnum(bond) == EditableMolecule::Direction::Up,
          "V2000 code 1 maps to Direction::Up");

    // 2) One undo restores BOTH fields in a single step (not two separate undos).
    doc.undo();
    CHECK(doc.molecule().bondOrder(bond) == 1, "undo restores order to 1");
    CHECK(doc.molecule().bondStereoDirectionEnum(bond) == EditableMolecule::Direction::None,
          "undo restores stereo to None in the SAME undo as the order restore");

    // 3) Redo re-applies both together.
    doc.redo();
    CHECK(doc.molecule().bondStereoDirectionEnum(bond) == EditableMolecule::Direction::Up,
          "redo re-applies Up");

    // 4) Changing ORDER away from 1 (stereo request 0/None) -- setBondStereo's own guard
    // (bondOrder(id) != 1) makes the stereo-clear call a safe no-op once order is 2, but the
    // wedge must still be fully restorable on undo (this is the whole reason execute/invert
    // both always set order FIRST, then attempt stereo -- see this method's own comment).
    doc.changeBondTypeAndStereo(bond, 2, 0);
    CHECK(doc.molecule().bondOrder(bond) == 2, "order changed to 2 (double bond)");

    doc.undo();
    CHECK(doc.molecule().bondOrder(bond) == 1, "undo restores order to 1");
    CHECK(doc.molecule().bondStereoDirectionEnum(bond) == EditableMolecule::Direction::Up,
          "undo ALSO restores the original Up wedge, even though order changed away from 1 and back");

    // 5) Guard: calling with identical order+stereo pushes no history entry.
    bool canUndoBefore = doc.canUndo();
    doc.changeBondTypeAndStereo(bond, 1, 1); // already order=1, stereo=Up (V2000 code 1) -- no-op
    CHECK(doc.canUndo() == canUndoBefore, "no-op call (unchanged order and stereo) pushes no history entry");
}

static void test_selectByRectAndLasso() {
    std::printf("--- Test: selectByRect / addSelectionByRect / selectByLasso ---\n");
    DocumentState doc;

    // Chain A-B-C, positioned so a rect covering only A and B still has its B-C bond
    // edge crossing the rect boundary (B is at x=10, C is at x=20; a rect ending at
    // x=15 clips the B-C bond segment without containing C).
    AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b = doc.addAtom(QStringLiteral("C"), 10, 0);
    AtomId c = doc.addAtom(QStringLiteral("C"), 20, 0);
    BondId abBond = doc.addBond(a, b, 1);
    BondId bcBond = doc.addBond(b, c, 1);

    // 1) selectByRect: rect [-5,-5]-[15,5] covers A and B fully; the B-C bond's segment
    // (10,0)-(20,0) crosses the rect's right edge at x=15 even though C (20,0) is outside.
    doc.selectByRect(-5, -5, 15, 5);
    CHECK(doc.selection().atoms.contains(a), "selectByRect: A is inside the rect");
    CHECK(doc.selection().atoms.contains(b), "selectByRect: B is inside the rect");
    CHECK(!doc.selection().atoms.contains(c), "selectByRect: C is outside the rect");
    CHECK(doc.selection().bonds.contains(abBond), "selectByRect: A-B bond selected (both endpoints inside)");
    CHECK(doc.selection().bonds.contains(bcBond),
          "selectByRect: B-C bond selected via edge-crossing even though C is outside -- the specific behavior that makes selectByRect different from a naive both-endpoints test");

    // 2) addSelectionByRect: same geometry, but starting from an EMPTY selection --
    // must NOT select the B-C bond (strict both-endpoints-inside, no edge-crossing test),
    // proving the two functions are genuinely asymmetric, not a copy-paste bug.
    DocumentState doc2;
    AtomId a2 = doc2.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b2 = doc2.addAtom(QStringLiteral("C"), 10, 0);
    AtomId c2 = doc2.addAtom(QStringLiteral("C"), 20, 0);
    BondId ab2 = doc2.addBond(a2, b2, 1);
    BondId bc2 = doc2.addBond(b2, c2, 1);
    doc2.addSelectionByRect(-5, -5, 15, 5);
    CHECK(doc2.selection().atoms.contains(a2), "addSelectionByRect: A is inside the rect");
    CHECK(doc2.selection().atoms.contains(b2), "addSelectionByRect: B is inside the rect");
    CHECK(!doc2.selection().atoms.contains(c2), "addSelectionByRect: C is outside the rect");
    CHECK(ab2 != bc2, "setup sanity: A-B and B-C are different bonds");
    CHECK(doc2.selection().bonds.contains(ab2), "addSelectionByRect: A-B bond selected (both endpoints inside)");
    CHECK(!doc2.selection().bonds.contains(bc2),
          "addSelectionByRect: B-C bond NOT selected -- strict both-endpoints-inside, no edge-crossing fallback (the asymmetry with selectByRect)");

    // 3) addSelectionByRect called twice with overlapping rects must not duplicate ids
    // (QSet already guarantees this, but confirm the ADDING behavior itself: a second,
    // non-overlapping rect call should ADD c2 without losing a2/b2).
    doc2.addSelectionByRect(15, -5, 25, 5);
    CHECK(doc2.selection().atoms.contains(a2) && doc2.selection().atoms.contains(b2) && doc2.selection().atoms.contains(c2),
          "addSelectionByRect accumulates across calls instead of replacing");

    // 4) selectByLasso: a triangle enclosing A and B but not C.
    QList<QPointF> triangle = { QPointF(-5, -5), QPointF(15, -5), QPointF(5, 10) };
    doc.selectByLasso(triangle);
    CHECK(doc.selection().atoms.contains(a), "selectByLasso: A is inside the triangle");
    CHECK(doc.selection().atoms.contains(b), "selectByLasso: B is inside the triangle");
    CHECK(!doc.selection().atoms.contains(c), "selectByLasso: C is outside the triangle");
    CHECK(doc.selection().bonds.contains(abBond), "selectByLasso: A-B bond selected (both endpoints inside)");
    CHECK(!doc.selection().bonds.contains(bcBond),
          "selectByLasso: B-C bond NOT selected -- lasso requires full enclosure, no edge-crossing fallback (matches ChemDraw semantics, not selectByRect's)");

    // 5) selectByLasso with fewer than 3 points clears the selection and does not crash.
    doc.selectByLasso({ QPointF(0, 0), QPointF(1, 1) });
    CHECK(doc.selection().isEmpty(), "selectByLasso with < 3 points clears the selection (matches the real < 6-flat-numbers guard)");

    // 6) rxnArrow/rxnPlus coverage: selectByRect and selectByLasso include them;
    // addSelectionByRect does not touch them at all.
    RxnArrowId arrow = doc.addRxnArrow(50, 50);
    RxnPlusId plus = doc.addRxnPlus(60, 60);
    doc.selectByRect(45, 45, 65, 65);
    CHECK(doc.selection().rxnArrows.contains(arrow), "selectByRect includes an rxnArrow fully inside the rect");
    CHECK(doc.selection().rxnPluses.contains(plus), "selectByRect includes an rxnPlus fully inside the rect");

    // Margin note (found during plan self-review): vertices chosen with comfortable
    // margin around both the rxnArrow segment ((50,50)-(~53.75,50)) and the rxnPlus point
    // (60,60) -- an earlier draft placed (60,60) exactly ON a triangle edge, which is
    // undefined/flaky for ray-casting point-in-polygon (boundary points can go either way
    // depending on floating-point rounding). At y=60 this triangle spans x=[47.5,72.5],
    // giving 60 a 12.5-unit margin on both sides.
    QList<QPointF> bigTriangle = { QPointF(35, 35), QPointF(85, 35), QPointF(60, 85) };
    doc.selectByLasso(bigTriangle);
    CHECK(doc.selection().rxnArrows.contains(arrow), "selectByLasso includes an rxnArrow fully inside the polygon");
    CHECK(doc.selection().rxnPluses.contains(plus), "selectByLasso includes an rxnPlus fully inside the polygon");

    DocumentState doc3;
    RxnArrowId arrow3 = doc3.addRxnArrow(50, 50);
    doc3.addSelectionByRect(45, 45, 65, 65);
    CHECK(!doc3.selection().rxnArrows.contains(arrow3),
          "addSelectionByRect never touches rxnArrows -- matches the real function exactly (it has no rxnArrow/rxnPlus code path at all)");
}

static void test_selectFragmentRingChain() {
    std::printf("--- Test: selectFragment / selectRing / selectChain ---\n");

    // ---- selectFragment: two disconnected fragments ----
    {
        DocumentState doc;
        AtomId a1 = doc.addAtom(QStringLiteral("C"), 0, 0);
        AtomId a2 = doc.addAtom(QStringLiteral("C"), 1, 0);
        BondId b1 = doc.addBond(a1, a2, 1);
        AtomId c1 = doc.addAtom(QStringLiteral("C"), 10, 10);
        AtomId c2 = doc.addAtom(QStringLiteral("C"), 11, 10);
        BondId b2 = doc.addBond(c1, c2, 1);

        doc.selectFragment(a1, -1);
        CHECK(doc.selection().atoms.contains(a1) && doc.selection().atoms.contains(a2),
              "selectFragment: both atoms of the clicked fragment are selected");
        CHECK(!doc.selection().atoms.contains(c1) && !doc.selection().atoms.contains(c2),
              "selectFragment: the OTHER fragment's atoms are not selected");
        CHECK(doc.selection().bonds.contains(b1) && !doc.selection().bonds.contains(b2),
              "selectFragment: only the clicked fragment's bond is selected");

        // Resolve from a bond instead of an atom (bond.begin is the resolved atomId).
        doc.selectFragment(-1, b2);
        CHECK(doc.selection().atoms.contains(c1) && doc.selection().atoms.contains(c2),
              "selectFragment(bond-resolved): the OTHER fragment is now selected");
        CHECK(!doc.selection().atoms.contains(a1),
              "selectFragment(bond-resolved): the first fragment is no longer selected (REPLACES)");
    }

    // ---- selectFragment: null/invalid-id edge cases (the exact asymmetry review found) ----
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
        doc.selectAtom(a);
        CHECK(!doc.selection().isEmpty(), "setup: something is selected");

        doc.selectFragment(-1, 9999); // invalid bondId, atomId null -> CLEARS (not a no-op)
        CHECK(doc.selection().isEmpty(),
              "selectFragment(-1, invalid bond id) clears the selection -- atomId stays null "
              "after the failed bond resolution, landing in the same clear-branch as both-null");

        doc.selectAtom(a);
        CHECK(!doc.selection().isEmpty(), "setup: something is selected again");

        doc.selectFragment(9999, -1); // atomId given directly but doesn't exist -> NO-OP
        CHECK(doc.selection().atoms.contains(a),
              "selectFragment(nonexistent atomId, -1) leaves a pre-existing selection "
              "untouched -- a DIFFERENT outcome than the invalid-bondId case above, despite "
              "both inputs being \"nothing findable\"");
    }

    // ---- selectRing: a 6-ring built by fusing a new hexagon onto a pre-existing bond
    // (same fixture shape as the existing test_addRing fusion-seam test). Topologically
    // this is a SINGLE 6-membered ring (6 atoms, 6 bonds -- Euler's formula gives exactly
    // 1 cycle), not two independent fused rings -- a0's two incident bonds (the original
    // seam, and its edge into the newly-fused hexagon) both lead to the SAME ring via
    // shortestRingThroughBond, so this exercises the atom->incident-bonds->shortest-ring
    // pipeline without needing a genuine multi-ring tie-break (that tie-breaking rule is
    // a documented "first found wins" simplification per the spec, not independently
    // asserted here). ----
    {
        DocumentState doc;
        AtomId a0 = doc.molecule().addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a1 = doc.molecule().addAtom(QStringLiteral("C"), 1.5, 0.0);
        doc.molecule().addBond(a0, a1, 1);
        QList<double> coords;
        coords << 0.0 << 0.0 << 1.5 << 0.0;
        coords << 2.25 << 1.3 << 1.5 << 2.6 << 0.0 << 2.6 << -0.75 << 1.3;
        doc.addRing(coords, true); // fuses onto a0/a1, completing the hexagon

        CHECK(doc.molecule().atomCount() == 6, "setup: the fused hexagon has 6 atoms total");

        doc.selectRing(a0, -1);
        CHECK(doc.selection().atoms.size() == 6,
              "selectRing from a0 finds the 6-membered ring it belongs to");
        CHECK(doc.selection().atoms.contains(a0),
              "selectRing: the clicked atom itself is part of the returned ring");
    }

    // ---- selectRing: isolated atom (no ring) leaves selection untouched ----
    {
        DocumentState doc;
        AtomId lone = doc.addAtom(QStringLiteral("C"), 0, 0);
        AtomId other = doc.addAtom(QStringLiteral("C"), 5, 5);
        doc.selectAtom(other);
        CHECK(doc.selection().atoms.contains(other), "setup: `other` is selected");

        doc.selectRing(lone, -1);
        CHECK(doc.selection().atoms.contains(other) && !doc.selection().atoms.contains(lone),
              "selectRing on an atom with no ring leaves the PRIOR selection untouched -- "
              "not a clear, matches the real code's early return exactly");
    }

    // ---- selectRing: both-null clears; invalid bondId also leaves selection untouched ----
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
        doc.selectAtom(a);
        doc.selectRing(-1, -1);
        CHECK(doc.selection().isEmpty(), "selectRing(-1, -1) clears the selection");

        doc.selectAtom(a);
        doc.selectRing(-1, 9999);
        CHECK(doc.selection().atoms.contains(a),
              "selectRing(-1, invalid bond id) leaves the selection untouched (NOT a clear -- "
              "asymmetric with the both-null case, matches selectRing's own early-return path)");
    }

    // ---- selectChain: substituent chain hanging off a ring ----
    {
        DocumentState doc;
        // Build a benzene-like 6-ring, then a 3-carbon chain off one ring atom.
        QList<double> ring;
        const double kPi = 3.14159265358979323846;
        for (int i = 0; i < 6; ++i) {
            double angle = (kPi / 3.0) * i;
            ring.append(std::cos(angle));
            ring.append(std::sin(angle));
        }
        doc.addRing(ring, true);
        QList<AtomId> ringAtoms = doc.molecule().atomIds();
        AtomId attachPoint = ringAtoms[0];

        AtomId chain1 = doc.molecule().addAtom(QStringLiteral("C"), 3.0, 0.0);
        doc.molecule().addBond(attachPoint, chain1, 1);
        AtomId chain2 = doc.molecule().addAtom(QStringLiteral("C"), 4.0, 0.0);
        doc.molecule().addBond(chain1, chain2, 1);
        AtomId chain3 = doc.molecule().addAtom(QStringLiteral("C"), 5.0, 0.0);
        doc.molecule().addBond(chain2, chain3, 1);

        doc.selectChain(chain3, -1);
        CHECK(doc.selection().atoms.contains(chain1) && doc.selection().atoms.contains(chain2)
              && doc.selection().atoms.contains(chain3),
              "selectChain from the terminal atom selects the whole 3-carbon substituent");
        CHECK(doc.selection().atoms.contains(attachPoint),
              "selectChain includes the ring atom it stops at (included but not expanded past)");
        int ringAtomsSelected = 0;
        for (AtomId id : ringAtoms) if (doc.selection().atoms.contains(id)) ++ringAtomsSelected;
        CHECK(ringAtomsSelected == 1,
              "selectChain does NOT cross into the ring -- only the single attachment atom is "
              "included, none of the other 5 ring atoms");
    }

    // ---- selectChain: stops at a branch point (3+ heavy substituents) ----
    {
        DocumentState doc;
        AtomId center = doc.addAtom(QStringLiteral("C"), 0, 0);
        AtomId branch1 = doc.addAtom(QStringLiteral("C"), 1, 0);
        AtomId branch2 = doc.addAtom(QStringLiteral("C"), -1, 1);
        AtomId branch3 = doc.addAtom(QStringLiteral("C"), -1, -1);
        doc.addBond(center, branch1, 1);
        doc.addBond(center, branch2, 1);
        doc.addBond(center, branch3, 1);
        AtomId beyond = doc.addAtom(QStringLiteral("C"), 2, 0);
        doc.addBond(branch1, beyond, 1);

        doc.selectChain(beyond, -1);
        CHECK(doc.selection().atoms.contains(branch1) && doc.selection().atoms.contains(beyond),
              "selectChain: the terminal atom and its immediate branch-arm neighbor are selected");
        CHECK(doc.selection().atoms.contains(center),
              "selectChain: the branch-point atom itself is included (added before the "
              ">2-heavy-neighbors check stops further expansion)");
        CHECK(!doc.selection().atoms.contains(branch2) && !doc.selection().atoms.contains(branch3),
              "selectChain: the OTHER two branches are NOT crossed into -- branch point stops "
              "expansion unconditionally, even though center is reached (not the seed itself)");
    }

    // ---- selectChain: both-null clears; invalid bondId leaves selection untouched ----
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0, 0);
        doc.selectAtom(a);
        doc.selectChain(-1, -1);
        CHECK(doc.selection().isEmpty(), "selectChain(-1, -1) clears the selection");

        doc.selectAtom(a);
        doc.selectChain(-1, 9999);
        CHECK(doc.selection().atoms.contains(a),
              "selectChain(-1, invalid bond id) leaves the selection untouched, matching the "
              "real code's `else { return }` inside the bond-lookup branch");
    }
}

static void test_templateLibraryNameLists() {
    std::printf("--- Test: TemplateLibrary name-order and <group> metadata ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    // functionalGroupNames(): file insertion order, unsorted -- fg.sdf's first
    // record is "Ac", last is "Ts" (confirmed by direct file inspection and by
    // what the live JS-engine picker already showed before this port).
    {
        QStringList names = lib.functionalGroupNames();
        CHECK(names.size() == 62, "functionalGroupNames() returns all 62 fg.sdf records");
        CHECK(!names.isEmpty() && names.first() == QStringLiteral("Ac"),
              "functionalGroupNames() first entry is \"Ac\" (real file order)");
        CHECK(!names.isEmpty() && names.last() == QStringLiteral("Ts"),
              "functionalGroupNames() last entry is \"Ts\" (real file order)");
    }

    // libraryTemplateNames(): library.sdf has 276 raw records, but only 250 distinct
    // names -- 26 records are genuine same-name duplicates in the real bundled file
    // (e.g. "alpha-D-Allopyranose", "Bicyclo[4-1-1]octane", "Ring5" each appear 2-4
    // times; confirmed by direct inspection, not a loading bug -- the real JS engine's
    // Object.keys()-based store would collapse them identically). Of those 250, one
    // more ("Chlorophyll A_dative") is a V3000-format record using a dative/coordinate
    // bond that this Indigo build's indigoLoadMoleculeFromString rejects outright even
    // in isolation -- a real, narrow Indigo parsing limitation unrelated to this
    // loader, out of scope to chase further for one exotic organometallic template.
    // 250 - 1 = 249 is the real, currently achievable count.
    {
        QStringList names = lib.libraryTemplateNames();
        CHECK(names.size() == 249, "libraryTemplateNames() returns 249 of library.sdf's 250 distinct names");
        CHECK(!names.isEmpty() && names.first() == QStringLiteral("alpha-D-Allopyranose"),
              "libraryTemplateNames() first entry is \"alpha-D-Allopyranose\" (real file order)");
        CHECK(!names.isEmpty() && names.last() == QStringLiteral("Phenylalanine mustard"),
              "libraryTemplateNames() last entry is \"Phenylalanine mustard\" (real file order)");
    }

    // saltOrSolventNames(): 135 records, real first/last by file order.
    {
        QStringList names = lib.saltOrSolventNames();
        CHECK(names.size() == 135, "saltOrSolventNames() returns all 135 salts-and-solvents.sdf records");
        CHECK(!names.isEmpty() && names.first() == QStringLiteral("acetic acid"),
              "saltOrSolventNames() first entry is \"acetic acid\" (real file order)");
        CHECK(!names.isEmpty() && names.last() == QStringLiteral("Ammonium-sulfate"),
              "saltOrSolventNames() last entry is \"Ammonium-sulfate\" (real file order)");
    }

    // functionalGroupGroup(): real fallback value when no <group> is present.
    CHECK(lib.functionalGroupGroup(QStringLiteral("Ac")) == QStringLiteral("Functional Groups"),
          "functionalGroupGroup(\"Ac\") falls back to \"Functional Groups\" (fg.sdf has no <group> field)");
    CHECK(lib.functionalGroupGroup(QStringLiteral("NotARealName")) == QStringLiteral("Functional Groups"),
          "functionalGroupGroup() on an unknown name also falls back to \"Functional Groups\"");

    // libraryTemplateGroup(): real, diverse <group> values, at least 3 exact pairs asserted.
    CHECK(lib.libraryTemplateGroup(QStringLiteral("alpha-D-Allopyranose")) == QStringLiteral("alpha-D-Sugars"),
          "libraryTemplateGroup(\"alpha-D-Allopyranose\") == \"alpha-D-Sugars\"");
    CHECK(lib.libraryTemplateGroup(QStringLiteral("Cyclopenta-1,3-diene")) == QStringLiteral("Aromatics"),
          "libraryTemplateGroup(\"Cyclopenta-1,3-diene\") == \"Aromatics\"");
    CHECK(lib.libraryTemplateGroup(QStringLiteral("Bicyclo[1-1-1]pentane")) == QStringLiteral("Bicycles"),
          "libraryTemplateGroup(\"Bicyclo[1-1-1]pentane\") == \"Bicycles\"");
    CHECK(lib.libraryTemplateGroup(QStringLiteral("NotARealName")) == QStringLiteral("Templates"),
          "libraryTemplateGroup() on an unknown name falls back to \"Templates\"");
}

static void test_alignAndDistributeAtoms() {
    std::printf("--- Test: alignAtoms / distributeAtoms ---\n");

    // alignAtoms("left"): target is the real minimum x; y untouched.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 1.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 1.0, 2.0);   // real minimum x
        AtomId c = doc.addAtom(QStringLiteral("C"), 3.0, 3.0);
        doc.selectAtom(a); doc.addAtomToSelection(b); doc.addAtomToSelection(c);
        doc.alignAtoms(QStringLiteral("left"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 1.0 && y == 1.0, "alignAtoms(left): atom a moved to real min x, y untouched");
        doc.molecule().atomPos(b, x, y);
        CHECK(x == 1.0 && y == 2.0, "alignAtoms(left): atom b (already at min) stays put");
        doc.molecule().atomPos(c, x, y);
        CHECK(x == 1.0 && y == 3.0, "alignAtoms(left): atom c moved to real min x, y untouched");
        doc.undo();
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 5.0 && y == 1.0, "alignAtoms(left): undo restores atom a's original position");
    }

    // alignAtoms("right"): target is the real maximum x.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 9.0, 0.0);   // real maximum x
        doc.selectAtom(a); doc.addAtomToSelection(b);
        doc.alignAtoms(QStringLiteral("right"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 9.0, "alignAtoms(right): atom a moved to real max x");
    }

    // alignAtoms("centerH"): target is the real average x.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 10.0, 0.0);
        doc.selectAtom(a); doc.addAtomToSelection(b);
        doc.alignAtoms(QStringLiteral("centerH"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 5.0, "alignAtoms(centerH): atom a moved to the real average x (5.0)");
        doc.molecule().atomPos(b, x, y);
        CHECK(x == 5.0, "alignAtoms(centerH): atom b moved to the real average x (5.0)");
    }

    // alignAtoms("top"): vertical axis, real minimum y.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 5.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 0.0, 1.0);   // real minimum y
        doc.selectAtom(a); doc.addAtomToSelection(b);
        doc.alignAtoms(QStringLiteral("top"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(y == 1.0, "alignAtoms(top): atom a moved to real min y");
    }

    // alignAtoms: guard, <2 selected atoms is a no-op.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 1.0);
        doc.selectAtom(a);
        doc.alignAtoms(QStringLiteral("left"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 5.0 && y == 1.0, "alignAtoms: <2 selected atoms is a no-op");
    }

    // alignAtoms: an atom present but NOT selected is untouched.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId unselected = doc.addAtom(QStringLiteral("C"), 42.0, 42.0);
        doc.selectAtom(a); doc.addAtomToSelection(b);
        doc.alignAtoms(QStringLiteral("left"));
        double x = 0, y = 0;
        doc.molecule().atomPos(unselected, x, y);
        CHECK(x == 42.0 && y == 42.0, "alignAtoms: an unselected atom is left completely untouched");
    }

    // distributeAtoms("horizontal"): 4 atoms at irregular x, evenly spaced after.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId c = doc.addAtom(QStringLiteral("C"), 8.0, 0.0);
        AtomId d = doc.addAtom(QStringLiteral("C"), 9.0, 0.0);
        doc.selectAtom(a); doc.addAtomToSelection(b); doc.addAtomToSelection(c); doc.addAtomToSelection(d);
        doc.distributeAtoms(QStringLiteral("horizontal"));
        // real step = (9-0)/(4-1) = 3.0 -> expected x: 0, 3, 6, 9
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 0.0, "distributeAtoms(horizontal): min-x atom stays at 0.0");
        doc.molecule().atomPos(b, x, y);
        CHECK(x == 3.0, "distributeAtoms(horizontal): 2nd-sorted atom lands at min+step (3.0)");
        doc.molecule().atomPos(c, x, y);
        CHECK(x == 6.0, "distributeAtoms(horizontal): 3rd-sorted atom lands at min+2*step (6.0)");
        doc.molecule().atomPos(d, x, y);
        CHECK(x == 9.0, "distributeAtoms(horizontal): max-x atom stays at 9.0");
        doc.undo();
        doc.molecule().atomPos(b, x, y);
        CHECK(x == 1.0, "distributeAtoms(horizontal): undo restores original x positions");
    }

    // distributeAtoms("vertical"): same, on the y axis.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 0.0, 2.0);
        AtomId c = doc.addAtom(QStringLiteral("C"), 0.0, 4.0);
        doc.selectAtom(a); doc.addAtomToSelection(b); doc.addAtomToSelection(c);
        doc.distributeAtoms(QStringLiteral("vertical"));
        double x = 0, y = 0;
        doc.molecule().atomPos(b, x, y);
        CHECK(y == 2.0, "distributeAtoms(vertical): already-evenly-spaced middle atom stays put");
    }

    // distributeAtoms: guard, <3 selected atoms is a no-op.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 1.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 9.0, 1.0);
        doc.selectAtom(a); doc.addAtomToSelection(b);
        doc.distributeAtoms(QStringLiteral("horizontal"));
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 5.0, "distributeAtoms: <3 selected atoms is a no-op");
    }
}

static void test_setShowExplicitH() {
    std::printf("--- Test: setShowExplicitH ---\n");
    DocumentState doc;
    CHECK(doc.showExplicitH() == false, "showExplicitH defaults to false");
    doc.setShowExplicitH(true);
    CHECK(doc.showExplicitH() == true, "setShowExplicitH(true) flips the flag");
    doc.setShowExplicitH(false);
    CHECK(doc.showExplicitH() == false, "setShowExplicitH(false) flips it back");
    CHECK(!doc.canUndo(), "toggling showExplicitH does not push an undo entry (not an edit)");
}

static void test_layoutSelectedChain() {
    std::printf("--- Test: layoutSelectedChain ---\n");
    const double kEps = 1e-6;

    // Happy path: Anchor-A-B-C, select {A,B,C} only.
    {
        DocumentState doc;
        AtomId anchor = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 2.0, 0.0);
        AtomId c = doc.addAtom(QStringLiteral("C"), 3.0, 0.0);
        doc.addBond(anchor, a, 1);
        doc.addBond(a, b, 1);
        doc.addBond(b, c, 1);
        doc.selectAtom(a);
        doc.addAtomToSelection(b);
        doc.addAtomToSelection(c);

        double anchorXBefore = 0, anchorYBefore = 0;
        doc.molecule().atomPos(anchor, anchorXBefore, anchorYBefore);

        doc.layoutSelectedChain();

        double anchorXAfter = 0, anchorYAfter = 0;
        doc.molecule().atomPos(anchor, anchorXAfter, anchorYAfter);
        CHECK(anchorXAfter == anchorXBefore && anchorYAfter == anchorYBefore,
              "layoutSelectedChain: the anchor atom's own position is completely untouched");

        double ax = 0, ay = 0, bx = 0, by = 0, cx = 0, cy = 0;
        doc.molecule().atomPos(a, ax, ay);
        doc.molecule().atomPos(b, bx, by);
        doc.molecule().atomPos(c, cx, cy);

        double distAnchorA = std::sqrt((ax - anchorXAfter) * (ax - anchorXAfter) + (ay - anchorYAfter) * (ay - anchorYAfter));
        double distAB = std::sqrt((bx - ax) * (bx - ax) + (by - ay) * (by - ay));
        double distBC = std::sqrt((cx - bx) * (cx - bx) + (cy - by) * (cy - by));
        CHECK(std::abs(distAnchorA - 1.5) < kEps, "layoutSelectedChain: anchor-to-A distance is the bond length (1.5)");
        CHECK(std::abs(distAB - 1.5) < kEps, "layoutSelectedChain: A-to-B distance is the bond length (1.5)");
        CHECK(std::abs(distBC - 1.5) < kEps, "layoutSelectedChain: B-to-C distance is the bond length (1.5)");
        CHECK(std::abs(ax - 1.0) > kEps || std::abs(ay - 0.0) > kEps,
              "layoutSelectedChain: atom A's position actually changed from its pre-layout position");

        doc.undo();
        double axUndone = 0, ayUndone = 0;
        doc.molecule().atomPos(a, axUndone, ayUndone);
        CHECK(std::abs(axUndone - 1.0) < kEps && std::abs(ayUndone - 0.0) < kEps,
              "layoutSelectedChain: undo restores atom A's original position");
    }

    // Guard: 0 attachment bonds (fully disconnected selected fragment) -> no-op.
    {
        DocumentState doc;
        AtomId a = doc.addAtom(QStringLiteral("C"), 5.0, 5.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 6.0, 5.0);
        doc.addBond(a, b, 1);
        doc.selectAtom(a);
        doc.addAtomToSelection(b);
        bool couldUndoBefore = doc.canUndo();
        doc.layoutSelectedChain();
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 5.0 && y == 5.0, "layoutSelectedChain: 0 attachment bonds is a no-op");
        CHECK(doc.canUndo() == couldUndoBefore, "layoutSelectedChain: no-op pushes no undo entry");
    }

    // Guard: 2 attachment bonds (selection attached to the rest of the structure at two
    // separate points) -> no-op.
    {
        DocumentState doc;
        AtomId anchor1 = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 2.0, 0.0);
        AtomId anchor2 = doc.addAtom(QStringLiteral("C"), 3.0, 0.0);
        doc.addBond(anchor1, a, 1);
        doc.addBond(a, b, 1);
        doc.addBond(b, anchor2, 1);
        doc.selectAtom(a);
        doc.addAtomToSelection(b);
        doc.layoutSelectedChain();
        double x = 0, y = 0;
        doc.molecule().atomPos(a, x, y);
        CHECK(x == 1.0 && y == 0.0, "layoutSelectedChain: 2 attachment bonds is a no-op");
    }

    // Guard: branching within the selection (an atom with 3 internal neighbors) -> no-op.
    {
        DocumentState doc;
        AtomId anchor = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId center = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId branch1 = doc.addAtom(QStringLiteral("C"), 2.0, 1.0);
        AtomId branch2 = doc.addAtom(QStringLiteral("C"), 2.0, -1.0);
        AtomId branch3 = doc.addAtom(QStringLiteral("C"), 2.0, 0.0);
        doc.addBond(anchor, center, 1);
        doc.addBond(center, branch1, 1);
        doc.addBond(center, branch2, 1);
        doc.addBond(center, branch3, 1);
        doc.selectAtom(center);
        doc.addAtomToSelection(branch1);
        doc.addAtomToSelection(branch2);
        doc.addAtomToSelection(branch3);
        doc.layoutSelectedChain();
        double x = 0, y = 0;
        doc.molecule().atomPos(center, x, y);
        CHECK(x == 1.0 && y == 0.0, "layoutSelectedChain: branching within the selection is a no-op");
    }

    // Guard: a disconnected sub-group within the selection -- a valid attached chain
    // (Anchor-A-B) plus a wholly separate, unconnected pair also selected (X-Y, bonded to
    // each other but to nothing else) -> no-op, since the walk from the anchor-adjacent atom
    // can never reach X/Y.
    {
        DocumentState doc;
        AtomId anchor = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 2.0, 0.0);
        AtomId x = doc.addAtom(QStringLiteral("C"), 10.0, 10.0);
        AtomId y = doc.addAtom(QStringLiteral("C"), 11.0, 10.0);
        doc.addBond(anchor, a, 1);
        doc.addBond(a, b, 1);
        doc.addBond(x, y, 1);
        doc.selectAtom(a);
        doc.addAtomToSelection(b);
        doc.addAtomToSelection(x);
        doc.addAtomToSelection(y);
        doc.layoutSelectedChain();
        double ax = 0, ay = 0, xx = 0, xy = 0;
        doc.molecule().atomPos(a, ax, ay);
        doc.molecule().atomPos(x, xx, xy);
        CHECK(ax == 1.0 && ay == 0.0, "layoutSelectedChain: disconnected sub-group is a no-op (chain atom untouched)");
        CHECK(xx == 10.0 && xy == 10.0, "layoutSelectedChain: disconnected sub-group is a no-op (orphan atom untouched)");
    }

    // Documented non-guard: a clean cycle (3-ring A-B-C-A, all selected, anchor bonded only
    // to A) is NOT rejected -- the walk reaches all 3 atoms before ever needing the ring-
    // closing bond. This must actually succeed (positions change), matching real behavior.
    {
        DocumentState doc;
        AtomId anchor = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId a = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId b = doc.addAtom(QStringLiteral("C"), 2.0, 0.5);
        AtomId c = doc.addAtom(QStringLiteral("C"), 2.0, -0.5);
        doc.addBond(anchor, a, 1);
        doc.addBond(a, b, 1);
        doc.addBond(b, c, 1);
        doc.addBond(c, a, 1);
        doc.selectAtom(a);
        doc.addAtomToSelection(b);
        doc.addAtomToSelection(c);
        doc.layoutSelectedChain();
        double ax = 0, ay = 0;
        doc.molecule().atomPos(a, ax, ay);
        CHECK(std::abs(ax - 1.0) > kEps || std::abs(ay - 0.0) > kEps,
              "layoutSelectedChain: a clean ring is NOT rejected -- it lays out like a simple chain");
    }
}

static void test_applyDrawingStyleWithRescaleAffectsExportSize() {
    std::printf("--- Test: applyDrawingStyleWithRescale affects physical/export size ---\n");
    const double kEps = 1e-6;

    DocumentState doc;
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = doc.addAtom(QStringLiteral("C"), 1.5, 0.0);

    DrawingStyle oldStyle = doc.drawingStyle();
    double oldModelToMmScale = oldStyle.modelToMmScale;
    double oldPhysicalBondLengthMm = oldStyle.bondLength * oldModelToMmScale;

    // New draft style: double the bond length, everything else (including
    // modelToMmScale) carried over unchanged, exactly as SettingsService's
    // draft (seeded from the document's real style) would hand it in.
    DrawingStyle newStyle = oldStyle;
    newStyle.bondLength = oldStyle.bondLength * 2.0;

    doc.applyDrawingStyleWithRescale(newStyle);

    double ax, ay, bx, by;
    doc.molecule().atomPos(a, ax, ay);
    doc.molecule().atomPos(b, bx, by);
    double newDist = std::hypot(bx - ax, by - ay);
    CHECK(std::abs(newDist - 3.0) < kEps, "atom coordinates rescaled 2x on bond-length change");

    CHECK(std::abs(doc.drawingStyle().modelToMmScale - oldModelToMmScale) < kEps,
          "modelToMmScale stays fixed across the rescale (not re-captured)");

    double newPhysicalBondLengthMm = doc.drawingStyle().bondLength * doc.drawingStyle().modelToMmScale;
    CHECK(std::abs(newPhysicalBondLengthMm - oldPhysicalBondLengthMm * 2.0) < kEps,
          "doubling bond length doubles the final physical/export bond length");

    CHECK(doc.canUndo(), "rescale + style change pushes exactly one undo entry");
    doc.undo();
    doc.molecule().atomPos(a, ax, ay);
    doc.molecule().atomPos(b, bx, by);
    double restoredDist = std::hypot(bx - ax, by - ay);
    CHECK(std::abs(restoredDist - 1.5) < kEps, "undo restores original atom coordinates");
    CHECK(std::abs(doc.drawingStyle().bondLength - oldStyle.bondLength) < kEps,
          "undo restores original bondLength");
    CHECK(std::abs(doc.drawingStyle().modelToMmScale - oldModelToMmScale) < kEps,
          "undo restores original modelToMmScale");
}

static void test_selectSubstructureMatches() {
    std::printf("--- Test: selectSubstructureMatches ---\n");

    DocumentState doc;
    AtomId a = doc.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = doc.addAtom(QStringLiteral("C"), 1.0, 0.0);
    AtomId c = doc.addAtom(QStringLiteral("C"), 2.0, 0.0);
    AtomId d = doc.addAtom(QStringLiteral("O"), 5.0, 5.0); // unconnected
    doc.addBond(a, b, 1);
    doc.addBond(b, c, 1);

    // Determine the real 1-indexed Indigo positions for a and b, mirroring how a genuine
    // caller would have gotten them from IndigoService::substructureSearch's own output.
    QList<AtomId> order = doc.molecule().atomIdsInIndigoOrder();
    int idxA = order.indexOf(a) + 1;
    int idxB = order.indexOf(b) + 1;
    int idxD = order.indexOf(d) + 1;

    // A single match selecting the two connected atoms -> both atoms + the bond between
    // them selected.
    QString matchesAB = QStringLiteral("{\"matches\":[[%1,%2]]}").arg(idxA).arg(idxB);
    doc.selectSubstructureMatches(matchesAB);
    CHECK(doc.selection().atoms.size() == 2 && doc.selection().atoms.contains(a) && doc.selection().atoms.contains(b),
          "selectSubstructureMatches: selects both matched atoms");
    CHECK(doc.selection().bonds.size() == 1, "selectSubstructureMatches: selects the bond between two matched atoms");

    // A match on the lone unconnected atom -> that atom selected, no spurious bond, and the
    // old selection is fully replaced (a/b no longer selected).
    QString matchesD = QStringLiteral("{\"matches\":[[%1]]}").arg(idxD);
    doc.selectSubstructureMatches(matchesD);
    CHECK(doc.selection().atoms.size() == 1 && doc.selection().atoms.contains(d),
          "selectSubstructureMatches: replaces the old selection with the new match");
    CHECK(doc.selection().bonds.isEmpty(), "selectSubstructureMatches: no bond for a lone unconnected match");

    // Out-of-range index -> silently ignored, no crash.
    doc.selectSubstructureMatches(QStringLiteral("{\"matches\":[[999]]}"));
    CHECK(doc.selection().atoms.isEmpty(), "selectSubstructureMatches: out-of-range index selects nothing");

    // Malformed JSON -> selection cleared, no crash.
    doc.selectSubstructureMatches(matchesAB);
    doc.selectSubstructureMatches(QStringLiteral("not json"));
    CHECK(doc.selection().isEmpty(), "selectSubstructureMatches: malformed JSON clears the selection, no crash");

    // Empty matches array -> selection cleared (the real "clear matches" UI case).
    doc.selectSubstructureMatches(matchesAB);
    doc.selectSubstructureMatches(QStringLiteral("{\"matches\":[]}"));
    CHECK(doc.selection().isEmpty(), "selectSubstructureMatches: empty matches array clears the selection");
}

int main() {
    test_selectionStateBasics();
    test_editCommandBasics();
    test_documentStateUndoRedo();
    test_executeCommandReentrancyGuard();
    test_undoReentrancyGuard();
    test_selectionReconciliationRemovesFullyStaleId();
    test_selectionReconciliationKeepsCollidingRealAtom();
    test_selectionReconciliationCleansUpAfterDeleteSelectionEntities();
    test_selectionReconciliationRemovesStaleBondId();
    test_selectionReconciliationRemovesStaleRxnArrowId();
    test_selectionReconciliationRemovesStaleRxnPlusId();
    test_selectionReconciliationRemovesStaleMultitailArrowId();
    test_documentStateHistoryCap();
    test_documentStateSelection();
    test_documentStateEditingOperations();
    test_discreteTransforms();
    test_imageTransforms();
    test_liveDragMove();
    test_liveDragRotate();
    test_liveDragScale();
    test_dragStateResetOnUndoRedo();
    test_addRing();
    test_addChain();
    test_importReaction();
    test_importReactionRxnPlusGapCentering();
    test_importReactionRxnPlusGapCenteringSingleAtomFloor();
    test_importReactionRecentersOnPageOrigin();
    test_importReactionSingleReactantSingleProductRecentering();
    test_importReactionWithCatalyst();
    test_importReactionWithMultiAtomCatalyst();
    test_importReactionWithTallReactant();
    test_insertFunctionalGroup();
    test_insertFunctionalGroupSGroupCleanupOnUndo();
    test_insertFunctionalGroupAllTemplatesConnectToTarget();
    test_insertFunctionalGroupConh2OneNeighborConnects();
    test_getFunctionalGroupPreviewFreeFloating();
    test_getFunctionalGroupPreviewGraft();
    test_graftAngleOrientation_oneNeighbor();
    test_graftAngleOrientation_twoNeighbors();
    test_graftAngleOrientation_ringAtomTarget();
    test_graftAngleOrientation_noRotationCase();
    test_graftAngleOrientation_noTemplateInternalBonds();
    test_graftAngleOrientation_threeNeighbors();
    test_graftAngleOrientation_multiBondAttachAtom();
    test_graftAngleOrientation_undoRedoRoundTrip();
    test_insertLibraryTemplateFused();
    test_toggleSgroupExpanded();
    test_rxnArrowLifecycle();
    test_rxnPlusCommands();
    test_multitailArrowCommands();
    test_setStereoFlags();
    test_rgroupCommands();
    test_textCommands();
    test_addBracketSelection();
    test_imageCommands();
    test_deserializeMol();
    test_deserializeMolCenterOnPage();
    test_insertFunctionalGroupLabel();
    test_renameSGroup();
    test_clearCanvas();
    test_loadBenzene();
    test_documentStateBioSequenceViewForwarding();
    test_setMoleculeName();
    test_copySelection();
    test_insertStructureAt();
    test_insertStructureAtClampsToPageBounds();
    test_copyPasteRoundTrip();
    test_deleteSelectionEntitiesBondedPair();
    test_deleteSelectionEntitiesWholePill();
    test_deleteSelectionEntitiesPartialMember();
    test_deleteSelectionEntitiesMixedTypes();
    test_deleteSelectionEntitiesRedoAfterUndo();
    test_deleteSelectionEntitiesWholePillRedoAfterUndo();
    test_cutSelection();
    test_setStereoDescriptorsValid();
    test_setStereoDescriptorsReuseOrderRegression();
    test_setStereoDescriptorsResetsOnEachCall();
    test_setStereoDescriptorsMalformedJsonIsNoOp();
    test_setStereoDescriptorsOutOfRangeIndexSkipped();
    test_setCheckIssuesValid();
    test_setCheckIssuesDefaultTargetIsAtom();
    test_setCheckIssuesResetsOnEachCall();
    test_setCheckIssuesMalformedJsonIsNoOp();
    test_setCheckIssuesReuseOrderRegression();
    test_changeBondTypeAndStereo();
    test_selectByRectAndLasso();
    test_selectFragmentRingChain();
    test_templateLibraryNameLists();
    test_alignAndDistributeAtoms();
    test_setShowExplicitH();
    test_layoutSelectedChain();
    test_selectSubstructureMatches();
    test_applyDrawingStyleWithRescaleAffectsExportSize();
    std::printf("Summary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}

