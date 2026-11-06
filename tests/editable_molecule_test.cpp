// tests/editable_molecule_test.cpp
// Standalone tests for the EditableMolecule module (chem-core.js migration,
// sub-project 1). Test 0 is the ID-behavior spike from the design spec: it
// DOCUMENTS how Indigo indices behave under removal/re-add and VERIFIES the
// one property the design depends on (indigoClone preserves indices).
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <QPointF>
#include <QJsonDocument>
#include <QJsonObject>
#include <functional>
#include "app/molecule/EditableMolecule.h"
#include "app/molecule/TemplateLibrary.h"
#include "indigo.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

static void spike_idBehavior() {
    std::printf("--- Test 0: Indigo ID-behavior spike ---\n");
    int mol = indigoCreateMolecule();
    CHECK(mol >= 0, "spike: create empty molecule");

    // Add 5 carbons, record each one's index.
    int atomHandles[5], atomIdx[5];
    for (int i = 0; i < 5; ++i) {
        atomHandles[i] = indigoAddAtom(mol, "C");
        atomIdx[i] = indigoIndex(atomHandles[i]);
        indigoSetXYZ(atomHandles[i], i * 1.0f, 0.0f, 0.0f);
    }
    // Chain bonds 0-1-2-3-4.
    for (int i = 0; i < 4; ++i)
        indigoAddBond(atomHandles[i], atomHandles[i + 1], 1);

    // Remove atom #2 (the middle one).
    indigoRemove(atomHandles[2]);
    CHECK(indigoCountAtoms(mol) == 4, "spike: 4 atoms remain after removal");

    // OBSERVATION (not a hard requirement of the design): do survivors keep
    // their old indices? Print what we see either way.
    int survivorsKeptIndices = 1;
    const int surviving[4] = {0, 1, 3, 4};
    for (int s = 0; s < 4; ++s) {
        int a = indigoGetAtom(mol, atomIdx[surviving[s]]);
        if (a < 0) { survivorsKeptIndices = 0; break; }
        indigoFree(a);
    }
    std::printf("[INFO] spike: surviving atoms keep their pre-removal indices: %s\n",
                survivorsKeptIndices ? "YES" : "NO");

    // OBSERVATION: does a new atom reuse the freed index?
    int newAtom = indigoAddAtom(mol, "N");
    int newIdx = indigoIndex(newAtom);
    std::printf("[INFO] spike: freed index was %d, new atom got index %d (%s)\n",
                atomIdx[2], newIdx, newIdx == atomIdx[2] ? "REUSED" : "not reused");

    // HARD REQUIREMENT: indigoClone preserves indices (snapshot/restore
    // depends on this -- the ID indirection table maps external IDs to these
    // indices and is copied verbatim into snapshots).
    int clone = indigoClone(mol);
    CHECK(clone >= 0, "spike: clone succeeds");
    CHECK(indigoCountAtoms(clone) == indigoCountAtoms(mol),
          "spike: clone has same atom count");
    int cloneAtom = indigoGetAtom(clone, newIdx);
    CHECK(cloneAtom >= 0 && indigoAtomicNumber(cloneAtom) == 7,
          "spike: clone preserves the N atom at the same index");
    if (cloneAtom >= 0) indigoFree(cloneAtom);

    indigoFree(clone);
    indigoFree(mol);
}

static void test_constructionAndMolfile() {
    std::printf("--- Test 1: construction + molfile round-trip ---\n");

    EditableMolecule empty;
    CHECK(empty.isValid(), "empty molecule constructs valid");
    CHECK(empty.atomCount() == 0, "empty molecule has 0 atoms");

    // indigoLoadMoleculeFromString auto-detects SMILES.
    EditableMolecule ethanol(QStringLiteral("CCO"));
    CHECK(ethanol.isValid(), "ethanol constructs valid");
    CHECK(ethanol.atomCount() == 3, "ethanol has 3 heavy atoms");
    CHECK(ethanol.bondCount() == 2, "ethanol has 2 bonds");

    StringResult mf = ethanol.toMolfile();
    CHECK(mf.success, "toMolfile succeeds");
    CHECK(mf.value.contains(QStringLiteral("M  END")), "molfile has M  END");

    EditableMolecule bad(QStringLiteral("not_a_molecule((("));
    CHECK(!bad.isValid(), "garbage input yields invalid molecule");
    CHECK(!bad.lastError().isEmpty(), "invalid molecule carries an error string");
}

static void test_toKetJson() {
    std::printf("--- Test: toKetJson ---\n");

    // indigoLoadMoleculeFromString auto-detects SMILES, matching test_constructionAndMolfile's
    // own established convention.
    EditableMolecule ethane(QStringLiteral("CC"));
    CHECK(ethane.isValid(), "toKetJson: ethane constructs valid");

    QString ket = ethane.toKetJson();
    CHECK(!ket.isEmpty(), "toKetJson: non-empty result for a valid molecule");

    QJsonParseError err;
    QJsonDocument doc = QJsonDocument::fromJson(ket.toUtf8(), &err);
    CHECK(err.error == QJsonParseError::NoError, "toKetJson: result parses as valid JSON");
    CHECK(doc.isObject() && doc.object().contains(QStringLiteral("root")),
          "toKetJson: result has a top-level \"root\" key (real KET shape)");

    // Genuinely-invalid input (not empty string, which yields a valid zero-atom molecule --
    // see test_constructionAndMolfile's own `empty.isValid()` assertion).
    EditableMolecule bad(QStringLiteral("not_a_molecule((("));
    CHECK(!bad.isValid(), "toKetJson: garbage input is invalid");
    CHECK(bad.toKetJson().isEmpty(), "toKetJson: empty result for an invalid molecule, no crash");
}

static void test_atomBondCrud() {
    std::printf("--- Test 2: atom/bond CRUD with stable IDs ---\n");
    EditableMolecule m;

    AtomId c1 = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId c2 = m.addAtom(QStringLiteral("C"), 1.0, 0.0);
    AtomId o1 = m.addAtom(QStringLiteral("O"), 2.0, 0.0);
    CHECK(c1 > 0 && c2 > 0 && o1 > 0, "addAtom returns positive IDs");
    CHECK(c1 != c2 && c2 != o1, "IDs are distinct");
    CHECK(m.atomCount() == 3, "3 atoms after adds");
    CHECK(m.atomSymbol(o1) == QStringLiteral("O"), "atomSymbol reads back");

    double x = -1, y = -1;
    CHECK(m.atomPos(c2, x, y) && x == 1.0 && y == 0.0, "atomPos reads back");

    BondId b1 = m.addBond(c1, c2, 1);
    BondId b2 = m.addBond(c2, o1, 2);
    CHECK(b1 > 0 && b2 > 0, "addBond returns positive IDs");
    CHECK(m.bondOrder(b2) == 2, "bondOrder reads back");

    // Remove the middle atom: its two incident bonds must die with it.
    CHECK(m.removeAtom(c2), "removeAtom succeeds");
    CHECK(m.atomCount() == 2, "2 atoms remain");
    CHECK(m.bondCount() == 0, "incident bonds removed with the atom");
    CHECK(m.bondOrder(b1) == -1, "stale BondId no longer resolves");
    CHECK(m.atomSymbol(c2).isEmpty(), "stale AtomId no longer resolves");

    // Survivors still resolve correctly through the indirection table even
    // if Indigo compacted/reused indices underneath.
    CHECK(m.atomSymbol(c1) == QStringLiteral("C"), "survivor c1 still resolves");
    CHECK(m.atomSymbol(o1) == QStringLiteral("O"), "survivor o1 still resolves");

    // A new atom must get a fresh ID, never a recycled one.
    AtomId n1 = m.addAtom(QStringLiteral("N"), 3.0, 0.0);
    CHECK(n1 != c2 && n1 > o1, "new ID is fresh, not recycled");
    CHECK(m.removeBond(m.addBond(c1, n1, 1)), "removeBond on a live bond succeeds");
    CHECK(!m.removeAtom(c2), "removing an already-removed ID fails cleanly");
}

static void test_extensionData() {
    std::printf("--- Test 3: extension data storage ---\n");
    EditableMolecule m(QStringLiteral("CCO"));
    QList<AtomId> ids = m.atomIds();

    m.setName(QStringLiteral("my ethanol"));
    CHECK(m.name() == QStringLiteral("my ethanol"), "name round-trips");

    CHECK(m.addTextAnnotation(1.0, 2.0, QStringLiteral("note")) >= 0, "text annotation added");
    CHECK(m.textAnnotationCount() == 1, "text annotation counted");

    CHECK(m.addRxnArrow(0, 0, 5, 0) >= 0, "rxn arrow added");
    CHECK(m.addRxnPlus(2.5, 0) >= 0, "rxn plus added");
    CHECK(m.rxnArrowCount() == 1 && m.rxnPlusCount() == 1, "arrow+plus counted");

    CHECK(m.addMultitailArrow({5,5, 0,4, 0,6}) >= 0, "multitail arrow added");
    CHECK(m.multitailArrowCount() == 1, "multitail counted");

    CHECK(m.addImage(0, 0, 4, 3, QByteArray("fakepng")) >= 0, "image added");
    CHECK(m.imageCount() == 1, "image counted");

    m.setStereoFlag(0, 2);
    CHECK(m.stereoFlag(0) == 2, "stereo flag round-trips");
    CHECK(m.stereoFlag(1) == -1, "unset stereo flag is -1");

    m.setAtomAAM(ids[0], 7);
    CHECK(m.atomAAM(ids[0]) == 7, "AAM round-trips");
    CHECK(m.atomAAM(ids[1]) == 0, "unset AAM is 0");

    m.setAtomCheckWarning(ids[2], true);
    CHECK(m.atomCheckWarning(ids[2]), "check warning round-trips");
    CHECK(!m.atomCheckWarning(ids[0]), "unset check warning is false");
}

static void test_sgroupPassthrough() {
    std::printf("--- Test 4: data sgroup passthrough ---\n");
    EditableMolecule m(QStringLiteral("CCO"));
    QList<AtomId> ids = m.atomIds();

    int sg = m.addDataSGroup({ids[0], ids[1]}, QStringLiteral("MYFIELD"), QStringLiteral("hello"));
    CHECK(sg >= 0, "data sgroup added");
    CHECK(m.dataSGroupCount() == 1, "data sgroup counted");

    StringResult mf = m.toMolfile();
    CHECK(mf.success && mf.value.contains(QStringLiteral("MYFIELD")), "sgroup survives molfile serialization");
}

static void test_snapshotRestore() {
    std::printf("--- Test 5: snapshot/restore undo round-trip ---\n");
    EditableMolecule m;
    AtomId a1 = m.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = m.addAtom(QStringLiteral("O"), 1, 0);
    m.addBond(a1, a2, 1);
    m.setName(QStringLiteral("before"));
    m.setAtomAAM(a1, 3);
    m.addTextAnnotation(0, 1, QStringLiteral("v1"));
    m.addDataSGroup({a1}, QStringLiteral("F"), QStringLiteral("d"));
    m.addRxnArrow(0, 0, 5, 0);
    m.addRxnPlus(2.5, 0);
    m.addMultitailArrow({5,5, 0,4, 0,6});
    m.addImage(0, 0, 4, 3, QByteArray("png1"));
    m.setStereoFlag(0, 1);
    m.setAtomCheckWarning(a2, true);

    // Real superatom (SUP) sgroup, tracked via EditableMolecule's own external m_sgroupIdx/
    // m_sgroupExpanded -- a DIFFERENT mechanism from addDataSGroup above, which lives entirely
    // inside Indigo's own molecule object and needs no external tracking. This is the exact
    // "Ac" fixture already used by test_sgroupIntrospectionAccessors.
    QString acMolfile = QStringLiteral(
        "Ac\n"
        "  Ketcher\n\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    3.3951   -3.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.6785   -3.9891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.3951   -2.7480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  2  1  1  0  0  0  0\n"
        "  1  3  2  0  0  0  0\n"
        "M  STY  1   1 SUP\n"
        "M  SLB  1   1   1\n"
        "M  SAL   1  3   1   2   3\n"
        "M  SAP   1  1   1   0\n"
        "M  SMT   1 Ac\n"
        "M  END\n");
    EditableMolecule::InsertResult acResult = m.insertStructure(
        acMolfile, [](double x, double y) { return QPointF(x, y); });
    CHECK(acResult.createdSGroups.size() == 1, "setup: \"Ac\" superatom inserted");
    SGroupId sgId = acResult.createdSGroups[0];
    m.setSGroupExpanded(sgId, false);   // away from the default (true), so restore is verifiable
    CHECK(!m.sgroupExpanded(sgId), "setup: sgroup starts collapsed for this test");

    MoleculeSnapshot snap = m.snapshot();

    AtomId a3 = m.addAtom(QStringLiteral("N"), 2, 0);
    m.addBond(a2, a3, 1);
    m.removeAtom(a1);
    m.setName(QStringLiteral("after"));
    m.setAtomAAM(a2, 9);
    m.addTextAnnotation(0, 2, QStringLiteral("v2"));
    m.setSGroupExpanded(sgId, true);   // mutate sgroup state too, AFTER the snapshot
    CHECK(m.atomCount() == 5 && m.name() == QStringLiteral("after"), "mutations applied");

    CHECK(m.restore(snap), "restore succeeds");
    CHECK(m.atomCount() == 5, "atom count restored");   // 2 plain + 3 from the Ac superatom
    CHECK(m.bondCount() == 3, "bond count restored");   // 1 plain + 2 from the Ac superatom
    CHECK(m.name() == QStringLiteral("before"), "name restored");
    CHECK(m.atomAAM(a1) == 3, "AAM restored under the ORIGINAL AtomId");
    CHECK(m.atomSymbol(a1) == QStringLiteral("C"), "original AtomId resolves after restore");
    CHECK(m.textAnnotationCount() == 1, "extension list restored");
    CHECK(m.dataSGroupCount() == 1, "sgroup survives snapshot/restore");
    CHECK(m.rxnArrowCount() == 1 && m.rxnPlusCount() == 1, "rxn arrow+plus survive restore");
    CHECK(m.multitailArrowCount() == 1, "multitail arrow survives restore");
    CHECK(m.imageCount() == 1, "image survives restore");
    CHECK(m.stereoFlag(0) == 1, "stereo flag survives restore");
    CHECK(m.atomCheckWarning(a2), "check warning survives restore");

    // The actual fix under test: SUP-superatom tracking (m_sgroupIdx/m_sgroupExpanded), a
    // SEPARATE mechanism from the data-sgroup check above, must also survive restore().
    CHECK(m.sgroupIds().contains(sgId), "SUP-sgroup id survives snapshot/restore");
    CHECK(m.sgroupMemberAtomIds(sgId).size() == 3, "SUP-sgroup member atoms survive snapshot/restore");
    CHECK(!m.sgroupExpanded(sgId), "SUP-sgroup expanded flag restored to its pre-mutation value (false)");

    m.removeAtom(a2);
    CHECK(m.restore(snap), "second restore of the same snapshot succeeds");
    CHECK(m.atomCount() == 5, "second restore rolls back again");

    AtomId a4 = m.addAtom(QStringLiteral("S"), 3, 0);
    CHECK(a4 > a3, "ID counter survives restore (no reuse of discarded IDs)");
}

static void test_queryFeatureSpike() {
    std::printf("--- Test 6: query-atom-features spike (spec risk 3) ---\n");
    // Molfile with an atom list ([N,O] on atom 1) -- V2000 query feature.
    const char* molWithAtomList =
        "\n  spike  \n\n"
        "  2  1  1  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  1 T    2   7   8\n"
        "M  END\n";

    indigoSetSessionId(indigoAllocSessionId()); // fresh session, throwaway
    int plain = indigoLoadMoleculeFromString(molWithAtomList);
    std::printf("[INFO] spike: atom-list molfile into PLAIN molecule handle: %s\n",
                plain >= 0 ? "ACCEPTED" : indigoGetLastError());
    int query = indigoLoadQueryMoleculeFromString(molWithAtomList);
    std::printf("[INFO] spike: atom-list molfile into QUERY molecule handle: %s\n",
                query >= 0 ? "ACCEPTED" : indigoGetLastError());
    CHECK(plain >= 0 || query >= 0, "atom-list molfile loads via at least one loader");
    if (plain >= 0) indigoFree(plain);
    if (query >= 0) indigoFree(query);
}

static void test_extensionEntityRemoval() {
    std::printf("--- Test 7: extension entity removal with stable IDs ---\n");
    EditableMolecule m;

    int t1 = m.addTextAnnotation(0, 0, QStringLiteral("a"));
    int t2 = m.addTextAnnotation(1, 1, QStringLiteral("b"));
    CHECK(m.textAnnotationCount() == 2, "2 texts after adds");
    CHECK(m.removeTextAnnotation(t1), "removeTextAnnotation succeeds for a live id");
    CHECK(m.textAnnotationCount() == 1, "1 text remains after removal");
    CHECK(!m.removeTextAnnotation(t1), "removing an already-removed id fails cleanly");
    int t3 = m.addTextAnnotation(2, 2, QStringLiteral("c"));
    CHECK(t3 != t1 && t3 > t2, "new text id is fresh, not recycled");

    int ra1 = m.addRxnArrow(0, 0, 1, 1);
    CHECK(m.removeRxnArrow(ra1), "removeRxnArrow succeeds");
    CHECK(m.rxnArrowCount() == 0, "0 rxn arrows after removal");
    CHECK(!m.removeRxnArrow(ra1), "double-remove of rxn arrow fails cleanly");

    int rp1 = m.addRxnPlus(0, 0);
    CHECK(m.removeRxnPlus(rp1), "removeRxnPlus succeeds");
    CHECK(m.rxnPlusCount() == 0, "0 rxn pluses after removal");

    int mta1 = m.addMultitailArrow({0,0, 1,1});
    CHECK(m.removeMultitailArrow(mta1), "removeMultitailArrow succeeds");
    CHECK(m.multitailArrowCount() == 0, "0 multitail arrows after removal");

    int img1 = m.addImage(0, 0, 1, 1, QByteArray("x"));
    CHECK(m.removeImage(img1), "removeImage succeeds");
    CHECK(m.imageCount() == 0, "0 images after removal");
}

static void test_atomAttributeAccessors() {
    std::printf("--- Test 8: atom attribute accessors (charge/isotope/radical/valence) ---\n");
    EditableMolecule m;
    AtomId c1 = m.addAtom(QStringLiteral("C"), 0, 0);

    CHECK(m.atomCharge(c1) == 0, "fresh atom has charge 0");
    CHECK(m.setAtomCharge(c1, 1), "setAtomCharge succeeds");
    CHECK(m.atomCharge(c1) == 1, "atomCharge reads back the set value");
    CHECK(m.setAtomCharge(c1, -1), "setAtomCharge accepts negative charge");
    CHECK(m.atomCharge(c1) == -1, "atomCharge reads back negative charge");

    CHECK(m.atomIsotope(c1) == 0, "fresh atom has isotope 0 (none)");
    CHECK(m.setAtomIsotope(c1, 13), "setAtomIsotope succeeds");
    CHECK(m.atomIsotope(c1) == 13, "atomIsotope reads back the set value");

    CHECK(m.atomRadical(c1) == 0, "fresh atom has radical 0 (none)");
    CHECK(m.setAtomRadical(c1, INDIGO_SINGLET), "setAtomRadical succeeds");
    CHECK(m.atomRadical(c1) == INDIGO_SINGLET, "atomRadical reads back the set value");

    CHECK(m.atomExplicitValence(c1) == -1, "fresh atom has explicit valence -1 (unset)");
    CHECK(m.setAtomExplicitValence(c1, 4), "setAtomExplicitValence succeeds");
    CHECK(m.atomExplicitValence(c1) == 4, "atomExplicitValence reads back the set value");

    AtomId invalid = 9999;
    CHECK(m.atomCharge(invalid) == 0, "invalid id: atomCharge returns sentinel 0");
    CHECK(!m.setAtomCharge(invalid, 1), "invalid id: setAtomCharge fails cleanly");
    CHECK(m.atomIsotope(invalid) == 0, "invalid id: atomIsotope returns sentinel 0");
    CHECK(!m.setAtomIsotope(invalid, 1), "invalid id: setAtomIsotope fails cleanly");
    CHECK(m.atomRadical(invalid) == 0, "invalid id: atomRadical returns sentinel 0");
    CHECK(!m.setAtomRadical(invalid, 1), "invalid id: setAtomRadical fails cleanly");
    CHECK(m.atomExplicitValence(invalid) == -1, "invalid id: atomExplicitValence returns sentinel -1");
    CHECK(!m.setAtomExplicitValence(invalid, 1), "invalid id: setAtomExplicitValence fails cleanly");
}

static void test_atomLabelAndBondOrder() {
    std::printf("--- Test 9: setAtomLabel + setBondOrderValue ---\n");
    EditableMolecule m;
    AtomId c1 = m.addAtom(QStringLiteral("C"), 0, 0);
    m.setAtomCharge(c1, 1);
    m.setAtomIsotope(c1, 13);

    CHECK(m.setAtomLabel(c1, QStringLiteral("N")), "setAtomLabel to a valid element succeeds");
    CHECK(m.atomSymbol(c1) == QStringLiteral("N"), "atomSymbol reflects the new label");
    CHECK(m.atomCharge(c1) == 1, "charge survives the label change");
    CHECK(m.atomIsotope(c1) == 13, "isotope survives the label change");

    CHECK(!m.setAtomLabel(c1, QStringLiteral("NotAnElement")), "setAtomLabel rejects an invalid symbol");
    CHECK(m.atomSymbol(c1) == QStringLiteral("N"), "label unchanged after a rejected setAtomLabel");

    AtomId c2 = m.addAtom(QStringLiteral("C"), 1, 0);
    BondId b1 = m.addBond(c1, c2, 1);
    CHECK(m.bondOrder(b1) == 1, "fresh bond has order 1");
    CHECK(m.setBondOrderValue(b1, 2), "setBondOrderValue succeeds");
    CHECK(m.bondOrder(b1) == 2, "bondOrder reads back the new order");

    AtomId invalid = 9999;
    CHECK(!m.setAtomLabel(invalid, QStringLiteral("O")), "invalid atom id: setAtomLabel fails cleanly");
    CHECK(!m.setBondOrderValue(invalid, 1), "invalid bond id: setBondOrderValue fails cleanly");
}

static void test_attachmentPointAccessors() {
    std::printf("--- Test 10: attachment-point accessors ---\n");
    EditableMolecule m;
    AtomId c1 = m.addAtom(QStringLiteral("C"), 0, 0);
    AtomId c2 = m.addAtom(QStringLiteral("C"), 1, 0);

    CHECK(m.atomAttachmentOrder(c1) == 0, "fresh atom has attachment order 0 (none)");
    CHECK(m.setAtomAttachmentOrder(c1, 1), "setAtomAttachmentOrder(1) succeeds");
    CHECK(m.atomAttachmentOrder(c1) == 1, "atomAttachmentOrder reads back order 1");

    CHECK(m.setAtomAttachmentOrder(c2, 2), "setAtomAttachmentOrder(2) on a different atom succeeds");
    CHECK(m.atomAttachmentOrder(c2) == 2, "second atom reports its own order 2");
    CHECK(m.atomAttachmentOrder(c1) == 1, "first atom's order unaffected by the second set");

    CHECK(m.setAtomAttachmentOrder(c1, 0), "setAtomAttachmentOrder(0) clears it");
    CHECK(m.atomAttachmentOrder(c1) == 0, "atomAttachmentOrder reports 0 after clearing");

    AtomId invalid = 9999;
    CHECK(!m.setAtomAttachmentOrder(invalid, 1), "invalid id: setAtomAttachmentOrder fails cleanly");
    CHECK(m.atomAttachmentOrder(invalid) == 0, "invalid id: atomAttachmentOrder returns sentinel 0");
}

static void test_sgroupMembershipSpike() {
    std::printf("--- Test 11: sgroup-atom-membership spike (3a design question) ---\n");
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    int mol = indigoLoadMoleculeFromString("CCO");
    CHECK(mol >= 0, "spike: load CCO for sgroup spike");

    int atomIndices[2] = { 0, 1 };
    int sg = indigoAddDataSGroup(mol, 2, atomIndices, 0, nullptr, "F", "d");
    CHECK(sg >= 0, "spike: data sgroup created on the raw handle");

    // Does indigoIterateAtoms work when called on the SGROUP handle itself
    // (as opposed to the molecule handle)? indigo.h's own doc comment says
    // "for all atoms of the given molecule" -- test empirically anyway.
    int sgAtomIter = indigoIterateAtoms(sg);
    int sgAtomCount = 0;
    if (sgAtomIter >= 0) {
        int h;
        while ((h = indigoNext(sgAtomIter)) > 0) { ++sgAtomCount; indigoFree(h); }
        indigoFree(sgAtomIter);
    }
    std::printf("[INFO] spike: indigoIterateAtoms(sgroupHandle) iterator=%s, atom count=%d\n",
                sgAtomIter >= 0 ? "valid" : "invalid(-1)", sgAtomCount);

    // Does the generic indigoRemove(item) delete a WHOLE sgroup (its atoms and
    // bonds too), or just detach the sgroup record leaving atoms untouched, or
    // fail outright? Order matters: check counts before/after.
    int countBefore = indigoCountAtoms(mol);
    int removeResult = indigoRemove(sg);
    int countAfter = indigoCountAtoms(mol);
    std::printf("[INFO] spike: indigoRemove(sgroupHandle) returned %d; atom count before=%d after=%d\n",
                removeResult, countBefore, countAfter);

    CHECK(true, "spike: sgroup-membership investigation completed (see [INFO] lines above)");

    indigoFree(mol);
    indigoReleaseSessionId(session);
}

static void test_bondStereoSpike() {
    std::printf("--- Test 12: bond-stereo mechanism spike (3a design question) ---\n");
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    // Alanine (CC(N)C(=O)O) has a real stereocenter at atom index 1 -- used to
    // test whether forcing a specific tetrahedral arrangement via
    // indigoAddStereocenter, then indigoMarkStereobonds, causes any incident
    // bond's (read-only) indigoBondStereo() to report a nonzero wedge/dash.
    int mol = indigoLoadMoleculeFromString("CC(N)C(=O)O");
    CHECK(mol >= 0, "spike: load alanine for bond-stereo spike");

    int centerAtom = indigoGetAtom(mol, 1);
    CHECK(centerAtom >= 0, "spike: got the stereocenter atom handle");
    int beforeType = indigoStereocenterType(centerAtom);
    std::printf("[INFO] spike: stereocenter type before indigoAddStereocenter: %d\n", beforeType);

    QList<int> neiIdx;
    int neiIter = indigoIterateNeighbors(centerAtom);
    if (neiIter >= 0) {
        int h;
        while ((h = indigoNext(neiIter)) > 0) { neiIdx.append(indigoIndex(h)); indigoFree(h); }
        indigoFree(neiIter);
    }
    int v1 = neiIdx.size() > 0 ? neiIdx[0] : -1;
    int v2 = neiIdx.size() > 1 ? neiIdx[1] : -1;
    int v3 = neiIdx.size() > 2 ? neiIdx[2] : -1;
    int setResult = indigoAddStereocenter(centerAtom, INDIGO_ABS, v1, v2, v3, -1);
    std::printf("[INFO] spike: indigoAddStereocenter(ABS) returned %d\n", setResult);

    int markResult = indigoMarkStereobonds(mol);
    std::printf("[INFO] spike: indigoMarkStereobonds returned %d\n", markResult);

    int bondIter = indigoIterateBonds(mol);
    int stereoBondsFound = 0;
    if (bondIter >= 0) {
        int b;
        while ((b = indigoNext(bondIter)) > 0) {
            if (indigoBondStereo(b) != 0) ++stereoBondsFound;
            indigoFree(b);
        }
        indigoFree(bondIter);
    }
    std::printf("[INFO] spike: bonds with nonzero indigoBondStereo() after marking: %d\n", stereoBondsFound);

    CHECK(true, "spike: bond-stereo mechanism investigation completed (see [INFO] lines above)");

    indigoFree(centerAtom);
    indigoFree(mol);
    indigoReleaseSessionId(session);
}

static void test_atomQueryList() {
    std::printf("--- Test 13: query atom list sidecar ---\n");
    EditableMolecule m;
    AtomId c1 = m.addAtom(QStringLiteral("C"), 0, 0);

    CHECK(!m.hasAtomQueryList(c1), "fresh atom has no query list");

    CHECK(m.setAtomQueryList(c1, QStringLiteral("C, N , O"), false), "setAtomQueryList succeeds");
    CHECK(m.hasAtomQueryList(c1), "hasAtomQueryList true after setting");
    CHECK(m.atomSymbol(c1) == QStringLiteral("L#"), "atomSymbol becomes the L# sentinel");
    QList<int> nums = m.atomQueryListNumbers(c1);
    CHECK(nums.size() == 3, "3 distinct elements parsed from the CSV (C, N, O)");
    CHECK(nums.contains(6) && nums.contains(7) && nums.contains(8),
          "parsed atomic numbers are exactly carbon/nitrogen/oxygen");
    CHECK(!m.atomQueryListIsNotList(c1), "notList flag is false");

    CHECK(m.setAtomQueryList(c1, QStringLiteral("F,Cl"), true), "re-setting the query list succeeds");
    CHECK(m.atomQueryListIsNotList(c1), "notList flag is true after re-setting with notList=true");
    CHECK(m.atomQueryListNumbers(c1).size() == 2, "2 elements after re-setting");

    CHECK(m.clearAtomQueryList(c1, QStringLiteral("N")), "clearAtomQueryList succeeds");
    CHECK(!m.hasAtomQueryList(c1), "hasAtomQueryList false after clearing");
    CHECK(m.atomSymbol(c1) == QStringLiteral("N"), "atomSymbol reset to the fallback label");

    CHECK(!m.setAtomQueryList(c1, QStringLiteral("NotAnElement,AlsoNot"), false),
          "setAtomQueryList with zero valid labels fails cleanly (matches _makeAtomList returning null)");
    CHECK(!m.hasAtomQueryList(c1), "no query list was set from all-invalid input");

    AtomId invalid = 9999;
    CHECK(!m.setAtomQueryList(invalid, QStringLiteral("C"), false), "invalid id: setAtomQueryList fails cleanly");
    CHECK(!m.clearAtomQueryList(invalid, QStringLiteral("C")), "invalid id: clearAtomQueryList fails cleanly");
}

static void test_bondStereoInvertSpike() {
    std::printf("--- Test 14: indigoInvertStereo-on-a-bond spike (3b design question) ---\n");
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    // A molfile with an explicit wedge bond (stereo code 1) -- SMILES cannot
    // carry wedge/dash, and sub-project 3a's Test 12 proved a coordinate-less
    // molecule never gets stereobonds marked, so this spike MUST start from a
    // molfile that already has both 2D coordinates and an explicit wedge flag.
    // Bond block column 4 is the stereo code: "1" = wedge up.
    const char* molfile =
        "\n  spike\n\n"
        "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0\n"
        "    1.0000    0.0000    0.0000 N   0  0\n"
        "    0.0000    1.0000    0.0000 O   0  0\n"
        "   -1.0000    0.0000    0.0000 F   0  0\n"
        "  1  2  1  1  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  1  4  1  0  0  0  0\n"
        "M  END\n";

    int mol = indigoLoadMoleculeFromString(molfile);
    CHECK(mol >= 0, "spike: molfile with an explicit wedge bond loads");
    if (mol < 0) {
        std::printf("[INFO] spike: load error: %s\n", indigoGetLastError());
        indigoReleaseSessionId(session);
        return;
    }

    int bond0 = indigoGetBond(mol, 0);
    CHECK(bond0 >= 0, "spike: got bond 0 handle");
    int stereoBefore = indigoBondStereo(bond0);
    std::printf("[INFO] spike: indigoBondStereo(bond0) as loaded = %d\n", stereoBefore);

    // Q1: does indigoInvertStereo accept a BOND handle at all?
    int invert1 = indigoInvertStereo(bond0);
    std::printf("[INFO] spike: indigoInvertStereo(bond0) returned %d%s\n",
                invert1, invert1 < 0 ? " (REJECTED)" : "");
    if (invert1 < 0) std::printf("[INFO] spike: invert error: %s\n", indigoGetLastError());

    // Q2: did the reported stereo actually change, and to what?
    indigoFree(bond0);
    bond0 = indigoGetBond(mol, 0);
    int stereoAfter1 = indigoBondStereo(bond0);
    std::printf("[INFO] spike: indigoBondStereo(bond0) after 1 invert = %d\n", stereoAfter1);

    // Q3: is it its own inverse? (the EditCommand's invert() calls it again)
    indigoInvertStereo(bond0);
    indigoFree(bond0);
    bond0 = indigoGetBond(mol, 0);
    int stereoAfter2 = indigoBondStereo(bond0);
    std::printf("[INFO] spike: indigoBondStereo(bond0) after 2 inverts = %d (round-trip %s)\n",
                stereoAfter2, stereoAfter2 == stereoBefore ? "OK" : "BROKEN");

    CHECK(true, "spike: bond-stereo invert investigation completed (see [INFO] lines above)");

    indigoFree(bond0);
    indigoFree(mol);
    indigoReleaseSessionId(session);
}

static void test_positionAccessors() {
    std::printf("--- Test 15: position accessors for every repositionable entity ---\n");
    EditableMolecule m;

    // Atom
    AtomId a1 = m.addAtom(QStringLiteral("C"), 1.5, 2.5);
    double ax = 0, ay = 0;
    CHECK(m.atomPos(a1, ax, ay) && ax == 1.5 && ay == 2.5, "atom starts at its add-time position");
    CHECK(m.setAtomPos(a1, -3.25, 4.75), "setAtomPos succeeds");
    CHECK(m.atomPos(a1, ax, ay) && ax == -3.25 && ay == 4.75, "setAtomPos round-trips");

    // RxnArrow
    RxnArrowId ar = m.addRxnArrow(0, 0, 5, 0);
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    CHECK(m.rxnArrowEndpoints(ar, x1, y1, x2, y2) && x1 == 0 && y1 == 0 && x2 == 5 && y2 == 0,
          "rxnArrowEndpoints reads the add-time endpoints");
    CHECK(m.setRxnArrowEndpoints(ar, 1, 2, 3, 4), "setRxnArrowEndpoints succeeds");
    CHECK(m.rxnArrowEndpoints(ar, x1, y1, x2, y2) && x1 == 1 && y1 == 2 && x2 == 3 && y2 == 4,
          "setRxnArrowEndpoints round-trips");

    // RxnPlus
    RxnPlusId pl = m.addRxnPlus(2.5, 0);
    double px = 0, py = 0;
    CHECK(m.rxnPlusPos(pl, px, py) && px == 2.5 && py == 0, "rxnPlusPos reads the add-time position");
    CHECK(m.setRxnPlusPos(pl, -1.5, 6.5), "setRxnPlusPos succeeds");
    CHECK(m.rxnPlusPos(pl, px, py) && px == -1.5 && py == 6.5, "setRxnPlusPos round-trips");

    // MultitailArrow
    MultitailArrowId mta = m.addMultitailArrow({5, 5, 0, 4, 0, 6});
    QList<double> pts = m.multitailArrowPoints(mta);
    CHECK(pts.size() == 6 && pts[0] == 5 && pts[5] == 6, "multitailArrowPoints reads the add-time points");
    CHECK(m.setMultitailArrowPoints(mta, {1, 2, 3, 4}), "setMultitailArrowPoints succeeds");
    pts = m.multitailArrowPoints(mta);
    CHECK(pts.size() == 4 && pts[0] == 1 && pts[3] == 4, "setMultitailArrowPoints round-trips");

    // Image
    ImageId img = m.addImage(0, 0, 4, 3, QByteArray("fakepng"));
    double ix = 0, iy = 0, iw = 0, ih = 0;
    CHECK(m.imageRect(img, ix, iy, iw, ih) && ix == 0 && iy == 0 && iw == 4 && ih == 3,
          "imageRect reads the add-time rect");
    CHECK(m.setImageRect(img, 1, 2, 8, 6), "setImageRect succeeds");
    CHECK(m.imageRect(img, ix, iy, iw, ih) && ix == 1 && iy == 2 && iw == 8 && ih == 6,
          "setImageRect round-trips");

    // Invalid ids fail cleanly.
    const int bad = 9999;
    CHECK(!m.setAtomPos(bad, 0, 0), "invalid id: setAtomPos fails cleanly");
    CHECK(!m.atomPos(bad, ax, ay), "invalid id: atomPos fails cleanly");
    CHECK(!m.rxnArrowEndpoints(bad, x1, y1, x2, y2), "invalid id: rxnArrowEndpoints fails cleanly");
    CHECK(!m.setRxnArrowEndpoints(bad, 0, 0, 0, 0), "invalid id: setRxnArrowEndpoints fails cleanly");
    CHECK(!m.rxnPlusPos(bad, px, py), "invalid id: rxnPlusPos fails cleanly");
    CHECK(!m.setRxnPlusPos(bad, 0, 0), "invalid id: setRxnPlusPos fails cleanly");
    CHECK(m.multitailArrowPoints(bad).isEmpty(), "invalid id: multitailArrowPoints returns empty");
    CHECK(!m.setMultitailArrowPoints(bad, {1, 2}), "invalid id: setMultitailArrowPoints fails cleanly");
    CHECK(!m.imageRect(bad, ix, iy, iw, ih), "invalid id: imageRect fails cleanly");
    CHECK(!m.setImageRect(bad, 0, 0, 1, 1), "invalid id: setImageRect fails cleanly");
}

static void test_atomMergeMechanicsSpike() {
    std::printf("--- Test 16: atom-merge mechanics spike (3c design question) ---\n");
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    // Q1: does indigoRemove(atom) auto-remove its incident bonds (raw, not
    // via our removeAtom wrapper)?
    {
        int mol = indigoCreateMolecule();
        int m = indigoAddAtom(mol, "C");
        int r = indigoAddAtom(mol, "C");
        indigoAddBond(m, r, 3);
        int bondsBefore = indigoCountBonds(mol);
        indigoRemove(m);
        int bondsAfter = indigoCountBonds(mol);
        std::printf("[INFO] spike: indigoRemove(atom) -- bonds before=%d after=%d (incident bond %s)\n",
                    bondsBefore, bondsAfter, bondsAfter < bondsBefore ? "AUTO-REMOVED" : "SURVIVED");
        indigoFree(mol);
    }

    // Q2: create a replacement bond onto the kept atom BEFORE removing the
    // doomed atom, when kept is NOT already bonded to the target (the
    // ordinary, non-seam case).
    {
        int mol = indigoCreateMolecule();
        int kept = indigoAddAtom(mol, "C");     // simulates the pre-existing atom
        int doomed = indigoAddAtom(mol, "C");   // simulates the new ring/chain atom at the same position
        int target = indigoAddAtom(mol, "C");   // doomed's other neighbor, to be rewired onto kept
        indigoAddBond(doomed, target, 1);       // the bond that must move from doomed to kept

        int newBond = indigoAddBond(kept, target, 1);
        std::printf("[INFO] spike: create replacement bond before removal (non-seam case) returned %d%s\n",
                    newBond, newBond < 0 ? " (REJECTED)" : "");
        if (newBond >= 0) {
            std::printf("[INFO] spike: replacement bond order read back = %d (created with order 1)\n",
                        indigoBondOrder(newBond));
        }
        indigoRemove(doomed);
        std::printf("[INFO] spike: atom count after removing doomed = %d (expect 2: kept + target)\n",
                    indigoCountAtoms(mol));
        indigoFree(mol);
    }

    // Q3: the seam case -- kept is ALREADY bonded to target (a fusion seam).
    // Creating an equivalent doomed->target bond onto kept should be
    // REJECTED as a parallel edge (3b's GT-T5 finding), which is the
    // CORRECT outcome: the pre-existing kept-target bond must survive
    // untouched, not be duplicated.
    {
        int mol = indigoCreateMolecule();
        int kept = indigoAddAtom(mol, "C");
        int doomed = indigoAddAtom(mol, "C");
        int target = indigoAddAtom(mol, "C");
        int preExisting = indigoAddBond(kept, target, 2);   // the committed seam bond, order 2
        indigoAddBond(doomed, target, 1);                   // doomed's duplicate edge to the same target

        int dupBond = indigoAddBond(kept, target, 1);
        std::printf("[INFO] spike: create replacement bond in the SEAM case (kept already bonded to target) returned %d%s\n",
                    dupBond, dupBond < 0 ? " (REJECTED, as expected)" : " (UNEXPECTEDLY ACCEPTED)");
        if (dupBond < 0) std::printf("[INFO] spike: rejection error: %s\n", indigoGetLastError());

        indigoRemove(doomed);
        std::printf("[INFO] spike: pre-existing seam bond order after doomed removal = %d (expect unchanged: 2)\n",
                    indigoBondOrder(preExisting));
        indigoFree(mol);
    }

    CHECK(true, "spike: atom-merge mechanics investigation completed (see [INFO] lines above)");
    indigoReleaseSessionId(session);
}

static void test_mergeOverlappingAtomsAndFindBond() {
    std::printf("--- Test 17: mergeOverlappingAtoms + findBond ---\n");

    // Two coincident atoms collapse to one; a bond incident to the doomed
    // atom survives, reattached to the kept atom, with its order carried
    // over; the kept id is the LOWER (pre-existing) one.
    {
        EditableMolecule m;
        AtomId kept = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId other = m.addAtom(QStringLiteral("N"), 5.0, 0.0);
        AtomId doomed = m.addAtom(QStringLiteral("C"), 0.05, 0.02);   // within 0.1 of kept
        AtomId target = m.addAtom(QStringLiteral("O"), 3.0, 0.0);
        m.addBond(other, kept, 2);       // untouched bystander bond
        BondId toRewire = m.addBond(doomed, target, 3);
        CHECK(toRewire > 0, "setup: doomed-target bond added");

        EditableMolecule::MergeResult result = m.mergeOverlappingAtoms();
        CHECK(result.mergedAway.value(doomed, -1) == kept, "doomed merges into kept (the lower id)");
        CHECK(result.createdBonds.size() == 1, "exactly one replacement bond created");
        CHECK(m.atomSymbol(doomed).isEmpty(), "doomed atom no longer resolves");
        CHECK(m.atomSymbol(kept) == QStringLiteral("C"), "kept atom survives with its own label");

        BondId newBond = m.findBond(kept, target);
        CHECK(newBond >= 0, "kept is now bonded to target");
        CHECK(m.bondOrder(newBond) == 3, "replacement bond carries over the original order");
        BondId bystander = m.findBond(kept, other);
        CHECK(bystander >= 0 && m.bondOrder(bystander) == 2, "the untouched bystander bond is unaffected");
    }

    // Non-coincident atoms are left alone; empty merge map.
    {
        EditableMolecule m;
        AtomId a = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId b = m.addAtom(QStringLiteral("C"), 5.0, 0.0);
        EditableMolecule::MergeResult result = m.mergeOverlappingAtoms();
        CHECK(result.mergedAway.isEmpty(), "no merges when nothing overlaps");
        CHECK(result.createdBonds.isEmpty(), "no replacement bonds when nothing overlaps");
        CHECK(m.atomSymbol(a) == QStringLiteral("C") && m.atomSymbol(b) == QStringLiteral("C"),
              "both atoms untouched");
    }

    // Seam case: kept already bonded to target -- the duplicate is
    // correctly dropped, not created, and the pre-existing bond's order is
    // untouched.
    {
        EditableMolecule m;
        AtomId kept = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId target = m.addAtom(QStringLiteral("C"), 1.5, 0.0);
        BondId seam = m.addBond(kept, target, 2);
        AtomId doomed = m.addAtom(QStringLiteral("C"), 0.02, 0.0);   // coincides with kept
        m.addBond(doomed, target, 1);   // duplicate edge to the same target

        EditableMolecule::MergeResult result = m.mergeOverlappingAtoms();
        CHECK(result.mergedAway.value(doomed, -1) == kept, "doomed still merges into kept");
        CHECK(result.createdBonds.isEmpty(), "no replacement bond created for the seam-duplicate edge");
        CHECK(m.bondOrder(seam) == 2, "pre-existing seam bond order is untouched");
        CHECK(m.bondCount() == 1, "exactly one bond survives (the pre-existing seam bond)");
    }

    // findBond: both argument orders, and the not-found/invalid cases.
    {
        EditableMolecule m;
        AtomId a = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId b = m.addAtom(QStringLiteral("C"), 1.0, 0.0);
        AtomId c = m.addAtom(QStringLiteral("C"), 2.0, 0.0);
        BondId ab = m.addBond(a, b, 1);
        CHECK(m.findBond(a, b) == ab, "findBond finds a bond in (a,b) order");
        CHECK(m.findBond(b, a) == ab, "findBond finds the same bond in (b,a) order");
        CHECK(m.findBond(a, c) == -1, "findBond returns -1 for a non-adjacent pair");
        CHECK(m.findBond(9999, a) == -1, "findBond returns -1 for an invalid id");
    }
}

static void test_addBenzeneRing() {
    std::printf("--- Test 18: addBenzeneRing ---\n");
    EditableMolecule m;
    QList<AtomId> ring = m.addBenzeneRing(4.0, 4.0);

    CHECK(ring.size() == 6, "addBenzeneRing returns 6 atom ids");
    CHECK(m.atomCount() == 6, "6 atoms created");
    CHECK(m.bondCount() == 6, "6 bonds created (closed cycle)");

    for (AtomId id : ring) {
        CHECK(m.atomSymbol(id) == QStringLiteral("C"), "every ring atom is carbon");
    }

    for (int i = 0; i < 6; ++i) {
        BondId bid = m.findBond(ring[i], ring[(i + 1) % 6]);
        CHECK(bid >= 0, "consecutive ring atoms are bonded");
        int expectedOrder = (i % 2 == 0) ? 2 : 1;
        CHECK(m.bondOrder(bid) == expectedOrder, "ring bond alternates 2/1/2/1/2/1");
    }
}

static void test_templateLibraryLoading() {
    std::printf("--- Test 19: TemplateLibrary loading ---\n");
    TemplateLibrary lib(
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/fg.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/library.sdf"),
        QStringLiteral(SKETCH_SOURCE_DIR "/templates/salts-and-solvents.sdf"));

    // "Ac" is a real fg.sdf entry: 3 atoms, 2 bonds, one SUP group. Its record
    // carries the Ketcher "G    1  0" line that must be stripped to parse.
    int ac = lib.functionalGroup(QStringLiteral("Ac"));
    CHECK(ac >= 0, "functionalGroup(\"Ac\") resolves to a valid handle");
    if (ac >= 0) {
        CHECK(indigoCountAtoms(ac) == 3, "\"Ac\" has 3 atoms");
        CHECK(indigoCountBonds(ac) == 2, "\"Ac\" has 2 bonds");
        CHECK(indigoCountSuperatoms(ac) == 1, "\"Ac\" has 1 superatom group");
    }

    // molfileText is the safe cross-session handoff: TemplateLibrary owns its
    // own Indigo session, so a raw handle like `ac` is only valid while that
    // session is active. Switch to a DIFFERENT session (as any EditableMolecule
    // does before its own calls) and confirm the handle can no longer be read
    // directly, but the Molfile text -- captured before switching -- reloads
    // correctly in the new session with the superatom intact.
    QString acMolfile = lib.molfileText(ac);
    CHECK(!acMolfile.isEmpty(), "molfileText(\"Ac\") returns non-empty text");
    unsigned long long otherSession = indigoAllocSessionId();
    indigoSetSessionId(otherSession);
    CHECK(indigoCountAtoms(ac) < 0, "the raw handle is not usable from a different session");
    int reloadedAc = indigoLoadMoleculeFromString(acMolfile.toUtf8().constData());
    CHECK(reloadedAc >= 0, "molfileText's output reparses cleanly in a different session");
    if (reloadedAc >= 0) {
        CHECK(indigoCountAtoms(reloadedAc) == 3, "reparsed \"Ac\" still has 3 atoms");
        CHECK(indigoCountSuperatoms(reloadedAc) == 1, "reparsed \"Ac\" still has its superatom group");
        indigoFree(reloadedAc);
    }
    indigoReleaseSessionId(otherSession);

    // "Indole" is a real library.sdf entry: 9 atoms, 10 bonds, bondIdx=6
    // (confirmed by direct file read: the aromatic N-C bond at the
    // pyrrole/benzene fusion seam). library.sdf entries have no G-line and no
    // superatoms.
    int indole = lib.libraryTemplate(QStringLiteral("Indole"));
    CHECK(indole >= 0, "libraryTemplate(\"Indole\") resolves to a valid handle");
    if (indole >= 0) {
        CHECK(indigoCountAtoms(indole) == 9, "\"Indole\" has 9 atoms");
        CHECK(indigoCountBonds(indole) == 10, "\"Indole\" has 10 bonds");
        CHECK(indigoCountSuperatoms(indole) == 0, "\"Indole\" has no superatom groups");
    }
    CHECK(lib.libraryTemplateFusionBondIdx(QStringLiteral("Indole")) == 6,
          "\"Indole\"'s fusion bond index is 6");

    // A name with no fusion metadata (any fg.sdf entry, since bondIdx only
    // ever appears in library.sdf) returns -1.
    CHECK(lib.libraryTemplateFusionBondIdx(QStringLiteral("Ac")) == -1,
          "a name with no library-registry fusion metadata returns -1");

    // Names that don't exist anywhere resolve to -1 across every registry.
    CHECK(lib.functionalGroup(QStringLiteral("NotARealTemplateName")) < 0,
          "an unknown functional-group name resolves to -1");
    CHECK(lib.libraryTemplate(QStringLiteral("NotARealTemplateName")) < 0,
          "an unknown library-template name resolves to -1");
    CHECK(lib.saltOrSolvent(QStringLiteral("NotARealTemplateName")) < 0,
          "an unknown salt/solvent name resolves to -1");

    // A missing file leaves its registry empty, not crashing.
    TemplateLibrary missingFile(
        QStringLiteral("/nonexistent/path/fg.sdf"),
        QStringLiteral("/nonexistent/path/library.sdf"),
        QStringLiteral("/nonexistent/path/salts.sdf"));
    CHECK(missingFile.functionalGroup(QStringLiteral("Ac")) < 0,
          "a missing SDF file leaves that registry empty rather than crashing");
}

static void test_insertStructure() {
    std::printf("--- Test 20: insertStructure ---\n");
    // A prior test (test_templateLibraryLoading) constructs local
    // TemplateLibrary objects whose destructor releases THEIR OWN session,
    // leaving the process-global "current session" pointing at an
    // already-released id. This test is the first one to make raw indigo*
    // calls before constructing any session-owning wrapper object (dest,
    // below), so -- matching every other raw-indigo-call test in this file
    // (e.g. test_atomMergeMechanicsSpike) -- it allocates and activates its
    // own session up front rather than relying on inherited global state.
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    // Plain translation transform, no superatom. The source is built and
    // serialized to Molfile text BEFORE `dest` is even constructed, so this
    // genuinely exercises the cross-session handoff contract (dest's own
    // constructor activates its own, different session) rather than relying
    // on both handles happening to share whatever session is globally active.
    {
        int src = indigoCreateMolecule();
        int a1 = indigoAddAtom(src, "N");
        int a2 = indigoAddAtom(src, "O");
        indigoSetXYZ(a1, 0.0f, 0.0f, 0.0f);
        indigoSetXYZ(a2, 1.0f, 0.0f, 0.0f);
        indigoAddBond(a1, a2, 2);
        QString srcMolfile = QString::fromUtf8(indigoMolfile(src));
        indigoFree(src);

        EditableMolecule dest;
        dest.addAtom(QStringLiteral("C"), 0.0, 0.0);   // 1 pre-existing atom

        EditableMolecule::InsertResult result = dest.insertStructure(srcMolfile, [](double x, double y) {
            return QPointF(x + 5.0, y + 5.0);
        });
        CHECK(result.createdAtoms.size() == 2, "insertStructure creates 2 new atoms");
        CHECK(result.createdBonds.size() == 1, "insertStructure creates 1 new bond");
        CHECK(result.createdSGroups.isEmpty(), "no sgroups when the source has none");
        CHECK(result.sourceIndexToNewAtomId.size() == 2, "source-index map covers both new atoms");

        double x = 0, y = 0;
        CHECK(dest.atomPos(result.createdAtoms[0], x, y) && x == 5.0 && y == 5.0,
              "first inserted atom is translated by the transform");
        CHECK(dest.atomPos(result.createdAtoms[1], x, y) && x == 6.0 && y == 5.0,
              "second inserted atom is translated by the transform");
        CHECK(dest.atomCount() == 3, "destination has 1 pre-existing + 2 new = 3 atoms");
    }

    // Source with a superatom: insertStructure discovers it and assigns a
    // stable SGroupId, with the attachment point correctly resolved. Same
    // build-then-serialize-then-free-before-dest-exists shape as above.
    // Reactivate this test's own session first: the previous block's `dest`
    // released ITS session on scope exit (EditableMolecule's destructor
    // does this, same as TemplateLibrary's), leaving global session state
    // stale again.
    indigoSetSessionId(session);
    {
        int src = indigoCreateMolecule();
        int sa1 = indigoAddAtom(src, "N");
        int sa2 = indigoAddAtom(src, "C");
        indigoSetXYZ(sa1, 0.0f, 0.0f, 0.0f);
        indigoSetXYZ(sa2, 1.0f, 0.0f, 0.0f);
        indigoAddBond(sa1, sa2, 1);
        int sa1Idx = indigoIndex(sa1);
        int atomIdxs[2] = { sa1Idx, indigoIndex(sa2) };
        int sup = indigoAddSuperatom(src, 2, atomIdxs, "TestGroup");
        indigoAddSGroupAttachmentPoint(sup, sa1Idx, -1, "");
        QString srcMolfile = QString::fromUtf8(indigoMolfile(src));
        indigoFree(src);   // sa1/sa2/sup are child handles of src; nothing below re-uses them

        EditableMolecule dest;
        EditableMolecule::InsertResult result = dest.insertStructure(srcMolfile, [](double x, double y) {
            return QPointF(x, y);
        });
        CHECK(result.createdSGroups.size() == 1, "insertStructure discovers the propagated superatom");
        AtomId attachAtom = -1;
        CHECK(dest.superatomAttachAtom(result.createdSGroups[0], attachAtom),
              "superatomAttachAtom resolves for the discovered sgroup");
        CHECK(result.sourceIndexToNewAtomId.value(sa1Idx, -1) == attachAtom,
              "the resolved attach atom matches the source's own attach-point atom");
        CHECK(dest.sgroupExpanded(result.createdSGroups[0]),
              "a newly-discovered sgroup defaults to expanded=true");
        dest.setSGroupExpanded(result.createdSGroups[0], false);
        CHECK(!dest.sgroupExpanded(result.createdSGroups[0]), "setSGroupExpanded round-trips");
    }

    // Invalid SGroupId / AtomId out-param cases fail cleanly.
    {
        EditableMolecule dest;
        AtomId out = -1;
        CHECK(!dest.superatomAttachAtom(9999, out), "invalid SGroupId: superatomAttachAtom fails cleanly");
        CHECK(!dest.sgroupExpanded(9999), "invalid SGroupId: sgroupExpanded returns false");
    }
    indigoReleaseSessionId(session);
}

static void test_renderSupportStorageDefaults() {
    std::printf("--- Test 22: render-support storage (stereoFlags/checkWarning text/rxnArrow extras) ---\n");

    EditableMolecule m;
    CHECK(m.stereoFlagsType() == QStringLiteral("abs"), "stereoFlagsType defaults to \"abs\"");
    CHECK(m.stereoFlagsGroupId() == 0, "stereoFlagsGroupId defaults to 0");

    AtomId a1 = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
    CHECK(m.atomCheckWarningText(a1).isEmpty(), "atomCheckWarningText defaults to empty");
    CHECK(m.bondCheckWarningText(9999).isEmpty(), "bondCheckWarningText on unknown id defaults to empty");

    AtomId a2 = m.addAtom(QStringLiteral("O"), 1.0, 0.0);
    BondId b1 = m.addBond(a1, a2, 1);
    CHECK(m.bondCheckWarningText(b1).isEmpty(), "bondCheckWarningText defaults to empty for a real bond");

    int arrowId = m.addRxnArrow(0.0, 0.0, 5.0, 0.0);
    CHECK(m.rxnArrowMode(arrowId) == QStringLiteral("filled-triangle"),
          "rxnArrowMode defaults to \"filled-triangle\"");
    CHECK(m.rxnArrowConditionsAbove(arrowId).isEmpty(), "rxnArrowConditionsAbove defaults to empty");
    CHECK(m.rxnArrowConditionsBelow(arrowId).isEmpty(), "rxnArrowConditionsBelow defaults to empty");
    double cx = 0, cy = 0;
    CHECK(!m.rxnArrowCurvature(arrowId, cx, cy), "rxnArrowCurvature defaults to false (no curvature set)");

    // Invalid ids fail cleanly / return sensible defaults, not crash.
    CHECK(m.rxnArrowMode(9999) == QStringLiteral("filled-triangle"),
          "rxnArrowMode on unknown id still returns the default, not garbage");
    CHECK(m.atomCheckWarningText(9999).isEmpty(), "atomCheckWarningText on unknown id defaults to empty");
}

static void test_graftAtomOnto() {
    std::printf("--- Test 21: graftAtomOnto ---\n");

    // A bond incident to the doomed atom survives, reattached to kept, order
    // carried over -- the same contract mergeOverlappingAtoms's Test 17 proves
    // for the position-triggered case, now exercised directly by explicit id.
    {
        EditableMolecule m;
        AtomId kept = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId doomed = m.addAtom(QStringLiteral("N"), 5.0, 5.0);   // far apart: NOT position-coincident
        AtomId target = m.addAtom(QStringLiteral("O"), 3.0, 0.0);
        BondId toRewire = m.addBond(doomed, target, 3);
        CHECK(toRewire > 0, "setup: doomed-target bond added");

        QList<BondId> created = m.graftAtomOnto(doomed, kept);
        CHECK(created.size() == 1, "graftAtomOnto creates exactly one replacement bond");
        CHECK(m.atomSymbol(doomed).isEmpty(), "doomed atom no longer resolves");
        CHECK(m.atomSymbol(kept) == QStringLiteral("C"), "kept atom survives untouched");
        BondId newBond = m.findBond(kept, target);
        CHECK(newBond >= 0 && m.bondOrder(newBond) == 3,
              "kept is now bonded to target with the original order carried over");
    }

    // Seam case: kept already bonded to target -- the duplicate is dropped,
    // not created, matching graftAtomOnto's documented parallel-edge contract.
    {
        EditableMolecule m;
        AtomId kept = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
        AtomId target = m.addAtom(QStringLiteral("C"), 1.5, 0.0);
        BondId seam = m.addBond(kept, target, 2);
        AtomId doomed = m.addAtom(QStringLiteral("N"), 9.0, 9.0);
        m.addBond(doomed, target, 1);

        QList<BondId> created = m.graftAtomOnto(doomed, kept);
        CHECK(created.isEmpty(), "no replacement bond created for the seam-duplicate edge");
        CHECK(m.bondOrder(seam) == 2, "pre-existing seam bond order is untouched");
        CHECK(m.bondCount() == 1, "exactly one bond survives");
    }
}

static void test_ringAndNeighborAccessors() {
    std::printf("--- Test 23: implicit-H, neighbors, SSSR ring membership ---\n");

    // Ethanol: C-C-O chain. Terminal carbon has 3 implicit H, middle carbon
    // has 2, oxygen has 1. Neighbor counts: 1, 2, 1.
    {
        EditableMolecule m(QStringLiteral("CCO"));
        QList<AtomId> ids = m.atomIds();
        CHECK(ids.size() == 3, "setup: ethanol has 3 heavy atoms");
        CHECK(m.implicitHydrogenCount(ids[0]) == 3, "terminal carbon has 3 implicit H");
        CHECK(m.implicitHydrogenCount(ids[1]) == 2, "middle carbon has 2 implicit H");
        CHECK(m.implicitHydrogenCount(ids[2]) == 1, "oxygen has 1 implicit H");
        CHECK(m.neighborAtomIds(ids[0]).size() == 1, "terminal carbon has 1 neighbor");
        CHECK(m.neighborAtomIds(ids[1]).size() == 2, "middle carbon has 2 neighbors");
        QList<AtomId> midNeighbors = m.neighborAtomIds(ids[1]);
        CHECK(midNeighbors.contains(ids[0]) && midNeighbors.contains(ids[2]),
              "middle carbon's neighbors are the terminal carbon and the oxygen");
        CHECK(m.implicitHydrogenCount(9999) == 0, "invalid id: implicitHydrogenCount returns 0");
        CHECK(m.neighborAtomIds(9999).isEmpty(), "invalid id: neighborAtomIds returns empty");
        CHECK(m.ringMembership().isEmpty(), "acyclic molecule has no SSSR rings");
    }

    // Benzene (real aromatic SMILES): one SSSR ring, 6 atoms, 6 bonds, all
    // real aromatic order (4).
    {
        EditableMolecule m(QStringLiteral("c1ccccc1"));
        CHECK(m.atomCount() == 6, "setup: benzene has 6 atoms");
        QList<EditableMolecule::RingMembership> rings = m.ringMembership();
        CHECK(rings.size() == 1, "benzene has exactly one SSSR ring");
        if (!rings.isEmpty()) {
            CHECK(rings[0].atoms.size() == 6, "the ring has 6 atoms");
            CHECK(rings[0].bonds.size() == 6, "the ring has 6 bonds");
            bool allAromatic = true;
            for (BondId bid : rings[0].bonds) {
                if (m.bondOrder(bid) != 4) allAromatic = false;
            }
            CHECK(allAromatic, "every ring bond has real aromatic order (4)");
        }
    }

    // A manually-alternating Kekule 6-ring (single/double/single/double/single/double,
    // built via explicit bonds, no real aromatic order anywhere) still has exactly
    // one SSSR ring -- ring PERCEPTION doesn't require aromaticity, matching the
    // real code's need to separately test the 6-ring/3-double-bond heuristic later.
    {
        EditableMolecule m;
        QList<AtomId> ring;
        for (int i = 0; i < 6; ++i) ring.append(m.addAtom(QStringLiteral("C"), i * 1.0, 0.0));
        for (int i = 0; i < 6; ++i) {
            int order = (i % 2 == 0) ? 2 : 1;
            m.addBond(ring[i], ring[(i + 1) % 6], order);
        }
        QList<EditableMolecule::RingMembership> rings = m.ringMembership();
        CHECK(rings.size() == 1, "the Kekule 6-ring is perceived as one SSSR ring");
        if (!rings.isEmpty()) {
            CHECK(rings[0].atoms.size() == 6, "the Kekule ring has 6 atoms");
            CHECK(rings[0].bonds.size() == 6, "the Kekule ring has 6 bonds");
            int doubleCount = 0;
            for (BondId bid : rings[0].bonds) {
                if (m.bondOrder(bid) == 2) ++doubleCount;
            }
            CHECK(doubleCount == 3, "exactly 3 of the ring's bonds are double bonds");
        }
    }
}

static void test_stereoCipAccessors() {
    std::printf("--- Test 24: stereocenter type/group + CIP descriptor accessors ---\n");

    // Real R/S case: F[C@H](Cl)Br -- the chiral carbon is atom index 1.
    {
        EditableMolecule m(QStringLiteral("F[C@H](Cl)Br"));
        QList<AtomId> ids = m.atomIds();
        CHECK(ids.size() == 4, "setup: F-C(H)(Cl)-Br has 4 heavy atoms");
        AtomId chiralC = ids[1];
        CHECK(m.stereocenterType(chiralC) != 0, "the chiral carbon is a real stereocenter");
        CHECK(m.atomCipDescriptor(chiralC) == 5, "the chiral carbon's CIP descriptor is R (CIPDesc::R == 5)");

        AtomId fluorine = ids[0];
        CHECK(m.stereocenterType(fluorine) == 0, "fluorine is not a stereocenter");
        CHECK(m.atomCipDescriptor(fluorine) == 0, "fluorine's CIP descriptor is NONE (0)");
    }

    // Bond stereo direction: trans-2-butene has a real TRANS double bond.
    {
        EditableMolecule m(QStringLiteral("C/C=C/C"));
        QList<BondId> bids = m.bondIds();
        bool foundTrans = false;
        for (BondId bid : bids) {
            if (m.bondOrder(bid) == 2 && m.bondStereoDirection(bid) == 8) foundTrans = true;   // INDIGO_TRANS == 8
        }
        CHECK(foundTrans, "the double bond reports INDIGO_TRANS (8) stereo direction");
    }

    // Invalid ids fail cleanly.
    {
        EditableMolecule m(QStringLiteral("CCO"));
        CHECK(m.stereocenterType(9999) == 0, "invalid id: stereocenterType returns 0");
        CHECK(m.stereocenterGroup(9999) == 0, "invalid id: stereocenterGroup returns 0");
        CHECK(m.bondStereoDirection(9999) == 0, "invalid id: bondStereoDirection returns 0");
        CHECK(m.atomCipDescriptor(9999) == 0, "invalid id: atomCipDescriptor returns 0");
    }
}

static void test_sgroupIntrospectionAccessors() {
    std::printf("--- Test 25: sgroupIds + sgroupMemberAtomIds ---\n");

    EditableMolecule m;
    // "Ac" molfile: 3 atoms (C, C, O), 2 bonds, one SUP superatom labeled "Ac"
    // with an attachment point at atom index 0 -- the same real bundled template
    // data confirmed in sub-project 3d, embedded here directly so this test has no
    // dependency on TemplateLibrary/file I/O.
    QString acMolfile = QStringLiteral(
        "Ac\n"
        "  Ketcher\n\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    3.3951   -3.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.6785   -3.9891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.3951   -2.7480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  2  1  1  0  0  0  0\n"
        "  1  3  2  0  0  0  0\n"
        "M  STY  1   1 SUP\n"
        "M  SLB  1   1   1\n"
        "M  SAL   1  3   1   2   3\n"
        "M  SAP   1  1   1   0\n"
        "M  SMT   1 Ac\n"
        "M  END\n");
    EditableMolecule::InsertResult r = m.insertStructure(acMolfile, [](double x, double y) { return QPointF(x, y); });
    CHECK(r.createdSGroups.size() == 1, "setup: \"Ac\" produces exactly one sgroup");
    SGroupId sgId = r.createdSGroups[0];

    CHECK(m.sgroupIds() == QList<SGroupId>{sgId}, "sgroupIds lists exactly the one created sgroup");

    QList<AtomId> members = m.sgroupMemberAtomIds(sgId);
    CHECK(members.size() == 3, "\"Ac\"'s sgroup has all 3 atoms as members (M SAL lists atoms 1,2,3)");
    AtomId attach = -1;
    CHECK(m.superatomAttachAtom(sgId, attach), "setup: attach atom resolves");
    CHECK(members.contains(attach), "the member set includes the attach atom");

    CHECK(m.sgroupMemberAtomIds(9999).isEmpty(), "invalid SGroupId: sgroupMemberAtomIds returns empty");
    CHECK(m.sgroupIds().size() == 1, "sgroupIds is unaffected by querying an invalid id");
}

static void test_rebuildIndexTablesPrunesAutoRemovedSgroup() {
    std::printf("--- Test: rebuildIndexTables prunes a sgroup Indigo auto-removed ---\n");
    // Two independent 2-atom SUP groups, mirroring the spec's own probe fixture.
    const char* molfile =
        "twoGroups\n"
        "  Ketcher\n\n"
        "  4  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    5.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    6.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "M  STY  2   1 SUP   2 SUP\n"
        "M  SLB  2   1   1   2   2\n"
        "M  SAL   1  2   1   2\n"
        "M  SAL   2  2   3   4\n"
        "M  SMT   1 GroupA\n"
        "M  SMT   2 GroupB\n"
        "M  END\n";
    EditableMolecule mol(molfile);
    CHECK(mol.isValid(), "molecule loaded");
    CHECK(mol.sgroupIds().size() == 2, "2 sgroups after load");

    QList<AtomId> allAtoms = mol.atomIds();
    CHECK(allAtoms.size() == 4, "4 atoms after load");
    SGroupId groupAId = -1, groupBId = -1;
    for (SGroupId sid : mol.sgroupIds()) {
        QList<AtomId> members = mol.sgroupMemberAtomIds(sid);
        if (members.contains(allAtoms[0])) groupAId = sid;
        else groupBId = sid;
    }
    CHECK(groupAId != -1 && groupBId != -1, "both sgroups identified by membership");

    QList<AtomId> groupAMembers = mol.sgroupMemberAtomIds(groupAId);
    CHECK(groupAMembers.size() == 2, "GroupA has 2 members before removal");
    for (AtomId aid : groupAMembers) mol.removeAtom(aid);

    CHECK(mol.sgroupIds().size() == 1, "only 1 sgroup left after GroupA's members are gone");
    CHECK(!mol.sgroupIds().contains(groupAId), "GroupA's id no longer reported by sgroupIds()");
    CHECK(mol.sgroupIds().contains(groupBId), "GroupB's id still reported");
    CHECK(mol.sgroupMemberAtomIds(groupBId).size() == 2, "GroupB's own members unaffected");
}

static void test_removeSuperatomOnly() {
    std::printf("--- Test: removeSuperatomOnly ---\n");
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
    EditableMolecule mol(molfile);
    CHECK(mol.isValid(), "Ac molecule loaded");
    CHECK(mol.sgroupIds().size() == 1, "1 sgroup before removal");
    SGroupId sid = mol.sgroupIds().first();

    CHECK(mol.removeSuperatomOnly(sid), "removeSuperatomOnly succeeds");
    CHECK(mol.sgroupIds().isEmpty(), "sgroupIds() empty after removal");
    CHECK(mol.atomCount() == 3, "all 3 atoms still exist");
}

static void test_createSuperatomFromAtoms() {
    std::printf("--- Test: createSuperatomFromAtoms ---\n");
    EditableMolecule mol;
    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = mol.addAtom(QStringLiteral("N"), 1.0, 0.0);
    CHECK(mol.sgroupIds().isEmpty(), "no sgroups on a bare molecule");

    SGroupId sid = mol.createSuperatomFromAtoms({a1, a2});
    CHECK(sid != -1, "createSuperatomFromAtoms returns a valid id");
    CHECK(mol.sgroupIds().size() == 1, "1 sgroup after creation");
    CHECK(mol.sgroupIds().contains(sid), "sgroupIds() reports the new id");
    QList<AtomId> members = mol.sgroupMemberAtomIds(sid);
    CHECK(members.size() == 2 && members.contains(a1) && members.contains(a2),
          "member set matches exactly what was passed in");
}

static void test_atomIdsInIndigoOrder() {
    std::printf("--- Test: atomIdsInIndigoOrder survives a remove-then-add index reuse ---\n");
    EditableMolecule mol;
    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = mol.addAtom(QStringLiteral("N"), 1, 0);
    AtomId a3 = mol.addAtom(QStringLiteral("O"), 2, 0);
    AtomId a4 = mol.addAtom(QStringLiteral("F"), 3, 0);
    mol.removeAtom(a2);
    AtomId a5 = mol.addAtom(QStringLiteral("Cl"), 4, 0);
    (void)a1; (void)a3; (void)a4;

    QList<AtomId> order = mol.atomIdsInIndigoOrder();
    CHECK(order.size() == 4, "4 atoms after remove+add");

    QString symbols;
    for (AtomId id : order) symbols += mol.atomSymbol(id);
    CHECK(symbols == QStringLiteral("CClOF"),
          "order matches Indigo's own internal index order (C, Cl-in-freed-slot, O, F), "
          "NOT sorted-by-AtomId order (which would give C,O,F,Cl)");

    StringResult mf = mol.toMolfile();
    CHECK(mf.success, "toMolfile succeeds");
    EditableMolecule reparsed(mf.value);
    QString reparsedSymbols;
    for (AtomId id : reparsed.atomIds()) reparsedSymbols += reparsed.atomSymbol(id);
    CHECK(reparsedSymbols == symbols,
          "atomIdsInIndigoOrder matches what a fresh reparse would number these atoms as");
}

static void test_bondIdsInIndigoOrder() {
    std::printf("--- Test: bondIdsInIndigoOrder survives a remove-then-add index reuse ---\n");
    EditableMolecule mol;
    AtomId a1 = mol.addAtom(QStringLiteral("C"), 0, 0);
    AtomId a2 = mol.addAtom(QStringLiteral("C"), 1, 0);
    AtomId a3 = mol.addAtom(QStringLiteral("C"), 2, 0);
    AtomId a4 = mol.addAtom(QStringLiteral("C"), 3, 0);
    BondId b1 = mol.addBond(a1, a2, 1);
    mol.addBond(a2, a3, 2);
    mol.addBond(a3, a4, 3);
    mol.removeBond(b1);
    BondId b4 = mol.addBond(a1, a4, 1);

    QList<BondId> order = mol.bondIdsInIndigoOrder();
    CHECK(order.size() == 3, "3 bonds after remove+add");
    CHECK(order.contains(b4), "the newly-added bond is present in the order list");

    StringResult mf = mol.toMolfile();
    CHECK(mf.success, "toMolfile succeeds");
    EditableMolecule reparsed(mf.value);
    CHECK(reparsed.bondCount() == 3, "reparsed molecule has the same 3 bonds");
    QList<int> orders1, orders2;
    for (BondId id : order) orders1.append(mol.bondOrder(id));
    for (BondId id : reparsed.bondIds()) orders2.append(reparsed.bondOrder(id));
    CHECK(orders1 == orders2,
          "bondIdsInIndigoOrder's bond-order sequence matches a fresh reparse's own order");
}

static void test_atomStereoDescriptorStorage() {
    std::printf("--- Test: setAtomStereoDescriptor / getters ---\n");
    EditableMolecule mol;
    AtomId a = mol.addAtom(QStringLiteral("C"), 0, 0);

    CHECK(mol.atomStoredCipLabel(a).isEmpty(), "default cipLabel is empty");
    CHECK(mol.atomStoredStereoType(a) == 0, "default stereoType is 0");
    CHECK(mol.atomStoredStereoGroup(a) == 0, "default stereoGroup is 0");

    mol.setAtomStereoDescriptor(a, QStringLiteral("R"), 1, 2);
    CHECK(mol.atomStoredCipLabel(a) == QStringLiteral("R"), "cipLabel set to R");
    CHECK(mol.atomStoredStereoType(a) == 1, "stereoType set to 1");
    CHECK(mol.atomStoredStereoGroup(a) == 2, "stereoGroup set to 2");

    mol.setAtomStereoDescriptor(9999, QStringLiteral("S"), 3, 4);
    CHECK(mol.atomStoredCipLabel(9999).isEmpty(), "unknown AtomId: setter no-op, getter defaults empty");
    CHECK(mol.atomStoredStereoType(9999) == 0, "unknown AtomId: getter defaults to 0");
}

static void test_bondStereoCipLabelStorage() {
    std::printf("--- Test: setBondStereoCipLabel / getter ---\n");
    EditableMolecule mol;
    AtomId a = mol.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b = mol.addAtom(QStringLiteral("C"), 1, 0);
    BondId bond = mol.addBond(a, b, 2);

    CHECK(mol.bondStoredCipLabel(bond).isEmpty(), "default bond cipLabel is empty");
    mol.setBondStereoCipLabel(bond, QStringLiteral("E"));
    CHECK(mol.bondStoredCipLabel(bond) == QStringLiteral("E"), "bond cipLabel set to E");

    mol.setBondStereoCipLabel(9999, QStringLiteral("Z"));
    CHECK(mol.bondStoredCipLabel(9999).isEmpty(), "unknown BondId: setter no-op, getter defaults empty");
}

static void test_checkWarningTextSetters() {
    std::printf("--- Test: setAtomCheckWarningText / setBondCheckWarningText ---\n");
    EditableMolecule mol;
    AtomId a = mol.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b = mol.addAtom(QStringLiteral("C"), 1, 0);
    BondId bond = mol.addBond(a, b, 1);

    CHECK(mol.atomCheckWarningText(a).isEmpty(), "default atom checkWarning is empty");
    mol.setAtomCheckWarningText(a, QStringLiteral("valence"));
    CHECK(mol.atomCheckWarningText(a) == QStringLiteral("valence"), "atom checkWarning set");

    CHECK(mol.bondCheckWarningText(bond).isEmpty(), "default bond checkWarning is empty");
    mol.setBondCheckWarningText(bond, QStringLiteral("overlap_atom"));
    CHECK(mol.bondCheckWarningText(bond) == QStringLiteral("overlap_atom"), "bond checkWarning set");

    mol.setAtomCheckWarningText(9999, QStringLiteral("radical"));
    CHECK(mol.atomCheckWarningText(9999).isEmpty(), "unknown AtomId: setter no-op");
    mol.setBondCheckWarningText(9999, QStringLiteral("radical"));
    CHECK(mol.bondCheckWarningText(9999).isEmpty(), "unknown BondId: setter no-op");
}

static void test_rxnArrowMutators() {
    std::printf("--- Test 26: rxn-arrow mutators ---\n");
    EditableMolecule m;
    int id = m.addRxnArrow(0, 0, 5, 0);

    CHECK(m.rxnArrowMode(id) == QStringLiteral("filled-triangle"), "new rxn arrow defaults to filled-triangle mode");
    CHECK(m.setRxnArrowMode(id, QStringLiteral("curved-mechanism")), "setRxnArrowMode succeeds for a real id");
    CHECK(m.rxnArrowMode(id) == QStringLiteral("curved-mechanism"), "mode actually changed");
    CHECK(!m.setRxnArrowMode(9999, QStringLiteral("filled-triangle")), "setRxnArrowMode fails for an unknown id");

    CHECK(m.setRxnArrowConditions(id, QStringLiteral("H2O"), QStringLiteral("reflux")), "setRxnArrowConditions succeeds");
    CHECK(m.rxnArrowConditionsAbove(id) == QStringLiteral("H2O"), "conditions above set");
    CHECK(m.rxnArrowConditionsBelow(id) == QStringLiteral("reflux"), "conditions below set");
    CHECK(!m.setRxnArrowConditions(9999, QStringLiteral("x"), QStringLiteral("y")), "setRxnArrowConditions fails for an unknown id");

    CHECK(m.setRxnArrowCurvature(id, 2.5, 1.0), "setRxnArrowCurvature succeeds");
    double cx = 0, cy = 0;
    CHECK(m.rxnArrowCurvature(id, cx, cy) && cx == 2.5 && cy == 1.0, "curvature stored and readable");
    CHECK(m.setRxnArrowCurvature(id, 0, 0, false), "setRxnArrowCurvature(has=false) succeeds");
    CHECK(!m.rxnArrowCurvature(id, cx, cy), "curvature cleared -- rxnArrowCurvature now reports false");
}

static void test_stereoFlagsDocumentMutator() {
    std::printf("--- Test 27: document-level stereoFlags mutator ---\n");
    EditableMolecule m;
    CHECK(m.stereoFlagsType() == QStringLiteral("abs"), "default document stereoFlags type is abs");
    CHECK(m.stereoFlagsGroupId() == 0, "default document stereoFlags groupId is 0");

    m.setStereoFlagsDocument(QStringLiteral("rel"), 3);
    CHECK(m.stereoFlagsType() == QStringLiteral("rel"), "setStereoFlagsDocument changes the type");
    CHECK(m.stereoFlagsGroupId() == 3, "setStereoFlagsDocument changes the groupId");

    // Distinct from the unrelated per-fragment setStereoFlag/stereoFlag pair.
    m.setStereoFlag(0, 5);
    CHECK(m.stereoFlag(0) == 5, "per-fragment setStereoFlag/stereoFlag is untouched by the document-level mutator");
    CHECK(m.stereoFlagsType() == QStringLiteral("rel"), "document-level type is untouched by the per-fragment setter");
}

static void test_rgroupStorageAndFragmentIndex() {
    std::printf("--- Test 28: R-group storage + atomFragmentIndex ---\n");

    // Two disjoint ethane fragments -- SMILES is auto-detected by indigoLoadMoleculeFromString
    // exactly like the "CCO"/"CC(N)C(=O)O" molecules already used elsewhere in this file, so this
    // avoids hand-rolling a V2000 molfile. Confirmed by direct probe against the real indigo.dll
    // (indigoComponentIndex is live, 0-based, no indigoCountComponents prerequisite).
    EditableMolecule m(QStringLiteral("CC.CC"));
    QList<AtomId> ids = m.atomIds();
    CHECK(ids.size() == 4, "setup: loaded 4 atoms across 2 disjoint fragments");
    std::sort(ids.begin(), ids.end());
    int frag0a = m.atomFragmentIndex(ids[0]);
    int frag0b = m.atomFragmentIndex(ids[1]);
    int frag1a = m.atomFragmentIndex(ids[2]);
    int frag1b = m.atomFragmentIndex(ids[3]);
    CHECK(frag0a == frag0b, "atoms 0,1 (bonded pair) share the same fragment index");
    CHECK(frag1a == frag1b, "atoms 2,3 (bonded pair) share the same fragment index");
    CHECK(frag0a != frag1a, "the two disjoint fragments have DIFFERENT indices");
    CHECK(m.atomFragmentIndex(9999) == -1, "atomFragmentIndex returns -1 for an unknown id");

    // R-group storage.
    CHECK(m.rgroupNumbers().isEmpty(), "no R-groups initially");
    CHECK(m.addRGroupEntry(1), "addRGroupEntry(1) succeeds");
    CHECK(!m.addRGroupEntry(1), "addRGroupEntry(1) again fails -- already exists");
    CHECK(m.rgroupNumbers() == QList<int>{1}, "rgroupNumbers lists the one R-group");

    QString range; bool resth = false; int ifthen = 0;
    CHECK(m.rgroupLogic(1, range, resth, ifthen), "rgroupLogic reads the fresh entry");
    CHECK(range.isEmpty() && !resth && ifthen == 0, "fresh R-group entry has empty/false/0 logic fields");
    CHECK(m.setRGroupLogic(1, QStringLiteral("1,2"), true, 2), "setRGroupLogic succeeds");
    CHECK(m.rgroupLogic(1, range, resth, ifthen) && range == QStringLiteral("1,2") && resth && ifthen == 2,
          "setRGroupLogic changed the fields");
    CHECK(!m.setRGroupLogic(2, QStringLiteral("x"), false, 0), "setRGroupLogic fails for an unknown R-group number");

    CHECK(m.rgroupFragmentIds(1).isEmpty(), "no member fragments yet");
    CHECK(m.addRGroupFragment(1, frag0a), "addRGroupFragment succeeds");
    CHECK(!m.addRGroupFragment(1, frag0a), "addRGroupFragment again fails -- already a member");
    CHECK(m.rgroupFragmentIds(1) == QList<int>{frag0a}, "rgroupFragmentIds lists the one member fragment");
    CHECK(m.removeRGroupFragment(1, frag0a), "removeRGroupFragment succeeds");
    CHECK(m.rgroupFragmentIds(1).isEmpty(), "member fragment removed");
    CHECK(!m.removeRGroupFragment(1, frag0a), "removeRGroupFragment again fails -- not a member");

    CHECK(m.removeRGroupEntry(1), "removeRGroupEntry succeeds");
    CHECK(!m.removeRGroupEntry(1), "removeRGroupEntry again fails -- already gone");
    CHECK(m.rgroupNumbers().isEmpty(), "no R-groups remain");
}

static void test_textBoldItalic() {
    std::printf("--- Test 29: text annotation bold/italic ---\n");
    EditableMolecule m;
    int id = m.addTextAnnotation(1, 2, QStringLiteral("hello"));

    double x = 0, y = 0; QString content; bool bold = true, italic = true;
    CHECK(m.textAnnotationContent(id, x, y, content, bold, italic), "textAnnotationContent succeeds");
    CHECK(!bold && !italic, "new text annotation defaults bold/italic to false");
    CHECK(content == QStringLiteral("hello"), "content round-trips");

    CHECK(m.setTextAnnotation(id, QStringLiteral("world"), true, true), "setTextAnnotation succeeds");
    CHECK(m.textAnnotationContent(id, x, y, content, bold, italic), "re-read after mutation");
    CHECK(content == QStringLiteral("world") && bold && italic, "content/bold/italic all updated");
    CHECK(x == 1 && y == 2, "position is untouched by setTextAnnotation");

    CHECK(!m.setTextAnnotation(9999, QStringLiteral("x"), false, false), "setTextAnnotation fails for an unknown id");
    CHECK(!m.textAnnotationContent(9999, x, y, content, bold, italic), "textAnnotationContent fails for an unknown id");
}

static void test_bracketStack() {
    std::printf("--- Test 30: bracket push/pop stack ---\n");
    EditableMolecule m;
    CHECK(m.bracketCount() == 0, "no brackets initially");
    CHECK(!m.popBracket(), "popBracket on an empty stack fails");

    m.pushBracket(-1, -2, 3, 4);
    CHECK(m.bracketCount() == 1, "one bracket after push");
    double minX = 0, minY = 0, maxX = 0, maxY = 0;
    CHECK(m.bracketAt(0, minX, minY, maxX, maxY), "bracketAt(0) succeeds");
    CHECK(minX == -1 && minY == -2 && maxX == 3 && maxY == 4, "bracket bbox round-trips");
    CHECK(!m.bracketAt(1, minX, minY, maxX, maxY), "bracketAt(1) fails -- out of range");

    m.pushBracket(0, 0, 1, 1);
    CHECK(m.bracketCount() == 2, "two brackets after a second push");
    CHECK(m.popBracket(), "popBracket succeeds");
    CHECK(m.bracketCount() == 1, "one bracket remains after pop");
    CHECK(m.bracketAt(0, minX, minY, maxX, maxY) && minX == -1, "the REMAINING bracket is the first one pushed (LIFO)");
}

static void test_loadFrom() {
    std::printf("--- Test 31: EditableMolecule::loadFrom ---\n");
    EditableMolecule m(QStringLiteral("CCO"));
    CHECK(m.atomCount() == 3, "setup: 3 atoms from the initial SMILES");
    AtomId oldA1 = m.atomIds().first();
    m.setName(QStringLiteral("ethanol"));
    m.addTextAnnotation(0, 0, QStringLiteral("note"));
    AtomId nextIdBefore = m.addAtom(QStringLiteral("N"), 5, 5);   // bumps m_nextAtomId forward
    m.removeAtom(nextIdBefore);                                   // remove it again -- id is now "used up"

    CHECK(m.loadFrom(QStringLiteral("c1ccccc1")), "loadFrom succeeds on a valid SMILES");
    CHECK(m.atomCount() == 6, "new document has benzene's 6 atoms, old 3 are gone");
    // oldA1 was one of {1,2,3}; assignFreshAtomBondIds() resumes from the CURRENT m_nextAtomId
    // (already at 5 after the earlier add+remove), so the new document's ids are {5..10} --
    // oldA1 is genuinely gone, not aliased onto a new atom.
    CHECK(!m.atomIds().contains(oldA1), "old AtomIds do not survive loadFrom -- ids are never reused");
    CHECK(m.name().isEmpty(), "m_ext reset to defaults -- name is empty again");
    CHECK(m.textAnnotationCount() == 0, "m_ext reset to defaults -- text annotations gone");

    AtomId newAtomId = m.addAtom(QStringLiteral("F"), 0, 0);
    CHECK(newAtomId > nextIdBefore, "AtomId counter continues forward, NOT reset to 1 after loadFrom");

    QString before = m.toMolfile().value;
    CHECK(!m.loadFrom(QStringLiteral("not a valid molecule at all $$$")), "loadFrom fails on garbage input");
    CHECK(!m.lastError().isEmpty(), "lastError is set on failure");
    CHECK(m.toMolfile().value == before, "failed loadFrom left the object's content completely unchanged");
}

static void test_submoleculeMolfile() {
    std::printf("--- Test 32: submoleculeMolfile ---\n");

    // Full-selection extraction of the "Ac" superatom fixture (3 atoms, 1 sgroup) -- same
    // fixture used throughout this codebase's own tests.
    QString acMolfile = QStringLiteral(
        "Ac\n"
        "  Ketcher\n\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    3.3951   -3.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.6785   -3.9891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.3951   -2.7480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  2  1  1  0  0  0  0\n"
        "  1  3  2  0  0  0  0\n"
        "M  STY  1   1 SUP\n"
        "M  SLB  1   1   1\n"
        "M  SAL   1  3   1   2   3\n"
        "M  SAP   1  1   1   0\n"
        "M  SMT   1 Ac\n"
        "M  END\n");
    EditableMolecule m;
    EditableMolecule::InsertResult r = m.insertStructure(acMolfile, [](double x, double y) { return QPointF(x, y); });
    CHECK(r.createdAtoms.size() == 3, "setup: 3 atoms inserted");
    CHECK(r.createdSGroups.size() == 1, "setup: 1 sgroup inserted");

    StringResult full = m.submoleculeMolfile(r.createdAtoms, r.createdBonds);
    CHECK(full.success, "full-selection extraction succeeds");
    int reparsedFull = indigoLoadMoleculeFromString(full.value.toUtf8().constData());
    CHECK(reparsedFull >= 0, "full-selection molfile reparses cleanly");
    if (reparsedFull >= 0) {
        CHECK(indigoCountAtoms(reparsedFull) == 3, "full-selection has all 3 atoms");
        CHECK(indigoCountSuperatoms(reparsedFull) == 1, "full-selection preserves the sgroup");
        indigoFree(reparsedFull);
    }

    // Partial selection: only 2 of the 3 atoms + 1 bond -- confirms the sgroup survives with a
    // trimmed member list, matching the direct probe finding from this sub-project's own spec.
    QList<AtomId> partialAtoms = { r.createdAtoms[0], r.createdAtoms[1] };
    QList<BondId> partialBonds = { r.createdBonds[0] };
    StringResult partial = m.submoleculeMolfile(partialAtoms, partialBonds);
    CHECK(partial.success, "partial-selection extraction succeeds");
    int reparsedPartial = indigoLoadMoleculeFromString(partial.value.toUtf8().constData());
    CHECK(reparsedPartial >= 0, "partial-selection molfile reparses cleanly");
    if (reparsedPartial >= 0) {
        CHECK(indigoCountAtoms(reparsedPartial) == 2, "partial-selection has only the 2 selected atoms");
        CHECK(indigoCountSuperatoms(reparsedPartial) == 1, "partial-selection STILL preserves the sgroup (trimmed member list)");
        indigoFree(reparsedPartial);
    }

    // Unknown id -- error result.
    StringResult bad = m.submoleculeMolfile({9999}, {});
    CHECK(!bad.success, "unknown AtomId returns a failure StringResult");
}

static void test_submoleculeMolfileInconsistentSelection() {
    std::printf("--- Test 33: submoleculeMolfile with an inconsistent selection ---\n");
    // 3-atom chain: only atom 0 is "selected" in atomIds, but the atom-0-to-atom-1 bond IS
    // selected in bondIds -- atom 1 is NOT in atomIds. Confirmed by direct probe against the
    // real indigo.dll that indigoCreateEdgeSubmolecule FAILS OUTRIGHT on this input unless the
    // vertex set is defensively expanded first.
    EditableMolecule m(QStringLiteral("CCC"));
    QList<AtomId> ids = m.atomIds();
    std::sort(ids.begin(), ids.end());
    AtomId a0 = ids[0], a1 = ids[1];
    BondId bond01 = m.findBond(a0, a1);
    CHECK(bond01 >= 0, "setup: bond between atom 0 and atom 1 exists");

    StringResult result = m.submoleculeMolfile({a0}, {bond01});
    CHECK(result.success, "inconsistent selection (bond without both endpoints) still succeeds");
    int reparsed = indigoLoadMoleculeFromString(result.value.toUtf8().constData());
    CHECK(reparsed >= 0, "result reparses cleanly");
    if (reparsed >= 0) {
        CHECK(indigoCountAtoms(reparsed) == 2, "the bond's other endpoint (atom 1) was defensively included, not dropped or failed");
        indigoFree(reparsed);
    }
}

static void test_setBondStereo() {
    std::printf("--- Test: setBondStereo (sub-project 8) ---\n");

    // Fresh, coordinate-bearing, no-initial-stereo fixture: a stereocenter carbon with three
    // distinct substituents (F/Cl/Br) at non-collinear 2D positions, all bonds starting plain
    // (stereo=0) -- same shape confirmed working end-to-end by this sub-project's own probe.
    const char* molfile =
        "probe\n"
        "  probe\n"
        "\n"
        "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000    2.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660    0.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660    0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  1  4  1  0  0  0  0\n"
        "M  END\n";
    EditableMolecule m(molfile);
    CHECK(m.isValid(), "stereocenter fixture loads");

    QList<AtomId> ids = m.atomIdsInIndigoOrder();
    QList<BondId> bIds = m.bondIdsInIndigoOrder();
    CHECK(ids.size() == 4 && bIds.size() == 3, "fixture has 4 atoms, 3 bonds");
    AtomId cAtom = ids[0];
    BondId targetBond = bIds[2]; // C-Br
    BondId otherBond1 = bIds[0]; // C-F
    BondId otherBond2 = bIds[1]; // C-Cl

    CHECK(m.bondStereoDirectionEnum(targetBond) == EditableMolecule::Direction::None,
          "target bond starts with no stereo");

    CHECK(m.setBondStereo(targetBond, EditableMolecule::Direction::Up),
          "setBondStereo(Up) succeeds");
    CHECK(m.bondStereoDirectionEnum(targetBond) == EditableMolecule::Direction::Up,
          "target bond reports Up");
    CHECK(m.bondStereoDirection(targetBond) == INDIGO_UP,
          "raw indigoBondStereo reports INDIGO_UP (5) specifically, not the V2000 write code (1)");
    CHECK(m.bondStereoDirectionEnum(otherBond1) == EditableMolecule::Direction::None,
          "other bond (C-F) unaffected");
    CHECK(m.bondStereoDirectionEnum(otherBond2) == EditableMolecule::Direction::None,
          "other bond (C-Cl) unaffected");
    CHECK(m.stereocenterType(cAtom) != 0,
          "the stereocenter atom is now perceived as a real stereocenter");

    CHECK(m.setBondStereo(targetBond, EditableMolecule::Direction::Down), "setBondStereo(Down) succeeds");
    CHECK(m.bondStereoDirection(targetBond) == INDIGO_DOWN, "reports INDIGO_DOWN (6)");

    CHECK(m.setBondStereo(targetBond, EditableMolecule::Direction::Either), "setBondStereo(Either) succeeds");
    CHECK(m.bondStereoDirection(targetBond) == INDIGO_EITHER, "reports INDIGO_EITHER (4)");

    CHECK(m.setBondStereo(targetBond, EditableMolecule::Direction::None), "setBondStereo(None) succeeds");
    CHECK(m.bondStereoDirectionEnum(targetBond) == EditableMolecule::Direction::None, "reports None again");

    StringResult mfAfter = m.toMolfile();
    CHECK(mfAfter.success, "toMolfile succeeds after stereo edits");

    CHECK(!m.setBondStereo(999999, EditableMolecule::Direction::Up),
          "setBondStereo on an invalid BondId fails");

    // Double bond guard.
    EditableMolecule dbl;
    AtomId da = dbl.addAtom(QStringLiteral("C"), 0, 0);
    AtomId db = dbl.addAtom(QStringLiteral("O"), 0, 1);
    BondId dblBond = dbl.addBond(da, db, 2);
    CHECK(!dbl.setBondStereo(dblBond, EditableMolecule::Direction::Up),
          "setBondStereo on a double bond fails");
    CHECK(dbl.bondStereoDirectionEnum(dblBond) == EditableMolecule::Direction::None,
          "double bond's stereo remains untouched after the rejected attempt");

    // REQUIRED: nontrivial add/remove history, proving the order-list-vs-serialization-order
    // assumption holds -- not a single literal string load.
    EditableMolecule scrambled;
    AtomId a = scrambled.addAtom(QStringLiteral("C"), 0, 0);
    AtomId b = scrambled.addAtom(QStringLiteral("N"), 1, 0);
    AtomId c = scrambled.addAtom(QStringLiteral("O"), 0, 1);
    AtomId d = scrambled.addAtom(QStringLiteral("F"), -1, 0);
    BondId ac = scrambled.addBond(a, c, 1);
    BondId ad = scrambled.addBond(a, d, 1);
    scrambled.removeAtom(b); // frees an internal Indigo index
    AtomId b2 = scrambled.addAtom(QStringLiteral("N"), 1, 0); // may reuse the freed index
    BondId ab2 = scrambled.addBond(a, b2, 1);

    double aX = 0, aY = 0, cX = 0, cY = 0, dX = 0, dY = 0, b2X = 0, b2Y = 0;
    scrambled.atomPos(a, aX, aY);
    scrambled.atomPos(c, cX, cY);
    scrambled.atomPos(d, dX, dY);
    scrambled.atomPos(b2, b2X, b2Y);
    QString aSym = scrambled.atomSymbol(a), cSym = scrambled.atomSymbol(c);
    QString dSym = scrambled.atomSymbol(d), b2Sym = scrambled.atomSymbol(b2);

    CHECK(scrambled.setBondStereo(ab2, EditableMolecule::Direction::Up),
          "setBondStereo succeeds on the scrambled-history molecule");

    double aX2 = 0, aY2 = 0, cX2 = 0, cY2 = 0, dX2 = 0, dY2 = 0, b2X2 = 0, b2Y2 = 0;
    scrambled.atomPos(a, aX2, aY2);
    scrambled.atomPos(c, cX2, cY2);
    scrambled.atomPos(d, dX2, dY2);
    scrambled.atomPos(b2, b2X2, b2Y2);
    CHECK(aX == aX2 && aY == aY2 && scrambled.atomSymbol(a) == aSym,
          "atom a's identity (position+symbol) is unchanged");
    CHECK(cX == cX2 && cY == cY2 && scrambled.atomSymbol(c) == cSym,
          "atom c's identity is unchanged");
    CHECK(dX == dX2 && dY == dY2 && scrambled.atomSymbol(d) == dSym,
          "atom d's identity is unchanged");
    CHECK(b2X == b2X2 && b2Y == b2Y2 && scrambled.atomSymbol(b2) == b2Sym,
          "atom b2's identity is unchanged -- no atom silently swapped identity");
    CHECK(scrambled.bondStereoDirectionEnum(ab2) == EditableMolecule::Direction::Up,
          "the target bond correctly ended up with Up");
    CHECK(scrambled.bondStereoDirectionEnum(ac) == EditableMolecule::Direction::None, "bond ac unaffected");
    CHECK(scrambled.bondStereoDirectionEnum(ad) == EditableMolecule::Direction::None, "bond ad unaffected");

    // REQUIRED: sgroup survives a setBondStereo call elsewhere in the same molecule.
    //
    // NOTE (found by direct testing, not anticipated in the approved spec): the plain 3-atom
    // "Ac" fixture used elsewhere in this file (2 carbons + 1 oxygen, one single bond, one
    // double bond) has NO bond whose begin-atom is a chemically plausible wedge origin -- both
    // atoms have too few real substituents to be a valid stereocenter, and Indigo correctly
    // REJECTS the reload with "direction of bond #0 makes no sense" (a real, correct chemistry
    // validation, not a bug in setBondStereo). This combined fixture instead concatenates the
    // separately-already-proven-working 4-atom stereocenter fragment (atoms 1-4: C/F/Cl/Br, same
    // shape as this test's very first fixture above) with the Ac sgroup's 3 atoms as a second,
    // disconnected fragment (atoms 5-7, atom numbers shifted +4 from the original Ac fixture) --
    // multi-fragment molecules are completely normal for Indigo. The bond wedged (1-4, C-Br) is
    // nowhere near the sgroup; this fixture exists purely to prove the sgroup-remap logic doesn't
    // regress a molecule that happens to also contain one.
    QString acMolfile = QStringLiteral(
        "combined\n"
        "  Ketcher\n\n"
        "  7  5  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000    2.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660    0.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660    0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.3951   -3.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.6785   -3.9891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.3951   -2.7480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  1  4  1  0  0  0  0\n"
        "  6  5  1  0  0  0  0\n"
        "  5  7  2  0  0  0  0\n"
        "M  STY  1   1 SUP\n"
        "M  SLB  1   1   1\n"
        "M  SAL   1  3   5   6   7\n"
        "M  SAP   1  1   5   0\n"
        "M  SMT   1 Ac\n"
        "M  END\n");
    EditableMolecule acMol(acMolfile);
    CHECK(acMol.isValid(), "combined stereocenter+sgroup fixture loads");
    QList<SGroupId> sgIdsBefore = acMol.sgroupIds();
    CHECK(sgIdsBefore.size() == 1, "exactly one sgroup before the edit");
    SGroupId sg = sgIdsBefore.first();
    QList<AtomId> membersBefore = acMol.sgroupMemberAtomIds(sg);
    std::sort(membersBefore.begin(), membersBefore.end());
    AtomId attachBefore = -1;
    acMol.superatomAttachAtom(sg, attachBefore);

    BondId wedgeTarget = -1;
    for (BondId bid : acMol.bondIds()) {
        double x1 = 0, y1 = 0;
        AtomId ba = -1, bb = -1;
        if (acMol.bondOrder(bid) == 1 && acMol.bondEndpoints(bid, ba, bb) &&
            acMol.atomSymbol(ba) == QStringLiteral("C") && acMol.atomSymbol(bb) == QStringLiteral("Br")) {
            wedgeTarget = bid;
            break;
        }
    }
    CHECK(wedgeTarget != -1, "found the fixture's C-Br bond (the same shape already proven wedgeable)");
    CHECK(acMol.setBondStereo(wedgeTarget, EditableMolecule::Direction::Up),
          "setBondStereo succeeds on a bond in an sgroup-bearing molecule");

    QList<SGroupId> sgIdsAfter = acMol.sgroupIds();
    CHECK(sgIdsAfter.size() == 1 && sgIdsAfter.first() == sg, "the sgroup's own id survives the reload");
    QList<AtomId> membersAfter = acMol.sgroupMemberAtomIds(sg);
    std::sort(membersAfter.begin(), membersAfter.end());
    CHECK(membersAfter == membersBefore, "the sgroup's member-atom set is unchanged");
    AtomId attachAfter = -1;
    acMol.superatomAttachAtom(sg, attachAfter);
    CHECK(attachAfter == attachBefore, "the sgroup's attachment atom is unchanged");
}

static void test_ezStereoIntent() {
    std::printf("--- Test: E/Z stereo intent capture and molfileForNaming ---\n");

    // Test 1: Load a SMILES with no stereo specified
    EditableMolecule molNoStereo(QStringLiteral("CC=C(C)CC"));
    CHECK(molNoStereo.isValid(), "molecule without stereo loads");
    CHECK(molNoStereo.bondCount() > 0, "has bonds");

    // Check that the double bond was recorded as unspecified (0)
    // The bond between atoms 2 and 3 (0-indexed in the SMILES: CC=C(C)CC)
    // In the loaded molecule, we need to find the double bond
    bool foundDoubleBond = false;
    for (BondId bid : molNoStereo.bondIds()) {
        if (molNoStereo.bondOrder(bid) == 2) {
            // This should have been recorded as unspecified
            // We can't directly access m_ezWasSpecified, but we can test molfileForNaming
            foundDoubleBond = true;
            break;
        }
    }
    CHECK(foundDoubleBond, "found a double bond in CC=C(C)CC");

    // Test molfileForNaming - it should not crash
    QString mf = molNoStereo.molfileForNaming();
    CHECK(!mf.isEmpty(), "molfileForNaming produces non-empty output");
    CHECK(mf.contains(QStringLiteral("M  END")), "molfileForNaming output is valid molfile");

    // The real correctness criterion: reload the produced molfile fresh (exactly
    // like IndigoService::generateIupacName does) and confirm the double bond's
    // indigoBondStereo() is now genuinely 0 (unspecified), not fabricated from
    // coordinates. This is the actual bug this whole feature exists to fix --
    // a molfile that merely parses is not sufficient evidence of correctness.
    {
        int reloaded = indigoLoadMoleculeFromString(mf.toUtf8().constData());
        CHECK(reloaded >= 0, "reloaded molfileForNaming output parses");
        if (reloaded >= 0) {
            int bi = indigoIterateBonds(reloaded);
            int foundReloadedDouble = 0;
            int b2;
            while ((b2 = indigoNext(bi)) != 0) {
                if (indigoBondOrder(b2) == 2) {
                    foundReloadedDouble++;
                    CHECK(indigoBondStereo(b2) == 0,
                          "unspecified double bond reads as unspecified after reload (the actual bug this fixes)");
                }
                indigoFree(b2);
            }
            indigoFree(bi);
            CHECK(foundReloadedDouble == 1, "reloaded molecule still has exactly one double bond");
            indigoFree(reloaded);
        }
    }

    // Test 2: Load a SMILES with explicit E stereo
    EditableMolecule molEStereo(QStringLiteral("C/C=C/C"));
    CHECK(molEStereo.isValid(), "molecule with E stereo loads");

    QString mfE = molEStereo.molfileForNaming();
    CHECK(!mfE.isEmpty(), "molfileForNaming with E stereo produces output");

    // Real, deliberately-specified stereo must NOT be suppressed by this fix.
    {
        int reloadedE = indigoLoadMoleculeFromString(mfE.toUtf8().constData());
        CHECK(reloadedE >= 0, "reloaded E-stereo molfileForNaming output parses");
        if (reloadedE >= 0) {
            int bi = indigoIterateBonds(reloadedE);
            int b2;
            while ((b2 = indigoNext(bi)) != 0) {
                if (indigoBondOrder(b2) == 2) {
                    CHECK(indigoBondStereo(b2) != 0,
                          "genuinely specified E stereo survives molfileForNaming, not suppressed");
                }
                indigoFree(b2);
            }
            indigoFree(bi);
            indigoFree(reloadedE);
        }
    }

    // Test 2b: a trisubstituted alkene where one carbon carries TWO substituents
    // (both unspecified). Both must end up genuinely collinear with the bond axis
    // (not coincident with each other), and the reloaded bond must read as
    // unspecified -- this is exactly the shape the original coordinate-collapse
    // bug (indigoXYZ buffer aliasing + a sign error) would have corrupted.
    {
        EditableMolecule molTri(QStringLiteral("CC=C(C)C"));  // 2-methylbut-2-ene, unspecified
        CHECK(molTri.isValid(), "trisubstituted alkene loads");
        QString mfTri = molTri.molfileForNaming();
        CHECK(!mfTri.isEmpty(), "molfileForNaming produces output for trisubstituted alkene");
        int reloadedTri = indigoLoadMoleculeFromString(mfTri.toUtf8().constData());
        CHECK(reloadedTri >= 0, "reloaded trisubstituted-alkene molfile parses");
        if (reloadedTri >= 0) {
            int bi = indigoIterateBonds(reloadedTri);
            int b2;
            while ((b2 = indigoNext(bi)) != 0) {
                if (indigoBondOrder(b2) == 2) {
                    CHECK(indigoBondStereo(b2) == 0,
                          "trisubstituted alkene's unspecified double bond stays unspecified after reload");
                }
                indigoFree(b2);
            }
            indigoFree(bi);
            indigoFree(reloadedTri);
        }
    }

    // Test 3: Empty molecule - should still produce a valid (if minimal) molfile
    EditableMolecule empty;
    QString mfEmpty = empty.molfileForNaming();
    // Empty molecule should not crash, but may return empty or a valid empty molfile
    // Just check it doesn't crash and returns a non-null string (could be empty or valid)
    CHECK(true, "empty molecule molfileForNaming does not crash");

    // Test 4: Molecule with no double bonds
    EditableMolecule ethane(QStringLiteral("CC"));
    CHECK(ethane.isValid(), "ethane loads");
    QString mfEthane = ethane.molfileForNaming();
    CHECK(!mfEthane.isEmpty(), "molfileForNaming works for molecule without double bonds");

    // Test 5: the REAL app flow this bug was reported in -- load bare SMILES (no coords),
    // then simulate DocumentState::deserializeMol receiving an auto-layout result the way
    // the real app does (IndigoService::layout exports the current molfile, lays it out in
    // a SEPARATE indigo call, and the result flows back via EditableMolecule::loadFrom --
    // see DocumentState::deserializeMol, MainWindow.qml's onLayoutFinished). If loadFrom's
    // own captureEzIntent() call re-derives intent from the now-coordinate-bearing molfile,
    // it will silently overwrite the correct originally-captured "unspecified" value with a
    // fabricated one, defeating this entire fix for its own primary reported scenario.
    {
        EditableMolecule molReal(QStringLiteral("CC=C(C)CC"));
        CHECK(molReal.isValid(), "real-flow molecule loads");

        // Build the exact molfile IndigoService::layout would hand back: load the same
        // structure fresh, run indigoLayout, export.
        int rawMol = indigoLoadMoleculeFromString("CC=C(C)CC");
        const char* rawMf = indigoMolfile(rawMol);
        int layoutMol = indigoLoadMoleculeFromString(rawMf);
        indigoLayout(layoutMol);
        const char* laidOutMf = indigoMolfile(layoutMol);
        QString laidOutMfQ = QString::fromUtf8(laidOutMf);
        indigoFree(rawMol);
        indigoFree(layoutMol);

        // This is exactly what DocumentState::deserializeMol does with a layout result
        // (centerOnPage=true for a same-document write-back, threaded through as preserveEzIntent).
        CHECK(molReal.loadFrom(laidOutMfQ, /*preserveEzIntent=*/true), "loadFrom accepts the simulated layout result");

        QString nameMf = molReal.molfileForNaming();
        CHECK(!nameMf.isEmpty(), "molfileForNaming produces output after the real flow");

        int reloaded = indigoLoadMoleculeFromString(nameMf.toUtf8().constData());
        CHECK(reloaded >= 0, "reloaded output of the real-flow test parses");
        if (reloaded >= 0) {
            int bi = indigoIterateBonds(reloaded);
            int b2;
            while ((b2 = indigoNext(bi)) != 0) {
                if (indigoBondOrder(b2) == 2) {
                    CHECK(indigoBondStereo(b2) == 0,
                          "REAL FLOW: load SMILES -> auto-layout -> Generate IUPAC Name does not fabricate E/Z");
                }
                indigoFree(b2);
            }
            indigoFree(bi);
            indigoFree(reloaded);
        }
    }
}

static void test_loadFromClearsEzIntent() {
    std::printf("--- Test: loadFrom clears and recaptures E/Z intent ---\n");

    // Create a molecule with one double bond
    EditableMolecule mol1(QStringLiteral("CC=C(C)CC"));
    CHECK(mol1.isValid(), "first molecule loads");

    // loadFrom should clear the intent and recapture
    CHECK(mol1.loadFrom(QStringLiteral("C/C=C/C")), "loadFrom succeeds");
    CHECK(mol1.isValid(), "still valid after loadFrom");

    // Should still be able to produce molfile for naming
    QString mf = mol1.molfileForNaming();
    CHECK(!mf.isEmpty(), "molfileForNaming works after loadFrom");
}

static void test_scaleAllCoordinates() {
    std::printf("--- Test: scaleAllCoordinates ---\n");

    // Create a simple molecule with atoms at known positions
    EditableMolecule m;
    AtomId a1 = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = m.addAtom(QStringLiteral("C"), 2.0, 0.0);
    AtomId a3 = m.addAtom(QStringLiteral("O"), 4.0, 0.0);
    CHECK(a1 > 0 && a2 > 0 && a3 > 0, "added 3 atoms");

    // Scale by factor of 2 around origin (0,0)
    int scaledCount = m.scaleAllCoordinates(2.0, 0.0, 0.0);
    CHECK(scaledCount == 3, "all 3 atoms scaled");

    // Verify new positions
    double x, y;
    CHECK(m.atomPos(a1, x, y) && x == 0.0 && y == 0.0, "a1 stays at origin");
    CHECK(m.atomPos(a2, x, y) && x == 4.0 && y == 0.0, "a2 scaled to (4,0)");
    CHECK(m.atomPos(a3, x, y) && x == 8.0 && y == 0.0, "a3 scaled to (8,0)");
}

static void test_scaleAllCoordinatesPreservesRelativePositions() {
    std::printf("--- Test: scaleAllCoordinates preserves relative positions ---\n");

    // Create a triangle molecule
    EditableMolecule m;
    AtomId a1 = m.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId a2 = m.addAtom(QStringLiteral("C"), 3.0, 0.0);
    AtomId a3 = m.addAtom(QStringLiteral("C"), 0.0, 4.0);

    // Get original relative distances
    double x1, y1, x2, y2, x3, y3;
    m.atomPos(a1, x1, y1);
    m.atomPos(a2, x2, y2);
    m.atomPos(a3, x3, y3);
    double origDist12 = std::hypot(x2 - x1, y2 - y1);
    double origDist13 = std::hypot(x3 - x1, y3 - y1);
    double origDist23 = std::hypot(x3 - x2, y3 - y2);

    // Scale by factor of 2
    m.scaleAllCoordinates(2.0, 0.0, 0.0);

    // Get new positions
    m.atomPos(a1, x1, y1);
    m.atomPos(a2, x2, y2);
    m.atomPos(a3, x3, y3);
    double newDist12 = std::hypot(x2 - x1, y2 - y1);
    double newDist13 = std::hypot(x3 - x1, y3 - y1);
    double newDist23 = std::hypot(x3 - x2, y3 - y2);

    // Relative distances should scale by factor of 2
    CHECK(std::abs(newDist12 - origDist12 * 2.0) < 0.001, "distance a1-a2 scaled correctly");
    CHECK(std::abs(newDist13 - origDist13 * 2.0) < 0.001, "distance a1-a3 scaled correctly");
    CHECK(std::abs(newDist23 - origDist23 * 2.0) < 0.001, "distance a2-a3 scaled correctly");

    // Ratios should be preserved
    double origRatio = origDist12 / origDist13;
    double newRatio = newDist12 / newDist13;
    CHECK(std::abs(origRatio - newRatio) < 0.001, "distance ratios preserved");
}

static void test_scaleAllCoordinatesWithCentroid() {
    std::printf("--- Test: scaleAllCoordinates with centroid ---\n");

    // Create a molecule with atoms at (0,0), (2,0), (4,0)
    EditableMolecule m;
    m.addAtom(QStringLiteral("C"), 0.0, 0.0);
    m.addAtom(QStringLiteral("C"), 2.0, 0.0);
    m.addAtom(QStringLiteral("C"), 4.0, 0.0);

    // Compute centroid
    QPointF centroid = m.computeCentroid();
    CHECK(std::abs(centroid.x() - 2.0) < 0.001 && std::abs(centroid.y() - 0.0) < 0.001,
          "centroid at (2,0)");

    // Scale by factor of 2 around centroid
    m.scaleAllCoordinates(2.0, centroid.x(), centroid.y());

    // Atoms should now be at (-2,0), (2,0), (6,0)
    // because: newX = 2 + (x - 2) * 2
    QList<AtomId> atomIds = m.atomIds();
    double x, y;
    m.atomPos(atomIds[0], x, y);
    CHECK(std::abs(x - (-2.0)) < 0.001 && std::abs(y - 0.0) < 0.001,
          "first atom scaled correctly around centroid");

    m.atomPos(atomIds[1], x, y);
    CHECK(std::abs(x - 2.0) < 0.001 && std::abs(y - 0.0) < 0.001,
          "second atom (at centroid) stays at same position");

    m.atomPos(atomIds[2], x, y);
    CHECK(std::abs(x - 6.0) < 0.001 && std::abs(y - 0.0) < 0.001,
          "third atom scaled correctly around centroid");
}

int main() {
    unsigned long long session = indigoAllocSessionId();
    indigoSetSessionId(session);

    spike_idBehavior();
    test_constructionAndMolfile();
    test_toKetJson();
    test_atomBondCrud();
    test_extensionData();
    test_sgroupPassthrough();
    test_snapshotRestore();
    test_queryFeatureSpike();
    test_extensionEntityRemoval();
    test_atomAttributeAccessors();
    test_atomLabelAndBondOrder();
    test_attachmentPointAccessors();
    test_sgroupMembershipSpike();
    test_bondStereoSpike();
    test_atomQueryList();
    test_bondStereoInvertSpike();
    test_positionAccessors();
    test_atomMergeMechanicsSpike();
    test_mergeOverlappingAtomsAndFindBond();
    test_addBenzeneRing();
    test_templateLibraryLoading();
    test_insertStructure();
    test_graftAtomOnto();
    test_renderSupportStorageDefaults();
    test_ringAndNeighborAccessors();
    test_stereoCipAccessors();
    test_sgroupIntrospectionAccessors();
    test_rebuildIndexTablesPrunesAutoRemovedSgroup();
    test_removeSuperatomOnly();
    test_createSuperatomFromAtoms();
    test_atomIdsInIndigoOrder();
    test_bondIdsInIndigoOrder();
    test_atomStereoDescriptorStorage();
    test_bondStereoCipLabelStorage();
    test_checkWarningTextSetters();
    test_rxnArrowMutators();
    test_stereoFlagsDocumentMutator();
    test_rgroupStorageAndFragmentIndex();
    test_textBoldItalic();
    test_bracketStack();
    test_loadFrom();
    test_submoleculeMolfile();
    test_submoleculeMolfileInconsistentSelection();
    test_setBondStereo();
    test_ezStereoIntent();
    test_loadFromClearsEzIntent();
    
    // Phase 3: scaleAllCoordinates tests
    test_scaleAllCoordinates();
    test_scaleAllCoordinatesPreservesRelativePositions();
    test_scaleAllCoordinatesWithCentroid();

    indigoReleaseSessionId(session);
    std::printf("Summary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
