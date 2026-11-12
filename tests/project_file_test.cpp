// tests/project_file_test.cpp
// Tests for ProjectFile: save/load round-trip, atomic-save failure safety,
// and load validation (regression coverage for the data-integrity bugs
// found by the five-phase QtSketch review and fixed in this pass).
#include <cstdio>
#include <cmath>
#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include "app/molecule/ProjectFile.h"
#include "app/molecule/EditableMolecule.h"
#include "app/molecule/DrawingStyle.h"
#include "app/molecule/PageSettings.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

static void test_saveLoadRoundTrip() {
    std::printf("--- Test: saveProject/loadProject round-trip ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("roundtrip.qtsketch");

    EditableMolecule mol("C");
    DrawingStyle style(2.0, 0.8, "Arial", 12, false);
    PageSettings page(200.0, 280.0, 5, 5, 5, 5);

    CHECK(ProjectFile::saveProject(path, mol, style, page), "save succeeds");
    CHECK(QFile::exists(path), "file exists after save");

    EditableMolecule loaded;
    DrawingStyle loadedStyle;
    PageSettings loadedPage;
    QString err = ProjectFile::loadProject(path, loaded, loadedStyle, loadedPage);
    CHECK(err.isEmpty(), "load succeeds with no error");
    CHECK(loadedStyle.bondLength == 2.0, "bondLength round-trips");
    CHECK(loadedStyle.fontFamily == "Arial", "fontFamily round-trips");
    CHECK(loadedPage.widthMm == 200.0, "page widthMm round-trips");
}

static void test_atomicSavePreservesOriginalOnFailure() {
    std::printf("--- Test: atomic save never destroys the original on failure ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("protected.qtsketch");

    EditableMolecule mol("C");
    DrawingStyle style;
    PageSettings page;
    CHECK(ProjectFile::saveProject(path, mol, style, page), "initial save succeeds");

    QFile f(path);
    qint64 originalSize = f.size();
    CHECK(originalSize > 0, "original file is non-empty");

    // Simulate a write-target that cannot be created (invalid directory) --
    // the original file at `path` must remain completely untouched.
    QString badPath = dir.filePath("does_not_exist_dir/nested.qtsketch");
    bool badSaveResult = ProjectFile::saveProject(badPath, mol, style, page);
    CHECK(!badSaveResult, "save to an unwritable path reports failure");

    QFile f2(path);
    CHECK(f2.exists(), "original file still exists after the unrelated failed save");
    CHECK(f2.size() == originalSize, "original file is byte-for-byte unchanged");
}

static void test_loadRejectsMissingFormatVersion() {
    std::printf("--- Test: load rejects a file missing formatVersion ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("noversion.qtsketch");
    QFile f(path);
    f.open(QIODevice::WriteOnly);
    f.write("{\"chemistry\": {\"molfile\": \"\"}}");
    f.close();

    EditableMolecule mol;
    DrawingStyle style;
    PageSettings page;
    QString err = ProjectFile::loadProject(path, mol, style, page);
    CHECK(!err.isEmpty(), "load reports an error for missing formatVersion");
}

static void test_loadRejectsMissingChemistry() {
    std::printf("--- Test: load rejects a file missing the chemistry field ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("nochemistry.qtsketch");
    QFile f(path);
    f.open(QIODevice::WriteOnly);
    f.write("{\"formatVersion\": 1}");
    f.close();

    EditableMolecule mol;
    DrawingStyle style;
    PageSettings page;
    QString err = ProjectFile::loadProject(path, mol, style, page);
    CHECK(!err.isEmpty(), "load reports an error for missing chemistry field");
}

static void test_loadRejectsCorruptMolfile() {
    std::printf("--- Test: load rejects a non-empty but unparseable molfile ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("corruptmol.qtsketch");
    QFile f(path);
    f.open(QIODevice::WriteOnly);
    f.write("{\"formatVersion\": 1, \"chemistry\": {\"molfile\": \"this is not a real molfile\"}}");
    f.close();

    EditableMolecule mol;
    DrawingStyle style;
    PageSettings page;
    QString err = ProjectFile::loadProject(path, mol, style, page);
    CHECK(!err.isEmpty(), "load reports an error instead of silently accepting corrupt chemistry");
}

static void test_loadRejectsUnsupportedVersion() {
    std::printf("--- Test: load rejects an unsupported format version ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("futureversion.qtsketch");
    QFile f(path);
    f.open(QIODevice::WriteOnly);
    f.write("{\"formatVersion\": 999, \"chemistry\": {\"molfile\": \"\"}}");
    f.close();

    EditableMolecule mol;
    DrawingStyle style;
    PageSettings page;
    QString err = ProjectFile::loadProject(path, mol, style, page);
    CHECK(!err.isEmpty(), "load reports an error for an unsupported future version");
}

static void test_loadIgnoresNonPositiveScaleValues() {
    std::printf("--- Test: load keeps usable style values when the file holds zero/negative ones ---\n");
    // a zero bondLength from a damaged file must be ignored
    const char* files[] = {
        "{\"formatVersion\": 1, \"chemistry\": {\"molfile\": \"\"}, \"drawingStyle\": {\"bondLength\": 0}}",
        "{\"formatVersion\": 1, \"chemistry\": {\"molfile\": \"\"}, \"drawingStyle\": {\"bondLength\": -2}, \"modelToMmScale\": 0}",
        "{\"formatVersion\": 1, \"chemistry\": {\"molfile\": \"\"}, \"drawingStyle\": {\"bondLength\": 1.5}, \"modelToMmScale\": -1}",
    };
    for (const char* content : files) {
        QTemporaryDir dir;
        QString path = dir.filePath("badscale.qtsketch");
        QFile f(path);
        f.open(QIODevice::WriteOnly);
        f.write(content);
        f.close();

        EditableMolecule mol;
        DrawingStyle style;
        PageSettings page;
        QString err = ProjectFile::loadProject(path, mol, style, page);
        CHECK(err.isEmpty(), "the file still loads");
        CHECK(style.bondLength > 0.0, "bondLength stays positive");
        CHECK(std::isfinite(style.modelToMmScale) && style.modelToMmScale > 0.0,
              "modelToMmScale stays finite and positive");
    }
}

static void test_rgroupsRoundTrip() {
    std::printf("--- Test: R-group entries survive save/load ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("rgroups.qtsketch");

    EditableMolecule mol("CC");
    DrawingStyle style;
    PageSettings page;
    CHECK(mol.addRGroupEntry(1), "add R1");
    CHECK(mol.setRGroupLogic(1, "1-3", true, 0), "set R1 logic");
    CHECK(mol.addRGroupFragment(1, 0), "add R1 fragment");
    CHECK(ProjectFile::saveProject(path, mol, style, page), "save succeeds");

    EditableMolecule loaded;
    DrawingStyle loadedStyle;
    PageSettings loadedPage;
    QString err = ProjectFile::loadProject(path, loaded, loadedStyle, loadedPage);
    CHECK(err.isEmpty(), "load succeeds");
    CHECK(loaded.rgroupNumbers() == QList<int>{1}, "R1 is back after load");
    QString range;
    bool resth = false;
    int ifthen = -1;
    CHECK(loaded.rgroupLogic(1, range, resth, ifthen) && range == "1-3" && resth && ifthen == 0,
          "R1 logic round-trips");
    CHECK(loaded.rgroupFragmentIds(1) == QList<int>{0}, "R1 fragment round-trips");
}

static const char* kTwoSuperatoms =
    "\n  test\n\n"
    "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n"
    "  2  3  1  0\n"
    "  2  4  1  0\n"
    "M  STY  2   1 SUP   2 SUP\n"
    "M  SAL   1  1   3\n"
    "M  SMT   1 OH\n"
    "M  SAL   2  1   4\n"
    "M  SMT   2 NH2\n"
    "M  END\n";

static void test_sgroupStateRoundTrip() {
    std::printf("--- Test: collapsed state and label of a group survive save/load ---\n");
    QTemporaryDir dir;
    QString path = dir.filePath("sgroups.qtsketch");

    // Reload once so the live ids (3, 4) differ from the ids a fresh load hands out (1, 2).
    EditableMolecule mol(kTwoSuperatoms);
    CHECK(mol.loadFrom(kTwoSuperatoms), "reload succeeds");
    QList<SGroupId> ids = mol.sgroupIdsInIndigoOrder();
    CHECK(ids.size() == 2, "two groups before save");
    if (ids.size() != 2) return;
    mol.setSGroupExpanded(ids[1], false);
    mol.setSGroupLabel(ids[1], "Am");

    DrawingStyle style;
    PageSettings page;
    CHECK(ProjectFile::saveProject(path, mol, style, page), "save succeeds");

    EditableMolecule loaded;
    DrawingStyle loadedStyle;
    PageSettings loadedPage;
    QString err = ProjectFile::loadProject(path, loaded, loadedStyle, loadedPage);
    CHECK(err.isEmpty(), "load succeeds");
    QList<SGroupId> loadedIds = loaded.sgroupIdsInIndigoOrder();
    CHECK(loadedIds.size() == 2, "two groups after load");
    if (loadedIds.size() != 2) return;
    CHECK(loaded.sgroupExpanded(loadedIds[0]), "first group stays expanded");
    CHECK(!loaded.sgroupExpanded(loadedIds[1]), "second group stays collapsed");
    CHECK(loaded.sgroupLabel(loadedIds[1]) == "Am", "second group keeps its label");
}

int main() {
    std::printf("Running ProjectFile tests...\n");
    test_rgroupsRoundTrip();
    test_sgroupStateRoundTrip();
    test_saveLoadRoundTrip();
    test_atomicSavePreservesOriginalOnFailure();
    test_loadRejectsMissingFormatVersion();
    test_loadRejectsMissingChemistry();
    test_loadRejectsCorruptMolfile();
    test_loadRejectsUnsupportedVersion();
    test_loadIgnoresNonPositiveScaleValues();
    std::printf("Summary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
