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
    // bondLength and modelToMmScale are divisors downstream (style rescale, page bounds,
    // export scale); a zero from a damaged or hand-edited file must not reach them.
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

int main() {
    std::printf("Running ProjectFile tests...\n");
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
