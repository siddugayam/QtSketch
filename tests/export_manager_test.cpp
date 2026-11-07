// tests/export_manager_test.cpp
// Tests for ExportManager class
#include <cstdio>
#include <QGuiApplication>
#include <QEventLoop>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QFile>
#include <QImage>
#include <QUrl>
#include <QVariantMap>
#include <QTimer>
#include "app/render/ExportManager.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

// Simple molfile for a benzene ring
static const QString benzeneMolfile =
    "\n"
    "  Marvin  02170812212D\n"
    "\n"
    "  6  6  0  0  0  0  0  0  0  0999 V2000\n"
    "   -0.6668    1.7792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.6668    1.7792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.3335    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.6668   -1.7792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.6668   -1.7792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.3335    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  3  4  1  0  0  0  0\n"
    "  4  5  1  0  0  0  0\n"
    "  5  6  1  0  0  0  0\n"
    "  6  1  1  0  0  0  0\n"
    "M  END\n";

// Helper to wait for export to complete with timeout
bool waitForExport(ExportManager& manager, QSignalSpy& spy, int timeoutMs = 10000) {
    QEventLoop loop;
    QTimer timeoutTimer;
    timeoutTimer.setSingleShot(true);
    timeoutTimer.setInterval(timeoutMs);
    
    QObject::connect(&manager, &ExportManager::exportFinished, &loop, &QEventLoop::quit);
    QObject::connect(&timeoutTimer, &QTimer::timeout, &loop, &QEventLoop::quit);
    
    timeoutTimer.start();
    loop.exec();
    
    return spy.count() > 0;
}

int main() {
    int argc = 0;
    char** argv = nullptr;
    QGuiApplication app(argc, argv);

    std::printf("--- ExportManager Tests ---\n\n");

    ExportManager manager;

    // Test 1: availableFormats
    std::printf("--- Test 1: availableFormats ---\n");
    QStringList formats = manager.availableFormats();
    CHECK(formats.contains("pdf"), "availableFormats includes pdf");
    CHECK(formats.contains("svg"), "availableFormats includes svg");
    CHECK(formats.contains("png"), "availableFormats includes png");

    // Test 2: defaultOptions
    std::printf("--- Test 2: defaultOptions ---\n");
    CHECK(!manager.defaultOptions("pdf").isEmpty(), "defaultOptions for pdf is non-empty");
    CHECK(!manager.defaultOptions("svg").isEmpty(), "defaultOptions for svg is non-empty");
    CHECK(!manager.defaultOptions("png").isEmpty(), "defaultOptions for png is non-empty");

    // Test 3-6: File-based tests
    QTemporaryDir tempDir;
    if (!tempDir.isValid()) {
        std::printf("[WARN] Could not create temp directory, skipping file-based tests\n");
    } else {
        // Test 3: PNG export
        std::printf("--- Test 3: exportToFile PNG ---\n");
        QString pngPath = tempDir.filePath("test_export.png");
        QUrl pngUrl = QUrl::fromLocalFile(pngPath);
        QSignalSpy spyPng(&manager, &ExportManager::exportFinished);
        manager.exportToFile(benzeneMolfile, pngUrl, "png");
        bool finished = waitForExport(manager, spyPng, 10000);
        CHECK(finished, "exportToFile PNG finished");
        CHECK(spyPng.count() == 1, "exportFinished signal emitted for PNG");
        QFile pngFile(pngPath);
        CHECK(pngFile.exists(), "PNG file was created");
        CHECK(pngFile.size() > 0, "PNG file is non-empty");
        QList<QVariant> args = spyPng.takeFirst();
        CHECK(args.size() == 2, "exportFinished has 2 arguments");
        CHECK(args[0].toBool() == true, "exportFinished success is true");
        CHECK(args[1].toString().isEmpty(), "exportFinished error is empty");
        pngFile.remove();

        // Test 4: SVG export
        std::printf("--- Test 4: exportToFile SVG ---\n");
        QString svgPath = tempDir.filePath("test_export.svg");
        QUrl svgUrl = QUrl::fromLocalFile(svgPath);
        QSignalSpy spySvg(&manager, &ExportManager::exportFinished);
        manager.exportToFile(benzeneMolfile, svgUrl, "svg");
        finished = waitForExport(manager, spySvg, 10000);
        CHECK(finished, "exportToFile SVG finished");
        CHECK(spySvg.count() == 1, "exportFinished signal emitted for SVG");
        QFile svgFile(svgPath);
        CHECK(svgFile.exists(), "SVG file was created");
        CHECK(svgFile.size() > 0, "SVG file is non-empty");
        args = spySvg.takeFirst();
        CHECK(args.size() == 2, "exportFinished has 2 arguments for SVG");
        CHECK(args[0].toBool() == true, "exportFinished success is true for SVG");
        svgFile.remove();

        // Test 5: PDF export
        std::printf("--- Test 5: exportToFile PDF ---\n");
        QString pdfPath = tempDir.filePath("test_export.pdf");
        QUrl pdfUrl = QUrl::fromLocalFile(pdfPath);
        QSignalSpy spyPdf(&manager, &ExportManager::exportFinished);
        manager.exportToFile(benzeneMolfile, pdfUrl, "pdf");
        finished = waitForExport(manager, spyPdf, 10000);
        CHECK(finished, "exportToFile PDF finished");
        CHECK(spyPdf.count() == 1, "exportFinished signal emitted for PDF");
        QFile pdfFile(pdfPath);
        CHECK(pdfFile.exists(), "PDF file was created");
        CHECK(pdfFile.size() > 0, "PDF file is non-empty");
        args = spyPdf.takeFirst();
        CHECK(args.size() == 2, "exportFinished has 2 arguments for PDF");
        CHECK(args[0].toBool() == true, "exportFinished success is true for PDF");
        pdfFile.remove();

        // Test 5b: document style/page settings actually flow through options
        // -- regression test for the old hardcoded-default-style bug (every
        // export used to ignore the real document's bondLength/page size
        // entirely, confirmed by direct code read before this fix).
        std::printf("--- Test 5b: custom drawingStyle/pageSettings via options ---\n");
        {
            QString customPath = tempDir.filePath("test_custom_style.png");
            QUrl customUrl = QUrl::fromLocalFile(customPath);
            QVariantMap drawingStyle;
            drawingStyle["bondLength"] = 1.5;
            drawingStyle["modelToMmScale"] = 10.0; // deliberately far from the 5.08/1.5 default
            QVariantMap opts;
            opts["drawingStyle"] = drawingStyle;
            QSignalSpy spyCustom(&manager, &ExportManager::exportFinished);
            manager.exportToFile(benzeneMolfile, customUrl, "png", opts);
            bool customFinished = waitForExport(manager, spyCustom, 10000);
            CHECK(customFinished, "exportToFile with custom options finished");
            CHECK(QFile::exists(customPath), "custom-style PNG was created");

            // Also export with defaults (no options) for comparison.
            QString defaultPath = tempDir.filePath("test_default_style.png");
            QUrl defaultUrl = QUrl::fromLocalFile(defaultPath);
            QSignalSpy spyDefault(&manager, &ExportManager::exportFinished);
            manager.exportToFile(benzeneMolfile, defaultUrl, "png");
            waitForExport(manager, spyDefault, 10000);

            QImage customImg(customPath);
            QImage defaultImg(defaultPath);
            CHECK(!customImg.isNull() && !defaultImg.isNull(), "both images load");
            // modelToMmScale=10.0 is roughly 3x the ~3.39 default, so the
            // custom-style PNG (same DPI) must come out noticeably larger --
            // confirms the style actually reached RenderScene, not just that
            // the call didn't crash.
            CHECK(customImg.width() > defaultImg.width() * 2,
                  "custom modelToMmScale produces a visibly larger export than the default");
            QFile::remove(customPath);
            QFile::remove(defaultPath);
        }

        // Test 6: unsupported format
        std::printf("--- Test 6: exportToFile with unsupported format ---\n");
        QString bogusPath = tempDir.filePath("test_export.bogus");
        QUrl bogusUrl = QUrl::fromLocalFile(bogusPath);
        QSignalSpy spyBogus(&manager, &ExportManager::exportFinished);
        manager.exportToFile(benzeneMolfile, bogusUrl, "bogus_format");
        finished = waitForExport(manager, spyBogus, 10000);
        CHECK(finished, "exportToFile with bogus format finished");
        CHECK(spyBogus.count() == 1, "exportFinished signal emitted for bogus format");
        args = spyBogus.takeFirst();
        CHECK(args.size() == 2, "exportFinished has 2 arguments for bogus format");
        CHECK(args[0].toBool() == false, "exportFinished success is false for bogus format");
        CHECK(!args[1].toString().isEmpty(), "exportFinished error is non-empty for bogus format");
    }

    std::printf("\nSummary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
