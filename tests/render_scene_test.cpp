// tests/render_scene_test.cpp
// Tests for RenderScene class
#include <cstdio>
#include <cmath>
#include <QGuiApplication>
#include <QImage>
#include <QPainter>
#include <QFontDatabase>
#include "app/render/RenderScene.h"
#include "app/render/Exporter.h"
#include "app/molecule/RenderPrimitives.h"
#include "app/molecule/DrawingStyle.h"
#include "app/molecule/PageSettings.h"
#include <QTemporaryDir>
#include <QSvgRenderer>

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

int main(int argc, char** argv) {
    // RenderScene constructs real QFont/QFontMetricsF and paints via QPainter --
    // this needs the GUI platform plugin initialized (QGuiApplication), not just
    // QCoreApplication or no Application object at all (confirmed real crash,
    // 0xC0000602 fail-fast, before this fix).
    QGuiApplication app(argc, argv);
    std::printf("--- RenderScene Tests ---\n\n");

    // Test 1: build with minimal RenderPrimitives returns true
    std::printf("--- Test 1: build with minimal primitives ---\n");
    RenderScene scene;
    
    // Create minimal RenderPrimitives with 2 atoms and 1 bond
    RenderPrimitives prims;
    
    AtomPrim atom1;
    atom1.id = 0;
    atom1.x = 0.0;
    atom1.y = 0.0;
    atom1.label = "C";
    atom1.element = "C";
    prims.atoms.append(atom1);
    
    AtomPrim atom2;
    atom2.id = 1;
    atom2.x = 1.5;
    atom2.y = 0.0;
    atom2.label = "C";
    atom2.element = "C";
    prims.atoms.append(atom2);
    
    BondPrim bond1;
    bond1.id = 0;
    bond1.begin = 0;
    bond1.end = 1;
    bond1.type = 1;
    prims.bonds.append(bond1);
    
    DrawingStyle style;
    style.modelToMmScale = 5.08 / 1.5;
    
    PageSettings page;
    
    bool built = scene.build(prims, style, page);
    CHECK(built, "build returns true with valid primitives");
    
    // Check model bounds. Both atoms sit on y=0 (a flat, collinear 2-atom
    // molecule), so the bounding box has zero height -- a legitimate shape,
    // not an error -- which makes QRectF::isEmpty() (width<=0 || height<=0)
    // the wrong check here. Width is the meaningful invariant for this fixture.
    const QRectF& bounds = scene.modelBounds();
    CHECK(bounds.width() > 0.0, "modelBounds has non-zero width");
    CHECK(bounds.left() <= 0.0, "modelBounds left is correct");
    CHECK(bounds.top() <= 0.0, "modelBounds top is correct");
    CHECK(bounds.right() >= 1.5, "modelBounds right is correct");
    CHECK(bounds.bottom() >= 0.0, "modelBounds bottom is correct");

    // Test 2: mmBounds scales modelBounds by modelToMmScale
    std::printf("--- Test 2: mmBounds scaling ---\n");
    QRectF mmBounds = scene.mmBounds();
    double scale = scene.modelToMmScale();
    
    // Check that mmBounds dimensions are scaled versions of modelBounds
    double expectedWidthMm = bounds.width() * scale;
    double expectedHeightMm = bounds.height() * scale;
    
    // Allow for floating point comparison with tolerance
    double tolerance = 0.001 * std::max(expectedWidthMm, expectedHeightMm);
    CHECK(std::abs(mmBounds.width() - expectedWidthMm) < tolerance, 
          "mmBounds width matches modelBounds width * scale");
    CHECK(std::abs(mmBounds.height() - expectedHeightMm) < tolerance, 
          "mmBounds height matches modelBounds height * scale");

    // Test 3: pxBounds scales mmBounds by DPI
    std::printf("--- Test 3: pxBounds scaling ---\n");
    double dpi = 72.0;
    QRectF pxBounds = scene.pxBounds(dpi);
    
    const double kMmPerInch = 25.4;
    const double mmToPx = dpi / kMmPerInch;
    
    expectedWidthMm = mmBounds.width() * mmToPx;
    expectedHeightMm = mmBounds.height() * mmToPx;
    
    tolerance = 0.001 * std::max(expectedWidthMm, expectedHeightMm);
    CHECK(std::abs(pxBounds.width() - expectedWidthMm) < tolerance, 
          "pxBounds width matches mmBounds width * mmToPx");
    CHECK(std::abs(pxBounds.height() - expectedHeightMm) < tolerance, 
          "pxBounds height matches mmBounds height * mmToPx");

    // Test 4: Font substitution detection
    std::printf("--- Test 4: Font substitution ---\n");
    RenderScene scene2;
    RenderPrimitives prims2;
    
    AtomPrim atom3;
    atom3.id = 0;
    atom3.x = 0.0;
    atom3.y = 0.0;
    atom3.label = "C";
    atom3.element = "C";
    prims2.atoms.append(atom3);
    
    DrawingStyle style2;
    style2.fontFamily = "ThisFontDoesNotExist12345";
    style2.fontSize = 10;
    style2.modelToMmScale = 5.08 / 1.5;
    
    scene2.build(prims2, style2, page);
    CHECK(scene2.fontSubstituted(), "fontSubstituted is true for nonexistent font");
    CHECK(scene2.requestedFontFamily() == "ThisFontDoesNotExist12345", 
          "requestedFontFamily matches what we set");
    CHECK(!scene2.resolvedFontFamily().isEmpty(), "resolvedFontFamily is not empty");
    CHECK(scene2.requestedFontFamily() != scene2.resolvedFontFamily(), 
          "requested and resolved font families differ");
    CHECK(!scene2.fontWarnings().isEmpty(), "fontWarnings is non-empty");

    // Test 5: No font substitution for real font
    std::printf("--- Test 5: No font substitution for real font ---\n");
    RenderScene scene3;
    RenderPrimitives prims3;
    
    AtomPrim atom4;
    atom4.id = 0;
    atom4.x = 0.0;
    atom4.y = 0.0;
    atom4.label = "C";
    atom4.element = "C";
    prims3.atoms.append(atom4);
    
    DrawingStyle style3;
    style3.fontFamily = "Arial";
    style3.fontSize = 10;
    style3.modelToMmScale = 5.08 / 1.5;
    
    scene3.build(prims3, style3, page);
    CHECK(!scene3.fontSubstituted(), "fontSubstituted is false for real font");

    // Test 6: paint on QImage does not crash
    std::printf("--- Test 6: paint on QImage ---\n");
    QImage image(100, 100, QImage::Format_ARGB32);
    image.fill(Qt::white);
    
    QPainter painter(&image);
    bool painterActive = painter.isActive();
    CHECK(painterActive, "QPainter is active on QImage");
    
    scene.paint(&painter);
    painter.end();
    
    CHECK(!image.isNull(), "QImage is still valid after painting");

    // Test 7: clear resets bounds and font tracking
    std::printf("--- Test 7: clear resets state ---\n");
    RenderScene scene4;
    scene4.build(prims, style, page);
    CHECK(scene4.modelBounds().width() > 0.0, "modelBounds has non-zero width before clear");
    
    scene4.clear();
    CHECK(scene4.modelBounds().isEmpty(), "modelBounds is empty after clear");
    CHECK(scene4.requestedFontFamily().isEmpty(), "requestedFontFamily is empty after clear");
    CHECK(scene4.resolvedFontFamily().isEmpty(), "resolvedFontFamily is empty after clear");
    CHECK(!scene4.fontSubstituted(), "fontSubstituted is false after clear");
    CHECK(scene4.fontWarnings().isEmpty(), "fontWarnings is empty after clear");

    // Test 8: build with more complex primitives
    std::printf("--- Test 8: build with complex primitives ---\n");
    RenderScene scene5;
    RenderPrimitives prims5;
    
    // Add 3 atoms in a triangle
    AtomPrim atom5;
    atom5.id = 0;
    atom5.x = 0.0;
    atom5.y = 0.0;
    atom5.label = "C";
    atom5.element = "C";
    prims5.atoms.append(atom5);
    
    AtomPrim atom6;
    atom6.id = 1;
    atom6.x = 1.5;
    atom6.y = 0.0;
    atom6.label = "C";
    atom6.element = "C";
    prims5.atoms.append(atom6);
    
    AtomPrim atom7;
    atom7.id = 2;
    atom7.x = 0.75;
    atom7.y = 1.3;
    atom7.label = "C";
    atom7.element = "C";
    prims5.atoms.append(atom7);
    
    // Add 3 bonds (triangle)
    BondPrim bond2;
    bond2.id = 0;
    bond2.begin = 0;
    bond2.end = 1;
    bond2.type = 1;
    prims5.bonds.append(bond2);
    
    BondPrim bond3;
    bond3.id = 1;
    bond3.begin = 1;
    bond3.end = 2;
    bond3.type = 1;
    prims5.bonds.append(bond3);
    
    BondPrim bond4;
    bond4.id = 2;
    bond4.begin = 2;
    bond4.end = 0;
    bond4.type = 1;
    prims5.bonds.append(bond4);
    
    built = scene5.build(prims5, style, page);
    CHECK(built, "build returns true with complex primitives");
    
    const QRectF& bounds5 = scene5.modelBounds();
    CHECK(!bounds5.isEmpty(), "modelBounds is non-empty for complex primitives");

    // Test 9: PngExporter actually applies modelToMmScale when painting,
    // not just when allocating the image. Regression test for a real bug:
    // the exporter's painter transform used to be a bare Y-flip with no
    // magnification at all, so painted content came out ~1/modelToMmScale
    // the size the allocated canvas expected (here, modelToMmScale is
    // 5.08/1.5 =~ 3.39x) -- everything exported miniature and offset
    // instead of filling the frame.
    std::printf("--- Test 9: PngExporter applies modelToMmScale to painted content ---\n");
    {
        QTemporaryDir dir;
        QString pngPath = dir.filePath("scale_check.png");
        PngExporter exporter;
        exporter.setResolution(300);
        bool exported = exporter.exportToFile(scene5, pngPath);
        CHECK(exported, "PngExporter::exportToFile succeeds");

        QImage img(pngPath);
        CHECK(!img.isNull(), "exported PNG loads back as a valid image");

        int minX = img.width(), maxX = -1, minY = img.height(), maxY = -1;
        for (int y = 0; y < img.height(); ++y) {
            for (int x = 0; x < img.width(); ++x) {
                if (img.pixelColor(x, y) != Qt::white) {
                    if (x < minX) minX = x;
                    if (x > maxX) maxX = x;
                    if (y < minY) minY = y;
                    if (y > maxY) maxY = y;
                }
            }
        }
        bool foundPaintedContent = maxX >= minX && maxY >= minY;
        CHECK(foundPaintedContent, "image contains non-background (painted) pixels");

        if (foundPaintedContent) {
            int paintedWidthPx = maxX - minX + 1;
            // Expected: modelBounds width (1.5 units) * modelToMmScale (~3.39)
            // * dpi/25.4 (300/25.4 =~ 11.8) =~ 60px. A bug that skips
            // modelToMmScale entirely would produce ~1.5 * 11.8 =~ 18px --
            // well below this floor.
            CHECK(paintedWidthPx > 40, "painted content spans the full scaled width, not ~1/3.39 of it");
        }
    }

    // Test 10: export stays inside the frame wherever the scene sits.
    std::printf("--- Test 10: PNG/SVG export frames a scene that is away from the model origin ---\n");
    {
        const double offsets[][2] = {{0.0, 0.0}, {0.0, -0.75}, {10.0, 7.0}, {-7.0, -9.0}};
        for (const auto& off : offsets) {
            RenderPrimitives p;
            const double pts[3][2] = {{off[0] + 3.0, off[1]}, {off[0], off[1]}, {off[0], off[1] + 1.5}};
            for (int i = 0; i < 3; ++i) {
                AtomPrim a; a.id = i; a.x = pts[i][0]; a.y = pts[i][1]; a.element = "C";
                p.atoms.append(a);
            }
            for (int i = 0; i < 2; ++i) {
                BondPrim b; b.id = i; b.begin = i; b.end = i + 1; b.type = 1;
                p.bonds.append(b);
            }
            RenderScene s;
            s.build(p, style, page);

            auto inkBox = [](const QImage& img, int& minX, int& maxX, int& minY, int& maxY) {
                minX = img.width(); maxX = -1; minY = img.height(); maxY = -1;
                for (int y = 0; y < img.height(); ++y)
                    for (int x = 0; x < img.width(); ++x) {
                        if (img.pixelColor(x, y) != Qt::white) {
                            minX = qMin(minX, x); maxX = qMax(maxX, x);
                            minY = qMin(minY, y); maxY = qMax(maxY, y);
                        }
                    }
            };

            QTemporaryDir dir;
            const QString pngPath = dir.filePath("offset.png");
            PngExporter png;
            png.setResolution(300);
            png.exportToFile(s, pngPath);
            QImage img(pngPath);
            int minX, maxX, minY, maxY;
            inkBox(img, minX, maxX, minY, maxY);
            // 3.0 x 1.5 model units at 300 dpi is about 120 x 60 px.
            CHECK(maxX - minX > img.width() * 0.9 && maxY - minY > img.height() * 0.9,
                  "PNG: both bonds are drawn across the full frame");

            const QString svgPath = dir.filePath("offset.svg");
            SvgExporter svg;
            svg.exportToFile(s, svgPath);
            QSvgRenderer renderer(svgPath);
            QImage svgImg(240, 120, QImage::Format_ARGB32);
            svgImg.fill(Qt::white);
            { QPainter sp(&svgImg); renderer.render(&sp); }
            inkBox(svgImg, minX, maxX, minY, maxY);
            CHECK(renderer.isValid() && maxX - minX > svgImg.width() * 0.9 && maxY - minY > svgImg.height() * 0.9,
                  "SVG: both bonds are visible across the full viewBox");
        }
    }

    // Test 11: atom labels reach the exported file at a readable size.
    std::printf("--- Test 11: PNG/SVG export draws atom labels ---\n");
    {
        RenderPrimitives p;
        AtomPrim o; o.id = 0; o.x = 4.0; o.y = 2.0; o.element = "O"; o.label = "O";
        p.atoms.append(o);
        RenderScene s;
        s.build(p, style, page);

        auto inkBox = [](const QImage& img, int& minX, int& maxX, int& minY, int& maxY) {
            minX = img.width(); maxX = -1; minY = img.height(); maxY = -1;
            for (int y = 0; y < img.height(); ++y)
                for (int x = 0; x < img.width(); ++x)
                    if (qGray(img.pixel(x, y)) < 128) {
                        minX = qMin(minX, x); maxX = qMax(maxX, x);
                        minY = qMin(minY, y); maxY = qMax(maxY, y);
                    }
        };

        QTemporaryDir dir;
        const QString pngPath = dir.filePath("label.png");
        PngExporter png;
        png.setResolution(300);
        CHECK(png.exportToFile(s, pngPath), "PNG export succeeds");
        QImage img(pngPath);
        int minX, maxX, minY, maxY;
        inkBox(img, minX, maxX, minY, maxY);
        CHECK(maxX >= minX, "PNG: the O label is drawn");
        // 10 pt at 300 dpi is about 42 px per em; a capital O is roughly 30 px tall.
        const int h = maxY - minY + 1;
        CHECK(h > 15 && h < 60, "PNG: the label is drawn at its point size");
        CHECK(minX > 0 && minY > 0 && maxX < img.width() - 1 && maxY < img.height() - 1,
              "PNG: the label is not clipped by the frame");

        const QString svgPath = dir.filePath("label.svg");
        SvgExporter svg;
        CHECK(svg.exportToFile(s, svgPath), "SVG export succeeds");
        QSvgRenderer renderer(svgPath);
        QImage svgImg(renderer.defaultSize() * 4, QImage::Format_ARGB32);
        svgImg.fill(Qt::white);
        { QPainter sp(&svgImg); renderer.render(&sp); }
        inkBox(svgImg, minX, maxX, minY, maxY);
        CHECK(renderer.isValid() && maxX >= minX, "SVG: the O label is drawn");
    }

    std::printf("\nSummary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
