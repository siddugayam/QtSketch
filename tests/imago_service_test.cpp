#include <cstdio>
#include <QGuiApplication>
#include <QEventLoop>
#include <QImage>
#include <QPainter>
#include <QPen>
#include <QPolygonF>
#include <QTemporaryDir>
#include <QTimer>
#include <QUrl>
#include <cmath>
#include "app/ImagoService.h"
#include "indigo.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

int main(int argc, char **argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QGuiApplication app(argc, argv);

    QTemporaryDir tmpDir;
    CHECK(tmpDir.isValid(), "temporary directory valid");

    // Draw a crisp black-on-white cyclohexane hexagon image for OCR verification.
    const QString pngPath = tmpDir.filePath("cyclohexane.png");
    {
        QImage img(600, 600, QImage::Format_RGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.setRenderHint(QPainter::Antialiasing, true);
        p.setPen(QPen(Qt::black, 8.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        QPolygonF hex;
        constexpr double kPi = 3.14159265358979323846;
        const double cx = 300.0, cy = 300.0, r = 140.0;
        for (int i = 0; i < 6; ++i) {
            const double a = (i * 60.0 - 90.0) * kPi / 180.0;
            hex << QPointF(cx + r * std::cos(a), cy + r * std::sin(a));
        }
        p.drawPolygon(hex);
        p.end();
        CHECK(img.save(pngPath, "PNG"), "saved synthetic cyclohexane PNG");
    }

    ImagoService svc;
    bool finished = false;
    QString gotMol, gotErr;
    QObject::connect(&svc, &ImagoService::imageRecognized, &svc,
                     [&](const QString &molfile, int, const QString &error) {
                         gotMol = molfile;
                         gotErr = error;
                         finished = true;
                     });

    svc.recognizeImage(QUrl::fromLocalFile(pngPath));

    QEventLoop loop;
    QTimer timer;
    timer.setSingleShot(true);
    QObject::connect(&timer, &QTimer::timeout, &loop, &QEventLoop::quit);
    timer.start(15000);
    while (!finished && timer.isActive()) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 25);
    }

    CHECK(finished, "ImagoService::recognizeImage emitted imageRecognized");
#if QTSKETCH_HAS_IMAGO
    CHECK(gotErr.isEmpty() && !gotMol.isEmpty(), "Imago OCR produced non-empty molfile without error");
    unsigned long long sid = indigoAllocSessionId();
    indigoSetSessionId(sid);
    int mol = indigoLoadMoleculeFromString(gotMol.toUtf8().constData());
    CHECK(mol >= 0, "Indigo loaded Imago-recognized molfile");
    if (mol >= 0) {
        CHECK(indigoCountAtoms(mol) == 6 && indigoCountBonds(mol) == 6,
              "Imago recognized 6-membered ring (6 atoms, 6 bonds)");
        indigoFree(mol);
    }
    indigoReleaseSessionId(sid);
#else
    CHECK(gotMol.isEmpty() && !gotErr.isEmpty(), "Stub ImagoService reports clean unavailable message");
#endif

    std::printf("\nSummary: %d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
