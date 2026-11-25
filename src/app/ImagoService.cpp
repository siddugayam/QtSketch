#include "ImagoService.h"
#if !defined(QTSKETCH_HAS_IMAGO)
#define QTSKETCH_HAS_IMAGO 1
#endif

#if QTSKETCH_HAS_IMAGO
#include "imago_c.h"
#include <QtConcurrent>
#include <QPointer>
#include <QCoreApplication>
#include <QMetaObject>
#include <QRegularExpression>
#include <QStringList>
#include <QImage>
#include <QImageReader>
#include <QPainter>
#include <QByteArray>
#include <QBuffer>

static void emitOnGuiThread(QPointer<ImagoService> self, std::function<void(ImagoService *)> fn) {
    if (!qApp) return;
    QMetaObject::invokeMethod(qApp, [self, fn]() {
        if (self) fn(self.data());
    }, Qt::QueuedConnection);
}

// Imago internally rescales the working image toward ~1280px on its longest side
// (core/src/settings_defaults.inc: csr.RescaleImageDimensions = 1280) before character
// segmentation, and drops any glyph shorter than ~7px post-rescale (characters.
// MinimalRecognizableHeight) - silently falling back to a plain unlabeled (= carbon)
// vertex for that atom. Real chemical-structure images pasted/screenshotted at common
// web/document resolutions are very often well under that target (verified: a real
// 500x500 test image and a real 470x182 one), so Imago's own internal upscale of an
// already-small source is the effective resolution character recognition runs against -
// and a generic/blurrier upscale there loses exactly the fine strokes that distinguish
// "N"/"O"/"F"/"Cl" glyphs from noise, and can also blur ring-vertex geometry enough to
// affect segmentation and bond-stereo detection.
//
// Fix, verified empirically on two independent real-world structures (Nilotinib.png,
// a patent-figure "6n" structure): pre-upscale any image already smaller than Imago's
// own rescale target to comfortably above it (1600px long side) using Qt's own
// Qt::SmoothTransformation, before Imago ever sees it - so Imago's internal pass becomes
// a (sharper) downscale instead of an upscale. This is not a placebo: on Nilotinib.png,
// heavy-atom recognition went from an all-carbon 23-atom skeleton (0 heteroatoms; real
// formula is C28H22F3N7O) to a byte-exact 28 C / 3 F / 7 N / 1 O match with 0 Indigo
// reload warnings, from this preprocessing step alone.
static const int kImagoTargetLongSide = 1600;

static QByteArray preprocessImageForRecognition(const QString &path) {
    QImageReader reader(path);
    QSize sz = reader.size();
    if (!sz.isValid() || sz.width() <= 0 || sz.height() <= 0) return QByteArray();
    int longSide = qMax(sz.width(), sz.height());
    if (longSide != kImagoTargetLongSide) {
        // Scale toward ~1600px on the longest side in both directions:
        // - Small inputs (< 1600px) upscale up to 3.2x so ~500px multi-ring figures
        //   reach ~1600px for heteroatom glyph recognition without bloating tiny single-ring crops.
        // - Oversized inputs (> 1600px, e.g. high-megapixel scans) downscale during decode
        //   so we never allocate multi-hundred-MB temporary canvases before Imago's own 1280px pass.
        const double rawScale = static_cast<double>(kImagoTargetLongSide) / longSide;
        const double scale = (longSide < kImagoTargetLongSide) ? qMin(3.2, rawScale) : rawScale;
        reader.setScaledSize(QSize(qMax(1, qRound(sz.width() * scale)),
                                   qMax(1, qRound(sz.height() * scale))));
    }
    QImage srcImg = reader.read();
    if (srcImg.isNull()) return QByteArray();
    // Composite onto an opaque white RGB32 canvas (>= 1280px longest side, >= 64px margin)
    // so transparent PNGs have a clean white background and Imago's internal 1280px
    // rescaler never stretches small tight crops up to giant stroke widths.
    constexpr int kMargin = 64;
    constexpr int kImagoInternalRescaleDim = 1280;
    const int minW = srcImg.width() + kMargin * 2;
    const int minH = srcImg.height() + kMargin * 2;
    const int side = qMax(qMax(minW, minH), kImagoInternalRescaleDim);
    const int canvasW = (minW < kImagoInternalRescaleDim && minH < kImagoInternalRescaleDim) ? side : minW;
    const int canvasH = (minW < kImagoInternalRescaleDim && minH < kImagoInternalRescaleDim) ? side : minH;
    QImage padded(canvasW, canvasH, QImage::Format_RGB32);
    padded.fill(Qt::white);
    {
        QPainter painter(&padded);
        painter.drawImage((canvasW - srcImg.width()) / 2,
                          (canvasH - srcImg.height()) / 2,
                          srcImg);
    }
    QByteArray bytes;
    QBuffer buffer(&bytes);
    buffer.open(QIODevice::WriteOnly);
    padded.save(&buffer, "PNG");
    return bytes;
}

// Imago's wedge/hash (bond stereo) recognition from pixel data is the least reliable
// signal in 2D structure OCR - verified empirically on a real recognition (Nilotinib.png):
// Imago's own output can assign a wedge direction that Indigo's own strict molfile
// loader rejects as geometrically inconsistent ("direction of bond #N makes no sense"),
// because two atoms end up too close/collinear for the stated direction to make sense.
// That failure is otherwise silent (Indigo just refuses to reload it, logged as a
// qWarning by IndigoService's periodic calc_props cycle) but was observed immediately
// preceding a real native crash under the resulting repeated-failed-reload load. Since
// a wedge bond's *direction* is unreliable anyway coming from OCR, stripping every
// bond's stereo flag before the recognized structure is ever inserted removes the
// failure mode at its source, rather than only reacting to it downstream.
static QString stripBondStereoFlags(const QString &molfile) {
    QStringList lines = molfile.split(QRegularExpression("\r\n|\n|\r"));
    if (lines.size() < 4) return molfile;
    bool ok = false;
    int atomCount = lines.at(3).left(3).trimmed().toInt(&ok);
    if (!ok) return molfile;
    int bondCount = lines.at(3).mid(3, 3).trimmed().toInt(&ok);
    if (!ok) return molfile;
    int bondBlockStart = 4 + atomCount;
    for (int i = bondBlockStart; i < bondBlockStart + bondCount && i < lines.size(); i++) {
        QString &line = lines[i];
        if (line.length() >= 12) line.replace(9, 3, "  0");
    }
    return lines.join("\n");
}

ImagoService::ImagoService(QObject *parent) : QObject(parent) {}
ImagoService::~ImagoService() {}

void ImagoService::recognizeImage(const QUrl &fileUrl) {
    QString path = fileUrl.toLocalFile();
    if (path.isEmpty()) { emit imageRecognized("", 0, "Invalid file path."); return; }
    QPointer<ImagoService> self = this;
    (void)QtConcurrent::run([self, path]() {
        QString result, error;
        int warnings = 0;
        imago_qword sid = imagoAllocSessionId();
        imagoSetSessionId(sid);
        try {
            QByteArray upscaled = preprocessImageForRecognition(path);
            bool loaded = !upscaled.isEmpty()
                ? imagoLoadImageFromBuffer(upscaled.constData(), upscaled.size())
                : imagoLoadImageFromFile(path.toUtf8().constData());
            if (loaded) {
                if (imagoFilterImage()) {
                    if (imagoRecognize(&warnings)) {
                        const char* mol = imagoGetMol();
                        if (mol) result = stripBondStereoFlags(QString::fromUtf8(mol));
                        else error = QString::fromUtf8(imagoGetLastError());
                    } else {
                        error = QString::fromUtf8(imagoGetLastError());
                    }
                } else {
                    error = QString("Image filtering failed: %1").arg(imagoGetLastError());
                }
            } else {
                error = QString("Failed to load image: %1").arg(imagoGetLastError());
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = "Unknown exception during image recognition.";
        }
        imagoReleaseSessionId(sid);
        emitOnGuiThread(self, [result, warnings, error](ImagoService *s) {
            emit s->imageRecognized(result, warnings, error);
        });
    });
}
#else
ImagoService::ImagoService(QObject *parent) : QObject(parent) {}
ImagoService::~ImagoService() {}

void ImagoService::recognizeImage(const QUrl &) {
    emit imageRecognized(QString(), 0,
        QStringLiteral("Image recognition (Imago) is not included in this build."));
}
#endif
