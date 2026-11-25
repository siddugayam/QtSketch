#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QQmlExpression>
#include <QQuickStyle>
#include <QQuickWindow>
#include <QStyleHints>
#include <QFile>
#include <QImage>
#include <QTextStream>
#include <QTimer>
#include <QStandardPaths>
#include <QDir>
#include <QFileInfo>
#include <QIcon>
#ifdef Q_OS_MACOS
#include <mach-o/dyld.h>
#endif

#include <QUrl>
#include "IndigoService.h"
#include "ImagoService.h"
#include "FileIO.h"
#include "AppController.h"

static QFile* g_logFile = nullptr;
static QTextStream* g_logStream = nullptr;

namespace {
class QaGrabber final : public QObject {
    Q_OBJECT
public:
    explicit QaGrabber(QObject *parent = nullptr) : QObject(parent) {}

    Q_INVOKABLE bool grabWindow(QQuickWindow *win, const QString &filePath) const {
        if (!win || filePath.isEmpty()) return false;
        const QImage img = win->grabWindow();
        if (img.isNull()) return false;
        QDir().mkpath(QFileInfo(filePath).absolutePath());
        return img.save(filePath);
    }

    Q_INVOKABLE void quitApp(int exitCode = 0) const {
        QCoreApplication::exit(exitCode);
    }
};
} // namespace

void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    if (g_logStream) {
        *g_logStream << msg << "\n";
        g_logStream->flush();
    }
    // Also print to stderr so we can see it in terminal!
    fprintf(stderr, "%s\n", msg.toLocal8Bit().constData());
    fflush(stderr);
}

int main(int argc, char *argv[])
{
    QString macBundleQmlDir;
#ifdef Q_OS_MACOS
    // Resolve bundled PlugIns and Resources/qml directly from the Mach-O executable
    // path so standalone .app launches always find libqcocoa.dylib and QML modules,
    // including when opened directly inside a mounted .dmg or under App Translocation.
    {
        uint32_t bufSize = 0;
        _NSGetExecutablePath(nullptr, &bufSize);
        QByteArray rawPath(static_cast<int>(bufSize), '\0');
        if (_NSGetExecutablePath(rawPath.data(), &bufSize) == 0) {
            const QFileInfo exeInfo(QString::fromUtf8(rawPath.constData()));
            const QString resolvedExe = exeInfo.canonicalFilePath().isEmpty()
                                            ? exeInfo.absoluteFilePath()
                                            : exeInfo.canonicalFilePath();
            const QDir macosDir = QFileInfo(resolvedExe).dir();
            const QString plugInsDir = QDir::cleanPath(macosDir.filePath(QStringLiteral("../PlugIns")));
            if (QDir(plugInsDir).exists()) {
                if (qEnvironmentVariableIsEmpty("QT_PLUGIN_PATH")) {
                    qputenv("QT_PLUGIN_PATH", plugInsDir.toUtf8());
                }
                if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM_PLUGIN_PATH")) {
                    const QString platformsDir = plugInsDir + QStringLiteral("/platforms");
                    if (QDir(platformsDir).exists()) {
                        qputenv("QT_QPA_PLATFORM_PLUGIN_PATH", platformsDir.toUtf8());
                    }
                }
                QCoreApplication::setLibraryPaths({plugInsDir});
            }
            const QString qmlDir = QDir::cleanPath(macosDir.filePath(QStringLiteral("../Resources/qml")));
            if (QDir(qmlDir).exists()) {
                macBundleQmlDir = qmlDir;
                if (qEnvironmentVariableIsEmpty("QML2_IMPORT_PATH")) {
                    qputenv("QML2_IMPORT_PATH", qmlDir.toUtf8());
                }
                if (qEnvironmentVariableIsEmpty("QML_IMPORT_PATH")) {
                    qputenv("QML_IMPORT_PATH", qmlDir.toUtf8());
                }
            }
        }
    }
#endif

    // FluentWinUI3 gives every stock control the Windows 11 Fluent metrics
    // (32px control height, 4px radii, Fluent focus rings) instead of the Basic
    // style's touch-sized 100x40 buttons and 140x40 combo boxes, which were the
    // single cause of the oversized caption buttons, page-size box and zoom
    // slider. Must be set before the QGuiApplication is constructed.
    QQuickStyle::setStyle("FluentWinUI3");

    QGuiApplication app(argc, argv);
    // The app ships a single light theme. Without pinning the colour scheme the
    // Fluent style follows the OS setting, so every stock control (dialogs,
    // text fields, combo boxes) rendered dark inside the light app on a machine
    // set to dark mode -- the defect that made the dialogs look foreign.
    app.styleHints()->setColorScheme(Qt::ColorScheme::Light);
    // Required for the QML Settings backing store (tool-panel section states)
    app.setOrganizationName("sketch");
    app.setApplicationName("sketch");

    QIcon appIcon;
    appIcon.addFile(":/qt/qml/Sketch/App/icons/app_icon_16.png", QSize(16, 16));
    appIcon.addFile(":/qt/qml/Sketch/App/icons/app_icon_32.png", QSize(32, 32));
    appIcon.addFile(":/qt/qml/Sketch/App/icons/app_icon_48.png", QSize(48, 48));
    appIcon.addFile(":/qt/qml/Sketch/App/icons/app_icon.png", QSize(256, 256));
    app.setWindowIcon(appIcon);

    // AppDataLocation depends on organizationName/applicationName above, so this
    // must run after they're set. A relative "qml_errors.log" resolved to whatever
    // the launch CWD happened to be, and Truncate wiped a crashing run's log
    // evidence on the very next launch -- use a stable writable path and append.
    QString logDir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    QDir().mkpath(logDir);
    g_logFile = new QFile(logDir + "/qml_errors.log");
    if (g_logFile->open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text)) {
        g_logStream = new QTextStream(g_logFile);
    }

    qInstallMessageHandler(myMessageOutput);

    QQmlApplicationEngine engine;
    if (!macBundleQmlDir.isEmpty()) {
        engine.addImportPath(macBundleQmlDir);
    }

    engine.loadFromModule("Sketch.App", "MainWindow");
    if (engine.rootObjects().isEmpty()) {
        qWarning() << "Error: rootObjects is empty! Failed to load QML.";
        return -1;
    }

    const QString qaScriptPath = qEnvironmentVariable("QTSKETCH_QA_SCRIPT");
    if (!qaScriptPath.isEmpty()) {
        QFile qaFile(qaScriptPath);
        if (qaFile.size() <= 1024 * 1024 && qaFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            auto *grabber = new QaGrabber(&engine);
            engine.rootContext()->setContextProperty(QStringLiteral("qaGrabber"), grabber);
            const QString script = QString::fromUtf8(qaFile.readAll());
            QObject *rootObj = engine.rootObjects().constFirst();
            QTimer::singleShot(600, rootObj, [&engine, rootObj, script, qaScriptPath]() {
                QQmlExpression expr(engine.rootContext(), rootObj, script);
                expr.setSourceLocation(qaScriptPath, 1);
                bool valueIsUndefined = false;
                expr.evaluate(&valueIsUndefined);
                if (expr.hasError()) {
                    qWarning() << "QA script error:" << expr.error().toString();
                    QCoreApplication::exit(2);
                }
            });
        }
    }

    return app.exec();
}

#include "main.moc"


