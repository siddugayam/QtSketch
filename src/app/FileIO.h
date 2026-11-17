#ifndef FILEIO_H
#define FILEIO_H

#include <QObject>
#include <QFile>
#include <QSaveFile>
#include <QTextStream>
#include <QUrl>
#include <QDebug>
#include <QFileInfo>
#include <QDir>
#include <QStandardPaths>
#include <QtQml/qqml.h>

class FileIO : public QObject {
    Q_OBJECT
    QML_ELEMENT
public:
    Q_INVOKABLE bool write(const QString& fileUrl, const QString& data) {
        QUrl url(fileUrl);
        QSaveFile file(url.isLocalFile() ? url.toLocalFile() : fileUrl);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            qWarning() << "FileIO write failed:" << file.errorString();
            return false;
        }
        QTextStream out(&file);
        out << data;
        out.flush();
        bool ok = (out.status() == QTextStream::Ok) && (file.error() == QFile::NoError);
        if (ok) ok = file.commit();
        if (!ok) {
            qWarning() << "FileIO write failed:" << file.errorString();
        }
        return ok;
    }

    Q_INVOKABLE QString read(const QString& fileUrl) {
        QUrl url(fileUrl);
        QFile file(url.isLocalFile() ? url.toLocalFile() : fileUrl);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qWarning() << "FileIO read failed:" << file.errorString();
            return "";
        }
        QTextStream in(&file);
        return in.readAll();
    }

    Q_INVOKABLE bool exists(const QString& fileUrl) {
        QUrl url(fileUrl);
        return QFile::exists(url.isLocalFile() ? url.toLocalFile() : fileUrl);
    }

    Q_INVOKABLE bool remove(const QString& fileUrl) {
        QUrl url(fileUrl);
        QString path = url.isLocalFile() ? url.toLocalFile() : fileUrl;

        // Every existing call site (MainWindow.qml/MessageDialogs.qml) only ever
        // deletes recovery/autosave files under this app's own AppDataLocation
        // tree. Guard against a compromised/buggy QML call passing an arbitrary
        // path by refusing to remove anything outside that tree.
        QString appData = QDir::cleanPath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation));
        QString cleanPath = QDir::cleanPath(path);
        if (appData.isEmpty() || !cleanPath.startsWith(appData + "/")) {
            qWarning() << "FileIO remove refused (outside app data directory):" << cleanPath;
            return false;
        }
        return QFile::remove(cleanPath);
    }

    Q_INVOKABLE bool mkpath(const QString& dirUrl) {
        QUrl url(dirUrl);
        QString path = url.isLocalFile() ? url.toLocalFile() : dirUrl;
        QDir dir;
        return dir.mkpath(path);
    }

    Q_INVOKABLE QStringList listFiles(const QString& dirUrl, const QString& filter) {
        QUrl url(dirUrl);
        QString path = url.isLocalFile() ? url.toLocalFile() : dirUrl;
        QDir dir(path);
        return dir.entryList(QStringList() << filter, QDir::Files);
    }

    Q_INVOKABLE QString readImageAsDataUri(const QUrl &fileUrl) {
        QFile file(fileUrl.toLocalFile());
        if (file.size() > 5 * 1024 * 1024) {
            return "";
        }
        if (!file.open(QIODevice::ReadOnly)) {
            return "";
        }
        QString suffix = QFileInfo(fileUrl.toLocalFile()).suffix().toLower();
        QString mimeType;
        if (suffix == "png") mimeType = "image/png";
        else if (suffix == "jpg" || suffix == "jpeg") mimeType = "image/jpeg";
        else if (suffix == "gif") mimeType = "image/gif";
        else if (suffix == "bmp") mimeType = "image/bmp";
        else return "";

        return "data:" + mimeType + ";base64," + QString(file.readAll().toBase64());
    }
};

#endif // FILEIO_H
