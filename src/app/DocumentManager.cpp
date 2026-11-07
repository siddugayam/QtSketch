#include "DocumentManager.h"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QTimer>
#include <QUuid>

namespace {
// Same "walk up from the binary, up to 5 levels" discovery pattern already proven in
// V8Process's own constructor for src/v8_worker.js (v8_process.cpp:23-36) -- templates/
// sits at the repo root, structurally parallel to src/, so the identical loop applies
// unchanged. Returns the templates/ directory path, or an empty string if not found
// within 5 levels (harmless either way -- see the comment at its call site below).
QString findTemplatesDir() {
    QString appDir = QCoreApplication::applicationDirPath();
    QDir dir(appDir);
    for (int i = 0; i < 5; ++i) {
        QString candidate = dir.filePath(QStringLiteral("templates/fg.sdf"));
        if (QFile::exists(candidate)) return dir.filePath(QStringLiteral("templates"));
        if (!dir.cdUp()) break;
    }
    return QString();
}
}

DocumentManager::DocumentManager(QObject *parent) : QObject(parent) {
    addDocument();
}

QVariantList DocumentManager::docIds() const {
    QVariantList out;
    for (int id : m_order) out.append(id);
    return out;
}

void DocumentManager::setActiveDocId(int docId) {
    if (m_activeDocId == docId || !m_documents.contains(docId)) return;
    m_activeDocId = docId;
    emit activeDocIdChanged();
}

int DocumentManager::addDocument() {
    int docId = m_nextDocId++;
    if (!m_templateLibrary) {
        // TemplateLibrary::loadSdf gracefully leaves a hash empty for a missing/
        // unreadable path (confirmed: TemplateLibrary.cpp:27-33, "missing/unreadable
        // file: leave target empty") -- constructing with a not-found directory (empty
        // findTemplatesDir() result, so these paths won't exist either) is safe and
        // never crashes. No existence check needed here; TemplateLibrary already does
        // the only check that matters, per file, internally.
        QString dir = findTemplatesDir();
        m_templateLibrary = std::make_unique<TemplateLibrary>(
            dir + QStringLiteral("/fg.sdf"),
            dir + QStringLiteral("/library.sdf"),
            dir + QStringLiteral("/salts-and-solvents.sdf"));
    }
    auto *proc = new V8Process(this, m_templateLibrary.get());
    QTimer::singleShot(0, proc, [proc]() { proc->init(); });
    m_documents.insert(docId, proc);
    m_recoveryIds.insert(docId, QUuid::createUuid().toString(QUuid::WithoutBraces));
    m_order.append(docId);
    emit docIdsChanged();
    setActiveDocId(docId);
    return docId;
}

void DocumentManager::closeDocument(int docId) {
    if (!m_documents.contains(docId)) return;
    V8Process *proc = m_documents.take(docId);
    m_recoveryIds.remove(docId);
    m_order.removeAll(docId);
    proc->deleteLater();
    emit docIdsChanged();
    if (m_activeDocId == docId) {
        m_activeDocId = -1;
        if (!m_order.isEmpty()) setActiveDocId(m_order.last());
        else emit activeDocIdChanged();
    }
}

QObject* DocumentManager::documentFor(int docId) const {
    return m_documents.value(docId, nullptr);
}

QString DocumentManager::recoveryIdFor(int docId) const {
    return m_recoveryIds.value(docId, QString::number(docId));
}
