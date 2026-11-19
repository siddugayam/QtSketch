// src/app/molecule/BiopolymerAsyncOperations.cpp
#include "BiopolymerAsyncOperations.h"
#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include "BiopolymerImportSession.h"
#include "BiopolymerExportSession.h"
#include <QThreadPool>
#include <QtConcurrent/QtConcurrent>
#include <QThread>
#include <QCoreApplication>

// ========== OperationId ==========

OperationId::OperationId() : m_id(0) {
}

OperationId::OperationId(quint64 id) : m_id(id) {
}

quint64 OperationId::value() const {
    return m_id;
}

bool OperationId::isValid() const {
    return m_id > 0;
}

bool OperationId::operator==(const OperationId& other) const {
    return m_id == other.m_id;
}

bool OperationId::operator!=(const OperationId& other) const {
    return m_id != other.m_id;
}

bool OperationId::operator<(const OperationId& other) const {
    return m_id < other.m_id;
}

OperationId::operator quint64() const {
    return m_id;
}

// ========== ProgressInfo ==========

bool ProgressInfo::isActive() const {
    return status == OperationStatus::Running || status == OperationStatus::Pending;
}

bool ProgressInfo::isFinished() const {
    return status == OperationStatus::Completed || 
           status == OperationStatus::Failed || 
           status == OperationStatus::Cancelled;
}

bool ProgressInfo::isSuccessful() const {
    return status == OperationStatus::Completed && errorMessage.isEmpty();
}

QString ProgressInfo::progressString() const {
    if (totalSteps > 0) {
        return QString("%1% (%2/%3)").arg(progress).arg(currentStep).arg(totalSteps);
    }
    return QString("%1%").arg(progress);
}

QString ProgressInfo::elapsedString() const {
    if (elapsedMilliseconds < 1000) {
        return QString("%1 ms").arg(elapsedMilliseconds);
    } else if (elapsedMilliseconds < 60000) {
        return QString("%1 s").arg(elapsedMilliseconds / 1000.0, 0, 'f', 2);
    } else {
        int minutes = elapsedMilliseconds / 60000;
        int seconds = (elapsedMilliseconds % 60000) / 1000;
        return QString("%1:%2 min").arg(minutes).arg(seconds, 2, 'f', 0, '0');
    }
}

// ========== CancellationToken ==========

CancellationToken::CancellationToken() : m_cancelled(0), m_parentCancelled(0) {
}

bool CancellationToken::isCancelled() const {
    return m_cancelled.loadAcquire() != 0 || 
           (m_parentCancelled.loadAcquire() != 0 && *reinterpret_cast<QAtomicInt*>(m_parentCancelled.loadAcquire()) != 0);
}

void CancellationToken::cancel() {
    m_cancelled.storeRelease(1);
}

void CancellationToken::reset() {
    m_cancelled.storeRelease(0);
}

CancellationToken CancellationToken::createChildToken() const {
    CancellationToken child;
    child.m_parentCancelled.storeRelease(reinterpret_cast<quintptr>(const_cast<QAtomicInt*>(&m_cancelled)));
    return child;
}

// ========== OperationContext ==========

OperationContext::OperationContext(OperationId id, AsyncOperationType type,
                                   const BiopolymerWorkspace& workspace, QObject* parent)
    : QObject(parent),
      m_id(id),
      m_type(type),
      m_workspace(workspace) {
    
    m_progress.id = id;
    m_progress.type = type;
    m_progress.status = OperationStatus::Pending;
    m_progress.startTime = QDateTime::currentDateTime();
}

OperationId OperationContext::id() const {
    return m_id;
}

AsyncOperationType OperationContext::type() const {
    return m_type;
}

const BiopolymerWorkspace& OperationContext::workspace() const {
    return m_workspace;
}

CancellationToken& OperationContext::cancellationToken() {
    return m_cancellationToken;
}

const CancellationToken& OperationContext::cancellationToken() const {
    return m_cancellationToken;
}

ProgressInfo& OperationContext::progress() {
    return m_progress;
}

const ProgressInfo& OperationContext::progress() const {
    return m_progress;
}

void OperationContext::setProgress(int progress, const QString& description) {
    m_progress.progress = progress;
    if (!description.isEmpty()) {
        m_progress.currentStepDescription = description;
    }
    emit progressChanged(progress, description);
}

void OperationContext::incrementProgress(const QString& description) {
    m_progress.progress += 1;
    if (m_progress.progress > 100) {
        m_progress.progress = 100;
    }
    if (!description.isEmpty()) {
        m_progress.currentStepDescription = description;
    }
    emit progressChanged(m_progress.progress, description);
}

void OperationContext::setStep(int step, int total, const QString& description) {
    m_progress.currentStep = step;
    m_progress.totalSteps = total;
    if (total > 0) {
        m_progress.progress = static_cast<int>((static_cast<qreal>(step) / total) * 100);
    }
    if (!description.isEmpty()) {
        m_progress.currentStepDescription = description;
    }
    emit stepChanged(step, total, description);
}

void OperationContext::setResult(const QString& data, const QMap<QString, QVariant>& metadata) {
    m_progress.resultData = data;
    m_progress.metadata = metadata;
    m_progress.status = OperationStatus::Completed;
    m_progress.endTime = QDateTime::currentDateTime();
    m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
    emit statusChanged(OperationStatus::Completed);
    emit finished(data, QString());
}

void OperationContext::setError(const QString& error) {
    m_progress.errorMessage = error;
    m_progress.status = OperationStatus::Failed;
    m_progress.endTime = QDateTime::currentDateTime();
    m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
    emit statusChanged(OperationStatus::Failed);
    emit finished(QString(), error);
}

bool OperationContext::isCancelled() const {
    return m_cancellationToken.isCancelled();
}

bool OperationContext::shouldStop() const {
    return isCancelled();
}

// ========== AsyncResult ==========

template <typename T>
AsyncResult<T>::AsyncResult() = default;

template <typename T>
AsyncResult<T>::AsyncResult(const T& value) : m_value(value) {
}

template <typename T>
AsyncResult<T>::AsyncResult(T&& value) : m_value(std::move(value)) {
}

template <typename T>
AsyncResult<T>::AsyncResult(const QString& error) : m_error(error) {
}

template <typename T>
bool AsyncResult<T>::isSuccess() const {
    return m_error.isEmpty() && !m_cancelled;
}

template <typename T>
bool AsyncResult<T>::isFailure() const {
    return !m_error.isEmpty() && !m_cancelled;
}

template <typename T>
bool AsyncResult<T>::isCancelled() const {
    return m_cancelled;
}

template <typename T>
const T& AsyncResult<T>::value() const {
    return m_value;
}

template <typename T>
T& AsyncResult<T>::value() {
    return m_value;
}

template <typename T>
const QString& AsyncResult<T>::error() const {
    return m_error;
}

template <typename T>
void AsyncResult<T>::setValue(const T& value) {
    m_value = value;
    m_error.clear();
    m_cancelled = false;
}

template <typename T>
void AsyncResult<T>::setValue(T&& value) {
    m_value = std::move(value);
    m_error.clear();
    m_cancelled = false;
}

template <typename T>
void AsyncResult<T>::setError(const QString& error) {
    m_error = error;
    m_cancelled = false;
}

template <typename T>
void AsyncResult<T>::setCancelled() {
    m_cancelled = true;
    m_error.clear();
}

template <typename T>
QString AsyncResult<T>::toString() const {
    if (m_cancelled) {
        return "Cancelled";
    } else if (!m_error.isEmpty()) {
        return QString("Error: %1").arg(m_error);
    } else {
        return "Success";
    }
}

// ========== AsyncTask ==========

AsyncTask::AsyncTask(OperationId id, AsyncOperationType type,
                     const BiopolymerWorkspace& workspace, QObject* parent)
    : QObject(parent),
      m_id(id),
      m_type(type),
      m_workspace(workspace) {
    
    m_progress.id = id;
    m_progress.type = type;
    m_progress.status = OperationStatus::Pending;
}

AsyncTask::~AsyncTask() {
}

OperationId AsyncTask::id() const {
    return m_id;
}

AsyncOperationType AsyncTask::type() const {
    return m_type;
}

const BiopolymerWorkspace& AsyncTask::workspace() const {
    return m_workspace;
}

ProgressInfo AsyncTask::progress() const {
    return m_progress;
}

const CancellationToken& AsyncTask::cancellationToken() const {
    return m_cancellationToken;
}

bool AsyncTask::isCancelled() const {
    return m_cancellationToken.isCancelled();
}

void AsyncTask::cancel() {
    m_cancellationToken.cancel();
    m_progress.status = OperationStatus::Cancelled;
    emit statusChanged(m_id, OperationStatus::Cancelled);
}

void AsyncTask::pause() {
    // Default implementation: do nothing
    // Subclasses can override to implement pausing
}

void AsyncTask::resume() {
    // Default implementation: do nothing
    // Subclasses can override to implement resuming
}

void AsyncTask::updateProgress(int progress, const QString& description) {
    m_progress.progress = progress;
    if (!description.isEmpty()) {
        m_progress.currentStepDescription = description;
    }
    emit progress(m_id, progress, description);
}

void AsyncTask::updateStep(int current, int total, const QString& description) {
    m_progress.currentStep = current;
    m_progress.totalSteps = total;
    if (total > 0) {
        m_progress.progress = static_cast<int>((static_cast<qreal>(current) / total) * 100);
    }
    if (!description.isEmpty()) {
        m_progress.currentStepDescription = description;
    }
    emit stepChanged(m_id, current, total, description);
}

void AsyncTask::reportSuccess(const QString& result, const QMap<QString, QVariant>& metadata) {
    m_progress.resultData = result;
    m_progress.metadata = metadata;
    m_progress.status = OperationStatus::Completed;
    m_progress.endTime = QDateTime::currentDateTime();
    m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
    emit statusChanged(m_id, OperationStatus::Completed);
    emit finished(m_id, result, QString());
}

void AsyncTask::reportError(const QString& error) {
    m_progress.errorMessage = error;
    m_progress.status = OperationStatus::Failed;
    m_progress.endTime = QDateTime::currentDateTime();
    m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
    emit statusChanged(m_id, OperationStatus::Failed);
    emit finished(m_id, QString(), error);
}

bool AsyncTask::shouldStop() const {
    return isCancelled();
}

void AsyncTask::sleepWithCancellation(int ms) {
    if (shouldStop()) {
        return;
    }
    
    int remaining = ms;
    const int checkInterval = 50; // Check every 50ms
    
    while (remaining > 0 && !shouldStop()) {
        int sleepTime = qMin(remaining, checkInterval);
        QThread::msleep(sleepTime);
        remaining -= sleepTime;
        
        // Process events to allow cancellation signal to be received
        QCoreApplication::processEvents(QEventLoop::AllEvents, checkInterval);
    }
}

// ========== BiopolymerAsyncOperations ==========

BiopolymerAsyncOperations::BiopolymerAsyncOperations(
    const BiopolymerWorkspace& workspace, QObject* parent)
    : QObject(parent),
      m_workspace(workspace),
      m_threadPool(nullptr),
      m_ownsThreadPool(false),
      m_maxConcurrentOperations(4),
      m_nextId(1) {
    
    ensureThreadPool();
}

BiopolymerAsyncOperations::~BiopolymerAsyncOperations() {
    cancelAllOperations();
    
    if (m_ownsThreadPool && m_threadPool) {
        m_threadPool->waitForDone();
        delete m_threadPool;
        m_threadPool = nullptr;
    }
}

void BiopolymerAsyncOperations::ensureThreadPool() {
    if (!m_threadPool) {
        m_threadPool = new QThreadPool(this);
        m_threadPool->setMaxThreadCount(m_maxConcurrentOperations);
        m_ownsThreadPool = true;
    }
}

OperationId BiopolymerAsyncOperations::generateOperationId() {
    return OperationId(m_nextId.fetchAndAddAcquire(1));
}

AsyncTask* BiopolymerAsyncOperations::getOperation(OperationId id) const {
    QMutexLocker locker(&m_mutex);
    return m_operations.value(id, nullptr);
}

void BiopolymerAsyncOperations::addOperation(AsyncTask* task) {
    QMutexLocker locker(&m_mutex);
    m_operations[task->id()] = task;
}

void BiopolymerAsyncOperations::removeOperation(OperationId id) {
    QMutexLocker locker(&m_mutex);
    m_operations.remove(id);
}

OperationId BiopolymerAsyncOperations::startImport(const QString& format, 
                                                   const QString& seqType, 
                                                   const QString& data) {
    OperationId id = generateOperationId();
    
    ImportTask* task = new ImportTask(id, m_workspace, format, seqType, data);
    
    // Connect signals
    connect(task, &AsyncTask::started, this, [this, id]() {
        emit operationStarted(id, task->type());
    });
    connect(task, &AsyncTask::progress, this, &BiopolymerAsyncOperations::operationProgress);
    connect(task, &AsyncTask::stepChanged, this, &BiopolymerAsyncOperations::operationStepChanged);
    connect(task, &AsyncTask::statusChanged, this, &BiopolymerAsyncOperations::onStatusChanged);
    connect(task, &AsyncTask::finished, this, &BiopolymerAsyncOperations::onTaskFinished);
    
    addOperation(task);
    
    QtConcurrent::run([task]() {
        task->execute();
    });
    
    return id;
}

OperationId BiopolymerAsyncOperations::startExport(const QString& format, 
                                                  const BiopolymerState& state) {
    OperationId id = generateOperationId();
    
    ExportTask* task = new ExportTask(id, m_workspace, format, state);
    
    connect(task, &AsyncTask::started, this, [this, id]() {
        emit operationStarted(id, task->type());
    });
    connect(task, &AsyncTask::progress, this, &BiopolymerAsyncOperations::operationProgress);
    connect(task, &AsyncTask::stepChanged, this, &BiopolymerAsyncOperations::operationStepChanged);
    connect(task, &AsyncTask::statusChanged, this, &BiopolymerAsyncOperations::onStatusChanged);
    connect(task, &AsyncTask::finished, this, &BiopolymerAsyncOperations::onTaskFinished);
    
    addOperation(task);
    
    QtConcurrent::run([task]() {
        task->execute();
    });
    
    return id;
}

OperationId BiopolymerAsyncOperations::startVectorExport(
    const QString& format, const BiopolymerState& state,
    const QMap<QString, QVariant>& options) {
    
    OperationId id = generateOperationId();
    
    VectorExportTask* task = new VectorExportTask(id, m_workspace, format, state, options);
    
    connect(task, &AsyncTask::started, this, [this, id]() {
        emit operationStarted(id, task->type());
    });
    connect(task, &AsyncTask::progress, this, &BiopolymerAsyncOperations::operationProgress);
    connect(task, &AsyncTask::stepChanged, this, &BiopolymerAsyncOperations::operationStepChanged);
    connect(task, &AsyncTask::statusChanged, this, &BiopolymerAsyncOperations::onStatusChanged);
    connect(task, &AsyncTask::finished, this, &BiopolymerAsyncOperations::onTaskFinished);
    
    addOperation(task);
    
    QtConcurrent::run([task]() {
        task->execute();
    });
    
    return id;
}

OperationId BiopolymerAsyncOperations::startValidation(const BiopolymerState& state) {
    OperationId id = generateOperationId();
    
    ValidationTask* task = new ValidationTask(id, m_workspace, state);
    
    connect(task, &AsyncTask::started, this, [this, id]() {
        emit operationStarted(id, task->type());
    });
    connect(task, &AsyncTask::progress, this, &BiopolymerAsyncOperations::operationProgress);
    connect(task, &AsyncTask::stepChanged, this, &BiopolymerAsyncOperations::operationStepChanged);
    connect(task, &AsyncTask::statusChanged, this, &BiopolymerAsyncOperations::onStatusChanged);
    connect(task, &AsyncTask::finished, this, &BiopolymerAsyncOperations::onTaskFinished);
    
    addOperation(task);
    
    QtConcurrent::run([task]() {
        task->execute();
    });
    
    return id;
}

bool BiopolymerAsyncOperations::cancelOperation(OperationId id) {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        task->cancel();
        emit operationCancelled(id);
        return true;
    }
    
    return false;
}

void BiopolymerAsyncOperations::cancelAllOperations() {
    QMutexLocker locker(&m_mutex);
    
    for (AsyncTask* task : m_operations.values()) {
        task->cancel();
    }
    
    // Emit cancelled for all
    for (OperationId id : m_operations.keys()) {
        emit operationCancelled(id);
    }
}

bool BiopolymerAsyncOperations::pauseOperation(OperationId id) {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        task->pause();
        return true;
    }
    
    return false;
}

bool BiopolymerAsyncOperations::resumeOperation(OperationId id) {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        task->resume();
        return true;
    }
    
    return false;
}

OperationStatus BiopolymerAsyncOperations::getOperationStatus(OperationId id) const {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        return task->progress().status;
    }
    
    return OperationStatus::Unknown;
}

ProgressInfo BiopolymerAsyncOperations::getOperationProgress(OperationId id) const {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        return task->progress();
    }
    
    return ProgressInfo();
}

bool BiopolymerAsyncOperations::isOperationRunning(OperationId id) const {
    OperationStatus status = getOperationStatus(id);
    return status == OperationStatus::Running || status == OperationStatus::Pending;
}

bool BiopolymerAsyncOperations::hasRunningOperations() const {
    QMutexLocker locker(&m_mutex);
    
    for (AsyncTask* task : m_operations) {
        if (task && (task->progress().status == OperationStatus::Running ||
                      task->progress().status == OperationStatus::Pending)) {
            return true;
        }
    }
    
    return false;
}

QList<OperationId> BiopolymerAsyncOperations::activeOperationIds() const {
    QMutexLocker locker(&m_mutex);
    QList<OperationId> ids;
    
    for (AsyncTask* task : m_operations) {
        if (task && task->progress().isActive()) {
            ids.append(task->id());
        }
    }
    
    return ids;
}

QList<OperationId> BiopolymerAsyncOperations::allOperationIds() const {
    QMutexLocker locker(&m_mutex);
    return m_operations.keys();
}

void BiopolymerAsyncOperations::cleanupFinishedOperations() {
    QMutexLocker locker(&m_mutex);
    
    QList<OperationId> toRemove;
    for (auto it = m_operations.constBegin(); it != m_operations.constEnd(); ++it) {
        if (it.value()->progress().isFinished()) {
            toRemove.append(it.key());
            delete it.value();
        }
    }
    
    for (OperationId id : toRemove) {
        m_operations.remove(id);
    }
}

void BiopolymerAsyncOperations::clearAllOperations() {
    cancelAllOperations();
    
    QMutexLocker locker(&m_mutex);
    qDeleteAll(m_operations.values());
    m_operations.clear();
}

void BiopolymerAsyncOperations::setMaxConcurrentOperations(int max) {
    m_maxConcurrentOperations = max;
    if (m_threadPool) {
        m_threadPool->setMaxThreadCount(max);
    }
}

int BiopolymerAsyncOperations::maxConcurrentOperations() const {
    return m_maxConcurrentOperations;
}

void BiopolymerAsyncOperations::setThreadPoolSize(int size) {
    setMaxConcurrentOperations(size);
}

int BiopolymerAsyncOperations::threadPoolSize() const {
    return maxConcurrentOperations();
}

void BiopolymerAsyncOperations::onTaskFinished(OperationId id, 
                                               const QString& result, 
                                               const QString& error) {
    QMutexLocker locker(&m_mutex);
    AsyncTask* task = m_operations.value(id, nullptr);
    
    if (task) {
        // Emit finished signal
        emit operationFinished(id, result, error);
        
        // Clean up the task
        task->deleteLater();
        m_operations.remove(id);
        
        // Check if all operations are finished
        if (!hasRunningOperations()) {
            emit allOperationsFinished();
        }
    }
}

void BiopolymerAsyncOperations::onStatusChanged(OperationId id, OperationStatus status) {
    emit operationStatusChanged(id, status);
    
    if (status == OperationStatus::Cancelled) {
        emit operationCancelled(id);
    }
}

// ========== ImportTask ==========

ImportTask::ImportTask(OperationId id, const BiopolymerWorkspace& workspace,
                       const QString& format, const QString& seqType, const QString& data,
                       QObject* parent)
    : AsyncTask(id, AsyncOperationType::ImportSequence, workspace, parent),
      m_format(format),
      m_seqType(seqType),
      m_data(data) {
}

void ImportTask::execute() {
    emit started(id());
    
    m_progress.status = OperationStatus::Running;
    m_progress.startTime = QDateTime::currentDateTime();
    emit statusChanged(id(), OperationStatus::Running);
    
    try {
        // Create import session
        BiopolymerImportSession importSession(workspace());
        
        // Set progress
        updateStep(1, 3, "Initializing import");
        
        // Import based on format
        BiopolymerState state;
        if (m_format.toUpper() == "SEQUENCE" || m_format.toUpper() == "FASTA") {
            state = importSession.loadSequence(m_seqType, m_data);
        } else {
            // For other formats, use generic import
            state = importSession.importFromText(m_format, m_data);
        }
        
        updateStep(2, 3, "Validating structure");
        
        // Validate the imported structure
        importSession.validate(state);
        
        updateStep(3, 3, "Finalizing");
        
        // Serialize to JSON for result
        QString result = importSession.serialize(state, "JSON");
        
        m_progress.endTime = QDateTime::currentDateTime();
        m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
        
        reportSuccess(result);
        
    } catch (const std::exception& e) {
        reportError(QString("Import failed: %1").arg(e.what()));
    } catch (...) {
        reportError("Import failed: Unknown error");
    }
}

// ========== ExportTask ==========

ExportTask::ExportTask(OperationId id, const BiopolymerWorkspace& workspace,
                       const QString& format, const BiopolymerState& state,
                       QObject* parent)
    : AsyncTask(id, AsyncOperationType::ExportSequence, workspace, parent),
      m_format(format),
      m_state(state) {
}

void ExportTask::execute() {
    emit started(id());
    
    m_progress.status = OperationStatus::Running;
    m_progress.startTime = QDateTime::currentDateTime();
    emit statusChanged(id(), OperationStatus::Running);
    
    try {
        // Create export session
        BiopolymerExportSession exportSession(m_state, workspace());
        
        // Set progress
        updateStep(1, 2, "Preparing export");
        
        // Determine export format
        ExportFormat fmt = exportFormatFromString(m_format);
        
        // Generate export
        ExportResult result = exportSession.generateExport(fmt);
        
        updateStep(2, 2, "Finalizing");
        
        if (result.isSuccess()) {
            m_progress.endTime = QDateTime::currentDateTime();
            m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
            
            reportSuccess(result.outputText);
        } else {
            reportError(result.errorMessage);
        }
        
    } catch (const std::exception& e) {
        reportError(QString("Export failed: %1").arg(e.what()));
    } catch (...) {
        reportError("Export failed: Unknown error");
    }
}

// ========== VectorExportTask ==========

VectorExportTask::VectorExportTask(OperationId id, const BiopolymerWorkspace& workspace,
                                   const QString& format, const BiopolymerState& state,
                                   const QMap<QString, QVariant>& options,
                                   QObject* parent)
    : AsyncTask(id, AsyncOperationType::ExportVector, workspace, parent),
      m_format(format),
      m_state(state),
      m_options(options) {
}

void VectorExportTask::execute() {
    emit started(id());
    
    m_progress.status = OperationStatus::Running;
    m_progress.startTime = QDateTime::currentDateTime();
    emit statusChanged(id(), OperationStatus::Running);
    
    try {
        // Note: This is a placeholder implementation
        // In a real implementation, this would use BiopolymerVectorExport
        
        updateStep(1, 3, "Processing residues");
        
        // Simulate processing
        for (int chainId : m_state.chainIds()) {
            if (shouldStop()) {
                return;
            }
            
            const BiopolymerChain& chain = m_state.chain(chainId);
            for (int residueId : chain.residueIds()) {
                if (shouldStop()) {
                    return;
                }
                
                const BiopolymerResidue& residue = m_state.residue(residueId);
                Q_UNUSED(residue);
            }
        }
        
        updateStep(2, 3, "Generating vector graphics");
        
        // Simulate vector generation
        sleepWithCancellation(100);
        
        updateStep(3, 3, "Finalizing");
        
        // Generate a simple SVG as result
        QString svg = QString("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                             "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"600\">\n"
                             "  <rect x=\"50\" y=\"50\" width=\"700\" height=\"500\" fill=\"white\" stroke=\"black\"/>\n"
                             "  <text x=\"400\" y=\"300\" text-anchor=\"middle\">Vector Export Placeholder</text>\n"
                             "</svg>");
        
        m_progress.endTime = QDateTime::currentDateTime();
        m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
        
        reportSuccess(svg);
        
    } catch (const std::exception& e) {
        reportError(QString("Vector export failed: %1").arg(e.what()));
    } catch (...) {
        reportError("Vector export failed: Unknown error");
    }
}

// ========== ValidationTask ==========

ValidationTask::ValidationTask(OperationId id, const BiopolymerWorkspace& workspace,
                               const BiopolymerState& state, QObject* parent)
    : AsyncTask(id, AsyncOperationType::ValidateStructure, workspace, parent),
      m_state(state) {
}

void ValidationTask::execute() {
    emit started(id());
    
    m_progress.status = OperationStatus::Running;
    m_progress.startTime = QDateTime::currentDateTime();
    emit statusChanged(id(), OperationStatus::Running);
    
    try {
        // Validate the state
        updateStep(1, 3, "Checking chain connectivity");
        
        bool allValid = true;
        QStringList errors;
        
        // Check each chain
        for (int chainId : m_state.chainIds()) {
            if (shouldStop()) {
                return;
            }
            
            const BiopolymerChain& chain = m_state.chain(chainId);
            
            // Check residues
            for (int residueId : chain.residueIds()) {
                if (shouldStop()) {
                    return;
                }
                
                const BiopolymerResidue& residue = m_state.residue(residueId);
                if (!residue.isValid()) {
                    allValid = false;
                    errors << QString("Invalid residue %1 in chain %2").arg(residueId).arg(chainId);
                }
            }
        }
        
        updateStep(2, 3, "Checking connections");
        
        // Check connections
        for (const BiopolymerConnection& conn : m_state.connections()) {
            if (shouldStop()) {
                return;
            }
            
            if (!m_state.hasResidue(conn.fromResidueId()) || 
                !m_state.hasResidue(conn.toResidueId())) {
                allValid = false;
                errors << QString("Invalid connection %1: missing residue").arg(conn.id());
            }
        }
        
        updateStep(3, 3, "Finalizing validation");
        
        QString result;
        if (allValid) {
            result = "Validation passed: All checks passed";
        } else {
            result = QString("Validation failed: %1 issues found").arg(errors.size());
            for (const QString& error : errors) {
                result += QString("\n  - %1").arg(error);
            }
        }
        
        m_progress.endTime = QDateTime::currentDateTime();
        m_progress.elapsedMilliseconds = m_progress.startTime.msecsTo(m_progress.endTime);
        
        reportSuccess(result);
        
    } catch (const std::exception& e) {
        reportError(QString("Validation failed: %1").arg(e.what()));
    } catch (...) {
        reportError("Validation failed: Unknown error");
    }
}

// ========== String Conversions ==========

QString asyncOperationTypeToString(AsyncOperationType type) {
    switch (type) {
    case AsyncOperationType::ImportSequence: return "ImportSequence";
    case AsyncOperationType::ImportFile: return "ImportFile";
    case AsyncOperationType::ExportSequence: return "ExportSequence";
    case AsyncOperationType::ExportVector: return "ExportVector";
    case AsyncOperationType::ExportImage: return "ExportImage";
    case AsyncOperationType::ConvertFormat: return "ConvertFormat";
    case AsyncOperationType::ValidateStructure: return "ValidateStructure";
    case AsyncOperationType::AnalyzeStructure: return "AnalyzeStructure";
    case AsyncOperationType::SearchSubstructure: return "SearchSubstructure";
    case AsyncOperationType::AlignSequences: return "AlignSequences";
    default: return "Unknown";
    }
}

AsyncOperationType asyncOperationTypeFromString(const QString& typeString) {
    const QString lower = typeString.toLower();
    if (lower == "importsequence" || lower == "import_sequence") return AsyncOperationType::ImportSequence;
    if (lower == "importfile" || lower == "import_file") return AsyncOperationType::ImportFile;
    if (lower == "exportsequence" || lower == "export_sequence") return AsyncOperationType::ExportSequence;
    if (lower == "exportvector" || lower == "export_vector") return AsyncOperationType::ExportVector;
    if (lower == "exportimage" || lower == "export_image") return AsyncOperationType::ExportImage;
    if (lower == "convertformat" || lower == "convert_format") return AsyncOperationType::ConvertFormat;
    if (lower == "validatestructure" || lower == "validate_structure") return AsyncOperationType::ValidateStructure;
    if (lower == "analyzestructure" || lower == "analyze_structure") return AsyncOperationType::AnalyzeStructure;
    if (lower == "searchsubstructure" || lower == "search_substructure") return AsyncOperationType::SearchSubstructure;
    if (lower == "alignsequences" || lower == "align_sequences") return AsyncOperationType::AlignSequences;
    return AsyncOperationType::Unknown;
}

QString operationStatusToString(OperationStatus status) {
    switch (status) {
    case OperationStatus::Pending: return "Pending";
    case OperationStatus::Running: return "Running";
    case OperationStatus::Completed: return "Completed";
    case OperationStatus::Failed: return "Failed";
    case OperationStatus::Cancelled: return "Cancelled";
    case OperationStatus::Paused: return "Paused";
    default: return "Unknown";
    }
}
