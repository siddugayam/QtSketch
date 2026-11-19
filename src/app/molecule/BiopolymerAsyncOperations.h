// src/app/molecule/BiopolymerAsyncOperations.h
#ifndef BIOPOLYMERASYNCOPERATIONS_H
#define BIOPOLYMERASYNCOPERATIONS_H

#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include <QString>
#include <QList>
#include <QMap>
#include <QFuture>
#include <QFutureWatcher>
#include <QMutex>
#include <QAtomicInt>
#include <QObject>
#include <functional>

class QThreadPool;

// ========== Async Operation Type ==========

enum class AsyncOperationType {
    ImportSequence,      // Import a sequence
    ImportFile,          // Import a file
    ExportSequence,      // Export to sequence format
    ExportVector,        // Export to vector format
    ExportImage,         // Export to image format
    ConvertFormat,       // Convert between formats
    ValidateStructure,  // Validate biopolymer structure
    AnalyzeStructure,   // Analyze biopolymer structure
    SearchSubstructure,  // Search for substructures
    AlignSequences,      // Sequence alignment
    Unknown
};

// Convert async operation type to string
QString asyncOperationTypeToString(AsyncOperationType type);

// Parse async operation type from string
AsyncOperationType asyncOperationTypeFromString(const QString& typeString);

// ========== Operation ID ==========

// Unique identifier for an async operation
class OperationId {
public:
    OperationId();
    explicit OperationId(quint64 id);
    
    quint64 value() const;
    bool isValid() const;
    
    bool operator==(const OperationId& other) const;
    bool operator!=(const OperationId& other) const;
    bool operator<(const OperationId& other) const;
    
    operator quint64() const;
    
private:
    quint64 m_id;
};

// ========== Operation Status ==========

enum class OperationStatus {
    Pending,     // Operation has not started
    Running,     // Operation is running
    Completed,   // Operation completed successfully
    Failed,      // Operation failed
    Cancelled,   // Operation was cancelled
    Paused       // Operation was paused
};

// Convert operation status to string
QString operationStatusToString(OperationStatus status);

// ========== Progress Info ==========

// Progress information for an async operation
struct ProgressInfo {
    OperationId id;
    AsyncOperationType type = AsyncOperationType::Unknown;
    OperationStatus status = OperationStatus::Pending;
    
    // Progress
    int progress = 0;           // 0-100
    int currentStep = 0;
    int totalSteps = 0;
    QString currentStepDescription;
    
    // Timing
    QDateTime startTime;
    QDateTime endTime;
    qint64 elapsedMilliseconds = 0;
    
    // Result
    QString resultData;
    QString errorMessage;
    QMap<QString, QVariant> metadata;
    
    // Check if operation is active
    bool isActive() const;
    
    // Check if operation is finished
    bool isFinished() const;
    
    // Check if operation is successful
    bool isSuccessful() const;
    
    // Get formatted progress string
    QString progressString() const;
    
    // Get formatted elapsed time string
    QString elapsedString() const;
};

// ========== Cancellation Token ==========

// Token for cooperative cancellation
class CancellationToken {
public:
    CancellationToken();
    
    // Check if cancellation was requested
    bool isCancelled() const;
    
    // Request cancellation
    void cancel();
    
    // Reset the token
    void reset();
    
    // Create a child token that inherits cancellation state
    CancellationToken createChildToken() const;
    
private:
    QAtomicInt m_cancelled;
    QAtomicInt m_parentCancelled;  // Pointer to parent token's cancelled flag (optional)
};

// ========== Operation Context ==========

// Context for an async operation
class OperationContext : public QObject {
    Q_OBJECT
    
public:
    OperationContext(OperationId id, AsyncOperationType type, 
                    const BiopolymerWorkspace& workspace, QObject* parent = nullptr);
    
    // Get operation information
    OperationId id() const;
    AsyncOperationType type() const;
    const BiopolymerWorkspace& workspace() const;
    
    // Get/set cancellation token
    CancellationToken& cancellationToken();
    const CancellationToken& cancellationToken() const;
    
    // Get/set progress
    ProgressInfo& progress();
    const ProgressInfo& progress() const;
    
    // Update progress
    void setProgress(int progress, const QString& description = QString());
    void incrementProgress(const QString& description = QString());
    void setStep(int step, int total, const QString& description = QString());
    
    // Report result
    void setResult(const QString& data, const QMap<QString, QVariant>& metadata = QMap<QString, QVariant>());
    void setError(const QString& error);
    
    // Check if cancelled
    bool isCancelled() const;
    
    // Check if should stop
    bool shouldStop() const;
    
signals:
    void progressChanged(int progress, const QString& description);
    void stepChanged(int current, int total, const QString& description);
    void statusChanged(OperationStatus status);
    void finished(const QString& result, const QString& error);
    
private:
    OperationId m_id;
    AsyncOperationType m_type;
    const BiopolymerWorkspace& m_workspace;
    CancellationToken m_cancellationToken;
    ProgressInfo m_progress;
};

// ========== Async Result ==========

// Result of an async operation
template <typename T = QString>
class AsyncResult {
public:
    AsyncResult();
    AsyncResult(const T& value);
    AsyncResult(T&& value);
    AsyncResult(const QString& error);
    
    // Check if operation succeeded
    bool isSuccess() const;
    
    // Check if operation failed
    bool isFailure() const;
    
    // Check if operation was cancelled
    bool isCancelled() const;
    
    // Get the result value
    const T& value() const;
    T& value();
    
    // Get the error message
    const QString& error() const;
    
    // Set the result
    void setValue(const T& value);
    void setValue(T&& value);
    
    // Set the error
    void setError(const QString& error);
    
    // Set cancelled
    void setCancelled();
    
    // Convert to string
    QString toString() const;
    
private:
    T m_value;
    QString m_error;
    bool m_cancelled = false;
};

// ========== Async Task ==========

// Base class for async tasks
class AsyncTask : public QObject {
    Q_OBJECT
    
public:
    AsyncTask(OperationId id, AsyncOperationType type, 
              const BiopolymerWorkspace& workspace, QObject* parent = nullptr);
    virtual ~AsyncTask();
    
    // Execute the task
    virtual void execute() = 0;
    
    // Get operation ID
    OperationId id() const;
    
    // Get operation type
    AsyncOperationType type() const;
    
    // Get workspace
    const BiopolymerWorkspace& workspace() const;
    
    // Get progress
    ProgressInfo progress() const;
    
    // Get cancellation token
    const CancellationToken& cancellationToken() const;
    
    // Check if cancelled
    bool isCancelled() const;
    
public slots:
    // Cancel the task
    virtual void cancel();
    
    // Pause the task
    virtual void pause();
    
    // Resume the task
    virtual void resume();
    
signals:
    void started(OperationId id);
    void progress(OperationId id, int progress, const QString& description);
    void stepChanged(OperationId id, int current, int total, const QString& description);
    void finished(OperationId id, const QString& result, const QString& error);
    void statusChanged(OperationId id, OperationStatus status);
    
protected:
    // Update progress
    void updateProgress(int progress, const QString& description = QString());
    void updateStep(int current, int total, const QString& description = QString());
    
    // Report completion
    void reportSuccess(const QString& result, const QMap<QString, QVariant>& metadata = QMap<QString, QVariant>());
    void reportError(const QString& error);
    
    // Check if should stop
    bool shouldStop() const;
    
    // Sleep while checking for cancellation
    void sleepWithCancellation(int ms);
    
private:
    OperationId m_id;
    AsyncOperationType m_type;
    const BiopolymerWorkspace& m_workspace;
    ProgressInfo m_progress;
    CancellationToken m_cancellationToken;
};

// ========== BiopolymerAsyncOperations ==========

// Manager for async biopolymer operations
//
// This class manages the execution of expensive operations asynchronously,
// including:
// - Sequence import/export
// - Vector export generation
// - Structure validation
// - Format conversion
//
// Design principles:
// - Operations run in thread pool to avoid blocking the UI
// - Support cooperative cancellation (operations check for cancellation)
// - Cancelled operations cannot publish results
// - Progress reporting for long-running operations
class BiopolymerAsyncOperations : public QObject {
    Q_OBJECT
    
public:
    explicit BiopolymerAsyncOperations(const BiopolymerWorkspace& workspace, 
                                       QObject* parent = nullptr);
    ~BiopolymerAsyncOperations();
    
    // ========== Operation Management ==========
    
    // Start an async operation
    template <typename Task, typename... Args>
    OperationId startOperation(Args&&... args);
    
    // Start an import operation
    OperationId startImport(const QString& format, const QString& seqType, 
                           const QString& data);
    
    // Start an export operation
    OperationId startExport(const QString& format, const BiopolymerState& state);
    
    // Start a vector export operation
    OperationId startVectorExport(const QString& format, const BiopolymerState& state,
                                  const QMap<QString, QVariant>& options);
    
    // Start a validation operation
    OperationId startValidation(const BiopolymerState& state);
    
    // ========== Operation Control ==========
    
    // Cancel an operation
    bool cancelOperation(OperationId id);
    
    // Cancel all operations
    void cancelAllOperations();
    
    // Pause an operation
    bool pauseOperation(OperationId id);
    
    // Resume an operation
    bool resumeOperation(OperationId id);
    
    // ========== Operation Status ==========
    
    // Get operation status
    OperationStatus getOperationStatus(OperationId id) const;
    
    // Get progress for an operation
    ProgressInfo getOperationProgress(OperationId id) const;
    
    // Check if an operation is running
    bool isOperationRunning(OperationId id) const;
    
    // Check if any operations are running
    bool hasRunningOperations() const;
    
    // Get all active operation IDs
    QList<OperationId> activeOperationIds() const;
    
    // Get all operation IDs
    QList<OperationId> allOperationIds() const;
    
    // ========== Cleanup ==========
    
    // Cleanup finished operations
    void cleanupFinishedOperations();
    
    // Clear all operations
    void clearAllOperations();
    
    // ========== Configuration ==========
    
    // Set maximum concurrent operations
    void setMaxConcurrentOperations(int max);
    
    // Get maximum concurrent operations
    int maxConcurrentOperations() const;
    
    // Set thread pool size
    void setThreadPoolSize(int size);
    
    // Get thread pool size
    int threadPoolSize() const;
    
signals:
    // Emitted when an operation starts
    void operationStarted(OperationId id, AsyncOperationType type);
    
    // Emitted when an operation progresses
    void operationProgress(OperationId id, int progress, const QString& description);
    
    // Emitted when an operation's step changes
    void operationStepChanged(OperationId id, int current, int total, const QString& description);
    
    // Emitted when an operation's status changes
    void operationStatusChanged(OperationId id, OperationStatus status);
    
    // Emitted when an operation finishes
    void operationFinished(OperationId id, const QString& result, const QString& error);
    
    // Emitted when an operation is cancelled
    void operationCancelled(OperationId id);
    
    // Emitted when all operations complete
    void allOperationsFinished();
    
private slots:
    // Handle task completion
    void onTaskFinished(OperationId id, const QString& result, const QString& error);
    
    // Handle status change
    void onStatusChanged(OperationId id, OperationStatus status);
    
private:
    // ========== Private Members ==========
    
    const BiopolymerWorkspace& m_workspace;
    QThreadPool* m_threadPool;
    QMap<OperationId, AsyncTask*> m_operations;
    QMutex m_mutex;
    QAtomicInt m_nextId;
    int m_maxConcurrentOperations;
    bool m_ownsThreadPool;
    
    // ========== Private Helper Methods ==========
    
    // Generate a new operation ID
    OperationId generateOperationId();
    
    // Create thread pool if needed
    void ensureThreadPool();
    
    // Get operation by ID
    AsyncTask* getOperation(OperationId id) const;
    
    // Add operation to map
    void addOperation(AsyncTask* task);
    
    // Remove operation from map
    void removeOperation(OperationId id);
};

// ========== Template Implementation ==========

template <typename Task, typename... Args>
OperationId BiopolymerAsyncOperations::startOperation(Args&&... args) {
    ensureThreadPool();
    
    OperationId id = generateOperationId();
    
    // Create the task
    Task* task = new Task(id, m_workspace, this);
    
    // Connect signals
    connect(task, &AsyncTask::started, this, [this, id]() {
        emit operationStarted(id, task->type());
    });
    connect(task, &AsyncTask::progress, this, [this](OperationId id, int progress, const QString& desc) {
        emit operationProgress(id, progress, desc);
    });
    connect(task, &AsyncTask::stepChanged, this, [this](OperationId id, int current, int total, const QString& desc) {
        emit operationStepChanged(id, current, total, desc);
    });
    connect(task, &AsyncTask::statusChanged, this, &BiopolymerAsyncOperations::onStatusChanged);
    connect(task, &AsyncTask::finished, this, &BiopolymerAsyncOperations::onTaskFinished);
    
    // Add to map
    addOperation(task);
    
    // Start the task in the thread pool
    QFuture<void> future = QtConcurrent::run([task]() {
        task->execute();
    });
    
    // Store future for potential cancellation
    Q_UNUSED(future);
    
    return id;
}

// ========== Specific Task Classes ==========

// Import task
class ImportTask : public AsyncTask {
    Q_OBJECT
    
public:
    ImportTask(OperationId id, const BiopolymerWorkspace& workspace, 
              const QString& format, const QString& seqType, const QString& data,
              QObject* parent = nullptr);
    
    void execute() override;
    
private:
    QString m_format;
    QString m_seqType;
    QString m_data;
};

// Export task
class ExportTask : public AsyncTask {
    Q_OBJECT
    
public:
    ExportTask(OperationId id, const BiopolymerWorkspace& workspace,
              const QString& format, const BiopolymerState& state,
              QObject* parent = nullptr);
    
    void execute() override;
    
private:
    QString m_format;
    BiopolymerState m_state;
};

// Vector export task
class VectorExportTask : public AsyncTask {
    Q_OBJECT
    
public:
    VectorExportTask(OperationId id, const BiopolymerWorkspace& workspace,
                    const QString& format, const BiopolymerState& state,
                    const QMap<QString, QVariant>& options,
                    QObject* parent = nullptr);
    
    void execute() override;
    
private:
    QString m_format;
    BiopolymerState m_state;
    QMap<QString, QVariant> m_options;
};

// Validation task
class ValidationTask : public AsyncTask {
    Q_OBJECT
    
public:
    ValidationTask(OperationId id, const BiopolymerWorkspace& workspace,
                  const BiopolymerState& state, QObject* parent = nullptr);
    
    void execute() override;
    
private:
    BiopolymerState m_state;
};

#endif // BIOPOLYMERASYNCOPERATIONS_H
