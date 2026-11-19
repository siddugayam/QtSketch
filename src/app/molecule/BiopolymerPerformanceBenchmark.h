// src/app/molecule/BiopolymerPerformanceBenchmark.h
#ifndef BIOPOLYMERPERFORMANCEBENCHMARK_H
#define BIOPOLYMERPERFORMANCEBENCHMARK_H

#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include <QString>
#include <QList>
#include <QMap>
#include <QElapsedTimer>
#include <QDateTime>

// ========== Benchmark Configuration ==========

// Configuration for performance benchmarks
struct BenchmarkConfig {
    // Test sizes (number of residues)
    QList<int> testSizes = {100, 1000, 10000};
    
    // Operations to benchmark
    bool testDocumentLoad = true;
    bool testSequenceRendering = true;
    bool testExportGeneration = true;
    bool testSelectionFeedback = true;
    bool testZoomPan = true;
    bool testUndoRedo = true;
    
    // Iteration counts (higher for smaller documents)
    int smallDocumentIterations = 100;    // < 1000 residues
    int mediumDocumentIterations = 10;    // 1000-5000 residues
    int largeDocumentIterations = 1;      // > 5000 residues
    
    // Warmup iterations
    int warmupIterations = 3;
    
    // Target performance thresholds (milliseconds)
    int targetSelectionFeedbackMs = 100;   // Target: < 100ms
    int targetFrameBudgetMs = 33;          // Target: < 33ms per frame
    int targetLoadTimeMs = 1000;           // Target: < 1s for 1000 residues
    
    // Memory tracking
    bool trackMemory = true;
    
    // Output options
    bool verboseOutput = false;
    bool saveResults = true;
    QString resultsFilePath;
    
    // Serialization
    QMap<QString, QVariant> toJson() const;
    static BenchmarkConfig fromJson(const QMap<QString, QVariant>& json);
};

// ========== Benchmark Result ==========

// Result of a single benchmark measurement
struct BenchmarkResult {
    QString testName;
    int testSize = 0;                    // Number of residues
    QString operation;
    int iteration = 0;
    
    // Timing
    qint64 elapsedNanoseconds = 0;
    qint64 elapsedMilliseconds = 0;
    double elapsedSeconds = 0.0;
    
    // Memory (if tracked)
    qint64 memoryBeforeBytes = 0;
    qint64 memoryAfterBytes = 0;
    qint64 memoryDeltaBytes = 0;
    
    // Statistics
    double averageMs = 0.0;
    double minMs = 0.0;
    double maxMs = 0.0;
    double stdDevMs = 0.0;
    
    // Pass/fail status
    bool passed = false;
    QStringList notes;
    
    // Formatted output
    QString formattedResult() const;
    QString csvLine() const;
    
    // Check if result meets target
    bool meetsTarget(int targetMs) const;
};

// ========== Benchmark Suite Result ==========

// Result of a complete benchmark suite
struct BenchmarkSuiteResult {
    QString suiteName;
    QDateTime timestamp;
    BenchmarkConfig config;
    
    // Results per test size
    QMap<int, QList<BenchmarkResult>> resultsBySize;
    
    // Summary statistics
    QMap<int, QMap<QString, BenchmarkResult>> summariesBySizeAndOp;
    
    // Overall pass/fail
    bool overallPass = false;
    QList<QString> failures;
    QList<QString> warnings;
    
    // System information
    QString systemInfo;
    
    // Formatted output
    QString formattedReport() const;
    QString csvReport() const;
    
    // Check if all benchmarks passed
    bool allPassed() const;
    
    // Get worst performing operations
    QList<BenchmarkResult> worstPerformers(int limit = 5) const;
};

// ========== Benchmark Operation ==========

// Represents a benchmarkable operation
struct BenchmarkOperation {
    QString name;
    QString description;
    std::function<BenchmarkResult(const BiopolymerState&, int)> executor;
    int targetMs = -1;  // -1 = no target
    
    // Run the benchmark
    QList<BenchmarkResult> run(const BiopolymerState& state, int size, 
                               int iterations, int warmup = 0) const;
};

// ========== BiopolymerPerformanceBenchmark ==========

// Class for running performance benchmarks on biopolymer operations
//
// This class provides comprehensive benchmarking for:
// - Document loading and parsing
// - Sequence rendering performance
// - Export generation performance
// - Selection feedback latency
// - Zoom/pan responsiveness
// - Undo/redo performance
//
// Design principle: Benchmarks should be isolated and repeatable.
// Each benchmark runs multiple iterations and reports statistics.
class BiopolymerPerformanceBenchmark {
public:
    // Constructor
    explicit BiopolymerPerformanceBenchmark(const BiopolymerWorkspace& workspace);
    
    // Initialize with configuration
    void initialize(const BenchmarkConfig& config);
    
    // ========== Benchmark Execution ==========
    
    // Run all benchmarks
    BenchmarkSuiteResult runAllBenchmarks();
    
    // Run benchmarks for specific document sizes
    BenchmarkSuiteResult runBenchmarksForSizes(const QList<int>& sizes);
    
    // Run a specific benchmark
    QList<BenchmarkResult> runBenchmark(const QString& benchmarkName, int size);
    
    // ========== Individual Benchmarks ==========
    
    // Benchmark document loading
    BenchmarkResult benchmarkDocumentLoad(int size);
    
    // Benchmark sequence rendering
    BenchmarkResult benchmarkSequenceRendering(int size);
    
    // Benchmark export generation
    BenchmarkResult benchmarkExportGeneration(int size);
    
    // Benchmark selection feedback
    BenchmarkResult benchmarkSelectionFeedback(int size);
    
    // Benchmark zoom/pan
    BenchmarkResult benchmarkZoomPan(int size);
    
    // Benchmark undo/redo
    BenchmarkResult benchmarkUndoRedo(int size);
    
    // ========== Benchmark Generation ==========
    
    // Generate a biopolymer of specific size for benchmarking
    BiopolymerState generateTestDocument(int size, PolymerType type = PolymerType::Protein);
    
    // Generate a random protein sequence
    QString generateRandomProtein(int length);
    
    // Generate a random DNA sequence
    QString generateRandomDna(int length);
    
    // Generate a random RNA sequence
    QString generateRandomRna(int length);
    
    // ========== Analysis ==========
    
    // Analyze results and identify bottlenecks
    QString analyzeResults(const BenchmarkSuiteResult& result);
    
    // Get performance recommendations
    QStringList getRecommendations(const BenchmarkSuiteResult& result);
    
    // ========== System Information ==========
    
    // Get system information for benchmark report
    QString getSystemInfo() const;
    
    // Get memory usage
    qint64 getMemoryUsage() const;
    
    // ========== Configuration ==========
    
    // Get/set configuration
    const BenchmarkConfig& config() const;
    void setConfig(const BenchmarkConfig& config);
    
    // Get workspace
    const BiopolymerWorkspace& workspace() const;
    
    // ========== Static Helpers ==========
    
    // Get default configuration
    static BenchmarkConfig defaultConfig();
    
    // Create a benchmark suite with custom operations
    static BenchmarkSuiteResult runCustomBenchmark(
        const QString& suiteName,
        const QList<BenchmarkOperation>& operations,
        const QList<int>& sizes,
        const BenchmarkConfig& config = BenchmarkConfig()
    );
    
    // Format time in milliseconds
    static QString formatTimeMs(qint64 ms);
    
    // Format memory size
    static QString formatMemorySize(qint64 bytes);

private:
    // ========== Private Members ==========
    
    const BiopolymerWorkspace& m_workspace;
    BenchmarkConfig m_config;
    BenchmarkSuiteResult m_lastResult;
    
    // ========== Private Helper Methods ==========
    
    // Run a benchmark operation multiple times
    QList<BenchmarkResult> runBenchmarkOperation(
        const QString& testName,
        const QString& operation,
        const std::function<void()>& operationFunc,
        int size,
        int iterations,
        int warmup
    );
    
    // Measure single operation execution
    BenchmarkResult measureOperation(
        const QString& testName,
        const QString& operation,
        const std::function<void()>& operationFunc,
        int size,
        int iteration
    );
    
    // Warmup the system
    void warmup(int iterations);
    
    // Calculate statistics from multiple results
    BenchmarkResult calculateStatistics(const QList<BenchmarkResult>& results) const;
    
    // Check if system is responsive
    bool checkSystemResponsiveness() const;
    
    // Create a test state with the specified size
    BiopolymerState createTestState(int size);
};

// ========== Inline Implementations ==========

inline bool BenchmarkResult::meetsTarget(int targetMs) const {
    return elapsedMilliseconds <= targetMs;
}

inline bool BenchmarkSuiteResult::allPassed() const {
    return overallPass && failures.isEmpty();
}

#endif // BIOPOLYMERPERFORMANCEBENCHMARK_H
