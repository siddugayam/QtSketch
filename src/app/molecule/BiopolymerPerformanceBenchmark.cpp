// src/app/molecule/BiopolymerPerformanceBenchmark.cpp
#include "BiopolymerPerformanceBenchmark.h"
#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include "BiopolymerImportSession.h"
#include <QDebug>
#include <QCoreApplication>
#include <QRandomGenerator>
#include <QProcess>

#ifdef Q_OS_WIN
#include <windows.h>
#include <psapi.h>
#endif

// ========== BenchmarkConfig ==========

QMap<QString, QVariant> BenchmarkConfig::toJson() const {
    QMap<QString, QVariant> json;
    // Add test sizes
    QStringList sizesStr;
    for (int size : testSizes) {
        sizesStr << QString::number(size);
    }
    json["testSizes"] = sizesStr;
    
    json["testDocumentLoad"] = testDocumentLoad;
    json["testSequenceRendering"] = testSequenceRendering;
    json["testExportGeneration"] = testExportGeneration;
    json["testSelectionFeedback"] = testSelectionFeedback;
    json["testZoomPan"] = testZoomPan;
    json["testUndoRedo"] = testUndoRedo;
    json["smallDocumentIterations"] = smallDocumentIterations;
    json["mediumDocumentIterations"] = mediumDocumentIterations;
    json["largeDocumentIterations"] = largeDocumentIterations;
    json["warmupIterations"] = warmupIterations;
    json["targetSelectionFeedbackMs"] = targetSelectionFeedbackMs;
    json["targetFrameBudgetMs"] = targetFrameBudgetMs;
    json["targetLoadTimeMs"] = targetLoadTimeMs;
    json["trackMemory"] = trackMemory;
    json["verboseOutput"] = verboseOutput;
    json["saveResults"] = saveResults;
    json["resultsFilePath"] = resultsFilePath;
    return json;
}

BenchmarkConfig BenchmarkConfig::fromJson(const QMap<QString, QVariant>& json) {
    BenchmarkConfig config;
    
    if (json.contains("testSizes")) {
        QStringList sizesStr = json["testSizes"].toStringList();
        config.testSizes.clear();
        for (const QString& sizeStr : sizesStr) {
            config.testSizes << sizeStr.toInt();
        }
    }
    
    if (json.contains("testDocumentLoad")) {
        config.testDocumentLoad = json["testDocumentLoad"].toBool();
    }
    if (json.contains("testSequenceRendering")) {
        config.testSequenceRendering = json["testSequenceRendering"].toBool();
    }
    if (json.contains("testExportGeneration")) {
        config.testExportGeneration = json["testExportGeneration"].toBool();
    }
    if (json.contains("testSelectionFeedback")) {
        config.testSelectionFeedback = json["testSelectionFeedback"].toBool();
    }
    if (json.contains("testZoomPan")) {
        config.testZoomPan = json["testZoomPan"].toBool();
    }
    if (json.contains("testUndoRedo")) {
        config.testUndoRedo = json["testUndoRedo"].toBool();
    }
    
    if (json.contains("smallDocumentIterations")) {
        config.smallDocumentIterations = json["smallDocumentIterations"].toInt();
    }
    if (json.contains("mediumDocumentIterations")) {
        config.mediumDocumentIterations = json["mediumDocumentIterations"].toInt();
    }
    if (json.contains("largeDocumentIterations")) {
        config.largeDocumentIterations = json["largeDocumentIterations"].toInt();
    }
    if (json.contains("warmupIterations")) {
        config.warmupIterations = json["warmupIterations"].toInt();
    }
    
    if (json.contains("targetSelectionFeedbackMs")) {
        config.targetSelectionFeedbackMs = json["targetSelectionFeedbackMs"].toInt();
    }
    if (json.contains("targetFrameBudgetMs")) {
        config.targetFrameBudgetMs = json["targetFrameBudgetMs"].toInt();
    }
    if (json.contains("targetLoadTimeMs")) {
        config.targetLoadTimeMs = json["targetLoadTimeMs"].toInt();
    }
    
    if (json.contains("trackMemory")) {
        config.trackMemory = json["trackMemory"].toBool();
    }
    if (json.contains("verboseOutput")) {
        config.verboseOutput = json["verboseOutput"].toBool();
    }
    if (json.contains("saveResults")) {
        config.saveResults = json["saveResults"].toBool();
    }
    if (json.contains("resultsFilePath")) {
        config.resultsFilePath = json["resultsFilePath"].toString();
    }
    
    return config;
}

// ========== BenchmarkResult ==========

QString BenchmarkResult::formattedResult() const {
    QStringList parts;
    parts << QString("Test: %1").arg(testName);
    parts << QString("Size: %1 residues").arg(testSize);
    parts << QString("Operation: %1").arg(operation);
    parts << QString("Elapsed: %1 ms").arg(elapsedMilliseconds);
    
    if (memoryDeltaBytes != 0) {
        parts << QString("Memory delta: %1").arg(formatMemorySize(memoryDeltaBytes));
    }
    
    if (passed) {
        parts << "Status: PASS";
    } else {
        parts << "Status: FAIL";
    }
    
    if (!notes.isEmpty()) {
        parts << QString("Notes: %1").arg(notes.join(", "));
    }
    
    return parts.join(" | ");
}

QString BenchmarkResult::csvLine() const {
    return QString("%1,%2,%3,%4,%5,%6")
        .arg(testName)
        .arg(testSize)
        .arg(operation)
        .arg(elapsedMilliseconds)
        .arg(memoryDeltaBytes)
        .arg(passed ? "PASS" : "FAIL");
}

// ========== BenchmarkSuiteResult ==========

QString BenchmarkSuiteResult::formattedReport() const {
    QStringList lines;
    
    lines << QString("== Performance Benchmark Report ==");
    lines << QString("Suite: %1").arg(suiteName);
    lines << QString("Timestamp: %1").arg(timestamp.toString(Qt::ISODate));
    lines << QString("System: %1").arg(systemInfo);
    lines << "";
    
    // Summary by size
    for (int size : resultsBySize.keys()) {
        lines << QString("--- Size: %1 residues ---").arg(size);
        
        if (summariesBySizeAndOp.contains(size)) {
            for (const BenchmarkResult& summary : summariesBySizeAndOp[size]) {
                lines << QString("  %1: %2 ms (avg)").arg(summary.operation).arg(summary.averageMs);
            }
        }
        lines << "";
    }
    
    // Overall status
    lines << QString("Overall: %1").arg(overallPass ? "PASS" : "FAIL");
    
    if (!failures.isEmpty()) {
        lines << "Failures:";
        for (const QString& failure : failures) {
            lines << QString("  - %1").arg(failure);
        }
    }
    
    if (!warnings.isEmpty()) {
        lines << "Warnings:";
        for (const QString& warning : warnings) {
            lines << QString("  - %1").arg(warning);
        }
    }
    
    return lines.join("\n");
}

QString BenchmarkSuiteResult::csvReport() const {
    QStringList lines;
    lines << "Suite,Timestamp,Size,Operation,ElapsedMs,MemoryDelta,Status";
    
    for (int size : resultsBySize.keys()) {
        for (const BenchmarkResult& result : resultsBySize[size]) {
            lines << result.csvLine();
        }
    }
    
    return lines.join("\n");
}

QList<BenchmarkResult> BenchmarkSuiteResult::worstPerformers(int limit) const {
    QList<BenchmarkResult> allResults;
    
    for (int size : resultsBySize.keys()) {
        allResults << resultsBySize[size];
    }
    
    // Sort by elapsed time (descending)
    std::sort(allResults.begin(), allResults.end(), 
        [](const BenchmarkResult& a, const BenchmarkResult& b) {
            return a.elapsedMilliseconds > b.elapsedMilliseconds;
        });
    
    return allResults.mid(0, limit);
}

// ========== BenchmarkOperation ==========

QList<BenchmarkResult> BenchmarkOperation::run(
    const BiopolymerState& state, int size, int iterations, int warmup) const {
    
    QList<BenchmarkResult> results;
    
    // Warmup
    for (int i = 0; i < warmup; ++i) {
        executor(state, size);
    }
    
    // Run iterations
    for (int i = 0; i < iterations; ++i) {
        BenchmarkResult result;
        result.testName = name;
        result.operation = description;
        result.testSize = size;
        result.iteration = i + 1;
        
        QElapsedTimer timer;
        timer.start();
        
        result = executor(state, size);
        
        timer.invalidate();
        result.elapsedNanoseconds = timer.nsecsElapsed();
        result.elapsedMilliseconds = timer.msecsElapsed();
        result.elapsedSeconds = timer.nsecsElapsed() / 1000000000.0;
        
        results.append(result);
    }
    
    return results;
}

// ========== BiopolymerPerformanceBenchmark ==========

BiopolymerPerformanceBenchmark::BiopolymerPerformanceBenchmark(
    const BiopolymerWorkspace& workspace)
    : m_workspace(workspace) {
    m_config = defaultConfig();
}

void BiopolymerPerformanceBenchmark::initialize(const BenchmarkConfig& config) {
    m_config = config;
}

BenchmarkConfig BiopolymerPerformanceBenchmark::defaultConfig() {
    BenchmarkConfig config;
    config.testSizes = {100, 1000, 10000};
    config.testDocumentLoad = true;
    config.testSequenceRendering = true;
    config.testExportGeneration = true;
    config.testSelectionFeedback = true;
    config.testZoomPan = true;
    config.testUndoRedo = true;
    config.smallDocumentIterations = 100;
    config.mediumDocumentIterations = 10;
    config.largeDocumentIterations = 1;
    config.warmupIterations = 3;
    config.targetSelectionFeedbackMs = 100;
    config.targetFrameBudgetMs = 33;
    config.targetLoadTimeMs = 1000;
    config.trackMemory = true;
    config.verboseOutput = false;
    config.saveResults = true;
    return config;
}

const BenchmarkConfig& BiopolymerPerformanceBenchmark::config() const {
    return m_config;
}

void BiopolymerPerformanceBenchmark::setConfig(const BenchmarkConfig& config) {
    m_config = config;
}

const BiopolymerWorkspace& BiopolymerPerformanceBenchmark::workspace() const {
    return m_workspace;
}

BenchmarkSuiteResult BiopolymerPerformanceBenchmark::runAllBenchmarks() {
    return runBenchmarksForSizes(m_config.testSizes);
}

BenchmarkSuiteResult BiopolymerPerformanceBenchmark::runBenchmarksForSizes(
    const QList<int>& sizes) {
    
    BenchmarkSuiteResult suiteResult;
    suiteResult.suiteName = "Biopolymer Performance Benchmark";
    suiteResult.timestamp = QDateTime::currentDateTime();
    suiteResult.config = m_config;
    suiteResult.systemInfo = getSystemInfo();
    
    for (int size : sizes) {
        QList<BenchmarkResult> results;
        
        // Create test document
        BiopolymerState testState = generateTestDocument(size);
        
        if (m_config.testDocumentLoad) {
            results << runBenchmark("Document Load", size);
        }
        
        if (m_config.testSequenceRendering) {
            results << runBenchmark("Sequence Rendering", size);
        }
        
        if (m_config.testExportGeneration) {
            results << runBenchmark("Export Generation", size);
        }
        
        if (m_config.testSelectionFeedback) {
            results << runBenchmark("Selection Feedback", size);
        }
        
        if (m_config.testZoomPan) {
            results << runBenchmark("Zoom/Pan", size);
        }
        
        if (m_config.testUndoRedo) {
            results << runBenchmark("Undo/Redo", size);
        }
        
        suiteResult.resultsBySize[size] = results;
        
        // Calculate summaries
        QMap<QString, BenchmarkResult> summaries;
        for (const BenchmarkResult& result : results) {
            if (!summaries.contains(result.operation)) {
                summaries[result.operation] = result;
            } else {
                // Accumulate
                BenchmarkResult& summary = summaries[result.operation];
                summary.elapsedMilliseconds += result.elapsedMilliseconds;
                summary.iteration++;
            }
        }
        
        // Calculate averages
        for (BenchmarkResult& summary : summaries) {
            if (summary.iteration > 0) {
                summary.averageMs = static_cast<double>(summary.elapsedMilliseconds) / summary.iteration;
            }
        }
        
        suiteResult.summariesBySizeAndOp[size] = summaries.values();
    }
    
    // Check overall pass/fail
    suiteResult.overallPass = true;
    for (int size : sizes) {
        if (suiteResult.summariesBySizeAndOp.contains(size)) {
            for (const BenchmarkResult& summary : suiteResult.summariesBySizeAndOp[size]) {
                if (summary.operation == "Selection Feedback" && 
                    !summary.meetsTarget(m_config.targetSelectionFeedbackMs)) {
                    suiteResult.overallPass = false;
                    suiteResult.failures << QString("Selection feedback too slow for %1 residues: %2 ms > %3 ms")
                        .arg(size).arg(summary.averageMs).arg(m_config.targetSelectionFeedbackMs);
                }
                if (summary.operation == "Sequence Rendering" &&
                    !summary.meetsTarget(m_config.targetFrameBudgetMs * 3)) { // Allow 3 frames
                    suiteResult.warnings << QString("Sequence rendering borderline for %1 residues: %2 ms")
                        .arg(size).arg(summary.averageMs);
                }
            }
        }
    }
    
    m_lastResult = suiteResult;
    return suiteResult;
}

QList<BenchmarkResult> BiopolymerPerformanceBenchmark::runBenchmark(
    const QString& benchmarkName, int size) {
    
    QList<BenchmarkResult> results;
    
    // Determine iteration count based on size
    int iterations = m_config.smallDocumentIterations;
    if (size >= 1000 && size < 5000) {
        iterations = m_config.mediumDocumentIterations;
    } else if (size >= 5000) {
        iterations = m_config.largeDocumentIterations;
    }
    
    // Create test state
    BiopolymerState testState = generateTestDocument(size);
    
    if (benchmarkName == "Document Load") {
        results << benchmarkDocumentLoad(size);
    } else if (benchmarkName == "Sequence Rendering") {
        results << benchmarkSequenceRendering(size);
    } else if (benchmarkName == "Export Generation") {
        results << benchmarkExportGeneration(size);
    } else if (benchmarkName == "Selection Feedback") {
        results << benchmarkSelectionFeedback(size);
    } else if (benchmarkName == "Zoom/Pan") {
        results << benchmarkZoomPan(size);
    } else if (benchmarkName == "Undo/Redo") {
        results << benchmarkUndoRedo(size);
    }
    
    return results;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkDocumentLoad(int size) {
    BenchmarkResult result;
    result.testName = "Document Load";
    result.operation = "Load from sequence";
    result.testSize = size;
    
    // Generate random sequence
    QString sequence = generateRandomProtein(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Create import session and load
    BiopolymerImportSession importSession(m_workspace);
    BiopolymerState state = importSession.loadSequence("PEPTIDE", sequence);
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.meetsTarget(m_config.targetLoadTimeMs);
    
    if (!result.passed) {
        result.notes << QString("Target: %1 ms").arg(m_config.targetLoadTimeMs);
    }
    
    return result;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkSequenceRendering(int size) {
    BenchmarkResult result;
    result.testName = "Sequence Rendering";
    result.operation = "Render to canvas";
    result.testSize = size;
    
    BiopolymerState testState = generateTestDocument(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Simulate rendering by iterating through residues
    // In a real implementation, this would call the actual rendering code
    for (int chainId : testState.chainIds()) {
        const BiopolymerChain& chain = testState.chain(chainId);
        for (int residueId : chain.residueIds()) {
            const BiopolymerResidue& residue = testState.residue(residueId);
            // Access residue properties to simulate rendering
            Q_UNUSED(residue);
            Q_UNUSED(residueId);
            Q_UNUSED(chain);
            Q_UNUSED(chainId);
        }
    }
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.meetsTarget(m_config.targetFrameBudgetMs * 3); // Allow multiple frames
    
    return result;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkExportGeneration(int size) {
    BenchmarkResult result;
    result.testName = "Export Generation";
    result.operation = "Generate FASTA export";
    result.testSize = size;
    
    BiopolymerState testState = generateTestDocument(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Simulate export generation
    QStringList sequences;
    for (int chainId : testState.chainIds()) {
        const BiopolymerChain& chain = testState.chain(chainId);
        QString seq;
        for (int residueId : chain.residueIds()) {
            const BiopolymerResidue& residue = testState.residue(residueId);
            seq += residue.label();
        }
        sequences << seq;
    }
    
    // Combine into FASTA
    QString fasta;
    for (int i = 0; i < sequences.size(); ++i) {
        fasta += QString(">Chain %1\n%2\n").arg(i + 1).arg(sequences[i]);
    }
    
    Q_UNUSED(fasta);
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.elapsedMilliseconds < 500; // Target: < 500ms
    
    return result;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkSelectionFeedback(int size) {
    BenchmarkResult result;
    result.testName = "Selection Feedback";
    result.operation = "Select residue";
    result.testSize = size;
    
    BiopolymerState testState = generateTestDocument(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Simulate selection by accessing a residue
    if (!testState.chainIds().isEmpty()) {
        int chainId = testState.chainIds().first();
        const BiopolymerChain& chain = testState.chain(chainId);
        if (!chain.residueIds().isEmpty()) {
            int residueId = chain.residueIds().first();
            const BiopolymerResidue& residue = testState.residue(residueId);
            Q_UNUSED(residue);
            Q_UNUSED(residueId);
            Q_UNUSED(chain);
            Q_UNUSED(chainId);
        }
    }
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.meetsTarget(m_config.targetSelectionFeedbackMs);
    
    if (!result.passed) {
        result.notes << QString("Target: %1 ms").arg(m_config.targetSelectionFeedbackMs);
    }
    
    return result;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkZoomPan(int size) {
    BenchmarkResult result;
    result.testName = "Zoom/Pan";
    result.operation = "Apply transform";
    result.testSize = size;
    
    BiopolymerState testState = generateTestDocument(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Simulate zoom/pan by recalculating positions
    for (int chainId : testState.chainIds()) {
        const BiopolymerChain& chain = testState.chain(chainId);
        for (int residueId : chain.residueIds()) {
            const BiopolymerResidue& residue = testState.residue(residueId);
            // Simulate position calculation
            qreal x = residue.x();
            qreal y = residue.y();
            qreal scaledX = x * 1.5; // Simulate zoom
            qreal scaledY = y * 1.5;
            Q_UNUSED(scaledX);
            Q_UNUSED(scaledY);
            Q_UNUSED(residue);
            Q_UNUSED(residueId);
            Q_UNUSED(chain);
            Q_UNUSED(chainId);
        }
    }
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.meetsTarget(m_config.targetFrameBudgetMs);
    
    return result;
}

BenchmarkResult BiopolymerPerformanceBenchmark::benchmarkUndoRedo(int size) {
    BenchmarkResult result;
    result.testName = "Undo/Redo";
    result.operation = "Undo last operation";
    result.testSize = size;
    
    BiopolymerState testState = generateTestDocument(size);
    
    QElapsedTimer timer;
    timer.start();
    
    // Simulate undo by copying state (in real implementation, this would be a proper undo)
    BiopolymerState copy = testState;
    Q_UNUSED(copy);
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.passed = result.elapsedMilliseconds < 100; // Target: < 100ms
    
    return result;
}

BiopolymerState BiopolymerPerformanceBenchmark::generateTestDocument(
    int size, PolymerType type) {
    
    BiopolymerWorkspace tempWorkspace;
    BiopolymerState state(tempWorkspace);
    
    // Create a chain with the specified number of residues
    QString sequence;
    switch (type) {
    case PolymerType::Protein:
        sequence = generateRandomProtein(size);
        break;
    case PolymerType::Dna:
        sequence = generateRandomDna(size);
        break;
    case PolymerType::Rna:
        sequence = generateRandomRna(size);
        break;
    default:
        sequence = generateRandomProtein(size);
    }
    
    // Import the sequence
    BiopolymerImportSession importSession(tempWorkspace);
    state = importSession.loadSequence("PEPTIDE", sequence);
    
    return state;
}

QString BiopolymerPerformanceBenchmark::generateRandomProtein(int length) {
    static const QString aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
    QString result;
    
    QRandomGenerator rng(QDateTime::currentDateTime().toMSecsSinceEpoch());
    
    for (int i = 0; i < length; ++i) {
        int index = rng.bounded(aminoAcids.length());
        result += aminoAcids[index];
    }
    
    return result;
}

QString BiopolymerPerformanceBenchmark::generateRandomDna(int length) {
    static const QString dnaBases = "ACGT";
    QString result;
    
    QRandomGenerator rng(QDateTime::currentDateTime().toMSecsSinceEpoch());
    
    for (int i = 0; i < length; ++i) {
        int index = rng.bounded(dnaBases.length());
        result += dnaBases[index];
    }
    
    return result;
}

QString BiopolymerPerformanceBenchmark::generateRandomRna(int length) {
    static const QString rnaBases = "ACGU";
    QString result;
    
    QRandomGenerator rng(QDateTime::currentDateTime().toMSecsSinceEpoch());
    
    for (int i = 0; i < length; ++i) {
        int index = rng.bounded(rnaBases.length());
        result += rnaBases[index];
    }
    
    return result;
}

QString BiopolymerPerformanceBenchmark::getSystemInfo() const {
    QStringList info;
    
    // OS
    info << QString("OS: %1").arg(QSysInfo::productType());
    
    // CPU
    info << QString("CPU cores: %1").arg(QThread::idealThreadCount());
    
    // Memory
    info << QString("Memory: %1").arg(formatMemorySize(getMemoryUsage()));
    
    // Qt version
    info << QString("Qt: %1").arg(QT_VERSION_STR);
    
    return info.join("; ");
}

qint64 BiopolymerPerformanceBenchmark::getMemoryUsage() const {
#ifdef Q_OS_WIN
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    return pmc.WorkingSetSize;
#elif defined(Q_OS_LINUX)
    // Read /proc/self/status
    QFile file("/proc/self/status");
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);
        while (!in.atEnd()) {
            QString line = in.readLine();
            if (line.startsWith("VmRSS:")) {
                QStringList parts = line.split(" ", Qt::SkipEmptyParts);
                if (parts.size() >= 2) {
                    return parts[1].toLongLong() * 1024; // Convert from KB to bytes
                }
            }
        }
    }
    return 0;
#else
    return 0;
#endif
}

BenchmarkSuiteResult BiopolymerPerformanceBenchmark::runCustomBenchmark(
    const QString& suiteName, const QList<BenchmarkOperation>& operations,
    const QList<int>& sizes, const BenchmarkConfig& config) {
    
    BiopolymerPerformanceBenchmark benchmark(BiopolymerWorkspace());
    benchmark.m_config = config;
    
    BenchmarkSuiteResult suiteResult;
    suiteResult.suiteName = suiteName;
    suiteResult.timestamp = QDateTime::currentDateTime();
    suiteResult.config = config;
    suiteResult.systemInfo = benchmark.getSystemInfo();
    
    for (int size : sizes) {
        BiopolymerState testState = benchmark.generateTestDocument(size);
        
        for (const BenchmarkOperation& op : operations) {
            QList<BenchmarkResult> results = op.run(testState, size, 
                benchmark.getIterationCount(size), config.warmupIterations);
            
            // Calculate statistics
            BenchmarkResult summary = benchmark.calculateStatistics(results);
            summary.testSize = size;
            summary.testName = op.name;
            summary.operation = op.description;
            
            suiteResult.resultsBySize[size].append(summary);
        }
    }
    
    return suiteResult;
}

int BiopolymerPerformanceBenchmark::getIterationCount(int size) const {
    if (size < 1000) {
        return m_config.smallDocumentIterations;
    } else if (size < 5000) {
        return m_config.mediumDocumentIterations;
    } else {
        return m_config.largeDocumentIterations;
    }
}

QList<BenchmarkResult> BiopolymerPerformanceBenchmark::runBenchmarkOperation(
    const QString& testName, const QString& operation, 
    const std::function<void()>& operationFunc, int size, 
    int iterations, int warmup) {
    
    QList<BenchmarkResult> results;
    
    // Warmup
    for (int i = 0; i < warmup; ++i) {
        operationFunc();
    }
    
    // Run iterations
    for (int i = 0; i < iterations; ++i) {
        results << measureOperation(testName, operation, operationFunc, size, i + 1);
    }
    
    return results;
}

BenchmarkResult BiopolymerPerformanceBenchmark::measureOperation(
    const QString& testName, const QString& operation,
    const std::function<void()>& operationFunc, int size, int iteration) {
    
    BenchmarkResult result;
    result.testName = testName;
    result.operation = operation;
    result.testSize = size;
    result.iteration = iteration;
    
    // Track memory if enabled
    qint64 memoryBefore = 0;
    if (m_config.trackMemory) {
        memoryBefore = getMemoryUsage();
    }
    
    QElapsedTimer timer;
    timer.start();
    
    operationFunc();
    
    timer.invalidate();
    
    result.elapsedNanoseconds = timer.nsecsElapsed();
    result.elapsedMilliseconds = timer.msecsElapsed();
    result.elapsedSeconds = timer.nsecsElapsed() / 1000000000.0;
    
    // Track memory if enabled
    if (m_config.trackMemory) {
        qint64 memoryAfter = getMemoryUsage();
        result.memoryBeforeBytes = memoryBefore;
        result.memoryAfterBytes = memoryAfter;
        result.memoryDeltaBytes = memoryAfter - memoryBefore;
    }
    
    return result;
}

void BiopolymerPerformanceBenchmark::warmup(int iterations) {
    for (int i = 0; i < iterations; ++i) {
        // Do some work to warm up the system
        QElapsedTimer timer;
        timer.start();
        volatile int sum = 0;
        for (int j = 0; j < 1000000; ++j) {
            sum += j;
        }
        Q_UNUSED(sum);
        timer.invalidate();
    }
}

BenchmarkResult BiopolymerPerformanceBenchmark::calculateStatistics(
    const QList<BenchmarkResult>& results) const {
    
    if (results.isEmpty()) {
        return BenchmarkResult();
    }
    
    BenchmarkResult stats = results.first();
    
    qint64 totalMs = 0;
    qint64 minMs = results.first().elapsedMilliseconds;
    qint64 maxMs = results.first().elapsedMilliseconds;
    
    for (const BenchmarkResult& result : results) {
        totalMs += result.elapsedMilliseconds;
        minMs = qMin(minMs, result.elapsedMilliseconds);
        maxMs = qMax(maxMs, result.elapsedMilliseconds);
    }
    
    stats.averageMs = static_cast<double>(totalMs) / results.size();
    stats.minMs = static_cast<double>(minMs);
    stats.maxMs = static_cast<double>(maxMs);
    
    // Calculate standard deviation
    double sumSquaredDiffs = 0;
    for (const BenchmarkResult& result : results) {
        double diff = result.elapsedMilliseconds - stats.averageMs;
        sumSquaredDiffs += diff * diff;
    }
    stats.stdDevMs = qSqrt(sumSquaredDiffs / results.size());
    
    return stats;
}

bool BiopolymerPerformanceBenchmark::checkSystemResponsiveness() const {
    QElapsedTimer timer;
    timer.start();
    
    // Perform a small computation
    volatile int sum = 0;
    for (int i = 0; i < 100000; ++i) {
        sum += i;
    }
    Q_UNUSED(sum);
    
    timer.invalidate();
    
    // System should respond within 100ms for a simple computation
    return timer.msecsElapsed() < 100;
}

BiopolymerState BiopolymerPerformanceBenchmark::createTestState(int size) {
    return generateTestDocument(size);
}

QString BiopolymerPerformanceBenchmark::analyzeResults(const BenchmarkSuiteResult& result) {
    QStringList analysis;
    
    analysis << "=== Performance Analysis ===\n";
    
    // Check selection feedback
    for (int size : result.resultsBySize.keys()) {
        for (const BenchmarkResult& res : result.resultsBySize[size]) {
            if (res.operation == "Select residue") {
                if (res.averageMs > m_config.targetSelectionFeedbackMs) {
                    analysis << QString("Selection feedback for %1 residues is too slow: %2 ms (target: %3 ms)")
                        .arg(size).arg(res.averageMs).arg(m_config.targetSelectionFeedbackMs);
                } else {
                    analysis << QString("Selection feedback for %1 residues is good: %2 ms (< %3 ms)")
                        .arg(size).arg(res.averageMs).arg(m_config.targetSelectionFeedbackMs);
                }
            }
        }
    }
    
    // Check frame budget
    for (int size : result.resultsBySize.keys()) {
        for (const BenchmarkResult& res : result.resultsBySize[size]) {
            if (res.operation == "Render to canvas") {
                if (res.averageMs > m_config.targetFrameBudgetMs) {
                    analysis << QString("Rendering for %1 residues exceeds frame budget: %2 ms (target: %3 ms)")
                        .arg(size).arg(res.averageMs).arg(m_config.targetFrameBudgetMs);
                } else {
                    analysis << QString("Rendering for %1 residues is within budget: %2 ms (< %3 ms)")
                        .arg(size).arg(res.averageMs).arg(m_config.targetFrameBudgetMs);
                }
            }
        }
    }
    
    return analysis.join("\n");
}

QStringList BiopolymerPerformanceBenchmark::getRecommendations(
    const BenchmarkSuiteResult& result) {
    
    QStringList recommendations;
    
    for (int size : result.resultsBySize.keys()) {
        for (const BenchmarkResult& res : result.resultsBySize[size]) {
            if (res.operation == "Select residue" && res.averageMs > m_config.targetSelectionFeedbackMs) {
                recommendations << QString("Optimize selection feedback for large documents (%1 residues)").arg(size);
            }
            if (res.operation == "Render to canvas" && res.averageMs > m_config.targetFrameBudgetMs * 2) {
                recommendations << QString("Implement lazy rendering or view culling for %1+ residue documents").arg(size);
            }
        }
    }
    
    return recommendations;
}

// ========== Static Helpers ==========

QString BiopolymerPerformanceBenchmark::formatTimeMs(qint64 ms) {
    if (ms < 1) {
        return QString("%1 us").arg(ms * 1000);
    } else if (ms < 1000) {
        return QString("%1 ms").arg(ms);
    } else {
        return QString("%1 s").arg(ms / 1000.0, 0, 'f', 2);
    }
}

QString BiopolymerPerformanceBenchmark::formatMemorySize(qint64 bytes) {
    const qint64 kb = 1024;
    const qint64 mb = 1024 * kb;
    const qint64 gb = 1024 * mb;
    
    if (bytes < kb) {
        return QString("%1 B").arg(bytes);
    } else if (bytes < mb) {
        return QString("%1 KB").arg(static_cast<double>(bytes) / kb, 0, 'f', 2);
    } else if (bytes < gb) {
        return QString("%1 MB").arg(static_cast<double>(bytes) / mb, 0, 'f', 2);
    } else {
        return QString("%1 GB").arg(static_cast<double>(bytes) / gb, 0, 'f', 2);
    }
}

// ========== String Conversions ==========

QString benchmarkFormatToString(VectorExportFormat format) {
    return vectorExportFormatToString(static_cast<VectorExportFormat>(format));
}
