// src/app/molecule/BiopolymerExportSession.h
#ifndef BIOPOLYMEREXPORTSESSION_H
#define BIOPOLYMEREXPORTSESSION_H

#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include <QString>
#include <QList>
#include <QMap>
#include <QVariant>
#include <QSet>

// ========== Export Format Types ==========

// Supported export formats for biopolymer data
enum class ExportFormat {
    Sequence,      // Plain sequence (single chain only)
    Fasta,         // FASTA format (multi-chain supported)
    Helm,          // HELM notation (Hierarchical Editing Language for Macromolecules)
    Idt,           // IDT format (Integrated DNA Technologies)
    AxoLabs,       // AxoLabs format
    Unknown        // Unknown/unsupported format
};

// Convert export format to string
QString exportFormatToString(ExportFormat format);

// Parse export format from string (case-insensitive)
ExportFormat exportFormatFromString(const QString& formatString);

// List of all supported export formats
QList<ExportFormat> allExportFormats();

// ========== Export Capability ==========

// Capability assessment for a biopolymer export
struct ExportCapability {
    ExportFormat format = ExportFormat::Unknown;
    bool isSupported = false;
    bool requiresExplicitConversion = false;
    QString statusMessage;  // "Success", "ConversionRequired", or "Unsupported"
    
    // Diagnostics for what's blocked or changed
    QList<QString> blockedReasons;
    QList<QString> changedFeatures;
    QList<QString> omittedFeatures;
    
    // Source document information
    int documentId = -1;
    quint64 documentRevision = 0;
    QString documentName;
    
    // Content assessment
    int chainCount = 0;
    int totalResidueCount = 0;
    QSet<PolymerType> polymerTypes;
    bool hasValidMappings = false;
    bool hasInvalidMappings = false;
    bool hasUnsupportedModifications = false;
    bool hasUnsupportedConnections = false;
    bool hasMultipleChains = false;
    bool hasTerminalGroups = false;
    bool hasBranchedTopology = false;
    
    // Check if export can proceed
    bool canExport() const;
    
    // Check if export requires user confirmation
    bool needsUserConfirmation() const;
    
    // Get formatted diagnostic message
    QString diagnosticMessage() const;
    
    // Get formatted change summary
    QString changeSummary() const;
    
    // Check if all mappings are valid
    bool allMappingsValid() const;
    
    // Check if format supports multiple chains
    static bool formatSupportsMultipleChains(ExportFormat format);
    
    // Check if format requires explicit conversion
    static bool formatRequiresConversion(ExportFormat format, const BiopolymerState& state);
};

// ========== Export Result ==========

// Result of a biopolymer export operation
struct ExportResult {
    ExportFormat format = ExportFormat::Unknown;
    bool success = false;
    QString outputText;
    QString errorMessage;
    
    // Source information
    int documentId = -1;
    quint64 documentRevision = 0;
    
    // Content information
    int chainCount = 0;
    int residueCount = 0;
    
    // Validation information
    bool mappingsValid = false;
    QList<QString> warnings;
    QList<QString> validationIssues;
    
    // Check if export succeeded
    bool isSuccess() const;
    
    // Check if export failed
    bool isFailure() const;
    
    // Get formatted output with warnings
    QString formattedOutput() const;
};

// ========== Export Options ==========

// Options for biopolymer export
struct ExportOptions {
    ExportFormat format = ExportFormat::Fasta;
    
    // FASTA-specific options
    QString fastaTitle;
    bool fastaIncludeDescription = true;
    int fastaLineWrap = 60;
    
    // Sequence-specific options
    bool sequenceSingleChainOnly = true;  // Block if multiple chains
    bool sequenceWarnOnModifications = true;
    
    // HELM-specific options
    bool helmIncludePolymerSection = true;
    bool helmIncludeAnnotationSection = true;
    
    // General options
    bool includeChainNames = true;
    bool includeResidueNumbers = false;
    bool includeTerminalGroups = true;
    bool forceSimplification = false;  // Require explicit user confirmation
    
    // Serialization
    QMap<QString, QVariant> toJson() const;
    static ExportOptions fromJson(const QMap<QString, QVariant>& json);
};

// ========== Export Diagnostics ==========

// Diagnostic information for export validation
struct ExportDiagnostic {
    QString category;  // "Mapping", "Chemistry", "Topology", "Format"
    QString severity;  // "Error", "Warning", "Info"
    QString message;
    int chainId = -1;
    int residueId = -1;
    QString feature;
    
    // Check if diagnostic is an error
    bool isError() const;
    
    // Check if diagnostic is a warning
    bool isWarning() const;
    
    // Get formatted message
    QString formattedMessage() const;
};

// ========== Export Validation ==========

// Validation result for export pre-check
struct ExportValidation {
    ExportFormat format = ExportFormat::Unknown;
    bool isValid = false;
    QList<ExportDiagnostic> diagnostics;
    
    // Quick checks
    bool hasErrors() const;
    bool hasWarnings() const;
    bool hasBlockers() const;
    
    // Get all error messages
    QList<QString> errorMessages() const;
    
    // Get all warning messages
    QList<QString> warningMessages() const;
    
    // Merge with another validation
    void merge(const ExportValidation& other);
};

// ========== BiopolymerExportSession ==========

// Session for managing biopolymer export operations
//
// This class manages the export workflow for biopolymers, including:
// - Capability assessment for different export formats
// - Validation of document state before export
// - Generation of export output
// - Management of export drafts and results
//
// Design principle: Export never modifies the document. All export operations
// are read-only and produce isolated results.
class BiopolymerExportSession {
public:
    // Request ID type for tracking export operations
    using RequestId = int;
    
    // Constructor
    explicit BiopolymerExportSession(const BiopolymerState& state, 
                                      const BiopolymerWorkspace& workspace);
    
    // Initialize session with document information
    void initialize(int documentId, quint64 documentRevision, const QString& documentName);
    
    // ========== Capability Assessment ==========
    
    // Assess export capability for a specific format
    ExportCapability assessCapability(ExportFormat format) const;
    
    // Assess capability for all formats
    QMap<ExportFormat, ExportCapability> assessAllCapabilities() const;
    
    // Get supported formats for current state
    QList<ExportFormat> supportedFormats() const;
    
    // Check if a format is supported
    bool isFormatSupported(ExportFormat format) const;
    
    // ========== Validation ==========
    
    // Validate document state for export
    ExportValidation validateForExport(ExportFormat format) const;
    
    // Validate all chains have valid mappings
    bool validateAllMappings() const;
    
    // Validate no unsupported modifications
    bool validateNoUnsupportedModifications() const;
    
    // Validate chain count for format
    bool validateChainCount(ExportFormat format) const;
    
    // ========== Export Generation ==========
    
    // Generate export output for a specific format
    ExportResult generateExport(ExportFormat format, const ExportOptions& options = ExportOptions()) const;
    
    // Generate sequence export (plain text)
    ExportResult generateSequenceExport(const ExportOptions& options = ExportOptions()) const;
    
    // Generate FASTA export
    ExportResult generateFastaExport(const ExportOptions& options = ExportOptions()) const;
    
    // Generate HELM export
    ExportResult generateHelmExport(const ExportOptions& options = ExportOptions()) const;
    
    // Generate IDT export
    ExportResult generateIdtExport(const ExportOptions& options = ExportOptions()) const;
    
    // Generate AxoLabs export
    ExportResult generateAxoLabsExport(const ExportOptions& options = ExportOptions()) const;
    
    // ========== Format-Specific Implementation ==========
    
    // Build sequence string for a chain
    QString buildSequenceString(int chainId, bool includeNumbers = false) const;
    
    // Build sequence string for all chains
    QMap<int, QString> buildAllSequenceStrings(bool includeNumbers = false) const;
    
    // Build FASTA record for a chain
    QString buildFastaRecord(int chainId, const QString& title, int lineWrap = 60) const;
    
    // Build FASTA for all chains
    QString buildFastaForAllChains(const QString& baseTitle, int lineWrap = 60) const;
    
    // Build HELM notation for the biopolymer
    QString buildHelmNotation() const;
    
    // Build IDT sequence
    QString buildIdtSequence() const;
    
    // Build AxoLabs format
    QString buildAxoLabsFormat() const;
    
    // ========== Utility Methods ==========
    
    // Get document information
    int documentId() const;
    quint64 documentRevision() const;
    QString documentName() const;
    
    // Get biopolymer state
    const BiopolymerState& state() const;
    const BiopolymerWorkspace& workspace() const;
    
    // Get summary information
    int chainCount() const;
    int totalResidueCount() const;
    QList<PolymerType> polymerTypes() const;
    bool hasValidMappings() const;
    bool hasInvalidMappings() const;
    bool hasMultipleChains() const;
    
    // Check if state supports a polymer type
    bool supportsPolymerType(PolymerType type) const;
    
    // Check if all residues in a chain are standard
    bool allResiduesStandard(int chainId) const;
    
    // Get unsupported modifications
    QList<QString> unsupportedModifications() const;
    
    // Get unsupported connections
    QList<QString> unsupportedConnections() const;
    
    // ========== Session State ==========
    
    // Reset session
    void reset();
    
    // Clear drafts
    void clearDrafts();
    
    // Get last export result
    ExportResult lastResult() const;
    
    // ========== Static Helpers ==========
    
    // Check if format supports polymer type
    static bool formatSupportsPolymerType(ExportFormat format, PolymerType polymerType);
    
    // Get format file extension
    static QString formatFileExtension(ExportFormat format);
    
    // Get format MIME type
    static QString formatMimeType(ExportFormat format);
    
    // Get format display name
    static QString formatDisplayName(ExportFormat format);
    
    // Get format description
    static QString formatDescription(ExportFormat format);
    
    // Check if format is text-based
    static bool formatIsTextBased(ExportFormat format);
    
    // Get default options for format
    static ExportOptions defaultOptions(ExportFormat format);

private:
    // ========== Private Members ==========
    
    const BiopolymerState& m_state;
    const BiopolymerWorkspace& m_workspace;
    
    int m_documentId = -1;
    quint64 m_documentRevision = 0;
    QString m_documentName;
    
    mutable ExportResult m_lastResult;
    
    // ========== Private Helper Methods ==========
    
    // Validate residue identity
    bool validateResidueIdentity(const BiopolymerResidue& residue) const;
    
    // Check if chain has valid sequence
    bool chainHasValidSequence(int chainId) const;
    
    // Get chain polymer type
    PolymerType chainPolymerType(int chainId) const;
    
    // Check if format supports the current polymer types
    bool formatSupportsCurrentPolymerTypes(ExportFormat format) const;
    
    // Build chain name for export
    QString buildChainName(int chainId) const;
    
    // Build residue identity string
    QString buildResidueIdentityString(const BiopolymerResidue& residue) const;
    
    // Check if simplification is required
    bool requiresSimplification(ExportFormat format) const;
    
    // Get required simplification type
    QString requiredSimplificationType(ExportFormat format) const;
};

#endif // BIOPOLYMEREXPORTSESSION_H
