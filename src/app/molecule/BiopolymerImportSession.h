// src/app/molecule/BiopolymerImportSession.h
#ifndef BIOPOLYMERIMPORTSESSION_H
#define BIOPOLYMERIMPORTSESSION_H

// Biopolymer import session for draft management and format-specific parsing.
//
// Manages the lifecycle of importing biopolymer data:
// - Separates input draft from generated output
// - Preserves drafts when changing formats
// - Provides format-specific parsing (FASTA, Sequence, HELM, IDT, etc.)
// - Detects unsupported features with precise location information
// - Validates input before chemical commit
// - Provides non-atomistic preview for ambiguous input
//
// Design principle: Import is a multi-step process:
// 1. Parse raw text into a draft structure
// 2. Validate draft and detect unsupported features
// 3. User reviews warnings and resolves ambiguities
// 4. Commit to actual biopolymer state (creates atoms/bonds)
//
// The draft is preserved across format changes and until explicit commit or abort.

#include <QString>
#include <QList>
#include <QVector>
#include <QHash>
#include <QMap>
#include <QSet>
#include <QDateTime>
#include "BiopolymerState.h"

// Forward declarations
class EditableMolecule;
class BiopolymerState;

// ========== Import Error and Warning Types ==========

// Severity levels for import issues
enum class ImportIssueSeverity {
    Error,       // Blocks commit, must be resolved
    Warning,     // Allows commit but user should be aware
    Info         // Informational only
};

// Location information for error reporting
struct SourceLocation {
    int line = -1;         // Line number (1-indexed)
    int column = -1;       // Column number (1-indexed)
    int recordIndex = -1;  // Record index for multi-record formats (FASTA)
    QString recordId;      // Record identifier (FASTA header, etc.)
    
    // Check if location is valid
    bool isValid() const { return line >= 1; }
    
    // String representation
    QString toString() const;
};

// Single import issue
struct ImportIssue {
    ImportIssueSeverity severity = ImportIssueSeverity::Error;
    QString code;           // Machine-readable code (e.g., "UNSUPPORTED_SYMBOL")
    QString message;        // Human-readable message
    SourceLocation location;
    QString suggestion;     // Suggested resolution
    
    // Check if this issue blocks commit
    bool blocksCommit() const { return severity == ImportIssueSeverity::Error; }
    
    // String representation
    QString toString() const;
};

// ========== Parse Result ==========

// Result of parsing a biopolymer input
struct ParseResult {
    bool success = false;
    QString errorMessage;
    
    // Parsed chains
    QList<ParsedChain> chains;
    
    // List of all issues found during parsing
    QList<ImportIssue> issues;
    
    // Unsupported features detected
    QSet<QString> unsupportedFeatures;
    
    // Whether the input contains ambiguous symbols
    bool hasAmbiguousSymbols = false;
    
    // Whether the input contains unsupported chemistry
    bool hasUnsupportedChemistry = false;
    
    // Check if parse was successful enough to attempt commit
    bool canAttemptCommit() const;
    
    // Get all blocking issues
    QList<ImportIssue> blockingIssues() const;
    
    // Get all warnings
    QList<ImportIssue> warnings() const;
};

// ========== Parsed Structures ==========

// Parsed monomer identity (from input text)
struct ParsedMonomer {
    QString symbol;              // Single-letter code from input
    QString originalSymbol;      // Original symbol (may be ambiguous)
    int position = -1;          // Position in chain (1-indexed)
    SourceLocation location;    // Source location of this monomer
    bool isAmbiguous = false;   // True for B, Z, X, N, etc.
    bool isSupported = true;    // True if we can map to a known monomer
    QString resolvedSymbol;     // Resolved symbol after disambiguation
};

// Parsed chain (from input text)
struct ParsedChain {
    QString name;                // Chain name (from FASTA header or default)
    QString originalName;        // Original chain name from input
    PolymerType polymerType = PolymerType::Unknown;
    QString polymerTypeStr;     // String representation of polymer type
    
    QList<ParsedMonomer> monomers;
    
    SourceLocation headerLocation;  // Location of chain header
    SourceLocation bodyLocation;    // Location of chain body
    
    int recordIndex = -1;       // Index in multi-record input
    QString recordId;            // Record identifier (FASTA defline)
    QString description;        // Description from header
    
    // Validation state
    bool isValid = true;
    QString invalidityReason;
};

// ========== Input Format Types ==========

// Supported input format types
enum class ImportFormat {
    Unknown,
    FASTA,           // FASTA format (multi-sequence)
    Sequence,        // Plain sequence (single chain)
    HELM,            // HELM notation
    IDT,             // Integrated DNA Technologies format
    AxoLabs,         // AxoLabs format
    Simple,          // Simple format (comma-separated chains)
};

// Format detection result
struct FormatDetection {
    ImportFormat format = ImportFormat::Unknown;
    double confidence = 0.0;    // 0.0 to 1.0
    QString reason;             // Why this format was detected
};

// ========== Import Options ==========

// Options for controlling import behavior
struct ImportOptions {
    // Whether to automatically resolve ambiguous symbols
    // If false, ambiguous symbols will block commit
    bool autoResolveAmbiguous = false;
    
    // Default polymer type if not specified in input
    PolymerType defaultPolymerType = PolymerType::Peptide;
    
    // Whether to create multiple chains from multi-record input
    bool createMultipleChains = true;
    
    // Whether to preserve FASTA record boundaries
    bool preserveRecordBoundaries = true;
    
    // Whether to validate stereochemistry
    bool validateStereochemistry = true;
    
    // Whether to detect and report unsupported features
    bool detectUnsupportedFeatures = true;
    
    // Whether to attempt partial import on errors
    bool allowPartialImport = false;
};

// ========== Import Preview ==========

// Preview of what will be created on commit
struct ImportPreview {
    // Number of chains that will be created
    int chainCount = 0;
    
    // Number of residues that will be created
    int residueCount = 0;
    
    // Estimated number of atoms
    int estimatedAtomCount = 0;
    
    // Estimated number of bonds
    int estimatedBondCount = 0;
    
    // Chain information
    QList<PreviewChain> chains;
    
    // Whether commit is blocked by errors
    bool commitBlocked = false;
    
    // Issues that block commit
    QList<ImportIssue> blockingIssues;
    
    // All warnings
    QList<ImportIssue> warnings;
};

// Preview chain information
struct PreviewChain {
    QString name;
    PolymerType polymerType = PolymerType::Unknown;
    int residueCount = 0;
    int estimatedAtomCount = 0;
    bool hasErrors = false;
    bool hasWarnings = false;
};

// ========== Main Import Session Class ==========

// Manages a single import session with draft state
class BiopolymerImportSession {
public:
    BiopolymerImportSession();
    explicit BiopolymerImportSession(const QString& sourceText, const ImportOptions& options = ImportOptions());
    ~BiopolymerImportSession();
    
    // ========== Source Management ==========
    
    // Set the source text to import
    void setSourceText(const QString& text);
    
    // Get the current source text
    QString sourceText() const;
    
    // Clear the source text
    void clearSource();
    
    // ========== Format Detection ==========
    
    // Detect the format of the source text
    FormatDetection detectFormat() const;
    
    // Get the detected or specified format
    ImportFormat format() const;
    
    // Set the format explicitly (overrides auto-detection)
    void setFormat(ImportFormat format);
    
    // ========== Parsing ==========
    
    // Parse the source text into a draft
    ParseResult parse();
    
    // Parse with specific format
    ParseResult parseAs(ImportFormat format);
    
    // Re-parse after options change
    ParseResult reparse();
    
    // Check if parsing has been attempted
    bool isParsed() const;
    
    // ========== Draft Access ==========
    
    // Get all parsed chains
    QList<ParsedChain> parsedChains() const;
    
    // Get parsed chain by index
    ParsedChain* parsedChain(int index);
    const ParsedChain* parsedChain(int index) const;
    
    // Get parsed chain by name
    ParsedChain* parsedChainByName(const QString& name);
    const ParsedChain* parsedChainByName(const QString& name) const;
    
    // Get total monomer count
    int monomerCount() const;
    
    // ========== Validation ==========
    
    // Validate the current draft
    QList<ImportIssue> validate() const;
    
    // Validate with specific options
    QList<ImportIssue> validate(const ImportOptions& options) const;
    
    // Check if draft can be committed
    bool canCommit() const;
    
    // Get all issues
    QList<ImportIssue> allIssues() const;
    
    // Get issues for a specific chain
    QList<ImportIssue> issuesForChain(int chainIndex) const;
    
    // Get issues for a specific monomer
    QList<ImportIssue> issuesForMonomer(int chainIndex, int monomerIndex) const;
    
    // ========== Ambiguity Resolution ==========
    
    // Get all ambiguous symbols
    QList<ParsedMonomer*> ambiguousSymbols();
    
    // Resolve an ambiguous symbol to a specific monomer
    bool resolveAmbiguousSymbol(int chainIndex, int monomerIndex, const QString& resolvedSymbol);
    
    // Auto-resolve all ambiguous symbols using default rules
    int autoResolveAllAmbiguous();
    
    // ========== Preview ==========
    
    // Generate a preview of what will be created
    ImportPreview generatePreview() const;
    
    // Estimate atom/bond counts for a chain
    static PreviewChain estimateChainPreview(const ParsedChain& chain);
    
    // ========== Commit ==========
    
    // Commit the draft to a BiopolymerState
    // Returns true if commit succeeded, false if blocked by errors
    bool commitToState(BiopolymerState& state, EditableMolecule& molecule) const;
    
    // Commit with progress callback
    // The callback receives (currentStep, totalSteps, description)
    using ProgressCallback = std::function<void(int, int, const QString&)>;
    bool commitToState(BiopolymerState& state, EditableMolecule& molecule, 
                       const ProgressCallback& progress) const;
    
    // ========== Options ==========
    
    // Get current options
    ImportOptions options() const;
    
    // Set options
    void setOptions(const ImportOptions& options);
    
    // ========== State ==========
    
    // Check if session has valid parsed data
    bool hasValidData() const;
    
    // Check if session has errors
    bool hasErrors() const;
    
    // Check if session has warnings
    bool hasWarnings() const;
    
    // Check if session has unsupported features
    bool hasUnsupportedFeatures() const;
    
    // Get unsupported features
    QSet<QString> unsupportedFeatures() const;
    
    // ========== Reset ==========
    
    // Reset the session (clear everything)
    void reset();
    
    // ========== Format-Specific Parsing ==========
    
    // Parse FASTA format
    static ParseResult parseFasta(const QString& text, const ImportOptions& options);
    
    // Parse plain sequence format
    static ParseResult parseSequence(const QString& text, const ImportOptions& options);
    
    // Parse simple format (comma-separated chains)
    static ParseResult parseSimple(const QString& text, const ImportOptions& options);
    
    // ========== Helper Methods ==========
    
    // Check if a symbol is a valid amino acid
    static bool isValidAminoAcid(const QString& symbol);
    
    // Check if a symbol is a valid nucleotide
    static bool isValidNucleotide(const QString& symbol);
    
    // Check if a symbol is ambiguous
    static bool isAmbiguousSymbol(const QString& symbol);
    
    // Get all valid amino acid symbols
    static QSet<QString> validAminoAcidSymbols();
    
    // Get all valid nucleotide symbols
    static QSet<QString> validNucleotideSymbols();
    
    // Get all ambiguous symbols
    static QSet<QString> ambiguousSymbols();
    
    // Resolve ambiguous symbol with default rules
    static QString resolveAmbiguousSymbol(const QString& symbol, PolymerType polymerType);
    
    // Normalize sequence (remove whitespace, convert to uppercase)
    static QString normalizeSequence(const QString& sequence);
    
    // ========== Serialization Support ==========
    
    // Convert to JSON for serialization
    QMap<QString, QVariant> toJson() const;
    
    // Load from JSON
    static BiopolymerImportSession fromJson(const QMap<QString, QVariant>& json);
    
private:
    // Source text
    QString m_sourceText;
    
    // Parsed data
    QList<ParsedChain> m_parsedChains;
    QList<ImportIssue> m_issues;
    QSet<QString> m_unsupportedFeatures;
    
    // Format
    ImportFormat m_format = ImportFormat::Unknown;
    bool m_formatExplicit = false;
    
    // Options
    ImportOptions m_options;
    
    // Parsing state
    bool m_isParsed = false;
    bool m_parseSuccess = false;
    
    // Helper methods
    FormatDetection detectFormatInternal() const;
    
    // Parse FASTA internally
    ParseResult parseFastaInternal(const QString& text);
    
    // Parse sequence internally
    ParseResult parseSequenceInternal(const QString& text);
    
    // Parse simple format internally
    ParseResult parseSimpleInternal(const QString& text);
    
    // Validate a parsed chain
    QList<ImportIssue> validateChain(const ParsedChain& chain, int chainIndex) const;
    
    // Validate a monomer
    QList<ImportIssue> validateMonomer(const ParsedMonomer& monomer, 
                                       int chainIndex, int monomerIndex) const;
    
    // Commit a single chain to state
    int commitChain(const ParsedChain& chain, BiopolymerState& state, 
                    EditableMolecule& molecule, const ProgressCallback& progress,
                    int& currentStep) const;
    
    // Create atoms and bonds for a monomer
    QList<AtomId> createMonomerAtoms(int residueId, const ParsedMonomer& monomer,
                                    EditableMolecule& molecule, double x, double y) const;
    
    // Connect two residues
    void connectResidues(int fromResidueId, int toResidueId, BondId bondId,
                         BiopolymerState& state, EditableMolecule& molecule) const;
    
    // Calculate position for a new monomer
    QPointF calculateMonomerPosition(int chainIndex, int monomerIndex, 
                                     const QList<QPointF>& existingPositions) const;
};

// ========== Format Registry ==========

// Registry of format parsers for extensibility
class ImportFormatRegistry {
public:
    // Register a format parser
    using FormatParser = std::function<ParseResult(const QString&, const ImportOptions&)>;
    
    static void registerParser(ImportFormat format, const QString& name, 
                               const FormatParser& parser);
    
    // Get parser for a format
    static FormatParser getParser(ImportFormat format);
    
    // Get all registered formats
    static QList<ImportFormat> allFormats();
    
    // Get format name
    static QString formatName(ImportFormat format);
    
    // Get format by name
    static ImportFormat formatByName(const QString& name);
    
private:
    struct FormatInfo {
        QString name;
        FormatParser parser;
    };
    
    static QMap<ImportFormat, FormatInfo>& registry();
};

#endif // BIOPOLYMERIMPORTSESSION_H
