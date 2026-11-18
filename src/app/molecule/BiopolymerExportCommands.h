// src/app/molecule/BiopolymerExportCommands.h
#ifndef BIOPOLYMEREXPORTCOMMANDS_H
#define BIOPOLYMEREXPORTCOMMANDS_H

#include "BiopolymerExportSession.h"
#include <QString>
#include <QObject>

// Forward declaration
class DocumentState;

// ========== Export Command Base Class ==========

// Base class for biopolymer export commands
//
// Export commands are special - they never modify the document.
// They only read from the document and produce output.
// This is different from edit commands which modify the document.
//
// Export commands can be:
// 1. Synchronous - for small documents
// 2. Asynchronous - for large documents (via worker threads)
//
// Design principle: Export commands must be read-only operations.
class BiopolymerExportCommand {
public:
    virtual ~BiopolymerExportCommand();
    
    // Execute the export
    virtual ExportResult execute(DocumentState& document) const = 0;
    
    // Get command name
    virtual QString name() const = 0;
    
    // Get export format
    virtual ExportFormat format() const = 0;
    
    // Get options
    virtual ExportOptions options() const = 0;
    
    // Check if command can be executed
    virtual bool canExecute(const DocumentState& document) const = 0;
    
    // Get validation result
    virtual ExportValidation validate(const DocumentState& document) const = 0;
    
    // Check if export requires user confirmation
    virtual bool requiresConfirmation() const;
    
    // Clone the command
    virtual BiopolymerExportCommand* clone() const = 0;
};

// ========== Concrete Export Commands ==========

// Command to export as Sequence
class ExportSequenceCommand : public BiopolymerExportCommand {
public:
    explicit ExportSequenceCommand(const ExportOptions& options = ExportOptions());
    ~ExportSequenceCommand() override;
    
    ExportResult execute(DocumentState& document) const override;
    QString name() const override;
    ExportFormat format() const override;
    ExportOptions options() const override;
    bool canExecute(const DocumentState& document) const override;
    ExportValidation validate(const DocumentState& document) const override;
    bool requiresConfirmation() const override;
    ExportSequenceCommand* clone() const override;
    
    // Set options
    void setOptions(const ExportOptions& options);
    
private:
    ExportOptions m_options;
};

// Command to export as FASTA
class ExportFastaCommand : public BiopolymerExportCommand {
public:
    explicit ExportFastaCommand(const ExportOptions& options = ExportOptions());
    ~ExportFastaCommand() override;
    
    ExportResult execute(DocumentState& document) const override;
    QString name() const override;
    ExportFormat format() const override;
    ExportOptions options() const override;
    bool canExecute(const DocumentState& document) const override;
    ExportValidation validate(const DocumentState& document) const override;
    bool requiresConfirmation() const override;
    ExportFastaCommand* clone() const override;
    
    // Set options
    void setOptions(const ExportOptions& options);
    
private:
    ExportOptions m_options;
};

// Command to export as HELM
class ExportHelmCommand : public BiopolymerExportCommand {
public:
    explicit ExportHelmCommand(const ExportOptions& options = ExportOptions());
    ~ExportHelmCommand() override;
    
    ExportResult execute(DocumentState& document) const override;
    QString name() const override;
    ExportFormat format() const override;
    ExportOptions options() const override;
    bool canExecute(const DocumentState& document) const override;
    ExportValidation validate(const DocumentState& document) const override;
    bool requiresConfirmation() const override;
    ExportHelmCommand* clone() const override;
    
    // Set options
    void setOptions(const ExportOptions& options);
    
private:
    ExportOptions m_options;
};

// Command to export as IDT
class ExportIdtCommand : public BiopolymerExportCommand {
public:
    explicit ExportIdtCommand(const ExportOptions& options = ExportOptions());
    ~ExportIdtCommand() override;
    
    ExportResult execute(DocumentState& document) const override;
    QString name() const override;
    ExportFormat format() const override;
    ExportOptions options() const override;
    bool canExecute(const DocumentState& document) const override;
    ExportValidation validate(const DocumentState& document) const override;
    bool requiresConfirmation() const override;
    ExportIdtCommand* clone() const override;
    
    // Set options
    void setOptions(const ExportOptions& options);
    
private:
    ExportOptions m_options;
};

// Command to export as AxoLabs
class ExportAxoLabsCommand : public BiopolymerExportCommand {
public:
    explicit ExportAxoLabsCommand(const ExportOptions& options = ExportOptions());
    ~ExportAxoLabsCommand() override;
    
    ExportResult execute(DocumentState& document) const override;
    QString name() const override;
    ExportFormat format() const override;
    ExportOptions options() const override;
    bool canExecute(const DocumentState& document) const override;
    ExportValidation validate(const DocumentState& document) const override;
    bool requiresConfirmation() const override;
    ExportAxoLabsCommand* clone() const override;
    
    // Set options
    void setOptions(const ExportOptions& options);
    
private:
    ExportOptions m_options;
};

// ========== Export Command Factory ==========

// Factory for creating export commands
class BiopolymerExportCommandFactory {
public:
    // Create an export command for a specific format
    static BiopolymerExportCommand* createCommand(ExportFormat format, 
                                                  const ExportOptions& options = ExportOptions());
    
    // Create export command from format string
    static BiopolymerExportCommand* createCommand(const QString& formatString,
                                                  const ExportOptions& options = ExportOptions());
    
    // Get all available export formats
    static QList<ExportFormat> availableFormats();
    
    // Get display name for a format
    static QString formatDisplayName(ExportFormat format);
    
    // Get description for a format
    static QString formatDescription(ExportFormat format);
};

// ========== Export Manager ==========

// Manages export operations and provides a high-level interface
//
// This class provides a simplified interface for exporting biopolymers
// without needing to create command objects directly.
class BiopolymerExportManager {
public:
    // Constructor
    explicit BiopolymerExportManager(DocumentState& document);
    
    // Export to a specific format
    ExportResult exportTo(ExportFormat format, const ExportOptions& options = ExportOptions());
    
    // Export to Sequence
    ExportResult exportToSequence(const ExportOptions& options = ExportOptions());
    
    // Export to FASTA
    ExportResult exportToFasta(const ExportOptions& options = ExportOptions());
    
    // Export to HELM
    ExportResult exportToHelm(const ExportOptions& options = ExportOptions());
    
    // Export to IDT
    ExportResult exportToIdt(const ExportOptions& options = ExportOptions());
    
    // Export to AxoLabs
    ExportResult exportToAxoLabs(const ExportOptions& options = ExportOptions());
    
    // Validate export for a format
    ExportValidation validateExport(ExportFormat format) const;
    
    // Check if export is supported for a format
    bool isExportSupported(ExportFormat format) const;
    
    // Get supported formats
    QList<ExportFormat> supportedFormats() const;
    
    // Get capability assessment
    QMap<ExportFormat, ExportCapability> assessCapabilities() const;
    
    // Get default options for a format
    static ExportOptions defaultOptions(ExportFormat format);
    
private:
    DocumentState& m_document;
    
    // Create export session
    BiopolymerExportSession createSession() const;
};

#endif // BIOPOLYMEREXPORTCOMMANDS_H
