// src/app/molecule/BiopolymerImportCommands.h
#ifndef BIOPOLYMERIMPORTCOMMANDS_H
#define BIOPOLYMERIMPORTCOMMANDS_H

// Biopolymer import commands for transactional import operations.
//
// These commands wrap import operations to make them undoable.
// Import can either create a new document or replace the current one.

#include <functional>
#include <QList>
#include <QMap>
#include <QSet>
#include <QString>
#include "EditCommand.h"
#include "BiopolymerImportSession.h"

using AtomId = int;
using BondId = int;

// Forward declarations
class EditableMolecule;
class BiopolymerState;
class BiopolymerImportSession;

// ========== Import Command Types ==========

// Base class for biopolymer import commands
class BiopolymerImportCommand {
public:
    virtual ~BiopolymerImportCommand();
    
    // Execute the import
    virtual void perform(EditableMolecule& molecule, BiopolymerState& state) = 0;
    
    // Undo the import (restore previous state)
    virtual void revert(EditableMolecule& molecule, BiopolymerState& state) = 0;
    
    // Get description
    virtual QString description() const = 0;
};

// ========== Open as New Document ==========

// Command to import biopolymer data as a new document
// This clears the current document and creates the biopolymer
class ImportAsNewDocumentCommand : public BiopolymerImportCommand {
public:
    // Source data for the import
    QString sourceText;
    QString format;
    QString sourceFormat;
    
    // Previous state (for undo)
    EditableMolecule* previousMolecule = nullptr;
    BiopolymerState* previousBiopolymerState = nullptr;
    bool owned = false;
    
    ImportAsNewDocumentCommand(const QString& sourceText, const QString& format,
                              const QString& sourceFormat);
    
    ~ImportAsNewDocumentCommand() override;
    
    void perform(EditableMolecule& molecule, BiopolymerState& state) override;
    void revert(EditableMolecule& molecule, BiopolymerState& state) override;
    QString description() const override;
    
private:
    void capturePreviousState(EditableMolecule& molecule, BiopolymerState& state);
};

// ========== Replace Current Document ==========

// Command to import biopolymer data replacing the current document
// Similar to ImportAsNewDocument but preserves more metadata
class ReplaceDocumentCommand : public BiopolymerImportCommand {
public:
    QString sourceText;
    QString format;
    QString sourceFormat;
    
    // Previous state
    EditableMolecule* previousMolecule = nullptr;
    BiopolymerState* previousBiopolymerState = nullptr;
    bool owned = false;
    
    ReplaceDocumentCommand(const QString& sourceText, const QString& format,
                          const QString& sourceFormat);
    
    ~ReplaceDocumentCommand() override;
    
    void perform(EditableMolecule& molecule, BiopolymerState& state) override;
    void revert(EditableMolecule& molecule, BiopolymerState& state) override;
    QString description() const override;
    
private:
    void capturePreviousState(EditableMolecule& molecule, BiopolymerState& state);
};

// ========== Import Session Wrapper ==========

// High-level wrapper for managing import sessions and creating commands
class BiopolymerImportManager {
public:
    // Import mode
    enum class ImportMode {
        OpenAsNew,    // Create new document with imported data
        Replace,       // Replace current document with imported data
        Merge         // Merge with current document (future)
    };
    
    // Import result
    struct ImportResult {
        bool success = false;
        QString errorMessage;
        int chainCount = 0;
        int residueCount = 0;
        int atomCount = 0;
        int bondCount = 0;
        QList<QString> warnings;
        QList<QString> validationIssues;
    };
    
    // Constructor
    BiopolymerImportManager();
    
    // Import biopolymer from text
    ImportResult importFromText(EditableMolecule& molecule, BiopolymerState& bioState,
                                const QString& text, const QString& format = QString(),
                                ImportMode mode = ImportMode::OpenAsNew);
    
    // Import biopolymer from text with options
    ImportResult importFromText(EditableMolecule& molecule, BiopolymerState& bioState,
                                const QString& text, const ImportOptions& options,
                                const QString& format = QString(), ImportMode mode = ImportMode::OpenAsNew);
    
    // Create import command (for undoable import)
    BiopolymerImportCommand* createImportCommand(const QString& sourceText,
                                                 const QString& format,
                                                 ImportMode mode);
    
    // Validate input before import
    static ParseResult validateInput(const QString& text,
                                                             const QString& format = QString());
    
    // Check if input can be imported
    static bool canImport(const QString& text, const QString& format = QString());
    
    // Get supported formats
    static QList<QString> supportedFormats();
    
    // Detect format from text
    static QString detectFormat(const QString& text);
    
    // ========== Import with Draft ==========
    
    // Start a new import session (for UI with preview)
    BiopolymerImportSession* startImportSession(const QString& text, 
                                               const ImportOptions& options = {});
    
    // Continue an existing session
    BiopolymerImportSession* continueImportSession(BiopolymerImportSession* session);
    
    // Commit import session
    ImportResult commitImportSession(BiopolymerImportSession* session,
                                    EditableMolecule& molecule, BiopolymerState& bioState,
                                    ImportMode mode);
    
    // Abort import session
    void abortImportSession(BiopolymerImportSession* session);
    
    // ========== Batch Import ==========
    
    // Import multiple sequences
    struct BatchImportResult {
        QList<ImportResult> results;
        int successCount = 0;
        int failureCount = 0;
    };
    
    BatchImportResult batchImport(const QList<QString>& texts, EditableMolecule& molecule,
                                  BiopolymerState& bioState, ImportMode mode = ImportMode::OpenAsNew);
    
private:
    // Current import session (for UI)
    BiopolymerImportSession* m_currentSession = nullptr;
    
    // Helper methods
    ImportResult doImport(EditableMolecule& molecule, BiopolymerState& bioState,
                         BiopolymerImportSession& session, ImportMode mode);
    
    void restorePreviousState(EditableMolecule& molecule, BiopolymerState& bioState,
                              EditableMolecule* previousMolecule, 
                              BiopolymerState* previousBiopolymerState) const;
};

#endif // BIOPOLYMERIMPORTCOMMANDS_H
