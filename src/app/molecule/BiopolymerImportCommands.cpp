// src/app/molecule/BiopolymerImportCommands.cpp
#include "BiopolymerImportCommands.h"
#include "BiopolymerImportSession.h"
#include "BiopolymerState.h"
#include "ProjectFile.h"
#include <QDebug>

// ========== ImportAsNewDocumentCommand Implementation ==========

ImportAsNewDocumentCommand::ImportAsNewDocumentCommand(const QString& sourceText,
                                                     const QString& format,
                                                     const QString& sourceFormat)
    : sourceText(sourceText), format(format), sourceFormat(sourceFormat) {
}

ImportAsNewDocumentCommand::~ImportAsNewDocumentCommand() {
    if (owned) {
        delete previousMolecule;
        delete previousBiopolymerState;
    }
}

void ImportAsNewDocumentCommand::perform(EditableMolecule& molecule) {
    capturePreviousState(molecule);
    
    // Clear current molecule
    molecule.clear();
    
    // Import biopolymer
    BiopolymerImportSession session(sourceText);
    session.setOptions(BiopolymerImportSession::ImportOptions());
    
    // Auto-detect format if not specified
    if (format.isEmpty()) {
        session.setFormat(BiopolymerImportSession::ImportFormat::Sequence);
    }
    
    // Parse and commit
    session.parse();
    if (session.canCommit()) {
        session.commitToState(/* state */ molecule.biopolymerState(), molecule);
    }
}

void ImportAsNewDocumentCommand::revert(EditableMolecule& molecule) {
    if (previousMolecule) {
        // Restore from previous state
        molecule = *previousMolecule;
    }
}

QString ImportAsNewDocumentCommand::description() const {
    return QString("Import biopolymer as new document from %1").arg(sourceFormat.isEmpty() ? "text" : sourceFormat);
}

void ImportAsNewDocumentCommand::capturePreviousState(EditableMolecule& molecule) {
    if (previousMolecule) {
        if (owned) {
            delete previousMolecule;
        }
        owned = false;
    }
    
    previousMolecule = new EditableMolecule(molecule);
    owned = true;
}

// ========== ReplaceDocumentCommand Implementation ==========

ReplaceDocumentCommand::ReplaceDocumentCommand(const QString& sourceText,
                                             const QString& format,
                                             const QString& sourceFormat)
    : sourceText(sourceText), format(format), sourceFormat(sourceFormat) {
}

ReplaceDocumentCommand::~ReplaceDocumentCommand() {
    if (owned) {
        delete previousMolecule;
        delete previousBiopolymerState;
    }
}

void ReplaceDocumentCommand::perform(EditableMolecule& molecule) {
    capturePreviousState(molecule);
    
    // Clear current molecule
    molecule.clear();
    
    // Import biopolymer
    BiopolymerImportSession session(sourceText);
    session.setOptions(BiopolymerImportSession::ImportOptions());
    
    // Parse and commit
    session.parse();
    if (session.canCommit()) {
        session.commitToState(/* state */ molecule.biopolymerState(), molecule);
    }
}

void ReplaceDocumentCommand::revert(EditableMolecule& molecule) {
    if (previousMolecule) {
        molecule = *previousMolecule;
    }
}

QString ReplaceDocumentCommand::description() const {
    return QString("Replace document with imported biopolymer from %1").arg(sourceFormat.isEmpty() ? "text" : sourceFormat);
}

void ReplaceDocumentCommand::capturePreviousState(EditableMolecule& molecule) {
    if (previousMolecule) {
        if (owned) {
            delete previousMolecule;
        }
        owned = false;
    }
    
    previousMolecule = new EditableMolecule(molecule);
    owned = true;
}

// ========== BiopolymerImportManager Implementation ==========

BiopolymerImportManager::BiopolymerImportManager() {
}

BiopolymerImportSession::ParseResult BiopolymerImportManager::validateInput(const QString& text,
                                                                         const QString& format) {
    BiopolymerImportSession session;
    session.setSourceText(text);
    
    if (!format.isEmpty()) {
        ImportFormat fmt = ImportFormatRegistry::formatByName(format);
        if (fmt != ImportFormat::Unknown) {
            session.setFormat(fmt);
        }
    }
    
    return session.parse();
}

bool BiopolymerImportManager::canImport(const QString& text, const QString& format) {
    ParseResult result = validateInput(text, format);
    return result.success && result.canAttemptCommit();
}

QList<QString> BiopolymerImportManager::supportedFormats() {
    QList<QString> formats;
    for (ImportFormat fmt : ImportFormatRegistry::allFormats()) {
        formats.append(ImportFormatRegistry::formatName(fmt));
    }
    return formats;
}

QString BiopolymerImportManager::detectFormat(const QString& text) {
    BiopolymerImportSession session;
    session.setSourceText(text);
    FormatDetection detection = session.detectFormat();
    return ImportFormatRegistry::formatName(detection.format);
}

BiopolymerImportManager::ImportResult BiopolymerImportManager::importFromText(
    EditableMolecule& molecule, BiopolymerState& bioState,
    const QString& text, const QString& format, ImportMode mode) {
    
    BiopolymerImportSession::ImportOptions options;
    return importFromText(molecule, bioState, text, options, format, mode);
}

BiopolymerImportManager::ImportResult BiopolymerImportManager::importFromText(
    EditableMolecule& molecule, BiopolymerState& bioState,
    const QString& text, const BiopolymerImportSession::ImportOptions& options,
    const QString& format, ImportMode mode) {
    
    BiopolymerImportSession session(text, options);
    
    if (!format.isEmpty()) {
        ImportFormat fmt = ImportFormatRegistry::formatByName(format);
        if (fmt != ImportFormat::Unknown) {
            session.setFormat(fmt);
        }
    }
    
    return doImport(molecule, bioState, session, mode);
}

BiopolymerImportCommand* BiopolymerImportManager::createImportCommand(
    const QString& sourceText, const QString& format, ImportMode mode) {
    
    switch (mode) {
        case ImportMode::OpenAsNew:
        case ImportMode::Replace:
            return new ReplaceDocumentCommand(sourceText, format, 
                                              format.isEmpty() ? "text" : format);
        default:
            return new ImportAsNewDocumentCommand(sourceText, format,
                                                   format.isEmpty() ? "text" : format);
    }
}

BiopolymerImportSession* BiopolymerImportManager::startImportSession(
    const QString& text, const BiopolymerImportSession::ImportOptions& options) {
    
    // Clean up previous session
    abortImportSession(m_currentSession);
    
    m_currentSession = new BiopolymerImportSession(text, options);
    return m_currentSession;
}

BiopolymerImportSession* BiopolymerImportManager::continueImportSession(
    BiopolymerImportSession* session) {
    
    // Just return the session (no state change)
    return session;
}

BiopolymerImportManager::ImportResult BiopolymerImportManager::commitImportSession(
    BiopolymerImportSession* session, EditableMolecule& molecule, 
    BiopolymerState& bioState, ImportMode mode) {
    
    if (!session) {
        return ImportResult{false, "No active import session", 0, 0, 0, 0, {}, {}};
    }
    
    ImportResult result = doImport(molecule, bioState, *session, mode);
    
    // Clean up session
    abortImportSession(session);
    
    return result;
}

void BiopolymerImportManager::abortImportSession(BiopolymerImportSession* session) {
    if (session == m_currentSession) {
        m_currentSession = nullptr;
    }
    delete session;
}

BiopolymerImportManager::BatchImportResult BiopolymerImportManager::batchImport(
    const QList<QString>& texts, EditableMolecule& molecule,
    BiopolymerState& bioState, ImportMode mode) {
    
    BatchImportResult result;
    
    for (const QString& text : texts) {
        ImportResult importResult = importFromText(molecule, bioState, text, QString(), mode);
        result.results.append(importResult);
        
        if (importResult.success) {
            result.successCount++;
        } else {
            result.failureCount++;
        }
    }
    
    return result;
}

BiopolymerImportManager::ImportResult BiopolymerImportManager::doImport(
    EditableMolecule& molecule, BiopolymerState& bioState,
    BiopolymerImportSession& session, ImportMode mode) {
    
    ImportResult result;
    
    // Parse the input
    BiopolymerImportSession::ParseResult parseResult = session.parse();
    
    if (!parseResult.success) {
        result.success = false;
        result.errorMessage = parseResult.errorMessage;
        
        for (const auto& issue : parseResult.issues) {
            if (issue.severity == ImportIssueSeverity::Error) {
                result.validationIssues.append(issue.toString());
            } else {
                result.warnings.append(issue.toString());
            }
        }
        
        return result;
    }
    
    // Check for blocking issues
    if (!session.canCommit()) {
        result.success = false;
        result.errorMessage = "Cannot commit due to validation errors";
        
        for (const auto& issue : parseResult.blockingIssues()) {
            result.validationIssues.append(issue.toString());
        }
        
        for (const auto& issue : parseResult.warnings()) {
            result.warnings.append(issue.toString());
        }
        
        return result;
    }
    
    // Get preview information
    BiopolymerImportSession::ImportPreview preview = session.generatePreview();
    
    // Save previous state for potential rollback
    EditableMolecule* previousMolecule = nullptr;
    BiopolymerState* previousBiopolymerState = nullptr;
    
    if (mode == ImportMode::Replace) {
        previousMolecule = new EditableMolecule(molecule);
        previousBiopolymerState = new BiopolymerState(bioState);
    }
    
    try {
        // Clear or prepare the document
        if (mode == ImportMode::OpenAsNew || mode == ImportMode::Replace) {
            molecule.clear();
            bioState = BiopolymerState(); // Reset biopolymer state
        }
        
        // Commit the import
        bool commitSuccess = session.commitToState(bioState, molecule,
            [&](int current, int total, const QString& desc) {
                // Progress callback (not used in this implementation)
                qDebug() << "Import progress:" << current << "/" << total << "-" << desc;
            });
        
        if (!commitSuccess) {
            // Restore previous state on failure
            restorePreviousState(molecule, bioState, previousMolecule, previousBiopolymerState);
            
            result.success = false;
            result.errorMessage = "Failed to commit import";
            return result;
        }
        
        // Update result with success information
        result.success = true;
        result.chainCount = preview.chainCount;
        result.residueCount = preview.residueCount;
        result.atomCount = preview.estimatedAtomCount;
        result.bondCount = preview.estimatedBondCount;
        
        // Collect warnings
        for (const auto& warning : preview.warnings) {
            result.warnings.append(warning.message);
        }
        
        // Clean up previous state (no longer needed)
        delete previousMolecule;
        delete previousBiopolymerState;
        
    } catch (...) {
        // Restore on exception
        restorePreviousState(molecule, bioState, previousMolecule, previousBiopolymerState);
        delete previousMolecule;
        delete previousBiopolymerState;
        
        result.success = false;
        result.errorMessage = "Exception during import";
        return result;
    }
    
    return result;
}

void BiopolymerImportManager::restorePreviousState(EditableMolecule& molecule,
                                                  BiopolymerState& bioState,
                                                  EditableMolecule* previousMolecule,
                                                  BiopolymerState* previousBiopolymerState) const {
    if (previousMolecule) {
        molecule = *previousMolecule;
    }
    if (previousBiopolymerState) {
        bioState = *previousBiopolymerState;
    }
}
