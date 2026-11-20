// src/app/molecule/BiopolymerExportCommands.cpp
#include "BiopolymerExportCommands.h"
#include "BiopolymerExportSession.h"
#include "DocumentState.h"

// ========== BiopolymerExportCommand Base Implementation ==========

BiopolymerExportCommand::~BiopolymerExportCommand() {
    // Virtual destructor
}

bool BiopolymerExportCommand::requiresConfirmation() const {
    return false;
}

// ========== ExportSequenceCommand Implementation ==========

ExportSequenceCommand::ExportSequenceCommand(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Sequence;
}

ExportSequenceCommand::~ExportSequenceCommand() {
    // Destructor
}

ExportResult ExportSequenceCommand::execute(DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    session.initialize(-1, document.editRevision(), "Document");
    return session.generateSequenceExport(m_options);
}

QString ExportSequenceCommand::name() const {
    return "Export as Sequence";
}

ExportFormat ExportSequenceCommand::format() const {
    return ExportFormat::Sequence;
}

ExportOptions ExportSequenceCommand::options() const {
    return m_options;
}

bool ExportSequenceCommand::canExecute(const DocumentState& document) const {
    return !document.biopolymerState().chainIds().isEmpty();
}

ExportValidation ExportSequenceCommand::validate(const DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    return session.validateForExport(ExportFormat::Sequence);
}

bool ExportSequenceCommand::requiresConfirmation() const {
    return false;
}

ExportSequenceCommand* ExportSequenceCommand::clone() const {
    return new ExportSequenceCommand(m_options);
}

void ExportSequenceCommand::setOptions(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Sequence;
}

// ========== ExportFastaCommand Implementation ==========

ExportFastaCommand::ExportFastaCommand(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Fasta;
}

ExportFastaCommand::~ExportFastaCommand() {
    // Destructor
}

ExportResult ExportFastaCommand::execute(DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    session.initialize(-1, document.editRevision(), "Document");
    return session.generateFastaExport(m_options);
}

QString ExportFastaCommand::name() const {
    return "Export as FASTA";
}

ExportFormat ExportFastaCommand::format() const {
    return ExportFormat::Fasta;
}

ExportOptions ExportFastaCommand::options() const {
    return m_options;
}

bool ExportFastaCommand::canExecute(const DocumentState& document) const {
    return !document.biopolymerState().chainIds().isEmpty();
}

ExportValidation ExportFastaCommand::validate(const DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    return session.validateForExport(ExportFormat::Fasta);
}

bool ExportFastaCommand::requiresConfirmation() const {
    return false;
}

ExportFastaCommand* ExportFastaCommand::clone() const {
    return new ExportFastaCommand(m_options);
}

void ExportFastaCommand::setOptions(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Fasta;
}

// ========== ExportHelmCommand Implementation ==========

ExportHelmCommand::ExportHelmCommand(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Helm;
}

ExportHelmCommand::~ExportHelmCommand() {
    // Destructor
}

ExportResult ExportHelmCommand::execute(DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    session.initialize(-1, document.editRevision(), "Document");
    return session.generateHelmExport(m_options);
}

QString ExportHelmCommand::name() const {
    return "Export as HELM";
}

ExportFormat ExportHelmCommand::format() const {
    return ExportFormat::Helm;
}

ExportOptions ExportHelmCommand::options() const {
    return m_options;
}

bool ExportHelmCommand::canExecute(const DocumentState& document) const {
    return !document.biopolymerState().chainIds().isEmpty();
}

ExportValidation ExportHelmCommand::validate(const DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    return session.validateForExport(ExportFormat::Helm);
}

bool ExportHelmCommand::requiresConfirmation() const {
    return false;
}

ExportHelmCommand* ExportHelmCommand::clone() const {
    return new ExportHelmCommand(m_options);
}

void ExportHelmCommand::setOptions(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Helm;
}

// ========== ExportIdtCommand Implementation ==========

ExportIdtCommand::ExportIdtCommand(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Idt;
}

ExportIdtCommand::~ExportIdtCommand() {
    // Destructor
}

ExportResult ExportIdtCommand::execute(DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    session.initialize(-1, document.editRevision(), "Document");
    return session.generateIdtExport(m_options);
}

QString ExportIdtCommand::name() const {
    return "Export as IDT";
}

ExportFormat ExportIdtCommand::format() const {
    return ExportFormat::Idt;
}

ExportOptions ExportIdtCommand::options() const {
    return m_options;
}

bool ExportIdtCommand::canExecute(const DocumentState& document) const {
    const BiopolymerState& bioState = document.biopolymerState();
    if (bioState.chainIds().isEmpty()) {
        return false;
    }
    
    // Check if all chains are DNA/RNA
    for (int chainId : bioState.chainIds()) {
        const BiopolymerChain* chain = bioState.chain(chainId);
        if (chain && chain->polymerType != PolymerType::DNA && 
            chain->polymerType != PolymerType::RNA) {
            return false;
        }
    }
    
    return true;
}

ExportValidation ExportIdtCommand::validate(const DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    return session.validateForExport(ExportFormat::Idt);
}

bool ExportIdtCommand::requiresConfirmation() const {
    return false;
}

ExportIdtCommand* ExportIdtCommand::clone() const {
    return new ExportIdtCommand(m_options);
}

void ExportIdtCommand::setOptions(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::Idt;
}

// ========== ExportAxoLabsCommand Implementation ==========

ExportAxoLabsCommand::ExportAxoLabsCommand(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::AxoLabs;
}

ExportAxoLabsCommand::~ExportAxoLabsCommand() {
    // Destructor
}

ExportResult ExportAxoLabsCommand::execute(DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    session.initialize(-1, document.editRevision(), "Document");
    return session.generateAxoLabsExport(m_options);
}

QString ExportAxoLabsCommand::name() const {
    return "Export as AxoLabs";
}

ExportFormat ExportAxoLabsCommand::format() const {
    return ExportFormat::AxoLabs;
}

ExportOptions ExportAxoLabsCommand::options() const {
    return m_options;
}

bool ExportAxoLabsCommand::canExecute(const DocumentState& document) const {
    const BiopolymerState& bioState = document.biopolymerState();
    if (bioState.chainIds().isEmpty()) {
        return false;
    }
    
    // Check if all chains are DNA/RNA
    for (int chainId : bioState.chainIds()) {
        const BiopolymerChain* chain = bioState.chain(chainId);
        if (chain && chain->polymerType != PolymerType::DNA && 
            chain->polymerType != PolymerType::RNA) {
            return false;
        }
    }
    
    return true;
}

ExportValidation ExportAxoLabsCommand::validate(const DocumentState& document) const {
    BiopolymerExportSession session(document.biopolymerState(), *document.bioWorkspace());
    return session.validateForExport(ExportFormat::AxoLabs);
}

bool ExportAxoLabsCommand::requiresConfirmation() const {
    return false;
}

ExportAxoLabsCommand* ExportAxoLabsCommand::clone() const {
    return new ExportAxoLabsCommand(m_options);
}

void ExportAxoLabsCommand::setOptions(const ExportOptions& options) {
    m_options = options;
    m_options.format = ExportFormat::AxoLabs;
}

// ========== ExportCommandFactory Implementation ==========

BiopolymerExportCommand* BiopolymerExportCommandFactory::createCommand(ExportFormat format,
                                                                        const ExportOptions& options) {
    switch (format) {
        case ExportFormat::Sequence:
            return new ExportSequenceCommand(options);
        case ExportFormat::Fasta:
            return new ExportFastaCommand(options);
        case ExportFormat::Helm:
            return new ExportHelmCommand(options);
        case ExportFormat::Idt:
            return new ExportIdtCommand(options);
        case ExportFormat::AxoLabs:
            return new ExportAxoLabsCommand(options);
        case ExportFormat::Unknown:
            return nullptr;
    }
    return nullptr;
}

BiopolymerExportCommand* BiopolymerExportCommandFactory::createCommand(const QString& formatString,
                                                                        const ExportOptions& options) {
    ExportFormat format = exportFormatFromString(formatString);
    return createCommand(format, options);
}

QList<ExportFormat> BiopolymerExportCommandFactory::availableFormats() {
    return allExportFormats();
}

QString BiopolymerExportCommandFactory::formatDisplayName(ExportFormat format) {
    return BiopolymerExportSession::formatDisplayName(format);
}

QString BiopolymerExportCommandFactory::formatDescription(ExportFormat format) {
    return BiopolymerExportSession::formatDescription(format);
}

// ========== BiopolymerExportManager Implementation ==========

BiopolymerExportManager::BiopolymerExportManager(DocumentState& document)
    : m_document(document) {
    // Initialize
}

ExportResult BiopolymerExportManager::exportTo(ExportFormat format, const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    return session.generateExport(format, options);
}

ExportResult BiopolymerExportManager::exportToSequence(const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    ExportOptions opts = options;
    opts.format = ExportFormat::Sequence;
    return session.generateSequenceExport(opts);
}

ExportResult BiopolymerExportManager::exportToFasta(const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    ExportOptions opts = options;
    opts.format = ExportFormat::Fasta;
    return session.generateFastaExport(opts);
}

ExportResult BiopolymerExportManager::exportToHelm(const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    ExportOptions opts = options;
    opts.format = ExportFormat::Helm;
    return session.generateHelmExport(opts);
}

ExportResult BiopolymerExportManager::exportToIdt(const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    ExportOptions opts = options;
    opts.format = ExportFormat::Idt;
    return session.generateIdtExport(opts);
}

ExportResult BiopolymerExportManager::exportToAxoLabs(const ExportOptions& options) {
    BiopolymerExportSession session = createSession();
    ExportOptions opts = options;
    opts.format = ExportFormat::AxoLabs;
    return session.generateAxoLabsExport(opts);
}

ExportValidation BiopolymerExportManager::validateExport(ExportFormat format) const {
    BiopolymerExportSession session = createSession();
    return session.validateForExport(format);
}

bool BiopolymerExportManager::isExportSupported(ExportFormat format) const {
    BiopolymerExportSession session = createSession();
    return session.isFormatSupported(format);
}

QList<ExportFormat> BiopolymerExportManager::supportedFormats() const {
    BiopolymerExportSession session = createSession();
    return session.supportedFormats();
}

QMap<ExportFormat, ExportCapability> BiopolymerExportManager::assessCapabilities() const {
    BiopolymerExportSession session = createSession();
    return session.assessAllCapabilities();
}

ExportOptions BiopolymerExportManager::defaultOptions(ExportFormat format) {
    return BiopolymerExportSession::defaultOptions(format);
}

BiopolymerExportSession BiopolymerExportManager::createSession() const {
    return BiopolymerExportSession(m_document.biopolymerState(), *m_document.bioWorkspace());
}
