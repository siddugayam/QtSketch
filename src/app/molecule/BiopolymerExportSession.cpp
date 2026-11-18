// src/app/molecule/BiopolymerExportSession.cpp
#include "BiopolymerExportSession.h"
#include "BiopolymerState.h"
#include "BiopolymerWorkspace.h"
#include <cmath>

// ========== Export Format Helpers ==========

QString exportFormatToString(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence: return "sequence";
        case ExportFormat::Fasta: return "fasta";
        case ExportFormat::Helm: return "helm";
        case ExportFormat::Idt: return "idt";
        case ExportFormat::AxoLabs: return "axolabs";
        case ExportFormat::Unknown: return "unknown";
    }
    return "unknown";
}

ExportFormat exportFormatFromString(const QString& formatString) {
    QString lower = formatString.toLower();
    if (lower == "sequence" || lower == "seq") return ExportFormat::Sequence;
    if (lower == "fasta" || lower == "fa") return ExportFormat::Fasta;
    if (lower == "helm") return ExportFormat::Helm;
    if (lower == "idt") return ExportFormat::Idt;
    if (lower == "axolabs" || lower == "axo") return ExportFormat::AxoLabs;
    return ExportFormat::Unknown;
}

QList<ExportFormat> allExportFormats() {
    return {
        ExportFormat::Sequence,
        ExportFormat::Fasta,
        ExportFormat::Helm,
        ExportFormat::Idt,
        ExportFormat::AxoLabs
    };
}

// ========== ExportCapability Implementation ==========

bool ExportCapability::canExport() const {
    return isSupported && !requiresExplicitConversion;
}

bool ExportCapability::needsUserConfirmation() const {
    return requiresExplicitConversion || !blockedReasons.isEmpty();
}

QString ExportCapability::diagnosticMessage() const {
    if (isSupported && !requiresExplicitConversion) {
        return QString("Export as %1 is fully supported").arg(exportFormatToString(format));
    }
    
    QStringList messages;
    for (const QString& reason : blockedReasons) {
        messages.append(QString("Blocked: %1").arg(reason));
    }
    for (const QString& change : changedFeatures) {
        messages.append(QString("Will be changed: %1").arg(change));
    }
    for (const QString& omitted : omittedFeatures) {
        messages.append(QString("Will be omitted: %1").arg(omitted));
    }
    
    if (messages.isEmpty()) {
        return "Export may require adjustments";
    }
    
    return messages.join("; ");
}

QString ExportCapability::changeSummary() const {
    QStringList changes;
    for (const QString& change : changedFeatures) {
        changes.append(change);
    }
    for (const QString& omitted : omittedFeatures) {
        changes.append(QString("%1 (omitted)").arg(omitted));
    }
    
    if (changes.isEmpty()) {
        return "No changes required";
    }
    
    return changes.join(", ");
}

bool ExportCapability::allMappingsValid() const {
    return !hasInvalidMappings;
}

bool ExportCapability::formatSupportsMultipleChains(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence: return false;  // Single chain only
        case ExportFormat::Fasta: return true;
        case ExportFormat::Helm: return true;
        case ExportFormat::Idt: return false;  // Typically single chain
        case ExportFormat::AxoLabs: return false;  // Typically single chain
        case ExportFormat::Unknown: return false;
    }
    return false;
}

bool ExportCapability::formatRequiresConversion(ExportFormat format, const BiopolymerState& state) {
    // Check if format requires conversion based on state
    if (state.chainIds().isEmpty()) {
        return false;
    }
    
    // Sequence format requires conversion if multiple chains
    if (format == ExportFormat::Sequence && state.chainIds().size() > 1) {
        return true;
    }
    
    // Check for unsupported modifications
    if (state.hasUnsupportedModifications()) {
        // Some formats may not support modifications
        switch (format) {
            case ExportFormat::Sequence:
            case ExportFormat::Fasta:
                return true;  // Plain sequence can't represent modifications
            case ExportFormat::Helm:
            case ExportFormat::Idt:
            case ExportFormat::AxoLabs:
                return false;  // These formats may support modifications
            case ExportFormat::Unknown:
                return true;
        }
    }
    
    // Check for unsupported connections
    if (state.hasUnsupportedConnections()) {
        // HELM can represent connections, others may not
        switch (format) {
            case ExportFormat::Sequence:
            case ExportFormat::Fasta:
            case ExportFormat::Idt:
            case ExportFormat::AxoLabs:
                return true;
            case ExportFormat::Helm:
                return false;
            case ExportFormat::Unknown:
                return true;
        }
    }
    
    return false;
}

// ========== ExportResult Implementation ==========

bool ExportResult::isSuccess() const {
    return success && !outputText.isEmpty();
}

bool ExportResult::isFailure() const {
    return !success || outputText.isEmpty();
}

QString ExportResult::formattedOutput() const {
    if (!isSuccess()) {
        return errorMessage;
    }
    
    QString result = outputText;
    
    if (!warnings.isEmpty()) {
        result += "\n\n// Warnings:\n";
        for (const QString& warning : warnings) {
            result += QString("// - %1\n").arg(warning);
        }
    }
    
    return result;
}

// ========== ExportOptions Implementation ==========

QMap<QString, QVariant> ExportOptions::toJson() const {
    QMap<QString, QVariant> json;
    json["format"] = exportFormatToString(format);
    json["fastaTitle"] = fastaTitle;
    json["fastaIncludeDescription"] = fastaIncludeDescription;
    json["fastaLineWrap"] = fastaLineWrap;
    json["sequenceSingleChainOnly"] = sequenceSingleChainOnly;
    json["sequenceWarnOnModifications"] = sequenceWarnOnModifications;
    json["helmIncludePolymerSection"] = helmIncludePolymerSection;
    json["helmIncludeAnnotationSection"] = helmIncludeAnnotationSection;
    json["includeChainNames"] = includeChainNames;
    json["includeResidueNumbers"] = includeResidueNumbers;
    json["includeTerminalGroups"] = includeTerminalGroups;
    json["forceSimplification"] = forceSimplification;
    return json;
}

ExportOptions ExportOptions::fromJson(const QMap<QString, QVariant>& json) {
    ExportOptions options;
    options.format = exportFormatFromString(json.value("format", "fasta").toString());
    options.fastaTitle = json.value("fastaTitle", "Biopolymer").toString();
    options.fastaIncludeDescription = json.value("fastaIncludeDescription", true).toBool();
    options.fastaLineWrap = json.value("fastaLineWrap", 60).toInt();
    options.sequenceSingleChainOnly = json.value("sequenceSingleChainOnly", true).toBool();
    options.sequenceWarnOnModifications = json.value("sequenceWarnOnModifications", true).toBool();
    options.helmIncludePolymerSection = json.value("helmIncludePolymerSection", true).toBool();
    options.helmIncludeAnnotationSection = json.value("helmIncludeAnnotationSection", true).toBool();
    options.includeChainNames = json.value("includeChainNames", true).toBool();
    options.includeResidueNumbers = json.value("includeResidueNumbers", false).toBool();
    options.includeTerminalGroups = json.value("includeTerminalGroups", true).toBool();
    options.forceSimplification = json.value("forceSimplification", false).toBool();
    return options;
}

// ========== ExportDiagnostic Implementation ==========

bool ExportDiagnostic::isError() const {
    return severity == "Error";
}

bool ExportDiagnostic::isWarning() const {
    return severity == "Warning";
}

QString ExportDiagnostic::formattedMessage() const {
    if (chainId >= 0 && residueId >= 0) {
        return QString("[%1] %2: %3 - %4").arg(category, severity, message, feature);
    } else if (chainId >= 0) {
        return QString("[%1] %2: %3").arg(category, severity, message);
    }
    return QString("[%1] %2: %3").arg(category, severity, message);
}

// ========== ExportValidation Implementation ==========

bool ExportValidation::hasErrors() const {
    for (const ExportDiagnostic& diag : diagnostics) {
        if (diag.isError()) {
            return true;
        }
    }
    return false;
}

bool ExportValidation::hasWarnings() const {
    for (const ExportDiagnostic& diag : diagnostics) {
        if (diag.isWarning()) {
            return true;
        }
    }
    return false;
}

bool ExportValidation::hasBlockers() const {
    return hasErrors();
}

QList<QString> ExportValidation::errorMessages() const {
    QList<QString> errors;
    for (const ExportDiagnostic& diag : diagnostics) {
        if (diag.isError()) {
            errors.append(diag.message);
        }
    }
    return errors;
}

QList<QString> ExportValidation::warningMessages() const {
    QList<QString> warnings;
    for (const ExportDiagnostic& diag : diagnostics) {
        if (diag.isWarning()) {
            warnings.append(diag.message);
        }
    }
    return warnings;
}

void ExportValidation::merge(const ExportValidation& other) {
    if (!other.isValid) {
        isValid = false;
    }
    diagnostics.append(other.diagnostics);
}

// ========== BiopolymerExportSession Implementation ==========

BiopolymerExportSession::BiopolymerExportSession(const BiopolymerState& state,
                                                 const BiopolymerWorkspace& workspace)
    : m_state(state), m_workspace(workspace) {
    // Initialize default document info
    m_documentName = "Untitled";
}

void BiopolymerExportSession::initialize(int documentId, quint64 documentRevision, 
                                        const QString& documentName) {
    m_documentId = documentId;
    m_documentRevision = documentRevision;
    m_documentName = documentName;
}

// ========== Capability Assessment ==========

ExportCapability BiopolymerExportSession::assessCapability(ExportFormat format) const {
    ExportCapability capability;
    capability.format = format;
    capability.documentId = m_documentId;
    capability.documentRevision = m_documentRevision;
    capability.documentName = m_documentName;
    
    // Basic state assessment
    capability.chainCount = m_state.chainIds().size();
    capability.totalResidueCount = m_state.allResidueIds().size();
    capability.hasValidMappings = m_state.hasValidMappings();
    capability.hasInvalidMappings = m_state.hasInvalidMappings();
    capability.hasUnsupportedModifications = m_state.hasUnsupportedModifications();
    capability.hasUnsupportedConnections = m_state.hasUnsupportedConnections();
    capability.hasMultipleChains = capability.chainCount > 1;
    capability.hasTerminalGroups = m_state.hasTerminalGroups();
    capability.hasBranchedTopology = m_state.hasBranchedTopology();
    
    // Collect polymer types
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (chain) {
            capability.polymerTypes.insert(chain->polymerType);
        }
    }
    
    // Check format-specific support
    switch (format) {
        case ExportFormat::Sequence: {
            // Plain sequence export
            // Block if multiple chains (unless overridden by options)
            if (capability.hasMultipleChains) {
                capability.isSupported = false;
                capability.blockedReasons.append("Multiple chains - plain sequence only supports single chain");
            }
            
            // Block if has invalid mappings
            if (capability.hasInvalidMappings) {
                capability.isSupported = false;
                capability.blockedReasons.append("Invalid residue-to-atom mappings");
            }
            
            // Check if modifications require conversion
            if (capability.hasUnsupportedModifications) {
                capability.requiresExplicitConversion = true;
                capability.changedFeatures.append("Modified residues will be simplified to standard");
            }
            
            capability.isSupported = capability.blockedReasons.isEmpty();
            capability.statusMessage = capability.isSupported ? "Success" : 
                                       (capability.requiresExplicitConversion ? "ConversionRequired" : "Unsupported");
            break;
        }
        
        case ExportFormat::Fasta: {
            // FASTA format supports multiple chains
            if (capability.hasInvalidMappings) {
                capability.isSupported = false;
                capability.blockedReasons.append("Invalid residue-to-atom mappings");
            }
            
            if (capability.hasUnsupportedModifications) {
                capability.requiresExplicitConversion = true;
                capability.changedFeatures.append("Modified residues will be represented as standard");
            }
            
            capability.isSupported = capability.blockedReasons.isEmpty();
            capability.statusMessage = capability.isSupported ? "Success" : 
                                       (capability.requiresExplicitConversion ? "ConversionRequired" : "Unsupported");
            break;
        }
        
        case ExportFormat::Helm: {
            // HELM can represent complex structures
            if (capability.hasInvalidMappings) {
                capability.isSupported = false;
                capability.blockedReasons.append("Invalid residue-to-atom mappings");
            }
            
            // HELM can represent modifications and connections
            if (capability.hasUnsupportedModifications) {
                capability.requiresExplicitConversion = true;
                capability.changedFeatures.append("Some modifications may not be fully represented");
            }
            
            if (capability.hasUnsupportedConnections) {
                capability.requiresExplicitConversion = true;
                capability.changedFeatures.append("Some connections may not be fully represented");
            }
            
            capability.isSupported = capability.blockedReasons.isEmpty();
            capability.statusMessage = capability.isSupported ? "Success" : 
                                       (capability.requiresExplicitConversion ? "ConversionRequired" : "Unsupported");
            break;
        }
        
        case ExportFormat::Idt: {
            // IDT format (DNA/RNA)
            if (capability.hasInvalidMappings) {
                capability.isSupported = false;
                capability.blockedReasons.append("Invalid residue-to-atom mappings");
            }
            
            // Check polymer type - IDT is for nucleotides
            bool hasOnlyNucleotides = true;
            for (PolymerType type : capability.polymerTypes) {
                if (type != PolymerType::DNA && type != PolymerType::RNA) {
                    hasOnlyNucleotides = false;
                    break;
                }
            }
            
            if (!hasOnlyNucleotides) {
                capability.isSupported = false;
                capability.blockedReasons.append("IDT format only supports DNA/RNA");
            }
            
            capability.isSupported = capability.blockedReasons.isEmpty();
            capability.statusMessage = capability.isSupported ? "Success" : "Unsupported";
            break;
        }
        
        case ExportFormat::AxoLabs: {
            // AxoLabs format
            if (capability.hasInvalidMappings) {
                capability.isSupported = false;
                capability.blockedReasons.append("Invalid residue-to-atom mappings");
            }
            
            capability.isSupported = capability.blockedReasons.isEmpty();
            capability.statusMessage = capability.isSupported ? "Success" : "Unsupported";
            break;
        }
        
        case ExportFormat::Unknown: {
            capability.isSupported = false;
            capability.statusMessage = "Unsupported";
            capability.blockedReasons.append("Unknown export format");
            break;
        }
    }
    
    return capability;
}

QMap<ExportFormat, ExportCapability> BiopolymerExportSession::assessAllCapabilities() const {
    QMap<ExportFormat, ExportCapability> capabilities;
    for (ExportFormat format : allExportFormats()) {
        capabilities[format] = assessCapability(format);
    }
    return capabilities;
}

QList<ExportFormat> BiopolymerExportSession::supportedFormats() const {
    QList<ExportFormat> formats;
    for (ExportFormat format : allExportFormats()) {
        if (isFormatSupported(format)) {
            formats.append(format);
        }
    }
    return formats;
}

bool BiopolymerExportSession::isFormatSupported(ExportFormat format) const {
    ExportCapability capability = assessCapability(format);
    return capability.isSupported;
}

// ========== Validation ==========

ExportValidation BiopolymerExportSession::validateForExport(ExportFormat format) const {
    ExportValidation validation;
    validation.format = format;
    validation.isValid = true;
    
    // Check mappings
    if (m_state.hasInvalidMappings()) {
        validation.isValid = false;
        validation.diagnostics.append({
            "Mapping", "Error", 
            "Document contains invalid residue-to-atom mappings",
            -1, -1, "mapping"
        });
    }
    
    // Check chain count for Sequence format
    if (format == ExportFormat::Sequence && m_state.chainIds().size() > 1) {
        validation.isValid = false;
        validation.diagnostics.append({
            "Format", "Error",
            QString("Sequence format only supports single chain, document has %1 chains")
                .arg(m_state.chainIds().size()),
            -1, -1, "chainCount"
        });
    }
    
    // Check for unsupported modifications
    if (m_state.hasUnsupportedModifications()) {
        validation.diagnostics.append({
            "Chemistry", "Warning",
            "Document contains unsupported modifications",
            -1, -1, "modifications"
        });
    }
    
    // Check for unsupported connections
    if (m_state.hasUnsupportedConnections()) {
        validation.diagnostics.append({
            "Topology", "Warning",
            "Document contains unsupported connections",
            -1, -1, "connections"
        });
    }
    
    // Check for terminal groups (may be a warning for some formats)
    if (m_state.hasTerminalGroups()) {
        validation.diagnostics.append({
            "Chemistry", "Info",
            "Document contains terminal groups",
            -1, -1, "terminalGroups"
        });
    }
    
    // Format-specific validation
    if (format == ExportFormat::Idt) {
        // IDT only supports nucleotides
        for (int chainId : m_state.chainIds()) {
            const BiopolymerChain* chain = m_state.chain(chainId);
            if (chain && chain->polymerType != PolymerType::DNA && 
                chain->polymerType != PolymerType::RNA) {
                validation.isValid = false;
                validation.diagnostics.append({
                    "Format", "Error",
                    QString("IDT format only supports DNA/RNA, chain %1 is %2")
                        .arg(chainId).arg(BiopolymerState::polymerTypeToString(chain->polymerType)),
                    chainId, -1, "polymerType"
                });
            }
        }
    }
    
    return validation;
}

bool BiopolymerExportSession::validateAllMappings() const {
    return m_state.hasValidMappings() && !m_state.hasInvalidMappings();
}

bool BiopolymerExportSession::validateNoUnsupportedModifications() const {
    return !m_state.hasUnsupportedModifications();
}

bool BiopolymerExportSession::validateChainCount(ExportFormat format) const {
    if (format == ExportFormat::Sequence) {
        return m_state.chainIds().size() <= 1;
    }
    return true;
}

// ========== Export Generation ==========

ExportResult BiopolymerExportSession::generateExport(ExportFormat format, 
                                                   const ExportOptions& options) const {
    ExportResult result;
    result.format = format;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    result.chainCount = m_state.chainIds().size();
    result.residueCount = m_state.allResidueIds().size();
    result.mappingsValid = validateAllMappings();
    
    // Generate export based on format
    switch (format) {
        case ExportFormat::Sequence:
            return generateSequenceExport(options);
        case ExportFormat::Fasta:
            return generateFastaExport(options);
        case ExportFormat::Helm:
            return generateHelmExport(options);
        case ExportFormat::Idt:
            return generateIdtExport(options);
        case ExportFormat::AxoLabs:
            return generateAxoLabsExport(options);
        case ExportFormat::Unknown:
            result.success = false;
            result.errorMessage = "Unknown export format";
            return result;
    }
    
    return result;
}

ExportResult BiopolymerExportSession::generateSequenceExport(const ExportOptions& options) const {
    ExportResult result;
    result.format = ExportFormat::Sequence;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    
    // Validate
    if (m_state.chainIds().size() > 1) {
        result.success = false;
        result.errorMessage = "Sequence export only supports single chain";
        result.validationIssues.append("Multiple chains not supported");
        return result;
    }
    
    if (!validateAllMappings()) {
        result.success = false;
        result.errorMessage = "Invalid residue-to-atom mappings";
        result.validationIssues.append("Invalid mappings");
        return result;
    }
    
    // Build sequence from first chain
    if (m_state.chainIds().isEmpty()) {
        result.success = false;
        result.errorMessage = "No chains in document";
        return result;
    }
    
    int chainId = m_state.chainIds().first();
    QString sequence = buildSequenceString(chainId, options.includeResidueNumbers);
    
    // Check for modifications
    if (m_state.hasUnsupportedModifications() && options.sequenceWarnOnModifications) {
        result.warnings.append("Document contains modifications that may not be fully represented");
    }
    
    result.success = true;
    result.outputText = sequence;
    result.residueCount = m_state.chainResidueCount(chainId);
    result.chainCount = 1;
    
    return result;
}

ExportResult BiopolymerExportSession::generateFastaExport(const ExportOptions& options) const {
    ExportResult result;
    result.format = ExportFormat::Fasta;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    
    // Validate
    if (!validateAllMappings()) {
        result.success = false;
        result.errorMessage = "Invalid residue-to-atom mappings";
        result.validationIssues.append("Invalid mappings");
        return result;
    }
    
    // Build FASTA for all chains
    QString fastaTitle = options.fastaTitle;
    if (fastaTitle.isEmpty()) {
        fastaTitle = m_documentName;
    }
    
    result.outputText = buildFastaForAllChains(fastaTitle, options.fastaLineWrap);
    result.success = true;
    result.residueCount = m_state.allResidueIds().size();
    result.chainCount = m_state.chainIds().size();
    
    // Check for modifications
    if (m_state.hasUnsupportedModifications()) {
        result.warnings.append("Document contains modifications that may not be fully represented in FASTA");
    }
    
    return result;
}

ExportResult BiopolymerExportSession::generateHelmExport(const ExportOptions& options) const {
    ExportResult result;
    result.format = ExportFormat::Helm;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    
    // Validate
    if (!validateAllMappings()) {
        result.success = false;
        result.errorMessage = "Invalid residue-to-atom mappings";
        result.validationIssues.append("Invalid mappings");
        return result;
    }
    
    // Build HELM notation
    result.outputText = buildHelmNotation();
    result.success = !result.outputText.isEmpty();
    result.residueCount = m_state.allResidueIds().size();
    result.chainCount = m_state.chainIds().size();
    
    if (m_state.hasUnsupportedModifications() || m_state.hasUnsupportedConnections()) {
        result.warnings.append("Some features may not be fully represented in HELM");
    }
    
    return result;
}

ExportResult BiopolymerExportSession::generateIdtExport(const ExportOptions& options) const {
    ExportResult result;
    result.format = ExportFormat::Idt;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    
    // Validate polymer types
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (chain && chain->polymerType != PolymerType::DNA && 
            chain->polymerType != PolymerType::RNA) {
            result.success = false;
            result.errorMessage = QString("IDT format only supports DNA/RNA, chain %1 is %2")
                .arg(chainId).arg(BiopolymerState::polymerTypeToString(chain->polymerType));
            return result;
        }
    }
    
    if (!validateAllMappings()) {
        result.success = false;
        result.errorMessage = "Invalid residue-to-atom mappings";
        result.validationIssues.append("Invalid mappings");
        return result;
    }
    
    // Build IDT sequence
    result.outputText = buildIdtSequence();
    result.success = !result.outputText.isEmpty();
    result.residueCount = m_state.allResidueIds().size();
    result.chainCount = m_state.chainIds().size();
    
    return result;
}

ExportResult BiopolymerExportSession::generateAxoLabsExport(const ExportOptions& options) const {
    ExportResult result;
    result.format = ExportFormat::AxoLabs;
    result.documentId = m_documentId;
    result.documentRevision = m_documentRevision;
    
    if (!validateAllMappings()) {
        result.success = false;
        result.errorMessage = "Invalid residue-to-atom mappings";
        result.validationIssues.append("Invalid mappings");
        return result;
    }
    
    // Build AxoLabs format
    result.outputText = buildAxoLabsFormat();
    result.success = !result.outputText.isEmpty();
    result.residueCount = m_state.allResidueIds().size();
    result.chainCount = m_state.chainIds().size();
    
    return result;
}

// ========== Format-Specific Builders ==========

QString BiopolymerExportSession::buildSequenceString(int chainId, bool includeNumbers) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) return "";
    
    QStringList symbols;
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = m_state.residue(residueId);
        if (residue) {
            QString symbol = residue->identity.symbol;
            if (includeNumbers) {
                symbol += QString("%1").arg(residue->sequencePosition);
            }
            symbols.append(symbol);
        }
    }
    
    return symbols.join("");
}

QMap<int, QString> BiopolymerExportSession::buildAllSequenceStrings(bool includeNumbers) const {
    QMap<int, QString> sequences;
    for (int chainId : m_state.chainIds()) {
        sequences[chainId] = buildSequenceString(chainId, includeNumbers);
    }
    return sequences;
}

QString BiopolymerExportSession::buildFastaRecord(int chainId, const QString& title, int lineWrap) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) return "";
    
    QString sequence = buildSequenceString(chainId, false);
    QString chainName = chain->name;
    
    // Build header
    QString header = ">";
    if (!title.isEmpty()) {
        header += title;
        if (!chainName.isEmpty()) {
            header += QString(" | Chain: %1").arg(chainName);
        }
    } else {
        header += QString("Chain_%1").arg(chainId);
    }
    
    // Wrap sequence
    QStringList lines;
    lines.append(header);
    
    for (int i = 0; i < sequence.size(); i += lineWrap) {
        int end = std::min(i + lineWrap, sequence.size());
        lines.append(sequence.mid(i, end - i));
    }
    
    return lines.join("\n");
}

QString BiopolymerExportSession::buildFastaForAllChains(const QString& baseTitle, int lineWrap) const {
    QStringList records;
    
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (chain) {
            QString title = baseTitle;
            if (!baseTitle.isEmpty() && chain->name.isEmpty()) {
                // Use base title for all chains
                records.append(buildFastaRecord(chainId, title, lineWrap));
            } else {
                // Include chain name in title
                QString chainTitle = baseTitle;
                if (!chain->name.isEmpty()) {
                    if (!chainTitle.isEmpty()) {
                        chainTitle += QString(" | %1").arg(chain->name);
                    } else {
                        chainTitle = chain->name;
                    }
                }
                records.append(buildFastaRecord(chainId, chainTitle, lineWrap));
            }
        }
    }
    
    return records.join("\n");
}

QString BiopolymerExportSession::buildHelmNotation() const {
    // HELM notation format: POLYMER1, POLYMER2, ... | CONNECTIONS | ANNOTATIONS
    
    QStringList polymerSections;
    QList<QString> chainNames;
    
    // Build polymer sections
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (!chain) continue;
        
        QString chainName = chain->name;
        if (chainName.isEmpty()) {
            chainName = QString("CHAIN%1").arg(chainId);
        }
        chainNames.append(chainName);
        
        QStringList residueSymbols;
        for (int residueId : chain->residueIds) {
            const BiopolymerResidue* residue = m_state.residue(residueId);
            if (residue) {
                // Use HELM monomer name
                QString helmSymbol = residue->identity.helmNotation;
                if (helmSymbol.isEmpty()) {
                    helmSymbol = residue->identity.symbol;
                }
                residueSymbols.append(helmSymbol);
            }
        }
        
        if (!residueSymbols.isEmpty()) {
            QString polymerSection = chainName + ":" + residueSymbols.join(",");
            polymerSections.append(polymerSection);
        }
    }
    
    if (polymerSections.isEmpty()) {
        return "";
    }
    
    // Build connections (simplified for now)
    QStringList connectionSections;
    for (const BiopolymerConnection& conn : m_state.connections()) {
        // HELM connection format: chain1.residue1:atom1-chain2.residue2:atom2
        QString connStr = QString("%1.%2:%3-%4.%5:%6")
            .arg(chainNames[conn.fromChainIndex])
            .arg(conn.fromResidueId)
            .arg(BiopolymerState::atomRoleToString(conn.fromRole))
            .arg(chainNames[conn.toChainIndex])
            .arg(conn.toResidueId)
            .arg(BiopolymerState::atomRoleToString(conn.toRole));
        connectionSections.append(connStr);
    }
    
    // Build the full HELM notation
    QString helm;
    helm += polymerSections.join(", ") + "$";
    
    if (!connectionSections.isEmpty()) {
        helm += connectionSections.join(", ") + "$";
    } else {
        helm += "$V2.0$";
    }
    
    return helm;
}

QString BiopolymerExportSession::buildIdtSequence() const {
    // IDT format: just the sequence with optional modifications
    // For simplicity, we'll just build the sequence string
    
    if (m_state.chainIds().isEmpty()) {
        return "";
    }
    
    // For IDT, we typically use the first chain or concatenate all
    QStringList sequences;
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (chain) {
            QString sequence = buildSequenceString(chainId, false);
            if (!sequence.isEmpty()) {
                sequences.append(sequence);
            }
        }
    }
    
    return sequences.join("");
}

QString BiopolymerExportSession::buildAxoLabsFormat() const {
    // AxoLabs format - simplified implementation
    // Similar to IDT, just the sequence
    return buildIdtSequence();
}

// ========== Utility Methods ==========

int BiopolymerExportSession::documentId() const {
    return m_documentId;
}

quint64 BiopolymerExportSession::documentRevision() const {
    return m_documentRevision;
}

QString BiopolymerExportSession::documentName() const {
    return m_documentName;
}

const BiopolymerState& BiopolymerExportSession::state() const {
    return m_state;
}

const BiopolymerWorkspace& BiopolymerExportSession::workspace() const {
    return m_workspace;
}

int BiopolymerExportSession::chainCount() const {
    return m_state.chainIds().size();
}

int BiopolymerExportSession::totalResidueCount() const {
    return m_state.allResidueIds().size();
}

QList<PolymerType> BiopolymerExportSession::polymerTypes() const {
    QList<PolymerType> types;
    for (int chainId : m_state.chainIds()) {
        const BiopolymerChain* chain = m_state.chain(chainId);
        if (chain) {
            if (!types.contains(chain->polymerType)) {
                types.append(chain->polymerType);
            }
        }
    }
    return types;
}

bool BiopolymerExportSession::hasValidMappings() const {
    return m_state.hasValidMappings();
}

bool BiopolymerExportSession::hasInvalidMappings() const {
    return m_state.hasInvalidMappings();
}

bool BiopolymerExportSession::hasMultipleChains() const {
    return m_state.chainIds().size() > 1;
}

bool BiopolymerExportSession::supportsPolymerType(PolymerType type) const {
    return type == PolymerType::PEPTIDE || 
           type == PolymerType::DNA || 
           type == PolymerType::RNA ||
           type == PolymerType::UNKNOWN;
}

bool BiopolymerExportSession::allResiduesStandard(int chainId) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) return false;
    
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = m_state.residue(residueId);
        if (residue && !residue->identity.isStandard) {
            return false;
        }
    }
    return true;
}

QList<QString> BiopolymerExportSession::unsupportedModifications() const {
    QList<QString> modifications;
    
    for (int residueId : m_state.allResidueIds()) {
        const BiopolymerResidue* residue = m_state.residue(residueId);
        if (residue && !residue->identity.modifications.isEmpty()) {
            for (const auto& mod : residue->identity.modifications) {
                if (!mod.isSupported) {
                    modifications.append(QString("%1 at position %2")
                        .arg(mod.name).arg(residue->sequencePosition));
                }
            }
        }
    }
    
    return modifications;
}

QList<QString> BiopolymerExportSession::unsupportedConnections() const {
    QList<QString> connections;
    
    for (const BiopolymerConnection& conn : m_state.connections()) {
        if (!conn.isSupported) {
            connections.append(QString("Connection from %1:%2 to %3:%4")
                .arg(conn.fromChainIndex).arg(conn.fromResidueId)
                .arg(conn.toChainIndex).arg(conn.toResidueId));
        }
    }
    
    return connections;
}

void BiopolymerExportSession::reset() {
    m_lastResult = ExportResult();
}

void BiopolymerExportSession::clearDrafts() {
    // Clear any draft data
    reset();
}

ExportResult BiopolymerExportSession::lastResult() const {
    return m_lastResult;
}

// ========== Static Helpers ==========

bool BiopolymerExportSession::formatSupportsPolymerType(ExportFormat format, PolymerType polymerType) {
    switch (format) {
        case ExportFormat::Sequence:
        case ExportFormat::Fasta:
            return polymerType == PolymerType::PEPTIDE || 
                   polymerType == PolymerType::DNA || 
                   polymerType == PolymerType::RNA;
        case ExportFormat::Helm:
            return true;  // HELM supports all polymer types
        case ExportFormat::Idt:
            return polymerType == PolymerType::DNA || polymerType == PolymerType::RNA;
        case ExportFormat::AxoLabs:
            return polymerType == PolymerType::DNA || polymerType == PolymerType::RNA;
        case ExportFormat::Unknown:
            return false;
    }
    return false;
}

QString BiopolymerExportSession::formatFileExtension(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence: return "seq";
        case ExportFormat::Fasta: return "fasta";
        case ExportFormat::Helm: return "helm";
        case ExportFormat::Idt: return "txt";
        case ExportFormat::AxoLabs: return "txt";
        case ExportFormat::Unknown: return "txt";
    }
    return "txt";
}

QString BiopolymerExportSession::formatMimeType(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence: return "text/plain";
        case ExportFormat::Fasta: return "text/plain";
        case ExportFormat::Helm: return "chemical/x-helm";
        case ExportFormat::Idt: return "text/plain";
        case ExportFormat::AxoLabs: return "text/plain";
        case ExportFormat::Unknown: return "text/plain";
    }
    return "text/plain";
}

QString BiopolymerExportSession::formatDisplayName(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence: return "Sequence";
        case ExportFormat::Fasta: return "FASTA";
        case ExportFormat::Helm: return "HELM";
        case ExportFormat::Idt: return "IDT";
        case ExportFormat::AxoLabs: return "AxoLabs";
        case ExportFormat::Unknown: return "Unknown";
    }
    return "Unknown";
}

QString BiopolymerExportSession::formatDescription(ExportFormat format) {
    switch (format) {
        case ExportFormat::Sequence:
            return "Plain sequence text (single chain only)";
        case ExportFormat::Fasta:
            return "FASTA format (supports multiple chains)";
        case ExportFormat::Helm:
            return "HELM notation (Hierarchical Editing Language for Macromolecules)";
        case ExportFormat::Idt:
            return "IDT format (Integrated DNA Technologies)";
        case ExportFormat::AxoLabs:
            return "AxoLabs format";
        case ExportFormat::Unknown:
            return "Unknown format";
    }
    return "Unknown format";
}

bool BiopolymerExportSession::formatIsTextBased(ExportFormat format) {
    return true;  // All current formats are text-based
}

ExportOptions BiopolymerExportSession::defaultOptions(ExportFormat format) {
    ExportOptions options;
    options.format = format;
    
    switch (format) {
        case ExportFormat::Sequence:
            options.sequenceSingleChainOnly = true;
            options.sequenceWarnOnModifications = true;
            break;
        case ExportFormat::Fasta:
            options.fastaTitle = "Biopolymer";
            options.fastaIncludeDescription = true;
            options.fastaLineWrap = 60;
            break;
        case ExportFormat::Helm:
            options.helmIncludePolymerSection = true;
            options.helmIncludeAnnotationSection = true;
            break;
        case ExportFormat::Idt:
        case ExportFormat::AxoLabs:
            break;
        case ExportFormat::Unknown:
            break;
    }
    
    return options;
}
