// src/app/molecule/BiopolymerImportSession.cpp
#include "BiopolymerImportSession.h"
#include "IndigoMonomerAdapter.h"
#include <QRegularExpression>
#include <QTextStream>
#include <algorithm>

// ========== SourceLocation Implementation ==========

QString SourceLocation::toString() const {
    QStringList parts;
    if (line >= 1) {
        parts << QString("line %1").arg(line);
        if (column >= 1) {
            parts << QString("col %1").arg(column);
        }
    }
    if (!recordId.isEmpty()) {
        parts << QString("record '%1'").arg(recordId);
    }
    if (recordIndex >= 0) {
        parts << QString("record #%1").arg(recordIndex);
    }
    return parts.join(", ");
}

// ========== ImportIssue Implementation ==========

QString ImportIssue::toString() const {
    QStringList parts;
    
    // Severity prefix
    switch (severity) {
        case ImportIssueSeverity::Error:
            parts << "ERROR";
            break;
        case ImportIssueSeverity::Warning:
            parts << "WARNING";
            break;
        case ImportIssueSeverity::Info:
            parts << "INFO";
            break;
    }
    
    // Code
    if (!code.isEmpty()) {
        parts << QString("[%1]").arg(code);
    }
    
    // Location
    if (location.isValid()) {
        parts << location.toString();
    }
    
    // Message
    parts << message;
    
    // Suggestion
    if (!suggestion.isEmpty()) {
        parts << QString("Suggestion: %1").arg(suggestion);
    }
    
    return parts.join(" ");
}

// ========== ParseResult Implementation ==========

bool ParseResult::canAttemptCommit() const {
    if (!success) {
        return false;
    }
    return !hasUnsupportedChemistry;
}

QList<ImportIssue> ParseResult::blockingIssues() const {
    QList<ImportIssue> blocking;
    for (const auto& issue : issues) {
        if (issue.blocksCommit()) {
            blocking.append(issue);
        }
    }
    return blocking;
}

QList<ImportIssue> ParseResult::warnings() const {
    QList<ImportIssue> result;
    for (const auto& issue : issues) {
        if (issue.severity == ImportIssueSeverity::Warning || 
            issue.severity == ImportIssueSeverity::Info) {
            result.append(issue);
        }
    }
    return result;
}

// ========== BiopolymerImportSession Implementation ==========

BiopolymerImportSession::BiopolymerImportSession() {
    // Initialize with default options
    m_options = ImportOptions();
}

BiopolymerImportSession::BiopolymerImportSession(const QString& sourceText, 
                                                   const ImportOptions& options)
    : m_sourceText(sourceText), m_options(options) {
    // Auto-detect format
    m_format = detectFormatInternal().format;
    
    // Parse the source
    parse();
}

BiopolymerImportSession::~BiopolymerImportSession() {
}

// ========== Source Management ==========

void BiopolymerImportSession::setSourceText(const QString& text) {
    m_sourceText = text;
    m_isParsed = false;
    m_parseSuccess = false;
    m_parsedChains.clear();
    m_issues.clear();
    m_unsupportedFeatures.clear();
    
    // Auto-detect format
    if (!m_formatExplicit) {
        m_format = detectFormatInternal().format;
    }
}

QString BiopolymerImportSession::sourceText() const {
    return m_sourceText;
}

void BiopolymerImportSession::clearSource() {
    m_sourceText.clear();
    reset();
}

// ========== Format Detection ==========

FormatDetection BiopolymerImportSession::detectFormat() const {
    return detectFormatInternal();
}

ImportFormat BiopolymerImportSession::format() const {
    return m_format;
}

void BiopolymerImportSession::setFormat(ImportFormat format) {
    m_format = format;
    m_formatExplicit = true;
}

FormatDetection BiopolymerImportSession::detectFormatInternal() const {
    if (m_sourceText.isEmpty()) {
        return {ImportFormat::Unknown, 0.0, "Empty input"};
    }
    
    QString text = m_sourceText.trimmed();
    
    // Check for FASTA format (starts with '>')
    if (text.startsWith('>')) {
        return {ImportFormat::FASTA, 0.9, "Starts with FASTA header"};
    }
    
    // Check for FASTA multi-record (contains '>' in middle)
    int fastaCount = text.count('>');
    if (fastaCount > 1) {
        return {ImportFormat::FASTA, 0.85, QString("Contains %1 FASTA headers").arg(fastaCount)};
    }
    
    // Check for HELM (starts with specific patterns)
    QString upper = text.toUpper();
    if (upper.startsWith("PEPTIDE") || upper.startsWith("RNA") || 
        upper.startsWith("DNA") || upper.contains("$$$$")) {
        return {ImportFormat::HELM, 0.8, "HELM notation pattern"};
    }
    
    // Check for IDT format (specific patterns)
    if (text.contains("IDT") || text.contains("Scale:") || 
        text.contains("Purification:")) {
        return {ImportFormat::IDT, 0.75, "IDT format pattern"};
    }
    
    // Check for simple format (comma-separated chains)
    if (text.contains(',')) {
        QList<QString> parts = text.split(',');
        bool allValid = true;
        for (const QString& part : parts) {
            QString trimmed = part.trimmed();
            if (trimmed.isEmpty()) continue;
            if (!isValidAminoAcid(trimmed) && !isValidNucleotide(trimmed) &&
                !isAmbiguousSymbol(trimmed) && trimmed != "-") {
                allValid = false;
                break;
            }
        }
        if (allValid && parts.size() > 1) {
            return {ImportFormat::Simple, 0.7, "Comma-separated chains"};
        }
    }
    
    // Default to Sequence format
    // Check if it looks like a sequence (only valid monomer symbols)
    bool looksLikeSequence = true;
    for (QChar c : text) {
        if (c.isSpace() || c == '-' || c == '_') continue;
        QString symbol(1, c.toUpper());
        if (!isValidAminoAcid(symbol) && !isValidNucleotide(symbol) && 
            !isAmbiguousSymbol(symbol)) {
            looksLikeSequence = false;
            break;
        }
    }
    
    if (looksLikeSequence) {
        return {ImportFormat::Sequence, 0.6, "Looks like a biopolymer sequence"};
    }
    
    return {ImportFormat::Unknown, 0.0, "Unknown format"};
}

// ========== Parsing ==========

ParseResult BiopolymerImportSession::parse() {
    if (m_sourceText.isEmpty()) {
        return ParseResult{false, "Empty source text", {}, 
                          {{ImportIssueSeverity::Error, "EMPTY_SOURCE", 
                            "Source text is empty", {}, "Provide input text"}}, {},
                          false, false};
    }
    
    m_isParsed = true;
    
    // Use explicit format if set
    if (m_formatExplicit && m_format != ImportFormat::Unknown) {
        return parseAs(m_format);
    }
    
    // Auto-detect and parse
    FormatDetection detection = detectFormatInternal();
    m_format = detection.format;
    
    switch (m_format) {
        case ImportFormat::FASTA:
            return parseFastaInternal(m_sourceText);
        case ImportFormat::Sequence:
            return parseSequenceInternal(m_sourceText);
        case ImportFormat::Simple:
            return parseSimpleInternal(m_sourceText);
        default:
            // Try common formats in order
            ParseResult result;
            
            // Try FASTA
            result = parseFastaInternal(m_sourceText);
            if (result.success) {
                m_format = ImportFormat::FASTA;
                return result;
            }
            
            // Try Sequence
            result = parseSequenceInternal(m_sourceText);
            if (result.success) {
                m_format = ImportFormat::Sequence;
                return result;
            }
            
            // Try Simple
            result = parseSimpleInternal(m_sourceText);
            if (result.success) {
                m_format = ImportFormat::Simple;
                return result;
            }
            
            // All failed
            return ParseResult{false, "Could not parse input in any known format", 
                              {}, 
                              {{ImportIssueSeverity::Error, "PARSE_FAILED", 
                                "Failed to parse input", {}, "Check input format"}}, 
                              {}, false, false};
    }
}

ParseResult BiopolymerImportSession::parseAs(ImportFormat format) {
    m_format = format;
    m_formatExplicit = true;
    
    switch (format) {
        case ImportFormat::FASTA:
            return parseFastaInternal(m_sourceText);
        case ImportFormat::Sequence:
            return parseSequenceInternal(m_sourceText);
        case ImportFormat::Simple:
            return parseSimpleInternal(m_sourceText);
        default:
            return ParseResult{false, QString("Format '%1' not yet implemented").arg(ImportFormatRegistry::formatName(format)),
                              {}, 
                              {{ImportIssueSeverity::Error, "UNSUPPORTED_FORMAT", 
                                QString("Format not supported: %1").arg(ImportFormatRegistry::formatName(format)), 
                                {}, QString("Use FASTA, Sequence, or Simple format")}}, 
                              {}, false, false};
    }
}

ParseResult BiopolymerImportSession::reparse() {
    m_isParsed = false;
    m_parseSuccess = false;
    m_parsedChains.clear();
    m_issues.clear();
    m_unsupportedFeatures.clear();
    return parse();
}

bool BiopolymerImportSession::isParsed() const {
    return m_isParsed;
}

// ========== FASTA Parsing ==========

ParseResult BiopolymerImportSession::parseFastaInternal(const QString& text) {
    ParseResult result;
    result.success = true;
    
    QStringList lines = text.split('\n');
    
    int currentLine = 0;
    int recordIndex = 0;
    ParsedChain* currentChain = nullptr;
    
    for (const QString& line : lines) {
        currentLine++;
        
        QString trimmed = line.trimmed();
        if (trimmed.isEmpty()) {
            continue;
        }
        
        // FASTA header line (starts with '>')
        if (trimmed.startsWith('>')) {
            // Finalize previous chain
            if (currentChain) {
                // Validate the chain
                QList<ImportIssue> chainIssues = validateChain(*currentChain, recordIndex - 1);
                result.issues.append(chainIssues);
                
                for (const auto& issue : chainIssues) {
                    if (issue.blocksCommit()) {
                        result.success = false;
                    }
                }
                
                m_parsedChains.append(*currentChain);
                delete currentChain;
            }
            
            // Parse header
            QString header = trimmed.mid(1).trimmed();
            
            // Split header: first word is ID, rest is description
            QStringList headerParts = header.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            QString recordId = headerParts.isEmpty() ? QString("record_%1").arg(recordIndex) : headerParts[0];
            QString description = headerParts.size() > 1 ? header.split(' ').mid(1).join(' ') : "";
            
            currentChain = new ParsedChain();
            currentChain->name = recordId;
            currentChain->originalName = recordId;
            currentChain->recordIndex = recordIndex;
            currentChain->recordId = recordId;
            currentChain->description = description;
            currentChain->headerLocation = SourceLocation{currentLine, 1, recordIndex, recordId};
            currentChain->polymerType = m_options.defaultPolymerType;
            currentChain->polymerTypeStr = BiopolymerState::polymerTypeToString(m_options.defaultPolymerType);
            
            recordIndex++;
        } else {
            // Sequence line
            if (!currentChain) {
                // Sequence without header
                result.issues.append({
                    ImportIssueSeverity::Error, "MISSING_HEADER",
                    "Sequence line without header",
                    SourceLocation{currentLine, 1, 0, ""},
                    "Add a FASTA header starting with '>'"
                });
                result.success = false;
                continue;
            }
            
            // Track body location
            if (currentChain->bodyLocation.line < 0) {
                currentChain->bodyLocation = SourceLocation{currentLine, 1, currentChain->recordIndex, currentChain->recordId};
            }
            
            // Parse sequence line
            QString sequence = normalizeSequence(line);
            int column = 1;
            
            for (QChar c : sequence) {
                QString symbol(1, c.toUpper());
                
                ParsedMonomer monomer;
                monomer.symbol = symbol;
                monomer.originalSymbol = symbol;
                monomer.position = currentChain->monomers.size() + 1;
                monomer.location = SourceLocation{currentLine, column, currentChain->recordIndex, currentChain->recordId};
                
                // Check if symbol is valid
                if (isValidAminoAcid(symbol)) {
                    monomer.isSupported = true;
                } else if (isValidNucleotide(symbol)) {
                    monomer.isSupported = true;
                    // Update polymer type based on nucleotides
                    if (currentChain->polymerType == PolymerType::Peptide) {
                        currentChain->polymerType = PolymerType::DNA;
                        currentChain->polymerTypeStr = "DNA";
                    }
                } else if (isAmbiguousSymbol(symbol)) {
                    monomer.isAmbiguous = true;
                    monomer.isSupported = m_options.autoResolveAmbiguous;
                    result.hasAmbiguousSymbols = true;
                    
                    result.issues.append({
                        ImportIssueSeverity::Warning, "AMBIGUOUS_SYMBOL",
                        QString("Ambiguous symbol '%1'").arg(symbol),
                        monomer.location,
                        QString("Resolve to a specific monomer or enable auto-resolve")
                    });
                } else if (symbol == "-" || symbol == "_") {
                    // Gap or separator, skip
                    column++;
                    continue;
                } else {
                    // Unsupported symbol
                    monomer.isSupported = false;
                    result.hasUnsupportedChemistry = true;
                    
                    result.issues.append({
                        ImportIssueSeverity::Error, "UNSUPPORTED_SYMBOL",
                        QString("Unsupported symbol '%1'").arg(symbol),
                        monomer.location,
                        QString("Use standard amino acid or nucleotide symbols")
                    });
                    result.success = false;
                    
                    m_unsupportedFeatures.insert(QString("SYMBOL_%1").arg(symbol));
                }
                
                currentChain->monomers.append(monomer);
                column++;
            }
        }
    }
    
    // Finalize last chain
    if (currentChain) {
        QList<ImportIssue> chainIssues = validateChain(*currentChain, recordIndex);
        result.issues.append(chainIssues);
        
        for (const auto& issue : chainIssues) {
            if (issue.blocksCommit()) {
                result.success = false;
            }
        }
        
        m_parsedChains.append(*currentChain);
        delete currentChain;
    }
    
    // Check for unsupported features
    if (m_options.detectUnsupportedFeatures) {
        detectUnsupportedFeaturesInChains(result);
    }
    
    m_parseSuccess = result.success;
    m_issues = result.issues;
    m_unsupportedFeatures = result.unsupportedFeatures;
    
    return result;
}

void BiopolymerImportSession::detectUnsupportedFeaturesInChains(ParseResult& result) const {
    // Check for features that are not yet supported
    for (const ParsedChain& chain : m_parsedChains) {
        // Check for mixed polymer types
        bool hasPeptide = false;
        bool hasNucleotide = false;
        
        for (const ParsedMonomer& monomer : chain.monomers) {
            if (isValidAminoAcid(monomer.symbol)) {
                hasPeptide = true;
            } else if (isValidNucleotide(monomer.symbol)) {
                hasNucleotide = true;
            }
        }
        
        if (hasPeptide && hasNucleotide) {
            result.unsupportedFeatures.insert("MIXED_POLYMER_TYPES");
            result.issues.append({
                ImportIssueSeverity::Warning, "MIXED_POLYMER_TYPES",
                QString("Chain '%1' contains both peptide and nucleotide symbols").arg(chain.name),
                chain.headerLocation,
                "Separate into different chains"
            });
        }
    }
}

// ========== Sequence Parsing ==========

ParseResult BiopolymerImportSession::parseSequenceInternal(const QString& text) {
    ParseResult result;
    result.success = true;
    
    QString sequence = normalizeSequence(text);
    
    if (sequence.isEmpty()) {
        result.success = false;
        result.errorMessage = "Empty sequence";
        result.issues.append({
            ImportIssueSeverity::Error, "EMPTY_SEQUENCE",
            "Sequence is empty",
            SourceLocation{1, 1, 0, ""},
            "Provide a valid biopolymer sequence"
        });
        return result;
    }
    
    ParsedChain chain;
    chain.name = "A";
    chain.originalName = "A";
    chain.polymerType = m_options.defaultPolymerType;
    chain.polymerTypeStr = BiopolymerState::polymerTypeToString(m_options.defaultPolymerType);
    chain.headerLocation = SourceLocation{1, 1, 0, ""};
    chain.bodyLocation = SourceLocation{1, 1, 0, ""};
    chain.recordIndex = 0;
    chain.recordId = "sequence";
    
    int column = 1;
    for (QChar c : sequence) {
        QString symbol(1, c.toUpper());
        
        ParsedMonomer monomer;
        monomer.symbol = symbol;
        monomer.originalSymbol = symbol;
        monomer.position = chain.monomers.size() + 1;
        monomer.location = SourceLocation{1, column, 0, ""};
        
        if (isValidAminoAcid(symbol)) {
            monomer.isSupported = true;
        } else if (isValidNucleotide(symbol)) {
            monomer.isSupported = true;
            if (chain.polymerType == PolymerType::Peptide) {
                chain.polymerType = PolymerType::DNA;
                chain.polymerTypeStr = "DNA";
            }
        } else if (isAmbiguousSymbol(symbol)) {
            monomer.isAmbiguous = true;
            monomer.isSupported = m_options.autoResolveAmbiguous;
            result.hasAmbiguousSymbols = true;
            
            result.issues.append({
                ImportIssueSeverity::Warning, "AMBIGUOUS_SYMBOL",
                QString("Ambiguous symbol '%1'").arg(symbol),
                monomer.location,
                "Resolve to a specific monomer or enable auto-resolve"
            });
        } else if (symbol == "-" || symbol == "_") {
            column++;
            continue;
        } else {
            monomer.isSupported = false;
            result.hasUnsupportedChemistry = true;
            
            result.issues.append({
                ImportIssueSeverity::Error, "UNSUPPORTED_SYMBOL",
                QString("Unsupported symbol '%1'").arg(symbol),
                monomer.location,
                "Use standard amino acid or nucleotide symbols"
            });
            result.success = false;
            
            m_unsupportedFeatures.insert(QString("SYMBOL_%1").arg(symbol));
        }
        
        chain.monomers.append(monomer);
        column++;
    }
    
    // Validate the chain
    QList<ImportIssue> chainIssues = validateChain(chain, 0);
    result.issues.append(chainIssues);
    
    for (const auto& issue : chainIssues) {
        if (issue.blocksCommit()) {
            result.success = false;
        }
    }
    
    m_parsedChains.append(chain);
    m_parseSuccess = result.success;
    m_issues = result.issues;
    m_unsupportedFeatures = result.unsupportedFeatures;
    
    return result;
}

// ========== Simple Format Parsing ==========

ParseResult BiopolymerImportSession::parseSimpleInternal(const QString& text) {
    ParseResult result;
    result.success = true;
    
    // Split by comma to get chains
    QStringList chainTexts = text.split(',', Qt::SkipEmptyParts);
    
    if (chainTexts.isEmpty()) {
        result.success = false;
        result.errorMessage = "No chains found";
        result.issues.append({
            ImportIssueSeverity::Error, "NO_CHAINS",
            "No chains found in input",
            SourceLocation{1, 1, 0, ""},
            "Provide at least one chain"
        });
        return result;
    }
    
    for (int i = 0; i < chainTexts.size(); i++) {
        QString chainText = chainTexts[i].trimmed();
        
        ParsedChain chain;
        chain.name = QString("%1").arg(char('A' + i));
        chain.originalName = chain.name;
        chain.polymerType = m_options.defaultPolymerType;
        chain.polymerTypeStr = BiopolymerState::polymerTypeToString(m_options.defaultPolymerType);
        chain.headerLocation = SourceLocation{1, 1, i, chain.name};
        chain.bodyLocation = SourceLocation{1, 1, i, chain.name};
        chain.recordIndex = i;
        chain.recordId = chain.name;
        
        QString sequence = normalizeSequence(chainText);
        int column = 1;
        
        for (QChar c : sequence) {
            QString symbol(1, c.toUpper());
            
            ParsedMonomer monomer;
            monomer.symbol = symbol;
            monomer.originalSymbol = symbol;
            monomer.position = chain.monomers.size() + 1;
            monomer.location = SourceLocation{1, column, i, chain.name};
            
            if (isValidAminoAcid(symbol)) {
                monomer.isSupported = true;
            } else if (isValidNucleotide(symbol)) {
                monomer.isSupported = true;
                if (chain.polymerType == PolymerType::Peptide) {
                    chain.polymerType = PolymerType::DNA;
                    chain.polymerTypeStr = "DNA";
                }
            } else if (isAmbiguousSymbol(symbol)) {
                monomer.isAmbiguous = true;
                monomer.isSupported = m_options.autoResolveAmbiguous;
                result.hasAmbiguousSymbols = true;
                
                result.issues.append({
                    ImportIssueSeverity::Warning, "AMBIGUOUS_SYMBOL",
                    QString("Ambiguous symbol '%1' in chain '%2'").arg(symbol, chain.name),
                    monomer.location,
                    "Resolve to a specific monomer or enable auto-resolve"
                });
            } else if (symbol == "-" || symbol == "_") {
                column++;
                continue;
            } else {
                monomer.isSupported = false;
                result.hasUnsupportedChemistry = true;
                
                result.issues.append({
                    ImportIssueSeverity::Error, "UNSUPPORTED_SYMBOL",
                    QString("Unsupported symbol '%1' in chain '%2'").arg(symbol, chain.name),
                    monomer.location,
                    "Use standard amino acid or nucleotide symbols"
                });
                result.success = false;
                
                m_unsupportedFeatures.insert(QString("SYMBOL_%1").arg(symbol));
            }
            
            chain.monomers.append(monomer);
            column++;
        }
        
        // Validate the chain
        QList<ImportIssue> chainIssues = validateChain(chain, i);
        result.issues.append(chainIssues);
        
        for (const auto& issue : chainIssues) {
            if (issue.blocksCommit()) {
                result.success = false;
            }
        }
        
        m_parsedChains.append(chain);
    }
    
    // Detect unsupported features
    if (m_options.detectUnsupportedFeatures) {
        ParseResult tempResult = result;
        detectUnsupportedFeaturesInChains(tempResult);
        result.unsupportedFeatures = tempResult.unsupportedFeatures;
        result.issues.append(tempResult.issues);
    }
    
    m_parseSuccess = result.success;
    m_issues = result.issues;
    m_unsupportedFeatures = result.unsupportedFeatures;
    
    return result;
}

// ========== Static Parse Methods ==========

ParseResult BiopolymerImportSession::parseFasta(const QString& text, const ImportOptions& options) {
    BiopolymerImportSession session;
    session.m_options = options;
    session.m_sourceText = text;
    session.m_format = ImportFormat::FASTA;
    session.m_formatExplicit = true;
    return session.parseFastaInternal(text);
}

ParseResult BiopolymerImportSession::parseSequence(const QString& text, const ImportOptions& options) {
    BiopolymerImportSession session;
    session.m_options = options;
    session.m_sourceText = text;
    session.m_format = ImportFormat::Sequence;
    session.m_formatExplicit = true;
    return session.parseSequenceInternal(text);
}

ParseResult BiopolymerImportSession::parseSimple(const QString& text, const ImportOptions& options) {
    BiopolymerImportSession session;
    session.m_options = options;
    session.m_sourceText = text;
    session.m_format = ImportFormat::Simple;
    session.m_formatExplicit = true;
    return session.parseSimpleInternal(text);
}

// ========== Validation ==========

QList<ImportIssue> BiopolymerImportSession::validateChain(const ParsedChain& chain, 
                                                         int chainIndex) const {
    QList<ImportIssue> issues;
    
    if (chain.monomers.isEmpty()) {
        issues.append({
            ImportIssueSeverity::Error, "EMPTY_CHAIN",
            QString("Chain '%1' has no monomers").arg(chain.name),
            chain.headerLocation,
            "Add at least one monomer to the chain"
        });
        return issues;
    }
    
    // Check for duplicate positions
    QSet<int> positions;
    for (const ParsedMonomer& monomer : chain.monomers) {
        if (positions.contains(monomer.position)) {
            issues.append({
                ImportIssueSeverity::Error, "DUPLICATE_POSITION",
                QString("Duplicate position %1 in chain '%2'").arg(monomer.position).arg(chain.name),
                monomer.location,
                "Each monomer must have a unique position"
            });
        }
        positions.insert(monomer.position);
    }
    
    // Check for gaps in positions
    int expectedPosition = 1;
    for (const ParsedMonomer& monomer : chain.monomers) {
        if (monomer.position != expectedPosition) {
            issues.append({
                ImportIssueSeverity::Warning, "POSITION_GAP",
                QString("Position gap in chain '%1': expected %2, found %3")
                    .arg(chain.name).arg(expectedPosition).arg(monomer.position),
                monomer.location,
                "Positions should be consecutive"
            });
        }
        expectedPosition = monomer.position + 1;
    }
    
    return issues;
}

QList<ImportIssue> BiopolymerImportSession::validateMonomer(const ParsedMonomer& monomer,
                                                          int chainIndex, 
                                                          int monomerIndex) const {
    QList<ImportIssue> issues;
    
    if (monomer.symbol.isEmpty()) {
        issues.append({
            ImportIssueSeverity::Error, "EMPTY_SYMBOL",
            "Monomer has empty symbol",
            monomer.location,
            "Provide a valid monomer symbol"
        });
    } else if (!monomer.isSupported && !monomer.isAmbiguous) {
        issues.append({
            ImportIssueSeverity::Error, "UNSUPPORTED_MONOMER",
            QString("Unsupported monomer symbol: '%1'").arg(monomer.symbol),
            monomer.location,
            "Use a standard amino acid or nucleotide symbol"
        });
    }
    
    return issues;
}

// ========== Validation (Public) ==========

QList<ImportIssue> BiopolymerImportSession::validate() const {
    return validate(m_options);
}

QList<ImportIssue> BiopolymerImportSession::validate(const ImportOptions& options) const {
    QList<ImportIssue> issues = m_issues;
    
    // Re-validate with new options if different
    if (options.autoResolveAmbiguous != m_options.autoResolveAmbiguous ||
        options.defaultPolymerType != m_options.defaultPolymerType) {
        
        for (int i = 0; i < m_parsedChains.size(); i++) {
            const ParsedChain& chain = m_parsedChains[i];
            QList<ImportIssue> chainIssues = validateChain(chain, i);
            issues.append(chainIssues);
            
            for (int j = 0; j < chain.monomers.size(); j++) {
                const ParsedMonomer& monomer = chain.monomers[j];
                
                // Check for ambiguous symbols with new settings
                if (isAmbiguousSymbol(monomer.symbol) && !options.autoResolveAmbiguous) {
                    // Find if we already have this issue
                    bool alreadyReported = false;
                    for (const ImportIssue& issue : issues) {
                        if (issue.code == "AMBIGUOUS_SYMBOL" &&
                            issue.location.line == monomer.location.line &&
                            issue.location.column == monomer.location.column) {
                            alreadyReported = true;
                            break;
                        }
                    }
                    
                    if (!alreadyReported) {
                        issues.append({
                            ImportIssueSeverity::Error, "AMBIGUOUS_SYMBOL",
                            QString("Ambiguous symbol '%1' requires resolution").arg(monomer.symbol),
                            monomer.location,
                            "Resolve to a specific monomer"
                        });
                    }
                }
            }
        }
    }
    
    return issues;
}

bool BiopolymerImportSession::canCommit() const {
    if (!m_isParsed || !m_parseSuccess) {
        return false;
    }
    
    // Check for blocking issues
    for (const ImportIssue& issue : m_issues) {
        if (issue.blocksCommit()) {
            return false;
        }
    }
    
    // Check for unsupported chemistry if not allowed
    if (!m_options.allowPartialImport && m_unsupportedFeatures.size() > 0) {
        return false;
    }
    
    return true;
}

// ========== Draft Access ==========

QList<ParsedChain> BiopolymerImportSession::parsedChains() const {
    return m_parsedChains;
}

ParsedChain* BiopolymerImportSession::parsedChain(int index) {
    if (index >= 0 && index < m_parsedChains.size()) {
        return &m_parsedChains[index];
    }
    return nullptr;
}

const ParsedChain* BiopolymerImportSession::parsedChain(int index) const {
    if (index >= 0 && index < m_parsedChains.size()) {
        return &m_parsedChains[index];
    }
    return nullptr;
}

ParsedChain* BiopolymerImportSession::parsedChainByName(const QString& name) {
    for (auto& chain : m_parsedChains) {
        if (chain.name == name || chain.originalName == name) {
            return &chain;
        }
    }
    return nullptr;
}

const ParsedChain* BiopolymerImportSession::parsedChainByName(const QString& name) const {
    for (const auto& chain : m_parsedChains) {
        if (chain.name == name || chain.originalName == name) {
            return &chain;
        }
    }
    return nullptr;
}

int BiopolymerImportSession::monomerCount() const {
    int count = 0;
    for (const auto& chain : m_parsedChains) {
        count += chain.monomers.size();
    }
    return count;
}

QList<ImportIssue> BiopolymerImportSession::allIssues() const {
    return m_issues;
}

QList<ImportIssue> BiopolymerImportSession::issuesForChain(int chainIndex) const {
    QList<ImportIssue> result;
    for (const ImportIssue& issue : m_issues) {
        if (issue.location.recordIndex == chainIndex) {
            result.append(issue);
        }
    }
    return result;
}

QList<ImportIssue> BiopolymerImportSession::issuesForMonomer(int chainIndex, int monomerIndex) const {
    QList<ImportIssue> result;
    const ParsedChain* chain = parsedChain(chainIndex);
    if (!chain || monomerIndex < 0 || monomerIndex >= chain->monomers.size()) {
        return result;
    }
    
    const ParsedMonomer& monomer = chain->monomers[monomerIndex];
    for (const ImportIssue& issue : m_issues) {
        if (issue.location.recordIndex == chainIndex &&
            issue.location.line == monomer.location.line &&
            issue.location.column == monomer.location.column) {
            result.append(issue);
        }
    }
    return result;
}

// ========== Ambiguity Resolution ==========

QList<ParsedMonomer*> BiopolymerImportSession::ambiguousSymbols() {
    QList<ParsedMonomer*> result;
    for (auto& chain : m_parsedChains) {
        for (auto& monomer : chain.monomers) {
            if (monomer.isAmbiguous && !monomer.isSupported) {
                result.append(&monomer);
            }
        }
    }
    return result;
}

bool BiopolymerImportSession::resolveAmbiguousSymbol(int chainIndex, int monomerIndex, 
                                                    const QString& resolvedSymbol) {
    if (chainIndex < 0 || chainIndex >= m_parsedChains.size()) {
        return false;
    }
    
    ParsedChain& chain = m_parsedChains[chainIndex];
    if (monomerIndex < 0 || monomerIndex >= chain.monomers.size()) {
        return false;
    }
    
    ParsedMonomer& monomer = chain.monomers[monomerIndex];
    
    // Check if the resolved symbol is valid
    if (!isValidAminoAcid(resolvedSymbol) && !isValidNucleotide(resolvedSymbol)) {
        return false;
    }
    
    // Update the monomer
    monomer.resolvedSymbol = resolvedSymbol;
    monomer.symbol = resolvedSymbol;
    monomer.isAmbiguous = false;
    monomer.isSupported = true;
    
    // Remove ambiguity issues for this monomer
    QList<ImportIssue> newIssues;
    for (const ImportIssue& issue : m_issues) {
        if (issue.location.recordIndex == chainIndex &&
            issue.location.line == monomer.location.line &&
            issue.location.column == monomer.location.column &&
            issue.code == "AMBIGUOUS_SYMBOL") {
            // Skip this issue
            continue;
        }
        newIssues.append(issue);
    }
    m_issues = newIssues;
    
    // Check if we need to update polymer type
    if (isValidNucleotide(resolvedSymbol)) {
        bool hasPeptide = false;
        for (const ParsedMonomer& m : chain.monomers) {
            if (isValidAminoAcid(m.symbol) || isValidAminoAcid(m.resolvedSymbol)) {
                hasPeptide = true;
                break;
            }
        }
        
        if (!hasPeptide) {
            chain.polymerType = PolymerType::DNA;
            chain.polymerTypeStr = "DNA";
        }
    }
    
    return true;
}

int BiopolymerImportSession::autoResolveAllAmbiguous() {
    int resolvedCount = 0;
    
    for (auto& chain : m_parsedChains) {
        for (auto& monomer : chain.monomers) {
            if (monomer.isAmbiguous && !monomer.isSupported) {
                QString resolved = resolveAmbiguousSymbol(monomer.originalSymbol, 
                                                         chain.polymerType);
                if (!resolved.isEmpty()) {
                    monomer.resolvedSymbol = resolved;
                    monomer.symbol = resolved;
                    monomer.isAmbiguous = false;
                    monomer.isSupported = true;
                    resolvedCount++;
                    
                    // Update polymer type if needed
                    if (isValidNucleotide(resolved)) {
                        bool hasPeptide = false;
                        for (const ParsedMonomer& m : chain.monomers) {
                            if (isValidAminoAcid(m.symbol) || isValidAminoAcid(m.resolvedSymbol)) {
                                hasPeptide = true;
                                break;
                            }
                        }
                        
                        if (!hasPeptide) {
                            chain.polymerType = PolymerType::DNA;
                            chain.polymerTypeStr = "DNA";
                        }
                    }
                }
            }
        }
    }
    
    // Rebuild issues
    m_issues.clear();
    m_parseSuccess = true;
    for (int i = 0; i < m_parsedChains.size(); i++) {
        QList<ImportIssue> chainIssues = validateChain(m_parsedChains[i], i);
        m_issues.append(chainIssues);
    }
    
    return resolvedCount;
}

// ========== Preview ==========

ImportPreview BiopolymerImportSession::generatePreview() const {
    ImportPreview preview;
    
    if (!m_isParsed || !m_parseSuccess) {
        preview.commitBlocked = true;
        return preview;
    }
    
    preview.chainCount = m_parsedChains.size();
    
    for (const ParsedChain& chain : m_parsedChains) {
        PreviewChain previewChain;
        previewChain.name = chain.name;
        previewChain.polymerType = chain.polymerType;
        previewChain.residueCount = chain.monomers.size();
        
        // Estimate atoms and bonds
        int atomCount = 0;
        int bondCount = 0;
        
        for (const ParsedMonomer& monomer : chain.monomers) {
            // Rough estimate: amino acids ~10-20 atoms, nucleotides ~20-30 atoms
            if (isValidAminoAcid(monomer.symbol) || isValidAminoAcid(monomer.resolvedSymbol)) {
                atomCount += 15;  // Average for amino acid
            } else if (isValidNucleotide(monomer.symbol) || isValidNucleotide(monomer.resolvedSymbol)) {
                atomCount += 25;  // Average for nucleotide
            } else {
                atomCount += 10;  // Unknown
            }
        }
        
        // Bonds: roughly atoms - 1 for backbone + side chain bonds
        bondCount = atomCount - chain.monomers.size();
        
        previewChain.estimatedAtomCount = atomCount;
        previewChain.hasErrors = false;
        previewChain.hasWarnings = false;
        
        // Check for issues
        for (const ImportIssue& issue : m_issues) {
            if (issue.location.recordId == chain.recordId) {
                if (issue.blocksCommit()) {
                    previewChain.hasErrors = true;
                } else {
                    previewChain.hasWarnings = true;
                }
            }
        }
        
        preview.chains.append(previewChain);
        preview.residueCount += previewChain.residueCount;
        preview.estimatedAtomCount += previewChain.estimatedAtomCount;
        preview.estimatedBondCount += bondCount;
    }
    
    // Check for blocking issues
    preview.commitBlocked = !canCommit();
    preview.blockingIssues = this->blockingIssues();
    preview.warnings = this->warnings();
    
    return preview;
}

PreviewChain BiopolymerImportSession::estimateChainPreview(const ParsedChain& chain) {
    PreviewChain previewChain;
    previewChain.name = chain.name;
    previewChain.polymerType = chain.polymerType;
    previewChain.residueCount = chain.monomers.size();
    
    int atomCount = 0;
    for (const ParsedMonomer& monomer : chain.monomers) {
        if (isValidAminoAcid(monomer.symbol) || isValidAminoAcid(monomer.resolvedSymbol)) {
            atomCount += 15;
        } else if (isValidNucleotide(monomer.symbol) || isValidNucleotide(monomer.resolvedSymbol)) {
            atomCount += 25;
        } else {
            atomCount += 10;
        }
    }
    
    previewChain.estimatedAtomCount = atomCount;
    previewChain.hasErrors = false;
    previewChain.hasWarnings = false;
    
    return previewChain;
}

// ========== Commit ==========

bool BiopolymerImportSession::commitToState(BiopolymerState& state, 
                                           EditableMolecule& molecule) const {
    return commitToState(state, molecule, nullptr);
}

bool BiopolymerImportSession::commitToState(BiopolymerState& state, 
                                           EditableMolecule& molecule,
                                           const ProgressCallback& progress) const {
    if (!canCommit()) {
        return false;
    }
    
    if (m_parsedChains.isEmpty()) {
        return false;
    }
    
    int totalSteps = m_parsedChains.size() * 2; // Create chains + create residues
    int currentStep = 0;
    
    // Clear existing biopolymer state if we're replacing
    // Note: This is handled by the caller for Replace mode
    
    for (const ParsedChain& chain : m_parsedChains) {
        if (progress) {
            progress(currentStep, totalSteps, QString("Creating chain '%1'").arg(chain.name));
        }
        
        int chainId = commitChain(chain, state, molecule, progress, currentStep);
        if (chainId < 0) {
            return false;
        }
    }
    
    // Set import provenance
    ImportProvenance provenance;
    provenance.sourceFormat = ImportFormatRegistry::formatName(m_format);
    provenance.sourceText = m_sourceText;
    provenance.importTime = QDateTime::currentDateTime();
    provenance.unsupportedFeatures = m_unsupportedFeatures;
    provenance.fullyMapped = true;
    
    state.setImportProvenance(provenance);
    
    return true;
}

int BiopolymerImportSession::commitChain(const ParsedChain& chain, 
                                        BiopolymerState& state, 
                                        EditableMolecule& molecule,
                                        const ProgressCallback& progress,
                                        int& currentStep) const {
    // Add chain to biopolymer state
    int chainId = state.addChain(chain.name, chain.polymerType);
    
    // Position tracking
    QList<QPointF> positions;
    double x = 0.0;
    double y = 0.0;
    
    // Create residues
    for (int i = 0; i < chain.monomers.size(); i++) {
        const ParsedMonomer& monomer = chain.monomers[i];
        
        if (progress) {
            currentStep++;
            progress(currentStep, -1, 
                    QString("Creating residue %1/%2 in chain '%3'")
                        .arg(i + 1).arg(chain.monomers.size()).arg(chain.name));
        }
        
        // Get the symbol to use
        QString symbol = monomer.resolvedSymbol.isEmpty() ? monomer.symbol : monomer.resolvedSymbol;
        
        // Create monomer identity
        MonomerIdentity identity;
        identity.symbol = symbol;
        identity.alias = IndigoMonomerAdapter::monomerAlias(symbol);
        identity.monomerClass = (chain.polymerType == PolymerType::Peptide) ? "AminoAcid" : "Base";
        identity.ambiguous = monomer.isAmbiguous;
        identity.stereochemistry = "L"; // Default for now
        
        // Add residue to chain
        int residueId = state.addResidue(chainId, identity, i + 1);
        
        // Create atoms for this residue
        QList<AtomId> atomIds = createMonomerAtoms(residueId, monomer, molecule, x, y);
        
        // Set residue atoms
        state.setResidueAtoms(residueId, atomIds);
        
        // Position for next residue
        x += 3.0; // Bond length
        positions.append(QPointF(x, y));
    }
    
    // Create connections between residues
    for (int i = 0; i < chain.monomers.size() - 1; i++) {
        int fromResidueId = state.residueIdsInChain(chainId)[i];
        int toResidueId = state.residueIdsInChain(chainId)[i + 1];
        
        // Find atoms to connect
        // For now, we'll just create a placeholder bond
        // In a real implementation, we'd find the correct atoms based on monomer templates
        
        if (!state.atomsOwnedByResidue(fromResidueId).isEmpty() && 
            !state.atomsOwnedByResidue(toResidueId).isEmpty()) {
            
            AtomId fromAtom = state.atomsOwnedByResidue(fromResidueId).first();
            AtomId toAtom = state.atomsOwnedByResidue(toResidueId).first();
            
            BondId bondId = molecule.addBond(fromAtom, toAtom, 1); // Single bond
            
            // Set connection in biopolymer state
            QString connectionType = (chain.polymerType == PolymerType::Peptide) ? "peptide" : "phosphodiester";
            state.setResidueConnection(fromResidueId, toResidueId, bondId, connectionType);
        }
    }
    
    return chainId;
}

QList<AtomId> BiopolymerImportSession::createMonomerAtoms(int residueId, 
                                                         const ParsedMonomer& monomer,
                                                         EditableMolecule& molecule,
                                                         double x, double y) const {
    QList<AtomId> atomIds;
    
    QString symbol = monomer.resolvedSymbol.isEmpty() ? monomer.symbol : monomer.resolvedSymbol;
    
    // Use IndigoMonomerAdapter to create atoms for this monomer
    // For now, we'll create a simplified representation
    
    // Add a central atom (alpha carbon or sugar carbon)
    AtomId centralAtom = molecule.addAtom(symbol, x, y);
    atomIds.append(centralAtom);
    
    // Add side chain atoms based on monomer type
    // This is a placeholder - real implementation would use monomer templates
    
    if (symbol == "Gly" || symbol == "A" || symbol == "G" || symbol == "C" || symbol == "T" || symbol == "U") {
        // Simple representation for glycine and nucleotides
        // For now, just add a few atoms in a simple pattern
        atomIds.append(molecule.addAtom("C", x + 1.0, y + 0.5));
        atomIds.append(molecule.addAtom("O", x + 1.0, y - 0.5));
        
        // Connect them
        molecule.addBond(centralAtom, atomIds[1], 1);
        molecule.addBond(centralAtom, atomIds[2], 1);
    } else {
        // Standard amino acid - add side chain
        atomIds.append(molecule.addAtom("C", x + 1.0, y + 0.5));
        atomIds.append(molecule.addAtom("C", x + 2.0, y + 0.5));
        atomIds.append(molecule.addAtom("O", x + 3.0, y + 0.5));
        
        // Connect them
        molecule.addBond(centralAtom, atomIds[1], 1);
        molecule.addBond(atomIds[1], atomIds[2], 1);
        molecule.addBond(atomIds[2], atomIds[3], 1);
    }
    
    return atomIds;
}

// ========== State Queries ==========

bool BiopolymerImportSession::hasValidData() const {
    return m_isParsed && m_parseSuccess && !m_parsedChains.isEmpty();
}

bool BiopolymerImportSession::hasErrors() const {
    for (const ImportIssue& issue : m_issues) {
        if (issue.severity == ImportIssueSeverity::Error) {
            return true;
        }
    }
    return false;
}

bool BiopolymerImportSession::hasWarnings() const {
    for (const ImportIssue& issue : m_issues) {
        if (issue.severity == ImportIssueSeverity::Warning) {
            return true;
        }
    }
    return false;
}

bool BiopolymerImportSession::hasUnsupportedFeatures() const {
    return !m_unsupportedFeatures.isEmpty();
}

QSet<QString> BiopolymerImportSession::unsupportedFeatures() const {
    return m_unsupportedFeatures;
}

// ========== Options ==========

ImportOptions BiopolymerImportSession::options() const {
    return m_options;
}

void BiopolymerImportSession::setOptions(const ImportOptions& options) {
    m_options = options;
}

// ========== Reset ==========

void BiopolymerImportSession::reset() {
    m_sourceText.clear();
    m_parsedChains.clear();
    m_issues.clear();
    m_unsupportedFeatures.clear();
    m_format = ImportFormat::Unknown;
    m_formatExplicit = false;
    m_isParsed = false;
    m_parseSuccess = false;
}

// ========== Helper Methods ==========

QList<ImportIssue> BiopolymerImportSession::blockingIssues() const {
    QList<ImportIssue> blocking;
    for (const ImportIssue& issue : m_issues) {
        if (issue.blocksCommit()) {
            blocking.append(issue);
        }
    }
    return blocking;
}

QList<ImportIssue> BiopolymerImportSession::warnings() const {
    QList<ImportIssue> result;
    for (const ImportIssue& issue : m_issues) {
        if (!issue.blocksCommit()) {
            result.append(issue);
        }
    }
    return result;
}

// ========== Static Validation Methods ==========

bool BiopolymerImportSession::isValidAminoAcid(const QString& symbol) {
    static const QSet<QString> validAminoAcids = {
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
        "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"
    };
    return validAminoAcids.contains(symbol.toUpper());
}

bool BiopolymerImportSession::isValidNucleotide(const QString& symbol) {
    static const QSet<QString> validNucleotides = {
        "A", "G", "C", "T", "U",
        "dA", "dG", "dC", "dT", "dU",
        "Ade", "Gua", "Cyt", "Thy", "Ura",
        "Adenine", "Guanine", "Cytosine", "Thymine", "Uracil"
    };
    return validNucleotides.contains(symbol.toUpper());
}

bool BiopolymerImportSession::isAmbiguousSymbol(const QString& symbol) {
    static const QSet<QString> ambiguousSymbols = {
        "B", "Z", "X", "N", "J", "O", "U"
    };
    return ambiguousSymbols.contains(symbol.toUpper());
}

QSet<QString> BiopolymerImportSession::validAminoAcidSymbols() {
    return {
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    };
}

QSet<QString> BiopolymerImportSession::validNucleotideSymbols() {
    return {"A", "G", "C", "T", "U"};
}

QSet<QString> BiopolymerImportSession::ambiguousSymbolNames() {
    return {"B", "Z", "X", "N", "J", "O", "U"};
}

QString BiopolymerImportSession::resolveAmbiguousSymbol(const QString& symbol, 
                                                        PolymerType polymerType) {
    // Default resolution rules
    // For peptides: X -> Ala (most common), B -> Asx (Asn/Asp), Z -> Glx (Gln/Glu)
    // For nucleotides: N -> A (adenine), etc.
    
    QString upper = symbol.toUpper();
    
    if (polymerType == PolymerType::Peptide || polymerType == PolymerType::Unknown) {
        if (upper == "X") return "Ala";
        if (upper == "B") return "Asn"; // Asn or Asp
        if (upper == "Z") return "Gln"; // Gln or Glu
    } else if (polymerType == PolymerType::DNA || polymerType == PolymerType::RNA) {
        if (upper == "N") return "A";
        if (upper == "X") return "A";
        if (upper == "B") return "G"; // In IUPAC, B = C/G/T (not A)
    }
    
    return ""; // Cannot resolve
}

QString BiopolymerImportSession::normalizeSequence(const QString& sequence) {
    QString result;
    for (QChar c : sequence) {
        if (!c.isSpace() && c != '-' && c != '_' && c != '.' && c != ',') {
            result.append(c.toUpper());
        }
    }
    return result;
}

// ========== Format Registry ==========

QMap<ImportFormat, ImportFormatRegistry::FormatInfo>& ImportFormatRegistry::registry() {
    static QMap<ImportFormat, FormatInfo> instance;
    static bool initialized = false;
    
    if (!initialized) {
        // Register default parsers
        registerParser(ImportFormat::FASTA, "FASTA", 
                      [](const QString& text, const ImportOptions& options) {
                          return BiopolymerImportSession::parseFasta(text, options);
                      });
        registerParser(ImportFormat::Sequence, "Sequence",
                      [](const QString& text, const ImportOptions& options) {
                          return BiopolymerImportSession::parseSequence(text, options);
                      });
        registerParser(ImportFormat::Simple, "Simple",
                      [](const QString& text, const ImportOptions& options) {
                          return BiopolymerImportSession::parseSimple(text, options);
                      });
        
        initialized = true;
    }
    
    return instance;
}

void ImportFormatRegistry::registerParser(ImportFormat format, const QString& name,
                                         const FormatParser& parser) {
    registry()[format] = {name, parser};
}

ImportFormatRegistry::FormatParser ImportFormatRegistry::getParser(ImportFormat format) {
    auto it = registry().find(format);
    if (it != registry().end()) {
        return it->parser;
    }
    return nullptr;
}

QList<ImportFormat> ImportFormatRegistry::allFormats() {
    return registry().keys();
}

QString ImportFormatRegistry::formatName(ImportFormat format) {
    auto it = registry().find(format);
    if (it != registry().end()) {
        return it->name;
    }
    return "Unknown";
}

ImportFormat ImportFormatRegistry::formatByName(const QString& name) {
    for (auto it = registry().begin(); it != registry().end(); ++it) {
        if (it->name == name) {
            return it.key();
        }
    }
    return ImportFormat::Unknown;
}

// ========== Serialization Support ==========

QMap<QString, QVariant> BiopolymerImportSession::toJson() const {
    QMap<QString, QVariant> json;
    json["sourceText"] = m_sourceText;
    json["format"] = static_cast<int>(m_format);
    json["formatExplicit"] = m_formatExplicit;
    json["isParsed"] = m_isParsed;
    json["parseSuccess"] = m_parseSuccess;
    
    // Serialize options
    QMap<QString, QVariant> optionsJson;
    optionsJson["autoResolveAmbiguous"] = m_options.autoResolveAmbiguous;
    optionsJson["defaultPolymerType"] = static_cast<int>(m_options.defaultPolymerType);
    optionsJson["createMultipleChains"] = m_options.createMultipleChains;
    optionsJson["preserveRecordBoundaries"] = m_options.preserveRecordBoundaries;
    optionsJson["validateStereochemistry"] = m_options.validateStereochemistry;
    optionsJson["detectUnsupportedFeatures"] = m_options.detectUnsupportedFeatures;
    optionsJson["allowPartialImport"] = m_options.allowPartialImport;
    json["options"] = optionsJson;
    
    // Serialize parsed chains
    QList<QVariant> chainsJson;
    for (const ParsedChain& chain : m_parsedChains) {
        QMap<QString, QVariant> chainJson;
        chainJson["name"] = chain.name;
        chainJson["originalName"] = chain.originalName;
        chainJson["polymerType"] = static_cast<int>(chain.polymerType);
        chainJson["polymerTypeStr"] = chain.polymerTypeStr;
        chainJson["recordIndex"] = chain.recordIndex;
        chainJson["recordId"] = chain.recordId;
        chainJson["description"] = chain.description;
        chainJson["isValid"] = chain.isValid;
        chainJson["invalidityReason"] = chain.invalidityReason;
        
        // Serialize header location
        QMap<QString, QVariant> headerLocJson;
        headerLocJson["line"] = chain.headerLocation.line;
        headerLocJson["column"] = chain.headerLocation.column;
        headerLocJson["recordIndex"] = chain.headerLocation.recordIndex;
        headerLocJson["recordId"] = chain.headerLocation.recordId;
        chainJson["headerLocation"] = headerLocJson;
        
        // Serialize body location
        QMap<QString, QVariant> bodyLocJson;
        bodyLocJson["line"] = chain.bodyLocation.line;
        bodyLocJson["column"] = chain.bodyLocation.column;
        bodyLocJson["recordIndex"] = chain.bodyLocation.recordIndex;
        bodyLocJson["recordId"] = chain.bodyLocation.recordId;
        chainJson["bodyLocation"] = bodyLocJson;
        
        // Serialize monomers
        QList<QVariant> monomersJson;
        for (const ParsedMonomer& monomer : chain.monomers) {
            QMap<QString, QVariant> monomerJson;
            monomerJson["symbol"] = monomer.symbol;
            monomerJson["originalSymbol"] = monomer.originalSymbol;
            monomerJson["position"] = monomer.position;
            monomerJson["isAmbiguous"] = monomer.isAmbiguous;
            monomerJson["isSupported"] = monomer.isSupported;
            monomerJson["resolvedSymbol"] = monomer.resolvedSymbol;
            
            // Serialize location
            QMap<QString, QVariant> locJson;
            locJson["line"] = monomer.location.line;
            locJson["column"] = monomer.location.column;
            locJson["recordIndex"] = monomer.location.recordIndex;
            locJson["recordId"] = monomer.location.recordId;
            monomerJson["location"] = locJson;
            
            monomersJson.append(monomerJson);
        }
        chainJson["monomers"] = monomersJson;
        
        chainsJson.append(chainJson);
    }
    json["chains"] = chainsJson;
    
    // Serialize issues
    QList<QVariant> issuesJson;
    for (const ImportIssue& issue : m_issues) {
        QMap<QString, QVariant> issueJson;
        issueJson["severity"] = static_cast<int>(issue.severity);
        issueJson["code"] = issue.code;
        issueJson["message"] = issue.message;
        issueJson["suggestion"] = issue.suggestion;
        
        QMap<QString, QVariant> locJson;
        locJson["line"] = issue.location.line;
        locJson["column"] = issue.location.column;
        locJson["recordIndex"] = issue.location.recordIndex;
        locJson["recordId"] = issue.location.recordId;
        issueJson["location"] = locJson;
        
        issuesJson.append(issueJson);
    }
    json["issues"] = issuesJson;
    
    // Serialize unsupported features
    QList<QVariant> unsupportedJson;
    for (const QString& feature : m_unsupportedFeatures) {
        unsupportedJson.append(feature);
    }
    json["unsupportedFeatures"] = unsupportedJson;
    
    return json;
}

BiopolymerImportSession BiopolymerImportSession::fromJson(const QMap<QString, QVariant>& json) {
    BiopolymerImportSession session;
    
    session.m_sourceText = json.value("sourceText").toString();
    session.m_format = static_cast<ImportFormat>(json.value("format").toInt());
    session.m_formatExplicit = json.value("formatExplicit").toBool();
    session.m_isParsed = json.value("isParsed").toBool();
    session.m_parseSuccess = json.value("parseSuccess").toBool();
    
    // Deserialize options
    QMap<QString, QVariant> optionsJson = json.value("options").toMap();
    session.m_options.autoResolveAmbiguous = optionsJson.value("autoResolveAmbiguous").toBool();
    session.m_options.defaultPolymerType = static_cast<PolymerType>(optionsJson.value("defaultPolymerType").toInt());
    session.m_options.createMultipleChains = optionsJson.value("createMultipleChains").toBool();
    session.m_options.preserveRecordBoundaries = optionsJson.value("preserveRecordBoundaries").toBool();
    session.m_options.validateStereochemistry = optionsJson.value("validateStereochemistry").toBool();
    session.m_options.detectUnsupportedFeatures = optionsJson.value("detectUnsupportedFeatures").toBool();
    session.m_options.allowPartialImport = optionsJson.value("allowPartialImport").toBool();
    
    // Deserialize chains
    QList<QVariant> chainsJson = json.value("chains").toList();
    for (const QVariant& chainVar : chainsJson) {
        QMap<QString, QVariant> chainJson = chainVar.toMap();
        ParsedChain chain;
        
        chain.name = chainJson.value("name").toString();
        chain.originalName = chainJson.value("originalName").toString();
        chain.polymerType = static_cast<PolymerType>(chainJson.value("polymerType").toInt());
        chain.polymerTypeStr = chainJson.value("polymerTypeStr").toString();
        chain.recordIndex = chainJson.value("recordIndex").toInt();
        chain.recordId = chainJson.value("recordId").toString();
        chain.description = chainJson.value("description").toString();
        chain.isValid = chainJson.value("isValid").toBool();
        chain.invalidityReason = chainJson.value("invalidityReason").toString();
        
        // Deserialize header location
        QMap<QString, QVariant> headerLocJson = chainJson.value("headerLocation").toMap();
        chain.headerLocation.line = headerLocJson.value("line").toInt();
        chain.headerLocation.column = headerLocJson.value("column").toInt();
        chain.headerLocation.recordIndex = headerLocJson.value("recordIndex").toInt();
        chain.headerLocation.recordId = headerLocJson.value("recordId").toString();
        
        // Deserialize body location
        QMap<QString, QVariant> bodyLocJson = chainJson.value("bodyLocation").toMap();
        chain.bodyLocation.line = bodyLocJson.value("line").toInt();
        chain.bodyLocation.column = bodyLocJson.value("column").toInt();
        chain.bodyLocation.recordIndex = bodyLocJson.value("recordIndex").toInt();
        chain.bodyLocation.recordId = bodyLocJson.value("recordId").toString();
        
        // Deserialize monomers
        QList<QVariant> monomersJson = chainJson.value("monomers").toList();
        for (const QVariant& monomerVar : monomersJson) {
            QMap<QString, QVariant> monomerJson = monomerVar.toMap();
            ParsedMonomer monomer;
            
            monomer.symbol = monomerJson.value("symbol").toString();
            monomer.originalSymbol = monomerJson.value("originalSymbol").toString();
            monomer.position = monomerJson.value("position").toInt();
            monomer.isAmbiguous = monomerJson.value("isAmbiguous").toBool();
            monomer.isSupported = monomerJson.value("isSupported").toBool();
            monomer.resolvedSymbol = monomerJson.value("resolvedSymbol").toString();
            
            // Deserialize location
            QMap<QString, QVariant> locJson = monomerJson.value("location").toMap();
            monomer.location.line = locJson.value("line").toInt();
            monomer.location.column = locJson.value("column").toInt();
            monomer.location.recordIndex = locJson.value("recordIndex").toInt();
            monomer.location.recordId = locJson.value("recordId").toString();
            
            chain.monomers.append(monomer);
        }
        
        session.m_parsedChains.append(chain);
    }
    
    // Deserialize issues
    QList<QVariant> issuesJson = json.value("issues").toList();
    for (const QVariant& issueVar : issuesJson) {
        QMap<QString, QVariant> issueJson = issueVar.toMap();
        ImportIssue issue;
        
        issue.severity = static_cast<ImportIssueSeverity>(issueJson.value("severity").toInt());
        issue.code = issueJson.value("code").toString();
        issue.message = issueJson.value("message").toString();
        issue.suggestion = issueJson.value("suggestion").toString();
        
        QMap<QString, QVariant> locJson = issueJson.value("location").toMap();
        issue.location.line = locJson.value("line").toInt();
        issue.location.column = locJson.value("column").toInt();
        issue.location.recordIndex = locJson.value("recordIndex").toInt();
        issue.location.recordId = locJson.value("recordId").toString();
        
        session.m_issues.append(issue);
    }
    
    // Deserialize unsupported features
    QList<QVariant> unsupportedJson = json.value("unsupportedFeatures").toList();
    for (const QVariant& featureVar : unsupportedJson) {
        session.m_unsupportedFeatures.insert(featureVar.toString());
    }
    
    return session;
}
