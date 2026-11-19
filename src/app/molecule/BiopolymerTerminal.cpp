// src/app/molecule/BiopolymerTerminal.cpp
#include "BiopolymerTerminal.h"

// ========== Terminal Type Helpers ==========

QString terminalTypeToString(TerminalType type) {
    switch (type) {
        case TerminalType::None: return "None";
        case TerminalType::N_Terminal: return "N_Terminal";
        case TerminalType::C_Terminal: return "C_Terminal";
        case TerminalType::Phosphate_5: return "5'_Phosphate";
        case TerminalType::Phosphate_3: return "3'_Phosphate";
        case TerminalType::Hydroxyl_5: return "5'_Hydroxyl";
        case TerminalType::Hydroxyl_3: return "3'_Hydroxyl";
        case TerminalType::Amine: return "Amine";
        case TerminalType::Acetyl: return "Acetyl";
        case TerminalType::Amide: return "Amide";
        case TerminalType::Carboxyl: return "Carboxyl";
        case TerminalType::Methyl: return "Methyl";
        case TerminalType::Hydroxyl: return "Hydroxyl";
        case TerminalType::Phosphonate: return "Phosphonate";
        case TerminalType::Biotin: return "Biotin";
        case TerminalType::Fluorescent: return "Fluorescent";
        case TerminalType::Custom: return "Custom";
        case TerminalType::Unknown: return "Unknown";
    }
    return "Unknown";
}

TerminalType terminalTypeFromString(const QString& typeString) {
    static const QMap<QString, TerminalType> typeMap = {
        {"none", TerminalType::None},
        {"n_terminal", TerminalType::N_Terminal},
        {"n-term", TerminalType::N_Terminal},
        {"c_terminal", TerminalType::C_Terminal},
        {"c-term", TerminalType::C_Terminal},
        {"5'_phosphate", TerminalType::Phosphate_5},
        {"5p", TerminalType::Phosphate_5},
        {"phosphate_5", TerminalType::Phosphate_5},
        {"3'_phosphate", TerminalType::Phosphate_3},
        {"3p", TerminalType::Phosphate_3},
        {"phosphate_3", TerminalType::Phosphate_3},
        {"5'_hydroxyl", TerminalType::Hydroxyl_5},
        {"5oh", TerminalType::Hydroxyl_5},
        {"3'_hydroxyl", TerminalType::Hydroxyl_3},
        {"3oh", TerminalType::Hydroxyl_3},
        {"amine", TerminalType::Amine},
        {"acetyl", TerminalType::Acetyl},
        {"ac", TerminalType::Acetyl},
        {"amide", TerminalType::Amide},
        {"carboxyl", TerminalType::Carboxyl},
        {"cooh", TerminalType::Carboxyl},
        {"methyl", TerminalType::Methyl},
        {"hydroxyl", TerminalType::Hydroxyl},
        {"phosphonate", TerminalType::Phosphonate},
        {"biotin", TerminalType::Biotin},
        {"fluorescent", TerminalType::Fluorescent},
        {"custom", TerminalType::Custom},
        {"unknown", TerminalType::Unknown}
    };
    
    QString lower = typeString.toLower();
    auto it = typeMap.find(lower);
    if (it != typeMap.end()) {
        return it.value();
    }
    return TerminalType::Unknown;
}

QList<TerminalType> allTerminalTypes() {
    return {
        TerminalType::None,
        TerminalType::N_Terminal,
        TerminalType::C_Terminal,
        TerminalType::Phosphate_5,
        TerminalType::Phosphate_3,
        TerminalType::Hydroxyl_5,
        TerminalType::Hydroxyl_3,
        TerminalType::Amine,
        TerminalType::Acetyl,
        TerminalType::Amide,
        TerminalType::Carboxyl,
        TerminalType::Methyl,
        TerminalType::Hydroxyl,
        TerminalType::Phosphonate,
        TerminalType::Biotin,
        TerminalType::Fluorescent,
        TerminalType::Custom,
        TerminalType::Unknown
    };
}

bool isNTerminal(TerminalType type) {
    return type == TerminalType::N_Terminal ||
           type == TerminalType::Amine ||
           type == TerminalType::Acetyl;
}

bool isCTerminal(TerminalType type) {
    return type == TerminalType::C_Terminal ||
           type == TerminalType::Amide ||
           type == TerminalType::Carboxyl;
}

bool is5PrimeTerminal(TerminalType type) {
    return type == TerminalType::Phosphate_5 ||
           type == TerminalType::Hydroxyl_5;
}

bool is3PrimeTerminal(TerminalType type) {
    return type == TerminalType::Phosphate_3 ||
           type == TerminalType::Hydroxyl_3;
}

// ========== Terminal Position Helpers ==========

QString terminalPositionToString(TerminalPosition position) {
    switch (position) {
        case TerminalPosition::First: return "First";
        case TerminalPosition::Last: return "Last";
        case TerminalPosition::Both: return "Both";
        case TerminalPosition::Internal: return "Internal";
        case TerminalPosition::Unknown: return "Unknown";
    }
    return "Unknown";
}

TerminalPosition terminalPositionFromString(const QString& positionString) {
    static const QMap<QString, TerminalPosition> positionMap = {
        {"first", TerminalPosition::First},
        {"last", TerminalPosition::Last},
        {"both", TerminalPosition::Both},
        {"internal", TerminalPosition::Internal},
        {"unknown", TerminalPosition::Unknown}
    };
    
    QString lower = positionString.toLower();
    auto it = positionMap.find(lower);
    if (it != positionMap.end()) {
        return it.value();
    }
    return TerminalPosition::Unknown;
}

// ========== BiopolymerTerminal Implementation ==========

BiopolymerTerminal::BiopolymerTerminal() {
    // Default constructor
}

BiopolymerTerminal::BiopolymerTerminal(int termId, int chId, TerminalType termType,
                                       TerminalPosition pos, const QString& termName)
    : terminalId(termId), chainId(chId), type(termType), position(pos), name(termName) {
    // Initialize based on type
    if (termType == TerminalType::Acetyl) {
        shortName = "Ac";
        chemicalFormula = "C2H3O";
        molecularWeight = 43.05;
        description = "Acetyl group";
    } else if (termType == TerminalType::Amine) {
        shortName = "NH2";
        chemicalFormula = "NH2";
        molecularWeight = 16.0;
        description = "Amine group";
    } else if (termType == TerminalType::Amide) {
        shortName = "CONH2";
        chemicalFormula = "CONH2";
        molecularWeight = 44.03;
        description = "Amide group";
    } else if (termType == TerminalType::Carboxyl) {
        shortName = "COOH";
        chemicalFormula = "CO2H";
        molecularWeight = 44.01;
        description = "Carboxyl group";
    } else if (termType == TerminalType::Phosphate_5) {
        shortName = "5p";
        chemicalFormula = "PO3H2";
        molecularWeight = 94.97;
        description = "5' phosphate group";
    } else if (termType == TerminalType::Phosphate_3) {
        shortName = "3p";
        chemicalFormula = "PO3H2";
        molecularWeight = 94.97;
        description = "3' phosphate group";
    } else if (termType == TerminalType::Hydroxyl_5) {
        shortName = "5OH";
        chemicalFormula = "OH";
        molecularWeight = 17.01;
        description = "5' hydroxyl group";
    } else if (termType == TerminalType::Hydroxyl_3) {
        shortName = "3OH";
        chemicalFormula = "OH";
        molecularWeight = 17.01;
        description = "3' hydroxyl group";
    }
}

bool BiopolymerTerminal::isValid() const {
    return terminalId >= 0 && chainId >= 0 && type != TerminalType::None;
}

bool BiopolymerTerminal::isNTerminal() const {
    return ::isNTerminal(type) || position == TerminalPosition::First;
}

bool BiopolymerTerminal::isCTerminal() const {
    return ::isCTerminal(type) || position == TerminalPosition::Last;
}

bool BiopolymerTerminal::is5Prime() const {
    return ::is5PrimeTerminal(type) || position == TerminalPosition::First;
}

bool BiopolymerTerminal::is3Prime() const {
    return ::is3PrimeTerminal(type) || position == TerminalPosition::Last;
}

QString BiopolymerTerminal::fullName() const {
    if (name.isEmpty()) {
        return terminalTypeToString(type);
    }
    return name;
}

QString BiopolymerTerminal::displayName() const {
    if (!shortName.isEmpty()) {
        return shortName;
    }
    return fullName();
}

// ========== Static Factory Methods ==========

BiopolymerTerminal BiopolymerTerminal::createNTerminalAcetyl(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Acetyl, TerminalPosition::First, "N-Terminal Acetyl");
    term.residueId = residueId;
    term.shortName = "Ac-N";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::createNTerminalAmine(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Amine, TerminalPosition::First, "N-Terminal Amine");
    term.residueId = residueId;
    term.shortName = "NH2-N";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::createCTerminalAmide(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Amide, TerminalPosition::Last, "C-Terminal Amide");
    term.residueId = residueId;
    term.shortName = "CONH2-C";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::createCTerminalCarboxyl(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Carboxyl, TerminalPosition::Last, "C-Terminal Carboxyl");
    term.residueId = residueId;
    term.shortName = "COOH-C";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::create5PrimePhosphate(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Phosphate_5, TerminalPosition::First, "5' Phosphate");
    term.residueId = residueId;
    term.shortName = "5p";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::create5PrimeHydroxyl(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Hydroxyl_5, TerminalPosition::First, "5' Hydroxyl");
    term.residueId = residueId;
    term.shortName = "5OH";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::create3PrimePhosphate(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Phosphate_3, TerminalPosition::Last, "3' Phosphate");
    term.residueId = residueId;
    term.shortName = "3p";
    term.isSupported = true;
    return term;
}

BiopolymerTerminal BiopolymerTerminal::create3PrimeHydroxyl(int termId, int chainId, int residueId) {
    BiopolymerTerminal term(termId, chainId, TerminalType::Hydroxyl_3, TerminalPosition::Last, "3' Hydroxyl");
    term.residueId = residueId;
    term.shortName = "3OH";
    term.isSupported = true;
    return term;
}

bool BiopolymerTerminal::operator==(const BiopolymerTerminal& other) const {
    return terminalId == other.terminalId &&
           chainId == other.chainId &&
           type == other.type &&
           position == other.position &&
           name == other.name;
}

bool BiopolymerTerminal::operator!=(const BiopolymerTerminal& other) const {
    return !(*this == other);
}

// ========== Terminal Database Static Initialization ==========

// Peptide N-terminals
const TerminalDatabase::TerminalRecord TerminalDatabase::N_Terminal_Amine = {
    TerminalType::Amine,
    "N-Terminal Amine",
    "NH2",
    "NH2",
    16.0,
    TerminalPosition::First,
    {PolymerType::PEPTIDE},
    true,
    "Free amine group at N-terminus"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::N_Terminal_Acetyl = {
    TerminalType::Acetyl,
    "N-Terminal Acetyl",
    "Ac",
    "C2H3O",
    43.05,
    TerminalPosition::First,
    {PolymerType::PEPTIDE},
    true,
    "Acetylated N-terminus"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::N_Terminal_Methyl = {
    TerminalType::Methyl,
    "N-Terminal Methyl",
    "Me",
    "CH3",
    15.04,
    TerminalPosition::First,
    {PolymerType::PEPTIDE},
    true,
    "Methylated N-terminus"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::N_Terminal_Formyl = {
    TerminalType::Acetyl,  // Formyl is similar to acetyl
    "N-Terminal Formyl",
    "For",
    "CHO",
    29.02,
    TerminalPosition::First,
    {PolymerType::PEPTIDE},
    true,
    "Formylated N-terminus"
};

// Peptide C-terminals
const TerminalDatabase::TerminalRecord TerminalDatabase::C_Terminal_Carboxyl = {
    TerminalType::Carboxyl,
    "C-Terminal Carboxyl",
    "COOH",
    "CO2H",
    44.01,
    TerminalPosition::Last,
    {PolymerType::PEPTIDE},
    true,
    "Free carboxyl group at C-terminus"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::C_Terminal_Amide = {
    TerminalType::Amide,
    "C-Terminal Amide",
    "NH2",
    "CONH2",
    44.03,
    TerminalPosition::Last,
    {PolymerType::PEPTIDE},
    true,
    "Amide group at C-terminus"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::C_Terminal_MethylEster = {
    TerminalType::Methyl,
    "C-Terminal Methyl Ester",
    "OMe",
    "CO2CH3",
    59.04,
    TerminalPosition::Last,
    {PolymerType::PEPTIDE},
    true,
    "Methyl ester at C-terminus"
};

// DNA/RNA 5' terminals
const TerminalDatabase::TerminalRecord TerminalDatabase::DNA_5Prime_Phosphate = {
    TerminalType::Phosphate_5,
    "5' Phosphate",
    "5p",
    "PO3H2",
    94.97,
    TerminalPosition::First,
    {PolymerType::DNA, PolymerType::RNA},
    true,
    "Phosphate group at 5' end"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::DNA_5Prime_Hydroxyl = {
    TerminalType::Hydroxyl_5,
    "5' Hydroxyl",
    "5OH",
    "OH",
    17.01,
    TerminalPosition::First,
    {PolymerType::DNA, PolymerType::RNA},
    true,
    "Hydroxyl group at 5' end"
};

// DNA/RNA 3' terminals
const TerminalDatabase::TerminalRecord TerminalDatabase::DNA_3Prime_Phosphate = {
    TerminalType::Phosphate_3,
    "3' Phosphate",
    "3p",
    "PO3H2",
    94.97,
    TerminalPosition::Last,
    {PolymerType::DNA, PolymerType::RNA},
    true,
    "Phosphate group at 3' end"
};

const TerminalDatabase::TerminalRecord TerminalDatabase::DNA_3Prime_Hydroxyl = {
    TerminalType::Hydroxyl_3,
    "3' Hydroxyl",
    "3OH",
    "OH",
    17.01,
    TerminalPosition::Last,
    {PolymerType::DNA, PolymerType::RNA},
    true,
    "Hydroxyl group at 3' end"
};

// ========== Terminal Database Implementation ==========

TerminalDatabase& TerminalDatabase::instance() {
    static TerminalDatabase db;
    return db;
}

TerminalDatabase::TerminalDatabase() {
    initialize();
}

TerminalDatabase::~TerminalDatabase() {
    // Cleanup
}

void TerminalDatabase::initialize() {
    // Register common terminals
    registerCustomTerminal(N_Terminal_Amine);
    registerCustomTerminal(N_Terminal_Acetyl);
    registerCustomTerminal(N_Terminal_Methyl);
    registerCustomTerminal(N_Terminal_Formyl);
    registerCustomTerminal(C_Terminal_Carboxyl);
    registerCustomTerminal(C_Terminal_Amide);
    registerCustomTerminal(C_Terminal_MethylEster);
    registerCustomTerminal(DNA_5Prime_Phosphate);
    registerCustomTerminal(DNA_5Prime_Hydroxyl);
    registerCustomTerminal(DNA_3Prime_Phosphate);
    registerCustomTerminal(DNA_3Prime_Hydroxyl);
}

TerminalDatabase::TerminalRecord TerminalDatabase::getTerminal(const QString& name) const {
    auto it = m_terminalsByName.find(name);
    if (it != m_terminalsByName.end()) {
        return it.value();
    }
    return TerminalRecord();
}

TerminalDatabase::TerminalRecord TerminalDatabase::getTerminal(TerminalType type) const {
    auto it = m_terminalsByType.find(type);
    if (it != m_terminalsByType.end()) {
        return it.value();
    }
    return TerminalRecord();
}

QList<TerminalDatabase::TerminalRecord> TerminalDatabase::allTerminals() const {
    QList<TerminalRecord> all;
    for (const auto& entry : m_terminalsByName) {
        all.append(entry);
    }
    return all;
}

QList<TerminalDatabase::TerminalRecord> TerminalDatabase::terminalsForPolymerType(
    PolymerType polymerType) const {
    QList<TerminalRecord> result;
    for (const auto& entry : m_terminalsByName) {
        const TerminalRecord& record = entry;
        if (record.applicableTo.contains(polymerType) || record.applicableTo.isEmpty()) {
            result.append(record);
        }
    }
    return result;
}

QList<TerminalDatabase::TerminalRecord> TerminalDatabase::nTerminalsForPolymerType(
    PolymerType polymerType) const {
    QList<TerminalRecord> result;
    for (const auto& entry : m_terminalsByName) {
        const TerminalRecord& record = entry;
        if ((record.position == TerminalPosition::First) &&
            (record.applicableTo.contains(polymerType) || record.applicableTo.isEmpty())) {
            result.append(record);
        }
    }
    return result;
}

QList<TerminalDatabase::TerminalRecord> TerminalDatabase::cTerminalsForPolymerType(
    PolymerType polymerType) const {
    QList<TerminalRecord> result;
    for (const auto& entry : m_terminalsByName) {
        const TerminalRecord& record = entry;
        if ((record.position == TerminalPosition::Last) &&
            (record.applicableTo.contains(polymerType) || record.applicableTo.isEmpty())) {
            result.append(record);
        }
    }
    return result;
}

bool TerminalDatabase::isKnownTerminal(const QString& name) const {
    return m_terminalsByName.contains(name);
}

void TerminalDatabase::registerCustomTerminal(const TerminalRecord& record) {
    if (!record.name.isEmpty()) {
        m_terminalsByName[record.name] = record;
        m_terminalsByType[record.type] = record;
    }
}

// ========== TerminalValidator Implementation ==========

TerminalValidator::TerminalValidator(const BiopolymerState& state)
    : m_state(state) {
    // Initialize validator
}

bool TerminalValidator::validateTerminal(const BiopolymerTerminal& terminal) const {
    // Check if terminal is valid
    if (!terminal.isValid()) {
        return false;
    }
    
    // Check chain exists
    const BiopolymerChain* chain = m_state.chain(terminal.chainId);
    if (!chain) {
        return false;
    }
    
    // Check compatibility
    if (!isCompatible(terminal.chainId, terminal)) {
        return false;
    }
    
    // Check for conflicts
    if (wouldConflict(terminal.chainId, terminal)) {
        return false;
    }
    
    return true;
}

QList<QString> TerminalValidator::validateAllTerminals() const {
    QList<QString> issues;
    
    for (const BiopolymerTerminal& term : m_state.terminals()) {
        if (!validateTerminal(term)) {
            issues.append(QString("Invalid terminal %1 on chain %2")
                .arg(term.displayName()).arg(term.chainId));
        }
    }
    
    return issues;
}

bool TerminalValidator::isCompatible(int chainId, const BiopolymerTerminal& terminal) const {
    // Check chain type
    if (!checkChainType(chainId, terminal.type)) {
        return false;
    }
    
    // Check position
    if (!checkPosition(chainId, terminal.position)) {
        return false;
    }
    
    return true;
}

bool TerminalValidator::canHaveTerminalType(int chainId, TerminalType type, 
                                            TerminalPosition position) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) {
        return false;
    }
    
    // Check if type is applicable to polymer type
    TerminalDatabase& db = TerminalDatabase::instance();
    QList<TerminalDatabase::TerminalRecord> records = 
        db.terminalsForPolymerType(chain->polymerType);
    
    for (const auto& record : records) {
        if (record.type == type) {
            return true;
        }
    }
    
    return false;
}

bool TerminalValidator::wouldConflict(int chainId, const BiopolymerTerminal& terminal) const {
    return !getConflictingTerminals(chainId, terminal).isEmpty();
}

QList<int> TerminalValidator::getConflictingTerminals(int chainId,
                                                      const BiopolymerTerminal& terminal) const {
    QList<int> conflicts;
    
    // Check for terminals at the same position
    for (const BiopolymerTerminal& existing : m_state.terminals()) {
        if (existing.chainId == chainId && existing.position == terminal.position) {
            if (existing.terminalId != terminal.terminalId) {
                conflicts.append(existing.terminalId);
            }
        }
    }
    
    return conflicts;
}

QList<BiopolymerTerminal> TerminalValidator::getSuggestedTerminals(int chainId, 
                                                                    TerminalPosition position) const {
    QList<BiopolymerTerminal> suggestions;
    
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) {
        return suggestions;
    }
    
    // Get applicable terminals from database
    TerminalDatabase& db = TerminalDatabase::instance();
    QList<TerminalDatabase::TerminalRecord> records;
    
    if (position == TerminalPosition::First) {
        records = db.nTerminalsForPolymerType(chain->polymerType);
    } else if (position == TerminalPosition::Last) {
        records = db.cTerminalsForPolymerType(chain->polymerType);
    } else {
        records = db.terminalsForPolymerType(chain->polymerType);
    }
    
    // Create suggestions
    int nextId = 1;  // Will be replaced with real ID when added
    for (const auto& record : records) {
        if (record.isSupported) {
            BiopolymerTerminal term(nextId++, chainId, record.type, position, record.name);
            term.shortName = record.shortName;
            term.chemicalFormula = record.chemicalFormula;
            term.molecularWeight = record.molecularWeight;
            term.description = record.description;
            term.isSupported = record.isSupported;
            suggestions.append(term);
        }
    }
    
    return suggestions;
}

bool TerminalValidator::checkChainType(int chainId, TerminalType type) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) {
        return false;
    }
    
    TerminalDatabase& db = TerminalDatabase::instance();
    QList<TerminalDatabase::TerminalRecord> records = 
        db.terminalsForPolymerType(chain->polymerType);
    
    for (const auto& record : records) {
        if (record.type == type) {
            return true;
        }
    }
    
    return false;
}

bool TerminalValidator::checkPosition(int chainId, TerminalPosition position) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) {
        return false;
    }
    
    // For terminal positions, we need at least one residue
    if (chain->residueIds.isEmpty()) {
        return false;
    }
    
    // First position requires at least one residue
    if (position == TerminalPosition::First) {
        return true;
    }
    
    // Last position requires at least one residue
    if (position == TerminalPosition::Last) {
        return true;
    }
    
    return true;
}

// ========== TerminalManager Implementation ==========

TerminalManager::TerminalManager(BiopolymerState& state)
    : m_state(state) {
    // Initialize manager
}

bool TerminalManager::addTerminal(const BiopolymerTerminal& terminal) {
    // Validate the terminal
    TerminalValidator validator(m_state);
    if (!validator.validateTerminal(terminal)) {
        return false;
    }
    
    // Add to state
    return m_state.addTerminal(terminal);
}

int TerminalManager::addTerminal(int chainId, TerminalType type, TerminalPosition position) {
    // Create a new terminal
    int termId = nextTerminalId();
    BiopolymerTerminal term(termId, chainId, type, position, "");
    
    // Get details from database
    TerminalDatabase& db = TerminalDatabase::instance();
    TerminalDatabase::TerminalRecord record = db.getTerminal(type);
    if (record.type != TerminalType::None) {
        term.name = record.name;
        term.shortName = record.shortName;
        term.chemicalFormula = record.chemicalFormula;
        term.molecularWeight = record.molecularWeight;
        term.description = record.description;
        term.isSupported = record.isSupported;
    }
    
    // Set residue ID based on position
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (chain && !chain->residueIds.isEmpty()) {
        if (position == TerminalPosition::First) {
            term.residueId = chain->residueIds.first();
        } else if (position == TerminalPosition::Last) {
            term.residueId = chain->residueIds.last();
        }
    }
    
    // Add the terminal
    if (addTerminal(term)) {
        return termId;
    }
    
    return -1;
}

bool TerminalManager::removeTerminal(int terminalId) {
    return m_state.removeTerminal(terminalId);
}

void TerminalManager::removeTerminalsFromChain(int chainId) {
    m_state.removeTerminalsFromChain(chainId);
}

void TerminalManager::removeTerminalsByType(TerminalType type) {
    m_state.removeTerminalsByType(type);
}

BiopolymerTerminal* TerminalManager::getTerminal(int terminalId) {
    return m_state.getTerminal(terminalId);
}

const BiopolymerTerminal* TerminalManager::getTerminal(int terminalId) const {
    return m_state.getTerminal(terminalId);
}

QList<BiopolymerTerminal> TerminalManager::allTerminals() const {
    return m_state.terminals();
}

QList<BiopolymerTerminal> TerminalManager::terminalsForChain(int chainId) const {
    return m_state.terminalsForChain(chainId);
}

QList<BiopolymerTerminal> TerminalManager::terminalsForPosition(TerminalPosition position) const {
    QList<BiopolymerTerminal> result;
    for (const BiopolymerTerminal& term : allTerminals()) {
        if (term.position == position) {
            result.append(term);
        }
    }
    return result;
}

QList<BiopolymerTerminal> TerminalManager::terminalsByType(TerminalType type) const {
    return m_state.terminalsByType(type);
}

BiopolymerTerminal* TerminalManager::nTerminal(int chainId) {
    for (BiopolymerTerminal& term : m_state.terminals()) {
        if (term.chainId == chainId && term.isNTerminal()) {
            return &term;
        }
    }
    return nullptr;
}

const BiopolymerTerminal* TerminalManager::nTerminal(int chainId) const {
    for (const BiopolymerTerminal& term : m_state.terminals()) {
        if (term.chainId == chainId && term.isNTerminal()) {
            return &term;
        }
    }
    return nullptr;
}

BiopolymerTerminal* TerminalManager::cTerminal(int chainId) {
    for (BiopolymerTerminal& term : m_state.terminals()) {
        if (term.chainId == chainId && term.isCTerminal()) {
            return &term;
        }
    }
    return nullptr;
}

const BiopolymerTerminal* TerminalManager::cTerminal(int chainId) const {
    for (const BiopolymerTerminal& term : m_state.terminals()) {
        if (term.chainId == chainId && term.isCTerminal()) {
            return &term;
        }
    }
    return nullptr;
}

bool TerminalManager::hasNTerminal(int chainId) const {
    return nTerminal(chainId) != nullptr;
}

bool TerminalManager::hasCTerminal(int chainId) const {
    return cTerminal(chainId) != nullptr;
}

bool TerminalManager::chainHasTerminalType(int chainId, TerminalType type) const {
    for (const BiopolymerTerminal& term : terminalsForChain(chainId)) {
        if (term.type == type) {
            return true;
        }
    }
    return false;
}

int TerminalManager::terminalCount() const {
    return allTerminals().size();
}

int TerminalManager::chainsWithTerminalsCount() const {
    QSet<int> chainIds;
    for (const BiopolymerTerminal& term : allTerminals()) {
        chainIds.insert(term.chainId);
    }
    return chainIds.size();
}

bool TerminalManager::hasUnsupportedTerminals() const {
    for (const BiopolymerTerminal& term : allTerminals()) {
        if (!term.isSupported) {
            return true;
        }
    }
    return false;
}

QList<BiopolymerTerminal> TerminalManager::unsupportedTerminals() const {
    QList<BiopolymerTerminal> result;
    for (const BiopolymerTerminal& term : allTerminals()) {
        if (!term.isSupported) {
            result.append(term);
        }
    }
    return result;
}

QList<QString> TerminalManager::validateAll() const {
    TerminalValidator validator(m_state);
    return validator.validateAllTerminals();
}

bool TerminalManager::setNTerminal(int chainId, TerminalType type) {
    // Remove existing N-terminal
    BiopolymerTerminal* existing = nTerminal(chainId);
    if (existing) {
        removeTerminal(existing->terminalId);
    }
    
    // Add new N-terminal
    int termId = addTerminal(chainId, type, TerminalPosition::First);
    return termId >= 0;
}

bool TerminalManager::setCTerminal(int chainId, TerminalType type) {
    // Remove existing C-terminal
    BiopolymerTerminal* existing = cTerminal(chainId);
    if (existing) {
        removeTerminal(existing->terminalId);
    }
    
    // Add new C-terminal
    int termId = addTerminal(chainId, type, TerminalPosition::Last);
    return termId >= 0;
}

TerminalType TerminalManager::defaultNTerminal(PolymerType polymerType) const {
    switch (polymerType) {
        case PolymerType::PEPTIDE:
            return TerminalType::Amine;
        case PolymerType::DNA:
        case PolymerType::RNA:
            return TerminalType::Phosphate_5;
        default:
            return TerminalType::None;
    }
}

TerminalType TerminalManager::defaultCTerminal(PolymerType polymerType) const {
    switch (polymerType) {
        case PolymerType::PEPTIDE:
            return TerminalType::Carboxyl;
        case PolymerType::DNA:
        case PolymerType::RNA:
            return TerminalType::Hydroxyl_3;
        default:
            return TerminalType::None;
    }
}

bool TerminalManager::allChainsHaveTerminals() const {
    for (int chainId : m_state.chainIds()) {
        if (!hasNTerminal(chainId) || !hasCTerminal(chainId)) {
            return false;
        }
    }
    return true;
}

QList<int> TerminalManager::chainsMissingNTerminals() const {
    QList<int> missing;
    for (int chainId : m_state.chainIds()) {
        if (!hasNTerminal(chainId)) {
            missing.append(chainId);
        }
    }
    return missing;
}

QList<int> TerminalManager::chainsMissingCTerminals() const {
    QList<int> missing;
    for (int chainId : m_state.chainIds()) {
        if (!hasCTerminal(chainId)) {
            missing.append(chainId);
        }
    }
    return missing;
}

int TerminalManager::nextTerminalId() const {
    int maxId = 0;
    for (const BiopolymerTerminal& term : allTerminals()) {
        maxId = std::max(maxId, term.terminalId);
    }
    return maxId + 1;
}
