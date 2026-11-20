// src/app/molecule/BiopolymerModification.cpp
#include "BiopolymerModification.h"
#include "BiopolymerState.h"
#include "EditableMolecule.h"

// ========== Modification Type Helpers ==========

QString modificationTypeToString(ModificationType type) {
    switch (type) {
        case ModificationType::None: return "None";
        case ModificationType::Phosphorylation: return "Phosphorylation";
        case ModificationType::Glycosylation: return "Glycosylation";
        case ModificationType::Acetylation: return "Acetylation";
        case ModificationType::Methylation: return "Methylation";
        case ModificationType::Sulfation: return "Sulfation";
        case ModificationType::Hydroxylation: return "Hydroxylation";
        case ModificationType::Carboxylation: return "Carboxylation";
        case ModificationType::Nitrosylation: return "Nitrosylation";
        case ModificationType::Ubiquitination: return "Ubiquitination";
        case ModificationType::Sumoylation: return "Sumoylation";
        case ModificationType::Prenylation: return "Prenylation";
        case ModificationType::Palmitoylation: return "Palmitoylation";
        case ModificationType::Myristoylation: return "Myristoylation";
        case ModificationType::Custom: return "Custom";
        case ModificationType::TerminalGroup: return "TerminalGroup";
        case ModificationType::Crosslink: return "Crosslink";
        case ModificationType::Cyclization: return "Cyclization";
        case ModificationType::BranchPoint: return "BranchPoint";
        case ModificationType::Unknown: return "Unknown";
    }
    return "Unknown";
}

ModificationType modificationTypeFromString(const QString& typeString) {
    static const QMap<QString, ModificationType> typeMap = {
        {"none", ModificationType::None},
        {"phosphorylation", ModificationType::Phosphorylation},
        {"p", ModificationType::Phosphorylation},
        {"glycosylation", ModificationType::Glycosylation},
        {"acetylation", ModificationType::Acetylation},
        {"ac", ModificationType::Acetylation},
        {"methylation", ModificationType::Methylation},
        {"sulfation", ModificationType::Sulfation},
        {"hydroxylation", ModificationType::Hydroxylation},
        {"carboxylation", ModificationType::Carboxylation},
        {"nitrosylation", ModificationType::Nitrosylation},
        {"ubiquitination", ModificationType::Ubiquitination},
        {"sumoylation", ModificationType::Sumoylation},
        {"prenylation", ModificationType::Prenylation},
        {"palmitoylation", ModificationType::Palmitoylation},
        {"myristoylation", ModificationType::Myristoylation},
        {"custom", ModificationType::Custom},
        {"terminalgroup", ModificationType::TerminalGroup},
        {"terminal", ModificationType::TerminalGroup},
        {"crosslink", ModificationType::Crosslink},
        {"cyclization", ModificationType::Cyclization},
        {"cycle", ModificationType::Cyclization},
        {"branchpoint", ModificationType::BranchPoint},
        {"branch", ModificationType::BranchPoint}
    };
    
    QString lower = typeString.toLower();
    auto it = typeMap.find(lower);
    if (it != typeMap.end()) {
        return it.value();
    }
    return ModificationType::Unknown;
}

QList<ModificationType> allModificationTypes() {
    return {
        ModificationType::None,
        ModificationType::Phosphorylation,
        ModificationType::Glycosylation,
        ModificationType::Acetylation,
        ModificationType::Methylation,
        ModificationType::Sulfation,
        ModificationType::Hydroxylation,
        ModificationType::Carboxylation,
        ModificationType::Nitrosylation,
        ModificationType::Ubiquitination,
        ModificationType::Sumoylation,
        ModificationType::Prenylation,
        ModificationType::Palmitoylation,
        ModificationType::Myristoylation,
        ModificationType::Custom,
        ModificationType::TerminalGroup,
        ModificationType::Crosslink,
        ModificationType::Cyclization,
        ModificationType::BranchPoint,
        ModificationType::Unknown
    };
}

bool isTerminalModification(ModificationType type) {
    return type == ModificationType::TerminalGroup;
}

bool isPostTranslationalModification(ModificationType type) {
    return type == ModificationType::Phosphorylation ||
           type == ModificationType::Glycosylation ||
           type == ModificationType::Acetylation ||
           type == ModificationType::Methylation ||
           type == ModificationType::Sulfation ||
           type == ModificationType::Hydroxylation ||
           type == ModificationType::Carboxylation ||
           type == ModificationType::Nitrosylation ||
           type == ModificationType::Ubiquitination ||
           type == ModificationType::Sumoylation ||
           type == ModificationType::Prenylation ||
           type == ModificationType::Palmitoylation ||
           type == ModificationType::Myristoylation;
}

// ========== Modification Position Helpers ==========

QString modificationPositionToString(ModificationPosition position) {
    switch (position) {
        case ModificationPosition::AlphaCarbon: return "AlphaCarbon";
        case ModificationPosition::NTerminal: return "NTerminal";
        case ModificationPosition::CTerminal: return "CTerminal";
        case ModificationPosition::SideChain: return "SideChain";
        case ModificationPosition::Backbone: return "Backbone";
        case ModificationPosition::Any: return "Any";
        case ModificationPosition::Unknown: return "Unknown";
    }
    return "Unknown";
}

ModificationPosition modificationPositionFromString(const QString& positionString) {
    static const QMap<QString, ModificationPosition> positionMap = {
        {"alphacarbon", ModificationPosition::AlphaCarbon},
        {"ca", ModificationPosition::AlphaCarbon},
        {"nterminal", ModificationPosition::NTerminal},
        {"n-term", ModificationPosition::NTerminal},
        {"n", ModificationPosition::NTerminal},
        {"cterminal", ModificationPosition::CTerminal},
        {"c-terminal", ModificationPosition::CTerminal},
        {"c", ModificationPosition::CTerminal},
        {"sidechain", ModificationPosition::SideChain},
        {"sc", ModificationPosition::SideChain},
        {"backbone", ModificationPosition::Backbone},
        {"bb", ModificationPosition::Backbone},
        {"any", ModificationPosition::Any},
        {"", ModificationPosition::Any}
    };
    
    QString lower = positionString.toLower();
    auto it = positionMap.find(lower);
    if (it != positionMap.end()) {
        return it.value();
    }
    return ModificationPosition::Unknown;
}

// ========== ModificationTarget Implementation ==========

bool ModificationTarget::isValid() const {
    return atomId >= 0 || role != AtomRole::None;
}

bool ModificationTarget::matchesRole(AtomRole targetRole) const {
    return role == targetRole;
}

// ========== BiopolymerModification Implementation ==========

BiopolymerModification::BiopolymerModification() {
    // Default constructor
}

BiopolymerModification::BiopolymerModification(int modId, int resId, ModificationType modType,
                                               const QString& modName, const QString& shortName)
    : modificationId(modId), residueId(resId), type(modType), name(modName), shortName(shortName) {
    // Initialize based on type
    if (modType == ModificationType::TerminalGroup) {
        position = ModificationPosition::Any;
    } else if (modType == ModificationType::Crosslink || modType == ModificationType::Cyclization) {
        position = ModificationPosition::Any;
    } else if (modType == ModificationType::BranchPoint) {
        position = ModificationPosition::AlphaCarbon;
        isBranchPoint = true;
    } else {
        position = ModificationPosition::SideChain;
    }
}

bool BiopolymerModification::isValid() const {
    return modificationId >= 0 && residueId >= 0 && type != ModificationType::None;
}

bool BiopolymerModification::isTerminal() const {
    return type == ModificationType::TerminalGroup || isTerminalModification(type);
}

bool BiopolymerModification::isPTM() const {
    return isPostTranslationalModification(type);
}

bool BiopolymerModification::isCrosslink() const {
    return type == ModificationType::Crosslink || 
           (type == ModificationType::TerminalGroup && name.contains("disulfide", Qt::CaseInsensitive));
}

bool BiopolymerModification::isCyclization() const {
    return type == ModificationType::Cyclization || isCyclic;
}

bool BiopolymerModification::isBranch() const {
    return type == ModificationType::BranchPoint || isBranchPoint;
}

QString BiopolymerModification::fullName() const {
    if (name.isEmpty()) {
        return modificationTypeToString(type);
    }
    return name;
}

QString BiopolymerModification::displayName() const {
    if (!shortName.isEmpty()) {
        return shortName;
    }
    return fullName();
}

// ========== Static Factory Methods ==========

BiopolymerModification BiopolymerModification::createPhosphorylation(int modId, int residueId,
                                                                    const QString& atomName) {
    BiopolymerModification mod(modId, residueId, ModificationType::Phosphorylation,
                               "Phosphorylation", "p");
    mod.chemicalFormula = "PO3H2";
    mod.molecularWeight = 94.97;  // PO3H2 molecular weight
    mod.description = "Phosphate group addition";
    mod.target.atomName = atomName;
    return mod;
}

BiopolymerModification BiopolymerModification::createGlycosylation(int modId, int residueId,
                                                                   const QString& atomName) {
    BiopolymerModification mod(modId, residueId, ModificationType::Glycosylation,
                               "N-Glycosylation", "N-Gly");
    mod.chemicalFormula = "C8H13N1O5";  // Typical N-acetylglucosamine
    mod.molecularWeight = 203.2;  // N-acetylglucosamine MW
    mod.description = "N-linked glycosylation";
    mod.target.atomName = atomName;
    return mod;
}

BiopolymerModification BiopolymerModification::createAcetylation(int modId, int residueId,
                                                                  ModificationPosition position) {
    BiopolymerModification mod(modId, residueId, ModificationType::Acetylation,
                               "Acetylation", "Ac");
    mod.chemicalFormula = "C2H3O";
    mod.molecularWeight = 43.05;
    mod.description = "Acetyl group addition";
    mod.position = position;
    mod.isSupported = true;
    return mod;
}

BiopolymerModification BiopolymerModification::createAmideTerminal(int modId, int residueId) {
    BiopolymerModification mod(modId, residueId, ModificationType::TerminalGroup,
                               "C-Terminal Amide", "NH2");
    mod.chemicalFormula = "NH2";
    mod.molecularWeight = 16.0;  // Approximate
    mod.description = "C-terminal amide group";
    mod.position = ModificationPosition::CTerminal;
    mod.isSupported = true;
    return mod;
}

BiopolymerModification BiopolymerModification::createDisulfide(int modId, int residueId1, int residueId2) {
    BiopolymerModification mod(modId, residueId1, ModificationType::Crosslink,
                               "Disulfide Bond", "SS");
    mod.chemicalFormula = "S2";
    mod.molecularWeight = 32.07 * 2;  // Two sulfur atoms
    mod.description = "Disulfide bond between cysteine residues";
    mod.targetResidueId = residueId2;
    mod.isSupported = true;
    return mod;
}

BiopolymerModification BiopolymerModification::createCyclization(int modId, int residueId1, int residueId2) {
    BiopolymerModification mod(modId, residueId1, ModificationType::Cyclization,
                               "Cyclization", "Cycle");
    mod.description = "Head-to-tail cyclization";
    mod.cyclePartnerResidueId = residueId2;
    mod.isCyclic = true;
    mod.isSupported = true;
    return mod;
}

BiopolymerModification BiopolymerModification::createBranchPoint(int modId, int residueId,
                                                                   const QList<int>& branchResidues) {
    BiopolymerModification mod(modId, residueId, ModificationType::BranchPoint,
                               "Branch Point", "Branch");
    mod.description = "Branch point in biopolymer topology";
    mod.branchResidues = branchResidues;
    mod.isBranchPoint = true;
    mod.position = ModificationPosition::AlphaCarbon;
    return mod;
}

bool BiopolymerModification::operator==(const BiopolymerModification& other) const {
    return modificationId == other.modificationId &&
           residueId == other.residueId &&
           type == other.type &&
           name == other.name &&
           shortName == other.shortName &&
           position == other.position;
}

bool BiopolymerModification::operator!=(const BiopolymerModification& other) const {
    return !(*this == other);
}

// ========== Modification Database Implementation ==========

// Static initialization of common modifications
const ModificationDatabase::ModificationRecord ModificationDatabase::Phosphorylation_Serine = {
    ModificationType::Phosphorylation,
    "Phosphorylation",
    "p",
    "PO3H2",
    94.97,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"S", "ser", "Ser", "SER"},
    true,
    "Phosphorylation of serine residue"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::Phosphorylation_Threonine = {
    ModificationType::Phosphorylation,
    "Phosphorylation",
    "p",
    "PO3H2",
    94.97,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"T", "thr", "Thr", "THR"},
    true,
    "Phosphorylation of threonine residue"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::Phosphorylation_Tyrosine = {
    ModificationType::Phosphorylation,
    "Phosphorylation",
    "p",
    "PO3H2",
    94.97,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"Y", "tyr", "Tyr", "TYR"},
    true,
    "Phosphorylation of tyrosine residue"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::N_Glycosylation = {
    ModificationType::Glycosylation,
    "N-Glycosylation",
    "N-Gly",
    "C8H13N1O5",
    203.2,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"N", "asn", "Asn", "ASN"},
    true,
    "N-linked glycosylation at asparagine"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::O_Glycosylation = {
    ModificationType::Glycosylation,
    "O-Glycosylation",
    "O-Gly",
    "C6H10O5",
    162.14,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"S", "ser", "Ser", "SER", "T", "thr", "Thr", "THR"},
    true,
    "O-linked glycosylation at serine/threonine"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::N_Terminal_Acetylation = {
    ModificationType::Acetylation,
    "N-Terminal Acetylation",
    "Ac",
    "C2H3O",
    43.05,
    ModificationPosition::NTerminal,
    {PolymerType::Peptide, PolymerType::DNA, PolymerType::RNA},
    {},
    true,
    "Acetylation of N-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::N_Terminal_Methylation = {
    ModificationType::Methylation,
    "N-Terminal Methylation",
    "Me",
    "CH3",
    15.04,
    ModificationPosition::NTerminal,
    {PolymerType::Peptide},
    {},
    true,
    "Methylation of N-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::Lysine_Methylation = {
    ModificationType::Methylation,
    "Lysine Methylation",
    "Me",
    "CH3",
    15.04,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"K", "lys", "Lys", "LYS"},
    true,
    "Methylation of lysine side chain"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::N_Terminal_Acetyl = {
    ModificationType::TerminalGroup,
    "N-Terminal Acetyl",
    "Ac",
    "C2H3O",
    43.05,
    ModificationPosition::NTerminal,
    {PolymerType::Peptide},
    {},
    true,
    "Acetyl group at N-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::N_Terminal_Amide = {
    ModificationType::TerminalGroup,
    "N-Terminal Amide",
    "NH2",
    "NH2",
    16.0,
    ModificationPosition::NTerminal,
    {PolymerType::Peptide},
    {},
    true,
    "Amide group at N-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::C_Terminal_Amide = {
    ModificationType::TerminalGroup,
    "C-Terminal Amide",
    "NH2",
    "NH2",
    16.0,
    ModificationPosition::CTerminal,
    {PolymerType::Peptide},
    {},
    true,
    "Amide group at C-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::C_Terminal_CarboxylicAcid = {
    ModificationType::TerminalGroup,
    "C-Terminal Carboxylic Acid",
    "COOH",
    "CO2H",
    44.01,
    ModificationPosition::CTerminal,
    {PolymerType::Peptide},
    {},
    true,
    "Carboxylic acid at C-terminus"
};

const ModificationDatabase::ModificationRecord ModificationDatabase::Disulfide_Bond = {
    ModificationType::Crosslink,
    "Disulfide Bond",
    "SS",
    "S2",
    64.14,
    ModificationPosition::SideChain,
    {PolymerType::Peptide},
    {"C", "cys", "Cys", "CYS"},
    true,
    "Disulfide bond between cysteine residues"
};

ModificationDatabase& ModificationDatabase::instance() {
    static ModificationDatabase db;
    return db;
}

ModificationDatabase::ModificationDatabase() {
    initialize();
}

ModificationDatabase::~ModificationDatabase() {
    // Cleanup
}

void ModificationDatabase::initialize() {
    // Register common modifications
    registerCustomModification(Phosphorylation_Serine);
    registerCustomModification(Phosphorylation_Threonine);
    registerCustomModification(Phosphorylation_Tyrosine);
    registerCustomModification(N_Glycosylation);
    registerCustomModification(O_Glycosylation);
    registerCustomModification(N_Terminal_Acetylation);
    registerCustomModification(N_Terminal_Methylation);
    registerCustomModification(Lysine_Methylation);
    registerCustomModification(N_Terminal_Acetyl);
    registerCustomModification(N_Terminal_Amide);
    registerCustomModification(C_Terminal_Amide);
    registerCustomModification(C_Terminal_CarboxylicAcid);
    registerCustomModification(Disulfide_Bond);
}

ModificationDatabase::ModificationRecord ModificationDatabase::getModification(const QString& name) const {
    auto it = m_modificationsByName.find(name);
    if (it != m_modificationsByName.end()) {
        return it.value();
    }
    return ModificationRecord();
}

ModificationDatabase::ModificationRecord ModificationDatabase::getModification(ModificationType type,
                                                                                const QString& residueSymbol) const {
    auto it = m_modificationsByType.find(type);
    if (it != m_modificationsByType.end()) {
        for (const ModificationRecord& record : it.value()) {
            if (record.applicableResidues.isEmpty() || 
                record.applicableResidues.contains(residueSymbol)) {
                return record;
            }
        }
    }
    return ModificationRecord();
}

QList<ModificationDatabase::ModificationRecord> ModificationDatabase::allModifications() const {
    QList<ModificationRecord> all;
    for (const auto& entry : m_modificationsByName) {
        all.append(entry);
    }
    return all;
}

QList<ModificationDatabase::ModificationRecord> ModificationDatabase::modificationsForResidue(
    const QString& residueSymbol, PolymerType polymerType) const {
    QList<ModificationRecord> result;
    for (const auto& entry : m_modificationsByName) {
        const ModificationRecord& record = entry;
        if (record.applicableTo.contains(polymerType) || record.applicableTo.isEmpty()) {
            if (record.applicableResidues.isEmpty() || 
                record.applicableResidues.contains(residueSymbol)) {
                result.append(record);
            }
        }
    }
    return result;
}

QList<ModificationDatabase::ModificationRecord> ModificationDatabase::modificationsForPolymerType(
    PolymerType polymerType) const {
    QList<ModificationRecord> result;
    for (const auto& entry : m_modificationsByName) {
        const ModificationRecord& record = entry;
        if (record.applicableTo.contains(polymerType) || record.applicableTo.isEmpty()) {
            result.append(record);
        }
    }
    return result;
}

bool ModificationDatabase::isKnownModification(const QString& name) const {
    return m_modificationsByName.contains(name);
}

void ModificationDatabase::registerCustomModification(const ModificationRecord& record) {
    if (!record.name.isEmpty()) {
        m_modificationsByName[record.name] = record;
        m_modificationsByType[record.type].append(record);
    }
}

// ========== ModificationValidator Implementation ==========

ModificationValidator::ModificationValidator(const BiopolymerState& state)
    : m_state(state) {
    // Initialize validator
}

bool ModificationValidator::validateModification(const BiopolymerModification& modification) const {
    // Check if modification is valid
    if (!modification.isValid()) {
        return false;
    }
    
    // Check residue exists
    const BiopolymerResidue* residue = m_state.residue(modification.residueId);
    if (!residue) {
        return false;
    }
    
    // Check compatibility
    if (!isCompatible(modification.residueId, modification)) {
        return false;
    }
    
    // Check for conflicts
    if (wouldConflict(modification.residueId, modification)) {
        return false;
    }
    
    return true;
}

QList<QString> ModificationValidator::validateAllModifications() const {
    QList<QString> issues;
    
    for (const BiopolymerModification& mod : m_state.modifications()) {
        if (!validateModification(mod)) {
            issues.append(QString("Invalid modification %1 on residue %2")
                .arg(mod.displayName()).arg(mod.residueId));
        }
    }
    
    return issues;
}

bool ModificationValidator::isCompatible(int residueId, const BiopolymerModification& modification) const {
    // Check residue type
    if (!checkResidueType(residueId, modification.type)) {
        return false;
    }
    
    // Check position
    if (!checkPosition(residueId, modification.position)) {
        return false;
    }
    
    // Check atom availability
    if (!checkAtomAvailability(residueId, modification.target)) {
        return false;
    }
    
    return true;
}

bool ModificationValidator::canHaveModificationType(int residueId, ModificationType type) const {
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return false;
    }
    
    // Get applicable modifications
    ModificationDatabase& db = ModificationDatabase::instance();
    const BiopolymerChain* chain = m_state.chain(residue->chainId);
    PolymerType polyType = chain ? chain->polymerType : PolymerType::Unknown;
    QList<ModificationDatabase::ModificationRecord> records = 
        db.modificationsForResidue(residue->identity.symbol, polyType);
    
    for (const auto& record : records) {
        if (record.type == type) {
            return true;
        }
    }
    
    return false;
}

bool ModificationValidator::wouldConflict(int residueId, const BiopolymerModification& modification) const {
    return !getConflictingModifications(residueId, modification).isEmpty();
}

QList<int> ModificationValidator::getConflictingModifications(int residueId,
                                                              const BiopolymerModification& modification) const {
    QList<int> conflicts;
    
    // Get existing modifications for this residue
    for (const BiopolymerModification& existing : m_state.modifications()) {
        if (existing.residueId == residueId) {
            // Same position conflicts
            if (existing.position == modification.position && 
                existing.modificationId != modification.modificationId) {
                conflicts.append(existing.modificationId);
            }
            
            // Same target atom conflicts
            if (existing.target.isValid() && modification.target.isValid() &&
                existing.target.atomId == modification.target.atomId &&
                existing.modificationId != modification.modificationId) {
                conflicts.append(existing.modificationId);
            }
        }
    }
    
    return conflicts;
}

QList<BiopolymerModification> ModificationValidator::getSuggestedModifications(int residueId) const {
    QList<BiopolymerModification> suggestions;
    
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return suggestions;
    }
    
    // Get applicable modifications from database
    ModificationDatabase& db = ModificationDatabase::instance();
    const BiopolymerChain* chain = m_state.chain(residue->chainId);
    PolymerType polyType = chain ? chain->polymerType : PolymerType::Unknown;
    QList<ModificationDatabase::ModificationRecord> records = 
        db.modificationsForResidue(residue->identity.symbol, polyType);
    
    // Create suggestions
    int nextId = 1;  // Will be replaced with real ID when added
    for (const auto& record : records) {
        if (record.isSupported) {
            BiopolymerModification mod(nextId++, residueId, record.type,
                                      record.name, record.shortName);
            mod.chemicalFormula = record.chemicalFormula;
            mod.molecularWeight = record.molecularWeight;
            mod.description = record.description;
            mod.position = record.defaultPosition;
            suggestions.append(mod);
        }
    }
    
    return suggestions;
}

bool ModificationValidator::checkResidueType(int residueId, ModificationType type) const {
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return false;
    }
    
    ModificationDatabase& db = ModificationDatabase::instance();
    const BiopolymerChain* chain = m_state.chain(residue->chainId);
    PolymerType polyType = chain ? chain->polymerType : PolymerType::Unknown;
    QList<ModificationDatabase::ModificationRecord> records = 
        db.modificationsForPolymerType(polyType);
    
    for (const auto& record : records) {
        if (record.type == type) {
            return true;
        }
    }
    
    return false;
}

bool ModificationValidator::checkPosition(int residueId, ModificationPosition position) const {
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return false;
    }
    
    // Check if position is valid for the residue type
    if (position == ModificationPosition::NTerminal || position == ModificationPosition::CTerminal) {
        // Check if residue is at terminal
        const BiopolymerChain* chain = m_state.chain(residue->chainId);
        if (chain) {
            if (position == ModificationPosition::NTerminal) {
                return chain->residueIds.first() == residueId;
            } else {
                return chain->residueIds.last() == residueId;
            }
        }
        return false;
    }
    
    return true;
}

bool ModificationValidator::checkAtomAvailability(int residueId, const ModificationTarget& target) const {
    if (!target.isValid()) {
        return true;  // No specific target, so it's available
    }
    
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return false;
    }
    
    // Check if the atom exists in the residue
    if (target.atomId >= 0) {
        // For now, we assume the atom exists if we have a valid ID
        // In a full implementation, we'd check against the actual atom list
        return true;
    }
    
    // Check if the role exists
    // Note: Atom role checking requires EditableMolecule access which is not available
    // For now, we assume the role is valid if it's not None
    if (target.role != AtomRole::None) {
        // Without access to EditableMolecule, we can't verify atom roles
        // Just return true to allow the modification
        return true;
    }
    
    return true;
}

// ========== ModificationManager Implementation ==========

ModificationManager::ModificationManager(BiopolymerState& state)
    : m_state(state) {
    // Initialize manager
}

bool ModificationManager::addModification(const BiopolymerModification& modification) {
    // Validate the modification
    ModificationValidator validator(m_state);
    if (!validator.validateModification(modification)) {
        return false;
    }
    
    // Add to state
    return m_state.addModification(modification);
}

int ModificationManager::addModification(int residueId, ModificationType type, const QString& name) {
    // Create a new modification
    int modId = nextModificationId();
    BiopolymerModification mod(modId, residueId, type, name, "");
    
    // Get details from database
    ModificationDatabase& db = ModificationDatabase::instance();
    ModificationDatabase::ModificationRecord record = db.getModification(name);
    if (!name.isEmpty() && !record.name.isEmpty()) {
        mod.name = record.name;
        mod.shortName = record.shortName;
        mod.chemicalFormula = record.chemicalFormula;
        mod.molecularWeight = record.molecularWeight;
        mod.description = record.description;
        mod.position = record.defaultPosition;
        mod.isSupported = record.isSupported;
    }
    
    // Add the modification
    if (addModification(mod)) {
        return modId;
    }
    
    return -1;
}

bool ModificationManager::removeModification(int modificationId) {
    return m_state.removeModification(modificationId);
}

void ModificationManager::removeModificationsFromResidue(int residueId) {
    m_state.removeModificationsFromResidue(residueId);
}

void ModificationManager::removeModificationsByType(ModificationType type) {
    m_state.removeModificationsByType(type);
}

BiopolymerModification* ModificationManager::getModification(int modificationId) {
    return m_state.getModification(modificationId);
}

const BiopolymerModification* ModificationManager::getModification(int modificationId) const {
    return m_state.getModification(modificationId);
}

QList<BiopolymerModification> ModificationManager::allModifications() const {
    return m_state.modifications();
}

QList<BiopolymerModification> ModificationManager::modificationsForResidue(int residueId) const {
    return m_state.modificationsForResidue(residueId);
}

QList<BiopolymerModification> ModificationManager::modificationsByType(ModificationType type) const {
    return m_state.modificationsByType(type);
}

bool ModificationManager::residueHasModifications(int residueId) const {
    return !modificationsForResidue(residueId).isEmpty();
}

bool ModificationManager::residueHasModificationType(int residueId, ModificationType type) const {
    for (const BiopolymerModification& mod : modificationsForResidue(residueId)) {
        if (mod.type == type) {
            return true;
        }
    }
    return false;
}

int ModificationManager::modificationCount() const {
    return allModifications().size();
}

int ModificationManager::modifiedResidueCount() const {
    QSet<int> residueIds;
    for (const BiopolymerModification& mod : allModifications()) {
        residueIds.insert(mod.residueId);
    }
    return residueIds.size();
}

QList<BiopolymerModification> ModificationManager::findModificationsAtPosition(
    ModificationPosition position) const {
    QList<BiopolymerModification> result;
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.position == position) {
            result.append(mod);
        }
    }
    return result;
}

QList<BiopolymerModification> ModificationManager::findTerminalModifications() const {
    return findModificationsAtPosition(ModificationPosition::NTerminal) +
           findModificationsAtPosition(ModificationPosition::CTerminal);
}

QList<BiopolymerModification> ModificationManager::findCrosslinks() const {
    QList<BiopolymerModification> result;
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isCrosslink()) {
            result.append(mod);
        }
    }
    return result;
}

QList<BiopolymerModification> ModificationManager::findCyclizations() const {
    QList<BiopolymerModification> result;
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isCyclization()) {
            result.append(mod);
        }
    }
    return result;
}

QList<BiopolymerModification> ModificationManager::findBranchPoints() const {
    QList<BiopolymerModification> result;
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isBranch()) {
            result.append(mod);
        }
    }
    return result;
}

bool ModificationManager::hasUnsupportedModifications() const {
    for (const BiopolymerModification& mod : allModifications()) {
        if (!mod.isSupported) {
            return true;
        }
    }
    return false;
}

QList<BiopolymerModification> ModificationManager::unsupportedModifications() const {
    QList<BiopolymerModification> result;
    for (const BiopolymerModification& mod : allModifications()) {
        if (!mod.isSupported) {
            result.append(mod);
        }
    }
    return result;
}

QList<QString> ModificationManager::validateAll() const {
    ModificationValidator validator(m_state);
    return validator.validateAllModifications();
}

// ========== Crosslink Management ==========

bool ModificationManager::addDisulfideBond(int residueId1, int residueId2) {
    const BiopolymerResidue* residue1 = m_state.residue(residueId1);
    const BiopolymerResidue* residue2 = m_state.residue(residueId2);
    
    if (!residue1 || !residue2) {
        return false;
    }
    
    // Check if both residues are cysteine
    if (residue1->identity.symbol != "C" && residue1->identity.symbol.toUpper() != "CYS") {
        return false;
    }
    if (residue2->identity.symbol != "C" && residue2->identity.symbol.toUpper() != "CYS") {
        return false;
    }
    
    // Check if already connected
    if (hasDisulfideBond(residueId1, residueId2)) {
        return false;
    }
    
    // Add the modification
    int modId = nextModificationId();
    BiopolymerModification mod = BiopolymerModification::createDisulfide(modId, residueId1, residueId2);
    
    return addModification(mod);
}

bool ModificationManager::removeDisulfideBond(int residueId1, int residueId2) {
    // Find and remove the disulfide modification
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.type == ModificationType::Crosslink && 
            ((mod.residueId == residueId1 && mod.targetResidueId == residueId2) ||
             (mod.residueId == residueId2 && mod.targetResidueId == residueId1))) {
            return removeModification(mod.modificationId);
        }
    }
    return false;
}

bool ModificationManager::hasDisulfideBond(int residueId1, int residueId2) const {
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.type == ModificationType::Crosslink && 
            ((mod.residueId == residueId1 && mod.targetResidueId == residueId2) ||
             (mod.residueId == residueId2 && mod.targetResidueId == residueId1))) {
            return true;
        }
    }
    return false;
}

QList<std::pair<int, int>> ModificationManager::allDisulfideBonds() const {
    QList<std::pair<int, int>> bonds;
    QSet<std::pair<int, int>> seen;
    
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.type == ModificationType::Crosslink && mod.targetResidueId >= 0) {
            int id1 = std::min(mod.residueId, mod.targetResidueId);
            int id2 = std::max(mod.residueId, mod.targetResidueId);
            std::pair<int, int> bond = {id1, id2};
            
            if (!seen.contains(bond)) {
                seen.insert(bond);
                bonds.append(bond);
            }
        }
    }
    
    return bonds;
}

// ========== Cyclization Management ==========

bool ModificationManager::addCyclization(int residueId1, int residueId2) {
    const BiopolymerResidue* residue1 = m_state.residue(residueId1);
    const BiopolymerResidue* residue2 = m_state.residue(residueId2);
    
    if (!residue1 || !residue2) {
        return false;
    }
    
    // Check if they're in the same chain
    if (residue1->chainId != residue2->chainId) {
        return false;
    }
    
    // Check if already cyclic
    if (isCyclic(residue1->chainId)) {
        return false;
    }
    
    // Add the cyclization
    int modId = nextModificationId();
    BiopolymerModification mod = BiopolymerModification::createCyclization(modId, residueId1, residueId2);
    
    return addModification(mod);
}

bool ModificationManager::isCyclic(int chainId) const {
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) {
        return false;
    }
    
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.type == ModificationType::Cyclization) {
            const BiopolymerResidue* residue = m_state.residue(mod.residueId);
            if (residue && residue->chainId == chainId) {
                return true;
            }
        }
    }
    
    return false;
}

QList<std::pair<int, int>> ModificationManager::getCyclizationPoints(int chainId) const {
    QList<std::pair<int, int>> points;
    QSet<std::pair<int, int>> seen;
    
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.type == ModificationType::Cyclization && mod.cyclePartnerResidueId >= 0) {
            const BiopolymerResidue* residue = m_state.residue(mod.residueId);
            if (residue && residue->chainId == chainId) {
                int id1 = std::min(mod.residueId, mod.cyclePartnerResidueId);
                int id2 = std::max(mod.residueId, mod.cyclePartnerResidueId);
                std::pair<int, int> point = {id1, id2};
                
                if (!seen.contains(point)) {
                    seen.insert(point);
                    points.append(point);
                }
            }
        }
    }
    
    return points;
}

// ========== Branch Management ==========

bool ModificationManager::addBranchPoint(int residueId, const QList<int>& branchResidues) {
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) {
        return false;
    }
    
    // Check if residue can be a branch point
    // (In a full implementation, we'd check the residue type)
    
    // Add the branch point
    int modId = nextModificationId();
    BiopolymerModification mod = BiopolymerModification::createBranchPoint(modId, residueId, branchResidues);
    
    return addModification(mod);
}

bool ModificationManager::isBranchPoint(int residueId) const {
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isBranch() && mod.residueId == residueId) {
            return true;
        }
    }
    return false;
}

QList<int> ModificationManager::getBranchResidues(int residueId) const {
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isBranch() && mod.residueId == residueId) {
            return mod.branchResidues;
        }
    }
    return QList<int>();
}

QList<int> ModificationManager::allBranchPoints() const {
    QList<int> points;
    for (const BiopolymerModification& mod : allModifications()) {
        if (mod.isBranch() && !points.contains(mod.residueId)) {
            points.append(mod.residueId);
        }
    }
    return points;
}

int ModificationManager::nextModificationId() const {
    int maxId = 0;
    for (const BiopolymerModification& mod : allModifications()) {
        maxId = std::max(maxId, mod.modificationId);
    }
    return maxId + 1;
}
