// src/app/molecule/BiopolymerModification.h
#ifndef BIOPOLYMERMODIFICATION_H
#define BIOPOLYMERMODIFICATION_H

#include "EditableMolecule.h"  // For AtomRole
#include <QString>
#include <QList>
#include <QMap>
#include <QSet>

// Forward declarations to break circular dependencies
enum class PolymerType;
class BiopolymerState;

// ========== Modification Type ==========

// Types of biopolymer modifications
enum class ModificationType {
    None,           // No modification
    Phosphorylation, // Phosphorylation (e.g., pS, pT, pY)
    Glycosylation,  // Glycosylation (N-linked, O-linked)
    Acetylation,    // Acetylation
    Methylation,    // Methylation
    Sulfation,      // Sulfation
    Hydroxylation,  // Hydroxylation
    Carboxylation,  // Carboxylation
    Nitrosylation,  // Nitrosylation
    Ubiquitination, // Ubiquitination
    Sumoylation,    // Sumoylation
    Prenylation,    // Prenylation
    Palmitoylation,  // Palmitoylation
    Myristoylation,  // Myristoylation
    Custom,         // Custom/modified residue
    TerminalGroup,  // Terminal modification (N-term, C-term)
    Crosslink,      // Crosslink between residues
    Cyclization,    // Cyclic structure
    BranchPoint,    // Branch point in topology
    Unknown         // Unknown modification type
};

// Convert modification type to string
QString modificationTypeToString(ModificationType type);

// Parse modification type from string
ModificationType modificationTypeFromString(const QString& typeString);

// Get all modification types
QList<ModificationType> allModificationTypes();

// Check if modification type is a terminal group
bool isTerminalModification(ModificationType type);

// Check if modification type is a post-translational modification (PTM)
bool isPostTranslationalModification(ModificationType type);

// ========== Modification Position ==========

// Position of a modification within a residue
enum class ModificationPosition {
    AlphaCarbon,    // C-alpha
    NTerminal,      // N-terminus
    CTerminal,      // C-terminus
    SideChain,      // Side chain
    Backbone,       // Backbone atom
    Any,            // Any position
    Unknown         // Unknown position
};

QString modificationPositionToString(ModificationPosition position);
ModificationPosition modificationPositionFromString(const QString& positionString);

// ========== Modification Target ==========

// Target atom/role for a modification
struct ModificationTarget {
    AtomId atomId = -1;           // Target atom ID
    AtomRole role = AtomRole::None; // Target atom role
    QString atomName;             // Target atom name
    
    // Check if target is valid
    bool isValid() const;
    
    // Check if target matches a role
    bool matchesRole(AtomRole targetRole) const;
};

// ========== Biopolymer Modification ==========

// Represents a modification to a biopolymer residue
//
// A modification can be:
// - A chemical modification (e.g., phosphorylation, glycosylation)
// - A terminal group (e.g., acetyl, amide)
// - A structural feature (e.g., branch point, cyclization)
//
// Modifications are associated with specific residues and may affect
// the residue's identity, connectivity, or chemical properties.
struct BiopolymerModification {
    int modificationId = -1;      // Unique modification ID
    int residueId = -1;           // Residue being modified
    ModificationType type = ModificationType::None;
    ModificationPosition position = ModificationPosition::SideChain;
    
    // Modification identity
    QString name;                 // Common name (e.g., "Phosphorylation")
    QString shortName;            // Short name (e.g., "p", "Ac")
    QString chemicalFormula;      // Chemical formula of the modification
    double molecularWeight = 0.0; // Molecular weight added
    
    // Modification target
    ModificationTarget target;     // Target atom/role
    
    // For crosslinks/Connections
    int targetResidueId = -1;     // Target residue for crosslinks
    int targetChainId = -1;       // Target chain for crosslinks
    AtomRole targetRole = AtomRole::None;  // Target role for crosslinks
    
    // For cyclization
    bool isCyclic = false;        // Whether this modification creates a cycle
    int cyclePartnerResidueId = -1; // Partner residue in cycle
    
    // For branched topology
    bool isBranchPoint = false;    // Whether this is a branch point
    QList<int> branchResidues;     // Residues branching from this point
    
    // Validation and support
    bool isSupported = true;       // Whether this modification is supported
    QStringList validationIssues; // Validation issues
    
    // Metadata
    QString description;          // Description of the modification
    QString reference;            // Reference (e.g., UniProt, database)
    QMap<QString, QString> properties; // Additional properties
    
    // ========== Constructors ==========
    
    BiopolymerModification();
    BiopolymerModification(int modId, int resId, ModificationType modType, 
                           const QString& modName, const QString& shortName);
    
    // ========== Accessors ==========
    
    // Check if modification is valid
    bool isValid() const;
    
    // Check if modification is a terminal group
    bool isTerminal() const;
    
    // Check if modification is a PTM
    bool isPTM() const;
    
    // Check if modification is a crosslink
    bool isCrosslink() const;
    
    // Check if modification is a cyclization
    bool isCyclization() const;
    
    // Check if modification is a branch point
    bool isBranch() const;
    
    // Get full name with residue
    QString fullName() const;
    
    // Get display name
    QString displayName() const;
    
    // ========== Static Factory Methods ==========
    
    // Create a phosphorylation modification
    static BiopolymerModification createPhosphorylation(int modId, int residueId, 
                                                        const QString& atomName = QString());
    
    // Create a glycosylation modification
    static BiopolymerModification createGlycosylation(int modId, int residueId, 
                                                       const QString& atomName = QString());
    
    // Create an acetylation modification
    static BiopolymerModification createAcetylation(int modId, int residueId, 
                                                  ModificationPosition position = ModificationPosition::NTerminal);
    
    // Create an amide terminal group
    static BiopolymerModification createAmideTerminal(int modId, int residueId);
    
    // Create a disulfide bond (crosslink)
    static BiopolymerModification createDisulfide(int modId, int residueId1, int residueId2);
    
    // Create a head-to-tail cyclization
    static BiopolymerModification createCyclization(int modId, int residueId1, int residueId2);
    
    // Create a branch point
    static BiopolymerModification createBranchPoint(int modId, int residueId, 
                                                   const QList<int>& branchResidues = QList<int>());
    
    // ========== Comparison ==========
    
    bool operator==(const BiopolymerModification& other) const;
    bool operator!=(const BiopolymerModification& other) const;
};

// ========== Modification Database ==========

// Database of known modifications
class ModificationDatabase {
public:
    // Modification record
    struct ModificationRecord {
        ModificationType type = ModificationType::None;
        QString name;
        QString shortName;
        QString chemicalFormula;
        double molecularWeight = 0.0;
        ModificationPosition defaultPosition = ModificationPosition::SideChain;
        QSet<PolymerType> applicableTo;  // Which polymer types this applies to
        QSet<QString> applicableResidues; // Which residue symbols this applies to
        bool isSupported = true;
        QString description;
    };
    
    // Get the singleton instance
    static ModificationDatabase& instance();
    
    // Get modification by name
    ModificationRecord getModification(const QString& name) const;
    
    // Get modification by type and residue
    ModificationRecord getModification(ModificationType type, const QString& residueSymbol) const;
    
    // Get all modifications
    QList<ModificationRecord> allModifications() const;
    
    // Get modifications for a residue
    QList<ModificationRecord> modificationsForResidue(const QString& residueSymbol, 
                                                       PolymerType polymerType) const;
    
    // Get modifications for a polymer type
    QList<ModificationRecord> modificationsForPolymerType(PolymerType polymerType) const;
    
    // Check if a modification is known
    bool isKnownModification(const QString& name) const;
    
    // Register a custom modification
    void registerCustomModification(const ModificationRecord& record);
    
    // ========== Common Modifications ==========
    
    // Phosphorylation
    static const ModificationRecord Phosphorylation_Serine;
    static const ModificationRecord Phosphorylation_Threonine;
    static const ModificationRecord Phosphorylation_Tyrosine;
    
    // Glycosylation
    static const ModificationRecord N_Glycosylation;
    static const ModificationRecord O_Glycosylation;
    
    // Acetylation
    static const ModificationRecord N_Terminal_Acetylation;
    
    // Methylation
    static const ModificationRecord N_Terminal_Methylation;
    static const ModificationRecord Lysine_Methylation;
    
    // Terminal groups
    static const ModificationRecord N_Terminal_Acetyl;
    static const ModificationRecord N_Terminal_Amide;
    static const ModificationRecord C_Terminal_Amide;
    static const ModificationRecord C_Terminal_CarboxylicAcid;
    
    // Crosslinks
    static const ModificationRecord Disulfide_Bond;
    
private:
    ModificationDatabase();
    ~ModificationDatabase();
    
    // Initialize the database
    void initialize();
    
    QMap<QString, ModificationRecord> m_modificationsByName;
    QMap<ModificationType, QList<ModificationRecord>> m_modificationsByType;
};

// ========== Modification Validation ==========

// Validates modifications for a biopolymer state
class ModificationValidator {
public:
    explicit ModificationValidator(const BiopolymerState& state);
    
    // Validate a modification
    bool validateModification(const BiopolymerModification& modification) const;
    
    // Validate all modifications
    QList<QString> validateAllModifications() const;
    
    // Check if a modification is compatible with the residue
    bool isCompatible(int residueId, const BiopolymerModification& modification) const;
    
    // Check if a residue can have a specific modification type
    bool canHaveModificationType(int residueId, ModificationType type) const;
    
    // Check if adding a modification would conflict
    bool wouldConflict(int residueId, const BiopolymerModification& modification) const;
    
    // Get conflicting modifications
    QList<int> getConflictingModifications(int residueId, 
                                           const BiopolymerModification& modification) const;
    
    // Get suggested modifications for a residue
    QList<BiopolymerModification> getSuggestedModifications(int residueId) const;

private:
    const BiopolymerState& m_state;
    
    // Check residue type compatibility
    bool checkResidueType(int residueId, ModificationType type) const;
    
    // Check position compatibility
    bool checkPosition(int residueId, ModificationPosition position) const;
    
    // Check atom availability
    bool checkAtomAvailability(int residueId, const ModificationTarget& target) const;
};

// ========== Modification Manager ==========

// Manages modifications within a biopolymer state
//
// This class provides high-level operations for working with modifications,
// including adding, removing, and querying modifications.
class ModificationManager {
public:
    explicit ModificationManager(BiopolymerState& state);
    
    // Add a modification to a residue
    bool addModification(const BiopolymerModification& modification);
    
    // Add a modification by type
    int addModification(int residueId, ModificationType type, 
                       const QString& name = QString());
    
    // Remove a modification
    bool removeModification(int modificationId);
    
    // Remove all modifications from a residue
    void removeModificationsFromResidue(int residueId);
    
    // Remove all modifications of a type
    void removeModificationsByType(ModificationType type);
    
    // Get modification by ID
    BiopolymerModification* getModification(int modificationId);
    const BiopolymerModification* getModification(int modificationId) const;
    
    // Get all modifications
    QList<BiopolymerModification> allModifications() const;
    
    // Get modifications for a residue
    QList<BiopolymerModification> modificationsForResidue(int residueId) const;
    
    // Get modifications of a type
    QList<BiopolymerModification> modificationsByType(ModificationType type) const;
    
    // Check if a residue has modifications
    bool residueHasModifications(int residueId) const;
    
    // Check if a residue has a specific modification type
    bool residueHasModificationType(int residueId, ModificationType type) const;
    
    // Get modification count
    int modificationCount() const;
    
    // Get residue count with modifications
    int modifiedResidueCount() const;
    
    // Find modifications at a position
    QList<BiopolymerModification> findModificationsAtPosition(ModificationPosition position) const;
    
    // Find terminal modifications
    QList<BiopolymerModification> findTerminalModifications() const;
    
    // Find crosslinks
    QList<BiopolymerModification> findCrosslinks() const;
    
    // Find cyclizations
    QList<BiopolymerModification> findCyclizations() const;
    
    // Find branch points
    QList<BiopolymerModification> findBranchPoints() const;
    
    // Check if state has unsupported modifications
    bool hasUnsupportedModifications() const;
    
    // Get unsupported modifications
    QList<BiopolymerModification> unsupportedModifications() const;
    
    // Validate all modifications
    QList<QString> validateAll() const;
    
    // ========== Crosslink/Connection Management ==========
    
    // Add a disulfide bond between two residues
    bool addDisulfideBond(int residueId1, int residueId2);
    
    // Remove a disulfide bond
    bool removeDisulfideBond(int residueId1, int residueId2);
    
    // Check if two residues are linked by a disulfide
    bool hasDisulfideBond(int residueId1, int residueId2) const;
    
    // Get all disulfide bonds
    QList<std::pair<int, int>> allDisulfideBonds() const;
    
    // ========== Cyclization Management ==========
    
    // Add a cyclization between two residues
    bool addCyclization(int residueId1, int residueId2);
    
    // Check if the chain is cyclic
    bool isCyclic(int chainId) const;
    
    // Get cyclization points for a chain
    QList<std::pair<int, int>> getCyclizationPoints(int chainId) const;
    
    // ========== Branch Management ==========
    
    // Add a branch point
    bool addBranchPoint(int residueId, const QList<int>& branchResidues);
    
    // Check if a residue is a branch point
    bool isBranchPoint(int residueId) const;
    
    // Get branch residues from a branch point
    QList<int> getBranchResidues(int residueId) const;
    
    // Get all branch points
    QList<int> allBranchPoints() const;

private:
    BiopolymerState& m_state;
    
    // Next modification ID
    int nextModificationId() const;
};

#endif // BIOPOLYMERMODIFICATION_H
