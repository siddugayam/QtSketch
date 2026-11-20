// src/app/molecule/BiopolymerState.h
#ifndef BIOPOLYMERSTATE_H
#define BIOPOLYMERSTATE_H

// Document-owned biopolymer state for QtSketch.
//
// Owns chain, residue, and connection IDs with stable, never-reused identifiers.
// Maintains residue-to-atom mappings using stable AtomId values from EditableMolecule.
// Tracks mapping validity and provides validation explanations for invalid regions.
// Preserves import provenance and declares unsupported features explicitly.
//
// Residue positions are derived from chain order. Renaming or inserting before a
// residue does not change its identity.
//
// Nucleotide representation: each displayed nucleotide is one residue, even when
// the chemistry uses several components (base, sugar, phosphate).
//
// Ownership rules:
// - A nucleotide owns its base, sugar, and incoming backbone phosphate when present.
// - Terminal groups are separately identified and associated with the chain end.
// - Each mapped atom has one primary residue or terminal-group owner.
// - Inter-residue bonds are represented as connections between owned atoms.

#include <QString>
#include <QList>
#include <QHash>
#include <QMap>
#include <QSet>
#include <QDateTime>
#include "EditableMolecule.h"
#include "BiopolymerModification.h"
#include "BiopolymerTerminal.h"

using AtomId = int;
using BondId = int;

// Polymer type enumeration
enum class PolymerType {
    Peptide,
    DNA,
    RNA,
    Unknown
};

// Monomer identity with stereochemistry and attachment info
struct MonomerIdentity {
    QString symbol;           // Single-letter code (e.g., "A", "C", "G", "T", "U")
    QString alias;           // Full name (e.g., "Ala", "dA", "A")
    QString monomerClass;    // "AminoAcid" or "Base"
    bool ambiguous = false;  // True for B, Z, X, N, etc.
    
    // For stereochemistry validation
    QString stereochemistry;  // "L" or "D" for peptides
    
    // For connection validation
    QList<int> attachmentPoints;  // Indices of attachment atoms in the monomer template
};

// Terminal group metadata
struct TerminalGroup {
    enum class Position { NTerminus, CTerminus, Phosphate5Prime, Phosphate3Prime };
    Position position;
    QString name;             // e.g., "NH2", "COOH", "PO4"
    AtomId ownerAtomId = -1;  // Atom that owns this terminal group
};

// Residue in a chain
struct BiopolymerResidue {
    int id = 0;                      // Stable, never-reused residue ID
    int chainId = 0;                 // Owning chain ID
    MonomerIdentity identity;
    int sequencePosition = 0;        // Position in chain (1-indexed)
    
    // Atom ownership: which atoms belong to this residue
    QList<AtomId> ownedAtomIds;
    
    // Connection to next residue (for peptides: peptide bond; for nucleotides: phosphodiester)
    BondId nextBondId = -1;
    
    // Validity tracking
    bool isValid = true;
    QString invalidityReason;       // If !isValid, explains why
};

// Chain of residues
struct BiopolymerChain {
    int id = 0;                      // Stable, never-reused chain ID
    QString name;                    // Chain name (e.g., "A", "B")
    PolymerType polymerType = PolymerType::Unknown;
    
    // Ordered list of residue IDs
    QList<int> residueIds;
    
    // Terminal groups
    TerminalGroup nTerminus;        // For peptides
    TerminalGroup cTerminus;        // For peptides
    TerminalGroup phosphate5Prime;  // For nucleotides
    TerminalGroup phosphate3Prime;  // For nucleotides
    
    // Validity tracking
    bool isValid = true;
    QString invalidityReason;
};

// Connection between residues (peptide bonds, disulfides, phosphodiester bonds)
struct BiopolymerConnection {
    int id = 0;                      // Stable, never-reused connection ID
    int fromResidueId = 0;
    int toResidueId = 0;
    BondId bondId = -1;              // The actual bond in EditableMolecule
    QString type;                   // "peptide", "phosphodiester", "disulfide", etc.
};

// Connection types for extended connections
enum class ConnectionType {
    Peptide,
    Phosphodiester,
    Disulfide,
    Custom,
    Unknown
};

// Mapping validity result
struct MappingValidity {
    bool isValid = true;
    QList<QString> issues;          // List of validation issues
    
    // Quick check methods
    bool hasIssues() const { return !issues.isEmpty(); }
    QString summary() const { return issues.join(", "); }
};

// Import provenance tracking
struct ImportProvenance {
    QString sourceFormat;           // "FASTA", "Sequence", "HELM", "IDT", "AxoLabs", "MOL", "SDF"
    QString sourceText;
    QDateTime importTime;
    
    // Unsupported features detected during import
    QSet<QString> unsupportedFeatures;
    
    // Whether the import was fully mapped
    bool fullyMapped = false;
};

// Main biopolymer state class
class BiopolymerState {
public:
    BiopolymerState();
    explicit BiopolymerState(const BiopolymerState& other);
    BiopolymerState& operator=(const BiopolymerState& other);
    ~BiopolymerState();

    // ========== Chain Management ==========
    
    // Create a new chain with a unique ID
    int addChain(const QString& name, PolymerType type);
    
    // Remove a chain and all its residues
    bool removeChain(int chainId);
    
    // Rename a chain (does not change identity)
    bool renameChain(int chainId, const QString& newName);
    
    // Get chain by ID
    BiopolymerChain* chain(int chainId);
    const BiopolymerChain* chain(int chainId) const;
    
    // Get all chain IDs
    QList<int> chainIds() const;
    
    // Get chain count
    int chainCount() const;
    
    // ========== Residue Management ==========
    
    // Add a residue to a chain
    int addResidue(int chainId, const MonomerIdentity& identity, int position = -1);
    
    // Insert a residue before an existing residue in a chain
    int insertResidueBefore(int chainId, int beforeResidueId, const MonomerIdentity& identity);
    
    // Insert a residue after an existing residue in a chain
    int insertResidueAfter(int chainId, int afterResidueId, const MonomerIdentity& identity);
    
    // Replace a residue
    bool replaceResidue(int residueId, const MonomerIdentity& identity);
    
    // Delete a residue from a chain
    bool deleteResidue(int residueId);
    
    // Delete a range of residues from a chain
    bool deleteResidueRange(int chainId, int startPosition, int endPosition);
    
    // Get residue by ID
    BiopolymerResidue* residue(int residueId);
    const BiopolymerResidue* residue(int residueId) const;
    
    // Get residue by chain and position
    BiopolymerResidue* residueInChain(int chainId, int position);
    const BiopolymerResidue* residueInChain(int chainId, int position) const;
    
    // Get all residue IDs in a chain
    QList<int> residueIdsInChain(int chainId) const;
    
    // Get all residues
    QList<int> allResidueIds() const;
    
    // ========== Atom Mapping ==========
    
    // Set the atoms owned by a residue
    void setResidueAtoms(int residueId, const QList<AtomId>& atomIds);
    
    // Add atoms to a residue's ownership
    void addResidueAtoms(int residueId, const QList<AtomId>& atomIds);
    
    // Remove atoms from a residue's ownership
    void removeResidueAtoms(int residueId, const QList<AtomId>& atomIds);
    
    // Get atoms owned by a residue
    QList<AtomId> atomsOwnedByResidue(int residueId) const;
    
    // Find which residue owns an atom (returns -1 if not mapped)
    int residueOwningAtom(AtomId atomId) const;
    
    // Set the bond connecting two residues
    void setResidueConnection(int fromResidueId, int toResidueId, BondId bondId, const QString& type);
    
    // Get the bond connecting two residues
    BondId connectionBond(int fromResidueId, int toResidueId) const;
    
    // Get all connections
    QList<BiopolymerConnection> allConnections() const;
    
    // ========== Terminal Groups ==========
    
    // Set N-terminus for a peptide chain
    void setNTerminus(int chainId, const TerminalGroup& terminus);
    
    // Set C-terminus for a peptide chain
    void setCTerminus(int chainId, const TerminalGroup& terminus);
    
    // Set 5' phosphate for a nucleotide chain
    void setPhosphate5Prime(int chainId, const TerminalGroup& phosphate);
    
    // Set 3' phosphate for a nucleotide chain
    void setPhosphate3Prime(int chainId, const TerminalGroup& phosphate);
    
    // Get terminal groups for a chain
    TerminalGroup getNTerminus(int chainId) const;
    TerminalGroup getCTerminus(int chainId) const;
    TerminalGroup getPhosphate5Prime(int chainId) const;
    TerminalGroup getPhosphate3Prime(int chainId) const;

    // ========== Biopolymer Terminals (Extended) ==========

    // Add a terminal to the state
    bool addTerminal(const BiopolymerTerminal& terminal);

    // Remove a terminal by ID
    bool removeTerminal(int terminalId);

    // Remove all terminals from a chain
    void removeTerminalsFromChain(int chainId);

    // Remove terminals by type
    void removeTerminalsByType(TerminalType type);

    // Get terminal by ID
    BiopolymerTerminal* getTerminal(int terminalId);
    const BiopolymerTerminal* getTerminal(int terminalId) const;

    // Get all terminals
    QList<BiopolymerTerminal> terminals() const;

    // Get terminals for a specific chain
    QList<BiopolymerTerminal> terminalsForChain(int chainId) const;

    // Get terminals by type
    QList<BiopolymerTerminal> terminalsByType(TerminalType type) const;

    // Allocate a new unique terminal ID
    int allocateTerminalId();

    // ========== Biopolymer Modifications ==========

    // Add a modification to the state
    bool addModification(const BiopolymerModification& modification);

    // Remove a modification by ID
    bool removeModification(int modificationId);

    // Remove all modifications from a residue
    void removeModificationsFromResidue(int residueId);

    // Remove all modifications from a chain
    void removeModificationsFromChain(int chainId);

    // Remove modifications by type
    void removeModificationsByType(ModificationType type);

    // Get modification by ID
    BiopolymerModification* getModification(int modificationId);
    const BiopolymerModification* getModification(int modificationId) const;

    // Get all modifications
    QList<BiopolymerModification> modifications() const;

    // Get modifications for a specific residue
    QList<BiopolymerModification> modificationsForResidue(int residueId) const;

    // Get modifications for a specific chain
    QList<BiopolymerModification> modificationsForChain(int chainId) const;

    // Get modifications by type
    QList<BiopolymerModification> modificationsByType(ModificationType type) const;

    // Check if a residue has modifications
    bool hasModifications(int residueId) const;

    // Check if a chain has modifications
    bool hasChainModifications(int chainId) const;

    // Allocate a new unique modification ID
    int allocateModificationId();

    // Get next ID counters including modifications and terminals
    struct AllIdCounters {
        int nextResidueId = 1;
        int nextChainId = 1;
        int nextConnectionId = 1;
        int nextTerminalId = 1;
        int nextModificationId = 1;
    };
    AllIdCounters allIdCounters() const;
    
    // ========== Validation ==========
    
    // Validate all mappings
    MappingValidity validateAllMappings() const;
    
    // Validate a specific residue's mapping
    MappingValidity validateResidueMapping(int residueId) const;
    
    // Validate a chain
    MappingValidity validateChain(int chainId) const;
    
    // Check if a residue is fully mapped
    bool isResidueFullyMapped(int residueId) const;
    
    // Check if a chain is fully mapped
    bool isChainFullyMapped(int chainId) const;
    
    // ========== Glue Methods for Export/Validation ==========
    
    // Check if all mappings are valid
    bool hasValidMappings() const;
    
    // Check if any mappings are invalid
    bool hasInvalidMappings() const;
    
    // Check if there are unsupported modifications
    bool hasUnsupportedModifications() const;
    
    // Check if there are unsupported connections
    bool hasUnsupportedConnections() const;
    
    // Check if there are terminal groups
    bool hasTerminalGroups() const;
    
    // Check if there is branched topology
    bool hasBranchedTopology() const;
    
    // Get total residue count
    int totalResidueCount() const;
    
    // Get polymer types present in the state
    QList<PolymerType> polymerTypes() const;
    
    // Get residue position in chain
    int residuePositionInChain(int chainId, int residueId) const;
    
    // Check if residue is N-terminal
    bool isNTerminal(int chainId, int residueId) const;
    
    // Check if residue is C-terminal
    bool isCTerminal(int chainId, int residueId) const;
    
    // Get chain residue count
    int chainResidueCount(int chainId) const;
    
    // ========== ID Management ==========
    
    // Allocate a new unique residue ID
    int allocateResidueId();
    
    // Allocate a new unique chain ID
    int allocateChainId();
    
    // Allocate a new unique connection ID
    int allocateConnectionId();
    
    // Reset ID counters (used for deserialization)
    void resetIdCounters(int nextResidueId, int nextChainId, int nextConnectionId,
                        int nextTerminalId = 1, int nextModificationId = 1);
    
    // Get current ID counters (for serialization)
    struct IdCounters {
        int nextResidueId = 1;
        int nextChainId = 1;
        int nextConnectionId = 1;
    };
    IdCounters idCounters() const;
    
    // ========== Import/Export ==========
    
    // Set import provenance
    void setImportProvenance(const ImportProvenance& provenance);
    
    // Get import provenance
    ImportProvenance getImportProvenance() const;
    
    // Clear import provenance
    void clearImportProvenance();
    
    // ========== State Query ==========
    
    // Check if any biopolymer data exists
    bool hasBiopolymerData() const;
    
    // Check if all chains are valid
    bool isFullyValid() const;
    
    // Get all validation issues
    QList<QString> allValidationIssues() const;
    
    // ========== Selection Support ==========
    
    // Get residues from a set of atom IDs
    QSet<int> residuesFromAtoms(const QSet<AtomId>& atomIds) const;
    
    // Get atoms from a set of residue IDs
    QSet<AtomId> atomsFromResidues(const QSet<int>& residueIds) const;
    
    // ========== Polymer Type Support ==========
    
    // Check if a polymer type is supported
    static bool isPolymerTypeSupported(PolymerType type);
    
    // Get string representation of polymer type
    static QString polymerTypeToString(PolymerType type);
    
    // Parse polymer type from string
    static PolymerType polymerTypeFromString(const QString& str);
    
    // ========== Serialization Support ==========
    
    // Prepare for serialization (clear any runtime-only data)
    void prepareForSerialization();
    
    // Finalize after deserialization (rebuild any runtime data)
    void finalizeAfterDeserialization();

private:
    // ID counters (monotonic, never reused)
    int m_nextResidueId = 1;
    int m_nextChainId = 1;
    int m_nextConnectionId = 1;
    int m_nextTerminalId = 1;
    int m_nextModificationId = 1;
    
    // Data storage
    QMap<int, BiopolymerChain> m_chains;
    QMap<int, BiopolymerResidue> m_residues;
    QList<BiopolymerConnection> m_connections;
    
    // Extended terminal storage (using BiopolymerTerminal)
    QMap<int, BiopolymerTerminal> m_terminals;
    
    // Modification storage
    QMap<int, BiopolymerModification> m_modifications;
    
    // Atom-to-residue mapping (for fast lookup)
    QHash<AtomId, int> m_atomToResidue;
    
    // Import provenance
    ImportProvenance m_importProvenance;
    
    // Validation cache
    mutable QHash<int, MappingValidity> m_validationCache;
    
    // Helper methods
    BiopolymerChain* findChain(int chainId);
    const BiopolymerChain* findChain(int chainId) const;
    
    void rebuildAtomToResidueMap();
    void invalidateValidationCache() const;
};

#endif // BIOPOLYMERSTATE_H
