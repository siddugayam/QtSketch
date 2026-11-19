// src/app/molecule/BiopolymerConnectionExtended.h
#ifndef BIOPOLYMERCONNECTIONEXTENDED_H
#define BIOPOLYMERCONNECTIONEXTENDED_H

#include <QString>
#include <QList>
#include <QHash>
#include <QMap>
#include <QSet>
#include "BiopolymerState.h"

// Connection class for extended biopolymer connections
// Handles disulfides, interchain connections, cyclization, branches

// Connection subtype enumeration
enum class ConnectionSubtype {
    Peptide,            // Standard peptide bond
    Phosphodiester,     // Standard phosphodiester bond
    Disulfide,          // Disulfide bond (cysteine-cysteine)
    Thioester,         // Thioester bond
    Lactam,            // Lactam bridge
    Custom,            // Custom connection type
    Cyclization,       // Cyclization connection
    Branch,            // Branch point connection
    Interchain,        // Interchain connection
    Unknown
};

// Disulfide bond representation
struct BiopolymerDisulfide {
    int disulfideId = 0;              // Unique disulfide ID
    int residueId1 = 0;               // First cysteine residue ID
    int residueId2 = 0;               // Second cysteine residue ID
    int chainId1 = 0;                 // Chain ID of first residue
    int chainId2 = 0;                 // Chain ID of second residue
    
    // Atom-level details
    int atomId1 = -1;                 // Sulfur atom ID in first residue
    int atomId2 = -1;                 // Sulfur atom ID in second residue
    BondId bondId = -1;               // The actual bond ID in EditableMolecule
    
    // Bond properties
    QString type = "Disulfide";      // "Disulfide", "Diselenide", etc.
    double bondLength = 2.05;         // Typical S-S bond length in Angstroms
    bool isIntrachain = false;        // Whether both residues are in the same chain
    bool isValid = true;              // Whether the disulfide is chemically valid
    QString validationError;         // If invalid, explains why
    
    // Metadata
    QString name;                     // Name of the disulfide
    QString description;              // Description
    bool isSelected = false;          // Whether this disulfide is selected
    bool isVisible = true;            // Whether this disulfide is visible
    
    // Timestamp for undo/redo
    quint64 creationTimestamp = 0;
    
    // Quick check methods
    bool isValidDisulfide() const { 
        return isValid && residueId1 > 0 && residueId2 > 0 && 
               residueId1 != residueId2 && bondId >= 0; 
    }
    
    bool isInterchain() const {
        return chainId1 != chainId2;
    }
    
    bool connects(int r1, int r2) const {
        return (residueId1 == r1 && residueId2 == r2) || 
               (residueId1 == r2 && residueId2 == r1);
    }
    
    bool involvesResidue(int residueId) const {
        return residueId1 == residueId || residueId2 == residueId;
    }
    
    bool involvesChain(int chainId) const {
        return chainId1 == chainId || chainId2 == chainId;
    }
    
    int otherResidue(int residueId) const {
        if (residueId1 == residueId) return residueId2;
        if (residueId2 == residueId) return residueId1;
        return -1;
    }
    
    int otherChain(int chainId) const {
        if (chainId1 == chainId) return chainId2;
        if (chainId2 == chainId) return chainId1;
        return -1;
    }
};

// Interchain connection representation (between different chains)
struct BiopolymerInterchainConnection {
    int connectionId = 0;             // Unique connection ID
    int chainId1 = 0;                 // First chain ID
    int chainId2 = 0;                 // Second chain ID
    int residueId1 = 0;               // Residue ID in first chain
    int residueId2 = 0;               // Residue ID in second chain
    
    // Atom-level details
    int atomId1 = -1;                 // Atom ID in first residue
    int atomId2 = -1;                 // Atom ID in second residue
    BondId bondId = -1;               // The actual bond ID in EditableMolecule
    
    // Connection properties
    ConnectionSubtype subtype = ConnectionSubtype::Interchain;
    QString type = "Interchain";      // Type of connection
    double bondLength = 0.0;          // Bond length
    bool isValid = true;              // Whether the connection is valid
    QString validationError;         // Validation error message
    
    // Metadata
    QString name;                     // Name of the connection
    QString description;              // Description
    bool isSelected = false;
    bool isVisible = true;
    
    // Timestamp for undo/redo
    quint64 creationTimestamp = 0;
    
    // Quick check methods
    bool isValidConnection() const {
        return isValid && chainId1 > 0 && chainId2 > 0 && 
               chainId1 != chainId2 && residueId1 > 0 && residueId2 > 0;
    }
    
    bool connectsChains(int c1, int c2) const {
        return (chainId1 == c1 && chainId2 == c2) || 
               (chainId1 == c2 && chainId2 == c1);
    }
    
    bool involvesChain(int chainId) const {
        return chainId1 == chainId || chainId2 == chainId;
    }
    
    bool involvesResidue(int residueId) const {
        return residueId1 == residueId || residueId2 == residueId;
    }
};

// Branch point representation (for branched biopolymers)
struct BiopolymerBranch {
    int branchId = 0;                 // Unique branch ID
    int chainId = 0;                 // Chain ID where branch originates
    int residueId = 0;               // Residue ID where branch originates
    int branchChainId = 0;           // Chain ID of the branch
    int branchResidueId = 0;         // First residue ID in the branch
    
    // Connection details
    int atomId = -1;                  // Atom ID where branch attaches
    BondId bondId = -1;               // Bond ID connecting branch
    
    // Branch properties
    QString type = "Branch";          // Type of branch
    QString name;                     // Name of the branch
    QString description;              // Description
    bool isValid = true;              // Whether the branch is valid
    QString validationError;         // Validation error
    
    // Geometry
    bool isStereocenter = false;     // Whether the branch point is a stereocenter
    QString stereochemistry;          // Stereochemistry if applicable
    
    // Metadata
    bool isSelected = false;
    bool isVisible = true;
    quint64 creationTimestamp = 0;
    
    // Quick check methods
    bool isValidBranch() const {
        return isValid && chainId > 0 && residueId > 0 && 
               branchChainId > 0 && branchResidueId > 0 &&
               bondId >= 0;
    }
    
    bool isAttachedTo(int residueId) const {
        return this->residueId == residueId || branchResidueId == residueId;
    }
    
    bool isInChain(int chainId) const {
        return this->chainId == chainId || branchChainId == chainId;
    }
};

// Cyclization representation
struct BiopolymerCyclization {
    int cyclizationId = 0;           // Unique cyclization ID
    QList<int> residueIds;            // Residues involved in cyclization
    QList<int> chainIds;              // Chains involved in cyclization
    
    // Connection details
    QList<BondId> bondIds;            // Bond IDs forming the cycle
    QList<int> atomIds;              // Atom IDs involved in the cycle
    
    // Cyclization properties
    ConnectionSubtype subtype = ConnectionSubtype::Cyclization;
    QString type = "Cyclization";    // Type of cyclization
    int ringSize = 0;                 // Number of atoms in the ring
    QString name;                     // Name of the cyclization
    QString description;              // Description
    bool isValid = true;              // Whether the cyclization is valid
    QString validationError;         // Validation error
    
    // Geometry
    bool isPlanar = false;            // Whether the ring is planar
    bool isAromatic = false;          // Whether the ring is aromatic
    
    // Metadata
    bool isSelected = false;
    bool isVisible = true;
    quint64 creationTimestamp = 0;
    
    // Quick check methods
    bool isValidCyclization() const {
        return isValid && residueIds.size() >= 3 && bondIds.size() >= 3;
    }
    
    bool involvesResidue(int residueId) const {
        return residueIds.contains(residueId);
    }
    
    bool involvesChain(int chainId) const {
        return chainIds.contains(chainId);
    }
};

// Unified extended connection class
class BiopolymerExtendedConnection {
public:
    int id = 0;                      // Unique connection ID
    ConnectionSubtype subtype = ConnectionSubtype::Unknown;
    
    // Connection data (one of these will be valid based on subtype)
    BiopolymerDisulfide disulfide;
    BiopolymerInterchainConnection interchain;
    BiopolymerBranch branch;
    BiopolymerCyclization cyclization;
    
    // Common properties
    QString type;                    // Connection type string
    QString name;                    // Name
    QString description;             // Description
    bool isValid = true;             // Whether the connection is valid
    QString validationError;        // Validation error
    bool isSelected = false;
    bool isVisible = true;
    quint64 creationTimestamp = 0;
    
    // Atom and bond references
    QList<int> residueIds;           // All residues involved
    QList<int> chainIds;             // All chains involved
    QList<int> atomIds;              // All atoms involved
    QList<BondId> bondIds;           // All bonds involved
    
    // Constructor
    BiopolymerExtendedConnection() = default;
    BiopolymerExtendedConnection(int connectionId, ConnectionSubtype subType)
        : id(connectionId), subtype(subType) {}
    
    // Factory methods
    static BiopolymerExtendedConnection createDisulfide(int connId, const BiopolymerDisulfide& ds);
    static BiopolymerExtendedConnection createInterchain(int connId, const BiopolymerInterchainConnection& ic);
    static BiopolymerExtendedConnection createBranch(int connId, const BiopolymerBranch& br);
    static BiopolymerExtendedConnection createCyclization(int connId, const BiopolymerCyclization& cyc);
    
    // Check methods
    bool isDisulfide() const { return subtype == ConnectionSubtype::Disulfide; }
    bool isInterchain() const { return subtype == ConnectionSubtype::Interchain; }
    bool isBranch() const { return subtype == ConnectionSubtype::Branch; }
    bool isCyclization() const { return subtype == ConnectionSubtype::Cyclization; }
    
    bool isValidConnection() const;
    bool involvesResidue(int residueId) const;
    bool involvesChain(int chainId) const;
    bool involvesAtom(int atomId) const;
    
    // Accessors
    const BiopolymerDisulfide& asDisulfide() const;
    BiopolymerDisulfide& asDisulfide();
    const BiopolymerInterchainConnection& asInterchain() const;
    BiopolymerInterchainConnection& asInterchain();
    const BiopolymerBranch& asBranch() const;
    BiopolymerBranch& asBranch();
    const BiopolymerCyclization& asCyclization() const;
    BiopolymerCyclization& asCyclization();
};

// Connection validator for extended connections
class ExtendedConnectionValidator {
public:
    explicit ExtendedConnectionValidator(const BiopolymerState& state);
    ~ExtendedConnectionValidator();
    
    // Validation
    bool validateDisulfide(const BiopolymerDisulfide& disulfide) const;
    bool validateInterchain(const BiopolymerInterchainConnection& connection) const;
    bool validateBranch(const BiopolymerBranch& branch) const;
    bool validateCyclization(const BiopolymerCyclization& cyclization) const;
    bool validateExtendedConnection(const BiopolymerExtendedConnection& connection) const;
    
    // Batch validation
    QList<QString> validateAllDisulfides() const;
    QList<QString> validateAllInterchainConnections() const;
    QList<QString> validateAllBranches() const;
    QList<QString> validateAllCyclizations() const;
    QList<QString> validateAllExtendedConnections() const;
    
    // Compatibility checks
    bool canFormDisulfide(int residueId1, int residueId2) const;
    bool canFormInterchainConnection(int residueId1, int residueId2) const;
    bool canCreateBranch(int chainId, int residueId, int branchChainId) const;
    bool canCreateCyclization(const QList<int>& residueIds) const;
    
    // Get suggestions
    QList<BiopolymerDisulfide> getPossibleDisulfides(int chainId) const;
    QList<BiopolymerInterchainConnection> getPossibleInterchainConnections() const;
    
    // Check for conflicts
    bool wouldConflict(const BiopolymerDisulfide& disulfide) const;
    bool wouldConflict(const BiopolymerInterchainConnection& connection) const;
    bool wouldConflict(const BiopolymerBranch& branch) const;
    bool wouldConflict(const BiopolymerCyclization& cyclization) const;
    
    // Get conflicting connections
    QList<int> getConflictingDisulfides(const BiopolymerDisulfide& disulfide) const;
    QList<int> getConflictingInterchainConnections(const BiopolymerInterchainConnection& connection) const;
    QList<int> getConflictingBranches(const BiopolymerBranch& branch) const;
    
    // Utility
    bool areResiduesCompatibleForDisulfide(int residueId1, int residueId2) const;
    bool areResiduesCompatibleForConnection(int residueId1, int residueId2) const;
    
private:
    const BiopolymerState& m_state;
    
    bool checkResidueExistence(int residueId) const;
    bool checkChainExistence(int chainId) const;
    bool checkAtomExistence(int atomId) const;
    bool checkBondExistence(BondId bondId) const;
    bool checkResidueType(int residueId, const QSet<QString>& allowedTypes) const;
};

// Connection manager for extended connections
class ExtendedConnectionManager {
public:
    explicit ExtendedConnectionManager(BiopolymerState& state);
    ~ExtendedConnectionManager();
    
    // ========== Disulfide Management ==========
    
    // Add a disulfide bond
    int addDisulfide(const BiopolymerDisulfide& disulfide);
    int addDisulfide(int residueId1, int residueId2, BondId bondId, int atomId1 = -1, int atomId2 = -1);
    
    // Remove a disulfide bond
    bool removeDisulfide(int disulfidId);
    
    // Remove disulfides involving a residue
    void removeDisulfidesFromResidue(int residueId);
    
    // Remove disulfides from a chain
    void removeDisulfidesFromChain(int chainId);
    
    // Get disulfide by ID
    BiopolymerDisulfide* getDisulfide(int disulfidId);
    const BiopolymerDisulfide* getDisulfide(int disulfidId) const;
    
    // Get all disulfides
    QList<BiopolymerDisulfide> allDisulfides() const;
    
    // Get disulfides for a residue
    QList<BiopolymerDisulfide> disulfidesForResidue(int residueId) const;
    
    // Get disulfides for a chain
    QList<BiopolymerDisulfide> disulfidesForChain(int chainId) const;
    
    // Check if a residue has disulfides
    bool hasDisulfides(int residueId) const;
    
    // Check if a chain has disulfides
    bool hasChainDisulfides(int chainId) const;
    
    // Get interchain disulfides
    QList<BiopolymerDisulfide> interchainDisulfides() const;
    
    // Get intrachain disulfides
    QList<BiopolymerDisulfide> intrachainDisulfides() const;
    
    // Find disulfide between two residues
    BiopolymerDisulfide* findDisulfide(int residueId1, int residueId2);
    const BiopolymerDisulfide* findDisulfide(int residueId1, int residueId2) const;
    
    // ========== Interchain Connection Management ==========
    
    // Add an interchain connection
    int addInterchainConnection(const BiopolymerInterchainConnection& connection);
    int addInterchainConnection(int chainId1, int residueId1, int chainId2, int residueId2, 
                                BondId bondId, int atomId1 = -1, int atomId2 = -1);
    
    // Remove an interchain connection
    bool removeInterchainConnection(int connectionId);
    
    // Remove interchain connections involving a chain
    void removeInterchainConnectionsFromChain(int chainId);
    
    // Remove interchain connections involving a residue
    void removeInterchainConnectionsFromResidue(int residueId);
    
    // Get interchain connection by ID
    BiopolymerInterchainConnection* getInterchainConnection(int connectionId);
    const BiopolymerInterchainConnection* getInterchainConnection(int connectionId) const;
    
    // Get all interchain connections
    QList<BiopolymerInterchainConnection> allInterchainConnections() const;
    
    // Get interchain connections for a chain
    QList<BiopolymerInterchainConnection> interchainConnectionsForChain(int chainId) const;
    
    // Get interchain connections between two chains
    QList<BiopolymerInterchainConnection> interchainConnectionsBetween(int chainId1, int chainId2) const;
    
    // Check if chains are connected
    bool areChainsConnected(int chainId1, int chainId2) const;
    
    // Find interchain connection between two residues
    BiopolymerInterchainConnection* findInterchainConnection(int residueId1, int residueId2);
    const BiopolymerInterchainConnection* findInterchainConnection(int residueId1, int residueId2) const;
    
    // ========== Branch Management ==========
    
    // Add a branch
    int addBranch(const BiopolymerBranch& branch);
    int addBranch(int chainId, int residueId, int branchChainId, int branchResidueId,
                  BondId bondId, int atomId = -1);
    
    // Remove a branch
    bool removeBranch(int branchId);
    
    // Remove branches from a chain
    void removeBranchesFromChain(int chainId);
    
    // Remove branches from a residue
    void removeBranchesFromResidue(int residueId);
    
    // Get branch by ID
    BiopolymerBranch* getBranch(int branchId);
    const BiopolymerBranch* getBranch(int branchId) const;
    
    // Get all branches
    QList<BiopolymerBranch> allBranches() const;
    
    // Get branches for a chain
    QList<BiopolymerBranch> branchesForChain(int chainId) const;
    
    // Get branches from a residue
    QList<BiopolymerBranch> branchesFromResidue(int residueId) const;
    
    // Check if a residue is a branch point
    bool isBranchPoint(int residueId) const;
    
    // Check if a chain has branches
    bool hasChainBranches(int chainId) const;
    
    // ========== Cyclization Management ==========
    
    // Add a cyclization
    int addCyclization(const BiopolymerCyclization& cyclization);
    
    // Remove a cyclization
    bool removeCyclization(int cyclizationId);
    
    // Remove cyclizations involving a residue
    void removeCyclizationsFromResidue(int residueId);
    
    // Get cyclization by ID
    BiopolymerCyclization* getCyclization(int cyclizationId);
    const BiopolymerCyclization* getCyclization(int cyclizationId) const;
    
    // Get all cyclizations
    QList<BiopolymerCyclization> allCyclizations() const;
    
    // Get cyclizations for a residue
    QList<BiopolymerCyclization> cyclizationsForResidue(int residueId) const;
    
    // Check if a residue is in a cyclization
    bool isInCyclization(int residueId) const;
    
    // ========== Extended Connection Management ==========
    
    // Add an extended connection
    int addExtendedConnection(const BiopolymerExtendedConnection& connection);
    
    // Remove an extended connection
    bool removeExtendedConnection(int connectionId);
    
    // Get extended connection by ID
    BiopolymerExtendedConnection* getExtendedConnection(int connectionId);
    const BiopolymerExtendedConnection* getExtendedConnection(int connectionId) const;
    
    // Get all extended connections
    QList<BiopolymerExtendedConnection> allExtendedConnections() const;
    
    // Get extended connections by type
    QList<BiopolymerExtendedConnection> extendedConnectionsByType(ConnectionSubtype subtype) const;
    
    // ========== Query Methods ==========
    
    // Get total counts
    int disulfideCount() const;
    int interchainConnectionCount() const;
    int branchCount() const;
    int cyclizationCount() const;
    int extendedConnectionCount() const;
    
    // Check for any extended connections
    bool hasExtendedConnections() const;
    bool hasDisulfides() const;
    bool hasInterchainConnections() const;
    bool hasBranches() const;
    bool hasCyclizations() const;
    
    // Check for interchain connections
    bool hasInterchainConnectionsBetweenChains() const;
    
    // ========== Validation ==========
    
    // Validate all extended connections
    QList<QString> validateAll() const;
    
    // Check if all disulfides are valid
    bool areAllDisulfidesValid() const;
    
    // Check if all interchain connections are valid
    bool areAllInterchainConnectionsValid() const;
    
    // ========== Utility ==========
    
    // Get all connections involving a residue
    QList<BiopolymerExtendedConnection> connectionsForResidue(int residueId) const;
    
    // Get all connections involving a chain
    QList<BiopolymerExtendedConnection> connectionsForChain(int chainId) const;
    
    // Check if residue is involved in any extended connection
    bool isResidueInExtendedConnection(int residueId) const;
    
    // ========== ID Management ==========
    
    // Allocate new IDs
    int allocateDisulfideId();
    int allocateInterchainConnectionId();
    int allocateBranchId();
    int allocateCyclizationId();
    int allocateExtendedConnectionId();
    
    // Get next ID values
    int nextDisulfideId() const;
    int nextInterchainConnectionId() const;
    int nextBranchId() const;
    int nextCyclizationId() const;
    int nextExtendedConnectionId() const;
    
private:
    BiopolymerState& m_state;
    
    // Storage for extended connections
    QMap<int, BiopolymerDisulfide> m_disulfides;
    QMap<int, BiopolymerInterchainConnection> m_interchainConnections;
    QMap<int, BiopolymerBranch> m_branches;
    QMap<int, BiopolymerCyclization> m_cyclizations;
    QMap<int, BiopolymerExtendedConnection> m_extendedConnections;
    
    // ID counters
    int m_nextDisulfideId = 1;
    int m_nextInterchainConnectionId = 1;
    int m_nextBranchId = 1;
    int m_nextCyclizationId = 1;
    int m_nextExtendedConnectionId = 1;
};


#endif // BIOPOLYMERCONNECTIONEXTENDED_H
