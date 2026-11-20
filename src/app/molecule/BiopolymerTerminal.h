// src/app/molecule/BiopolymerTerminal.h
#ifndef BIOPOLYMERTERMINAL_H
#define BIOPOLYMERTERMINAL_H

#include <QString>
#include <QList>
#include <QMap>
#include <QSet>

// Forward declarations to break circular dependencies
enum class PolymerType;
class BiopolymerState;
class BiopolymerModification;

// ========== Terminal Type ==========

// Types of chain terminals
enum class TerminalType {
    None,           // No terminal group
    N_Terminal,     // N-terminus (peptides)
    C_Terminal,     // C-terminus (peptides)
    Phosphate_5,    // 5' phosphate (DNA/RNA)
    Phosphate_3,    // 3' phosphate (DNA/RNA)
    Hydroxyl_5,     // 5' hydroxyl (DNA/RNA)
    Hydroxyl_3,     // 3' hydroxyl (DNA/RNA)
    Amine,          // Amine group (NH2)
    Acetyl,         // Acetyl group (COCH3)
    Amide,          // Amide group (CONH2)
    Carboxyl,       // Carboxyl group (COOH)
    Methyl,         // Methyl group (CH3)
    Hydroxyl,       // Hydroxyl group (OH)
    Phosphonate,    // Phosphonate group
    Biotin,         // Biotin label
    Fluorescent,    // Fluorescent label
    Custom,         // Custom terminal group
    Unknown         // Unknown terminal type
};

// Convert terminal type to string
QString terminalTypeToString(TerminalType type);

// Parse terminal type from string
TerminalType terminalTypeFromString(const QString& typeString);

// Get all terminal types
QList<TerminalType> allTerminalTypes();

// Check if terminal type is for N-terminus
bool isNTerminal(TerminalType type);

// Check if terminal type is for C-terminus
bool isCTerminal(TerminalType type);

// Check if terminal type is for 5' end
bool is5PrimeTerminal(TerminalType type);

// Check if terminal type is for 3' end
bool is3PrimeTerminal(TerminalType type);

// ========== Terminal Position ==========

// Position of a terminal group
enum class TerminalPosition {
    First,      // First residue in chain (N-term / 5')
    Last,       // Last residue in chain (C-term / 3')
    Both,       // Both ends
    Internal,   // Internal position (not typical for terminals)
    Unknown
};

QString terminalPositionToString(TerminalPosition position);
TerminalPosition terminalPositionFromString(const QString& positionString);

// ========== Biopolymer Terminal Group ==========

// Represents a terminal group on a biopolymer chain
//
// Terminal groups are chemical groups at the ends of biopolymer chains.
// They can affect the chemistry, charge, and reactivity of the biopolymer.
struct BiopolymerTerminal {
    int terminalId = -1;          // Unique terminal ID
    int chainId = -1;             // Chain the terminal belongs to
    TerminalType type = TerminalType::None;
    TerminalPosition position = TerminalPosition::First;
    
    // Terminal identity
    QString name;                 // Common name
    QString shortName;            // Short name/abbreviation
    QString chemicalFormula;      // Chemical formula
    double molecularWeight = 0.0; // Molecular weight added
    
    // Residue association
    int residueId = -1;           // Residue at the terminal (if applicable)
    
    // Validation and support
    bool isSupported = true;       // Whether this terminal is supported
    QStringList validationIssues; // Validation issues
    
    // Metadata
    QString description;          // Description
    QMap<QString, QString> properties; // Additional properties
    
    // ========== Constructors ==========
    
    BiopolymerTerminal();
    BiopolymerTerminal(int termId, int chId, TerminalType termType,
                       TerminalPosition pos, const QString& termName);
    
    // ========== Accessors ==========
    
    // Check if terminal is valid
    bool isValid() const;
    
    // Check if terminal is at N-terminus
    bool isNTerminal() const;
    
    // Check if terminal is at C-terminus
    bool isCTerminal() const;
    
    // Check if terminal is at 5' end
    bool is5Prime() const;
    
    // Check if terminal is at 3' end
    bool is3Prime() const;
    
    // Get full name
    QString fullName() const;
    
    // Get display name
    QString displayName() const;
    
    // ========== Static Factory Methods ==========
    
    // Create an N-terminal acetyl group
    static BiopolymerTerminal createNTerminalAcetyl(int termId, int chainId, int residueId);
    
    // Create an N-terminal amine group
    static BiopolymerTerminal createNTerminalAmine(int termId, int chainId, int residueId);
    
    // Create a C-terminal amide group
    static BiopolymerTerminal createCTerminalAmide(int termId, int chainId, int residueId);
    
    // Create a C-terminal carboxyl group
    static BiopolymerTerminal createCTerminalCarboxyl(int termId, int chainId, int residueId);
    
    // Create a 5' phosphate group
    static BiopolymerTerminal create5PrimePhosphate(int termId, int chainId, int residueId);
    
    // Create a 5' hydroxyl group
    static BiopolymerTerminal create5PrimeHydroxyl(int termId, int chainId, int residueId);
    
    // Create a 3' phosphate group
    static BiopolymerTerminal create3PrimePhosphate(int termId, int chainId, int residueId);
    
    // Create a 3' hydroxyl group
    static BiopolymerTerminal create3PrimeHydroxyl(int termId, int chainId, int residueId);
    
    // ========== Comparison ==========
    
    bool operator==(const BiopolymerTerminal& other) const;
    bool operator!=(const BiopolymerTerminal& other) const;
};

// ========== Terminal Database ==========

// Database of known terminal groups
class TerminalDatabase {
public:
    // Terminal record
    struct TerminalRecord {
        TerminalType type = TerminalType::None;
        QString name;
        QString shortName;
        QString chemicalFormula;
        double molecularWeight = 0.0;
        TerminalPosition position = TerminalPosition::First;
        QSet<PolymerType> applicableTo;  // Which polymer types this applies to
        bool isSupported = true;
        QString description;
    };
    
    // Get the singleton instance
    static TerminalDatabase& instance();
    
    // Get terminal by name
    TerminalRecord getTerminal(const QString& name) const;
    
    // Get terminal by type
    TerminalRecord getTerminal(TerminalType type) const;
    
    // Get all terminals
    QList<TerminalRecord> allTerminals() const;
    
    // Get terminals for a polymer type
    QList<TerminalRecord> terminalsForPolymerType(PolymerType polymerType) const;
    
    // Get N-terminals for a polymer type
    QList<TerminalRecord> nTerminalsForPolymerType(PolymerType polymerType) const;
    
    // Get C-terminals for a polymer type
    QList<TerminalRecord> cTerminalsForPolymerType(PolymerType polymerType) const;
    
    // Check if a terminal is known
    bool isKnownTerminal(const QString& name) const;
    
    // Register a custom terminal
    void registerCustomTerminal(const TerminalRecord& record);
    
    // ========== Common Terminals ==========
    
    // Peptide N-terminals
    static const TerminalRecord N_Terminal_Amine;
    static const TerminalRecord N_Terminal_Acetyl;
    static const TerminalRecord N_Terminal_Methyl;
    static const TerminalRecord N_Terminal_Formyl;
    
    // Peptide C-terminals
    static const TerminalRecord C_Terminal_Carboxyl;
    static const TerminalRecord C_Terminal_Amide;
    static const TerminalRecord C_Terminal_MethylEster;
    
    // DNA/RNA 5' terminals
    static const TerminalRecord DNA_5Prime_Phosphate;
    static const TerminalRecord DNA_5Prime_Hydroxyl;
    
    // DNA/RNA 3' terminals
    static const TerminalRecord DNA_3Prime_Phosphate;
    static const TerminalRecord DNA_3Prime_Hydroxyl;

private:
    TerminalDatabase();
    ~TerminalDatabase();
    
    // Initialize the database
    void initialize();
    
    QMap<QString, TerminalRecord> m_terminalsByName;
    QMap<TerminalType, TerminalRecord> m_terminalsByType;
};

// ========== Terminal Validator ==========

// Validates terminal groups for a biopolymer state
class TerminalValidator {
public:
    explicit TerminalValidator(const BiopolymerState& state);
    
    // Validate a terminal group
    bool validateTerminal(const BiopolymerTerminal& terminal) const;
    
    // Validate all terminal groups
    QList<QString> validateAllTerminals() const;
    
    // Check if a terminal is compatible with the chain
    bool isCompatible(int chainId, const BiopolymerTerminal& terminal) const;
    
    // Check if a chain can have a specific terminal type
    bool canHaveTerminalType(int chainId, TerminalType type, TerminalPosition position) const;
    
    // Check if adding a terminal would conflict
    bool wouldConflict(int chainId, const BiopolymerTerminal& terminal) const;
    
    // Get conflicting terminals
    QList<int> getConflictingTerminals(int chainId, const BiopolymerTerminal& terminal) const;
    
    // Get suggested terminals for a chain
    QList<BiopolymerTerminal> getSuggestedTerminals(int chainId, TerminalPosition position) const;

private:
    const BiopolymerState& m_state;
    
    // Check chain type compatibility
    bool checkChainType(int chainId, TerminalType type) const;
    
    // Check position compatibility
    bool checkPosition(int chainId, TerminalPosition position) const;
};

// ========== Terminal Manager ==========

// Manages terminal groups within a biopolymer state
//
// This class provides high-level operations for working with terminal groups,
// including adding, removing, and querying terminals.
class TerminalManager {
public:
    explicit TerminalManager(BiopolymerState& state);
    
    // Add a terminal group
    bool addTerminal(const BiopolymerTerminal& terminal);
    
    // Add a terminal by type
    int addTerminal(int chainId, TerminalType type, TerminalPosition position = TerminalPosition::First);
    
    // Remove a terminal
    bool removeTerminal(int terminalId);
    
    // Remove all terminals from a chain
    void removeTerminalsFromChain(int chainId);
    
    // Remove all terminals of a type
    void removeTerminalsByType(TerminalType type);
    
    // Get terminal by ID
    BiopolymerTerminal* getTerminal(int terminalId);
    const BiopolymerTerminal* getTerminal(int terminalId) const;
    
    // Get all terminals
    QList<BiopolymerTerminal> allTerminals() const;
    
    // Get terminals for a chain
    QList<BiopolymerTerminal> terminalsForChain(int chainId) const;
    
    // Get terminals for a position
    QList<BiopolymerTerminal> terminalsForPosition(TerminalPosition position) const;
    
    // Get terminals of a type
    QList<BiopolymerTerminal> terminalsByType(TerminalType type) const;
    
    // Get N-terminal for a chain
    BiopolymerTerminal* nTerminal(int chainId);
    const BiopolymerTerminal* nTerminal(int chainId) const;
    
    // Get C-terminal for a chain
    BiopolymerTerminal* cTerminal(int chainId);
    const BiopolymerTerminal* cTerminal(int chainId) const;
    
    // Check if a chain has N-terminal
    bool hasNTerminal(int chainId) const;
    
    // Check if a chain has C-terminal
    bool hasCTerminal(int chainId) const;
    
    // Check if a chain has a specific terminal type
    bool chainHasTerminalType(int chainId, TerminalType type) const;
    
    // Get terminal count
    int terminalCount() const;
    
    // Get chain count with terminals
    int chainsWithTerminalsCount() const;
    
    // Check if state has unsupported terminals
    bool hasUnsupportedTerminals() const;
    
    // Get unsupported terminals
    QList<BiopolymerTerminal> unsupportedTerminals() const;
    
    // Validate all terminals
    QList<QString> validateAll() const;
    
    // ========== Convenience Methods ==========
    
    // Set N-terminal for a chain
    bool setNTerminal(int chainId, TerminalType type);
    
    // Set C-terminal for a chain
    bool setCTerminal(int chainId, TerminalType type);
    
    // Get the default N-terminal for a polymer type
    TerminalType defaultNTerminal(PolymerType polymerType) const;
    
    // Get the default C-terminal for a polymer type
    TerminalType defaultCTerminal(PolymerType polymerType) const;
    
    // Check if all chains have terminals defined
    bool allChainsHaveTerminals() const;
    
    // Get chains missing N-terminals
    QList<int> chainsMissingNTerminals() const;
    
    // Get chains missing C-terminals
    QList<int> chainsMissingCTerminals() const;

private:
    BiopolymerState& m_state;
    
    // Next terminal ID
    int nextTerminalId() const;
};

// ========== Terminal Edit Commands ==========

// For edit commands that integrate with the undo/redo system,
// see BiopolymerTerminalCommands.h which contains commands that
// work with the BiopolymerEditCommand base class.

#endif // BIOPOLYMERTERMINAL_H
