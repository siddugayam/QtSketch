// src/app/molecule/BiopolymerEditCommandWrapper.h
#ifndef BIOPOLYMEREDITCOMMANDWRAPPER_H
#define BIOPOLYMEREDITCOMMANDWRAPPER_H

// Wrapper to adapt BiopolymerEditCommand to DocumentState's EditCommand system.
//
// BiopolymerEditCommand operates on (EditableMolecule&, BiopolymerState&) pairs,
// while DocumentState's EditCommand uses capture-by-reference lambdas.
// This wrapper bridges the two systems.

#include "BiopolymerEditCommands.h"
#include "EditCommand.h"
#include <memory>

// Forward declarations
class DocumentState;

// Wrapper class that adapts a BiopolymerEditCommand to an EditCommand
class BiopolymerEditCommandWrapper {
public:
    // Create a wrapper for a biopolymer edit command
    static EditCommand wrap(BiopolymerEditCommand* cmd, DocumentState& doc);
    
    // Create and wrap a command in one step
    template <typename T, typename... Args>
    static EditCommand createAndWrap(DocumentState& doc, Args&&... args) {
        return wrap(new T(std::forward<Args>(args)...), doc);
    }

    // Constructor
    explicit BiopolymerEditCommandWrapper(BiopolymerEditCommand* command, DocumentState& document);

    // Destructor
    ~BiopolymerEditCommandWrapper();

    // EditCommand interface
    void execute();
    void unexecute();
    bool canUndo() const;
    QString description() const;
    bool isSignificant() const;
    void setSignificant(bool significant);
    BiopolymerEditCommand* command() const;
    BiopolymerEditCommand* takeCommand();
    
private:
    std::shared_ptr<BiopolymerEditCommand> m_command;
    DocumentState& m_document;
    
    // Helper methods for common operations
    
    // Add a chain
    static EditCommand makeAddChainCommand(DocumentState& doc, int chainId, 
                                          const QString& name, const QString& type);
    
    // Remove a chain
    static EditCommand makeRemoveChainCommand(DocumentState& doc, int chainId, 
                                            const QString& name, const QString& type,
                                            const QList<int>& residues, 
                                            const QList<AtomId>& atoms, 
                                            const QList<BondId>& bonds);
    
    // Rename a chain
    static EditCommand makeRenameChainCommand(DocumentState& doc, int chainId,
                                            const QString& oldName, const QString& newName);
    
    // Add a residue
    static EditCommand makeAddResidueCommand(DocumentState& doc, int residueId, int chainId,
                                           int position, const QString& symbol,
                                           const QString& type,
                                           const QList<AtomId>& atoms, 
                                           const QList<BondId>& bonds);
    
    // Insert a residue
    static EditCommand makeInsertResidueCommand(DocumentState& doc, int residueId, int chainId,
                                              int beforeResidueId, const QString& symbol,
                                              const QString& type,
                                              const QList<AtomId>& atoms, 
                                              const QList<BondId>& bonds,
                                              BondId newBond, BondId oldBond);
    
    // Replace a residue
    static EditCommand makeReplaceResidueCommand(DocumentState& doc, int residueId,
                                               const QString& oldSymbol,
                                               const QString& newSymbol, const QString& type,
                                               const QList<AtomId>& oldAtoms, 
                                               const QList<AtomId>& newAtoms,
                                               const QList<BondId>& oldBonds, 
                                               const QList<BondId>& newBonds,
                                               const QMap<AtomId, AtomId>& atomMap);
    
    // Delete a residue
    static EditCommand makeDeleteResidueCommand(DocumentState& doc, int residueId, int chainId,
                                              int position, const QString& symbol,
                                              const QString& type,
                                              const QList<AtomId>& atoms, 
                                              const QList<BondId>& bonds,
                                              BondId bondToRecreate, int prevRes, int nextRes);
    
    // Delete a residue range
    static EditCommand makeDeleteResidueRangeCommand(DocumentState& doc, int chainId, int start,
                                                    int end,
                                                    const QList<int>& residueIds,
                                                    const QList<QString>& symbols,
                                                    const QString& type,
                                                    const QList<AtomId>& atoms,
                                                    const QList<BondId>& bonds,
                                                    BondId bondToRecreate);
    
    // Connect residues
    static EditCommand makeConnectResiduesCommand(DocumentState& doc, int fromRes, int toRes,
                                                  BondId bond, AtomId fromAtom, AtomId toAtom,
                                                  const QString& type);
    
    // Set residue atoms
    static EditCommand makeSetResidueAtomsCommand(DocumentState& doc, int residueId,
                                                  const QList<AtomId>& oldAtoms,
                                                  const QList<AtomId>& newAtoms);
    
    // Internal wrapper function
    static void executeWrapper(BiopolymerEditCommand* cmd, DocumentState& doc);
    static void invertWrapper(BiopolymerEditCommand* cmd, DocumentState& doc);
};

// ========== Biopolymer State Mutator ==========

// High-level interface for creating biopolymer edit commands.
// This class provides methods that create properly-captured commands for common operations.
class BiopolymerStateMutator {
public:
    BiopolymerStateMutator(DocumentState& doc);
    ~BiopolymerStateMutator();
    
    // Chain operations
    
    // Add a new chain
    // Returns the command and the new chain ID
    std::pair<EditCommand, int> addChain(const QString& name, PolymerType type);
    
    // Remove a chain
    EditCommand removeChain(int chainId);
    
    // Rename a chain
    EditCommand renameChain(int chainId, const QString& newName);
    
    // Residue operations
    
    // Add a residue to a chain (at end)
    std::pair<EditCommand, int> addResidue(int chainId, const QString& symbol, 
                                         const QList<AtomId>& atomIds = {});
    
    // Insert a residue before another residue
    std::pair<EditCommand, int> insertResidueBefore(int chainId, int beforeResidueId,
                                                   const QString& symbol,
                                                   const QList<AtomId>& atomIds = {});
    
    // Insert a residue after another residue
    std::pair<EditCommand, int> insertResidueAfter(int chainId, int afterResidueId,
                                                  const QString& symbol,
                                                  const QList<AtomId>& atomIds = {});
    
    // Replace a residue
    EditCommand replaceResidue(int residueId, const QString& newSymbol);
    
    // Delete a residue
    EditCommand deleteResidue(int residueId);
    
    // Delete a range of residues
    EditCommand deleteResidueRange(int chainId, int startPosition, int endPosition);
    
    // Connection operations
    
    // Connect two residues
    EditCommand connectResidues(int fromResidueId, int toResidueId, 
                                 BondId bondId, const QString& type = "peptide");
    
    // Set residue atoms (for mapping)
    EditCommand setResidueAtoms(int residueId, const QList<AtomId>& atomIds);
    
    // Build a residue from a monomer template
    // Creates atoms and adds residue to biopolymer state
    std::pair<EditCommand, int> buildResidueFromTemplate(int chainId, const QString& symbol,
                                                         double x, double y);
    
private:
    std::shared_ptr<BiopolymerEditCommand> m_command;
    DocumentState& m_document;
    
    // Helper to capture current state for a residue
    struct ResidueState {
        int residueId;
        QString symbol;
        int chainId;
        int position;
        QList<AtomId> atomIds;
        QList<BondId> bondIds;
    };
    
    ResidueState captureResidueState(int residueId) const;
    
    // Helper to capture chain state
    struct ChainState {
        int chainId;
        QString name;
        PolymerType type;
        QList<int> residueIds;
        QList<AtomId> atomIds;
        QList<BondId> bondIds;
    };
    
    ChainState captureChainState(int chainId) const;
};

#endif // BIOPOLYMEREDITCOMMANDWRAPPER_H
