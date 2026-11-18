// src/app/molecule/BiopolymerEditCommands.h
#ifndef BIOPOLYMEREDITCOMMANDS_H
#define BIOPOLYMEREDITCOMMANDS_H

// Biopolymer-specific edit commands for transactional residue editing.
//
// Each command captures the complete state change for a biopolymer edit operation,
// including chemistry changes (atoms, bonds), biopolymer metadata (residues, chains, mappings),
// and validation state. Commands are created by BiopolymerStateMutator and executed
// through DocumentState's executeCommand/undo/redo mechanism.
//
// Design principle: each edit operation creates ONE undo entry, even if it involves
// multiple internal steps (atom creation, bond formation, residue mapping, etc.).
// This is achieved by capturing the complete before/after state in the command.

#include <functional>
#include <QList>
#include <QMap>
#include <QSet>
#include <QString>

using AtomId = int;
using BondId = int;

// Forward declarations
class EditableMolecule;
class BiopolymerState;

// Base class for biopolymer edit commands
class BiopolymerEditCommand {
public:
    virtual ~BiopolymerEditCommand();
    
    // Execute the command (apply changes)
    virtual void execute(EditableMolecule& molecule, BiopolymerState& bioState) = 0;
    
    // Invert the command (undo changes)
    virtual void invert(EditableMolecule& molecule, BiopolymerState& bioState) = 0;
    
    // Get a human-readable description
    virtual QString description() const = 0;
    
    // Check if the command is valid (can be executed)
    virtual bool isValid() const = 0;
    
    // Create a copy of this command (for wrapper use)
    virtual BiopolymerEditCommand* clone() const = 0;
};

// ========== Concrete Command Types ==========

// Add a new chain
class AddChainCommand : public BiopolymerEditCommand {
public:
    int chainId = -1;
    QString chainName;
    QString polymerTypeStr;
    
    AddChainCommand(int chainId, const QString& name, const QString& type)
        : chainId(chainId), chainName(name), polymerTypeStr(type) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new AddChainCommand(*this); }
};

// Remove a chain
class RemoveChainCommand : public BiopolymerEditCommand {
public:
    int chainId = -1;
    QString chainName;
    QString polymerTypeStr;
    QList<int> residueIds;
    QList<AtomId> atomIds;
    QList<BondId> bondIds;
    
    RemoveChainCommand(int chainId, const QString& name, const QString& type,
                       const QList<int>& residues, const QList<AtomId>& atoms,
                       const QList<BondId>& bonds)
        : chainId(chainId), chainName(name), polymerTypeStr(type),
          residueIds(residues), atomIds(atoms), bondIds(bonds) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new RemoveChainCommand(*this); }
};

// Rename a chain
class RenameChainCommand : public BiopolymerEditCommand {
public:
    int chainId = -1;
    QString oldName;
    QString newName;
    
    RenameChainCommand(int chainId, const QString& oldName, const QString& newName)
        : chainId(chainId), oldName(oldName), newName(newName) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new RenameChainCommand(*this); }
};

// Add a residue to a chain
class AddResidueCommand : public BiopolymerEditCommand {
public:
    int residueId = -1;
    int chainId = -1;
    int position = -1;
    QString symbol;
    QString polymerTypeStr;
    QList<AtomId> createdAtomIds;
    QList<BondId> createdBondIds;
    
    AddResidueCommand(int residueId, int chainId, int position,
                       const QString& symbol, const QString& type,
                       const QList<AtomId>& atoms, const QList<BondId>& bonds)
        : residueId(residueId), chainId(chainId), position(position),
          symbol(symbol), polymerTypeStr(type),
          createdAtomIds(atoms), createdBondIds(bonds) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new AddResidueCommand(*this); }
};

// Insert a residue before another residue
class InsertResidueCommand : public BiopolymerEditCommand {
public:
    int residueId = -1;
    int chainId = -1;
    int beforeResidueId = -1;
    QString symbol;
    QString polymerTypeStr;
    QList<AtomId> createdAtomIds;
    QList<BondId> createdBondIds;
    BondId newBondToNext = -1;
    BondId oldBondBetweenNeighbors = -1;
    
    InsertResidueCommand(int residueId, int chainId, int beforeResidueId,
                         const QString& symbol, const QString& type,
                         const QList<AtomId>& atoms, const QList<BondId>& bonds,
                         BondId newBond, BondId oldBond)
        : residueId(residueId), chainId(chainId), beforeResidueId(beforeResidueId),
          symbol(symbol), polymerTypeStr(type),
          createdAtomIds(atoms), createdBondIds(bonds),
          newBondToNext(newBond), oldBondBetweenNeighbors(oldBond) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new InsertResidueCommand(*this); }
};

// Replace a residue
class ReplaceResidueCommand : public BiopolymerEditCommand {
public:
    int residueId = -1;
    QString oldSymbol;
    QString newSymbol;
    QString polymerTypeStr;
    QList<AtomId> oldAtomIds;
    QList<AtomId> newAtomIds;
    QList<BondId> oldBondIds;
    QList<BondId> newBondIds;
    QMap<AtomId, AtomId> oldToNewAtomMap;
    
    ReplaceResidueCommand(int residueId, const QString& oldSymbol, const QString& newSymbol,
                         const QString& type,
                         const QList<AtomId>& oldAtoms, const QList<AtomId>& newAtoms,
                         const QList<BondId>& oldBonds, const QList<BondId>& newBonds,
                         const QMap<AtomId, AtomId>& atomMap)
        : residueId(residueId), oldSymbol(oldSymbol), newSymbol(newSymbol),
          polymerTypeStr(type),
          oldAtomIds(oldAtoms), newAtomIds(newAtoms),
          oldBondIds(oldBonds), newBondIds(newBonds),
          oldToNewAtomMap(atomMap) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new ReplaceResidueCommand(*this); }
};

// Delete a residue
class DeleteResidueCommand : public BiopolymerEditCommand {
public:
    int residueId = -1;
    int chainId = -1;
    int position = -1;
    QString symbol;
    QString polymerTypeStr;
    QList<AtomId> deletedAtomIds;
    QList<BondId> deletedBondIds;
    BondId bondToRecreate = -1;  // Bond between neighbors that needs to be recreated
    int prevResidueId = -1;
    int nextResidueId = -1;
    
    DeleteResidueCommand(int residueId, int chainId, int position,
                        const QString& symbol, const QString& type,
                        const QList<AtomId>& atoms, const QList<BondId>& bonds,
                        BondId bondToRecreate, int prevRes, int nextRes)
        : residueId(residueId), chainId(chainId), position(position),
          symbol(symbol), polymerTypeStr(type),
          deletedAtomIds(atoms), deletedBondIds(bonds),
          bondToRecreate(bondToRecreate),
          prevResidueId(prevRes), nextResidueId(nextRes) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new DeleteResidueCommand(*this); }
};

// Delete a range of residues
class DeleteResidueRangeCommand : public BiopolymerEditCommand {
public:
    int chainId = -1;
    int startPosition = -1;
    int endPosition = -1;
    QList<int> deletedResidueIds;
    QList<QString> deletedSymbols;
    QString polymerTypeStr;
    QList<AtomId> deletedAtomIds;
    QList<BondId> deletedBondIds;
    BondId bondToRecreate = -1;
    
    DeleteResidueRangeCommand(int chainId, int start, int end,
                             const QList<int>& residueIds, const QList<QString>& symbols,
                             const QString& type,
                             const QList<AtomId>& atoms, const QList<BondId>& bonds,
                             BondId bondToRecreate)
        : chainId(chainId), startPosition(start), endPosition(end),
          deletedResidueIds(residueIds), deletedSymbols(symbols),
          polymerTypeStr(type),
          deletedAtomIds(atoms), deletedBondIds(bonds),
          bondToRecreate(bondToRecreate) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new DeleteResidueRangeCommand(*this); }
};

// Connect two residues
class ConnectResiduesCommand : public BiopolymerEditCommand {
public:
    int fromResidueId = -1;
    int toResidueId = -1;
    BondId bondId = -1;
    AtomId fromAtomId = -1;
    AtomId toAtomId = -1;
    QString connectionType;
    
    ConnectResiduesCommand(int fromRes, int toRes, BondId bond,
                          AtomId fromAtom, AtomId toAtom, const QString& type)
        : fromResidueId(fromRes), toResidueId(toRes), bondId(bond),
          fromAtomId(fromAtom), toAtomId(toAtom), connectionType(type) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new ConnectResiduesCommand(*this); }
};

// Set residue atoms (for mapping)
class SetResidueAtomsCommand : public BiopolymerEditCommand {
public:
    int residueId = -1;
    QList<AtomId> oldAtomIds;
    QList<AtomId> newAtomIds;
    
    SetResidueAtomsCommand(int residueId, const QList<AtomId>& oldAtoms, const QList<AtomId>& newAtoms)
        : residueId(residueId), oldAtomIds(oldAtoms), newAtomIds(newAtoms) {}
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override { return new SetResidueAtomsCommand(*this); }
};

// ========== Command Builder ==========

// Helper class to build biopolymer edit commands with proper state capture
class BiopolymerEditCommandBuilder {
public:
    // Build an add chain command
    static AddChainCommand buildAddChain(int chainId, const QString& name, const QString& type);
    
    // Build a remove chain command
    static RemoveChainCommand buildRemoveChain(int chainId, const QString& name, const QString& type,
                                              const QList<int>& residues, const QList<AtomId>& atoms,
                                              const QList<BondId>& bonds);
    
    // Build a rename chain command
    static RenameChainCommand buildRenameChain(int chainId, const QString& oldName, const QString& newName);
    
    // Build an add residue command
    static AddResidueCommand buildAddResidue(int residueId, int chainId, int position,
                                            const QString& symbol, const QString& type,
                                            const QList<AtomId>& atoms, const QList<BondId>& bonds);
    
    // Build an insert residue command
    static InsertResidueCommand buildInsertResidue(int residueId, int chainId, int beforeResidueId,
                                                   const QString& symbol, const QString& type,
                                                   const QList<AtomId>& atoms, const QList<BondId>& bonds,
                                                   BondId newBond, BondId oldBond);
    
    // Build a replace residue command
    static ReplaceResidueCommand buildReplaceResidue(int residueId, const QString& oldSymbol,
                                                     const QString& newSymbol, const QString& type,
                                                     const QList<AtomId>& oldAtoms, const QList<AtomId>& newAtoms,
                                                     const QList<BondId>& oldBonds, const QList<BondId>& newBonds,
                                                     const QMap<AtomId, AtomId>& atomMap);
    
    // Build a delete residue command
    static DeleteResidueCommand buildDeleteResidue(int residueId, int chainId, int position,
                                                  const QString& symbol, const QString& type,
                                                  const QList<AtomId>& atoms, const QList<BondId>& bonds,
                                                  BondId bondToRecreate, int prevRes, int nextRes);
    
    // Build a delete residue range command
    static DeleteResidueRangeCommand buildDeleteResidueRange(int chainId, int start, int end,
                                                             const QList<int>& residueIds,
                                                             const QList<QString>& symbols,
                                                             const QString& type,
                                                             const QList<AtomId>& atoms,
                                                             const QList<BondId>& bonds,
                                                             BondId bondToRecreate);
    
    // Build a connect residues command
    static ConnectResiduesCommand buildConnectResidues(int fromRes, int toRes, BondId bond,
                                                       AtomId fromAtom, AtomId toAtom, const QString& type);
    
    // Build a set residue atoms command
    static SetResidueAtomsCommand buildSetResidueAtoms(int residueId, const QList<AtomId>& oldAtoms,
                                                      const QList<AtomId>& newAtoms);
};

#endif // BIOPOLYMEREDITCOMMANDS_H
