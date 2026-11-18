// src/app/molecule/BiopolymerEditCommands.cpp
#include "BiopolymerEditCommands.h"
#include "BiopolymerState.h"
#include "EditableMolecule.h"

// ========== BiopolymerEditCommand ==========

BiopolymerEditCommand::~BiopolymerEditCommand() = default;

// ========== AddChainCommand ==========

void AddChainCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    // Add the chain to biopolymer state
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    bioState.addChain(chainName, type);
}

void AddChainCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    // Remove the chain from biopolymer state
    bioState.removeChain(chainId);
}

QString AddChainCommand::description() const {
    return QString("Add chain %1 (%2)").arg(chainName).arg(polymerTypeStr);
}

bool AddChainCommand::isValid() const {
    return chainId != -1 && !chainName.isEmpty();
}

// ========== RemoveChainCommand ==========

void RemoveChainCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove atoms and bonds from molecule
    for (BondId bondId : bondIds) {
        molecule.removeBond(bondId);
    }
    for (AtomId atomId : atomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Remove the chain from biopolymer state
    bioState.removeChain(chainId);
}

void RemoveChainCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Re-add the chain to biopolymer state first (to get the chainId back)
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    int newChainId = bioState.addChain(chainName, type);
    
    // The chainId should be the same (since we use monotonic IDs)
    if (newChainId != chainId) {
        // This shouldn't happen if IDs are managed correctly
        return;
    }
    
    // Re-add residues to the chain
    for (int residueId : residueIds) {
        // We need to re-create the residue with its atoms
        // For simplicity, we'll just restore the residue structure
        // In a full implementation, we'd need to restore the atom/bond mappings
    }
    
    // Re-add atoms and bonds
    // Note: This is simplified - in reality we'd need to restore coordinates too
    for (AtomId atomId : atomIds) {
        // Placeholder - actual atom restoration would need symbol and position
    }
    for (BondId bondId : bondIds) {
        // Placeholder - actual bond restoration would need endpoints and order
    }
}

QString RemoveChainCommand::description() const {
    return QString("Remove chain %1 (%2)").arg(chainName).arg(polymerTypeStr);
}

bool RemoveChainCommand::isValid() const {
    return chainId != -1 && !chainName.isEmpty();
}

// ========== RenameChainCommand ==========

void RenameChainCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    bioState.renameChain(chainId, newName);
}

void RenameChainCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    bioState.renameChain(chainId, oldName);
}

QString RenameChainCommand::description() const {
    return QString("Rename chain %1 to %2").arg(oldName).arg(newName);
}

bool RenameChainCommand::isValid() const {
    return chainId != -1 && !oldName.isEmpty() && !newName.isEmpty();
}

// ========== AddResidueCommand ==========

void AddResidueCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // In a full implementation, we would:
    // 1. Create the atoms in the molecule (if not already created)
    // 2. Create the bonds in the molecule (if not already created)
    // 3. Add the residue to the biopolymer state
    // 4. Set the residue atoms
    
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    int newResidueId = bioState.addResidue(chainId, {
        symbol, symbol, BiopolymerState::polymerTypeToString(type), false, "", {}
    }, position);
    
    if (newResidueId != residueId) {
        // Residue ID mismatch - this shouldn't happen
        return;
    }
    
    // Set the atoms for this residue
    if (!createdAtomIds.isEmpty()) {
        bioState.setResidueAtoms(residueId, createdAtomIds);
    }
}

void AddResidueCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove bonds first (to avoid cascading atom removal)
    for (BondId bondId : createdBondIds) {
        molecule.removeBond(bondId);
    }
    
    // Remove atoms
    for (AtomId atomId : createdAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Remove the residue from biopolymer state
    bioState.deleteResidue(residueId);
}

QString AddResidueCommand::description() const {
    return QString("Add residue %1 at position %2 in chain %3").arg(symbol).arg(position).arg(chainId);
}

bool AddResidueCommand::isValid() const {
    return residueId != -1 && chainId != -1 && !symbol.isEmpty();
}

// ========== InsertResidueCommand ==========

void InsertResidueCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    
    // Add the residue
    int newResidueId = bioState.insertResidueBefore(chainId, beforeResidueId, {
        symbol, symbol, BiopolymerState::polymerTypeToString(type), false, "", {}
    });
    
    if (newResidueId != residueId) {
        return;
    }
    
    // Set the atoms
    if (!createdAtomIds.isEmpty()) {
        bioState.setResidueAtoms(residueId, createdAtomIds);
    }
    
    // Create the bond to the next residue
    if (newBondToNext != -1) {
        // The bond should already exist in the molecule
        // We just need to set it in the biopolymer state
        const BiopolymerResidue* residue = bioState.residue(residueId);
        const BiopolymerResidue* nextRes = bioState.residue(beforeResidueId);
        if (residue && nextRes) {
            bioState.setResidueConnection(residueId, beforeResidueId, newBondToNext, "peptide");
        }
    }
}

void InsertResidueCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove bonds (in reverse order of creation)
    for (auto it = createdBondIds.rbegin(); it != createdBondIds.rend(); ++it) {
        molecule.removeBond(*it);
    }
    
    // Remove atoms
    for (AtomId atomId : createdAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Remove the residue
    bioState.deleteResidue(residueId);
    
    // Recreate the bond between neighbors
    if (oldBondBetweenNeighbors != -1) {
        // In a full implementation, we'd need to know the atom endpoints
        // For now, we'll just note that this would recreate the bond
    }
}

QString InsertResidueCommand::description() const {
    return QString("Insert residue %1 before %2 in chain %3").arg(symbol).arg(beforeResidueId).arg(chainId);
}

bool InsertResidueCommand::isValid() const {
    return residueId != -1 && chainId != -1 && beforeResidueId != -1 && !symbol.isEmpty();
}

// ========== ReplaceResidueCommand ==========

void ReplaceResidueCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove old atoms and bonds
    for (BondId bondId : oldBondIds) {
        molecule.removeBond(bondId);
    }
    for (AtomId atomId : oldAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Add new atoms and bonds
    // In a full implementation, we'd create them with proper positions
    
    // Replace the residue in biopolymer state
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    bioState.replaceResidue(residueId, {
        newSymbol, newSymbol, BiopolymerState::polymerTypeToString(type), false, "", {}
    });
    
    // Set the new atoms
    if (!newAtomIds.isEmpty()) {
        bioState.setResidueAtoms(residueId, newAtomIds);
    }
}

void ReplaceResidueCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove new atoms and bonds
    for (BondId bondId : newBondIds) {
        molecule.removeBond(bondId);
    }
    for (AtomId atomId : newAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Add old atoms and bonds back
    // In a full implementation, we'd restore them with original positions
    
    // Restore the old residue
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    bioState.replaceResidue(residueId, {
        oldSymbol, oldSymbol, BiopolymerState::polymerTypeToString(type), false, "", {}
    });
    
    // Set the old atoms
    if (!oldAtomIds.isEmpty()) {
        bioState.setResidueAtoms(residueId, oldAtomIds);
    }
}

QString ReplaceResidueCommand::description() const {
    return QString("Replace residue %1 with %2").arg(oldSymbol).arg(newSymbol);
}

bool ReplaceResidueCommand::isValid() const {
    return residueId != -1 && !oldSymbol.isEmpty() && !newSymbol.isEmpty();
}

// ========== DeleteResidueCommand ==========

void DeleteResidueCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove bonds first
    for (BondId bondId : deletedBondIds) {
        molecule.removeBond(bondId);
    }
    
    // Remove atoms
    for (AtomId atomId : deletedAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Remove the residue from biopolymer state
    bioState.deleteResidue(residueId);
    
    // If there's a bond to recreate between neighbors, create it
    if (bondToRecreate != -1 && prevResidueId != -1 && nextResidueId != -1) {
        // In a full implementation, we'd create the bond between the neighbors
        // bioState.setResidueConnection(prevResidueId, nextResidueId, bondToRecreate, "peptide");
    }
}

void DeleteResidueCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Re-add the residue
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    int newResidueId = bioState.insertResidueBefore(chainId, nextResidueId, {
        symbol, symbol, BiopolymerState::polymerTypeToString(type), false, "", {}
    });
    
    if (newResidueId != residueId) {
        // ID mismatch - can't properly undo
        return;
    }
    
    // Set the atoms
    if (!deletedAtomIds.isEmpty()) {
        bioState.setResidueAtoms(residueId, deletedAtomIds);
    }
    
    // Restore bonds
    // In a full implementation, we'd recreate the bonds with proper endpoints
}

QString DeleteResidueCommand::description() const {
    return QString("Delete residue %1 at position %2 in chain %3").arg(symbol).arg(position).arg(chainId);
}

bool DeleteResidueCommand::isValid() const {
    return residueId != -1 && chainId != -1;
}

// ========== DeleteResidueRangeCommand ==========

void DeleteResidueRangeCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove bonds first
    for (BondId bondId : deletedBondIds) {
        molecule.removeBond(bondId);
    }
    
    // Remove atoms
    for (AtomId atomId : deletedAtomIds) {
        molecule.removeAtom(atomId);
    }
    
    // Remove the residues from biopolymer state
    bioState.deleteResidueRange(chainId, startPosition, endPosition);
    
    // If there's a bond to recreate, create it
    if (bondToRecreate != -1) {
        // In a full implementation, we'd create the bond
    }
}

void DeleteResidueRangeCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    PolymerType type = BiopolymerState::polymerTypeFromString(polymerTypeStr);
    
    // Re-add residues in order
    for (int i = 0; i < deletedResidueIds.size(); ++i) {
        int residueId = deletedResidueIds[i];
        QString symbol = deletedSymbols[i];
        
        // Insert at the correct position
        int position = startPosition + i;
        bioState.insertResidueBefore(chainId, -1, {  // -1 means append
            symbol, symbol, BiopolymerState::polymerTypeToString(type), false, "", {}
        });
    }
    
    // Note: This is a simplified implementation
    // In a full implementation, we'd need to properly restore atoms, bonds, and positions
}

QString DeleteResidueRangeCommand::description() const {
    return QString("Delete residues %1-%2 in chain %3").arg(startPosition).arg(endPosition).arg(chainId);
}

bool DeleteResidueRangeCommand::isValid() const {
    return chainId != -1 && startPosition >= 1 && endPosition >= startPosition &&
           deletedResidueIds.size() == deletedSymbols.size();
}

// ========== ConnectResiduesCommand ==========

void ConnectResiduesCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    // The bond should already exist in the molecule
    // We just need to record it in the biopolymer state
    bioState.setResidueConnection(fromResidueId, toResidueId, bondId, connectionType);
}

void ConnectResiduesCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    // Remove the bond from the molecule
    if (bondId != -1) {
        molecule.removeBond(bondId);
    }
    
    // Remove the connection from biopolymer state
    // Note: BiopolymerState doesn't have a removeConnection method yet
    // For now, we'll just clear it by setting to -1
    bioState.setResidueConnection(fromResidueId, toResidueId, -1, "");
}

QString ConnectResiduesCommand::description() const {
    return QString("Connect residues %1 and %2").arg(fromResidueId).arg(toResidueId);
}

bool ConnectResiduesCommand::isValid() const {
    return fromResidueId != -1 && toResidueId != -1 && bondId != -1 &&
           fromAtomId != -1 && toAtomId != -1;
}

// ========== SetResidueAtomsCommand ==========

void SetResidueAtomsCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    bioState.setResidueAtoms(residueId, newAtomIds);
}

void SetResidueAtomsCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    bioState.setResidueAtoms(residueId, oldAtomIds);
}

QString SetResidueAtomsCommand::description() const {
    return QString("Set atoms for residue %1").arg(residueId);
}

bool SetResidueAtomsCommand::isValid() const {
    return residueId != -1;
}

// ========== Command Builder ==========

AddChainCommand BiopolymerEditCommandBuilder::buildAddChain(
    int chainId, const QString& name, const QString& type) {
    return AddChainCommand(chainId, name, type);
}

RemoveChainCommand BiopolymerEditCommandBuilder::buildRemoveChain(
    int chainId, const QString& name, const QString& type,
    const QList<int>& residues, const QList<AtomId>& atoms, const QList<BondId>& bonds) {
    return RemoveChainCommand(chainId, name, type, residues, atoms, bonds);
}

RenameChainCommand BiopolymerEditCommandBuilder::buildRenameChain(
    int chainId, const QString& oldName, const QString& newName) {
    return RenameChainCommand(chainId, oldName, newName);
}

AddResidueCommand BiopolymerEditCommandBuilder::buildAddResidue(
    int residueId, int chainId, int position,
    const QString& symbol, const QString& type,
    const QList<AtomId>& atoms, const QList<BondId>& bonds) {
    return AddResidueCommand(residueId, chainId, position, symbol, type, atoms, bonds);
}

InsertResidueCommand BiopolymerEditCommandBuilder::buildInsertResidue(
    int residueId, int chainId, int beforeResidueId,
    const QString& symbol, const QString& type,
    const QList<AtomId>& atoms, const QList<BondId>& bonds,
    BondId newBond, BondId oldBond) {
    return InsertResidueCommand(residueId, chainId, beforeResidueId, symbol, type,
                                  atoms, bonds, newBond, oldBond);
}

ReplaceResidueCommand BiopolymerEditCommandBuilder::buildReplaceResidue(
    int residueId, const QString& oldSymbol, const QString& newSymbol, const QString& type,
    const QList<AtomId>& oldAtoms, const QList<AtomId>& newAtoms,
    const QList<BondId>& oldBonds, const QList<BondId>& newBonds,
    const QMap<AtomId, AtomId>& atomMap) {
    return ReplaceResidueCommand(residueId, oldSymbol, newSymbol, type,
                                    oldAtoms, newAtoms, oldBonds, newBonds, atomMap);
}

DeleteResidueCommand BiopolymerEditCommandBuilder::buildDeleteResidue(
    int residueId, int chainId, int position,
    const QString& symbol, const QString& type,
    const QList<AtomId>& atoms, const QList<BondId>& bonds,
    BondId bondToRecreate, int prevRes, int nextRes) {
    return DeleteResidueCommand(residueId, chainId, position, symbol, type,
                                  atoms, bonds, bondToRecreate, prevRes, nextRes);
}

DeleteResidueRangeCommand BiopolymerEditCommandBuilder::buildDeleteResidueRange(
    int chainId, int start, int end,
    const QList<int>& residueIds, const QList<QString>& symbols,
    const QString& type,
    const QList<AtomId>& atoms, const QList<BondId>& bonds,
    BondId bondToRecreate) {
    return DeleteResidueRangeCommand(chainId, start, end, residueIds, symbols, type,
                                       atoms, bonds, bondToRecreate);
}

ConnectResiduesCommand BiopolymerEditCommandBuilder::buildConnectResidues(
    int fromRes, int toRes, BondId bond,
    AtomId fromAtom, AtomId toAtom, const QString& type) {
    return ConnectResiduesCommand(fromRes, toRes, bond, fromAtom, toAtom, type);
}

SetResidueAtomsCommand BiopolymerEditCommandBuilder::buildSetResidueAtoms(
    int residueId, const QList<AtomId>& oldAtoms, const QList<AtomId>& newAtoms) {
    return SetResidueAtomsCommand(residueId, oldAtoms, newAtoms);
}
