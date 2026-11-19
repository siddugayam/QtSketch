// src/app/molecule/BiopolymerTerminalCommands.cpp
#include "BiopolymerTerminalCommands.h"
#include "BiopolymerTerminal.h"

// ========== Terminal Edit Commands ==========

AddTerminalCommand::AddTerminalCommand(int chainId, TerminalType type, TerminalPosition position)
    : m_chainId(chainId), m_type(type), m_position(position) {}

AddTerminalCommand::~AddTerminalCommand() = default;

void AddTerminalCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    TerminalManager manager(bioState);
    m_terminalId = manager.addTerminal(m_chainId, m_type, m_position);
}

void AddTerminalCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_terminalId >= 0) {
        TerminalManager manager(bioState);
        manager.removeTerminal(m_terminalId);
    }
}

QString AddTerminalCommand::description() const {
    return QString("Add %1 terminal to chain %2")
        .arg(terminalTypeToString(m_type)).arg(m_chainId);
}

bool AddTerminalCommand::isValid() const {
    return m_chainId >= 0 && m_type != TerminalType::None;
}

BiopolymerEditCommand* AddTerminalCommand::clone() const {
    return new AddTerminalCommand(m_chainId, m_type, m_position);
}

// ==========

RemoveTerminalCommand::RemoveTerminalCommand(int terminalId)
    : m_terminalId(terminalId) {}

RemoveTerminalCommand::~RemoveTerminalCommand() = default;

void RemoveTerminalCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    TerminalManager manager(bioState);
    const BiopolymerTerminal* term = manager.getTerminal(m_terminalId);
    if (term) {
        m_savedTerminal = *term;
        m_hadTerminal = true;
    }
    manager.removeTerminal(m_terminalId);
}

void RemoveTerminalCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_hadTerminal) {
        TerminalManager manager(bioState);
        manager.addTerminal(m_savedTerminal);
    }
}

QString RemoveTerminalCommand::description() const {
    return QString("Remove terminal %1").arg(m_terminalId);
}

bool RemoveTerminalCommand::isValid() const {
    return m_terminalId >= 0;
}

BiopolymerEditCommand* RemoveTerminalCommand::clone() const {
    return new RemoveTerminalCommand(m_terminalId);
}

// ==========

SetTerminalCommand::SetTerminalCommand(int chainId, TerminalPosition position, TerminalType newType)
    : m_chainId(chainId), m_position(position), m_newType(newType) {}

SetTerminalCommand::~SetTerminalCommand() = default;

void SetTerminalCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    TerminalManager manager(bioState);
    
    // Get current terminal
    const BiopolymerTerminal* current = nullptr;
    if (m_position == TerminalPosition::First) {
        current = manager.nTerminal(m_chainId);
    } else {
        current = manager.cTerminal(m_chainId);
    }
    
    if (current) {
        m_oldType = current->type;
        m_hadTerminal = true;
    }
    
    // Set new terminal
    if (m_position == TerminalPosition::First) {
        manager.setNTerminal(m_chainId, m_newType);
    } else {
        manager.setCTerminal(m_chainId, m_newType);
    }
}

void SetTerminalCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    TerminalManager manager(bioState);
    
    if (m_hadTerminal) {
        if (m_position == TerminalPosition::First) {
            manager.setNTerminal(m_chainId, m_oldType);
        } else {
            manager.setCTerminal(m_chainId, m_oldType);
        }
    } else {
        // Remove the terminal
        if (m_position == TerminalPosition::First) {
            BiopolymerTerminal* term = manager.nTerminal(m_chainId);
            if (term) {
                manager.removeTerminal(term->terminalId);
            }
        } else {
            BiopolymerTerminal* term = manager.cTerminal(m_chainId);
            if (term) {
                manager.removeTerminal(term->terminalId);
            }
        }
    }
}

QString SetTerminalCommand::description() const {
    return QString("Set %1 terminal type for chain %2")
        .arg(terminalPositionToString(m_position)).arg(m_chainId);
}

bool SetTerminalCommand::isValid() const {
    return m_chainId >= 0 && m_newType != TerminalType::None;
}

BiopolymerEditCommand* SetTerminalCommand::clone() const {
    return new SetTerminalCommand(m_chainId, m_position, m_newType);
}

// ========== Modification Edit Commands ==========

AddModificationCommand::AddModificationCommand(int residueId, ModificationType type, const QString& name)
    : m_residueId(residueId), m_type(type), m_name(name) {}

AddModificationCommand::~AddModificationCommand() = default;

void AddModificationCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ModificationManager manager(bioState);
    BiopolymerModification mod;
    mod.residueId = m_residueId;
    mod.type = m_type;
    mod.name = m_name.isEmpty() ? modificationTypeToString(m_type) : m_name;
    mod.isSupported = true;
    
    m_modificationId = manager.addModification(mod);
}

void AddModificationCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_modificationId >= 0) {
        ModificationManager manager(bioState);
        manager.removeModification(m_modificationId);
    }
}

QString AddModificationCommand::description() const {
    return QString("Add %1 modification to residue %2")
        .arg(m_name.isEmpty() ? modificationTypeToString(m_type) : m_name).arg(m_residueId);
}

bool AddModificationCommand::isValid() const {
    return m_residueId >= 0 && m_type != ModificationType::None;
}

BiopolymerEditCommand* AddModificationCommand::clone() const {
    return new AddModificationCommand(m_residueId, m_type, m_name);
}

// ==========

RemoveModificationCommand::RemoveModificationCommand(int modificationId)
    : m_modificationId(modificationId) {}

RemoveModificationCommand::~RemoveModificationCommand() = default;

void RemoveModificationCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ModificationManager manager(bioState);
    const BiopolymerModification* mod = manager.getModification(m_modificationId);
    if (mod) {
        m_savedModification = *mod;
        m_hadModification = true;
    }
    manager.removeModification(m_modificationId);
}

void RemoveModificationCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_hadModification) {
        ModificationManager manager(bioState);
        manager.addModification(m_savedModification);
    }
}

QString RemoveModificationCommand::description() const {
    return QString("Remove modification %1").arg(m_modificationId);
}

bool RemoveModificationCommand::isValid() const {
    return m_modificationId >= 0;
}

BiopolymerEditCommand* RemoveModificationCommand::clone() const {
    return new RemoveModificationCommand(m_modificationId);
}

// ==========

SetModificationCommand::SetModificationCommand(int modificationId, ModificationType newType, const QString& newName)
    : m_modificationId(modificationId), m_newType(newType), m_newName(newName) {}

SetModificationCommand::~SetModificationCommand() = default;

void SetModificationCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ModificationManager manager(bioState);
    const BiopolymerModification* mod = manager.getModification(m_modificationId);
    
    if (mod) {
        m_oldType = mod->type;
        m_oldName = mod->name;
        m_hadModification = true;
        
        // Update the modification
        BiopolymerModification updated = *mod;
        updated.type = m_newType;
        if (!m_newName.isEmpty()) {
            updated.name = m_newName;
        }
        // Remove old and add updated
        manager.removeModification(m_modificationId);
        manager.addModification(updated);
    }
}

void SetModificationCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_hadModification) {
        ModificationManager manager(bioState);
        const BiopolymerModification* mod = manager.getModification(m_modificationId);
        if (mod) {
            BiopolymerModification reverted = *mod;
            reverted.type = m_oldType;
            reverted.name = m_oldName;
            manager.removeModification(m_modificationId);
            manager.addModification(reverted);
        }
    }
}

QString SetModificationCommand::description() const {
    return QString("Set modification %1 type to %2")
        .arg(m_modificationId).arg(modificationTypeToString(m_newType));
}

bool SetModificationCommand::isValid() const {
    return m_modificationId >= 0 && m_newType != ModificationType::None;
}

BiopolymerEditCommand* SetModificationCommand::clone() const {
    return new SetModificationCommand(m_modificationId, m_newType, m_newName);
}

// ========== Extended Connection Edit Commands ==========

AddDisulfideBondCommand::AddDisulfideBondCommand(int residueId1, int residueId2, BondId bondId,
                                               int atomId1, int atomId2)
    : m_residueId1(residueId1), m_residueId2(residueId2), m_bondId(bondId),
      m_atomId1(atomId1), m_atomId2(atomId2) {}

AddDisulfideBondCommand::~AddDisulfideBondCommand() = default;

void AddDisulfideBondCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ExtendedConnectionManager manager(bioState);
    m_disulfideId = manager.addDisulfide(m_residueId1, m_residueId2, m_bondId, m_atomId1, m_atomId2);
}

void AddDisulfideBondCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_disulfideId >= 0) {
        ExtendedConnectionManager manager(bioState);
        manager.removeDisulfide(m_disulfideId);
    }
}

QString AddDisulfideBondCommand::description() const {
    return QString("Add disulfide bond between residues %1 and %2")
        .arg(m_residueId1).arg(m_residueId2);
}

bool AddDisulfideBondCommand::isValid() const {
    return m_residueId1 >= 0 && m_residueId2 >= 0 && m_residueId1 != m_residueId2;
}

BiopolymerEditCommand* AddDisulfideBondCommand::clone() const {
    return new AddDisulfideBondCommand(m_residueId1, m_residueId2, m_bondId, m_atomId1, m_atomId2);
}

// ==========

RemoveDisulfideBondCommand::RemoveDisulfideBondCommand(int disulfideId)
    : m_disulfideId(disulfideId) {}

RemoveDisulfideBondCommand::~RemoveDisulfideBondCommand() = default;

void RemoveDisulfideBondCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ExtendedConnectionManager manager(bioState);
    const BiopolymerDisulfide* ds = manager.getDisulfide(m_disulfideId);
    if (ds) {
        m_savedDisulfide = *ds;
        m_hadDisulfide = true;
    }
    manager.removeDisulfide(m_disulfideId);
}

void RemoveDisulfideBondCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_hadDisulfide) {
        ExtendedConnectionManager manager(bioState);
        manager.addDisulfide(m_savedDisulfide);
    }
}

QString RemoveDisulfideBondCommand::description() const {
    return QString("Remove disulfide bond %1").arg(m_disulfideId);
}

bool RemoveDisulfideBondCommand::isValid() const {
    return m_disulfideId >= 0;
}

BiopolymerEditCommand* RemoveDisulfideBondCommand::clone() const {
    return new RemoveDisulfideBondCommand(m_disulfideId);
}

// ==========

AddInterchainConnectionCommand::AddInterchainConnectionCommand(
    int chainId1, int residueId1, int chainId2, int residueId2,
    BondId bondId, int atomId1, int atomId2)
    : m_chainId1(chainId1), m_residueId1(residueId1),
      m_chainId2(chainId2), m_residueId2(residueId2),
      m_bondId(bondId), m_atomId1(atomId1), m_atomId2(atomId2) {}

AddInterchainConnectionCommand::~AddInterchainConnectionCommand() = default;

void AddInterchainConnectionCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ExtendedConnectionManager manager(bioState);
    m_connectionId = manager.addInterchainConnection(
        m_chainId1, m_residueId1, m_chainId2, m_residueId2, m_bondId, m_atomId1, m_atomId2);
}

void AddInterchainConnectionCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_connectionId >= 0) {
        ExtendedConnectionManager manager(bioState);
        manager.removeInterchainConnection(m_connectionId);
    }
}

QString AddInterchainConnectionCommand::description() const {
    return QString("Add interchain connection between chains %1:%2 and %3:%4")
        .arg(m_chainId1).arg(m_residueId1).arg(m_chainId2).arg(m_residueId2);
}

bool AddInterchainConnectionCommand::isValid() const {
    return m_chainId1 >= 0 && m_residueId1 >= 0 && 
           m_chainId2 >= 0 && m_residueId2 >= 0 &&
           m_chainId1 != m_chainId2;
}

BiopolymerEditCommand* AddInterchainConnectionCommand::clone() const {
    return new AddInterchainConnectionCommand(
        m_chainId1, m_residueId1, m_chainId2, m_residueId2, m_bondId, m_atomId1, m_atomId2);
}

// ==========

RemoveInterchainConnectionBondCommand::RemoveInterchainConnectionBondCommand(int connectionId)
    : m_connectionId(connectionId) {}

RemoveInterchainConnectionBondCommand::~RemoveInterchainConnectionBondCommand() = default;

void RemoveInterchainConnectionBondCommand::execute(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    ExtendedConnectionManager manager(bioState);
    const BiopolymerInterchainConnection* ic = manager.getInterchainConnection(m_connectionId);
    if (ic) {
        m_savedConnection = *ic;
        m_hadConnection = true;
    }
    manager.removeInterchainConnection(m_connectionId);
}

void RemoveInterchainConnectionBondCommand::invert(EditableMolecule& molecule, BiopolymerState& bioState) {
    Q_UNUSED(molecule);
    if (m_hadConnection) {
        ExtendedConnectionManager manager(bioState);
        manager.addInterchainConnection(m_savedConnection);
    }
}

QString RemoveInterchainConnectionBondCommand::description() const {
    return QString("Remove interchain connection %1").arg(m_connectionId);
}

bool RemoveInterchainConnectionBondCommand::isValid() const {
    return m_connectionId >= 0;
}

BiopolymerEditCommand* RemoveInterchainConnectionBondCommand::clone() const {
    return new RemoveInterchainConnectionBondCommand(m_connectionId);
}
