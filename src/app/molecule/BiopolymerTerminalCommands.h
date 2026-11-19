// src/app/molecule/BiopolymerTerminalCommands.h
#ifndef BIOPOLYMERTERMINALCOMMANDS_H
#define BIOPOLYMERTERMINALCOMMANDS_H

#include "BiopolymerTerminal.h"
#include "BiopolymerModification.h"
#include "BiopolymerConnectionExtended.h"
#include "BiopolymerEditCommands.h"

// Terminal edit commands that work with BiopolymerEditCommand base class

// Add terminal command
class AddTerminalCommand : public BiopolymerEditCommand {
public:
    AddTerminalCommand(int chainId, TerminalType type, TerminalPosition position);
    ~AddTerminalCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_chainId;
    TerminalType m_type;
    TerminalPosition m_position;
    int m_terminalId = -1;
};

// Remove terminal command
class RemoveTerminalCommand : public BiopolymerEditCommand {
public:
    explicit RemoveTerminalCommand(int terminalId);
    ~RemoveTerminalCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_terminalId;
    BiopolymerTerminal m_savedTerminal;
    bool m_hadTerminal = false;
};

// Set terminal command
class SetTerminalCommand : public BiopolymerEditCommand {
public:
    SetTerminalCommand(int chainId, TerminalPosition position, TerminalType newType);
    ~SetTerminalCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_chainId;
    TerminalPosition m_position;
    TerminalType m_newType;
    TerminalType m_oldType = TerminalType::None;
    bool m_hadTerminal = false;
};

// ========== Modification Edit Commands ==========

// Add modification command
class AddModificationCommand : public BiopolymerEditCommand {
public:
    AddModificationCommand(int residueId, ModificationType type, const QString& name = "");
    ~AddModificationCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_residueId;
    ModificationType m_type;
    QString m_name;
    int m_modificationId = -1;
};

// Remove modification command
class RemoveModificationCommand : public BiopolymerEditCommand {
public:
    explicit RemoveModificationCommand(int modificationId);
    ~RemoveModificationCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_modificationId;
    BiopolymerModification m_savedModification;
    bool m_hadModification = false;
};

// Set modification command
class SetModificationCommand : public BiopolymerEditCommand {
public:
    SetModificationCommand(int modificationId, ModificationType newType, const QString& newName = "");
    ~SetModificationCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_modificationId;
    ModificationType m_newType;
    QString m_newName;
    ModificationType m_oldType = ModificationType::None;
    QString m_oldName;
    bool m_hadModification = false;
};

// ========== Extended Connection Edit Commands ==========

// Add disulfide command
class AddDisulfideBondCommand : public BiopolymerEditCommand {
public:
    AddDisulfideBondCommand(int residueId1, int residueId2, BondId bondId, 
                           int atomId1 = -1, int atomId2 = -1);
    ~AddDisulfideBondCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_residueId1;
    int m_residueId2;
    BondId m_bondId;
    int m_atomId1;
    int m_atomId2;
    int m_disulfideId = -1;
};

// Remove disulfide command
class RemoveDisulfideBondCommand : public BiopolymerEditCommand {
public:
    explicit RemoveDisulfideBondCommand(int disulfideId);
    ~RemoveDisulfideBondCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_disulfideId;
    BiopolymerDisulfide m_savedDisulfide;
    bool m_hadDisulfide = false;
};

// Add interchain connection command
class AddInterchainConnectionCommand : public BiopolymerEditCommand {
public:
    AddInterchainConnectionCommand(int chainId1, int residueId1, int chainId2, int residueId2,
                                   BondId bondId, int atomId1 = -1, int atomId2 = -1);
    ~AddInterchainConnectionCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_chainId1;
    int m_residueId1;
    int m_chainId2;
    int m_residueId2;
    BondId m_bondId;
    int m_atomId1;
    int m_atomId2;
    int m_connectionId = -1;
};

// Remove interchain connection command
class RemoveInterchainConnectionBondCommand : public BiopolymerEditCommand {
public:
    explicit RemoveInterchainConnectionBondCommand(int connectionId);
    ~RemoveInterchainConnectionBondCommand() override;
    
    void execute(EditableMolecule& molecule, BiopolymerState& bioState) override;
    void invert(EditableMolecule& molecule, BiopolymerState& bioState) override;
    QString description() const override;
    bool isValid() const override;
    BiopolymerEditCommand* clone() const override;
    
private:
    int m_connectionId;
    BiopolymerInterchainConnection m_savedConnection;
    bool m_hadConnection = false;
};

#endif // BIOPOLYMERTERMINALCOMMANDS_H
