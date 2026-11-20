// src/app/molecule/BiopolymerEditCommandWrapper.cpp
#include "BiopolymerEditCommandWrapper.h"
#include "BiopolymerEditCommands.h"
#include "DocumentState.h"

// Memory management for shared_ptr
#include <memory>

// Static wrap method
EditCommand BiopolymerEditCommandWrapper::wrap(BiopolymerEditCommand* cmd, DocumentState& doc) {
    // Create a wrapper and use its execute/ undo lambdas
    // Use shared_ptr so the lambdas can share ownership
    auto wrapper = std::make_shared<BiopolymerEditCommandWrapper>(cmd, doc);
    EditCommand editCmd;
    editCmd.execute = [wrapper]() {
        wrapper->execute();
    };
    editCmd.invert = [wrapper]() {
        wrapper->unexecute();
    };
    return editCmd;
}

BiopolymerEditCommandWrapper::BiopolymerEditCommandWrapper(
        BiopolymerEditCommand* command, DocumentState& document)
    : m_command(command), m_document(document) {
    Q_ASSERT(command != nullptr);
}

BiopolymerEditCommandWrapper::~BiopolymerEditCommandWrapper() = default;

void BiopolymerEditCommandWrapper::execute() {
    // Translate from document editing to biopolymer-specific
    m_command->execute(m_document.molecule(), m_document.biopolymerState());
}

void BiopolymerEditCommandWrapper::unexecute() {
    // Undo the edit by reversing the biopolymer-specific change
    m_command->invert(m_document.molecule(), m_document.biopolymerState());
}

bool BiopolymerEditCommandWrapper::canUndo() const {
    return m_command->isValid();
}

QString BiopolymerEditCommandWrapper::description() const {
    return m_command->description();
}

bool BiopolymerEditCommandWrapper::isSignificant() const {
    return true; // Biopolymer commands are always significant
}

void BiopolymerEditCommandWrapper::setSignificant(bool significant) {
    // Not supported by BiopolymerEditCommand, but we need this for the interface
    Q_UNUSED(significant);
}

// -- Access to wrapped command
BiopolymerEditCommand* BiopolymerEditCommandWrapper::command() const {
    return m_command.get();
}

// -- Access to wrapped command
BiopolymerEditCommand* BiopolymerEditCommandWrapper::takeCommand() {
    return m_command.get();
}

// ========== BiopolymerStateMutator ==========

BiopolymerStateMutator::BiopolymerStateMutator(DocumentState& doc)
    : m_document(doc) {
}

BiopolymerStateMutator::~BiopolymerStateMutator() = default;

std::pair<EditCommand, int> BiopolymerStateMutator::addChain(const QString& name, PolymerType type) {
    Q_UNUSED(name); Q_UNUSED(type);
    return {EditCommand{}, -1};
}

EditCommand BiopolymerStateMutator::removeChain(int chainId) {
    Q_UNUSED(chainId);
    return EditCommand{};
}

EditCommand BiopolymerStateMutator::renameChain(int chainId, const QString& newName) {
    Q_UNUSED(chainId); Q_UNUSED(newName);
    return EditCommand{};
}

std::pair<EditCommand, int> BiopolymerStateMutator::addResidue(int chainId, const QString& symbol,
                                                               const QList<AtomId>& atomIds) {
    Q_UNUSED(chainId); Q_UNUSED(symbol); Q_UNUSED(atomIds);
    return {EditCommand{}, -1};
}

std::pair<EditCommand, int> BiopolymerStateMutator::insertResidueBefore(int chainId, int beforeResidueId,
                                                                         const QString& symbol,
                                                                         const QList<AtomId>& atomIds) {
    Q_UNUSED(chainId); Q_UNUSED(beforeResidueId); Q_UNUSED(symbol); Q_UNUSED(atomIds);
    return {EditCommand{}, -1};
}

std::pair<EditCommand, int> BiopolymerStateMutator::insertResidueAfter(int chainId, int afterResidueId,
                                                                        const QString& symbol,
                                                                        const QList<AtomId>& atomIds) {
    Q_UNUSED(chainId); Q_UNUSED(afterResidueId); Q_UNUSED(symbol); Q_UNUSED(atomIds);
    return {EditCommand{}, -1};
}

EditCommand BiopolymerStateMutator::replaceResidue(int residueId, const QString& newSymbol) {
    Q_UNUSED(residueId); Q_UNUSED(newSymbol);
    return EditCommand{};
}

EditCommand BiopolymerStateMutator::deleteResidue(int residueId) {
    Q_UNUSED(residueId);
    return EditCommand{};
}

EditCommand BiopolymerStateMutator::deleteResidueRange(int chainId, int startPosition, int endPosition) {
    Q_UNUSED(chainId); Q_UNUSED(startPosition); Q_UNUSED(endPosition);
    return EditCommand{};
}

EditCommand BiopolymerStateMutator::connectResidues(int fromResidueId, int toResidueId,
                                                    BondId bondId, const QString& type) {
    Q_UNUSED(fromResidueId); Q_UNUSED(toResidueId); Q_UNUSED(bondId); Q_UNUSED(type);
    return EditCommand{};
}

EditCommand BiopolymerStateMutator::setResidueAtoms(int residueId, const QList<AtomId>& atomIds) {
    Q_UNUSED(residueId); Q_UNUSED(atomIds);
    return EditCommand{};
}

std::pair<EditCommand, int> BiopolymerStateMutator::buildResidueFromTemplate(int chainId, const QString& symbol,
                                                                             double x, double y) {
    Q_UNUSED(chainId); Q_UNUSED(symbol); Q_UNUSED(x); Q_UNUSED(y);
    return {EditCommand{}, -1};
}

BiopolymerStateMutator::ResidueState BiopolymerStateMutator::captureResidueState(int residueId) const {
    Q_UNUSED(residueId);
    return ResidueState{};
}

BiopolymerStateMutator::ChainState BiopolymerStateMutator::captureChainState(int chainId) const {
    Q_UNUSED(chainId);
    return ChainState{};
}
