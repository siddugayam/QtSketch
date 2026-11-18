// src/app/molecule/BiopolymerEditCommandWrapper.cpp
#include "BiopolymerEditCommandWrapper.h"
#include "BiopolymerEditCommands.h"
#include "DocumentState.h"

BiopolymerEditCommandWrapper::BiopolymerEditCommandWrapper(
        BiopolymerEditCommand* command, DocumentState& document)
    : m_command(command), m_document(document) {
    Q_ASSERT(command != nullptr);
}

BiopolymerEditCommandWrapper::~BiopolymerEditCommandWrapper() {
    // Delegate ownership to caller via command()
    // Caller is responsible for cleaning up the underlying command
}

bool BiopolymerEditCommandWrapper::execute() {
    // Translate from document editing to biopolymer-specific
    return m_command->performEdit(m_document);
}

bool BiopolymerEditCommandWrapper::unexecute() {
    // Undo the edit by reversing the biopolymer-specific change
    return m_command->undoEdit(m_document);
}

bool BiopolymerEditCommandWrapper::canUndo() const {
    return m_command->canUndo(m_document);
}

QString BiopolymerEditCommandWrapper::description() const {
    return m_command->name();
}

bool BiopolymerEditCommandWrapper::isSignificant() const {
    return m_command->isSignificant();
}

void BiopolymerEditCommandWrapper::setSignificant(bool significant) {
    m_command->setSignificant(significant);
}

// -- factory for main EditCommand subclass
EditCommand* BiopolymerEditCommandWrapper::clone() const {
    return new BiopolymerEditCommandWrapper(m_command->clone(), m_document);
}

BiopolymerEditCommand* BiopolymerEditCommandWrapper::command() const {
    return m_command.get();
}

// -- Access to wrapped command
BiopolymerEditCommand* BiopolymerEditCommandWrapper::takeCommand() {
    return m_command.release();
}
