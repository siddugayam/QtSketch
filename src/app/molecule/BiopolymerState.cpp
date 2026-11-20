// src/app/molecule/BiopolymerState.cpp
#include "BiopolymerState.h"

BiopolymerState::BiopolymerState() = default;

BiopolymerState::BiopolymerState(const BiopolymerState& other)
    : m_nextResidueId(other.m_nextResidueId),
      m_nextChainId(other.m_nextChainId),
      m_nextConnectionId(other.m_nextConnectionId),
      m_nextTerminalId(other.m_nextTerminalId),
      m_nextModificationId(other.m_nextModificationId),
      m_chains(other.m_chains),
      m_residues(other.m_residues),
      m_connections(other.m_connections),
      m_terminals(other.m_terminals),
      m_modifications(other.m_modifications),
      m_atomToResidue(other.m_atomToResidue),
      m_importProvenance(other.m_importProvenance)
{
    // Validation cache is intentionally not copied (will be rebuilt on demand)
}

BiopolymerState& BiopolymerState::operator=(const BiopolymerState& other) {
    if (this != &other) {
        m_nextResidueId = other.m_nextResidueId;
        m_nextChainId = other.m_nextChainId;
        m_nextConnectionId = other.m_nextConnectionId;
        m_nextTerminalId = other.m_nextTerminalId;
        m_nextModificationId = other.m_nextModificationId;
        m_chains = other.m_chains;
        m_residues = other.m_residues;
        m_connections = other.m_connections;
        m_terminals = other.m_terminals;
        m_modifications = other.m_modifications;
        m_atomToResidue = other.m_atomToResidue;
        m_importProvenance = other.m_importProvenance;
        // Validation cache is intentionally not copied (will be rebuilt on demand)
        invalidateValidationCache();
    }
    return *this;
}

BiopolymerState::~BiopolymerState() = default;

// ========== Chain Management ==========

int BiopolymerState::addChain(const QString& name, PolymerType type) {
    int chainId = allocateChainId();
    BiopolymerChain chain;
    chain.id = chainId;
    chain.name = name;
    chain.polymerType = type;
    m_chains[chainId] = chain;
    invalidateValidationCache();
    return chainId;
}

bool BiopolymerState::removeChain(int chainId) {
    auto it = m_chains.find(chainId);
    if (it == m_chains.end()) {
        return false;
    }
    
    // Remove all residues in this chain
    for (int residueId : it->residueIds) {
        m_residues.remove(residueId);
    }
    
    // Remove connections involving this chain's residues
    QList<BiopolymerConnection> keptConnections;
    for (const auto& conn : m_connections) {
        if (it->residueIds.contains(conn.fromResidueId) || 
            it->residueIds.contains(conn.toResidueId)) {
            continue;
        }
        keptConnections.append(conn);
    }
    m_connections = keptConnections;
    
    // Remove from atom-to-residue map
    rebuildAtomToResidueMap();
    
    m_chains.erase(it);
    invalidateValidationCache();
    return true;
}

bool BiopolymerState::renameChain(int chainId, const QString& newName) {
    auto it = m_chains.find(chainId);
    if (it == m_chains.end()) {
        return false;
    }
    it->name = newName;
    return true;
}

BiopolymerChain* BiopolymerState::chain(int chainId) {
    return findChain(chainId);
}

const BiopolymerChain* BiopolymerState::chain(int chainId) const {
    return findChain(chainId);
}

QList<int> BiopolymerState::chainIds() const {
    return m_chains.keys();
}

int BiopolymerState::chainCount() const {
    return m_chains.size();
}

// ========== Residue Management ==========

int BiopolymerState::addResidue(int chainId, const MonomerIdentity& identity, int position) {
    auto chain = findChain(chainId);
    if (!chain) {
        return -1;
    }
    
    int residueId = allocateResidueId();
    BiopolymerResidue residue;
    residue.id = residueId;
    residue.chainId = chainId;
    residue.identity = identity;
    
    if (position == -1) {
        // Append at end
        residue.sequencePosition = chain->residueIds.size() + 1;
        chain->residueIds.append(residueId);
    } else {
        // Insert at specific position (1-indexed)
        residue.sequencePosition = position;
        // Adjust positions of existing residues
        for (int i = position - 1; i < chain->residueIds.size(); ++i) {
            int existingId = chain->residueIds[i];
            m_residues[existingId].sequencePosition++;
        }
        chain->residueIds.insert(position - 1, residueId);
    }
    
    m_residues[residueId] = residue;
    invalidateValidationCache();
    return residueId;
}

int BiopolymerState::insertResidueBefore(int chainId, int beforeResidueId, const MonomerIdentity& identity) {
    auto chain = findChain(chainId);
    if (!chain) {
        return -1;
    }
    
    auto residueIt = m_residues.find(beforeResidueId);
    if (residueIt == m_residues.end() || residueIt->chainId != chainId) {
        return -1;
    }
    
    int position = residueIt->sequencePosition;
    return addResidue(chainId, identity, position);
}

int BiopolymerState::insertResidueAfter(int chainId, int afterResidueId, const MonomerIdentity& identity) {
    auto chain = findChain(chainId);
    if (!chain) {
        return -1;
    }
    
    auto residueIt = m_residues.find(afterResidueId);
    if (residueIt == m_residues.end() || residueIt->chainId != chainId) {
        return -1;
    }
    
    int position = residueIt->sequencePosition + 1;
    return addResidue(chainId, identity, position);
}

bool BiopolymerState::replaceResidue(int residueId, const MonomerIdentity& identity) {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return false;
    }
    it->identity = identity;
    invalidateValidationCache();
    return true;
}

bool BiopolymerState::deleteResidue(int residueId) {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return false;
    }
    
    int chainId = it->chainId;
    auto chain = findChain(chainId);
    if (!chain) {
        return false;
    }
    
    // Remove from chain
    int position = it->sequencePosition;
    chain->residueIds.removeAll(residueId);
    
    // Adjust positions of subsequent residues
    for (int i = position; i < chain->residueIds.size(); ++i) {
        int existingId = chain->residueIds[i];
        m_residues[existingId].sequencePosition--;
    }
    
    // Remove residue
    m_residues.erase(it);
    
    // Remove connections involving this residue
    QList<BiopolymerConnection> keptConnections;
    for (const auto& conn : m_connections) {
        if (conn.fromResidueId == residueId || conn.toResidueId == residueId) {
            continue;
        }
        keptConnections.append(conn);
    }
    m_connections = keptConnections;
    
    // Rebuild atom-to-residue map
    rebuildAtomToResidueMap();
    
    invalidateValidationCache();
    return true;
}

bool BiopolymerState::deleteResidueRange(int chainId, int startPosition, int endPosition) {
    auto chain = findChain(chainId);
    if (!chain) {
        return false;
    }
    
    if (startPosition < 1 || endPosition > chain->residueIds.size() || startPosition > endPosition) {
        return false;
    }
    
    // Collect residue IDs to delete
    QList<int> residuesToDelete;
    for (int pos = startPosition; pos <= endPosition; ++pos) {
        int residueId = chain->residueIds[pos - 1];
        residuesToDelete.append(residueId);
    }
    
    // Delete each residue
    for (int residueId : residuesToDelete) {
        deleteResidue(residueId);
    }
    
    invalidateValidationCache();
    return true;
}

BiopolymerResidue* BiopolymerState::residue(int residueId) {
    auto it = m_residues.find(residueId);
    return (it != m_residues.end()) ? &it.value() : nullptr;
}

const BiopolymerResidue* BiopolymerState::residue(int residueId) const {
    auto it = m_residues.find(residueId);
    return (it != m_residues.end()) ? &it.value() : nullptr;
}

BiopolymerResidue* BiopolymerState::residueInChain(int chainId, int position) {
    auto chain = findChain(chainId);
    if (!chain || position < 1 || position > chain->residueIds.size()) {
        return nullptr;
    }
    return residue(chain->residueIds[position - 1]);
}

const BiopolymerResidue* BiopolymerState::residueInChain(int chainId, int position) const {
    auto chain = findChain(chainId);
    if (!chain || position < 1 || position > chain->residueIds.size()) {
        return nullptr;
    }
    return residue(chain->residueIds[position - 1]);
}

QList<int> BiopolymerState::residueIdsInChain(int chainId) const {
    auto chain = findChain(chainId);
    return chain ? chain->residueIds : QList<int>();
}

QList<int> BiopolymerState::allResidueIds() const {
    return m_residues.keys();
}

// ========== Atom Mapping ==========

void BiopolymerState::setResidueAtoms(int residueId, const QList<AtomId>& atomIds) {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return;
    }
    it->ownedAtomIds = atomIds;
    rebuildAtomToResidueMap();
    invalidateValidationCache();
}

void BiopolymerState::addResidueAtoms(int residueId, const QList<AtomId>& atomIds) {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return;
    }
    for (AtomId atomId : atomIds) {
        if (!it->ownedAtomIds.contains(atomId)) {
            it->ownedAtomIds.append(atomId);
        }
    }
    rebuildAtomToResidueMap();
    invalidateValidationCache();
}

void BiopolymerState::removeResidueAtoms(int residueId, const QList<AtomId>& atomIds) {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return;
    }
    for (AtomId atomId : atomIds) {
        it->ownedAtomIds.removeAll(atomId);
    }
    rebuildAtomToResidueMap();
    invalidateValidationCache();
}

QList<AtomId> BiopolymerState::atomsOwnedByResidue(int residueId) const {
    auto it = m_residues.find(residueId);
    return (it != m_residues.end()) ? it->ownedAtomIds : QList<AtomId>();
}

int BiopolymerState::residueOwningAtom(AtomId atomId) const {
    auto it = m_atomToResidue.find(atomId);
    return (it != m_atomToResidue.end()) ? it.value() : -1;
}

void BiopolymerState::setResidueConnection(int fromResidueId, int toResidueId, BondId bondId, const QString& type) {
    // Check if connection already exists
    for (auto& conn : m_connections) {
        if (conn.fromResidueId == fromResidueId && conn.toResidueId == toResidueId) {
            conn.bondId = bondId;
            conn.type = type;
            invalidateValidationCache();
            return;
        }
    }
    
    // Create new connection
    BiopolymerConnection conn;
    conn.id = allocateConnectionId();
    conn.fromResidueId = fromResidueId;
    conn.toResidueId = toResidueId;
    conn.bondId = bondId;
    conn.type = type;
    m_connections.append(conn);
    invalidateValidationCache();
}

BondId BiopolymerState::connectionBond(int fromResidueId, int toResidueId) const {
    for (const auto& conn : m_connections) {
        if (conn.fromResidueId == fromResidueId && conn.toResidueId == toResidueId) {
            return conn.bondId;
        }
    }
    return -1;
}

QList<BiopolymerConnection> BiopolymerState::allConnections() const {
    return m_connections;
}

// ========== Terminal Groups ==========

void BiopolymerState::setNTerminus(int chainId, const TerminalGroup& terminus) {
    auto chain = findChain(chainId);
    if (chain) {
        chain->nTerminus = terminus;
        invalidateValidationCache();
    }
}

void BiopolymerState::setCTerminus(int chainId, const TerminalGroup& terminus) {
    auto chain = findChain(chainId);
    if (chain) {
        chain->cTerminus = terminus;
        invalidateValidationCache();
    }
}

void BiopolymerState::setPhosphate5Prime(int chainId, const TerminalGroup& phosphate) {
    auto chain = findChain(chainId);
    if (chain) {
        chain->phosphate5Prime = phosphate;
        invalidateValidationCache();
    }
}

void BiopolymerState::setPhosphate3Prime(int chainId, const TerminalGroup& phosphate) {
    auto chain = findChain(chainId);
    if (chain) {
        chain->phosphate3Prime = phosphate;
        invalidateValidationCache();
    }
}

TerminalGroup BiopolymerState::getNTerminus(int chainId) const {
    auto chain = findChain(chainId);
    return chain ? chain->nTerminus : TerminalGroup{};
}

TerminalGroup BiopolymerState::getCTerminus(int chainId) const {
    auto chain = findChain(chainId);
    return chain ? chain->cTerminus : TerminalGroup{};
}

TerminalGroup BiopolymerState::getPhosphate5Prime(int chainId) const {
    auto chain = findChain(chainId);
    return chain ? chain->phosphate5Prime : TerminalGroup{};
}

TerminalGroup BiopolymerState::getPhosphate3Prime(int chainId) const {
    auto chain = findChain(chainId);
    return chain ? chain->phosphate3Prime : TerminalGroup{};
}

// ========== Biopolymer Terminals (Extended) ==========

bool BiopolymerState::addTerminal(const BiopolymerTerminal& terminal) {
    // Check if terminal ID is valid
    if (terminal.terminalId < 0) {
        return false;
    }
    
    // Check if chain exists
    if (!m_chains.contains(terminal.chainId)) {
        return false;
    }
    
    // Check if terminal with same ID already exists
    if (m_terminals.contains(terminal.terminalId)) {
        return false;
    }
    
    // Add the terminal
    m_terminals[terminal.terminalId] = terminal;
    
    // Update next terminal ID if needed
    m_nextTerminalId = std::max(m_nextTerminalId, terminal.terminalId + 1);
    
    invalidateValidationCache();
    return true;
}

bool BiopolymerState::removeTerminal(int terminalId) {
    if (!m_terminals.contains(terminalId)) {
        return false;
    }
    
    m_terminals.remove(terminalId);
    invalidateValidationCache();
    return true;
}

void BiopolymerState::removeTerminalsFromChain(int chainId) {
    QList<int> toRemove;
    for (auto it = m_terminals.begin(); it != m_terminals.end(); ++it) {
        if (it->chainId == chainId) {
            toRemove.append(it.key());
        }
    }
    
    for (int termId : toRemove) {
        m_terminals.remove(termId);
    }
    
    invalidateValidationCache();
}

void BiopolymerState::removeTerminalsByType(TerminalType type) {
    QList<int> toRemove;
    for (auto it = m_terminals.begin(); it != m_terminals.end(); ++it) {
        if (it->type == type) {
            toRemove.append(it.key());
        }
    }
    
    for (int termId : toRemove) {
        m_terminals.remove(termId);
    }
    
    invalidateValidationCache();
}

BiopolymerTerminal* BiopolymerState::getTerminal(int terminalId) {
    auto it = m_terminals.find(terminalId);
    return (it != m_terminals.end()) ? &it.value() : nullptr;
}

const BiopolymerTerminal* BiopolymerState::getTerminal(int terminalId) const {
    auto it = m_terminals.find(terminalId);
    return (it != m_terminals.end()) ? &it.value() : nullptr;
}

QList<BiopolymerTerminal> BiopolymerState::terminals() const {
    return m_terminals.values();
}

QList<BiopolymerTerminal> BiopolymerState::terminalsForChain(int chainId) const {
    QList<BiopolymerTerminal> result;
    for (const auto& term : m_terminals) {
        if (term.chainId == chainId) {
            result.append(term);
        }
    }
    return result;
}

QList<BiopolymerTerminal> BiopolymerState::terminalsByType(TerminalType type) const {
    QList<BiopolymerTerminal> result;
    for (const auto& term : m_terminals) {
        if (term.type == type) {
            result.append(term);
        }
    }
    return result;
}

int BiopolymerState::allocateTerminalId() {
    return m_nextTerminalId++;
}

// ========== Biopolymer Modifications ==========

bool BiopolymerState::addModification(const BiopolymerModification& modification) {
    // Check if modification ID is valid
    if (modification.modificationId < 0) {
        return false;
    }
    
    // Check if residue exists
    if (!m_residues.contains(modification.residueId)) {
        return false;
    }
    
    // Check if modification with same ID already exists
    if (m_modifications.contains(modification.modificationId)) {
        return false;
    }
    
    // Add the modification
    m_modifications[modification.modificationId] = modification;
    
    // Update next modification ID if needed
    m_nextModificationId = std::max(m_nextModificationId, modification.modificationId + 1);
    
    invalidateValidationCache();
    return true;
}

bool BiopolymerState::removeModification(int modificationId) {
    if (!m_modifications.contains(modificationId)) {
        return false;
    }
    
    m_modifications.remove(modificationId);
    invalidateValidationCache();
    return true;
}

void BiopolymerState::removeModificationsFromResidue(int residueId) {
    QList<int> toRemove;
    for (auto it = m_modifications.begin(); it != m_modifications.end(); ++it) {
        if (it->residueId == residueId) {
            toRemove.append(it.key());
        }
    }
    
    for (int modId : toRemove) {
        m_modifications.remove(modId);
    }
    
    invalidateValidationCache();
}

void BiopolymerState::removeModificationsFromChain(int chainId) {
    QList<int> toRemove;
    for (auto it = m_modifications.begin(); it != m_modifications.end(); ++it) {
        const BiopolymerResidue* residue = this->residue(it->residueId);
        if (residue && residue->chainId == chainId) {
            toRemove.append(it.key());
        }
    }
    
    for (int modId : toRemove) {
        m_modifications.remove(modId);
    }
    
    invalidateValidationCache();
}

void BiopolymerState::removeModificationsByType(ModificationType type) {
    QList<int> toRemove;
    for (auto it = m_modifications.begin(); it != m_modifications.end(); ++it) {
        if (it->type == type) {
            toRemove.append(it.key());
        }
    }
    
    for (int modId : toRemove) {
        m_modifications.remove(modId);
    }
    
    invalidateValidationCache();
}

BiopolymerModification* BiopolymerState::getModification(int modificationId) {
    auto it = m_modifications.find(modificationId);
    return (it != m_modifications.end()) ? &it.value() : nullptr;
}

const BiopolymerModification* BiopolymerState::getModification(int modificationId) const {
    auto it = m_modifications.find(modificationId);
    return (it != m_modifications.end()) ? &it.value() : nullptr;
}

QList<BiopolymerModification> BiopolymerState::modifications() const {
    return m_modifications.values();
}

QList<BiopolymerModification> BiopolymerState::modificationsForResidue(int residueId) const {
    QList<BiopolymerModification> result;
    for (const auto& mod : m_modifications) {
        if (mod.residueId == residueId) {
            result.append(mod);
        }
    }
    return result;
}

QList<BiopolymerModification> BiopolymerState::modificationsForChain(int chainId) const {
    QList<BiopolymerModification> result;
    for (const auto& mod : m_modifications) {
        const BiopolymerResidue* residue = this->residue(mod.residueId);
        if (residue && residue->chainId == chainId) {
            result.append(mod);
        }
    }
    return result;
}

QList<BiopolymerModification> BiopolymerState::modificationsByType(ModificationType type) const {
    QList<BiopolymerModification> result;
    for (const auto& mod : m_modifications) {
        if (mod.type == type) {
            result.append(mod);
        }
    }
    return result;
}

bool BiopolymerState::hasModifications(int residueId) const {
    for (const auto& mod : m_modifications) {
        if (mod.residueId == residueId) {
            return true;
        }
    }
    return false;
}

bool BiopolymerState::hasChainModifications(int chainId) const {
    for (const auto& mod : m_modifications) {
        const BiopolymerResidue* residue = this->residue(mod.residueId);
        if (residue && residue->chainId == chainId) {
            return true;
        }
    }
    return false;
}

int BiopolymerState::allocateModificationId() {
    return m_nextModificationId++;
}

BiopolymerState::AllIdCounters BiopolymerState::allIdCounters() const {
    return {
        m_nextResidueId,
        m_nextChainId,
        m_nextConnectionId,
        m_nextTerminalId,
        m_nextModificationId
    };
}

// ========== Validation ==========

MappingValidity BiopolymerState::validateAllMappings() const {
    MappingValidity validity;
    
    for (int chainId : m_chains.keys()) {
        MappingValidity chainValidity = validateChain(chainId);
        if (!chainValidity.isValid) {
            validity.isValid = false;
            validity.issues.append(chainValidity.issues);
        }
    }
    
    return validity;
}

MappingValidity BiopolymerState::validateResidueMapping(int residueId) const {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return MappingValidity{false, {"Residue not found"}};
    }
    
    const auto& residue = it.value();
    MappingValidity validity;
    
    // Check if residue has atoms
    if (residue.ownedAtomIds.isEmpty()) {
        validity.isValid = false;
        validity.issues.append(QString("Residue %1 has no mapped atoms").arg(residueId));
    }
    
    // Check if residue is in a valid chain
    if (!m_chains.contains(residue.chainId)) {
        validity.isValid = false;
        validity.issues.append(QString("Residue %1 references non-existent chain %2").arg(residueId).arg(residue.chainId));
    }
    
    // Check for duplicate atoms
    QSet<AtomId> uniqueAtoms(residue.ownedAtomIds.begin(), residue.ownedAtomIds.end());
    if (uniqueAtoms.size() != residue.ownedAtomIds.size()) {
        validity.isValid = false;
        validity.issues.append(QString("Residue %1 has duplicate atom IDs").arg(residueId));
    }
    
    return validity;
}

MappingValidity BiopolymerState::validateChain(int chainId) const {
    auto chainIt = m_chains.find(chainId);
    if (chainIt == m_chains.end()) {
        return MappingValidity{false, {"Chain not found"}};
    }
    
    const auto& chain = chainIt.value();
    MappingValidity validity;
    
    // Check each residue in the chain
    for (int residueId : chain.residueIds) {
        MappingValidity residueValidity = validateResidueMapping(residueId);
        if (!residueValidity.isValid) {
            validity.isValid = false;
            validity.issues.append(residueValidity.issues);
        }
    }
    
    // Check connections between residues
    for (int i = 0; i < chain.residueIds.size() - 1; ++i) {
        int fromResidueId = chain.residueIds[i];
        int toResidueId = chain.residueIds[i + 1];
        
        bool hasConnection = false;
        for (const auto& conn : m_connections) {
            if (conn.fromResidueId == fromResidueId && conn.toResidueId == toResidueId) {
                hasConnection = true;
                if (conn.bondId == -1) {
                    validity.isValid = false;
                    validity.issues.append(QString("Missing bond between residues %1 and %2").arg(fromResidueId).arg(toResidueId));
                }
                break;
            }
        }
        
        if (!hasConnection) {
            validity.isValid = false;
            validity.issues.append(QString("No connection between residues %1 and %2").arg(fromResidueId).arg(toResidueId));
        }
    }
    
    return validity;
}

bool BiopolymerState::isResidueFullyMapped(int residueId) const {
    auto it = m_residues.find(residueId);
    if (it == m_residues.end()) {
        return false;
    }
    return !it->ownedAtomIds.isEmpty();
}

bool BiopolymerState::isChainFullyMapped(int chainId) const {
    auto chain = findChain(chainId);
    if (!chain) {
        return false;
    }
    
    for (int residueId : chain->residueIds) {
        if (!isResidueFullyMapped(residueId)) {
            return false;
        }
    }
    return true;
}

// ========== ID Management ==========

int BiopolymerState::allocateResidueId() {
    return m_nextResidueId++;
}

int BiopolymerState::allocateChainId() {
    return m_nextChainId++;
}

int BiopolymerState::allocateConnectionId() {
    return m_nextConnectionId++;
}

void BiopolymerState::resetIdCounters(int nextResidueId, int nextChainId, int nextConnectionId,
                                        int nextTerminalId, int nextModificationId) {
    m_nextResidueId = nextResidueId;
    m_nextChainId = nextChainId;
    m_nextConnectionId = nextConnectionId;
    m_nextTerminalId = nextTerminalId;
    m_nextModificationId = nextModificationId;
}

BiopolymerState::IdCounters BiopolymerState::idCounters() const {
    return {m_nextResidueId, m_nextChainId, m_nextConnectionId};
}

// ========== Import/Export ==========

void BiopolymerState::setImportProvenance(const ImportProvenance& provenance) {
    m_importProvenance = provenance;
}

ImportProvenance BiopolymerState::getImportProvenance() const {
    return m_importProvenance;
}

void BiopolymerState::clearImportProvenance() {
    m_importProvenance = ImportProvenance{};
}

// ========== State Query ==========

bool BiopolymerState::hasBiopolymerData() const {
    return !m_chains.isEmpty() || !m_residues.isEmpty();
}

bool BiopolymerState::isFullyValid() const {
    for (int chainId : m_chains.keys()) {
        if (!isChainFullyMapped(chainId)) {
            return false;
        }
    }
    return true;
}

QList<QString> BiopolymerState::allValidationIssues() const {
    QList<QString> issues;
    for (int chainId : m_chains.keys()) {
        MappingValidity chainValidity = validateChain(chainId);
        if (!chainValidity.isValid) {
            issues.append(chainValidity.issues);
        }
    }
    return issues;
}

// ========== Selection Support ==========

QSet<int> BiopolymerState::residuesFromAtoms(const QSet<AtomId>& atomIds) const {
    QSet<int> residueIds;
    for (AtomId atomId : atomIds) {
        int residueId = residueOwningAtom(atomId);
        if (residueId != -1) {
            residueIds.insert(residueId);
        }
    }
    return residueIds;
}

QSet<AtomId> BiopolymerState::atomsFromResidues(const QSet<int>& residueIds) const {
    QSet<AtomId> atomIds;
    for (int residueId : residueIds) {
        QList<AtomId> atoms = atomsOwnedByResidue(residueId);
        for (AtomId atomId : atoms) {
            atomIds.insert(atomId);
        }
    }
    return atomIds;
}

// ========== Polymer Type Support ==========

bool BiopolymerState::isPolymerTypeSupported(PolymerType type) {
    return type == PolymerType::Peptide || type == PolymerType::DNA || type == PolymerType::RNA;
}

QString BiopolymerState::polymerTypeToString(PolymerType type) {
    switch (type) {
    case PolymerType::Peptide: return "PEPTIDE";
    case PolymerType::DNA: return "DNA";
    case PolymerType::RNA: return "RNA";
    default: return "UNKNOWN";
    }
}

PolymerType BiopolymerState::polymerTypeFromString(const QString& str) {
    QString upper = str.toUpper();
    if (upper == "PEPTIDE" || upper == "PROTEIN") return PolymerType::Peptide;
    if (upper == "DNA") return PolymerType::DNA;
    if (upper == "RNA") return PolymerType::RNA;
    return PolymerType::Unknown;
}

// ========== Serialization Support ==========

void BiopolymerState::prepareForSerialization() {
    // Clear runtime-only data
    m_atomToResidue.clear();
    m_validationCache.clear();
}

void BiopolymerState::finalizeAfterDeserialization() {
    // Rebuild runtime data
    rebuildAtomToResidueMap();
    invalidateValidationCache();
}

// ========== Helper Methods ==========

BiopolymerChain* BiopolymerState::findChain(int chainId) {
    auto it = m_chains.find(chainId);
    return (it != m_chains.end()) ? &it.value() : nullptr;
}

const BiopolymerChain* BiopolymerState::findChain(int chainId) const {
    auto it = m_chains.find(chainId);
    return (it != m_chains.end()) ? &it.value() : nullptr;
}

void BiopolymerState::rebuildAtomToResidueMap() {
    m_atomToResidue.clear();
    for (const auto& residue : m_residues) {
        for (AtomId atomId : residue.ownedAtomIds) {
            m_atomToResidue[atomId] = residue.id;
        }
    }
}

void BiopolymerState::invalidateValidationCache() const {
    m_validationCache.clear();
}

// ========== Glue Methods for Export/Validation ==========

bool BiopolymerState::hasValidMappings() const {
    return validateAllMappings().isValid;
}

bool BiopolymerState::hasInvalidMappings() const {
    MappingValidity result = validateAllMappings();
    return !result.isValid && !result.issues.isEmpty();
}

bool BiopolymerState::hasUnsupportedModifications() const {
    for (const auto& mod : m_modifications) {
        if (!mod.isSupported) {
            return true;
        }
    }
    return false;
}

bool BiopolymerState::hasUnsupportedConnections() const {
    for (const auto& conn : m_connections) {
        if (conn.type != "peptide" && conn.type != "phosphodiester") {
            return true;
        }
    }
    return false;
}

bool BiopolymerState::hasTerminalGroups() const {
    return !m_terminals.isEmpty();
}

bool BiopolymerState::hasBranchedTopology() const {
    for (const auto& mod : m_modifications) {
        if (mod.type == ModificationType::BranchPoint) {
            return true;
        }
    }
    return false;
}

int BiopolymerState::totalResidueCount() const {
    return m_residues.size();
}

QList<PolymerType> BiopolymerState::polymerTypes() const {
    QSet<PolymerType> types;
    for (const auto& chain : m_chains) {
        types.insert(chain.polymerType);
    }
    return types.values();
}

int BiopolymerState::residuePositionInChain(int chainId, int residueId) const {
    const BiopolymerChain* chain = this->chain(chainId);
    if (!chain) return -1;
    
    const BiopolymerResidue* residue = this->residue(residueId);
    if (!residue) return -1;
    
    return residue->sequencePosition;
}

bool BiopolymerState::isNTerminal(int chainId, int residueId) const {
    const BiopolymerChain* chain = this->chain(chainId);
    if (!chain) return false;
    if (chain->residueIds.isEmpty()) return false;
    
    return chain->residueIds.first() == residueId;
}

bool BiopolymerState::isCTerminal(int chainId, int residueId) const {
    const BiopolymerChain* chain = this->chain(chainId);
    if (!chain) return false;
    if (chain->residueIds.isEmpty()) return false;
    
    return chain->residueIds.last() == residueId;
}

int BiopolymerState::chainResidueCount(int chainId) const {
    const BiopolymerChain* chain = this->chain(chainId);
    if (!chain) return 0;
    return chain->residueIds.size();
}
