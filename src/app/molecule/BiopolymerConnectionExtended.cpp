// src/app/molecule/BiopolymerConnectionExtended.cpp
#include "BiopolymerConnectionExtended.h"

// ========== Helper Functions ==========

QString connectionSubtypeToString(ConnectionSubtype subtype) {
    switch (subtype) {
        case ConnectionSubtype::Peptide: return "Peptide";
        case ConnectionSubtype::Phosphodiester: return "Phosphodiester";
        case ConnectionSubtype::Disulfide: return "Disulfide";
        case ConnectionSubtype::Thioester: return "Thioester";
        case ConnectionSubtype::Lactam: return "Lactam";
        case ConnectionSubtype::Custom: return "Custom";
        case ConnectionSubtype::Cyclization: return "Cyclization";
        case ConnectionSubtype::Branch: return "Branch";
        case ConnectionSubtype::Interchain: return "Interchain";
        default: return "Unknown";
    }
}

ConnectionSubtype connectionSubtypeFromString(const QString& str) {
    static const QMap<QString, ConnectionSubtype> map = {
        {"peptide", ConnectionSubtype::Peptide},
        {"phosphodiester", ConnectionSubtype::Phosphodiester},
        {"disulfide", ConnectionSubtype::Disulfide},
        {"thioester", ConnectionSubtype::Thioester},
        {"lactam", ConnectionSubtype::Lactam},
        {"custom", ConnectionSubtype::Custom},
        {"cyclization", ConnectionSubtype::Cyclization},
        {"branch", ConnectionSubtype::Branch},
        {"interchain", ConnectionSubtype::Interchain}
    };
    QString lower = str.toLower();
    auto it = map.find(lower);
    return (it != map.end()) ? it.value() : ConnectionSubtype::Unknown;
}

// ========== BiopolymerExtendedConnection Implementation ==========

BiopolymerExtendedConnection BiopolymerExtendedConnection::createDisulfide(
    int connId, const BiopolymerDisulfide& ds) {
    BiopolymerExtendedConnection conn(connId, ConnectionSubtype::Disulfide);
    conn.disulfide = ds;
    conn.type = "Disulfide";
    conn.name = ds.name;
    conn.description = ds.description;
    conn.isValid = ds.isValid;
    conn.validationError = ds.validationError;
    conn.creationTimestamp = ds.creationTimestamp;
    conn.residueIds = {ds.residueId1, ds.residueId2};
    conn.chainIds = {ds.chainId1, ds.chainId2};
    if (ds.atomId1 >= 0) conn.atomIds.append(ds.atomId1);
    if (ds.atomId2 >= 0) conn.atomIds.append(ds.atomId2);
    if (ds.bondId >= 0) conn.bondIds.append(ds.bondId);
    conn.isSelected = ds.isSelected;
    conn.isVisible = ds.isVisible;
    return conn;
}

BiopolymerExtendedConnection BiopolymerExtendedConnection::createInterchain(
    int connId, const BiopolymerInterchainConnection& ic) {
    BiopolymerExtendedConnection conn(connId, ConnectionSubtype::Interchain);
    conn.interchain = ic;
    conn.type = ic.type;
    conn.name = ic.name;
    conn.description = ic.description;
    conn.isValid = ic.isValid;
    conn.validationError = ic.validationError;
    conn.creationTimestamp = ic.creationTimestamp;
    conn.residueIds = {ic.residueId1, ic.residueId2};
    conn.chainIds = {ic.chainId1, ic.chainId2};
    if (ic.atomId1 >= 0) conn.atomIds.append(ic.atomId1);
    if (ic.atomId2 >= 0) conn.atomIds.append(ic.atomId2);
    if (ic.bondId >= 0) conn.bondIds.append(ic.bondId);
    conn.isSelected = ic.isSelected;
    conn.isVisible = ic.isVisible;
    return conn;
}

BiopolymerExtendedConnection BiopolymerExtendedConnection::createBranch(
    int connId, const BiopolymerBranch& br) {
    BiopolymerExtendedConnection conn(connId, ConnectionSubtype::Branch);
    conn.branch = br;
    conn.type = br.type;
    conn.name = br.name;
    conn.description = br.description;
    conn.isValid = br.isValid;
    conn.validationError = br.validationError;
    conn.creationTimestamp = br.creationTimestamp;
    conn.residueIds = {br.residueId, br.branchResidueId};
    conn.chainIds = {br.chainId, br.branchChainId};
    if (br.atomId >= 0) conn.atomIds.append(br.atomId);
    if (br.bondId >= 0) conn.bondIds.append(br.bondId);
    conn.isSelected = br.isSelected;
    conn.isVisible = br.isVisible;
    return conn;
}

BiopolymerExtendedConnection BiopolymerExtendedConnection::createCyclization(
    int connId, const BiopolymerCyclization& cyc) {
    BiopolymerExtendedConnection conn(connId, ConnectionSubtype::Cyclization);
    conn.cyclization = cyc;
    conn.type = cyc.type;
    conn.name = cyc.name;
    conn.description = cyc.description;
    conn.isValid = cyc.isValid;
    conn.validationError = cyc.validationError;
    conn.creationTimestamp = cyc.creationTimestamp;
    conn.residueIds = cyc.residueIds;
    conn.chainIds = cyc.chainIds;
    conn.atomIds = cyc.atomIds;
    conn.bondIds = cyc.bondIds;
    conn.isSelected = cyc.isSelected;
    conn.isVisible = cyc.isVisible;
    return conn;
}

bool BiopolymerExtendedConnection::isValidConnection() const {
    if (!isValid) return false;
    
    switch (subtype) {
        case ConnectionSubtype::Disulfide:
            return disulfide.isValidDisulfide();
        case ConnectionSubtype::Interchain:
            return interchain.isValidConnection();
        case ConnectionSubtype::Branch:
            return branch.isValidBranch();
        case ConnectionSubtype::Cyclization:
            return cyclization.isValidCyclization();
        default:
            return false;
    }
}

bool BiopolymerExtendedConnection::involvesResidue(int residueId) const {
    return residueIds.contains(residueId);
}

bool BiopolymerExtendedConnection::involvesChain(int chainId) const {
    return chainIds.contains(chainId);
}

bool BiopolymerExtendedConnection::involvesAtom(int atomId) const {
    return atomIds.contains(atomId);
}

const BiopolymerDisulfide& BiopolymerExtendedConnection::asDisulfide() const {
    if (subtype != ConnectionSubtype::Disulfide) {
        static BiopolymerDisulfide empty;
        return empty;
    }
    return disulfide;
}

BiopolymerDisulfide& BiopolymerExtendedConnection::asDisulfide() {
    if (subtype != ConnectionSubtype::Disulfide) {
        static BiopolymerDisulfide empty;
        return empty;
    }
    return disulfide;
}

const BiopolymerInterchainConnection& BiopolymerExtendedConnection::asInterchain() const {
    if (subtype != ConnectionSubtype::Interchain) {
        static BiopolymerInterchainConnection empty;
        return empty;
    }
    return interchain;
}

BiopolymerInterchainConnection& BiopolymerExtendedConnection::asInterchain() {
    if (subtype != ConnectionSubtype::Interchain) {
        static BiopolymerInterchainConnection empty;
        return empty;
    }
    return interchain;
}

const BiopolymerBranch& BiopolymerExtendedConnection::asBranch() const {
    if (subtype != ConnectionSubtype::Branch) {
        static BiopolymerBranch empty;
        return empty;
    }
    return branch;
}

BiopolymerBranch& BiopolymerExtendedConnection::asBranch() {
    if (subtype != ConnectionSubtype::Branch) {
        static BiopolymerBranch empty;
        return empty;
    }
    return branch;
}

const BiopolymerCyclization& BiopolymerExtendedConnection::asCyclization() const {
    if (subtype != ConnectionSubtype::Cyclization) {
        static BiopolymerCyclization empty;
        return empty;
    }
    return cyclization;
}

BiopolymerCyclization& BiopolymerExtendedConnection::asCyclization() {
    if (subtype != ConnectionSubtype::Cyclization) {
        static BiopolymerCyclization empty;
        return empty;
    }
    return cyclization;
}

// ========== ExtendedConnectionValidator Implementation ==========

ExtendedConnectionValidator::ExtendedConnectionValidator(const BiopolymerState& state)
    : m_state(state) {}

ExtendedConnectionValidator::~ExtendedConnectionValidator() = default;

bool ExtendedConnectionValidator::validateDisulfide(const BiopolymerDisulfide& disulfide) const {
    if (!checkResidueExistence(disulfide.residueId1) || 
        !checkResidueExistence(disulfide.residueId2)) {
        return false;
    }
    
    if (disulfide.residueId1 == disulfide.residueId2) {
        return false;  // Can't have disulfide with same residue
    }
    
    if (!checkBondExistence(disulfide.bondId)) {
        return false;
    }
    
    if (!areResiduesCompatibleForDisulfide(disulfide.residueId1, disulfide.residueId2)) {
        return false;
    }
    
    if (wouldConflict(disulfide)) {
        return false;
    }
    
    return true;
}

bool ExtendedConnectionValidator::validateInterchain(const BiopolymerInterchainConnection& connection) const {
    if (!checkResidueExistence(connection.residueId1) || 
        !checkResidueExistence(connection.residueId2)) {
        return false;
    }
    
    if (!checkChainExistence(connection.chainId1) || 
        !checkChainExistence(connection.chainId2)) {
        return false;
    }
    
    if (connection.chainId1 == connection.chainId2) {
        return false;  // Interchain should be between different chains
    }
    
    if (!checkBondExistence(connection.bondId)) {
        return false;
    }
    
    if (!areResiduesCompatibleForConnection(connection.residueId1, connection.residueId2)) {
        return false;
    }
    
    if (wouldConflict(connection)) {
        return false;
    }
    
    return true;
}

bool ExtendedConnectionValidator::validateBranch(const BiopolymerBranch& branch) const {
    if (!checkResidueExistence(branch.residueId) || 
        !checkResidueExistence(branch.branchResidueId)) {
        return false;
    }
    
    if (!checkChainExistence(branch.chainId) || 
        !checkChainExistence(branch.branchChainId)) {
        return false;
    }
    
    if (!checkBondExistence(branch.bondId)) {
        return false;
    }
    
    if (wouldConflict(branch)) {
        return false;
    }
    
    return true;
}

bool ExtendedConnectionValidator::validateCyclization(const BiopolymerCyclization& cyclization) const {
    if (cyclization.residueIds.size() < 3) {
        return false;  // Need at least 3 residues for a ring
    }
    
    for (int residueId : cyclization.residueIds) {
        if (!checkResidueExistence(residueId)) {
            return false;
        }
    }
    
    for (BondId bondId : cyclization.bondIds) {
        if (!checkBondExistence(bondId)) {
            return false;
        }
    }
    
    if (wouldConflict(cyclization)) {
        return false;
    }
    
    return true;
}

bool ExtendedConnectionValidator::validateExtendedConnection(const BiopolymerExtendedConnection& connection) const {
    switch (connection.subtype) {
        case ConnectionSubtype::Disulfide:
            return validateDisulfide(connection.asDisulfide());
        case ConnectionSubtype::Interchain:
            return validateInterchain(connection.asInterchain());
        case ConnectionSubtype::Branch:
            return validateBranch(connection.asBranch());
        case ConnectionSubtype::Cyclization:
            return validateCyclization(connection.asCyclization());
        default:
            return false;
    }
}

QList<QString> ExtendedConnectionValidator::validateAllDisulfides() const {
    QList<QString> issues;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerDisulfide& ds : manager.allDisulfides()) {
        if (!validateDisulfide(ds)) {
            issues.append(QString("Invalid disulfide %1 between residues %2 and %3")
                .arg(ds.disulfideId).arg(ds.residueId1).arg(ds.residueId2));
        }
    }
    
    return issues;
}

QList<QString> ExtendedConnectionValidator::validateAllInterchainConnections() const {
    QList<QString> issues;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerInterchainConnection& ic : manager.allInterchainConnections()) {
        if (!validateInterchain(ic)) {
            issues.append(QString("Invalid interchain connection %1 between chains %2 and %3")
                .arg(ic.connectionId).arg(ic.chainId1).arg(ic.chainId2));
        }
    }
    
    return issues;
}

QList<QString> ExtendedConnectionValidator::validateAllBranches() const {
    QList<QString> issues;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerBranch& br : manager.allBranches()) {
        if (!validateBranch(br)) {
            issues.append(QString("Invalid branch %1 at residue %2")
                .arg(br.branchId).arg(br.residueId));
        }
    }
    
    return issues;
}

QList<QString> ExtendedConnectionValidator::validateAllCyclizations() const {
    QList<QString> issues;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerCyclization& cyc : manager.allCyclizations()) {
        if (!validateCyclization(cyc)) {
            issues.append(QString("Invalid cyclization %1 with %2 residues")
                .arg(cyc.cyclizationId).arg(cyc.residueIds.size()));
        }
    }
    
    return issues;
}

QList<QString> ExtendedConnectionValidator::validateAllExtendedConnections() const {
    QList<QString> issues;
    issues.append(validateAllDisulfides());
    issues.append(validateAllInterchainConnections());
    issues.append(validateAllBranches());
    issues.append(validateAllCyclizations());
    return issues;
}

bool ExtendedConnectionValidator::canFormDisulfide(int residueId1, int residueId2) const {
    return areResiduesCompatibleForDisulfide(residueId1, residueId2);
}

bool ExtendedConnectionValidator::canFormInterchainConnection(int residueId1, int residueId2) const {
    return areResiduesCompatibleForConnection(residueId1, residueId2);
}

bool ExtendedConnectionValidator::canCreateBranch(int chainId, int residueId, int branchChainId) const {
    return checkChainExistence(chainId) && 
           checkChainExistence(branchChainId) &&
           checkResidueExistence(residueId);
}

bool ExtendedConnectionValidator::canCreateCyclization(const QList<int>& residueIds) const {
    if (residueIds.size() < 3) return false;
    
    for (int residueId : residueIds) {
        if (!checkResidueExistence(residueId)) {
            return false;
        }
    }
    
    return true;
}

bool ExtendedConnectionValidator::wouldConflict(const BiopolymerDisulfide& disulfide) const {
    return !getConflictingDisulfides(disulfide).isEmpty();
}

bool ExtendedConnectionValidator::wouldConflict(const BiopolymerInterchainConnection& connection) const {
    return !getConflictingInterchainConnections(connection).isEmpty();
}

bool ExtendedConnectionValidator::wouldConflict(const BiopolymerBranch& branch) const {
    return !getConflictingBranches(branch).isEmpty();
}

bool ExtendedConnectionValidator::wouldConflict(const BiopolymerCyclization& cyclization) const {
    // Check for conflicts with existing connections
    for (int residueId : cyclization.residueIds) {
        const BiopolymerResidue* residue = m_state.residue(residueId);
        if (!residue) continue;
        
        // Check if residue is already in too many connections
        // This is a simplified check
    }
    return false;
}

QList<int> ExtendedConnectionValidator::getConflictingDisulfides(const BiopolymerDisulfide& disulfide) const {
    QList<int> conflicts;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerDisulfide& existing : manager.allDisulfides()) {
        // Check for same residue involvement
        if ((existing.residueId1 == disulfide.residueId1 || existing.residueId1 == disulfide.residueId2 ||
             existing.residueId2 == disulfide.residueId1 || existing.residueId2 == disulfide.residueId2) &&
            existing.disulfideId != disulfide.disulfideId) {
            conflicts.append(existing.disulfideId);
        }
    }
    
    return conflicts;
}

QList<int> ExtendedConnectionValidator::getConflictingInterchainConnections(
    const BiopolymerInterchainConnection& connection) const {
    QList<int> conflicts;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerInterchainConnection& existing : manager.allInterchainConnections()) {
        // Check for same chain pair
        if (((existing.chainId1 == connection.chainId1 && existing.chainId2 == connection.chainId2) ||
             (existing.chainId1 == connection.chainId2 && existing.chainId2 == connection.chainId1)) &&
            existing.connectionId != connection.connectionId) {
            conflicts.append(existing.connectionId);
        }
    }
    
    return conflicts;
}

QList<int> ExtendedConnectionValidator::getConflictingBranches(const BiopolymerBranch& branch) const {
    QList<int> conflicts;
    ExtendedConnectionManager manager(const_cast<BiopolymerState&>(m_state));
    
    for (const BiopolymerBranch& existing : manager.allBranches()) {
        // Check for same residue
        if (existing.residueId == branch.residueId && existing.branchId != branch.branchId) {
            conflicts.append(existing.branchId);
        }
    }
    
    return conflicts;
}

QList<BiopolymerDisulfide> ExtendedConnectionValidator::getPossibleDisulfides(int chainId) const {
    QList<BiopolymerDisulfide> possible;
    
    const BiopolymerChain* chain = m_state.chain(chainId);
    if (!chain) return possible;
    
    // Find all cysteine residues in the chain
    QList<int> cysteineResidues;
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = m_state.residue(residueId);
        if (residue && residue->identity.symbol == "C") {  // Cysteine
            cysteineResidues.append(residueId);
        }
    }
    
    // Generate possible disulfide pairs
    for (int i = 0; i < cysteineResidues.size(); ++i) {
        for (int j = i + 1; j < cysteineResidues.size(); ++j) {
            int resId1 = cysteineResidues[i];
            int resId2 = cysteineResidues[j];
            
            if (canFormDisulfide(resId1, resId2)) {
                BiopolymerDisulfide ds;
                ds.residueId1 = resId1;
                ds.residueId2 = resId2;
                ds.chainId1 = chainId;
                ds.chainId2 = chainId;
                ds.isIntrachain = true;
                ds.isValid = true;
                ds.name = QString("Disulfide_%1_%2").arg(resId1).arg(resId2);
                possible.append(ds);
            }
        }
    }
    
    return possible;
}

QList<BiopolymerInterchainConnection> ExtendedConnectionValidator::getPossibleInterchainConnections() const {
    QList<BiopolymerInterchainConnection> possible;
    
    QList<int> chainIds = m_state.chainIds();
    for (int i = 0; i < chainIds.size(); ++i) {
        for (int j = i + 1; j < chainIds.size(); ++j) {
            int chainId1 = chainIds[i];
            int chainId2 = chainIds[j];
            
            const BiopolymerChain* chain1 = m_state.chain(chainId1);
            const BiopolymerChain* chain2 = m_state.chain(chainId2);
            
            if (!chain1 || !chain2) continue;
            
            // Try to find compatible residues
            for (int residueId1 : chain1->residueIds) {
                for (int residueId2 : chain2->residueIds) {
                    if (canFormInterchainConnection(residueId1, residueId2)) {
                        BiopolymerInterchainConnection ic;
                        ic.chainId1 = chainId1;
                        ic.residueId1 = residueId1;
                        ic.chainId2 = chainId2;
                        ic.residueId2 = residueId2;
                        ic.isValid = true;
                        ic.name = QString("Interchain_%1_%2").arg(chainId1).arg(chainId2);
                        possible.append(ic);
                    }
                }
            }
        }
    }
    
    return possible;
}

bool ExtendedConnectionValidator::areResiduesCompatibleForDisulfide(int residueId1, int residueId2) const {
    const BiopolymerResidue* res1 = m_state.residue(residueId1);
    const BiopolymerResidue* res2 = m_state.residue(residueId2);
    
    if (!res1 || !res2) return false;
    
    // For now, just check if both are cysteine (C)
    // In reality, this would check for sulfur-containing residues
    return res1->identity.symbol == "C" && res2->identity.symbol == "C";
}

bool ExtendedConnectionValidator::areResiduesCompatibleForConnection(int residueId1, int residueId2) const {
    const BiopolymerResidue* res1 = m_state.residue(residueId1);
    const BiopolymerResidue* res2 = m_state.residue(residueId2);
    
    if (!res1 || !res2) return false;
    
    // Check if residues can form a bond
    // This would depend on the chemistry
    return true;  // Simplified for now
}

bool ExtendedConnectionValidator::checkResidueExistence(int residueId) const {
    return m_state.residue(residueId) != nullptr;
}

bool ExtendedConnectionValidator::checkChainExistence(int chainId) const {
    return m_state.chain(chainId) != nullptr;
}

bool ExtendedConnectionValidator::checkAtomExistence(int atomId) const {
    // For now, we can't check atom existence directly from BiopolymerState
    // This would need access to the underlying molecule
    return atomId >= 0;
}

bool ExtendedConnectionValidator::checkBondExistence(BondId bondId) const {
    // For now, we can't check bond existence directly from BiopolymerState
    // This would need access to the underlying molecule
    return bondId >= 0;
}

bool ExtendedConnectionValidator::checkResidueType(int residueId, const QSet<QString>& allowedTypes) const {
    const BiopolymerResidue* residue = m_state.residue(residueId);
    if (!residue) return false;
    return allowedTypes.contains(residue->identity.symbol);
}

// ========== ExtendedConnectionManager Implementation ==========

ExtendedConnectionManager::ExtendedConnectionManager(BiopolymerState& state)
    : m_state(state) {}

ExtendedConnectionManager::~ExtendedConnectionManager() = default;

// ========== Disulfide Management ==========

int ExtendedConnectionManager::addDisulfide(const BiopolymerDisulfide& disulfide) {
    // Check if ID is valid
    int disulfidId = disulfide.disulfideId;
    if (disulfidId < 0) {
        disulfidId = allocateDisulfideId();
    }
    
    // Validate
    ExtendedConnectionValidator validator(m_state);
    if (!validator.validateDisulfide(disulfide)) {
        return -1;
    }
    
    // Add to storage
    BiopolymerDisulfide newDisulfide = disulfide;
    newDisulfide.disulfideId = disulfidId;
    m_disulfides[disulfidId] = newDisulfide;
    
    // Add to extended connections
    BiopolymerExtendedConnection extended = BiopolymerExtendedConnection::createDisulfide(disulfidId, newDisulfide);
    m_extendedConnections[disulfidId] = extended;
    
    return disulfidId;
}

int ExtendedConnectionManager::addDisulfide(int residueId1, int residueId2, BondId bondId, 
                                           int atomId1, int atomId2) {
    BiopolymerDisulfide disulfide;
    disulfide.residueId1 = residueId1;
    disulfide.residueId2 = residueId2;
    
    const BiopolymerResidue* res1 = m_state.residue(residueId1);
    const BiopolymerResidue* res2 = m_state.residue(residueId2);
    
    if (res1) {
        disulfide.chainId1 = res1->chainId;
    }
    if (res2) {
        disulfide.chainId2 = res2->chainId;
    }
    
    disulfide.isIntrachain = (disulfide.chainId1 == disulfide.chainId2);
    disulfide.atomId1 = atomId1;
    disulfide.atomId2 = atomId2;
    disulfide.bondId = bondId;
    disulfide.isValid = true;
    disulfide.name = QString("Disulfide_%1_%2").arg(residueId1).arg(residueId2);
    
    return addDisulfide(disulfide);
}

bool ExtendedConnectionManager::removeDisulfide(int disulfidId) {
    auto it = m_disulfides.find(disulfidId);
    if (it == m_disulfides.end()) {
        return false;
    }
    
    m_disulfides.erase(it);
    m_extendedConnections.remove(disulfidId);
    return true;
}

void ExtendedConnectionManager::removeDisulfidesFromResidue(int residueId) {
    QList<int> toRemove;
    for (auto it = m_disulfides.begin(); it != m_disulfides.end(); ++it) {
        if (it->residueId1 == residueId || it->residueId2 == residueId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_disulfides.remove(id);
        m_extendedConnections.remove(id);
    }
}

void ExtendedConnectionManager::removeDisulfidesFromChain(int chainId) {
    QList<int> toRemove;
    for (auto it = m_disulfides.begin(); it != m_disulfides.end(); ++it) {
        if (it->chainId1 == chainId || it->chainId2 == chainId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_disulfides.remove(id);
        m_extendedConnections.remove(id);
    }
}

BiopolymerDisulfide* ExtendedConnectionManager::getDisulfide(int disulfidId) {
    auto it = m_disulfides.find(disulfidId);
    return (it != m_disulfides.end()) ? &it.value() : nullptr;
}

const BiopolymerDisulfide* ExtendedConnectionManager::getDisulfide(int disulfidId) const {
    auto it = m_disulfides.find(disulfidId);
    return (it != m_disulfides.end()) ? &it.value() : nullptr;
}

QList<BiopolymerDisulfide> ExtendedConnectionManager::allDisulfides() const {
    return m_disulfides.values();
}

QList<BiopolymerDisulfide> ExtendedConnectionManager::disulfidesForResidue(int residueId) const {
    QList<BiopolymerDisulfide> result;
    for (const auto& ds : m_disulfides) {
        if (ds.residueId1 == residueId || ds.residueId2 == residueId) {
            result.append(ds);
        }
    }
    return result;
}

QList<BiopolymerDisulfide> ExtendedConnectionManager::disulfidesForChain(int chainId) const {
    QList<BiopolymerDisulfide> result;
    for (const auto& ds : m_disulfides) {
        if (ds.chainId1 == chainId || ds.chainId2 == chainId) {
            result.append(ds);
        }
    }
    return result;
}

bool ExtendedConnectionManager::hasDisulfides(int residueId) const {
    for (const auto& ds : m_disulfides) {
        if (ds.residueId1 == residueId || ds.residueId2 == residueId) {
            return true;
        }
    }
    return false;
}

bool ExtendedConnectionManager::hasChainDisulfides(int chainId) const {
    for (const auto& ds : m_disulfides) {
        if (ds.chainId1 == chainId || ds.chainId2 == chainId) {
            return true;
        }
    }
    return false;
}

QList<BiopolymerDisulfide> ExtendedConnectionManager::interchainDisulfides() const {
    QList<BiopolymerDisulfide> result;
    for (const auto& ds : m_disulfides) {
        if (!ds.isIntrachain) {
            result.append(ds);
        }
    }
    return result;
}

QList<BiopolymerDisulfide> ExtendedConnectionManager::intrachainDisulfides() const {
    QList<BiopolymerDisulfide> result;
    for (const auto& ds : m_disulfides) {
        if (ds.isIntrachain) {
            result.append(ds);
        }
    }
    return result;
}

BiopolymerDisulfide* ExtendedConnectionManager::findDisulfide(int residueId1, int residueId2) {
    for (auto& ds : m_disulfides) {
        if (ds.connects(residueId1, residueId2)) {
            return &ds;
        }
    }
    return nullptr;
}

const BiopolymerDisulfide* ExtendedConnectionManager::findDisulfide(int residueId1, int residueId2) const {
    for (const auto& ds : m_disulfides) {
        if (ds.connects(residueId1, residueId2)) {
            return &ds;
        }
    }
    return nullptr;
}

// ========== Interchain Connection Management ==========

int ExtendedConnectionManager::addInterchainConnection(const BiopolymerInterchainConnection& connection) {
    int connectionId = connection.connectionId;
    if (connectionId < 0) {
        connectionId = allocateInterchainConnectionId();
    }
    
    ExtendedConnectionValidator validator(m_state);
    if (!validator.validateInterchain(connection)) {
        return -1;
    }
    
    BiopolymerInterchainConnection newConnection = connection;
    newConnection.connectionId = connectionId;
    m_interchainConnections[connectionId] = newConnection;
    
    BiopolymerExtendedConnection extended = BiopolymerExtendedConnection::createInterchain(connectionId, newConnection);
    m_extendedConnections[connectionId] = extended;
    
    return connectionId;
}

int ExtendedConnectionManager::addInterchainConnection(int chainId1, int residueId1, 
                                                        int chainId2, int residueId2,
                                                        BondId bondId, int atomId1, int atomId2) {
    BiopolymerInterchainConnection connection;
    connection.chainId1 = chainId1;
    connection.residueId1 = residueId1;
    connection.chainId2 = chainId2;
    connection.residueId2 = residueId2;
    connection.bondId = bondId;
    connection.atomId1 = atomId1;
    connection.atomId2 = atomId2;
    connection.isValid = true;
    connection.name = QString("Interchain_%1_%2").arg(chainId1).arg(chainId2);
    
    return addInterchainConnection(connection);
}

bool ExtendedConnectionManager::removeInterchainConnection(int connectionId) {
    auto it = m_interchainConnections.find(connectionId);
    if (it == m_interchainConnections.end()) {
        return false;
    }
    
    m_interchainConnections.erase(it);
    m_extendedConnections.remove(connectionId);
    return true;
}

void ExtendedConnectionManager::removeInterchainConnectionsFromChain(int chainId) {
    QList<int> toRemove;
    for (auto it = m_interchainConnections.begin(); it != m_interchainConnections.end(); ++it) {
        if (it->chainId1 == chainId || it->chainId2 == chainId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_interchainConnections.remove(id);
        m_extendedConnections.remove(id);
    }
}

void ExtendedConnectionManager::removeInterchainConnectionsFromResidue(int residueId) {
    QList<int> toRemove;
    for (auto it = m_interchainConnections.begin(); it != m_interchainConnections.end(); ++it) {
        if (it->residueId1 == residueId || it->residueId2 == residueId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_interchainConnections.remove(id);
        m_extendedConnections.remove(id);
    }
}

BiopolymerInterchainConnection* ExtendedConnectionManager::getInterchainConnection(int connectionId) {
    auto it = m_interchainConnections.find(connectionId);
    return (it != m_interchainConnections.end()) ? &it.value() : nullptr;
}

const BiopolymerInterchainConnection* ExtendedConnectionManager::getInterchainConnection(int connectionId) const {
    auto it = m_interchainConnections.find(connectionId);
    return (it != m_interchainConnections.end()) ? &it.value() : nullptr;
}

QList<BiopolymerInterchainConnection> ExtendedConnectionManager::allInterchainConnections() const {
    return m_interchainConnections.values();
}

QList<BiopolymerInterchainConnection> ExtendedConnectionManager::interchainConnectionsForChain(int chainId) const {
    QList<BiopolymerInterchainConnection> result;
    for (const auto& ic : m_interchainConnections) {
        if (ic.chainId1 == chainId || ic.chainId2 == chainId) {
            result.append(ic);
        }
    }
    return result;
}

QList<BiopolymerInterchainConnection> ExtendedConnectionManager::interchainConnectionsBetween(int chainId1, int chainId2) const {
    QList<BiopolymerInterchainConnection> result;
    for (const auto& ic : m_interchainConnections) {
        if ((ic.chainId1 == chainId1 && ic.chainId2 == chainId2) ||
            (ic.chainId1 == chainId2 && ic.chainId2 == chainId1)) {
            result.append(ic);
        }
    }
    return result;
}

bool ExtendedConnectionManager::areChainsConnected(int chainId1, int chainId2) const {
    for (const auto& ic : m_interchainConnections) {
        if ((ic.chainId1 == chainId1 && ic.chainId2 == chainId2) ||
            (ic.chainId1 == chainId2 && ic.chainId2 == chainId1)) {
            return true;
        }
    }
    return false;
}

BiopolymerInterchainConnection* ExtendedConnectionManager::findInterchainConnection(int residueId1, int residueId2) {
    for (auto& ic : m_interchainConnections) {
        if (ic.connects(residueId1, residueId2)) {
            return &ic;
        }
    }
    return nullptr;
}

const BiopolymerInterchainConnection* ExtendedConnectionManager::findInterchainConnection(int residueId1, int residueId2) const {
    for (const auto& ic : m_interchainConnections) {
        if (ic.connects(residueId1, residueId2)) {
            return &ic;
        }
    }
    return nullptr;
}

// ========== Branch Management ==========

int ExtendedConnectionManager::addBranch(const BiopolymerBranch& branch) {
    int branchId = branch.branchId;
    if (branchId < 0) {
        branchId = allocateBranchId();
    }
    
    ExtendedConnectionValidator validator(m_state);
    if (!validator.validateBranch(branch)) {
        return -1;
    }
    
    BiopolymerBranch newBranch = branch;
    newBranch.branchId = branchId;
    m_branches[branchId] = newBranch;
    
    BiopolymerExtendedConnection extended = BiopolymerExtendedConnection::createBranch(branchId, newBranch);
    m_extendedConnections[branchId] = extended;
    
    return branchId;
}

int ExtendedConnectionManager::addBranch(int chainId, int residueId, int branchChainId, 
                                        int branchResidueId, BondId bondId, int atomId) {
    BiopolymerBranch branch;
    branch.chainId = chainId;
    branch.residueId = residueId;
    branch.branchChainId = branchChainId;
    branch.branchResidueId = branchResidueId;
    branch.bondId = bondId;
    branch.atomId = atomId;
    branch.isValid = true;
    branch.name = QString("Branch_%1_%2").arg(chainId).arg(residueId);
    
    return addBranch(branch);
}

bool ExtendedConnectionManager::removeBranch(int branchId) {
    auto it = m_branches.find(branchId);
    if (it == m_branches.end()) {
        return false;
    }
    
    m_branches.erase(it);
    m_extendedConnections.remove(branchId);
    return true;
}

void ExtendedConnectionManager::removeBranchesFromChain(int chainId) {
    QList<int> toRemove;
    for (auto it = m_branches.begin(); it != m_branches.end(); ++it) {
        if (it->chainId == chainId || it->branchChainId == chainId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_branches.remove(id);
        m_extendedConnections.remove(id);
    }
}

void ExtendedConnectionManager::removeBranchesFromResidue(int residueId) {
    QList<int> toRemove;
    for (auto it = m_branches.begin(); it != m_branches.end(); ++it) {
        if (it->residueId == residueId || it->branchResidueId == residueId) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_branches.remove(id);
        m_extendedConnections.remove(id);
    }
}

BiopolymerBranch* ExtendedConnectionManager::getBranch(int branchId) {
    auto it = m_branches.find(branchId);
    return (it != m_branches.end()) ? &it.value() : nullptr;
}

const BiopolymerBranch* ExtendedConnectionManager::getBranch(int branchId) const {
    auto it = m_branches.find(branchId);
    return (it != m_branches.end()) ? &it.value() : nullptr;
}

QList<BiopolymerBranch> ExtendedConnectionManager::allBranches() const {
    return m_branches.values();
}

QList<BiopolymerBranch> ExtendedConnectionManager::branchesForChain(int chainId) const {
    QList<BiopolymerBranch> result;
    for (const auto& br : m_branches) {
        if (br.chainId == chainId || br.branchChainId == chainId) {
            result.append(br);
        }
    }
    return result;
}

QList<BiopolymerBranch> ExtendedConnectionManager::branchesFromResidue(int residueId) const {
    QList<BiopolymerBranch> result;
    for (const auto& br : m_branches) {
        if (br.residueId == residueId) {
            result.append(br);
        }
    }
    return result;
}

bool ExtendedConnectionManager::isBranchPoint(int residueId) const {
    for (const auto& br : m_branches) {
        if (br.residueId == residueId) {
            return true;
        }
    }
    return false;
}

bool ExtendedConnectionManager::hasChainBranches(int chainId) const {
    for (const auto& br : m_branches) {
        if (br.chainId == chainId || br.branchChainId == chainId) {
            return true;
        }
    }
    return false;
}

// ========== Cyclization Management ==========

int ExtendedConnectionManager::addCyclization(const BiopolymerCyclization& cyclization) {
    int cyclizationId = cyclization.cyclizationId;
    if (cyclizationId < 0) {
        cyclizationId = allocateCyclizationId();
    }
    
    ExtendedConnectionValidator validator(m_state);
    if (!validator.validateCyclization(cyclization)) {
        return -1;
    }
    
    BiopolymerCyclization newCyclization = cyclization;
    newCyclization.cyclizationId = cyclizationId;
    m_cyclizations[cyclizationId] = newCyclization;
    
    BiopolymerExtendedConnection extended = BiopolymerExtendedConnection::createCyclization(cyclizationId, newCyclization);
    m_extendedConnections[cyclizationId] = extended;
    
    return cyclizationId;
}

bool ExtendedConnectionManager::removeCyclization(int cyclizationId) {
    auto it = m_cyclizations.find(cyclizationId);
    if (it == m_cyclizations.end()) {
        return false;
    }
    
    m_cyclizations.erase(it);
    m_extendedConnections.remove(cyclizationId);
    return true;
}

void ExtendedConnectionManager::removeCyclizationsFromResidue(int residueId) {
    QList<int> toRemove;
    for (auto it = m_cyclizations.begin(); it != m_cyclizations.end(); ++it) {
        if (it->residueIds.contains(residueId)) {
            toRemove.append(it.key());
        }
    }
    
    for (int id : toRemove) {
        m_cyclizations.remove(id);
        m_extendedConnections.remove(id);
    }
}

BiopolymerCyclization* ExtendedConnectionManager::getCyclization(int cyclizationId) {
    auto it = m_cyclizations.find(cyclizationId);
    return (it != m_cyclizations.end()) ? &it.value() : nullptr;
}

const BiopolymerCyclization* ExtendedConnectionManager::getCyclization(int cyclizationId) const {
    auto it = m_cyclizations.find(cyclizationId);
    return (it != m_cyclizations.end()) ? &it.value() : nullptr;
}

QList<BiopolymerCyclization> ExtendedConnectionManager::allCyclizations() const {
    return m_cyclizations.values();
}

QList<BiopolymerCyclization> ExtendedConnectionManager::cyclizationsForResidue(int residueId) const {
    QList<BiopolymerCyclization> result;
    for (const auto& cyc : m_cyclizations) {
        if (cyc.residueIds.contains(residueId)) {
            result.append(cyc);
        }
    }
    return result;
}

bool ExtendedConnectionManager::isInCyclization(int residueId) const {
    for (const auto& cyc : m_cyclizations) {
        if (cyc.residueIds.contains(residueId)) {
            return true;
        }
    }
    return false;
}

// ========== Extended Connection Management ==========

int ExtendedConnectionManager::addExtendedConnection(const BiopolymerExtendedConnection& connection) {
    int connectionId = connection.id;
    if (connectionId < 0) {
        connectionId = allocateExtendedConnectionId();
    }
    
    ExtendedConnectionValidator validator(m_state);
    if (!validator.validateExtendedConnection(connection)) {
        return -1;
    }
    
    BiopolymerExtendedConnection newConnection = connection;
    newConnection.id = connectionId;
    m_extendedConnections[connectionId] = newConnection;
    
    // Also add to specific storage based on type
    switch (connection.subtype) {
        case ConnectionSubtype::Disulfide:
            m_disulfides[connectionId] = connection.asDisulfide();
            break;
        case ConnectionSubtype::Interchain:
            m_interchainConnections[connectionId] = connection.asInterchain();
            break;
        case ConnectionSubtype::Branch:
            m_branches[connectionId] = connection.asBranch();
            break;
        case ConnectionSubtype::Cyclization:
            m_cyclizations[connectionId] = connection.asCyclization();
            break;
        default:
            break;
    }
    
    return connectionId;
}

bool ExtendedConnectionManager::removeExtendedConnection(int connectionId) {
    auto it = m_extendedConnections.find(connectionId);
    if (it == m_extendedConnections.end()) {
        return false;
    }
    
    // Remove from specific storage
    switch (it->subtype) {
        case ConnectionSubtype::Disulfide:
            m_disulfides.remove(connectionId);
            break;
        case ConnectionSubtype::Interchain:
            m_interchainConnections.remove(connectionId);
            break;
        case ConnectionSubtype::Branch:
            m_branches.remove(connectionId);
            break;
        case ConnectionSubtype::Cyclization:
            m_cyclizations.remove(connectionId);
            break;
        default:
            break;
    }
    
    m_extendedConnections.erase(it);
    return true;
}

BiopolymerExtendedConnection* ExtendedConnectionManager::getExtendedConnection(int connectionId) {
    auto it = m_extendedConnections.find(connectionId);
    return (it != m_extendedConnections.end()) ? &it.value() : nullptr;
}

const BiopolymerExtendedConnection* ExtendedConnectionManager::getExtendedConnection(int connectionId) const {
    auto it = m_extendedConnections.find(connectionId);
    return (it != m_extendedConnections.end()) ? &it.value() : nullptr;
}

QList<BiopolymerExtendedConnection> ExtendedConnectionManager::allExtendedConnections() const {
    return m_extendedConnections.values();
}

QList<BiopolymerExtendedConnection> ExtendedConnectionManager::extendedConnectionsByType(
    ConnectionSubtype subtype) const {
    QList<BiopolymerExtendedConnection> result;
    for (const auto& conn : m_extendedConnections) {
        if (conn.subtype == subtype) {
            result.append(conn);
        }
    }
    return result;
}

// ========== Query Methods ==========

int ExtendedConnectionManager::disulfideCount() const {
    return m_disulfides.size();
}

int ExtendedConnectionManager::interchainConnectionCount() const {
    return m_interchainConnections.size();
}

int ExtendedConnectionManager::branchCount() const {
    return m_branches.size();
}

int ExtendedConnectionManager::cyclizationCount() const {
    return m_cyclizations.size();
}

int ExtendedConnectionManager::extendedConnectionCount() const {
    return m_extendedConnections.size();
}

bool ExtendedConnectionManager::hasExtendedConnections() const {
    return !m_extendedConnections.isEmpty();
}

bool ExtendedConnectionManager::hasDisulfides() const {
    return !m_disulfides.isEmpty();
}

bool ExtendedConnectionManager::hasInterchainConnections() const {
    return !m_interchainConnections.isEmpty();
}

bool ExtendedConnectionManager::hasBranches() const {
    return !m_branches.isEmpty();
}

bool ExtendedConnectionManager::hasCyclizations() const {
    return !m_cyclizations.isEmpty();
}

bool ExtendedConnectionManager::hasInterchainConnectionsBetweenChains() const {
    return !m_interchainConnections.isEmpty();
}

// ========== Validation ==========

QList<QString> ExtendedConnectionManager::validateAll() const {
    ExtendedConnectionValidator validator(m_state);
    return validator.validateAllExtendedConnections();
}

bool ExtendedConnectionManager::areAllDisulfidesValid() const {
    ExtendedConnectionValidator validator(m_state);
    return validator.validateAllDisulfides().isEmpty();
}

bool ExtendedConnectionManager::areAllInterchainConnectionsValid() const {
    ExtendedConnectionValidator validator(m_state);
    return validator.validateAllInterchainConnections().isEmpty();
}

// ========== Utility ==========

QList<BiopolymerExtendedConnection> ExtendedConnectionManager::connectionsForResidue(int residueId) const {
    QList<BiopolymerExtendedConnection> result;
    for (const auto& conn : m_extendedConnections) {
        if (conn.involvesResidue(residueId)) {
            result.append(conn);
        }
    }
    return result;
}

QList<BiopolymerExtendedConnection> ExtendedConnectionManager::connectionsForChain(int chainId) const {
    QList<BiopolymerExtendedConnection> result;
    for (const auto& conn : m_extendedConnections) {
        if (conn.involvesChain(chainId)) {
            result.append(conn);
        }
    }
    return result;
}

bool ExtendedConnectionManager::isResidueInExtendedConnection(int residueId) const {
    for (const auto& conn : m_extendedConnections) {
        if (conn.involvesResidue(residueId)) {
            return true;
        }
    }
    return false;
}

// ========== ID Management ==========

int ExtendedConnectionManager::allocateDisulfideId() {
    return m_nextDisulfideId++;
}

int ExtendedConnectionManager::allocateInterchainConnectionId() {
    return m_nextInterchainConnectionId++;
}

int ExtendedConnectionManager::allocateBranchId() {
    return m_nextBranchId++;
}

int ExtendedConnectionManager::allocateCyclizationId() {
    return m_nextCyclizationId++;
}

int ExtendedConnectionManager::allocateExtendedConnectionId() {
    return m_nextExtendedConnectionId++;
}

int ExtendedConnectionManager::nextDisulfideId() const {
    return m_nextDisulfideId;
}

int ExtendedConnectionManager::nextInterchainConnectionId() const {
    return m_nextInterchainConnectionId;
}

int ExtendedConnectionManager::nextBranchId() const {
    return m_nextBranchId;
}

int ExtendedConnectionManager::nextCyclizationId() const {
    return m_nextCyclizationId;
}

int ExtendedConnectionManager::nextExtendedConnectionId() const {
    return m_nextExtendedConnectionId;
}

