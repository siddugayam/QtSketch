// src/app/molecule/BiopolymerWorkspace.cpp
#include "BiopolymerWorkspace.h"
#include "DocumentState.h"
#include "EditableMolecule.h"
#include <cmath>
#include <algorithm>

// ========== BiopolymerSelectionEvent Implementation ==========

bool BiopolymerSelectionEvent::isEmpty() const {
    return atomIds.isEmpty() && bondIds.isEmpty() && 
           residueIds.isEmpty() && chainIds.isEmpty();
}

bool BiopolymerSelectionEvent::containsAtom(AtomId atomId) const {
    return atomIds.contains(atomId);
}

bool BiopolymerSelectionEvent::containsResidue(int residueId) const {
    return residueIds.contains(residueId);
}

bool BiopolymerSelectionEvent::containsChain(int chainId) const {
    return chainIds.contains(chainId);
}

QSet<AtomId> BiopolymerSelectionEvent::allAtoms() const {
    return atomIds;
}

QSet<int> BiopolymerSelectionEvent::allResidues() const {
    return residueIds;
}

QSet<int> BiopolymerSelectionEvent::allChains() const {
    return chainIds;
}

// ========== BiopolymerVisibilityState Implementation ==========

bool BiopolymerVisibilityState::isChainVisible(int chainId) const {
    auto it = chainVisible.find(chainId);
    if (it != chainVisible.end()) {
        return it.value();
    }
    return true; // Default to visible
}

bool BiopolymerVisibilityState::isResidueVisible(int residueId) const {
    auto it = residueVisible.find(residueId);
    if (it != residueVisible.end()) {
        return it.value();
    }
    return true; // Default to visible
}

void BiopolymerVisibilityState::setChainVisible(int chainId, bool visible) {
    chainVisible[chainId] = visible;
}

void BiopolymerVisibilityState::setResidueVisible(int residueId, bool visible) {
    residueVisible[residueId] = visible;
}

void BiopolymerVisibilityState::reset() {
    chainVisible.clear();
    residueVisible.clear();
    showSideChains = true;
    showBackbone = true;
    showLabels = true;
}

// ========== BiopolymerWorkspaceState Implementation ==========

QMap<QString, QVariant> BiopolymerWorkspaceState::toJson() const {
    QMap<QString, QVariant> json;
    
    // Serialize selection
    QMap<QString, QVariant> selectionJson;
    selectionJson["scope"] = static_cast<int>(selection.scope);
    selectionJson["isProgrammatic"] = selection.isProgrammatic;
    
    // Serialize atom IDs
    QList<QVariant> atomIdsList;
    for (AtomId id : selection.atomIds) {
        atomIdsList.append(id);
    }
    selectionJson["atomIds"] = atomIdsList;
    
    // Serialize bond IDs
    QList<QVariant> bondIdsList;
    for (BondId id : selection.bondIds) {
        bondIdsList.append(id);
    }
    selectionJson["bondIds"] = bondIdsList;
    
    // Serialize residue IDs
    QList<QVariant> residueIdsList;
    for (int id : selection.residueIds) {
        residueIdsList.append(id);
    }
    selectionJson["residueIds"] = residueIdsList;
    
    // Serialize chain IDs
    QList<QVariant> chainIdsList;
    for (int id : selection.chainIds) {
        chainIdsList.append(id);
    }
    selectionJson["chainIds"] = chainIdsList;
    
    json["selection"] = selectionJson;
    
    // Serialize visibility
    QMap<QString, QVariant> visibilityJson;
    
    QMap<QString, QVariant> chainVisibleJson;
    for (auto it = visibility.chainVisible.begin(); it != visibility.chainVisible.end(); ++it) {
        chainVisibleJson[QString::number(it.key())] = it.value();
    }
    visibilityJson["chainVisible"] = chainVisibleJson;
    
    QMap<QString, QVariant> residueVisibleJson;
    for (auto it = visibility.residueVisible.begin(); it != visibility.residueVisible.end(); ++it) {
        residueVisibleJson[QString::number(it.key())] = it.value();
    }
    visibilityJson["residueVisible"] = residueVisibleJson;
    
    visibilityJson["showSideChains"] = visibility.showSideChains;
    visibilityJson["showBackbone"] = visibility.showBackbone;
    visibilityJson["showLabels"] = visibility.showLabels;
    
    json["visibility"] = visibilityJson;
    
    // Serialize layout parameters
    QMap<QString, QVariant> layoutJson;
    layoutJson["strategy"] = static_cast<int>(layoutParams.strategy);
    layoutJson["bondLength"] = layoutParams.bondLength;
    layoutJson["bondAngle"] = layoutParams.bondAngle;
    layoutJson["helixRadius"] = layoutParams.helixRadius;
    layoutJson["helixPitch"] = layoutParams.helixPitch;
    layoutJson["helixStartAngle"] = layoutParams.helixStartAngle;
    layoutJson["sheetSpacing"] = layoutParams.sheetSpacing;
    layoutJson["circleRadius"] = layoutParams.circleRadius;
    layoutJson["circleStartAngle"] = layoutParams.circleStartAngle;
    layoutJson["margin"] = layoutParams.margin;
    layoutJson["autoCenter"] = layoutParams.autoCenter;
    
    json["layoutParams"] = layoutJson;
    
    // Serialize sync options
    QMap<QString, QVariant> syncJson;
    syncJson["syncSelection"] = syncOptions.syncSelection;
    syncJson["syncVisibility"] = syncOptions.syncVisibility;
    syncJson["syncLayout"] = syncOptions.syncLayout;
    syncJson["autoCenter"] = syncOptions.autoCenter;
    syncJson["maintainSelectionOnLayout"] = syncOptions.maintainSelectionOnLayout;
    syncJson["showSyncIndicators"] = syncOptions.showSyncIndicators;
    
    json["syncOptions"] = syncJson;
    
    return json;
}

BiopolymerWorkspaceState BiopolymerWorkspaceState::fromJson(const QMap<QString, QVariant>& json) {
    BiopolymerWorkspaceState state;
    
    // Deserialize selection
    QMap<QString, QVariant> selectionJson = json.value("selection").toMap();
    state.selection.scope = static_cast<SelectionScope>(selectionJson.value("scope").toInt());
    state.selection.isProgrammatic = selectionJson.value("isProgrammatic").toBool();
    
    for (const QVariant& id : selectionJson.value("atomIds").toList()) {
        state.selection.atomIds.insert(id.toInt());
    }
    for (const QVariant& id : selectionJson.value("bondIds").toList()) {
        state.selection.bondIds.insert(id.toInt());
    }
    for (const QVariant& id : selectionJson.value("residueIds").toList()) {
        state.selection.residueIds.insert(id.toInt());
    }
    for (const QVariant& id : selectionJson.value("chainIds").toList()) {
        state.selection.chainIds.insert(id.toInt());
    }
    
    // Deserialize visibility
    QMap<QString, QVariant> visibilityJson = json.value("visibility").toMap();
    
    QMap<QString, QVariant> chainVisibleJson = visibilityJson.value("chainVisible").toMap();
    for (auto it = chainVisibleJson.begin(); it != chainVisibleJson.end(); ++it) {
        state.visibility.chainVisible[it.key().toInt()] = it.value().toBool();
    }
    
    QMap<QString, QVariant> residueVisibleJson = visibilityJson.value("residueVisible").toMap();
    for (auto it = residueVisibleJson.begin(); it != residueVisibleJson.end(); ++it) {
        state.visibility.residueVisible[it.key().toInt()] = it.value().toBool();
    }
    
    state.visibility.showSideChains = visibilityJson.value("showSideChains").toBool();
    state.visibility.showBackbone = visibilityJson.value("showBackbone").toBool();
    state.visibility.showLabels = visibilityJson.value("showLabels").toBool();
    
    // Deserialize layout parameters
    QMap<QString, QVariant> layoutJson = json.value("layoutParams").toMap();
    state.layoutParams.strategy = static_cast<PlacementStrategy>(layoutJson.value("strategy").toInt());
    state.layoutParams.bondLength = layoutJson.value("bondLength").toDouble();
    state.layoutParams.bondAngle = layoutJson.value("bondAngle").toDouble();
    state.layoutParams.helixRadius = layoutJson.value("helixRadius").toDouble();
    state.layoutParams.helixPitch = layoutJson.value("helixPitch").toDouble();
    state.layoutParams.helixStartAngle = layoutJson.value("helixStartAngle").toDouble();
    state.layoutParams.sheetSpacing = layoutJson.value("sheetSpacing").toDouble();
    state.layoutParams.circleRadius = layoutJson.value("circleRadius").toDouble();
    state.layoutParams.circleStartAngle = layoutJson.value("circleStartAngle").toDouble();
    state.layoutParams.margin = layoutJson.value("margin").toDouble();
    state.layoutParams.autoCenter = layoutJson.value("autoCenter").toBool();
    
    // Deserialize sync options
    QMap<QString, QVariant> syncJson = json.value("syncOptions").toMap();
    state.syncOptions.syncSelection = syncJson.value("syncSelection").toBool();
    state.syncOptions.syncVisibility = syncJson.value("syncVisibility").toBool();
    state.syncOptions.syncLayout = syncJson.value("syncLayout").toBool();
    state.syncOptions.autoCenter = syncJson.value("autoCenter").toBool();
    state.syncOptions.maintainSelectionOnLayout = syncJson.value("maintainSelectionOnLayout").toBool();
    state.syncOptions.showSyncIndicators = syncJson.value("showSyncIndicators").toBool();
    
    return state;
}

// ========== BiopolymerWorkspace Implementation ==========

BiopolymerWorkspace::BiopolymerWorkspace(DocumentState& document, QObject* parent)
    : QObject(parent), m_document(document) {
    initialize();
}

BiopolymerWorkspace::~BiopolymerWorkspace() {
    unregisterAllCallbacks();
}

void BiopolymerWorkspace::initialize() {
    // Set default layout parameters
    m_layoutParams = LayoutParameters();
    
    // Set default sync options
    m_syncOptions = WorkspaceSyncOptions();
    
    // Reset state
    reset();
}

void BiopolymerWorkspace::reset() {
    // Clear selection
    m_selection = BiopolymerSelectionEvent();
    
    // Reset visibility
    m_visibility = BiopolymerVisibilityState();
    
    // Clear custom placements
    m_residuePlacements.clear();
    
    // Update views
    updateSequenceViewFromBiopolymer();
    
    emit dataChanged();
}

// ========== Synchronization Control ==========

void BiopolymerWorkspace::setSynchronizationEnabled(bool enabled) {
    m_syncEnabled = enabled;
}

bool BiopolymerWorkspace::isSynchronizationEnabled() const {
    return m_syncEnabled;
}

void BiopolymerWorkspace::setSyncOptions(const WorkspaceSyncOptions& options) {
    m_syncOptions = options;
}

WorkspaceSyncOptions BiopolymerWorkspace::syncOptions() const {
    return m_syncOptions;
}

void BiopolymerWorkspace::synchronize() {
    if (!m_syncEnabled) return;
    
    if (m_syncOptions.syncSelection) {
        syncSelection();
    }
    
    if (m_syncOptions.syncVisibility) {
        syncVisibility();
    }
    
    if (m_syncOptions.syncLayout) {
        applyLayout();
    }
    
    emit synchronizationComplete();
}

// ========== Selection Management ==========

void BiopolymerWorkspace::setSelection(const BiopolymerSelectionEvent& selection) {
    m_selection = selection;
    
    if (m_syncEnabled && m_syncOptions.syncSelection) {
        syncSelectionToChemistry();
    }
    
    emit selectionChanged(m_selection);
    notifyEvent(WorkspaceEvent{WorkspaceEventType::SelectionChanged, m_selection, m_visibility, m_layoutParams, "BiopolymerWorkspace"});
}

BiopolymerSelectionEvent BiopolymerWorkspace::selection() const {
    return m_selection;
}

void BiopolymerWorkspace::selectResidue(int residueId) {
    BiopolymerSelectionEvent selection;
    selection.scope = SelectionScope::Residue;
    selection.residueIds.insert(residueId);
    
    // Also select the atoms
    selection.atomIds = residueIdsToAtomIds({residueId});
    
    setSelection(selection);
}

void BiopolymerWorkspace::selectResidues(const QSet<int>& residueIds) {
    BiopolymerSelectionEvent selection;
    selection.scope = SelectionScope::MultipleResidues;
    selection.residueIds = residueIds;
    
    // Also select the atoms
    selection.atomIds = residueIdsToAtomIds(residueIds);
    
    setSelection(selection);
}

void BiopolymerWorkspace::selectChain(int chainId) {
    BiopolymerSelectionEvent selection;
    selection.scope = SelectionScope::Chain;
    selection.chainIds.insert(chainId);
    
    // Select all residues in the chain
    QList<int> residueIdsInChain = m_document.biopolymerState().residueIdsInChain(chainId);
    for (int residueId : residueIdsInChain) {
        selection.residueIds.insert(residueId);
    }
    
    // Also select the atoms
    selection.atomIds = residueIdsToAtomIds(selection.residueIds);
    
    setSelection(selection);
}

void BiopolymerWorkspace::selectAll() {
    BiopolymerSelectionEvent selection;
    selection.scope = SelectionScope::MultipleResidues;
    
    // Select all residues
    for (int residueId : m_document.biopolymerState().allResidueIds()) {
        selection.residueIds.insert(residueId);
    }
    
    // Also select all atoms
    selection.atomIds = residueIdsToAtomIds(selection.residueIds);
    
    setSelection(selection);
}

void BiopolymerWorkspace::clearSelection() {
    BiopolymerSelectionEvent selection;
    setSelection(selection);
}

void BiopolymerWorkspace::selectResiduesFromAtoms(const QSet<AtomId>& atomIds) {
    QSet<int> residueIds = atomIdsToResidueIds(atomIds);
    selectResidues(residueIds);
}

QSet<AtomId> BiopolymerWorkspace::selectAtomsFromResidues(const QSet<int>& residueIds) const {
    return residueIdsToAtomIds(residueIds);
}

// ========== Selection Synchronization ==========

void BiopolymerWorkspace::syncSelectionToChemistry() {
    if (!m_syncEnabled || !m_syncOptions.syncSelection) return;
    
    // Update document selection based on workspace selection
    DocumentState& doc = m_document;
    SelectionState& selection = doc.selection();
    
    // Clear current selection
    selection.clear();
    
    // Select atoms from residues
    for (AtomId atomId : m_selection.atomIds) {
        selection.addAtomToSelection(atomId);
    }
    
    // Also select bonds that are between selected atoms
    for (BondId bondId : doc.molecule().bondIds()) {
        AtomId a1 = -1, a2 = -1;
        if (doc.molecule().bondEndpoints(bondId, a1, a2)) {
            if (m_selection.atomIds.contains(a1) && m_selection.atomIds.contains(a2)) {
                selection.addBondToSelection(bondId);
            }
        }
    }
}

void BiopolymerWorkspace::syncSelectionFromChemistry() {
    if (!m_syncEnabled || !m_syncOptions.syncSelection) return;
    
    // Update workspace selection based on document selection
    const SelectionState& selection = m_document.selection();
    
    BiopolymerSelectionEvent event;
    
    // Get selected atoms from document
    for (AtomId atomId : selection.atomIds()) {
        event.atomIds.insert(atomId);
    }
    
    // Get selected bonds from document
    for (BondId bondId : selection.bondIds()) {
        event.bondIds.insert(bondId);
    }
    
    // Convert atoms to residues
    event.residueIds = atomIdsToResidueIds(event.atomIds);
    
    // Determine scope
    if (event.residueIds.size() == 1) {
        event.scope = SelectionScope::Residue;
    } else if (event.residueIds.size() > 1) {
        event.scope = SelectionScope::MultipleResidues;
    } else if (event.atomIds.size() == 1) {
        event.scope = SelectionScope::Atom;
    } else if (event.bondIds.size() == 1) {
        event.scope = SelectionScope::Bond;
    }
    
    // Only update if selection changed
    if (event != m_selection) {
        m_selection = event;
        emit selectionChanged(m_selection);
    }
}

void BiopolymerWorkspace::syncSelection() {
    // Determine which direction to sync based on which view was last changed
    // For now, sync from chemistry to biopolymer
    syncSelectionFromChemistry();
    // Then sync back to ensure consistency
    syncSelectionToChemistry();
}

// ========== Visibility Management ==========

void BiopolymerWorkspace::setVisibility(const BiopolymerVisibilityState& visibility) {
    m_visibility = visibility;
    
    if (m_syncEnabled && m_syncOptions.syncVisibility) {
        syncVisibilityToChemistry();
    }
    
    emit visibilityChanged(m_visibility);
    notifyEvent(WorkspaceEvent{WorkspaceEventType::VisibilityChanged, m_selection, m_visibility, m_layoutParams, "BiopolymerWorkspace"});
}

BiopolymerVisibilityState BiopolymerWorkspace::visibility() const {
    return m_visibility;
}

void BiopolymerWorkspace::setChainVisible(int chainId, bool visible) {
    m_visibility.setChainVisible(chainId, visible);
    
    if (m_syncEnabled && m_syncOptions.syncVisibility) {
        syncVisibilityToChemistry();
    }
    
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::setResidueVisible(int residueId, bool visible) {
    m_visibility.setResidueVisible(residueId, visible);
    
    if (m_syncEnabled && m_syncOptions.syncVisibility) {
        syncVisibilityToChemistry();
    }
    
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::setAllChainsVisible(bool visible) {
    for (int chainId : m_document.biopolymerState().chainIds()) {
        m_visibility.setChainVisible(chainId, visible);
    }
    
    if (m_syncEnabled && m_syncOptions.syncVisibility) {
        syncVisibilityToChemistry();
    }
    
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::setAllResiduesVisible(bool visible) {
    for (int residueId : m_document.biopolymerState().allResidueIds()) {
        m_visibility.setResidueVisible(residueId, visible);
    }
    
    if (m_syncEnabled && m_syncOptions.syncVisibility) {
        syncVisibilityToChemistry();
    }
    
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::toggleSideChains() {
    m_visibility.showSideChains = !m_visibility.showSideChains;
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::toggleBackbone() {
    m_visibility.showBackbone = !m_visibility.showBackbone;
    emit visibilityChanged(m_visibility);
}

void BiopolymerWorkspace::toggleLabels() {
    m_visibility.showLabels = !m_visibility.showLabels;
    emit visibilityChanged(m_visibility);
}

// ========== Visibility Synchronization ==========

void BiopolymerWorkspace::syncVisibilityToChemistry() {
    if (!m_syncEnabled || !m_syncOptions.syncVisibility) return;
    
    DocumentState& doc = m_document;
    const BiopolymerState& bioState = doc.biopolymerState();
    
    // For now, visibility sync is a placeholder
    // In a full implementation, we would hide/show atoms based on visibility state
    
    // This would involve:
    // 1. Hiding atoms of invisible residues
    // 2. Hiding bonds between invisible atoms
    // 3. Updating the render state
    
    // For this implementation, we'll just emit the signal
    // The actual visibility sync would be handled by the rendering layer
}

void BiopolymerWorkspace::syncVisibilityFromChemistry() {
    if (!m_syncEnabled || !m_syncOptions.syncVisibility) return;
    
    // For now, this is a no-op
    // In a full implementation, we would:
    // 1. Check which atoms are visible in the chemistry view
    // 2. Update the visibility state accordingly
    
    // This would require integration with the rendering system
}

void BiopolymerWorkspace::syncVisibility() {
    syncVisibilityFromChemistry();
    syncVisibilityToChemistry();
}

// ========== Layout Management ==========

void BiopolymerWorkspace::setLayoutParameters(const LayoutParameters& params) {
    m_layoutParams = params;
}

LayoutParameters BiopolymerWorkspace::layoutParameters() const {
    return m_layoutParams;
}

void BiopolymerWorkspace::setPlacementStrategy(PlacementStrategy strategy) {
    m_layoutParams.strategy = strategy;
}

PlacementStrategy BiopolymerWorkspace::placementStrategy() const {
    return m_layoutParams.strategy;
}

void BiopolymerWorkspace::applyLayout() {
    const BiopolymerState& bioState = m_document.biopolymerState();
    
    for (int chainId : bioState.chainIds()) {
        applyLayoutToChain(chainId);
    }
    
    emit layoutChanged();
}

void BiopolymerWorkspace::applyLayoutToChain(int chainId) {
    switch (m_layoutParams.strategy) {
        case PlacementStrategy::Snake:
            applySnakeLayoutToChain(chainId);
            break;
        case PlacementStrategy::Helix:
            applyHelixLayoutToChain(chainId);
            break;
        case PlacementStrategy::Sheet:
            applySheetLayoutToChain(chainId);
            break;
        case PlacementStrategy::Circle:
            applyCircleLayoutToChain(chainId);
            break;
        case PlacementStrategy::Custom:
            // Custom layout is handled by individual residue placements
            break;
        case PlacementStrategy::FromTemplate:
            applyTemplateLayout();
            break;
    }
    
    emit layoutChanged();
}

void BiopolymerWorkspace::layoutResidue(int residueId) {
    const BiopolymerResidue* residue = m_document.biopolymerState().residue(residueId);
    if (!residue) return;
    
    // If we have a custom placement, use it
    auto it = m_residuePlacements.find(residueId);
    if (it != m_residuePlacements.end()) {
        // Apply custom placement
        const ResiduePlacement& placement = it.value();
        
        // For now, we'll just update the atoms' positions
        // In a full implementation, we'd need to move the residue's atoms
        
        return;
    }
    
    // Otherwise, layout as part of the chain
    applyLayoutToChain(residue->chainId);
}

void BiopolymerWorkspace::layoutResidueRange(int chainId, int startPosition, int endPosition) {
    const BiopolymerChain* chain = m_document.biopolymerState().chain(chainId);
    if (!chain) return;
    
    // For now, just re-layout the entire chain
    applyLayoutToChain(chainId);
}

// ========== Residue Placement ==========

void BiopolymerWorkspace::setResiduePlacement(int residueId, const ResiduePlacement& placement) {
    m_residuePlacements[residueId] = placement;
}

ResiduePlacement BiopolymerWorkspace::residuePlacement(int residueId) const {
    auto it = m_residuePlacements.find(residueId);
    if (it != m_residuePlacements.end()) {
        return it.value();
    }
    return ResiduePlacement();
}

void BiopolymerWorkspace::resetResiduePlacement(int residueId) {
    m_residuePlacements.remove(residueId);
}

void BiopolymerWorkspace::setResiduePosition(int residueId, double x, double y) {
    ResiduePlacement placement = residuePlacement(residueId);
    placement.residueId = residueId;
    placement.position = QPointF(x, y);
    placement.hasCustomPlacement = true;
    m_residuePlacements[residueId] = placement;
}

void BiopolymerWorkspace::setResiduePosition3D(int residueId, double x, double y, double z) {
    ResiduePlacement placement = residuePlacement(residueId);
    placement.residueId = residueId;
    placement.position3D = QPointF(x, y); // Note: QPointF only has x,y - z would need extension
    placement.hasCustomPlacement = true;
    m_residuePlacements[residueId] = placement;
}

void BiopolymerWorkspace::setResidueOrientation(int residueId, const QQuaternion& orientation) {
    ResiduePlacement placement = residuePlacement(residueId);
    placement.residueId = residueId;
    placement.orientation = orientation;
    placement.hasCustomPlacement = true;
    m_residuePlacements[residueId] = placement;
}

// ========== Workspace State ==========

BiopolymerWorkspaceState BiopolymerWorkspace::state() const {
    BiopolymerWorkspaceState state;
    state.selection = m_selection;
    state.visibility = m_visibility;
    state.layoutParams = m_layoutParams;
    state.residuePlacements = m_residuePlacements;
    state.syncOptions = m_syncOptions;
    return state;
}

void BiopolymerWorkspace::restoreState(const BiopolymerWorkspaceState& state) {
    m_selection = state.selection;
    m_visibility = state.visibility;
    m_layoutParams = state.layoutParams;
    m_residuePlacements = state.residuePlacements;
    m_syncOptions = state.syncOptions;
    
    // Update views
    updateSequenceViewFromBiopolymer();
    
    emit dataChanged();
    emit selectionChanged(m_selection);
    emit visibilityChanged(m_visibility);
    emit layoutChanged();
}

QMap<QString, QVariant> BiopolymerWorkspace::toJson() const {
    return state().toJson();
}

void BiopolymerWorkspace::fromJson(const QMap<QString, QVariant>& json) {
    BiopolymerWorkspaceState state = BiopolymerWorkspaceState::fromJson(json);
    restoreState(state);
}

// ========== Event Handling ==========

void BiopolymerWorkspace::registerEventCallback(const WorkspaceEventCallback& callback) {
    m_eventCallbacks.append(callback);
}

void BiopolymerWorkspace::registerSelectionChangedCallback(const SelectionChangedCallback& callback) {
    m_selectionCallbacks.append(callback);
}

void BiopolymerWorkspace::registerVisibilityChangedCallback(const VisibilityChangedCallback& callback) {
    m_visibilityCallbacks.append(callback);
}

void BiopolymerWorkspace::unregisterAllCallbacks() {
    m_eventCallbacks.clear();
    m_selectionCallbacks.clear();
    m_visibilityCallbacks.clear();
}

void BiopolymerWorkspace::emitSelectionChanged(const BiopolymerSelectionEvent& selection) {
    for (const auto& callback : m_selectionCallbacks) {
        callback(selection);
    }
}

void BiopolymerWorkspace::emitVisibilityChanged(const BiopolymerVisibilityState& visibility) {
    for (const auto& callback : m_visibilityCallbacks) {
        callback(visibility);
    }
}

void BiopolymerWorkspace::notifyEvent(const WorkspaceEvent& event) {
    for (const auto& callback : m_eventCallbacks) {
        callback(event);
    }
}

// ========== Utility Methods ==========

bool BiopolymerWorkspace::hasBiopolymerData() const {
    return m_document.hasBiopolymerData();
}

DocumentState& BiopolymerWorkspace::document() {
    return m_document;
}

const DocumentState& BiopolymerWorkspace::document() const {
    return m_document;
}

BiopolymerState& BiopolymerWorkspace::biopolymerState() {
    return m_document.biopolymerState();
}

const BiopolymerState& BiopolymerWorkspace::biopolymerState() const {
    return m_document.biopolymerState();
}

EditableMolecule& BiopolymerWorkspace::molecule() {
    return m_document.molecule();
}

const EditableMolecule& BiopolymerWorkspace::molecule() const {
    return m_document.molecule();
}

BiopolymerSequenceView& BiopolymerWorkspace::sequenceView() {
    return m_document.bioView();
}

const BiopolymerSequenceView& BiopolymerWorkspace::sequenceView() const {
    return m_document.bioView();
}

void BiopolymerWorkspace::updateSequenceView() {
    updateSequenceViewFromBiopolymer();
}

QRectF BiopolymerWorkspace::selectionBoundingBox() const {
    if (m_selection.atomIds.isEmpty()) {
        return QRectF();
    }
    
    const EditableMolecule& mol = m_document.molecule();
    
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    
    for (AtomId atomId : m_selection.atomIds) {
        double x = 0, y = 0;
        if (mol.atomPos(atomId, x, y)) {
            minX = std::min(minX, x);
            minY = std::min(minY, y);
            maxX = std::max(maxX, x);
            maxY = std::max(maxY, y);
        }
    }
    
    if (minX > maxX || minY > maxY) {
        return QRectF();
    }
    
    return QRectF(minX, minY, maxX - minX, maxY - minY);
}

QRectF BiopolymerWorkspace::biopolymerBoundingBox() const {
    const EditableMolecule& mol = m_document.molecule();
    const BiopolymerState& bioState = m_document.biopolymerState();
    
    QSet<AtomId> allAtoms;
    for (int residueId : bioState.allResidueIds()) {
        QList<AtomId> resAtoms = bioState.atomsOwnedByResidue(residueId);
        for (AtomId atomId : resAtoms) {
            allAtoms.insert(atomId);
        }
    }
    
    if (allAtoms.isEmpty()) {
        return QRectF();
    }
    
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    
    for (AtomId atomId : allAtoms) {
        double x = 0, y = 0;
        if (mol.atomPos(atomId, x, y)) {
            minX = std::min(minX, x);
            minY = std::min(minY, y);
            maxX = std::max(maxX, x);
            maxY = std::max(maxY, y);
        }
    }
    
    if (minX > maxX || minY > maxY) {
        return QRectF();
    }
    
    return QRectF(minX, minY, maxX - minX, maxY - minY);
}

void BiopolymerWorkspace::centerOnSelection() {
    QRectF bbox = selectionBoundingBox();
    if (bbox.isValid()) {
        // In a full implementation, this would center the view on the bbox
        // For now, we just emit the signal
        emit layoutChanged();
    }
}

void BiopolymerWorkspace::centerOnBiopolymer() {
    QRectF bbox = biopolymerBoundingBox();
    if (bbox.isValid()) {
        // In a full implementation, this would center the view on the bbox
        emit layoutChanged();
    }
}

// ========== Helper Methods ==========

QSet<AtomId> BiopolymerWorkspace::residueIdsToAtomIds(const QSet<int>& residueIds) const {
    QSet<AtomId> atomIds;
    const BiopolymerState& bioState = m_document.biopolymerState();
    
    for (int residueId : residueIds) {
        QList<AtomId> resAtoms = bioState.atomsOwnedByResidue(residueId);
        for (AtomId atomId : resAtoms) {
            atomIds.insert(atomId);
        }
    }
    
    return atomIds;
}

QSet<int> BiopolymerWorkspace::atomIdsToResidueIds(const QSet<AtomId>& atomIds) const {
    QSet<int> residueIds;
    const BiopolymerState& bioState = m_document.biopolymerState();
    
    for (AtomId atomId : atomIds) {
        int residueId = bioState.residueOwningAtom(atomId);
        if (residueId != -1) {
            residueIds.insert(residueId);
        }
    }
    
    return residueIds;
}

// ========== Layout Strategies ==========

void BiopolymerWorkspace::applySnakeLayout() {
    m_layoutParams.strategy = PlacementStrategy::Snake;
    applyLayout();
}

void BiopolymerWorkspace::applyHelixLayout() {
    m_layoutParams.strategy = PlacementStrategy::Helix;
    applyLayout();
}

void BiopolymerWorkspace::applySheetLayout() {
    m_layoutParams.strategy = PlacementStrategy::Sheet;
    applyLayout();
}

void BiopolymerWorkspace::applyCircleLayout() {
    m_layoutParams.strategy = PlacementStrategy::Circle;
    applyLayout();
}

void BiopolymerWorkspace::applyCustomLayout(const QList<ResiduePlacement>& placements) {
    m_layoutParams.strategy = PlacementStrategy::Custom;
    
    for (const ResiduePlacement& placement : placements) {
        m_residuePlacements[placement.residueId] = placement;
    }
    
    emit layoutChanged();
}

void BiopolymerWorkspace::applyTemplateLayout() {
    m_layoutParams.strategy = PlacementStrategy::FromTemplate;
    
    // In a full implementation, this would use monomer template positions
    // For now, just apply a default layout
    applyLayout();
}

// ========== Snake Layout ==========

void BiopolymerWorkspace::applySnakeLayoutToChain(int chainId) {
    const BiopolymerChain* chain = m_document.biopolymerState().chain(chainId);
    if (!chain) return;
    
    const BiopolymerState& bioState = m_document.biopolymerState();
    EditableMolecule& mol = m_document.molecule();
    
    double x = 0.0;
    double y = 0.0;
    double angle = 0.0;
    
    // Convert bond angle to radians
    double bondAngleRad = m_layoutParams.bondAngle * M_PI / 180.0;
    
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = bioState.residue(residueId);
        if (!residue) continue;
        
        // Get atoms for this residue
        QList<AtomId> atomIds = bioState.atomsOwnedByResidue(residueId);
        
        if (atomIds.isEmpty()) {
            // Residue has no atoms - use a placeholder position
            x += m_layoutParams.bondLength * std::cos(angle);
            y += m_layoutParams.bondLength * std::sin(angle);
            continue;
        }
        
        // For simplicity, we'll just move the first atom (central atom)
        // In a full implementation, we'd move all atoms relative to the central atom
        AtomId centralAtom = atomIds.first();
        mol.setAtomPos(centralAtom, x, y);
        
        // Move other atoms relative to central atom
        // This is a simplified approach
        for (int i = 1; i < atomIds.size(); ++i) {
            AtomId atomId = atomIds[i];
            // Get current relative position
            double rx = 0, ry = 0;
            if (mol.atomPos(atomId, rx, ry)) {
                // For now, keep relative positions
                // In a full implementation, we'd use monomer template offsets
            }
        }
        
        // Move to next position
        x += m_layoutParams.bondLength * std::cos(angle);
        y += m_layoutParams.bondLength * std::sin(angle);
        
        // Alternate angle for snake layout
        angle = -angle;
    }
}

// ========== Helix Layout ==========

void BiopolymerWorkspace::applyHelixLayoutToChain(int chainId) {
    const BiopolymerChain* chain = m_document.biopolymerState().chain(chainId);
    if (!chain) return;
    
    const BiopolymerState& bioState = m_document.biopolymerState();
    EditableMolecule& mol = m_document.molecule();
    
    double centerX = 0.0;
    double centerY = 0.0;
    double currentAngle = m_layoutParams.helixStartAngle * M_PI / 180.0;
    double currentZ = 0.0;
    
    int residueCount = chain->residueIds.size();
    double angleIncrement = 2.0 * M_PI / std::max(1, residueCount);
    double zIncrement = m_layoutParams.helixPitch / std::max(1, residueCount);
    
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = bioState.residue(residueId);
        if (!residue) continue;
        
        QList<AtomId> atomIds = bioState.atomsOwnedByResidue(residueId);
        
        if (atomIds.isEmpty()) continue;
        
        // Calculate position in 3D (projected to 2D for now)
        double radius = m_layoutParams.helixRadius;
        double x = centerX + radius * std::cos(currentAngle);
        double y = centerY + radius * std::sin(currentAngle);
        
        // Move central atom
        AtomId centralAtom = atomIds.first();
        mol.setAtomPos(centralAtom, x, y);
        
        // Move other atoms relative to central atom
        for (int i = 1; i < atomIds.size(); ++i) {
            AtomId atomId = atomIds[i];
            // Keep relative positions for now
        }
        
        // Increment angle and Z
        currentAngle += angleIncrement;
        currentZ += zIncrement;
    }
}

// ========== Sheet Layout ==========

void BiopolymerWorkspace::applySheetLayoutToChain(int chainId) {
    const BiopolymerChain* chain = m_document.biopolymerState().chain(chainId);
    if (!chain) return;
    
    const BiopolymerState& bioState = m_document.biopolymerState();
    EditableMolecule& mol = m_document.molecule();
    
    double x = 0.0;
    double y = 0.0;
    
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = bioState.residue(residueId);
        if (!residue) continue;
        
        QList<AtomId> atomIds = bioState.atomsOwnedByResidue(residueId);
        
        if (atomIds.isEmpty()) continue;
        
        // Zig-zag pattern for sheet
        double offset = (residue->sequencePosition % 2 == 0) ? m_layoutParams.sheetSpacing : 0.0;
        
        // Move central atom
        AtomId centralAtom = atomIds.first();
        mol.setAtomPos(centralAtom, x, y + offset);
        
        x += m_layoutParams.bondLength;
    }
}

// ========== Circle Layout ==========

void BiopolymerWorkspace::applyCircleLayoutToChain(int chainId) {
    const BiopolymerChain* chain = m_document.biopolymerState().chain(chainId);
    if (!chain) return;
    
    const BiopolymerState& bioState = m_document.biopolymerState();
    EditableMolecule& mol = m_document.molecule();
    
    double centerX = 0.0;
    double centerY = 0.0;
    double currentAngle = m_layoutParams.circleStartAngle * M_PI / 180.0;
    
    int residueCount = chain->residueIds.size();
    double angleIncrement = 2.0 * M_PI / std::max(1, residueCount);
    
    for (int residueId : chain->residueIds) {
        const BiopolymerResidue* residue = bioState.residue(residueId);
        if (!residue) continue;
        
        QList<AtomId> atomIds = bioState.atomsOwnedByResidue(residueId);
        
        if (atomIds.isEmpty()) continue;
        
        // Calculate position on circle
        double x = centerX + m_layoutParams.circleRadius * std::cos(currentAngle);
        double y = centerY + m_layoutParams.circleRadius * std::sin(currentAngle);
        
        // Move central atom
        AtomId centralAtom = atomIds.first();
        mol.setAtomPos(centralAtom, x, y);
        
        currentAngle += angleIncrement;
    }
}

// ========== Sequence View Synchronization ==========

void BiopolymerWorkspace::updateSequenceViewFromBiopolymer() {
    const BiopolymerState& bioState = m_document.biopolymerState();
    BiopolymerSequenceView& view = m_document.bioView();
    
    // Build sequence text from biopolymer state
    QString sequence;
    
    for (int chainId : bioState.chainIds()) {
        const BiopolymerChain* chain = bioState.chain(chainId);
        if (!chain) continue;
        
        for (int residueId : chain->residueIds) {
            const BiopolymerResidue* residue = bioState.residue(residueId);
            if (residue) {
                sequence += residue->identity.symbol;
            }
        }
        
        // Add separator between chains
        sequence += "|";
    }
    
    // Remove trailing separator
    if (!sequence.isEmpty() && sequence.endsWith('|')) {
        sequence.chop(1);
    }
    
    // Determine sequence type
    QString seqType;
    if (!bioState.chainIds().isEmpty()) {
        const BiopolymerChain* firstChain = bioState.chain(bioState.chainIds().first());
        if (firstChain) {
            seqType = BiopolymerState::polymerTypeToString(firstChain->polymerType);
        }
    }
    
    if (seqType.isEmpty()) {
        seqType = "PEPTIDE";
    }
    
    // Build the sequence view
    view.buildSequenceView(sequence, seqType);
}

void BiopolymerWorkspace::updateBiopolymerFromSequenceView() {
    // This is intentionally limited
    // The sequence view is read-only, so we don't typically update from it
    // However, if there are changes in the sequence view, we could update the biopolymer state
    
    // For now, this is a no-op
}

// ========== Document Change Handling ==========

void BiopolymerWorkspace::onDocumentChanged() {
    // Update sequence view when document changes
    updateSequenceViewFromBiopolymer();
    
    // Re-synchronize selection and visibility
    if (m_syncEnabled) {
        if (m_syncOptions.syncSelection) {
            syncSelectionFromChemistry();
        }
        
        if (m_syncOptions.syncVisibility) {
            syncVisibilityFromChemistry();
        }
    }
    
    emit dataChanged();
}
