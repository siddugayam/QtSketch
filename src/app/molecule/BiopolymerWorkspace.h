// src/app/molecule/BiopolymerWorkspace.h
#ifndef BIOPOLYMERWORKSPACE_H
#define BIOPOLYMERWORKSPACE_H

// Biopolymer workspace management for synchronized sequence and chemistry views.
//
// Manages the relationship between:
// 1. BiopolymerState - The semantic biopolymer model (chains, residues, mappings)
// 2. EditableMolecule - The chemical structure (atoms, bonds)
// 3. BiopolymerSequenceView - The sequence display (read-only visualization)
//
// Provides:
// - Synchronized selection between sequence and chemistry views
// - Visibility synchronization (show/hide residues in both views)
// - Placement and layout strategies for imported biopolymers
// - Residue placement editing and 3D orientation support
// - Automatic synchronization when either view changes
//
// Design principle: The workspace maintains consistency between the semantic
// biopolymer model and the chemical structure. When one changes, the other is
// updated to match.

#include <QObject>
#include <QList>
#include <QSet>
#include <QHash>
#include <QMap>
#include <QPointF>
#include <QRectF>
#include <QVariant>
#include <QQuaternion>
#include "BiopolymerState.h"
#include "BiopolymerSequenceView.h"

// Forward declarations
class EditableMolecule;
class DocumentState;

// ========== Selection Types ==========

// Selection scope for synchronization
enum class SelectionScope {
    Atom,           // Single atom
    Bond,           // Single bond
    Residue,        // Single residue (and its atoms)
    Chain,          // Entire chain
    MultipleResidues, // Multiple residues
    Custom          // Custom selection
};

// Selection event data
struct BiopolymerSelectionEvent {
    SelectionScope scope = SelectionScope::Custom;
    
    // For atom selection
    QSet<AtomId> atomIds;
    
    // For bond selection
    QSet<BondId> bondIds;
    
    // For residue selection
    QSet<int> residueIds;
    
    // For chain selection
    QSet<int> chainIds;
    
    // Whether this is a programmatic selection (not user-initiated)
    bool isProgrammatic = false;
    
    // Check if selection is empty
    bool isEmpty() const;
    
    // Check if selection contains a specific atom
    bool containsAtom(AtomId atomId) const;
    
    // Check if selection contains a specific residue
    bool containsResidue(int residueId) const;
    
    // Check if selection contains a specific chain
    bool containsChain(int chainId) const;
    
    // Get all selected atoms
    QSet<AtomId> allAtoms() const;
    
    // Get all selected residues
    QSet<int> allResidues() const;
    
    // Get all selected chains
    QSet<int> allChains() const;
};

// ========== Visibility Types ==========

// Visibility state for residues/chains
struct BiopolymerVisibilityState {
    // Per-chain visibility
    QMap<int, bool> chainVisible;
    
    // Per-residue visibility (overrides chain)
    QMap<int, bool> residueVisible;
    
    // Whether to show side chains
    bool showSideChains = true;
    
    // Whether to show backbone
    bool showBackbone = true;
    
    // Whether to show labels
    bool showLabels = true;
    
    // Check if a chain is visible
    bool isChainVisible(int chainId) const;
    
    // Check if a residue is visible
    bool isResidueVisible(int residueId) const;
    
    // Set chain visibility
    void setChainVisible(int chainId, bool visible);
    
    // Set residue visibility
    void setResidueVisible(int residueId, bool visible);
    
    // Reset to default (all visible)
    void reset();
};

// ========== Placement Strategies ==========

// Placement strategy for imported biopolymers
enum class PlacementStrategy {
    Snake,           // Linear snake layout (default)
    Helix,           // Alpha-helix layout
    Sheet,           // Beta-sheet layout
    Circle,          // Circular layout
    Custom,          // Custom user-defined layout
    FromTemplate     // Use template positions
};

// Layout parameters
struct LayoutParameters {
    PlacementStrategy strategy = PlacementStrategy::Snake;
    
    // For snake layout
    double bondLength = 1.5;
    double bondAngle = 120.0; // degrees
    
    // For helix layout
    double helixRadius = 2.0;
    double helixPitch = 3.6; // rise per turn
    double helixStartAngle = 0.0;
    
    // For sheet layout
    double sheetSpacing = 5.0;
    
    // For circle layout
    double circleRadius = 10.0;
    double circleStartAngle = 0.0;
    
    // General parameters
    double margin = 1.0;
    bool autoCenter = true;
};

// Residue placement information
struct ResiduePlacement {
    int residueId = -1;
    QPointF position;          // 2D position in the view
    QPointF position3D;        // 3D position (if available)
    QQuaternion orientation;   // 3D orientation
    bool hasCustomPlacement = false;
    
    // Check if placement is valid
    bool isValid() const { return residueId != -1; }
};

// ========== Synchronization Options ==========

// Options for workspace synchronization
struct WorkspaceSyncOptions {
    // Whether to synchronize selection
    bool syncSelection = true;
    
    // Whether to synchronize visibility
    bool syncVisibility = true;
    
    // Whether to synchronize layout
    bool syncLayout = true;
    
    // Whether to auto-center on changes
    bool autoCenter = true;
    
    // Whether to maintain selection during layout changes
    bool maintainSelectionOnLayout = true;
    
    // Whether to show synchronization indicators
    bool showSyncIndicators = false;
};

// ========== Workspace State ==========

// Complete workspace state for serialization
struct BiopolymerWorkspaceState {
    // Selection state
    BiopolymerSelectionEvent selection;
    
    // Visibility state
    BiopolymerVisibilityState visibility;
    
    // Layout parameters
    LayoutParameters layoutParams;
    
    // Placement overrides
    QMap<int, ResiduePlacement> residuePlacements;
    
    // Synchronization options
    WorkspaceSyncOptions syncOptions;
    
    // Convert to JSON for serialization
    QMap<QString, QVariant> toJson() const;
    
    // Load from JSON
    static BiopolymerWorkspaceState fromJson(const QMap<QString, QVariant>& json);
};

// ========== Synchronization Event Types ==========

// Event types for workspace changes
enum class WorkspaceEventType {
    SelectionChanged,    // Selection changed in one view
    VisibilityChanged,   // Visibility changed
    LayoutChanged,       // Layout changed
    DataChanged,         // Underlying data changed (import, edit, etc.)
    SyncRequested,       // Manual sync requested
    ResetRequested       // Reset to default state
};

// Workspace event data
struct WorkspaceEvent {
    WorkspaceEventType type = WorkspaceEventType::DataChanged;
    BiopolymerSelectionEvent selection;
    BiopolymerVisibilityState visibility;
    LayoutParameters layoutParams;
    
    // Source of the event
    QString source;
    
    // Timestamp
    quint64 timestamp = 0;
};

// ========== Main Workspace Class ==========

// Callback types for workspace events
using WorkspaceEventCallback = std::function<void(const WorkspaceEvent&)>;
using SelectionChangedCallback = std::function<void(const BiopolymerSelectionEvent&)>;
using VisibilityChangedCallback = std::function<void(const BiopolymerVisibilityState&)>;

// Main biopolymer workspace class
class BiopolymerWorkspace : public QObject {
    Q_OBJECT

public:
    // Constructor
    BiopolymerWorkspace(DocumentState& document, QObject* parent = nullptr);
    ~BiopolymerWorkspace();
    
    // ========== Initialization ==========
    
    // Initialize the workspace
    void initialize();
    
    // Reset the workspace to default state
    void reset();
    
    // ========== Synchronization Control ==========
    
    // Enable/disable synchronization
    void setSynchronizationEnabled(bool enabled);
    bool isSynchronizationEnabled() const;
    
    // Set synchronization options
    void setSyncOptions(const WorkspaceSyncOptions& options);
    WorkspaceSyncOptions syncOptions() const;
    
    // Force a full synchronization
    void synchronize();
    
    // ========== Selection Management ==========
    
    // Set selection in the workspace
    void setSelection(const BiopolymerSelectionEvent& selection);
    
    // Get current selection
    BiopolymerSelectionEvent selection() const;
    
    // Select a residue
    void selectResidue(int residueId);
    
    // Select multiple residues
    void selectResidues(const QSet<int>& residueIds);
    
    // Select a chain
    void selectChain(int chainId);
    
    // Select all residues
    void selectAll();
    
    // Clear selection
    void clearSelection();
    
    // Select residues corresponding to atom selection
    void selectResiduesFromAtoms(const QSet<AtomId>& atomIds);
    
    // Select atoms corresponding to residue selection
    QSet<AtomId> selectAtomsFromResidues(const QSet<int>& residueIds) const;
    
    // ========== Selection Synchronization ==========
    
    // Synchronize selection from biopolymer to chemistry
    void syncSelectionToChemistry();
    
    // Synchronize selection from chemistry to biopolymer
    void syncSelectionFromChemistry();
    
    // Synchronize selection bidirectionally
    void syncSelection();
    
    // ========== Visibility Management ==========
    
    // Set visibility state
    void setVisibility(const BiopolymerVisibilityState& visibility);
    
    // Get current visibility state
    BiopolymerVisibilityState visibility() const;
    
    // Show/hide a chain
    void setChainVisible(int chainId, bool visible);
    
    // Show/hide a residue
    void setResidueVisible(int residueId, bool visible);
    
    // Show/hide all chains
    void setAllChainsVisible(bool visible);
    
    // Show/hide all residues
    void setAllResiduesVisible(bool visible);
    
    // Toggle side chain visibility
    void toggleSideChains();
    
    // Toggle backbone visibility
    void toggleBackbone();
    
    // Toggle labels
    void toggleLabels();
    
    // ========== Visibility Synchronization ==========
    
    // Synchronize visibility from biopolymer to chemistry
    void syncVisibilityToChemistry();
    
    // Synchronize visibility from chemistry to biopolymer
    void syncVisibilityFromChemistry();
    
    // Synchronize visibility bidirectionally
    void syncVisibility();
    
    // ========== Layout Management ==========
    
    // Set layout parameters
    void setLayoutParameters(const LayoutParameters& params);
    
    // Get current layout parameters
    LayoutParameters layoutParameters() const;
    
    // Set placement strategy
    void setPlacementStrategy(PlacementStrategy strategy);
    
    // Get current placement strategy
    PlacementStrategy placementStrategy() const;
    
    // Apply layout to the current biopolymer
    void applyLayout();
    
    // Apply layout to a specific chain
    void applyLayoutToChain(int chainId);
    
    // Layout a single residue
    void layoutResidue(int residueId);
    
    // Layout a range of residues
    void layoutResidueRange(int chainId, int startPosition, int endPosition);
    
    // ========== Residue Placement ==========
    
    // Set custom placement for a residue
    void setResiduePlacement(int residueId, const ResiduePlacement& placement);
    
    // Get residue placement
    ResiduePlacement residuePlacement(int residueId) const;
    
    // Reset residue placement to default
    void resetResiduePlacement(int residueId);
    
    // Set 2D position for a residue
    void setResiduePosition(int residueId, double x, double y);
    
    // Set 3D position for a residue
    void setResiduePosition3D(int residueId, double x, double y, double z);
    
    // Set 3D orientation for a residue
    void setResidueOrientation(int residueId, const QQuaternion& orientation);
    
    // ========== Workspace State ==========
    
    // Get complete workspace state
    BiopolymerWorkspaceState state() const;
    
    // Restore workspace state
    void restoreState(const BiopolymerWorkspaceState& state);
    
    // Save state to JSON
    QMap<QString, QVariant> toJson() const;
    
    // Load state from JSON
    void fromJson(const QMap<QString, QVariant>& json);
    
    // ========== Event Handling ==========
    
    // Register a callback for workspace events
    void registerEventCallback(const WorkspaceEventCallback& callback);
    
    // Register a callback for selection changes
    void registerSelectionChangedCallback(const SelectionChangedCallback& callback);
    
    // Register a callback for visibility changes
    void registerVisibilityChangedCallback(const VisibilityChangedCallback& callback);
    
    // Unregister all callbacks
    void unregisterAllCallbacks();
    
    // ========== Utility Methods ==========
    
    // Check if workspace has biopolymer data
    bool hasBiopolymerData() const;
    
    // Get the document
    DocumentState& document();
    const DocumentState& document() const;
    
    // Get the biopolymer state
    BiopolymerState& biopolymerState();
    const BiopolymerState& biopolymerState() const;
    
    // Get the molecule
    EditableMolecule& molecule();
    const EditableMolecule& molecule() const;
    
    // Get the sequence view
    BiopolymerSequenceView& sequenceView();
    const BiopolymerSequenceView& sequenceView() const;
    
    // Update sequence view from biopolymer state
    void updateSequenceView();
    
    // Get bounding box of selected residues
    QRectF selectionBoundingBox() const;
    
    // Get bounding box of all biopolymer atoms
    QRectF biopolymerBoundingBox() const;
    
    // Center view on selection
    void centerOnSelection();
    
    // Center view on biopolymer
    void centerOnBiopolymer();
    
    // ========== Layout Strategies ==========
    
    // Apply snake layout
    void applySnakeLayout();
    
    // Apply helix layout
    void applyHelixLayout();
    
    // Apply sheet layout
    void applySheetLayout();
    
    // Apply circle layout
    void applyCircleLayout();
    
    // Apply custom layout
    void applyCustomLayout(const QList<ResiduePlacement>& placements);
    
    // Apply layout from template positions
    void applyTemplateLayout();
    
signals:
    // Emitted when selection changes
    void selectionChanged(const BiopolymerSelectionEvent& selection);
    
    // Emitted when visibility changes
    void visibilityChanged(const BiopolymerVisibilityState& visibility);
    
    // Emitted when layout changes
    void layoutChanged();
    
    // Emitted when workspace data changes
    void dataChanged();
    
    // Emitted when synchronization completes
    void synchronizationComplete();
    
private slots:
    // Handle document changes
    void onDocumentChanged();
    
private:
    // Reference to the document
    DocumentState& m_document;
    
    // Current selection
    BiopolymerSelectionEvent m_selection;
    
    // Current visibility state
    BiopolymerVisibilityState m_visibility;
    
    // Current layout parameters
    LayoutParameters m_layoutParams;
    
    // Custom residue placements
    QMap<int, ResiduePlacement> m_residuePlacements;
    
    // Synchronization options
    WorkspaceSyncOptions m_syncOptions;
    
    // Whether synchronization is enabled
    bool m_syncEnabled = true;
    
    // Event callbacks
    QList<WorkspaceEventCallback> m_eventCallbacks;
    QList<SelectionChangedCallback> m_selectionCallbacks;
    QList<VisibilityChangedCallback> m_visibilityCallbacks;
    
    // Helper methods
    
    // Update sequence view from biopolymer state
    void updateSequenceViewFromBiopolymer();
    
    // Update biopolymer state from sequence view (limited)
    void updateBiopolymerFromSequenceView();
    
    // Convert residue IDs to atom IDs
    QSet<AtomId> residueIdsToAtomIds(const QSet<int>& residueIds) const;
    
    // Convert atom IDs to residue IDs
    QSet<int> atomIdsToResidueIds(const QSet<AtomId>& atomIds) const;
    
    // Apply snake layout to a chain
    void applySnakeLayoutToChain(int chainId);
    
    // Apply helix layout to a chain
    void applyHelixLayoutToChain(int chainId);
    
    // Apply sheet layout to a chain
    void applySheetLayoutToChain(int chainId);
    
    // Apply circle layout to a chain
    void applyCircleLayoutToChain(int chainId);
    
    // Emit selection changed signal
    void emitSelectionChanged(const BiopolymerSelectionEvent& selection);
    
    // Emit visibility changed signal
    void emitVisibilityChanged(const BiopolymerVisibilityState& visibility);
    
    // Notify all event callbacks
    void notifyEvent(const WorkspaceEvent& event);
};

#endif // BIOPOLYMERWORKSPACE_H
