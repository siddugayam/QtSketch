// src/app/PlacementCleanup.h
#pragma once
#include <QList>
#include <QSet>
#include "molecule/EditableMolecule.h"

// Post-placement cleanup utilities
// Used to fix up geometry issues after multi-atom placements

class DocumentState;

class PlacementCleanup {
public:
    // Constants for cleanup thresholds
    static constexpr double MIN_ATOM_SEPARATION = 0.3;      // Minimum distance between atoms (model units)
    static constexpr double BOND_LENGTH_TOLERANCE = 0.1;    // Tolerance for bond length deviation

    // Run cleanup on newly placed atoms
    // Only modifies atoms in newlyPlacedAtoms, not any pre-existing atoms
    // Returns the list of atoms that were actually moved
    static QList<AtomId> cleanupPlacement(EditableMolecule& mol, 
                                         const QList<AtomId>& newlyPlacedAtoms,
                                         double targetBondLength);

private:
    // Helper: find pairs of atoms that are too close together
    static QList<std::pair<AtomId, AtomId>> findNearCoincidentPairs(
        EditableMolecule& mol, const QList<AtomId>& newlyPlacedAtoms, double minSeparation);

    // Helper: nudge an atom to achieve target bond length
    static bool nudgeToBondLength(EditableMolecule& mol, AtomId atomToMove, AtomId fixedAtom,
                                  double targetLength, double tolerance);
};
