// src/app/PlacementCleanup.cpp
#include "PlacementCleanup.h"
#include <cmath>
#include <algorithm>
#include <QDebug>

// Find pairs of atoms that are too close together
QList<std::pair<AtomId, AtomId>> PlacementCleanup::findNearCoincidentPairs(
    EditableMolecule& mol, const QList<AtomId>& newlyPlacedAtoms, double minSeparation) {
    QList<std::pair<AtomId, AtomId>> closePairs;
    QSet<AtomId> newAtomSet(newlyPlacedAtoms.begin(), newlyPlacedAtoms.end());

    // Only check pairs among newly placed atoms
    for (int i = 0; i < newlyPlacedAtoms.size(); ++i) {
        AtomId id1 = newlyPlacedAtoms[i];
        double x1 = 0, y1 = 0;
        if (!mol.atomPos(id1, x1, y1)) continue;

        for (int j = i + 1; j < newlyPlacedAtoms.size(); ++j) {
            AtomId id2 = newlyPlacedAtoms[j];
            double x2 = 0, y2 = 0;
            if (!mol.atomPos(id2, x2, y2)) continue;

            double dx = x2 - x1;
            double dy = y2 - y1;
            double dist = std::sqrt(dx * dx + dy * dy);

            if (dist < minSeparation) {
                closePairs.append({id1, id2});
            }
        }
    }

    return closePairs;
}

// Nudge an atom to achieve target bond length
// atomToMove: the atom to move
// fixedAtom: the atom to stay fixed
// targetLength: the desired bond length
// tolerance: acceptable deviation
// Returns true if the atom was moved
bool PlacementCleanup::nudgeToBondLength(EditableMolecule& mol, AtomId atomToMove, AtomId fixedAtom,
                                        double targetLength, double tolerance) {
    double fx = 0, fy = 0;
    if (!mol.atomPos(fixedAtom, fx, fy)) return false;

    double mx = 0, my = 0;
    if (!mol.atomPos(atomToMove, mx, my)) return false;

    double dx = mx - fx;
    double dy = my - fy;
    double currentLength = std::sqrt(dx * dx + dy * dy);

    if (currentLength == 0) {
        // Atoms are coincident, move atomToMove to a default position
        mol.setAtomPos(atomToMove, fx + targetLength, fy);
        return true;
    }

    if (std::abs(currentLength - targetLength) < tolerance) {
        return false;  // Already close enough
    }

    // Scale the vector to the target length
    double scale = targetLength / currentLength;
    double newX = fx + dx * scale;
    double newY = fy + dy * scale;

    mol.setAtomPos(atomToMove, newX, newY);
    return true;
}

// Run cleanup on newly placed atoms
// Only modifies atoms in newlyPlacedAtoms, not any pre-existing atoms
QList<AtomId> PlacementCleanup::cleanupPlacement(EditableMolecule& mol, 
                                                const QList<AtomId>& newlyPlacedAtoms,
                                                double targetBondLength) {
    QList<AtomId> movedAtoms;
    if (newlyPlacedAtoms.size() <= 1) {
        // Single atom placements don't need cleanup
        return movedAtoms;
    }

    // Capture pre-existing atoms before merge to check invariants
    QSet<AtomId> preExistingAtoms;
    for (AtomId id : mol.atomIds()) {
        if (!newlyPlacedAtoms.contains(id)) {
            preExistingAtoms.insert(id);
        }
    }

    // Step 1: Fix near-coincident atoms, scoped to newly placed atoms so
    // pre-existing atoms in the document are never merged with each other
    EditableMolecule::MergeResult mergeResult = mol.mergeOverlappingAtoms(0.1, newlyPlacedAtoms);

    // Chemistry-invariant verification: assert/warn if any pre-existing atom was merged away
    for (auto it = mergeResult.mergedAway.constBegin(); it != mergeResult.mergedAway.constEnd(); ++it) {
        AtomId doomed = it.key();
        if (preExistingAtoms.contains(doomed)) {
            qWarning() << "PlacementCleanup: chemistry invariant violated, pre-existing atom"
                       << doomed << "was merged away";
        }
    }

    // Filter newlyPlacedAtoms to surviving atoms (excluding those merged away)
    QList<AtomId> survivingNewAtoms;
    for (AtomId id : newlyPlacedAtoms) {
        if (!mergeResult.mergedAway.contains(id)) {
            survivingNewAtoms.append(id);
        }
    }

    // Step 2: Fix bond lengths for surviving newly placed atoms
    for (AtomId newAtom : survivingNewAtoms) {
        // Get all bonds connected to this atom
        QList<BondId> atomBonds = mol.bondIds();
        for (BondId bid : atomBonds) {
            AtomId a1 = -1, a2 = -1;
            if (!mol.bondEndpoints(bid, a1, a2)) continue;

            // Check if this bond is connected to the new atom
            if (a1 != newAtom && a2 != newAtom) continue;

            // Get the other atom (the one that's not newAtom)
            AtomId otherAtom = (a1 == newAtom) ? a2 : a1;

            // Only normalize the attachment bond -- the single connection between
            // the newly-placed structure and the pre-existing document. A bond
            // where BOTH endpoints are newly placed belongs to the tool/template
            // that built it (ring/chain geometry already uses targetBondLength by
            // construction; a template graft's internal bonds are deliberately its
            // OWN lengths, e.g. a 1.0-unit template bond next to a 1.5-unit
            // document bondLength) and must not be rewritten here.
            if (!preExistingAtoms.contains(otherAtom)) {
                continue;
            }

            // Nudge to target bond length
            if (nudgeToBondLength(mol, newAtom, otherAtom, targetBondLength, BOND_LENGTH_TOLERANCE)) {
                if (!movedAtoms.contains(newAtom)) {
                    movedAtoms.append(newAtom);
                }
            }
        }
    }

    return movedAtoms;
}
