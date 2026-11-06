// tests/placement_engines_test.cpp
// Tests for PlacementEngines and PlacementCleanup
#include <cstdio>
#include <cmath>
#include <QList>
#include <QPointF>
#include "app/PlacementEngines.h"
#include "app/PlacementUtils.h"
#include "app/PlacementCleanup.h"
#include "app/molecule/EditableMolecule.h"

static int g_pass = 0, g_fail = 0;
#define CHECK(cond, name) do { \
    if (cond) { ++g_pass; std::printf("[PASS] %s\n", name); } \
    else { ++g_fail; std::printf("[FAIL] %s (line %d)\n", name, __LINE__); } \
} while (0)

static const double PI = 3.14159265358979323846;

static void test_snapAngle_noExistingAngles() {
    std::printf("--- Test: PlacementUtils::snapAngle with no existing angles ---\n");
    double snapped = PlacementUtils::snapAngle(0.1, {});
    CHECK(std::abs(snapped - 0.0) < 0.001, "snaps to nearest 30-degree increment (0)");

    double snapped2 = PlacementUtils::snapAngle(PI / 6.0 + 0.05, {});
    CHECK(std::abs(snapped2 - PI / 6.0) < 0.001, "snaps to nearest 30-degree increment (30deg)");
}

static void test_snapAngle_singleStericAvoidance() {
    std::printf("--- Test: PlacementUtils::snapAngle avoids one close existing angle ---\n");
    // Raw angle snaps to 0; an existing bond already sits at angle 0 -- must nudge away.
    double snapped = PlacementUtils::snapAngle(0.02, {0.0});
    CHECK(std::abs(snapped - 0.0) > 0.01, "nudges away from the single existing angle");
    CHECK(!PlacementUtils::isAngleTooClose(snapped, {0.0}), "nudged result clears the existing angle");
}

static void test_snapAngle_iterativeStericAvoidance() {
    std::printf("--- Test: PlacementUtils::snapAngle re-checks after nudging (iterative, not one-shot) ---\n");
    // Existing bonds at 0 and at +30deg (one snap interval away). A naive single-nudge
    // implementation would snap to 0, nudge once to +30deg, and stop there WITHOUT
    // noticing +30deg is ALSO occupied -- this is exactly the bug the consolidation
    // was supposed to fix (the old AtomPlacementEngine::snapAngle "simplistic approach").
    QList<double> existing = {0.0, PI / 6.0};
    double snapped = PlacementUtils::snapAngle(0.02, existing);
    CHECK(!PlacementUtils::isAngleTooClose(snapped, existing),
          "iterative avoidance clears ALL existing angles, not just the first nudge target");
}

static void test_snapAngle_givesUpGracefullyWhenSurrounded() {
    std::printf("--- Test: PlacementUtils::snapAngle falls back to best candidate when every nearby slot is occupied ---\n");
    // Occupy every 30-degree slot around the full circle -- no clear angle can exist.
    QList<double> existing;
    for (int i = 0; i < 12; ++i) existing.append(i * PI / 6.0);
    double snapped = PlacementUtils::snapAngle(0.02, existing);
    // Must not crash/hang and must return SOME finite angle (the MAX_NUDGES-bounded fallback).
    CHECK(std::isfinite(snapped), "returns a finite fallback angle instead of hanging or crashing");
}

static void test_cleanupPlacement_skipsInternalBonds() {
    std::printf("--- Test: PlacementCleanup::cleanupPlacement leaves template/ring-internal bonds alone ---\n");
    // Regression test: a bond between two NEWLY PLACED atoms must be left exactly as
    // the tool/template built it, even if its length differs from targetBondLength --
    // e.g. a template-internal bond of length 1.0 next to a document bondLength of 1.5.
    EditableMolecule mol;
    AtomId target = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId attach = mol.addAtom(QStringLiteral("C"), 1.5, 0.0);   // correct attach distance
    AtomId far = mol.addAtom(QStringLiteral("C"), 2.5, 0.0);      // 1.0 further out (template-internal)
    mol.addBond(target, attach, 1);
    mol.addBond(attach, far, 1);

    QList<AtomId> newlyPlaced = {attach, far};
    PlacementCleanup::cleanupPlacement(mol, newlyPlaced, 1.5);

    double fx = 0, fy = 0;
    mol.atomPos(far, fx, fy);
    double ax = 0, ay = 0;
    mol.atomPos(attach, ax, ay);
    double internalBondLength = std::hypot(fx - ax, fy - ay);
    CHECK(std::abs(internalBondLength - 1.0) < 0.001,
          "new-to-new (template-internal) bond length is untouched, stays 1.0");
}

static void test_cleanupPlacement_nudgesAttachmentBond() {
    std::printf("--- Test: PlacementCleanup::cleanupPlacement normalizes the attachment bond ---\n");
    // A bond between a newly placed atom and a PRE-EXISTING atom with the wrong
    // length must be corrected to targetBondLength.
    EditableMolecule mol;
    AtomId target = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId newAtom = mol.addAtom(QStringLiteral("C"), 3.0, 0.0);  // wrong length (should be 1.5)
    mol.addBond(target, newAtom, 1);

    QList<AtomId> newlyPlaced = {newAtom};
    // Single-atom placements are a documented no-op in cleanupPlacement; use two atoms
    // so the normalization path actually runs.
    AtomId newAtom2 = mol.addAtom(QStringLiteral("C"), 3.0, 1.0);
    mol.addBond(newAtom, newAtom2, 1);
    newlyPlaced.append(newAtom2);

    PlacementCleanup::cleanupPlacement(mol, newlyPlaced, 1.5);

    double tx = 0, ty = 0, nx = 0, ny = 0;
    mol.atomPos(target, tx, ty);
    mol.atomPos(newAtom, nx, ny);
    double attachLength = std::hypot(nx - tx, ny - ty);
    CHECK(std::abs(attachLength - 1.5) < 0.01, "attachment bond is nudged to targetBondLength");
}

static void test_cleanupPlacement_mergesNearCoincidentAtoms() {
    std::printf("--- Test: PlacementCleanup::cleanupPlacement merges near-coincident newly-placed atoms ---\n");
    EditableMolecule mol;
    AtomId a = mol.addAtom(QStringLiteral("C"), 0.0, 0.0);
    AtomId b = mol.addAtom(QStringLiteral("C"), 0.01, 0.0);  // effectively coincident with a
    mol.addBond(a, b, 1);

    int countBefore = mol.atomIds().size();
    QList<AtomId> newlyPlaced = {a, b};
    PlacementCleanup::cleanupPlacement(mol, newlyPlaced, 1.5);
    int countAfter = mol.atomIds().size();
    CHECK(countAfter < countBefore, "near-coincident atoms get merged, atom count drops");
}

int main() {
    std::printf("Running Placement Engines tests...\n");
    test_snapAngle_noExistingAngles();
    test_snapAngle_singleStericAvoidance();
    test_snapAngle_iterativeStericAvoidance();
    test_snapAngle_givesUpGracefullyWhenSurrounded();
    test_cleanupPlacement_skipsInternalBonds();
    test_cleanupPlacement_nudgesAttachmentBond();
    test_cleanupPlacement_mergesNearCoincidentAtoms();
    std::printf("Summary: %d passed, %d failed.\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
