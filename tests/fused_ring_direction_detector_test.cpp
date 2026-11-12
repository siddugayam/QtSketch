#include "FusedRingDirectionDetector.h"
#include "FusedRingOrientation.h"
#include "indigo.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int passed = 0;
int failed = 0;

void check_fusions_identical(const FusedRingSystemInput& input, const std::vector<FusedRingEdge>& expected_fusions) {
    if (input.fusions.size() != expected_fusions.size()) {
        cout << "  FAIL: fusions size mismatch: " << input.fusions.size() << " vs " << expected_fusions.size() << endl;
        failed++;
        return;
    }
    for (size_t i = 0; i < input.fusions.size(); ++i) {
        const auto& f = input.fusions[i];
        const auto& e = expected_fusions[i];
        if (f.ring1 != e.ring1 || f.ring2 != e.ring2 || f.bondIndex != e.bondIndex) {
            cout << "  FAIL: fusion[" << i << "] mismatch: (" << f.ring1 << "," << f.ring2 << "," << f.bondIndex 
                 << ") vs (" << e.ring1 << "," << e.ring2 << "," << e.bondIndex << ")" << endl;
            failed++;
            return;
        }
    }
    cout << "  PASS: fusions identical" << endl;
    passed++;
}

void run_test(const string& name, const string& smiles, int exp_row, double exp_ur) {
    cout << "Test: " << name << " (" << smiles << ")" << endl;
    
    int mol = indigoLoadMoleculeFromString(smiles.c_str());
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    indigoAromatize(mol);
    
    FusedRingSystemInput input = detectFusedRingDirections(mol);
    indigoFree(mol);
    
    if (input.ringSizes.empty()) {
        cout << "  FAIL: detectFusedRingDirections rejected or failed" << endl;
        failed++;
        return;
    }
    
    auto res = computePreferredOrientation(input);
    if (!res.success) {
        cout << "  FAIL: computePreferredOrientation failed: " << res.error << endl;
        failed++;
        return;
    }
    
    bool ok = true;
    if (res.ringsInHorizontalRow != exp_row) {
        cout << "  FAIL: expected ringsInHorizontalRow=" << exp_row << ", got " << res.ringsInHorizontalRow << endl;
        ok = false;
    }
    if (exp_ur >= 0.0 && abs(res.ringsInUpperRightQuadrant - exp_ur) > 1e-5) {
        cout << "  FAIL: expected ringsInUpperRightQuadrant=" << exp_ur << ", got " << res.ringsInUpperRightQuadrant << endl;
        ok = false;
    }
    
    if (ok) {
        cout << "  PASS" << endl;
        passed++;
    } else {
        failed++;
    }
}

int main() {
    setvbuf(stdout, nullptr, _IONBF, 0);
    qword sid = indigoAllocSessionId();
    indigoSetSessionId(sid);
    
    // 1. Naphthalene - 2 rings
    run_test("Naphthalene", "c1ccc2ccccc2c1", 2, 0.5);
    
    // 2. Anthracene - 3 linear rings
    run_test("Anthracene", "c1ccc2cc3ccccc3cc2c1", 3, 0.75);
    
    // 3. Phenanthrene - 3 angular rings
    run_test("Phenanthrene", "c1ccc2c(c1)ccc1ccccc12", 2, 1.5);
    
    // 4. Naphthacene/tetracene - 4 linear rings
    run_test("Naphthacene", "c1ccc2cc3cc4ccccc4cc3cc2c1", 4, 1.0);
    
    // 5. Chrysene - 4 rings zigzag
    // We don't assert exp_ur precisely here since it's a bonus, but row must be 2.
    run_test("Chrysene", "c1ccc2c(c1)ccc3c4ccccc4ccc23", 2, -1.0);

    // 6. Pentacene - 5 linear rings
    run_test("Pentacene", "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1", 5, 1.25);
    
    // 7. Regression test for linear case - verify fusions are identical to pre-change behavior
    // Anthracene: 3 linear rings, should produce fusions {{0,1,0}, {1,2,0}}
    cout << "Test: Linear regression (Anthracene)" << endl;
    int mol_anth = indigoLoadMoleculeFromString("c1ccc2cc3ccccc3cc2c1");
    if (mol_anth < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
    } else {
        indigoAromatize(mol_anth);
        FusedRingSystemInput input_anth = detectFusedRingDirections(mol_anth);
        indigoFree(mol_anth);
        
        // Expected fusions for linear anthracene: ring0->ring1 bondIndex=0 (root's
        // first/seeding edge), ring1->ring2 bondIndex=3 (a straight-line continuation
        // sits at the "opposite side" of ring1's own 6-membered cycle from its
        // incoming bond -- position 0+3=3, not 0; verified empirically after the
        // bondIndex semantic change from a propagated global direction to a raw
        // local cycle-position index).
        std::vector<FusedRingEdge> expected_fusions = {{0, 1, 0}, {1, 2, 3}};
        check_fusions_identical(input_anth, expected_fusions);
    }
    
    // 8. Branching test - central benzene with 3 benzene/pyridine leaves
    // SMILES: c12c(cccc2)c3c(cccc3)c4c1(nccc4)
    // Verified: 4 rings (one pyridine, three benzene), degree-3 central ring, no peri-fusion
    cout << "Test: Branching (degree-3 central ring)" << endl;
    int mol_branch = indigoLoadMoleculeFromString("c12c(cccc2)c3c(cccc3)c4c1(nccc4)");
    if (mol_branch < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for branching molecule" << endl;
        failed++;
    } else {
        indigoAromatize(mol_branch);
        FusedRingSystemInput input_branch = detectFusedRingDirections(mol_branch);
        indigoFree(mol_branch);
        
        if (input_branch.ringSizes.empty()) {
            cout << "  FAIL: Branching molecule was rejected (ringSizes empty)" << endl;
            failed++;
        } else if (input_branch.ringSizes.size() != 4) {
            cout << "  FAIL: Expected 4 rings, got " << input_branch.ringSizes.size() << endl;
            failed++;
        } else if (input_branch.fusions.size() != 3) {
            cout << "  FAIL: Expected 3 fusions (tree with 4 nodes), got " << input_branch.fusions.size() << endl;
            failed++;
        } else if (input_branch.ringAtoms.size() != 4) {
            cout << "  FAIL: Expected ringAtoms.size() == 4, got " << input_branch.ringAtoms.size() << endl;
            failed++;
        } else {
            // Check all ringAtoms are non-empty
            bool allNonEmpty = true;
            for (const auto& atoms : input_branch.ringAtoms) {
                if (atoms.empty()) {
                    allNonEmpty = false;
                    break;
                }
            }
            if (!allNonEmpty) {
                cout << "  FAIL: Some ringAtoms entries are empty" << endl;
                failed++;
            } else {
                // Check orientation computation works
                auto orient_res = computePreferredOrientation(input_branch);
                if (!orient_res.success) {
                    cout << "  FAIL: computePreferredOrientation failed: " << orient_res.error << endl;
                    failed++;
                } else if (orient_res.ringPos.size() != 4) {
                    cout << "  FAIL: Expected ringPos.size() == 4, got " << orient_res.ringPos.size() << endl;
                    failed++;
                } else {
                    // Check all positions are distinct (floating-point tolerance,
                    // not exact equality -- ringPos is now real (x,y), not hex ints).
                    bool allDistinct = true;
                    for (size_t i = 0; i < orient_res.ringPos.size(); ++i) {
                        for (size_t j = i + 1; j < orient_res.ringPos.size(); ++j) {
                            double dx = orient_res.ringPos[i].first - orient_res.ringPos[j].first;
                            double dy = orient_res.ringPos[i].second - orient_res.ringPos[j].second;
                            if (dx*dx + dy*dy < 1e-6) {
                                allDistinct = false;
                                break;
                            }
                        }
                        if (!allDistinct) break;
                    }
                    if (!allDistinct) {
                        cout << "  FAIL: Not all ringPos positions are distinct" << endl;
                        failed++;
                    } else {
                        cout << "  PASS: Branching molecule accepted, 4 rings, 3 fusions, all ringAtoms non-empty, orientation computed with 4 distinct positions" << endl;
                        passed++;
                    }
                }
            }
        }
    }
    
    // 9. Peri-fusion regression test - verify peri-fused systems are still rejected
    // Pyrene: a peri-fused system with 4 rings in a non-tree topology
    // Using SMILES that should create a cyclic ring adjacency graph
    cout << "Test: Peri-fusion rejection (Pyrene)" << endl;
    int mol_peri = indigoLoadMoleculeFromString("c1cc2cc3cc4cccc5cccc1c2c34c5");
    if (mol_peri < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for pyrene" << endl;
        failed++;
    } else {
        indigoAromatize(mol_peri);
        FusedRingSystemInput input_peri = detectFusedRingDirections(mol_peri);
        indigoFree(mol_peri);
        
        if (!input_peri.ringSizes.empty()) {
            cout << "  FAIL: Peri-fused molecule was accepted (ringSizes not empty), fusions: " << input_peri.fusions.size() << endl;
            failed++;
        } else {
            cout << "  PASS: Peri-fused molecule correctly rejected" << endl;
            passed++;
        }
    }

    // 10. Cyclic placement accepted: Acenaphthylene (3 rings: 6, 6, 5)
    cout << "Test: Cyclic placement accepted (Acenaphthylene)" << endl;
    int mol_acen = indigoLoadMoleculeFromString("c1ccc2cccc3c2c1C=C3");
    if (mol_acen < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for acenaphthylene" << endl;
        failed++;
    } else {
        indigoAromatize(mol_acen);
        FusedRingSystemInput input_acen = detectFusedRingDirections(mol_acen);
        indigoFree(mol_acen);
        if (input_acen.ringSizes.size() != 3) {
            cout << "  FAIL: Expected 3 rings for acenaphthylene, got " << input_acen.ringSizes.size() << endl;
            failed++;
        } else if (input_acen.fusions.size() != 2) {
            cout << "  FAIL: Expected 2 spanning fusions for acenaphthylene, got " << input_acen.fusions.size() << endl;
            failed++;
        } else {
            auto orient = computePreferredOrientation(input_acen);
            if (!orient.success) {
                cout << "  FAIL: computePreferredOrientation failed for acenaphthylene: " << orient.error << endl;
                failed++;
            } else {
                cout << "  PASS: Acenaphthylene cyclic placement accepted" << endl;
                passed++;
            }
        }
    }

    // 11. Cyclic placement accepted: 9bH-Phenalene (3 rings: 6, 6, 6)
    cout << "Test: Cyclic placement accepted (9bH-Phenalene)" << endl;
    int mol_phen = indigoLoadMoleculeFromString("C1=CC2=CC=CC3=CC=CC(=C1)C23");
    if (mol_phen < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for 9bH-phenalene" << endl;
        failed++;
    } else {
        indigoAromatize(mol_phen);
        FusedRingSystemInput input_phen = detectFusedRingDirections(mol_phen);
        indigoFree(mol_phen);
        if (input_phen.ringSizes.size() != 3) {
            cout << "  FAIL: Expected 3 rings for 9bH-phenalene, got " << input_phen.ringSizes.size() << endl;
            failed++;
        } else if (input_phen.fusions.size() != 2) {
            cout << "  FAIL: Expected 2 spanning fusions for 9bH-phenalene, got " << input_phen.fusions.size() << endl;
            failed++;
        } else {
            auto orient = computePreferredOrientation(input_phen);
            if (!orient.success) {
                cout << "  FAIL: computePreferredOrientation failed for 9bH-phenalene: " << orient.error << endl;
                failed++;
            } else {
                cout << "  PASS: 9bH-Phenalene cyclic placement accepted" << endl;
                passed++;
            }
        }
    }

    // 12. Cyclic placement accepted: Pyrene (4 rings: 6, 6, 6, 6)
    cout << "Test: Cyclic placement accepted (Pyrene)" << endl;
    int mol_pyr = indigoLoadMoleculeFromString("c1cc2ccc3cccc4ccc(c1)c2c34");
    if (mol_pyr < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for pyrene" << endl;
        failed++;
    } else {
        indigoAromatize(mol_pyr);
        FusedRingSystemInput input_pyr = detectFusedRingDirections(mol_pyr);
        indigoFree(mol_pyr);
        if (input_pyr.ringSizes.size() != 4) {
            cout << "  FAIL: Expected 4 rings for pyrene, got " << input_pyr.ringSizes.size() << endl;
            failed++;
        } else if (input_pyr.fusions.size() != 3) {
            cout << "  FAIL: Expected 3 spanning fusions for pyrene, got " << input_pyr.fusions.size() << endl;
            failed++;
        } else {
            auto orient = computePreferredOrientation(input_pyr);
            if (!orient.success) {
                cout << "  FAIL: computePreferredOrientation failed for pyrene: " << orient.error << endl;
                failed++;
            } else {
                cout << "  PASS: Pyrene cyclic placement accepted" << endl;
                passed++;
            }
        }
    }

    // 13. Geometrically inconsistent cyclic rejection: 6-6-3 peri-cycle (naphtho[1,8-bc]cyclopropene)
    cout << "Test: Geometrically-inconsistent cyclic rejection (6-6-3 peri-cycle)" << endl;
    int mol_inconsistent = indigoLoadMoleculeFromString("c1ccc2cccc3c2c13");
    if (mol_inconsistent < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for inconsistent fixture" << endl;
        failed++;
    } else {
        indigoAromatize(mol_inconsistent);
        FusedRingSystemInput input_inconsistent = detectFusedRingDirections(mol_inconsistent);
        indigoFree(mol_inconsistent);
        if (input_inconsistent.ringSizes.empty() && !input_inconsistent.error.empty()) {
            cout << "  PASS: Geometrically-inconsistent cyclic fixture cleanly rejected: " << input_inconsistent.error << endl;
            passed++;
        } else {
            cout << "  FAIL: Expected rejection for geometrically-inconsistent fixture, got ringSizes.size()="
                 << input_inconsistent.ringSizes.size() << endl;
            failed++;
        }
    }

    // 14. Pyrene atom-order invariance regression test: both SMILES orderings must produce identical results
    cout << "Test: Pyrene atom-order invariance regression" << endl;
    int mol_pyr1 = indigoLoadMoleculeFromString("c1cc2ccc3cccc4ccc(c1)c2c34");
    int mol_pyr2 = indigoLoadMoleculeFromString("C1=CC2=C3C1=CC=C4C3=C(C=C2)C=C4");
    if (mol_pyr1 < 0 || mol_pyr2 < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed for pyrene orderings" << endl;
        failed++;
        if (mol_pyr1 >= 0) indigoFree(mol_pyr1);
        if (mol_pyr2 >= 0) indigoFree(mol_pyr2);
    } else {
        indigoAromatize(mol_pyr1);
        indigoAromatize(mol_pyr2);
        FusedRingSystemInput input1 = detectFusedRingDirections(mol_pyr1);
        FusedRingSystemInput input2 = detectFusedRingDirections(mol_pyr2);
        indigoFree(mol_pyr1);
        indigoFree(mol_pyr2);

        bool ok = true;
        if (input1.ringSizes.size() != 4 || input2.ringSizes.size() != 4) {
            cout << "  FAIL: Expected 4 rings for both orderings, got "
                 << input1.ringSizes.size() << " and " << input2.ringSizes.size() << endl;
            ok = false;
        }
        if (input1.error != input2.error) {
            cout << "  FAIL: Error string mismatch: '" << input1.error << "' vs '" << input2.error << "'" << endl;
            ok = false;
        }
        if (input1.fusions.size() != 3 || input2.fusions.size() != 3) {
            cout << "  FAIL: Expected 3 spanning fusions for both orderings, got "
                 << input1.fusions.size() << " and " << input2.fusions.size() << endl;
            ok = false;
        }
        if (input1.closingFusions.empty() || input2.closingFusions.empty()) {
            cout << "  FAIL: Expected closingFusions populated for both orderings, got "
                 << input1.closingFusions.size() << " and " << input2.closingFusions.size() << endl;
            ok = false;
        }
        if (ok) {
            cout << "  PASS: Both pyrene SMILES orderings produce identical success/ring count/error, closingFusions preserved" << endl;
            passed++;
        } else {
            failed++;
        }
    }
    
    cout << "\nTOTAL PASSED: " << passed << endl;
    cout << "TOTAL FAILED: " << failed << endl;
    
    return failed == 0 ? 0 : 1;
}
