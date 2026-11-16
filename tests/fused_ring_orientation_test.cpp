#include "../src/app/FusedRingOrientation.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

int passed = 0;
int failed = 0;

void run_test(const string& name, const FusedRingSystemInput& input, int exp_row, double exp_ur) {
    cout << "Test: " << name << endl;
    auto res = computePreferredOrientation(input);
    if (!res.success) {
        cout << "  FAIL: returned error: " << res.error << endl;
        failed++;
        return;
    }
    
    bool ok = true;
    if (res.ringsInHorizontalRow != exp_row) {
        cout << "  FAIL: expected ringsInHorizontalRow=" << exp_row << ", got " << res.ringsInHorizontalRow << endl;
        ok = false;
    }
    if (abs(res.ringsInUpperRightQuadrant - exp_ur) > 1e-5) {
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
    // Case 1: Naphthalene
    FusedRingSystemInput naph;
    naph.ringSizes = {6, 6};
    naph.fusions = {{0, 1, 0}};
    run_test("Naphthalene", naph, 2, 0.5); 

    // Case 2: Anthracene
    FusedRingSystemInput anthra;
    anthra.ringSizes = {6, 6, 6};
    anthra.fusions = {{0, 1, 0}, {1, 2, 3}};
    run_test("Anthracene", anthra, 3, 0.75);

    // Case 3: Phenanthrene
    FusedRingSystemInput phenan;
    phenan.ringSizes = {6, 6, 6};
    phenan.fusions = {{0, 1, 0}, {1, 2, 2}};
    run_test("Phenanthrene", phenan, 2, 1.5);

    // Case 4: Naphthacene (Tetracene)
    FusedRingSystemInput tetra;
    tetra.ringSizes = {6, 6, 6, 6};
    tetra.fusions = {{0, 1, 0}, {1, 2, 3}, {2, 3, 3}};
    run_test("Naphthacene", tetra, 4, 1.0);
    
    // Case 5: Branching / Triple junction cycle (Phenalene non-tree system)
    FusedRingSystemInput phenal;
    phenal.ringSizes = {6, 6, 6};
    phenal.fusions = {{0, 1, 0}, {1, 2, 5}, {2, 0, 5}};
    cout << "Test: Phenalene" << endl;
    auto res = computePreferredOrientation(phenal);
    if (res.success && res.error.empty()) {
        cout << "  PASS: Phenalene non-tree system is oriented" << endl;
        passed++;
    } else {
        cout << "  FAIL: Expected phenalene orientation success, got error: " << res.error << endl;
        failed++;
    }
    
    // Case 6: Branching molecule (degree-3 central ring) - verify new ringPos field
    cout << "Test: Branching molecule (ringPos)" << endl;
    FusedRingSystemInput branch_input;
    branch_input.ringSizes = {6, 6, 6, 6};
    branch_input.fusions = {{0, 1, 0}, {0, 2, 2}, {0, 3, 4}};
    
    auto res_branch = computePreferredOrientation(branch_input);
    if (!res_branch.success) {
        cout << "  FAIL: computePreferredOrientation failed: " << res_branch.error << endl;
        failed++;
    } else if (res_branch.ringPos.size() != 4) {
        cout << "  FAIL: Expected ringPos.size() == 4, got " << res_branch.ringPos.size() << endl;
        failed++;
    } else {
        bool allDistinct = true;
        for (size_t i = 0; i < res_branch.ringPos.size(); ++i) {
            for (size_t j = i + 1; j < res_branch.ringPos.size(); ++j) {
                double dx = res_branch.ringPos[i].first - res_branch.ringPos[j].first;
                double dy = res_branch.ringPos[i].second - res_branch.ringPos[j].second;
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
            cout << "  PASS: Branching molecule accepted, ringPos has 4 distinct positions" << endl;
            passed++;
        }
    }

    // Case 7: Pyrene (4 rings in 2x2 grid with cyclic closure)
    {
        FusedRingSystemInput pyrene;
        pyrene.ringSizes = {6, 6, 6, 6};
        pyrene.fusions = {{0, 1, 0}, {0, 2, 1}, {2, 3, 2}, {1, 3, 4}, {1, 2, 5}};
        cout << "Test: Pyrene non-tree system" << endl;
        auto res_pyrene = computePreferredOrientation(pyrene);
        if (res_pyrene.success && res_pyrene.error.empty()) {
            cout << "  PASS: Pyrene non-tree system is oriented" << endl;
            passed++;
        } else {
            cout << "  FAIL: Expected pyrene orientation success, got error: " << res_pyrene.error << endl;
            failed++;
        }
    }

    // Case 8: Geometrically impossible cyclic closure (mismatched sizes)
    {
        FusedRingSystemInput bad_cycle;
        bad_cycle.ringSizes = {5, 5, 5};
        bad_cycle.fusions = {{0, 1, 0}, {1, 2, 2}, {2, 0, 2}};
        cout << "Test: Geometrically impossible cyclic closure" << endl;
        auto res_bad = computePreferredOrientation(bad_cycle);
        if (!res_bad.success && !res_bad.error.empty() && res_bad.error.find("Non-tree") == string::npos) {
            cout << "  PASS: Impossible closure rejected with geometric error: " << res_bad.error << endl;
            passed++;
        } else {
            cout << "  FAIL: Expected geometric contradiction error, got: " << res_bad.error << endl;
            failed++;
        }
    }
    
    cout << endl;
    cout << "TOTAL PASSED: " << passed << endl;
    cout << "TOTAL FAILED: " << failed << endl;
    
    return failed == 0 ? 0 : 1;
}
