#include "FusedSystemTopology.h"
#include "indigo.h"
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

static int passed = 0;
static int failed = 0;

struct IndigoMolContext {
    Graph g;
    std::vector<std::set<int>> sssrRings;
    std::map<int, int> indigoToGraphIdx;
};

static IndigoMolContext loadMoleculeContext(int mol) {
    IndigoMolContext ctx;
    int atomIter = indigoIterateAtoms(mol);
    int idx = 0;
    while (int atom = indigoNext(atomIter)) {
        int origIdx = indigoIndex(atom);
        if (indigoAtomicNumber(atom) > 1) { // heavy atom
            ctx.indigoToGraphIdx[origIdx] = idx;
            GraphNode node;
            node.id = idx;
            node.indigoIdx = origIdx;
            node.atomicNumber = indigoAtomicNumber(atom);
            ctx.g.nodes.push_back(node);
            idx++;
        }
        indigoFree(atom);
    }
    indigoFree(atomIter);

    int bondIter = indigoIterateBonds(mol);
    while (int bond = indigoNext(bondIter)) {
        int src = indigoIndex(indigoSource(bond));
        int dst = indigoIndex(indigoDestination(bond));
        if (ctx.indigoToGraphIdx.count(src) && ctx.indigoToGraphIdx.count(dst)) {
            int u = ctx.indigoToGraphIdx[src];
            int v = ctx.indigoToGraphIdx[dst];
            int order = indigoBondOrder(bond);
            ctx.g.bonds.push_back({u, v, order});
            ctx.g.nodes[u].neighbors.push_back(v);
            ctx.g.nodes[u].bondOrders.push_back(order);
            ctx.g.nodes[v].neighbors.push_back(u);
            ctx.g.nodes[v].bondOrders.push_back(order);
        }
        indigoFree(bond);
    }
    indigoFree(bondIter);

    int sssrIter = indigoIterateSSSR(mol);
    if (sssrIter >= 0) {
        int subMol = 0;
        while ((subMol = indigoNext(sssrIter)) != 0) {
            std::set<int> rNodes;
            int ringAtomIter = indigoIterateAtoms(subMol);
            if (ringAtomIter >= 0) {
                int atomHandle = 0;
                while ((atomHandle = indigoNext(ringAtomIter)) != 0) {
                    int origIdx = indigoIndex(atomHandle);
                    if (ctx.indigoToGraphIdx.count(origIdx)) {
                        rNodes.insert(ctx.indigoToGraphIdx[origIdx]);
                    }
                    indigoFree(atomHandle);
                }
                indigoFree(ringAtomIter);
            }
            ctx.sssrRings.push_back(rNodes);
            indigoFree(subMol);
        }
        indigoFree(sssrIter);
    }

    return ctx;
}

static void test_naphthalene() {
    cout << "Test 1: Naphthalene (c1ccc2ccccc2c1)" << endl;
    int mol = indigoLoadMoleculeFromString("c1ccc2ccccc2c1");
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    indigoAromatize(mol);
    IndigoMolContext ctx = loadMoleculeContext(mol);
    indigoFree(mol);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);

    bool ok = true;
    if (topo.ringCount() != 2) {
        cout << "  FAIL: expected 2 rings, got " << topo.ringCount() << endl;
        ok = false;
    }

    // 0 peri centers, 0 spiro, 1 ortho
    if (!topo.periCenters().empty()) {
        cout << "  FAIL: expected 0 peri centers, got " << topo.periCenters().size() << endl;
        ok = false;
    }
    if (!topo.interiorAtoms().empty()) {
        cout << "  FAIL: expected 0 interior atoms, got " << topo.interiorAtoms().size() << endl;
        ok = false;
    }
    if (topo.peripheralAtoms().size() != 10) {
        cout << "  FAIL: expected 10 peripheral atoms, got " << topo.peripheralAtoms().size() << endl;
        ok = false;
    }
    if (topo.spiroJunctions().size() != 0) {
        cout << "  FAIL: expected 0 spiro junctions, got " << topo.spiroJunctions().size() << endl;
        ok = false;
    }
    if (topo.orthoJunctions().size() != 1) {
        cout << "  FAIL: expected 1 ortho junction, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    } else {
        const auto &oj = topo.orthoJunctions()[0];
        if (oj.ring1 != 0 || oj.ring2 != 1) {
            cout << "  FAIL: ortho junction rings mismatch: " << oj.ring1 << ", " << oj.ring2 << endl;
            ok = false;
        }
        // Hand-derived shared atoms in SMILES 'c1ccc2ccccc2c1': atoms 3 and 8
        if (oj.atom1 != 3 || oj.atom2 != 8) {
            cout << "  FAIL: ortho junction atoms mismatch: (" << oj.atom1 << ", " << oj.atom2 
                 << ") vs expected (3, 8)" << endl;
            ok = false;
        }
    }

    if (topo.junctionType(0, 1) != FusedBondType::Ortho) {
        cout << "  FAIL: junctionType(0, 1) != Ortho" << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: Naphthalene (2 rings, 1 ortho [3,8], 0 peri, 0 spiro)" << endl;
        passed++;
    } else {
        failed++;
    }
}

static void test_acenaphthylene() {
    cout << "Test 2: Acenaphthylene (c1ccc2cccc3c2c1C=C3)" << endl;
    int mol = indigoLoadMoleculeFromString("c1ccc2cccc3c2c1C=C3");
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    indigoAromatize(mol);
    IndigoMolContext ctx = loadMoleculeContext(mol);
    indigoFree(mol);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);

    bool ok = true;
    if (topo.ringCount() != 3) {
        cout << "  FAIL: expected 3 rings, got " << topo.ringCount() << endl;
        ok = false;
    }

    // Hand-derived peri center: atom 8 shared by all 3 rings
    std::set<int> expPeri = {8};
    if (topo.periCenters() != expPeri) {
        cout << "  FAIL: peri centers mismatch. Expected {8}, got: ";
        for (int a : topo.periCenters()) cout << a << " ";
        cout << endl;
        ok = false;
    }
    if (topo.interiorAtoms() != expPeri) {
        cout << "  FAIL: interior atoms mismatch" << endl;
        ok = false;
    }
    if (topo.peripheralAtoms().size() != 11 || topo.peripheralAtoms().count(8)) {
        cout << "  FAIL: peripheral atoms mismatch (expected 11 atoms excluding 8)" << endl;
        ok = false;
    }

    // Check peri junction details
    if (topo.periJunctions().size() != 1) {
        cout << "  FAIL: expected 1 peri junction, got " << topo.periJunctions().size() << endl;
        ok = false;
    } else {
        const auto &pj = topo.periJunctions()[0];
        if (pj.atom != 8) {
            cout << "  FAIL: peri junction atom expected 8, got " << pj.atom << endl;
            ok = false;
        }
        std::vector<int> expRings = {0, 1, 2};
        if (pj.rings != expRings) {
            cout << "  FAIL: peri junction rings mismatch" << endl;
            ok = false;
        }
    }

    // Zero spiro junctions
    if (!topo.spiroJunctions().empty()) {
        cout << "  FAIL: expected 0 spiro junctions, got " << topo.spiroJunctions().size() << endl;
        ok = false;
    }

    // Rest ortho junctions: exactly 3 ortho junctions
    if (topo.orthoJunctions().size() != 3) {
        cout << "  FAIL: expected 3 ortho junctions, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    }

    // Hand-derived pairs: (0,1) shares {3,8}, (0,2) shares {8,9}, (1,2) shares {7,8}
    // Check junctionType between all pairs
    if (topo.junctionType(0, 1) != FusedBondType::Ortho ||
        topo.junctionType(0, 2) != FusedBondType::Ortho ||
        topo.junctionType(1, 2) != FusedBondType::Ortho)
    {
        cout << "  FAIL: expected Ortho for all ring pairs" << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: Acenaphthylene (3 rings, 1 peri center [atom 8], 3 ortho junctions, 0 spiro)" << endl;
        passed++;
    } else {
        failed++;
    }
}

static void test_phenalene() {
    cout << "Test 3: 9bH-phenalene (C1=CC2=CC=CC3=CC=CC(=C1)C23)" << endl;
    int mol = indigoLoadMoleculeFromString("C1=CC2=CC=CC3=CC=CC(=C1)C23");
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    indigoAromatize(mol);
    IndigoMolContext ctx = loadMoleculeContext(mol);
    indigoFree(mol);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);

    bool ok = true;
    if (topo.ringCount() != 3) {
        cout << "  FAIL: expected 3 rings, got " << topo.ringCount() << endl;
        ok = false;
    }

    // Hand-derived peri center: atom 12 (central carbon shared by all 3 rings)
    // Matches Stage A1 detectPeriFusion fixture
    std::set<int> expPeri = {12};
    if (topo.periCenters() != expPeri) {
        cout << "  FAIL: peri centers mismatch. Expected {12}, got: ";
        for (int a : topo.periCenters()) cout << a << " ";
        cout << endl;
        ok = false;
    }
    if (topo.interiorAtoms() != expPeri) {
        cout << "  FAIL: interior atoms mismatch" << endl;
        ok = false;
    }
    if (topo.peripheralAtoms().size() != 12 || topo.peripheralAtoms().count(12)) {
        cout << "  FAIL: peripheral atoms mismatch (expected 12 atoms excluding 12)" << endl;
        ok = false;
    }

    if (topo.periJunctions().size() != 1) {
        cout << "  FAIL: expected 1 peri junction, got " << topo.periJunctions().size() << endl;
        ok = false;
    } else {
        const auto &pj = topo.periJunctions()[0];
        if (pj.atom != 12) {
            cout << "  FAIL: peri junction atom expected 12, got " << pj.atom << endl;
            ok = false;
        }
        std::vector<int> expRings = {0, 1, 2};
        if (pj.rings != expRings) {
            cout << "  FAIL: peri junction rings mismatch" << endl;
            ok = false;
        }
    }

    if (!topo.spiroJunctions().empty()) {
        cout << "  FAIL: expected 0 spiro junctions, got " << topo.spiroJunctions().size() << endl;
        ok = false;
    }

    if (topo.orthoJunctions().size() != 3) {
        cout << "  FAIL: expected 3 ortho junctions, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: 9bH-phenalene (3 rings, 1 central peri center [atom 12], 3 ortho junctions, 0 spiro)" << endl;
        passed++;
    } else {
        failed++;
    }
}

static void test_pyrene() {
    cout << "Test 4: Pyrene (c1cc2ccc3cccc4ccc(c1)c2c34)" << endl;
    int mol = indigoLoadMoleculeFromString("c1cc2ccc3cccc4ccc(c1)c2c34");
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    indigoAromatize(mol);
    IndigoMolContext ctx = loadMoleculeContext(mol);
    indigoFree(mol);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);

    bool ok = true;
    if (topo.ringCount() != 4) {
        cout << "  FAIL: expected 4 rings, got " << topo.ringCount() << endl;
        ok = false;
    }

    // Hand-derived peri centers: exactly 2 peri centers, atoms {14, 15}
    std::set<int> expPeri = {14, 15};
    if (topo.periCenters() != expPeri) {
        cout << "  FAIL: peri centers mismatch. Expected {14, 15}, got: ";
        for (int a : topo.periCenters()) cout << a << " ";
        cout << endl;
        ok = false;
    }
    if (topo.interiorAtoms() != expPeri) {
        cout << "  FAIL: interior atoms mismatch" << endl;
        ok = false;
    }
    if (topo.peripheralAtoms().size() != 14 || topo.peripheralAtoms().count(14) || topo.peripheralAtoms().count(15)) {
        cout << "  FAIL: peripheral atoms mismatch (expected 14 atoms)" << endl;
        ok = false;
    }

    if (topo.periJunctions().size() != 2) {
        cout << "  FAIL: expected 2 peri junctions, got " << topo.periJunctions().size() << endl;
        ok = false;
    } else {
        // Both peri centers must meet exactly 3 rings each
        for (const auto &pj : topo.periJunctions()) {
            if (pj.rings.size() != 3) {
                cout << "  FAIL: peri junction atom " << pj.atom << " does not meet 3 rings (got " 
                     << pj.rings.size() << ")" << endl;
                ok = false;
            }
        }
    }

    if (!topo.spiroJunctions().empty()) {
        cout << "  FAIL: expected 0 spiro junctions, got " << topo.spiroJunctions().size() << endl;
        ok = false;
    }

    // 5 ortho fusion bonds in pyrene
    if (topo.orthoJunctions().size() != 5) {
        cout << "  FAIL: expected 5 ortho junctions, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: Pyrene (4 rings, 2 peri centers [atoms 14, 15], 5 ortho junctions, 0 spiro)" << endl;
        passed++;
    } else {
        failed++;
    }
}

static void test_spiro_decane() {
    cout << "Test 5: Spiro[4.5]decane (C1CCCC12CCCCC2)" << endl;
    int mol = indigoLoadMoleculeFromString("C1CCCC12CCCCC2");
    if (mol < 0) {
        cout << "  FAIL: indigoLoadMoleculeFromString failed" << endl;
        failed++;
        return;
    }
    IndigoMolContext ctx = loadMoleculeContext(mol);
    indigoFree(mol);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);

    bool ok = true;
    if (topo.ringCount() != 2) {
        cout << "  FAIL: expected 2 rings, got " << topo.ringCount() << endl;
        ok = false;
    }

    // 0 peri centers, 0 ortho junctions
    if (!topo.periCenters().empty()) {
        cout << "  FAIL: expected 0 peri centers, got " << topo.periCenters().size() << endl;
        ok = false;
    }
    if (!topo.interiorAtoms().empty()) {
        cout << "  FAIL: expected 0 interior atoms, got " << topo.interiorAtoms().size() << endl;
        ok = false;
    }
    if (topo.peripheralAtoms().size() != 10) {
        cout << "  FAIL: expected 10 peripheral atoms, got " << topo.peripheralAtoms().size() << endl;
        ok = false;
    }
    if (!topo.orthoJunctions().empty()) {
        cout << "  FAIL: expected 0 ortho junctions, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    }

    // Exactly 1 spiro junction at atom 4
    if (topo.spiroJunctions().size() != 1) {
        cout << "  FAIL: expected 1 spiro junction, got " << topo.spiroJunctions().size() << endl;
        ok = false;
    } else {
        const auto &sj = topo.spiroJunctions()[0];
        if (sj.ring1 != 0 || sj.ring2 != 1) {
            cout << "  FAIL: spiro rings mismatch: " << sj.ring1 << ", " << sj.ring2 << endl;
            ok = false;
        }
        if (sj.atom != 4) {
            cout << "  FAIL: spiro atom expected 4, got " << sj.atom << endl;
            ok = false;
        }
    }

    if (topo.junctionType(0, 1) != FusedBondType::Spiro) {
        cout << "  FAIL: junctionType(0, 1) != Spiro" << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: Spiro[4.5]decane (2 rings, 1 spiro atom [atom 4], 0 ortho, 0 peri)" << endl;
        passed++;
    } else {
        failed++;
    }
}

static void test_hand_crafted() {
    cout << "Test 6: Hand-crafted topology edge cases" << endl;
    // 3 rings: Ring 0={0,1,2,3}, Ring 1={2,3,4,5}, Ring 2={4,5,6,7}
    // Bonds: (0,1), (1,2), (2,3), (3,0), (3,4), (4,5), (5,2), (5,6), (6,7), (7,4)
    // Here ring 0 and ring 1 share {2, 3} (bonded) -> Ortho
    // Ring 1 and ring 2 share {4, 5} (bonded) -> Ortho
    // Ring 0 and ring 2 are disjoint -> None
    std::vector<std::set<int>> rings = {
        {0, 1, 2, 3},
        {2, 3, 4, 5},
        {4, 5, 6, 7}
    };
    std::vector<GraphBond> bonds = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},
        {3, 4, 1}, {4, 5, 1}, {5, 2, 1},
        {5, 6, 1}, {6, 7, 1}, {7, 4, 1}
    };

    FusedSystemTopology topo(rings, bonds);
    bool ok = true;
    if (topo.orthoJunctions().size() != 2) {
        cout << "  FAIL: expected 2 ortho junctions, got " << topo.orthoJunctions().size() << endl;
        ok = false;
    }
    if (topo.periCenters().size() != 0) {
        cout << "  FAIL: expected 0 peri centers, got " << topo.periCenters().size() << endl;
        ok = false;
    }
    if (topo.adjacentRings(0).size() != 1 || topo.adjacentRings(0)[0].neighborRing != 1) {
        cout << "  FAIL: ring 0 adjacency mismatch" << endl;
        ok = false;
    }
    if (topo.adjacentRings(1).size() != 2) {
        cout << "  FAIL: ring 1 adjacency mismatch" << endl;
        ok = false;
    }

    if (ok) {
        cout << "  PASS: Hand-crafted linear ring chain" << endl;
        passed++;
    } else {
        failed++;
    }
}

int main() {
    setvbuf(stdout, nullptr, _IONBF, 0);
    qword sid = indigoAllocSessionId();
    indigoSetSessionId(sid);

    test_naphthalene();
    test_acenaphthylene();
    test_phenalene();
    test_pyrene();
    test_spiro_decane();
    test_hand_crafted();

    indigoReleaseSessionId(sid);

    cout << "\nTopology test results: " << passed << " passed, " << failed << " failed." << endl;
    return failed == 0 ? 0 : 1;
}
