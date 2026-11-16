#include "PeripheralNumbering.h"
#include "SystemLocant.h"
#include "FusedSystemTopology.h"
#include "indigo.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <algorithm>
#include <QString>

using namespace std;

static int passed = 0;
static int failed = 0;

#define EXPECT_EQ(actual, expected) do { \
    auto actual_val_ = (actual); \
    auto expected_val_ = (expected); \
    if (actual_val_ != expected_val_) { \
        std::cerr << "FAIL at " << __LINE__ << ": expected " << #expected << ", got " << #actual << "\n"; \
        failed++; \
    } else { \
        passed++; \
    } \
} while (0)

#define ASSERT_EQ(actual, expected) do { \
    auto actual_val_ = (actual); \
    auto expected_val_ = (expected); \
    if (actual_val_ != expected_val_) { \
        std::cerr << "ASSERT_EQ FAIL at " << __LINE__ << ": expected " << #expected << ", got " << #actual << "\n"; \
        failed++; \
        return; \
    } else { \
        passed++; \
    } \
} while (0)

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

static void test_pyrene_interior_iupac_and_cas() {
    cout << "Test: Pyrene interior atoms IUPAC (3a1, 5a1) and CAS (10b, 10c)" << endl;
    int m = indigoLoadMoleculeFromString("c1cc2ccc3cccc4ccc(c1)c2c34");
    if (m < 0) {
        cerr << "FAIL: Failed to load pyrene SMILES\n";
        failed++;
        return;
    }
    indigoAromatize(m);
    auto ctx = loadMoleculeContext(m);
    indigoFree(m);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);
    auto perim = topo.perimeterCycle();
    auto interior = topo.interiorAtoms();

    ASSERT_EQ(perim.size(), 14u);
    ASSERT_EQ(interior.size(), 2u);

    // Standard pyrene numbering: [3 non-BH, 1 BH, 2 non-BH, 1 BH, 3 non-BH, 1 BH, 2 non-BH, 1 BH]
    // Orient perim to start with 3 non-BH atoms (1, 2, 3), followed by 3a, 4, 5, 5a, 6, 7, 8, 8a, 9, 10, 10a
    auto isBh = [&](int atom) {
        int deg = 0;
        for (int nei : ctx.g.nodes[atom].neighbors) {
            if (topo.peripheralAtoms().count(nei) || interior.count(nei)) deg++;
        }
        return deg > 2;
    };

    std::vector<int> orientedPerim;
    int n = static_cast<int>(perim.size());
    for (int dir : {1, -1}) {
        for (int i = 0; i < n; ++i) {
            std::vector<int> cand;
            for (int k = 0; k < n; ++k) {
                int idx = (i + dir * k + n * 10) % n;
                cand.push_back(perim[idx]);
            }
            if (!isBh(cand[0]) && !isBh(cand[1]) && !isBh(cand[2]) && isBh(cand[3]) &&
                !isBh(cand[4]) && !isBh(cand[5]) && isBh(cand[6]) &&
                isBh(cand[n - 1])) {
                orientedPerim = cand;
                break;
            }
        }
        if (!orientedPerim.empty()) break;
    }

    ASSERT_EQ(orientedPerim.size(), 14u);

    // Call numberInteriorAtoms
    auto locants = numberInteriorAtoms(orientedPerim, ctx.g, interior);
    ASSERT_EQ(locants.size(), 2u);

    // IUPAC: 3a1, 5a1
    std::vector<QString> renderedIupac = {renderIupac(locants[0]), renderIupac(locants[1])};
    std::sort(renderedIupac.begin(), renderedIupac.end());
    EXPECT_EQ(renderedIupac, (std::vector<QString>{QString::fromUtf8("3a\u00B9"), QString::fromUtf8("5a\u00B9")}));

    // CAS: 10b, 10c
    std::vector<QString> renderedCas = {renderCas(locants[0]), renderCas(locants[1])};
    std::sort(renderedCas.begin(), renderedCas.end());
    EXPECT_EQ(renderedCas, (std::vector<QString>{"10b", "10c"}));
}

static void test_phenalene_interior_both_styles_agree() {
    cout << "Test: Phenalene single interior atom agrees in both styles (9b)" << endl;
    int m = indigoLoadMoleculeFromString("C1=CC2=CC=CC3=CC=CC(=C1)C23");
    if (m < 0) {
        cerr << "FAIL: Failed to load phenalene SMILES\n";
        failed++;
        return;
    }
    indigoAromatize(m);
    auto ctx = loadMoleculeContext(m);
    indigoFree(m);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);
    auto perim = topo.perimeterCycle();
    auto interior = topo.interiorAtoms();

    ASSERT_EQ(perim.size(), 12u);
    ASSERT_EQ(interior.size(), 1u);

    // Orient perim to start with 3 non-BH atoms (1, 2, 3), followed by 3a, 4, 5, 6, 6a, 7, 8, 9, 9a
    auto isBh = [&](int atom) {
        int deg = 0;
        for (int nei : ctx.g.nodes[atom].neighbors) {
            if (topo.peripheralAtoms().count(nei) || interior.count(nei)) deg++;
        }
        return deg > 2;
    };

    std::vector<int> orientedPerim;
    int n = static_cast<int>(perim.size());
    for (int dir : {1, -1}) {
        for (int i = 0; i < n; ++i) {
            std::vector<int> cand;
            for (int k = 0; k < n; ++k) {
                int idx = (i + dir * k + n * 10) % n;
                cand.push_back(perim[idx]);
            }
            if (!isBh(cand[0]) && !isBh(cand[1]) && !isBh(cand[2]) && isBh(cand[3]) &&
                isBh(cand[n - 1])) {
                orientedPerim = cand;
                break;
            }
        }
        if (!orientedPerim.empty()) break;
    }

    ASSERT_EQ(orientedPerim.size(), 12u);

    auto locants = numberInteriorAtoms(orientedPerim, ctx.g, interior);
    ASSERT_EQ(locants.size(), 1u);

    EXPECT_EQ(renderIupac(locants[0]), QString("9b"));
    EXPECT_EQ(renderCas(locants[0]), QString("9b"));
}

static void test_coronene_interior_locants() {
    cout << "Test: Coronene interior atoms IUPAC and CAS" << endl;
    int m = indigoLoadMoleculeFromString("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67");
    if (m < 0) {
        cerr << "FAIL: Failed to load coronene SMILES\n";
        failed++;
        return;
    }
    indigoAromatize(m);
    auto ctx = loadMoleculeContext(m);
    indigoFree(m);

    FusedSystemTopology topo(ctx.sssrRings, ctx.g);
    auto perim = topo.perimeterCycle();
    auto interior = topo.interiorAtoms();

    ASSERT_EQ(perim.size(), 18u);
    ASSERT_EQ(interior.size(), 6u);

    auto isBh = [&](int atom) {
        int deg = 0;
        for (int nei : ctx.g.nodes[atom].neighbors) {
            if (topo.peripheralAtoms().count(nei) || interior.count(nei)) deg++;
        }
        return deg > 2;
    };

    std::vector<int> orientedPerim;
    int n = static_cast<int>(perim.size());
    for (int dir : {1, -1}) {
        for (int i = 0; i < n; ++i) {
            std::vector<int> cand;
            for (int k = 0; k < n; ++k) {
                int idx = (i + dir * k + n * 10) % n;
                cand.push_back(perim[idx]);
            }
            if (!isBh(cand[0]) && !isBh(cand[1]) && isBh(cand[2]) && isBh(cand[n - 1])) {
                orientedPerim = cand;
                break;
            }
        }
        if (!orientedPerim.empty()) break;
    }

    ASSERT_EQ(orientedPerim.size(), 18u);

    auto locants = numberInteriorAtoms(orientedPerim, ctx.g, interior);
    ASSERT_EQ(locants.size(), 6u);

    std::vector<QString> renderedIupac;
    std::vector<QString> renderedCas;
    for (const auto &loc : locants) {
        renderedIupac.push_back(renderIupac(loc));
        renderedCas.push_back(renderCas(loc));
    }
    std::sort(renderedIupac.begin(), renderedIupac.end());
    std::sort(renderedCas.begin(), renderedCas.end());

    std::vector<QString> expectedIupac = {
        QString::fromUtf8("10a\u00B9"),
        QString::fromUtf8("12a\u00B9"),
        QString::fromUtf8("2a\u00B9"),
        QString::fromUtf8("4a\u00B9"),
        QString::fromUtf8("6a\u00B9"),
        QString::fromUtf8("8a\u00B9")
    };
    std::sort(expectedIupac.begin(), expectedIupac.end());
    EXPECT_EQ(renderedIupac, expectedIupac);

    std::vector<QString> expectedCas = {"12b", "12c", "12d", "12e", "12f", "12g"};
    EXPECT_EQ(renderedCas, expectedCas);
}

int main() {
    setvbuf(stdout, nullptr, _IONBF, 0);
    qword sid = indigoAllocSessionId();
    indigoSetSessionId(sid);

    test_pyrene_interior_iupac_and_cas();
    test_phenalene_interior_both_styles_agree();
    test_coronene_interior_locants();

    indigoReleaseSessionId(sid);

    cout << endl;
    cout << "TOTAL PASSED: " << passed << endl;
    cout << "TOTAL FAILED: " << failed << endl;

    return failed == 0 ? 0 : 1;
}
