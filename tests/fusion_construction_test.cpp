#include "FusionConstruction.h"
#include "FusedRingDirectionDetector.h"
#include "FusedSystemTopology.h"
#include "PeripheralNumbering.h"
#include "indigo.h"
#include "indigo-inchi.h"
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>

static int passed = 0;
static int failed = 0;

#define EXPECT_TRUE(cond) do { \
    if (!(cond)) { \
        std::cerr << "FAIL at line " << __LINE__ << ": expected true for " << #cond << "\n"; \
        failed++; \
    } else { \
        passed++; \
    } \
} while(0)

#define EXPECT_FALSE(cond) do { \
    if (cond) { \
        std::cerr << "FAIL at line " << __LINE__ << ": expected false for " << #cond << "\n"; \
        failed++; \
    } else { \
        passed++; \
    } \
} while(0)

#define EXPECT_EQ(actual, expected) do { \
    auto a_ = (actual); \
    auto e_ = (expected); \
    if (a_ != e_) { \
        std::cerr << "FAIL at line " << __LINE__ << ": expected " << #expected << ", got " << a_ << "\n"; \
        failed++; \
    } else { \
        passed++; \
    } \
} while(0)

struct TestContext {
    Graph g;
    std::vector<std::set<int>> sssrRings;
    std::map<int, int> indigoToGraphIdx;
};

static TestContext loadContext(int mol) {
    TestContext ctx;
    int atomIter = indigoIterateAtoms(mol);
    if (atomIter >= 0) {
        int a;
        while ((a = indigoNext(atomIter)) != 0) {
            int origIdx = indigoIndex(a);
            int gIdx = static_cast<int>(ctx.g.nodes.size());
            ctx.indigoToGraphIdx[origIdx] = gIdx;
            GraphNode node;
            node.id = gIdx;
            node.atomicNumber = indigoAtomicNumber(a);
            ctx.g.nodes.push_back(node);
            indigoFree(a);
        }
        indigoFree(atomIter);
    }
    int bondIter = indigoIterateBonds(mol);
    if (bondIter >= 0) {
        int b;
        while ((b = indigoNext(bondIter)) != 0) {
            int src = indigoIndex(indigoSource(b));
            int dst = indigoIndex(indigoDestination(b));
            int order = indigoBondOrder(b);
            if (ctx.indigoToGraphIdx.count(src) && ctx.indigoToGraphIdx.count(dst)) {
                int u = ctx.indigoToGraphIdx[src];
                int v = ctx.indigoToGraphIdx[dst];
                ctx.g.bonds.push_back({u, v, order});
                ctx.g.nodes[u].neighbors.push_back(v);
                ctx.g.nodes[v].neighbors.push_back(u);
            }
            indigoFree(b);
        }
        indigoFree(bondIter);
    }
    int sssrIter = indigoIterateSSSR(mol);
    if (sssrIter >= 0) {
        int subMol;
        while ((subMol = indigoNext(sssrIter)) != 0) {
            std::set<int> rNodes;
            int ringAtomIter = indigoIterateAtoms(subMol);
            if (ringAtomIter >= 0) {
                int a;
                while ((a = indigoNext(ringAtomIter)) != 0) {
                    int origIdx = indigoIndex(a);
                    if (ctx.indigoToGraphIdx.count(origIdx)) {
                        rNodes.insert(ctx.indigoToGraphIdx[origIdx]);
                    }
                    indigoFree(a);
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

static void test_hits_fr32_synthetic() {
    std::cout << "Test: hitsFr32Limitation synthetic tree tests\n";
    // 3 components in a triangle (0-1, 1-2, 2-0)
    ComponentTree triangleTree;
    triangleTree.components.resize(3);
    triangleTree.junctions = {
        {0, 1, {}, {}},
        {1, 2, {}, {}},
        {2, 0, {}, {}}
    };
    EXPECT_TRUE(hitsFr32Limitation(triangleTree));

    // 3 components in a line (0-1, 0-2)
    ComponentTree chainTree;
    chainTree.components.resize(3);
    chainTree.junctions = {
        {0, 1, {}, {}},
        {0, 2, {}, {}}
    };
    EXPECT_FALSE(hitsFr32Limitation(chainTree));
}

static void test_decomposes_clean_three_component_chain() {
    std::cout << "Test: DecomposesCleanThreeComponentChain (dibenzo[c,g]phenanthrene)\n";
    // dibenzo[c,g]phenanthrene ([5]helicene)
    int m = indigoLoadMoleculeFromString("C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=CC=CC=C54");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput input = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, input, m);
    indigoFree(m);

    for (const auto &c : tree.components) {
        std::cout << "  Comp: " << c.name.toStdString() << " (rings=" << c.ringCount << ")\n";
    }

    EXPECT_TRUE(tree.valid);
    EXPECT_EQ(tree.components.size(), 3u);
    EXPECT_FALSE(hitsFr32Limitation(tree));

    // Confirm parent is phenanthrene (3 rings)
    EXPECT_TRUE(tree.parentIndex >= 0);
    if (tree.parentIndex >= 0) {
        EXPECT_EQ(tree.components[tree.parentIndex].name.toStdString(), "phenanthrene");
        EXPECT_EQ(tree.components[tree.parentIndex].ringCount, 3);
    }
}

static void test_detects_fr32_limitation_on_mutually_fused_pair() {
    std::cout << "Test: DetectsFr32LimitationOnMutuallyFusedPair (acepentalene)\n";
    // acepentalene / cyclopenta[cd]pentalene: 3 mutually-fused 5-rings
    int m = indigoLoadMoleculeFromString("C1=CC2=CC=C3C2=C1C=C3");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput input = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, input, m);
    indigoFree(m);

    for (const auto &c : tree.components) {
        std::cout << "  Comp: " << c.name.toStdString() << " (rings=" << c.ringCount << ")\n";
    }

    EXPECT_TRUE(tree.valid);
    EXPECT_EQ(tree.components.size(), 3u);
    EXPECT_TRUE(hitsFr32Limitation(tree));
}

static void test_parent_selection_determinism() {
    std::cout << "Test: ParentSelectionDeterministicAcrossAtomOrder\n";
    const char *s1 = "C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=CC=CC=C54";
    const char *s2 = "c1ccc2c(c1)ccc3c2c4c(cc3)ccc5ccccc45";

    int m1 = indigoLoadMoleculeFromString(s1);
    indigoAromatize(m1);
    auto ctx1 = loadContext(m1);
    FusedRingSystemInput in1 = detectFusedRingDirections(m1);
    ComponentTree tree1 = decomposeComponentTree(ctx1.g, ctx1.sssrRings, in1, m1);
    indigoFree(m1);

    int m2 = indigoLoadMoleculeFromString(s2);
    indigoAromatize(m2);
    auto ctx2 = loadContext(m2);
    FusedRingSystemInput in2 = detectFusedRingDirections(m2);
    ComponentTree tree2 = decomposeComponentTree(ctx2.g, ctx2.sssrRings, in2, m2);
    indigoFree(m2);

    EXPECT_TRUE(tree1.valid);
    EXPECT_TRUE(tree2.valid);
    EXPECT_EQ(tree1.components.size(), tree2.components.size());
    EXPECT_EQ(tree1.components[tree1.parentIndex].name.toStdString(), tree2.components[tree2.parentIndex].name.toStdString());
    EXPECT_EQ(tree1.components[tree1.parentIndex].ringCount, tree2.components[tree2.parentIndex].ringCount);
}

static void test_citation_order_fr36() {
    std::cout << "Test: CitationOrderFR36 (furo[2',3':4,5]pyrrolo[2,3-b]imidazo[4,5-e]pyrazine)\n";
    ComponentTree tree;
    tree.valid = true;
    tree.parentIndex = 0;
    FusionComponent c0; c0.name = "pyrazine"; c0.prefix = "pyrazino"; c0.order = 0;
    FusionComponent c1; c1.name = "pyrrole"; c1.prefix = "pyrrolo"; c1.order = 1;
    FusionComponent c2; c2.name = "imidazole"; c2.prefix = "imidazo"; c2.order = 1;
    FusionComponent c3; c3.name = "furan"; c3.prefix = "furo"; c3.order = 2;
    tree.components = {c0, c1, c2, c3};
    tree.junctions = {
        {0, 1, {}, {}},
        {0, 2, {}, {}},
        {1, 3, {}, {}}
    };

    QStringList order = orderCitation(tree);
    EXPECT_EQ(order.size(), 4);
    if (order.size() == 4) {
        EXPECT_EQ(order[0].toStdString(), "furo");
        EXPECT_EQ(order[1].toStdString(), "pyrrolo");
        EXPECT_EQ(order[2].toStdString(), "imidazo");
        EXPECT_EQ(order[3].toStdString(), "pyrazine");
    }
}

static void test_benzoheterocycle_merging_fr35() {
    std::cout << "Test: BenzoheterocycleMergingFR35 (thieno[3,2-f][2,1]benzothiazole)\n";

    int m = indigoLoadMoleculeFromString("c1cc2c(s1)cc3nscc3c2");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    applyBenzoheterocycleMerging(tree, ctx.g, m);
    QString constructed = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(constructed.toStdString(), "thieno[3,2-f][2,1]benzothiazole");
    indigoFree(m);

    EXPECT_TRUE(tree.valid);
    EXPECT_EQ(tree.components.size(), 2u);
    QStringList citation = orderCitation(tree);
    EXPECT_EQ(citation.size(), 2);
    if (citation.size() == 2) {
        EXPECT_EQ(citation[0].toStdString(), "thieno");
        EXPECT_TRUE(citation[1].contains("benzothiazole"));
    }
}

static void test_benzoheterocycle_merging_skipped_on_multiplicative() {
    std::cout << "Test: BenzoheterocycleMergingSkippedOnMultiplicative (dibenzo[c,g]phenanthrene)\n";
    int m = indigoLoadMoleculeFromString("C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=CC=CC=C54");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    for (const auto &c : tree.components) {
        std::cout << "  Before merge: " << c.name.toStdString() << " (rings=" << c.ringCount << ", atoms=" << c.atomIndices.size() << ")\n";
    }

    applyBenzoheterocycleMerging(tree, ctx.g, m);
    QString constructed = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(constructed.toStdString(), "dibenzo[c,g]phenanthrene");
    indigoFree(m);

    EXPECT_TRUE(tree.valid);
    EXPECT_EQ(tree.components.size(), 3u);
    QStringList citation = orderCitation(tree);
    EXPECT_EQ(citation.size(), 3);
    if (citation.size() == 3) {
        EXPECT_EQ(citation[0].toStdString(), "benzo");
        EXPECT_EQ(citation[1].toStdString(), "benzo");
        EXPECT_EQ(citation[2].toStdString(), "phenanthrene");
    }
}

static void test_descriptor_encoding_dipyrido_pyrazine() {
    std::cout << "Test: DescriptorEncoding (dipyrido[1,2-a:2',1'-c]pyrazine)\n";
    int m = indigoCreateMolecule();
    int a[14];
    for (int i = 0; i < 14; ++i) {
        if (i == 0 || i == 3) a[i] = indigoAddAtom(m, "N");
        else a[i] = indigoAddAtom(m, "C");
    }
    // Pyrazine ring
    indigoAddBond(a[0], a[1], 4);
    indigoAddBond(a[1], a[2], 4);
    indigoAddBond(a[2], a[3], 4);
    indigoAddBond(a[3], a[4], 4);
    indigoAddBond(a[4], a[5], 4);
    indigoAddBond(a[5], a[0], 4);
    // Pyridine 1
    indigoAddBond(a[0], a[9], 4);
    indigoAddBond(a[9], a[8], 4);
    indigoAddBond(a[8], a[7], 4);
    indigoAddBond(a[7], a[6], 4);
    indigoAddBond(a[6], a[1], 4);
    // Pyridine 2
    indigoAddBond(a[2], a[10], 4);
    indigoAddBond(a[10], a[11], 4);
    indigoAddBond(a[11], a[12], 4);
    indigoAddBond(a[12], a[13], 4);
    indigoAddBond(a[13], a[3], 4);

    std::cout << "  Constructed smiles: " << indigoSmiles(m) << "\n";
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);
    for (const auto &c : tree.components) {
        std::cout << "  Candidate 4 comp: " << c.name.toStdString() << " (rings=" << c.ringCount << ")\n";
    }

    EXPECT_TRUE(tree.valid);
    EXPECT_EQ(tree.components.size(), 3u);
    if (tree.parentIndex >= 0 && tree.parentIndex < static_cast<int>(tree.components.size())) {
        EXPECT_EQ(tree.components[tree.parentIndex].name.toStdString(), "pyrazine");
    }
    QString desc = encodeFusionDescriptor(tree, ctx.g, m);
    EXPECT_EQ(desc.toStdString(), "[1,2-a:2',1'-c]");
    indigoFree(m);
}

static void test_descriptor_encoding_naphtho_difuran() {
    std::cout << "Test: DescriptorEncoding (naphtho[1,2-c:7,8-c']difuran)\n";
    int m = indigoCreateMolecule();
    int at[16];
    for (int i = 0; i < 16; ++i) {
        if (i == 11 || i == 14) at[i] = indigoAddAtom(m, "O");
        else at[i] = indigoAddAtom(m, "C");
    }
    // Naphthalene ring A: 0,1,2,3,4,9
    indigoAddBond(at[0], at[1], 4);
    indigoAddBond(at[1], at[2], 4);
    indigoAddBond(at[2], at[3], 4);
    indigoAddBond(at[3], at[4], 4);
    indigoAddBond(at[4], at[9], 4);
    indigoAddBond(at[9], at[0], 4);
    // Naphthalene ring B: 4,5,6,7,8,9
    indigoAddBond(at[4], at[5], 4);
    indigoAddBond(at[5], at[6], 4);
    indigoAddBond(at[6], at[7], 4);
    indigoAddBond(at[7], at[8], 4);
    indigoAddBond(at[8], at[9], 4);
    // Furan 1 on bond 0-1 (1,2 of naphthalene, opposite O 11)
    indigoAddBond(at[0], at[10], 4);
    indigoAddBond(at[10], at[11], 4);
    indigoAddBond(at[11], at[12], 4);
    indigoAddBond(at[12], at[1], 4);
    // Furan 2 on bond 7-8 (7,8 of naphthalene, opposite O 14)
    indigoAddBond(at[7], at[13], 4);
    indigoAddBond(at[13], at[14], 4);
    indigoAddBond(at[14], at[15], 4);
    indigoAddBond(at[15], at[8], 4);

    std::cout << "  Constructed naphtho smiles: " << indigoSmiles(m) << "\n";

    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString desc = encodeFusionDescriptor(tree, ctx.g, m);
    EXPECT_EQ(desc.toStdString(), "[1,2-c:7,8-c']");
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "naphtho[1,2-c:7,8-c']difuran");
    indigoFree(m);
}

static void test_locant_omission_dibenzo_phenanthrene() {
    std::cout << "Test: LocantOmissionFR46 (dibenzo[c,g]phenanthrene)\n";
    int m = indigoLoadMoleculeFromString("C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=CC=CC=C54");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString desc = encodeFusionDescriptor(tree, ctx.g, m);
    EXPECT_EQ(desc.toStdString(), "[c,g]");
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "dibenzo[c,g]phenanthrene");
    indigoFree(m);
}

static void test_identical_components_difuro_pyridine() {
    std::cout << "Test: IdenticalComponentsDifuroPyridine (difuro[3,2-b:3',4'-e]pyridine)\n";
    int m = indigoCreateMolecule();
    int p[6];
    p[0] = indigoAddAtom(m, "N");
    for (int i = 1; i < 6; ++i) p[i] = indigoAddAtom(m, "C");
    // Pyridine ring: p0-p1-p2-p3-p4-p5-p0
    for (int i = 0; i < 6; ++i) indigoAddBond(p[i], p[(i + 1) % 6], 4);

    // Furan 1 on p1-p2 (side b), fusion 3,2-b
    // Furan atoms: 1=O, 2=C, 3=C, 4=C, 5=C.
    // p1 is 3, p2 is 2.
    int f1_o = indigoAddAtom(m, "O");
    int f1_c5 = indigoAddAtom(m, "C");
    int f1_c4 = indigoAddAtom(m, "C");
    indigoAddBond(p[2], f1_o, 4);    // 2-1
    indigoAddBond(f1_o, f1_c5, 4);   // 1-5
    indigoAddBond(f1_c5, f1_c4, 4);  // 5-4
    indigoAddBond(f1_c4, p[1], 4);   // 4-3

    // Furan 2 on p4-p5 (side e), fusion 3',4'-e
    // Furan atoms: 1'=O, 2'=C, 3'=C, 4'=C, 5'=C.
    // p4 is 3', p5 is 4'.
    int f2_c2 = indigoAddAtom(m, "C");
    int f2_o = indigoAddAtom(m, "O");
    int f2_c5 = indigoAddAtom(m, "C");
    indigoAddBond(p[4], f2_c2, 4);   // 3-2
    indigoAddBond(f2_c2, f2_o, 4);   // 2-1
    indigoAddBond(f2_o, f2_c5, 4);   // 1-5
    indigoAddBond(f2_c5, p[5], 4);   // 5-4

    std::cout << "  Constructed difuro smiles: " << indigoSmiles(m) << "\n";

    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "difuro[3,2-b:3',4'-e]pyridine");
    indigoFree(m);
}

static void test_higher_order_punctuation_dithieno_cyclohepta_imidazole() {
    std::cout << "Test: HigherOrderPunctuationFR49 (dithieno[2',3':3,4;2\",3\":6,7]cyclohepta[1,2-d]imidazole)\n";
    int m = indigoCreateMolecule();
    // Imidazole
    int im[5];
    im[0] = indigoAddAtom(m, "N");
    im[1] = indigoAddAtom(m, "C");
    im[2] = indigoAddAtom(m, "N");
    im[3] = indigoAddAtom(m, "C");
    im[4] = indigoAddAtom(m, "C");
    for (int i = 0; i < 5; ++i) indigoAddBond(im[i], im[(i + 1) % 5], 1);

    // Cyclohepta (7 atoms: c1=im[3], c2=im[4], c3..c7)
    int c[8];
    c[1] = im[3];
    c[2] = im[4];
    for (int i = 3; i <= 7; ++i) c[i] = indigoAddAtom(m, "C");
    for (int i = 2; i <= 6; ++i) indigoAddBond(c[i], c[i + 1], 1);
    indigoAddBond(c[7], c[1], 1);

    // Thiophene 1 on c[3]-c[4] (fusion 2',3':3,4)
    // S is 1', c[3] is 2', c[4] is 3'
    int th1_s = indigoAddAtom(m, "S");
    int th1_c5 = indigoAddAtom(m, "C");
    int th1_c4 = indigoAddAtom(m, "C");
    indigoAddBond(c[3], th1_s, 1);   // 2'-1'
    indigoAddBond(th1_s, th1_c5, 1); // 1'-5'
    indigoAddBond(th1_c5, th1_c4, 1);// 5'-4'
    indigoAddBond(th1_c4, c[4], 1);  // 4'-3'

    // Thiophene 2 on c[6]-c[7] (fusion 2",3":6,7)
    // S is 1", c[6] is 2", c[7] is 3"
    int th2_s = indigoAddAtom(m, "S");
    int th2_c5 = indigoAddAtom(m, "C");
    int th2_c4 = indigoAddAtom(m, "C");
    indigoAddBond(c[6], th2_s, 1);   // 2"-1"
    indigoAddBond(th2_s, th2_c5, 1); // 1"-5"
    indigoAddBond(th2_c5, th2_c4, 1);// 5"-4"
    indigoAddBond(th2_c4, c[7], 1);  // 4"-3"

    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "dithieno[2',3':3,4;2'',3'':6,7]cyclohepta[1,2-d]imidazole");
    indigoFree(m);
}

static void test_separate_citation_fr491() {
    std::cout << "Test: SeparateCitationFR491 (furo[3,4-b]furo[3',2':4,5]furo[2,3-e]pyridine)\n";
    int m = indigoCreateMolecule();
    int p[6];
    p[0] = indigoAddAtom(m, "N");
    for (int i = 1; i < 6; ++i) p[i] = indigoAddAtom(m, "C");
    for (int i = 0; i < 6; ++i) indigoAddBond(p[i], p[(i + 1) % 6], 1);

    // Furan 1 on p1-p2 (side b): furo[3,4-b]
    // 3 is p1, 4 is p2
    int f1_c2 = indigoAddAtom(m, "C");
    int f1_o = indigoAddAtom(m, "O");
    int f1_c5 = indigoAddAtom(m, "C");
    indigoAddBond(p[1], f1_c2, 1);  // 3-2
    indigoAddBond(f1_c2, f1_o, 1);  // 2-1
    indigoAddBond(f1_o, f1_c5, 1);  // 1-5
    indigoAddBond(f1_c5, p[2], 1);  // 5-4

    // Furan 2 on p4-p5 (side e): furo[2,3-e]
    // 2 is p4, 3 is p5
    int f2_o = indigoAddAtom(m, "O");
    int f2_c5 = indigoAddAtom(m, "C");
    int f2_c4 = indigoAddAtom(m, "C");
    indigoAddBond(p[4], f2_o, 1);   // 2-1
    indigoAddBond(f2_o, f2_c5, 1);  // 1-5
    indigoAddBond(f2_c5, f2_c4, 1); // 5-4
    indigoAddBond(f2_c4, p[5], 1);  // 4-3

    // Furan 3 on f2_c4 - f2_c5 (bond 4-5 of Furan 2): furo[3',2':4,5]
    // 3' is f2_c4, 2' is f2_c5
    int f3_o = indigoAddAtom(m, "O");
    int f3_c5 = indigoAddAtom(m, "C");
    int f3_c4 = indigoAddAtom(m, "C");
    indigoAddBond(f2_c5, f3_o, 1);  // 2'-1'
    indigoAddBond(f3_o, f3_c5, 1);  // 1'-5'
    indigoAddBond(f3_c5, f3_c4, 1); // 5'-4'
    indigoAddBond(f3_c4, f2_c4, 1); // 4'-3'

    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "furo[3,4-b]furo[3',2':4,5]furo[2,3-e]pyridine");
    indigoFree(m);
}

static void test_composite_grouping_bis_fr492() {
    std::cout << "Test: CompositeGroupingBisFR492 (bis(pyrimido[5',4':4,5]pyrrolo)[2,3-c:3',2'-e]pyridazine)\n";
    int m = indigoCreateMolecule();
    int pd[6];
    pd[0] = indigoAddAtom(m, "N");
    pd[1] = indigoAddAtom(m, "N");
    for (int i = 2; i < 6; ++i) pd[i] = indigoAddAtom(m, "C");
    for (int i = 0; i < 6; ++i) indigoAddBond(pd[i], pd[(i + 1) % 6], 1);

    // Branch 1 on pd2-pd3 (side c): pyrimido[5',4':4,5]pyrrolo[2,3-c]
    // Pyrrole 1 on pd2-pd3 (atoms 2,3 of pyrrole)
    int py1_n = indigoAddAtom(m, "N");
    int py1_c5 = indigoAddAtom(m, "C");
    int py1_c4 = indigoAddAtom(m, "C");
    indigoAddBond(pd[2], py1_n, 1);   // 2-1
    indigoAddBond(py1_n, py1_c5, 1);  // 1-5
    indigoAddBond(py1_c5, py1_c4, 1); // 5-4
    indigoAddBond(py1_c4, pd[3], 1);  // 4-3

    // Pyrimidine 1 on py1_c4 - py1_c5 (atoms 5',4' of pyrimidine)
    int pm1_c6 = indigoAddAtom(m, "C");
    int pm1_n1 = indigoAddAtom(m, "N");
    int pm1_c2 = indigoAddAtom(m, "C");
    int pm1_n3 = indigoAddAtom(m, "N");
    indigoAddBond(py1_c4, pm1_c6, 1); // 5'-6'
    indigoAddBond(pm1_c6, pm1_n1, 1); // 6'-1'
    indigoAddBond(pm1_n1, pm1_c2, 1); // 1'-2'
    indigoAddBond(pm1_c2, pm1_n3, 1); // 2'-3'
    indigoAddBond(pm1_n3, py1_c5, 1); // 3'-4'

    // Branch 2 on pd4-pd5 (side e): pyrimido[5',4':4,5]pyrrolo[3',2'-e]
    // Pyrrole 2 on pd4-pd5 (atoms 3',2' of pyrrole)
    int py2_n = indigoAddAtom(m, "N");
    int py2_c5 = indigoAddAtom(m, "C");
    int py2_c4 = indigoAddAtom(m, "C");
    indigoAddBond(pd[5], py2_n, 1);   // 2'-1'
    indigoAddBond(py2_n, py2_c5, 1);  // 1'-5'
    indigoAddBond(py2_c5, py2_c4, 1); // 5'-4'
    indigoAddBond(py2_c4, pd[4], 1);  // 4'-3'

    // Pyrimidine 2 on py2_c4 - py2_c5
    int pm2_c6 = indigoAddAtom(m, "C");
    int pm2_n1 = indigoAddAtom(m, "N");
    int pm2_c2 = indigoAddAtom(m, "C");
    int pm2_n3 = indigoAddAtom(m, "N");
    indigoAddBond(py2_c4, pm2_c6, 1);
    indigoAddBond(pm2_c6, pm2_n1, 1);
    indigoAddBond(pm2_n1, pm2_c2, 1);
    indigoAddBond(pm2_c2, pm2_n3, 1);
    indigoAddBond(pm2_n3, py2_c5, 1);

    indigoAromatize(m);
    std::cout << "  dithieno smiles: " << indigoSmiles(m) << "\n";
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);

    EXPECT_TRUE(tree.valid);
    QString name = constructFusionName(tree, ctx.g, m);
    EXPECT_EQ(name.toStdString(), "bis(pyrimido[5',4':4,5]pyrrolo)[2,3-c:3',2'-e]pyridazine");
    indigoFree(m);
}

static void test_skeletal_replacement_hexaaza_acepentalene() {
    std::cout << "Test: SkeletalReplacementAlternative1 (1,2,3,4,5,6-hexaazacyclopenta[cd]pentalene)\n";
    int m = indigoLoadMoleculeFromString("N1=NC2=NN=C3C2=C1N=N3");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);
    EXPECT_TRUE(hitsFr32Limitation(tree));

    QString name = nameBySkeletalReplacementFallback(tree, ctx.g, ctx.sssrRings, m);
    EXPECT_EQ(name.toStdString(), "1,2,3,4,5,6-hexaazacyclopenta[cd]pentalene");
    indigoFree(m);
}

static void test_skeletal_replacement_tetraaza_phenalene() {
    std::cout << "Test: SkeletalReplacementAlternative1 (1,4,9-tetraazaphenalene)\n";
    int m = indigoLoadMoleculeFromString("n1ccc2nccc3ccnc1[n]23");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);
    EXPECT_TRUE(hitsFr32Limitation(tree));

    QString name = nameBySkeletalReplacementFallback(tree, ctx.g, ctx.sssrRings, m);
    EXPECT_TRUE(name == "1,4,9-tetraazaphenalene" || name == "1,4,9,9b-tetraazaphenalene");
    indigoFree(m);
}

static void test_alternative_2_cyclobuta_indeno_naphthalene() {
    std::cout << "Test: Alternative 2 (cyclobuta[1,7]indeno[5,6-b]naphthalene)\n";
    int m = indigoLoadMoleculeFromString("c5cccc4c5cc3c(c4)c2c1c(c3)ccc1c2");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);
    EXPECT_TRUE(hitsFr32Limitation(tree));

    QString name = nameByLessPreferredParentFallback(ctx.g, ctx.sssrRings, in, m);
    std::cout << "ALTERNATIVE 2 NAME: '" << name.toStdString() << "'" << std::endl;
    EXPECT_EQ(name.toStdString(), "cyclobuta[1,7]indeno[5,6-b]naphthalene");
    indigoFree(m);
}

static void test_alternative_3_etheno_benzo_cyclohepta_naphthalene() {
    std::cout << "Test: Alternative 3 (1,12-ethenobenzo[4,5]cyclohepta[1,2,3-de]naphthalene)\n";
    int m = indigoLoadMoleculeFromString("C1c2cccc3ccc4ccc5c(c4c23)c(ccc5)C=1");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);

    QString name = nameByBridgedFusedFallback(ctx.g, ctx.sssrRings, in, m);
    std::cout << "ALTERNATIVE 3 NAME: '" << name.toStdString() << "'" << std::endl;
    EXPECT_EQ(name.toStdString(), "1,12-ethenobenzo[4,5]cyclohepta[1,2,3-de]naphthalene");
    indigoFree(m);

    // Test seed 3 reordering
    int m3 = indigoLoadMoleculeFromString("c1cc2c3c(c1)C=Cc1cccc4c1c3c(cc4)cc2");
    if (m3 >= 0) {
        indigoAromatize(m3);
        auto ctx3 = loadContext(m3);
        FusedRingSystemInput in3 = detectFusedRingDirections(m3);
        QString name3 = nameByBridgedFusedFallback(ctx3.g, ctx3.sssrRings, in3, m3);
        EXPECT_EQ(name3.toStdString(), "1,12-ethenobenzo[4,5]cyclohepta[1,2,3-de]naphthalene");
        indigoFree(m3);
    }
}

static void test_strict_ordering_alt1_preferred_over_alt2_and_alt3() {
    std::cout << "Test: Strict ordering: Alternative 1 preferred over Alternative 2 and 3\n";
    int m = indigoLoadMoleculeFromString("N1=NC2=NN=C3C2=C1N=N3");
    EXPECT_TRUE(m >= 0);
    if (m < 0) return;
    indigoAromatize(m);
    auto ctx = loadContext(m);
    FusedRingSystemInput in = detectFusedRingDirections(m);
    ComponentTree tree = decomposeComponentTree(ctx.g, ctx.sssrRings, in, m);
    EXPECT_TRUE(hitsFr32Limitation(tree));

    // Alternative 1 matches and must take priority
    QString name1 = nameBySkeletalReplacementFallback(tree, ctx.g, ctx.sssrRings, m);
    EXPECT_FALSE(name1.isEmpty());
    EXPECT_EQ(name1.toStdString(), "1,2,3,4,5,6-hexaazacyclopenta[cd]pentalene");
    indigoFree(m);
}

int main() {
    unsigned long long sid = indigoAllocSessionId();
    indigoSetSessionId(sid);

    test_hits_fr32_synthetic();
    test_decomposes_clean_three_component_chain();
    test_detects_fr32_limitation_on_mutually_fused_pair();
    test_parent_selection_determinism();
    test_citation_order_fr36();
    test_benzoheterocycle_merging_fr35();
    test_benzoheterocycle_merging_skipped_on_multiplicative();
    test_descriptor_encoding_dipyrido_pyrazine();
    test_descriptor_encoding_naphtho_difuran();
    test_locant_omission_dibenzo_phenanthrene();
    test_identical_components_difuro_pyridine();
    test_higher_order_punctuation_dithieno_cyclohepta_imidazole();
    test_separate_citation_fr491();
    test_composite_grouping_bis_fr492();
    test_skeletal_replacement_hexaaza_acepentalene();
    test_skeletal_replacement_tetraaza_phenalene();
    test_alternative_2_cyclobuta_indeno_naphthalene();
    test_alternative_3_etheno_benzo_cyclohepta_naphthalene();
    test_strict_ordering_alt1_preferred_over_alt2_and_alt3();

    indigoReleaseSessionId(sid);


    std::cout << "\nTOTAL PASSED: " << passed << "\n";
    std::cout << "TOTAL FAILED: " << failed << "\n";
    return (failed == 0) ? 0 : 1;
}
