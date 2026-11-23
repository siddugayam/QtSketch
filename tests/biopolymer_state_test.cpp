// tests/biopolymer_state_test.cpp
#include <QtTest/QtTest>
#include <QCoreApplication>

#include "BiopolymerState.h"
#include "EditableMolecule.h"

class BiopolymerStateTest : public QObject {
    Q_OBJECT

private slots:
    void testCreateChain() {
        BiopolymerState state;
        
        // Test adding a chain
        int chainId = state.addChain("A", PolymerType::Peptide);
        QCOMPARE(chainId, 1);
        QCOMPARE(state.chainCount(), 1);
        QCOMPARE(state.chainIds(), QList<int>{1});
        
        const BiopolymerChain* chain = state.chain(chainId);
        QVERIFY(chain != nullptr);
        QCOMPARE(chain->name, QString("A"));
        QCOMPARE(chain->polymerType, PolymerType::Peptide);
    }
    
    void testAddResidue() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";
        ala.monomerClass = "AminoAcid";
        ala.ambiguous = false;
        
        int residueId = state.addResidue(chainId, ala);
        QCOMPARE(residueId, 1);
        QCOMPARE(state.allResidueIds(), QList<int>{1});
        
        const BiopolymerResidue* residue = state.residue(residueId);
        QVERIFY(residue != nullptr);
        QCOMPARE(residue->identity.symbol, QString("A"));
        QCOMPARE(residue->identity.alias, QString("Ala"));
        QCOMPARE(residue->sequencePosition, 1);
    }
    
    void testResidueAtoms() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";
        
        int residueId = state.addResidue(chainId, ala);
        
        // Set atom ownership
        QList<AtomId> atomIds = {1, 2, 3, 4, 5};
        state.setResidueAtoms(residueId, atomIds);
        
        QCOMPARE(state.atomsOwnedByResidue(residueId), atomIds);
        
        // Check atom-to-residue mapping
        for (AtomId atomId : atomIds) {
            QCOMPARE(state.residueOwningAtom(atomId), residueId);
        }
        
        // Check non-existent atom
        QCOMPARE(state.residueOwningAtom(999), -1);
    }
    
    void testChainValidation() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";
        
        int residueId1 = state.addResidue(chainId, ala);
        int residueId2 = state.addResidue(chainId, ala);
        
        // Set atoms for both residues
        state.setResidueAtoms(residueId1, {1, 2, 3});
        state.setResidueAtoms(residueId2, {4, 5, 6});
        
        // Without connections, chain should not be fully valid
        MappingValidity validity = state.validateChain(chainId);
        QVERIFY(!validity.isValid);
        QVERIFY(validity.issues.join("\n").contains("No connection between residues"));
        
        // Add connection
        state.setResidueConnection(residueId1, residueId2, 1, "peptide");
        
        // Now should be valid (assuming atoms are correct)
        validity = state.validateChain(chainId);
        QVERIFY(validity.isValid);
    }
    
    void testDeleteResidue() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";
        
        int residueId1 = state.addResidue(chainId, ala);
        int residueId2 = state.addResidue(chainId, ala);
        
        QCOMPARE(state.allResidueIds().size(), 2);
        
        // Delete first residue
        QVERIFY(state.deleteResidue(residueId1));
        
        QCOMPARE(state.allResidueIds().size(), 1);
        QCOMPARE(state.allResidueIds().first(), residueId2);
        
        // Check that the remaining residue has the correct position
        const BiopolymerResidue* residue = state.residue(residueId2);
        QCOMPARE(residue->sequencePosition, 1);
    }
    
    void testInsertResidue() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";
        
        int residueId1 = state.addResidue(chainId, ala);
        int residueId3 = state.addResidue(chainId, ala);
        
        // Insert before residue 3
        MonomerIdentity val;
        val.symbol = "V";
        val.alias = "Val";
        
        int residueId2 = state.insertResidueBefore(chainId, residueId3, val);
        QCOMPARE(residueId2, 3);
        
        // Check positions
        QCOMPARE(state.residue(residueId1)->sequencePosition, 1);
        QCOMPARE(state.residue(residueId2)->sequencePosition, 2);
        QCOMPARE(state.residue(residueId3)->sequencePosition, 3);
        
        // Check order in chain
        QCOMPARE(state.residueIdsInChain(chainId), (QList<int>{residueId1, residueId2, residueId3}));
    }

    void testAddResidueBoundaryClamping() {
        BiopolymerState state;
        int chainId = state.addChain("A", PolymerType::Peptide);

        MonomerIdentity ala;
        ala.symbol = "A";
        ala.alias = "Ala";

        MonomerIdentity gly;
        gly.symbol = "G";
        gly.alias = "Gly";

        MonomerIdentity val;
        val.symbol = "V";
        val.alias = "Val";

        MonomerIdentity leu;
        leu.symbol = "L";
        leu.alias = "Leu";

        // Seed initial residue [r1]
        int r1 = state.addResidue(chainId, ala);
        // Position 0 clamps to front -> [rZero, r1]
        int rZero = state.addResidue(chainId, gly, 0);
        // Negative position clamps to front -> [rNeg, rZero, r1]
        int rNeg = state.addResidue(chainId, val, -5);
        // Out-of-range high position clamps to end -> [rNeg, rZero, r1, rBeyond]
        int rBeyond = state.addResidue(chainId, leu, 999);

        const QList<int> expectedOrder{rNeg, rZero, r1, rBeyond};
        QCOMPARE(state.residueIdsInChain(chainId), expectedOrder);
        for (int i = 0; i < expectedOrder.size(); ++i) {
            const BiopolymerResidue *res = state.residue(expectedOrder[i]);
            QVERIFY(res != nullptr);
            QCOMPARE(res->sequencePosition, i + 1);
        }
    }
    
    void testRenameChain() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        QCOMPARE(state.chain(chainId)->name, QString("A"));
        
        QVERIFY(state.renameChain(chainId, "B"));
        QCOMPARE(state.chain(chainId)->name, QString("B"));
        
        // Rename non-existent chain should fail
        QVERIFY(!state.renameChain(999, "C"));
    }
    
    void testRemoveChain() {
        BiopolymerState state;
        
        int chainId1 = state.addChain("A", PolymerType::Peptide);
        int chainId2 = state.addChain("B", PolymerType::Peptide);
        
        QCOMPARE(state.chainCount(), 2);
        
        QVERIFY(state.removeChain(chainId1));
        QCOMPARE(state.chainCount(), 1);
        QCOMPARE(state.chainIds().first(), chainId2);
        
        // Remove non-existent chain should fail
        QVERIFY(!state.removeChain(999));
    }
    
    void testPolymerTypeSupport() {
        QVERIFY(BiopolymerState::isPolymerTypeSupported(PolymerType::Peptide));
        QVERIFY(BiopolymerState::isPolymerTypeSupported(PolymerType::DNA));
        QVERIFY(BiopolymerState::isPolymerTypeSupported(PolymerType::RNA));
        QVERIFY(!BiopolymerState::isPolymerTypeSupported(PolymerType::Unknown));
    }
    
    void testPolymerTypeConversion() {
        QCOMPARE(BiopolymerState::polymerTypeToString(PolymerType::Peptide), QString("PEPTIDE"));
        QCOMPARE(BiopolymerState::polymerTypeToString(PolymerType::DNA), QString("DNA"));
        QCOMPARE(BiopolymerState::polymerTypeToString(PolymerType::RNA), QString("RNA"));
        QCOMPARE(BiopolymerState::polymerTypeToString(PolymerType::Unknown), QString("UNKNOWN"));
        
        QCOMPARE(BiopolymerState::polymerTypeFromString("PEPTIDE"), PolymerType::Peptide);
        QCOMPARE(BiopolymerState::polymerTypeFromString("PROTEIN"), PolymerType::Peptide);
        QCOMPARE(BiopolymerState::polymerTypeFromString("DNA"), PolymerType::DNA);
        QCOMPARE(BiopolymerState::polymerTypeFromString("RNA"), PolymerType::RNA);
        QCOMPARE(BiopolymerState::polymerTypeFromString("UNKNOWN"), PolymerType::Unknown);
    }
    
    void testHasBiopolymerData() {
        BiopolymerState state1;
        QVERIFY(!state1.hasBiopolymerData());
        
        BiopolymerState state2;
        state2.addChain("A", PolymerType::Peptide);
        QVERIFY(state2.hasBiopolymerData());
    }
    
    void testSelectionSupport() {
        BiopolymerState state;
        
        int chainId = state.addChain("A", PolymerType::Peptide);
        
        MonomerIdentity ala;
        ala.symbol = "A";
        
        int residueId1 = state.addResidue(chainId, ala);
        int residueId2 = state.addResidue(chainId, ala);
        
        state.setResidueAtoms(residueId1, {1, 2, 3});
        state.setResidueAtoms(residueId2, {4, 5, 6});
        
        // Test residuesFromAtoms
        QSet<AtomId> atomSet = {1, 2, 4, 5};
        QSet<int> residueSet = state.residuesFromAtoms(atomSet);
        QCOMPARE(residueSet.size(), 2);
        QVERIFY(residueSet.contains(residueId1));
        QVERIFY(residueSet.contains(residueId2));
        
        // Test atomsFromResidues
        QSet<int> residueSet2 = {residueId1};
        QSet<AtomId> atomSet2 = state.atomsFromResidues(residueSet2);
        QCOMPARE(atomSet2.size(), 3);
        QVERIFY(atomSet2.contains(1));
        QVERIFY(atomSet2.contains(2));
        QVERIFY(atomSet2.contains(3));
    }
    
    void testIdCounters() {
        BiopolymerState state;
        
        // Initial counters should be at 1
        BiopolymerState::IdCounters counters = state.idCounters();
        QCOMPARE(counters.nextResidueId, 1);
        QCOMPARE(counters.nextChainId, 1);
        QCOMPARE(counters.nextConnectionId, 1);
        
        // Add a chain
        state.addChain("A", PolymerType::Peptide);
        counters = state.idCounters();
        QCOMPARE(counters.nextChainId, 2);
        
        // Add a residue
        MonomerIdentity ala;
        state.addResidue(1, ala);
        counters = state.idCounters();
        QCOMPARE(counters.nextResidueId, 2);
        
        // Add a connection
        state.setResidueConnection(1, 1, 1, "test");
        counters = state.idCounters();
        QCOMPARE(counters.nextConnectionId, 2);
    }
};

QTEST_MAIN(BiopolymerStateTest)
#include "biopolymer_state_test.moc"
