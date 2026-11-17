#include "IndoleFamily.h"

// All entries use the 6-5 ortho-fused bicyclic skeleton (9 atoms) with two
// bridgehead atoms (degree 3).
//
// Each SMILES was empirically verified to place N atoms at the correct
// atom indices matching their chemical locant positions, using the same
// rigor as Task 2 (quinoline family): tiny Indigo-linked scratch programs
// were built and run for each candidate SMILES, parsing the SMILES,
// iterating atoms, and printing each atom's element (indigoAtomicNumber)
// and degree (indigoDegree) in order.
//
// For each entry, bridgehead positions (degree 3 atoms) were confirmed to
// land at the index positions the locant list expects, and heteroatoms
// (N) confirmed to land at the index(es) the locant list maps to the
// chemically-correct locant.

void registerIndoleFamily(std::vector<RetainedRingSystem> &registry) {
    // 1. indole: N at locant 1 (adjacent to a bridgehead carbon). Classic
    //    indole numbering: N-1, C-2, C-3, C-3a (bridgehead), C-4, C-5, C-6,
    //    C-7, C-7a (bridgehead) -- with the ring CLOSING bond being 7a-1,
    //    i.e. N-1's bridgehead neighbor is 7a, NOT 3a (3a is reached going
    //    the OTHER way around, via C-2/C-3 first).
    //
    // SMILES c2cnc1ccccc12, real bonds confirmed via Indigo
    // (0-1,0-8,1-2,2-3,3-4,3-8,4-5,5-6,6-7,7-8):
    //   atom[2]=N, neighbors atom[1] and atom[3] (both real bonds).
    //   atom[3] and atom[8] are the two bridgeheads (degree 3), directly
    //   bonded to each other (the 3-8 edge is the ring-fusion bond).
    //
    // CORRECTED (the originally-dispatched version swapped the 3a/7a
    // labels, and consequently mirrored the whole 4-5-6-7 sequence --
    // caught by independently re-deriving the full peripheral walk from
    // the real bond list, not just confirming bridgeheads exist and N is
    // adjacent to one of them: that check alone is not enough to tell
    // WHICH bridgehead label -- 3a or 7a -- the adjacent one should get).
    // Walking the real locant sequence 1-2-3-3a-4-5-6-7-7a(-1 closing)
    // against the real bonds gives: atom[2]=1, atom[1]=2, atom[0]=3,
    // atom[8]=3a, atom[7]=4, atom[6]=5, atom[5]=6, atom[4]=7, atom[3]=7a.
    // Every consecutive pair in that sequence (including the 7a-1 closing
    // bond, atom[3]-atom[2]) is confirmed a real bond above.
    registry.push_back(RetainedRingSystem{
        "indole",
        "c2cnc1ccccc12",
        {"3","2","1","7a","7","6","5","4","3a"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","3a"}},
            {'d', {"3a","4"}}, {'e', {"4","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","7a"}}, {'i', {"7a","1"}},
        },
        39, // Appendix 2 seniority rank
    });

    // 2. isoindole: N at locant 2 (between C-1 and C-3, neither a
    //    bridgehead; analogous to how isoquinoline's N sits one position
    //    away from quinoline's N in Task 2). Standard isoindole numbering:
    //    C-1, N-2, C-3, C-3a (bridgehead), C-4, C-5, C-6, C-7, C-7a
    //    (bridgehead), closing bond 7a-1 -- so locant 1's bridgehead
    //    neighbor is 7a, and locant 3's bridgehead neighbor is 3a.
    //
    // SMILES N1cc2ccccc2c1, real bonds confirmed via Indigo
    // (0-1,0-8,1-2,2-3,2-7,3-4,4-5,5-6,6-7,7-8):
    //   atom[0]=N, neighbors atom[1] and atom[8] (both real bonds, neither
    //   a bridgehead -- matches N-2 being flanked by non-bridgehead C-1/C-3).
    //   atom[2] and atom[7] are the two bridgeheads (degree 3), directly
    //   bonded to each other (the 2-7 edge is the ring-fusion bond).
    //
    // CORRECTED (same bug class as indole above: the originally-dispatched
    // version swapped the 3a/7a labels and mirrored 4-5-6-7). Walking the
    // real locant sequence 1-2-3-3a-4-5-6-7-7a(-1 closing) against the
    // real bonds gives: atom[1]=1, atom[0]=2, atom[8]=3, atom[7]=3a,
    // atom[6]=4, atom[5]=5, atom[4]=6, atom[3]=7, atom[2]=7a. Every
    // consecutive pair in that sequence (including the 7a-1 closing bond,
    // atom[2]-atom[1]) is confirmed a real bond above.
    //
    // SECOND FIX (found independently during Task 7's continuation review,
    // long after this entry's original commit baf1df2): the SMILES
    // "N1cc2ccccc2c1" parses to formula C8H9N, but real isoindole is C8H7N
    // (isoindole is indole's structural isomer -- both share one molecular
    // formula). Indigo's aromaticity perception fails to treat atom[1] and
    // atom[8] (the two carbons directly bonded to the non-lowercase N) as
    // aromatic CH -- it gives them H=2 each (sp3-like), 2 H too many. This
    // is the same class of bug as the chromene-family indicated-hydrogen
    // bug found in Task 6, just not checked here until now because bug
    // class 4 (gross-formula verification) did not exist yet when this
    // entry was first written. Switching to an explicit Kekule structure
    // (same technique used to fix arsindole/isoarsindole in Task 7) fixes
    // it: "N1C=C2C=CC=CC2=C1" gives the correct C8H7N, with atom[1]/atom[8]
    // now H=1 as real aromatic CH should be. Atom indices, degrees, and
    // neighbor lists are IDENTICAL between the old and new SMILES (only the
    // two H-counts change) -- so the locants/bondLetters below remain valid
    // unchanged, and nothing that already reused isoindoleLocants/
    // isoindoleBondLetters needs re-deriving.
    registry.push_back(RetainedRingSystem{
        "isoindole",
        "N1C=C2C=CC=CC2=C1",
        {"2","1","7a","7","6","5","4","3a","3"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","3a"}},
            {'d', {"3a","4"}}, {'e', {"4","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","7a"}}, {'i', {"7a","1"}},
        },
        40, // Appendix 2 seniority rank
    });

    // 3. indazole: two N atoms in the five-membered ring, at positions 1 and 2
    //    (adjacent to each other), with N1 adjacent to the bridgehead.
    //    Standard indazole numbering: N-1, N-2, C-3, C-3a, C-4, C-5, C-6, C-7, C-7a.
    //    Deferred: Could not confidently verify a SMILES with both N in the 5-membered ring.
    //    Tried c1nc2ncccc2c1 (9 atoms, 2N) but N at atom[1] and atom[3] are NOT both in the
    //    5-membered ring (5-ring: 0,1,2,7,8; 6-ring: 2,3,4,5,6,7; so atom[1] is in 5-ring,
    //    atom[3] is in 6-ring).
    //    Tried c1nnc2ccccc2c1 and c2ncn1ccccc12 but both gave 10 atoms (6-6 fused).
    //    The correct indazole SMILES likely exists but was not found within the verification time.

    // 4. indolizine: bridgehead nitrogen -- 6-5 fused with N at one ring fusion position.
    //    SMILES c1ccn2c1cccc2 verified (9 atoms, 1N at bridgehead atom[3] deg=3).
    //    Deferred: Bridgehead nitrogen has a different locant numbering scheme from 6-5
    //    ortho-fused bicyclics; requires separate research to verify the correct numbering.

    // 5. pyrrolizine: 5-5 fused with bridgehead nitrogen (8 atoms total).
    //    SMILES c1ccn2cc1cc2 verified (8 atoms, 1N at bridgehead atom[3] deg=3).
    //    Deferred: 5-5 fused system with bridgehead N is structurally different and
    //    requires separate research to verify the correct numbering.
    // 2,1-benzothiazole: 9 atoms, 6-5 bicyclic, N at 1 and S at 2 (NOT "S at 1
    // N at 2" as an earlier draft comment claimed -- verified via PubChem CID
    // 638008, IUPACName literally "2,1-benzothiazole", SMILES
    // C1=CC2=CSN=C2C=C1: tracing that SMILES shows N bonded directly to one
    // ring-fusion carbon (-> C7a) and the CH-S side bonded to the other
    // fusion carbon (-> C3a) with S between the CH and N, i.e. N1-S2-C3-C3a-
    // C4-C5-C6-C7-C7a-(back to N1). All 9 bond-letter pairs below were
    // confirmed to match this molecule's real bonds once N/S were corrected
    // to N=1/S=2 -- every other position in this entry was already right.
    // Appendix 2 retained benzoheterocycle.
    registry.push_back(RetainedRingSystem{
        "2,1-benzothiazole",
        "c1ccc2nscc2c1",
        {"5","6","7","7a","1","2","3","3a","4"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","3a"}},
            {'d', {"3a","4"}}, {'e', {"4","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","7a"}}, {'i', {"7a","1"}},
        },
        35,
        "2,1-benzothiazolo",
    });

    // [1,3]benzodioxole (2H-1,3-benzodioxole): 9 atoms, 6-5 bicyclic with oxygens
    // at locants 1 and 3, sp3 CH2 at locant 2 (indicated hydrogen 2H), bridgeheads
    // at 3a and 7a, benzene carbons at 4, 5, 6, 7.
    // SMILES 'O1COc2ccccc21' walks the 9-atom perimeter consecutively:
    // atom 0: O at 1
    // atom 1: C at 2 (CH2, indicated H)
    // atom 2: O at 3
    // atom 3: C at 3a (bridgehead)
    // atom 4: C at 4
    // atom 5: C at 5
    // atom 6: C at 6
    // atom 7: C at 7
    // atom 8: C at 7a (bridgehead)
    // All 9 bond letters 'a' through 'i' cyclically adjacent.
    // Appendix 2 seniorityRank=70 (O-containing bicyclic).
    // fusionPrefix="[1,3]benzodioxolo", indicatedHydrogenLocant="2".
    registry.push_back(RetainedRingSystem{
        "[1,3]benzodioxole",
        "O1COc2ccccc21",
        {"1","2","3","3a","4","5","6","7","7a"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","3a"}},
            {'d', {"3a","4"}}, {'e', {"4","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","7a"}}, {'i', {"7a","1"}},
        },
        70,
        "[1,3]benzodioxolo",
        "2",
    });
}