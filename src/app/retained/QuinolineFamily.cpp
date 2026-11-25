#include "QuinolineFamily.h"

// All entries use the naphthalene 6-6 ortho-fused bicyclic skeleton with the
// same peripheral cycle order: 1-2-3-4-4a-5-6-7-8-8a-1.
//
// Reference naphthalene (from HydrocarbonFamily.cpp):
//   SMILES: c1ccc2ccccc2c1
//   atom[0] -> "2", atom[1] -> "3", atom[2] -> "4",
//   atom[3] -> "4a", atom[4] -> "5", atom[5] -> "6",
//   atom[6] -> "7", atom[7] -> "8", atom[8] -> "8a", atom[9] -> "1"
//   Bridgeheads at atom[3] ("4a") and atom[8] ("8a"), both deg=3.
//
// All heterobicyclics below reuse this same atom->locant mapping.
// Each SMILES was empirically verified to place N atoms at the correct
// atom indices matching their chemical locant positions.
//
// RETROACTIVE FIX (found independently during Task 8's matching-engine
// verification, long after this file's original commit 35f1a4c): every
// entry below except 1,5-naphthyridine used an UPPERCASE "N" for its
// pyridine-type (no N-H) nitrogen atoms mixed into an otherwise lowercase-
// aromatic SMILES string (e.g. isoquinoline's "N1ccc2ccccc2c1"). Indigo's
// aromaticity perception does not treat that uppercase N as part of the
// aromatic ring even after an explicit indigoAromatize() call -- it leaves
// the two ring bonds adjacent to it as non-aromatic/open, which silently
// gives the WRONG gross molecular formula (confirmed via indigoGrossFormula:
// e.g. isoquinoline parsed to C9H11N instead of the real C9H7N; quinazoline/
// quinoxaline/cinnoline/phthalazine/4 of the 5 naphthyridines all showed the
// same class of error). This predates bug class 4's discovery (Task 6, long
// after this file was committed in Task 2) so gross formula was never
// checked against these entries until now. 1,5-naphthyridine was ALREADY
// using lowercase "n" for both its nitrogens and already gave the correct
// C8H6N2 -- that one entry is the reason this bug wasn't caught by a
// same-file comparison sooner.
//
// FIX: every pyridine-type N below (none of these systems have an N-H --
// all are 2-ring-bond, lone-pair nitrogens, same valence class as the
// carbons they replace) is now written as lowercase "n", matching
// 1,5-naphthyridine's already-correct pattern. Confirmed via fresh Indigo
// runs: atom order, element, and degree at every index are IDENTICAL
// between the old uppercase and new lowercase SMILES (case-only change) --
// only the computed H-count differs, now correctly C9H7N (mono-aza) or
// C8H6N2 (diaza) for every entry, matching the real compounds' known
// formulas. No locant/bondLetter changes needed for any entry.

void registerQuinolineFamily(std::vector<RetainedRingSystem> &registry) {
    // 1. quinoline: N at locant 1 (atom[9])
    // SMILES c1ccc2ccccc2N1 verified: atom[9]=N, all others C, 10 atoms
    // Locants: atom[0]="2", atom[1]="3", atom[2]="4", atom[3]="4a",
    //          atom[4]="5", atom[5]="6", atom[6]="7", atom[7]="8",
    //          atom[8]="8a", atom[9]="1" (N)
    registry.push_back(RetainedRingSystem{
        "quinoline",
        "c1ccc2ccccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        34, // Appendix 2 seniority rank
    });

    // 2. isoquinoline: N at locant 2 (atom[0])
    // SMILES N1ccc2ccccc2c1 verified: atom[0]=N, all others C, 10 atoms
    // Locants: same mapping as quinoline
    registry.push_back(RetainedRingSystem{
        "isoquinoline",
        "n1ccc2ccccc2c1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        35, // Appendix 2 seniority rank
    });

    // 3. cinnoline: N at locants 1 and 2 (atoms[9] and atom[0])
    // SMILES N1ccc2ccccc2N1 verified: atom[0]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "cinnoline",
        "n1ccc2ccccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        24, // Appendix 2 seniority rank
    });

    // 4. quinazoline: N at locants 1 and 3.
    // CORRECTED from the originally-dispatched "c1Nc2ccccc2cN1", which
    // opened ring2 one token position too early (at index 2 instead of the
    // reference skeleton's index 3), shifting the whole backbone and
    // actually placing N at locants 6 and 8 instead of the required 1 and
    // 3 -- caught by independently re-tracing the real bond topology
    // (ring-closure token positions), not by validateRetainedRingSystem
    // (which only checks internal list self-consistency, not that the
    // atom order matches the intended chemistry). Verified: SMILES
    // "c1Ncc2ccccc2N1" has 10 atoms, bridgeheads at indices 3 and 8
    // (matching the shared naphthalene-skeleton reference exactly), N at
    // index1 (-> locant "3") and index9 (-> locant "1"): N(1)-C(2)-N(3)-
    // C(4)-C(4a)-... -- the textbook quinazoline arrangement.
    registry.push_back(RetainedRingSystem{
        "quinazoline",
        "c1ncc2ccccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        25, // Appendix 2 seniority rank
    });

    // 5. quinoxaline: N at locants 1 and 4 (atoms[9] and atom[2])
    // SMILES c1cNc2ccccc2N1 verified: atom[2]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "quinoxaline",
        "c1cnc2ccccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        26, // Appendix 2 seniority rank
    });

    // 6. phthalazine: N at locants 2 and 3 (atoms[0] and atom[1])
    // SMILES N1Ncc2ccccc2c1 verified: atom[0]=N, atom[1]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "phthalazine",
        "n1ncc2ccccc2c1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        31, // Appendix 2 seniority rank
    });

    // 7. 1,5-naphthyridine: N at locants 1 and 5 (atoms[9] and atom[4])
    // SMILES c1ccc2ncccc2n1 verified: atom[4]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "1,5-naphthyridine",
        "c1ccc2ncccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        27, // Appendix 2 seniority rank
    });

    // 8. 1,6-naphthyridine: N at locants 1 and 6 (atoms[9] and atom[5])
    // SMILES c1ccc2cNccc2N1 verified: atom[5]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "1,6-naphthyridine",
        "c1ccc2cnccc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        28, // Appendix 2 seniority rank
    });

    // 9. 1,7-naphthyridine: N at locants 1 and 7 (atoms[9] and atom[6])
    // SMILES c1ccc2ccNcc2N1 verified: atom[6]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "1,7-naphthyridine",
        "c1ccc2ccncc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        29, // Appendix 2 seniority rank
    });

    // 10. 1,8-naphthyridine: N at locants 1 and 8 (atoms[9] and atom[7])
    // SMILES c1ccc2cccNc2N1 verified: atom[7]=N, atom[9]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "1,8-naphthyridine",
        "c1ccc2cccnc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        30, // Appendix 2 seniority rank
    });

    // 11. 2,6-naphthyridine: N at locants 2 and 6 (atoms[0] and atom[5])
    // SMILES N1ccc2cNccc2c1 verified: atom[0]=N, atom[5]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "2,6-naphthyridine",
        "n1ccc2cnccc2c1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        32, // Appendix 2 seniority rank
    });

    // 12. 2,7-naphthyridine: N at locants 2 and 7 (atoms[0] and atom[6])
    // SMILES N1ccc2ccNcc2c1 verified: atom[0]=N, atom[6]=N, 10 atoms
    registry.push_back(RetainedRingSystem{
        "2,7-naphthyridine",
        "n1ccc2ccncc2c1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        33, // Appendix 2 seniority rank
    });

    // 13. pteridine: BICYCLIC (pyrazino[2,3-d]pyrimidine, rank 23)
    // Moved from TricyclicNFamily.cpp (where it was incorrectly categorized)
    // to here, since pteridine is a bicyclic 6-6 ortho-fused system, not
    // tricyclic. Appendix 2 entry 23 (https://iupac.qmul.ac.uk/fusedring/app2.html).
    //
    // CORRECTED: a prior pass placed N at locants 1,2,4,7 and explicitly
    // noted in its own comment that this did NOT match the real pteridine
    // positions (1,3,5,8) -- it passed validateRetainedRingSystem's cyclic-
    // adjacency self-test (which only checks internal array consistency,
    // not real chemistry) and was committed anyway. Real pteridine N
    // positions (N1,C2,N3,C4,C4a,N5,C6,C7,N8,C8a -- pyrimidine ring N1/N3,
    // pyrazine ring N5/N8) confirmed via biochemistry-literature naming
    // convention ("2-aminopteridin-4(3H)-one" places substituents at C2/C4,
    // confirming those are carbon, not nitrogen, ruling out the prior
    // pass's N-at-2 placement).
    //
    // Using this file's own established naphthalene atom->locant mapping
    // (see header comment: atom[0]="2", atom[1]="3", atom[2]="4",
    // atom[3]="4a", atom[4]="5", atom[5]="6", atom[6]="7", atom[7]="8",
    // atom[8]="8a", atom[9]="1"), N at real locants {1,3,5,8} maps to atom
    // indices {9,1,4,7}. SMILES "c1ncc2nccnc2n1" (lowercase 'n' throughout,
    // per this file's own RETROACTIVE FIX note above -- uppercase N breaks
    // aromaticity perception). Verified via dump_rings.exe: 10 atoms, N at
    // exactly atoms 1,4,7,9, bridgeheads atom[3]/atom[8] remain C degree=3,
    // same two SSSR ring sets as every other entry in this file. Gross
    // formula C6H4N4, matching real pteridine (MW 132.12).
    registry.push_back(RetainedRingSystem{
        "pteridine",
        "c1ncc2nccnc2n1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        23, // Appendix 2 seniority rank
    });

    // 14. purine: BICYCLIC (9H-purine, rank 21)
    // 5-6 ortho-fused imidazole[4,5-d]pyrimidine system.
    // Appendix 2 entry 21. Standard IUPAC purine peripheral locants:
    // N1, C2, N3, C4(bridgehead), N9(indicated H), C8, N7, C5(bridgehead), C6.
    // SMILES: n1cnc2[nH]cnc2c1
    // Atom mapping:
    // atom[0]=1 (N), atom[1]=2 (C), atom[2]=3 (N), atom[3]=4 (C, bh)
    // atom[4]=9 (N), atom[5]=8 (C), atom[6]=7 (N), atom[7]=5 (C, bh), atom[8]=6 (C)
    // Cyclic perimeter sequence: {"1","2","3","4","9","8","7","5","6"}
    // Bond letters: 'a':1-2, 'b':2-3, 'c':3-4, 'd':4-9, 'e':9-8, 'f':8-7, 'g':7-5, 'h':5-6, 'i':6-1
    registry.push_back(RetainedRingSystem{
        "purine",
        "n1cnc2[nH]cnc2c1",
        {"1","2","3","4","9","8","7","5","6"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","9"}}, {'e', {"9","8"}}, {'f', {"8","7"}},
            {'g', {"7","5"}}, {'h', {"5","6"}}, {'i', {"6","1"}},
        },
        21, // Appendix 2 seniority rank
        "purino", // fusion prefix
        "", // indicated hydrogen locant: IUPAC PIN bare name is "purine"
    });

    // 15. quinolizine: Bridgehead N - SKIPPED
}
