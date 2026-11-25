#include "HydrocarbonFamily.h"

void registerHydrocarbonFamily(std::vector<RetainedRingSystem> &registry) {
    // Naphthalene -- THE motivating entry (naphtho[2,3-g]isoquinoline).
    // SMILES 'c1ccc2ccccc2c1' parsed by Indigo produces atoms in order:
    // atom 0: C beta (neighbors 1, 9) -> locant "2"
    // atom 1: C beta (neighbors 0, 2) -> locant "3"
    // atom 2: C alpha (neighbors 1, 3) -> locant "4"
    // atom 3: C bridgehead (neighbors 2, 4, 8) -> locant "4a"
    // atom 4: C alpha (neighbors 3, 5) -> locant "5"
    // atom 5: C beta (neighbors 4, 6) -> locant "6"
    // atom 6: C beta (neighbors 5, 7) -> locant "7"
    // atom 7: C alpha (neighbors 6, 8) -> locant "8"
    // atom 8: C bridgehead (neighbors 3, 7, 9) -> locant "8a"
    // atom 9: C alpha (neighbors 0, 8) -> locant "1"
    registry.push_back(RetainedRingSystem{
        "naphthalene",
        "c1ccc2ccccc2c1",
        {"2","3","4","4a","5","6","7","8","8a","1"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","4"}},
            {'d', {"4","4a"}}, {'e', {"4a","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","8"}}, {'i', {"8","8a"}},
            {'j', {"8a","1"}},
        },
        56, // Appendix 1 seniority position
        "naphtho", // Task 10: fusion prefix, for citing naphthalene as an
                   // attached component (e.g. "naphtho[2,3-g]isoquinoline")
    });

    // Anthracene -- 14 carbons; 9 and 10 are the two unique "meso" positions.
    // SMILES 'c1ccc2cc3ccccc3cc2c1' parsed by Indigo produces atoms in order:
    // atom 0: C beta (neighbors 1, 13) -> locant "2"
    // atom 1: C beta (neighbors 0, 2) -> locant "3"
    // atom 2: C alpha (neighbors 1, 3) -> locant "4"
    // atom 3: C bridgehead (neighbors 2, 4, 12) -> locant "4a"
    // atom 4: C meso (neighbors 3, 5) -> locant "10"
    // atom 5: C bridgehead (neighbors 4, 6, 10) -> locant "10a"
    // atom 6: C alpha (neighbors 5, 7) -> locant "5"
    // atom 7: C beta (neighbors 6, 8) -> locant "6"
    // atom 8: C beta (neighbors 7, 9) -> locant "7"
    // atom 9: C alpha (neighbors 8, 10) -> locant "8"
    // atom 10: C bridgehead (neighbors 5, 9, 11) -> locant "8a"
    // atom 11: C meso (neighbors 10, 12) -> locant "9"
    // atom 12: C bridgehead (neighbors 3, 11, 13) -> locant "9a"
    // atom 13: C alpha (neighbors 0, 12) -> locant "1"
    registry.push_back(RetainedRingSystem{
        "anthracene",
        "c1ccc2cc3ccccc3cc2c1",
        {"2","3","4","4a","10","10a","5","6","7","8","8a","9","9a","1"},
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","4"}},
            {'d', {"4","4a"}},  {'e', {"4a","10"}}, {'f', {"10","10a"}},
            {'g', {"10a","5"}}, {'h', {"5","6"}},   {'i', {"6","7"}},
            {'j', {"7","8"}},   {'k', {"8","8a"}},  {'l', {"8a","9"}},
            {'m', {"9","9a"}},  {'n', {"9a","1"}},
        },
        46,
        "anthra", // Blue Book FR-2.1.9 irregular stem: "anthra-" (not "anthraceneo-")
    });

    // Tetracene (also called naphthacene) -- 4 linearly-fused all-carbon
    // aromatic 6-membered rings, 18 carbons. Locant numbering independently
    // triangulated from three real external facts (no single source gave
    // the full locant list directly, so this was cross-checked rather than
    // taken from one place):
    //   1. Rubrene's real IUPAC name is "5,6,11,12-tetraphenyltetracene"
    //      (confirmed via a fresh web search of supplier/catalog listings
    //      giving this exact name) -- fixes the two meso-pair locants.
    //   2. "Tetracene-5,12-dione" is the standard/major tetracene quinone
    //      (ChEBI:51287, analogous to anthracene-9,10-dione) -- confirms
    //      5 and 12 are the two meso positions of the SAME inner ring
    //      (diagonally opposite within it), not of different rings.
    //   3. PubChem CID 12446083's IUPACName is
    //      "1,2,3,4,7,8,9,10-octahydrotetracene" (already independently
    //      confirmed this session for the ground-truth molecule) --
    //      confirms the two OUTER rings carry locants 1,2,3,4 and 7,8,9,10.
    //   4. IUPAC's own fused-ring-nomenclature Appendix 1 listing
    //      (iupac.qmul.ac.uk/fusedring/app1.html) places tetracene at
    //      table position 39, anthracene at 47, and naphthalene at 57 --
    //      consistently ONE MORE than this registry's own seniorityRank
    //      values for anthracene (46) and naphthalene (56), confirming
    //      the registry's seniorityRank = (that table position - 1), so
    //      tetracene's seniorityRank = 39 - 1 = 38 (more senior than
    //      anthracene/naphthalene, consistent with FR-2.3(b) "greater
    //      number of rings" -- 4 rings outranks 3 or 2).
    // These three facts, plus the standard ortho-fused-acene periphery-walk
    // pattern already used for naphthalene/anthracene above (each bridgehead
    // gets the previous number's own letter suffix), together fully pin
    // down the numbering below -- independently re-derived from first
    // principles AND cross-checked against all three facts simultaneously
    // (they are only mutually consistent with ONE locant assignment).
    //
    // SMILES 'c1ccc2cc3cc4ccccc4cc3cc2c1' parsed by Indigo (confirmed via a
    // fresh scratch Indigo run) produces atoms in order:
    // atom 0: C (neighbors 1, 17) -> locant "2"
    // atom 1: C (neighbors 0, 2) -> locant "3"
    // atom 2: C (neighbors 1, 3) -> locant "4"
    // atom 3: C bridgehead (neighbors 2, 4, 16) -> locant "4a"
    // atom 4: C meso, ring-B top (neighbors 3, 5) -> locant "5"
    // atom 5: C bridgehead (neighbors 4, 6, 14) -> locant "5a"
    // atom 6: C meso, ring-C top (neighbors 5, 7) -> locant "6"
    // atom 7: C bridgehead (neighbors 6, 8, 12) -> locant "6a"
    // atom 8: C (neighbors 7, 9) -> locant "7"
    // atom 9: C (neighbors 8, 10) -> locant "8"
    // atom 10: C (neighbors 9, 11) -> locant "9"
    // atom 11: C (neighbors 10, 12) -> locant "10"
    // atom 12: C bridgehead (neighbors 11, 7, 13) -> locant "10a"
    // atom 13: C meso, ring-C bottom (neighbors 12, 14) -> locant "11"
    // atom 14: C bridgehead (neighbors 13, 5, 15) -> locant "11a"
    // atom 15: C meso, ring-B bottom (neighbors 14, 16) -> locant "12"
    // atom 16: C bridgehead (neighbors 15, 3, 17) -> locant "12a"
    // atom 17: C (neighbors 16, 0) -> locant "1"
    registry.push_back(RetainedRingSystem{
        "tetracene",
        "c1ccc2cc3cc4ccccc4cc3cc2c1",
        {"2","3","4","4a","5","5a","6","6a","7","8","9","10","10a","11","11a","12","12a","1"},
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","4"}},
            {'d', {"4","4a"}},  {'e', {"4a","5"}},  {'f', {"5","5a"}},
            {'g', {"5a","6"}},  {'h', {"6","6a"}},  {'i', {"6a","7"}},
            {'j', {"7","8"}},   {'k', {"8","9"}},   {'l', {"9","10"}},
            {'m', {"10","10a"}},{'n', {"10a","11"}},{'o', {"11","11a"}},
            {'p', {"11a","12"}},{'q', {"12","12a"}},{'r', {"12a","1"}},
        },
        38,
        // fusionPrefix: verified directly against Blue Book FR-2.1.2
        // (fetched iupac.qmul.ac.uk/fusedring/FR21.html this session):
        // "In an attached component prefix the ending '-acene' is changed
        // to '-aceno'." Tetracene is NOT one of FR-2.1.9's four irregular
        // stems (naphtho/anthra/benzo/phenanthro only cover
        // naphthalene/anthracene/benzene/phenanthrene), so the regular
        // rule applies: tetracene -> tetraceno. (A prior guess of
        // "naphthaceno" here, by loose analogy with naphtho/anthra, was
        // never independently verified and is superseded by this check.)
        "tetraceno",
    });

    // Phenanthrene -- 14 carbons; angular 6-6-6 ring system.
    // seniorityRank=47 (lower number = higher seniority, per
    // FusionConstruction.h): one step less senior than anthracene (46),
    // far more senior than the generic cyclo-prefixed rings (>=100). The
    // specific Blue Book Appendix 1 table-position this was originally
    // cited against was not independently verified -- do not trust an
    // exact table-position claim here. The rank's practical correctness is
    // established empirically (dibenzo[c,g]phenanthrene parent-selection
    // fixture, tests/fusion_construction_test.cpp), not by that citation.
    registry.push_back(RetainedRingSystem{
        "phenanthrene",
        "c1ccc2ccc3ccccc3c2c1",
        {"3","2","1","10a","10","9","8a","8","7","6","5","4b","4a","4"},
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","4"}},
            {'d', {"4","4a"}},  {'e', {"4a","4b"}}, {'f', {"4b","5"}},
            {'g', {"5","6"}},   {'h', {"6","7"}},   {'i', {"7","8"}},
            {'j', {"8","8a"}},  {'k', {"8a","9"}},  {'l', {"9","10"}},
            {'m', {"10","10a"}},{'n', {"10a","1"}},
        },
        47,
        "phenanthro",
    });

    // Indene -- 9 carbons; 6-5 ortho-fused ring system.
    // seniorityRank=57 (lower number = higher seniority): one step less
    // senior than phenanthrene (47), still far more senior than the
    // generic cyclo-prefixed rings (>=100). Blue Book Appendix 1
    // table-position not independently verified -- see phenanthrene's
    // comment above for the same caveat; correctness here rests on the
    // passing fixture tests, not on a verified table citation.
    // SMILES 'C1C=Cc2ccccc12' producing atoms in order:
    // atom 0: C1 -> "1", atom 1: C2 -> "2", atom 2: C3 -> "3", atom 3: C3a -> "3a",
    // atom 4: C4 -> "4", atom 5: C5 -> "5", atom 6: C6 -> "6", atom 7: C7 -> "7",
    // atom 8: C7a -> "7a"
    registry.push_back(RetainedRingSystem{
        "indene",
        "C1C=Cc2ccccc12",
        {"1","2","3","3a","4","5","6","7","7a"},
        {
            {'a', {"1","2"}}, {'b', {"2","3"}}, {'c', {"3","3a"}},
            {'d', {"3a","4"}}, {'e', {"4","5"}}, {'f', {"5","6"}},
            {'g', {"6","7"}}, {'h', {"7","7a"}}, {'i', {"7a","1"}},
        },
        57,
        "indeno",
    });

    // Picene -- 22 carbons; 5-ring angular PAH (phenacene series).
    // seniorityRank=37 (Appendix 1 entry 37).
    // Fusion prefix "piceno" (regular -ene -> -eno stem per FR-2.1.2).
    //
    // Locant-to-atom correspondence independently re-derived and verified via
    // OPSIN (https://www.ebi.ac.uk/opsin/ws/<name>.json), NOT assumed from the
    // SMILES atom-walk order. Fetched real SMILES for 1,2-dihydropicene,
    // 3,4-dihydropicene, 5,6-dihydropicene, 13,14-dihydropicene,
    // 9,10-dihydropicene and 11,12-dihydropicene, then matched each against
    // this reference SMILES via an any-bond carbon-skeleton SMARTS (confirmed
    // atom-order-identical to this SMILES's own parse order) to find which
    // atom pair each name's sp3 CH2-CH2 saturation lands on. A prior version
    // of this entry assumed locants "1,2,3,4,..." simply followed the SMILES
    // atom order directly (atom 0 = "2", ..., atom 21 = "1") without this
    // check; that assumption was wrong -- it put the real "1,2" and "3,4"
    // edges, and the real "5,6" and "13,14" edges, on each other's atoms
    // (OPSIN's 1,2-dihydropicene lands on what that version called locants
    // "3","4", not "1","2"), which produced the wrong fusion bond letter for
    // dinaphtho[1,2-c:2',1'-m]picene (engine said bond 'a', should be 'c').
    // The real numbering walks this SMILES's atom indices in the OPPOSITE
    // direction, starting at atom 2: atom2=1, atom1=2, atom0=3, atom21=4,
    // atom20=4a, atom19=5, atom18=6, atom17=6a, atom16=6b, atom15=7,
    // atom14=8, atom13=8a, atom12=9, atom11=10, atom10=11, atom9=12,
    // atom8=12a, atom7=12b, atom6=13, atom5=14, atom4=14a, atom3=14b
    // (back to atom2=1). All six OPSIN-derived locant pairs above are
    // consistent with this single walk (verified, not just curve-fit to a
    // subset), including the fixture's own bond 'c'={"3","4"} and bond
    // 'm'={"9","10"}.
    registry.push_back(RetainedRingSystem{
        "picene",
        "c1ccc2c3ccc4c5ccccc5ccc4c3ccc2c1",
        {
            "3", "2", "1", "14b", "14a", "14", "13", "12b", "12a", "12",
            "11", "10", "9", "8a", "8", "7", "6b", "6a", "6", "5",
            "4a", "4"
        },
        {
            {'a', {"1", "2"}},
            {'b', {"2", "3"}},
            {'c', {"3", "4"}},
            {'d', {"4", "4a"}},
            {'e', {"4a", "5"}},
            {'f', {"5", "6"}},
            {'g', {"6", "6a"}},
            {'h', {"6a", "6b"}},
            {'i', {"6b", "7"}},
            {'j', {"7", "8"}},
            {'k', {"8", "8a"}},
            {'l', {"8a", "9"}},
            {'m', {"9", "10"}},
            {'n', {"10", "11"}},
            {'o', {"11", "12"}},
            {'p', {"12", "12a"}},
            {'q', {"12a", "12b"}},
            {'r', {"12b", "13"}},
            {'s', {"13", "14"}},
            {'t', {"14", "14a"}},
            {'u', {"14a", "14b"}},
            {'v', {"14b", "1"}},
        },
        37,
        "piceno",
    });

    registry.push_back(RetainedRingSystem{
        "fluorene",
        "c1cccc2c3ccccc3Cc12",
        {
            "1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8",
            "8a", "9", "9a"
        },
        {
            {'a', {"1", "2"}},
            {'b', {"2", "3"}},
            {'c', {"3", "4"}},
            {'d', {"4", "4a"}},
            {'e', {"4a", "4b"}},
            {'f', {"4b", "5"}},
            {'g', {"5", "6"}},
            {'h', {"6", "7"}},
            {'i', {"7", "8"}},
            {'j', {"8", "8a"}},
            {'k', {"8a", "9"}},
            {'l', {"9", "9a"}},
            {'m', {"9a", "1"}},
        },
        10,
        "",
        "9",
    });

    // Pentacene -- 22 carbons; 5 linearly-fused aromatic 6-rings.
    registry.push_back(RetainedRingSystem{
        "pentacene",
        "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1",
        {
            "2","3","4","4a","5","5a","6","6a","7","7a","8",
            "9","10","11","11a","12","12a","13","13a","14","14a","1"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","4"}},
            {'d', {"4","4a"}},  {'e', {"4a","5"}},  {'f', {"5","5a"}},
            {'g', {"5a","6"}},  {'h', {"6","6a"}},  {'i', {"6a","7"}},
            {'j', {"7","7a"}},  {'k', {"7a","8"}},  {'l', {"8","9"}},
            {'m', {"9","10"}},  {'n', {"10","11"}}, {'o', {"11","11a"}},
            {'p', {"11a","12"}},{'q', {"12","12a"}},{'r', {"12a","13"}},
            {'s', {"13","13a"}},{'t', {"13a","14"}},{'u', {"14","14a"}},
            {'v', {"14a","1"}},
        },
        36,
        "pentaceno",
    });

    // Chrysene -- 18 carbons; 4 angular zigzag 6-rings.
    registry.push_back(RetainedRingSystem{
        "chrysene",
        "c1cccc2c3ccc4ccccc4c3ccc21",
        {
            "1","2","3","4","4a","4b","5","6","6a","7","8","9","10","10a","10b","11","12","12a"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","4"}},
            {'d', {"4","4a"}},  {'e', {"4a","4b"}}, {'f', {"4b","5"}},
            {'g', {"5","6"}},   {'h', {"6","6a"}},  {'i', {"6a","7"}},
            {'j', {"7","8"}},   {'k', {"8","9"}},   {'l', {"9","10"}},
            {'m', {"10","10a"}},{'n', {"10a","10b"}},{'o', {"10b","11"}},
            {'p', {"11","12"}}, {'q', {"12","12a"}},{'r', {"12a","1"}},
        },
        42,
        "chryseno",
    });

    // Fluoranthene -- 16 carbons; 4 rings (naphthalene + benzene + 5-ring), 1 interior atom 10c.
    registry.push_back(RetainedRingSystem{
        "fluoranthene",
        "c1ccc2cccc3c4ccccc4-c1c23",
        {
            "1","2","3","3a","4","5","6","6a","6b","7","8","9","10","10a","10b","10c"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","3a"}},
            {'d', {"3a","4"}},  {'e', {"4","5"}},   {'f', {"5","6"}},
            {'g', {"6","6a"}},  {'h', {"6a","6b"}}, {'i', {"6b","7"}},
            {'j', {"7","8"}},   {'k', {"8","9"}},   {'l', {"9","10"}},
            {'m', {"10","10a"}},{'n', {"10a","10b"}},{'o', {"10b","1"}},
        },
        44,
        "fluorantheno",
        "",
        {"10c"},
    });

    // Perylene -- 20 carbons; 5 rings (two naphthalenes joined peri), 2 interior atoms 12c, 12d.
    registry.push_back(RetainedRingSystem{
        "perylene",
        "c1ccc2cccc3c4cccc5cccc(c45)c1c23",
        {
            "1","2","3","3a","4","5","6","6a","6b","7","8","9","9a",
            "10","11","12","12a","12d","12b","12c"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","3a"}},
            {'d', {"3a","4"}},  {'e', {"4","5"}},   {'f', {"5","6"}},
            {'g', {"6","6a"}},  {'h', {"6a","6b"}}, {'i', {"6b","7"}},
            {'j', {"7","8"}},   {'k', {"8","9"}},   {'l', {"9","9a"}},
            {'m', {"9a","10"}}, {'n', {"10","11"}}, {'o', {"11","12"}},
            {'p', {"12","12a"}},{'q', {"12a","12b"}},{'r', {"12b","1"}},
        },
        39,
        "peryleno",
        "",
        {"12c", "12d"},
    });

    // Acenaphthylene -- 12 carbons; 3 peri-fused rings (6,6,5), 1 interior atom 8b.
    registry.push_back(RetainedRingSystem{
        "acenaphthylene",
        "C1=Cc2cccc3cccc1c23",
        {
            "1","2","2a","3","4","5","5a","6","7","8","8a","8b"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","2a"}},  {'c', {"2a","3"}},
            {'d', {"3","4"}},   {'e', {"4","5"}},   {'f', {"5","5a"}},
            {'g', {"5a","6"}},  {'h', {"6","7"}},   {'i', {"7","8"}},
            {'j', {"8","8a"}},  {'k', {"8a","1"}},
        },
        48,
        "acenaphthyleno",
        "",
        {"8b"},
    });

    // Pyrene -- 16 carbons; 4 peri-fused 6-rings, 2 interior atoms 10b, 10c.
    registry.push_back(RetainedRingSystem{
        "pyrene",
        "c1ccc2ccc3cccc4ccc1c2c34",
        {
            "1","2","3","3a","4","5","5a","6","7","8","8a","9","10","10a","10b","10c"
        },
        {
            {'a', {"1","2"}},   {'b', {"2","3"}},   {'c', {"3","3a"}},
            {'d', {"3a","4"}},  {'e', {"4","5"}},   {'f', {"5","5a"}},
            {'g', {"5a","6"}},  {'h', {"6","7"}},   {'i', {"7","8"}},
            {'j', {"8","8a"}},  {'k', {"8a","9"}},  {'l', {"9","10"}},
            {'m', {"10","10a"}},{'n', {"10a","1"}},
        },
        40,
        "pyreno",
        "",
        {"10b", "10c"},
    });
}