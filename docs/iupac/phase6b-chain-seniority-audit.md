# Phase 6B Audit: Chain-as-Parent Seniority (P-44.1.1) Gates

Audit of the `38` FDA drugs (`36` + `2` carboxylate anion wrappers) rejected with `"A chain-based principal group outranks the ring in this structure; chain-as-parent seniority (P-44.1.1) for this ring/class combination is not yet supported."` in the v100/v101 3,311-drug sweep (`docs/iupac/phase6-stereo-seniority-generalization-plan.md`, Task 3).

**Primary code paths audited (`src/app/IupacNamer.cpp`):**
- `isChainParentWithRingSubstituentSupported` (`L8900-8906`) — already includes 30 `GroupType` values (`ACID`, `AMIDE`, `THIOAMIDE`, `HYDRAZIDE`, `AMIDINE`, `NITRILE`, `ALDEHYDE`, `KETONE`, `ALCOHOL`, `THIOL`, `SELENOL`, `TELLUROL`, `HYDROPEROXIDE`, `AMINE`, `IMINE`, `THIAL`, `THIONE`, `SULFONIC_ACID`, `SULFINIC_ACID`, `SULFINYL_HALIDE`, `PHOSPHONIC_ACID`, `PHOSPHONIC_DIHALIDE`, `ARSONIC_ACID`, `ARSONIC_DIHALIDE`, `ESTER`, `ACYL_HALIDE`, `SULFONYL_HALIDE`, `SULFONAMIDE`, `SULFINAMIDE`, `STIBONIC_ACID`). Only `PHOSPHINIC_ACID`, `ARSINIC_ACID`, `STIBINIC_ACID`, `BORONIC_ACID`, `BORINIC_ACID`, `PHOSPHINE` are omitted (`PHOSPHINIC_ACID`/`ARSINIC_ACID`/`STIBINIC_ACID`/`BORINIC_ACID` already short-circuit before `runChainParentSeniorDispatch`; `BORONIC_ACID`/`PHOSPHINE` only fail downstream).
- `nameChainParentWithRingSubstituent` (`L9270-9592`) — internal structural gates (`return "";`).
- `runChainParentSeniorDispatch` (`L23732-23858`) — caller + `cSidePureLinear` fallback (`L23780-23855`) and error emission (`L23856`).
- `selectSeniorParent` (`L15321-15403`) and `classifyFunctionalGroups` (`L14699-15254`).
- Pre-Phase-2 single-ring demotion (`ringSubstituentInfos`, `L23170-23518`).

None of the `38` drugs fail because of `isChainParentWithRingSubstituentSupported` missing a major class — **all 38 fail inside internal structural gates** grouped below into 7 shared root causes.

---

## Gate Summary Table

| Gate | Shared Internal Failure Point | Count (Rows / Unique) | Representative Drugs |
| :--- | :--- | :---: | :--- |
| **G1** | Heteroatom/acyl/amide/sulfonamide/polyether linker between `chainSet` and `ringNodeSet` (`L9328` carbon-only BFS + `L9351-9461` direct/single-bridge attachment check) | 17 / 11 | `AMINOHIPPURIC ACID`, `CILASTATIN`, `PROGLUMIDE`, `IODOHIPPURATE`, `METHENAMINE HIPPURATE`, `IOCETAMIC ACID`, `LODOXAMIDE`, `TROFINETIDE`, `ERDOSTEINE`, `SULFADICRAMIDE`, `PIFLUFOLASTAT`, `NONOXYNOL 9` |
| **G2** | Ring-demotion / ester-to-ring bridge & dead `mainChainExoCount == 1` guard (`L23176-23213`, `L9385`, `L9538`) | 6 / 4 | `METIPRANOLOL`, `MOXISYLYTE`, `DIPIVEFRIN`, `VITAMIN A PALMITATE` (`ESTER` on O-aryl or alkenyl chain) |
| **G3** | Disconnected principal chains separated by ring (`L9285-9292`, `L9361` `attachCount > 1` on `rNode != ipsoRingNode`) | 4 / 2 | `ALBUTEROL`, `LEVOSALBUTAMOL` |
| **G4** | Demoted ring/branch stereo & `E`/`Z` double-bond failure inside `nameChainParentWithRingSubstituent` -> `nameRingAsSubstituent` (`L9469-9491`) | 2 / 1 | `CARBOPROST`, `CARBOPROST TROMETHAMINE` |
| **G5** | Second ring in multi-ring heterocycle misclassified as acyclic chain `AMINE` (`L15078-15105`, `L15324`, `L34263`) | 1 / 1 | `NICOTINE` |
| **G6** | `AMINE` on acyclic side chain + carbamate on ring + chiral amine C (`L23784-23855` / `L9351`) | 2 / 1 | `RIVASTIGMINE`, `RIVASTIGMINE TARTRATE` |
| **G7** | Rootless acyl/guanidine/oxime linkers & ring-attachment gaps (`L14938`, `L15007`, `L23213`) | 6 / 3 | `DIETHYLCARBAMAZINE`, `PROGUANIL`, `FLUVOXAMINE` |

*(Total = `17 + 6 + 4 + 2 + 1 + 2 + 6 = 38` rows.)*

---

## Detailed Root-Cause Breakdown

### G1 — Non-Carbon / Functional Linker Between Principal Chain and Ring (17 rows)
- **Code location:** `src/app/IupacNamer.cpp:9320-9461` (`nameChainParentWithRingSubstituent` Step 2–3).
- **Mechanism:**
  1. Step 2 (`L9328`) grows `chainSet` exclusively along carbon-carbon bonds (`g.nodes[nei].atomicNumber == 6`). Any amide `-C(=O)-NH-`, sulfonamide `-S(=O)2-NH-`, thioether `-S-`, or repeating ether `-(CH2CH2O)n-` between the principal chain and the ring halts traversal.
  2. Step 3 (`L9351-9461`) only recognizes:
     - Direct `ring - C(chainSet)` bond (`L9353-9366`),
     - Single divalent `O`/`S` bridge `ring - O/S - C(chainSet)` (`L9367-9400`), or
     - Anilide fallback ONLY when `winningType == GroupType::AMIDE` and the amide `N` is directly bonded to `ringNodeSet` (`L9403-9459`).
  3. Consequently, whenever an `ACID`/`ALCOHOL`/`AMIDE` chain connects to a ring via `-NH-C(=O)-ring`, `-N(Ac)-ring`, `-NH-S(=O)2-ring`, `-C(=O)-NH-ring`, or multiple `-O-` / `-S-` bridges, `ipsoRingNode < 0` and `nameChainParentWithRingSubstituent` returns `""`.
- **Sub-patterns inside G1:**
  - **G1a (`chain(ACID) - NH - C(=O) - ring` / reversed carboxamido substituent on chain):** `AMINOHIPPURIC ACID`, `AMINOHIPPURATE SODIUM`, `IODOHIPPURATE`, `METHENAMINE HIPPURATE`, `PROGLUMIDE`, `TROFINETIDE`, `PIFLUFOLASTAT`, `PIFLUFOLASTAT F18`, `CILASTATIN`, `CILASTATIN SODIUM`.
    - Analog probe `CC(=O)NCC(O)=O` already names as `2-ethanamidoethanoic acid`, showing `nameAcyclicChainParentWithSubstituents` can name acylamido substituents on an acid chain once the `ring-C(=O)-NH-` or `ring-NH-C(=O)-` substituent string is supplied as an `extraSubstituents` entry (or via `nameBranchGraph` when crossing heteroatom bridges from `chainSet`)!
    - Note secondary sub-gates in `CILASTATIN` (chain split by internal `-S-` between heptenoic acid and cysteine acid + `(Z)`/`(R)`/`(S)`), `TROFINETIDE` (`ring` = proline with quaternary stereocenter + acylamide), `PIFLUFOLASTAT` (`urea` splitting the two carboxyl chains + `[18F]`).
  - **G1b (`ring - NH - C(=O) - chain(ACID)` / `ring - N(Ac) - chain(ACID)`):** `IOCETAMIC ACID`, `ERDOSTEINE`, `LODOXAMIDE`, `LODOXAMIDE TROMETHAMINE`. Note `LODOXAMIDE` has two identical oxamic-acid chains separated by the benzene ring (`L9362` analog / multiplicative `P-52`).
  - **G1c (`chain(AMIDE) - S(=O)2 - ring`):** `SULFADICRAMIDE`. Amide N-substituent is benzenesulfonyl (`-S(=O)2-ring`), rejected by `L9412` (`ringNodeSet.count(nNei)` requires direct `N-ring` bond without `-S(=O)2-`).
  - **G1d (`chain(ALCOHOL) - (O - CH2CH2)_n - ring`):** `NONOXYNOL 9`. Single-ether `NONOXYNOL_SHORT` (`CCc1ccc(OCCO)cc1`) succeeds via `L9377`; longer polyether fails because `chainSet` stops at the distal `-O-`.

### G2 — Ring-as-Ester-Alkyl vs Ring Substituent Bridge & Dead `mainChainExoCount == 1` Guard (6 rows)
- **Code location:** `src/app/IupacNamer.cpp:23176-23213` (pre-Phase-2 `ringSubstituentInfos`) and `L9385`, `L9538` (`nameChainParentWithRingSubstituent` / `selectSeniorParent`).
- **Mechanism:**
  1. For O-aryl esters (`METIPRANOLOL`, `METIPRANOLOL HYDROCHLORIDE`, `MOXISYLYTE`, `DIPIVEFRIN`, `DIPIVEFRIN HYDROCHLORIDE`), `combinedWinner = GroupType::ESTER` outranks the chain `AMINE`/`ALCOHOL`. Because `ESTER` acyl carbon (`CH3-C(=O)-`) is outside `ringNodeSet`, `selectSeniorParent` picks `CHAIN_PARENT`.
  2. Inside `nameChainParentWithRingSubstituent`, `L9385` explicitly blocks `esterOxygen` from acting as a substituent bridge (correctly: in `AcO-ring`, the ring IS the ester's O-aryl group, not an alkoxy substituent on acetic acid!).
  3. Meanwhile, pre-Phase-2 ring demotion (`L23176-23202`) already has code selecting the primary ester/chain attachment when `esterOCandidates.size() == 1 && !chainCCandidates.empty()` (`METIPRANOLOL`, `MOXISYLYTE`) or `chainCCandidates.size() > 1` (`FLUVOXAMINE_AMINE`), **but `L23213` checks `(mainChainExoCount == 1 || twoRingsAreDisjoint)`**! Because every multi-substituent ring has `mainChainExoCount >= 2`, `L23213` renders `L23176-23202` unreachable (`MTP_ANISOLE_AC` passes only because `-OCH3` isn't counted in `mainChainExoCount`, while `MTP_TOLUENE_AC` fails with `mainChainExoCount == 2`).
  4. `DIPIVEFRIN` has **two** pivalate `ESTER` groups attached to one benzene ring (diester of catechol -> requires multiplicative diester O-arylene naming).
  5. `VITAMIN A PALMITATE` fails `ESTER` O-alkyl naming when `nameBranchGraph` hits the alkenyl-substituted unsaturated ring (`CYCLOHEXENYL_ETH_ACETATE` -> `"Unsupported ester alkyl group."` / `"A chain-based..."`).

### G3 — Disconnected Principal Chains Separated by Ring (`ALBUTEROL`, 4 rows)
- **Code location:** `src/app/IupacNamer.cpp:9285-9292` and `L9351-9366`.
- **Mechanism:**
  - In `ALBUTEROL` / `LEVOSALBUTAMOL` (`CC(NCC(c1cc(CO)c(O)cc1)O)(C)C`), `winningType = GroupType::ALCOHOL`.
  - Phenol `-OH` is on `ringNodeSet` (`ringCount = 1`); benzylic `-CH2OH` (exocyclic C at ring C3) and `-CH(OH)CH2NHtBu` (C1) are both non-ring carbons -> `selectSeniorParent` gets `ringCount = 1, chainCount = 2` -> `CHAIN_PARENT`.
  - Step 1 (`L9285-9292`) collects **both** `-CH2OH` and `-CH(OH)CH2NHtBu` into `principalCarbons`.
  - Step 2 seeds DFS from both components; Step 3 hits two separate `ringNodeSet` attachments (`C3` and `C1`) and returns `""` at `L9362` (`else if (rNode != ipsoRingNode) return "";`).
  - Probes prove both halves work cleanly in isolation (`CC(NCC(c1ccccc1)O)(C)C` -> `2-[(1,1-dimethylethyl)amino]-1-phenylethanol`; `CC(NCCc1cc(CO)c(O)cc1)(C)C` -> `{5-(2-[(1,1-dimethylethyl)amino]ethyl)-2-hydroxyphenyl}methanol`).
  - Under Blue Book P-44.4.1, when one non-symmetric chain component strictly outranks another (e.g., 2-carbon ethanol `-CH(OH)CH2NHtBu` vs 1-carbon hydroxymethyl `-CH2OH`), selecting the senior single connected chain component leaves `-CH2OH` as a ring substituent `(hydroxymethyl)` inside `nameRingAsSubstituent`! Symmetric ties (such as `baseline.csv:1064` `COC(=O)CCc1ccc2cc(CCC(=O)OC)ccc2c1` and `LODOXAMIDE`) remain distinct (multiplicative parent).

### G4 — Demoted Ring/Alkenyl Branch Stereo / E-Z Drop Inside `nameChainParentWithRingSubstituent` (`CARBOPROST`, 2 rows)
- **Code location:** `src/app/IupacNamer.cpp:9469-9491` and `nameRingAsSubstituent` alkenyl/E-Z branch stereo handling.
- **Mechanism:**
  - Probes show `CCCCC(C)(O)C=CC1CCCC1CCCCCCC(O)=O` (flat, no E/Z, no ring OH) succeeds (`7-[2-(3-hydroxy-3-methylhept-1-en-1-yl)cyclopentyl]heptanoic acid`), `OC1CC(O)CC1CCCCCCC(O)=O` succeeds (`7-(2,4-dihydroxycyclopentyl)heptanoic acid`), and `CCCCC[C@@](O)(C)CCC1CCCC1CCCCCCC(O)=O` succeeds!
  - Adding `/C=C/` E/Z geometry (`CARBOPROST_RING_2OH_CHAIN`) or `/C=C/[C@H]1` (`CARBOPROST_STEREO_SIMPLE`) makes `nameRingAsSubstituent` (called from `L9469` with `branchBondStereo = nullptr` / stereocenter checks) return `""`, causing `L9491` to fail cleanly back to `"A chain-based principal group outranks..."`.

### G5 — Multi-Ring Heterocycle Substituent Classified as Acyclic `AMINE` (`NICOTINE`, 1 row)
- **Code location:** `src/app/IupacNamer.cpp:15078-15105`, `L15321-15403`, and `L34263-34290`.
- **Mechanism:**
  - For `CN1[C@H](c2cnccc2)CCC1` (`NICOTINE`), `ringCount == 2` chooses `pyridine` (`ringNodeSet`) as parent ring and names `1-methylpyrrolidin-2-yl` as `phase2ForcedRingSubstituentNames`, then calls `namePhase2Monocyclic()`.
  - Inside `namePhase2Monocyclic()`, `classifyFunctionalGroups` (`L15081`) only skips ring-internal `C-N` bonds if `ringNodeSet.count(i) > 0` (`pyridine` only). Carbons of the `pyrrolidine` ring are thus misclassified as acyclic `GroupType::AMINE`, and `selectSeniorParent` (`L15388`) checks `ringHasHeteroatom` only when `choice.ringCount > 0` (false for pyridine), erroneously selecting `CHAIN_PARENT`.
  - Analog `CN1CCCC1c1ccccc1` succeeds (`1-methyl-2-phenylpyrrolidine`) only because pyrrolidine outranks benzene (`ring1Eligible && !ring2Eligible`) so `ringNodeSet = pyrrolidine`. (`NICOTINE` is also tracked as `pinned=false` in `baseline.csv:740`.)

### G6 — Acyclic Chiral Amine + Masked Ring Carbamate (`RIVASTIGMINE`, 2 rows)
- **Code location:** `src/app/IupacNamer.cpp:23780-23856` and `L9285-9366`.
- **Mechanism:**
  - `CCN(C(Oc1cc([C@@H](N(C)C)C)ccc1)=O)C` combines a rootless phenyl carbamate (skipped via `SubordinateSkip` at `L14664` because an amine exists outside) and a tertiary benzylic `AMINE` (`[C@@H](N(C)C)C`).
  - `selectSeniorParent` picks `CHAIN_PARENT` (`AMINE`). Notice `principalCarbons` includes **all 3** N-bonded carbons (`CH(CH3)ring`, `CH3`, `CH3`), which breaks both `nameChainParentWithRingSubstituent` (`attachCount` / multi-seed) and causes the flat variant (`RIVASTIGMINE_FLAT`) to loop/duplicate N-substituents in `cSidePureLinear`; with `@`, stereo resolution rejects the fallback and hits `L23856`.

### G7 — Rootless Hetero-Acyl Linkers (`DIETHYLCARBAMAZINE`, `PROGUANIL`) & Substituted Oxime + `mainChainExoCount` (`FLUVOXAMINE`) (6 rows)
- **Code location:** `src/app/IupacNamer.cpp:14938`, `L15007`, `L23213`, and oxime naming.
- **Mechanism:**
  - **`DIETHYLCARBAMAZINE`** (`CCN(C(N1CCN(C)CC1)=O)CC`): Rootless urea `N-C(=O)-N(Et)2` attached to piperazine N1 is skipped at `L14938` (`!allSSSRNodesScan.count(c) && isOtherFuncWithRings`). Amine carbons on `N(Et)2` (`CH2CH3`) trigger `CHAIN_PARENT` (`ethanamine`), whose `chainSet` stops at the urea `N` (`attachCount == 0`) instead of recognizing `4-methylpiperazine-1-carboxamide` (Blue Book PIN `N,N-diethyl-4-methylpiperazine-1-carboxamide`).
  - **`PROGUANIL`** (`CC(NC(NC(Nc1ccc(Cl)cc1)=N)=N)C`): Rootless biguanide carbons are skipped at `L15008`, leaving the two external `AMINE` carbons (`iPr-NH` and `4-chlorophenyl-NH`). `selectSeniorParent` picks `CHAIN_PARENT` (`propan-2-amine`), whose `chainSet` stops at the biguanide `N` (`attachCount == 0`, while simpler `NC(=N)NC(=N)Nc1ccccc1` passes on the ring!).
  - **`FLUVOXAMINE`** (`COCCCC/C(/c1ccc(C(F)(F)F)cc1)=N\OCCN`): Fails two layered gates: (1) `L23213` dead `mainChainExoCount == 1` guard blocks ring demotion when `-CF3` coexists on benzene (`CCCC(c1ccccc1)OCCN` passes while `CCCC(c1ccc(C(F)(F)F)cc1)OCCN` fails), and (2) `-C(=NOCCN)-` oxime O-alkylamino bridge is rejected downstream once ring demotion is bypassed (`CC(=NOCCN)c1ccccc1` -> `"Substituted amine/hydrazine substituents are not supported in this phase."`).

---

## Generalization Priorities for Task 4 (Phase 6B)

Fixing the **shared structural mechanisms** resolves whole drug classes cleanly without per-drug special cases:
1. **Multi-ring / secondary-ring N-C exclusion in `classifyFunctionalGroups` / `selectSeniorParent` (`G5`):** Do not classify carbons whose C-N bond lies inside another `allSSSRRings` ring (or when a multi-ring forced substituent ring already owns that fragment) as competing chain `AMINE`s. Fixes `NICOTINE` + secondary hetero-ring cases (verify `baseline.csv:740` pin state: `pinned=false`).
2. **Dead `mainChainExoCount == 1` guard + `ringSubstituentInfos` (`G2` / `G7`):** Replace `mainChainExoCount == 1` at `L23213` with checking whether a unique principal attachment (`attachRingNode != -1 && foundAttachChainNode != -1`) was already disambiguated by `L23176-23202`! Immediately unlocks `METIPRANOLOL`, `MOXISYLYTE`, and `4-CF3` analogs.
3. **Senior single connected `chainSet` selection when `principalCarbons` spans multiple disconnected chain components (`G3` / `G1b`):** When `principalCarbons` contains multiple components separated by `ringNodeSet` (e.g., `ALBUTEROL`'s `-CH(OH)CH2NHtBu` vs `-CH2OH`, NOT symmetric identical chains like `baseline.csv:1064` or `LODOXAMIDE`), pick the strictly senior chain component per P-44.4.1 (count of principal groups -> longest carbon chain / highest substituent count -> strict tie-break), leaving the subordinate component to be named cleanly as a ring substituent by `nameRingAsSubstituent`!
4. **Heteroatom-linker attachment bridge in `nameChainParentWithRingSubstituent` (`G1a`, `G1b`, `G1c`, `G1d`):** When `attachCount == 0` and no single `O`/`S` bridge exists, search for an acyclic chain component attached to the ring through a non-carbon functional linker (amide `-NH-C(=O)-ring`, anilide/acylamino `-C(=O)-NH-ring` or `-N(R)-ring`, sulfonamide `-NH-S(=O)2-ring`, polyether `-(O-C-C)n-ring`) and inject `nameBranchGraph(g, linkerRootFromChain, chainAttachNode, ...)` at `chainAttachNode` as `extraSubstituents` (or handle `-SO2-ring` in the `AMIDE` N-substituent fallback `L9411` + sulfonamide/urea/carbamate). Unlocks `AMINOHIPPURIC ACID`, `IODOHIPPURATE`, `METHENAMINE HIPPURATE`, `PROGLUMIDE`, `IOCETAMIC ACID`, `SULFADICRAMIDE`, `NONOXYNOL 9`, `PIFLUFOLASTAT` (mono-amide & urea chain), `TROFINETIDE`, `ERDOSTEINE`, `CILASTATIN`.
5. **Pass `branchBondStereo` + ring stereocenters cleanly in `nameChainParentWithRingSubstituent` (`G4` + Review Focus #4/#5):** Thread `branchBondStereo` / `handledBranchBonds` through `nameRingAsSubstituent` and `nameAcyclicChainParentWithSubstituents` at `L9469` & `L9584` so `CARBOPROST` + ring-stereo + branch-stereo P-44.1.1 structures emit full `(R)/(S)/(E)/(Z)` prefixes cleanly.
6. **`RIVASTIGMINE` / `DIETHYLCARBAMAZINE` / `PROGUANIL` / `FLUVOXAMINE` (`G6`, `G7`):** Handle single-parent amide/amine N-linked urea/carbamate/oxime/biguanide linkers (`nameRingAsSubstituent` carbamoyloxy / carboxamide / oxime-ether / biguanide prefixes, and single `principalCarbons` seed for tertiary amines in `nameChainParentWithRingSubstituent`) where Blue Book P-44 makes the ring (e.g. piperazine-1-carboxamide `DIETHYLCARBAMAZINE`) or chain (`PROGUANIL`, `FLUVOXAMINE`, `RIVASTIGMINE`) senior.
