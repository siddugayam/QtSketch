# Phase 6C Ring-Coverage Scoping Report (`v102` Baseline)

**Scope:** Investigation-only taxonomy of the two largest ring-coverage rejection buckets in the 3,311-drug desalted FDA sweep after Phase 6B (`namer_results_desalted_v102.tsv`: `1,979 / 3,311` named = `59.77%`, `0` crashes, `0` hangs):
1. `"Fused, bridged, spiro, or multiple ring systems are not supported in Phase 2."`: **611 rows** (`448` unique SMILES)
2. `"Unrecognized or unsupported substituent on ring."`: **82 rows** (`63` unique SMILES)

Combined, these two buckets account for **693 / 1,332 rejections (52.0% of all remaining failures)**.

> **Related note:** The remaining P-44.1.1 `"A chain-based principal group outranks the ring..."` bucket dropped from `38` rows in `v101` to **7 rows** in `v102` (such as downstream alkenyl-stereo prostaglandin branches like `CARBOPROST`). Note also that `44` carboxylate-salt rows fail upstream with `"Cannot name the neutral acid for anion: Fused, bridged, spiro, or multiple ring systems are not supported in Phase 2."` and `11` with `"Cannot name the neutral acid for anion: Unrecognized or unsupported substituent on ring."`, bringing total impacted rows including carboxylate anions to `655` and `93`.

---

## 1. Bucket 1: `"Unrecognized or unsupported substituent on ring."` (82 rows, 63 unique SMILES)

### Structural Breakdown (Full Bucket Census)

| Sub-bucket | Rows | Unique SMILES | Primary `IupacNamer.cpp` Gate(s) | Representative Drugs |
|---|---:|---:|---|---|
| **U1. Monocyclic P-heterocycle / phosphoramide** | 2 | 2 | `L23740` / `L24839` (`Section 5` / monocyclic `nameSubstituentSubtree`) | `IFOSFAMIDE`, `MIFAMURTIDE ACID` |
| **U2. Monocyclic with complex functionalized alkyl/amide substituent** | 3 | 2 | `L23740` / `L24839` (tertiary amine / N-acyl on saturated ring / complex chain) | `CAPTOPRIL`, `TAPENTADOL` (+ `HCl`) |
| **U3. Multi-ring beta-lactams (cephalosporins / carbapenems)** | 28 | ~22 | `L23736` (`deferredSection5Error`) / `L26774` surfaced via `L34427` / `L34510` / `L34522` | `CEFTIZOXIME`, `CEFTRIAXONE`, `CEFPROZIL`, `CEFDINIR`, `AMOXICILLIN`, `MOXALACTAM` |
| **U4. Macrocycles / cyclic peptides / taxanes** | 6 | 5 | `L23736` / `L26774` (`deferredSection5Error`) | `PLERIXAFOR`, `ANIDULAFUNGIN`, `MICAFUNGIN`, `CASPOFUNGIN`, `QUINUPRISTIN` |
| **U5. Fused / disjoint multi-ring with unnamable branch/linker** | 43 | ~32 | `L23736` / `L26774` (`deferredSection5Error`) / `L30665` (fused ring sub) | `PACLITAXEL`, `DOCETAXEL`, `LAPATINIB`, `SAQUINAVIR`, `INDINAVIR`, `TERFENADINE`, `NAPHAZOLINE`, `SERTRALINE`, `LUMEFANTRINE`, `DROPERIDOL` |

### Key Finding on the 82-Row Bucket
- **77 of 82 rows (93.9%) have 2 or more SSSR rings** (`rc` distribution across all 82 rows: `rc=1: 5`, `rc=2: 1`, `rc=3: 36`, `rc=4: 20`, `rc=5: 10`, `rc=6: 6`, `rc=7: 3`, `rc=9: 1`). They report `"Unrecognized or unsupported substituent on ring."` rather than `"Fused, bridged, spiro..."` because `deferredSection5Error` (`src/app/IupacNamer.cpp:23736`, `26774`) captures the first failure in `ringSubstituentInfos` / Section 5 completeness and overrides the fallback Phase 2 error at `src/app/IupacNamer.cpp:34427`, `34510`, `34522`.
- **Only 5 of 82 rows (6.1%) are genuinely monocyclic** (`IFOSFAMIDE`, `CAPTOPRIL`, `TAPENTADOL` + `HCl`, `MIFAMURTIDE ACID`): solving monocyclic substituent gaps alone recovers at most 3-5 drugs; unlocking the rest of this bucket requires the multi-ring and saturated-heterocycle capabilities overlapping Bucket 2.

---

## 2. Bucket 2: `"Fused, bridged, spiro, or multiple ring systems are not supported in Phase 2."` (611 rows, 448 unique SMILES)

### Topology Taxonomy Across All 611 Rows
*(Computed via RDKit SSSR atom-ring connectivity on all 611 rows; 0 parse failures. Note: all rows have >= 2 SSSR rings, distributed as `rc=2: 70`, `rc=3: 177`, `rc=4: 172`, `rc=5: 111`, `rc=6: 43`, `rc=7: 18`, `rc=8: 5`, `rc=9: 13`, `rc=10: 2`.)*

#### Group A: Disjoint Rings Only (All SSSR rings isolated, `124` rows = `20.3%`)
| Sub-bucket | Rows | Key Structural Sub-types | Representative Drugs |
|---|---:|---|---|
| **A1. Two isolated rings (`rc=2`, `comps=2`)** | **70** | Separated by heterocycle linker (`29`), two separated benzenes / multiplicative-symmetric linkers (`28`), two non-aromatic carbocycles (`8`), direct ring-ring bond where senior-chain/sub dispatch stops (`5`) | `DIETHYLSTILBESTROL`, `OLSALAZINE`, `PENTAMIDINE`, `HEXESTROL`, `HEXACHLOROPHENE`, `FLURBIPROFEN`, `FELBINAC`, `NITROFURANTOIN`, `CLONIDINE`, `APRACLONIDINE`, `STAVUDINE`, `MUPIROCIN`, `DAPSONE` |
| **A2. Three or more isolated rings (`rc>=3`, `comps=rc`)** | **54** | Chains/heteroatom linkers carrying 3-5 separated benzene/heterocycle rings (`33` with `4+` rings, `21` with `3` rings); no `ringCount>=3` disjoint assembly loop exists | `RIMONABANT`, `KETOCONAZOLE`, `LEVOKETOCONAZOLE`, `POSACONAZOLE`, `PHENTOLAMINE`, `RILPIVIRINE`, `SULFINPYRAZONE`, `LETROZOLE`, `TERCONAZOLE` |

#### Group B: Single Polycyclic Ring Component (`83` rows = `13.6%`)
| Sub-bucket | Rows | Key Structural Sub-types | Representative Drugs |
|---|---:|---|---|
| **B1. Single fused polycyclic component (`comps=1`, no spiro/bridged)** | **61** | `[5,6,6,6]` steroid-like carbocycles (`35`), fused 3+ ring heterocycles (`19`), tricyclic `[6,x,6]` carbocycles like dibenzocycloheptenes (`7`) | `PREDNISONE`, `ESTRADIOL`, `ETHINYL ESTRADIOL`, `FINASTERIDE`, `METHOXSALEN`, `NEDOCROMIL`, `ANAGRELIDE`, `NORTRIPTYLINE`, `AMITRIPTYLINE`, `CYCLOBENZAPRINE` |
| **B2. Single polycyclic with bridged fusion** | **13** | Morphinans, bridged heterocycles/anhydrides, tricyclic azabicycles | `LEVORPHANOL`, `DEXTROMETHORPHAN`, `PENTAZOCINE`, `CANTHARIDIN`, `VARENICLINE` |
| **B3. Single polycyclic with spiro junction** | **9** | Spirolactones, griseofulvin, endoperoxides | `SPIRONOLACTONE`, `DROSPIRENONE`, `GRISEOFULVIN`, `FLUORESCEIN`, `ARTENIMOL` |

#### Group C: Multi-Component Ring Systems (`404` rows = `66.1%`)
| Sub-bucket | Rows | Key Structural Sub-types | Representative Drugs |
|---|---:|---|---|
| **C1. Fused system + >=1 isolated/fused companion ring (no spiro/bridged/macro)** | **330** | Bicyclic heterocycle (quinoline, indole, purine, benzimidazole, beta-lactam, etc.) + companion ring(s) (`222`), 3-ring fused + companion (`44`), 4-ring fused + companion (`22`), non-`6,6` bicyclic carbocycle + companion (`19`), `6,6` naphthalene/tetralin + companion (`18`), 5-ring fused + companion (`5`) | `EPIRUBICIN`, `PENICILLIN G`, `TROVAFLOXACIN`, `MOXIFLOXACIN`, `IRINOTECAN`, `PODOFILOX`, `KETOROLAC`, `LANSOPRAZOLE`, `RISPERIDONE`, `MONTELUKAST`, `TELMISARTAN`, `ALPELISIB`, `TUCATINIB` |
| **C2. Multi-component with bridged or spiro junction** | **57** | Tropane/morphinan/azabicyclic/spiro system + companion aryl/hetero rings; complex vinca alkaloids | `BUTORPHANOL`, `GRANISETRON`, `SPIRAPRIL`, `VINBLASTINE`, `LURASIDONE`, `SOLIFENACIN`, `SPARSENTAN`, `IBREXAFUNGERP` |
| **C3. Multi-component with macrocycle (>=12 ring atoms)** | **17** | Macrolides, macrocyclic HCV protease inhibitors, porphyrins, ansamycins | `SIROLIMUS`, `TACROLIMUS`, `EVEROLIMUS`, `SIMEPREVIR`, `TEMOPORFIN`, `PACRITINIB`, `IXABEPILONE` |

---

## 3. Exact `IupacNamer.cpp` Rejection / Bottleneck Sites

- **`ringSubstituentInfos` / Section 5 (`src/app/IupacNamer.cpp:23736`, `26774`, `26778`):** Sets `deferredSection5Error = "Unrecognized or unsupported substituent on ring."` when `ringCount >= 2` and substituent naming on a single ring fails, masking downstream Phase 2 multi-ring topology messages (`L34427`, `L34510`, `L34522`).
- **`namePhase2Monocyclic` edge gate (`L24839`):** Emits `"Unrecognized or unsupported substituent on ring."` when a monocyclic parent's exocyclic substituent cannot be named by `nameSubstituentSubtree` (e.g., complex polyfunctional chains/amines/phosphoramides on `U1`/`U2`).
- **Two disjoint rings handler (`L34395`–`L34575`):** Handles only direct C-C / direct ring-ring attachment when neither ring needs chain-seniority dispatch (`A1` gap: separated 2-ring linkers, P-44.2/P-45 multiplicative assemblies like `DIETHYLSTILBESTROL`, `OLSALAZINE`, `PENTAMIDINE`, and `direct_bond_2_rings` like `FLURBIPROFEN` / `FELBINAC` where the chain acid outranks the benzene ring).
- **Fused ring branching / non-6-membered fused orientation gate (`L13515`):** `detectFusedRingDirections` only handles 6-membered benzenoid/polycyclic ortho-fusion; anything with 5/7-membered or non-6,6 fusion hits `L13515` (`"Fused, bridged, spiro, or multiple ring systems are not supported in Phase 2."`).
- **Fused ring substituent fallback (`L30665`):** Emits `"Unrecognized or unsupported substituent on ring."` when `nameFusedRingSubstituent` and `nameBranchGraph` both return empty on a fused core.
- **Top-level `ringCount >= 3` / unhandled multi-component fallthrough (`L34732`):** Rejects any remaining 3+ ring or unsupported 2-ring fused/bridged/spiro/multi-component graph (`A2`, `B1`, `B2`, `B3`, `C1`, `C2`, `C3` = ~541 rows of Bucket 2).

---

## 4. Recommended Phase 6C Workstream Breakdown & Yield Estimates

Ordered by ROI (`FDA drugs recoverable / implementation risk`):

### Workstream 1 — Separated 2-Ring & Chain-Over-Multiple-Isolated-Rings Assemblies (`A1` + subset of `A2`: est. **~75–95 drugs**)
- **Target:** `A1` (`70` rows) + `A2` (`54` rows, especially common diphenylmethyl/dibenzyl/hetero-linked drugs like `RIMONABANT`, `CETIRIZINE`-like separated rings, `FLURBIPROFEN`, `FELBINAC`, `DIETHYLSTILBESTROL`, `HEXESTROL`, `PENTAMIDINE`, `OLSALAZINE`, `DAPSONE`).
- **Mechanism:**
  1. Generalize `nameChainParentWithRingSubstituent` beyond single-ring-component exclusion (`L9374` / `L9435`) so a senior acyclic chain can attach `>= 2` isolated rings (or a biphenyl `ring-ring` substituent via `nameBranchGraph`, recovering `FLURBIPROFEN` and `FELBINAC` immediately).
  2. Support general `nameBranchGraph` naming across heteroatom/chain bridges connecting 2 disjoint rings (fixing `L34539` fallthrough when rings are separated by a short chain/azo/amide/ester/ether bridge rather than directly bonded), plus symmetric P-15.3 multiplicative assembly for identical rings (`DIETHYLSTILBESTROL`, `OLSALAZINE`, `PENTAMIDINE`, `HEXACHLOROPHENE`, `PROBUCOL`, `DAPSONE`).
- **Risk:** Low–Medium (reuses existing monocyclic and acyclic namers without new ring topology algorithms).

### Workstream 2 — Retained & General Bicyclic Heterocycles as Substituents and Parents (`C1` 2-ring heterocycle subset: est. **~120–160 drugs**; highest volume)
- **Target:** The largest single sub-bucket in the entire sweep: **`C1` bicyclic heterocycle + companion ring(s) (`222` rows)** plus matching rows in `U5`/`A1` (indoles, quinolines, isoquinolines, purines, benzimidazoles, benzoxazoles, beta-lactam `[4,5]`/`[4,6]` penam/cephem cores when combined with `U3`).
- **Mechanism:**
  1. Ensure retained bicyclic heterocycle cores (and Hantzsch–Widman + benzo fused bicyclic cores via peripheral numbering) work cleanly both as substituents inside `nameRingSystemAsSubstituent` / `nameBranchGraph` AND as parent cores inside `ringCount >= 2` multi-component dispatch (`L34732`), instead of funneling only through `ringCount == 2` naphthalene/retained gates.
  2. Fix `U3`/`U5` (`L23736`, `L26774`) interaction so multi-ring parents don't fail prematurely in single-ring `ringSubstituentInfos` before fused-heterocycle / retained-core dispatch runs.
- **Risk:** Medium (requires clean peripheral locant numbering on common bicyclic heterocycles and unambiguous `-yl` attachment locants).

### Workstream 3 — Retained & Polycyclic Carbocycles + Tricyclics (`B1` `61` rows + `C1` `3+/non-66` carbocycle subset `~85` rows: est. **~65–90 drugs**)
- **Target:** Steroids with `5,6,6,6` tetracyclic skeletons (`35` in `B1`, e.g. `PREDNISONE`, `ESTRADIOL`, `ETHINYL ESTRADIOL`, `FINASTERIDE` whose ketone/amide/alcohol functional patterns differ from already-named testosterone esters), tricyclic `6,7,6` / `6,6,6` scaffolds (`NORTRIPTYLINE`, `AMITRIPTYLINE`, `CYCLOBENZAPRINE`), and `[6,6]` naphthalene/tetralin + companion ring drugs (`18` in `C1` + `SERTRALINE` in `U5`).
- **Mechanism:**
  1. Expand retained-ring registry / von Baeyer vs retained-core functional-group coverage (`3-oxo`, `17-one`, `17-ol`, `lactam` aza-steroids) and non-6-membered ortho-fusion support in `detectFusedRingDirections` (`L13515`).
  2. Route `naphthalene/tetralin + companion ring` (`18` + `SERTRALINE`) through `nameRingSystemAsSubstituent` as a substituent or parent paired with existing phenyl/saturated-ring substituent naming.
- **Risk:** Medium–High (steroid and tricyclic retained-name numbering/stereochemistry vs von Baeyer rules require careful golden-set PIN calibration).

### Workstream 4 — Bridged, Spiro, & Saturated-Heterocycle Extensions (`B2`/`B3`/`C2`: `79` rows + `U1`/`U2`/`U3`: est. **~30–50 drugs**; defer `C3`/`U4` macrocycles)
- **Target:** General von Baeyer bridged substituents (`quinuclidine`/`tropane`/`morphinan` scaffolds), oxa/aza spiro ring systems (`SPIRONOLACTONE`, `GRISEOFULVIN`, `IRBESARTAN`), and remaining monocyclic gaps (`TAPENTADOL` tertiary-amine alkyl branch, `CAPTOPRIL` N-acylpyrrolidine).
- **Explicit Deferral:** Macrocycles & cyclic peptides (`C3` = `17` + `U4` = `6`: `ANIDULAFUNGIN`, `CASPOFUNGIN`, `SIROLIMUS`, `PACLITAXEL` taxanes, `VINBLASTINE` bis-indole alkaloids) should be **explicitly deferred** beyond Phase 6C due to high combinatorial complexity and low ROI (`~23` rows total).

---

## 5. Summary Yield Projection

| Workstream | Target Pool (`v102`) | Conservative Yield | Cumulative Projected Coverage (`vs 1,979 = 59.77%`) |
|---|---:|---:|---:|
| WS1: Separated 2-ring & isolated multi-ring chains (`A1`+`A2`) | 124 rows | +75 to +95 | ~62.0% – 62.6% |
| WS2: Bicyclic heterocycles as substituents & multi-comp parents (`C1`+`U3`/`U5`) | ~250 rows | +120 to +160 | ~65.7% – 67.5% |
| WS3: Polycyclic/steroid/tricyclic carbocycles + tetralins (`B1`+`C1`) | ~145 rows | +65 to +90 | ~67.6% – 70.2% |
| WS4: Bridged/spiro + monocyclic sub cleanup (`B2`/`B3`/`C2`+`U1`/`U2`) | ~84 rows | +30 to +50 | **~68.5% – 71.7%** |
