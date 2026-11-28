# Phase 6: Branch Stereocenters & Chain-Seniority Generalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the two largest-ROI rejection buckets left after Tasks 1-11 (`v100`: `1835/3311` named, `1476` rejected) — substituent-branch stereocenters (`153` drugs) and the remaining P-44.1.1 chain-outranks-ring seniority gaps (`38` drugs) — by fixing the single shared mechanism behind each bucket, not by special-casing the named drugs.

**Architecture:** Both buckets reject through a single, already-general dispatch point reached from many call sites — not N independent bugs. Phase 6A threads the existing (but under-used) `formatBranchStereoPrefix` mechanism through every branch-naming call path that currently skips it, so a stereocenter inside a substituent branch gets consumed there instead of falling through to `formatStereoPrefix`'s top-level reject. Phase 6B investigates, then fixes, why `nameChainParentWithRingSubstituent` itself returns empty for 7 known drugs even though the one gating whitelist (`isChainParentWithRingSubstituentSupported`) already admits their functional-group classes. Phase 6C (multi-ring parents / unrecognized ring substituents, `722` drugs combined) is **not implemented in this plan** — its two sample drugs don't characterize 6+ scattered rejection sites well enough to scope tasks; Task 5 produces the investigation needed to write that plan next.

**Tech Stack:** C++17, Qt6 (QString/containers only, no UI), Indigo cheminformatics library, CTest, CMake/MinGW.

**Spec:** This plan's problem statements are grounded directly against `src/app/IupacNamer.cpp` on branch `dev` (post-commit `e646d024`) and the FDA 3,311-drug sweep (`.claude_scratch_fdafull/v100_sweep_output.txt`, `namer_results_desalted_v100.tsv`). No separate design doc exists for Phase 6; this plan **is** the spec, verified against code and data in the same session that wrote it.

## Global Constraints

- Build: PowerShell, from `build\`, `C:\mingw64\bin\mingw32-make.exe -j4 all` (full rebuild after every `IupacNamer.cpp` change — partial targets leave test binaries stale).
- Run any Qt6-linked test binary or `ctest` via **Bash**, not PowerShell, with `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin"` **appended** (never prepended) — PowerShell's base PATH shadows the project's MinGW runtime with conda's bundled one and makes every test binary exit `0xc0000139`.
- Zero regressions gate on every task: `iupac_namer_test.exe`, `iupac_golden_test.exe`, `retained_ring_registry_test.exe`, `retained_ring_matcher_test.exe`, `ctest -E imago_service_test` must all stay at current pass counts or better (currently `1642/0`, `809/0, 0 pin mismatches`, `91/0`, `13/0`, `27/27`).
- FDA sweep re-run at the end of each task (copy the newest `run_fda_sweep_v*.py` in `C:\Users\jp18b\.gemini\antigravity\brain\1c192a82-12a9-451a-a61a-93d5b7b711c8\scratch\`, bump the version number throughout, redirect stdout straight to a file — never pipe through `tail` on a long-running background run, it truncates the script's own summary lines out of the captured file).
- No drug name may appear in production code (`IupacNamer.cpp`) as a literal special case. Named drugs in this plan are **test fixtures only** — each one is a stand-in for a structural pattern class; if a task's fix only makes the named drug pass without generalizing to its pattern class, the task is not done.
- Never commit/push without the user's explicit instruction (standing rule this session).

## Review Focus

- A substituent branch with **two or more independent stereocenters** (not just one) — `formatBranchStereoPrefix`'s locant/descriptor ordering must stay correct when more than one center is consumed from the same branch call.
- A stereocenter that sits **exactly on the ring/branch boundary atom itself** (the attachment carbon), not strictly inside the branch — must not be double-counted by both `formatStereoPrefix` (parent-level) and `formatBranchStereoPrefix` (branch-level), and must not be silently dropped by neither.
- A chain-parent dispatch where the **ring substituent itself also carries a stereocenter** (e.g. a stereocenter on the ring that's demoted to substituent by Phase 6B's fix) — confirm Phase 6B's fix doesn't bypass Phase 6A's stereo wiring for the newly-demoted ring.
- A molecule that would hit **both** the 6A and 6B gates at once (chain-parent-wins class with a stereocenter in the demoted ring's substituent chain) — not represented by either task's named-drug list; add a combined fixture.
- Re-running `isChainParentWithRingSubstituentSupported`'s whitelist against the *current* `GroupType` enum before touching `nameChainParentWithRingSubstituent`'s body — a new `GroupType` value added since this plan was written would silently stay unsupported with no test catching it.

---

### Task 1: Audit stereo-wiring gaps across branch-naming call paths (6A investigation)

**Files:**
- Read only: `src/app/IupacNamer.cpp` (`nameRingSystemAsSubstituent` ~2805-2920, `nameBranchGraph` ~3220-3600 and its recursive calls, `nameSubstituentSubtree` ~2769-3220, `formatBranchStereoPrefix` definition at 9587, its 2 current call sites at 2518 and 3498, `formatStereoPrefix` definition at 9546)
- Create: `docs/iupac/phase6a-stereo-wiring-audit.md` (the task's deliverable — a table, not code)

**Interfaces:**
- Consumes: nothing (first task).
- Produces: a table of every branch-producing call path, each row = {call path name, line}, {currently calls `formatBranchStereoPrefix`? yes/no}, {if no: which structural pattern class it represents (ring-as-substituent / N-acyl branch / ester branch / plain alkyl branch / other)}. Task 2 consumes this table directly — it is the list of fix sites.

- [ ] **Step 1: Re-run the FDA v100 sweep's rejection list filtered to `"Stereocenters on substituent branches are not supported in this phase."` and extract every `(DRUG_NAME, SMILES)` pair**

Run (Bash): `grep -F "Stereocenters on substituent branches are not supported in this phase." .claude_scratch_fdafull/namer_results_desalted_v100.tsv | cut -f1,2,3`
Expected: ~153 rows (confirm the exact count against the live file; `153` is the report's figure, re-measure rather than trust it).

- [ ] **Step 2: For each of the 9 named fixtures (`DILTIAZEM`, `LOVASTATIN`, `SIMVASTATIN`, `COLCHICINE`, `LINEZOLID`, `CHLORAMPHENICOL`, `SELEGILINE`, `MEDROXYPROGESTERONE ACETATE`, `HYDROCORTISONE ACETATE`), run `smiles_driver_v72.exe <SMILES>` individually and capture the exact reject message and a GDB backtrace at the `formatStereoPrefix` reject line (9565) to find which caller walked that stereocenter's branch without going through `formatBranchStereoPrefix`**

Run (Bash, per drug): `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin" && gdb -batch -ex "break IupacNamer.cpp:9565" -ex run -ex bt --args build/smiles_driver_v72.exe "<SMILES>"`
Expected: a backtrace naming the specific branch-producing function (e.g. `nameBranchGraph` at a specific recursive call, or `nameRingSystemAsSubstituent`) that produced the branch name without registering its stereocenter in `handledBranchStereoIds`/`graphIdToLocant`.

- [ ] **Step 3: Group the 9 backtraces (plus any additional pattern classes the ~153-row list in Step 1 suggests by inspection) into distinct structural pattern classes, and write the audit table**

Write `docs/iupac/phase6a-stereo-wiring-audit.md` with one row per distinct call path found, per the Produces block above. Do not list the 9 drug names as the table's rows — list the call paths; attach drug names only as the "confirmed by" column.

- [ ] **Step 4: Commit the audit document**

```bash
git add docs/iupac/phase6a-stereo-wiring-audit.md
git commit -m "docs: audit branch-stereo wiring gaps for Phase 6A"
```

---

### Task 2: Wire `formatBranchStereoPrefix` through the gaps found in Task 1 (6A implementation)

**Files:**
- Modify: `src/app/IupacNamer.cpp` (exact call paths named in Task 1's audit table — not predetermined here, since Task 1 is what determines them)
- Test: `tests/iupac_namer_test.cpp`

**Interfaces:**
- Consumes: Task 1's audit table (`docs/iupac/phase6a-stereo-wiring-audit.md`) — the list of call paths to fix.
- Consumes: existing `formatBranchStereoPrefix(const Graph &g, int rootIdx, int parentIdx, const std::map<int, QChar> &stereoByGraphId, std::set<int> &handledBranchStereoIds, const std::vector<std::set<int>> &allRings = {}, const std::map<std::pair<int,int>, QChar> &branchBondStereo = {}, std::set<std::pair<int,int>> *handledBranchBonds = nullptr)` (line 9587) — do not change its signature; call it from each new site the same way the 2 existing call sites (2518, 3498) do.
- Produces: no new public symbols — this task only adds call sites to the existing function.

- [ ] **Step 1: Write 11 failing tests in `tests/iupac_namer_test.cpp`: the 9 named fixtures, plus 2 Review-Focus fixtures**

```cpp
TEST(verify "Task 6A branch stereocenter (DILTIAZEM)" produces a name containing
     the correct stereo descriptor instead of REJECTED, exact expected string
     determined by hand-deriving the IUPAC name from the SMILES before writing
     the assertion)
```
One such test per named fixture; use the project's existing `EXPECT_EQ`/golden-style pattern already used for other `[PASS]`-tagged cases in this file (see nearby Task 10/11 tests for the exact macro/harness style). Do not guess the expected name string — derive it by hand against Blue Book stereo-citation rules (P-92 et seq.) before writing the assertion, since a wrong hard-coded expectation would pass trivially once the branch stops rejecting but still be the wrong name.

Add 2 more fixtures for Review Focus items 1-2: (a) a substituent branch containing **two independent stereocenters** (pick or construct a SMILES from the Step 1-of-Task-1 153-row list that has 2 stereo-marked atoms in the same branch; if none exists in that list, construct a minimal synthetic SMILES), asserting both descriptors appear with correct relative locants; (b) a molecule where the stereocenter sits **on the ring/branch attachment atom itself**, asserting it is cited exactly once (not duplicated between `formatStereoPrefix` and `formatBranchStereoPrefix`, not dropped by both).

- [ ] **Step 2: Build and run; verify all 9 new tests fail with the current `"Stereocenters on substituent branches are not supported in this phase."` message**

Run: full rebuild per Global Constraints, then `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin" && ./build/iupac_namer_test.exe`
Expected: 9 new FAILs with that exact reject string in the failure output; all pre-existing tests still pass.

- [ ] **Step 3: For each call path in Task 1's audit table, add a `formatBranchStereoPrefix(...)` call matching the pattern at line 2518/3498 (compute the branch name first via the existing `nameBranchGraph`/`nameRingSystemAsSubstituent` call, then call `formatBranchStereoPrefix` with the same root/parent indices and splice its returned prefix onto the branch name exactly as the existing 2 sites do)**

No new code pattern to invent — copy the existing splice logic (lines 2517-2524 and 3495-3510) at each new site, adapting only the local variable names.

- [ ] **Step 4: Rebuild and run the full test suite + full FDA sweep (bump version number per Global Constraints)**

Run: per Global Constraints zero-regressions gate, plus `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin" && ./build/iupac_namer_test.exe` (expect the 9 new tests now PASS)
Expected: zero regressions on the 4 test binaries + ctest; FDA sweep gained count close to but not necessarily exactly 153 (some of the 153 may need Task 2's fix plus something else — note any remaining REJECTED rows with the same message for later follow-up, don't force every one to pass in this task).

- [ ] **Step 5: Commit**

```bash
git add src/app/IupacNamer.cpp tests/iupac_namer_test.cpp
git commit -m "feat: wire formatBranchStereoPrefix through remaining branch-naming call paths (Phase 6A)"
```

---

### Task 3: Trace why `nameChainParentWithRingSubstituent` returns empty for the 7 named drugs (6B investigation)

**Files:**
- Read only: `src/app/IupacNamer.cpp` (`nameChainParentWithRingSubstituent` definition starting 9222, its single gating whitelist `isChainParentWithRingSubstituentSupported` at 8860, the terminal reject at ~23770 and its surrounding lambda `runChainParentSeniorDispatch` starting ~23649)
- Create: `docs/iupac/phase6b-chain-seniority-audit.md`

**Interfaces:**
- Consumes: nothing (independent of Task 1/2).
- Produces: a table of {drug, SMILES, which internal check inside `nameChainParentWithRingSubstituent` returns `""`/fails first, and what structural condition that check is testing for}. Task 4 consumes this table as its fix list.

- [ ] **Step 1: Confirm `isChainParentWithRingSubstituentSupported` already admits each of the 7 drugs' winning `GroupType`**

Run (Bash, per drug, with a temporary debug print or GDB watch on `combinedWinner` at line 23688): confirm each of `AMINOHIPPURIC ACID` (AMIDE/ACID), `ALBUTEROL` (ALCOHOL/AMINE), `RIVASTIGMINE` (ESTER), `DIETHYLCARBAMAZINE` (AMIDE), `PROGUANIL` (AMIDINE-like), `FLUVOXAMINE` (AMINE), `CILASTATIN` (ACID) maps to a `GroupType` already in the whitelist at line 8860-8865.
Expected: all 7 pass the whitelist (per the plan author's read of the enum values already present) — if any does NOT pass, that is itself the fix (add it to the whitelist) and Task 4 is simpler than assumed for that drug; record either outcome in the table.

- [ ] **Step 2: For each drug whose `GroupType` passes the whitelist, GDB-trace into `nameChainParentWithRingSubstituent` (9222) to find the first `return ""` reached**

Run (Bash): `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin" && gdb -batch -ex "break IupacNamer.cpp:9222" -ex run -ex "next 200" --args build/smiles_driver_v72.exe "<SMILES>"` (adjust step count to reach the function's end; alternative: add a temporary `qDebug()` at each `return ""` site inside the function and rebuild once).
Expected: a specific line number and the structural condition it tests (e.g. "chain reachable set includes a second ring attachment" or "principal carbon count != 1") for each of the 7 drugs.

- [ ] **Step 3: Write the audit table to `docs/iupac/phase6b-chain-seniority-audit.md`**

One row per drug (or grouped by shared failing condition, if more than one drug shares the same internal `return ""` site — group them, since that is the shape of the real fix).

- [ ] **Step 4: Commit**

```bash
git add docs/iupac/phase6b-chain-seniority-audit.md
git commit -m "docs: audit chain-parent-with-ring-substituent dispatch gaps for Phase 6B"
```

---

### Task 4: Fix the internal gates found in Task 3 (6B implementation)

**Files:**
- Modify: `src/app/IupacNamer.cpp` (exact lines inside `nameChainParentWithRingSubstituent`, 9222-~9550, named by Task 3's audit table — not predetermined here)
- Test: `tests/iupac_namer_test.cpp`

**Interfaces:**
- Consumes: Task 3's audit table (`docs/iupac/phase6b-chain-seniority-audit.md`).
- Produces: no new public symbols — modifies the body of the existing function only; its signature (consumed at the 6 call sites: 23693, 27605, 27946, 28217, 33928, 33992) does not change.

- [ ] **Step 1: Write 9 failing tests in `tests/iupac_namer_test.cpp`: the 7 named fixtures, plus 2 Review-Focus fixtures**

Same pattern as Task 2 Step 1 — hand-derive each expected IUPAC name against Blue Book P-44.1.1 before writing the assertion.

Add 2 more fixtures for Review Focus items 3-4: (a) a chain-parent-wins molecule where the **demoted ring substituent itself carries a stereocenter**, asserting Task 2's `formatBranchStereoPrefix` wiring still fires on the newly-demoted ring (confirms Task 4 doesn't bypass Task 2's fix path); (b) a **combined 6A+6B fixture** — a chain-parent-wins molecule whose demoted ring's substituent chain (not the ring itself) carries the stereocenter, asserting both mechanisms compose correctly. Construct both synthetically if no FDA sweep row provides them.

- [ ] **Step 2: Build and run; verify all 7 fail with `"A chain-based principal group outranks the ring in this structure; chain-as-parent seniority (P-44.1.1) for this ring/class combination is not yet supported."`**

Run: full rebuild, then `export PATH="$PATH:/c/Qt/6.11.1/mingw_64/bin" && ./build/iupac_namer_test.exe`

- [ ] **Step 3: Fix each distinct internal gate found in Task 3's table**

For each grouped condition in the audit table, relax or correct the specific `return ""` guard so it admits the structural pattern the 7 drugs represent, without weakening it to accept the cases it's correctly rejecting elsewhere (check this against the full regression suite in Step 4, not by inspection alone).

- [ ] **Step 4: Rebuild, run full test suite + ctest, re-run FDA sweep (bump version)**

Expected: zero regressions; the 7 new tests pass; FDA sweep gained count reported honestly (don't force to exactly 38 — some of the 38 may share a condition not among the 7 named fixtures and need follow-up, same caveat as Task 2 Step 4).

- [ ] **Step 5: Commit**

```bash
git add src/app/IupacNamer.cpp tests/iupac_namer_test.cpp
git commit -m "fix: generalize nameChainParentWithRingSubstituent dispatch gates (Phase 6B)"
```

---

### Task 5: Scope Phase 6C as its own plan (investigation only — no code in this task)

**Files:**
- Read only: `src/app/IupacNamer.cpp` (the 6 `"Unrecognized or unsupported substituent on ring."` sites: 23429/23433, 24493, 25378, 26418/26422, 30297; the 2 `"Fused, bridged, spiro, or multiple ring systems are not supported in Phase 2."` sites: 13236, 34059)
- Create: `docs/iupac/phase6c-ring-coverage-scoping.md`

**Interfaces:**
- Consumes: nothing.
- Produces: a scoping document that becomes the **Spec** input for a future Phase 6C plan (not written in this plan — `722` combined drugs across 2 message types and 8 scattered code sites is too large and too uncharacterized to task-break here; one sample drug per message type does not establish how many distinct root causes exist behind either message).

- [ ] **Step 1: Pull every `(DRUG, SMILES)` row from `namer_results_desalted_v100.tsv` matching either message, and bucket them by rough structural shape (symmetric multi-aryl parent / fused polycyclic / spiro / bridged / other) by inspecting a sample of at least 40 SMILES per message**

Run (Bash): `grep -F "Unrecognized or unsupported substituent on ring." .claude_scratch_fdafull/namer_results_desalted_v100.tsv | cut -f1,2,3 > phase6c_unrecognized_sample.tsv` (similarly for the fused/bridged/spiro message).

- [ ] **Step 2: For each bucket, identify which of the 6 (or 2) code sites each sample row actually reaches (GDB or temporary debug print), to determine whether the 6 sites are 6 independent gates or 1 shared helper called from 6 places**

This mirrors Task 1/3's investigation pattern — do not assume groupings, confirm them.

- [ ] **Step 3: Write `docs/iupac/phase6c-ring-coverage-scoping.md` with the bucket breakdown, each bucket's real code site(s), and a rough drug-count-per-bucket estimate, suitable as the Spec for a dedicated Phase 6C plan written afterward**

- [ ] **Step 4: Commit**

```bash
git add docs/iupac/phase6c-ring-coverage-scoping.md
git commit -m "docs: scope Phase 6C ring-coverage gaps for a future plan"
```

---

## Execution Order Summary

1. Task 1 before Task 2 (audit before fix); Task 3 before Task 4 (audit before fix). Tasks {1,2} and {3,4} are independent of each other and may run in either order or in parallel.
2. Task 5 has no dependency on 1-4 and may run at any point; it produces no production-code change.
3. Phase 6C itself (the fix) is explicitly out of scope for this plan — Task 5's output is the next plan's spec, not a task list here.
