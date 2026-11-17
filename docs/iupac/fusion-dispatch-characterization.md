# Fusion-Dispatch Call Site Characterization

This document characterizes all three call sites for the fusion-naming functions
`nameLinearChainFusion` and `nameBranchingFusion` in the IUPAC namer codebase
(as of commit `0a2f2d5` on branch `dev`).

Both functions share the same signature:

```cpp
std::vector<RingSystemNameResult> nameLinearChainFusion(
    const Graph &g,
    int mol,
    const std::map<int, int> &indigoToGraphIdx,
    int ringCount,
    int freeValenceAtom = -1,
    const std::vector<std::set<int>> *ringSubsetOverride = nullptr,
    const std::vector<std::pair<int, QString>> &extraSubstituents = {});

std::vector<RingSystemNameResult> nameBranchingFusion(
    const Graph &g,
    int mol,
    const std::map<int, int> &indigoToGraphIdx,
    int ringCount,
    int freeValenceAtom = -1,
    const std::vector<std::set<int>> *ringSubsetOverride = nullptr,
    const std::vector<std::pair<int, QString>> &extraSubstituents = {});
```

Both functions start by creating a `FusionGraphContext` (introduced in an earlier task),
which wraps the ~51-line setup block that was previously duplicated at the top of
each function. The context computes SSSR rings (either from `ringSubsetOverride` or
via Indigo) and checks pairwise disjointness, providing `sssrRings` and `allDisjoint`
to the rest of the function logic.

---

## Call Site 1: `nameRingSystemAsSubstituent` (ring-as-substituent context)

**Location:** `IupacNamer.cpp:2870` and `2887`

**Function:** `nameBranchGraph` (line 3098, NOT its forward declaration at 1790) calls
`nameRingSystemAsSubstituent` (line 3122), which contains the actual fusion-dispatch
calls.

**Call expressions:**
```cpp
// Line 2870
std::vector<RingSystemNameResult> chainResults = 
    nameLinearChainFusion(g, mol, indigoToGraphIdx, ringCount, rootAtom);

// Line 2887
std::vector<RingSystemNameResult> branchResults = 
    nameBranchingFusion(g, mol, indigoToGraphIdx, ringCount, rootAtom);
```

### Site 1 Input Construction

- **`g`**: The `Graph` object passed as parameter to `nameRingSystemAsSubstituent`
- **`mol`**: The molecule handle passed as parameter to `nameRingSystemAsSubstituent`
- **`indigoToGraphIdx`**: Built locally from `g.nodes` (lines 2812-2815):
  ```cpp
  std::map<int, int> indigoToGraphIdx;
  for (const auto &node : g.nodes) {
      indigoToGraphIdx[node.indigoIdx] = node.id;
  }
  ```
- **`ringCount`**: Computed via BFS traversal of ortho-fused rings (lines 2827-2849).
  Starts from the ring directly containing `rootAtom` (`ringSystemNodes`), then walks
  adjacency (two rings sharing exactly 2 atoms) to find the full transitively-connected
  fused system. Result is the size of `fusedSystemRingIdx` (line 2850).
- **`freeValenceAtom` (5th param)**: `rootAtom` (the attachment point of this ring
  system as a substituent)
- **`ringSubsetOverride` (6th param)**: NOT passed (defaults to `nullptr`)
- **`extraSubstituents` (7th param)**: NOT passed (defaults to `{}`)

### Site 1 Fusion Function Assumptions

Both functions use `FusionGraphContext` to obtain SSSR rings and check disjointness.
With `ringSubsetOverride == nullptr`, they compute SSSR rings from `mol` via Indigo
(locally, not using any caller-provided ring information).

The functions' early logic checks:
- `nameLinearChainFusion`: If `ringCount < 3`, returns "Linear chain: no match" sentinel
- `nameBranchingFusion`: No explicit ringCount check, but internally checks for
  branching topology (degree > 2)

Both functions then build the `shared` matrix and `sharedNodesPairs` to identify
fusion connections, mirroring the logic in Site 1's own BFS traversal.

### Site 1 Result Handling

Results from both functions are collected into a `candidates` vector (lines 2871-2873,
2888-2890). Each result is checked:
- If `result.success`: added to `candidates`
- If `result.error` is a "no match" sentinel: skipped (fall through to next mechanism)
- If `result.error` is a specific rejection message (e.g., "All-carbon fused ring
  systems..." for linear chain): skipped (fall through)
- Otherwise: returns empty string `""` immediately (hard failure)

After both mechanisms are tried, if `candidates` is non-empty, the candidate
minimizing `rootAtom`'s locant is selected (lines 3052-3074):
```cpp
RingSystemNameResult *bestCandidate = nullptr;
int bestLocant = INT_MAX;
for (auto &candidate : candidates) {
    // Check if rootAtom has a valid numeric locant
    auto it = candidate.locantMap.find(rootAtom);
    if (it == candidate.locantMap.end()) continue;
    QString locantStr = it->second;
    // Reject letter-suffixed locants (fusion carbons can't be attachment points)
    if (!locantStr.isEmpty() && locantStr.back().isLetter()) continue;
    bool ok = false;
    int locantValue = locantStr.toInt(&ok);
    if (ok && locantValue < bestLocant) {
        bestLocant = locantValue;
        bestCandidate = &candidate;
    }
}
```

If no valid candidate found, returns empty string.

The selected name is then formatted (lines 3081-3095) and returned wrapped in parentheses
with locant and "-yl" suffix: `"(" + stem + "-" + locantStr + "-yl" + ")"`

### Site 1 Failure Handling

- Hard failures (non-sentinel errors) cause immediate `return ""`
- "No match" sentinels allow fall-through to next mechanism
- Empty candidates after all mechanisms: `return ""` (line 3049)

---

## Call Site 2: `generateNameImpl` (multi-component handling)

**Location:** `IupacNamer.cpp:23343` and `23349`

**Context:** Deep inside multi-component molecule handling, specifically for
molecules with exactly one real, non-trivial ring system (fused pair or ortho-fused
chain/branch) plus one or more separate, plain single rings elsewhere.

**Call expressions:**
```cpp
// Line 23343
auto chainRes = nameLinearChainFusion(g, mol, indigoToGraphIdx, NC, freeValAtom, 
                                       &realComponentIndigo, extraSubs);

// Line 23349
auto branchRes = nameBranchingFusion(g, mol, indigoToGraphIdx, NC, freeValAtom, 
                                       &realComponentIndigo, extraSubs);
```

### Site 2 Input Construction

- **`g`**: `Graph g`, local to `generateNameImpl` (built at line 12012 from the
  molecule's heavy atoms), NOT a function parameter -- `generateNameImpl` takes
  only `int mol` (line 8739).
- **`mol`**: the `int mol` parameter `generateNameImpl` itself receives.
- **`indigoToGraphIdx`**: `std::map<int, int> indigoToGraphIdx`, local to
  `generateNameImpl` (built at line 11784), NOT a parameter.
- **`NC`**: Set to `realComponentSize` (line 23310), which is the size of
  `realComponentIndigo` vector (line 22825: `int realComponentSize = 
  static_cast<int>(realComponentIndigo.size());`)
- **`freeValenceAtom`**: Computed as first element of `attachNodesSeen` or -1 (line 23337):
  ```cpp
  int freeValAtom = attachNodesSeen.empty() ? -1 : *attachNodesSeen.begin();
  ```
  Where `attachNodesSeen` contains the attachment points from isolated rings to the
  real component (lines 22896-22900).
- **`ringSubsetOverride`**: `&realComponentIndigo` (pointer to the vector of ring sets
  for the real component)
- **`extraSubstituents`**: `extraSubs` (line 22903: `std::vector<std::pair<int, QString>> 
  extraSubs;`), populated by naming each directly-attached isolated ring via
  `nameRingAsSubstituent` (lines 22905-22910):
  ```cpp
  for (size_t k = 0; k < isolatedAttachPairs.size(); ++k) {
      int attachReal = isolatedAttachPairs[k].first;
      int attachIso = isolatedAttachPairs[k].second;
      QString subName = nameRingAsSubstituent(g, isolatedGraphRings[k], attachIso, 
                                               attachReal, allSSSRRings, 
                                               realComponentGraphUnion, ...);
      if (subName.isEmpty()) { namesOk = false; break; }
      extraSubs.push_back({attachReal, subName});
  }
  ```

### Site 2 Fusion Function Assumptions

With `ringSubsetOverride == &realComponentIndigo`, the `FusionGraphContext` uses the
provided ring subset instead of computing SSSR from `mol`. This is critical: the
multi-component context wants to name only the REAL component (the fused system),
not the entire molecule (which includes isolated rings).

The `extraSubstituents` parameter provides pre-named substituent information for
the isolated rings attached to the real component. This allows the fusion functions
to incorporate these substituents into their numbering without having to re-detect
and re-name them.

### Site 2 Result Handling

Results are checked directly for success (lines 23344-23346, 23350-23352):
```cpp
// For linear chain
if (!chainRes.empty() && chainRes[0].success) {
    mcSuccess = true;
    mcName = chainRes[0].baseName;
}

// For branching
if (!branchRes.empty() && branchRes[0].success) {
    mcSuccess = true;
    mcName = branchRes[0].baseName;
}
```

Only the first element of the result vector is used (`chainRes[0]` or `branchRes[0]`).
The `baseName` is used directly as the component name.

If either mechanism succeeds, `mcSuccess` is set to true and `mcName` to the result,
and the function returns early (line 23361: `return {true, mcName, ""};`).

### Site 2 Failure Handling

- If neither mechanism succeeds, falls through to later checks (line 23354-23362)
- The calls are guarded by topology checks (lines 23339-23354):
  - Linear chain: only called if `validFusionC && endsC == 2 && centersC == NC - 2 && NC == 3`
  - Branching: only called if `validFusionC && hasBranchingC`
- No sentinel error checking - uses direct success/failure

---

## Call Site 3: `generateNameImpl` (retained-registry pre-check)

**Location:** `IupacNamer.cpp:23566` and `23634`

**Context:** Inside `generateNameImpl`, near the retained-registry pre-check, immediately
after the Task C InChIKey-override return falls through.

**Call expressions:**
```cpp
// Line 23566
chainResults = nameLinearChainFusion(g, mol, indigoToGraphIdx, ringCount);

// Line 23634
branchResults = nameBranchingFusion(g, mol, indigoToGraphIdx, ringCount);
```

### Site 3 Input Construction

- **`g`**: `Graph g`, local to `generateNameImpl` (built at line 12012), NOT a
  parameter -- same local as Site 2 above, since both sites are in the same
  function.
- **`mol`**: the `int mol` parameter `generateNameImpl` itself receives.
- **`indigoToGraphIdx`**: `std::map<int, int> indigoToGraphIdx`, local to
  `generateNameImpl` (built at line 11784), NOT a parameter.
- **`ringCount`**: already in scope by this point in the function -- defined
  at line 8922 as `int ringCount = indigoCountSSSR(mol);`, not recomputed at
  this call site. The call passes this existing `ringCount`; the guards at
  this site (`if (ringCount >= 3)` at line 23372) use it directly, and the
  SSSR-ring SET itself (`sssrRings`, producing `N = sssrRings.size()` at line
  23409) is a separate local used by the fusion-detection logic inside this
  block, not by the `nameLinearChainFusion`/`nameBranchingFusion` calls'
  `ringCount` argument.
- **`freeValenceAtom` (5th param)**: NOT passed (defaults to -1)
- **`ringSubsetOverride` (6th param)**: NOT passed (defaults to `nullptr`)
- **`extraSubstituents` (7th param)**: NOT passed (defaults to `{}`)

### Site 3 Fusion Function Assumptions

With all optional parameters at defaults, the functions compute SSSR from `mol`
via `FusionGraphContext` and proceed with standard logic.

These calls are guarded by topology checks computed in the same scope (lines
23372-23410):
```cpp
int N = sssrRings.size();
std::vector<std::vector<int>> shared(N, std::vector<int>(N, 0));
// ... build shared matrix ...
std::vector<int> degree(N, 0);
bool validFusion = true;
for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
        if (shared[i][j] == 2) degree[i]++;
        else if (shared[i][j] != 0) validFusion = false;
    }
}
int ends = 0, centers = 0;
bool hasBranching = false;
for (int i = 0; i < N; ++i) {
    if (degree[i] == 1) ends++;
    else if (degree[i] == 2) centers++;
    else if (degree[i] > 2) hasBranching = true;
}
```

### Site 3 Result Handling

Results are assigned to local variables and checked (lines 23566-23568, 23634-23636):
```cpp
// For linear chain
std::vector<RingSystemNameResult> chainResults;
if (validFusion && ends == 2 && centers == N - 2) {
    chainResults = nameLinearChainFusion(g, mol, indigoToGraphIdx, ringCount);
    linearChainApplies = !chainResults.empty() && 
                        chainResults[0].error != "Linear chain: no match";
}

// For branching
std::vector<RingSystemNameResult> branchResults;
if (validFusion && hasBranching) {
    branchResults = nameBranchingFusion(g, mol, indigoToGraphIdx, ringCount);
    branchingApplies = !branchResults.empty() && 
                      branchResults[0].error != "Branching fusion: no match";
}
```

Note the sentinel checking: `chainResults[0].error != "Linear chain: no match"` and
`branchResults[0].error != "Branching fusion: no match"`. This differs from Site 1,
which checks each element in the vector.

If linear chain applies and succeeds (line 23621-23623):
```cpp
// linearChainApplies guarantees chainResults[0] exists and isn't "no match"
// the only remaining case is success.
return {true, chainResults[0].baseName, ""};
```

For branching, if applies and succeeds (lines 23637-23649), returns the baseName.

Special error handling for exocyclic substituents (lines 23573-23619 for linear,
23641-23659 for branching): for the linear-chain path, if the error is
"Substituents on fused ring systems named via this mechanism are not yet
supported."; for the branching path, the equivalent check is against the
branching-specific string "Substituents on branching fused ring systems are
not yet supported." (line 23641) -- the two error strings differ, it is not
the same literal text reused for both mechanisms. Either way it performs chain-vs-ring
seniority check and may fall back to chain naming.

Other errors cause immediate return (line 23619: `return {false, "", chainResults[0].error};`)

### Site 3 Failure Handling

- Sentinel "no match" errors allow fall-through
- Specific rejection errors may trigger fallback logic
- Other errors cause immediate rejection
- Unhandled cases fall through to later mechanism checks

---

## Contract Comparison

### Parameter Passing Summary

| Parameter | Site 1 | Site 2 | Site 3 |
|-----------|--------|--------|--------|
| `g` | From function param (`nameRingSystemAsSubstituent`) | Local to `generateNameImpl` (line 12012) | Local to `generateNameImpl` (line 12012) |
| `mol` | From function param | From function param (`generateNameImpl`'s own `int mol`) | From function param (`generateNameImpl`'s own `int mol`) |
| `indigoToGraphIdx` | Built locally from `g.nodes` | Local to `generateNameImpl` (line 11784) | Local to `generateNameImpl` (line 11784) |
| `ringCount` | BFS-computed fused system size | `realComponentSize` | existing local `ringCount` (set earlier, line 8922) |
| `freeValenceAtom` | `rootAtom` (explicit) | `freeValAtom` (explicit) | -1 (default) |
| `ringSubsetOverride` | `nullptr` (default) | `&realComponentIndigo` (explicit) | `nullptr` (default) |
| `extraSubstituents` | `{}` (default) | `extraSubs` (explicit) | `{}` (default) |

### Key Differences

1. **`indigoToGraphIdx` construction**:
   - Site 1: `g` is a `const Graph &` parameter `nameRingSystemAsSubstituent`
     itself receives (line 2803); `indigoToGraphIdx` is built locally inside
     it from `g.nodes` (lines 2812-2815).
   - Site 2 & 3: `generateNameImpl` receives only `int mol` (line 8739) --
     both `g` (line 12012) and `indigoToGraphIdx` (line 11784) are locals it
     builds itself, not parameters. Sites 2 and 3 share the SAME already-
     built `g`/`indigoToGraphIdx` locals since both calls are inside the
     same function body.
   - **Impact**: Site 1's `g` arrives from its caller; Site 1 only builds its
     own `indigoToGraphIdx` locally. Sites 2 & 3 build (and share) both `g`
     and `indigoToGraphIdx` themselves, earlier in the same function.

2. **`ringSubsetOverride` usage**:
   - Site 1: Uses `nullptr` → functions compute SSSR from `mol`
   - Site 2: Uses `&realComponentIndigo` → functions use provided ring subset
   - Site 3: Uses `nullptr` → functions compute SSSR from `mol`
   - **Impact**: Site 2 constrains the functions to a specific ring subset (the real
     component); Sites 1 & 3 let functions detect all rings in `mol`

3. **`extraSubstituents` usage**:
   - Site 1: Uses `{}` → functions detect all substituents internally
   - Site 2: Uses `extraSubs` → functions incorporate pre-named isolated ring substituents
   - Site 3: Uses `{}` → functions detect all substituents internally
   - **Impact**: Site 2 provides pre-computed substituent names; Sites 1 & 3 require
     functions to name substituents themselves

4. **`freeValenceAtom` usage**:
   - Site 1: Uses `rootAtom` (the substituent attachment point)
   - Site 2: Uses first attachment point or -1
   - Site 3: Uses -1 (no specific free valence atom)
   - **Impact**: Sites 1 & 2 provide a free valence atom for the functions to skip
     in substituent detection; Site 3 does not

5. **Result handling**:
   - Site 1: Collects ALL results, filters sentinels, picks best by locant minimization
   - Site 2: Uses only `result[0]`, checks `.success` directly
   - Site 3: Uses only `result[0]`, checks both `.success` and sentinel errors

### Structural Similarities

- Sites 1 and 3 both use default `ringSubsetOverride` and `extraSubstituents`
- Sites 1 and 3 both let functions compute SSSR from `mol`
- Sites 1 and 2 both provide explicit `freeValenceAtom` values
- Sites 2 and 3 are both in `generateNameImpl` with similar topology-checking logic

### Blocking Differences

**Site 2 differs from Sites 1 and 3 in two critical ways:**

1. **`ringSubsetOverride`**: Site 2 passes `&realComponentIndigo` to constrain the
   fusion functions to only the real component's rings. Sites 1 and 3 pass `nullptr`,
   causing the functions to compute SSSR from the entire `mol` and consider all rings.
   This is a fundamental difference in scope: Site 2 is naming a specific subset of
   the molecule's rings, while Sites 1 and 3 name all rings in the molecule (or the
   provided graph).

2. **`extraSubstituents`**: Site 2 passes pre-computed substituent names for isolated
   rings attached to the real component. Sites 1 and 3 pass empty defaults, requiring
   the fusion functions to detect and name all substituents internally.

These two differences mean Site 2's contract is **fundamentally different** from
Sites 1 and 3. Site 2 is operating in a multi-component naming context where:
- Only a subset of rings (the "real component") should be named
- Isolated ring substituents have already been named and must be incorporated

In contrast, Sites 1 and 3 operate in single-ring-system contexts where:
- All rings in the relevant scope (entire graph or entire molecule) should be considered
- All substituents must be detected and named by the fusion functions themselves

### Sites 1 and 3: Can They Be Unified?

Sites 1 and 3 do NOT have the same parameter pattern -- Site 1 explicitly
passes `rootAtom` as the 5th parameter (`freeValenceAtom`); Site 3 passes
no 5th parameter and relies on the default (-1). They do agree on params
6-7 (`ringSubsetOverride`/`extraSubstituents`, both default). `g`/
`indigoToGraphIdx` are both already in scope at the call either way --
Site 1's `g` is a parameter its enclosing function (`nameRingSystemAsSubstituent`)
receives, with `indigoToGraphIdx` built locally from it; Site 3's enclosing
function (`generateNameImpl`) builds both as its own locals (see "Key
Differences" item 1 above).

**The `freeValenceAtom` difference**:
- Site 1: Passes `rootAtom` -- the attachment point for this ring system as
  a substituent, which `nameLinearChainFusion`/`nameBranchingFusion` then
  skip during substituent detection so the attachment point itself isn't
  mistaken for a substituent.
- Site 3: Passes -1 -- there is no single attachment point to exclude,
  since Site 3 is naming the whole molecule's ring system as the PARENT,
  not as a substituent on something else.

This argument-value difference is real but, on its own, is not a blocking
structural difference: a single shared entry point could simply take
`freeValenceAtom` as a required parameter and let each caller supply its
own value (`rootAtom` or -1) -- that's a cosmetic unification, not a
functional obstacle.

**The REAL blocking difference is result handling (see "Key Differences"
item 5 above), not the argument value:**
- Site 1 collects the ENTIRE result vector, filters out multiple sentinel
  strings, and picks the best candidate by lowest-locant minimization
  across however many results come back.
- Site 3 only ever looks at `result[0]`, checking both `.success` and a
  specific sentinel error string, then falls through to a separate
  chain-vs-ring seniority tiebreak and a retained-registry pre-check that
  Site 1 has no equivalent of.

A shared dispatcher would need to either (a) always do Site 1's full
collect-and-minimize handling, which Site 3's simpler `result[0]`-plus-
fallback logic was never designed to need and might silently change
behavior for, or (b) take a result-handling STRATEGY as a parameter,
which is no longer "one shared entry point" in any meaningful sense --
it's two call sites with a thin, pointless common wrapper around them.

**Conclusion**: keep Sites 1 and 3 separate. The blocking difference is
each site's own post-call result-handling logic, not the call's argument
values.

---

## Final Conclusion

**Can the call sites be unified?**

**NO** - Site 2 has a fundamentally different contract from Sites 1 and 3:

1. Site 2 passes `ringSubsetOverride` to constrain the fusion to a specific ring subset
2. Site 2 passes `extraSubstituents` with pre-named isolated ring substituents

These parameters change the semantic behavior of the fusion functions in ways
that cannot be reconciled without breaking Site 2's multi-component naming logic.

**Sites 1 and 3 have matching contracts in terms of ring detection and substituent
handling** (both use defaults for `ringSubsetOverride` and `extraSubstituents`),
but they differ in:
- The value of `freeValenceAtom` (Site 1: `rootAtom`, Site 3: -1) -- a real but
  non-blocking argument-value difference; a shared entry point could simply take
  this as a required parameter.
- Where `indigoToGraphIdx` is held: Site 1's enclosing function
  (`nameRingSystemAsSubstituent`) builds its own local from a `g` parameter it
  receives; Sites 2 and 3 share the SAME local `generateNameImpl` already built
  for itself (line 11784) -- both are locals, neither is "passed in" from
  outside; this is not itself a blocker either.
- **Their RESULT-HANDLING logic (the actual blocker, see below).**

These first two differences are **context-appropriate** rather than arbitrary:
- Site 1 needs to exclude the attachment point (`rootAtom`) from substituent detection
- Site 3 has no attachment point to exclude

But neither is what actually blocks unification. A forced unification of just
`freeValenceAtom`/`indigoToGraphIdx` handling would be straightforward (pass
`freeValenceAtom` explicitly everywhere; both sites already share the same
`indigoToGraphIdx` local shape). **The real blocker is result handling**: Site 1
collects the entire result vector, filters multiple sentinel strings, and picks
the best candidate by lowest-locant minimization; Site 3 only ever reads
`result[0]`, checks one sentinel string, and then runs a separate chain-vs-ring
seniority fallback that Site 1 has no equivalent of. A shared dispatcher would
either have to always do Site 1's full collect-and-minimize handling (which
Site 3 was never designed to need, and might silently change its behavior) or
take a result-handling strategy as a parameter -- at which point it's two call
sites with a thin, pointless wrapper, not one real shared entry point.

**Recommendation**: Keep all three call sites separate and documented as-is.

The specific blocking differences are:
- Site 2's use of `ringSubsetOverride` and `extraSubstituents` creates a fundamentally
  different naming context (multi-component with pre-named substituents) that
  cannot be unified with the single-ring-system contexts of Sites 1 and 3 without
  breaking the multi-component naming functionality.
- Between Sites 1 and 3, the different RESULT-HANDLING logic (collect-and-
  minimize vs. `result[0]`-plus-seniority-fallback) is the real obstacle --
  not the `freeValenceAtom` argument value, which is cosmetic on its own.
