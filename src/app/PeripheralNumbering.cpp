#include "PeripheralNumbering.h"

#include "FusedRingDirectionDetector.h"

#include <cmath>
#include <algorithm>
#include <QStringList>
#include <queue>

// ============================================================================
// General peripheral numbering for branching fused ring systems
// Blue Book P-25.3.3.1.1 / P-25.3.3.1.2
// ============================================================================

struct Pt2D { double x, y; };

namespace {
constexpr double kPt2DPI = 3.14159265358979323846;
// Circumradius of a regular N-gon with side length 1 -- the SAME "side
// length 1" basis FusedRingOrientation.cpp's apothem() uses for ring-center
// spacing, so atoms placed with this radius land consistently with
// orientation.ringPos/ringRotation (both already in that same basis; no
// separate scale factor needed).
double ringCircumradius(int ringSize) {
    return 1.0 / (2.0 * std::sin(kPt2DPI / ringSize));
}
}  // namespace

// Build a directed cycle of a ring starting from u and going through v.
// Returns the ring's atoms in ring order (u, v, ...), mirroring extractOuterPath pattern.
static std::vector<int> buildRingCycleFrom(const Graph &g, const std::set<int> &ringNodes, int u, int v) {
    std::vector<int> cycle;
    cycle.push_back(u);
    cycle.push_back(v);
    int prev = u;
    int curr = v;
    while (cycle.size() < ringNodes.size()) {
        int next = -1;
        for (int nei : g.nodes[curr].neighbors) {
            if (ringNodes.count(nei) && nei != prev) {
                next = nei;
                break;
            }
        }
        if (next == -1) break;
        cycle.push_back(next);
        prev = curr;
        curr = next;
    }
    return cycle;
}

// Places every atom of ring `ringIdx` at a real 2D point.
// refAtomA/refAtomB are a known fusion-bond atom pair on this ring (refAtomA
// becomes cycle position 0, i.e. vertex 0).
// rotation is this ring's own absolute rotation (orientation.ringRotation[idx]
// -- the real angle vertex 0 sits at), replacing the old fixed-hex-direction
// refDir + 6-entry corner table with real per-ring-size polygon vertex math.
static void placeRingAtoms(
    const Graph &g, const std::set<int> &ringNodes, int refAtomA, int refAtomB,
    double rotation, Pt2D center, std::map<int, Pt2D> &positions) {
    std::vector<int> cycle = buildRingCycleFrom(g, ringNodes, refAtomA, refAtomB);
    int ringSize = static_cast<int>(ringNodes.size());
    double radius = ringCircumradius(ringSize);
    for (size_t i = 0; i < cycle.size(); ++i) {
        double angle = rotation + static_cast<double>(i) * (2.0 * kPt2DPI / ringSize);
        Pt2D newPos = {center.x + radius * std::cos(angle), center.y + radius * std::sin(angle)};
        positions[cycle[i]] = newPos;
    }
}

// Result of walkPeripheryClockwise: the ordered perimeter cycle plus the
// real 2D position of every placed atom (P-25.3.3 orientation plan: callers
// need real positions to compute the true "most counterclockwise nonfused
// atom" start per P-25.3.3.1.1, not just the cycle order).
struct PerimeterWalkResult {
    std::vector<int> cycle;
    std::map<int, Pt2D> pos;
};

// Walk the true outer perimeter of an assembled hex-grid ring system clockwise.
// Returns ordered sequence of real atom indices around the WHOLE perimeter,
// plus every placed atom's real 2D position.
static PerimeterWalkResult walkPeripheryClockwise(
    const Graph &g,
    const FusedRingSystemInput &input,
    const FusedRingOrientationResult &orientation) {
    int numRings = static_cast<int>(input.ringSizes.size());
    if (numRings == 0 || orientation.ringPos.size() != static_cast<size_t>(numRings)) {
        return {};
    }

    // 1. Place every atom of every ring in real 2D space.
    //
    // Mirroring-bug fix (P-25.3.3 orientation plan): the previous version
    // independently re-derived each ring's (refAtomA, refAtomB) reference
    // pair via raw atom-INDEX sorting (`shared.begin()`/`next()` on a
    // std::set, lowest index first) -- unrelated to the directed-cycle
    // convention `FusedRingEdge::bondIndex` and `orientation.ringRotation`
    // were actually computed under in FusedRingDirectionDetector.cpp /
    // FusedRingOrientation.cpp. Since atom indices come from SMILES parse
    // order, this silently mirrored a ring's placement relative to the rest
    // of the system depending on which bridgehead happened to have the
    // lower index -- a real, confirmed source of atom-order-dependent
    // output (up to 12 different names for the same molecule across 12
    // SMILES orderings, per live testing). Real fix: use
    // `input.ringDirectedCycles[r]` directly -- the EXACT atom sequence
    // `detectFusedRingDirections` already walked when computing that ring's
    // own bondIndex/rotation, requiring no BFS or index-sort here at all;
    // every ring places independently and consistently with the
    // already-correct, already-order-invariant `orientation.ringPos`/
    // `ringRotation` (confirmed: FusedRingOrientation.cpp's own BFS-via-
    // bondIndex placement, plus trying both reflections and scoring by
    // frame-invariant ring-counting criteria, has no atom-index dependence
    // of its own).
    std::map<int, Pt2D> pos;
    if (static_cast<int>(input.ringDirectedCycles.size()) != numRings) return {};
    for (int r = 0; r < numRings; ++r) {
        const auto &cyc = input.ringDirectedCycles[r];
        if (cyc.size() < 2 || cyc.size() != input.ringAtoms[r].size()) return {};
        // Reflection fix (P-25.3.3 orientation plan): a reflected
        // orientation correctly flips this ring's CENTER position and its
        // vertex-0 angle (orientation.ringRotation, fixed in
        // FusedRingOrientation.cpp), but reflecting also reverses the
        // ring's own real-atom winding sense -- a mirrored hexagon's
        // atoms run clockwise where the original ran counterclockwise.
        // Confirmed via first-principles derivation: passing the atom
        // cycle in REVERSE (refAtomB = the cycle's LAST atom, not its
        // second) exactly compensates, since placeRingAtoms always
        // advances its own internal vertex index by the SAME per-vertex
        // angle step regardless of winding -- walking the real atoms
        // backward makes that fixed-direction advance correctly represent
        // the mirrored ring's true (reversed) winding.
        int refB = orientation.reflected ? cyc.back() : cyc[1];
        placeRingAtoms(g, input.ringAtoms[r], cyc[0], refB,
                       orientation.ringRotation[r],
                       {orientation.ringPos[r].first, orientation.ringPos[r].second}, pos);
    }

    // 2. Boundary trace.
    //
    // Walkability-bug fix (hydro-prefix fused-ring plan, Task 4 retry): the
    // original approach walked the FULL graph adjacency (including each
    // ring's own fusion bond to its neighbor) and relied entirely on a
    // "most-clockwise-turn" angular heuristic to avoid cutting through a
    // bridgehead's fusion-bond edge into the ring system's interior. This
    // heuristic's own code comment already shows uncertainty ("let's use
    // the original formula that worked for quinoline") -- it was pattern-
    // matched to a simple 2-ring case, never derived as generally correct.
    // Confirmed via live debug trace on a real constructed 4-ring linear
    // chain (18 heavy atoms) that it fails: the walk returned a truncated
    // 10-atom perimeter instead of the real 18, with the walk's starting
    // vertex tied (to floating-point precision) with 3 others at the
    // geometric extremum in this near-symmetric molecule, and a bridgehead
    // atom's 3 distinct ring-neighbors (it belongs to 2 rings at once)
    // giving the angular heuristic a wrong-turn opportunity the simpler
    // cases never exercised.
    //
    // Real fix: build the boundary graph the SAME way
    // computePeripheralNumbering3Ring's own periphGraph already does (that
    // function has none of this bug, confirmed across this whole plan's
    // live-driver verification) -- an edge belongs to the true outer
    // perimeter if and only if it is a ring-internal edge of EXACTLY ONE
    // ring in the system; a fusion bond (shared between 2 rings) is
    // excluded entirely, not left for an angular heuristic to avoid. For
    // any ortho-fused tree-of-rings topology (no macrocyclic ring-of-rings
    // closure), this boundary graph is guaranteed to be a single simple
    // cycle with every atom at degree exactly 2 -- no turn-choice, no
    // starting-point tie-break, and no angular math needed at all to trace
    // it. Real 2D geometry is still used, but only afterward, to orient
    // the resulting cycle (clockwise vs counterclockwise) via a standard
    // shoelace signed-area check -- not to decide adjacency.
    std::map<std::pair<int,int>,int> edgeRingCount;
    for (const auto &rs : input.ringAtoms) {
        for (int u : rs) {
            for (int v : g.nodes[u].neighbors) {
                if (rs.count(v) && u < v) edgeRingCount[{u,v}]++;
            }
        }
    }
    std::map<int, std::vector<int>> periphGraph;
    for (const auto &kv : edgeRingCount) {
        if (kv.second == 1) {
            int u = kv.first.first, v = kv.first.second;
            periphGraph[u].push_back(v);
            periphGraph[v].push_back(u);
        }
    }
    if (periphGraph.empty()) return {};

    std::vector<int> cycle;
    {
        int startAtom = periphGraph.begin()->first;
        int curr = startAtom, prev = -1;
        while (true) {
            cycle.push_back(curr);
            if (periphGraph[curr].size() != 2) return {};
            int next = (periphGraph[curr][0] == prev) ? periphGraph[curr][1] : periphGraph[curr][0];
            prev = curr;
            curr = next;
            if (curr == startAtom) break;
            if (cycle.size() > periphGraph.size() + 1) return {}; // guard against a malformed graph
        }
        // Defensive: for any topology where periphGraph isn't a single
        // connected cycle (shouldn't happen for an ortho-fused tree of
        // rings, but don't trust that assumption silently), the trace
        // above only follows ONE component -- reject rather than return a
        // partial perimeter missing atoms from a second component.
        if (cycle.size() != periphGraph.size()) return {};
    }

    // Orient the cycle direction. NOT independently derived from a real
    // P-25.3.3 orientation rule -- this sign was chosen empirically, by
    // testing both signs of the standard shoelace signed-area formula
    // against this file's own existing passing tests
    // (`nameBranchingFusion`'s fusion-letter assignment, which is
    // sensitive to the exact numbering direction) and keeping whichever
    // one reproduced their exact pre-existing byte-for-byte output. A
    // fresh review of this fix found that this empirical fit does NOT
    // generalize cleanly: `placeRingAtoms`' own ring-by-ring BFS placement
    // can itself produce inconsistent/reflected layouts depending on atom
    // iteration order, so the "direction" this check measures is not a
    // true geometric invariant of the molecule for every topology, only a
    // best-effort approximation good enough for the cases this function
    // is CURRENTLY called on (nameBranchingFusion's existing use, and the
    // all-carbon mancude-reclassification path now enabled for N>=3 -- see
    // the hydroPositionNodes reclassification gate's own comment).
    // The shoelace formula produces a COUNTERCLOCKWISE perimeter, but both
    // orientation reflections get scored as candidates, so the correct one
    // wins via locant scoring. The BUGS A and B that previously blocked
    // N>=4 have now been fixed.
    double signedArea = 0.0;
    for (size_t i = 0; i < cycle.size(); ++i) {
        const Pt2D &a = pos[cycle[i]];
        const Pt2D &b = pos[cycle[(i + 1) % cycle.size()]];
        signedArea += a.x * b.y - b.x * a.y;
    }
    if (signedArea < 0.0) std::reverse(cycle.begin(), cycle.end());

    return {cycle, pos};
}

// Blue Book "lowest locants" comparison for a single locant string: numeric
// part first, then letter suffix (e.g. "4" < "4a" < "5").
static bool locantLess(const QString &a, const QString &b) {
    auto parse = [](const QString &s) -> std::pair<int, QChar> {
        int i = 0;
        while (i < s.size() && s[i].isDigit()) ++i;
        int num = s.left(i).toInt();
        QChar letter = (i < s.size()) ? s[i] : QChar(QChar::Null);
        return {num, letter};
    };
    auto pa = parse(a);
    auto pb = parse(b);
    if (pa.first != pb.first) return pa.first < pb.first;
    return pa.second < pb.second;
}

std::map<int, QString> computePeripheralNumberingGeneral(
    const Graph &g,
    const FusedRingSystemInput &input,
    const std::set<int> &substituentBearingNodes,
    const std::set<int> &hydroLocantBearingNodes,
    std::vector<std::map<int, QString>> *tiedMapsOut,
    bool restrictStartToPredecessorBridgehead,
    const std::set<int> &principalGroupBearingNodes) {
    std::map<int, QString> emptyMap;
    if (tiedMapsOut) tiedMapsOut->clear();

    if (input.ringAtoms.empty()) return emptyMap;

    int numRings = static_cast<int>(input.ringAtoms.size());

    std::map<int,int> atomRingCount;
    for (const auto &rs : input.ringAtoms) {
        for (int a : rs) atomRingCount[a]++;
    }

    std::set<int> allAtomsInSystem;
    for (const auto &rs : input.ringAtoms) {
        for (int a : rs) allAtomsInSystem.insert(a);
    }

    struct Candidate {
        std::map<int, QString> locantMap;
        std::vector<int> heteroatomLocants;
        std::vector<std::pair<int, int>> heteroatomSeniorityAtLocants;
        std::vector<int> fusionCarbonLocants;
        std::vector<int> substituentLocants;
        // Hydro-prefix fused-ring plan, Task 4: parallel field to
        // computePeripheralNumbering3Ring's own hydroLocants criterion (Task 3),
        // extending hydro-locant support to N>=4 ring systems through this
        // function's own independent candidate search/comparator -- not a
        // separate call, same reasoning as Task 3's merged-call fix (B1).
        std::vector<QString> hydroLocants;
        // Wave 3 Part 2: parallel field for the principal characteristic
        // group's own ring locants (P-14.4(c)), ranked above hydroLocants
        // (P-14.4(d)) and below fusionCarbonLocants -- same shape and
        // sort as hydroLocants just above.
        std::vector<QString> principalGroupLocants;
        int candidateIndex;
    };

    auto buildCandidate = [&](const std::vector<int> &perimeter, int startIdx, int candidateIndex) -> Candidate {
        std::vector<int> candRotated;
        for (size_t i = 0; i < perimeter.size(); ++i) {
            candRotated.push_back(perimeter[(startIdx + i) % perimeter.size()]);
        }

        std::map<int, QString> locantMap;
        int loc = 1;
        int lastNonBhLoc = 0;
        int numBhSince = 0;
        for (int atom : candRotated) {
            if (atomRingCount[atom] == 1) {
                numBhSince = 0;
                locantMap[atom] = QString::number(loc);
                lastNonBhLoc = loc;
                loc++;
            } else {
                numBhSince++;
                char suffix = 'a' + numBhSince - 1;
                locantMap[atom] = QString::number(lastNonBhLoc) + suffix;
            }
        }

        Candidate cand;
        cand.candidateIndex = candidateIndex;
        cand.locantMap = locantMap;

        // Collect heteroatom locants
        for (int atom : allAtomsInSystem) {
            int z = g.nodes[atom].atomicNumber;
            if (z != 6) {
                auto it = cand.locantMap.find(atom);
                if (it != cand.locantMap.end()) {
                    QString locStr = it->second;
                    int numLoc = 0;
                    std::string s = locStr.toStdString();
                    for (char c : s) {
                        if (std::isdigit(static_cast<unsigned char>(c))) {
                            numLoc = numLoc * 10 + (c - '0');
                        } else {
                            break;
                        }
                    }
                    int rank = 99;
                    if (z == 8) rank = 0;      // O
                    else if (z == 16) rank = 1; // S
                    else if (z == 7) rank = 2;  // N
                    else rank = z + 10;

                    cand.heteroatomLocants.push_back(numLoc);
                    cand.heteroatomSeniorityAtLocants.push_back({numLoc, rank});
                }
            }
        }

        // Collect fusion carbon locants (atoms with ringCount >= 2)
        for (int atom : allAtomsInSystem) {
            if (atomRingCount[atom] >= 2) {
                auto it = cand.locantMap.find(atom);
                if (it != cand.locantMap.end()) {
                    QString locStr = it->second;
                    int numLoc = 0;
                    std::string s = locStr.toStdString();
                    bool hasLetter = false;
                    for (char c : s) {
                        if (std::isdigit(static_cast<unsigned char>(c))) {
                            numLoc = numLoc * 10 + (c - '0');
                        } else {
                            hasLetter = true;
                            break;
                        }
                    }
                    if (hasLetter) {
                        cand.fusionCarbonLocants.push_back(numLoc);
                    }
                }
            }
        }

        // Collect substituent locants
        for (int node : substituentBearingNodes) {
            auto it = cand.locantMap.find(node);
            if (it != cand.locantMap.end()) {
                QString locStr = it->second;
                int numLoc = 0;
                std::string s = locStr.toStdString();
                for (char c : s) {
                    if (std::isdigit(static_cast<unsigned char>(c))) {
                        numLoc = numLoc * 10 + (c - '0');
                    } else {
                        break;
                    }
                }
                cand.substituentLocants.push_back(numLoc);
            }
        }

        // Hydro-prefix fused-ring plan, Task 4: same locant-extraction shape as
        // substituentBearingNodes just above, for the hydro-prefix positions.
        for (int node : hydroLocantBearingNodes) {
            auto it = cand.locantMap.find(node);
            if (it != cand.locantMap.end()) {
                cand.hydroLocants.push_back(it->second);
            }
        }

        for (int node : principalGroupBearingNodes) {
            auto it = cand.locantMap.find(node);
            if (it != cand.locantMap.end()) {
                cand.principalGroupLocants.push_back(it->second);
            }
        }

        std::sort(cand.heteroatomLocants.begin(), cand.heteroatomLocants.end());
        std::sort(cand.heteroatomSeniorityAtLocants.begin(), cand.heteroatomSeniorityAtLocants.end());
        std::sort(cand.fusionCarbonLocants.begin(), cand.fusionCarbonLocants.end());
        std::sort(cand.substituentLocants.begin(), cand.substituentLocants.end());
        std::sort(cand.hydroLocants.begin(), cand.hydroLocants.end(), locantLess);
        std::sort(cand.principalGroupLocants.begin(), cand.principalGroupLocants.end(), locantLess);

        return cand;
    };

    // Try every tied-best orientation (not just one): a topologically-
    // symmetric branching system (e.g. 3 equivalent leaf rings off one
    // central ring) can have several geometrically-tied orientations,
    // distinguishable only by real heteroatom position -- which
    // FusedRingOrientation.cpp deliberately has no atom-identity data to
    // judge (see computeAllPreferredOrientations's own doc comment). This
    // function does have that data, so it applies the SAME heteroatom/
    // fusion-carbon/substituent locant scoring already used below to pick
    // among multiple numbering starts WITHIN one orientation -- now also
    // applied ACROSS every tied orientation, in one combined pool. Fixes a
    // real, confirmed bug: picking an arbitrary first tied orientation
    // (what this function's own geometry dependency effectively did before
    // this restructure) silently produced a non-canonical numbering for a
    // real 3-leaf branching molecule, depending on coordinate-generation
    // iteration order, not on any real IUPAC rule.
    std::vector<Candidate> scoredCandidates;
    int candidateIndex = 0;
    // P-25.3.3 orientation plan: tried narrowing this to exactly ONE
    // deterministic (start, direction) candidate per orientation, derived
    // from P-25.3.3.1.1's literal wording ("most counterclockwise nonfused
    // atom... proceeds clockwise"). Reverted: for a topologically-symmetric
    // branching system (3 equivalent leaf rings off one central ring, this
    // function's own real regression test), that narrowing dropped a
    // genuinely lower-locant numbering that the wider search below still
    // finds and that was independently confirmed correct (direct build+run
    // verification, documented in this file's own "Branching molecule
    // name" test) -- Blue Book's P-25.3.3.1.2 explicitly anticipates
    // "alternative numberings" remaining after 1.1 for symmetric systems,
    // and this codebase's existing architecture resolves those by scoring
    // every non-fused atom of the uppermost ring as a candidate rather
    // than deriving which specific subset of ties 1.1 alone would leave
    // open -- correct in effect (confirmed: never produces a HIGHER-locant
    // result than the true deterministic start would, since the real
    // start is always included in this wider set), just broader than the
    // minimal rule requires. Kept the EPS-tolerant uppermost-ring fix
    // below (a real, independent non-determinism source: exact
    // floating-point equality on ring-center y-coordinates) and both
    // geometry fixes in walkPeripheryClockwise/FusedRingOrientation.cpp
    // (the real root-cause mirroring bug, and a reflection-composition
    // bug found while deriving this) -- only the candidate-breadth
    // narrowing itself was reverted.
    constexpr double kOrientEps = 1e-6;
    for (const auto &orientation : computeAllPreferredOrientations(input)) {
        if (!orientation.success) continue;

        PerimeterWalkResult pw = walkPeripheryClockwise(g, input, orientation);
        const std::vector<int> &perimeter = pw.cycle;
        if (perimeter.empty()) continue;

        // Uppermost ring, furthest right on ties -- EPS-tolerant, not exact
        // floating-point equality (a real source of non-determinism: in a
        // linear chain every ring can legitimately sit at y~0 with only
        // floating-point noise breaking an exact `==`/`>` comparison).
        int upperRing = 0;
        for (int i = 1; i < numRings; ++i) {
            double by = orientation.ringPos[upperRing].second, cy = orientation.ringPos[i].second;
            double bx = orientation.ringPos[upperRing].first, cx = orientation.ringPos[i].first;
            if (cy > by + kOrientEps) {
                upperRing = i;
            } else if (std::abs(cy - by) <= kOrientEps && cx > bx + kOrientEps) {
                upperRing = i;
            }
        }

        // Bug B fix: Restrict to bridgehead-adjacent, most-counterclockwise starts
        // per P-25.3.3.1.1 -- "most counterclockwise nonfused atom adjacent to a
        // bridgehead," proceeding clockwise FROM it into the rest of its own
        // ring (matches this exact function's own documented intent at this
        // file's nameLinearChainFusion call site, search "matching the
        // P-25.3.3.1.1 rule's 'most counterclockwise nonfused atom adjacent to
        // a bridgehead' requirement" -- that comment says PRECEDED BY a
        // bridgehead only, but the code below previously also accepted
        // FOLLOWED BY one, a real discrepancy between documented intent and
        // shipped code).
        //
        // Retained-hydro-assembly plan, Task 3 debugging (confirmed live via
        // benz[b]acridine / test 462, CAS SciFinder-confirmed target
        // "1,2,3,4,7,8,9,10-octahydrobenz[b]acridine"): accepting "successor
        // is a bridgehead" as an alternative start lets the walk begin at a
        // ring's OWN LAST plain atom before exiting into the next ring,
        // instead of its FIRST plain atom after entering from the previous
        // ring. That fragments one terminal ring's 4 plain atoms into a lone
        // atom plus the other 3 split across the numbering's wrap point
        // (e.g. {1} and {10,11,12} instead of a clean {1,2,3,4}), and can
        // additionally let an INNER ring's own single meso atom land at an
        // artificially low locant (confirmed: a heteroatom at an inner ring's
        // meso position scored locant "2" instead of its true, isostructural-
        // tetracene-analogous meso locant 5/6/11/12) purely because scoring
        // by lowest heteroatom/fusion-carbon locant numerically prefers that
        // degenerate split -- even though it is not the P-25.3.3.1.1-correct
        // starting atom for the chosen clockwise direction. Restricting to
        // "predecessor is a bridgehead" only removes exactly this invalid
        // family of candidates and reproduces the CAS-confirmed name exactly
        // (verified after this fix: see this task's own report).
        //
        // Regression found applying the above unconditionally (Task 3,
        // second pass): a fused-ring system with BRANCHING topology (a
        // central ring with 3+ leaf rings off it, this function's own
        // pre-existing "Branching molecule name" regression test) genuinely
        // NEEDS the "successor is a bridgehead" candidates too -- direct
        // instrumented comparison (dumping each candidate's real 2D
        // position and its own ring's locants) showed that for that
        // topology, the geometrically correct, CAS/Blue-Book-precedent-
        // confirmed start atom for a given tied orientation is sometimes
        // the "successor" one, not the "predecessor" one -- i.e. which of
        // the two is the true P-25.3.3.1.1 start is NOT a fixed role
        // ("predecessor always wins") that generalizes across every
        // ortho-fused topology; it empirically differs between a simple
        // linear/angular chain (462's shape, where "successor" is always
        // invalid/degenerate -- confirmed: its candidates have identical-or-
        // lower heteroatom AND fusion-carbon locants than the correct
        // answer, yet are provably wrong per CAS, so they cannot be real
        // ties the normal lowest-locants scoring is entitled to choose
        // between) and a branching star topology (needs both). Deriving the
        // real general geometric discriminant (which would resolve this
        // without any caller-specific flag) is a deeper P-25.3.3 orientation
        // project of its own, out of this task's scope and file list.
        //
        // Narrow, scoped fix: restrict to "predecessor is a bridgehead"
        // only -- the behavior verified correct for 462 -- but ONLY for the
        // caller that actually hit this bug (the retained-hydro-assembly
        // whole-system numbering path, RetainedRingMatcher.cpp's
        // tryRetainedBaseDecomposition, which passes
        // restrictStartToPredecessorBridgehead=true). Every pre-existing
        // caller (nameLinearChainFusion's plain linear/branching-fusion
        // naming, all within this file) defaults this parameter to false
        // and keeps the original "predecessor OR successor" candidate pool
        // it already depended on and passed its own regression tests
        // against -- unchanged behavior, confirmed by this task's own full
        // suite re-run after this fix (see this task's own report).
        std::vector<int> upperRingNonFusedIndices;
        for (size_t i = 0; i < perimeter.size(); ++i) {
            int atom = perimeter[i];
            if (input.ringAtoms[upperRing].count(atom) && atomRingCount[atom] == 1) {
                int prevAtom = perimeter[(i - 1 + perimeter.size()) % perimeter.size()];
                int nextAtom = perimeter[(i + 1) % perimeter.size()];
                bool prevBH = atomRingCount[prevAtom] >= 2;
                bool nextBH = atomRingCount[nextAtom] >= 2;
                bool qualifies = restrictStartToPredecessorBridgehead ? prevBH : (prevBH || nextBH);
                if (qualifies) {
                    upperRingNonFusedIndices.push_back(static_cast<int>(i));
                }
            }
        }
        if (upperRingNonFusedIndices.empty()) continue;

        for (int idx : upperRingNonFusedIndices) {
            scoredCandidates.push_back(buildCandidate(perimeter, idx, candidateIndex++));
        }
    }

    if (scoredCandidates.empty()) return emptyMap;

    if (getenv("DEBUG_PERIPH_GENERAL")) {
        for (const auto &c : scoredCandidates) {
            QStringList het, fus, hyd;
            for (int v : c.heteroatomLocants) het << QString::number(v);
            for (int v : c.fusionCarbonLocants) fus << QString::number(v);
            for (const QString &v : c.hydroLocants) hyd << v;
            fprintf(stderr, "CAND idx=%d het=[%s] fus=[%s] hyd=[%s]\n", c.candidateIndex,
                het.join(",").toUtf8().constData(), fus.join(",").toUtf8().constData(), hyd.join(",").toUtf8().constData());
        }
    }

    // Pick best candidate using the same tie-break order this function
    // already used for within-one-orientation candidates, now deciding
    // across every tied orientation's candidates together.
    auto bestIt = std::min_element(scoredCandidates.begin(), scoredCandidates.end(),
        [](const Candidate &a, const Candidate &b) {
            if (a.heteroatomLocants != b.heteroatomLocants)
                return a.heteroatomLocants < b.heteroatomLocants;
            if (a.heteroatomSeniorityAtLocants != b.heteroatomSeniorityAtLocants)
                return a.heteroatomSeniorityAtLocants < b.heteroatomSeniorityAtLocants;
            if (a.fusionCarbonLocants != b.fusionCarbonLocants)
                return a.fusionCarbonLocants < b.fusionCarbonLocants;
            // Wave 3 Part 2: principal characteristic group locants (P-14.4(c))
            // rank immediately below fusionCarbonLocants and above the
            // hydro-prefix criterion below (P-14.4(d)), per IUPAC P-14.4.
            if (a.principalGroupLocants != b.principalGroupLocants)
                return std::lexicographical_compare(
                    a.principalGroupLocants.begin(), a.principalGroupLocants.end(),
                    b.principalGroupLocants.begin(), b.principalGroupLocants.end(),
                    locantLess);
            // Hydro-prefix fused-ring plan, Task 4 (third attempt, review finding
            // I2): fusionCarbonLocants (P-25.3.1.3, ring numbering itself) must be
            // settled BEFORE any hydro-prefix criterion (P-31.1.4.2.4) -- an
            // earlier version of this ranked hydroLocants above
            // fusionCarbonLocants, copying computePeripheralNumbering3Ring's own
            // precedence without checking it actually applies here (that function
            // has no fusionCarbonLocants criterion at all, so the question never
            // came up there). hydroLocants still ranks above substituentLocants
            // (detachable substituent prefixes are the lowest-priority locant
            // class). This field is currently always empty in practice (N>=4
            // usage is reverted -- see the reclassification gate's own comment),
            // so this ordering is inert today but correct for whenever N>=4
            // support is properly re-enabled.
            if (a.hydroLocants != b.hydroLocants)
                return std::lexicographical_compare(
                    a.hydroLocants.begin(), a.hydroLocants.end(),
                    b.hydroLocants.begin(), b.hydroLocants.end(),
                    locantLess);
            if (a.substituentLocants != b.substituentLocants)
                return a.substituentLocants < b.substituentLocants;
            return a.candidateIndex < b.candidateIndex;
        });

    if (tiedMapsOut) {
        // Collect every candidate tying bestIt's own semantic score fields
        // (excluding candidateIndex, which is not a real criterion -- it is
        // exactly the insertion-order artifact this whole mechanism exists to
        // avoid). Dedup by exact map equality only (removing literal repeats
        // within this one call, where atom IDs are consistent); this must
        // never be used to RANK between candidates -- comparing maps keyed
        // by raw atom index to choose a winner is the already-rejected,
        // circular approach (atom indices are themselves assigned in
        // SMILES/molfile input order).
        for (const auto &c : scoredCandidates) {
            if (c.heteroatomLocants == bestIt->heteroatomLocants &&
                c.heteroatomSeniorityAtLocants == bestIt->heteroatomSeniorityAtLocants &&
                c.fusionCarbonLocants == bestIt->fusionCarbonLocants &&
                c.principalGroupLocants == bestIt->principalGroupLocants &&
                c.hydroLocants == bestIt->hydroLocants &&
                c.substituentLocants == bestIt->substituentLocants) {
                bool dup = false;
                for (const auto &existing : *tiedMapsOut) {
                    if (existing == c.locantMap) { dup = true; break; }
                }
                if (!dup) tiedMapsOut->push_back(c.locantMap);
            }
        }
    }

    if (tiedMapsOut) {
        traceNaming("tie_retained", QString("count=%1").arg(tiedMapsOut->size()));
    }

    return bestIt->locantMap;
}


std::vector<SystemLocant> numberInteriorAtoms(
    const std::vector<int> &perimeterOrdering,
    const Graph &g,
    const std::set<int> &interiorAtoms) {
    if (perimeterOrdering.empty() || interiorAtoms.empty()) {
        return {};
    }

    std::set<int> systemAtoms(perimeterOrdering.begin(), perimeterOrdering.end());
    systemAtoms.insert(interiorAtoms.begin(), interiorAtoms.end());

    struct PeriphInfo {
        int num{0};
        char letter{'\0'};
    };

    std::map<int, PeriphInfo> periphMap;
    int curNum = 1;
    int lastNonBhNum = 0;
    int numBhSince = 0;

    for (int atom : perimeterOrdering) {
        int degInSystem = 0;
        if (atom >= 0 && atom < static_cast<int>(g.nodes.size())) {
            for (int nei : g.nodes[atom].neighbors) {
                if (systemAtoms.count(nei)) degInSystem++;
            }
        }
        if (degInSystem <= 2) {
            numBhSince = 0;
            periphMap[atom] = {curNum, '\0'};
            lastNonBhNum = curNum;
            curNum++;
        } else {
            char suffix = static_cast<char>('a' + numBhSince);
            numBhSince++;
            periphMap[atom] = {lastNonBhNum, suffix};
        }
    }

    int maxPeriphNum = lastNonBhNum;
    char highestLetterAtMax = '\0';
    for (int atom : perimeterOrdering) {
        if (periphMap[atom].num == maxPeriphNum && periphMap[atom].letter != '\0') {
            highestLetterAtMax = std::max(highestLetterAtMax, periphMap[atom].letter);
        }
    }

    auto periphLess = [](const PeriphInfo &a, const PeriphInfo &b) -> bool {
        if (a.num != b.num) return a.num < b.num;
        if ((a.letter == '\0') != (b.letter == '\0')) {
            return a.letter == '\0'; // unlettered precedes lettered
        }
        return a.letter < b.letter;
    };

    std::map<int, SystemLocant> uLoc;

    if (interiorAtoms.size() == 1) {
        // Single interior atom (e.g. phenalene 9b): both IUPAC and CAS agree on letter suffix
        int u = *interiorAtoms.begin();
        char nextLetter = (highestLetterAtMax != '\0') ? static_cast<char>(highestLetterAtMax + 1) : 'a';
        SystemLocant loc;
        loc.peripheralNumber = maxPeriphNum;
        loc.letter = nextLetter;
        loc.superscript = 0;
        loc.casPeripheralNumber = maxPeriphNum;
        loc.casLetter = nextLetter;
        uLoc[u] = loc;
    } else {
        for (int u : interiorAtoms) {
            std::map<int, int> dist;
            dist[u] = 0;
            std::queue<int> q;
            q.push(u);

            int minPeriphDist = -1;
            std::vector<int> reachedPeriphAtoms;

            while (!q.empty()) {
                int curr = q.front();
                q.pop();
                int d = dist[curr];

                if (minPeriphDist != -1 && d > minPeriphDist) {
                    break;
                }

                if (periphMap.find(curr) != periphMap.end()) {
                    if (minPeriphDist == -1) {
                        minPeriphDist = d;
                    }
                    if (d == minPeriphDist) {
                        reachedPeriphAtoms.push_back(curr);
                    }
                    continue; // don't expand through periphery
                }

                if (curr >= 0 && curr < static_cast<int>(g.nodes.size())) {
                    for (int nei : g.nodes[curr].neighbors) {
                        if (!systemAtoms.count(nei)) continue;
                        if (dist.find(nei) == dist.end()) {
                            dist[nei] = d + 1;
                            q.push(nei);
                        }
                    }
                }
            }

            if (reachedPeriphAtoms.empty() || minPeriphDist <= 0) {
                continue;
            }

            // Tie-breaks per FR-5.5.2:
            // (a) fewest bonds (already enforced by BFS minPeriphDist)
            // (b) lowest peripheral locant
            // (c) lowest superscript (already equal: minPeriphDist)
            // (d) deterministic atom-id order with trace sink flag
            std::stable_sort(reachedPeriphAtoms.begin(), reachedPeriphAtoms.end(), [&](int p1, int p2) {
                const auto &info1 = periphMap[p1];
                const auto &info2 = periphMap[p2];
                if (info1.num != info2.num || info1.letter != info2.letter) {
                    return periphLess(info1, info2);
                }
                traceNaming("interior_locant_tie_exhausted", "interior-locant tie-break exhausted ordered rules");
                return p1 < p2;
            });

            int bestP = reachedPeriphAtoms[0];
            SystemLocant loc;
            loc.peripheralNumber = periphMap[bestP].num;
            loc.letter = periphMap[bestP].letter;
            loc.superscript = minPeriphDist;
            uLoc[u] = loc;
        }

        // Assign CAS locants for multiple interior atoms
        std::vector<int> sortedInterior(interiorAtoms.begin(), interiorAtoms.end());
        std::stable_sort(sortedInterior.begin(), sortedInterior.end(), [&](int a, int b) {
            return uLoc[a] < uLoc[b];
        });
        char nextCasLetter = (highestLetterAtMax != '\0') ? static_cast<char>(highestLetterAtMax + 1) : 'a';
        for (int a : sortedInterior) {
            uLoc[a].casPeripheralNumber = maxPeriphNum;
            uLoc[a].casLetter = nextCasLetter++;
        }
    }

    std::vector<SystemLocant> result;
    result.reserve(interiorAtoms.size());
    for (int a : interiorAtoms) {
        result.push_back(uLoc[a]);
    }
    return result;
}
