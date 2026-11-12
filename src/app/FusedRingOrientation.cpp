#include "FusedRingOrientation.h"
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>

namespace {
constexpr double PI = 3.14159265358979323846;

// Apothem (center-to-edge-midpoint distance) of a regular N-gon with side
// length 1. Circumradius isn't needed -- only apothem-sum placement is used.
double apothem(int ringSize) {
    return 1.0 / (2.0 * std::tan(PI / ringSize));
}

// The absolute angle (radians, standard math convention: 0 = +x axis,
// counterclockwise positive) of the midpoint of a ring's bondIndex-th edge,
// given the ring's own absolute rotation (the angle its own vertex 0 sits
// at). Edge k's midpoint sits exactly between vertex k and vertex k+1, each
// vertex spaced 2*PI/ringSize apart.
double edgeMidpointAngle(double ringRotation, int ringSize, CycleEdgeIndex bondIndex) {
    double vertexSpacing = 2.0 * PI / ringSize;
    return ringRotation + (bondIndex.value + 0.5) * vertexSpacing;
}
}  // namespace

struct PlacedRing {
    double x = 0.0, y = 0.0;
    double rotation = 0.0;  // absolute angle of this ring's own vertex 0
    bool placed = false;
};

struct ScoreTuple {
    int maxRow;
    double ur;
    double inv_ll; 
    double above;
    
    bool operator<(const ScoreTuple& o) const {
        if (maxRow != o.maxRow) return maxRow < o.maxRow;
        if (ur != o.ur) return ur < o.ur;
        if (inv_ll != o.inv_ll) return inv_ll < o.inv_ll;
        return above < o.above;
    }
};

std::vector<FusedRingOrientationResult> computeAllPreferredOrientations(const FusedRingSystemInput &input) {
    FusedRingOrientationResult res;
    if (!input.error.empty()) {
        res.success = false;
        res.error = input.error;
        return {res};
    }
    if (input.ringSizes.empty()) {
        res.success = true;
        return {res};
    }

    int num_rings = static_cast<int>(input.ringSizes.size());
    if (input.fusions.size() >= static_cast<size_t>(num_rings)) {
        res.success = false;
        res.error = "Non-tree ring systems (ortho-and-peri-fused like phenalene/pyrene) are not yet supported.";
        return {res};
    }

    // adj[ring] = list of (neighborRing, bondIndexWithinRing, isOutgoing)
    // -- isOutgoing true means THIS ring's bondIndex is the one recorded in
    // input.fusions (ring1==this ring); false means we're looking at it from
    // ring2's side and need the mirrored angle (edge midpoint + PI).
    struct AdjEntry { int neighbor; CycleEdgeIndex bondIndex; bool isOutgoing; };
    std::vector<std::vector<AdjEntry>> adj(num_rings);
    for (const auto &edge : input.fusions) {
        adj[edge.ring1].push_back({edge.ring2, edge.bondIndex, true});
        adj[edge.ring2].push_back({edge.ring1, edge.bondIndex, false});
    }

    std::vector<int> inDegree(num_rings, 0);
    for (const auto &edge : input.fusions) {
        if (edge.ring2 >= 0 && edge.ring2 < num_rings) {
            inDegree[edge.ring2]++;
        }
    }
    int rootRing = 0;
    for (int i = 0; i < num_rings; ++i) {
        if (inDegree[i] == 0) {
            rootRing = i;
            break;
        }
    }

    std::vector<PlacedRing> placed(num_rings);
    placed[rootRing] = {0.0, 0.0, 0.0, true};

    std::queue<int> q;
    q.push(rootRing);
    int visited_count = 1;

    const double TOLERANCE = 1e-6;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (const auto &nbr : adj[u]) {
            int v = nbr.neighbor;
            if (placed[v].placed) continue;
            double edgeAngle;
            if (nbr.isOutgoing) {
                edgeAngle = edgeMidpointAngle(placed[u].rotation, input.ringSizes[u], nbr.bondIndex);
            } else {
                edgeAngle = edgeMidpointAngle(placed[u].rotation, input.ringSizes[u], nbr.bondIndex) + PI;
            }
            double dist = apothem(input.ringSizes[u]) + apothem(input.ringSizes[v]);
            double vx = placed[u].x + dist * std::cos(edgeAngle);
            double vy = placed[u].y + dist * std::sin(edgeAngle);
            double vRotation = edgeAngle + PI - 0.5 * (2.0 * PI / input.ringSizes[v]);
            placed[v] = {vx, vy, vRotation, true};
            visited_count++;
            q.push(v);
        }
    }

    if (visited_count != num_rings) {
        res.success = false;
        res.error = "Disconnected ring system";
        return {res};
    }

    for (int i = 0; i < num_rings; i++) {
        for (int j = i + 1; j < num_rings; j++) {
            double dx = placed[i].x - placed[j].x;
            double dy = placed[i].y - placed[j].y;
            double distSq = dx * dx + dy * dy;
            double minExpectedDist = std::min(apothem(input.ringSizes[i]), apothem(input.ringSizes[j])) * 0.5;
            if (distSq < minExpectedDist * minExpectedDist) {
                res.success = false;
                res.error = "Geometric contradiction: Rings overlap";
                return {res};
            }
        }
    }

    // Candidate orientations: for each real fusion bond, the rotation that
    // aligns that bond's direction with the horizontal axis, times 2 for a
    // reflection (flip y), times 2 more for a 180-degree rotation -- a small,
    // chemically-meaningful candidate set (not continuous/arbitrary rotation
    // search), generalizing the old code's 12 discrete hex symmetries.
    //
    // NOTE: a broader "every vertex of every ring" candidate set was tried
    // here and made things WORSE (regressed naphthalene/anthracene/
    // naphthacene, which were correct with the narrower set) -- the
    // branching-ground-truth regression's real cause is NOT candidate
    // breadth; it's elsewhere (scoring criteria or tie-break behavior).
    // Reverted to this narrower, empirically-correct-for-linear-chains set.
    //
    // Fix for Bug A: Each fusion bond's base rotation now gets TWO partners:
    // (1) its Y-reflection (already present), and (2) its 180-degree X-mirror
    // rotation. For a linear chain, this ensures BOTH terminal rings can be
    // "uppermost" regardless of BFS root choice, eliminating atom-order
    // dependence on which end gets locant 1.
    std::vector<double> candidateRotations;
    for (const auto &edge : input.fusions) {
        double a = edgeMidpointAngle(placed[edge.ring1].rotation, input.ringSizes[edge.ring1], edge.bondIndex);
        candidateRotations.push_back(-a);
        // Add 180-degree (PI) rotation partner for this fusion bond
        candidateRotations.push_back(-a + PI);
    }
    if (candidateRotations.empty()) candidateRotations.push_back(0.0);

    struct CandidateResult {
        ScoreTuple score;
        std::vector<std::pair<double,double>> pos;
        std::vector<double> rotation;
        double y_main, x_mid;
        bool reflected;
    };
    std::vector<CandidateResult> allCandidates;

    for (double rot : candidateRotations) {
        for (int reflect = 0; reflect < 2; ++reflect) {
            std::vector<std::pair<double,double>> cur_pos(num_rings);
            std::vector<double> cur_rotation(num_rings);
            for (int i = 0; i < num_rings; i++) {
                double x = placed[i].x * std::cos(rot) - placed[i].y * std::sin(rot);
                double y = placed[i].x * std::sin(rot) + placed[i].y * std::cos(rot);
                if (reflect) y = -y;
                cur_pos[i] = {x, y};
                // Each ring's own vertex-0 angle must undergo the EXACT
                // same point transform (rotate by rot, then reflect) as the
                // ring's own relative position vector would -- NOT "negate
                // first, then add rot" (that was the previous formula here,
                // confirmed wrong by direct derivation during the P-25.3.3
                // orientation plan): reflecting a point at angle theta maps
                // it to angle -theta, so reflect(rotate(angle, rot)) =
                // -(angle + rot) = -angle - rot, not -angle + rot. Verified
                // the non-reflected branch (rot alone, no sign issue) was
                // already correct; only the reflected branch's sign on
                // `rot` was wrong.
                cur_rotation[i] = reflect ? (-placed[i].rotation - rot) : (placed[i].rotation + rot);
            }

            // Group rings into horizontal rows (same y within tolerance) to
            // find the largest contiguous row -- mirrors the old code's
            // r_groups-by-hex-row logic, generalized to real y with a
            // tolerance instead of exact integer r equality.
            std::vector<std::pair<double,double>> sortedByY = cur_pos;
            std::sort(sortedByY.begin(), sortedByY.end(), [](auto &a, auto &b){ return a.second < b.second; });

            int max_block_size = 0;
            double block_y = 0.0;
            for (size_t i = 0; i < sortedByY.size(); ) {
                size_t j = i;
                while (j < sortedByY.size() && std::abs(sortedByY[j].second - sortedByY[i].second) < 0.1) j++;
                int blockSize = static_cast<int>(j - i);
                if (blockSize > max_block_size) {
                    max_block_size = blockSize;
                    block_y = sortedByY[i].second;
                }
                i = j;
            }

            double y_main = block_y;
            double x_mid = 0.0;
            {
                double sum = 0.0; int cnt = 0;
                for (auto &p : cur_pos) {
                    if (std::abs(p.second - y_main) < 0.1) { sum += p.first; cnt++; }
                }
                if (cnt > 0) x_mid = sum / cnt;
            }

            double ur = 0, ll = 0, above = 0;
            for (int i = 0; i < num_rings; i++) {
                double y_rel = cur_pos[i].second - y_main;
                double x_rel = cur_pos[i].first - x_mid;
                bool above_axis = (y_rel > 0.1), below_axis = (y_rel < -0.1), bisect_h = (std::abs(y_rel) <= 0.1);
                bool right_axis = (x_rel > 0.1), left_axis = (x_rel < -0.1), bisect_v = (std::abs(x_rel) <= 0.1);

                if (above_axis && right_axis) ur += 1.0;
                else if (bisect_h && right_axis) ur += 0.5;
                else if (above_axis && bisect_v) ur += 0.5;
                else if (bisect_h && bisect_v) ur += 0.25;

                if (below_axis && left_axis) ll += 1.0;
                else if (bisect_h && left_axis) ll += 0.5;
                else if (below_axis && bisect_v) ll += 0.5;
                else if (bisect_h && bisect_v) ll += 0.25;

                if (above_axis) above += 1.0;
                else if (bisect_h) above += 0.5;
            }

            ScoreTuple cand_score = {max_block_size, ur, -ll, above};
            allCandidates.push_back({cand_score, cur_pos, cur_rotation, y_main, x_mid, reflect != 0});
        }
    }

    // Find the best score, then collect EVERY candidate that ties it (within
    // floating-point tolerance) -- not just the first one found. For a
    // topologically-symmetric system (e.g. 3 equivalent leaf rings off one
    // central ring), multiple genuinely-tied orientations exist and only a
    // heteroatom-aware tie-break (which this purely-geometric module has no
    // atom-identity data to perform) can correctly choose among them. Real,
    // confirmed bug this fixes: picking an arbitrary "first" tied candidate
    // here silently chose a different (and wrong, by the already-hand-
    // verified ground truth) orientation than before this refactor, for a
    // real 3-leaf branching molecule -- not a hex-grid-specific issue, a
    // pre-existing gap this is the first change to actually expose. Callers
    // that have real atom/heteroatom data (IupacNamer.cpp) use
    // computeAllPreferredOrientations to get every tied candidate and pick
    // correctly; computePreferredOrientation (below) stays as a thin
    // backward-compatible wrapper returning just the first for callers that
    // don't need tie-aware behavior.
    ScoreTuple best_score = {-1, -1e9, -1e9, -1e9};
    for (const auto &c : allCandidates) {
        if (best_score < c.score) best_score = c.score;
    }

    auto scoreTies = [](const ScoreTuple &a, const ScoreTuple &b) {
        const double eps = 1e-6;
        return a.maxRow == b.maxRow && std::abs(a.ur - b.ur) < eps &&
               std::abs(a.inv_ll - b.inv_ll) < eps && std::abs(a.above - b.above) < eps;
    };

    std::vector<FusedRingOrientationResult> results;
    for (const auto &c : allCandidates) {
        if (!scoreTies(c.score, best_score)) continue;

        FusedRingOrientationResult r;
        r.success = true;
        r.ringsInHorizontalRow = best_score.maxRow;
        r.ringsInUpperRightQuadrant = best_score.ur;
        r.ringsInLowerLeftQuadrant = -best_score.inv_ll;
        r.ringsAboveHorizontalRow = best_score.above;

        std::vector<int> sorted_rings(num_rings);
        for (int i = 0; i < num_rings; i++) sorted_rings[i] = i;
        std::sort(sorted_rings.begin(), sorted_rings.end(), [&](int a, int b) {
            double yA = c.pos[a].second, yB = c.pos[b].second;
            double xA = c.pos[a].first, xB = c.pos[b].first;
            bool a_in_main = (std::abs(yA - c.y_main) <= 0.1);
            bool b_in_main = (std::abs(yB - c.y_main) <= 0.1);
            if (a_in_main != b_in_main) return a_in_main > b_in_main;
            if (a_in_main) return xA < xB;
            if (std::abs(yA - yB) > 0.1) return yA > yB;
            return xA < xB;
        });

        r.ringOrder = sorted_rings;
        r.ringPos = c.pos;
        r.ringRotation = c.rotation;
        r.reflected = c.reflected;
        results.push_back(r);
    }
    return results;
}

FusedRingOrientationResult computePreferredOrientation(const FusedRingSystemInput &input) {
    std::vector<FusedRingOrientationResult> all = computeAllPreferredOrientations(input);
    if (all.empty()) {
        FusedRingOrientationResult res;
        if (input.ringSizes.empty()) { res.success = true; return res; }
        res.success = false;
        res.error = "Non-tree ring systems (ortho-and-peri-fused like phenalene/pyrene) are not yet supported.";
        return res;
    }
    return all[0];
}
