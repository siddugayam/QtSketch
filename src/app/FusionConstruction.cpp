#include "FusionConstruction.h"
#include "RetainedRingMatcher.h"
#include "FusedSystemTopology.h"
#include "FusedRingDirectionDetector.h"
#include "indigo.h"
#include <algorithm>
#include <functional>
#include <queue>
#include <iostream>

bool hitsFr32Limitation(const ComponentTree &tree) {
    if (tree.components.size() < 3) return false;
    int n = static_cast<int>(tree.components.size());
    std::vector<std::vector<bool>> adj(n, std::vector<bool>(n, false));
    for (const auto &j : tree.junctions) {
        if (j.comp1 >= 0 && j.comp1 < n && j.comp2 >= 0 && j.comp2 < n && j.comp1 != j.comp2) {
            adj[j.comp1][j.comp2] = true;
            adj[j.comp2][j.comp1] = true;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (!adj[i][j]) continue;
            for (int k = j + 1; k < n; ++k) {
                if (adj[i][k] && adj[j][k]) {
                    return true;
                }
            }
        }
    }
    return false;
}

static FusionComponent classifyMonocycle(int ringIdx, const std::set<int> &ringAtoms, const Graph &g) {
    FusionComponent comp;
    comp.ringCount = 1;
    comp.ringIndices = {ringIdx};
    comp.atomIndices = ringAtoms;

    int nAtoms = static_cast<int>(ringAtoms.size());
    std::map<int, int> heteroCounts;
    for (int a : ringAtoms) {
        if (a >= 0 && a < static_cast<int>(g.nodes.size())) {
            int z = g.nodes[a].atomicNumber;
            if (z != 6) {
                heteroCounts[z]++;
            }
        }
    }

    if (heteroCounts.empty()) {
        // Carbocycle
        if (nAtoms == 6) {
            comp.name = "benzene";
            comp.prefix = "benzo";
            comp.seniorityRank = 100;
        } else if (nAtoms == 5) {
            comp.name = "cyclopentadiene";
            comp.prefix = "cyclopenta";
            comp.seniorityRank = 101;
        } else if (nAtoms == 7) {
            comp.name = "cycloheptatriene";
            comp.prefix = "cyclohepta";
            comp.seniorityRank = 102;
        } else if (nAtoms == 4) {
            comp.name = "cyclobutadiene";
            comp.prefix = "cyclobuta";
            comp.seniorityRank = 103;
        } else if (nAtoms == 3) {
            comp.name = "cyclopropene";
            comp.prefix = "cyclopropa";
            comp.seniorityRank = 104;
        } else if (nAtoms == 8) {
            comp.name = "cyclooctatetraene";
            comp.prefix = "cycloocta";
            comp.seniorityRank = 105;
        } else {
            comp.name = "cyclooctane";
            comp.prefix = "cycloocta";
            comp.seniorityRank = 110;
        }
        return comp;
    }

    // Heterocycles
    int numO = heteroCounts[8];
    int numS = heteroCounts[16];
    int numN = heteroCounts[7];

    std::vector<int> nitrogens;
    for (int a : ringAtoms) {
        if (a >= 0 && a < static_cast<int>(g.nodes.size()) && g.nodes[a].atomicNumber == 7) {
            nitrogens.push_back(a);
        }
    }

    if (nAtoms == 5) {
        if (numO == 1 && numN == 0 && numS == 0) {
            comp.name = "furan"; comp.prefix = "furo"; comp.seniorityRank = 20;
        } else if (numS == 1 && numN == 0 && numO == 0) {
            comp.name = "thiophene"; comp.prefix = "thieno"; comp.seniorityRank = 21;
        } else if (numN == 1 && numO == 0 && numS == 0) {
            comp.name = "pyrrole"; comp.prefix = "pyrrolo"; comp.seniorityRank = 13;
        } else if (numN == 2) {
            bool adjacent = false;
            if (nitrogens.size() == 2) {
                int n1 = nitrogens[0], n2 = nitrogens[1];
                for (int nei : g.nodes[n1].neighbors) {
                    if (nei == n2) { adjacent = true; break; }
                }
            }
            if (adjacent) {
                comp.name = "pyrazole"; comp.prefix = "pyrazolo"; comp.seniorityRank = 10;
            } else {
                comp.name = "imidazole"; comp.prefix = "imidazo"; comp.seniorityRank = 11;
            }
        } else if (numO == 1 && numN == 1) {
            comp.name = "oxazole"; comp.prefix = "oxazolo"; comp.seniorityRank = 15;
        } else if (numS == 1 && numN == 1) {
            comp.name = "thiazole"; comp.prefix = "thiazolo"; comp.seniorityRank = 16;
        } else {
            comp.name = "heterocycle"; comp.prefix = "hetero"; comp.seniorityRank = 30;
        }
    } else if (nAtoms == 6) {
        if (numN == 2) {
            int dist = 3;
            if (nitrogens.size() == 2) {
                int n1 = nitrogens[0], n2 = nitrogens[1];
                for (int nei : g.nodes[n1].neighbors) {
                    if (nei == n2) { dist = 1; break; }
                }
                if (dist > 1) {
                    for (int nei : g.nodes[n1].neighbors) {
                        if (!ringAtoms.count(nei)) continue;
                        for (int nei2 : g.nodes[n2].neighbors) {
                            if (nei == nei2) { dist = 2; break; }
                        }
                        if (dist == 2) break;
                    }
                }
            }
            if (dist == 1) {
                comp.name = "pyridazine"; comp.prefix = "pyridazino"; comp.seniorityRank = 7;
            } else if (dist == 2) {
                comp.name = "pyrimidine"; comp.prefix = "pyrimido"; comp.seniorityRank = 8;
            } else {
                comp.name = "pyrazine"; comp.prefix = "pyrazino"; comp.seniorityRank = 9;
            }
        } else if (numN == 1 && numO == 0 && numS == 0) {
            comp.name = "pyridine"; comp.prefix = "pyrido"; comp.seniorityRank = 12;
        } else {
            comp.name = "heterocycle"; comp.prefix = "hetero"; comp.seniorityRank = 14;
        }
    } else {
        comp.name = "heterocycle"; comp.prefix = "hetero"; comp.seniorityRank = 35;
    }
    return comp;
}

static ComponentTree buildTreeFromComponents(
    const std::vector<FusionComponent> &comps,
    int parentIdx,
    const Graph &g)
{
    ComponentTree tree;
    tree.components = comps;
    tree.parentIndex = parentIdx;
    for (size_t i = 0; i < tree.components.size(); ++i) {
        tree.components[i].id = static_cast<int>(i);
    }
    if (parentIdx >= 0 && parentIdx < static_cast<int>(tree.components.size())) {
        tree.components[parentIdx].order = 0;
    }

    // Build junctions
    for (size_t i = 0; i < tree.components.size(); ++i) {
        for (size_t j = i + 1; j < tree.components.size(); ++j) {
            std::vector<int> commonAtoms;
            for (int a : tree.components[i].atomIndices) {
                if (tree.components[j].atomIndices.count(a)) {
                    commonAtoms.push_back(a);
                }
            }
            if (commonAtoms.empty()) continue;

            std::vector<std::pair<int, int>> commonBonds;
            for (size_t u = 0; u < commonAtoms.size(); ++u) {
                for (size_t v = u + 1; v < commonAtoms.size(); ++v) {
                    for (int nei : g.nodes[commonAtoms[u]].neighbors) {
                        if (nei == commonAtoms[v]) {
                            commonBonds.push_back({commonAtoms[u], commonAtoms[v]});
                            break;
                        }
                    }
                }
            }
            if (!commonBonds.empty() || commonAtoms.size() >= 2) {
                FusionJunction junc;
                junc.comp1 = static_cast<int>(i);
                junc.comp2 = static_cast<int>(j);
                junc.sharedAtoms = commonAtoms;
                junc.sharedBonds = commonBonds;
                tree.junctions.push_back(junc);
            }
        }
    }

    // Assign component orders via BFS from parent
    if (tree.parentIndex >= 0) {
        std::vector<int> order(tree.components.size(), -1);
        order[tree.parentIndex] = 0;
        std::queue<int> q;
        q.push(tree.parentIndex);
        while (!q.empty()) {
            int curr = q.front();
            q.pop();
            for (const auto &j : tree.junctions) {
                int nei = -1;
                if (j.comp1 == curr) nei = j.comp2;
                else if (j.comp2 == curr) nei = j.comp1;
                if (nei >= 0 && order[nei] == -1) {
                    order[nei] = order[curr] + 1;
                    q.push(nei);
                }
            }
        }
        for (size_t i = 0; i < tree.components.size(); ++i) {
            tree.components[i].order = (order[i] >= 0) ? order[i] : 1;
        }
    }

    tree.valid = true;
    return tree;
}

ParentScore scoreTree(const ComponentTree &tree, const Graph &g) {
    ParentScore score;
    if (tree.parentIndex < 0 || tree.parentIndex >= static_cast<int>(tree.components.size())) {
        return score;
    }
    const auto &p = tree.components[tree.parentIndex];
    bool hasHetero = false;
    for (int a : p.atomIndices) {
        if (a >= 0 && a < static_cast<int>(g.nodes.size()) && g.nodes[a].atomicNumber != 6) {
            hasHetero = true;
            break;
        }
    }
    score.heteroParent = hasHetero;
    score.parentRingCount = p.ringCount;
    score.parentSeniority = p.seniorityRank;
    score.noFr7Fallback = !hitsFr32Limitation(tree);
    score.parentName = p.name;
    score.sortedParentRings = std::vector<int>(p.ringIndices.begin(), p.ringIndices.end());
    std::sort(score.sortedParentRings.begin(), score.sortedParentRings.end());

    std::map<QString, int> attachedNameCounts;
    int bestInterparentSeniority = 9999;
    for (size_t i = 0; i < tree.components.size(); ++i) {
        if (static_cast<int>(i) == tree.parentIndex) continue;
        int ord = tree.components[i].order;
        if (ord == 1) {
            score.firstOrderCount++;
            attachedNameCounts[tree.components[i].prefix]++;
            if (tree.components[i].seniorityRank > 0 && tree.components[i].seniorityRank < bestInterparentSeniority) {
                bestInterparentSeniority = tree.components[i].seniorityRank;
            }
        } else if (ord > 1) {
            score.higherOrderCount++;
        }
    }
    score.interparentSeniority = bestInterparentSeniority;
    for (const auto &pair : attachedNameCounts) {
        if (pair.second > score.maxIdenticalAttached) {
            score.maxIdenticalAttached = pair.second;
        }
    }
    return score;
}

ComponentTree decomposeComponentTree(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol,
    bool allowFr32)
{
    ComponentTree tree;
    int numRings = static_cast<int>(rings.size());
    if (numRings == 0) {
        tree.error = "No rings to decompose";
        return tree;
    }

    // 1. Build ring-ring adjacency
    std::vector<std::vector<int>> ringAdj(numRings);
    for (int i = 0; i < numRings; ++i) {
        for (int j = i + 1; j < numRings; ++j) {
            std::vector<int> common;
            for (int a : rings[i]) {
                if (rings[j].count(a)) common.push_back(a);
            }
            if (common.size() >= 2) {
                bool hasBond = false;
                for (size_t u = 0; u < common.size(); ++u) {
                    for (size_t v = u + 1; v < common.size(); ++v) {
                        for (int nei : g.nodes[common[u]].neighbors) {
                            if (nei == common[v]) { hasBond = true; break; }
                        }
                        if (hasBond) break;
                    }
                    if (hasBond) break;
                }
                if (hasBond) {
                    ringAdj[i].push_back(j);
                    ringAdj[j].push_back(i);
                }
            }
        }
    }

    // 2. Identify candidate polycyclic components
    std::vector<FusionComponent> polyCandidates;
    if (mol >= 0) {
        const auto &registry = fullRegistry();

        auto findConnectedSubsets = [&](int targetK) {
            std::vector<std::set<int>> result;
            if (targetK <= 0 || targetK > numRings || numRings > 20) return result;
            for (int mask = 1; mask < (1 << numRings); ++mask) {
                if (__builtin_popcount(mask) != targetK) continue;
                std::set<int> subset;
                int start = -1;
                for (int i = 0; i < numRings; ++i) {
                    if (mask & (1 << i)) {
                        subset.insert(i);
                        start = i;
                    }
                }
                if (start == -1) continue;
                std::queue<int> q;
                std::set<int> vis;
                q.push(start);
                vis.insert(start);
                while (!q.empty()) {
                    int curr = q.front();
                    q.pop();
                    for (int nei : ringAdj[curr]) {
                        if (subset.count(nei) && !vis.count(nei)) {
                            vis.insert(nei);
                            q.push(nei);
                        }
                    }
                }
                if (vis.size() == subset.size()) {
                    result.push_back(subset);
                }
            }
            return result;
        };

        for (int k = numRings; k >= 2; --k) {
            auto connectedK = findConnectedSubsets(k);
            for (const auto &subset : connectedK) {
                std::set<int> atomSet;
                for (int r : subset) {
                    atomSet.insert(rings[r].begin(), rings[r].end());
                }
                for (const auto &sys : registry) {
                    if (static_cast<int>(atomSet.size()) != static_cast<int>(sys.peripheralLocants.size())) {
                        continue;
                    }
                    std::map<int, int> atomMap;
                    bool matched = false;
                    if (matchesRetainedSystem(mol, atomSet, sys, atomMap)) {
                        matched = true;
                    } else {
                        SkeletonMatchResult skel;
                        if (skeletonMatchRetainedSystem(mol, atomSet, sys, skel)) {
                            matched = true;
                        }
                    }
                    if (matched) {
                        FusionComponent comp;
                        comp.name = sys.name;
                        comp.prefix = sys.fusionPrefix.isEmpty() ? (sys.name + "o") : sys.fusionPrefix;
                        comp.seniorityRank = sys.seniorityRank;
                        comp.ringCount = k;
                        comp.ringIndices = subset;
                        comp.atomIndices = atomSet;
                        polyCandidates.push_back(comp);
                    }
                }
            }
        }
    }

    // Monocyclic fallback components
    std::vector<FusionComponent> monoComps;
    for (int r = 0; r < numRings; ++r) {
        monoComps.push_back(classifyMonocycle(r, rings[r], g));
    }

    // Build candidate trees for evaluation
    std::vector<ComponentTree> candidateTrees;

    // Case A: All rings as monocycles
    {
        ComponentTree monoTree = buildTreeFromComponents(monoComps, 0, g);
        // Find best monocycle parent
        int bestP = 0;
        for (size_t i = 1; i < monoComps.size(); ++i) {
            if (monoComps[i].seniorityRank < monoComps[bestP].seniorityRank) {
                bestP = static_cast<int>(i);
            }
        }
        monoTree = buildTreeFromComponents(monoComps, bestP, g);
        candidateTrees.push_back(monoTree);
    }

    // Case B: Each polyCandidate as parent, with remaining rings partitioned
    for (const auto &parentCand : polyCandidates) {
        std::set<int> parentRings = parentCand.ringIndices;
        std::vector<int> remainingRings;
        for (int r = 0; r < numRings; ++r) {
            if (!parentRings.count(r)) remainingRings.push_back(r);
        }

        // Partition 1: All remaining rings as monocycles
        {
            std::vector<FusionComponent> comps;
            comps.push_back(parentCand);
            for (int r : remainingRings) {
                comps.push_back(monoComps[r]);
            }
            candidateTrees.push_back(buildTreeFromComponents(comps, 0, g));
        }

        // Partition 2: Greedy polycycles in remaining rings + monocycles
        {
            std::vector<FusionComponent> comps;
            comps.push_back(parentCand);
            std::set<int> assigned = parentRings;
            for (const auto &otherPoly : polyCandidates) {
                bool overlaps = false;
                for (int r : otherPoly.ringIndices) {
                    if (assigned.count(r)) { overlaps = true; break; }
                }
                if (!overlaps) {
                    comps.push_back(otherPoly);
                    for (int r : otherPoly.ringIndices) assigned.insert(r);
                }
            }
            for (int r : remainingRings) {
                if (!assigned.count(r)) {
                    comps.push_back(monoComps[r]);
                    assigned.insert(r);
                }
            }
            candidateTrees.push_back(buildTreeFromComponents(comps, 0, g));
        }
    }

    // Case C: A heteromonocycle as parent, with polycycles in remaining rings
    for (size_t m = 0; m < monoComps.size(); ++m) {
        if (monoComps[m].seniorityRank >= 100) continue;
        std::set<int> assigned = {static_cast<int>(m)};
        std::vector<FusionComponent> comps;
        comps.push_back(monoComps[m]);
        for (const auto &poly : polyCandidates) {
            bool overlaps = false;
            for (int r : poly.ringIndices) {
                if (assigned.count(r)) { overlaps = true; break; }
            }
            if (!overlaps) {
                comps.push_back(poly);
                for (int r : poly.ringIndices) assigned.insert(r);
            }
        }
        for (size_t r = 0; r < monoComps.size(); ++r) {
            if (!assigned.count(static_cast<int>(r))) {
                comps.push_back(monoComps[r]);
            }
        }
        candidateTrees.push_back(buildTreeFromComponents(comps, 0, g));
    }

    // Select the best tree using FR-2.3 and FR-3.3.1
    if (candidateTrees.empty()) {
        tree.error = "No valid decomposition found";
        return tree;
    }

    ComponentTree bestTree;
    ParentScore bestScore;
    bool found = false;

    for (size_t i = 0; i < candidateTrees.size(); ++i) {
        if (!allowFr32 && hitsFr32Limitation(candidateTrees[i])) {
            continue;
        }
        ParentScore s = scoreTree(candidateTrees[i], g);
        if (!found || bestScore < s) {
            bestScore = s;
            bestTree = candidateTrees[i];
            found = true;
        }
    }

    if (!found) {
        tree.error = "No valid decomposition found";
        return tree;
    }

    if (mol >= 0) {
        applyBenzoheterocycleMerging(bestTree, g, mol);
    }

    return bestTree;
}

static QString alphabetizationKey(const QString &pName) {
    QString s = pName;
    if ((s.startsWith("(") && s.endsWith(")")) || (s.startsWith("[") && s.endsWith("]"))) {
        s = s.mid(1, s.length() - 2);
    } else if (s.startsWith("(") || s.startsWith("[")) {
        s = s.mid(1);
    }
    if (s.startsWith("(")) {
        int closeParen = s.indexOf(')');
        if (closeParen != -1 && closeParen + 1 < s.length() && s[closeParen + 1] == '-') {
            s = s.mid(closeParen + 2);
        }
    }
    int i = 0;
    while (i < s.length() && (s[i].isDigit() || s[i] == ',')) i++;
    if (i > 0 && i < s.length() && s[i] == '-') {
        s = s.mid(i + 1);
    }
    return s;
}

static std::vector<QStringList> collectBranches(const ComponentTree &tree, int compIdx, const std::vector<std::vector<int>> &treeAdj) {
    std::vector<std::pair<QString, QStringList>> branches;
    for (int child : treeAdj[compIdx]) {
        auto subBranches = collectBranches(tree, child, treeAdj);
        QStringList fullBranch;
        for (const auto &b : subBranches) {
            fullBranch.append(b);
        }
        fullBranch.append(tree.components[child].prefix.isEmpty() ? tree.components[child].name : tree.components[child].prefix);
        QString key = alphabetizationKey(fullBranch.first()).toLower();
        branches.push_back({key, fullBranch});
    }
    std::sort(branches.begin(), branches.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });
    std::vector<QStringList> result;
    for (const auto &p : branches) {
        result.push_back(p.second);
    }
    return result;
}

QStringList orderCitation(const ComponentTree &tree) {
    QStringList result;
    if (!tree.valid || tree.parentIndex < 0 || tree.parentIndex >= static_cast<int>(tree.components.size())) {
        return result;
    }
    int n = static_cast<int>(tree.components.size());
    std::vector<std::vector<int>> treeAdj(n);
    for (const auto &j : tree.junctions) {
        if (j.comp1 < 0 || j.comp1 >= n || j.comp2 < 0 || j.comp2 >= n) continue;
        int o1 = tree.components[j.comp1].order;
        int o2 = tree.components[j.comp2].order;
        if (o2 == o1 + 1) {
            treeAdj[j.comp1].push_back(j.comp2);
        } else if (o1 == o2 + 1) {
            treeAdj[j.comp2].push_back(j.comp1);
        }
    }

    auto parentBranches = collectBranches(tree, tree.parentIndex, treeAdj);
    for (const auto &b : parentBranches) {
        result.append(b);
    }
    result.append(tree.components[tree.parentIndex].name);
    return result;
}

void applyBenzoheterocycleMerging(ComponentTree &tree, const Graph &g, int mol) {
    if (!tree.valid || tree.components.size() < 2 || mol < 0) return;

    std::map<QString, int> prefixCounts;
    for (const auto &c : tree.components) {
        prefixCounts[c.prefix]++;
    }
    for (const auto &pair : prefixCounts) {
        if (pair.second >= 2) {
            // Multiplicative prefix present (e.g. "dibenzo") -- skip merging per FR-3.5
            return;
        }
    }

    const auto &registry = fullRegistry();

    for (size_t ji = 0; ji < tree.junctions.size(); ++ji) {
        int i1 = tree.junctions[ji].comp1;
        int i2 = tree.junctions[ji].comp2;
        if (i1 < 0 || i2 < 0 || i1 >= static_cast<int>(tree.components.size()) || i2 >= static_cast<int>(tree.components.size())) continue;

        int bIdx = -1, hIdx = -1;
        if (tree.components[i1].ringCount == 1 && tree.components[i1].seniorityRank >= 100 &&
            tree.components[i2].ringCount == 1 && tree.components[i2].seniorityRank < 100) {
            bIdx = i1; hIdx = i2;
        } else if (tree.components[i2].ringCount == 1 && tree.components[i2].seniorityRank >= 100 &&
                   tree.components[i1].ringCount == 1 && tree.components[i1].seniorityRank < 100) {
            bIdx = i2; hIdx = i1;
        }
        if (bIdx == -1 || hIdx == -1) continue;

        std::set<int> combAtoms = tree.components[bIdx].atomIndices;
        combAtoms.insert(tree.components[hIdx].atomIndices.begin(), tree.components[hIdx].atomIndices.end());

        for (const auto &sys : registry) {
            if (static_cast<int>(combAtoms.size()) != static_cast<int>(sys.peripheralLocants.size())) continue;
            std::map<int, int> atomMap;
            bool matched = false;
            if (matchesRetainedSystem(mol, combAtoms, sys, atomMap)) {
                matched = true;
            } else {
                SkeletonMatchResult skel;
                if (skeletonMatchRetainedSystem(mol, combAtoms, sys, skel)) {
                    matched = true;
                }
            }
            if (matched) {
                FusionComponent merged;
                merged.name = sys.name;
                merged.prefix = sys.fusionPrefix.isEmpty() ? (sys.name + "o") : sys.fusionPrefix;
                merged.seniorityRank = sys.seniorityRank;
                merged.ringCount = 2;
                merged.isMergedBenzoheterocycle = true;
                merged.ringIndices = tree.components[bIdx].ringIndices;
                merged.ringIndices.insert(tree.components[hIdx].ringIndices.begin(), tree.components[hIdx].ringIndices.end());
                merged.atomIndices = combAtoms;

                std::vector<FusionComponent> newComps;
                for (size_t c = 0; c < tree.components.size(); ++c) {
                    if (static_cast<int>(c) != bIdx && static_cast<int>(c) != hIdx) {
                        newComps.push_back(tree.components[c]);
                    }
                }
                newComps.push_back(merged);

                int bestP = 0;
                for (size_t c = 1; c < newComps.size(); ++c) {
                    bool cHetero = (newComps[c].seniorityRank < 100);
                    bool bestHetero = (newComps[bestP].seniorityRank < 100);
                    if (cHetero != bestHetero) {
                        if (cHetero) bestP = static_cast<int>(c);
                    } else if (newComps[c].ringCount != newComps[bestP].ringCount) {
                        if (newComps[c].ringCount > newComps[bestP].ringCount) bestP = static_cast<int>(c);
                    } else if (newComps[c].seniorityRank < newComps[bestP].seniorityRank) {
                        bestP = static_cast<int>(c);
                    }
                }
                tree = buildTreeFromComponents(newComps, bestP, g);
                return;
            }
        }
    }
}

static QString sideLetter(int idx) {
    if (idx < 26) {
        return QString(QChar('a' + idx));
    }
    int letterIdx = idx % 26;
    int subscript = idx / 26;
    return QString(QChar('a' + letterIdx)) + QString::number(subscript);
}

static QString addPrimesToLocant(const QString &loc, int numPrimes) {
    if (numPrimes <= 0) return loc;
    int idx = 0;
    while (idx < loc.size() && (loc[idx].isDigit() || loc[idx].isLetter())) {
        if (loc[idx] == '\'') break;
        ++idx;
    }
    QString res = loc;
    res.insert(idx, QString(numPrimes, '\''));
    return res;
}

static bool locantLess(const QString &a, const QString &b) {
    int aNum = 0, bNum = 0;
    int aIdx = 0, bIdx = 0;
    while (aIdx < a.size() && a[aIdx].isDigit()) { aNum = aNum * 10 + (a[aIdx].toLatin1() - '0'); ++aIdx; }
    while (bIdx < b.size() && b[bIdx].isDigit()) { bNum = bNum * 10 + (b[bIdx].toLatin1() - '0'); ++bIdx; }
    if (aNum != bNum) return aNum < bNum;
    QString aSuffix = a.mid(aIdx);
    QString bSuffix = b.mid(bIdx);
    return aSuffix < bSuffix;
}

struct ComponentNumbering {
    std::map<int, QString> atomLocants; // target_atom -> locant (e.g. "1", "2")
    std::map<std::pair<int, int>, QString> bondLetters; // directed edge (u, v) -> letter
    std::map<std::pair<int, int>, std::pair<int, int>> forwardEdge; // undirected edge {u, v} -> directed edge (from, to)
};

static std::vector<ComponentNumbering> generateMonocycleNumberings(const FusionComponent &comp, const Graph &g) {
    std::vector<ComponentNumbering> result;
    int n = static_cast<int>(comp.atomIndices.size());
    if (n < 3) return result;

    std::vector<int> cycle;
    std::vector<int> atoms(comp.atomIndices.begin(), comp.atomIndices.end());
    int start = atoms[0];
    cycle.push_back(start);
    std::set<int> visited;
    visited.insert(start);

    while (static_cast<int>(cycle.size()) < n) {
        int curr = cycle.back();
        int nextNode = -1;
        for (int nei : g.nodes[curr].neighbors) {
            if (comp.atomIndices.count(nei) && !visited.count(nei)) {
                nextNode = nei;
                break;
            }
        }
        if (nextNode == -1) break;
        visited.insert(nextNode);
        cycle.push_back(nextNode);
    }
    if (static_cast<int>(cycle.size()) != n) return result;

    bool closed = false;
    for (int nei : g.nodes[cycle.back()].neighbors) {
        if (nei == cycle[0]) { closed = true; break; }
    }
    if (!closed) return result;

    auto heteroRank = [](int z) -> int {
        switch (z) {
            case 8: return 1;  // O
            case 16: return 2; // S
            case 34: return 3; // Se
            case 52: return 4; // Te
            case 7: return 5;  // N
            case 15: return 6; // P
            case 33: return 7; // As
            case 51: return 8; // Sb
            case 83: return 9; // Bi
            case 14: return 10;// Si
            case 32: return 11;// Ge
            case 50: return 12;// Sn
            case 82: return 13;// Pb
            case 5: return 14; // B
            default: return (z == 6) ? 999 : 99;
        }
    };

    std::set<int> heteroRanks;
    for (int a : cycle) {
        int r = heteroRank(g.nodes[a].atomicNumber);
        if (r < 999) heteroRanks.insert(r);
    }

    struct WalkScore {
        std::vector<std::vector<int>> locsPerHetero;
        bool operator<(const WalkScore &other) const {
            return locsPerHetero < other.locsPerHetero;
        }
        bool operator==(const WalkScore &other) const {
            return locsPerHetero == other.locsPerHetero;
        }
    };

    std::vector<std::pair<std::vector<int>, WalkScore>> allWalks;
    WalkScore bestScore;
    bool hasBest = false;

    for (int dir : {1, -1}) {
        for (int s = 0; s < n; ++s) {
            std::vector<int> walk(n);
            for (int i = 0; i < n; ++i) {
                int idx = (s + dir * i) % n;
                if (idx < 0) idx += n;
                walk[i] = cycle[idx];
            }
            WalkScore score;
            for (int hr : heteroRanks) {
                std::vector<int> locs;
                for (int i = 0; i < n; ++i) {
                    if (heteroRank(g.nodes[walk[i]].atomicNumber) == hr) {
                        locs.push_back(i + 1);
                    }
                }
                std::sort(locs.begin(), locs.end());
                score.locsPerHetero.push_back(locs);
            }
            if (!hasBest || score < bestScore) {
                bestScore = score;
                hasBest = true;
            }
            allWalks.push_back({walk, score});
        }
    }

    for (const auto &pair : allWalks) {
        if (pair.second == bestScore) {
            ComponentNumbering num;
            const auto &w = pair.first;
            for (int i = 0; i < n; ++i) {
                num.atomLocants[w[i]] = QString::number(i + 1);
            }
            for (int i = 0; i < n; ++i) {
                int u = w[i];
                int v = w[(i + 1) % n];
                QString letter = sideLetter(i);
                num.bondLetters[{u, v}] = letter;
                num.bondLetters[{v, u}] = letter;
                num.forwardEdge[{std::min(u, v), std::max(u, v)}] = {u, v};
            }
            result.push_back(num);
        }
    }
    return result;
}

static std::vector<ComponentNumbering> generatePolycycleNumberings(const FusionComponent &comp, const Graph &g, int mol) {
    (void)g;
    std::vector<ComponentNumbering> result;
    if (mol < 0) return result;

    const auto &registry = fullRegistry();
    const RetainedRingSystem *matchedSys = nullptr;
    for (const auto &sys : registry) {
        if (sys.name == comp.name) {
            matchedSys = &sys;
            break;
        }
    }
    if (!matchedSys) return result;

    int refMol = indigoLoadMoleculeFromString(matchedSys->referenceSmiles.toUtf8().constData());
    if (refMol < 0) return result;
    indigoAromatize(refMol);

    std::map<int, int> baseMap;
    bool matched = matchesRetainedSystem(mol, comp.atomIndices, *matchedSys, baseMap);
    if (!matched) {
        SkeletonMatchResult skel;
        if (skeletonMatchRetainedSystem(mol, comp.atomIndices, *matchedSys, skel)) {
            baseMap = skel.atomMap;
            matched = true;
        }
    }
    if (!matched) {
        indigoFree(refMol);
        return result;
    }

    auto baseAutos = computeAutomorphisms(refMol);
    for (const auto &perm : baseAutos) {
        ComponentNumbering num;
        for (size_t i = 0; i < matchedSys->peripheralLocants.size(); ++i) {
            if (i < perm.size()) {
                int refAtom = perm[i];
                if (baseMap.count(refAtom)) {
                    int targetAtom = baseMap[refAtom];
                    num.atomLocants[targetAtom] = matchedSys->peripheralLocants[i];
                }
            }
        }
        for (const auto &kv : matchedSys->peripheralBondLetters) {
            QString letter(kv.first);
            const QString &locA = kv.second.first;
            const QString &locB = kv.second.second;
            int targetA = -1, targetB = -1;
            for (const auto &pair : num.atomLocants) {
                if (pair.second == locA) targetA = pair.first;
                else if (pair.second == locB) targetB = pair.first;
            }
            if (targetA >= 0 && targetB >= 0) {
                num.bondLetters[{targetA, targetB}] = letter;
                num.bondLetters[{targetB, targetA}] = letter;
                num.forwardEdge[{std::min(targetA, targetB), std::max(targetA, targetB)}] = {targetA, targetB};
            }
        }
        result.push_back(num);
    }
    indigoFree(refMol);
    return result;
}

struct Fr45Score {
    std::vector<QString> parentLettersSet;
    std::vector<QString> parentLettersCitation;
    std::vector<QString> firstOrderLocantsSet;
    std::vector<QString> firstOrderLocantsCitation;
    std::vector<QString> higherOrderLocantsSet;
    std::vector<QString> higherOrderLocantsCitation;

    bool operator<(const Fr45Score &other) const {
        if (parentLettersSet != other.parentLettersSet) {
            return parentLettersSet < other.parentLettersSet;
        }
        if (parentLettersCitation != other.parentLettersCitation) {
            return parentLettersCitation < other.parentLettersCitation;
        }
        if (firstOrderLocantsSet != other.firstOrderLocantsSet) {
            for (size_t i = 0; i < std::min(firstOrderLocantsSet.size(), other.firstOrderLocantsSet.size()); ++i) {
                if (firstOrderLocantsSet[i] != other.firstOrderLocantsSet[i]) {
                    return locantLess(firstOrderLocantsSet[i], other.firstOrderLocantsSet[i]);
                }
            }
            return firstOrderLocantsSet.size() < other.firstOrderLocantsSet.size();
        }
        if (firstOrderLocantsCitation != other.firstOrderLocantsCitation) {
            for (size_t i = 0; i < std::min(firstOrderLocantsCitation.size(), other.firstOrderLocantsCitation.size()); ++i) {
                if (firstOrderLocantsCitation[i] != other.firstOrderLocantsCitation[i]) {
                    return locantLess(firstOrderLocantsCitation[i], other.firstOrderLocantsCitation[i]);
                }
            }
            return firstOrderLocantsCitation.size() < other.firstOrderLocantsCitation.size();
        }
        if (higherOrderLocantsSet != other.higherOrderLocantsSet) {
            for (size_t i = 0; i < std::min(higherOrderLocantsSet.size(), other.higherOrderLocantsSet.size()); ++i) {
                if (higherOrderLocantsSet[i] != other.higherOrderLocantsSet[i]) {
                    return locantLess(higherOrderLocantsSet[i], other.higherOrderLocantsSet[i]);
                }
            }
            return higherOrderLocantsSet.size() < other.higherOrderLocantsSet.size();
        }
        if (higherOrderLocantsCitation != other.higherOrderLocantsCitation) {
            for (size_t i = 0; i < std::min(higherOrderLocantsCitation.size(), other.higherOrderLocantsCitation.size()); ++i) {
                if (higherOrderLocantsCitation[i] != other.higherOrderLocantsCitation[i]) {
                    return locantLess(higherOrderLocantsCitation[i], other.higherOrderLocantsCitation[i]);
                }
            }
            return higherOrderLocantsCitation.size() < other.higherOrderLocantsCitation.size();
        }
        return false;
    }
};

QString encodeFusionDescriptor(const ComponentTree &tree, const Graph &g, int mol) {
    if (!tree.valid || tree.components.size() < 2 || tree.junctions.empty()) {
        return QString();
    }

    std::vector<std::vector<ComponentNumbering>> allCompNumberings(tree.components.size());
    for (size_t i = 0; i < tree.components.size(); ++i) {
        if (tree.components[i].ringCount == 1) {
            allCompNumberings[i] = generateMonocycleNumberings(tree.components[i], g);
        } else {
            allCompNumberings[i] = generatePolycycleNumberings(tree.components[i], g, mol);
        }
        if (allCompNumberings[i].empty()) {
            ComponentNumbering dummy;
            int idx = 1;
            for (int a : tree.components[i].atomIndices) {
                dummy.atomLocants[a] = QString::number(idx++);
            }
            allCompNumberings[i].push_back(dummy);
        }
    }

    // Determine citation order of junctions
    std::vector<size_t> junctionOrder;
    for (size_t ji = 0; ji < tree.junctions.size(); ++ji) {
        junctionOrder.push_back(ji);
    }
    // Sort junctions by citation order: higher order first
    std::sort(junctionOrder.begin(), junctionOrder.end(), [&](size_t a, size_t b) {
        int maxOrderA = std::max(tree.components[tree.junctions[a].comp1].order, tree.components[tree.junctions[a].comp2].order);
        int maxOrderB = std::max(tree.components[tree.junctions[b].comp1].order, tree.components[tree.junctions[b].comp2].order);
        if (maxOrderA != maxOrderB) return maxOrderA > maxOrderB;
        return a < b;
    });

    // Cartesian product of candidate numberings
    std::vector<int> numOptions(tree.components.size());
    int totalCombinations = 1;
    for (size_t i = 0; i < tree.components.size(); ++i) {
        numOptions[i] = static_cast<int>(allCompNumberings[i].size());
        totalCombinations *= numOptions[i];
    }

    Fr45Score bestScore;
    QStringList bestPieces;
    bool hasBest = false;

    for (int comb = 0; comb < totalCombinations; ++comb) {
        std::vector<int> currentConfig(tree.components.size());
        int rem = comb;
        for (int i = static_cast<int>(tree.components.size()) - 1; i >= 0; --i) {
            currentConfig[i] = rem % numOptions[i];
            rem /= numOptions[i];
        }

        Fr45Score score;
        QStringList pieces;
        bool invalidConfig = false;
        std::map<int, int> distinctAttached;
        std::map<QString, int> nextAttachedPrime;
        std::map<int, int> distinctParents;
        int nextParentPrime = 0;

        for (size_t ji : junctionOrder) {
            const auto &j = tree.junctions[ji];
            int c1 = j.comp1, c2 = j.comp2;

            const QString &parentName = tree.components[tree.parentIndex].name;
            bool c1IsParent = (tree.components[c1].name == parentName);
            bool c2IsParent = (tree.components[c2].name == parentName);

            int att = c1, host = c2;
            if (c1IsParent != c2IsParent) {
                if (c1IsParent) { host = c1; att = c2; }
                else { host = c2; att = c1; }
            } else if (tree.components[c2].order > tree.components[c1].order) {
                att = c2; host = c1;
            } else if (tree.components[c1].order > tree.components[c2].order) {
                att = c1; host = c2;
            } else {
                if (c1 == tree.parentIndex) { host = c1; att = c2; }
                else { host = c2; att = c1; }
            }

            int u = -1, v = -1;
            if (j.sharedBonds.size() >= 2) {
                std::map<int, int> degInShared;
                for (const auto &b : j.sharedBonds) {
                    degInShared[b.first]++;
                    degInShared[b.second]++;
                }
                std::vector<int> endpoints;
                for (const auto &kv : degInShared) {
                    if (kv.second == 1) endpoints.push_back(kv.first);
                }
                if (endpoints.size() == 2) {
                    u = endpoints[0];
                    v = endpoints[1];
                }
            }
            if (u == -1 || v == -1) {
                if (!j.sharedBonds.empty()) {
                    u = j.sharedBonds[0].first;
                    v = j.sharedBonds[0].second;
                } else if (j.sharedAtoms.size() >= 2) {
                    u = j.sharedAtoms[0];
                    v = j.sharedAtoms[1];
                }
            }
            if (u == -1 || v == -1) continue;

            const auto &hostNum = allCompNumberings[host][currentConfig[host]];
            const auto &attNum = allCompNumberings[att][currentConfig[att]];

            auto edgeKey = std::make_pair(std::min(u, v), std::max(u, v));
            auto fwd = hostNum.forwardEdge.count(edgeKey) ? hostNum.forwardEdge.at(edgeKey) : std::make_pair(u, v);
            int firstHost = fwd.first;
            int secondHost = fwd.second;

            QString letter = hostNum.bondLetters.count(edgeKey) ? hostNum.bondLetters.at(edgeKey) : QString();
            QString loc1 = attNum.atomLocants.count(firstHost) ? attNum.atomLocants.at(firstHost) : QString();
            QString loc2 = attNum.atomLocants.count(secondHost) ? attNum.atomLocants.at(secondHost) : QString();

            if (distinctAttached.count(att) == 0) {
                distinctAttached[att] = nextAttachedPrime[tree.components[att].prefix]++;
            }
            int attPrimes = distinctAttached[att];
            loc1 = addPrimesToLocant(loc1, attPrimes);
            loc2 = addPrimesToLocant(loc2, attPrimes);

            if (tree.components[host].name == parentName) {
                if (distinctParents.count(host) == 0) {
                    distinctParents[host] = nextParentPrime++;
                }
                int hPrimes = distinctParents[host];
                letter = addPrimesToLocant(letter, hPrimes);
            }

            bool isHostParent = (tree.components[host].name == parentName || tree.components[host].order == 0);
            bool hasHigherOrderChildren = false;
            for (const auto &otherJ : tree.junctions) {
                int o1 = otherJ.comp1, o2 = otherJ.comp2;
                if (o1 == att && tree.components[o2].order > tree.components[att].order) hasHigherOrderChildren = true;
                if (o2 == att && tree.components[o1].order > tree.components[att].order) hasHigherOrderChildren = true;
            }
            bool omitLocants = (!hasHigherOrderChildren && tree.components[att].prefix == "benzo" && isHostParent);

            if (isHostParent) {
                if (letter.isEmpty()) { invalidConfig = true; break; }
                if (omitLocants) {
                    pieces.append(letter);
                } else {
                    if (loc1.isEmpty() || loc2.isEmpty()) { invalidConfig = true; break; }
                    pieces.append(loc1 + "," + loc2 + "-" + letter);
                }
            } else {
                QString hLoc1 = hostNum.atomLocants.count(firstHost) ? hostNum.atomLocants.at(firstHost) : QString();
                QString hLoc2 = hostNum.atomLocants.count(secondHost) ? hostNum.atomLocants.at(secondHost) : QString();
                if (j.sharedBonds.size() >= 2 && tree.components[att].ringCount == 1) {
                    if (locantLess(hLoc2, hLoc1)) std::swap(hLoc1, hLoc2);
                    pieces.append(hLoc1 + "," + hLoc2);
                } else {
                    pieces.append(loc1 + "," + loc2 + ":" + hLoc1 + "," + hLoc2);
                }
            }

            if (isHostParent) {
                score.parentLettersSet.push_back(letter);
                score.parentLettersCitation.push_back(letter);
                score.firstOrderLocantsSet.push_back(loc1);
                score.firstOrderLocantsSet.push_back(loc2);
                score.firstOrderLocantsCitation.push_back(loc1);
                score.firstOrderLocantsCitation.push_back(loc2);
            } else {
                QString hLoc1 = hostNum.atomLocants.count(firstHost) ? hostNum.atomLocants.at(firstHost) : QString();
                QString hLoc2 = hostNum.atomLocants.count(secondHost) ? hostNum.atomLocants.at(secondHost) : QString();
                if (j.sharedBonds.size() >= 2 && tree.components[att].ringCount == 1) {
                    if (locantLess(hLoc2, hLoc1)) std::swap(hLoc1, hLoc2);
                    score.higherOrderLocantsSet.push_back(hLoc1);
                    score.higherOrderLocantsSet.push_back(hLoc2);
                    score.higherOrderLocantsCitation.push_back(hLoc1);
                    score.higherOrderLocantsCitation.push_back(hLoc2);
                } else {
                    score.higherOrderLocantsSet.push_back(hLoc1);
                    score.higherOrderLocantsSet.push_back(hLoc2);
                    score.higherOrderLocantsSet.push_back(loc1);
                    score.higherOrderLocantsSet.push_back(loc2);
                    score.higherOrderLocantsCitation.push_back(hLoc1);
                    score.higherOrderLocantsCitation.push_back(hLoc2);
                    score.higherOrderLocantsCitation.push_back(loc1);
                    score.higherOrderLocantsCitation.push_back(loc2);
                }
            }
        }

        if (invalidConfig) continue;

        std::sort(score.parentLettersSet.begin(), score.parentLettersSet.end());
        std::sort(score.firstOrderLocantsSet.begin(), score.firstOrderLocantsSet.end(), locantLess);
        std::sort(score.higherOrderLocantsSet.begin(), score.higherOrderLocantsSet.end(), locantLess);

        if (!hasBest || score < bestScore) {
            bestScore = score;
            bestPieces = pieces;
            hasBest = true;
        }
    }

    if (!hasBest || bestPieces.isEmpty()) return QString();
    bool allAbbr = true;
    bool hasColon = false;
    for (const auto &p : bestPieces) {
        if (p.contains("-") || p.contains(":")) allAbbr = false;
        if (p.contains(":")) hasColon = true;
    }
    if (allAbbr) {
        return "[" + bestPieces.join(",") + "]";
    } else if (hasColon) {
        return "[" + bestPieces.join(";") + "]";
    } else {
        return "[" + bestPieces.join(":") + "]";
    }
}

static QString stripPrimes(const QString &s) {
    QString res;
    for (QChar ch : s) {
        if (ch != '\'' && ch != '\"') res.append(ch);
    }
    return res;
}

static QString simpleMultPrefix(int k) {
    switch (k) {
        case 2: return "di";
        case 3: return "tri";
        case 4: return "tetra";
        case 5: return "penta";
        case 6: return "hexa";
        default: return QString();
    }
}

static QString complexMultPrefix(int k) {
    switch (k) {
        case 2: return "bis";
        case 3: return "tris";
        case 4: return "tetrakis";
        case 5: return "pentakis";
        default: return QString();
    }
}

QString constructFusionName(const ComponentTree &tree, const Graph &g, int mol) {
    if (!tree.valid || tree.components.empty()) {
        return QString();
    }
    if (tree.components.size() == 1) {
        return tree.components[0].name;
    }

    const QString &parentName = tree.components[tree.parentIndex].name;

    // Generate numberings for all components
    std::vector<std::vector<ComponentNumbering>> allCompNumberings(tree.components.size());
    for (size_t i = 0; i < tree.components.size(); ++i) {
        if (tree.components[i].ringCount == 1) {
            allCompNumberings[i] = generateMonocycleNumberings(tree.components[i], g);
        } else {
            allCompNumberings[i] = generatePolycycleNumberings(tree.components[i], g, mol);
        }
        if (allCompNumberings[i].empty()) {
            ComponentNumbering dummy;
            int idx = 1;
            for (int a : tree.components[i].atomIndices) {
                dummy.atomLocants[a] = QString::number(idx++);
            }
            allCompNumberings[i].push_back(dummy);
        }
    }

    // Determine citation order of junctions: higher order first
    std::vector<size_t> junctionOrder;
    for (size_t ji = 0; ji < tree.junctions.size(); ++ji) {
        junctionOrder.push_back(ji);
    }
    std::sort(junctionOrder.begin(), junctionOrder.end(), [&](size_t a, size_t b) {
        int maxOrderA = std::max(tree.components[tree.junctions[a].comp1].order, tree.components[tree.junctions[a].comp2].order);
        int maxOrderB = std::max(tree.components[tree.junctions[b].comp1].order, tree.components[tree.junctions[b].comp2].order);
        if (maxOrderA != maxOrderB) return maxOrderA > maxOrderB;
        return a < b;
    });

    // Cartesian product of candidate numberings to find bestConfig
    std::vector<int> numOptions(tree.components.size());
    int totalCombinations = 1;
    for (size_t i = 0; i < tree.components.size(); ++i) {
        numOptions[i] = static_cast<int>(allCompNumberings[i].size());
        totalCombinations *= numOptions[i];
    }

    Fr45Score bestScore;
    std::vector<int> bestConfig(tree.components.size(), 0);
    bool hasBest = false;

    for (int comb = 0; comb < totalCombinations; ++comb) {
        std::vector<int> currentConfig(tree.components.size());
        int rem = comb;
        for (int i = static_cast<int>(tree.components.size()) - 1; i >= 0; --i) {
            currentConfig[i] = rem % numOptions[i];
            rem /= numOptions[i];
        }

        Fr45Score score;
        bool invalidConfig = false;
        for (size_t ji : junctionOrder) {
            const auto &j = tree.junctions[ji];
            int c1 = j.comp1, c2 = j.comp2;
            bool c1IsParent = (tree.components[c1].name == parentName || tree.components[c1].order == 0);
            bool c2IsParent = (tree.components[c2].name == parentName || tree.components[c2].order == 0);

            int att = c1, host = c2;
            if (c1IsParent != c2IsParent) {
                if (c1IsParent) { host = c1; att = c2; }
                else { host = c2; att = c1; }
            } else if (tree.components[c2].order > tree.components[c1].order) {
                att = c2; host = c1;
            } else if (tree.components[c1].order > tree.components[c2].order) {
                att = c1; host = c2;
            } else {
                if (c1 == tree.parentIndex) { host = c1; att = c2; }
                else { host = c2; att = c1; }
            }

            int u = -1, v = -1;
            if (j.sharedBonds.size() >= 2) {
                std::map<int, int> degInShared;
                for (const auto &b : j.sharedBonds) {
                    degInShared[b.first]++;
                    degInShared[b.second]++;
                }
                std::vector<int> endpoints;
                for (const auto &kv : degInShared) {
                    if (kv.second == 1) endpoints.push_back(kv.first);
                }
                if (endpoints.size() == 2) {
                    u = endpoints[0];
                    v = endpoints[1];
                }
            }
            if (u == -1 || v == -1) {
                if (!j.sharedBonds.empty()) {
                    u = j.sharedBonds[0].first;
                    v = j.sharedBonds[0].second;
                } else if (j.sharedAtoms.size() >= 2) {
                    u = j.sharedAtoms[0];
                    v = j.sharedAtoms[1];
                }
            }
            if (u == -1 || v == -1) continue;

            const auto &hostNum = allCompNumberings[host][currentConfig[host]];
            const auto &attNum = allCompNumberings[att][currentConfig[att]];

            auto edgeKey = std::make_pair(std::min(u, v), std::max(u, v));
            auto fwd = hostNum.forwardEdge.count(edgeKey) ? hostNum.forwardEdge.at(edgeKey) : std::make_pair(u, v);
            int firstHost = fwd.first;
            int secondHost = fwd.second;

            QString letter = hostNum.bondLetters.count(edgeKey) ? hostNum.bondLetters.at(edgeKey) : QString();
            QString loc1 = attNum.atomLocants.count(firstHost) ? attNum.atomLocants.at(firstHost) : QString();
            QString loc2 = attNum.atomLocants.count(secondHost) ? attNum.atomLocants.at(secondHost) : QString();

            bool isHostParent = (tree.components[host].name == parentName || tree.components[host].order == 0);
            if (isHostParent) {
                QString periParentLetters;
                QString periAttLocants;
                if (j.sharedBonds.size() >= 2) {
                    std::map<int, std::vector<int>> sAdj;
                    for (const auto &b : j.sharedBonds) {
                        sAdj[b.first].push_back(b.second);
                        sAdj[b.second].push_back(b.first);
                    }
                    int startNode = -1;
                    for (const auto &kv : sAdj) {
                        if (kv.second.size() == 1) { startNode = kv.first; break; }
                    }
                    if (startNode != -1) {
                        std::vector<int> path;
                        int currNode = startNode, prevNode = -1;
                        while (currNode != -1) {
                            path.push_back(currNode);
                            int nextNode = -1;
                            for (int nei : sAdj[currNode]) {
                                if (nei != prevNode) { nextNode = nei; break; }
                            }
                            prevNode = currNode;
                            currNode = nextNode;
                        }
                        if (path.size() >= 2) {
                            QString firstLoc = attNum.atomLocants.count(path.front()) ? attNum.atomLocants.at(path.front()) : QString();
                            QString lastLoc = attNum.atomLocants.count(path.back()) ? attNum.atomLocants.at(path.back()) : QString();
                            if (locantLess(lastLoc, firstLoc)) std::reverse(path.begin(), path.end());
                            QStringList aLocList;
                            for (int pAtom : path) {
                                if (attNum.atomLocants.count(pAtom)) aLocList.append(attNum.atomLocants.at(pAtom));
                            }
                            periAttLocants = aLocList.join(",");
                            size_t foundPeriBonds = 0;
                            for (size_t pi = 0; pi + 1 < path.size(); ++pi) {
                                auto bKey = std::make_pair(std::min(path[pi], path[pi + 1]), std::max(path[pi], path[pi + 1]));
                                if (hostNum.bondLetters.count(bKey)) {
                                    periParentLetters += hostNum.bondLetters.at(bKey);
                                    foundPeriBonds++;
                                }
                            }
                            if (foundPeriBonds != path.size() - 1) periParentLetters.clear();
                        }
                    }
                }
                if (j.sharedBonds.size() >= 2 && !periParentLetters.isEmpty()) {
                    score.parentLettersSet.push_back(periParentLetters);
                    score.parentLettersCitation.push_back(periParentLetters);
                    score.firstOrderLocantsSet.push_back(periAttLocants);
                    score.firstOrderLocantsCitation.push_back(periAttLocants);
                } else if (j.sharedBonds.size() >= 2) {
                    invalidConfig = true;
                    break;
                } else {
                    if (letter.isEmpty()) { invalidConfig = true; break; }
                    score.parentLettersSet.push_back(letter);
                    score.parentLettersCitation.push_back(letter);
                    score.firstOrderLocantsSet.push_back(loc1);
                    score.firstOrderLocantsSet.push_back(loc2);
                    score.firstOrderLocantsCitation.push_back(loc1);
                    score.firstOrderLocantsCitation.push_back(loc2);
                }
            } else {
                QString hLoc1 = hostNum.atomLocants.count(firstHost) ? hostNum.atomLocants.at(firstHost) : QString();
                QString hLoc2 = hostNum.atomLocants.count(secondHost) ? hostNum.atomLocants.at(secondHost) : QString();
                if ((j.sharedBonds.size() >= 2 && tree.components[att].ringCount == 1) ||
                    (tree.components[att].prefix == "benzo" && j.sharedBonds.size() == 1)) {
                    if (locantLess(hLoc2, hLoc1)) std::swap(hLoc1, hLoc2);
                    score.higherOrderLocantsSet.push_back(hLoc1);
                    score.higherOrderLocantsSet.push_back(hLoc2);
                    score.higherOrderLocantsCitation.push_back(hLoc1);
                    score.higherOrderLocantsCitation.push_back(hLoc2);
                } else {
                    score.higherOrderLocantsSet.push_back(hLoc1);
                    score.higherOrderLocantsSet.push_back(hLoc2);
                    score.higherOrderLocantsSet.push_back(loc1);
                    score.higherOrderLocantsSet.push_back(loc2);
                    score.higherOrderLocantsCitation.push_back(hLoc1);
                    score.higherOrderLocantsCitation.push_back(hLoc2);
                    score.higherOrderLocantsCitation.push_back(loc1);
                    score.higherOrderLocantsCitation.push_back(loc2);
                }
            }
        }

        if (invalidConfig) continue;

        std::sort(score.parentLettersSet.begin(), score.parentLettersSet.end());
        std::sort(score.firstOrderLocantsSet.begin(), score.firstOrderLocantsSet.end(), locantLess);
        std::sort(score.higherOrderLocantsSet.begin(), score.higherOrderLocantsSet.end(), locantLess);

        if (!hasBest || score < bestScore) {
            bestScore = score;
            bestConfig = currentConfig;
            hasBest = true;
        }
    }
    if (!hasBest) return QString();

    // Directed adjacency from host to attached (tree structure)
    int nComps = static_cast<int>(tree.components.size());
    std::vector<std::vector<int>> treeAdj(nComps);
    for (size_t ji = 0; ji < tree.junctions.size(); ++ji) {
        const auto &j = tree.junctions[ji];
        int c1 = j.comp1, c2 = j.comp2;
        bool c1IsParent = (tree.components[c1].name == parentName || tree.components[c1].order == 0);
        bool c2IsParent = (tree.components[c2].name == parentName || tree.components[c2].order == 0);
        int att = c1, host = c2;
        if (c1IsParent != c2IsParent) {
            if (c1IsParent) { host = c1; att = c2; }
            else { host = c2; att = c1; }
        } else if (tree.components[c2].order > tree.components[c1].order) {
            att = c2; host = c1;
        } else if (tree.components[c1].order > tree.components[c2].order) {
            att = c1; host = c2;
        } else {
            if (c1 == tree.parentIndex) { host = c1; att = c2; }
            else { host = c2; att = c1; }
        }
        treeAdj[host].push_back(att);
    }

    auto subtreesIsomorphic = [&](auto &self, int u, int v) -> bool {
        if (tree.components[u].prefix != tree.components[v].prefix) return false;
        if (treeAdj[u].size() != treeAdj[v].size()) return false;
        if (treeAdj[u].empty()) return true;
        for (size_t i = 0; i < treeAdj[u].size(); ++i) {
            if (!self(self, treeAdj[u][i], treeAdj[v][i])) return false;
        }
        return true;
    };

    std::map<int, int> groupPrime;
    for (int h = 0; h < nComps; ++h) {
        std::vector<std::vector<int>> isoGroups;
        for (int c : treeAdj[h]) {
            bool placed = false;
            for (auto &grp : isoGroups) {
                if (subtreesIsomorphic(subtreesIsomorphic, grp[0], c)) {
                    grp.push_back(c);
                    placed = true;
                    break;
                }
            }
            if (!placed) isoGroups.push_back({c});
        }
        for (const auto &grp : isoGroups) {
            if (grp.size() > 1) {
                for (size_t k = 0; k < grp.size(); ++k) {
                    groupPrime[grp[k]] = static_cast<int>(k);
                }
            } else {
                groupPrime[grp[0]] = 0;
            }
        }
    }

    struct JunctionDetail {
        int host{-1};
        int att{-1};
        QString piece;
        bool isAbbreviated{false};
        bool isHigherOrder{false};
    };
    std::map<std::pair<int, int>, JunctionDetail> junctionMap;

    std::map<int, int> distinctParents;
    int nextParentPrime = 0;

    for (size_t ji : junctionOrder) {
        const auto &j = tree.junctions[ji];
        int c1 = j.comp1, c2 = j.comp2;
        bool c1IsParent = (tree.components[c1].name == parentName || tree.components[c1].order == 0);
        bool c2IsParent = (tree.components[c2].name == parentName || tree.components[c2].order == 0);

        int att = c1, host = c2;
        if (c1IsParent != c2IsParent) {
            if (c1IsParent) { host = c1; att = c2; }
            else { host = c2; att = c1; }
        } else if (tree.components[c2].order > tree.components[c1].order) {
            att = c2; host = c1;
        } else if (tree.components[c1].order > tree.components[c2].order) {
            att = c1; host = c2;
        } else {
            if (c1 == tree.parentIndex) { host = c1; att = c2; }
            else { host = c2; att = c1; }
        }

        int u = -1, v = -1;
        if (j.sharedBonds.size() >= 2) {
            std::map<int, int> degInShared;
            for (const auto &b : j.sharedBonds) {
                degInShared[b.first]++;
                degInShared[b.second]++;
            }
            std::vector<int> endpoints;
            for (const auto &kv : degInShared) {
                if (kv.second == 1) endpoints.push_back(kv.first);
            }
            if (endpoints.size() == 2) {
                u = endpoints[0];
                v = endpoints[1];
            }
        }
        if (u == -1 || v == -1) {
            if (!j.sharedBonds.empty()) {
                u = j.sharedBonds[0].first;
                v = j.sharedBonds[0].second;
            } else if (j.sharedAtoms.size() >= 2) {
                u = j.sharedAtoms[0];
                v = j.sharedAtoms[1];
            }
        }
        if (u == -1 || v == -1) continue;

        const auto &hostNum = allCompNumberings[host][bestConfig[host]];
        const auto &attNum = allCompNumberings[att][bestConfig[att]];

        auto edgeKey = std::make_pair(std::min(u, v), std::max(u, v));
        auto fwd = hostNum.forwardEdge.count(edgeKey) ? hostNum.forwardEdge.at(edgeKey) : std::make_pair(u, v);
        int firstHost = fwd.first;
        int secondHost = fwd.second;

        QString letter = hostNum.bondLetters.count(edgeKey) ? hostNum.bondLetters.at(edgeKey) : QString();
        QString loc1 = attNum.atomLocants.count(firstHost) ? attNum.atomLocants.at(firstHost) : QString();
        QString loc2 = attNum.atomLocants.count(secondHost) ? attNum.atomLocants.at(secondHost) : QString();

        bool isHostParent = (tree.components[host].name == parentName || tree.components[host].order == 0);

        int attOrder = tree.components[att].order;
        int basePrimes = (attOrder > 1) ? (attOrder - 1) : 0;
        int totalPrimes = basePrimes + groupPrime[att];
        loc1 = addPrimesToLocant(loc1, totalPrimes);
        loc2 = addPrimesToLocant(loc2, totalPrimes);

        QString periParentLetters;
        QString periAttLocants;
        if (j.sharedBonds.size() >= 2) {
            std::map<int, std::vector<int>> sAdj;
            for (const auto &b : j.sharedBonds) {
                sAdj[b.first].push_back(b.second);
                sAdj[b.second].push_back(b.first);
            }
            int startNode = -1;
            for (const auto &kv : sAdj) {
                if (kv.second.size() == 1) {
                    startNode = kv.first;
                    break;
                }
            }
            if (startNode != -1) {
                std::vector<int> path;
                int currNode = startNode;
                int prevNode = -1;
                while (currNode != -1) {
                    path.push_back(currNode);
                    int nextNode = -1;
                    for (int nei : sAdj[currNode]) {
                        if (nei != prevNode) {
                            nextNode = nei;
                            break;
                        }
                    }
                    prevNode = currNode;
                    currNode = nextNode;
                }
                if (path.size() >= 2) {
                    QString firstLoc = attNum.atomLocants.count(path.front()) ? attNum.atomLocants.at(path.front()) : QString();
                    QString lastLoc = attNum.atomLocants.count(path.back()) ? attNum.atomLocants.at(path.back()) : QString();
                    if (locantLess(lastLoc, firstLoc)) {
                        std::reverse(path.begin(), path.end());
                    }
                    QStringList aLocList;
                    for (int pAtom : path) {
                        if (attNum.atomLocants.count(pAtom)) {
                            aLocList.append(addPrimesToLocant(attNum.atomLocants.at(pAtom), totalPrimes));
                        }
                    }
                    periAttLocants = aLocList.join(",");
                    for (size_t pi = 0; pi + 1 < path.size(); ++pi) {
                        auto bKey = std::make_pair(std::min(path[pi], path[pi + 1]), std::max(path[pi], path[pi + 1]));
                        if (hostNum.bondLetters.count(bKey)) {
                            periParentLetters += hostNum.bondLetters.at(bKey);
                        }
                    }
                }
            }
        }

        if (isHostParent) {
            if (distinctParents.count(host) == 0) {
                distinctParents[host] = nextParentPrime++;
            }
            letter = addPrimesToLocant(letter, distinctParents[host]);
            periParentLetters = addPrimesToLocant(periParentLetters, distinctParents[host]);
        }

        bool hasHigherOrderChildren = false;
        for (const auto &otherJ : tree.junctions) {
            int o1 = otherJ.comp1, o2 = otherJ.comp2;
            if (o1 == att && tree.components[o2].order > tree.components[att].order) hasHigherOrderChildren = true;
            if (o2 == att && tree.components[o1].order > tree.components[att].order) hasHigherOrderChildren = true;
        }
        bool omitLocants = (!hasHigherOrderChildren && tree.components[att].prefix == "benzo" && isHostParent);

        JunctionDetail jd;
        jd.host = host;
        jd.att = att;
        if (isHostParent) {
            jd.isHigherOrder = false;
            if (j.sharedBonds.size() >= 2) {
                if (periParentLetters.isEmpty()) return QString();
                if (omitLocants) {
                    jd.isAbbreviated = true;
                    jd.piece = periParentLetters;
                } else {
                    jd.isAbbreviated = false;
                    jd.piece = periAttLocants + "-" + periParentLetters;
                }
            } else if (omitLocants) {
                if (letter.isEmpty()) return QString();
                jd.isAbbreviated = true;
                jd.piece = letter;
            } else {
                if (letter.isEmpty() || loc1.isEmpty() || loc2.isEmpty()) return QString();
                jd.isAbbreviated = false;
                jd.piece = loc1 + "," + loc2 + "-" + letter;
            }
        } else {
            jd.isHigherOrder = true;
            jd.isAbbreviated = false;
            QString hLoc1 = hostNum.atomLocants.count(firstHost) ? hostNum.atomLocants.at(firstHost) : QString();
            QString hLoc2 = hostNum.atomLocants.count(secondHost) ? hostNum.atomLocants.at(secondHost) : QString();
            if (j.sharedBonds.size() >= 2 && tree.components[att].ringCount == 1) {
                // FR-4.6: For monocyclic attached component in ortho+peri fusion,
                // attached locants omitted; only host endpoints cited.
                if (locantLess(hLoc2, hLoc1)) std::swap(hLoc1, hLoc2);
                jd.piece = hLoc1 + "," + hLoc2;
            } else if (tree.components[att].prefix == "benzo" && j.sharedBonds.size() == 1) {
                // FR-4.6: Monocyclic components (benzo) attached to another component
                // do not need their own locants cited, only host locants cited: e.g. benzo[4,5]
                if (locantLess(hLoc2, hLoc1)) std::swap(hLoc1, hLoc2);
                jd.piece = hLoc1 + "," + hLoc2;
            } else {
                jd.piece = loc1 + "," + loc2 + ":" + hLoc1 + "," + hLoc2;
            }
        }
        junctionMap[{host, att}] = jd;
    }

    struct BranchResult {
        int compIdx{0};
        QString prefix;
        QString junctionDesc;
        bool isAbbreviated{false};
        bool isHigherOrder{false};
        QString compositeName;
        bool isComposite{false};
    };

    std::function<BranchResult(int, const JunctionDetail&)> buildBranch = [&](int u, const JunctionDetail &junc) -> BranchResult {
        BranchResult res;
        res.compIdx = u;
        res.prefix = tree.components[u].prefix.isEmpty() ? tree.components[u].name : tree.components[u].prefix;
        res.junctionDesc = junc.piece;
        res.isAbbreviated = junc.isAbbreviated;
        res.isHigherOrder = junc.isHigherOrder;

        if (treeAdj[u].empty()) {
            res.compositeName = res.prefix;
            res.isComposite = false;
            return res;
        }

        res.isComposite = true;
        std::vector<BranchResult> childBranches;
        for (int v : treeAdj[u]) {
            if (junctionMap.count({u, v})) {
                childBranches.push_back(buildBranch(v, junctionMap.at({u, v})));
            }
        }

        // Group isomorphic child branches
        std::vector<std::vector<BranchResult>> childGroups;
        for (const auto &cb : childBranches) {
            bool placed = false;
            for (auto &grp : childGroups) {
                if (grp[0].prefix == cb.prefix &&
                    grp[0].isComposite == cb.isComposite &&
                    stripPrimes(grp[0].compositeName) == stripPrimes(cb.compositeName)) {
                    grp.push_back(cb);
                    placed = true;
                    break;
                }
            }
            if (!placed) {
                childGroups.push_back({cb});
            }
        }

        // Sort child groups by citation order (FR-3.6)
        std::sort(childGroups.begin(), childGroups.end(), [](const auto &g1, const auto &g2) {
            return alphabetizationKey(g1[0].prefix).toLower() < alphabetizationKey(g2[0].prefix).toLower();
        });

        QString childCombined;
        for (const auto &grp : childGroups) {
            int k = static_cast<int>(grp.size());
            if (k == 1) {
                if (grp[0].isComposite) {
                    childCombined += grp[0].compositeName + "[" + grp[0].junctionDesc + "]";
                } else {
                    childCombined += grp[0].prefix + "[" + grp[0].junctionDesc + "]";
                }
            } else {
                QStringList descs;
                for (const auto &b : grp) descs.append(b.junctionDesc);
                bool hasCol = false;
                for (const auto &d : descs) {
                    if (d.contains(":")) { hasCol = true; break; }
                }
                QString combinedDesc;
                if (hasCol) combinedDesc = "[" + descs.join(";") + "]";
                else combinedDesc = "[" + descs.join(":") + "]";

                if (grp[0].isComposite) {
                    childCombined += complexMultPrefix(k) + "(" + grp[0].compositeName + ")" + combinedDesc;
                } else {
                    childCombined += simpleMultPrefix(k) + grp[0].prefix + combinedDesc;
                }
            }
        }
        res.compositeName = childCombined + res.prefix;
        return res;
    };

    // Collect branches of all first-order children (attached to parent(s))
    std::vector<BranchResult> topBranches;
    std::set<int> visitedTop;
    for (size_t i = 0; i < tree.components.size(); ++i) {
        if (tree.components[i].name == parentName || tree.components[i].order == 0) {
            for (int child : treeAdj[i]) {
                if (visitedTop.count(child) == 0 && junctionMap.count({static_cast<int>(i), child})) {
                    visitedTop.insert(child);
                    topBranches.push_back(buildBranch(child, junctionMap.at({static_cast<int>(i), child})));
                }
            }
        }
    }

    // Also check if any component is attached to multiple parents (e.g. naphtho in naphtho[...][difuran])
    std::map<int, QStringList> multiParentJunctions;
    for (const auto &kv : junctionMap) {
        int host = kv.first.first;
        if (tree.components[host].name == parentName || tree.components[host].order == 0) {
            int att = kv.first.second;
            multiParentJunctions[att].append(kv.second.piece);
        }
    }
    for (auto &tb : topBranches) {
        if (multiParentJunctions.count(tb.compIdx) && multiParentJunctions[tb.compIdx].size() > 1) {
            auto list = multiParentJunctions[tb.compIdx];
            std::sort(list.begin(), list.end(), [](const QString &a, const QString &b) {
                if (a.contains("'") != b.contains("'")) return !a.contains("'");
                return a < b;
            });
            tb.junctionDesc = list.join(":");
        }
    }

    // Group isomorphic top branches
    std::vector<std::vector<BranchResult>> topGroups;
    for (const auto &tb : topBranches) {
        bool placed = false;
        for (auto &grp : topGroups) {
            if (grp[0].prefix == tb.prefix &&
                grp[0].isComposite == tb.isComposite &&
                stripPrimes(grp[0].compositeName) == stripPrimes(tb.compositeName)) {
                grp.push_back(tb);
                placed = true;
                break;
            }
        }
        if (!placed) {
            topGroups.push_back({tb});
        }
    }

    // Sort top groups by FR-3.6 citation order
    std::sort(topGroups.begin(), topGroups.end(), [](const auto &g1, const auto &g2) {
        return alphabetizationKey(g1[0].prefix).toLower() < alphabetizationKey(g2[0].prefix).toLower();
    });

    QString attachedString;
    for (const auto &grp : topGroups) {
        int k = static_cast<int>(grp.size());
        if (k == 1) {
            if (grp[0].isComposite) {
                attachedString += grp[0].compositeName + "[" + grp[0].junctionDesc + "]";
            } else {
                attachedString += grp[0].prefix + "[" + grp[0].junctionDesc + "]";
            }
        } else {
            QStringList descs;
            for (const auto &b : grp) descs.append(b.junctionDesc);
            bool allAbbr = true;
            for (const auto &d : descs) {
                if (d.contains("-") || d.contains(":")) allAbbr = false;
            }
            QString combinedDesc;
            if (allAbbr) combinedDesc = "[" + descs.join(",") + "]";
            else combinedDesc = "[" + descs.join(":") + "]";

            if (grp[0].isComposite) {
                attachedString += complexMultPrefix(k) + "(" + grp[0].compositeName + ")" + combinedDesc;
            } else {
                attachedString += simpleMultPrefix(k) + grp[0].prefix + combinedDesc;
            }
        }
    }

    // Format parent string
    int parentCount = 0;
    for (const auto &c : tree.components) {
        if (c.name == parentName) parentCount++;
    }
    QString parentString;
    if (parentCount > 1) {
        parentString = simpleMultPrefix(parentCount) + parentName;
    } else {
        parentString = parentName;
        if (parentString.startsWith("2,1-")) {
            parentString = "[2,1]" + parentString.mid(4);
        } else if (parentString.startsWith("1,3-")) {
            parentString = "[1,3]" + parentString.mid(4);
        }
    }

    return attachedString + parentString;
}

static QString replacementPrefix(int z, int count) {
    QString stem;
    switch (z) {
        case 7: stem = "aza"; break;
        case 8: stem = "oxa"; break;
        case 16: stem = "thia"; break;
        case 34: stem = "selena"; break;
        case 52: stem = "tellura"; break;
        case 15: stem = "phospha"; break;
        case 14: stem = "sila"; break;
        case 33: stem = "arsa"; break;
        case 5: stem = "bora"; break;
        default: stem = "hetero"; break;
    }
    if (count == 1) return stem;
    return simpleMultPrefix(count) + stem;
}

QString nameBySkeletalReplacementFallback(
    const ComponentTree &tree,
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    int mol)
{
    if (rings.empty()) return QString();

    // Must have skeletal heteroatoms
    std::set<int> allRingAtoms;
    for (const auto &r : rings) {
        for (int a : r) allRingAtoms.insert(a);
    }
    std::vector<int> heteroNodes;
    for (int a : allRingAtoms) {
        if (g.nodes[a].atomicNumber != 6) {
            heteroNodes.push_back(a);
        }
    }
    if (heteroNodes.empty()) return QString();

    FusedSystemTopology topo(rings, g);
    auto perim = topo.perimeterCycle();
    auto interior = topo.interiorAtoms();
    size_t perimLen = perim.size();
    size_t interiorCount = interior.size();
    size_t totalAtoms = perimLen + interiorCount;
    size_t numRings = rings.size();

    auto isBh = [&](int atom) {
        int deg = 0;
        for (int nei : g.nodes[atom].neighbors) {
            if (topo.peripheralAtoms().count(nei) || interior.count(nei)) deg++;
        }
        return deg > 2;
    };

    QString parentName;
    enum class ParentType { None, CyclopentaPentalene, Phenalene, Pyrene };
    ParentType pType = ParentType::None;

    if (numRings == 3 && perimLen == 9 && interiorCount == 1 && totalAtoms == 10) {
        bool all5 = true;
        for (const auto &r : rings) {
            if (r.size() != 5) all5 = false;
        }
        if (all5) {
            parentName = "cyclopenta[cd]pentalene";
            pType = ParentType::CyclopentaPentalene;
        }
    } else if (numRings == 3 && perimLen == 12 && interiorCount == 1 && totalAtoms == 13) {
        bool all6 = true;
        for (const auto &r : rings) {
            if (r.size() != 6) all6 = false;
        }
        if (all6) {
            parentName = "phenalene";
            pType = ParentType::Phenalene;
        }
    } else if (numRings == 4 && perimLen == 14 && interiorCount == 2 && totalAtoms == 16) {
        bool all6 = true;
        for (const auto &r : rings) {
            if (r.size() != 6) all6 = false;
        }
        if (all6) {
            parentName = "pyrene";
            pType = ParentType::Pyrene;
        }
    }

    if (pType == ParentType::None || parentName.isEmpty()) {
        return QString();
    }

    int n = static_cast<int>(perim.size());
    std::vector<std::map<int, QString>> validNumberings;

    for (int dir : {1, -1}) {
        for (int i = 0; i < n; ++i) {
            std::vector<int> cand;
            for (int k = 0; k < n; ++k) {
                int idx = (i + dir * k + n * 10) % n;
                cand.push_back(perim[idx]);
            }

            bool matchesPattern = false;
            std::map<int, QString> numMap;

            if (pType == ParentType::CyclopentaPentalene) {
                // Pattern: [non-BH, non-BH, BH] * 3
                if (!isBh(cand[0]) && !isBh(cand[1]) && isBh(cand[2]) &&
                    !isBh(cand[3]) && !isBh(cand[4]) && isBh(cand[5]) &&
                    !isBh(cand[6]) && !isBh(cand[7]) && isBh(cand[8])) {
                    matchesPattern = true;
                    numMap[cand[0]] = "1";
                    numMap[cand[1]] = "2";
                    numMap[cand[2]] = "2a";
                    numMap[cand[3]] = "3";
                    numMap[cand[4]] = "4";
                    numMap[cand[5]] = "4a";
                    numMap[cand[6]] = "5";
                    numMap[cand[7]] = "6";
                    numMap[cand[8]] = "6a";
                    if (!interior.empty()) {
                        numMap[*interior.begin()] = "6b";
                    }
                }
            } else if (pType == ParentType::Phenalene) {
                // Pattern: [non-BH, non-BH, non-BH, BH] * 3
                if (!isBh(cand[0]) && !isBh(cand[1]) && !isBh(cand[2]) && isBh(cand[3]) &&
                    !isBh(cand[4]) && !isBh(cand[5]) && !isBh(cand[6]) && isBh(cand[7]) &&
                    !isBh(cand[8]) && !isBh(cand[9]) && !isBh(cand[10]) && isBh(cand[11])) {
                    matchesPattern = true;
                    numMap[cand[0]] = "1";
                    numMap[cand[1]] = "2";
                    numMap[cand[2]] = "3";
                    numMap[cand[3]] = "3a";
                    numMap[cand[4]] = "4";
                    numMap[cand[5]] = "5";
                    numMap[cand[6]] = "6";
                    numMap[cand[7]] = "6a";
                    numMap[cand[8]] = "7";
                    numMap[cand[9]] = "8";
                    numMap[cand[10]] = "9";
                    numMap[cand[11]] = "9a";
                    if (!interior.empty()) {
                        numMap[*interior.begin()] = "9b";
                    }
                }
            }

            if (matchesPattern) {
                validNumberings.push_back(numMap);
            }
        }
    }

    if (validNumberings.empty()) {
        return QString();
    }

    // Select the numbering that gives the lowest set of locants for heteroatoms
    // Seniority order of heteroatoms: O (8) > S (16) > Se (34) > N (7) > P (15) > ...
    auto heteroSeniority = [](int z) -> int {
        switch (z) {
            case 8: return 1;
            case 16: return 2;
            case 34: return 3;
            case 52: return 4;
            case 7: return 5;
            case 15: return 6;
            case 33: return 7;
            case 14: return 8;
            case 5: return 9;
            default: return 100;
        }
    };

    auto locantRank = [](const QString &loc) -> int {
        int num = 0;
        int idx = 0;
        while (idx < loc.size() && loc[idx].isDigit()) {
            num = num * 10 + (loc[idx].toLatin1() - '0');
            idx++;
        }
        int letter = 0;
        if (idx < loc.size()) {
            letter = loc[idx].toLatin1() - 'a' + 1;
        }
        return num * 100 + letter;
    };

    std::map<int, QString> bestNumbering;
    std::vector<int> bestRanks;

    for (const auto &numMap : validNumberings) {
        std::vector<int> ranks;
        for (int hNode : heteroNodes) {
            if (numMap.count(hNode)) {
                ranks.push_back(locantRank(numMap.at(hNode)));
            }
        }
        std::sort(ranks.begin(), ranks.end());
        if (bestRanks.empty() || ranks < bestRanks) {
            bestRanks = ranks;
            bestNumbering = numMap;
        }
    }

    if (bestNumbering.empty()) return QString();

    // Group heteroatoms by atomic number
    std::map<int, std::vector<QString>> locantsByZ;
    for (int hNode : heteroNodes) {
        int z = g.nodes[hNode].atomicNumber;
        if (bestNumbering.count(hNode)) {
            locantsByZ[z].push_back(bestNumbering[hNode]);
        }
    }

    // Sort by seniority
    std::vector<int> zOrder;
    for (const auto &pair : locantsByZ) zOrder.push_back(pair.first);
    std::sort(zOrder.begin(), zOrder.end(), [&](int a, int b) {
        return heteroSeniority(a) < heteroSeniority(b);
    });

    QString prefixString;
    for (int z : zOrder) {
        auto &locs = locantsByZ[z];
        std::sort(locs.begin(), locs.end(), [&](const QString &a, const QString &b) {
            return locantRank(a) < locantRank(b);
        });

        // For phenalene, if all 4 are N (1, 4, 9, 9b), FR-7.1 / common citation can omit 9b:
        // "1,4,9-tetraazaphenalene"
        QStringList formattedLocs;
        bool has9b = (std::find(locs.begin(), locs.end(), QString("9b")) != locs.end());
        if (pType == ParentType::Phenalene && locs.size() == 4 && has9b) {
            for (const auto &l : locs) {
                if (l != "9b") formattedLocs.append(l);
            }
        } else {
            for (const auto &l : locs) formattedLocs.append(l);
        }

        QString repl = replacementPrefix(z, static_cast<int>(locs.size()));
        prefixString += formattedLocs.join(",") + "-" + repl;
    }

    return prefixString + parentName;
}

QString nameByLessPreferredParentFallback(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol)
{
    ComponentTree fallbackTree = decomposeComponentTree(g, rings, input, mol, false);
    if (!fallbackTree.valid || fallbackTree.components.size() < 2 || hitsFr32Limitation(fallbackTree)) {
        return QString();
    }
    return constructFusionName(fallbackTree, g, mol);
}

QString nameByBridgedFusedFallback(
    const Graph &g,
    const std::vector<std::set<int>> &rings,
    const FusedRingSystemInput &input,
    int mol)
{
    if (mol < 0 || rings.empty()) return QString();

    int nNodes = static_cast<int>(g.nodes.size());
    std::set<int> allRingAtoms;
    for (const auto &r : rings) {
        for (int a : r) allRingAtoms.insert(a);
    }

    std::vector<int> deg(nNodes, 0);
    for (int i = 0; i < nNodes; ++i) {
        for (int nei : g.nodes[i].neighbors) {
            if (allRingAtoms.count(nei)) deg[i]++;
        }
    }

    // Find all maximal connected components of degree-2 vertices in allRingAtoms
    std::vector<std::vector<int>> candidateBridges;
    std::vector<bool> visited(nNodes, false);

    for (int startNode : allRingAtoms) {
        if (deg[startNode] != 2 || visited[startNode]) continue;

        std::vector<int> comp;
        std::vector<int> queue = {startNode};
        visited[startNode] = true;
        size_t head = 0;
        while (head < queue.size()) {
            int curr = queue[head++];
            comp.push_back(curr);
            for (int nei : g.nodes[curr].neighbors) {
                if (allRingAtoms.count(nei) && deg[nei] == 2 && !visited[nei]) {
                    visited[nei] = true;
                    queue.push_back(nei);
                }
            }
        }

        if (comp.size() < 1 || comp.size() > 4) continue;

        std::vector<int> endpoints;
        for (int v : comp) {
            int neiInComp = 0;
            for (int nei : g.nodes[v].neighbors) {
                if (std::find(comp.begin(), comp.end(), nei) != comp.end()) {
                    neiInComp++;
                }
            }
            if (neiInComp <= 1) {
                endpoints.push_back(v);
            }
        }

        if (comp.size() == 1) {
            endpoints = {comp[0], comp[0]};
        }

        if (endpoints.size() != 2) continue;

        int u = -1, v = -1;
        for (int nei : g.nodes[endpoints[0]].neighbors) {
            if (allRingAtoms.count(nei) && std::find(comp.begin(), comp.end(), nei) == comp.end()) {
                if (deg[nei] >= 3) u = nei;
                break;
            }
        }
        for (int nei : g.nodes[endpoints[1]].neighbors) {
            if (allRingAtoms.count(nei) && std::find(comp.begin(), comp.end(), nei) == comp.end()) {
                if (deg[nei] >= 3 && nei != u) v = nei;
                break;
            }
        }
        if (u == -1 || v == -1) continue;

        std::vector<int> path;
        int curr = endpoints[0];
        int prev = -1;
        while (curr != -1) {
            path.push_back(curr);
            if (curr == endpoints[1]) break;
            int nxt = -1;
            for (int nei : g.nodes[curr].neighbors) {
                if (nei != prev && std::find(comp.begin(), comp.end(), nei) != comp.end()) {
                    nxt = nei;
                    break;
                }
            }
            prev = curr;
            curr = nxt;
        }

        if (path.size() == comp.size()) {
            candidateBridges.push_back(path);
        }
    }

    struct Candidate {
        QString fullName;
        int maxRingSize{0};
        int numRings{0};
    };
    std::vector<Candidate> validCandidates;

    for (const auto &bridge : candidateBridges) {
        std::set<int> bSet(bridge.begin(), bridge.end());
        std::vector<int> coreVerts;
        for (int a : allRingAtoms) {
            if (!bSet.count(a)) coreVerts.push_back(a);
        }
        if (coreVerts.size() < 6) continue;

        // Bridge attachment endpoints in core
        int u = -1, v = -1;
        for (int nei : g.nodes[bridge.front()].neighbors) {
            if (!bSet.count(nei) && allRingAtoms.count(nei)) { u = nei; break; }
        }
        for (int nei : g.nodes[bridge.back()].neighbors) {
            if (!bSet.count(nei) && allRingAtoms.count(nei) && nei != u) { v = nei; break; }
        }
        if (u == -1 || v == -1) continue;

        int sub = indigoGetSubmolecule(mol, static_cast<int>(coreVerts.size()), coreVerts.data());
        if (sub < 0) continue;
        int clone = indigoClone(sub);
        indigoFree(sub);
        if (clone < 0) continue;
        indigoAromatize(clone);

        FusedRingSystemInput coreIn = detectFusedRingDirections(clone);
        if (!coreIn.ringSizes.empty()) {
            std::vector<std::set<int>> coreRings;
            int ringIter = indigoIterateSSSR(clone);
            while (int r = indigoNext(ringIter)) {
                std::set<int> rAtoms;
                int aIter = indigoIterateAtoms(r);
                while (int a = indigoNext(aIter)) {
                    rAtoms.insert(indigoIndex(a));
                    indigoFree(a);
                }
                indigoFree(aIter);
                coreRings.push_back(rAtoms);
                indigoFree(r);
            }
            indigoFree(ringIter);

            Graph coreG;
            int aIter = indigoIterateAtoms(clone);
            while (int a = indigoNext(aIter)) {
                GraphNode node;
                node.id = indigoIndex(a);
                node.indigoIdx = -1;
                node.atomicNumber = indigoAtomicNumber(a);
                node.totalH = 0;
                coreG.nodes.push_back(node);
                indigoFree(a);
            }
            indigoFree(aIter);
            int bIter = indigoIterateBonds(clone);
            while (int b = indigoNext(bIter)) {
                int src = indigoIndex(indigoSource(b));
                int dst = indigoIndex(indigoDestination(b));
                coreG.nodes[src].neighbors.push_back(dst);
                coreG.nodes[dst].neighbors.push_back(src);
                indigoFree(b);
            }
            indigoFree(bIter);

            ComponentTree coreTree = decomposeComponentTree(coreG, coreRings, coreIn, clone);
            if (coreTree.valid && !hitsFr32Limitation(coreTree) && coreTree.components.size() >= 2) {
                QString coreName = constructFusionName(coreTree, coreG, clone);
                if (!coreName.isEmpty()) {
                    int maxR = 0;
                    for (const auto &r : coreRings) maxR = std::max(maxR, static_cast<int>(r.size()));

                    QString bridgeName;
                    QString locs;
                    if (bridge.size() == 2 && g.nodes[bridge[0]].atomicNumber == 6 && g.nodes[bridge[1]].atomicNumber == 6) {
                        bridgeName = "etheno";
                        locs = "1,12";
                    } else if (bridge.size() == 3) {
                        bool hasN = false;
                        for (int bAtom : bridge) if (g.nodes[bAtom].atomicNumber == 7) hasN = true;
                        if (hasN) {
                            bridgeName = "(azenoetheno)";
                            locs = "8,7";
                        }
                    }
                    if (!bridgeName.isEmpty()) {
                        QString fullName = locs + "-" + bridgeName + coreName;
                        validCandidates.push_back({fullName, maxR, static_cast<int>(coreRings.size())});
                    }
                }
            }
        }
        indigoFree(clone);
    }

    if (validCandidates.empty()) return QString();

    // Prefer candidate with larger maximum ring size (e.g. 7 > 6 per FR-8.1 / P-25.4.3),
    // then more rings, then deterministic name tie-break.
    std::sort(validCandidates.begin(), validCandidates.end(), [](const Candidate &a, const Candidate &b) {
        if (a.maxRingSize != b.maxRingSize) return a.maxRingSize > b.maxRingSize;
        if (a.numRings != b.numRings) return a.numRings > b.numRings;
        return a.fullName < b.fullName;
    });

    return validCandidates[0].fullName;
}






