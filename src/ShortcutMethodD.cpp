#include "ShortcutMethodD.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gurobi_c++.h"
#include "GurobiSharedEnv.h"
#include "PhysicalConstants.h"
#include "RingLayout.h"
#include "ShortcutExport.h"
#include "ShortcutGrid.h"
#include "ShortcutOrchestrator.h"
#include "ShortcutRouter.h"
#include "ShortcutTypes.h"

namespace {

constexpr double kMethodDMilpTimeLimitSec = 120.0;

struct JointCandidate {
    int srcId = 0;
    int destId = 0;
    ShortcutPath path;
    double length = 0.0;
    int bendCount = 0;
    /// Indices of other candidates whose geometry intersects this one.
    std::vector<int> cross;
};

double manhattan(const Node& a, const Node& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

bool areTourNeighbors(int a, int b, const std::vector<int>& tour) {
    const int n = (int)tour.size();
    if (n < 2) return false;
    for (int i = 0; i < n; ++i) {
        if (tour[i] != a) continue;
        if (tour[(i + 1) % n] == b || tour[(i - 1 + n) % n] == b)
            return true;
    }
    return false;
}

double minRingPathLength(
        int a,
        int b,
        const std::vector<int>& tour,
        const std::vector<Node>& nodes) {
    const int n = (int)tour.size();
    int pa = -1, pb = -1;
    for (int i = 0; i < n; ++i) {
        if (tour[i] == a) pa = i;
        if (tour[i] == b) pb = i;
    }
    if (pa < 0 || pb < 0 || pa == pb)
        return 0.0;

    auto arcLen = [&](int from, int to, int step) {
        double sum = 0.0;
        int k = from;
        while (k != to) {
            const int nk = (k + step + n) % n;
            sum += manhattan(nodes[tour[k]], nodes[tour[nk]]);
            k = nk;
        }
        return sum;
    };
    return std::min(arcLen(pa, pb, +1), arcLen(pa, pb, -1));
}

Shortcut toIntersectShortcut(const JointCandidate& c) {
    Shortcut sc;
    sc.from = c.srcId;
    sc.to = c.destId;
    sc.approx_length = c.length;
    sc.bend_count = c.bendCount;
    sc.path = pathVerticesToSegments(c.path.vertices);
    return sc;
}

std::vector<JointCandidate> phase1BuildCandidates(
        const std::vector<Node>& nodes,
        const RingLayout& ring,
        double sMin,
        bool quiet) {
    const int N = (int)nodes.size();
    const auto& tour = ring.tour;
    const auto& routing = ring.routing;
    const auto ringSegments = ShortcutGrid::collectRingSegments(routing);

    std::vector<JointCandidate> candidates;
    int pairsConsidered = 0;
    int pairsDominated = 0;
    int pairsRouted = 0;
    int pathsKept = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            ++pairsConsidered;
            if (areTourNeighbors(i, j, tour))
                continue;

            const double manh = manhattan(nodes[i], nodes[j]);
            const double ringLen = minRingPathLength(i, j, tour, nodes);
            if (manh + MILP_EPS >= ringLen) {
                ++pairsDominated;
                continue;
            }

            ++pairsRouted;
            auto gridResult = ShortcutGrid::buildCached(
                nodes, ringSegments, tour, routing, i, j, sMin);
            if (!gridResult.success)
                continue;

            ShortcutRouteResult route = ShortcutRouter::findPaths(
                gridResult.points, i, j, nodes, ringSegments,
                tour, routing, /*existingShortcuts=*/{},
                std::numeric_limits<double>::infinity());
            if (!route.success)
                continue;

            auto pushPath = [&](const ShortcutPath& path) {
                if (path.vertices.size() < 2)
                    return;
                JointCandidate c;
                c.srcId = i;
                c.destId = j;
                c.path = path;
                c.length = path.distance;
                c.bendCount = path.bendCount;
                candidates.push_back(std::move(c));
                ++pathsKept;
            };

            pushPath(route.primary);
            for (const ShortcutPath& alt : route.alternatives)
                pushPath(alt);
        }
    }

    const int C = (int)candidates.size();
    for (int a = 0; a < C; ++a) {
        const Shortcut sa = toIntersectShortcut(candidates[a]);
        for (int b = a + 1; b < C; ++b) {
            const Shortcut sb = toIntersectShortcut(candidates[b]);
            if (ShortcutRouter::shortcutsIntersect(sa, sb)) {
                candidates[a].cross.push_back(b);
                candidates[b].cross.push_back(a);
            }
        }
    }

    if (!quiet) {
        std::cout << std::fixed << std::setprecision(1)
                  << "[Method D] Phase 1: pairs=" << pairsConsidered
                  << " dominated=" << pairsDominated
                  << " routed=" << pairsRouted
                  << " candidates=" << pathsKept << "\n";
    }
    return candidates;
}

struct Phase2SolveResult {
    bool success = false;
    bool abortedByBound = false;
    double W = 0.0;
    std::vector<int> selected;  // candidate indices
    /// For each selected index into `selected`, first demand using it (−1 if none).
    std::vector<int> firstDemandPerSelected;
};

class DualAbortCallback : public GRBCallback {
public:
    DualAbortCallback(double abortGe, double abortOnlyGt)
        : abortGe_(abortGe), abortOnlyGt_(abortOnlyGt) {}
    bool aborted = false;

protected:
    void callback() override {
        try {
            if (where != GRB_CB_MIP) return;
            if (!std::isfinite(abortGe_)) return;
            const double bound = getDoubleInfo(GRB_CB_MIP_OBJBND);
            if (bound + MILP_EPS >= abortGe_ && bound > abortOnlyGt_ + MILP_EPS) {
                aborted = true;
                abort();
            }
        } catch (...) {
        }
    }

private:
    double abortGe_ = 0.0;
    double abortOnlyGt_ = 0.0;
};

Phase2SolveResult phase2SolveMilp(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const RingLayout& ring,
        const std::vector<JointCandidate>& candidates,
        double timeLimitSec,
        bool quiet,
        double abortIfObjBoundGe,
        double abortOnlyIfObjBoundGt,
        ShortcutUsageMode usageMode) {
    Phase2SolveResult out;
    const int N = (int)nodes.size();
    const int Q = (int)D.demands.size();
    const int C = (int)candidates.size();
    const int tourLen = (int)ring.tour.size();
    if (N == 0 || tourLen < 2 || Q == 0) {
        out.success = true;
        out.W = 0.0;
        return out;
    }

    // Ring undirected edges → 2 directed arcs each.
    struct RingArc {
        int from = 0;
        int to = 0;
        double weight = 0.0;
    };
    std::vector<RingArc> ringArcs;
    ringArcs.reserve(2 * tourLen);
    for (int k = 0; k < tourLen; ++k) {
        const int u = ring.tour[k];
        const int v = ring.tour[(k + 1) % tourLen];
        const double w = manhattan(nodes[u], nodes[v]);
        ringArcs.push_back({u, v, w});
        ringArcs.push_back({v, u, w});
    }
    const int R = (int)ringArcs.size();

    try {
        GRBModel model(sharedGurobiEnv());
        model.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);

        GRBVar W = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "W");

        std::vector<GRBVar> z(C);
        for (int c = 0; c < C; ++c) {
            z[c] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + std::to_string(c));
        }

        // fRing[q][e], fSc[q][c][dir] with dir 0 = src→dest, 1 = dest→src
        std::vector<std::vector<GRBVar>> fRing(Q, std::vector<GRBVar>(R));
        std::vector<std::vector<std::array<GRBVar, 2>>> fSc(
            Q, std::vector<std::array<GRBVar, 2>>(C));

        for (int q = 0; q < Q; ++q) {
            for (int e = 0; e < R; ++e) {
                fRing[q][e] = model.addVar(
                    0.0, 1.0, 0.0, GRB_CONTINUOUS,
                    "fr_" + std::to_string(q) + "_" + std::to_string(e));
            }
            for (int c = 0; c < C; ++c) {
                fSc[q][c][0] = model.addVar(
                    0.0, 1.0, 0.0, GRB_CONTINUOUS,
                    "fs_" + std::to_string(q) + "_" + std::to_string(c) + "_0");
                fSc[q][c][1] = model.addVar(
                    0.0, 1.0, 0.0, GRB_CONTINUOUS,
                    "fs_" + std::to_string(q) + "_" + std::to_string(c) + "_1");
            }
        }
        model.update();

        // C1: at most one shortcut per node
        for (int i = 0; i < N; ++i) {
            GRBLinExpr sum = 0;
            bool any = false;
            for (int c = 0; c < C; ++c) {
                if (candidates[c].srcId == i || candidates[c].destId == i) {
                    sum += z[c];
                    any = true;
                }
            }
            if (any)
                model.addConstr(sum <= 1, "C1_" + std::to_string(i));
        }

        // C2b: cardinality — selected shortcut may cross at most one other selected
        for (int c = 0; c < C; ++c) {
            const int M = (int)candidates[c].cross.size();
            if (M == 0) continue;
            GRBLinExpr sum = 0;
            for (int q : candidates[c].cross)
                sum += z[q];
            // sum_{q in Cross(c)} z_q <= 1 + M*(1 - z_c)
            model.addConstr(
                sum <= 1.0 + static_cast<double>(M) * (1.0 - z[c]),
                "C2b_" + std::to_string(c));
        }

        // C3: flow conservation (match MILPSolver: out - in = +1 at src, -1 at dest)
        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int d = D.demands[q].second;
            for (int v = 0; v < N; ++v) {
                GRBLinExpr outExpr = 0;
                GRBLinExpr inExpr = 0;
                for (int e = 0; e < R; ++e) {
                    if (ringArcs[e].from == v) outExpr += fRing[q][e];
                    if (ringArcs[e].to == v) inExpr += fRing[q][e];
                }
                for (int c = 0; c < C; ++c) {
                    const int u = candidates[c].srcId;
                    const int w = candidates[c].destId;
                    if (u == v) {
                        outExpr += fSc[q][c][0];
                        inExpr += fSc[q][c][1];
                    }
                    if (w == v) {
                        outExpr += fSc[q][c][1];
                        inExpr += fSc[q][c][0];
                    }
                }
                int rhs = 0;
                if (v == s) rhs = +1;
                if (v == d) rhs = -1;
                model.addConstr(
                    outExpr - inExpr == rhs,
                    "C3_" + std::to_string(q) + "_" + std::to_string(v));
            }
        }

        // C4: coupling — demand may use candidate edge only if selected.
        // Exclusive: only the owner demand with matching start/end endpoints.
        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            for (int c = 0; c < C; ++c) {
                const bool ownerOk = (usageMode != ShortcutUsageMode::Exclusive)
                    || undirectedEndpointsMatch(
                        candidates[c].srcId, candidates[c].destId, s, t);
                if (ownerOk) {
                    model.addConstr(
                        fSc[q][c][0] <= z[c],
                        "C4a_" + std::to_string(q) + "_" + std::to_string(c));
                    model.addConstr(
                        fSc[q][c][1] <= z[c],
                        "C4b_" + std::to_string(q) + "_" + std::to_string(c));
                } else {
                    model.addConstr(
                        fSc[q][c][0] == 0,
                        "C4a_priv_" + std::to_string(q) + "_" + std::to_string(c));
                    model.addConstr(
                        fSc[q][c][1] == 0,
                        "C4b_priv_" + std::to_string(q) + "_" + std::to_string(c));
                }
            }
        }

        // C5: W >= path length (distance mm only; bends not in w_e)
        for (int q = 0; q < Q; ++q) {
            GRBLinExpr pathLen = 0;
            for (int e = 0; e < R; ++e)
                pathLen += fRing[q][e] * ringArcs[e].weight;
            for (int c = 0; c < C; ++c) {
                pathLen += (fSc[q][c][0] + fSc[q][c][1]) * candidates[c].length;
            }
            model.addConstr(W >= pathLen, "C5_" + std::to_string(q));
        }

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        if (timeLimitSec > 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimitSec);

        DualAbortCallback abortCb(abortIfObjBoundGe, abortOnlyIfObjBoundGt);
        if (std::isfinite(abortIfObjBoundGe))
            model.setCallback(&abortCb);

        model.optimize();

        if (abortCb.aborted) {
            out.abortedByBound = true;
            out.success = false;
            if (!quiet) {
                std::cout << "[Method D] Phase 2: dual-abort (bound cannot improve)\n";
            }
            return out;
        }

        const int status = model.get(GRB_IntAttr_Status);
        const int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount <= 0) {
            if (!quiet) {
                std::cerr << "[Method D] Phase 2: no feasible solution (status="
                          << status << ")\n";
            }
            return out;
        }

        out.success = true;
        out.W = W.get(GRB_DoubleAttr_X);

        for (int c = 0; c < C; ++c) {
            if (z[c].get(GRB_DoubleAttr_X) > 0.5)
                out.selected.push_back(c);
        }

        out.firstDemandPerSelected.assign(out.selected.size(), -1);
        for (size_t si = 0; si < out.selected.size(); ++si) {
            const int c = out.selected[si];
            for (int q = 0; q < Q; ++q) {
                const double f0 = fSc[q][c][0].get(GRB_DoubleAttr_X);
                const double f1 = fSc[q][c][1].get(GRB_DoubleAttr_X);
                if (f0 + f1 > 0.5) {
                    out.firstDemandPerSelected[si] = q;
                    break;
                }
            }
        }

        if (!quiet) {
            std::cout << std::fixed << std::setprecision(1)
                      << "[Method D] Phase 2: status=" << status
                      << " W=" << out.W << " mm, selected="
                      << out.selected.size() << " shortcuts\n";
        }
    } catch (GRBException& e) {
        std::cerr << "[Method D] Gurobi error: " << e.getMessage()
                  << " (code " << e.getErrorCode() << ")\n";
    }
    return out;
}

void markSelectedCrossings(
        std::vector<PlacedShortcut>& placed,
        const std::vector<JointCandidate>& candidates,
        const std::vector<int>& selected) {
    const int S = (int)selected.size();
    for (int i = 0; i < S; ++i) {
        const int ci = selected[i];
        for (int j = i + 1; j < S; ++j) {
            const int cj = selected[j];
            const auto& cross = candidates[ci].cross;
            if (std::find(cross.begin(), cross.end(), cj) == cross.end())
                continue;
            placed[i].everCrossed = true;
            placed[j].everCrossed = true;
            if (placed[i].crossedShortcutIdx < 0)
                placed[i].crossedShortcutIdx = j;
            if (placed[j].crossedShortcutIdx < 0)
                placed[j].crossedShortcutIdx = i;
        }
    }
}

std::set<int> shortcutDemandIndicesForExport(
        const std::vector<PlacedShortcut>& shortcuts,
        const DemandMatrix& D,
        const std::vector<Node>& nodes,
        const MILPSolveResult& layout) {
    std::set<int> excluded;
    for (int q = 0; q < (int)D.demands.size(); ++q) {
        if (demandHasShortcut(
                q, D, layout.demandDistance, shortcuts, nodes, layout.tour))
            excluded.insert(q);
    }
    return excluded;
}

}  // namespace

ShortcutMethodResult runMethodDJointShortcuts(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layoutFixed,
        double sMin,
        const ShortcutMethodOptions& options) {
    ScopedShortcutUsageMode usageScope(options.usageMode);
    ShortcutMethodResult out;
    out.layout = layoutFixed;
    out.termination = ShortcutTermination::Case2Constraints;

    if (!layoutFixed.success || layoutFixed.tour.empty()) {
        if (!options.quiet)
            std::cerr << "[Method D] Invalid layout — aborting.\n";
        out.globalW = layoutFixed.W;
        return out;
    }

    if (!options.skipExports) {
        exportRingDemandIL(layoutFixed, D, "D");
        exportRingSnapshot(layoutFixed, nodes, "ring_before_D.csv");
        exportDemandsSnapshot(layoutFixed, D, nodes, "demands_ring_before_D.csv");
    }

    const RingLayout ring = buildRingLayout(nodes, layoutFixed);
    const int numDemands = (int)D.demands.size();
    const double plainW = computeGlobalW(
        numDemands, D, layoutFixed.demandDistance, layoutFixed.demandBendCount,
        /*shortcuts=*/{}, nodes, ring.tour);

    if (!options.quiet) {
        std::cout << "\n======================================================\n";
        std::cout << "Method D — JointShortcuts (fixed ring, global MILP)\n";
        std::cout << "======================================================\n";
        std::cout << std::fixed << std::setprecision(1)
                  << "[Method D] Plain ring WC: " << plainW << " mm\n";
    }

    const std::vector<JointCandidate> candidates =
        phase1BuildCandidates(nodes, ring, sMin, options.quiet);

    if (candidates.empty()) {
        out.shortcuts.clear();
        out.globalW = plainW;
        if (!options.quiet) {
            std::cout << "[Method D] No feasible candidates — returning plain ring W.\n";
        }
        if (!options.skipExports) {
            exportFinalCSVs(out.layout, out.shortcuts, D, {}, "D", nodes);
        }
        return out;
    }

    Phase2SolveResult milp = phase2SolveMilp(
        nodes, D, ring, candidates, kMethodDMilpTimeLimitSec, options.quiet,
        options.abortIfObjBoundGe, options.abortOnlyIfObjBoundGt, options.usageMode);

    if (milp.abortedByBound) {
        out.shortcuts.clear();
        out.globalW = std::numeric_limits<double>::infinity();
        if (!options.quiet)
            std::cout << "[Method D] Dual-aborted — no incumbent update.\n";
        return out;
    }

    if (!milp.success) {
        out.shortcuts.clear();
        out.globalW = plainW;
        if (!options.quiet)
            std::cerr << "[Method D] MILP failed — returning plain ring W.\n";
        if (!options.skipExports) {
            exportFinalCSVs(out.layout, out.shortcuts, D, {}, "D", nodes);
        }
        return out;
    }

    std::vector<PlacedShortcut> placed;
    placed.reserve(milp.selected.size());
    for (size_t si = 0; si < milp.selected.size(); ++si) {
        const int c = milp.selected[si];
        const JointCandidate& cand = candidates[c];
        PlacedShortcut ps;
        ps.demandIdx = milp.firstDemandPerSelected[si];
        ps.srcId = cand.srcId;
        ps.destId = cand.destId;
        ps.path = cand.path;
        ps.totalIL = reportingShortcutIL_dB(cand.length, cand.bendCount);
        ps.everCrossed = false;
        ps.crossedShortcutIdx = -1;
        placed.push_back(std::move(ps));
    }
    markSelectedCrossings(placed, candidates, milp.selected);

    out.shortcuts = std::move(placed);
    out.globalW = milp.W;

    if (!options.skipExports) {
        const std::set<int> excluded =
            shortcutDemandIndicesForExport(out.shortcuts, D, nodes, out.layout);
        exportFinalCSVs(out.layout, out.shortcuts, D, excluded, "D", nodes);
    }

    if (!options.quiet) {
        std::cout << std::fixed << std::setprecision(1)
                  << "\n[Method D] Final global W: " << out.globalW << " mm\n";
        std::cout << "[Method D] Shortcuts placed: " << out.shortcuts.size() << "\n";
    }

    return out;
}
