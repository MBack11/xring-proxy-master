#include "ProxyMasterMilp.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "gurobi_c++.h"
#include "GurobiSharedEnv.h"
#include "PhysicalConstants.h"
#include "RingDemandRouting.h"
#include "RingGeometry.h"

namespace {

constexpr int kMtzThreshold = 20;

class SubtourLazyCallback : public GRBCallback {
public:
    SubtourLazyCallback(int n, std::vector<GRBVar>& b) : N(n), b_(b) {}

protected:
    void callback() override {
        try {
            if (where != GRB_CB_MIPSOL) return;
            std::vector<int> next(N, -1);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    if (getSolution(b_[idx(i, j)]) > 0.5) {
                        next[i] = j;
                        break;
                    }
                }
            std::vector<bool> visited(N, false);
            for (int start = 0; start < N; ++start) {
                if (visited[start]) continue;
                std::vector<int> tour;
                int cur = start;
                while (cur != -1 && !visited[cur]) {
                    visited[cur] = true;
                    tour.push_back(cur);
                    cur = next[cur];
                }
                if (!tour.empty() && (int)tour.size() < N) {
                    GRBLinExpr expr = 0;
                    for (int k = 0; k < (int)tour.size(); ++k) {
                        const int from = tour[k];
                        const int to = tour[(k + 1) % tour.size()];
                        expr += b_[idx(from, to)];
                    }
                    addLazy(expr <= (int)tour.size() - 1);
                }
            }
        } catch (GRBException& e) {
            std::cerr << "Subtour callback: " << e.getMessage() << "\n";
        }
    }

private:
    int N = 0;
    std::vector<GRBVar>& b_;
    int idx(int i, int j) const { return i * N + j; }
};

MasterSolution extractSolution(
        GRBModel& model,
        const std::vector<GRBVar>& b,
        const std::vector<GRBVar>& r,
        const std::vector<GRBVar>& z,
        const std::vector<ShortcutPairIndex>& pairs,
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        int N) {
    MasterSolution sol;
    sol.Wproxy = model.get(GRB_DoubleAttr_ObjVal);

    auto bIdx = [N](int i, int j) { return i * N + j; };

    std::vector<int> succ(N, -1);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            if (b[bIdx(i, j)].get(GRB_DoubleAttr_X) > 0.5)
                succ[i] = j;
        }

    sol.layout.success = true;
    sol.layout.tour.clear();
    int cur = 0;
    for (int k = 0; k < N; ++k) {
        sol.layout.tour.push_back(cur);
        cur = succ[cur];
        if (cur < 0) break;
    }

    sol.layout.tourEdges.clear();
    for (int k = 0; k < (int)sol.layout.tour.size(); ++k) {
        const int u = sol.layout.tour[k];
        const int v = sol.layout.tour[(k + 1) % sol.layout.tour.size()];
        const int opt = static_cast<int>(r[bIdx(u, v)].get(GRB_DoubleAttr_X) + 0.5);
        sol.layout.tourEdges.push_back({u, v, opt});
    }

    for (int c = 0; c < (int)z.size(); ++c) {
        if (z[c].get(GRB_DoubleAttr_X) > 0.5)
            sol.selectedShortcutIndices.push_back(c);
    }

    applyShortestArcDemandRouting(nodes, D, sol.layout);
    sol.layout.W = sol.Wproxy;
    return sol;
}

}  // namespace

std::vector<std::pair<int, int>> ringUndirectedEdges(
        const MILPSolveResult& layout) {
    std::set<std::pair<int, int>> edges;
    for (const auto& te : layout.tourEdges) {
        const int a = std::min(te.from, te.to);
        const int b = std::max(te.from, te.to);
        edges.insert({a, b});
    }
    return {edges.begin(), edges.end()};
}

std::string ringTourKey(const MILPSolveResult& layout) {
    std::string s;
    for (const auto& [a, b] : ringUndirectedEdges(layout))
        s += std::to_string(a) + "-" + std::to_string(b) + ";";
    return s;
}

FixedRingProxyResult solveProxyMasterFixedRing(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layoutFixed,
        double timeLimitSec,
        bool quiet,
        ShortcutUsageMode usageMode) {
    FixedRingProxyResult out;
    const int N = (int)nodes.size();
    const int Q = (int)D.demands.size();
    if (N < 3 || Q == 0 || !layoutFixed.success || layoutFixed.tourEdges.empty()) {
        out.message = "invalid layout";
        return out;
    }

    const std::vector<ShortcutPairIndex> pairs = buildMasterShortcutPairs(nodes);
    const std::vector<std::pair<int, int>> xPairs =
        computeShortcutCrossingPairs(nodes, pairs);
    const RingGeometry geom(nodes);

    auto bIdx = [N](int i, int j) { return i * N + j; };

    std::vector<std::vector<int>> bFixed(N, std::vector<int>(N, 0));
    for (const auto& te : layoutFixed.tourEdges)
        bFixed[te.from][te.to] = 1;

    std::set<std::pair<int, int>> tourUndir;
    for (const auto& e : ringUndirectedEdges(layoutFixed))
        tourUndir.insert(e);

    try {
        GRBModel model(sharedGurobiEnv());
        model.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);
        if (timeLimitSec > 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimitSec);

        GRBVar W = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "W");

        const int C = (int)pairs.size();
        std::vector<GRBVar> z(C);
        for (int c = 0; c < C; ++c)
            z[c] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + std::to_string(c));

        std::vector<std::vector<GRBVar>> fRing(Q, std::vector<GRBVar>(N * N));
        std::vector<std::vector<std::array<GRBVar, 2>>> fSc(
            Q, std::vector<std::array<GRBVar, 2>>(C));
        for (int q = 0; q < Q; ++q) {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    fRing[q][bIdx(i, j)] = model.addVar(
                        0.0, 1.0, 0.0, GRB_CONTINUOUS,
                        "fr_" + std::to_string(q) + "_" + std::to_string(i)
                            + "_" + std::to_string(j));
                }
            for (int c = 0; c < C; ++c)
                for (int d = 0; d < 2; ++d)
                    fSc[q][c][d] = model.addVar(
                        0.0, 1.0, 0.0, GRB_CONTINUOUS,
                        "fs_" + std::to_string(q) + "_" + std::to_string(c)
                            + "_" + std::to_string(d));
        }
        model.update();

        for (int v = 0; v < N; ++v) {
            GRBLinExpr occ = 0;
            for (int c = 0; c < C; ++c)
                if (pairs[c].i == v || pairs[c].j == v)
                    occ += z[c];
            model.addConstr(occ <= 1, "M3_" + std::to_string(v));
        }

        for (int c = 0; c < C; ++c) {
            const int i = pairs[c].i, j = pairs[c].j;
            const int a = std::min(i, j), b = std::max(i, j);
            if (tourUndir.count({a, b}))
                model.addConstr(z[c] == 0, "M2_" + std::to_string(c));
        }

        for (const auto& [p, q] : xPairs)
            model.addConstr(z[p] + z[q] <= 1, "M2b");

        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    model.addConstr(
                        fRing[q][bIdx(i, j)] <= bFixed[i][j] + bFixed[j][i],
                        "M4r");
                }
            for (int c = 0; c < C; ++c) {
                const bool ownerOk = (usageMode != ShortcutUsageMode::Exclusive)
                    || undirectedEndpointsMatch(pairs[c].i, pairs[c].j, s, t);
                for (int d = 0; d < 2; ++d) {
                    if (ownerOk)
                        model.addConstr(fSc[q][c][d] <= z[c], "M4s");
                    else
                        model.addConstr(fSc[q][c][d] == 0, "M4s_priv");
                }
            }
        }

        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            for (int v = 0; v < N; ++v) {
                GRBLinExpr delta = 0;
                for (int j = 0; j < N; ++j) {
                    if (v == j) continue;
                    delta += fRing[q][bIdx(v, j)];
                    delta -= fRing[q][bIdx(j, v)];
                }
                for (int c = 0; c < C; ++c) {
                    const int i = pairs[c].i, j = pairs[c].j;
                    if (v == i) {
                        delta += fSc[q][c][0];
                        delta -= fSc[q][c][1];
                    } else if (v == j) {
                        delta += fSc[q][c][1];
                        delta -= fSc[q][c][0];
                    }
                }
                double rhs = 0.0;
                if (v == s) rhs = 1.0;
                else if (v == t) rhs = -1.0;
                model.addConstr(delta == rhs, "M5");
            }
        }

        for (int q = 0; q < Q; ++q) {
            GRBLinExpr len = 0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    len += geom.manhattanDist(i, j) * fRing[q][bIdx(i, j)];
                }
            for (int c = 0; c < C; ++c)
                for (int d = 0; d < 2; ++d)
                    len += pairs[c].delta * fSc[q][c][d];
            model.addConstr(W >= len, "M6");
        }

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.optimize();

        if (model.get(GRB_IntAttr_SolCount) <= 0) {
            out.message = "no solution, status="
                + std::to_string(model.get(GRB_IntAttr_Status));
            return out;
        }

        out.success = true;
        out.Wproxy = W.get(GRB_DoubleAttr_X);
        for (int c = 0; c < C; ++c)
            if (z[c].get(GRB_DoubleAttr_X) > 0.5)
                out.selectedShortcutIndices.push_back(c);
        out.message = "ok";
    } catch (GRBException& e) {
        out.message = std::string("GRB: ") + e.getMessage();
    }
    return out;
}

struct ProxyMasterMilpSession::Impl {
    const std::vector<Node>& nodes;
    const DemandMatrix& D;
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared;
    int N = 0;
    int Q = 0;
    std::vector<ShortcutPairIndex> pairs;
    std::vector<std::pair<int, int>> xPairs;
    RingGeometry geom;
    bool useLazySubtour = false;
    bool built = false;
    size_t cutsApplied = 0;

    std::unique_ptr<GRBModel> model;
    GRBVar W;
    std::vector<GRBVar> b;
    std::vector<GRBVar> r;
    std::vector<GRBVar> z;
    std::unique_ptr<SubtourLazyCallback> subtourCb;

    Impl(const std::vector<Node>& n, const DemandMatrix& d, ShortcutUsageMode mode)
        : nodes(n), D(d), usageMode(mode), N((int)n.size()),
          Q((int)d.demands.size()), geom(n) {
        pairs = buildMasterShortcutPairs(nodes);
        xPairs = computeShortcutCrossingPairs(nodes, pairs);
        useLazySubtour = N > kMtzThreshold;
    }

    int bIdx(int i, int j) const { return i * N + j; }

    void ensureBuilt() {
        if (built) return;
        if (N < 3 || Q == 0)
            throw std::runtime_error("empty instance");

        model = std::make_unique<GRBModel>(sharedGurobiEnv());
        model->set(GRB_IntParam_OutputFlag, 0);
        if (useLazySubtour)
            model->set(GRB_IntParam_LazyConstraints, 1);

        W = model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "W");

        b.assign(N * N, GRBVar());
        r.assign(N * N, GRBVar());
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                b[bIdx(i, j)] = model->addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "b_" + std::to_string(i) + "_" + std::to_string(j));
                r[bIdx(i, j)] = model->addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "r_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        const int C = (int)pairs.size();
        z.resize(C);
        for (int c = 0; c < C; ++c)
            z[c] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + std::to_string(c));

        std::vector<std::vector<GRBVar>> fRing(Q, std::vector<GRBVar>(N * N));
        std::vector<std::vector<std::array<GRBVar, 2>>> fSc(
            Q, std::vector<std::array<GRBVar, 2>>(C));
        for (int q = 0; q < Q; ++q) {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    fRing[q][bIdx(i, j)] = model->addVar(
                        0.0, 1.0, 0.0, GRB_CONTINUOUS,
                        "fr_" + std::to_string(q) + "_" + std::to_string(i)
                            + "_" + std::to_string(j));
                }
            for (int c = 0; c < C; ++c)
                for (int d = 0; d < 2; ++d)
                    fSc[q][c][d] = model->addVar(
                        0.0, 1.0, 0.0, GRB_CONTINUOUS,
                        "fs_" + std::to_string(q) + "_" + std::to_string(c)
                            + "_" + std::to_string(d));
        }
        model->update();

        for (int i = 0; i < N; ++i) {
            GRBLinExpr outE = 0, inE = 0;
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                outE += b[bIdx(i, j)];
                inE += b[bIdx(j, i)];
            }
            model->addConstr(outE == 1, "C1_" + std::to_string(i));
            model->addConstr(inE == 1, "C2_" + std::to_string(i));
        }

        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                model->addConstr(
                    b[bIdx(i, j)] + b[bIdx(j, i)] <= 1,
                    "C3_" + std::to_string(i) + "_" + std::to_string(j));

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                model->addConstr(r[bIdx(i, j)] <= b[bIdx(i, j)],
                    "r_le_b_" + std::to_string(i) + "_" + std::to_string(j));
            }

        if (!useLazySubtour) {
            std::vector<GRBVar> u(N);
            for (int i = 1; i < N; ++i)
                u[i] = model->addVar(1.0, N - 1, 0.0, GRB_CONTINUOUS,
                    "u_" + std::to_string(i));
            model->update();
            for (int i = 1; i < N; ++i)
                for (int j = 1; j < N; ++j) {
                    if (i == j) continue;
                    model->addConstr(
                        u[i] - u[j] + (N - 1) * b[bIdx(i, j)] <= N - 2,
                        "MTZ_" + std::to_string(i) + "_" + std::to_string(j));
                }
        }

        for (const RingCrossingPair& cp : geom.crossingPairs()) {
            GRBLinExpr lhs = b[bIdx(cp.i, cp.j)] + b[bIdx(cp.k, cp.l)];
            if (cp.opt1 == 0) lhs += 1 - r[bIdx(cp.i, cp.j)];
            else              lhs += r[bIdx(cp.i, cp.j)];
            if (cp.opt2 == 0) lhs += 1 - r[bIdx(cp.k, cp.l)];
            else              lhs += r[bIdx(cp.k, cp.l)];
            model->addConstr(lhs <= 3, "C5");
        }

        for (const RingPassThroughTriple& pt : geom.passThroughTriples()) {
            GRBLinExpr lhs = b[bIdx(pt.i, pt.j)];
            if (pt.opt == 0) lhs += 1 - r[bIdx(pt.i, pt.j)];
            else             lhs += r[bIdx(pt.i, pt.j)];
            model->addConstr(lhs <= 1, "C5b");
        }

        for (int v = 0; v < N; ++v) {
            GRBLinExpr occ = 0;
            for (int c = 0; c < C; ++c)
                if (pairs[c].i == v || pairs[c].j == v)
                    occ += z[c];
            model->addConstr(occ <= 1, "M3_" + std::to_string(v));
        }

        for (int c = 0; c < C; ++c) {
            const int i = pairs[c].i, j = pairs[c].j;
            model->addConstr(
                z[c] + b[bIdx(i, j)] + b[bIdx(j, i)] <= 1,
                "M2_" + std::to_string(c));
        }

        for (const auto& [p, q] : xPairs)
            model->addConstr(z[p] + z[q] <= 1, "M2b_" + std::to_string(p)
                + "_" + std::to_string(q));

        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    model->addConstr(
                        fRing[q][bIdx(i, j)]
                            <= b[bIdx(i, j)] + b[bIdx(j, i)],
                        "M4r_" + std::to_string(q) + "_" + std::to_string(i)
                            + "_" + std::to_string(j));
                }
            for (int c = 0; c < C; ++c) {
                const bool ownerOk = (usageMode != ShortcutUsageMode::Exclusive)
                    || undirectedEndpointsMatch(pairs[c].i, pairs[c].j, s, t);
                for (int d = 0; d < 2; ++d) {
                    if (ownerOk) {
                        model->addConstr(
                            fSc[q][c][d] <= z[c],
                            "M4s_" + std::to_string(q) + "_" + std::to_string(c)
                                + "_" + std::to_string(d));
                    } else {
                        model->addConstr(
                            fSc[q][c][d] == 0,
                            "M4s_priv_" + std::to_string(q) + "_"
                                + std::to_string(c) + "_" + std::to_string(d));
                    }
                }
            }
        }

        for (int q = 0; q < Q; ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            for (int v = 0; v < N; ++v) {
                GRBLinExpr delta = 0;
                for (int j = 0; j < N; ++j) {
                    if (v == j) continue;
                    delta += fRing[q][bIdx(v, j)];
                    delta -= fRing[q][bIdx(j, v)];
                }
                for (int c = 0; c < C; ++c) {
                    const int i = pairs[c].i, j = pairs[c].j;
                    if (v == i) {
                        delta += fSc[q][c][0];
                        delta -= fSc[q][c][1];
                    } else if (v == j) {
                        delta += fSc[q][c][1];
                        delta -= fSc[q][c][0];
                    }
                }
                double rhs = 0.0;
                if (v == s) rhs = 1.0;
                else if (v == t) rhs = -1.0;
                model->addConstr(delta == rhs, "M5_" + std::to_string(q)
                    + "_" + std::to_string(v));
            }
        }

        for (int q = 0; q < Q; ++q) {
            GRBLinExpr len = 0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    len += geom.manhattanDist(i, j) * fRing[q][bIdx(i, j)];
                }
            for (int c = 0; c < C; ++c)
                for (int d = 0; d < 2; ++d)
                    len += pairs[c].delta * fSc[q][c][d];
            model->addConstr(W >= len, "M6_" + std::to_string(q));
        }

        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        subtourCb = std::make_unique<SubtourLazyCallback>(N, b);
        built = true;
    }

    void appendCuts(const std::vector<MasterNogoodCut>& nogoodCuts) {
        while (cutsApplied < nogoodCuts.size()) {
            const size_t k = cutsApplied;
            GRBLinExpr lhsForward = 0;
            for (const auto& [i, j] : nogoodCuts[k].directedArcs)
                lhsForward += b[bIdx(i, j)];
            model->addConstr(lhsForward <= N - 1,
                "nogood_fwd_" + std::to_string(k));

            GRBLinExpr lhsReverse = 0;
            for (const auto& [i, j] : nogoodCuts[k].directedArcs)
                lhsReverse += b[bIdx(j, i)];
            model->addConstr(lhsReverse <= N - 1,
                "nogood_rev_" + std::to_string(k));
            ++cutsApplied;
        }
        model->update();
    }

    ProxyMasterMilpResult solve(
            const std::vector<MasterNogoodCut>& nogoodCuts,
            const ProxyMasterMilpOptions& options) {
        ProxyMasterMilpResult out;
        out.numXPairs = (int)xPairs.size();
        if (N < 3 || Q == 0) {
            out.message = "empty instance";
            return out;
        }

        try {
            ensureBuilt();
            appendCuts(nogoodCuts);

            model->set(GRB_IntParam_OutputFlag, options.quiet ? 0 : 1);
            if (options.timeLimitSec > 0.0)
                model->set(GRB_DoubleParam_TimeLimit, options.timeLimitSec);
            else
                model->set(GRB_DoubleParam_TimeLimit, GRB_INFINITY);
            if (options.mipGap >= 0.0)
                model->set(GRB_DoubleParam_MIPGap, options.mipGap);

            const bool usePool = options.poolSolutions > 1;
            if (usePool) {
                model->set(GRB_IntParam_PoolSearchMode, 2);
                model->set(GRB_IntParam_PoolSolutions, options.poolSolutions);
                model->set(GRB_DoubleParam_PoolGap, 0.0);
            } else {
                model->set(GRB_IntParam_PoolSearchMode, 0);
                model->set(GRB_IntParam_PoolSolutions, 1);
            }

            if (useLazySubtour)
                model->setCallback(subtourCb.get());
            else
                model->setCallback(nullptr);

            model->optimize();

            if (model->get(GRB_IntAttr_SolCount) <= 0) {
                out.message = "no solution, status="
                    + std::to_string(model->get(GRB_IntAttr_Status));
                return out;
            }

            out.success = true;
            out.Wproxy = model->get(GRB_DoubleAttr_ObjVal);
            out.solverStatus = model->get(GRB_IntAttr_Status);
            out.objBound = model->get(GRB_DoubleAttr_ObjBound);
            out.mipGapAtStop = model->get(GRB_DoubleAttr_MIPGap);
            out.runtimeSec = model->get(GRB_DoubleAttr_Runtime);
            out.message = "ok";

            if (usePool) {
                const int nSol = model->get(GRB_IntAttr_SolCount);
                for (int s = 0; s < nSol; ++s) {
                    model->set(GRB_IntParam_SolutionNumber, s);
                    const double poolObj = model->get(GRB_DoubleAttr_PoolObjVal);
                    if (poolObj > out.Wproxy + MILP_EPS)
                        continue;
                    out.pool.push_back(extractSolution(
                        *model, b, r, z, pairs, nodes, D, N));
                }
            }

            if (out.pool.empty()) {
                model->set(GRB_IntParam_SolutionNumber, 0);
                out.pool.push_back(extractSolution(
                    *model, b, r, z, pairs, nodes, D, N));
            }

            if (!options.quiet) {
                std::cout << std::fixed << std::setprecision(2)
                          << "[ProxyMaster] W_proxy*=" << out.Wproxy
                          << " mm, pool=" << out.pool.size()
                          << ", XPairs(M2b)=" << out.numXPairs
                          << ", cuts=" << cutsApplied << "\n";
            }
        } catch (GRBException& e) {
            out.message = std::string("GRB: ") + e.getMessage();
            if (!options.quiet)
                std::cerr << "[ProxyMaster] " << out.message << "\n";
        } catch (const std::exception& e) {
            out.message = e.what();
        } catch (...) {
            out.message = "unknown exception";
        }
        return out;
    }
};

ProxyMasterMilpSession::ProxyMasterMilpSession(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        ShortcutUsageMode usageMode)
    : impl_(std::make_unique<Impl>(nodes, D, usageMode)) {}

ProxyMasterMilpSession::~ProxyMasterMilpSession() = default;

ProxyMasterMilpResult ProxyMasterMilpSession::solve(
        const std::vector<MasterNogoodCut>& nogoodCuts,
        const ProxyMasterMilpOptions& options) {
    return impl_->solve(nogoodCuts, options);
}

ProxyMasterMilpResult solveProxyMasterMilp(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const std::vector<MasterNogoodCut>& nogoodCuts,
        const ProxyMasterMilpOptions& options) {
    ProxyMasterMilpSession session(nodes, D, options.usageMode);
    return session.solve(nogoodCuts, options);
}
