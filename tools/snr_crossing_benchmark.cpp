// SNR-oriented shortcut crossing audit — Proxy Master vs W_base.
//
// Study baseline:
//   - Ring: Method B (WC) with 90s cap + warm-start from A; fallback A
//   - Shared W_base: Method D joint shortcuts
//   - Private W_base: WC-greedy Method A+Shortcuts (no Method-D MILP)
//
// Runs:
//   - N=6 paper case (d) DSP, seeds 1..20
//   - N=8 synth density=1, seeds 1..20
//   - N=12 curated 10 fast cases from (a)/(b)/(c)
//   - N=14/16 synth density=1, seeds 1..10
//
// Combined wall budget per (seed x mode): N6=240s, N8=360s, N>=12=600s.
//
// Usage:
//   ./snr_crossing_benchmark              — study batch (N6+N8+N12+N14+N16)
//   ./snr_crossing_benchmark d|n8|n12|n14|n16
//   ./snr_crossing_benchmark a 10         — single case/seed smoke test
//   ./snr_crossing_benchmark --private-shortcuts d
//       — private start→end-only shortcut flow; writes *_private.csv

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "DemandMatrix.h"
#include "DemandMatrices.h"
#include "InstanceSetup.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "PhysicalConstants.h"
#include "ProxyMasterLoop.h"
#include "ShortcutGrid.h"
#include "ShortcutMethodD.h"
#include "ShortcutMethods.h"
#include "ShortcutOrchestrator.h"
#include "ShortcutRouter.h"
#include "WavelengthLoadBalance.h"

namespace {

using Clock = std::chrono::steady_clock;

struct BenchCase {
    char id = 'a';
    std::string slug;
    std::string title;
    int N = 0;
    DemandMatrix D;
};

struct CrossPairDetail {
    int scI = -1;
    int scJ = -1;
    int loadI = 0;
    int loadJ = 0;
    int signalsAffected = 0;  // loadI + loadJ
};

struct ShortcutDetail {
    int scIdx = -1;
    int srcId = -1;
    int destId = -1;
    int load = 0;
    int crossesWith = 0;  // #other shortcuts this one geometrically intersects
};

struct LayoutSnrStats {
    int numShortcuts = 0;
    int signalsOnShortcuts = 0;
    int crossPairs = 0;
    std::vector<int> loadPerSc;
    std::vector<ShortcutDetail> shortcuts;
    std::vector<CrossPairDetail> pairs;
};

struct Job {
    char caseId = 'a';
    int seed = 0;
};

struct SeedRow {
    char caseId = '?';
    std::string slug;
    int seed = 0;
    int N = 0;
    double wStar = std::numeric_limits<double>::quiet_NaN();
    double gap = std::numeric_limits<double>::quiet_NaN();
    bool proven = false;
    std::string stopReason;
    double proxySec = 0.0;
    LayoutSnrStats proxySnr;
    bool proxySnrOk = false;
    double wBase = std::numeric_limits<double>::quiet_NaN();
    double baseSec = 0.0;
    bool baseOk = false;
    std::string baseRingMethod;       // B | A | B_timeout
    std::string baseShortcutMethod;   // D | A_greedy
    LayoutSnrStats baseSnr;
    bool baseSnrOk = false;
};

constexpr double kMethodBRingTimeLimitSec = 90.0;

double wallBudgetForN(int N) {
    if (N <= 6) return 240.0;
    if (N <= 8) return 360.0;
    return 600.0;
}

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r')) --b;
    return s.substr(a, b - a);
}

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == ',') {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(trim(cur));
    return out;
}

bool loadCase(
        const std::string& root,
        char id,
        const std::string& slug,
        const std::string& title,
        BenchCase& out) {
    out.id = id;
    out.slug = slug;
    out.title = title;
    const std::string nodesPath = root + "/benchmarks/parsed/" + slug + "_nodes.csv";
    const std::string demPath = root + "/benchmarks/parsed/" + slug + "_demands.csv";

    std::ifstream nf(nodesPath);
    if (!nf) {
        std::cerr << "missing " << nodesPath << "\n";
        return false;
    }
    std::string line;
    std::getline(nf, line);
    int nNodes = 0;
    while (std::getline(nf, line)) {
        if (trim(line).empty()) continue;
        ++nNodes;
    }
    out.N = nNodes;

    std::ifstream df(demPath);
    if (!df) {
        std::cerr << "missing " << demPath << "\n";
        return false;
    }
    std::getline(df, line);
    out.D.demands.clear();
    while (std::getline(df, line)) {
        if (trim(line).empty()) continue;
        auto cols = splitCsv(line);
        if (cols.size() < 2) continue;
        out.D.add(std::stoi(cols[0]), std::stoi(cols[1]));
    }
    return out.N >= 3 && !out.D.demands.empty();
}

std::string mapStop(const ProxyMasterLoopResult& r) {
    if (r.stoppedFine) return "fine";
    if (r.stoppedCoarse) return "coarse";
    if (r.stoppedBudget) return "wall-time-limit";
    if (r.stopReason.find("max rounds") != std::string::npos) return "round-limit";
    if (!r.stopReason.empty()) return r.stopReason;
    return "unknown";
}

// Local WC-aware SP that records ALL shortcut indices on the path (not only first).
struct Edge {
    int to = -1;
    double w = 0.0;
    int sc = -1;
};

std::vector<int> shortcutIndicesOnSp(
        int src,
        int dest,
        const std::vector<std::vector<Edge>>& adj) {
    const int N = (int)adj.size();
    std::vector<int> empty;
    if (src < 0 || dest < 0 || src >= N || dest >= N) return empty;
    if (src == dest) return empty;

    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(N, INF);
    std::vector<int> parent(N, -1);
    std::vector<int> parentSc(N, -1);
    using State = std::pair<double, int>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    dist[src] = 0.0;
    pq.push({0.0, src});
    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u] + MILP_EPS) continue;
        if (u == dest) break;
        for (const Edge& e : adj[u]) {
            const double nd = d + e.w;
            if (nd + MILP_EPS < dist[e.to]) {
                dist[e.to] = nd;
                parent[e.to] = u;
                parentSc[e.to] = e.sc;
                pq.push({nd, e.to});
            }
        }
    }
    if (!std::isfinite(dist[dest])) return empty;

    std::vector<int> used;
    for (int v = dest; v != src && v >= 0; v = parent[v]) {
        if (parentSc[v] >= 0) used.push_back(parentSc[v]);
        if (parent[v] < 0) break;
    }
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    return used;
}

LayoutSnrStats computeSnrStats(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const ShortcutMethodResult& res,
        ShortcutUsageMode usageMode) {
    LayoutSnrStats st;
    if (!res.layout.success) return st;
    const int S = (int)res.shortcuts.size();
    st.numShortcuts = S;
    st.loadPerSc.assign(S, 0);
    if (S == 0) return st;

    const std::vector<double> ringDist = allRingDemandDistances(nodes, D, res.layout);
    const int N = (int)nodes.size();
    const bool privateOnly = (usageMode == ShortcutUsageMode::Exclusive);

    auto buildRingAdj = [&]() {
        std::vector<std::vector<Edge>> adj(N);
        const auto& tour = res.layout.tour;
        for (int i = 0; i < (int)tour.size(); ++i) {
            const int u = tour[i];
            const int v = tour[(i + 1) % (int)tour.size()];
            if (u < 0 || v < 0 || u >= N || v >= N) continue;
            const double w = std::abs(nodes[u].x - nodes[v].x)
                + std::abs(nodes[u].y - nodes[v].y);
            adj[u].push_back({v, w, -1});
            adj[v].push_back({u, w, -1});
        }
        return adj;
    };

    auto spDistOnly = [&](const std::vector<std::vector<Edge>>& adj, int s, int t) {
        const double INF = std::numeric_limits<double>::infinity();
        std::vector<double> dist(N, INF);
        using State = std::pair<double, int>;
        std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
        dist[s] = 0.0;
        pq.push({0.0, s});
        while (!pq.empty()) {
            const auto [d, u] = pq.top();
            pq.pop();
            if (d > dist[u] + MILP_EPS) continue;
            if (u == t) break;
            for (const Edge& e : adj[u]) {
                const double nd = d + e.w;
                if (nd + MILP_EPS < dist[e.to]) {
                    dist[e.to] = nd;
                    pq.push({nd, e.to});
                }
            }
        }
        return dist[t];
    };

    std::vector<double> spDist((int)D.demands.size(), 0.0);
    double Wstar = 0.0;
    std::set<int> demandsOnSc;

    if (privateOnly) {
        // Per-demand graph: ring + only owner start↔end shortcuts.
        for (int q = 0; q < (int)D.demands.size(); ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            auto adj = buildRingAdj();
            for (int si = 0; si < S; ++si) {
                const auto& sc = res.shortcuts[si];
                if (sc.srcId < 0 || sc.destId < 0 || sc.srcId >= N || sc.destId >= N)
                    continue;
                if (!undirectedEndpointsMatch(sc.srcId, sc.destId, s, t)) continue;
                adj[sc.srcId].push_back({sc.destId, sc.path.distance, si});
                adj[sc.destId].push_back({sc.srcId, sc.path.distance, si});
            }
            spDist[q] = spDistOnly(adj, s, t);
            if (std::isfinite(spDist[q])) Wstar = std::max(Wstar, spDist[q]);
        }
        for (int q = 0; q < (int)D.demands.size(); ++q) {
            const double ring = (q < (int)ringDist.size())
                ? ringDist[q]
                : std::numeric_limits<double>::infinity();
            if (std::isfinite(ring) && ring <= Wstar + MILP_EPS) continue;

            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            auto adj = buildRingAdj();
            for (int si = 0; si < S; ++si) {
                const auto& sc = res.shortcuts[si];
                if (sc.srcId < 0 || sc.destId < 0 || sc.srcId >= N || sc.destId >= N)
                    continue;
                if (!undirectedEndpointsMatch(sc.srcId, sc.destId, s, t)) continue;
                adj[sc.srcId].push_back({sc.destId, sc.path.distance, si});
                adj[sc.destId].push_back({sc.srcId, sc.path.distance, si});
            }
            const auto used = shortcutIndicesOnSp(s, t, adj);
            if (used.empty()) continue;
            demandsOnSc.insert(q);
            for (int si : used) {
                if (si >= 0 && si < S) ++st.loadPerSc[si];
            }
        }
    } else {
        auto adj = buildRingAdj();
        for (int si = 0; si < S; ++si) {
            const auto& sc = res.shortcuts[si];
            if (sc.srcId < 0 || sc.destId < 0 || sc.srcId >= N || sc.destId >= N)
                continue;
            adj[sc.srcId].push_back({sc.destId, sc.path.distance, si});
            adj[sc.destId].push_back({sc.srcId, sc.path.distance, si});
        }

        // W* from free SP over ring+shortcuts (same WC gate as demandPathUsesShortcut).
        for (int q = 0; q < (int)D.demands.size(); ++q) {
            const int s = D.demands[q].first;
            const int t = D.demands[q].second;
            spDist[q] = spDistOnly(adj, s, t);
            if (std::isfinite(spDist[q])) Wstar = std::max(Wstar, spDist[q]);
        }

        for (int q = 0; q < (int)D.demands.size(); ++q) {
            const double ring = (q < (int)ringDist.size())
                ? ringDist[q]
                : std::numeric_limits<double>::infinity();
            // WC gate: stay on ring if ring(d) <= W*+eps
            if (std::isfinite(ring) && ring <= Wstar + MILP_EPS) continue;

            const auto used = shortcutIndicesOnSp(
                D.demands[q].first, D.demands[q].second, adj);
            if (used.empty()) continue;
            demandsOnSc.insert(q);
            for (int si : used) {
                if (si >= 0 && si < S) ++st.loadPerSc[si];
            }
        }
    }
    st.signalsOnShortcuts = (int)demandsOnSc.size();

    const std::vector<Shortcut> routerSc = placedShortcutsToRouterFormat(res.shortcuts);
    std::vector<int> crossCount(S, 0);
    for (int i = 0; i < S; ++i) {
        for (int j = i + 1; j < S; ++j) {
            if (!ShortcutRouter::shortcutsIntersect(routerSc[i], routerSc[j]))
                continue;
            CrossPairDetail p;
            p.scI = i;
            p.scJ = j;
            p.loadI = st.loadPerSc[i];
            p.loadJ = st.loadPerSc[j];
            p.signalsAffected = p.loadI + p.loadJ;
            st.pairs.push_back(p);
            ++crossCount[i];
            ++crossCount[j];
        }
    }
    st.crossPairs = (int)st.pairs.size();

    for (int i = 0; i < S; ++i) {
        ShortcutDetail d;
        d.scIdx = i;
        d.srcId = res.shortcuts[i].srcId;
        d.destId = res.shortcuts[i].destId;
        d.load = st.loadPerSc[i];
        d.crossesWith = crossCount[i];
        st.shortcuts.push_back(d);
    }
    return st;
}

ShortcutMethodResult runBaselineFull(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        double sMin,
        bool* ok,
        ShortcutUsageMode usageMode,
        double wallRemainingSec,
        std::string* ringMethodOut,
        std::string* shortcutMethodOut) {
    ShortcutMethodResult empty;
    *ok = false;
    if (ringMethodOut) *ringMethodOut = "";
    if (shortcutMethodOut) *shortcutMethodOut = "";
    if (wallRemainingSec <= 0.5) {
        std::cerr << "[baseline] abort: wallRemainingSec=" << wallRemainingSec << "\n";
        return empty;
    }

    const auto t0 = Clock::now();
    MILPSolver solverA(nodes, D);
    MILPSolveResult layoutA = solverA.solve(false, "", {}, -1.0, nullptr, -1.0, true);
    if (!layoutA.success) {
        std::cerr << "[baseline] Method A ring failed status=" << layoutA.status
                  << " — trying Method B cold\n";
    }

    const double elapsedA = std::chrono::duration<double>(Clock::now() - t0).count();
    const double remainAfterA = wallRemainingSec - elapsedA;
    if (remainAfterA <= 0.5) {
        if (!layoutA.success) return empty;
        if (ringMethodOut) *ringMethodOut = "A";
        ShortcutMethodOptions opt;
        opt.quiet = true;
        opt.skipExports = true;
        opt.usageMode = usageMode;
        opt.wallTimeLimitSec = 0.5;
        ShortcutMethodResult res;
        if (usageMode == ShortcutUsageMode::Exclusive) {
            if (shortcutMethodOut) *shortcutMethodOut = "A_greedy";
            res = runMethodAWithShortcuts(nodes, D, layoutA, sMin, opt);
        } else {
            if (shortcutMethodOut) *shortcutMethodOut = "D";
            res = runMethodDJointShortcuts(nodes, D, layoutA, sMin, opt);
        }
        *ok = std::isfinite(res.globalW);
        return res;
    }

    const double bCap = std::min(kMethodBRingTimeLimitSec, remainAfterA);
    MILPSolver solverB(nodes, D);
    const MILPSolveResult* warm = layoutA.success ? &layoutA : nullptr;
    MILPSolveResult layoutB =
        solverB.solve(true, "", {}, -1.0, warm, bCap, true);

    const MILPSolveResult* ring = nullptr;
    std::string ringMethod;
    if (layoutB.success) {
        ring = &layoutB;
        ringMethod = (layoutB.status == GRB_TIME_LIMIT) ? "B_timeout" : "B";
    } else if (layoutA.success) {
        ring = &layoutA;
        ringMethod = "A";
    } else {
        std::cerr << "[baseline] Method B also failed status=" << layoutB.status << "\n";
        return empty;
    }
    if (ringMethodOut) *ringMethodOut = ringMethod;

    const double elapsedRing = std::chrono::duration<double>(Clock::now() - t0).count();
    const double remainForSc = wallRemainingSec - elapsedRing;

    ShortcutMethodOptions opt;
    opt.quiet = true;
    opt.skipExports = true;
    opt.usageMode = usageMode;
    if (remainForSc > 0.0)
        opt.wallTimeLimitSec = remainForSc;

    ShortcutMethodResult res;
    if (usageMode == ShortcutUsageMode::Exclusive) {
        if (shortcutMethodOut) *shortcutMethodOut = "A_greedy";
        res = runMethodAWithShortcuts(nodes, D, *ring, sMin, opt);
    } else {
        if (shortcutMethodOut) *shortcutMethodOut = "D";
        res = runMethodDJointShortcuts(nodes, D, *ring, sMin, opt);
    }
    *ok = std::isfinite(res.globalW);
    return res;
}

SeedRow runJob(
        const BenchCase& bc,
        int seed,
        double wallSec,
        ShortcutUsageMode usageMode) {
    SeedRow row;
    row.caseId = bc.id;
    row.slug = bc.slug;
    row.seed = seed;
    row.N = bc.N;

    constexpr double spacing = 5.0 * ShortcutGrid::DEFAULT_S_MIN;
    const std::vector<Node> nodes =
        generateNodes(bc.N, static_cast<unsigned>(seed), spacing);
    const double sMin = ShortcutGrid::DEFAULT_S_MIN;

    // Reserve wall time for W_base so proxy overruns cannot starve the baseline.
    const double baseReserve = (bc.N >= 12)
        ? std::min(180.0, wallSec * 0.35)
        : std::min(90.0, wallSec * 0.35);
    const double proxyBudget = std::max(30.0, wallSec - baseReserve);

    ProxyMasterLoopOptions opt;
    opt.quiet = true;
    opt.poolSolutions = 10;
    opt.maxRounds = 20;
    opt.wallTimeLimitSec = proxyBudget;
    opt.traceRounds = false;
    opt.incrementalMaster = true;
    opt.sMin = sMin;
    opt.usageMode = usageMode;

    const auto t0 = Clock::now();
    const ProxyMasterLoopResult loop = runProxyMasterLoop(nodes, bc.D, opt);
    row.proxySec = std::chrono::duration<double>(Clock::now() - t0).count();
    row.wStar = loop.Wstar;
    row.gap = loop.gap;
    row.proven = loop.provenOptimal;
    row.stopReason = mapStop(loop);

    if (std::isfinite(loop.Wstar) && loop.best.layout.success) {
        ScopedShortcutUsageMode scope(usageMode);
        (void)runWavelengthLoadBalance(nodes, bc.D, loop.best, loop.Wstar);
        row.proxySnr = computeSnrStats(nodes, bc.D, loop.best, usageMode);
        row.proxySnrOk = true;
    }

    // Give baseline the reserved slice plus any unused proxy budget.
    const double unusedProxy = std::max(0.0, proxyBudget - row.proxySec);
    const double baseBudget = baseReserve + unusedProxy;

    const auto t1 = Clock::now();
    bool baseOk = false;
    ShortcutMethodResult baseRes = runBaselineFull(
        nodes, bc.D, sMin, &baseOk, usageMode, baseBudget,
        &row.baseRingMethod, &row.baseShortcutMethod);

    // If standalone A/B ring MILP is infeasible, fall back to the proxy incumbent
    // ring and still build W_base shortcuts (tagged proxy_fallback).
    if (!baseOk && loop.best.layout.success) {
        std::cerr << "[baseline] A/B infeasible — proxy_fallback ring\n";
        const double elapsed = std::chrono::duration<double>(Clock::now() - t1).count();
        const double remain = std::max(1.0, baseBudget - elapsed);
        ShortcutMethodOptions scOpt;
        scOpt.quiet = true;
        scOpt.skipExports = true;
        scOpt.usageMode = usageMode;
        scOpt.wallTimeLimitSec = remain;
        if (usageMode == ShortcutUsageMode::Exclusive) {
            row.baseShortcutMethod = "A_greedy";
            baseRes = runMethodAWithShortcuts(
                nodes, bc.D, loop.best.layout, sMin, scOpt);
        } else {
            row.baseShortcutMethod = "D";
            baseRes = runMethodDJointShortcuts(
                nodes, bc.D, loop.best.layout, sMin, scOpt);
        }
        baseOk = std::isfinite(baseRes.globalW);
        if (baseOk) row.baseRingMethod = "proxy_fallback";
    }

    row.baseSec = std::chrono::duration<double>(Clock::now() - t1).count();
    row.baseOk = baseOk;
    if (baseOk) {
        row.wBase = baseRes.globalW;
        ScopedShortcutUsageMode scope(usageMode);
        (void)runWavelengthLoadBalance(nodes, bc.D, baseRes, baseRes.globalW);
        row.baseSnr = computeSnrStats(nodes, bc.D, baseRes, usageMode);
        row.baseSnrOk = true;
    }
    return row;
}

// Curated N=12: first 10 of overnight fast proven set.
const std::vector<Job> kCuratedN12 = {
    {'c', 13}, {'b', 10}, {'c', 10}, {'c', 6}, {'b', 4},
    {'c', 19}, {'b', 6}, {'c', 5}, {'c', 7}, {'b', 8},
};

void meanPairStats(const LayoutSnrStats& snr, double* meanLoadCross, double* meanAff) {
    *meanLoadCross = 0.0;
    *meanAff = 0.0;
    if (snr.pairs.empty()) return;
    double sumLoad = 0.0;
    double sumAff = 0.0;
    for (const auto& p : snr.pairs) {
        sumLoad += 0.5 * (p.loadI + p.loadJ);
        sumAff += p.signalsAffected;
    }
    *meanLoadCross = sumLoad / snr.pairs.size();
    *meanAff = sumAff / snr.pairs.size();
}

void writeSeedHeader(std::ostream& os) {
    os << "case,slug,seed,N,W_star,W_base,gap,proven,stop_reason,"
          "proxy_time_s,base_time_s,base_ring_method,base_shortcut_method,"
          "proxy_num_shortcuts,proxy_signals_on_shortcuts,proxy_cross_pairs,"
          "proxy_mean_load_on_crossing_scs,proxy_mean_signals_affected_per_pair,"
          "base_num_shortcuts,base_signals_on_shortcuts,base_cross_pairs,"
          "base_mean_load_on_crossing_scs,base_mean_signals_affected_per_pair\n";
}

void writeSeedRow(std::ostream& os, const SeedRow& r) {
    auto num = [](double x) -> std::string {
        if (!std::isfinite(x)) return "";
        std::ostringstream s;
        s << std::fixed << std::setprecision(6) << x;
        return s.str();
    };
    auto writeSnrCols = [&](bool ok, const LayoutSnrStats& snr) {
        if (!ok) {
            os << ",,,,";
            return;
        }
        double meanLoad = 0.0, meanAff = 0.0;
        meanPairStats(snr, &meanLoad, &meanAff);
        os << snr.numShortcuts << "," << snr.signalsOnShortcuts << ","
           << snr.crossPairs << ","
           << std::fixed << std::setprecision(4) << meanLoad << "," << meanAff;
    };

    os << r.caseId << "," << r.slug << "," << r.seed << "," << r.N << ","
       << num(r.wStar) << "," << (r.baseOk ? num(r.wBase) : "") << ","
       << num(r.gap) << ","
       << (r.proven ? "Y" : "N") << "," << r.stopReason << ","
       << num(r.proxySec) << "," << num(r.baseSec) << ","
       << r.baseRingMethod << "," << r.baseShortcutMethod << ",";
    writeSnrCols(r.proxySnrOk, r.proxySnr);
    os << ",";
    writeSnrCols(r.baseSnrOk, r.baseSnr);
    os << "\n";
}

void writePairHeader(std::ostream& os) {
    os << "source,case,slug,seed,N,sc_i,sc_j,load_i,load_j,signals_affected\n";
}

void writeScHeader(std::ostream& os) {
    os << "source,case,slug,seed,N,sc_idx,src,dest,load,crosses_with\n";
}

void appendSnrDetails(
        std::ostream& pairs,
        std::ostream& scs,
        const std::string& source,
        const SeedRow& r,
        bool ok,
        const LayoutSnrStats& snr) {
    if (!ok) return;
    for (const auto& p : snr.pairs) {
        pairs << source << "," << r.caseId << "," << r.slug << "," << r.seed << ","
              << r.N << "," << p.scI << "," << p.scJ << ","
              << p.loadI << "," << p.loadJ << "," << p.signalsAffected << "\n";
    }
    for (const auto& s : snr.shortcuts) {
        scs << source << "," << r.caseId << "," << r.slug << "," << r.seed << ","
            << r.N << "," << s.scIdx << "," << s.srcId << ","
            << s.destId << "," << s.load << "," << s.crossesWith << "\n";
    }
}

void printAggregateSource(
        const std::string& label,
        const std::string& source,
        const std::vector<SeedRow>& rows,
        bool useProxy) {
    int n = (int)rows.size();
    int withSnr = 0;
    int withCross = 0;
    double sumPairs = 0.0;
    double sumAffAllPairs = 0.0;
    int nPairs = 0;
    double sumPairsGivenCross = 0.0;

    for (const auto& r : rows) {
        const bool ok = useProxy ? r.proxySnrOk : r.baseSnrOk;
        const LayoutSnrStats& snr = useProxy ? r.proxySnr : r.baseSnr;
        if (!ok) continue;
        ++withSnr;
        sumPairs += snr.crossPairs;
        if (snr.crossPairs > 0) {
            ++withCross;
            sumPairsGivenCross += snr.crossPairs;
        }
        for (const auto& p : snr.pairs) {
            sumAffAllPairs += p.signalsAffected;
            ++nPairs;
        }
    }

    std::cout << "\n=== SNR aggregate [" << source << "]: " << label << " ===\n"
              << "instances: " << n << "  with_layout: " << withSnr
              << "  with_≥1_crossing: " << withCross << "\n";
    if (withSnr > 0) {
        std::cout << std::fixed << std::setprecision(3)
                  << "mean cross_pairs / instance: "
                  << (sumPairs / withSnr) << "\n";
    }
    if (withCross > 0) {
        std::cout << std::fixed << std::setprecision(3)
                  << "mean cross_pairs | instance has crossing: "
                  << (sumPairsGivenCross / withCross) << "\n";
    }
    if (nPairs > 0) {
        std::cout << std::fixed << std::setprecision(3)
                  << "crossing pairs total: " << nPairs << "\n"
                  << "mean signals_affected / crossing pair (load_i+load_j): "
                  << (sumAffAllPairs / nPairs) << "\n";
    } else {
        std::cout << "no shortcut crossings observed\n";
    }
}

void printAggregate(const std::string& label, const std::vector<SeedRow>& rows) {
    printAggregateSource(label, "proxy", rows, true);
    printAggregateSource(label, "base", rows, false);
}

}  // namespace

int main(int argc, char* argv[]) {
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);

    const std::string root =
        ".";
    const std::string outDir = root + "/benchmarks/results";

    std::map<char, BenchCase> cases;
    struct Meta { char id; const char* slug; const char* title; };
    const Meta metas[] = {
        {'a', "a_mpeg4", "mpeg4"},
        {'b', "b_vopd", "vopd"},
        {'c', "c_mwd", "mwd"},
        {'d', "d_dsp", "dsp"},
    };
    for (const auto& m : metas) {
        BenchCase bc;
        if (!loadCase(root, m.id, m.slug, m.title, bc)) {
            std::cerr << "failed to load case " << m.id << "\n";
            return 1;
        }
        cases[m.id] = std::move(bc);
        std::cout << "loaded (" << m.id << ") " << m.slug
                  << " N=" << cases[m.id].N
                  << " Q=" << cases[m.id].D.demands.size() << "\n";
    }

    // Synth N=8/14/16 (density=1).
    auto addSynth = [&](char id, int N, const char* slug, const char* title) {
        BenchCase synth;
        synth.id = id;
        synth.slug = slug;
        synth.title = title;
        synth.N = N;
        if (!initDemandMatrix(N, 1, synth.D)) {
            std::cerr << "initDemandMatrix(" << N << ",1) failed\n";
            std::exit(1);
        }
        cases[id] = std::move(synth);
        std::cout << "loaded (" << id << ") " << slug << " N=" << N
                  << " Q=" << cases[id].D.demands.size() << "\n";
    };
    addSynth('s', 8, "synth_n8_d1", "synth N=8 dens=1");
    addSynth('t', 14, "synth_n14_d1", "synth N=14 dens=1");
    addSynth('u', 16, "synth_n16_d1", "synth N=16 dens=1");

    enum Mode { All, OnlyD, OnlyN8, OnlyN12, OnlyN14, OnlyN16, Single };
    Mode mode = All;
    char singleCase = 0;
    int singleSeed = -1;
    bool privateShortcuts = false;

    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--private-shortcuts" || a == "--private") {
            privateShortcuts = true;
            continue;
        }
        args.push_back(a);
    }

    if (!args.empty()) {
        const std::string& a1 = args[0];
        // Prefer "case seed" single-job form when a seed is provided.
        if (a1.size() == 1
                && ((a1[0] >= 'a' && a1[0] <= 'd')
                    || a1[0] == 's' || a1[0] == 't' || a1[0] == 'u')
                && args.size() >= 2) {
            mode = Single;
            singleCase = a1[0];
            singleSeed = std::stoi(args[1]);
        } else if (a1 == "d") {
            mode = OnlyD;
        } else if (a1 == "n8") {
            mode = OnlyN8;
        } else if (a1 == "n12") {
            mode = OnlyN12;
        } else if (a1 == "n14") {
            mode = OnlyN14;
        } else if (a1 == "n16") {
            mode = OnlyN16;
        } else {
            std::cerr << "Usage: ./snr_crossing_benchmark [--private-shortcuts] "
                         "[d|n8|n12|n14|n16|<a|b|c|d|s|t|u> SEED]\n";
            return 1;
        }
    }

    const ShortcutUsageMode usageMode = privateShortcuts
        ? ShortcutUsageMode::Exclusive
        : ShortcutUsageMode::Shared;

    std::vector<Job> jobs;
    if (mode == All || mode == OnlyD) {
        for (int s = 1; s <= 20; ++s) jobs.push_back({'d', s});
    }
    if (mode == All || mode == OnlyN8) {
        for (int s = 1; s <= 20; ++s) jobs.push_back({'s', s});
    }
    if (mode == All || mode == OnlyN12) {
        jobs.insert(jobs.end(), kCuratedN12.begin(), kCuratedN12.end());
    }
    if (mode == All || mode == OnlyN14) {
        for (int s = 1; s <= 10; ++s) jobs.push_back({'t', s});
    }
    if (mode == All || mode == OnlyN16) {
        for (int s = 1; s <= 10; ++s) jobs.push_back({'u', s});
    }
    if (mode == Single) jobs.push_back({singleCase, singleSeed});

    // Study outputs live under snr_study/ so prior Method-D-private CSVs stay intact.
    const std::string studyDir = outDir + "/snr_study";
    {
        const std::string mkdirCmd = "mkdir -p \"" + studyDir + "\"";
        if (std::system(mkdirCmd.c_str()) != 0) {
            std::cerr << "failed to create " << studyDir << "\n";
            return 1;
        }
    }
    const std::string priv = privateShortcuts ? "private" : "shared";
    const std::string seedCsv = studyDir + "/seeds_" + priv + ".csv";
    const std::string pairCsv = studyDir + "/pairs_" + priv + ".csv";
    const std::string scCsv = studyDir + "/shortcuts_" + priv + ".csv";

    std::set<std::string> doneKeys;
    {
        std::ifstream in(seedCsv);
        std::string line;
        if (in && std::getline(in, line)
                && line.find("base_ring_method") != std::string::npos
                && line.find("proxy_cross_pairs") != std::string::npos) {
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                auto cols = splitCsv(line);
                // case,slug,seed,...,W_base at col 5 (0-based index 5)
                if (cols.size() < 13 || cols[0].empty()) continue;
                // Skip incomplete rows (empty W_base) so resume re-runs them.
                if (cols[5].empty()) continue;
                doneKeys.insert(cols[0] + ":" + cols[2]);
            }
        } else if (!line.empty()) {
            const std::string bak = seedCsv + ".old.bak";
            std::rename(seedCsv.c_str(), bak.c_str());
            std::rename(pairCsv.c_str(), (pairCsv + ".old.bak").c_str());
            std::rename(scCsv.c_str(), (scCsv + ".old.bak").c_str());
            std::cout << "archived incompatible CSVs → *.old.bak\n";
            doneKeys.clear();
        }
    }
    const bool resume = !doneKeys.empty();
    if (!resume) {
        std::ofstream ofs(seedCsv); writeSeedHeader(ofs);
        std::ofstream ofp(pairCsv); writePairHeader(ofp);
        std::ofstream ofc(scCsv); writeScHeader(ofc);
    } else {
        // Drop incomplete seed rows (empty W_base) so a re-run appends cleanly.
        {
            std::ifstream in(seedCsv);
            std::string line;
            std::vector<std::string> keep;
            if (in && std::getline(in, line)) {
                keep.push_back(line);
                while (std::getline(in, line)) {
                    if (line.empty()) continue;
                    auto cols = splitCsv(line);
                    if (cols.size() < 6 || cols[5].empty()) continue;
                    keep.push_back(line);
                }
            }
            std::ofstream ofs(seedCsv);
            for (const auto& l : keep) ofs << l << "\n";
        }
    }

    std::cout << "\n=== snr_crossing_benchmark study (proxy + W_base"
              << (privateShortcuts
                      ? ", PRIVATE / W_base=A_greedy"
                      : ", SHARED / W_base=MethodD")
              << "; ring=B@90s) ===\n"
              << "jobs=" << jobs.size()
              << (resume ? (" resume_skip=" + std::to_string(doneKeys.size())) : "")
              << "\nout seeds: " << seedCsv
              << "\nout pairs: " << pairCsv
              << "\nout shortcuts: " << scCsv << "\n\n";

    std::vector<SeedRow> allRows;
    std::map<int, std::vector<SeedRow>> rowsByN;

    for (const Job& job : jobs) {
        const std::string key = std::string(1, job.caseId) + ":" + std::to_string(job.seed);
        if (doneKeys.count(key)) {
            std::cout << "  (" << job.caseId << ") seed " << job.seed << " ... SKIP (resume)\n";
            continue;
        }
        const BenchCase& bc = cases.at(job.caseId);
        const double wall = wallBudgetForN(bc.N);
        std::cout << "  (" << job.caseId << ") seed " << job.seed
                  << " N=" << bc.N << " wall=" << wall << "s ..." << std::flush;
        SeedRow row;
        try {
            row = runJob(bc, job.seed, wall, usageMode);
        } catch (const std::exception& e) {
            std::cerr << " EXCEPTION: " << e.what() << "\n";
            row.caseId = job.caseId;
            row.slug = bc.slug;
            row.seed = job.seed;
            row.N = bc.N;
            row.stopReason = std::string("exception:") + e.what();
        }
        std::cout << " W*=" << row.wStar
                  << " W_base=" << (row.baseOk ? std::to_string(row.wBase) : "FAIL")
                  << " ring=" << row.baseRingMethod
                  << " sc=" << row.baseShortcutMethod
                  << " proven=" << (row.proven ? "Y" : "N")
                  << " tP=" << std::fixed << std::setprecision(1) << row.proxySec << "s"
                  << " tB=" << row.baseSec << "s";
        auto dumpSnr = [](const char* tag, bool ok, const LayoutSnrStats& snr) {
            if (!ok) return;
            std::cout << " " << tag << "(sc=" << snr.numShortcuts
                      << ",sig=" << snr.signalsOnShortcuts
                      << ",xp=" << snr.crossPairs;
            if (!snr.pairs.empty()) {
                double meanAff = 0.0;
                for (const auto& p : snr.pairs) meanAff += p.signalsAffected;
                meanAff /= snr.pairs.size();
                std::cout << ",aff=" << std::setprecision(2) << meanAff;
            }
            std::cout << ")";
        };
        dumpSnr("P", row.proxySnrOk, row.proxySnr);
        dumpSnr("B", row.baseSnrOk, row.baseSnr);
        std::cout << "\n";

        {
            std::ofstream ofs(seedCsv, std::ios::app);
            writeSeedRow(ofs, row);
            std::ofstream ofp(pairCsv, std::ios::app);
            std::ofstream ofc(scCsv, std::ios::app);
            appendSnrDetails(ofp, ofc, "proxy", row, row.proxySnrOk, row.proxySnr);
            appendSnrDetails(ofp, ofc, "base", row, row.baseSnrOk, row.baseSnr);
        }

        allRows.push_back(row);
        rowsByN[row.N].push_back(row);
    }

    for (const auto& kv : rowsByN) {
        printAggregate("N=" + std::to_string(kv.first), kv.second);
    }
    if (allRows.size() > 1) printAggregate("all jobs this run", allRows);

    std::cout << "\nDone.\n";
    return 0;
}
