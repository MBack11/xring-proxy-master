#include "WavelengthLoadBalance.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "PhysicalConstants.h"
#include "RingLayout.h"
#include "ShortcutRouter.h"
#include "ShortcutTypes.h"

namespace {

Shortcut placedToRouter(const PlacedShortcut& placed) {
    Shortcut sc;
    sc.from = placed.srcId;
    sc.to = placed.destId;
    sc.demandIdx = placed.demandIdx;
    sc.approx_length = placed.path.distance;
    sc.bend_count = placed.path.bendCount;
    sc.everCrossed = placed.everCrossed;
    sc.path = pathVerticesToSegments(placed.path.vertices);
    return sc;
}

double manhattan(const Node& a, const Node& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

struct TourRing {
    std::vector<int> tour;
    std::vector<int> pos;
    int nTour = 0;

    TourRing(const std::vector<int>& t, int nNodes)
        : tour(t), pos(nNodes, -1), nTour((int)t.size()) {
        for (int i = 0; i < nTour; ++i)
            if (tour[i] >= 0 && tour[i] < nNodes)
                pos[tour[i]] = i;
    }

    WlbRoute arc(
            int s,
            int d,
            int step,
            const std::vector<Node>& nodes,
            const char* label) const {
        WlbRoute p;
        p.label = label;
        if (s == d) {
            p.length = 0.0;
            return p;
        }
        if (s < 0 || d < 0 || s >= (int)pos.size() || d >= (int)pos.size())
            return p;
        if (pos[s] < 0 || pos[d] < 0) return p;

        int k = pos[s];
        const int kd = pos[d];
        int bends = 0;
        double dist = 0.0;
        while (k != kd) {
            const int nk = (k + step + nTour) % nTour;
            const int u = tour[k];
            const int v = tour[nk];
            p.ringArcs.push_back({u, v});
            dist += manhattan(nodes[u], nodes[v]);
            if (ringEdgeNeedsBend(
                    nodes[u].x, nodes[u].y, nodes[v].x, nodes[v].y))
                ++bends;
            k = nk;
        }
        p.length = dist;
        p.bendCount = bends;
        return p;
    }

    WlbRoute cwArc(int s, int d, const std::vector<Node>& nodes) const {
        return arc(s, d, +1, nodes, "CW");
    }
    WlbRoute ccwArc(int s, int d, const std::vector<Node>& nodes) const {
        return arc(s, d, -1, nodes, "CCW");
    }
};

WlbRoute concatRoutes(const WlbRoute& a, const WlbRoute& b) {
    WlbRoute out;
    out.ringArcs = a.ringArcs;
    out.ringArcs.insert(out.ringArcs.end(), b.ringArcs.begin(), b.ringArcs.end());
    out.shortcutHops = a.shortcutHops;
    out.shortcutHops.insert(
        out.shortcutHops.end(), b.shortcutHops.begin(), b.shortcutHops.end());
    out.length = a.length + b.length;
    out.bendCount = a.bendCount + b.bendCount;
    return out;
}

std::string makeLabel(const WlbRoute& r) {
    if (r.shortcutHops.empty()) {
        // Infer CW vs CCW from first arc if possible — label already set.
        return r.label.empty() ? "ring" : r.label;
    }
    std::ostringstream oss;
    if (!r.ringArcs.empty() && !r.label.empty() && r.label != "SC")
        oss << r.label << "+";
    oss << "SC";
    for (size_t i = 0; i < r.shortcutHops.size(); ++i) {
        if (i) oss << ",";
        oss << "#" << r.shortcutHops[i].shortcutIdx;
    }
    return oss.str();
}

using ScDirKey = std::pair<int, int>;  // (shortcutIdx, dirSign) dirSign = ±1

int shortcutDirSign(
        const std::vector<PlacedShortcut>& shortcuts,
        int si,
        int from,
        int to) {
    if (si < 0 || si >= (int)shortcuts.size()) return 0;
    const auto& sc = shortcuts[si];
    if (from == sc.srcId && to == sc.destId) return +1;
    if (from == sc.destId && to == sc.srcId) return -1;
    return 0;
}

ScDirKey hopKey(
        const std::vector<PlacedShortcut>& shortcuts,
        const WlbShortcutHop& h) {
    return {h.shortcutIdx, shortcutDirSign(shortcuts, h.shortcutIdx, h.from, h.to)};
}

bool sameConflictSignature(const WlbRoute& a, const WlbRoute& b) {
    if (a.ringArcs.size() != b.ringArcs.size()) return false;
    if (a.shortcutHops.size() != b.shortcutHops.size()) return false;
    std::set<std::pair<int, int>> ea(a.ringArcs.begin(), a.ringArcs.end());
    std::set<std::pair<int, int>> eb(b.ringArcs.begin(), b.ringArcs.end());
    if (ea != eb) return false;
    std::set<std::tuple<int, int, int>> sa;
    std::set<std::tuple<int, int, int>> sb;
    for (const auto& h : a.shortcutHops)
        sa.insert({h.shortcutIdx, h.from, h.to});
    for (const auto& h : b.shortcutHops)
        sb.insert({h.shortcutIdx, h.from, h.to});
    return sa == sb;
}

void addUnique(
        std::vector<WlbRoute>& paths,
        WlbRoute cand,
        double Wstar) {
    if (!std::isfinite(cand.length) || cand.length > Wstar + MILP_EPS)
        return;
    cand.label = makeLabel(cand);
    for (const auto& p : paths)
        if (sameConflictSignature(p, cand))
            return;
    paths.push_back(std::move(cand));
}

WlbRoute dijkstraSp(
        int src,
        int dest,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour,
        const std::vector<PlacedShortcut>& shortcuts) {
    WlbRoute info;
    const int N = (int)nodes.size();
    if (src < 0 || dest < 0 || src >= N || dest >= N) return info;
    if (src == dest) {
        info.length = 0.0;
        info.label = "SP";
        return info;
    }

    struct GEdge {
        int to = -1;
        double weight = 0.0;
        int shortcutIdx = -1;
        int bends = 0;
    };
    std::vector<std::vector<GEdge>> adj(N);
    const int tourLen = (int)tour.size();
    for (int i = 0; i < tourLen; ++i) {
        const int u = tour[i];
        const int v = tour[(i + 1) % tourLen];
        if (u < 0 || v < 0 || u >= N || v >= N) continue;
        const double w = manhattan(nodes[u], nodes[v]);
        const int b = ringEdgeNeedsBend(
            nodes[u].x, nodes[u].y, nodes[v].x, nodes[v].y) ? 1 : 0;
        adj[u].push_back({v, w, -1, b});
        adj[v].push_back({u, w, -1, b});
    }
    for (int si = 0; si < (int)shortcuts.size(); ++si) {
        const auto& sc = shortcuts[si];
        if (sc.srcId < 0 || sc.destId < 0) continue;
        adj[sc.srcId].push_back({sc.destId, sc.path.distance, si, sc.path.bendCount});
        adj[sc.destId].push_back({sc.srcId, sc.path.distance, si, sc.path.bendCount});
    }

    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(N, INF);
    std::vector<int> parent(N, -1);
    std::vector<int> parentSc(N, -1);
    std::vector<int> bendAt(N, 0);
    using State = std::pair<double, int>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    dist[src] = 0.0;
    pq.push({0.0, src});
    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u] + MILP_EPS) continue;
        if (u == dest) break;
        for (const GEdge& e : adj[u]) {
            const double nd = d + e.weight;
            const int nb = bendAt[u] + e.bends;
            if (nd + MILP_EPS < dist[e.to]
                    || (std::abs(nd - dist[e.to]) <= MILP_EPS && nb < bendAt[e.to])) {
                dist[e.to] = nd;
                bendAt[e.to] = nb;
                parent[e.to] = u;
                parentSc[e.to] = e.shortcutIdx;
                pq.push({nd, e.to});
            }
        }
    }
    if (!std::isfinite(dist[dest])) return info;

    info.length = dist[dest];
    info.bendCount = bendAt[dest];
    info.label = "SP";
    for (int v = dest; v != src && v >= 0; v = parent[v]) {
        const int u = parent[v];
        if (u < 0) break;
        if (parentSc[v] >= 0)
            info.shortcutHops.push_back({parentSc[v], u, v});
        else
            info.ringArcs.push_back({u, v});
    }
    std::reverse(info.ringArcs.begin(), info.ringArcs.end());
    std::reverse(info.shortcutHops.begin(), info.shortcutHops.end());
    return info;
}

std::vector<WlbRoute> enumerateAdmissible(
        int demandIdx,
        int src,
        int dest,
        const MILPSolveResult& layout,
        const std::vector<Node>& nodes,
        const TourRing& ring,
        const std::vector<PlacedShortcut>& shortcuts,
        double Wstar) {
    std::vector<WlbRoute> paths;
    addUnique(paths, ring.cwArc(src, dest, nodes), Wstar);
    addUnique(paths, ring.ccwArc(src, dest, nodes), Wstar);

    // Layout's original Method-D / ring routing (always admissible if W* correct).
    if (demandIdx >= 0 && demandIdx < (int)layout.demandFlowEdges.size()) {
        WlbRoute layoutRoute;
        layoutRoute.label = "layout";
        layoutRoute.ringArcs = layout.demandFlowEdges[demandIdx];
        if (demandIdx < (int)layout.demandDistance.size())
            layoutRoute.length = layout.demandDistance[demandIdx];
        if (demandIdx < (int)layout.demandBendCount.size())
            layoutRoute.bendCount = layout.demandBendCount[demandIdx];
        // If layout used shortcuts, distance already includes them; arcs may be
        // ring-only in demandFlowEdges. Prefer SP for shortcut-using layout.
        addUnique(paths, layoutRoute, Wstar);
    }

    addUnique(paths, dijkstraSp(src, dest, nodes, layout.tour, shortcuts), Wstar);

    for (int si = 0; si < (int)shortcuts.size(); ++si) {
        const auto& sc = shortcuts[si];
        const int u = sc.srcId;
        const int v = sc.destId;
        if (u < 0 || v < 0) continue;

        const int ends[2][2] = {{u, v}, {v, u}};
        for (const auto& endsPair : ends) {
            const int entry = endsPair[0];
            const int exit = endsPair[1];
            WlbRoute scHop;
            scHop.label = "SC";
            scHop.length = sc.path.distance;
            scHop.bendCount = sc.path.bendCount;
            scHop.shortcutHops = {{si, entry, exit}};

            for (int modeIn = 0; modeIn < 2; ++modeIn) {
                WlbRoute toEntry = (modeIn == 0)
                    ? ring.cwArc(src, entry, nodes)
                    : ring.ccwArc(src, entry, nodes);
                if (!std::isfinite(toEntry.length) && src != entry) continue;
                for (int modeOut = 0; modeOut < 2; ++modeOut) {
                    WlbRoute fromExit = (modeOut == 0)
                        ? ring.cwArc(exit, dest, nodes)
                        : ring.ccwArc(exit, dest, nodes);
                    if (!std::isfinite(fromExit.length) && exit != dest) continue;
                    WlbRoute full = concatRoutes(concatRoutes(toEntry, scHop), fromExit);
                    full.label = (modeIn == 0 ? "CW" : "CCW");
                    if (!scHop.shortcutHops.empty())
                        full.label += "+SC#" + std::to_string(si) + "+"
                            + (modeOut == 0 ? "CW" : "CCW");
                    addUnique(paths, full, Wstar);
                }
            }
        }
    }
    return paths;
}

// Segment keys for load counters.
using RingKey = std::pair<int, int>;  // directed

int maxDirShortcutLoad(
        int si,
        const std::map<ScDirKey, int>& scDirLoad) {
    int mf = 0, mr = 0;
    auto itF = scDirLoad.find({si, +1});
    auto itR = scDirLoad.find({si, -1});
    if (itF != scDirLoad.end()) mf = itF->second;
    if (itR != scDirLoad.end()) mr = itR->second;
    return std::max(mf, mr);
}

int conflictLoadOfRoute(
        const WlbRoute& r,
        const std::map<RingKey, int>& ringLoad,
        const std::map<ScDirKey, int>& scDirLoad,
        const std::vector<std::vector<int>>& scCross,
        const std::vector<PlacedShortcut>& shortcuts) {
    int mx = 0;
    for (const auto& e : r.ringArcs) {
        auto it = ringLoad.find(e);
        const int load = (it == ringLoad.end()) ? 0 : it->second;
        mx = std::max(mx, load);
    }
    for (const auto& h : r.shortcutHops) {
        if (h.shortcutIdx < 0) continue;
        const ScDirKey key = hopKey(shortcuts, h);
        auto it = scDirLoad.find(key);
        int cl = (it == scDirLoad.end()) ? 0 : it->second;
        // Crossing clique: this direction on si vs the denser direction on cj
        // (opposite dirs on cj do not form a clique with each other).
        if (h.shortcutIdx < (int)scCross.size()) {
            for (int cj : scCross[h.shortcutIdx])
                cl += maxDirShortcutLoad(cj, scDirLoad);
        }
        mx = std::max(mx, cl);
    }
    return mx;
}

void applyLoad(
        const WlbRoute& r,
        std::map<RingKey, int>& ringLoad,
        std::map<ScDirKey, int>& scDirLoad,
        const std::vector<PlacedShortcut>& shortcuts) {
    for (const auto& e : r.ringArcs)
        ringLoad[e] += 1;
    for (const auto& h : r.shortcutHops) {
        if (h.shortcutIdx < 0) continue;
        scDirLoad[hopKey(shortcuts, h)] += 1;
    }
}

bool routesConflict(
        const WlbRoute& a,
        const WlbRoute& b,
        const std::vector<Shortcut>& routerSc,
        const std::vector<PlacedShortcut>& shortcuts) {
    // Same directed ring arc → conflict (CW vs CCW already separated by direction).
    std::set<RingKey> ea(a.ringArcs.begin(), a.ringArcs.end());
    for (const auto& e : b.ringArcs)
        if (ea.count(e)) return true;

    // Same shortcut, same direction → conflict; opposite directions do NOT.
    std::set<ScDirKey> sa;
    for (const auto& h : a.shortcutHops)
        sa.insert(hopKey(shortcuts, h));
    for (const auto& h : b.shortcutHops)
        if (sa.count(hopKey(shortcuts, h))) return true;

    // Different shortcuts that cross → conflict regardless of direction.
    for (const auto& ha : a.shortcutHops) {
        for (const auto& hb : b.shortcutHops) {
            if (ha.shortcutIdx == hb.shortcutIdx) continue;
            if (ha.shortcutIdx < 0 || hb.shortcutIdx < 0
                    || ha.shortcutIdx >= (int)routerSc.size()
                    || hb.shortcutIdx >= (int)routerSc.size())
                continue;
            if (ShortcutRouter::shortcutsIntersect(
                    routerSc[ha.shortcutIdx], routerSc[hb.shortcutIdx]))
                return true;
        }
    }
    return false;
}

/// DSATUR: repeatedly color the uncolored vertex with highest saturation
/// (ties: highest degree, then lowest index).
std::vector<int> dsaturColor(
        int n,
        const std::vector<std::vector<char>>& conflict) {
    std::vector<int> color(n, 0);  // 1-indexed when assigned
    std::vector<char> colored(n, 0);
    std::vector<int> degree(n, 0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && conflict[i][j])
                ++degree[i];

    int remaining = n;
    while (remaining > 0) {
        int best = -1;
        int bestSat = -1;
        int bestDeg = -1;
        for (int i = 0; i < n; ++i) {
            if (colored[i]) continue;
            std::set<int> neighColors;
            for (int j = 0; j < n; ++j)
                if (conflict[i][j] && colored[j])
                    neighColors.insert(color[j]);
            const int sat = (int)neighColors.size();
            if (best < 0
                    || sat > bestSat
                    || (sat == bestSat && degree[i] > bestDeg)
                    || (sat == bestSat && degree[i] == bestDeg && i < best)) {
                best = i;
                bestSat = sat;
                bestDeg = degree[i];
            }
        }
        std::set<int> used;
        for (int j = 0; j < n; ++j)
            if (conflict[best][j] && colored[j])
                used.insert(color[j]);
        int c = 1;
        while (used.count(c)) ++c;
        color[best] = c;
        colored[best] = 1;
        --remaining;
    }
    return color;
}

std::string segmentsString(const WlbRoute& r) {
    std::ostringstream oss;
    bool first = true;
    for (const auto& [u, v] : r.ringArcs) {
        if (!first) oss << " ";
        first = false;
        oss << u << "->" << v;
    }
    for (const auto& h : r.shortcutHops) {
        if (!first) oss << " ";
        first = false;
        oss << "SC#" << h.shortcutIdx << ":" << h.from << "->" << h.to;
    }
    return oss.str();
}

}  // namespace

WavelengthLoadBalanceResult runWavelengthLoadBalance(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const ShortcutMethodResult& loopResult,
        double Wstar) {
    WavelengthLoadBalanceResult out;
    out.Wstar = Wstar;

    if (!loopResult.layout.success || loopResult.layout.tour.empty()) {
        out.message = "invalid loop result layout";
        return out;
    }
    if (!std::isfinite(Wstar)) {
        out.message = "non-finite W*";
        return out;
    }

    const MILPSolveResult& layout = loopResult.layout;
    const auto& shortcuts = loopResult.shortcuts;
    const int Q = (int)D.demands.size();
    const int N = (int)nodes.size();
    TourRing ring(layout.tour, N);

    std::vector<Shortcut> routerSc;
    routerSc.reserve(shortcuts.size());
    for (const auto& ps : shortcuts)
        routerSc.push_back(placedToRouter(ps));

    const int S = (int)shortcuts.size();
    std::vector<std::vector<int>> scCross(S);
    for (int i = 0; i < S; ++i) {
        for (int j = i + 1; j < S; ++j) {
            if (ShortcutRouter::shortcutsIntersect(routerSc[i], routerSc[j])) {
                scCross[i].push_back(j);
                scCross[j].push_back(i);
            }
        }
    }

    // Step 1: Adm(q)
    std::vector<std::vector<WlbRoute>> adm(Q);
    for (int q = 0; q < Q; ++q) {
        const int s = D.demands[q].first;
        const int t = D.demands[q].second;
        adm[q] = enumerateAdmissible(
            q, s, t, layout, nodes, ring, shortcuts, Wstar);
        if (adm[q].empty()) {
            out.message = "empty Adm(q) for demand " + std::to_string(q)
                + " — W* likely inconsistent with layout";
            out.allRoutesRespectWstar = false;
            return out;
        }
    }

    // Step 2–3: fixed demands first, then ascending |Adm|
    std::map<RingKey, int> ringLoad;
    std::map<ScDirKey, int> scDirLoad;
    std::vector<WlbRoute> chosen(Q);
    std::vector<char> assigned(Q, 0);

    auto assignOne = [&](int q, const WlbRoute& r) {
        chosen[q] = r;
        assigned[q] = 1;
        applyLoad(r, ringLoad, scDirLoad, shortcuts);
        if (r.length > Wstar + MILP_EPS)
            out.allRoutesRespectWstar = false;
    };

    for (int q = 0; q < Q; ++q) {
        if ((int)adm[q].size() == 1)
            assignOne(q, adm[q][0]);
    }

    std::vector<int> open;
    for (int q = 0; q < Q; ++q)
        if (!assigned[q]) open.push_back(q);
    std::sort(open.begin(), open.end(), [&](int a, int b) {
        if (adm[a].size() != adm[b].size())
            return adm[a].size() < adm[b].size();
        return a < b;
    });

    for (int q : open) {
        int bestIdx = 0;
        int bestCost = conflictLoadOfRoute(
            adm[q][0], ringLoad, scDirLoad, scCross, shortcuts);
        int bestBends = adm[q][0].bendCount;
        for (int i = 1; i < (int)adm[q].size(); ++i) {
            const int c = conflictLoadOfRoute(
                adm[q][i], ringLoad, scDirLoad, scCross, shortcuts);
            const int b = adm[q][i].bendCount;
            if (c < bestCost
                    || (c == bestCost && b < bestBends)
                    || (c == bestCost && b == bestBends
                        && adm[q][i].label < adm[q][bestIdx].label)) {
                bestIdx = i;
                bestCost = c;
                bestBends = b;
            }
        }
        assignOne(q, adm[q][bestIdx]);
    }

    // Final conflictLoad lower bound (max over directed segments after assignments).
    int maxConflictLoad = 0;
    for (const auto& [e, load] : ringLoad)
        maxConflictLoad = std::max(maxConflictLoad, load);
    for (const auto& [key, load] : scDirLoad) {
        const int si = key.first;
        int cl = load;
        if (si >= 0 && si < S) {
            for (int cj : scCross[si])
                cl += maxDirShortcutLoad(cj, scDirLoad);
        }
        maxConflictLoad = std::max(maxConflictLoad, cl);
    }

    // Step 4: conflict graph + DSATUR
    std::vector<std::vector<char>> conflict(Q, std::vector<char>(Q, 0));
    for (int i = 0; i < Q; ++i) {
        for (int j = i + 1; j < Q; ++j) {
            if (routesConflict(chosen[i], chosen[j], routerSc, shortcuts)) {
                conflict[i][j] = 1;
                conflict[j][i] = 1;
            }
        }
    }
    const std::vector<int> colors = dsaturColor(Q, conflict);
    int lambda = 0;
    for (int c : colors) lambda = std::max(lambda, c);

    out.assignments.resize(Q);
    for (int q = 0; q < Q; ++q) {
        WlbDemandAssignment a;
        a.demandIdx = q;
        a.src = D.demands[q].first;
        a.dest = D.demands[q].second;
        a.route = chosen[q];
        a.wavelength = colors[q];
        a.admCount = (int)adm[q].size();
        a.wasFixed = (a.admCount == 1);
        out.assignments[q] = std::move(a);
    }

    out.lambdaNeeded = lambda;
    out.conflictLoadLowerBound = maxConflictLoad;
    out.success = true;
    out.message = "ok";
    return out;
}

void printWavelengthLoadBalanceReport(
        const WavelengthLoadBalanceResult& res,
        const DemandMatrix& D) {
    (void)D;
    std::cout << std::fixed;
    std::cout.precision(3);
    if (!res.success) {
        std::cout << "[WLB] FAILED: " << res.message << "\n";
        return;
    }

    std::cout << "\ndemand | route | segments | λ | |Adm| | fixed?\n";
    std::cout << "-------|-------|----------|---|------|--------\n";
    for (const auto& a : res.assignments) {
        std::cout << "  " << a.demandIdx
                  << " (" << a.src << "→" << a.dest << ")"
                  << " | " << a.route.label
                  << " | " << segmentsString(a.route)
                  << " | " << a.wavelength
                  << " | " << a.admCount
                  << " | " << (a.wasFixed ? "Y" : "N")
                  << "\n";
    }

    std::cout << "\n--- WLB summary ---\n";
    std::cout << "W*                  : " << res.Wstar << "\n";
    std::cout << "λ_needed (DSATUR)   : " << res.lambdaNeeded << "\n";
    std::cout << "conflictLoad LB     : " << res.conflictLoadLowerBound << "\n";
    std::cout << "all routes ≤ W*     : "
              << (res.allRoutesRespectWstar ? "yes" : "NO") << "\n";
}
