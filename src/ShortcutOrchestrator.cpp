#include "ShortcutOrchestrator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <utility>

#include "RingLayout.h"
#include "ShortcutGrid.h"

namespace {

bool segmentsNearlyEqual(const Segment& a, const Segment& b, double eps = 1e-6) {
    auto close = [&](double x, double y) { return std::abs(x - y) <= eps; };
    return (close(a.x1, b.x1) && close(a.y1, b.y1) && close(a.x2, b.x2) && close(a.y2, b.y2))
        || (close(a.x1, b.x2) && close(a.y1, b.y2) && close(a.x2, b.x1) && close(a.y2, b.y1));
}

bool isExemptRingSegment(
        const Segment& rs,
        int srcId,
        int destId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing) {
    for (const Segment& s : ShortcutGrid::incidentExemptSegmentsForNode(srcId, tour, routing)) {
        if (segmentsNearlyEqual(rs, s))
            return true;
    }
    for (const Segment& s : ShortcutGrid::incidentExemptSegmentsForNode(destId, tour, routing)) {
        if (segmentsNearlyEqual(rs, s))
            return true;
    }
    return false;
}

bool pathCrossesRing(
        const std::vector<Segment>& pathSegs,
        const std::vector<Segment>& ringSegments,
        int srcId,
        int destId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing) {
    for (const Segment& ps : pathSegs) {
        for (const Segment& rs : ringSegments) {
            if (isExemptRingSegment(rs, srcId, destId, tour, routing))
                continue;
            if (ShortcutRouter::segmentsGeometricallyIntersect(ps, rs))
                return true;
        }
    }
    return false;
}

PlacedShortcut makePlacedShortcut(
        int demandIdx,
        int srcId,
        int destId,
        const ShortcutRouteResult& route) {
    PlacedShortcut placed;
    placed.demandIdx = demandIdx;
    placed.srcId = srcId;
    placed.destId = destId;
    placed.path = route.primary;
    placed.alternatives = route.alternatives;
    placed.totalIL = reportingShortcutIL_dB(route.primary.distance, route.primary.bendCount);
    placed.everCrossed = route.crossedShortcutIdx >= 0;
    placed.crossedShortcutIdx = route.crossedShortcutIdx;
    return placed;
}

Shortcut toRouterShortcut(const PlacedShortcut& placed) {
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

std::optional<ShortcutRouteResult> findFreshRoute(
        int srcId,
        int destId,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        const std::vector<Shortcut>& existingShortcuts,
        double sMin,
        double maxDistance) {
    auto ringSegments = ShortcutGrid::collectRingSegments(routing);
    auto gridResult = ShortcutGrid::buildCached(
        nodes, ringSegments, tour, routing, srcId, destId, sMin);
    if (!gridResult.success)
        return std::nullopt;

    ShortcutRouteResult route = ShortcutRouter::findPaths(
        gridResult.points, srcId, destId, nodes, ringSegments,
        tour, routing, existingShortcuts, maxDistance);
    if (!route.success)
        return std::nullopt;
    return route;
}

bool isWorseDemand(
        double distance,
        int bendCount,
        double bestDistance,
        int bestBendCount) {
    if (distance > bestDistance + MILP_EPS)
        return true;
    if (distance < bestDistance - MILP_EPS)
        return false;
    return bendCount > bestBendCount;
}

bool endpointsMatch(int srcA, int destA, int srcB, int destB) {
    return (srcA == srcB && destA == destB) || (srcA == destB && destA == srcB);
}

const PlacedShortcut* findEndpointMatchedShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts) {
    if (demandIdx < 0 || demandIdx >= (int)D.demands.size())
        return nullptr;
    // Prefer demandIdx association (relaxed endpoints may differ from demand).
    for (const PlacedShortcut& sc : shortcuts) {
        if (sc.demandIdx == demandIdx)
            return &sc;
    }
    const int srcId = D.demands[demandIdx].first;
    const int destId = D.demands[demandIdx].second;
    for (const PlacedShortcut& sc : shortcuts) {
        if (endpointsMatch(srcId, destId, sc.srcId, sc.destId))
            return &sc;
    }
    return nullptr;
}

struct GraphEdge {
    int to = -1;
    double weight = 0.0;
    int shortcutIdx = -1;  // >=0 if this arc is a placed shortcut
    int bendContribution = 0;
};

struct GraphRouteResult {
    double distance = std::numeric_limits<double>::infinity();
    int bendCount = 0;
    int firstShortcutIdx = -1;
};

std::vector<std::vector<GraphEdge>> buildBidirectionalAdj(
        const std::vector<Node>& nodes,
        const std::vector<int>& tour,
        const std::vector<PlacedShortcut>& shortcuts) {
    const int N = (int)nodes.size();
    std::vector<std::vector<GraphEdge>> adj(N);
    const int tourLen = (int)tour.size();
    for (int i = 0; i < tourLen; ++i) {
        const int u = tour[i];
        const int v = tour[(i + 1) % tourLen];
        if (u < 0 || v < 0 || u >= N || v >= N)
            continue;
        const double w = std::abs(nodes[u].x - nodes[v].x)
            + std::abs(nodes[u].y - nodes[v].y);
        const int bends = ringEdgeNeedsBend(
            nodes[u].x, nodes[u].y, nodes[v].x, nodes[v].y) ? 1 : 0;
        adj[u].push_back({v, w, -1, bends});
        adj[v].push_back({u, w, -1, bends});
    }
    for (int si = 0; si < (int)shortcuts.size(); ++si) {
        const PlacedShortcut& sc = shortcuts[si];
        if (sc.srcId < 0 || sc.destId < 0 || sc.srcId >= N || sc.destId >= N)
            continue;
        const double w = sc.path.distance;
        const int bends = sc.path.bendCount;
        adj[sc.srcId].push_back({sc.destId, w, si, bends});
        adj[sc.destId].push_back({sc.srcId, w, si, bends});
    }
    return adj;
}

GraphRouteResult dijkstraDemandRoute(
        int src,
        int dest,
        const std::vector<std::vector<GraphEdge>>& adj) {
    GraphRouteResult out;
    const int N = (int)adj.size();
    if (src < 0 || dest < 0 || src >= N || dest >= N)
        return out;
    if (src == dest) {
        out.distance = 0.0;
        return out;
    }

    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(N, INF);
    std::vector<int> bendAt(N, 0);
    std::vector<int> parent(N, -1);
    std::vector<int> parentEdgeSc(N, -1);
    using State = std::pair<double, int>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    dist[src] = 0.0;
    pq.push({0.0, src});

    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u] + MILP_EPS) continue;
        if (u == dest) break;
        for (const GraphEdge& e : adj[u]) {
            const double nd = d + e.weight;
            const int nb = bendAt[u] + e.bendContribution;
            if (nd + MILP_EPS < dist[e.to]
                || (std::abs(nd - dist[e.to]) <= MILP_EPS && nb < bendAt[e.to])) {
                dist[e.to] = nd;
                bendAt[e.to] = nb;
                parent[e.to] = u;
                parentEdgeSc[e.to] = e.shortcutIdx;
                pq.push({nd, e.to});
            }
        }
    }

    if (!std::isfinite(dist[dest]))
        return out;

    out.distance = dist[dest];
    out.bendCount = bendAt[dest];
    for (int v = dest; v != src && v >= 0; v = parent[v]) {
        if (parentEdgeSc[v] >= 0) {
            out.firstShortcutIdx = parentEdgeSc[v];
            break;
        }
        if (parent[v] < 0) break;
    }
    return out;
}

/// W* = max_d sp(d). Demands with ring(d) <= W*+eps stay on the ring;
/// only those that would otherwise exceed W* may use shortcuts.
struct WcAwareSpBundle {
    std::vector<GraphRouteResult> sp;
    double Wstar = 0.0;
};

WcAwareSpBundle computeAllSpAndWstar(
        const DemandMatrix& D,
        const std::vector<std::vector<GraphEdge>>& adj) {
    WcAwareSpBundle out;
    const int Q = (int)D.demands.size();
    out.sp.resize(Q);
    for (int q = 0; q < Q; ++q) {
        out.sp[q] = dijkstraDemandRoute(
            D.demands[q].first, D.demands[q].second, adj);
        if (std::isfinite(out.sp[q].distance))
            out.Wstar = std::max(out.Wstar, out.sp[q].distance);
    }
    return out;
}

bool demandNeedsShortcutForWc(double ringDist, double Wstar) {
    return ringDist > Wstar + MILP_EPS;
}

EffectiveDemandMetrics metricsFromRingOnly(
        int demandIdx,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<double>& ringDemandIL) {
    EffectiveDemandMetrics m;
    if (demandIdx >= 0 && demandIdx < (int)ringDemandDistance.size())
        m.distance = ringDemandDistance[demandIdx];
    if (demandIdx >= 0 && demandIdx < (int)ringDemandBendCount.size())
        m.bendCount = ringDemandBendCount[demandIdx];
    if (demandIdx >= 0 && demandIdx < (int)ringDemandIL.size())
        m.reportingIL = ringDemandIL[demandIdx];
    else
        m.reportingIL = reportingShortcutIL_dB(m.distance, m.bendCount);
    return m;
}

EffectiveDemandMetrics metricsFromSpRoute(
        int demandIdx,
        int srcId,
        int destId,
        const GraphRouteResult& route,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<double>& ringDemandIL) {
    EffectiveDemandMetrics m;
    m.distance = route.distance;
    m.bendCount = route.bendCount;
    if (route.firstShortcutIdx >= 0
            && route.firstShortcutIdx < (int)shortcuts.size()) {
        const PlacedShortcut& sc = shortcuts[route.firstShortcutIdx];
        // Exclusive shortcut hop: use stored reporting IL when path is just that edge.
        if (endpointsMatch(srcId, destId, sc.srcId, sc.destId)
                && std::abs(route.distance - sc.path.distance) <= MILP_EPS) {
            m.reportingIL = sc.totalIL;
        } else {
            m.reportingIL = reportingShortcutIL_dB(m.distance, m.bendCount);
        }
    } else if (demandIdx >= 0 && demandIdx < (int)ringDemandIL.size()) {
        m.reportingIL = ringDemandIL[demandIdx];
    } else {
        m.reportingIL = reportingShortcutIL_dB(m.distance, m.bendCount);
    }
    return m;
}

/// WC-aware route choice for one demand given precomputed sp(d) and W*.
EffectiveDemandMetrics routeDemandWcAware(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<double>& ringDemandIL,
        const std::vector<PlacedShortcut>& shortcuts,
        const GraphRouteResult& spRoute,
        double Wstar) {
    const double ringDist =
        (demandIdx >= 0 && demandIdx < (int)ringDemandDistance.size())
            ? ringDemandDistance[demandIdx]
            : std::numeric_limits<double>::infinity();

    if (!demandNeedsShortcutForWc(ringDist, Wstar))
        return metricsFromRingOnly(
            demandIdx, ringDemandDistance, ringDemandBendCount, ringDemandIL);

    if (std::isfinite(spRoute.distance) && demandIdx >= 0
            && demandIdx < (int)D.demands.size()) {
        return metricsFromSpRoute(
            demandIdx, D.demands[demandIdx].first, D.demands[demandIdx].second,
            spRoute, shortcuts, ringDemandIL);
    }
    return metricsFromRingOnly(
        demandIdx, ringDemandDistance, ringDemandBendCount, ringDemandIL);
}

EffectiveDemandMetrics metricsFromEndpointMatch(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<double>& ringDemandIL,
        const std::vector<PlacedShortcut>& shortcuts) {
    EffectiveDemandMetrics m = metricsFromRingOnly(
        demandIdx, ringDemandDistance, ringDemandBendCount, ringDemandIL);

    if (const PlacedShortcut* sc = findEndpointMatchedShortcut(demandIdx, D, shortcuts)) {
        m.distance = sc->path.distance;
        m.bendCount = sc->path.bendCount;
        m.reportingIL = sc->totalIL;
    }
    return m;
}

}  // namespace

namespace {
thread_local ShortcutUsageMode g_shortcutUsageMode = ShortcutUsageMode::Shared;
}  // namespace

ShortcutUsageMode getShortcutUsageMode() {
    return g_shortcutUsageMode;
}

void setShortcutUsageMode(ShortcutUsageMode mode) {
    g_shortcutUsageMode = mode;
}

bool shortcutPathCrossesRing(
        const ShortcutPath& path,
        const std::vector<Segment>& ringSegments,
        int srcId,
        int destId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing) {
    return pathCrossesRing(
        pathVerticesToSegments(path.vertices), ringSegments,
        srcId, destId, tour, routing);
}

std::vector<Shortcut> placedShortcutsToRouterFormat(
        const std::vector<PlacedShortcut>& shortcuts) {
    std::vector<Shortcut> out;
    out.reserve(shortcuts.size());
    for (const PlacedShortcut& placed : shortcuts)
        out.push_back(toRouterShortcut(placed));
    return out;
}

const PlacedShortcut* findShortcutForDemand(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts) {
    return findEndpointMatchedShortcut(demandIdx, D, shortcuts);
}

const PlacedShortcut* findShortcutForDemand(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return findShortcutForDemand(
        demandIdx, D, /*ringDemandDistance=*/{}, shortcuts, nodes, tour);
}

const PlacedShortcut* findShortcutForDemand(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    if (tour.empty() || nodes.empty())
        return findEndpointMatchedShortcut(demandIdx, D, shortcuts);
    if (demandIdx < 0 || demandIdx >= (int)D.demands.size())
        return nullptr;
    if (shortcuts.empty())
        return nullptr;

    const auto adj = buildBidirectionalAdj(nodes, tour, shortcuts);
    const WcAwareSpBundle bundle = computeAllSpAndWstar(D, adj);
    const double ringDist =
        (demandIdx < (int)ringDemandDistance.size())
            ? ringDemandDistance[demandIdx]
            : std::numeric_limits<double>::infinity();
    if (!demandNeedsShortcutForWc(ringDist, bundle.Wstar))
        return nullptr;

    const GraphRouteResult& route = bundle.sp[demandIdx];
    if (route.firstShortcutIdx >= 0
            && route.firstShortcutIdx < (int)shortcuts.size())
        return &shortcuts[route.firstShortcutIdx];

    return findEndpointMatchedShortcut(demandIdx, D, shortcuts);
}

bool demandHasShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts) {
    return findEndpointMatchedShortcut(demandIdx, D, shortcuts) != nullptr;
}

bool demandPathUsesShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return demandPathUsesShortcut(
        demandIdx, D, /*ringDemandDistance=*/{}, shortcuts, nodes, tour);
}

bool demandPathUsesShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    if (getShortcutUsageMode() == ShortcutUsageMode::Exclusive)
        return demandHasShortcut(demandIdx, D, shortcuts);
    if (tour.empty() || nodes.empty() || shortcuts.empty())
        return demandHasShortcut(demandIdx, D, shortcuts);
    if (demandIdx < 0 || demandIdx >= (int)D.demands.size())
        return false;

    const auto adj = buildBidirectionalAdj(nodes, tour, shortcuts);
    const WcAwareSpBundle bundle = computeAllSpAndWstar(D, adj);
    const double ringDist =
        (demandIdx < (int)ringDemandDistance.size())
            ? ringDemandDistance[demandIdx]
            : std::numeric_limits<double>::infinity();
    if (!demandNeedsShortcutForWc(ringDist, bundle.Wstar))
        return false;
    return bundle.sp[demandIdx].firstShortcutIdx >= 0;
}

bool demandHasShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return demandPathUsesShortcut(demandIdx, D, shortcuts, nodes, tour);
}

bool demandHasShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return demandPathUsesShortcut(
        demandIdx, D, ringDemandDistance, shortcuts, nodes, tour);
}

EffectiveDemandMetrics effectiveDemandMetrics(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<double>& ringDemandIL,
        const std::vector<PlacedShortcut>& shortcuts) {
    return metricsFromEndpointMatch(
        demandIdx, D, ringDemandDistance, ringDemandBendCount, ringDemandIL, shortcuts);
}

EffectiveDemandMetrics effectiveDemandMetrics(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<double>& ringDemandIL,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    if (tour.empty() || nodes.empty()
            || getShortcutUsageMode() == ShortcutUsageMode::Exclusive) {
        return metricsFromEndpointMatch(
            demandIdx, D, ringDemandDistance, ringDemandBendCount,
            ringDemandIL, shortcuts);
    }
    if (demandIdx < 0 || demandIdx >= (int)D.demands.size())
        return {};

    const auto adj = buildBidirectionalAdj(nodes, tour, shortcuts);
    const WcAwareSpBundle bundle = computeAllSpAndWstar(D, adj);
    return routeDemandWcAware(
        demandIdx, D, ringDemandDistance, ringDemandBendCount, ringDemandIL,
        shortcuts, bundle.sp[demandIdx], bundle.Wstar);
}

int findWorstCaseDemand(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts) {
    return findWorstCaseDemandEffective(
        D, ringDemandDistance, ringDemandBendCount, shortcuts, {});
}

int findWorstCaseDemand(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return findWorstCaseDemandEffective(
        D, ringDemandDistance, ringDemandBendCount, shortcuts, nodes, tour, {});
}

int findWorstCaseDemandEffective(
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::set<int>& excludedDemands) {
    const int numDemands = (int)D.demands.size();
    int worstQ = -1;
    double worstDistance = -1.0;
    int worstBends = -1;
    for (int q = 0; q < numDemands; ++q) {
        if (excludedDemands.count(q)) continue;
        const EffectiveDemandMetrics m = effectiveDemandMetrics(
            q, D, ringDemandDistance, ringDemandBendCount, {}, shortcuts);
        if (worstQ < 0 || isWorseDemand(m.distance, m.bendCount, worstDistance, worstBends)) {
            worstDistance = m.distance;
            worstBends = m.bendCount;
            worstQ = q;
        }
    }
    return worstQ;
}

int findWorstCaseDemandEffective(
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour,
        const std::set<int>& excludedDemands) {
    if (tour.empty() || nodes.empty()
            || getShortcutUsageMode() == ShortcutUsageMode::Exclusive) {
        return findWorstCaseDemandEffective(
            D, ringDemandDistance, ringDemandBendCount, shortcuts, excludedDemands);
    }

    const auto adj = buildBidirectionalAdj(nodes, tour, shortcuts);
    const WcAwareSpBundle bundle = computeAllSpAndWstar(D, adj);

    const int numDemands = (int)D.demands.size();
    int worstQ = -1;
    double worstDistance = -1.0;
    int worstBends = -1;
    for (int q = 0; q < numDemands; ++q) {
        if (excludedDemands.count(q)) continue;
        const EffectiveDemandMetrics m = routeDemandWcAware(
            q, D, ringDemandDistance, ringDemandBendCount, {},
            shortcuts, bundle.sp[q], bundle.Wstar);
        if (worstQ < 0 || isWorseDemand(m.distance, m.bendCount, worstDistance, worstBends)) {
            worstDistance = m.distance;
            worstBends = m.bendCount;
            worstQ = q;
        }
    }
    return worstQ;
}

double effectiveDemandDistance(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts) {
    return effectiveDemandMetrics(
        demandIdx, D, ringDemandDistance, {}, {}, shortcuts).distance;
}

double effectiveDemandDistance(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return effectiveDemandMetrics(
        demandIdx, D, ringDemandDistance, {}, {}, shortcuts, nodes, tour).distance;
}

double computeGlobalW(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts) {
    double W = 0.0;
    for (int q = 0; q < numDemands; ++q) {
        const EffectiveDemandMetrics m = effectiveDemandMetrics(
            q, D, ringDemandDistance, ringDemandBendCount, {}, shortcuts);
        W = std::max(W, m.distance);
    }
    return W;
}

double computeGlobalW(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    if (tour.empty() || nodes.empty()
            || getShortcutUsageMode() == ShortcutUsageMode::Exclusive) {
        return computeGlobalW(
            numDemands, D, ringDemandDistance, ringDemandBendCount, shortcuts);
    }

    const auto adj = buildBidirectionalAdj(nodes, tour, shortcuts);
    const WcAwareSpBundle bundle = computeAllSpAndWstar(D, adj);

    double W = 0.0;
    for (int q = 0; q < numDemands; ++q) {
        const EffectiveDemandMetrics m = routeDemandWcAware(
            q, D, ringDemandDistance, ringDemandBendCount, {},
            shortcuts, bundle.sp[q], bundle.Wstar);
        W = std::max(W, m.distance);
    }
    return W;
}

bool strictlyImprovesGlobalW(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        double oldGlobalW) {
    return computeGlobalW(
        numDemands, D, ringDemandDistance, ringDemandBendCount, shortcuts)
        + MILP_EPS < oldGlobalW;
}

bool strictlyImprovesGlobalW(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const std::vector<PlacedShortcut>& shortcuts,
        double oldGlobalW,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour) {
    return computeGlobalW(
        numDemands, D, ringDemandDistance, ringDemandBendCount, shortcuts, nodes, tour)
        + MILP_EPS < oldGlobalW;
}

double worstRingRoutedDistance(
        int numDemands,
        const DemandMatrix& D,
        const std::vector<double>& ringDemandDistance,
        const std::vector<PlacedShortcut>& shortcuts) {
    double worst = -1.0;
    for (int q = 0; q < numDemands; ++q) {
        if (demandHasShortcut(q, D, shortcuts))
            continue;
        if (q < (int)ringDemandDistance.size())
            worst = std::max(worst, ringDemandDistance[q]);
    }
    return worst;
}

void applyShortcutCrossingSideEffects(
        std::vector<PlacedShortcut>& shortcuts,
        const PlacedShortcut& placed) {
    if (placed.crossedShortcutIdx < 0)
        return;
    if (placed.crossedShortcutIdx < (int)shortcuts.size())
        shortcuts[placed.crossedShortcutIdx].everCrossed = true;
}

std::optional<PlacedShortcut> tryPlaceShortcut(
        int demandIdx,
        int srcId,
        int destId,
        const std::vector<Node>& nodes,
        const std::vector<int>& currentTour,
        const std::vector<EdgeOption>& currentRouting,
        const std::vector<Shortcut>& existingShortcutsForRouter,
        const std::set<int>& usedNodes,
        double sMin,
        double maxDistance) {
    if (usedNodes.count(srcId) || usedNodes.count(destId))
        return std::nullopt;

    auto ringSegments = ShortcutGrid::collectRingSegments(currentRouting);
    auto gridResult = ShortcutGrid::buildCached(
        nodes, ringSegments, currentTour, currentRouting, srcId, destId, sMin);
    if (!gridResult.success)
        return std::nullopt;

    ShortcutRouteResult route = ShortcutRouter::findPaths(
        gridResult.points, srcId, destId, nodes, ringSegments,
        currentTour, currentRouting, existingShortcutsForRouter, maxDistance);
    if (!route.success)
        return std::nullopt;

    return makePlacedShortcut(demandIdx, srcId, destId, route);
}

std::optional<PlacedShortcut> tryPlaceShortcut(
        int demandIdx,
        const DemandMatrix& D,
        const std::vector<Node>& nodes,
        const std::vector<int>& currentTour,
        const std::vector<EdgeOption>& currentRouting,
        const std::vector<Shortcut>& existingShortcutsForRouter,
        const std::set<int>& usedNodes,
        double sMin,
        double maxDistance,
        double currentRingDistance) {
    const int srcId = D.demands[demandIdx].first;
    const int destId = D.demands[demandIdx].second;

    if (usedNodes.count(srcId) || usedNodes.count(destId)) {
        std::cerr << "[tryPlaceShortcut] demand " << demandIdx << " ("
                  << srcId + 1 << "->" << destId + 1
                  << "): endpoint node already used by another shortcut\n";
        return std::nullopt;
    }

    auto placed = tryPlaceShortcut(
        demandIdx, srcId, destId, nodes, currentTour, currentRouting,
        existingShortcutsForRouter, usedNodes, sMin, maxDistance);
    if (!placed) {
        std::cerr << "[tryPlaceShortcut] demand " << demandIdx << " ("
                  << srcId + 1 << "->" << destId + 1
                  << "): no valid geometric route\n";
        return std::nullopt;
    }

    // Strict raw-distance improvement required; bend count must not override equal distance.
    if (placed->path.distance >= currentRingDistance - MILP_EPS) {
        std::cerr << "[tryPlaceShortcut] demand " << demandIdx << " ("
                  << srcId + 1 << "->" << destId + 1
                  << "): route distance " << placed->path.distance
                  << " mm not shorter than ring " << currentRingDistance << " mm\n";
        return std::nullopt;
    }

    return placed;
}

std::vector<double> allRingDemandDistances(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layout) {
    MILPSolver solver(nodes, D);
    MILPSolveResult tmp = layout;
    solver.applyShortestArcDemandRouting(tmp, {});
    return tmp.demandDistance;
}

ShortcutRevalidationResult revalidateAllShortcuts(
        std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::vector<int>& newTour,
        const std::vector<EdgeOption>& newRouting,
        double sMin,
        const std::vector<double>& ringDemandDistance,
        double wcCap) {
    ShortcutRevalidationResult result;
    const auto ringSegments = ShortcutGrid::collectRingSegments(newRouting);

    auto ownRingDistance = [&](int demandIdx) -> double {
        if (demandIdx < 0 || demandIdx >= (int)ringDemandDistance.size())
            return std::numeric_limits<double>::infinity();
        return ringDemandDistance[demandIdx];
    };

    // Accept a path if it is strictly shorter than the demand's ring distance
    // and strictly shorter than the new WC cap (user criterion).
    auto withinCap = [&](double dist, double ringDist) -> bool {
        const double maxAllowed = std::min(ringDist, wcCap);
        return dist + MILP_EPS < maxAllowed;
    };

    auto pathOk = [&](const ShortcutPath& path, int srcId, int destId, double ringDist) {
        if (shortcutPathCrossesRing(path, ringSegments, srcId, destId, newTour, newRouting))
            return false;
        return withinCap(path.distance, ringDist);
    };

    std::vector<PlacedShortcut> kept;
    kept.reserve(shortcuts.size());

    for (size_t i = 0; i < shortcuts.size(); ++i) {
        PlacedShortcut placed = shortcuts[i];
        const double ringDist = ownRingDistance(placed.demandIdx);

        if (pathOk(placed.path, placed.srcId, placed.destId, ringDist)) {
            kept.push_back(std::move(placed));
            ++result.kept;
            continue;
        }

        bool repaired = false;

        for (const ShortcutPath& alt : placed.alternatives) {
            if (!pathOk(alt, placed.srcId, placed.destId, ringDist))
                continue;
            placed.path = alt;
            placed.totalIL = reportingShortcutIL_dB(alt.distance, alt.bendCount);
            repaired = true;
            break;
        }

        if (!repaired) {
            std::vector<Shortcut> others;
            others.reserve(shortcuts.size());
            for (size_t j = 0; j < shortcuts.size(); ++j) {
                if (j == i) continue;
                // Prefer already-accepted peers; fall back to original list.
                others.push_back(toRouterShortcut(shortcuts[j]));
            }
            for (const PlacedShortcut& k : kept)
                others.push_back(toRouterShortcut(k));

            const double maxAllowed = std::min(ringDist, wcCap);
            auto fresh = findFreshRoute(
                placed.srcId, placed.destId, nodes, newTour, newRouting,
                others, sMin, maxAllowed);

            if (fresh
                && pathOk(fresh->primary, placed.srcId, placed.destId, ringDist)) {
                placed.path = fresh->primary;
                placed.alternatives = fresh->alternatives;
                placed.totalIL = reportingShortcutIL_dB(
                    fresh->primary.distance, fresh->primary.bendCount);
                placed.everCrossed = fresh->crossedShortcutIdx >= 0;
                repaired = true;
            }
        }

        if (repaired) {
            kept.push_back(std::move(placed));
            ++result.repaired;
        } else {
            std::cerr << "[revalidateAllShortcuts] Dropping shortcut for demand "
                      << placed.demandIdx
                      << " (no route shorter than min(ring=" << ringDist
                      << ", WC=" << wcCap << ") mm after ring change).\n";
            ++result.dropped;
        }
    }

    shortcuts = std::move(kept);
    return result;
}
