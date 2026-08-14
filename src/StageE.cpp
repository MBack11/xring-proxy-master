#include "StageE.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#include "PhysicalConstants.h"
#include "RingGeometry.h"

namespace {

struct Rect {
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
};

Rect boundingRect(const std::vector<Node>& nodes, int i, int j) {
    const double xi = nodes[i].x, xj = nodes[j].x;
    const double yi = nodes[i].y, yj = nodes[j].y;
    return {std::min(xi, xj), std::max(xi, xj),
            std::min(yi, yj), std::max(yi, yj)};
}

bool onBoundaryThroughEndpoints(
        const RingSegment& seg,
        const Rect& R,
        const std::vector<Node>& nodes,
        int i,
        int j) {
    constexpr double eps = 1e-9;
    const double xi = nodes[i].x, yi = nodes[i].y;
    const double xj = nodes[j].x, yj = nodes[j].y;

    if (seg.isHorizontal) {
        if (std::abs(seg.fixed - R.ymin) > eps && std::abs(seg.fixed - R.ymax) > eps)
            return false;
        const bool atI = seg.lo <= xi + eps && seg.hi >= xi - eps
            && std::abs(seg.fixed - yi) < eps;
        const bool atJ = seg.lo <= xj + eps && seg.hi >= xj - eps
            && std::abs(seg.fixed - yj) < eps;
        return atI || atJ;
    }
    if (std::abs(seg.fixed - R.xmin) > eps && std::abs(seg.fixed - R.xmax) > eps)
        return false;
    const bool atI = seg.lo <= yi + eps && seg.hi >= yi - eps
        && std::abs(seg.fixed - xi) < eps;
    const bool atJ = seg.lo <= yj + eps && seg.hi >= yj - eps
        && std::abs(seg.fixed - xj) < eps;
    return atI || atJ;
}

bool segmentIrrelevant(
        const RingSegment& seg,
        const Rect& R,
        const std::vector<Node>& nodes,
        int i,
        int j) {
    constexpr double eps = 1e-9;
    if (onBoundaryThroughEndpoints(seg, R, nodes, i, j))
        return true;

    if (seg.isHorizontal) {
        if (seg.fixed < R.ymin - eps || seg.fixed > R.ymax + eps) return true;
        if (seg.hi < R.xmin - eps || seg.lo > R.xmax + eps) return true;
        return false;
    }
    if (seg.fixed < R.xmin - eps || seg.fixed > R.xmax + eps) return true;
    if (seg.hi < R.ymin - eps || seg.lo > R.ymax + eps) return true;
    return false;
}

bool segmentFullySevering(const RingSegment& seg, const Rect& R) {
    constexpr double eps = 1e-9;
    if (seg.isHorizontal) {
        return seg.fixed > R.ymin + eps && seg.fixed < R.ymax - eps
            && seg.lo <= R.xmin + eps && seg.hi >= R.xmax - eps;
    }
    return seg.fixed > R.xmin + eps && seg.fixed < R.xmax - eps
        && seg.lo <= R.ymin + eps && seg.hi >= R.ymax - eps;
}

double overhangPenalty(const RingSegment& seg, const Rect& R, double sMin) {
    if (seg.isHorizontal) {
        const double left = 2.0 * ((R.xmin - seg.lo) + sMin);
        const double right = 2.0 * ((seg.hi - R.xmax) + sMin);
        return std::max(0.0, std::min(left, right));
    }
    const double bottom = 2.0 * ((R.ymin - seg.lo) + sMin);
    const double top = 2.0 * ((seg.hi - R.ymax) + sMin);
    return std::max(0.0, std::min(bottom, top));
}

double hatShortcutLength(
        const std::vector<Node>& nodes,
        const MILPSolveResult& layout,
        int i,
        int j,
        double sMin) {
    const Rect R = boundingRect(nodes, i, j);
    const double delta = std::abs(nodes[i].x - nodes[j].x)
        + std::abs(nodes[i].y - nodes[j].y);

    RingGeometry geom(nodes);
    double maxPen = 0.0;

    for (const auto& te : layout.tourEdges) {
        for (int opt = 0; opt <= 1; ++opt) {
            for (const RingSegment& seg : geom.getSegments(te.from, te.to, opt)) {
                if (segmentIrrelevant(seg, R, nodes, i, j))
                    continue;
                if (!segmentFullySevering(seg, R))
                    continue;
                maxPen = std::max(maxPen, overhangPenalty(seg, R, sMin));
            }
        }
    }
    return delta + maxPen;
}

double shortestPathWc(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layout,
        const std::vector<std::pair<int, int>>& scEdges,
        const std::vector<double>& scWeights,
        ShortcutUsageMode usageMode) {
    const int N = (int)nodes.size();
    auto manhattan = [&](int a, int b) {
        return std::abs(nodes[a].x - nodes[b].x)
             + std::abs(nodes[a].y - nodes[b].y);
    };

    std::vector<std::vector<std::pair<int, double>>> ringAdj(N);
    for (const auto& te : layout.tourEdges) {
        const double w = manhattan(te.from, te.to);
        ringAdj[te.from].push_back({te.to, w});
        ringAdj[te.to].push_back({te.from, w});
    }

    auto spOnAdj = [&](const std::vector<std::vector<std::pair<int, double>>>& adj,
                       int src, int dst) {
        std::vector<double> dist(N, std::numeric_limits<double>::infinity());
        dist[src] = 0.0;
        using Q = std::pair<double, int>;
        std::priority_queue<Q, std::vector<Q>, std::greater<Q>> pq;
        pq.push({0.0, src});
        while (!pq.empty()) {
            const auto [d, u] = pq.top();
            pq.pop();
            if (d > dist[u] + MILP_EPS) continue;
            if (u == dst) return d;
            for (const auto& [v, w] : adj[u]) {
                const double nd = d + w;
                if (nd + MILP_EPS < dist[v]) {
                    dist[v] = nd;
                    pq.push({nd, v});
                }
            }
        }
        return dist[dst];
    };

    if (usageMode == ShortcutUsageMode::Exclusive) {
        double wc = 0.0;
        for (const auto& [s, t] : D.demands) {
            auto adj = ringAdj;
            for (size_t k = 0; k < scEdges.size(); ++k) {
                const int a = scEdges[k].first, b = scEdges[k].second;
                if (!undirectedEndpointsMatch(a, b, s, t)) continue;
                adj[a].push_back({b, scWeights[k]});
                adj[b].push_back({a, scWeights[k]});
            }
            wc = std::max(wc, spOnAdj(adj, s, t));
        }
        return wc;
    }

    auto adj = ringAdj;
    for (size_t k = 0; k < scEdges.size(); ++k) {
        const int a = scEdges[k].first, b = scEdges[k].second;
        adj[a].push_back({b, scWeights[k]});
        adj[b].push_back({a, scWeights[k]});
    }

    double wc = 0.0;
    for (const auto& [s, t] : D.demands)
        wc = std::max(wc, spOnAdj(adj, s, t));
    return wc;
}

}  // namespace

double computeStageE(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layout,
        const std::vector<int>& selectedShortcutIndices,
        const std::vector<ShortcutPairIndex>& pairs,
        double sMin,
        ShortcutUsageMode usageMode) {
    std::vector<std::pair<int, int>> scEdges;
    std::vector<double> scWeights;
    for (int c : selectedShortcutIndices) {
        const int i = pairs[c].i, j = pairs[c].j;
        scEdges.push_back({i, j});
        scWeights.push_back(hatShortcutLength(nodes, layout, i, j, sMin));
    }
    return shortestPathWc(nodes, D, layout, scEdges, scWeights, usageMode);
}
