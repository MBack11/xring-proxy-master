#include "ShortcutGrid.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>

namespace {

constexpr double kEps = 1e-9;

void mergeCoordImpl(std::vector<double>& coords, double v, double minSep) {
    for (double c : coords) {
        if (std::abs(c - v) < minSep - kEps) return;
    }
    coords.push_back(v);
    std::sort(coords.begin(), coords.end());
}

bool segIsHorizontal(const Segment& s) {
    return std::abs(s.y1 - s.y2) < kEps;
}

bool segIsVertical(const Segment& s) {
    return std::abs(s.x1 - s.x2) < kEps;
}

void offsetSegmentInward(const Segment& s, bool ccw, double sm, Segment& out) {
    if (segIsHorizontal(s)) {
        double dy = (s.x2 > s.x1) ? (ccw ? sm : -sm) : (ccw ? -sm : sm);
        out = {s.x1, s.y1 + dy, s.x2, s.y2 + dy};
    } else {
        double dx = (s.y2 > s.y1) ? (ccw ? -sm : sm) : (ccw ? sm : -sm);
        out = {s.x1 + dx, s.y1, s.x2 + dx, s.y2};
    }
}

}  // namespace

bool ShortcutGrid::isHorizontal(const Segment& s) {
    return segIsHorizontal(s);
}

bool ShortcutGrid::isVertical(const Segment& s) {
    return segIsVertical(s);
}

bool ShortcutGrid::pointOnSegment(double px, double py, const Segment& s) {
    if (isHorizontal(s)) {
        if (std::abs(py - s.y1) > EPS) return false;
        double lo = std::min(s.x1, s.x2), hi = std::max(s.x1, s.x2);
        return px >= lo - EPS && px <= hi + EPS;
    }
    if (isVertical(s)) {
        if (std::abs(px - s.x1) > EPS) return false;
        double lo = std::min(s.y1, s.y2), hi = std::max(s.y1, s.y2);
        return py >= lo - EPS && py <= hi + EPS;
    }
    return false;
}

std::vector<Segment> ShortcutGrid::collectRingSegments(
        const std::vector<EdgeOption>& routing) {
    std::vector<Segment> segs;
    for (const EdgeOption& eo : routing)
        for (const Segment& s : eo.option->segments)
            segs.push_back(s);
    return segs;
}

static int tourPosOf(const std::vector<int>& tour, int nodeId) {
    for (int i = 0; i < (int)tour.size(); ++i)
        if (tour[i] == nodeId) return i;
    return -1;
}

std::vector<Segment> ShortcutGrid::incidentExemptSegmentsForNode(
        int nodeId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing) {
    int N = (int)tour.size();
    int i = tourPosOf(tour, nodeId);
    if (i < 0) return {};

    int eiOut = i;
    int eiIn  = (i - 1 + N) % N;

    const auto& outSegs = routing[eiOut].option->segments;
    const auto& inSegs  = routing[eiIn].option->segments;
    if (outSegs.empty() || inSegs.empty()) return {};

    // Outgoing: first segment leaving the node.
    // Incoming: last segment arriving at the node (the ring leg touching this node).
    return {outSegs.front(), inSegs.back()};
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 3.1 — minimum parallel-edge gap d
// ─────────────────────────────────────────────────────────────────────────────

bool ShortcutGrid::computeParallelGap(
        const std::vector<Segment>& ringSegments,
        const std::vector<Node>& nodes,
        ParallelGapInfo& out) {
    std::vector<Segment> horiz, vert;
    for (const Segment& s : ringSegments) {
        if (isHorizontal(s)) horiz.push_back(s);
        else if (isVertical(s)) vert.push_back(s);
    }

    double bestD = std::numeric_limits<double>::infinity();
    Segment bestA{}, bestB{};
    bool bestHoriz = true;

    auto considerHoriz = [&](const Segment& a, const Segment& b) {
        double dy = std::abs(a.y1 - b.y1);
        if (dy < EPS) return;
        if (dy < bestD) {
            bestD = dy;
            bestA = a;
            bestB = b;
            bestHoriz = true;
        }
    };

    auto considerVert = [&](const Segment& a, const Segment& b) {
        double dx = std::abs(a.x1 - b.x1);
        if (dx < EPS) return;
        if (dx < bestD) {
            bestD = dx;
            bestA = a;
            bestB = b;
            bestHoriz = false;
        }
    };

    for (int i = 0; i < (int)horiz.size(); ++i)
        for (int j = i + 1; j < (int)horiz.size(); ++j)
            considerHoriz(horiz[i], horiz[j]);

    for (int i = 0; i < (int)vert.size(); ++i)
        for (int j = i + 1; j < (int)vert.size(); ++j)
            considerVert(vert[i], vert[j]);

    if (!std::isfinite(bestD)) return false;

    out.d = bestD;
    out.horizontalGap = bestHoriz;
    out.edgeA = bestA;
    out.edgeB = bestB;

    double cx = 0.0, cy = 0.0;
    for (const Node& n : nodes) { cx += n.x; cy += n.y; }
    cx /= nodes.size();
    cy /= nodes.size();

    if (bestHoriz) {
        double yA = bestA.y1, yB = bestB.y1;
        if (yA > yB) std::swap(yA, yB);
        out.anchorY = 0.5 * (yA + yB);

        auto xRange = [](const Segment& s) {
            return std::make_pair(std::min(s.x1, s.x2), std::max(s.x1, s.x2));
        };
        auto [aLo, aHi] = xRange(bestA);
        auto [bLo, bHi] = xRange(bestB);
        double oLo = std::max(aLo, bLo), oHi = std::min(aHi, bHi);
        out.anchorX = (oLo <= oHi + EPS) ? 0.5 * (oLo + oHi) : cx;
    } else {
        double xA = bestA.x1, xB = bestB.x1;
        if (xA > xB) std::swap(xA, xB);
        out.anchorX = 0.5 * (xA + xB);

        auto yRange = [](const Segment& s) {
            return std::make_pair(std::min(s.y1, s.y2), std::max(s.y1, s.y2));
        };
        auto [aLo, aHi] = yRange(bestA);
        auto [bLo, bHi] = yRange(bestB);
        double oLo = std::max(aLo, bLo), oHi = std::min(aHi, bHi);
        out.anchorY = (oLo <= oHi + EPS) ? 0.5 * (oLo + oHi) : cy;
    }

    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Coordinate propagation at uniform step
// ─────────────────────────────────────────────────────────────────────────────

std::vector<double> ShortcutGrid::spacedCoords(
        double anchor, double lo, double hi, double step) {
    std::set<double> coords;
    coords.insert(anchor);
    for (double c = anchor - step; c >= lo - EPS; c -= step) coords.insert(c);
    for (double c = anchor + step; c <= hi + EPS; c += step) coords.insert(c);
    return std::vector<double>(coords.begin(), coords.end());
}

void ShortcutGrid::mergeCoord(std::vector<double>& coords, double v, double minSep) {
    mergeCoordImpl(coords, v, minSep);
}

namespace {

void addBottleneckCoords(
        std::vector<double>& xs,
        std::vector<double>& ys,
        const ParallelGapInfo& gap,
        double sMin,
        double laneStep) {
    auto xRange = [](const Segment& s) {
        return std::make_pair(std::min(s.x1, s.x2), std::max(s.x1, s.x2));
    };
    auto yRange = [](const Segment& s) {
        return std::make_pair(std::min(s.y1, s.y2), std::max(s.y1, s.y2));
    };
    const double fineStep = std::max(sMin, laneStep);

    if (gap.horizontalGap) {
        const double yLo = std::min(gap.edgeA.y1, gap.edgeB.y1);
        const double yHi = std::max(gap.edgeA.y1, gap.edgeB.y1);
        for (double y = yLo; y <= yHi + kEps; y += fineStep)
            mergeCoordImpl(ys, y, sMin);
        mergeCoordImpl(ys, gap.anchorY, sMin);

        const auto [aLo, aHi] = xRange(gap.edgeA);
        const auto [bLo, bHi] = xRange(gap.edgeB);
        const double xLo = std::max(aLo, bLo);
        const double xHi = std::min(aHi, bHi);
        if (xLo <= xHi + kEps) {
            for (double x = xLo; x <= xHi + kEps; x += fineStep)
                mergeCoordImpl(xs, x, sMin);
        }
    } else {
        const double xLo = std::min(gap.edgeA.x1, gap.edgeB.x1);
        const double xHi = std::max(gap.edgeA.x1, gap.edgeB.x1);
        for (double x = xLo; x <= xHi + kEps; x += fineStep)
            mergeCoordImpl(xs, x, sMin);
        mergeCoordImpl(xs, gap.anchorX, sMin);

        const auto [aLo, aHi] = yRange(gap.edgeA);
        const auto [bLo, bHi] = yRange(gap.edgeB);
        const double yLo = std::max(aLo, bLo);
        const double yHi = std::min(aHi, bHi);
        if (yLo <= yHi + kEps) {
            for (double y = yLo; y <= yHi + kEps; y += fineStep)
                mergeCoordImpl(ys, y, sMin);
        }
    }
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Global clearance projection lines (s_min inward from ring boundaries)
// ─────────────────────────────────────────────────────────────────────────────

void ShortcutGrid::collectGlobalClearanceLines(
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<Node>& nodes,
        double sMin,
        std::set<double>& clearanceXs,
        std::set<double>& clearanceYs) {
    clearanceXs.clear();
    clearanceYs.clear();
    if (ringSegments.empty() || tour.empty()) return;

    double area = 0.0;
    int N = (int)tour.size();
    std::map<int, std::pair<double, double>> pos;
    for (const Node& n : nodes) pos[n.id] = {n.x, n.y};
    for (int i = 0; i < N; ++i) {
        auto [x1, y1] = pos[tour[i]];
        auto [x2, y2] = pos[tour[(i + 1) % N]];
        area += x1 * y2 - x2 * y1;
    }
    const bool ccw = area > 0.0;

    for (const Segment& s : ringSegments) {
        Segment inset{};
        offsetSegmentInward(s, ccw, sMin, inset);
        if (segIsHorizontal(s))
            clearanceYs.insert(inset.y1);
        else if (segIsVertical(s))
            clearanceXs.insert(inset.x1);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Usability classification (Step 3.7)
// ─────────────────────────────────────────────────────────────────────────────

static bool otherEndpointOf(const Segment& s, double nx, double ny, double& ox, double& oy) {
    if (std::abs(s.x1 - nx) < kEps && std::abs(s.y1 - ny) < kEps) {
        ox = s.x2; oy = s.y2; return true;
    }
    if (std::abs(s.x2 - nx) < kEps && std::abs(s.y2 - ny) < kEps) {
        ox = s.x1; oy = s.y1; return true;
    }
    return false;
}

static int nodeIdAtCoord(const std::vector<Node>& nodes, double x, double y) {
    for (const Node& n : nodes) {
        if (std::abs(n.x - x) < kEps && std::abs(n.y - y) < kEps)
            return n.id;
    }
    return -1;
}

static bool pointOnClippedExemptSegment(
        double px, double py,
        const Segment& s,
        double nx, double ny,
        int srcId, int destId,
        const std::vector<Node>& nodes) {
    if (!segIsHorizontal(s) && !segIsVertical(s)) return false;

    if (segIsHorizontal(s)) {
        if (std::abs(py - s.y1) > kEps) return false;
        double lo = std::min(s.x1, s.x2), hi = std::max(s.x1, s.x2);
        if (px < lo - kEps || px > hi + kEps) return false;
    } else {
        if (std::abs(px - s.x1) > kEps) return false;
        double lo = std::min(s.y1, s.y2), hi = std::max(s.y1, s.y2);
        if (py < lo - kEps || py > hi + kEps) return false;
    }

    if (std::abs(px - nx) < kEps && std::abs(py - ny) < kEps)
        return true;

    double ox = 0.0, oy = 0.0;
    if (!otherEndpointOf(s, nx, ny, ox, oy)) return false;

    if (std::abs(px - ox) < kEps && std::abs(py - oy) < kEps) {
        int oid = nodeIdAtCoord(nodes, ox, oy);
        return oid == srcId || oid == destId;
    }

    return true;
}

bool ShortcutGrid::isOnTourNode(
        double x,
        double y,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes) {
    for (const Node& n : nodes) {
        if (n.id == srcId || n.id == destId) continue;
        bool onTour = false;
        for (int tid : tour) {
            if (tid == n.id) { onTour = true; break; }
        }
        if (!onTour) continue;
        if (std::abs(x - n.x) < EPS && std::abs(y - n.y) < EPS)
            return true;
    }
    return false;
}

bool ShortcutGrid::classifyUsability(
        double x,
        double y,
        const std::vector<Segment>& ringSegments,
        int srcId,
        int destId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        const std::vector<Node>& nodes) {
    // Ring tour nodes (except src/dest) cannot be used as shortcut waypoints.
    if (isOnTourNode(x, y, tour, srcId, destId, nodes))
        return false;

    // Off-ring grid points: normal (finite) cost.
    bool onRing = false;
    for (const Segment& s : ringSegments) {
        if (pointOnSegment(x, y, s)) { onRing = true; break; }
    }
    if (!onRing) return true;

    // On-ring points: infinite cost unless on an exempt first segment at src/dest.
    std::map<int, std::pair<double, double>> pos;
    for (const Node& n : nodes) pos[n.id] = {n.x, n.y};

    for (int nodeId : {srcId, destId}) {
        auto [nx, ny] = pos[nodeId];
        for (const Segment& s : incidentExemptSegmentsForNode(nodeId, tour, routing)) {
            if (pointOnClippedExemptSegment(x, y, s, nx, ny, srcId, destId, nodes))
                return true;
        }
    }
    return false;
}

void ShortcutGrid::addPoint(
        std::vector<GridPoint>& pts,
        double x,
        double y,
        const std::string& source,
        bool usable) {
    for (const GridPoint& p : pts) {
        if (std::abs(p.x - x) < EPS && std::abs(p.y - y) < EPS) {
            return;
        }
    }
    pts.push_back({x, y, source, usable});
}

// ─────────────────────────────────────────────────────────────────────────────
// Main build
// ─────────────────────────────────────────────────────────────────────────────

ShortcutGridResult ShortcutGrid::build(
        const std::vector<Node>& nodes,
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        int srcId,
        int destId,
        double sMin) {
    ShortcutGridResult result;
    result.sMin = sMin;

    if (!computeParallelGap(ringSegments, nodes, result.gapInfo)) {
        result.message = "Could not compute parallel-edge gap d (no parallel ring-edge pairs).";
        return result;
    }

    result.d = result.gapInfo.d;
    const double half = result.d * 0.5;

    if (half < sMin - EPS) {
        result.success = false;
        result.preconditionOk = false;
        result.message = "No shortcut possible: d/2 (" + std::to_string(half)
                       + " mm) is smaller than s_min (" + std::to_string(sMin) + " mm).";
        return result;
    }
    result.preconditionOk = true;

    double minX = nodes[0].x, maxX = nodes[0].x;
    double minY = nodes[0].y, maxY = nodes[0].y;
    for (const Node& n : nodes) {
        minX = std::min(minX, n.x); maxX = std::max(maxX, n.x);
        minY = std::min(minY, n.y); maxY = std::max(maxY, n.y);
    }
    const double margin = result.d;

    const double anchorX = result.gapInfo.anchorX;
    const double anchorY = result.gapInfo.anchorY;

    const double spanX = (maxX - minX) + 2.0 * margin;
    const double spanY = (maxY - minY) + 2.0 * margin;
    const double coarseStep = std::max(
        sMin,
        std::max(spanX, spanY) / std::max(1, MAX_GRID_LINES - 1));
    const double laneStep = std::max(sMin, half);
    const double step = (coarseStep > laneStep + EPS) ? coarseStep : laneStep;

    std::vector<double> xs = spacedCoords(anchorX, minX - margin, maxX + margin, step);
    std::vector<double> ys = spacedCoords(anchorY, minY - margin, maxY + margin, step);
    result.gridStep = step;
    result.xsCount = (int)xs.size();
    result.ysCount = (int)ys.size();

    for (const Node& n : nodes) {
        mergeCoord(xs, n.x, sMin);
        mergeCoord(ys, n.y, sMin);
    }
    if (step > laneStep + EPS)
        addBottleneckCoords(xs, ys, result.gapInfo, sMin, half);
    result.xsCount = (int)xs.size();
    result.ysCount = (int)ys.size();

    std::vector<GridPoint> pts;

    auto usableAt = [&](double x, double y) {
        return classifyUsability(x, y, ringSegments, srcId, destId, tour, routing, nodes);
    };

    // 3.4 — main grid
    for (double x : xs)
        for (double y : ys)
            addPoint(pts, x, y, "main_grid", usableAt(x, y));

    // 3.5 — node grid-lines (vertical x=node.x and horizontal y=node.y)
    for (const Node& n : nodes) {
        for (double y : ys) {
            if (std::abs(n.y - y) < EPS && n.id != srcId && n.id != destId)
                continue;
            addPoint(pts, n.x, y, "node_line", usableAt(n.x, y));
        }
        for (double x : xs) {
            if (std::abs(n.x - x) < EPS && n.id != srcId && n.id != destId)
                continue;
            addPoint(pts, x, n.y, "node_line", usableAt(x, n.y));
        }
    }

    // 3.6 — global clearance projection lines at s_min (full-layout lanes)
    std::set<double> clearanceXs, clearanceYs;
    collectGlobalClearanceLines(ringSegments, tour, nodes, sMin, clearanceXs, clearanceYs);

    std::vector<double> clearanceXv(clearanceXs.begin(), clearanceXs.end());
    std::vector<double> clearanceYv(clearanceYs.begin(), clearanceYs.end());

    // All intersections of global clearance horizontal × vertical lines
    for (double cx : clearanceXv)
        for (double cy : clearanceYv)
            addPoint(pts, cx, cy, "inner_ring", usableAt(cx, cy));

    // Continuous lanes: each clearance line spans the full layout at d/2 spacing
    for (double cx : clearanceXv)
        for (double y : ys)
            addPoint(pts, cx, y, "inner_ring", usableAt(cx, y));
    for (double cy : clearanceYv)
        for (double x : xs)
            addPoint(pts, x, cy, "inner_ring", usableAt(x, cy));

    for (const Node& n : nodes) {
        bool forceUsable = (n.id == srcId || n.id == destId);
        addPoint(pts, n.x, n.y, "node_anchor",
                 forceUsable || usableAt(n.x, n.y));
    }

    result.points = std::move(pts);
    result.success = true;
    if (result.message.empty())
        result.message = "Grid built successfully.";
    return result;
}

namespace {

struct GridCacheKey {
    std::vector<int> tour;
    std::vector<int> routingOpts;
    int srcId = -1;
    int destId = -1;
    long long sMinKey = 0;

    bool operator==(const GridCacheKey& other) const {
        return srcId == other.srcId
            && destId == other.destId
            && sMinKey == other.sMinKey
            && tour == other.tour
            && routingOpts == other.routingOpts;
    }
};

struct GridCacheKeyHash {
    size_t operator()(const GridCacheKey& key) const {
        size_t h = std::hash<int>()(key.srcId);
        h ^= std::hash<int>()(key.destId) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<long long>()(key.sMinKey) + 0x9e3779b9 + (h << 6) + (h >> 2);
        for (int v : key.tour)
            h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        for (int v : key.routingOpts)
            h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

long long quantizeSMin(double sMin) {
    return static_cast<long long>(std::llround(sMin / 1e-9));
}

GridCacheKey makeGridCacheKey(
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        int srcId,
        int destId,
        double sMin) {
    GridCacheKey key;
    key.tour = tour;
    key.routingOpts.reserve(routing.size());
    for (const EdgeOption& eo : routing)
        key.routingOpts.push_back(eo.optionIndex);
    key.srcId = srcId;
    key.destId = destId;
    key.sMinKey = quantizeSMin(sMin);
    return key;
}

std::unordered_map<GridCacheKey, ShortcutGridResult, GridCacheKeyHash>& gridCache() {
    static std::unordered_map<GridCacheKey, ShortcutGridResult, GridCacheKeyHash> cache;
    return cache;
}

ShortcutGrid::BuildCacheStats& gridCacheStatsMutable() {
    static ShortcutGrid::BuildCacheStats stats;
    return stats;
}

}  // namespace

ShortcutGrid::BuildCacheStats ShortcutGrid::buildCacheStats() {
    BuildCacheStats stats = gridCacheStatsMutable();
    stats.entries = gridCache().size();
    return stats;
}

void ShortcutGrid::resetBuildCacheStats() {
    gridCacheStatsMutable() = {};
    gridCache().clear();
}

ShortcutGridResult ShortcutGrid::buildCached(
        const std::vector<Node>& nodes,
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        int srcId,
        int destId,
        double sMin) {
    const GridCacheKey key = makeGridCacheKey(tour, routing, srcId, destId, sMin);
    auto& cache = gridCache();
    auto it = cache.find(key);
    if (it != cache.end()) {
        ++gridCacheStatsMutable().hits;
        return it->second;
    }

    ++gridCacheStatsMutable().misses;
    ShortcutGridResult result = build(
        nodes, ringSegments, tour, routing, srcId, destId, sMin);
    cache.emplace(key, result);
    return result;
}
