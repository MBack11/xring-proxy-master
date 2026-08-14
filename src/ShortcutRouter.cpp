#include "ShortcutRouter.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace {

long long coordKey(double v) {
    return static_cast<long long>(std::llround(v / 1e-6));
}

bool pointOnSegment(double px, double py, const Segment& s) {
    static constexpr double EPS = 1e-9;
    if (std::abs(s.x1 - s.x2) < EPS) {
        if (std::abs(px - s.x1) > EPS) return false;
        double lo = std::min(s.y1, s.y2), hi = std::max(s.y1, s.y2);
        return py >= lo - EPS && py <= hi + EPS;
    }
    if (std::abs(s.y1 - s.y2) < EPS) {
        if (std::abs(py - s.y1) > EPS) return false;
        double lo = std::min(s.x1, s.x2), hi = std::max(s.x1, s.x2);
        return px >= lo - EPS && px <= hi + EPS;
    }
    return false;
}

}  // namespace

namespace {
// Count axis-aligned direction changes on a polyline, ignoring sub-s_min stubs
// (grid artifacts such as a tiny lateral step from a node onto a clearance lane).
int countBendsIgnoringShortStubs(
        const std::vector<std::pair<double, double>>& verts,
        double minSegmentLength) {
    constexpr double kEps = 1e-9;
    struct Move {
        int axis;  // 0=H, 1=V
        int sign;
        double len;
    };
    std::vector<Move> moves;
    for (size_t i = 0; i + 1 < verts.size(); ++i) {
        const double dx = verts[i + 1].first - verts[i].first;
        const double dy = verts[i + 1].second - verts[i].second;
        const double len = std::abs(dx) + std::abs(dy);
        if (len < kEps) continue;
        Move m;
        if (std::abs(dx) >= std::abs(dy)) {
            m.axis = 0;
            m.sign = (dx > 0) ? 1 : -1;
        } else {
            m.axis = 1;
            m.sign = (dy > 0) ? 1 : -1;
        }
        m.len = len;
        if (!moves.empty() && moves.back().axis == m.axis && moves.back().sign == m.sign)
            moves.back().len += m.len;
        else
            moves.push_back(m);
    }

    std::vector<Move> kept;
    for (const Move& m : moves) {
        if (m.len + kEps < minSegmentLength) continue;
        if (!kept.empty() && kept.back().axis == m.axis && kept.back().sign == m.sign)
            kept.back().len += m.len;
        else
            kept.push_back(m);
    }

    int bends = 0;
    for (size_t i = 1; i < kept.size(); ++i) {
        if (kept[i].axis != kept[i - 1].axis || kept[i].sign != kept[i - 1].sign)
            ++bends;
    }
    return bends;
}

}  // namespace

namespace {

long long ringCellKey(int ix, int iy) {
    return (static_cast<long long>(ix) << 32) ^ static_cast<unsigned int>(iy);
}

}  // namespace

void ShortcutRouter::RingSegmentIndex::build(
        const std::vector<Segment>& ring, double cell) {
    segments = &ring;
    cellSize = std::max(cell, 1e-6);
    buckets.clear();
    for (int i = 0; i < (int)ring.size(); ++i) {
        const Segment& s = ring[i];
        const double xLo = std::min(s.x1, s.x2);
        const double xHi = std::max(s.x1, s.x2);
        const double yLo = std::min(s.y1, s.y2);
        const double yHi = std::max(s.y1, s.y2);
        const int ix0 = static_cast<int>(std::floor(xLo / cellSize));
        const int ix1 = static_cast<int>(std::floor(xHi / cellSize));
        const int iy0 = static_cast<int>(std::floor(yLo / cellSize));
        const int iy1 = static_cast<int>(std::floor(yHi / cellSize));
        for (int ix = ix0; ix <= ix1; ++ix)
            for (int iy = iy0; iy <= iy1; ++iy)
                buckets[ringCellKey(ix, iy)].push_back(i);
    }
}

void ShortcutRouter::RingSegmentIndex::forEachCandidate(
        const Segment& edge, const std::function<void(int)>& fn) const {
    if (!segments) return;
    const double xLo = std::min(edge.x1, edge.x2);
    const double xHi = std::max(edge.x1, edge.x2);
    const double yLo = std::min(edge.y1, edge.y2);
    const double yHi = std::max(edge.y1, edge.y2);
    const int ix0 = static_cast<int>(std::floor(xLo / cellSize)) - 1;
    const int ix1 = static_cast<int>(std::floor(xHi / cellSize)) + 1;
    const int iy0 = static_cast<int>(std::floor(yLo / cellSize)) - 1;
    const int iy1 = static_cast<int>(std::floor(yHi / cellSize)) + 1;
    std::unordered_set<int> seen;
    for (int ix = ix0; ix <= ix1; ++ix) {
        for (int iy = iy0; iy <= iy1; ++iy) {
            auto it = buckets.find(ringCellKey(ix, iy));
            if (it == buckets.end()) continue;
            for (int idx : it->second) {
                if (seen.insert(idx).second) fn(idx);
            }
        }
    }
}

static int gShortcutCrossDim = 2;

int ShortcutRouter::stateId(int pointIdx, GridDir dir) {
    return stateId(pointIdx, dir, 0);
}

int ShortcutRouter::stateId(int pointIdx, GridDir dir, int crossSlot) {
    return (pointIdx * NUM_DIRS + static_cast<int>(dir)) * gShortcutCrossDim + crossSlot;
}

std::pair<int, GridDir> ShortcutRouter::decodeState(int sid) {
    auto [pIdx, dir, cross] = decodeStateEx(sid);
    (void)cross;
    return {pIdx, dir};
}

std::tuple<int, GridDir, int> ShortcutRouter::decodeStateEx(int sid) {
    int crossSlot = sid % gShortcutCrossDim;
    int rest = sid / gShortcutCrossDim;
    int dir = rest % NUM_DIRS;
    int pIdx = rest / NUM_DIRS;
    return {pIdx, static_cast<GridDir>(dir), crossSlot};
}

GridDir ShortcutRouter::directionBetween(double x1, double y1, double x2, double y2) {
    if (std::abs(x2 - x1) > EPS) {
        return (x2 > x1) ? GridDir::RIGHT : GridDir::LEFT;
    }
    return (y2 > y1) ? GridDir::UP : GridDir::DOWN;
}

bool ShortcutRouter::lexLess(double d1, int b1, double d2, int b2) {
    if (std::abs(d1 - d2) > EPS) return d1 < d2;
    return b1 < b2;
}

bool ShortcutRouter::lexEq(double d1, int b1, double d2, int b2) {
    return std::abs(d1 - d2) <= EPS && b1 == b2;
}

bool ShortcutRouter::isParallelOverlap(const Segment& s1, const Segment& s2) {
    bool isH1 = std::abs(s1.y1 - s1.y2) < EPS;
    bool isH2 = std::abs(s2.y1 - s2.y2) < EPS;
    bool isV1 = std::abs(s1.x1 - s1.x2) < EPS;
    bool isV2 = std::abs(s2.x1 - s2.x2) < EPS;

    if (isH1 && isH2 && std::abs(s1.y1 - s2.y1) < EPS) {
        double lo1 = std::min(s1.x1, s1.x2), hi1 = std::max(s1.x1, s1.x2);
        double lo2 = std::min(s2.x1, s2.x2), hi2 = std::max(s2.x1, s2.x2);
        return lo1 < hi2 - EPS && lo2 < hi1 - EPS;
    }
    if (isV1 && isV2 && std::abs(s1.x1 - s2.x1) < EPS) {
        double lo1 = std::min(s1.y1, s1.y2), hi1 = std::max(s1.y1, s1.y2);
        double lo2 = std::min(s2.y1, s2.y2), hi2 = std::max(s2.y1, s2.y2);
        return lo1 < hi2 - EPS && lo2 < hi1 - EPS;
    }
    return false;
}

bool ShortcutRouter::isTJunction(double px, double py, const Segment& obs, bool gridIsH) {
    bool obsIsH = std::abs(obs.y1 - obs.y2) < EPS;
    bool obsIsV = std::abs(obs.x1 - obs.x2) < EPS;
    if (gridIsH == obsIsH) return false;
    if (obsIsH) {
        if (std::abs(py - obs.y1) > EPS) return false;
        double lo = std::min(obs.x1, obs.x2), hi = std::max(obs.x1, obs.x2);
        return px > lo + EPS && px < hi - EPS;
    }
    if (obsIsV) {
        if (std::abs(px - obs.x1) > EPS) return false;
        double lo = std::min(obs.y1, obs.y2), hi = std::max(obs.y1, obs.y2);
        return py > lo + EPS && py < hi - EPS;
    }
    return false;
}

bool ShortcutRouter::segmentsGeometricallyIntersect(const Segment& a, const Segment& b) {
    if (isStrictCrossing(a, b) || isParallelOverlap(a, b)) return true;

    bool aH = std::abs(a.y1 - a.y2) < EPS;
    bool aV = std::abs(a.x1 - a.x2) < EPS;
    bool bH = std::abs(b.y1 - b.y2) < EPS;
    bool bV = std::abs(b.x1 - b.x2) < EPS;

    if (aH || aV) {
        if (isTJunction(a.x1, a.y1, b, aH) || isTJunction(a.x2, a.y2, b, aH)) return true;
    }
    if (bH || bV) {
        if (isTJunction(b.x1, b.y1, a, bH) || isTJunction(b.x2, b.y2, a, bH)) return true;
    }

    auto onSeg = [](double px, double py, const Segment& s) {
        return pointOnSegment(px, py, s);
    };
    if (onSeg(a.x1, a.y1, b) || onSeg(a.x2, a.y2, b)) return true;
    if (onSeg(b.x1, b.y1, a) || onSeg(b.x2, b.y2, a)) return true;
    return false;
}

bool ShortcutRouter::edgeIntersectsShortcut(const Segment& edge, const Shortcut& sc) {
    for (const Segment& s : sc.path) {
        if (segmentsGeometricallyIntersect(edge, s)) return true;
    }
    return false;
}

bool ShortcutRouter::shortcutsIntersect(const Shortcut& a, const Shortcut& b) {
    for (const Segment& sa : a.path) {
        for (const Segment& sb : b.path) {
            if (segmentsGeometricallyIntersect(sa, sb)) return true;
        }
    }
    return false;
}

int ShortcutRouter::countShortcutIntersections(
        const std::vector<Shortcut>& shortcuts, int index) {
    int count = 0;
    for (int j = 0; j < (int)shortcuts.size(); ++j) {
        if (j == index) continue;
        if (shortcutsIntersect(shortcuts[index], shortcuts[j])) ++count;
    }
    return count;
}

bool ShortcutRouter::isStrictCrossing(const Segment& s1, const Segment& s2) {
    bool isH1 = std::abs(s1.y1 - s1.y2) < EPS;
    bool isV1 = std::abs(s1.x1 - s1.x2) < EPS;
    bool isH2 = std::abs(s2.y1 - s2.y2) < EPS;
    bool isV2 = std::abs(s2.x1 - s2.x2) < EPS;

    if ((isH1 && isH2) || (isV1 && isV2)) return false;

    if (isH1 && isV2) {
        double minX = std::min(s1.x1, s1.x2), maxX = std::max(s1.x1, s1.x2);
        double minY = std::min(s2.y1, s2.y2), maxY = std::max(s2.y1, s2.y2);
        return (s2.x1 > minX + EPS && s2.x1 < maxX - EPS)
            && (s1.y1 > minY + EPS && s1.y1 < maxY - EPS);
    }
    if (isV1 && isH2) {
        double minY = std::min(s1.y1, s1.y2), maxY = std::max(s1.y1, s1.y2);
        double minX = std::min(s2.x1, s2.x2), maxX = std::max(s2.x1, s2.x2);
        return (s1.x1 > minX + EPS && s1.x1 < maxX - EPS)
            && (s2.y1 > minY + EPS && s2.y1 < maxY - EPS);
    }
    return false;
}

bool ShortcutRouter::edgeIsOnExemptSegment(
        const Segment& edge,
        const std::vector<Segment>& exemptSegments) {
    auto contained = [&](const Segment& container) {
        bool eH = std::abs(edge.y1 - edge.y2) < EPS;
        bool cH = std::abs(container.y1 - container.y2) < EPS;
        bool eV = std::abs(edge.x1 - edge.x2) < EPS;
        bool cV = std::abs(container.x1 - container.x2) < EPS;
        if (eH && cH && std::abs(edge.y1 - container.y1) < EPS) {
            double eLo = std::min(edge.x1, edge.x2), eHi = std::max(edge.x1, edge.x2);
            double cLo = std::min(container.x1, container.x2), cHi = std::max(container.x1, container.x2);
            return eLo >= cLo - EPS && eHi <= cHi + EPS;
        }
        if (eV && cV && std::abs(edge.x1 - container.x1) < EPS) {
            double eLo = std::min(edge.y1, edge.y2), eHi = std::max(edge.y1, edge.y2);
            double cLo = std::min(container.y1, container.y2), cHi = std::max(container.y1, container.y2);
            return eLo >= cLo - EPS && eHi <= cHi + EPS;
        }
        return false;
    };

    for (const Segment& ex : exemptSegments) {
        if (contained(ex)) return true;
    }
    return false;
}

bool ShortcutRouter::edgePassesThroughTourNode(
        double x1, double y1, double x2, double y2,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes) {
    Segment edge{x1, y1, x2, y2};
    std::map<int, std::pair<double, double>> pos;
    for (const Node& n : nodes) pos[n.id] = {n.x, n.y};

    for (int nid : tour) {
        if (nid == srcId || nid == destId) continue;
        auto it = pos.find(nid);
        if (it == pos.end()) continue;
        auto [nx, ny] = it->second;
        if (pointOnSegment(nx, ny, edge))
            return true;
    }
    return false;
}

bool ShortcutRouter::gridEdgeAllowed(
        double x1, double y1, double x2, double y2,
        const std::vector<Segment>& ringSegments,
        const std::vector<Segment>& exemptSegments,
        const std::vector<Segment>& hardShortcutSegments,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes) {
    return !classifyGridEdge(
        x1, y1, x2, y2, ringSegments, exemptSegments, hardShortcutSegments,
        {}, tour, srcId, destId, nodes).hardBlocked;
}

ShortcutRouter::EdgeCrossInfo ShortcutRouter::classifyGridEdge(
        double x1, double y1, double x2, double y2,
        const std::vector<Segment>& ringSegments,
        const std::vector<Segment>& exemptSegments,
        const std::vector<Segment>& hardShortcutSegments,
        const std::vector<std::pair<int, std::vector<Segment>>>& softShortcuts,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes,
        const RingSegmentIndex* ringIndex) {
    EdgeCrossInfo info;
    Segment edge{x1, y1, x2, y2};

    if (edgeIsOnExemptSegment(edge, exemptSegments)) {
        info.hardBlocked = false;
        info.softShortcutIdx = -1;
        return info;
    }

    if (edgePassesThroughTourNode(x1, y1, x2, y2, tour, srcId, destId, nodes)) {
        return info;
    }

    const double exLo = std::min(x1, x2);
    const double exHi = std::max(x1, x2);
    const double eyLo = std::min(y1, y2);
    const double eyHi = std::max(y1, y2);

    auto bboxOverlap = [&](const Segment& ring) {
        const double rxLo = std::min(ring.x1, ring.x2);
        const double rxHi = std::max(ring.x1, ring.x2);
        const double ryLo = std::min(ring.y1, ring.y2);
        const double ryHi = std::max(ring.y1, ring.y2);
        return exHi >= rxLo - EPS && rxHi >= exLo - EPS
            && eyHi >= ryLo - EPS && ryHi >= eyLo - EPS;
    };

    bool ringCrossing = false;
    if (ringIndex) {
        ringIndex->forEachCandidate(edge, [&](int idx) {
            if (ringCrossing) return;
            const Segment& ring = (*ringIndex->segments)[idx];
            if (isStrictCrossing(edge, ring)) ringCrossing = true;
        });
    } else {
        for (const Segment& ring : ringSegments) {
            if (!bboxOverlap(ring)) continue;
            if (isStrictCrossing(edge, ring)) {
                ringCrossing = true;
                break;
            }
        }
    }
    if (ringCrossing) return info;

    for (const Segment& ring : ringSegments) {
        if (!bboxOverlap(ring)) continue;
        if (isParallelOverlap(edge, ring)) return info;
    }

    for (const Segment& obs : hardShortcutSegments) {
        if (segmentsGeometricallyIntersect(edge, obs)) return info;
    }

    info.hardBlocked = false;
    info.softShortcutIdx = -1;

    int crossCount = 0;
    for (const auto& [scIdx, segs] : softShortcuts) {
        bool intersects = false;
        for (const Segment& obs : segs) {
            if (segmentsGeometricallyIntersect(edge, obs)) {
                intersects = true;
                break;
            }
        }
        if (intersects) {
            ++crossCount;
            info.softShortcutIdx = (crossCount == 1) ? scIdx : -2;
            if (crossCount > 1) {
                info.hardBlocked = true;
                return info;
            }
        }
    }

    return info;
}

int ShortcutRouter::findNodeIndex(const Graph& g, double x, double y) {
    for (int i = 0; i < (int)g.xs.size(); ++i) {
        if (std::abs(g.xs[i] - x) <= EPS && std::abs(g.ys[i] - y) <= EPS)
            return i;
    }
    return -1;
}

void ShortcutRouter::connectLine(
        Graph& g,
        const std::vector<int>& indices,
        bool sortByY) {
    std::vector<int> order;
    order.reserve(indices.size());
    for (int i : indices)
        if (g.usable[i]) order.push_back(i);
    if (order.size() < 2) return;

    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (sortByY) {
            if (std::abs(g.ys[a] - g.ys[b]) > EPS) return g.ys[a] < g.ys[b];
            return g.xs[a] < g.xs[b];
        }
        if (std::abs(g.xs[a] - g.xs[b]) > EPS) return g.xs[a] < g.xs[b];
        return g.ys[a] < g.ys[b];
    });

    for (int k = 0; k + 1 < (int)order.size(); ++k) {
        int u = order[k], v = order[k + 1];
        g.adj[u].push_back(v);
        g.adj[v].push_back(u);
    }
}

std::optional<ShortcutRouter::Graph> ShortcutRouter::buildGraph(
        const std::vector<GridPoint>& gridPoints) {
    Graph g;
    g.xs.reserve(gridPoints.size());
    g.ys.reserve(gridPoints.size());
    g.usable.reserve(gridPoints.size());
    for (const GridPoint& p : gridPoints) {
        g.xs.push_back(p.x);
        g.ys.push_back(p.y);
        g.usable.push_back(p.usable);
    }
    g.adj.assign(gridPoints.size(), {});

    std::map<long long, std::vector<int>> byX;
    std::map<long long, std::vector<int>> byY;
    for (int i = 0; i < (int)gridPoints.size(); ++i) {
        byX[coordKey(g.xs[i])].push_back(i);
        byY[coordKey(g.ys[i])].push_back(i);
    }

    for (auto& [_, idxs] : byX) connectLine(g, idxs, true);
    for (auto& [_, idxs] : byY) connectLine(g, idxs, false);

    for (auto& nbrs : g.adj) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    return g;
}

void ShortcutRouter::enumeratePaths(
        int srcState,
        const std::vector<std::vector<int>>& preds,
        const std::vector<double>& dist,
        const std::vector<int>& bends,
        double optDist,
        int optBends,
        int destPoint,
        int maxPaths,
        std::vector<std::vector<int>>& outStatePaths) {
    std::vector<int> endStates;
    for (int cross = 0; cross < gShortcutCrossDim; ++cross) {
        for (int d = 0; d < NUM_DIRS; ++d) {
            int sid = stateId(destPoint, static_cast<GridDir>(d), cross);
            if (lexEq(dist[sid], bends[sid], optDist, optBends))
                endStates.push_back(sid);
        }
    }
    std::sort(endStates.begin(), endStates.end());

    std::vector<int> pathStates;
    std::set<std::vector<int>> seen;

    std::function<void(int)> dfs = [&](int sid) {
        if ((int)outStatePaths.size() >= maxPaths) return;

        pathStates.push_back(sid);
        if (sid == srcState) {
            std::vector<int> rev = pathStates;
            std::reverse(rev.begin(), rev.end());
            if (!seen.count(rev)) {
                seen.insert(rev);
                outStatePaths.push_back(rev);
            }
            pathStates.pop_back();
            return;
        }

        auto predList = preds[sid];
        std::sort(predList.begin(), predList.end());
        for (int ps : predList) {
            dfs(ps);
            if ((int)outStatePaths.size() >= maxPaths) break;
        }
        pathStates.pop_back();
    };

    for (int es : endStates) dfs(es);
}

ShortcutRouteResult ShortcutRouter::findPaths(
        const std::vector<GridPoint>& gridPoints,
        int srcId,
        int destId,
        const std::vector<Node>& nodes,
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        const std::vector<Shortcut>& existingShortcuts,
        double maxDistance) {
    ShortcutRouteResult result;

    auto graphOpt = buildGraph(gridPoints);
    if (!graphOpt) {
        result.message = "Failed to build routing graph.";
        return result;
    }
    Graph g = std::move(*graphOpt);

    double sx = 0.0, sy = 0.0, dx = 0.0, dy = 0.0;
    bool foundSrc = false, foundDest = false;
    for (const Node& n : nodes) {
        if (n.id == srcId) { sx = n.x; sy = n.y; foundSrc = true; }
        if (n.id == destId) { dx = n.x; dy = n.y; foundDest = true; }
    }
    if (!foundSrc || !foundDest) {
        result.message = "Source or destination node not found.";
        return result;
    }

    int srcIdx = findNodeIndex(g, sx, sy);
    int destIdx = findNodeIndex(g, dx, dy);
    if (srcIdx < 0 || destIdx < 0) {
        result.message = "Source/dest coordinates missing from grid.";
        return result;
    }
    if (!g.usable[srcIdx] || !g.usable[destIdx]) {
        result.message = "Source or destination grid point is not usable.";
        return result;
    }

    std::vector<Segment> exemptSegments;
    for (int nodeId : {srcId, destId}) {
        auto segs = ShortcutGrid::incidentExemptSegmentsForNode(nodeId, tour, routing);
        exemptSegments.insert(exemptSegments.end(), segs.begin(), segs.end());
    }

    RingSegmentIndex ringIndex;
    {
        double minX = nodes[0].x, maxX = nodes[0].x;
        double minY = nodes[0].y, maxY = nodes[0].y;
        for (const Node& n : nodes) {
            minX = std::min(minX, n.x); maxX = std::max(maxX, n.x);
            minY = std::min(minY, n.y); maxY = std::max(maxY, n.y);
        }
        const double cell = std::max(
            std::max(maxX - minX, maxY - minY) / 8.0, 0.3);
        ringIndex.build(ringSegments, cell);
    }

    if (srcIdx == destIdx) {
        result.success = true;
        result.primary.vertices = {{sx, sy}};
        result.message = "Trivial path (src == dest).";
        return result;
    }

    std::vector<Segment> hardShortcutSegments;
    std::vector<std::pair<int, std::vector<Segment>>> softShortcuts;
    for (int si = 0; si < (int)existingShortcuts.size(); ++si) {
        const Shortcut& sc = existingShortcuts[si];
        if (countShortcutIntersections(existingShortcuts, si) >= 1) {
            for (const Segment& s : sc.path)
                hardShortcutSegments.push_back(s);
        } else {
            softShortcuts.push_back({si, sc.path});
        }
    }

    const int nPts = (int)g.xs.size();
    gShortcutCrossDim = std::max(2, 1 + (int)existingShortcuts.size());
    const int nStates = nPts * NUM_DIRS * gShortcutCrossDim;
    const double INF = std::numeric_limits<double>::infinity();

    struct SearchOutcome {
        double optDist = std::numeric_limits<double>::infinity();
        int optBends = INT_MAX;
        int crossedIdx = -1;
        std::vector<std::vector<int>> statePaths;
    };

    auto lexLessDistBendCross = [](double d1, int b1, int c1, double d2, int b2, int c2) {
        if (std::abs(d1 - d2) > EPS) return d1 < d2;
        if (b1 != b2) return b1 < b2;
        return c1 < c2;
    };

    auto doSearch = [&]() -> std::optional<SearchOutcome> {
        std::vector<double> dist(nStates, INF);
        std::vector<int> bendCount(nStates, INT_MAX);
        std::vector<std::vector<int>> preds(nStates);

        using PQEntry = std::tuple<double, int, int>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        int srcState = stateId(srcIdx, GridDir::NONE, 0);
        dist[srcState] = 0.0;
        bendCount[srcState] = 0;
        pq.push({0.0, 0, srcState});

        while (!pq.empty()) {
            auto [dCur, bCur, sid] = pq.top();
            pq.pop();
            if (std::abs(dCur - dist[sid]) > EPS || bCur != bendCount[sid]) continue;

            auto [pIdx, dirIn, crossSlot] = decodeStateEx(sid);
            int crossedIdx = (crossSlot == 0) ? -1 : (crossSlot - 1);
            for (int nIdx : g.adj[pIdx]) {
                double x1 = g.xs[pIdx], y1 = g.ys[pIdx];
                double x2 = g.xs[nIdx], y2 = g.ys[nIdx];

                EdgeCrossInfo edgeInfo = classifyGridEdge(
                    x1, y1, x2, y2, ringSegments, exemptSegments, hardShortcutSegments,
                    softShortcuts, tour, srcId, destId, nodes, &ringIndex);
                if (edgeInfo.hardBlocked) continue;

                int newCrossSlot = crossSlot;
                if (edgeInfo.softShortcutIdx == -2) continue;
                if (edgeInfo.softShortcutIdx >= 0) {
                    int si = edgeInfo.softShortcutIdx;
                    if (crossedIdx == -1) newCrossSlot = si + 1;
                    else if (crossedIdx == si) newCrossSlot = crossSlot;
                    else continue;
                }

                GridDir dirOut = directionBetween(x1, y1, x2, y2);
                double edgeLen = std::abs(g.xs[nIdx] - g.xs[pIdx])
                               + std::abs(g.ys[nIdx] - g.ys[pIdx]);
                // Sub-s_min edges are clearance/grid stubs: do not commit direction
                // and do not count a bend (avoids phantom bend on tiny lateral jogs).
                GridDir dirForState = dirOut;
                int addedBend = 0;
                if (edgeLen + EPS < ShortcutGrid::DEFAULT_S_MIN) {
                    dirForState = dirIn;
                    addedBend = 0;
                } else {
                    addedBend = (dirIn != GridDir::NONE && dirIn != dirOut) ? 1 : 0;
                }
                double nd = dCur + edgeLen;
                if (nd > maxDistance + EPS) continue;
                int nb = bCur + addedBend;
                int nsid = stateId(nIdx, dirForState, newCrossSlot);

                if (lexLess(nd, nb, dist[nsid], bendCount[nsid])) {
                    dist[nsid] = nd;
                    bendCount[nsid] = nb;
                    preds[nsid] = {sid};
                    pq.push({nd, nb, nsid});
                } else if (lexEq(nd, nb, dist[nsid], bendCount[nsid])) {
                    preds[nsid].push_back(sid);
                }
            }
        }

        SearchOutcome out;
        int endState = -1;
        int bestCrossSlot = INT_MAX;
        for (int cross = 0; cross < gShortcutCrossDim; ++cross) {
            for (int d = 0; d < NUM_DIRS; ++d) {
                int sid = stateId(destIdx, static_cast<GridDir>(d), cross);
                if (lexLessDistBendCross(dist[sid], bendCount[sid], cross,
                                         out.optDist, out.optBends, bestCrossSlot)) {
                    out.optDist = dist[sid];
                    out.optBends = bendCount[sid];
                    bestCrossSlot = cross;
                    endState = sid;
                }
            }
        }

        if (!std::isfinite(out.optDist) || endState < 0) return std::nullopt;

        enumeratePaths(srcState, preds, dist, bendCount, out.optDist, out.optBends,
                       destIdx, 1 + MAX_ALTERNATIVES, out.statePaths);
        if (out.statePaths.empty()) return std::nullopt;

        int endCrossSlot = std::get<2>(decodeStateEx(endState));
        if (endCrossSlot > 0) {
            out.crossedIdx = endCrossSlot - 1;
        }

        return out;
    };

    std::optional<SearchOutcome> chosen = doSearch();

    if (!chosen) {
        result.message = "No path found (respecting one-crossing limit and ring arc length).";
        return result;
    }

    auto statePathToShortcut = [&](const std::vector<int>& sp) {
        ShortcutPath path;
        path.distance = chosen->optDist;
        for (int sid : sp) {
            int pIdx = std::get<0>(decodeStateEx(sid));
            path.vertices.emplace_back(g.xs[pIdx], g.ys[pIdx]);
        }
        // Authoritative bend count from simplified geometry (ignores sub-s_min stubs).
        path.bendCount = countBendsIgnoringShortStubs(
            path.vertices, ShortcutGrid::DEFAULT_S_MIN);
        return path;
    };

    result.success = true;
    result.crossedShortcutIdx = chosen->crossedIdx;
    result.primary = statePathToShortcut(chosen->statePaths[0]);
    for (int i = 1; i < (int)chosen->statePaths.size(); ++i)
        result.alternatives.push_back(statePathToShortcut(chosen->statePaths[i]));

    result.message = (chosen->crossedIdx >= 0)
        ? "Found path with one shortcut intersection (index "
          + std::to_string(chosen->crossedIdx) + ")."
        : "Found path with no shortcut intersections.";
    return result;
}
