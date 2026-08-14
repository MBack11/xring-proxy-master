#include "RingGeometry.h"

#include <cmath>

RingGeometry::RingGeometry(const std::vector<Node>& nodes)
    : nodes(nodes), N((int)nodes.size()) {
    computeOnce();
}

double RingGeometry::manhattanDist(int i, int j) const {
    return std::abs(nodes[i].x - nodes[j].x)
         + std::abs(nodes[i].y - nodes[j].y);
}

std::vector<RingSegment> RingGeometry::getSegments(int i, int j, int opt) const {
    const double xi = nodes[i].x, yi = nodes[i].y;
    const double xj = nodes[j].x, yj = nodes[j].y;
    std::vector<RingSegment> segs;
    if (opt == 0) {
        segs.push_back({true, yi, std::min(xi, xj), std::max(xi, xj)});
        segs.push_back({false, xj, std::min(yi, yj), std::max(yi, yj)});
    } else {
        segs.push_back({false, xi, std::min(yi, yj), std::max(yi, yj)});
        segs.push_back({true, yj, std::min(xi, xj), std::max(xi, xj)});
    }
    return segs;
}

bool RingGeometry::segmentsCross(const RingSegment& s1, const RingSegment& s2) const {
    if (s1.isHorizontal != s2.isHorizontal) {
        const RingSegment& h = s1.isHorizontal ? s1 : s2;
        const RingSegment& v = s1.isHorizontal ? s2 : s1;
        return (h.lo < v.fixed && v.fixed < h.hi)
            && (v.lo < h.fixed && h.fixed < v.hi);
    }
    if (std::abs(s1.fixed - s2.fixed) > 1e-9) return false;
    return s1.lo < s2.hi && s2.lo < s1.hi;
}

bool RingGeometry::segmentsTouchPerpendicular(
        const RingSegment& s1,
        const RingSegment& s2) const {
    if (s1.isHorizontal == s2.isHorizontal) return false;
    const RingSegment& h = s1.isHorizontal ? s1 : s2;
    const RingSegment& v = s1.isHorizontal ? s2 : s1;
    constexpr double eps = 1e-9;
    return (h.lo <= v.fixed + eps && v.fixed <= h.hi + eps)
        && (v.lo <= h.fixed + eps && h.fixed <= v.hi + eps);
}

bool RingGeometry::edgesShareEndpoint(int i, int j, int k, int l) const {
    return i == k || i == l || j == k || j == l;
}

bool RingGeometry::collinearOverlapBeyondShared(
        int i,
        int j,
        int opt1,
        int k,
        int l,
        int opt2,
        int sharedNode) const {
    const double px = nodes[sharedNode].x;
    const double py = nodes[sharedNode].y;
    constexpr double eps = 1e-9;

    const auto segs1 = getSegments(i, j, opt1);
    const auto segs2 = getSegments(k, l, opt2);
    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (s1.isHorizontal != s2.isHorizontal) continue;
            if (std::abs(s1.fixed - s2.fixed) > eps) continue;

            const double olLo = std::max(s1.lo, s2.lo);
            const double olHi = std::min(s1.hi, s2.hi);
            if (olHi - olLo <= eps) continue;

            const double sharedCoord = s1.isHorizontal ? px : py;
            if (olHi > sharedCoord + eps || olLo < sharedCoord - eps)
                return true;
        }
    }
    return false;
}

bool RingGeometry::doesCross(
        int i,
        int j,
        int opt1,
        int k,
        int l,
        int opt2) const {
    const auto segs1 = getSegments(i, j, opt1);
    const auto segs2 = getSegments(k, l, opt2);
    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (segmentsCross(s1, s2)) return true;
        }
    }

    if (edgesShareEndpoint(i, j, k, l)) {
        const int sharedNodes[] = {i, j, k, l};
        for (int s : sharedNodes) {
            const bool isShared = (s == i || s == j) && (s == k || s == l);
            if (!isShared) continue;
            if (collinearOverlapBeyondShared(i, j, opt1, k, l, opt2, s))
                return true;
        }
        return false;
    }

    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (segmentsTouchPerpendicular(s1, s2)) return true;
        }
    }
    return false;
}

bool RingGeometry::pathPassesThroughNode(int i, int j, int opt, int node) const {
    if (node == i || node == j) return false;

    const double xk = nodes[node].x;
    const double yk = nodes[node].y;
    const double xi = nodes[i].x, yi = nodes[i].y;
    const double xj = nodes[j].x, yj = nodes[j].y;
    constexpr double eps = 1e-9;

    auto onClosedSegment = [&](const RingSegment& s) {
        if (s.isHorizontal) {
            if (std::abs(yk - s.fixed) > eps) return false;
            return xk >= s.lo - eps && xk <= s.hi + eps;
        }
        if (std::abs(xk - s.fixed) > eps) return false;
        return yk >= s.lo - eps && yk <= s.hi + eps;
    };

    for (const RingSegment& s : getSegments(i, j, opt)) {
        if (!onClosedSegment(s)) continue;
        const bool atI = std::abs(xk - xi) < eps && std::abs(yk - yi) < eps;
        const bool atJ = std::abs(xk - xj) < eps && std::abs(yk - yj) < eps;
        if (!atI && !atJ) return true;
    }
    return false;
}

void RingGeometry::computeOnce() {
    if (computed_) return;
    computed_ = true;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            for (int k = 0; k < N; ++k)
                for (int l = 0; l < N; ++l) {
                    if (k == l) continue;
                    if (edgeIdx(i, j) >= edgeIdx(k, l)) continue;
                    for (int opt1 = 0; opt1 <= 1; ++opt1)
                        for (int opt2 = 0; opt2 <= 1; ++opt2) {
                            if (doesCross(i, j, opt1, k, l, opt2))
                                cachedCrossingPairs_.push_back(
                                    {i, j, opt1, k, l, opt2});
                        }
                }
        }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            for (int opt = 0; opt <= 1; ++opt)
                for (int node = 0; node < N; ++node) {
                    if (pathPassesThroughNode(i, j, opt, node))
                        cachedPassThroughTriples_.push_back({i, j, opt, node});
                }
        }
}
