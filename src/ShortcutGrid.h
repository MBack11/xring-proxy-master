#pragma once

#include <set>
#include <string>
#include <vector>

#include "Nodes.h"
#include "RingTypes.h"

struct GridPoint {
    double x;
    double y;
    std::string source;  // "main_grid", "node_line", "inner_ring"
    bool usable;         // per (src, dest) search
};

struct ParallelGapInfo {
    double d = 0.0;
    bool horizontalGap = true;  // true: measured in y between horizontal edges
    Segment edgeA{};
    Segment edgeB{};
    double anchorX = 0.0;
    double anchorY = 0.0;
};

struct ShortcutGridResult {
    bool success = false;
    std::string message;
    double d = 0.0;
    double sMin = 0.0;
    bool preconditionOk = false;
    ParallelGapInfo gapInfo{};
    std::vector<GridPoint> points;
    double gridStep = 0.0;
    int xsCount = 0;
    int ysCount = 0;
};

class ShortcutGrid {
public:
    static constexpr double DEFAULT_S_MIN = 0.2;
    static constexpr int MAX_GRID_LINES = 30;

    static ShortcutGridResult build(
        const std::vector<Node>& nodes,
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        int srcId,
        int destId,
        double sMin = DEFAULT_S_MIN);

    static ShortcutGridResult buildCached(
        const std::vector<Node>& nodes,
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        int srcId,
        int destId,
        double sMin = DEFAULT_S_MIN);

    struct BuildCacheStats {
        size_t hits = 0;
        size_t misses = 0;
        size_t entries = 0;
    };
    static BuildCacheStats buildCacheStats();
    static void resetBuildCacheStats();

    static std::vector<Segment> collectRingSegments(
        const std::vector<EdgeOption>& routing);

    // First segment of outgoing edge + last segment of incoming edge (ring leg at node).
    static std::vector<Segment> incidentExemptSegmentsForNode(
        int nodeId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing);

    static bool computeParallelGap(
        const std::vector<Segment>& ringSegments,
        const std::vector<Node>& nodes,
        ParallelGapInfo& out);

private:
    static constexpr double EPS = 1e-9;

    static bool isHorizontal(const Segment& s);
    static bool isVertical(const Segment& s);
    static bool pointOnSegment(double px, double py, const Segment& s);
    static std::vector<double> spacedCoords(double anchor, double lo, double hi, double step);
    static void mergeCoord(std::vector<double>& coords, double v, double minSep);
    static void collectGlobalClearanceLines(
        const std::vector<Segment>& ringSegments,
        const std::vector<int>& tour,
        const std::vector<Node>& nodes,
        double sMin,
        std::set<double>& clearanceXs,
        std::set<double>& clearanceYs);
    static void addPoint(
        std::vector<GridPoint>& pts,
        double x,
        double y,
        const std::string& source,
        bool usable);
    static bool classifyUsability(
        double x,
        double y,
        const std::vector<Segment>& ringSegments,
        int srcId,
        int destId,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        const std::vector<Node>& nodes);

    static bool isOnTourNode(
        double x,
        double y,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes);
};
