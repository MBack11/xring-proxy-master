#pragma once

#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ShortcutTypes.h"
#include "ShortcutGrid.h"

enum class GridDir : int {
    NONE  = 0,
    RIGHT = 1,
    LEFT  = 2,
    UP    = 3,
    DOWN  = 4
};

struct ShortcutPath {
    std::vector<std::pair<double, double>> vertices;
    double distance = 0.0;
    int bendCount = 0;
};

struct ShortcutRouteResult {
    bool success = false;
    std::string message;
    ShortcutPath primary;
    std::vector<ShortcutPath> alternatives;
    int crossedShortcutIdx = -1;  // index into existingShortcuts, if one crossing used
};

class ShortcutRouter {
public:
  static constexpr int MAX_ALTERNATIVES = 10;

  // True if two shortcut paths share any geometric intersection.
  static bool shortcutsIntersect(const Shortcut& a, const Shortcut& b);
  // How many other shortcuts shortcut at index i intersects.
  static int countShortcutIntersections(
      const std::vector<Shortcut>& shortcuts, int index);

  static ShortcutRouteResult findPaths(
      const std::vector<GridPoint>& gridPoints,
      int srcId,
      int destId,
      const std::vector<Node>& nodes,
      const std::vector<Segment>& ringSegments,
      const std::vector<int>& tour,
      const std::vector<EdgeOption>& routing,
      const std::vector<Shortcut>& existingShortcuts = {},
      double maxDistance = std::numeric_limits<double>::infinity());

  static bool segmentsGeometricallyIntersect(const Segment& a, const Segment& b);

private:
    struct RingSegmentIndex {
        double cellSize = 1.0;
        const std::vector<Segment>* segments = nullptr;
        std::unordered_map<long long, std::vector<int>> buckets;

        void build(const std::vector<Segment>& ring, double cell);
        void forEachCandidate(const Segment& edge, const std::function<void(int)>& fn) const;
    };

    static constexpr double EPS = 1e-9;
    static constexpr int NUM_DIRS = 5;

    struct Graph {
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<bool> usable;
        std::vector<std::vector<int>> adj;
    };

    static int stateId(int pointIdx, GridDir dir);
    static std::pair<int, GridDir> decodeState(int sid);
    static GridDir directionBetween(double x1, double y1, double x2, double y2);
    static bool lexLess(double d1, int b1, double d2, int b2);
    static bool lexEq(double d1, int b1, double d2, int b2);

    static bool isParallelOverlap(const Segment& s1, const Segment& s2);
    static bool isStrictCrossing(const Segment& s1, const Segment& s2);
    static bool isTJunction(double px, double py, const Segment& obs, bool gridIsH);
    static bool edgeIntersectsShortcut(const Segment& edge, const Shortcut& sc);
    static bool edgeIsOnExemptSegment(
        const Segment& edge,
        const std::vector<Segment>& exemptSegments);
    static bool gridEdgeAllowed(
        double x1, double y1, double x2, double y2,
        const std::vector<Segment>& ringSegments,
        const std::vector<Segment>& exemptSegments,
        const std::vector<Segment>& hardShortcutSegments,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes);

    struct EdgeCrossInfo {
        bool hardBlocked = true;
        int softShortcutIdx = -1;  // -1 none, -2 multiple, >=0 into existingShortcuts
    };
    static EdgeCrossInfo classifyGridEdge(
        double x1, double y1, double x2, double y2,
        const std::vector<Segment>& ringSegments,
        const std::vector<Segment>& exemptSegments,
        const std::vector<Segment>& hardShortcutSegments,
        const std::vector<std::pair<int, std::vector<Segment>>>& softShortcuts,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes,
        const RingSegmentIndex* ringIndex = nullptr);

    static int stateId(int pointIdx, GridDir dir, int crossUsed);
    static std::tuple<int, GridDir, int> decodeStateEx(int sid);

    static bool edgePassesThroughTourNode(
        double x1, double y1, double x2, double y2,
        const std::vector<int>& tour,
        int srcId,
        int destId,
        const std::vector<Node>& nodes);

    static std::optional<Graph> buildGraph(const std::vector<GridPoint>& gridPoints);
    static int findNodeIndex(const Graph& g, double x, double y);
    static void connectLine(
        Graph& g,
        const std::vector<int>& indices,
        bool sortByY);
    static void enumeratePaths(
        int srcState,
        const std::vector<std::vector<int>>& preds,
        const std::vector<double>& dist,
        const std::vector<int>& bends,
        double optDist,
        int optBends,
        int destPoint,
        int maxPaths,
        std::vector<std::vector<int>>& outStatePaths);
};
