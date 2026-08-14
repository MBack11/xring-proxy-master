#include "LShapeCrossing.h"

#include <cmath>

#include "ShortcutRouter.h"

namespace {

double manhattan(const Node& a, const Node& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

}  // namespace

std::vector<ShortcutPairIndex> buildMasterShortcutPairs(
        const std::vector<Node>& nodes) {
    const int N = (int)nodes.size();
    std::vector<ShortcutPairIndex> pairs;
    pairs.reserve(N * (N - 1) / 2);
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            pairs.push_back({i, j, manhattan(nodes[i], nodes[j])});
    return pairs;
}

Shortcut makeLShapeShortcut(
        int i,
        int j,
        bool vhFirst,
        const std::vector<Node>& nodes) {
    Shortcut sc;
    sc.from = i;
    sc.to = j;
    sc.approx_length = manhattan(nodes[i], nodes[j]);
    sc.bend_count = 1;

    const double xi = nodes[i].x, yi = nodes[i].y;
    const double xj = nodes[j].x, yj = nodes[j].y;

    if (vhFirst) {
        sc.path.push_back({xi, yi, xi, yj});
        sc.path.push_back({xi, yj, xj, yj});
    } else {
        sc.path.push_back({xi, yi, xj, yi});
        sc.path.push_back({xj, yi, xj, yj});
    }
    return sc;
}

bool lShapePairAlwaysCrosses(
        int i,
        int j,
        int k,
        int l,
        const std::vector<Node>& nodes) {
    const bool opts[2] = {false, true};  // HV, VH
    for (bool pOpt : opts) {
        const Shortcut sp = makeLShapeShortcut(i, j, pOpt, nodes);
        for (bool qOpt : opts) {
            const Shortcut sq = makeLShapeShortcut(k, l, qOpt, nodes);
            if (!ShortcutRouter::shortcutsIntersect(sp, sq))
                return false;
        }
    }
    return true;
}

std::vector<std::pair<int, int>> computeShortcutCrossingPairs(
        const std::vector<Node>& nodes,
        const std::vector<ShortcutPairIndex>& pairs) {
    std::vector<std::pair<int, int>> xPairs;
    const int C = (int)pairs.size();
    for (int p = 0; p < C; ++p) {
        const int i = pairs[p].i, j = pairs[p].j;
        for (int q = p + 1; q < C; ++q) {
            const int k = pairs[q].i, l = pairs[q].j;
            if (lShapePairAlwaysCrosses(i, j, k, l, nodes))
                xPairs.emplace_back(p, q);
        }
    }
    return xPairs;
}
