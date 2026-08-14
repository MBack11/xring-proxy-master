#pragma once

#include <utility>
#include <vector>

#include "Nodes.h"
#include "ShortcutTypes.h"

/// One master shortcut candidate p = (i, j), i < j.
struct ShortcutPairIndex {
    int i = 0;
    int j = 0;
    double delta = 0.0;
};

/// All unordered node pairs (master candidate set P).
std::vector<ShortcutPairIndex> buildMasterShortcutPairs(
    const std::vector<Node>& nodes);

/// Build an L-shape shortcut at Manhattan length δ_ij.
/// vhFirst=true → VH (vertical then horizontal); false → HV.
Shortcut makeLShapeShortcut(
    int i,
    int j,
    bool vhFirst,
    const std::vector<Node>& nodes);

/// True iff all four {VH,HV}×{VH,HV} L-shape combinations cross.
bool lShapePairAlwaysCrosses(
    int i,
    int j,
    int k,
    int l,
    const std::vector<Node>& nodes);

/// XPairs: indices (p, q) with p < q into `pairs`, where lShapePairAlwaysCrosses.
std::vector<std::pair<int, int>> computeShortcutCrossingPairs(
    const std::vector<Node>& nodes,
    const std::vector<ShortcutPairIndex>& pairs);
