#pragma once

#include <vector>

#include "MILPSolver.h"
#include "RingTypes.h"

struct RingLayout {
    std::vector<Edge> edges;
    std::vector<int> tour;
    std::vector<EdgeOption> routing;
};

RingLayout buildRingLayout(
    const std::vector<Node>& nodes,
    const MILPSolveResult& layout);

std::vector<Segment> pathVerticesToSegments(
    const std::vector<std::pair<double, double>>& vertices);
