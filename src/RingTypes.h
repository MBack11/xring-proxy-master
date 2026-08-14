#pragma once

#include <vector>

#include "Nodes.h"

struct Segment {
    double x1, y1, x2, y2;
};

struct RoutingOption {
    std::vector<Segment> segments;
};

struct Edge {
    int from;
    int to;
    double distance;
    RoutingOption option1;
    RoutingOption option2;
};

struct EdgeOption {
    const Edge* edge;
    int optionIndex;  // 1 = HV, 2 = VH
    const RoutingOption* option;
};
