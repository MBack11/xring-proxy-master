#pragma once

#include <vector>

#include "RingTypes.h"

struct Shortcut {
    int from = 0;
    int to = 0;
    double approx_length = 0.0;
    int bend_count = 0;
    bool everCrossed = false;
    int demandIdx = -1;
    std::vector<Segment> path;
};
