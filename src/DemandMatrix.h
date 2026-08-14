#pragma once

#include <vector>
#include <utility>

struct DemandMatrix {
    std::vector<std::pair<int,int>> demands;

    void add(int sender, int receiver) {
        demands.push_back({sender, receiver});
    }
};
