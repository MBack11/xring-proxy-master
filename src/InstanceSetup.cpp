#include "InstanceSetup.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <sstream>

std::vector<Node> generateNodes(int n, unsigned seed, double spacing) {
    const int COLS = 20;
    const int ROWS = 20;

    std::vector<std::pair<double, double>> cells;
    cells.reserve(COLS * ROWS);
    for (int r = 0; r < ROWS; ++r)
        for (int c = 0; c < COLS; ++c) {
            const double x = c * spacing + (r % 2) * (spacing * 0.5);
            const double y = r * spacing;
            cells.push_back({x, y});
        }

    std::mt19937 rng(seed);
    std::shuffle(cells.begin(), cells.end(), rng);

    std::vector<Node> nodes;
    nodes.reserve(n);
    for (int i = 0; i < n; ++i)
        nodes.push_back({i, cells[i].first, cells[i].second});
    return nodes;
}

std::string formatTour(const std::vector<int>& tour) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < tour.size(); ++i) {
        if (i) oss << ",";
        oss << tour[i];
    }
    oss << "]";
    return oss.str();
}

double theoreticalLowerBoundW(
        const std::vector<Node>& nodes,
        const DemandMatrix& D) {
    double maxDirect = 0.0;
    for (const auto& [s, d] : D.demands) {
        const double dist = std::abs(nodes[s].x - nodes[d].x)
                          + std::abs(nodes[s].y - nodes[d].y);
        maxDirect = std::max(maxDirect, dist);
    }
    return maxDirect;
}
