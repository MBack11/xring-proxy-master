#include "RingDemandRouting.h"

#include <cmath>
#include <tuple>
#include <utility>

#include "PhysicalConstants.h"

void applyShortestArcDemandRouting(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        MILPSolveResult& result) {
    const int N = (int)nodes.size();
    const int tourLen = (int)result.tour.size();
    if (tourLen == 0) return;

    std::vector<int> tourPos(N, -1);
    for (int i = 0; i < tourLen; ++i)
        tourPos[result.tour[i]] = i;

    auto manhattan = [&](int i, int j) {
        return std::abs(nodes[i].x - nodes[j].x) + std::abs(nodes[i].y - nodes[j].y);
    };

    auto collectArc = [&](int ps, int pd, int step) {
        std::vector<std::pair<int, int>> edges;
        int k = ps;
        while (k != pd) {
            const int nk = (k + step + tourLen) % tourLen;
            edges.emplace_back(result.tour[k], result.tour[nk]);
            k = nk;
        }
        return edges;
    };

    auto arcMetrics = [&](const std::vector<std::pair<int, int>>& edges) {
        double dist = 0.0;
        double il = 0.0;
        int bends = 0;
        for (const auto& [i, j] : edges) {
            dist += manhattan(i, j);
            il += reportingRingEdgeIL_dB(
                nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y);
            if (ringEdgeNeedsBend(nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y))
                ++bends;
        }
        return std::tuple<double, double, int>{dist, il, bends};
    };

    const int numDemands = (int)D.demands.size();
    result.demandFlowEdges.assign(numDemands, {});
    result.demandIL.assign(numDemands, 0.0);
    result.demandDistance.assign(numDemands, 0.0);
    result.demandBendCount.assign(numDemands, 0);

    double maxW = 0.0;
    for (int q = 0; q < numDemands; ++q) {
        const int s = D.demands[q].first;
        const int d = D.demands[q].second;
        if (tourPos[s] < 0 || tourPos[d] < 0) continue;

        const int ps = tourPos[s];
        const int pd = tourPos[d];
        const auto cw = collectArc(ps, pd, +1);
        const auto ccw = collectArc(ps, pd, -1);
        const auto cwM = arcMetrics(cw);
        const auto ccwM = arcMetrics(ccw);

        const double ccwDist = std::get<0>(ccwM);
        const double cwDist = std::get<0>(cwM);
        const bool takeCcw =
            (ccwDist < cwDist - MILP_EPS) ||
            (std::abs(ccwDist - cwDist) <= MILP_EPS &&
             std::get<2>(ccwM) < std::get<2>(cwM));

        const auto& chosen = takeCcw ? ccw : cw;
        const auto& chosenM = takeCcw ? ccwM : cwM;

        result.demandFlowEdges[q] = chosen;
        result.demandDistance[q] = std::get<0>(chosenM);
        result.demandIL[q] = std::get<1>(chosenM);
        result.demandBendCount[q] = std::get<2>(chosenM);
        maxW = std::max(maxW, result.demandDistance[q]);
    }
    result.W = maxW;
}
