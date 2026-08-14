//
// Created by Mika Back on 02.06.26.
//

#pragma once

#include <set>
#include <string>
#include <vector>
#include "gurobi_c++.h"
#include "Nodes.h"
#include "DemandMatrix.h"

struct RingSegment {
    bool isHorizontal;
    double fixed;
    double lo, hi;
};

struct MILPSolveOptions {
    /// Omit static C5 cuts and add crossing constraints lazily at incumbents.
    bool lazyC5 = false;
    /// Gurobi MIPFocus; negative = use solve()-mode defaults.
    int mipFocus = -1;
    /// Gurobi ImproveStartTime (seconds); negative = use solve()-mode defaults.
    double improveStartTime = -1.0;
    /// If true, omit b[0,1]==1 symmetry break (full rotation-invariant search).
    bool skipSymmetryBreak = false;
};

struct MILPSolveResult {
    bool success = false;
    int status = 0;
    double W = 0.0;
    std::vector<int> tour;
    struct TourEdge {
        int from;
        int to;
        int routing;  // 0 = HV, 1 = VH
    };
    std::vector<TourEdge> tourEdges;
    std::vector<double> demandIL;
    std::vector<double> demandDistance;
    std::vector<int> demandBendCount;
    std::vector<std::vector<std::pair<int, int>>> demandFlowEdges;
};

class MilpLazyCallback;

class MILPSolver {
public:
    MILPSolver(const std::vector<Node>& nodes, const DemandMatrix& D);

    MILPSolveResult solve(
        bool optimizeWorstCase,
        const std::string& fileSuffix,
        const std::set<int>& excludedDemands = {},
        double bestObjStop = -1.0,
        const MILPSolveResult* warmStart = nullptr,
        double timeLimitSeconds = -1.0,
        bool quiet = false,
        const MILPSolveOptions& options = {});

    void applyShortestArcDemandRouting(
        MILPSolveResult& result,
        const std::set<int>& excludedDemands = {}) const;

private:
    friend class MilpLazyCallback;

    struct CrossingPair {
        int i, j, opt1, k, l, opt2;
    };

    struct PassThroughTriple {
        int i, j, opt, node;
    };

    const std::vector<Node>& nodes;
    const DemandMatrix& D;
    int N;

    std::vector<CrossingPair> cachedCrossingPairs;
    std::vector<PassThroughTriple> cachedPassThroughTriples;
    bool geometryConstraintsComputed = false;

    void computeGeometryConstraintsOnce();

    double manhattanDist(int i, int j) const;
    double reportingDemandEdgeIL(int i, int j) const;
    double insertionLoss(int i, int j) const;
    inline int edge_idx(int i, int j) const { return i * N + j; }

    std::vector<RingSegment> getSegments(int i, int j, int opt) const;
    bool segmentsCross(const RingSegment& s1, const RingSegment& s2) const;
    bool segmentsTouchPerpendicular(const RingSegment& s1, const RingSegment& s2) const;
    bool edgesShareEndpoint(int i, int j, int k, int l) const;
    bool collinearOverlapBeyondShared(
        int i, int j, int opt1, int k, int l, int opt2, int sharedNode) const;
    bool doesCross(int i, int j, int opt1, int k, int l, int opt2) const;
    bool pathPassesThroughNode(int i, int j, int opt, int node) const;

    double computeBigM() const;

    void applyWarmStartPositionModel(
        const MILPSolveResult& warmStart,
        const std::set<int>& excludedDemands,
        std::vector<GRBVar>& pos,
        std::vector<GRBVar>& wrapVar,
        std::vector<GRBVar>& yVar,
        std::vector<GRBVar>& forwardVar,
        std::vector<GRBVar>& demandDistVar,
        std::vector<GRBVar>& dirVar);

    void verifyPositionModelSolution(
        const MILPSolveResult& result,
        const std::vector<GRBVar>& demandDistVar,
        const std::vector<GRBVar>& forwardVar,
        const std::vector<GRBVar>& pos,
        double L_total_val,
        const std::set<int>& excludedDemands) const;

    void selfCheckRingCrossings(
        const std::vector<GRBVar>& b,
        const std::vector<GRBVar>& r,
        const std::string& context) const;

    void exportCSV(const MILPSolveResult& result,
                   const std::string& fileSuffix,
                   const std::set<int>& excludedDemands) const;

    GRBEnv env;
};
