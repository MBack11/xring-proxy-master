#pragma once

#include <vector>

#include "MILPSolver.h"
#include "Nodes.h"

struct RingCrossingPair {
    int i = 0, j = 0, opt1 = 0, k = 0, l = 0, opt2 = 0;
};

struct RingPassThroughTriple {
    int i = 0, j = 0, opt = 0, node = 0;
};

/// Ring-edge geometry for C5 / C5b (from Idea-1 MILPSolver, ring-only).
class RingGeometry {
public:
    explicit RingGeometry(const std::vector<Node>& nodes);

    int n() const { return N; }
    int edgeIdx(int i, int j) const { return i * N + j; }

    double manhattanDist(int i, int j) const;

    const std::vector<RingCrossingPair>& crossingPairs() const {
        return cachedCrossingPairs_;
    }
    const std::vector<RingPassThroughTriple>& passThroughTriples() const {
        return cachedPassThroughTriples_;
    }

    bool doesCross(int i, int j, int opt1, int k, int l, int opt2) const;
    std::vector<RingSegment> getSegments(int i, int j, int opt) const;

private:
    const std::vector<Node>& nodes;
    int N = 0;

    std::vector<RingCrossingPair> cachedCrossingPairs_;
    std::vector<RingPassThroughTriple> cachedPassThroughTriples_;
    bool computed_ = false;

    bool segmentsCross(const RingSegment& s1, const RingSegment& s2) const;
    bool segmentsTouchPerpendicular(const RingSegment& s1, const RingSegment& s2) const;
    bool edgesShareEndpoint(int i, int j, int k, int l) const;
    bool collinearOverlapBeyondShared(
        int i, int j, int opt1, int k, int l, int opt2, int sharedNode) const;
    bool pathPassesThroughNode(int i, int j, int opt, int node) const;
    void computeOnce();
};
