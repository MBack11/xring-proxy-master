#pragma once

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "DemandMatrix.h"
#include "LShapeCrossing.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ShortcutOrchestrator.h"

/// No-good cut on an evaluated ring (directed tour arcs as found).
struct MasterNogoodCut {
    std::vector<std::pair<int, int>> directedArcs;
};

struct MasterSolution {
    MILPSolveResult layout;
    std::vector<int> selectedShortcutIndices;
    double Wproxy = 0.0;
};

struct ProxyMasterMilpOptions {
    double timeLimitSec = 120.0;
    int poolSolutions = 10;
    bool quiet = false;
    /// If >= 0, set Gurobi MIPGap; otherwise use Gurobi default.
    double mipGap = -1.0;
    /// Shared: any demand may flow on selected z[c]. Exclusive: only owner
    /// demand with endpoints {s_q,t_q} matching pair c.
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared;
};

struct ProxyMasterMilpResult {
    bool success = false;
    double Wproxy = 0.0;
    std::vector<MasterSolution> pool;
    int numXPairs = 0;
    std::string message;
    int solverStatus = 0;
    double objBound = 0.0;
    double mipGapAtStop = 0.0;
    double runtimeSec = 0.0;
};

/// Proxy optimum with ring fixed (only z + flow + W).
struct FixedRingProxyResult {
    bool success = false;
    double Wproxy = 0.0;
    std::vector<int> selectedShortcutIndices;
    std::string message;
};

/// Persistent master MILP: build once, append nogood cuts, reoptimize.
class ProxyMasterMilpSession {
public:
    ProxyMasterMilpSession(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        ShortcutUsageMode usageMode = ShortcutUsageMode::Shared);
    ~ProxyMasterMilpSession();
    ProxyMasterMilpSession(const ProxyMasterMilpSession&) = delete;
    ProxyMasterMilpSession& operator=(const ProxyMasterMilpSession&) = delete;

    ProxyMasterMilpResult solve(
        const std::vector<MasterNogoodCut>& nogoodCuts,
        const ProxyMasterMilpOptions& options = {});

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// One-shot wrapper (builds a fresh session). Prefer ProxyMasterMilpSession in loops.
ProxyMasterMilpResult solveProxyMasterMilp(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const std::vector<MasterNogoodCut>& nogoodCuts = {},
    const ProxyMasterMilpOptions& options = {});

/// Per-ring proxy LB: optimize abstract shortcuts + flow on a fixed layout.
FixedRingProxyResult solveProxyMasterFixedRing(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const MILPSolveResult& layoutFixed,
    double timeLimitSec = 60.0,
    bool quiet = true,
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared);

/// Undirected edge signature for deduplicating evaluated rings.
std::string ringTourKey(const MILPSolveResult& layout);

/// Undirected edge list of a ring layout.
std::vector<std::pair<int, int>> ringUndirectedEdges(const MILPSolveResult& layout);
