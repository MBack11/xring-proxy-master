#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "Nodes.h"
#include "ProxyMasterMilp.h"
#include "ShortcutMethodD.h"
#include "ShortcutMethods.h"

struct ProxyMasterLoopOptions {
    double eps = 1e-3;
    int maxRounds = 50;
    int poolSolutions = 10;
    double masterTimeLimitSec = 120.0;
    double wallTimeLimitSec = 600.0;
    double sMin = 0.2;
    bool quiet = false;
    /// Per-round trace for plateau diagnosis (does not change loop logic).
    bool traceRounds = false;
    /// Keep one Gurobi master model and append nogood cuts (same semantics).
    bool incrementalMaster = true;
    /// Shared (default) vs private start→end-only shortcut flow in M / E / V.
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared;
};

struct ProxyMasterLoopStageVTrace {
    std::string ringKey;
    double wTrue = std::numeric_limits<double>::quiet_NaN();
    bool fineTestFired = false;
};

struct ProxyMasterLoopRoundTrace {
    int round = 0;
    double LBk = 0.0;
    double LBnext = std::numeric_limits<double>::quiet_NaN();
    int poolSize = 0;
    int poolNew = 0;
    int poolDuplicate = 0;
    int poolResurfacedCut = 0;
    int cutsAdded = 0;
    bool stageVInvoked = false;
    std::string stageVRingKey;
    double stageVWtrue = std::numeric_limits<double>::quiet_NaN();
    std::vector<ProxyMasterLoopStageVTrace> stageVEvals;
    double incumbentBefore = std::numeric_limits<double>::infinity();
    double incumbentAfter = std::numeric_limits<double>::infinity();
    bool incumbentChanged = false;
    std::vector<std::string> poolKeys;
    std::vector<std::string> resurfacedKeys;
    bool fineStop = false;
    bool coarseStop = false;
};

struct ProxyMasterLoopResult {
    double Wstar = std::numeric_limits<double>::infinity();
    double LB = 0.0;
    double gap = std::numeric_limits<double>::infinity();
    bool provenOptimal = false;
    bool stoppedFine = false;
    bool stoppedCoarse = false;
    bool stoppedBudget = false;
    int rounds = 0;
    /// 1-based outer-loop round when W* was last improved (0 = never).
    int bestFoundRound = 0;
    int ringsEvaluated = 0;
    std::string poolModeNote;
    ShortcutMethodResult best;
    std::string stopReason;
    std::vector<ProxyMasterLoopRoundTrace> roundTraces;
};

/// Full M → E → V → C loop (Steps 3–4).
ProxyMasterLoopResult runProxyMasterLoop(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const ProxyMasterLoopOptions& options = {});

/// Step 1 checkpoint helper: Method D on baseline ring (B if ok, else A).
double runMethodDBaselineW(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    double sMin,
    bool quiet = true,
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared);
