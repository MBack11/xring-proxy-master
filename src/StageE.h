#pragma once

#include <vector>

#include "DemandMatrix.h"
#include "LShapeCrossing.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ProxyMasterMilp.h"
#include "ShortcutOrchestrator.h"

/// Stage E fast realizability estimate W_hat(R) — ranking only, not a bound.
double computeStageE(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const MILPSolveResult& layout,
    const std::vector<int>& selectedShortcutIndices,
    const std::vector<ShortcutPairIndex>& pairs,
    double sMin,
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared);
