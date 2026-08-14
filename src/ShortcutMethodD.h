#pragma once

#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ShortcutGrid.h"
#include "ShortcutMethods.h"

/// Method D — JointShortcuts: global shortcut selection via MILP on a fixed ring.
/// Phase 1 routes all non-adjacent node pairs in isolation; Phase 2 selects a
/// feasible subset minimizing worst-case demand distance (mm). Crossing rule
/// is C2b (at most one soft crossing per selected shortcut).
ShortcutMethodResult runMethodDJointShortcuts(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const MILPSolveResult& layoutFixed,
    double sMin = ShortcutGrid::DEFAULT_S_MIN,
    const ShortcutMethodOptions& options = {});
