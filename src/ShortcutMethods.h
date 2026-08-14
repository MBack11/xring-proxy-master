#pragma once

#include <limits>
#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ShortcutOrchestrator.h"
#include "ShortcutGrid.h"

struct ResolveAttemptRecord {
    int iteration = 0;
    std::string outcome;
    double wBefore = 0.0;
    double wAfter = 0.0;
};

struct ShortcutMethodProfile {
    double shortcutSearchSec = 0.0;
    double resolveSec = 0.0;
    int wImprovementEvents = 0;
    std::vector<ResolveAttemptRecord> resolveAttempts;
};

struct ShortcutMethodOptions {
    bool quiet = false;
    bool skipExports = false;
    double resolveTimeLimitSec = 20.0;
    /// Total wall-clock budget for this method run (shortcut loop). <0 = unlimited.
    double wallTimeLimitSec = -1.0;
    ShortcutUsageMode usageMode = ShortcutUsageMode::Shared;
    ShortcutMethodProfile* profile = nullptr;
    /// Dual-safe early abort for Method D Phase 2 (quality-preserving):
    /// stop when ObjBound >= abortIfObjBoundGe AND ObjBound > abortOnlyIfObjBoundGt.
    /// Default: disabled (abortIfObjBoundGe = +inf).
    double abortIfObjBoundGe = std::numeric_limits<double>::infinity();
    double abortOnlyIfObjBoundGt = -std::numeric_limits<double>::infinity();
};

enum class ShortcutTermination {
    Case2Constraints,
    Case3NoFurtherImprovement,
    WallTimeLimit
};

struct ShortcutMethodResult {
    MILPSolveResult layout;
    std::vector<PlacedShortcut> shortcuts;
    double globalW = 0.0;
    ShortcutTermination termination = ShortcutTermination::Case2Constraints;
};

ShortcutMethodResult runMethodAWithShortcuts(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const MILPSolveResult& layoutFixed,
    double sMin = ShortcutGrid::DEFAULT_S_MIN,
    const ShortcutMethodOptions& options = {});

ShortcutMethodResult runMethodBWithShortcuts(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    MILPSolver& solver,
    const MILPSolveResult& initialLayout,
    double sMin = ShortcutGrid::DEFAULT_S_MIN,
    const std::string& exportSuffix = "B",
    const char* logLabel = "Method B+Shortcuts",
    const ShortcutMethodOptions& options = {});

const char* terminationToString(ShortcutTermination reason);
