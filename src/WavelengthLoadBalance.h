#pragma once

#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ShortcutMethods.h"
#include "ShortcutOrchestrator.h"

/// Directed use of one shortcut waveguide (src→dest or dest→src).
struct WlbShortcutHop {
    int shortcutIdx = -1;
    int from = -1;
    int to = -1;
};

/// One admissible route for a demand (ring arcs + optional shortcut hops).
struct WlbRoute {
    std::string label;  // "CW", "CCW", "SC#k", "CW+SC#k+CCW", ...
    std::vector<std::pair<int, int>> ringArcs;  // directed (from,to)
    std::vector<WlbShortcutHop> shortcutHops;   // directed shortcut uses
    double length = 0.0;
    int bendCount = 0;
};

struct WlbDemandAssignment {
    int demandIdx = -1;
    int src = -1;
    int dest = -1;
    WlbRoute route;
    int wavelength = 0;  // 1-indexed
    int admCount = 0;
    bool wasFixed = false;  // |Adm|=1
};

struct WavelengthLoadBalanceResult {
    bool success = false;
    std::string message;
    double Wstar = 0.0;
    int lambdaNeeded = 0;
    int conflictLoadLowerBound = 0;
    bool allRoutesRespectWstar = true;
    std::vector<WlbDemandAssignment> assignments;
};

/// Post-process after the proxy master: admissible CW/CCW/shortcut routes
/// under W*, load-greedy assignment, then DSATUR wavelength coloring.
///
/// Physical model (matches ring CW/CCW):
/// - CW and CCW are independent waveguides (opposite ring directions do not conflict).
/// - Each shortcut is two independent directional waveguides; opposite directions
///   on the same shortcut do NOT conflict for λ.
/// - Shortcut crossings conflict fully (any direction on p vs any direction on q).
WavelengthLoadBalanceResult runWavelengthLoadBalance(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const ShortcutMethodResult& loopResult,
    double Wstar);

void printWavelengthLoadBalanceReport(
    const WavelengthLoadBalanceResult& res,
    const DemandMatrix& D);
