#pragma once

#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"

void applyShortestArcDemandRouting(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    MILPSolveResult& result);
