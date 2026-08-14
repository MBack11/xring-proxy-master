#pragma once

#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "Nodes.h"

std::vector<Node> generateNodes(int n, unsigned seed, double spacing);

std::string formatTour(const std::vector<int>& tour);

double theoreticalLowerBoundW(
    const std::vector<Node>& nodes,
    const DemandMatrix& D);
