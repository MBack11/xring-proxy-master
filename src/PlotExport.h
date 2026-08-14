#pragma once

#include <string>

#include "DemandMatrix.h"
#include "Nodes.h"
#include "ShortcutMethods.h"
#include "WavelengthLoadBalance.h"

/// Write nodes/ring/demands/shortcuts/wavelengths CSVs for visualize.py.
/// suffix typically "PM" (proxy master) or "D".
void exportProxyMasterPlotCsvs(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const ShortcutMethodResult& result,
    const WavelengthLoadBalanceResult* wlb,
    const std::string& suffix);
