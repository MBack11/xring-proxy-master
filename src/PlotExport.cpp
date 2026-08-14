#include "PlotExport.h"

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "PhysicalConstants.h"
#include "ShortcutOrchestrator.h"

namespace {

void writeRingCsv(
        const MILPSolveResult& layout,
        const std::vector<Node>& nodes,
        const std::string& filename) {
    std::ofstream f(filename);
    f << "from_id,to_id,from_x,from_y,to_x,to_y,routing\n";
    for (const auto& te : layout.tourEdges) {
        f << te.from + 1 << "," << te.to + 1 << ","
          << nodes[te.from].x << "," << nodes[te.from].y << ","
          << nodes[te.to].x << "," << nodes[te.to].y << ","
          << (te.routing == 0 ? "HV" : "VH") << "\n";
    }
}

void writeDemandsCsvFromWlb(
        const WavelengthLoadBalanceResult& wlb,
        const DemandMatrix& D,
        const std::string& filename) {
    std::ofstream f(filename);
    f << "demand,sender,receiver,from_id,to_id,direction\n";
    for (const auto& a : wlb.assignments) {
        const int s = a.src;
        const int t = a.dest;
        for (const auto& [u, v] : a.route.ringArcs) {
            f << a.demandIdx << "," << s + 1 << "," << t + 1 << ","
              << u + 1 << "," << v + 1 << ",flow\n";
        }
    }
}

void writeDemandsCsv(
        const MILPSolveResult& layout,
        const DemandMatrix& D,
        const std::string& filename) {
    std::ofstream f(filename);
    f << "demand,sender,receiver,from_id,to_id,direction\n";
    for (int q = 0; q < (int)D.demands.size(); ++q) {
        const int s = D.demands[q].first;
        const int t = D.demands[q].second;
        if (q >= (int)layout.demandFlowEdges.size()) continue;
        for (const auto& [u, v] : layout.demandFlowEdges[q]) {
            f << q << "," << s + 1 << "," << t + 1 << ","
              << u + 1 << "," << v + 1 << ",flow\n";
        }
    }
}

void writeShortcutsCsv(
        const std::vector<PlacedShortcut>& shortcuts,
        const std::string& filename) {
    std::ofstream f(filename);
    f << "shortcut_id,from_id,to_id,leg_index,leg_type,"
         "phys_x1,phys_y1,phys_x2,phys_y2,"
         "logical_x1,logical_y1,logical_x2,logical_y2,"
         "approx_length,bend_count,everCrossed\n";
    for (int si = 0; si < (int)shortcuts.size(); ++si) {
        const PlacedShortcut& sc = shortcuts[si];
        const auto& verts = sc.path.vertices;
        for (size_t leg = 0; leg + 1 < verts.size(); ++leg) {
            const double x1 = verts[leg].first;
            const double y1 = verts[leg].second;
            const double x2 = verts[leg + 1].first;
            const double y2 = verts[leg + 1].second;
            f << si << ","
              << sc.srcId + 1 << "," << sc.destId + 1 << ","
              << leg << ",normal,"
              << x1 << "," << y1 << "," << x2 << "," << y2 << ","
              << x1 << "," << y1 << "," << x2 << "," << y2 << ","
              << sc.path.distance << "," << sc.path.bendCount << ","
              << (sc.everCrossed ? "true" : "false") << "\n";
        }
    }
}

void writeEffectiveIlCsv(
        const MILPSolveResult& layout,
        const DemandMatrix& D,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::vector<Node>& nodes,
        const std::string& suffix) {
    std::ofstream f("effective_demand_il_" + suffix + ".csv");
    f << "demand,sender,receiver,effective_distance,effective_bend_count,effective_il\n";
    for (int q = 0; q < (int)D.demands.size(); ++q) {
        const EffectiveDemandMetrics m = effectiveDemandMetrics(
            q, D, layout.demandDistance, layout.demandBendCount, layout.demandIL,
            shortcuts, nodes, layout.tour);
        f << q << ","
          << D.demands[q].first + 1 << ","
          << D.demands[q].second + 1 << ","
          << m.distance << ","
          << m.bendCount << ","
          << m.reportingIL << "\n";
    }
}

void writeRingIlCsv(
        const MILPSolveResult& layout,
        const DemandMatrix& D,
        const std::string& suffix) {
    std::ofstream f("demand_ring_il_" + suffix + ".csv");
    f << "demand,sender,receiver,ring_il\n";
    for (int q = 0; q < (int)D.demands.size(); ++q) {
        const double il = (q < (int)layout.demandIL.size())
            ? layout.demandIL[q] : 0.0;
        f << q << ","
          << D.demands[q].first + 1 << ","
          << D.demands[q].second + 1 << ","
          << il << "\n";
    }
}

void writeWavelengthsCsv(
        const WavelengthLoadBalanceResult* wlb,
        const ShortcutMethodResult& result,
        const DemandMatrix& D,
        const std::string& suffix) {
    std::ofstream f("wavelengths_" + suffix + ".csv");
    f << "kind,id,wavelength\n";
    if (wlb && wlb->success) {
        for (const auto& a : wlb->assignments) {
            f << "demand," << a.demandIdx << "," << a.wavelength << "\n";
            for (const auto& h : a.route.shortcutHops)
                f << "uses_shortcut," << a.demandIdx << "," << h.shortcutIdx << "\n";
        }
        return;
    }
    // Fallback: no WLB — leave empty (visualize.py uses greedy fallback).
    (void)result;
    (void)D;
}

}  // namespace

void exportProxyMasterPlotCsvs(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const ShortcutMethodResult& result,
        const WavelengthLoadBalanceResult* wlb,
        const std::string& suffix) {
    std::ofstream nodesFile("nodes.csv");
    nodesFile << "id,x,y\n";
    for (const auto& n : nodes)
        nodesFile << n.id + 1 << "," << n.x << "," << n.y << "\n";

    writeRingCsv(result.layout, nodes, "ring_" + suffix + ".csv");
    if (wlb && wlb->success)
        writeDemandsCsvFromWlb(*wlb, D, "demands_" + suffix + ".csv");
    else
        writeDemandsCsv(result.layout, D, "demands_" + suffix + ".csv");
    writeShortcutsCsv(result.shortcuts, "shortcuts_" + suffix + ".csv");
    writeEffectiveIlCsv(result.layout, D, result.shortcuts, nodes, suffix);
    writeRingIlCsv(result.layout, D, suffix);
    writeWavelengthsCsv(wlb, result, D, suffix);
}
