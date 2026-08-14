#pragma once

#include <set>
#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ShortcutOrchestrator.h"

struct WcIterationRecord {
    int iteration = 0;
    int wcDemand = -1;
    int sender = 0;
    int receiver = 0;
    double ringIL = 0.0;
    bool shortcutPlaced = false;
    double shortcutIL = 0.0;
};

struct WavelengthUsage {
    int totalWavelengths = 0;
    /// 1-indexed λ per demand (end-to-end; every demand gets a channel).
    std::vector<int> demandWavelength;
    /// Optional: 1-indexed λ tag per shortcut (max/first user); may be empty.
    std::vector<int> shortcutWavelength;
    /// Shortcuts used on each demand's WC-aware route.
    std::vector<std::set<int>> demandShortcutIndices;
};

WavelengthUsage computeWavelengthUsage(
    const MILPSolveResult& layout,
    const DemandMatrix& D,
    const std::vector<Node>& nodes,
    const std::vector<PlacedShortcut>& shortcuts = {});

/// Shared-shortcut / end-to-end wavelength assignment: ONE channel per demand.
/// Same directed ring edge → conflict. Opposite directions on the same physical
/// edge do NOT conflict. Same or intersecting shortcuts → conflict.
/// Greedy first-fit = next free wavelength.
WavelengthUsage computeWavelengthUsageShared(
    const MILPSolveResult& layout,
    const DemandMatrix& D,
    const std::vector<Node>& nodes,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<std::set<int>>& shortcutIndicesPerDemand);

void logWavelengthUsage(
    const char* label,
    const WavelengthUsage& usage);

/// Export demand/shortcut → λ map for visualize.py (wavelengths_{suffix}.csv).
void exportWavelengthAssignment(
    const WavelengthUsage& usage,
    const std::string& suffix);

void exportRingDemandIL(
    const MILPSolveResult& layout,
    const DemandMatrix& D,
    const std::string& suffix);

void exportRingSnapshot(
    const MILPSolveResult& layout,
    const std::vector<Node>& nodes,
    const std::string& filename);

void exportDemandsSnapshot(
    const MILPSolveResult& layout,
    const DemandMatrix& D,
    const std::vector<Node>& nodes,
    const std::string& filename);

void exportShortcutHistory(
    const std::vector<WcIterationRecord>& history,
    const std::string& suffix);

void exportFinalCSVs(
    const MILPSolveResult& layout,
    const std::vector<PlacedShortcut>& shortcuts,
    const DemandMatrix& D,
    const std::set<int>& excludedDemands,
    const std::string& suffix,
    const std::vector<Node>& nodes);

void exportEffectiveDemandIL(
    const MILPSolveResult& layout,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::string& suffix,
    const std::vector<Node>& nodes);

// Diagnostic: pairwise ring-segment crossing check from exported ring_*.csv.
bool checkRingCsvCrossings(const std::string& csvPath);
