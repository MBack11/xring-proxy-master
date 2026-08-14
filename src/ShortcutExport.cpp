#include "ShortcutExport.h"

void exportRingDemandIL(
        const MILPSolveResult&,
        const DemandMatrix&,
        const std::string&) {}

void exportRingSnapshot(
        const MILPSolveResult&,
        const std::vector<Node>&,
        const std::string&) {}

void exportDemandsSnapshot(
        const MILPSolveResult&,
        const DemandMatrix&,
        const std::vector<Node>&,
        const std::string&) {}

void exportFinalCSVs(
        const MILPSolveResult&,
        const std::vector<PlacedShortcut>&,
        const DemandMatrix&,
        const std::set<int>&,
        const std::string&,
        const std::vector<Node>&) {}

void exportShortcutHistory(
        const std::vector<WcIterationRecord>&,
        const std::string&) {}

void exportEffectiveDemandIL(
        const MILPSolveResult&,
        const DemandMatrix&,
        const std::vector<PlacedShortcut>&,
        const std::string&,
        const std::vector<Node>&) {}

void exportWavelengthAssignment(
        const WavelengthUsage&,
        const std::string&) {}

WavelengthUsage computeWavelengthUsage(
        const MILPSolveResult&,
        const DemandMatrix&,
        const std::vector<Node>&,
        const std::vector<PlacedShortcut>&) {
    return {};
}

WavelengthUsage computeWavelengthUsageShared(
        const MILPSolveResult&,
        const DemandMatrix&,
        const std::vector<Node>&,
        const std::vector<PlacedShortcut>&,
        const std::vector<std::set<int>>&) {
    return {};
}

void logWavelengthUsage(const char*, const WavelengthUsage&) {}

bool checkRingCsvCrossings(const std::string&) { return false; }
