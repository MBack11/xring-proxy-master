#pragma once

#include <optional>
#include <set>
#include <vector>

#include "DemandMatrix.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "PhysicalConstants.h"
#include "RingTypes.h"
#include "ShortcutRouter.h"
#include "ShortcutTypes.h"

struct PlacedShortcut {
    int demandIdx = -1;
    int srcId = 0;
    int destId = 0;
    double totalIL = 0.0;
    bool everCrossed = false;
    int crossedShortcutIdx = -1;
    ShortcutPath path;
    std::vector<ShortcutPath> alternatives;
};

std::optional<PlacedShortcut> tryPlaceShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<Node>& nodes,
    const std::vector<int>& currentTour,
    const std::vector<EdgeOption>& currentRouting,
    const std::vector<Shortcut>& existingShortcutsForRouter,
    const std::set<int>& usedNodes,
    double sMin,
    double maxDistance,
    double currentRingDistance);

/// Route a shortcut between explicit endpoints while keeping `demandIdx` for
/// bookkeeping (used by WC endpoint relaxation). Does not require the raw
/// route distance to beat the demand's ring distance — caller decides
/// admissibility from the trial ring+shortcut graph.
std::optional<PlacedShortcut> tryPlaceShortcut(
    int demandIdx,
    int srcId,
    int destId,
    const std::vector<Node>& nodes,
    const std::vector<int>& currentTour,
    const std::vector<EdgeOption>& currentRouting,
    const std::vector<Shortcut>& existingShortcutsForRouter,
    const std::set<int>& usedNodes,
    double sMin,
    double maxDistance);

struct ShortcutRevalidationResult {
    int kept = 0;
    int repaired = 0;
    int dropped = 0;
};

/// After a ring re-solve: keep/repair shortcuts that stay strictly shorter than
/// both the demand's new ring distance and `wcCap` (typically the new WC).
/// Irreparable shortcuts are dropped (soft-fail) instead of rejecting the ring.
ShortcutRevalidationResult revalidateAllShortcuts(
    std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& newTour,
    const std::vector<EdgeOption>& newRouting,
    double sMin,
    const std::vector<double>& ringDemandDistance,
    double wcCap);

// Ring arc distances for every demand (ignores shortcut exclusions).
std::vector<double> allRingDemandDistances(
    const std::vector<Node>& nodes,
    const DemandMatrix& D,
    const MILPSolveResult& layout);

double worstRingRoutedDistance(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts);

bool shortcutPathCrossesRing(
    const ShortcutPath& path,
    const std::vector<Segment>& ringSegments,
    int srcId,
    int destId,
    const std::vector<int>& tour,
    const std::vector<EdgeOption>& routing);

std::vector<Shortcut> placedShortcutsToRouterFormat(
    const std::vector<PlacedShortcut>& shortcuts);

/// Exclusive (private): demand q may use a shortcut only if endpoints match
/// {s_q, t_q} (start→end owner). Applies to MIP flow coupling when wired
/// through, and to post-solve metrics / WC routing.
/// Shared: WC-aware multi-demand use + RWA λ assignment (current default).
enum class ShortcutUsageMode {
    Exclusive = 0,
    Shared = 1
};

/// Undirected endpoint match: (a,b) equals demand (s,t) in either orientation.
inline bool undirectedEndpointsMatch(int a, int b, int s, int t) {
    return (a == s && b == t) || (a == t && b == s);
}

ShortcutUsageMode getShortcutUsageMode();
void setShortcutUsageMode(ShortcutUsageMode mode);

struct ScopedShortcutUsageMode {
    explicit ScopedShortcutUsageMode(ShortcutUsageMode mode)
        : prev_(getShortcutUsageMode()) {
        setShortcutUsageMode(mode);
    }
    ~ScopedShortcutUsageMode() { setShortcutUsageMode(prev_); }
    ScopedShortcutUsageMode(const ScopedShortcutUsageMode&) = delete;
    ScopedShortcutUsageMode& operator=(const ScopedShortcutUsageMode&) = delete;
private:
    ShortcutUsageMode prev_;
};

/// Effective demand metrics under WC-aware routing:
/// W* = max_d sp(d) on the bidirectional ring+shortcut graph;
/// if ring(d) <= W*+eps stay on the ring, else use sp(d).
/// Endpoint-match fallback when tour is empty.
/// When ShortcutUsageMode::Exclusive: owner-only endpoint match.
struct EffectiveDemandMetrics {
    double distance = 0.0;
    int bendCount = 0;
    double reportingIL = 0.0;
};

EffectiveDemandMetrics effectiveDemandMetrics(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<double>& ringDemandIL,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// Endpoint-match fallback when tour is empty.
EffectiveDemandMetrics effectiveDemandMetrics(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<double>& ringDemandIL,
    const std::vector<PlacedShortcut>& shortcuts);

int findWorstCaseDemand(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

int findWorstCaseDemand(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts);

int findWorstCaseDemandEffective(
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour,
    const std::set<int>& excludedDemands = {});

int findWorstCaseDemandEffective(
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::set<int>& excludedDemands = {});

double effectiveDemandDistance(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

double effectiveDemandDistance(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts);

double computeGlobalW(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

double computeGlobalW(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts);

bool strictlyImprovesGlobalW(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    double oldGlobalW,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

bool strictlyImprovesGlobalW(
    int numDemands,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<int>& ringDemandBendCount,
    const std::vector<PlacedShortcut>& shortcuts,
    double oldGlobalW);

/// Exclusive endpoint match (placement / export / MILP exclusion).
bool demandHasShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts);

/// True if the WC-aware route uses any shortcut edge.
/// Without ring distances, ring(d) is treated as +∞ (free-SP usage check).
bool demandPathUsesShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// WC-aware: shortcut used only when ring(d) > W*+eps and SP uses a shortcut.
bool demandPathUsesShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// Case-3 overload: path uses a shortcut when nodes+tour provided.
bool demandHasShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// Case-3 overload with ring distances for WC-aware usage.
bool demandHasShortcut(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// Exclusive endpoint match.
const PlacedShortcut* findShortcutForDemand(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts);

/// First shortcut on the WC-aware route (else endpoint match).
const PlacedShortcut* findShortcutForDemand(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

/// WC-aware shortcut lookup with ring distances.
const PlacedShortcut* findShortcutForDemand(
    int demandIdx,
    const DemandMatrix& D,
    const std::vector<double>& ringDemandDistance,
    const std::vector<PlacedShortcut>& shortcuts,
    const std::vector<Node>& nodes,
    const std::vector<int>& tour);

void applyShortcutCrossingSideEffects(
    std::vector<PlacedShortcut>& shortcuts,
    const PlacedShortcut& placed);
