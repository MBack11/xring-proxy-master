#include "ShortcutMethods.h"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gurobi_c++.h"
#include "RingLayout.h"
#include "ShortcutExport.h"
#include "ShortcutGrid.h"

namespace {

using Clock = std::chrono::steady_clock;

double elapsedSec(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

const char* gurobiStatusName(int status) {
    switch (status) {
        case GRB_OPTIMAL: return "OPTIMAL";
        case GRB_TIME_LIMIT: return "TIME_LIMIT";
        case GRB_USER_OBJ_LIMIT: return "USER_OBJ_LIMIT";
        case GRB_INFEASIBLE: return "INFEASIBLE";
        default: return "OTHER";
    }
}

std::set<int> shortcutDemandIndices(
        const std::vector<PlacedShortcut>& shortcuts,
        const DemandMatrix& D) {
    std::set<int> excluded;
    for (int q = 0; q < (int)D.demands.size(); ++q) {
        if (demandHasShortcut(q, D, shortcuts))
            excluded.insert(q);
    }
    return excluded;
}

RingLayout layoutToRing(const std::vector<Node>& nodes, const MILPSolveResult& layout) {
    return buildRingLayout(nodes, layout);
}

void recordWImprovement(ShortcutMethodProfile* profile, double before, double after) {
    if (!profile) return;
    if (after + MILP_EPS < before)
        ++profile->wImprovementEvents;
}

/// Tour neighbors of `nodeId`: itself, left, right (deduped). A node is
/// admissible only if it is not already a shortcut endpoint in `usedNodes`.
std::vector<int> admissibleEndpointChoices(
        int nodeId,
        const std::vector<int>& tour,
        const std::set<int>& usedNodes) {
    std::vector<int> out;
    const int n = (int)tour.size();
    int idx = -1;
    for (int i = 0; i < n; ++i) {
        if (tour[i] == nodeId) {
            idx = i;
            break;
        }
    }
    if (idx < 0)
        return out;

    const int candidates[3] = {
        nodeId,
        tour[(idx - 1 + n) % n],
        tour[(idx + 1) % n],
    };
    std::set<int> seen;
    for (int id : candidates) {
        if (usedNodes.count(id))
            continue;
        if (seen.insert(id).second)
            out.push_back(id);
    }
    return out;
}

struct RelaxedWcPlaceResult {
    PlacedShortcut placed;
    double globalWAfter = 0.0;
};

/// Endpoint relaxation for one WC demand: try up to 3×3 neighbor pairs,
/// accept only if q* ceases to be WC and global W improves strictly.
std::optional<RelaxedWcPlaceResult> tryPlaceRelaxedWcShortcut(
        int wcDemand,
        const DemandMatrix& D,
        const std::vector<Node>& nodes,
        const std::vector<int>& tour,
        const std::vector<EdgeOption>& routing,
        const std::vector<PlacedShortcut>& shortcuts,
        const std::set<int>& usedNodes,
        double sMin,
        double maxDistance,
        double ringDistance,
        double globalWBefore,
        const std::vector<double>& ringDemandDistance,
        const std::vector<int>& ringDemandBendCount,
        const ShortcutMethodOptions& options) {
    if (wcDemand < 0 || wcDemand >= (int)D.demands.size())
        return std::nullopt;

    const int sStar = D.demands[wcDemand].first;
    const int dStar = D.demands[wcDemand].second;
    const auto srcChoices = admissibleEndpointChoices(sStar, tour, usedNodes);
    const auto destChoices = admissibleEndpointChoices(dStar, tour, usedNodes);
    const auto existingRouter = placedShortcutsToRouterFormat(shortcuts);
    const int numDemands = (int)D.demands.size();

    std::optional<RelaxedWcPlaceResult> best;
    int bestBends = std::numeric_limits<int>::max();
    bool bestIsBaseline = false;

    for (int sPrime : srcChoices) {
        for (int dPrime : destChoices) {
            if (sPrime == dPrime)
                continue;

            auto placed = tryPlaceShortcut(
                wcDemand, sPrime, dPrime, nodes, tour, routing,
                existingRouter, usedNodes, sMin, maxDistance);
            if (!placed)
                continue;

            std::vector<PlacedShortcut> trial = shortcuts;
            trial.push_back(*placed);
            applyShortcutCrossingSideEffects(trial, *placed);

            const double W_cand = computeGlobalW(
                numDemands, D, ringDemandDistance, ringDemandBendCount,
                trial, nodes, tour);
            const int argmaxCand = findWorstCaseDemandEffective(
                D, ringDemandDistance, ringDemandBendCount, trial,
                nodes, tour, {});

            // Admissible iff q* is no longer WC and global W improves strictly.
            if (argmaxCand == wcDemand)
                continue;
            if (!(W_cand + MILP_EPS < globalWBefore))
                continue;

            const EffectiveDemandMetrics qMetrics = effectiveDemandMetrics(
                wcDemand, D, ringDemandDistance, ringDemandBendCount, {},
                trial, nodes, tour);
            const int bendsQ = qMetrics.bendCount;
            const bool isBaseline = (sPrime == sStar && dPrime == dStar);

            bool take = false;
            if (!best) {
                take = true;
            } else if (W_cand + MILP_EPS < best->globalWAfter) {
                take = true;
            } else if (std::abs(W_cand - best->globalWAfter) <= MILP_EPS) {
                if (bendsQ < bestBends) {
                    take = true;
                } else if (bendsQ == bestBends && isBaseline && !bestIsBaseline) {
                    take = true;
                }
            }

            if (take) {
                best = RelaxedWcPlaceResult{*placed, W_cand};
                bestBends = bendsQ;
                bestIsBaseline = isBaseline;
            }
        }
    }

    if (!best)
        return std::nullopt;

    recordWImprovement(options.profile, globalWBefore, best->globalWAfter);
    if (!options.quiet) {
        const bool relaxed =
            best->placed.srcId != sStar || best->placed.destId != dStar;
        std::cout << std::fixed << std::setprecision(1)
                  << "[Step C] Placed shortcut for demand " << wcDemand
                  << " (" << best->placed.path.distance << " mm vs ring "
                  << ringDistance << " mm), global W "
                  << globalWBefore << " -> " << best->globalWAfter;
        if (relaxed) {
            std::cout << " [endpoint relaxation "
                      << sStar + 1 << "->" << dStar + 1 << " → "
                      << best->placed.srcId + 1 << "->"
                      << best->placed.destId + 1 << "]";
        }
        std::cout << "\n";
    }
    return best;
}

void logResolveAttempt(
        ShortcutMethodProfile* profile,
        int iteration,
        const std::string& outcome,
        double wBefore,
        double wAfter) {
    if (!profile) return;
    profile->resolveAttempts.push_back({iteration, outcome, wBefore, wAfter});
}

double remainingWallSec(
        const ShortcutMethodOptions& options,
        Clock::time_point loopStart) {
    if (options.wallTimeLimitSec < 0.0)
        return options.resolveTimeLimitSec;
    return std::max(0.0, options.wallTimeLimitSec - elapsedSec(loopStart));
}

bool wallBudgetExhausted(
        const ShortcutMethodOptions& options,
        Clock::time_point loopStart,
        double minSec = 0.5) {
    if (options.wallTimeLimitSec < 0.0)
        return false;
    return remainingWallSec(options, loopStart) < minSec;
}

bool stopForWallTime(
        ShortcutMethodResult& out,
        const ShortcutMethodOptions& options,
        Clock::time_point loopStart,
        const char* logLabel) {
    if (!wallBudgetExhausted(options, loopStart))
        return false;
    out.termination = ShortcutTermination::WallTimeLimit;
    if (!options.quiet) {
        std::cout << "[" << logLabel << "] Terminated: wall-time limit ("
                  << options.wallTimeLimitSec << " s) reached\n";
    }
    return true;
}

}  // namespace

const char* terminationToString(ShortcutTermination reason) {
    switch (reason) {
        case ShortcutTermination::Case2Constraints:
            return "Case 2: shortcut constraints violated";
        case ShortcutTermination::Case3NoFurtherImprovement:
            return "Case 3: WC demand already has shortest shortcut route";
        case ShortcutTermination::WallTimeLimit:
            return "Wall-time limit reached";
        default:
            return "unknown";
    }
}

ShortcutMethodResult runMethodAWithShortcuts(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layoutFixed,
        double sMin,
        const ShortcutMethodOptions& options) {
    ScopedShortcutUsageMode usageScope(options.usageMode);
    ShortcutMethodResult out;
    out.layout = layoutFixed;

    if (!options.skipExports) {
        exportRingDemandIL(layoutFixed, D, "A");
        exportRingSnapshot(layoutFixed, nodes, "ring_before_A.csv");
        exportDemandsSnapshot(layoutFixed, D, nodes, "demands_ring_before_A.csv");
    }

    RingLayout ring = layoutToRing(nodes, layoutFixed);
    std::vector<PlacedShortcut> shortcuts;
    std::vector<WcIterationRecord> wcHistory;
    std::set<int> usedNodes;
    const int numDemands = (int)D.demands.size();
    const int MAX_ITERATIONS = 50;
    int iteration = 0;

    double W_current = computeGlobalW(
        numDemands, D, layoutFixed.demandDistance, layoutFixed.demandBendCount,
        shortcuts, nodes, ring.tour);

    if (!options.quiet) {
        std::cout << "\n======================================================\n";
        std::cout << "INNER LOOP: Method A+Shortcuts\n";
        std::cout << "======================================================\n";
    }

    const auto loopStart = Clock::now();

    while (true) {
        iteration++;
        if (iteration > MAX_ITERATIONS) {
            if (!options.quiet)
                std::cerr << "[InnerLoop] MAX_ITERATIONS safety cap reached.\n";
            break;
        }
        if (stopForWallTime(out, options, loopStart, "Method A+Shortcuts"))
            break;

        const int wcDemand = findWorstCaseDemandEffective(
            D, layoutFixed.demandDistance, layoutFixed.demandBendCount, shortcuts,
            nodes, ring.tour, {});
        if (wcDemand < 0) break;

        const auto& wcPair = D.demands[wcDemand];
        const double wcEffective = effectiveDemandDistance(
            wcDemand, D, layoutFixed.demandDistance, shortcuts, nodes, ring.tour);
        if (!options.quiet) {
            std::cout << std::fixed << std::setprecision(1)
                      << "[InnerLoop] WC demand " << wcDemand << " ("
                      << wcPair.first + 1 << "->" << wcPair.second + 1
                      << "), effective distance=" << wcEffective << " mm\n";
        }

        if (demandHasShortcut(
                wcDemand, D, layoutFixed.demandDistance, shortcuts, nodes, ring.tour)) {
            out.termination = ShortcutTermination::Case3NoFurtherImprovement;
            if (!options.quiet) {
                std::cout << "[InnerLoop] Terminated: no improvement for WC demand "
                          << wcDemand << "\n";
            }
            break;
        }

        const double ringDistance = layoutFixed.demandDistance[wcDemand];
        const double ringReportingIL = layoutFixed.demandIL[wcDemand];
        const double globalWBefore = W_current;
        const double maxDistance = ringDistance;

        WcIterationRecord rec;
        rec.iteration = iteration;
        rec.wcDemand = wcDemand;
        rec.sender = wcPair.first + 1;
        rec.receiver = wcPair.second + 1;
        rec.ringIL = ringReportingIL;

        const auto searchStart = Clock::now();
        auto candidate = tryPlaceRelaxedWcShortcut(
            wcDemand, D, nodes, ring.tour, ring.routing, shortcuts, usedNodes,
            sMin, maxDistance, ringDistance, globalWBefore,
            layoutFixed.demandDistance, layoutFixed.demandBendCount, options);
        if (options.profile)
            options.profile->shortcutSearchSec += elapsedSec(searchStart);
        if (stopForWallTime(out, options, loopStart, "Method A+Shortcuts"))
            break;

        if (!candidate) {
            rec.shortcutPlaced = false;
            rec.shortcutIL = 0.0;
            wcHistory.push_back(rec);
            out.termination = ShortcutTermination::Case2Constraints;
            if (!options.quiet) {
                std::cout << "[InnerLoop] Terminated: no improvement for WC demand "
                          << wcDemand << "\n";
            }
            break;
        }

        std::vector<PlacedShortcut> trialShortcuts = shortcuts;
        trialShortcuts.push_back(candidate->placed);
        applyShortcutCrossingSideEffects(trialShortcuts, candidate->placed);

        rec.shortcutPlaced = true;
        rec.shortcutIL = candidate->placed.totalIL;
        wcHistory.push_back(rec);

        shortcuts = std::move(trialShortcuts);
        usedNodes.insert(candidate->placed.srcId);
        usedNodes.insert(candidate->placed.destId);
        W_current = candidate->globalWAfter;
    }

    out.shortcuts = shortcuts;
    out.globalW = W_current;

    if (!options.skipExports) {
        std::set<int> excluded = shortcutDemandIndices(shortcuts, D);
        exportShortcutHistory(wcHistory, "A");
        exportFinalCSVs(out.layout, shortcuts, D, excluded, "A", nodes);
    }

    if (!options.quiet) {
        std::cout << std::fixed << std::setprecision(1)
                  << "\n[Method A+Shortcuts] Final global W: " << out.globalW << " mm\n";
        std::cout << "[Method A+Shortcuts] Shortcuts placed: " << shortcuts.size() << "\n";
    }

    return out;
}

ShortcutMethodResult runMethodBWithShortcuts(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        MILPSolver& solver,
        const MILPSolveResult& initialLayout,
        double sMin,
        const std::string& exportSuffix,
        const char* logLabel,
        const ShortcutMethodOptions& options) {
    ScopedShortcutUsageMode usageScope(options.usageMode);
    ShortcutMethodResult out;

    MILPSolveResult layoutCurrent = initialLayout;
    if (!layoutCurrent.success) {
        if (!options.quiet)
            std::cerr << "[" << logLabel << "] Initial MILP layout invalid.\n";
        return out;
    }

    if (!options.skipExports) {
        exportRingDemandIL(layoutCurrent, D, exportSuffix);
        exportRingSnapshot(layoutCurrent, nodes, "ring_before_" + exportSuffix + ".csv");
        exportDemandsSnapshot(layoutCurrent, D, nodes,
                              "demands_ring_before_" + exportSuffix + ".csv");
    }

    double W_current = computeGlobalW(
        (int)D.demands.size(), D, layoutCurrent.demandDistance, layoutCurrent.demandBendCount,
        {}, nodes, layoutCurrent.tour);
    std::vector<PlacedShortcut> shortcuts;
    std::vector<WcIterationRecord> wcHistory;
    std::set<int> usedNodes;
    const int numDemands = (int)D.demands.size();
    const int MAX_ITERATIONS = 50;
    int iteration = 0;

    if (!options.quiet) {
        std::cout << "\n======================================================\n";
        std::cout << "INNER LOOP: " << logLabel << "\n";
        std::cout << "======================================================\n";
    }

    const auto loopStart = Clock::now();

    while (true) {
        iteration++;
        if (iteration > MAX_ITERATIONS) {
            if (!options.quiet)
                std::cerr << "[InnerLoop] MAX_ITERATIONS safety cap reached.\n";
            break;
        }
        if (stopForWallTime(out, options, loopStart, logLabel))
            break;

        const int wcDemand = findWorstCaseDemandEffective(
            D, layoutCurrent.demandDistance, layoutCurrent.demandBendCount, shortcuts,
            nodes, layoutCurrent.tour, {});

        if (wcDemand < 0) break;

        RingLayout ring = layoutToRing(nodes, layoutCurrent);

        if (demandHasShortcut(
                wcDemand, D, layoutCurrent.demandDistance, shortcuts, nodes,
                ring.tour)) {
            out.termination = ShortcutTermination::Case3NoFurtherImprovement;
            if (!options.quiet) {
                std::cout << "[InnerLoop] Terminated: no improvement for WC demand "
                          << wcDemand << "\n";
            }
            break;
        }

        const auto& wcPair = D.demands[wcDemand];
        const double wcEffective = effectiveDemandDistance(
            wcDemand, D, layoutCurrent.demandDistance, shortcuts, nodes, ring.tour);
        if (!options.quiet) {
            std::cout << std::fixed << std::setprecision(1)
                      << "[InnerLoop] WC demand " << wcDemand << " ("
                      << wcPair.first + 1 << "->" << wcPair.second + 1
                      << "), effective distance=" << wcEffective << " mm\n";
        }

        const double ringDistance = layoutCurrent.demandDistance[wcDemand];
        const double ringReportingIL = layoutCurrent.demandIL[wcDemand];
        const double globalWBefore = W_current;
        const double maxDistance = ringDistance;

        WcIterationRecord rec;
        rec.iteration = iteration;
        rec.wcDemand = wcDemand;
        rec.sender = wcPair.first + 1;
        rec.receiver = wcPair.second + 1;
        rec.ringIL = ringReportingIL;

        const auto searchStart = Clock::now();
        auto shortcutResult = tryPlaceRelaxedWcShortcut(
            wcDemand, D, nodes, ring.tour, ring.routing, shortcuts, usedNodes,
            sMin, maxDistance, ringDistance, globalWBefore,
            layoutCurrent.demandDistance, layoutCurrent.demandBendCount, options);
        if (options.profile)
            options.profile->shortcutSearchSec += elapsedSec(searchStart);
        if (stopForWallTime(out, options, loopStart, logLabel))
            break;

        if (!shortcutResult) {
            rec.shortcutPlaced = false;
            rec.shortcutIL = 0.0;
            wcHistory.push_back(rec);
            out.termination = ShortcutTermination::Case2Constraints;
            if (!options.quiet) {
                std::cout << "[InnerLoop] Terminated: no improvement for WC demand "
                          << wcDemand << "\n";
            }
            break;
        }

        std::vector<PlacedShortcut> trialShortcuts = shortcuts;
        trialShortcuts.push_back(shortcutResult->placed);
        applyShortcutCrossingSideEffects(trialShortcuts, shortcutResult->placed);

        rec.shortcutPlaced = true;
        rec.shortcutIL = shortcutResult->placed.totalIL;
        wcHistory.push_back(rec);

        shortcuts = std::move(trialShortcuts);
        usedNodes.insert(shortcutResult->placed.srcId);
        usedNodes.insert(shortcutResult->placed.destId);
        W_current = shortcutResult->globalWAfter;

        int wcNew = findWorstCaseDemandEffective(
            D, layoutCurrent.demandDistance, layoutCurrent.demandBendCount, shortcuts,
            nodes, ring.tour, {});
        if (!demandHasShortcut(
                wcNew, D, layoutCurrent.demandDistance, shortcuts, nodes,
                ring.tour)) {
            std::set<int> excludedDemands = shortcutDemandIndices(shortcuts, D);

            const double resolveLimit = std::min(
                options.resolveTimeLimitSec, remainingWallSec(options, loopStart));
            if (resolveLimit < 0.5) {
                logResolveAttempt(options.profile, iteration, "WALL_TIME_LIMIT",
                                  W_current, W_current);
                out.termination = ShortcutTermination::WallTimeLimit;
                if (!options.quiet) {
                    std::cout << "[" << logLabel << "] Terminated: wall-time limit ("
                              << options.wallTimeLimitSec << " s) reached before re-solve\n";
                }
                break;
            }

            const auto resolveStart = Clock::now();
            MILPSolveOptions reSolveOpts;
            MILPSolveResult candidate = solver.solve(
                true, "", excludedDemands, W_current,
                &layoutCurrent, resolveLimit, true, reSolveOpts);
            const double resolveElapsed = elapsedSec(resolveStart);
            if (options.profile)
                options.profile->resolveSec += resolveElapsed;
            if (stopForWallTime(out, options, loopStart, logLabel))
                break;

            const double wBeforeResolve = W_current;
            if (candidate.status == GRB_TIME_LIMIT && !candidate.success) {
                logResolveAttempt(options.profile, iteration, "TIMEOUT",
                                  wBeforeResolve, wBeforeResolve);
            } else if (!candidate.success) {
                logResolveAttempt(options.profile, iteration, "SOLVE_FAILED",
                                  wBeforeResolve, wBeforeResolve);
            } else {
                // Revalidate on a copy: keep/repair shortcuts shorter than the new WC;
                // drop irreparable ones instead of rejecting the whole ring.
                RingLayout candidateRing = layoutToRing(nodes, candidate);
                const std::vector<double> ringDistForReval =
                    allRingDemandDistances(nodes, D, candidate);
                if (ringDistForReval.empty()) {
                    logResolveAttempt(options.profile, iteration, "REJECTED_REVALIDATION",
                                      wBeforeResolve, wBeforeResolve);
                } else {
                    std::vector<PlacedShortcut> trialShortcuts = shortcuts;
                    const auto ringSegs =
                        ShortcutGrid::collectRingSegments(candidateRing.routing);

                    // Provisional WC under candidate ring: bidirectional effective
                    // distance using non-crossing shortcuts as graph edges.
                    std::vector<PlacedShortcut> usable;
                    usable.reserve(trialShortcuts.size());
                    for (const PlacedShortcut& sc : trialShortcuts) {
                        if (!shortcutPathCrossesRing(
                                sc.path, ringSegs, sc.srcId, sc.destId,
                                candidateRing.tour, candidateRing.routing))
                            usable.push_back(sc);
                    }
                    double wcCap = 0.0;
                    for (int q = 0; q < numDemands; ++q) {
                        const double d = effectiveDemandDistance(
                            q, D, ringDistForReval, usable,
                            nodes, candidateRing.tour);
                        wcCap = std::max(wcCap, d);
                    }

                    const ShortcutRevalidationResult rev = revalidateAllShortcuts(
                        trialShortcuts, nodes, candidateRing.tour,
                        candidateRing.routing, sMin, ringDistForReval, wcCap);

                    // Candidate layout still excludes shortcut demands in demandDistance;
                    // recompute global W with repaired/dropped shortcuts + full ring dists.
                    MILPSolveResult routed = candidate;
                    if ((int)routed.demandDistance.size() != numDemands
                        || (int)routed.demandBendCount.size() != numDemands) {
                        routed.demandDistance = ringDistForReval;
                        routed.demandBendCount.assign(numDemands, 0);
                    } else {
                        for (int q = 0; q < numDemands; ++q) {
                            if (q < (int)ringDistForReval.size())
                                routed.demandDistance[q] = ringDistForReval[q];
                        }
                    }

                    const double newGlobalW = computeGlobalW(
                        numDemands, D, routed.demandDistance, routed.demandBendCount,
                        trialShortcuts, nodes, candidateRing.tour);

                    if (newGlobalW < W_current - 1e-4) {
                        if (!options.quiet) {
                            std::cout << std::fixed << std::setprecision(1)
                                      << "[" << logLabel << "] Iteration " << iteration
                                      << ": re-solve status=" << gurobiStatusName(candidate.status)
                                      << ", W: " << W_current << " -> " << newGlobalW
                                      << " mm, WC_cap=" << wcCap
                                      << " mm, shortcuts kept/repaired/dropped="
                                      << rev.kept << "/" << rev.repaired << "/"
                                      << rev.dropped
                                      << ", excluded=" << excludedDemands.size()
                                      << ", targeting WC demand " << wcNew << "\n";
                        }
                        logResolveAttempt(options.profile, iteration, "ACCEPTED",
                                          wBeforeResolve, newGlobalW);
                        recordWImprovement(options.profile, wBeforeResolve, newGlobalW);
                        layoutCurrent = candidate;
                        // Keep full ring distances on layout for subsequent WC search.
                        layoutCurrent.demandDistance = routed.demandDistance;
                        layoutCurrent.demandBendCount = routed.demandBendCount;
                        shortcuts = std::move(trialShortcuts);
                        usedNodes.clear();
                        for (const PlacedShortcut& sc : shortcuts) {
                            usedNodes.insert(sc.srcId);
                            usedNodes.insert(sc.destId);
                        }
                        W_current = newGlobalW;
                    } else {
                        logResolveAttempt(options.profile, iteration, "REJECTED_NO_IMPROVEMENT",
                                          wBeforeResolve, newGlobalW);
                        if (!options.quiet) {
                            std::cout << "[" << logLabel << "] Iteration " << iteration
                                      << ": re-solve status=" << gurobiStatusName(candidate.status)
                                      << ", no improvement after revalidation"
                                      << " (dropped=" << rev.dropped << "), continuing\n";
                        }
                    }
                }
            }
        }
    }

    out.layout = layoutCurrent;
    out.shortcuts = shortcuts;
    out.globalW = W_current;

    if (!options.skipExports) {
        std::set<int> excluded = shortcutDemandIndices(shortcuts, D);
        exportShortcutHistory(wcHistory, exportSuffix);
        exportFinalCSVs(out.layout, shortcuts, D, excluded, exportSuffix, nodes);
    }

    if (!options.quiet) {
        std::cout << std::fixed << std::setprecision(1)
                  << "\n[" << logLabel << "] Final global W: " << out.globalW << " mm\n";
        std::cout << "[" << logLabel << "] Shortcuts placed: " << shortcuts.size() << "\n";
    }

    return out;
}
