#include "ProxyMasterLoop.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <utility>

#include "LShapeCrossing.h"
#include "MILPSolver.h"
#include "PhysicalConstants.h"
#include "StageE.h"

namespace {

using Clock = std::chrono::steady_clock;

std::string tourKey(const MILPSolveResult& layout) {
    std::vector<std::pair<int, int>> edges = ringUndirectedEdges(layout);
    std::string s;
    for (const auto& [a, b] : edges)
        s += std::to_string(a) + "-" + std::to_string(b) + ";";
    return s;
}

double runMethodDOnLayout(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const MILPSolveResult& layout,
        double sMin,
        bool quiet,
        ShortcutMethodResult* out,
        ShortcutUsageMode usageMode,
        double abortIfObjBoundGe = std::numeric_limits<double>::infinity(),
        double abortOnlyIfObjBoundGt = -std::numeric_limits<double>::infinity()) {
    ShortcutMethodOptions opt;
    opt.quiet = quiet;
    opt.skipExports = true;
    opt.usageMode = usageMode;
    opt.abortIfObjBoundGe = abortIfObjBoundGe;
    opt.abortOnlyIfObjBoundGt = abortOnlyIfObjBoundGt;
    ShortcutMethodResult res = runMethodDJointShortcuts(nodes, D, layout, sMin, opt);
    if (out) *out = res;
    return res.globalW;
}

}  // namespace

double runMethodDBaselineW(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        double sMin,
        bool quiet,
        ShortcutUsageMode usageMode) {
    MILPSolver solverA(nodes, D);
    MILPSolveResult layoutA = solverA.solve(false, "", {}, -1.0, nullptr, -1.0, quiet);
    if (!layoutA.success) return std::numeric_limits<double>::infinity();

    MILPSolver solverB(nodes, D);
    MILPSolveResult layoutB = solverB.solve(true, "", {}, -1.0, nullptr, -1.0, quiet);
    const MILPSolveResult& ring = layoutB.success ? layoutB : layoutA;
    return runMethodDOnLayout(nodes, D, ring, sMin, quiet, nullptr, usageMode);
}

ProxyMasterLoopResult runProxyMasterLoop(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        const ProxyMasterLoopOptions& options) {
    ScopedShortcutUsageMode usageScope(options.usageMode);
    ProxyMasterLoopResult result;
    result.poolModeNote =
        "Gurobi PoolSearchMode=2, PoolGap=0.0 (proven optimal pool)";

    const std::vector<ShortcutPairIndex> pairs = buildMasterShortcutPairs(nodes);
    std::vector<MasterNogoodCut> cuts;
    std::set<std::string> seenRings;

    // Incremental master: one model, append cuts. Non-incremental rebuilds each call.
    std::unique_ptr<ProxyMasterMilpSession> masterSession;
    if (options.incrementalMaster)
        masterSession = std::make_unique<ProxyMasterMilpSession>(
            nodes, D, options.usageMode);

    auto solveMaster = [&](const ProxyMasterMilpOptions& mOpt) {
        if (masterSession)
            return masterSession->solve(cuts, mOpt);
        return solveProxyMasterMilp(nodes, D, cuts, mOpt);
    };

    const auto t0 = Clock::now();

    for (int round = 0; round < options.maxRounds; ++round) {
        if (options.wallTimeLimitSec > 0.0) {
            const double elapsed = std::chrono::duration<double>(
                Clock::now() - t0).count();
            if (elapsed >= options.wallTimeLimitSec) {
                result.stoppedBudget = true;
                result.stopReason = "wall-time budget";
                break;
            }
        }

        ProxyMasterMilpOptions mOpt;
        mOpt.quiet = options.quiet;
        mOpt.timeLimitSec = options.masterTimeLimitSec;
        mOpt.poolSolutions = options.poolSolutions;
        mOpt.usageMode = options.usageMode;

        const ProxyMasterMilpResult master = solveMaster(mOpt);
        if (!master.success || master.pool.empty()) {
            result.stopReason = "master infeasible or empty pool";
            break;
        }

        result.LB = master.Wproxy;
        result.rounds = round + 1;
        const double LBk = master.Wproxy;

        ProxyMasterLoopRoundTrace trace;
        if (options.traceRounds) {
            trace.round = round + 1;
            trace.LBk = LBk;
            trace.incumbentBefore = result.Wstar;
            trace.poolSize = (int)master.pool.size();
            for (const MasterSolution& sol : master.pool) {
                const std::string key = tourKey(sol.layout);
                trace.poolKeys.push_back(key);
                if (seenRings.count(key))
                    ++trace.poolDuplicate;
                else
                    ++trace.poolNew;
                if (seenRings.count(key)) {
                    trace.resurfacedKeys.push_back(key);
                    ++trace.poolResurfacedCut;
                }
            }
        }

        struct Ranked {
            MasterSolution sol;
            double what = 0.0;
        };
        std::vector<Ranked> ranked;
        ranked.reserve(master.pool.size());
        for (const MasterSolution& sol : master.pool) {
            const std::string key = tourKey(sol.layout);
            if (seenRings.count(key)) continue;  // skip Stage E on already-evaluated rings
            Ranked r;
            r.sol = sol;
            r.what = computeStageE(
                nodes, D, sol.layout, sol.selectedShortcutIndices, pairs,
                options.sMin, options.usageMode);
            ranked.push_back(r);
        }
        std::sort(ranked.begin(), ranked.end(),
            [](const Ranked& a, const Ranked& b) { return a.what < b.what; });

        bool fineStop = false;
        for (const Ranked& r : ranked) {
            const std::string key = tourKey(r.sol.layout);
            if (seenRings.count(key)) continue;
            seenRings.insert(key);

            ShortcutMethodResult dRes;
            // Dual-safe early abort: cannot improve W* and cannot fine-stop
            const double abortGe = result.Wstar;
            const double abortOnlyGt = LBk;
            const double wTrue = runMethodDOnLayout(
                nodes, D, r.sol.layout, options.sMin, options.quiet, &dRes,
                options.usageMode, abortGe, abortOnlyGt);
            ++result.ringsEvaluated;

            if (options.traceRounds) {
                if (!trace.stageVInvoked) {
                    trace.stageVInvoked = true;
                    trace.stageVRingKey = key;
                    trace.stageVWtrue = wTrue;
                }
                ProxyMasterLoopStageVTrace sv;
                sv.ringKey = key;
                sv.wTrue = wTrue;
                sv.fineTestFired = (std::isfinite(wTrue) && wTrue <= LBk + MILP_EPS);
                trace.stageVEvals.push_back(sv);
                ++trace.cutsAdded;
            }

            if (std::isfinite(wTrue) && wTrue + MILP_EPS < result.Wstar) {
                result.Wstar = wTrue;
                result.best = dRes;
                result.bestFoundRound = round + 1;
            }

            MasterNogoodCut cut;
            for (const auto& te : r.sol.layout.tourEdges)
                cut.directedArcs.push_back({te.from, te.to});
            cuts.push_back(cut);

            if (std::isfinite(wTrue) && wTrue <= LBk + MILP_EPS) {
                fineStop = true;
                result.stoppedFine = true;
                const double fineGap = result.Wstar - LBk;
                if (fineGap >= -MILP_EPS)
                    result.provenOptimal = true;
                result.stopReason = "fine test: W_true(R) <= LB_k + MILP_EPS";
                if (options.traceRounds) {
                    trace.fineStop = true;
                    trace.incumbentAfter = result.Wstar;
                    trace.incumbentChanged =
                        (trace.incumbentAfter + MILP_EPS < trace.incumbentBefore);
                    result.roundTraces.push_back(trace);
                }
                break;
            }
        }

        if (fineStop) break;

        ProxyMasterMilpOptions mOpt2 = mOpt;
        mOpt2.poolSolutions = 1;
        const ProxyMasterMilpResult rescan = solveMaster(mOpt2);
        if (!rescan.success) {
            result.stopReason = "re-solve after cuts failed";
            if (options.traceRounds) {
                trace.incumbentAfter = result.Wstar;
                trace.incumbentChanged =
                    (trace.incumbentAfter + MILP_EPS < trace.incumbentBefore);
                result.roundTraces.push_back(trace);
            }
            break;
        }

        const double LBnext = rescan.Wproxy;
        result.LB = LBnext;

        if (options.traceRounds) {
            trace.LBnext = LBnext;
            trace.incumbentAfter = result.Wstar;
            trace.incumbentChanged =
                (trace.incumbentAfter + MILP_EPS < trace.incumbentBefore);
            result.roundTraces.push_back(trace);
        }

        if (std::isfinite(result.Wstar)) {
            const double gapCoarse = result.Wstar - LBnext;
            if (gapCoarse >= -MILP_EPS && gapCoarse <= MILP_EPS) {
                result.stoppedCoarse = true;
                result.provenOptimal = true;
                result.stopReason = "coarse test: W* - LB_{k+1} <= MILP_EPS";
                if (options.traceRounds && !result.roundTraces.empty())
                    result.roundTraces.back().coarseStop = true;
                break;
            }
            if (gapCoarse < -MILP_EPS) {
                result.stopReason = "invalid LB: W* < LB (gap="
                    + std::to_string(gapCoarse) + "), not proven";
            }
        }
    }

    if (!result.stopReason.empty() && !result.provenOptimal && !result.stoppedBudget
            && result.stopReason.find("invalid LB") == std::string::npos)
        result.stopReason = "max rounds exhausted";

    if (std::isfinite(result.Wstar) && std::isfinite(result.LB)) {
        result.gap = result.Wstar - result.LB;
        if (result.gap < -MILP_EPS && result.provenOptimal) {
            result.provenOptimal = false;
            result.stopReason += " [retracted: negative gap]";
        }
    }

    return result;
}
