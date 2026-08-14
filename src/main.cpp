// Proxy Master — benchmark runner (N=8, density=1).
// Usage:
//   ./proxy_master [step] [seed]
//   ./proxy_master batch10          — seeds 1..10 summary table

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "DemandMatrices.h"
#include "InstanceSetup.h"
#include "LShapeCrossing.h"
#include "MILPSolver.h"
#include "PhysicalConstants.h"
#include "ProxyMasterLoop.h"
#include "ProxyMasterMilp.h"
#include "ShortcutGrid.h"
#include "ShortcutMethodD.h"
#include "ShortcutMethods.h"
#include "StageE.h"
#include "WavelengthLoadBalance.h"
#include "PlotExport.h"

namespace {

constexpr int N = 8;
constexpr int DENSITY = 1;

struct Instance {
    std::vector<Node> nodes;
    DemandMatrix D;
    bool ok = false;
};

Instance loadInstance(int n, int seed) {
    Instance inst;
    if (!initDemandMatrix(n, DENSITY, inst.D))
        return inst;
    constexpr double spacing = 5.0 * ShortcutGrid::DEFAULT_S_MIN;
    inst.nodes = generateNodes(n, static_cast<unsigned>(seed), spacing);
    inst.ok = !inst.nodes.empty();
    return inst;
}

Instance loadInstance(int seed) {
    return loadInstance(N, seed);
}

double methodDOnLayout(
        const Instance& inst,
        const MILPSolveResult& layout) {
    ShortcutMethodOptions opt;
    opt.quiet = true;
    opt.skipExports = true;
    return runMethodDJointShortcuts(
        inst.nodes, inst.D, layout, ShortcutGrid::DEFAULT_S_MIN, opt).globalW;
}

struct Step1Result {
    bool ok = false;
    double LB1 = 0.0;
    double wTrueR1 = 0.0;
    double wBaseline = 0.0;
    double wHatR1 = 0.0;
    int xPairs = 0;
    int scR1 = 0;
};

Step1Result runStep1(const Instance& inst) {
    Step1Result r;
    ProxyMasterMilpOptions mOpt;
    mOpt.quiet = true;
    mOpt.poolSolutions = 1;
    const ProxyMasterMilpResult master = solveProxyMasterMilp(
        inst.nodes, inst.D, {}, mOpt);
    if (!master.success || master.pool.empty())
        return r;

    r.ok = true;
    r.LB1 = master.Wproxy;
    r.xPairs = master.numXPairs;
    r.scR1 = (int)master.pool.front().selectedShortcutIndices.size();
    r.wTrueR1 = methodDOnLayout(inst, master.pool.front().layout);
    r.wBaseline = runMethodDBaselineW(
        inst.nodes, inst.D, ShortcutGrid::DEFAULT_S_MIN, true);

    const std::vector<ShortcutPairIndex> pairs =
        buildMasterShortcutPairs(inst.nodes);
    r.wHatR1 = computeStageE(
        inst.nodes, inst.D, master.pool.front().layout,
        master.pool.front().selectedShortcutIndices, pairs,
        ShortcutGrid::DEFAULT_S_MIN);
    return r;
}

struct BatchRow {
    int seed = 0;
    bool ok = false;
    std::string note;
    double LB1 = 0.0;
    double wTrueR1 = 0.0;
    double wBaseline = 0.0;
    double wHatR1 = 0.0;
    double wStar = 0.0;
    double gap = 0.0;
    int poolSize = 0;
    int ringsEval = 0;
    int rounds = 0;
    bool proven = false;
    std::string stopReason;
};

struct StandaloneRow {
    int seed = 0;
    bool ok = false;
    std::string note;
    double wStar = 0.0;
    double lbFinal = 0.0;
    double gap = 0.0;
    int rounds = 0;
    int ringsEval = 0;
    bool proven = false;
    bool stoppedFine = false;
    bool stoppedCoarse = false;
    std::string stopReason;
    double wBaselinePost = 0.0;  // computed after loop, for comparison only
};

BatchRow runBatchSeed(int seed) {
    BatchRow row;
    row.seed = seed;
    const Instance inst = loadInstance(seed);
    if (!inst.ok) {
        row.note = "instance";
        return row;
    }

    const Step1Result s1 = runStep1(inst);
    if (!s1.ok) {
        row.note = "master";
        return row;
    }
    if (!std::isfinite(s1.wBaseline)) {
        row.note = "baseline ring";
        return row;
    }

    row.ok = true;
    row.LB1 = s1.LB1;
    row.wTrueR1 = s1.wTrueR1;
    row.wBaseline = s1.wBaseline;
    row.wHatR1 = s1.wHatR1;

    ProxyMasterMilpOptions poolOpt;
    poolOpt.quiet = true;
    poolOpt.poolSolutions = 10;
    const ProxyMasterMilpResult pool = solveProxyMasterMilp(
        inst.nodes, inst.D, {}, poolOpt);
    row.poolSize = pool.success ? (int)pool.pool.size() : 0;

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    row.wStar = loop.Wstar;
    row.gap = loop.gap;
    row.ringsEval = loop.ringsEvaluated;
    row.rounds = loop.rounds;
    row.proven = loop.provenOptimal;
    row.stopReason = loop.stopReason;
    return row;
}

StandaloneRow runStandaloneSeed(int seed) {
    StandaloneRow row;
    row.seed = seed;
    const Instance inst = loadInstance(seed);
    if (!inst.ok) {
        row.note = "instance";
        return row;
    }

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    row.ok = std::isfinite(loop.Wstar);
    if (!row.ok) {
        row.note = loop.stopReason.empty() ? "no incumbent" : loop.stopReason;
        return row;
    }

    row.wStar = loop.Wstar;
    row.lbFinal = loop.LB;
    row.gap = loop.gap;
    row.rounds = loop.rounds;
    row.ringsEval = loop.ringsEvaluated;
    row.proven = loop.provenOptimal;
    row.stoppedFine = loop.stoppedFine;
    row.stoppedCoarse = loop.stoppedCoarse;
    row.stopReason = loop.stopReason;
    row.wBaselinePost = runMethodDBaselineW(
        inst.nodes, inst.D, ShortcutGrid::DEFAULT_S_MIN, true);
    return row;
}

void reportBatch10Standalone() {
    std::cout << "=== Standalone loop batch (no W_base access during solve) ===\n";
    std::cout << "N=" << N << ", density=" << DENSITY
              << ", seeds 1..10 (skip seed 6 if instance-only)\n\n";
    std::cout << "seed | ok | W* | W_base(post) | match | rounds | D-eval | gap | proven | stop\n";
    std::cout << "-----|----|----|--------------|-------|--------|--------|-----|--------|----\n";

    std::vector<StandaloneRow> rows;
    for (int seed = 1; seed <= 10; ++seed) {
        if (seed == 6) continue;
        std::cout << "Running seed " << seed << "..." << std::endl;
        rows.push_back(runStandaloneSeed(seed));
    }

    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(2);
    int okCount = 0;
    int matchCount = 0;
    for (const StandaloneRow& r : rows) {
        std::cout << std::setw(4) << r.seed << " | "
                  << (r.ok ? " Y " : " N ") << " | ";
        if (!r.ok) {
            std::cout << "  -  |      -       |   -   |   -   |   -   |  -  |   N   | "
                      << r.note << "\n";
            continue;
        }
        ++okCount;
        const bool match = std::isfinite(r.wBaselinePost)
            && std::abs(r.wStar - r.wBaselinePost) <= MILP_EPS;
        if (match) ++matchCount;
        std::cout << std::setw(4) << r.wStar << " | "
                  << std::setw(12) << r.wBaselinePost << " | "
                  << (match ? "  Y  " : "  N  ") << " | "
                  << std::setw(6) << r.rounds << " | "
                  << std::setw(6) << r.ringsEval << " | "
                  << std::setw(4) << r.gap << " | "
                  << (r.proven ? "  Y   " : "  N   ") << " | "
                  << r.stopReason << "\n";
    }

    std::cout << "\n--- summary (" << okCount << " successful) ---\n";
    std::cout << "W* matches W_base (post-hoc): " << matchCount << "/" << okCount << "\n";
}

void reportBatch20Standalone() {
    std::cout << "=== Standalone loop batch seeds 1..20 ===\n";
    std::cout << "N=" << N << ", density=" << DENSITY << "\n\n";
    std::cout << "seed | ok | W* | W_base(post) | match | rounds | D-eval | gap | proven | stop\n";
    std::cout << "-----|----|----|--------------|-------|--------|--------|-----|--------|----\n";

    std::vector<StandaloneRow> rows;
    for (int seed = 1; seed <= 20; ++seed) {
        std::cout << "Running seed " << seed << "..." << std::endl;
        rows.push_back(runStandaloneSeed(seed));
    }

    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(2);
    int okCount = 0, provenCount = 0, matchCount = 0;
    int hitRoundLimit = 0;
    for (const StandaloneRow& r : rows) {
        std::cout << std::setw(4) << r.seed << " | "
                  << (r.ok ? " Y " : " N ") << " | ";
        if (!r.ok) {
            std::cout << "  -  |      -       |   -   |   -   |   -   |  -  |   N   | "
                      << r.note << "\n";
            continue;
        }
        ++okCount;
        if (r.proven) ++provenCount;
        const bool match = std::isfinite(r.wBaselinePost)
            && std::abs(r.wStar - r.wBaselinePost) <= MILP_EPS;
        if (match) ++matchCount;
        if (r.stopReason.find("max rounds") != std::string::npos
                || r.stopReason.find("exhausted") != std::string::npos)
            ++hitRoundLimit;
        std::cout << std::setw(4) << r.wStar << " | "
                  << std::setw(12) << (std::isfinite(r.wBaselinePost)
                      ? std::to_string(r.wBaselinePost) : "-") << " | "
                  << std::setw(5) << (match ? "Y" : (std::isfinite(r.wBaselinePost) ? "N" : "n/a")) << " | "
                  << std::setw(6) << r.rounds << " | "
                  << std::setw(6) << r.ringsEval << " | "
                  << std::setw(4) << r.gap << " | "
                  << (r.proven ? "  Y   " : "  N   ") << " | "
                  << r.stopReason << "\n";
    }
    std::cout << "\n--- summary ---\n";
    std::cout << "successful: " << okCount << "/20, proven=Y: " << provenCount
              << ", W*==W_base: " << matchCount << "/" << okCount
              << ", hit round limit: " << hitRoundLimit << "\n";
}

struct N12BatchRow {
    int seed = 0;
    bool proxyOk = false;
    bool baseOk = false;
    double LB1 = 0.0;
    double wStar = 0.0;
    double wBase = std::numeric_limits<double>::quiet_NaN();
    double gap = 0.0;
    bool proven = false;
    int rounds = 0;
    int ringsEval = 0;
    std::string stopReason;
    std::string note;
};

N12BatchRow runN12Seed(int seed) {
    constexpr int nNodes = 12;
    N12BatchRow row;
    row.seed = seed;
    const Instance inst = loadInstance(nNodes, seed);
    if (!inst.ok) {
        row.note = "instance";
        return row;
    }

    ProxyMasterMilpOptions mOpt;
    mOpt.quiet = true;
    mOpt.poolSolutions = 1;
    const ProxyMasterMilpResult master0 = solveProxyMasterMilp(
        inst.nodes, inst.D, {}, mOpt);
    if (master0.success)
        row.LB1 = master0.Wproxy;

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    row.proxyOk = std::isfinite(loop.Wstar);
    if (!row.proxyOk) {
        row.note = loop.stopReason.empty() ? "no incumbent" : loop.stopReason;
        return row;
    }

    row.wStar = loop.Wstar;
    row.gap = loop.gap;
    row.rounds = loop.rounds;
    row.ringsEval = loop.ringsEvaluated;
    row.proven = loop.provenOptimal;
    row.stopReason = loop.stopReason;

    row.wBase = runMethodDBaselineW(
        inst.nodes, inst.D, ShortcutGrid::DEFAULT_S_MIN, true);
    row.baseOk = std::isfinite(row.wBase);
    return row;
}

double medianInt(const std::vector<int>& v) {
    if (v.empty()) return 0.0;
    std::vector<int> s = v;
    std::sort(s.begin(), s.end());
    const size_t n = s.size();
    if (n % 2 == 1) return static_cast<double>(s[n / 2]);
    return 0.5 * (s[n / 2 - 1] + s[n / 2]);
}

void reportBatch12Standalone() {
    std::cout << "=== N=12 standalone proxy master batch ===\n";
    std::cout << "N=12, density=" << DENSITY << ", seeds 1..10\n\n";
    std::cout << "seed | LB_1 | W* | W_base | gap | proven | rounds | D-evals | stop\n";
    std::cout << "-----|------|----|--------|-----|--------|--------|---------|----\n";

    std::vector<N12BatchRow> rows;
    for (int seed = 1; seed <= 10; ++seed) {
        std::cout << "Running seed " << seed << "..." << std::endl;
        rows.push_back(runN12Seed(seed));
    }

    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(2);

    int proxyOk = 0, provenCount = 0, roundLimit = 0;
    int better = 0, equal = 0, worse = 0;
    std::vector<int> roundsList, evalsList;
    std::vector<std::pair<int, double>> positiveGaps;

    for (const N12BatchRow& r : rows) {
        std::cout << std::setw(4) << r.seed << " | ";
        if (!r.proxyOk) {
            std::cout << "  -  |  -  |   -   |   -  |   N   |   -   |    -   | "
                      << r.note << "\n";
            continue;
        }
        ++proxyOk;
        if (r.proven) ++provenCount;
        if (r.stopReason.find("exhausted") != std::string::npos
                || r.stopReason.find("max rounds") != std::string::npos)
            ++roundLimit;
        roundsList.push_back(r.rounds);
        evalsList.push_back(r.ringsEval);
        if (r.gap > MILP_EPS) positiveGaps.push_back({r.seed, r.gap});

        if (r.baseOk) {
            if (r.wStar < r.wBase - MILP_EPS) ++better;
            else if (r.wStar > r.wBase + MILP_EPS) ++worse;
            else ++equal;
        }

        std::cout << std::setw(4) << r.LB1 << " | "
                  << std::setw(4) << r.wStar << " | "
                  << (r.baseOk ? std::to_string(r.wBase) : std::string("-")) << " | "
                  << std::setw(4) << r.gap << " | "
                  << (r.proven ? "  Y   " : "  N   ") << " | "
                  << std::setw(6) << r.rounds << " | "
                  << std::setw(7) << r.ringsEval << " | "
                  << r.stopReason;
        if (!r.baseOk) std::cout << " [baseline-fail]";
        if (r.proven && r.baseOk && r.wStar > r.wBase + MILP_EPS)
            std::cout << " [W*>W_base while proven!]";
        std::cout << "\n";
    }

    std::cout << "\n--- aggregate (seeds 1-10, N=12) ---\n";
    std::cout << "proxy successful: " << proxyOk << "/10\n";
    std::cout << "proven=Y: " << provenCount << "/10, hit round limit: " << roundLimit << "\n";
    std::cout << "W* vs W_base (where baseline ok): better=" << better
              << ", equal=" << equal << ", worse=" << worse << "\n";
    if (!roundsList.empty()) {
        std::cout << "rounds: min=" << *std::min_element(roundsList.begin(), roundsList.end())
                  << ", median=" << medianInt(roundsList)
                  << ", max=" << *std::max_element(roundsList.begin(), roundsList.end())
                  << " (N=8 after cut fix: min=1, median=1, max=5)\n";
        std::cout << "D-evals: min=" << *std::min_element(evalsList.begin(), evalsList.end())
                  << ", median=" << medianInt(evalsList)
                  << ", max=" << *std::max_element(evalsList.begin(), evalsList.end())
                  << " (N=8 after cut fix: min=1, median=1, max=5)\n";
    }
    if (positiveGaps.empty())
        std::cout << "seeds with gap > 0: none\n";
    else {
        std::cout << "seeds with gap > 0:\n";
        for (const auto& [s, g] : positiveGaps)
            std::cout << "  seed " << s << ": gap=" << g << "\n";
    }
}

struct ExtendedRow {
    int seed = 0;
    bool proxyOk = false;
    double wStar = 0.0;
    int rounds = 0;
    int ringsEval = 0;
    double gap = 0.0;
    bool proven = false;
    std::string stopKind;
    std::string stopReason;
    double wMethodB = std::numeric_limits<double>::quiet_NaN();
    double wBase = std::numeric_limits<double>::quiet_NaN();
    bool methodBOk = false;
    bool baseOk = false;
    std::string note;
};

double runMethodBW(
        const Instance& inst,
        bool* okOut = nullptr) {
    MILPSolver solverB(inst.nodes, inst.D);
    MILPSolveResult layoutB = solverB.solve(
        true, "", {}, -1.0, nullptr, -1.0, true);
    if (!layoutB.success) {
        if (okOut) *okOut = false;
        return std::numeric_limits<double>::infinity();
    }
    ShortcutMethodOptions opt;
    opt.quiet = true;
    opt.skipExports = true;
    const ShortcutMethodResult res = runMethodBWithShortcuts(
        inst.nodes, inst.D, solverB, layoutB,
        ShortcutGrid::DEFAULT_S_MIN, "B", "Method B+Shortcuts", opt);
    if (okOut) *okOut = res.layout.success && std::isfinite(res.globalW);
    return res.globalW;
}

std::string compareVsMethodB(double wStar, double wB) {
    if (!std::isfinite(wB)) return "n/a";
    if (wStar < wB - MILP_EPS) return "better";
    if (wStar > wB + MILP_EPS) return "worse";
    return "equal";
}

double medianOf(std::vector<int> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2 == 1) return static_cast<double>(v[n / 2]);
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

ExtendedRow runExtendedSeed(int seed) {
    ExtendedRow row;
    row.seed = seed;
    const Instance inst = loadInstance(seed);
    if (!inst.ok) {
        row.note = "instance";
        return row;
    }

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    row.proxyOk = std::isfinite(loop.Wstar);
    if (row.proxyOk) {
        row.wStar = loop.Wstar;
        row.rounds = loop.rounds;
        row.ringsEval = loop.ringsEvaluated;
        row.gap = loop.gap;
        row.proven = loop.provenOptimal;
        row.stopReason = loop.stopReason;
        if (loop.stoppedFine) row.stopKind = "fine";
        else if (loop.stoppedCoarse) row.stopKind = "coarse";
        else if (loop.stoppedBudget) row.stopKind = "budget";
        else row.stopKind = "other";
    } else {
        row.note = loop.stopReason.empty() ? "proxy failed" : loop.stopReason;
    }

    row.wBase = runMethodDBaselineW(
        inst.nodes, inst.D, ShortcutGrid::DEFAULT_S_MIN, true);
    row.baseOk = std::isfinite(row.wBase);
    row.wMethodB = runMethodBW(inst, &row.methodBOk);
    return row;
}

void reportBatch1120Extended() {
    std::cout << "=== Extended batch: standalone proxy master + Method B ===\n";
    std::cout << "N=" << N << ", density=" << DENSITY << ", seeds 11..20\n\n";
    std::cout << "seed | W* | rounds | D-evals | gap | proven | stop | W_B | W_base | proxy==W_base | proxy vs Method B\n";
    std::cout << "-----|----|--------|---------|-----|--------|------|-----|--------|---------------|------------------\n";

    std::vector<ExtendedRow> rows;
    for (int seed = 11; seed <= 20; ++seed) {
        std::cout << "Running seed " << seed << "..." << std::endl;
        rows.push_back(runExtendedSeed(seed));
    }

    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(2);

    int proxyOkCount = 0;
    int provenCount = 0;
    int matchBaseCount = 0;
    int betterThanB = 0, equalToB = 0, worseThanB = 0;
    std::vector<int> roundsList;
    std::vector<int> evalsList;
    std::vector<std::pair<int, double>> positiveGaps;

    for (const ExtendedRow& r : rows) {
        std::cout << std::setw(4) << r.seed << " | ";
        if (!r.proxyOk) {
            std::cout << "  -  |   -   |    -   |   -  |   N   |  -   |   -   |   -   |    -     |      -        | "
                      << r.note << "\n";
            continue;
        }
        ++proxyOkCount;
        if (r.proven) ++provenCount;
        roundsList.push_back(r.rounds);
        evalsList.push_back(r.ringsEval);
        if (r.gap > MILP_EPS)
            positiveGaps.push_back({r.seed, r.gap});

        const bool matchBase = r.baseOk
            && std::abs(r.wStar - r.wBase) <= MILP_EPS;
        if (matchBase) ++matchBaseCount;

        const std::string vsB = compareVsMethodB(r.wStar, r.wMethodB);
        if (vsB == "better") ++betterThanB;
        else if (vsB == "equal") ++equalToB;
        else if (vsB == "worse") ++worseThanB;

        std::cout << std::setw(4) << r.wStar << " | "
                  << std::setw(6) << r.rounds << " | "
                  << std::setw(7) << r.ringsEval << " | "
                  << std::setw(4) << r.gap << " | "
                  << (r.proven ? "  Y   " : "  N   ") << " | "
                  << std::setw(4) << r.stopKind << " | "
                  << (r.methodBOk ? std::to_string(r.wMethodB) : std::string("-")) << " | "
                  << (r.baseOk ? std::to_string(r.wBase) : std::string("-")) << " | "
                  << std::setw(13) << (matchBase ? "Y" : (r.baseOk ? "N" : "n/a")) << " | "
                  << vsB;
        if (!r.baseOk || !r.methodBOk)
            std::cout << " [" << (r.baseOk ? "" : "base-fail ")
                      << (r.methodBOk ? "" : "methodB-fail") << "]";
        std::cout << "\n";
    }

    std::cout << "\n--- aggregate summary (seeds 11-20) ---\n";
    std::cout << "proxy runs successful: " << proxyOkCount << "/10\n";
    std::cout << "proven=Y: " << provenCount << "/" << proxyOkCount << "\n";
    std::cout << "proxy == W_base: " << matchBaseCount << "/" << proxyOkCount << "\n";
    std::cout << "proxy vs Method B: better=" << betterThanB
              << ", equal=" << equalToB << ", worse=" << worseThanB << "\n";

    if (positiveGaps.empty()) {
        std::cout << "seeds with gap > 0: none\n";
    } else {
        std::cout << "seeds with gap > 0:\n";
        for (const auto& [seed, gap] : positiveGaps)
            std::cout << "  seed " << seed << ": gap=" << gap << "\n";
    }

    if (!roundsList.empty()) {
        std::cout << "rounds: min=" << *std::min_element(roundsList.begin(), roundsList.end())
                  << ", median=" << medianOf(roundsList)
                  << ", max=" << *std::max_element(roundsList.begin(), roundsList.end()) << "\n";
        std::cout << "D-evals: min=" << *std::min_element(evalsList.begin(), evalsList.end())
                  << ", median=" << medianOf(evalsList)
                  << ", max=" << *std::max_element(evalsList.begin(), evalsList.end()) << "\n";
    }
}

void reportStep1(const Instance& inst, int seed) {
    std::cout << "=== Step 1 checkpoint (N=" << N << ", density=" << DENSITY
              << ", seed=" << seed << ") ===\n";
    const Step1Result s1 = runStep1(inst);
    if (!s1.ok) {
        std::cout << "Master FAILED\n";
        return;
    }
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "LB_1 = W_proxy*     : " << s1.LB1 << " mm\n";
    std::cout << "W_true(R_1)        : " << s1.wTrueR1 << " mm\n";
    std::cout << "W_true(baseline B) : " << s1.wBaseline << " mm\n";
    std::cout << "W_hat(R_1)         : " << s1.wHatR1 << " mm\n";
    std::cout << "XPairs (M2b)       : " << s1.xPairs << "\n";
    std::cout << "Shortcuts in R_1   : " << s1.scR1 << "\n";
}

void reportWavelengthLoadBalance(int nNodes, int seed) {
    std::cout << "=== Wavelength load-balance post-process ===\n";
    std::cout << "N=" << nNodes << ", density=" << DENSITY
              << ", seed=" << seed << "\n\n";

    const Instance inst = loadInstance(nNodes, seed);
    if (!inst.ok) {
        std::cerr << "instance load failed\n";
        return;
    }

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = (nNodes >= 12) ? 600.0 : 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Loop: W*=" << loop.Wstar
              << " LB=" << loop.LB
              << " gap=" << loop.gap
              << " proven=" << (loop.provenOptimal ? "Y" : "N")
              << " rounds=" << loop.rounds
              << " D-evals=" << loop.ringsEvaluated
              << "\n  stop: " << loop.stopReason << "\n";

    if (!std::isfinite(loop.Wstar) || !loop.best.layout.success) {
        std::cerr << "No usable loop incumbent — skip WLB\n";
        return;
    }

    // Prefer the loop's best ShortcutMethodResult; ensure tour is present.
    ShortcutMethodResult best = loop.best;
    if (best.layout.tour.empty() && !loop.best.layout.tourEdges.empty()) {
        // Reconstruct tour order from edges if needed — normally set by Method D.
    }

    const WavelengthLoadBalanceResult wlb = runWavelengthLoadBalance(
        inst.nodes, inst.D, best, loop.Wstar);
    printWavelengthLoadBalanceReport(wlb, inst.D);

    exportProxyMasterPlotCsvs(inst.nodes, inst.D, best, &wlb, "PM");
    std::cout << "\nCSV exports written (nodes.csv, ring_PM.csv, ...)\n";
    std::cout << "Plot with: python3 tools/plot_final.py --suffix PM --open\n";
}

void reportPlot(int nNodes, int seed) {
    std::cout << "=== Final plot export: N=" << nNodes
              << " density=" << DENSITY << " seed=" << seed << " ===\n\n";
    const Instance inst = loadInstance(nNodes, seed);
    if (!inst.ok) {
        std::cerr << "instance load failed\n";
        return;
    }

    ProxyMasterLoopOptions loopOpt;
    loopOpt.quiet = true;
    loopOpt.poolSolutions = 10;
    loopOpt.maxRounds = 20;
    loopOpt.wallTimeLimitSec = (nNodes >= 12) ? 600.0 : 300.0;
    const ProxyMasterLoopResult loop = runProxyMasterLoop(
        inst.nodes, inst.D, loopOpt);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Loop: W*=" << loop.Wstar
              << " proven=" << (loop.provenOptimal ? "Y" : "N")
              << " stop=" << loop.stopReason << "\n";

    if (!std::isfinite(loop.Wstar) || !loop.best.layout.success) {
        std::cerr << "No usable incumbent\n";
        return;
    }

    ShortcutMethodResult best = loop.best;
    const WavelengthLoadBalanceResult wlb = runWavelengthLoadBalance(
        inst.nodes, inst.D, best, loop.Wstar);
    exportProxyMasterPlotCsvs(inst.nodes, inst.D, best, &wlb, "PM");
    std::cout << "Exported CSVs for suffix PM.\n";
    std::cout << "Shortcuts: " << best.shortcuts.size() << "\n";
}

void reportStep4(const Instance& inst) {
    std::cout << "\n=== Step 4 full loop ===\n";
    ProxyMasterLoopOptions opt;
    opt.quiet = true;
    opt.poolSolutions = 10;
    opt.maxRounds = 20;
    opt.wallTimeLimitSec = 300.0;
    const ProxyMasterLoopResult res = runProxyMasterLoop(inst.nodes, inst.D, opt);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "W* (incumbent)  : " << res.Wstar << " mm\n";
    std::cout << "LB (final)      : " << res.LB << " mm\n";
    std::cout << "gap             : " << res.gap << " mm\n";
    std::cout << "rounds          : " << res.rounds << "\n";
    std::cout << "rings evaluated : " << res.ringsEvaluated << "\n";
    std::cout << "proven optimal  : " << (res.provenOptimal ? "yes" : "no") << "\n";
    std::cout << "stop reason     : " << res.stopReason << "\n";
}

void reportBenchInc(int nNodes, int seed) {
    std::cout << "=== Bench incremental master: N=" << nNodes
              << " dens=" << DENSITY << " seed=" << seed << " ===\n\n";
    const Instance inst = loadInstance(nNodes, seed);
    if (!inst.ok) {
        std::cerr << "instance failed\n";
        return;
    }

    struct Row {
        bool incremental = false;
        double sec = 0.0;
        double wStar = 0.0;
        double lb = 0.0;
        int rounds = 0;
        int evals = 0;
        bool proven = false;
        std::string stop;
    };
    std::vector<Row> rows;
    for (bool inc : {false, true}) {
        ProxyMasterLoopOptions opt;
        opt.quiet = true;
        opt.poolSolutions = 10;
        opt.maxRounds = 20;
        opt.wallTimeLimitSec = 300.0;
        opt.incrementalMaster = inc;
        const auto t0 = std::chrono::steady_clock::now();
        const ProxyMasterLoopResult res = runProxyMasterLoop(inst.nodes, inst.D, opt);
        const double sec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
        Row r;
        r.incremental = inc;
        r.sec = sec;
        r.wStar = res.Wstar;
        r.lb = res.LB;
        r.rounds = res.rounds;
        r.evals = res.ringsEvaluated;
        r.proven = res.provenOptimal;
        r.stop = res.stopReason;
        rows.push_back(r);
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "mode         | sec    | W*     | LB     | rounds | D-eval | proven | stop\n";
    std::cout << "-------------|--------|--------|--------|--------|--------|--------|-----\n";
    for (const Row& r : rows) {
        std::cout << (r.incremental ? "incremental " : "rebuild     ")
                  << " | " << std::setw(6) << r.sec
                  << " | " << std::setw(6) << r.wStar
                  << " | " << std::setw(6) << r.lb
                  << " | " << std::setw(6) << r.rounds
                  << " | " << std::setw(6) << r.evals
                  << " | " << (r.proven ? "Y" : "N")
                  << "      | " << r.stop << "\n";
    }
    if (rows.size() == 2) {
        const bool sameW = std::abs(rows[0].wStar - rows[1].wStar) <= 1e-6
            && std::abs(rows[0].lb - rows[1].lb) <= 1e-6
            && rows[0].proven == rows[1].proven;
        std::cout << "\nSame result: " << (sameW ? "YES" : "NO") << "\n";
        if (rows[0].sec > 1e-9) {
            std::cout << "Speedup (rebuild/inc): "
                      << (rows[0].sec / rows[1].sec) << "x\n";
        }
    }
}

void reportBatch10() {
    std::cout << "=== Batch: N=" << N << ", density=" << DENSITY
              << ", seeds 1..10 ===\n\n";
    std::cout << "seed | ok | LB_1 | W_true(R1) | W_base | W_hat(R1) | W* | gap | pool | rings | proven | note\n";
    std::cout << "-----|----|------|------------|--------|-----------|----|-----|------|-------|--------|-----\n";

    std::vector<BatchRow> rows;
    for (int seed = 1; seed <= 10; ++seed) {
        std::cout << "Running seed " << seed << "..." << std::endl;
        rows.push_back(runBatchSeed(seed));
    }

    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(2);
    int okCount = 0;
    int r1Better = 0, loopBetter = 0, tieBase = 0;
    for (const BatchRow& r : rows) {
        std::cout << std::setw(4) << r.seed << " | "
                  << (r.ok ? " Y " : " N ") << " | ";
        if (!r.ok) {
            std::cout << "  -   |     -      |   -    |     -     |  - |  -  |  -   |   -   | "
                      << r.note << "\n";
            continue;
        }
        ++okCount;
        if (r.wTrueR1 < r.wBaseline - MILP_EPS) ++r1Better;
        else if (r.wTrueR1 > r.wBaseline + MILP_EPS) {}
        else ++tieBase;
        if (r.wStar < r.wBaseline - MILP_EPS) ++loopBetter;

        std::cout << std::setw(5) << r.LB1 << " | "
                  << std::setw(10) << r.wTrueR1 << " | "
                  << std::setw(6) << r.wBaseline << " | "
                  << std::setw(9) << r.wHatR1 << " | "
                  << std::setw(4) << r.wStar << " | "
                  << std::setw(4) << r.gap << " | "
                  << std::setw(4) << r.poolSize << " | "
                  << std::setw(5) << r.ringsEval << " | "
                  << (r.proven ? "  Y   " : "  N   ") << " | \n";
    }

    std::cout << "\n--- aggregated (" << okCount << " successful) ---\n";
    std::cout << "R1 better than baseline: " << r1Better << "/" << okCount << "\n";
    std::cout << "Loop W* better than baseline: " << loopBetter << "/" << okCount << "\n";
    std::cout << "R1 tie with baseline: " << tieBase << "/" << okCount << "\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string step = "1";
    int seed = 1;
    if (argc > 1) step = argv[1];
    if (argc > 2) seed = std::stoi(argv[2]);

    if (step == "batch10") {
        reportBatch10();
        return 0;
    }

    if (step == "batch10_standalone") {
        reportBatch10Standalone();
        return 0;
    }

    if (step == "batch20_standalone") {
        reportBatch20Standalone();
        return 0;
    }

    if (step == "batch12_standalone") {
        reportBatch12Standalone();
        return 0;
    }

    if (step == "plot") {
        int nNodes = N;
        int plotSeed = 1;
        if (argc == 3) {
            plotSeed = std::stoi(argv[2]);
        } else if (argc >= 4) {
            nNodes = std::stoi(argv[2]);
            plotSeed = std::stoi(argv[3]);
        }
        reportPlot(nNodes, plotSeed);
        return 0;
    }

    if (step == "wlb") {
        int nNodes = N;
        if (argc > 2) nNodes = std::stoi(argv[2]);
        int wlbSeed = 1;
        if (argc > 3) wlbSeed = std::stoi(argv[3]);
        else if (argc > 2 && nNodes <= 32) {
            // allow: ./proxy_master wlb <seed>  with default N=8
            // or:    ./proxy_master wlb <N> <seed>
        }
        // Parse: wlb [N] [seed]  — if only one int after wlb and it's small seed-like
        // Prefer: wlb N seed
        if (argc == 3) {
            // single arg: treat as seed with default N=8
            wlbSeed = std::stoi(argv[2]);
            nNodes = N;
        } else if (argc >= 4) {
            nNodes = std::stoi(argv[2]);
            wlbSeed = std::stoi(argv[3]);
        }
        reportWavelengthLoadBalance(nNodes, wlbSeed);
        return 0;
    }

    if (step == "batch1120") {
        reportBatch1120Extended();
        return 0;
    }

    if (step == "bench_inc") {
        int nNodes = N;
        int benchSeed = seed;
        if (argc >= 4) {
            nNodes = std::stoi(argv[2]);
            benchSeed = std::stoi(argv[3]);
        } else if (argc == 3) {
            benchSeed = std::stoi(argv[2]);
        }
        reportBenchInc(nNodes, benchSeed);
        return 0;
    }

    const Instance inst = loadInstance(seed);
    if (!inst.ok) return 1;

    if (step == "1" || step == "all") reportStep1(inst, seed);
    if (step == "4" || step == "all") reportStep4(inst);

    if (step != "1" && step != "4" && step != "all") {
        std::cerr << "Unknown step: " << step
                  << " (use 1, 4, all, batch10, batch10_standalone, batch12_standalone,"
                  << " batch20_standalone, batch1120, bench_inc [N] [seed],"
                  << " wlb [N] [seed])\n";
        return 1;
    }
    return 0;
}
