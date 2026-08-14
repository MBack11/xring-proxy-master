// Overnight benchmark batch: proxy master vs W_base on demand matrices (a)–(d).
// Usage:
//   ./benchmark_batch              — full run seeds 1..20, all cases
//   ./benchmark_batch estimate     — one seed (1) per case, time estimate only
//   ./benchmark_batch a 1 3        — case a, seeds 1..3
//
// Seed range: 1..20 (placement via generateNodes; demands fixed per case).
// Wavelength: post-process WLB only (directional shortcut conflict fix required).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "DemandMatrix.h"
#include "InstanceSetup.h"
#include "MILPSolver.h"
#include "Nodes.h"
#include "ProxyMasterLoop.h"
#include "ProxyMasterMilp.h"
#include "ShortcutGrid.h"
#include "ShortcutMethodD.h"
#include "ShortcutMethods.h"
#include "WavelengthLoadBalance.h"

namespace {

using Clock = std::chrono::steady_clock;

struct BenchCase {
    char id = 'a';
    std::string slug;
    std::string title;
    int N = 0;
    std::vector<std::string> nodeNames;
    DemandMatrix D;
};

struct SeedRow {
    int seed = 0;
    int N = 0;
    double wStar = std::numeric_limits<double>::quiet_NaN();
    double lb1 = std::numeric_limits<double>::quiet_NaN();
    int rounds = 0;
    int dEvals = 0;
    double gap = std::numeric_limits<double>::quiet_NaN();
    bool proven = false;
    std::string stopReason;
    double proxySec = 0.0;
    int lambdaProxy = -1;
    int lambdaLbProxy = -1;
    bool wlbProxyOk = false;
    double wBase = std::numeric_limits<double>::quiet_NaN();
    double baseSec = 0.0;
    bool baseOk = false;
    int lambdaBase = -1;
    int lambdaLbBase = -1;
    bool wlbBaseOk = false;
    std::string wCmp;   // better / equal / worse / n/a
    std::string lamCmp; // better / equal / worse / n/a
    bool negativeGap = false;
    bool provenWorseThanBase = false;
};

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r')) --b;
    return s.substr(a, b - a);
}

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == ',') {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(trim(cur));
    return out;
}

bool loadCase(
        const std::string& root,
        char id,
        const std::string& slug,
        const std::string& title,
        BenchCase& out) {
    out.id = id;
    out.slug = slug;
    out.title = title;
    const std::string nodesPath = root + "/benchmarks/parsed/" + slug + "_nodes.csv";
    const std::string demPath = root + "/benchmarks/parsed/" + slug + "_demands.csv";

    std::ifstream nf(nodesPath);
    if (!nf) {
        std::cerr << "missing " << nodesPath << "\n";
        return false;
    }
    std::string line;
    std::getline(nf, line);  // header
    out.nodeNames.clear();
    while (std::getline(nf, line)) {
        if (trim(line).empty()) continue;
        auto cols = splitCsv(line);
        if (cols.size() < 2) continue;
        const int idNum = std::stoi(cols[0]);
        if (idNum != (int)out.nodeNames.size()) {
            std::cerr << "node id gap in " << nodesPath << "\n";
            return false;
        }
        out.nodeNames.push_back(cols[1]);
    }
    out.N = (int)out.nodeNames.size();

    std::ifstream df(demPath);
    if (!df) {
        std::cerr << "missing " << demPath << "\n";
        return false;
    }
    std::getline(df, line);
    out.D.demands.clear();
    while (std::getline(df, line)) {
        if (trim(line).empty()) continue;
        auto cols = splitCsv(line);
        if (cols.size() < 2) continue;
        const int s = std::stoi(cols[0]);
        const int t = std::stoi(cols[1]);
        out.D.add(s, t);
    }
    return out.N >= 3 && !out.D.demands.empty();
}

std::string mapStop(const ProxyMasterLoopResult& r) {
    if (r.stoppedFine) return "fine";
    if (r.stoppedCoarse) return "coarse";
    if (r.stoppedBudget) return "wall-time-limit";
    if (r.stopReason.find("max rounds") != std::string::npos) return "round-limit";
    if (!r.stopReason.empty()) return r.stopReason;
    return "unknown";
}

std::string cmpW(double wStar, double wBase, bool baseOk) {
    if (!baseOk || !std::isfinite(wBase) || !std::isfinite(wStar)) return "n/a";
    constexpr double eps = 1e-6;
    if (wStar + eps < wBase) return "better";
    if (wBase + eps < wStar) return "worse";
    return "equal";
}

std::string cmpLam(int a, int b, bool okA, bool okB) {
    if (!okA || !okB || a < 0 || b < 0) return "n/a";
    if (a < b) return "better";
    if (a > b) return "worse";
    return "equal";
}

ShortcutMethodResult runBaselineFull(
        const std::vector<Node>& nodes,
        const DemandMatrix& D,
        double sMin,
        bool* ok) {
    ShortcutMethodResult empty;
    *ok = false;
    MILPSolver solverA(nodes, D);
    MILPSolveResult layoutA = solverA.solve(false, "", {}, -1.0, nullptr, -1.0, true);
    if (!layoutA.success) return empty;
    MILPSolver solverB(nodes, D);
    MILPSolveResult layoutB = solverB.solve(true, "", {}, -1.0, nullptr, -1.0, true);
    const MILPSolveResult& ring = layoutB.success ? layoutB : layoutA;
    ShortcutMethodOptions opt;
    opt.quiet = true;
    opt.skipExports = true;
    ShortcutMethodResult res = runMethodDJointShortcuts(nodes, D, ring, sMin, opt);
    *ok = std::isfinite(res.globalW);
    return res;
}

SeedRow runSeed(const BenchCase& bc, int seed, const std::string& /*root*/) {
    SeedRow row;
    row.seed = seed;
    row.N = bc.N;

    constexpr double spacing = 5.0 * ShortcutGrid::DEFAULT_S_MIN;
    const std::vector<Node> nodes =
        generateNodes(bc.N, static_cast<unsigned>(seed), spacing);
    const DemandMatrix& D = bc.D;
    const double sMin = ShortcutGrid::DEFAULT_S_MIN;

    // --- Proxy master ---
    ProxyMasterLoopOptions opt;
    opt.quiet = true;
    opt.poolSolutions = 10;
    opt.maxRounds = 20;
    opt.wallTimeLimitSec = (bc.N >= 12) ? 600.0 : 300.0;
    opt.traceRounds = true;
    opt.incrementalMaster = true;
    opt.sMin = sMin;

    const auto t0 = Clock::now();
    const ProxyMasterLoopResult loop = runProxyMasterLoop(nodes, D, opt);
    row.proxySec = std::chrono::duration<double>(Clock::now() - t0).count();

    row.wStar = loop.Wstar;
    row.rounds = loop.rounds;
    row.dEvals = loop.ringsEvaluated;
    row.gap = loop.gap;
    row.proven = loop.provenOptimal;
    row.stopReason = mapStop(loop);
    row.lb1 = loop.roundTraces.empty() ? loop.LB : loop.roundTraces.front().LBk;
    if (std::isfinite(row.gap) && row.gap < -1e-6) row.negativeGap = true;

    if (opt.traceRounds && !loop.roundTraces.empty()) {
        std::cout << "  --- round trace ---\n";
        for (const auto& tr : loop.roundTraces) {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "  round " << tr.round
                      << " LB_k=" << tr.LBk;
            if (std::isfinite(tr.LBnext))
                std::cout << " LB_{k+1}=" << tr.LBnext;
            std::cout << " W* before=" << tr.incumbentBefore
                      << " W* after=" << tr.incumbentAfter
                      << " cuts+" << tr.cutsAdded
                      << (tr.fineStop ? " [fine]" : "")
                      << (tr.coarseStop ? " [coarse]" : "")
                      << "\n";
            for (const auto& sv : tr.stageVEvals) {
                std::cout << "    D-eval W_true=" << sv.wTrue
                          << (sv.fineTestFired ? " (fine-test)" : "")
                          << " ring=" << sv.ringKey << "\n";
            }
        }
        std::cout << "  final LB=" << loop.LB
                  << " W*=" << loop.Wstar
                  << " gap(W*-LB)=" << loop.gap << "\n";
    }

    if (std::isfinite(loop.Wstar) && loop.best.layout.success) {
        const WavelengthLoadBalanceResult wlb =
            runWavelengthLoadBalance(nodes, D, loop.best, loop.Wstar);
        row.wlbProxyOk = wlb.success && wlb.allRoutesRespectWstar;
        if (wlb.success) {
            row.lambdaProxy = wlb.lambdaNeeded;
            row.lambdaLbProxy = wlb.conflictLoadLowerBound;
        }
    }

    // --- W_base ---
    const auto t1 = Clock::now();
    bool baseOk = false;
    ShortcutMethodResult baseRes = runBaselineFull(nodes, D, sMin, &baseOk);
    row.baseSec = std::chrono::duration<double>(Clock::now() - t1).count();
    row.baseOk = baseOk;
    if (baseOk) {
        row.wBase = baseRes.globalW;
        const WavelengthLoadBalanceResult wlbB =
            runWavelengthLoadBalance(nodes, D, baseRes, baseRes.globalW);
        row.wlbBaseOk = wlbB.success && wlbB.allRoutesRespectWstar;
        if (wlbB.success) {
            row.lambdaBase = wlbB.lambdaNeeded;
            row.lambdaLbBase = wlbB.conflictLoadLowerBound;
        }
    }

    row.wCmp = cmpW(row.wStar, row.wBase, row.baseOk);
    row.lamCmp = cmpLam(row.lambdaProxy, row.lambdaBase, row.wlbProxyOk, row.wlbBaseOk);
    if (row.proven && row.baseOk && std::isfinite(row.wStar) && std::isfinite(row.wBase)
            && row.wStar > row.wBase + 1e-6) {
        row.provenWorseThanBase = true;
    }
    return row;
}

void writeCsvHeader(std::ostream& os) {
    os << "seed,N,W_star_proxy,LB1,rounds,D_evals,gap,proven,stop_reason,"
          "proxy_time_s,lambda_needed_proxy,lambda_LB_proxy,"
          "W_base,W_base_time_s,lambda_needed_base,lambda_LB_base,"
          "W_star_vs_W_base,lambda_proxy_vs_lambda_base\n";
}

void writeCsvRow(std::ostream& os, const SeedRow& r) {
    auto num = [](double x) -> std::string {
        if (!std::isfinite(x)) return "";
        std::ostringstream s;
        s << std::fixed << std::setprecision(6) << x;
        return s.str();
    };
    os << r.seed << "," << r.N << ","
       << num(r.wStar) << "," << num(r.lb1) << ","
       << r.rounds << "," << r.dEvals << "," << num(r.gap) << ","
       << (r.proven ? "Y" : "N") << "," << r.stopReason << ","
       << num(r.proxySec) << ","
       << (r.lambdaProxy >= 0 ? std::to_string(r.lambdaProxy) : "") << ","
       << (r.lambdaLbProxy >= 0 ? std::to_string(r.lambdaLbProxy) : "") << ","
       << (r.baseOk ? num(r.wBase) : "") << "," << num(r.baseSec) << ","
       << (r.lambdaBase >= 0 ? std::to_string(r.lambdaBase) : "") << ","
       << (r.lambdaLbBase >= 0 ? std::to_string(r.lambdaLbBase) : "") << ","
       << r.wCmp << "," << r.lamCmp << "\n";
}

double medianOf(std::vector<double> v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

void summarizeCase(
        std::ostream& os,
        const BenchCase& bc,
        const std::vector<SeedRow>& rows) {
    int provenY = 0, limitHit = 0;
    int better = 0, equal = 0, worse = 0, baseFail = 0;
    int lamBetter = 0, lamEqual = 0, lamWorse = 0;
    int lamTightP = 0, lamTightB = 0, lamOkP = 0, lamOkB = 0;
    int negGap = 0, provenWorse = 0;
    std::vector<double> rounds, evals, proxyT, baseT, lamP, lamB;

    os << "\n=== Case (" << bc.id << ") " << bc.title
       << "  N=" << bc.N << " Q=" << bc.D.demands.size()
       << "  seeds=" << rows.size() << " ===\n";

    for (const SeedRow& r : rows) {
        if (r.proven) ++provenY;
        if (r.stopReason == "round-limit" || r.stopReason == "wall-time-limit")
            ++limitHit;
        if (!r.baseOk) {
            ++baseFail;
            os << "  [baseline-fail] seed " << r.seed << "\n";
        } else if (r.wCmp == "better") ++better;
        else if (r.wCmp == "equal") ++equal;
        else if (r.wCmp == "worse") ++worse;

        if (r.lamCmp == "better") ++lamBetter;
        else if (r.lamCmp == "equal") ++lamEqual;
        else if (r.lamCmp == "worse") ++lamWorse;

        if (r.wlbProxyOk) {
            ++lamOkP;
            lamP.push_back((double)r.lambdaProxy);
            if (r.lambdaProxy == r.lambdaLbProxy) ++lamTightP;
        }
        if (r.wlbBaseOk) {
            ++lamOkB;
            lamB.push_back((double)r.lambdaBase);
            if (r.lambdaBase == r.lambdaLbBase) ++lamTightB;
        }
        rounds.push_back((double)r.rounds);
        evals.push_back((double)r.dEvals);
        proxyT.push_back(r.proxySec);
        if (r.baseOk) baseT.push_back(r.baseSec);
        if (r.negativeGap) {
            ++negGap;
            os << "  [NEGATIVE GAP] seed " << r.seed << " gap=" << r.gap << "\n";
        }
        if (r.provenWorseThanBase) {
            ++provenWorse;
            os << "  [PROVEN W*>W_base] seed " << r.seed
               << " W*=" << r.wStar << " W_base=" << r.wBase << "\n";
        }
    }

    auto mm = [](const std::vector<double>& v) {
        if (v.empty()) return std::string("n/a");
        const double mn = *std::min_element(v.begin(), v.end());
        const double mx = *std::max_element(v.begin(), v.end());
        const double md = medianOf(v);
        std::ostringstream s;
        s << std::fixed << std::setprecision(3)
          << "min=" << mn << " med=" << md << " max=" << mx;
        return s.str();
    };

    os << "proven=Y: " << provenY << "/" << rows.size()
       << "  limit-hit: " << limitHit << "/" << rows.size() << "\n";
    os << "W* vs W_base: better=" << better << " equal=" << equal
       << " worse=" << worse << " base-fail=" << baseFail << "\n";
    os << "λ proxy vs λ base: fewer=" << lamBetter << " equal=" << lamEqual
       << " more=" << lamWorse << "\n";
    os << "rounds: " << mm(rounds) << "\n";
    os << "D-evals: " << mm(evals) << "\n";
    os << "proxy time (s): " << mm(proxyT) << "\n";
    os << "W_base time (s): " << mm(baseT) << "\n";
    os << "λ_needed proxy: " << mm(lamP)
       << "  tight(λ=LB): " << lamTightP << "/" << lamOkP << "\n";
    os << "λ_needed base:  " << mm(lamB)
       << "  tight(λ=LB): " << lamTightB << "/" << lamOkB << "\n";
    os << "negative gaps: " << negGap
       << "  proven W*>W_base: " << provenWorse << "\n";
}

std::set<int> loadDoneSeeds(const std::string& csvPath) {
    std::set<int> done;
    std::ifstream in(csvPath);
    if (!in) return done;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        auto cols = splitCsv(line);
        if (!cols.empty()) done.insert(std::stoi(cols[0]));
    }
    return done;
}

SeedRow parseCsvRow(const std::string& line) {
    SeedRow r;
    auto c = splitCsv(line);
    if (c.size() < 18) return r;
    auto dbl = [](const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        return std::stod(s);
    };
    auto ior = [](const std::string& s) {
        if (s.empty()) return -1;
        return std::stoi(s);
    };
    r.seed = std::stoi(c[0]);
    r.N = std::stoi(c[1]);
    r.wStar = dbl(c[2]);
    r.lb1 = dbl(c[3]);
    r.rounds = std::stoi(c[4]);
    r.dEvals = std::stoi(c[5]);
    r.gap = dbl(c[6]);
    r.proven = (c[7] == "Y");
    r.stopReason = c[8];
    r.proxySec = dbl(c[9]);
    r.lambdaProxy = ior(c[10]);
    r.lambdaLbProxy = ior(c[11]);
    r.wlbProxyOk = r.lambdaProxy >= 0;
    r.wBase = dbl(c[12]);
    r.baseSec = dbl(c[13]);
    r.baseOk = !c[12].empty() && std::isfinite(r.wBase);
    r.lambdaBase = ior(c[14]);
    r.lambdaLbBase = ior(c[15]);
    r.wlbBaseOk = r.lambdaBase >= 0;
    r.wCmp = c[16];
    r.lamCmp = c[17];
    if (std::isfinite(r.gap) && r.gap < -1e-6) r.negativeGap = true;
    if (r.proven && r.baseOk && r.wStar > r.wBase + 1e-6) r.provenWorseThanBase = true;
    return r;
}

std::vector<SeedRow> loadExistingRows(const std::string& csvPath) {
    std::vector<SeedRow> rows;
    std::ifstream in(csvPath);
    if (!in) return rows;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        rows.push_back(parseCsvRow(line));
    }
    return rows;
}

}  // namespace

int main(int argc, char* argv[]) {
    // Line-buffer stdout so long overnight logs flush per seed.
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);

    const std::string root =
        ".";

    const bool estimateOnly =
        (argc >= 2 && std::string(argv[1]) == "estimate");

    std::vector<BenchCase> cases(4);
    if (!loadCase(root, 'a', "a_mpeg4", "MPEG4 decoder", cases[0])) return 1;
    if (!loadCase(root, 'b', "b_vopd", "vopd", cases[1])) return 1;
    if (!loadCase(root, 'c', "c_mwd", "mwd", cases[2])) return 1;
    if (!loadCase(root, 'd', "d_dsp", "DSP", cases[3])) return 1;

    std::cout << "Loaded cases:\n";
    for (const auto& c : cases) {
        std::cout << "  (" << c.id << ") " << c.title
                  << " N=" << c.N << " Q=" << c.D.demands.size() << "\n";
    }
    std::cout << "Seeds: 1..20 (placement generateNodes; fixed demands)\n";
    std::cout << "Total: 4 cases × 20 seeds × 2 methods = 160 invocations\n\n";

    if (argc >= 4 && std::string(argv[1]) == "diag") {
        const char want = argv[2][0];
        const int diagSeed = std::stoi(argv[3]);
        for (const BenchCase& bc : cases) {
            if (bc.id != want) continue;
            std::cout << "=== DIAG (" << bc.id << ") " << bc.title
                      << " seed=" << diagSeed << " ===\n";
            std::cout << "LB = proxy-master bound (W_proxy); W* = incumbent WC after Method D\n\n";
            SeedRow row = runSeed(bc, diagSeed, root);
            std::cout << "\n--- summary ---\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "LB_1 (round-1 proxy LB): " << row.lb1 << " mm\n";
            std::cout << "W* (final incumbent):    " << row.wStar << " mm\n";
            std::cout << "gap = W* - LB_final:     " << row.gap << " mm\n";
            std::cout << "proven: " << (row.proven ? "Y" : "N")
                      << "  stop: " << row.stopReason << "\n";
            return 0;
        }
        std::cerr << "unknown case " << argv[2] << "\n";
        return 1;
    }

    int seedLo = 1, seedHi = 20;
    std::vector<BenchCase> runCases = cases;
    if (estimateOnly) {
        seedLo = 1;
        seedHi = 1;
        std::cout << "=== ESTIMATE MODE: seed 1 only, all 4 cases ===\n\n";
    } else if (argc >= 4) {
        // ./benchmark_batch a 1 3
        const char want = argv[1][0];
        seedLo = std::stoi(argv[2]);
        seedHi = std::stoi(argv[3]);
        runCases.clear();
        for (const auto& c : cases)
            if (c.id == want) runCases.push_back(c);
        if (runCases.empty()) {
            std::cerr << "unknown case\n";
            return 1;
        }
    }

    std::vector<SeedRow> allRows;
    std::map<char, std::vector<SeedRow>> byCase;

    for (const BenchCase& bc : runCases) {
        std::cout << "---- Case (" << bc.id << ") " << bc.title << " ----\n";
        const std::string outPath =
            root + "/benchmarks/results/results_" + bc.slug + ".csv";

        std::vector<SeedRow> rows;
        std::set<int> done;
        if (!estimateOnly) {
            rows = loadExistingRows(outPath);
            done = loadDoneSeeds(outPath);
            if (!done.empty()) {
                std::cout << "  resume: " << done.size()
                          << " seeds already in " << outPath << "\n";
            }
            if (rows.empty()) {
                std::ofstream ofs(outPath);
                writeCsvHeader(ofs);
                ofs.flush();
            }
        }

        for (int seed = seedLo; seed <= seedHi; ++seed) {
            if (!estimateOnly && done.count(seed)) {
                std::cout << "  seed " << seed << " ... SKIP (resume)\n";
                continue;
            }
            std::cout << "  seed " << seed << " ..." << std::flush;
            SeedRow row;
            try {
                row = runSeed(bc, seed, root);
            } catch (const std::exception& e) {
                std::cerr << " EXCEPTION: " << e.what() << " — continuing\n";
                row.seed = seed;
                row.N = bc.N;
                row.stopReason = std::string("exception:") + e.what();
            } catch (...) {
                std::cerr << " UNKNOWN EXCEPTION — continuing\n";
                row.seed = seed;
                row.N = bc.N;
                row.stopReason = "exception:unknown";
            }
            std::cout << " W*=" << row.wStar
                      << " W_base=" << (row.baseOk ? std::to_string(row.wBase) : "FAIL")
                      << " proxy=" << std::fixed << std::setprecision(1) << row.proxySec
                      << "s base=" << row.baseSec << "s"
                      << " proven=" << (row.proven ? "Y" : "N")
                      << " λp=" << row.lambdaProxy << " λb=" << row.lambdaBase
                      << "\n";
            if (row.negativeGap || row.provenWorseThanBase) {
                std::cerr << "  !! anomaly on (" << bc.id << ", seed " << seed
                          << ") — continuing\n";
            }
            rows.push_back(row);
            allRows.push_back(row);
            if (!estimateOnly) {
                std::ofstream ofs(outPath, std::ios::app);
                writeCsvRow(ofs, row);
                ofs.flush();
            }
        }
        byCase[bc.id] = rows;
        if (!estimateOnly)
            std::cout << "Wrote/updated " << outPath << "\n";
    }

    // Summary
    std::ostringstream summary;
    for (const BenchCase& bc : runCases) {
        if (!byCase.count(bc.id)) continue;
        summarizeCase(summary, bc, byCase[bc.id]);
    }

    // Overall
    if (!estimateOnly && runCases.size() > 1) {
        summary << "\n=== OVERALL (all cases combined) ===\n";
        BenchCase all;
        all.id = '*';
        all.title = "all (a)–(d)";
        all.N = -1;
        summarizeCase(summary, all, allRows);
    }

    if (estimateOnly) {
        double totalProxy = 0, totalBase = 0;
        for (const SeedRow& r : allRows) {
            totalProxy += r.proxySec;
            totalBase += r.baseSec;
        }
        double est = 0;
        for (const SeedRow& r : allRows)
            est += 20.0 * (r.proxySec + r.baseSec);
        summary << "\n=== TIME ESTIMATE (from seed 1 × 20) ===\n";
        summary << std::fixed << std::setprecision(1);
        summary << "Measured seed-1 wall (proxy+base) sum over 4 cases: "
                << (totalProxy + totalBase) << " s\n";
        summary << "Rough full-batch estimate (×20 seeds): " << est << " s ("
                << (est / 3600.0) << " h)\n";
        summary << "Note: multi-seed variance can be large; treat as order-of-magnitude.\n";
        std::cout << summary.str();
        std::ofstream estf(root + "/benchmarks/results/estimate_seed1.txt");
        estf << summary.str();
        return 0;
    }

    std::cout << summary.str();
    {
        std::ofstream sf(root + "/benchmarks/results/summary_all.txt");
        sf << summary.str();
        std::cout << "Wrote benchmarks/results/summary_all.txt\n";
    }
    return 0;
}
