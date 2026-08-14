/*
#include "MILPSolver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

// ── Physikalische Konstanten (Ortín-Obón 2017) ───────────────
const double ALPHA = 0.0274;  // dB/mm Propagation Loss
const double BETA  = 0.010;   // dB Through Loss + Bend Loss pro Kante
const double EPS   = 1e-6;    // Tie-Breaker für nicht-kritische Demands

// ── Subtour Callback (DFJ Lazy Constraints) ──────────────────
// Wird nur für große N verwendet — für kleine N übernimmt MTZ
class SubtourCallback : public GRBCallback {
private:
    int N;
    std::vector<GRBVar>& b;

    int edge_idx(int i, int j) const { return i * N + j; }

public:
    SubtourCallback(int N, std::vector<GRBVar>& b) : N(N), b(b) {}

protected:
    void callback() override {
        try {
            if (where != GRB_CB_MIPSOL) return;

            std::vector<int> next(N, -1);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    if (getSolution(b[edge_idx(i,j)]) > 0.5) {
                        next[i] = j;
                        break;
                    }
                }

            std::vector<bool> visited(N, false);
            for (int start = 0; start < N; start++) {
                if (visited[start]) continue;
                std::vector<int> tour;
                int cur = start;
                while (cur != -1 && !visited[cur]) {
                    visited[cur] = true;
                    tour.push_back(cur);
                    cur = next[cur];
                }
                if (tour.size() > 0 && (int)tour.size() < N) {
                    GRBLinExpr expr = 0;
                    for (int k = 0; k < (int)tour.size(); k++) {
                        int from = tour[k];
                        int to   = tour[(k + 1) % tour.size()];
                        expr += b[edge_idx(from, to)];
                    }
                    addLazy(expr <= (int)tour.size() - 1);
                }
            }
        } catch (GRBException& e) {
            std::cerr << "Callback Fehler: " << e.getMessage() << "\n";
        }
    }
};

// ─────────────────────────────────────────────────────────────

MILPSolver::MILPSolver(const std::vector<Node>& nodes, const DemandMatrix& D)
    : nodes(nodes), D(D), N(nodes.size()) {}

double MILPSolver::manhattanDist(int i, int j) const {
    return std::abs(nodes[i].x - nodes[j].x)
         + std::abs(nodes[i].y - nodes[j].y);
}

double MILPSolver::insertionLoss(int i, int j) const {
    return ALPHA * manhattanDist(i, j) + BETA;
}

std::vector<RingSegment> MILPSolver::getSegments(int i, int j, int opt) const {
    double xi = nodes[i].x, yi = nodes[i].y;
    double xj = nodes[j].x, yj = nodes[j].y;
    std::vector<RingSegment> segs;
    if (opt == 0) {
        segs.push_back({true,  yi, std::min(xi,xj), std::max(xi,xj)});
        segs.push_back({false, xj, std::min(yi,yj), std::max(yi,yj)});
    } else {
        segs.push_back({false, xi, std::min(yi,yj), std::max(yi,yj)});
        segs.push_back({true,  yj, std::min(xi,xj), std::max(xi,xj)});
    }
    return segs;
}

bool MILPSolver::segmentsCross(const RingSegment& s1, const RingSegment& s2) const {
    if (s1.isHorizontal != s2.isHorizontal) {
        const RingSegment& h = s1.isHorizontal ? s1 : s2;
        const RingSegment& v = s1.isHorizontal ? s2 : s1;
        return (h.lo < v.fixed && v.fixed < h.hi)
            && (v.lo < h.fixed && h.fixed < v.hi);
    } else {
        if (std::abs(s1.fixed - s2.fixed) > 1e-9) return false;
        return s1.lo < s2.hi && s2.lo < s1.hi;
    }
}

bool MILPSolver::doesCross(int i, int j, int opt1, int k, int l, int opt2) const {
    auto segs1 = getSegments(i, j, opt1);
    auto segs2 = getSegments(k, l, opt2);
    for (auto& s1 : segs1)
        for (auto& s2 : segs2)
            if (segmentsCross(s1, s2)) return true;
    return false;
}

void MILPSolver::solve() {
    try {
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_OutputFlag, 1);
        env.start();
        GRBModel model = GRBModel(env);

        // Schwellwert: MTZ für kleine N, Lazy DFJ für große N
        const int MTZ_THRESHOLD = 20;
        const bool useLazy = (N > MTZ_THRESHOLD);

        if (useLazy) {
            model.set(GRB_IntParam_LazyConstraints, 1);
            std::cout << "Subtour-Elimination: Lazy DFJ (N=" << N << ")\n";
        } else {
            std::cout << "Subtour-Elimination: MTZ (N=" << N << ")\n";
        }

        // ── Variable W: Worst-Case IL (Hauptziel) ────────────
        GRBVar W = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "W");

        // ── Variablen b[i→j]: Kante im Ring? ────────────────
        std::vector<GRBVar> b(N * N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                b[edge_idx(i,j)] = model.addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "b_" + std::to_string(i) + "_" + std::to_string(j)
                );
            }

        // ── Variablen r[i→j]: Routing Option (0=HV, 1=VH) ──
        std::vector<GRBVar> r(N * N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                r[edge_idx(i,j)] = model.addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "r_" + std::to_string(i) + "_" + std::to_string(j)
                );
            }

        // ── Variablen f[q][i→j]: Flow für jeden Demand ───────
        int numDemands = D.demands.size();
        std::vector<std::vector<GRBVar>> f(numDemands,
            std::vector<GRBVar>(N * N));

        for (int q = 0; q < numDemands; q++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    double w = insertionLoss(i, j);
                    f[q][edge_idx(i,j)] = model.addVar(
                        0.0, 1.0, MILP_EPS * w, GRB_CONTINUOUS,
                        "f_" + std::to_string(q) + "_"
                        + std::to_string(i) + "_" + std::to_string(j)
                    );
                }
        model.update();

        // ── C1: Jeder Node hat genau einen Nachfolger ────────
        for (int i = 0; i < N; i++) {
            GRBLinExpr out = 0;
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                out += b[edge_idx(i,j)];
            }
            model.addConstr(out == 1, "C1_" + std::to_string(i));
        }

        // ── C2: Jeder Node hat genau einen Vorgänger ─────────
        for (int j = 0; j < N; j++) {
            GRBLinExpr in = 0;
            for (int i = 0; i < N; i++) {
                if (i == j) continue;
                in += b[edge_idx(i,j)];
            }
            model.addConstr(in == 1, "C2_" + std::to_string(j));
        }

        // ── C3: Keine Gegenrichtung gleichzeitig ─────────────
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++)
                model.addConstr(
                    b[edge_idx(i,j)] + b[edge_idx(j,i)] <= 1,
                    "C3_" + std::to_string(i) + "_" + std::to_string(j)
                );

        // ── C4: Subtour-Elimination ───────────────────────────
        // MTZ für N ≤ 20: starke LP-Relaxation, schnell für kleine N
        // Lazy DFJ für N > 20: skaliert besser für große N
        if (!useLazy) {
            std::vector<GRBVar> u(N);
            for (int i = 1; i < N; i++)
                u[i] = model.addVar(
                    1.0, N - 1, 0.0, GRB_CONTINUOUS,
                    "u_" + std::to_string(i)
                );
            model.update();

            for (int i = 1; i < N; i++)
                for (int j = 1; j < N; j++) {
                    if (i == j) continue;
                    model.addConstr(
                        u[i] - u[j] + N * b[edge_idx(i,j)] <= N - 1,
                        "C4_" + std::to_string(i) + "_" + std::to_string(j)
                    );
                }
        }
        // Lazy DFJ wird über den Callback zur Laufzeit eingeworfen:
        // Σ b[i→j] für alle Kanten im Teilring S ≤ |S| - 1

        // ── C5: Crossing-free Constraint ─────────────────────
        int crossCount = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                for (int k = 0; k < N; k++)
                    for (int l = 0; l < N; l++) {
                        if (k == l) continue;
                        if (edge_idx(i,j) >= edge_idx(k,l)) continue;
                        for (int opt1 = 0; opt1 <= 1; opt1++)
                            for (int opt2 = 0; opt2 <= 1; opt2++) {
                                if (doesCross(i, j, opt1, k, l, opt2)) {
                                    GRBLinExpr lhs = b[edge_idx(i,j)]
                                                   + b[edge_idx(k,l)];
                                    if (opt1 == 0) lhs += 1 - r[edge_idx(i,j)];
                                    else           lhs += r[edge_idx(i,j)];
                                    if (opt2 == 0) lhs += 1 - r[edge_idx(k,l)];
                                    else           lhs += r[edge_idx(k,l)];
                                    model.addConstr(lhs <= 3,
                                        "C5_" + std::to_string(i) + "_"
                                        + std::to_string(j) + "_"
                                        + std::to_string(k) + "_"
                                        + std::to_string(l) + "_"
                                        + std::to_string(opt1)
                                        + std::to_string(opt2));
                                    crossCount++;
                                }
                            }
                    }
            }
        std::cout << "Crossing-Constraints: " << crossCount << "\n";

        // ── C6: Flow-Kopplung an Ring ─────────────────────────
        for (int q = 0; q < numDemands; q++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    model.addConstr(
                        f[q][edge_idx(i,j)] <= b[edge_idx(i,j)] + b[edge_idx(j,i)],
                        "C6_" + std::to_string(q) + "_"
                        + std::to_string(i) + "_" + std::to_string(j)
                    );
                }

        // ── C7: Flow-Erhaltung ────────────────────────────────
        for (int q = 0; q < numDemands; q++) {
            int s = D.demands[q].first;
            int d = D.demands[q].second;
            for (int v = 0; v < N; v++) {
                GRBLinExpr out = 0, in = 0;
                for (int j = 0; j < N; j++) {
                    if (j == v) continue;
                    out += f[q][edge_idx(v,j)];
                    in  += f[q][edge_idx(j,v)];
                }
                int rhs = 0;
                if (v == s) rhs = +1;
                if (v == d) rhs = -1;
                model.addConstr(out - in == rhs,
                    "C7_" + std::to_string(q) + "_" + std::to_string(v));
            }
        }

        // ── C8: Worst-case Bound ──────────────────────────────
        for (int q = 0; q < numDemands; q++) {
            GRBLinExpr il = 0;
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    il += f[q][edge_idx(i,j)] * insertionLoss(i,j);
                }
            model.addConstr(W >= il, "C8_" + std::to_string(q));
        }

        // ── Callback anbinden + Optimieren ───────────────────
        SubtourCallback cb(N, b);
        if (useLazy) model.setCallback(&cb);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.optimize();

        // ── Ergebnis ─────────────────────────────────────────
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            std::cout << "\n========================================\n";
            std::cout << "Optimaler Worst-Case IL (W): "
                      << W.get(GRB_DoubleAttr_X) << " dB\n";
            std::cout << "========================================\n";

            // Ring-Reihenfolge
            std::cout << "\nRing-Reihenfolge:\n  Node 1";
            int current = 0;
            for (int step = 0; step < N - 1; step++)
                for (int j = 0; j < N; j++) {
                    if (j == current) continue;
                    if (b[edge_idx(current,j)].get(GRB_DoubleAttr_X) > 0.5) {
                        std::cout << " → " << j + 1;
                        current = j;
                        break;
                    }
                }
            std::cout << " → 1\n";

            // IL pro Demand
            std::cout << "\nInsertion Loss pro Demand:\n";
            double worstIL = 0.0;
            int worstQ = -1;
            for (int q = 0; q < numDemands; q++) {
                double il = 0.0;
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        il += f[q][edge_idx(i,j)].get(GRB_DoubleAttr_X)
                            * insertionLoss(i,j);
                    }
                std::cout << "  Demand " << D.demands[q].first + 1
                          << " → " << D.demands[q].second + 1
                          << ": " << il << " dB\n";
                if (il > worstIL) { worstIL = il; worstQ = q; }
            }
            std::cout << "\nWorst-case: Demand "
                      << D.demands[worstQ].first + 1 << " → "
                      << D.demands[worstQ].second + 1
                      << " | " << worstIL << " dB\n";

            // CSV Ring
            std::ofstream csv("../ring.csv");
            csv << "from_id,to_id,from_x,from_y,to_x,to_y,routing\n";
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    if (b[edge_idx(i,j)].get(GRB_DoubleAttr_X) > 0.5) {
                        int opt = (int)(r[edge_idx(i,j)].get(GRB_DoubleAttr_X) + 0.5);
                        csv << i+1 << "," << j+1 << ","
                            << nodes[i].x << "," << nodes[i].y << ","
                            << nodes[j].x << "," << nodes[j].y << ","
                            << (opt == 0 ? "HV" : "VH") << "\n";
                    }
                }
            csv.close();

            // CSV Nodes
            std::ofstream ncsv("../nodes.csv");
            ncsv << "id,x,y\n";
            for (auto& n : nodes)
                ncsv << n.id+1 << "," << n.x << "," << n.y << "\n";
            ncsv.close();

            // CSV Demands
            std::ofstream dcsv("../demands.csv");
            dcsv << "demand,sender,receiver,from_id,to_id,direction\n";
            for (int q = 0; q < numDemands; q++) {
                int s = D.demands[q].first;
                int d = D.demands[q].second;

                bool isCW = false;
                for (int i = 0; i < N && !isCW; i++)
                    for (int j = 0; j < N && !isCW; j++) {
                        if (i == j) continue;
                        if (b[edge_idx(i,j)].get(GRB_DoubleAttr_X) > 0.5 &&
                            f[q][edge_idx(i,j)].get(GRB_DoubleAttr_X) > 0.5)
                            isCW = true;
                    }

                std::string dir = isCW ? "CW" : "CCW";
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        if (f[q][edge_idx(i,j)].get(GRB_DoubleAttr_X) > 0.5) {
                            dcsv << q << ","
                                 << s+1 << "," << d+1 << ","
                                 << i+1 << "," << j+1 << ","
                                 << dir << "\n";
                        }
                    }
            }
            dcsv.close();

            std::cout << "\nCSV exportiert: ring.csv, nodes.csv, demands.csv\n";

        } else {
            std::cout << "Keine optimale Lösung. Status: "
                      << model.get(GRB_IntAttr_Status) << "\n";
        }

    } catch (GRBException& e) {
        std::cerr << "Gurobi Fehler: " << e.getMessage() << "\n";
    }
}
*/

#include "MILPSolver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "PhysicalConstants.h"

class MilpLazyCallback : public GRBCallback {
private:
    int N;
    std::vector<GRBVar>& b;
    std::vector<GRBVar>& r;
    const MILPSolver* solver;
    bool subtour;
    bool crossing;
    int edge_idx(int i, int j) const { return i * N + j; }

public:
    MilpLazyCallback(
            int N,
            std::vector<GRBVar>& b,
            std::vector<GRBVar>& r,
            const MILPSolver* solver,
            bool subtour,
            bool crossing)
        : N(N), b(b), r(r), solver(solver), subtour(subtour), crossing(crossing) {}

protected:
    void callback() override {
        try {
            if (where != GRB_CB_MIPSOL) return;

            if (subtour) {
                std::vector<int> next(N, -1);
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        if (getSolution(b[edge_idx(i, j)]) > 0.5) {
                            next[i] = j;
                            break;
                        }
                    }

                std::vector<bool> visited(N, false);
                for (int start = 0; start < N; start++) {
                    if (visited[start]) continue;
                    std::vector<int> tour;
                    int cur = start;
                    while (cur != -1 && !visited[cur]) {
                        visited[cur] = true;
                        tour.push_back(cur);
                        cur = next[cur];
                    }
                    if (!tour.empty() && (int)tour.size() < N) {
                        GRBLinExpr expr = 0;
                        for (int k = 0; k < (int)tour.size(); k++) {
                            int from = tour[k];
                            int to = tour[(k + 1) % tour.size()];
                            expr += b[edge_idx(from, to)];
                        }
                        addLazy(expr <= (int)tour.size() - 1);
                    }
                }
            }

            if (crossing) {
                struct ActiveEdge { int i; int j; int opt; };
                std::vector<ActiveEdge> active;
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        if (i == j) continue;
                        if (getSolution(b[edge_idx(i, j)]) > 0.5) {
                            const int opt = static_cast<int>(
                                getSolution(r[edge_idx(i, j)]) + 0.5);
                            active.push_back({i, j, opt});
                        }
                    }
                }

                for (size_t a = 0; a < active.size(); ++a) {
                    for (size_t bIdx = a + 1; bIdx < active.size(); ++bIdx) {
                        const ActiveEdge& e1 = active[a];
                        const ActiveEdge& e2 = active[bIdx];
                        if (!solver->doesCross(
                                e1.i, e1.j, e1.opt, e2.i, e2.j, e2.opt))
                            continue;
                        GRBLinExpr lhs = b[edge_idx(e1.i, e1.j)] + b[edge_idx(e2.i, e2.j)];
                        if (e1.opt == 0) lhs += 1 - r[edge_idx(e1.i, e1.j)];
                        else             lhs += r[edge_idx(e1.i, e1.j)];
                        if (e2.opt == 0) lhs += 1 - r[edge_idx(e2.i, e2.j)];
                        else             lhs += r[edge_idx(e2.i, e2.j)];
                        addLazy(lhs <= 3);
                    }
                }
            }
        } catch (GRBException& e) {
            std::cerr << "Callback Fehler: " << e.getMessage() << "\n";
        }
    }
};

MILPSolver::MILPSolver(const std::vector<Node>& nodes, const DemandMatrix& D)
    : nodes(nodes), D(D), N(nodes.size()), env(true) {
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
}

double MILPSolver::computeBigM() const {
    double maxEdge = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (i != j)
                maxEdge = std::max(maxEdge, manhattanDist(i, j));
    return static_cast<double>(N) * maxEdge;
}

void MILPSolver::applyWarmStartPositionModel(
        const MILPSolveResult& warmStart,
        const std::set<int>& excludedDemands,
        std::vector<GRBVar>& pos,
        std::vector<GRBVar>& wrapVar,
        std::vector<GRBVar>& yVar,
        std::vector<GRBVar>& forwardVar,
        std::vector<GRBVar>& demandDistVar,
        std::vector<GRBVar>& dirVar) {
    if (warmStart.tour.size() < 2 || pos.empty())
        return;

    std::vector<double> posVal(N, 0.0);
    posVal[warmStart.tour.front()] = 0.0;
    double cum = 0.0;
    for (size_t k = 0; k + 1 < warmStart.tour.size(); ++k) {
        const int from = warmStart.tour[k];
        const int to = warmStart.tour[k + 1];
        cum += manhattanDist(from, to);
        posVal[to] = cum;
    }

    for (int i = 0; i < N; ++i) {
        if (i == 0) continue;
        pos[i].set(GRB_DoubleAttr_Start, posVal[i]);
    }

    const double Ltot = cum;
    const int numDemands = (int)D.demands.size();

    MILPSolveResult routed = warmStart;
    applyShortestArcDemandRouting(routed, excludedDemands);

    for (int q = 0; q < numDemands; ++q) {
        if (excludedDemands.count(q)) continue;
        const int s = D.demands[q].first;
        const int d = D.demands[q].second;

        double yStart = 0.0;
        double wrapStart = 0.0;
        if (posVal[d] + MILP_EPS < posVal[s]) {
            yStart = Ltot;
            wrapStart = 1.0;
        }
        const double fwd = posVal[d] - posVal[s] + yStart;
        const double bwd = Ltot - fwd;
        const double distStart = std::min(fwd, bwd);
        const bool takeForward = (fwd <= bwd + MILP_EPS);

        wrapVar[q].set(GRB_DoubleAttr_Start, wrapStart);
        yVar[q].set(GRB_DoubleAttr_Start, yStart);
        forwardVar[q].set(GRB_DoubleAttr_Start, fwd);
        dirVar[q].set(GRB_DoubleAttr_Start, takeForward ? 1.0 : 0.0);
        demandDistVar[q].set(GRB_DoubleAttr_Start, distStart);

        if (q < (int)routed.demandDistance.size()
                && std::abs(routed.demandDistance[q] - distStart) > 1e-3) {
            demandDistVar[q].set(GRB_DoubleAttr_Start, routed.demandDistance[q]);
        }
    }
}

void MILPSolver::verifyPositionModelSolution(
        const MILPSolveResult& result,
        const std::vector<GRBVar>& demandDistVar,
        const std::vector<GRBVar>& forwardVar,
        const std::vector<GRBVar>& pos,
        double L_total_val,
        const std::set<int>& excludedDemands) const {
    constexpr double tol = 1e-4;
    const int numDemands = (int)D.demands.size();

    for (int q = 0; q < numDemands; ++q) {
        if (excludedDemands.count(q)) continue;
        if (q >= (int)demandDistVar.size()) continue;
        const double fwd = forwardVar[q].get(GRB_DoubleAttr_X);
        const double bwd = L_total_val - fwd;
        const double expected = std::min(fwd, bwd);
        const double actual = demandDistVar[q].get(GRB_DoubleAttr_X);
        if (std::abs(actual - expected) > tol) {
            std::cerr << "[PositionModel] WARNING: demand " << q
                      << " demandDistance=" << actual
                      << " but min(forward,backward)=" << expected << "\n";
        }
        if (q < (int)result.demandDistance.size()
                && std::abs(result.demandDistance[q] - actual) > tol) {
            std::cerr << "[PositionModel] WARNING: demand " << q
                      << " arc routing distance=" << result.demandDistance[q]
                      << " != model demandDistance=" << actual << "\n";
        }
    }

    if (result.tour.size() < 2) return;

    const int tourLen = (int)result.tour.size();
    for (int step = 0; step < tourLen - 1; ++step) {
        const int cur = result.tour[step];
        const int nxt = result.tour[step + 1];
        const double pCur = pos[cur].get(GRB_DoubleAttr_X);
        const double pNxt = pos[nxt].get(GRB_DoubleAttr_X);
        if (pNxt <= pCur + tol) {
            std::cerr << "[PositionModel] WARNING: pos not increasing along tour at "
                      << cur << "->" << nxt << " (" << pCur << " -> " << pNxt << ")\n";
        }
    }
}

double MILPSolver::manhattanDist(int i, int j) const {
    return std::abs(nodes[i].x - nodes[j].x) + std::abs(nodes[i].y - nodes[j].y);
}

double MILPSolver::insertionLoss(int i, int j) const {
    return ALPHA * manhattanDist(i, j) + BETA_THROUGH;
}

double MILPSolver::reportingDemandEdgeIL(int i, int j) const {
    return reportingRingEdgeIL_dB(
        nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y);
}

std::vector<RingSegment> MILPSolver::getSegments(int i, int j, int opt) const {
    double xi = nodes[i].x, yi = nodes[i].y;
    double xj = nodes[j].x, yj = nodes[j].y;
    std::vector<RingSegment> segs;
    if (opt == 0) {
        segs.push_back({true,  yi, std::min(xi,xj), std::max(xi,xj)});
        segs.push_back({false, xj, std::min(yi,yj), std::max(yi,yj)});
    } else {
        segs.push_back({false, xi, std::min(yi,yj), std::max(yi,yj)});
        segs.push_back({true,  yj, std::min(xi,xj), std::max(xi,xj)});
    }
    return segs;
}

bool MILPSolver::segmentsCross(const RingSegment& s1, const RingSegment& s2) const {
    if (s1.isHorizontal != s2.isHorizontal) {
        const RingSegment& h = s1.isHorizontal ? s1 : s2;
        const RingSegment& v = s1.isHorizontal ? s2 : s1;
        return (h.lo < v.fixed && v.fixed < h.hi)
            && (v.lo < h.fixed && h.fixed < v.hi);
    }
    if (std::abs(s1.fixed - s2.fixed) > 1e-9) return false;
    return s1.lo < s2.hi && s2.lo < s1.hi;
}

bool MILPSolver::segmentsTouchPerpendicular(const RingSegment& s1, const RingSegment& s2) const {
    if (s1.isHorizontal == s2.isHorizontal) return false;
    const RingSegment& h = s1.isHorizontal ? s1 : s2;
    const RingSegment& v = s1.isHorizontal ? s2 : s1;
    constexpr double eps = 1e-9;
    return (h.lo <= v.fixed + eps && v.fixed <= h.hi + eps)
        && (v.lo <= h.fixed + eps && h.fixed <= v.hi + eps);
}

bool MILPSolver::edgesShareEndpoint(int i, int j, int k, int l) const {
    return i == k || i == l || j == k || j == l;
}

bool MILPSolver::collinearOverlapBeyondShared(
        int i, int j, int opt1, int k, int l, int opt2, int sharedNode) const {
    const double px = nodes[sharedNode].x;
    const double py = nodes[sharedNode].y;
    constexpr double eps = 1e-9;

    const auto segs1 = getSegments(i, j, opt1);
    const auto segs2 = getSegments(k, l, opt2);
    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (s1.isHorizontal != s2.isHorizontal) continue;
            if (std::abs(s1.fixed - s2.fixed) > eps) continue;

            const double olLo = std::max(s1.lo, s2.lo);
            const double olHi = std::min(s1.hi, s2.hi);
            if (olHi - olLo <= eps) continue;

            const double sharedCoord = s1.isHorizontal ? px : py;
            if (olHi > sharedCoord + eps || olLo < sharedCoord - eps)
                return true;
        }
    }
    return false;
}

bool MILPSolver::doesCross(int i, int j, int opt1, int k, int l, int opt2) const {
    const auto segs1 = getSegments(i, j, opt1);
    const auto segs2 = getSegments(k, l, opt2);
    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (segmentsCross(s1, s2)) return true;
        }
    }

    if (edgesShareEndpoint(i, j, k, l)) {
        const int sharedNodes[] = {i, j, k, l};
        for (int s : sharedNodes) {
            const bool isShared = (s == i || s == j) && (s == k || s == l);
            if (!isShared) continue;
            if (collinearOverlapBeyondShared(i, j, opt1, k, l, opt2, s))
                return true;
        }
        return false;
    }

    for (const RingSegment& s1 : segs1) {
        for (const RingSegment& s2 : segs2) {
            if (segmentsTouchPerpendicular(s1, s2)) return true;
        }
    }
    return false;
}

bool MILPSolver::pathPassesThroughNode(int i, int j, int opt, int node) const {
    if (node == i || node == j) return false;

    const double xk = nodes[node].x;
    const double yk = nodes[node].y;
    const double xi = nodes[i].x, yi = nodes[i].y;
    const double xj = nodes[j].x, yj = nodes[j].y;
    constexpr double eps = 1e-9;

    auto onClosedSegment = [&](const RingSegment& s) {
        if (s.isHorizontal) {
            if (std::abs(yk - s.fixed) > eps) return false;
            return xk >= s.lo - eps && xk <= s.hi + eps;
        }
        if (std::abs(xk - s.fixed) > eps) return false;
        return yk >= s.lo - eps && yk <= s.hi + eps;
    };

    for (const RingSegment& s : getSegments(i, j, opt)) {
        if (!onClosedSegment(s)) continue;
        const bool atI = std::abs(xk - xi) < eps && std::abs(yk - yi) < eps;
        const bool atJ = std::abs(xk - xj) < eps && std::abs(yk - yj) < eps;
        if (!atI && !atJ) return true;
    }
    return false;
}

void MILPSolver::computeGeometryConstraintsOnce() {
    cachedCrossingPairs.clear();
    cachedPassThroughTriples.clear();

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++) {
                    if (k == l) continue;
                    if (edge_idx(i,j) >= edge_idx(k,l)) continue;
                    for (int opt1 = 0; opt1 <= 1; opt1++)
                        for (int opt2 = 0; opt2 <= 1; opt2++) {
                            if (doesCross(i, j, opt1, k, l, opt2))
                                cachedCrossingPairs.push_back({i, j, opt1, k, l, opt2});
                        }
                }
        }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            for (int opt = 0; opt <= 1; opt++)
                for (int node = 0; node < N; node++) {
                    if (pathPassesThroughNode(i, j, opt, node))
                        cachedPassThroughTriples.push_back({i, j, opt, node});
                }
        }

    geometryConstraintsComputed = true;
    std::cout << "Crossing-Constraints: " << cachedCrossingPairs.size()
              << ", Pass-through-Constraints: " << cachedPassThroughTriples.size()
              << std::endl;
}

void MILPSolver::applyShortestArcDemandRouting(
        MILPSolveResult& result,
        const std::set<int>& excludedDemands) const {
    const int tourLen = (int)result.tour.size();
    if (tourLen == 0) return;

    std::vector<int> tourPos(N, -1);
    for (int i = 0; i < tourLen; ++i)
        tourPos[result.tour[i]] = i;

    auto collectArc = [&](int ps, int pd, int step) {
        std::vector<std::pair<int, int>> edges;
        int k = ps;
        while (k != pd) {
            const int nk = (k + step + tourLen) % tourLen;
            edges.emplace_back(result.tour[k], result.tour[nk]);
            k = nk;
        }
        return edges;
    };

    auto arcMetrics = [&](const std::vector<std::pair<int, int>>& edges) {
        double dist = 0.0;
        double il = 0.0;
        int bends = 0;
        for (const auto& [i, j] : edges) {
            dist += manhattanDist(i, j);
            il += reportingDemandEdgeIL(i, j);
            if (ringEdgeNeedsBend(nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y))
                ++bends;
        }
        return std::tuple<double, double, int>{dist, il, bends};
    };

    const int numDemands = (int)D.demands.size();
    result.demandFlowEdges.assign(numDemands, {});
    result.demandIL.assign(numDemands, 0.0);
    result.demandDistance.assign(numDemands, 0.0);
    result.demandBendCount.assign(numDemands, 0);

    double maxW = 0.0;
    for (int q = 0; q < numDemands; ++q) {
        if (excludedDemands.count(q)) continue;

        const int s = D.demands[q].first;
        const int d = D.demands[q].second;
        if (tourPos[s] < 0 || tourPos[d] < 0) continue;

        const int ps = tourPos[s];
        const int pd = tourPos[d];
        const auto cw = collectArc(ps, pd, +1);
        const auto ccw = collectArc(ps, pd, -1);
        const auto cwM = arcMetrics(cw);
        const auto ccwM = arcMetrics(ccw);

        const double ccwDist = std::get<0>(ccwM);
        const double cwDist = std::get<0>(cwM);
        const bool takeCcw =
            (ccwDist < cwDist - MILP_EPS) ||
            (std::abs(ccwDist - cwDist) <= MILP_EPS &&
             std::get<2>(ccwM) < std::get<2>(cwM));
        const auto& chosen = takeCcw ? ccw : cw;
        const auto& chosenM = takeCcw ? ccwM : cwM;

        result.demandFlowEdges[q] = chosen;
        result.demandDistance[q] = std::get<0>(chosenM);
        result.demandIL[q] = std::get<1>(chosenM);
        result.demandBendCount[q] = std::get<2>(chosenM);
        maxW = std::max(maxW, result.demandDistance[q]);
    }
    result.W = maxW;
}

void MILPSolver::selfCheckRingCrossings(
        const std::vector<GRBVar>& b,
        const std::vector<GRBVar>& r,
        const std::string& context) const {
    struct ActiveEdge {
        int i;
        int j;
        int opt;
    };

    std::vector<ActiveEdge> active;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            if (b[edge_idx(i, j)].get(GRB_DoubleAttr_X) > 0.5) {
                const int opt = static_cast<int>(r[edge_idx(i, j)].get(GRB_DoubleAttr_X) + 0.5);
                active.push_back({i, j, opt});
            }
        }
    }

    bool found = false;
    for (size_t a = 0; a < active.size(); ++a) {
        for (size_t bIdx = a + 1; bIdx < active.size(); ++bIdx) {
            const ActiveEdge& e1 = active[a];
            const ActiveEdge& e2 = active[bIdx];
            if (!doesCross(e1.i, e1.j, e1.opt, e2.i, e2.j, e2.opt)) continue;
            std::cerr << "[SelfCheck] WARNING: solved model contains illegal ring edge geometry! ("
                      << e1.i << "," << e1.j << ") opt=" << e1.opt
                      << " x (" << e2.i << "," << e2.j << ") opt=" << e2.opt
                      << " [" << context << "]\n";
            found = true;
        }
    }

    if (!found) {
        std::cout << "[SelfCheck] No crossing ring edges in solved model [" << context << "]\n";
    }
}

// ── solve with optional demand exclusion and BestObjStop ─────────────────
MILPSolveResult MILPSolver::solve(
        bool optimizeWorstCase,
        const std::string& fileSuffix,
        const std::set<int>& excludedDemands,
        double bestObjStop,
        const MILPSolveResult* warmStart,
        double timeLimitSeconds,
        bool quiet,
        const MILPSolveOptions& options) {
    MILPSolveResult result;
    try {
        env.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);
        GRBModel model = GRBModel(env);

        const int MTZ_THRESHOLD = 20;
        const bool useLazySubtour = (N > MTZ_THRESHOLD);
        const bool useLazyC5 = options.lazyC5;
        if (useLazySubtour || useLazyC5)
            model.set(GRB_IntParam_LazyConstraints, 1);

        const double BIGM = computeBigM();
        const int numDemands = (int)D.demands.size();

        double w_cost = optimizeWorstCase ? 1.0 : 0.0;
        GRBVar W = model.addVar(0.0, GRB_INFINITY, w_cost, GRB_CONTINUOUS, "W");

        std::vector<GRBVar> b(N * N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                double b_cost = optimizeWorstCase ? 0.0 : insertionLoss(i, j);
                b[edge_idx(i,j)] = model.addVar(
                    0.0, 1.0, b_cost, GRB_BINARY,
                    "b_" + std::to_string(i) + "_" + std::to_string(j)
                );
            }

        std::vector<GRBVar> r(N * N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                r[edge_idx(i,j)] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_" + std::to_string(i) + "_" + std::to_string(j));
            }

        if (warmStart != nullptr && warmStart->success) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    b[edge_idx(i,j)].set(GRB_DoubleAttr_Start, 0.0);
                    r[edge_idx(i,j)].set(GRB_DoubleAttr_Start, 0.0);
                }
            for (const auto& te : warmStart->tourEdges) {
                b[edge_idx(te.from, te.to)].set(GRB_DoubleAttr_Start, 1.0);
                r[edge_idx(te.from, te.to)].set(GRB_DoubleAttr_Start, static_cast<double>(te.routing));
            }
        }

        std::vector<std::vector<GRBVar>> f;
        std::vector<GRBVar> pos;
        std::vector<GRBVar> demandDistVar(numDemands);
        std::vector<GRBVar> forwardVar(numDemands);
        std::vector<GRBVar> wrapVar(numDemands);
        std::vector<GRBVar> yVar(numDemands);
        std::vector<GRBVar> dirVar(numDemands);

        if (!optimizeWorstCase) {
            f.assign(numDemands, std::vector<GRBVar>(N * N));
            for (int q = 0; q < numDemands; q++) {
                if (excludedDemands.count(q)) continue;
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        f[q][edge_idx(i,j)] = model.addVar(
                            0.0, 1.0, 0.0, GRB_CONTINUOUS,
                            "f_" + std::to_string(q) + "_" + std::to_string(i) + "_" + std::to_string(j));
                    }
            }
        } else {
            pos.resize(N);
            for (int i = 0; i < N; ++i) {
                if (i == 0)
                    pos[i] = model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, "pos_0");
                else
                    pos[i] = model.addVar(0.0, BIGM, 0.0, GRB_CONTINUOUS, "pos_" + std::to_string(i));
            }
        }
        model.update();

        for (int i = 0; i < N; i++) {
            GRBLinExpr out = 0;
            for (int j = 0; j < N; j++) if (i != j) out += b[edge_idx(i,j)];
            model.addConstr(out == 1, "C1_" + std::to_string(i));
        }

        for (int j = 0; j < N; j++) {
            GRBLinExpr in = 0;
            for (int i = 0; i < N; i++) if (i != j) in += b[edge_idx(i,j)];
            model.addConstr(in == 1, "C2_" + std::to_string(j));
        }

        // C3 omitted: redundant with C1+C2+C4 for |S|=2 (length-2 subtours already excluded).

        if (!useLazySubtour) {
            std::vector<GRBVar> u(N);
            for (int i = 1; i < N; i++)
                u[i] = model.addVar(1.0, N - 1, 0.0, GRB_CONTINUOUS, "u_" + std::to_string(i));
            model.update();
            for (int i = 1; i < N; i++)
                for (int j = 1; j < N; j++) {
                    if (i == j) continue;
                    model.addConstr(u[i] - u[j] + N * b[edge_idx(i,j)] <= N - 1, "C4_" + std::to_string(i) + "_" + std::to_string(j));
                }
            if (!options.skipSymmetryBreak) {
                // Break rotation × direction symmetry (may exclude tours without edge 0→1).
                model.addConstr(b[edge_idx(0,1)] == 1, "SymmetryBreak_0_1");
            }
        }

        if (!geometryConstraintsComputed)
            computeGeometryConstraintsOnce();

        if (!useLazyC5) {
            for (const CrossingPair& cp : cachedCrossingPairs) {
                GRBLinExpr lhs = b[edge_idx(cp.i, cp.j)] + b[edge_idx(cp.k, cp.l)];
                if (cp.opt1 == 0) lhs += 1 - r[edge_idx(cp.i, cp.j)];
                else              lhs += r[edge_idx(cp.i, cp.j)];
                if (cp.opt2 == 0) lhs += 1 - r[edge_idx(cp.k, cp.l)];
                else              lhs += r[edge_idx(cp.k, cp.l)];
                model.addConstr(lhs <= 3, "C5_...");
            }
        } else if (!quiet) {
            std::cout << "C5: lazy mode (" << cachedCrossingPairs.size()
                      << " candidate pairs, cuts on demand)\n";
        }

        for (const PassThroughTriple& pt : cachedPassThroughTriples) {
            GRBLinExpr lhs = b[edge_idx(pt.i, pt.j)];
            if (pt.opt == 0) lhs += 1 - r[edge_idx(pt.i, pt.j)];
            else             lhs += r[edge_idx(pt.i, pt.j)];
            model.addConstr(lhs <= 1, "C5b_...");
        }

        GRBLinExpr L_total = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                L_total += b[edge_idx(i,j)] * manhattanDist(i, j);
            }

        if (optimizeWorstCase) {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j || j == 0) continue;
                    const double w_ij = manhattanDist(i, j);
                    model.addConstr(
                        pos[j] - pos[i] - w_ij <= BIGM * (1 - b[edge_idx(i,j)]),
                        "PosLink_ub_" + std::to_string(i) + "_" + std::to_string(j));
                    model.addConstr(
                        pos[j] - pos[i] - w_ij >= -BIGM * (1 - b[edge_idx(i,j)]),
                        "PosLink_lb_" + std::to_string(i) + "_" + std::to_string(j));
                }

            for (int q = 0; q < numDemands; ++q) {
                if (excludedDemands.count(q)) continue;
                const int s = D.demands[q].first;
                const int d = D.demands[q].second;

                wrapVar[q] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "wrap_" + std::to_string(q));
                yVar[q] = model.addVar(0.0, BIGM, 0.0, GRB_CONTINUOUS, "y_" + std::to_string(q));
                forwardVar[q] = model.addVar(0.0, BIGM, 0.0, GRB_CONTINUOUS, "forward_" + std::to_string(q));
                demandDistVar[q] = model.addVar(0.0, BIGM, 0.0, GRB_CONTINUOUS, "demandDist_" + std::to_string(q));
                dirVar[q] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "dir_" + std::to_string(q));

                model.addConstr(yVar[q] <= L_total, "WrapY_L_" + std::to_string(q));
                model.addConstr(yVar[q] <= BIGM * wrapVar[q], "WrapY_U1_" + std::to_string(q));
                model.addConstr(yVar[q] >= L_total - BIGM * (1 - wrapVar[q]), "WrapY_U2_" + std::to_string(q));
                model.addConstr(yVar[q] >= 0.0, "WrapY_LB_" + std::to_string(q));

                model.addConstr(forwardVar[q] == pos[d] - pos[s] + yVar[q], "Forward_" + std::to_string(q));
                model.addConstr(forwardVar[q] >= 0.0, "Forward_LB_" + std::to_string(q));
                model.addConstr(forwardVar[q] <= L_total, "Forward_UB_" + std::to_string(q));

                model.addConstr(demandDistVar[q] <= forwardVar[q], "DemDist_F_" + std::to_string(q));
                model.addConstr(demandDistVar[q] <= L_total - forwardVar[q], "DemDist_B_" + std::to_string(q));
                model.addConstr(
                    demandDistVar[q] >= forwardVar[q] - BIGM * (1 - dirVar[q]),
                    "DemDist_selF_" + std::to_string(q));
                model.addConstr(
                    demandDistVar[q] >= (L_total - forwardVar[q]) - BIGM * dirVar[q],
                    "DemDist_selB_" + std::to_string(q));

                model.addConstr(W >= demandDistVar[q], "C8_" + std::to_string(q));
            }
            if (warmStart != nullptr && warmStart->success) {
                applyWarmStartPositionModel(
                    *warmStart, excludedDemands, pos, wrapVar, yVar,
                    forwardVar, demandDistVar, dirVar);
            }
            model.update();
        } else {
            for (int q = 0; q < numDemands; q++) {
                if (excludedDemands.count(q)) continue;
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        model.addConstr(
                            f[q][edge_idx(i,j)] <= b[edge_idx(i,j)] + b[edge_idx(j,i)],
                            "C6_" + std::to_string(q) + "_" + std::to_string(i) + "_" + std::to_string(j));
                    }
            }

            for (int q = 0; q < numDemands; q++) {
                if (excludedDemands.count(q)) continue;
                int s = D.demands[q].first;
                int d = D.demands[q].second;
                for (int v = 0; v < N; v++) {
                    GRBLinExpr out = 0, in = 0;
                    for (int j = 0; j < N; j++) {
                        if (j == v) continue;
                        out += f[q][edge_idx(v,j)];
                        in  += f[q][edge_idx(j,v)];
                    }
                    int rhs = 0;
                    if (v == s) rhs = +1;
                    if (v == d) rhs = -1;
                    model.addConstr(out - in == rhs,
                        "C7_" + std::to_string(q) + "_" + std::to_string(v));
                }
            }

            for (int q = 0; q < numDemands; q++) {
                if (excludedDemands.count(q)) continue;
                GRBLinExpr pathDist = 0;
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        if (i == j) continue;
                        pathDist += f[q][edge_idx(i,j)] * manhattanDist(i,j);
                    }
                model.addConstr(W >= pathDist, "C8_" + std::to_string(q));
            }
        }

        MilpLazyCallback cb(N, b, r, this, useLazySubtour, useLazyC5);
        if (useLazySubtour || useLazyC5)
            model.setCallback(&cb);

        const bool warmStarted =
            warmStart != nullptr && warmStart->success;
        const bool isReSolve = bestObjStop >= 0.0;

        if (options.mipFocus >= 0)
            model.set(GRB_IntParam_MIPFocus, options.mipFocus);
        else if (optimizeWorstCase && warmStarted)
            model.set(GRB_IntParam_MIPFocus, 1);

        if (options.improveStartTime >= 0.0)
            model.set(GRB_DoubleParam_ImproveStartTime, options.improveStartTime);
        else if (optimizeWorstCase && isReSolve)
            model.set(GRB_DoubleParam_ImproveStartTime, 10.0);

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        if (bestObjStop >= 0.0)
            model.set(GRB_DoubleParam_BestObjStop, bestObjStop - 1e-4);
        if (timeLimitSeconds >= 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimitSeconds);
        model.optimize();

        result.status = model.get(GRB_IntAttr_Status);
        const bool hasIncumbent = model.get(GRB_IntAttr_SolCount) > 0;
        const bool validStatus =
            result.status == GRB_OPTIMAL ||
            result.status == GRB_USER_OBJ_LIMIT ||
            (result.status == GRB_TIME_LIMIT && hasIncumbent);

        if (validStatus) {
            result.success = true;
            result.W = W.get(GRB_DoubleAttr_X);

            {
                std::ostringstream ctx;
                ctx << "suffix=" << (fileSuffix.empty() ? "-" : fileSuffix)
                    << " excluded=" << excludedDemands.size()
                    << (warmStart != nullptr ? " re-solve" : " first-solve");
                selfCheckRingCrossings(b, r, ctx.str());
            }

            result.demandIL.assign(numDemands, 0.0);
            result.demandDistance.assign(numDemands, 0.0);
            result.demandBendCount.assign(numDemands, 0);

            int current = 0;
            result.tour.push_back(current);
            for (int step = 0; step < N - 1; step++) {
                bool found = false;
                for (int j = 0; j < N; j++) {
                    if (j == current) continue;
                    if (b[edge_idx(current,j)].get(GRB_DoubleAttr_X) > 0.5) {
                        int opt = (int)(r[edge_idx(current,j)].get(GRB_DoubleAttr_X) + 0.5);
                        result.tourEdges.push_back({current, j, opt});
                        current = j;
                        result.tour.push_back(current);
                        found = true;
                        break;
                    }
                }
                if (!found) break;
            }
            if (!result.tourEdges.empty()) {
                int first = result.tourEdges.front().from;
                int last = result.tourEdges.back().to;
                if (first != last) {
                    int opt = (int)(r[edge_idx(last, first)].get(GRB_DoubleAttr_X) + 0.5);
                    result.tourEdges.push_back({last, first, opt});
                }
            }

            result.demandFlowEdges.assign(numDemands, {});
            if (!optimizeWorstCase) {
                applyShortestArcDemandRouting(result, excludedDemands);
            } else {
                applyShortestArcDemandRouting(result, excludedDemands);

                double L_total_val = 0.0;
                for (const auto& te : result.tourEdges)
                    L_total_val += manhattanDist(te.from, te.to);

                verifyPositionModelSolution(
                    result, demandDistVar, forwardVar, pos,
                    L_total_val, excludedDemands);

                result.W = 0.0;
                for (int q = 0; q < numDemands; ++q) {
                    if (excludedDemands.count(q)) continue;
                    const double modelDist = demandDistVar[q].get(GRB_DoubleAttr_X);
                    result.demandDistance[q] = modelDist;
                    result.W = std::max(result.W, modelDist);
                }
            }

            if (!quiet) {
                std::cout << "Gefundener Worst-Case Distanz (W): " << result.W << " mm\n";
            }

            if (!fileSuffix.empty())
                exportCSV(result, fileSuffix, excludedDemands);
        } else if (!quiet) {
            std::cout << "Keine gültige Lösung. Status: " << result.status << "\n";
        }
    } catch (GRBException& e) {
        std::cerr << "Gurobi Fehler: " << e.getMessage() << "\n";
    }
    return result;
}

void MILPSolver::exportCSV(
        const MILPSolveResult& result,
        const std::string& fileSuffix,
        const std::set<int>& excludedDemands) const {
    int numDemands = (int)D.demands.size();

    std::ofstream csv("ring_" + fileSuffix + ".csv");
    csv << "from_id,to_id,from_x,from_y,to_x,to_y,routing\n";
    for (const auto& te : result.tourEdges) {
        csv << te.from + 1 << "," << te.to + 1 << ","
            << nodes[te.from].x << "," << nodes[te.from].y << ","
            << nodes[te.to].x << "," << nodes[te.to].y << ","
            << (te.routing == 0 ? "HV" : "VH") << "\n";
    }
    csv.close();

    std::ofstream ncsv("nodes.csv");
    ncsv << "id,x,y\n";
    for (const auto& n : nodes)
        ncsv << n.id + 1 << "," << n.x << "," << n.y << "\n";
    ncsv.close();

    std::ofstream dcsv("demands_" + fileSuffix + ".csv");
    dcsv << "demand,sender,receiver,from_id,to_id,direction\n";

    std::vector<int> tourPos(N, -1);
    const int tourLen = (int)result.tour.size();
    for (int i = 0; i < tourLen; ++i)
        tourPos[result.tour[i]] = i;

    auto tourStepIsForward = [&](int from, int to) {
        if (from < 0 || from >= N || to < 0 || to >= N) return true;
        if (tourPos[from] < 0 || tourPos[to] < 0) return true;
        return result.tour[(tourPos[from] + 1) % tourLen] == to;
    };

    for (int q = 0; q < numDemands; q++) {
        if (excludedDemands.count(q)) continue;
        int s = D.demands[q].first;
        int d = D.demands[q].second;

        if (q < (int)result.demandFlowEdges.size() && !result.demandFlowEdges[q].empty()) {
            const bool forward = tourStepIsForward(
                result.demandFlowEdges[q].front().first,
                result.demandFlowEdges[q].front().second);
            const std::string dir = forward ? "CW" : "CCW";
            for (const auto& [i, j] : result.demandFlowEdges[q]) {
                dcsv << q << ","
                     << s + 1 << "," << d + 1 << ","
                     << i + 1 << "," << j + 1 << ","
                     << dir << "\n";
            }
            continue;
        }

        std::cerr << "[exportCSV] demand " << q << " has no flow edges in result.\n";
    }
    dcsv.close();
    std::cout << "CSV exportiert für Methode " << fileSuffix << "\n";
}