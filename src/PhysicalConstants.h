#pragma once

#include <cmath>

constexpr double ALPHA = 0.0274;
constexpr double BETA_THROUGH = 0.010;
constexpr double BETA_BEND = 0.010;
/// Waveguide Add/Drop coupling loss when a demand uses a shortcut endpoint it
/// does not own. Distinct from bend/through loss — TODO: calibrate.
constexpr double BETA_DROP = 0.5;
constexpr double MILP_EPS = 1e-6;

// Legacy alias (reporting / backward compat in non-decision paths).
constexpr double BETA = BETA_THROUGH;

inline bool ringEdgeNeedsBend(double xi, double yi, double xj, double yj) {
    return std::abs(xi - xj) > MILP_EPS && std::abs(yi - yj) > MILP_EPS;
}

// Reporting-only dB for one ring edge (not used for MILP W / worst-case decisions).
inline double reportingRingEdgeIL_dB(double xi, double yi, double xj, double yj) {
    const double dist = std::abs(xi - xj) + std::abs(yi - yj);
    double il = ALPHA * dist + BETA_THROUGH;
    if (ringEdgeNeedsBend(xi, yi, xj, yj))
        il += BETA_BEND;
    return il;
}

// Reporting-only dB for a shortcut path (not used for shortcut-vs-ring acceptance).
inline double reportingShortcutIL_dB(double distance, int bendCount) {
    return ALPHA * distance + BETA_THROUGH + BETA_BEND * static_cast<double>(bendCount);
}

// Legacy name kept for call sites that only need reporting IL.
inline double computeShortcutIL(double distance, int bendCount) {
    return reportingShortcutIL_dB(distance, bendCount);
}
