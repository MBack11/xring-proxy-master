#pragma once

#include "gurobi_c++.h"

/// Process-wide Gurobi environment (avoids repeated env.start() cost).
inline GRBEnv& sharedGurobiEnv() {
    static GRBEnv* env = nullptr;
    if (!env) {
        env = new GRBEnv(true);
        env->set(GRB_IntParam_OutputFlag, 0);
        env->start();
    }
    return *env;
}
