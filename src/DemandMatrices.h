#pragma once

#include "DemandMatrix.h"

// density: 1 = low (~1 demand per node), 2 = high (~2 demands per node)
// numNodes: 8, 12, 14, 16, 20, 24, 28, or 32
bool initDemandMatrix(int numNodes, int density, DemandMatrix& D);

void initDemandMatrix_8_Density_1(DemandMatrix& D);
void initDemandMatrix_8_Density_2(DemandMatrix& D);
void initDemandMatrix_12_Density_1(DemandMatrix& D);
void initDemandMatrix_12_Density_2(DemandMatrix& D);
void initDemandMatrix_14_Density_1(DemandMatrix& D);
void initDemandMatrix_16_Density_1(DemandMatrix& D);
void initDemandMatrix_16_Density_2(DemandMatrix& D);
void initDemandMatrix_20_Density_1(DemandMatrix& D);
void initDemandMatrix_20_Density_2(DemandMatrix& D);
void initDemandMatrix_24_Density_1(DemandMatrix& D);
void initDemandMatrix_24_Density_2(DemandMatrix& D);
void initDemandMatrix_28_Density_1(DemandMatrix& D);
void initDemandMatrix_28_Density_2(DemandMatrix& D);
void initDemandMatrix_32_Density_1(DemandMatrix& D);
void initDemandMatrix_32_Density_2(DemandMatrix& D);
