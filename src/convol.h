#ifndef MATRIXMUL_H
#define MATRIXMUL_H

#include "utils.h"
#include <stdlib.h>

void convolN(mat *dstMat, mat *srcMat, size_t row, size_t col, mat *h, int divFactor, size_t n);
void binomFilConvolN(mat *dstMat, mat *srcMat, size_t row, size_t col, size_t n);

#endif // MATRIXMUL_H
