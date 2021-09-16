#ifndef HIST_H
#define HIST_H

#include "utils.h"
#include <stdlib.h>

mat *buildHist(mat *buf, size_t row, size_t col);
void stretchHist(mat *dstMat, mat *srcMat, size_t row, size_t col);
void eqHist(mat *dstMat, mat *srcMat, size_t row, size_t col);

#endif // HIST_H
