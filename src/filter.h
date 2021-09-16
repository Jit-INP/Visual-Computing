#ifndef FILTER_H
#define FILTER_H


#include "utils.h"
#include <stdlib.h>

void revFilter(mat *dst, mat *src);
void prodWFilter(mat *x, mat *h, mat *y);
void *allocFilterMem( void );
mat sumBuff(mat *buff, size_t limit);
mat summateFilter(mat *m);
void buildFilterSizMatrix(mat *op, mat *buf, size_t x, size_t y, size_t row, size_t col);
void genGausianKernel(mat *gKernel, int *divFactor);
void getScharrOperators(mat *hx, mat *hy, int *divFactor);
void getSobelOperators(mat *hx, mat *hy, int *divFactor);

#endif // FILTER_H
