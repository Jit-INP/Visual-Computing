#include "filter.h"
#include <string.h>
#include "matrix.h"

void revFilter(mat *dst, mat *src){
    reverseMatrix(dst, src, FILTER_SIZ, FILTER_SIZ);
}

void prodWFilter(mat *x, mat *h, mat *y){
    prodMatrix(x, h, FILTER_SIZ, FILTER_SIZ, y);
}

void *allocFilterMem( void ){
    return allocMatMem(FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
}

mat summateFilter(mat *m){
    NULL_PTR_CHK(m);

    mat sum = 0;
    for (size_t i = 0; i < FILTER_SIZ; ++i) {
        for (size_t j = 0; j < FILTER_SIZ; ++j) {
            sum += m[ (i * FILTER_SIZ) + j];
        }
    }
    return sum;
}

mat sumBuff(mat *buff, size_t limit){
    NULL_PTR_CHK(buff);

    mat sum = 0;
    for (size_t i = 0; i <= limit; ++i) {
        sum += buff[i];
    }
    return sum;
}

void buildFilterSizMatrix(mat *op, mat *buf, size_t x, size_t y, size_t row, size_t col){
    chkMatrixValidity(buf, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(op);
    SAME_PTR_CHK(buf, op);

    int mid = FILTER_SIZ / 2;
    for (int i = -mid; i <= mid; ++i) {
        for (int j = -mid; j <= mid; ++j) {
            if( ( ((signed)x + i) < 0 ) || ( ((signed)y + j) < 0 ) || ( ((signed)x + i) > (signed)(row - 1) ) || ( ((signed)y + j) > (signed)(col - 1) ) ){
                op[( (i + mid) * FILTER_SIZ) + (j + mid)] = 0;
            } else {
                op[( (i + mid) * FILTER_SIZ) + (j + mid)] = buf [ ( ((signed)x + i) * (signed)col ) + ((signed)y + j)];
            }
        }
    }
}

void genPascalNbs(mat *pascBuf){
    mat *pascMat = allocFilterMem();
    mat initVal = 1;
    memsetMatrix(pascMat, &initVal, FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    for (size_t i = 1; i < FILTER_SIZ; ++i) {
        for (size_t j = 0; j < FILTER_SIZ; ++j) {
            pascMat[(i * FILTER_SIZ) + j] = (mat) sumBuff(pascMat + ( (i - 1) * FILTER_SIZ), j);
        }
    }
    for (size_t k = 0; k < FILTER_SIZ; ++k){
        size_t i = FILTER_SIZ - 1 - k;
        size_t j = k;
        pascBuf[k] = pascMat[(i * FILTER_SIZ) + j];
    }
    free(pascMat);
}

void genGausianKernel(mat *gKernel, int *divFactor){
    NULL_PTR_CHK(gKernel);
    NULL_PTR_CHK(divFactor);

    mat *pascNb = malloc(FILTER_SIZ * sizeof (mat));
    genPascalNbs(pascNb);
    mat *buf;
    size_t tmp;
    matrixMul(pascNb, FILTER_SIZ, 1, pascNb, 1, FILTER_SIZ, &buf, &tmp, &tmp);
    for (size_t i = 0; i < FILTER_SIZ; ++i) {
        for (size_t j = 0; j < FILTER_SIZ; ++j) {
            gKernel[ (i * FILTER_SIZ) + j] = buf[ (i * FILTER_SIZ) + j];
        }
    }
    free(buf);
    int sum = sumBuff(pascNb, FILTER_SIZ - 1);
    *divFactor = sum * sum;
    free(pascNb);
}

void getScharrOperators(mat *hx, mat *hy, int *divFactor){
    NULL_PTR_CHK(hx);
    NULL_PTR_CHK(hy);
    NULL_PTR_CHK(divFactor);
    mat scharrX[9] = {-3, 0, 3, -10, 0, 10, -3, 0, 3};
    mat scharrY[9] = {-3, -10, -3, 0, 0, 0, 3, 10, 3};
    cpyMatrix(hx, &scharrX[0], FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    cpyMatrix(hy, &scharrY[0], FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    *divFactor = 16;
}

void getSobelOperators(mat *hx, mat *hy, int *divFactor){
    NULL_PTR_CHK(hx);
    NULL_PTR_CHK(hy);
    NULL_PTR_CHK(divFactor);
    mat sobelX[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    mat sobelY[9] = {1, 2, 1, 0, 0, 0, -1, -2, -1};
    cpyMatrix(hx, &sobelX[0], FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    cpyMatrix(hy, &sobelY[0], FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    *divFactor = 4;
}
