#include <string.h>
#include "convol.h"
#include "filter.h"
#include "matrix.h"

mat getConvVal(mat *x, mat *revH, mat divFactor){
    mat *y = (void *)allocFilterMem();
    prodWFilter( x, revH, y);
    int sum = summateFilter(y);
    mat convPix = (mat) ((double) sum / divFactor);
    free(y);
    return convPix;
}

mat convXY(mat *srcMat, size_t x, size_t y, size_t row, size_t col, mat *revH, int divFactor){
    mat *tmp = allocFilterMem();
    buildFilterSizMatrix(tmp, srcMat, x, y, row, col);
    mat conv = getConvVal(tmp, revH, divFactor);
    free(tmp);
    return conv;
}

void convol(mat *dstMat, mat *srcMat, size_t row, size_t col, mat *h, int divFactor){
    mat *revH = allocFilterMem();
    revFilter(revH, h);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            *(dstMat + ( (i * col) + j) ) = convXY(srcMat, i, j, row, col, revH, divFactor);
        }
    }
    free(revH);
    return;
}

void convolN(mat *dstMat, mat *srcMat, size_t row, size_t col, mat *h, int divFactor, size_t n){
    chkMatrixValidity(srcMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(dstMat);
    ABOVE_ZERO_CHK(divFactor);
    ABOVE_ZERO_CHK(n);
    SAME_PTR_CHK(dstMat, srcMat);

    mat *tmpBuff = allocMatMem(row, col, GRAYSCALE_DIM);
    cpyMatrix(tmpBuff, srcMat, row, col, GRAYSCALE_DIM);
    for (size_t i = 0; i < n; ++i) {
        convol(dstMat, tmpBuff, row, col, h, divFactor);
        cpyMatrix(tmpBuff, dstMat, row, col, GRAYSCALE_DIM);
    }
    free(tmpBuff);
    return;
}

void binomFilConvolN(mat *dstMat, mat *srcMat, size_t row, size_t col, size_t n){
    mat *h = allocFilterMem();
    int divFactor;
    genGausianKernel( h, &divFactor);
    convolN(dstMat, srcMat, row, col, h, divFactor, n);
}
