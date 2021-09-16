#include "hist.h"
#include <string.h>
#include "matrix.h"
#include "filter.h"

mat *buildHist(mat *buf, size_t row, size_t col){
    chkMatrixValidity(buf, row, col, GRAYSCALE_DIM);
    mat *hist = malloc(sizeof (mat) * (unsigned)(MAX_PIX_VAL + 1));
    mat initVal = 0;
    memsetMatrix(hist, &initVal, 1, (unsigned)(MAX_PIX_VAL + 1), GRAYSCALE_DIM);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            hist[buf[(i * col) + j]]++;
        }
    }
    return hist;
}

int findMinIndex(mat *buf, size_t nbEle){
    int minI = (signed)nbEle;
    for (size_t i = 0; i < nbEle; ++i) {
        if(buf[i] != 0){
            minI = (signed)i;
            break;
        }
    }
    assert(minI != (signed)nbEle);
    return minI;
}

int findMaxIndex(mat *buf, size_t nbEle){
    int maxI = (signed)nbEle;
    for (size_t i = nbEle - 1; i != 0; --i) {
        if(buf[i] != 0){
            maxI = (signed)i;
            break;
        }
    }
    assert(maxI != (signed)nbEle);
    return maxI;
}

void stretchHist(mat *dstMat, mat *srcMat, size_t row, size_t col){
    chkMatrixValidity(srcMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(dstMat);
    SAME_PTR_CHK(dstMat, srcMat);

    mat *hist = buildHist(srcMat, row, col);
    int iMin = findMinIndex(hist, (unsigned)(MAX_PIX_VAL + 1));
    int iMax = findMaxIndex(hist, (unsigned)(MAX_PIX_VAL + 1));
    int iDel = iMax - iMin;
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            mat iCur = srcMat[ (i * col) + j ];
            signed int val = (int)( ( (iCur - iMin) * MAX_PIX_VAL ) / (double) iDel );
            dstMat[ (i * col) + j ] = val;
        }
    }
    free(hist);
}

mat histEqXY(mat *histBuf, mat *imgBuf, size_t x, size_t y, size_t row, size_t col ){
    mat val = imgBuf[(x * col) + y];
    assert(val > -1);
    int sum = sumBuff(histBuf, (unsigned)val);
    return (mat) ( ( (double)(sum * MAX_PIX_VAL) ) / (row * col) );
}

void eqHist(mat *dstMat, mat *srcMat, size_t row, size_t col){
    chkMatrixValidity(srcMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(dstMat);
    SAME_PTR_CHK(dstMat, srcMat);

    mat *hist = buildHist(srcMat, row, col);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            dstMat[(i * col) + j] = histEqXY(hist, srcMat, i, j, row, col);
        }
    }

    free(hist);
    return;
}
