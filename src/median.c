#include "median.h"
#include "filter.h"
#include "matrix.h"

int medComp(const void *ele1, const void *ele2){
    mat f = *((mat*)ele1);
    mat s = *((mat*)ele2);
    if (f > s){
        return  1;
    }
    if (f < s){
        return -1;
    }
    return 0;
}

void sortMedianMatrix(mat *y){
    qsort(y, FILTER_SIZ * FILTER_SIZ, sizeof (mat), medComp);
}

mat getMedEle(mat *x){
    mat *y = (void *)allocFilterMem();
    cpyMatrix(y, x, FILTER_SIZ, FILTER_SIZ, GRAYSCALE_DIM);
    sortMedianMatrix(y);
    mat medEle = y[(FILTER_SIZ * FILTER_SIZ) / 2];
    free(y);
    return medEle;
}

mat medianXY(mat *buf, size_t x, size_t y, size_t row, size_t col){
    mat *tmp = allocFilterMem();
    buildFilterSizMatrix(tmp, buf, x, y, row, col);
    mat med = getMedEle(tmp);
    free(tmp);
    return med;
}

void median(mat *dstMat, mat *srcMat, size_t row, size_t col){
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            *(dstMat + ( (i * col) + j) ) = medianXY(srcMat, i, j, row, col);
        }
    }
    return;
}

void medianN(mat *dstMat, mat *srcMat, size_t row, size_t col, size_t n){
    chkMatrixValidity(srcMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(dstMat);
    ABOVE_ZERO_CHK(n);

    mat *tmpBuff = allocMatMem(row, col, GRAYSCALE_DIM);
    cpyMatrix(tmpBuff, srcMat, row, col, GRAYSCALE_DIM);
    for (size_t i = 0; i < n; ++i) {
        median(dstMat, tmpBuff, row, col);
        cpyMatrix(tmpBuff, dstMat, row, col, GRAYSCALE_DIM);
    }
    free(tmpBuff);
    return;
}
