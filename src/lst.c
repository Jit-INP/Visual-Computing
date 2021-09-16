#include "lst.h"
#include "matrix.h"

void chkLstValidity(mat *lst, size_t nbEle, size_t dim){
    chkMatrixValidity(lst, 1, nbEle, dim);
}

mat *allocLstMem(size_t nbEle, size_t dim){
    return allocMatMem(1, nbEle, dim);
}

bool cmpLst(mat *a, mat *b, size_t nbEle, size_t dim){
    return matCmp(a, b, 1, nbEle, dim);
}

void printLst(mat *lst, size_t nbEle, size_t dim){
    printMatrix(lst, 1, nbEle, dim);
}

void cpyLst(mat *dstLst, mat *srcLst, size_t nbEle, size_t dim){
    cpyMatrix(dstLst, srcLst, 1, nbEle, dim);
}

void memsetLst(mat *lst, mat *val, size_t nbEle, size_t dim){
    memsetMatrix(lst, val, 1, nbEle, dim);
}


mat* genRandLst(size_t nbEle, size_t dim, mat *maxVal){
    return genRandMat(1, nbEle, dim, maxVal);
}

void lstAvg(mat *lst, size_t nbEle, size_t dim, mat *res){
    matAvg(lst, 1, nbEle, dim, res);
}
