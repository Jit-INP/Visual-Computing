#ifndef LST_H
#define LST_H

#include "utils.h"
#include <stdbool.h>

void chkLstValidity(mat *lst, size_t nbEle, size_t dim);
mat *allocLstMem(size_t nbEle, size_t dim);
bool cmpLst(mat *a, mat *b, size_t nbEle, size_t dim);
void printLst(mat *lst, size_t nbEle, size_t dim);
void cpyLst(mat *dstLst, mat *srcLst, size_t nbEle, size_t dim);
void memsetLst(mat *lst, mat *val, size_t nbEle, size_t dim);
mat* genRandLst(size_t nbEle, size_t dim, mat *maxVal);
void lstAvg(mat *lst, size_t nbEle, size_t dim, mat *res);

#endif // LST_H
