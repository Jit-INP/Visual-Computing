#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"
#include <stdlib.h>
#include <stdbool.h>

void chkMatrixValidity(mat *matr, size_t row, size_t col, size_t dim);
void printMatrix(mat *buf, size_t row, size_t col, size_t dim);
mat *allocMatMem(size_t r, size_t c, size_t dim);
void memsetMatrix(mat *matr, mat *val, size_t row, size_t col, size_t dim);
void cpyMatrix(mat *dst, mat *src, size_t row, size_t col, size_t dim);
void reverseMatrix(mat *revH, mat *h, size_t r, size_t c);
void prodMatrix(mat *x, mat *h, size_t r, size_t c, mat *y);
void getMagnitudeMat(mat *xComp, mat *yComp, size_t row, size_t col, mat *mag);
void matrixMul(mat *matA, size_t rowA, size_t colA, mat *matB, size_t rowB, size_t colB, mat **matC, size_t *rowC, size_t *colC);
void matrixThresholding(mat *x, size_t r, size_t c, size_t th, mat *y);
void getDirecMat(mat *xComp, mat *yComp, size_t row, size_t col, mat *direc);
void matAvg(mat *matr, size_t row, size_t col, size_t dim, mat* res);
bool matCmp(mat *mat1, mat* mat2, size_t row, size_t col, size_t dim);
mat* genRandMat(size_t row, size_t col, size_t dim, mat *maxVal);

#endif // MATRIX_H
