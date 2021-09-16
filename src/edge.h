#ifndef EDGE_H
#define EDGE_H

#include "utils.h"
#include <stdlib.h>

void sobelEdgeDet(mat *imgMat, mat *xComp, mat *yComp, mat *edgeImg, size_t gradTh, size_t row, size_t col);
void scharrEdgeDet(mat *imgMat, mat *xComp, mat *yComp, mat *edgeImg, size_t gradTh, size_t row, size_t col);
void cannyEdgeDet(mat *imgMat, mat *filMat, mat *xComp, mat *yComp, mat *gradThComp, mat *nmsImg, mat *hysThImg, mat *edgeImg, size_t row, size_t col);

#endif // EDGE_H
