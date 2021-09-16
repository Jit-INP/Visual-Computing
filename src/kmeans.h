#ifndef KMEANS_H
#define KMEANS_H

#include "utils.h"
#include <stdlib.h>

void applyKMeans(mat* opMat, mat *imgMat, size_t row, size_t col, size_t dim, size_t K);
void applySpatialKMeans(mat* opMat, mat *imgMat, size_t row, size_t col, size_t dim, size_t K);

#endif // KMEANS_H
