#include "kmeans.h"
#include<math.h>
#include<string.h>
#include "arr.h"
#include "matrix.h"
#include "lst.h"

unsigned int calcDist(mat* a, mat* b, size_t dim){
    unsigned int sqSum = 0;
    for (size_t k = 0; k < dim; ++k) {
        int d = a[k] - b[k];
        sqSum = sqSum + ((unsigned)(d * d));
    }
    unsigned int dist = (unsigned int) +sqrt(sqSum);
    return dist;
}

size_t getNearestCentroidIndex(mat *centLst, mat* curPix, size_t K, size_t dim){
    chkLstValidity(centLst, K, dim);
    mat *distLst = allocLstMem(K, 1);

    mat initVal = 0;
    memsetLst(distLst, &initVal, K, 1);

    for (size_t i = 0; i < 1; ++i) {
        for (size_t j = 0; j < K; ++j) {
            distLst[j] = (signed) calcDist(&centLst[(i * K * dim) + (j * dim)], curPix, dim);
        }
    }
    size_t idx = getArrMinEleIndex((mat*)distLst, K);

    free(distLst);
    return idx;
}

mat **allocClusterLst(size_t clustSiz, size_t nbClust){
    mat **cluster = malloc(sizeof (mat*) * nbClust);
    for (size_t i = 0; i < nbClust; ++i) {
        cluster[i] = allocLstMem(clustSiz, 1);
    }
    return cluster;
}

void deallocClusterLst(mat **clusterLst, size_t nbClust){
    for (size_t i = 0; i < nbClust; ++i) {
        free(clusterLst[i]);
    }
    free(clusterLst);
}

void buildClusters(mat** clustLst, mat* clustLstCurIdx, mat* centLst, mat* imgMat, size_t row, size_t col, size_t dim, size_t K){
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            mat curPixIdx = (signed)(i * col * dim) + (signed)(j * dim); // index of start of current pix in imgMat array
            mat *curPix = &imgMat[curPixIdx];
            size_t nCentIdx = getNearestCentroidIndex(centLst, curPix, K, dim);
            addToArray( clustLst[nCentIdx], row * col, &clustLstCurIdx[nCentIdx], (signed)curPixIdx );
        }
    }
}

void buildMatForAvg(mat *matr, mat* imgMat, mat *clust, size_t nbEle, size_t dim){
    for (size_t i = 0; i < nbEle; ++i) {
        cpyLst( &matr[i * dim], &imgMat[clust[i]], 1, dim);
    }
}

void getNewCentroids(mat** clustLst, mat* clustLstCurIdx, mat* imgMat, size_t dim, mat *newCentLst, size_t nbClust){
    for (size_t i = 0; i < nbClust; ++i) {
        mat *curClust = clustLst[i];
        if(clustLstCurIdx[i] != 0){
            mat *avgMat = allocLstMem((unsigned) clustLstCurIdx[i], dim);
            buildMatForAvg(avgMat, imgMat, curClust, (unsigned)clustLstCurIdx[i], dim);
            lstAvg(avgMat, (unsigned) clustLstCurIdx[i], dim, &newCentLst[i * dim]);
            free(avgMat);
        }
    }
}

void buildOpMat(mat *opMat, mat **clustLst, mat *clustLstCurIdx, mat *centLst, size_t dim, size_t nbClust){
    for (size_t i = 0; i < nbClust; ++i) {
        size_t curClustSiz = (unsigned) clustLstCurIdx[i];
        if(curClustSiz == 0){
            continue;
        }
        mat *curClust = clustLst[i];
        mat *curCent = &centLst[i * dim];
        for (size_t j = 0; j < curClustSiz; ++j) {
            mat curIdx = curClust[j];
            cpyLst(&opMat[(unsigned)curIdx], curCent, 1, dim);
        }
    }
}

void applyKMeans(mat* opMat, mat *imgMat, size_t row, size_t col, size_t dim, size_t K){
    chkMatrixValidity(imgMat, row, col, dim);
    NULL_PTR_CHK(opMat);
    SAME_PTR_CHK(imgMat, opMat);
    ABOVE_ZERO_CHK(K);

//    opMat = allocMatMem(row, col, dim);

//    imgMat = genRandLst(row * col, dim, MAX_PIX_VAL);
//    printMatrix(imgMat, row, col, dim);
    mat *maxLst = allocLstMem(1, dim);
    mat initVal = MAX_PIX_VAL;
    memsetLst(maxLst, &initVal, dim, 1);
    mat *centLst = genRandLst(K, dim, maxLst);
    free(maxLst);
    mat *newCentLst = allocLstMem(K, dim);

    bool stopItr = false;
    while(stopItr == false){
        mat **clustLst = allocClusterLst(row * col, K);

        mat *clustLstCurIdx = allocLstMem(K, 1);
        mat initClustIdxVal = 0;
        memsetLst(clustLstCurIdx, &initClustIdxVal, K, 1);

//        printLst(centLst, K, dim);
        buildClusters(clustLst, clustLstCurIdx, centLst, imgMat, row, col, dim, K);

//        for (size_t i = 0; i < K; ++i) {
//            printLst(clustLst[i], (size_t)clustLstCurIdx[i], 1);
//        }

        getNewCentroids(clustLst, clustLstCurIdx, imgMat, dim, newCentLst, K);

        if( ( stopItr = cmpLst(centLst, newCentLst, K, dim) ) == true){
            buildOpMat(opMat, clustLst, clustLstCurIdx, centLst, dim, K);
        }

        free(clustLstCurIdx);
        deallocClusterLst(clustLst, K);
        cpyLst(centLst, newCentLst, K, dim);
    }

//    printMatrix(opMat, row, col, dim);
    free(newCentLst);
    free(centLst);
}


mat *imgMatToImgSpatialMat(mat *imgMat, size_t row, size_t col, size_t *dim){
    size_t colDim = *dim;
    *dim += 2;
    mat *spatImgMat = allocMatMem(row, col, *dim);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            for (size_t k = 0; k < *dim; ++k) {
                size_t spatImgMatIdx = (i * col * *dim) + (j * *dim) + k;
                if(k < colDim){
                    size_t imgMatIdx = spatImgMatIdx - (2 * i * col) - (2 * j); // - 2 ( ( i * col) + j)
                    spatImgMat[spatImgMatIdx] = imgMat[imgMatIdx];
                } else if (k == colDim){
                    spatImgMat[spatImgMatIdx] = (signed)i;
                } else {
                    spatImgMat[spatImgMatIdx] = (signed)j;
                }
            }
        }
    }
    return spatImgMat;
}


void imgSpatialMatToImgMat(mat *imgMat, mat *spatImgMat, size_t row, size_t col, size_t *dim){
    *dim -= 2;
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            for (size_t k = 0; k < *dim; ++k) {
                size_t imgMatIdx = (i * col * *dim) + (j * *dim) + k;
                size_t spatImgMatIdx = imgMatIdx + (2 * i * col) + (2 * j); // + 2 ( ( i * col) + j)
                imgMat[imgMatIdx] = spatImgMat[spatImgMatIdx];
            }
        }
    }
}

void applySpatialKMeans(mat* opMat, mat *imgMat, size_t row, size_t col, size_t dim, size_t K){
//        row = 2;
//        col = 2;
//        dim = 1;
//        opMat = allocMatMem(row, col, dim);

//        imgMat = genRandLst(row * col, dim, MAX_PIX_VAL);
//        printMatrix(imgMat, row, col, dim);
//        mat *spatImgMat = imgMatToImgSpatialMat(imgMat, row, col, &dim);
//        printMatrix(spatImgMat, row, col, dim);
//        imgSpatialMatToImgMat(opMat, spatImgMat, row, col, &dim);
//        printMatrix(opMat, row, col, dim);

    chkMatrixValidity(imgMat, row, col, dim);
    NULL_PTR_CHK(opMat);
    SAME_PTR_CHK(imgMat, opMat);
    ABOVE_ZERO_CHK(K);
    size_t colDim = dim;
    mat *spatImgMat = imgMatToImgSpatialMat(imgMat, row, col, &dim);
    mat *spatImgOpMat = allocMatMem(row, col, dim);

    mat *maxLst = allocLstMem(1, dim);
    mat initVal = MAX_PIX_VAL;
    memsetLst(maxLst, &initVal, colDim, 1);
    maxLst[colDim] = (signed)row;
    maxLst[colDim + 1] = (signed)col;
    mat *centLst = genRandLst(K, dim, maxLst);
    free(maxLst);
    mat *newCentLst = allocLstMem(K, dim);

    bool stopItr = false;
    while(stopItr == false){
        mat **clustLst = allocClusterLst(row * col, K);

        mat *clustLstCurIdx = allocLstMem(K, 1);
        mat initClustIdxVal = 0;
        memsetLst(clustLstCurIdx, &initClustIdxVal, K, 1);

//        printLst(centLst, K, dim);
        buildClusters(clustLst, clustLstCurIdx, centLst, spatImgMat, row, col, dim, K);

//        for (size_t i = 0; i < K; ++i) {
//            printLst(clustLst[i], (size_t)clustLstCurIdx[i], 1);
//        }

        getNewCentroids(clustLst, clustLstCurIdx, spatImgMat, dim, newCentLst, K);

        if( ( stopItr = cmpLst(centLst, newCentLst, K, dim) ) == true){
            buildOpMat(spatImgOpMat, clustLst, clustLstCurIdx, centLst, dim, K);
            imgSpatialMatToImgMat(opMat, spatImgOpMat, row, col, &dim);
            free(spatImgMat);
            free(spatImgOpMat);
        }
        free(clustLstCurIdx);
        deallocClusterLst(clustLst, K);
        cpyLst(centLst, newCentLst, K, dim);
    }
//    printMatrix(opMat, row, col, dim);
    free(newCentLst);
    free(centLst);
}
