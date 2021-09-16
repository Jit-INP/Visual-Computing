#include "edge.h"
#include "arr.h"
#include "matrix.h"
#include "filter.h"
#include "convol.h"

void edgeDet(mat *imgMat, mat *hx, mat *hy, int divFactor, mat *xComp, mat *yComp, mat *edgeImg, size_t gradTh, size_t row, size_t col){

    convolN(xComp, imgMat, row, col, hx, divFactor, 1);
    convolN(yComp, imgMat, row, col, hy, divFactor, 1);

    mat *gradMag = allocMatMem(row, col, GRAYSCALE_DIM);
    getMagnitudeMat(xComp, yComp, row, col, gradMag);

    matrixThresholding(gradMag, row, col, gradTh, edgeImg);
    free(gradMag);
}

void sobelEdgeDet(mat *imgMat, mat *xComp, mat *yComp, mat *edgeImg, size_t gradTh, size_t row, size_t col){
    chkMatrixValidity(imgMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(xComp);
    NULL_PTR_CHK(yComp);
    NULL_PTR_CHK(edgeImg);
    SAME_PTR_CHK(imgMat, xComp);
    SAME_PTR_CHK(imgMat, yComp);
    SAME_PTR_CHK(imgMat, edgeImg);
    assert(gradTh <= (unsigned)MAX_PIX_VAL);

    mat *hy = allocFilterMem();
    mat *hx = allocFilterMem();
    int divFactor;

    getSobelOperators( hx, hy, &divFactor );

    edgeDet(imgMat, hx, hy, divFactor, xComp, yComp, edgeImg, gradTh, row, col);

    free(hx);
    free(hy);
}

void scharrEdgeDet(mat *imgMat, mat *xComp, mat *yComp, mat *edgeImg, size_t gradTh, size_t row, size_t col){
    chkMatrixValidity(imgMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(xComp);
    NULL_PTR_CHK(yComp);
    NULL_PTR_CHK(edgeImg);
    SAME_PTR_CHK(imgMat, xComp);
    SAME_PTR_CHK(imgMat, yComp);
    SAME_PTR_CHK(imgMat, edgeImg);
    assert(gradTh <= (unsigned)MAX_PIX_VAL);

    mat *hy = allocFilterMem();
    mat *hx = allocFilterMem();
    int divFactor;

    getScharrOperators( hx, hy, &divFactor );

    edgeDet(imgMat, hx, hy, divFactor, xComp, yComp, edgeImg, gradTh, row, col);

    free(hx);
    free(hy);

}


void nonMaxSuppression(mat *imgMat, mat *imgDirec, mat *nmsMat, size_t row, size_t col){
    chkMatrixValidity(imgMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(imgDirec);
    SAME_PTR_CHK(imgDirec, imgMat);
    SAME_PTR_CHK(imgDirec, nmsMat);
    SAME_PTR_CHK(imgMat, nmsMat);

    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            mat q, r;
            mat curImgMat = imgMat[(i * col) + j];
            if( ( ((signed)i - 1) < 0 ) || ( ((signed)j - 1) < 0 ) || ( ((signed)i + 1) > (signed)(row - 1) ) || ( ((signed)j + 1) > (signed)(col - 1) ) ){
                q = 255;
                r = 255;
            } else {
                mat curAng = imgDirec[(i * col) + j];
                if((curAng >= 0 && curAng < 22.5) || (curAng >= 157.5 && curAng <= 180)){
                    q = imgMat[(i * col) + (j + 1)];
                    r = imgMat[(i * col) + (j - 1)];
                } else if(curAng >= 22.5 && curAng < 67.5) {
                    q = imgMat[((i + 1) * col) + (j - 1)];
                    r = imgMat[((i - 1) * col) + (j + 1)];
                } else if(curAng >= 67.5 && curAng < 112.5) {
                    q = imgMat[((i + 1) * col) + j];
                    r = imgMat[((i - 1) * col) + j];
                } else if(curAng >= 112.5 && curAng < 157.5) {
                    q = imgMat[((i - 1) * col) + (j - 1)];
                    r = imgMat[((i + 1) * col) + (j + 1)];
                } else {
                    assert(0);
                }
            }
            if((curImgMat >= q) && (curImgMat >= r)){
                nmsMat[(i * col) + j] = curImgMat;
            } else {
                nmsMat[(i * col) + j] = 0;
            }
        }
    }
}

#define STRONG_PIX_VAL MAX_PIX_VAL
#define WEAK_PIX_VAL 25

void hysteresisThresholding(mat *nmsImg, mat *hysThImg, size_t row, size_t col, double lowThRatio, double highThRatio){
    chkMatrixValidity(nmsImg, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(hysThImg);
    SAME_PTR_CHK(nmsImg, hysThImg);
    ABOVE_ZERO_CHK(lowThRatio);
    ABOVE_ZERO_CHK(highThRatio);

    double highTh = getArrMax(nmsImg, row * col) * highThRatio;
    double lowTh = highTh * lowThRatio;
    assert(highTh > lowTh);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            mat curVal = nmsImg[(i * col) + j];
            if(curVal > highTh){
                hysThImg[(i * col) + j] = STRONG_PIX_VAL;
            } else if(curVal < lowTh){
                hysThImg[(i * col) + j] = 0;
            } else {
                hysThImg[(i * col) + j] = WEAK_PIX_VAL;
            }
        }
    }
}

void trackEdges(mat *hysThImg, mat *resImg, size_t row, size_t col){
    chkMatrixValidity(hysThImg, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(resImg);
    SAME_PTR_CHK(hysThImg, resImg);

    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            if( ( ((signed)i - 1) < 0 ) || ( ((signed)j - 1) < 0 ) || ( ((signed)i + 1) > (signed)(row - 1) ) || ( ((signed)j + 1) > (signed)(col - 1) ) ){
                resImg[(i * col) + j] = 0;
            } else {
                if(hysThImg[(i * col) + j] == WEAK_PIX_VAL){
                    if( (hysThImg[((i + 1) * col) + (j - 1)] == STRONG_PIX_VAL) || (hysThImg[( (i + 1) * col) + j] == STRONG_PIX_VAL) || (hysThImg[((i + 1) * col) + (j + 1)] == STRONG_PIX_VAL) ||
                            (hysThImg[(i * col) + (j - 1)] == STRONG_PIX_VAL) || (hysThImg[(i * col) + (j + 1)] == STRONG_PIX_VAL) ||
                            (hysThImg[((i - 1) * col) + (j - 1)] == STRONG_PIX_VAL) || (hysThImg[((i - 1) * col) + j] == STRONG_PIX_VAL) || (hysThImg[((i - 1) * col) + (j + 1)] == STRONG_PIX_VAL) ){
                        resImg[(i * col) + j] = STRONG_PIX_VAL;
                    } else {
                        resImg[(i * col) + j] = 0;
                    }
                } else if(hysThImg[(i * col) + j] == STRONG_PIX_VAL){
                    resImg[(i * col) + j] = STRONG_PIX_VAL;
                } else if(hysThImg[(i * col) + j] == 0){
                    resImg[(i * col) + j] = 0;
                }
            }
        }
    }
}

void cannyEdgeDet(mat *imgMat, mat *filMat, mat *xComp, mat *yComp, mat *gradThComp, mat *nmsImg, mat* hysThImg, mat *edgeImg, size_t row, size_t col){
    chkMatrixValidity(imgMat, row, col, GRAYSCALE_DIM);
    NULL_PTR_CHK(filMat);
    SAME_PTR_CHK(imgMat, filMat);
    NULL_PTR_CHK(xComp);
    SAME_PTR_CHK(imgMat, xComp);
    NULL_PTR_CHK(yComp);
    SAME_PTR_CHK(imgMat, yComp);
    NULL_PTR_CHK(gradThComp);
    SAME_PTR_CHK(imgMat, gradThComp);
    NULL_PTR_CHK(nmsImg);
    SAME_PTR_CHK(imgMat, nmsImg);
    NULL_PTR_CHK(hysThImg);
    SAME_PTR_CHK(imgMat, hysThImg);
    NULL_PTR_CHK(edgeImg);
    SAME_PTR_CHK(imgMat, edgeImg);

    binomFilConvolN(filMat, imgMat, row, col, 1);

    sobelEdgeDet(filMat, xComp, yComp, gradThComp, 0, row, col);

    mat *direc = allocMatMem(row, col, GRAYSCALE_DIM);
    getDirecMat(xComp, yComp, row, col, direc);

    nonMaxSuppression(gradThComp, direc, nmsImg, row, col);

    hysteresisThresholding(nmsImg, hysThImg, row, col, 0.05, 0.09);

    trackEdges(hysThImg, edgeImg, row, col);
}
