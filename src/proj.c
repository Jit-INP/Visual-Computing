#include "proj.h"

void centerThePCL(point3d *points, int N){
    int i = 0;
    double xMin = 100000, xMax = -1000000, yMin = 1000000, yMax = -10000000, zMin = 100000000, zMax = -10000000;
    for(i = 0; i < N; i++){
        if(points[i].x > xMax)
            xMax = points[i].x;
        else if (points[i].x < xMin)
            xMin = points[i].x;

        if(points[i].y > yMax)
            yMax = points[i].y;
        else if (points[i].y < yMin)
            yMin = points[i].y;

        if(points[i].z > zMax)
            zMax = points[i].z;
        else if (points[i].z < zMin)
            zMin = points[i].z;
    }

    double xMid = (xMin+xMax)/2, yMid = (yMin+yMax)/2, zMid = (zMin+zMax)/2;
    for(i = 0; i < N; i++){
        points[i].x = points[i].x-xMid;
        points[i].y = points[i].y-yMid;
        points[i].z = points[i].z-zMid;
    }
}


void Rx(double alpha, double * R){

    int i = 0, j = 0;
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            R[i*4+j] = 0;
    // 1      0     0    0
    //   0     cosA -sinA  0
    //   0     sinA  cosA  0
    //   0      0     0    1

    R[0] = 1;

    R[5] = cos(alpha);
    R[6] = -sin(alpha);

    R[9] = sin(alpha);
    R[10] = cos(alpha);

    R[15] = 1;

    return;
}

void Ry(double alpha, double *R){

    int i = 0, j = 0;
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            R[i*4+j] = 0;
    /* cosA    0    sinA   0
      0      1     0     0
    -sinA    0    cosA   0
      0      0     0     1
      */
    R[0] = cos(alpha);
    R[2] = sin(alpha);

    R[5] = 1;

    R[8] = -sin(alpha);
    R[10] = cos(alpha);

    R[15] = 1;

    return;
}

void Rz(double alpha, double *R){

    int i = 0, j = 0;
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            R[i*4+j] = 0;
    /* cosA  -sinA   0   0
     sinA   cosA   0   0
      0      0     1   0
      0      0     0   1
      */
    R[0] = cos(alpha);
    R[1] = -sin(alpha);

    R[4] = sin(alpha);
    R[5] = cos(alpha);

    R[10] = 1;
    R[15] = 1;

    return;
}

void getMatAEleD(double *matr, size_t col, size_t rowNo, double *arr){
    for (size_t j = 0; j < col; ++j) {
        arr[j] = matr[ (rowNo * col) + j];
    }
}

void getMatBEleD(double *matr, size_t nbRow, size_t nbCol, size_t colNo, double *arr){
    for (size_t i = 0; i < nbRow; ++i) {
        arr[i] = matr[ (i * nbCol) + colNo];
    }
}

double multEleD(double *matA, double *matB, size_t nbEle){
    double sum = 0;
    for (size_t i = 0; i < nbEle; ++i) {
        sum += (matA[i] * matB[i]);
    }
    return sum;
}

void matMul(double *matA, size_t rowA, size_t colA, double *matB, size_t rowB, size_t colB, double **matC, size_t *rowC, size_t *colC){
    chkMatrixValidity((void *)matA, rowA, colA, GRAYSCALE_DIM);
    chkMatrixValidity((void *)matB, rowB, colB, GRAYSCALE_DIM);
    NULL_PTR_CHK(matC);
    NULL_PTR_CHK(rowC);
    NULL_PTR_CHK(colC);
    SAME_PTR_CHK(*matC,(void*)matA);
    SAME_PTR_CHK(*matC,(void*)matB);
    assert(colA == rowB);

    *rowC = rowA;
    *colC = colB;
    *matC = (void *)malloc(rowA * colB * sizeof(double));

    double *p = *matC;

    for (size_t i = 0; i < rowA; ++i) {
        for (size_t j = 0; j < colB; ++j) {
            double *arrA = malloc(1 * colA * sizeof(double));
            getMatAEleD(matA, colA, i, arrA);

            double *arrB = malloc(rowB * 1 * sizeof(double));
            getMatBEleD(matB, rowB, colB, j, arrB);

            p[(i * colB) + j] = multEleD(arrA, arrB, colA);

            free(arrA);
            free(arrB);
        }
    }

}

void computeTrans(double gama, double beta, double alpha, double T_x, double T_y, double T_z, double **result){
    double R_x[16], R_y[16], R_z[16];
    Rx(gama, R_x);
    Ry(beta, R_y);
    Rz(alpha, R_z);

    double *temp;
    size_t tempR, tempC;
    matMul(R_y, 4, 4, R_z, 4, 4, &temp, &tempR, &tempC);
    matMul(R_x, 4, 4, temp, 4, 4, result, &tempR, &tempC);
    free(temp);
    *( (*result) + 3 ) = T_x;
    *( (*result) + 7 ) = T_y;
    *( (*result) + 11 ) = T_z;
    return;
}

typedef struct{
    double x;
    double y;
    double z;
    double d;
    mat r;
    mat g;
    mat b;
}RigTrans;

RigTrans *rigTransIn3d(point3d *pt, size_t N, double gama, double beta, double alpha, double T_x, double T_y, double T_z){
    RigTrans *t = malloc(N * sizeof (RigTrans));
    double *transMat;
    computeTrans(gama, beta, alpha, T_x, T_y, T_z, &transMat);
    for (size_t n = 0; n < N; ++n) {
        double *ptMat = malloc(4 * 1 * sizeof (double));
        ptMat[0] = pt[n].x;
        ptMat[1] = pt[n].y;
        ptMat[2] = pt[n].z;
        ptMat[3] = 1;
        double *resMat = 0;
        size_t tempR, tempC;
        matMul(transMat, 4, 4, ptMat, 4, 1, &resMat, &tempR, &tempC);
        t[n].x = resMat[0];
        t[n].y = resMat[1];
        t[n].z = resMat[2];
        t[n].d = resMat[3];
        t[n].r = pt[n].r;
        t[n].g = pt[n].g;
        t[n].b = pt[n].b;
        free(ptMat);
        free(resMat);
    }
    free(transMat);
    return t;
}

RigTrans *rigTransToPerspec(RigTrans *r, size_t N, double f){
    RigTrans *t = malloc(N * sizeof (RigTrans));
    for (size_t n = 0; n < N; ++n) {
        t[n].x = r[n].x;
        t[n].y = r[n].y;
        t[n].z = 0;
        t[n].d = 1 + (r[n].z / f);
        t[n].r = r[n].r;
        t[n].g = r[n].g;
        t[n].b = r[n].b;
    }
    return t;
}

typedef struct{
    double x;
    double y;
    mat r;
    mat g;
    mat b;
}PerspecImg;

PerspecImg *perspecToImage(RigTrans *r, size_t N){
    PerspecImg *p = malloc(N * sizeof (PerspecImg));
    for (size_t n = 0; n < N; ++n) {
        p[n].x = r[n].x / r[n].d;
        p[n].y = r[n].y / r[n].d;
        p[n].r = r[n].r;
        p[n].g = r[n].g;
        p[n].b = r[n].b;
    }
    return p;
}

PerspecImg *scalePerspecImage(PerspecImg *p, size_t N, size_t row, size_t col){
    PerspecImg *s = malloc(N * sizeof (PerspecImg));
    for (size_t n = 0; n < N; ++n) {
        s[n].x = (p[n].x * row) + row / 2;
        s[n].y = (p[n].y * col) + col / 2;
        s[n].r = p[n].r;
        s[n].g = p[n].g;
        s[n].b = p[n].b;
    }
    return s;
}

void prespecToImg(PerspecImg *pt, size_t N, mat *img, size_t row, size_t col){
    size_t dim = 3;
    mat val[3] = {0, 0, 0};
    memsetMatrix(img, val, row, col, dim);
    for (size_t n = 0; n < N; ++n) {
        if((pt[n].x < 0) || (pt[n].x > row)){
            continue;
        }
        if((pt[n].y < 0) || (pt[n].y > col)){
            continue;
        }
        size_t rowNo = (unsigned)pt[n].x;
        size_t colNo = (unsigned)pt[n].y;
        img[ ( rowNo * col * dim ) + ( colNo * dim ) + 0 ] = pt[n].r;
        img[ ( rowNo * col * dim ) + ( colNo * dim ) + 1 ] = pt[n].g;
        img[ ( rowNo * col * dim ) + ( colNo * dim ) + 2 ] = pt[n].b;
    }
}

PerspecImg *perspecToUV(PerspecImg *pt, size_t N, size_t resolX, size_t resolY, size_t x0, size_t y0){
    PerspecImg *uv = malloc(N * sizeof (PerspecImg));
    for (size_t n = 0; n < N; ++n) {
        uv[n].x = (pt[n].x / resolX) + x0;
        uv[n].y = (pt[n].y / resolY) + y0;
        uv[n].r = pt[n].r;
        uv[n].g = pt[n].g;
        uv[n].b = pt[n].b;
    }
    return uv;
}
#include <string.h>
RigTrans *occludeRigTrans(RigTrans *r, size_t N, size_t *nN, double del){
    RigTrans *p = malloc(N * sizeof (RigTrans));
    memset(p, 0, N * sizeof (RigTrans));
    size_t cnt = 0;
    for (size_t n = 0; n < N; ++n) {
        size_t c = 0;
        for (; c < N; ++c) {
            if( ( r[n].x > (r[c].x - del) ) && ( r[n].x < r[c].x + del ) && ( r[n].y > (r[c].y - del) ) && ( r[n].y < r[c].y + del ) ){
                if(r[c].z > r[n].z){
                    break;
                }
            }
        }
        if(c == N){
            p[cnt].x = r[n].x;
            p[cnt].y = r[n].y;
            p[cnt].z = r[n].z;
            p[cnt].d = r[n].d;
            p[cnt].r = r[n].r;
            p[cnt].g = r[n].g;
            p[cnt].b = r[n].b;
            ++cnt;
        }
    }
    *nN = cnt;
    RigTrans *o = malloc(cnt * sizeof (RigTrans));
    for (size_t n = 0; n < cnt; ++n) {
        o[n].x = p[n].x;
        o[n].y = p[n].y;
        o[n].z = p[n].z;
        o[n].d = p[n].d;
        o[n].r = p[n].r;
        o[n].g = p[n].g;
        o[n].b = p[n].b;
    }
    free(p);
    return o;
}

#define ORIGINAL_IMG_F 1

void from3dtoPinholeProj(point3d *pt, size_t N, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, 0, 0, 0, 0, 0, 0);
    RigTrans *rp = rigTransToPerspec(r, N, ORIGINAL_IMG_F);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    prespecToImg(s, N, img, row, col);
    free(s);
}

void from3dtoUVProj(point3d *pt, size_t N, size_t resolX, size_t resolY, size_t x0, size_t y0, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, 0, 0, 0, 0, 0, 0);
    RigTrans *rp = rigTransToPerspec(r, N, ORIGINAL_IMG_F);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    PerspecImg *uv = perspecToUV(s, N, resolX, resolY, x0, y0);
    free(s);
    prespecToImg(uv, N, img, row, col);
    free(uv);
}

void from3dtoOrthProj(point3d *pt, size_t N, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, 0, 0, 0, 0, 0, 0);
    RigTrans *rp = rigTransToPerspec(r, N, 20);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    PerspecImg *uv = perspecToUV(s, N, 1, 1, 0, 0);
    free(s);
    prespecToImg(uv, N, img, row, col);
    free(uv);
}

void from3dTransformToPinholeProj(point3d *pt, size_t N, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, gama, beta, alpha, T_x, T_y, T_z);
    RigTrans *rp = rigTransToPerspec(r, N, f);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    prespecToImg(s, N, img, row, col);
    free(s);
}

void from3dTransformtoUVProj(point3d *pt, size_t N, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z, size_t resolX, size_t resolY, size_t x0, size_t y0, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, gama, beta, alpha, T_x, T_y, T_z);
    RigTrans *rp = rigTransToPerspec(r, N, f);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    PerspecImg *uv = perspecToUV(s, N, resolX, resolY, x0, y0);
    free(s);
    prespecToImg(uv, N, img, row, col);
    free(uv);
}

void from3dTransformtoOrthProj(point3d *pt, size_t N, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, gama, beta, alpha, T_x, T_y, T_z);
    RigTrans *rp = rigTransToPerspec(r, N, 20);
    free(r);
    PerspecImg *p = perspecToImage(rp, N);
    free(rp);
    PerspecImg *s = scalePerspecImage(p, N, row, col);
    free(p);
    PerspecImg *uv = perspecToUV(s, N, 1, 1, 0, 0);
    free(s);
    prespecToImg(uv, N, img, row, col);
    free(uv);
}

void from3dtoPinholeProjWOcclusion(point3d *pt, size_t N, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, 0, 0, 0, 0, 0, 0);
    size_t nN;
    RigTrans *p = occludeRigTrans(r, N, &nN, 0.002);
    free(r);
    RigTrans *rp = rigTransToPerspec(p, nN, ORIGINAL_IMG_F);
    free(p);
    PerspecImg *pp = perspecToImage(rp, nN);
    free(rp);
    PerspecImg *s = scalePerspecImage(pp, nN, row, col);
    free(pp);
    prespecToImg(s, nN, img, row, col);
    free(s);
}

void from3dTransformToPinholeProjWOcclusion(point3d *pt, size_t N, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col){
    RigTrans *r = rigTransIn3d(pt, N, gama, beta, alpha, T_x, T_y, T_z);
    size_t nN;
    RigTrans *p = occludeRigTrans(r, N, &nN, 0.002);
    free(r);
    RigTrans *rp = rigTransToPerspec(p, nN, ORIGINAL_IMG_F);
    free(p);
    PerspecImg *pp = perspecToImage(rp, nN);
    free(rp);
    PerspecImg *s = scalePerspecImage(pp, nN, row, col);
    free(pp);
    prespecToImg(s, nN, img, row, col);
    free(s);
}
