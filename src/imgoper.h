#ifndef IMGOPER_H
#define IMGOPER_H

#include <stdlib.h>

void doBinConv(const char* ipFilNam, size_t n);
void doMedianFilter(const char* ipFilNam, size_t n);
void doHistStch(const char* ipFilNam);
void doHistEq(const char* ipFilNam);
void doScharrEdgeDet(const char* ipFilNam, size_t gradTh);
void doSobelEdgeDet(const char* ipFilNam, size_t gradTh);
void doCannyEdgeDet(const char* ipFilNam);
void doKMeans(const char* ipFilNam, size_t K);
void doSpatialKMeans(const char* ipFilNam, size_t K);
void doAll(const char *ipFilNam);
void doPinholeProjection(const char *ipFilNam);
void doUVProjection(const char* ipFilNam, size_t resolX, size_t resolY, size_t x0, size_t y0);
void doOrthProjection(const char* ipFilNam);
void doTransformNPinholeProjection(const char* ipFilNam, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z);
void doTransformNUVProjection(const char* ipFilNam, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z, size_t resolX, size_t resolY, size_t x0, size_t y0);
void doTransormNOrthProjection(const char* ipFilNam, double gama, double beta, double alpha, double T_x, double T_y, double T_z);
void doOcclusion(const char* ipFilNam);
void doTransormNOcclusion(const char* ipFilNam, double gama, double beta, double alpha, double T_x, double T_y, double T_z);

#endif // IMGOPER_H
