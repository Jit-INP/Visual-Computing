#ifndef PROJ_H
#define PROJ_H

#include "offfil.h"
#include "matrix.h"

void from3dtoPinholeProj(point3d *pt, size_t N, mat *img, size_t row, size_t col);
void from3dtoUVProj(point3d *pt, size_t N, size_t resolX, size_t resolY, size_t x0, size_t y0, mat *img, size_t row, size_t col);
void from3dtoOrthProj(point3d *pt, size_t N, mat *img, size_t row, size_t col);
void from3dTransformToPinholeProj(point3d *pt, size_t N, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col);
void from3dTransformtoUVProj(point3d *pt, size_t N, double f, double gama, double beta, double alpha, double T_x, double T_y, double T_z, size_t resolX, size_t resolY, size_t x0, size_t y0, mat *img, size_t row, size_t col);
void from3dTransformtoOrthProj(point3d *pt, size_t N, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col);
void from3dtoPinholeProjWOcclusion(point3d *pt, size_t N, mat *img, size_t row, size_t col);
void from3dTransformToPinholeProjWOcclusion(point3d *pt, size_t N, double gama, double beta, double alpha, double T_x, double T_y, double T_z, mat *img, size_t row, size_t col);

#endif // PROJ_H
