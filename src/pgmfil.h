#ifndef PGMFIL_H
#define PGMFIL_H

#include <stdlib.h>
#include "utils.h"

typedef enum{
    eP2 = 0,
    eP3,
    eP5,
    eP6,
    eNbPGMTypes
} PGMType;

void pgmToImgMat(const char *filNam, mat **imgMat, size_t *row, size_t *col, size_t *ch, size_t *maxVal, PGMType *type);
void imgMatToPgm(const char *filNam, mat *imgMat, size_t row, size_t col, size_t ch, size_t maxVal, PGMType type);

#endif // PGMFIL_H
