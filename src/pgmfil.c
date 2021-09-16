#include "pgmfil.h"
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "matrix.h"

// to do:
// current max pix value -> max possible pix value
// max pix value, one that's read from file
// filter size not hardcoded
// clean up menu
// option to enable or disable intermediates
// check why name not appended with ./dir/img.pgm or is it always
// check processing with p2


int fileGetC(FILE *fileH){
    int ich = getc( fileH );
    assert(ich != EOF );
    return ich;
}

// get char skipping comments
char pm_getc(FILE* file){
    int ich;
    char ch;

    ich = fileGetC(file);
    ch = (char) ich;

    if ( ch == '#' ){
        do{
            ich = fileGetC(file);
            ch = (char) ich;
        } while ( ch != '\n' && ch != '\r' ); // skip the comment and yield newline char
    }
    return ch;
}

// p1 -> p4
bit pm_getbit(FILE* file){
    char ch;

    do {
        ch = pm_getc( file );
    } while ( ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' ); // trying to find 1st non-whitespace char

    switch (ch) {
    case '0':
        return 0;
    case '1':
        return 1;
    default:
        assert(0);
    }
}


unsigned char pm_getrawbyte(FILE* file){
    return  (unsigned char) fileGetC(file);
}

int pm_getint( FILE* file){
    char ch;
    do {
        ch = pm_getc( file );
    } while ( ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' );// get 1st non-whitespace char

    if ( ch < '0' || ch > '9' ){ // check if it's valid ascii dec
        assert(0);
    }

    int i = 0;
    do {
        i = i * 10 + ch - '0';
        ch = pm_getc( file );
    } while ( ch >= '0' && ch <= '9' ); // until next whitespace convert chars in ascii decimal to numbers

    return i;
}

pix *allocImgMem(size_t r, size_t c, size_t nbCh){
    ABOVE_ZERO_CHK(r);
    ABOVE_ZERO_CHK(c);
    ABOVE_ZERO_CHK(nbCh);
    return malloc(r * c * nbCh * sizeof (pix));
}

void imgMatToPixMap(pix *pMap, mat *matr, size_t row, size_t col, size_t ch){
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            for (size_t k = 0; k < ch; ++k) {
                size_t idx = (i * col * ch) + (j * ch) + k;
                pix pVal;
                mat mVal = matr[idx];
                if(mVal < 0){
                    pVal = 0;
                } else {
                    if(mVal > MAX_PIX_VAL){
                        pVal = (pix) MAX_PIX_VAL;
                    } else {
                        pVal = (pix) mVal;
                    }
                }
                pMap[idx] = pVal;
            }
        }
    }
}

void pixMapToImgMat(mat *matr, pix *pMap, size_t row, size_t col, size_t ch){
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            for (size_t k = 0; k < ch; ++k) {
                matr[(i * col * ch) + (j * ch) + k] = (mat) pMap[(i * col * ch) + (j * ch) + k];
            }
        }
    }
}

void pgmToImgMat(const char *filNam, mat **imgMat, size_t *row, size_t *col, size_t *ch, size_t *maxVal, PGMType *type){
    NULL_PTR_CHK(filNam);
    NULL_PTR_CHK(imgMat);
    NULL_PTR_CHK(row);
    NULL_PTR_CHK(col);
    NULL_PTR_CHK(type);

    FILE* ifp = fopen(filNam, "r");
    NULL_PTR_CHK(ifp);

    fileGetC( ifp );
    int ich2 = fileGetC( ifp );

    switch (ich2) {
    case '2':
        *type = eP2;
        *ch = 1;
        break;
    case '3':
        *type = eP3;
        *ch = 3;
        break;
    case '5':
        *type = eP5;
        *ch = 1;
        break;
    case '6':
        *type = eP6;
        *ch = 3;
        break;
    default:
        assert(0);
        break;
    }

    *col = (unsigned)pm_getint( ifp );
    *row = (unsigned)pm_getint( ifp );
    *maxVal = (unsigned)pm_getint( ifp );

    pix *pixMap = allocImgMem(*row, *col, *ch);

    for(size_t i=0; i < *row; ++i){
        for(size_t j=0; j < *col; ++j){
            for (size_t k=0; k < *ch; ++k) {
                pix val;
                if((*type == eP5 )||(*type == eP6 )){
                    val = pm_getrawbyte(ifp);
                } else {
                    val = (pix)pm_getint(ifp) - '0';
                }
                pixMap[(i * *col * *ch) + (j * *ch) + k] = val;
            }
        }
    }
    fclose(ifp);

    *imgMat = allocMatMem(*row, *col, *ch);
    pixMapToImgMat(*imgMat, pixMap, *row, *col, *ch);
    free(pixMap);
}

void imgMatToPgm(const char *filNam, mat *imgMat, size_t row, size_t col, size_t ch, size_t maxVal, PGMType type){
    FILE* ifp = fopen(filNam, "w");
    NULL_PTR_CHK(ifp);

    char *imgType = strdup("P2");
    switch(type){
    case eP2:
        strcpy(imgType, "P2");
        break;
    case eP3:
        strcpy(imgType, "P3");
        break;
    case eP5:
        strcpy(imgType, "P5");
        break;
    case eP6:
        strcpy(imgType, "P6");
        break;
    default:
        assert(0);
        break;
    }
    fprintf(ifp, "%s\n", imgType);
    free(imgType);
    fprintf(ifp, "%d %d\n", (signed)col, (signed)row);
    fprintf(ifp, "%d\n", (signed)maxVal);

    pix *pixMap = allocImgMem(row, col, ch);
    imgMatToPixMap(pixMap, imgMat, row, col, ch);
    for(size_t i=0; i < row; ++i){
        for(size_t j=0; j < col; ++j){
            for (size_t k=0; k < ch; ++k) {
                pix val = pixMap[(i * col * ch) + (j * ch) + k];
                if((type == eP5 )||(type == eP6 )){
                    fprintf(ifp, "%c", val);
                } else {
                    fprintf(ifp, "%d ", val + '0');
                }
            }
        }
    }
    free(pixMap);
    fclose(ifp);
}
