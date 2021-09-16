#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define UNUSED(...) (void)(__VA_ARGS__)

//#define FILTER_SIZ 5
#define MAX_FIL_NAM_SIZ 64

//typedef unsigned char pix;
//typedef char filter;
//typedef short pixprod;
typedef int mat;

typedef unsigned char bit;
typedef unsigned char pix;

typedef struct{
   double  x; // If you want to access x, you need to do p.x, suppose you decleared p as "struct point3d p;"
   double  y; // Similar to above
   double  z;
   int r;
   int g;
   int b;
}point3d;

#define ASSERT_MSG(cond, msg) assert((cond) && (msg))

#define NULL_PTR_CHK(x) ASSERT_MSG(x != NULL, "null ptr chk failed")
#define SAME_PTR_CHK(x, y) ASSERT_MSG(x != y, "same ptr passed for ip and op")
#define ABOVE_ZERO_CHK(x) ASSERT_MSG(x > 0.0, "ip not above 0")
// this check is added since we want to generalize the library for colour images while supporting grayscale still
#define GRAY_IMG_CHK(nbCh) ASSERT_MSG(nbCh == 1, "not a grayscale image")

#define GRAYSCALE_DIM (pix)1

#define MAX_PIX_VAL ( (mat) ( pow( 2, 8 * sizeof(pix) ) - 1 ) )

#endif // UTILS_H
