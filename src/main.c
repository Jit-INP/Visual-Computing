#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imgoper.h"
#include "utils.h"

#include "convol.h"
#include "filter.h"
#include "hist.h"
#include "matrix.h"
#include "median.h"
#include "kmeans.h"
#include "proj.h"


typedef enum{
    eBinConv = 0,
    eMedFilter,
    eHistStch,
    eHistEq,
    eScharrEdge,
    eCannyEdge,
    eKMeans,
    eSpatialKMeans,
    eNbImgOper
}ImgOper;

int main(int argc, char **argv){
    UNUSED(argc, argv);
//    doPinholeProjection("frameCube.off");
//    doPinholeProjection("human.off");
//    doUVProjection("frameCube.off", 2, 1, 50, 50);
//    doUVProjection("human.off", 1, 1, 0, 0);
//    doOrthProjection("frameCube.off");
//    doTransformNUVProjection("frameCube.off", 0.5, 0.5236, 0.5236, 0, 0, 0, 0, 2, 1, 50, 50);
//    doTransformNPinholeProjection("frameCube.off", 0.5, 0.5236, 0.5236, 0, 0, 0, 0);
//    doTransormNOrthProjection("frameCube.off", 0.262, 0.262, 0.262, 0, 0, 0);
//    doOcclusion("human.off");
//    doTransormNOcclusion("human.off", 3.14, 0, 0, 0, 0, 0);
//    doTransformNPinholeProjection("human.off", 1, 3.14, 0, 0, 0, 0, 0);
    doPinholeProjection("frameCube.off");
    doPinholeProjection("human.off");
    doTransformNPinholeProjection("frameCube.off", 1, 0.5236, 0.5236, 0, 0, 0, 0);
    doUVProjection("frameCube.off", 1, 2, 50, 50);
    doTransformNUVProjection("frameCube.off", 0.5, 0.5236, 0.5236, 0, 0, 0, 0, 2, 3, 50, 50);
    doOrthProjection("frameCube.off");
    doTransormNOrthProjection("frameCube.off", 0.262, 0.262, 0.262, 0, 0, 0);
    //doOcclusion("human.off");
    //doTransormNOcclusion("human.off", 0, 0, -1.57, 0, 0, 0);
    return 0;
}
