/* Naive off reader for point cloud from .off file */
/* Only works for TP5 */
/* Point cloud is stored in points */
/* The size of points is N */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "offfil.h"


point3d *readOff(const char *fileName, int *N){
    point3d *points;

	FILE *ifp = fopen(fileName, "r");
	if (ifp == NULL) {
  		fprintf(stderr, "Can't open input file in.list!\n");
  		exit(1);
	}

	char fooChar[16];
	int N_v, fooInt;
	fscanf(ifp, "%s", fooChar);
	//printf("%s \n", fooChar);
	fscanf(ifp, "%d %d %d", &N_v, &fooInt, &fooInt);
	*N = N_v;
    points = (point3d *) malloc((unsigned)N_v * sizeof(point3d));
	//printf("%d %d %d\n", N_v, fooInt, fooInt);
	int i = 0;
	for(i = 0; i < N_v; i++){
        fscanf(ifp, "%lf %lf %lf %d %d %d %d", &(points[i].x), &(points[i].y), &(points[i].z), &(points[i].r), &(points[i].g), &(points[i].b), &fooInt);
		//printf("%f %f %f %d %d %d\n", points[i].x, points[i].y, points[i].z, points[i].r, points[i].g, points[i].b);
	}

	return points;
}
