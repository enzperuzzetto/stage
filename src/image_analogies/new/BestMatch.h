#ifndef __BESTMATCH_H__

#define __BESTMATCH_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pyramid.h"


float**
initKernel(int size);

void
freeKernel(float** kernel, int size);

float
dist(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel);

int
Ashikhmin(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel);

int
BestMatch(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float w, float levelweight, int L, int K, float** kernel);

#endif
