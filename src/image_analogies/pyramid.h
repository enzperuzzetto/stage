#ifndef __PYRAMID__
#define __PYRAMID__

/**
 * @file pyramid.h
 * @brief 
 **/

#define A 0.40

float GAUSS[5] = { 0.25-A/2.0, 0.25, A, 0.25, 0.25-A/2.0 };


pnm
pyramidGaussien(int cols, int rows, float* input, int level);// ou stocker??

#endif
