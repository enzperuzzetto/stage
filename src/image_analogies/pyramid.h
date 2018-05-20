#ifndef __PYRAMID_H__

#define __PYRAMID_H__


typedef struct Pyramid
{
  int cols;
  int rows;
  float* data;
  float* mean;
  float* sd;
  int* s;
}Pyramid;

/**
 * @brief
 *
 * @Param
 **/
int*
pyramid_level_dim(int cols, int rows, int level);

/**
 * @brief
 *
 * @Param
 **/
Pyramid*
init_pyramid(int cols, int rows, int s);

/**
 * @brief
 *
 * @Param
 **/
void
free_pyramid(Pyramid* pyramid, int s);

/**
 * @brief
 *
 * @Param
 **/
float*
pyramid_level(int cols, int rows, float* data, int level);

#endif
