#ifndef __FEATURES_H__

#define __FEATURES_H__

/**
 * @brief
 *
 * @Param
 **/
float
mean(int cols, int rows, float* in, int index, int size, int half);


/**
 * @brief
 *
 * @Param
 **/
float
variance(int cols, int rows, float* in, float mean, int index, int size, int half);


/**
 * @brief
 *
 * @Param
 **/
float
standard_deviation(int cols, int rows, float* in, float mean, int index, int size, int half);


/**
 * @brief
 *
 * @Param
 **/
void
mean_img(int cols, int rows, float* in, float* means, int size, int half);


/**
 * @brief
 *
 * @Param
 **/
void
sd_img(int cols, int rows, float* in, float* mean, float* sd, int size, int half);


/**
 * @brief
 *
 * @Param
 **/
float
dist( int p, Pyramid* source, Pyramid* source_filter, int q, Pyramid* target, Pyramid* target_filter, int l);

/**
 * @brief
 *
 * @Param
 **/
int
max(int a, int b);

#endif
