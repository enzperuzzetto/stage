#ifndef __FEATURES__
#define __FEATURES__

/**
 * @file features.h
 * @brief Header of features.c that implements functions
 *        for convert color channel.
 **/


/**
 * @brief Matrix to convert channel 
 **/

float RGB2YIQ[3][3] = {
  { 0.299, 0.587, 0.114 },
  { 0.595716, -0.274453, -0.321263 },
  { 0.211456, -0.522591, 0.311135 }
};

float YIQ2RGB[3][3] = {
  { 1.0, 0.9563, 0.6210 },
  { 1.0, -0.2721, -0.6474 },
  { 1.0, -1.1070, 1.7046 }
};


/**
 * @brief compute the multiplication of matrix
 *
 **/

void
multiple(int cols, int rows, float M[3][3], float* B, float* output);


void
convertRGB2YIQ(int cols, int rows, float* input, float* output);

void
convertYIQ2RGB(int cols, int rows, float* input, float* output);

#endif
