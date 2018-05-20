#ifndef __CONVERSION_H__

#define __CONVERSION_H__

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
 * @brief
 *
 * @Param
 **/
float*
convertUnsignedShort2Float(unsigned short* in, int cols, int rows);

/**
 * @brief
 *
 * @Param
 **/
void
convertFloat2UnsignedShort(float* in, unsigned short* out, int cols, int rows);

/**
 * @brief
 *
 * @Param
 **/
void
multiple(int cols, int rows, float A[3][3], float* B, float* output);

/**
 * @brief
 *
 * @Param
 **/
float*
convertRGB2YIQ(int cols, int rows, float* in);

/**
 * @brief
 *
 * @Param
 **/
float*
convertYIQ2RGB(int cols, int rows, float* in);

/**
 * @brief
 *
 * @Param
 **/
float*
channel(int cols, int rows, float* yiq, int channel);

/**
 * @brief
 *
 * @Param
 **/
float*
putChannel(int cols, int rows, float* y, float* I, float* q);

#endif
