#include <stdlib.h>

#include "conversion.h"


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
 * @brief Conversion de type du tableau de unsigned short à float
 *
 * @fn float* convertUnsignedShort2Float(unsigned short* in, int cols, int rows)
 * @param cols,rows: dimension du tableau.
 * @param in: tableau a convertir.
 * @return tableau converti en float.
 **/
float*
convertUnsignedShort2Float(unsigned short* in, int cols, int rows)
{
  float* out = malloc(sizeof(float) * cols * rows);

  for(int j=0; j<rows; j++){
    for(int i=0; i<cols; i++){
      out[j*cols+i] = (float)in[j*cols+i];
    }
  }

  return out;
}

/**
 * @brief Converti un tableau de float en unsigned short.
 *
 * @fn void convertFloat2UnsignedShort(float* in, unsigned short* out, int cols, int rows)
 * @param in tableau a convertir.
 * @param out tableau converti. Il faut donc lui donner un tableau de même dimension.
 * @param cols,rows dimension des tableaux.
 **/
void
convertFloat2UnsignedShort(float* in, unsigned short* out, int cols, int rows)
{
  float val = 0.0;
  for(int j=0; j<cols; j++){
    for(int i=0; i<rows; i++){
      val = in[i*cols+j];
      if(val > 255)
	val = 255;
      else if(val < 0)
	val = 0.0;

      out[i*cols+j] =(unsigned short)val;
    }
  }
}

/**
 * @brief Permet la multiplication d'un tableau avec une Matrice de conversion de taille 3*3.
 *
 * @fn void multiple(int cols, int rows, float A[3][3], float* B, float* output)
 * @param cols,rows dimension du tableau.
 * @param A Matrice de conversion.
 * @parma B tableau d'entré.
 * @param output tableau résultant de la multiplication entre la matrice et le tableau.
 **/
void
multiple(int cols, int rows, float A[3][3], float* B, float* output)
{
  for(int j=0; j<3 * cols; j += 3){
    for(int i=0; i<rows; i++){
      output[i*3*cols + j  ] = A[0][0] * B[i*3*cols + j] + A[0][1] * B[i*3*cols + j+1] + A[0][2] * B[i*3*cols + j+2];
      output[i*3*cols + j+1] = A[1][0] * B[i*3*cols + j] + A[1][1] * B[i*3*cols + j+1] + A[1][2] * B[i*3*cols + j+2];
      output[i*3*cols + j+2] = A[2][0] * B[i*3*cols + j] + A[2][1] * B[i*3*cols + j+1] + A[2][2] * B[i*3*cols + j+2];
    }
  }
}

/**
 * @brief Permet de convertir un tableau de RGB en tablea de YIQ. [r1g1b1, r2g2b2, ...] -> [y1q1i1, y2q2i2,...]
 *
 * @fn float* convertRGB2YIQ(int cols, int rows, float* in)
 * @param cols,rows dimensions du tableau.
 * @param in tableau RGB.
 * @return un tableau de flaot de même dimension que le tableau d'entré.
 **/
float*
convertRGB2YIQ(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, RGB2YIQ, in, out);

  return out;
}

/**
 * @brief Permet de convertir un tableau de YIQ en tablea de RGB. [y1q1i1, y2q2i2,...] -> [r1g1b1, r2g2b2, ...]
 *
 * @fn float* convertYIQ2RGB(int cols, int rows, float* in)
 * @param cols,rows dimensions du tableau.
 * @param in tableau YIQ.
 * @return un tableau de flaot de même dimension que le tableau d'entré.
 **/
float*
convertYIQ2RGB(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, YIQ2RGB, in, out);

  return out;
}

/**
 * @brief Récupère un seul canal d'un tableau RGB ou YIQ.
 *
 * @fn float* channel(int cols, int rows, float* yiq, int channel)
 * @param cols,rows dimension de l'image. (Attention pas dimension du tableau).
 * @param yiq tableau de RGB/YIQ.
 * @param channel Canal a récupérer.
 * @return un tableau de dimension de l'image contenant un canal.
 **/
float*
channel(int cols, int rows, float* yiq, int channel)
{
  float* out = malloc(sizeof(float) *cols * rows);
  int k=0;
  for(int i=0; i<rows; i++){
    for(int j=channel; j<3*cols; j+=3){
      out[k] = yiq[i*3*cols+j];
      k++;
    }
  }

  return out;
}

/**
 * @brief Permet de mettre dans un tableau trois canaux.
 *
 * @fn float* putChannel(int cols, int rows, float* y, float* I, float* q)
 * @param cols,rows dimensions d'un canal.
 * @parma y,i,q trois tableaux de même dimension représentent les trois canaux.
 * @return un tableau de type float composé de trois canaux.
 **/
float*
putChannel(int cols, int rows, float* y, float* I, float* q)
{
  float* yiq = malloc(sizeof(float)*cols*3*rows);

  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      yiq[i*3*cols+(j*3)] = y[i*cols+j];
      yiq[i*3*cols+(j*3)+1] = I[i*cols+j];
      yiq[i*3*cols+(j*3)+2] = q[i*cols+j];
    }
  }

  return yiq;
}
