#include "stats.h"

/**
 * @brief Calcule la moyenne des pixels de l'image.
 *
 * @fn float mean(int cols, int cols, float* in)
 * @param cols, rows: dimensions de l'image.
 * @param in: tableau de la luminance de l'image.
 * @return float: moyenne de l'image.
 **/
float
mean(int cols, int rows, float* in)
{
  float m = 0.0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      m += in[i*cols+j];
    }
  }

  return m/(float)(rows*cols);
}

/**
 * @brief Calcule l'écart-type des pixels de l'image.
 *
 * @fn float sd(int cols, int rows, float* in, float mean)
 * @param cols,rows dimensions de l'image.
 * @param in tableau de la luminance de l'image.
 * @param mean moyenne de l'image.
 * @return float ecart-type de l'image.
 **/
float
sd(int cols, int rows, float* in, float mean)
{
  float var = 0.0, tmp;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      tmp = (in[i*cols+j]-mean) * (in[i*cols+j]-mean);
      var += tmp;
    }
  }

  return sqrtf(var/(float)(rows*cols));
}

/**
 * @brief retourne la plus grande valeur.
 *
 * @fn int max(int a, int b)
 * @param a,b valeur a comparer
 * @return valeur la plus grande
 **/
int
max(int a, int b)
{
  return (a<b)?b:a;
}

/**
 * @brief Permet de mettre la luminance de l'image A dans le ton de la luminance de l'image B
 *
 * @fn float* luminanceRemapping(int Aw, int Ah, float* A, int Bw, int Bh, float* B)
 * @param Aw,Ah dimensions de l'image A.
 * @param A luminance de l'image A.
 * @param Bw,Bh dimensions de l'image B.
 * @param B luminance de l'image B.
 * @return Un tableau de la luminance de l'image A remapping
 **/
float*
luminanceRemapping(int Aw, int Ah, float* A, int Bw, int Bh, float* B)
{
  float sdA, sdB, meanA, meanB, quot;
  meanA = mean(Aw, Ah, A);
  sdA = sd(Aw, Ah, A, meanA);
  meanB = mean(Bw, Bh, B);
  sdB = sd(Bw, Bh, B, meanB);

  quot = sdB / sdA;

  float* out = malloc(sizeof(float) * Aw * Ah);

  int pixel;
  for(int i=0; i<Ah; i++){
    for(int j=0; j<Aw; j++){
      pixel = i*Aw+j;
      out[pixel] = quot * (A[pixel] - meanA) + meanB;
    }
  }

  return out;
}
