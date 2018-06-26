#include <stdlib.h>
#include <math.h>

#include "pyramid.h"

/**
 * @def GAUSS
 * Tableau de taille 5 représentant un noyeau gaussien permettant de faire la convolution pour calculer une image
 * à un certain niveau de la pyramide gaussienne.
 **/

float GAUSS[5] = { 0.25-0.40/2.0, 0.25, 0.40, 0.25, 0.25-0.40/2.0 };

/**
 * @brief Initialise un tableau de taille L de structure Pyramid (allocation du tableau et de la structure).
 *        Cela représente la pyramide gaussienne où le niveau 0 représente le plus bas niveau (l'image de plus basse résolution) et L-1 l'image de base (l'image de plus grande résolution).
 *
 * @fn Pyramid* init_pyramid(int cols, int rows, int s, int L)
 * @param cols,rows dimensions de l'image.
 * @param s entier servant de booleen permettant de savoir si l'allocation du tableau s est nécéssaire ou non.
 * @return Le tableau de la pyramide gaussienne
 **/
Pyramid*
init_pyramid(int cols, int rows, int s, int L)
{
  Pyramid* pyramid = malloc(sizeof(Pyramid) * L);

  for(int l=0; l<L; l++){
    int* dim = pyramid_level_dim(cols, rows, L-1-l);
    pyramid[l].cols = dim[0];
    pyramid[l].rows = dim[1];
    pyramid[l].lum = malloc(sizeof(float) * dim[0] * dim[1]);
    if( s )
      pyramid[l].s = malloc(sizeof(int) * dim[0] * dim[1]);
    else
      pyramid[l].s = NULL;
    free(dim);
  }
  return pyramid;
}

/**
 * @brief Libère la mèmoire du tableau de la pyramide gaussienne.
 *
 * @fn void free_pyramid(Pyramid* pyramid, int s, int L)
 * @param pyramid tableau à désallouer.
 * @param s booleen permettant de savoir si oui ou non le tableau s doit être désallouer.
 * @param L taille de la pyramid.
 **/
void
free_pyramid(Pyramid* pyramid, int s, int L)
{
  for(int l=0; l<L; l++){
    free(pyramid[l].lum);
    if( s )
      free(pyramid[l].s);
  }

  free(pyramid);
}

/**
 * @brief Calcule les dimensions de l'image à un certain niveau de la pyramide.
 *
 * @fn int* pyramid_level_dim(int cols, int rows, int level)
 * @param cols,rows dimension de l'image de base.
 * @param level niveau de la pyramide ou l'on souhaite calculer les dimensions. (Attention ici 0 représente l'image de base donc faire L-level)
 * @return un tableau de deux int où tab[0] = cols et tab[1] = rows
 **/
int*
pyramid_level_dim(int cols, int rows, int level)
{
  int* dim = malloc(sizeof(int) * 2);
  if(level == 0){
    dim[0] = cols;
    dim[1] = rows;
  }else{
    int N =(int) powf(2.0, (float)level);
    dim[0] =(int)(cols/(float)N);
    dim[1] =(int)(rows/(float)N);
  }
  return dim;
}


/**
 * @brief Permet de calculer l'image au niveau level de la pyramide.
 *
 * @fn float* pyramid_level(int cols, int rows, float* data, int level)
 * @param cols,rows dimension de l'image de base.
 * @param data tableau de l'image de base.
 * @param level niveau de la pyramide ou l'on souhaite calculer les dimensions. (Attention ici 0 représente l'image de base donc faire L-level)
 * @return Un tableau de l'image au niveau level de la pyramide.
 **/
float*
pyramid_level(int cols, int rows, float* data, int level)
{
  int* dim = pyramid_level_dim(cols, rows, level);
  int cols_t = dim[0];
  int rows_t = dim[1];
  int N =(int) powf(2.0, (float)level);

  float* tmp1= malloc(sizeof(float)* rows * cols_t);
  float* out = malloc(sizeof(float)* rows_t * cols_t);

  //convolution rows
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols_t; j++){
	if(j*N-1 < 0)
	  tmp1[i*cols_t+j] = GAUSS[2] * data[i*cols + j*N] + GAUSS[3] * data[i*cols + j*N+1] + GAUSS[4] * data[i*cols + j*N+2];
	else if(j*N+1 == cols-1)
	  tmp1[i*cols_t+j] = GAUSS[0] * data[i*cols+ j*N-2] + GAUSS[1] * data[i*cols+ j*N-1] + GAUSS[2] * data[i*cols + j*N] + GAUSS[3] * data[i*cols + j*N+1];
	else
	  tmp1[i*cols_t+j] = GAUSS[0] * data[i*cols+ j*N-2] + GAUSS[1] * data[i*cols+ j*N-1] + GAUSS[2] * data[i*cols + j*N] + GAUSS[3] * data[i*cols + j*N+1] + GAUSS[4] * data[i*cols + j*N+2];
      }
    }

    //convolution cols
    for(int j=0; j<cols_t; j++){
      for(int i=0; i<rows_t; i++){
	if(i*N-1 < 0)
	  out[i*cols_t+j] = GAUSS[2] * tmp1[i*N*cols_t + j] + GAUSS[3] * tmp1[(i*N+1)*cols_t + j] + GAUSS[4] * tmp1[(i*N+2)*cols_t + j];
	else if(i*N+1 == rows-1)
	  out[i*cols_t+j] = GAUSS[0] * tmp1[(i*N-2)*cols_t + j] + GAUSS[1] * tmp1[(i*N-1)*cols_t + j] + GAUSS[2] * tmp1[i*N*cols_t + j] + GAUSS[3] * tmp1[(i*N+1)*cols_t + j];
	else
	  out[i*cols_t+j] = GAUSS[0] * tmp1[(i*N-2)*cols_t + j] + GAUSS[1] * tmp1[(i*N-1)*cols_t + j] + GAUSS[2] * tmp1[i*N*cols_t + j] + GAUSS[3] * tmp1[(i*N+1)*cols_t + j] + GAUSS[4] * tmp1[(i*N+2)*cols_t + j];
      }
    }

    free(dim);
    free(tmp1);

    return out;
}
