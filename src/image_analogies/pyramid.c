#include <stdlib.h>
#include <math.h>

#include "pyramid.h"


float GAUSS[5] = { 0.25-0.40/2.0, 0.25, 0.40, 0.25, 0.25-0.40/2.0 };

int*
pyramid_level_dim(int cols, int rows, int level)
{
  int* dim = malloc(sizeof(int) * 2);

  int N =(int) powf(2.0, (float)level);
  dim[0] =(int)(cols/(float)N);
  dim[1] =(int)(rows/(float)N);

  return dim;
}

Pyramid*
init_pyramid(int cols, int rows, int s)
{
  Pyramid* pyramid = malloc(sizeof(Pyramid) * L);

  for(int l=0; l<L; l++){
    int* dim = pyramid_level_dim(cols, rows, L-1-l);
    pyramid[l].cols = dim[0];
    pyramid[l].rows = dim[1];
    pyramid[l].data = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].mean = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].sd = malloc(sizeof(float) * dim[0] * dim[1]);
    if( s )
      pyramid[l].s = malloc(sizeof(int) * dim[0] * dim[1]);
    else
      pyramid[l].s = NULL;
    free(dim);
  }

  return pyramid;
}

void
free_pyramid(Pyramid* pyramid, int s)
{
  for(int l=0; l<L; l++){
    free(pyramid[l].data);
    free(pyramid[l].mean);
    free(pyramid[l].sd);
    if( s )
      free(pyramid[l].s);
  }

  free(pyramid);
}
	     
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

float
dist( int p, Pyramid* source, Pyramid* source_filter, int q, Pyramid* target, Pyramid* target_filter, int l)
{
  float luminance = 0.0, mean = 0.0, sd = 0.0, dist = 0.0;
  int q2 = q/4.0, p2 = p/4.0;
  if( l == 0 ){
    
    luminance = (source[l].data[p] - target[l].data[q]) * (source[l].data[p] - target[l].data[q]);    
    mean = (source_filter[l].mean[p] - target_filter[l].mean[q]) * (source_filter[l].mean[p] - target_filter[l].mean[q]) + (source[l].mean[p] - target[l].mean[q]) * (source[l].mean[p] - target[l].mean[q]);    
    sd = (source_filter[l].sd[p] - target_filter[l].sd[q]) * (source_filter[l].sd[p] - target_filter[l].sd[q]) + (source[l].sd[p] - target[l].sd[q]) * (source[l].sd[p] - target[l].sd[q]); 
    
  
  }else{

    luminance = (source[l].data[p] - target[l].data[q]) * (source[l].data[p] - target[l].data[q]);
    luminance += (source[l-1].data[p2] - target[l-1].data[q2]) * (source[l-1].data[p2] - target[l-1].data[q2]) + (source_filter[l-1].data[p2] - target_filter[l-1].data[q2]) * (source_filter[l-1].data[p2] - target_filter[l-1].data[q2]);
    
    mean = (source_filter[l].mean[p] - target_filter[l].mean[q]) * (source_filter[l].mean[p] - target_filter[l].mean[q]) + (source[l].mean[p] - target[l].mean[q]) * (source[l].mean[p] - target[l].mean[q]);
    mean += (source_filter[l-1].mean[p2] - target_filter[l-1].mean[q2]) * (source_filter[l-1].mean[p2] - target_filter[l-1].mean[q2]) + (source[l-1].mean[p2] - target[l-1].mean[q2]) * (source[l-1].mean[p2] - target[l-1].mean[q2]);
    
    sd = (source_filter[l].sd[p] - target_filter[l].sd[q]) * (source_filter[l].sd[p] - target_filter[l].sd[q]) + (source[l].sd[p] - target[l].sd[q]) * (source[l].sd[p] - target[l].sd[q]);
    sd += (source_filter[l-1].sd[p2] - target_filter[l-1].sd[q2]) * (source_filter[l-1].sd[p2] - target_filter[l-1].sd[q2]) + (source[l-1].sd[p2] - target[l-1].sd[q2]) * (source[l-1].sd[p2] - target[l-1].sd[q2]);    
  }
  dist = sqrtf(luminance + mean + sd);

  return dist;
}