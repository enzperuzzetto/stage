#include <stdlib.h>
#include <math.h>
#include<stdio.h>

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

//Calcul la position du  pixel au niveau l-1
int
pixelL_1(int p, int cols)
{
  int i = p/cols;
  int j = p - i*cols;
  int new_i = i/2.0;
  int new_j = j/2.0;;

  return new_i * cols/2.0 + new_j;
}

/*
float
dist( int p, Pyramid* source, Pyramid* source_filter, int q, Pyramid* target, Pyramid* target_filter, int l)
{
  float d_prime=0.0, d=0.0;
  float luminance, mean, sd, tmp;
  int k = 0, i1, j1, p1, q1;
  if( l == 0){

    tmp = source[l].data[p] - target[l].data[q];
    luminance = tmp * tmp;

    tmp = source[l].mean[p] - target[l].mean[q];
    mean = tmp * tmp;

    tmp = source[l].sd[p] - target[l].sd[q];
    sd = tmp * tmp;

    d = luminance + mean + sd;

    luminance = 0.0;
    for(int i=0; i<NEIGHBOOR_SIZE_FINER; i++){
      for(int j=0; j<NEIGHBOOR_SIZE_FINER; j++){
	i1 = i-NEIGHBOOR_SIZE_FINER/2.0;
	j1 = j-NEIGHBOOR_SIZE_FINER/2.0;

	p1 = p + (i1*source_filter[l].cols+j1);
	q1 = q + (i1*target_filter[l].cols+j1);

	if(p1 < 0 || q1 < 0)
	  break;
	if(p1 == p || q1 == q)
	  break;

	tmp = (source_filter[l].data[p1] - target_filter[l].data[q1]);
	luminance += tmp * tmp;
	k++;
      }
      if(p1 == p || q1 == q)
	break;
    }

    if( k == 0)
      luminance = 99999999999999999;
    
    tmp = source_filter[l].mean[p] - target_filter[l].mean[q];
    mean = tmp * tmp;

    tmp = source_filter[l].sd[p] - target_filter[l].sd[q];
    sd = tmp * tmp;

    d_prime = luminance + mean + sd;

  }else{

    int q2 = 0.0, p2 = 0.0;
    q2 = pixelL_1(q, target[l].cols);
    p2 = pixelL_1(p, source[l].cols);
    
    tmp = source[l].data[p] - target[l].data[q];
    luminance = tmp * tmp;
    tmp = source[l-1].data[p2] - target[l-1].data[q2];
    luminance += tmp * tmp;

    tmp = source[l].mean[p] - target[l].mean[q];
    mean = tmp * tmp;
    tmp = source[l-1].mean[p2] - target[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = source[l].sd[p] - target[l].sd[q];
    sd = tmp * tmp;
    tmp = source[l-1].sd[p2] - target[l-1].sd[q2];
    sd += tmp * tmp;

    d = luminance + mean + sd;

    luminance = 0.0;
    for(int i=0; i<NEIGHBOOR_SIZE_FINER; i++){
      for(int j=0; j<NEIGHBOOR_SIZE_FINER; j++){
	i1 = i-NEIGHBOOR_SIZE_FINER/2.0;
	j1 = j-NEIGHBOOR_SIZE_FINER/2.0;

	p1 = p + (i1*source_filter[l].cols+j1);
	q1 = q + (i1*target_filter[l].cols+j1);

	if(p1 < 0 || q1 < 0)
	  break;
	if(p1 == p || q1 == q)
	  break;

	tmp = (source_filter[l].data[p1] - target_filter[l].data[q1]);
	luminance += tmp * tmp;
	k++;
      }
      if(p1 == p || q1 == q)
	break;
    }

    if (k == 0)
      luminance = 99999999999999999;
    k=0;

    for(int i=0; i<NEIGHBOOR_SIZE_COARSER; i++){
      for(int j=0; j<NEIGHBOOR_SIZE_COARSER; j++){
	i1 = i-NEIGHBOOR_SIZE_COARSER/2.0;
	j1 = j-NEIGHBOOR_SIZE_COARSER/2.0;

	p1 = p2 + (i1*source_filter[l-1].cols+j1);
	q1 = q2 + (i1*target_filter[l-1].cols+j1);

	if(p1 < 0 || q1 < 0)
	  break;
	
	tmp = (source_filter[l-1].data[p1] - target_filter[l-1].data[q1]);
	luminance += tmp * tmp;
	k++;
      }
    }

    if (k == 0)
      luminance = 99999999999999999;

    tmp = source_filter[l].mean[p] - target_filter[l].mean[q];
    mean = tmp * tmp;
    tmp = source_filter[l-1].mean[p2] - target_filter[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = source_filter[l].sd[p] - target_filter[l].sd[q];
    sd = tmp * tmp;
    tmp = source_filter[l-1].sd[p2] - target[l-1].sd[q2];
    sd += tmp * tmp;

    d_prime = luminance + mean + sd;

  }
  
  return  W * d + d_prime;
}
*/

float
dist( int xs, int ys, Pyramid* source, Pyramid* source_filter, int xt, int yt, Pyramid* target, Pyramid* target_filter, int l)
{
  int begin_t, end_t, istart, iend, jstart, jend;
  int s = (int)((NEIGHBOOR_SIZE_FINER-1)/2.0);
  int cols_s = source_filter[l].cols;
  int cols_t = target_filter[l].cols;
  int rows_t = target_filter[l].rows;
  int p = xs * cols_s + ys, q = xt * cols_t + yt, p1, q1 ,p2 = (xs/2.0) * (cols_s/2.0) + (ys/2.0), q2 = (xt/2.0) * (cols_t/2.0) + (yt/2.0), k = 0;
  double lum, mean, sd, tmp, d, d_prime;

  // A et B
  tmp = source[l].data[p] - target[l].data[q];
  lum = tmp * tmp;

  tmp = source[l].mean[p] - target[l].mean[q];
  mean = tmp * tmp;

  tmp = source[l].sd[p] - target[l].sd[q];
  sd = tmp * tmp;

  if(l > 0){
    
    tmp = source[l-1].data[p2] - target[l-1].data[q2];
    lum += tmp * tmp;

    tmp = source[l-1].mean[p2] - target[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = source[l-1].sd[p2] - target[l-1].sd[q2];
    sd += tmp * tmp;

  }

  d = lum + mean + sd;

  //A' et B'

  lum = 0.0;
  mean = 0.0;
  sd = 0.0;
  
  if(xt < s)
    begin_t = -s + (s - xt);
  else
    begin_t = -s;
  if(xt > rows_t-1 - s)
    end_t = s - (xt - (rows_t-1 - s));
  else
    end_t = s;

  istart = begin_t;
  iend = end_t;

  if(yt < s)
    begin_t = -s + (s - yt);
  else
    begin_t = -s;
  if(yt > cols_t-1 - s)
    end_t = s -(yt - (cols_t-1 - s));
  else
    end_t = s;

  jstart = begin_t;
  jend = end_t;
  
  for(int i=istart; i<iend+1; i++){
    for(int j=jstart; j<jend+1; j++){

      p1 = p + (i*cols_s+j);
      q1 = q + (i*cols_t+j);
           
      if(p1 == p || q1 == q)
	break;

      tmp = source_filter[l].data[p1] - target_filter[l].data[q1];
      lum += tmp * tmp;
      k++;
    }
    if(p1 == p || q1 == q)
	break;
  }

  if( k == 0)
    lum = 99999999999999999;

  tmp = source_filter[l].mean[p] - target_filter[l].mean[q];
  mean = tmp * tmp;

  tmp = source_filter[l].sd[p] - target_filter[l].sd[q];
  sd = tmp * tmp;

  if( l > 0){

    xs /= 2.0;
    ys /= 2.0;
    xt /= 2.0;
    yt /= 2.0;
    s /= 2.0;
    cols_s /= 2.0;
    cols_t /= 2.0;
    rows_t /= 2.0;
    k = 0;

    if(xt < s)
      begin_t = -s + (s - xt);
    else
      begin_t = -s;
    if(xt > rows_t-1 - s)
      end_t = s - (xt - (rows_t-1 - s));
    else
      end_t = s;

    istart = begin_t;
    iend = end_t;

    if(yt < s)
      begin_t = -s + (s - yt);
    else
      begin_t = -s;
    if(yt > cols_t-1 - s)
      end_t = s -(yt - (cols_t-1 - s));
    else
      end_t = s;

    jstart = begin_t;
    jend = end_t;
  
    for(int i=istart; i<iend+1; i++){
      for(int j=jstart; j<jend+1; j++){

	p1 = p2 + (i*cols_s+j);
	q1 = q2 + (i*cols_t+j);

	tmp = source_filter[l-1].data[p1] - target_filter[l-1].data[q1];
	lum += tmp * tmp;
	k++;
      }
    }

    tmp = source_filter[l-1].mean[p2] - target_filter[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = source_filter[l-1].sd[p2] - target_filter[l-1].sd[q2];
    sd += tmp * tmp;

  }

  d_prime = lum + mean + sd;

  return  W * d + d_prime;
}
