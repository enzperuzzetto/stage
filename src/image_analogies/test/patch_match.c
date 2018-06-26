#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <bcl.h>

#define PATCH 5
#define NB_MAX_ITER 5
#define ALPHA 0.5

int
max(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}

float
dist( int p, float* A, int q, float* B)
{
  float r = (A[p] - B[q]) * (A[p] - B[q]);
  return sqrt(r);
}

float
Dist(int p, int cols_s, int rows_s, float* A, int q, int cols_t, int rows_t, float* B)
{
  float d = 0.0;
  int i_p1, j_p1, p_patch, q_patch;
  for(int i=0; i<PATCH; i++){
    for(int j=0; j<PATCH; j++){
      i_p1 = i - PATCH/2.0;
      j_p1 = j - PATCH/2.0;

      p_patch = p + i_p1*cols_s +j_p1;
      q_patch = q + i_p1*cols_t + j_p1;

      if(p_patch < 0)
	p_patch = 0;
      else if( p_patch > cols_s*rows_s -1)
	p_patch = cols_s*rows_s -1;

      if(q_patch < 0)
	q_patch = 0;
      else if( q_patch > cols_t*rows_t -1)
	q_patch = cols_t*rows_t -1;

      d += dist(p_patch, A, q_patch, B);
    }
  }

  return d;
}

int*
initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
{
  int* nnf = malloc(sizeof(int) * cols_s * rows_s);
  srand(time(NULL));
  for(int k=0; k<cols_s*rows_s; k++){
    nnf[k] = rand()%(cols_t*rows_t);
  }

  return nnf;
}

int
propagation(int cols_s, int rows_s, float* A, int cols_t, int rows_t, float* B, int* nnf, int p)
{
  int pixel_min = p;
  float dist_min = Dist(p, cols_s, rows_s, A, nnf[p], cols_t, rows_t, B);
  //left
  int p1 = p-1;

  if(p1 <0)
    p1=0;
  else if(p1 > cols_s*rows_s-1)
    p1 = cols_s*rows_s-1;

  float d = Dist(p1, cols_s, rows_s, A, nnf[p1], cols_t, rows_t, B);

  if( d < dist_min ){
    dist_min = d;
    pixel_min = p1;
  }

  //top
  p1 = p - cols_s;
  //###
  if(p1 <0)
    p1=0;
  else if(p1 > cols_s*rows_s-1)
    p1 = cols_s*rows_s-1;

  d = Dist(p1, cols_s, rows_s, A, nnf[p1], cols_t, rows_t, B);

  if( d < dist_min ){
    dist_min = d;
    pixel_min = p1;
  }

  return pixel_min;
}

int
search(int cols_s, int rows_s, float* A, int cols_t, int rows_t, float* B, int* nnf, int p)
{
  int w = max(cols_t, rows_t), pixel_x = 0, pixel_y = 0, pixel = 0, pixel_min = nnf[p], i=0;
  float dist_min = Dist(p, cols_s, rows_s, A, pixel_min, cols_t, rows_t, B), d = 0.0, alpha = powf(ALPHA, i), x = 0.0, y = 0.0;

  //for(int i=0; i< NB_MAX_ITER; i++){
  while( w*alpha > 1){

    x = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;
    y = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;

    /*
    if(w*alpha < 1)
      break;
    */
    pixel_x = w*alpha*x;
    pixel_y = w*alpha*y;
    pixel = pixel_min + (pixel_x*cols_t + pixel_y);


    if(pixel <0)
      pixel=0;
    else if(pixel > cols_t*rows_t-1)
      pixel = cols_t*rows_t-1;

    d = Dist(p, cols_s, rows_s, A, pixel, cols_t, rows_t, B);

    if(d < dist_min){
      dist_min = d;
      pixel_min = pixel;
    }

    i++;
    alpha = powf(ALPHA, i);
  }

  return pixel_min;
}

int
iteration(int cols_s, int rows_s, float* A, int cols_t, int rows_t, float* B, int p, int* nnf)
{
  int pixel = propagation(cols_s, rows_s, A, cols_t, rows_t, B, nnf, p);
  pixel = search(cols_s, rows_s, A, cols_t, rows_t, B, nnf, pixel);

  return pixel;
}

void
process(char *ims_name, char *imt_name )
{
  pnm ims = pnm_load(ims_name);
  int cols_s = pnm_get_width(ims);
  int rows_s = pnm_get_height(ims);

  pnm imt = pnm_load(imt_name);
  int cols_t = pnm_get_width(imt);
  int rows_t = pnm_get_height(imt);


  pnm imd = pnm_new(cols_s, rows_s, PnmRawPpm);


  unsigned short* data_s = pnm_get_channel(ims, NULL, 0);
  unsigned short* data_t = pnm_get_channel(imt, NULL, 0);
  unsigned short* data_out_r = malloc(sizeof(unsigned short) * cols_s * rows_s);

  float* A_r = malloc(sizeof(float) * cols_s * rows_s);
  float* B_r = malloc(sizeof(float) * cols_t * rows_t);
  float* out_r = malloc(sizeof(float) * cols_s * rows_s);

  for(int i=0; i< rows_s; i++){
    for(int j=0; j<cols_s; j++){
      A_r[i*cols_s+j] = (float)data_s[i*cols_s+j];
    }
  }


  for(int i=0; i< rows_t; i++){
    for(int j=0; j<cols_t; j++){
      B_r[i*cols_t+j] = (float)data_t[i*cols_t+j];

    }
  }

  int* nnf = initialisation(cols_s, rows_s, cols_t, rows_t);

  for(int i=0; i<NB_MAX_ITER; i++){
    //initialisation
    int p, q;
    if(i%2 == 0){//top to bot
      for(int x=0; x<rows_s; x++){
	for(int y=0; y<cols_s; y++){
	  p = x*cols_s+y;
	  q = iteration(cols_s, rows_s, A_r, cols_t, rows_t, B_r, p, nnf);
	  if( i > 0){
	    if( dist(p, A_r, q, B_r) < dist(p, A_r, p, out_r)){
	      out_r[p] = B_r[q];
	      nnf[p] = q;
	    }


	  }
	  else{
	    out_r[p] = B_r[q];
	    nnf[p] = q;
	  }
	}
      }
    }else{//bot to top
      for(int x=rows_s-1; x>=0; x--){
	for(int y=cols_s-1; y>=0; y--){
	  p = x*cols_s+y;
	  q = iteration(cols_s, rows_s, A_r, cols_t, rows_t, B_r, p, nnf);
	  if( i > 0){
	    if( dist(p, A_r, q, B_r) < dist(p, A_r, p, out_r)){
	      out_r[p] = B_r[q];
	      nnf[p] = q;
	    }

	  }
	  else{
	    out_r[p] = B_r[q];
	    nnf[p] = q;
	  }
	}
      }
    }
  }

  free(nnf);
  float val;
  for(int j=0; j<cols_s; j++){
    for(int i=0; i<rows_s; i++){
      val = out_r[i*cols_t+j];
      if(val > 255)
	val = 255;
      else if(val < 0)
	val = 0.0;

      data_out_r[i*cols_t+j] =(unsigned short)val;
    }
  }

  pnm_set_channel(imd, data_out_r, 0);

  pnm_save(imd, PnmRawPpm, "PatchMatch.ppm");

  free(data_s);
  free(data_t);
  free(data_out_r);
  free(A_r);
  free(B_r);
  free(out_r);
  pnm_free(ims);
  pnm_free(imt);
  pnm_free(imd);
}



void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <imt> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 2
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1], argv[2]);
  return EXIT_SUCCESS;
}
