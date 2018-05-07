#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bcl.h>

#define A 0.40

float GAUSS[5] = { 0.25-A/2.0, 0.25, A, 0.25, 0.25-A/2.0 };

void
process(char *ims_name, int level)
{
  pnm ims = pnm_load(ims_name);
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);

  int N =(int) powf(2.0, (float)level);
  int cols_t = cols/(float)N;
  int rows_t = rows/(float)N;

  pnm imd = pnm_new(cols_t, rows_t, PnmRawPpm);

  unsigned short* output = malloc(sizeof(unsigned short) * cols_t * rows_t);
  float* tmp = malloc(sizeof(float)* rows * cols);
  float* tmp1= malloc(sizeof(float)* rows * cols_t);
  float* out = malloc(sizeof(float)* rows_t * cols_t);
  
  for(int k=0; k<3; k++){
    unsigned short* in = pnm_get_channel(ims, NULL, k);

    for(int j=0; j<rows; j++){
      for(int i=0; i<cols; i++){
	tmp[j*cols+i] = (float)in[j*cols+i];
      }
    }

 
    //convolution rows
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols_t; j++){
	if(j*N-1 < 0)
	  tmp1[i*cols_t+j] = GAUSS[2] * tmp[i*cols + j*N] + GAUSS[3] * tmp[i*cols + j*N+1] + GAUSS[4] * tmp[i*cols + j*N+2];
	else if(j*N+1 == cols-1)
	  tmp1[i*cols_t+j] = GAUSS[0] * tmp[i*cols+ j*N-2] + GAUSS[1] * tmp[i*cols+ j*N-1] + GAUSS[2] * tmp[i*cols + j*N] + GAUSS[3] * tmp[i*cols + j*N+1];
	else
	  tmp1[i*cols_t+j] = GAUSS[0] * tmp[i*cols+ j*N-2] + GAUSS[1] * tmp[i*cols+ j*N-1] + GAUSS[2] * tmp[i*cols + j*N] + GAUSS[3] * tmp[i*cols + j*N+1] + GAUSS[4] * tmp[i*cols + j*N+2];
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

    

    float val;
    for(int j=0; j<cols_t; j++){
      for(int i=0; i<rows_t; i++){
	val = out[i*cols_t+j];
	if(val > 255)
	  val = 255;
	else if(val < 0)
	  val = 0.0;
      
	output[i*cols_t+j] =(unsigned short)val;
      }
    }

    pnm_set_channel(imd, output, k);

    free(in);
  }

  pnm_save(imd, PnmRawPpm, "reduce.ppm");  

  free(output);
  free(out);
  free(tmp1);
  free(tmp);
  pnm_free(ims);
  pnm_free(imd);
}
  



void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <level> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 2
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1], atoi(argv[2]));
  return EXIT_SUCCESS;
}
