#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bcl.h>

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

void
process(char *ims_name)
{
  pnm ims = pnm_load(ims_name);
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);

  unsigned short* in = pnm_get_image(ims);
  
  pnm imd = pnm_new(cols, rows, PnmRawPpm);
  

  unsigned short* output = malloc(sizeof(unsigned short) * 3 * cols * rows);
  float* tmp = malloc(sizeof(float)* rows * 3 * cols);
  float* out = malloc(sizeof(float)* rows * 3 * cols);


  for(int j=0; j<3*cols; j++){
    for(int i=0; i<rows; i++){
      tmp[i*3*cols+j] = (float)in[i*3*cols+j];
    }
  }

  multiple(cols, rows, RGB2YIQ, tmp, out);

  multiple(cols, rows, YIQ2RGB, out, tmp);

  float value;
  for(int j=0; j<3*cols; j++){
    for(int i=0; i<rows; i++){
      value = tmp[i*3*cols+j];
      if( value > 255)
	value = 255;
      else if(value < 0)
	value = 0;
      output[i*3*cols+j] = (unsigned short)value;
    }
  }

  for(int j=0; j<cols; j++){
    for(int i=0; i<rows; i++){
      for(int k=0; k<3; k++){
	pnm_set_component(imd, i, j, k, output[i*3*cols + (j*3)+k]);
      }
    }
  }

  pnm_save(imd, PnmRawPpm, "YIQ.ppm");
  
  free(output);
  free(out);
  free(tmp);
  pnm_free(ims);
  pnm_free(imd);
}
  



void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 1
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1]);
  return EXIT_SUCCESS;
}
