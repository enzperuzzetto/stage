#include <stdio.h>
#include <stdlib.h>

#include <bcl.h>

float RGB2YUV[3][3] = {
  { 0.299, 0.587, 0.114 },
  { -0.14713, -0.28886, 0.436 },
  { 0.615, -0.51498, -0.10001 }
};

float YUV2RGB[3][3] = {
  { 1.0, 0.0, 1.13983 },
  { 1.0, -0.39465, -0.58060 },
  { 1.0, 2.03211, 0.0 }
};

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

void
convertFloat2UnsignedShort(float* in, unsigned short* out, int cols, int rows)
{
  float val = 0.0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      val = in[i*cols+j];
      if(val > 255)
	val = 255;
      else if(val < 0)
	val = 0.0;
      
      out[i*cols+j] =(unsigned short)val;
    }
  }
}

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

float*
convertRGB2YUV(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, RGB2YUV, in, out);

  return out;
}

float*
convertYUV2RGB(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, YUV2RGB, in, out);

  return out;
}

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

float*
putChannel(int cols, int rows, float* y, float* u, float* v)
{
  float* yuv = malloc(sizeof(float)*cols*3*rows);

  //int k=0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      yuv[i*3*cols+(j*3)] = y[i*cols+j];
      yuv[i*3*cols+(j*3)+1] = u[i*cols+j];
      yuv[i*3*cols+(j*3)+2] = v[i*cols+j];
    }
  }

  return yuv;
}

void
process(char* ims_name)
{
  pnm ims = pnm_load(ims_name);
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);

  pnm imd = pnm_new(cols, rows, PnmRawPpm);
  
  unsigned short* in = pnm_get_image(ims);

  float* tmp = convertUnsignedShort2Float(in, 3 * cols, rows);


  float* yuv = convertRGB2YUV(cols, rows, tmp);

  float* y = channel(cols, rows, yuv, 0);
  float* u = channel(cols, rows, yuv, 1);
  float* v = channel(cols, rows, yuv, 2);

  yuv = putChannel(cols, rows, y, u, v);
  
  float* out = convertYUV2RGB(cols, rows, yuv);
  
  convertFloat2UnsignedShort(out, in, 3*cols, rows);  

  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(imd, i, j, k, in[i*3*cols+(j*3)+k]);
      }
    }
  }
  
  pnm_save(imd, PnmRawPpm, "yuv.ppm");
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
