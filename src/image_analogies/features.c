#include <features.h>


void
multiple(int cols, int rows, float M[3][3], float* B, float* output)
{
   for(int j=0; j<3 * cols; j += 3){
    for(int i=0; i<rows; i++){
      output[i*3*cols + j  ] = M[0][0] * B[i*3*cols + j] + M[0][1] * B[i*3*cols + j+1] + M[0][2] * B[i*3*cols + j+2];
      output[i*3*cols + j+1] = M[1][0] * B[i*3*cols + j] + M[1][1] * B[i*3*cols + j+1] + M[1][2] * B[i*3*cols + j+2];
      output[i*3*cols + j+2] = M[2][0] * B[i*3*cols + j] + M[2][1] * B[i*3*cols + j+1] + M[2][2] * B[i*3*cols + j+2];
    }
  }
}

void
convertRGB2YIQ(int cols, int rows, float* input, float* output)
{
   multiple(cols, rows, RGB2YIQ, input, output);
}

void
convertYIQ2RGB(int cols, int rows, float* input, float* output)
{
   multiple(cols, rows, YIQ2RGB, input, output);
}
