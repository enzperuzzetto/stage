#include <math.h>

#include "features.h"


float
mean(int cols, int rows, float* in, int index, int size, int half)
{
  float mean = 0.0;
  int k = 0, i_p1 = 0, j_p1 = 0, pixel = 0;
  
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      i_p1 = i - size/2.0;
      j_p1 = j - size/2.0;
      pixel = index + (i_p1*cols+j_p1);

      if(half && (pixel == index))
	break;

      if(pixel >= 0 && pixel < rows*cols){
	mean += in[pixel];
	k++;
      }
    }

    if(half && (pixel == index))
	break;
  }
  if(k==0)
    mean = 0.0;
  else
    mean /= (float)k;
  
  return mean;
}

float
variance(int cols, int rows, float* in, float mean, int index, int size, int half)
{
  float var = 0.0;
  int k=0, i_p1 = 0, j_p1 = 0, pixel = 0;

  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      i_p1 = i - size/2.0;
      j_p1 = j - size/2.0;
      pixel = index + (i_p1*cols+j_p1);
      
      if(half && pixel == index)
	break;

      if(pixel >= 0 && pixel < rows*cols){
	var += (in[pixel] - mean) * (in[pixel] - mean);
	k++;
      }
    }
    if(half && (pixel == index))
	break;
  }

  if(k == 0)
    var = 0.0;
  else
    var /= (float)k;
  
  return var;
}

float
standard_deviation(int cols, int rows, float* in, float mean, int index, int size, int half)
{
  float var=0.0, st_d=0.0;
  var = variance(cols, rows, in, mean, index, size, half);
  st_d = sqrt(var);
  return st_d;
}

void
mean_img(int cols, int rows, float* in, float* means, int size, int half)
{
  int pixel = 0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = i*cols+j;
      means[pixel] = mean(cols, rows, in, pixel, size, half);
    }
  }
}


void
sd_img(int cols, int rows, float* in, float* mean, float* sd, int size, int half)
{
  int pixel = 0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = i*cols+j;
      sd[pixel] = standard_deviation(cols, rows, in, mean[pixel], pixel, size, half);
    }
  }
}


int
max(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}

float
mean_all_img(int cols, int rows, float* in)
{
  float m = 0.0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      m += in[i*cols+j];
    }
  }

  return m/(float)(rows*cols);
}

float
sd_all_img(int cols, int rows, float* in, float mean)
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
