#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <bcl.h>

#define NEIGHBOOR_SIZE_COARSER 3
#define NEIGHBOOR_SIZE_FINER 5
#define K 3
#define L 3
#define NB_MAX_ITER 7
#define ALPHA 0.5



//************* Structure pyramid ******************************//

typedef struct Pyramid
{
  int cols;
  int rows;
  float* data;
  float* mean;
  float* sd;
  int* s;
}Pyramid;

//************** Convertion function ***************************//

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

float GAUSS[5] = { 0.25-0.40/2.0, 0.25, 0.40, 0.25, 0.25-0.40/2.0 };

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
  for(int j=0; j<cols; j++){
    for(int i=0; i<rows; i++){
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
convertRGB2YIQ(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, RGB2YIQ, in, out);

  return out;
}

float*
convertYIQ2RGB(int cols, int rows, float* in)
{
  float* out = malloc(sizeof(float)* rows * 3 * cols);

  multiple(cols, rows, YIQ2RGB, in, out);

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
putChannel(int cols, int rows, float* y, float* I, float* q)
{
  float* yiq = malloc(sizeof(float)*cols*3*rows);

  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      yiq[i*3*cols+(j*3)] = y[i*cols+j];
      yiq[i*3*cols+(j*3)+1] = I[i*cols+j];
      yiq[i*3*cols+(j*3)+2] = q[i*cols+j];
    }
  }

  return yiq;
}
//************************************************************************************//

//*************** Pyramid function  **************************************************//

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
  
  
//*************** Features function **************************************************//

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

      if(half && pixel == index)
	break;

      if(pixel >= 0 && pixel < rows*cols){
	mean += in[pixel];
	k++;
      }
    }
  }

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
  }
   
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

int
max(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}

//****************************************************************************************************//

//********************* Match function ***************************************************************//

int*
initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
{
  int* nnf = malloc(sizeof(int) * cols_t * rows_t);
  srand(time(NULL));
  for(int i=0; i<cols_t*rows_t; i++){
    nnf[i] =  rand()%(cols_s*rows_s);
  }

  return nnf;
}

int
propagation(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
  int pixel_min = q;
  float dist_min = dist(nnf[pixel_min], source, source_filter, q, target, target_filter, l);

  int left = q - 1;

  if(left < 0)
    left = 0;

  float d = dist(nnf[left], source, source_filter, q, target, target_filter, l);

  if(d < dist_min){
    dist_min = d;
    pixel_min = left;
  }

  int top = q - target[l].cols;

  if(top < 0)
    top = 0;

  d = dist(nnf[top], source, source_filter, q, target, target_filter, l);

  if(d < dist_min){
    dist_min = d;
    pixel_min = top;
  }

  return pixel_min;
}

int
search(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
  int w = max(source[l].cols, source[l].rows), pixel_x = 0, pixel_y = 0, pixel = 0, pixel_min = nnf[q], i = 0;
  float dist_min = dist(pixel_min, source, source_filter, q, target, target_filter, l), d = 0.0, alpha = powf(ALPHA, i), x = 0.0, y = 0.0;

  while( w*alpha > 1){
    x = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;
    y = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;

    pixel_x = w*alpha*x;
    pixel_y = w*alpha*y;
    pixel = pixel_min + (pixel_x*source[l].cols + pixel_y);

    if(pixel < 0)
      pixel = 0;
    else if(pixel > source[l].cols*source[l].rows-1)
      pixel = source[l].cols*source[l].rows-1;

    d = dist(pixel, source, source_filter, q, target, target_filter, l);

    if( d < dist_min){
      dist_min = d;
      pixel_min = pixel;
    }

    i++;
    alpha = powf(ALPHA, i);
  }

  return pixel_min;
}

int
iteration(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
  int pixelTarget = propagation(source, source_filter, target, target_filter, l, q, nnf);
  int pixelSource = search(source, source_filter, target, target_filter, l, pixelTarget, nnf);

  return pixelSource;
}

int
patchMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int pixelMatch = 0, tmp = 0;
  for(int k=0; k<NB_MAX_ITER; k++){
    int* nnf = initialisation(source[l].cols, source[l].rows, target[l].cols, target[l].rows);

    tmp = iteration(source, source_filter, target, target_filter, l, q, nnf);

    if(k > 0){
      if( dist(pixelMatch, source, source_filter, q, target, target_filter, l) > dist(tmp, source, source_filter, q, target, target_filter, l) ){
	pixelMatch = tmp;
      }
    }
    else{
      pixelMatch = tmp;
    }

    free(nnf);
  }

  return pixelMatch;
    
}


int
bruteForceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int index_min = 0;
  int pixel = 0;
  float dist_min = dist(0, source, source_filter, q, target, target_filter, l);
  float d = 0.0;
  for(int i=0; i<source[l].rows; i++){
    for(int j=0; j<source[l].cols; j++){
      pixel = i*source[l].cols+j;
      d = dist(pixel, source, source_filter, q, target, target_filter, l);
      if(d < dist_min){
	dist_min = d;
	index_min = pixel;
      }
    }
  }
  
  return index_min;
}

int
bestCoherenceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int pixel = 0, p = 0;
  int r = 0, i_p1 = 0, j_p1 = 0;
  float dist_min = dist(0, source, source_filter, q, target, target_filter, l);
  float d = 0.0;
  for(int i=0; i<NEIGHBOOR_SIZE_FINER; i++){
    for(int j=0; j<NEIGHBOOR_SIZE_FINER; j++){
      i_p1 = i - NEIGHBOOR_SIZE_FINER/2.0;
      j_p1 = j - NEIGHBOOR_SIZE_FINER/2.0;
      pixel = q + (i_p1*target_filter[l].cols+j_p1);

      if(pixel == q)
	break;
      if( pixel < 0 )
	pixel = 0;
      else if( pixel > target_filter[l].cols * target_filter[l].rows -1 )
	pixel = source[l].cols * source[l].rows -1;

      
      p = target_filter[l].s[pixel] + (q-pixel);
      
      if(p < 0)
	p = 0;
      else if(p > source[l].cols * source[l].rows -1)
	p = source[l].cols * source[l].rows -1;
      
      d = dist(p, source, source_filter, q, target, target_filter, l);
   
      if(d < dist_min){
	r = pixel;
	dist_min = d;
      }
    }
  }

  p = target_filter[l].s[r] + (q - r);
  
  if(p < 0)
    p = 0;
  else if(p > source[l].cols * source[l].rows -1)
    p = source[l].cols * source[l].rows -1;
  
  return p;
}

int
bestMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  //int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
  int p_patch = patchMatch(source, source_filter, target, target_filter, l, q);
  int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);
  

  //float dist_brute = dist(p_brute, source, source_filter, q, target, target_filter, l);
  float dist_patch = dist(p_patch, source, source_filter, q, target, target_filter, l);
  dist_patch *= dist_patch;
  float dist_coh = dist(p_coh, source, source_filter, q, target, target_filter, l);
  dist_coh *= dist_coh;
  
  if(dist_coh < dist_patch * ( 1 + powf(2, l-L)*K) )
    return p_coh;
  else
  return p_patch;
  //return p_brute;
}

//**********************************************************************************************//

//**************** Main ***********************************************************************//

void
process(char* ims_name, char* ims_filter_name, char* imt_name)
{
  pnm imsRGB = pnm_load(ims_name);
  pnm ims_filterRGB = pnm_load(ims_filter_name);
  int cols_s = pnm_get_width(imsRGB);
  int rows_s = pnm_get_height(imsRGB);

  pnm imtRGB = pnm_load(imt_name);
  int cols_t = pnm_get_width(imtRGB);
  int rows_t = pnm_get_height(imtRGB);

  pnm imt_filter = pnm_new(cols_t, rows_t, PnmRawPpm);

  //*********** covert to YIQ ***************//
  
  printf("Convertion en YIQ\n");
  unsigned short* sRGB = pnm_get_image(imsRGB);
  unsigned short* s_fRGB = pnm_get_image(ims_filterRGB);
  unsigned short* tRGB = pnm_get_image(imtRGB);

  float* tmp_s = convertUnsignedShort2Float(sRGB, 3 * cols_s, rows_s);
  float* tmp_sf = convertUnsignedShort2Float(s_fRGB, 3 * cols_s, rows_s);
  float* tmp_t = convertUnsignedShort2Float(tRGB, 3 * cols_t, rows_t);
  
  float* s_yiq = convertRGB2YIQ(cols_s, rows_s, tmp_s);
  float* sf_yiq = convertRGB2YIQ(cols_s, rows_s, tmp_sf);
  float* t_yiq = convertRGB2YIQ(cols_t, rows_t, tmp_t);

  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * cols_t * rows_t);
 
  float* data_source = channel(cols_s, rows_s, s_yiq, 0);
  float* data_source_filter = channel(cols_s, rows_s, sf_yiq, 0);
  float* data_target = channel(cols_t, rows_t, t_yiq, 0);
  float* data_target_filter = malloc(sizeof(float) * cols_t * rows_t);

  printf("fait\n");

  //*****************************************//

  
  //*********** Compute L level of pyramid *****************************//
  printf("Démarrage algo \n");
    
  Pyramid* source = init_pyramid(cols_s, rows_s, 0);
  Pyramid* source_filter = init_pyramid(cols_s, rows_s, 0);
  Pyramid* target = init_pyramid(cols_t, rows_t, 0);
  Pyramid* target_filter = init_pyramid(cols_t, rows_t, 1);
    
  for(int l=0; l<L; l++){
    printf("Début level %d\n",l);
     
    if( l == L-1) {
      free(source[l].data);
      free(source_filter[l].data);
      free(target[l].data);
      
      source[l].data = data_source;
      source_filter[l].data = data_source_filter;
      target[l].data = data_target;

    }else{
      free(source[l].data);
      free(source_filter[l].data);
      free(target[l].data);
      free(target_filter[l].data);
      
      source[l].data = pyramid_level(cols_s, rows_s, data_source, L-1-l);
      source_filter[l].data = pyramid_level(cols_s, rows_s, data_source_filter, L-1-l);
           
      target[l].data = pyramid_level(cols_t, rows_t, data_target, L-1-l);
      target_filter[l].data = pyramid_level(cols_t, rows_t, data_target_filter, L-1-l);
    }
    //************* Compute Features Means & standard deviation **************************
      
    mean_img(source[l].cols, source[l].rows, source[l].data, source[l].mean, NEIGHBOOR_SIZE_FINER, 0);
    mean_img(target[l].cols, target[l].rows, target[l].data, target[l].mean, NEIGHBOOR_SIZE_FINER, 0);      
      
    sd_img(source[l].cols, source[l].rows, source[l].data, source[l].mean, source[l].sd, NEIGHBOOR_SIZE_FINER, 0);
    sd_img(target[l].cols, target[l].rows, target[l].data, target[l].mean, target[l].sd, NEIGHBOOR_SIZE_FINER, 0);
      
    //*********** Compute feature l-1 with neighboor size coarser ***********************// 
    if( l > 0){
      mean_img(source[l-1].cols, source[l-1].rows, source[l-1].data, source[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(source_filter[l-1].cols, source_filter[l-1].rows, source_filter[l-1].data, source_filter[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(target[l-1].cols, target[l-1].rows, target[l-1].data, target[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(target_filter[l-1].cols, target_filter[l-1].rows, target_filter[l-1].data, target_filter[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
	
      sd_img(source[l-1].cols, source[l-1].rows, source[l-1].data, source[l-1].mean, source[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(source_filter[l-1].cols, source_filter[l-1].rows, source_filter[l-1].data, source_filter[l-1].mean, source_filter[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(target[l-1].cols, target[l-1].rows, target[l-1].data, target[l-1].mean, target[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(target_filter[l-1].cols, target_filter[l-1].rows, target_filter[l-1].data, target_filter[l-1].mean, target_filter[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
    }
  
    //*************************************************************************************

    //************* Best Match   **********************************************************

    int p = 0, q = 0;

    for(int i=0; i<target_filter[l].rows; i++){
      for(int j=0; j<target_filter[l].cols; j++){

	//compute mean & sd to q in A' and B'
	q = i*target_filter[l].cols+j;

	source_filter[l].mean[q] = mean(source_filter[l].cols, source_filter[l].rows, source_filter[l].data, q, NEIGHBOOR_SIZE_FINER, 1);
	source_filter[l].sd[q] = standard_deviation(source_filter[l].cols, source_filter[l].rows, source_filter[l].data, source_filter[l].mean[q], q, NEIGHBOOR_SIZE_FINER, 1);
	        
	target_filter[l].mean[q] = mean(target_filter[l].cols, target_filter[l].rows, target_filter[l].data, q, NEIGHBOOR_SIZE_FINER, 1);
	target_filter[l].sd[q] = standard_deviation(target_filter[l].cols, target_filter[l].rows, target_filter[l].data, target_filter[l].mean[q], q, NEIGHBOOR_SIZE_FINER, 1);
	      
	p = bestMatch(source, source_filter, target, target_filter, l, q);
	 
	target_filter[l].data[q] = source_filter[l].data[p];
	target_filter[l].s[q] = p;
      }
    }
      
    printf("Fait \n");
  }

  //*************************************************************************************

  printf("Convertion en RGB\n");
       
  //************* Copy I & Q B to B' ************************************//
    
  float* it = channel(cols_t, rows_t, t_yiq, 1);
  float* qt = channel(cols_t, rows_t, t_yiq, 2);

  float* yiq = putChannel(cols_t, rows_t, target_filter[L-1].data, it, qt);

  float* rgb = convertYIQ2RGB(cols_t, rows_t, yiq);
  
  convertFloat2UnsignedShort(rgb, data_out, 3*cols_t, rows_t);

  //*********************************************************************//

  free_pyramid(source, 0);
  free_pyramid(source_filter, 0);
  free_pyramid(target, 0);
  free_pyramid(target_filter, 1);

  for(int i=0; i<rows_t; i++){
    for(int j=0; j<cols_t; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(imt_filter, i, j, k, data_out[i*3*cols_t+(j*3)+k]);
      }
    }
  }
  

  printf("fait\n");
 
  printf("Sauvegarde\n");
  pnm_save(imt_filter, PnmRawPpm, "filter.ppm");
  printf("Fait\n");


  free(imsRGB);
  free(ims_filterRGB);
  free(imtRGB);
  pnm_free(imt_filter);
  free(sRGB);
  free(s_fRGB);
  free(tRGB);
  free(tmp_s);
  free(tmp_sf);
  free(tmp_t);
  free(s_yiq);
  free(sf_yiq);
  free(t_yiq);
  free(data_out);
  free(data_target_filter);
  free(it);
  free(qt);
  free(yiq);
  free(rgb);
}

void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <ims_filter> <imt> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 3
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1], argv[2], argv[3]);
  return EXIT_SUCCESS;
}

//**********************************************************************************//
