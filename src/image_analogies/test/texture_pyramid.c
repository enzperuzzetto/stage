#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"

#define L 1


typedef struct Pyramid
{
  int cols;
  int rows;
  float* r;
  float* g;
  float* b;
  int* s;
}Pyramid;

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
    pyramid[l].r = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].g = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].b = malloc(sizeof(float) * dim[0] * dim[1]);
    if( s )
      pyramid[l].s = malloc(sizeof(int) * dim[0] * dim[1]);
    else
      pyramid[l].s = NULL;
    free(dim);
  }
  return pyramid;
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

int
pixelL_1(int p, int cols)
{
  int i = p/cols;
  int j = p - i*cols;
  int new_i = i/2.0;
  int new_j = j/2.0;;

  return new_i * cols/2.0 + new_j;
}

float
dist(int tex_cols,int p, float* Rtex, float* Gtex, float* Btex, int cols, int q, float* Rout, float* Gout, float* Bout, int size, int l)
{
  int i1, j1, p1, q1, k=0;
  float sum=0.0,r, g, b;
  for(int i=0; i< size; i++){
    for(int j=0; j<size; j++){
      i1 = i-size/2.0;
      j1 = j-size/2.0;

      p1 = p + (i1*tex_cols+j1);
      q1 = q + (i1*cols+j1);

      if(p1 < 0 || q1 < 0)
	break;
      if(p1 == p || q1 == q)
	break;
      r = (Rtex[p1] - Rout[q1]) * (Rtex[p1] - Rout[q1]);
      g = (Gtex[p1] - Gout[q1]) * (Gtex[p1] - Gout[q1]);
      b = (Btex[p1] - Bout[q1]) * (Btex[p1] - Bout[q1]);
      sum += r + g + b;
      k++;
    }
    if(p1 == p || q1 == q)
      break;
  }
  if(k==0)
    sum = 99999999999999999;
  
  if(l > 0){
    float sum2 =0.0;
    int q2 = 0.0, p2 = 0.0;
    q2 = pixelL_1(q, cols);
    p2 = pixelL_1(p, tex_cols);
    for(int i=0; i< size; i++){
      for(int j=0; j<size; j++){
	i1 = i-size/2.0;
	j1 = j-size/2.0;

	p1 = p2 + (i1*tex_cols+j1);
	q1 = q2 + (i1*cols+j1);

	if(p1 < 0 || q1 < 0)
	  break;
	if(p1 == p || q1 == q)
	  break;
	r = (Rtex[p1] - Rout[q1]) * (Rtex[p1] - Rout[q1]);
	g = (Gtex[p1] - Gout[q1]) * (Gtex[p1] - Gout[q1]);
	b = (Btex[p1] - Bout[q1]) * (Btex[p1] - Bout[q1]);
	sum2 += r + g + b;
	k++;
      }
      if(p1 == p || q1 == q)
	break;
    }
    if(k==0)
      sum2 = 99999999999999999;
    sum += sum2;
  }
  
  return sqrtf(sum);
}

void
process(char* texture_name, int cols, int rows, int windowSize)
{
  pnm ims = pnm_load(texture_name);
  int tex_cols = pnm_get_width(ims);
  int tex_rows = pnm_get_height(ims);

 

  unsigned short* Redtex = pnm_get_channel(ims, NULL, 0);
  unsigned short* Greentex = pnm_get_channel(ims, NULL, 1);
  unsigned short* Bluetex = pnm_get_channel(ims, NULL, 2);

  float* Rtex = malloc(sizeof(float) * tex_cols * tex_rows);
  float* Gtex = malloc(sizeof(float) * tex_cols * tex_rows);
  float* Btex = malloc(sizeof(float) * tex_cols * tex_rows);

  for(int i=0; i<tex_rows; i++){
    for(int j=0; j<tex_cols; j++){
      Rtex[i*tex_cols+j] = (float)Redtex[i*tex_cols+j];
      Gtex[i*tex_cols+j] = (float)Greentex[i*tex_cols+j];
      Btex[i*tex_cols+j] = (float)Bluetex[i*tex_cols+j];
    }
  }

  /*
  float* Rout = malloc(sizeof(float) * cols * rows);
  float* Gout = malloc(sizeof(float) * cols * rows);
  float* Bout = malloc(sizeof(float) * cols * rows);
  
  
  
  int* s = malloc(sizeof(int)* cols * rows);
  */
  Pyramid* texture = init_pyramid(tex_cols, tex_rows, 1);
  Pyramid* result = init_pyramid(cols, rows, 1);
  int cols1, cols2, rows1, rows2, pixel, pmin, i1,j1;
  //compute pyramid Tex
  
  for(int l=0; l<L;l++){
    free(texture[l].r);
    free(texture[l].g);
    free(texture[l].b);
    if(l==L-1){
      texture[l].r = Rtex;
      texture[l].g = Gtex;
      texture[l].b = Btex;
    }else{
      texture[l].r = pyramid_level(tex_cols, tex_rows, Rtex, L-1-l);
      texture[l].g = pyramid_level(tex_cols, tex_rows, Gtex, L-1-l);
      texture[l].b = pyramid_level(tex_cols, tex_rows, Btex, L-1-l);       
    }
    if(l==0){
      //int r,g,b;
      for(int i=0; i<result[l].rows; i++){
	for(int j=0; j<result[l].cols; j++){
	  pixel = rand()%(texture[l].cols*texture[l].rows);
	  //r = rand()%256;
	  //g = rand()%256;
	  //b = rand()%256;
	  result[l].r[i*result[l].cols+j] = Rtex[pixel];
	  result[l].g[i*result[l].cols+j] = Gtex[pixel];
	  result[l].b[i*result[l].cols+j] = Btex[pixel];
	  result[l].s[i*result[l].cols+j] = pixel;
	}
      }
    }else{
      int q2;
      for(int i=0; i<result[l].rows; i++){
	for(int j=0; j<result[l].cols; j++){
	  pixel = i*result[l].cols+j;
	  q2 = pixelL_1(pixel, result[l].cols);
	  result[l].r[pixel] = result[l-1].r[q2];
	  result[l].g[pixel] = result[l-1].g[q2];
	  result[l].b[pixel] = result[l-1].b[q2];
	  result[l].s[pixel] = result[l-1].s[q2];
	}
      }
    }
  
    
    printf("LEVEL %d\n",l);
    cols1 = texture[l].cols;
    rows1 = texture[l].rows;
    cols2 = result[l].cols;
    rows2 = result[l].rows;
    float d, dist_min;
    int q, p;// pixel;
    
    for(int i=0; i<rows2; i++){
      for(int j=0; j<cols2; j++){
	q = i*cols2+j;
#if 0
	pmin = 0;
	dist_min = dist(cols1, pmin, texture[l].r, texture[l].g, texture[l].b, cols2, q, result[l].r, result[l].g, result[l].b, windowSize, l);
	for(int x=0; x<rows1; x++){
	  for(int y=0; y<cols1; y++){
	    pixel = x*cols1+y;
	    //printf(" %f ", texture[l].b[pixel]);
	    d =  dist(cols1, pixel, texture[l].r, texture[l].g, texture[l].b, cols2, q, result[l].r, result[l].g, result[l].b, windowSize, l);

	    //printf(" d%f dmin%f p%d pmin%d\n",d, dist_min, pixel, p);
	    if(d < dist_min){
	      dist_min = d;
	      pmin = pixel;
	    }

	  }
	}
#endif	
#if 1
	pmin = result[l].s[q];
	dist_min = dist(cols1, pmin, texture[l].r, texture[l].g, texture[l].b, cols2, q, result[l].r, result[l].g, result[l].b, windowSize, l);
      
	for(int x=0; x<windowSize; x++){
	  for(int y=0; y<windowSize; y++){
	    i1 = x - windowSize/2.0;
	    j1 = y - windowSize/2.0;
	  
	    pixel = q + (i1*cols2+j1);

	    if(pixel == q  || pixel < 0 )
	      break;
	      
	    p = result[l].s[pixel] - i1*cols1-j1;
	    if( p >= tex_cols * tex_rows || p < 0)
	      p =  rand()%(cols1*rows1);
	      
	    d = dist(cols1, p, texture[l].r, texture[l].g, texture[l].b, cols2, q, result[l].r, result[l].g, result[l].b, windowSize, l);
	  
	    if(d < dist_min){
	      dist_min = d;
	      //r = pixel;
	      pmin = p;
	    }

	  }
	  if(pixel == q)
	    break;
	}
	#endif
	result[l].r[q] = texture[l].r[pmin];
	result[l].g[q] = texture[l].g[pmin];
	result[l].b[q] = texture[l].b[pmin];
	result[l].s[q] = pmin;
      }
    }
    pnm imt = pnm_new(cols2, rows2, PnmRawPpm);
  
    float valr, valb, valg;
    for(int i=0; i<rows2; i++){
      for(int j=0; j<cols2; j++){
	valr = result[l].r[i*cols2+j];
	valg = result[l].g[i*cols2+j];
	valb = result[l].b[i*cols1+j];
	//printf(" %f %f %f \n", valr, valb, valg);
	if(valr > 255)
	  valr = 255;
	else if ( valr <0)
	  valr = 0;

	if(valg > 255)
	  valg = 255;
	else if ( valg <0)
	  valg = 0;

	if(valb > 255)
	  valb = 255;
	else if ( valb <0)
	  valb = 0;
      
	pnm_set_component(imt, i, j, 0, valr);
	pnm_set_component(imt, i, j, 1, valg);
	pnm_set_component(imt, i, j, 2, valb);
      }
    }
    //printf("\n");
    char name[10];
    sprintf(name, "level%d.ppm",l);
    pnm_save(imt, PnmRawPpm, name);
  }
  /*
    #if 0 // Wei and Levoy texture synthesis
  srand(time(NULL));
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      Rout[i*cols+j] = Rtex[pixel];
      Gout[i*cols+j] = Gtex[pixel];
      Bout[i*cols+j] = Btex[pixel];
      s[i*cols+j] = pixel;
    }
  }
  for(int k=0; k<2;k++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;
	
	p = 0;
	dist_min = dist(tex_cols, p, Rtex, Gtex, Btex,cols, q, Rout, Gout, Bout, windowSize);
	for(int x=0; x<tex_rows; x++){
	  for(int y=0; y<tex_cols; y++){
	    pixel = x*tex_cols+y;

	    d =  dist(tex_cols, pixel, Rtex, Gtex, Btex,cols, q, Rout, Gout, Bout, windowSize);

	    //printf(" d%f dmin%f p%d pmin%d\n",d, dist_min, pixel, p);
	    if(d < dist_min){
	      dist_min = d;
	      p = pixel;
	    }

	  }
	}

	Rout[q] = Rtex[p];
	Gout[q] = Gtex[p];
	Bout[q] = Btex[p];
	s[q] = p;
      }
    }
  }
#endif

#if 0 // Ashikhmin texture synthesis
   
  srand(time(NULL));
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      Rout[i*cols+j] = Rtex[pixel];
      Gout[i*cols+j] = Gtex[pixel];
      Bout[i*cols+j] = Btex[pixel];
      s[i*cols+j] = pixel;
    }
  }
  int i1, j1, pmin;// r, rprim, 
  for(int k=0; k<2; k++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;


	pmin = s[q];
	dist_min = dist(tex_cols, pmin, Rtex, Gtex, Btex,cols, q, Rout, Gout, Bout, windowSize);
      
	for(int x=0; x<windowSize; x++){
	  for(int y=0; y<windowSize; y++){
	    i1 = x - windowSize/2.0;
	    j1 = y - windowSize/2.0;
	  
	    pixel = q + (i1*cols+j1);

	    if((pixel == q && k==0) || pixel < 0 )
	      break;
	      
	    p = s[pixel] - i1*tex_cols-j1;
	    if( p >= tex_cols * tex_rows || p < 0)
	      p =  rand()%(tex_cols*tex_rows);
	      
	    d =  dist(tex_cols, p, Rtex, Gtex, Btex,cols, q, Rout, Gout, Bout, windowSize);
	  
	    if(d < dist_min){
	      dist_min = d;
	      //r = pixel;
	      pmin = p;
	    }

	  }
	  if(pixel == q && k==0)
	    break;
	}

	Rout[q] = Rtex[pmin];
	Gout[q] = Gtex[pmin];
	Bout[q] = Btex[pmin];
	s[q] = pmin;
      }
    }
  }

  
#endif
  
  float valr, valb, valg;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      valr = Rout[i*cols+j];
      valg = Gout[i*cols+j];
      valb = Bout[i*cols+j];
      if(valr > 255)
	valr = 255;
      else if ( valr <0)
	valr = 0;

      if(valg > 255)
	valg = 255;
      else if ( valg <0)
	valg = 0;

      if(valb > 255)
	valb = 255;
      else if ( valb <0)
	valb = 0;
      
      pnm_set_component(imt, i, j, 0, valr);
      pnm_set_component(imt, i, j, 1, valg);
      pnm_set_component(imt, i, j, 2, valb);
    }
  }
  

  pnm_save(imt, PnmRawPpm, "syn.ppm");
  */
  
}



void
usage (char *s){
  fprintf(stderr, "Usage: %s <texture> <width> <height> <sizewindows> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 4
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
    
  process(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  
  return EXIT_SUCCESS;
}
