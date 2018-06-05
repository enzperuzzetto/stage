#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"

#define L 2


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
dist(int tex_cols, int xs, int ys, Pyramid* t, int cols, int rows, int xt, int yt, Pyramid* res, int size, int l)
{  
  int begin_t, end_t, istart, iend, jstart, jend;
  int s = (int)((size-1)/2.0);
  int p = xs * tex_cols + ys, q = xt * cols + yt, p1, q1 ,k = 0;
  double r, g, b, sum = 0.0;
  
  if(xt < s)
    begin_t = -s + (s - xt);
  else
    begin_t = -s;
  if(xt > rows-1 - s)
    end_t = s - (xt - (rows-1 - s));
  else
    end_t = s;

  istart = begin_t;
  iend = end_t;

  if(yt < s)
    begin_t = -s + (s - yt);
  else
    begin_t = -s;
  if(yt > cols-1 - s)
    end_t = s -(yt - (cols-1 - s));
  else
    end_t = s;

  jstart = begin_t;
  jend = end_t;
  
  for(int i=istart; i<iend+1; i++){
    for(int j=jstart; j<jend+1; j++){

      p1 = p + (i*tex_cols+j);
      q1 = q + (i*cols+j);
           
      if(p1 == p || q1 == q)
	break;

      r = (t[l].r[p1] - res[l].r[q1]) * (t[l].r[p1] - res[l].r[q1]);
      g = (t[l].g[p1] - res[l].g[q1]) * (t[l].g[p1] - res[l].g[q1]);
      b = (t[l].b[p1] - res[l].b[q1]) * (t[l].b[p1] - res[l].b[q1]);
      sum += r + g + b;
      k++;
    }
    if(p1 == p || q1 == q)
	break;
  }

  if(k==0)
    sum = 99999999999999999;

  if(l>0){ // niveau l-1
    xs /= 2.0;
    ys /= 2.0;
    xt /= 2.0;
    yt /= 2.0;
    s /= 2.0;
    tex_cols /= 2.0;
    cols /= 2.0;
    rows /= 2.0;
    p = xs * tex_cols +ys;
    q = xt * cols + yt;
    k = 0;
    
    if(xt < s)
      begin_t = -s + (s - xt);
    else
      begin_t = -s;
    if(xt > rows-1 - s)
      end_t = s - (xt - (rows-1 - s));
    else
      end_t = s;

    istart = begin_t;
    iend = end_t;

    if(yt < s)
      begin_t = -s + (s - yt);
    else
      begin_t = -s;
    if(yt > cols-1 - s)
      end_t = s -(yt - (cols-1 - s));
    else
      end_t = s;

    jstart = begin_t;
    jend = end_t;
  
    for(int i=istart; i<iend+1; i++){
      for(int j=jstart; j<jend+1; j++){

	p1 = p + (i*tex_cols+j);
	q1 = q + (i*cols+j);

	r = (t[l-1].r[p1] - res[l-1].r[q1]) * (t[l-1].r[p1] - res[l-1].r[q1]);
	g = (t[l-1].g[p1] - res[l-1].g[q1]) * (t[l-1].g[p1] - res[l-1].g[q1]);
	b = (t[l-1].b[p1] - res[l-1].b[q1]) * (t[l-1].b[p1] - res[l-1].b[q1]);
	sum += r + g + b;
	k++;
      }
    }

    if( k==0 )
      sum = 99999999999999999;
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
  int cols1, cols2, rows1, rows2, pixel, pmin;
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
    
    float val;
    for(int i=0; i<result[l].rows; i++){
      for(int j=0; j<result[l].cols; j++){
	pixel = rand()%(texture[l].cols*texture[l].rows);
	val = rand()%256;
	result[l].r[i*result[l].cols+j] = val;
	result[l].g[i*result[l].cols+j] = val;
	result[l].b[i*result[l].cols+j] = val;
	result[l].s[i*result[l].cols+j] = pixel;
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
#if 1 //Wei Levoy
	pmin = 0;
	
	dist_min = dist(cols1, 0, 0, texture, cols2, rows2, i, j, result, windowSize, l);
	for(int x=0; x<rows1; x++){
	  for(int y=0; y<cols1; y++){
	    pixel = x*cols1+y;
	    //printf(" %f ", texture[l].b[pixel]);
	    d =  dist(cols1, x, y, texture, cols2, rows2, i, j, result, windowSize, l);
       
	    //printf(" d%f dmin%f p%d pmin%d\n",d, dist_min, pixel, p);
	    if(d <= dist_min){
	      dist_min = d;
	      pmin = pixel;
	    }

	  }
	}
#endif	
#if 0  // Ashikhmin
	int xs, ys, start =-(int)(windowSize/2.0), end = (int)(windowSize/2.0);
	
	pmin = result[l].s[q];
	xs = pmin/cols1;
	ys = pmin - xs * cols1;
	dist_min = dist(cols1, xs, ys, texture, cols2, rows2, i, j, result, windowSize, l);
      
	for(int x=start; x<end+1; x++){
	  for(int y=start; y<end+1; y++){
		  
	    pixel = q + (x*cols2+y);

	    if(pixel == q)
	      break;
	    if ( pixel < 0 )
	      continue;
	    
	    p = result[l].s[pixel] - x*cols1-y;
	    
	    if( p >= cols1 * rows1 || p < 0)
	      p =  rand()%(cols1*rows1);

	    xs = p/cols1;
	    ys = p - xs * cols1;
	    
	    d = dist(cols1, xs, ys, texture, cols2, rows2, i, j, result, windowSize, l);
	  
	    if(d < dist_min){
	      dist_min = d;
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
