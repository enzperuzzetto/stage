#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"
#include "conversion.h"
#include "features.h"


float
dist(int p, float* mean_tex, float* sd_tex, int q, float* meant, float* sdt)
{

  float  sd, mean;
  mean = (mean_tex[p] - meant[q]) *  (mean_tex[p] - meant[q]);
  sd = (sd_tex[p] - sdt[q]) *  (sd_tex[p] - sdt[q]);
  return sqrtf( mean + sd);
}

float
dist1(int tex_cols,int p, float* lumtex, float* mean_tex, float* sd_tex, int cols, int q, float* lumt, float* meant, float* sdt, int size)
{
  int i1, j1, p1, q1;
  float lum = 0.0, sd, mean;
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

      lum += (lumtex[p1] -lumt[q1]) * (lumtex[p1] - lumt[q1]);
    }
    if(p1 == p || q1 == q)
      break;
  }

  mean = (mean_tex[p] - meant[q]) *  (mean_tex[p] - meant[q]);
  sd = (sd_tex[p] - sdt[q]) *  (sd_tex[p] - sdt[q]);

  return sqrt(lum + mean + sd);
}
      

void
process(char* texture_name, int cols, int rows, int windowSize)
{
  pnm tex = pnm_load(texture_name);
  int tex_cols = pnm_get_width(tex);
  int tex_rows = pnm_get_height(tex);

  pnm imt = pnm_new(cols, rows, PnmRawPpm);

  unsigned short* RGBtex = pnm_get_image(tex);
  float* tmp = convertUnsignedShort2Float(RGBtex, 3 * tex_cols, tex_rows);

  float* YIQtex = convertRGB2YIQ(tex_cols, tex_rows, tmp);

  float* lumtex = channel(tex_cols, tex_rows, YIQtex, 0);

  float* meanTex = malloc(sizeof(float) * tex_cols * tex_rows);
  float* sdTex = malloc(sizeof(float) * tex_cols * tex_rows);

  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * cols * rows);
  float* lum_t = malloc(sizeof(float) * cols * rows);
  float* meant = malloc(sizeof(float) * cols * rows);
  float* sdt = malloc(sizeof(float) * cols * rows);
  int* s = malloc(sizeof(int) * cols * rows);
  
  mean_img(tex_cols, tex_rows, lumtex, meanTex, windowSize, 1);
  sd_img(tex_cols, tex_rows, lumtex, meanTex, sdTex, windowSize, 1);

  int q, p,pixel;
  /*for(int i=0; i<tex_rows; i++){
    for(int j=0; j<tex_cols; j++){
      lum_t[i*cols+j] = lumtex[i*tex_cols+j];
      s[i*cols+j] = i*tex_cols+j;
    }
    }*/
  
  float d, dist_min;
#if 0 // Wei and Levoy texture synthesis
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      q = i*cols+j;

      //if(lum_t[q] == 0){

      meant[q] = mean(cols, rows, lum_t, q, windowSize, 1);
      sdt[q] =  standard_deviation(cols, rows, lum_t, meant[q], q, windowSize, 1);
	
      p = 0;
      dist_min = dist(p, meanTex, sdTex, q,  meant, sdt);
      for(int x=0; x<tex_rows; x++){
	for(int y=0; y<tex_cols; y++){
	  pixel = x*tex_cols+y;

	  d = dist(pixel, meanTex, sdTex, q,  meant, sdt);

	  if(d < dist_min){
	    dist_min = d;
	    p = pixel;
	  }

	}
      }

      //printf( " %f %f %f %f\n",meant[q], sdt[q], meanTex[p], sdTex[p]);
      lum_t[q] = lumtex[p];
      s[q] = p;
      // }
    }
  }
#endif

#if 1 // Ashikhmin texture synthesis
   
  srand(time(NULL));
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      s[i*cols+j] = rand()%(tex_cols*tex_rows);
    }
  }
  
  int i1, j1, r=-1, rprime;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      q = i*cols+j;

      meant[q] = mean(cols, rows, lum_t, q, windowSize, 1);
      sdt[q] =  standard_deviation(cols, rows, lum_t, meant[q], q, windowSize, 1);

      
      r = q + ((-windowSize/2.0)*cols-(windowSize/2.0));
      if( r < 0 ){
	for(int x=0; x<windowSize; x++){
	  for(int y=0; y<windowSize; y++){
	    i1 = x - windowSize/2.0;
	    j1 = y - windowSize/2.0;
	  
	    r = q + (i1*cols+j1);

	    if( r >= 0 || r == q)
	      break;
	  }
	  if( r >= 0 || r == q)
	    break;
	}

	if(r == q)
	  r = -1;
      }
     
      if(r >= 0){
	dist_min = dist1(tex_cols, r, lumtex, meanTex, sdTex, cols, q, lum_t, meant, sdt, windowSize);
      
	for(int x=0; x<windowSize; x++){
	  for(int y=0; y<windowSize; y++){
	    i1 = x - windowSize/2.0;
	    j1 = y - windowSize/2.0;
	  
	    pixel = q + (i1*cols+j1);

	    if(pixel == q || pixel < 0 )
	      break;

	    if( cols > tex_cols)
	      p = s[pixel] + (q-pixel) + (cols - tex_cols);
	    else if( cols < tex_cols)
	      p = s[pixel] + (q-pixel) + (tex_cols - cols);
	    else
	      p = s[pixel] + (q-pixel);

	    if( p >= tex_cols * tex_rows)
	      p = s[pixel];

	    //d = dist(p, meanTex, sdTex, q,  meant, sdt);
	    d = dist1(tex_cols, p, lumtex, meanTex, sdTex, cols, q, lum_t, meant, sdt, windowSize);
	  
	    if(d < dist_min){
	      dist_min = d;
	      r = pixel;
	      rprime = -i1*tex_cols-j1;
	    }

	  }
	  if(pixel == q)
	    break;
	}

	p = s[r] + rprime;

	if( p > tex_cols*tex_rows-1)
	  p = s[r];
      
	lum_t[q] = lumtex[p];
	s[q] = p;
      }

      else
	lum_t[q] = lumtex[s[q]];
       
    }
  }
  
#endif
  float* isf = channel(tex_cols, tex_rows, YIQtex, 1);
  float* qsf = channel(tex_cols, tex_rows, YIQtex, 2);

  float* it = malloc(sizeof(float) * cols * rows);
  float* qt = malloc(sizeof(float) * cols * rows);

  pixel = 0;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = s[i*cols+j];
      it[i*cols+j] = isf[pixel];
      qt[i*cols+j] = qsf[pixel];
    }
  }

  float* yiq = putChannel(cols, rows, lum_t, it, qt);

  float* rgb = convertYIQ2RGB(cols, rows, yiq);

  convertFloat2UnsignedShort(rgb, data_out, 3*cols, rows);

  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(imt, i, j, k, data_out[i*3*cols+(j*3)+k]);
      }
    }
  }

  pnm_save(imt, PnmRawPpm, "syn.ppm");
  
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
