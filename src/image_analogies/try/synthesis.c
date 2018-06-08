#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"
#include "conversion.h"
#include "features.h"
#include "pyramid.h"


float
dist(int tex_cols,int xs, int ys, float* lumtex, float* mean_tex, float* sd_tex, int cols, int rows, int xt, int yt, float* lumt, float* meant, float* sdt, int size, int nb_pass)
{
  int begin_t, end_t, istart, iend, jstart, jend;
  int s = (int)((size-1)/2.0);
  int p = xs * tex_cols + ys, q = xt * cols + yt, p1, q1 ,k = 0;
  double lum = 0.0, tmp;//mean = 0.0, sd = 0.0, tmp;
  (void)mean_tex;
  (void)meant;
  (void)sd_tex;
  (void)sdt;
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
           
      if((p1 == p || q1 == q) && nb_pass == 0)
	break;
     
      tmp = lumtex[p1] - lumt[q1];
      lum += tmp * tmp;
      k++;
    }
    if((p1 == p || q1 == q) && nb_pass == 0)
      break;
  }
  
  if(k==0)
    lum = 99999999999999999;

  /*tmp = mean_tex[p] - meant[q];
  mean = tmp * tmp;

  tmp = sd_tex[p] - sdt[q];
  sd = tmp * tmp;*/
  
  return sqrtf(lum); //+ mean + sd);
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
  float d, dist_min;
  
#if 0 // Wei and Levoy texture synthesis
  //bruit blanc
  srand(time(NULL));
  int val;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      val = rand()%256;
      lum_t[i*cols+j] = val;
      s[i*cols+j] = pixel;
    }
  }
 
  for(int k=0; k<1; k++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;

	meant[q] = mean(cols, rows, lum_t, q, windowSize, 1);
	sdt[q] =  standard_deviation(cols, rows, lum_t, meant[q], q, windowSize, 1);
	
	p = 0;
	
	dist_min = dist(tex_cols, 0, 0, lumtex, meanTex, sdTex, cols, rows, i, j, lum_t, meant, sdt, windowSize,k);
	for(int x=0; x<tex_rows; x++){
	  for(int y=0; y<tex_cols; y++){
	    pixel = x*tex_cols+y;

	    d = dist(tex_cols, x, y, lumtex, meanTex, sdTex, cols, rows, i, j, lum_t, meant, sdt, windowSize,k);
	    if(d < dist_min){
	      dist_min = d;
	      p = pixel;
	    }

	  }
	}

	lum_t[q] = lumtex[p];
	s[q] = p;
      }
    }
  }
#endif

#if 1 // Ashikhmin texture synthesis
   
  srand(time(NULL));
  int val;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      val = rand()%256;
      lum_t[i*cols+j] = val;
      s[i*cols+j] = pixel;
    }
  }
  
  int pmin, xs, ys, start =-(int)(windowSize/2.0), end = (int)(windowSize/2.0);
  for(int k=0; k<1; k++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;

	meant[q] = mean(cols, rows, lum_t, q, windowSize, 1);
	sdt[q] =  standard_deviation(cols, rows, lum_t, meant[q], q, windowSize, 1);

	pmin = s[q];
	xs = pmin/tex_cols;
	ys = pmin - xs * tex_cols;
	
	dist_min = dist(tex_cols, xs, ys, lumtex, meanTex, sdTex, cols, rows, i, j, lum_t, meant, sdt, windowSize,k);
      
	for(int x=start; x<end+1; x++){
	  for(int y=start; y<end+1; y++){
	    
	    pixel = q + (x*cols+y);

	    if(pixel == q && k==0) 
	      break;
	    if(pixel < 0)
	      continue;

	    p = s[pixel] - x*tex_cols-y;
	    
	    if( p >= tex_cols * tex_rows || p < 0)
	      p =  rand()%(tex_cols*tex_rows);
		
	    xs = p/tex_cols;
	    ys = p - xs * tex_cols;
	    
	    d = dist(tex_cols, xs, ys, lumtex, meanTex, sdTex, cols, rows, i, j, lum_t, meant, sdt, windowSize,k);
	  
	    if(d < dist_min){
	      dist_min = d;
	      pmin = p;
	    }

	  }
	  if(pixel == q && k==0)
	    break;
	}
	  
	lum_t[q] = lumtex[pmin];
	s[q] = pmin;
      }	
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
