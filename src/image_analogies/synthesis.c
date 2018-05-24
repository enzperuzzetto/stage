#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bcl.h"
#include "conversion.h"
#include "features.h"


float
dist(int p,  float* mean_tex, float* sd_tex, int q, float* meant, float* sdt)
{

  float  sd, mean;
  mean = (mean_tex[p] - meant[q]) *  (mean_tex[p] - meant[q]);
  sd = (sd_tex[p] - sdt[q]) *  (sd_tex[p] - sdt[q]);
  return sqrtf( mean + sd);
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
  for(int i=0; i<tex_rows; i++){
    for(int j=0; j<tex_cols; j++){
      lum_t[i*cols+j] = lumtex[i*tex_cols+j];
      s[i*cols+j] = i*tex_cols+j;
    }
  }
  
  float d, dist_min;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      q = i*cols+j;

      if(lum_t[q] == 0){

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

	printf( " %f %f %f %f\n",meant[q], sdt[q], meanTex[p], sdTex[p]);
	lum_t[q] = lumtex[p];
	s[q] = p;
      }
    }
  }
  
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
