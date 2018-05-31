#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"


float
dist(int tex_cols,int p, float* Rtex, float* Gtex, float* Btex, int cols, int q, float* Rout, float* Gout, float* Bout, int size)
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
    return 99999999999999999;
  return sqrtf(sum);
}

void
process(char* texture_name, int cols, int rows, int windowSize)
{
  pnm ims = pnm_load(texture_name);
  int tex_cols = pnm_get_width(ims);
  int tex_rows = pnm_get_height(ims);

  pnm imt = pnm_new(cols, rows, PnmRawPpm);

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

  float* Rout = malloc(sizeof(float) * cols * rows);
  float* Gout = malloc(sizeof(float) * cols * rows);
  float* Bout = malloc(sizeof(float) * cols * rows);
  
  
  
  int* s = malloc(sizeof(int)* cols * rows);
  
  float d, dist_min;
  int q, p, pixel;
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

#if 1 // Ashikhmin texture synthesis
   
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
