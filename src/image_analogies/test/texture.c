#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "bcl.h"

float
dist(int tex_cols, int xs, int ys, float* Rtex, float* Gtex, float* Btex, int cols, int rows, int xt, int yt, float* Rout, float* Gout, float* Bout, int size, int nb_pass)
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
           
      if((p1 == p || q1 == q) && nb_pass == 0)
	break;
      if(p1 <0)
	continue;
     
      r = (Rtex[p1] - Rout[q1]) * (Rtex[p1] - Rout[q1]);
      g = (Gtex[p1] - Gout[q1]) * (Gtex[p1] - Gout[q1]);
      b = (Btex[p1] - Bout[q1]) * (Btex[p1] - Bout[q1]);
      sum += r + g + b;
      k++;
    }
    if((p1 == p || q1 == q) && nb_pass == 0)
      break;
  }
  
  if(k==0)
    sum = 99999999999999999;
  
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
  
  int pixel, q, p;
  float dist_min, d;
#if 0 // Wei and Levoy texture synthesis

  //bruit blanc
  srand(time(NULL));
  int val;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      val = rand()%256;
      Rout[i*cols+j] = val;
      Gout[i*cols+j] = val;
      Bout[i*cols+j] = val;
      s[i*cols+j] = pixel;
    }
  }

  for(int k=0; k<1;k++){ // nb passes
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;

       
	p = 0;
	
	dist_min = dist(tex_cols,  0, 0, Rtex, Gtex, Btex, cols, rows, i, j, Rout, Gout, Bout, windowSize,k);
	for(int x=0; x<tex_rows; x++){
	  for(int y=0; y<tex_cols; y++){
	    pixel = x*tex_cols+y;

	    d = dist(tex_cols, x, y, Rtex, Gtex, Btex, cols, rows, i, j, Rout, Gout, Bout, windowSize,k);
	
	    if(d <= dist_min){
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
  int val;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = rand()%(tex_cols*tex_rows);
      val = rand()%256;
      Rout[i*cols+j] = val;
      Gout[i*cols+j] = val;
      Bout[i*cols+j] = val;
      s[i*cols+j] = pixel;
    }
  }
  
  int pmin, xs, ys, start =-(int)(windowSize/2.0), end = (int)(windowSize/2.0);
  for(int k=0; k<1; k++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
	q = i*cols+j;

	pmin = s[q];
	xs = pmin/tex_cols;
	ys = pmin - xs * tex_cols;
 
	dist_min = dist(tex_cols, xs, ys, Rtex, Gtex, Btex, cols, rows, i, j, Rout, Gout, Bout, windowSize,k);
      
	for(int x=start; x<end+1; x++){
	  for(int y=start; y<end+1; y++){
	    
	    pixel = q + (x*cols+y);

	    if(pixel == q && k==0)
	      break;
	    if( pixel < 0 )
	      continue;
	    
	    p = s[pixel] - x*tex_cols-y;
	    
	    if( p >= tex_cols * tex_rows || p < 0)
	      p =  rand()%(tex_cols*tex_rows);

	    xs = p/tex_cols;
	    ys = p - xs * tex_cols;
	    
	    d = dist(tex_cols,  xs, ys, Rtex, Gtex, Btex, cols, rows, i, j, Rout, Gout, Bout, windowSize,k);
	    
	    if(d < dist_min){
	      dist_min = d;
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
  

  free(s);
  free(Bout);
  free(Gout);
  free(Rout);
  free(Rtex);
  free(Gtex);
  free(Btex);
  free(Redtex);
  free(Greentex);
  free(Bluetex);
  pnm_free(imt);
  pnm_free(ims);
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
