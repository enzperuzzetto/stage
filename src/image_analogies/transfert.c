#include "transfert.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bcl.h"
#include "pyramid.h"
#include "features.h"
#include "conversion.h"


typedef struct Pyramid
{
  int cols;
  int rows;
  float* lum;
  float* mean;
  float* sd;
  int* s;
}Pyramid;

Pyramid*
init_pyramid_transfert(int cols, int rows)
{
  Pyramid* pyramid = malloc(sizeof(Pyramid) * L);

  for(int l=0; l<L; l++){
    int* dim = pyramid_level_dim(cols, rows, L-1-l);
    pyramid[l].cols = dim[0];
    pyramid[l].rows = dim[1];
    pyramid[l].lum = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].mean = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].sd = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].s = malloc(sizeof(int) * dim[0] * dim[1]);

    free(dim);
  }

  return pyramid;
}

void
free_pyramid_transfert(Pyramid* pyramid)
{
  for(int l=0; l<L; l++){
    free(pyramid[l].lum);
    free(pyramid[l].mean);
    free(pyramid[l].sd);
    free(pyramid[l].s);
  }

  free(pyramid);
}

//####################################################################################################################
float
dist_transfert(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, int w)
{
   int s = (int)(size/2.0);
   int Acols, Arows, Bcols, Brows, ip, jp, iq, jq, p1, q1, q, k=0;
   float lum1 = 0.0, lum2 = 0.0, d = 0.0, dprim = 0.0, tmp;

   Acols = A[l].cols;
   Arows = A[l].rows;
   Bcols = B[l].cols;
   Brows = B[l].rows;
   q = xt*Bcols+yt;

   for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){
      ip = xs + i;
      jp = ys + j;
      iq = xt + i;
      jq = yt + j;

      if(ip < 0)
	ip = Arows - ip;
      if(jp < 0)
	jp = Acols - jp;
      if(ip > Arows-1)
	ip -= Arows;
      if(jp > Acols-1)
	jp-= Acols;

      if(iq < 0)
	iq = Brows - iq;
      if(jq < 0)
	jq = Bcols - jq;
      if(iq > Brows-1)
	iq -= Brows;
      if(jq > Bcols-1)
	jq -= Bcols;

      p1 = ip*Acols+jp;
      q1 = iq*Bcols+jq;
      /*
      if((q1 == q) && nb_pass == 0)
	break;
      */
      tmp = A[l].lum[p1] - B[l].lum[q1];
      lum1 += tmp * tmp;
      k++;
    }
   }

   d = lum1 / (float)k;

   k=0;

   for(int i=-s; i<s+1; i++){
    for(int j=-s; j<s+1; j++){
      ip = xs + i;
      jp = ys + j;
      iq = xt + i;
      jq = yt + j;

      if(ip < 0)
	ip = Arows - ip;
      if(jp < 0)
	jp = Acols - jp;
      if(ip > Arows-1)
	ip -= Arows;
      if(jp > Acols-1)
	jp-= Acols;

      if(iq < 0)
	iq = Brows - iq;
      if(jq < 0)
	jq = Bcols - jq;
      if(iq > Brows-1)
	iq -= Brows;
      if(jq > Bcols-1)
	jq -= Bcols;

      p1 = ip*Acols+jp;
      q1 = iq*Bcols+jq;
      
      if(q1 == q)// && nb_pass == 0)
	break;

      tmp = Aprim[l].lum[p1] - Bprim[l].lum[q1];
      lum2 += tmp * tmp;
      k++;
    }
     if(q1 == q)// && nb_pass == 0)
	break;
   }

   dprim = lum2 /(float)k;
   
   if(l>0){
     xs /= 2.0;
     ys /= 2.0;
     xt /= 2.0;
     yt /= 2.0;
     s  /= 2.0;
     Acols /= 2.0;
     Arows /= 2.0;
     Bcols /= 2.0;
     Brows /= 2.0;
     k = 0;
     lum1 = 0.0;
     lum2 = 0.0;

     for(int i=-s; i<s+1; i++){
       for(int j=-s; j<s+1; j++){
	 ip = xs + i;
	 jp = ys + j;
	 iq = xt + i;
	 jq = yt + j;

	 if(ip < 0)
	   ip = Arows - ip;
	 if(jp < 0)
	   jp = Acols - jp;
	 if(ip > Arows-1)
	   ip -= Arows;
	 if(jp > Acols-1)
	   jp-= Acols;

	 if(iq < 0)
	   iq = Brows - iq;
	 if(jq < 0)
	   jq = Bcols - jq;
	 if(iq > Brows-1)
	   iq -= Brows;
	 if(jq > Bcols-1)
	   jq -= Bcols;

	 p1 = ip*Acols+jp;
	 q1 = iq*Bcols+jq;

	 tmp = A[l-1].lum[p1] - B[l-1].lum[q1];
	 lum1 += tmp * tmp;
	 tmp = Aprim[l-1].lum[p1] - Bprim[l-1].lum[q1];
	 lum2 += tmp * tmp;
	 k++;
       }
     }

     d += lum1/(float)k;
     dprim += lum2/(float)k;
   }
   //printf(" %f %f\n ",d ,dprim);
   return w*d + (1-w)*dprim;
}
/*
float
dist_transfert(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, int w)
{
  int istart, iend, jstart, jend, Acols, Arows, Bcols, Brows, p, q, p1, q1, p2, q2, k, s;
  //int Acols, Arows, Bcols, Brows, p, q, p2, q2, s;
  double d, dprim, lum = 0.0, mean = 0.0, sd = 0.0, tmp;
  Acols = A[l].cols;
  Arows = A[l].rows;
  Bcols = B[l].cols;
  Brows = B[l].rows;
  s = (int)((size-1)/2.0);
  p = xs * Acols + ys;
  q = xt * Bcols + yt;
  p2 = (xs/2.0) * (Acols/2.0) + (ys/2.0);
  q2 = (xt/2.0) * (Bcols/2.0) + (yt/2.0);
  k = 0;

  tmp = A[l].lum[p] - B[l].lum[q];
  lum = tmp * tmp;
  
  tmp = A[l].mean[p] - B[l].mean[q];
  mean = tmp * tmp;

  tmp = A[l].sd[p] - B[l].sd[q];
  sd = tmp * tmp;
  
  if(l > 0){

    tmp = A[l-1].lum[p2] - B[l-1].lum[q2];
    lum += tmp * tmp;
    
    tmp = A[l-1].mean[p2] - B[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = A[l-1].sd[p2] - B[l-1].sd[q2];
    sd += tmp * tmp;
  }

  d = lum + mean + sd;

  lum = 0.0;
  mean = 0.0;
  sd = 0.0;
  
  if(xt < s)
    istart = -s + (s - xt);
  else
    istart = -s;
  if(xt > Brows-1 - s)
    iend = s - (xt - (Brows-1 - s));
  else
    iend = s;

  if(yt < s)
    jstart = -s + (s - yt);
  else
    jstart = -s;
  if(yt > Bcols-1 - s)
    jend = s -(yt - (Bcols-1 - s));
  else
    jend = s;

  for(int i=istart; i<iend+1; i++){
    for(int j=jstart; j<jend+1; j++){

      p1 = p + (i*Acols+j);
      q1 = q + (i*Bcols+j);
           
      if(p1 == p || q1 == q)
	break;
      if(p1 < 0 ){
	p1 = Arows*Acols + p1;
	//continue;
      }
      if(p1 >= Arows*Acols)
	p1 = p1 - Arows*Acols;
      
      tmp = Aprim[l].lum[p1] - Bprim[l].lum[q1];
      lum += tmp * tmp;
      k++;
    }
    if(p1 == p || q1 == q)
      break;
  }
  
  if(k == 0)
    lum = 99999999999999999;
  
  /*tmp = Aprim[l].lum[p] - Bprim[l].lum[q];
    lum = tmp * tmp;
  //
  tmp = Aprim[l].mean[p] - Bprim[l].mean[q];
  mean = tmp * tmp;

  tmp = Aprim[l].sd[p] - Bprim[l].sd[q];
  sd = tmp * tmp;
  
  if(l > 0){ //niveau l-1
    xs /= 2.0;
    ys /= 2.0;
    xt /= 2.0;
    yt /= 2.0;
    s  /= 2.0;
    Acols /= 2.0;
    Arows /= 2.0;
    Bcols /= 2.0;
    Brows /= 2.0;
    
    if(xt < s)
      istart = -s + (s - xt);
    else
      istart = -s;
    if(xt > Brows-1 - s)
      iend = s - (xt - (Brows-1 - s));
    else
      iend = s;

    if(yt < s)
      jstart = -s + (s - yt);
    else
      jstart = -s;
    if(yt > Bcols-1 - s)
      jend = s -(yt - (Bcols-1 - s));
    else
      jend = s;

    for(int i=istart; i<iend+1; i++){
      for(int j=jstart; j<jend+1; j++){

	p1 = p2 + (i*Acols+j);
	q1 = q2 + (i*Bcols+j);

	if(p1 < 0)
	  p1 = Arows*Acols + p1;
	//continue;
	if(p1 >= Arows*Acols)
	  p1 = p1 - Arows*Acols;
	
	tmp = Aprim[l-1].lum[p1] - Bprim[l-1].lum[q1];
	lum += tmp * tmp;
	k++;
      }
    }

    /*if( k==0)
      lum =  99999999999999999;
    */
    /* tmp = Aprim[l-1].lum[p2] - Bprim[l-1].lum[q2];
       lum += tmp * tmp;
    //
    tmp = Aprim[l-1].mean[p2] - Bprim[l-1].mean[q2];
    mean += tmp * tmp;

    tmp = Aprim[l-1].sd[p2] - Bprim[l-1].sd[q2];
    sd += tmp * tmp;
  }

  dprim = lum + mean + sd;
      
  return  w * d + dprim;
}

*/
int
WeiLevoy_transfert(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, int w)
{
  float dist_min, d;
  int pmin=0, p, Arows, Acols;
  Arows = A[l].rows;
  Acols = A[l].cols;
  
  dist_min = dist_transfert(0, 0, A, Aprim, xt, yt, B, Bprim, l, size, w);
  
  for(int i=0; i<Arows; i++){
    for(int j=0; j<Acols; j++){
      p = i*Acols+j;

      d = dist_transfert(i, j, A, Aprim, xt, yt, B, Bprim, l, size, w);
      if(d <= dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}

int
Ashikhmin_transfert(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, int w)
{
  float dist_min, d;
  int pmin, p, pixel, q, Acols, Arows, Bcols, xs, ys, start, end;
  Acols = A[l].cols;
  Arows = A[l].rows;
  Bcols = B[l].cols;
  start = -(int)(size/2.0);
  end = -start;
  q = xt * Bcols + yt;

  pmin = Bprim[l].s[xt*Bcols+yt];
  xs = pmin/Acols;
  ys = pmin - xs * Acols;

  dist_min = dist_transfert(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, w);
  
  for(int i=start; i<end+1; i++){
    for(int j=start; j<end+1; j++){
      pixel = q + i*Bcols+j;

      if(pixel == q)
	break;
      if(pixel < 0)
	continue;

      p = Bprim[l].s[pixel] - i*Acols-j;

      if(p < 0 || p >= Acols*Arows )
	p = rand() % (Acols*Arows);

      xs = p/Acols;
      ys = p - xs * Acols;

      d = dist_transfert(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, w);

      if(d <= dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}

int
BestMatch_transfert(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, int w)
{
  int p_brute, p_coh, x, y;
  float dist_brute, dist_coh;
  
  p_brute = WeiLevoy_transfert(A, Aprim, xt, yt, B, Bprim, l, size, w);
  x = p_brute/A[l].cols;
  y = p_brute - x * A[l].cols;
  dist_brute = dist_transfert(x, y, A, Aprim, xt, yt, B, Bprim, l, size, w);
  
  p_coh = Ashikhmin_transfert(A, Aprim, xt, yt, B, Bprim, l, size, w);
  x = p_coh/A[l].cols;
  y = p_coh - x * A[l].cols;
  dist_coh = dist_transfert(x, y, A, Aprim, xt, yt, B, Bprim, l, size, w);

  if(dist_coh < dist_brute * (1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_brute;
}


void
transfert(char* texture_name, char* image_name, float w, int size)
{
   printf("INITIALISATION \n");
  pnm ims = pnm_load(texture_name);
  int Acols = pnm_get_width(ims);
  int Arows = pnm_get_height(ims);

  pnm imt = pnm_load(image_name);
  int Bcols = pnm_get_width(imt);
  int Brows = pnm_get_height(imt);

  pnm out = pnm_new(Bcols, Brows, PnmRawPpm);

  unsigned short* sRGB = pnm_get_image(ims);
  unsigned short* tRGB = pnm_get_image(imt);
  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * Bcols * Brows);

  float* tmp_s = convertUnsignedShort2Float(sRGB, 3 * Acols, Arows);
  float* tmp_t = convertUnsignedShort2Float(tRGB, 3 * Bcols, Brows);
  
  float* s_yiq = convertRGB2YIQ(Acols, Arows, tmp_s);
  float* t_yiq = convertRGB2YIQ(Bcols, Brows, tmp_t);
 
  float* data_source = channel(Acols, Arows, s_yiq, 0);
  float* data_source_filtre = channel(Acols, Arows, s_yiq, 0);
  float* data_target = channel(Bcols, Brows, t_yiq, 0);

  Pyramid* A = init_pyramid_transfert(Acols, Arows);
  Pyramid* Aprim = init_pyramid_transfert(Acols, Arows);
  Pyramid* B = init_pyramid_transfert(Bcols, Brows);
  Pyramid* Bprim = init_pyramid_transfert(Bcols, Brows);

  int pixel, p, q;
  float val;

  for(int l=0; l<L; l++){
    printf("LEVEL %d\n",l);

    Acols = A[l].cols;
    Arows = A[l].rows;
    Bcols = B[l].cols;
    Brows = B[l].rows;

    free(A[l].lum);
    free(Aprim[l].lum);
    free(B[l].lum);

    if( l == L-1){
      A[l].lum = data_source;
      Aprim[l].lum = data_source_filtre;
      B[l].lum = data_target;
    }
    else{
      A[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source, L-1-l);
      Aprim[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source_filtre, L-1-l);
      B[l].lum = pyramid_level(B[L-1].cols, B[L-1].rows, data_target, L-1-l);
    }
    printf("CALCUL DESCRIPTEURS\n");
    mean_img(Acols, Arows, A[l].lum, A[l].mean, 2, 0);
    mean_img(Acols, Arows, Aprim[l].lum, Aprim[l].mean, size, 1);
    mean_img(Bcols, Brows, B[l].lum, B[l].mean, 2, 0);

    sd_img(Acols, Arows, A[l].lum, A[l].mean, A[l].sd, 2, 0);
    sd_img(Acols, Arows, Aprim[l].lum, Aprim[l].mean, Aprim[l].sd, size, 1);
    sd_img(Bcols, Brows, B[l].lum, B[l].mean, B[l].sd, 2, 0);

    if(l > 0){
      int acols, arows, bcols, brows, s;
      acols = Acols/2.0;
      arows = Arows/2.0;
      bcols = Bcols/2.0;
      brows = Brows/2.0;
      s = size/2.0;
      
      mean_img(acols, arows, A[l-1].lum, A[l-1].mean, s, 0);
      mean_img(acols, arows, Aprim[l-1].lum, Aprim[l-1].mean, s, 0);
      mean_img(bcols, brows, B[l-1].lum, B[l-1].mean, s, 0);
      mean_img(bcols, brows, Bprim[l-1].lum, Bprim[l-1].mean, s, 0);

      sd_img(acols, arows, A[l-1].lum, A[l-1].mean, A[l-1].sd, s, 0);
      sd_img(acols, arows, Aprim[l-1].lum, Aprim[l-1].mean, Aprim[l-1].sd, s, 0);
      sd_img(bcols, brows, B[l-1].lum, B[l-1].mean, B[l-1].sd, s, 0);
      sd_img(bcols, brows, Bprim[l-1].lum, Bprim[l-1].mean, Bprim[l-1].sd, s, 0);
    }

    for(int i=0; i<Brows; i++){
      for(int j=0; j<Bcols; j++){
	pixel = rand() % (Acols * Arows);
	val = rand() % 256;
	Bprim[l].lum[i*Bcols+j] = val;
	Bprim[l].s[i*Bcols+j] = pixel;
      }
    }
    printf("FAIT\n");
    for(int i=0; i<Brows; i++){
      for(int j=0; j<Bcols; j++){
	q = i*Bcols+j;

	Bprim[l].mean[q] = mean(Bcols, Brows, Bprim[l].lum, q, size, 1);
	Bprim[l].sd[q] = standard_deviation(Bcols, Brows, Bprim[l].lum, Bprim[l].mean[q], q, size, 1);
	
	p = BestMatch_transfert(A, Aprim, i, j, B, Bprim, l, size, w);

	Bprim[l].lum[q] = Aprim[l].lum[p];
	Bprim[l].s[q] = p;
      }
    }
    printf("FAIT\n");
  }

  printf("CONVERSION\n");
  Bcols = Bprim[L-1].cols;
  Brows = Bprim[L-1].rows;
  /*
  float* it = channel(Bcols, Brows, t_yiq, 1);
  float* qt = channel(Bcols, Brows, t_yiq, 2);
  */
  float* isf = channel(A[L-1].cols, A[L-1].rows, s_yiq, 1);
  float* qsf = channel(A[L-1].cols, A[L-1].rows, s_yiq, 2);
  
  float* it = malloc(sizeof(float) * Bcols * Brows);
  float* qt = malloc(sizeof(float) * Bcols * Brows);

  for(int i=0; i<Brows; i++){
    for(int j=0; j<Bcols; j++){
      pixel = Bprim[L-1].s[i*Bcols+j];
      //printf("%d ",pixel);
      it[i*Bcols+j] = isf[pixel];
      qt[i*Bcols+j] = qsf[pixel];
    }
  }
  
  float* yiq = putChannel(Bcols, Brows, Bprim[L-1].lum, it, qt);

  float* rgb = convertYIQ2RGB(Bcols, Brows, yiq);

  convertFloat2UnsignedShort(rgb, data_out, 3*Bcols, Brows);

  printf("SAUVEGARDE\n");
  for(int i=0; i<Brows; i++){
    for(int j=0; j<Bcols; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(out, i, j, k, data_out[i*3*Bcols+(j*3)+k]);
      }
    }
  }

  pnm_save(out, PnmRawPpm, "transfert.ppm");

  printf("LIBERATION MEMOIRE \n");

  free(data_out);
  free(rgb);
  free(yiq);
  //free(isf);
  //free(qsf);
  free(qt);
  free(it);
  free_pyramid_transfert(A);
  free_pyramid_transfert(B);
  free_pyramid_transfert(Bprim);
  free(t_yiq);
  free(s_yiq);
  free(tmp_t);
  free(tmp_s);
  free(tRGB);
  free(sRGB);
  free(ims);
  free(imt);
  pnm_free(out);
}
