#include "filtre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
init_pyramid_filtre(int cols, int rows)
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
free_pyramid_filtre(Pyramid* pyramid)
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
/*float
dist_filtre(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
{
   int s = (int)(size/2.0);
   int Acols, Arows, Bcols, Brows, ip, jp, iq, jq, p1, q1, q, k=0;
   float lum1 = 0.0, lum2 = 0.0, d = 0.0, dprim = 0.0, tmp;

   Acols = A[l].cols;
   Arows = A[l].rows;
   Bcols = B[l].cols;
   Brows = B[l].rows;
   q = xt*Bcols+yt;

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
      
      if((q1 == q) && nb_pass == 0)
	break;
      
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

   return d + dprim;
}
*/	 
     
float
dist_filtre(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
{
  int istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,Acols, Arows, Bcols, Brows, p, q, p1=0, q1=0, p2, q2, k, s;
  double d, dprim, lum1 = 0.0, lum2=0.0, tmp;

  Acols = A[l].cols;
  Arows = A[l].rows;
  Bcols = B[l].cols;
  Brows = B[l].rows;
  s = (int)((size+1)/2.0);
  p = xs * Acols + ys;
  q = xt * Bcols + yt;
  p2 = (xs/2.0) * (Acols/2.0) + (ys/2.0);
  q2 = (xt/2.0) * (Bcols/2.0) + (yt/2.0);
  k = 0;

  if(xt < s)
    istart1 = -s + (s - xt);
  else
    istart1 = -s;
  if(xt > Brows-1 - s)
    iend1 = s - (xt - (Brows-1 - s));
  else
    iend1 = s;

  if(yt < s)
    jstart1 = -s + (s - yt);
  else
    jstart1 = -s;
  if(yt > Bcols-1 - s)
    jend1 = s -(yt - (Bcols-1 - s));
  else
    jend1 = s;

  for(int i=istart1; i<iend1+1; i++){
    for(int j=jstart1; j<jend1+1; j++){

      p1 = p + (i*Acols+j);
      q1 = q + (i*Bcols+j);
           
      if(p1 < 0)
	p1 = Acols*Arows + p1;
      if( p1 >= Acols*Arows)
	p1 = p1 - Arows*Acols;
     
      tmp = A[l].lum[p1] - B[l].lum[q1];
      lum1 += tmp * tmp;
      k++;
    }
  }
  
  if(k == 0)
    lum1 = 99999999999999999;
 
  d = lum1 / (float)k;

  for(int i=istart1; i<iend1+1; i++){
    for(int j=jstart1; j<jend1+1; j++){

      p1 = p + (i*Acols+j);
      q1 = q + (i*Bcols+j);
           
      if(p1 == p || q1 == q)
	break;
      if(p1 < 0)
	p1 = Acols*Arows + p1;
      if( p1 >= Acols*Arows)
	p1 = p1 - Arows*Acols;
     
      tmp = Aprim[l].lum[p1] - Bprim[l].lum[q1];
      lum2 += tmp * tmp;
      k++;
    }
    if(p1 == p || q1 == q)
      break;
  }
  
  if(k == 0)
    lum2 = 99999999999999999;
 
  dprim = lum2 / (float)k;
  
  if(l > 0){
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
    
    if(xt < s)
      istart2 = -s + (s - xt);
    else
      istart2 = -s;
    if(xt > Brows-1 - s)
      iend2 = s - (xt - (Brows-1 - s));
    else
      iend2 = s;

    if(yt < s)
      jstart2 = -s + (s - yt);
    else
      jstart2 = -s;
    if(yt > Bcols-1 - s)
      jend2 = s -(yt - (Bcols-1 - s));
    else
      jend2 = s;

    for(int i=istart2; i<iend2+1; i++){
      for(int j=jstart2; j<jend2+1; j++){

	p1 = p2 + (i*Acols+j);
	q1 = q2 + (i*Bcols+j);

	if(p1 < 0)
	  p1 = Acols*Arows + p1;
	if( p1 >= Acols*Arows)
	  p1 = p1 - Arows*Acols;
	
	tmp = A[l-1].lum[p1] - B[l-1].lum[q1];
	lum1 += tmp * tmp;
	tmp = Aprim[l-1].lum[p1] - Bprim[l-1].lum[q1];
	lum2 += tmp * tmp;
	k++;
      }
    }

    if( k==0){
      lum1 =  99999999999999999;
      lum2 = 99999999999999999;
    }
  
    d += lum1 / (float)k;
    dprim += lum2 / (float)k;
  }

  return d + dprim;
}

/*
float
dist_filtre(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
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
      if(p1 < 0)
	p1 = Acols*Arows + p1;
      if( p1 >= Acols*Arows)
	p1 = p1 - Arows*Acols;
      //continue;
     
      tmp = Aprim[l].lum[p1] - Bprim[l].lum[q1];
      lum += tmp * tmp;
      k++;
    }
    if(p1 == p || q1 == q)
      break;
  }
  
  if(k == 0)
    lum = 99999999999999999;
  
  tmp = Aprim[l].lum[p] - Bprim[l].lum[q];
    lum = tmp * tmp;
  
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
    k = 0;
    p = xs * Acols + ys;
    q = xt * Bcols + yt;
    
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
	
	if(p1 < 0 || p1 >= Acols*Arows)
	  continue;
	
	if(p1 < 0)
	  p1 = Acols*Arows + p1;
	if( p1 >= Acols*Arows)
	  p1 = p1 - Arows*Acols;
	
	tmp = Aprim[l-1].lum[p1] - Bprim[l-1].lum[q1];
	lum += tmp * tmp;
	k++;
      }
    }

    if( k==0)
      lum =  99999999999999999;

    tmp = Aprim[l-1].lum[p2] - Bprim[l-1].lum[q2];
       lum += tmp * tmp;

    tmp = Aprim[l-1].mean[p] - Bprim[l-1].mean[q];
    mean += tmp * tmp;

    tmp = Aprim[l-1].sd[p] - Bprim[l-1].sd[q];
    sd += tmp * tmp;
  }

  dprim = lum + mean + sd;
      
  return d + dprim;
}
*/

int
WeiLevoy_filtre(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
{
  float dist_min, d;
  int pmin=0, p, Arows, Acols;
  Arows = A[l].rows;
  Acols = A[l].cols;
  
  dist_min = dist_filtre(0, 0, A, Aprim, xt, yt, B, Bprim, l, size);
  
  for(int i=0; i<Arows; i++){
    for(int j=0; j<Acols; j++){
      p = i*Acols+j;

      d = dist_filtre(i, j, A, Aprim, xt, yt, B, Bprim, l, size);
      if(d < dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}

int
Ashikhmin_filtre(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
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

  dist_min = dist_filtre(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size);
  
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

      d = dist_filtre(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size);

      if(d < dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}


int
BestMatch_filtre(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size)
{
  int p_brute, p_coh, x, y;
  float dist_brute, dist_coh;
  
  p_brute = WeiLevoy_filtre(A, Aprim, xt, yt, B, Bprim, l, size);
  x = p_brute/A[l].cols;
  y = p_brute - x * A[l].cols;
  dist_brute = dist_filtre(x, y, A, Aprim, xt, yt, B, Bprim, l, size);
  
  p_coh = Ashikhmin_filtre(A, Aprim, xt, yt, B, Bprim, l, size);
  x = p_coh/A[l].cols;
  y = p_coh - x * A[l].cols;
  dist_coh = dist_filtre(x, y, A, Aprim, xt, yt, B, Bprim, l, size);

  
  if(dist_coh < dist_brute * (1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_brute;
  
}



void
filtre(char* source_name, char* source_filtre_name, char* target_name, int size)
{
  printf("INITIALISATION \n");
  pnm ims = pnm_load(source_name);
  pnm ims_filtre = pnm_load(source_filtre_name);
  int Acols = pnm_get_width(ims);
  int Arows = pnm_get_height(ims);

  pnm imt = pnm_load(target_name);
  int Bcols = pnm_get_width(imt);
  int Brows = pnm_get_height(imt);

  pnm out = pnm_new(Bcols, Brows, PnmRawPpm);

  unsigned short* sRGB = pnm_get_image(ims);
  unsigned short* s_fRGB = pnm_get_image(ims_filtre);
  unsigned short* tRGB = pnm_get_image(imt);
  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * Bcols * Brows);

  float* tmp_s = convertUnsignedShort2Float(sRGB, 3 * Acols, Arows);
  float* tmp_sf = convertUnsignedShort2Float(s_fRGB, 3 * Acols, Arows);
  float* tmp_t = convertUnsignedShort2Float(tRGB, 3 * Bcols, Brows);
  
  float* s_yiq = convertRGB2YIQ(Acols, Arows, tmp_s);
  float* sf_yiq = convertRGB2YIQ(Acols, Arows, tmp_sf);
  float* t_yiq = convertRGB2YIQ(Bcols, Brows, tmp_t);
 
  float* data_source = channel(Acols, Arows, s_yiq, 0);
  float* data_source_filter = channel(Acols, Arows, sf_yiq, 0);
  float* data_target = channel(Bcols, Brows, t_yiq, 0);

  Pyramid* A = init_pyramid_filtre(Acols, Arows);
  Pyramid* Aprim = init_pyramid_filtre(Acols, Arows);
  Pyramid* B = init_pyramid_filtre(Bcols, Brows);
  Pyramid* Bprim = init_pyramid_filtre(Bcols, Brows);

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
      Aprim[l].lum = data_source_filter;
      B[l].lum = data_target;
    }
    else{
      A[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source, L-1-l);
      Aprim[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source_filter, L-1-l);
      B[l].lum = pyramid_level(B[L-1].cols, B[L-1].rows, data_target, L-1-l);
    }
    /* printf("CALCUL DESCRIPTEURS\n");
    mean_img(Acols, Arows, A[l].lum, A[l].mean, size, 0);
    mean_img(Acols, Arows, Aprim[l].lum, Aprim[l].mean, size, 1);
    mean_img(Bcols, Brows, B[l].lum, B[l].mean, size, 0);

    sd_img(Acols, Arows, A[l].lum, A[l].mean, A[l].sd, size, 0);
    sd_img(Acols, Arows, Aprim[l].lum, Aprim[l].mean, Aprim[l].sd, size, 1);
    sd_img(Bcols, Brows, B[l].lum, B[l].mean, B[l].sd, size, 0);

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
      }*/

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

	/*	Bprim[l].mean[q] = mean(Bcols, Brows, Bprim[l].lum, q, size, 1);
	Bprim[l].sd[q] = standard_deviation(Bcols, Brows, Bprim[l].lum, Bprim[l].mean[q], q, size, 1);
	*/
	p = BestMatch_filtre(A, Aprim, i, j, B, Bprim, l, size);

	Bprim[l].lum[q] = B[l].lum[q] + (Aprim[l].lum[p] - A[l].lum[p]);
	Bprim[l].s[q] = p;
      }
    }
    printf("FAIT\n");
  }
  
  printf("CONVERSION\n");
  Bcols = Bprim[L-1].cols;
  Brows = Bprim[L-1].rows;
#if 1
  float* it = channel(Bcols, Brows, t_yiq, 1);
  float* qt = channel(Bcols, Brows, t_yiq, 2);
#endif
#if 0
  float* isf = channel(A[L-1].cols, A[L-1].rows, sf_yiq, 1);
  float* qsf = channel(A[L-1].cols, A[L-1].rows, sf_yiq, 2);
  
  float* it = malloc(sizeof(float) * Bcols * Brows);
  float* qt = malloc(sizeof(float) * Bcols * Brows);

  for(int i=0; i<Brows; i++){
    for(int j=0; j<Bcols; j++){
      pixel = Bprim[L-1].s[i*Bcols+j];
      it[i*Bcols+j] = isf[pixel];
      qt[i*Bcols+j] = qsf[pixel];
    }
  }
#endif
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

  pnm_save(out, PnmRawPpm, "filtre.ppm");
  
  printf("LIBERATION MEMOIRE \n");

  free(data_out);
  free(rgb);
  free(yiq);
  free(qt);
  free(it);
  free_pyramid_filtre(A);
  free_pyramid_filtre(Aprim);
  free_pyramid_filtre(B);
  free_pyramid_filtre(Bprim);
  free(t_yiq);
  free(sf_yiq);
  free(s_yiq);
  free(tmp_t);
  free(tmp_sf);
  free(tmp_s);
  free(tRGB);
  free(s_fRGB);
  free(sRGB);
  free(ims);
  free(ims_filtre);
  free(imt);
  pnm_free(out);
}
