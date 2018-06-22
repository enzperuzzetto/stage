#include "synthese.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bcl.h"
#include "pyramid.h"
#include "conversion.h"


typedef struct Pyramid
{
  int cols;
  int rows;
  float* lum;
  int* s;
}Pyramid;


Pyramid*
init_pyramid(int cols, int rows)
{
  Pyramid* pyramid = malloc(sizeof(Pyramid) * L);

  for(int l=0; l<L; l++){
    int* dim = pyramid_level_dim(cols, rows, L-1-l);
    pyramid[l].cols = dim[0];
    pyramid[l].rows = dim[1];
    pyramid[l].lum = malloc(sizeof(float) * dim[0] * dim[1]);
    pyramid[l].s = malloc(sizeof(int) * dim[0] * dim[1]);

     free(dim);
  }

  return pyramid;
}

void
free_pyramid(Pyramid* pyramid)
{
  for(int l=0; l<L; l++){
    free(pyramid[l].lum);
    free(pyramid[l].s);
  }

  free(pyramid);
}

//###############################################################################################################################
float
dist(int xs, int ys, Pyramid* texture, int xt, int yt, Pyramid* result, int l, int size, int nbpass)
{
  int s = (int)(size/2.0);
  int Acols, Arows, Bcols, Brows, ip, jp, iq, jq, p1, q1, q, k=0;
  float lum1 = 0.0, lum2=0.0, tmp;

  Acols = texture[l].cols;
  Arows = texture[l].rows;
  Bcols = result[l].cols;
  Brows = result[l].rows;
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
	jp -= Acols;

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

      if((q1 == q) && nbpass == 0)
	break;
      
      tmp = texture[l].lum[p1] - result[l].lum[q1];
      lum1 += tmp * tmp;
      k++;
    }
    if((q1 == q) && nbpass == 0)
	break;
  }

  //lum1 /= (float)k;
   
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

	tmp = texture[l-1].lum[p1] - result[l-1].lum[q1];
	lum2 += tmp * tmp;
	k++;
      }
    }

    //lum2 /= (float)k;
  }

  return lum1 + lum2;
}
	 
/*
float
dist(int xs, int ys, Pyramid* texture, int xt, int yt, Pyramid* result, int l, int size)
{
  int istart, iend, jstart, jend, ctex, rtex, cout, rout, p, q, p1, q1, k, s;
  double lum = 0.0, tmp;
  ctex = texture[l].cols;
  rtex = texture[l].rows;
  cout = result[l].cols;
  rout = result[l].rows;
  s = (int)((size-1)/2.0);
  p = xs * ctex + ys;
  q = xt * cout + yt;
  k = 0;
  
 
  if(xt < s)
    istart = -s + (s - xt);
  else
    istart = -s;
  if(xt > rout-1 - s)
    iend = s - (xt - (rout-1 - s));
  else
    iend = s;

  if(yt < s)
    jstart = -s + (s - yt);
  else
    jstart = -s;
  if(yt > cout-1 - s)
    jend = s -(yt - (cout-1 - s));
  else
    jend = s;

  for(int i=istart; i<iend+1; i++){
    for(int j=jstart; j<jend+1; j++){

      p1 = p + (i*ctex+j);
      q1 = q + (i*cout+j);
           
      if(p1 == p || q1 == q)
	break;

      if(p1 < 0)
	continue;
      
      if(p1 < 0 ){
	//p1 = -p1;
	
	p1 = ctex*rtex + p1;
	  //printf("%d ",p1);
      }
      if(p1 >= ctex*rtex)
	p1 = p1 - ctex*rtex;
      tmp = texture[l].lum[p1] - result[l].lum[q1];
      lum += tmp * tmp;
      k++;
    }
    if(p1 == p || q1 == q)
      break;
  }
  
  if(k == 0)
    lum = 99999999999999999;

  if(l > 0){ //niveau l-1
    xs /= 2.0;
    ys /= 2.0;
    xt /= 2.0;
    yt /= 2.0;
    s  /= 2.0;
    ctex /= 2.0;
    rtex /= 2.0;
    cout /= 2.0;
    rout /= 2.0;

    p = xs * ctex + ys;
    q = xt * cout + yt;

    if(xt < s)
      istart = -s + (s - xt);
    else
      istart = -s;
    if(xt > rout-1 - s)
      iend = s - (xt - (rout-1 - s));
    else
      iend = s;

    if(yt < s)
      jstart = -s + (s - yt);
    else
      jstart = -s;
    if(yt > cout-1 - s)
      jend = s -(yt - (cout-1 - s));
    else
      jend = s;

    for(int i=istart; i<iend+1; i++){
      for(int j=jstart; j<jend+1; j++){

	p1 = p + (i*ctex+j);
	q1 = q + (i*cout+j);
	
	if(p1 < 0 ){
	  //p1 = p1%(ctex*rtex);
	  //p1 = -p1;
	  p1 = ctex*rtex + p1;
	  // printf("%d ", p1);
	}
	if(p1 >= ctex*rtex)
	  p1 = p1 - ctex*rtex;
	
	
	tmp = texture[l-1].lum[p1] - result[l-1].lum[q1];
	lum += tmp * tmp;
      }
    }
  }
  
      
  return sqrtf(lum);
}
*/
int
WeiLevoy(Pyramid* texture, int xt, int yt, Pyramid* result, int l, int size, int nbpass)
{
  float dist_min, d;
  int pmin=0, p, rtex, ctex;
  rtex = texture[l].rows;
  ctex = texture[l].cols;
  
  dist_min = dist(0, 0, texture, xt, yt, result, l, size, nbpass);
  
  for(int i=0; i<rtex; i++){
    for(int j=0; j<ctex; j++){
      p = i*ctex+j;

      d = dist(i, j, texture, xt, yt, result, l, size, nbpass);
      if(d <= dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}

int
Ashikhmin(Pyramid* texture, int xt, int yt, Pyramid* result, int l, int size, int nbpass)
{
  float dist_min, d;
  int pmin, p, pixel, q, ctex,rtex, cout, xs, ys, start, end;
  ctex = texture[l].cols;
  rtex = texture[l].rows;
  cout = result[l].cols;
  start = -(int)(size/2.0);
  end = -start;
  q = xt * cout + yt;

  pmin = result[l].s[xt*cout+yt];
  xs = pmin/ctex;
  ys = pmin - xs * ctex;

  dist_min = dist(xs, ys, texture, xt, yt, result, l, size, nbpass);
  
  for(int i=start; i<end+1; i++){
    for(int j=start; j<end+1; j++){
      pixel = q + i*cout+j;

      if(pixel == q)
	break;
      if(pixel < 0)
	continue;

      p = result[l].s[pixel] - i*ctex-j;

      if(p < 0 || p >= ctex*rtex )
	p = rand() % (ctex*rtex);

      xs = p/ctex;
      ys = p - xs * ctex;

      d = dist(xs, ys, texture, xt, yt, result, l, size,nbpass);

      if(d <= dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}
	
	


int
BestMatch(Pyramid* texture, int xt, int yt, Pyramid* result, int l, int size,int nbpass)
{
  int p_brute, p_coh, x, y;
  float dist_brute, dist_coh;
  /*
  p_brute = WeiLevoy(texture, xt, yt, result, l, size, nbpass);
  x = p_brute/texture[l].cols;
  y = p_brute - x * texture[l].cols;
  dist_brute = dist(x, y, texture, xt, yt, result, l, size,nbpass);
  */
  p_coh = Ashikhmin(texture, xt, yt, result, l, size,nbpass);
  x = p_coh/texture[l].cols;
  y = p_coh - x * texture[l].cols;
  dist_coh = dist(x, y, texture, xt, yt, result, l, size,nbpass);
  
  // if(dist_coh < dist_brute * (1 + powf(2, l-L)*K) )
  return p_coh;
    //else
  //return p_brute;
  
}





void
synthese(char* texture_name, int cols, int rows, int size)
{
  printf("INITIALISATION \n");
  
  pnm tex = pnm_load(texture_name);
  int tex_cols = pnm_get_width(tex);
  int tex_rows = pnm_get_height(tex);

  pnm imt = pnm_new(cols, rows, PnmRawPpm);

  unsigned short* RGBtex = pnm_get_image(tex);
  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * cols * rows);

  float* tmp =  convertUnsignedShort2Float(RGBtex, 3 * tex_cols, tex_rows);
  float* YIQtex = convertRGB2YIQ(tex_cols, tex_rows, tmp);

  float* lumtex = channel(tex_cols, tex_rows, YIQtex, 0);
  float* itex = channel(tex_cols, tex_rows, YIQtex, 1);
  float* qtex = channel(tex_cols, tex_rows, YIQtex, 2);

  Pyramid* texture = init_pyramid(tex_cols, tex_rows);
  Pyramid* result  = init_pyramid(cols, rows);

  int ctex, rtex, cout, rout, pixel, q, p;
  float val;

  for(int l=0; l<L; l++){
    printf("LEVEL %d\n",l);
    
    ctex = texture[l].cols;
    rtex = texture[l].rows;
    cout = result[l].cols;
    rout = result[l].rows;

    free(texture[l].lum);
    if(l == L-1)
      texture[l].lum = lumtex;
    else
      texture[l].lum = pyramid_level(tex_cols, tex_rows, lumtex, L-1-l);

    for(int i=0; i<rout; i++){
      for(int j=0; j<cout; j++){
	pixel = rand() % (ctex * rtex);
	val = rand() % 256;
	result[l].lum[i*cout+j] = val;
	result[l].s[i*cout+j] = pixel;
      }
    }

    for(int nbpass = 0; nbpass<2; nbpass++){
      for(int i=0; i<rout; i++){
	for(int j=0; j<cout; j++){
	  q = i*cout+j;

	  p = BestMatch(texture, i, j, result, l, size, nbpass);

	  result[l].lum[q] = texture[l].lum[p];
	  result[l].s[q] = p;
	}
      }
    }
    printf("FAIT \n");
  }
  printf("CONVERSION \n");
  float* it = malloc(sizeof(float) * cols * rows);
  float* qt = malloc(sizeof(float) * cols * rows);

  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      pixel = result[L-1].s[i*cols+j];
      it[i*cols+j] = itex[pixel];
      qt[i*cols+j] = qtex[pixel];
    }
  }

  float* yiq = putChannel(cols, rows, result[L-1].lum, it, qt);
  float* rgb =  convertYIQ2RGB(cols, rows, yiq);

  convertFloat2UnsignedShort(rgb, data_out, 3*cols, rows);

  printf("SAUVEGARDE \n");
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(imt, i, j, k, data_out[i*3*cols+(j*3)+k]);
      }
    }
  }

  pnm_save(imt, PnmRawPpm, "syn.ppm");

  printf("LIBERATION MEMOIRE \n");
  
  free(data_out);
  free(rgb);
  free(yiq);
  free(qt);
  free(it);
  free_pyramid(result);
  free_pyramid(texture);
  free(qtex);
  free(itex);
  free(YIQtex);
  free(tmp);
  free(RGBtex);
  pnm_free(imt);
  free(tex);
  
}
