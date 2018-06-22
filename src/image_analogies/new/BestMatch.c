#include "BestMatch.h"

float**
initKernel(int size)
{
    float** kernel = malloc(sizeof(float*) * size+1);
    for(int i=0; i<=size+1; i++){
      kernel[i] = malloc(sizeof(float)*size);
    }

    for(int width=1; width<=size; width+=2){
      float denom = 1<<(width-1);
      kernel[width][0] = 1/denom;
      printf("kernel[%d]: %f ",width, kernel[width][0]);
      for(int i=1; i<width; i++){
        kernel[width][i] = kernel[width][i-1] * (width-i) /i;
        printf("%f ",kernel[width][i]);
      }
      printf("\n");
    }
    printf("\n");
    return kernel;
}

void
freeKernel(float** kernel, int size)
{
  for(int i=0; i<size+1; i++)
    free(kernel[i]);
  free(kernel);
}


float
dist_level(int xs, int ys, Pyramid* A, int xt, int yt, Pyramid* B, int l, int size, float W, int weighting, int skipCenter, float** kernel)
{
  int offset = size/2, dstX, dstY, srcX, srcY;
  float totalweight = 0., d = 0., tmp;

  for(int i=0; i<size; i++){
    srcX = xs + i - offset;
    dstX = xt + i - offset;

    for(int j=0; j<size; j++){
      srcY = ys + j - offset;
      dstY = yt + j - offset;

      if( skipCenter && dstY == yt && dstX == xt)
        break;

      float weight = kernel[size][i] * kernel[size][j];

      //B' no valid
      if(dstX*B[l].cols+dstY < 0 || dstX*B[l].cols+dstY > B[l].cols*B[l].rows-1)
        continue;


      if(srcX < 0)
        srcX = A[l].rows + srcX;
      else if(srcX >= A[l].rows)
        srcX -= A[l].rows;
      if(srcY < 0)
        srcY = A[l].cols + srcX;
      else if(srcY >= A[l].cols)
        srcY -= A[l].cols;

      //if(srcX <0 || srcX >= A[l].rows || srcY<0 || srcY>= A[l].rows)
        //continue;

      tmp = A[l].lum[srcX*A[l].cols+srcY] - B[l].lum[dstX*B[l].cols+dstY];
      float dp = tmp * tmp;

      if( weighting ){
        totalweight += weight * W;
        weight *= W;
      }
      else{
        totalweight += weight;
      }

      d += weight * weight * dp;

      //totalweight += weight;
    }
    if( skipCenter && dstY == yt && dstX == xt)
      break;
  }

  if(totalweight > 0)
    d /= (totalweight * totalweight);

  return d;
}

float
dist(int xs, int ys, Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel)
{
  float lweight = 1;
  float d = 0., d2 = 0.;

  //lweight *= levelWeight * levelWeight;

  d2 = dist_level(xs, ys, Aprim, xt, yt, Bprim, l, size, W, 0, 1, kernel);

  if( W > 0.0 )//&& W < 1.0)
    d2 += dist_level(xs, ys, A, xt, yt, B, l, size, W, 1, 0, kernel);
  /*else if( W >= 1.0 )
    d2 += dist_level(xs, ys, A, xt, yt, B, l, size, W, 1, 0, kernel);
*/
  d2 *= lweight * lweight;
  d += d2;

  if(l > 0){
    d2 = 0.;
    lweight *= levelWeight * levelWeight;

    d2 = dist_level(xs/2, ys/2, Aprim, xt/2, yt/2, Bprim, l-1, size/2+1, W, 0, 0, kernel);

    if( W > 0.0)
      d2 += dist_level(xs/2, ys/2, A, xt/2, yt/2, B, l-1, size/2+1, W, 1, 0, kernel);

    d2 *= lweight * lweight;

    d += d2;
  }

  return d;
}

int
WeiLevoy(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel)
{
  float dist_min, d;
  int pmin=0, p, Arows, Acols;
  Arows = A[l].rows;
  Acols = A[l].cols;

  dist_min = dist(0, 0, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);

  for(int i=0; i<Arows; i++){
    for(int j=0; j<Acols; j++){
      p = i*Acols+j;

      d = dist(i, j, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);
      if(d < dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}

int
Ashikhmin(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel)
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

  dist_min = dist(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);

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

      d = dist(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);

      if(d < dist_min){
	dist_min = d;
	pmin = p;
      }

    }
  }

  return pmin;
}





int
BestMatch(Pyramid* A, Pyramid* Aprim, int xt, int yt, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, int L, int K, float** kernel)
{
  int p_brute, p_coh, x, y;
  float dist_brute, dist_coh;

  p_brute = WeiLevoy(A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);
  x = p_brute/A[l].cols;
  y = p_brute - x * A[l].cols;
  dist_brute = dist(x, y, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);

  p_coh = Ashikhmin(A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);
  x = p_coh/A[l].cols;
  y = p_coh - x * A[l].cols;
  dist_coh = dist(x, y, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, kernel);


  if(dist_coh < dist_brute * (1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_brute;

}
