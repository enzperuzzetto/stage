#include "BestMatch.h"
#include<stdio.h>


int
bruteForceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int pmin = 0;
  int pixel = 0;
  int xt = q/target[l].cols;
  int yt = q - xt*target[l].cols;
  //float dist_min = dist(index_min, source, source_filter, q, target, target_filter, l);
  float dist_min = dist(0, 0, source, source_filter, xt, yt, target, target_filter, l);
  float d = 0.0;
  
  for(int i=0; i<source[l].rows; i++){
    for(int j=0; j<source[l].cols; j++){
      pixel = i*source[l].cols+j;
      //d = dist(pixel, source, source_filter, q, target, target_filter, l);
      d = dist(i, j, source, source_filter, xt, yt, target, target_filter, l);
      
      if(d <= dist_min){
	dist_min = d;
	pmin = pixel;
      }
    }
  }
  //printf("%d ", index_min);
  return pmin;
}

int
bestCoherenceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  //srand(time(NULL));
  int pixel, pmin = target_filter[l].s[q], p;
  int xt = q/target[l].cols;
  int yt = q - xt*target[l].cols;
  int xs = pmin/source_filter[l].cols;
  int ys = pmin - xs * source_filter[l].cols;
  int start = -(int)(NEIGHBOOR_SIZE_FINER/2.0);
  int end = (int)(NEIGHBOOR_SIZE_FINER/2.0);
  //float dist_min = dist(pmin, source, source_filter, q, target, target_filter, l);
  float dist_min = dist(xs, ys, source, source_filter, xt, yt, target, target_filter, l);
  
  float d = 0.0;
  for(int i=start; i<end+1; i++){
    for(int j=start; j<end+1; j++){

      pixel = q + (i*target_filter[l].cols+j);

      if(pixel == q )
	break;
      if(pixel < 0)
	continue;
      
      p = target_filter[l].s[pixel] -i*source_filter[l].cols-j;
      
      if(p < 0 || p > source[l].cols * source[l].rows -1){
	p =  rand()%(source_filter[l].cols*source_filter[l].rows);
	//printf(" %d ", p);
      }
      
      //d = dist(p, source, source_filter, q, target, target_filter, l);
      d = dist(xs, ys, source, source_filter, xt, yt, target, target_filter, l);
   
      if(d < dist_min){
	pmin = p;
	dist_min = d;
      }
    }
     if(pixel == q )
	break;
  }

  return pmin;
}

int
bestMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  /*#if 0
  int p_patch = patchMatch(source, source_filter, target, target_filter, l, q);
  int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);
  
  float dist_patch = dist(p_patch, source, source_filter, q, target, target_filter, l);
  dist_patch *= dist_patch;
  float dist_coh = dist(p_coh, source, source_filter, q, target, target_filter, l);
  dist_coh *= dist_coh;
  
  if(dist_coh < dist_patch * ( 1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_patch;
#endif
  */
#if 1
  int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
  int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);

  int xs = p_brute/source[l].cols;
  int ys = p_brute - xs * source[l].cols;
  int xt = q/target[l].cols;
  int yt = q - xt * target[l].cols;
  float dist_brute = dist(xs, ys, source, source_filter, xt, yt, target, target_filter, l);
  //dist_brute *= dist_brute;

  xs = p_coh/source[l].cols;
  ys = p_coh - xs * source[l].cols;
  float dist_coh = dist(xs, ys, source, source_filter, xt, yt, target, target_filter, l);
  //dist_coh *= dist_coh;

 
  if(dist_coh < dist_brute * ( 1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_brute;
#endif
  /*
#if 0
   int p_patch = patchMatch(source, source_filter, target, target_filter, l, q);
   return p_patch;
#endif

#if 0
   int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);
   return p_coh;
#endif

#if 0
   int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
   return p_brute;
   #endif*/
}
