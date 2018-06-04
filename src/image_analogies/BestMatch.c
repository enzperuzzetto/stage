#include "BestMatch.h"
#include<stdio.h>


int
bruteForceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int index_min = 0;
  int pixel = 0;
  float dist_min = dist(index_min, source, source_filter, q, target, target_filter, l);  
  
  float d = 0.0;
  for(int i=0; i<source[l].rows; i++){
    for(int j=0; j<source[l].cols; j++){
      pixel = i*source[l].cols+j;
      d = dist(pixel, source, source_filter, q, target, target_filter, l);
    
      if(d < dist_min){
	dist_min = d;
	index_min = pixel;
      }
    }
  }
  printf("%d ", index_min);
  return index_min;
}

int
bestCoherenceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  //srand(time(NULL));
  int pixel, pmin = target_filter[l].s[q], p;
  int i_p1, j_p1;
  float dist_min = dist(pmin, source, source_filter, q, target, target_filter, l);
  
  float d = 0.0;
  for(int i=0; i<NEIGHBOOR_SIZE_FINER; i++){
    for(int j=0; j<NEIGHBOOR_SIZE_FINER; j++){
      i_p1 = i - NEIGHBOOR_SIZE_FINER/2.0;
      j_p1 = j - NEIGHBOOR_SIZE_FINER/2.0;
      pixel = q + (i_p1*target_filter[l].cols+j_p1);

      if(pixel == q || pixel < 0)
	break;
      
      p = target_filter[l].s[pixel] -i_p1*source_filter[l].cols-j_p1;
      
      if(p < 0 || p > source[l].cols * source[l].rows -1){
	p =  rand()%(source_filter[l].cols*source_filter[l].rows);
	//printf(" %d ", p);
      }
      
      d = dist(p, source, source_filter, q, target, target_filter, l);
   
      if(d < dist_min){
	pmin = p;
	dist_min = d;
      }
    }
  }

  return pmin;
}

int
bestMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
#if 0
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
  
#if 0
  int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
  int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);

  float dist_brute = dist(p_brute, source, source_filter, q, target, target_filter, l);
  dist_brute *= dist_brute;
  float dist_coh = dist(p_coh, source, source_filter, q, target, target_filter, l);
  dist_coh *= dist_coh;

 
  if(dist_coh < dist_brute * ( 1 + powf(2, l-L)*K) )
    return p_coh;
  else
    return p_brute;
#endif

#if 0
   int p_patch = patchMatch(source, source_filter, target, target_filter, l, q);
   return p_patch;
#endif

#if 1
   int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);
   return p_coh;
#endif

#if 0
   int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
   return p_brute;
#endif
}
