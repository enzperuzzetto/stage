#include <BestMatch.h>


int
bruteForceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int index_min = 0;
  int pixel = 0;
  float dist_min = dist(0, source, source_filter, q, target, target_filter, l);
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
  
  return index_min;
}

int
bestCoherenceMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  int pixel = 0, p = 0;
  int r = 0, i_p1 = 0, j_p1 = 0;
  float dist_min = dist(0, source, source_filter, q, target, target_filter, l);
  float d = 0.0;
  for(int i=0; i<NEIGHBOOR_SIZE_FINER; i++){
    for(int j=0; j<NEIGHBOOR_SIZE_FINER; j++){
      i_p1 = i - NEIGHBOOR_SIZE_FINER/2.0;
      j_p1 = j - NEIGHBOOR_SIZE_FINER/2.0;
      pixel = q + (i_p1*target_filter[l].cols+j_p1);

      if(pixel == q)
	break;
      if( pixel < 0 )
	pixel = 0;
      else if( pixel > target_filter[l].cols * target_filter[l].rows -1 )
	pixel = source[l].cols * source[l].rows -1;

      
      p = target_filter[l].s[pixel] + (q-pixel);
      
      if(p < 0)
	p = 0;
      else if(p > source[l].cols * source[l].rows -1)
	p = source[l].cols * source[l].rows -1;
      
      d = dist(p, source, source_filter, q, target, target_filter, l);
   
      if(d < dist_min){
	r = pixel;
	dist_min = d;
      }
    }
  }

  p = target_filter[l].s[r] + (q - r);
  
  if(p < 0)
    p = 0;
  else if(p > source[l].cols * source[l].rows -1)
    p = source[l].cols * source[l].rows -1;
  
  return p;
}

int
bestMatch(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q)
{
  //int p_brute = bruteForceMatch(source, source_filter, target, target_filter, l, q);
  int p_patch = patchMatch(source, source_filter, target, target_filter, l, q);
  int p_coh = bestCoherenceMatch(source, source_filter, target, target_filter, l, q);
  

  //float dist_brute = dist(p_brute, source, source_filter, q, target, target_filter, l);
  float dist_patch = dist(p_patch, source, source_filter, q, target, target_filter, l);
  dist_patch *= dist_patch;
  float dist_coh = dist(p_coh, source, source_filter, q, target, target_filter, l);
  dist_coh *= dist_coh;
  
  if(dist_coh < dist_patch * ( 1 + powf(2, l-L)*K) )
    return p_coh;
  else
  return p_patch;
  //return p_brute;
}
