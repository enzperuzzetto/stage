#include <patchMatch.h>

int
initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
{
  int* nnf = malloc(sizeof(int) * cols_t * rows_t);
  srand(time(NULL));
  for(int i=0; i<cols_t*rows_t; i++){
    nnf[i] =  rand()%(cols_s*rows_s);
  }

  return nnf;
}

int
iteration(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
   int pixelTarget = propagation(source, source_filter, target, target_filter, l, q, nnf);
  int pixelSource = search(source, source_filter, target, target_filter, l, pixelTarget, nnf);

  return pixelSource;
}

int
search(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
   int w = max(source[l].cols, source[l].rows), pixel_x = 0, pixel_y = 0, pixel = 0, pixel_min = nnf[q], i = 0;
  float dist_min = dist(pixel_min, source, source_filter, q, target, target_filter, l), d = 0.0, alpha = powf(ALPHA, i), x = 0.0, y = 0.0;

  while( w*alpha > 1){
    x = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;
    y = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;

    pixel_x = w*alpha*x;
    pixel_y = w*alpha*y;
    pixel = pixel_min + (pixel_x*source[l].cols + pixel_y);

    if(pixel < 0)
      pixel = 0;
    else if(pixel > source[l].cols*source[l].rows-1)
      pixel = source[l].cols*source[l].rows-1;

    d = dist(pixel, source, source_filter, q, target, target_filter, l);

    if( d < dist_min){
      dist_min = d;
      pixel_min = pixel;
    }

    i++;
    alpha = powf(ALPHA, i);
  }

  return pixel_min;
}

int
propagation(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf)
{
   int pixel_min = q;
  float dist_min = dist(nnf[pixel_min], source, source_filter, q, target, target_filter, l);

  int left = q - 1;

  if(left < 0)
    left = 0;

  float d = dist(nnf[left], source, source_filter, q, target, target_filter, l);

  if(d < dist_min){
    dist_min = d;
    pixel_min = left;
  }

  int top = q - target[l].cols;

  if(top < 0)
    top = 0;

  d = dist(nnf[top], source, source_filter, q, target, target_filter, l);

  if(d < dist_min){
    dist_min = d;
    pixel_min = top;
  }

  return pixel_min;
}


int
patchMatch( Pyramid* source,
	    Pyramid* source_filter,
	    Pyramid* target,
	    Pyramid* target_filter,
	    int l,
	    int q)
{
  int pixelMatch = 0, tmp = 0;
  for(int k=0; k<NB_MAX_ITER; k++){
    int* nnf = initialisation(source[l].cols, source[l].rows, target[l].cols, target[l].rows);

    tmp = iteration(source, source_filter, target, target_filter, l, q, nnf);

    if(k > 0){
      if( dist(pixelMatch, source, source_filter, q, target, target_filter, l) > dist(tmp, source, source_filter, q, target, target_filter, l) ){
	pixelMatch = tmp;
      }
    }
    else{
      pixelMatch = tmp;
    }

    free(nnf);
  }

  return pixelMatch;
}

