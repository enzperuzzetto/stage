#include "patchMatch.h"

/**
 * @brief Permet d'initailiser un tableau d'indices random.
 *
 * @fn int* initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
 * @param cols_s,row_s dimensions du tableau de sortie.
 * @param cos_t,rows_t dimensions de l'image où prendre les indices.
 * @return un tableau d'indices random d'une image.
 **/
int*
initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
{
  int* nnf = malloc(sizeof(int) * cols_t * rows_t);
  srand(time(NULL));
  for(int i=0; i<cols_t*rows_t; i++){
    nnf[i] =  rand()%(cols_s*rows_s);
  }

  return nnf;
}

/**
 * @brief Etape de la recherche de l'algorithme de pathcMatch.
 *        Recherche aléatoire dans une fenêtre de recherche qui diminue à chaque itération de boucle.
 *
 * @fn int search(Pyramid* A, Pyramid* Aprim, int q, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel, int* nnf)
 * @param A,A' pyramide gaussienne de l'image source A et A'.
 * @param q indice du pixel des images B et B'.
 * @param B,B' pyramide gaussienne des images B et B'.
 * @parma l niveau de la pyramide gaussienne.
 * @param size taille du patch/voisinnage.
 * @param float W: poid soulignant la similarité entre (A,B) et (A',B').
 * @param int weighting: booleen qui détermine l'utilisation de W.
 * @param int skipCenter: booleen qui détermine si l'on fait sur le L-shape du patch ou sur le patch entier.
 * @param onePixel: booleen qui prend que le pixel centrale au niveau l.
 * @param float** kernel: matrice de noyaux gaussiens.
 * @param nnf  tableau d'indice de A et A' aléatoire.
 * @return le pixel ayant la plus petiite distance à la fin des itérations.
 **/
int
search(Pyramid* A, Pyramid* Aprim, int q, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, int onePixel, float** kernel, int* nnf)
{
  int w = max(A[l].cols, A[l].rows), pixel_x = 0, pixel_y = 0, pixel = 0, pixel_min = nnf[q], i = 0;
  int xt = q/B[l].cols;
  int yt = q - xt*B[l].cols;
  int xs = pixel_min/A[l].cols;
  int ys = pixel_min - xs*A[l].cols;

  float dist_min = dist(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, onePixel, kernel), d = 0.0, alpha = powf(ALPHA, i), x = 0.0, y = 0.0;

  while( w*alpha > 1){
    x = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;
    y = ((float)rand()/ (float)(RAND_MAX/3.0)) -1;

    pixel_x = w*alpha*x;
    pixel_y = w*alpha*y;
    pixel = pixel_min + (pixel_x*A[l].cols + pixel_y);

    if(pixel < 0)
      pixel = 0;
    else if(pixel > A[l].cols*A[l].rows-1)
      pixel = A[l].cols*A[l].rows-1;

    xs = pixel/A[l].cols;
    ys = pixel - xs*A[l].cols;
    d = dist(xs, ys, A, Aprim, xt, yt, B, Bprim, l, size, W, levelWeight, onePixel, kernel);

    if( d < dist_min){
      dist_min = d;
      pixel_min = pixel;
    }

    i++;
    alpha = powf(ALPHA, i);
  }

  return pixel_min;
}

/**
 * @brief Etapde de propagation de l'algorithme de patchMatch.
 *        Compare et retourne le pixel qui à la plus petite distance entre le pixel q, son voisin de gauche et son voisin du dessus.
 *
 * @fn int propagation(Pyramid* A, Pyramid* Aprim, int x, int y, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel, int* nnf)
 * @param A,A' pyramide gaussienne de l'image source A et A'.
 * @param q indice du pixel des images B et B'.
 * @param B,B' pyramide gaussienne des images B et B'.
 * @parma l niveau de la pyramide gaussienne.
 * @param size taille du patch/voisinnage.
 * @param float W: poid soulignant la similarité entre (A,B) et (A',B').
 * @param int weighting: booleen qui détermine l'utilisation de W.
 * @param int skipCenter: booleen qui détermine si l'on fait sur le L-shape du patch ou sur le patch entier.
 * @param onePixel: booleen qui prend que le pixel centrale au niveau l.
 * @param float** kernel: matrice de noyaux gaussiens.
 * @param nnf  tableau d'indice de A et A' aléatoire.
 * @return Soit le pixel q soit son pixel de gauche soit celui du dessus.
 **/
int
propagation(Pyramid* A, Pyramid* Aprim, int x, int y, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, int onePixel, float** kernel, int* nnf)
{
  int q = x*B[l].cols+y;
  int pixel_min = q;
  int p = nnf[pixel_min];
  int xs = p/A[l].cols;
  int ys = p - xs*A[l].cols;
  float dist_min = dist(xs, ys, A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel);

  int left = q - 1;
  if(left < 0)
    left = 0;

  p = nnf[left];
  xs = p/A[l].cols;
  ys = p - xs*A[l].cols;

  float d = dist(xs, ys, A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel);

  if(d < dist_min){
    dist_min = d;
    pixel_min = left;
  }

  int top = q - B[l].cols;
  if(top < 0)
    top = 0;

  p = nnf[top];
  xs = p/A[l].cols;
  ys = p - xs*A[l].cols;

  d = dist(xs, ys, A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel);

  if(d < dist_min){
    dist_min = d;
    pixel_min = top;
  }

  return pixel_min;
}

/**
 * @brief Appel des fonctions propagation puis search
 **/
int
iteration(Pyramid* A, Pyramid* Aprim, int x, int y, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, int onePixel, float** kernel, int* nnf)
{
  int pixelTarget = propagation(A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel, nnf);
  int pixelSource = search(A, Aprim, pixelTarget, B, Bprim, l, size, W, levelWeight, onePixel, kernel, nnf);

  return pixelSource;
}

/**
 * @brief Une itération de l'algorithme patchMatch
 **/
int
patchMatch( Pyramid* A, Pyramid* Aprim, int x, int y, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, int onePixel, float** kernel, int* nnf, int k)
{
  int pixelMatch = nnf[x*B[l].cols+y], tmp = 0;
  int xs = pixelMatch/A[l].cols;
  int ys = pixelMatch - xs*A[l].cols;

  tmp = iteration(A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel, nnf);
  int xxs = tmp/A[l].cols;
  int yys = tmp - xxs*A[l].cols;

  if(k > 0){
    if( dist(xs, ys, A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel) > dist(xxs, yys, A, Aprim, x, y, B, Bprim, l, size, W, levelWeight, onePixel, kernel) ){
       pixelMatch = tmp;
    }
  }
  else{
    pixelMatch = tmp;
  }

  return pixelMatch;
}
