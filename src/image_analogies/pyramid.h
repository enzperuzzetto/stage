#ifndef __PYRAMID_H__

#define __PYRAMID_H__

#define NEIGHBOOR_SIZE_COARSER 3
#define NEIGHBOOR_SIZE_FINER 5
#define K 5
#define L 3
#define W 1.0

typedef struct Pyramid
{
  int cols;
  int rows;
  float* data;
  float* mean;
  float* sd;
  int* s;
}Pyramid;

/**
 * @brief Retourne dans un tableau les dimensions de l'image au niveau level de la pyramide gaussienne
 *
 * @Param cols: dimension en largeur de l'image d'origine
 *        rows: dimension en hauteur de l'image d'origine
 *        level: niveau de la pyramid gaussienne
 **/
int*
pyramid_level_dim(int cols, int rows, int level);

/**
 * @brief Retourne un tableau de taille L représentent chaque niveau de la pyramid. Permet juste l'allocation en mémoire de la pyramid
 *
 * @Param cols: dimension en largeur de l'image d'origine
 *        rows: dimension en hauteur de l'image d'origine
 *        s   : 0 si on doit allouer le tableau s 1 si on doit l'allouer
 **/
Pyramid*
init_pyramid(int cols, int rows, int s);

/**
 * @brief Libère la mémoire de la pyramid
 *
 * @Param pyramid: le tableau a libérer
 *        s      : 1 si il possède le tableau s, 0 sinon
 **/
void
free_pyramid(Pyramid* pyramid, int s);

/**
 * @brief Calcul l'image au niveau level de la pyramid
 *
 * @Param cols , rows : dimension de l'image de base
 *        data: donnée de l'image de base (luminance)
 *        level: niveau de la pyramide
 **/
float*
pyramid_level(int cols, int rows, float* data, int level);

/**
 * @brief Calcul la distance euclidienne de A B A' B' aux pixels p et q et au niveau l et l-1
 *
 * @Param p: pixel de A et A'
 *        source, source_filter: données de A et A'
 *        q: pixel de B et B'
 *        target, target_filter: données de B et B'
 *        l: niveau de la pyramide
 **/
float
dist( int p, Pyramid* source, Pyramid* source_filter, int q, Pyramid* target, Pyramid* target_filter, int l);

#endif
