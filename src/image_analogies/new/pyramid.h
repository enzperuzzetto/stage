#ifndef __PYRAMID_H__

#define __PYRAMID_H__

typedef struct Pyramid
{
  int cols;
  int rows;
  float* lum;
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
init_pyramid(int cols, int rows, int s, int L);

/**
 * @brief Libère la mémoire de la pyramidf
 *
 * @Param pyramid: le tableau a libérer
 *        s      : 1 si il possède le tableau s, 0 sinon
 **/
void
free_pyramid(Pyramid* pyramid, int s, int L);

/**
 * @brief Calcul l'image au niveau level de la pyramid
 *
 * @Param cols , rows : dimension de l'image de base
 *        data: donnée de l'image de base (luminance)
 *        level: niveau de la pyramide
 **/
float*
pyramid_level(int cols, int rows, float* data, int level);

#endif
