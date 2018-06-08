#ifndef __PYRAMID_H__
#define __PYRAMID_H__


#define K 0.5
#define L 3


/**
 * @brief Retourne dans un tableau les dimensions de l'image au niveau level de la pyramide gaussienne
 *
 * @Param cols: dimension en largeur de l'image d'origine
 *        rows: dimension en hauteur de l'image d'origine
 *   
 **/
int*
pyramid_level_dim(int cols, int rows, int level);


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
