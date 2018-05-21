#ifndef __FEATURES_H__

#define __FEATURES_H__



/**
 * @brief Calcule la moyenne du voisinnage du pixel
 *
 * @Param cols, rows: dimension de l'image
 *        in: tableau de données
 *        index: pixel a calculer
 *        size: taille du voisinnage
 *        half: 1 si l'on calcule que la moitier (A' B'), 0 sinon
 **/
float
mean(int cols, int rows, float* in, int index, int size, int half);


/**
 * @brief Calcule la variance du voisinnage du pixel
 *
 * @Param cols, rows: dimension de l'image
 *        in: tableau de données
 *        mean: moyenne du voisinnage à se pixel
 *        index: pixel a calculer
 *        size: taille du voisinnage
 *        half: 1 si l'on calcule que la moitier (A' B'), 0 sinon
 **/
float
variance(int cols, int rows, float* in, float mean, int index, int size, int half);


/**
 * @brief Calcule l'écart-type du voisinnage du pixel
 *
 * @Param cols, rows: dimension de l'image
 *        in: tableau de données
 *        mean: moyenne du voisinnage à se pixel
 *        index: pixel a calculer
 *        size: taille du voisinnage
 *        half: 1 si l'on calcule que la moitier (A' B'), 0 sinon
 **/
float
standard_deviation(int cols, int rows, float* in, float mean, int index, int size, int half);


/**
 * @brief Calcule la moyenne de chaque pixel en fonction de leur voisinnage
 *
 * @Param cols, rows: dimension de l'image
 *        in: tableau de données
 *        mean: tableau où l'on stocke les moyennes 
 *        size: taille du voisinnage
 *        half: 1 si l'on calcule que la moitier (A' B'), 0 sinon
 **/
void
mean_img(int cols, int rows, float* in, float* means, int size, int half);


/**
 * @brief Calcule l'écart-type de chaque pixel en fonction de leur voisinnage 
 *
 * @Param cols, rows: dimension de l'image
 *        in: tableau de données
 *        mean: tableau des moyennes
 *        sd  : tableau où l'on stocke les écart-types
 *        size: taille du voisinnage
 *        half: 1 si l'on calcule que la moitier (A' B'), 0 sinon
 **/
void
sd_img(int cols, int rows, float* in, float* mean, float* sd, int size, int half);


/**
 * @brief Retourne la valeur max entre deux entiers
 *
 * @Param a, b: deux entiers
 **/
int
max(int a, int b);

#endif
