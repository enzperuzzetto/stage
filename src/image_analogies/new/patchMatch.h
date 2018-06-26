#ifndef __PATCH_MATCH_H__

#define __PATCH_MATCH_H__

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "pyramid.h"
#include "BestMatch.h"
#include "stats.h"

#define ALPHA 0.5


/**
 * @brief Phase d'initialisation
 *        Créer un tableau aléatoire des indices de A de taille de B
 *
 * @Param cols_s, rows_s: dimension de A
 *        cols_t, rows_t: dimension de B
 **/
//int*
//initialisation(int cols_s, int rows_s, int cols_t, int rows_t);

/**
 * @brief Phase d'itération
 *        propagation puis recherche aléatoire
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramide
 *        q: pixel de B et B'
 *        nnf: tableau aléatoire de la phase d'initialisation
 **/
//int
//iteration(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief Phase de recherche aléatoire
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramide
 *        q: pixel de B et B'
 *        nnf: tableau aléatoire de la phase d'initialisation
 **/
//int
//search(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief Phase de propagation
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramide
 *        q: pixel de B et B'
 *        nnf: tableau aléatoire de la phase d'initialisation
 **/
//int
//propagation(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief algorithme patchMatch regroupant toutes les fonctions décrites au-dessus
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramide
 *        q: pixel de B et B'
 **/
int
patchMatch( Pyramid* A, Pyramid* Aprim, int x, int y, Pyramid* B, Pyramid* Bprim, int l, int size, float W, float levelWeight, float** kernel, int* nnf, int k);

#endif
