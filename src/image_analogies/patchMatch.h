#ifndef __PATCH_MATCH_H__

#define __PATCH_MATCH_H__

#include <time.h>

#include <pyramid.h>


/**
 * @brief
 *
 * @Param
 **/
int
initialisation(int cols_s, int rows_s, int cols_t, int rows_t);

/**
 * @brief
 *
 * @Param
 **/
int
iteration(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief
 *
 * @Param
 **/
int
search(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief
 *
 * @Param
 **/
int
propagation(Pyramid* source, Pyramid* source_filter, Pyramid* target, Pyramid* target_filter, int l, int q, int* nnf);

/**
 * @brief
 *
 * @Param
 **/
int
patchMatch( Pyramid* source,
	    Pyramid* source_filter,
	    Pyramid* target,
	    Pyramid* target_filter,
	    int l,
	    int q);


#endif
