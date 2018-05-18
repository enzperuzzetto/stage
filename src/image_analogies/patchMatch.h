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
iteration();

/**
 * @brief
 *
 * @Param
 **/
int
search();

/**
 * @brief
 *
 * @Param
 **/
int
propagation();

/**
 * @brief
 *
 * @Param
 **/
int
patchMatch( int cols_s,
	    int cols_t,
	    Pyramid* source,
	    Pyramid* source_filter,
	    int cols_t,
	    int rows_t,
	    Pyramid* target,
	    Pyramid* target_filter,
	    int q);
