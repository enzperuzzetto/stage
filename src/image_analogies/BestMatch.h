#ifndef __BESTMATCH_H__

#define __BESTMATCH_H__

#include <patchMatch.h>

/**
 * @brief
 *
 * @Param
 **/
int
bruteForceMatch(Pyramid* source,
		Pyramid* source_filter,
		Pyramid* target,
		Pyramid* target_filter,
		int l,
		int q);

/**
 * @brief
 *
 * @Param
 **/
int
bestCoherenceMatch(Pyramid* source,
		   Pyramid* source_filter,
		   Pyramid* target,
		   Pyramid* target_filter,
		   int l,
		   int q);

/**
 * @brief
 *
 * @Param
 **/
int
bestMatch(Pyramid* source,
	  Pyramid* source_filter,
	  Pyramid* target,
	  Pyramid* target_filter,
	  int l,
	  int q);

#endif
