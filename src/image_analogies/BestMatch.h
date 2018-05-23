#ifndef __BESTMATCH_H__

#define __BESTMATCH_H__

#include "patchMatch.h"

#define NEIGHBOOR_SIZE_COARSER 3
#define NEIGHBOOR_SIZE_FINER 5
#define K 0.5
#define L 3


/**
 * @brief Recherche correspondance du pixel q par brute force
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramid
 *        q: pixel de B et B'
 **/
int
bruteForceMatch(Pyramid* source,
		Pyramid* source_filter,
		Pyramid* target,
		Pyramid* target_filter,
		int l,
		int q);

/**
 * @brief Recherche de la meilleur coherence du pixel q
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramid
 *        q: pixel de B et B'
 **/
int
bestCoherenceMatch(Pyramid* source,
		   Pyramid* source_filter,
		   Pyramid* target,
		   Pyramid* target_filter,
		   int l,
		   int q);

/**
 * @brief Permet de trouver la meilleur correspondance du pixel q
 *        PatchMatch/Coherence soit BruteForce/Coherence
 *
 * @Param source, source_filter, target, target_filter: A A' B B'
 *        l: niveau de la pyramid
 *        q: pixel de B et B'
 **/
int
bestMatch(Pyramid* source,
	  Pyramid* source_filter,
	  Pyramid* target,
	  Pyramid* target_filter,
	  int l,
	  int q);

#endif
