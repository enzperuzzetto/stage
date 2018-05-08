#ifndef __SEARCH_STRUCTURES__
#define __SEARCH_STRUCURES__

/**
 * @file search.h
 * @brief 
 **/

int
BestMatch(float* ex, float* ex_prim, float* input, float* output, int* s, int l, int q);

int
BestApproximateMatch(float* ex, float* ex_prim, float* input, float* output, int l, int q);

int
BestCoherenceMatch(float* ex, float* ex_prim, float* input, float* output, int* s, int l, int q);

#endif
