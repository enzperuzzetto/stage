#include <search_structures.h>

//#include <math.h>


int
BestMatch(float* ex, float* ex_prim, float* input, float* output, int* s, int l, int q)
{
  (void)ex;
  (void)ex_prim;
  (void)input;
  (void)output;
  (void)l;
  (void)s;
  (void)q;
  /*int p_app = BestApproximateMatch(ex, ex_prim, input, output, l, q);
  int p_coh = BestCoherenceMatch(ex, ex_prim, input, output, s, l, q);
  float dapp = 0.0;
  float dcoh = 0.0;

    if(dcoh <= dapp*(1+powf(2, l-L)*k) )
    return p_coh;
  else
    return p_app;
  */
  return -1;
}

int
BestApproximateMatch(float* ex, float* ex_prim, float* input, float* output, int l, int q)
{
  (void)ex;
  (void)ex_prim;
  (void)input;
  (void)output;
  (void)l;
  (void)q;

  return -1;
}

int
BestCoherenceMatch(float* ex, float* ex_prim, float* input, float* output, int* s, int l, int q)
{
  (void)ex;
  (void)ex_prim;
  (void)input;
  (void)output;
  (void)l;
  

  int r = 0;
  return s[r] + (q-r);
}
