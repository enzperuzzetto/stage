#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int*
initialisation(int cols_s, int rows_s, int cols_t, int rows_t)
{
  int* nnf = malloc(sizeof(int) * cols_s * rows_s);
  srand(time(NULL));
  for(int k=0; k<cols_s*rows_s; k++){
    nnf[k] = rand()%(cols_t*rows_t);
  }

  return nnf;  
}

int
main(void)
{

  int* nnf = initialisation(5, 5, 10, 10);

  for(int i=0; i<25; i++){
    printf("%d ",nnf[i]);
  }

  return 1;
}
