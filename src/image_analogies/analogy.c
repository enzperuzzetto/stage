#include <stdio.h>
#include <stdlib.h>

#include "pyramid.h"
#include "features.h"
#include "search_structures.h"


void
process(char *ims_name, char *ims_filter_name, char *imd_name)
{
  (void)ims_name;
  (void)ims_filter_name;
  (void)imd_name;
}
  



void
usage (char *s){
  fprintf(stderr, "Usage: %s <A> <A'> <B>  \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 3
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1], argv[2], argv[3]);
  return EXIT_SUCCESS;
}
