#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "transfert.h"
#include "synthese.h"
#include "filtre.h"


void
Usage(char* s)
{
  fprintf(stderr, "Usage: %s traitement fichier option\n",s);
  fprintf(stderr, "traitement: filtre / synthese / transfert\n");
  fprintf(stderr, "\nPour filtre:\nA A' B et size(taille voisinage)\n");
  fprintf(stderr, "\nPour synthese:\nA' et dimension (w, h), size(taille voisinage)\n");
  fprintf(stderr, "\nPour transfert:\nA A B et w ]0,1[, size(taille voisinage)\n");
  exit(EXIT_FAILURE);
}

int
main(int argc, char** argv)
{
  if(argc != 6) Usage(argv[0]);

  int size = atoi(argv[5]);
  if(size < 0) Usage(argv[0]);
  
  if(strcmp(argv[1], "filtre") == 0)
    filtre(argv[2], argv[3], argv[4], size);

  else if(strcmp(argv[1], "synthese") == 0){
    int width = atoi(argv[3]);
    int height = atoi(argv[4]);
    synthese(argv[2], width, height, size);
  }
  else if(strcmp(argv[1], "transfert") == 0){
    float w = atof(argv[4]);
    //if(w<=0) //|| w>1)
	//Usage(argv[0]);
    transfert(argv[2], argv[3], w, size);
  }
  else
    Usage(argv[0]);

  return EXIT_SUCCESS;
}
    
