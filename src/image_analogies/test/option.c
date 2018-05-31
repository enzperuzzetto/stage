#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>


int
main(int argc, char** argv)
{
  char c;
  char* optstring = "s:t:f:h";
  while((c = getopt(argc, argv, optstring)) !=EOF ){
    switch (c){
    case 'f':
      printf("-filtre");
      break;
    case 's':
      printf("-synthèse");
      break;
    case 't':
      printf( "-transfert");
      break;
    case 'h':
      fprintf(stderr, "Usage: %s [option] images\n"
	      "Options:\n"
	      "-filtre        applique le filtre de A A' sur B\n"
	      "-synthése      crée une texture à partir de A'au dimension de l'image B\n"
	      "-transfert     transfert de style de l'image A sur l'image B\n"
	      "-h/--help      message d'aide\n"
	      , argv[0]);
      return 0;
    }
  }

  if(optind == 1 && argc != 4){
     fprintf(stderr, "Usage: %s [option] images\n"
	      "Options:\n"
	      "-filtre        applique le filtre de A A' sur B\n"
	      "-synthése      crée une texture à partir de A'au dimension de l'image B\n"
	      "-transfert     transfert de style de l'image A sur l'image B\n"
	      "-h/--help      message d'aide\n"
	      , argv[0]);
      return 0;
  }
    
  while(optind < argc){
    printf(" %s ", argv[optind]);
    optind++;
  }


  return 1;

}
