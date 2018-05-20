//----------------------------------------------------------------------
//      File:           main.cpp
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Fichier pricipal pour l'analogie d'image
//----------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "loadSaveTiffClass.h"
#include "Point.h"
#include "BestCoherence.h"
#include "tiffio.h"
#include "annUse.h"
#include <vector>
#include <math.h>

static float K=1;
int gk=2;
int pk=1;
//nombre de niveaux de la pyramide gaussienne
int level=1;
bool DEBUG = false;

//true si on sauvegarde le résultat en mosaique avec les images sources
bool mosaique=false;
//true si on remappe les luminances après calcul
bool remapping=false;
//true si on effectue une texture-by-numbers
bool bynumbers=false;

using namespace std;

//Retourne une image crée par analogie avec les imagesA,A' et B.
//L'algorithme utilisera une pyramide gaussienne a level niveaux
Image *CreateImageAnalogy(Image *A,Image *Ap,Image * B, FeatureType Type)
{
  Image *Bp = new Image(B->getHauteur(0),B->getLargeur(0));

  //calcul des pyramides gaussiennes  
  if (DEBUG) {
    cout << "Calcul des pyramides gaussiennes..." << endl;
  }
  A->calculerLesPyramides(level);
  Ap->calculerLesPyramides(level);
  B->calculerLesPyramides(level);
  Bp->calculerLesPyramides(level);
  if (DEBUG) {
    cout << "Calcul terminé!" << endl;
  }
  
  //initialises des structures
  if (DEBUG) {
    cout << "Initialisation: niveau " << level << "  (taille de l'image: " << B->getLargeur(level) << "x" << B->getHauteur(level) << ")" << endl;
  }
  
  Point * q;
  Point * p;
  
  //calcul de l'image B' au niveau level en utilisant la force brute
  float progression = 0;
  for(uint32 y=0 ; y<B->getHauteur(level) ; y++) {
    for(uint32 x=0 ; x<B->getLargeur(level) ; x++) {
      if (DEBUG) {
	progression = ((float) (y*B->getLargeur(level) + x)) / ((float) (B->getHauteur(level)*B->getLargeur(level)));
	cout << "\r  Progression: " << 100*progression << "%    ";
      }
      q = new Point(x,y);
      p = A->ComputeBruteForce(Ap, B, Bp, q, level, gk, Type);
      Bp->setFeature(x,y,level,Ap->getFeature(p->getX(),p->getY(),level));
      delete p;
      delete q;
    }
  }

  if (DEBUG) {
    cout << "\r  Progression: 100%        " << endl ;
  }

  //calcul des images de B' aux niveau level-1 à 0
  Point **s;
  
  for(int l=level-1;l>=0;l--)
    {
      if (DEBUG) {
	cout << " * level " << l << " (taille: " << B->getLargeur(l) << "x" << B->getHauteur(l) << ")" << endl;
      }
      //initialisation de s
      s = new Point* [(B->getLargeur(l))*(B->getHauteur(l))];
      //initialisation de la structure de recherche pour ANN
      annUse *ann=new annUse(A,Ap,l,pk,gk,Type);
      float progression = 0;
      //calcul de B' au niveau l
      for(uint32 y=0 ; y<B->getHauteur(l) ; y++)
	for(uint32 x=0 ; x<B->getLargeur(l) ; x++) {
	  if (DEBUG) {
	    progression = ((float) (y*B->getLargeur(l) + x)) / ((float) (B->getHauteur(l)*B->getLargeur(l)));
	    cout << "\r     Progression: " << 100*progression << "%    ";
	  }
	  q = new Point(x,y);
	  p = BestMatch(A,Ap,B,Bp,l,q,s,ann,pk,gk,level,K,Type);
	  Bp->setFeature(x,y,l,Ap->getFeature(p->getX(),p->getY(),l));
	  s[indS(l,x,y,B)]=new Point(*p);
	  delete p;
	  delete q;
	}
      if (DEBUG) {
	cout << "\r     Progression: 100%       " <<endl;
      }
      
      delete ann;
      for (uint32 tmp=0 ; tmp<(B->getLargeur(l)*B->getHauteur(l)) ; tmp++) {
	delete s[tmp];
      }
      delete [] s;
    }
  
  if(remapping)
    Bp->luminanceRemapping(Ap);

  //calcul des couleurs de B' en fonction de B
  if(!bynumbers)
    Bp->computeRGB(B);

  return Bp;
}

//Sauvegarde l'image de nom res comme mosaique des images A, Ap, B et Bp
//i et j représente l'eppaisseur des bandes noires qui séparent les images
void Mosaique(Image *A, Image *Ap, Image *B, Image *Bp,uint32 i,uint32 j,char* res){

  int a1=A->getLargeur(0);
  int b1=B->getLargeur(0);
  int max;
  if(a1>b1){
    max=2*a1;
  }
  else {
    max=2*b1;
  }
  Image *resultat = new Image(A->getHauteur(0)+B->getHauteur(0)+j,max+i);

  //recopie de A  

  for(uint32 y=0;y<A->getHauteur(0);y++)
    {
      for(uint32 x=0;x<A->getLargeur(0);x++)
	{
	  resultat->setFeature(x,y,0,A->getFeature(x,y,0));
	}
    }

 //recopie de Ap  

  for(uint32 y=0;y<Ap->getHauteur(0);y++)
    {
      for(uint32 x=0;x<Ap->getLargeur(0);x++)
	{
	  resultat->setFeature(max/2+x+i,y,0,Ap->getFeature(x,y,0));
	}
    }

 //recopie de B

  for(uint32 y=0;y<B->getHauteur(0);y++)
    {
      for(uint32 x=0;x<B->getLargeur(0);x++)
	{
	  resultat->setFeature(x,y+j+A->getHauteur(0),0,B->getFeature(x,y,0));
	}
    }

  //recopie de Bp

  for(uint32 y=0;y<Bp->getHauteur(0);y++)
    {
      for(uint32 x=0;x<Bp->getLargeur(0);x++)
	{
	  resultat->setFeature(max/2+x+i,j+y+A->getHauteur(0),0,Bp->getFeature(x,y,0));
	}
    }

  resultat->sauvegarderImageCouleur(res,0);
  delete resultat;
}

void usage() 
{
  cout << endl;
  cout << "ANALOGIMAGE : Usage" << endl;
  cout << "     AnalogImage <A> <A'> <B> <B'> [options]" << endl;
  cout << "     Options:" << endl;
  cout << "      -h/-help : Affiche l'aide" << endl;
  cout << "      -d/-debug : Affiche des messages de débuggage" << endl;
  cout << "      -t/-type <LUM/RGB/RGBLUM/IQ/DEPLUM/DEPRGB> : informations utilisées pour le calcul (par défaut: LUM)" << endl;
  cout << "      -Ad <Ad> : Carte de profondeur de A" << endl;
  cout << "      -Bd <Bd> : Carte de profondeur de B" << endl;
  cout << "         (uniquement si type = DEPLUM ou DEPRGB)" << endl;
  cout << "      -k <K> : Valeur de K (par défaut: 1.0)" << endl;
  cout << "      -l/-level <L> : Nombre de niveaux des pyramides gaussiennes (par défaut: 1)" << endl;
  cout << "      -pk <pk> : Taille du voisinage au niveau l-1 (par défaut: 1)" << endl;
  cout << "      -gk <gk> : Taille du voisinage au niveau l (par défaut: 2)" << endl;
  cout << "      -r/-remapping : Applique un remapping des luminances avant l'algo" << endl;
  cout << "      -m/-mosaic : L'image de sortie contient A, A', B, et B'" << endl;
  cout << "      -c/-color <I> : Récupère les couleurs de l'image I = Ap ou B (par défaut, de l'image B)" << endl;
  cout << endl;
}


void help()
{
  usage();
  cout << endl;
  cout << "TYPES DES PARAMETRES:" << endl;
  cout << "    les A, A', Ad, B, B', et Bd sont des noms d'images .tif" << endl;
  cout << "    K est un flottant" << endl;
  cout << "    pk, gk, et l, sont des entiers" << endl;
  cout << "    I est soit la chaine \"Ap\" soit la chaine \"B\"" << endl;
  cout << "UTILITE DES PARAMETRES:" << endl;
  cout << "    K détermine si on utilise le coherence match ou l'approximate match." << endl;
  cout << "       Plus K est grand, plus les textures de l'image B' sont cohérentes avec celles de l'image A'" << endl;
  cout << "    pk détermine la taille du voisinage considéré dans les images de niveau l-1." << endl;
  cout << "       On considère le pixel p +/-pk pixels autour." << endl;
  cout << "       pk=1: voisinage de 3x3 pixels (valeur par défaut)" << endl;
  cout << "    gk détermine la taille du voisinage considéfé dans les images de niveau l." << endl;
  cout << "       gk=2: voisinage de 5x5 pixels (valeur par défaut)" << endl;
  cout << "    l (level) détermine le nombre de niveaux des pyramides gaussiennes calculées." << endl;
  cout << "       chaque niveau correspond à l'image du niveau précédent réduite de moitié (en utilisant un flou gaussien)." << endl;
  cout << "       si à l=0, taille=400x400, alors à l=1, taille=200x200, puis à l=2, taille=100x100, etc..." << endl;
  cout << "    le remapping permet d'aligner les luminances des images." << endl;
  cout << "       Il faut l'utiliser si l'image B est plus claire ou plus foncée que l'image A." << endl;
  cout << "    la mosaic permet d'obtenir une image regroupant les trois images en entrée de l'algo et celle obtenue en sortie." << endl;
  cout << "       Par défaut, la sortie de l'algo est la seule image B'." << endl;
  cout << "    color permet de récupérer la couleur soit de l'image A' (cas de la synthèse de textures)," << endl;
  cout << "       soit de l'image B (cas des filtres artistiques)" << endl;
  cout << "    Ad et Bd sont les cartes de profondeur," << endl;
  cout << "       on ne les utilise que si le type est soit DEPLUM, soit DEPRGB." << endl;
  cout << "QUELQUES CONSEILS:" << endl;
  cout << "    Choisissez bien vos paramètres en fonction de ce que vous voulez faire:" << endl;
  cout << "       Transfert de textures                   : L>1 , 0.01<K<0.50 , pk>0 , gk>1 , -color Ap" << endl;
  cout << "       Filtres artistiques                     :       0.50<K<2.00 , pk>0 , gk>1 , -color B" << endl;
  cout << "       Brouillard ou flou (avec la profondeur) : L>0 ,      K<0.10 , pk<1 , gk<2 , -color B" << endl;
  cout << endl;
}

//////////////////////////////////////////////////
// Fonction principale
//////////////////////////////////////////////////
int main(int argc, char** argv)
{  
  if (argc < 5) {
    for(int i=0; i<argc ; i++) {
      if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"-help")) {
	help();
	exit(0);
      }
    }
    usage();
    exit(-1);
  }
  if (!strstr (argv[1], ".tif")) {
    cout << "Argument n°" << 1 << ": " << argv[1] << "; Erreur: Image .tif valide attendue" << endl;
    exit(-1);
  }
  if (!strstr (argv[2], ".tif")) {
    cout << "Argument n°" << 2 << ": " << argv[2] << "; Erreur: Image .tif valide attendue" << endl;
    exit(-1);
  }
  if (!strstr (argv[3], ".tif")) {
    cout << "Argument n°" << 3 << ": " << argv[3] << "; Erreur: Image .tif valide attendue" << endl;
    exit(-1);
  }
  if (!strstr (argv[4], ".tif")) {
    cout << "Argument n°" << 4 << ": " << argv[4] << "; Erreur: Image .tif valide attendue" << endl;
    exit(-1);
  }

  char* nomA = argv[1];
  char* nomAp = argv[2];
  char* nomB = argv[3];
  char* nomBp = argv[4];
  char* nomAd = NULL;
  char* nomBd = NULL;

  FeatureType Type=LUMIN;
  //récupération des arguments de la commande
  for(int i=5; i<argc ; i++)
    {
      if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"-help")) {
	help();
	exit(0);
      } else if(!strcmp(argv[i],"-d") || !strcmp(argv[i],"-debug")) {
	DEBUG = true;
      } else if(!strcmp(argv[i],"-t") || !strcmp(argv[i],"-type")) {
	i++; if (i==argc) {exit(-1);}
	if (!strcmp(argv[i],"LUM")) {
	  Type = LUMIN;
	} else if (!strcmp(argv[i],"RGB")) {
	  Type = RGB;
	} else if (!strcmp(argv[i],"RGBLUM")) {
	  Type = RGBLUM;
	} else if (!strcmp(argv[i],"IQ")) {
	  Type = IQ;
	} else if (!strcmp(argv[i],"DEPRGB")) {
	  Type = DEPRGB;
	} else if (!strcmp(argv[i],"DEPLUM")) {
	  Type = DEPLUM;
	} else {
	  cout << "Argument n°" << i << ": " << argv[i] << "; Attention: Après -t ou -type, il faut mettre \"LUM\", \"RGB\", \"RGBLUM\", \"IQ\", \"DEPLUM\", ou \"DEPRGB\" seulement" << endl;
	}
      } else if(!strcmp(argv[i],"-m") || !strcmp(argv[i],"-mosaic")) {
	mosaique=true;
      } else if(!strcmp(argv[i],"-r") || !strcmp(argv[i],"-remapping")) {
	remapping=true;
      } else if(!strcmp(argv[i],"-Ad") || !strcmp(argv[i],"-AD")) {
	i++; if (i==argc) {exit(-1);}
	nomAd = argv[i];
      } else if(!strcmp(argv[i],"-Bd") || !strcmp(argv[i],"-BD")) {
	i++; if (i==argc) {exit(-1);}
	nomBd = argv[i];
      } else if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-K")) {
	i++; if (i==argc) {exit(-1);}
	K = atof(argv[i]);
      } else if(!strcmp(argv[i],"-pk")) {
	i++; if (i==argc) {exit(-1);}
	pk= atoi(argv[i]);
      } else if(!strcmp(argv[i],"-gk")) {
	i++; if (i==argc) {exit(-1);}
	gk= atoi(argv[i]);
      } else if(!strcmp(argv[i],"-l") || !strcmp(argv[i],"-L") || !strcmp(argv[i],"-level")) {
	i++; if (i==argc) {exit(-1);}
	level = atoi(argv[i]);
      } else if(!strcmp(argv[i],"-c") || !strcmp(argv[i],"-color")) {
	i++; if (i==argc) {exit(-1);}
	if (!strcmp(argv[i],"Ap")) {
	  bynumbers=true;
	} else if (!strcmp(argv[i],"B")) {
	  bynumbers=false;
	} else {
	  cout << "Argument n°" << i << ": " << argv[i] << "; Attention: Après -c ou -color, il faut mettre \"A'\" ou \"B\" seulement" << endl;
	}
      } else { 
	cout << "Argument n°" << i << ": " << argv[i] << " : Argument incorrect, ignoré" << endl;
      }
    }

  if (DEBUG) {
    cout << " -- Lancement de l'algorithme --" << endl;
  }

  // Lecture des images
  vector<uint8> f1;
  Image *A = new Image(nomA, nomAd);  
  Image *Ap= new Image(nomAp);
  Image *B = new Image(nomB, nomBd);
  
  if (DEBUG) {
    cout << "Images A, A', et B chargées" << endl;
    cout << "Création de B'..." << endl;
  }

  //création de l'image par analogie
  Image *Bp = CreateImageAnalogy(A, Ap, B, Type);
  if (DEBUG) {
    cout << "Images B' terminée" << endl;
  }

  //sauvegarde de l'image
  if(!mosaique)
    Bp->sauvegarderImageCouleur(nomBp, 0);
  else
    Mosaique(A,Ap,B,Bp,3,3,nomBp);

  delete A;
  delete Ap;
  delete B;
  delete Bp;

  if (DEBUG) {
    cout << " -- Fin du calcul --" << endl;
  }

  return 0;
}
