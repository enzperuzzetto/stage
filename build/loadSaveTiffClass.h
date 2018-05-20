#ifndef __LOADSAVETIFF__
#define __LOADSAVETIFF__

#include "classeFeatures.h"

#include <stdio.h>
#include "tiffio.h"

class Features;
class Point;

class Image {
  
 public:
  Image(int hauteur, int largeur);
  Image(TIFF* fichier, TIFF* fichierProf=NULL);
  Image(char* nomFichier, char* nomFichierProf=NULL);
  Image(Image* source, int ratio=1);
  ~Image();

  Features* getFeature(int x, int y,int l);
  uint32 getHauteur (int l);
  uint32 getLargeur (int l);
  int calculerLesPyramides(int nombreDeNiveaux);
  Point* ComputeBruteForce(Image* imageAp, Image* imageB, Image* imageBp, Point* q, int level, int gk, FeatureType Type);
  
  int sauvegarderImageCouleur(TIFF* fichier,int l);
  int sauvegarderImageCouleur(char* nomFichier,int l);

  void setFeature(int x,int y,int l,Features* f);

  void luminanceRemapping(Image* B);

  int computeRGB(Image* colorImage);

  int moyenneLum();
  int varianceLum(int moy);
 private:
  Features** img;
  uint32 hauteur;
  uint32 largeur;
  // pyramides gaussiennes
  Image* precedente;
  Image* suivante;

  int LoadPixels (TIFF* file, uint32** memory, uint32* width, uint32* height);
  int SwapBuffer (uint32** memoryD, uint32* memoryS, uint32 width, uint32 height);
  int SavePixels (TIFF* file, Features** img, uint32 width, uint32 height);
  
};

#endif
