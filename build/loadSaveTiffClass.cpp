/*
 * loadSaveTiff.cpp
 * 
 * Le code est le même que le fichier loadSaveTiff du TP3 
 * qui utilisait des images en niveau de gris.
 */

#include "loadSaveTiffClass.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>


using namespace std;
typedef unsigned char uint8;


Image::Image(int hauteur, int largeur) {
  this->hauteur = hauteur;
  this->largeur = largeur;
  img = new Features *[hauteur * largeur];
  for (int i = 0 ; i<hauteur*largeur ; i++) 
    img[i] = new Features();
  precedente = NULL; suivante = NULL;
}

Image::Image(TIFF* fichier,  TIFF* fichierProfondeur) {
  uint32** bufferSource=(uint32**)malloc(sizeof(uint32*));
  LoadPixels (fichier, bufferSource, &largeur, &hauteur);
  uint32** bufferDest=(uint32**)malloc(sizeof(uint32*));
  SwapBuffer (bufferDest, *bufferSource, largeur, hauteur);

  uint32** bufferSourceProf = NULL;
  uint32** bufferDestProf = NULL;
  if (fichierProfondeur != NULL) {
    bufferSourceProf = (uint32**)malloc(sizeof(uint32*));
    LoadPixels (fichierProfondeur, bufferSourceProf, &largeur, &hauteur);
    bufferDestProf = (uint32**)malloc(sizeof(uint32*));
    SwapBuffer (bufferDestProf, *bufferSourceProf, largeur, hauteur);
  }

  img = new Features *[hauteur * largeur];

  for (uint32 i = 0 ; i<hauteur*largeur ; i++) 
    if (fichierProfondeur != NULL) {
      uint8 D =0.3*(((*bufferDestProf)[i] & 0x000000FF))+0.59*(((*bufferDestProf)[i] & 0x0000FF00)>> 8)+0.11*(((*bufferDestProf)[i] & 0x00FF0000)>> 16);
      img[i] = new Features((*bufferDest)[i], D);
    } else {
      img[i] = new Features((*bufferDest)[i]);
    }

  if (fichierProfondeur != NULL) {
    _TIFFfree ((tdata_t) *bufferSourceProf);
    _TIFFfree ((tdata_t) *bufferDestProf);
    free (bufferSourceProf);
    free (bufferDestProf);
  }
  _TIFFfree ((tdata_t) *bufferSource);
  _TIFFfree ((tdata_t) *bufferDest);
  free (bufferSource);
  free (bufferDest);
  precedente = NULL; suivante = NULL;
}

Image::Image(char* nomFichier, char* nomFichierProfondeur) {
  TIFF* fichier = TIFFOpen (nomFichier, "r");
  TIFF* fichierProfondeur = NULL;
  if (nomFichierProfondeur != NULL) {
    fichierProfondeur = TIFFOpen (nomFichierProfondeur, "r");
  }

  uint32** bufferSource=(uint32**)malloc(sizeof(uint32*));
  LoadPixels (fichier, bufferSource, &largeur, &hauteur);
  uint32** bufferDest=(uint32**)malloc(sizeof(uint32*));
  SwapBuffer (bufferDest, *bufferSource, largeur, hauteur);

  uint32** bufferSourceProf = NULL;
  uint32** bufferDestProf = NULL;
  if (fichierProfondeur != NULL) {
    bufferSourceProf = (uint32**)malloc(sizeof(uint32*));
    LoadPixels (fichierProfondeur, bufferSourceProf, &largeur, &hauteur);
    bufferDestProf = (uint32**)malloc(sizeof(uint32*));
    SwapBuffer (bufferDestProf, *bufferSourceProf, largeur, hauteur);
  }

  img = new Features *[hauteur * largeur];

  for (uint32 i = 0 ; i<hauteur*largeur ; i++) 
    if (fichierProfondeur != NULL) {
      uint8 D =0.3*(((*bufferDestProf)[i] & 0x000000FF))+0.59*(((*bufferDestProf)[i] & 0x0000FF00)>> 8)+0.11*(((*bufferDestProf)[i] & 0x00FF0000)>> 16);
      img[i] = new Features((*bufferDest)[i], D);
    } else {
      img[i] = new Features((*bufferDest)[i]);
    }

  if (fichierProfondeur != NULL) {
    _TIFFfree ((tdata_t) *bufferSourceProf);
    _TIFFfree ((tdata_t) *bufferDestProf);
    free (bufferSourceProf);
    free (bufferDestProf);
    TIFFClose (fichierProfondeur);
  }
  _TIFFfree ((tdata_t) *bufferSource);
  _TIFFfree ((tdata_t) *bufferDest);
  free (bufferSource);
  free (bufferDest);
  precedente = NULL; suivante = NULL;

  TIFFClose (fichier);
}

// ça sert à quoi ça ?
Image::Image(Image* source, int ratio) {
  img = NULL;
  hauteur = (source->getHauteur(0)-1)/2;
  largeur = (source->getLargeur(0)-1)/2;
  precedente = source;
  suivante = NULL;
}

Image::~Image () {
  if (suivante != NULL) {
    delete suivante;
    suivante = NULL;
  }
  for (uint32 i=0 ; i<hauteur*largeur ; i++) {
    delete img[i];
  }
  delete [] img;
}


uint32 Image::getHauteur (int l) {
  if(l==0)
    return hauteur;
  else
    return suivante->getHauteur(l-1);
}

uint32 Image::getLargeur (int l) {
  if(l==0)
    return largeur;
  else
    return suivante->getLargeur(l-1);
}

Features* Image::getFeature(int x, int y,int l) {
  if(l==0)
    return img[y*largeur+x];
  else
    return suivante->getFeature(x,y,l-1);
}

int SubSize (int originalSize, int level) {
  int res = originalSize;
  int i;
  for (i=level ; i>0 ; i--) {
    res = (res-1)/2;
  }
  return res;
}

int Image::sauvegarderImageCouleur(TIFF* fichier,int l) {
  if(l==0)
    return SavePixels (fichier, img, largeur, hauteur);
  else 
    return suivante->sauvegarderImageCouleur(fichier,l-1);
}

int Image::sauvegarderImageCouleur(char* nomFichier,int l) {
  if(l==0)
    {
      TIFF* fichier = TIFFOpen (nomFichier, "w");
      int res;
      res = SavePixels (fichier, img, largeur, hauteur);
      TIFFClose (fichier);
      return res;
    }
  else
    return suivante->sauvegarderImageCouleur(nomFichier,l-1);
}




int Image::LoadPixels (TIFF* file, uint32** memory, uint32* width, uint32* height)
{
  int res;
  
  /* Get the width & height of the image */

  res = TIFFGetField (file, TIFFTAG_IMAGEWIDTH, width);
  if (res != 1) {
    printf ("Tag IMAGEWIDTH was not found in tiff file\n");
    return 1;
  }
  res = TIFFGetField (file, TIFFTAG_IMAGELENGTH, height);
  if (res != 1) {
    printf ("Tag IMAGEHEIGHT was not found in tiff file\n");
    return 1;
  }
  //printf ("TIFF file has an image of size %ix%i.\n", (int)(*width), (int)(*height));

  /* Allocate memory for the pixels and load them from file */

  *memory = (uint32*)_TIFFmalloc ((*width) * (*height) * sizeof(uint32));

  if (*memory == (uint32*)NULL) {
    printf ("Could not allocate memory for image pixels.\n");
    return 1;
  }
  res = TIFFReadRGBAImage (file, *width, *height, *memory, 0);
  if (res != 1) {
    printf ("Could not load image pixels into memory\n");
    _TIFFfree ((tdata_t)(*memory));
    *memory = (uint32*)NULL;
    return 1;
  }

  return 0;
}


int Image::SwapBuffer (uint32** memoryD, uint32* memoryS, uint32 width, uint32 height)
{
  int        i;
  uint32*    swap;
  uint32*    cursS;
  uint32*    cursD;
  uint32     pixel;
  uint32     x, y;
	
  swap = (uint32*)_TIFFmalloc (width * height * sizeof(uint32));
  if (swap == (uint32*)NULL) {
    printf ("Could not allocate swap buffer\n");
    return 1;
  }

  cursS = memoryS;
  cursD = swap + (height - 1) * width;
	
  /* Il faut tester le format de la machine : little ou big endian. */
  i = 0;
  *((char *)(&i)) = 1;

  if (i == 1)
  {
    /* En little endian, il ne faut pas echanger les canaux avant la sauvegarde
       car en memoire (PAS dans la variable pixel, mais en representation interne),
       les informations sont rangees dans l'ordre RVBX. 
    */
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        *cursD = *cursS;
        cursS++;
        cursD++;
      }
      cursD -= 2 * width;
    }
  }
  else 
  {
    /* En big endian, il faut echanger les canaux avant la sauvegarde car en memoire,
       les informations sont rangees dans l'ordre XBVR. 
    */
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        pixel = *cursS;
        *cursD = ((pixel & 0xFF000000) >> 24) |
                 ((pixel & 0x00FF0000) >>  8) |
                 ((pixel & 0x0000FF00) <<  8) |
                 ((pixel & 0x000000FF) << 24);
        cursS++;
        cursD++;
      }
      cursD -= 2 * width;
    }
  }

  *memoryD = swap;
  return 0;
}



int Image::SavePixels (TIFF* file, Features** img, uint32 width, uint32 height)
{
  int     res;
  int     dataSize;



  uint32* memory = (uint32*)malloc(height*width*sizeof(uint32));
  for (uint32 i = 0 ; i<height*width ; i++) 
    memory[i] = img[i]->toPixel();
  
  res =   TIFFSetField (file, TIFFTAG_IMAGEWIDTH,     width);
  if (res == 1)
    res = TIFFSetField (file, TIFFTAG_IMAGELENGTH,    height);
  if (res == 1)
    res = TIFFSetField (file, TIFFTAG_BITSPERSAMPLE,  (uint16)8);
  if (res == 1)
    res = TIFFSetField (file, TIFFTAG_PHOTOMETRIC,    (uint16)2);
  if (res == 1)
    res = TIFFSetField (file, TIFFTAG_SAMPLESPERPIXEL,(uint16)4);
  if (res == 1)
    res = TIFFSetField (file, TIFFTAG_PLANARCONFIG,
                                              (uint16)PLANARCONFIG_CONTIG);
	if (res == 1) {
    dataSize = sizeof(uint32) * width * height;
    res = (int)TIFFWriteRawStrip (file, 0, memory, dataSize);
    if (res == -1)
      res = 0;
  }
  if (res == 0) {
    printf ("An error occured while trying to save TIFF file\n");
    return 1;
  }

  free(memory);

  return 0;
}


void Image::setFeature(int x,int y,int l,Features* f)
{  
  if(l==0) {
    if (img[y*largeur+x] == NULL) {
      img[y*largeur+x]=new Features(*f);
    } else {
      img[y*largeur+x]->setR(f->getR());
      img[y*largeur+x]->setG(f->getG());
      img[y*largeur+x]->setB(f->getB());
      img[y*largeur+x]->setX(f->getX());
      img[y*largeur+x]->setD(f->getD());
      img[y*largeur+x]->setLUM(f->getLUM());
      img[y*largeur+x]->setI(f->getI());
      img[y*largeur+x]->setQ(f->getQ());
    }
  } else {
    suivante->setFeature(x,y,l-1,f);  
  }
}

int Image::calculerLesPyramides(int nombreDeNiveaux) {
  int res;
  if (suivante != NULL) {
    delete suivante;
  }
  if (nombreDeNiveaux > 0) {
    int demieLargeur = (this->largeur-1)/2;
    int demieHauteur = (this->hauteur-1)/2;
    suivante = new Image (demieHauteur, demieLargeur);
    for (int xCourant = 0 ; xCourant < demieLargeur ; xCourant++) {
      for (int yCourant = 0 ; yCourant < demieHauteur ; yCourant++) {
	uint8 r=(this->getFeature(2*xCourant, 2*yCourant, 0)->getR())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant, 0)->getR())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant, 0)->getR())/16
	  +(this->getFeature(2*xCourant, 2*yCourant+1, 0)->getR())/8
	  +(this->getFeature(2*xCourant+1, 2*yCourant+1, 0)->getR())/4 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+1, 0)->getR())/8
	  +(this->getFeature(2*xCourant, 2*yCourant+2, 0)->getR())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant+2, 0)->getR())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+2, 0)->getR())/16;
	uint8 g=(this->getFeature(2*xCourant, 2*yCourant, 0)->getG())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant, 0)->getG())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant, 0)->getG())/16
	  +(this->getFeature(2*xCourant, 2*yCourant+1, 0)->getG())/8
	  +(this->getFeature(2*xCourant+1, 2*yCourant+1, 0)->getG())/4 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+1, 0)->getG())/8
	  +(this->getFeature(2*xCourant, 2*yCourant+2, 0)->getG())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant+2, 0)->getG())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+2, 0)->getG())/16;
	uint8 b=(this->getFeature(2*xCourant, 2*yCourant, 0)->getB())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant, 0)->getB())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant, 0)->getB())/16
	  +(this->getFeature(2*xCourant, 2*yCourant+1, 0)->getB())/8
	  +(this->getFeature(2*xCourant+1, 2*yCourant+1, 0)->getB())/4 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+1, 0)->getB())/8
	  +(this->getFeature(2*xCourant, 2*yCourant+2, 0)->getB())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant+2, 0)->getB())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+2, 0)->getB())/16;
	uint8 x=(this->getFeature(2*xCourant, 2*yCourant, 0)->getX())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant, 0)->getX())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant, 0)->getX())/16
	  +(this->getFeature(2*xCourant, 2*yCourant+1, 0)->getX())/8
	  +(this->getFeature(2*xCourant+1, 2*yCourant+1, 0)->getX())/4 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+1, 0)->getX())/8
	  +(this->getFeature(2*xCourant, 2*yCourant+2, 0)->getX())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant+2, 0)->getX())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+2, 0)->getX())/16;
	uint8 d=(this->getFeature(2*xCourant, 2*yCourant, 0)->getD())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant, 0)->getD())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant, 0)->getD())/16
	  +(this->getFeature(2*xCourant, 2*yCourant+1, 0)->getD())/8
	  +(this->getFeature(2*xCourant+1, 2*yCourant+1, 0)->getD())/4 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+1, 0)->getD())/8
	  +(this->getFeature(2*xCourant, 2*yCourant+2, 0)->getD())/16
	  +(this->getFeature(2*xCourant+1, 2*yCourant+2, 0)->getD())/8 
	  +(this->getFeature(2*xCourant+2, 2*yCourant+2, 0)->getD())/16;

	Features* tmp = new Features(r,g,b,x,d);
	suivante->setFeature(xCourant,yCourant,0,tmp);
	delete tmp;
      }
    }
    suivante->precedente = this;
    res = suivante->calculerLesPyramides (nombreDeNiveaux-1);
  } else {
    res = 0;
  }
  return res;
}



int Image::computeRGB(Image* colorImage)
{
  if (this->hauteur != colorImage->hauteur || this->largeur != colorImage->largeur) {
    cout << "Impossible d'utiliser computeRGB sur deux images de taille différentes" << endl;
    return -1;
  } else {
    for (uint32 x=0 ; x<largeur ; x++) {
      for (uint32 y=0 ; y<hauteur ; y++) {
	img[y*largeur+x]->computeRGB(colorImage->getFeature(x,y,0));
      }
    }
  }
  return 0;
}

void  Image::luminanceRemapping(Image* B)
{
  int ma,mb,va,vb;
  ma=moyenneLum();
  va=varianceLum(ma);
  mb=B->moyenneLum();
  vb=B->varianceLum(mb);
  
  for (uint32 x=0 ; x<largeur ; x++)
    for (uint32 y=0 ; y<hauteur ; y++)
      img[y*largeur+x]->setLUM((vb/va)*(img[y*largeur+x]->getLUM()-ma)+mb);
}

int  Image::moyenneLum()
{
  int s=0;
  for (uint32 x=0 ; x<largeur ; x++)
    for (uint32 y=0 ; y<hauteur ; y++)
      s+=img[y*largeur+x]->getLUM();
  return s/(largeur*hauteur);   
}

int  Image::varianceLum(int moy)
{
  int s=0;
  for (uint32 x=0 ; x<largeur ; x++)
    for (uint32 y=0 ; y<hauteur ; y++)
      s+=pow(img[y*largeur+x]->getLUM()-moy,2);
  return sqrt(s/(largeur*hauteur));   

}


Point* Image::ComputeBruteForce (Image* imageAp, Image* imageB, Image* imageBp, Point* q, int level, int gk, FeatureType type) 
{
  if (level > 0) {
    Point* p = this->suivante->ComputeBruteForce(imageAp->suivante, imageB->suivante, imageBp->suivante, q, level-1, gk, type);
    return p;
  } else {
    Point* p = new Point(0,0);
    int min = INT_MAX;
    int xB = q->getX(); int yB = q->getY();
    int largeurB = imageB->getLargeur(level);
    int hauteurB = imageB->getHauteur(level);
    int largeurA = this->getLargeur(level);
    int hauteurA = this->getHauteur(level);
    Features** imgB = imageB->img;
    Features** imgBp = imageBp->img;
    int Xleft = -gk; if (xB+Xleft < 0) {Xleft = -xB;}
    int Yup = -gk; if (yB+Yup < 0) {Yup = -yB;}
    int Xright = gk; if (xB+Xright > largeurB-1) {Xright = largeurB-1-xB;}
    int Ydown = gk; if (yB+Ydown > hauteurB-1) {Ydown = hauteurB-1-yB;}
    Features** imgAp = imageAp->img;
   
    switch (type){
    case RGBLUM :
      // on utilise la luminance ET les canaux RGB
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 rp;
	  uint8 rq;
	  uint8 gp;
	  uint8 gq;
	  uint8 bp;
	  uint8 bq;
	  uint8 lump;
	  uint8 lumq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = this->img[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgB[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = this->img[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgB[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = this->img[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgB[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	      lump = this->img[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgB[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	      lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	    rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	    rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	    res += (rp-rq)*(rp-rq);
	    gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	    gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	    res += (gp-gq)*(gp-gq);
	    bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	    bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	    res += (bp-bq)*(bp-bq);
	    lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	    lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	    res += (lump-lumq)*(lump-lumq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;

    case RGB :
      // on utilise les canaux RGB
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 rp;
	  uint8 rq;
	  uint8 gp;
	  uint8 gq;
	  uint8 bp;
	  uint8 bq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = this->img[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgB[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = this->img[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgB[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = this->img[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgB[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	    rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	    rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	    res += (rp-rq)*(rp-rq);
	    gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	    gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	    res += (gp-gq)*(gp-gq);
	    bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	    bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	    res += (bp-bq)*(bp-bq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;
      
    case LUMIN :
      // on utilise seulement la luminance
      
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 lump;
	  uint8 lumq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      lump = this->img[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgB[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	      lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;

    case IQ :
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 ip;
	  uint8 iq;
	  uint8 qp;
	  uint8 qq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      ip = this->img[(yA+y)*largeurA+(xA+x)]->getI();
	      iq = imgB[(yB+y)*largeurB+(xB+x)]->getI();
	      res += (ip-iq)*(ip-iq);
	      qp = this->img[(yA+y)*largeurA+(xA+x)]->getQ();
	      qq = imgB[(yB+y)*largeurB+(xB+x)]->getQ();
	      res += (qp-qq)*(qp-qq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      ip = imgAp[(yA+y)*largeurA+(xA+x)]->getI();
	      iq = imgBp[(yB+y)*largeurB+(xB+x)]->getI();
	      res += (ip-iq)*(ip-iq);
	      qp = imgAp[(yA+y)*largeurA+(xA+x)]->getQ();
	      qq = imgBp[(yB+y)*largeurB+(xB+x)]->getQ();
	      res += (qp-qq)*(qp-qq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	    ip = imgAp[(yA+y)*largeurA+(xA+x)]->getI();
	    iq = imgBp[(yB+y)*largeurB+(xB+x)]->getI();
	    res += (ip-iq)*(ip-iq);
	    qp = imgAp[(yA+y)*largeurA+(xA+x)]->getQ();
	    qq = imgBp[(yB+y)*largeurB+(xB+x)]->getQ();
	    res += (qp-qq)*(qp-qq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;

    case DEPRGB:
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 rp;
	  uint8 rq;
	  uint8 gp;
	  uint8 gq;
	  uint8 bp;
	  uint8 bq;
	  uint8 dp;
	  uint8 dq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = this->img[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgB[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = this->img[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgB[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = this->img[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgB[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	      dp = this->img[(yA+y)*largeurA+(xA+x)]->getD();
	      dq = imgB[(yB+y)*largeurB+(xB+x)]->getD();
	      res += (dp-dq)*(dp-dq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	      rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	      res += (rp-rq)*(rp-rq);
	      gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	      gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	      res += (gp-gq)*(gp-gq);
	      bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	      bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	      res += (bp-bq)*(bp-bq);
	      dp = imgAp[(yA+y)*largeurA+(xA+x)]->getD();
	      dq = imgBp[(yB+y)*largeurB+(xB+x)]->getD();
	      res += (dp-dq)*(dp-dq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	    rp = imgAp[(yA+y)*largeurA+(xA+x)]->getR();
	    rq = imgBp[(yB+y)*largeurB+(xB+x)]->getR();
	    res += (rp-rq)*(rp-rq);
	    gp = imgAp[(yA+y)*largeurA+(xA+x)]->getG();
	    gq = imgBp[(yB+y)*largeurB+(xB+x)]->getG();
	    res += (gp-gq)*(gp-gq);
	    bp = imgAp[(yA+y)*largeurA+(xA+x)]->getB();
	    bq = imgBp[(yB+y)*largeurB+(xB+x)]->getB();
	    res += (bp-bq)*(bp-bq);
	    dp = imgAp[(yA+y)*largeurA+(xA+x)]->getD();
	    dq = imgBp[(yB+y)*largeurB+(xB+x)]->getD();
	    res += (dp-dq)*(dp-dq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;
      
    case DEPLUM:
      for (int yA = -Yup ; yA < hauteurA-Ydown ; yA++) {
	for (int xA = -Xleft ; xA < largeurA-Xright ; xA++) {
	  int res = 0;
	  uint8 lump;
	  uint8 lumq;
	  uint8 dp;
	  uint8 dq;
	  int x; int y;
	  for (y = Yup ; y <= Ydown ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      lump = this->img[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgB[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	      dp = this->img[(yA+y)*largeurA+(xA+x)]->getD();
	      dq = imgB[(yB+y)*largeurB+(xB+x)]->getD();
	      res += (dp-dq)*(dp-dq);
	    }
	  }
	  if (res > min) {continue;}
	  for (y = Yup ; y < 0 ; y++) {
	    for (x = Xleft ; x <= Xright ; x++) {
	      lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	      lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	      res += (lump-lumq)*(lump-lumq);
	      dp = imgAp[(yA+y)*largeurA+(xA+x)]->getD();
	      dq = imgBp[(yB+y)*largeurB+(xB+x)]->getD();
	      res += (dp-dq)*(dp-dq);
	    }
	  }
	  if (res > min) {continue;}
	  y = 0;
	  for (x = Xleft ; x < 0; x++) {
	    lump = imgAp[(yA+y)*largeurA+(xA+x)]->getLUM();
	    lumq = imgBp[(yB+y)*largeurB+(xB+x)]->getLUM();
	    res += (lump-lumq)*(lump-lumq);
	    dp = imgAp[(yA+y)*largeurA+(xA+x)]->getD();
	    dq = imgBp[(yB+y)*largeurB+(xB+x)]->getD();
	    res += (dp-dq)*(dp-dq);
	  }
	  if (res < min) {
	    min = res;
	    p->setX(xA);
	    p->setY(yA);
	    if (res == 0) {
	      return p;
	    }
	  }
	}
      }
      break;
      
    default:
      break;
      
    }
    
    return p;
  }
}
