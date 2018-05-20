//----------------------------------------------------------------------
//      File:           annUse.h
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    classe permettant d'utiliser l'algorithme ANN
//                      pour l'analogie d'image
//----------------------------------------------------------------------
#ifndef __ANNUSE__
#define __ANNUSE__

#include "tiffio.h"
#include "loadSaveTiffClass.h"
#include <stdio.h>
#include "Point.h"
#include "ann/include/ANN/ANN.h"

class Image;
class Point;

class annUse {
 public:
  //construit une instance de la classe a partir des images sources A et Ap
  //à un level l, on cherchera une correspondance en utilisant un vosinage mgk
  //au niveau l et mpk au niveau inferrieur
  annUse(Image* A,Image* App,int l, int mpk,int mgk, FeatureType Type);

  //retourne le point de p correspondant au point q des images B et Bp dans A et Ap
  //au niveau l
  Point * bestApproximateMatch(Image* B, Image* Bp,int l, Point* q,Image* A, FeatureType Type);

  //destruction d'une instance de annUse
  ~annUse();

 private:

  //taille des ANNpoints
  int dim;
  //marge d'erreur tolérée
  double eps;
  //taille des voisinnages
  int gk,pk;

  //nombre de points de l'espace de recherche
  int nPts;

  //données de l'ANN
  ANNpointArray	dataPts;
  ANNidxArray nnIdx;
  ANNdistArray dists;
  ANNkd_tree* kdTree;	

  //construit un ANNpoint pour le point (x,y) des images I et Ip
  ANNpoint makeANNPoint(Image* I,Image* Ip,int l,uint32 x, uint32 y, FeatureType Type, ANNpoint p);
};


#endif
