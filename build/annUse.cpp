//----------------------------------------------------------------------
//      File:           annUse.cpp
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    classe permettant d'utiliser l'algorithme ANN
//                      pour l'analogie d'image
//----------------------------------------------------------------------

#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include "annUse.h"
#include "Point.h"


//construit une instance de la classe a partir des images sources A et Ap
//à un level l, on cherchera une correspondance en utilisant un vosinage mgk
//au niveau l et mpk au niveau inferrieur
annUse::annUse(Image* A,Image* Ap,int l, int mpk,int mgk, FeatureType Type)
{
  pk=mpk;
  gk=mgk;
  int m,M;

  uint32 largeurA = A->getLargeur(l);
  if (largeurA/2 != (largeurA-1)/2) {
    largeurA--;
  }
  uint32 hauteurA = A->getHauteur(l);
  if (hauteurA/2 != (hauteurA-1)/2) {
    hauteurA--;
  }

  if(gk>=2*pk) {
    m=gk;
    M=gk+1;
  } else {
    m=pk*2;
    M=pk*2+1;
  }
  int gK=gk*2+1;
  int pK=pk*2+1;
  Features f;
  //calcul de la taille d'un point en fonction de la taille du voisinnage
  dim=(gK*gK+pK*pK*2+(gK*gK)/2)*f.getNBChamps(Type);
  eps=0;
  
  //calcul du nombre de points
  nPts= (largeurA-m-M)*(hauteurA-m-M);

  //allocation des ressources
  dataPts = annAllocPts(nPts, dim);
  nnIdx = new ANNidx[1];
  dists = new ANNdist[1];

  //remplissage du tableau dataPts
  int i=0;
  for(uint32 y=m; y<hauteurA-M  ; y++)
    for(uint32 x=m; x<largeurA-M; x++)
      {
	makeANNPoint(A,Ap,l,x,y,Type,dataPts[i]);
	i++;
      }
	  
  //construction de l'arbre de recherche
  kdTree = new ANNkd_tree( dataPts, nPts, dim);    
}

//destruction d'une instance de annUse
annUse::~annUse()
{
  //destruction des ressources utilisées par l'ann
  annDeallocPts(dataPts);
  delete [] nnIdx;
  delete [] dists;
  delete kdTree;
  annClose();
}

//retourne le point de p correspondant au point q des images B et Bp dans A et Ap
//au niveau l
Point* annUse::bestApproximateMatch(Image* B, Image* Bp,int l, Point* q,Image* A, FeatureType Type)
{
  int marge;
  if(gk==0 && pk==0)
    marge=2;
  else
    marge=(gk>=pk*2)?gk*3:pk*4+2;
  //construction de l'ANNpoint correspondant à q et son voisinnage
  ANNpoint queryPt = annAllocPt(dim);
  makeANNPoint(B,Bp,l,q->getX(),q->getY(), Type, queryPt);

  //recherche du point correspondant dans A et Ap
  kdTree->annkSearch(queryPt, 1, nnIdx, dists, eps);

  //calcul et retour des coordonnées du point
  uint32 x=nnIdx[0]%(A->getLargeur(l)-marge);
  uint32 y=nnIdx[0]/(A->getLargeur(l)-marge);

  //libération de l'ANNpoint
  annDeallocPt(queryPt);

  return new Point(x+gk,y+gk);
}


//construit un ANNpoint pour le point (x,y) des images I et Ip
ANNpoint annUse::makeANNPoint(Image* I,Image* Ip,int l,uint32 x, uint32 y, FeatureType Type, ANNpoint p)
{
  Features f;
  
  //ajout des valeurs de I au niveau l
  int j=0;
  for(uint32 y1=y-gk ; y1<=y+gk ; y1++)
    for(uint32 x1=x-gk ; x1<=x+gk ; x1++)
      {
	I->getFeature(x1,y1,l)->add(p,j,Type);
	j+=f.getNBChamps(Type);
      }

  //ajout des valeurs de Ip au niveau l
  for(uint32 y1=y-gk;y1<=y+gk;y1++) {
    uint32 x1;
    for(x1=x-gk;x1<=x+gk;x1++) {
      if(x1==x && y1==y)
	break;
      Ip->getFeature(x1,y1,l)->add(p,j,Type);
      j+=f.getNBChamps(Type);
    }
    if(x1==x && y1==y)
      break;
  }
  
  uint xp=x/2;
  uint yp=y/2;
  //ajout des valeurs de I au niveau l-1
  for(uint32 y1=yp-pk;y1<=yp+pk;y1++)
    for(uint32 x1=xp-pk;x1<=xp+pk;x1++) {
      I->getFeature(x1,y1,l+1)->add(p,j,Type);
      j+=f.getNBChamps(Type);
    }
  
  //ajout des valeurs de Ip au niveau l-1
  for(uint32 y1=yp-pk;y1<=yp+pk;y1++)
    for(uint32 x1=xp-pk;x1<=xp+pk;x1++) {
      Ip->getFeature(x1,y1,l+1)->add(p,j,Type);
      j+=f.getNBChamps(Type);
    }
  return p;
}

