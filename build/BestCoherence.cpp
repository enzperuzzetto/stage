//----------------------------------------------------------------------
//      File:           BestCoherence.cpp
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Calcul du BestMatch
//----------------------------------------------------------------------
#include "BestCoherence.h"
#include <stdio.h>
#include <stdlib.h>
#include <climits>

using namespace std;

static int eps=0;

/*return la norme de la difference de deux vecteur, 
  -1 si les deux vectereturnur n'ont pas la même dimension*/
int Norme(vector<uint8> f1,vector<uint8> f2){
  int res=0;
  if(f1.size()!=f2.size()){
    return -1;
  }else{
    vector<uint8>::iterator it1 = f1.begin(); 
    vector<uint8>::iterator it2 = f2.begin();
    while(it1 != f1.end()){
      res+=((*it1)-(*it2))*((*it1)-(*it2));
      it1++;
      it2++;	
    }
    return res;
  }
}



void makeF(Image* I, Image* Ip, int l, uint32 x, uint32 y, vector<uint8> * F, int gk, int pk, FeatureType Type)
{
  F->clear();

  for(uint32 y1=y-gk;y1<=y+gk;y1++)
    for(uint32 x1=x-gk;x1<=x+gk;x1++)
      I->getFeature(x1,y1,l)->add(F, Type);
  

  for(uint32 y1=y-gk;y1<=y+gk;y1++)
    {
      uint32 x1;
      for(x1=x-gk;x1<=x+gk;x1++)
	{
	  if(x1==x && y1==y)
	    break;
	  Ip->getFeature(x1,y1,l)->add(F, Type);
	}
      if(x1==x && y1==y)
	break;
    }

  uint xp=x/2;
  uint yp=y/2;

  for(uint32 y1=yp-pk ; y1<=yp+pk ; y1++)
    for(uint32 x1=xp-pk ; x1<=xp+pk ; x1++) {
      I->getFeature(x1,y1,l+1)->add(F, Type);
    }
	      

  for(uint32 y1=yp-pk ; y1<=yp+pk ; y1++)
    for(uint32 x1=xp-pk ; x1<=xp+pk ; x1++)
      Ip->getFeature(x1,y1,l+1)->add(F, Type);

  
}




int indS(int l, int x,int y,Image *I)
{
  return y*I->getLargeur(l)+x;
}

Point* BestCoherence(Image* A1,Image* A2, Image* B1, Image* B2, Point** s, int l, Point* q, int gk, int pk, FeatureType Type)
{
  uint32 x=q->getX();
  uint32 y=q->getY();
  int n;
  int min=INT_MAX;
  Point* minr=NULL;
  vector<uint8> F1;
  vector<uint8> F2;
  
  for(uint32 y1=y-gk;y1<=y+gk;y1++)
    {
      uint32 x1;
      for(x1=x-gk;x1<=x+gk;x1++)
	{
	  if(x1==x && y1==y)
	    break;
	  Point* r = new Point(x1,y1);
	  Point* a1=new Point(s[indS(l,x1,y1,B2)]->getX()+(q->getX()-r->getX()),
			      s[indS(l,x1,y1,B2)]->getY()+(q->getY()-r->getY()));

	  if(a1->getX()<(uint32)gk || a1->getY()<(uint32)gk) {
	    delete r;
	    delete a1;	    
	    continue;
	  }
	  if(a1->getX()+gk>=A1->getLargeur(l) || a1->getY()+gk>=A1->getHauteur(l)) {
	    delete r;
	    delete a1;
	    continue;	
	  }  
	  if(a1->getX()/2+pk>=A1->getLargeur(l+1) || a1->getY()/2+pk>=A1->getHauteur(l+1)) {
	    delete r;
	    delete a1;
	    continue;
	  }

	  makeF(A1,A2,l,a1->getX(),a1->getY(),&F2,gk,pk, Type);
	  makeF(B1,B2,l,x,y,&F1,gk,pk, Type);
	  
	  n=Norme(F1,F2);
	  if (n<min) {
	    if (minr != NULL) {
	      delete minr;
	    }
	    min=n;
	    minr=new Point(*a1);
	  }
	  delete r;
	  delete a1;
	}
      if(x1==x && y1==y)
	break;
    }
  return minr;

}



void makeF(Image* I,Image* Ip,int l,uint32 x, uint32 y,vector<uint8> * F,Point* pMin,Point* pMax , FeatureType Type)
{
  F->clear();

  for(uint32 y1=y-pMin->getY();y1<=y+pMax->getY();y1++)
    for(uint32 x1=x-pMin->getX();x1<=x+pMax->getX();x1++)
      I->getFeature(x1,y1,l)->add(F, Type);

  for(uint32 y1=y-pMin->getY();y1<=y;y1++)	
    for(uint32 x1=x-pMin->getX();x1<=x+pMax->getX();x1++)
      {
	if(x1==x && y1==y)
	  break;
	Ip->getFeature(x1,y1,l)->add(F,Type);
      }  
}

Point* BruteForceSearch(Image* A,Image* Ap, Image* B, Image* Bp,int l, Point* q, FeatureType Type)
{
  uint32 k=2;
  Point* pMin=new Point(0,0);
  Point* pMax=new Point(0,0);

  if(k <= q->getX()) {
    pMin->setX(k);
  } else {
    pMin->setX(q->getX());
  }

  if(k <= q->getY()) {
    pMin->setY(k);
  } else {
    pMin->setY(q->getY());
  }

  if(q->getX() < B->getLargeur(l)-k) {
    pMax->setX(k);
  } else {
    pMax->setX(B->getLargeur(l)-q->getX()-(k-1));
  }

  if(q->getY() < B->getHauteur(l)-k) {
    pMax->setY(k);
  } else {
    pMax->setY(B->getHauteur(l)-q->getY()-(k-1));
  }

  vector<uint8> Vq;
  vector<uint8> Vp;
  makeF(B,Bp,l,q->getX(),q->getY(),&Vq,pMin,pMax, Type);

  int min = INT_MAX;
  Point* pres = new Point(0,0);
  int nc;
  
  for(uint32 y=0;y<A->getHauteur(l);y++) {
    for(uint32 x=0;x<A->getLargeur(l);x++) {
	if(x >= pMin->getX() && y>=pMin->getY())
	  if(x+pMax->getX() < A->getLargeur(l) && y+pMax->getY()<A->getHauteur(l)) {
	    makeF(A,Ap,l,x,y,&Vp,pMin,pMax,Type);
	    nc = Norme(Vp,Vq);
	    if(nc < min && nc != -1) {
	      min=nc;
	      pres->setX(x);
	      pres->setY(y);
	      if(nc<eps) {
		delete pMin;
		delete pMax;
		return pres;
	      }
	    }
	  }
    }
  }
  delete pMin;
  delete pMax;
  return pres;
}


Point* BestMatch(Image *A, Image *Ap, Image * B, Image * Bp, int l, Point *q, Point** s, annUse *ann, int pk, int gk, int level, float K, FeatureType Type)
{  
  if((q->getX()<(uint32)gk || q->getY()<(uint32)gk) ||
     (q->getX()+gk>=B->getLargeur(l) || q->getY()+gk>=B->getHauteur(l)) ||
     (q->getX()/2<(uint32)pk || q->getY()/2<(uint32)pk) ||
     (q->getX()/2+pk>=B->getLargeur(l+1) || q->getY()/2+pk>=B->getHauteur(l+1)))
    return A->ComputeBruteForce(Ap, B, Bp, q, l, gk, Type);

  Point* papp=ann->bestApproximateMatch(B,Bp,l,q,A,Type);
  Point* pcoh=BestCoherence(A,Ap,B,Bp,s,l,q,gk,pk,Type);

  if(pcoh==NULL)
    return papp;
  
  if(pcoh->getX() < (uint32)gk || pcoh->getY() < (uint32)gk || pcoh->getX()+gk >= A->getLargeur(l) || pcoh->getY()+gk >= A->getHauteur(l) || pcoh->getX()/2 < (uint32)pk || pcoh->getY()/2 < (uint32)pk || pcoh->getX()/2+pk >= A->getLargeur(l+1) || pcoh->getY()/2+pk >= A->getHauteur(l+1)) {
    delete pcoh;
    return papp;
  }

  if(papp->getX() < (uint32)gk || papp->getY() < (uint32)gk || papp->getX()+gk >= A->getLargeur(l) || papp->getY()+gk >= A->getHauteur(l) || papp->getX()/2 < (uint32)pk || papp->getY()/2 < (uint32)pk || papp->getX()/2+pk >= A->getLargeur(l+1) || papp->getY()/2+pk >= A->getHauteur(l+1)) {
    delete papp;
    return pcoh;
  }
  
  vector<uint8> Fapp;
  vector<uint8> Fcoh;
  vector<uint8> Fq;
  makeF(A,Ap,l,papp->getX(),papp->getY(),&Fapp,gk,pk,Type);
  makeF(A,Ap,l,pcoh->getX(),pcoh->getY(),&Fcoh,gk,pk,Type);
  makeF(B,Bp,l,q->getX(),q->getY(),&Fq,gk,pk,Type);
  int dapp=Norme(Fapp,Fq);
  int dcoh=Norme(Fcoh,Fq);

  if(dcoh<=dapp*(1+K*pow(2,1+level-l))) {
    delete papp;
    return pcoh;
  } else {
    delete pcoh;
    return papp;
    }
}
