//----------------------------------------------------------------------
//      File:           BestCoherence.h
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Calcul du BestMatch
//----------------------------------------------------------------------

#include "classeFeatures.h"
#include "loadSaveTiffClass.h"
#include "Point.h"
#include "tiffio.h"
#include "annUse.h"
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

//return la norme de la difference de deux vecteur, 
//-1 si les deux vecteur n'ont pas la même dimension
int Norme(vector<uint8> f1,vector<uint8> f2);


Point* BestCoherence(Image* A1,Image* A2, Image* B1, Image* B2, Point** s, int l, Point* q, int pk,int gk, FeatureType Type);

Point* BruteForceSearch(Image* A,Image* Ap, Image* B, Image* Bp,int l, Point* q);

int indS(int l,int x,int y,Image *I);

Point* BestMatch(Image *A,Image *Ap,Image * B,Image * Bp,int l,Point *q,Point** s,annUse *ann, int pk,int gk,int level, float K, FeatureType Type);
