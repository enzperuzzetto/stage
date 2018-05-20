//----------------------------------------------------------------------
//      File:           Point.h
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Définnition d'un point à 2 dimensions
//----------------------------------------------------------------------

#ifndef __POINT__
#define __POINT__

#include "tiffio.h"
#include <stdio.h>

class Point {

 public:
  //contructeur du point de coordonnées (0,0)
  Point();
  //constructeur d'un point a partir de coordonnées
  Point(uint32 mx,uint32 my);
  //contructeur d'un point par copie
  Point(Point & p);
  uint32 getX ();
  uint32 getY ();
  void setX (uint32 mx);
  void setY (uint32 my);


  Point &operator+ (const Point &p);
  Point &operator- (const Point &p);
 private:
  uint32 x;
  uint32 y;

};


#endif
