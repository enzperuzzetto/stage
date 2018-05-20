//----------------------------------------------------------------------
//      File:           Point.cpp
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Définnition d'un point à 2 dimensions
//----------------------------------------------------------------------

#include "Point.h"
#include<iostream>

using namespace std;

//contructeur du point de coordonnées (0,0)
Point::Point()
{
  x=0;
  y=0;
}

//contructeur d'un point par copie
Point::Point(Point & p)
{
  x=p.getX();
  y=p.getY();
}

//constructeur d'un point a partir de coordonnées
Point::Point(uint32 mx,uint32 my)
{
  x=mx;
  y=my;
}

uint32 Point::getX ()
{
  return x;
}

uint32 Point::getY ()
{
  return y;
}

Point &Point::operator+ (const Point &p)
{
  x += p.x;
  y += p.y;
  return *this;
  
}
Point &Point::operator- (const Point &p)
{   
  x -= p.x;
  y -= p.y;
  return *this;  
}

void Point::setX (uint32 mx)
{
  x=mx;
}

void Point::setY (uint32 my)
{
  y=my;
}
