//----------------------------------------------------------------------
//      File:           classeFeatures.h
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Deffinition d'un feature associé à un point d'une image
//----------------------------------------------------------------------

#ifndef __CLASSEFEATURES__
#define __CLASSEFEATURES__

enum FeatureType {RGB, LUMIN, RGBLUM, IQ, DEPRGB, DEPLUM };

#include "tiffio.h"
#include "annUse.h"
#include "ann/include/ANN/ANN.h"
#include <stdio.h>
#include <iostream>
#include <vector>


using namespace std;

class Features {

 public:
  Features();
  Features(uint32 p, uint8 d=0);
  Features(Features & f);
  Features(uint8 mR,uint8 mG,uint8 mB, uint8 mX, uint8 mD=0);
  uint32 toPixel();
  uint8 getR ();
  uint8 getG ();
  uint8 getB ();
  uint8 getX ();
  uint8 getD ();
  uint8 getLUM ();
  uint8 getI();
  uint8 getQ();
  static uint8 getNBChamps (FeatureType type);
  void setR (uint8 mr);
  void setG (uint8 mg);
  void setB (uint8 mb);
  void setX (uint8 mx);
  void setD (uint8 md);
  void setLUM (uint8 ml);
  void setI(uint8 mi);
  void setQ(uint8 mq);

  void computeRGB(Features* colorFeature);
  void computeRGB(uint8 Rb, uint8 Gb, uint8 Bb, uint8 Xb);

  void add(vector <uint8>*  v, FeatureType Type);
  void add(ANNpoint p, int j, FeatureType Type);

 private:
  uint8 R;
  uint8 G;
  uint8 B;
  uint8 X;
  uint8 D;
  uint8 LUM;
  uint8 I;
  uint8 Q;

  void computeYIQ();

};

#endif
