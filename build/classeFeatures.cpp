//----------------------------------------------------------------------
//      File:           classeFeatures.cpp 
//      Programmers:    Lamranih Aboud El Assad Hanane, Haboubi Jihene
//                      Siret Orianne, Fusade Rémi, Ducout Romain
//      Description:    Deffinition d'un feature associé à un point d'une image
//----------------------------------------------------------------------

#include "classeFeatures.h"
#include<iostream>

using namespace std;

Features::Features(){
  R=0;
  G=0;
  B=0;
  X=255;
  D=0;
  LUM=0;
  I=0;
  Q=0;
}

Features::Features(uint8 mR,uint8 mG,uint8 mB, uint8 mX, uint8 mD){
  R=mR;
  G=mG;
  B=mB;
  X=mX;
  D=mD;
  computeYIQ();
}

Features:: Features(uint32 p, uint8 d){
  
  X = (uint8)((p & 0xFF000000)>> 24);
  B = (uint8)((p & 0x00FF0000)>> 16);
  G = (uint8)((p & 0x0000FF00)>> 8);
  R = (uint8)((p & 0x000000FF));
  D = d;
  computeYIQ();
}

Features::Features(Features & f){
  R=f.getR();
  G=f.getG();
  B=f.getB();
  X=f.getX();
  D=f.getD();
  LUM=f.getLUM();
  I=f.getI();
  Q=f.getQ();
}

uint8 Features::getNBChamps (FeatureType Type)
{
  switch( Type ) {
  case RGB:
    return 3;
    break;
  case RGBLUM:
    return 4;
    break;
  case IQ:
    return 2;
    break;
  case LUMIN:
    return 1;
    break;
  case DEPRGB:
    return 4;
    break;
  case DEPLUM:
    return 2;
    break;
  default:
    return 0;
  }
}

uint8  Features::getR (){
  return R;  
}

uint8 Features::getG (){
  return G; 
}

uint8 Features::getB (){
  return B; 
}

uint8 Features::getX (){
  return X; 
}

uint8 Features::getD (){
  return D; 
}

uint8 Features::getLUM (){
  return LUM; 
}

uint8 Features::getI(){
  return I;
}

uint8 Features::getQ(){
  return Q;
}

void Features::setR (uint8 mr){
  R=mr; 
}

void Features::setG (uint8 mg){
  G=mg; 
}

void Features::setB (uint8 mb){
  B=mb; 
}

void Features::setX (uint8 mx){
  X=mx; 
}

void Features::setD (uint8 md){
  D=md; 
}

void Features::setLUM (uint8 ml){
  LUM=ml;
}

void Features::setI(uint8 mi){
  I=mi;
}

void Features::setQ(uint8 mq){
  Q=mq;
}

uint32 Features::toPixel()
{  
  return (((uint32 )X << 24) |((uint32) B << 16)|((uint32)G << 8) |(uint32) R) ;
}

void Features::computeYIQ()
{
  LUM=0.299*R+0.587*G+0.114*B;
  I=0.596*R-0.274*G-0.322*B;
  Q=0.212*R-0.523*G+0.311*B;
}

void Features::add(vector<uint8> * v, FeatureType Type)
{
  switch( Type ) {
      case RGB:
	(*v).push_back(R);
	(*v).push_back(G);
	(*v).push_back(B);
	break;
      case RGBLUM:
	(*v).push_back(R);
	(*v).push_back(G);
	(*v).push_back(B);
	(*v).push_back(LUM);
	break;
      case IQ:
	(*v).push_back(I);
	(*v).push_back(Q);
	break;
      case LUMIN:
	(*v).push_back(LUM);
	break;
      case DEPRGB:
	(*v).push_back(R);
	(*v).push_back(G);
	(*v).push_back(B);
	(*v).push_back(D);
	break;
      case DEPLUM:
	(*v).push_back(LUM);
	(*v).push_back(D);
	break;
   }
}

void Features::add(ANNpoint p,int j, FeatureType Type)
{
  switch( Type ) {
      case RGB:
	p[j]=R;
	p[j+1]=G;
	p[j+2]=B;
	break;
      case RGBLUM:
	p[j]=R;
	p[j+1]=G;
	p[j+2]=B;
	p[j+3]=LUM;
	break;
      case IQ:
	p[j]=I;
	p[j+1]=Q;
	break;
      case LUMIN:
	p[j]=LUM;
	break;
      case DEPRGB:
	p[j]=R;
	p[j+1]=G;
	p[j+2]=B;
	p[j+3]=D;
	break;
      case DEPLUM:
	p[j]=LUM;
	p[j+1]=D;
	break;
   }
}

void Features::computeRGB(Features* colorFeature) 
{
  uint8 Rb = colorFeature->getR();
  uint8 Gb = colorFeature->getG();
  uint8 Bb = colorFeature->getB();
  this->computeRGB(Rb, Gb, Bb, colorFeature->getX());
}

void Features::computeRGB(uint8 Rb, uint8 Gb, uint8 Bb, uint8 Xb) 
{
  // changement de base pour l'image b
  float Ib = 0.596*Rb - 0.275*Gb - 0.321*Bb;
  float Qb = 0.212*Rb - 0.523*Gb + 0.311*Bb;
  // changement des RGB de l'image b'
  int tmp = (1.702*Qb + this->LUM - 1.105*Ib);
  if (tmp > 255) {tmp = 255;}
  if (tmp < 0) {tmp = 0;}
  this->B = (uint8) tmp;
  tmp = (-0.647*Qb + this->LUM - 0.272*Ib);
  if (tmp > 255) {tmp = 255;}
  if (tmp < 0) {tmp = 0;}
  this->G = (uint8) tmp;
  tmp = (0.621*Qb + this->LUM + 0.956*Ib);
  if (tmp > 255) {tmp = 255;}
  if (tmp < 0) {tmp = 0;}
  this->R = (uint8) tmp;
  this->X = Xb;
}
