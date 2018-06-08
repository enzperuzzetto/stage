#ifndef RGB_H
#define RGB_H

#include <gl/image.h>

#define R 0
#define G 1
#define B 2
#define a(x,y,W)        (3*((y)*(W)+(x)))
#define aa(x,y) (2*((y)*data.widthout+(x)))
typedef float pixelvalue;

#endif
