#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rgb.h"

typedef struct{
 int localx, localy, localz;
 int widthin, widthout;
 int heightin, heightout; 
 int nfin, nfout;
}params;
#define OUT 0
#define IN 1

// Image handling definitions and functions
extern "C" void write_rgb(int, int, char *, pixelvalue *);
extern "C" void read_rgb(char *, int*, int *,pixelvalue **);
extern "C" void convert_pixel(pixelvalue red, pixelvalue green, pixelvalue blue,
                   unsigned char *out);
extern "C" void convert_pixel_from_rgb(pixelvalue *red, pixelvalue *green, pixelvalue *blue,
                   unsigned char *in);
void read_ppm(char *fn, int &X, int &Y, pixelvalue * &result);
void write_ppm(int X, int Y, char *out_fn, pixelvalue *result);

void create_texture(pixelvalue *image, pixelvalue *result, params &data);

void parse_flags(int argc, char ** argv);
void init_params(params &data);
#endif
