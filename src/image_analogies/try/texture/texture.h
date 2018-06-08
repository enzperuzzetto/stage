#include "common.h"

typedef struct{
 double sign, diff;
 int x,y;
 int secondx, secondy;
}signature;
void *SIGNATURES;

void init(pixelvalue *result, pixelvalue *image);
double findbestmatch(pixelvalue *image, int canx, int cany,
                   pixelvalue *result, int curx, int cury,
                   int &bestx, int &besty, params &data);
//void findbestmatch(pixelvalue *image,double signa,int &bestx, int &besty, params &data);
double char_neighb(pixelvalue *image,int x, int y, params &data, int flag);
int comparesign (const void *, const void *);
double findbestlocal(pixelvalue *image,int x,int y, int x1, int y1,
                    int &bestx, int &besty, params &data);
double compare_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1);
double compare_full_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1);
double compare_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1,params &data);
double compare_neighb(pixelvalue *image,int x, int y,
                      int x1, int y1, params &data);

int create_candidates(int x,int y);
int create_all_candidates(int x,int y);
void postprocess(pixelvalue *result);
