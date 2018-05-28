#include "common.h"
#include "texture.h"
#include <time.h>

int vrstartx, vrfinishx, vrstarty, vrfinishy;
int *candlistx, *candlisty;
extern params data;
extern int *atlas;
int anotherpass=0,maxcand = 40;
pixelvalue *target;
int *xloopout, *yloopout;
int *xloopin, *yloopin;

double compare_rest(pixelvalue *image,int x, int y,
                      pixelvalue *tar,int x1, int y1);

///////////////////////////////////////////////////////////////////
// This is the main texture synthesis function. Called just once
// from main to generate image from 'image' into 'result'.
// Synthesis parameters (image and neighborhood sizes) are in global
// 'data' structure.
///////////////////////////////////////////////////////////////////

void create_texture(pixelvalue *image, pixelvalue *result, params &dumb){
 int i,j,k, ncand, bestx,besty;
 double diff,curdiff;
 int tsx,tsy;
 //read_ppm("14.t.ppm",tsx,tsy,target);
 //if(tsx != data.widthout || tsx != data.heightout){
 // printf("Target image is of incorrect size. Exiting\n");
 // exit(1);
 //}
 // Atlas is to hold information about where in the input image did a
 // pixel came from. In principle, only need a small atlas, but
 // updating it would complicate the code, so make it of the size
 // of the output image.

 //atlas = (int *)malloc(2*data.widthout*data.heightout*sizeof(int));
 //atlasy = (int *)malloc(data.widthout*data.heightout*sizeof(int));

 // This is to hold a list of candidates - very small arrays
 // (25 or so integers each)
 candlistx = (int*)malloc((data.localx*data.localy+1)*sizeof(int));
 candlisty = (int*)malloc((data.localx*data.localy+1)*sizeof(int));

 if(!anotherpass) init(result, image);

 printf("Finished initialization\n");

///////////////////////////////////////////////////////////////////////
// This is the start of the main texture synthesis loop. If there is
// anything worth optimizing, it's here (and, obviously, in functions
// it calls): initialization ('init' function)
// and edge handling (a loop below this one) consume currently just a
// fraction of total time.
//////////////////////////////////////////////////////////////////////

 for(i=0;i<data.heightout-data.localy/2;i++){
  for(j=0;j<data.widthout;j++){
   // First, create a list of candidates for particular pixel.
   if(anotherpass) ncand = create_all_candidates(j,i);
   else  ncand = create_candidates(j,i);

   // If there are multiple candidates, choose the best based on L_2
   // norm

   if(ncand > 1){
   diff = 1e10;
   for(k=0;k<ncand;k++){
    curdiff = compare_neighb(image,candlistx[k],candlisty[k],result,j,i);
    curdiff += compare_rest(image,candlistx[k],candlisty[k],target,j,i);
    if(curdiff < diff){
     diff = curdiff;
     bestx = candlistx[k];
     besty = candlisty[k];
    }
   }
   }else{
    bestx = candlistx[0];
    besty = candlisty[0];
   }

   // Copy the best candidate to the output image and record its position
   // in the atlas (atlas is used to create candidates)

//bestx=besty = 0;
   result[a(j,i,data.widthout)+R] = image[a(bestx,besty,data.widthin)+R];
   result[a(j,i,data.widthout)+G] = image[a(bestx,besty,data.widthin)+G];
   result[a(j,i,data.widthout)+B] = image[a(bestx,besty,data.widthin)+B];
   atlas[aa(j,i)] = bestx;
   atlas[aa(j,i)+1] = besty;
  }
 }

 // Use full neighborhoods for the last few rows. This is a small
 // fraction of total area - can be ignored for optimization purposes.

 for(;i<data.heightout;i++){
  for(j=0;j<data.widthout;j++){
   ncand = create_all_candidates(j,i);
   //ncand = create_candidates(j,i);
   if(ncand > 1){
    diff = 1e10;
    for(k=0;k<ncand;k++){
     curdiff = compare_full_neighb(image,candlistx[k],candlisty[k],result,j,i);
     if(curdiff < diff){
      diff = curdiff;
      bestx = candlistx[k];
      besty = candlisty[k];
     }
    } // for k
   }else{
    bestx = candlistx[0];
    besty = candlisty[0];
   }
   result[a(j,i,data.widthout)+R] = image[a(bestx,besty,data.widthin)+R];
   result[a(j,i,data.widthout)+G] = image[a(bestx,besty,data.widthin)+G];
   result[a(j,i,data.widthout)+B] = image[a(bestx,besty,data.widthin)+B];
   atlas[aa(j,i)] = bestx;
   atlas[aa(j,i)+1] = besty;
  }
 }
 printf("Finished first pass\n");

 ///////////////////////////////////////////////////////////////////////
 // End of main texture synthesis loop
 ///////////////////////////////////////////////////////////////////////

 // Redo first couple of rows for better vertical tilability
 // This takes a small fraction of total time, so it can be ignored for now

 //for(i=0;i<data.heightout;i++){
 for(i=0;i<data.localy/2;i++){
  for(j=0;j<data.widthout;j++){
   ncand = create_all_candidates(j,i);
   if(ncand > 1){
    diff = 1e10;
    for(k=0;k<ncand;k++){
     curdiff = compare_full_neighb(image,candlistx[k],candlisty[k],result,j,i);
     if(curdiff < diff){
      diff = curdiff;
      bestx = candlistx[k];
      besty = candlisty[k];
     }
    }
   }else{
    bestx = candlistx[0];
    besty = candlisty[0];
   }
   result[a(j,i,data.widthout)+R] = image[a(bestx,besty,data.widthin)+R];
   result[a(j,i,data.widthout)+G] = image[a(bestx,besty,data.widthin)+G];
   result[a(j,i,data.widthout)+B] = image[a(bestx,besty,data.widthin)+B];
   atlas[aa(j,i)] = bestx;
   atlas[aa(j,i)+1] = besty;
  }
 }
 //postprocess(result);
 printf("Finished frame\n");
 //write_ppm(data.widthout, data.heightout, "zad.ppm",result);
 //printf("Finished wrtiting frame\n");
 //free(atlas);
}

// Marks boundaries of texture pieces. Not necessary - ignore it now.

void postprocess(pixelvalue *result){
 int i,j,sx,sy;
 for(i=1;i<data.heightout-1;i++){
  for(j=1;j<data.widthout-1;j++){
   sx = atlas[aa(j,i)]; sy = atlas[aa(j,i)+1];
   if(!(sx == atlas[aa(j+1,i)]-1 && sx == atlas[aa(j-1,i)]+1 &&
      sy == atlas[aa(j,i+1)+1]-1 && sy == atlas[aa(j,i-1)+1]+1) ){
    result[a(j,i,data.widthout)+R] = 1.0;
    result[a(j,i,data.widthout)+G] = 1.0;
    result[a(j,i,data.widthout)+B] = 1.0;
    /*
    result[a(j,i,data.widthout)+R] = (result[a(j-1,i,data.widthout)+R]+
     result[a(j+1,i,data.widthout)+R]+result[a(j,i-1,data.widthout)+R]+
     result[a(j,i+1,data.widthout)+R]+result[a(j-1,i-1,data.widthout)+R]+
     result[a(j-1,i+1,data.widthout)+R]+result[a(j+1,i-1,data.widthout)+R]+
     result[a(j+1,i+1,data.widthout)+R])/8;
    result[a(j,i,data.widthout)+G] = (result[a(j-1,i,data.widthout)+G]+
     result[a(j+1,i,data.widthout)+G]+result[a(j,i-1,data.widthout)+G]+
     result[a(j,i+1,data.widthout)+G]+result[a(j-1,i-1,data.widthout)+G]+
     result[a(j-1,i+1,data.widthout)+G]+result[a(j+1,i-1,data.widthout)+G]+
     result[a(j+1,i+1,data.widthout)+G])/8;
    result[a(j,i,data.widthout)+B] = (result[a(j-1,i,data.widthout)+B]+
     result[a(j+1,i,data.widthout)+B]+result[a(j,i-1,data.widthout)+B]+
     result[a(j,i+1,data.widthout)+B]+result[a(j-1,i-1,data.widthout)+B]+
     result[a(j-1,i+1,data.widthout)+B]+result[a(j+1,i-1,data.widthout)+B]+
     result[a(j+1,i+1,data.widthout)+B])/8;
    */
   }
  }
 }
}

// Creates a list of valid candidates for given pixel using only L-shaped
// causal area

int create_candidates(int x,int y){
 int address,i,j,k,n = 0;
 for(i=0;i<=data.localy/2;i++){
  for(j=-data.localx/2;j<=data.localx/2;j++){
   if(i==0 && j>=0) continue;
   //address = aa((data.widthout+x+j)%data.widthout,
   //                      (data.heightout+y-i)%data.heightout);
   address = aa(xloopout[x+j],yloopout[y-i]);
   candlistx[n] = atlas[address]-j;
   candlisty[n] = atlas[address+1]+i;

   if(candlistx[n] >= vrfinishx || candlistx[n] < vrstartx) {
     candlistx[n] = vrstartx+(int)(drand48()*(vrfinishx-vrstartx));
     candlisty[n] = vrstarty+(int)(drand48()*(vrfinishy-vrstarty));
     n++;
     continue;
   }

   //double tmp = drand48(); tmp *= tmp; tmp *=tmp;
   if(candlisty[n] >= vrfinishy ){
     // || tmp > ((double)vrfinishy-candlisty[n])/(vrfinishy-vrstarty) ){ 
     candlisty[n] = vrstarty+(int)(drand48()*(vrfinishy-vrstarty));
     candlistx[n] = vrstartx+(int)(drand48()*(vrfinishx-vrstartx));
     n++;
     continue;
   }
   for(k=0;k<n;k++){
    if(candlistx[n] == candlistx[k] && candlisty[n] == candlisty[k]){ n--; break;}
   }
   n++;
   //if(n == maxcand) return n;
  }
 }
/*
 for(i=-data.localy/2;i<0;i++){
  for(j=-data.localx/2;j<=data.localx/2;j++){
   address = a((data.widthout+x+j)%data.widthout,
                         (data.heightout+y-i)%data.heightout,data.widthout);
   if(target[address+B] != 1.0){
    address = aa((data.widthout+x+j)%data.widthout,
                         (data.heightout+y-i)%data.heightout);
    candlistx[n] = atlas[address];
    candlisty[n] = atlas[address+1];
    n++;
    return n;
   }
  }
 }
*/
 return n;
}

// Created a list of candidates using the complete square around the pixel

int create_all_candidates(int x,int y){
 int address,i,j,k,n = 0;
 for(i=-data.localy/2;i<=data.localy/2;i++){
  for(j=-data.localx/2;j<=data.localx/2;j++){
   if(i==0 && j>=0) continue;
   //address = aa((data.widthout+x+j)%data.widthout,
   //                      (data.heightout+y-i)%data.heightout);
   address = aa(xloopout[x+j],yloopout[y-i]);
   candlistx[n] = atlas[address]-j;
   candlisty[n] = atlas[address+1]+i;

   if(candlistx[n] >= vrfinishx || candlistx[n] < vrstartx) {
     candlistx[n] = vrstartx+(int)(drand48()*(vrfinishx-vrstartx));
     candlisty[n] = vrstarty+(int)(drand48()*(vrfinishy-vrstarty));
     n++;
     continue;
   }

   //tmp = drand48(); tmp *= tmp; tmp *=tmp;
   if(candlisty[n] >= vrfinishy || candlisty[n] < vrstarty){
     // || tmp > ((double)vrfinishy-candlisty[n])/(vrfinishy-vrstarty) ){ 
     candlisty[n] = vrstarty+(int)(drand48()*(vrfinishy-vrstarty));
     candlistx[n] = vrstartx+(int)(drand48()*(vrfinishx-vrstartx));
     n++;
     continue;
   }
   for(k=0;k<n;k++){
    if(candlistx[n] == candlistx[k] && candlisty[n] == candlisty[k]){ n--; break;}
   }
   n++;
   //if(n == maxcand) return n;
  }
 }
 return n;
}

// Initializes the output image and atlases to a random collection of pixels

void init(pixelvalue *result, pixelvalue *image){
 int i,j,tmpx,tmpy;
 vrstartx = data.localx/2; vrstarty = data.localy/2;
 vrfinishx = data.widthin-data.localx/2;
 vrfinishy = data.heightin-data.localy/2;
 for(i=0;i<data.heightout;i++){
  for(j=0;j<data.widthout;j++){
   if(target[a(j,i,data.widthout)+R] == 1.0 &&
       target[a(j,i,data.widthout)+G] == 1.0 &&
       target[a(j,i,data.widthout)+B] == 1.0){
    tmpx = vrstartx+(int)(drand48()*(vrfinishx-vrstartx));
    tmpy = vrstarty+(int)(drand48()*(vrfinishy-vrstarty));
    if(!anotherpass){
    atlas[aa(j,i)] = tmpx; atlas[aa(j,i)+1] = tmpy;
    result[a(j,i,data.widthout)+R] = image[a(tmpx,tmpy,data.widthin)+R];
    result[a(j,i,data.widthout)+G] = image[a(tmpx,tmpy,data.widthin)+G];;
    result[a(j,i,data.widthout)+B] = image[a(tmpx,tmpy,data.widthin)+B];;
    }
   }
  }
 }
 return;
}


// Compares two square neighborhoods, returns L_2 difference

double compare_full_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1){
 //printf("Comparing %d %d with %d %d\n",x,y,x1,y1);
 double tmp,res = 0;
 int i,j,addr,addr1;
 for(i=-(data.localy/2);i<=data.localy/2;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    if( !( i > 0 && y1 > data.localy && y1+i < data.heightout) ){
    addr = a(x+j,y+i,data.widthin);
    //addr1 = a((data.widthout+x1+j)%data.widthout,
    //              (data.heightout+y1+i)%data.heightout,data.widthout);
    addr1 = a(xloopout[x1+j],yloopout[y1+i],data.widthout);

    tmp = image[addr+R] - image1[addr1+R];
    res += tmp*tmp;
    tmp = image[addr+G] - image1[addr1+G];
    res += tmp*tmp;
    tmp = image[addr+B] - image1[addr1+B];
    res += tmp*tmp;
    }
   }
 }
 return res;
}

// Compares two L-shaped neighborhoods, returns L_2 difference

double compare_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1){
 double tmp,res = 0;
 int i,j,addr,addr1;
 for(i=-(data.localy/2);i<0;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    addr = a(x+j,y+i,data.widthin);
    //addr1 = a((data.widthout+x1+j)%data.widthout,
    //              (data.heightout+y1+i)%data.heightout,data.widthout);
    addr1 = a(xloopout[x1+j],yloopout[y1+i],data.widthout);

    tmp = image[addr+R] - image1[addr1+R];
    res += tmp*tmp;
    tmp = image[addr+G] - image1[addr1+G];
    res += tmp*tmp;
    tmp = image[addr+B] - image1[addr1+B];
    res += tmp*tmp;
   }
 }
 
//return res;
 for(j=-(data.localx/2);j<0;j++){
 //for(j=-(data.localx/2);j<=data.localx/2;j++){

   addr = a(x+j,y,data.widthin);
   //addr1 = a((data.widthout+x1+j)%data.widthout, y1,data.widthout);
   addr1 = a(xloopout[x1+j],y1,data.widthout);

   tmp = image[addr+R] - image1[addr1+R];
   res += tmp*tmp;
   tmp = image[addr+G] - image1[addr1+G];
   res += tmp*tmp;
   tmp = image[addr+B] - image1[addr1+B];
   res += tmp*tmp;

  }

 return res;
}

double compare_rest(pixelvalue *image,int x, int y,
                      pixelvalue *tar,int x1, int y1){
 double tmp,res = 0;
 int i,j,addr,addr1;
 for(i=(data.localy/2);i>0;i--){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    addr = a(x+j,y+i,data.widthin);
    //addr1 = a((data.widthout+x1+j)%data.widthout,
    //              (data.heightout+y1+i)%data.heightout,data.widthout);
    addr1 = a(xloopout[x1+j],yloopout[y1+i],data.widthout);
    if(tar[addr1+B] != 1.0){
     tmp = image[addr+R] - tar[addr1+R];
     res += tmp*tmp;
     tmp = image[addr+G] - tar[addr1+G];
     res += tmp*tmp;
     tmp = image[addr+B] - tar[addr1+B];
     res += tmp*tmp;
    }
   }
 }
 
//return res;
 for(j=(data.localx/2);j>0;j--){
 //for(j=-(data.localx/2);j<=data.localx/2;j++){

   addr = a(x+j,y,data.widthin);
   //addr1 = a((data.widthout+x1+j)%data.widthout, y1,data.widthout);
   addr1 = a(xloopout[x1+j],y1,data.widthout);
   if(tar[addr1+B] != 1.0){
    tmp = image[addr+R] - tar[addr1+R];
    res += tmp*tmp;
    tmp = image[addr+G] - tar[addr1+G];
    res += tmp*tmp;
    tmp = image[addr+B] - tar[addr1+B];
    res += tmp*tmp;
   }
  }

 return res;
}


////////////////////////////////////////////////////////////////////
// IGNORE THE REST !!! THESE FUNCTIONS ARE HERE DUE TO
// MY EARLIER EXPERIMENTATION BUT NONE ARE CURRENTLY USED
// I wanted to keep them around in case I need to come back to some
// of the ideas I've tried before. Plus I'm not sure if removing some of them 
// won't cause compiler to complain.
///////////////////////////////////////////////////////////////////



int comparesign(const void *s1, const void *s2){
 if( ((signature *)s1)->sign < ((signature *)s2)->sign) return -1;
 if( ((signature *)s1)->sign > ((signature *)s2)->sign) return 1;
 return 0;
}
double findbestmatch(pixelvalue *image, int candidatex, int candidatey,
                   pixelvalue *result, int currentx, int currenty,
                   int &bestx, int &besty, params &data){

//void findbestmatch(pixelvalue *image,double signa,
//                  int &bestx, int &besty, params &data){
/*
int i,j;
double curdiff, diff = 1e10;
for(i=data.localy;i<data.heightin;i++){
  for(j=data.localx/2;j<data.widthin-data.localx/2;j++){
   curdiff = compare_neighb(image, j,i,result,currentx,currenty, data);
   if(curdiff < diff){
     diff = curdiff;
     bestx = j; besty = i;
    }
  }
}
*/
int addr = (candidatey-data.localy/2)*(data.widthin-data.localx+1)+
       (candidatex-data.localx/2);
int secondchoicex = 
  ((signature *)SIGNATURES)[addr].secondx;
int secondchoicey = 
  ((signature *)SIGNATURES)[addr].secondy;
double dist1 =
 compare_neighb(image,candidatex, candidatey, result,currentx, currenty, data);
double dist2 =
 compare_neighb(image,secondchoicex, secondchoicey, result,currentx, currenty, data);
//if(dist1 < dist2){
 bestx = candidatex; besty = candidatey;return dist1;
//}else{
// bestx = secondchoicex; besty = secondchoicey;return dist2;
//}

/*
int index;
int index_above = (data.heightin-data.localy/2)*(data.widthin-data.localx+1)-1;
int index_below = 0;
 double tmp;
 while(index_above - index_below > 1){
  index = (index_above+index_below)/2;
  if(signa > ((signature *)SIGNATURES)[index].sign) index_below = index;
  else index_above = index;
 }
 if( fabs(((signature *)SIGNATURES)[index_above].sign - signa) <
     fabs(((signature *)SIGNATURES)[index_below].sign - signa) ){
  bestx = ((signature *)SIGNATURES)[index_above].x;
  besty = ((signature *)SIGNATURES)[index_above].y;
  tmp = ((signature *)SIGNATURES)[index_above].sign;
 }else{
  bestx = ((signature *)SIGNATURES)[index_below].x;
  besty = ((signature *)SIGNATURES)[index_below].y;
  tmp = ((signature *)SIGNATURES)[index_below].sign;
 }
*/
 //printf("Best match to %f is %f with index_b = %d, index_a = %d\n",
  //      signa, tmp, index_below,index_above );
/*
 double curdiff, diff = 10e10;
 for(i=0;i<data.heightin;i++){
  for(j=0;j<data.widthin;j++){
   curdiff = fabs(char_neighb(image,j,i, data, IN) - signature);
   if(curdiff < diff){
    diff = curdiff;
    bestx = j; besty = i;
   }
  }
 }
*/
}

double char_neighb(pixelvalue *image,int x,int y, params &data, int flag){
 // Assume the neighbourhood is X by Y where both X and Y are odd.
 double res = 0;
 int i,j;
 pixelvalue tmp;
 if(flag == IN){
  for(i=-(data.localy/2);i<0;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    tmp = image[a((data.widthin+x+j)%data.widthin,(data.heightin+y+i)%data.heightin,data.widthin)+R];
    res += tmp;
   }
  }
  // Last row (incomplete)
  for(j=-(data.localx/2);j<0;j++){
   tmp = image[a((data.widthin+x+j)%data.widthin,y,data.widthin)+R];
   res += tmp;
  }
  res /= data.localx*data.localy/2;
 }else if(flag == OUT){
  for(i=-(data.localy/2);i<0;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    tmp = image[a((data.widthout+x+j)%data.widthout,(data.heightout+y+i)%data.heightout,data.widthout)+R];
    res += tmp;
   }
  }
  // Last row (incomplete)
  for(j=-(data.localx/2);j<0;j++){
   tmp = image[a((data.widthout+x+j)%data.widthout,y,data.widthout)+R];
   res += tmp;
  }
  res /= data.localx*data.localy/2;
 }
 return res;
}

double compare_neighb(pixelvalue *image,int x, int y,
                      pixelvalue *image1,int x1, int y1, params &data){
 return compare_neighb(image,x, y,image1,x1,y1);
}
double compare_neighb(pixelvalue *image,int x, int y,
                      int x1, int y1, params &data){
 double tmp,res = 0;
 int i,j;
 for(i=-(data.localy/2);i<0;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    tmp = image[a((data.widthin+x+j)%data.widthin,
                  (data.heightin+y+i)%data.heightin,data.widthin)+R] -
          image[a((data.widthin+x1+j)%data.widthin,
                  (data.heightin+y1+i)%data.heightin,data.widthin)+R];
    res += tmp*tmp;
   }
 }
 for(j=-(data.localx/2);j<0;j++){
   tmp = image[a((data.widthin+x+j)%data.widthin,y,data.widthin)+R] -
         image[a((data.widthin+x1+j)%data.widthin,y1,data.widthin)+R];
   res += tmp*tmp;
  }
 return res;
}
double findbestlocal(pixelvalue *image,int x,int y, int x1, int y1,
                    int &bestx, int &besty, params &data){
 int i,j;
 double curdiff,diff = 1e10;
 for(i=-(data.localy/2);i<=data.localy/2;i++){
   for(j=-(data.localx/2);j<=data.localx/2;j++){
    curdiff = compare_neighb(image,x,y,
        x1+j,y1+i,data);
    if(curdiff < diff){
     diff = curdiff;
     bestx = x1+j; besty = y1+i;
    }
   }
 }
 return diff;
}
