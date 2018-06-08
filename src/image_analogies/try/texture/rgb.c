#include "rgb.h"

/* Performs a conversion from pixelvalue type to unsigned char.
   Code for color space conversions should go here, in case pixelvalues
   represent something other that simple RGB (YCrCb, for example).
   Copy values to out[] otherwise.
   out should be at least 3*sizeof(unsigned char).
*/

void convert_pixel(pixelvalue red, pixelvalue green, pixelvalue blue,
                   unsigned char *out){

   /* Simply copy values to out[], performing boundary check */
   
   red = ( red > 0 ) ? red : 0;                  
   out[R] = (unsigned char) (( red < 255 ) ? red : 255);
   green = ( green > 0 ) ? green : 0;
   out[G] = (unsigned char) (( green < 255 ) ? green : 255);
   blue = ( blue > 0 ) ? blue : 0;
   out[B] = (unsigned char) (( blue < 255 ) ? blue : 255);

   /* Commented out is YCrCb --> RGB conversion

   tmp = 1.164*red+1.596*green;
   tmp = ( tmp > 0 ) ? tmp : 0;
   out[R] = ( tmp < 255 ) ? tmp : 255;
   tmp = 1.164*red-0.813*green-0.392*blue;
   tmp = ( tmp > 0 ) ? tmp : 0;
   out[G] = ( tmp < 255 ) ? tmp : 255;
   tmp = 1.164*red+2.017*blue;
   tmp = ( tmp > 0 ) ? tmp : 0;
   out[B] = ( tmp < 255 ) ? tmp : 255;
   */

}

void convert_pixel_from_rgb(pixelvalue *red, pixelvalue *green, pixelvalue *blue,
                   unsigned char *in){

   /* Simply copy values from in[] */
   *red = (pixelvalue)(in[R]);
   *green = (pixelvalue)(in[G]);
   *blue = (pixelvalue)(in[B]);

}


void write_rgb(int X, int Y, char *out_fn, pixelvalue *result){
 //FILE *out_file;
 pixelvalue r,g,b;
 int i,j;
 unsigned char tmp[3];
 unsigned short *bufferr,*bufferg,*bufferb;
 IMAGE *image = iopen(out_fn,"w",RLE(1),3,X,Y,3);
 if(image == NULL){
  printf("Can't open output rgb file\n");
  exit(1);
 }
 bufferr = (unsigned short *)malloc(X*sizeof(unsigned short));
 bufferg = (unsigned short *)malloc(X*sizeof(unsigned short));
 bufferb = (unsigned short *)malloc(X*sizeof(unsigned short));
 for(i=0;i<Y;i++){
  for(j=0;j<X;j++){
   r = result[a(j,i,X)+R]*255;
   g = result[a(j,i,X)+G]*255;
   b = result[a(j,i,X)+B]*255;
   convert_pixel(r,g,b,tmp);
   bufferr[j] = tmp[R];
   bufferg[j] = tmp[G];
   bufferb[j] = tmp[B];
  }
  putrow(image,bufferr,Y-i-1,0);
  putrow(image,bufferg,Y-i-1,1);
  putrow(image,bufferb,Y-i-1,2);
 }
 free(bufferr);
 free(bufferg);
 free(bufferb);
 iclose(image);
}

void read_rgb(char *in_fn, int *Xout, int *Yout,pixelvalue **resultout){
 //FILE *out_file;
 pixelvalue r,g,b;
 pixelvalue *result;
 int x,y,X,Y;
 unsigned char tmp[3];
 unsigned short *bufferr,*bufferg,*bufferb;
 IMAGE *image = iopen(in_fn,"r");
 if(image == NULL){
  printf("Can't open input rgb file\n");
  exit(1);
 }
 X = image->xsize; *Xout = X;
 Y = image->ysize; *Yout = Y;
 result = (pixelvalue *)malloc(3*X*Y*sizeof(pixelvalue));
 if(result == NULL){
  printf("Can't allocate image of size %dx%d. Exiting\n",X,Y);
  exit(1);
 }else{
  printf("Reading image %s of size %dx%d\n",in_fn,X,Y);
 }
 bufferr = (unsigned short *)malloc(X*sizeof(unsigned short));
 bufferg = (unsigned short *)malloc(X*sizeof(unsigned short));
 bufferb = (unsigned short *)malloc(X*sizeof(unsigned short));
 if(image->zsize == 1) {
   printf("Warning: This is a black and write image\n");
   for(y=0; y<image->ysize; y++) {
    getrow(image,bufferr,y,0);
    for(x=0; x<image->xsize; x++){
       tmp[R] = tmp[G] = tmp[B] = bufferr[x];
       convert_pixel_from_rgb(&r,&g,&b,tmp);
       result[a(x,image->ysize-y-1,X)+R] = r/255;
       result[a(x,image->ysize-y-1,X)+G] = g/255;
       result[a(x,image->ysize-y-1,X)+B] = b/255;
    }
   }
  }else if(image->zsize >= 3) {  /* if the image has alpha zsize is 4 */
   printf("This is a rgb image\n");
   for(y=0; y<image->ysize; y++) {
    getrow(image,bufferr,y,0);
    getrow(image,bufferg,y,1);
    getrow(image,bufferb,y,2);
    for(x=0; x<image->xsize; x++){
     tmp[R] = bufferr[x]; tmp[G] = bufferg[x]; tmp[B] = bufferb[x];
     convert_pixel_from_rgb(&r,&g,&b,tmp);
     result[a(x,image->ysize-y-1,X)+R] = r/255;
     result[a(x,image->ysize-y-1,X)+G] = g/255;
     result[a(x,image->ysize-y-1,X)+B] = b/255;
     //if(x == 0)
     //printf("Pixel %d %d: %d %d %d\n",x,y,tmp[R],tmp[G],tmp[B]);
    }
   }
  }
 free(bufferr);
 free(bufferg);
 free(bufferb);
 iclose(image);
 *resultout = result;
 printf("Done reading\n");
}
