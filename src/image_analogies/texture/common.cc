#include "common.h"

/* Reads in a binary .ppm file. Allocates an array of necessary size
   for the image. Returns image size in X and Y
*/
void read_ppm(char *fn, int &X, int &Y, pixelvalue * &result){
 int i,j,tmpint;
 FILE *in_file = fopen(fn,"r");
 if(in_file == NULL){
  printf("Can't open input file %s. Exiting.\n",fn);
  exit(1);
 }
 fscanf(in_file,"P%d\n", &tmpint);
 if(tmpint != 6){
  printf("Input file is not binary ppm. Exiting.\n");
  exit(1);
 }
 char ch;
 fscanf(in_file, "%c", &ch);
 while(ch == '#'){
  while(ch != '\n') fscanf(in_file,"%c", &ch); fscanf(in_file, "%c", &ch);
 }
 fseek (in_file, -1, 1);
 fscanf(in_file,"%d %d\n%d\n", &X, &Y, &tmpint);
 if(tmpint != 255) printf("Warning: maxvalue is not 255 in ppm file\n");
 result = (pixelvalue *)malloc(3*X*Y*sizeof(pixelvalue));
 if(result == NULL){
  printf("Can't allocate image of size %dx%d. Exiting\n",X,Y);
  exit(1);
 }else{
  printf("Reading image %s of size %dx%d\n",fn,X,Y);
 }
 unsigned char tmp[3];
 pixelvalue r,g,b;
 for(i=0;i<Y;i++){
  for(j=0;j<X;j++){
   fscanf(in_file,"%c",&tmp[R]);
   fscanf(in_file,"%c",&tmp[G]);
   fscanf(in_file,"%c",&tmp[B]);
   convert_pixel_from_rgb(&r,&g,&b,tmp);
   result[a(j,Y-1-i,X)+R] = r/tmpint;
   result[a(j,Y-1-i,X)+G] = g/tmpint;
   result[a(j,Y-1-i,X)+B] = b/tmpint;
  }
 }
 fclose(in_file);
 printf("Done reading image.\n");
}

/* Writes output image of size XxY to a raw ppm file with name out_fn.
   Assumes row-major order for result. result is 3*X*Y array of
   pixel values to be written out. It  calls function
   convert_pixel on each value befor writing it to a file.
*/
void write_ppm(int X, int Y, char *out_fn, pixelvalue *result){
 printf("Inside write out, %s: %d %d\n",out_fn,X,Y);
 FILE *out_file;
 int i,j;
 unsigned char tmp[3];
 pixelvalue r,g,b;
 out_file = fopen(out_fn,"w");
 //printf("File opened\n");
 if(out_file == NULL){
  printf("Can't open output file. Exiting.\n");
  exit(1);
 }
 //printf("P6\n%d %d\n255\n",X,Y);
 fprintf(out_file,"P6\n%d %d\n255\n",X,Y);
 printf("Starting write out\n");
 for(i=0;i<Y;i++){
  for(j=0;j<X;j++){
   r = result[a(j,Y-1-i,X)+R]*255;
   g = result[a(j,Y-1-i,X)+G]*255;
   b = result[a(j,Y-1-i,X)+B]*255;
   convert_pixel(r,g,b,tmp);
   fprintf(out_file,"%c",tmp[R]);
   fprintf(out_file,"%c",tmp[G]);
   fprintf(out_file,"%c",tmp[B]);
  }
  //printf("Finished line %d\n",i);
 }
 fclose(out_file);
}

