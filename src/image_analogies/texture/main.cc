#include "common.h"
#include "ui.h"

pixelvalue *image;
int WIDTHin,HEIGHTin;
pixelvalue *result;
extern pixelvalue *target;
int WIDTH,HEIGHT;
int localx, localy,targetin=0;
params data;
int *atlas;
extern int *xloopout, *yloopout;
extern int *xloopin, *yloopin;

int  main(int argc, char **argv){
 if(argc < 7){
  printf("Usage: texture in_filename out_filename width height x y [target_name]\n");
  printf(" where width and height are the output image size;\n");
  printf("       x and y are the size of locally similar neighbourhood;\n");
  printf("       target_name is optional name of the target image\n");
  printf("       All images expected to be in raw ppm format\n");
  printf(" WARNING: RESEARCH CODE, NO EXTENSIVE TESTING, YOU ARE ON YOUR OWN !!!\n");
  exit(1);
 }
 parse_flags(argc, argv);
 init_params(data);
 printf("Ready to init UI\n");
 create_texture(image, result, data);

 glutInit(&argc, argv);
 modes[RGBA] = 1;
 modes[DOUBLEBUFFER] = 0;
 modes[DEPTH] = 0;
 setInitDisplayMode();
 makeMenus();
 makeWindow(0);
 makeWindow(1);
 makeWindow(2);
 glutIdleFunc(idleFunc);
 glutMenuStateFunc(menuStateFunc);
 glutMainLoop();
 return(0);

 write_ppm(data.widthout, data.heightout, argv[2],result);
 printf("Wrote frame\n");
 exit(0);
}

void parse_flags(int argc, char ** argv){
 int i, tsx,tsy;
 read_ppm(argv[1],WIDTHin,HEIGHTin,image);
 xloopin = (int *)malloc(2*WIDTHin*sizeof(int));
 yloopin = (int *)malloc(2*HEIGHTin*sizeof(int));
 for(i=-WIDTHin/2;i<WIDTHin+WIDTHin/2;i++){
  xloopin[i+WIDTHin/2] = (WIDTHin+i)%WIDTHin;
 }
 for(i=-HEIGHTin/2;i<HEIGHTin+HEIGHTin/2;i++){
  yloopin[i+HEIGHTin/2] = (HEIGHTin+i)%HEIGHTin;
 }
 xloopin += WIDTHin/2; yloopin += HEIGHTin/2;

 size[0][0] = WIDTHin; size[0][1] = HEIGHTin;
 pos[1][0] = pos[0][0]+WIDTHin+20, pos[1][1] = pos[0][1];
 WIDTH = atoi(argv[3]);
 HEIGHT = atoi(argv[4]);
 size[1][0] = WIDTH; size[1][1] = HEIGHT;
 pos[2][0] = pos[1][0]+WIDTH+20, pos[2][1] = pos[0][1];
 size[2][0] = WIDTH; size[2][1] = HEIGHT;

 result = (pixelvalue *)malloc(3*WIDTH*HEIGHT*sizeof(pixelvalue));
 if(argc>7){
  read_ppm(argv[7],tsx,tsy,target);
  if(tsx != WIDTH || tsy != HEIGHT){
   printf("Wrong target size. Exiting.\n");
   exit(1);
  }
  targetin = 1;
 }else target = (pixelvalue *)malloc(3*WIDTH*HEIGHT*sizeof(pixelvalue));
 atlas = (int *)malloc(2*WIDTH*HEIGHT*sizeof(int));
 xloopout = (int *)malloc(2*WIDTH*sizeof(int));
 yloopout = (int *)malloc(2*HEIGHT*sizeof(int));
 for(i=-WIDTH/2;i<WIDTH+WIDTH/2;i++){
  xloopout[i+WIDTH/2] = (WIDTH+i)%WIDTH;
 }
 for(i=-HEIGHT/2;i<HEIGHT+HEIGHT/2;i++){
  yloopout[i+HEIGHT/2] = (HEIGHT+i)%HEIGHT;
 }
 xloopout += WIDTH/2; yloopout += HEIGHT/2;

 if (result == NULL){
  printf("Can't allocate %dx%d image. Exiting.\n",WIDTH,HEIGHT);
  exit(1);
 }
 localx = atoi(argv[5]);
 localy = atoi(argv[6]);
}

void init_params(params &data){
 int i,j;
 data.localx = localx; data.localy = localy;
 data.widthin = WIDTHin; data.widthout = WIDTH;
 data.heightin = HEIGHTin; data.heightout = HEIGHT;
 if(!targetin){
  for(i=0;i<data.heightout;i++){
   for(j=0;j<data.widthout;j++){
    target[a(j,i,data.widthout)+R] = 1.0;
    target[a(j,i,data.widthout)+G] = 1.0;
    target[a(j,i,data.widthout)+B] = 1.0;
   }
  }
 }
 for(i=0;i<data.heightout;i++){
  for(j=0;j<data.widthout;j++){
   result[a(j,i,data.widthout)+R]  = 1.0;
   result[a(j,i,data.widthout)+G]  = 1.0;
   result[a(j,i,data.widthout)+B]  = 1.0;
  }
 }
}
