#include "common.h"
#include "ui.h"
#include <time.h>
extern params data;
extern pixelvalue *target,*result,*image;
extern int *atlas, anotherpass;
int cursourcex,cursourcey,brushsizex,brushsizey;
int updateneeded, exactcolor, niter;
pixelvalue *brush;

void makeWindow(int index)
{
  char str[99];
  switch (index) {

  case 0:              /* ordinary RGB windows  */
    brushsizex = 11;//data.localx;
    brushsizey = 11; //data.localy;
    exactcolor = 0;
    niter = 1;
    brush = (pixelvalue *)malloc(3*brushsizex*brushsizey*sizeof(pixelvalue));
  case 1:
  case 2:
  case 3:

    setInitDisplayMode();
    glutInitWindowPosition(pos[index][0], pos[index][1]);
    glutInitWindowSize(size[index][0], size[index][1]);
    winId[index] = glutCreateWindow(" ");
    PR("Window %d id = %d \n", index, winId[index]);
    gfxInit(index);

    addCallbacks();

    sprintf(str, "window %d (RGB)", index);
    glutSetWindowTitle(str);
    sprintf(str, "icon %d", index);
    glutSetIconTitle(str);
    glutSetMenu(menu1);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    break;
  }
}

void idleFunc(void)
{
  int i;
/*
  for (i = 0; i < MAXWIN; i++) {
    if (winId[i] && winVis[i] ) {
      glutSetWindow(winId[i]);
      glutPostRedisplay();
    }
  }
*/
}
void makeMenus(void)
{

/* General control / debug */

  menu2 = glutCreateMenu(menuFunc);
  glutAddMenuEntry("3x3", 303);
  glutAddMenuEntry("5x5", 305);
  glutAddMenuEntry("7x7", 307);
  glutAddMenuEntry("9x9", 309);
  glutAddMenuEntry("11x11", 311);

/* Shapes */

  menu3 = glutCreateMenu(menuFunc);
  glutAddMenuEntry("1", 201);
  glutAddMenuEntry("3", 203);
  glutAddMenuEntry("5", 205);
  glutAddMenuEntry("10", 210);
  glutAddMenuEntry("input", 211);

/* Open/close windows */

  menu4 = glutCreateMenu(menuFunc);
  glutAddMenuEntry("open all windows", 450);
  glutAddMenuEntry("close all windows", 451);
  glutAddMenuEntry(" ", 9999);
  glutAddMenuEntry("create win 0", 400);
  glutAddMenuEntry("create win 1", 401);
  glutAddMenuEntry("create win 2", 402);
  glutAddMenuEntry("create win 3", 403);
  glutAddMenuEntry("create sub window", 404);
  glutAddMenuEntry("create color index win 6", 406);
  glutAddMenuEntry("create color index win 7", 407);
  glutAddMenuEntry(" ", 9999);
  glutAddMenuEntry("destroy win 0", 410);
  glutAddMenuEntry("destroy win 1", 411);
  glutAddMenuEntry("destroy win 2", 412);
  glutAddMenuEntry("destroy win 3", 413);
  glutAddMenuEntry("destroy sub window", 414);
  glutAddMenuEntry("destroy color index win 6", 416);
  glutAddMenuEntry("destroy color index win 7", 417);

/* Window manager stuff */

  menu5 = glutCreateMenu(menuFunc);
  glutAddMenuEntry("move current win", 430);
  glutAddMenuEntry("resize current win", 431);
  glutAddMenuEntry("iconify current win", 432);
  glutAddMenuEntry("show current win", 433);
  glutAddMenuEntry("hide current win", 434);
  glutAddMenuEntry("push current win", 435);
  glutAddMenuEntry("pop current win", 436);
  glutAddMenuEntry(" ", 9999);
  glutAddMenuEntry("move win 1", 420);
  glutAddMenuEntry("resize win 1", 421);
  glutAddMenuEntry("iconify win 1", 422);
  glutAddMenuEntry("show win 1", 423);
  glutAddMenuEntry("hide win 1", 424);
  glutAddMenuEntry("push win 1", 425);
  glutAddMenuEntry("pop win 1", 426);


/* Main menu */

  menu1 = glutCreateMenu(menuFunc);
  glutAddMenuEntry("synthesize", 101);
  glutAddMenuEntry("toggle anotherpass", 103);
  glutAddMenuEntry("toggle exactcolor", 106);
  glutAddMenuEntry("toggle eraser", 107);
  glutAddSubMenu("brush size", menu2);
  glutAddSubMenu("number of passes", menu3);
  glutAddMenuEntry("write result", 102);
  glutAddMenuEntry("write target", 104);
  glutAddMenuEntry("reset target", 105);
  glutAddMenuEntry("quit", 100);
}

void menuFunc(int value){
 char outfname[100];
 int i,j;
 switch (value) {
  case 100:
   exit(0);
  case 101:
   int starttime = time(NULL);
   for(i=0;i<niter;i++){
    create_texture(image, result, data);
    printf("Iteration %d finished\n",i);
   }
   printf("Time per iteration: %f seconds",((float)(time(NULL)-starttime))/niter);
   printf("Begin display %d %d\n", data.widthout, data.heightout);
   glutSetWindow(winId[2]);
   glDrawPixels( data.widthout, data.heightout, GL_RGB, GL_FLOAT,
                     (const GLvoid *) result );
   printf("Done display\n");
   break;
  case 102:
   printf("File name to store the result:");
   scanf("%s",outfname);
   write_ppm(data.widthout, data.heightout, outfname,result);
   break;
  case 104:
   printf("File name to store the target:");
   scanf("%s",outfname);
   write_ppm(data.widthout, data.heightout, outfname,target);
   break;
  case 103:
   anotherpass = !anotherpass;
   printf("anotherpass is now %d\n", anotherpass);
   break;   
  case 105:
   anotherpass = 0;
   printf("anotherpass is now %d\n", anotherpass);
   for(i=0;i<data.heightout;i++){
    for(j=0;j<data.widthout;j++){
     target[a(j,i,data.widthout)+R] = 1.0;
     target[a(j,i,data.widthout)+G] = 1.0;
     target[a(j,i,data.widthout)+B] = 1.0;
    }
   }
   glutSetWindow(winId[1]);
   glRasterPos2i(0,0);
   glDrawPixels( data.widthout, data.heightout, GL_RGB, GL_FLOAT,
                     (const GLvoid *) target );
   break;
  case 106:
   exactcolor = !exactcolor;
   printf("exactcolor is now %d\n", exactcolor);
   break;
  case 107:
   for(i=-brushsizey/2;i<=brushsizey/2;i++){
    for(j=-brushsizex/2;j<=brushsizex/2;j++){
     brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+R] = 1.0;
     brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+G] = 1.0;
     brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+B] = 1.0;
    }
   }
   break;
  case 303:
   brushsizex = brushsizey = 3;
   break;
  case 305:
   brushsizex = brushsizey = 5;
   break;
  case 307:
   brushsizex = brushsizey = 7;
   break;
  case 311:
   brushsizex = brushsizey = 11;
   break;
  case 201:
   niter = 1;
   break;
  case 203:
   niter = 3;
   break;
  case 205:
   niter = 5;
   break;
  case 210:
   niter = 10;
   break;
  case 211:
   printf("Number of iterations to run before displaying the result:");
   scanf("%d",&niter);
   break;
 }
}

void menuStateFunc(int state)
{
}

void setInitDisplayMode(void)
{
  int i;

  displayMode = 0;

  for (i = 0; i < MODES; i++) {
    if (modes[i]) {
      /* printf("Requesting %s \n", modeNames[i]);  */
      displayMode |= glutMode[i];
    }
  }

  glutInitDisplayMode(displayMode);

  if (!glutGet(GLUT_DISPLAY_MODE_POSSIBLE))
    printf("Warning:This display mode not supported\n");
}

void gfxInit(int index)
{

#define XX  3
#define YY  3
#define ZZ  -2.5

  float vertex[][3] =
  {
    {-XX, -YY, ZZ},
    {+XX, -YY, ZZ},
    {+XX, +YY, ZZ},
    {-XX, +YY, ZZ}
  };

/* warning: This func mixes RGBA and CMAP calls in an ugly
   fashion */

  //redefineShapes(currentShape);  /* set up display lists  */
  glutSetWindow(winId[index]);  /* hack - redefineShapes
                                   changes glut win */
  int winwidth,winheight;
  if(index == 0){
   winwidth = data.widthin; winheight = data.heightin;
  }else{
   winwidth = data.widthout; winheight = data.heightout;
  }
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, (GLfloat)winwidth, 0.0, (GLfloat)winheight);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

/* Set proj+view */
/*
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, 1.0, 1.0, 20.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.);
  glTranslatef(0.0, 0.0, -1.0);
*/

/* Set basic material, lighting for RGB windows */

/*
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
*/
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_DITHER);
  glDisable(GL_ALPHA_TEST);
  return;
}

void addCallbacks(void)
{
  glutDisplayFunc(drawScene);
  glutVisibilityFunc(visible);
  glutReshapeFunc(reshapeFunc);
  glutKeyboardFunc(keyFunc);
  glutSpecialFunc(specialFunc);
  glutMouseFunc(mouseFunc);
  glutMotionFunc(motionFunc);
  glutEntryFunc(entryFunc);

/* Callbacks for exotic input devices. Must get my dials &
   buttons back. */
 /*
  glutSpaceballMotionFunc(spaceballMotionCB);
  glutSpaceballRotateFunc(spaceballRotateCB);
  glutSpaceballButtonFunc(spaceballButtonCB);

  glutButtonBoxFunc(buttonBoxCB);
  glutDialsFunc(dialsCB);

  glutTabletMotionFunc(tabletMotionCB);
  glutTabletButtonFunc(tabletButtonCB);
 */
}

void drawScene(void)
{
  int winIndex;
  glClear(GL_COLOR_BUFFER_BIT);
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  winIndex = idToIndex(glutGetWindow());
  switch(winIndex){
   case 0:
    glDrawPixels( data.widthin, data.heightin, GL_RGB, GL_FLOAT,
                     (const GLvoid *) image );
    break;
   case 1:
    glRasterPos2i(0,0);
    glDrawPixels( data.widthout, data.heightout, GL_RGB, GL_FLOAT,
                     (const GLvoid *) target );
    break;
   case 2:
    glDrawPixels( data.widthout, data.heightout, GL_RGB, GL_FLOAT,
                     (const GLvoid *) result );
    break;
  }
  //glutSwapBuffers();
}

void reshapeFunc(int width, int height) {}

void visible(int state) {
  int winid,i;
  static GLboolean someVisible = GL_TRUE;

  winid = glutGetWindow();
  /* printf("visible: state = %d \n", state);  */

  if (state == GLUT_VISIBLE) {
    PR("Window id %d visible \n", winid);
    winVis[idToIndex(winid)] = GL_TRUE;
  } else {
    PR("Window %d not visible \n", winId);
    winVis[idToIndex(winid)] = GL_FALSE;
  }

  if ((winVis[0] == GL_FALSE) && (winVis[1] == GL_FALSE) &&
      (winVis[2] == GL_FALSE)) {
    glutIdleFunc(NULL);
    PR("All windows not visible; idle func disabled\n");
    someVisible = GL_FALSE;
  } else {
      PR("Some windows now visible; idle func enabled\n");
      someVisible = GL_TRUE;
      for (i = 0; i < MAXWIN; i++) {
       if (winId[i] && winVis[i] ) {
        glutSetWindow(winId[i]);
        glutPostRedisplay();
       }
      }
  }
}

void keyFunc(unsigned char key, int x, int y) {}
void mouseFunc(int button, int state, int x, int y){
 int wind,i,j;
 if(button == GLUT_LEFT_BUTTON){
  wind = idToIndex(glutGetWindow());
 if(wind == 0){
  if(x > data.widthin-data.localx/2-1)
    cursourcex = data.widthin-data.localx/2-1;
  else if(x < data.localx/2)
    cursourcex = data.localx/2;
  else cursourcex = x;
  if(y > data.heightin-data.localy/2-1)
    cursourcey = data.localy/2;
  else if(y < data.localy/2)
    cursourcey = data.heightin-data.localy/2-1;
  else cursourcey = data.heightin - 1 - y;
  pixelvalue resultr=0,resultg=0,resultb=0;
  if(exactcolor){
   resultr = image[a(cursourcex,cursourcey,data.widthin)+R];
   resultg = image[a(cursourcex,cursourcey,data.widthin)+G];
   resultb = image[a(cursourcex,cursourcey,data.widthin)+B];
  }else{
   for(i=-brushsizey/2;i<=brushsizey/2;i++){
    for(j=-brushsizex/2;j<=brushsizex/2;j++){
     resultr += 
      image[a(cursourcex+j,cursourcey+i,data.widthin)+R];
     resultg +=
      image[a(cursourcex+j,cursourcey+i,data.widthin)+G];
     resultb +=
      image[a(cursourcex+j,cursourcey+i,data.widthin)+B];
    }
   }
   resultr /= brushsizex*brushsizey;
   resultg /= brushsizex*brushsizey;
   resultb /= brushsizex*brushsizey;
  }
  for(i=-brushsizey/2;i<=brushsizey/2;i++){
   for(j=-brushsizex/2;j<=brushsizex/2;j++){
    brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+R] = resultr;
    brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+G] = resultg;
    brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+B] = resultb;
   }
  }
  return;
 }else if(wind == 1){
  motionFunc(x, y);
 }
 }
}
void motionFunc(int x, int y) {
 int wind = idToIndex(glutGetWindow());
 int i,j;
 if(wind == 1){
  for(i=-brushsizey/2;i<=brushsizey/2;i++){
   for(j=-brushsizex/2;j<=brushsizex/2;j++){
    target[a(x+j,data.heightout-y+i,data.widthout)+R] = brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+R];
    target[a(x+j,data.heightout-y+i,data.widthout)+G] = brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+G];
    target[a(x+j,data.heightout-y+i,data.widthout)+B] = brush[a(brushsizex/2+j,brushsizey/2+i,brushsizex)+B];
    //atlas[aa(x+j,data.heightout-y+i)] = cursourcex;   //+j;
    //atlas[aa(x+j,data.heightout-y+i)+1] = cursourcey; //+i;
   }
  }
  glRasterPos2i(x-brushsizex/2,data.heightout-y-brushsizey/2);
  glDrawPixels( brushsizex, brushsizey, GL_RGB, GL_FLOAT,
                     (const GLvoid *) brush );
  
 }
}
void passiveMotionFunc(int x, int y) {}
void entryFunc(int state) {}
void specialFunc(int key, int x, int y) {}
int idToIndex(int id)
{
  int i;
  for (i = 0; i < MAXWIN; i++) {
    if (winId[i] == id)
      return i;
  }
  fprintf(stderr, "error: id %d not found \n", id);
  return (-1);
}
