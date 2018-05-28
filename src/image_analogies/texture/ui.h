#ifndef UI_H
#define UI_H

#include <GL/glut.h>
#define MAXWIN     9    /* max no. of windows */
GLboolean debug = GL_FALSE;  /* dump all events  */
int menu1, menu2, menu3, menu4, menu5, menu6 = 0, menu7, menu8;
int pos[MAXWIN][2] =
{
  {50, 150},            /* win 0  */
  {450, 150},           /* win 1  */
  {50, 600},            /* win 2  */
  {450, 600},           /* win 3  */
  {10, 10},             /* subwin 4 (relative to parent win 0) */
  {300, 400},           /* help win 5  */
  {850, 150},           /* cmap win 6  */
  {850, 600},           /* cmap win 7  */
  {250, 450}            /* text win 8  */
};

int size[MAXWIN][2] =
{
  {350, 350},           /* win 0  */
  {350, 350},           /* win 1  */
  {350, 350},           /* win 2  */
  {350, 350},           /* win 3  */
  {200, 200},           /* subwin 4  */
  {700, 300},           /* help win 5  */
  {350, 350},           /* cmap win 6  */
  {350, 350},           /* cmap win 7  */
  {800, 450}            /* text win 8  */
};
#define PR     if(debug)printf
int winId[MAXWIN] =
{0};                    /* table of glut window id's  */
GLboolean winVis[MAXWIN] =
{GL_FALSE};                /* is window visible  */

enum {
  RGBA, INDEX, SINGLE, DOUBLEBUFFER, DEPTH, ACCUM, ALPHA, STENCIL, MULTISAMPLE,
  STEREO, MODES
};
int glutMode[] =
{GLUT_RGBA, GLUT_INDEX, GLUT_SINGLE, GLUT_DOUBLE, GLUT_DEPTH,
  GLUT_ACCUM, GLUT_ALPHA, GLUT_STENCIL, GLUT_MULTISAMPLE, GLUT_STEREO};
int modes[MODES] = {0};
int displayMode = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH;

void menuStateFunc(int state);
void idleFunc(void);
void setInitDisplayMode(void);
void makeMenus(void);
void makeWindow(int index);
void gfxInit(int index);
void addCallbacks(void);
void menuFunc(int value);
int idToIndex(int id);

void drawScene(void);
void reshapeFunc(int width, int height);
void visible(int state);
void keyFunc(unsigned char key, int x, int y);
void mouseFunc(int button, int state, int x, int y);
void motionFunc(int x, int y);
void passiveMotionFunc(int x, int y);
void entryFunc(int state);
void specialFunc(int key, int x, int y);

#endif
