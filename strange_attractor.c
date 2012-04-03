/*============================================================================*
 * strange_attactor.c
 * Time-stamp: <2011-12-18 23:50:49 (mkmccjr)>
 *
 * Compile: "gcc -W -Wall -std=c99 -pedantic -O3 \
 *               -o fractal.x fractal.c -lm"
 *
 *============================================================================*/
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.c"
#include "arrays.c"
#include "color_constants.h"
#include "color_conversion.c"
#include "par.c"



/* Macro definitions to for efficiency.
 * -------------------------------------------------------------------------- */
#define MAX(x, y) (x > y) ? x : y
#define MIN(x, y) (x < y) ? x : y


/* Constants for the fractal iteration
 * -------------------------------------------------------------------------- */
static int iter_max;
static double aa, bb, cc, dd;
static double expo, cut;
static char *method;


/* Image parameters
 * -------------------------------------------------------------------------- */
static int Nx, Ny;
static double xmin, xmax, ymin, ymax, Lx, Ly, dx, dy;
static double center[2];


/* Other
 * -------------------------------------------------------------------------- */
static char *filename;
static void process_input(char *athinput, int argc, char *argv[]);

static inline double (*xnew) (double, double);
static inline double (*ynew) (double, double);

static inline double xnew_clifford(double x, double y);
static inline double ynew_clifford(double x, double y);
static inline double xnew_svensson(double x, double y);
static inline double ynew_svensson(double x, double y);
static inline double xnew_peter(double x, double y);
static inline double ynew_peter(double x, double y);

static void render_fractal(int **sample_data);




/* Main function
 * -------------------------------------------------------------------------- */
int main (int argc, char *argv[])
{
  int i, j;
  int max;

  /* image stuff */
  FILE *fp;
  char *comment  = "# ";

  int **array;

  char *definput = "input.frac";  /* default input filename */
  char *athinput = definput;


  /* parse command line options */
  for (i=1; i<argc; i++) {
    switch(*(argv[i]+1)) {
    case 'i':                      /* -i <file>   */
      athinput = argv[++i];
      break;
    default:
      break;
    }
  }

  /* process the input file */
  process_input(athinput, argc, argv);


  /* render the fractal */
  array = (int**) calloc_2d_array(Ny, Nx, sizeof(int));

  render_fractal(array);

  max = 0;
  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      max = MAX(max, array[j][i]);
    }
  }

  double val;

  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      val = cut * (double) array[j][i] / max;
      val = MAX(0,val);  val = MIN(1,val);
      val = pow(val, expo);

      array[j][i] = 255 * (1.0-val);
    }
  }

  /* Write to a file */
  fp = fopen(filename, "wb");
  /* write header to the file */
  fprintf(fp, "P5\n %s\n %d\n %d\n %d\n", comment, Nx, Ny, 255);

  /* write image data bytes to the file */
  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      fwrite(&array[j][i], sizeof(char), 1, fp);
    }
  }

  /* clean up */
  fclose(fp);
  free_2d_array((void**) array);

  return 0;
}
/* End main
 * -------------------------------------------------------------------------- */



/* Process the input file
 * -------------------------------------------------------------------------- */
static void process_input(char *athinput, int argc, char *argv[])
{
  /* first comes first
     -------------------------------------------------- */
  par_open(athinput);
  par_cmdline(argc,argv);

  filename = par_gets_def("image", "file", "attractor.ppm");

  Nx = par_geti("image", "Nx");
  Ny = par_geti("image", "Ny");

  iter_max = (int) par_getd("image", "npts");


  /* fractal properties
     -------------------------------------------------- */
  center[0] = par_getd("fractal", "center_x");
  center[1] = par_getd("fractal", "center_y");

  Lx = par_getd("fractal", "Lx");
  Ly = Lx * Ny / Nx;

  xmin = center[0] - Lx/2;  xmax = center[0] + Lx/2;
  ymin = center[1] - Ly/2;  ymax = center[1] + Ly/2;

  dx = Lx/Nx;
  dy = Ly/Ny;



  /* image properties
     -------------------------------------------------- */
  aa = par_getd("fractal", "a");
  bb = par_getd("fractal", "b");
  cc = par_getd("fractal", "c");
  dd = par_getd("fractal", "d");

  cut  = par_getd_def("image", "cut", 10.0);
  expo = par_getd_def("image", "exp",  0.5);

  method = par_gets_def("fractal", "method", "clifford");

  if (strcmp(method, "peter") == 0) {
    xnew = &xnew_peter;
    ynew = &ynew_peter;
  } else if (strcmp(method, "svensson") == 0) {
    xnew = &xnew_svensson;
    ynew = &ynew_svensson;
  } else{
    xnew = &xnew_clifford;
    ynew = &ynew_clifford;
  }

  par_close();
  return;
}
/* End process input
 * -------------------------------------------------------------------------- */



/* Different equations
 * -------------------------------------------------------------------------- */
static inline double xnew_clifford(double x, double y) {
  return (sin(aa*y) + cc * cos(aa*x));
}

static inline double ynew_clifford(double x, double y) {
  return (sin(bb*x) + dd * cos(bb*y));
}

static inline double xnew_peter(double x, double y) {
  return (sin(aa*y) - cos(bb*x));
}

static inline double ynew_peter(double x, double y) {
  return (sin(cc*x) - cos(dd*y));
}

static inline double xnew_svensson(double x, double y) {
  return (dd * sin(aa*x) - sin(bb*y));
}

static inline double ynew_svensson(double x, double y) {
  return (cc * cos(aa*x) + cos(bb*y));
}
/* End equations
 * -------------------------------------------------------------------------- */



/* Render the fractal
 * -------------------------------------------------------------------------- */
static void render_fractal(int **sample_data)
{
  int i, j ,c;
  double x, y, x1, y1;

  x = y = 0.1;
  /* "burn in" the orbit. */
  for (c = 0; c < 50; c++){
    x1 = xnew(x, y);  y1 = ynew(x, y);
    x  = x1;          y  = y1;
  }

  /* now, iterate the fractal */
  for (c = 0; c < iter_max; c++){
    x1 = xnew(x, y);  y1 = ynew(x, y);
    x  = x1;          y  = y1;

    i = (int) ((x-xmin)/dx);
    j = (int) ((y-ymin)/dy);

    i = MIN(Nx-1, i);  i = MAX(0, i);
    j = MIN(Ny-1, j);  j = MAX(0, j);

    sample_data[j][i]++;
  }

  return;
}
/* End render fractal
 * -------------------------------------------------------------------------- */

/*
Local Variables:
  compile-command:
    "gcc -W -Wall -std=c99 -pedantic -O3 \
         -o strange_attractor.x strange_attractor.c -lm"
End:
*/
