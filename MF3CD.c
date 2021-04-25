/* Dey, Subhadip & Bhattacharya, Avik & Frery, Alejandro & López-Martínez, Carlos. (2020).
A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data.

To compile and run, e.g.:
   gcc MF3CD.c -o MF3CD.exe -lm
   ./MF3CD.exe T2

C impl. 20210423 by Ash Richardson, Senior Data Scientist, BC Wildfire Service */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include<string.h>
#define for0(i,n) for(i = 0; i < n; i++) /* for loop shorthand */

#define N_IN 4 /* number of input files */
#define T11 0 /* T3 matrix input data indexing */
#define T12_re 1
#define T12_im 2
#define T22 3

const char* T_fn[] = {"T11.bin", /* T2 matrix input filenames */
                      "T12_real.bin",
                      "T12_imag.bin",
                      "T22.bin"};

float ** T; /* buffers for T matrix elements */
float ** out_d; /* output data buffers */

#define N_OUT 5 /* number of output bands */
#define _dop 0 /* output data indexing */
#define _theta 1
#define _pd 2
#define _pv 3
#define _ps 4

const char * out_fn[] = {"DoP_DP.bin", /* output filenames */
                         "theta_DP.bin",
                         "pd_DP.bin",
                         "pv_DP.bin",
                         "ps_DP.bin"};
char sep(){
  #ifdef _WIN32
    return '\\'; /* windows path separator */
  #else
    return '/'; /* posix path sep */
  #endif
}

void err(const char * m){
  printf("Error: %s\n", m);
  exit(1); /* print message and bail */
}

#define MAX_ARRAYS 1024 /* large number */
int n_arrays = 0; /* count arrays initialized */
void ** arrays; /* track mallocs, free at end */

void * alloc(size_t n){
  void * d = malloc(n); /* create array */
  if(!d) err("failed to allocate memory");
  memset(d, '\0', n); /* must touch memory on windows */
  arrays[n_arrays ++] = d; /* save pointer to free later */
  return d;
}

float * falloc(size_t n){
  return (float *)alloc(n * sizeof(float)); /* float32 array */
}

#define READ 1
#define WRITE 0

FILE * open(const char * fn, int mode){
  printf("+%s %s\n", mode?"r":"w", fn); /* open a file for read or write */
  FILE * f = fopen(fn, mode?"rb":"wb");
  if(!f) err("file access failed");
  return f;
}

void read_config(char * file_name, int * nrow, int * ncol){
  size_t x;
  char tmp[4096]; /* based on PolSARPro by Eric POTTIER and Laurent FERRO-FAMIL */
  FILE * f = open(file_name, READ);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp); // number of rows
  *nrow = atoi(tmp);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp); // number of cols
  *ncol = atoi(tmp);
  fclose(f);
  printf("nrow %d ncol %d\n", *nrow, *ncol);
}

float * read(const char * file_name, size_t n_float){
  FILE * f = open(file_name, READ);
  float * d = falloc(n_float);
  size_t nr = fread(d, sizeof(float), n_float, f);
  if(nr != n_float){
    printf("Expected number of floats: %zu\n", n_float);
    printf("Number of floats read: %zu\n", nr);
    err("Unexpected float read count");
  }
  fclose(f);
  return d; /* return array of floats we read in */
}

#define STR_MAX 4096 /* string variable length */
void hwrite(char * bfn, size_t nrow, size_t ncol, size_t nband){
  size_t i;
  char hfn[STR_MAX];
  size_t L = strlen(bfn);

  strcpy(hfn, bfn); /* change ext from bin to hdr*/
  hfn[L - 3] = 'h';
  hfn[L - 2] = 'd';
  hfn[L - 1] = 'r';

  FILE * f = open(hfn, WRITE);
  fprintf(f, "ENVI\n");
  fprintf(f, "samples = %zu\n", ncol);
  fprintf(f, "lines = %zu\n", nrow);
  fprintf(f,"bands = %zu\n", nband);
  fprintf(f, "header offset = 0\n");
  fprintf(f, "file type = ENVI Standard\n");
  fprintf(f, "data type = 4\n");
  fprintf(f, "interleave = bsq\n");
  fprintf(f, "byte order = 0\n");
  fprintf(f, "band names = {band 1");
  for0(i, nband - 1) fprintf(f, ",\nband %zu", i + 2);
  fprintf(f, "}\n");
  fclose(f);
}

int main(int argc, char ** argv){

  if(argc < 2) err("M4FC.exe [input T3 directory]");
  char * path = argv[1]; /* T3 matrix data path */
  int i, j, k, np, nrow, ncol, di, dj, ii, jj, x, ix, jx, nw;

  char fn[STR_MAX];
  strcpy(fn, path);
  fn[strlen(path)] = sep();
  strcpy(fn + strlen(path) + 1, "config.txt"); /* path to config.txt */
  read_config(fn, &nrow, &ncol); /* read image dimensions */
  np = nrow * ncol; /* number of px */

  arrays = (void *) malloc(sizeof(void *) * MAX_ARRAYS); /* array of pointers to free later */
  memset(arrays, '\0', sizeof(void *) * MAX_ARRAYS); /* always touch memory on win OS */
  n_arrays = 0; /* start from the beginning */

  T = (float **) alloc(sizeof(float *) * N_IN); /* input file buffers */
  for0(k, N_IN){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, T_fn[k]); /* [path][sep][filename] e.g. T3/T11.bin */
    T[k] = read(fn, np); /* read each input data band */
  }

  out_d = (float **) alloc(sizeof(float *) * N_OUT); /* output bands buffers */
  for0(i, N_OUT) out_d[i] = falloc(np); /* allocate output space */

  double t11, t12_r, t12_i, t22; /* intermediary variables */
  double dop_d, theta_d, pd_d, pv_d, ps_d;
  double det, m1, trace, val_d, trace2, h, g;
  float * out_d_dop, * out_d_theta, * out_d_pd, * out_d_pv, * out_d_ps;

  float * t11_p = T[T11];
  float * t12_r_p = T[T12_re];
  float * t12_i_p = T[T12_im];
  float * t22_p = T[T22];

  out_d_dop = out_d[_dop];
  out_d_theta = out_d[_theta];
  out_d_pd = out_d[_pd];
  out_d_pv = out_d[_pv];
  out_d_ps = out_d[_ps];

  for0(i, np){
    t11   = (double)t11_p[i];
    t12_r = (double)t12_r_p[i];
    t12_i = (double)t12_i_p[i];
    t22   = (double)t22_p[i];

    trace = t11 + t22; // s0_d
    trace2 = trace * trace;
    det = t11 * t22 - (t12_r * t12_r) - (t12_i * t12_i);
    dop_d = sqrt(1. - (4. * det / trace2));
    pv_d = (1. - dop_d) * trace;
    h = t11 - t22;
    g = t11;

    val_d = (dop_d * trace * h) / (t22 * g + dop_d * dop_d * trace2);
    theta_d = atan(val_d);

    pd_d = dop_d * (trace / 2.) * (1. - sin(2. * theta_d));
    ps_d = dop_d * (trace / 2.) * (1. + sin(2. * theta_d));

    out_d_dop[i] = (float)dop_d;
    out_d_theta[i] = (float)(theta_d * 180. / M_PI);
    out_d_pd[i] = (float)pd_d;
    out_d_pv[i] = (float)pv_d;
    out_d_ps[i] = (float)ps_d;
  }

  FILE * out_f[N_OUT];
  for0(i, N_OUT){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, out_fn[i]);
    out_f[i] = open(fn, WRITE);
    hwrite(fn, nrow, ncol, 1); /* write envi header */
    nw = fwrite(out_d[i], sizeof(float), np, out_f[i]);
    if(nw != np) err("failed to write expected number of floats");
    fclose(out_f[i]);
  }

  for0(k, n_arrays) free(arrays[n_arrays]); /* free anything we malloc'ed */
  free(arrays);
  return 0;
}