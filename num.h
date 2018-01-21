#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#include "mklwrapper.h"

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define SGN(a)    ((a) > 0 ? 1.0 : -1.0)

/* typedef struct FCOMPLEX {double r, i; } fcomplex; */


double pow_di( double x, int n );
int pow_i( int x, int n );
void gauleg(double x1, double x2, double x[], double w[], int n);

void nrerror(char error_text[]);
void gaulag(double x[], double w[], int n, double alf);
void gaujac(double x[], double w[], int n, double alf, double bet);
double gammln(double xx);
double plgndr(int l, int m, double x);
double qlgndr(int l, int m, double x);
void free_c3tensor(double _Complex ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);
double _Complex ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);
char **chmatrix(long nrl, long nrh, long ncl, long nch);
double *vector(long nl, long nh);
double _Complex *cvector(long nl, long nh);
void free_cvector(double _Complex *v, long nl, long nh);
void free_vector(double *v, long nl, long nh);
void free_cvector(double _Complex *v, long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
double _Complex **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_cmatrix(double _Complex **m, long nrl, long nrh, long ncl, long nch);
void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch);

double _Complex c_conv_int2(int p, int N, double H, double _Complex *A, double _Complex *B);
double _Complex rc_conv_int2(int p, int N, double H, double *A, double _Complex *B);
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);
void fourn(double data[], int nn[], int ndim, int isign);
/*fcomplex Cdiv(fcomplex a, fcomplex b);
 fcomplex Cmul(fcomplex a, fcomplex b);
 fcomplex RCmul(double a, fcomplex b);
 fcomplex Csub(fcomplex a, fcomplex b);
 fcomplex Cadd(fcomplex a, fcomplex b);
 fcomplex Cexp(double a);
 fcomplex Conj(fcomplex a);
 double Cabs(fcomplex a);
 */
void set_integration_constants();
double integr(double *fg, int a, int b, double H);
double _Complex cintegr(double _Complex *fg, int a, int b, double H);
double emlint1(double *fint, int nxmin, int nxmax, double HH);
double _Complex cemlint1(double _Complex *fint, int nxmin, int nxmax, double HH);
double emlint(double **fint, int nxmin, int nxmax, int nymin, int nymax,
              double HH);


void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void splint_der(double xa[], double ya[], double y2a[], int n, double x, double *y);



void  printvec(char *fnm, int n, double *x, double *y);
void  printcol(char *fnm, int n, int row, double *x, double **y);
void time_deriv(int tpts, double *tgrid, double *f, double yd1, 
                double ydn, double *df, double *ddf);
double **rallocation(int M, int *mdim);
double _Complex **callocation(int M, int *mdim);
double _Complex ***cmat_alloc(int M, int *mdim, int tdim);
double _Complex **cvec_alloc(int M, int *mdim);
void free_cmat(double _Complex ***mat, int M, int *mdim, int tdim);
void swap(double *v, int i, int j);
void opp_sort(double *v, int left, int right);
void sort(double *v, int left, int right);
double qromb(double (*func)(double, int, double *, double *, double *), double a, double b, int N, double *x, double *y, double *y2);
double trapzd(double (*func)(double, int, double *, double *, double *), double a, double b, int n, int N, double *x, double *y, double *y2);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void conv_int(int tpts, double *tgrid, double *ya, double *yb, double *yc);
double conv_int2(int p, int tpts, double *tgrid, double *ya, double *yb);
double spline_int(int n1, int n2, double *x, double *y);
double _Complex cspline_int(int n1, int n2, double *x, double _Complex *y);
double _Complex cxspline_int(int n1, int n2, double *x, double _Complex *y);

void cr_conv_int(int tpts, double *tgrid, double _Complex *ya, double *yb,double _Complex *yc);
void cr_conv_int_W(int tpts, double *tgrid, double _Complex *ya, double *yb, double _Complex *yc);
void rc_conv_int_W(int tpts, double *tgrid, double *ya, double _Complex *yb, double _Complex *yc);
