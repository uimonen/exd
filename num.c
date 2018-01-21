#include "num.h"

#define NR_END 1
#define FREE_ARG char*

static double *_c;
static char EML_INIT = 1;

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

double _Complex *cvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double _Complex *v;

    v = (double _Complex *) malloc((size_t) ((nh - nl + 1 + NR_END)
            * sizeof(double _Complex)));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    free((FREE_ARG) (v + nl - NR_END));
    v = NULL;
}

void free_cvector(double _Complex *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    free((FREE_ARG) (v + nl - NR_END));
    v = NULL;
}

int *ivector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    int *v;

    v = (int *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free a int vector allocated with vector() */
{
    free((FREE_ARG) (v + nl - NR_END));
    v = NULL;
}

long *lvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    long *v;

    v = (long *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_lvector(long *v, long nl, long nh)
/* free a int vector allocated with vector() */
{
    free((FREE_ARG) (v + nl - NR_END));
    v = NULL;
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

char **chmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    char **m;

    /* allocate pointers to rows */
    m -= nrl;
    m = (char **) malloc((size_t) ((nrow + NR_END) * sizeof(char*)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(char)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m -= nrl;
    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double*)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *) malloc((size_t) ((nrow * ncol + NR_END)
            * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    m = NULL;
}


double _Complex **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a fcomplex matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double _Complex **m;

    /* allocate pointers to rows */
    m -= nrl;
    m = (double _Complex **) malloc((size_t) ((nrow + NR_END)
            * sizeof(double _Complex*)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double _Complex *) malloc((size_t) ((nrow * ncol + NR_END)
            * sizeof(double _Complex)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_cmatrix(double _Complex **m, long nrl, long nrh, long ncl, long nch)
/* free a double _Complex matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    m = NULL;
}

void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    m = NULL;
}

double _Complex ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
        long ndh)
/* allocate a double _Complex 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double _Complex ***t;

    /* allocate pointers to pointers to rows */
    t = (double _Complex ***) malloc((size_t) ((nrow + NR_END)
            * sizeof(double _Complex**)));
    if (!t)
        nrerror("allocation failure 1 in d3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double _Complex **) malloc((size_t) ((nrow * ncol + NR_END)
            * sizeof(double _Complex*)));
    if (!t[nrl])
        nrerror("allocation failure 2 in d3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double _Complex *) malloc((size_t) ((nrow * ncol * ndep
            + NR_END) * sizeof(double _Complex)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in d3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

void free_c3tensor(double _Complex ***t, long nrl, long nrh, long ncl,
        long nch, long ndl, long ndh)
/* free a double _Complex d3tensor allocated by d3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double ***t;

    /* allocate pointers to pointers to rows */
    t = (double ***) malloc((size_t) ((nrow + NR_END) * sizeof(double**)));
    if (!t)
        nrerror("allocation failure 1 in d3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double **) malloc((size_t) ((nrow * ncol + NR_END)
            * sizeof(double*)));
    if (!t[nrl])
        nrerror("allocation failure 2 in d3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double *) malloc((size_t) ((nrow * ncol * ndep + NR_END)
            * sizeof(double)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in d3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int **m;

    /* allocate pointers to rows */
    m -= nrl;
    m = (int **) malloc((size_t) ((nrow + NR_END) * sizeof(int*)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (int *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free a int matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
  m = NULL;
}


int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  int ***t;
  
  /* allocate pointers to pointers to rows */
  t = (int ***) malloc((size_t) ((nrow + NR_END) * sizeof(int**)));
  if (!t)
    nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (int **) malloc((size_t) ((nrow * ncol + NR_END)
                                        * sizeof(int*)));
  if (!t[nrl])
    nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (int *) malloc((size_t) ((nrow * ncol * ndep + NR_END)
                                            * sizeof(int)));
  if (!t[nrl][ncl])
    nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  
  /* return pointer to array of pointers to rows */
  return t;
}


void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG) (t[nrl] + ncl - NR_END));
  free((FREE_ARG) (t + nrl - NR_END));

} 
double qgaus(double(*func)(double), double a, double b, double *xf, double *wf,
        int n) {
    int j;
    double s;

    s = 0;
    for (j = 1; j <= n; j++)
        s += wf[j] * ((*func)(xf[j]));

    return s;
}

#define EPS 3.0e-11

void gauleg(double x1, double x2, double x[], double w[], int n) {
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n + 1) >> 1;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for (i = 1; i <= m; i++) {
        z = cos(3.141592654 * (i - 0.25) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > EPS);
        x[i] = xm - xl * z;
        x[n + 1 - i] = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n + 1 - i] = w[i];
    }
}
#undef EPS

#define EPS 3.0e-14
#define MAXIT 50

void gaulag(double x[], double w[], int n, double alf) {
    void nrerror(char error_text[]);
    int i, its, j;
    double ai;
    double p1, p2, p3, pp, z, z1;

    for (i = 1; i <= n; i++) {
        if (i == 1) {
            z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
        } else if (i == 2) {
            z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
        } else {
            ai = i - 2;
            z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf / (1.0 + 3.5
                    * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * alf);
        }
        for (its = 1; its <= MAXIT; its++) {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
            }
            pp = (n * p1 - (n + alf) * p2) / z;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS)
                break;
        }
        if (its > MAXIT)
            nrerror("too many iterations in gaulag");
        x[i] = z;
        w[i] = -exp(gammln(alf + n) - gammln((double) n)) / (pp * n * p2);
    }
}
#undef EPS
#undef MAXIT

double gammln(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = { 76.18009172947146, -86.50532032941677,
            24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
            -0.5395239384953e-5 };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
        ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

#define EPS 3.0e-14
#define MAXIT 10

void gaujac(double x[], double w[], int n, double alf, double bet) {
    double gammln(double xx);
    void nrerror(char error_text[]);
    int i, its, j;
    double alfbet, an, bn, r1, r2, r3;
    double a, b, c, p1, p2, p3, pp, temp, z, z1;

    for (i = 1; i <= n; i++) {
        if (i == 1) {
            an = alf / n;
            bn = bet / n;
            r1 = (1.0 + alf) * (2.78 / (4.0 + n * n) + 0.768 * an / n);
            r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an * bn;
            z = 1.0 - r1 / r2;
        } else if (i == 2) {
            r1 = (4.1 + alf) / ((1.0 + alf) * (1.0 + 0.156 * alf));
            r2 = 1.0 + 0.06 * (n - 8.0) * (1.0 + 0.12 * alf) / n;
            r3 = 1.0 + 0.012 * bet * (1.0 + 0.25 * fabs(alf)) / n;
            z -= (1.0 - z) * r1 * r2 * r3;
        } else if (i == 3) {
            r1 = (1.67 + 0.28 * alf) / (1.0 + 0.37 * alf);
            r2 = 1.0 + 0.22 * (n - 8.0) / n;
            r3 = 1.0 + 8.0 * bet / ((6.28 + bet) * n * n);
            z -= (x[1] - z) * r1 * r2 * r3;
        } else if (i == n - 1) {
            r1 = (1.0 + 0.235 * bet) / (0.766 + 0.119 * bet);
            r2 = 1.0 / (1.0 + 0.639 * (n - 4.0) / (1.0 + 0.71 * (n - 4.0)));
            r3 = 1.0 / (1.0 + 20.0 * alf / ((7.5 + alf) * n * n));
            z += (z - x[n - 3]) * r1 * r2 * r3;
        } else if (i == n) {
            r1 = (1.0 + 0.37 * bet) / (1.67 + 0.28 * bet);
            r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
            r3 = 1.0 / (1.0 + 8.0 * alf / ((6.28 + alf) * n * n));
            z += (z - x[n - 2]) * r1 * r2 * r3;
        } else {
            z = 3.0 * x[i - 1] - 3.0 * x[i - 2] + x[i - 3];
        }
        alfbet = alf + bet;
        for (its = 1; its <= MAXIT; its++) {
            temp = 2.0 + alfbet;
            p1 = (alf - bet + temp * z) / 2.0;
            p2 = 1.0;
            for (j = 2; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                temp = 2 * j + alfbet;
                a = 2 * j * (j + alfbet) * (temp - 2.0);
                b = (temp - 1.0) * (alf * alf - bet * bet + temp * (temp - 2.0)
                        * z);
                c = 2.0 * (j - 1 + alf) * (j - 1 + bet) * temp;
                p1 = (b * p2 - c * p3) / a;
            }
            pp = (n * (alf - bet - temp * z) * p1 + 2.0 * (n + alf) * (n + bet)
                    * p2) / (temp * (1.0 - z * z));
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS)
                break;
        }
        if (its > MAXIT)
            nrerror("too many iterations in gaujac");
        x[i] = z;
        w[i] = exp(gammln(alf + n) + gammln(bet + n) - gammln(n + 1.0)
                - gammln(n + alfbet + 1.0)) * temp * pow(2.0, alfbet) / (pp
                * p2);
    }
}
#undef EPS
#undef MAXIT

#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = (double *) malloc((size_t) ((n + 1) * sizeof(double)));

    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = TINY;
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }
    free(vv);
}
#undef TINY

double plgndr(int l, int m, double x) {
    void nrerror(char error_text[]);
    double fact, pll, pmm, pmmp1, somx2;
    int i, ll;

    if (m < 0 || m > l || fabs(x) > 1.0)
        nrerror("Bad arguments in routine plgndr");
    pmm = 1.0;
    if (m > 0) {
        somx2 = sqrt((1.0 - x) * (1.0 + x));
        fact = 1.0;
        for (i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pow_di(-1.0, m) * pmm;
    else {
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1))
            return pow_di(-1.0, m) * pmmp1;
        else {
            for (ll = m + 2; ll <= l; ll++) {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm)
                        / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pow_di(-1.0, m) * pll;
        }
    }
}

void set_integration_constants() {
    _c[1] = 0.3298611111;
    _c[2] = 1.3208333333;
    _c[3] = 0.7666666667;
    _c[4] = 1.1013888889;
    _c[5] = 0.9812500000;

}

double _Complex cemlint1(double _Complex *fint, int nxmin, int nxmax, double HH) {
    int i;
    double _Complex sum;

    if (EML_INIT) {
        _c = (double *) malloc((size_t) ((5 + 1) * sizeof(double)));
        set_integration_constants();
        EML_INIT = 0;
    }

    sum = 0.0;

    for (i = nxmin; i <= (nxmin + 4); i++)
        sum += _c[i - nxmin + 1] * (fint[i] + fint[nxmax + nxmin - i]);

    for (i = (nxmin + 5); i <= (nxmax - 5); i++)
        sum += fint[i];

    return (HH * sum);

}

double emlint1(double *fint, int nxmin, int nxmax, double HH) {
    int i;
    double sum;

    if (EML_INIT) {
        _c = (double *) malloc((size_t) ((5 + 1) * sizeof(double)));
        set_integration_constants();
        EML_INIT = 0;
    }

    sum = 0.0;

    for (i = nxmin; i <= (nxmin + 4); i++)
        sum += _c[i - nxmin + 1] * (fint[i] + fint[nxmax + nxmin - i]);

    for (i = (nxmin + 5); i <= (nxmax - 5); i++)
        sum += fint[i];

    return (HH * sum);

}

double emlint(double **fint, int nxmin, int nxmax, int nymin, int nymax,
        double HH) {
    int i, j;
    double sum;

    sum = 0.0;

    for (i = nxmin; i <= (nxmin + 4); i++) {
        for (j = nymin; j <= (nymin + 4); j++)
            sum += _c[i - nxmin + 1] * _c[j - nymin + 1] * (fint[i][j]
                    + fint[nxmax - i + nxmin][j] + fint[i][nymax - j + nymin]
                    + fint[nxmax - i + nxmin][nymax - j + nymin]);
        for (j = (nymin + 5); j <= (nymax - 5); j++)
            sum += _c[i - nxmin + 1]
                    * (fint[i][j] + fint[nxmax - i + nxmin][j]);
    }

    for (i = (nxmin + 5); i <= (nxmax - 5); i++) {
        for (j = nymin; j <= (nymin + 4); j++)
            sum += _c[j - nymin + 1]
                    * (fint[i][j] + fint[i][nymax - j + nymin]);
        for (j = (nymin + 5); j <= (nymax - 5); j++)
            sum += fint[i][j];
    }

    return (HH * HH * sum);

}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]) {
    int i, k;
    double p, qn, sig, un, *u;

    u = vector(1, n - 1);
    if (yp1 > 0.99e30)
        y2[1] = u[1] = 0.0;
    else {
        y2[1] = -0.5;
        u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
    for (i = 2; i <= n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
                / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n
                - 1]));
    }
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
    for (k = n - 1; k >= 1; k--)
        y2[k] = y2[k] * y2[k + 1] + u[k];
    free_vector(u, 1, n - 1);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
    void nrerror(char error_text[]);
    int klo, khi, k;
    double h, b, a;

    klo = 1;
    khi = n;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0)
        nrerror("Bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b
            - b) * y2a[khi]) * (h * h) / 6.0;
}

void splint_der(double xa[], double ya[], double y2a[], int n, double x,
        double *y) {
    void nrerror(char error_text[]);
    int klo, khi, k;
    double h, b, a;

    klo = 1;
    khi = n;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0)
        nrerror("Bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = (ya[khi] - ya[klo]) / h + (-(3.0 * a * a - 1.0) * y2a[klo] + (3.0 * b
            * b - 1.0) * y2a[khi]) * h / 6.0;
}

void lubksb(double **a, int n, int *indx, double b[]) {
    int i, ii = 0, ip, j;
    double sum;

    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

/*
 
 
 
 
 double Cabs(fcomplex a)
 {
 double c;
 
 c=sqrt(a.r*a.r+a.i*a.i);
 
 return c;
 }
 
 
 fcomplex Conj(fcomplex a)
 {
 fcomplex c;
 c.r=a.r;
 c.i=-a.i;
 return c;
 }
 
 
 fcomplex Cadd(fcomplex a, fcomplex b)
 {
 fcomplex c;
 c.r=a.r+b.r;
 c.i=a.i+b.i;
 return c;
 }
 
 fcomplex Csub(fcomplex a, fcomplex b)
 {
 fcomplex c;
 c.r=a.r-b.r;
 c.i=a.i-b.i;
 return c;
 }
 
 
 fcomplex Cmul(fcomplex a, fcomplex b)
 {
 fcomplex c;
 c.r=a.r*b.r-a.i*b.i;
 c.i=a.i*b.r+a.r*b.i;
 return c;
 }
 
 fcomplex RCmul(double a, fcomplex b)
 {
 fcomplex c;
 c.r=a*b.r;
 c.i=a*b.i;
 
 return c;
 }
 
 fcomplex Cdiv(fcomplex a, fcomplex b)
 {
 fcomplex c;
 double r,den;
 if (fabs(b.r) >= fabs(b.i)) {
 r=b.i/b.r;
 den=b.r+r*b.i;
 c.r=(a.r+r*a.i)/den;
 c.i=(a.i-r*a.r)/den;
 } else {
 r=b.r/b.i;
 den=b.i+r*b.r;
 c.r=(a.r*r+a.i)/den;
 c.i=(a.i*r-a.r)/den;
 }
 return c;
 }
 
 
 fcomplex Cexp(double a)
 {
 // c=e(i a) 
 fcomplex c;
 
 c.r=cos(a);
 c.i=sin(a);
 
 return c;
 
 }
 
 */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(double data[], int nn[], int ndim, int isign)

{
    int idim;
    unsigned long i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
    unsigned long ibit, k1, k2, n, nprev, nrem, ntot;
    double tempi, tempr;
    double theta, wi, wpi, wpr, wr, wtemp;

    for (ntot = 1, idim = 1; idim <= ndim; idim++)
        ntot *= nn[idim];
    nprev = 1;
    for (idim = ndim; idim >= 1; idim--) {
        n = nn[idim];
        nrem = ntot / (n * nprev);
        ip1 = nprev << 1;
        ip2 = ip1 * n;
        ip3 = ip2 * nrem;
        i2rev = 1;
        for (i2 = 1; i2 <= ip2; i2 += ip1) {
            if (i2 < i2rev) {
                for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2) {
                    for (i3 = i1; i3 <= ip3; i3 += ip2) {
                        i3rev = i2rev + i3 - i2;
                        SWAP(data[i3],data[i3rev]);
                        SWAP(data[i3+1],data[i3rev+1])
                        ;
                    }
                }
            }
            ibit = ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit) {
                i2rev -= ibit;
                ibit >>= 1;
            }
            i2rev += ibit;
        }
        ifp1 = ip1;
        while (ifp1 < ip2) {
            ifp2 = ifp1 << 1;
            theta = isign * 6.28318530717959 / (ifp2 / ip1);
            wtemp = sin(0.5 * theta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (i3 = 1; i3 <= ifp1; i3 += ip1) {
                for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2) {
                    for (i2 = i1; i2 <= ip3; i2 += ifp2) {
                        k1 = i2;
                        k2 = k1 + ifp1;
                        tempr = (double) wr * data[k2] - (double) wi * data[k2
                                + 1];
                        tempi = (double) wr * data[k2 + 1] + (double) wi
                                * data[k2];
                        data[k2] = data[k1] - tempr;
                        data[k2 + 1] = data[k1 + 1] - tempi;
                        data[k1] += tempr;
                        data[k1 + 1] += tempi;
                    }
                }
                wr = (wtemp = wr) * wpr - wi * wpi + wr;
                wi = wi * wpr + wtemp * wpi + wi;
            }
            ifp1 = ifp2;
        }
        nprev *= n;
    }
}

void opp_sort(double *v, int left, int right) {
    int i, last;

    if (left >= right)
        return;

    swap(v, left, (left + right) / 2);
    last = left;
    for (i = left + 1; i <= right; i++)
        if (v[i] > v[left])
            swap(v, ++last, i);

    swap(v, left, last);
    opp_sort(v, left, last - 1);
    opp_sort(v, last + 1, right);

}

void sort(double *v, int left, int right) {
    int i, last;

    if (left >= right)
        return;

    swap(v, left, (left + right) / 2);
    last = left;
    for (i = left + 1; i <= right; i++)
        if (v[i] < v[left])
            swap(v, ++last, i);

    swap(v, left, last);
    sort(v, left, last - 1);
    sort(v, last + 1, right);

}

void swap(double *v, int i, int j) {
    double temp;

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}

#define EPS 1.0e-6
#define JMAX 20 
#define JMAXP (JMAX+1)
#define K 7 

double qromb(double(*func)(double, int, double *, double *, double *),
        double a, double b, int n, double *a2, double *b2, double *c2) {
    double ss, dss;
    double s[JMAXP + 1], h[JMAXP + 1];
    int j;

    h[1] = 1.0;
    for (j = 1; j <= JMAX; j++) {
        s[j] = trapzd(func, a, b, j, n, a2, b2, c2);
        if (j >= K) {
            polint(&h[j - K], &s[j - K], K, 0.0, &ss, &dss);
            if (fabs(dss) < EPS * fabs(ss))
                return ss;
        }
        s[j + 1] = s[j];
        h[j + 1] = 0.25 * h[j];
    }
    nrerror("Too many steps in routine qromb");
    return ss;
    /*	return 0.0; */
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/*#define FUNC(x) ((*func)(x))*/

double trapzd(double(*func)(double, int, double *, double *, double *),
        double a, double b, int n, int N, double *a2, double *b2, double *c2) {
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1) {
        return (s = 0.5 * (b - a) * (func(a, N, a2, b2, c2) + func(b, N, a2,
                b2, c2)));
    } else {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 1; j <= it; j++, x += del)
            sum += func(x, N, a2, b2, c2);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}
/*#undef FUNC*/

void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs(x - xa[1]);
    c = vector(1, n);
    d = vector(1, n);
    for (i = 1; i <= n; i++) {
        if ((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            if ((den = ho - hp) == 0.0)
                nrerror("Error in routine polint");
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
    free_vector(d, 1, n);
    free_vector(c, 1, n);
}

double _Complex cspline_int(int n1, int n2, double *x, double _Complex *y) {
    const int N = n2 - n1 + 1;
    int i;
    double yrp1, yrpn, yip1, yipn, *yr, *yi, *yr2, *yi2, *xx, d, d1, dn;
    double _Complex J = 0.0;

    yr2 = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yi2 = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yr = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yi = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    xx = (double *) malloc((size_t) ((N + 1) * sizeof(double)));

    for (i = 1; i <= N; i++) {
        xx[i] = x[i + n1 - 1];
        yr[i] = creal(y[i + n1 - 1]);
        yi[i] = cimag(y[i + n1 - 1]);
    }

    d1 = xx[2] - xx[1];
    dn = xx[N] - xx[N - 1];
    yrp1 = (yr[2] - yr[1]) / d1;
    yrpn = (yr[N] - yr[N - 1]) / dn;
    yip1 = (yi[2] - yi[1]) / d1;
    yipn = (yi[N] - yi[N - 1]) / dn;

    spline(xx, yr, N, yrp1, yrpn, yr2);
    spline(xx, yi, N, yip1, yipn, yi2);

    for (i = 2; i <= N; i++) {
        d = xx[i] - xx[i - 1];
        J += (d * (yr[i - 1] + yr[i] - d * d / 12.0 * (yr2[i] + yr2[i - 1])));
        J += (I * d
                * (yi[i - 1] + yi[i] - d * d / 12.0 * (yi2[i] + yi2[i - 1])));
    }
    J *= 0.5;

    free(yr2);
    free(yi2);
    free(xx);
    free(yr);
    free(yi);

    return J;

}

double _Complex cxspline_int(int n1, int n2, double *x, double _Complex *y) {
    const int N = n2 - n1 + 1;
    int i;
    double yrp1, yrpn, yip1, yipn, *yr, *yi, *yr2, *yi2, *xx, d, d1, dn;
    double _Complex J = 0.0;

    yr2 = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yi2 = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yr = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yi = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    xx = (double *) malloc((size_t) ((N + 1) * sizeof(double)));

    for (i = 1; i <= N; i++) {
        xx[i] = x[i + n1 - 1];
        yr[i] = creal(y[i + n1 - 1]);
        yi[i] = cimag(y[i + n1 - 1]);
    }

    d1 = xx[2] - xx[1];
    dn = xx[N] - xx[N - 1];
    yrp1 = (yr[2] - yr[1]) / d1;
    yrpn = (yr[N] - yr[N - 1]) / dn;
    yip1 = (yi[2] - yi[1]) / d1;
    yipn = (yi[N] - yi[N - 1]) / dn;

    spline(xx, yr, N, yrp1, yrpn, yr2);
    spline(xx, yi, N, yip1, yipn, yi2);

    for (i = 2; i <= N; i++) {
        d = xx[i] - xx[i - 1];
        J += (d * (yr[i - 1] + yr[i] - d * d / 12.0 * (yr2[i] + yr2[i - 1])));
        J -= (I * d
                * (yi[i - 1] + yi[i] + d * d / 12.0 * (yi2[i] + yi2[i - 1])));
    }
    J *= 0.5;

    free(yr2);
    free(yi2);
    free(xx);
    free(yr);
    free(yi);
    return J;
}

double spline_int(int n1, int n2, double *x, double *y) {
    const int N = n2 - n1 + 1;
    int i;
    double yp1, ypn, *yy, *y2, *xx, d, J = 0.0;

    y2 = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    yy = (double *) malloc((size_t) ((N + 1) * sizeof(double)));
    xx = (double *) malloc((size_t) ((N + 1) * sizeof(double)));

    for (i = 1; i <= N; i++) {
        xx[i] = x[i + n1 - 1];
        yy[i] = y[i + n1 - 1];
    }

    yp1 = (yy[2] - yy[1]) / (xx[2] - xx[1]);
    ypn = (yy[N] - yy[N - 1]) / (xx[N] - xx[N - 1]);

    spline(xx, yy, N, yp1, ypn, y2);

    for (i = 2; i <= N; i++) {
        d = xx[i] - xx[i - 1];
        J += (d * (yy[i - 1] + yy[i] - d * d / 12.0 * (y2[i] + y2[i - 1])));
    }
    J *= 0.5;

    free(y2);
    free(xx);
    free(yy);

    return J;

}

double pow_di(double x, int n) {
    double pow = 1.0;

    if (n != 0) {
        if (n < 0) {
            if (x == 0)
                return (pow);
            n = -n;
            x = 1 / x;
        }
        for (;;) {
            if (n & 01)
                pow *= x;
            if (n >>= 1)
                x *= x;
            else
                break;
        }
    }
    return (pow);
}

int pow_i(int x, int n) {
    int pow = 1.0;

    if (n != 0) {
        if (n < 0) {
            if (x == 0)
                return (pow);
            n = -n;
            x = 1 / x;
        }
        for (;;) {
            if (n & 01)
                pow *= x;
            if (n >>= 1)
                x *= x;
            else
                break;
        }
    }
    return (pow);
}
