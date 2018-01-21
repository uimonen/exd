#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#include "num.h"

#define INPUT "input.in"
#define OUTPUT "output.out"
#define cdouble double _Complex

/********STRUCTURES*****************/

typedef struct PPP_par{
    int kick;
    int nb;                     /*  Number of sites.           */
    int up;                     /*  Number of up electrons.    */
    int part;                   /*  Total number of particles. */
    int down;                   /*  Number of down electrons.  */
    double U;                   /*  On-site electron-electron strength sometimes called Hubbard U. */
    double E0;                  /*  Electric field strenght     */
    double as;                  /*      */
    double al;                  /*      */
    double a_c;                 /* On-site potentia.l            */
    double b_c;                 /* Hopping between the sites.    */        
    int int_type;               /* Interaction type 1 = long range, 0 = Hubbard. */
    int size;                   /* Size of the Hilbert space.     */
    int max1;                   /* Number of configurations for the up electrons */
    int max2;                   /* Number of configurations for the down electrons. */
    double alpha;               /* Phase factor of the wavefunction (when the magnetic field is applied).     */
    double tstep;               /* Time-step length.     */
    int NT;                     /* Number of time-steps. */
    int chemical_structure;     /* Used ppp-model. */
    int do_spectrum;            /* Parameter to determine is the spectral function calculated. */
    double acc;                 /* Accuracy of the specral function caluculation.  */
    double lim;                 /* Limis of the plot for the spectral function.   */
} PPP_par;


typedef struct PPP_ints{
    int nsites;                     /*   Number of basis functions                  */
    double **psi;                   /*   Molecular orbital expansion coefficients   */
    double **hmat;                  /*   Kinetic energy matrix                      */
    int **con;                      /*   Connection matrix                          */
    double **zmat;                  /*   Dipole matrix                              */
    double **v_ijkl;                /*   Two-electron integrals                     */
    double **dist;          	      /*   Distance matrix                            */
    double **basis;                 /*               */
    int ***ps;                      /*   State vectros of the Hilbert space for one e removed  */
    int ***ps_1;                    /*   State vectors the Hilbert space for one e added       */
    int ***ps_2;                    /*                    */
    double **zzz;                   /*   Dipole matrix in normal basis              */
    double _Complex **ham;          /*   Hamiltonian for the initial state          */
    double _Complex **ham_fin;      /*   Hamiltonian for the final state            */
    double _Complex **UH;           /*   Wave functions for the initial state       */
    double _Complex **UH_t;
    double _Complex **UH_fin;       /*   Wave functions for the final state         */
    double _Complex **UH_fin_t;     /*   Wave functions for the final state         */
    double _Complex **UH_added;     /*   Wave functions for the one electron added to the initial state           */
    double _Complex **UH_removed;   /*   Wave functions for the one electron removed from the initial state       */
    double *eh;                     /*   Eigenvalues of the intial state            */
    double *eh_fin;                 /*   Eigenvalues of the final state             */
    double *eh_added;               /*   Eigenvalues of the state one electron added      */
    double *eh_removed;             /*   Eigenvalues of the state one electron removed    */
    double **ext;
    double **vpot;                  /*   Spin dependent potential for the intial state    */
    double **vpot_fin;              /*   Spin dependent potential for the final state     */
    double **dens;                  /*   Site densities    */
    double _Complex ***g;
    double _Complex ***f;
    double *H;
    double *H1;
    double *H2;    
    double *s;    
} PPP_ints;


/*************************************/


/* ioprogs.c */
void read_input(char *input_file, PPP_ints *pints, PPP_par *ppp);
void print_parameters(char *output_file, PPP_par *ppp, PPP_ints *pints);
void printbasis(char *filename, PPP_ints *pints, PPP_par *ppp, int i);
void printhamiltonian(char *filename, int max2, double _Complex **ham);
void print_matrices(PPP_ints *pints);
void printdensities(char *filename, double _Complex **UH, PPP_ints *pints, PPP_par *ppp);
void printwf(char *filename, double _Complex **UH, PPP_ints *pints, PPP_par *ppp);

/* basis.c */
void construct_ham(PPP_ints *pints, PPP_par *ppp);
double tmn(PPP_par *par, double **V, int ns, int i, int j);
int connections(PPP_par *par, int ns, int i, int j);
double zmn(PPP_par *par,int i, int j);
double vijkl(PPP_par *par,int i, int j);
void zeroham_mat(PPP_ints *pints);
void makebasis(PPP_ints *pints, PPP_par *ppp, int i);
void makehamiltonian(PPP_ints *pints, PPP_par *ppp, int i);
void confnumber(PPP_par *ppp, int *max1, int *max2);

/* propagate.c*/
void propagate(PPP_ints *pints, PPP_par *ppp);
void kick_matrix(PPP_ints *pints, PPP_par *ppp, double _Complex **kck);
void kick(PPP_ints *pints, PPP_par *ppp);

/* outputs.c*/
int max(int a, int b);
void writefactors(int maxm, int nb, double _Complex ***f, double _Complex ***g);
void makefactors(int max2, int max3, int max4, int maxm, PPP_ints *pints, PPP_par *ppp);
void makespectral(int maxm, PPP_ints *pints, PPP_par *ppp);
void makenorm(PPP_ints *pints, PPP_par *ppp);
void cft(int n1, int N, double _Complex *y, double _Complex *xx);
void GreenF(int maxm, PPP_ints *pints, PPP_par *ppp);

/* diag.c */
void diag(int nb, double **M, double **U, double *e);
void cdiag(int N, double _Complex **A, double _Complex **u2, double *d);
void cdiag_old(int N, double _Complex **A, double _Complex **u2, double *d);
void gdiag(int nb, double **M, double **rho, double **U, double *e);


/* LAPACK ROUTINES */
void dgeev(char *jobvl, char *jobvr, int *N, double *A, int *LDA, double *WR,
           double *WI, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
void dsyevx(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA,
            double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W,
            double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *IFAIL, int *INFO );
double dlamch(char *cmach);
void dsygvx_(int *itype,char *jobz,char *range,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,
             double *vl,double *vu,int *il,int *iu,double *abstol,int *m,double *w,double *z,int *ldz,
             double *work,int *lwork,int *iwork,int *ifail,int *info);

void zgetrf(int *m, int *n, double *A, int *LDA, int *ipiv, int *info);
void dgetrf(int *m, int *n, double *A, int *LDA, int *ipiv, int *info);
void zheev(char *jobz, char *uplo, int *n, double *a, int *lda, double *x,
           double *work, int *lwork, double *rwork, int *info);
void zheevx(char *jobz, char *range, char *uplo, int *n, double *a, int *lda,
            double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
            double *w, double *z, int *ldz, double *work, int *lwork,
            double *rwork, int *iwork, int *ifail, int *info);

