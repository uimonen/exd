#include "exactdiag.h"
//#include </opt/intel/Compiler/11.0/059/Frameworks/mkl/include/mkl_lapack.h>
//#include </opt/intel/Compiler/11.0/059/Frameworks/mkl/include/mkl_blas.h>
#include "num.h"


void diag(int nb, double **M, double **U, double *e)
{
    /* DIAGONALIZE A REAL SYMMETRIC MATRIX A OF DIMENSION nb, RETURNING */
    /* THE EIGENVECTORS IN U AND THE EIGENVALUES IN e                   */
    /* USING THE LAPACK ROUTINE dsyevx_                                 */
    char jobz='V', range='A', uplo='u', cmach='s';
    int n,lda,ldz;
    int i,j,  lwork=8*nb, il=1, iu=nb, m, *iwork, *ifail, info;
    double *f, *A, vl=0.0, vu=1.0e15, abstol, *q, *work;
    
    f=(double *) malloc((size_t) ((nb+1)*sizeof(double)));
    A=(double *) malloc((size_t) ((nb*nb+1)*sizeof(double)));
    q=(double *) malloc((size_t) ((nb*nb+1)*sizeof(double)));
    work=(double *) malloc((size_t) ((8*nb+1)*sizeof(double)));
    iwork=(int *) malloc((size_t) ((8*nb+1)*sizeof(int)));
    ifail=(int *) malloc((size_t) ((nb+1)*sizeof(int)));
    
    n=lda=ldz=nb;
    for(i=1;i<=nb;i++) for(j=i;j<=nb;j++) A[j+nb*(i-1)]=A[i+nb*(j-1)]=M[i][j];
    
    abstol=2.0*dlamch(&(cmach));
    dsyevx(&(jobz), &(range), &(uplo), &(n), &(A[1]), &(lda), &(vl), &(vu), 
           &(il), &(iu), &(abstol), &(m), &(e[1]), &(q[1]), &(ldz), &(work[1]),
           &(lwork), &(iwork[1]), &(ifail[1]), &(info));
    
    
    if(info)
    {
        if(info<0) printf("The %dth argument had an illegal value\n", -info);
        else
        {
            printf("**************\n%d eigenvalues failed to converge\n", info);
            for(i=1;i<=info;i++) printf("%d\n", ifail[i]);
            printf("**************\n");
        }
        
    }
    
    for(i=1;i<=nb;i++) for(j=1;j<=nb;j++) U[i][j]=q[j+nb*(i-1)];
    
    free(f);
    free(A);
    free(q);
    free(work);
    free(iwork);
    free(ifail);
    
}

/**************************************************************************/

void cdiag(int N, double _Complex **A, double _Complex **u2, double *d)
{
    /* Given a hermitian complex matrix u, the eigenvalues (d) and   */
    /* eigenvectors (q)  are caluclated                              */
    
    int ind, ind2, i, j, n=N, m, lda, ldz, info, lwork, il=1, iu=N, *iwork, *ifail;
    double *work, *rwork, vl=0.0, vu=1e20, abstol, *q, *u;
    char jobz, uplo, range, cmach='s';
    
    
    lwork=2*N;
    u=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    q=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    work=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    rwork=(double *) malloc((size_t) ((7*N*N+1)*sizeof(double)));
    iwork=(int *) malloc((size_t) ((5*N*N+1)*sizeof(int)));
    ifail=(int *) malloc((size_t) ((N+1)*sizeof(int)));
    
    for(i=1;i<=N;i++) 
        for(j=1;j<=N;j++)
        {
           // ind=(j+N*(i-1))<<1;
            //ind2=(i+N*(j-1))<<1;
            ind = (j+N*(i-1))*2;
            u[ind-1]=creal(A[i][j]);
            u[ind]=cimag(A[i][j]);
           // u[ind2]=-cimag(A[i][j]);
        }
    
    jobz='v';
    uplo='l';  // was u
    range='a';
    
    n=m=lda=ldz=N;
    abstol=dlamch(&(cmach)); /* the absolute error tolerance for the eigenvals */
    
    zheevx(&jobz, &range, &uplo, &n, &(u[1]), &lda, &vl, &vu, &il, &iu, &abstol,
           &m, &(d[1]), &(q[1]), &ldz, &(work[1]), &lwork, &(rwork[1]),
           &(iwork[1]), &(ifail[1]), &info);
 
    
    if(info) printf("INFO %d\n", info);
    
    for(i=1;i<=N;i++)
        for(j=1;j<=N;j++)
        {
            //ind=(j+N*(i-1))<<1;
          ind = (j+N*(i-1))*2;
            //u2[i][j]=q[ind-1]+I*q[ind];
            u2[i][j]=q[ind-1]+I*q[ind];
        }
    
    free(u);
    free(q);
    free(work);
    free(rwork);
    free(iwork);
    free(ifail);
    
    return;
    
}
/**************************************************************************/

void cdiag_old(int N, double _Complex **A, double _Complex **u2, double *d)
{
    /* Given a hermitian complex matrix u, the eigenvalues (d) and   */
    /* eigenvectors (q)  are caluclated                              */
    
    int ind, ind2, i, j, n=N, m, lda, ldz, info, lwork, il=1, iu=N, *iwork, *ifail;
    double *work, *rwork, vl=0.0, vu=1e15, abstol, *q, *u;
    char jobz, uplo, range, cmach='s';
    
    
    lwork=2*N;
    u=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    q=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    work=(double *) malloc((size_t) ((4*N*N+1)*sizeof(double)));
    rwork=(double *) malloc((size_t) ((7*N*N+1)*sizeof(double)));
    iwork=(int *) malloc((size_t) ((5*N*N+1)*sizeof(int)));
    ifail=(int *) malloc((size_t) ((N+1)*sizeof(int)));
    
    for(i=1;i<=N;i++) 
        for(j=1;j<=i;j++)
        {
            ind=(j+N*(i-1))<<1;
            ind2=(i+N*(j-1))<<1;
            u[ind-1]=u[ind2-1]=creal(A[i][j]);
            u[ind]=cimag(A[i][j]);
            u[ind2]=-cimag(A[i][j]);
        }
    
    jobz='v';
    uplo='l';  // was u
    range='a';
    
    n=m=lda=ldz=N;
    abstol=dlamch(&(cmach)); /* the absolute error tolerance for the eigenvals */
    
    zheevx(&jobz, &range, &uplo, &n, &(u[1]), &lda, &vl, &vu, &il, &iu, &abstol,
           &m, &(d[1]), &(q[1]), &ldz, &(work[1]), &lwork, &(rwork[1]),
           &(iwork[1]), &(ifail[1]), &info);
 
    
    if(info) printf("INFO %d\n", info);
    
    for(i=1;i<=N;i++)
        for(j=1;j<=N;j++)
        {
            ind=(j+N*(i-1))<<1;
            //u2[i][j]=q[ind-1]+I*q[ind];
            u2[i][j]=q[ind-1]+I*q[ind];
        }
    
    free(u);
    free(q);
    free(work);
    free(rwork);
    free(iwork);
    free(ifail);
    
    return;
    
}

/**************************************************************************/

void gdiag(int nb, double **M, double **rho, double **U, double *e)
{
    /* DIAGONALIZE A REAL SYMMETRIC MATRIX A OF DIMENSION nb, RETURNING */
    /* THE EIGENVECTORS IN U AND THE EIGENVALUES IN e                   */
    /* USING THE LAPACK ROUTINE dsyevx_                                 */
    char jobz='V', range='A', uplo='u', cmach='s';
    int n,lda,ldb,ldz, itype=1;
    int i,j, lwork=8*nb,  il=1, iu=nb, m, *iwork, *ifail, info;
    double *f, *A, *B, vl=0.0, vu=1.0e15, abstol, *q, *work;
    
    f=(double *) malloc((size_t) ((nb+1)*sizeof(double)));
    A=(double *) malloc((size_t) ((nb*nb+1)*sizeof(double)));
    B=(double *) malloc((size_t) ((nb*nb+1)*sizeof(double)));
    q=(double *) malloc((size_t) ((nb*nb+1)*sizeof(double)));
    work=(double *) malloc((size_t) ((8*nb+1)*sizeof(double)));
    iwork=(int *) malloc((size_t) ((8*nb+1)*sizeof(int)));
    ifail=(int *) malloc((size_t) ((nb+1)*sizeof(int)));
    
    n=lda=ldb=ldz=nb;
    for(i=1;i<=nb;i++) for(j=i;j<=nb;j++) A[j+nb*(i-1)]=A[i+nb*(j-1)]=M[i][j];
    for(i=1;i<=nb;i++) for(j=i;j<=nb;j++) B[j+nb*(i-1)]=B[i+nb*(j-1)]=rho[i][j];
    
    abstol=2.0*dlamch(&(cmach));
    
    dsygvx_(&(itype),&(jobz),&(range), &(uplo),&(n), &(A[1]),&(lda),&(B[1]),
            &(ldb),&(vl),&(vu),&(il),&(iu),&(abstol),&(m),&(e[1]),&(q[1]),&(ldz),
            &(work[1]), &(lwork),&(iwork[1]),&(ifail[1]),&(info));
    
    if(info)
    {
        if(info<0) printf("The %dth argument had an illegal value\n", -info);
        else if(info<=nb)
        {
            printf("**************\n%d eigenvalues failed to converge\n", info);
            for(i=1;i<=info;i++) printf("%d\n", ifail[i]);
            printf("**************\n");
        }
        else if(info>nb) printf("info %d \n", info);
        
    }
    
    for(i=1;i<=nb;i++) for(j=1;j<=nb;j++) U[i][j]=q[j+nb*(i-1)];
    
    free(f);
    free(A);
    free(B);
    free(q);
    free(work);
    free(iwork);
    free(ifail);
    
}


void ctranpose(int N, double _Complex **A, double _Complex **B){

  int i,j;
  
  for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        B[j][i] = A[j][i];
      }
  }
  
}
