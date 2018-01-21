#include "exactdiag.h"


void dospctrum(PPP_ints *pints, PPP_par *ppp){
    
  printf("calculating the spectral function...\n");  
  
  int max2, max3, max4, maxm;
  int size = ppp->size;
  double lim = ppp->lim;
  double acc = ppp->acc;
  int nb = ppp->nb;
  int spin = 2;

  int n = (int)(2*lim*acc);
  
  FILE *fd;
           	
  pints->s = vector(0,n);
  pints->H = vector(0,n);
 
  max2 = ppp->max2;

  // One electron removed
  printf("calculate the eigenvector and eigenvalues for state one electron removed...\n");
        
  ppp->down = ppp->down - 1; // one electron removed

  confnumber(ppp, &(ppp->max1), &(ppp->max2));
  max3 = ppp->max2 =(ppp->max2) * (ppp->max1);
  ppp->size = max3;

  pints->H1 = vector(0,n);
  pints->UH_removed = cmatrix(1,max3,1,max3);
  pints->eh_removed = vector(1,max3);

  printf("\nConstruct basis for the case one eelctron removed...\n");
  makebasis(pints, ppp, 1);
  printbasis("basis_removed.dat", pints, ppp, 1);
  printf("\nConstruct Hamiltonian for the case one electron removed...\n");
  makehamiltonian(pints, ppp,1);
  printhamiltonian("hamiltonian_removed.dat",size, pints->ham);
  cdiag(max3,pints->ham,pints->UH_removed,pints->eh_removed);

  // One electron added
  printf("\nClculate the eigenvector and eigenvalues for state one electron added...\n");
       
  ppp->down = ppp->down+2; // back to to the intial state +1 and one added +1
 
  confnumber(ppp, &(ppp->max1), &(ppp->max2));
  max4 = ppp->max2 = (ppp->max2) * (ppp->max1);
  ppp->size = max4;  

  pints->H2 = vector(0,n);
  pints->UH_added = cmatrix(1,max4,1,max4);
  pints->eh_added = vector(1,max4);

  printf("\nConstruct basis for the case one electron added...\n");
  makebasis(pints, ppp, 2);      
  printbasis("basis_added.dat", pints, ppp, 2);
  printf("\nConstruct Hamiltonian for the case one electron added...\n");
  makehamiltonian(pints, ppp, 2);
  printhamiltonian("hamiltonian_added.dat",size, pints->ham);

  cdiag(max4,pints->ham,pints->UH_added,pints->eh_added);

  maxm = max(max2,max(max3,max4));
  pints->f =c3tensor(1, maxm, 1, nb, 1, spin);
  pints->g =c3tensor(1, maxm, 1, nb, 1, spin);

  printf("\nMake factors...\n");
  makefactors(max2, max3, max4, maxm, pints, ppp);
  printf("\nPrint factors...\n");
  writefactors(maxm, ppp->nb, pints->f, pints->g);
  printf("\nMake specral...\n");
  makespectral(maxm, pints, ppp);
  printf("\nMake norm...\n");
  makenorm(pints,ppp);
  printf("\nGreeF...\n");
  GreenF(maxm, pints, ppp);

  free_c3tensor(pints->f, 1, maxm, 1, nb, 1, spin);
  free_c3tensor(pints->g, 1, maxm, 1, nb, 1, spin);
  free_cmatrix(pints->UH_added, 1, max4, 1, max4);
  free_cmatrix(pints->UH_removed, 1, max3, 1, max3);
  free_vector(pints->eh_added, 1, max4);
  free_vector(pints->eh_removed, 1, max3);
  free_vector(pints->H1, 0, n);
  free_vector(pints->H2, 0, n);
  free_vector(pints->H, 0, n);
  free_vector(pints->s, 0, n);

}    


// determines the maximum of two integers a and b
int max(int a, int b){
  return ((a) > (b)) ? (a) : (b); 
}


//make factors f(k,i,is) and g(k,j,is')
void makefactors(int max2, int max3, int max4, int maxm, PPP_ints *pints, PPP_par *ppp){

  FILE *fd;
  int nb = ppp->nb;
  double _Complex ***f = pints->f;        
  double _Complex ***g = pints->g;        
  int ***ps = pints->ps;     // State vectors for N-particle system
  int ***ps_1 = pints->ps_1; // State vectors for N-1 particle system
  int ***ps_2 = pints->ps_2; // State vectors for N+1 particle system
  int part = ppp->part;     // Particle number
  double _Complex **UH = pints->UH;
  double _Complex **UH_added = pints->UH_added;
  double _Complex **UH_removed = pints->UH_removed;
  double *eh = pints->eh;
  double *eh_added = pints->eh_added;
  double *eh_removed = pints->eh_removed;
  int test;
  
  int i, j, k, is, i7, l, n;
	double sumf = 0.0, sumg = 0.0;        
  int sk = 0;
  
  for (k = 1; k <= maxm; k++){
    for (i = 1; i <= nb; i++){
      for (is = 1; is <= 2; is++){
        f[k][i][is] = 0.0;
        g[k][i][is] = 0.0;
      }
    }   
  }
  
  // f        
	for (k = 1; k <= max3; k++){
	  for (i = 1; i <= nb; i++){
	    for (is = 1; is <= 2; is++){
        //sumf = 0.0;
	      for(j = 1; j <= max2; j++){
          if (ps[j][i][is] == 1){
						sk  = 0;                      // commmute the annihilation operator to its place
						for(i7 = 1; i7 <= i-1; i7++) sk += ps[j][i7][is];
            if(is==2) sk = sk + pow_i(-1,part/2);
						ps[j][i][is] = 0;
            for (l = 1; l <= max3; l++){
              test = pow_i(-1, sk);
              for(n = 1; n <= nb; n++){
                if(ps_1[l][n][1] != ps[j][n][1]) test = 0;
                if(ps_1[l][n][2] != ps[j][n][2]) test = 0; 
              }
              f[k][i][is] = f[k][i][is] + UH_removed[k][l]*conj(UH[1][j])*test;
              //sumf += f[k][i][is] + UH_removed[k][l]*conj(UH[1][j])*test;
              //sumf +=  UH_removed[k][l]*conj(UH[1][j])*test;
              //sumf = f[k][i][is] + UH_removed[k][l]*conj(UH[1][j])*test;
            }
            ps[j][i][is] = 1;
					} 
				}
        //f[k][i][is] = sumf;
			}
		}
	}  
  

  //sk = 0;
  // g
	for (k = 1; k <= max4; k++){
	  for (i = 1; i <= nb; i++){
	    for (is = 1; is <= 2; is++){
        //sumg = 0.0;
	      for(j = 1; j <= max2; j++){
          if (ps[j][i][is] == 0){
						sk  = 0;                      // commmute the creation operator to its place
						for(i7 = 1; i7 <= i-1; i7++) sk += ps[j][i7][is];
            if(is==2) sk = sk + pow_i(-1,part/2);
						ps[j][i][is] = 1;
            for (l = 1; l <= max4; l++){
              test = pow_i(-1, sk);
              for(n = 1; n <= nb; n++){
                if(ps_2[l][n][1] != ps[j][n][1]) test = 0;
                if(ps_2[l][n][2] != ps[j][n][2]) test = 0; 
              }
               g[k][i][is] = g[k][i][is] + UH_added[k][l]*conj(UH[1][j])*test;
               //sumg = sumg + UH_added[k][l]*conj(UH[1][j])*test;
              //sumg = sumg + UH_added[k][l]*conj(UH[1][j])*test;
              // sumg = g[k][i][is] + UH_added[k][l]*conj(UH[1][j])*test;
            	
            }
            ps[j][i][is] = 0;
					} 
				}
        //g[k][i][is] = sumg;
       	 //  printf("%f\n", sumg);	
      }
		}
	}          
  
                                                             
  // These should match with the densities at initial state.
   fd = fopen("sums.dat", "w");
   for( i = 1; i <= nb; i++){
    for(is = 1; is <= 2; is++){
      sumf = 0.0;
      sumg = 0.0;
      for(k = 1; k <= maxm; k++){
        sumf += f[k][i][is]*conj(f[k][i][is]);
        sumg += g[k][i][is]*conj(g[k][i][is]);
      }
      fprintf(fd, "%d %d %f %f %f\n", i, is, sumf, sumg, sumf + sumg);
    }
  }
  fclose(fd);

}



// Writes the factors f and g to files f.dat and g.dat
void writefactors(int maxm, int nb, double _Complex ***f, double _Complex ***g){
	
	int is, j, i;
	FILE *fd;

	fd = fopen("f.dat", "w");
	for (is = 1; is <= 2; is++){
		for(j = 1; j <= maxm; j++){
			for(i = 1; i <= nb; i++){
				fprintf(fd,"%f %f", creal(f[j][i][is]), cimag(f[j][i][is]));				
			}
			fprintf(fd, "\n");
		}
	}                              
	fclose(fd);

	fd = fopen("g.dat", "w");
	for (is = 1; is <= 2; is++){
		for(j = 1; j <= maxm; j++){
			for(i = 1; i <= nb; i++){
				fprintf(fd,"%f %f", creal(g[j][i][is]), cimag(g[j][i][is]) );				
			}
			fprintf(fd, "\n");
		}
	}                              
	fclose(fd);
                              
}


// calculates the spectral function A(w)
void makespectral(int maxm, PPP_ints *pints, PPP_par *ppp){

  double lim = ppp->lim;
  double acc = ppp->acc;
  int nb = ppp->nb;

  int n = (int)(2*lim*acc);
  int l, i, k, j, is;
  double A, A1, A2, d1, d2, w;

  double *H=pints->H;
  double *H1=pints->H1;
  double *H2=pints->H2;
  double *s=pints->s;

  double _Complex ***f = pints->f;        
  double _Complex ***g = pints->g;        

  double *eh = pints->eh;
  double *eh_removed = pints->eh_removed;
  double *eh_added = pints->eh_added;

  double pi=4.0*atan(1.0);
  double c=0.01/pi;

  for(l=0; l<= n; l++){
    w = - lim + l / acc;
    A1 = 0.0;
    A2 = 0.0;
    d1 = 0.0;
    d2 = 0.0;

    for(i = 1; i <= nb; i++){
      for(is = 1; is <= 2; is++){
        for(k = 1; k <= maxm; k++){
          d1 = c/((w-eh[1]+eh_removed[k])*(w-eh[1]+eh_removed[k])  + 0.01*0.01);
          d2 = c/((w+eh[1]-eh_added[k])*(w+eh[1]-eh_added[k]) + 0.01*0.01);
          A1 = A1 + 2*pi*f[k][i][is]*conj(f[k][i][is])*d1; 
          A2 = A2 + 2*pi*g[k][i][is]*conj(g[k][i][is])*d2;	  
        }
      }
    }
    s[l] = w;
    H1[l] = A1;
    H2[l] = A2;
    H[l] = A1 + A2;
  }

  pints->H = H;                              
  pints->s = s;                              
  pints->H1 = H1;                              
  pints->H2 = H2;                              
                              
}

// calculates the spectral function A(w)
void GreenF(int maxm, PPP_ints *pints, PPP_par *ppp){

  double lim = ppp->lim;
  double acc = ppp->acc;
  int nb = ppp->nb;

  int n = (int)(2*lim*acc);
  int l, i, k, j, is;
  double _Complex A, A1, A2, d1, d2;
  double w;

  double *H  = pints->H;
  double *H1 = pints->H1;
  double *H2 = pints->H2;
  double *s  = pints->s;

  double _Complex ***f = pints->f;        
  double _Complex ***g = pints->g;        

  double *eh = pints->eh;
  double *eh_removed = pints->eh_removed;
  double *eh_added = pints->eh_added;

  double pi=4.0*atan(1.0);
  double eta = 1e-9;
  double c=eta/pi;

  FILE *fd;
  fd = fopen("GreenF.dat","w");
  
  for(l=0; l<= n; l++){
    w = - lim + l / acc;
    A1 = 0.0;
    A2 = 0.0;
    d1 = 0.0;
    d2 = 0.0;

    for(i = 1; i <= nb; i++){
      for(is = 1; is <= 2; is++){
        for(k = 1; k <= maxm; k++){
          d1 = 1.0/((w-eh[1]+eh_removed[k])  + I*eta);
          d2 = 1.0/((w+eh[1]-eh_added[k]) + I*eta);
          A1 = A1 + f[k][i][is]*conj(f[k][i][is])*d1; 
          A2 = A2 + g[k][i][is]*conj(g[k][i][is])*d2;	  
        }
      }
    }
    fprintf(fd, "%f %f %f \n", w, creal(A1+A2), cimag(A1+A2));
  }
  fclose(fd);
}



//normalization for spectral function, if not correctly normalized.
void makenorm(PPP_ints *pints, PPP_par *ppp){

  int nb = ppp->nb;
  double lim = ppp->lim;
  double acc = ppp->acc;

  int n = (int)(2*lim*acc);
  int l, i, k, j, is;
  double A  = 0.0;
  double *s=pints->s;
  double *H=pints->H;
  double *H1=pints->H1;
  double *H2=pints->H2;
  double In = 0.0;
  double e = 0.0;

  double pi=4.0*atan(1.0);
  double c=0.2/pi;

  FILE *fd;  
  FILE *fd2;  
  fd = fopen("testing.dat","w");  
  fd2 = fopen("spectrum.dat","w");  

  for(k = 0; k <= n-1; k++){
    A = ( H[k] + H[k+1]) / (2*acc);
    In += A;
  }                                
  
  e = ( (double) (nb) ) / In;
  
  fprintf(fd, "%f %f %d", In/pi, e, nb); 

  for( i = 0; i <= n; i++){
    H[i] = e*H[i];
    H1[i] = e*H1[i];
    H2[i] = e*H2[i];
    fprintf(fd2, "%f %f %f %f \n", s[i], H1[i], H2[i], H[i]);
  }                            

  fclose(fd);                            
  fclose(fd2);                            

}



//-------------------------------------------------------
/*void Construct_G(int TSTEPS, int NL, int NTOT, int N_S, double mu, double beta, 
double *eigval, double _Complex **t_orbs, double _Complex **G)
{
  int p,q,i,j,k;
  FILE* fd;

  double _Complex n;

  // Initialize the Green function matrix to zero content
	for(i = 1; i <= N_S; i++){
		for(q = 1; q <= TSTEPS; q++){
			for(p = 1; p <= TSTEPS ;p++){
				G[ i+N_S*(i-1) ][ q+TSTEPS*(p-1) ] = 0.0;
			}
		}
	}

	for(q = 1; q <= TSTEPS; q++){
		for(p = 1; p <= q; p++){
			for(i = 1;i <= N_S; i++)
			{
				for(k = 1; k <= NTOT; k++)
				{
		      // LESSER PART
			     G[i+N_S*(i-1)][q+TSTEPS*(p-1)]+=I*t_orbs[k+NTOT*(i-1)][p]*conj(t_orbs[k+NTOT*(i-1)][q]);
					//G[i+N_S*(i-1)][q+TSTEPS*(p-1)] += I*fdist(eigval[k],beta,mu)*t_orbs[k+NTOT*(i-1)][p]*conj(t_orbs[k+NTOT*(i-1)][q]);
				}
			}
		}
	}

	for(q = 1; q <= TSTEPS; q++){
		for(p = 1 ;p < q; p++){ 
			for(i = 1; i <= N_S; i++){
				for(k = 1; k <= NTOT; k++)
				{
					// GREATER PART
					 G[i+N_S*(i-1)][p+TSTEPS*(q-1)]+=I*t_orbs[k+NTOT*(i-1)][q]*conj(t_orbs[k+NTOT*(i-1)][p]);
					//G[i+N_S*(i-1)][p+TSTEPS*(q-1)] += I*(fdist(eigval[k],beta,mu)-1.0)*t_orbs[k+NTOT*(i-1)][q]*conj(t_orbs[k+NTOT*(i-1)][p]);
				}
			}
		}
	}

  fd = fopen("GF","w");
	for(i = 1; i <= N_S; i++){
		for(p = 1; p <= TSTEPS; p++)
		{
			for(q = 1; q <= TSTEPS; q++)
			{
				fprintf(fd,"%d %d %d %10.10f %10.10f\n",i,p,q,
					creal(G[i+N_S*(i-1)][q+TSTEPS*(p-1)]), cimag(G[i+N_S*(i-1)][q+TSTEPS*(p-1)]));
			}
			fprintf(fd,"\n");
		}
	}
	fclose(fd);

  n = 0.0;
  for(i = 1; i <= N_S; i++){
    n += -I*G[i+N_S*(i-1)][1];
	}
  
	printf("n = %10.10f + I * %10.10f\n", creal(n), cimag(n));

}
*/
//-------------------------------------------------------

/*
* FOR A GIVEN CENTRE-OF-TIME COORDINATE T, PLOT THE TRACE OF THE SPECTRAL FUNCTION 
*
*
*/
void print_Atrace(char *fn, int T, int _TT, double _D, int nb, double _Complex **G, double mu)
{
  int N=9000-1; /* THE NUMBER OF POINTS USED IN THE FOURIER TRANSFORM */	
  int ind;
  int p, q, T1, T2, i,  N1;
  double _Complex *xx, *y;
  double _Complex su, **A;
  FILE *fname;
  double sign;
  double pi2 = 2.0*M_PI;

  /* FOR A GIVEN CENTER OF TIME VARIABLE T, THE TIME-COORDINATES t and t'  */
  /* ARE IN THE RANGE SPECIFIED BY T1 AND T2 */
  T1=MAX(1,2*T-_TT);
  T2=MIN(2*T-1,_TT);
  N1=T2-T1+1;  /* THE NUMBER OF POINTS ON THE RELATIVE-TIME GRID */

  xx=(double _Complex*) malloc((size_t) ((N+4)*sizeof(double _Complex)));
  y=(double _Complex *) malloc((size_t) ((2*N1+5)*sizeof(double _Complex)));

  A=cmatrix(1, N1, 1, nb);

  /* 1. Calculate A(t_p-t_q) = i[ G^>(t_p, t_q) -G^<(t_p,t_q) ] */
  for(i=1;i<=nb;i++)
  {
    ind=i+nb*(i-1);
    for(p=T1;p<=T2;p++)
    {
      q=2*T-p;
      /* IF t_p > t_q THEN G^>_ij(t_p,t_q) = G_ij(p,q) AND G^<=-G^*_ji(q,p) */
      /* IF t_p < t_q THEN G^>_ij(t_p,t_q) = -G^*_ji(q,p) AND G^<=G_ij(p,q) */
      /* IF t_p = t_q THEN G^>_ii=-i+G_ii^<  and A_ii(0) = 1.0 */
      if(p>q) A[p+1-T1][i]=I*(G[ind][q+_TT*(p-1)]+conj(G[ind][p+_TT*(q-1)]));
      else if(p==q) A[p+1-T1][i]=1.0;
      else A[p+1-T1][i]=-I*(G[ind][q+_TT*(p-1)]+conj(G[ind][p+_TT*(q-1)]));
    }
  }

  for(p=1;p<=N1;p++)
  {
    /* Calculate y(t) = Trace A(t) */
    for(su=0.0,i=1;i<=nb;i++) su+=A[p][i];
    y[p-1]=su;
  }

  /* FOURIER TRANFORM y */
  cft(N1, N, y, xx);

  free_cmatrix(A, 1, N1, 1, nb);

  /* PLOT THE SPECTRAL FUNCTION. A SHIFT OF mu (CHEMICAL POTENTIAL) IS */
  /* NECESSARY                                                         */
  /* THE FT ARRAY IS REVERSED, DUE TO THE OPPOSITE SIGN IN THE DEFINITION */
  /* OF THE FOURIER TRANSFORM IN THE INTEL MKL LIBRARY                    */
  fname=fopen(fn, "w");

  for(p=N/2;p<=N-1;p++)
  {
    sign = pow(-1.0,p);
    fprintf(fname,"%f %e %e %e\n", pi2*(p-N)/(N*2.0*_D)+mu,2.0*_D*cabs(xx[N-p-1]), 2.0*sign*_D*creal(xx[N-p-1]), 2.0*sign*_D*cimag(xx[N-p-1]));
  }
  for(p=0;p<=N/2-1;p++)
  {
    sign = pow(-1.0,p+1);
    fprintf(fname,"%f %e %e %e\n", pi2*p/(N*2.0*_D)+mu, 2.0*_D*cabs(xx[N-p-1]), 2.0*sign*_D*creal(xx[N-p-1]), 2.0*sign*_D*cimag(xx[N-p-1]) );
  }

  fclose(fname);

  free(xx);
  free(y);
}

//-------------------------------------------------------

void cft(int n1, int N, double _Complex *y, double _Complex *xx)
{
  int i;
  long status;
  DFTI_DESCRIPTOR *my_desk_handle1;

  for(i = 0; i <= N/2-n1/2-1; i++) 			xx[i] = 0.0;
  for(i = N/2-n1/2; i <= N/2+n1/2; i++) xx[i] = y[i-N/2+n1/2];
  for(i = N/2+n1/2+1; i <= N-1; i++) 		xx[i] = 0.0;

  status = DftiCreateDescriptor(&my_desk_handle1, DFTI_DOUBLE, DFTI_COMPLEX,1,N);
  if(status != DFTI_NO_ERROR)
  {
    printf("status %ld\n",status);
    nrerror("error in create descriptor");
  }
  status = DftiCommitDescriptor(my_desk_handle1);
  status = DftiComputeForward(my_desk_handle1,xx);
  status = DftiFreeDescriptor(&my_desk_handle1);
}

//-------------------------------------------------------
void store_torbs(int p, int NL, int NTOT, int N_S, double _Complex **orbs, double _Complex **t_orbs)
{
  int ind;
  int i,k;

  for(i=1;i<=N_S;i++)
  {
    ind=NL+i;
    for(k=1;k<=NTOT;k++)
    {
      t_orbs[k+NTOT*(i-1)][p] = orbs[ind][k];
    }
  }
}
//-------------------------------------------------------




