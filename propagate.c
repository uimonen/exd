#include "exactdiag.h"

void propagate(PPP_ints *pints, PPP_par *ppp){
  
  
  printf("propagation...\n");
  
  double test_norm;
  int max2 = ppp->max2;
  int size = ppp->size;
  int nb = ppp->nb;
  double tstep = ppp->tstep;
  int NT = ppp->NT; 
  int i, j, k, it, ispin, is;  
  double e0, t;
  double E0 = ppp->E0;
  
  double _Complex ss, sss, ww, wr2;
  double _Complex *w, *wc;
  double **dens;
  double _Complex *pi, **pf;
  
  FILE *fd;
  FILE *fd2;
  FILE *fd_up;
  FILE *fd_down;
  
  fd = fopen("td_densities.dat", "w");
  fprintf(fd,"# time site densities\n" );
  fclose(fd);
  
  fd = fopen("td_densities_up.dat", "w");
  fprintf(fd,"# time site densities\n" );
  fclose(fd);
  
  fd = fopen("td_densities_down.dat", "w");
  fprintf(fd,"# time site densities\n" );
  fclose(fd);
  
  w = cvector(1, max2);
  wc = cvector(1, max2);
  pi = cvector(1,max2);
  pf = cmatrix(1, max2, 1,max2);
  
  dens = matrix(1, nb, 1, 2);
  test_norm = 0.0;
  
  
  if(ppp->kick == 1){
    printf("\n\nPropagation with a kick....\n\n\n");
    kick(pints, ppp);
    
    
    for(i = 1; i <= max2; i++) {
      w[i] = 0.0;
      wc[i] = 0.0;
      if(i%2!=0)pi[i] = -pints->UH_t[1][i];
      else pi[i] = -pints->UH_t[1][i];
      for(j = 1; j <= max2; j++){
        if(j%2!=0)pf[i][j] = pints->UH_fin[j][i];
        else pf[i][j] = -pints->UH_fin[j][i];
      }
    }

  }
else{
 
  for(i = 1; i <= max2; i++) {
    w[i] = 0.0;
    wc[i] = 0.0;
    if(i%2!=0)pi[i] = -pints->UH[1][i];
    else pi[i] = -pints->UH[1][i];
    for(j = 1; j <= max2; j++){
      if(j%2!=0)pf[i][j] = pints->UH_fin[j][i];
      else pf[i][j] = -pints->UH_fin[j][i];
    }
  }

}
    
  fd2 = fopen("test.dat","w"); 
  for(i = 1; i <= max2; i++) {
    fprintf(fd2,"%14.14f \n",pi[i]);
  }
  fprintf(fd2,"\n");
  
  for(i = 1; i <= max2; i++) {
      for(j = 1; j <= max2; j++){
        fprintf(fd2,"%14.14f   ", pf[i][j]);
    }
    fprintf(fd2,"\n");
  }
  fclose(fd2);
  
  fd = fopen("basis_propagation.dat","w");
  
  /* Expand the initial state with final state wave functions*/
  for(i = 1; i <= max2; i++){
    for(j = 1; j <= max2; j++){
      w[i] +=  pi[j] * pf[j][i];
    }
    fprintf(fd, "%10.30f\n", w[i]);
    test_norm +=  w[i] * conj(w[i]);
   
  }
  printf("test_norm = %20.20f \n", test_norm);
  
  fclose(fd);
  
  
  /* an additional test??*/
  sss = 0.0;
  for(i = 1; i <= max2; i++){
    ss = 0.0;
    for(j = 1; j <= max2; j++){
      ss += w[j]*pf[j][i];
    }
    sss = sss + pow(ss,2);
  }
  printf("test = %20.20f\n", sss); /* should also be 1 with high accuracy */
  
  
  
  /* do the dynamics */
  e0 = pints->eh_fin[1];
  for(i = 1; i <= max2; i++){
    pints->eh[i] = pints->eh_fin[i] - e0;
  }
  
  for(it = 0; it <= NT; it++){
    t = it*tstep;
    for(j = 1; j <= max2; j++){
      wc[j] = w[j]*cexp(-I*(pints->eh_fin[j])*t);
    }
    printf("%d\n", it);
   
    for(i = 1; i <= nb; i++){
      dens[i][1] = 0.0;
      dens[i][2] = 0.0;
    }
    
    for(i = 1; i <= max2; i++){
       ww = 0.0;
      for (j = 1; j <= max2; j++) {
        ww += wc[j]*pf[i][j];
      }
      wr2 = ww*conj(ww);
     // printf("wr2 = %f ww = %f\n",wr2, ww);
      for (ispin = 1; ispin <= 2; ispin++) {
        for (is = 1; is <= nb; is++) {
          dens[is][ispin] += pints->ps[i][is][ispin]*wr2; 
        }
      }
    }
    
    
    double sum = 0.0;
    double sup = 0.0; double sdown = 0.0;
    fd = fopen("td_densities.dat","a");
    fd_up = fopen("td_densities_up.dat","a");
    fd_down = fopen("td_densities_down.dat","a");
    
    fprintf(fd,"%13.13f   ", t );
    fprintf(fd_up,"%13.13f   ", t );
    fprintf(fd_down,"%13.13f   ", t );
    
    for (is = 1; is <= nb; is++) {
      fprintf(fd, "%13.13f   ", dens[is][1]+dens[is][2]);
      fprintf(fd_up, "%13.13f   ", dens[is][1]);
      fprintf(fd_down, "%13.13f   ", dens[is][2]);
      
      sum += dens[is][1]+dens[is][2];
      sup += dens[is][1];
      sdown += dens[is][2];
    }
    fprintf(fd,"%13.13f \n", sum);
    fprintf(fd_up,"%13.13f \n", sup);
    fprintf(fd_down,"%13.13f \n", sdown);
    
    fclose(fd);
    fclose(fd_up);
    fclose(fd_down);  
    
  }
  
  free_cvector(w, 1, max2);    
  free_cvector(wc, 1, max2);
  free_cvector(pi, 1, max2);
  free_cmatrix(pf, 1, max2, 1, max2);
  free_matrix(dens, 1, nb, 1, 2);   
  
}


void kick(PPP_ints *pints, PPP_par *ppp){

  int size = ppp->size;
  double _Complex **kck = cmatrix(1,size,1, size);
  int i, j, k, l;
  double _Complex dum, sum;
  kick_matrix(pints, ppp, kck);
  
  for(i = 1; i <= size; i++) {
    for(j = 1; j <= size; j++){
      for( sum = 0.0, k = 1; k <= size; k++){
        for(l = 1;l <= size; l++){

          dum = conj(kck[i][k])*kck[j][l];
          sum += dum*pints->UH[k][l];
       
        }
      }
       pints->UH_t[i][j] = sum;
     
    }
  }
 
}


void kick_matrix(PPP_ints *pints, PPP_par *ppp, double _Complex **kck)
{
  int i, j, k;
  double **U, *d;
  double _Complex su;
  int nb = ppp->nb;
  int size = ppp->size;
  double E0 = ppp->E0;
  
  
  U=matrix(1, size, 1, size);
  d=(double *) malloc((size_t) ((size+1)*sizeof(double )));
  
  
  diag(size, pints->zzz, U, d);
  
  for(i=1;i<=size;i++)  
    for(j=1;j<=size;j++)
    { 
      for(su=0.0,k=1;k<=size;k++) su+=(U[k][i]*U[k][j]*cexp(-I*E0*d[k]));
      
      kck[i][j]=su;
   //   printf("RE(su) = %f IM(su) = %f\n", creal(su), cimag(su));
    }
  
  free_matrix(U,1,size,1,size);
  free(d);
  
  
}








