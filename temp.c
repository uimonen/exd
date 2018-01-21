/*
 * Propagation with  a kick
 *
 */

void propagate_kick(PPP_ints *pints, PPP_par *ppp){
  
  
  printf("propagation with the kick...\n");
  
  double test_norm;
  int max2 = ppp->max2;
  int nb = ppp->nb;
  double tstep = ppp->tstep;
  int NT = ppp->NT; 
  int i, j, k, it, ispin, is;  
  double e0, t;
  int size = ppp->size;
  double E0 = ppp->E0;
  
  
  double supa, sdowna;
  
  
  double _Complex ss, sss, ww, wr2;
  double _Complex *w, *wc;
  double **dens;
  double _Complex *pi, **pf;
  
  
  pints->UH_t = cmatrix(1,size,1,size);
  pints->eh_t = vector(1,size);
  pints->UH_fin_t = cmatrix(1,size,1,size);
  pints->eh_fin_t = vector(1,size);
  
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
  
  
  int first = 1;
  double first_step = 1E-12;  
  
  
  /* do the dynamics */
  for(it = 0; it <= NT; it++){
    
    for(i = 1; i<= nb; i++){
      for(j = 1; j<=2; j++){
        pints->vpot[i][j] = 0.0;
        pints->vpot_fin[i][j] = 0.0;
        if(i == 1 && first == 1) pints->vpot[i][j] = E0/first_step; 
      }
    }      
    
    if(first == 1) t = it*first_step;
    else t = it*tstep;
    
    makebasis(pints, ppp, 3);
    makehamiltonian_t(pints, ppp, 3);
    
    cdiag(size,pints->ham,pints->UH_t,pints->eh_t);
    
    if(first == 1){
      for(i = 1; i <= max2; i++) {
        w[i] = 0.0;
        if(i%2!=0)pi[i] = -pints->UH[1][i];
        else pi[i] = pints->UH[1][i];
        for(j = 1; j <= max2; j++){
          if(j%2!=0)pf[i][j] = pints->UH_t[j][i];
          else pf[i][j] = -pints->UH_t[j][i];
        }
      }
      first = 0;
    }
    else{
      for(i = 1; i <= max2; i++) {
        for(j = 1; j <= max2; j++){
          if(j%2!=0)pf[i][j] = pints->UH_t[j][i];
          else pf[i][j] = -pints->UH_t[j][i];
        }
      }
    }
    
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
    
    //    for(i = 1; i <= max2; i++){
    //     for(j = 1; j <= max2; j++){
    //       w[i] = pf[j][i]*cexp(-I*(pints->eh_t[j])*t)*pi[j];
    //     }
    //   }
    
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
        ww += w[j]*pf[i][j];
      }
      wr2 = ww*conj(ww);
      //   printf("wr2 = %f ww = %f\n",wr2, ww);
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
    
    
    free_cvector(w, 1, max2);    
    free_cvector(wc, 1, max2);
    free_cvector(pi, 1, max2);
    free_cmatrix(pf, 1, max2, 1, max2);
    free_matrix(dens, 1, nb, 1, 2);   
    
  }
  
