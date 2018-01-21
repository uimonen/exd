/*
==============================================================================
 * 	Exactdiagonalization of the Hilber space for Hubbard and 
 *  Coulombic systems.w
 *  FILES:  exactdiag.c basis.c diag.c ioprogs.c num.c outputs.c propagate.c      
==============================================================================
 */

#include <stdio.h>
#include "exactdiag.h"


int main()
{
    FILE *fd;
    printf("\n\n\nexact digonalization. \n\n\n");
    
    fd = fopen(OUTPUT, "w");
    fprintf(fd, "exact digonalization. \n");
    fclose(fd);

    int nb,size, i, j;
    double **U;
    double *e;
    PPP_ints *pints;
    PPP_par *ppp;

    pints = (PPP_ints *) malloc(sizeof(PPP_ints));
    ppp = (PPP_par *) malloc(sizeof(PPP_par));
    
    read_input(INPUT, pints, ppp);
    print_parameters(OUTPUT, ppp, pints);
    zeroham_mat(pints);
    construct_ham(pints, ppp);
    print_matrices(pints);
    
    nb = ppp->nb;
    U = matrix(1,nb,1,nb);
    e = vector(1,nb);    
 
    diag(nb,pints->hmat,U,e);
    confnumber(ppp, &(ppp->max1), &(ppp->max2));
    ppp->max2 =(ppp->max2) * (ppp->max1);

    fd = fopen(OUTPUT,"a");
    fprintf(fd,"Huckel eigenvalues:\n");
    for (i=1; i<=nb; i++) {
      fprintf(fd,"%d %f\n", i, e[i]);
    }
    fclose(fd);
  
    size = ppp->max2;
    ppp->size = size;

    
    pints->UH = cmatrix(1,size,1,size);
    pints->UH_t = cmatrix(1,size,1,size);
    pints->eh = vector(1,size);
    pints->UH_fin = cmatrix(1,size,1,size);
    pints->eh_fin = vector(1,size);

   
    makebasis(pints, ppp, 0);
    printbasis("basis.dat", pints, ppp, 0);
    makehamiltonian(pints, ppp, 0);
    printhamiltonian("hamiltonian.dat",size, pints->ham);
    printhamiltonian("hamiltonian_fin.dat",size, pints->ham_fin);
    printzzz("zzz.dat",size, pints->zzz);

    if(ppp->kick){          // Make sure no external potential is applied.
      for(i = 1; i<= nb; i++){
        for(j = 1; j<=2; j++){
          pints->vpot[i][j] = 0.0;
          pints->vpot_fin[i][j] = 0.0;
        }
      }
    }    

    // Determine the initial and final state eigenvectors and eigenvalues.
    cdiag(size,pints->ham,pints->UH,pints->eh);
    cdiag(size,pints->ham_fin,pints->UH_fin,pints->eh_fin);

    fd = fopen(OUTPUT,"a");
    fprintf(fd,"\nEigenvalues initial state:\n");
    for (i = 1; i <= size; i++) {
      fprintf(fd,"%d %f\n", i, pints->eh[i]);
    }
    fprintf(fd,"\nEigenvalues final state:\n");
    for (i = 1; i <= size; i++) {
      fprintf(fd,"%d %f\n", i, pints->eh_fin[i]);
    }
    fclose(fd);
 
    printdensities("densities.dat",pints->UH,pints, ppp);
    printdensities("densities_fin.dat",pints->UH_fin, pints,ppp);
    printwf("wavefunction.dat",pints->UH,pints, ppp);
    printwf("wavefunction_fin.dat",pints->UH_fin, pints,ppp);

    propagate(pints, ppp);
    if(ppp->do_spectrum == 1 ) dospctrum(pints, ppp);
 
    free_vector(e, 1, nb);
    free_matrix(U, 1, nb, 1, nb);
    free_matrix(pints->psi, 1, nb, 1, nb);
    free_matrix(pints->zmat, 1, nb, 1, nb);
    free_matrix(pints->v_ijkl, 1, nb, 1, nb);
    free_matrix(pints->hmat, 1, nb, 1, nb);
    free_matrix(pints->vpot, 1, nb, 1, 2);
    free_cmatrix(pints->ham, 1, size, 1, size);
    free_i3tensor(pints->ps, 1, size, 1, nb, 1, 2);
    free_cmatrix(pints->UH, 1, size, 1, size);
    free_vector(pints->eh, 1, size);

    printf("\nEND!\n\n");

 return 0;
}


