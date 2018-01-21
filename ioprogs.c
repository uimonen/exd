
#include "exactdiag.h"
#include "num.h"


/******************************************************************/

void read_input(char *fn, PPP_ints *pints, PPP_par *ppp)
{
    int c, i, nb;
    double d;
    FILE *fname;
    
    printf("read_input...\n");  
    fname=fopen(fn, "r");
    
    fscanf(fname,"%d", &(ppp->nb));                /* Number of sites */
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    fscanf(fname,"%d", &(ppp->up));               /* Number of spin-ups */
    c=getc(fname); while(c!='\n') c=getc(fname);
    if(ppp->up > 10){
      printf("Number of spin ups is too big...\n");
       exit(1);
    }    
    
    fscanf(fname,"%d", &(ppp->down));             /* Number of spin-downs */
    c=getc(fname); while(c!='\n') c=getc(fname);
    if(ppp->up > 6){
      printf("Number of spin downs is too big...\n");
       exit(1);
    }     
 
    ppp->part = ppp->up + ppp->down;
   
    fscanf(fname,"%lf", &(ppp->U));               /* Hubbard term */
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    fscanf(fname,"%lf", &(ppp->a_c));             /* On-site parameter for the system sites */
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    fscanf(fname,"%lf", &(ppp->b_c));             /* Hopping energy for the system sites */
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    fscanf(fname,"%d",  &(ppp->int_type));        /* Type of the interaction */
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    fscanf(fname,"%d",  &(ppp->chemical_structure));        /* Type of the chemical structure */
    c=getc(fname); while(c!='\n') c=getc(fname);
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    nb = ppp->nb;
    double temp[nb+1][3];
    pints->vpot = matrix(1, nb, 1, 2);
    
    for(i = 1; i <= nb; i++){
        fscanf(fname, "%lf", &temp[i][1]);
        pints->vpot[i][1]=temp[i][1];
    }
    printf("\n");
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    for(i = 1; i <= nb; i++){
        fscanf(fname, "%lf", &temp[i][2]);
        pints->vpot[i][2]=temp[i][2];
    }
    
    c=getc(fname); while(c!='\n') c=getc(fname);
    c=getc(fname); while(c!='\n') c=getc(fname);
   
    pints->vpot_fin = matrix(1, nb, 1, 2);
    
    for(i = 1; i <= nb; i++){
        fscanf(fname, "%lf", &temp[i][1]);
        pints->vpot_fin[i][1]=temp[i][1];
    }
    printf("\n");
    c=getc(fname); while(c!='\n') c=getc(fname);
    
    for(i = 1; i <= nb; i++){
        fscanf(fname, "%lf", &temp[i][2]);
        pints->vpot_fin[i][2]=temp[i][2];
    }
        
    c=getc(fname); while(c!='\n') c=getc(fname);
    c=getc(fname); while(c!='\n') c=getc(fname);
   
    fscanf(fname,"%d", &(ppp->kick));             /* kick  */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%lf", &(ppp->E0));             /* Electric field strenght  */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%d", &(ppp->NT));             /* number of time-steps */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%lf", &(ppp->tstep));             /* time-step */
    c=getc(fname); while(c!='\n') c=getc(fname);
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%d", &(ppp->do_spectrum));             /* is the spectral function calculated */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%lf", &(ppp->lim));             /* is the spectral function calculated */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fscanf(fname,"%lf", &(ppp->acc));             /* is the spectral function calculated */
    c=getc(fname); while(c!='\n') c=getc(fname);

    fclose(fname);
    
}

/******************************************************************/
void print_matrices(PPP_ints *pints){
    
    int i,j;
    int ns = pints->nsites;
    
    printf("\n Hamiltonian in site basis:\n");
    for(i = 1; i <= ns; i++){
        for(j = 1; j <= ns; j++){
            printf("%f ",pints->hmat[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("\n Connection matrix:\n");
    for(i = 1; i <= ns; i++){
        for(j = 1; j <= ns; j++){
            printf("%d ",pints->con[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    
    printf("\n vijkl:\n");
    for(i = 1; i <= ns; i++){
        for(j = 1; j <= ns; j++){
            printf("%f ",pints->v_ijkl[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    
    printf("\n zmat (zmn):\n");
    for(i = 1; i <= ns; i++){
        for(j = 1; j <= ns; j++){
            printf("%f ",pints->zmat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    
}

/******************************************************************/

void print_parameters(char *output_file, PPP_par *ppp, PPP_ints *pints){
    FILE *fd;
    int i,j;
    
    fd = fopen(output_file,"a");
    fprintf(fd,"\n");
    fprintf(fd, "Number of sites: %d\n", ppp->nb);
    fprintf(fd, "Number of spin-ups: %d\n", ppp->up);
    fprintf(fd, "Number of spin-downs: %d\n", ppp->down);
    fprintf(fd, "Hubbard U: %f\n", ppp->U);
    fprintf(fd, "On-site parameter for the system sites: %f\n", ppp->a_c);
    fprintf(fd, "Hopping energy for the system sites: %f\n", ppp->b_c);
    fprintf(fd, "Type of interaction: %d\n", ppp->int_type);
    fprintf(fd, "Chemical structure: %d\n", ppp->chemical_structure);
    fprintf(fd,"\n");
    
    fprintf(fd, "Spin-dependent potentials for initial state: \n");
    for(i=1;i<=2; i++){
        for(j=1;j<=ppp->nb; j++){
            fprintf(fd, "%f  ", pints->vpot[j][i]);
        }
        fprintf(fd,"\n");
    }
    fprintf(fd,"\n");

    fprintf(fd, "Spin-dependent potentials for final state: \n");
    for(i=1;i<=2; i++){
        for(j=1;j<=ppp->nb; j++){
            fprintf(fd, "%f  ", pints->vpot_fin[j][i]);
        }
        fprintf(fd,"\n");
    }
    fprintf(fd,"\n");

    fprintf(fd, "Time propagation parameter: %d\n", ppp->kick);
    fprintf(fd, "Electric field strength: %f\n", ppp->E0);   
    fprintf(fd, "Number of time-steps: %d\n", ppp->NT);
    fprintf(fd, "Length of time-step: %f\n", ppp->tstep);
    fprintf(fd,"\n");

    fprintf(fd, "Is the specral function calculated: %d ", ppp->do_spectrum);
    if(ppp->do_spectrum == 1) fprintf(fd, "Yes.\n");
    else fprintf(fd, "No.\n");
    fprintf(fd, "Limits of the plot for the spectral function: %f \n", ppp->lim);
    fprintf(fd, "Accuracy of the specral function calculation: %f, 1/ acc = %f\ ", ppp->acc, 1/(ppp->acc));
    fprintf(fd,"\n");
    fprintf(fd,"\n");
    fprintf(fd,"\n");


    fclose(fd);
    
}

/******************************************************************/

void printbasis(char *filename, PPP_ints *pints, PPP_par *ppp, int type){
    
    int nb = ppp->nb;
    int max2 = ppp->max2;
    int hub = ppp->int_type;
    int size = ppp->size;
    int wp[size];
    double sum;
    int ***ps;    

    int i,j, i1, i2, i3, k, l;
    FILE *fd;

    if(type == 0) ps = pints->ps;
    if(type == 1) ps = pints->ps_1;
    if(type == 2) ps = pints->ps_2;

    
    for(i=1;i<=size;i++){
        wp[i] = 0.0;
    }
    
    fd = fopen(filename,"w");
    for(i=1;i<=max2;i++){
        for(i2=1; i2<=nb; i2++){
            for(i3=1; i3<=nb; i3++){
                if(pints->hmat[i2][i3]==-1){
                    if(ps[i][i2][1]==1 && ps[i][i3][2]==1 && ps[i][i2][2] == 0 && ps[i][i3][1] == 0){
                        wp[i]=wp[i]+1;
                    }
                }
                
            }
        }
        fprintf(fd,"%d  %d  ",i, wp[i]);
        for(i2=1;i2<=nb;i2++){
            fprintf(fd,"%d",ps[i][i2][1]);
        }
        fprintf(fd,"\n");
        fprintf(fd,"      ");
        for(i3=1;i3<=nb;i3++){
            fprintf(fd,"%d", ps[i][i3][2]);
        }
        fprintf(fd,"\n");
        
    }
    
    fclose(fd);
    
}

/******************************************************************/

void printhamiltonian(char *filename, int max2, double _Complex **ham){
    
    int i, j;
    FILE *fd;
    fd = fopen(filename,"w");
    
    for(i=1;i<=max2;i++){
        for(j=1;j<=max2;j++){
            fprintf(fd,"%5.5f  %5.5f ", creal(ham[i][j]), cimag(ham[i][j]));
        }
        fprintf(fd,"\n");
    }
    
    fclose(fd);
    
}

/******************************************************************/

void printdensities(char *filename, double _Complex **UH, PPP_ints *pints, PPP_par *ppp){
    
    int i,j,k;
    double sup, sdown;
    
    int max2 = ppp->max2;
    int nb = ppp->nb;
    FILE *fd;
    fd = fopen(filename,"w");
    
    
    for(k = 1; k <= max2; k++){
        fprintf(fd,"%d\n",k);
        for(j = 1; j <= nb; j++){
            sup = 0.0;
            sdown = 0.0;
            for(i = 1; i <= max2; i++){
                if(pints->ps[i][j][1] == 1) sup = sup + UH[k][i]*conj(UH[k][i]);
                if(pints->ps[i][j][2] == 1) sdown = sdown + UH[k][i]*conj(UH[k][i]);
            }
            
            fprintf(fd,"%d  %10.10f  %10.10f  %10.10f\n", j, sup+sdown, sup, sdown);
        }
    }
    
    
    
}

/******************************************************************/

void printwf(char *filename, double _Complex **UH, PPP_ints *pints, PPP_par *ppp){

    int i1, i2;
    FILE *fd;
    int max2 = ppp->max2;   
    
    fd = fopen(filename,"w");
    
    for(i1 = 1; i1 <= max2; i1++){
        fprintf(fd,"%d  state \n", i1);
        for(i2 = 1; i2 <= max2; i2++){
            if(i1%2 != 0)fprintf(fd,"%d  %20.20e \n",i2, -UH[i1][i2]);  
            else fprintf(fd,"%d  %20.20e \n",i2, UH[i1][i2]); 
        }
        
    }
    fclose(fd);


}




/******************************************************************/

void printzzz(char *filename, int max2, double  **ham){
    
    int i, j;
    FILE *fd;
    fd = fopen(filename,"w");
    
    for(i=1;i<=max2;i++){
        for(j=1;j<=max2;j++){
            fprintf(fd,"%5.5f ", ham[i][j]);
        }
        fprintf(fd,"\n");
    }
    
    fclose(fd);
    
}
