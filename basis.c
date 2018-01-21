
#include "exactdiag.h"

/******************************************************************/

void construct_ham(PPP_ints *pints, PPP_par *ppp) {
    int i, j, ns;
    
    ns = pints->nsites = ppp->nb;
    
    printf("constructing Hamiltonian in site basis...\n");    
    pints->psi = matrix(1, ns, 1, ns);
    pints->hmat = matrix(1, ns, 1, ns);
    pints->con = imatrix(1, ns, 1, ns);
    pints->zmat = matrix(1, ns, 1, ns);
    pints->v_ijkl = matrix(1, ns, 1, ns);
    
    for (i = 1; i <= ns; i++)
        for (j = 1; j <= ns; j++){
            pints->psi[j][i] = 0.0;
            pints->con[i][j]  = connections(ppp, ns, i, j);
      }
    for (i = 1; i <= ns; i++)
        for (j = i; j <= ns; j++)
            pints->zmat[i][j] = pints->zmat[j][i] = zmn(ppp, i, j);
    
    for (i = 1; i <= ns; i++)
        for (j = 1; j <= ns; j++)
            pints->v_ijkl[i][j] = vijkl(ppp, i, j);
    for (i = 1; i <= ns; i++){
        for (j = i; j <= ns; j++){
            pints->hmat[i][j] = pints->hmat[j][i] = tmn(ppp, pints->v_ijkl, ns, i, j);
          
        }
    }
}

/******************************************************************/

double tmn(PPP_par *par, double **V, int ns, int i, int j) {
    double t = 0.0;
    double a_c, b_c;
    int k;
      
    a_c = par->a_c;
    b_c = par->b_c;
    
    if(par->chemical_structure == 0){
        
        if (abs(i - j) > 1) return 0;
        
        if (i == j) {
            t = a_c;
            /*     for(k=1;k<=ns;k++)       
             {
             if(i!=k) t-=V[i][k];
             }
             */
        } else {
            t = b_c;
        }
        
        return t;
        
    } 
    
    if(par->chemical_structure==1) //1=benzene
    {
        double t=0.0;
        //   printf ("U=%lf\n",par->U);
        //   printf ("Vij=%lf\n",V[i][j]);
        //printf("ns=%d\n",ns);
        
        if(i==j)
        {
            t=-0.5*par->U;
            for(k=1;k<=ns;k++)
            {
                if(i!=k) t-=V[i][k]; // *basis->zmat[i][j];      //-V[i][1]-V[i][2]-V[i][3]-V[i][4]-V[i][5]-V[i][6];//i=k
            }
        }
        
        else if(abs(i-j)==1 || abs(i-j)%5==0 && i!=j) t=-0.093; //0.5*par->U-V[i][j]; //t=-0.3539;//t=V[i][j]; //t=-2.539;
        
        else t=0.0;
        
        return t;
    }    
    
    
    
	if(par->chemical_structure==5) //5=thiophene
	{
		int a, b;		
		double tx=0.093, t=0.0;
		double t_s =0.1102; 
		double e_n=0.2866;   
		
		double pos [5][3] = {{-0.86109   ,     0.29405   ,     0.00000 },
			{-1.16796   ,    -1.03680   ,     0.00000 },
			{-2.58569   ,    -1.26697   ,     0.00000 },
			{-3.29801   ,    -0.10168   ,     0.00000 },
			{-2.27263   ,     1.28530   ,     0.00000 }};
        
        //printf("Input i=%d j=%d\n", i, j);
        
        
        
        if ((i==1 && j==2) || (i==2 && j==3) || (i==3 && j==4))
		{
            
            double d = sqrt(pow(pos[i-1][0]-pos[j-1][0],2) + 
                            pow(pos[i-1][1]-pos[j-1][1],2) + 
                            pow(pos[i-1][2]-pos[j-1][2],2));
            
            //			printf ("t -= tx + 0.1175*(d-1.397) i=%d j=%d\n",i,j);
			t = -tx + 0.1175*(d-1.397);//0.1175
			return t;
			//return t = -tx + 3.2*((d-1.397));
		}
		
		if ((i==4 && j==5) || (i==1 && j==5))
		{
            //			printf ("t -= t_s i=%d j=%d\n",i,j);
			t = -t_s;// - V[i][j];
			return t;
            
		}
        if (i==1 && j==1)
        {
            t=-1.057518;
            return t;
        }
        
        if (i==2 && j==2)
        {
            t=-1.057336;
            return t;
        }
        
        if (i==3 && j==3)
        {
            t=-1.052350;
            return t;
        }
        
        if (i==4 && j==4)
        {
            t=-1.057509;
            return t;
        }
        
        if (i==5 && j==5)
        {
            t=-1.027810;
            return t;
        }
        else
        {
			printf("Is there a last else? i=%d j=%d\n", i,j);
            return 0.0;
        }
    }
    
}

/******************************************************************/

int connections(PPP_par *par, int ns, int i, int j) {
    int t;
    
    if(par->chemical_structure==0){
        
        if (abs(i - j) > 1) return 0;;
        
        if (i == j) {
            t = 0;
        } else {
            t = 1;
        }
        return t;
    }
    
    if (par->chemical_structure == 1) {
        
        if( i == j ) t = 0;
        else if (abs(i-j) == 1) t = 1;
        
        else if(i==1 && j==6) t=1;
        else if(i==6 && j==1) t=1;
        
        else t = 0;
        
        return t;
    }
    
    if (par->chemical_structure == 5) {
        
        if( i == j ) t=0;
        else if (abs(i-j) == 1) t = 1;
        
        else if(i==1 && j==5) t=1;
        else if(i==5 && j==1) t=1;
       
        else t = 0;
        
        return t;
    }
    
}

/******************************************************************/

double zmn(PPP_par *par, int i, int j) {
    double z = 0.0;
    
    if (i - j) return 0.0;
    else if (i == j && i == 1) z = 1.0;
    else if (i == j && i == 4) z = 0.0;
    
    return z;
}

/******************************************************************/

double vijkl(PPP_par *par, int i, int j) {
    
    int k, l, nl, ns;
    double as=par->as, al=par->al, t, r, U=(par->U);///0.5; //, avl=2.6645;
    double avl, u0=0.0472; /* Bond length deviation of +/- u0 Angstrom */
    
    
    if(par->chemical_structure == 0){
        
    if (par->int_type == 0) {
        if (i == j) return U;
        else return 0.0; /* HUBBARD */
    }
    if (par->int_type == 1) {
        if (i == j) return U;
        else return U / (2.0 * fabs(i - j)); /* NORMAL COULOMB POTENTIAL */
    } else return 0.0;
        
    }
    
    if(par->chemical_structure==1) //1=benzene
    {
        double dist [] = {0.0,1.40,2.214,2.80,2.214,1.40,0};
        double alpha=0.4881;
        if (i==j) return U;
        
        else {
            double d = dist[abs((i-1)-(j-1))];
            t = U/(sqrt(1.0+alpha*(d*d)));
            return t;
        }
        
        //return t;
    } 
    
    
    if (par->chemical_structure==5) //5=thyophene
    {
        double U_s=5.0/27.2113;
        
        double pos [5][3] = {{-0.86109   ,     0.29405   ,     0.00000 },
            {-1.16796   ,    -1.03680   ,     0.00000 },
            {-2.58569   ,    -1.26697   ,     0.00000 },
            {-3.29801   ,    -0.10168   ,     0.00000 },
            {-2.27263   ,     1.28530   ,     0.00000 }};
        
        double alpha=0.48826;
        
        if (i==j)
        {
            if (i!=5)
            {
                return U;
            }
            
            else if (i==5)
            { 
                return U_s;
            }
        }
        
        else {
            
            double dd = (pow(pos[i-1][0]-pos[j-1][0],2) +
                         pow(pos[i-1][1]-pos[j-1][1],2) +
                         pow(pos[i-1][2]-pos[j-1][2],2));
            
            if (i!=5 && j!=5)
            {
                //printf ("if i=%d j=%d\n",i,j);
                return t = U/sqrt(1.0+alpha*dd);
            }
            
            else
            {
                //printf ("else i=%d j=%d\n",i,j);
                //return t = (14.397)/sqrt(pow((28.794)/(U+U_s),2)+dd);
                return t = (14.397/27.2113)/sqrt(pow((28.794/27.2113)/(U+U_s),2)+dd);
            }
            //return t;
            
        }
        
    }
    
    
    
}

/******************************************************************/

void zeroham_mat(PPP_ints *pints) {
    int i, j;
    
    printf("\nput vijkl to zero!\n");
    
    for (i = 1; i <= pints->nsites; i++) {
        for (j = 1; j <= pints->nsites; j++) {
            pints->v_ijkl[i][j] = 0.0;
            pints->hmat[i][j] = 0.0;
            pints->zmat[i][j] = 0.0;
        }
    }
    
}

/******************************************************************/
void confnumber(PPP_par *ppp, int *max1, int *max2) {
    
    int nb = ppp->nb;
    int up = ppp->up;
    int down = ppp->down;
    int i;
    
    *max1 = 1;
    for (i = nb - up + 1; i <= nb; i++) {
        *max1 = *max1 * i;
    }
    for (i = 1; i <= up; i++) {
        *max1 = *max1 / i;
    }
    if (down > 0) {
        *max2 = 1;
        for (i = nb - down + 1; i <= nb; i++) {
            *max2 = *max2 * i;
        }
        for (i = 1; i <= down; i++) {
            *max2 = *max2 / i;
        }
    } else {
        *max2 = 1;
    }
    
}

/******************************************************************/

/* make basis assuming 10>up>0, 5>down>-1 up>down or up=down */
void makebasis(PPP_ints *pints, PPP_par *ppp, int type) {
    
    int nb = ppp->nb;
    int up = ppp->up;
    int down = ppp->down;
    int size = ppp->size;
    int spin = 2;
    int i, j, k;
    int ***ps;
    
    if(type == 0) ps = pints->ps = i3tensor(1, size, 1, nb, 1, spin);
    if(type == 1) ps = pints->ps_1 = i3tensor(1, size, 1, nb, 1, spin);
    if(type == 2) ps = pints->ps_2 = i3tensor(1, size, 1, nb, 1, spin);    


    int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;
    int j1, j2, j3, j4, j5, j6, j7, j8, j9, j10;
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    int m1, m2, m3, m4, m5, m6, m7, m8, m9, m10;
    
    printf("\nMaking basis...\n");
    
    for (i = 1; i <= size; i++) {
        for (j = 1; j <= nb; j++) {
            for (k = 1; k <= spin; k++) {
                  ps[i][j][k] = 0.0;
            }
        }
    }
    
    i = 0;
   
    for (i1 = 1; i1 <= nb; i1++) {
        if (up > 1) n2 = i1 + 1;
        else n2 = nb;
        
        for (i2 = n2; i2 <= nb; i2++) {
            if (up > 2) n3 = i2 + 1;
            else n3 = nb;
            
            for (i3 = n3; i3 <= nb; i3++) {
                if (up > 3) n4 = i3 + 1;
                else n4 = nb;
                
                for (i4 = n4; i4 <= nb; i4++) {
                    if (up > 4) n5 = i4 + 1;
                    else n5 = nb;
                    
                    for (i5 = n5; i5 <= nb; i5++) {
                        if (up > 5) n6 = i5 + 1;
                        else n6 = nb;
                        
                        for (i6 = n6; i6 <= nb; i6++) {
                            if (up > 6) n7 = i6 + 1;
                            else n7 = nb;
                            
                            for (i7 = n7; i7 <= nb; i7++) {
                                if (up > 7) n8 = i7 + 1;
                                else n8 = nb;
                                
                                for (i8 = n8; i8 <= nb; i8++) {
                                    if (up > 8) n9 = i8 + 1;
                                    else n9 = nb;
                                    
                                    for (i9 = n9; i9 <= nb; i9++){
                                        if (down > 0) m1 = 1;
                                        else m1 = nb;
                                        
                                        for (j1 = m1; j1 <= nb; j1++) {
                                            if (down > 1) m2 = j1 + 1;
                                            else m2 = nb;
                                            
                                            for (j2 = m2; j2 <= nb; j2++) {
                                                if (down > 2) m3 = j2 + 1;
                                                else m3 = nb;
                                                
                                                for (j3 = m3; j3 <= nb; j3++) {
                                                    if (down > 3) m4 = j3 + 1;
                                                    else m4 = nb;
                                                    
                                                    for (j4 = m4; j4 <= nb; j4++) {
                                                        if(down > 5) m5 = j4 +1;
                                                        else m5 = nb;
                                                   
                                                    for(j5 = m5; j5 <= nb; j5++ ){
                                                        i = i + 1;
                                                        ps[i][i1][1] = 1;
                                                        
                                                        if (up >= 2) ps[i][i2][1] = 1;
                                                        if (up >= 3) ps[i][i3][1] = 1;
                                                        if (up >= 4) ps[i][i4][1] = 1;
                                                        if (up >= 5) ps[i][i5][1] = 1;
                                                        if (up >= 6) ps[i][i6][1] = 1;
                                                        if (up >= 7) ps[i][i7][1] = 1;
                                                        if (up >= 8) ps[i][i8][1] = 1;
                                                        if (up >= 9) ps[i][i9][1] = 1;
                                                        
                                                        if (down >= 1) ps[i][j1][2] = 1;
                                                        if (down >= 2) ps[i][j2][2] = 1;
                                                        if (down >= 3) ps[i][j3][2] = 1;
                                                        if (down >= 4) ps[i][j4][2] = 1;
                                                        if (down >= 5) ps[i][j5][2] = 1;
                                                      }
                                                    }
                                                    
                                                }
                                            }
                                            
                                        } // end j1
                                    } //end i9 
                                }
                                
                            }
                            
                        } // end i6
                        
                    }
                }
            }
        }
    } // end i1

    
    if(type == 0) pints->ps = ps; 
    if(type == 1) pints->ps_1 = ps; 
    if(type == 2) pints->ps_2 = ps; 
    

}
/******************************************************************/

void makehamiltonian(PPP_ints *pints, PPP_par *ppp, int type) {
    
    int nb = ppp->nb;
    int up = ppp->up;
    int down = ppp->down;
    int size = ppp->size;
    int max1 = ppp->max1;
    int max2 = ppp->max2;
    double alpha = ppp->alpha = 0.0;
    int spin = 2;
    int i, j, k, l, is, test;
    int i1, i2, i3, i4, i7, sk, nt, sj;
    double Z, expal = 1.0;
    int ***ps;
    
    if(type == 0){ 
      pints->zzz = matrix(1, size, 1, size);  
      pints->ham = cmatrix(1, size, 1, size);  
      pints->ham_fin = cmatrix(1, size, 1, size);  
    }
    
    if(type == 0) ps = pints->ps;
    if(type == 1) ps = pints->ps_1;
    if(type == 2) ps = pints->ps_2;    

    
    printf("\nConstructing hamiltonian....\n");
    
    for(i = 1; i <= size; i++){
        for(j = 1; j <= size; j++){
            pints->ham[i][j] = 0.0;
            pints->ham_fin[i][j] = 0.0;
         //if(i == 1 && j == 1 ) pints->zzz[i][j] = 1.0; 
          //else if(i == 1 && j == 2 ) pints->zzz[i][j] = 1.0; 
          //else if(i == 2 && j == 1 ) pints->zzz[i][j] = 1.0; 
          pints->zzz[i][j] = 0.0;
        }
    }  
   
    
    
    Z=(up+down)/nb;
    
    for(i1 = 1; i1 <= max2; i1++){
        for(i2 = 1; i2 <= max2; i2++){
            if(i1 == i2){
                if(ppp->int_type == 1){
                    for(k = 1; k <= nb; k++){
                        for(l = 1; l <= nb; l++){
                            if(k != l){
                                for(is = 1; is <= 2; is++){
                                    test = ps[i1][k][is];
                                    if(ps[i1][k][is] == ps[i2][l][is]){
                                      if (ppp->int_type == 0) {  // Hubbard, no background.
                                        pints->ham[i1][i2]=pints->ham[i1][i2]+0.5*pints->v_ijkl[k][l]*(test);
                                        pints->ham_fin[i1][i2]=pints->ham_fin[i1][i2]+0.5*pints->v_ijkl[k][l]*(test);
                                       
                                       }
                                       else{
                                          pints->ham[i1][i2]=pints->ham[i1][i2]+0.5*pints->v_ijkl[k][l]*(ps[i1][k][is]-Z)*(ps[i2][l][is]-Z);
                                          pints->ham_fin[i1][i2]=pints->ham_fin[i1][i2]+0.5*pints->v_ijkl[k][l]*(ps[i1][k][is]-Z)*(ps[i2][l][is]-Z);
                                    
                                        }
                                     } 
                                }
                                test=ps[i1][k][1];
                                if(ps[i1][k][1] == ps[i2][l][2]){
                                    pints->ham[i1][i2]=pints->ham[i1][i2]+pints->v_ijkl[k][l]*test;
                                    pints->ham_fin[i1][i2]=pints->ham_fin[i1][i2]+pints->v_ijkl[k][l]*test;
                                }
                            }
                            if(ps[i1][k][1]==ps[i2][l][2] && k==l){
                                pints->ham[i1][i2]=pints->ham[i1][i2]+(pints->v_ijkl[k][l])*(ps[i1][k][1]);
                                pints->ham_fin[i1][i2]=pints->ham_fin[i1][i2]+(pints->v_ijkl[k][l])*(ps[i1][k][1]);
                            }                        
                        }
                    }
                }
          
                else{
                    for(i = 1; i <= nb; i++){
                        if(ps[i1][i][1] == ps[i2][i][2]){
                            pints->ham[i1][i2]=pints->ham[i1][i2]+(ppp->U)*ps[i1][i][1];
                            pints->ham_fin[i1][i2]=pints->ham_fin[i1][i2]+(ppp->U)*ps[i1][i][1];
                        }
                    }
                }
            }
            else{
                for(is = 1; is <= 2; is++){
                    for(j = 1; j <= nb; j++){
                        for(k = 1; k <= nb; k++){
                            if(j<k) expal = cexp(I*alpha);
                            else expal = cexp(-I*alpha);
                            nt = ps[i2][k][is]-ps[i2][j][is];
                            if(pints->con[j][k] == 1 && nt == 1){
                                sk = 0;
                                for(i7=1;i7<=k-1;i7++){
                                    sk=sk+ps[i2][i7][is];
                                }
                                ps[i2][k][is] = 0;
                                ps[i2][j][is] = 1;
                                sj = 0;
                                for(i7 = 1; i7 <= j-1; i7++){
                                    sj=sj+ps[i2][i7][is];
                                }
                                test = pow(-1,sj+sk);
                                for(i3 = 1; i3 <= nb; i3++){
                                    for(i4 = 1; i4 <= 2; i4++){
                                        if(ps[i1][i3][i4] != ps[i2][i3][i4]) test = 0;
                                        
                                    } 
                                }
                               
                                pints->ham[i1][i2] +=  test*expal*(pints->hmat[j][k]);
                                //pints->zzz[i1][i2] +=  test*(pints->zmat[j][k]);
                                pints->ham_fin[i1][i2] +=  test*expal*(pints->hmat[j][k]);
                                ps[i2][k][is] = 1;
                                ps[i2][j][is] = 0;
                            }
                            
                        }
                    }
                }
            }
        }
    }
   
    for(i1 = 1; i1 <= max2; i1++){
        for(i2 = 1; i2 <= nb; i2++){
            if(ps[i1][i2][1] == 1) pints->ham[i1][i1] = pints->ham[i1][i1] + pints->vpot[i2][1];
            if(ps[i1][i2][2] == 1) pints->ham[i1][i1] = pints->ham[i1][i1] + pints->vpot[i2][2];
            
            if(ps[i1][i2][2] == 1) pints->zzz[i1][i1] = pints->zzz[i1][i1] + pints->zmat[i2][i2];
            if(ps[i1][i2][1] == 1) pints->zzz[i1][i1] = pints->zzz[i1][i1] + pints->zmat[i2][i2];
            
            if(ps[i1][i2][1] == 1) pints->ham_fin[i1][i1] = pints->ham_fin[i1][i1] + pints->vpot_fin[i2][1];
            if(ps[i1][i2][2] == 1) pints->ham_fin[i1][i1] = pints->ham_fin[i1][i1] + pints->vpot_fin[i2][2];
            
        }
    }

}


/******************************************************************/

int permutations(int x, int y) {
    
    int divider, p, permut;
    permut = 1;
    divider = 1;
    
    for (p = 1; p <= 2 * x; p++) {
        permut = permut * p;
        if (p <= y) divider = divider * p;
        if (p <= ((2 * x) - y)) divider = divider * p;
    }
    return permut = permut / divider;
    
}

/******************************************************************/


