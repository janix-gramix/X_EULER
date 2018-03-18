//
//  PB051_100.c
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "PB051_100.h"

#define PB051_MAXP 1000000
#define PB051_SQMAXP 3200
#define PB051_MAXGDIG   7
#define PB051_MINLGSUI  2
#define PB051_MAXLGSUI  8

typedef struct ListDeltas {
    int nb;         // nombre de delta possibles
    int dt[1<<(PB051_MAXGDIG-1)] ;  // valeurs des deltas
} ListDeltas ;

int PB051(PB_RESULT *pbR) {
    ListDeltas  L_delta[1<<(PB051_MAXGDIG-1)] ;
    CTX_PRIMETABLE * ctxP  ;
    int minSuit = PB051_MINLGSUI ;
    T_prime p ;
    pbR->nbClock = clock()  ;
    if( PB051_SQMAXP*PB051_SQMAXP< 10*PB051_MAXP) {
        fprintf(stdout,"\t PB%d Need more Prime %d < %lld \n",pbR->pbNum,PB051_SQMAXP,Sqrt64(10*PB051_MAXP));
        return 0 ;
        
    }
    if((ctxP = Gen_tablePrime(PB051_SQMAXP)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    {
        int i ;
        int deltas[1<<(PB051_MAXGDIG-1)] ;
        int  lg = 0 , pow10 = 1 ;
        L_delta[0].nb = 0 ;
        deltas[lg++] = 0 ;
        for(i=1; i<PB051_MAXGDIG;i++ , lg *= 2) {
            int j;
            pow10 *= 10 ;
            for(j=0;j<lg;j++) {
                int nPat = lg+j ;
                deltas[nPat] = deltas[j] +pow10 ;
                int k ;
                L_delta[nPat].nb = 0 ;
                for(k=1;k<=nPat;k++) {
                    if((k & nPat) == k) { // k est-il contenu dans le pattern nPat
                        L_delta[nPat].dt[L_delta[nPat].nb++] = deltas[k] ;
                    }
                }
                
            }
        }
    }
    
    int incP ;
    for(p=5,incP=2;p<PB051_MAXP;p+=incP , incP = 6 - incP) {
        u_int8_t occurDig[10] ;
        int k ;
        if(!Is_Prime(p,tbPrime)) continue ;
        memset(occurDig,0,10 - minSuit) ;
        {
            int p1 ;
            u_int8_t bit = 1; // on extrait les digits en sautant le premier (poids faible)
            for(p1=p/10;p1 != 0; p1 /= 10 , bit <<= 1 ) {
                u_int8_t dg = p1 % 10 ;
                if(dg < 10 - minSuit) { // on construi le pattern pour de digit
                    occurDig[dg] |= bit ;
                }
            }
        }
        for(k=0;k<10 - minSuit;k++){ // on traite les patterns des digits 0, 1, ... pour lesquels le minimum est valide
            if(occurDig[k]) {
                ListDeltas * L_d = L_delta +  occurDig[k] ; // liste des deltas pour le pattern
                int tbP[10] ;
                int id ;
                for(id=0;id<L_d->nb;id++) { // listes des deltas .
                    int j,nbP ;
                    tbP[0]= p ;
                    nbP = 1;
                    int delta = L_d->dt[id]  ;
                    for(j=1;j<10-k;j++) {
                        if(Is_Prime(p+j*delta,tbPrime)) tbP[nbP++] = p+j*delta ;
                    }
                    if(nbP >= minSuit) {
                        if(pbR->isVerbose) {
                            int i; T_prime p1 ;
                            char pattern[PB051_MAXGDIG];
                            for(i=0, p1=p;p1!=0;i++ , p1 /= 10, delta /= 10) {
                                if((delta % 10) == 1) {
                                    pattern[i] = '*' ;
                                } else {
                                    pattern[i] = (p1 % 10) + '0' ;
                                }
                            }
                            fprintf(stdout,"\t PB%d ",pbR->pbNum);
                            while(i-- > 0) fprintf(stdout,"%c",pattern[i]);
                            fprintf(stdout,"->[%d] ",nbP);
                            for(i=0;i<nbP;i++) fprintf(stdout,"%d%c",tbP[i],(i==nbP-1) ? '\n' : ',');
                        }
                        minSuit = nbP + 1 ;
                        if(minSuit > PB051_MAXLGSUI ) {
                            Free_tablePrime(ctxP);
                            sprintf(pbR->strRes,"%d",p) ;
                            pbR->nbClock = clock() - pbR->nbClock ;
                            return 1 ;
                        }
                    }
                }
            }
        }
    }
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}

#define PB052_MAXD  10
// pour que x et 6x ait le m nb de digit x commence par 1
// de ce fait il y a au moins 6 digits differents
// donc x est >= 123456
int PB052(PB_RESULT *pbR) {
    int i;
    pbR->nbClock = clock()  ;
    int nbDig = 6 ;
    u_int8_t digX[PB052_MAXD+1] = { 6,5,4,3,2,1,0,0,0,0,0} ;
    while(nbDig < PB052_MAXD) {
        int k ;
        u_int8_t digXsort[PB052_MAXD] ;
        u_int8_t digXxk[PB052_MAXD+1] ;
        memcpy(digXsort,digX,nbDig);
        HeapSortUint8(digXsort,nbDig) ;
        for(k=6;k>1;k--) {
            int i ;
            for(i=0;i<nbDig;i++) {
            }
            // on ajoute digXxk + digX
            digXxk[0] = 0 ;
            for(i=0;i < nbDig ;i++) {
                digXxk[i] += k * digX[i];
                digXxk[i+1] = 0 ;
                if(digXxk[i] >= 10) {
                    digXxk[i+1] =  digXxk[i] / 10;
                    digXxk[i] = digXxk[i] % 10 ;
                }
            }
            if(digXxk[nbDig])  break ;
            HeapSortUint8(digXxk,nbDig) ;
            if(memcmp(digXsort,digXxk,nbDig) != 0)  break ;
        }
        if(k==1) {
            int i ;
            u_int32_t n =0 ;
            for(i=1;i<=nbDig;i++) { n = n*10 + digX[nbDig-i] ; }
            if(pbR->isVerbose)  fprintf(stdout,"\t PB%d %d,%d,%d,%d,%d,%d \n",pbR->pbNum,n,2*n,3*n,4*n,5*n,6*n) ;
            sprintf(pbR->strRes,"%d",n);
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
        if(digXxk[nbDig]) {
            nbDig++ ;
            for(i=0;i<=6;i++) {digX[i] = 6 -i ;}
            for(;i<nbDig;i++) digX[i] = 1 ;
            printf(" nbDig=%d",nbDig) ;
            
        } else {
            digX[0] ++ ;
            for(i=0;digX[i]>=10;i++) {
                digX[i] -= 10 ;
                digX[i+1]++ ;
            }
        }
        
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}

int PB052a(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nbDig = 6 ;
    u_int32_t n = 123456 ;
    while(nbDig < PB052_MAXD) {
        int i,k ;
        u_int32_t N1 ;
        u_int8_t digXsort[PB052_MAXD] ;
        u_int8_t digXxk[PB052_MAXD+1] ;
        for(i=0,N1=n;N1 != 0; i++) {
            digXsort[i] = N1 % 10 ; N1 /= 10 ;
        }
        HeapSortUint8(digXsort,nbDig) ;
        for(k=6;k>1;k--) {
            for(i=0,N1=n*k;N1 != 0; i++) {
                digXxk[i] = N1 % 10 ; N1 /= 10 ;
            }
            if(digXxk[nbDig])  break ;
            HeapSortUint8(digXxk,nbDig) ;
            if(memcmp(digXsort,digXxk,nbDig) != 0)  break ;
        }
        if(k==1) {
            if(pbR->isVerbose)  fprintf(stdout,"\t PB%da %d,%d,%d,%d,%d,%d \n",pbR->pbNum,n,2*n,3*n,4*n,5*n,6*n) ;
            sprintf(pbR->strRes,"%d",n);
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
        if(digXxk[nbDig]) {
            n = 0 ;
            nbDig++ ;
            for(i=nbDig;i>=6;i--) { n = 10*n + 1 ; }
            for(;i>=0;i--) {n = 10*n + 6 - i ;}
            printf(" nbDig=%d",nbDig) ;
        } else {
            n++ ;
        }
        
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB053_MAXN  100
#define PB053_MINV  1000000
int PB053(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int32_t n, n2,r,Cr  ;
    int32_t nbSup = 0;
    for(n=2,n2=0,r=0,Cr=1;n<=PB053_MAXN;n++) {
        if((n & 1) == 0 ){
            n2++ ;
        }
        // calculate C(r,n)= (C(r,n-1) * n) / (n-r)
        Cr = (Cr * n) / (n-r) ;
        if( Cr > PB053_MINV ) {
            do {
                // calculate C(r-1,n) = (C(r,n) * r)/(n-r+1)
                Cr = (Cr * r)/(n-r+1) ;
                r-- ;
            } while(Cr > PB053_MINV);
        } else {
            while(r < n2) {
                // calculate C(r+1,n) = (C(r,n) * (n-r))/(r+1)
                int32_t nCr = (Cr * (n-r)) /(r+1) ;
                if(nCr > PB053_MINV) break ;
                Cr = nCr ;
                r++ ;
            }
        }
        if(r < n2) {
            nbSup += n-1 - 2*r ;
        }
    }
    sprintf(pbR->strRes,"%d",nbSup);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB055_MAXN  10000
#define PB055_MAXITER   50
#define PB055_MAXDIGIT  PB055_MAXITER+5

int PB055(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int32_t n  ;
    int32_t nbLychrel = 0;
    for(n=1;n<PB055_MAXN;n++) {
        u_int8_t Dig0[PB055_MAXDIGIT] ;
        u_int8_t Dig1[PB055_MAXDIGIT] ;
        u_int8_t *pDigCur ;
        u_int8_t *pDigAnt ;
        
        int nbDig,k,n1  ;
        for(nbDig=0,n1=n;n1 != 0;) {
            Dig0[nbDig++] = n1 % 10 ;
            n1 /= 10 ;
        }
        pDigCur = Dig0 ;
        pDigAnt = Dig1 ;
        for(k=0;k<PB055_MAXITER;k++) {
            int i ;
            u_int8_t carry = 0 ;
            {
                u_int8_t *tmp = pDigCur;
                pDigCur = pDigAnt ;
                pDigAnt = tmp ;
            }
            for(i=0;i<nbDig;i++) {
                pDigCur[i] = pDigAnt[i]+pDigAnt[nbDig-i-1] + carry ;
                if(pDigCur[i] >= 10) {
                    pDigCur[i] -= 10 ; carry = 1 ;
                } else {
                    carry = 0 ;
                }
            }
            if( carry) pDigCur[nbDig++] = 1 ;
            for(i=0;2*i<nbDig;i++) {
                if(pDigCur[i] != pDigCur[nbDig-i-1]) break ;
            }
            if(2*i>=nbDig) break ;
        }
        if(k==PB055_MAXITER) {
            nbLychrel++ ;
        }
    }
    sprintf(pbR->strRes,"%d",nbLychrel);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



#define PB057_N  1000
#define PB057_MAXDIGIT  PB057_N+5

typedef struct FractCont {
    u_int8_t *pt1 ;
    u_int8_t *pt ;
    int nbDig ;
}  FractCont ;
int PB057(PB_RESULT *pbR) {
    int n ;
    int nbNsupD = 0 ;
    pbR->nbClock = clock()  ;
    u_int8_t DigN1[PB057_MAXDIGIT] ;
    u_int8_t DigN0[PB057_MAXDIGIT] ;
    u_int8_t DigD1[PB057_MAXDIGIT] ;
    u_int8_t DigD0[PB057_MAXDIGIT] ;
    FractCont FC[2] ; // numerateur indice 0, denominateur indice 1
    FC[0].nbDig = 1 ;
    FC[0].pt = DigN0 ;
    FC[0].pt1 = DigN1 ;
    FC[0].pt[0] = 1 ;
    FC[0].pt1[0] = 1 ;
    
    FC[1].nbDig = 1 ;
    FC[1].pt = DigD0 ;
    FC[1].pt1 = DigD1 ;
    FC[1].pt[0] = 1 ;
    FC[1].pt1[0] = 0 ;
    for(n=1;n<=PB057_N;n++) {
        int i ;
        int k ;
        for(k=0;k<2;k++) {
            u_int8_t *tmp = FC[k].pt1 ;
            FC[k].pt1  = FC[k].pt ; FC[k].pt = tmp ;
            // N(n+1) = 2*N(n) + N(n-1)
            int carry = 0 ;
            for(i=0;i<FC[k].nbDig;i++) {
                FC[k].pt[i] += 2 * FC[k].pt1[i] + carry ;
                if(FC[k].pt[i] >= 10) {
                    carry = FC[k].pt[i] / 10 ;
                    FC[k].pt[i] = FC[k].pt[i] % 10 ;
                } else {
                    carry = 0 ;
                }
            }
            if(carry) {
                FC[k].pt1[FC[k].nbDig] = 0 ;
                FC[k].pt[FC[k].nbDig++] = carry ;
            }
        }
        if(FC[0].nbDig > FC[1].nbDig) {
            nbNsupD++ ;
        }
        
    }
    sprintf(pbR->strRes,"%d",nbNsupD);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

// N = 2n+1
// DownRight N*N
// DownLeft  N*(N-1)+1
// UpLeft    N*(N-2)+2
// UpRight   N*(N-3)+3
// NbPoint = 4n+1
#define PB058_MAXN      100000
int PB058(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    int N ;
    int nbPoint, nbPrime , DL  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB058_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    nbPoint = 1 ;
    nbPrime = 0 ;
    DL = 1 ;
    const T_prime * tbPrime =GetTbPrime(ctxP);
    for(N=3;N< PB058_MAXN;N += 2 ) {
        nbPoint += 4 ;
        DL += 4*N - 6 ;
        if(Is_Prime(DL,tbPrime)) {
            nbPrime++ ;
        }
        if(Is_Prime(DL-N+1,tbPrime)) {
            nbPrime++ ;
        }
        if(Is_Prime(DL-2*N+2,tbPrime)) {
            nbPrime++ ;
        }
        if(10*nbPrime < nbPoint) break ;
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(N<PB058_MAXN) {
        sprintf(pbR->strRes,"%d",N);
        return 1 ;
    } else {
        return 0 ;
    }
}

static u_int8_t PB59_encrypted[] = {
    79,59,12,2,79,35,8,28,20,2,3,68,8,9,68,45,0,12,9,67,68,4,7,5,23,27,1,21,79,85,78,79,85,71,38,10,71,27,12,2,79,6,2,8,13,9,1,13,9,8,68,
    19,7,1,71,56,11,21,11,68,6,3,22,2,14,0,30,79,1,31,6,23,19,10,0,73,79,44,2,79,19,6,28,68,16,6,16,15,79,35,8,11,72,71,14,10,3,79,12,2,79,
    19,6,28,68,32,0,0,73,79,86,71,39,1,71,24,5,20,79,13,9,79,16,15,10,68,5,10,3,14,1,10,14,1,3,71,24,13,19,7,68,32,0,0,73,79,87,71,39,1,71,
    12,22,2,14,16,2,11,68,2,25,1,21,22,16,15,6,10,0,79,16,15,10,22,2,79,13,20,65,68,41,0,16,15,6,10,0,79,1,31,6,23,19,28,68,19,7,5,19,79,12,
    2,79,0,14,11,10,64,27,68,10,14,15,2,65,68,83,79,40,14,9,1,71,6,16,20,10,8,1,79,19,6,28,68,14,1,68,15,6,9,75,79,5,9,11,68,19,7,13,20,79,8,
    14,9,1,71,8,13,17,10,23,71,3,13,0,7,16,71,27,11,71,10,18,2,29,29,8,1,1,73,79,81,71,59,12,2,79,8,14,8,12,19,79,23,15,6,10,2,28,68,19,7,22,
    8,26,3,15,79,16,15,10,68,3,14,22,12,1,1,20,28,72,71,14,10,3,79,16,15,10,68,3,14,22,12,1,1,20,28,68,4,14,10,71,1,1,17,10,22,71,10,28,19,6,
    10,0,26,13,20,7,68,14,27,74,71,89,68,32,0,0,71,28,1,9,27,68,45,0,12,9,79,16,15,10,68,37,14,20,19,6,23,19,79,83,71,27,11,71,27,1,11,3,68,2,
    25,1,21,22,11,9,10,68,6,13,11,18,27,68,19,7,1,71,3,13,0,7,16,71,28,11,71,27,12,6,27,68,2,25,1,21,22,11,9,10,68,10,6,3,15,27,68,5,10,8,14,
    10,18,2,79,6,2,12,5,18,28,1,71,0,2,71,7,13,20,79,16,2,28,16,14,2,11,9,22,74,71,87,68,45,0,12,9,79,12,14,2,23,2,3,2,71,24,5,20,79,10,8,27,
    68,19,7,1,71,3,13,0,7,16,92,79,12,2,79,19,6,28,68,8,1,8,30,79,5,71,24,13,19,1,1,20,28,68,19,0,68,19,7,1,71,3,13,0,7,16,73,79,93,71,59,12,
    2,79,11,9,10,68,16,7,11,71,6,23,71,27,12,2,79,16,21,26,1,71,3,13,0,7,16,75,79,19,15,0,68,0,6,18,2,28,68,11,6,3,15,27,68,19,0,68,2,25,1,21,
    22,11,9,10,72,71,24,5,20,79,3,8,6,10,0,79,16,8,79,7,8,2,1,71,6,10,19,0,68,19,7,1,71,24,11,21,3,0,73,79,85,87,79,38,18,27,68,6,3,16,15,0,17,
    0,7,68,19,7,1,71,24,11,21,3,0,71,24,5,20,79,9,6,11,1,71,27,12,21,0,17,0,7,68,15,6,9,75,79,16,15,10,68,16,0,22,11,11,68,3,6,0,9,72,16,71,29,
    1,4,0,3,9,6,30,2,79,12,14,2,68,16,7,1,9,79,12,2,79,7,6,2,1,73,79,85,86,79,33,17,10,10,71,6,10,71,7,13,20,79,11,16,1,68,11,14,10,3,79,5,9,11,
    68,6,2,11,9,8,68,15,6,23,71,0,19,9,79,20,2,0,20,11,10,72,71,7,1,71,24,5,20,79,10,8,27,68,6,12,7,2,31,16,2,11,74,71,94,86,71,45,17,19,79,16,
    8,79,5,11,3,68,16,7,11,71,13,1,11,6,1,17,10,0,71,7,13,10,79,5,9,11,68,6,12,7,2,31,16,2,11,68,15,6,9,75,79,12,2,79,3,6,25,1,71,27,12,2,79,
    22,14,8,12,19,79,16,8,79,6,2,12,11,10,10,68,4,7,13,11,11,22,2,1,68,8,9,68,32,0,0,73,79,85,84,79,48,15,10,29,71,14,22,2,79,22,2,13,11,21,
    1,69,71,59,12,14,28,68,14,28,68,9,0,16,71,14,68,23,7,29,20,6,7,6,3,68,5,6,22,19,7,68,21,10,23,18,3,16,14,1,3,71,9,22,8,2,68,15,26,9,6,1,
    68,23,14,23,20,6,11,9,79,11,21,79,20,11,14,10,75,79,16,15,6,23,71,29,1,5,6,22,19,7,68,4,0,9,2,28,68,1,29,11,10,79,35,8,11,74,86,91,68,52,
    0,68,19,7,1,71,56,11,21,11,68,5,10,7,6,2,1,71,7,17,10,14,10,71,14,10,3,79,8,14,25,1,3,79,12,2,29,1,71,0,10,71,10,5,21,27,12,71,14,9,8,1,3,
    71,26,23,73,79,44,2,79,19,6,28,68,1,26,8,11,79,11,1,79,17,9,9,5,14,3,13,9,8,68,11,0,18,2,79,5,9,11,68,1,14,13,19,7,2,18,3,10,2,28,23,73,79,
    37,9,11,68,16,10,68,15,14,18,2,79,23,2,10,10,71,7,13,20,79,3,11,0,22,30,67,68,19,7,1,71,8,8,8,29,29,71,0,2,71,27,12,2,79,11,9,3,29,71,60,11,
    9,79,11,1,79,16,15,10,68,33,14,16,15,10,22,73
    
} ;

#define PB059_MAXASCII  128
int PB059(PB_RESULT *pbR) {
    u_int16_t HIST[PB059_MAXASCII*3] ;
    int i;
    int Sum = 0 ;
    pbR->nbClock = clock()  ;
    memset(HIST,0,sizeof(HIST));
    for(i=0;i<sizeof(PB59_encrypted);i++) {
        HIST[ (i % 3) * PB059_MAXASCII + PB59_encrypted[i]]++ ;
    }
    for(i=0;i<3;i++) {
        int j, jmax = 0 ;
        int max = 0 ;
        for(j=0;j<PB059_MAXASCII;j++) {
            if(HIST[j+i*PB059_MAXASCII]>max) {
                max = HIST[j+i*PB059_MAXASCII] ;
                jmax = j ;
            }
        }
        jmax ^= ' ' ;
        for(j=i;j<sizeof(PB59_encrypted);j += 3) {
            PB59_encrypted[j] ^= jmax ;
            Sum += PB59_encrypted[j] ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d \"%30.30s ...\"\n",pbR->pbNum,PB59_encrypted) ;
    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum);
    return 1 ;
}


#define PB060_MAXP      100000



int PB060(PB_RESULT *pbR) {
    static u_int8_t *isP1P2 ;
    
    static int NPMAX,NP ;
    CTX_PRIMETABLE * ctxP  ;
    int32_t *pow10 ;
    int32_t maxS, minS = 0 ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB060_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    u_int32_t nbPrime = GetNbPrime(ctxP) ;
    NPMAX = GetNbPrime(ctxP)-1 ;
    int32_t maxP = tbPrime[nbPrime-1] ;
    {   // on limite la valeur max des nb premiers a tester pour que la concat de 2 d'entre eux
        // soit testable avec la table calculee
        int32_t pow10 = 10 ;
        while(pow10 < PB060_MAXP) pow10 *= 10 ;
        u_int64_t maxP2 = maxP * (u_int64_t) maxP ;
        while(tbPrime[NP] * (u_int64_t) (pow10 + 1) >= maxP2 ) NPMAX-- ;
    }
    
    NP = NPMAX ;
    maxS = tbPrime[NP-1] ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d NPMAX=%d,maxS=%d,maxP=%d\n",pbR->pbNum,NP,maxS,maxP);
    pow10 = calloc(NP, sizeof(pow10[0])) ;
    {
        int i ;
        for(i=0;i<NP;i++) {
            int32_t Pi = tbPrime[i] ;
            int32_t ip10 = 10 ;
            while(ip10< Pi) ip10 *= 10 ;
            pow10[i] = ip10 ;
        }
    }
    {
        isP1P2 = calloc(NP*NP,sizeof(isP1P2[0])) ;
        // parcours arborescent jusqu'a atteindre la profondeur 5
        int index[5],i ;
        u_int64_t P[5],SP ;
        if(minS==0) maxS = tbPrime[NP-1] ;
        
        for(i=0,index[i]=0,SP=0;i>=0;) {
            int j , isOK ;
            int ip = index[i] ;
            isOK = 1 ;
            // borne sur l'index et test que la somme P[k] < maxS avec P[k] croissant
            if(ip >=NP || SP+(5-i)*(P[i]=tbPrime[ip]) > maxS ) {
                if(--i >= 0) { SP -= P[i] ; index[i]++ ; }
                continue ;
            } else {
                for(j=0;j<i;j++) {
                    int jp = index[j] ; // test compatibilite P[j] et le nouveau P[i]
                    // test non calcule (zero , sinon 1 pour non compat et 2 pour compat
                    if(isP1P2[NP*jp+ip]==0)  { isP1P2[NP*jp+ip] = Is_Prime2(P[i]+P[j]*pow10[ip],P[j]+P[i]*pow10[jp],tbPrime) ? 2 : 1 ; }
                    if(isP1P2[NP*jp+ip]==1)  {
                        isOK = 0; break ;
                    }
                }
            }
            if(isOK) {
                if(i == 4) {
                    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d  %lld = %lld + %lld + %lld + %lld + %lld\n",pbR->pbNum,SP+P[4],P[0],P[1],P[2],P[3],P[4]);
                    if(SP+P[4] < maxS) {
                        minS = maxS =(int32_t) (SP+P[4]) ;
                    }
                    index[i]++ ;
                } else {
                    SP += P[i] ;
                    index[i+1] = index[i]+1 ;
                    i++ ;
                }
            } else {
                index[i]++ ;
            }
        }
        free(isP1P2) ;
    }
    free(pow10);
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(minS) {
        sprintf(pbR->strRes,"%d",minS) ;
        return 1 ;
    } else {
        return 0 ;
    }
}



// maximum compris entre 1010 et 9999
// P(s,n) = (s-2)n(n-1)/2 + n
// P(s+1,n) = P(s,n) + n(n-1)/ 2 ;

//  nombre de polygones differents
#define PB061_NS       6

// max de valeurs pour le gros polygone (choisi comme point de depart de la boucle)
#define PB061_MAXPX0     100
// nombre de prefixe differents (en fait entre 10 et 99) Majorant
#define PB061_NBPREF    100
// majorant stric du nombre de polygone d'un meme type pouvent tomber sur un prefixe
#define PB061_MAXBYPREF 4

// hors le premier type de polygone (permTyp=PB061_NS-1 choisi,s=0 où toutes les valeurs sont stockees a la suite)
// pour les autres on indexe les polygones par prefixe,typePoly et si polygones plusieurs ont
// le meme index on les mets a la suite terminee par un zero (Hyp : max conflits < PB061_MAXBYPREF)
// de ce fait on adresse directement la tete de liste par la macro IND(pref,s)
#define IND(pref,s) (PB061_MAXPX0+ ((pref)*PB061_NS+(s)) * PB061_MAXBYPREF)

typedef u_int16_t  T_Polygonal;
int PB061(PB_RESULT *pbR) {
    int n,s ;
    int nbSol = 0 ;
    int isOnlyFirst = 0 ;
    int Smin ;
    T_Polygonal  *T_values;
    ;
    char *     nameT[6] = { "Tria","Squa","Pent","Hexa","Hept","Octo" } ;
    
    pbR->nbClock = clock()  ;
    
    {
        u_int32_t T0 = 0, Ts ;
        T_values = calloc(PB061_MAXPX0+(PB061_NS*PB061_MAXBYPREF*PB061_NBPREF),sizeof(T_values[0])) ;
        for(n=0;T0 < 10000 ;n++ , T0 += n) {
            for(s=0,Ts=T0;s<PB061_NS;s++ ,Ts += (n*(n-1))/2) {
                if(Ts>=1010 && Ts<10000 && ((Ts % 100) >= 10) ) {
                    // index par, prefixe et type=s
                    int ind = (s==PB061_NS-1) ? 0 : IND(Ts/100,s);
                    while(T_values[ind]) ind++ ; // on cherche la premiere valeur non nulle
                    if(s ==PB061_NS-1) {
                        T_values[ind] = Ts ;
                    }else {
                        if(ind >= IND(Ts/100,s)+ PB061_MAXBYPREF-1) {
                            fprintf(stdout,"\t PB%0.3d Too many identical prefix(%d) pour Type %d\n",pbR->pbNum, Ts/100,s) ;
                            return 0 ;
                        }
                        T_values[ind] = Ts % 100 ; // on ne stocke que le suffixe
                    }
                }
            }
            
        }
    }
    Smin = 9999 * PB061_NS;
    {
        int inds[PB061_NS] ;
        u_int8_t permTyp[PB061_NS] ;
        T_Polygonal T[PB061_NS] ;
        T_Polygonal  pref0 = 0 ;
        permTyp[0] = PB061_NS-1 ;
        for(s=1;s<PB061_NS;s++) permTyp[s] = s-1 ;
        inds[0] = 0 ;
        for(s=0;s>=0;){
            int isOKs = 0 ;
            if((T[s] = T_values[inds[s]++])) {
                if(s==0) {
                    pref0 = T[s] / 100 ;
                    T[s] = T[s] % 100 ;
                }
                isOKs = 1;
            }
            if(isOKs) {
                if(s < PB061_NS-2) {
                    s++ ;
                    inds[s] = IND(T[s-1],permTyp[s]) ;
                    
                    continue ;
                } else { // fin de boucle on verifie le rebouclage
                    int indf ;
                    u_int8_t typf = permTyp[PB061_NS-1] ;
                    for(indf=IND(T[s],typf);T_values[indf] ; indf++) {
                        if(T_values[indf] == pref0) {
                            T[PB061_NS-1] = pref0 ;
                            {
                                int32_t S =0 ; int i ;
                                T_Polygonal Tc[PB061_NS] ;
                                for(i=0;i<PB061_NS;i++) {
                                    Tc[i] = (i==0) ? (100*pref0+T[0]) : (100*T[i-1]+T[i]) ;
                                    S += Tc[i] ;
                                }
                                nbSol++ ;
                                if(S < Smin) {
                                    Smin = S ;
                                    if(pbR->isVerbose){
                                        fprintf(stdout,"\t PB%0.3d S=%d,Nb=%d ",pbR->pbNum,S,nbSol);
                                        if(PB061_NS>6) {
                                            for(i=0;i<PB061_NS;i++) fprintf(stdout,"%d(%d)%c",permTyp[i],Tc[i],(i==PB061_NS-1) ? '\n' :' ') ;
                                        }else {
                                            for(i=0;i<PB061_NS;i++) fprintf(stdout,"%s(%d)%c",nameT[permTyp[i]],Tc[i],(i==PB061_NS-1) ? '\n' :' ') ;
                                            
                                        }
                                        
                                    }
                                }
                            }
                        }
                    }
                    isOKs= 0 ;
                }
            }
            if(!isOKs) { // on doit decrementer s
                if(s==0 || (isOnlyFirst && nbSol)) {
                    break ;
                }
                int i ;
                u_int8_t tmp ;
                for(i=s+1;i<PB061_NS;i++) {
                    if(permTyp[i] > permTyp[s]) {
                        tmp = permTyp[s] ;
                        permTyp[s] = permTyp[i] ;
                        permTyp[i] = tmp ;
                        inds[s] = IND(T[s-1],permTyp[s]) ;
                        break ;
                    }
                }
                if(i==PB061_NS) {
                    tmp = permTyp[s] ;
                    for(i=s+1;i<PB061_NS;i++) {
                        permTyp[i-1] = permTyp[i] ;
                    }
                    permTyp[PB061_NS-1] = tmp ;
                    s-- ;
                }
                
            }
        }
    }
    free(T_values);
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3d S=%d,Nb=%d\n",pbR->pbNum,Smin,nbSol);
    }
    if(nbSol > 0) sprintf(pbR->strRes,"%d",Smin);
    pbR->nbClock = clock() - pbR->nbClock ;
    return (nbSol > 0) ;
}

static u_int8_t * digCube = NULL ;

int CmpIndexPB62(const void *pt1,const void *pt2) {
    int cmp ;
    const u_int8_t * st1 = digCube + ((int32_t *)pt1)[0] ;
    const u_int8_t * st2 = digCube + ((int32_t *)pt2)[0] ;
    cmp = strcmp((char *)st1,(char *)st2) ;
    if(cmp != 0) { // pou les garder dans l'ordre initial en cas d'egalite
        return cmp ;
    } else return (int) (st1 - st2) ;
}

#define PB062_NBASK 5
#define PB062_MAX   100000
#define PB062_MAXD  18+1 // 3x5
// (racine cubique de 10 x 100000)+1
#define PB062_S1     215443
// (racine cubique de 100 x 100000)+1
#define PB062_S2     464157
#define PB062_S3     999999
#define PB062_SMIN  6   // 100
#define PB062_SMAX  12  // 10000

int PB062(PB_RESULT *pbR) {
    
    u_int64_t n , cube , bestCub=0 ;
    int i,pow, nbDig, nbEqualMax = 1 ;
    int32_t * index = NULL ;
    pbR->nbClock = clock() ;
    u_int32_t seuil[] = {
        0   , 3     ,  5
        , 10    , 22    , 47
        ,100    ,216    ,465
        ,1000   ,2155   ,4642
        ,10000  ,21545  ,46416
        ,100000 } ;
    {
        for(i=1,pow=100000;pow>1;pow /= 10) {
            seuil[i++] = (PB062_S1 / pow) + 1;
            seuil[i++] = (PB062_S2 / pow) + 1;
            seuil[i++] = (PB062_S3 / pow) + 1 ;
        }
    }
    for(nbDig=0,i=0;i<seuil[PB062_SMIN];i++) {
        if(i == seuil[nbDig]) nbDig++ ;
    }
    for(n=seuil[PB062_SMIN], cube=1000000 ;n<seuil[PB062_SMAX];n++) {
        u_int8_t * ptDig ;
        if(n==seuil[nbDig]) {
            nbDig++ ;
            digCube = malloc((nbDig+1)* (seuil[nbDig]-seuil[nbDig-1])) ;
            index = malloc((seuil[nbDig]-seuil[nbDig-1])*sizeof(index[0])) ;
        }
        index[n - seuil[nbDig-1]] = (int32_t) (n - seuil[nbDig-1])* (nbDig+1) ;
        ptDig = digCube + (n - seuil[nbDig-1]) * (nbDig+1) ;
        sprintf((char *)ptDig,"%lld",cube);
        HeapSortUint8(ptDig,nbDig);
        if(n+1==seuil[nbDig]) {
            u_int8_t *ant,*cur ;
            int nbSim = 1;
            qsort(index,seuil[nbDig]-seuil[nbDig-1],sizeof(index[0]),CmpIndexPB62);
            for(i=1,cur=digCube+index[0];i<seuil[nbDig]-seuil[nbDig-1];i++) {
                ant = cur ;
                cur = digCube+index[i] ;
                if(strcmp((char *)ant,(char *)cur)== 0) {
                    nbSim++;
                }else {
                    if(nbSim>nbEqualMax) {
                        int k ;
                        nbEqualMax = nbSim ;
                        bestCub = seuil[nbDig-1]+index[i-nbSim]/(nbDig+1) ;
                        bestCub = bestCub *bestCub * bestCub ;
                        if(pbR->isVerbose){
                            fprintf(stdout,"\t PB%d %d %lld ",pbR->pbNum, nbEqualMax,bestCub);
                            for(k=0;k<nbSim;k++)fprintf(stdout,"%d%c",seuil[nbDig-1]+index[i-nbSim+k]/(nbDig+1) ,(k==nbSim-1) ? '\n' : ' ' );
                        }
                    }
                    nbSim = 1 ;
                }
                
            }
            free(digCube) ;
            if(nbEqualMax == PB062_NBASK) break ;
        }
        cube += 3*n*(n+1)+1 ;
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    if(nbEqualMax == PB062_NBASK){
        sprintf(pbR->strRes,"%lld",bestCub) ;
        return 1 ;
    } else   return 0 ;
}


// si  10**(n-1) <= a**n < 10 **n
// en passant au log10 :  n-1 <= log10(a) < n
// en divisant par n : 1-1/n <= log10(a) < 1
// donc : a < 10 et 1 -log10(a) <= 1/n <=> 1 / (1 - log10(a)) >= n
int PB063(PB_RESULT *pbR) {
    u_int16_t nb = 0 ;
    int a;
    pbR->nbClock = clock()  ;
    for(a=1;a<10; a++) {
        nb += (int) ( 1 / ( 1 - log10((double)a) ) ) ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nb);
    return 1 ;
}


int PB064(PB_RESULT *pbR) {
    int32_t N , k0, k02, n , d ;
    int NbImpair = 0 ;
    pbR->nbClock = clock()  ;
    for(N=2,k0=1,k02=4;N<10000;N++) {
        int i ;
        if(N == k02) { // k02 = (k0+1)*(k0+1)
            k0++ ;
            k02 += 2*k0 + 1 ; continue ;
        }
        n = k0 ; d=1 ; i = 0 ; // so k0 =(int) srqt(N)
        do {
            d = (N - n * n) / d ; // quite easy to recursively show that division is exact
            n = k0 - ( (k0 + n) % d ) ;
            i++ ;
        } while(d !=1 || n != k0) ; // test loop on (n,d) = (k0,1)first couple
        if(i&1) NbImpair++ ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",NbImpair);
    return 1 ;
}

#define PB065_NITER         100
#define PB065_MAXDIGIT  (PB065_NITER+10)

#define PB065_N0    1
#define PB065_N1    2

#define PB065_D0    0
#define PB065_D1    1

typedef struct FractCont16 {
    u_int16_t *pt1 ;
    u_int16_t *pt ;
    int nbDig ;
}  FractCont16 ;
int PB065(PB_RESULT *pbR) {
    int n ;
    pbR->nbClock = clock()  ;
    u_int16_t DigN1[PB065_MAXDIGIT] ;
    u_int16_t DigN0[PB065_MAXDIGIT] ;
    u_int16_t DigD1[PB065_MAXDIGIT] ;
    u_int16_t DigD0[PB065_MAXDIGIT] ;
    FractCont16 FC[2] ; // numerateur indice 0, denominateur indice 1
    FC[0].nbDig = 1 ;
    FC[0].pt = DigN0 ;
    FC[0].pt1 = DigN1 ;
    FC[0].pt[0] = PB065_N1 ;
    FC[0].pt1[0] = PB065_N0 ;
    
    FC[1].nbDig = 1 ;
    FC[1].pt = DigD0 ;
    FC[1].pt1 = DigD1 ;
    FC[1].pt[0] = PB065_D1 ;
    FC[1].pt1[0] = PB065_D0 ;
    for(n=2;n<=PB065_NITER;n++) {
        int i ;
        int k ;
        int a = (n % 3) ? 1 : 2*(n/3) ;
        for(k=0;k<2;k++) { // numerateur et denominateur
            u_int16_t *tmp = FC[k].pt1 ;
            // pt->N(n) pt1->N(n->1)
            FC[k].pt1  = FC[k].pt ; FC[k].pt = tmp ;
            // pt1->N(n) pt->N(n-1) => N(n+1)
            // N(n+1) = a*N(n) + N(n-1)
            int carry = 0 ;
            for(i=0;i<FC[k].nbDig;i++) {
                FC[k].pt[i] += a * FC[k].pt1[i] + carry ;
                if(FC[k].pt[i] >= 10) {
                    carry = FC[k].pt[i] / 10 ;
                    FC[k].pt[i] = FC[k].pt[i] % 10 ;
                } else {
                    carry = 0 ;
                }
            }
            while(carry) {
                FC[k].pt1[FC[k].nbDig] = 0 ;
                FC[k].pt[FC[k].nbDig++] = carry % 10 ;
                carry /= 10 ;
            }
        }
    }
    {
        int k,SumD = 0;
        for(k=0;k<FC[0].nbDig;k++) SumD += FC[0].pt[k] ;
        sprintf(pbR->strRes,"%d",SumD);
        
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#include "p067_data.h"
#define PB067_SIZE  100
#define PB067_TSIZE ((PB067_SIZE*(PB067_SIZE+1))/2)
int PB067(PB_RESULT *pbR) {
    
    int ic,ir ;
    pbR->nbClock = clock() ;
    const u_int8_t * p067_data = P067_GetData() ;
    int vals[PB067_TSIZE] ;
    /* on recopie le tableau */
    {
        for(ic=0; ic <PB067_TSIZE;ic++) vals[ic] = p067_data[ic] ;
    }
    // on commernce a l'avant derniere ligne
    for(ir=PB067_SIZE-2;ir>=0;ir--) {
        int ic0 = (ir*(ir+1))/2 ;
        for(ic=0;ic<=ir;ic++) {
            int icnr = ic0+ir+1+ic ;
            vals[ic0+ic] += (vals[icnr] > vals[icnr+1]) ? vals[icnr] : vals[icnr+1] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",vals[0]) ;
    return 1 ;
}

// pour avoir 16 digits il faut que le 10 soit exterieur (sinon compte 2 fois donc 17 digits)
// comme on commence par le plus petit exterieur les chiffres exterieur sont 6,7,8,9,10
// et donc les chiffres interieurs 1,2,3,4,5
// Le total = sum(ext) +2 sum(int) = 40 + 2x15 = 70, donc le total sur une ligne est 70/5=14
// Pour completer le 6 seul couple possible 3,5.
// Pour completer le 10 seul couple possible 1,3
// Donc 3 partege entre 6 et 10, donc 6 et 10 sont voisin.
// 6,5,3 > 6,3,5 donc le 10 est le voisin apres le 6. Et l'on a place 6,5,3,10,1
// Pour completer le 9, puisque le 3 est pris par 6 et 10 il ne reste que 1,4
// Comme le 1 est place , le 9 est apres le 10.
// Il ne reste a placer en int que le 2, qui conduit donc en externe a placer le 8 puis le 4.
// Les lignes sont donc :
// 6,5,3 , 10,3,1, 9,1,4, 8,4,2, 7,2,5

int PB068(PB_RESULT *pbR) {
    
    
    pbR->nbClock = clock() ;
    strcpy(pbR->strRes,"6531031914842725") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB069_MAXN  1000000
int PB069a(PB_RESULT *pbR) {
    int32_t N = PB069_MAXN ;
    int32_t *phi = malloc(N * sizeof(phi[0])) ;
    int i, nBest = 2 ;
    pbR->nbClock = clock() ;
    for(i=0;i<N;i++) phi[i]=i ;
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3da ",pbR->pbNum) ;
    }
    for(i=2;i<N;i++) {
        if(phi[i] == i) { // nouveau nombre premier
            int np ;
            for(np=i;np<N;np+= i) {
                phi[np] = phi[np]/i * (i-1) ;
            }
        }
        if(i*(u_int64_t)phi[nBest] > phi[i]*(u_int64_t)nBest) {
            if(pbR->isVerbose)fprintf(stdout,"(%d,%d)",i,phi[i]);
            nBest = i ;
        }
    }
    if(pbR->isVerbose){
        fprintf(stdout,"\n\t PB%0.3da best n/phi %.2f for %d\n",pbR->pbNum,((float) nBest)/phi[nBest] ,nBest) ;
    }
    
    sprintf(pbR->strRes,"%d",nBest);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


// comme n/phi(n) = prod (Pi/(Pi-1) ) = prod(1/ (1 -1/Pi))
// les plus grande valeurs sont pour le maximum de Pi, et pour chacun la plus petite valeur
// Donc on cherche 2x3x5x7x11... xPi < PB069_MAXN
int PB069b(PB_RESULT *pbR) {
    int32_t N = PB069_MAXN ;
    int32_t P = 1 ;
    int32_t phi = 1 ;
    pbR->nbClock = clock() ;
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3db ",pbR->pbNum) ;
    }
    // ss optimiser les nbre premiers
    int p ;
    for(p=2;p*P < N ;p++) {
        int d,d2,isPrime ;
        for(d=2,d2=4,isPrime=1; d2<=p;d++) {
            if((p % d )== 0) { isPrime = 0 ; break ; }
            d2 += 2*d+1 ;
        }
        if(isPrime) { P *= p ; phi *= p-1 ; }
    }
    
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3db best n/phi %.2f for %d\n",pbR->pbNum,((float) P)/phi ,P) ;
    }
    
    sprintf(pbR->strRes,"%d",P);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB070_MAXN  10000000
int PB070(PB_RESULT *pbR) {
    int32_t N = PB070_MAXN ;
    int32_t *phi = malloc(N * sizeof(phi[0])) ;
    int i, nBest = 2 ;
    pbR->nbClock = clock() ;
    for(i=0;i<N;i++) phi[i]=i ;
    for(i=2;i<N;i++) {
        if(phi[i] == i) { // nouveau nombre premier
            int np ;
            for(np=i;np<N;np+= i) {
                phi[np] = phi[np]/i * (i-1) ;
            }
        }
        if(i*(u_int64_t)phi[nBest] < phi[i]*(u_int64_t)nBest) {
            unsigned char str_i[10], str_phi[10] ;
            int lg = sprintf((char *)str_i,"%d",i) ;
            HeapSortUint8(str_i,lg);
            lg=sprintf((char *) str_phi,"%d",phi[i]);
            HeapSortUint8(str_phi,lg);
            if(strcmp((char *)str_i,(char *)str_phi)== 0) {
                if(pbR->isVerbose)fprintf(stdout,"(%d,%d)",i,phi[i]);
                nBest = i ;
            }
        }
    }
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3d best n/phi %.6f for %d phi=%d\n",pbR->pbNum,((float) nBest)/phi[nBest] ,nBest,phi[nBest]) ;
    }
    
    sprintf(pbR->strRes,"%d",nBest);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
//
// on utilise le fait que la solution est
// Pi*Pj < N  avec (1-1/Pi)*(1-1/Pj) maximum <=> Pi/(Pi-1)*Pj/(Pj-1) minimum
// La solution est proche de sqrt(N) donc on filtre mieux en partant de la
int PB070a(PB_RESULT *pbR) {
    int32_t N = PB070_MAXN ;
    int32_t nSqrt = (int32_t) Sqrt64(N) ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(N/2+1)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    int i,j,jmax, nBest = 2 , phiBest = 1 ;
    pbR->nbClock = clock() ;
    u_int32_t nbPrime = GetNbPrime(ctxP) ;
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    for(i=0;i<nbPrime && (tbPrime[i] < nSqrt) ;i++) ;
    // on continue la boucle sur les Pi (croissants) et jmax est la valeur maxi pour que Pi*Pj < N)
    for(jmax=i-1;i<nbPrime;i++){
        int Pi = tbPrime[i] ;
        int Pj ;
        int n,phi ;
        for(j=jmax;j>=0; j--) {
            if( (n=Pi*(Pj=tbPrime[j])) >= N) { jmax-- ; continue ; }
            if(n*(u_int64_t) phiBest < (phi=(Pi-1)*(Pj-1)) * (u_int64_t) nBest) {
                unsigned char str_n[10], str_phi[10] ;
                int lg = sprintf((char *)str_n,"%d",Pi*Pj) ;
                HeapSortUint8(str_n,lg);
                lg=sprintf((char *) str_phi,"%d",phi);
                HeapSortUint8(str_phi,lg);
                if(strcmp((char *)str_n,(char *)str_phi)== 0) {
                    nBest = n ;
                    phiBest = (Pi-1)*(Pj-1);
                    //                    if(pbR->isVerbose)fprintf(stdout,"(%d,%d)",nBest,phiBest);
                }
            } else {
                // a Pi fixe, on ne peut que degrader en decroissant j
                break ;
            }
        }
    }
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3d best n/phi %.6f for %d phi=%d\n",pbR->pbNum,((float) nBest)/phiBest ,nBest,phiBest) ;
    }
    
    sprintf(pbR->strRes,"%d",nBest);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB070_MAX 100000
// facile de voir que l'on a n(k)*d(k-1)-d(k)*n(k-1) = 1
// donc si n(k) = 3 et d(k) =7
// 3 d(k-1) - 7 n(k-1) = 1
// 3x5-7x2=1
// donc d(k-1) = 5 + 7xi et n(k-1) = 2 + 3xi
// i = (1000000 - 5) / 7
// et n = 2 + 3*i
int PB071(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i = (1000000 - 5) / 7 ;
    int d = 5 + 7 * i ;
    int n = 2 + 3 * i ;
    if(pbR->isVerbose){
        fprintf(stdout,"\t PB%0.3d The fraction before 3/7 is %d/%d\n",pbR->pbNum,n,d) ;
    }
    sprintf(pbR->strRes,"%d",n);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB072_MAXN  1000000
int PB072(PB_RESULT *pbR) {
    int32_t N = PB072_MAXN ;
    int32_t *phi = malloc((N +1)* sizeof(phi[0])) ;
    u_int64_t S = 0 ;
    int i ;
    pbR->nbClock = clock() ;
    for(i=0;i<=N;i++) phi[i]=i ;
    for(i=2;i<=N;i++) {
        if(phi[i] == i) { // nouveau nombre premier
            int np ;
            for(np=i;np<=N;np+= i) {
                phi[np] = phi[np]/i * (i-1) ;
            }
        }
        S += phi[i] ;
    }
    
    sprintf(pbR->strRes,"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB073_MAXN  12000

#define PB073_N     1
#define PB073_D     3

#define PB073_Nend     1
#define PB073_Dend     2


int PB073(PB_RESULT *pbR) {
    int32_t N = PB073_MAXN ;
    int nb = 0 ;
    int d ;
    int n ;
    int d_end=PB073_Dend ;
    int n_end=PB073_Nend ;
    // on va chercher d0 et n0 tel que
    // on ait besout n x d0 - d * n0 = 1
    int d0, n0 ;
    d=PB073_D ;
    n=PB073_N ;
    { // solve besout
        int s0 = 1, s1 = 0;
        int t0 = 0, t1 = -1 ;
        int n1 = n ;
        int d1 = d ;
        do {
            int q = d1 / n1 ;
            int tmp = d1 - q * n1 ;
            d1 = n1 ;  n1 = tmp ;
            
            tmp = s0 + q * s1 ;
            s0 = s1 ; s1 = tmp ;
            
            tmp = t0  + q * t1 ;
            t0 = t1 ; t1 = tmp ;
            
        } while ( n1 ) ;
        d0 = -t0 ;
        n0 = s0 ;
        if(n*d0-d*n0 == -1) { // on inverse le signe
            int q = n0/n+1 ;
            n0 = -n0 + q*n;
            d0 = -d0 + q*d ;
        }
    }
    
    //    int d0=4 , n0 = 1 ; // satisfait besout n x d0 - d * n0 = 1
    do {
        int a = (N+d0)/d ; // on cherche d = a * d - d0 le plus grand possible
        int tmp = d ;
        d = a * d - d0 ;
        d0 = tmp ;
        tmp = n ;
        n = a * n - n0 ; // n = a * n - n0 ;
        n0 = tmp ;
        nb++ ;
        // on garde le couple (n,d) comme (n0,d0) car
        // besout  n x d0 - d * n0 = 1 est toujours satisfait
        // (a *n - n0) *d - (a*d - d0) * n = n x d0 - d * n0 = 1;
    } while(d != d_end && n != n_end ) ;
    
    
    
    pbR->nbClock = clock() ;
    
    sprintf(pbR->strRes,"%d",nb-1);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



typedef u_int32_t TY_FACT   ;


static TY_FACT FactDig[10] = { 1,1,2,6,24,120,720,5040,40320,362880 } ;

void FactVal(TY_FACT *valCur) {
    TY_FACT val = valCur[0] ;
    if(val == 169 ) {
        valCur[1] = FactDig[1] + FactDig[6] + FactDig[9] ; // 363601
        valCur[2] = 2 * FactDig[3] + 2 * FactDig[6] + FactDig[0] + FactDig[1] ;
        return ;
        
    }else if( val == 871 || val == 872) {
        valCur[1] = FactDig[8] + FactDig[7] + FactDig[val % 10] ;
        return ;
    } else {
        TY_FACT q;
        TY_FACT newVal , curval = val ;
        for(newVal=0, q=curval/10 ;curval ;curval = q, q=curval /10 ){
            newVal += FactDig[curval - 10 * q] ;
        }
        if(newVal == val) {
            return ;
        }
        valCur[1] = newVal ;
        FactVal(valCur+1);
    }
}

// ne retourne que la longueur, pas les valeurs
#define FACT_MEM    1
#if FACT_MEM
#define FACT_MEM_SIZE (1<<19)
static u_int16_t lgF[FACT_MEM_SIZE] ;
static u_int16_t privFactLgMem(u_int16_t lg,TY_FACT valInit,TY_FACT val) {
    if(val < FACT_MEM_SIZE ) {
        if(lgF[val]) {
            if(valInit < FACT_MEM_SIZE) lgF[valInit] = lg+lgF[val] ;
            return lg+lgF[val] ;
        } else if (val == 169) {
            lgF[val] = 3 ;
            if(valInit < FACT_MEM_SIZE) lgF[valInit] = lg+lgF[val] ;
            return lg+lgF[val] ;
        } else if ( val == 871 || val == 872) {
            lgF[val] = 2 ;
            if(valInit < FACT_MEM_SIZE) lgF[valInit] = lg+lgF[val] ;
            return lg+lgF[val] ;
        }
    }
    {
        TY_FACT q;
        TY_FACT newVal, valcur = val ;
        for(newVal=0  ;valcur > 9 ;valcur = q ){
            q=valcur/10 ;
            newVal += FactDig[valcur - 10 * q] ;
        }
        newVal +=FactDig[valcur] ;
        if(newVal == val) {
            if(val < FACT_MEM_SIZE ) lgF[val] = 1 ;
            return lg+1 ;
        }
        return privFactLgMem(lg+1,valInit,newVal);
    }
    
}
#else

static u_int16_t privFactLgH3(u_int16_t lg,TY_FACT val2,TY_FACT val1,TY_FACT val) {
    TY_FACT q;
    TY_FACT newVal, valcur = val ;
    for(newVal=0 ;valcur > 9 ;valcur = q ){
        q=valcur/10 ;
        newVal += FactDig[valcur - 10 * q] ;
    }
    newVal += FactDig[valcur] ;
    if(val == newVal || newVal == val1 || newVal == val2 ) {
        return lg+1 ;
    }
    return privFactLgH3(lg+1,val1,val,newVal);
}

static u_int16_t privFactLgH2(u_int16_t lg,TY_FACT val1,TY_FACT val) {
    TY_FACT q;
    TY_FACT newVal, valcur = val ;
    if(val == 169) {
        return lg +3 ;
    }
    for(newVal=0 ;valcur > 9 ;valcur = q ){
        q=valcur/10 ;
        newVal += FactDig[valcur - 10 * q] ;
    }
    newVal += FactDig[valcur] ;
    if(val == newVal || newVal == val1  ) {
        return lg+1 ;
    }
    return privFactLgH2(lg+1,val,newVal);
}

#endif
static u_int16_t privFactLg(u_int16_t lg,TY_FACT val) {
    if(val == 169 ) {
        return lg + 3 ;
    } else if ( val == 871 || val == 872) {
        return lg + 2 ;
    } else {
        TY_FACT q;
        TY_FACT newVal, valcur = val ;
        for(newVal=0 ;valcur > 9 ;valcur = q ){
            q=valcur/10 ;
            newVal += FactDig[valcur - 10 * q] ;
        }
        newVal += FactDig[valcur] ;
        if(val == newVal) {
            return lg+1 ;
        }
        return privFactLg(lg+1,newVal);
    }
}


static u_int16_t FactLg(u_int32_t val) {
#if FACT_MEM
    return privFactLgMem(0,val,val) ;
#else
    return privFactLg(0,val) ;
    //    return privFactLgH3(0,0,0,val) ;
    //    return privFactLgH2(0,0,val) ;
#endif
}

#define PB074_MAX_VALUE 1000000
#define PB074_PRINT 0
#define PB074_LGCOUNT   60
int PB074(PB_RESULT *pbR) {
    TY_FACT k ;
    TY_FACT kBest = 1 ;
    int bestLg = 0 ;
    int nbCount = 0 ;
    pbR->nbClock = clock() ;
    for(k=1;k<PB074_MAX_VALUE;k++) {
        u_int16_t lg = FactLg(k) ;
        if(lg > bestLg ) {
            bestLg = lg ;
            kBest = k ;
        }
        if(lg==PB074_LGCOUNT ) {
            nbCount++ ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d  lg=%d seen %d first for %d\n"
                              ,pbR->pbNum
                              ,bestLg,nbCount,kBest);
    sprintf(pbR->strRes,"%d",nbCount);
#if PB074_PRINT
    {
        TY_FACT *chainVal = malloc((bestLg+1)*sizeof(chainVal[0])) ;
        chainVal[0] = kBest ;
        FactVal(chainVal);
        
        {
            int i ;
            for(i=0 ; i<bestLg ; i++) {
                printf("%d%s",chainVal[i], (i==bestLg-1) ? "\n" : "->" );
            }
        }
        free(chainVal);
    }
#endif
    return 1 ;
}

//
// methode differente on remarque que l'on peut permuter les digits pour obtenir la meme longueur
// donc on cherche le calcul des longueurs uniquement pour le representant a digit croissant
int PB074a(PB_RESULT *pbR) {
    int nbCount = 0 ;
    pbR->nbClock = clock() ;
    int d1,d2,d3,d4,d5,d6 ;
    int dhist[10] = {0,0,0,0,0,0,0,0,0,0} ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3da ",pbR->pbNum);
    for(d1=0;d1<=9;d1++) {
        dhist[d1]++ ;
        for(d2=d1;d2<=9;d2++) {
            dhist[d2]++ ;
            for(d3=d2;d3<=9;d3++) {
                dhist[d3]++ ;
                for(d4=d3;d4<=9;d4++) {
                    dhist[d4]++ ;
                    for(d5=d4;d5<=9;d5++) {
                        dhist[d5]++ ;
                        for(d6=d5;d6<=9;d6++) {
                            dhist[d6]++ ;
                            int k = d1*100000+d2*10000+d3*1000+d4*100+d5*10+d6 ;
                            if(privFactLg(0,k) == PB074_LGCOUNT) {
                                int i,npermut = 1;
                                for(i=2;i< 6-dhist[0];i++) { // (nbdig-1)!
                                    npermut *= i ;
                                }
                                if(dhist[1]) { // rempacement des 1 par des zeros sauf le 1 leader
                                    npermut = (npermut * (1 << (dhist[1]-1)) )* (12-2*dhist[0]-dhist[1]) ;
                                } else {
                                    npermut *= (6 -dhist[0])  ;
                                }
                                for(i=2;i<=9;i++) {
                                    if(dhist[i] > 1) npermut /= FactDig[dhist[i]] ;
                                }
                                nbCount += npermut ;
                                if(pbR->isVerbose)fprintf(stdout," +Permut(%d)=%d ",k,npermut);
                            }
                            dhist[d6]-- ;
                        }
                        dhist[d5]-- ;
                    }
                    dhist[d4]-- ;
                }
                dhist[d3]-- ;
            }
            dhist[d2]-- ;
        }
        dhist[d1]-- ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d seen %d \n"
                              ,pbR->pbNum,nbCount);
    sprintf(pbR->strRes,"%d",nbCount);
#if PB074_PRINT
    {
        TY_FACT *chainVal = malloc((bestLg+1)*sizeof(chainVal[0])) ;
        chainVal[0] = kBest ;
        FactVal(chainVal);
        
        {
            int i ;
            for(i=0 ; i<bestLg ; i++) {
                printf("%d%s",chainVal[i], (i==bestLg-1) ? "\n" : "->" );
            }
        }
        free(chainVal);
    }
#endif
    return 1 ;
}


#define PB075_MAXN  1500000
// pour m > n > 0
// a = (m**2 - n**2)*k ; b = 2*m*n*k ; c = (m**2 + n**2)*k ;
// donc S/2 = (m**2 + m*n)*k = m * (m+n) * k
int PB075(PB_RESULT *pbR) {
    int32_t N = PB075_MAXN/2 ;
    pbR->nbClock = clock() ;
    int32_t m,mSqrt = (int32_t) Sqrt64(N);
    int32_t * nbPyth = calloc(N+1 ,sizeof(nbPyth[0])) ;
    u_int32_t nb = 0 ;
    
    for(m=2;m<mSqrt;m++) {
        int32_t n ;
        int32_t L = m*m ;
        int32_t H = L ;
        for(n=1;n<m;n++) {
            L += m ;
            H += 2*n - 1 ;
            if( ((n&1) ^(m& 1)) == 0) continue ; // on saute si parite identique (pour l'unicite)
            {
                int kH,kL ;
                for(kL=L,kH=H;kL <=N;kL+= L, kH += H) {
                    if(nbPyth[kL] == 0) {
                        nbPyth[kL] = kH  ;
                    } else if(nbPyth[kL] != kH ) {
                        nbPyth[kL] = -1 ;
                    }
                }
            }
        }
    }
    int i ;
    for(i=1;i<=N;i++) {
        if(nbPyth[i]>0) {
            nb++ ;
        }
    }
    sprintf(pbR->strRes,"%d",nb);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB076_MAXN  100
#define IND76(k,n) ((PB076_MAXN)*(n-1)+(k-1))
int PB076(PB_RESULT *pbR) {
    int32_t Ekn[(PB076_MAXN)*(PB076_MAXN)] ;
    int n ;
    pbR->nbClock = clock() ;
    
    Ekn[IND76(1,1)] = 1 ;
    for(n=2;n<=PB076_MAXN;n++) {
        int k ;
        Ekn[IND76(1,n)] = 1 ;
        for(k=2;k<=n/2;k++) {
            // decompostion de n en k nombre.
            // soit la decomposition contient au mois un 1, et donc en l'enlevant on tombe sur Ekn(k-1,n-1)
            // soit la decomposition ne contient pas de 1, et l'on peut enlever 1 a chaque element Ekn(k,n-k)
            // Le deuxieme cas ne peut se produire que si k <= n/2
            Ekn[IND76(k,n)] = Ekn[IND76(k-1,n-1)] + Ekn[IND76(k,n-k)] ;
        }
        for(;k<=n;k++) {
            Ekn[IND76(k,n)] = Ekn[IND76(k-1,n-1)]  ;
        }
    }
    int k,S ; // on saute k=1 car somme en au moins 2 elements demande
    for(k=2,S=0;k<=PB076_MAXN;k++) S += Ekn[IND76(k,PB076_MAXN)] ;
    sprintf(pbR->strRes,"%d",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

// base sur le theoreme pentagonal d'euler
// P(n) = Sum(k){ (-1)**(k+1) P(n - (k(3k-1))/2) }
int PB076a(PB_RESULT *pbR) {
    int32_t Pn[PB076_MAXN+1] ;
    int n ;
    pbR->nbClock = clock() ;
    Pn[0] = Pn[1] = 1 ;
    for(n=2;n<=PB076_MAXN;n++) {
        int k , P = 0 ,Pk = 1  ;
        for(k=1; Pk <= n ; Pk += 3*k+1, k++ ) {
            if(k&1) {
                P += Pn[n - Pk] ;
                // car P-k = Pk+k
                if(Pk + k <= n ) P += Pn[n - Pk - k ] ;
            } else {
                P -= Pn[n - Pk] ;
                if(Pk + k <= n ) P -= Pn[n - Pk - k ] ;
            }
        }
        Pn[n] = P ;
    }
    sprintf(pbR->strRes,"%d",Pn[PB076_MAXN]-1);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB077_MAXN  100
#define PB077_NB_ASK    5000
// calculate E(k,n) où k represente le kime nombre premier Pk  et 0 <= n < Maxn
// E(kn) represente le nombre de decomposition n'utilisant que des premiers <= Pk
// On calcule par recurrence sur k  ( E(0,n) vaut 1 si n pair 0 sinon )
// Ensuite E(k,n) = Sum (i=0, ...) E(k-1, n - i * Pk).
int PB077(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE *ctxP ;
    if((ctxP = Gen_tablePrime(PB077_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    u_int32_t nbPrime = GetNbPrime(ctxP) ;
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    
    int32_t *Epn = calloc(nbPrime*PB077_MAXN,sizeof(Epn[0])) ;
    int ip,k ;
    int p = 2 ;
    for(k=0;2*k<PB077_MAXN;k++) {
        Epn[2*k] = 1 ; Epn[2*k+1] = 0 ;
    }
    
    for(ip=1;ip<nbPrime;ip++) {
        p = tbPrime[ip] ;
        int np ;
        for(np=0; np < PB077_MAXN; np+=p) {
            for(k=np;k < PB077_MAXN;k++) {
                Epn[ip * PB077_MAXN + k] += Epn[(ip-1)*PB077_MAXN+k-np] ;
            }
        }
    }
    int kmin ;
    for(kmin=0;kmin<PB077_MAXN;kmin++) {
        //        printf("%d ",Epn[(ctxP->nbPrime-1)*PB077_MAXN + kmin]) ;
        if(Epn[(nbPrime-1)*PB077_MAXN + kmin]-1 > PB077_NB_ASK) {
            if(pbR->isVerbose) fprintf(stdout,"\t PB%0.3d %d a %d decompositions en premiers\n"
                                       ,pbR->pbNum,kmin,Epn[(nbPrime-1)*PB077_MAXN + kmin]-1);
            break ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    free(Epn) ;
    Free_tablePrime(ctxP) ;
    if(kmin<PB077_MAXN){
        sprintf(pbR->strRes,"%d",kmin);
        return 1 ;
    } else {
        return 0 ;
    }
}



// voir PB076a
// meme calcul avec la valeur exacte de Pn[n] (sans -1)
#define PB078_MAX   100000
#define PB078_DIV   1000000
int PB078(PB_RESULT *pbR) {
    int32_t * Pn= malloc((PB078_MAX+1) *sizeof(Pn[0])) ;
    int n ;
    pbR->nbClock = clock() ;
    Pn[0] = Pn[1] = 1 ;
    for(n=2;n<=PB078_MAX;n++) {
        int k ;
        int32_t Pk = 1  ;
        int32_t P = 0 ;
        for(k=1; Pk <= n ; Pk += 3*k+1, k++ ) {
            if(k&1) {
                P += Pn[n - Pk] ;
                // car P-k = Pk+k
                if(Pk + k <= n ) P += Pn[n - Pk - k ] ;
            } else {
                P -= Pn[n - Pk] ;
                if(Pk + k <= n ) P -= Pn[n - Pk - k ] ;
            }
        }
        P =  P % PB078_DIV ;
        //        if (P < 0) P += PB078_DIV ;
        Pn[n] =  P ;
        if(P == 0) break ;
    }
    sprintf(pbR->strRes,"%d",n);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#include "p081_data.h"



#define MATDYN(i,j)  ((i)*(sizeM+2)+(j))
#define MATDYN_PC  1000000000
#define MATDYN_PD  1000000000

int ProgDyn_Mat(int sizeM,const u_int16_t * cout,const int *startCase, const int *endCase,const int *Neighbour ) {
    u_int32_t   *antDist = malloc((sizeM+2)*(sizeM+2)*sizeof(antDist[0])) ;
    u_int32_t   *newDist = malloc((sizeM+2)*(sizeM+2)*sizeof(newDist[0])) ;
    
    u_int32_t   *Cout = malloc((sizeM+2)*(sizeM+2)*sizeof(Cout[0])) ; ;
    
    int i,j ;
    for(j=0;j<=sizeM+1;j++) {
        Cout[MATDYN(0,j)] = Cout[MATDYN(sizeM+1,j)] = MATDYN_PC ;
        newDist[MATDYN(0,j)] = newDist[MATDYN(sizeM+1,j)] = MATDYN_PD ;
        antDist[MATDYN(0,j)] = antDist[MATDYN(sizeM+1,j)] = MATDYN_PD ;
    }
    for(i=1;i<=sizeM;i++) {
        Cout[MATDYN(i,0)] = Cout[MATDYN(i,sizeM+1)] = MATDYN_PC ;
        newDist[MATDYN(i,0)] = newDist[MATDYN(i,sizeM+1)] = MATDYN_PD ;
        antDist[MATDYN(i,0)] = antDist[MATDYN(i,sizeM+1)] = MATDYN_PD ;
        for(j=1;j<=sizeM;j++) {
            Cout[MATDYN(i,j)] = cout[i*sizeM+j-sizeM-1] ;
            newDist[MATDYN(i,j)] = MATDYN_PD ;
        }
    }
    // debut
    while(*startCase) {
        newDist[*startCase] = Cout[*startCase] ;
        startCase++ ;
    }
    int isMod = 0;
    do {
        isMod = 0 ;
        u_int32_t   * tmp = antDist ;
        antDist = newDist ;
        newDist = tmp ;
        for(i=1;i<=sizeM;i++) {
            int ic = MATDYN(i,1) ;
            for(j=1;j<=sizeM;j++) {
                int k ;
                u_int32_t bestD ;
                int cout_ic = Cout[ic] ;
                bestD = antDist[ic] - cout_ic ;
                for(k=0;Neighbour[k];k++) {
                    if(antDist[ic+Neighbour[k]] < bestD) { bestD = antDist[ic+Neighbour[k]] ; isMod = 1 ; }
                }
                newDist[ic++] = bestD + cout_ic ;
            }
        }
    } while(isMod) ;
    u_int32_t minD = newDist[*endCase];
    
    while(*++endCase){
        if(newDist[*endCase] < minD) {
            minD = newDist[*endCase] ;
        }
    }
    free(Cout); free(antDist); free(newDist);
    return minD;
}


#define P081_SIZE  80

int PB081(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int neigh[] = { -1 , -P081_SIZE-2 , 0  } ;
    int start[] = { P081_SIZE+3 , 0 } ;
    int end[] = { (P081_SIZE+2)*P081_SIZE+ P081_SIZE  , 0 } ;
    const u_int16_t * p81_data = P081_GetData() ;
    int best = ProgDyn_Mat(P081_SIZE,p81_data,start,end,neigh) ;
    sprintf(pbR->strRes,"%d",best);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB082(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int neigh[] = { -1 , -P081_SIZE-2 , +P081_SIZE+2, 0  } ;
    int i ;
    int start[P081_SIZE+1] ;
    int end[P081_SIZE+1] ;
    const u_int16_t * p81_data = P081_GetData() ;
    for(i=1;i<=P081_SIZE;i++){
        start[i-1] = i*(P081_SIZE+2) + 1 ;
        end[i-1] = i*(P081_SIZE+2) + P081_SIZE ;
    }
    start[P081_SIZE] = end[P081_SIZE] = 0 ;
    int best = ProgDyn_Mat(P081_SIZE,p81_data,start,end,neigh) ;
    sprintf(pbR->strRes,"%d",best);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB083(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int neigh[] = { 1, -1 , -P081_SIZE-2 , +P081_SIZE+2, 0  } ;
    int start[] = { P081_SIZE+3 , 0 } ;
    int end[] = { (P081_SIZE+2)*P081_SIZE+ P081_SIZE  , 0 } ;
    const u_int16_t * p81_data = P081_GetData() ;
    int best = ProgDyn_Mat(P081_SIZE,p81_data,start,end,neigh) ;
    sprintf(pbR->strRes,"%d",best);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



// number of rectangles is NB(n,p) = n(n+1)/2 x p(p+1)/2
// search n,p nearest 2000000
// Loop on Tn and Tp so Tn*Tp <= 2000000 < Tn(Tp+1)
// use recursion to compute Tn*tp - 2000000

#define PB085_NB    2000000
int PB085(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n,p,Tn,Tp ;
    int nBest = 0 ;
    int pBest = 0 ;
    int deltaNBBest =  PB085_NB ;
    int D,D1 ;
    for(n=1,Tn=1,p=(int)Sqrt64(2*PB085_NB),Tp=(p*(p+1))/2,D1=0, D = Tp-PB085_NB; n<=p;n++,Tn += n , D += n * Tp ) {
        while(D > 0 ) {
            D1 = D ;
            D -= Tn * p ;
            Tp -= p-- ;
        }
        if(-D < deltaNBBest ) {
            deltaNBBest = - D ;
            nBest = n ;
            pBest = p ;
        }
        if (D1 < deltaNBBest ) {
            deltaNBBest = D1 ;
            nBest = n ;
            pBest = p+1 ;
        }
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%0.3d For %d,%d np=%d Nb rectangles=%d\n",pbR->pbNum,nBest,pBest,nBest*pBest,(nBest*(nBest+1)*pBest*(pBest+1))/4) ;
    sprintf(pbR->strRes,"%d",nBest*pBest);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


// en mettant le parallepipede a <= b <= c a plat deplie on voir rapidement
// que la plus coure distance au carre est (a+b)**2 + c**2
// donc a+b,c sont les petites valeurs d'un triangle pythagoricien
// donc on cherche a decomposer un triangle Pyth primaire (m,n) m<n premiers m**m - n**2 , 2mn
// on compte les cas : a+b = m**2 - n**2 , c =2mn et a<=b<=c et c<=M
// et les cas a+b =2mn c = m**2-n**2 et c <=M
// puis on mutliplie par k tant que l'on ne depasse pas M

#define PB086_MAX    5000
#define PB086_NBASK 1000000

int PB086(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i,n,m,c,k ;
    int M = PB086_MAX ;
    int minM = 1 ;
    u_int32_t *histoM = calloc(PB086_MAX,sizeof(histoM[0])) ;
    for(m=2;;m++) {
        // c >= 1/3 (a+b+c) as c>=b>=a
        // (a+b=c) = m**2 - n**2 + 2*m*n has a minima for n in [0,m] when n=0 => m**2
        // so minc <= m*m / 3
        int minC = (m*m)/3  ;
        for(n=(m&1) + 1;n<m;n += 2) {
            int apb,tmp,permut ;
            if(PGCD(m,n) != 1)  continue ;
            // first case a+b = m**2 - n**2 , c = 2*m*n
            // second case a+b = 2*m*n  c = m**2 - n**2
            c = 2 * m * n ;
            apb = m*m - n*n ;
            for(permut=0;permut<2;permut++) { // loop for permutation c <=> apb
                if (c<=M && apb <= 2*c) {
                    k = M / c ;
                    if( apb <= c ) {
                        for(i=1;i<=k;i++) histoM[i*c] += (i*apb)/2 ;
                    } else {
                        for(i=1;i<=k;i++) histoM[i*c] += i*c-(i*apb-1)/2 ;
                    }
                }
                tmp = apb ;
                apb = c ;
                c = tmp ;
            }
        }
        for(;minM+1<minC;) { // we know that histo under minc will remain unchanged
            // so we cand cumulate histogram
            minM++ ;
            histoM[minM] += histoM[minM-1] ;
            if(histoM[minM]>PB086_NBASK) {
                break ;
            }
        }
        if(histoM[minM]>PB086_NBASK) break ;
    }
    free(histoM) ;
    if(pbR->isVerbose)printf("\t PB%0.3d Nb[%d]=%d => Nb[%d]=%d\n",pbR->pbNum,minM-1,histoM[minM-1],minM,histoM[minM]) ;
    sprintf(pbR->strRes,"%d",minM);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB087_MAXN  50000000
int PB087(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE *ctxP ;
    u_int32_t n_sqr2 = (u_int32_t)Sqrt64(50000000) ;
    u_int32_t n_sqr4 = (u_int32_t)Sqrt64(n_sqr2+1) ;
    u_int32_t n_sqr3 = 800 ;
    while(n_sqr3*n_sqr3*n_sqr3 > PB087_MAXN) n_sqr3-- ;
    u_int32_t  * pow2 , * pow3 ,*pow4 ;
    pow2 = malloc(n_sqr2*sizeof(pow2[0])) ;
    pow3 = malloc(n_sqr3*sizeof(pow3[0])) ;
    pow4 = malloc(n_sqr4*sizeof(pow4[0])) ;
    if((ctxP = Gen_tablePrime(n_sqr2)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    
    int i,i2,i3,i4 ;
    u_int32_t nbPrime = GetNbPrime(ctxP) ;
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    
    for(i=i2=i3=i4=0;i<nbPrime;i++) {
        u_int32_t p,p2 ;
        p = tbPrime[i] ;
        pow2[i2++] = p2 = p*p ;
        if(p< n_sqr3) {
            pow3[i3++] = p*p2 ;
            if(p < n_sqr4) {
                pow4[i4++] = p2*p2 ;
            }
        }
    }
    n_sqr2 = i2 ;
    n_sqr3 = i3 ;
    n_sqr4 = i4 ;
    int nbLoop = 4;
    int nb = 0 ;
    int sizeLoop = (((PB087_MAXN/16+1) / nbLoop) +1) * 16 ;
    u_int16_t * isDec = calloc(sizeLoop/16,sizeof(isDec[0])) ;
    int nl ;
    for(nl=0;nl<nbLoop;nl++) {
        int infLoop = sizeLoop * nl ;
        for(i4=0;i4<n_sqr4;i4++) {
            int32_t S4 = pow4[i4] - infLoop ;
            if(S4 >= sizeLoop) break ;
            for(i3=0;i3<n_sqr3;i3++) {
                int32_t S3 = S4 + pow3[i3] ;
                if(S3 >= sizeLoop) break ;
                for(i2=0;i2<n_sqr2;i2++) {
                    int32_t S2 = S3 + pow2[i2] ;
                    if(S2 < 0) continue ;
                    if( S2 >= sizeLoop) break ;
                    if((isDec[S2>>4] & ( 1 << ( S2 & 0xf ))) ==0) {
                        nb++; isDec[S2>>4] |= 1 << ( S2 & 0xf ) ;
                    }
                }
            }
        }
        if(nl < nbLoop-1) memset(isDec,0, (sizeLoop/16) * sizeof(isDec[0]));
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nb) ;
    free(isDec); free(pow2); free(pow3); free(pow4) ;
    Free_tablePrime(ctxP) ;
    return 1 ;
}


#define PB088_MAXK  12000

static int CmpKminN(const void *el1,const void *el2) {
    return ((u_int32_t *)el1)[0] - ((u_int32_t *)el2)[0] ;
}
int PB088(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int32_t *kminN = calloc(PB088_MAXK+1,sizeof(kminN[0])) ;
    int nv ;
    u_int32_t V[32] /* values */ ,P[32] /* product */, S[32] /* sum */ ;
    int Pmax = 2 * PB088_MAXK ;
    int isNotAll2 = 1;
    for(nv=2;isNotAll2;nv++) {
        int iP , iPbad ;
        u_int32_t k = 1 ;
        V[nv] = 1;
        P[nv] = 1 ;
        S[nv] = - nv ;
        for(iPbad = 0, iP=nv-1, V[iP]=1 ; iP<nv ; iP++) {
            u_int32_t newV = V[iP] + 1;
            for(; iP>=0 ; iP--) {
                V[iP] = newV ;
                P[iP] = P[iP+1]*newV ;
                S[iP] = S[iP+1]+newV ;
                k = P[iP] - S[iP]  ;
                if(P[iP] > Pmax ||  k >  PB088_MAXK ) break ;
            }
            if(iP >= 0) {
                if(iP == 0 && V[iP] == 2) {
                    isNotAll2 = 0 ; break ;// END, only V==2 and too big
                }
                iP = iPbad++ ;
                continue ;
            } else {
                iPbad = 0 ; // find a new value
                if(kminN[k]) {
                    if(P[0] < kminN[k]) {
                        kminN[k] = P[0];
                    }
                } else {
                    kminN[k] = P[0] ;
                }
            }
            
        }
    }
    {   // compute the Sum
        int k ; u_int64_t S=0;
        qsort(kminN+1,PB088_MAXK,sizeof(kminN[0]),CmpKminN) ;
        u_int32_t ant = kminN[2] ;
        S += kminN[2] ;
        for(k=3;k<=PB088_MAXK;k++) {
            if(kminN[k] == ant) continue ;
            S += kminN[k] ;
            ant = kminN[k] ;
        }
        sprintf(pbR->strRes,"%lld",S) ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define isPresent(v0,v1)    ( (( D1 & (1<<v0)) && (D2 & (1<<v1)) ) || (( D2 & (1<<v0)) && (D1 & (1<<v1)) ) )
#define isPresent6(v0)     ( (( D1 & (1<<v0)) && ( D2 & ( (1<<6) + (1<<9) ) ) ) || (( D2 & (1<<v0)) && (D1 & ( (1<<6) + (1<<9) ) ) ) )
int PB090(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nb = 0 ;
    // as number of 6 bits between 10
    u_int16_t diffPosition[210] ;
    // by complementary place 4 zero's in ten position
    {
        int i = 0 ;
        int i1,i2,i3,i4 ;
        for(i1=9;i1>2;i1--) {
            for(i2=i1-1;i2>1;i2--) {
                for(i3=i2-1;i3>0;i3--) {
                    for(i4=i3-1;i4>=0;i4--) {
                        u_int16_t v = (1<<i1) + (1<<i2) + (1<<i3) + (1<<i4)  ;
                        diffPosition[i++] = 1023 ^ v ;
                    }
                }
            }
        }
    }
    int i1,i2 ;
    for(i1=0;i1<210;i1++) {
        u_int16_t D1 = diffPosition[i1] ;
        for(i2=i1;i2<210;i2++) {
            u_int16_t D2 = diffPosition[i2] ;
            if(!isPresent(0,1)) continue ;
            if(!isPresent(0,4)) continue ;
            if(!isPresent6(0)) continue ; // 09
            if(!isPresent6(1)) continue ; // 16
            if(!isPresent(2,5)) continue ;
            if(!isPresent6(3)) continue ; // 36
            if(!isPresent6(4)) continue ; // 64
            if(!isPresent(8,1)) continue ;
            nb++ ;
        }
    }
    
    sprintf(pbR->strRes,"%d",nb);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
#define PB091_MAXY  50

int PB091(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nbO=0, nbp=0,nbq = 0 ;
    int x1,y1,x2,y2 ;
    for(x1=0;x1<=PB091_MAXY;x1++) {
        for(y1=0;y1<=PB091_MAXY;y1++){
            if(x1+y1 == 0) continue ;
            for(y2=0;y2<y1;y2++) {
                for(x2=0;x2<=PB091_MAXY;x2++) {
                    if(x2+y2 == 0) continue ;
                    if(x1*(x2-x1)+y1*(y2-y1) == 0) { nbp++ ;}
                    else if(x1*x2+y1*y2 == 0) {  nbO++; }
                    else if(x2*(x1-x2)+y2*(y1-y2)== 0 ) { nbq++; }
                    
                }
            }
            for(y2=y1,x2=0;x2<x1;x2++) {
                if(x2+y2==0) continue ;
                if(x1*(x2-x1)+y1*(y2-y1) == 0) {  nbp++ ;}
                else if(x1*x2+y1*y2 == 0) {  nbO++; }
                else if(x2*(x1-x2)+y2*(y1-y2)== 0 ) {  nbq++; }
                
            }
            
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d Nb=%d (O:%d,P:%d,Q:%d)\n",pbR->pbNum,nbO+nbp+nbq, nbO,nbp,nbq) ;
    sprintf(pbR->strRes,"%d",nbO+nbp+nbq);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB092_NBDIG     19
#define PB092_NBASK     7
#define PB092_T89       2

int PB092(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int8_t dig[PB092_NBDIG] ;
    int16_t sum[PB092_NBDIG] ;
    int maxValue = 81*PB092_NBDIG ;
    int nbT[PB092_T89+1] ;
    nbT[1] = 0 ;
    nbT[PB092_T89] = 0 ;
    int i,n,lgBack  ;
    u_int8_t *terminal = calloc(maxValue+1,sizeof(terminal[0]));
    int16_t *backTrace = malloc((maxValue+1)*sizeof(backTrace[0])) ;
    terminal[1] = 1 ; nbT[1]++ ;
    terminal[89] = PB092_T89 ; nbT[PB092_T89]++ ;
    for(i=1;i<=maxValue;i++) {
        n = i; lgBack = 0 ;
        while(n>maxValue || terminal[n]==0 ) {
            if(n<=maxValue) backTrace[lgBack++] = n ;
            int nxt = 0 ;
            while(n>0) {
                nxt += (n % 10) * (n % 10) ;
                n /= 10 ;
            }
            n = nxt ;
        }
        int k ;
        for(k=0;k<lgBack;k++) {
            terminal[backTrace[k]] = terminal[n] ;
        }
        nbT[terminal[n]] += lgBack ;
    }
    // on parcours tous les nombres au dela de max value
    memset(sum,0,PB092_NBDIG*sizeof(sum[0])) ;
    memset(dig,0,PB092_NBDIG*sizeof(dig[0])) ;
    for(n=maxValue,i=PB092_NBDIG-1;n!=0;i--) {
        dig[i] = n % 10 ;
        sum[i] = dig[i]*dig[i] ;
        n /= 10 ;
    }
    while(++i<=PB092_NBDIG-1) { sum[i] += sum[i-1] ; }
    int is = PB092_NBDIG-1 ;
    while(is >=0) {
        if(dig[is] < 9) {
            sum[is] += 2*dig[is]+1 ;
            dig[is]++ ;
            while(is <PB092_NBDIG-1){
                is++ ;
                sum[is] = sum[is-1] ;
                dig[is] = 0 ;
            }
            nbT[terminal[sum[is]]]++ ;
        } else {
            is-- ;
        }
    }
    
    
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d Term=1[%d] Term=89[%d]\n",pbR->pbNum,nbT[1],nbT[PB092_T89]) ;
    sprintf(pbR->strRes,"%d",nbT[PB092_T89]);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB092a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int8_t dig[PB092_NBDIG] ;
    int16_t sum[PB092_NBDIG] ;
    int maxValue = 81*PB092_NBDIG ;
    u_int64_t nbT[PB092_T89+1] ;
    nbT[1] = 0 ;
    nbT[PB092_T89] = 0 ;
    int i,n,lgBack  ;
    u_int8_t *terminal = calloc(maxValue+1,sizeof(terminal[0]));
    int16_t *backTrace = malloc((maxValue+1)*sizeof(backTrace[0])) ;
    terminal[1] = 1 ; // nbT[1]++ ;
    terminal[89] = PB092_T89 ; // nbT[PB092_T89]++ ;
    for(i=1;i<=maxValue;i++) {
        n = i; lgBack = 0 ;
        while(n>maxValue || terminal[n]==0 ) {
            if(n<=maxValue) backTrace[lgBack++] = n ;
            int nxt ;
            for(nxt = 0 ; n > 0; n /= 10 ) {
                nxt += (n % 10) * (n % 10) ;
            }
            n = nxt ;
        }
        int k ;
        for(k=0;k<lgBack;k++) {
            terminal[backTrace[k]] = terminal[n] ;
        }
//        nbT[terminal[n]] += lgBack ;
    }
    {
        u_int64_t * histoSum[PB092_NBDIG] ,sum ;
        int i,nb,is ;
        for(nb=0;nb<PB092_NBDIG;nb++) histoSum[nb] = calloc((81+1)*(nb+1),sizeof(histoSum[nb][0])) ;
        for(i=0;i<10;i++) histoSum[0][i*i]++  ; // premier // digit
        for(nb=1;nb<PB092_NBDIG;nb++) {
            for(is=0;is<=81*nb;is++) {
                for(i=0;i<10;i++) {
                    if(histoSum[nb-1][is]) histoSum[nb][is+i*i] += histoSum[nb-1][is] ;
                }
            }
            nbT[1] = nbT[PB092_T89] = 0 ;
            for(is=1;is<=81*nb;is++) {
                nbT[terminal[is]] += histoSum[nb][is] ;
            }
            if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d %cTerm=89[%lld]\n",pbR->pbNum,(nb==PB092_NBASK-1) ? '*':' ',nbT[PB092_T89]) ;
            if(nb==PB092_NBASK-1 ) sprintf(pbR->strRes,"%lld",nbT[PB092_T89]); ;
        }
        for(i=0;i<PB092_NBDIG;i++) free(histoSum[i]) ;
    }
    
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



typedef struct Fract_093 {
    int32_t N ;
    int32_t D ;
} Fract_093 ;
Fract_093 Oper(int nop, Fract_093 d0, Fract_093 d1) {
    Fract_093 fr ;
        switch(nop){
            case 0: fr.D=d0.D*d1.D; fr.N=d0.N*d1.D+d1.N*d0.D;  break ; // +
            case 1: fr.D=d0.D*d1.D ; fr.N=d0.N*d1.N ; break ; // *
            case 2: fr.D=d0.D*d1.D; fr.N=d0.N*d1.D-d1.N*d0.D; if(fr.N<0) fr.N = -fr.N ; break ; // -
            case 3: fr.D=d0.D*d1.N ; fr.N=d0.N*d1.D ; break ; // d0/d1
            case 4: fr.D=d1.D*d0.N ; fr.N=d1.N*d0.D ; break ; // d1/d0
        }
       if(fr.D > 1) {
            int32_t gd = PGCD(fr.N,fr.D) ;
            if(gd > 1) { fr.N /= gd ; fr.D /= gd ; }
        }
    return fr ;
}


Fract_093 FractDig(u_int8_t dig) {
    Fract_093 fr ;
    fr.D = 1 ; fr.N = dig ;
    return fr ;
}

#define PB095_NBD       4
#define PB095_NBO       5
#define PB095_EXPMAX   100

// EXPMAX 1000 pour 5 digits ou  6 digits

void AddVal(u_int8_t *isValFind, Fract_093 res) {
    if(res.D != 1 || res.N == 0 || res.N > PB095_EXPMAX ) return ;
    isValFind[res.N-1] = 1 ;
}
typedef struct niv_093 {
    Fract_093 vals[PB095_NBD] ;
    int16_t op1 ;
    int16_t op2 ;
    int16_t oper ;
} niv_093 ;

// version avec parcours recursif
int PB093a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int8_t isValFind[PB095_EXPMAX] ;
    int maxCons =0 ;
    char maxABCD[PB095_NBD+1] ;
    maxABCD[PB095_NBD] = 0 ;
    niv_093 niv[PB095_NBD] ;
    int i , indDig ;
    for(i=0;i<PB095_NBD;i++) {
        niv[0].vals[i] = FractDig(i+1) ;
        niv[i].op1 = 0 ; niv[i].op2 = 1 ;
        niv[i].oper = 0 ;
    }
    indDig = PB095_NBD - 1 ;
    do {
        if(indDig == PB095_NBD - 1) {
            int curNiv = 0 ;
            memset(isValFind,0,PB095_EXPMAX) ;
            do {
              int lgNiv = PB095_NBD - curNiv ;
              if(niv[curNiv].oper >=PB095_NBO){
                  niv[curNiv].oper = 0 ;
                    if(niv[curNiv].op2>= lgNiv-1 ) {
                        if(niv[curNiv].op1>= lgNiv-2 ) {
                            niv[curNiv].op1 = 0 ; niv[curNiv].op2 = 1 ;
                            curNiv-- ;
                            continue ;
                        } else {
                            niv[curNiv].op1++ ;
                            niv[curNiv].op2 = niv[curNiv].op1 + 1 ;
                        }
                    } else {
                        niv[curNiv].op2++ ;
                    }
                }
                Fract_093 newVal = Oper(niv[curNiv].oper,niv[curNiv].vals[niv[curNiv].op1],niv[curNiv].vals[niv[curNiv].op2]);
                niv[curNiv].oper++ ;
                if(newVal.D != 0) { // no division par zero
                    if(lgNiv == 2) {
                        AddVal(isValFind,newVal);
                        continue ;
                    } else {
                        // on avance d'un cran
                        int i,j;
                        for(i=0;i<niv[curNiv].op1;i++) {
                            niv[curNiv+1].vals[i] = niv[curNiv].vals[i] ;
                        }
                        niv[curNiv+1].vals[i++] = newVal ;
                        for(j=niv[curNiv].op1+1;j<lgNiv;j++) {
                            if(j != niv[curNiv].op2 ) {
                                niv[curNiv+1].vals[i++] = niv[curNiv].vals[j] ;
                            }
                        }
                        curNiv++ ;
                    }
                }
            } while(curNiv >= 0) ;
            {
                int i ;
                for(i=0;i<PB095_EXPMAX && isValFind[i] ; i++) ;
                if(i>=maxCons){
                    maxCons = i ;
                    int k ;
                    for(k=0;k<PB095_NBD;k++) maxABCD[k] = niv[0].vals[k].N + '0' ;
                }
            }
        }
        if(niv[0].vals[indDig].N < 10 - PB095_NBD + indDig ){
            niv[0].vals[indDig].N++ ;
            for(;indDig<PB095_NBD-1;){
                indDig++ ;
                niv[0].vals[indDig].N = niv[0].vals[indDig-1].N+1 ;
            }
        } else {
            indDig-- ;
        }
    } while(indDig >=0) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(maxCons < PB095_EXPMAX) {
        sprintf(pbR->strRes,"%s",maxABCD);
        if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d Nb=%d (%s)\n",pbR->pbNum,maxCons,maxABCD) ;
        return 1 ;
    } else {
        if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d ERROR MAX_PRESUMED=%d reached for (%s)\n",pbR->pbNum,maxCons,maxABCD) ;
        return 0 ;
        
    }
}

int PB093(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int a,b,c,d ;
    // 2 parenthesage possible ((a@b)@c)@d) et (a@b)@(c@d)
    // ou @ designe un operateur commutatif (+,x) ou un operateur non commutatif a un seul sens.
    // Pour le premeir choix a,b puis c, et 3 operations. Pour le second choix a,b et " operateurs
    u_int8_t isValFind[PB095_EXPMAX] ;
    //   int tbValues[3000] ;
    int maxCons =0 ;
    int maxABCD = 0 ;
    u_int8_t dig[4] ;
    for(a=1;a<7;a++) {
        dig[0] = a ;
        for(b=a+1;b<8;b++) {
            dig[1] = b ;
            for(c=b+1;c<9;c++) {
                dig[2] = c ;
                for(d=c+1;d<10;d++) {
                    dig[3]= d ;
                    memset(isValFind,0,PB095_EXPMAX) ;
                    do {
                        int op1 ;
                        if(dig[0]>dig[1]) continue ;
                        for(op1=0;op1<PB095_NBO;op1++) {
                            Fract_093 v1 = Oper(op1 ,FractDig(dig[0]),FractDig(dig[1]) );
                            if(v1.D==0) continue ;
                            int op2 ;
                            for(op2=0;op2<PB095_NBO;op2++) {
                                Fract_093 v2 = Oper(op2 ,FractDig(dig[2]),FractDig(dig[3]) );
                                if(v2.D != 0) {
                                    int op3 ;
                                    for(op3=0;op3<PB095_NBO;op3++) {
                                        Fract_093 v3 =  Oper(op3 ,v1,v2);
                                        AddVal(isValFind,v3);
                                    }
                                }
                                v2 = Oper(op2,v1,FractDig(dig[2])) ;
                                if(v2.D != 0) {
                                    int op3 ;
                                    for(op3=0;op3<PB095_NBO;op3++) {
                                        Fract_093 v3 =  Oper(op3 ,v2,FractDig(dig[3]));
                                        AddVal(isValFind,v3);
                                    }
                                }
                                v2 = Oper(op2,v1,FractDig(dig[3])) ;
                                if(v2.D != 0) {
                                    int op3 ;
                                    for(op3=0;op3<PB095_NBO;op3++) {
                                        Fract_093 v3 =  Oper(op3 ,v2,FractDig(dig[2]));
                                        AddVal(isValFind,v3);
                                    }
                                }
                            }
                        }
                    } while(NextPermutRg(dig,4,1)>=0) ;
                    HeapSortUint8(dig,4) ;
                    int i ;
                    for(i=0;i<PB095_EXPMAX && isValFind[i] ; i++) ;
                    if(i>=maxCons){
                        maxCons = i ;
                        maxABCD= a*1000+b*100+c*10+d ;
                    }
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(maxCons < PB095_EXPMAX) {
        sprintf(pbR->strRes,"%d",maxABCD);
        if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d Nb=%d (%d)\n",pbR->pbNum,maxCons,maxABCD) ;
        return 1 ;
    } else {
        if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d ERROR MAX_PRESUMED=%d reached for (%d)\n",pbR->pbNum,maxCons,maxABCD) ;
        return 0 ;
        
    }
}



#define PB100_DEBUG 1
#define PB100_MIN_N 1000000000000
#define INIT_FROM_MIN_N 1

/* RESULT
 delta=-567209652218	; p=707106781186	; n=1000000000000	; d=1279794257938 ;
 delta=0	; p=756872327472	; n=1070379110496	; d=2584123765442 ;
 PB100(72.711029s) blue=756872327473 total=1070379110497 1070379110496x1070379110497=2x756872327472x756872327473
 */
int PB100(PB_RESULT *pbR) {
    //
    // p/n doit approximer 1/sqrt(2).
    // 2 *(p+1)**2 > n**2 > 2* p**2
    // d = 2 *(p+1)**2 - n**2
    // delta = p*(p+1) - (n*(n+1))/2
    // si delta==0 gagne
    
    int64_t n,p,p2,d, delta ;
    pbR->nbClock = clock() ;
    
    n = 1 ; p = 0 ;
    /*
     A decommenter pour partir d'une solution pour trouver la suivante
     delta=0	; p=2	; n=3	; d=9 ;
     delta=0	; p=14	; n=20	; d=50 ;
     delta=0	; p=84	; n=119	; d=289 ;
     delta=0	; p=492	; n=696	; d=1682 ;
     delta=0	; p=2870	; n=4059	; d=9801 ;
     delta=0	; p=16730	; n=23660	; d=57122 ;
     delta=0	; p=97512	; n=137903	; d=332929 ;
     delta=0	; p=568344	; n=803760	; d=1940450 ;
     delta=0	; p=3312554	; n=4684659	; d=11309769 ;
     delta=0	; p=19306982	; n=27304196	; d=65918162 ;
     delta=0	; p=112529340	; n=159140519	; d=384199201 ;
     delta=0	; p=655869060	; n=927538920	; d=2239277042 ;
     delta=0	; p=3822685022	; n=5406093003	; d=13051463049 ;
     delta=0	; p=22280241074	; n=31509019100	; d=76069501250 ;
     delta=0	; p=129858761424	; n=183648021599	; d=443365544449 ;
     delta=0	; p=756872327472	; n=1070379110496	; d=2584123765442 ;
     delta=0	; p=129858761424	; n=183648021599	; d=443365544449 ;
     
     delta=-56711258	;       p= 707106781	; n= 1000000000	;    d=2300791048 ;
     delta=-10168600468 ;    p= 7071067811    ;n= 10000000000   ; d=3804934688 ;
     delta=-71885299958 ;    p= 70710678118  ; n= 100000000000  ; d=97650756322 ;
     delta=-567209652218;    p= 707106781186 ; n= 1000000000000 ; d=127979425793 ;
     */
    // formule pour des petites valeurs d'initialisation
    d = 2*(p+1)*(p+1) - n * n ;
    delta = p * (p + 1) - (n * (n + 1))/2 ;
    
#if defined(INIT_FROM_MIN_N)
    n = PB100_MIN_N ;
    { // on cherche 2*p*p = n*n
        int64_t r = n ;
        int nd ;
        for(nd=0; (r > 10) && ((r % 10) == 0) ; nd++) {  r /= 10 ; }
        r =  r * r;
        p = 0 ;
        while(nd-- > 0 ) {
            int i ;
            p = 10 * p ;
            r = 100 * r ;
            for(i=0; r >= 4 * p + 4 * i + 2 ;i++) {
                r -= 4 * p + 4 * i + 2 ;
            }
            p += i ;
        } // n * n = 2 * p * p + r
        d = 4*p+2-r ; // 4*p+2-r = 4*p+2 + 2*p*p - n*n= 2*(p+1)(p+1) - n*n
        delta = (d-n)/2-p-1 ; // (d-n)/2-p-1 =p*(p+1)-(n*(n+1))/2
    }
#endif
    p2 = 2*p ;
    
#if PB100_DEBUG
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d delta=%lld\t; p=%lld\t; n=%lld\t; d=%lld ;\n",pbR->pbNum,delta ,p2 >> 1,n,d) ;
#endif
    while(1) {
        d -= (n * 2) + 1 ;
        delta -= n + 1 ;
        n++ ;
        if(d < 0)  {
            delta += p2 + 2 ;
            d += (p2 * 2) + 6 ;
            p2 += 2 ;
            if(delta == 0) {
#if PB100_DEBUG
                if(pbR->isVerbose)fprintf(stdout,"\t PB%d delta=%lld\t; p=%lld\t; n=%lld\t; d=%lld ;\n",pbR->pbNum, delta ,p2 >> 1,n,d) ;
#endif
                if(n > PB100_MIN_N ) break ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d blue=%lld total=%lld %lldx%lld=2x%lldx%lld \n"
                              ,pbR->pbNum
                              ,(p2>>1)+1, n+1, n,n+1,(p2>>1),(p2>>1)+1
                              );
    sprintf(pbR->strRes,"%lld",(p2>>1)+1);
    return 1;
    
}

