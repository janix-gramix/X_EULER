//
//  PB201_250.c
//  X_euler
//
//  Created by Jeannot on 14/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "PB201_250.h"
#define PB206_NBD   9



int PB206(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int64_t N,N2 ;
    int i ;
    int is=0;
    int64_t digits[PB206_NBD] ;
    int64_t pow10[2*PB206_NBD] ;
    int64_t maxDigit[PB206_NBD] ;
    pow10[0]=1 ;
    for(i=0;i<2*PB206_NBD-2;i++) pow10[i+1] = pow10[i]*10 ;


    for(i=0;i<PB206_NBD-1;i++) {
        digits[i] = 0 ;
         maxDigit[i]= 9*pow10[i] ;
    }
    
    maxDigit[PB206_NBD-2] = 3*pow10[PB206_NBD-2] ;
    digits[PB206_NBD-1] = pow10[PB206_NBD-1] ;
     maxDigit[PB206_NBD-1] = pow10[PB206_NBD-1] ;
     is=0;
    N=pow10[PB206_NBD-1] ;
    while(1){
        if(is & 1) {
            is++ ; continue ;
        } else {
             if(is<PB206_NBD-1) { // verification of added digit
                N2 = N*N ;
                 if( ( (N2/pow10[is] ) % 10) == 9-(is/2) ) {
                     is++ ; continue ;
                }
            } else { // verification for high digits (low are OK by construction)
                N2 = N*N ;
                for(i=PB206_NBD-2;i>=PB206_NBD/2;i--){
                    if(((N2/pow10[2*i]) % 10) != 9-i) break ;
                }
                if(i<PB206_NBD/2) {
                    break ;
                }
                is-- ;
            }
        }

        while(is>=0 && digits[is]==maxDigit[is]) {
            N -= digits[is] ;
            digits[is--]=0 ;
        }
        if(is< 0) break ;
        else { digits[is] += pow10[is] ; N +=pow10[is] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(is<PB206_NBD-1) {
        return 0 ;
    } else {
        N *= 10 ;
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s %lld **2 = %lld\n",pbR->ident,N,N*N) ;
        sprintf(pbR->strRes,"%lld",N);
        return 1 ;
    }
}



#define U 10203040506070809LL
#define A  1000000000000000LL
#define B    10000000000000LL
#define C      100000000000LL
#define V        9090909090LL


int PB206a(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int a,b,c,i;
    int64_t u,m,mm;
    for(a=0;a<10;a++) {
        for(b=0;b<10;b++) {
            for(c=0;c<10;c++) {
                u=U+a*A+b*B+c*C;
                m=Sqrt64(u); m |= 1;
                for(;(mm=m*m)<u+V;m+=2){
                    for(i=9;i>=1;i--) {
                        if(mm%10!=i) break;
                        mm/=100;
                    }
                    if(i==0) {
                        m *= 10 ;
                        pbR->nbClock = clock() - pbR->nbClock ;
                        if(pbR->isVerbose) fprintf(stdout,"\tPB%s %lld **2 = %lld\n",pbR->ident,m,m*m) ;
                        sprintf(pbR->strRes,"%lld",m);
                        return 1 ;
                        
                    }
                }
            }
        }
    }
    return 0 ;
}

static int Chk357(int nbFact,int *factP,int n,const T_prime * tbPrime) {
    if(nbFact==2) {
        return Is_Prime32(factP[1]+2,tbPrime) && Is_Prime32(n+1, tbPrime) ;
    }
    int nbDiv2 = 1<< (nbFact -1) ;
    int mask ;
    for(mask=0;mask<nbDiv2;mask++) {
        int d1 = 1 ;
        int d2 = 2 ;
        int i ;
        for(i=0;i<nbFact-1;i++ ) {
            if((1<<i) & mask) d1 *= factP[i+1] ;
            else d2 *= factP[i+1] ;
        }
     if(!Is_Prime32(d1+d2,tbPrime)) return 0 ;
    }
    return 1 ;
}

#define PB357_MAXP  100000000
// #define PB357_MAXP  1000
int PB357(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;

    pbR->nbClock = clock();
    if((ctxP = Gen_tablePrime(PB357_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }

    int64_t Sum = 1 ;
    int32_t nb = 1 ;
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int factP[40] ;
    int nbFact = 0 ;
    factP[0] = 2 ;

    int i,j;
    for(i=1;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        int n  = p-1 ;
        if((p & 3)==1) continue ;
        if(!Is_Prime32((p+3)/2, tbPrime)) continue ;
       int n1 = n ;
        int nbFact = 0 ;
        for(j=0;n>1 && j<i;j++){
            int pj = tbPrime[j] ;
             if((n % pj) == 0) {
                n /= pj ;
                if((n % pj)== 0 ) break ;
                factP[nbFact++] = pj ;
            }
        }
        if(n > 1) continue ;
        if(Chk357(nbFact, factP,n1,tbPrime)) {
            nb++; Sum += n1 ;
        }

        
    }
    Free_tablePrime(ctxP) ;

    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;

    sprintf(pbR->strRes,"%lld",Sum);
    return 1 ;
    
}
int PB357a(PB_RESULT *pbR) {
        CTX_PRIMETABLE * ctxP  ;
        
        pbR->nbClock = clock();
        if((ctxP = Gen_tablePrime(PB357_MAXP)) == NULL) {
            fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
            return 0 ;
        }
        
        const T_prime * tbPrime = GetTbPrime(ctxP);
        int nbPrime = GetNbPrime(ctxP) ;
        u_int8_t *isPG = calloc(PB357_MAXP+1,sizeof(isPG[0]));
        int i,d,k;
        for(i=0;i<nbPrime;i++) {
            isPG[tbPrime[i]-1] = 1 ;
        }
    int ip = 0 ;
    int p = tbPrime[ip] ;
    int imax = Sqrt32(PB357_MAXP);
        for(i=2;i<=imax;i++) {
            while(p<i*2) { p=tbPrime[++ip] ; }
            int kMax = PB357_MAXP/i ;
            int ip2 = ip ;
            int p2 = p ;
            k = i ;
            while(p2<=kMax+i) {
                for(;k<p2-i;k++) {
                    if(isPG[k*i]) isPG[k*i] = 0 ;
                 }
                p2 = tbPrime[++ip2] ; k++ ;
           }
            for(;k<=kMax;k++) {
                isPG[k*i] = 0 ;
            }
         }
        int64_t Sum = 0 ;
        int32_t nb = 0 ;
    for(i=1;i<PB357_MAXP;i++) {
        if(isPG[i]) { nb++ ; Sum += i ; }
    }
    Free_tablePrime(ctxP) ;
            
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
        
        sprintf(pbR->strRes,"%lld",Sum);
        return 1 ;

}
