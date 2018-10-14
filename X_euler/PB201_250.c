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
    int64_t digits[PB206_NBD] ;
    int64_t pow10[PB206_NBD*2-1] ;
    int nbDtoCheck[PB206_NBD] ;
    int64_t maxDigit[PB206_NBD] ;
    pow10[0]=1 ;
    for(i=0;i<PB206_NBD*2-2;i++) pow10[i+1] = pow10[i]*10 ;

    for(i=0;i<PB206_NBD-1;i++) {
        digits[i] = 0 ;
        nbDtoCheck[i] = i/2+1 ;
        maxDigit[i]= 9*pow10[i] ;
    }
    maxDigit[PB206_NBD-2] = 3*pow10[PB206_NBD-2] ;
    digits[PB206_NBD-1] = pow10[PB206_NBD-1] ;
    nbDtoCheck[PB206_NBD-1] = PB206_NBD ;
    maxDigit[PB206_NBD-1] = pow10[PB206_NBD-1] ;
    
    int is=0;
    N=pow10[PB206_NBD-1] ;
    while(is < PB206_NBD) {
        N2 = N*N ;
        for(i=0;i<nbDtoCheck[is];i++){
            if((N2 % 10)  != 9-i) break ;
            N2 /= 100 ;
        }
        if(i==nbDtoCheck[is]) {
            if(is==PB206_NBD-1) break ;
            is++ ; continue ;
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
                    for(i=9;i>=3;i--) {
                        if(mm%10!=i) break;
                        mm/=100;
                    }
                    if(i==2) {
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
