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
    int i,k ;
    int64_t digits[PB206_NBD] ;
    int64_t pow10[PB206_NBD*2-1] ;
    int nbDtoCheck[PB206_NBD] ;
    int maxDigit[PB206_NBD] ;
    for(i=0;i<PB206_NBD-1;i++) {
        digits[i] = 0 ;
        nbDtoCheck[i] = i/2+1 ;
        maxDigit[i]= 9 ;
    }
    digits[PB206_NBD-1] = 0 ;
    nbDtoCheck[PB206_NBD-1] = PB206_NBD ;
    maxDigit[PB206_NBD-1] = 1 ;
    
    pow10[0]=1 ;
    for(i=0;i<PB206_NBD*2-2;i++) pow10[i+1] = pow10[i]*10 ;
    int is=0;
    N=0 ;
    while(is < PB206_NBD) {
        N2 = N*N ;
        for(k=0;k<nbDtoCheck[is];k++){
            if((N2 % 10)  != 9-k) break ;
            N2 /= 100 ;
        }
        if(k==nbDtoCheck[is]) {
            if(is==PB206_NBD-1) break ;
            is++ ; continue ;
        }
        while(is>=0 && digits[is]==maxDigit[is]) {
            digits[is]=0 ;
            N -= maxDigit[is]*pow10[is] ;
            is-- ;
        }
        if(is< 0) break ;
        else {
            digits[is]++ ;
            N +=pow10[is] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(is<PB206_NBD-1) {
        return 0 ;
    } else {
        N *= 10 ;
    }
    sprintf(pbR->strRes,"%lld",N);
    return 1 ;
}
