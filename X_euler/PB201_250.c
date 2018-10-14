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
#define PB206_NBD   8
int PB206(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int64_t N,N2 ;
    int i,k ;
    int digits[PB206_NBD+1] ;
    int64_t pow10[PB206_NBD*2+1] ;
    for(i=0;i<PB206_NBD+1;i++) digits[i] = 0 ;
    pow10[0]=1 ;
    for(i=0;i<PB206_NBD*2;i++) pow10[i+1] = pow10[i]*10 ;
    int is=0;
    N=0 ;
    while(is <= PB206_NBD) {
        N2 = N*N ;
        for(k=0;2*k<=is;k++){
            int dk = N2/pow10[2*k] % 10  ;
            if(dk != 9-k) break ;
        }
        if(is==PB206_NBD) {
            for(;k<=is;k++){
                int dk = N2/pow10[2*k] % 10  ;
                if(dk != 9-k) break ;
            }
            if(k>is) break ;
        } else if(2*k>is) {
 //           printf("%lld ",N);
            is++ ;
            continue ;
        }
        if(is==PB206_NBD) {
            if(digits[is] == 1) {
                digits[is]=0 ;
                N -= pow10[is] ;
                is-- ;
            } else {
                digits[is]++ ;
                N +=pow10[is] ;
                continue ;
            }
        }
        while(is>=0 && digits[is]==9) {
            digits[is]=0 ;
            N -= 9*pow10[is] ;
            is-- ;
        }
        if(is< 0) break ;
        else {
            digits[is]++ ;
            N +=pow10[is] ;
            
        }
    
    }
    
    sprintf(pbR->strRes,"%d",N);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
