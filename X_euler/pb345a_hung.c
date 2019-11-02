//
//  pb345a_hung.c
//  X_euler
//
//  Created by Jeannot on 02/11/2019.
//  Copyright © 2019 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "euler_utils.h"
#include "hungarian.h"
#include "PB_other.h"
#include "p345_data.h"




int PB345a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int32_t S = 0 ;
    const int32_t * M = P345_GetData()  ;
    const int32_t * MH[PB345_DIM] ;
    int i,j ;
    for(i=0;i<PB345_DIM;i++)MH[i] = M+PB345_DIM*i ;
    hungarian_problem_t p;
    
    int matrix_size = hungarian_init(&p,MH,PB345_DIM,PB345_DIM, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;
    hungarian_solve(&p);
    for(i=0;i<PB345_DIM;i++) {
        for(j=0;j<PB345_DIM;j++) {
            if(p.assignment[i][j]) S += M[i*PB345_DIM+j] ;
        }
    }
//        hungarian_print_status(&p);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}
