//
//  gmp_utils.h
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef gmp_utils_h
#define gmp_utils_h

#include <stdio.h>
#include <gmp.h>
#include "euler_utils.h"
typedef struct FractContG FractContG ;
FractContG * FCG_alloc(void) ;
FractContG * FCG_free(FractContG * FCG);
void FCG_init(FractContG * FCG,uint32_t a0) ;
void FCG_NextCoef(FractContG * FCG,uint32_t a) ;

void FCG_GetNum(const FractContG * FCG,mpz_t op) ;

int FCG_CmpNum0(const FractContG * FCG,const mpz_t op) ;

void FCG_GetNum0(const FractContG * FCG,mpz_t op) ;

int FCG_CmpDen(FractContG * FCG,const mpz_t op) ;
void FCG_GetDen(const FractContG * FCG,mpz_t op) ;

int FCG_CmpDen0(FractContG * FCG,const mpz_t op) ;
void FCG_GetDen0(const FractContG * FCG,mpz_t op) ;



#endif /* gmp_utils_h */
