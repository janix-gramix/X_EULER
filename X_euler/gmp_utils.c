//
//  gmp_utils.c
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdlib.h>
#include <gmp.h>
#include "gmp_utils.h"


struct  FractContG {
    mpz_t   n1 ;
    mpz_t   n0 ;
    mpz_t   d1 ;
    mpz_t   d0 ;
}  ;

FractContG * FCG_alloc(void) {
    FractContG * FCG = calloc(1,sizeof(FCG[0])) ;
    mpz_init(FCG->d0) ;
    mpz_init(FCG->d1) ;
    mpz_init(FCG->n0) ;
    mpz_init(FCG->n1) ;
    return  FCG ;
}

FractContG * FCG_free(FractContG * FCG){
    free(FCG);
    return NULL ;
}
void FCG_init(FractContG * FCG,u_int32_t a0) {
    mpz_set_ui (FCG->d0, 0) ;
    mpz_set_ui (FCG->d1, 1) ;
    mpz_set_ui (FCG->n0, 1) ;
    mpz_set_ui (FCG->n1, a0) ;
}

void FCG_NextCoef(FractContG * FCG,u_int32_t a) {
    mpz_addmul_ui(FCG->d0,FCG->d1,a) ;
    mpz_swap (FCG->d0,FCG->d1) ;
    mpz_addmul_ui(FCG->n0,FCG->n1,a) ;
    mpz_swap (FCG->n0,FCG->n1) ;
}

int FCG_CmpNum(const FractContG * FCG,const mpz_t op) {
    return mpz_cmp(FCG->n1,op) ;
}

void FCG_GetNum(const FractContG * FCG,mpz_t op) {
    mpz_set (op, FCG->n1) ;
}

int FCG_CmpNum0(const FractContG * FCG,const mpz_t op) {
    return mpz_cmp(FCG->n0,op) ;
}

void FCG_GetNum0(const FractContG * FCG,mpz_t op) {
    mpz_set (op, FCG->n0) ;
}

int FCG_CmpDen(FractContG * FCG,const mpz_t op) {
    return mpz_cmp(FCG->d1,op) ;
}
void FCG_GetDen(const FractContG * FCG,mpz_t op) {
    mpz_set (op, FCG->d1) ;
}

int FCG_CmpDen0(FractContG * FCG,const mpz_t op) {
    return mpz_cmp(FCG->d0,op) ;
}
void FCG_GetDen0(const FractContG * FCG,mpz_t op) {
    mpz_set (op, FCG->d0) ;
}


