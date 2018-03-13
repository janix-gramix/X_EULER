//
//  euler_utils.h
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef euler_utils_h
#define euler_utils_h

#include <stdio.h>

#define PB_MAX_STRLEN   100

typedef struct PB_RESULT {
    int     pbNum ;
    int     isVerbose ;
    char    strRes[PB_MAX_STRLEN] ;
    clock_t nbClock ;
} PB_RESULT ;

typedef int(*PB_FIND)(PB_RESULT *pbR);


typedef struct PB_CALL {
    int             pbNum ;
    PB_FIND         pbSolve ;
    const   char *  Solution ;
    const   char *  name ;
    
} PB_CALL ;

#define PROTO_PB(NUM)   int PB##NUM(PB_RESULT *pbR)

// PGCD en 32 ou 64 bits
u_int32_t PGCD(u_int32_t n1,u_int32_t n2 ) ;
u_int64_t PGCD64(u_int64_t n1,u_int64_t n2 ) ;
u_int64_t Sqrt64(u_int64_t val) ;

void HeapSortUint8(u_int8_t *H,int n) ;
void HeapSortUint8Rev(u_int8_t *H,int n) ;


typedef  u_int32_t T_prime ;

// calcul de nombre premier recursif (sans stockage)
// routine de completion
// doir retourner 0 pour arreter
typedef int(*TY_CPL_nxtPrime)(void *ctx,T_prime nxtPrime);
u_int32_t FindPrime(T_prime maxValue,void *ctx,TY_CPL_nxtPrime nxtPrime) ;

// generation d'une table de nombre premiers
typedef struct CTX_PRIMETABLE CTX_PRIMETABLE;
// par maxValue du plus grand nombre premier
CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) ;
// par nombre de premiers
CTX_PRIMETABLE * Gen_tablePrimeNb(T_prime maxNb) ;
// recupere la table
const T_prime * GetTbPrime(CTX_PRIMETABLE * ctx) ;
//
u_int32_t GetNbPrime(CTX_PRIMETABLE * ctx) ;
CTX_PRIMETABLE * Free_tablePrime(CTX_PRIMETABLE * ctx) ;
// cherche n dans la table, return true if n is prime
int Search_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) ;
// retourne le rang (2 == rang 0, prmeir nbre premeir
// retourne (-1) si pas premier
int SearchRg_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) ;


u_int32_t FindNbDiv(u_int64_t N, const T_prime *tbPrime) ;
u_int32_t FindNbDivPrime(u_int64_t N, const T_prime *tbPrime) ;
int Is_Prime(u_int64_t N, const T_prime *tbPrime) ;

// return true if P1 and P2 are prime
int Is_Prime2(u_int64_t N1,u_int64_t N2,const T_prime *tbPrime) ;





#endif /* euler_utils_h */
