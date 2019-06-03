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
#include <stdint.h>
#include <time.h>

#define PB_MAX_STRLEN   100

typedef struct PB_RESULT {
    const char    *ident ;
    int     isVerbose ;
    char    strRes[PB_MAX_STRLEN] ;
    clock_t nbClock ;
} PB_RESULT ;

typedef int(*PB_FIND)(PB_RESULT *pbR);


typedef struct PB_CALL {
    const char      *ident  ;
    PB_FIND         pbSolve ;
    const   char *  Solution ;
    const   char *  name ;
    
} PB_CALL ;

#define PROTO_PB(NUM)   int PB##NUM(PB_RESULT *pbR)

// PGCD en 32 ou 64 bits
uint32_t PGCD(uint32_t n1,uint32_t n2 ) ;
uint64_t PGCD64(uint64_t n1,uint64_t n2 ) ;
uint64_t Sqrt64(uint64_t val) ;
uint32_t Sqrt32(uint32_t val) ;
uint16_t Sqrt16(uint16_t val) ;


void HeapSortUint8(uint8_t *H,int n) ;
void HeapSortUint8Rev(uint8_t *H,int n) ;

int NextArrangement(uint8_t *arr,int k, int n) ;
int NextSub(uint8_t *sub,int k, int n) ;
int NextSub16(uint16_t *sub,int k, int n) ;
int NextSub32(uint32_t *sub,int k, int n) ;

int NextPermut(uint8_t *perm,int lg) ;
int NextPermutRg(uint8_t *perm,int lg,int rg) ;
int ChkPermutRg(uint8_t *perm,int lg,int rg) ;
int NextPermutRev(uint8_t *perm,int lg) ;
int NextPermutRgRev(uint8_t *perm,int lg,int rg) ;

typedef struct Decomp {
    uint16_t  Sum ;
    uint16_t * val ;
    uint16_t nbVal ;
} Decomp ;

Decomp  * DecompAlloc(uint16_t Sum) ;
int DecompNext(Decomp  * DeC );
void DecompRewind(Decomp  * DeC );
Decomp * DecompFree(Decomp  * DeC );
/* pour decomposer Sum en elements (ordre alpha decroissant)
Ex ; Sum = 9
 9
 8.1
 7.2
 7.1.1
 6.3
 6.2.1
 6.1.1.1
 5.4
 5.3.1
 5.2.2
 5.2.1.1
 5.1.1.1.1
 4.4.1
 4.3.2
 4.3.1.1
 4.2.2.1
 4.2.1.1.1
 4.1.1.1.1.1
 3.3.3
 3.3.2.1
 3.3.1.1.1
 ...
 */

typedef  uint32_t T_prime ;

// calcul de nombre premier recursif (sans stockage)
// routine de completion
// doir retourner 0 pour arreter
typedef int(*TY_CPL_nxtPrime)(void *ctx,T_prime nxtPrime);
uint32_t FindPrime(T_prime maxValue,void *ctx,TY_CPL_nxtPrime nxtPrime) ;

// generation d'une table de nombre premiers
typedef struct CTX_PRIMETABLE CTX_PRIMETABLE;
// par maxValue du plus grand nombre premier
CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) ;
// par nombre de premiers
CTX_PRIMETABLE * Gen_tablePrimeNb(T_prime maxNb) ;
// recupere la table
const T_prime * GetTbPrime(CTX_PRIMETABLE * ctx) ;
//
uint32_t GetNbPrime(CTX_PRIMETABLE * ctx) ;
CTX_PRIMETABLE * Free_tablePrime(CTX_PRIMETABLE * ctx) ;
// cherche n dans la table, return true if n is prime
int Search_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) ;
// retourne le rang (2 == rang 0, prmeir nbre premeir
// retourne (-1) si pas premier
int SearchRg_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) ;

int MRP_isPrime(uint64_t n) ;


uint32_t FindNbDiv(uint64_t N, const T_prime *tbPrime) ;
uint32_t FindNbDivPrime(uint64_t N, const T_prime *tbPrime) ;
int Is_Prime(uint64_t N, const T_prime *tbPrime) ;

int Is_Prime32(uint32_t N, const T_prime *tbPrime) ;

// return true if P1 and P2 are prime
int Is_Prime2(uint64_t N1,uint64_t N2,const T_prime *tbPrime) ;


typedef struct FractCont64 {
    uint64_t N0 ;
    uint64_t D0 ;
    uint64_t N1 ;
    uint64_t D1 ;
    
} FractCont64 ;

void NextFract(FractCont64 * F, int a) ; // compute next continued fraction


#endif /* euler_utils_h */
