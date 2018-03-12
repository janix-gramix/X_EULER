//
//  main.c
//  EulerProject
//  Created by Jeannot on 26/01/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

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


int PB001(PB_RESULT *pbR) {
    int32_t i ;
    int64_t S = 0 ;
    int32_t i3 = 0;
    int32_t i5 = 0 ;
    pbR->nbClock = clock()  ;
    for(i=0;i<1000;i++) {
        if(i3 == 0  || i5 ==0) S += i ;
        if(++i3==3) i3 =0 ;
        if(++i5==5) i5 =0 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",S);
    return 1 ;
}

int PB002(PB_RESULT *pbR) {
    int64_t S = 0 ;
    int32_t u0 = 1 ;
    int32_t u1 = 2 ;
    pbR->nbClock = clock()  ;

    while ( u1 < 4000000) {
        if((u1 & 1) == 0) S += u1 ;
        u1 += u0 ;
        u0 = u1 -u0 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",S);
    return 1 ;
}

int PB003(PB_RESULT *pbR) {
    u_int64_t n = 600851475143 ;
    u_int64_t d = 2 ;
    u_int64_t d2 = 2*2 ;
    u_int64_t bp = 1 ;
    pbR->nbClock = clock() ;
    while ( d2 < n) {
        while ( (n % d) == 0) {
            bp = d ;
            n /= d ;
        }
        d++ ;
        d2 += d * 2 - 1 ;
    }
    if (n > bp) bp = n ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",bp);
    return 1 ;
    
}

int ChkPalyn(int n) {
    // on suppose que n > 100000
    if (   (((n % 100001) % 10010) % 1100) ==0  ){
        return 1 ;
    }
    return 0 ;
}

int PB004(PB_RESULT *pbR) {
    // on va assurer n1 >= n2
    // on balaye a somme constante (en partant du max S/2 * S/2)
    // on decremente la somme des que le produit est <  bestP
    // on termine quand le max S/2 * S/2 est < bestP
    int32_t n1 = 999 ;
    int32_t n2 = 999 ;
    int32_t P = n1*n2 ;
    int32_t S = n1+n2 ;
    int32_t bestP = 100000 ;
    int32_t bestn1 = 0 ;
    pbR->nbClock = clock() ;
    while ( 1 ) {
        //       printf("%d=%dx%d ",P,n1,n2);
        if ( (P > bestP) && ChkPalyn(P) ) {
            bestP = P ;
            bestn1 =n1 ;
        }
        if((P > bestP) && (n1 < 999) ) {
            P += n2 -n1 -1 ;
            n1++ ;
            n2-- ;
        } else {
            S-- ;
            n1 =(S+1)/2 ;
            n2 = S-n1 ;
            P = n1*n2 ;
            if(P <= bestP ) {
                break ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%d P=%d = %d x %d \n",pbR->pbNum,bestP,bestn1,bestP/ bestn1);
    if ( bestn1 ) {
        sprintf(pbR->strRes,"%d",bestP);
        return 1 ;
    } else {
        return 0 ;
    }
}

u_int32_t PGCD(u_int32_t n1,u_int32_t n2 ) {
    u_int32_t r ;
    if (n1 > n2) {
        r = n2 ;
        n2 = n1 ;
    } else {
        r = n1 ;
    }
    while ( r > 0) {
        n1 = n2 ;
        n2 = r ;
        r = n1 % n2 ;
    }
    return n2 ;
}

u_int64_t PGCD64(u_int64_t n1,u_int64_t n2 ) {
    u_int64_t r ;
    if (n1 > n2) {
        r = n2 ;
        n2 = n1 ;
    } else {
        r = n1 ;
    }
    while ( r > 0) {
        n1 = n2 ;
        n2 = r ;
        r = n1 % n2 ;
    }
    return n2 ;
}


int PB005(PB_RESULT *pbR) {
    int PPCM = 1;
    int i ;
    for(i=2;i<20;i++) {
        PPCM *= i/ PGCD(PPCM,i);
    }
    sprintf(pbR->strRes,"%d",PPCM);
    return 1 ;
}

int PB006(PB_RESULT *pbR) {
    int64_t n = 100 ;
    int64_t Sn = (n *(n+1)) / 2 ;
    int64_t Sn2 = (n*(n+1)*(2*n+1)) / 6 ;
    int64_t diff = Sn * Sn - Sn2 ;
    sprintf(pbR->strRes,"%lld",diff);
    return 1;
}

// #define PB007_NB_PRIME_ASK 100000000
// #define PB007_NB_MAX 2000000000

// #define PB007_NB_PRIME_ASK 100
// #define PB007_NB_MAX 542

// #define PB007_NB_PRIME_ASK 400001
#define PB007_NB_PRIME_ASK 10001
#define PB007_NB_MAX 2000000

#define PB007_GENTABLE  0
#define PB007_PRINTABLE 0

#define T_prime u_int32_t

// routine de completion
// doir retourner 0 pour arreter
typedef int(*TY_CPL_nxtPrime)(void *ctx,T_prime nxtPrime);


typedef struct CTX_PB007 {
    u_int32_t nbPrime ;
    T_prime lastPrime ;
    u_int32_t nbAsk ;
#if PB007_GENTABLE
    T_prime *tbPrime ;
#endif
    
} CTX_PB007 ;

int PB0007_nxtPrime(void *ptCtx,T_prime nxtPrime) {
    CTX_PB007 * ctx = (CTX_PB007 *) ptCtx ;
    ctx->lastPrime = nxtPrime ;
#if PB007_GENTABLE
    ctx->tbPrime[ctx->nbPrime] = nxtPrime ;
#endif
    ctx->nbPrime++ ;
    if(ctx->nbPrime >= ctx->nbAsk) {
        return 0 ;
    } else {
        return 1;
    }
}


#define IS_Composed(p)  (isComposed[(p)/8] & (1 << ((p) & 0x7)) )
#define SET_Composed(p)  (isComposed[(p)/8] |=  (1 << ((p)  & 0x7)) )


u_int32_t PB007a_FindPrime(T_prime nbMax,void *ctx,TY_CPL_nxtPrime nxtPrime) {
    int32_t nSqrt = 1+ (int32_t)sqrt( (double) nbMax ) ;
    T_prime sizeTable = nbMax ;
    T_prime sizeTable2 = nbMax >> 1 ;
    u_int8_t *isComposed = calloc( (sizeTable+15) /16,  sizeof(isComposed[0])) ;
    u_int32_t nbPrime = 0 ;
    T_prime curPrime = 0 ;
    T_prime lastPrime = 0 ;

    // on traite le cas 2 a part, car apres on ne traite que les nombres impairs
    nbPrime ++ ;
    lastPrime = 2 ;
    if(nxtPrime(ctx,2) == 0) {
        return nbPrime ;
    }
    while(++curPrime < sizeTable2) {
        if( ! IS_Composed(curPrime)) {
            lastPrime = 2*curPrime+1 ;
            nbPrime++ ;
            if(nxtPrime(ctx,lastPrime) == 0) {
                break ; // demande d'arret
            }
            if(lastPrime < nSqrt ) {
                T_prime icp2 ;
                for(icp2 = curPrime+lastPrime ; icp2 < sizeTable2; icp2 += lastPrime ) {
                    SET_Composed(icp2) ;
                }
            }
        }
    }
    free(isComposed);
    return nbPrime ;
}


//
// On replie la table.
// la taille de la table peut être quelconque
// Plus rapide pour taille nSqrt ou 32368 si trop grande valeurs
//
u_int32_t PB007b_FindPrime(T_prime nbMax,void *ctx,TY_CPL_nxtPrime nxtPrime) {
    int32_t nSqrt = 1+ (int32_t)sqrt( (double) nbMax ) ;
    int isEnd = 0;
    T_prime *tbPrime = malloc(nSqrt * sizeof(tbPrime[0])) ;
    int32_t *offSet = malloc(nSqrt * sizeof(offSet[0])) ;
    int32_t sizeTable =  (nSqrt < 32768) ? nSqrt : 32768 ;
    int8_t *isComposed = calloc( sizeTable , sizeof(isComposed[0])) ;
    T_prime nbPrime = 0 ;
    T_prime lastPrime = 0 ;
    T_prime offSetTable = 0 ;
    T_prime nbPrimeSqrt = 0 ;
    
    nbPrime ++ ;
    lastPrime = 2 ;
    if(nxtPrime(ctx,2) == 0) {
        return nbPrime ;
    }
    // remarque le nb premier 2 n'est pas stocke dans tbPrime car pas utilise dans le crible erasto.
    // pour commencer a 3 donc indPrime = (3>>1)
    T_prime indPrime = 1 ;
    
    while ( 1) {
        T_prime icp ;
        for(icp=indPrime - offSetTable ;icp < sizeTable; icp++ ) {
            if(!isComposed[icp] ) {
                lastPrime = 2*(icp+offSetTable)+1 ;

                if(lastPrime < nSqrt) {
                    T_prime icp2 ;
                    tbPrime[nbPrimeSqrt] = lastPrime ;
                    for(icp2 = icp + lastPrime; icp2 < sizeTable ; icp2 += lastPrime ) {
                        isComposed[icp2] = 1 ;
                    }
                    offSet[nbPrimeSqrt++] = icp2 - sizeTable ;
                }
                nbPrime++ ;
                if(nxtPrime(ctx,lastPrime) == 0) {
                    isEnd = 1;
                    break ; // demande d'arret
                }
            }
        }
        offSetTable += sizeTable ;
        if(isEnd || offSetTable >= nbMax) break ;
        indPrime = offSetTable ;
        if ( offSetTable + sizeTable > nbMax) { sizeTable = nbMax - offSetTable ; }
        memset(isComposed,0,sizeTable) ;
        {
            int np ;
            for(np=0;np<nbPrimeSqrt;np++) {
                int32_t p = tbPrime[np] ;
                int32_t indPrime = offSet[np] ;
                while ( indPrime < sizeTable) {
                    isComposed[indPrime] = 1 ;
                    indPrime += p ;
                }
                offSet[np] = indPrime - sizeTable ;
            }
        }
    }
    free(tbPrime);
    free(offSet) ;
    free(isComposed);
    return nbPrime ;
}



int PB007(PB_RESULT *pbR) {
    CTX_PB007 ctx ;
    T_prime nbMax = PB007_NB_MAX ;
    ctx.nbPrime = 0 ;
    ctx.lastPrime = 0 ;
    ctx.nbAsk = PB007_NB_PRIME_ASK ;
    pbR->nbClock = clock() ;
#if PB007_GENTABLE
    ctx.tbPrime = malloc(ctx.nbAsk * sizeof(ctx.tbPrime[0])) ;
#endif
    PB007a_FindPrime(nbMax,&ctx,PB0007_nxtPrime) ;
#if  PB007_PRINTABLE
     {
        int i;
        for(i=0;i<ctx.nbPrime ;i++) {
            printf("%d ",ctx.tbPrime[i]) ;
        }
     }
     printf("\n");
#endif
    if(pbR->isVerbose) {
        pbR->nbClock = clock() -  pbR->nbClock ;
        printf("\t PB%da(%.6fs)  prime n°%d  = %d\n",pbR->pbNum,(float)pbR->nbClock / CLOCKS_PER_SEC ,ctx.nbPrime,ctx.lastPrime) ;
    }
    ctx.nbPrime = 0 ;
    ctx.lastPrime = 0 ;
    pbR->nbClock = clock();
    PB007b_FindPrime(nbMax,&ctx,PB0007_nxtPrime) ;
#if  PB007_PRINTABLE
    {
        int i;
        for(i=0;i<ctx.nbPrime ;i++) {
            printf("%d ",ctx.tbPrime[i]) ;
        }
    }
    printf("\n");
#endif
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",ctx.lastPrime);
    if(pbR->isVerbose) fprintf(stdout,"\t PB%db  prime n°%d  = %d\n",pbR->pbNum,ctx.nbPrime,ctx.lastPrime) ;
    return 1 ;
    
    
}

#define PB008_NB    13
int PB008(PB_RESULT *pbR) {
    char strDig[] =
    "73167176531330624919225119674426574742355349194934"
    "96983520312774506326239578318016984801869478851843"
    "85861560789112949495459501737958331952853208805511"
    "12540698747158523863050715693290963295227443043557"
    "66896648950445244523161731856403098711121722383113"
    "62229893423380308135336276614282806444486645238749"
    "30358907296290491560440772390713810515859307960866"
    "70172427121883998797908792274921901699720888093776"
    "65727333001053367881220235421809751254540594752243"
    "52584907711670556013604839586446706324415722155397"
    "53697817977846174064955149290862569321978468622482"
    "83972241375657056057490261407972968652414535100474"
    "82166370484403199890008895243450658541227588666881"
    "16427171479924442928230863465674813919123162824586"
    "17866458359124566529476545682848912883142607690042"
    "24219022671055626321111109370544217506941658960408"
    "07198403850962455444362981230987879927244284909188"
    "84580156166097919133875499200524063689912560717606"
    "05886116467109405077541002256983155200055935729725"
    "71636269561882670428252483600823257530420752963450" ;
    int lg = (int) strlen(strDig) ;
    u_int64_t bestProd = 0 ;
    int32_t bestIndex = 0  ;
    int32_t indNxt ;
    pbR->nbClock = clock() ;

    for(indNxt=PB008_NB;indNxt<lg;indNxt++) {
        int i ;
        u_int64_t prod = 1 ;
        for(i=0;i<PB008_NB;i++) {
            prod *= ( strDig[indNxt -PB008_NB+i] - '0' ) ;
            if(prod == 0) {
                indNxt += i ;
                break ;
            }
        }
        if(prod > bestProd) {
            bestProd = prod ;
            bestIndex = indNxt - PB008_NB ;
        }
    }
    strDig[bestIndex+PB008_NB] = 0 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%d Prod(%s)  =%lld\n",pbR->pbNum,strDig+bestIndex,bestProd) ;
    sprintf(pbR->strRes,"%lld",bestProd) ;
    return 1 ;
}

#define PB009_SUM   1000
int PB009(PB_RESULT *pbR) {
// pour m > n > 0
// a = m**2 - n**2 ; b = 2*m*n ; c = m**2 + n**2 ;
// donc S = 2 * (m**2 + m*n) = 2 * m * (m+n)
// 1000 = 2 * 25 * 20
// m = 20 n=5
    int m = 20 ;
    int n = 5 ;
    int a = m*m - n*n ;
    int b = 2*m*n ;
    int c = m*m + n*n ;
   if(pbR->isVerbose)fprintf(stdout,"\t PB%d %dx%dx%d=%d %d+%d+%d=%d %d**2+%d**2=%d %d**2=%d\n",pbR->pbNum, a,b,c,a*b*c,a,b,c,a+b+c,a,b,a*a+b*b,c,c*c) ;
   sprintf(pbR->strRes,"%d",a*b*c) ;
    return 1;
}

typedef struct CTX_PB010 {
    u_int64_t sum ;
} CTX_PB010 ;

#define PB010_NB_MAX    2000000
int PB010_nxtPrime(void *ptCtx,T_prime nxtPrime) {
    CTX_PB010 * ctx = (CTX_PB010 *) ptCtx ;
    if(nxtPrime > PB010_NB_MAX) return 0 ;
    ctx->sum += nxtPrime ;
    return 1;
}

int PB010(PB_RESULT *pbR) {    CTX_PB010 ctx ;
    T_prime nbMax = PB010_NB_MAX ;
    ctx.sum = 0 ;
    pbR->nbClock = clock() ;
    PB007b_FindPrime(nbMax,&ctx,PB010_nxtPrime) ;
    pbR->nbClock = clock() - pbR->nbClock ;

    if(pbR->isVerbose)fprintf(stdout,"\t PB%d SUM(prime<%d)  = %lld\n",pbR->pbNum,nbMax,ctx.sum) ;
    sprintf(pbR->strRes,"%lld",ctx.sum) ;
    return 1 ;
}

#define PB11_SIZE   20
int PB011(PB_RESULT *pbR) {
    static int tbIn[PB11_SIZE*PB11_SIZE] = {
         8,02,22,97,38,15,00,40,00,75,04,05,07,78,52,12,50,77,91, 8,
        49,49,99,40,17,81,18,57,60,87,17,40,98,43,69,48,04,56,62,00,
        81,49,31,73,55,79,14,29,93,71,40,67,53,88,30,03,49,13,36,65,
        52,70,95,23,04,60,11,42,69,24,68,56,01,32,56,71,37,02,36,91,
        22,31,16,71,51,67,63,89,41,92,36,54,22,40,40,28,66,33,13,80,
        24,47,32,60,99,03,45,02,44,75,33,53,78,36,84,20,35,17,12,50,
        32,98,81,28,64,23,67,10,26,38,40,67,59,54,70,66,18,38,64,70,
        67,26,20,68,02,62,12,20,95,63,94,39,63, 8,40,91,66,49,94,21,
        24,55,58,05,66,73,99,26,97,17,78,78,96,83,14,88,34,89,63,72,
        21,36,23, 9,75,00,76,44,20,45,35,14,00,61,33,97,34,31,33,95,
        78,17,53,28,22,75,31,67,15,94, 3,80,04,62,16,14, 9,53,56,92,
        16,39,05,42,96,35,31,47,55,58,88,24,00,17,54,24,36,29,85,57,
        86,56,00,48,35,71,89,07,05,44,44,37,44,60,21,58,51,54,17,58,
        19,80,81,68,05,94,47,69,28,73,92,13,86,52,17,77,04,89,55,40,
        04,52, 8,83,97,35,99,16,07,97,57,32,16,26,26,79,33,27,98,66,
        88,36,68,87,57,62,20,72,03,46,33,67,46,55,12,32,63,93,53,69,
        04,42,16,73,38,25,39,11,24,94,72,18, 8,46,29,32,40,62,76,36,
        20,69,36,41,72,30,23,88,34,62,99,69,82,67,59,85,74,04,36,16,
        20,73,35,29,78,31,90,01,74,31,49,71,48,86,81,16,23,57,05,54,
        01,70,54,71,83,51,54,69,16,92,33,48,61,43,52,01,89,19,67,48 } ;
    int ic,ir ;
    int bestProd = 0 ;
    int bestInd = 0 ;
    int bestDelta = 0;
    int prod ;
    pbR->nbClock = clock()  ;
    for(ir=0;ir<PB11_SIZE;ir++) {
        for(ic=0;ic<PB11_SIZE;ic++) {
            int ind = ir * PB11_SIZE + ic ;
            if(ic <= PB11_SIZE - 4) {
                prod = tbIn[ind] * tbIn[ind+1] * tbIn[ind+2] * tbIn[ind+3] ; // RIGHT
                if(prod > bestProd) {
                    bestProd = prod ;
                    bestInd = ind ;
                    bestDelta = 1;
                }
                if(ir <=  PB11_SIZE - 4) {
                    prod = tbIn[ind] * tbIn[ind+1+PB11_SIZE] * tbIn[ind+2+2*PB11_SIZE] * tbIn[ind+3+3*PB11_SIZE] ; // DIAG down-right
                    if(prod > bestProd) {
                        bestProd = prod ;
                        bestInd = ind ;
                        bestDelta = 1+PB11_SIZE;
                    }
                }
            }
            if(ir <=  PB11_SIZE - 4) {
                prod = tbIn[ind] * tbIn[ind+PB11_SIZE] * tbIn[ind+2*PB11_SIZE] * tbIn[ind+3*PB11_SIZE] ; // down
                if(prod > bestProd) {
                    bestProd = prod ;
                    bestInd = ind ;
                    bestDelta = PB11_SIZE;
                }
                if(ir >=  3) {
                    prod = tbIn[ind] * tbIn[ind-1+PB11_SIZE] * tbIn[ind-2+2*PB11_SIZE] * tbIn[ind-3+3*PB11_SIZE] ; // DIAG down-left
                    if(prod > bestProd) {
                        bestProd = prod ;
                        bestInd = ind ;
                        bestDelta = PB11_SIZE-1;
                    }
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d  [%d,%d]x...x[%d,%d]=%dx%dx%dx%d  = %d\n"
           ,pbR->pbNum
           , bestInd / PB11_SIZE ,bestInd % PB11_SIZE
           , (bestInd+3*bestDelta) / PB11_SIZE ,(bestInd+3*bestDelta) % PB11_SIZE
           , tbIn[bestInd],tbIn[bestInd+bestDelta],tbIn[bestInd+2*bestDelta],tbIn[bestInd+3*bestDelta]
           ,bestProd ) ;
    sprintf(pbR->strRes,"%d",bestProd);
    return 1 ;
}

typedef struct CTX_PRIMETABLE {
    u_int32_t   nbPrime ;
    T_prime     maxValue ;
    u_int32_t   maxNbPrime ;
    T_prime     *tbPrime ;
} CTX_PRIMETABLE ;

int CPL_tablePrime(void *ptCtx,T_prime nxtPrime) {
    CTX_PRIMETABLE * ctx = (CTX_PRIMETABLE *) ptCtx ;
    if(nxtPrime > ctx->maxValue) return 0 ;
    ctx->tbPrime[ctx->nbPrime] = nxtPrime ;
    ctx->nbPrime++ ;
    if(ctx->nbPrime >= ctx->maxNbPrime) {
       ctx->maxValue = nxtPrime * nxtPrime ;
       return 0 ;
    }
    return 1;
}

CTX_PRIMETABLE * Free_tablePrime(CTX_PRIMETABLE * ctx) {
    if(ctx != NULL) free(ctx->tbPrime) ;
    free(ctx);
    return NULL ;
}
#define USE_PB007b  1

CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) {
    CTX_PRIMETABLE *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    ctx->maxValue = maxValue ;
    if(maxValue > 100) {
        ctx->maxNbPrime = (T_prime) (1+ maxValue / (log((double)maxValue) - 4)) ;
    } else {
        ctx->maxNbPrime = 30 ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime(ctx) ; }
#if USE_PB007b
    PB007b_FindPrime(maxValue,ctx,CPL_tablePrime) ;
#else
    PB007a_FindPrime(maxValue,ctx,CPL_tablePrime) ;
#endif
    return ctx ;
}

CTX_PRIMETABLE * Gen_tablePrimeNb(T_prime maxNb) {
    CTX_PRIMETABLE *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    if(maxNb < 30)  {
        ctx->maxNbPrime = 30 ;
    } else {
        ctx->maxNbPrime = maxNb ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime(ctx) ; }
#if USE_PB007b
    ctx->maxValue = 0x7fffffff ;
    PB007b_FindPrime(ctx->maxValue,ctx,CPL_tablePrime) ;
#else
    ctx->maxValue = 0xfffffff ;
    PB007a_FindPrime(ctx->maxValue,ctx,CPL_tablePrime) ;
#endif

    return ctx ;
}


u_int64_t Sqrt64(u_int64_t val) {
    // on utilse sqrt beaucoup plus efficace (30 fois)
#if 0
    return (u_int64_t) sqrt((double)val);
#else
    // Sylvain racine carre a la main classique
    uint64_t sq = 0, b = 1LL << (sizeof(val) * 8 - 2);
    while (b > val) {
        b >>= 2;
    }
    while (b) {
        uint64_t c = sq | b;
        sq >>= 1 ;
        if (val >= c) {
            val -= c;
            sq |= b;
        }
        b >>= 2;
    }
    return sq ;
#endif
}

u_int32_t FindNbDiv(u_int64_t N, T_prime *tbPrime) {
    u_int64_t sqr = Sqrt64(N) ;
    T_prime p ;
    u_int32_t nDiv = 1;
    for(p= *tbPrime++ ; p<= sqr ;p=*tbPrime++) {
        u_int32_t exp = 0 ;
        while((N % p) == 0) {
            N /= p ;
            exp++ ;
        }
        nDiv *= (exp+1) ;
    }
    if(N > 1) nDiv *= 2 ;
    return nDiv ;
}

u_int32_t FindNbDivPrime(u_int64_t N, T_prime *tbPrime) {
    u_int64_t sqr = Sqrt64(N) ;
    T_prime p ;
    u_int32_t nDivP = 0;
    for(p= *tbPrime++ ; p<= sqr ;p=*tbPrime++) {
        u_int32_t exp = 0 ;
        while((N % p) == 0) {
            N /= p ;
            exp = 1 ;
        }
        nDivP += exp ;
    }
    if(N > 1) nDivP++ ;
    return nDivP ;
}


static int Is_Prime(u_int64_t N, T_prime *tbPrime) {
    u_int64_t sqr = Sqrt64(N) ;
    T_prime p ;
    for(p= *tbPrime++ ; p<= sqr ;p=*tbPrime++) {
        if((N % p) ==0) return 0 ;
    }
    return 1 ;
}


// return true if P1 and P2 are prime
static int Is_Prime2(u_int64_t N1,u_int64_t N2, T_prime *tbPrime) {
    u_int64_t sqr = Sqrt64((N1>N2) ? N1 : N2 ) ;
    
    T_prime p ;
    for(p= *tbPrime++ ; p<= sqr ;p=*tbPrime++) {
        if((N1 % p) ==0 || (N2 % p) ==0) return 0 ;
    }
    return 1 ;
}


#define PB012_MAX_SQPRIME   1000000 // permet de factoriser jusqu'a 10**12
int PB012(PB_RESULT *pbR) {
   CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB012_MAX_SQPRIME)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    } else {
        int64_t n = 0 ;
        int64_t ntr;
        int64_t nbDiv ;
        do {
            n++ ;
            ntr =(n*(n+1))/2 ;
            nbDiv=FindNbDiv( ntr,ctxP->tbPrime) ;
            
        } while (nbDiv <= 500) ;
        pbR->nbClock = clock() - pbR->nbClock ;
        if(pbR->isVerbose)fprintf(stdout,"\t PB%d triangle(%lld)=%lld a %lld diviseurs\n",pbR->pbNum,n,ntr,nbDiv);
        sprintf(pbR->strRes,"%lld",ntr);
    }
    Free_tablePrime(ctxP);
    return 1;
}

#define PB013_NBNUM     100
#define PB013_LENNUM    50
int PB013(PB_RESULT *pbR) {
    static char strNumber[PB013_NBNUM*PB013_LENNUM] =
"37107287533902102798797998220837590246510135740250"
"46376937677490009712648124896970078050417018260538"
"74324986199524741059474233309513058123726617309629"
"91942213363574161572522430563301811072406154908250"
"23067588207539346171171980310421047513778063246676"
"89261670696623633820136378418383684178734361726757"
"28112879812849979408065481931592621691275889832738"
"44274228917432520321923589422876796487670272189318"
"47451445736001306439091167216856844588711603153276"
"70386486105843025439939619828917593665686757934951"
"62176457141856560629502157223196586755079324193331"
"64906352462741904929101432445813822663347944758178"
"92575867718337217661963751590579239728245598838407"
"58203565325359399008402633568948830189458628227828"
"80181199384826282014278194139940567587151170094390"
"35398664372827112653829987240784473053190104293586"
"86515506006295864861532075273371959191420517255829"
"71693888707715466499115593487603532921714970056938"
"54370070576826684624621495650076471787294438377604"
"53282654108756828443191190634694037855217779295145"
"36123272525000296071075082563815656710885258350721"
"45876576172410976447339110607218265236877223636045"
"17423706905851860660448207621209813287860733969412"
"81142660418086830619328460811191061556940512689692"
"51934325451728388641918047049293215058642563049483"
"62467221648435076201727918039944693004732956340691"
"15732444386908125794514089057706229429197107928209"
"55037687525678773091862540744969844508330393682126"
"18336384825330154686196124348767681297534375946515"
"80386287592878490201521685554828717201219257766954"
"78182833757993103614740356856449095527097864797581"
"16726320100436897842553539920931837441497806860984"
"48403098129077791799088218795327364475675590848030"
"87086987551392711854517078544161852424320693150332"
"59959406895756536782107074926966537676326235447210"
"69793950679652694742597709739166693763042633987085"
"41052684708299085211399427365734116182760315001271"
"65378607361501080857009149939512557028198746004375"
"35829035317434717326932123578154982629742552737307"
"94953759765105305946966067683156574377167401875275"
"88902802571733229619176668713819931811048770190271"
"25267680276078003013678680992525463401061632866526"
"36270218540497705585629946580636237993140746255962"
"24074486908231174977792365466257246923322810917141"
"91430288197103288597806669760892938638285025333403"
"34413065578016127815921815005561868836468420090470"
"23053081172816430487623791969842487255036638784583"
"11487696932154902810424020138335124462181441773470"
"63783299490636259666498587618221225225512486764533"
"67720186971698544312419572409913959008952310058822"
"95548255300263520781532296796249481641953868218774"
"76085327132285723110424803456124867697064507995236"
"37774242535411291684276865538926205024910326572967"
"23701913275725675285653248258265463092207058596522"
"29798860272258331913126375147341994889534765745501"
"18495701454879288984856827726077713721403798879715"
"38298203783031473527721580348144513491373226651381"
"34829543829199918180278916522431027392251122869539"
"40957953066405232632538044100059654939159879593635"
"29746152185502371307642255121183693803580388584903"
"41698116222072977186158236678424689157993532961922"
"62467957194401269043877107275048102390895523597457"
"23189706772547915061505504953922979530901129967519"
"86188088225875314529584099251203829009407770775672"
"11306739708304724483816533873502340845647058077308"
"82959174767140363198008187129011875491310547126581"
"97623331044818386269515456334926366572897563400500"
"42846280183517070527831839425882145521227251250327"
"55121603546981200581762165212827652751691296897789"
"32238195734329339946437501907836945765883352399886"
"75506164965184775180738168837861091527357929701337"
"62177842752192623401942399639168044983993173312731"
"32924185707147349566916674687634660915035914677504"
"99518671430235219628894890102423325116913619626622"
"73267460800591547471830798392868535206946944540724"
"76841822524674417161514036427982273348055556214818"
"97142617910342598647204516893989422179826088076852"
"87783646182799346313767754307809363333018982642090"
"10848802521674670883215120185883543223812876952786"
"71329612474782464538636993009049310363619763878039"
"62184073572399794223406235393808339651327408011116"
"66627891981488087797941876876144230030984490851411"
"60661826293682836764744779239180335110989069790714"
"85786944089552990653640447425576083659976645795096"
"66024396409905389607120198219976047599490197230297"
"64913982680032973156037120041377903785566085089252"
"16730939319872750275468906903707539413042652315011"
"94809377245048795150954100921645863754710598436791"
"78639167021187492431995700641917969777599028300699"
"15368713711936614952811305876380278410754449733078"
"40789923115535562561142322423255033685442488917353"
"44889911501440648020369068063960672322193204149535"
"41503128880339536053299340368006977710650566631954"
"81234880673210146739058568557934581403627822703280"
"82616570773948327592232845941706525094512325230608"
"22918802058777319719839450180888072429661980811197"
"77158542502016545090413245809786882778948721859617"
"72107838435069186155435662884062257473692284509516"
"20849603980134001723930671666823555245252804609722"
"53503534226472524250874054075591789781264330331690";
    u_int8_t digitNumber[PB013_NBNUM*PB013_LENNUM] ;
    u_int32_t sum[PB013_LENNUM] ;
    int i ;
    pbR->nbClock = clock()  ;
    for(i=0;i<PB013_NBNUM*PB013_LENNUM;i++) {
        digitNumber[i] = strNumber[i] - '0' ;
    }
    for(i=0;i<PB013_LENNUM;i++){
        int j ;
        sum[i] = 0 ;
        for(j=0;j<PB013_NBNUM;j++) {
            sum[i] += digitNumber[i+j*PB013_LENNUM] ;
        }
    }
    for(i=PB013_LENNUM-1;i>0;i--) {
        sum[i-1] += sum[i]/10 ;
        sum[i] = sum[i] % 10 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose){
        fprintf(stdout,"\t bP%d  Sum=",pbR->pbNum);
        for(i=0;i<PB013_LENNUM;i++) fprintf(stdout,"%d",sum[i]);
        fprintf(stdout,"\n");
    }
    {
        int lg  = 0 ;
        for(i=0;lg < 10;i++) {
            lg+= sprintf(pbR->strRes+lg,"%d",sum[i]);
        }
        pbR->strRes[10] = '\0' ;
    }
    return 1 ;
}

typedef u_int32_t TY_SYR   ;


void SyracuseVal(TY_SYR *valCur) {
    TY_SYR val = valCur[0] ;
    if(val == 1) {
        return ;
    } else {
        if(val & 1) {
            valCur[1] = (3 * val + 1) ;
        }else {
            valCur[1] = val / 2 ;
        }
        SyracuseVal(valCur+1);
    }
}

// ne retourne que la longueur, pas les valeurs
#define SYRACUSE_MEM    1
#if SYRACUSE_MEM
#define SYRACUSE_MEM_SIZE (1<<19)
static u_int16_t lgS[SYRACUSE_MEM_SIZE] ;
static u_int16_t privSyracuseLgMem(u_int16_t lg,TY_SYR valInit,TY_SYR val) {
    if(val < SYRACUSE_MEM_SIZE ) {
        if(lgS[val]) {
            if(valInit < SYRACUSE_MEM_SIZE) lgS[valInit] = lg+lgS[val] ;
            return lg+lgS[val] ;
        } else if (val == 1) {
            lgS[1] = 1 ;
            if(valInit < SYRACUSE_MEM_SIZE) lgS[valInit] = lg+lgS[val] ;
            return lg+lgS[val] ;
        }
    }
    if(val == 1) {
        return lg;
    } else
    {
        if(val & 1) {
            val = (3 * val + 1) ;
        }else {
            val = val / 2 ;
        }
        return privSyracuseLgMem(lg+1,valInit,val);
    }
}
#else 
static u_int16_t privSyracuseLg(u_int16_t lg,TY_SYR val) {
    if(val == 1) {
        return lg+1;
    } else {
        if(val & 1) {
            val = (3 * val + 1) ;
        }else {
            val = val / 2 ;
        }
        return privSyracuseLg(lg+1,val);
    }
}
#endif

static u_int16_t SyracuseLg(TY_SYR val) {
#if SYRACUSE_MEM
    return privSyracuseLgMem(0,val,val) ;
#else
    return privSyracuseLg(0,val) ;
#endif
}

#define PB014_MAX_VALUE 1000000
#define PB014_PRINT 0
int PB014(PB_RESULT *pbR) {
    TY_SYR k ;
    TY_SYR kBest = 1 ;
    int bestLg = 0 ;
    pbR->nbClock = clock() ;
    for(k=1;k<PB014_MAX_VALUE;k++) {
        u_int16_t lg = SyracuseLg(k) ;
        if(lg > bestLg ) {
            bestLg = lg ;
            kBest = k ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d  lg=%d for %d\n"
            ,pbR->pbNum
            ,bestLg,kBest);
    sprintf(pbR->strRes,"%d",kBest);
#if PB014_PRINT
    {
        TY_SYR *chainVal = malloc((bestLg+1)*sizeof(chainVal[0])) ;
        chainVal[0] = kBest ;
        SyracuseVal(chainVal);
        
        {
            int i,lg ;
            for(i=0,lg=0 ; chainVal[i] != 1 ; i++) {
                printf("%d->",chainVal[i] );
                lg++ ;
            }
            printf("%d\n",chainVal[i]);
        }
        free(chainVal);
    }
#endif
    return 1 ;
}

int PB015(PB_RESULT *pbR) {
// answer = C(20,40) = 40! / ( 20! * 20!) = 21x22x...x40 / 1x2x3...x20
// on apparie les indices pairs en haut et les 10 derniers en bas (haut =2xbas)
//    = 21x23x25x27x...x39x 2**10 / 1x2x3x4...x10
// On fait les simplifications 25x2 = 5x10  ; 27=3x9 ; 21x2 = 6x7 ; 2**6=2x4x8 ;
// Il reste en haut 23x29x31x33x35x37x39x2**2
    sprintf(pbR->strRes,"%lld", ((u_int64_t) 23) * 29 * 31 * 33 * 35 *37 * 39 * 4 );
    return 1 ;
}



#define PB016_MAXL  1000/3
#define PB016_EXP   1000
int PB016(PB_RESULT *pbR) {
    u_int8_t digLarge[PB016_MAXL] ;
    int i,ie ,len = 0 ;
    u_int32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB016_MAXL;i++) {digLarge[i] = 0 ; }
    digLarge[0] = 1 ; len = 1 ;
    for(ie=0;ie<PB016_EXP;ie++) {
        for(i=0;i<len;i++) {
            digLarge[i] *= 2 ;
        }
        for(i=0;i<len;i++) {
            if(digLarge[i] >= 10) {
                digLarge[i] -= 10 ;
                digLarge[i+1]++ ;
            }
        }
        if(digLarge[len]) len++ ;
    }
    for(i=0;i<len;i++){
        S += digLarge[i] ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%d  Sumdig(2**%d)=%d\n"
            ,pbR->pbNum
            ,PB016_EXP,S) ;
    sprintf(pbR->strRes,"%d",S);
    return 1 ;
}

int PB016_gmp(PB_RESULT *pbR) {
    char  digLarge[PB016_MAXL] ;
    u_int32_t S = 0 ;
    int i ;
    pbR->nbClock = clock() ;
    mpz_t pow2 ;
    mpz_init(pow2);
    mpz_ui_pow_ui (pow2, 2, PB016_EXP);
    gmp_sprintf (digLarge,"%Zd",pow2);
    for(i=0;digLarge[i] != 0; i++) {
        S += digLarge[i] - '0' ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp Sumdig(2**%d)=%d\n"
            ,pbR->pbNum
            ,PB016_EXP,S) ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
}

#define PB018_SIZE  15
#define PB018_TSIZE ((PB018_SIZE*(PB018_SIZE+1))/2)
int PB018(PB_RESULT *pbR) {
    int vals[PB018_TSIZE] = {
    75
    ,95 ,64
    ,17 ,47 ,82
    ,18 ,35 ,87 ,10
    ,20 ,04 ,82 ,47 ,65
    ,19 ,01 ,23 ,75 ,03 ,34
    ,88 ,02 ,77 ,73 ,07 ,63 ,67
    ,99 ,65 ,04 ,28 ,06 ,16 ,70 ,92
    ,41 ,41 ,26 ,56 ,83 ,40 ,80 ,70 ,33
    ,41 ,48 ,72 ,33 ,47 ,32 ,37 ,16 ,94 ,29
    ,53 ,71 ,44 ,65 ,25 ,43 ,91 ,52 ,97 ,51 ,14
    ,70 ,11 ,33 ,28 ,77 ,73 ,17 ,78 ,39 ,68 ,17 ,57
    ,91 ,71 ,52 ,38 ,17 ,14 ,91 ,43 ,58 ,50 ,27 ,29 ,48
    ,63 ,66 ,04 ,68 ,89 ,53 ,67 ,30 ,73 ,16 ,69 ,87 ,40 ,31
    ,04 ,62 ,98 ,27 ,23 ,9 ,70 ,98 ,73 ,93 ,38 ,53 ,60 ,04 ,23 } ;
    int ic,ir ;
    pbR->nbClock = clock() ;
    // on commernce a l'avant derniere ligne
    for(ir=PB018_SIZE-2;ir>=0;ir--) {
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

#define PB020_MAXL  (100*2)
#define PB020_FACT   100
int PB020(PB_RESULT *pbR) {
    u_int32_t digLarge[PB020_MAXL] ;
    int i,ie ,len = 0 ;
    u_int32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB020_MAXL;i++) {digLarge[i] = 0 ; }
    digLarge[0] = 1 ; len = 1 ;
    for(ie=1;ie<=PB020_FACT;ie++) {
        for(i=0;i<len;i++) {
            digLarge[i] *= ie ;
        }
        for(i=0;i<len;i++) {
            if(digLarge[i] >= 10) {
                digLarge[i+1] += digLarge[i] /10 ;
                digLarge[i] = digLarge[i] % 10 ;
            }
        }
        while(digLarge[len]) {
            if(digLarge[len] >= 10) {
                digLarge[len+1] = digLarge[len] /10 ;
                digLarge[len] = digLarge[len] % 10 ;
            }
            len++ ;
        }
        
    }
    for(i=0;i<len;i++){
        S += digLarge[i] ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d  Sumdig(%d!)=%d\n",pbR->pbNum, PB020_FACT,S) ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
}


#define PB021_MAXM  10000
int PB021(PB_RESULT *pbR) {
    int SumDiv[PB021_MAXM+1] ;
    int i ;
    u_int32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=1;i<=PB021_MAXM;i++) {
        SumDiv[i] = 0 ;
    }
    for(i=1;2*i<PB021_MAXM;i++) {
        int mxi ;
        for(mxi=2*i;mxi <= PB021_MAXM ; mxi += i) {
            SumDiv[mxi] += i ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
    for(i=1;i<PB021_MAXM;i++){
        if((SumDiv[i] > i) && (SumDiv[i] <= PB021_MAXM)) {
            if(SumDiv[SumDiv[i]] == i ) {
                fprintf(stdout,"%d<->%d ",i,SumDiv[i]) ;
                S += i + SumDiv[i] ;
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
}

#define PB023_MAXM  28123
int PB023(PB_RESULT *pbR) {
    int SumDiv[PB023_MAXM+1] ;
    int *isSumAbun ;
    int abundant[PB023_MAXM] ;
    int nbAbund =0 ;
    int i ;
    u_int32_t SumNo2a = 0 ;
    int nbNoSum2a = 0 ;
    pbR->nbClock = clock()  ;
    for(i=1;i<=PB023_MAXM;i++) {
        SumDiv[i] = 0 ;
    }
    for(i=1;2*i<=PB023_MAXM;i++) {
        int mxi ;
        for(mxi=2*i;mxi <= PB023_MAXM ; mxi += i) {
            SumDiv[mxi] += i ;
        }
    }
    isSumAbun = SumDiv ;
    for(i=1;i<=PB023_MAXM;i++){
        if(SumDiv[i] > i) {
            abundant[nbAbund++] = i ;
        }
        isSumAbun[i] = 0 ;
    }
    for(i=0;i<nbAbund;i++) {
        int j ;
        for(j=i;j<nbAbund;j++) {
            int S = abundant[i]+abundant[j] ;
            if(S<= PB023_MAXM) {
                isSumAbun[S] = 1;
            } else {
                break ;
            }
        }
    }
    for(i=1;i<=PB023_MAXM;i++){
        if(isSumAbun[i] == 0) {
            SumNo2a += i ;
            nbNoSum2a++ ;
        }
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d  nbNo2ab=%d ,Sum=%d\n",pbR->pbNum,nbNoSum2a,SumNo2a) ;
    sprintf(pbR->strRes,"%d",SumNo2a);
    return 1 ;
}


#define PB024_SIZE  10
#define PB024_RANG   1000000
int PB024(PB_RESULT *pbR) {
    int i  ;
    u_int64_t factoriels[PB024_SIZE] ;
    u_int8_t digits[PB024_SIZE] ;
    int rang = PB024_RANG - 1 ;
    pbR->nbClock = clock() ;
    factoriels[0] = 1 ;
    for(i=1;i<PB024_SIZE;i++) {
        factoriels[i] = factoriels[i-1] *i ;
    }
    for(i=0;i<PB024_SIZE;i++) { digits[i] = i ; }
    for(i=0;i<PB024_SIZE;i++) {
        int rgUnused = rang / factoriels[PB024_SIZE-1-i] ;
        rang -= rgUnused * factoriels[PB024_SIZE-1-i] ;
        if(rgUnused>0){
            int j ;
            u_int8_t saveDig = digits[i+rgUnused] ;
            for(j=rgUnused+i; j > i ; j--) {
                digits[j] = digits[j-1] ;
            }
            digits[i] = saveDig ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    for(i=0;i<PB024_SIZE;i++)pbR->strRes[i] = '0' + digits[i];
    pbR->strRes[10] = 0 ;
    return 1 ;

 }

#define PB025_MAXDIGIT 1000
int PB025(PB_RESULT *pbR) {
    u_int8_t FA[PB025_MAXDIGIT+1] ;
    u_int8_t FB[PB025_MAXDIGIT+1] ;
    int i,n,len ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB025_MAXDIGIT+1;i++) {
        FA[i]=FB[i]= 0 ;
    }
    FA[0] = FB[0] = 1 ;
    len = 1 ;
    n=2 ;
    do {
        u_int8_t *FD = FA ;
        if(++n & 1) FD =FB ;
        for(i=0;i<len;i++) {
           if(FA[i]+FB[i]< 10) {
               FD[i] = FA[i]+FB[i] ;
           } else {
               FD[i] = FA[i]+FB[i]-10 ;
               FD[i+1]++ ;
           }
        }
        if(FD[len]) len++ ;
    } while(len < PB025_MAXDIGIT) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d len>=%d pour Fibonnaci n°%d \n",pbR->pbNum , PB025_MAXDIGIT,n) ;
    sprintf(pbR->strRes,"%d",n);
    return 1 ;
}

#define PB026_MAX 1000
int PB026(PB_RESULT *pbR) {
    int maxCycle = 0 ;
    int dmax = 0 ;
    int d ;
    u_int16_t isReste[PB026_MAX] ;
    pbR->nbClock = clock() ;
    for(d=2;d<PB026_MAX;d++) {
        int r = 1;
        int lc = 1 ;
        memset(isReste,0,d*sizeof(isReste[0]));
        isReste[1] = 1 ;
        while(1) {
            r *= 10 ; lc++ ;
            if(r >= d) {
                r = r % d ;
                if(r==0) break ;
                if( isReste[r]) {
                    if(lc - isReste[r] > maxCycle) {
                        maxCycle = lc - isReste[r] ;
                        dmax = d ;
                    }
                    break ;
                } else {
                    isReste[r] = lc ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d Pour 1/%d cycle de %d \n",pbR->pbNum, dmax,maxCycle) ;
    sprintf(pbR->strRes,"%d",dmax) ;
    return 1 ;
}

#define PB027_MAX_A 1000 // doit etre pair
#define PB027_MAX_B 1000

int32_t CmpPrime(const void *p1,const void *p2) {
    return *((int32_t *)p1) - *((int32_t *)p2) ;
}

int PB027(PB_RESULT *pbR) {
    // b est positif et est un nbre premier (P(0)=B)
    // a est impair (P(1)-P(0) =  1 + a
    // P(b) n'est pas premier donc n<b
    // si P(b-a)=(b-a)(b-a)+a)+b n'est pas premier. Donc si a<b, n<b-a
    // n<b => P(n) < (b**2 + b*a+b < Max_b**2 + Max_b*Max_a
    pbR->nbClock = clock()  ;
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB027_MAX_B*(PB027_MAX_A+PB027_MAX_B)) ;
    if(ctxP == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    {
        int32_t ib = 0;
        int32_t maxNbprime = 0 ;
        int32_t  a , b ;
        int32_t aMax = 0 , bMax = 0 ;
        while ( (b=ctxP->tbPrime[ib++]) < PB027_MAX_B) {
            for ( a= (- PB027_MAX_A+1 ); a < PB027_MAX_A ;  a+=2) {
                int nbPrime = 0 ;
                int32_t n = 0 ;
                while(1) {
                    int32_t pn = n*(n+a)+b ;
                    if(pn<=0) { break ;}
                    if(bsearch(&pn,ctxP->tbPrime,ctxP->nbPrime,sizeof(pn),CmpPrime) == NULL) break ;
                    nbPrime ++ ;
                    n++ ;
                }
                if(nbPrime > maxNbprime) {
                    maxNbprime = nbPrime ;
                    aMax = a;
                    bMax = b ;
                 }
            }
        }
        pbR->nbClock = clock() - pbR->nbClock ;

        if(pbR->isVerbose)fprintf(stdout,"\t PB%d %d nxn%c%dxn+%d produce %d premiers\n"
                ,pbR->pbNum
                ,aMax*bMax
                ,(aMax > 0) ? '+' : '-'
                ,(aMax > 0) ? aMax : -aMax
                , bMax, maxNbprime) ;
        sprintf(pbR->strRes,"%d",aMax*bMax) ;
        
    }
    Free_tablePrime(ctxP);
    return 1;
}

int PB028(PB_RESULT *pbR) {
    /*
     43 44 45 46 47 48 49
     42 21 22 23 24 25 26
     41 20  7  8  9 10 27
     40 19  6  1  2 11 28
     39 18  5  4  3 12 29
     38 17 16 15 14 13 30
     37 36 35 34 33 32 31
     
     Quite easy to demonstrate that the answer is a polynom of degree 3 of the half size
     (as the diagonal is odd squares (1 9 25 49 ...) and the other deduce by addition of a * p + b.
     So if size is 2k+1
     P(0) = 1
     P(1) = P(0)+3+5+7+9 = 25
     P(2) = P(1)+13+17+21+25 = 101
     P(3) = P(2)+31+37+43+49 = 261
     
     P(k) = A k**3+ B k**2+  C *p + D
     D = 1 ;
     (1) A+B+C = 24
     (2) 8A+4B+2C= 100 <=> (2') 4A + 2B + C = 50
     (3) 27A+9B+3C = 260
     
     (3)-4*(2)+(1) <=> 6A=260 - 3 * 100 + 3 * 24 = 32
     A = 16/3
     (2')-(1) 3A + B = 26 => B =10
     C =26/3
     
     donc P(k) = (16 k**3 + 30 k**2 + 26 k + 3) / 3
     */
    int k=500 ;
   
    sprintf(pbR->strRes,"%d",(16*k*k*k + 30*k*k + 26*k + 3) / 3) ;
    return 1 ;
    
}



int PB029(PB_RESULT *pbR) {
    // les seules valeurs qui peuvent avoir des puissances communes sont 2, 3, 5, 6, 7 , 10
    // pour 2 il y a les puissances de 2, 4, 8,16,32,64 <100 soit exp multiple de 1,2,3,4,5,6
    // pour 3 il y a les puissances de 3 9 27 81 soit exp multiple de 1,2,3,4
    // pour 5 => 5 et 25 exp mult de 1,2
    // pour 6 => 6 et 36 exp mult de 1,2
    // pour 7 => 7 et 49 exp mult de 1,2
    // pour 10 => 10 et 100 exp muult de 1,2
    // 4,8,9 sont deja traites des 2 ou 3
    // on va compter les collisions et le total sera 99x99 - nbColl
    int nbCol = 0 ;
    int p ;
    pbR->nbClock = clock();
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
    for(p=2;p<=3;p++) {
        int expM ; int n ; // p**expM <= 100
        for(n=p,expM=0; n < 100; expM++) {
            n *= p ;
        }
        {
            int i, exp ;
            u_int8_t *isOccupied = calloc(expM*101,sizeof(isOccupied[0])) ;
            for(i=2;i<=100;i++) { isOccupied[i] = 1 ; }
            for(exp=2;exp <= expM; exp++) {
                for(i=2;i<=100;i++) {
                    if(isOccupied[i*exp]) nbCol++ ;
                    else isOccupied[i*exp] = 1 ;
                }
            }
            free(isOccupied) ;
            if(pbR->isVerbose)fprintf(stdout,"(%d,%d,%d)",p,expM,nbCol) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n");
    // pour p = 5 ,6 ,7 , 10 nbCol = 49 ;
    nbCol += 49 * 4 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",99*99 - nbCol) ;
    return 1 ;
}

int PB030(PB_RESULT *pbR) {
    int32_t pow5[10] ;
    int i, n ;
    int SumRes = 0;
    int Maxn = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<10;i++) {
        pow5[i] = i*i*i*i*i ;
    }
    Maxn = 6 * pow5[9] +1 ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d sum[ ",pbR->pbNum) ;
    for(n=2;n<Maxn;n++) {
        if(n == pow5[n % 10] + pow5[ (n/10) % 10] + pow5[(n/100) % 10] + pow5[(n/1000) % 10] + pow5[(n/10000) % 10] + + pow5[(n/100000)]) {
            SumRes += n ;
            if(pbR->isVerbose) fprintf(stdout,"%d ",n);
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"] = %d\n",SumRes) ;
    sprintf(pbR->strRes,"%d",SumRes);
    return 1 ;
}

#define PB031_NBV   8
int PB031(PB_RESULT *pbR) {
    int value[PB031_NBV] = { 200 , 100 , 50, 20 , 10 ,5 , 2 ,1};
    int nbPieces[PB031_NBV] ;
    int nbDecomp = 0 ;
    int np , sum ;
    pbR->nbClock = clock() ;
    np = -1 ;
    sum = 200 ;
    // on parcourt la decomposition dans l'ordre lexicographique
    while(1) {
        if(sum == 0) {
            nbDecomp++ ;
//          { int i ;    for(i=0;i<PB031_NBV;i++) { printf("%c%dx%d",(i==0) ? '\n' : '+' ,nbPieces[i],value[i]) ; }   }
            while(sum < value[np]) {
                sum += nbPieces[np] * value[np] ;
                nbPieces[np] = 0 ;
                if(--np < 0 ) break ;
            }
            if(np < 0) {
                break ;
            }
            nbPieces[np]++ ;
            sum -= value[np] ;
        }
        while(++np < PB031_NBV-1) {
            nbPieces[np] = 0 ;
        }
        // derniere piece de 1
        nbPieces[np] = sum ;
        sum = 0 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbDecomp) ;
    return 1 ;
}

void HeapSortUint8(u_int8_t *H,int n) {
    int i;
    for(i=n-1;i>=0;i--) {
        int i1 ;
        for(i1=(i-1)/2;i1>=0;i1--) {
            int k = i1 ;
            u_int8_t heap = 0 , v=H[k] ;
            while(heap==0 && (2*k)< i) {
                int j=2*k+1;
                if( (j<i) && (H[j]<H[j+1]) ) { j++;  }
                if(v>=H[j]) {
                    heap=1;
                } else {
                    H[k]=H[j];
                    k=j;
                }
            }
            H[k]=v;
        }
        { u_int8_t  t=H[0];  H[0]=H[i];   H[i]=t; }
    }
}


void HeapSortUint8Rev(u_int8_t *H,int n) {
    int i;
    for(i=n-1;i>=0;i--) {
        int i1 ;
        for(i1=(i-1)/2;i1>=0;i1--) {
            int k = i1 ;
            u_int8_t heap = 0 , v=H[k] ;
            while(heap==0 && (2*k)< i) {
                int j=2*k+1;
                if( (j<i) && (H[j]>H[j+1]) ) { j++;  }
                if(v<=H[j]) {
                    heap=1;
                } else {
                    H[k]=H[j];
                    k=j;
                }
            }
            H[k]=v;
        }
        { u_int8_t  t=H[0];  H[0]=H[i];   H[i]=t; }
    }
}



// return -1 if last permutation, or the lower index of modification
static int NextPermut(u_int8_t *perm,int lg) {
    int i ;
    for(i=(--lg - 1);i>=0 && perm[i]>perm[i+1];i--) ;
    if(i<0) return i ;
    { int j ;
        u_int8_t tmp = perm[i] ;
        for(j=lg;perm[j]<tmp;j--) ;
        perm[i] =perm[j] ;
        perm[j] = tmp ;
        for(j=i+1;j<lg;j++,lg--) {
            tmp = perm[j] ; perm[j] =perm[lg] ; perm[lg] = tmp ;
        }
        return i ;
    }
}

// force the change of rg (next value possible), return -1 if impossible
static int NextPermutRg(u_int8_t *perm,int lg,int rg) {
    if(rg+1 >=lg || rg == 0) return -1 ;
    HeapSortUint8Rev(perm+rg+1,lg-rg-1) ;
    return NextPermut(perm,lg);
}


int ChkPermutRg(u_int8_t *perm,int lg,int rg) {
    int i ;
    for(i=(lg - 1);i>rg && perm[i]<perm[rg];i--) ;
    return (i > rg) ;
}


// idem ordre reverse
// return -1 if last permutation, or the lower index of modification
static int NextPermutRev(u_int8_t *perm,int lg) {
    int i ;
    for(i=(--lg - 1);i>=0 && perm[i]<perm[i+1];i--) ;
    if(i<0) return i ;
    { int j ;
        u_int8_t tmp = perm[i] ;
        for(j=lg;perm[j]>tmp;j--) ;
        perm[i] =perm[j] ;
        perm[j] = tmp ;
        for(j=i+1;j<lg;j++,lg--) {
            tmp = perm[j] ; perm[j] =perm[lg] ; perm[lg] = tmp ;
        }
        return i ;
    }
}

int NextPermutRgRev(u_int8_t *perm,int lg,int rg) {
    if(rg+1 >=lg || rg == 0) return -1 ;
    HeapSortUint8(perm+rg+1,lg-rg-1) ;
    return NextPermutRev(perm,lg);
}



int PB032(PB_RESULT *pbR) {
    u_int8_t dg[9] ;
    int Sum = 0 ;
    pbR->nbClock = clock() ;
    { int i ; for(i=0;i<9;i++) {  dg[i] = i+1 ; } }
    int rg = 0 ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
    do {
        int nb4 = 1000*dg[0]+100*dg[1]+10*dg[2]+dg[3]  ;
        do {
        // on verifie si decomposition OK ; c'est forcement les 4 premiers = dddd = d x dddd ou dddd = dd x ddd
            if(nb4 == dg[4] * (1000*dg[5]+100*dg[6]+10*dg[7]+dg[8]) ) {
               Sum += nb4 ;
               if(pbR->isVerbose)fprintf(stdout," %d=%dx%d%d%d%d ",nb4,dg[4],dg[5],dg[6],dg[7],dg[8]) ;
               while((rg = NextPermut(dg,9) ) >= 4) ;
                
            } else if(nb4 == (10*dg[4]+dg[5]) * (100*dg[6]+10*dg[7]+dg[8]) ) {
               Sum += nb4 ;
               if(pbR->isVerbose)fprintf(stdout," %d=%d%dx%d%d%d ",nb4,dg[4],dg[5],dg[6],dg[7],dg[8]) ;
               while((rg = NextPermut(dg,9) ) >= 4) ;
           } else {
               rg=NextPermut(dg,9) ;
            }
        } while ( rg >= 4) ;
    } while(rg >= 0 ) ;
    if(pbR->isVerbose)fprintf(stdout,"\n");
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum) ;
    return 1 ;
}

int PB033(PB_RESULT *pbR) {
    int a,n,d ;
    int Num = 1;
    int Den = 1;
    pbR->nbClock = clock();
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
    for(n=1;n<9;n++) {
        for(d=n+1;d<10;d++) {
            // 2 cas (an/da = n/d ou na/ad = n/d)
            for(a=1;a<10;a++) {
                if ( (10*a+n)*d == (10*d+a)*n ) {
                    Num *= n ; Den *= d ;
                    if(pbR->isVerbose)printf(" %d%d/%d%d=%d/%d",a,n,d,a,n,d) ;
                } else if ( (10*n+a)*d == (10*a+d)*n ) {
                    Num *= n ; Den *= d ;
                    if(pbR->isVerbose)printf(" %d%d/%d%d=%d/%d",n,a,a,d,n,d) ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    {
        int pgcd = PGCD(Num,Den);
        Num /= pgcd ; Den /=pgcd ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d Prod=%d/%d\n",pbR->pbNum,Num,Den) ;
    sprintf(pbR->strRes,"%d",Den);
    return 1 ;
}

// max value = 9! x 7 ;
#define PB034_MAX 2540160
#define PB034_MAXDIG    7
int PB034(PB_RESULT *pbR) {
    int fact[10] , diffFact[10] ;
    int dig[PB034_MAXDIG] ;
    int i, n , k , nbDig ;
    int  Diff,Sum = 0 ;
    pbR->nbClock = clock()  ;
    fact[0] = 1 ;
    for(i=1;i<=9;i++) {
        fact[i] =fact[i-1] * i ;
    }
    for(i=1;i<=9;i++) {
        diffFact[i] = fact[i] - fact[i-1] ;
    }
    diffFact[0] = fact[0] - fact[9] ;
    for(i=0;i<PB034_MAXDIG;i++) { dig[i] = 0 ; }
    n = 3 ;
    dig[0] = n ;
    Diff = fact[dig[0]] - n ;
    nbDig = 1 ;
    k = 0 ;
   if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
   while (n < PB034_MAX){
         if(Diff == 0) {
            Sum += n ;
             if(pbR->isVerbose) {
                 fprintf(stdout," %d",n); for(i=nbDig-1;i>=0;i--) { fprintf(stdout,"%c%d!",(i==nbDig-1) ? '=' : '+',dig[i]) ; }
             }
        }
        while(dig[k] == 9) {
            dig[k++] = 0;
            Diff += diffFact[0] ;
        }
        if(k >= nbDig) {
            dig[k] = 1 ;
            nbDig++ ;
        } else {
            dig[k]++ ; Diff += diffFact[dig[k]]  - 1;
        }
        k = 0 ;
        n++ ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum) ;
    return 1 ;
    
}

#define PB035_MAXN  1000000
int PB035(PB_RESULT *pbR) {
    int i , puis10 ;
    int nbFind = 0 ;
    u_int8_t    *isPrimeUsed = NULL ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB035_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    } else {
        isPrimeUsed = calloc(ctxP->nbPrime,sizeof(isPrimeUsed[0])) ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d",pbR->pbNum) ;
    for(i=0,puis10=1;i<ctxP->nbPrime;i++) {
        T_prime n = ctxP->tbPrime[i] ;
        if(isPrimeUsed[i] == 0) { // on met a zero ceux deja trouve
            int lgCyclePrime = 1 ;
            isPrimeUsed[i] = 1 ;
            T_prime nr = n ;
            while(10*puis10 <= n ) {
                puis10 *= 10 ;
            }
            // si dans le cycle on aboutit a une valeur inferieure c'est deja crame
            while ((nr = (nr % 10) * puis10  + nr /10) > n) {
                T_prime *ptNr ;
                // il faur nr premier
                if( (ptNr=bsearch(&nr,ctxP->tbPrime,ctxP->nbPrime,sizeof(nr),CmpPrime)) != NULL) {
                        lgCyclePrime++ ;
                        isPrimeUsed[ptNr-ctxP->tbPrime] = 1 ;
                } else {
                   lgCyclePrime = 0 ;
                   break ;
                }
            }
            
            if(nr == n) {
                if(pbR->isVerbose)fprintf(stdout," %d",n) ;
                nbFind += lgCyclePrime ;
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbFind);
        free(isPrimeUsed) ;
    return 1 ;
}

static int IsPal2(u_int32_t n) {
    u_int32_t B=1, H = 1 << 31 ;
    if(n<=1) return 1 ;
    while(H>n) { H >>= 1 ; }
    while(H>B) {
        u_int32_t cmp = (H+B) & n ;
        if(cmp == H || cmp == B) return 0  ;
        H >>= 1;
        B <<= 1 ;
    }
    return 1 ;
}

int PB036(PB_RESULT *pbR) {
    int i1,i2,i3 ;
    int Sum = 0 ;
    int nbFind = 0 ;
    pbR->nbClock = clock()  ;
      if(pbR->isVerbose)fprintf(stdout,"\t PB%d ",pbR->pbNum) ;
    for(i1=1;i1<10;i1++) {
        int n=i1;
        if(IsPal2(n)) {
            if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
            nbFind++ ;
            Sum += n;
        }
        n=i1*11 ;
        if(IsPal2(n)) {
            if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
            nbFind++ ;
            Sum += n;
        }
    }
    for(i1=1;i1<10;i1++) {
        for(i2=0;i2<10;i2++) {
            int n=i1*101 + i2*10 ;
            if(IsPal2(n)) {
                nbFind++ ;
                if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                Sum += n;
            }
            n= i1*1001 + i2 * 110 ;
            if(IsPal2(n)) {
                nbFind++ ;
                if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                Sum += n;
            }
        }
    }
    for(i1=1;i1<10;i1++) {
        for(i2=0;i2<10;i2++) {
            for(i3=0;i3<10;i3++) {
                int n=i1*10001 + i2*1010 + i3 * 100 ;
                if(IsPal2(n)) {
                    if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                    nbFind++ ;
                    Sum += n;
                }
                n= i1*100001 + i2 * 10010 + i3 * 1100 ;
                if(IsPal2(n)) {
                    if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                    nbFind++ ;
                    Sum += n;
                }
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d Nbfind=%d\n",pbR->pbNum,nbFind) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum);
    return 1 ;
}

#define PB037_MAXP  1000000
#define PB037_MAX_NB    11
int PB037(PB_RESULT *pbR) {
    int  truncPrime[PB037_MAX_NB] ;
    int nbAlloc = 1024 ;
    int *truncRight = malloc(nbAlloc*sizeof(truncRight[0])) ;
    int *truncLeft = malloc(nbAlloc*sizeof(truncLeft[0])) ;
    int nbRight = 0 , nbLeft = 0;
    int i,k,indm,indm1 = 0 , indrm, indrm1 = 0 , indlm , indlm1 = 0;
    int nbFind = 0 ;
    int puis10 = 10 ;
    int Sum = 0 ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB037_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    // on tabule ceux a 2 digits (composition de 2,3,5,7)
    truncPrime[nbFind] = 23 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 37 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 53 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 73 ; Sum += truncPrime[nbFind++] ;
    // on tbule les right a 2 digits
    truncRight[nbRight++] = 23 ;
    truncRight[nbRight++] = 29 ;
    truncRight[nbRight++] = 31 ;
    truncRight[nbRight++] = 37 ;
    truncRight[nbRight++] = 53 ;
    truncRight[nbRight++] = 59 ;
    truncRight[nbRight++] = 71 ;
    truncRight[nbRight++] = 73 ;
    truncRight[nbRight++] = 79 ;
    // on tabule les left a 2 digits
    truncLeft[nbLeft++] = 13 ;
    truncLeft[nbLeft++] = 17 ;
    truncLeft[nbLeft++] = 23 ;
    truncLeft[nbLeft++] = 37 ;
    truncLeft[nbLeft++] = 43 ;
    truncLeft[nbLeft++] = 47 ;
    truncLeft[nbLeft++] = 53 ;
    truncLeft[nbLeft++] = 73 ;
    truncLeft[nbLeft++] = 83 ;
    truncLeft[nbLeft++] = 97 ;
    do {
        indm = nbFind ;
        indrm = nbRight ;
        indlm = nbLeft ;
        puis10 *= 10 ;
        if(puis10*10 > PB037_MAXP) {
            free(truncRight); free(truncLeft) ;
            if(pbR->isVerbose) {
                fprintf(stdout,"\n\t PB%d ",pbR->pbNum) ;
                for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
            }
            return 0 ;
        }
        if(4*(indrm-indrm1)+ nbRight >= nbAlloc || 4*(indlm-indlm1)+ nbLeft >= nbAlloc) {
            nbAlloc *= 2 ;
            truncRight = realloc(truncRight,nbAlloc*sizeof(truncRight[0])) ;
            truncLeft = realloc(truncLeft,nbAlloc*sizeof(truncLeft[0])) ;
        }
        for(i=indrm1;i<indrm;i++) {
            // on va le completer avec des 1, 3, 7, 9
            u_int32_t p = 10*truncRight[i]+1 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p += 3-1 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p+= 7-3 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p += 9 - 7 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
        }
        for(k=1;k<=9;k++) {
            for(i=indlm1;i<indlm;i++) {
            // on va le completer avec de 1 a 9
                u_int32_t p = truncLeft[i] + k * puis10  ;
                if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncLeft[nbLeft++] = p ; }
            }
        }
        for(i=indrm;i<nbRight;i++) {
            u_int32_t p = truncRight[i] ;
            if(bsearch(&p,truncLeft+indlm,nbLeft-indlm,sizeof(p),CmpPrime)) {
                 truncPrime[nbFind++] = truncRight[i] ;
                Sum += truncRight[i]  ;
            }
        }
        indm1 = indm ;
        indrm1 = indrm ;
        indlm1 = indlm ;
    } while(nbFind < PB037_MAX_NB ) ;
    if(pbR->isVerbose) {
        fprintf(stdout,"\t PB%d right=%d left=%d RL=%d\n",pbR->pbNum,nbRight,nbLeft,nbFind) ;
        fprintf(stdout,"\t PB%d ",pbR->pbNum) ;
        for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum);
    free(truncRight);
    return 1 ;
}

int PB037r(PB_RESULT *pbR) {
    int  truncPrime[PB037_MAX_NB] ;
    int nbAlloc = 1024 ;
    int *truncRight = malloc(nbAlloc*sizeof(truncRight[0])) ;
    int nbRight = 0 ;
    int i,indm,indm1 = 0 , indrm, indrm1 = 0 ;
    int nbFind = 0 ;
    int puis10 = 10 ;
    int Sum = 0 ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB037_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    // on tabule ceux a 2 digits (composition de 2,3,5,7)
    truncPrime[nbFind] = 23 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 37 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 53 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 73 ; Sum += truncPrime[nbFind++] ;
    // on tbule les right a 2 digits
    truncRight[nbRight++] = 23 ;
    truncRight[nbRight++] = 29 ;
    truncRight[nbRight++] = 31 ;
    truncRight[nbRight++] = 37 ;
    truncRight[nbRight++] = 53 ;
    truncRight[nbRight++] = 59 ;
    truncRight[nbRight++] = 71 ;
    truncRight[nbRight++] = 73 ;
    truncRight[nbRight++] = 79 ;
    do {
        indm = nbFind ;
        indrm = nbRight ;
        puis10 *= 10 ;
        if(puis10*10 > PB037_MAXP) {
            free(truncRight);
            if(pbR->isVerbose) {
                fprintf(stdout,"\n\t PB%d ",pbR->pbNum) ;
                for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
            }
            return 0 ;
        }
        if(4*(indrm-indrm1)+ nbRight >= nbAlloc ) {
            nbAlloc *= 2 ;
            truncRight = realloc(truncRight,nbAlloc*sizeof(truncRight[0])) ;
        }
        for(i=indrm1;i<indrm;i++) {
            // on va le completer avec des 1, 3, 7, 9
            u_int32_t p = 10*truncRight[i]+1 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p += 3-1 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p+= 7-3 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
            p += 9 - 7 ;
            if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime)) { truncRight[nbRight++] = p ; }
        }
        for(i=indrm;i<nbRight;i++) {
            u_int32_t k = puis10 ;
            u_int32_t p = truncRight[i] ;
            while( k >= 10) {
                p = p % k ;
                if(bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime) == NULL)  break ;
                k /= 10 ;
            }
            if ( k == 1) {
                truncPrime[nbFind++] = truncRight[i] ;
                Sum += truncRight[i]  ;
            }
        }
        indm1 = indm ;
        indrm1 = indrm ;
    } while(nbFind < PB037_MAX_NB ) ;
    if(pbR->isVerbose) {
        fprintf(stdout,"\t PB%d right=%d RL=%d\n",pbR->pbNum,nbRight,nbFind) ;
        fprintf(stdout,"\t PB%d ",pbR->pbNum) ;
        for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",Sum);
    free(truncRight);
    return 1 ;
}


// max value = 9! x 7 ;

int PB038(PB_RESULT *pbR) {
    u_int8_t dig[10] ;
    int i , k ;
    int maxP = 0 ;
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%d (1,2)",pbR->pbNum) ;
    // XXXX par 1,2 XXXX =>4  XXXXx2 =>5
    for(i=5000;i<=10000;i++) {
        int n1 = i * 100002 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout," %d->%d",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d (1,2,3)",pbR->pbNum) ;
    // XXX par 1 , 2 , 3  XXX => 3c XXXx2=>3c XXXx3=>3c
    for(i=100;i<=333;i++) {
        int n1 = i * 1002003 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d (1,2,3,4)",pbR->pbNum) ;
    // XX par 1,2,3,4  XX=>2c XXx2=>2c XXx3=>2c XXx4=>3c
    for(i=25;i<=33;i++) {
        int n1 = i * 10203004 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%d (1,2,3,4,5)",pbR->pbNum) ;
    // X par 1,2,3,4,5  X=>1c  Xx2=>2c Xx3=>2c Xx4=>2c Xx4=>2c
    for(i=1;i<=9;i++) {
        int n1 = i * 102030405 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }

    
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",maxP) ;
    return 1 ;
    
}

#define PB039_MAX   1000

int PB039(PB_RESULT *pbR) {
    // pour m > n > 0
    // a = (m**2 - n**2)*k ; b = 2*m*n*k ; c = (m**2 + n**2)*k ;
    // donc S = 2 * (m**2 + m*n)*k = 2 * m * (m+n) * k
    int32_t m , n , P , nMax =0 , Pbest=1 ;
    int32_t *nbSol = calloc(PB039_MAX/2+1,sizeof(nbSol[0])) ;
    pbR->nbClock = clock()  ;

    for(n=1;n*(2*n+1)<=PB039_MAX/2;n++) {
        for(m=n+1;(P=m*(m+n))<=PB039_MAX/2;m++) {
            if( PGCD(m,n)== 1 && (((m & n) & 1) == 0 )) {
                do {
                    nbSol[P]++ ;
                    P += m*(m+n) ;
                } while(P <= PB039_MAX/2 ) ;
            }
        }
    }
    for(n=1;n<PB039_MAX/2;n++) {
        if(nbSol[n]> nMax) {
            nMax = nbSol[n] ;
            Pbest = n ;
        }
    }
    if(pbR->isVerbose) {
        int mmn ;
        fprintf(stdout,"\t PB%03d P=%d",pbR->pbNum,2*Pbest);
        for(mmn=2; mmn<Pbest;mmn++) {
            if(nbSol[mmn] == 1 && (Pbest % mmn) == 0) {
                int k = Pbest / mmn ;
                for(m=2;m*m< mmn; m++) {
                    if((mmn % m ) == 0) {
                        n = mmn / m  - m ;
                        if(m>n && PGCD(m,n)== 1 && (((m & n) & 1) == 0 )) {
                            fprintf(stdout," (%d,%d,%d)",k*(m*m - n*n) ,k*2*m*n,k*(m*m + n*n));
                        }
                    }
                }
            }
        }
        fprintf(stdout,"\n");
    }
    pbR->nbClock = clock() -  pbR->nbClock ;

    sprintf(pbR->strRes,"%d",2*Pbest) ;
    free(nbSol);
    return 1;
}


#define PB040_MAXN  7
#define PB040_IND_MAX   1000000
int PB040(PB_RESULT *pbR) {
    int i,nb,P=1,ind;
    u_int32_t maxByDig[PB040_MAXN+1] ;
    u_int32_t offsByDig[PB040_MAXN+1] ;
    pbR->nbClock = clock()  ;
    maxByDig[0] = 1 ; nb = 9 ; offsByDig[0] = -1 ;
    for(i=1;i<=PB040_MAXN;i++) {
        maxByDig[i] = maxByDig[i-1] + nb * i ;
        nb *= 10 ;
        offsByDig[i] = (offsByDig[i-1]+1)*10 ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%03d P",pbR->pbNum);
        
    for(ind=1;ind<=PB040_IND_MAX;ind *=10) {
        int n, r ;
        for(i=1;maxByDig[i]<= ind;i++) ;
        n = (ind+offsByDig[i]) / i ;
        r = (ind+offsByDig[i]) % i ;
        while(++r < i) {
            n /= 10 ;
        }
        P *= n % 10 ;
        if(pbR->isVerbose) fprintf(stdout,"%c%d",(ind==1) ? '=' : 'x' , n % 10) ;
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\n");
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",P);
    return 1 ;
}


// on commence a 7
// car permut de 8 et 9 sont multiples de 3
#define PB041_MAXN  7654321
int PB041(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    u_int8_t permut[7] = {7,6,5,4,3,2,1} ;
    T_prime pBest = 0 ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB041_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    do {
        if((permut[6] & 1) && (permut[6] != 5)) { // pas pair ni multiple de 5
            T_prime p = permut[0]*1000000+permut[1]*100000+permut[2]*10000+permut[3]*1000+permut[4]*100+permut[5]*10+permut[6] ;
            if( bsearch(&p,ctxP->tbPrime,ctxP->nbPrime,sizeof(p),CmpPrime) != NULL) {
                pBest = p ;
                break ;
            }
        }
    } while(NextPermutRev(permut,7) >= 0 ) ;
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pBest) {
        sprintf(pbR->strRes,"%d",pBest);
        return 1 ;
    } else {
        return 0 ;
        
    }
}

#include "p042_words.h"
int PB042(PB_RESULT *pbR) {
    const char ** P042_W = P042_GetData() ;
    const char *word ;
    int nbTriangle = 0 ;
    pbR->nbClock = clock()  ;
    while((word = *P042_W++) != NULL ) {
        int32_t val=0 ;
        int32_t n ;
        u_int8_t c ;
        while((c=*word++)) {
            val += c - 'A' + 1 ;
        }
        val *= 2 ;
        n = (int32_t)Sqrt64(val) ;
        if(val == n*(n+1)) {
            nbTriangle++ ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbTriangle);
    return 1 ;
}






int PB043(PB_RESULT *pbR) {
    // attention 0 est le digit de poids faible, 9 celui de poids fort
    u_int8_t dig[10] = {0,1,2,3,4,5,6,7,8,9} ;
    int rg  ;
    int nbL = 0 ;
    u_int64_t Sum = 0 ;
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%03d ",pbR->pbNum);
    do  {
        if( ( (dig[0]+10*dig[1]+100*dig[2]) % 17 ) != 0 ) {
            // passer au suivant pour les 3 digits de poids faible
            rg = 2 ;
        } else if ( ( (dig[1]+10*dig[2]+100*dig[3]) % 13 ) != 0 ) {
            rg = 3 ;
        } else if ( ( (dig[2]+10*dig[3]+100*dig[4]) % 11 ) != 0 ) {
            rg = 4 ;
        } else if ( ( dig[4] % 5) != 0 ) {
            // test juste sur le dig4 car 10 et 100 divisible par 5
            rg = 4 ;
        } else if ( ( (dig[3]+10*dig[4]+100*dig[5]) % 7 ) != 0 ) {
            rg = 5 ;
        } else if ( ( dig[6] % 2 ) != 0 ) {
            // doit etre pair
            rg = 6 ;
        } else if ( ( (dig[5]+dig[6]+dig[7]) % 3 ) != 0 ) {
            rg = 7 ;
        } else {
            u_int64_t N = dig[0]+10*(dig[1]+10*(dig[2]+10*(dig[3]+10*(dig[4]+10*(dig[5]+10*(dig[6]+10*(dig[7]+10*(dig[8]+10LL*dig[9])))))))) ;
            if(pbR->isVerbose) fprintf(stdout,"%lld ",N);
            rg = 8 ;
            Sum += N ;
        }
        nbL++ ;
        
    } while(NextPermutRg(dig, 10, rg) >= 0) ;
    if(pbR->isVerbose) fprintf(stdout,"nbLoop=%d\n",nbL);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",Sum);
    return 1 ;
}

u_int32_t IsPentagonal ( u_int64_t P) {
    u_int32_t n = (u_int32_t) Sqrt64((2*P)/3) + 1 ;
    if(P*2 == ((u_int64_t)n )*(3*n -1)) {
        return n ;
    } else {
        return 0 ;
    }
}

typedef  int64_t PENT_T  ;
int PB044(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%03d ",pbR->pbNum);
    {
        PENT_T Pi,Pj ;
        int32_t  i,j,k,r ;
        for(k=2;;k++) {
            int32_t d ;
            int32_t Pk2 = k*(3*k-1) ;
            // on utilise l'idee que si Pk = Pi - Pj on a
            // k(3k-1) = (i-j) * (3i + 3j - 1)
            // donc on peut chercher les factorisation de k(3k-1) avec :
            // i-j = d ; 3i+3j-1 = Pk2/d <=> j = (Pk2/d-3*d+1)/6 et j > 0
             for(d=1; 0 < Pk2 - 3*d*d + d ;d++)  { // boucle sur les diviseurs
                if((Pk2 % d)==0 ) {
                    j = (Pk2/d-3*d+1)/6 ;
                    i = j + d ;
                    Pi = (i*(PENT_T)(3*i-1))/2 ;
                    Pj = (j*(PENT_T)(3*j-1))/2 ;
                    if((r=IsPentagonal(Pi+Pj)) && IsPentagonal(Pi-Pj) ) {
                        pbR->nbClock = clock() - pbR->nbClock ;
                        if(pbR->isVerbose)fprintf(stdout,"P(%d)=%d P(%d)=%lld P(%d)=%lld=%lld+%lld P(%d)=%lld=%lld+%lld\n",k,Pk2/2,j,Pj,i,Pi,Pi-Pj,Pj,r,Pi+Pj,Pi,Pj) ;
                        sprintf(pbR->strRes,"%d",Pk2/2);
                        return 1 ;

                    }
                }
            }
        }
    }
}

// Tn = n(n+1)/2
// Hn = n(2n-1) = 2n(2n-1)/2
// Pn = n(3n-1)/2
// les hexagonaux sont les triangles de rang impair
// donc on a va chercher un pentagonal qui est aussi un trinagle impair
#define PB045_NBS   2
int PB045(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%03d ",pbR->pbNum);
    {
        PENT_T Pk2 ;
        int32_t  k, nbSol=PB045_NBS;
        for(k=2;nbSol>0;k++) {
            Pk2 = k*(PENT_T)(3*k-1) ;
            int32_t d = (int32_t) Sqrt64(Pk2) ;
            if((d&1) && (d*(PENT_T)(d+1) == Pk2) ) {
                if(pbR->isVerbose)fprintf(stdout,"P(%d)=%lld=T(%d)=H(%d) ",k,Pk2/2,d,(d+1)/2) ;
                nbSol-- ;
                if(nbSol == PB045_NBS - 2) sprintf(pbR->strRes,"%lld",Pk2/2);
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\n");
    return 1 ;
}

#define PB046_MAXN  100000
int PB046(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    u_int8_t * oddDec = calloc(PB046_MAXN/2, sizeof(oddDec[0])) ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB046_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    {
        int i ;
        int curOdd = 1;
        for(i=1;i<ctxP->nbPrime;i++) {
            int n ;
            T_prime p = ctxP->tbPrime[i] ;
            int indP = p/2 ;
            while(curOdd < indP) {
                if(!oddDec[curOdd]) {
                    Free_tablePrime(ctxP) ;
                    sprintf(pbR->strRes,"%d",2*curOdd+1) ;
                    pbR->nbClock = clock() - pbR->nbClock ;
                    return 1 ;
                }
                curOdd++ ;
            }
            curOdd++ ;
            indP += 1 ; // (2* 1*1 /2)
            for(n=1;indP<PB046_MAXN/2;n++) {
                oddDec[indP] = 1 ;
                indP += 2*n+1 ;
            }
        }
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}

// donc test jusqu'a 10000*10000
#define PB047_MAXN  10000
int PB047(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB047_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    {
        int n=2 ;
        for(n=2;n<PB047_MAXN*PB047_MAXN;) {
            if(FindNbDivPrime(n,ctxP->tbPrime) < 4 ) {
                n++ ; continue ;
            }
            if(FindNbDivPrime(n+1,ctxP->tbPrime) < 4 ) {
                n += 2 ; continue ;
            }
            if(FindNbDivPrime(n+2,ctxP->tbPrime) < 4 ) {
                n += 3 ; continue ;
            }
            if(FindNbDivPrime(n+3,ctxP->tbPrime) < 4 ) {
                n += 4 ; continue ;
            }
            sprintf(pbR->strRes,"%d",n);
            Free_tablePrime(ctxP) ;
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}


#define PB048_MAXN  1000
#define PB048_MOD   10000000000LL // dix derniet digits
#define PB048_MODHALF   100000
u_int64_t mult10d(u_int64_t a, u_int64_t b) {
    u_int64_t al = a % PB048_MODHALF ;
    u_int64_t ah = ( a / PB048_MODHALF ) % PB048_MODHALF ;
    u_int64_t bl = b % PB048_MODHALF ;
    u_int64_t bh = ( b / PB048_MODHALF ) % PB048_MODHALF ;
    return (al * bl + (al * bh + ah * bl) * PB048_MODHALF) % PB048_MOD ;
}
int PB048(PB_RESULT *pbR) {
    u_int64_t Sum = 0 ;
    int n;
    pbR->nbClock = clock()  ;
    for(n=1;n<=PB048_MAXN;n++) {
        int exp ;
        u_int64_t pow = n ;
        u_int64_t Npow= 1 ;
        for(exp=n ;exp != 0; exp >>=1) {
            if(exp & 1) Npow = mult10d(Npow,pow)  ;
            pow = mult10d(pow,pow)  ;
        }
        Sum = (Sum + Npow) % PB048_MOD ;
    }
    sprintf(pbR->strRes,"%lld",Sum);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

typedef struct perm4d {
    u_int16_t nb ;
    u_int16_t  perm[24] ;
} perm4d;

#define PB049_MAXN  10000
#define PB049_SOL1  1487
int PB049(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    perm4d *Permut ;
    int i;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB047_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    Permut = calloc(PB049_MAXN,sizeof(Permut[0]));
    for(i=0;ctxP->tbPrime[i]<1000;i++) ;
    for(;i<ctxP->nbPrime;i++) {
        int num ;
        u_int8_t dig[4];
        int p = ctxP->tbPrime[i] ;
        dig[0] = p / 1000 ;
        dig[1] = (p/100) % 10 ;
        dig[2] = (p/10) % 10 ;
        dig[3] = p % 10 ;
        HeapSortUint8(dig,4);
        num = dig[0] * 1000 + dig[1] * 100 + dig[2]*10 + dig[3] ;
        Permut[num].perm[Permut[num].nb++] = p;
    }
    for(i=1000;i<PB049_MAXN;i++) {
        int nb = Permut[i].nb ;
        int i1,i2,i3 ;
        if(nb < 3 || Permut[i].perm[0] == PB049_SOL1) continue ;
        for(i1=0;i1<nb-2;i1++) {
            for(i2=i1+1;i2<nb-1; i2++) {
                int e3 = 2*Permut[i].perm[i2] - Permut[i].perm[i1] ;
                for(i3=i2+1;i3<nb; i3++) {
                    if( Permut[i].perm[i3] == e3) {
                        sprintf(pbR->strRes,"%d%d%d",Permut[i].perm[i1],Permut[i].perm[i2],Permut[i].perm[i3]) ;
                        i = PB049_MAXN ; i3= i2 = i1 = nb ;
                    }
                }
            }
        }
        
    }
//    sprintf(pbR->strRes,"%lld",Sum);
    free(Permut);
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB050_MAXN  1000000
int PB050(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    int32_t indP0,indP1 ;
    int32_t S ;
    u_int32_t bestIndP0 = 0 ;
    u_int32_t bestLg = 0 ;
    u_int32_t bestP = 0 ;
    pbR->nbClock = clock()  ;
    // we precompute N = sqrt(Max) * log2(Max)  primes
    u_int32_t sizeP = (u_int32_t)Sqrt64(PB050_MAXN) ;
    {   int i=1, nb = 1;
        while((i *= 2) < PB050_MAXN ) { nb++ ; }
        sizeP *= nb ;
    }
    if((ctxP = Gen_tablePrime(sizeP)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    indP0 = 0; //first prime of the sum is index=0 (2)
    indP1 = indP0 ; S=0;
    do { // loop on the first prime
        int32_t indP2 ;
        int32_t S1 ;
        // on cherche lindex indP1 S > AMX
        while(S<PB050_MAXN) {
            S += ctxP->tbPrime[indP1++] ;
        }
        // decrease the chain by the end until we reach a prime
        indP2 = indP1;
        S1 = 0 ;
        while(indP2>=indP0+bestLg) { //check not poor result
            S1 += ctxP->tbPrime[--indP2] ;
            if(Is_Prime(S-S1,ctxP->tbPrime)){
                bestLg =indP2 - indP0 ;
                bestIndP0 = indP0 ;
                bestP = S-S1 ;
                break ;
            }
        }
        // increment the index of the first prime in the sum
        S -= ctxP->tbPrime[indP0++] ;
    //exit if the chain beginning at indP0 and length=bestLg
    // has a SUM exeeding MAX. It will be worse for superior indP0
    } while (indP0+bestLg < indP1) ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%03d lg=%d %d=sum(%d..%d)\n",pbR->pbNum,bestLg,bestP,
                               ctxP->tbPrime[bestIndP0],ctxP->tbPrime[bestIndP0+bestLg-1]);

    pbR->nbClock = clock() - pbR->nbClock ;
    Free_tablePrime(ctxP);
    sprintf(pbR->strRes,"%d",bestP);
    return 1 ;
 }

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
        if(!Is_Prime(p,ctxP->tbPrime)) continue ;
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
                        if(Is_Prime(p+j*delta,ctxP->tbPrime)) tbP[nbP++] = p+j*delta ;
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

#define PB056_MAX   100
#define PB056_MAXL  201
int PB056_gmp(PB_RESULT *pbR) {
    char  digLarge[PB056_MAXL] ;
    u_int32_t Max = 0 ;
    int a,b ;
    int aBest =0, bBest = 0 ;
    pbR->nbClock = clock() ;
    mpz_t pow ;
    mpz_init(pow);
    for(b=PB056_MAX-1;b>1;b--) {
        if(18*b <= Max) break ; // car max 2*b digits
        for(a=PB056_MAX-1;a>1;a--) {
            int i , S = 0 ;
            mpz_ui_pow_ui (pow, a, b);
            gmp_sprintf (digLarge,"%Zd",pow);
            if(strlen(digLarge)*9 <= Max ) continue ;
            for(i=0;digLarge[i] != 0; i++) {
                S += digLarge[i] - '0' ;
            }
            if(S > Max) {
                Max = S ;
                aBest = a;
                bBest = b ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp DigSum(%d**%d)=%d\n"
                              ,pbR->pbNum
                              ,aBest,bBest,Max) ;
    sprintf(pbR->strRes,"%d",Max) ;
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
    for(N=3;N< PB058_MAXN;N += 2 ) {
        nbPoint += 4 ;
        DL += 4*N - 6 ;
        if(Is_Prime(DL,ctxP->tbPrime)) {
            nbPrime++ ;
        }
        if(Is_Prime(DL-N+1,ctxP->tbPrime)) {
            nbPrime++ ;
        }
        if(Is_Prime(DL-2*N+2,ctxP->tbPrime)) {
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
    
    NPMAX = ctxP->nbPrime-1 ;
    {   // on limite la valeur max des nb premiers a tester pour que la concat de 2 d'entre eux
        // soit testable avec la table calculee
        int32_t pow10 = 10 ;
        int32_t maxP = ctxP->tbPrime[ctxP->nbPrime-1] ;
        while(pow10 < PB060_MAXP) pow10 *= 10 ;
        u_int64_t maxP2 = maxP * (u_int64_t) maxP ;
        while(ctxP->tbPrime[NP] * (u_int64_t) (pow10 + 1) >= maxP2 ) NPMAX-- ;
    }
  
    NP = NPMAX ;
     maxS = ctxP->tbPrime[NP-1] ;
     if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d NPMAX=%d,maxS=%d,maxP=%d\n",pbR->pbNum,NP,maxS,ctxP->maxValue);
     pow10 = calloc(NP, sizeof(pow10[0])) ;
     {
         int i ;
         for(i=0;i<NP;i++) {
             int32_t Pi = ctxP->tbPrime[i] ;
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
        if(minS==0) maxS = ctxP->tbPrime[NP-1] ;
         
        for(i=0,index[i]=0,SP=0;i>=0;) {
            int j , isOK ;
            int ip = index[i] ;
            isOK = 1 ;
            // borne sur l'index et test que la somme P[k] < maxS avec P[k] croissant
            if(ip >=NP || SP+(5-i)*(P[i]=ctxP->tbPrime[ip]) > maxS ) {
                if(--i >= 0) { SP -= P[i] ; index[i]++ ; }
                continue ;
          } else {
                for(j=0;j<i;j++) {
                    int jp = index[j] ; // test compatibilite P[j] et le nouveau P[i]
                    // test non calcule (zero , sinon 1 pour non compat et 2 pour compat
                    if(isP1P2[NP*jp+ip]==0)  { isP1P2[NP*jp+ip] = Is_Prime2(P[i]+P[j]*pow10[ip],P[j]+P[i]*pow10[jp],ctxP->tbPrime) ? 2 : 1 ; }
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

typedef struct FractContG {
    mpz_t   n1 ;
    mpz_t   n0 ;
    mpz_t   d1 ;
    mpz_t   d0 ;
} FractContG ;

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


int PB066(PB_RESULT *pbR) {
    int32_t N , n , d ;
    int32_t k0, k02 ;
    pbR->nbClock = clock()  ;
    mpz_t max_x ;
    mpz_init_set_ui (max_x, 0) ;
    int bestN = 0 ;
    FractContG * FCG = FCG_alloc() ;
    for(N=2,k0=1,k02=4;N<=1000;N++) {
        int32_t i,a;
        if(N == k02) { // k02 = (k0+1)*(k0+1)
            k0++ ;
            k02 += 2*k0 + 1 ; continue ;
        }
        n = k0 ; d=1 ; i = 0 ; // so k0 =(int) srqt(N)
        FCG_init(FCG,k0) ;
        do {
            d = (N - n * n) / d ; // quite easy to recursively show that division is exact
            a = (k0+n)/d ;
            n = k0 - ( (k0 + n) % d ) ;
            i++ ;
            FCG_NextCoef(FCG,a) ;
       }while(d !=1 || n != k0 || (i&1)) ; // test loop on (n,d) = (k0,1) et i even (solution of PELL equation)
       if(FCG_CmpNum0(FCG,max_x) > 0 ) { // ask for n0 because one step ahead
            FCG_GetNum0(FCG,max_x) ;
            bestN = N ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) {
        char * str_x = mpz_get_str(NULL,10,max_x) ;
        fprintf(stdout,"\t PB%0.3d x(%d)=%s\n",pbR->pbNum, bestN,str_x);
        free (str_x) ;
    }
    sprintf(pbR->strRes,"%d",bestN);
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
    for(i=0;i<ctxP->nbPrime && (ctxP->tbPrime[i] < nSqrt) ;i++) ;
    // on continue la boucle sur les Pi (croissants) et jmax est la valeur maxi pour que Pi*Pj < N)
    for(jmax=i-1;i<ctxP->nbPrime;i++){
        int Pi = ctxP->tbPrime[i] ;
        int Pj ;
        int n,phi ;
        for(j=jmax;j>=0; j--) {
            if( (n=Pi*(Pj=ctxP->tbPrime[j])) >= N) { jmax-- ; continue ; }
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


static u_int16_t FactLg(TY_SYR val) {
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
    int32_t *Epn = calloc(ctxP->nbPrime*PB077_MAXN,sizeof(Epn[0])) ;
    int ip,k ;
    int p = 2 ;
    for(k=0;2*k<PB077_MAXN;k++) {
        Epn[2*k] = 1 ; Epn[2*k+1] = 0 ;
    }
    
    for(ip=1;ip<ctxP->nbPrime;ip++) {
        p = ctxP->tbPrime[ip] ;
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
        if(Epn[(ctxP->nbPrime-1)*PB077_MAXN + kmin]-1 > PB077_NB_ASK) {
            if(pbR->isVerbose) fprintf(stdout,"\t PB%0.3d %d a %d decompositions en premiers\n"
                                       ,pbR->pbNum,kmin,Epn[(ctxP->nbPrime-1)*PB077_MAXN + kmin]-1);
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

#define PB080_NBD   100
#define PB080_N     100

int PB080_gmp(PB_RESULT *pbR) {
    char  digDecimal[PB080_NBD+10] ;
   int S=0,i,k ,nxt_k2 ;
    pbR->nbClock = clock() ;
    mpz_t n2 ;
    mpz_init(n2);
    for(i=1,k=1,nxt_k2=1;i<PB080_N;i++) {
        if(i == nxt_k2) {
            nxt_k2 += 2*k + 1 ;
            k++ ;
            continue ;
        }
        mpz_ui_pow_ui (n2,10,(PB080_NBD-1)*2);
        mpz_mul_si (n2,n2,i) ;
        mpz_sqrt (n2,n2) ;
        gmp_sprintf (digDecimal,"%Zd",n2);
        {
            int j ;
            for(j=0;digDecimal[j] != 0; j++) {
                S += digDecimal[j] - '0' ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp SumdDecimal[1..%d]=%d\n"
                              ,pbR->pbNum,PB080_N,S) ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
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
    for(i=i2=i3=i4=0;i<ctxP->nbPrime;i++) {
        u_int32_t p,p2 ;
        p = ctxP->tbPrime[i] ;
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


#define PB597_LG    1800
#define PB597_NB    13
#define PB597_SIZEB 40
#define PB597_DEBUG 0

typedef struct devTorpids_a {
    mpq_t   P ;                   // proba
    u_int32_t   bitEnv ;
    u_int16_t   isEven ;
    
} devTorpids_a ;

#if PB597_DEBUG
static void PB597_printTorpids_a(FILE *fout,devTorpids_a * dT) {
    gmp_fprintf(fout,"%Qd %c:",dT->P,dT->isEven ? 'E' : 'O') ;
    int bitTest,ndec ;
    for(bitTest=1,ndec=0;bitTest<=dT->bitEnv;bitTest <<=1 , ndec++) { if(bitTest & dT->bitEnv) fprintf(fout,"%c%c", (bitTest==1) ? '{' : ',','A' + ndec) ; }
    fprintf(fout,"}\n");
}
#endif

#define MPQ_INIT_ND(mpq,N,D)   mpq_set_ui(mpq,N,D) ; mpq_canonicalize(mpq)


int PB597_a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    u_int32_t nb, curBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    // calcul factoriel n!
    for(i=1,nb=1;i<=PB597_NB;i++) nb *= i ;
    devTorpids_a * antDT = malloc(nb*sizeof(antDT[0])) ;
    devTorpids_a * curDT = malloc(nb*sizeof(curDT[0])) ;
    u_int32_t antNb, curNb ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    {
        u_int64_t pg = PGCD64(curLg,sizeB);
        curLg /= pg ;
        sizeB /= pg ;
    }
    curNb= 0 ;
    mpq_t gcd, ND,ONE ;
    mpq_init(gcd);
    mpq_init(ND);
    mpq_init(ONE) ;
    mpq_set_ui(ONE,1,1);

    mpq_init(curDT[curNb].P) ;
    mpq_set_ui(curDT[curNb].P,1,1) ;
    curDT[curNb].isEven = 1;
    curDT[curNb].bitEnv = 1 ;
    sumLg[0] = Lg[0] = curLg ;
    curNb++ ;
     for(curBoat=1;curBoat<PB597_NB;curBoat++) {
        devTorpids_a *tmp = antDT ;
        antDT = curDT ;
        curDT = tmp ;
        antNb = curNb ;
        curNb = 0 ;
        curLg += sizeB ;
        Lg[curBoat] = curLg ;
        sumLg[curBoat] = sumLg[curBoat-1] + curLg ;
#if PB597_DEBUG
      fprintf(stdout,"==%d\n",curBoat);
#endif
        for(i=0;i<antNb;i++) {
            devTorpids_a * ptAntDT = antDT + i ;
            MPQ_INIT_ND(ND,curLg,sumLg[curBoat]) ;
            // on rajoute le bateau qui ne bumpe personne
            mpq_init (curDT[curNb].P) ;
            mpq_mul (curDT[curNb].P, ptAntDT->P, ND) ;
            curDT[curNb].bitEnv = 1 << curBoat ;
            curDT[curNb].isEven = ptAntDT->isEven ;
#if PB597_DEBUG
            printf("*");
            PB597_printTorpids_a(stdout,curDT+curNb);
#endif
            curNb++ ;
            int place ;
            u_int16_t nboat ;
            // on va maintenant rajouter les bumps
            for(place=0,nboat=0;(1<< nboat) <= ptAntDT->bitEnv ;place++) {
                while(((1<< nboat) & ptAntDT->bitEnv) == 0) {   nboat++ ; }
                mpq_init(curDT[curNb].P) ;
               if(place == 0) {
                   MPQ_INIT_ND(ND,sumLg[curBoat] - curLg,sumLg[curBoat]) ;
                   mpq_mul(curDT[curNb].P,ptAntDT->P,ND);
               } else {
                   mpq_inv(ND,ND) ; // ND = D/N
                   mpq_sub(ND,ND,ONE); // ND = D/N -1 = (D-N)/N
                   mpq_mul(curDT[curNb].P,curDT[curNb-1].P,ND);
               }
               MPQ_INIT_ND(ND,(curBoat - nboat) * sizeB,sumLg[curBoat]-sumLg[nboat] - (curBoat - nboat) * Lg[nboat]) ;
               mpq_mul(curDT[curNb].P,curDT[curNb].P,ND);
               nboat++ ;
               curDT[curNb].bitEnv =  (ptAntDT->bitEnv & ((1 << nboat) -1)) |  (1 << curBoat)  ;
               curDT[curNb].isEven = (place & 1) ? ptAntDT->isEven : !ptAntDT->isEven ;
#if PB597_DEBUG
               PB597_printTorpids_a(stdout,curDT+curNb);
#endif
               curNb++ ;
           }
        }
    }
    {
        mpq_t ND ;
        mpq_init(ND);
        mpq_set_ui(ND,0,1);
        for(i=0;i<curNb;i++) {
            if(curDT[i].isEven) {
                mpq_add(ND,ND,curDT[i].P) ;
            }
        }
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%da PEVEN=%Qd\n",pbR->pbNum,ND );
        mpz_mul_ui(mpq_numref(ND) ,mpq_numref(ND),100000000000) ;
        mpz_tdiv_q (mpq_numref(ND),mpq_numref(ND),mpq_denref(ND));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(ND) );
    }
    free(curDT); free(antDT) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


typedef struct devTorpids {
    u_int16_t   isUsed ;
    mpq_t   E_P ;                   // proba for even permutations
    mpq_t   O_P ;                   // proba for odd permutations
} devTorpids ;

#if PB597_DEBUG
static void PB597_printTorpids(FILE *fout,devTorpids * dT,u_int32_t index) {
    gmp_fprintf(fout,"[%d]E:%Qd O:%Qd",index,dT->E_P,dT->O_P) ;
    int bitTest,kboat ;
    for(bitTest=1,kboat=0;bitTest<=index;bitTest <<=1 , kboat++) { if(bitTest & index) fprintf(fout,"%c%c", (bitTest==1) ? '{' : ',','A' + kboat) ; }
    
    fprintf(fout,"}\n");
}
#endif
#define MPQ_ADD_EO(dst,src1,src2)   mpq_add((dst)->E_P,(src1)->E_P,(src2)->E_P) ; mpq_add((dst)->O_P,(src1)->O_P,(src2)->O_P)
#define MPQ_SET_ND_EO(dst,ND)       mpq_set((dst)->E_P,ND) ; mpq_set((dst)->O_P,ND)
#define MPQ_INIT_SET_EO(dst,src)       mpq_init((dst)->E_P); mpq_set((dst)->E_P,(src)->E_P) ;mpq_init((dst)->O_P); mpq_set((dst)->O_P,(src)->O_P)
#define MPQ_MUL_ND_EO(dst,src1,ND)   mpq_mul((dst)->E_P,(src1)->E_P,(ND)) ; mpq_mul((dst)->O_P,(src1)->O_P,(ND))
// crossing odd and even
#define MPQ_MUL_ND_ExO(dst,src1,ND)   mpq_mul((dst)->E_P,(src1)->O_P,(ND)) ; mpq_mul((dst)->O_P,(src1)->E_P,(ND))
#define MPQ_MUL_EO(dst,src1,src2)   mpq_mul((dst)->E_P,(src1)->E_P,(src2)->E_P) ; mpq_mul((dst)->O_P,(src1)->O_P,(src2)->O_P)
#define MPQ_SELFMUL_ND_ExO(dst,ND,tmp)   mpq_mul(tmp,(dst)->O_P,(ND)) ; mpq_mul((dst)->O_P,(dst)->E_P,(ND)) ; mpq_set((dst)->E_P,tmp)

#define MPQ_INIT_MUL_ND_EO(dst,src1,ND)   mpq_init((dst)->E_P); mpq_mul((dst)->E_P,(src1)->E_P,(ND)) ;mpq_init((dst)->O_P); mpq_mul((dst)->O_P,(src1)->O_P,(ND))
#define MPQ_ADD_MUL_ND_EO(dst,src1,ND,tmp)    mpq_mul(tmp,(src1)->E_P,(ND)) ; mpq_add((dst)->E_P,(dst)->E_P,tmp) ; mpq_mul(tmp,(src1)->O_P,(ND)) ; mpq_add((dst)->O_P,(dst)->O_P,tmp)


int PB597(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    u_int32_t  curBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    mpq_t tmpz , ND ,ONE ;
    devTorpids * DT = calloc(1 << (PB597_NB+1),sizeof(DT[0])) ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    devTorpids newDT ;
    mpq_init(ND);    mpq_init(tmpz);
    mpq_init(ONE) ;  mpq_set_ui(ONE,1,1);
    mpq_init(newDT.E_P) ; mpq_init(newDT.O_P) ;
    {
        u_int64_t pg = PGCD64(curLg,sizeB);
        curLg /= pg ;
        sizeB /= pg ;
    }
    {
        u_int32_t indNb= 1 ; // 1 <<0
        DT[indNb].isUsed = 1 ;
        mpq_init(DT[indNb].E_P);   mpq_set_ui(DT[indNb].E_P,1,1);
        mpq_init(DT[indNb].O_P);   mpq_set_ui(DT[indNb].O_P,0,1);
    }
    sumLg[0] = Lg[0] = curLg ;
    for(curBoat=1;curBoat<PB597_NB;curBoat++) {
        u_int32_t indAnt ;
        curLg += sizeB ;
        Lg[curBoat] = curLg ;
        sumLg[curBoat] = sumLg[curBoat-1] + curLg ;
#if PB597_DEBUG
        fprintf(stdout,"==%d\n",curBoat);
#endif
        for(indAnt=1<<(curBoat-1);indAnt< (1<<curBoat);indAnt++) {
            devTorpids * ptAntDT = DT + indAnt ;
            if(ptAntDT->isUsed == 0) continue ;
            MPQ_INIT_ND(ND,curLg,sumLg[curBoat]) ;
            devTorpids * ptCurDT = DT + (1 << curBoat) ;
            // on rajoute le bateau qui ne bumpe personne
            if(ptCurDT->isUsed) {
                MPQ_ADD_MUL_ND_EO(ptCurDT,ptAntDT,ND,tmpz) ; // ptCurDT += ptAntDT * N/D
            } else {
                ptCurDT->isUsed = 1 ;
                MPQ_INIT_MUL_ND_EO(ptCurDT,ptAntDT,ND) ; // ptCurDT = ptAntDT * N/D
            }
#if PB597_DEBUG
            printf("*");
            PB597_printTorpids(stdout,ptCurDT,1 << curBoat);
#endif
            int isFirst ;
            u_int16_t nboat ;
            //
            for(isFirst=1,nboat=0;(1<< nboat) <= indAnt ;) {
                while(((1<< nboat) & indAnt) == 0) {   nboat++ ; }
                if(isFirst) { // Begin ?
                    MPQ_INIT_ND(ND,sumLg[curBoat] - curLg,sumLg[curBoat]) ;
                    MPQ_MUL_ND_ExO(&newDT,ptAntDT,ND) ; // newDT = ptAndtDT * N/D with crossing Even/Odd
                    isFirst = 0 ;
                } else { // on croise even et odd
                    mpq_inv(ND,ND) ; // ND = D/N
                    mpq_sub(ND,ND,ONE); // ND = D/N -1 = (D-N)/N
                    MPQ_SELFMUL_ND_ExO(&newDT,ND,tmpz) ; // newtDT *= N/D with crossing Even/Odd
                }
                MPQ_INIT_ND(ND,(curBoat - nboat) * sizeB,sumLg[curBoat]-sumLg[nboat] - (curBoat - nboat) * Lg[nboat]) ;
                MPQ_MUL_ND_EO(&newDT,&newDT,ND) ; // newDT *= N/D
                nboat++ ;
                ptCurDT = DT + ( (indAnt & ((1 << nboat) -1)) |  (1 << curBoat) ) ;
                if(ptCurDT->isUsed) {
                    MPQ_ADD_EO(ptCurDT,ptCurDT,&newDT);  // ptCurDT += newDT
                } else {
                    ptCurDT->isUsed = 1 ;
                    MPQ_INIT_SET_EO(ptCurDT,&newDT) ; // ptCurDT = newDT
                }
#if PB597_DEBUG
                PB597_printTorpids(stdout,ptCurDT,  ( (indAnt & ((1 << nboat) -1)) |  (1 << curBoat) ) );
#endif
            }
        }
    }
    {
        mpq_t ND ;
        mpq_init(ND);
        mpq_set_ui(ND,0,1);
        for(i= (1 << (PB597_NB-1));i< ( 1<< (PB597_NB)) ;i++) {
            if(DT[i].isUsed) {
                mpq_add(ND,ND,DT[i].E_P) ;
           }
        }
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%d PEVEN=%Qd\n",pbR->pbNum,ND );
        mpz_mul_ui(mpq_numref(ND) ,mpq_numref(ND),100000000000) ;
        mpz_tdiv_q (mpq_numref(ND),mpq_numref(ND),mpq_denref(ND));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(ND) );
    }
    free(DT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB597_x(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int32_t  nbBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    mpq_t tmpz , ND ,ONE, EVEN , FACT ;
    devTorpids * DT = calloc(1 << (PB597_NB),sizeof(DT[0])) ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    devTorpids  tmpDT ;
    mpq_init(ND);      mpq_init(tmpz); mpq_init(FACT) ;
    mpq_init(ONE) ;    mpq_set_ui(ONE,1,1);
    mpq_init(EVEN);    mpq_set_ui(EVEN,0,1) ;
    mpq_init(tmpDT.E_P) ; mpq_init(tmpDT.O_P) ;
    u_int64_t pg = PGCD64(curLg,sizeB);
    curLg /= pg ;  sizeB /= pg ; // to avoid big numbers
    DT[1].isUsed = 1 ; // init for race with one boat
    mpq_init(DT[1].E_P);   mpq_set_ui(DT[1].E_P,1,1);
    mpq_init(DT[1].O_P);   mpq_set_ui(DT[1].O_P,0,1);
    sumLg[0] = Lg[0] = curLg ;
    for(nbBoat=1;nbBoat<PB597_NB;nbBoat++) { // add boat by boat
        u_int32_t indAnt ;
        curLg += sizeB ;
        Lg[nbBoat] = curLg ; // length to run for new boat
        sumLg[nbBoat] = sumLg[nbBoat-1] + curLg ; // sigma length (usefull for coefficients)
#if PB597_DEBUG
        if(nbBoat != PB597_NB-1) fprintf(stdout,"==%d\n",nbBoat);
#endif
        // loop on states for race with curBoat - 1
        for(indAnt=1<<(nbBoat-1);indAnt< (1<<nbBoat);indAnt++) {
            devTorpids * ptAntDT = DT + indAnt ;
            if(ptAntDT->isUsed == 0) continue ;
            MPQ_INIT_ND(ND,curLg,sumLg[nbBoat]) ;
            devTorpids * ptCurDT = DT + (1 << nbBoat) ;
            if(nbBoat == PB597_NB-1) {// last loop cumulate only even
                mpq_mul(tmpz,ptAntDT->E_P,ND) ;
                mpq_add(EVEN,EVEN,tmpz) ;
            } else {
                if(ptCurDT->isUsed) {
                    MPQ_ADD_MUL_ND_EO(ptCurDT,ptAntDT,ND,tmpz) ; // ptCurDT += ptAntDT * N/D
                } else {
                    ptCurDT->isUsed = 1 ;
                    MPQ_INIT_MUL_ND_EO(ptCurDT,ptAntDT,ND) ; // ptCurDT = ptAntDT * N/D
                }
            }
#if PB597_DEBUG
            if(nbBoat != PB597_NB-1) {
                printf("*");
                PB597_printTorpids(stdout,ptCurDT,1 << nbBoat);
            }
#endif
            int nplace ;
            int indDT ;
            u_int64_t  nFACT, dFACT, pg ;
            nFACT = sumLg[nbBoat] - curLg ;
            dFACT = sumLg[nbBoat] ;
            u_int16_t kboat ;
            // loop on kboat for the precedent state
            // newboat takes his place
            for(nplace=0,kboat=0; ;nplace++) {
                while(((1<< kboat) & indAnt) == 0) {   kboat++ ; }
                u_int32_t N=(nbBoat - kboat) * sizeB ;
                u_int32_t D= sumLg[nbBoat]-sumLg[kboat] - (nbBoat - kboat) * Lg[kboat] ;
                kboat++ ;
                indDT =  (indAnt & ((1 << kboat) -1)) |  (1 << nbBoat) ;
                nFACT *= N ;
                dFACT *= D ;
                pg = PGCD64(nFACT,dFACT) ;
                if(pg > 1) { nFACT /= pg ; dFACT /= pg ; }
                ptCurDT = DT + indDT  ;
                MPQ_INIT_ND(FACT,nFACT,dFACT);
                if(nbBoat == PB597_NB-1) {// last loop cumulate only even
                    if(nplace & 1) {
                        mpq_mul(tmpz,ptAntDT->E_P,FACT);
                    }else {
                        mpq_mul(tmpz,ptAntDT->O_P,FACT);
                    }
                    mpq_add(EVEN,EVEN,tmpz) ;
                } else {
                    if(nplace & 1) {
                        MPQ_MUL_ND_EO(&tmpDT,ptAntDT,FACT) ;
                    } else {
                        MPQ_MUL_ND_ExO(&tmpDT,ptAntDT,FACT) ;
                    }
                    if(ptCurDT->isUsed) {
                        MPQ_ADD_EO(ptCurDT,ptCurDT,&tmpDT);
                    } else {
                        ptCurDT->isUsed = 1 ;
                        MPQ_INIT_SET_EO(ptCurDT,&tmpDT) ;
                    }
                }
#if PB597_DEBUG
                if(nbBoat != PB597_NB-1) PB597_printTorpids(stdout,ptCurDT,   indDT  );
#endif
                if((1<< kboat) <= indAnt) {
                    nFACT *= (D-N) ;
                    dFACT *= N ;
                } else {
                    break ;
                }
            }
        }
    }
    {
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%d %d boats PEVEN=%Qd\n",pbR->pbNum,PB597_NB,EVEN );
        mpz_mul_ui(mpq_numref(EVEN) ,mpq_numref(EVEN),100000000000) ;
        mpz_tdiv_q (mpq_numref(EVEN),mpq_numref(EVEN),mpq_denref(EVEN));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(EVEN) );
    }
    free(DT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
typedef struct CTX_597 {
    mpq_t   * backEven ;    // for even(k,k) k=0..NB
    mpq_t   * frontEven ;   // for even(n,L-NB+n) n=0..NB
} CTX_597 ;

mpq_t * getEven(CTX_597 * ctx,int n,int lg) {
    return (n==lg) ? ctx->backEven : ctx->frontEven ;
}
void computePeven(CTX_597 * ctx,int n,int lg) {
    static int isInit=0 ;
    mpq_t   tmpS, tmpP ,tmpQ, ONE ;
    if(!isInit) { // first init
        mpq_init(ONE) ; mpq_set_ui(ONE,1,1); mpq_init(tmpS) ;
        mpq_init(tmpP); mpq_init(tmpQ) ; isInit = 1 ;
    }
    mpq_t  * even = getEven(ctx,n,lg) ;
    if(n<2) {
        mpq_init(even[n]) ; // assume first time
        mpq_set_ui(even[n],1,1);
        return ;
    } else {
        int k ;
        mpq_set_ui(tmpS,0,1) ; // tmpS = 0
        for(k=0;k<n;k++) {
            mpq_set(tmpP,even[n-k-1]) ;// tmpP = (Pa=EV(n-k-1,l-k-1))
            mpq_add(tmpP,tmpP,ctx->backEven[k]) ; // tmpP = Pa + (Pb=EV(k,k))
            mpq_mul(tmpQ,even[n-k-1],ctx->backEven[k]) ; // tmpQ = Pa * Pb
            mpq_sub(tmpP,tmpP,tmpQ) ;// tmpP = Pa + Pb - Pa * Pb
            mpq_sub(tmpP,tmpP,tmpQ); // tmpP = Pa + Pb - 2 * Pa * Pb
            if(((k) & 1) == 0) mpq_sub(tmpP,ONE,tmpP); // ? 1 - Pa - Pb + 2 * Pa * Pb
            mpq_set_ui(tmpQ,lg-k,1) ; // tmpQ = lg-k
            mpq_mul(tmpP,tmpP,tmpQ) ;  // tmpP *= lg-k
            mpq_add(tmpS,tmpS,tmpP) ; // tmpS += tmpP
        }
        mpq_set_ui(tmpQ,1,(n*(2*lg-n+1))/2); // n*(2*lg-n+1)
        mpq_init(even[n]) ; // assume first time
        mpq_mul(even[n],tmpQ,tmpS);
        return ;
    }
}

int PB597_y(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int32_t n, L, sizeB ;
    CTX_597 ctx ;
    L = PB597_LG ;
    sizeB = PB597_SIZEB ;
    u_int64_t pg = PGCD64(L,sizeB);
    L /= pg ;  sizeB /= pg ; // GCD to avoid big numbers
    // Back : NB+1 for even(k,k) k=0..NB ; Front : NB+1 for even(n,L-NB+n) n=0..NB
    ctx.backEven = malloc((PB597_NB+1)*sizeof(ctx.backEven[0]));
    ctx.frontEven = malloc((PB597_NB+1)*sizeof(ctx.frontEven[0]));
    // Loop on number of boats
    for(n=0 ;n<=PB597_NB;n++) {
        computePeven(&ctx,n,n) ;
        computePeven(&ctx,n,L- PB597_NB + n ) ;
    }
    mpz_t Num1000 ;    mpz_init(Num1000);
    mpz_mul_ui(Num1000 ,mpq_numref(getEven(&ctx,PB597_NB,L)[PB597_NB]),100000000000) ;
    mpz_tdiv_q (Num1000,Num1000,mpq_denref(getEven(&ctx,PB597_NB,L)[PB597_NB]));
    if(pbR->isVerbose){
        gmp_fprintf (stdout,"\t PB%dy %2d boats PEVEN=0.%Zd=%Qd\n",pbR->pbNum,PB597_NB,Num1000,getEven(&ctx,PB597_NB,L)[PB597_NB]);
    }
    gmp_sprintf (pbR->strRes,"%Zd",Num1000);
    pbR->nbClock = clock() - pbR->nbClock ;
    free(ctx.backEven) ; free(ctx.frontEven) ;
    return 1 ;
}

typedef struct P_EXP {
    u_int32_t   p;
    u_int32_t   exp ;
    
} P_EXP ;
#define PB622_NBP   11 // number of prime factors for 2**60-1
int PB622(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i, jp,k;
    u_int64_t S=0;
    P_EXP fact2exp60_m1[PB622_NBP] = {{3, 2}, {5, 2}, {7, 1}, {11, 1}, {13, 1}, {31, 1}, {41, 1}, {61,1}, {151, 1}, {331, 1}, {1321, 1}};
    P_EXP fact2div60byP[3][PB622_NBP] ;
    u_int32_t quotPrime60[3] = { (1<<(60/2))-1 , (1<<(60/3)) -1 , (1<<(60/5)) -1 } ;
    {
        for(i=0;i<PB622_NBP;i++) {
            u_int32_t p = fact2exp60_m1[i].p ;
            for(jp=0;jp<3;jp++) {
                int exp, N ;
                for(exp=0, N=quotPrime60[jp];(N % p )== 0; exp++ , N /= p) ;
                fact2div60byP[jp][i].p = p; fact2div60byP[jp][i].exp = exp ;
            }
        }
    }
    {  // loop to test all factor of 2**60-1
        int32_t exp[PB622_NBP] , ip=PB622_NBP-1 ;
        for(i=0;i<PB622_NBP;i++) exp[i] = 0 ;
        while(ip >= 0) {
            if(exp[ip]++ < fact2exp60_m1[ip].exp) {
                while(ip<PB622_NBP-1) exp[++ip] = 0 ;
                for(jp=0;jp<3;jp++) { // test not divisor of element quotPrime60
                    for(k=0;k<PB622_NBP;k++) {
                        if(fact2div60byP[jp][k].exp < exp[k] ) break ; // not divisor
                    }
                    if(k==PB622_NBP) break ; // divisor of quotPrime60[jp]
                }
                if(jp==3) { // 60 is the minimum value, good candidate
                    u_int64_t N ;
                    for(k=0,N=1;k<PB622_NBP;k++) {
                        int ie ;
                        for(ie=exp[k];ie>0;ie--) N *= fact2exp60_m1[k].p ;
                    }
                    S += N+1 ;
                }
            } else {  ip-- ; }
        }
    }
    sprintf (pbR->strRes,"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB1000_NUM  10000000
int PB1000(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrimeNb(PB1000_NUM)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    fprintf(stdout,"\t PB%0.4d prime[%d] = %d\n",pbR->pbNum, PB1000_NUM,ctxP->tbPrime[PB1000_NUM-1]) ;
    sprintf(pbR->strRes,"%d",ctxP->tbPrime[PB1000_NUM-1]);
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



typedef struct TotalRun {
    clock_t TotalClock ;
    int isVerbose ;
    int nbPBOK ;
    int nbPBerror ;
} TotalRun ;


void Execute(TotalRun *ttr, PB_CALL *pbCall) {
    PB_RESULT pbR ;
    memset(&pbR,0,sizeof(pbR)) ;
    pbR.isVerbose = ttr->isVerbose ;
    if( pbCall->pbSolve != NULL) {
        pbR.pbNum = pbCall->pbNum ;
        pbR.nbClock = 0 ;
        if(pbCall->pbSolve(&pbR)) {
            ttr->TotalClock += pbR.nbClock ;
            if(strcmp(pbCall->Solution,pbR.strRes)==0) {
                ttr->nbPBOK++ ;
                fprintf(stdout,"OK\tPB%03d(%.06fs) Sol=%s '%s'\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,(pbCall->name) ? pbCall->name : "") ;
            } else {
                ttr->nbPBerror++ ;
                fprintf(stdout,"ERROR\tPB%03d(%.06fs) Find=%s != Exp=%s\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,pbCall->Solution) ;
            }
        } else {
            ttr->TotalClock += pbR.nbClock ;
            ttr->nbPBerror++ ;
            fprintf(stdout," BAD EXECUTION PB%03d(%.06fs)\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC) ;
        }
        
    } else {
        fprintf(stdout,"SKIPPED\tPB%3d\n",pbCall->pbNum) ;
    }

    
}

static PB_CALL ALL_calls[] = {
    {  1,  PB001,  "233168"}
    ,{  2,  PB002,  "4613732"}
    ,{  3,  PB003,  "6857"}
    ,{  4,  PB004,  "906609"}
    ,{  5,  PB005,  "232792560"}
    ,{  6,  PB006,  "25164150"}
    ,{  7,  PB007,  "104743"}
    ,{  8,  PB008,  "23514624000"}
    ,{  9,  PB009,  "31875000"}
    ,{ 10,  PB010,  "142913828922","Summation of primes"}
    ,{ 11,  PB011,  "70600674"}
    ,{ 12,  PB012,  "76576500"}
    ,{ 13,  PB013,  "5537376230"}
    ,{ 14,  PB014,  "837799"}
    ,{ 15,  PB015,  "137846528820"}
    ,{ 16,  PB016,  "1366"}
    ,{ 16,  PB016_gmp,  "1366"}
    //        ,{ 17,  PB017,  ""}
    ,{ 18,  PB018,  "1074"}
    //        ,{ 19,  PB019,  ""}
    ,{ 20,  PB020,  "648"}
    ,{ 21,  PB021,  "31626"}
    //        ,{ 22,  PB022,  ""}
    ,{ 23,  PB023,  "4179871"}
    ,{ 24,  PB024,  "2783915460"}
    ,{ 25,  PB025,  "4782"}
    ,{ 26,  PB026,  "983"}
    ,{ 27,  PB027,  "-59231"}
    ,{ 28,  PB028,  "669171001"}
    ,{ 29,  PB029,  "9183"}
    ,{ 30,  PB030,  "443839"}
    ,{ 31,  PB031,  "73682"}
    ,{ 32,  PB032,  "45228"}
    ,{ 33,  PB033,  "100"}
    ,{ 34,  PB034,  "40730"}
    ,{ 35,  PB035,  "55"}
    ,{ 36,  PB036,  "872187"}
    ,{ 37,  PB037,  "748317"}
    ,{ 37,  PB037r, "748317"}
    ,{ 38,  PB038,  "932718654"}
    ,{ 39,  PB039,  "840","Integer right triangles"}
    ,{ 40,  PB040,  "210"}
    ,{ 41,  PB041,  "7652413"}
    ,{ 42,  PB042,  "162"}
    ,{ 43,  PB043,  "16695334890" , "Sub-string divisibility"}
    ,{ 44,  PB044,  "5482660" ,"Pentagon numbers"}
    ,{ 45,  PB045,  "1533776805" ,"Triangular, pentagonal, and hexagonal"}
    ,{ 46,  PB046,  "5777" ,    "Goldbach's other conjecture" }
    ,{ 47,  PB047,  "134043" ,  "Distinct primes factors" }
    ,{ 48,  PB048,  "9110846700" ,  "Self powers" }
    ,{ 49,  PB049,  "296962999629" ,"Prime permutations" }
    ,{ 50,  PB050,  "997651" ,  "Consecutive prime sum" }
    ,{ 51,  PB051,  "121313" ,   "Prime digit replacements" }
    ,{ 52,  PB052,  "142857", "Permuted multiples"}
    ,{ 52,  PB052a,  "142857", "Permuted multiples"}
    ,{ 53,  PB053,  "4075", "Combinatoric selections" }
    ,{ 55,  PB055,  "249" , "Lychrel numbers" }
    ,{ 56,  PB056_gmp,  "972",  "Powerful digit sum" }
    ,{ 57,  PB057,  "153",  "Square root convergents" }
    ,{ 58,  PB058,  "26241" ,   "Spiral primes" }
    ,{ 59,  PB059,  "107359" ,  "XOR decryption" }
    ,{ 60,  PB060,  "26033" ,  "Prime pair sets" }
    ,{ 61,  PB061,  "28684" , "Cyclical figurate numbers"}
    ,{ 62,  PB062,  "127035954683" , "Cubic permutations"}
    ,{ 63,  PB063,  "49" , "Powerful digit counts"}
    ,{ 64,  PB064,  "1322" , "Odd period square roots"}
    ,{ 65,  PB065,  "272" , "Convergents of e"}
    ,{ 66,  PB066,  "661", "Diophantine equation"}
    ,{ 67,  PB067,  "7273", "Maximum path sum II"}
    ,{ 68,  PB068,  "6531031914842725", "Magic 5-gon ring"}
//    ,{ 69,  PB069a,  "510510", "Totient maximum"},
    ,{ 69,  PB069b,  "510510", "Totient maximum"}
//    ,{ 70,  PB070,  "8319823", "Totient permutation"}
    ,{ 70,  PB070a,  "8319823", "Totient permutation"}
    ,{ 71,  PB071,  "428570", "Ordered fractions"}
    ,{ 72,  PB072,  "303963552391", "Counting fractions"}
    ,{ 73,  PB073,  "7295372", "Counting fractions in a range"}
//    ,{ 74,  PB074,  "402", "Digit factorial chains"}
    ,{ 74,  PB074a,  "402", "Digit factorial chains"}
    ,{ 75,  PB075,  "161667", "Singular integer right triangles"}
//    ,{ 76,  PB076,  "190569291", "Counting summations"}
    ,{ 76,  PB076a,  "190569291", "Counting summations"}
    ,{ 77,  PB077,  "71", "Prime summations"}
    ,{ 78,  PB078,  "55374", "Counting summations"}
    ,{ 80,  PB080_gmp,  "40886", "Square root digital expansion"}
    ,{ 81,  PB081,  "427337", "Path sum: two ways"}
    ,{ 82,  PB082,  "260324", "Path sum: three ways"}
    ,{ 83,  PB083,  "425185", "Path sum: four ways"}
    ,{ 85,  PB085,  "2772", "Counting rectangles"}
    ,{ 86,  PB086,  "1818", "Cuboid route"}
    ,{ 87,  PB087,  "1097343", "Prime power triples"}
    ,{ 88,  PB088,  "7587457", "Product-sum numbers"}
    ,{ 90,  PB090,  "1217", "Cube digit pairs"}
    ,{ 91,  PB091,  "14234", "Right triangles with integer coordinates"}

//    ,{100,  PB100,  "756872327473"}
//    ,{ 1000,  PB1000,  "179424673", "Test 10 000 000 prime numbers"}
// version moins rapide pour grande valeurs
//    ,{ 597,  PB597_a,  "50018178282", "Torpids"}
// version rapide mais calculs complets (pour debug)
//    ,{ 597,  PB597,  "50018178282", "Torpids"}
// version la plus rapide avec calcul FACT en u_int64_t et optimisation last loop
//    ,{ 597,  PB597_x,  "50018178282", "Torpids"}
// version beaucoup plus rapide en O(n**2)
    ,{ 597,  PB597_y,  "50018178282", "Torpids"}
    ,{ 622,  PB622,  "3010983666182123972", "Riffle Shuffles"}
    
    ,{  0,  NULL,   ""}
};


static PB_CALL CUR_calls[] = {
 //    { 51,  PB051,  "121313" ,   "Prime digit replacements" },
    { 622,  PB622,  "3010983666182123972", "Riffle Shuffles"},


    {  0,  NULL,   ""}
} ;


int main(int argc, const char * argv[]) {
    PB_CALL *ptCall ;
    TotalRun ttr ;
    

    ttr.isVerbose = 1 ;
    ttr.nbPBerror  = ttr.nbPBOK = 0 ;
    ttr.TotalClock = 0 ;
    int isCur = 1 ;
    for(ptCall = (isCur) ? CUR_calls : ALL_calls ; ptCall->pbNum != 0 ; ptCall++) {
        Execute(&ttr,ptCall);
    }
    fprintf(stdout,"### Execution of %d PB en %.06fs (%d en erreur)\n",ttr.nbPBOK+ttr.nbPBerror,(float) ttr.TotalClock / CLOCKS_PER_SEC,ttr.nbPBerror) ;
    return 0;
}
