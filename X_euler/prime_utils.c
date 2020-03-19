#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef  uint32_t T_prime ;
typedef struct CTX_PRIMETABLE CTX_PRIMETABLE;
CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) ;
const T_prime * GetTbPrime(CTX_PRIMETABLE * ctx) ;
uint32_t GetNbPrime(CTX_PRIMETABLE * ctx) ;

typedef int(*TY_CPL_nxtPrime)(void *ctx,T_prime nxtPrime);
uint32_t FindPrime(T_prime maxValue,void *ctx,TY_CPL_nxtPrime nxtPrime) ;

struct  CTX_PRIMETABLE {
    uint32_t   nbPrime ;
    T_prime     maxValue ;
    uint32_t   maxNbPrime ;
    T_prime     *tbPrime ;
}  ;

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

const T_prime * GetTbPrime(CTX_PRIMETABLE * ctx) {
    return ctx->tbPrime ;
}

uint32_t GetNbPrime(CTX_PRIMETABLE * ctx) {
    return ctx->nbPrime ;
}

CTX_PRIMETABLE * Free_tablePrime(CTX_PRIMETABLE * ctx) {
    if(ctx != NULL) free(ctx->tbPrime) ;
    free(ctx);
    return NULL ;
}

#define USE_PRIMEb  1



#define IS_Composed(p)  (isComposed[(p)/8] & (1 << ((p) & 0x7)) )
#define SET_Composed(p)  (isComposed[(p)/8] |=  (1 << ((p)  & 0x7)) )


uint32_t FindPrime_a(T_prime nbMax,void *ctx,TY_CPL_nxtPrime nxtPrime) {
    int32_t nSqrt = 1+ (int32_t)sqrt( (double) nbMax ) ;
    T_prime sizeTable = nbMax ;
    T_prime sizeTable2 = nbMax >> 1 ;
    uint8_t *isComposed = calloc( (sizeTable+15) /16,  sizeof(isComposed[0])) ;
    uint32_t nbPrime = 0 ;
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
uint32_t FindPrime_b(T_prime nbMax,void *ctx,TY_CPL_nxtPrime nxtPrime) {
    int32_t nSqrt = 1+ (int32_t)sqrt( (double) nbMax ) ;
    int isEnd = 0;
    T_prime *tbPrime = malloc(nSqrt * sizeof(tbPrime[0])) ;
    int32_t *offSet = malloc(nSqrt * sizeof(offSet[0])) ;
    int32_t sizeTable =  (nSqrt < 32768) ? nSqrt : 32768 ;
    int8_t *isComposed = calloc( sizeTable , sizeof(isComposed[0])) ;
    uint32_t nbPrime = 0 ;
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
                    offSet[nbPrimeSqrt++] = (int32_t)(icp2 - sizeTable) ;
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
        if ( offSetTable + sizeTable > nbMax) { sizeTable = (int32_t)(nbMax - offSetTable) ; }
        memset(isComposed,0,sizeTable) ;
        {
            int np ;
            for(np=0;np<nbPrimeSqrt;np++) {
                T_prime p = tbPrime[np] ;
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

uint32_t FindPrime(T_prime maxValue,void *ctx,TY_CPL_nxtPrime nxtPrime) {
#if USE_PRIMEb
    return FindPrime_b(maxValue,ctx,nxtPrime) ;
#else
    return FindPrime_a(maxValue,ctx,nxtPrime) ;
#endif
}

CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) {
    CTX_PRIMETABLE *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    ctx->maxValue = maxValue ;
    if(maxValue > 100) {
        ctx->maxNbPrime = (uint32_t) (1+ maxValue / (log((double)maxValue) - 4)) ;
    } else {
        ctx->maxNbPrime = 30 ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime(ctx) ; }
    FindPrime(maxValue,ctx,CPL_tablePrime) ;
    return ctx ;
}

CTX_PRIMETABLE * Gen_tablePrimeNb(T_prime maxNb) {
    CTX_PRIMETABLE *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    if(maxNb < 30)  {
        ctx->maxNbPrime = 30 ;
    } else {
        ctx->maxNbPrime = (uint32_t) maxNb ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime(ctx) ; }
//   ctx->maxValue = 0x7fffffff ;
    ctx->maxValue = ctx->maxNbPrime * (4 + log (ctx->maxNbPrime))  ;
    FindPrime(ctx->maxValue,ctx,CPL_tablePrime) ;
    return ctx ;
}

static int32_t CmpPrime(const void *p1,const void *p2) {
    return *((int32_t *)p1) - *((int32_t *)p2) ;
}

int Search_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) {
    return (bsearch(&n,ctxP->tbPrime,ctxP->nbPrime,sizeof(n),CmpPrime) != NULL) ;
}

int SearchRg_TablePrime(CTX_PRIMETABLE *ctxP, T_prime n) {
    T_prime * pt= bsearch(&n,ctxP->tbPrime,ctxP->nbPrime,sizeof(n),CmpPrime) ;
    if(pt != NULL) return (int)(pt-ctxP->tbPrime) ;
    else return -1 ;
}
