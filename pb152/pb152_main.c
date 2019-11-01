#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PB_MAX_STRLEN   100

typedef struct PB_RESULT {
    const char    *ident ;
    int     isVerbose ;
    char    strRes[PB_MAX_STRLEN] ;
    clock_t nbClock ;
} PB_RESULT ;

int PB152_n(PB_RESULT *pbR, int n) ;

int main(int argc, char**argv) {
    PB_RESULT pbR ;
    pbR.ident = "PB152" ;
    pbR.isVerbose = 1 ;
    int n = 80 ;
    if(argc > 1) {
        n = atoi(argv[1]);
    }
    if(n > 320) n = 320 ;
    printf("Execution for k <= %d\n",n);
    PB152_n(&pbR,n);
    
    fprintf(stdout,"\tPB%s(%.06fs) Sol=%s \n",pbR.ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes) ;
    return 0 ;
    
}


#define PB152_MAXN 320

typedef int64_t    sumType ;
typedef int32_t     constraintType ;
typedef int         indexType ;


// Power of p (prime) admissible for a development
// np design the max n  with p**n admissible
typedef struct Powp152 {
    constraintType p ;  // prime p
    constraintType powp ;   // p**n
    constraintType ppcmPowp ; // p**(np+1-n) usefull to cumulate ppcm
} Powp152 ;

// list of Powp for the maxValue of k
typedef struct ListPowp {
    int      max_n ;     // max value for ceofficient
    Powp152 Powp[100];
    int     nbPowp ;    // number of Powp
    int64_t den_ppcm ;  // ppcm of the whole powp (not square)
    int curPowp ;       // current powp to treat.
} listPowp ;


typedef struct Node152 {
    indexType elem ;              // element in the level for the node
    sumType sum ;                 // cumulated sum from the beginning
} Node152 ;

typedef struct Element152 {
    int     val ;           // k
    sumType weight ;        // contribution to 1/2 pondered by ppcm
    sumType cumToEnd ;
} Element152 ;

typedef struct Element152_a {
    int val ;                           // int for denominator
    sumType weight ;                    // contribution to 1/2
    constraintType weightConstraint ;   // weight modulus the constraint
} Element152_a ;

// a context of development for a powp
typedef struct level_a {
    int         numLevel ;  // index of treated powp (with nb elements > 0)
    Element152_a  Elem[50] ;  // coefficients for current powp
    int         nbElem ;    // number of coefficients
    constraintType  constraint ;    // current constraint for the level (= p*p)
    Powp152     PP ;    // current Powp
} level_a ;


// whole set of development of powp
typedef struct levelsGlobal {
    level_a  LS[100] ;
    int     nbLevel ;
    int     computVal[PB152_MAXN] ; // history of k precedently used (to filter coefficients and use only once)
    int     nbComputVal ;   // number of k precedently used
    
} levelsGlobal ;


// construct the list of admissible powp's for maxvalue for k
int ComputeListPowP(listPowp *LP, int maxValue){
    // maxValue must be <= 300
    int Primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,119,127,131,137,139,149,151,157,161} ;
    LP->den_ppcm = 1 ;
    LP->max_n = maxValue ;
    int ip,p ;
#define MAX152_NP   16
    int64_t invSqare[MAX152_NP], sumInv[1<<MAX152_NP] ;
    for(ip=0;(p=Primes[ip])<=maxValue/2;ip++) ; // seach first prime above maxValue/2 (other not admissible)
    int nbPowp = 0 ;
    for(;--ip >=0;) {
        p = Primes[ip] ;
        int powp,np ;
        for(powp=p;powp*p<=maxValue/2;) powp = powp * p ; // search max power of p candidate
        int notFound = 1 ;
        for(;powp>1 && notFound; powp /= p) {// search by decreasing power of the first admissible
            np = maxValue/powp ;
            int nbInv = 0 ;
            int sqp = p*p ;
            if(p != 2 && np>= powp-1) np = powp-1 ;
            if(p==2 && np>= 2*powp-1) np =2*powp -1 ;
            if(np>MAX152_NP) { printf("FATAL ERROR np=%d too big for powp=%d\n",np,powp) ; }
            int i ;
            for(i=1;i<=np;i++) { // for each i (k=i) search inverse for k*k mod[p*p]
                if((i % p) == 0) continue ;
                int i2 = (i*i) % sqp ;
                int inv_i2, ninv ; // search 1/(i*i) with if a power of i*i
                for(inv_i2 = 1; (ninv = (inv_i2*i2) % sqp )!=1;inv_i2=ninv) ;
                invSqare[nbInv++] = inv_i2 ; // store 1/(i*i) mod[p*p]
            }
            sumInv[0] = 0 ;
            for(i=0;i<nbInv && notFound ;i++) { // search a combinaison sum of inverse == 0 mod[p*p]
                int j,jmax = 1 << i ;
                for(j=0;j<jmax && notFound;j++) {
                    sumInv[j+jmax] = (sumInv[j]+invSqare[i]) % sqp ;
                    if(sumInv[j+jmax]==0) notFound = 0 ;
                }
            }
            
        }
        if(!notFound) { // existing combinaison sum
            powp *= p ; // add p*powp and all divisor to list of admissible powp
            LP->den_ppcm *= powp ; printf("x%d",powp) ;
            constraintType ppcmPowp = p ;
            for(;powp > 1;powp /= p) {
                LP->Powp[nbPowp].powp = powp ;
                LP->Powp[nbPowp].p = p ;
                LP->Powp[nbPowp++].ppcmPowp = ppcmPowp ; ppcmPowp *= p ;
            }
        }
    }
    LP->curPowp = 0 ;
    LP->nbPowp = nbPowp ;
    printf("=%llu [%llu] \n",(u_int64_t)LP->den_ppcm,(u_int64_t)(LP->den_ppcm *  LP->den_ppcm/2) );
    return LP->nbPowp ;
}



int PB152_n(PB_RESULT *pbR,int max_n) {
    int nbSol = 0 ;
    pbR->nbClock = clock() ;
    listPowp LP ; // list of power of p admissible
    levelsGlobal LG ; // context for levels
    ComputeListPowP(&LP, max_n) ; // compute powp
    u_int64_t den_ppcm = LP.den_ppcm  ;
    
    Node152 nod[PB152_MAXN] ;
    Element152 Elem[PB152_MAXN] ;
    int k ;
    int nbElem = 0 ;
    for(k=2;k<=max_n;k++) { // compute the weight for each k retained
        if((den_ppcm % k) == 0 ) {
            Elem[nbElem].weight = (den_ppcm / k)*(den_ppcm / k) ;
            Elem[nbElem++].val = k ;
        }
    }
    sumType cum = 0 ;   // cumulative to the end, so we can filter the partail sums too low.
    for(k=nbElem-1;k>=0;k--) {
        cum += Elem[k].weight ;
        Elem[k].cumToEnd = cum ;
    }
    printf("=%lld [%lld] nbNum=%d\n",den_ppcm,den_ppcm*den_ppcm/2,nbElem );
    int is ;
    is =0 ;
    nod[is].sum = den_ppcm*den_ppcm/2 - Elem[0].weight ;
    nod[is].elem = 0 ;
    is++ ;
    nod[is].sum = nod[is-1].sum - Elem[1].weight ;
    nod[is].elem = 1 ;
    while(is>=1) { // tree walk do develop the sums
        if(nod[is].sum > 0) {
            is++ ;
            nod[is].elem = nod[is-1].elem + 1;
            nod[is].sum = nod[is-1].sum - Elem[nod[is].elem].weight ;
        } else {
            if(nod[is].sum == 0) {
                nbSol++ ;
            }
            while (is > 0) {
                nod[is].sum += Elem[nod[is].elem].weight ;
                int ie ;
                for(ie=nod[is].elem+1;Elem[ie].cumToEnd >= nod[is].sum &&  ie<nbElem ;ie++) {
                    if( Elem[ie].weight <= nod[is].sum) break ;
                }
                if(ie<nbElem && Elem[ie].cumToEnd>= nod[is].sum ) {
                    nod[is].elem = ie ; nod[is].sum -= Elem[ie].weight ;
                    break ;
                }
                is-- ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol) ;
    return 1 ;
}

