

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

int PB152a_n(PB_RESULT *pbR, int n) ;

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
    PB152a_n(&pbR,n);
    
    fprintf(stdout,"\tPB%s(%.06fs) Sol=%s \n",pbR.ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes) ;
    return 0 ;
    
}

#define PB152_MAXN 320



typedef int64_t    sumType ;
typedef int32_t     constraintType ;
typedef int         indexType ;


#if !defined(INT128)
typedef struct countType {
    u_int64_t   h ;
    u_int64_t l ;
} countType ;
#define addCount(dst,src)   if( ( (dst).l = (dst).l + (src).l ) >= (src).l ) (dst).h += (src).h ; else (dst).h += (src).h + 1
#define int2count(low)   ((countType){0,(low)})
#define l_count(c)       ((c).l)
#define h_count(c)       ((c).h)
#define isCountNull(c)     ((c).l==0  && (c).h==0 )
#  else
typedef __int128   countType ;
#define addCount(dst,src)  (dst) += (src)
#define int2count(low)   ((countType)(low))
#define l_count(c)       ((u_int64_t)(c))
#define h_count(c)       ((int64_t)((c)>>64))
#define isCountNull(c)     ((c)==0)
#endif

typedef struct Node152_a {
    int  level ;        // level for the node
    indexType elem ;              // element in the level for the node
    constraintType sumConstraint   ;   // cumulates sum modulus the constraint
    sumType sum ;                       // cumulated sum from the beginning
} Node152_a ;


typedef struct hist152 {
    sumType     sum ;   // sigma(1/k**2) pondered
    countType count ;   // times sum has been observed
} hist152 ;


typedef void (*CBsol152) (void *ctx,hist152 hist) ;


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



int PB152a_n(PB_RESULT *pbR,int max_n) ;



// update context levelsGlobal for the next Powp with nbElement (coef) > 0
// return 0 if last level
int nextLevelsGlobal(levelsGlobal *LG,listPowp *LP) {
    int k ;
    int oldLevel = LG->nbLevel ;
    level_a * LS = LG->LS+oldLevel ;
    int nbElem=0 ;
    while(nbElem==0 && LP->curPowp < LP->nbPowp) { // loop on remaining powp
        Powp152 * newPowp = LP->Powp + LP->curPowp++ ;
        int64_t powp = newPowp->powp ;
        int64_t ppcm = LP->den_ppcm ;
        for(k=(int)powp;k<=LP->max_n;k+=powp) { // search element k multiple of powp
            int l ;
            for(l=0;l<LG->nbComputVal;l++) { // verify if k not precendently used.
                if(LG->computVal[l] == k) break ;
            }
            if(l<LG->nbComputVal) continue ;
            if((LP->den_ppcm % k) == 0) {   // check if k is admissible
                LS->Elem[nbElem++].val = k ; // store k
            }
        }
        if(nbElem > 0) { // powp has one admissible element
            // update the levelState
            LS->constraint = newPowp->ppcmPowp * newPowp->ppcmPowp ;  // constraint
            
            LS->nbElem = nbElem ;
            LS->PP = *newPowp ; // valid powp
            int ie ;
            for(ie=0;ie<nbElem;ie++) {  // compute weight of each element with ppcm
                u_int64_t weight = (ppcm / LS->Elem[ie].val)*(ppcm / LS->Elem[ie].val)  ;
                
                LS->Elem[ie].weight = weight  ;
                LS->Elem[ie].weightConstraint = weight % LS->constraint ;
                
                LG->computVal[LG->nbComputVal++] = LS->Elem[ie].val ;
            }
            LS->numLevel = LG->nbLevel ; // increment level number (can be decorralated qith powp number, due to powp with no new lement
            LG->nbLevel++ ;
        }
    }  ;
    return LG->nbLevel-oldLevel ;
}

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


int PB152a_n(PB_RESULT *pbR,int max_n) {
    u_int64_t den_ppcm = 1  ;
    int nbSol = 0 ;
    pbR->nbClock = clock() ;
    listPowp LP ; // list of power of p admissible
    levelsGlobal LG ; // context for levels
    ComputeListPowP(&LP, max_n) ; // compute powp
    den_ppcm = LP.den_ppcm ;
    LG.nbLevel =0 ;
    LG.nbComputVal=0 ;
    LG.computVal[LG.nbComputVal++] = 2 ;
    // compute all the level
    while(nextLevelsGlobal(&LG,&LP) ) {} //compute coefficients for all the powp
    
    printf("=%lld [%lld] nbLevel%d\n",LP.den_ppcm,LP.den_ppcm*LP.den_ppcm/2,LG.nbLevel );
    Node152_a nod[PB152_MAXN] ; // nodes for tree walk
    int is ;
    //   is =0 ; nod[is].sum = den_ppcm*den_ppcm/2 - Elem[0].weight ;  // car 2 obligatoire
    is =0 ; nod[is].sum = den_ppcm*den_ppcm/4 ;  // car 2 obligatoire
    nod[is].sumConstraint = 0 ;
    nod[is].level = -1 ; // special case for 1/2**2
    is++ ;
    nod[is].level = 0 ;
    nod[is].elem = 0 ;
    
    while(is>=1) {
        if(nod[is-1].sum == 0) {
            nbSol++ ;
        }
        while (is > 0) { //loop : try to advance by developing current node
            if(nod[is-1].sum > 0 && nod[is].level < LG.nbLevel) { // precedent node OK ?
                int ie ;
                level_a *curLevel = LG.LS+nod[is].level ;
                for(ie=nod[is].elem; ie<curLevel->nbElem ; ie++) { // search a valid element
                    if( curLevel->Elem[ie].weight <= nod[is-1].sum ) {
                        nod[is].elem = ie ;
                        nod[is].sum = nod[is-1].sum-curLevel->Elem[ie].weight ;
                        if(nod[is].level == nod[is-1].level) {
                            nod[is].sumConstraint = nod[is-1].sumConstraint - curLevel->Elem[ie]. weightConstraint;
                            if(nod[is].sumConstraint < 0) nod[is].sumConstraint += curLevel->constraint ;
                        } else {
                            nod[is].sumConstraint = nod[is].sum % curLevel->constraint  ;
                        }
                        break ;
                    }
                }
                if(ie < curLevel->nbElem) { // found an admissible element
                    is++ ; nod[is].elem = nod[is-1].elem+1 ; nod[is].level=nod[is-1].level ;
                    break ; // continue on current level
                } else if(nod[is-1].sumConstraint == 0){
                    nod[is].level++ ; // current level exhausted, try next level
                    nod[is].elem = 0 ;
                    break ; // continue on next level
                }
            } // end ob branch (last level, or sum< 0, must rewind
            
            while( --is > 0) { // rewind the tree
                if(nod[is].elem+1 <  LG.LS[nod[is].level].nbElem) {
                    nod[is].elem++ ; break ;
                } else if ( nod[is-1].sumConstraint == 0) {
                    nod[is].level++ ; nod[is].elem = 0 ;  break ;
                }
            }
        } // end loop try to advance
    } // end global loop
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s 1/2=sigma(1/n**2 <1<n<=%d has %d sol \n",pbR->ident,max_n,nbSol);
    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol) ;
    return 1 ;
}


