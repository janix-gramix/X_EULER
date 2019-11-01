
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

int PB152c_n(PB_RESULT *pbR, int n) ;

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
    PB152c_n(&pbR,n);
    
    fprintf(stdout,"\tPB%s(%.06fs) Sol=%s \n",pbR.ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes) ;
    return 0 ;
    
}


uint64_t PGCD64(uint64_t n1,uint64_t n2 ) {
    uint64_t r ;
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
// ppcm*ppcm/k**2
typedef struct Element152_c {
    int     val ;           // k
    sumType weight ;        // contribution to 1/2 pondered by ppcm
} Element152_c ;


typedef struct hist152 {
    sumType     sum ;   // sigma(1/k**2) pondered
    countType count ;   // times sum has been observed
} hist152 ;

// comparison on sum
int CmpHisto(const void *el1,const void *el2){
    hist152 * h1 = (hist152 *)el1 ;
    hist152 * h2 = (hist152 *)el2 ;
    if(h1->sum > h2->sum) {
        return 1 ;
    } else if (h1->sum < h2->sum) {
        return -1 ;
    } else return 0 ;
}

// compress array of histos by merging histo with same sum (add counts)
int CompressHisto(hist152 *Histo, int nbHisto) {
    int nbHistoMerge = 0 ;
    qsort(Histo,nbHisto,sizeof(Histo[0]),CmpHisto);
    int ih,j;
    hist152 *sortHisto = Histo ;
    for(ih=0,nbHistoMerge=0;ih<nbHisto;) {
        sumType sum = Histo[ih].sum ;
        sortHisto[nbHistoMerge] = Histo[ih] ;
        for(j=ih+1;j<nbHisto && Histo[j].sum == sum;j++) {
            //          sortHisto[nbHistoMerge].count += Histo[j].count ;
            addCount(sortHisto[nbHistoMerge].count , Histo[j].count) ;
        }
        ih = j ;
        nbHistoMerge++ ;
    }
    return nbHistoMerge ;
}

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

// for development for the set of elements(coefficients) corresponding to the same Powp
// As the sum must comply to modulus, the sum is decomposed in quotient (intSum) and remainder (modSum)
typedef struct sumLevel {
    sumType        intSum ;
    constraintType modSum ;
} sumLevel ;



// results of the development for  Powp
typedef struct LevelDev {
    indexType           nbSum ;     // number of resulting sums
    constraintType      constraint ;    // constraint (modulus = p*p)
    indexType           nbDiffMod ;     // number of different remainder
    sumLevel *          sumL ;      // array of resulting sums
    indexType *         indMod ;    // for each possible remainder index in sumL for the first
    // sumL with this remainder.
    sumType             maxSum ;    // max intSum obtain for the development
} LevelDev ;

// current context of development
typedef struct levelState {
    int         numLevel ;  // index of treated powp (with nb elements > 0)
    Element152_c  Elem[50] ;  // coefficients for current powp
    int         nbElem ;    // number of coefficients
    constraintType  constraint ;    // current constraint for the level (= p*p)
    int64_t     ppcm ;  // ppcm of elements for the level
    Powp152     PP ;    // current Powp
    int32_t     fact ;  // coefficient factor to convert sum from the precedent level (must apply fact*fact)
    int     computVal[PB152_MAXN] ; // history of k precedently used (to filter coefficients and use only once)
    int     nbComputVal ;   // number of k precedently used
    int64_t cumPpcm ;       // ppcm of precedent level ppcm's. (used to compute ppcm for the current level)
} levelState ;

// comparison of sumL, first by remainder modSum, secondly by increasing quotient
int CmpSumL(const void *el1,const void *el2){
    sumLevel * sumL1 = (sumLevel *)el1 ;
    sumLevel * sumL2 = (sumLevel *)el2 ;
    if(sumL1->modSum > sumL2->modSum) {
        return 1 ;
    } else if (sumL1->modSum < sumL2->modSum) {
        return -1 ;
    } else if (sumL1->intSum > sumL2->intSum) {
        return 1 ;
    } else if (sumL1->intSum < sumL2->intSum) {
        return -1 ;
    } else {
        return 0 ;
    }
}

LevelDev * Free_levelDev(LevelDev * Lv) {
    if(Lv) {
        free(Lv->sumL) ;
        free(Lv->indMod) ;
        free(Lv) ;
    }
    return NULL ;
}
// Compute all the partial sums for a powp.
// For nbElem (coefficients) => 2**nbElem partial sums.
// The sum are sorted and indexed by remainder (mod[constraint]) for faster access.
LevelDev * Cmp52_level (const Element152_c * Elem, int nbElem, int32_t modConstraint) {
    indexType nbSum = 1 << nbElem ;
    LevelDev * Lv = calloc(1,sizeof(Lv[0])) ;
    Lv->nbSum = nbSum ;
    Lv->constraint = modConstraint ;
    Lv->sumL = malloc((nbSum+1) * sizeof(Lv->sumL[0])) ;
    Lv->indMod = calloc(modConstraint,sizeof(Lv->indMod[0])) ;
    int i ;
    Lv->sumL[0].intSum = 0 ;
    Lv->sumL[0].modSum = 0 ;
    sumType MaxSum = 0 ;
    sumLevel * sumL = Lv->sumL ;
    for(i=0;i<nbElem;i++) { // add element by Element
        sumType intSum = Elem[i].weight / modConstraint ;
        constraintType modSum = (constraintType) (Elem[i].weight - intSum * modConstraint) ;
        int j,jmax = 1 << i ;
        for(j=0;j<jmax;j++) {
            sumL[j+jmax].intSum = sumL[j].intSum+intSum ;
            sumL[j+jmax].modSum = sumL[j].modSum+modSum ;
            if(sumL[j+jmax].modSum-modConstraint >= 0 ) { // avoid division
                sumL[j+jmax].modSum = sumL[j+jmax].modSum-modConstraint ;
                sumL[j+jmax].intSum++ ;
            }
            if(sumL[j+jmax].intSum  > MaxSum) MaxSum = sumL[j+jmax].intSum ;
        }
    }
    Lv->sumL[nbSum].intSum = 0 ;
    Lv->sumL[nbSum].modSum = modConstraint ;
    Lv->maxSum = MaxSum ;
    qsort(sumL,nbSum,sizeof(sumL[0]),CmpSumL) ; // sort by remainder,quotient
    Lv->nbDiffMod = 0 ;
    for(i=0;i<nbSum;) { // compute access index for each remainder present
        int j , nb ;
        constraintType mod = sumL[i].modSum ;
        Lv->indMod[mod] = i+1 ;
        Lv->nbDiffMod++ ;
        for(j=i+1;j<nbSum && sumL[j].modSum== mod;j++) ;
        i = j ;
    }
    return Lv ;
}

// update context levelState for the next Powp with nbElement (coef) > 0
// return 0 if last level
int nextLevelState(levelState *LS,listPowp *LP) {
    int k ;
    int oldLevel = LS->numLevel ;
    int nbElem=0 ;
    while(nbElem==0 && LP->curPowp < LP->nbPowp) { // loop on remaining powp
        Powp152 * newPowp = LP->Powp + LP->curPowp++ ;
        int64_t powp = newPowp->powp ;
        int64_t ppcm = LS->cumPpcm ;
        for(k=(int)powp;k<=LP->max_n;k+=powp) { // search element k multiple of powp
            int l ;
            for(l=0;l<LS->nbComputVal;l++) { // verify if k not precendently used.
                if(LS->computVal[l] == k) break ;
            }
            if(l<LS->nbComputVal) continue ;
            if((LP->den_ppcm % k) == 0) {   // check if k is admissible
                ppcm = ppcm * k  / PGCD64(ppcm,k) ; // update ppcm
                LS->Elem[nbElem++].val = k ; // store k
            }
        }
        if(nbElem > 0) { // powp has one admissible element
            // update the levelState
            LS->constraint = newPowp->p * newPowp->p ;  // constraint
            
            LS->fact = (int32_t)(ppcm * LS->PP.p)/ LS->ppcm;    // fact from precedent level
            LS->cumPpcm = ppcm/PGCD64(newPowp->ppcmPowp,ppcm) ; // update cumulative ppcm
            LS->nbElem = nbElem ;
            LS->ppcm = ppcm ;   // ppcm for the level
            LS->PP = *newPowp ; // valid powp
            int ie ;
            for(ie=0;ie<nbElem;ie++) {  // compute weight of each element with ppcm
                LS->Elem[ie].weight = (ppcm / LS->Elem[ie].val)*(ppcm / LS->Elem[ie].val) ;
                LS->computVal[LS->nbComputVal++] = LS->Elem[ie].val ;
            }
            LS->numLevel++ ; // increment level number (can be decorralated qith powp number, due to powp with no new lement
        }
    }  ;
    return LS->numLevel-oldLevel ;
}

// decreasing order
int CmpOrderByPowp(const void *el1,const void *el2) {
    return (int) ((Powp152 *)el2)[0].powp - (int) ((Powp152 *)el1)[0].powp ;
}
// reoder powp by decreasing powp (and not p,pwop )
void REORDER152 (listPowp *LP) {
    qsort(LP->Powp+LP->curPowp,LP->nbPowp-LP->curPowp,sizeof(LP->Powp[0]),CmpOrderByPowp) ;
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

#define PB152_SIZE_HISTO    20480000



int PB152c_n(PB_RESULT *pbR,int max_n) {
    pbR->nbClock = clock() ;
    listPowp LP ;
    ComputeListPowP(&LP, max_n) ;
    levelState LS ;
    LS.nbComputVal = 0 ; LS.numLevel = 0 ;
    LS.computVal[LS.nbComputVal++] = 2 ; // 1/2**2 is mandatory, so remove it
    LS.constraint = 1; LS.ppcm = 2 ;   LS.PP.p = 1;    LS.PP.powp = 2;
    LS.PP.ppcmPowp = 2;  LS.cumPpcm = 1; LS.fact = 1 ;
    countType nbSolTot = int2count(0) ;
    sumType sumFinal = 1 ;
    
    int nbHisto0 = 0 ;
    hist152 * Histo0 = malloc(1*sizeof(Histo0[0])) ;
    Histo0[nbHisto0].sum = sumFinal ;
    Histo0[nbHisto0++].count = int2count(1);
    int64_t nbCountS0 = 0 ;
    sumType sumMax ;
    LevelDev * Lv ;
    countType * countS0 = NULL ;
    sumType offsetS0 = 0 ;
#define M_IN_COUNT   1
#define M_OUT_COUNT  2
#define M_OUT_HISTO  0
    int mode = 0 ; // in and out by histo
    
    while(nextLevelState(&LS, &LP)) {
        int curLevel = LS.numLevel ;
        countType * countS1 = NULL ;
        hist152 * Histo1 = NULL ;
        constraintType constraint = LS.constraint ;
        int nbElem = LS.nbElem ;
        Lv = Cmp52_level(LS.Elem,nbElem,constraint) ;
        int64_t factS = LS.fact*LS.fact ;
        sumType maxLv = Lv->maxSum ;
        int nbHisto1 = 0 , nbNul = 0 ;
        int ih,ihMax ;
        if(mode & M_OUT_COUNT ) {
            ihMax = (indexType) nbCountS0 ;
        } else {
            Histo1 = malloc(PB152_SIZE_HISTO*sizeof(Histo1[0])) ;
            offsetS0 = Histo0[0].sum ;
            ihMax = nbHisto0 ;
            sumMax =  Histo0[nbHisto0-1].sum ;
        }
        sumMax =  (sumMax*factS)/constraint  ;
        int64_t nbCountS1 =(sumMax-offsetS0*factS/constraint+Lv->maxSum+1) ;
        if(nbCountS1-1 > sumMax) nbCountS1 =  sumMax+1 ;
        sumType offsetS1 = sumMax - (nbCountS1-1) ;
        if((mode == M_OUT_HISTO)  && nbCountS1 < 2 * nbHisto0 - 10  ) {
            mode = M_OUT_COUNT ; printf("Swap to mode OUT_COUNT at level %d\n",curLevel) ;
        }
        printf("%.3f Level=%d,Elems=[%d...%d] Mod[%d],ppcm=%lld,nbEl=%d,MaxLv=%lld,fact=%d HMax=%lld ExpHout=%lld Hin=%d"
               ,(clock() - pbR->nbClock) / (float) CLOCKS_PER_SEC
               ,curLevel,LS.Elem[0].val,LS.Elem[nbElem-1].val,constraint
               ,LS.ppcm,nbElem,maxLv,LS.fact ,sumMax,nbCountS1 , ihMax) ;
        if(mode & M_OUT_COUNT ) {
            countS1 = calloc(nbCountS1,sizeof(countS1[0])) ;
        }
        u_int64_t nout = 0 ;
        for(ih=0;ih<ihMax;ih++) {
            sumType sum  ;  countType count ;
            if(mode & M_IN_COUNT ) {
                count = countS0[ih] ;
                sum = (offsetS0+ih) ;
            } else {
                count = Histo0[ih].count ;
                sum = Histo0[ih].sum ;
            }
            sum *= factS ;
            sumType intSum = sum / constraint ;
            constraintType modSum = (constraintType) (sum - intSum * constraint) ;
            if(Lv->indMod[modSum]>0) {
                int32_t ind ;
                for(ind = Lv->indMod[modSum]-1 ; Lv->sumL[ind].modSum == modSum;ind++) {
                    sumType newSum = (intSum - Lv->sumL[ind].intSum) ;
                    //                   if(newSum < 0) { nbNul++ ;break ;} // sum exeeds 1/2
                    nout++ ;
                    if(mode & M_OUT_COUNT ) {
                        int nc = (indexType)(newSum - offsetS1) ;
                        addCount(countS1[nc],count) ;
                    } else {
                        Histo1[nbHisto1].sum = newSum ;
                        Histo1[nbHisto1++].count = count ;
                    }
                }
            }
            if(mode == M_OUT_HISTO ) {
                if(nout > nbCountS1 ) { // test swap to count mode ?
                    mode = M_OUT_COUNT ; ih=-1; nout=0 ;
                    printf("***");
                    countS1 = calloc(nbCountS1,sizeof(countS1[0])) ;
                    continue ;
                } else if( nbHisto1 > PB152_SIZE_HISTO - Lv->nbSum || ih==nbHisto0-1) {
                    // need compression ? end of level, or histo buffer out almost full
                    nbHisto1 = CompressHisto(Histo1,nbHisto1) ;
                }
            }
        } // end loop ih on precedent active sums
        Lv = Free_levelDev(Lv);
        if(mode == M_OUT_HISTO ) {
            nbHisto0 = nbHisto1 ;
            free(Histo0) ;
            Histo0 = Histo1 ;
        } else {
            if(mode == M_OUT_COUNT) {
                REORDER152 (&LP);
                mode |= M_IN_COUNT ;
            }
            nbCountS0 = nbCountS1 ;
            offsetS0 = offsetS1 ;
            free(countS0);
            countS0 = countS1 ;
        }
        if(mode == M_OUT_HISTO ) {
            printf(" Hout=%lld/->%d \n",nout,nbHisto0) ;
            
#if defined(PB152_DEBUG)
            printf("(%lld,%lld)...(%lld,%lld)\n",(u_int64_t)(Histo0[0].sum),l_count(Histo0[0].count),(u_int64_t)Histo0[nbHisto0-1].sum,l_count(Histo0[nbHisto0-1].count));
            if(nbHisto0 < 50) {
                for(ih=0;ih<nbHisto0;ih++) printf("(%lld,%lld)",(u_int64_t)(Histo0[ih].sum),l_count(Histo0[ih].count));
                printf("\n");
            }
#endif
        } else {
            printf("=%d+%d(0)  Hout=%lld/%d\n",(indexType)(ihMax-nbNul),nbNul,nout,(indexType)nbCountS0);
#if defined(PB152_DEBUG)
            if(nbCountS0 < 50) {
                int i,i0 ;
                for(i0=0,i=0;i<nbCountS0 && i0 < 100;i++) { if(!isCountNull(countS0[i]))  { printf("[%d,%lld]",i,l_count(countS0[i])); i0++ ; } }
                printf("\n");
            }
#endif
        }
    }
    if(mode & M_OUT_COUNT) {
        int i ;
        //       nbSolTot = countS0[(int64_t)-offsetS0] ;
        nbSolTot = countS0[0] ;
        free(countS0);
    } else {
        int i ;
        for(i=0;i<nbHisto0;i++) {
            if(Histo0[i].sum == 0 ) { addCount(nbSolTot,Histo0[i].count) ; }
        }
        free(Histo0) ;
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
# define POW10_9  1000000000LL
    if(h_count(nbSolTot) || l_count(nbSolTot) >= POW10_9 * POW10_9) {
        u_int64_t nbH10 = h_count(nbSolTot) / POW10_9 ;
        u_int64_t nbHr = h_count(nbSolTot) % POW10_9 ;
        u_int64_t nbL10 = l_count(nbSolTot) / POW10_9 ;
        u_int64_t nbLr = l_count(nbSolTot) % POW10_9 ;
        u_int64_t p64_10 = (1LL << (64-9)) / (POW10_9/(1<<9)) ;
        u_int64_t p64r = ((1LL << (64-9)) - p64_10 * (POW10_9/(1<<9))) << 9 ;
        // (nbH10 * POW10_9 + nbHr) * (P64_10 * POW10_9 + p64r) + (nbL10 * POW10_9 + nbLr)
        u_int64_t nbLow = (nbHr * p64r + nbLr) ;
        u_int64_t nbHig = (nbLow / POW10_9) + nbH10 * p64r + nbHr * p64_10 + POW10_9 * nbH10 * p64_10 + nbL10;
        nbLow -= (nbLow / POW10_9) * POW10_9 ;
        
        if(pbR->isVerbose) fprintf(stdout,"\t PB%s 1/2=sigma(1/n**2 <1<n<=%d has %lld%lld sol\n",pbR->ident,max_n,nbHig,nbLow);
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld%lld",nbHig,nbLow) ;
    } else {
        if(pbR->isVerbose) fprintf(stdout,"\t PB%s 1/2=sigma(1/n**2 <1<n<=%d has %lld sol \n",pbR->ident,max_n,l_count(nbSolTot));
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",l_count(nbSolTot)) ;
    }
    return 1 ;
}


