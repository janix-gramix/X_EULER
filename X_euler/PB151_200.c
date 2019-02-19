//
//  PB151_200.c
//  X_euler
//
//  Created by Jeannot on 01/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "euler_utils.h"
#include "faray_utils.h"
#include "PB151_200.h"

//#define PB152_MAXN  45
#define PB152_MAXN 80

typedef __int128    sumType ;
typedef u_int64_t   countType ;
typedef int32_t     constraintType ;
typedef int         indexType ;

typedef struct Element152 {
    int val ;                           // int for denominator
//   int level ;                         // level for modulus constraint
//    indexType deltaLevelEnd  ;          // increment to end level element (0 <=> last)
    sumType weight ;                    // contribution to 1/2
    sumType cumWtoEnd ;                 // cumulated weight por remaning elments
    constraintType constraint ;         // modulus constraint
} Element152 ;


typedef struct hist152 {
 //   int64_t sum;
    sumType sum ;
    countType count ;
} hist152 ;

int CmpHisto(const void *el1,const void *el2){
    hist152 * h1 = (hist152 *)el1 ;
    hist152 * h2 = (hist152 *)el2 ;
    if(h1->sum > h2->sum) {
        return 1 ;
    } else if (h1->sum < h2->sum) {
        return -1 ;
    } else return 0 ;
}


typedef struct sumLevel {
    sumType intSum ;
    constraintType modSum ;
} sumLevel ;

typedef struct LevelDev {
    indexType           nbSum ;
    constraintType      constraint ;
    indexType           nbDiffMod ;
    sumLevel *          sumL ;
    indexType *         indMod ;
} LevelDev ;

typedef struct Level {
    indexType   firstElem ;
    indexType   lastElem ;
    constraintType  constraint ;
} Level ;
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

LevelDev * Cmp52_level (const Element152 * Elem, int nbElem, int32_t modConstraint) {
    indexType nbSum = 1 << nbElem ;
    LevelDev * Lv = calloc(1,sizeof(Lv[0])) ;
    Lv->nbSum = nbSum ;
    Lv->constraint = modConstraint ;
    Lv->sumL = malloc(nbSum * sizeof(Lv->sumL[0])) ;
    Lv->indMod = calloc(modConstraint,sizeof(Lv->indMod[0])) ;
    int i ;
    Lv->sumL[0].intSum = 0 ;
    Lv->sumL[0].modSum = 0 ;
    sumLevel * sumL = Lv->sumL ;
    for(i=0;i<nbElem;i++) {
        sumType intSum = Elem[i].weight / modConstraint ;
        constraintType modSum = (constraintType) (Elem[i].weight - intSum * modConstraint) ;
        int j,jmax = 1 << i ;
        for(j=0;j<jmax;j++) {
            sumL[j+jmax].intSum = sumL[j].intSum+intSum ;
            sumL[j+jmax].modSum = sumL[j].modSum+modSum ;
            if(sumL[j+jmax].modSum-modConstraint >= 0 ) {
                sumL[j+jmax].modSum = sumL[j+jmax].modSum-modConstraint ;
                sumL[j+jmax].intSum++ ;
            }
         }
    }
    qsort(sumL,nbSum,sizeof(sumL[0]),CmpSumL) ;
    Lv->nbDiffMod = 0 ;
    for(i=0;i<nbSum;) {
        int j ;
        constraintType mod = sumL[i].modSum ;
        Lv->indMod[mod] = i+1 ;
        Lv->nbDiffMod++ ;
        for(j=i+1;j<nbSum && sumL[j].modSum== mod;j++) ;
        i = j ;

    }
    return Lv ;
}


int PB152c(PB_RESULT *pbR) {
    int Primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107} ;
    sumType den_ppcm = 1  ;
    int nbSol = 0 ;
    pbR->nbClock = clock() ;
    int ip,p ;
    int64_t invSqare[200],candidate[200],powCand[200], sumInv[65536] ;
    int nbCand = 0 ;
    for(ip=0;(p=Primes[ip])<=PB152_MAXN/2;ip++) ;
    for(;--ip >=0;) {
        p = Primes[ip] ;
        int powp ;
        for(powp=p;powp*p<=PB152_MAXN/2;) powp = powp * p ;
        int nbS0 = 0 ;
        int np ;
        while(powp>1) {
            np = PB152_MAXN/powp ;
            int nbInv = 0 ;
            int sqp = p*p ;
            if(p != 2 && np>= powp-1) np = powp-1 ;
            if(p==2 && np>= 2*powp-1) np =2*powp -1 ;
            int i ;
            for(i=1;i<=np;i++) {
                if(i != 0 && (i % p) == 0) continue ;
                int i2 = (i*i) % sqp ;
                int inv_i2, ninv ; ;
                for(inv_i2 = 1; (ninv = (inv_i2*i2) % sqp )!=1;inv_i2=ninv) ;
                invSqare[nbInv++] = inv_i2 ;
            }
            sumInv[0] = 0 ;
            for(i=0;i<nbInv;i++) {
                int j,jmax = 1 << i ;
                for(j=0;j<jmax;j++) {
                    sumInv[j+jmax] = (sumInv[j]+invSqare[i]) % sqp ;
                    if(sumInv[j+jmax]==0){  nbS0++ ;     break ;  }
                }
                if(nbS0) {  break ; }
            }
            if(nbS0) break ;
            else powp /= p ;
        }
        if(nbS0) {
            den_ppcm *= powp ; printf("x%d",powp) ;
            powCand[nbCand] = powp ;
            candidate[nbCand++] = p ;
        }
    }
    Element152 Elem[PB152_MAXN] ;
    hist152 * histo ;
 //   int lastElemByLevel[200] ;
 //   constraintType  constraintLevel[200] ;
    Level   level[200] ;
    int k ;
    int nbElem = 0 ;
    int ic ;
    //    int lastElemByLevel[100] ;
    int nbLevel = 0 ;
    u_int64_t nbSolTot = 0 ;
    int nbHisto = 0 ;
    histo = malloc(100000000*sizeof(histo[0])) ;
//    den_ppcm /= 2 ;
    level[nbLevel].firstElem = nbElem ;
    Elem[nbElem].val = 2;
    Elem[nbElem].weight = (den_ppcm / 2)*(den_ppcm / 2)  ;
//    Elem[nbElem].level = nbLevel ; // mandatory
    Elem[nbElem].constraint  = 1 ; // no constraint
//    constraintLevel[nbLevel] = 1 ;
    level[nbLevel].constraint = 1 ;
 //   Elem[nbElem].deltaLevelEnd = 0 ; // no other element
//    lastElemByLevel[nbLevel++] = nbElem ;
    level[nbLevel++].lastElem =nbElem ;
    nbElem++ ;
 
    for(ic=0;ic<nbCand;ic++) {
        //  for(ic=nbCand;--ic>=0;) {
        int p = candidate[ic] ;
        int powp = powCand[ic] ;
        //      if(p >= 5) continue ;
        int askPowp = p;
        while(powp > 1) {
            level[nbLevel].firstElem = nbElem ;
            for(k=powp;k<=PB152_MAXN;k+=powp) {
                int l ;
                for(l=0;l<nbElem;l++) {
                    if(Elem[l].val == k) break ;
                }
                if(l<nbElem) continue ;
                if((den_ppcm % k) == 0) {
                    Elem[nbElem].weight = (den_ppcm / k)*(den_ppcm / k) ;
                    Elem[nbElem].constraint = askPowp*askPowp  ;
   //                 Elem[nbElem].level = nbLevel ;
                    Elem[nbElem++].val = k ; printf("+%d ",k) ;
                }
            }
//            constraintLevel[nbLevel] = askPowp*askPowp  ;
            level[nbLevel].constraint = askPowp*askPowp ;
//            lastElemByLevel[nbLevel] = nbElem-1 ;
            level[nbLevel].lastElem = nbElem-1 ;
            nbLevel++ ;
            askPowp *= p ;
            powp /= p ;
        }
    }
    sumType cum = 0 ;
    for(k=nbElem-1;k>=0;k--) {
        cum += Elem[k].weight ;
        Elem[k].cumWtoEnd = cum ;
  //      Elem[k].deltaLevelEnd = lastElemByLevel[Elem[k].level]-k ;
  //      Elem[k].deltaLevelEnd = level[Elem[k].level].lastElem-k ;
    }
//    lastElemByLevel[nbLevel] = nbElem ;
    level[nbLevel].lastElem = nbElem ;
    Elem[nbElem].val = 0 ;
    Elem[nbElem].cumWtoEnd = 0 ;
    Elem[nbElem].weight = 0 ;
//    Elem[nbElem].level = 0 ;
    Elem[nbElem++].val = 0 ;
    
    sumType sumFinal = den_ppcm*den_ppcm/2 ;
    printf("=%llu [%llu] nbNum=%d\n",(u_int64_t)den_ppcm,(u_int64_t)(den_ppcm*den_ppcm/2),nbElem );
    int curLevel ;
    nbHisto = 0 ;
    histo[nbHisto].sum = sumFinal ;
    histo[nbHisto++].count = 1 ;
 //   int begNextLevel = 0 ;
  //  for(curLevel=1;curLevel<nbLevel;curLevel++) {
    for(curLevel=0;curLevel<nbLevel-1;curLevel++) {
        constraintType constraint = level[curLevel].constraint;
//       constraintType constraint = (constraintType) Elem[begNextLevel].constraint ;
        if(constraint != (constraintType) Elem[level[curLevel].firstElem].constraint) {
            printf("Level=%d %d=>%d)******\n",curLevel,constraint,(constraintType) Elem[level[curLevel].firstElem].constraint);
        }
//        int nbElem = Elem[begNextLevel].deltaLevelEnd+1 ;
       int nbElem = level[curLevel].lastElem - level[curLevel].firstElem +1 ;
//        LevelDev * Lv = Cmp52_level(Elem+begNextLevel,nbElem,constraint) ;
       LevelDev * Lv = Cmp52_level(Elem+level[curLevel].firstElem,nbElem,constraint) ;
        int nbInd = Lv->nbSum ;
        hist152 * nxtHisto =  histo + nbHisto ;
        printf("Level=%d,FisrtElem=[%d,%lld] const=%d, nbEl=%d, beg=%d nbSum=%d nbHisto=%d "
               ,curLevel,Elem[level[curLevel].firstElem].val,Elem[level[curLevel].firstElem].weight ,constraint,nbElem,level[curLevel].firstElem,Lv->nbSum ,nbHisto) ;
        int nxtNbHisto = 0 ;
//        sumType cumLevel = Elem[begNextLevel+nbElem].cumWtoEnd ;
        sumType cumLevel = Elem[level[curLevel+1].firstElem].cumWtoEnd ;
        int ih ;
        for(ih=0;ih<nbHisto;ih++) {
            sumType sum = histo[ih].sum ;
            countType count = histo[ih].count ;
            sumType intSum = sum / constraint ;
            constraintType modSum = (constraintType) (sum - intSum * constraint) ;
            int modComp = modSum ;
            if(Lv->indMod[modComp]>0) {
                int32_t ind ;
                for(ind = Lv->indMod[modComp] - 1 ;ind < nbInd && Lv->sumL[ind].modSum == modComp;ind++) {
                    nxtHisto[nxtNbHisto].sum = (intSum - Lv->sumL[ind].intSum) * constraint ;
                    if(nxtHisto[nxtNbHisto].sum < 0) break ; // on exede la valeur
                    nxtHisto[nxtNbHisto].count = count ;
                    if( cumLevel < nxtHisto[nxtNbHisto].sum) continue ;
                    nxtNbHisto++ ;
                }
            }
            if(nxtNbHisto > 2000000) {
                qsort(nxtHisto,nxtNbHisto,sizeof(nxtHisto[0]),CmpHisto);
                int ih1,i2;
                hist152 *sortHisto = nxtHisto ;
                for(ih1=0,i2=0;ih1<nxtNbHisto;) {
                    int j ;
                    sumType sum = nxtHisto[ih1].sum ;
                    sortHisto[i2] = nxtHisto[ih1] ;
                    for(j=ih1+1;j<nxtNbHisto && nxtHisto[j].sum == sum;j++) {
                        sortHisto[i2].count += nxtHisto[j].count ;
                    }
                    ih1 = j ;
                    i2++ ;
                }
                printf("CMP(%d->%d)",nxtNbHisto,i2);
                nxtNbHisto = i2 ;
           }
        }
        
        qsort(nxtHisto,nxtNbHisto,sizeof(nxtHisto[0]),CmpHisto);
        int i2;
        hist152 *sortHisto = nxtHisto ;
        for(ih=0,i2=0;ih<nxtNbHisto;) {
            int j ;
            sumType sum = nxtHisto[ih].sum ;
            sortHisto[i2] = nxtHisto[ih] ;
            for(j=ih+1;j<nxtNbHisto && nxtHisto[j].sum == sum;j++) {
                sortHisto[i2].count += nxtHisto[j].count ;
            }
            ih = j ;
            i2++ ;
        }
        printf(" %d => nbHisto=%d \n",nxtNbHisto,i2) ;
  //      begNextLevel += nbElem ;
        nbHisto = i2 ;
        histo =sortHisto ;
        for(ih=0;ih<nbHisto ;ih++) printf("(%llx%llx,%lld)",(u_int64_t)(histo[ih].sum>>64),(u_int64_t)(histo[ih].sum),histo[ih].count) ;
        printf("\n");
    }
    int i ;
    for(i=0;i<nbHisto;i++) {
        if(histo[i].sum == 0 ) nbSolTot += histo[i].count ;
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s 1/2=sigma(1/n**2 <1<n<=%d has %d ,%lld sol \n",pbR->ident,PB152_MAXN,nbSol,nbSolTot);
    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbSolTot) ;
    return 1 ;
}


typedef struct Node152 {
    indexType elem ;              // element for the node
    constraintType sumConstraint   ;   // cumulates sum modulus the constraint
    sumType sum ;                       // cumulated sum from the beginning
} Node152 ;

typedef void (*CBsol152) (void *ctx,hist152 hist) ;

typedef struct Element152a {
    int val ;                           // int for denominator
    int level ;                         // level for modulus constraint
    indexType deltaLevelEnd  ;          // increment to end level element (0 <=> last)
    sumType weight ;                    // contribution to 1/2
    sumType cumWtoEnd ;                 // cumulated weight por remaning elments
    constraintType weightConstraint ;   // weight modulus the constraint
    constraintType constraint ;         // modulus constraint
} Element152a ;

int PB152a(PB_RESULT *pbR) {
    int Primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97} ;
    u_int64_t den_ppcm = 1  ;
    int nbSol = 0 ;
    pbR->nbClock = clock() ;
    int ip,p ;
    int64_t invSqare[100],candidate[100],powCand[100], sumInv[2048] ;
    int nbCand = 0 ;
    for(ip=0;(p=Primes[ip])<=PB152_MAXN/2;ip++) ;
    for(;--ip >=0;) {
        p = Primes[ip] ;
        int powp ;
        for(powp=p;powp*p<=PB152_MAXN/2;) powp = powp * p ;
        int nbS0 = 0 ;
        int np ;
        while(powp>1) {
            np = PB152_MAXN/powp ;
            int nbInv = 0 ;
            int sqp = p*p ;
            if(p != 2 && np>= powp-1) np = powp-1 ;
            if(p==2 && np>= 2*powp-1) np =2*powp -1 ;
            int i ;
            for(i=1;i<=np;i++) {
                if(i != 0 && (i % p) == 0) continue ;
                int i2 = (i*i) % sqp ;
                int inv_i2, ninv ; ;
                for(inv_i2 = 1; (ninv = (inv_i2*i2) % sqp )!=1;inv_i2=ninv) ;
                invSqare[nbInv++] = inv_i2 ;
            }
            sumInv[0] = 0 ;
            for(i=0;i<nbInv;i++) {
                int j,jmax = 1 << i ; ;
                for(j=0;j<jmax;j++) {
                    sumInv[j+jmax] = (sumInv[j]+invSqare[i]) % sqp ;
                    if(sumInv[j+jmax]==0){  nbS0++ ;     break ;  }
                }
                if(nbS0) {  break ; }
            }
            if(nbS0) break ;
            else powp /= p ;
        }
        if(nbS0) {
            den_ppcm *= powp ; printf("x%d",powp) ;
            powCand[nbCand] = powp ;
            candidate[nbCand++] = p ;
        }
    }
    Element152a Elem[PB152_MAXN] ;
    int k ;
    int nbElem = 0 ;
    int ic ;
    int lastElemByLevel[100] ;
    int nbLevel = 0 ;
    Elem[nbElem].val = 2;
    Elem[nbElem].weight = (den_ppcm / 2)*(den_ppcm / 2)  ;
    Elem[nbElem].level = nbLevel ; // mandatory
    Elem[nbElem].constraint  = 1 ; // no constraint
//    Elem[nbElem].LevelEndElem = nbElem ; // no other element
   Elem[nbElem].deltaLevelEnd = 0 ; // no other element
    nbElem++ ;
   lastElemByLevel[nbLevel++] = 1 ;
    for(ic=0;ic<nbCand;ic++) {
//  for(ic=nbCand;--ic>=0;) {
        int p = candidate[ic] ;
        int powp = powCand[ic] ;
       int askPowp = p;
        while(powp > 1) {
            for(k=powp;k<=PB152_MAXN;k+=powp) {
                int l ;
                for(l=0;l<nbElem;l++) {
                    if(Elem[l].val == k) break ;
                }
                if(l<nbElem) continue ;
                if((den_ppcm % k) == 0 ) {
                    Elem[nbElem].weight = (den_ppcm / k)*(den_ppcm / k) ;
                    Elem[nbElem].constraint = askPowp*askPowp  ;
                    Elem[nbElem].weightConstraint = Elem[nbElem].weight % Elem[nbElem].constraint ;
                    Elem[nbElem].level = nbLevel ;
                    Elem[nbElem++].val = k ; printf("+%d ",k) ;
                }
            }
            lastElemByLevel[nbLevel] = nbElem-1 ;
            nbLevel++ ;
            askPowp *= p ;
            powp /= p ;
        }
    }
    lastElemByLevel[nbLevel] = nbElem ;
    int64_t cum = 0 ;
    for(k=nbElem-1;k>=0;k--) {
        cum += Elem[k].weight ;
        Elem[k].cumWtoEnd = cum ;
//        Elem[k].LevelEndElem = lastElemByLevel[Elem[k].level] ;
       Elem[k].deltaLevelEnd = lastElemByLevel[Elem[k].level]-k ;
    }
    Elem[nbElem].val = 0 ;
    Elem[nbElem].cumWtoEnd = -1 ;
    Elem[nbElem].weight = 0 ;
    Elem[nbElem].level = 0 ;
    Elem[nbElem++].val = 0 ;
   printf("=%lld [%lld] nbNum=%d\n",den_ppcm,den_ppcm*den_ppcm/2,nbElem );
    Node152 nod[PB152_MAXN] ;
    int is ;
    is =0 ; nod[is].sum = den_ppcm*den_ppcm/2 - Elem[0].weight ;  // car 2 obligatoire
    nod[is].sumConstraint = 0 ;
    nod[is++].elem = 0 ;
    nod[is].sum = nod[is-1].sum - Elem[1].weight ;
    nod[is].sumConstraint = nod[is-1].sumConstraint - Elem[1].weightConstraint;
    if(nod[is].sumConstraint < 0) nod[is].sumConstraint += Elem[1].constraint ;
    nod[is].elem = 1  ;
    nod[is+1].elem = nod[is].elem+1 ;
    is++ ;
    while(is>=1) {
        if(nod[is-1].sum == 0) {
 //               int j ; for(j=0;j<is;j++) printf("%d ",Elem[nod[j].elem].val);  printf("\n") ;
            nbSol++ ;
        }
 //     int j ; for(j=0;j<is;j++) printf("%d ",Elem[nod[j].elem].val); printf("[%d] .\n",Elem[nod[is].elem].level) ;
        while (is > 0) {
//          int j ; for(j=0;j<=is;j++) printf("%d ",Elem[nod[j].elem].val); printf(" *\n") ;
            if(nod[is-1].sum > 0) {
                int ie ;
                int curLevel = Elem[nod[is].elem].level ;
                int isOk = 0 ;
                int isCumOk = 1 ;
                for(ie=nod[is].elem; Elem[ie].level == curLevel ; ie++) {
                    if(Elem[ie].cumWtoEnd < nod[is-1].sum) {
                         isCumOk= 0 ; break ;
                    }
                   if( Elem[ie].weight <= nod[is-1].sum ) {
                        nod[is].elem = ie ;
                        nod[is].sum = nod[is-1].sum-Elem[ie].weight ;
                        if(Elem[ie].level == Elem[nod[is-1].elem].level) {
                           nod[is].sumConstraint = nod[is-1].sumConstraint - Elem[ie].weightConstraint;
                           if(nod[is].sumConstraint < 0) nod[is].sumConstraint += Elem[ie].constraint ;
                        } else {
                            nod[is].sumConstraint = nod[is].sum % Elem[ie].constraint  ;
                        }
                        isOk = 1 ;
                        break ;
                    }
                }
                if(isOk) {
//                    if(ie < Elem[ie].LevelEndElem || (nod[is].sumConstraint) == 0)  {
                  if(Elem[ie].deltaLevelEnd >0 || (nod[is].sumConstraint) == 0)  {
                        is++ ; nod[is].elem = nod[is-1].elem+1 ;   break ;
     //               } else if((nod[is-1].sum % Elem[nod[is-1].elem].constraint) == 0){
                    } else if(nod[is-1].sumConstraint == 0){
                        nod[is].elem += Elem[nod[is].elem].deltaLevelEnd+1;   break ;
                    }
    //            } else if(isCumOk && (nod[is-1].sum % Elem[nod[is-1].elem].constraint) == 0 )  {
                } else if(isCumOk && nod[is-1].sumConstraint == 0 )  {
                    nod[is].elem += Elem[nod[is].elem].deltaLevelEnd+1;
                    break ; // next lvel
                }
            }
            while(--is > 0) {
                if(Elem[nod[is].elem].deltaLevelEnd > 0 || nod[is-1].sumConstraint == 0) {
                    nod[is].elem++ ; break ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s 1/2=sigma(1/n**2 <1<n<=%d has %d sol \n",pbR->ident,PB152_MAXN,nbSol);

    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol) ;
    return 1 ;
}



int PB152(PB_RESULT *pbR) {
    int Primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97} ;
    u_int64_t den_ppcm = 1  ;
    int nbSol = 0 ;
    pbR->nbClock = clock() ;
    int i,ip,p ;
    int64_t invSqare[100], sumInv[2048] ;
    for(ip=0;(p=Primes[ip])<=PB152_MAXN/2;ip++) {
        int powp ;
        for(powp=p;powp*p<=PB152_MAXN/2;) {
            powp = powp * p ;
        }
      int nbS0 = 0 ;
      while(powp>1) {
            int np = PB152_MAXN/powp ;
            int nbInv = 0 ;
           int sqp = p*p ;
            for(i=1;i<=np;i++) {
                if((i % p) == 0) continue ;
                int64_t i2 = (i*i) % sqp ;
                int64_t inv_i2, ninv ; ;
                for(inv_i2 = 1; (ninv = (inv_i2*i2) % sqp ) !=1;inv_i2=ninv) ;
                invSqare[nbInv++] = inv_i2 ;
            }
            sumInv[0] = 0 ;
              for(i=0;i<nbInv;i++) {
                int j,jmax = 1 << i ; ;
                for(j=0;j<jmax;j++) {
                    sumInv[j+jmax] = (sumInv[j]+invSqare[i]) % sqp ;
                    if(sumInv[j+jmax]==0){ nbS0++ ; break ; }
                }
                if(nbS0) break ;
            }
            if(nbS0) break ;
            else powp /= p ;
        }
        if(nbS0) {
            den_ppcm *= powp ;
            printf("x%d",powp) ;
        }
    }
    int64_t Num[PB152_MAXN] ;
    int indNum[PB152_MAXN] ;
    int64_t cumNum[PB152_MAXN] ;
    int k ;
    int nbNum = 0 ;
    for(k=2;k<=PB152_MAXN;k++) {
        if((den_ppcm % k) == 0 ) {
            Num[nbNum] = (den_ppcm / k)*(den_ppcm / k) ;
            indNum[nbNum++] = k ;
        }
    }
    int64_t cum = 0 ;
    for(k=nbNum-1;k>=0;k--) {
        cum += Num[k] ;
        cumNum[k] = cum ;
    }
    printf("=%lld [%lld] nbNum=%d\n",den_ppcm,den_ppcm*den_ppcm/2,nbNum );
    int64_t Sum[PB152_MAXN] ;
    int Ind[PB152_MAXN] ;
    int is ;
    is =0 ; Sum[is] = den_ppcm*den_ppcm/2 - Num[0] ;  // car 2 obligatoire
    Ind[is] = 0 ;
    is++ ;
    Sum[is] = Sum[is-1] - Num[1] ;
    Ind[1] = 1  ;
    while(is>=1) {
        if (Sum[is] > 0) {
            is++ ;
            Ind[is] = Ind[is-1]+1 ;
            Sum[is] = Sum[is-1] - Num[Ind[is]] ;
        } else {
            if(Sum[is] == 0) {
 //               int j ; for(j=0;j<=is;j++) printf("%d ",indNum[Ind[j]]);
 //               printf("\n") ;
                nbSol++ ;
            }
            while (is > 0) {
                Sum[is] += Num[Ind[is]] ;
                int i ;
                for(i=Ind[is]+1;cumNum[i]>= Sum[is] &&  i<nbNum;i++) {
                    if( Num[i] <= Sum[is]) break ;
                }
                if(i<nbNum && cumNum[i]>= Sum[is] ) {
                    Ind[is] = i ;
                    Sum[is] -= Num[i] ;
                    break ;
                } else {
                    is-- ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol) ;
    return 1 ;
}


// nb(tiles) = 4*n*p with n>p
// so count multiples k*p from k=p+1 to k<= N/p
#define PB173_NB    1000000
int PB173(PB_RESULT *pbR) {
    int N = PB173_NB/4 ;
    pbR->nbClock = clock() ;
    int n ;
    int64_t S = 0 ;
    S += N-1 ; // product 1xn
    for(n=2;n*(n+1)<=N;n++) {
        S += N/n - n ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S) ;
    return 1 ;
}

#define PB174_NB    1000000
// #define PB174_NB    100

int PB174(PB_RESULT *pbR) {
    int N = PB174_NB/4 ;
    pbR->nbClock = clock() ;
    uint8_t *nbProd=malloc((N+1)*sizeof(nbProd[0])) ;
    int n ;
    for(n=2;n<=N;n++) { nbProd[n] = 1 ; }
    for(n=2;n*n<N;n++) {
        int np ;
        for(np=n*(n+1);np<=N;np+=n) { nbProd[np]++ ; }
    }
    int S = 0 ;
    for(n=2;n<=N;n++) { if(nbProd[n]<=10)S++ ;  }
    free(nbProd) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}
// #define PB179_NB    1000
#define PB179_NB    10000000

int PB179(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    uint16_t *nbDiv=malloc((N+2)*sizeof(nbDiv[0])) ;
    int n,d ;
    nbDiv[2] = nbDiv[3] = 2 ;
    for(n=4;n<=N;n+=2) { // d=2
        nbDiv[n] = 4 ;
        nbDiv[n+1] = 2 ;
    }
    int dMax = Sqrt32(N);
    for(d=3;d<=dMax;d++) { // d > 2
        int np = d*d ;
        nbDiv[np]++ ;
        for(np += d ;np<=N;np+=d) {
            nbDiv[np] += 2 ;
            
        }
    }
    int S = 0 ;
    for(n=2;n<N;n++) {
        if(nbDiv[n]==nbDiv[n+1] ) {
            S++ ;
        }
        
    }
    free(nbDiv) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}

int PB179c(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    int dMax = Sqrt32(N) ;
    int sizeCache = 16384 ;
    sizeCache += sizeCache & 1 ; // sizeCache even
    uint16_t *nbDiv=malloc((sizeCache)*sizeof(nbDiv[0])) ; // number of divisor
    int32_t *indNxtFactor = malloc((dMax+1)*sizeof(indNxtFactor[0])) ; // indice of nxt multiple
    int n,d ;
    int offset = 6 ; // offset for the current cache pof nbDiv
    for(d=3;d<=dMax;d++) { // if n = d1*d2 d1<=d2
        indNxtFactor[d] = d*d - offset ;
    }
    int nbDivAnt = 2 ; // n=5
    int S = 1 ; // (2,3)
    for(;offset < N;offset += sizeCache ) { // loop, by moving the cache
         for(n=0;n<sizeCache;n+=2) { // special case for d=2
            nbDiv[n] = 4 ; // n even
            nbDiv[n+1] = 2 ; // n odd
        }
        int indCacheMax = sizeCache ;
        if(offset+sizeCache >N+1) { // check if last cache move
            indCacheMax = N+1 - offset ;
        }
        for(d=3;d<=dMax;d++) { // d=3,...dMax
            int np = indNxtFactor[d] ;
            if(np < indCacheMax) { //
                if(np == d*d - offset ) { // square
                    nbDiv[np]++ ; np += d ;
                }
                for(;np<indCacheMax;np += d ) {
                    nbDiv[np] += 2 ; // add product, 2 divisors
                }
            }
            indNxtFactor[d] = np - sizeCache  ;
        }
        for(n=0;n< indCacheMax ; n++) { // check nbDiv[n] == nbDiv[n-1]
            if(nbDiv[n]==nbDivAnt) S++ ;
            nbDivAnt = nbDiv[n] ;
        }
    }
    free(nbDiv) ; free(indNxtFactor) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}


int PB179a(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    uint16_t *nbDiv=malloc((N+1)*sizeof(nbDiv[0])) ;
    int n ;
    // pow of 2
    for(n=2;n<=N;n++) {
        nbDiv[n] = 2 - (n&1);
    }
    int pow2 ;
    for(pow2=4;pow2<=N;pow2 *=2){
        for(n=pow2;n<=N;n+=pow2) nbDiv[n]++ ;
    }
    int p ;
    int maxP = Sqrt32(N) ;
    for(p=3;p<=N;p+=2) {
        if(nbDiv[p]>1) continue ;
        int k ;
        for(n=p,k=1;n<=N;n+=p, k++) {
            if(k<p) {
                nbDiv[n] *= 2 ;
            } else {
                k = 0 ;
            }
        }
        if(p<=maxP) {
            int powp ,exp1 ;
            int maxPowp = N/ p ;
            for(powp=p*p,exp1=3;;exp1++,powp *=p) {
                for(n=(int)powp,k=1;n<=N;n+=(int)powp, k++) {
                    if(k<p) {
                        nbDiv[n] *= exp1 ;
                    } else {
                        k = 0 ;
                    }
                }
                if(powp>maxPowp) break ;
            }
        }
    }
    int S = 0 ;
    for(n=2;n<N;n++) {
       if(nbDiv[n]==nbDiv[n+1] ) {
            S++ ;
        }
        
    }
    free(nbDiv) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}

int PB179b(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    uint16_t *pDivMax=malloc((N+1)*sizeof(pDivMax[0])) ;
    int n ;
    
    pDivMax[1] = 1 ;
    for(n=2;n<=N;n++) { // case d= 2
        pDivMax[n] = (n&1) ? 0 : 2 ;
    }
    int p ;
    int maxP = Sqrt32(N) ;
    // search max prime divisor for each n
    for(p=3;p<=maxP;p+=2) { // loop for p even
        if(pDivMax[p]) continue ;
        for(n=p;n<=N;n+=p) { // multiple of
            pDivMax[n] = p ;
        }
    }
    for(n=3;n<=N;n++) {
        if(n == pDivMax[n] || (pDivMax[n] == 0) ) {
            // n , case n< maxP and n > maxP
            pDivMax[n] = 2 ;
        } else {
            p = pDivMax[n] ;
            int d,d1,exp1 ; // search max pow of p
            for(d=n/p,exp1=2;d1=d/p,d==d1*p;d = d1) {
                exp1++ ;
            }
            pDivMax[n] = pDivMax[d]*exp1 ; // as d and p prime
        }
    }
    int S = 0 ;
    for(n=2;n<N;n++) {
        if(pDivMax[n]==pDivMax[n+1] ) {
            S++ ;
        }
    }
    free(pDivMax) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}

#define PB181_NB 60
#define PB181_NW 40
#define indBW(nb,nw)    ((nb)*(PB181_NW+1)+(nw))

#define P2(nb,nw)   xP2[indBW((nb),(nw))]
#define PW2(nb,nw)   xPW2[indBW((nb),(nw))]
#define M2(nb,nw)   xM2[indBW((nb),(nw))]

// first compute number of decomposition for one color
// (developpment for Prod(1/(1-x**n)) = (1+x+x**2+x**3+...) (1+x**2+x**4+x**6...) ...(1+x**n+x**2n+x**3n...)
// with Pentagonal number theorem (Euler)
// For 2 colors it is the developpement P2(x,y) = Prod(1/1-x**m . y**n) (n,m != (0,0) )
// the trick is to remark that :
// the elements of a decomposition are in 2 separate class.
//  a) nb(B) >= nb(W)
//  b) nb(B) < nb(W)
// Un element of class a) can be decomposed in BW and B, so the class is generated by P2(x,xy)
// Class B is generated by PW2(xy,y) = Prod(1/1-x**m . y**n) (n >0 , m >= 0 )
// so first compute PW2(x,y) by DP (see PB181a) and then
// compute P2(x,y) = P2(x,xy) * PW2(xy,y)
// the gain
int PB181(PB_RESULT *pbR) {
    u_int64_t xP1[PB181_NB+1] ;
    u_int64_t xP2[(1+PB181_NB)*(1+PB181_NW)] ;
    u_int64_t xPW2[(1+PB181_NW)*(1+PB181_NW)] ;

    int m,n ;
    pbR->nbClock = clock() ;
    xP1[0] = xP1[1] = 1 ;
    for(n=2;n<=PB181_NB;n++) {
        int k  ,Pk = 1  ;
        u_int64_t P = 0 ;
        for(k=1; Pk <= n ; Pk += 3*k+1, k++ ) {
            if(k&1) {
                P += xP1[n - Pk] ;
                // car P-k = Pk+k
                if(Pk + k <= n ) P += xP1[n - Pk - k ] ;
            } else {
                P -= xP1[n - Pk] ;
                if(Pk + k <= n ) P -= xP1[n - Pk - k ] ;
            }
        }
        xP1[n] = P ;
    }
    
    for (m = 0; m <= PB181_NW; m++) {
        for (n = 0; n <= PB181_NW-m; n++) {
            PW2(m,n)=0 ;
        }
    }
    
    PW2(0,0)=1 ;
    for (m = 0; m < PB181_NW; m++) {
        // n=1 to force a W
        for (n = 1 ; n <= PB181_NW-m; n++) {
            int i,j ;
            PW2(m,n)++;
            for (i = m; i < PB181_NW; i++) {
                // n+1 for because j-èn must be strictly positive
                for (j = n+1; j <= PB181_NW-i; j++) {
                    PW2(i,j) += PW2(i - m,j - n);
                }
            }
        }
    }
    int nb ;
    for(nb=1;nb<=PB181_NB;nb++) {
        P2(nb,0) = xP1[nb] ;
        if(nb<=PB181_NW) P2(0,nb)=xP1[nb] ;
    }
    P2(0,0) = 1 ;

    for(m=1;m<=PB181_NB;m++) {
        int maxN = (m < PB181_NW ) ? m : PB181_NW ;
        for(n=1;n<=maxN;n++) {
            int m2,n1;
            u_int64_t P = P2(m-n,n) ;
            for(n1=0;n1<n;n1++) {
               for(m2=0;m2+n1<n;m2++) {
                  P += P2(m-m2-n1,n1)*PW2(m2,n-n1-m2) ;
               }
            }
            P2(m,n) = P ;
            if(m<=PB181_NW) P2(n,m) = P ;
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%llu",P2(PB181_NB,PB181_NW)) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


// Partition of multiset
// par recurrence sur la taille des sous-ensembles de (m+n) elements (mB,nW)
// que l'on ajoute aux denombrement des partitions ne contenant que des sous-ensembles de taille inferieure.
int PB181a(PB_RESULT *pbR) {
    u_int64_t xP2[(1+PB181_NB)*(1+PB181_NW)] ;
    int m,n;
    pbR->nbClock = clock() ;
    for (m = 0; m <= PB181_NB; m++) {
        for (n = 0; n <= PB181_NW; n++) {
            P2(m,n)=0 ;
        }
    }
    P2(0,0)=1 ;
    for (n = 0; n <= PB181_NW; n++) {
        for (m = (n) ? 0 : 1 ; m <= PB181_NB; m++) {
            int i,j ;
            for (j = n; j <= PB181_NW; j++) {
                for (i = m; i <= PB181_NB; i++) {
                    P2(i,j) += P2(i - m,j - n);
                }
            }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%llu",P2(PB181_NB,PB181_NW)) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

//#define PB187_MAX   2000000000
#define PB187_MAX   100000000
int PB187(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB187_MAX/2)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbFind =0 ;
    int i,j;
    nbFind += nbPrime ; // 2 * pi
    for(i=1;i<nbPrime;i++) {
        int maxPj = PB187_MAX/ tbPrime[i] ;
        for(j=i;tbPrime[j] <= maxPj;j++)  ;
        nbFind += j-i ;
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbFind) ;
    return 1 ;
}
int PB187a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB187_MAX/2)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbFind =0 ;
    int i;
    int maxPi = Sqrt32(PB187_MAX) ;
    int pi ;
    nbFind += nbPrime ; // 2 * pi
 //   printf("%d ",nbPrime);

    int j = nbPrime - 1 ;
    for(i=1;(pi=tbPrime[i])<= maxPi ;i++) {
        int maxPj = PB187_MAX/ tbPrime[i] ;
        for(;tbPrime[j] > maxPj;j--)  ;
 //       printf("%d ",j+1);
        nbFind += j-i+1 ;
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbFind) ;
    return 1 ;
}
//
// implantation de l'algo d'euler pour caluler PI(x)  S(x)= x - Sigma(X/pi) +  Sigma(X/pi pj) ...
// S(x) = PI(x) - PI(sqrt(x))+ 1 ; les pi <= sqrt(x)
// n'est rentable qua partir de 10**9
int PB187b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    int maxPi = Sqrt32(PB187_MAX)+1 ;
    if((ctxP = Gen_tablePrime(maxPi)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbFind =0 ;
    int i;
    int Pi ;
 //   nbFind += nbPrime ; // 2 * pi,
    int nbM = nbPrime -1 ;
    for(i=0;i<nbPrime ;i++) {
        Pi=tbPrime[i] ;
        int invPi = PB187_MAX/ Pi ;
        // on veut calculer pi(invPi)
        int PIinvPi = invPi ;
        int sqInvPi = Sqrt32(invPi) ;
        for(;tbPrime[nbM]>sqInvPi;nbM--) ;
        const T_prime *pt1,*pt2,*pt3,*pt4,*pt5,*pt6,*pt7,*pt8 ;
        const T_prime *ptnbM = tbPrime+nbM ;
        int32_t p1,p2,p3,p4,p5,p6,p7 ;
        
        for(pt1=tbPrime;pt1<=ptnbM && (p1=  invPi / *pt1);pt1++) {
            PIinvPi -= p1 ;
            for(pt2=pt1+1;pt2<=ptnbM &&(p1 >= *pt2);pt2++) {
                PIinvPi += (p2 = p1 / *pt2 )  ;
                for(pt3=pt2+1;pt3<=ptnbM &&(p2 >= *pt3);pt3++) {
                    PIinvPi -=  (p3  = p2 / *pt3)  ;
                    for(pt4=pt3+1;pt4<=ptnbM &&(p3 >= *pt4);pt4++) {
                        PIinvPi +=  (p4 = p3 / *pt4) ;
                        for(pt5=pt4+1;pt5<=ptnbM && (p4 >= *pt5);pt5++) {
                            PIinvPi -=  (p5 = p4 / *pt5 ) ;
                            for(pt6=pt5+1;pt6<=ptnbM && (p5 >= *pt6);pt6++) {
                                PIinvPi += (p6 = p5 / *pt6 )  ;
                                for(pt7=pt6+1;pt7<=ptnbM && (p6 >= *pt7);pt7++) {
                                    PIinvPi -=  (p7 = p6 / *pt7 ) ;
                                    for(pt8=pt7+1;pt7<=ptnbM && (p7 >= *pt8);pt8++) {
                                        PIinvPi +=  p7 / *pt8 ;                                     }
                                }
                            }
                        }
                    }
                }
            }
        }
 //       printf("PI(%d)=%d\n",invPi,PIinvPi+nbM);
        nbFind += PIinvPi+nbM - i  ; // -i +i PI(i)
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbFind) ;
    return 1 ;
}

#define PB188_VAL   1777
#define PB188_EXP   1855
#define PB188_EXP10 8

// use Euler's totient function
// if a^n = 1  a**(EPhi(n)) =1 Mod[n]
//EPhi(n) = n Prod (1-1/p)  p divide n.
// For m0= n = 10**8 Ephi(10**8) = 10**8 x 1/2 x 4/5 = 2**8 * 5**8 x 2/5 = 2**9 * 5**7 = 20000000
// Donc le plus petit m1 tel que a**m1=1 [m0] est un diviseur de 20000000, on trouve m1=1250000
// En suite on cherche le plus petit des diviseur de Ephi[125000] = 1250000 * 4/5 = 1000000
// tel que a**m2=1 [m1]
//
u_int64_t modPow(int64_t exp,u_int64_t mod) {
    u_int64_t modPOW2 = PB188_VAL ;
    u_int64_t answer = (exp & 1) ? (PB188_VAL % mod) : 1 ;
    int i ;
    for(i=1;exp >= (1LL<<i) ;i++) {
        modPOW2 = (modPOW2 * modPOW2) % mod ;
        if( (1LL<<i) & exp) {
            answer = (answer * modPOW2 ) % mod  ;
        }
    }
    return answer ;
}
int PB188(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nbMod =0 ;
    u_int64_t Pn,Dn ;
    u_int64_t Pn1,Dn1,rop ;
    int nbPow2 = 2 ;
    int nbPow5 = 0 ;
    Pn1= 0 ;
    Pn = 1<<nbPow2 ;
    Dn1=1 ;
    Dn= modPow (Dn1,Pn) ;
    rop= modPow (Pn1,Pn) ;
    fprintf(stdout,"(%d,%d) 1777**%lld %% %lld =%lld  ; 1777**%lld %% %lld = %lld \n",nbPow2,nbPow5,Pn1,Pn,rop,Dn1,Pn,Dn);

    while(nbPow2 < PB188_EXP10 || nbPow5 < PB188_EXP10) {
        Dn1=Dn ;
        Pn1=Pn ;
        nbPow5++ ;
        int deltaPow2 = (nbPow2+4 < PB188_EXP10) ? 4 : PB188_EXP10-nbPow2 ;
        nbPow2 += deltaPow2 ;
        Pn = Pn * 5 * (1<<deltaPow2) ;
        Dn=modPow(Dn,Pn) ;
        nbMod++ ;
        rop=modPow(Pn1,Pn) ;
        fprintf(stdout,"(%d,%d) 1777**%lld %% %lld =%lld  ; 1777**%lld %% %lld = %lld \n",nbPow2,nbPow5,Pn1,Pn,rop,Dn1,Pn,Dn);
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Exp=%d (%d)digits=%lld \n",pbR->ident,nbMod,PB188_EXP10,Dn);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%llu",Dn % 100000000 ) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}






#define PB191_LEN   30
#define PB191_NBLL  12
typedef enum LL2 {
  OO_0,OA_0,AO_0,AA_0,OO_1,OA_1,AO_1,AA_1,LA,LO,AL,OL
} LL2 ;


int PB191(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbChains[PB191_NBLL]  ,newNb[PB191_NBLL];
    int i ;
    nbChains[OO_0]=1; nbChains[OA_0]=1;nbChains[AO_0]=1;nbChains[AA_0]=1;
    nbChains[OO_1]=0; nbChains[OA_1]=0;nbChains[AO_1]=0;nbChains[AA_1]=0;
    nbChains[LA]=1; nbChains[AL]=1;nbChains[LO]=1;nbChains[OL]=1;
    for(i=2;i<PB191_LEN;i++) {
        newNb[OO_0] = nbChains[OO_0] + nbChains[AO_0] ;
        newNb[OO_1] = nbChains[OO_1] + nbChains[AO_1] + nbChains[LO] ;

        newNb[AA_0] = nbChains[OA_0] ;
        newNb[AA_1] = nbChains[OA_1] + nbChains[LA] ;

        newNb[OA_0] = nbChains[OO_0] + nbChains[AO_0] ;
        newNb[OA_1] = nbChains[OO_1] + nbChains[AO_1] + nbChains[LO] ;

        newNb[AO_0] = nbChains[AA_0] + nbChains[OA_0] ;
        newNb[AO_1] = nbChains[AA_1] + nbChains[OA_1] + nbChains[LA] ;
        
        newNb[LA] = nbChains[AL] + nbChains[OL];
        newNb[LO] = nbChains[AL] + nbChains[OL] ;

        newNb[AL] = nbChains[OA_0] + nbChains[AA_0];
        newNb[OL] = nbChains[OO_0] + nbChains[AO_0];

        memcpy(nbChains,newNb,sizeof(nbChains)) ;
    }
    int64_t nbTot = 0;
    for(i=0;i<PB191_NBLL;i++) {
        nbTot += nbChains[i] ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NbChains=%lld\n",pbR->ident,nbTot) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbTot) ;
    return 1 ;
}

typedef enum LL2a {
    A0_0,A0_1,A1_0,A1_1,A2_0,A2_1
} LL2a ;

#define PB191_NBLLa  6

int PB191a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbChains[PB191_NBLLa]  ,newNb[PB191_NBLLa];
    int i ;
    for(i=0;i<PB191_NBLLa;i++) nbChains[i] = 0 ;
    nbChains[A0_0]=1;
    for(i=0;i<PB191_LEN;i++) {
        newNb[A0_0] = nbChains[A0_0] + nbChains[A1_0] + nbChains[A2_0] ; //  rajout O
        newNb[A0_1] = nbChains[A0_1] + nbChains[A1_1] + nbChains[A2_1]    // rejout O
            +  nbChains[A0_0] + nbChains[A1_0] + nbChains[A2_0] ; // rajout L
        newNb[A1_0] = nbChains[A0_0] ;
        newNb[A1_1] = nbChains[A0_1] ;
        newNb[A2_0] = nbChains[A1_0] ;
        newNb[A2_1] = nbChains[A1_1] ;
        memcpy(nbChains,newNb,sizeof(nbChains)) ;
    }
    int64_t nbTot = 0;
    for(i=0;i<PB191_NBLLa;i++) {
        nbTot += nbChains[i] ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NbChains=%lld\n",pbR->ident,nbTot) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbTot) ;
    return 1 ;
}


#define PB192_MAXN  100000
#define PB192_PREC    1000000000000LL
//#define PB192_PREC  1000000000000LL
#define PB192_HALF  1000LL
//#define PB192_HALF  1000000LL
//#define PB192_PREC  100LL

#define HIGH_PREC_192   0


int PB192(PB_RESULT *pbR) {
    int32_t N , a0, a2;
    pbR->nbClock = clock()  ;
    int64_t Sum = 0 ;
    for(N=2,a0=1,a2=4;N<=PB192_MAXN;N++) {
        int32_t n , d ;
        int64_t p0,q0,p1,q1,p2,q2,pk,qk ;
        int32_t a ;
        if(N == a2) { // a2 = (a0+1)*(a0+1)
            a0++ ;
            a2 += 2*a0 + 1 ; continue ;
        }
        // compute the convergent for sqrt(N)
        // in place with 3 consecutives p0/q0 p1/q1 p2/q2
        a = a0 ; d=1 ;  n = 0 ; // so k0 =(int) srqt(N)
        p1=1 ; q1=0;
        p2=a ; q2 = 1 ;
        do {
            n = d * a - n ;
            d = (N - n*n) / d ;
            a = (a0+n) / d ;
            p0 = p1 ;  p1 = p2 ;  p2 = a*p1 + p0  ;
            q0 = q1 ; q1 = q2 ;   q2 =a*q1 + q0 ;
        } while( q2 <= PB192_PREC) ;
        // p2/q2 exceed precision. p1/q1 is the last convergent OK
        // must test if pk/qk = (p0+k*p1)/(q0+k*q1) is better
        // with k max value not ot excced precision
       int64_t k = (PB192_PREC - q0) / q1 ;
        if(k == 0) {
            pk = p1 ;
            qk = q1 ;
        } else {
            pk = p0 + k * p1 ;
            qk = q0 + k * q1 ;
            if(2*k < a) {
                pk = p1 ;
                qk = q1 ;
            } else if( 2*k==a){
                // compute the remaining convergent np1/nq1 , np2/nq2
                // and compare to q1/q0 ( q(n)/q(n-1) )
                // comparaison depend on parity as convergents alternate.
                n = d * a - n ;
                d = (N - n*n) / d ;
                a = (a0+n) / d ;
                int64_t np1 = 1 , nq1 = 0 ;
                int64_t np2 = a , nq2 = 1 ;
                int is = 1 ;
                do {
                    n = d * a - n ;
                    d = (N- n*n) / d ;
                    a = (a0+n) / d ;
                    
                    int64_t tmp = np1 ;
                    np1 = np2 ;
                    np2 = a*np2 + tmp ;
                    
                    tmp = nq1 ;
                    nq1 = nq2 ;
                    nq2 =a*nq2 + tmp ;
                    
                    is = -is ; // parity
                } while((np2*q0-nq2*q1) * is < 0 ) ; // test loop on (n,d) = (k0,1)first couple
                if(is < 0) {
                    pk = p1 ;
                    qk = q1 ;
                }
            }
        }
        Sum += qk ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
}

 #define PB193_MAX   1125899906842624LL
// #define PB193_MAX   64LL

int PB193(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    int maxP = (int) Sqrt64(PB193_MAX)+1 ;
    if((ctxP = Gen_tablePrime(maxP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbM = nbPrime ;
    const T_prime *ptnbM = tbPrime+nbM ;
    const T_prime *pt1,*pt2,*pt3,*pt4,*pt5,*pt6,*pt7,*pt8 ;
    int64_t nbSqrMult = 0 ;
    int64_t P2,P3,P4,P5,P6,P7,P8 ;
    for(pt1=tbPrime;pt1<ptnbM ;pt1++) {
        int64_t P1 = PB193_MAX / ( *pt1 * (int64_t) *pt1)  ;
        nbSqrMult += P1 ;
        for(pt2=pt1+1;pt2<ptnbM && (P2=*pt2 * (int64_t) *pt2 ) <= P1 ;pt2++) {
            P2 = P1 / P2 ;
            nbSqrMult -= P2 ;
            for(pt3=pt2+1;pt3<ptnbM && (P3=*pt3 * (int64_t) *pt3 ) <= P2 ;pt3++) {
                P3 = P2 / P3 ;
                nbSqrMult += P3 ;
                for(pt4=pt3+1;pt4<ptnbM && (P4=*pt4 * (int64_t) *pt4 ) <= P3 ;pt4++) {
                    P4 = P3 / P4 ;
                    nbSqrMult -= P4 ;
                    for(pt5=pt4+1;pt5<ptnbM && (P5=*pt5 * (int64_t) *pt5 ) <= P4 ;pt5++) {
                        P5 = P4 / P5 ;
                        nbSqrMult += P5 ;
                        for(pt6=pt5+1;pt6<ptnbM && (P6=*pt6 * (int64_t) *pt6 ) <= P5;pt6++) {
                            P6 = P5 / P6 ;
                            nbSqrMult -= P6 ;
                            for(pt7=pt6+1;pt7<ptnbM  && (P7=*pt7 * (int64_t) *pt7 ) <= P6 ;pt7++) {
                                P7 = P6 / P7 ;
                                nbSqrMult += P7 ;
                                for(pt8=pt7+1;pt8<ptnbM && (P8=*pt8 * (int64_t) *pt8 ) <= P7 ;pt8++) {
                                    P8 = P7 / P8 ;
                                    nbSqrMult -= P8 ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int64_t nbFind = PB193_MAX - nbSqrMult ;
 
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbFind) ;
    return 1 ;
}

int PB193a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    int maxP = (int) Sqrt64(PB193_MAX) ;
    if((ctxP = Gen_tablePrime(maxP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int64_t *tbPrime2 = malloc((nbPrime+1)*sizeof(tbPrime2[0])) ;
    int i;
    for(i=0;i<nbPrime;i++) {
        tbPrime2[i] = tbPrime[i]*(int64_t)tbPrime[i] ;
    }
    tbPrime2[i] = tbPrime2[i-1] ; // on duplique le derneir element pour eviter un test
    Free_tablePrime(ctxP) ;

    int64_t *ptnbM = tbPrime2+nbPrime ;
    int64_t *pt1,*pt2,*pt3,*pt4,*pt5,*pt6,*pt7,*pt8 ;
    int64_t nbSqrMult = 0 ;
    int64_t P2,P3,P4,P5,P6,P7,P8 ;
    for(pt1=tbPrime2;pt1<ptnbM ;pt1++) {
        int64_t P1 = PB193_MAX / *pt1  ;
        nbSqrMult += P1 ;
        for(pt2=pt1+1; *pt2 <= P1 ;pt2++) {
            P2 = P1 / *pt2 ;
            nbSqrMult -= P2 ;
            for(pt3=pt2+1; *pt3  <= P2 ;pt3++) {
                P3 = P2 / *pt3 ;
                nbSqrMult += P3 ;
                for(pt4=pt3+1;*pt4  <= P3 ;pt4++) {
                    P4 = P3 / *pt4 ;
                    nbSqrMult -= P4 ;
                    for(pt5=pt4+1;*pt5  <= P4 ;pt5++) {
                        P5 = P4 / *pt5 ;
                        nbSqrMult += P5 ;
                        for(pt6=pt5+1;*pt6  <= P5;pt6++) {
                            P6 = P5 / *pt6 ;
                            nbSqrMult -= P6 ;
                            for(pt7=pt6+1;*pt7  <= P6 ;pt7++) {
                                P7 = P6 / *pt7 ;
                                nbSqrMult += P7 ;
                                for(pt8=pt7+1;*pt8  <= P7 ;pt8++) {
                                    P8 = P7 / *pt8 ;
                                    nbSqrMult -= P8 ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int64_t nbFind = PB193_MAX - nbSqrMult ;
    free(tbPrime2) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbFind) ;
    return 1 ;
}

 // #define PB195_MAXR   100


// #define PB195_MAXR   10000

#define PB195_MAXR  1053779

//#define PB195_MAXR 10000000
// parametrage triangle primitif m^n (premiers)
// a=m(3m+2n) ; b = (m+n)("m+n)
// c = sqrt(a*a+b*b-axb) = 3m*m+3*mxn+n*n
// R = sqrt(3)/2 * m * (m+n) si n%3 != 0
// si n= 3*p
// a/3 , b/3 est primitif
// R = 1/(2*sqrt(3)) * m * (m+3p)

int PB195(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    double R =  PB195_MAXR * 2 / sqrtl(3.0);
    double R3 = PB195_MAXR * 2 * sqrtl(3.0) ;
    int64_t nbSol = 0 ;
    int m,mr ;
    int mMax = sqrt(R3)+1 ;
     for(m=1;m<=mMax;m++) {
         int n ;
          for(n=1; (mr=m*(n+m))<=R3 ;n++) {
             if(PGCD(m,n) > 1) continue ;
             if((n % 3) != 0) {
                 nbSol += R / mr ;
             } else {
                 nbSol += R3 / mr ;
             }
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nbsol=%lld\n",pbR->ident ,nbSol) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbSol) ;
    return 1 ;
}

int PB195a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    double R =  PB195_MAXR * 2 / sqrtl(3);
    double R3 = PB195_MAXR * 2 * sqrtl(3) ;
    int64_t nbSol = 0 ;
    int m,n,mr ;
    int mMax = sqrt(R3)+1 ;
    int nMax = R3/4 + 2 ; // for m>=4
    u_int8_t *isNotPrime_mn = calloc(nMax,sizeof(isNotPrime_mn[0])) ;
    for(n=1;(mr=n+1)<=R3;n++) { nbSol += (n % 3) ? R/mr : R3/mr ; } // m=1
    for(n=1;(mr=2*n+4)<=R3;n +=2) { nbSol += (n % 3) ? R/mr : R3/mr ; } // m=2
    for(n=1;(mr=3*n+9)<=R;n++) { if( (n%3) != 0 ) nbSol += R/mr ; } // m=3
    for(m=4;m<=mMax;m++) { // m >= 4
        nMax =R3/m - m ; // m*(m+n)=R3
        nbSol += R/(m*(m+1)) ; // n=1
        for(n=2;n<=nMax;n++){
            if(isNotPrime_mn[n]) { isNotPrime_mn[n] = 0 ; continue ; }
            if(n <= m) {
                if( (m % n) == 0 ) { // n divisor of m
                    int np ; // invalidate multiple on n
                    for(np = 2*n; np<=nMax;np+=n) isNotPrime_mn[np] = 1 ;
                    continue ;
                }
            }
            // n is prime with m
            mr = m * (m+n) ;
            nbSol += (n % 3) ? R/mr : R3/mr ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nbsol=%lld\n",pbR->ident ,nbSol) ;
    free(isNotPrime_mn) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbSol) ;
    return 1 ;
}



#define PB198_MAXQ  100000000
#define PB198_MIND  100

int PB198e(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    int d0Max = (int32_t)Sqrt64(PB198_MAXQ/2) ;
    FRACTRED fr0 ;
    nbA += (PB198_MAXQ/2 - PB198_MIND) / (PB198_MIND*PB198_MIND) ;
    int nbLoop = 0 ;
    for(fr0.d=PB198_MIND+1;fr0.d<=d0Max;fr0.d++) {
        for(fr0.n=1;PB198_MIND*fr0.n<=fr0.d;fr0.n++) {
            nbLoop++ ;
            FRACTRED fr1 =Besout(fr0);
            if(fr1.d*fr0.n-fr1.n*fr0.d != -1) continue ;
            int diff = PB198_MAXQ/2/fr0.d - fr1.d ;
            nbA += diff / fr0.d  ;
 //           {  int nb = diff / fr0.d  ; int i ;  for(i=1;i<=nb;i++) printf("%d/%d->%d/%d\n",fr0.n,fr0.d,fr1.n+i*fr0.n,fr1.d+i*fr0.d);       }
            int d = -fr1.d + fr0.d ;
               diff = PB198_MAXQ/2/fr0.d - d ;
            nbA += diff / fr0.d  ;
//          {   int nb = diff / fr0.d  ; int i ; for(i=1;i<=nb;i++) printf("%d/%d->%d/%d\n",fr0.n,fr0.d,n+i*fr0.n,d+i*fr0.d); }
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with Besout d0<%d n0/d0<%d loops=%d\n",pbR->ident,nbA,d0Max,PB198_MIND,nbLoop) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}

int PB198f(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    int d0Max = (int32_t)Sqrt64(PB198_MAXQ/2) ;
    printf("D0max=%d\n",d0Max);
    nbA += (PB198_MAXQ/2 - PB198_MIND) / (PB198_MIND*PB198_MIND) ;
    
    SBTree *sbt = SBT_alloc() ;
    FRACTRED fr0 = {0,1} ;
    FRACTRED fr1= {1,PB198_MIND} ;
    SBT_init(sbt,fr0,fr1) ;
    int nbLoop = 0 ;
    while(sbt->indS > 0  ) {
        nbLoop++ ;
        if(sbt->fr0.d + sbt->fr1.d <=d0Max ) {
            SBT_ValidNxt(sbt,1) ;
        } else {
            if(sbt->fr0.d <=d0Max && sbt->fr0.n ) {
 //               printf("%d/%d->%d/%d :",sbt->fr0.n,sbt->fr0.d,sbt->fr1.n,sbt->fr1.d);
                int d0 = sbt->fr0.d ;  int n0 = sbt->fr0.n ;  int d1 = sbt->fr1.d ;    int n1 = sbt->fr1.n ;
                if(d0<d1) {
                    int q = n1/n0 ;
                    d1 -= q * d0 ;
                    n1 -= q * n0 ;
                }
                int diff = PB198_MAXQ/2/d0 -  d1 ;
                nbA += diff / d0 ;
//              {  int nb = diff /  d0  ; int i ;  for(i=1;i<=nb;i++) printf("%d/%d ",n1+i*n0,d1+i*d0); printf(";") ; }
                d1 = -d1 + d0 ;
                n1 = -n1 + n0 ;
                diff = PB198_MAXQ/2/d0 - d1 ;
                nbA += diff / d0 ;
  //            { int nb = diff / d0 ; int i ; for(i=1;i<=nb;i++) printf("%d/%d ",n1+i*n0,d1+i*d0); printf("\n");}
            }
            SBT_ValidNxt(sbt,0) ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with stack(%d) d0<%d loops=%d\n",pbR->ident,nbA,sbt->sizeStack,d0Max,nbLoop) ;
    SBT_free(sbt);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}

#define PB198_Nend     1
#define PB198_Dend     PB198_MIND

int PB198g(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    
    int64_t nbA = 0 ;
    
    nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; // add 1/2k for k=51,52,...,99,100,... 50000000
    int N = (int32_t)Sqrt64(PB198_MAXQ/2) ;
    nbA += (PB198_MAXQ/2 - PB198_MIND) / (PB198_MIND*PB198_MIND) ; // add (1/100 + k/(1+k*100) for (100*(1+k*100) <= 50000000
 
    
    int nbLoop = 0 ;
    int n0, n, d0, d ;
    n0 = n = 1 ;
    d0 = N+1 ;  d = N ;
//    d0 = Sqrt32(PB198_MAXQ/4) + 1 ; d =d0-1 ;
    int n_end = 1 ;
    int d_end = PB198_MIND ;
    //  satisfait besout n x d0 - d * n0 = 1
    do {
        nbLoop++ ;
        int a = (N+d0)/d ; // on cherche d = a * d - d0 le plus grand possible
        int tmp = d ;
        d = a * d - d0 ;
        d0 = tmp ;
        tmp = n ;
        n = a * n - n0 ; // n = a * n - n0 ;
        n0 = tmp ;
//               printf("%d/%d->%d/%d :",n0,d0,n,d);

        {
            int d1 = d ;    int n1 = n ;
            int diff = (PB198_MAXQ/2)/d0 ;
            if(d0<d1) {
                int q = n1/n0 ;
                d1 -= q * d0 ;
                n1 -= q * n0 ;
            }
            
            if(diff >= d0+d1) nbA += (diff-d1) / d0 ;
//              {  int nb = diff /  d0  ; int i ;  for(i=1;i<=nb;i++) printf("%d/%d ",n1+i*n0,d1+i*d0); printf(";") ; }
            d1 = -d1 + d0 ;
            n1 = -n1 + n0 ;
            
            if(diff >= d0+d1) nbA += (diff-d1) / d0 ;
//            { int nb = diff / d0 ; int i ; for(i=1;i<=nb;i++) printf("%d/%d ",n1+i*n0,d1+i*d0); printf("\n");}
            
        }
 //      } while(d >= n* PB198_MIND) ;
    } while(d != d_end || n != n_end ) ;
     if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version  direct with d0<%d loops=%d\n",pbR->ident,nbA,N,nbLoop) ;
    
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}


/*
static inline int PB198CB(int d0,int d1) {
    if ( d1 <= PB198_MAXQ/2*PB198_MIND  && d0*(int64_t)d1 <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}
*/


int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += PB198_MIND/2+ (PB198_MAXQ)/2 - PB198_MAXQ/(2*PB198_MIND) ;

    
    SBdTree *sbdt = SBdT_alloc() ;
    SBdT_init(sbdt,1, PB198_MIND) ;
    int nbLoop = 0 ;
    while(sbdt->indS > 0 ) {
        nbLoop++ ;
         if(sbdt->d1 <=PB198_MAXQ/(2*PB198_MIND) && sbdt->d0*(int64_t)sbdt->d1 <= PB198_MAXQ/2) {
 //            printf("%d->%d ",sbdt->d0,sbdt->d1) ;
            nbA++ ;
            if(SBdT_ValidNxt(sbdt,1)==0) {
                if(pbR->isVerbose) fprintf(stdout,"\tPB%s ERROR REALLOC SBT(%d)\n",pbR->ident,sbdt->sizeStack) ;
                SBdT_free(sbdt);
                return 0 ;
            }
        } else {
           SBdT_ValidNxt(sbdt,0) ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with stack(%d) (den only) loops=%d\n",pbR->ident,nbA-1,sbdt->sizeStack,nbLoop) ;
    SBdT_free(sbdt);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA-1) ;
    return 1 ;
}





// Impementation of Stern-Brocot Tree

int PB198a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; //
    FRACTRED fr0 = {1,2} , fr1 = {1,1} ; ;
    
    SBTree *sbt = SBT_alloc() ;
    SBT_init(sbt,fr0, fr1) ;
    int nbLoop = 0 ;
    while(sbt->indS > 0 ) {
        nbLoop++ ;
        int64_t dd = 2 * (int64_t) sbt->fr0.n * sbt->fr1.n ;
        int k = (int)((int64_t)(Sqrt64(dd*PB198_MAXQ+1)-sbt->fr0.n * (int64_t)sbt->fr1.d - sbt->fr1.n * (int64_t)sbt->fr0.d)/dd) - PB198_MIND +2 ;
        if(k > 0) {
            nbA += k ;
            if(SBT_ValidNxt(sbt,1)==0) {
                if(pbR->isVerbose) fprintf(stdout,"\tPB%s ERROR REALLOC SBT(%d)\n",pbR->ident,sbt->sizeStack) ;
                SBT_free(sbt);
                return 0 ;
            }
        } else {
            SBT_ValidNxt(sbt,0) ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with stack(%d) (Fract[1/2,1]) loops=%d\n",pbR->ident,nbA-1,sbt->sizeStack,nbLoop) ;

    SBT_free(sbt);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}


static int loopPB198b = 0 ;
int PB198bCB(FRACTRED fr0, FRACTRED fr1) {
    loopPB198b++ ;
    int64_t dd = 2 * (int64_t) fr0.n * fr1.n ;
    int k = (int)((int64_t)(Sqrt64(dd*PB198_MAXQ+1)-fr0.n * (int64_t)fr1.d - fr1.n * (int64_t)fr0.d)/dd) - PB198_MIND +2 ;
    if(k>0) return k ;
    else return 0 ;
}

int PB198b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; //
    FRACTRED fr0 = {1,2} , fr1 = {1,1} ;
    nbA += STBrcv(fr0,fr1,PB198bCB) ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version recursive (Fract[1/2,1]) loops=%d\n",pbR->ident,nbA-1,loopPB198b) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}





int PB198c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    SBdTree *sbdt = SBdT_alloc() ;
    int nbLoop = 0;
    for(i=PB198_MIND;i<iMax;i++) {
        SBdT_init(sbdt,i,i+1) ;
        while(sbdt->indS > 0) {
            nbLoop++ ;
            if (sbdt->d0*(int64_t)sbdt->d1 <= PB198_MAXQ/2) {
                nbA++ ;
                SBdT_ValidNxt(sbdt,1) ;
            } else {
                SBdT_ValidNxt(sbdt,0) ;
            }
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version stack(%d) Den[%d %d] loops=%d\n",pbR->ident,nbA-1,sbdt->sizeStack,PB198_MIND,iMax,nbLoop) ;
    SBdT_free(sbdt);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}


static int loopPB198d = 0 ;
static inline int PB198dCB(int d0,int d) {
    loopPB198d++ ;
    if (d0*(int64_t)d <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}

int PB198d(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    for(i=PB198_MIND;i<iMax;i++) {
        nbA += STBrcvDen(i,i+1,PB198dCB) ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version recursive Den[%d %d] loops=%d\n",pbR->ident,nbA-1,PB198_MIND,iMax,loopPB198d) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
    return 1 ;
}



#define PB199_NBITER    10

typedef struct T199_CIRCLE {
    int nb ;
    int type ;
    double  q1 ;
    double  q2 ;
    double  q3 ;
} T199_CIRCLE ;

#define T199_T111 1
#define T199_T112 2
#define T199_T123 3


int PB199(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n,na ;
    int nbMax = na = 2 ; // T111 + T112
    for(n=0;n<PB199_NBITER;n++) {
        na *= 3 ;
        nbMax += na ;
    }
    double S = 1 ;
    int indTbyNiv[PB199_NBITER+1] ;
    int nbTbyNiv[PB199_NBITER+1] ;
    T199_CIRCLE *TC199 = malloc(nbMax*sizeof(TC199[0])) ;
    na = 0 ;
    indTbyNiv[0] = na ;
    TC199[na].type = T199_T111 ;
    TC199[na].nb = 1 ;
    TC199[na].q1 = TC199[na].q2 = TC199[na].q3 = 1/(2*sqrt(3.0)-3) ;
    S -= 3 * (2*sqrt(3.0)-3)*(2*sqrt(3.0)-3) ;
    na++ ;
    TC199[na].type = T199_T112 ;
    TC199[na].nb = 3 ;
    TC199[na].q1 = TC199[na].q2 =  1/(2*sqrt(3.0)-3) ;
    TC199[na].q3 = -1.0 ;
    na++ ;
    nbTbyNiv[0] = na - indTbyNiv[0] ;
    for(n=0;n<PB199_NBITER;n++) {
        int no = indTbyNiv[n] ;
        int na = no + nbTbyNiv[n] ;
        indTbyNiv[n+1] = na ;
        while(no<indTbyNiv[n+1]){
            double q ;
            if(TC199[no].type== T199_T123) {
                q = TC199[no].q1 + TC199[no].q2 + TC199[no].q3 + 2*sqrt( TC199[no].q1 * TC199[no].q2 + TC199[no].q1 * TC199[no].q3 + TC199[no].q2 * TC199[no].q3) ;
                TC199[na].type = T199_T123 ;
                TC199[na].nb = TC199[no].nb ;
                TC199[na].q1  = TC199[no].q1 ;
                TC199[na].q2  = TC199[no].q2 ;
                TC199[na].q3 = q ;
                na++ ;

                TC199[na].type = T199_T123 ;
                TC199[na].nb = TC199[no].nb ;
                TC199[na].q1  = TC199[no].q1 ;
                TC199[na].q2  = TC199[no].q3 ;
                TC199[na].q3 = q ;
                na++ ;

                TC199[na].type = T199_T123 ;
                TC199[na].nb = TC199[no].nb ;
                TC199[na].q1  = TC199[no].q2 ;
                TC199[na].q2  = TC199[no].q3 ;
                TC199[na].q3 = q ;
                na++ ;

                
            } else if(TC199[no].type== T199_T112) {
                q = 2*TC199[no].q1 + TC199[no].q3 + 2*sqrt( TC199[no].q1 * (TC199[no].q1+2*TC199[no].q3) ) ;
                TC199[na].type = T199_T112 ;
                TC199[na].nb = TC199[no].nb ;
                TC199[na].q1 = TC199[na].q2 = TC199[no].q1 ;
                TC199[na].q3 = q ;
                na++ ;

                TC199[na].type = T199_T123 ;
                TC199[na].nb = 2*TC199[no].nb ;
                TC199[na].q1  = TC199[no].q1 ;
                TC199[na].q2  = TC199[no].q3 ;
                TC199[na].q3 = q ;
                na++ ;

            } else { // T199_T111
                q = TC199[no].q1 * (3 + 2 *sqrt(3)) ;
                 TC199[na].type = T199_T112 ;
                TC199[na].nb = 3 * TC199[no].nb ;
                TC199[na].q1 = TC199[na].q2 = TC199[no].q1 ;
                TC199[na].q3 = q ;
                na++ ;
            }
            S -= TC199[no++].nb / (q*q) ;
       }
        nbTbyNiv[n+1] = na - indTbyNiv[n+1] ;
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s %d -> S=%.8f\n",pbR->ident,n,S) ;
     }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.8f",S) ;
    return 1 ;
}






