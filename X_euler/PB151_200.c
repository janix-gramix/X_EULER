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

#include "faray_utils.h"
#include "PB151_200.h"

#define PB198_MAXQ  100000000
#define PB198_N     1
#define PB198_D     9

#define PB198_Nend     1
#define PB198_Dend     100


int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    // il faut rajouter 2 fois toutes les fractions p/q < 1/100 , irreductibles et q <= 10**8
    
    
    int32_t N = PB198_MAXQ / 200 ;
    FRACTRED * fr =malloc((N*(int64_t)(N+1))/200*sizeof(fr[0])) ;
//   int * den =malloc((N*(int64_t)(N+1))/200*sizeof(den[0])) ;
    int64_t nb = 0 ;
    int d ;
    int n ;
    int d_end=PB198_Dend ;
    int n_end=PB198_Nend ;
    // on va chercher d0 et n0 tel que
    // on ait besout n x d0 - d * n0 = 1
    int d0, n0 ;
    d=Sqrt32(PB198_MAXQ/2)+1 ;
    n=1 ;
    d0=d+1 ;
    n0=1 ;
    do { // voir PB073 pour l'algorithme
        int a = (N+d0)/d ; // on cherche d = a * d - d0 le plus grand possible
        int tmp = d ;
        d = a * d - d0 ;
        d0 = tmp ;
        tmp = n ;
        n = a * n - n0 ; // n = a * n - n0 ;
        n0 = tmp ;
        fr[nb].n= n;
        fr[nb++].d=d ;
    } while(d != d_end || n != n_end ) ;
    printf("\nEND Farey\n");
    int i ;
    for(i=0;i<nb-1;i++) {
        d = fr[i].d ;
         int j ;
        int64_t dMin ;
        for(j=i+1, dMin = fr[j].d ;fr[j].d>d;j++) {
            if(j==i+1) {
                if(dMin*d <= PB198_MAXQ/2) nbA++ ;
            } else if( fr[j].d < dMin ) {
                dMin = fr[j].d ;
                if(dMin*d <= PB198_MAXQ/2) nbA++ ;
            }
        }
        if((int64_t)fr[j].d*d <= PB198_MAXQ/2) {
            nbA++ ;
        }
    }
    

    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}

/*

int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    FRACTRED fr0 = {0,1} , fr1 = {1,100} ; ;

    
    SBTree *sbt = SBT_alloc() ;
    SBT_init(sbt,fr0, fr1) ;
    while(sbt->indS > 0 ) {
        nbA++ ;
        if((int64_t)fr0.d*fr1.d <= PB198_MAXQ/2) {
            if(SBT_ValidNxt(sbt,1)==0) {
                if(pbR->isVerbose) fprintf(stdout,"\tPB%s ERROR REALLOC SBT(%d)\n",pbR->ident,sbt->sizeStack) ;
                SBT_free(sbt);
                return 0 ;
            }
        } else {
            SBT_ValidNxt(sbt,0) ;
        }
    }
    SBT_free(sbt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}
*/




// Impementation of Stern-Brocot Tree

int PB198a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    FRACTRED fr0 = {1,2} , fr1 = {1,1} ; ;
    
    SBTree *sbt = SBT_alloc() ;
    SBT_init(sbt,fr0, fr1) ;
    while(sbt->indS > 0 ) {
        int64_t dd = 2 * (int64_t) sbt->fr0.n * sbt->fr1.n ;
        int k = (int)((int64_t)(Sqrt64(dd*PB198_MAXQ+1)-sbt->fr0.n * (int64_t)sbt->fr1.d - sbt->fr1.n * (int64_t)sbt->fr0.d)/dd) - 100 +2 ;
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
    SBT_free(sbt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}



int PB198bCB(FRACTRED fr0, FRACTRED fr1) {
    int64_t dd = 2 * (int64_t) fr0.n * fr1.n ;
    int k = (int)((int64_t)(Sqrt64(dd*PB198_MAXQ+1)-fr0.n * (int64_t)fr1.d - fr1.n * (int64_t)fr0.d)/dd) - 100 +2 ;
    if(k>0) return k ;
    else return 0 ;
}

int PB198b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    FRACTRED fr0 = {1,2} , fr1 = {1,1} ;
    nbA += STBrcv(fr0,fr1,PB198bCB) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}
static inline int PB198dCB(int d0,int d) {
    if (d0*(int64_t)d <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}


int PB198c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    SBdTree *sbdt = SBdT_alloc() ;
    for(i=100;i<iMax;i++) {
        SBdT_init(sbdt,i,i+1) ;
        while(sbdt->indS > 0) {
            if (sbdt->d0*(int64_t)sbdt->d1 <= PB198_MAXQ/2) {
                nbA++ ;
                SBdT_ValidNxt(sbdt,1) ;
            } else {
                SBdT_ValidNxt(sbdt,0) ;
            }
        }
    }
    
    SBdT_free(sbdt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}



int PB198d(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    for(i=100;i<iMax;i++) {
        nbA += STBrcvDen(i,i+1,PB198dCB) ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}










