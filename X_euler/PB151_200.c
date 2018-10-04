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
//#define PB198_MAXQ  1000000
#define PB198_N     1
#define PB198_D     9

#define PB198_Nend     1
#define PB198_Dend     100

/*
int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    // il faut rajouter 2 fois toutes les fractions p/q < 1/100 , irreductibles et q <= 10**8
   int32_t N = PB198_MAXQ / 200 ;
   int * den =malloc((N*(int64_t)(N+1))/200*sizeof(den[0])) ;
    int32_t nb = 0 ;
    int d_end=PB198_Dend ;
    int n_end=PB198_Nend ;
    // on va chercher d0 et n0 tel que
    // on ait besout n x d0 - d * n0 = 1
    int d=Sqrt32(PB198_MAXQ/2)+1 ;
    int d0=d+1 ;
    int n=1 ;
    int n0=1 ;
    do { // voir PB073 pour l'algorithme
        int a = (N+d0)/d ; // on cherche d = a * d - d0 le plus grand possible
        int tmp = d ;
        d = a * d - d0 ;
        d0 = tmp ;
        tmp = n ;
        n = a * n - n0 ; // n = a * n - n0 ;
        n0 = tmp ;
        den[nb++] = d ;
    } while(d != d_end || n != n_end ) ;
    printf("\nEND Farey(%d)\n",nb);
    int i ;
    for(i=0;i<nb-1;i++) {
        d = den[i] ;
         int j ;
        int64_t dMin ;
        int isMin = 0 ;
        for(j=i+1, dMin = den[j] ;den[j]>d;j++) {
            if(j==i+1) {
                if(dMin*d <= PB198_MAXQ/2) { nbA++ ; isMin = 1 ; }
            } else if( den[j] < dMin ) {
                dMin = den[j] ;
                if(isMin) nbA++;
                else if(dMin*d <= PB198_MAXQ/2) { nbA++ ; isMin=1 ; }
            }
        }
        if((int64_t)den[j]*d <= PB198_MAXQ/2) {
            nbA++ ;
        }
    }

    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}
*/
FRACTRED Besout(FRACTRED fr1) { // solve besout
    int s0 = 1, s1 = 0;
    int t0 = 0, t1 = -1 ;
    int n1 = fr1.n ;
    int d1 = fr1.d ;
    do {
        int q = d1 / n1 ;
        int tmp = d1 - q * n1 ;
        d1 = n1 ;  n1 = tmp ;
        
        tmp = s0 + q * s1 ;
        s0 = s1 ; s1 = tmp ;
        
        tmp = t0  + q * t1 ;
        t0 = t1 ; t1 = tmp ;
        
    } while ( n1 ) ;
    FRACTRED fr2 ;
    fr2.d = -t0 ;
    fr2.n = s0 ;
    if(fr1.d*fr2.n-fr1.n*fr2.d == -1) { // on inverse le signe
        int q = fr2.n/fr1.n+1 ;
        fr2.n = -fr2.n + q*fr1.n;
        fr2.d = -fr2.d + q*fr1.d ;
    }
    return fr2 ;
}
/*

int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    // il faut rajouter 2 fois toutes les fractions p/q < 1/100 , irreductibles et q <= 10**8
    int32_t N = PB198_MAXQ / 200 ;
//    int * den =malloc((N*(int64_t)(N+1))/200*sizeof(den[0])) ;
    int32_t nb = 0 ;
    int d_end=PB198_Dend ;
    int n_end=PB198_Nend ;
    // on va chercher d0 et n0 tel que
    // on ait besout n x d0 - d * n0 = 1
    int d=Sqrt32(PB198_MAXQ/2) ;
    int d0=d+1 ;
    int n=1 ;
    int n0=1 ;
    do { // voir PB073 pour l'algorithme
 //       int a = (N+d0)/d ; // on cherche d = a * d - d0 with
        printf("%d/%d->%d/%d ",n0,d0,n,d) ;
        int a = (PB198_MAXQ/2 + d0* (int64_t) d)/(d*(int64_t)d) ;
        if(a*d-d0 <= 0) {
            a= d0/d+1 ;
            d0 = a * d - d0 ;
            n0 = a * n - n0 ;
            FRACTRED fr0, fr1 ;
            fr0.d =d0 ;
            fr0.n = n0 ;
            fr1= Besout(fr0) ;
            d = fr1.d ;
            n = fr1.n ;
            
            a = (PB198_MAXQ/2 + d0* (int64_t) d)/(d0*(int64_t)d0) ;
            int tmp = d ;
            d = d + a * d0 ;
            d0 = tmp ;
            tmp = n ;
            n = n + a * n0 ;
            n0 = tmp ;

            
            
        } else {
            int tmp = d ;
            d = a * d - d0 ;
            if(d*(int64_t)tmp <= PB198_MAXQ/2) nbA++ ;
            d0 = tmp ;
            tmp = n ;
            n = a * n - n0 ; // n = a * n - n0 ;
            n0 = tmp ;
        }
 //       den[nb++] = d ;
    } while(d != d_end || n != n_end ) ;
    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}
*/

static inline int PB198CB(int d0,int d1) {
    if ( d1 <= PB198_MAXQ/200  && d0*(int64_t)d1 <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}



int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += 50+ (PB198_MAXQ)/2 - PB198_MAXQ/200 ;

    
    SBdTree *sbdt = SBdT_alloc() ;
    SBdT_init(sbdt,1, 100) ;
    int nbLoop = 0 ;
    while(sbdt->indS > 0 ) {
        nbLoop++ ;
         if(sbdt->d1 <=PB198_MAXQ/200 && sbdt->d0*(int64_t)sbdt->d1 <= PB198_MAXQ/2) {
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
    sprintf(pbR->strRes,"%lld",nbA-1) ;
    return 1 ;
}





// Impementation of Stern-Brocot Tree

int PB198a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    FRACTRED fr0 = {1,2} , fr1 = {1,1} ; ;
    
    SBTree *sbt = SBT_alloc() ;
    SBT_init(sbt,fr0, fr1) ;
    int nbLoop = 0 ;
    while(sbt->indS > 0 ) {
        nbLoop++ ;
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
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with stack(%d) (Fract[1/2,1]) loops=%d\n",pbR->ident,nbA-1,sbt->sizeStack,nbLoop) ;

    SBT_free(sbt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}


static int loopPB198b = 0 ;
int PB198bCB(FRACTRED fr0, FRACTRED fr1) {
    loopPB198b++ ;
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
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version recursive (Fract[1/2,1]) loops=%d\n",pbR->ident,nbA-1,loopPB198b) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}



int PB198c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    SBdTree *sbdt = SBdT_alloc() ;
    int nbLoop = 0;
    for(i=100;i<iMax;i++) {
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
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version stack(%d) Den[100 %d] loops=%d\n",pbR->ident,nbA-1,sbdt->sizeStack,iMax,nbLoop) ;
    SBdT_free(sbdt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
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
    nbA += (PB198_MAXQ-100)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    for(i=100;i<iMax;i++) {
        nbA += STBrcvDen(i,i+1,PB198dCB) ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version recursive Den[100 %d] loops=%d\n",pbR->ident,nbA-1,iMax,loopPB198d) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}










