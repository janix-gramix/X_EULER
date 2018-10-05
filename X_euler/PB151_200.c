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


int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    int d0Max = Sqrt32(PB198_MAXQ/2) ;
    FRACTRED fr0 ;
    nbA += (PB198_MAXQ/2 - 100) / (100*100) ;
    for(fr0.d=101;fr0.d<=d0Max;fr0.d++) {
        for(fr0.n=1;100*fr0.n<=fr0.d;fr0.n++) {
            FRACTRED fr1 =Besout(fr0);
            if(fr1.d*fr0.n-fr1.n*fr0.d != -1) continue ;
            int diff = PB198_MAXQ/2 - fr0.d * fr1.d ;
            if(diff >0) {
                int nb = diff / (fr0.d * fr0.d) ;
                int i ;
//                for(i=1;i<=nb;i++) printf("%d/%d->%d/%d\n",fr0.n,fr0.d,fr1.n+i*fr0.n,fr1.d+i*fr0.d);
                nbA += nb  ;
            }
            int d = -fr1.d + fr0.d ;
            int n = -fr1.n + fr0.n ;
            diff = PB198_MAXQ/2 - fr0.d * d ;

            if(diff > 0) {
                int nb = diff / (fr0.d * fr0.d) ;
                int i ;
//                for(i=1;i<=nb;i++) printf("%d/%d->%d/%d\n",fr0.n,fr0.d,n+i*fr0.n,d+i*fr0.d);
                nbA += nb  ;
            }
            
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}



static inline int PB198CB(int d0,int d1) {
    if ( d1 <= PB198_MAXQ/200  && d0*(int64_t)d1 <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}


/*
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
*/




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
    SBTree *sbt = SBT_alloc() ;
    int nbLoop = 0;
    for(i=100;i<iMax;i++) {
        FRACTRED fr0 = {1,i} ;
        FRACTRED fr1= {1,i+1} ;
        SBT_init(sbt,fr0,fr1) ;
        while(sbt->indS > 0) {
            nbLoop++ ;
            if (sbt->fr0.d*(int64_t)sbt->fr1.d <= PB198_MAXQ/2) {
//                printf("%d/%d->%d/%d\n",sbt->fr0.n,sbt->fr0.d,sbt->fr1.n,sbt->fr1.d);
                nbA++ ;
                SBT_ValidNxt(sbt,1) ;
            } else {
                SBT_ValidNxt(sbt,0) ;
            }
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version stack(%d) Den[100 %d] loops=%d\n",pbR->ident,nbA-1,sbt->sizeStack,iMax,nbLoop) ;
    SBT_free(sbt);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}
/*

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
*/

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










