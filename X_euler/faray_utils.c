//
//  faray_utils.c
//  X_euler
//
//  Created by Jeannot on 03/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "faray_utils.h"



static FRACTRED SignBesout(FRACTRED fr1,int sign) { // solve besout
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
    if((fr2.n*fr1.d-fr2.d*fr1.n)*sign < 0) { // ochange sign
        int q = fr2.n/fr1.n+1 ;
        fr2.n = -fr2.n + q*fr1.n;
        fr2.d = -fr2.d + q*fr1.d ;
    }
    return fr2 ;
}

FRACTRED Besout(FRACTRED fr1) {
    return SignBesout(fr1,1) ;
}
FRACTRED IBesout(FRACTRED fr1) {
    return SignBesout(fr1,-1) ;
}


void FRC_init(FRACTrec *FRrec,int maxDen , FRACTRED frStart, FRACTRED frEnd) {
    FRrec->maxDen = maxDen ;
    FRrec->frEnd = frEnd ;
    if(frStart.n) { // pas 0
        FRrec->fr1 = frStart ;
        FRrec->fr0 = IBesout(frStart) ;
    } else  {
        FRrec->fr0 = frStart ; // fraction null
        FRrec->fr1 = (FRACTRED){1,maxDen} ;
    }
}
// obtient (dans FRrec->fr1) la fraction suivante entre frStart et frEnd
// retourne 0 quand le parcours est termine
int FRC_getNext(FRACTrec *FRrec) {
    
    int a = (FRrec->maxDen+FRrec->fr0.d)/FRrec->fr1.d ; // searcd d = a * d - d0 biggest
    int tmp = FRrec->fr1.d ;
    FRrec->fr1.d = a * FRrec->fr1.d - FRrec->fr0.d ;
    FRrec->fr0.d = tmp ;
    tmp = FRrec->fr1.n ;
    FRrec->fr1.n = a * FRrec->fr1.n - FRrec->fr0.n ; // n = a * n - n0 ;
    FRrec->fr0.n = tmp ;
    if( FRrec->fr1.d == FRrec->frEnd.d && FRrec->fr1.n == FRrec->frEnd.n) return 0 ;
    else return 1 ;

}




void SBT_init(SBTree * sbt,FRACTRED fr0, FRACTRED fr1) {
    sbt->fr0 = fr0 ;
    sbt->indS = 0 ;
    sbt->stack[sbt->indS++] = sbt->fr1 = fr1 ;
}
#define SBT_SIZE    16384
SBTree * SBT_alloc() {
    SBTree *sbt = calloc(1,sizeof(sbt[0]));
    sbt->sizeStack = SBT_SIZE ;
    sbt->stack = malloc(sbt->sizeStack*sizeof(sbt->stack[0]));
    //    SBT_init(sbt,fr0,fr1) ;
    return sbt ;
}



int SBT_ValidNxt(SBTree * sbt, int isOK) {
    if(isOK ) {
        if(sbt->indS == sbt->sizeStack) {
            FRACTRED * nstack = realloc(sbt->stack,sbt->sizeStack*2*sizeof(nstack[0])) ;
            if(nstack) {
                sbt->sizeStack *= 2 ;
                sbt->stack = nstack ;
            } else {
                return 0 ;
            }
        }
        sbt->fr1.d += sbt->fr0.d ;
        sbt->fr1.n += sbt->fr0.n ;
        sbt->stack[sbt->indS++] = sbt->fr1 ;
        return sbt->indS ;
    } else {
        sbt->fr0 = sbt->fr1 ;
        sbt->indS-- ;
        if(sbt->indS>0) sbt->fr1 = sbt->stack[sbt->indS-1] ;
        return(sbt->indS);
    }
}

int SBT_getNext(SBTree * sbt,int maxDen) {
    while(sbt->indS > 0) {
        if(sbt->fr0.d + sbt->fr1.d <= maxDen) {
            return SBT_ValidNxt(sbt,1) ;
        } else {
            SBT_ValidNxt(sbt,0) ;
        }
    }
    return sbt->indS ;
}


SBTree * SBT_free(SBTree * sbt) {
    if(sbt) {
        free(sbt->stack);
        free(sbt);
    }
    return NULL ;
}





void SBdT_init(SBdTree * sbdt,int32_t d0, int32_t d1) {
    sbdt->d0 = d0 ;
    sbdt->indS = 0 ;
    sbdt->stack[sbdt->indS++] = sbdt->d1 = d1 ;
}
SBdTree * SBdT_alloc() {
    SBdTree *sbdt = calloc(1,sizeof(sbdt[0]));
    sbdt->sizeStack = SBT_SIZE ;
    sbdt->stack = malloc(sbdt->sizeStack*sizeof(sbdt->stack[0]));
    //    SBT_init(sbt,fr0,fr1) ;
    return sbdt ;
}


int SBdT_ValidNxt(SBdTree * sbdt, int isOK) {
    if(isOK ) {
        if(sbdt->indS == sbdt->sizeStack) {
            int32_t * nstack = realloc(sbdt->stack,sbdt->sizeStack*2*sizeof(nstack[0])) ;
            if(nstack) {
                sbdt->sizeStack *= 2 ;
                sbdt->stack = nstack ;
            } else {
                return 0 ;
            }
        }
        sbdt->d1 += sbdt->d0 ;
        sbdt->stack[sbdt->indS++] = sbdt->d1 ;
        return sbdt->indS ;
    } else {
        sbdt->d0 = sbdt->d1 ;
        sbdt->indS-- ;
        if(sbdt->indS>0) sbdt->d1 = sbdt->stack[sbdt->indS-1] ;
        return(sbdt->indS);
    }
}
SBdTree * SBdT_free(SBdTree * sbdt) {
    if(sbdt) {
        free(sbdt->stack);
        free(sbdt);
    }
    return NULL ;
}



int STBrcv(FRACTRED fr0, FRACTRED fr1,STB_CB stb_CB )
{
    int sum=0;
    int k ;
    if ((k=stb_CB(fr0,fr1)) != 0)
    {
        FRACTRED frs ;
        frs.d = fr1.d + fr0.d ; frs.n = fr0.n + fr1.n ;
        sum+=k+STBrcv(fr0, frs,stb_CB);
        sum+=STBrcv(frs, fr1,stb_CB);
    }
    return sum ;
}






int STBrcvDen(int d0, int d1,STBDen_CB stbDen_CB )
{
    int sum=0;
    int k ;
    if ((k=stbDen_CB(d0,d1)) != 0)
    {
        sum+=k+STBrcvDen(d0, d0+d1,stbDen_CB);
        sum+=STBrcvDen(d0+d1, d1,stbDen_CB);
    }
    return sum ;
}


