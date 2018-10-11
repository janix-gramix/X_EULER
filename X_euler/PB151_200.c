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
#include <gmp.h>

#include "faray_utils.h"
#include "PB151_200.h"


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
    sprintf(pbR->strRes,"%lld",Sum);
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
            int n = -fr1.n + fr0.n ;
            diff = PB198_MAXQ/2/fr0.d - d ;
            nbA += diff / fr0.d  ;
//          {   int nb = diff / fr0.d  ; int i ; for(i=1;i<=nb;i++) printf("%d/%d->%d/%d\n",fr0.n,fr0.d,n+i*fr0.n,d+i*fr0.d); }
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version with Besout d0<%d n0/d0<%d loops=%d\n",pbR->ident,nbA,d0Max,PB198_MIND,nbLoop) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
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
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}

#define PB198_Nend     1
#define PB198_Dend     PB198_MIND

int PB198g(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    
    int64_t nbA = 0 ;
    int bcl ;
    
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
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}



static inline int PB198CB(int d0,int d1) {
    if ( d1 <= PB198_MAXQ/2*PB198_MIND  && d0*(int64_t)d1 <= PB198_MAXQ/2) return 1 ;
    else return 0 ;
}



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
    sprintf(pbR->strRes,"%lld",nbA-1) ;
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
    sprintf(pbR->strRes,"%lld",nbA) ;
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
    sprintf(pbR->strRes,"%lld",nbA) ;
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
    nbA += (PB198_MAXQ-PB198_MIND)/2 ; //
    int i ;
    int iMax= Sqrt32(PB198_MAXQ) ;
    for(i=PB198_MIND;i<iMax;i++) {
        nbA += STBrcvDen(i,i+1,PB198dCB) ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s S=%lld,Version recursive Den[%d %d] loops=%d\n",pbR->ident,nbA-1,PB198_MIND,iMax,loopPB198d) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}










