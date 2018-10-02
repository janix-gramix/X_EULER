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


#include "PB151_200.h"

#define PB198_MAXQ  100000000
#define PB198_N     1
#define PB198_D     9

#define PB198_Nend     1
#define PB198_Dend     100

typedef struct FRACTRED {
    int n ;
    int d ;
} FRACTRED ;
int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; // ajout de 1/2k pour k=51,52,...,99,100
    // il faut rajouter 2 fois toutes les fractions p/q < 1/100 , irreductibles et q <= 10**8
    
    
    int32_t N = PB198_MAXQ / 200 ;
    FRACTRED * fr =malloc((N*(int64_t)(N+1))/200*sizeof(fr[0])) ;
    int64_t nb = 0 ;
    int d ;
    int n ;
    int d_end=PB198_Dend ;
    int n_end=PB198_Nend ;
    // on va chercher d0 et n0 tel que
    // on ait besout n x d0 - d * n0 = 1
    int d0, n0 ;
    d=N ;
    n=1 ;
    d0=N+1 ;
    n0=1 ;
//    fr[0].n=0 ;
//    fr[0].d=1 ;
//    nb=1 ;
//    printf("%d/%d ",n,d);
//    nbA += 2 ;
    //    int d0=4 , n0 = 1 ; // satisfait besout n x d0 - d * n0 = 1
    do {
        int a = (N+d0)/d ; // on cherche d = a * d - d0 le plus grand possible
        int tmp = d ;
        d = a * d - d0 ;
        d0 = tmp ;
        tmp = n ;
        n = a * n - n0 ; // n = a * n - n0 ;
        n0 = tmp ;
        fr[nb].n= n;
        fr[nb++].d=d ;
//        printf("%d/%d ",n,d);
        // on garde le couple (n,d) comme (n0,d0) car
        // besout  n x d0 - d * n0 = 1 est toujours satisfait
        // (a *n - n0) *d - (a*d - d0) * n = n x d0 - d * n0 = 1;
    } while(d != d_end || n != n_end ) ;
    int i ;
    FRACTRED maxFr = {0,1} ;
    for(i=0;i<nb-1;i++) {
        d = fr[i].d ;
        int j ;
        if(fr[i].n > 1)
        {
            for(j=i-1;j>=0 &&  fr[j].d>d;j--) ;
            if(j>=0 && 2*(int64_t)fr[j].d*d <= PB198_MAXQ) {
//                if(d>maxFr.d)   { maxFr = fr[i] ; printf("%d/%d<-%d/%d ",fr[j].n,fr[j].d,fr[i].n,fr[i].d) ; }
                nbA++ ;
            }
        }

        for(j=i+1;fr[j].d>d;j++) ;
        if(2*(int64_t)fr[j].d*d <= PB198_MAXQ) {
    // if(d>maxFr.d) */  { maxFr = fr[i] ; printf("%d/%d->%d/%d ",fr[i].n,fr[i].d,fr[j].n,fr[j].d) ; }
            nbA++ ;
        }
    }
    

    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}

#define d0Max 7072
static int64_t P(int n1, int d1, int n0, int d0) ;
static int64_t P(int n1, int d1, int n0, int d0)
{
    int64_t all = 0 ;
    int b,c ;
    if (d0 > d0Max) return 0;
    
    int64_t e = 2*(int64_t)d0*d1;
    if (PB198_MAXQ >= e) {
        all += (PB198_MAXQ-e)/(2*d0*d0);
    } else {
        all = 0;
    }
    
    b = (PB198_MAXQ-d1)/d0;
    for (c=1; c<=b; c++) {
        all += P(n0, d0, c*n0+n1, c*d0+d1);
    }
    return all;
}


int PB198a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbA = 0 ;
    nbA += (PB198_MAXQ-100)/2 ; //
    int i ;
    for(i=100;i<PB198_MAXQ+100;i++) {
        nbA += P(0, 1, 1, i);
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbA) ;
    return 1 ;
}

