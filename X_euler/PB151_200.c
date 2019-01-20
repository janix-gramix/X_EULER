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

#define PB179_NB    10000000
//#define PB179_NB    1000
int PB179(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    uint16_t *nbProd=malloc((N+2)*sizeof(nbProd[0])) ;
    int n ;
    nbProd[2] = 2 ;
    for(n=4;n<=N;n+=2) {
        nbProd[n] = 4 ;
        nbProd[n+1] = 2 ;
    }
    int nMax = Sqrt32(N);
    for(n=3;n<=nMax;n++) {
        int np = n*n ;
        nbProd[np]++ ;
        for(np += n ;np<=N;np+=n) {
            nbProd[np] += 2 ;
            
        }
    }
    int S = 0 ;
    for(n=2;n<N;n++) {
//        printf("(%d,%d)",n,nbProd[n]) ;
        if(nbProd[n]==nbProd[n+1] ) {
            S++ ;
        }
        
    }
    free(nbProd) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}

int PB179c(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    int dMax = Sqrt32(N) ;
    int nMax = dMax+1;
    nMax += nMax & 1 ; // nMax even
    uint16_t *nbProd=malloc((nMax)*sizeof(nbProd[0])) ;
    int32_t *nxtFactor = malloc((dMax+1)*sizeof(nxtFactor[0])) ;
    int n,d ;
    for(d=3;d<=dMax;d++) {
        nxtFactor[d] = d ;
    }
    int offset = 4 ;
//    nbProd[2-offset] = 2 ;
    int nbProdAnt = 2 ; // n=3
    int S = 1 ; // (2,3)
    for(;offset < N;offset += nMax ) {
         for(n=0;n<nMax;n+=2) {
            nbProd[n] = 4 ;
            nbProd[n+1] = 2 ;
        }
        int ndMax = nMax ;
        if(offset+nMax >N+1) {
            ndMax = N+1 - offset ;
        }
        for(d=3;d<=dMax;d++) {
           int fd = nxtFactor[d] ;
           int np = d*fd - offset ;
           if(np >= ndMax) continue ;
           if(fd == d) { nbProd[np]++ ; fd++ ; np += d ; }
           for(;np<ndMax;np += d ) {
               nbProd[np] += 2 ; fd++ ;
           }
           nxtFactor[d] = fd ;
        }
 //       printf("\n [%d->%d[ ",offset,offset+ndMax) ;
        for(n=0;n< ndMax ; n++) {
 //           printf("(%d,%d)",n+offset,nbProd[n]) ;
            if(nbProd[n]==nbProdAnt) S++ ;
            nbProdAnt = nbProd[n] ;
        }
    }
    free(nbProd) ; free(nxtFactor) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}


int PB179a(PB_RESULT *pbR) {
    int N =  PB179_NB ;
    pbR->nbClock = clock() ;
    uint16_t *nbDiv=malloc((N+1)*sizeof(nbDiv[0])) ;
    int n ;
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
            int powp;
            int exp1 ;
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
//        printf("(%d,%d)",n,nbProd[n]) ;
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
    for(n=2;n<=N;n++) {
        pDivMax[n] = (n&1) ? 0 : 2 ;
    }

    int p ;
    int maxP = Sqrt32(N) ;
    
    for(p=3;p<=maxP;p+=2) {
        if(pDivMax[p]) continue ;
        for(n=p;n<=N;n+=p) {
            pDivMax[n] = p ;
        }
    }
    for(n=3;n<=N;n++) {
        if(n == pDivMax[n] || (pDivMax[n] == 0) ) {
            pDivMax[n] = 2 ;
        } else {
            p = pDivMax[n] ;
            int d,d1,exp1 ;
            for(d=n/p,exp1=2;d1=d/p,d==d1*p;d = d1) {
                exp1++ ;
            }
            pDivMax[n] = pDivMax[d]*exp1 ;
        }
    }
    int S = 0 ;
    for(n=2;n<N;n++) {
//                printf("(%d,%d)",n,pDivMax[n]) ;
        if(pDivMax[n]==pDivMax[n+1] ) {
            S++ ;
        }
        
    }
    free(pDivMax) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
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
    int i;
    int P ;
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
            int n = -fr1.n + fr0.n ;
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
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",nbA) ;
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






