//
//  PB201_250.c
//  X_euler
//
//  Created by Jeannot on 14/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "PB201_250.h"

#define PB204_MAXP  100
#define PB204_MAXN  2000000000

int PB204(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int32_t tbPrime[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    int nbPrime = sizeof(tbPrime)/sizeof(tbPrime[0]) ;
    int32_t *ptnbM = tbPrime+nbPrime ;
    int32_t *pt1,*pt2,*pt3,*pt4,*pt5,*pt6,*pt7,*pt8,*pt9 ;
    int32_t nbHamm = 1 ; // for 1
    int32_t p1,p2,p3,p4,p5,p6,p7,p8,p9 ;
    int32_t P1,P2,P3,P4,P5,P6,P7,P8,P9 ;
    for(pt1=tbPrime;pt1<ptnbM ;pt1++) {
        p1 = *pt1 ; P1 = PB204_MAXN / p1 ;
        for(P1 = PB204_MAXN / p1 ;P1>=1 ; P1 /= p1) {
          nbHamm++ ;
            for(pt2=pt1+1;pt2<ptnbM && (p2=*pt2) <= P1 ;pt2++) {
              for(P2 = P1/ p2 ; P2>=1 ; P2 /= p2) {
                nbHamm++ ;
                for(pt3=pt2+1;pt3<ptnbM && (p3=*pt3) <= P2 ;pt3++) {
                  for(P3 = P2/ p3 ; P3>=1 ; P3 /= p3) {
                    nbHamm++ ;
                    for(pt4=pt3+1;pt4<ptnbM && (p4=*pt4) <= P3 ;pt4++) {
                      for(P4 = P3/ p4 ; P4>=1 ; P4 /= p4) {
                        nbHamm++ ;
                        for(pt5=pt4+1;pt5<ptnbM && (p5=*pt5) <= P4 ;pt5++) {
                          for(P5 = P4/ p5 ; P5>=1 ; P5 /= p5) {
                            nbHamm++ ;
                            for(pt6=pt5+1;pt6<ptnbM && (p6=*pt6) <= P5 ;pt6++) {
                              for(P6 = P5/p6 ; P6>=1 ; P6 /= p6) {
                                nbHamm++ ;
                                for(pt7=pt6+1;pt7<ptnbM && (p7=*pt7) <= P6 ;pt7++) {
                                  for(P7 = P6/p7 ; P7>=1 ; P7 /= p7) {
                                    nbHamm++ ;
                                    for(pt8=pt7+1;pt8<ptnbM && (p8=*pt8) <= P7 ;pt8++) {
                                      for(P8 = P7/p8 ; P8>=1 ; P8 /= p8) {
                                        nbHamm++ ;
                                        for(pt9=pt8+1;pt9<ptnbM && (p9=*pt9) <= P8 ;pt9++) {
                                          for(P9 = P8/p9 ; P9>=1 ; P9 /= p9) {
                                            nbHamm++ ;
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
        }
    }
    
    
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s \n",pbR->ident) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbHamm);
    return 1 ;
 }

int CmpHam(const void *el1,const void *el2) {
    return ((int *)el1)[0] - ((int *)el2)[0] ;
}
int PB204b(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int32_t tbPrime[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    int nbPrime = sizeof(tbPrime)/sizeof(tbPrime[0]) ;
    int32_t *ptnbM = tbPrime+nbPrime ;
    int iP[20] ;
    int32_t P[20] ;
    int32_t Px[20] ;
    int32_t nbHamm = 1 ; // for 1
    int is,ip ;
    for(ip=0;ip<nbPrime;ip++) {
        int32_t p = tbPrime[ip] ;
        for(is=0,iP[is]=ip,P[is]=p, Px[is] = PB204_MAXN/p;is>=0;) {
            if(Px[is]>=1) {
                nbHamm++ ;
                if(iP[is]<nbPrime-1 && tbPrime[iP[is]+1] <= Px[is] ) {
                    iP[is+1] = iP[is]+1 ; is++ ;
                    P[is] = tbPrime[iP[is]] ; Px[is] = Px[is-1]/P[is] ;
                  continue ;
                }
                Px[is] /= P[is] ;
                continue ;
            }
            ++iP[is] ;
            if(is > 0 && iP[is]<nbPrime && tbPrime[iP[is]] <= Px[is-1]) {
                P[is] = tbPrime[iP[is]] ;
                Px[is] = Px[is-1]/P[is] ;
                continue ;
            }
            if(--is >= 0) Px[is] /= P[is] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s \n",pbR->ident) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbHamm);
    return 1 ;
}


int PB204a(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int32_t tbPrime[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    int nbPrime = sizeof(tbPrime)/sizeof(tbPrime[0]) ;
    double * logPrime = malloc(nbPrime*sizeof(logPrime[0])) ;
    int i ;
    for(i=0;i<nbPrime;i++) logPrime[i] = log(tbPrime[i]) ;
    double *ptnbM = logPrime+nbPrime ;
    double *pt1,*pt2,*pt3,*pt4,*pt5,*pt6,*pt7,*pt8,*pt9 ;
    int32_t nbHamm = 1 ; // for 1
    double p1,p2,p3,p4,p5,p6,p7,p8,p9 ;
    double P1,P2,P3,P4,P5,P6,P7,P8,P9 ;
    for(pt1=logPrime;pt1<ptnbM ;pt1++) {
      p1 = *pt1 ;
      for(P1 = log(PB204_MAXN+1) - p1 ;P1>0 ; P1 -= p1) {
        nbHamm++ ;
        for(pt2=pt1+1;pt2<ptnbM && (p2=*pt2) < P1 ;pt2++) {
          for(P2 = P1 - p2 ; P2> 0 ; P2 -= p2) {
            nbHamm++ ;
            for(pt3=pt2+1;pt3<ptnbM && (p3=*pt3) < P2 ;pt3++) {
              for(P3 = P2 - p3 ; P3>0 ; P3 -= p3) {
                nbHamm++ ;
                for(pt4=pt3+1;pt4<ptnbM && (p4=*pt4) < P3 ;pt4++) {
                  for(P4 = P3 - p4 ; P4>0 ; P4 -= p4) {
                    nbHamm++ ;
                    for(pt5=pt4+1;pt5<ptnbM && (p5=*pt5) < P4 ;pt5++) {
                      for(P5 = P4 - p5 ; P5> 0 ; P5 -= p5) {
                        nbHamm++ ;
                        for(pt6=pt5+1;pt6<ptnbM && (p6=*pt6) < P5 ;pt6++) {
                          for(P6 = P5 - p6 ; P6> 0 ; P6 -= p6) {
                            nbHamm++ ;
                            for(pt7=pt6+1;pt7<ptnbM && (p7=*pt7) < P6 ;pt7++) {
                              for(P7 = P6 - p7 ; P7>0 ; P7 -= p7) {
                                nbHamm++ ;
                                for(pt8=pt7+1;pt8<ptnbM && (p8=*pt8) < P7 ;pt8++) {
                                  for(P8 = P7 -p8 ; P8> 0 ; P8 -= p8) {
                                    nbHamm++ ;
                                    for(pt9=pt8+1;pt9<ptnbM && (p9=*pt9) < P8 ;pt9++) {
                                      for(P9 = P8 - p9 ; P9>0 ; P9 -= p9) {
                                        nbHamm++ ;
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s \n",pbR->ident) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbHamm);
    return 1 ;
}





#define PB206_NBD   9

int PB206(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int64_t N,N2 ;
    int i ;
    int is=0;
    int64_t digits[PB206_NBD] ;
    int64_t pow10[2*PB206_NBD] ;
    int64_t maxDigit[PB206_NBD] ;
    pow10[0]=1 ;
    for(i=0;i<2*PB206_NBD-2;i++) pow10[i+1] = pow10[i]*10 ;


    for(i=0;i<PB206_NBD-1;i++) {
        digits[i] = 0 ;
         maxDigit[i]= 9*pow10[i] ;
    }
    
    maxDigit[PB206_NBD-2] = 3*pow10[PB206_NBD-2] ;
    digits[PB206_NBD-1] = pow10[PB206_NBD-1] ;
     maxDigit[PB206_NBD-1] = pow10[PB206_NBD-1] ;
     is=0;
    N=pow10[PB206_NBD-1] ;
    while(1){
        if(is & 1) {
            is++ ; continue ;
        } else {
             if(is<PB206_NBD-1) { // verification of added digit
                N2 = N*N ;
                 if( ( (N2/pow10[is] ) % 10) == 9-(is/2) ) {
                     is++ ; continue ;
                }
            } else { // verification for high digits (low are OK by construction)
                N2 = N*N ;
                for(i=PB206_NBD-2;i>=PB206_NBD/2;i--){
                    if(((N2/pow10[2*i]) % 10) != 9-i) break ;
                }
                if(i<PB206_NBD/2) {
                    break ;
                }
                is-- ;
            }
        }

        while(is>=0 && digits[is]==maxDigit[is]) {
            N -= digits[is] ;
            digits[is--]=0 ;
        }
        if(is< 0) break ;
        else { digits[is] += pow10[is] ; N +=pow10[is] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(is<PB206_NBD-1) {
        return 0 ;
    } else {
        N *= 10 ;
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s %lld **2 = %lld\n",pbR->ident,N,N*N) ;
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",N);
        return 1 ;
    }
}



#define U 10203040506070809LL
#define A  1000000000000000LL
#define B    10000000000000LL
#define C      100000000000LL
#define V        9090909090LL


int PB206a(PB_RESULT *pbR) {
    pbR->nbClock = clock();
    int a,b,c,i;
    int64_t u,m,mm;
    for(a=0;a<10;a++) {
        for(b=0;b<10;b++) {
            for(c=0;c<10;c++) {
                u=U+a*A+b*B+c*C;
                m=Sqrt64(u); m |= 1;
                for(;(mm=m*m)<u+V;m+=2){
                    for(i=9;i>=1;i--) {
                        if(mm%10!=i) break;
                        mm/=100;
                    }
                    if(i==0) {
                        m *= 10 ;
                        pbR->nbClock = clock() - pbR->nbClock ;
                        if(pbR->isVerbose) fprintf(stdout,"\tPB%s %lld **2 = %lld\n",pbR->ident,m,m*m) ;
                        snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",m);
                        return 1 ;
                        
                    }
                }
            }
        }
    }
    return 0 ;
}


 #define PB357_MAXP  100000000
//#define PB357_MAXP  100000000



// brut force version
// phase 1 : Sieve for primes (isP)
// Is Prime generator by Sieve like
// phase 2 : loop on divisor d to eliminate all k*d with k+d not prime
// phase 3 : sum for remianing n from phase 2.
int PB357(PB_RESULT *pbR) {
    int64_t Sum = 0 ;
    int32_t nb = 0 ;
    pbR->nbClock = clock();
    uint8_t  *isP = malloc((PB357_MAXP+1)*sizeof(isP[0]));
    uint8_t  *isPG = malloc((PB357_MAXP+1)*sizeof(isPG[0]));
    memset(isP,1,(PB357_MAXP+1)*sizeof(isP[0])) ;
    memset(isPG,1,(PB357_MAXP+1)*sizeof(isPG[0])) ;
    isP[0] = isP[1] = 0;
    int i, j ;
    for ( i = 2; i <=PB357_MAXP; i++) { // Sieve for primes
        if (isP[i]) {
            for (j = 2*i; j <= PB357_MAXP; j += i) {
                if(isP[j]) isP[j] = 0;
            }
        }
    }
    int d, k ; // Sieve for Prime generator
    for (d = 1; d <= PB357_MAXP; d++) {
        for (k = 1; k*d <= PB357_MAXP; k++) {
            if ( !isP[d+k] && isPG[k*d]) { // not prime and candidate ?
                isPG[k*d] = 0 ; // kill
            }
        }
    }
    for (i = 1; i <= PB357_MAXP; i++) { // sum the remaining candidates
        if (isPG[i]) {
            Sum += i;
            nb++ ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
}
// Phase 1 : compute list of prime < N=100 000 000
// Phase 2: add to Prime generator candidate n=p-1 with p prime
// Phase 3 : for each divisor d < sqrt(N)
//  for each d2*d, { d2*d<N and d2+d not prime } kill d2*d as candidate
//  The test d2+d not prime is done by :
//     search first prime p2 with p2>=2*d (for d2=2)
//     and increment d2 and p2 in parallel
//
int PB357a(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock();
    if((ctxP = Gen_tablePrime(PB357_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    uint8_t *isPG = calloc(PB357_MAXP+1,sizeof(isPG[0]));
    int i;
    for(i=0;i<nbPrime;i++) {
        int n = tbPrime[i]-1 ; // add p-1 as candidate
        isPG[n] = 1 ;
    }
    int ip = 0 ;
    int p = tbPrime[ip] ;
    int imax = Sqrt32(PB357_MAXP);
    int d,d2 ;
    for(d=2;d<=imax;d++) { // loop on divisor
        while(p<2*d) { p=tbPrime[++ip] ; } // search p >=2*i
        int d2Max = PB357_MAXP/d ;
        int ip2 = ip ;
        int p2 = p ;
        d2 = d ;
        while(p2<=d2Max+d) { // p2 next prime
            for(;d2+d<p2;d2++) { // kill all multiple < p2
                if(isPG[d2*d]) isPG[d2*d] = 0 ;
            }
            p2 = tbPrime[++ip2] ; d2++ ; // skip p2, and next p2
        }
        for(;d2<=d2Max;d2++) { // end loop
            isPG[d2*d] = 0 ;
        }
    }
    int64_t Sum = 0 ;
    int32_t nb = 0 ;
    for(i=1;i<PB357_MAXP;i++) {
        if(isPG[i]) { nb++ ; Sum += i ; }
    }
    Free_tablePrime(ctxP) ;
    free(isPG) ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
}

// assume nbFact > 2
// we know that 2 is a prime factor of n, so decompose n=d1*d2
// with d1 odd and d2 even.
static int Chk357(int nbFact,int *factP,int n,const T_prime * tbPrime) {
    int nbDiv2 = 1<< (nbFact -1) ;
    int mask ;
    for(mask=0;mask<nbDiv2;mask++) { // all partition of prime factors.
        int d1 = 1 ;
        int d2 = 2 ;
        int i ;
        for(i=0;i<nbFact-1;i++ ) {
            if((1<<i) & mask) d1 *= factP[i+1] ;
            else d2 *= factP[i+1] ;
        }
        if(!Is_Prime32(d1+d2,tbPrime)) return 0 ;
    }
    return 1 ;
}


// brut force using decompostion of n in product of prime factors
//   and using the fact that n mstt have no square factor.
// Phase 1 : compute list of prime < N=100 000 000
// Phase 2 :
//   for each prime p , p==3 % 4 and (p+3)/3 prime
//    so d=1 and d=2 are OK and 2 is not a square factor
//    decompose n in prime factor n=p-1 = p1 * p2 * p3 ... pk (check no square)
//    then check for each divisor d of n that d=p+n/d is prime.
//    The list of divisor d is build using n=p1 * p2 * ... * pk (2**k divisor <=> 2**(k-1) divisor pair)
int PB357b(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock();
    if((ctxP = Gen_tablePrime(PB357_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    int64_t Sum = 1 ;
    int32_t nb = 1 ;
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int factP[40] ;
    factP[0] = 2 ;
    int i,j;
    for(i=1;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        int n  = p-1 ; // so d=1 is OK
        if((p & 3)==1) continue ; // 2 is square
        if(!Is_Prime32((p+3)/2, tbPrime)) continue ; // so d=2 is OK
        int n1 = n ; // save n
        int nbFact = 0 ; // search prime factors of n
        int pj;
        for(j=0;pj=tbPrime[j],pj*pj<=n;j++){
            if((n % pj) == 0) {
                n /= pj ;
                if((n % pj)== 0 ) { n=0 ; break ; } // square
                factP[nbFact++] = pj ;
            }
        }
        if(n==0) continue ;
        factP[nbFact++] = n ; // the last factor is prime
        if(nbFact == 2 || Chk357(nbFact, factP,n1,tbPrime)) {
            nb++; Sum += n1 ;
        }
    }
    Free_tablePrime(ctxP) ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
    
}

// Phase 1 : compute list of prime < N=100 000 000
// Phase 2 :
//   for each prime p , p==3 % 4 and (p+3)/3 prime
//    so d=1 and d=2 are OK and 2 is not a square factor
//    for n=p-1 initialise decomposition count nbPG(n) to one (d=1 as 1+n/1 prime)
//    use the fact that p==3 % 4 (to avoid 2 square divisor)  to store only N/4 values
// Phase 3 :
//  for each prime p > 2
//      a) for each d=2 d<p/2 increment decompostion count nbPG(n) for n=d*(p-d) so d+n/d is prime
//      b) for each multiple kp = k*p increment the count of prime divisor nbDivP(kp)
//          if k is a multipe of p (p square prime divisor) kill kp by setting nbPG(kp)=0
// Phase 4 :
//   Sum for n with nbPG(n) == 2**nbDivP(n)
//   as if n=p1*p2*..pk with k==nbDivP(n) , n as 2**k divisor.
typedef struct IS_PG {
    uint8_t nbPG ; //
    uint8_t nbDivP ;
} IS_PG ;
int PB357c(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    int64_t Sum = 0 ;
    int32_t nb = 0 ;
    pbR->nbClock = clock();
    if((ctxP = Gen_tablePrime(PB357_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    IS_PG  *pg = calloc(PB357_MAXP/4+1,sizeof(pg[0]));
    int i,k,n;
    for(i=0;i<nbPrime;i++) { // n = p-1 , 2 not square divisor of n
        int n = tbPrime[i] -1 ;
        if((n&3) == 2) pg[n/4].nbPG = 1 ;
    }
    for(i=1;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        int d ; // increment nbPG((p-d)*d)
        for(d=2;(n=(p-d)*d)<PB357_MAXP && 2*d<=p;d++){
            if( (n/2&1) && pg[n/4].nbPG) pg[n/4].nbPG++ ;
        }
        int kp ; // increment nbDivP for kp.
        for (k=1,kp=p;kp<=PB357_MAXP/2;k++,kp+= p) {
            if(k!=p) {
                if((kp&1) && pg[kp/2].nbPG) pg[kp/2].nbDivP++ ;
            } else { // kill kp as k multiple of p, so p square divisor.
                if((kp&1) && pg[kp/2].nbPG) pg[kp/2].nbPG=0 ;
                k = 0 ;
            }
        }
     }
     for(i=1;i<nbPrime;i++) { // sum for nbPG == 2**nbDivP
        int p = tbPrime[i] ;
         if((p&3) != 3) continue ;
        int n  = (p-1)/4 ;
        if(pg[n].nbPG && pg[n].nbPG == ( 1<< (pg[n].nbDivP) ) ){
            nb++; Sum += n ;
        }
    }
    Sum = 4*Sum+ 2*nb + 1 ; // correction 2*nb for ((4n1+2)/4)*4 = 4n1+2
    nb++ ;
    Free_tablePrime(ctxP) ;
    free(pg) ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Nb=%d Sum= %lld\n",pbR->ident,nb,Sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
}



