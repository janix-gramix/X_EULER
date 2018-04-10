//
//  PB101_150.c
//  X_euler
//
//  Created by Jeannot on 19/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



#include "PB101_150.h"
#include "p102_data.h"
#include "p105_data.h"
#include "p107_data.h"

#define PB101_DEG   10
int64_t P(int32_t n) {
    // P(n) = (n**11+1)/(n+1)
    int64_t N = n ;
    int64_t N2 = N*N ;
    int64_t N4=N2*N2 ;
    int64_t N8 = N4*N4 ;
    int64_t Pn = (N8*N2*N+1)/(N+1) ;
    return Pn ;
}
int PB101(PB_RESULT *pbR) {
    int64_t P_[PB101_DEG+3] ;
    //Q[n]=Qn(n+1) = OP(n,n+1) ; D[n] = P[n+1] - OP(n,n+1)
    int64_t Q[PB101_DEG+2] ;
    int64_t D[PB101_DEG+2] ;

    pbR->nbClock = clock() ;
    Q[0] = 0 ;
    P_[1] =P(1) ;
    D[0] = P_[1]-Q[0];
    
    Q[1] = D[0] ;
    P_[2] =P(2) ;
    D[1] = P_[2] - Q[1] ;

    //Q2(2) = P(2)
    // Q2(n) = Q1(n)+ (n-1)/1! D[1] => Q2(3)=Q[2]
    Q[2] = D[0]
    + (3-1) * D[1] / 1 ;
    P_[3] =P(3) ;
    D[2] = P_[3]-Q [2] ;

    // Q3(n) = Q2(n) + (n-1)(n-2))/2! D[2]  ;  Q3(4)=Q[3]
    Q[3] = D[0]
    +  (4-1)*D[1]
    + ((4-1)*(4-2))/2 * D[2] ;
    P_[4] =P(4) ;
    D[3] = P_[4]-Q[3] ;
 
    // Q4(n) = Q3(n) + (n-1)(n-2)(n-3))/3! D[3]  Q4(5) = Q[4]
    Q[4] = D[0]
    +  (5-1)*D[1]
    + ((5-1)*(5-2))/2 * D[2]
    + ((5-1)*(5-2)*(5-3))/(2*3)  *D[3] ;
    P_[5] =P(5) ;
    D[4] = P_[5] - Q[4] ;
   
    // Q5(n) = Q4(n) + (n-1)(n-2)(n-3)(n-4) D[4] / 4! Q5(6) = Q[5]
    Q[5] = D[0]
    +  (6-1)*D[1]
    + ((6-1)*(6-2))/2 * D[2]
    + ((6-1)*(6-2)*(6-3))/(2*3) * D[3]
    + ((6-1)*(6-2)*(6-3)*(6-4))/(2*3*4) *D[4] ;
    P_[6] =P(6) ;
    D[5] = P_[6] - Q[5] ;
   
    // Q6(n) = Q5(n) + ((n-1)(n-2)(n-3)(n-4)(n-5))/5!  D[5]  Q6(7) = Q[6]
    Q[6] = D[0]
    +  (7-1)*D[1]
    + ((7-1)*(7-2))/2 * D[2]
    + ((7-1)*(7-2)*(7-3))/(2*3) * D[3]
    + ((7-1)*(7-2)*(7-3)*(7-4))/(2*3*4) *D[4]
    + ((7-1)*(7-2)*(7-3)*(7-4)*(7-5))/(2*3*4*5) *D[5] ;
    P_[7] =P(7) ;
    D[6] = P_[7] - Q[6] ;
 
    // Q7(n) = Q6(n) + ((n-1)(n-2)(n-3)(n-4)(n-5)(n-6))/6!  D[6]  Q7(8) = Q[7]
    Q[7] = D[0]
    +  (8-1)*D[1]
    + ((8-1)*(8-2))/2 * D[2]
    + ((8-1)*(8-2)*(8-3))/(2*3) * D[3]
    + ((8-1)*(8-2)*(8-3)*(8-4))/(2*3*4) *D[4]
    + ((8-1)*(8-2)*(8-3)*(8-4)*(8-5))/(2*3*4*5) *D[5]
    + ((8-1)*(8-2)*(8-3)*(8-4)*(8-5)*(8-6))/(2*3*4*5*6) *D[6] ;
    P_[8] =P(8) ;
    D[7] = P_[8] - Q[7] ;

    // Q8(n) = Q7(n) + ((n-1)(n-2)(n-3)(n-4)(n-5)(n-6)(n-7))/7!  D[7]  Q8(9) = Q[8]
    Q[8] = D[0]
    +  (9-1)*D[1]
    + ((9-1)*(9-2))/2 * D[2]
    + ((9-1)*(9-2)*(9-3))/(2*3) * D[3]
    + ((9-1)*(9-2)*(9-3)*(9-4))/(2*3*4) *D[4]
    + ((9-1)*(9-2)*(9-3)*(9-4)*(9-5))/(2*3*4*5) *D[5]
    + ((9-1)*(9-2)*(9-3)*(9-4)*(9-5)*(9-6))/(2*3*4*5*6) *D[6]
    + ((9-1)*(9-2)*(9-3)*(9-4)*(9-5)*(9-6)*(9-7))/(2*3*4*5*6*7) * D[7] ;
    P_[9] =P(9) ;
    D[8] = P_[9] - Q[8] ;

    // Q9(n) = Q8(n) + ((n-1)(n-2)(n-3)(n-4)(n-5)(n-6)(n-7)(n-8))/8!  D[8]  Q9(10) = Q[9]
    Q[9] = D[0]
    + (10-1)*D[1]
    + ((10-1)*(10-2))/2 * D[2]
    + ((10-1)*(10-2)*(10-3))/(2*3) * D[3]
    + ((10-1)*(10-2)*(10-3)*(10-4))/(2*3*4) *D[4]
    + ((10-1)*(10-2)*(10-3)*(10-4)*(10-5))/(2*3*4*5) *D[5]
    + ((10-1)*(10-2)*(10-3)*(10-4)*(10-5)*(10-6))/(2*3*4*5*6) *D[6]
    + ((10-1)*(10-2)*(10-3)*(10-4)*(10-5)*(10-6)*(10-7))/(2*3*4*5*6*7) * D[7]
    + ((10-1)*(10-2)*(10-3)*(10-4)*(10-5)*(10-6)*(10-7)*(10-8))/(2*3*4*5*6*7*8) * D[8] ;
    P_[10] =P(10) ;
    D[9] = P_[10] - Q[9] ;
    
    // Q10(n) = Q9(n) + ((n-1)(n-2)(n-3)(n-4)(n-5)(n-6)(n-7)(n-8)(n-9))/9!  D[9]  Q10(11) = Q[10]
    Q[10] = D[0]
    + (11-1)*D[1]
    + ((11-1)*(11-2))/2 * D[2]
    + ((11-1)*(11-2)*(11-3))/(2*3) * D[3]
    + ((11-1)*(11-2)*(11-3)*(11-4))/(2*3*4) *D[4]
    + ((11-1)*(11-2)*(11-3)*(11-4)*(11-5))/(2*3*4*5) *D[5]
    + ((11-1)*(11-2)*(11-3)*(11-4)*(11-5)*(11-6))/(2*3*4*5*6) *D[6]
    + ((11-1)*(11-2)*(11-3)*(11-4)*(11-5)*(11-6)*(11-7))/(2*3*4*5*6*7) * D[7]
    + ((11-1)*(11-2)*(11-3)*(11-4)*(11-5)*(11-6)*(11-7)*(11-8))/(2*3*4*5*6*7*8) * D[8]
    + ((11-1)*(11-2)*(11-3)*(11-4)*(11-5)*(11-6)*(11-7)*(11-8)*(11-9))/(2*3*4*5*6*7*8*9) * D[9];
    P_[11] =P(11) ;
    D[10] = P_[11] - Q[10] ;

    // Q11(n) = Q10(n) + ((n-1)(n-2)(n-3)(n-4)(n-5)(n-6)(n-7)(n-8)(n-9)(n-10))/10!  D[10]  Q11(12) = Q[11]
    Q[11] = D[0]
    + (12-1)*D[1]
    + ((12-1)*(12-2))/2 * D[2]
    + ((12-1)*(12-2)*(12-3))/(2*3) * D[3]
    + ((12-1)*(12-2)*(12-3)*(12-4))/(2*3*4) *D[4]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5))/(2*3*4*5) *D[5]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5)*(12-6))/(2*3*4*5*6) *D[6]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5)*(12-6)*(12-7))/(2*3*4*5*6*7) * D[7]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5)*(12-6)*(12-7)*(12-8))/(2*3*4*5*6*7*8) * D[8]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5)*(12-6)*(12-7)*(12-8)*(12-9))/(2*3*4*5*6*7*8*9) * D[9]
    + ((12-1)*(12-2)*(12-3)*(12-4)*(12-5)*(12-6)*(12-7)*(12-8)*(12-9)*(12-10))/(2*3*4*5*6*7*8*9*10) * D[10];
    P_[12] =P(12) ;
    D[11] = P_[12] - Q[11] ;
   
    int64_t S=0 ;
    int k ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s S",pbR->ident);
    for(k=1;k<PB101_DEG+2;k++) {
        if(D[k]) {
            S += Q[k] ;
            if(pbR->isVerbose) fprintf(stdout,"%c%lld",(k==1) ? '=' : '+',Q[k]);
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"=%lld\n",S);
    sprintf(pbR->strRes,"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB102(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    
    const Triangle * TR =  P102_GetData() ;
    int i ;
    int nbOK = 0;
    for(i=0;i<PB102_NBT;i++) {
        int32_t A = - TR[i].yb*TR[i].xc +TR[i].yc*TR[i].xb ;
        int32_t B = - TR[i].yc*TR[i].xa +TR[i].ya*TR[i].xc ;
        if(A > 0 && B > 0) {
            int32_t C = - TR[i].ya*TR[i].xb +TR[i].yb*TR[i].xa ;
            if(C > 0) nbOK++ ;
        } else if(A < 0 && B < 0) {
            int32_t C = - TR[i].ya*TR[i].xb +TR[i].yb*TR[i].xa ;
            if(C < 0) nbOK++ ;
        }
    }
    sprintf(pbR->strRes,"%d",nbOK) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

/*
double t = r * log10((1+sqrt(5))/2) + log10(1/sqrt(5))
long L = (long)Math.pow(10, t - (long)t + 8)
*/
static int CmpUint32(const void *el1, const void *el2) {
    return ((int32_t *)el1)[0] - ((int32_t *)el2)[0] ;
}
#define FACT9   362880
int PB104(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int k ;
     u_int32_t PanDigital[FACT9] ;
    u_int8_t perm[9] = {1,2,3,4,5,6,7,8,9} ;
    int is=0 ;
    do {
        PanDigital[is++] = 10*(10*(10*(10*(10*(10*(10*(10*perm[0]+perm[1])+perm[2])+perm[3])+perm[4])+perm[5])+perm[6])+perm[7])+perm[8] ;
    } while (NextPermut(perm,9) >= 0) ;
    u_int32_t F0=1 ;
    u_int32_t F1=1 ;
    k = 2 ;
    while(1) {
        k++ ;
        F0=(F0+F1) % 1000000000 ;
        if((F0 % 9 ) == 0 ) {
           if(bsearch(&F0, PanDigital, FACT9,sizeof(u_int32_t), CmpUint32) != NULL) {
                // approximation by PHY**k / sqrt(5)
                double logFk = k * log10((1+sqrt(5))/2) - log10(sqrt(5)) ;
                u_int32_t highDigits = (u_int32_t)pow(10,logFk - (u_int32_t)logFk +8 ) ;
                if(bsearch(&highDigits, PanDigital, FACT9,sizeof(u_int32_t), CmpUint32) != NULL) {
                    printf(" F(%d) is Double Pandigit\n",k) ;
                    break ;
                }
            }
        }
        u_int32_t tmp = F0 ;
        F0=F1 ;
        F1=tmp ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",k) ;
    return 1 ;
}




typedef struct Edge {
    u_int16_t   cost;
    u_int16_t   Vbeg ;
    u_int16_t   Vend ;
} Edge ;

typedef struct vertexTree { // tree for a vertex
    u_int16_t numTree ; // tree number.
    u_int16_t nxtVertex ; // link to nxt vertex in the tree
    
} nodeTree ;


int CmpEdge( const void *edg1, const void *edg2) {
    return ((Edge *)edg1)[0].cost - ((Edge *)edg2)[0].cost ;
}
int PB107(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    Edge EDG[PB107_SIZE*PB107_SIZE] ;
    u_int16_t treeLength[PB107_SIZE+1] ;
    nodeTree VtoT[PB107_SIZE] ; // association node -> tree
    u_int16_t nbTree = 0 ;
    int i, j ;
    const u_int16_t * cost = P107_GetData();
    u_int16_t nbEdge = 0 ;
    for(i=0;i<PB107_SIZE;i++) { // get the non null edge
        for(j=i+1;j<PB107_SIZE;j++) {
            if(cost[i*PB107_SIZE+j]){
                EDG[nbEdge].Vbeg = i ;
                EDG[nbEdge].Vend = j ;
                EDG[nbEdge].cost = cost[i*PB107_SIZE+j] ;
                nbEdge++ ;
            }
        }
    }
    treeLength[nbTree++] = 1 ; // special tree with one node
    for(i=0;i<PB107_SIZE;i++) { //each node in is own tree
        VtoT[i].numTree = 0 ; // can be the same tree as this tree is never used
        VtoT[i].nxtVertex = i ;
    }
    heapsort(EDG,nbEdge,sizeof(EDG[0]),CmpEdge) ; // sort edges by cost
    int indSortEdg = 0 ;
    int maxLength = 0 ;
    u_int32_t savingCost = 0 ;
    do {    // loop to build tree , adding the min cost edge
        Edge curEdg = EDG[indSortEdg++] ;
        u_int16_t begT= VtoT[curEdg.Vbeg].numTree ;
        u_int16_t endT =VtoT[curEdg.Vend].numTree ;
        if(begT && endT && begT == endT) {
            savingCost += curEdg.cost ; // the edge address the same tree => loop, forbiddeen
            continue ;
        }
        if(begT == 0 && endT == 0 ) { // nodes not in e tree => simple tree with 2 nodees
                treeLength[nbTree] = 2 ;
                if(maxLength < 2) maxLength = 2 ;
                VtoT[curEdg.Vbeg].numTree = VtoT[curEdg.Vend].numTree = nbTree ;
                VtoT[curEdg.Vbeg].nxtVertex = curEdg.Vend ;
                VtoT[curEdg.Vend].nxtVertex = curEdg.Vbeg ;
            nbTree++ ;
                continue;
        } else { // two trees, rattach the shortest to the longuest
            if(treeLength[begT] > treeLength[endT] ) { // permutation end <->beg
                u_int16_t tmp = endT ;
                endT = begT ;
                begT= tmp ;
                tmp = curEdg.Vbeg ;
                curEdg.Vbeg = curEdg.Vend ;
                curEdg.Vend = tmp ;
            }
          // rattach beg tree to end tree
            u_int16_t nxtEnd = VtoT[curEdg.Vend].nxtVertex ;
            u_int16_t nxtBeg = curEdg.Vbeg ;
            VtoT[curEdg.Vend].nxtVertex = nxtBeg ;
            u_int16_t antBeg ;
            do { // loop on the vertexes of beg tree
                antBeg= nxtBeg ;
                VtoT[nxtBeg].numTree = endT ;
                nxtBeg = VtoT[nxtBeg].nxtVertex ;
            } while(nxtBeg != curEdg.Vbeg) ;
            VtoT[curEdg.Vbeg].numTree = endT ;
            VtoT[antBeg].nxtVertex = nxtEnd ;
            treeLength[endT] += treeLength[begT] ;
            if(maxLength < treeLength[endT]) maxLength = treeLength[endT] ;
            continue ;
        }
    } while(maxLength != PB107_SIZE) ;
    while(indSortEdg<nbEdge) { // add the remaining edge to cost saving
        savingCost += EDG[indSortEdg++].cost ;
    }
    sprintf(pbR->strRes,"%d",savingCost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

// les solutions sont du type <=> xy = n (x+y)
// n = d1 x d2 x d3 avec  d2 <= d3 et d2^d3 pour assurer l'unicite
// x= d1 x d3 x (d2+d3) ;  y = x= d1 x d2 x (d2+d3)
// Ex n=4
// d1=1 d2=1 d3=4 => x=20 y=5  ;100 = 4 x(20+5)
// d1=2 d2=1 d3=2 => x=12 y=6   ;72 = 2 x (6+12)
// d1=4 d2=1 d3=1 => x=8  y=8  ; 64 = 4 x (8+8)
//
// plus astucieux
// 1/(n+a) + 1/(n+b) = 1/n <=> (n+a)(n+b)=n(2n+a+b) <=>ab = n2
// donc c'est le (nombre de diviseurs de n**2)/2 + 1. 


int PB108(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int d1,d2,d3 ;
    int n ;
    int nbSol = 0 ;
    int nbSolM=0 ;
    for(n=4;;n++) {
        nbSol = 0 ;
        for(d2=1;d2*d2<=n;d2++) {
            int n1 = n / d2 ;
            if(n != d2*n1) continue ;
            int g =PGCD(d2,n1);
            for(d1=g;d1*d2 <= n1 ;d1 += g) {
                d3 = n1/d1 ;
                if(n1 == d1*d3) {
                    if( PGCD(d2,d3) == 1)
                    {
//                           printf("(%d,%d)",d1*d2*(d2 + d3),d1*d3*(d2 + d3)) ;
                        nbSol++ ;
                    }
                }
            }
        }
        if(nbSol > nbSolM) {
            printf("%d=>%d \n",n,nbSol );
            nbSolM = nbSol ;
            if(nbSol > 1000) break ;
        }
    }
 
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",n) ;
    return 1 ;
}

#define PB108_ALPHAM    60
#define PB110_MINS      4000000
#define PB108_MINS      1000
int PB110(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int alpha[PB108_ALPHAM+1] ;
    int primes[PB108_ALPHAM] = { 2,3,5,7,11,13,17,19,23,29
                                ,31,37,41,43,47,53,59,61,67,71
                                ,73,79,83,89,97,101,103,107,109,113
                                ,127,131,137,139,149,151,157,163,167,173
                                ,179,181,191,193,197,199,211,223,227,229
                                ,233,239,241,251,257,263,269,271,277,281} ;
    int minS = (strcmp(pbR->ident,"110") == 0) ? PB110_MINS : PB108_MINS ;
    int sumA ;
     u_int64_t nMin = 1L << PB108_ALPHAM ;
    for(sumA=2;(1L<<sumA) <= nMin ;sumA++) {
        //      printf("\nSum(alpha)=%d ",sumA) ;
        int nbA,nbAmax ;
         u_int64_t n=1 ;
        for(nbAmax=0;nbAmax<=sumA;nbAmax++) {
            n *= primes[nbAmax] ;
            if(n> nMin) break ;
        }
        u_int64_t nbS = 1 ;
        for(nbA=1;nbA<nbAmax;nbA++) {
            alpha[nbA] = 1 ;
            nbS *= 3 ;
            if((2*(sumA-nbA)+1)*nbS > 2*minS) break ;
        }
        alpha[--nbA] = 0 ;
        alpha[0] = sumA-nbA+1 ;
        while(1) {
            int j ;
            int ia;
            u_int64_t nbS = 1 ;
            for(ia=0;ia<nbA;ia++) {
                nbS *= 2*alpha[ia]+1 ;
            }
            u_int64_t nbSol = (nbS+1)/2 ;
//            for(ia=0;ia<nbA;ia++) printf("%c%d",(ia==0)? '\n' : '.' ,alpha[ia]) ; printf(" nbS=%lld",nbSol) ;
           if(nbSol > minS ) {
                u_int64_t n=1;
                for(ia=nbA-1;ia>=0;ia--) {
                    for(j=alpha[ia];j>0;j--) {
                        n *= primes[ia] ;
                    }
                    if(n > nMin) break ;
                }
                if(n < nMin) {
                    nMin = n ;
                    if(pbR->isVerbose) {
                        fprintf(stdout,"\t PB%s SumA=%d Alpha=",pbR->ident,sumA);
                        for(ia=0;ia<nbA;ia++) fprintf(stdout,"%c%d",(ia==0)? ' ' : '.' ,alpha[ia]) ;
                        printf(" NbSol=%lld,n=%lld\n",nbSol,nMin) ;
                    }
                }
                break ;
            } else {
                int sumR ;
                for(sumR=1,j=nbA-1;j>=0;j--) {
                    sumR += alpha[j] ; // sumR to dispatch betwwen the remaining alpha's
                    if(sumR <= (nbAmax-j)*(alpha[j]-1)) {
                        sumR -= alpha[j] ;
                        alpha[j]-- ;
                        int k= j+1 ;
                        while(sumR > 0) {
                            alpha[k] = (sumR > alpha[j]) ? alpha[j] : sumR ;
                            sumR -= alpha[k++] ;
                        }
                        nbA = k ;
                        alpha[k] = 0 ;
                        break ;
                    }
                }
                if(j<0)    break ; // no successor
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nMin) ;
    return 1 ;
}

// the searched proportion is PB112_PERCENT / (PB112_PERCENT + 1)
// can be 999999999 for a proportion af 1/one billion no bouncy
#define PB112_PERCENT   99
#define PB112_MAXDIG    21
// #define PB112_DEBUG
// struct to count the different number categories
// for increasing and decreasing number differentiation of the numbers
// by the leading digit.
// Also, the "constant numbers" are not included in incr[] and decr[]
// to avoid multiples counts.
// Leading '0' is only counted for recursion.
typedef struct CountB {
    int64_t Incr[10] ;  // increasing by leading digit (constant excluded)
    int64_t Decr[10] ;  // increasing by leading digit (constant excluded)
    int64_t Bouncy0 ;   // Bouncy with '0' as leading digit
    int64_t BouncyN0 ;  // Bouncy with no '0' as leading digit
    int64_t Const ;     // Constant numbers (same digit)
// cumulative counters for recursion
    int64_t SumI ;
    int64_t SumD ;
    int64_t SumB ;
    int64_t SumC ;
} CountB ;

// info is only a comment
// nbDigJoker is only a commodity to print the last digits
void CountB_print(CountB *CB,char *info, int nbDigJoker) {
    int64_t SI = CB->SumI;
    int64_t SD = CB->SumD ;
    int64_t SB = CB->SumB ;
    int64_t SC = CB->SumC ;
    int64_t ST = SI+SD+SB -SC ; // verification
    char * joker=".................9" ;
#if defined(PB112_DEBUG)
    int j ;
    int64_t I = CB->Const  ; for(j=1;j<10;j++) I += CB->Incr[j] ;
    int64_t D = CB->Const  ; for(j=1;j<10;j++) D += CB->Decr[j] ;
    int64_t B = CB->BouncyN0 ;
    int64_t T = I+D+B - CB->Const ;
    printf(" %s%s => %6lld I +%6lld D + %10lld B T=%10lld Cum: %6lld I +%6lld D +%10lld B T=%10lld Perc=%.6lf\n"
           ,info,joker+17-nbDigJoker,I,D,B,T,SI,SD,SB,ST, ((double)100.0*SB)/ST) ;
#else
    printf(" %s%s => %6lld I +%6lld D +%10lld B T=%10lld Perc=%.8lf\n"
           ,info,joker+17-nbDigJoker,SI,SD,SB,ST, (((double)100.0)*SB)/ST);
#endif
}

// initialise newC with oldC
// reset the current counters, tansmit cumulatives counters
void CountB_Init(CountB *oldC, CountB *newC) {
    newC->SumI = oldC->SumI ;
    newC->SumD = oldC->SumD ;
    newC->SumB = oldC->SumB ;
    newC->SumC = oldC->SumC ;
    memset(newC->Incr,0,sizeof(newC->Incr));
    memset(newC->Decr,0,sizeof(newC->Decr));
    newC->Bouncy0 = newC->BouncyN0 = newC->Const = 0 ;
}

#define BIT_INC 1   // increasing
#define BIT_DEC 2   // decreasing
#define BIT_CONST 4 // increasing and decreasing

// add leading digits chain to statistics
//
void CountB_AddHead(char *digits, CountB *oldC, CountB *newC) {
    // a chain can be increasing and decreasing (constant chain)
    // check the status (increasing, decreasing of the leading digits chain)
    int status = BIT_INC | BIT_DEC | BIT_CONST ;
    int idH = digits[0] - '0' ; // leading digit of the added digit chain
    int idT = idH ;             // last digit of the chain, to compute butting
    int is ;
    for(is=0;digits[is] != 0; is++) {
        idT = digits[is] - '0' ;
        if(is > 0) {
            if(digits[is] > digits[is-1]) status &= ~(BIT_DEC|BIT_CONST) ;
            else if(digits[is] < digits[is-1]) status &= ~(BIT_INC|BIT_CONST) ;
        }
    }
    int ida; // leader digit of the precedent stats
    // deltaBouncy will be attributed to Boncy0 or BouncyN0 depending on idH (leading digt of the added chain)
    int64_t deltaBouncy = oldC->Bouncy0 + oldC->BouncyN0 ;
    if(oldC->Const) {
        for(ida=0;ida<10;ida++) {
            if(ida < idT) {  // T > A
                if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[ida] + 1;
                else  deltaBouncy += oldC->Decr[ida] + 1 ;
                deltaBouncy += oldC->Incr[ida] ;
            }else if (ida > idT) {  // T < A
                if(status & BIT_INC)  newC->Incr[idH] += oldC->Incr[ida]  + 1 ;
                else deltaBouncy += oldC->Incr[ida] + 1 ;
                deltaBouncy += oldC->Decr[ida] ;
            } else { // T == A
                if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[ida] + ((status & BIT_CONST) ? 0 : 1)  ;
                else deltaBouncy += oldC->Decr[ida]   ;
                if(status & BIT_INC) newC->Incr[idH] += oldC->Incr[ida]  + ((status & BIT_CONST) ? 0 : 1) ;
                else  deltaBouncy += oldC->Incr[ida]   ;
                if(! (status & (BIT_INC|BIT_DEC)))  deltaBouncy++ ;
            }
        }
    } else { // special case for one digit numbers  to avoid problems with constant numbers
        if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[idT] + ((status & BIT_CONST) ? 0 : 1)  ;
        else deltaBouncy += oldC->Decr[idT]   ;
        if(status & BIT_INC) newC->Incr[idH] += oldC->Incr[idT]  + ((status & BIT_CONST) ? 0 : 1) ;
        else  deltaBouncy += oldC->Incr[idT]   ;
        if(! (status & (BIT_INC|BIT_DEC)))  deltaBouncy++ ;
    }
    if(idH==0) newC->Bouncy0 += deltaBouncy ;
    else { // accumulation
        newC->BouncyN0 += deltaBouncy ;
        newC->SumB += deltaBouncy ;
        newC->SumI += newC->Incr[idH] + ((status & BIT_CONST) ? 1  : 0) ;
        newC->SumD += newC->Decr[idH] + ((status & BIT_CONST) ? 1  : 0) ;
        if(status & BIT_CONST) {
            newC->Const++ ;
            newC->SumC++ ;
        }
    }
}


int PB112(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CountB CB[PB112_MAXDIG+1] ; // to store statistic for numbers  10**(nd-1) : 9, 99, 999, 9999,
    char prefix[PB112_MAXDIG+1] ;
    int id=1,nd,j ;
    int64_t n ;
    CountB CBwork ; // 2 counts so when depassing the purpose, rolling back one step
    CountB CBnext ;
    memset(&CB[0],0,sizeof(CB[0])) ;
    CountB_Init(CB,&CBwork) ;
    // search the power of 10 reaching a superior proportion of bouncy
    for(nd=1;nd<PB112_MAXDIG-1;nd++) {
         CountB_AddHead("0",CB+nd-1,&CBwork) ; // no print, only for recursion
         CBnext=CBwork ;
         for(id=1;id<10;id++) { // add the
            sprintf(prefix,"%d",id) ;
            CountB_AddHead(prefix,CB+nd-1,&CBnext) ;
            if(PB112_PERCENT * (CBnext.SumI+CBnext.SumD-CBnext.SumC) <= CBnext.SumB) break ;
            CountB_print(&CBnext,prefix,nd-2) ;
            CBwork = CBnext ; // valid power, next turn.
        }
        printf("\n");
        if(id != 10) {  break ; // upper bound reached
        } else {
            CB[nd] = CBwork ;
            CountB_Init(CB+nd,&CBwork) ;
        }
     }
    n = id ; // leading digit for the search, CBwork contains statistic for numbers <= [id]9999999
    while(--nd > 0) { // loop to search lower digits
        n *=10 ;
         for(j=0;j<10;j++) { // search the digit (upper bound)
             int isLess = 0 ;
            sprintf(prefix,"%lld",n+j) ;
             CountB_Init(&CBwork,&CBnext) ;
             CountB_AddHead(prefix, CB+nd-1,&CBnext) ;
             if(PB112_PERCENT * (CBnext.SumI+CBnext.SumD-CBnext.SumC) < CBnext.SumB) {
                 isLess = 1; // upper bound found, => next lower dgit
             }
             printf(" n %c ",isLess ? '<' : ' ') ; CountB_print(&CBnext,prefix,nd-2) ;
             if(isLess) break ;
             CBwork = CBnext ; // valid digit, next turn.
        }
        n += j ; // add the valid digit
    }
    n = CBwork.SumI+CBwork.SumD + CBwork.SumB - CBwork.SumC ;
    int64_t nbBouncy = CBwork.SumB ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Under n=%lld exactly only %lld are not bouncy numbers\n",pbR->ident,n,n-nbBouncy) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",n) ;
    return 1 ;
}
