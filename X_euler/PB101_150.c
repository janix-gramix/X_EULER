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
#include "p122_data.h"

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
     uint32_t PanDigital[FACT9] ;
    uint8_t perm[9] = {1,2,3,4,5,6,7,8,9} ;
    int is=0 ;
    do {
        PanDigital[is++] = 10*(10*(10*(10*(10*(10*(10*(10*perm[0]+perm[1])+perm[2])+perm[3])+perm[4])+perm[5])+perm[6])+perm[7])+perm[8] ;
    } while (NextPermut(perm,9) >= 0) ;
    uint32_t F0=1 ;
    uint32_t F1=1 ;
    k = 2 ;
    while(1) {
        k++ ;
        F0=(F0+F1) % 1000000000 ;
        if((F0 % 9 ) == 0 ) {
           if(bsearch(&F0, PanDigital, FACT9,sizeof(uint32_t), CmpUint32) != NULL) {
                // approximation by PHY**k / sqrt(5)
                double logFk = k * log10((1+sqrt(5))/2) - log10(sqrt(5)) ;
                uint32_t highDigits = (uint32_t)pow(10,logFk - (uint32_t)logFk +8 ) ;
                if(bsearch(&highDigits, PanDigital, FACT9,sizeof(uint32_t), CmpUint32) != NULL) {
                    printf(" F(%d) is Double Pandigit\n",k) ;
                    break ;
                }
            }
        }
        uint32_t tmp = F0 ;
        F0=F1 ;
        F1=tmp ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",k) ;
    return 1 ;
}




typedef struct Edge {
    uint16_t   cost;
    uint16_t   Vbeg ;
    uint16_t   Vend ;
} Edge ;

typedef struct vertexTree { // tree for a vertex
    uint16_t numTree ; // tree number.
    uint16_t nxtVertex ; // link to nxt vertex in the tree
    
} nodeTree ;


int CmpEdge( const void *edg1, const void *edg2) {
    return ((Edge *)edg1)[0].cost - ((Edge *)edg2)[0].cost ;
}
int PB107(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    Edge EDG[PB107_SIZE*PB107_SIZE] ;
    uint16_t treeLength[PB107_SIZE+1] ;
    nodeTree VtoT[PB107_SIZE] ; // association node -> tree
    uint16_t nbTree = 0 ;
    int i, j ;
    const uint16_t * cost = P107_GetData();
    uint16_t nbEdge = 0 ;
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
    qsort(EDG,nbEdge,sizeof(EDG[0]),CmpEdge) ; // sort edges by cost
    int indSortEdg = 0 ;
    int maxLength = 0 ;
    uint32_t savingCost = 0 ;
    do {    // loop to build tree , adding the min cost edge
        Edge curEdg = EDG[indSortEdg++] ;
        uint16_t begT= VtoT[curEdg.Vbeg].numTree ;
        uint16_t endT =VtoT[curEdg.Vend].numTree ;
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
                uint16_t tmp = endT ;
                endT = begT ;
                begT= tmp ;
                tmp = curEdg.Vbeg ;
                curEdg.Vbeg = curEdg.Vend ;
                curEdg.Vend = tmp ;
            }
          // rattach beg tree to end tree
            uint16_t nxtEnd = VtoT[curEdg.Vend].nxtVertex ;
            uint16_t nxtBeg = curEdg.Vbeg ;
            VtoT[curEdg.Vend].nxtVertex = nxtBeg ;
            uint16_t antBeg ;
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
     uint64_t nMin = 1L << PB108_ALPHAM ;
    for(sumA=2;(1L<<sumA) <= nMin ;sumA++) {
        //      printf("\nSum(alpha)=%d ",sumA) ;
        int nbA,nbAmax ;
         uint64_t n=1 ;
        for(nbAmax=0;nbAmax<=sumA;nbAmax++) {
            n *= primes[nbAmax] ;
            if(n> nMin) break ;
        }
        uint64_t nbS = 1 ;
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
            uint64_t nbS = 1 ;
            for(ia=0;ia<nbA;ia++) {
                nbS *= 2*alpha[ia]+1 ;
            }
            uint64_t nbSol = (nbS+1)/2 ;
//            for(ia=0;ia<nbA;ia++) printf("%c%d",(ia==0)? '\n' : '.' ,alpha[ia]) ; printf(" nbS=%lld",nbSol) ;
           if(nbSol > minS ) {
                uint64_t n=1;
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

// le parcours des solutions peut se faire de la façon suivante
// Pour le dernier coup c'est un parmi les 21 doubles (D1..D20 + D25) => 21
// Pour les cas a 2 ou 1 coups : il faut pour le coup supllementaire un aprmi les 62 cas (S1..S20 S25,  D1..D20 +D25 T1..T20) + le cas vide
// => 21 x (62+1)
// Pour les cas a trois coup, pour ne pas doublonner les coups et 2 (Parmi les 62) il suffir de les classer donc 62+61+...+1 = (62x63)/2
// => 21 x 62x63)/2
// donc un total de 21 ( 63+62x63/2) = 21x63x(62+2)/2= 21x63x32.

int cmpVal(const void *e1,const void *e2) {
    return ((int *)e1)[0] - ((int *)e2)[0] ;
}
#define PB109_MAX   100
int PB109(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nb =0;
    int val[63],D[21];
    val[0] = 0 ;
    int i,j,k, is=1 ;
    for(i=1;i<=20;i++) {
        D[i-1] = 2*i ;
        val[is++]= i;
        val[is++] = 2*i ;
        val[is++] = 3*i ;
    }
    val[is++] = 25 ; val[is++] = 50 ; D[20] = 50 ;
    qsort(val,63,sizeof(val[0]),cmpVal) ;
    for(i=0;i<21;i++) {
        for(j=0;j<63;j++) {
            if(D[i]+val[j] >= PB109_MAX) break ;
            for(k=j;k<63;k++) {
                if(D[i]+val[j]+val[k] < PB109_MAX) nb++ ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nb) ;
    return 1 ;
}

#define PB111_MAXP  400000
#define PB111_NBD   10
int PB111(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB111_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    uint32_t pow10[PB111_NBD+1] = { 1 } ;
    int nd,k,id;
    pow10[PB111_NBD] = pow10[0] ;
    for(k=1;k<PB111_NBD;k++) {
        pow10[k] = pow10[k-1] * 10 ;
        pow10[PB111_NBD] += pow10[k] ;
    }
    uint64_t S ;
    int nbDOk = 0 ;
    S = 0 ;
    int nbD[10] = {0,0,0,0,0,0,0,0,0,0} ;
    for(nd=1;nbDOk < 10 ;nd++) {
        uint8_t noD[PB111_NBD] ;
        uint8_t valNoD[PB111_NBD] ;
        for(k=0;k<nd;k++) {noD[k]=k ; valNoD[k]= 0 ;}
        do {
            uint32_t Ndig = pow10[PB111_NBD] ;
            for(k=0;k<nd;k++) {
                Ndig -= pow10[noD[k]] ;
            }
            int is = nd-1 ;
            while (is >=0) {
                uint64_t NnoDig = 0 ;
                for(k=0;k<nd;k++) NnoDig += (uint64_t) pow10[noD[k]] * valNoD[k] ;
                {
                    for(id=0;id<10;id++) {
                        if(nbD[id] && nbD[id] < nd) continue ;
                        uint64_t N = (uint64_t) Ndig * id + NnoDig ;
                        if(N>pow10[PB111_NBD-1] && Is_Prime(N,tbPrime)) {
                            if(nbD[id] == 0) {
                                nbD[id] = nd ;
                                nbDOk++ ;
                            }
                            S += N ;
                        }
                    }
                }
                while(is >= 0 && valNoD[is] == 9) {
                    valNoD[is]= 0 ; is-- ;
                }
                if(is >= 0) { valNoD[is]++ ; is = nd-1 ; }
            }
        } while(NextSub(noD,nd,PB111_NBD)>=0 );
    }
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",S) ;
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
    printf(" %s%s => %6lld I +%6lld D + %10lld B T=%10lld Cum: %6lld C %6lld I +%6lld D +%10lld B T=%10lld Perc=%.6lf\n"
           ,info,joker+17-nbDigJoker,I,D,B,T,SC,SI,SD,SB,ST, ((double)100.0*SB)/ST) ;
#else
    printf(" %s%s =>  %6lld C %6lld I +%6lld D +%10lld B T=%10lld Perc=%.8lf\n"
           ,info,joker+17-nbDigJoker,SC,SI,SD,SB,ST, (((double)100.0)*SB)/ST);
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


#define PB113_NBDIG    100
// #define PB113_DEBUG
// struct to count the different number categories
// for increasing and decreasing number differentiation of the numbers
// by the leading digit.
// Also, the "constant numbers" are not included in incr[] and decr[]
// to avoid multiples counts.
// Leading '0' is only counted for recursion.
typedef struct CountNoB {
    int64_t Incr[10] ;  // increasing by leading digit (constant excluded)
    int64_t Decr[10] ;  // increasing by leading digit (constant excluded)
    int64_t Const ;     // Constant numbers (same digit)
    // cumulative counters for recursion
    int64_t SumI ;
    int64_t SumD ;
    int64_t SumC ;
} CountNoB ;

// info is only a comment
void CountNB_print(CountNoB *CNB,char *info) {
    int64_t SI = CNB->SumI;
    int64_t SD = CNB->SumD ;
    int64_t SC = CNB->SumC ;
 #if defined(PB113_DEBUG)
    int j ;
    int64_t I = CNB->Const  ; for(j=1;j<10;j++) I += CNB->Incr[j] ;
    int64_t D = CNB->Const  ; for(j=1;j<10;j++) D += CNB->Decr[j] ;
    int64_t C = CNB->Const ;
    
    printf(" %s =>  %6lld C %6lld I +%6lld D T=%6lld Cum: %6lld C %6lld I +%6lld D  NB=%10lld\n"
           ,info,C,I,D,I+D-C,SC,SI,SD,SI+SD-SC) ;
#else
    printf(" %s =>  %6lld C %6lld I +%6lld D\n"
           ,info,SC,SI,SD);
#endif
}

// initialise newC with oldC
// reset the current counters, tansmit cumulatives counters
void CountNB_Init(CountNoB *oldC, CountNoB *newC) {
    newC->SumI = oldC->SumI ;
    newC->SumD = oldC->SumD ;
    newC->SumC = oldC->SumC ;
    memset(newC->Incr,0,sizeof(newC->Incr));
    memset(newC->Decr,0,sizeof(newC->Decr));
    newC->Const = 0 ;
}

#define BIT_INC 1   // increasing
#define BIT_DEC 2   // decreasing
#define BIT_CONST 4 // increasing and decreasing

// add leading digits chain to statistics
//
void CountNB_AddHead(char *digits, CountNoB *oldC, CountNoB *newC) {
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
    if(oldC->Const) {
        for(ida=0;ida<10;ida++) {
            if(ida < idT) {  // T > A
                if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[ida] + 1;
            }else if (ida > idT) {  // T < A
                if(status & BIT_INC)  newC->Incr[idH] += oldC->Incr[ida]  + 1 ;
            } else { // T == A
                if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[ida] + ((status & BIT_CONST) ? 0 : 1)  ;
                if(status & BIT_INC) newC->Incr[idH] += oldC->Incr[ida]  + ((status & BIT_CONST) ? 0 : 1) ;
            }
        }
    } else { // special case for one digit numbers  to avoid problems with constant numbers
        if(status & BIT_DEC) newC->Decr[idH] += oldC->Decr[idT] + ((status & BIT_CONST) ? 0 : 1)  ;
        if(status & BIT_INC) newC->Incr[idH] += oldC->Incr[idT]  + ((status & BIT_CONST) ? 0 : 1) ;
    }
    if(idH!=0) { // accumulation
        newC->SumI += newC->Incr[idH] + ((status & BIT_CONST) ? 1  : 0) ;
        newC->SumD += newC->Decr[idH] + ((status & BIT_CONST) ? 1  : 0) ;
        if(status & BIT_CONST) {
            newC->Const++ ;
            newC->SumC++ ;
        }
    }
}


int PB113(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    char prefix[10] ;
    int id=1,nd ;
    CountNoB CNBwork ; // 2 counts so when depassing the purpose, rolling back one step
    CountNoB CNBnext ;
    memset(&CNBwork,0,sizeof(CNBwork)) ;
    CountNB_Init(&CNBwork,&CNBnext) ;
    // search the power of 10 reaching a superior proportion of bouncy
    for(nd=1;nd<=PB113_NBDIG;nd++) {
        for(id=0;id<10;id++) { // add the
            sprintf(prefix,"%d",id) ;
            CountNB_AddHead(prefix,&CNBwork,&CNBnext) ;
        }
#if defined(PB113_DEBUG)
        char info[30] ;
        sprintf(info,"under 10**%d",nd) ;
        CountNB_print(&CNBnext,info) ;
#endif
        
        CNBwork = CNBnext ;
        CountNB_Init(&CNBwork,&CNBnext) ;
    }
    nd-- ;
    int64_t noBouncy = CNBwork.SumI+CNBwork.SumD-CNBwork.SumC ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Under 10**%d not bouncy=%lld (+I=%lld +D=%lld -C=%lld)\n",
                               pbR->ident,nd,noBouncy,CNBwork.SumI,CNBwork.SumD,CNBwork.SumC) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",noBouncy) ;
    return 1 ;
}


#define PB114_LEN   50
typedef struct State114 {
    uint64_t B ;
    uint64_t R1 ;
    uint64_t R2 ;
    uint64_t R3 ;
} State114;


int PB114(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    State114 St[PB114_LEN+1] ;
    // init
    St[0].B = 1 ;
    St[0].R1 = St[0].R2= St[0].R3 = 0 ;
    int nt;
    for(nt=1;nt<PB114_LEN+1 ;nt++ ) {
        St[nt].B = St[nt-1].B + St[nt-1].R3 ;
        St[nt].R2 = St[nt-1].R1 ;
        St[nt].R3 = St[nt-1].R3 + St[nt-1].R2 ;
        St[nt].R1 = St[nt-1].B ;
    }
    nt-- ;
    uint64_t nbStates = St[nt].B+St[nt].R3 ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NnState=%lld=(%lld-B %lld %lld %lld-R3)\n",
                               pbR->ident, nbStates,St[nt].B,   St[nt].R1,St[nt].R2,St[nt].R3) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbStates) ;
    return 1 ;
}

/* On consider la matrice M (matrice de transition)
 1 0 0 1
 1 0 0 0
 0 1 0 0
 0 0 1 1
 
                                 1
                                 0
                                 0
 Le resultat est (1 0 0 1 ) M**n 0
 
 le produit final est pour ne retenir que les B er R3
 
 */

#define PB114_SIZE  4
#define PB114_MAXPOW2    6

#define I2(i,j)   ((i)*PB114_SIZE+(j))
#define I3(ip,i,j)   ((ip)*PB114_SIZE*PB114_SIZE+(i)*PB114_SIZE+(j))

int PB114a(PB_RESULT *pbR) {
    int i,j,k,ip ;
    uint64_t M[PB114_SIZE*PB114_SIZE] ;
    uint64_t Mpow2[PB114_MAXPOW2*PB114_SIZE*PB114_SIZE] ;
    uint64_t V[PB114_SIZE] , antV[PB114_SIZE] ;
    memset(M,0,sizeof(M)) ;
    memset(V,0,sizeof(V)) ;
    V[0] = 1 ;
    M[I2(0,0)] = 1 ; M[I2(0,PB114_SIZE-1)] = 1 ;
    M[I2(PB114_SIZE-1,PB114_SIZE-1)] = 1 ;
    for(i=1;i<PB114_SIZE;i++)M[I2(i,i-1)] = 1 ;
    memcpy(Mpow2,M,sizeof(M)) ;
    for(ip=1;ip<PB114_MAXPOW2;ip++) {
        for(i=0;i<PB114_SIZE;i++) {
            for(j=0;j<PB114_SIZE;j++) {
                uint64_t S =0 ;
                for(k=0;k<PB114_SIZE;k++) {
                    S += Mpow2[I3(ip-1,i,k)] * Mpow2[I3(ip-1,k,j)] ;
                }
                Mpow2[I3(ip,i,j)] = S ;
            }
        }
    }
    for(ip=0;ip<PB114_MAXPOW2;ip++) {
        if((1<<ip) & PB114_LEN) {
            // puissance presente, on mutliplie V par M**(1<<ip))
            memcpy(antV,V,sizeof(V)) ;
            for(i=0;i<PB114_SIZE;i++) {
                uint64_t S =0 ;
                for(k=0;k<PB114_SIZE;k++) {
                    S += Mpow2[I3(ip,i,k)] * antV[k] ;
                }
                V[i] = S ;
            }
            
        }
    }
    uint64_t nbStates = V[0]+V[PB114_SIZE-1] ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NnState=%lld=(%lld-B  %lld-R3)\n",
                               pbR->ident, nbStates,V[0],V[PB114_SIZE-1]) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbStates) ;
    return 1 ;
    
    
}


#define PB115_SIZE  51
#define PB115_MAX   1000000
int PB115(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    uint64_t St0[PB115_SIZE],St1[PB115_SIZE] ;
    uint64_t * antSt, *curSt ;
    // init
    memset(St0,0,sizeof(St0)) ;
    St0[0] =  1 ;
    antSt = St0 ;
    curSt = St1 ;
    int nt;
    for(nt=1;;nt++) {
        int i ;
        curSt[0] = antSt[0] + antSt[PB115_SIZE-1] ;
        for(i=1;i<PB115_SIZE-1;i++) {
            curSt[i] = antSt[i-1] ;
        }
        curSt[PB115_SIZE-1] = antSt[PB115_SIZE-2] + antSt[PB115_SIZE-1] ;
        uint64_t nbStates = curSt[0] + curSt[PB115_SIZE-1] ;
        if(nbStates > PB115_MAX) break ;
        uint64_t *tmp = antSt ;
        antSt = curSt ;
        curSt = tmp ;
     }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s n=%d\n",
                               pbR->ident, nt) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nt) ;
    return 1 ;
}


typedef struct State116 {
    uint64_t B ;
    uint64_t T1 ;
    uint64_t T2 ;
    uint64_t T3 ;
    uint64_t T4 ;
} State116;


#define PB116_LEN 50
int PB116(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    State116 StR0 ,StR1 , *ptRcur = &StR0 , *ptRant = &StR1 ;
    State116 StG0 ,StG1 , *ptGcur = &StG0  , *ptGant = &StG1 ;
    State116 StB0 ,StB1 , *ptBcur = &StB0 , *ptBant = &StB1 ;
    // init
    StR0.B = 1 ;
    StR0.T1 = StR0.T2 = 0 ;
    StG0.B = 1 ;
    StG0.T1 = StG0.T2 = StG0.T3 = 0 ;
    StB0.B = 1 ;
    StB0.T1 = StB0.T2 = StB0.T3 = StB0.T4 = 0 ;
    
    int nt;
    for(nt=0;nt<PB116_LEN ;nt++ ) {
        State116 * tmp  = ptRant ;
        ptRant = ptRcur ;
        ptRcur = tmp ;
  
        tmp  = ptGant ;
        ptGant = ptGcur ;
        ptGcur = tmp ;

        tmp  = ptBant ;
        ptBant = ptBcur ;
        ptBcur = tmp ;

        
        ptRcur->B = ptRant->B + ptRant->T2 ;
        ptRcur->T1 = ptRant->B + ptRant->T2 ;
        ptRcur->T2 = ptRant->T1;

        ptGcur->B = ptGant->B + ptGant->T3 ;
        ptGcur->T1 = ptGant->B + ptGant->T3 ;
        ptGcur->T2 = ptGant->T1;
        ptGcur->T3 = ptGant->T2;

        ptBcur->B = ptBant->B + ptBant->T4 ;
        ptBcur->T1 = ptBant->B + ptBant->T4 ;
        ptBcur->T2 = ptBant->T1;
        ptBcur->T3 = ptBant->T2;
        ptBcur->T4 = ptBant->T3;
    }
    uint64_t nbStates = ptRcur->B  + ptRcur->T2 + ptGcur->B  + ptGcur->T3 + ptBcur->B  + ptBcur->T4 -3 ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NnState=%lld R=%lld G=%lld B=%lld\n",
                               pbR->ident, nbStates,ptRcur->B  + ptRcur->T2 -1 ,ptGcur->B  + ptGcur->T3 -1 , ptBcur->B  + ptBcur->T4 -1) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbStates) ;
    return 1 ;
}


typedef struct State117 {
    uint64_t B ;
    uint64_t R1 ;
    uint64_t R2 ;
    uint64_t G1 ;
    uint64_t G2 ;
    uint64_t G3 ;
    uint64_t B1 ;
    uint64_t B2 ;
    uint64_t B3 ;
    uint64_t B4 ;
} State117;


#define PB117_LEN 50
int PB117(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    State117 St0 ,St1 , *ptCur = &St0 , *ptAnt = &St1 ;
    // init
    memset(&St0,0,sizeof(St0)) ;
    St0.B = 1 ;
    
    int nt;
    for(nt=0;nt<PB117_LEN ;nt++ ) {
        State117 * tmp  = ptAnt ;
        ptAnt = ptCur ;
        ptCur = tmp ;
       
        
        ptCur->B = ptAnt->B + ptAnt->R2 + ptAnt->G3 + ptAnt->B4 ;
        ptCur->R1 = ptCur->G1 = ptCur->B1 = ptAnt->B + ptAnt->R2 + ptAnt->G3 + ptAnt->B4 ;
        ptCur->R2 = ptAnt->R1;
        ptCur->G2 = ptAnt->G1;
        ptCur->G3 = ptAnt->G2;
        ptCur->B2 = ptAnt->B1;
        ptCur->B3 = ptAnt->B2;
        ptCur->B4 = ptAnt->B3;
    }
    uint64_t nbStates = ptCur->B  + ptCur->R2 + ptCur->G3 + ptCur->B4  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s NnState=%lld \n",
                               pbR->ident, nbStates) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbStates) ;
    return 1 ;
}
#define PB118_MAXDPRIME 50000
#define PB118_MAXP  100000000
typedef struct DPrime {
    uint32_t P ;       // value
    uint32_t mask ;    // binary mask for used digits
} DPrime ;

int PB118(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;

    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(Sqrt32(PB118_MAXP))) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);

    DPrime *tbDigPrime = malloc(PB118_MAXDPRIME*sizeof(tbDigPrime[0])) ;
    Decomp * DeC = DecompAlloc(9) ;
    int nbSub ;
    int indDPrime[10] ;
    int nbDPrime[10] ;
    int id = 1 ;
    int ip = 0 ;
    // primes with one digit
    indDPrime[id] = ip ;
    tbDigPrime[ip].P = 2 ;
    tbDigPrime[ip++].mask = 1<<(2-1) ;
    tbDigPrime[ip].P = 3 ;
    tbDigPrime[ip++].mask = 1<<(3-1) ;
    tbDigPrime[ip].P = 5 ;
    tbDigPrime[ip++].mask = 1<<(5-1) ;
    tbDigPrime[ip].P = 7 ;
    tbDigPrime[ip++].mask = 1<<(7-1) ;
    nbDPrime[id] = ip - indDPrime[id] ;
    for(id=2;id<8;id++) {
        // primes with (id+1) different digits
        uint8_t sub[9] ;
        int k ;
        indDPrime[id] = ip ;
        for(k=0;k<9;k++) sub[k] = k ;
        do {
            if(sub[id-1] & 1) continue ; // avoid even number
            int N=0 ;
            for(k=0;k<id;k++) N = 10*N + sub[k] + 1 ;
            if(Is_Prime32(N,tbPrime)) {
                tbDigPrime[ip].P = N ;
                int mask = 0 ;
                for(k=0;k<id;k++) mask |= 1 << sub[k] ;
                tbDigPrime[ip++].mask = mask  ;
            }
        } while (NextArrangement(sub,id,9) >= 0) ;
        nbDPrime[id] = ip - indDPrime[id] ;
    }
    { // case 8 digits (id=7) : imposes to the excluded digit to be prime (so a solution)
        uint8_t sub[9] ;
        int k ;
        indDPrime[id] = ip ;
        for(k=0;k<8;k++) sub[k] = k ;
        do {
            int ipExclus = 1 ; // minus one, 2,3,5,7
            while(ipExclus) {
                int N=0 ;
                for(k=0;k<id;k++) {
                    N = 10*N + ((sub[k] == ipExclus) ? 9 :  (sub[k] + 1) )  ;
                }
                if(Is_Prime32(N,tbPrime)) {
                    tbDigPrime[ip].P = N ;
                    int mask = 0 ;
                    for(k=0;k<id;k++) mask |= 1 << ((sub[k] == ipExclus) ? 8 :  (sub[k]) )  ;
                    tbDigPrime[ip++].mask = mask  ;
                }
                if(ipExclus == 1) ipExclus = 2;
                else if (ipExclus == 2) ipExclus = 4;
                else if (ipExclus == 4) ipExclus = 6 ;
                else ipExclus = 0 ;
            }
        } while (NextPermut(sub,8) >= 0) ;
        nbDPrime[id] = ip - indDPrime[id] ;
        id++ ;
    }
    
 
    //no solution as permutation of nine digits is a multiple of 3
    indDPrime[id] = ip ;
    nbDPrime[id] = 0 ;
    nbSub = nbDPrime[8] ; // as 8 digits are solution
    if(pbR->isVerbose) {
            printf("\t PB%s +%-5d form:8.1\n",pbR->ident,nbSub) ;
    }
    do { // loop on decomposition of 9 as a sum (number of digits of eeach number)
        int j ;
        int indFree[9] ;
        indFree[0] = 0 ;
        int is = 0 ;
        int mask = 0 ;
        int oldNbSub = nbSub ;
        if(DeC->val[0] > 7) continue ; // 9 and 8 digits
        while(is>=0) { // loop to fill the decomposition (tree exploration)
            if(indFree[is] < nbDPrime[DeC->val[is]]) {
                if(mask & tbDigPrime[indDPrime[DeC->val[is]]+indFree[is]].mask ) {
                    indFree[is]++ ; continue ; // mask conflit ?
                }
                mask |= tbDigPrime[indDPrime[DeC->val[is]]+indFree[is]].mask ;
                indFree[is]++ ;
                if(is == DeC->nbVal - 1) {
                    nbSub++ ; // remove last digit in the mask
                    mask ^= tbDigPrime[indDPrime[DeC->val[is]]+indFree[is]-1].mask ;
                } else {
                    is++ ;
                    if(DeC->val[is] == DeC->val[is-1]) { // same length, must be different
                        indFree[is] = indFree[is-1] ;
                    } else {
                        indFree[is] = 0 ;
                    }
                }
            } else { // no more candidate
                indFree[is] = 0 ;
                is-- ;
                if(is>=0) mask ^= tbDigPrime[indDPrime[DeC->val[is]]+indFree[is]-1].mask ;
            }
        }
        if(pbR->isVerbose) {
            if(nbSub > oldNbSub) {
                printf("\t PB%s +%-5d form:",pbR->ident,nbSub-oldNbSub) ;
                for(j=0;j<DeC->nbVal;j++) printf("%d%c",DeC->val[j],(j== DeC->nbVal -1) ? '\n' : '.');
            }
        }
    } while(DecompNext(DeC) >= 0) ;
    DecompFree(DeC) ;
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbSub) ;
    return 1 ;
}

#define PB119_MAX   30
#define PB119_NBD_MAX 20
int PB119(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nd,nb = 0 ;
    uint64_t s ;
    uint64_t Pow[PB119_NBD_MAX*9];
    uint64_t Pow10 ;
    uint64_t Pfind[2*PB119_MAX] ;
    for(s=2;s<PB119_NBD_MAX*9;s++){
        Pow[s] = s*s ;
    }
    int antNb = 0 ;
    for(nd=3,Pow10=1000;nb<PB119_MAX;nd++ , Pow10 *= 10) {
        antNb = nb ;
        for(s=2;s<=9*PB119_NBD_MAX;s++){
            uint64_t pow  ;
            for(pow=Pow[s];pow < Pow10 ; pow *= s){
                uint64_t p1 = pow  ;
                int s1 = 0 , np =0 ;
                while(p1) {
                    s1 += p1 % 10 ;
                    np++ ;
                    if(s1 > s) break ;
                    p1 /= 10 ;
                }
                if(s1 == s) {
                    Pfind[++nb] = pow  ;
                    if(pbR->isVerbose) fprintf(stdout,"\t PB%s %lld\t=\t%lld**%d\n",pbR->ident,pow,s,np) ;
                }
            }
            Pow[s] = pow ;
            
        }
    }
    while(antNb < PB119_MAX) {
        int in ;
        antNb++ ;
        for(in=antNb+1; in<=nb;in++) {
            if(Pfind[antNb] > Pfind[in]) {
                uint64_t tmp = Pfind[antNb] ;
                Pfind[antNb] = Pfind[in] ;
                Pfind[in] = tmp ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s The %dth is %lld\n",pbR->ident,antNb,Pfind[antNb]) ;
    sprintf(pbR->strRes,"%lld",Pfind[antNb]) ;
    return 1 ;
}


int PB120(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int SumR = 0 ;   
    SumR += 2*3 ;
    int n4 ;
    for(n4=4;n4<1000;n4 += 4) {
        SumR += n4 *(n4-2) + (n4+1)*(n4) + (n4+2)*(n4) + (n4+3)*(n4+2) ;
    }
    SumR += n4 * (n4-2) ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum of rest=%d\n",pbR->ident,SumR) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",SumR) ;
    return 1 ;
}

#define PB121_NBTURN    15
// On calcule par recurrence Pn(x) = Sigma(P(k,n) x**k ) avec P(k,n) = nb cas avec k bleu
// on a la recurrence Pn(x) = Pn-1(x) * (x+n)
// le nombre total  de cas est Factoriel(n+1)
// Il restec ensuite a sommer les cas favorables.
int PB121(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    
    uint64_t Cf[PB121_NBTURN+1] ;
    uint64_t fact = 2;
    Cf[0] = 1 ;
    Cf[1] = 1 ;
    int n,k ;
    // Pn(x) = Pn-1(x) * (x+n)
    for(n=2;n<=PB121_NBTURN;n++) {
        fact *= n+1 ;
        uint64_t ck = Cf[0] ;
        Cf[0] = ck * n ;
        Cf[n]= 0 ;
        for(k=1;k<=n;k++) {
            uint64_t tmp = ck + n * Cf[k] ;
            ck = Cf[k] ;
            Cf[k] = tmp ;
        }
    }
    uint64_t sumCkMin = 0 ;
    for(k=PB121_NBTURN;k>PB121_NBTURN/2;k--) {
        sumCkMin += Cf[k] ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum=%lld (%lld/%lld)\n",pbR->ident,fact/sumCkMin,fact,sumCkMin) ;
    sprintf(pbR->strRes,"%d",(int) (fact/sumCkMin)) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB122_MAX   100000
// #define PB122_MAX   15000
#define PB122_NBCHAIN   2000000000
#define PB122_CHK 1
#define PB122_MAXLEVEL  64
typedef struct Chain122 {
    uint32_t   n ;
    uint32_t   antCh ;
} Chain122 ;
int PB122(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int SumR = 0 ;
    
    Chain122 * tbCH = malloc(sizeof(tbCH[0])*PB122_NBCHAIN) ;
    uint8_t MinM[PB122_MAX+1] ;
    int nxtCH = 0 ;
    int antCH = 0 ;
    int nbMulMax = 1 ;
    int nbFind = 0 ;
    memset(MinM,0,sizeof(MinM)) ;
    while((1<<nbMulMax) <= PB122_MAX) nbMulMax++ ;
    nbMulMax = 2 * (nbMulMax - 1 ) ; // decomposition binaire (calculer les puissances, les sommer)
    MinM[2] = 1 ;
    nbFind = 2 ;
    int nm ;
    tbCH[nxtCH].n = 1 ;
    tbCH[nxtCH].antCh = nxtCH ;
    nxtCH++ ;
    for(nm=1;nbFind < PB122_MAX;nm++) {
        int k ;
        int curCH = nxtCH ;
        for(k=antCH;k<curCH;k++) {
            Chain122 endCH = tbCH[k] ;
            int  iCh = k ;
            while(1) {
                 int n = endCH.n + tbCH[iCh].n ;
                if ((n <= PB122_MAX) && (MinM[n] == 0 || MinM[n] >= nm)) {
 //                   printf("(%d->%d)",n,endCH.n) ;
                    tbCH[nxtCH].n = n ;
                    tbCH[nxtCH++].antCh = k ;
                    if(MinM[n] == 0) {
                        MinM[n] = nm ;
                        nbFind++ ;
                    }
                    if(iCh==0 && endCH.antCh > 0) {
                        tbCH[nxtCH].n = endCH.n - tbCH[endCH.antCh].n +1 ;
                        tbCH[nxtCH++].antCh =  tbCH[endCH.antCh].antCh ;
                        tbCH[nxtCH].n = n ;
                        tbCH[nxtCH].antCh = nxtCH-1 ; nxtCH++ ;
                    }
                }
                if(iCh) {
                    iCh = tbCH[iCh].antCh ;
                } else {
                    break ;
                }
            } ;
        }
        printf("%d,%d ",nxtCH,nbFind) ;
        antCH = curCH ;
    }
    int i,sumM = 0 ;
    for(i=1;i<=PB122_MAX;i++) sumM +=  MinM[i] ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum([1 %d]=%d, %d nodes Max=%d\n",pbR->ident,PB122_MAX,sumM,nxtCH,nm) ;
    sprintf(pbR->strRes,"%d",sumM) ;
    free(tbCH);
    pbR->nbClock = clock() - pbR->nbClock ;
#if PB122_CHK
    {
        const int32_t * MinRef = P122_GetData() ;
        int nbErr = 0 ;
        int iMax = ( PB122_MAX < P122_REFLG ) ? PB122_MAX : P122_REFLG ;
        for(i=1;i<=iMax;i++) {
            if(MinM[i] != MinRef[i]) {
                nbErr ++;
                printf("([%d]=%d exp=%d)\n",i,MinM[i],MinRef[i]);
            }
                
        }
        printf("CHECK again http://oeis.org/A003313/b003313.txt %d errors\n",nbErr);
    }
#endif
    
    return 1 ;
}


int PB122a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    Chain122 tbCH[30] ;
    int isDeltaOne[30] ;
    uint8_t MinM[PB122_MAX+1] ;
    int nbMulMax = 1 ;
    int nbFind = 0 ;
    int iOut= 0 ;
 //   memset(MinM,0,sizeof(MinM)) ;
    int i ;
    for(i=0;i<PB122_MAX+1;i++) MinM[i] = 0 ;
    while((1<<nbMulMax) <= PB122_MAX) nbMulMax++ ;
    nbMulMax += 2 ;
    //    nbMulMax = 2 * (nbMulMax - 1 )-1 ; // decomposition binaire (calculer les puissances, les sommer)
    MinM[1] = 0 ;
    MinM[2] = 1 ;
    nbFind = 2 ;
    for(;nbFind<PB122_MAX;nbMulMax++) {
        int is = 1 ;
        tbCH[is].n = 1 ;
        tbCH[is].antCh = 1 ;
        while(is > 0) {
           if(is>=nbMulMax) { is-- ; continue ;}
           if(tbCH[is].antCh < 1) {
               is-- ; continue;
           }
          int n = tbCH[is].n + tbCH[tbCH[is].antCh].n ;
//          iOut++ ; printf("%d(%d)=%d+%d%c",n,is,tbCH[is].n,tbCH[tbCH[is].antCh].n,(iOut & 7) ? ' ' : '\n') ;
          tbCH[is].antCh-- ;
            if(n <= PB122_MAX) {
               if(MinM[n] == 0 || MinM[n] >= is ) {
                    if(MinM[n] == 0) {
                        nbFind++ ;
                        if(nbFind >= PB122_MAX) {
                            MinM[n] = is ;
                            break ;
                        }
                    }
                    MinM[n] = is ;
                    if(n < PB122_MAX){
                        tbCH[++is].n = n ;
                        tbCH[is].antCh = is ;
                    }
                }
            }
        }
        printf("%d->%d ",nbMulMax,nbFind);
    }
    int sumM = 0 ;
    for(i=1;i<=PB122_MAX;i++) sumM +=  MinM[i] ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum([1 %d]=%d \n",pbR->ident,PB122_MAX,sumM) ;
    sprintf(pbR->strRes,"%d",sumM) ;
    pbR->nbClock = clock() - pbR->nbClock ;
#if PB122_CHK
    {
        const int32_t * MinRef = P122_GetData() ;
        int nbErr = 0 ;
        int iMax = ( PB122_MAX < P122_REFLG ) ? PB122_MAX : P122_REFLG ;
        for(i=1;i<=iMax;i++) {
            if(MinM[i] != MinRef[i]) {
                nbErr ++;
                printf("([%d]=%d exp=%d)\n",i,MinM[i],MinRef[i]);
            }
            
        }
        printf("CHECK again http://oeis.org/A003313/b003313.txt %d errors\n",nbErr);
    }
#endif
    return 1 ;
}

typedef struct Chain122b {
    uint32_t   n ;
//    uint32_t   n2 ;
    int32_t   antCh1 ;
    int32_t   antCh2 ;
} Chain122b ;

int PB122b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    Chain122b tbCH[30] ;
    uint8_t MinM[PB122_MAX+1] ;
    int nbMulMax = 1 ;
    int nbFind = 0 ;
    int iOut= 0 ;
    //   memset(MinM,0,sizeof(MinM)) ;
    int i ;
    for(i=0;i<PB122_MAX+1;i++) MinM[i] = 0 ;
    while((1<<nbMulMax) <= PB122_MAX) nbMulMax++ ;
    nbMulMax += 1 ;
    //    nbMulMax = 2 * (nbMulMax - 1 )-1 ; // decomposition binaire (calculer les puissances, les sommer)
    MinM[1] = 0 ;
    MinM[2] = 1 ;
    nbFind = 2 ;
    
    MinM[1] = 0 ;
    MinM[2] = 1 ;
    nbFind = 2 ;
    for(;nbFind<PB122_MAX;nbMulMax++) {
        int is = 1 ;
        tbCH[is].n = 1 ;
        tbCH[is].antCh1 = 1 ;
        while(is > 0) {
            if(is>=nbMulMax) { is-- ; continue ;}
            if(tbCH[is].antCh1 < 1) {
                is-- ; continue;
            }
            int n = tbCH[is].n + tbCH[tbCH[is].antCh1].n ;
            //          iOut++ ; printf("%d(%d)=%d+%d%c",n,is,tbCH[is].n,tbCH[tbCH[is].antCh].n,(iOut & 7) ? ' ' : '\n') ;
            tbCH[is].antCh1-- ;
            if(n <= PB122_MAX) {
                if(MinM[n] == 0 || MinM[n] >= is ) {
                    if(MinM[n] == 0) {
                        nbFind++ ;
                        if(nbFind >= PB122_MAX) {
                            MinM[n] = is ;
                            break ;
                        }
                    }
                    MinM[n] = is ;
                    if(n < PB122_MAX){
                        tbCH[++is].n = n ;
                        tbCH[is].antCh1 = is ;
                    }
                }
            }
        }
        printf("%d->%d ",nbMulMax,nbFind);
    }
    
    
    nbMulMax -=2 ;
    
//    for(;nbFind<PB122_MAX;nbMulMax++)
    {
        tbCH[0].n = 0 ;
        tbCH[0].antCh1 = tbCH[0].antCh2 = 0 ;
        int is = 1 ;
        tbCH[is].n = 1 ;
        tbCH[is].antCh1 = is ;
        tbCH[is].antCh2 = 0 ;
        int deltaIs = 0 ;
//        tbCH[is].n2 = tbCH[is].n ;
        while(is > 0) {
            if(is>=nbMulMax) { is-- ; continue ;}
            if(tbCH[is].antCh1 < tbCH[is].antCh2 ) {
                if(tbCH[is].antCh2 < is-1 ) {
//                if(tbCH[is].antCh2 < is-1 && tbCH[is].antCh2 < 2 ) {
//                if(tbCH[is].antCh2 < is-1 && tbCH[is].antCh2 < 100 ) {
                    tbCH[is].antCh1 = is ;
                    tbCH[is].antCh2++ ;
//                    tbCH[is].n2 = tbCH[is].n + tbCH[tbCH[is].antCh2].n  ;
                } else {
                    is-- ; continue;
                }
            }
//            int n = tbCH[is].n+tbCH[tbCH[is].antCh1].n + tbCH[tbCH[is].antCh2].n ;
            int n1 = tbCH[tbCH[is].antCh2].n + tbCH[tbCH[is].antCh1].n  ;
            deltaIs = tbCH[is].antCh2 ? 1 : 0 ;
            tbCH[is].antCh1-- ;
            int n = tbCH[is].n + n1  ;
//            iOut++ ; if(n<=PB122_MAX)printf("%d(%d)=%d+%d+%d%c",n,is,tbCH[is].n,tbCH[tbCH[is].antCh1+1].n,tbCH[tbCH[is].antCh2].n,(iOut & 7) ? ' ' : '\n') ;
            if(n <= PB122_MAX) {
//                if((MinM[n] >= is && tbCH[is].antCh2 == 0) || MinM[n] >= is+1  ) {
                if(MinM[n] >= is + deltaIs ) {
//                    if(tbCH[is].antCh2 == 0) {
                    if(deltaIs == 0) {
                        if(is < MinM[n]) {
                            MinM[n] = is ;
                            printf("\n%d(%d)->",n,is);
                            int j ;
                            for(j=is;j>0;j--) printf("%d[%d] ",tbCH[j].n,tbCH[j].antCh1+1) ;
                        }
                    } else {
                        if(is+1 < MinM[n]) {
                            MinM[n] = is+1 ;
                            printf("\n%d(%d)->*%d %d ",n,is+1,tbCH[tbCH[is].antCh1+1].n, tbCH[tbCH[is].antCh2].n );
                            int j ;
                            for(j=is;j>0;j--) printf("%d(%d)[%d,%d] ",tbCH[j].n,MinM[tbCH[j].n],tbCH[is].antCh1+1,tbCH[is].antCh2) ;
                        }
                    }
                    if(n < PB122_MAX){
//                        if(tbCH[is].antCh2 == 0) {
                        if(deltaIs == 0) {
                            tbCH[++is].n = n ;
                            tbCH[is].antCh1 = is ;
                            tbCH[is].antCh2 = 0 ;
//                            tbCH[is].n2 = n ;
                        } else {
                            tbCH[is+1].n = tbCH[tbCH[is].antCh1+1].n + tbCH[tbCH[is].antCh2].n ;
                            tbCH[is+1].antCh1 = is-1 ;
                            tbCH[is+1].antCh2 = is ;
                            is++ ;
  
                            tbCH[++is].n = n ;
                            tbCH[is].antCh1 = is ;
//                            tbCH[is].n2 = n ;
                            tbCH[is].antCh2 = 0 ;
                        }
                    }

                }
            }
        }
        printf("%d->%d ",nbMulMax,nbFind);
    }
    int sumM = 0 ;
    for(i=1;i<=PB122_MAX;i++) sumM +=  MinM[i] ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum([1 %d]=%d \n",pbR->ident,PB122_MAX,sumM) ;
    sprintf(pbR->strRes,"%d",sumM) ;
    pbR->nbClock = clock() - pbR->nbClock ;
#if PB122_CHK
    {
        const int32_t * MinRef = P122_GetData() ;
        int nbErr = 0 ;
        int iMax = ( PB122_MAX < P122_REFLG ) ? PB122_MAX : P122_REFLG ;
        for(i=1;i<=iMax;i++) {
            if(MinM[i] != MinRef[i]) {
                nbErr ++;
                printf("([%d]=%d exp=%d)\n",i,MinM[i],MinRef[i]);
            }
            
        }
        printf("CHECK again http://oeis.org/A003313/b003313.txt %d errors\n",nbErr);
    }
#endif
    return 1 ;
}

int PB122c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    Chain122b tbCH[30] ;
    uint8_t MinM[PB122_MAX+1] ;
    int nbMulMax = 1 ;
    int nbFind = 0 ;
    int iOut= 0 ;
    //   memset(MinM,0,sizeof(MinM)) ;
    int i ;
    for(i=0;i<PB122_MAX+1;i++) MinM[i] = 0 ;
    while((1<<nbMulMax) <= PB122_MAX) nbMulMax++ ;
    nbMulMax += 1 ;
    //    nbMulMax = 2 * (nbMulMax - 1 )-1 ; // decomposition binaire (calculer les puissances, les sommer)
    MinM[1] = 0 ;
    MinM[2] = 1 ;
    nbFind = 2 ;
    
    MinM[1] = 0 ;
    MinM[2] = 1 ;
    nbFind = 2 ;
    for(;nbFind<PB122_MAX;nbMulMax++) {
        int is = 1 ;
        tbCH[is].n = 1 ;
        tbCH[is].antCh1 = 1 ;
        while(is > 0) {
            if(is>=nbMulMax) { is-- ; continue ;}
            if(tbCH[is].antCh1 < 1) {
                is-- ; continue;
            }
            int n = tbCH[is].n + tbCH[tbCH[is].antCh1].n ;
            //          iOut++ ; printf("%d(%d)=%d+%d%c",n,is,tbCH[is].n,tbCH[tbCH[is].antCh].n,(iOut & 7) ? ' ' : '\n') ;
            tbCH[is].antCh1-- ;
            if(n <= PB122_MAX) {
                if(MinM[n] == 0 || MinM[n] >= is ) {
                    if(MinM[n] == 0) {
                        nbFind++ ;
                        if(nbFind >= PB122_MAX) {
                            MinM[n] = is ;
                            break ;
                        }
                    }
                    MinM[n] = is ;
                    if(n < PB122_MAX){
                        tbCH[++is].n = n ;
                        tbCH[is].antCh1 = is ;
                    }
                }
            }
        }
        printf("%d->%d ",nbMulMax,nbFind);
    }
    
    if(74195 < PB122_MAX){int Chain74195[20] = {2,4,8,16,17,33,66,83,149,215,430,579,1158,2316,4632,9264,18528,37056,74112,74195} ;
        int n ; for(n=0;n<20;n++)printf("%cMin[%d]=%d",(n & 7 ) ? ',' : '\n',Chain74195[n],MinM[Chain74195[n]] );
                                printf("\n");
    }
    
    
    
    nbMulMax -=2 ;
    
    
    {
        int is = 1 ;
        tbCH[is].n = 1 ;
        tbCH[is].antCh1 = 1 ;
        while(is > 0) {
            if(is>=nbMulMax) { is-- ; continue ;}
            if(tbCH[is].antCh1 < 1) {
                is-- ; continue;
            }
            int n = tbCH[is].n + tbCH[tbCH[is].antCh1].n ;
            //          iOut++ ; printf("%d(%d)=%d+%d%c",n,is,tbCH[is].n,tbCH[tbCH[is].antCh].n,(iOut & 7) ? ' ' : '\n') ;
            tbCH[is].antCh1-- ;
            if(n <= PB122_MAX) {
                if(MinM[n] == 0 || MinM[n] >= is-1 ) {
                    if(MinM[n] == 0) {
                        nbFind++ ;
                        if(nbFind >= PB122_MAX) {
                            MinM[n] = is ;
                            break ;
                        }
                    }
                    if(MinM[n] == 0 || MinM[n] > is) {
                        MinM[n] = is ;
                        printf("\n%d(%d)->",n,is);
                        int j ;
                        for(j=is;j>0;j--) printf("%d(%d) ",tbCH[j].n,MinM[tbCH[j].n]) ;
                    }
                    if(n < PB122_MAX){
                        tbCH[++is].n = n ;
                        tbCH[is].antCh1 = is ;
                    }
                }
            }
        }
       
    }
    if(74195 < PB122_MAX){int Chain74195[20] = {2,4,8,16,17,33,66,83,149,215,430,579,1158,2316,4632,9264,18528,37056,74112,74195} ;
        int n ; for(n=0;n<20;n++)printf("%cMin[%d]=%d",(n & 7 ) ? ',' : '\n',Chain74195[n],MinM[Chain74195[n]] );
        printf("\n");
    }

    
    int sumM = 0 ;
    for(i=1;i<=PB122_MAX;i++) sumM +=  MinM[i] ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum([1 %d]=%d \n",pbR->ident,PB122_MAX,sumM) ;
    sprintf(pbR->strRes,"%d",sumM) ;
    pbR->nbClock = clock() - pbR->nbClock ;
#if PB122_CHK
    {
        const int32_t * MinRef = P122_GetData() ;
        int nbErr = 0 ;
        int iMax = ( PB122_MAX < P122_REFLG ) ? PB122_MAX : P122_REFLG ;
        for(i=1;i<=iMax;i++) {
            if(MinM[i] != MinRef[i]) {
                nbErr ++;
                printf("([%d]=%d exp=%d)\n",i,MinM[i],MinRef[i]);
            }
            
        }
        printf("CHECK again http://oeis.org/A003313/b003313.txt %d errors\n",nbErr);
    }
#endif
    return 1 ;
}

#define PB123_MAX 10000000000LL


typedef int(*TY_CPL_nxtPrime)(void *ctx,T_prime nxtPrime);

int CPL123_nxtPrime(void *ctx,T_prime nxtPrime) {
    int32_t * nbPrime = (int32_t *) ctx ;
    *nbPrime += 1 ;
    if(*nbPrime & 1) {
        return (2LL * *nbPrime * nxtPrime > PB123_MAX)  ? 0 : 1 ;
    } else return 1 ;
}


int PB123(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int32_t nbPrime = 0 ;
    int32_t Pmax = 1+(int32_t) (sqrt((double) PB123_MAX) * log((double) PB123_MAX)) ;
    FindPrime(Pmax,&nbPrime,CPL123_nxtPrime) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbPrime) ;
    return 1 ;
}

// for PB124 and PB127
int RadCmp(const void * e1,const  void * e2) {
    return ((int *)e1)[0] - ((int *)e2)[0] ;
}
typedef struct RAD2N {
    int32_t * rad2n ;
    int32_t * n2rad ;
} RAD2N ;

// maxn : max value for n
// maxr : max rad(n) orderer (inclusive rad2n[maxr] is OK)
RAD2N * Rad2nAlloc(int maxn,int maxr,int isOrdered) {
    RAD2N * rdn = calloc(1,sizeof(rdn[0])) ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(Sqrt32(maxn))) == NULL) {
        fprintf(stdout,"\tFail to alloc prime table\n");
        return NULL ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    rdn->rad2n = malloc((maxn+1)*sizeof(rdn->rad2n[0]));
    rdn->n2rad = calloc(maxn,sizeof(rdn->n2rad[0]));
    int i ;
    int rgN = 0 ;
    rdn->rad2n[rgN++] = 1 ;
    rdn->n2rad[1] = 1 ;
    int indRad = 2 ;
    //    n2rad[1] = 1 ;
    while (rgN <= maxr) {
        while(rdn->n2rad[indRad]) indRad++ ;
        int curRad = indRad++ ;
        rdn->n2rad[curRad] = curRad ;
        int fact[20] ;
        int nbFact = 0 ;
        int val ;
        for(i=0,val=curRad;val>1 && i < nbPrime ;i++) {
            if((val % tbPrime[i]) == 0) {
                fact[nbFact++] = tbPrime[i] ;
                val /= tbPrime[i] ;
            }
        }
        if(val > 1) fact[nbFact++] = val ;
        int nbSup = 0 ;
        rdn->rad2n[rgN + nbSup++] = curRad ;
        for(i=0;i<nbFact;i++) {
            int p = fact[i] ;
            int powp ;
            int pSup = 0 ;
            for(powp=p; curRad*powp <= maxn; powp *=p) {
                int j ;
                int maxRad = maxn / powp ;
                for(j=0;j<nbSup;j++) {
                    if(rdn->rad2n[rgN+j] <= maxRad) {
                        int n = rdn->rad2n[rgN+j]*powp ;
                        rdn->rad2n[rgN+nbSup+pSup] = n ;
                        rdn->n2rad[n] = curRad ;  // printf("+%d ",n) ;
                        pSup++ ;
                    }
                }
            }
            nbSup += pSup ;
        }
        
        if(isOrdered) {
            qsort(rdn->rad2n,nbSup,sizeof(rdn->rad2n[0]),RadCmp) ;
        }
        if(rgN+nbSup < maxr) {
            rgN += nbSup ;
            continue ;
        } else {
            break ;
        }
    }
    Free_tablePrime(ctxP);
    return rdn ;
}

RAD2N * Rad2nFree(RAD2N * r2n) {
    if(r2n != NULL) {
        free(r2n->n2rad) ;
        free(r2n->rad2n);
    }
    free(r2n) ;
    return NULL ;
}

#define PB124_RG    10000
#define PB124_MAXP   100000
#define PB124_MAXF  20
int cmpRad(const void * e1,const  void * e2) {
    return ((int *)e1)[0] - ((int *)e2)[0] ;
}
int PB124(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nAsk = 0 ;
    RAD2N * rdn = Rad2nAlloc(PB124_MAXP,PB124_RG, 1) ;
    nAsk = rdn->rad2n[PB124_RG] ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nAsk) ;
    return 1 ;
}

#define PB125_MAXD  8
int PB125(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int lg ;
    int pow10[9]={1,10,100,1000,10000,100000,1000000,10000000,100000000} ;
    int incr[4] ;
    int i0,i1,i2,i3 ;
    int nbSol  ;
    int64_t Sum  ;
    nbSol = 0 ;
    Sum = 0 ;
    for(lg=1;lg<=PB125_MAXD;lg++) {
        int hlg = lg/2 ;
        int i ;
        for(i=0;i<hlg;i++) {
            incr[i]=6*(pow10[i]+pow10[lg-i-1]) ;
        }
        if(lg&1) { incr[i] = 6*pow10[i] ; hlg++ ; }
        for(i=hlg;i<4;i++) incr[i] = 0 ;
        for(i0=1;i0<10;i0++) {
            for(i1= 0 ; i1<10 ; i1++) {
                for(i2= 0 ; i2<10 ; i2++) {
                    for(i3= 0 ; i3<10 ; i3++) {
                        int N6=incr[0]*i0+incr[1]*i1+incr[2]*i2+incr[3]*i3 ;
                        int n,c ;
                            for(n=2;c=(n+1)*(2*n+1),n*c<=N6 ;n++) {
                            if(N6 % n) continue ;
                            int nd = N6 / n ;
                            int b = 6*(n+1) ;
                            int delta = b*b - 4 * 6 * (c-nd) ;
                            int sqDelta = Sqrt32(delta) ;
                            if(sqDelta*sqDelta == delta) {
//                                int k = (-b+sqDelta)/12 ;
//                                printf("%d=[%d..%d] ",N6/6,k+1,k+n) ;
                                Sum += N6/6 ;
                                nbSol++ ;
                                break ;
                            }
                        }
                        if(hlg <4) break ;
                    }
                    if(hlg<3) break ;
                }
                if (hlg<2) break ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",Sum) ;
    return 1 ;
}


#define PB126_NB    1000
#define PB126_MAXV 10000

//
// la Taille de la couche n est
// L(n) = 2 * S2 + 4 * (n-1)(S1 + n-2)
// Avec S1 = a+b+c et S2 = a*b + b*c + a*c

int PB126(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int minLayerNB = PB126_MAXV ;
    int * histLayer = calloc(PB126_MAXV,sizeof(histLayer[0])) ;
    int a=1,b=1,c=1 ;
    int layer ;

    int S2 ;
    // parcours a>=b>=c
    // on controle a chaque etape que S2=L(1) est inferieur a PB126_MAXV
    for(a=1;2*a+1<PB126_MAXV;a++) {
        for(b=1;b<=a && (a*b+a+b < PB126_MAXV);b++) {
            for(c=1;c<=b && ((S2=a*b+b*c+a*c)<PB126_MAXV);c++){
                int S1 = a+b+c ;
                int n ;
                for(n=1;(layer = S2 + 2*(n-1)*(S1+n-2) ) <PB126_MAXV;n++) {
                    histLayer[layer]++ ;
                }                
            }
        }
    }
    for(layer=6;layer < PB126_MAXV;layer++) {
        if(histLayer[layer]== PB126_NB) break ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
/*    {
        int i ;
        for(i=10;i<PB126_MAXV;i++) printf("%c%d",((i % 10) == 0) ? '\n':' ',histLayer[i]);
    }
*/    if(layer < PB126_MAXV){
        sprintf(pbR->strRes,"%d",layer*2) ;
        return 1 ;
    } else {
        return 0 ;
    }
}

#define PB127_MAXP   120000
// #define PB127_MAXP   50
#define PB127_MAXF  20
int PB127(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    RAD2N * rdn=Rad2nAlloc(PB127_MAXP,PB127_MAXP,0) ;
    //boucle externe sur b
    int a,b,c ;
    int sum = 0 ;
    for(b=3;b<PB127_MAXP;b++) {
        int rb = rdn->n2rad[b] ;
        if(rb == b) continue ;
        int maxRa = PB127_MAXP / rb ;
        for(a=1;a<b;a++) {
           c = a+b ;
            if(c>PB127_MAXP) break ;
            int ra = rdn->n2rad[a] ;
            if(ra >= maxRa) continue ;
            int rc = rdn->n2rad[c] ;
            if((int64_t)ra*rb*rc < c ) {
                if(rdn->n2rad[ra*rb*rc]==ra*rb*rc){
//                       printf("rad(%dx%dx%d)=%d<%d\n",a,b,c,ra*rb*rc,c) ;
                    sum += c ;
                }
            }
        }
    }
    Rad2nFree(rdn) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",sum) ;
    return 1 ;
}

// variante more efficient
// extern loop for c
// change the condition a<b by rad(a)<rad(b) (equivalant as ra!=rb to remove double count)
// limit exploration as ra < sqrt(c/rc)
// explore a by rad(a) order (easy as obtain naturally by rad(n) calculation in rad2n
int PB127a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    RAD2N * rdn=Rad2nAlloc(PB127_MAXP,PB127_MAXP,0) ;
    int a,b,c ;
    int sum = 0 ;
    for(c=3;c<PB127_MAXP;c++) {
        int rc = rdn->n2rad[c] ;
        if(rc == c) continue ;
        int maxRab = c/rc ;
        int maxRa = Sqrt32(maxRab) ;
        int ia,ra ;
        for(ia=0;a=rdn->rad2n[ia], (ra=rdn->n2rad[a])<=maxRa;ia++) {
            b=c-a ;
            int rb = rdn->n2rad[b] ;
            if( ra*rb < maxRab && ra<rb && PGCD(ra,rb) == 1) {
//              printf("rad(%dx%dx%d)=%d<%d\n",a,b,c,ra*rb*rc,c) ;
                sum += c ;
            }
        }
     }
    Rad2nFree(rdn) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",sum) ;
    return 1 ;
}

// exceot for first and last tile of each ring, each tile has +1 and -1 for difference (ring neighbours)
// plus 2 neighbours of the same parity => maximum 2 prime differences.
// Only the first and the last tile of each ring can be a candidate
// For first tile of ring k, n=3*(int64_t)k*(k-1)+2 and TEST <=> IsPrime(6*k-1) && IsPrime(6*k+1) && IsPrime(12*k+5)
// For last tile of ring k, n = 3*(int64_t)k*(k+1)+1  and TEST <=> IsPrime(6*k-1) && IsPrime(6*k+5) && && IsPrime(12*k+7)
#define PB128_P  5000
#define PB128_ASK 2000
int PB128(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB128_P)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    uint8_t *PD = calloc(PB128_P,sizeof(PD[0])) ;
    int k,nFound=2;
    int64_t n=0 ;
    for(k=2;nFound<PB128_ASK;k++) {
        if(Is_Prime32(6*k-1,tbPrime)) {
            if(Is_Prime32(12*k+5,tbPrime) && Is_Prime32(6*k+1,tbPrime) ) {
                nFound++;
                if(nFound>=PB128_ASK) {
                    n = 3*(int64_t)k*(k-1)+2 ;
                    break ;
                }
            }
            if( Is_Prime32(12*k-7,tbPrime) && Is_Prime32(6*k+5,tbPrime) ) {
                nFound++; if(nFound>=PB128_ASK) {
                    n = 3*(int64_t)k*(k+1)+1 ;
                    break ;
                }
            }
        }
    }
    sprintf(pbR->strRes,"%lld",n) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
 }

//
// little more efficient version by checking that
// p = 6*k-1 prime First tile TST <=> IsPrime(p+2) && IsPrime(2*p+7)
//                 Last tile TST <=> IsPrime(p+6) && IsPrime(2*p-5)
#define PB128_Pa 500000
int PB128a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB128_Pa)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int64_t n = 0 ;
    int i,nFound=2;
    for(i=4;nFound<PB128_ASK && i<nbPrime-2;i++) {
        int p = tbPrime[i] ; // p = 6*k-1
        if((p % 6) == 5) {
            // Is_Prime32(12*k+5,tbPrime) && Is_Prime32(6*k+1,tbPrime)
            if(tbPrime[i+1] == p+2 && Is_Prime32(2*p+7,tbPrime)) {
                nFound++;
                if(nFound>=PB128_ASK) {
                    int k = (p+1)/6 ;
                    n = 3*(int64_t)k*(k-1)+2 ; break ;
                }
            }
//           (Is_Prime32(p+6,tbPrime) && Is_Prime32(2*p-5,tbPrime))
            if((tbPrime[i+1] == p+6 || tbPrime[i+2] == p+6) && Is_Prime32(2*p-5,tbPrime)) {
                nFound++;
                if(nFound>=PB128_ASK) {
                    int k = (p+1)/6 ;
                    n = 3*(int64_t)k*(k-1)+2 + 6*k - 1 ; break ;
                }
            }
        }
    }
    sprintf(pbR->strRes,"%lld",n) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



#define PB129_NB    1000000
#define PB129_PRIME 1100
// on joue sur le fait que R(n) = (10**n -1) / 9
//  si n % 3 != 0  et n premier avec 10.
//      Si n=p premier On a 10**(p-1)=1 Mod[p]
//      Il reste a verifier qu'il ny a pas un diviseur d de p-1 tel que 10**d=1 Mod[p]
//      Si n pas premier n=p1*p2*..
//      Alors 10**((p1-1)*(p2-2).. ) =1 Mod[n] et (p1-1)*(p2-1)... < 10000000 donc ne marche pas
// si n % 3 == 0
//      Si n/3-1 est divisible par 3 alors 10**(n/3-1)=1 Mod[3* (n/3)] ne marche pas
//      si n/3 pas premier ramene au cas 10**(3*(p1-1)*(p2-2).. ) =1 Mod[3*(n/3)]
//      Si n/3 est  premier on est ramene a peu pres au cas précédent
//          avac 10**[3*(n/3-1)] = 1 Mod[3*n/3] qui convient car 3*(n/3-1) tres proche de n

int PB129(PB_RESULT *pbR) {
    int64_t tbPow[32] ;
    int maxPow ;
    int isMult3=0 ;
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB129_PRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int n,n0;
    for(n0=(PB129_NB & 0x7ffffffe) + 1 ;;n0 += 2)  { // n impair
        int i,p ;
        isMult3=0 ;
        if((n0 % 3) == 0) {
            isMult3= 1 ;
            n = n0/3 ;
            if(((n-1) % 3) == 0 ) continue ;
        } else {
            n = n0 ;
        }
        for(i=0;p=tbPrime[i], p*p <=n ;i++) {
            if( (n%p)==0) break ;
        }
        if(p*p <=n ) continue ;
        maxPow = 0 ;
        tbPow[maxPow] = 10 ;
        if((n % 5) == 0  ) continue ;
        int d1, n1=n-1 ;
        int64_t pd1,pd2 ;
        for(d1=2;d1*d1<=n1;d1++) {
            if((n1 % d1) == 0) {
                int d2 = n1 / d1 ;
                int k ;
                for(k=0;(1<<k) <= d2;k++) ;
                k-- ;
                while(maxPow++< k) {
                    tbPow[maxPow] =  (tbPow[maxPow-1]* tbPow[maxPow-1]) % n ;
                }
                pd1=1 ;
                pd2=1 ;
                while(k>=0) {
                    if((1<<k) & d1) pd1 = (pd1 * tbPow[k]) % n ;
                    if((1<<k) & d2) pd2 = (pd2 * tbPow[k]) % n ;
                    k--;
                }
                if(pd1==1 || pd2==1) {
                 break ;
                }
            }
        }
        if(d1*d1 > n1) {
            break ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(n0 <  PB129_PRIME * PB129_PRIME) {
        if(isMult3) n *= 3  ;
        sprintf(pbR->strRes,"%d",n) ;
        return 1 ;
    } else {
        return 0 ;
    }
}

// return 10**pow mod[n]
//
int pow10_modn(int pow,int n) {
    int k,ln2 ;
    int64_t tbPow[32] ;
    for(ln2=0;(1<<ln2) <= pow;ln2++) ;
    ln2-- ;
    k = 0 ;
    tbPow[k] = 10 ;
    int p10 = (pow & 1) ? 10 : 1 ; //
    while(k++ < ln2) {
        tbPow[k] =  (tbPow[k-1]* tbPow[k-1]) % n ;
        if( (pow & (1<<k)) != 0) p10 = ( p10 * tbPow[k] ) % n ;
    }
    return p10 ;
}

#define PB130_PRIME 1000
#define PB130_NB  25

//#define PB130_PRIME 10000
//#define PB130_NB  1000

// on peut eliminer les multiples de 3
// si n = 3 * k alors n-1 pas multiple de 3 hors A[n] est multiple de 3.
// Donc on prend n=p1**k1 *p2**k2 *.. pj**kj (non premier avec pi != 2,3,5)
// A[n] est un diviseur de Q = PPCM { (p1-1)**k1 , (p2-2)**k2 , ... (pj-1)**kj }
// puis il faut on verifie que A[n] divise n-1, pour cela il suffit de verifier
// que 10** PGCD(Q,n-1) = 1 mod[n] ou encore pour 1..j que 10 * PGCD( (pi-1)**ki , n-1 ) = 1 mod[pi**ki]
int PB130(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB130_PRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int n ,nbfound=0 ;
    uint64_t sum = 0 ;
// debut par 49 car avant multiple de 3, ou 5 ou premier
    for(n=49 ;;n += 2)  { // n impair
        int i,p ;
        if( ((n % 3) == 0) || ((n % 5) == 0)) continue ;
        int n0 = n ;
        // on saute 2,3,5
        for(i=3;p=tbPrime[i], p*p <=n ;i++) {
           
            if((n0 % p) == 0) {
                int pe,p1 ;
                for(p1= p-1 ,pe= p; (n0 % (pe*p))== 0; p1 *= p-1 , pe *=p  ) ;
                int g = PGCD(n-1,p1) ;
                if( pow10_modn(g,pe)!=1) {
                    n0 = n ; break ;
                } else if((n0=n0/pe) ==1 ) break ;
                
            }
            
        }
        if(n0==1) {
            nbfound++ ; sum += n ;
//            if(pbR->isVerbose) fprintf(stdout," %d",n)  ;
        } else if(n0 < n ) { // dernier diviseur
        
            int g = PGCD(n-1,n0-1) ;
            if( pow10_modn(g,n0) == 1) {
                nbfound++ ;
                sum += n ;
 //               if(pbR->isVerbose) fprintf(stdout," %d",n)  ;
            }
        }
        if(nbfound>=PB130_NB) break ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Repunit[%d]=%d Sum=%lld\n",pbR->ident,PB130_NB,n,sum) ;

    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",sum) ;
    return 1 ;
}


// n**3+n**2 * p = m**3
// si pi**ai divise n (pi premier exposant ai)
// n**2 (n+p) est divisible par pi**(2xai) car n+p premier avec pi (sauf** si pi=p)
// donc ai doit etre multiple de 3 (puisque decomposition de m**3)
// donc n=q**3
// n**2 (n+p) = q**6 (q**3+p)
// donc q**3+p = r**3 (cube) <=> p = r**3-q**3 = (r-q) (r**2+rq+q**2)
// dc r=q+1 en remplacant => p = 3*q*(q+1) + 1
// donc les p sont tels que : p=1 mod[3] (p-1)/3 = q*(q+1) (unicite de q)
// **(pi=p) n = p**a0 ;   n**3+n**2 * p = q**3 p**(2.a0) ( p**a0 q**3 + p)
// = q**3 x p**(2.a0+1)  (p**(a0-1) x q**3 + 1)
// si a0=1 il faudrait q**3+1 egal a un cube. Impossible (cubes consecutifs)
// si a0=1+3a il faudrait (p**3a x q**3 +1) egal a un cube.  Impossible (cubes consecutifs)
// Les solutions sont p premier p-1 =  3 * p1 et p1 = q * (q+1)
// Comme les n/ln(n) > sqrt(n) on va parcourir les q(q+1) et regarder si 3*q*(q+1)+1 est premier
//#define PB131_MAX       1000000000000LL
//#define PB131_PRIME     1000000

#define PB131_MAX       1000000
#define PB131_PRIME     1000
int PB131(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB131_PRIME+1)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int64_t q,p1 ;
    int nb = 0 ;
    for(q=1;( p1=3*q*(q+1)+1) <= PB131_MAX ;q++)  { // n impair
        int i;
        int64_t p ;
        // on saute 2
        for(i=1; (i< nbPrime) && ( p=tbPrime[i], (p*p) <= p1 ) ;i++) {
            if((p1 % p) == 0) break ;
        }
        if(i==nbPrime ||  p*p > p1) {
            nb++ ;
            //    if(pbR->isVerbose) fprintf(stdout," %d",p1)  ;
        }
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Found %d\n",pbR->ident,nb) ;
    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nb) ;
    return 1 ;
}


#define PB132_POW10     1000000000
#define PB132_MAXP      40
#define PB132_PRIME     1000000
int PB132(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB132_PRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int32_t pow = PB132_POW10 ;
    int32_t nbFind = 0 ;
    int32_t sum = 0 ;
    int i ;
    // on saute 2,3,5
    for(i=3;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        if( pow10_modn(pow,p)==1) {
            nbFind++ ;
            sum += p ;
//            printf("%d ",p);
            if(nbFind >= PB132_MAXP) break ;
        }
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Sum[%d prime factors]=%d\n",pbR->ident,nbFind,sum) ;
    
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",sum) ;
    return 1 ;
}

#define PB133_PRIME 100000

int PB133(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB133_PRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbFind =0 ;     int i;
    int sum  = 2+3+5; // qui ne sont pas possible 2 et 5 trivial, 3 car ne nombre de 1 est 10**n non multiple de 3
    // il faut sauter 2,3,5
    for(i=3;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        int n0 = p-1 ;
        int q ;
        int d1 ;
        int dmin = n0 ;
        for(d1=2;d1*d1<=n0;d1++) {
            if((n0 % d1) == 0) {
                if(pow10_modn(d1,p)==1) {
                    dmin = d1; break ; // c'est forcement la plus petite valeur
                }
                int d2 = n0 / d1 ;
                if(pow10_modn(d2,p)==1) {
                    dmin = d2 ;
                }
            }
        }
        // on teste maintenant si dmin admet comme facteur seulement 2 et 5
        while((dmin & 1) == 0) dmin /= 2;
        while((dmin % 5)==0) dmin /=  5;
        if(dmin>1) {
            sum += p ;
            nbFind++ ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",sum) ;
    return 1 ;
}

int PB133a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB133_PRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int nbFind =0 ;     int i;
    int sum  = 2+3+5; // 2 and 5 trivial, 3 because 10**n is not a multiple of 3
    // begin by 3 to skip 2,3,5
    // loop for primes
    for(i=3;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        int n0 = p-1 ;
        int g = 1 ;
        while((n0 & 1) == 0 ) { n0 /= 2; g *= 2 ; }
        while((n0 % 5) == 0 ) { n0 /= 5; g *= 5 ; }
        if(pow10_modn(g,p) !=1) {
            sum += p ;
            nbFind++ ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",sum) ;
    return 1 ;
}

#define PB134_MAXP  1000000
// n * 10**k1 + p1 = 0 mod(p2)
// <=> n = (p2-p1) * 1/ 10**k1 mod(p2)
// comme 10**p2-1 = 1 mod(p2) 10**(p2-1-k1) = 1/ 10**k1

int PB134(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB134_MAXP+10000)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int  i;
    int64_t sum  = 0;
    // loop for primes, begin by 5
    int pow10 = 1 , k1 = 0 ;
    int p1 ;
    for(i=2;(p1=tbPrime[i])<PB134_MAXP;i++) {
        int p1 = tbPrime[i] ;
        int p2 = tbPrime[i+1] ;
        while(pow10<p1) { k1++ ; pow10 *= 10 ; }
        int n =  ((p2-p1) * pow10_modn(p2-1-k1,p2)) % p2 ;
        sum += n*(int64_t)pow10+p1 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",sum) ;
    return 1 ;
}

int inverse(int a, int n) {
    int t = 0;     int nt = 1;
    int r = n;     int nr = a;
    while(nr != 0) {
        int q = r / nr ;
        int tmp = nt ;
        nt = t -q * nt ;
        t = tmp ;
        tmp = nr ;
        nr = r - q * nr ;
        r = tmp ;
    }
    if ( r > 1 )  return 0 ;
    else if (t < 0 ) t +=  n ;
    return t ;
}
// on inverse 10**k1 mod(p2) par extended gcg
int PB134a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB134_MAXP+10000)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int  i;
    int64_t sum  = 0;
    // loop for primes, begin by 5
    int pow10 = 1 , k1 = 0 ;
    int p1 ;
    for(i=2;(p1=tbPrime[i])<PB134_MAXP;i++) {
        int p1 = tbPrime[i] ;
        int p2 = tbPrime[i+1] ;
        while(pow10<p1) { k1++ ; pow10 *= 10 ; }
        int n = ((p2-p1) * inverse(pow10 % p2,p2)) % p2 ;
        sum += n*(int64_t)pow10+p1 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",sum) ;
    return 1 ;
}


#define PB135_MAX 1000000
// x = y+d , z=y-d
// x**2 - y**2 -z**2 = (4*d-y)* y
// so n=d1*d2
// d2=y ; 4*d-y=d1  <=> d = (d2+d1)/4
// z>0 <=> d<y <=> d2 < 3 * d1
int PB135(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n ;
    int num  =0 ;
    for(n=2;n<PB135_MAX;n++) {
        if((n & 3) ==1 || (n & 3) == 2) continue ;
        int d1,d2 ;
        int nb= 0 ;
        for(d1=1;d1*d1<=n;d1++) {
            d2 = n / d1 ;
            if((d1+d2) & 3 ) continue ;
            if(n==d1*d2) {
                nb++ ;
                if(d2< 3*d1 && d1 != d2) nb++ ;
                if(nb > 10) break ;
            }
        }
        if(nb==10) num++;
     }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",num) ;
    return 1 ;
}
// replace loop on n, by loop on d1 ,d2
// and  as n % 3 = 0 or 3 so we can hash for n/2
//
int PB135a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    uint8_t *nbSol = calloc(PB135_MAX/2,sizeof(nbSol[0])) ;
    int n ;
    int d1, d2 ;
    for(d1=1;d1<PB135_MAX;d1++) {
        int d2max = 3*d1-1 ;
        if(d2max > PB135_MAX/d1 ) d2max = PB135_MAX/d1 ;
        for(d2= 4 - (d1 & 3) , n = (d1*d2)/2 ;d2<= d2max ;d2+=4 , n+= 2*d1) {
            if(nbSol[n]<=10) nbSol[n]++ ;
        }
    }
    int num =0 ;
    for(n=2;n<PB135_MAX/2;n++) {
        if(nbSol[n]==10)  num++ ;
    }
    free(nbSol);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",num) ;
    return 1 ;
}

#define PB136_MAX   50000000
// if n % 4 == 0 n = 4 * p with p prime odd
// if n multiple of 8 no decomposition (d1+d2)%4 == 0
// if n multiple of 16 n= 16 * p , p prime odd
// if n multiple of 32 two decompositions : 4 x n/4 ; 8 x n/8
// if n % 4 == 3
//     n odd, and any decomposition d1*d2 gives a solution so n must be prime
//
int PB136(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB136_MAX)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    int nbPrime = GetNbPrime(ctxP) ;
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbSol = 2 ; // for n= 4 and n=16
    int i ;
    for(i=1;i<nbPrime;i++) {
        int p = tbPrime[i] ;
        if((p & 3) == 3) nbSol++ ;
        if(p <= PB136_MAX/4) {
            nbSol++ ;
            if(p <= PB136_MAX/16) nbSol++ ;
        }
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%d",nbSol) ;
    return 1 ;
}

#define PB137_RG    15
int PB137(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int rg ;
    int64_t x1 = 3 , y1 = 1 ;
    int64_t xn=x1 , yn= y1 ;
    int64_t N=((xn+yn)*y1)/2 ;
    printf("%lld=((%lld+%lld)x%lld)/2\t\t %lld**2 - 5 x %lld**2 = %d\n",N,xn,yn,yn,xn,yn,(int)(xn*xn-5*yn*yn));
   for(rg=2;rg<=PB137_RG;rg++) {
        int64_t xn1 = (x1*xn+5*y1*yn)/2 ;
        int64_t yn1= (x1*yn+xn*y1) /2 ;
        N = ((xn1+yn1)*yn1)/2 ;
        printf("%lld=((%lld+%lld)x%lld)/2\t\t %lld**2 - 5 x %lld**2 = %d\n",N,xn1,yn1,yn1,xn1,yn1,(int)(xn1*xn1-5*yn1*yn1));
        xn=xn1 ;
        yn = yn1 ;
        
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",N) ;
    return 1 ;
}



// donc on cherche a decomposer un triangle Pythagoricien primaire (m,n) m>n premiers m**m - n**2 , 2mn


#define PB139_MAX 100000000
//  brute force
int PB139(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i,n,m,a,b,c,p ;
    int32_t nbT = 0 ;
    for(m=2;;m++) {
        int nb = 0 ;
        for(n=(m&1) + 1;n<m;n += 2) {
            if(PGCD(m,n) != 1)  continue ;
             a = m*m - n*n ; b = 2*m*n ; c = m*m + n*n ;
            p = a+b+c ;
            if(p> PB139_MAX) break ;
            int dab = (a > b) ? (a-b) : (b-a) ;
            if((c % dab)==0) {
                nbT += PB139_MAX/p ;
                if(pbR->isVerbose)printf("\t PB%s [(%d,%d)%d,%d,%d,%d]\n",pbR->ident, m,n,a+b+c,a,b,c) ;
            }
            nb++ ;
            
        }
        if(nb==0) break ;
  //      printf("\n");
    }
//    if(pbR->isVerbose)printf("\t PB%s Nb[%d]=%d => Nb[%d]=%d\n",pbR->ident,minM-1,histoM[minM-1],minM,histoM[minM]) ;
    sprintf(pbR->strRes,"%d",nbT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
//
// on peut demonter (a-b) (ou (b-a)à divide c entraine a-b= +-1
// les triangles pythagoriciens sont tels m**2-n**2 -2mn = +-1
// <=> (m-n)**2 - 2n**2 = =-1 (equation de Pell
// Donc on cherche x**2 - 2 y**2 = 1 et x**2 - 2 y**2 = -1
// m = x+y ; n = y ;

int PB139a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int xp1 = 3 , yp1 = 2 ;
    int xm1 = 1 , ym1 = 1 ;
    int32_t xp = 1 , yp = 0 ;
    int32_t xm = 1 , ym = 0 ;
    int32_t nbT = 0 ;
    int nb  ;
    do {
        int32_t m,n,a,b,c;
        uint32_t p ;
        nb = 0 ;
        xm = xp*xm1 + 2*yp*ym1 ;  ym = xp*ym1 + yp*xm1 ;
        m = xm+ym ; n = ym ;
        a = m*m -n*n ; b = 2*m*n ; c= m*m +n*n ;
        p= a+b+c ;
        if(p<= PB139_MAX) {
            nb += PB139_MAX/p ;
            if(pbR->isVerbose)printf("\t PB%s [(%d,%d(%d)->%d,%d)%d,%d,%d,%d]\n",pbR->ident,xm,ym,xm*xm-2*ym*ym,m,n,p,a,b,c) ;
        }
        int32_t tmp = xp ;
        xp = xp*xp1 + 2*yp*yp1 ;    yp = tmp*yp1 + yp*xp1 ;
        m = xp+yp ; n = yp ;
        a = m*m -n*n ; b = 2*m*n ; c= m*m +n*n ;
        p= a+b+c ;
        if(p<= PB139_MAX) {
            nb += PB139_MAX/p ;
            if(pbR->isVerbose)printf("\t PB%s [(%d,%d(%d)->%d,%d)%d,%d,%d,%d]\n",pbR->ident,xp,yp,xp*xp-2*yp*yp,m,n,p,a,b,c) ;
        }
        nbT += nb ;
    } while(nb > 0) ;
    sprintf(pbR->strRes,"%d",nbT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



#define PB140_RG    30
int PB140(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int rg ;
    int64_t Sum = 0 ;
    int64_t x44e = 7 , y44e= 1 ;
    int64_t x44o = 8 , y44o= 2 ;
    int64_t x1 = 3 , y1 = 1 ;
    int64_t xn=2 , yn= 0 ;
    int64_t N ,xn44,yn44;
    for(rg=1;rg<=PB140_RG;rg++) {
        int64_t xn1 = (x1*xn+5*y1*yn)/2 ;
        int64_t yn1= (x1*yn+xn*y1) /2 ;
        xn = xn1 ;
        yn = yn1 ;
        
        if(rg & 1) {
            xn44 = (x44o*xn+5*y44o*yn)/2 ;
            yn44 = (x44o*yn+y44o*xn)/2 ;
            N=(-7+xn44)/5 ;
            
        }else {
            xn44 = (x44e*xn+5*y44e*yn)/2 ;
            yn44 = (x44e*yn+y44e*xn)/2 ;
            N=(-7+xn44)/5 ;
        }
        Sum += N ;
        printf("%lld=(%lld-7)/5\t\t %lld**2 - 5 x %lld**2 = %d\n",N,xn44,xn44,yn44,(int)(xn44*xn44-5*yn44*yn44));
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",Sum) ;
    return 1 ;
}



#define PB143_MAX   120000
// brute force
// find all (p,q) tel que p**2+p*q+q**2 est un carre parfait
// puis trouver (p,q), (p,r) tel que (q,r) OK

int PB143(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int32_t p,r ;
    int32_t M[1000] ;
    uint64_t Sum = 0 ;
    uint8_t *isFound = calloc(PB143_MAX,sizeof(isFound[0])) ;
    int totNb = 0 ;
    for(p=1;p<PB143_MAX;p++) {
        int nb =0 ;
        int M0 = 0 ;
        for(r=p;r+p+M0<PB143_MAX;r++) {
//            if(PGCD(p,r) != 1) continue ;
            int64_t S =p*(int64_t)p+r*(int64_t)r+p*(int64_t)r ;
            int64_t sq =Sqrt64(S);
            if(S==sq*sq) {
//                printf("(%d,%d)->%d ",p,r,sq);
                if(nb==0) M0 = r ;
                M[nb++] = r ;
            }
            
        }
        if(nb>1) {
            int i , j ;
            totNb += nb ;
            for(i=1;i<nb;i++) {
                for(j=0;j<i && p+M[i]+M[j]< PB143_MAX;j++) {
                    int64_t S = M[i]*(int64_t)M[i]+M[j]*(int64_t)M[j]+M[i]*(int64_t)M[j] ;
                    int64_t sq =Sqrt64(S);
                    if(S==sq*sq) {
                        int32_t pqr = p+M[i]+M[j] ;
                        if(isFound[pqr] == 0) {
                            Sum += pqr ;
                            isFound[pqr] = 1 ;
                        }
//                        printf("%lld+%lld+%lld=%lld ",p,M[i],M[j],pqr) ;
                    }
                }
            }
        }
    }
    printf("TotNb=%d\n",totNb);
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",Sum) ;
    return 1 ;
}

typedef struct DIOPH143 {
    int32_t p ;
    int32_t q ;
    
} DIOPH143 ;

int CmpDioph43(const void *e1, const void *e2) {
    const DIOPH143 * d1 = (DIOPH143 *) e1 ;
    const DIOPH143 * d2 = (DIOPH143 *) e2 ;
    int diff = d1->p - d2->p  ;
    if(diff) return diff ;
    else return d1->q - d2->q ;
}

#define PB143_MAXDIOPH  150000
//
// on genere les p,q par p =m*m-n*n ; q=n*(2*m+n)
// on trie puis idem
int PB143a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    uint64_t Sum = 0 ;
    uint8_t *isFound = calloc(PB143_MAX,sizeof(isFound[0])) ;
    DIOPH143 *dioph = malloc(PB143_MAXDIOPH*sizeof(dioph[0]));
    int nbDioph = 0 ;
    int m,n ;
    int mMax = Sqrt32(PB143_MAX) ;
    for(m=2;m<mMax;m++) {
        for(n=1;n<m;n++) {
            int p = m*m -n*n ;
            int q = n*(2*m+n) ;
            int s = p+q ;
            if(s >= PB143_MAX) continue ;
            if(PGCD(p,q) != 1) continue ;
            if(p < q) {
                int k ;
                for(k=1;k*s<PB143_MAX;k++) {
                        dioph[nbDioph].p = p*k ;
                        dioph[nbDioph++].q = q*k ;
                }
            } else {
                int k ;
                for(k=1;k*s<PB143_MAX;k++) {
                        dioph[nbDioph].p = q*k ;
                        dioph[nbDioph++].q = p*k ;
                }
                
            }
        }
    }
    qsort(dioph,nbDioph,sizeof(dioph[0]),CmpDioph43) ;
    memset(isFound,0,PB143_MAX);
    dioph[nbDioph].p = 0 ; // terminator
    int i;
    for(i=0;i<nbDioph-1;i++) {
        int j;
        for(j=i;dioph[i].p == dioph[j+1].p ; j++) ;
        if(j>i) {
            int p = dioph[i].p ;
            int k1,k2 ;
            for(k1=i+1;k1<=j;k1++) {
                int32_t pqr ;
                for(k2=i; k2<k1 && (pqr=p+dioph[k1].q+dioph[k2].q)< PB143_MAX; k2++) {
                    int64_t S = dioph[k1].q*(int64_t)dioph[k1].q+dioph[k2].q*(int64_t)dioph[k2].q+dioph[k1].q*(int64_t)dioph[k2].q ;
                    int64_t sq =Sqrt64(S);
                    if(S==sq*sq) {
                        if(isFound[pqr] == 0) {
                            Sum += pqr ;
                            isFound[pqr] = 1 ;
                        }
                    }
                }
            }            
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",Sum) ;
    return 1 ;
}

#define PB145_MAXD  9
int PB145a(PB_RESULT *pbR) {
    uint64_t nbSol = 0 ;
    pbR->nbClock = clock() ;
    uint64_t nbNivSol[4] ;
    nbNivSol[2]=20 ;
    nbNivSol[3]=5*20 ;
    nbNivSol[0] = 20*30 ;
    nbNivSol[1] = 0 ;

    int lg ;
    for(lg=2;lg<=PB145_MAXD;lg++) {
        if(lg >= 6) {
            if(lg & 1) nbNivSol[lg%4] *= 500 ;
            else nbNivSol[lg%4] *= 900 ; ;
        }
        nbSol += nbNivSol[lg % 4] ;
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s for %d digits NB=%lld\n",pbR->ident,lg,nbSol) ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbSol) ;
    return 1 ;
}

int PB145(PB_RESULT *pbR) {
    int64_t nbSol = 0 ;
    int digits[(PB145_MAXD+1)/2] ;
    int nbDecomp[19] =  { 1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1} ;
    int nbDecomp0[19] = { 0,0,1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1} ;
    pbR->nbClock = clock() ;
    int lg ;
    int i,c ;
    for(lg=2;lg<=PB145_MAXD;lg++) {
        if((lg % 4)==1) continue ;
        digits[0]=2 ;
        int hlg = lg/ 2 ;
        for(i=1;i<hlg;i++) { digits[i] = 0; }
        if(lg & 1) digits[i++] = 0 ;
        int is =lg;
        while(is>= 0) {
            for(i=0,c=0;i<hlg;i++) {  c += digits[i] ; if((c&1) == 0) break ;  c = (c>=10) ? 1 : 0 ;  }
            if(i== hlg) {
                for(;i<lg;i++) { c += digits[lg-i-1] ; if((c&1) == 0) break ; c = (c>=10) ? 1 : 0 ; }
                if(i==lg) {
                    int64_t isup = nbDecomp0[digits[0]] ;
                    for(i=1;i<hlg;i++)  isup *= nbDecomp[digits[i]] ;
                    nbSol += isup ;
                }
            }
            if(lg & 1) {
                if(digits[hlg] < 18) { digits[hlg] += 2 ; continue ;
                } else {  digits[hlg] = 0 ;   }
            }
            for(is=hlg-1;is>=0;is--) {
                if(digits[is]<18) {  digits[is]++ ;  break ;
                } else {
                    if(is==0) { is=-1 ; break ; }
                    digits[is] = 0 ;
                }
            }
        }
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s for %d digits NB=%lld\n",pbR->ident,lg,nbSol) ;
     }
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",nbSol) ;
    return 1 ;
}




#define PB147_MD    10
void SumInPlace147(int nbRec[PB147_MD*PB147_MD]) {
    int i , j ;
    for(i=0;i<PB147_MD;i++) {
        for(j=0;j<=i;j++) {
            if(i > 0 && j > 0) {
                int S = nbRec[PB147_MD*i+j] + nbRec[PB147_MD*(i-1)+j] + nbRec[PB147_MD*i+j-1]-nbRec[PB147_MD*(i-1)+j-1] ;
                nbRec[PB147_MD*i+j] = nbRec[PB147_MD*j+i]= S ;
            } else if(i > 0) {
                int S = nbRec[PB147_MD*i+j] + nbRec[PB147_MD*(i-1)+j] ;
                nbRec[PB147_MD*i+j] = nbRec[PB147_MD*j+i]= S ;
            } else if(j > 0) {
                int S = nbRec[PB147_MD*i+j] + nbRec[PB147_MD*i+j-1] ;
                nbRec[PB147_MD*i+j] = nbRec[PB147_MD*j+i]= S ;
            }
//            printf("%8d%c",nbRec[PB147_MD*i+j],(nbRec[PB147_MD*i+j] == ((i+1)*(i+2)*(i+3))/6 * ((j+1)*(j+2)*(j+3))/6) ? '=' : '!' ) ;
        }
//        printf("\n");
    }
}
int64_t NBRec(int64_t i, int64_t j) {
    return (i*(i+1))/2 * (j*(j+1))/2 ;
}
int64_t NBSumRec(int64_t i,int64_t j) {
    return (i*(i+1)*(i+2))/6 * (j*(j+1)*(j+2))/6 ;
}
int Check147(int nbRec[PB147_MD*PB147_MD], int64_t(*Fcheck)(int64_t,int64_t),FILE * fout) {
    int i,j,nbErr = 0 ;
    for(i=0;i<PB147_MD;i++) {
        for(j=0;j<PB147_MD;j++) {
            int err = (nbRec[PB147_MD*i+j] != Fcheck(i+1,j+1));
            if(err)nbErr++ ;
            if(fout)fprintf(fout,"%8d%c",nbRec[PB147_MD*i+j],(err) ? '*' : ' ' ) ;
        }
        if(fout)printf("\n");
    }
    return nbErr ;
}
int64_t NBRec2(int64_t i, int64_t j) {
    if(j<i) { // permutation
        j -= i ;
        i += j ;
        j = i - j ;
    }
    // i <= j
    return ((i-1)*i*(4*i*i+4*i+3))/6 + (j-i) * (i * (4*i*i - 1))/3 ;
}
int64_t NBSumRec2(int64_t i, int64_t j) {
    if(j<i) { // permutation
        j -= i ;
        i += j ;
        j = i - j ;
    }
    // i <= j
    return ((i-1)*i*(2*i*i*i*i+8*i*i*i+13*i*i+8*i+1))/30
        + (j-i) * ((i-1)*i*(3*i*i*i+8*i*i+8*i+3))/15
        + ((j-i)*(j-i+1))/2 * (i*(2*i*i*i+4*i*i+i-1) )/6;
}


int PB147(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;

    int  i,j,k;
    int64_t sum  = 0;
    int nbRec[PB147_MD*PB147_MD],nbRec2[PB147_MD*PB147_MD] ;
    int xMin[2*PB147_MD];
    int xMax[2*PB147_MD];
    nbRec[0] = 1 ;
    int nbErr = 0 ;
    for(i=0;i<PB147_MD;i++) {
        int ni  = ((i+1)*(i+2))/2 ;
        for(j=0;j<=i;j++) {
            nbRec[PB147_MD*i+j] = nbRec[PB147_MD*j+i]= ni * nbRec[j*PB147_MD] ;
        }
    }
    nbErr+=Check147(nbRec,NBRec,NULL);
    SumInPlace147(nbRec) ;
    nbErr+=Check147(nbRec,NBSumRec,NULL);

    for(i=1;i<=PB147_MD;i++) {
        for(j=i;j<=PB147_MD;j++) {
            for(k=0;k<=i;k++) {  xMin[k] = -k ; xMax[k] = k ;  }
            for(   ;k<=j;k++) { xMin[k] = k- 2 * i  ; xMax[k] = k ;  }
            for( ;k<=i+j;k++) { xMin[k] = k- 2 * i  ; xMax[k] = 2*j - k  ;   }
            //
            int ay,ax ;
            int S = 0 ;
            for(ay=i+j;ay>0;ay--) {
                for(ax=xMin[ay];ax<xMax[ay];ax++) {
                    int by ;
                    for(by=ay-1;by>0;by--) {
                        int bxMin = xMin[by] ;
                        if(ax < bxMin) continue ;
                        if(bxMin < ax+1) bxMin = ax+1 ;
                        int bxMax = xMax[by] ;
                        if(bxMax > xMax[ay]) bxMax = xMax[ay] ;
                        if(bxMax>=bxMin) S += bxMax - bxMin +1 ;
                    }
                }
                
            }
            nbRec2[PB147_MD*(i-1)+(j-1)] = nbRec2[PB147_MD*(j-1)+(i-1)] = S ;
        }
    }
    nbErr+=Check147(nbRec2,NBRec2,NULL);
    SumInPlace147(nbRec2) ;
    nbErr+=Check147(nbRec2,NBSumRec2,NULL);
    
    int64_t S = NBSumRec(43,47)+NBSumRec2(43,47) ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s NBerr=%d Found=%lld\n",pbR->ident,nbErr,S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    sprintf(pbR->strRes,"%lld",S) ;
    return 1 ;
}
