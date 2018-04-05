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
                // aprroximation by
                double logFk = k * log10((1+sqrt(5))/2) + log10(1/sqrt(5)) ;
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



