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

typedef struct Edge {
    u_int16_t   cost;
    u_int16_t   Nbeg ;
    u_int16_t   Nend ;
} Edge ;

typedef struct nodeTree { // tree for a node
    u_int16_t numTree ;
    u_int16_t nxtNode ; // link to nxt node in the tree
    
} nodeTree ;

typedef struct Tree {
    u_int16_t length ;
    u_int16_t firstNode ;
} Tree ;

int CmpEdge( const void *edg1, const void *edg2) {
    return ((Edge *)edg1)[0].cost - ((Edge *)edg2)[0].cost ;
}
int PB107(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    Edge EDG[PB107_SIZE*PB107_SIZE] ;
    Tree TR[PB107_SIZE+1] ;
    nodeTree nT[PB107_SIZE] ; // association node -> tree
    u_int16_t nbTree = 0 ;
    int i, j ;
    const u_int16_t * cost = P107_GetData();
    u_int16_t nbEdge = 0 ;
    for(i=0;i<PB107_SIZE;i++) { // get the non null edge
        for(j=i+1;j<PB107_SIZE;j++) {
            if(cost[i*PB107_SIZE+j]){
                EDG[nbEdge].Nbeg = i ;
                EDG[nbEdge].Nend = j ;
                EDG[nbEdge].cost = cost[i*PB107_SIZE+j] ;
                nbEdge++ ;
            }
        }
    }
    memset(nT,0,sizeof(nT)) ;
    heapsort(EDG,nbEdge,sizeof(EDG[0]),CmpEdge) ; // sort edges by cost
    int indSortEdg = 0 ;
    TR[0].length = 0 ;
    int maxLength = 0 ;
    u_int32_t savingCost = 0 ;
    do {    // loop to build tree , adding the min cost edge
        Edge curEdg = EDG[indSortEdg++] ;
        u_int16_t begT= nT[curEdg.Nbeg].numTree ;
        u_int16_t endT =nT[curEdg.Nend].numTree ;
        if(begT && endT && begT == endT) {
            savingCost += curEdg.cost ; // the edge address the same tree => loop, forbiddeen
            continue ;
        }
        if(begT == 0 && endT == 0 ) { // nodes not in e tree => simple tree with 2 nodees
               nbTree++ ;
                TR[nbTree].length = 2 ;
                if(maxLength < 2) maxLength = 2 ;
                TR[nbTree].firstNode = curEdg.Nbeg ;
                nT[curEdg.Nbeg].numTree = nT[curEdg.Nend].numTree = nbTree ;
                nT[curEdg.Nbeg].nxtNode = curEdg.Nend ;
                nT[curEdg.Nend].nxtNode = curEdg.Nbeg ;
                continue;
        } else if( begT==0 || endT == 0) { // One node not in a tree, rattach to the tree
            if(endT == 0) { // permutation  end <-> beg to shorten code
                endT = begT ;
                u_int16_t tmp = curEdg.Nbeg ;
                curEdg.Nbeg = curEdg.Nend ;
                curEdg.Nend = tmp ;
            }
            nT[curEdg.Nbeg].numTree = endT ;
            nT[curEdg.Nbeg].nxtNode = nT[curEdg.Nend].nxtNode  ;
            nT[curEdg.Nend].nxtNode = curEdg.Nbeg ;
            TR[endT].length++ ;
            if(maxLength < TR[endT].length) maxLength = TR[endT].length ;
            continue ;
        } else { // two trees, rattach the shortest to the longuest
            if(TR[begT].length > TR[endT].length ) { // permutation end <->beg
                u_int16_t tmp = endT ;
                endT = begT ;
                begT= tmp ;
                tmp = curEdg.Nbeg ;
                curEdg.Nbeg = curEdg.Nend ;
                curEdg.Nend = tmp ;
            }
          // rattach beg tree to end tree
            u_int16_t nxtEnd = nT[curEdg.Nend].nxtNode ;
            u_int16_t nxtBeg = curEdg.Nbeg ;
            nT[curEdg.Nend].nxtNode = nxtBeg ;
            u_int16_t antBeg ;
            do {
                antBeg= nxtBeg ;
                nT[nxtBeg].numTree = endT ;
                nxtBeg = nT[nxtBeg].nxtNode ;
            } while(nxtBeg != curEdg.Nbeg) ;
            nT[curEdg.Nbeg].numTree = endT ;
            nT[antBeg].nxtNode = nxtEnd ;
            TR[endT].length += TR[begT].length ;
            if(maxLength < TR[endT].length) maxLength = TR[endT].length ;
            continue ;
        }
    } while(maxLength != PB107_SIZE) ;
    while(indSortEdg<nbEdge) {
        savingCost += EDG[indSortEdg++].cost ;
    }
    sprintf(pbR->strRes,"%d",savingCost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



