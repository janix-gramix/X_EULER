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

#define PB103_MAXNB     16
#define PB103_NB  9
#define PB103_MAX_DELTA   300

typedef int32_t sum103_t ;

int CheckEquality(sum103_t *v,int lg) {
    int k ;
    for(k=2;2*k<=lg;k++) {
        int j ;
        u_int8_t perm2[PB103_MAXNB] ;
        for(j=0;j<2*k;j++)perm2[j] = j ;
        do {
            sum103_t S = 0 ;
            for(j=0;j<2*k;j++){
                S += v[perm2[j]] ;
            }
            if (S & 1) continue ; // pas divisible par 2
            S  /= 2 ;
            u_int8_t perm[PB103_MAXNB] ;
            for(j=0;j<k;j++)perm[j] = j ;
            do {
                int D = S ;
                for(j=0;j<k;j++){
                    D -= v[perm2[perm[j]]] ;
                    if( D == 0) {
                        return 0 ;
                    }
                }
            } while(NextSub(perm,k,2*k) >= 0 ) ;
        } while(NextSub(perm2,2*k,lg) >= 0) ;
    }
    return 1;
}

typedef struct AlterPaths {
    int maxK ;  // k max for the path length (k in 2k)
    int * nbK ; // number of path for length k  (nb=1 por k=2, nk=5, for k= 3
    //   nb(k) = 1/2 k(-1)/(k+1) x C[2k,k]
    int * begK ;// index for the first path of lengh k
    // a path is k indexes  in [0, 2k[
    int16_t *path ;
} AlterPaths ;

typedef struct CheckPaths {
    int maxN ; // N max, = length of ensemble to check
    int * npermK ; // number of permutation og lengh k
    int * nsumK ; // number of sum (or permutations to check) of length
                  // one more than hte number of sum (for global sum)
    int *begK ; // index for the first sum to calculate (of length k)
    int16_t *indSum ; // list on index for sum
    
} CheckPaths ;

typedef struct HalfPaths {
    int maxN ; // N max, = length of ensemble to check
    int npermk ; // number of permutation of lengh k
    int nsumk ; // number of sum (or permutations to check) of length
    int16_t * indSum ; // list on index for sum
    
} HalfPaths ;


AlterPaths * FreeAlterPath(AlterPaths * AltP) {
    if(AltP != NULL) {
        free(AltP->begK) ;
        free(AltP->nbK);
        free(AltP->path);
        free(AltP);
    }
    return NULL ;
}
AlterPaths * GetAlterPath(int maxk) {
    AlterPaths * AltP = calloc(1,sizeof(AltP[0])) ;
    AltP->nbK = malloc(maxk*sizeof(AltP->nbK[0])) ;
    AltP->begK = malloc(maxk*sizeof(AltP->begK[0])) ;
    int k ;
    AltP->nbK[0] = 0 ;
    AltP->begK[0] = 0 ;
    AltP->maxK = maxk ;
    for(k=1;k<maxk;k++) {
        AltP->begK[k] = AltP->begK[k-1]+k*AltP->nbK[k-1] ; //AltP->nbK[k-1] de longueur k
        int j;
        int k1= k+1 ; // pour tenir compte du decalage d'indice
        u_int64_t N = k1-1;
        for(j=k1+2;j<=2*k1;j++) { N *= j ;} ; // on par de k1+2 car 1/(k1+1)
        for(j=2;j<=k1;j++) N /= j ;
        AltP->nbK[k] = (int) N/2 ;
    }
    AltP->path = malloc((AltP->begK[maxk-1]+ maxk * AltP->nbK[maxk-1]) * sizeof(AltP->path[0])) ;
    u_int8_t ind[PB103_MAXNB] ;
    for(k=1;k<maxk;k++) { // on va generer les path
        int k1= k+1 ; // pour tenir compte du decalage d'indice
        int is = AltP->begK[k] ;
        int j ;
        for(j=0;j<k1-1;j++) ind[j]=j ; // choix de k1-1 parmis 2k1
        //       printf("*****Gen %d, beg=%d \n",k1,is) ;
        do {
            // on recopie ind[0..2k1] en appliquant la reverse path(2n,-2) => path(2n,0)
            int ij,  s = 1 ; // car premier element force a 1
            AltP->path[is++] = ind[0];
            for(j=1,ij=1;j<2*k1;j++) {
                if(ij < k1-1 && j==ind[ij]) { // +1
                    s++ ;
                    AltP->path[is++] = j ;
                    ij++ ;
                } else { // -1
                    s-- ;
                    if(s==-1) { j++ ; break ;} // on arrive en mode reverse
                }
            }
            for(;j<2*k1;j++) {
                if(ij < k1-1 && j==ind[ij]) { // +1
                    ij++ ;
                    s-- ;  // normalement inutile de compter
                } else { // -1
                    s++ ;
                    AltP->path[is++] = j ;
                }
            }
            //           for(j=0;j<k1;j++) printf("%d%c",AltP->path[is-k1+j],j==k1-1 ? '\n' : ',') ;
        } while(NextSub(ind+1, k1-2, 2*k1)>=0) ; // on contraint le premier a etre +1 pour gerer la symetrie
        //       printf("*****End %d, beg=%d \n",k1,is) ;
    }
    return AltP ;
}

CheckPaths * FreeCheckPath(CheckPaths * chkP) {
    if(chkP != NULL) {
        free(chkP->begK) ;
        free(chkP->nsumK);
        free(chkP->indSum);
        free(chkP->npermK);
        free(chkP);
    }
    return NULL ;
}

CheckPaths * GetCheckPath(int N, AlterPaths * altP,int isSup) {
    CheckPaths * chkP = calloc(1,sizeof(chkP[0])) ;
    int maxk = N/2 ;
    int dN1 = isSup ? 1 : 0 ;
    chkP->nsumK = malloc(maxk*sizeof(chkP->nsumK[0])) ;
    chkP->npermK = malloc(maxk*sizeof(chkP->npermK[0])) ;
    chkP->begK = malloc(maxk*sizeof(chkP->begK[0])) ;
    int k ;
    chkP->nsumK[0] = 0 ;
    chkP->begK[0] = 0 ;
    chkP->npermK[0] = 0 ;
    chkP->maxN = N ;
    for(k=1;k<maxk;k++) {
        int j ;
        int k1= k+1 ; // pour tenir compte du decalage d'indice
        chkP->begK[k] = chkP->begK[k-1]+k*chkP->nsumK[k-1]*chkP->npermK[k-1] ; //chkP->nbK[k-1] de longueur k
        u_int64_t CN2k = 1;
        for(j=0;j<2*k1-dN1;j++) { CN2k *= N-dN1-j ; }
        for(j=2;j<=2*k1-dN1;j++) { CN2k /= j ; }
        chkP->npermK[k] = (int) CN2k ;
        chkP->nsumK[k] = (altP->nbK[k] + 1)  ;
    }
    chkP->indSum = malloc((chkP->begK[maxk-1]+ maxk * chkP->nsumK[maxk-1])* chkP->npermK[maxk-1] * sizeof(chkP->indSum[0])) ;
    u_int8_t perm2[PB103_MAXNB] ;
    for(k=1;k<maxk;k++) { // on va generer les index
        int k1= k+1 ; // pour tenir compte du decalage d'indice
        int is = chkP->begK[k] ;
        int j ;
        if(dN1) perm2[2*k1-1] = N-1 ;
        for(j=0;j<2*k1-dN1;j++)perm2[j] = j ;
 //       printf("*****Gen %d, beg=%d \n",k1,is) ;
        do {
            int16_t * ind  = altP->path + altP->begK[k] ;
            int nb = altP->nbK[k] ;
            // add of the complement of the first sum
            int ij ;
            for(ij=0 ;ij<2*k1;) {
                while(*ind != ij) {
                    chkP->indSum[is++] = perm2[ij++] ;
                }
                ind++; ij++ ;
            }
 //           for(j=0;j<k1;j++) printf("%d%c",chkP->indSum[is-k1+j],j==k1-1 ? ' ' : '.') ;
            ind -= k1 ;
            while(nb-- > 0) {
                for(j=0;j<k1;j++) {
                    chkP->indSum[is++] = perm2[*ind++] ;
                }
  //              for(j=0;j<k1;j++) printf("%d%c",chkP->indSum[is-k1+j],j==k1-1 ? '\n' : ',') ;
            }
        } while(NextSub(perm2,2*k1-dN1,N-dN1) >= 0) ; //
  //      printf("*****End %d, beg=%d \n",k1,is) ;
    }
    return chkP ;

    
}

HalfPaths * FreeHalfPath(HalfPaths * hlfP) {
    if(hlfP != NULL) {
        free(hlfP->indSum);
        free(hlfP);
    }
    return NULL ;
}

HalfPaths * GetHalfPath(int N) {
    HalfPaths * hlfP = calloc(1,sizeof(hlfP[0])) ;
    int k = N/2 ;
    hlfP->maxN = N ;
    {
        int j ;
        u_int64_t CNk = 1;
        for(j=0;j<k;j++) { CNk *= N-j ; }
        for(j=2;j<=k;j++) { CNk /= j ; }
        hlfP->npermk = (int) CNk ;
    }
    hlfP->indSum = malloc(k * hlfP->npermk * sizeof(hlfP->indSum[0])) ;
    u_int8_t perm2[PB103_MAXNB] ;
    int is = 0 ;
    int j ;
    for(j=0;j<k;j++)perm2[j] = j ;
    do {
        for(j=0;j<k;j++) hlfP->indSum[is++] = perm2[j] ;
    } while(NextSub(perm2,k,N) >= 0) ; //
    
    return hlfP ;
    
    
}


// Comparaison avec parcours des sous-ensembles a 2*k elements du total a N elements
// Les indices des sous-ensembles a k elements sont precalcules dans AltP
// Ce precalcul ne depend pas de N
// si isSup est mis on suppose que le test a N-1 element est OK
// donc on ne prend que les sous-ensemble 2*k contenant N-1
static int CheckEqualityAlt(sum103_t *v,int N,AlterPaths * AltP,int isSup) {
    int k ;
    int dN = (isSup) ? 1 : 0 ;
    for(k=2;2*k<=N;k++) {
        int j ;
        u_int8_t perm2[PB103_MAXNB] ;
        for(j=0;j<2*k-dN;j++)perm2[j] = j ;
        if(dN) perm2[2*k-1] = N-1 ;
        do {
            sum103_t S = 0 ;
            for(j=0;j<2*k;j++){
                S += v[perm2[j]] ;
            }
            if (S & 1) continue ; // pas divisible par 2
            S  /= 2 ;
            int16_t * ind  = AltP->path + AltP->begK[k-1] ;
            int nb = AltP->nbK[k-1] ;
            while(nb-- > 0) {
                int D = S ;
                for(j=0;j<k;j++) {
                    D -=v[perm2[*ind++]] ;
                }
                if( D == 0) {
                    return 0 ;
                }
            }
        } while(NextSub(perm2,2*k-dN,N-dN) >= 0) ;
    }
    return 1;
}

// Comparaison avec precalcul total des sous-ensembles a 2*k elements parmi N
// ainsi que des sous-ensemble a k elements a traite
// Ce precalcul DEPEND de N
// Le precalcul peut integrer le fait que N-1 est verifie
static int CheckEqualityPreN(sum103_t *v,int N,CheckPaths * chkP) {
    int k ;
    for(k=(N>>1);k>1;k--) {
        int j ;
        int nbPerm = chkP->npermK[k-1] ;
        int16_t * ind  = chkP->indSum + chkP->begK[k-1] ;
        int nbSum = chkP->nsumK[k-1] -2 ;
        while (nbPerm-- > 0) {
            int16_t * indk = ind+k ;
            sum103_t S0 = v[*ind++] , S = v[*indk++] ;
            for(j=1;j<k;j++) { S0 += v[*ind++] ; S += v[*indk++] ; }
            ind = indk ;
            if(S==S0) return 0 ;
            S += S0 ; // sum total
            if(S & 1) { ind += k * (nbSum) ; continue ; }
            S >>= 1;
            int nb = nbSum ;
            while(nb-- > 0) {
                int16_t D = v[*ind++] ;
                for(j=1;j<k;j++) {
                    D +=v[*ind++] ;
                }
                if( D == S) {
                    return 0 ;
                }
                
            }
        }
    }
    return 1;
}


static u_int8_t isSum[100000] ;

static int CheckEqualityHalf(sum103_t *v,int N) {
    int k = N/2 ;
    if(k==1) return 1;
    int j ;
    sum103_t ST = v[N-1] ;
    for(j=2;j<=k;j++) ST += v[N-j] ;
    memset(isSum,0,ST+1) ;
    
    u_int8_t perm2[PB103_MAXNB] ;
    for(j=0;j<k;j++)perm2[j] = j ;
        int il=0  ;
        sum103_t S=0 ;
        do {
            if(il < k-1) {
                S=v[perm2[0]] ;
                for(j=1;j<k-1;j++) S += v[perm2[j]] ;
            }
            if(isSum[S+v[perm2[k-1]]]) {
                return 0 ;
            }
            isSum[S+v[perm2[k-1]]] = 1 ;
        } while( (il=NextSub(perm2,k,N)) >= 0) ;
        return 1;
}

static int CheckEqualityPreH(sum103_t *v,HalfPaths * hlfP) {
    int N = hlfP->maxN ;
    int k = hlfP->maxN / 2  ;
    if(k==1) return 1;
    int j ;
    sum103_t ST = v[N-1] ;
    for(j=2;j<=k;j++) ST += v[N-j] ;
    memset(isSum,0,ST+1) ;
    int nbPerm = hlfP->npermk ;
    int16_t *ind = hlfP->indSum ;
    while(nbPerm-- > 0) {
        sum103_t S = v[*ind++] ;
        for(j=1;j<k;j++) S += v[*ind++] ;
        if(isSum[S]) {
            return 0 ;
        }
        isSum[S] = 1 ;
    }
    return 1;
}


// calcule v0 minimum (mandatory) pour que la condition ii soit verifiee
sum103_t MinCheck(sum103_t * v,int lg) {
    sum103_t minV0 = 1 ;
    int j ;
    int lg2 = (lg-1)/2 ;
    for(j=0;j<lg2;j++){
        minV0 += v[lg-1-j] - v[j+1] ;
    }
    return minV0 ;
}

// calcul le v[N-1] a partir duquel plus besoin de faire de test
// en incremental si N-1 a ete teste.
sum103_t MaxCheck(sum103_t * v,int lg) {
    sum103_t maxV0 = 1+ v[lg-2] ;
    int j ;
    int lg2 = lg/2-1 ;
    for(j=0;j<lg2;j++){
        maxV0 += v[lg-3-j] - v[j] ;
    }
    return maxV0 ;
}


int Check(sum103_t *v, int lg) {
    if(v[0] < MinCheck(v,lg)) return 0 ;
    return CheckEquality(v,lg);
}

// version sans ponderation, on parcours par deltaS=SIGMA(delta[is]) croissant
// peu performant car la contribution de v0 est preponderante et
// n'est pris en compte qu'apres coup.
// Critere d'arret avec une pauvre estimation de la valeur min de v0
// Utilise CheckEqualityAlt
int PB103(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    sum103_t values[PB103_NB] ;
    sum103_t minV0 = 32000 ;
    static int MinV0[] ={1,1,1,2,3,6,11,20,39} ;
    int16_t Smin = 32000 ;
    int AntMinV0 ;
    if(PB103_NB*sizeof(MinV0[0]) <= sizeof(MinV0)) {
        AntMinV0 = MinV0[PB103_NB-1] * PB103_NB ;
    } else {
        AntMinV0 = (((PB103_NB-1)/2) * ((PB103_NB-1)/2)) *  PB103_NB ;
    }
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    CheckPaths *chkP =GetCheckPath(PB103_NB,AltP,0) ;
        values[0] = 0 ;
        int j ;
        int is = 1 ;
        int deltaS = 0 ;
        sum103_t delta = deltaS ;
        sum103_t deltaMax = Smin - AntMinV0 ;
        sum103_t Delta[PB103_NB] ;
        for(j=0;j<PB103_NB-1;j++) Delta[j] = 0;
        Delta[PB103_NB-1] = deltaS ;
        sum103_t offsetDeltaS = (PB103_NB*(PB103_NB-1))/2  ;
        while(deltaS <= deltaMax) {
        for(;is<PB103_NB;is++) {
            values[is] = values[is-1]+Delta[is]+1 ;
        }
        sum103_t v0 = MinCheck(values,PB103_NB) ;
        sum103_t S = deltaS + PB103_NB * v0 ;
        if(S <= Smin && CheckEqualityPreN(values,PB103_NB,chkP))  {
            int j ;
            Smin = S ;
            deltaMax = Smin - AntMinV0 ;
            if(v0 < minV0) minV0= v0 ;
        
            if(pbR->isVerbose) {
                fprintf(stdout,"\t PB%s S=%d Delta=%d,minv0=%d,DeltaMax=%d ",pbR->ident,S+offsetDeltaS,deltaS,minV0,deltaMax ) ;
                for(j=0;j<PB103_NB;j++){
                    fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                }
                for(j=1;j<PB103_NB;j++){
                    fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                }
            }
            int lg = 0 ;
            for(j=0;j<PB103_NB;j++){
                lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
            }
        }
        is = PB103_NB-2 ;
        while(PB103_NB-is > delta && is){
            delta += (PB103_NB-is) * Delta[is] ;
            Delta[is-- ] = 0;
        }
        if(is == 0) {
            deltaS++ ;
            delta = deltaS  ;
            is = 1 ;
        } else {
            Delta[is]++ ; delta -= PB103_NB-is ;
        }
        Delta[PB103_NB-1] = delta ;

    }
    FreeCheckPath(chkP) ;
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

//
// calcule les ponderations des deltas[is] ds la somme
// la relation entre les v[is] et les deltas est la suivante :
// v[is] = v[is-1]+ delta[is]+1
// Les pondérations tiennent compte de 2 facteurs :
// - la recurrence ci-dessus (les delta[is] ont un poids de PB103_NB - i
// - le fait qu'il vaur rajouter a toutes les valeurs v0 = 1+ v[N-1] +v[N-2]+v[N-lg2]+..+ v[(N-1-lg2] -v[lg2]-v[lg2-1]...-v[2]
//  donc un poids x N.
// On stocke dans pond[0] l'offset constant( donc ss influence sur les comparaisons)
// entre la somme des deltas ponderes (SIGMA(pond[is]*delta[is]) et la vrai somme
sum103_t InitPondSum(sum103_t * pondDelta,int lgAsk , AlterPaths *AltP) {
    sum103_t v[PB103_NB] ;
    sum103_t offsetSum = 1 ;
    int i ;
    for(i=1;i<PB103_NB;i++) {
        pondDelta[i] = PB103_NB - i ;
        offsetSum += i+1 ;
    }
    int lg2 = (PB103_NB-1)/2 ;
    for(i=0;i<lg2;i++) {
        int k ;
        offsetSum +=(PB103_NB-1-i)  * PB103_NB ;
        for(k=1;k<=PB103_NB-1-i;k++) {
            pondDelta[k]+= PB103_NB ;
        }
        offsetSum -= (i+1) * PB103_NB ;
        for(k=1;k<=i+1;k++) {
            pondDelta[k]-= PB103_NB ;
        }
    }
    pondDelta[0] = offsetSum ; // inutilise
    //  for(i=0;i<PB103_NB;i++) { printf("%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    
    
    //on va initialiser avec une valeur pas trop loin de l'optimum
    v[0] = 0 ;  v[1] = 2 ; v[2] = 3 ;     v[3] = 4 ;
    switch(lgAsk){
        case 0: return 0;
        case 1: return v[0] ;
        case 2: return v[0] + v[1]  ;
        case 3: return v[0] + v[1] + v[2];
        case 4: return v[0] + v[1] + v[2]+ v[3];
    }
    int j, lg ;
    for(lg=3;;lg++) {
        sum103_t v0 = MinCheck(v,lg+1) ;
        for(j=0;j<=lg;j++) v[j] += v0 ;
        if(!CheckEqualityAlt(v,lg,AltP,0) ) {
            for(j=0;j<=lg;j++){ printf("%d%c",v[j],(j==lg) ? '\n' : ',' ) ; }
            return 0 ;
        }
        if(lg < lgAsk-1) {
            for(j=lg;j>=0;j--) v[j+1] = v[j] ;
            v[0] = 0 ;
        } else {
            break ;
        }
    }
    sum103_t Smin = 0 ;
    for(j=0;j<lgAsk;j++) Smin += v[j] ;
    return Smin ;
    
}


// version avec ponderation, mais le parcours se fait en tenant compte
// pour l'ordre de toutes les ponderations sauf la derniere valeur
// Du coup nombre de cas excede largement Smin et sont rejetes au dernier moment
// Utilise CheckEqualityPreN

int PB103a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    sum103_t values[PB103_NB] ;
    sum103_t Smin = 32000 ;
    sum103_t deltaMax = 20000 ;
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    CheckPaths *chkP =GetCheckPath(PB103_NB,AltP,0) ;
    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;
    
    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%3s Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }
    values[0] = 0 ;
    int j ;
    int is = 1 ;
    int deltaS = 0 ;
    sum103_t delta = deltaS ;
    sum103_t Delta[PB103_NB] ;
    for(j=0;j<PB103_NB-1;j++) Delta[j] = 0;
    Delta[PB103_NB-1] = deltaS ;
    while(deltaS<=deltaMax) {
        for(;is<PB103_NB;is++) {
            values[is] = values[is-1]+Delta[is]+1 ; /*S+= values[is] ;*/
        }
        sum103_t S = deltaS + Delta[PB103_NB-1] * ( pondDelta[PB103_NB-1]-1) ;
        if( S <= Smin &&  CheckEqualityPreN(values,PB103_NB,chkP))  {
            sum103_t v0 = MinCheck(values,PB103_NB) ;
            int j ;
            Smin = S ;
            deltaMax = S ;
            if(pbR->isVerbose) {
                fprintf(stdout,"\t PB%s S=%d Delta=%d DeltaMax=%d ",pbR->ident,Smin+offsetSum,deltaS,deltaMax ) ;
                for(j=0;j<PB103_NB;j++){
                    fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                }
                for(j=1;j<PB103_NB;j++){
                    fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                }
            }
            int lg = 0 ;
            for(j=0;j<PB103_NB;j++){
                lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
            }
        }
        is = PB103_NB-2 ;
        while(pondDelta[is] > delta && is){
            delta += pondDelta[is] * Delta[is];
            Delta[is-- ] = 0;
        }
        if(is == 0) {
            deltaS++ ;
            delta = deltaS  ;
            is = 1 ;
        } else {
            Delta[is]++ ; delta -= pondDelta[is] ;
        }
        Delta[PB103_NB-1] = delta ;
    }
    FreeCheckPath(chkP) ;
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



// version avec ponderation, mais le parcours se fait en tenant compte
// pour l'ordre de toutes les ponderations sauf la derniere valeur
// Mais on ne teste que si la derniere valeur est un multiple de la ponderation
// Donc economie de test par rapport a la version a)
// Ne sort que les solutions correctes (puisque parcourrues ds l'ordre croissant
// Utilise CheckEqualityPreN
int PB103b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    
    sum103_t values[PB103_NB] ;
    sum103_t Smin = 32000 ;
    sum103_t deltaMax = 20000 ;
    
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    CheckPaths *chkP =GetCheckPath(PB103_NB,AltP,0) ;
 
    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;
    
    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%sb Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }


    int j ;
    int is = 1 ;
    int deltaS = 0 ;
    values[0] = 0 ;
    sum103_t delta = deltaS ;
    sum103_t Delta[PB103_NB] ;
    for(j=0;j<PB103_NB-1;j++) Delta[j] = 0;
    Delta[PB103_NB-1] = deltaS ;
    while(deltaS<=deltaMax) {
        if(( Delta[PB103_NB-1] * pondDelta[PB103_NB-1])  == delta ){
            for(is=1;is<PB103_NB;is++) {
                values[is] = values[is-1]+Delta[is]+1 ; /*S+= values[is] ;*/
            }
            sum103_t S = deltaS  ;
            if( S <= Smin &&  CheckEqualityPreN(values,PB103_NB,chkP))  {
                sum103_t v0 = MinCheck(values,PB103_NB) ;
                int j ;
                Smin = S ;
                deltaMax = S ;
                if(pbR->isVerbose) {
                    fprintf(stdout,"\t PB%s S=%d Delta=%d DeltaMax=%d ",pbR->ident,Smin+offsetSum,deltaS,deltaMax ) ;
                    for(j=0;j<PB103_NB;j++){
                        fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                    }
                    for(j=1;j<PB103_NB;j++){
                        fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                    }
                }
                int lg = 0 ;
                for(j=0;j<PB103_NB;j++){
                    lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
                }
            }
        }
        is = PB103_NB-2 ;
        while(pondDelta[is] > delta && is){
            delta += pondDelta[is] * Delta[is];
            Delta[is-- ] = 0;
        }
        if(is == 0) {
            deltaS++ ;
            delta = deltaS  ;
            is = 1 ;
        } else {
            Delta[is]++ ; delta -= pondDelta[is] ;
        }
        Delta[PB103_NB-1] = delta / pondDelta[PB103_NB-1] ;
    }
    FreeCheckPath(chkP) ;
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
// version avec ponderation, avec parcours en tenant compte
// de l'ordre de toutes les ponderations a l'ordre PB103_NB-1
// Si c'est correct, on teste otus les cas possible en balayant les valeurs
// possible pour le dernier delta.
// Utilise CheckEqualityAlt pour PB103_NB-1,
// CheckEqualityPreN pour PB103_NB (avec isSup=1)


int PB103c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    sum103_t values[PB103_NB] ;
    sum103_t Smin = 32000 ;
    sum103_t deltaMax = 20000 ;
    
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    CheckPaths *chkP =GetCheckPath(PB103_NB,AltP,1) ;

    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;
    
    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%sc Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }
    int j, is = 1 ;
    int deltaS = 0 ;
    
    values[0] = 0 ;
    sum103_t delta = deltaS ;
    sum103_t Delta[PB103_NB] ;
    for(j=0;j<PB103_NB-1;j++) Delta[j] = 0;
    Delta[PB103_NB-1] = deltaS ;
    while(deltaS<=deltaMax) {
        if(delta == 0 ){
            for(is=1;is<PB103_NB-1;is++) {
                values[is] = values[is-1]+Delta[is]+1 ;
            }
            if(CheckEqualityAlt(values,PB103_NB-1 ,AltP,0)) {
                values[PB103_NB-1] = values[PB103_NB-2]+1 ;
                sum103_t S = deltaS  ;
                while(S<=Smin) {
                    if(CheckEqualityPreN(values,PB103_NB,chkP))  {
                        sum103_t v0 = MinCheck(values,PB103_NB) ;
                        int j ;
                        Smin = S ;
                        deltaMax = S ;
                        if(pbR->isVerbose) {
                            fprintf(stdout,"\t PB%s S=%d Delta=%d DeltaMax=%d ",pbR->ident,Smin+offsetSum,deltaS,deltaMax ) ;
                            for(j=0;j<PB103_NB;j++){
                                fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                            }
                            for(j=1;j<PB103_NB;j++){
                                fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                            }
                        }
                        int lg = 0 ;
                        for(j=0;j<PB103_NB;j++){
                            lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
                        }
                        break ;
                    }
                    values[PB103_NB-1] ++ ; S += pondDelta[PB103_NB-1] ;
                }
            }
        }
        is = PB103_NB-2 ;
        while(pondDelta[is] > delta && is){
            delta += pondDelta[is] * Delta[is];
            Delta[is-- ] = 0;
        }
        if(is == 0) {
            deltaS++ ;
            delta = deltaS  ;
            is = 1 ;
        } else {
            Delta[is]++ ; delta -= pondDelta[is] ;
        }
    }
    FreeCheckPath(chkP) ;
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
//
//
// Parcours en incrementant les delta[is) tt que somme pondere inferieur a dmin
// Et l'on ne developpe que si les ordres inferieurs sont valides
// Pas dans l'ordre de dmin croissant, mais a chque nouvelle valeur on tient compte de l'amelioration de dmin
// utilise CheckEqualityAlt pour is =3
// utilise CheckEqualityPreN aves isSup pour les ordres superieurs
//
int PB103d(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    
    sum103_t values[PB103_NB] ;
    sum103_t Smin = 32000 ;
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    CheckPaths * chkPS[PB103_NB] ;
    { int i ; for(i=4;i<PB103_NB;i++) chkPS[i] = GetCheckPath(i+1,AltP,1) ; }
  
    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;

    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%s Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }
        sum103_t S[PB103_NB] ;
        values[0] = 0 ; values[1] = 1 ; S[1] = 0 ;
        int is = 1 ;
        while(is>0) {
             if(S[is] <= Smin) {
                 if(is <= 3) {
                    if(is== 3 && !CheckEqualityAlt(values,is+1,AltP,0) ) {
                        values[is]++ ; S[is] += pondDelta[is] ;
                        continue ;
                    }
                } else if(is < PB103_NB) {
                    if(!CheckEqualityPreN(values,is+1,chkPS[is]) ) {
                        values[is]++ ; S[is] += pondDelta[is] ; continue;
                    } else if(is==PB103_NB-1) {
                        // new solution
                        sum103_t v0 = MinCheck(values,PB103_NB) ;
                        int j ;
                        Smin = S[is] ;
                        if(pbR->isVerbose) {
                            fprintf(stdout,"\t PB%s S=%d D=%d  ",pbR->ident,Smin+offsetSum,Smin ) ;
                            for(j=0;j<PB103_NB;j++){
                                fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                            }
                            for(j=1;j<PB103_NB;j++){
                                fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                            }
                        }
                        int lg = 0 ;
                        for(j=0;j<PB103_NB;j++){
                            lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
                        }
                        values[is]++ ; S[is] += pondDelta[is] ; continue;
                    }
                }
                is++ ;
                S[is] = S[is-1] ;  values[is] = values[is-1] + 1 ;
            } else {
                is-- ;
                values[is]++ ; S[is] += pondDelta[is] ;
            }
        }
    { int i ; for(i=4;i<PB103_NB;i++) FreeCheckPath(chkPS[i]);  }
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}





int PB103e(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    sum103_t values[PB103_NB] ;
    sum103_t Smin = 32000 ;
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    HalfPaths * hlfP[PB103_NB] ;
    { int i ; for(i=4;i<PB103_NB;i++) hlfP[i] = GetHalfPath(i+1) ; }
    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;
    //    Smin = 2773 - offsetSum ; // valeur optimum pour 10
    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%s Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }
    sum103_t S[PB103_NB] ;
    sum103_t SRev=0 ;
    values[0] = 0 ; values[1] = 1 ; S[1] = 0 ;
    int is = 1 ;
    while(is>0) {
        if(S[is] <= Smin || SRev <= Smin ) {
            if(is <= 3) {
                if(is== 3 && !CheckEqualityAlt(values,is+1,AltP,0) ) {
                    values[is]++ ; S[is] += pondDelta[is] ;
                    SRev += pondDelta[PB103_NB-is] ;
                    continue ;
                }
            } else if(is < PB103_NB) {
                if(!CheckEqualityPreH(values,hlfP[is]) ) {
                    values[is]++ ; S[is] += pondDelta[is] ;  SRev += pondDelta[PB103_NB-is] ;  continue;
                } else if(is==PB103_NB-1) {
                    int j ;
                    {// new solution
                        sum103_t v0 = MinCheck(values,PB103_NB) ;
                        if(S[is] <= Smin)  {
                            Smin = S[is] ;
                            if(pbR->isVerbose) {
                                fprintf(stdout,"\t PB%s S=%d D=%d  ",pbR->ident,Smin+offsetSum,Smin ) ;
                                for(j=0;j<PB103_NB;j++){
                                    fprintf(stdout,"%d%c",values[j]+v0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                                }
                                for(j=1;j<PB103_NB;j++){
                                    fprintf(stdout,"%d%c",values[j]-values[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                                }
                            }
                        }
                        // on va tester le reverse
                        if(SRev < Smin) {
                            sum103_t vRev[PB103_NB] ;
                            vRev[0]  = 0 ;
                            for(j=1;j<PB103_NB;j++){
                                vRev[j] = vRev[j-1] + values[PB103_NB-j]-values[PB103_NB-j-1] ;
                            }
                            {
                                sum103_t vR0 = MinCheck(vRev,PB103_NB) ;
                                if(CheckEqualityPreH(vRev,hlfP[is])) {
                                    fprintf(stdout,"\t PB%s [R] S=%d D=%d  ",pbR->ident,SRev+offsetSum,SRev ) ;
                                    for(j=0;j<PB103_NB;j++){
                                        fprintf(stdout,"%d%c",vRev[j]+vR0,(j==PB103_NB-1) ? ' ' : ',' ) ;
                                    }
                                    for(j=1;j<PB103_NB;j++){
                                        fprintf(stdout,"%d%c",vRev[j]-vRev[j-1]-1,(j==PB103_NB-1) ? '\n' : '.' ) ;
                                    }
  
                                    Smin = SRev ;
                                }
                            }
                        }

                        int lg = 0 ;
                        for(j=0;j<PB103_NB;j++){
                            lg+=sprintf(pbR->strRes+lg,"%2.2d",values[j]+v0) ;
                        }
                    }
                    values[is]++ ; S[is] += pondDelta[is] ; SRev += pondDelta[PB103_NB-is] ; continue;
                }
            }
            is++ ;
            S[is] = S[is-1] ;  values[is] = values[is-1] + 1 ;
        } else {
            SRev -= (values[is] - values[is-1] - 1) * pondDelta[PB103_NB-is] ;
            is-- ;
            values[is]++ ; S[is] += pondDelta[is] ; SRev += pondDelta[PB103_NB-is] ;
        }
    }
    { int i ; for(i=4;i<PB103_NB;i++) FreeHalfPath(hlfP[i]) ;  }
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int Cmp103f(const void *el1,const void *el2) {
    static int size = 0;
    if(el2==NULL) {
        size = ((int *)el1)[0] ;
        return size ;
    }
    return memcmp(el1,el2,size) ;
}

typedef struct DevDeltas {
    sum103_t val[PB103_NB] ;
    sum103_t deltas[PB103_NB] ;
    sum103_t R_val[PB103_NB] ;
    sum103_t R_deltas[PB103_NB] ;
    sum103_t S ;
    sum103_t Srev ;
    sum103_t R_S ;
    sum103_t R_Srev ;
    int     nbDev ;
    int     R_nbDev ;
    int     isR ;
} DevDeltas ;

void InitLevel(DevDeltas  *DD, int is, sum103_t pondDelta[PB103_NB]) {
    int j ;
    DD[is].val[0] = DD[is].R_val[0] = 0 ;
    DD[is].S=0;
    DD[is].Srev=0 ;
    DD[is].R_S=0;
    DD[is].R_Srev=0 ;
    DD[is].nbDev = DD[is].R_nbDev = 0 ;
 
 //       DD[is].deltas[is-1] = 0 ;
    DD[is].deltas[is-1] = DD[is].deltas[is-2] ;
//   DD[is].val[is] = DD[is].val[is-1] + DD[is].deltas[is-1] + 1 ;
    
//    DD[is].R_deltas[is-1] = DD[is].R_deltas[is-2] ;
//    DD[is].R_val[is] = DD[is].R_val[is-1] +DD[is].R_deltas[is-1] + 1 ;
    
    
    for(j=1;j<=is;j++) {
        DD[is].S += DD[is].deltas[j-1] *pondDelta[j] ;
        DD[is].Srev += DD[is].deltas[j-1] *pondDelta[is+1-j] ;
        DD[is].val[j] = DD[is].val[j-1] + DD[is].deltas[j-1] + 1 ;
        
        if(j<is) DD[is].R_deltas[j-1] = DD[is].deltas[is-1-j] ;
        else DD[is].R_deltas[j-1] = DD[is].R_deltas[j-2] ;
        DD[is].R_S += DD[is].R_deltas[j-1] *pondDelta[j] ;
        DD[is].R_Srev += DD[is].R_deltas[j-1] *pondDelta[is+1-j] ;
        DD[is].R_val[j] = DD[is].R_val[j-1] + DD[is].R_deltas[j-1] + 1 ;
    }
    
    DD[is].isR= 0 ;
    return ;
}

#define MAX_DEV 100
int PB103f(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
//    sum103_t values[PB103_NB] ;
    sum103_t vRev[PB103_NB] ;
    
    sum103_t Smin = 32000 ;
    sum103_t pondDelta[PB103_NB] ;
    AlterPaths *AltP =GetAlterPath(PB103_NB/2) ;
    HalfPaths * hlfP[PB103_NB] ;
    { int i ; for(i=4;i<PB103_NB;i++) hlfP[i] = GetHalfPath(i+1) ; }
    Smin = InitPondSum(pondDelta,PB103_NB ,AltP) ;
    sum103_t offsetSum = pondDelta[0]  ;
    if(pbR->isVerbose) {
        int i ;
        fprintf(stdout,"\t PB%s Smin=%d,Dmin=%d Pond=",pbR->ident, Smin+offsetSum,Smin) ;
        for(i=0;i<PB103_NB;i++) { fprintf(stdout,"%d%c",pondDelta[i],(i==PB103_NB-1) ? '\n' : ' '); }
    }
    
    sum103_t deltas0[PB103_NB] = {0,0,1,2} ;
    
    DevDeltas  DD [PB103_NB] ;
    int is,j ;
    is = 4 ;
     for(j=0;j<is-1;j++) {
        DD[is].deltas[j] = deltas0[j] ;
    }
    InitLevel(DD,is,pondDelta) ;
    
    while(is>3) {
        // tritement d'un element
/*        printf("%d,%d %d,%d,%d,%d V=",is,DD[is].isR,DD[is].S,DD[is].Srev,DD[is].R_S,DD[is].R_Srev) ;
        for(j=0;j<=is;j++) printf("%d%c",DD[is].val[j],(j==is) ? ' ': ',') ;
        for(j=0;j<is;j++) printf("%d%c",DD[is].deltas[j],(j==is-1) ? ' ': '.') ;
        printf("R_V=");
        for(j=0;j<=is;j++) printf("%d%c",DD[is].R_val[j],(j==is) ? ' ': ',') ;
        for(j=0;j<is;j++) printf("%d%c",DD[is].R_deltas[j],(j==is-1) ? '\n': '.') ;
 */
    if(DD[is].isR== 0) {
 //           DD[is].isR= 1 ;
            if(/* DD[is].nbDev < MAX_DEV && */ (DD[is].S <= Smin || DD[is].Srev <= Smin)) {
                if(DD[is].nbDev || CheckEqualityPreH(DD[is].val,hlfP[is])) {
                    DD[is].nbDev++ ;
                    if(is == PB103_NB -1 ) {
                        // nouvelle solution a base de val
                        if(  DD[is].S <= Smin) {
                            Smin = DD[is].S ;
                            sum103_t v0 = MinCheck(DD[is].val,PB103_NB) ;
                            printf("S=%d,%d ",DD[is].S+offsetSum,DD[is].S) ;
                            for(j=0;j<=is;j++) printf("%d%c",DD[is].val[j]+v0,(j==is) ? ' ' : ',') ;
                            for(j=1;j<=is;j++) printf("%d%c",DD[is].deltas[j-1],(j==is) ? '\n' : '.') ;
                        }
                        if( DD[is].Srev <= Smin) {
                            Smin = DD[is].Srev ;
                            vRev[0] = 0;
                            for(j=1;j<=is;j++) vRev[j] = vRev[j-1] + DD[is].deltas[is-j]+ 1 ;
                            sum103_t v0 = MinCheck(vRev,PB103_NB) ;
                            printf("S[r]=%d,%d ",DD[is].Srev+offsetSum,DD[is].Srev) ;
                            for(j=0;j<=is;j++) printf("%d%c",vRev[j]+v0,(j==is) ? ' ' : ',') ;
                            for(j=1;j<=is;j++) printf("%d%c",DD[is].deltas[is-j],(j==is) ? '\n' : '.') ;
                        }
                        DD[is].isR= 1 ;
                        continue ;
                        
                    } else {
                        memcpy(DD[is+1].deltas,DD[is].deltas,is*sizeof(DD[0].deltas[0]) ) ;
                        DD[is].val[is]++ ; DD[is].deltas[is-1]++ ; DD[is].S += pondDelta[is] ; DD[is].Srev += pondDelta[1] ;
                        is ++ ; // on passer au niveau superieur
                        InitLevel(DD, is, pondDelta) ;
                        continue ;
                    }
                }
                // on
                DD[is].val[is]++ ; DD[is].deltas[is-1]++ ; DD[is].S += pondDelta[is] ; DD[is].Srev += pondDelta[1] ;
                continue ;
            }
            DD[is].isR= 1 ;
        }
        // on est dans le cas isR==1
//        DD[is].isR= 0 ;
        if(/* DD[is].R_nbDev < MAX_DEV && */ (DD[is].R_S <= Smin || DD[is].R_Srev <= Smin)){
            if(DD[is].R_nbDev || CheckEqualityPreH(DD[is].R_val,hlfP[is])) {
                DD[is].R_nbDev++ ;
                if(is == PB103_NB -1 ) {
                    // nouvelle solution a base de R_val
                    if(  DD[is].R_S <= Smin) {
                        Smin = DD[is].R_S ;
                        sum103_t v0 = MinCheck(DD[is].R_val,PB103_NB) ;
                        printf("*S=%d,%d ",DD[is].R_S+offsetSum,DD[is].R_S) ;
                        for(j=0;j<=is;j++) printf("%d%c",DD[is].R_val[j]+v0,(j==is) ? ' ' : ',') ;
                        for(j=1;j<=is;j++) printf("%d%c",DD[is].R_deltas[j-1],(j==is) ? '\n' : '.') ;
                    }
                    if( DD[is].R_Srev <= Smin) {
                        Smin = DD[is].R_Srev ;
                        vRev[0] = 0;
                        for(j=1;j<=is;j++) vRev[j] = vRev[j-1] + DD[is].R_deltas[is-j]+ 1 ;
                        sum103_t v0 = MinCheck(vRev,PB103_NB) ;
                        printf("*S[r]=%d,%d ",DD[is].R_Srev+offsetSum,DD[is].R_Srev) ;
                        for(j=0;j<=is;j++) printf("%d%c",vRev[j]+v0,(j==is) ? ' ' : ',') ;
                        for(j=1;j<=is;j++) printf("%d%c",DD[is].R_deltas[is-j],(j==is) ? '\n' : '.') ;
                    }
                    is-- ;
                    continue ;
                } else {
                    memcpy(DD[is+1].deltas,DD[is].R_deltas,is*sizeof(DD[0].deltas[0]) ) ;
                    DD[is].R_val[is]++ ; DD[is].R_deltas[is-1]++ ; DD[is].R_S += pondDelta[is] ; DD[is].R_Srev += pondDelta[1] ;
                    is ++ ; // on passe au niveau superieur
                    InitLevel(DD, is, pondDelta) ;
                    continue ;

                }
             }
            // score admissible on reste au meme niveau
            DD[is].R_val[is]++ ; DD[is].R_deltas[is-1]++ ; DD[is].R_S += pondDelta[is] ; DD[is].R_Srev += pondDelta[1] ;
            continue ;
            
        } else /* if (DD[is].S > Smin && DD[is].Srev > Smin) */ {
            // on est  dans le cas ou aucunes conditions de score n'est satisfaite
            is-- ;
        }
    }
    
    
    { int i ; for(i=4;i<PB103_NB;i++) FreeHalfPath(hlfP[i]) ;  }
    FreeAlterPath(AltP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int Cmpint16(const void *el1, const void *el2) {
    return ((int16_t *)el1)[0] - ((int16_t *)el2)[0] ;
}

int PB105(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int32_t * * tbSub = P105_GetData() ;
    int nt ;
    u_int32_t S =0 ;
    for(nt=0;tbSub[nt] != 0; nt++){
        int32_t * sub = tbSub[nt] ;
        int lg ;
        for(lg=0;sub[lg]!=0;lg++) ;
        heapsort(sub,lg,sizeof(sub[0]),Cmpint16) ;
        if(sub[0] >= MinCheck(sub,lg)  && CheckEquality(sub,lg) ) {
            int i;
            for(i=0;i<lg;i++) S += sub[i] ;
//            printf("%d ",nt);
        }
    }
    sprintf(pbR->strRes,"%d",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB106_ASK   12
int PB106(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    AlterPaths * altP = GetAlterPath(PB106_ASK/2) ;
    CheckPaths * chkP = GetCheckPath(PB106_ASK,altP,0) ;
    int S =0;
    int k ;
    for(k=2 ; 2*k <= PB106_ASK ; k++) {
        S += (chkP->nsumK[k-1]-1)*chkP->npermK[k-1];
    }
    sprintf(pbR->strRes,"%d",S) ;
    FreeAlterPath(altP);
    FreeCheckPath(chkP);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

