//
//  PB_gmp.c
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "PB_gmp.h"

#define PB016_MAXL  1000/3
#define PB016_EXP   1000

int PB016_gmp(PB_RESULT *pbR) {
    char  digLarge[PB016_MAXL] ;
    u_int32_t S = 0 ;
    int i ;
    pbR->nbClock = clock() ;
    mpz_t pow2 ;
    mpz_init(pow2);
    mpz_ui_pow_ui (pow2, 2, PB016_EXP);
    gmp_sprintf (digLarge,"%Zd",pow2);
    for(i=0;digLarge[i] != 0; i++) {
        S += digLarge[i] - '0' ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp Sumdig(2**%d)=%d\n"
                              ,pbR->pbNum
                              ,PB016_EXP,S) ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
}

#define PB056_MAX   100
#define PB056_MAXL  201
int PB056_gmp(PB_RESULT *pbR) {
    char  digLarge[PB056_MAXL] ;
    u_int32_t Max = 0 ;
    int a,b ;
    int aBest =0, bBest = 0 ;
    pbR->nbClock = clock() ;
    mpz_t pow ;
    mpz_init(pow);
    for(b=PB056_MAX-1;b>1;b--) {
        if(18*b <= Max) break ; // car max 2*b digits
        for(a=PB056_MAX-1;a>1;a--) {
            int i , S = 0 ;
            mpz_ui_pow_ui (pow, a, b);
            gmp_sprintf (digLarge,"%Zd",pow);
            if(strlen(digLarge)*9 <= Max ) continue ;
            for(i=0;digLarge[i] != 0; i++) {
                S += digLarge[i] - '0' ;
            }
            if(S > Max) {
                Max = S ;
                aBest = a;
                bBest = b ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp DigSum(%d**%d)=%d\n"
                              ,pbR->pbNum
                              ,aBest,bBest,Max) ;
    sprintf(pbR->strRes,"%d",Max) ;
    return 1 ;
}

int PB066(PB_RESULT *pbR) {
    int32_t N , n , d ;
    int32_t k0, k02 ;
    pbR->nbClock = clock()  ;
    mpz_t max_x ;
    mpz_init_set_ui (max_x, 0) ;
    int bestN = 0 ;
    FractContG * FCG = FCG_alloc() ;
    for(N=2,k0=1,k02=4;N<=1000;N++) {
        int32_t i,a;
        if(N == k02) { // k02 = (k0+1)*(k0+1)
            k0++ ;
            k02 += 2*k0 + 1 ; continue ;
        }
        n = k0 ; d=1 ; i = 0 ; // so k0 =(int) srqt(N)
        FCG_init(FCG,k0) ;
        do {
            d = (N - n * n) / d ; // quite easy to recursively show that division is exact
            a = (k0+n)/d ;
            n = k0 - ( (k0 + n) % d ) ;
            i++ ;
            FCG_NextCoef(FCG,a) ;
        }while(d !=1 || n != k0 || (i&1)) ; // test loop on (n,d) = (k0,1) et i even (solution of PELL equation)
        if(FCG_CmpNum0(FCG,max_x) > 0 ) { // ask for n0 because one step ahead
            FCG_GetNum0(FCG,max_x) ;
            bestN = N ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) {
        char * str_x = mpz_get_str(NULL,10,max_x) ;
        fprintf(stdout,"\t PB%0.3d x(%d)=%s\n",pbR->pbNum, bestN,str_x);
        free (str_x) ;
    }
    sprintf(pbR->strRes,"%d",bestN);
    return 1 ;
}



#define PB080_NBD   100
#define PB080_N     100

int PB080_gmp(PB_RESULT *pbR) {
    char  digDecimal[PB080_NBD+10] ;
    int S=0,i,k ,nxt_k2 ;
    pbR->nbClock = clock() ;
    mpz_t n2 ;
    mpz_init(n2);
    for(i=1,k=1,nxt_k2=1;i<PB080_N;i++) {
        if(i == nxt_k2) {
            nxt_k2 += 2*k + 1 ;
            k++ ;
            continue ;
        }
        mpz_ui_pow_ui (n2,10,(PB080_NBD-1)*2);
        mpz_mul_si (n2,n2,i) ;
        mpz_sqrt (n2,n2) ;
        gmp_sprintf (digDecimal,"%Zd",n2);
        {
            int j ;
            for(j=0;digDecimal[j] != 0; j++) {
                S += digDecimal[j] - '0' ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%dgmp SumdDecimal[1..%d]=%d\n"
                              ,pbR->pbNum,PB080_N,S) ;
    sprintf(pbR->strRes,"%d",S) ;
    return 1 ;
}



#define PB597_LG    1800
#define PB597_NB    13
#define PB597_SIZEB 40
#define PB597_DEBUG 0

typedef struct devTorpids_a {
    mpq_t   P ;                   // proba
    u_int32_t   bitEnv ;
    u_int16_t   isEven ;
    
} devTorpids_a ;

#if PB597_DEBUG
static void PB597_printTorpids_a(FILE *fout,devTorpids_a * dT) {
    gmp_fprintf(fout,"%Qd %c:",dT->P,dT->isEven ? 'E' : 'O') ;
    int bitTest,ndec ;
    for(bitTest=1,ndec=0;bitTest<=dT->bitEnv;bitTest <<=1 , ndec++) { if(bitTest & dT->bitEnv) fprintf(fout,"%c%c", (bitTest==1) ? '{' : ',','A' + ndec) ; }
    fprintf(fout,"}\n");
}
#endif

#define MPQ_INIT_ND(mpq,N,D)   mpq_set_ui(mpq,N,D) ; mpq_canonicalize(mpq)


int PB597_gmpa(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    u_int32_t nb, curBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    // calcul factoriel n!
    for(i=1,nb=1;i<=PB597_NB;i++) nb *= i ;
    devTorpids_a * antDT = malloc(nb*sizeof(antDT[0])) ;
    devTorpids_a * curDT = malloc(nb*sizeof(curDT[0])) ;
    u_int32_t antNb, curNb ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    {
        u_int64_t pg = PGCD64(curLg,sizeB);
        curLg /= pg ;
        sizeB /= pg ;
    }
    curNb= 0 ;
    mpq_t gcd, ND,ONE ;
    mpq_init(gcd);
    mpq_init(ND);
    mpq_init(ONE) ;
    mpq_set_ui(ONE,1,1);
    
    mpq_init(curDT[curNb].P) ;
    mpq_set_ui(curDT[curNb].P,1,1) ;
    curDT[curNb].isEven = 1;
    curDT[curNb].bitEnv = 1 ;
    sumLg[0] = Lg[0] = curLg ;
    curNb++ ;
    for(curBoat=1;curBoat<PB597_NB;curBoat++) {
        devTorpids_a *tmp = antDT ;
        antDT = curDT ;
        curDT = tmp ;
        antNb = curNb ;
        curNb = 0 ;
        curLg += sizeB ;
        Lg[curBoat] = curLg ;
        sumLg[curBoat] = sumLg[curBoat-1] + curLg ;
#if PB597_DEBUG
        fprintf(stdout,"==%d\n",curBoat);
#endif
        for(i=0;i<antNb;i++) {
            devTorpids_a * ptAntDT = antDT + i ;
            MPQ_INIT_ND(ND,curLg,sumLg[curBoat]) ;
            // on rajoute le bateau qui ne bumpe personne
            mpq_init (curDT[curNb].P) ;
            mpq_mul (curDT[curNb].P, ptAntDT->P, ND) ;
            curDT[curNb].bitEnv = 1 << curBoat ;
            curDT[curNb].isEven = ptAntDT->isEven ;
#if PB597_DEBUG
            printf("*");
            PB597_printTorpids_a(stdout,curDT+curNb);
#endif
            curNb++ ;
            int place ;
            u_int16_t nboat ;
            // on va maintenant rajouter les bumps
            for(place=0,nboat=0;(1<< nboat) <= ptAntDT->bitEnv ;place++) {
                while(((1<< nboat) & ptAntDT->bitEnv) == 0) {   nboat++ ; }
                mpq_init(curDT[curNb].P) ;
                if(place == 0) {
                    MPQ_INIT_ND(ND,sumLg[curBoat] - curLg,sumLg[curBoat]) ;
                    mpq_mul(curDT[curNb].P,ptAntDT->P,ND);
                } else {
                    mpq_inv(ND,ND) ; // ND = D/N
                    mpq_sub(ND,ND,ONE); // ND = D/N -1 = (D-N)/N
                    mpq_mul(curDT[curNb].P,curDT[curNb-1].P,ND);
                }
                MPQ_INIT_ND(ND,(curBoat - nboat) * sizeB,sumLg[curBoat]-sumLg[nboat] - (curBoat - nboat) * Lg[nboat]) ;
                mpq_mul(curDT[curNb].P,curDT[curNb].P,ND);
                nboat++ ;
                curDT[curNb].bitEnv =  (ptAntDT->bitEnv & ((1 << nboat) -1)) |  (1 << curBoat)  ;
                curDT[curNb].isEven = (place & 1) ? ptAntDT->isEven : !ptAntDT->isEven ;
#if PB597_DEBUG
                PB597_printTorpids_a(stdout,curDT+curNb);
#endif
                curNb++ ;
            }
        }
    }
    {
        mpq_t ND ;
        mpq_init(ND);
        mpq_set_ui(ND,0,1);
        for(i=0;i<curNb;i++) {
            if(curDT[i].isEven) {
                mpq_add(ND,ND,curDT[i].P) ;
            }
        }
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%da PEVEN=%Qd\n",pbR->pbNum,ND );
        mpz_mul_ui(mpq_numref(ND) ,mpq_numref(ND),100000000000) ;
        mpz_tdiv_q (mpq_numref(ND),mpq_numref(ND),mpq_denref(ND));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(ND) );
    }
    free(curDT); free(antDT) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


typedef struct devTorpids {
    u_int16_t   isUsed ;
    mpq_t   E_P ;                   // proba for even permutations
    mpq_t   O_P ;                   // proba for odd permutations
} devTorpids ;

#if PB597_DEBUG
static void PB597_printTorpids(FILE *fout,devTorpids * dT,u_int32_t index) {
    gmp_fprintf(fout,"[%d]E:%Qd O:%Qd",index,dT->E_P,dT->O_P) ;
    int bitTest,kboat ;
    for(bitTest=1,kboat=0;bitTest<=index;bitTest <<=1 , kboat++) { if(bitTest & index) fprintf(fout,"%c%c", (bitTest==1) ? '{' : ',','A' + kboat) ; }
    
    fprintf(fout,"}\n");
}
#endif
#define MPQ_ADD_EO(dst,src1,src2)   mpq_add((dst)->E_P,(src1)->E_P,(src2)->E_P) ; mpq_add((dst)->O_P,(src1)->O_P,(src2)->O_P)
#define MPQ_SET_ND_EO(dst,ND)       mpq_set((dst)->E_P,ND) ; mpq_set((dst)->O_P,ND)
#define MPQ_INIT_SET_EO(dst,src)       mpq_init((dst)->E_P); mpq_set((dst)->E_P,(src)->E_P) ;mpq_init((dst)->O_P); mpq_set((dst)->O_P,(src)->O_P)
#define MPQ_MUL_ND_EO(dst,src1,ND)   mpq_mul((dst)->E_P,(src1)->E_P,(ND)) ; mpq_mul((dst)->O_P,(src1)->O_P,(ND))
// crossing odd and even
#define MPQ_MUL_ND_ExO(dst,src1,ND)   mpq_mul((dst)->E_P,(src1)->O_P,(ND)) ; mpq_mul((dst)->O_P,(src1)->E_P,(ND))
#define MPQ_MUL_EO(dst,src1,src2)   mpq_mul((dst)->E_P,(src1)->E_P,(src2)->E_P) ; mpq_mul((dst)->O_P,(src1)->O_P,(src2)->O_P)
#define MPQ_SELFMUL_ND_ExO(dst,ND,tmp)   mpq_mul(tmp,(dst)->O_P,(ND)) ; mpq_mul((dst)->O_P,(dst)->E_P,(ND)) ; mpq_set((dst)->E_P,tmp)

#define MPQ_INIT_MUL_ND_EO(dst,src1,ND)   mpq_init((dst)->E_P); mpq_mul((dst)->E_P,(src1)->E_P,(ND)) ;mpq_init((dst)->O_P); mpq_mul((dst)->O_P,(src1)->O_P,(ND))
#define MPQ_ADD_MUL_ND_EO(dst,src1,ND,tmp)    mpq_mul(tmp,(src1)->E_P,(ND)) ; mpq_add((dst)->E_P,(dst)->E_P,tmp) ; mpq_mul(tmp,(src1)->O_P,(ND)) ; mpq_add((dst)->O_P,(dst)->O_P,tmp)


int PB597_gmp(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    u_int32_t  curBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    mpq_t tmpz , ND ,ONE ;
    devTorpids * DT = calloc(1 << (PB597_NB+1),sizeof(DT[0])) ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    devTorpids newDT ;
    mpq_init(ND);    mpq_init(tmpz);
    mpq_init(ONE) ;  mpq_set_ui(ONE,1,1);
    mpq_init(newDT.E_P) ; mpq_init(newDT.O_P) ;
    {
        u_int64_t pg = PGCD64(curLg,sizeB);
        curLg /= pg ;
        sizeB /= pg ;
    }
    {
        u_int32_t indNb= 1 ; // 1 <<0
        DT[indNb].isUsed = 1 ;
        mpq_init(DT[indNb].E_P);   mpq_set_ui(DT[indNb].E_P,1,1);
        mpq_init(DT[indNb].O_P);   mpq_set_ui(DT[indNb].O_P,0,1);
    }
    sumLg[0] = Lg[0] = curLg ;
    for(curBoat=1;curBoat<PB597_NB;curBoat++) {
        u_int32_t indAnt ;
        curLg += sizeB ;
        Lg[curBoat] = curLg ;
        sumLg[curBoat] = sumLg[curBoat-1] + curLg ;
#if PB597_DEBUG
        fprintf(stdout,"==%d\n",curBoat);
#endif
        for(indAnt=1<<(curBoat-1);indAnt< (1<<curBoat);indAnt++) {
            devTorpids * ptAntDT = DT + indAnt ;
            if(ptAntDT->isUsed == 0) continue ;
            MPQ_INIT_ND(ND,curLg,sumLg[curBoat]) ;
            devTorpids * ptCurDT = DT + (1 << curBoat) ;
            // on rajoute le bateau qui ne bumpe personne
            if(ptCurDT->isUsed) {
                MPQ_ADD_MUL_ND_EO(ptCurDT,ptAntDT,ND,tmpz) ; // ptCurDT += ptAntDT * N/D
            } else {
                ptCurDT->isUsed = 1 ;
                MPQ_INIT_MUL_ND_EO(ptCurDT,ptAntDT,ND) ; // ptCurDT = ptAntDT * N/D
            }
#if PB597_DEBUG
            printf("*");
            PB597_printTorpids(stdout,ptCurDT,1 << curBoat);
#endif
            int isFirst ;
            u_int16_t nboat ;
            //
            for(isFirst=1,nboat=0;(1<< nboat) <= indAnt ;) {
                while(((1<< nboat) & indAnt) == 0) {   nboat++ ; }
                if(isFirst) { // Begin ?
                    MPQ_INIT_ND(ND,sumLg[curBoat] - curLg,sumLg[curBoat]) ;
                    MPQ_MUL_ND_ExO(&newDT,ptAntDT,ND) ; // newDT = ptAndtDT * N/D with crossing Even/Odd
                    isFirst = 0 ;
                } else { // on croise even et odd
                    mpq_inv(ND,ND) ; // ND = D/N
                    mpq_sub(ND,ND,ONE); // ND = D/N -1 = (D-N)/N
                    MPQ_SELFMUL_ND_ExO(&newDT,ND,tmpz) ; // newtDT *= N/D with crossing Even/Odd
                }
                MPQ_INIT_ND(ND,(curBoat - nboat) * sizeB,sumLg[curBoat]-sumLg[nboat] - (curBoat - nboat) * Lg[nboat]) ;
                MPQ_MUL_ND_EO(&newDT,&newDT,ND) ; // newDT *= N/D
                nboat++ ;
                ptCurDT = DT + ( (indAnt & ((1 << nboat) -1)) |  (1 << curBoat) ) ;
                if(ptCurDT->isUsed) {
                    MPQ_ADD_EO(ptCurDT,ptCurDT,&newDT);  // ptCurDT += newDT
                } else {
                    ptCurDT->isUsed = 1 ;
                    MPQ_INIT_SET_EO(ptCurDT,&newDT) ; // ptCurDT = newDT
                }
#if PB597_DEBUG
                PB597_printTorpids(stdout,ptCurDT,  ( (indAnt & ((1 << nboat) -1)) |  (1 << curBoat) ) );
#endif
            }
        }
    }
    {
        mpq_t ND ;
        mpq_init(ND);
        mpq_set_ui(ND,0,1);
        for(i= (1 << (PB597_NB-1));i< ( 1<< (PB597_NB)) ;i++) {
            if(DT[i].isUsed) {
                mpq_add(ND,ND,DT[i].E_P) ;
            }
        }
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%d PEVEN=%Qd\n",pbR->pbNum,ND );
        mpz_mul_ui(mpq_numref(ND) ,mpq_numref(ND),100000000000) ;
        mpz_tdiv_q (mpq_numref(ND),mpq_numref(ND),mpq_denref(ND));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(ND) );
    }
    free(DT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB597_gmpx(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int32_t  nbBoat  ;
    u_int32_t curLg, sizeB ;
    u_int32_t Lg[PB597_NB] ;
    u_int32_t sumLg[PB597_NB] ;
    mpq_t tmpz , ND ,ONE, EVEN , FACT ;
    devTorpids * DT = calloc(1 << (PB597_NB),sizeof(DT[0])) ;
    curLg = PB597_LG - (PB597_NB -1) * PB597_SIZEB ;
    sizeB = PB597_SIZEB ;
    devTorpids  tmpDT ;
    mpq_init(ND);      mpq_init(tmpz); mpq_init(FACT) ;
    mpq_init(ONE) ;    mpq_set_ui(ONE,1,1);
    mpq_init(EVEN);    mpq_set_ui(EVEN,0,1) ;
    mpq_init(tmpDT.E_P) ; mpq_init(tmpDT.O_P) ;
    u_int64_t pg = PGCD64(curLg,sizeB);
    curLg /= pg ;  sizeB /= pg ; // to avoid big numbers
    DT[1].isUsed = 1 ; // init for race with one boat
    mpq_init(DT[1].E_P);   mpq_set_ui(DT[1].E_P,1,1);
    mpq_init(DT[1].O_P);   mpq_set_ui(DT[1].O_P,0,1);
    sumLg[0] = Lg[0] = curLg ;
    for(nbBoat=1;nbBoat<PB597_NB;nbBoat++) { // add boat by boat
        u_int32_t indAnt ;
        curLg += sizeB ;
        Lg[nbBoat] = curLg ; // length to run for new boat
        sumLg[nbBoat] = sumLg[nbBoat-1] + curLg ; // sigma length (usefull for coefficients)
#if PB597_DEBUG
        if(nbBoat != PB597_NB-1) fprintf(stdout,"==%d\n",nbBoat);
#endif
        // loop on states for race with curBoat - 1
        for(indAnt=1<<(nbBoat-1);indAnt< (1<<nbBoat);indAnt++) {
            devTorpids * ptAntDT = DT + indAnt ;
            if(ptAntDT->isUsed == 0) continue ;
            MPQ_INIT_ND(ND,curLg,sumLg[nbBoat]) ;
            devTorpids * ptCurDT = DT + (1 << nbBoat) ;
            if(nbBoat == PB597_NB-1) {// last loop cumulate only even
                mpq_mul(tmpz,ptAntDT->E_P,ND) ;
                mpq_add(EVEN,EVEN,tmpz) ;
            } else {
                if(ptCurDT->isUsed) {
                    MPQ_ADD_MUL_ND_EO(ptCurDT,ptAntDT,ND,tmpz) ; // ptCurDT += ptAntDT * N/D
                } else {
                    ptCurDT->isUsed = 1 ;
                    MPQ_INIT_MUL_ND_EO(ptCurDT,ptAntDT,ND) ; // ptCurDT = ptAntDT * N/D
                }
            }
#if PB597_DEBUG
            if(nbBoat != PB597_NB-1) {
                printf("*");
                PB597_printTorpids(stdout,ptCurDT,1 << nbBoat);
            }
#endif
            int nplace ;
            int indDT ;
            u_int64_t  nFACT, dFACT, pg ;
            nFACT = sumLg[nbBoat] - curLg ;
            dFACT = sumLg[nbBoat] ;
            u_int16_t kboat ;
            // loop on kboat for the precedent state
            // newboat takes his place
            for(nplace=0,kboat=0; ;nplace++) {
                while(((1<< kboat) & indAnt) == 0) {   kboat++ ; }
                u_int32_t N=(nbBoat - kboat) * sizeB ;
                u_int32_t D= sumLg[nbBoat]-sumLg[kboat] - (nbBoat - kboat) * Lg[kboat] ;
                kboat++ ;
                indDT =  (indAnt & ((1 << kboat) -1)) |  (1 << nbBoat) ;
                nFACT *= N ;
                dFACT *= D ;
                pg = PGCD64(nFACT,dFACT) ;
                if(pg > 1) { nFACT /= pg ; dFACT /= pg ; }
                ptCurDT = DT + indDT  ;
                MPQ_INIT_ND(FACT,nFACT,dFACT);
                if(nbBoat == PB597_NB-1) {// last loop cumulate only even
                    if(nplace & 1) {
                        mpq_mul(tmpz,ptAntDT->E_P,FACT);
                    }else {
                        mpq_mul(tmpz,ptAntDT->O_P,FACT);
                    }
                    mpq_add(EVEN,EVEN,tmpz) ;
                } else {
                    if(nplace & 1) {
                        MPQ_MUL_ND_EO(&tmpDT,ptAntDT,FACT) ;
                    } else {
                        MPQ_MUL_ND_ExO(&tmpDT,ptAntDT,FACT) ;
                    }
                    if(ptCurDT->isUsed) {
                        MPQ_ADD_EO(ptCurDT,ptCurDT,&tmpDT);
                    } else {
                        ptCurDT->isUsed = 1 ;
                        MPQ_INIT_SET_EO(ptCurDT,&tmpDT) ;
                    }
                }
#if PB597_DEBUG
                if(nbBoat != PB597_NB-1) PB597_printTorpids(stdout,ptCurDT,   indDT  );
#endif
                if((1<< kboat) <= indAnt) {
                    nFACT *= (D-N) ;
                    dFACT *= N ;
                } else {
                    break ;
                }
            }
        }
    }
    {
        if(pbR->isVerbose)gmp_fprintf (stdout,"\t PB%d %d boats PEVEN=%Qd\n",pbR->pbNum,PB597_NB,EVEN );
        mpz_mul_ui(mpq_numref(EVEN) ,mpq_numref(EVEN),100000000000) ;
        mpz_tdiv_q (mpq_numref(EVEN),mpq_numref(EVEN),mpq_denref(EVEN));
        gmp_sprintf (pbR->strRes,"%Zd",mpq_numref(EVEN) );
    }
    free(DT);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
typedef struct CTX_597 {
    mpq_t   * backEven ;    // for even(k,k) k=0..NB
    mpq_t   * frontEven ;   // for even(n,L-NB+n) n=0..NB
} CTX_597 ;

mpq_t * getEven(CTX_597 * ctx,int n,int lg) {
    return (n==lg) ? ctx->backEven : ctx->frontEven ;
}
void computePeven(CTX_597 * ctx,int n,int lg) {
    static int isInit=0 ;
    mpq_t   tmpS, tmpP ,tmpQ, ONE ;
    if(!isInit) { // first init
        mpq_init(ONE) ; mpq_set_ui(ONE,1,1); mpq_init(tmpS) ;
        mpq_init(tmpP); mpq_init(tmpQ) ; isInit = 1 ;
    }
    mpq_t  * even = getEven(ctx,n,lg) ;
    if(n<2) {
        mpq_init(even[n]) ; // assume first time
        mpq_set_ui(even[n],1,1);
        return ;
    } else {
        int k ;
        mpq_set_ui(tmpS,0,1) ; // tmpS = 0
        for(k=0;k<n;k++) {
            mpq_set(tmpP,even[n-k-1]) ;// tmpP = (Pa=EV(n-k-1,l-k-1))
            mpq_add(tmpP,tmpP,ctx->backEven[k]) ; // tmpP = Pa + (Pb=EV(k,k))
            mpq_mul(tmpQ,even[n-k-1],ctx->backEven[k]) ; // tmpQ = Pa * Pb
            mpq_sub(tmpP,tmpP,tmpQ) ;// tmpP = Pa + Pb - Pa * Pb
            mpq_sub(tmpP,tmpP,tmpQ); // tmpP = Pa + Pb - 2 * Pa * Pb
            if(((k) & 1) == 0) mpq_sub(tmpP,ONE,tmpP); // ? 1 - Pa - Pb + 2 * Pa * Pb
            mpq_set_ui(tmpQ,lg-k,1) ; // tmpQ = lg-k
            mpq_mul(tmpP,tmpP,tmpQ) ;  // tmpP *= lg-k
            mpq_add(tmpS,tmpS,tmpP) ; // tmpS += tmpP
        }
        mpq_set_ui(tmpQ,1,(n*(2*lg-n+1))/2); // n*(2*lg-n+1)
        mpq_init(even[n]) ; // assume first time
        mpq_mul(even[n],tmpQ,tmpS);
        return ;
    }
}

int PB597_gmpy(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    u_int32_t n, L, sizeB ;
    CTX_597 ctx ;
    L = PB597_LG ;
    sizeB = PB597_SIZEB ;
    u_int64_t pg = PGCD64(L,sizeB);
    L /= pg ;  sizeB /= pg ; // GCD to avoid big numbers
    // Back : NB+1 for even(k,k) k=0..NB ; Front : NB+1 for even(n,L-NB+n) n=0..NB
    ctx.backEven = malloc((PB597_NB+1)*sizeof(ctx.backEven[0]));
    ctx.frontEven = malloc((PB597_NB+1)*sizeof(ctx.frontEven[0]));
    // Loop on number of boats
    for(n=0 ;n<=PB597_NB;n++) {
        computePeven(&ctx,n,n) ;
        computePeven(&ctx,n,L- PB597_NB + n ) ;
    }
    mpz_t Num1000 ;    mpz_init(Num1000);
    mpz_mul_ui(Num1000 ,mpq_numref(getEven(&ctx,PB597_NB,L)[PB597_NB]),100000000000) ;
    mpz_tdiv_q (Num1000,Num1000,mpq_denref(getEven(&ctx,PB597_NB,L)[PB597_NB]));
    if(pbR->isVerbose){
        gmp_fprintf (stdout,"\t PB%dy %2d boats PEVEN=0.%Zd=%Qd\n",pbR->pbNum,PB597_NB,Num1000,getEven(&ctx,PB597_NB,L)[PB597_NB]);
    }
    gmp_sprintf (pbR->strRes,"%Zd",Num1000);
    pbR->nbClock = clock() - pbR->nbClock ;
    free(ctx.backEven) ; free(ctx.frontEven) ;
    return 1 ;
}

