
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define M_PI 3.14159265358979323846264338327950288

#include "euler_utils.h"
#include "DC3.h"
#include "PB_other.h"
#include "ukkonen.h"

#define PB696_NBS   3
#define PB696_NBT   4
#define PB696_NBN   30



//#define PB696_NBS   67108864
/*
 #define PB696_NBS   1000
 #define PB696_NBT   5
 #define PB696_NBN   1000
 */


#define PB696_NBS   100000000
#define PB696_NBT   400
#define PB696_NBN   100000000


#define PB696_MOD   1000000007

#define isSerie 1
#define isT3    2
#define isPaire 4

#define NB_T696a 15     // 44 33 34 22 23 24 11 12 13 14 00 01 02 03 04

#define NB_STP696a   8
static inline int nb2ntya(int nb1,int nb2) {
    static int nb2nty[]={10,6,3,1,0} ;
    //    if(nb1>nb2) nb1 = nb2 ;
    return nb2nty[nb1]+nb2-nb1 ;
}

static u_int64_t modPow(int64_t a,int64_t exp,u_int64_t mod) {
    u_int64_t aPow2 = a % mod ;
    u_int64_t aPowExp = (exp & 1) ? aPow2 : 1 ;
    int i ;
    for(i=1;exp >= (1LL<<i) ;i++) {
        aPow2 = (aPow2 * aPow2) % mod ;
        if( (1LL<<i) & exp) {
            aPowExp = (aPowExp * aPow2 ) % mod  ;
        }
    }
    return aPowExp ;
}


enum PB696ST { ST_Start = 0 ,   ST_C2=2 ,  ST_no2S = 3, ST_1T3=4 , ST_2T3=5, ST_NoT3 = 6 , ST_End=7} ;
static inline enum PB696ST nxtSta (enum PB696ST st,int flag) {
    enum PB696ST nxtSt ;
    switch(st) {
        case ST_Start :
            if(flag & isPaire ) nxtSt = (flag & isSerie) ? ST_no2S : ST_1T3 ;
            else nxtSt = ST_Start ;
            break ;
        case ST_no2S: nxtSt =  ST_C2 ; break ;
        case ST_C2:   nxtSt = ST_NoT3; break ;
        case ST_NoT3: nxtSt = (flag & isSerie) ? ST_no2S : ST_End ; break ;
        case ST_1T3: nxtSt = (flag & isT3) ? ST_2T3 : ST_End ; break ;
        case ST_2T3: nxtSt = (flag & isT3) ? ST_NoT3 : ST_End ; break ;
        case ST_End: nxtSt = ST_End; break ;
    }
    return nxtSt ;
}


typedef struct NB696 {
    int64_t noP ;
    int64_t withP ;
} NB696 ;

static void Conv(int ntMax1,NB696 *V1 , NB696 *V2 ,NB696 *Res ) {
    for(int t1=0;t1<ntMax1;t1++) {
        int64_t S = 0 , SP = 0 ;
        for(int t2=0;t2<=t1;t2++) {
            int64_t delta = V1[t2].noP * V2[t1-t2].noP ;
            delta %= PB696_MOD ;
            S += delta ; S %= PB696_MOD ;
            delta = V1[t2].withP * V2[t1-t2].noP + V1[t2].noP * V2[t1-t2].withP  ;
            SP += delta ; SP %= PB696_MOD ;
        }
        Res[t1].noP = S ;   Res[t1].withP = SP ;
    }
    return ;
}

#define INDXITIN(it,in) (((it)-ib)*maxI1+(in))
#define INDCOUNT(it,nty,st)    (((it)*NB_T696a+(nty))*NB_STP696a+(st))
#define SIZEBREAK  ((ntMax1)*NB_T696a*NB_STP696a)
#define SIZE_ITxIN  (maxI1*ntMax1)
#define NXT(nbSerie,flags,nbIt) (ibCur+INDCOUNT(it+(nbIt),nb2ntya(nb2-(nbSerie),4-(nbSerie)),nxtSta(st,(flags))))
#define SET_NXT tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD
#define UPNXT(nxt,newSt)    (deltaBreak+((nxt) & 0xfffffff8) +(newSt) )

int PB696b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB696_NBN ;
    int ntMax1 = PB696_NBT+1 ;
    int ns = PB696_NBS ;
    int maxI = ntMax1*3-2 ; // max length for connexe pattern
    if(maxI>n)maxI = n ;
    int maxI1 = maxI+1 ;
    // number of free cases in the 2 next columns
    // 44 33 34 22 23 24 11 12 13 14 00 01 02 03 04
    int nty2nb1[NB_T696a] = { 4,3,3,2,2,2,1,1,1,1,0,0,0,0,0} ;
    int nty2nb2[NB_T696a] = { 4,3,4,2,3,4,1,2,3,4,0,1,2,3,4} ;
    
    int64_t *tb1 = malloc(ntMax1*SIZEBREAK*sizeof(tb1[0]));
    int64_t *tb2 = malloc(ntMax1*SIZEBREAK*sizeof(tb2[0]));
    NB696 ** pt_ITxIN = malloc(ntMax1*sizeof(pt_ITxIN[0])) ;
    for(int ib=0;ib<ntMax1;ib++)pt_ITxIN[ib] = calloc((maxI1)*(ntMax1-ib),sizeof(pt_ITxIN[ib][0])) ;
    int deltaBreak = SIZEBREAK;
    int64_t *tbCur= tb1 ;
    int64_t *tbAnt = tb2 ;
    memset(tbCur,0,ntMax1*SIZEBREAK*sizeof(tbCur[0])) ;
    tbCur[0] = 1 ;
// Compute the nb of compact hands (no missing tile numbers) by
//   - IT (nb tiles)
//   - IN (nb of tile numbers)
//   - IB (nb of connected parts)
    for(int in=0;in<=maxI;in++) { // loop on numbers
        int maxIb = ntMax1 ;
        if(maxIb > in+1) maxIb = in+1 ;
        int isSerieAllowed = ( in<n-2) ? 1 : 0 ;
        int64_t *tmp = tbCur ; tbCur = tbAnt ;  tbAnt = tmp ;
          for(int ib=0;ib<maxIb;ib++){ // loop on breaks (frontiers between connected dev.)
            NB696 *nb_ITxIN = pt_ITxIN[ib] ;
            int ibCur = ib * deltaBreak ;
            for(int nty=0;nty<NB_T696a;nty++) {
                int nb2 = nty2nb2[nty] ;
                int nb1 = nty2nb1[nty]   ;
                for(int st=0;st<NB_STP696a;st++){ // states to avoid with multiple representation
                    for(int it=ib;it<ntMax1;it++) {
                        int old = ibCur + INDCOUNT(it,nty,st);
                        int64_t na = tbAnt[ old] ;
                        if(na==0) continue ;
                        else tbAnt[old] = 0 ;
                        int isBreak = 0 ;
                        if(in && nb1==4){ // break possible (if not in==0 and no remaining serial<=> nb1==4)
                            if(it < ntMax1 ){
                                isBreak = 1 ;
                                if(st==ST_Start ) {
                                    nb_ITxIN[INDXITIN(it,in)].noP += na ; nb_ITxIN[INDXITIN(it,in)].noP %= PB696_MOD  ;
                                }  else  {
                                    nb_ITxIN[INDXITIN(it,in)].withP += na; nb_ITxIN[INDXITIN(it,in)].withP %= PB696_MOD   ;
                                }
                            } else {   continue ; }
                        }
                        if( isSerieAllowed  && nb1>=2 && (it < ntMax1-2)  && ( ((st !=ST_no2S) && (st != ST_1T3) ) || isBreak)) {
                            // add 2 serials
                            int nxt = NXT(2,isSerie,2) ;
                            if( (st != ST_no2S) && (st != ST_1T3)  ) { SET_NXT ; }
                            if(isBreak) {  nxt = (st) ? UPNXT(nxt,ST_End) : nxt+deltaBreak ; SET_NXT ; }
                            if(nb1>=4 && (st==ST_Start)) { // add 2 serials + Pair
                                int nxt = NXT(2,isPaire|isSerie,2) ; SET_NXT ;
                                if(isBreak) { nxt += deltaBreak ; SET_NXT ; }
                            }
                        }
                        if( isSerieAllowed && nb1 && it < ntMax1-1){ // add 1xSerial ?
                            int nxt = NXT(1,isSerie,1) ; SET_NXT ;
                            if(isBreak) { nxt = (st) ? UPNXT(nxt,ST_End) : nxt+deltaBreak  ; SET_NXT ; }
                            if(nb1>=3 && (st==ST_Start) ) { // add 1xSerial + Pair
                                int nxt = NXT(1,isSerie|isPaire,1) ; SET_NXT ;
                                if(isBreak) {  nxt += deltaBreak ; SET_NXT ; }
                            }
                            if(nb1>=4 && (it < ntMax1 - 2) && ((st != ST_NoT3 ) || isBreak) ) { // add T3 ?
                                int nxt = NXT(1,isSerie|isT3,2) ;
                                if(st != ST_NoT3 ) { SET_NXT ; }
                                if(isBreak) { nxt = (st) ? UPNXT(nxt,ST_End) : nxt+deltaBreak ; SET_NXT ; }
                            }
                        }
                        if((st==ST_Start) && nb1>=2) { // add Pair ?
                            int nxt = NXT(0,isPaire,0) ;   SET_NXT ;
                            if(isBreak) {  nxt += deltaBreak ; SET_NXT ; }
                        }
                        if(nb1>=3 && it < ntMax1-1 && ((st != ST_NoT3 ) || isBreak) ) { // add T3 ?
                            int nxt = NXT(0,isT3,1) ;
                            if(st != ST_NoT3 ){ SET_NXT ;}
                            if(isBreak) { nxt = (st) ? UPNXT(nxt,ST_End) : nxt+deltaBreak ; SET_NXT ; }
                        }
                        if(nb1!=4) { // add nothing
                            int nxt = NXT(0,0,0) ; SET_NXT ;
                        }
                    }
                 } // end bcl it
            } // end bcl nty
        } // end bcl ib
    } // end bcl in
    free(tbCur); tbCur = NULL ;
    free(tbAnt); tbAnt = NULL ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s End phase 0 : compact hands %.6fs\n",pbR->ident,(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC);
    
    
#define INDXIBIN(ib,in) ((ib)*maxI1+(in))
/*    {
        for(int ib=0;ib<ntMax1;ib++) {
            printf("**** ib=%d ****\n",ib) ;
            for(int it=ib;it<ntMax1;it++) {
                printf("it=%d ",it) ;
                for(int in=0;in<=maxI;in++) printf("%lld ",pt_ITxIN[ib][INDXITIN(it,in)].noP);
                printf("\n   w ");
                for(int in=0;in<=maxI;in++) printf("%lld ",pt_ITxIN[ib][INDXITIN(it,in)].withP);
                printf("\n");
           }
        }
    }
*/
// precompute combination Cnp(ib,in) =
//  ( n-in+1 )
//  (   ib   )

    NB696 *nb_IT = calloc(ntMax1,sizeof(nb_IT[0])) ;
    int64_t invI[ntMax1+2] ;
    invI[1] =1 ;
    for(int ii=2;ii<=ntMax1+1;ii++) invI[ii] = modPow(ii,PB696_MOD-2 , PB696_MOD) ;
    int64_t *Cnp = malloc(ntMax1*maxI1*sizeof(Cnp[0]));
    for(int in=1;in<=maxI;in++) {
        int64_t fact = 1 ;
        for(int ib=0;ib<ntMax1;ib++) {
             fact *= n-in+1-ib ; fact %= PB696_MOD ;  fact *= invI[ib+1] ; fact %= PB696_MOD ;
            Cnp[INDXIBIN(ib,in)] = fact ;
        }
    }
// From the number of hands (multiconnected, and compact(no missing number)
// by IB (nb of connected parts), by IT(nb tiles) and IN(nb used numbers) computes:
// the number of hands of length n
// by inserting missing numbers <=> multiply by Cnp
    for(int ib=0;ib<ntMax1;ib++) {
        NB696 * nb_ITxIN = pt_ITxIN[ib] ;
        for(int in=ib;in<=maxI;in++) {
            int64_t fact = Cnp[INDXIBIN(ib,in)] ;
            for(int it= ib;it<ntMax1;it++) {
                int ind = INDXITIN(it,in) ;
                if(nb_ITxIN[ind].noP) {
                    nb_IT[it].noP += nb_ITxIN[ind].noP * fact ; nb_IT[it].noP %= PB696_MOD ;
                 }
                nb_IT[it].withP += nb_ITxIN[ind].withP * fact ; nb_IT[it].withP %= PB696_MOD ;
            }
        }
    }
    
    nb_IT[0].noP = 1 ;
//    printf("C_Cn=%d :",ntMax1-1) ; for(int it=0;it<ntMax1;it++) printf("%lld ",nb_IT[it].noP);
//    printf("\n\tw:"); for(int it=0;it<ntMax1;it++) printf("%lld ",nb_IT[it].withP); printf("\n");
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s End phase 1 : hands for one suit %.6fs\n",pbR->ident,(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC);
    int64_t Sum = 0 ;
    
    {
        int exp2 ;
        for(exp2=0; (1<<exp2) <= ns ;exp2++) ;
        
        NB696 *nb_pow2 = malloc(exp2*(ntMax1)*sizeof(nb_pow2[0]));
        for(int it=0;it<ntMax1;it++) {
            nb_pow2[it].noP = nb_IT[it].noP;
            nb_pow2[it].withP = nb_IT[it].withP;
        }
        for(int is=1;is<exp2;is++) {
             Conv(ntMax1,nb_pow2+(is-1)*ntMax1 , nb_pow2+(is-1)*ntMax1 , nb_pow2+is*ntMax1 );
        }
        int is ;
        NB696 *nb_cur = malloc(ntMax1*sizeof(nb_cur[0]));
        NB696 *nb_ant = malloc(ntMax1*sizeof(nb_ant[0]));
        for(is=0;is<exp2;is++) {
            if((1<<is) & ns) break ;
        }
        nb_cur = nb_pow2+is*ntMax1 ;
        for(++is;is<exp2;is++) {
            if((1<<is) & ns) {
                NB696 * tmp = nb_ant ; nb_ant = nb_cur; nb_cur = tmp ;
                Conv(ntMax1, nb_ant, nb_pow2+is*ntMax1 ,nb_cur );
            }
        }
        printf("%.6fs SP=%lld \n",(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC,nb_cur[ntMax1-1].withP);
        Sum = nb_cur[ntMax1-1].withP ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s %.6fs t=%d s=%d n=%d Hands=%lld(mod %d)\n"
                               ,pbR->ident,(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC,
                               ntMax1-1,ns,n,Sum,PB696_MOD);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

