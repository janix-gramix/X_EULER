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
#define PB696_NBN   9




#define PB696_MOD   1000000007

#define NB_T696 12  // 00 01 02 03 04 11 12 13 14 22 23 24
#define NB_STP696   11
static inline int nb2nty(int nb2,int nb1) {
    static int nb2nty[]={0,5,9} ;
    if(nb2>nb1) nb2 = nb1 ;
    if(nb2 > 2) nb2 =2 ;
    return nb2nty[nb2]+nb1-nb2 ;
}
// #define DBG
#if defined(DBG)
#define DBGSDIG    10
#define DBGSH2  100
#define DBGD2   11
#define DBGT3   111
#define DBGP3   (DBGD2+1)
#define DBGM3   (DBGSH2-1)
#define DBGSH3  1000
#define DBGSH5  100000
#define DBGSH6  1000000
#define DBGSH8  100000000
#define DBGFMT  "%lld->%lld\n"

#endif

//
#define isSerie 1
#define isT3    2
#define isPaire 4
static inline int nxtIp (int ip,int flag) {
    int nxtIp ;
    if(ip==0) nxtIp = 0 ;
    else if(ip==1 && (flag & isT3) )nxtIp= 8 ;
    else if(ip==8 ) {
        if(flag & isT3) {
            if(flag & isSerie) nxtIp = 9 ;
            else nxtIp = 10 ;
        } else {
            if(flag & isSerie) nxtIp = 3 ;
            else nxtIp = 6 ;
        }
    } else if(ip==9) nxtIp = 4 ;
    else if(ip==10) {
        nxtIp = 4 ;
    } else if(ip==2 && !(flag & isSerie) ) nxtIp = 6 ;
    else if(ip==5 && (flag & isSerie)) nxtIp = 3 ;
    else if(ip==5 && !(flag & isSerie)) nxtIp = 7 ;
    else nxtIp = (ip==7) ? 7 : ip+1 ;
    return nxtIp ;
}
int PB696b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB696_NBN ;
    int ntMax = PB696_NBT ;
    int ns = PB696_NBS ;
    //    int nty2nb2[NB_T696] = { 0,0,0,1,1,2} ;
    //    int nty2nb1[NB_T696] = { 0,1,2,1,2,2} ;
    int nty2nb2[NB_T696] = { 0,0,0,0,0,1,1,1,1,2,2,2} ;
    int nty2nb1[NB_T696] = { 0,1,2,3,4,1,2,3,4,2,3,4} ;
    int32_t *tb1 = malloc(NB_T696*(ntMax+1)*NB_STP696*sizeof(tb1[0]));
    int32_t *tb2 = malloc(NB_T696*(ntMax+1)*NB_STP696*sizeof(tb2[0]));
#if defined(DBG)
    int64_t *val1 = calloc(NB_T696*(ntMax+1)*NB_STP696*10000,sizeof(val1[0]));
    int64_t *val2 = calloc(NB_T696*(ntMax+1)*NB_STP696*10000,sizeof(val2[0]));
    int64_t *valCur= val1 ;
    int64_t *valAnt = val2 ;
#endif
    
    int64_t nbNoP[ntMax+1] ;
    int64_t nbWithP[ntMax+1] ;
    int nt ;
    for(nt=0;nt<=ntMax;nt++) {
        int32_t *tbCur= tb1 ;
        int32_t *tbAnt = tb2 ;
        memset(tbCur,0,NB_T696*(nt+1)*NB_STP696*sizeof(tbCur[0])) ;
        tbCur[0] = 1 ;
        
        
        for(int i=0;i<n;i++) {
#if defined(DBG)
            int64_t *val = valCur ;
            valCur = valAnt ;
            valAnt = val ;
            printf("\n************ %d **********\n",i);
#endif
            int32_t *tmp = tbCur ;
            tbCur = tbAnt ;
            tbAnt = tmp ;
            memset(tbCur,0,NB_T696*(nt+1)*NB_STP696*sizeof(tbCur[0])) ;
            for(int ip=0;ip<NB_STP696;ip++) { // paire ou pas paire
                for(int nty=0;nty<NB_T696;nty++) {
                    int nb2 = nty2nb2[nty] ;
                    int nb1 = nty2nb1[nty]   ;
                    int nb0 = 4 ;
                    for(int it=0;it<=nt;it++) {
                        int nxt, old = (it*NB_T696+nty)*NB_STP696+ip ;
                        int na = tbAnt[ old] ;
                        if(na==0) continue ;
                        //                    printf(" %d,%d,%d=%d->",it,nty,ip,na);
                        //                    printf(" %lld->",valAnt[old*10000+na-1]);
                        //                   if(nb2>=2 && it < nt-1 ) {
                        //                    if(nb2>=2 && (it < nt-1) && ip != 6 && ip != 3) {
                        if(nb2>=2 && (it < nt-1)  && (ip != 3) && (ip != 6) && (ip !=9 && (ip != 10))) {
                            // on rajoute 2xfois la serie
                            int nxty = nb2nty(nb1-2,2);
                            int nxtp =  nxtIp (ip,isSerie);
                            int nxit = it+2 ;
                            nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH6 + (DBGT3*i-DBGM3) * (DBGSH3+1) ;
                                printf("(%d,%d,%d=%d=>%d,%d,%d)" DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            if(ip==0) {
                                int nxty = 0 ;
                                int nxtp =  1 ;
                                int nxit = it+2 ;
                                // on rajoute 2xfois la serie + la paire
                                nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                                for(int iv=0;iv<na;iv++) {
                                    valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH8 + (DBGT3*i-DBGM3) * (DBGSH2*(DBGSH3+1)) + DBGD2*(i+1) ;
                                    printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                                }
#endif
                                tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                            }
                        }
                        
                        if(nb2 && it < nt){ // Serie possible
                            // on rajoute la serie seulement
                            int nxty = nb2nty(nb1-1,3) ;
                            int nxtp =  nxtIp (ip,isSerie);
                            int nxit = it+1 ;
                            nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH3 + (DBGT3*i-DBGM3) ;
                                printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                            if(ip==0) {
                                // on rajoute la serie +P0
                                int nxty = nb2nty(nb1-1,1) ;
                                int nxtp = 1 ;
                                int nxit = it+1 ;
                                nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                                for(int iv=0;iv<na;iv++) {
                                    valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH5 + (DBGT3*i-DBGM3)*DBGSH2 + DBGD2*(i+1) ;
                                    printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                                }
#endif
                                tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                            }
                            if((it < nt - 1) && ip != 3 && ip != 9 && ip != 10) {
                                // on rajoute la serie +T0 // nxty = 0
                                int nxty = 0 ;
                                int nxtp =  nxtIp (ip,isSerie|isT3);
                                int nxit = it + 2 ;
                                nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                                for(int iv=0;iv<na;iv++) {
                                    valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH6 + (DBGT3*i-DBGM3)*DBGSH3 + DBGT3*(i+1) ;
                                    printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                                }
#endif
                                tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD  ;
                            }
                        }
                        if(ip==0) {
                            // on rajoute P0
                            int nxty = nb2nty(nb1,2)  ; //
                            int nxtp = 1 ;
                            int nxit = it ;
                            nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH2 + DBGD2*(i+1) ;
                                printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        }
                        if(it < nt && ip != 3 && ip != 9 && ip !=10 ) {
                            // on rajoute T0
                            int nxty = nb2nty(nb1,1) ; //
                            int nxtp =  nxtIp (ip,isT3);
                            int nxit = it + 1 ;
                            nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH3 + DBGT3*(i+1) ;
                                printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        }
                        int nxty = nb2nty(nb1,4) ;
                        int nxtp =  nxtIp (ip,0);
                        int nxit = it ;
                        nxt = (nxit*NB_T696+nxty)*NB_STP696+nxtp ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] ;
                            printf("(%d,%d,%d=%d=>%d,%d,%d)"  DBGFMT ,it,nty,ip,na,nxit,nxty,nxtp,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                    }
                }
            }
            
        }
        int nbV=0 ;
        int noIP = 0, withIP = 0 ;
        for(int ii=0;ii<NB_T696;ii++) {
            noIP += tbCur[(nt*NB_T696+ii)*NB_STP696] ; noIP %= PB696_MOD  ;
            for(int iip=1;iip<NB_STP696;iip++) {
                int ind = (nt*NB_T696+ii)*NB_STP696+iip ;
                withIP += tbCur[ind] ; withIP %= PB696_MOD ;
            }
        }
        nbNoP[nt] = noIP ;
        nbWithP[nt] = withIP ;
        
        
        printf("nt=%d noIP=%lld , withIP=%lld\n",nt,nbNoP[nt],nbWithP[nt]);
    }
    int64_t Sum = 0 ;
    
    int64_t CnbNoP[ntMax+1] ;
    int64_t CnbWithP[ntMax+1] ;
    for(nt=0;nt<=ntMax;nt++) {
        CnbNoP[nt] = nbNoP[nt] ;
        CnbWithP[nt] = nbWithP[nt] ;
    }
    for(int is=2;is<=ns;is++) {
        int64_t OnbNoP[ntMax+1] ;
        int64_t OnbWithP[ntMax+1] ;
        memcpy(OnbNoP,CnbNoP,sizeof(OnbNoP)) ;
        memcpy(OnbWithP,CnbWithP,sizeof(OnbWithP)) ;
        memset(CnbNoP,0,sizeof(CnbNoP)) ;
        memset(CnbWithP,0,sizeof(CnbWithP)) ;
        for(nt=0;nt<=ntMax;nt++) {
            for(int nt1=0;nt1<=nt;nt1++) {
                CnbNoP[nt] += OnbNoP[nt1] * nbNoP[nt-nt1] ; CnbNoP[nt] %= PB696_MOD ;
                CnbWithP[nt] += OnbWithP[nt1] * nbNoP[nt-nt1] + OnbNoP[nt1] * nbWithP[nt-nt1] ;
                CnbWithP[nt] %= PB696_MOD ;
            }
        }
        
    }
    printf("nt=%d ns=%d n=%d Sol=%lld\n",ntMax,ns,n,CnbWithP[ntMax]);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",CnbWithP[ntMax]);
    free(tb1) ; free(tb2) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


