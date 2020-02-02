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

#define PB696_NBS   100000000
#define PB696_NBT   30
#define PB696_NBN   100000000


#define PB696_NBS   3
#define PB696_NBT   4
#define PB696_NBN   9

/*
#define PB696_NBS   1000
#define PB696_NBT   5
#define PB696_NBN   1000
*/


#define PB696_MOD   1000000007


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
#if defined(DBG)
typedef struct VAL2 {
    int64_t val ;
    int64_t vald ;
} VAL2 ;
int cmpVal2(const void * el1, const void *el2) {
    VAL2 v1 = ((VAL2 *)el1)[0] ;
    VAL2 v2 = ((VAL2 *)el2)[0] ;
    if(v1.vald>v2.vald) return 1;
    else if (v1.vald<v2.vald) return -1;
    return 0 ;
}


int cmpVal696(const void * el1, const void *el2) {
    int64_t i1 = ((int64_t *)el1)[0] ;
    int64_t i2 = ((int64_t *)el2)[0] ;
    if(i1>i2) return 1;
    else if (i1<i2) return -1;
    if(el1 > el2) return 1;
    else return -1 ;
    //    return 0 ;
}
#endif

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


#if 1
#define NB_T696a 15     // 44 33 34 22 23 24 11 12 13 14 00 01 02 03 04

#define NB_STP696a   11
static inline int nb2ntya(int nb1,int nb2) {
    static int nb2nty[]={10,6,3,1,0} ;
    //    if(nb1>nb2) nb1 = nb2 ;
    return nb2nty[nb1]+nb2-nb1 ;
}


//
static inline int nxtSta (int st,int flag) {
    int nxtSt ;
    if(st==0) {
        if(flag & isPaire ) nxtSt = (flag & isSerie) ? 1 : 6 ;
        else nxtSt = 0 ;
    } else if(st==1 && (flag & isT3) )nxtSt= 4 ;
    else if(st==6 && !(flag & isT3) ) nxtSt = 10 ;
    else if(st==7 && !(flag & isT3) ) nxtSt = 10 ;
    else if (st==4 && !(flag & isT3) ) nxtSt = 3;
    else if(st==3) nxtSt = (flag & isSerie) ? 1 : 10 ;
    else if(st==5) nxtSt = (flag & isSerie) ? 1 : 10 ;
    else if(st==8) nxtSt = (flag & isSerie) ? 1 : 10  ;
    else nxtSt = (st==10) ? 10 : st+1 ;
    return nxtSt ;
}
int PB696b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB696_NBN ;
    int ntMax = PB696_NBT ;
    int ns = PB696_NBS ;
    //    // 44 33 34 22 23 24 11 12 13 14 00 01 02 03 04
    int nty2nb1[NB_T696a] = { 4,3,3,2,2,2,1,1,1,1,0,0,0,0,0} ;
    int nty2nb2[NB_T696a] = { 4,3,4,2,3,4,1,2,3,4,0,1,2,3,4} ;
    
    int64_t *tb1 = malloc(NB_T696a*(ntMax+1)*NB_STP696a*sizeof(tb1[0]));
    int64_t *tb2 = malloc(NB_T696a*(ntMax+1)*NB_STP696a*sizeof(tb2[0]));
    
#if defined(DBG)
    int64_t *val1 = calloc(NB_T696a*(ntMax+1)*NB_STP696a*10000,sizeof(val1[0]));
    int64_t *val2 = calloc(NB_T696a*(ntMax+1)*NB_STP696a*10000,sizeof(val2[0]));
    int64_t *valCur= val1 ;
    int64_t *valAnt = val2 ;
#endif
    
    int nt ;
    //    int maxI = ntMax*4+5 ;
    int maxI = ntMax*3+1 ;
    //    int maxI=n;
    if(maxI>n)maxI = n ;
    int maxI1 = maxI+1 ;
#define INDXITIN(it,in) ((it)*maxI1+(in))
#define INDCOUNT(it,nty,st)    (((it)*NB_T696a+(nty))*NB_STP696a+(st))
    
    printf("Maxi=%d n=%d\n",maxI,n);
    int64_t *nbConnextNoP = calloc((maxI+1)*(ntMax+1),sizeof(nbConnextNoP[0]));
    int64_t *nbConnextWithP = calloc((maxI+1)*(ntMax+1),sizeof(nbConnextWithP[0]));
    nt = ntMax ;
    int64_t *tbCur= tb1 ;
    int64_t *tbAnt = tb2 ;
    memset(tbCur,0,NB_T696a*(nt+1)*NB_STP696a*sizeof(tbCur[0])) ;
    tbCur[0] = 1 ;
    for(int in=0;in<=maxI;in++) {
#if defined(DBG)
        int64_t *val = valCur ;
        valCur = valAnt ;
        valAnt = val ;
        printf("\n************ %d **********\n",in);
#endif
        int64_t *tmp = tbCur ;
        tbCur = tbAnt ;
        tbAnt = tmp ;
        memset(tbCur,0,NB_T696a*(nt+1)*NB_STP696a*sizeof(tbCur[0])) ;
        for(int st=0;st<NB_STP696a;st++) { // state to avoid double patterns
            for(int nty=0;nty<NB_T696a;nty++) {
                int nb2 = nty2nb2[nty] ;
                int nb1 = nty2nb1[nty]   ;
                int isSerieAllowed = ( in<n-2) ? 1 : 0 ;
                for(int it=0;it<=nt;it++) {
                    int nxt, old = INDCOUNT(it,nty,st);
                    int64_t na = tbAnt[ old] ;
                    if(na==0) continue ;
                    //                        if( isSerieAllowed  && nb1>=2 && (it < nt-1)  && (ip != 3) && (ip != 5) && (ip !=9 && (ip != 10))) {
                    if( isSerieAllowed  && nb1>=2 && (it < nt-1)  && (st != 1) && (st != 6)) {
                        // on rajoute 2xfois la serie
                        int nxty = nb2ntya(nb2-2,2);
                        int nxst =  nxtSta (st,isSerie);
                        int nxit = it+2 ;
                        nxt = INDCOUNT(nxit,nxty,nxst)  ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH6 + (DBGT3*(in+1)+DBGP3) * (DBGSH3+1) ;
                            printf("(%d,%d,%d=%lld=>%d,%d,%d)" DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        if(nb1>=4 && st==0) {
                            int nxty = nb2ntya(nb2-2,2) ;
                            int nxst = nxtSta(0,isPaire| isSerie) ;
                            int nxit = it+2 ;
                            // on rajoute 2xfois la serie + la paire
                            nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH8 + (DBGT3*(in+1)+DBGP3) * (DBGSH2*(DBGSH3+1)) + DBGD2*(in+1) ;
                                printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        }
                    }
                    
                    if( isSerieAllowed && nb1 && it < nt){ // Serie possible
                        // on rajoute la serie
                        int nxty = nb2ntya(nb2-1,3) ;
                        int nxst =  nxtSta (st,isSerie);
                        int nxit = it+1 ;
                        nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH3 + (DBGT3*(in+1)+DBGP3) ;
                            printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        if(nb1>=3 && st==0) {
                            // on rajoute la serie +P0
                            int nxty = nb2ntya(nb2-1,3) ;
                            int nxst = nxtSta(0,isPaire|isSerie) ;
                            int nxit = it+1 ;
                            nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH5 + (DBGT3*(in+1)+DBGP3)*DBGSH2 + DBGD2*(in+1) ;
                                printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                        }
                        if(nb1>=4 && (it < nt - 1) && st != 3 && st != 5 && st != 8) {
                            // on rajoute la serie +T0 // nxty = 0
                            int nxty =  nb2ntya(nb2-1,3) ;
                            int nxst =  nxtSta (st,isSerie|isT3);
                            int nxit = it + 2 ;
                            nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                            for(int iv=0;iv<na;iv++) {
                                valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH6 + (DBGT3*(in+1)+DBGP3)*DBGSH3 + DBGT3*(in+1) ;
                                printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                            }
#endif
                            tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD  ;
                        }
                        
                    }
                    if(st==0 && nb1>=2) {
                        // on rajoute P0
                        int nxty = nb2ntya(nb2,4)  ; //
                        int nxst = nxtSta(0,isPaire) ;
                        int nxit = it ;
                        nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH2 + DBGD2*(in+1) ;
                            printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                    }
                    if(nb1>=3 && it < nt && st != 3 && st != 5 && st !=8 ) {
                        // on rajoute T0
                        int nxty = nb2ntya(nb2,4) ; //
                        int nxst =  nxtSta (st,isT3);
                        int nxit = it + 1 ;
                        nxt = (nxit*NB_T696a+nxty)*NB_STP696a+nxst ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] * DBGSH3 + DBGT3*(in+1) ;
                            printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                    }
                    if(nb1==4) {
                        // if(it==nt)
                        {
                            if(st==0) {
                                //                                   printf("+%d(%d,%d) ",na,it,i);
                                nbConnextNoP[INDXITIN(it,in)] += na ; nbConnextNoP[INDXITIN(it,in)] %= PB696_MOD  ;
                            }else {
                                //  printf("+w%d(%d,%d) ",na,it,i);
                                nbConnextWithP[INDXITIN(it,in)] += na ; nbConnextWithP[INDXITIN(it,in)] %= PB696_MOD  ;
                            }
                        }
                    } else {
                        
                        int nxty = nb2ntya(nb2,4) ;
                        int nxst =  nxtSta (st,0);
                        int nxit = it ;
                        nxt = INDCOUNT(nxit,nxty,nxst) ;
#if defined(DBG)
                        for(int iv=0;iv<na;iv++) {
                            valCur[nxt*10000+tbCur[nxt]+iv] = valAnt[old*10000+iv] ;
                            printf("(%d,%d,%d=%lld=>%d,%d,%d)"  DBGFMT ,it,nty,st,na,nxit,nxty,nxst,valAnt[old*10000+iv],valCur[nxt*10000+tbCur[nxt]+iv]);
                        }
#endif
                        tbCur[nxt] += na ; tbCur[nxt] %= PB696_MOD ;
                    }
                }
                
            }
        }
    }
    
#if defined(DBG)
    int64_t *tbSortVal=malloc(100000*sizeof(tbSortVal[0]));
    VAL2 * tbVal2 = malloc(100000*sizeof(tbVal2[0]));
#endif
    
    int nbV=0 ;
    int noIP = 0, withIP = 0 ;
    for(int ii=0;ii<NB_T696a;ii++) {
        noIP += tbCur[(nt*NB_T696a+ii)*NB_STP696a] ; noIP %= PB696_MOD  ;
        for(int iip=1;iip<NB_STP696a;iip++) {
            int ind = (nt*NB_T696a+ii)*NB_STP696a+iip ;
            withIP += tbCur[ind] ; withIP %= PB696_MOD ;
#if defined(DBG)
            for(int iv=0;iv<tbCur[ind];iv++) tbSortVal[nbV++] = valCur[10000*ind+iv] ;
            
#endif
        }
    }
    //        nbNoP[nt] = noIP ;
    //        nbWithP[nt] = withIP ;
#if defined(DBG)
    int nbId = 0 ;
    printf("\n nbV=%d\n",nbV);
    qsort(tbSortVal,nbV,sizeof(tbSortVal[0]),cmpVal696) ;
    for(int iv=1;iv<nbV;iv++) {
        if(tbSortVal[iv-1]==tbSortVal[iv])
            printf("?%lld\n",tbSortVal[iv]);
    }
    int dig[20] ;
    int histo[DBGSDIG] ;
    for(int iv=0;iv<nbV;iv++) {
        memset(histo,0,sizeof(histo));
        int64_t val = tbSortVal[iv] ;
        for(int i=0;i<nt*3+2;i++) {
            dig[i] = val % DBGSDIG ;
            val /= DBGSDIG ;
            histo[dig[i]]++ ;
            if(histo[dig[i]]>4) printf("PBBBB %llx\n",tbSortVal[iv]);
        }
        tbVal2[iv].val = tbSortVal[iv] ;
        int64_t vald = 0 ;
        for(int id=1;id<DBGSDIG;id++) {
            for(int ih=0;ih<histo[id];ih++) vald = DBGSDIG *vald + id ;
        }
        tbVal2[iv].vald = vald ;
        int nbt=0, nbp=0 ;
        int isOK = 0;
        int ip ;
        for(ip=0;ip<nt*3+2;ip+=3) {
            if(dig[ip] != dig[ip+1]) continue ;
            int it ;
            for(it=0;it<ip;it+=3) {
                if(dig[it] == dig[it+1] && dig[it] == dig[it+2]) continue ;
                if(dig[it] == dig[it+1]+1 && dig[it] == dig[it+2]+2) continue ;
                break ;
            }
            if(it != ip) continue ;
            for(it=ip+2;it<nt*3+2;it+=3) {
                if(dig[it] == dig[it+1] && dig[it] == dig[it+2]) continue ;
                if(dig[it] == dig[it+1]+1 && dig[it] == dig[it+2]+2) continue ;
                break ;
            }
            if(it==nt*3+2) { isOK = 1; break ;}
        }
        if(!isOK) {
            printf("*%lld\n",tbSortVal[iv]);
        } else {
            /*           for(int it=0;it<ip;it+=3) printf("%d%d%d.",dig[it],dig[it+1],dig[it+2]) ;
             printf("%d%d.",dig[ip],dig[ip+1]) ;
             for(int it=ip+2;it<nt*3+2;it+=3) printf("%d%d%d.",dig[it],dig[it+1],dig[it+2]) ;
             printf("\n");
             */       }
        
    }
    qsort(tbVal2 ,nbV,sizeof(tbVal2[0]),cmpVal2) ;
    for(int iv=1;iv<nbV;iv++) {
        if(tbVal2[iv-1].vald==tbVal2[iv].vald) {
            nbId++ ;
            if(tbVal2[iv-1].val>tbVal2[iv].val) printf(DBGFMT,tbVal2[iv-1].val,tbVal2[iv].val);
            else printf(DBGFMT,tbVal2[iv].val,tbVal2[iv-1].val);
        }
    }
#endif
    
    //        printf("nt=%d noIP=%lld , withIP=%lld\n",nt,nbNoP[nt],nbWithP[nt]);
    
    int64_t *nbNoP = calloc(ntMax+1,sizeof(nbNoP[0])) ;
    int64_t *nbWithP = calloc(ntMax+1,sizeof(nbWithP[0])) ;
    
    
    
    int ntMax1= ntMax+1 ;
    //#define INDXITIN(it,in) ((in)*ntMax1+(it))
    
    int64_t *nbN_NoP = calloc((maxI+1)*(ntMax+1),sizeof(nbN_NoP[0]));
    int64_t *nbN_WithP = calloc((maxI+1)*(ntMax+1),sizeof(nbN_WithP[0]));
    
    int64_t *nbT_NoP1 = calloc((maxI+1)*(ntMax+1),sizeof(nbT_NoP1[0]));
    int64_t *nbT_WithP1 = calloc((maxI+1)*(ntMax+1),sizeof(nbT_WithP1[0]));
    int64_t *nbT_NoP2 = calloc((maxI+1)*(ntMax+1),sizeof(nbT_NoP2[0]));
    int64_t *nbT_WithP2 = calloc((maxI+1)*(ntMax+1),sizeof(nbT_WithP2[0]));
    int64_t *cur_nbT_NoP = nbT_NoP1 ;
    int64_t *ant_nbT_NoP = nbT_NoP2 ;
    int64_t *cur_nbT_WithP = nbT_WithP1 ;
    int64_t *ant_nbT_WithP = nbT_WithP2 ;
    int i ;
    
    for(i=0;i<=maxI;i++) {
        int it ;
        for(it=0;it<=ntMax;it++) {
            cur_nbT_NoP[INDXITIN(it,i)] = nbConnextNoP[INDXITIN(it,i)] ;
            cur_nbT_WithP[INDXITIN(it,i)] = nbConnextWithP[INDXITIN(it,i)] ;
        }
    }
    /*
     {
     int it ;
     for(it=0;it<=ntMax;it++) {
     printf("itC=%d: ",it) ;
     for(i=0;i<=maxI;i++) printf("%lld ",cur_nbT_NoP[i*(ntMax+1)+it]);
     printf("\n\tw: ") ;
     for(i=0;i<=maxI;i++) printf("%lld ",cur_nbT_WithP[i*(ntMax+1)+it]);
     printf("\n");
     }
     
     }
     */
    int64_t invI[ntMax+2] ;
    invI[1] =1 ;
    for(int ii=2;ii<=ntMax+1;ii++) invI[ii] = modPow(ii,PB696_MOD-2 , PB696_MOD) ;
    for(int it=0;it<=ntMax;it++) {
        for(int in=1;in<=maxI;in++) {
            nbNoP[it] += cur_nbT_NoP[INDXITIN(it,in)] * (n-in+1) ; nbNoP[it] %= PB696_MOD ;
            //            if(it==2 && cur_nbT_WithP[in*(ntMax+1)+it]) printf("C1-%lld+=(%d,%d)%lldx%lld ",nbWithP[it],in,it,cur_nbT_WithP[in*(ntMax+1)+it], n-in+1);
            nbWithP[it] += cur_nbT_WithP[INDXITIN(it,in)]*(n-in+1) ; nbWithP[it] %= PB696_MOD ;
            int64_t fact = (n-in)*(int64_t)(n-in-1) ; fact %= PB696_MOD ;
            //           if(it==2 && cur_nbT_NoP[in*(ntMax+1)+it]) printf("C1-%lld+=no(%d,%d)%lldx%lld ",nbWithP[it],in,it,cur_nbT_NoP[in*(ntMax+1)+it],fact);
            nbWithP[it] += cur_nbT_NoP[INDXITIN(it,in)] * fact ; nbWithP[it] %= PB696_MOD ;
        }
    }
    /*
     printf("Cn=%d :",1) ;
     for(int it=0;it<=ntMax;it++) printf("%lld ",nbNoP[it]);
     printf("\n\tw:");
     for(int it=0;it<=ntMax;it++) printf("%lld ",nbWithP[it]);
     printf("\n");
     */
    int ic ;
    for(ic=2;ic<=ntMax;ic++) {
        int64_t *tmp = ant_nbT_NoP ;
        ant_nbT_NoP = cur_nbT_NoP ;
        cur_nbT_NoP = tmp ;
        tmp = ant_nbT_WithP ;
        ant_nbT_WithP = cur_nbT_WithP ;
        cur_nbT_WithP = tmp ;
        int t1,t2 ;
        int i1,i2 ;
        for(t1=0;t1<=ntMax;t1++) {
            int maxi1 = 3*t1+1 ;
            if(maxi1>maxI) maxi1=maxI ;
            for(i1=1;i1<=maxi1;i1++) {
                if(i1*4 < t1*3) continue ;
                int64_t Sno = 0;
                int64_t Swith = 0 ;
                for(t2 = 1 ; t2 <t1;t2++) {
                    int maxi2 = 3*t2+2 ;
                    if(maxi2>i1)maxi2=i1 ;
                    for(i2=1;i2<maxi2;i2++) {
                        int64_t delta = ant_nbT_NoP[INDXITIN(t2,i2)] * nbConnextNoP[INDXITIN(t1-t2,i1-i2)] ;
                        if(delta){ Sno += delta ; Sno %= PB696_MOD ;  }
                        delta = ant_nbT_WithP[INDXITIN(t2,i2)] * nbConnextNoP[INDXITIN(t1-t2,i1-i2)]
                        + ant_nbT_NoP[INDXITIN(t2,i2)] * nbConnextWithP[INDXITIN(t1-t2,i1-i2)] ;
                        if(delta)  {  Swith += delta ; Swith %= PB696_MOD ;   }
                    }
                }
                //                    printf(" =>%lld\n",Sno);
                cur_nbT_NoP[INDXITIN(t1,i1)] = Sno ;
                cur_nbT_WithP[INDXITIN(t1,i1)] = Swith ;
            }
        }
        /*        {
         printf("**** ic=%d ****\n",ic);
         int it,i ;
         for(it=0;it<=ntMax;it++) {
         printf("it=%d:",it);
         for(i=0;i<=maxI;i++) printf("%lld ",cur_nbT_NoP[i*(ntMax+1)+it]);
         printf("\n   w ");
         for(i=0;i<=maxI;i++) printf("%lld ",cur_nbT_WithP[i*(ntMax+1)+it]);
         printf("\n");
         }
         }
         */
        for(int it=ic;it<=ntMax;it++) {
            int in ;
            for(in=2;in<=maxI;in++) {
                int64_t fact = n-in+1 ;
                for (int ij=2;ij<=ic;ij++) {
                    fact *=  (n-in+2-ij) ; fact %= PB696_MOD ;
                    fact *= invI[ij] ; fact %= PB696_MOD ;
                }
                nbNoP[it] += cur_nbT_NoP[INDXITIN(it,in)] * fact ; nbNoP[it] %= PB696_MOD ;
                nbWithP[it] += cur_nbT_WithP[INDXITIN(it,in)] * fact ; nbWithP[it] %= PB696_MOD ;
                fact *= (n-in+1-ic) ; fact %= PB696_MOD ;
                nbWithP[it] += cur_nbT_NoP[INDXITIN(it,in-1)] * fact ; nbWithP[it] %= PB696_MOD ;
            }
        }
    }
    
    free(nbN_NoP) ; nbN_NoP= NULL ;
    free(nbN_WithP) ; nbN_WithP= NULL ;
    free(cur_nbT_NoP); cur_nbT_NoP = NULL ;
    free(ant_nbT_NoP); ant_nbT_NoP = NULL ;
    free(cur_nbT_WithP); cur_nbT_WithP = NULL ;
    free(ant_nbT_WithP); ant_nbT_WithP = NULL ;
    
    
    nbNoP[0] = 1 ;
    printf("Cn=%d :",ntMax) ;
    for(int it=0;it<=ntMax;it++) printf("%lld ",nbNoP[it]);
    printf("\n\tw:");
    for(int it=0;it<=ntMax;it++) printf("%lld ",nbWithP[it]);
    printf("\n");
    
    printf(" End phase 1 %.6fs\n", (float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC);
    int64_t Sum = 0 ;
    
    int64_t CnbNoP[ntMax+1] ;
    int64_t CnbWithP[ntMax+1] ;
    for(nt=0;nt<=ntMax;nt++) {
        CnbNoP[nt] = nbNoP[nt] ;
        CnbWithP[nt] = nbWithP[nt] ;
    }
    
    Decomp  * DC = DecompAlloc(ntMax) ;
    int64_t ntBySeries[ntMax] ;
    int64_t ntMultiplicy[ntMax] ;
    int nbSeriesGroup  ;
    int64_t invFactorial[ntMax+2] ;
    int64_t arrangement[ntMax+2] ;
    int64_t *nbNoP_Mult = malloc(ntMax*(ntMax+1)*sizeof(nbNoP_Mult[0]));
    
    {
        int64_t fact = 1 ;
        int64_t factA = 1 ;
        invFactorial[0] = invFactorial[1] = fact ;
        arrangement[0] = factA ;
        factA = ns ;
        arrangement[1] = factA ;
        int64_t n ;
        for(n=2;n<=ntMax+1;n++) {
            fact *= invI[n] ; fact %= PB696_MOD ;
            invFactorial[n] = fact ;
            factA *= ns -n +1 ; factA %= PB696_MOD ;
            arrangement[n] = factA ;
        }
        for(int it=1;it<=ntMax;it++) {
            fact = 1 ;
            nbNoP_Mult[(it-1)*(ntMax+1)] = fact ;
            int ie , iemax = ntMax / it ;
            for(ie=1;ie<=iemax;ie++) {
                fact *= nbNoP[it] ; fact %= PB696_MOD ; fact *= invI[ie] ; fact %= PB696_MOD ;
                nbNoP_Mult[(it-1)*(ntMax+1)+ie] = fact ;
            }
        }
    }
    
    {
        int exps ;
        for(exps=0; (1<<exps) <= ns ;exps++) ;
        
        int64_t *nbNoP_exps = malloc(exps*(ntMax+1)*sizeof(nbNoP_exps[0]));
        int64_t *nbWithP_exps = malloc(exps*(ntMax+1)*sizeof(nbWithP_exps[0]));
        for(int it=0;it<=ntMax;it++) {
            nbNoP_exps[it] = nbNoP[it];
            nbWithP_exps[it] = nbWithP[it];
        }
        for(int is=1;is<exps;is++) {
            for(int t1=0;t1<=ntMax;t1++) {
                int64_t S = 0 , SP = 0 ;
                for(int t2=0;t2<=t1;t2++) {
                    int64_t delta = nbNoP_exps[(is-1)*(ntMax+1)+t2] * nbNoP_exps[(is-1)*(ntMax+1)+t1-t2] ;
                    delta %= PB696_MOD ;
                    S += delta ; S %= PB696_MOD ;
                    delta = nbWithP_exps[(is-1)*(ntMax+1)+t2] * nbNoP_exps[(is-1)*(ntMax+1)+t1-t2]
                    + nbNoP_exps[(is-1)*(ntMax+1)+t2] * nbWithP_exps[(is-1)*(ntMax+1)+t1-t2]  ;
                    SP += delta ; SP %= PB696_MOD ;
                }
                nbNoP_exps[is*(ntMax+1)+t1] = S ;
                nbWithP_exps[is*(ntMax+1)+t1] = SP ;
            }
        }
        int is,nis = 0 ;
        for(is=0;is<exps;is++) {
            if((1<<is) & ns) nis++ ;
        }
        int64_t *nbNoP_nis = malloc(nis*(ntMax+1)*sizeof(nbNoP_nis[0]));
        int64_t *nbWithP_nis = malloc(nis*(ntMax+1)*sizeof(nbWithP_nis[0]));
        for(is=0;is<exps;is++) {
            if((1<<is) & ns) break ;
        }
        memcpy(nbNoP_nis,nbNoP_exps+is*(ntMax+1),(ntMax+1)*sizeof(nbNoP_exps[0]));
        memcpy(nbWithP_nis,nbWithP_exps+is*(ntMax+1),(ntMax+1)*sizeof(nbWithP_exps[0]));
        for(nis=0,++is;is<exps;is++) {
            if((1<<is) & ns) {
                for(int t1=0;t1<=ntMax;t1++) {
                    int64_t S = 0 , SP = 0 ;
                    for(int t2=0;t2<=t1;t2++) {
                        int64_t delta = nbNoP_nis[nis*(ntMax+1)+t2] * nbNoP_exps[is*(ntMax+1)+t1-t2] ;
                        delta %= PB696_MOD ;
                        S += delta ; S %= PB696_MOD ;
                        delta = nbWithP_nis[nis*(ntMax+1)+t2] * nbNoP_exps[is*(ntMax+1)+t1-t2]
                        + nbNoP_nis[nis*(ntMax+1)+t2] * nbWithP_exps[is*(ntMax+1)+t1-t2]  ;
                        SP += delta ; SP %= PB696_MOD ;
                    }
                    nbNoP_nis[(nis+1)*(ntMax+1)+t1] = S ;
                    nbWithP_nis[(nis+1)*(ntMax+1)+t1] = SP ;
                }
                nis++ ;
            }
        }
        printf("%.6fs nt=%d ns=%d n=%d SP=%lld \n",(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC,nt,ns,n,nbWithP_nis[nis*(ntMax+1)+ntMax]);
        Sum = nbWithP_nis[nis*(ntMax+1)+ntMax] ;
    }
    printf("Sum=%lld \n",Sum);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    free(tb1) ; free(tb2) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#endif
