//
//  pb696
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PB696_NBS   100000000
#define PB696_NBT   30
#define PB696_NBN   100000000


#define PB696_MOD   1000000007

#define isSerie 1
#define isT3    2
#define isPaire 4
//free case for the 2 next column
#define NB_T696a 9 // 44 33 34 22 23 24 12 13 02
// conversion nb1 x nb2 => nty (compact index)
static inline int nb2ntya(int nb1,int nb2) {
    switch(nb1){
        case 4: return 0 ;
        case 3: return nb2 - 2 ;
        case 2: return nb2 + 1 ;
        case 1: return nb2 + 4 ;
        case 0:
        default: return 8 ;
    }
}
// fast exponention
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

#define NB_STP696a   7
// states (7) and automate to avoid multiple decompositions
enum PB696ST { ST_Start = 0 ,   ST_C2=1 ,  ST_no2S = 2, ST_1T3=3 , ST_2T3=4, ST_NoT3 = 5 , ST_End=6} ;
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
// convolution of Res = V1 * V2 
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
#define INDCOUNT(it,nty,st)    (((it)*NB_T696a+(nty))*NB_STP696a+(st))
#define SIZEBREAK  ((ntMax1)*NB_T696a*NB_STP696a)
#define SIZE_ITxIN  (maxI1*ntMax1)
#define NXT(nbSerie,nbIt,flags) int nxty = nb2ntya(nb2-(nbSerie),4-(nbSerie)) ; int nxst = nxtSta(st,(flags)) ; int nxit = it+(nbIt); int nxt = (ibCur+INDCOUNT(nxit,nxty,nxst))
#define SET_NXT tbCur[nxt] = (int32_t)((na + tbCur[nxt]) % PB696_MOD)
#define BRKNXT(newSt) (ibCur+deltaBreak+INDCOUNT(nxit,nxty,(newSt)))


int main(int argc, const char * argv[]) {
    clock_t nbClock = clock() ;
    int n = PB696_NBN ;
    int ntMax1 = PB696_NBT+1 ;
    int ns = PB696_NBS ;
    if(argc > 1){
        ntMax1 =atoi(argv[1])+1 ;
        if(ntMax1 > 2001) ntMax1=2001 ;
        if(argc > 2) {
            n = atoi(argv[2]);
            if(argc > 3)ns = atoi(argv[3]);
        }
    } else {
        printf("Args : [t (%d)] [n (%d) [s (%d)]]\n",ntMax1-1,n,ns);
    }
    int maxI = ntMax1*3-2 ; // max length for connexe pattern
    if(maxI>n)maxI = n ;
    int maxI1 = maxI+1 ;
 
    // number of free cases in the 2 next columns
    // 44 33 34 22 23 24 12 13 02
    int nty2nb1[NB_T696a] = { 4,3,3,2,2,2,1,1,0} ;
    int nty2nb2[NB_T696a] = { 4,3,4,2,3,4,2,3,2} ;

    int32_t *tb1 = malloc(ntMax1*SIZEBREAK*sizeof(tb1[0]));
    int32_t *tb2 = malloc(ntMax1*SIZEBREAK*sizeof(tb2[0]));

    // precompute combination Cnp(ib,in) =
    //  ( n-in+1 )
    //  (   ib   )
#define INDXIBIN(ib,in) ((ib)*maxI1+(in))
    NB696 *nb_IT = calloc(ntMax1,sizeof(nb_IT[0])) ;
    int64_t invI[ntMax1+1] ;
    invI[1] =1 ;
    for(int ii=2;ii<=ntMax1;ii++) invI[ii] = modPow(ii,PB696_MOD-2 , PB696_MOD) ;
    int64_t *Cnp = malloc(ntMax1*maxI1*sizeof(Cnp[0]));
    for(int in=1;in<=maxI;in++) {
        int64_t fact = 1 ;
        for(int ib=0;ib<ntMax1;ib++) {
            fact *= n-in+1-ib ; fact %= PB696_MOD ;  fact *= invI[ib+1] ; fact %= PB696_MOD ;
            Cnp[INDXIBIN(ib,in)] = fact ;
        }
    }

    int deltaBreak = SIZEBREAK;
    int32_t *tbCur= tb1 ;
    int32_t *tbAnt = tb2 ;
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
        int32_t *tmp = tbCur ; tbCur = tbAnt ;  tbAnt = tmp ;
        for(int ib=0;ib<maxIb;ib++){ // loop on breaks (frontiers between connected dev.)
             int64_t factCnp = Cnp[INDXIBIN(ib,in)] ;
             int ibCur = ib * deltaBreak ;
              for(int nty=0;nty<NB_T696a;nty++) {
                 int nb2 = nty2nb2[nty] , nb1 = nty2nb1[nty]   ;
                 for(int st=0;st<NB_STP696a;st++){ // states to avoid with multiple representation
                     for(int it= ib ;it<ntMax1;it++) {
                        int old = ibCur + INDCOUNT(it,nty,st);
                        int64_t na = tbAnt[ old] ;
                        if(na==0) continue ;
                        else tbAnt[old] = 0 ;
                        int isBreak = 0 ;
                        if(in && nb1==4){ // break possible (if not in==0 and no remaining serial<=> nb1==4)
                            if(it < ntMax1 ){
                                isBreak = 1 ;
                                if(st==ST_Start ) {
                                    nb_IT[it].noP += na * factCnp ; nb_IT[it].noP %= PB696_MOD ;
                                }  else  {
                                    nb_IT[it].withP += na * factCnp ; nb_IT[it].withP %= PB696_MOD ;
                                }
                            } else {   continue ;  } // no prolongation maximum tiles
                        }
                        if( isSerieAllowed  && nb1>=2 && (it < ntMax1-2)  && ( ((st !=ST_no2S) && (st != ST_1T3) ) || isBreak)) {
                            // add 2 serials
                            NXT(2,2,isSerie) ;
                            if( (st != ST_no2S) && (st != ST_1T3)  ) { SET_NXT ;  }
                            if(isBreak) {  nxt = (st) ? BRKNXT(ST_End) : nxt+deltaBreak ; SET_NXT ;   }
                            if(nb1>=4 && (st==ST_Start)) { // add 2 serials + Pair
                                NXT(2,2,isPaire|isSerie) ; SET_NXT ;
                                if(isBreak) { nxt += deltaBreak ; SET_NXT ; }
                            }
                        }
                        if( isSerieAllowed && nb1 && it < ntMax1-1){ // add 1xSerial ?
                            NXT(1,1,isSerie) ; SET_NXT ;
                            if(isBreak) { nxt = (st) ? BRKNXT(ST_End) : nxt+deltaBreak  ; SET_NXT ; }
                            if(nb1>=3 && (st==ST_Start) ) { // add 1xSerial + Pair
                                NXT(1,1,isSerie|isPaire) ; SET_NXT ;
                                if(isBreak) {  nxt += deltaBreak ; SET_NXT ; }
                            }
                            if(nb1>=4 && (it < ntMax1 - 2) && ((st != ST_NoT3 ) || isBreak) ) { // add T3 ?
                                NXT(1,2,isSerie|isT3) ;
                                if(st != ST_NoT3 ) { SET_NXT ;   }
                                if(isBreak) { nxt = (st) ? BRKNXT(ST_End) : nxt+deltaBreak ; SET_NXT ; }
                            }
                        }
                        if((st==ST_Start) && nb1>=2) { // add Pair ?
                            NXT(0,0,isPaire) ;   SET_NXT ;
                            if(isBreak) {  nxt += deltaBreak ; SET_NXT ; }
                        }
                        if(nb1>=3 && it < ntMax1-1 && ((st != ST_NoT3 ) || isBreak) ) { // add T3 ?
                           NXT(0,1,isT3) ;
                            if(st != ST_NoT3 ){ SET_NXT ;  }
                            if(isBreak) { nxt = (st) ? BRKNXT(ST_End) : nxt+deltaBreak ; SET_NXT ; }
                        }
                        if(nb1!=4) { // add nothing
                            NXT(0,0,0) ; SET_NXT ;
                        }
                    } // end bcl it
                } //end bcl st
            } // end bcl nty
        } // end bcl ib
    } // end bcl in
    free(tbCur); tbCur = NULL ;
    free(tbAnt); tbAnt = NULL ;
    nb_IT[0].noP = 1 ;
    //    printf("C_Cn=%d :",ntMax1-1) ; for(int it=0;it<ntMax1;it++) printf("%lld ",nb_IT[it].noP);
    //    printf("\n\tw:"); for(int it=0;it<ntMax1;it++) printf("%lld ",nb_IT[it].withP); printf("\n");
    

    printf("\t %.6fs End phase 1 : hands for one suit\n",(float)(clock()-nbClock) / CLOCKS_PER_SEC);
    int64_t Sum = 0 ;
    
    {   // compute number of suits
        int exp2 ;
        for(exp2=0; (1<<exp2) <= ns ;exp2++) ;
        
        NB696 *nb_pow2 = malloc(exp2*(ntMax1)*sizeof(nb_pow2[0]));
        for(int it=0;it<ntMax1;it++) { // init 2**0
            nb_pow2[it].noP = nb_IT[it].noP;
            nb_pow2[it].withP = nb_IT[it].withP;
        }
        for(int is=1;is<exp2;is++) { // compute 2**exp2 suits
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
        Sum = nb_cur[ntMax1-1].withP ;
    }
    printf("\t %.6fs t=%d s=%d n=%d Hands=%lld(mod %d)\n"
                               ,(float)(clock()-nbClock) / CLOCKS_PER_SEC,
                               ntMax1-1,ns,n,Sum,PB696_MOD);
    return 0 ;
}

