#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

typedef  uint32_t T_prime ;
typedef struct CTX_PRIMETABLE CTX_PRIMETABLE;
CTX_PRIMETABLE * Gen_tablePrime(T_prime maxValue) ;
const T_prime * GetTbPrime(CTX_PRIMETABLE * ctx) ;
uint32_t GetNbPrime(CTX_PRIMETABLE * ctx) ;

#define PB705_MAX 100000000
#define PB705_MOD  1000000007LL


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


int main(int argc,char**argv) {
    int32_t N=PB705_MAX;
    if(argc > 1) N =atoi(argv[1]);
    if(N > 1000000000) N= 1000000000 ;
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(N) ;
    clock_t nbClock = clock() ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int64_t D1=0,D2=0,D3=0,D4=0,D5=0,D6=0,D7=0,D8=0,D9=0 ;
    int64_t D21=0,D31=0,D41=0,D51=0,D61=0,D71=0,D81=0,D91=0 ;
    int64_t D32=0,D42=0,D52=0,D62=0,D72=0,D82=0,D92=0 ;
    int64_t D43=0,D53=0,D63=0,D73=0,D83=0,D93=0 ;
    int64_t D54=0,D64=0,D74=0,D84=0,D94=0 ;
    int64_t D65=0,D75=0,D85=0,D95=0 ;
    int64_t D76=0,D86=0,D96=0 ;
    int64_t D87=0,D97=0 ;
    int64_t D98=0 ;
    
    for(int i=0;i<nbPrime;i++) {
        int32_t p = tbPrime[i] ;
        int dig[10] ;
        int nbDig = 0 ;
        while(p >= 10) {
            int d = p % 10 ;
            p = p/10 ;
            if(d)dig[nbDig++] = d ;
        }
        if(p)dig[nbDig++] = p ;
         for(int nd=nbDig-1;nd>=0;nd--) {
            int d = dig[nd] ;
            switch(d) {
                case 1 : D1++; D21+=D2; D31+=D3; D41+=D4; D51+=D5; D61+=D6; D71+=D7; D81+=D8; D91+=D9; break ;
                case 2 : D2++; D32+=D3; D42+=D4; D52+=D5; D62+=D6; D72+=D7; D82+=D8; D92+=D9; break ;
                case 3 : D3++; D43+=D4; D53+=D5; D63+=D6; D73+=D7; D83+=D8; D93+=D9; break ;
                case 4 : D4++; D54+=D5; D64+=D6; D74+=D7; D84+=D8; D94+=D9; break ;
                case 5 : D5++; D65+=D6; D75+=D7; D85+=D8; D95+=D9; break ;
                case 6 : D6++; D76+=D7; D86+=D8; D96+=D9; break ;
                case 7 : D7++; D87+=D8; D97+=D9; break ;
                case 8 : D8++; D98+=D9; break ;
                case 9 : D9++; break ;
           }
        }
    }
    
    

    int64_t cost = 0 ;
    cost +=  (D21+D31+D51+D71)*36 + (D41+D91)*48 + (D61+D81)*54  ;
   {    int64_t D22=D2*(D2-1)/2 ;
        int64_t D23=D2*D3-D32; int64_t D24=D2*D4-D42; int64_t D25=D2*D5-D52 ;
        int64_t D26=D2*D6-D62; int64_t D27=D2*D7-D72; int64_t D28=D2*D8-D82 ;
        int64_t D29=D2*D9-D92 ;
        cost += D22*18 +(D32+D52+D72)*36 +(D23+D25+D27)*18 +D42*36 +D24*12 +D92*48 +D29*12  +D62*45 +D26*9 + D82*45 + D28*9 ;
    }
    {   int64_t D33=D3*(D3-1)/2 ;
        int64_t D34=D3*D4-D43; int64_t D35=D3*D5-D53 ; int64_t D36=D3*D6-D63 ;
        int64_t D37=D3*D7-D73; int64_t D38=D3*D8-D83 ; int64_t D39=D3*D9-D93 ;
        cost += D33*18 +(D53+D73)*36 +(D35+D37)*18 +D43*36 +D34*24 +D93*36 +D39*12  +D63*36 +D36*18+ D83*45 + D38*18 ;
    }
    {   int64_t D44=D4*(D4-1)/2 ;
        int64_t D45=D4*D5-D54 ; int64_t D46=D4*D6-D64 ; int64_t D47=D4*D7-D74 ;
        int64_t D48=D4*D8-D84 ; int64_t D49=D4*D9-D94 ;
        cost += D44*24 +(D54+D74)*36 +(D45+D47)*24 +D94*40 +D49*24  +D64*36 +D46*24 +D84*36 + D48*18 ;
    }
    {   int64_t D55=D5*(D5-1)/2 ;
        int64_t D56=D5*D6-D65; int64_t D57=D5*D7-D75 ;
        int64_t D58=D5*D8-D85; int64_t D59=D5*D9-D95 ;
        cost += D55*18 +D75*36 +D57*18 +D95*36 +D59*24  +D65*36 +D56*27 +D85*36 + D58*27 ;
    }
    {   int64_t D66=D6*(D6-1)/2;
        int64_t D67=D6*D7-D76; int64_t D68=D6*D8-D86; int64_t D69=D6*D9-D96;
        cost += D66*27 +D76*36 +D67*27 +D96*36 +D69*24  +D86*36 + D68*27;
    }
    {   int64_t D77 =D7*(D7-1)/2 ;
        int64_t D78=D7*D8-D87; int64_t D79= D7*D9-D97 ;
        cost += D77*18 +D97*36 +D79*24  +D87*36 + D78*27 ;
    }
    {   int64_t D88=D8*(D8-1)/2 ;
        int64_t D89=D8*D9-D98 ;
        cost += D88*27+D98*36+D89*30   ;
    }
    int64_t D99 = D9 * (D9-1)/2 ;
    cost += D99*24    ;
    int64_t costMod = cost % PB705_MOD ;
    int64_t nbPow2 = D2 + D3 + D5 + D7 + 2*D6 + 2*D8 ;
    int64_t nbPow3 = D4 + D9 ;
    costMod = (costMod * modPow(2,nbPow2-3,PB705_MOD)) % PB705_MOD ;
    if(nbPow3==0) costMod/=9 ;
    else if(nbPow3==1) costMod/=3 ;
    else costMod = (costMod * modPow(3,nbPow3-2,PB705_MOD)) % PB705_MOD ;
    nbClock = clock() - nbClock ;
    printf("\t%.3fs for Max=%d Nbinv=%lld (%%%lld) = %lld x 2**%lld x 3**%lld\n",
                             (float)nbClock/CLOCKS_PER_SEC,
                              PB705_MAX,costMod,PB705_MOD,cost,nbPow2-3,nbPow3-2);
    return 1 ;
}
