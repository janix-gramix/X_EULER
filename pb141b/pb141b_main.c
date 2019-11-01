#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PB_MAX_STRLEN   100

typedef struct PB_RESULT {
    const char    *ident ;
    int     isVerbose ;
    char    strRes[PB_MAX_STRLEN] ;
    clock_t nbClock ;
} PB_RESULT ;

int PB141b_e(PB_RESULT *pbR, int exp) ;

int main(int argc, char**argv) {
    PB_RESULT pbR ;
    pbR.ident = "PB141b" ;
    pbR.isVerbose = 1 ;
    int exmax = 12 ;
    if(argc > 1) {
        exmax = atoi(argv[1]);
    }
    if(exmax > 26) exmax = 26 ;
    printf("Execution for n <= 10**%d\n",exmax);
    PB141b_e(&pbR,exmax);
    
    fprintf(stdout,"\tPB%s(%.06fs) Sol=%s \n",pbR.ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes) ;
    return 0 ;
    
}


#define P141_INT128
#if defined(P141_INT128)
typedef __int128  bigInt141 ;
#else
typedef  uint64_t bigInt141 ;
#endif


typedef struct SOL141 {
    bigInt141   n; // n = r+q*d with r=k*b**2 , q = k*a*b , d =k*b**2
    int32_t     a ; // a/b goemetric ratio with a^b=1
    int32_t     b ;
    int64_t     k ;
} SOL141 ;
int AddSol141(int nbSol,SOL141 *sols,bigInt141 n, int32_t a, int32_t b, int64_t k) {
    sols[nbSol].n = n ; sols[nbSol].a = a ; sols[nbSol].b = b ; sols[nbSol].k = k ;
    //    printf("%d/%dx%lld\t%llu \n",a,b,k,(int64_t)n);
    return ++nbSol ;
}
// goodies to sort solutions
int CmpSol(const void *el1,const void *el2) {
    SOL141 * sol1 = (SOL141 *)el1 ;
    SOL141 * sol2 = (SOL141 *)el2 ;
    // key : n
    if(sol1->n > sol2->n) return 1;
    else if(sol1->n < sol2->n) return -1;
    int diff ;
    // key : a/b, k
    //    diff = sol1->a*sol2->b - sol1->b*sol2->a ;
    //    if(diff) return diff ;
    //    return (sol1->k - sol2->k) ;
    // keys : a,b,k
    diff = sol1->a - sol2->a ;
    if(diff) return diff ;
    diff = sol1->b - sol2->b ;
    if(diff) return diff ;
    if(sol1->k > sol2->k) return 1;
    else if (sol1->k < sol2->k) return -1;
    else return 0 ;
}
void print_bigInt141(bigInt141 val,char str[42]) {
    uint64_t hval = val / 1000000000000000000ULL ;
    if(hval) {
        sprintf(str,"%llu%0.18llu",hval,(uint64_t)(val-hval*1000000000000000000ULL));
    } else {
        sprintf(str,"%llu",(uint64_t)val) ;
    }
    return ;
}

#define PB141_MAX_ASK   1000000000000LL
#define EXP_PB141_MAX   26

uint64_t Sqrt64(uint64_t val) {
    return sqrtl(val);
}

uint32_t PGCD(uint32_t n1,uint32_t n2 ) {
    uint32_t r ;
    if (n1 > n2) {
        r = n2 ;
        n2 = n1 ;
    } else {
        r = n1 ;
    }
    while ( r > 0) {
        n1 = n2 ;
        n2 = r ;
        r = n1 % n2 ;
    }
    return n2 ;
}


// couple d=divisor of sf (square free number)
// coef = d * (sf/d)**2 so coef has the same divisors as sf
// and a subpart d square free
typedef struct DIV {
    int32_t coef ;
    int32_t d ;
} DIV ;

typedef struct SF {
    int32_t sf ;        // square free number
    int32_t i0 ;        // first indice in div
    int32_t inext ;     // last+1 indice in div
} SF ;

// precompute square free number sf
// for each sf decompose in sf=dfxds (w
int GetSF(SF *sf, DIV *dv, int nbMax) {
    int i,isf,id,p;
    int prime[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,0}  ;
    for(i=0;i<=nbMax;i++) sf[i].sf= i ;
    for(i=0;p=prime[i],p*p <= nbMax;i++) {
        int np2,p2=p*p ; // invalidate multiples of square
        for(np2=p2;np2<=nbMax;np2 += p2) sf[np2].sf=0 ;
    }
    isf=0 ;
    sf[isf].sf = 1 ;   sf[isf].i0 = 0 ; sf[isf].inext=1; isf++ ;
    id = 0 ;  dv[id].coef=dv[id].d=1 ; id++ ;
    for(i=2;i<=nbMax;i++) {
        if(sf[i].sf) { // loop on squarefree number sf
            int d , j ,nbd = 0 ;
            sf[isf] = sf[i] ;
            int s_f = sf[i].sf;
            sf[isf].i0 = id ;
            for(d=1;d*d<s_f;d++) {// loop on divisor of sf
                if((s_f % d)== 0) {
                    // store d x s_f/d ; d < s_f/d
                    dv[id].coef = d ; dv[id++].d = s_f/d ;   nbd++ ;
                }
            }
            for(j=0;j<nbd;j++) { // duplicate couples for d > s_f/d
                // and replace coef by coef**2 * d So foreach divisor d: (d,coef) =  (d, (sf/d)**2 * d )
                dv[id].coef =  dv[id-2*j-1].d*dv[id-2*j-1].d*dv[id-2*j-1].coef ;
                dv[id].d =  dv[id-2*j-1].coef ;
                dv[id-2*j-1].coef *= dv[id-2*j-1].coef * dv[id-2*j-1].d ;
                id++ ;
            }
            sf[isf].inext = sf[isf].i0 + 2*nbd ;
            isf++ ;
        }
    }
    sf[isf].sf = 0;  sf[isf].i0 = id ; sf[isf++].inext = id ; // terminator
    return isf-1 ;
}

// the test is base on n = m *m with m=IntPart(sqab)+1
// to assure the unicity must check a^b == 1
// as PGCD is longer , second test
int TestSol141(double sqab, bigInt141 n, int32_t a, int32_t b) {
    uint64_t m = (uint64_t)sqab + 1 ;
    return (n==m*(bigInt141)m && PGCD(a,b)==1) ;
}


int IsSquare(int64_t n) {
    static u_int8_t isSq[256] ;
    if(n < 0) {
        int i ; for(i=0;i<256;i++) isSq[i]=0 ;
        for(i=0;i<256;i++)isSq[0xff & (i*i)] = 1 ;
        return 0 ;
    }    else if(isSq[n & 0xff] == 0) {
        return 0 ;
    }
    int64_t m = Sqrt64(n) ;
    if( n==m*m ) {
        return (int) m ;
    } else {
        return 0 ;
    }
}


int PB141b_e(PB_RESULT *pbR,int ex141_max) {
    pbR->nbClock = clock() ;
    uint64_t nbTry = 0 ;
    SOL141  sol[400] ;
#if defined(P141_INT128)
    if(ex141_max > 26) ex141_max = 26 ;
#else
    if(ex141_max > 18) ex141_max = 18 ;
#endif
    
    bigInt141 pb141_max = 1;
    {  int i;  for(i=0;i<ex141_max;i++) pb141_max *= 10 ; }
    int nbSol = 0 ;
    int32_t maxB0 = (int) pow(pb141_max,1/6.0) + 1 ;
    SF *Sf = malloc((maxB0+1)*sizeof(Sf[0]));
    DIV *Div = malloc((maxB0+1)*30 *sizeof(Sf[0]));
    int nbSf = GetSF(Sf,Div,maxB0) ;
    IsSquare(-1);
    int32_t jsf ;
    SF  sf ; // current squarefree to test
    int32_t sfBK ;
    for(jsf=0;sf=Sf[jsf],sfBK = sf.sf, jsf < nbSf ;jsf++) { // loop on LCM(bf,kf) squarefree parts of b, k
        int32_t b0,k0 ;
        int32_t jb ;
        DIV divb ;
        for(jb=sf.i0;divb=Div[jb],b0=divb.coef, jb<sf.inext  ;jb++) {
            int32_t bf = divb.d ;
            uint64_t k02 ;
            DIV divk ;
            int32_t jk ; //k0 = coef = (sf/d)**2 * d ; b0s = d  (square free)
            for(jk=sf.i0 ;divk=Div[jk],k0=divk.coef , k02=k0*(uint64_t)k0, jk<sf.inext ;jk++) {
                if((divb.d * divk.d) % sfBK != 0) continue ; // PPCM(b0squareFree, k0squarefree) = sfBK
                int32_t kf = divk.d ;
                bigInt141 a3max = pb141_max / (b0*(bigInt141) k02) ;
                bigInt141 a3 ;
                int32_t a ;
                for (a = b0+1;a3 = (bigInt141)a*a*a, a3 < a3max; a++){
                    if( ( a|sfBK)&1 ) { // as sfBK is a divisor of b , must
                        long double sq_bfa = a*sqrtl(a*bf) ; // solution without additionnal square
                        uint32_t ks;
                        uint64_t ks2 ; // additional square to k => k=k*ks*ks
                        bigInt141 ks4 ;
                        for (ks=1;ks2=ks*(uint64_t)ks, ks4=ks2*(bigInt141) ks2, a3*ks4< a3max ; ks++){
                            uint64_t b12 = (ks*sq_bfa+1) ;
                            //                          uint64_t b12 = Sqrt64(ks*a3*bf);
                            uint64_t delta = b12*(bigInt141) b12 - bf*ks2*(bigInt141)a3 ;
                            nbTry++ ;
                            if( ( b0 * delta < a * (int64_t)kf ) && (delta % kf)  == 0  ) {
                                delta /= kf ;
                                int32_t bs = IsSquare(delta) ;
                                if(bs ) {
                                    uint64_t  k=ks2*k0 ;
                                    uint32_t b = (uint32_t) delta*b0 ;
                                    //                                    bigInt141 n = k02*ks4*a3*b+k*b02*delta*delta ;
                                    bigInt141 n = k*b*(a3*k+b);
                                    if( n < pb141_max && PGCD(bs*b0,a)==1 ) {
                                        nbSol=AddSol141(nbSol,sol,n,a,b,k) ;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bigInt141 Sum = 0 ;
    uint64_t SumAsk = 0 ;
    qsort(sol,nbSol,sizeof(sol[0]),CmpSol) ;
    int i ;
    char str[42];
    for(i=0;i<nbSol;i++) {
        Sum += sol[i].n  ;
        if(sol[i].n <= PB141_MAX_ASK ) SumAsk += sol[i].n ;
        if(pbR->isVerbose) {
            print_bigInt141(sol[i].n,str);
            fprintf(stdout,"\t PB%s (%d/%d)x%lld\t\t%s\n",pbR->ident,sol[i].a,sol[i].b,sol[i].k,str) ;
        }
    }
    print_bigInt141(Sum,str);
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Sum(n)=%s ; n<10**%d\n",pbR->ident,str,ex141_max);
    
    
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",SumAsk) ;
    return 1 ;
}

