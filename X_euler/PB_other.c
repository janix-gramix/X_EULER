//
//  PB_other.c
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#define M_PI 3.14159265358979323846264338327950288

#include "euler_utils.h"
#include "DC3.h"
#include "PB_other.h"
#include "ukkonen.h"

#define PB495_NCOLOR    2
#define PB495_NELBYCOL  1
#define PB495_NPACKET  2

typedef union CONF495 {
    u_int8_t  P[PB495_NPACKET] ;
    u_int64_t N ;
} CONF495 ;

int Cmp_P(const void *el1,const void *el2) {
    return (int)((u_int8_t *)el1)[0] - (int)((u_int8_t *)el2)[0] ;
}
int Cmp_N(const void *el1,const void *el2) {
    if(((u_int64_t *)el1)[0] > ((u_int64_t *)el2)[0]) {
        return 1 ;
    } else if(((u_int64_t *)el1)[0] == ((u_int64_t *)el2)[0]) {
        return 0 ;
    } else return -1 ;
}
int PB495(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nbConf = 1;
    int i ;
    for(i=0;i<PB495_NCOLOR;i++) {
        nbConf *= PB495_NPACKET ;
    }
    CONF495 *CF = calloc(nbConf,sizeof(CF[0])) ;
    int i0,i1,i2,i3,i4,i5 ;
    CONF495 C0,C1,C2,C3,C4,C5 ;
    int N0,N1,N2,N3,N4,N5 ;
    for(i0=0;i0<PB495_NPACKET;i0++) {
        C0.N=0 ; C0.P[i0] = 1 ;
        N0 = i0 ;
        for(i1=0;i1<PB495_NPACKET;i1++) {
            C1=C0 ; C1.P[i1] |=2 ;
            N1=PB495_NPACKET*N0+i1 ;
            if(PB495_NCOLOR == 2) {
                qsort(C1.P,PB495_NPACKET,sizeof(C1.P[0]),Cmp_P) ;
                CF[N1] = C1 ;
                continue ;
            }
           for(i2=0;i2<PB495_NPACKET;i2++) {
                C2=C1 ; C2.P[i2] |= 4 ;
                N2=PB495_NPACKET*N1+i2 ;
                if(PB495_NCOLOR == 3) {
                    qsort(C2.P,PB495_NPACKET,sizeof(C2.P[0]),Cmp_P) ;
                    CF[N2] = C2 ;
                    continue ;
                }
                for(i3=0;i3<PB495_NPACKET;i3++) {
                    C3=C2 ; C3.P[i3] |= 8 ;
                    N3=PB495_NPACKET*N2+i3 ;
                    if(PB495_NCOLOR == 4) {
                        qsort(C3.P,PB495_NPACKET,sizeof(C3.P[0]),Cmp_P) ;
                        CF[N3] = C3 ;
                        continue ;
                    }
                    for(i4=0;i4<PB495_NPACKET;i4++) {
                        C4=C3 ; C4.P[i4] |= 16 ;
                        N4=PB495_NPACKET*N3+i4 ;
                        if(PB495_NCOLOR == 5) {
                            qsort(C4.P,PB495_NPACKET,sizeof(C4.P[0]),Cmp_P) ;
                            CF[N4] = C4 ;
                            continue ;
                        }
                        for(i5=0;i5<PB495_NPACKET;i5++) {
                            C5=C4 ; C5.P[i5] |= 32 ;
                            N5=PB495_NPACKET*N4+i5 ;
                            if(PB495_NCOLOR == 6) {
                                qsort(C5.P,PB495_NPACKET,sizeof(C5.P[0]),Cmp_P) ;
                                CF[N5] = C5 ;
                                continue ;
                            }
                        }
                    }
                }
            }
        }
    }
    qsort(CF,nbConf,sizeof(CF[0]),Cmp_N) ;
    u_int64_t ant = CF[0].N  ;
    int nb = 1 ;
    int nbDiff = 1 ;
    for(i=1;i<nbConf;i++) {
        if(CF[i].N == ant) nb++ ;
        else {
            printf(" %d x %llx -",nb,ant) ;
            nb = 1 ;
            ant = CF[i].N ;
            nbDiff++ ;
        }
    }
    printf(" %d x %llx \n NbDiff=%d\n",nb,ant,nbDiff) ;

    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbConf);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB500_MAX   500500
#define PB500_MOD 500500507
int Cmp500(const void *e1,const void *e2){
    return ((int *)e1)[0] - ((int *)e2)[0] ;
}
int PB500(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrimeNb(PB500_MAX)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    int lastP = tbPrime[nbPrime-1] ;
    int lastP2 = Sqrt32(lastP)+1 ;
    int *powP2=malloc(lastP2*sizeof(powP2[0])) ;
    int nbPow = 0 ;
    int i,p;
    for(i=0;(p=tbPrime[i]) < lastP2;i++) {
        do {
            p = p*p ;
            powP2[nbPow++]=p ;
        } while(p<lastP2) ;
    }
    powP2[nbPow] = lastP ;
    qsort(powP2, nbPow,sizeof(powP2[0]), Cmp500);
    int ip,ipp,is;
    int64_t res=1 ;
    for(ip=0,ipp=0,is=0;is<PB500_MAX;is++) {
        if(tbPrime[ip]<powP2[ipp]) {
            res = (res * tbPrime[ip++]) % PB500_MOD ;
        } else {
            res = (res * powP2[ipp++]) % PB500_MOD ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s last prime (n°%d)=%d last p^n(N°%d)=%d\n",pbR->ident, ip,tbPrime[ip-1],ipp,powP2[ipp-1] );
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",res) ;
    return 1 ;
}


#define PB579_SIZE  5000
#define PB579_MOD   1000000000

typedef enum ClassV {
    Vabc = 0 ,/* Vaab = 1 , */ V0ab = 2 , V00c = 3
} ClassV ;
typedef struct V3 {
    int32_t x ;
    int32_t y ;
    int32_t z ;
} V3 ;

typedef struct V3GCD {
    V3 T ;
    int32_t gcd ;
    ClassV cl ;
    
} V3GCD ;

typedef struct cubeO {
    V3 A ;
    V3 B ;
    V3 C ;
} cubeO ;

#define PB579_MAXT  20000

typedef struct CTX_Cube {
    int nbC ;
    int nbCand ;
    uint64_t nbVertice ;
    uint64_t nbCube ;
} CTX_Cube ;

static inline V3 PVect( V3 T1,V3 T2,int N) {
    V3 V ;
    V.x = (T1.y*T2.z - T1.z*T2.y)/N ;
    V.y = (T1.z*T2.x - T1.x*T2.z)/N ;
    V.z = (T1.x*T2.y - T1.y*T2.x)/N ;
    return V ;
}
// vectoriel product, then abs and x<= y <= z
static inline V3 PVectAbs( V3 T1,V3 T2,int N) {
    V3 V ;
    V.x = (T1.y*T2.z - T1.z*T2.y)/N ; if(V.x < 0) V.x = - V.x ;
    V.y = (T1.z*T2.x - T1.x*T2.z)/N ; if(V.y < 0) V.y = - V.y ;
    V.z = (T1.x*T2.y - T1.y*T2.x)/N ; if(V.z < 0) V.z = - V.z ;
    if(V.x > V.y) { int tmp = V.y ; V.y = V.x ;  V.x = tmp ;  }
    if(V.y > V.z) { int tmp = V.z ; V.z = V.y ;  V.y = tmp ;
        if(V.x > V.y) { int tmp = V.y ; V.y = V.x ; V.x = tmp ; }
    }
    return V ;
}

#define IABS(x) ((x)>=0 ? (x) : -(x))
void AddCube(CTX_Cube * CC,V3 T1,V3 T2,int32_t N,int fact) {
    cubeO cube ;
    cube.A = T1 ;
    cube.B = T2 ;
    cube.C = PVect(T1,T2,N) ;
    CC->nbCand++ ;
    int offx, offy,offz ;
    int sizex, sizey,sizez ;
    {
        int min, max ;
        min=0 ; max=0 ;
        if(cube.A.x< 0) min += cube.A.x ;
        else if(cube.A.x>0) max += cube.A.x ;
        if(cube.B.x< 0) min += cube.B.x ;
        else if(cube.B.x>0) max += cube.B.x ;
        if(cube.C.x< 0) min += cube.C.x ;
        else if(cube.C.x>0) max += cube.C.x ;
        offx = -min ;
        sizex = max + offx ;
        
        min=0 ; max=0 ;
        if(cube.A.y< 0) min += cube.A.y ;
        else if(cube.A.y>0) max += cube.A.y ;
        if(cube.B.y< 0) min += cube.B.y ;
        else if(cube.B.y>0) max += cube.B.y ;
        if(cube.C.y< 0) min += cube.C.y ;
        else if(cube.C.y>0) max += cube.C.y ;
        offy = -min ;
        sizey = max + offy ;
        
        min=0 ; max=0 ;
        if(cube.A.z< 0) min += cube.A.z ;
        else if(cube.A.z>0) max += cube.A.z ;
        if(cube.B.z< 0) min += cube.B.z ;
        else if(cube.B.z>0) max += cube.B.z ;
        if(cube.C.z< 0) min += cube.C.z ;
        else if(cube.C.z>0) max += cube.C.z ;
        offz = -min ;
        sizez = max + offz ;
        
    }
    int sizexyz = (sizex  > sizey ) ? sizex : sizey ;
    if(sizez > sizexyz) sizexyz = sizez ;
    
    if(sizexyz <= PB579_SIZE) {
        int64_t k ;
        int d1 = PGCD(PGCD(IABS(cube.A.x),IABS(cube.A.y)),IABS(cube.A.z)) ;
        int d2 = PGCD(PGCD(IABS(cube.B.x),IABS(cube.B.y)),IABS(cube.B.z)) ;
        int d3 = PGCD(PGCD(IABS(cube.C.x),IABS(cube.C.y)),IABS(cube.C.z)) ;
        CC->nbC++ ;
        for(k=1;k*sizexyz <= PB579_SIZE;k++) {
#if defined(PB579_MOD)
            int64_t nbCubes = ( (((PB579_SIZE-k*sizex+1)*(PB579_SIZE-k*sizey+1)) % PB579_MOD)
                               *( ((PB579_SIZE-k*sizez+1)* fact) % PB579_MOD)
                              ) %  PB579_MOD ;
            CC->nbCube = (CC->nbCube + nbCubes) % PB579_MOD ;
            int64_t nbVertices = ((k*N+1) * ( ((k*k)% PB579_MOD ) * ((N*N)% PB579_MOD ) + (( (d1+d2+d3-N)*k )  % PB579_MOD) + 1)) % PB579_MOD ;
            CC->nbVertice = ( CC->nbVertice+((nbCubes * nbVertices) % PB579_MOD)) % PB579_MOD ;
#else
            int64_t nbCubes = (PB579_SIZE-k*sizex+1)*(PB579_SIZE-k*sizey+1)*(PB579_SIZE-k*sizez+1) * fact ;
            CC->nbCube += nbCubes;
            // (k*N+1) (k*k*N*N+ (d1+d2+d3-N)*k +1 ) pour k=1
            int64_t nbVertices = (k*N+1) * (k*k*N*N + (d1+d2+d3-N)* k +1 ) ;
            CC->nbVertice += nbVertices * nbCubes  ;
#endif
        }
    }
   
}

static inline V3 R1(V3 T) {  V3 Tn   ; Tn.x = T.y    ; Tn.y = T.z    ; Tn.z = T.x    ; return Tn ; }
static inline V3 R2(V3 T) {  V3 Tn   ; Tn.x = T.z    ; Tn.y = T.x    ; Tn.z = T.y    ; return Tn ; }
static inline V3 Pxy(V3 T){  V3 Tn   ; Tn.x = T.y    ; Tn.y = T.x    ; Tn.z = T.z    ; return Tn ; }
static inline V3 Pyz(V3 T){  V3 Tn   ; Tn.x = T.x    ; Tn.y = T.z    ; Tn.z = T.y    ; return Tn ; }
static inline V3 Pxz(V3 T){  V3 Tn   ; Tn.x = T.z    ; Tn.y = T.y    ; Tn.z = T.x    ; return Tn ; }

#define Sx(T)   (T).x = - (T).x
#define Sy(T)   (T).y = - (T).y
#define Sz(T)   (T).z = - (T).z


#define SC_Sx(S,T)  (-(S).x*(T).x + (S).y*(T).y + (S).z*(T).z)
#define SC_Sy(S,T)  ((S).x*(T).x - (S).y*(T).y + (S).z*(T).z)
#define SC_Sz(S,T)  ((S).x*(T).x + (S).y*(T).y - (S).z*(T).z)

void VerifCube(CTX_Cube *CC,V3GCD *tbTrip,int nbT, int i1, int i2,int N) {
    if(tbTrip[i1].gcd > 1 && tbTrip[i2].gcd > 1 && PGCD(tbTrip[i1].gcd,tbTrip[i2].gcd) > 1) return ;
    V3 T1 = tbTrip[i1].T ;
    V3 T2[6] ;
    T2[0] = tbTrip[i2].T ;
    if(tbTrip[i1].cl == V00c) {
        if(tbTrip[i2].cl == V0ab ) {
            AddCube(CC,R1(T2[0]),T1,N,6);
        }
        return ;
    } else if(tbTrip[i2].cl == V00c) {
        if(tbTrip[i1].cl == V0ab ) {
            AddCube(CC,R1(T1),T2[0],N,6); // pas de test a faire cela est OK
        }
        return ;
    }
    T2[1] = R1(tbTrip[i2].T) ; T2[2] = R2(tbTrip[i2].T) ; T2[3] = Pxy(tbTrip[i2].T) ;
    T2[4] = Pyz(tbTrip[i2].T); T2[5] = Pxz(tbTrip[i2].T) ;
    int ip;
    for(ip=0;ip<6;ip++){
         if(SC_Sy(T2[ip],T1)==0) { Sy(T2[ip]) ;  break ; }
         if(SC_Sz(T2[ip],T1)==0) { Sz(T2[ip]) ;  break ; }
    }
    if(ip == 6) return ;
    else {
        V3 T3 = PVectAbs(T1,T2[ip],N) ;;
        int i3 ;
        for(i3=i2;i3<nbT;i3++) {
            if(tbTrip[i3].T.x < T3.x) continue ;
            if(tbTrip[i3].T.x > T3.x) break ;
            if(T3.y==tbTrip[i3].T.y ){
               int fact ;
                if(i1==i2 || i2==i3 ) {
                    if(i3==i1) { // i1=i2=13
                        fact = (N==3) ? 4 : 8 ;
                    } else { //
                        fact = 12 ;
                    }
                } else {
                    fact = 24 ;
                }
                AddCube(CC,T1,T2[ip], N, fact);
               break ;
            }
        }
    }
}

int PB579(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_Cube CC ;
    V3GCD tbTrip[PB579_MAXT] ;
    int32_t n1,n2,n3;
    int32_t N ;
    CC.nbC = 0 ;
    // on ajoute les cubes alignes sur le reseau
    N = PB579_SIZE ;
//
    CC.nbVertice = 0 ;
#if defined(PB579_MOD)
    int64_t k ;
    for(k=1;k<=N;k++) {
        CC.nbVertice = (CC.nbVertice + ( ( (((N+2-k)*(N+2-k)*(N+2-k)) % PB579_MOD) * ((k*k*k) % PB579_MOD) ) % PB579_MOD)) % PB579_MOD ;
    }
    CC.nbCube =  ((((int64_t)N)*N * (N+1) * (N+1)) / 4 ) % PB579_MOD ;
#else
    CC.nbCube =  ((((int64_t)N)*N * (N+1) * (N+1)) / 4 ) ;
    CC.nbVertice =(((uint64_t)(N))*(N+1)*(3 * N*N*N*N*N + 39 * N*N*N*N + 213 * N*N*N + 627 * N*N + 640 * N + 158 ))/ 420 ;
#endif
    for(N=3;N<=PB579_SIZE;N+=2) {
        int32_t S1 = N*N ;
        int nbT = 0 ;
        int n1Max = Sqrt32(S1/3) ;
        for(n1=0;n1<=n1Max;n1++) {
            int32_t S2 = S1 - n1*n1 ;
            int n2Max =Sqrt32(S2>>1) ;
            n2=n1 ;
            int Diff = S2-n2*n2 ;
            n3=Sqrt32(S2-n2*n2) ;
            Diff -= n3*n3 ;
            for(; n2<=n2Max ; ) {
                while(Diff < 0) {
                     Diff += 2*n3-1 ; n3-- ;
                }
                if(Diff==0) { //  n3 =Sqrt32(S2-n2*n2) ;
                    int gcd  = PGCD(n1,n2);
                    tbTrip[nbT].T.x = n1 ;
                    tbTrip[nbT].T.y = n2 ;
                    tbTrip[nbT].T.z = n3 ;
                    tbTrip[nbT].gcd = PGCD(gcd,n3); ;
                    tbTrip[nbT].cl = (n1 != 0 ) ? Vabc : ( (n2==0) ? V00c : V0ab  ) ;
                    nbT++;
                }
                do {
                    Diff -= 2*n2+1 ; n2++  ;
                } while(Diff > 0) ;
            }
        }
        if(nbT > 0) {
            int i1,i2 ;
            if((N % 512) == 1)printf("%d->%d\n",N,nbT);
            for(i1=0;i1<nbT;i1++){
                for(i2=i1;i2<nbT;i2++) {
                    VerifCube(&CC,tbTrip,nbT,i1,i2,N) ;
                }
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s NbCubeType=%d NbCubes=%lld NBVertices=%lld \n"
                              ,pbR->ident,CC.nbC,CC.nbCube,CC.nbVertice) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",CC.nbVertice);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB620_MAX   500 //
int Count620(int s, int p,int q) {
    int k ;
    double a = -(double)(s+p) / (double) (s+q) ;
    double piInv = 1 / M_PI ;
    double factSin = ((double) (s+p)) / (s+q) ;
    double factSin2 = factSin*factSin ;
    int kmin = 1, kmax = s+p-1 ;
    while(kmax-kmin > 0) {
//        printf("(%d-%d)",kmin,kmax);
        k = (kmin + kmax+1) / 2 ;
        int i ;
        double b = (double) k / ((double) (s+q)) ;
        double x1 = 0 ;
        double y1 ,dy1 ,dx1 ;
        for(i=0;i<10;i++) {
            y1 = piInv * asin(factSin * sin(M_PI*x1)) - a*x1 - b ;
            dy1 = cos(M_PI*x1) * factSin / sqrt(1-factSin2*sin(M_PI*x1)*sin(M_PI*x1)) - a ;
            dx1 = - y1/dy1 ;
            x1 += dx1 ;
        } ;
        double distx2 = p+q - (s+q) * cos(M_PI * (a*x1 + b )) + (s+p) * cos(M_PI * x1) ;
        if(distx2 > M_PI*2) {
            kmin =k ;
        } else {
            if(k==kmax) {
                kmax-- ;
            } else {
                kmax = k ;
            }
        }
    }
    return kmin ;
}

int PB620(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nbSol = 0 ;
    int s,p,q,N ;
    for(N=16;N<=PB620_MAX;N++) {
        for(s=5;s<=N-11;s++) {
            
            for(p=5;2*p+1<=N-s;p++) {
                q = N-s-p ;
                nbSol+= Count620(s,p,q) ;
            }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB620a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    double piInv = 1 / M_PI ;
    int nbSol = 0 ;
    int s,p,q,N ;
    for(N=16;N<=PB620_MAX;N++) {
        for(s=5;s<=N-11;s++) {
            
            for(p=5;2*p+1<=N-s;p++) {
                q = N-s-p ;
                double Smax2 = ((p+q)-2*M_PI)  ;
                double x = acos(((s+q)*(s+q)-Smax2*Smax2-(s+p)*(s+p))/(2*Smax2 *(s+p)))*piInv ;
                double y = 1-acos(((s+p)*(s+p)-Smax2*Smax2-(s+q)*(s+q))/(2*Smax2 *(s+q)))*piInv ;
                int k = ((s+q)*y + (s+p)*x) ;
                nbSol+= k ;
            }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSol);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



typedef struct P_EXP {
    uint32_t   p;
    uint32_t   exp ;
    
} P_EXP ;
#define PB622_NBP   11 // number of prime factors for 2**60-1
int PB622(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i, jp,k;
    uint64_t S=0;
    P_EXP fact2exp60_m1[PB622_NBP] = {{3, 2}, {5, 2}, {7, 1}, {11, 1}, {13, 1}, {31, 1}, {41, 1}, {61,1}, {151, 1}, {331, 1}, {1321, 1}};
    P_EXP fact2div60byP[3][PB622_NBP] ;
    uint32_t quotPrime60[3] = { (1<<(60/2))-1 , (1<<(60/3)) -1 , (1<<(60/5)) -1 } ;
    {
        for(i=0;i<PB622_NBP;i++) {
            uint32_t p = fact2exp60_m1[i].p ;
            for(jp=0;jp<3;jp++) {
                int exp, N ;
                for(exp=0, N=quotPrime60[jp];(N % p )== 0; exp++ , N /= p) ;
                fact2div60byP[jp][i].p = p; fact2div60byP[jp][i].exp = exp ;
            }
        }
    }
    {  // loop to test all factor of 2**60-1
        int32_t exp[PB622_NBP] , ip=PB622_NBP-1 ;
        for(i=0;i<PB622_NBP;i++) exp[i] = 0 ;
        while(ip >= 0) {
            if(exp[ip]++ < fact2exp60_m1[ip].exp) {
                while(ip<PB622_NBP-1) exp[++ip] = 0 ;
                for(jp=0;jp<3;jp++) { // test not divisor of element quotPrime60
                    for(k=0;k<PB622_NBP;k++) {
                        if(fact2div60byP[jp][k].exp < exp[k] ) break ; // not divisor
                    }
                    if(k==PB622_NBP) break ; // divisor of quotPrime60[jp]
                }
                if(jp==3) { // 60 is the minimum value, good candidate
                    uint64_t N ;
                    for(k=0,N=1;k<PB622_NBP;k++) {
                        int ie ;
                        for(ie=exp[k];ie>0;ie--) N *= fact2exp60_m1[k].p ;
                    }
                    S += N+1 ;
                }
            } else {  ip-- ; }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB625_MAX   100000000000LL
// #define PB625_MAX   10000000LL //
// #define PB625_MAX   100000000LL //
// #define PB625_MAX   1000000L
#define PB625_MOD 998244353
// #define PB625_MOD 100000




int PB625(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ,p ,j , N = PB625_MAX ;
#if 1
    if(PB625_MAX > 1000000000L) {
        printf("Too memory for %lld\n",PB625_MAX) ;
        return 0 ;
    }
    uint64_t *sumPhi=calloc(PB625_MAX+1,sizeof(sumPhi[0])) ;
    uint64_t nbP = 1 ;
    for(p=2;p<=PB625_MAX;p++)  {
        if(sumPhi[p]) {
            nbP += sumPhi[p] ;
            nbP =  nbP % PB625_MOD ;
            sumPhi[p] = nbP ;
//            printf("(%d,%d)",p,sumPhi[p]) ;
            continue ;
        } else {
            int pm = p-1 ;
            nbP += pm ;
            nbP =  nbP % PB625_MOD ;
            sumPhi[p] = nbP ;
            int k ;
//            printf("(%d,%d)",p,sumPhi[p]) ;
            for(j=2*p,k=2;j<=PB625_MAX;j+=p,k++) {
                if(sumPhi[j]) {
                    sumPhi[j] = (sumPhi[j]/p)*(p-1)  ;
                } else {
                    sumPhi[j] = k*(p-1)   ;
                }
                
            }
        }
        
    }
    
    uint64_t isum =  (N+(N/2)+1)  ;
    uint64_t idif =  (N-(N/2))  ;
    if(isum & 1) idif /=2 ;
    else isum /=2 ;
    idif = idif % PB625_MOD ;
    isum = isum % PB625_MOD ;

    uint64_t S = (isum*idif) % PB625_MOD ;
    //    printf("SN=%lld ",S );
    for(i=1;2*i<=N;i++) {
//        if(i==plim) printf("(i=%d S=%lld)",i,S);
        S += ((uint64_t)(i % PB625_MOD ) *  sumPhi[N/i] ) % PB625_MOD;
        S = S % PB625_MOD ;
//        printf("(%d,%lld,%lld)",i,sumPhi[N/i],S);
    }
    if(pbR->isVerbose)fprintf(stdout, "\t PB%s MAX=%lld[%d] S=%lld\n",pbR->ident,PB625_MAX,PB625_MOD,S );
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    free(sumPhi) ;
#else
    uint64_t S = 0 ;
    for(i=1;i<=N;i++) {
        for(j=1;j<=i;j++) {
            S += PGCD(i,j);
        }
    }
    if(pbR->isVerbose)fprintf(stdout, "\t PB%s MAX=%lld[%d] S=%lld\n",pbR->ident,PB625_MAX,PB625_MOD,S );
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
#endif
     pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
typedef struct P625  {
    uint64_t p ;
    uint64_t q ;
    uint64_t ioffset ;
} P625 ;


typedef struct N625 {
    uint64_t   gcd ;
    uint64_t   index ;
} N625 ;

int PB625a(PB_RESULT *pbR) {
    int sizeD = pow(PB625_MAX, 0.625) ;
    pbR->nbClock = clock() ;
    int64_t  N = PB625_MAX ;
    uint32_t *sumPhi=calloc(sizeD,sizeof(sumPhi[0])) ;
    uint64_t *sumPhi_Ext=malloc(sizeD*sizeof(sumPhi_Ext[0])) ;
    
    { // comp sumPhi[n] = Sigma[1..n] Phi[n] ,Phi[n] = #[1..n] prime with n
        sumPhi[1] = 1 ;
        int32_t nbP = 1 ;
        int j,p; // limited by sizeD, int32 is sufficient.
        for(p=2;p<sizeD;p++)  {
            if(sumPhi[p]) {
                nbP += sumPhi[p] ;
                nbP = nbP % PB625_MOD ;
                sumPhi[p] = (int32_t) nbP ;
                continue ;
            } else {
                uint32_t pm = p-1 ;
                nbP += pm ;
                nbP = nbP % PB625_MOD ;
                sumPhi[p] = (int32_t) nbP ;
                int q ;
                for(j=2*p,q=2;j<sizeD;j+=p,q++) {
                    if(sumPhi[j]) {
                        sumPhi[j] = (sumPhi[j]/p)*(p-1)  ;
                    } else {
                        sumPhi[j] = q*(p-1)   ;
                    }
                }
            }
        }
    }
    uint64_t S1 = 0 , S2 = 0 ,S=0;
    uint64_t i0,i1 ;
    {
        int32_t d ;
        // calcul for i in [2,SizeD] sum { sumPhi[i] * (T(N/i) - T(N/(i+1)) }
        for(d=1;d<=sizeD;d++) {
            i0 = N/d ;
            i1 = N/(d+1) ;
            uint64_t isum = ((i0+i1+1) & 1 ) ? ((i0+i1+1) % PB625_MOD) : (((i0+i1+1)/2) % PB625_MOD) ;
            uint64_t idif = ((i0+i1+1) & 1 ) ? (((i0-i1)/2) % PB625_MOD) : ((i0-i1) % PB625_MOD) ;
            uint64_t nb= (isum*idif ) % PB625_MOD ;
            S1 += (nb * sumPhi[d]) % PB625_MOD ;
//            if(d==1) printf("Sn=%lld",S1);
        }
    }
//    printf(" S1=%lld ",S1);
    uint64_t id = N/sizeD ;
    while(id > 0) {
        uint64_t Nd = N / id ;
        uint64_t k ;
        i0 = (Nd & 1) ? (Nd % PB625_MOD) : ((Nd/2) % PB625_MOD) ;
        i1 = (Nd & 1) ? (((Nd+1)/2) % PB625_MOD) : ((Nd+1) % PB625_MOD) ;
        // calcul de sumPhi[Nd] stored in sumPhi_ext[id]
        // (Nd*(Nd+1))/2 = SumPhi[Nd] + sigma{ k in [2..Nd] sumPhi[Nd/k]}
        // SumPhi[Nd] = (Nd*(Nd+1))/2
        // - sigma{ 1<k<=Nd/sizeD ; sumPhi[N/(k*id)] := sumPhi_Ext[k*id] }
        // - sigma{ Nd/sizeD < k <= Nd; sumPhi[Nd/k } . This sum is cut in two parts
        // to take in account that Nd/k takes many time the same value for big k.
        uint64_t Si = (i0*i1) % PB625_MOD ; // Nd*((Nd+1)/2
        for(k=2;k<=Nd/sizeD;k++) {
            Si += PB625_MOD - sumPhi_Ext[k*id] ; //  Si -= sumPhi_Ext[k*id] ;
        }
        for(;(k+1)*(k+1) < Nd ;k++) {
            Si += PB625_MOD -  sumPhi[Nd/k] ; //   Si -= sumPhi[Nd/k] ;
        }
        for(k=Nd/k;k>0 ;k--) { //   Si -= sumPhi[k] * (Nd/k - Nd/(k+1)) ;
            Si += PB625_MOD -  ((sumPhi[k] * (Nd/k - Nd/(k+1))) % PB625_MOD) ;
        }
        sumPhi_Ext[id] = Si % PB625_MOD ;
        S2 += (id * sumPhi_Ext[id]) % PB625_MOD ;
        S2 = S2 % PB625_MOD ;
        id -- ;
    }
    S =(S1+S2) % PB625_MOD ;
    if(pbR->isVerbose)fprintf(stdout, "\t PB%s MAX=%lld[%d] S=%lld\n",pbR->ident,PB625_MAX,PB625_MOD,S );
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    free(sumPhi) ; free(sumPhi_Ext) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


// this version dont use the recursion for sumPhi and use Sieve to the end so
// is very slow for big values.
int PB625b(PB_RESULT *pbR) {
    int pLim = (Sqrt64(PB625_MAX)+1) ;
    int sizeD = pLim ;
    pbR->nbClock = clock() ;
    int32_t j ;
    uint64_t p , N = PB625_MAX ;
    uint64_t *sumPhi=calloc(sizeD,sizeof(sumPhi[0])) ;
    uint64_t *sumPhi2=malloc(sizeD*sizeof(sumPhi2[0])) ;
    N625 *nxt = malloc(sizeD*sizeof(nxt[0])) ;
    int64_t nbPrimeAlloc = sizeD ;
    P625 * tbPrime = malloc(nbPrimeAlloc * sizeof(tbPrime[0])) ;
    int nbPrime = 0 ;
    
    uint64_t nbP = 1 ;
    sumPhi[1] = 1 ;
    {
        int j,p; // limited by sizeD, int3é is sufficient.
        for(p=2;p<sizeD;p++)  {
            if(sumPhi[p]) {
                nbP += sumPhi[p] ;
                nbP = nbP % PB625_MOD ;
                sumPhi[p] = nbP ;
                continue ;
            } else {
                tbPrime[nbPrime].p = p ;
                uint64_t pm = p-1 ;
                nbP += pm ;
                nbP = nbP % PB625_MOD ;
                sumPhi[p] = nbP ;
                int q ;
                for(j=2*p,q=2;j<sizeD;j+=p,q++) {
                    if(sumPhi[j]) {
                        sumPhi[j] = (sumPhi[j]/p)*(p-1)  ;
                    } else {
                        sumPhi[j] = q*(p-1)   ;
                    }
                    
                }
                tbPrime[nbPrime].q = q ;
                tbPrime[nbPrime++].ioffset = j - sizeD ;
            }
        }
    }
    uint64_t isum =  (N+(N/2)+1) ;
    uint64_t idif =  (N-(N/2))  ;
    if(isum & 1) idif /=2 ;
    else isum /=2 ;
    idif = idif % PB625_MOD ;
    isum = isum % PB625_MOD ;
    uint64_t S = (isum*idif) % PB625_MOD ;
//    printf("SN=%lld ",S );
    
    uint64_t pl ;
    uint64_t S1 = 0 , S2 = 0 ;
    uint64_t i0,i1 ;

    {
        int32_t i ;
        for(i=2;i<=sizeD;i++) {
            i0 = N/i ;
            i1 = N/(i+1) ;
            uint64_t isum = ((i0+i1+1) & 1 ) ? ((i0+i1+1) % PB625_MOD) : (((i0+i1+1)/2) % PB625_MOD) ;
            uint64_t idif = ((i0+i1+1) & 1 ) ? (((i0-i1)/2) % PB625_MOD) : ((i0-i1) % PB625_MOD) ;
            uint64_t nb= (isum*idif ) % PB625_MOD ;
            S1 += (nb * sumPhi[i]) % PB625_MOD ;
            S1 = S1 % PB625_MOD ;
        }
    }
//    printf(" S1=%lld ",S1);
    uint64_t pgNext = i0 ;

    uint64_t invPgNext = N /pgNext ;
    while(invPgNext < sizeD) {
        S1 += ((uint64_t)pgNext * sumPhi[invPgNext]) % PB625_MOD;
        S1 = S1 % PB625_MOD ;
        pgNext-- ;
        invPgNext = N /pgNext ;
    }
    int nbNxt ;
    for(nbNxt=0;pgNext>0;nbNxt++ ) {
        nxt[nbNxt].gcd = pgNext ;
        nxt[nbNxt].index = N / pgNext ;
        pgNext-- ;
    }
    int indNxt = 0 ;
    for(pl=sizeD;pl<=PB625_MAX;pl+=sizeD ){
        memset(sumPhi2,0,sizeD*sizeof(sumPhi2[0])) ;
        int32_t ip ;
        for(ip=0;ip<nbPrime;ip++) {
            if(tbPrime[ip].ioffset < sizeD) {
                p = tbPrime[ip].p ;
                uint64_t q = tbPrime[ip].q ;
                uint64_t j ;
                for(j=tbPrime[ip].ioffset;j<sizeD;j+=p,q++) {
                    
                    if(sumPhi2[j]) {
                        sumPhi2[j] = ((sumPhi2[j])/p)*(p-1)  ;
                    } else {
                        sumPhi2[j] = q*(p-1)   ;
                    }
                }
                tbPrime[ip].q = q ;
                tbPrime[ip].ioffset = j - sizeD ;
            } else {
                tbPrime[ip].ioffset -= sizeD ;
            }
        }
        for(j=0;j<sizeD;j++) {
            if(sumPhi2[j]) {
                nbP += sumPhi2[j] % PB625_MOD ;
            } else {
                p = pl+j ;
                nbP +=  ( p -1 ) % PB625_MOD;
                int k ;
                for(k=indNxt;k<nbNxt;k++) {
                    int xp = nxt[k].index /p ;
                    S2 += ((sumPhi[xp]-1)*nxt[k].gcd ) % PB625_MOD;
                    S2 = S2 % PB625_MOD ;
                }
            }
            nbP = nbP % PB625_MOD ;
            if(pl+j == nxt[indNxt].index){
                S += ( nbP * nxt[indNxt].gcd ) % PB625_MOD ;
                S = S % PB625_MOD ;
                indNxt++ ;
                if(indNxt >= nbNxt) break ;
            }
        }
        if(indNxt >= nbNxt)  break ;
    }
    S =(S+S1-S2) % PB625_MOD ;
    if(pbR->isVerbose)fprintf(stdout, "\t PB%s MAX=%lld[%d] S=%lld\n",pbR->ident,PB625_MAX,PB625_MOD,S );
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB626_DIM   3

int Cmp_626Mat(const void * el1, const void * el2) {
    int64_t * mat1 = (int64_t *) el1 ;
    int64_t * mat2 = (int64_t *) el2 ;
    if(mat1[0]>mat2[0]) return 1 ;
    else return (mat1[0]- mat2[0]) ? -1 : 0 ;
}


int PB626(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nbMat = (1LL << (PB626_DIM * PB626_DIM)) ;
    int64_t *classMat = calloc(nbMat,sizeof(classMat[0])) ;
    int64_t *classMatPerm = calloc(nbMat,sizeof(classMatPerm[0])) ;
    int64_t imaskL = (1 << PB626_DIM) -1 ;
    int64_t imaskC = 1 ;
    int64_t i ;
    for(i=1;i<PB626_DIM;i++) {
        imaskC = (imaskC << PB626_DIM) | 1 ;
    }
    int nbClass =0 ;
    int nbBitMax =0 ;
    for(i=0;i<nbMat;i++) {
        if(classMat[i]) continue ;
        printf("\n%o->",i);
        int k ;
        nbClass ++ ;
        int nbBitMin = PB626_DIM*PB626_DIM ;
        for(k=0;k< (1 << (2*PB626_DIM-1));k++) {
            int ik,k1 ;
            int64_t imat= i ;
            for(ik=0,k1=k;k1;ik++) {
                if(k1&1) {
                    if(ik < PB626_DIM) {
                        imat ^= imaskL << (ik * PB626_DIM) ;
                    } else {
                        imat ^= imaskC << (ik-PB626_DIM) ;
                    }
                }
                k1 >>= 1 ;
            }
            int ib,nbB=0 ; ;
            for(ib=0;ib<PB626_DIM * PB626_DIM;ib++) {
                if((1LL<< ib) & imat) nbB++ ;
            }
          if(imat < nbMat) classMat[imat] = i+1 ;
            printf("%o,",imat);
            if(nbB< nbBitMin) {
                nbBitMin = nbB ;
            }
        }
        if(nbBitMin >= nbBitMax) { nbBitMax = nbBitMin ; printf("*") ; }
    }
    printf("\nNbClass=%d, nbBitMax=%d\n",nbClass,nbBitMax) ;
    int nbClassPerm=0 ;
    for(i=0;i<nbMat;i++) {
        int k ;
        if(classMatPerm[i]) continue ;
        printf("\n%o->",i);
        nbClassPerm ++ ;
        uint8_t permL[PB626_DIM] ;
        for(k=0;k<PB626_DIM;k++) permL[k] = k ;
        do {
            uint8_t permC[PB626_DIM] ;
            for(k=0;k<PB626_DIM;k++) permC[k] = k ;
            do {
                int64_t iMat = 0 ;
                int ib ;
                for(ib=0;ib<PB626_DIM * PB626_DIM;ib++) {
                    int il = ib / PB626_DIM ;
                    int ic = ib - il * PB626_DIM ;
                    if(i & ( 1LL << ( permL[il] * PB626_DIM + permC[ic])) ) iMat |= (1LL<< ib) ;
                }
                if(classMatPerm[iMat] == 0) {
                    printf("%o(%o),",iMat,classMat[iMat]) ;
                    classMatPerm[iMat] = i+1 ;
                }
                
            } while(NextPermut(permC,PB626_DIM)>=0) ;
        } while(NextPermut(permL,PB626_DIM)>=0) ;
        
 
    }
    printf("\nNbClassPerm=%d\n",nbClassPerm) ;
/*
    int nbBitMax = ((PB626_DIM-2)*(PB626_DIM-1))/2 + 1 ;
    int nbClass = 0 ;
    int nbB ;
    for(nbB=1;nbB<= nbBitMax;nbB++) {
        nbB += 
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbClass);
 */
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



#define PB1000_NUM  10000000
int PB1000(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrimeNb(PB1000_NUM)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    
    fprintf(stdout,"\t PB%s prime[%d] = %u\n",pbR->ident, PB1000_NUM,tbPrime[PB1000_NUM-1]) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%u",tbPrime[PB1000_NUM-1]);
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
//#define PB681_MAXA  1000000LL
#define PB681_MAXA  100000LL
#define PB681_MAXDIV    2000

#define MAX_DIVISOR 20
typedef struct Divisors {
    int P[MAX_DIVISOR];
    int EXP[MAX_DIVISOR] ;
    int nbP ;
    int curE[MAX_DIVISOR] ;
    int64_t curD[MAX_DIVISOR] ;
    int maxD ;
} Divisors ;

int64_t DIV_Init(Divisors * div, int maxD) {
    int i ;
    for(i=0;i<div->nbP;i++) { div->curE[i] = 0 ; div->curD[i] = 1 ; }
    div->maxD = maxD ;
    return 1 ;
}
int64_t DIV_getnext(Divisors * div) {
    int i ;
    for(i=0;i<div->nbP;i++) {
        if(div->curE[i]<div->EXP[i]) {
            int64_t curD=div->curD[i] * div->P[i] ;
            if(curD<=div->maxD) {
               int j  = i ;
               div->curD[j] = curD ;
               div->curE[j--]++ ;
                for(;j>=0;j--) {
                    div->curE[j] = 0 ;
                    div->curD[j] = curD ;
               }
               return curD ;
            }
        }
    }
    div->curD[0] = 0 ;
    return 0 ;
}
int CmpDiv(const void *e1,const void *e2) {
    int64_t i1=((int64_t *)e1)[0];
    int64_t i2 =((int64_t *)e2)[0] ;
    if(i1>i2) return 1;
    else if(i1<i2) return -1 ;
    return 0 ;
}
int PB681a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB681_MAXA+1)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;
    
    int64_t  sumS = 0 ;
    int64_t n ;
    Divisors divL ;
    sumS += 4 ;
    for(n=2;n<=PB681_MAXA;n++) {
        int64_t div_n[10000] ;
        int nbD_n = 0 ;
        div_n[nbD_n++]= 1 ;
        int64_t n1 ;
        int ip ;
        for(ip=0,n1=n;ip<nbPrime;ip++) {
            int64_t p = tbPrime[ip] ;
            if(p*p>n) break ;
            int imin=0 ;
            while ((n1 % p ) == 0 ) {
                int i ;
                int64_t p2 = p * p ;
                int is = imin+nbD_n ;
                for(i=imin;i<imin+nbD_n;i++) {
                    div_n[i+nbD_n]=div_n[i]*p ;
                    div_n[i+2*nbD_n]=div_n[i]*p2 ;
                }
                imin += 2*nbD_n ;
                n1 /= p ;
            }
            if(imin) nbD_n += imin ;
        }
        if(n1 > 1) {
            int i ;
            for(i=0;i<nbD_n;i++) {
                div_n[i+nbD_n]=div_n[i]*n1 ;
                div_n[i+2*nbD_n]=div_n[i]*n1*n1 ;
            }
            nbD_n *= 3 ;
        }
        qsort(div_n,nbD_n,sizeof(div_n[0]),CmpDiv);
        int i ;
 //       printf("\n%lld->",n);
 //       for(i=0;i<nbD_n;i++) printf("%lld,",div_n[i]);
        int il ;
        int64_t nl,nh ;
        for(il=0;il<nbD_n;il++) {
            nl=div_n[il];
            if(nl >n ) break ;
            nh = div_n[nbD_n-1-il] ;
            int ia ;
            int64_t a,b,c,d ;
            int minId ;
         int64_t sqnh = Sqrt64(nh-1)+1 ;
            for(minId=nbD_n-1-il;/*minId >=0 && */ div_n[minId]>=sqnh;minId--) ;
            minId++ ;
            for(ia=0;ia<nbD_n;ia++) {
                a=div_n[ia];
                if(a*a>nl) break ;
                b = nl / a ;
                if(nl != a*b) continue ;
                int64_t Sab = a+b ;
                int id ;
               for(id=minId;id<nbD_n;id++) {
                    d = div_n[id] ;
 //                  if(d*d<nh) {  printf("\n%lld+%lld+%lld+%lld %lld",a,b,c,d,n ) ; continue ; }
                    c = nh/d ;
                   if(c < b) break ;
                    if(d >= Sab+c) break ;
                    if(nh != c*d) continue ;
                    int64_t S =Sab+c+d ;
                    if((S&1)==0) {
                        sumS += S ;
  //                     S >>= 1 ;
  //                      printf("\n%lld+%lld+%lld+%lld S=%lld A=%lld sumS=%lld id=%d"
  //                             ,S-d,S-c,S-b,S-a,S,n,sumS,id) ;
                        

                    }
                }
            }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",sumS);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int PB681(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_PRIMETABLE * ctxP  ;
    if((ctxP = Gen_tablePrime(PB681_MAXA)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP);
    int nbPrime = GetNbPrime(ctxP) ;

    int64_t  sumS = 0 ;
    int64_t n ;
    Divisors divL ;
    for(n=1;n<=PB681_MAXA;n++) {
        divL.nbP = 0 ;
        int64_t n1 ;
        int ip ;
        for(ip=0,n1=n;ip<nbPrime;ip++) {
            int p = tbPrime[ip] ;
            if(p*p>n) break ;
            if ((n1 % p ) == 0 ) {
                divL.P[divL.nbP] = p ;
                divL.EXP[divL.nbP] = 1 ;
                n1 /= p ;
                while ((n1 % p ) == 0 ) {
                    divL.EXP[divL.nbP]++ ;
                    n1 /= p ;
                }
                divL.nbP++ ;
                if(n1==1) break ;
            }
        }
        if(n1 > 1) {
            divL.P[divL.nbP] = n1 ;
            divL.EXP[divL.nbP++] = 1 ;
        }
        int i ;
        for(i=0;i<divL.nbP;i++)divL.EXP[i] <<= 1 ;
        int64_t dl,dh;
        dl = DIV_Init(&divL,n) ;
 //       printf("\n%lld->",n) ;
        do {
            int64_t n2 = n*n ;
            Divisors divH,diva ;
            dh = n2/dl ;
//            if(dl*dh != n2) continue ;
            int64_t a,b,c,d ;
            int i;
            divH.nbP = 0;
            diva.nbP = 0 ;
            for(i=0;i<divL.nbP;i++) {
                if(divL.curE[i] < divL.EXP[i]) {
                    divH.P[divH.nbP] = divL.P[i] ;
                    divH.EXP[divH.nbP++] = divL.EXP[i] - divL.curE[i] ;
                }
                if(divL.curE[i]) {
                    diva.P[diva.nbP] = divL.P[i] ;
                    diva.EXP[diva.nbP++] = divL.curE[i] ;
                }
            }
            int cMax = Sqrt64(dh) ;
            int aMax = Sqrt32(dl) ;
            a = DIV_Init(&diva,aMax);
            do {
                b= dl/a ;
  //              if(a*b != dl) continue ;
                int64_t Sab = a+b ;
                c = DIV_Init(&divH,cMax);
                do {
                    if(c>=b) {
                        d = dh/c ;
 //                   if(c*d != dh) continue ;
                        int64_t S = Sab+c ;
                        if(S<=d) continue ;
                        S +=d ;
                        if((S&1)==0) {
                            sumS += S ;
//                           S >>= 1 ;
//                               printf("%lld+%lld+%lld+%lld S=%lld A=%lld sumS=%lld\n"
//                                      ,S-d,S-c,S-b,S-a,S,n,sumS) ;
                        }
                    }
                    
                }while((c=DIV_getnext(&divH)) != 0) ;
            }while((a=DIV_getnext(&diva)) != 0) ;
        
        }while((dl=DIV_getnext(&divL)) != 0) ;
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",sumS);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int cmpPlace( const void *el1, const void *el2) {
    return ((uint8_t *)el1)[0] - ((uint8_t *)el2)[0] ;
}

static  uint64_t Cnp(int n,int p) {
    if(p==0) return 1 ;
    if(n<p) return 0 ;
    if(p==1) return n ;
    uint64_t cnp = n*(n-1)/2 ;
    if(p>2) {
        int i  ;
        for(i=2;i<p;i++) {
            cnp *= n-i ;
            cnp /= i+1 ;
        }
    }
    return cnp ;
}

static int Index5(int S,int i0,int i1,int i2,int i3) {
    static int Cn3[] = {0,1,5,15,35,70,126,210,330,495,715,1001,1365,1820} ;
    static int Cn2[] = {0,1,4,10,20,35,56,84,120,165,220,286,364,455} ;
    static int Cn1[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91} ;
    int ind = Cn3[S-i0]+Cn2[S-i0-i1]+Cn1[S-i0-i1-i2]+S-i0-i1-i2-i3 ;
    return ind ;
}



#define FLOAT687
#if defined(FLOAT687)
typedef  double PROB687 ;
#define SetIndex5(NR,N0,N1,N2,N3,ant,val)  cur5[Index5(NR,N0,N1,N2,N3)] += ant*(val)
#endif
#   define MAX_NR   13


// insertion ****
#define INS4_N0(ant) SetIndex5(nr,n0,n1,n2,n3,ant,nb0)
#define INS4_N1(ant) if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3,ant,n1)
#define INS4_N2(ant) if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3,ant,2*n2)
#define INS4_N3(ant) if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3-1,ant,2*n3)
#define INS4_N4_2(ant) if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3,ant,n4)
#define INS4_N4_3(ant) if(n4>0)SetIndex5(nr,n0,n1,n2,n3+1,ant,2*n4)

int PB687(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nr ;
    
    PROB687 buf0[2380],buf1[2380] ;
    PROB687 ant, *ant5 = buf0, *cur5=buf1 ;
    
    ant5[Index5(1,0,0,0,0)] = 1 ;
    ant5[Index5(1,1,0,0,0)] = 0 ;
    ant5[Index5(1,0,1,0,0)] = 0 ;
    ant5[Index5(1,0,0,1,0)] = 0 ;
    ant5[Index5(1,0,0,0,1)] = 0 ;
    
    int n0,n1,n2,n3;
    // n0 * * * *
    // n1 ** * *
    // n2 ** **
    // n3 *** *
    // n4 ****
    n1 = 0 ;
    
    for(nr=2;nr<=MAX_NR;nr++) {
        int i ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) cur5[i] = 0 ;
        // pour rajouter un rang : cela peut être : *,*,* || **,* || ***
        for(n0=0;n0<nr;n0++) {
            for(n1=0;n1<nr-n0;n1++) {
                for(n2=0;n2<nr-n0-n1;n2++) {
                    for(n3=0;n3<nr-n0-n1-n2;n3++) {
                        int n4=nr-1-n0-n1-n2-n3  ;
                        int nb0 = 4*n0+3*n1+2*n2+2*n3+n4 ; // not +1, by symetry on ranks we can impose not to insert at the end
                        ant = ant5[Index5(nr-1,n0,n1,n2,n3)] ;
                        if(ant <= 0) continue ;
                        // insert ****
                        INS4_N0(ant) ; INS4_N1(ant) ; INS4_N2(ant) ; INS4_N3(ant) ; INS4_N4_2(ant) ; INS4_N4_3(ant) ;
                        
                        int pow3[4] = { 1,3,9,27} ;
                        int i0,i1,i21,i22,i31,i32,i42,i43,nbp ;
                        int d,d0,d1,d2,d3 ;
                        // ** * *  || * * * *
                        for(d=0;d<4;d++) {
                            d0=d1=d2=d3=0 ;
                            switch(d) {
                                case 0 : d0=1; nbp = 4; break ; // * * * *
                                case 1 : d1=1; nbp = 3; break ; // ** * *
                                case 2 : d2=1; nbp = 2 ;break ; // ** **
                                case 3 : d3=1; nbp = 2 ;break ; // *** *
                            }
                            for(i0=0;i0<=nbp && i0 <= nb0 ;i0++){
                                int64_t delta0 =  Cnp(nb0,i0) ;
                                if(d1 ) delta0 *= 3 ;
                                else if(d3) delta0 *= 2 ;
                                for(i1=0; i1<=nbp-i0 && i1<=n1;i1++){
                                    int64_t delta1 = delta0 * Cnp(n1,i1) ;
                                    int i2max = nbp-i0-i1 ;
                                    if(i2max==0) { SetIndex5( nr,n0+i1+d0,n1-i1+d1,n2+d2,n3+d3,ant,delta1) ; break ;}
                                    for(i22=0; i22*2 <= i2max && i22<=n2;i22++){
                                        for(i21=0;i21<=(i2max-2*i22) && i21+i22 <= n2;i21++){
                                            int64_t delta2 = delta1 * Cnp(n2,i22) * Cnp(n2-i22,i21)*(1<<i21);
                                            int i3max = i2max-2*i22 - i21 ;
                                            if(i3max==0) {
                                                SetIndex5( nr,n0+d0+i1+i22,n1+d1-i1+i21,n2+d2-i22-i21,n3+d3,ant,delta2) ;
                                                break ;
                                            }
                                            for(i32=0; i32*2<=i3max && i32 <= n3;i32++){
                                                for(i31=0; i31 <= i3max-2*i32 && i31+i32 <=n3;i31++){
                                                    int64_t delta3 = delta2 * Cnp(n3,i32) * Cnp(n3-i32,i31)*(1<<i31) ;
                                                    int i4max = i3max -2 * i32 - i31 ;
                                                    if(i4max == 0 ) { SetIndex5( nr,n0+d0+i1+i22+i32,n1+d1-i1+i21+i31,n2+d2-i22-i21,n3+d3-i32-i31,ant, delta3) ;break ; }
                                                    for(i43=0;i43*3<= i4max && i43 <=n4 ;i43++){
                                                        for(i42=0;i42*2<=(i4max-3*i43) && i43+i42<=n4;i42++){
                                                            int i41 = i4max-3*i43 -2*i42 ;
                                                            if(i41+i42+i43 >n4) continue ;
                                                            int64_t delta = delta3 ;
                                                            if(i43) delta *= n4 ;
                                                            if(i42) delta *= Cnp(n4-i43,i42)*pow3[i42]  ;
                                                            int j3 ;
                                                            for(j3=0;j3<=i41;j3++) {
                                                                SetIndex5( nr,n0+d0+i1+i22+i32+i43 ,n1+d1-i1+i21+i31+i42  ,n2+d2-i22-i21+i41-j3 ,n3+d3-i32-i31+j3,
                                                                          ant , delta * Cnp(n4-i43-i42,i41) * (1<<j3) * Cnp(i41,j3));
                                                            }
                                                            delta *= Cnp(n4-i43-i42,i41)*pow3[i41]  ;
                                                        } // end i42
                                                    } // end i43
                                                } // end i31
                                            } // end i32
                                        } // end i21
                                    } // end i22
                                } // end i1
                            } // end i0
                        } // end nbp ** * *  or * * * *
                        
                    } // end n3
                } // end n2
            } // end n1
        } // end n0
        
        PROB687 den,prob[20];
        
        for(i=0;i<=nr;i++)prob[i] = 0 ;
        den = 0 ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) den += cur5[i] ;
        for(n0=0;n0<=nr;n0++) {
            for(n1=0;n1<=nr-n0;n1++) {
                for(n2=0;n2<=nr-n0-n1;n2++) {
                    for(n3=0;n3<=nr-n0-n1-n2;n3++) {
                        prob[n0]+= cur5[Index5(nr,n0,n1,n2,n3)] ;
                    }
                }
            }
        }
        
        fprintf(stdout,"\t PB%s 4 color, %d rank : ",pbR->ident,nr);
        for(i=0;i<=nr;i++) { printf("%.11f ",prob[i]/(double)den) ; } ;
        printf("\n");
        if(nr==MAX_NR) {
            snprintf(pbR->strRes, sizeof(pbR->strRes),"%.10f",((double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13]))/den);
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
        PROB687 * tmp = ant5 ;
        ant5 = cur5 ;
        cur5 = tmp ;
    } // end loop nr
    printf("\n*********\n") ;
    return 0;
}

#undef FLOAT687
#undef SetIndex5
#if !defined(FLOAT687)
#include "gmp_utils.h"
typedef mpz_t PROB687m ;
#define SetIndex5(NR,N0,N1,N2,N3,ant,val)  mpz_addmul_ui(cur5[Index5(NR,N0,N1,N2,N3)],ant,val)
#   define MAX_NR   13
#endif


int PB687a(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nr ;
 
    PROB687m buf0[2380],buf1[2380] ;
    PROB687m ant, *ant5 = buf0, *cur5=buf1 ;


    int i ;
    for(i=0;i<2380;i++) { mpz_init(buf0[i]) ; mpz_init(buf1[i]) ; }
    mpz_init(ant);
    mpz_set_ui(ant5[Index5(1,0,0,0,0)],1);

    int n0,n1,n2,n3;
    // n0 * * * *
    // n1 ** * *
    // n2 ** **
    // n3 *** *
    // n4 ****
    n1 = 0 ;
    
    for(nr=2;nr<=MAX_NR;nr++) {
        int i ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) mpz_set_ui(cur5[i], 0) ;
        // pour rajouter un rang : cela peut être : *,*,* || **,* || ***
        for(n0=0;n0<nr;n0++) {
            for(n1=0;n1<nr-n0;n1++) {
                for(n2=0;n2<nr-n0-n1;n2++) {
                    for(n3=0;n3<nr-n0-n1-n2;n3++) {
                         int n4=nr-1-n0-n1-n2-n3  ;
                        //                          int nb0 = 4*n0+3*n1+2*n2+2*n3+n4+1 ;//
                        int nb0 = 4*n0+3*n1+2*n2+2*n3+n4 ;
                        mpz_set(ant,ant5[Index5(nr-1,n0,n1,n2,n3)]) ;
                        if(ant <= 0) continue ;
                           // insert ****
                        INS4_N0(ant) ; INS4_N1(ant) ; INS4_N2(ant) ; INS4_N3(ant) ; INS4_N4_2(ant) ; INS4_N4_3(ant) ;
                        
                        int pow3[4] = { 1,3,9,27} ;
                        int i0,i1,i21,i22,i31,i32,i42,i43,nbp ;
                        int d,d0,d1,d2,d3 ;
                        // ** * *  || * * * *
                        for(d=0;d<4;d++) {
                            d0=d1=d2=d3=0 ;
                            switch(d) {
                                case 0 : d0=1; nbp = 4; break ; // * * * *
                                case 1 : d1=1; nbp = 3; break ; // ** * *
                                case 2 : d2=1; nbp = 2 ;break ; // ** **
                                case 3 : d3=1; nbp = 2 ;break ; // *** *
                            }
                            for(i0=0;i0<=nbp && i0 <= nb0 ;i0++){
                                int64_t delta0 =  Cnp(nb0,i0) ;
                                if(d1 ) delta0 *= 3 ;
                                else if(d3) delta0 *= 2 ;
                                for(i1=0; i1<=nbp-i0 && i1<=n1;i1++){
                                    int64_t delta1 = delta0 * Cnp(n1,i1) ;
                                    int i2max = nbp-i0-i1 ;
                                    if(i2max==0) { SetIndex5( nr,n0+i1+d0,n1-i1+d1,n2+d2,n3+d3,ant,delta1) ; break ;}
                                    for(i22=0; i22*2 <= i2max && i22<=n2;i22++){
                                        for(i21=0;i21<=(i2max-2*i22) && i21+i22 <= n2;i21++){
                                            int64_t delta2 = delta1 * Cnp(n2,i22) * Cnp(n2-i22,i21)*(1<<i21);
                                            int i3max = i2max-2*i22 - i21 ;
                                            if(i3max==0) {
                                                SetIndex5( nr,n0+d0+i1+i22,n1+d1-i1+i21,n2+d2-i22-i21,n3+d3,ant,delta2) ;
                                                break ;
                                            }
                                            for(i32=0; i32*2<=i3max && i32 <= n3;i32++){
                                                for(i31=0; i31 <= i3max-2*i32 && i31+i32 <=n3;i31++){
                                                    int64_t delta3 = delta2 * Cnp(n3,i32) * Cnp(n3-i32,i31)*(1<<i31) ;
                                                    int i4max = i3max -2 * i32 - i31 ;
                                                    if(i4max == 0 ) { SetIndex5( nr,n0+d0+i1+i22+i32,n1+d1-i1+i21+i31,n2+d2-i22-i21,n3+d3-i32-i31,ant, delta3) ;break ; }
                                                    for(i43=0;i43*3<= (i4max) && i43 <=n4 ;i43++){
                                                        for(i42=0;i42*2<=(i4max-3*i43) && i43+i42<=n4;i42++){
                                                            int i41 = i4max-3*i43 -2*i42 ;
                                                            if(i41+i42+i43 >n4) continue ;
                                                            int64_t delta = delta3 ;
                                                            if(i43) delta *= n4 ;
                                                            if(i42) delta *= Cnp(n4-i43,i42)*pow3[i42]  ;
                                                            int j3 ;
                                                            for(j3=0;j3<=i41;j3++) {
                                                                SetIndex5( nr,n0+d0+i1+i22+i32+i43 ,n1+d1-i1+i21+i31+i42  ,n2+d2-i22-i21+i41-j3 ,n3+d3-i32-i31+j3,
                                                                          ant , delta * Cnp(n4-i43-i42,i41) * (1<<j3) * Cnp(i41,j3));
                                                            }
                                                            delta *= Cnp(n4-i43-i42,i41)*pow3[i41]  ;
                                                        } // end i42
                                                    } // end i43
                                                } // end i31
                                            } // end i32
                                        } // end i21
                                    } // end i22
                                } // end i1
                            } // end i0
                        } // end nbp ** * *  or * * * *
  
                    } // end n3
                } // end n2
            } // end n1
        } // end n0
        
        PROB687m den,denavg,prob[20];

        PROB687m g,q,p100,avg ;
        mpz_init (den) ; mpz_init (denavg) ;mpz_init (avg) ;mpz_init(g); mpz_init(q);mpz_init(p100);
        for(i=0;i<=nr;i++)mpz_init(prob[i]) ;
        mpz_set(g,cur5[0]) ;
        for(i=1;i<=Index5(nr,0,0,0,0);i++) {  mpz_gcd (g, g,cur5[i]); if(mpz_cmp_ui(g,1)==0) break ; } ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) {  mpz_divexact(cur5[i],cur5[i],g) ; mpz_add(den,den,cur5[i]); }
        for(n0=0;n0<=nr;n0++) {
            for(n1=0;n1<=nr-n0;n1++) {
                for(n2=0;n2<=nr-n0-n1;n2++) {
                    for(n3=0;n3<=nr-n0-n1-n2;n3++) {
                        mpz_add(prob[n0],prob[n0],cur5[Index5(nr,n0,n1,n2,n3)]) ;
                    }
                }
            }
        }
        
        for(i=0;i<=nr;i++){  mpz_addmul_ui (avg,prob[i],i); } ;
        mpz_gcd (g, avg,den); mpz_divexact(avg,avg,g) ; mpz_divexact(denavg,den,g) ;
        gmp_printf("\t PB%s 4 color, %d rank, avg=%Zd/%Zd :",pbR->ident,nr,avg,denavg);
        mpz_set(g,prob[0]) ;
        for(i=0;i<=nr;i++){  mpz_gcd (g, g,prob[i]); if(mpz_cmp_ui(g,1)==0) break ; } ;
        for(i=0;i<=nr;i++) { mpz_divexact(prob[i],prob[i],g) ;
            mpz_mul_ui(p100,prob[i],10000000000LL) ; mpz_tdiv_q(q,p100,den);
            printf("%.10f ",mpz_get_ui(q)/10000000000.0) ;
        } ;
       printf("\n");
        if(nr==MAX_NR) {
            PROB687m sum  ;
            mpz_init(sum); mpz_init(avg); mpz_add(sum,prob[2],prob[3]) ; mpz_add(sum,sum,prob[5]) ; mpz_add(sum,sum,prob[7]) ;
            mpz_add(sum,sum,prob[11]) ;mpz_add(sum,sum,prob[13]) ;
            mpz_gcd (g, sum,den); mpz_divexact(sum,sum,g) ; mpz_divexact(den,den,g) ;
             gmp_printf("\t PB%s Prob(Prime)= %Zd / %Zd\n",pbR->ident, sum,den);
            mpz_mul_ui(sum,sum,10000000000LL) ;  mpz_tdiv_q(q,sum,den);
            snprintf(pbR->strRes, sizeof(pbR->strRes),"%.10f",mpz_get_ui(q)/10000000000.0);
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
        PROB687m * tmp = ant5 ;
        ant5 = cur5 ;
        cur5 = tmp ;
    } // end loop nr
    return 0;
}

#define PB690_NBV   10
int64_t MULT(int k, int n) {
    int i ;
    int64_t M = 1 ;
    for(i=1;i<k;i++) {
       M = (M * (n+i))/i ;
    }
    return M ;
}

int PB690(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    int n = PB690_NBV ;
    int64_t *nbT = malloc((n+1)*sizeof(nbT[0])) ;
    nbT[1] = 1; nbT[2] = 1 ; nbT[3] = 1 ;
    
    for(i=4;i<=n;i++) {
        nbT[i] = ( 1 << (i-4) )  +  ( 1 << ((i-4)/2) )
        + ((i< 7) ? 0 : 1 );
        
    }
    for(i=1;i<=n;i++) printf("%lld ",nbT[i]) ;
    printf("\n");
    int64_t N = MULT(nbT[7],1) ;
    Decomp  * Dec  = DecompAlloc(n) ;
    while(DecompNext(Dec)>=0) {
        int j ;
        int nb = 1 ;
        int64_t P = 1 ;
        for(j=1;j<Dec->nbVal;j++) {
            if(Dec->val[j] != Dec->val[j-1]) {
                P *= MULT(nbT[Dec->val[j-1]],nb) ;
                nb = 1 ;
            } else {
                nb++ ;
            }
        }
        P *= MULT(nbT[Dec->val[j-1]],nb);
        N += P ;
        for(j=0;j<Dec->nbVal;j++) printf("%d ",Dec->val[j]);
        printf(" => %lld\n",P);
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",N);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

//#define PB691_LG    5000000
#define PB691_LG    5000000

uint8_t * InitSuit(int len) {
    uint8_t *suit = malloc(PB691_LG*sizeof(suit[0])) ;
    int64_t Fn = 1, Fd = 1; // approximation of 1/phi by rationnal fraction with fibonacci
    int i ;
    for(i=0;i<38;i++) {    int64_t tmp = Fn ; Fn += Fd ; Fd = tmp ; }
    suit[0] = 0 ;
    int is = 1;
    for(i=1;is*2 < PB691_LG;is *=2) {
        for(i=is;i<2*is;i++) suit[i] = 1 - suit[i-is] ;
    }
    for(i=is;i<PB691_LG;i++)suit[i] = 1 - suit[i-is] ;
    int64_t antPhi = 0 ;
    for(i=0;i<PB691_LG;i++) {
        antPhi += Fd ;
        if( antPhi >= Fn) { suit[i] = 1 - suit[i] ;    antPhi -= Fn ; }
    }
    return suit ;
}

DC3_INDEX * GetMaxDupByLength(int * ptMaxLcp,const DC3_INDEX *lcp,int n) {
    int maxLcp = 0 ;
    int i ;
    DC3_TYPE  *B, *A;// before and after
    B= malloc(n*sizeof(B[0])); A= malloc(n*sizeof(A[0]));
    B[0] = -1;
    for (i = 1; i < n; ++i) { // search suffixes before with same common part
        int li = lcp[i];
        if(li > maxLcp) maxLcp = li ;
        int id = i - 1 ;
        while (id != -1 && lcp[id] >= li) { ; id = B[id]; }
        B[i] = id;
    }
    A[n - 1] = n;
    for (i = n - 2; i >= 0; --i) { // search suffixes after with same common part
        int li = lcp[i];
        int id = i + 1;
        while (id != n && lcp[id] >= li) id = A[id];
        A[i] = id;
    }
    DC3_INDEX *maxDupBylengh = calloc(maxLcp+1,sizeof(maxDupBylengh[0])) ;
    for (int i = 0; i < n; ++i) {
        int li = lcp[i];
        int cnt = A[i] - B[i] - 1 ; //  ] A[i] ... B[i] [
        if(cnt > maxDupBylengh[li]) {
            maxDupBylengh[li] = cnt ;
        }
    }
    *ptMaxLcp = maxLcp ;
    free(A); free(B) ;
    return maxDupBylengh ;
}

int PB691a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    const int n = PB691_LG ;
    uint8_t * suit = InitSuit(PB691_LG) ;
    DC3_TYPE *suit_dc3 = malloc((PB691_LG)*sizeof(suit_dc3[0])) ;
    for(i=0;i<n;i++) suit_dc3[i] = suit[i]  ;
    DC3_INDEX * sa= malloc(n*sizeof(sa[0]));
    SuffixSort(suit_dc3, sa, n, 2);
//    for(i=0;i<n;i++) { for(int j=sa[i];j<n;j++) printf("%d ",suit[j]); printf("\n");}
    free(suit_dc3 ); // not necessary after
    DC3_INDEX * lcp = GetLCP(suit, sa, n,n);
    free(sa); free(suit) ; // not necessary after
    DC3_INDEX maxLcp  = 0 ;
    DC3_INDEX * maxDupBylengh = GetMaxDupByLength(&maxLcp,lcp,n) ;
    free(lcp) ;
    int32_t antdup = 1 ;
    int64_t S =0 ;
    for(i=1;i<=maxLcp;i++) {
//        printf("%d %d\n",i,maxDupBylengh[i]) ;
        if(maxDupBylengh[i]) {
            antdup = maxDupBylengh[i] ;
        }
        S += antdup ;
    }
    S += n ;
    free(maxDupBylengh);
    printf("MaxLcp=%d ans=%lld\n",maxLcp,S);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PCKTYPE uint16_t
#define PCK691  15  // sizeof(PCKTYPE)-1


#define PCKSET  (1<<PCK691)



PCKTYPE *suitPck ;
// #define SORTSIMPLE
#if defined(SORTSIMPLE)
#   define OFFSET0  0
#else
#   define OFFSET0  1
#endif
int cmpSa(const void * el1, const void *el2) {
    int d12 = ((int *)el2)[0] - ((int *)el1)[0] ;
    if(d12 == 0) return 0 ;
    PCKTYPE *s1 = suitPck + ((int *)el1)[0] + OFFSET0 ;
    while( (*s1 & PCKSET) && *s1 == *(s1+d12) ) { s1++ ;  }
    if(*s1 & PCKSET) {
        if(s1[d12] & PCKSET) {
            return *s1 - s1[d12] ;
        } else {
            return 1 ;
        }
    } else if(s1[d12] & PCKSET) {
        return -1 ;
    } else {
        return *s1 - s1[d12] ;
    }
}

PCKTYPE * GetSuitPck(const uint8_t * suit,int debI0[PCK691+1],int n) {
    int lgPck = n+PCK691 ;
    PCKTYPE *  suitPck = malloc(lgPck*sizeof(suitPck[0])) ;
    int i,i0 ;
    int nbPck =0 ;
    for(i0=0;i0<PCK691;i0++) {
        int lg1 = i0+((n-i0)/PCK691)*PCK691;
        debI0[i0] = nbPck ;
        for(i=i0;i<lg1;i += PCK691) {
            PCKTYPE cur = PCKSET ; // pack PCK691 bits + bit PCKSET
            for(int ni =0; ni<PCK691;ni++) cur |= suit[i+ni] <<  (PCK691 -ni -1) ;
            suitPck[nbPck++] = cur ;
        }
        PCKTYPE last = PCKSET ; int iEnd = i ;
        if(iEnd < n) { // pack remaining bits + bit PCKSET
            for(;i<n;i++) {last |= suit[i] << (PCK691-1 - i + iEnd) ; }
            suitPck[nbPck++] = last ;
        }
        suitPck[nbPck++] = n - iEnd ; // cout of bit in last packet, no PCKSET
    }
    debI0[PCK691] = nbPck ;
    return suitPck ;
}

int PB691b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    uint8_t * suit = InitSuit(PB691_LG) ;
    const int n = PB691_LG ;
     int debI0[PCK691+1] ;
    suitPck = GetSuitPck(suit,debI0,n) ;
    int nbPck = debI0[PCK691] ;
    int32_t *saPck = malloc(n*sizeof(saPck[0]));
    for(i=0;i<n;i++) {
        int i0 = i % PCK691 ;
        saPck[i] = debI0[i0] + (i-i0)/PCK691 ;
    }
    int32_t  *orderPck = malloc(nbPck * sizeof(orderPck[0]))  ;
    for(i=0;i<n;i++)orderPck[saPck[i]]= i ;
    DC3_INDEX * sa = malloc(n*sizeof(sa[0])) ;

#if defined(SORTSIMPLE)
    heapsort(saPck,n,sizeof(saPck[0]),cmpSa) ;
    for(i=0;i<n;i++) { sa[i] = orderPck[saPck[i]] ; }
    
#else
    // one pass radixSort and
    int32_t *sao = malloc(n*sizeof(sao[0])) ;
    DC3_INDEX *  count = CountRadixSort_uint16(saPck,sao,suitPck,n, 2*PCKSET) ;
    int is =0;
    for(int ic=0;ic <= 2*PCKSET+1 ; ic++) { // loop on radix , sort by group of firms elem equal
        if(count[ic]-is > 1) {
            heapsort(sao+is,count[ic]-is,sizeof(sao[0]),cmpSa) ;
        }
        is = count[ic] ;
    }
    for(i=0;i<n;i++) { sa[i] = orderPck[sao[i]] ; }
    free(sao); free(count);
#endif
    free(saPck); free(orderPck) ; free(suitPck); suitPck = NULL ;
    DC3_INDEX * lcp = GetLCP(suit, sa, n,n);
    free(sa); free(suit) ; // not necessary after
    DC3_INDEX maxLcp  = 0 ;
    DC3_INDEX * maxDupBylengh = GetMaxDupByLength(&maxLcp,lcp,n) ;
    free(lcp) ;
    int32_t antdup = 1 ;
    int64_t S =0 ;
    for(i=1;i<=maxLcp;i++) {
         if(maxDupBylengh[i]) {
            antdup = maxDupBylengh[i] ;
        }
        S += antdup ;
    }
    S += PB691_LG ;
    free(maxDupBylengh);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

typedef struct DupPat {
    int nbDup ;
    int start ;
    int end ;
} DupPat ;

#define NBITS691    18

uint8_t * suitInt8 ;
int suitLen ;

int cmpSaUint8(const void * el1, const void *el2) {
    int i1 = ((int *)el1)[0] ;
    int i2 = ((int *)el2)[0] ;
    int d12 = i2 - i1 ;
    if(d12==0) return 0 ;
    int len = suitLen- NBITS691 -( (d12 > 0) ? i2 : i1) ;
    if(len <= 0) return d12 ;
    uint8_t *s1 = suitInt8 + NBITS691 + i1 ;
    int diff = memcmp(s1,s1+d12,len) ;
    if(diff) return diff ;
    else return d12 ;
}

int PB691(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    uint8_t * suit = InitSuit(PB691_LG) ;
    int n = PB691_LG ;
    suitLen = n ;
    suitInt8 = suit ;
    int nbHisto = (1 << NBITS691) ;
    int maskStart = nbHisto - 1 ;
    int * iStart = malloc(PB691_LG*sizeof(iStart[0]));
    int * sa = malloc(PB691_LG*sizeof(sa[0]));
// radix sort on NBITS691 bits
    int * histoStart = calloc(nbHisto,sizeof(histoStart[0]));
    int curStart = 0 ;
    int i ;
    for(i=0;i<NBITS691-1;i++) curStart = (curStart << 1) | suit[i] ;
    for(i=0;i<n-NBITS691+1;i++) {
        curStart = ((curStart << 1) & maskStart) | suit[i+NBITS691-1] ;
        iStart[i] = curStart ;
        histoStart[curStart]++ ;
    }
    for(;i<n;i++) { curStart = (curStart << 1) & maskStart  ; iStart[i] = curStart ; histoStart[curStart]++ ; }
    int sumH = histoStart[0] ;
    histoStart[0] = 0 ;
    int ih ;
    for(ih=1;ih<nbHisto;ih++) {
        int tmp = histoStart[ih] ;
        histoStart[ih] = sumH ;
        sumH += tmp ;
    }
    for(i=0;i<n;i++) {
        sa[histoStart[iStart[i]]++] =  i ;
    }
    int nb_sa = 0;
    for(ih=0;ih<nbHisto;ih++) {
        if(histoStart[ih]-nb_sa > 1) {
            qsort(sa+nb_sa,histoStart[ih]-nb_sa,sizeof(sa[0]),cmpSaUint8);
            nb_sa = histoStart[ih] ;
        } else if(histoStart[ih]-nb_sa == 1) {
            nb_sa++ ;
        }
    }
    free (histoStart) ; free(iStart);
//    printf("nb_sa=%d\n",nb_sa);
//    for(i=0;i<nb_sa;i++) {printf("\n-"); for(int j=sa[i];j<n;j++) printf("%d",suit[j]); }

    DC3_INDEX * lcp = GetLCP(suit, sa, nb_sa,n);
    free(sa); free(suit) ; // not necessary after
    DC3_INDEX maxLcp  = 0 ;
    DC3_INDEX * maxDupBylengh = GetMaxDupByLength(&maxLcp,lcp,nb_sa) ;
    free(lcp) ;
    int antdup = 1 ;
    int64_t S =0 ;
    for(i=0;i<=maxLcp;i++) {
        if(maxDupBylengh[i]) {
            antdup = maxDupBylengh[i] ;
        }
        S += antdup ;
    }
    free(maxDupBylengh);
     snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
     pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB691c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    uint8_t * suit = InitSuit(PB691_LG+1) ;
    int n = PB691_LG ;
//    n = 50 ;
    for(int i=0;i<n;i++) suit[i] += 1 ;
    suit[n] = 0 ;
    STree *ST = X_CompSuffixTree(suit,3,n);
    for(i=0;i<n;i++) suit[i] -= 1 ;

    int *sa = NULL ; // not necessary ,compute directly lcp
    int *lcp =(int*) malloc(n*sizeof(lcp[0]));
    CompSALCP(ST,sa,lcp)  ;
    FreeSuffixTree(ST); free(sa); free(suit) ; // not necessary after
    DC3_INDEX maxLcp  = 0 ;
    DC3_INDEX * maxDupBylengh = GetMaxDupByLength(&maxLcp,lcp,n) ;
    free(lcp) ;
    int antdup = 1 ;
    int64_t S =0 ;
    for(i=0;i<=maxLcp;i++) {
        if(maxDupBylengh[i]) {
            antdup = maxDupBylengh[i] ;
        }
        S += antdup ;
    }
    free(maxDupBylengh);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


int PB691d(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    uint8_t * suit = InitSuit(PB691_LG) ;
    int n = PB691_LG ;
    STree *ST = ComputeSuffixTree(suit,2,n);
    int *sa = NULL ; // not necessary ,compute directly lcp
    int *lcp =(int*) malloc(n*sizeof(lcp[0]));
    CompSALCP(ST,sa,lcp)  ;
    
    printf("Nbnode=%d\n",NbNodeUsed(ST)); 
    FreeSuffixTree(ST);  free(sa); free(suit) ; // not necessary after
    DC3_INDEX maxLcp  = 0 ;
    DC3_INDEX * maxDupBylengh = GetMaxDupByLength(&maxLcp,lcp,n) ;
    free(lcp) ;
    int antdup = 1 ;
    int64_t S =0 ;
    for(i=0;i<=maxLcp;i++) {
        if(maxDupBylengh[i]) {
            antdup = maxDupBylengh[i] ;
        }
        S += antdup ;
    }
    free(maxDupBylengh);
    
    
//    freeSuffixTreeByPostOrder(root);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB693_MAX  3000000

int Calc693(int32_t x, int32_t *curModx, int32_t *nextModx,int32_t *curMax,int32_t *nextMax) {
    int isModeSparse = 0 ;
    int xStart = x ;
    int modx,y,nbY ;
    memset(curMax,0,((x+1)/2+1)*sizeof(curMax[0])) ;
    // initialise squares
    for(y=2,modx=1;2*y<=x+1;y++) {
        modx += 2*y-1 ; if(modx>=x) { modx -= x ;  }
        int modx2 = x+1 -modx ;
        if(modx < modx2) modx2 = modx ;
        if(modx2>1) {
            curMax[modx2] = 1 ;
        }
    }
    for(++x;;x++) { // loop on x until no y alive.
        if(isModeSparse == 0) { // many y alive, compute all vlaues
            memset(nextMax,0,((x+1)/2+1)*sizeof(nextMax[0])) ;
            nbY = 0 ;
            for(y=2,modx=1;2*y<=x+1;y++) {
                modx += 2*y-1 ; if(modx>=x) { modx -= x ;}
                int curMy= curMax[y]+1 ;
                if(curMy == 1) continue ; // y not alive
                int modx2 = x+1 -modx ;
                if(modx < modx2) modx2 = modx ; // restrains modx to [0 x/2]
                if (modx2 >1 && curMy >nextMax[modx2]) {
                    nbY++ ;
                    nextMax[modx2] = curMy ;
                }
            }
            if(nbY < x/8)
            {
                // swithc to sparse mode
                int iy = 0 ;
                isModeSparse = 1 ;
                for(y=2,modx=1;2*y<=x+1;y++) {
                    if(nextMax[y]) { // search alive y
                        curModx[iy++] = y ;
                     }
                }
                nbY = iy ;
                 memset(curMax,0,((x+1)/2+1)*sizeof(curMax[0])) ;
                 nextMax[(x+1)/2+1] = 0 ;
             }
          } else { // sparse mode
            if(nbY > 1 ) {
                nextMax[(x+1)/2+1] = curMax[(x+1)/2+1] = 0 ;
                int32_t modx ;
                int iy ;
                int nextNbY = 0 ;
                for(iy=0;iy<nbY;iy++) {
                    int32_t y =curModx[iy] ;
                    int yMax = curMax[y]+1  ;
                    curMax[y] = 0 ;
                    modx = ((int64_t)y * y ) % x ;
                    if(x+1 - modx < modx ) modx = x+1 - modx ;
                    if(modx > 1) {
                        if(nextMax[modx] == 0) {
                            nextMax[modx] = yMax ;
                            nextModx[nextNbY++] = modx ;
                        } else if(yMax > nextMax[modx]) {
                            nextMax[modx] = yMax ;
                        }
                    }
                }
                nbY = nextNbY ;
                if(nbY == 0) {  break ; }
                int32_t *tmp = nextModx ;
                nextModx = curModx ;
                curModx = tmp ;
            } else if (nbY==1) { // only one y alive
                modx = curModx[0] ;
                int yMax = curMax[modx]+1 ;
                while(1) {
                    modx = ((int64_t) modx * modx) % x ;
                    if(modx <= 1 ) break ;
                    yMax++ ; x++ ;
                }
                break ;
            }
        }
        int32_t * tmp = nextMax;
        nextMax = curMax;
        curMax = tmp ;
    }
    return x-xStart+1 ;
}
//
// search f(x) where we init values only for one y=p and when x=0 mod(freq)
int Calc693OneValue(int32_t p,int32_t n,int32_t freq, int32_t *curModx, int32_t *nextModx,int32_t *curMax,int32_t *nextMax) {
    int xStart = p+1 ;
    int p2 = p*p ;
    int x = xStart ;
    int lgMax = 0 ;
    memset(nextMax,0,((x+1)/2+1)*sizeof(nextMax[0])) ;
    memset(curMax,0,((x+1)/2+1)*sizeof(curMax[0])) ;
    int nbY = 0 ;
    int ifreq = 0 ;
    for(++x;;x++) {
        nextMax[(x+1)/2+1] =curMax[(x+1)/2+1] = 0 ;
        if( ifreq-- == 0 && x <= n ) {
            ifreq = freq-1 ;
            int modx = p2 % x ;
            if(x+1 - modx < modx ) modx = x+1 - modx ;
            if(modx > 1 && curMax[modx] == 0) {
                curMax[modx] = 1 ;
                curModx[nbY++] = modx ;
            }
        }
        int iy ;
        int nextNbY = 0 ;
        if(x>=n && nbY == 0 ) break ;
        for(iy=0;iy<nbY;iy++) {
            int32_t y =curModx[iy] ;
            int32_t yMax = curMax[y]+1  ;
            curMax[y] = 0 ;
            int32_t modx = ((int64_t)y * y ) % x ;
            if(x+1 - modx < modx ) modx = x+1 - modx ;
            if(modx > 1) {
                if(nextMax[modx] == 0) {
                    nextMax[modx] = yMax ;
                    nextModx[nextNbY++] = modx ;
                } else if(yMax > nextMax[modx]) {
                    nextMax[modx] = yMax ;
                }
            } else if(yMax > lgMax ) {
                    lgMax = yMax ;
            }
        }
        nbY = nextNbY ;
        int32_t *tmp = nextModx ;
        nextModx = curModx ;
        curModx = tmp ;
        tmp = nextMax;
        nextMax = curMax;
        curMax = tmp ;
    }
    return lgMax+1 ;
}


// #define PB693_TRACE
int PB693(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB693_MAX ;
    int32_t *tb1 = malloc(2*n*sizeof(tb1[0]));
    int32_t *tbModx1 = malloc(2*n*sizeof(tbModx1[0]));
    int32_t *curModx = tbModx1 ;
    int32_t *nextModx = tbModx1+n ;
    int32_t *curMax = tb1 ;
    int32_t *nextMax = tb1+n ;
    int x;
    int lgMax = Calc693OneValue(2,n,n/40,curModx,nextModx,curMax,nextMax);
    for(x=n;x>5;) {
        int32_t delta = Calc693(x,curModx,nextModx,curMax,nextMax) ;
#if !defined(PB693_TRACE)
        printf("\n\t%.4fs [%d %d]->%d ",(float)(clock() - pbR->nbClock)/ CLOCKS_PER_SEC,x,x+lgMax,
               x+delta-lgMax+1);
#endif
        if(delta > lgMax) {
            lgMax = delta ;
          }
        if(delta >= lgMax-1) {
            int xHigh = x ;
            int xLow = x - lgMax ;
            if(xLow < 5) xLow = 5 ;
            int end0 = x+delta ;
            int deltaLow = delta ;
            while(xHigh-xLow>=2) {
               int xMed = (xHigh+xLow)/2 ;
               int deltaMed = Calc693(xMed,curModx,nextModx,curMax,nextMax) ;
                if(deltaMed > lgMax) lgMax = deltaMed ;
                if(xMed+deltaMed >= end0) {
                    xHigh = xMed ;
                } else {
                    xLow = xMed ; deltaLow = deltaMed ;
                }
            }
            x = xLow - lgMax + deltaLow    ;
#if !defined(PB693_TRACE)
            printf("->[%d %d]",x,xHigh);
#else
            printf("\nH %d %d\nL %d %d\nX %d %d",xHigh,lgMax+1,xLow,deltaLow+1,x,lgMax+1) ;
#endif
        } else {
#if defined(PB693_TRACE)
            printf("\nH %d %d\nL %d %d",x,delta+1,x-lgMax+delta,lgMax+1) ;
#endif
            x -= lgMax-delta-1 ;
        }
    }
    printf("\nFor n=%d lgMax=%d \n",n,lgMax+1);
    free(tb1); free(tbModx1) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",lgMax+1);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

// version to compute detailed f(x)
int PB693b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB693_MAX ;
    int32_t *tb1 = malloc(2*n*sizeof(tb1[0]));
    int32_t *tbModx1 = malloc(2*n*sizeof(tbModx1[0]));
    int32_t *curModx = tbModx1 ;
    int32_t *nextModx = tbModx1+n ;
    int32_t *curMax = tb1 ;
    int32_t *nextMax = tb1+n ;
    int x;
 /*   int d=0,d1=0,d2=0 ;
    for(x=5;x<=n;x++) {
        d2=d1; d1= d ;
        d = Calc693(x,curModx,nextModx,curMax,nextMax);
        if(d1<d2 && d1 < d) printf("v%d(%d) ",x-1,d1+1) ;
        else if(d1>d2 && d1 >d)printf("^%d(%d) \n",x-1,d1+1) ;
//        else printf("%d->%d ",x-1,d1+1) ;
    }
  */
    int lgMax = 1 ;
    int deltaAnt = 1 ;
    for(x=n;x>5;) {
        int32_t delta = Calc693(x,curModx,nextModx,curMax,nextMax) ;
        printf("\n\t%.4fs [%d %d] ",(float)(clock() - pbR->nbClock)/ CLOCKS_PER_SEC,x,x+delta);
        if(delta >= deltaAnt-1) {
            int xHigh = x ;
            int deltaHigh = delta ;
            int endH = x+deltaHigh ;
            int xLow = x - lgMax ;
            if(xLow < 5) xLow = 5 ;
            int deltaLow =  Calc693(xLow,curModx,nextModx,curMax,nextMax);
            printf("Dich[%d %d] ",xLow,xHigh);
            while(xHigh-xLow>=2) {
                int xMed = (xHigh+xLow)/2 ;
                int deltaMed = Calc693(xMed,curModx,nextModx,curMax,nextMax) ;
                 if(xMed+deltaMed >= endH) {
                    xHigh = xMed ;  deltaHigh = deltaMed ;
               } else {
                    xLow = xMed ; deltaLow = deltaMed ;
                }
            }
            x = xLow-1 ;
            if(deltaHigh> lgMax)lgMax=deltaHigh ;

            printf("->[%d %d]",xLow,xHigh);

            printf("\nH %d %d \nL %d %d",xLow+1,deltaHigh+1,xLow,deltaLow+1);

        } else {
            x-- ;
        }
        deltaAnt = delta ;
    }
    free(tb1); free(tbModx1) ;
    printf("\nFor n=%d lgMax=%d \n",n,lgMax+1);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",lgMax+1);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int PB693a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    int n = PB693_MAX ;
    int32_t *tb1 = malloc(4*n*sizeof(tb1[0]));
    int32_t *tb2 = malloc(4*n*sizeof(tb2[0]));
    
    int32_t *curMax = tb1 ;
    int32_t *nextMax = tb2 ;
    int x,y;
    int lgMax = 0;
    curMax[2] = 1 ;
    memset(nextMax,0,3*sizeof(nextMax[0])) ;
    //   curMax[3] = 1 ;
    //   curMax[4] = 0 ;
    for(x=5;x<=n;x++) {
        int modx ;
        memset(nextMax,0,((x+1)/2+1)*sizeof(nextMax[0])) ;
        
        for(y=2,modx=1;2*y<=x;y++) {
            modx += 2*y-1 ; if(modx>=x) { modx -= x ; /*if(modx>=x) modx -= x ;*/ }
            int modx2 = x+1 -modx ;
            if(modx < modx2) modx2 = modx ;
            int curMy  ;
            if((curMy=curMax[y]+1) == 1) continue ;
            if(modx2<=1) {
                if(curMy>lgMax) {
                    lgMax = curMy ; printf("M=%d->%d ",lgMax,x) ;
                }
            } else if (curMy >nextMax[modx2]) {
                nextMax[modx2] = curMy ;
            }
        }
        // printf("\n");
        for(y=2,modx=1;2*y<=x+1;y++) {
            modx += 2*y-1 ; if(modx>=x) { modx -= x ; /*if(modx>=x) modx -= x ;*/ }
            int modx2 = x+1 -modx ;
            if(modx < modx2) modx2 = modx ;
            if(modx2>1 && nextMax[modx2] == 0) {
                nextMax[modx2] = 1 ;
            }
        }
        //     curMax[y] = 0 ;
        int32_t * tmp = nextMax;
        nextMax = curMax;
        curMax = tmp ;
    }
    int isChange  ;
    do {
        isChange = 0 ;
        memset(nextMax,0,((x+1)/2+1)*sizeof(nextMax[0])) ;
        int modx ;
        for(y=2,modx=1;2*y<=x;y++) {
            modx += 2*y-1 ; if(modx>=x) { modx -= x ; if(modx>=x) modx -= x ; }
            int modx2 = x+1 -modx ;
            if(modx < modx2) modx2 = modx ;
            int curMy ;
            if((curMy=curMax[y]+1) == 1) continue ;
            if(modx2<=1) {
                if(curMy>lgMax) { lgMax = curMy ;  }
            } else if (curMy >nextMax[modx2]) {
                isChange = 1 ;
                nextMax[modx2] = curMy ;
            }
        }
        int32_t * tmp = nextMax;
        nextMax = curMax;
        curMax = tmp ;
        
        
        x++ ;
    } while (isChange );
    printf("For n=%d lgMax=%d \n",n,lgMax+1);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",lgMax+1);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


//**************************************

/*
 #define PB696_NBS   3
#define PB696_NBT   4
#define PB696_NBN   30
*/





#define PB696_NBS   100000000
#define PB696_NBT   30
#define PB696_NBN   100000000


#define PB696_MOD   1000000007

#define isSerie 1
#define isT3    2
#define isPaire 4
#define NB_T696a 9 // 44 33 34 22 23 24 12 13 02
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

int PB696(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB696_NBN ;
    int ntMax1 = PB696_NBT+1 ;
    int ns = PB696_NBS ;
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
        int isSerieAllowed = ( in<n-2) ? 1 : 0 ; //  for small values for n
        int32_t *tmp = tbCur ; tbCur = tbAnt ;  tbAnt = tmp ;
        for(int ib=0;ib<maxIb;ib++){ // loop on breaks (frontiers between connected dev.)
            int64_t factCnp = Cnp[INDXIBIN(ib,in)] ;
            int ibCur = ib * deltaBreak ;
            for(int nty=0;nty<NB_T696a;nty++) {
                int nb2 = nty2nb2[nty] ;
                int nb1 = nty2nb1[nty]   ;
                for(int st=0;st<NB_STP696a;st++){ // states to avoid multiple representations
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
                                    nb_IT[it].noP += na * factCnp ; nb_IT[it].noP %= PB696_MOD ;
                                }  else  {
                                    nb_IT[it].withP += na * factCnp ; nb_IT[it].withP %= PB696_MOD ;
                                }
                            } else {   continue ; }
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
                    }
                } // end bcl it
            } // end bcl nty
        } // end bcl ib
    } // end bcl in
    free(tbCur); tbCur = NULL ;
    free(tbAnt); tbAnt = NULL ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s End phase 0 : compact hands %.6fs\n",pbR->ident,(float)(clock()-pbR->nbClock) / CLOCKS_PER_SEC);
    
    
    
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



#define PB697_N 10000000

// we must solve
// f(x) = Sum[x**i/i!] (i [0..n]) / exp(x) - 1/4
// f(x) = 0
// f'(x) = - x**n / n! / exp(x).
// Solve by newton
// x(k+1) = x(k) - f[x(k)]/f'[(x(k)] = x(k) + delta(k)
// delta(k) = Sum[ n!/ i! / x(k)**(n-i) ] - 1/4 x n! / x(k)**n x exp(x(k))
// Stirling approximation for n! ~ n**n / exp(n) * sqrt(2xpixn) * (1+1/12*n)
int PB697(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB697_N-1 ;
    double s = n ;
    double delta ;
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s",pbR->ident) ;
    do {
        double fact = 1 ;
        delta = 1 ;
        int i;
        double s1 = 1/ s ;
        for(i=n;i*fact > 0.01;i--){
            fact *=  i * s1 ;
            delta += fact ;
        }
         delta -=  exp (n * (log(n) - 1  - log(s)) + s + log(n)/2 + log(M_PI*2)/2)* (1+1/12.0/n) * 0.25 ;
        s += delta ;
        if(pbR->isVerbose) fprintf(stdout," +%.4f=>%.5f",delta,s) ;

    } while (delta > 0.01 || delta < -0.01) ;
    if(pbR->isVerbose) fprintf(stdout,"\n") ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.2f",s / log(10.0));
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

// by diccotomie
//
int PB697a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int n = PB697_N ;
    double vLow, vHigh ;
    vLow = n ;
    vHigh = n*2.0 ;
    while(vHigh - vLow > 0.01 ) {
        double vMed = (vHigh+vLow)/2 ;
        double exp_vMed_n = exp(-vMed /n) ;
        int in = 0 ;
        int i0 = n - 20000 ;
        double sumF, fact = 1;
        if(i0 < 0) {
            i0 = 0 ;
            sumF = fact  ;
        } else {
            sumF = 0 ;
        }
        for(int i=i0+1;i<n;i++){
            fact *= vMed / i ;
            while (sumF > 1.0) {
                sumF *= exp_vMed_n ;
                fact *= exp_vMed_n ;
                in ++ ;
            }
            sumF += fact ;
        }
        if(i0) {
            double lnFactI0 = (log(vMed) - log(i0)+1) * i0 ;
            int in0 = lnFactI0 / ( vMed/(n)) ;
            in += in0 ;
            sumF *= exp( i0 * ( log(vMed) - log(i0) + 1 ) -in0 * (vMed/n) - log(i0) /2 -log(M_PI*2)/2 ) ;
        }
        sumF /= exp(vMed*(n-in)/n) ;
        if(sumF < .25) {
            vHigh = vMed ;
        } else {
            vLow = vMed ;
        }
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\tPB%s Ln=%.4f ",pbR->ident,(vHigh+vLow)/2) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.2f",(vHigh+vLow) / log(10.0)/2);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB701_L 7
#define PB701_H 7
#define PB701_MAXCONNECT    5




typedef struct PB701_PP {
    int num ;
    int globalContact ;
    int nbSet ;
    uint16_t   Set[PB701_MAXCONNECT]  ;
} PB701_PP ;

typedef struct PB701_CONN {
    int globalContact ;
    int nbSet ;
    uint16_t Set[PB701_MAXCONNECT]  ;
    uint8_t sizeC[PB701_MAXCONNECT] ;
}PB701_CONN ;

typedef struct PB701_TRANSF {
    int     numPP ;
    uint8_t numCv[PB701_MAXCONNECT] ;
    uint8_t delta[PB701_MAXCONNECT] ;
} PB701_TRANSF ;

typedef struct PB701_INTER {
    uint8_t surfIncrement ;
    uint16_t newContact ;
} PB701_INTER ;


#define PB701_BIG
#if defined(PB701_BIG)
typedef __int128 Prob701 ;
#else
typedef uint64_t Prob701 ;
#endif

typedef struct PB701_STATE {
    uint16_t numPP ;    // PP corresponding
    NsumVal  Surf[PB701_MAXCONNECT+1] ; // surface for each set (including hidden set)
    Prob701  prob ;     // weight of the state
} PB701_STATE ;





typedef struct PB701_DATA {
    int nbCol ;   // number of differents column ( 1<<sizeColumn)
    int nbPP ;  //  number of column partitions in connected parts
    int *isSymetric  ; // 0 <=> mask symetric, +1 im > sym(im) , -1 im < sym(im)
    PB701_PP *PP  ; // description for each column partition in connected masks
    int *firstPPByMask ; // index of first PP by column.
    PB701_TRANSF *PT ; // description of transformation to add a column to a partition
} PB701_DATA ;

void PB701_Init(PB701_DATA *data701, int sizeColumn) {
    data701->nbCol = 1 << sizeColumn ;
    int *nbBit ;
    nbBit = malloc(data701->nbCol*sizeof(nbBit[0]));
    // compute nb bits by mask
    for(int ic=0;ic<data701->nbCol;ic++) {
        int nb = 0 ;
        for(int j=0;j < sizeColumn ;j++) {
            if(ic & (1 << j)) nb++ ;
        }
        nbBit[ic] = nb ;
    }
    // compute symetry of mask
    data701->isSymetric = malloc(data701->nbCol*sizeof(data701->isSymetric[0]));
    for(int ic=0;ic<data701->nbCol;ic++) {
        int iSym =0 ;
        for(int ib=0;ib<sizeColumn;ib++) {
            iSym |= ((ic >> ib) & 1) << (sizeColumn -ib -1) ;
        }
        if(iSym == ic) data701->isSymetric[ic] = 0 ;
        else if(ic > iSym) data701->isSymetric[ic] = 1;
        else data701->isSymetric[ic] = -1;
    }
    // compute intersection of 2 mask ,extended to connected parts for each mask.
    PB701_INTER *interSec ; // intersection of 2 masks extended to connexion parts
    interSec=malloc(data701->nbCol*data701->nbCol*sizeof(interSec[0]));
    for(int ic=0;ic<data701->nbCol;ic++) {
        for(int jc=0;jc<data701->nbCol;jc++) {
            int km = ic & jc ;
            if(km==0) {
                interSec[ic*data701->nbCol+jc].newContact = 0 ;
                interSec[ic*data701->nbCol+jc].surfIncrement = 0 ;
                continue ;
            }
            // extension to connected parts
            for(int m=0;m<sizeColumn;m++) {
                if(m>0  && ( (1<<(m-1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
                if(m<sizeColumn-1  && ( (1<<(m+1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
            }
            for(int m=sizeColumn-1;m>=0;m--) {
                if(m>0  && ( (1<<(m-1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
                if(m<sizeColumn-1  && ( (1<<(m+1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
            }
            interSec[ic*data701->nbCol+jc].surfIncrement = nbBit[km] ;
            interSec[ic*data701->nbCol+jc].newContact = (uint16_t) km ;
        }
    }
    // compute connected partition for each column
    PB701_CONN *PC = malloc(data701->nbCol*sizeof(PC[0])) ;
    int nbPPexp = 0 ;
    for(int ic=0;ic<data701->nbCol;ic++) {
        uint16_t Set[PB701_MAXCONNECT]  ;
        memset(Set,0,sizeof(Set));
        int nbConnect = 0 ;
        int n = ic & 1 ;
        for(int j=1;j<=sizeColumn;j++) {
            if(ic & (1 << j)) {
                n++ ;
            } else {
                if(n) {
                    Set[nbConnect] = (uint16_t)((1 << j) - (1<<(j-n))) ;
                    nbConnect++ ;
                    n=0 ;
                }
            }
        }
        assert(nbConnect <= sizeColumn);
        int ppByConnected[PB701_MAXCONNECT+1] = { 1,1,2,5,14,42} ;
        nbPPexp += ppByConnected[nbConnect] ;
        PC[ic].globalContact = ic ;
        PC[ic].nbSet = nbConnect ;
        memcpy(PC[ic].Set,Set,sizeof(Set)) ;
    }
    data701->PP = malloc(nbPPexp*sizeof(data701->PP[0]));
    data701->firstPPByMask = malloc((data701->nbCol+1)*sizeof(data701->firstPPByMask[0])) ;
    // compute different partitions of column in connected parts
    // dont add impossible configurations as aaabbbaabbb (entrelacement not planar)
    data701->nbPP = 0 ;
    data701->firstPPByMask[0] = data701->nbPP ;
    for(int ic=0;ic<data701->nbCol;ic++) {
        PB701_PP *curPP = data701->PP + data701->nbPP ;
        uint16_t Set[PB701_MAXCONNECT]  ;
        memcpy(Set,PC[ic].Set,sizeof(Set)) ;
        int nbConnect = PC[ic].nbSet ;
#define INIT_PART(n)   curPP = data701->PP + data701->nbPP ; curPP->num = data701->nbPP++ ;  curPP->globalContact = ic ; memset(curPP->Set,0,sizeof(curPP->Set)) ; curPP->nbSet = (n)
        
#define DST(i)  curPP->Set[i]
#define SRC(i)  Set[i]
        
#define PART_A(i0)         INIT_PART(1) ; DST(0)=SRC(i0)
        
#define PART_AB(i0,i1)     INIT_PART(1) ; DST(0)=SRC(i0)|SRC(i1)
#define PART_A_B(i0,i1)    INIT_PART(2) ; DST(0)=SRC(i0) ; DST(1)=SRC(i1)
        
#define PART_ABC(i0,i1,i2)     INIT_PART(1) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2)
#define PART_AB_C(i0,i1,i2)    INIT_PART(2) ; DST(0)=SRC(i0)|SRC(i1) ; DST(1)=SRC(i2)
#define PART_A_B_C(i0,i1,i2)   INIT_PART(3) ; DST(0)=SRC(i0) ; DST(1)=SRC(i1) ; DST(2)=SRC(i2)
        
#define PART_ABCD(i0,i1,i2,i3)     INIT_PART(1) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3)
#define PART_AB_CD(i0,i1,i2,i3)    INIT_PART(2) ; DST(0)=SRC(i0)|SRC(i1) ; DST(1)=SRC(i2)|SRC(i3)
#define PART_ABC_D(i0,i1,i2,i3)    INIT_PART(2) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2) ; DST(1)=SRC(i3)
#define PART_AB_C_D(i0,i1,i2,i3)   INIT_PART(3) ; DST(0)=SRC(i0)|SRC(i1) ; DST(1)=SRC(i2) ; DST(2)=SRC(i3)
#define PART_A_B_C_D(i0,i1,i2,i3)  INIT_PART(4) ; DST(0)=SRC(i0) ; DST(1)=SRC(i1) ; DST(2)=SRC(i2) ; DST(3)=SRC(i3)
        
#define PART_ABCDE(i0,i1,i2,i3,i4)     INIT_PART(1) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3)|SRC(i4)
#define PART_ABCD_E(i0,i1,i2,i3,i4)    INIT_PART(2) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3) ; DST(1) = SRC(i4)
#define PART_ABC_DE(i0,i1,i2,i3,i4)    INIT_PART(2) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2) ; DST(1)=SRC(i3)|SRC(i4)
#define PART_ABC_D_E(i0,i1,i2,i3,i4)   INIT_PART(3) ; DST(0)=SRC(i0)|SRC(i1)|SRC(i2) ; DST(1)=SRC(i3) ; DST(2)=SRC(i4)
#define PART_AB_CD_E(i0,i1,i2,i3,i4)   INIT_PART(3) ; DST(0)=SRC(i0)|SRC(i1) ; DST(1)=SRC(i2)|SRC(i3) ; DST(2)=SRC(i4)
#define PART_AB_C_D_E(i0,i1,i2,i3,i4)  INIT_PART(4) ; DST(0)=SRC(i0)|SRC(i1) ; DST(1)=SRC(i2) ; DST(2)=SRC(i3) ; DST(3)=SRC(i4)
#define PART_A_B_C_D_E(i0,i1,i2,i3,i4) INIT_PART(5) ; DST(0)=SRC(i0) ; DST(1)=SRC(i1) ; DST(2)=SRC(i2) ; DST(3)=SRC(i3) ; DST(4)=SRC(i4)
        switch(nbConnect) {
            case 0 :
                INIT_PART(0) ;
                break ;
            case 1 :
                PART_A(0) ;
                break ;
            case 2:
                PART_AB(0,1) ;PART_A_B(0,1);
                break ;
            case 3:
                PART_ABC(0,1,2);
                PART_AB_C(0,1,2) ; PART_AB_C(0,2,1) ; PART_AB_C(1,2,0) ;
                PART_A_B_C(0,1,2);
                break ;
            case 4:
                PART_ABCD(0,1,2,3);
                PART_ABC_D(0,1,2,3)  ; PART_ABC_D(0,1,3,2) ;PART_ABC_D(0,2,3,1) ;PART_ABC_D(1,2,3,0) ;
                PART_AB_CD(0,1,2,3)  ; PART_AB_CD(0,3,1,2) ;
                PART_AB_C_D(0,1,2,3) ; PART_AB_C_D(0,2,1,3);PART_AB_C_D(0,3,1,2);PART_AB_C_D(1,2,0,3);PART_AB_C_D(1,3,0,2);PART_AB_C_D(2,3,0,1) ;
                PART_A_B_C_D(0,1,2,3);
                break ;
            case 5:
                PART_ABCDE(0,1,2,3,4);
                PART_ABCD_E(0,1,2,3,4)  ;PART_ABCD_E(0,1,2,4,3)  ;PART_ABCD_E(0,1,3,4,2)  ;PART_ABCD_E(0,2,3,4,1)  ;PART_ABCD_E(1,2,3,4,0) ;
                PART_ABC_DE(0,1,2,3,4)  ;PART_ABC_DE(0,1,4,2,3)  ;PART_ABC_DE(0,3,4,1,2)  ;PART_ABC_DE(2,3,4,0,1)  ;PART_ABC_DE(1,2,3,0,4) ;
                PART_ABC_D_E(0,1,2,3,4) ;PART_ABC_D_E(0,1,3,2,4) ;PART_ABC_D_E(0,2,3,1,4) ;PART_ABC_D_E(1,2,3,0,4) ;PART_ABC_D_E(0,1,4,2,3);
                PART_ABC_D_E(0,2,4,1,3) ;PART_ABC_D_E(1,2,4,0,3) ;PART_ABC_D_E(0,3,4,1,2) ;PART_ABC_D_E(1,3,4,0,2) ;PART_ABC_D_E(2,3,4,0,1);
                PART_AB_CD_E(0,1,2,3,4) ;PART_AB_CD_E(0,3,1,2,4) ;PART_AB_CD_E(0,1,2,4,3) ;PART_AB_CD_E(0,4,1,2,3) ;PART_AB_CD_E(0,1,3,4,2);
                PART_AB_CD_E(0,4,1,3,2) ;PART_AB_CD_E(0,2,3,4,1) ;PART_AB_CD_E(0,4,2,3,1) ;PART_AB_CD_E(1,2,3,4,0) ;PART_AB_CD_E(1,4,2,3,0);
                PART_AB_C_D_E(0,1,2,3,4);PART_AB_C_D_E(0,2,1,3,4);PART_AB_C_D_E(0,3,1,2,4);PART_AB_C_D_E(0,4,1,2,3);PART_AB_C_D_E(1,2,0,3,4);
                PART_AB_C_D_E(1,3,0,2,4);PART_AB_C_D_E(1,4,0,2,3);PART_AB_C_D_E(2,3,0,1,4);PART_AB_C_D_E(2,4,0,1,3);PART_AB_C_D_E(3,4,0,1,2) ;
                PART_A_B_C_D_E(0,1,2,3,4);
                break ;
        }
        // sort each partition set for future comparison between partitions
        data701->firstPPByMask[ic+1] = data701->nbPP ;
        for(int ip = data701->firstPPByMask[ic];ip<data701->firstPPByMask[ic+1];ip++ ) {
            int nbs = data701->PP[ip].nbSet ;
            for(int is = 1;is<nbs;is++) {
                for(int js=is;js<nbs;js++) {
                    if(data701->PP[ip].Set[is-1] > data701->PP[ip].Set[js]) {
                        uint16_t tmp = data701->PP[ip].Set[is-1] ;
                        data701->PP[ip].Set[is-1] = data701->PP[ip].Set[js] ;
                        data701->PP[ip].Set[js] = tmp ;
                    }
                }
            }
        }
        data701->firstPPByMask[ic+1] = data701->nbPP ;
    }
    // compute for each partition every transformation to add a new column
    data701->PT = calloc(data701->nbPP*data701->nbCol,sizeof(data701->PT[0]));
    for(int ip=0;ip<data701->nbPP;ip++) { // loop on partition
        for(int jc=0;jc<data701->nbCol;jc++) { // lopp on supplementary partition
            int nbSet = 0 ;
            uint16_t Set[PB701_MAXCONNECT]  ;
            memset(Set,0,sizeof(Set));
            int globalContact = 0 ;
            for(int is=0;is<data701->PP[ip].nbSet;is++) {
                // for each connected part search the new set of connected parts in new column
                // remove all intersecting precedent sets and add the new set.
                // respect ascending order in sets for future comparison
                uint16_t newContact = interSec[data701->PP[ip].Set[is]*data701->nbCol+jc].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) { // remove intersecting sets
                        if((newContact & Set[jj])) {
                            newContact = interSec[(newContact | Set[jj]) * data701->nbCol+jc].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                Set[l-1] = Set[l] ;
                            }
                            nbSet--; jj-- ;
                            Set[nbSet] = 0 ;
                        }
                    }
                    int l ; // add new set in good place
                    for(l=nbSet-1;l>=0 && Set[l] > newContact  ;l--) {
                        Set[l+1] = Set[l] ;
                    }
                    Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }
            if(globalContact != jc) { // add missing connected parts of new column
                globalContact ^= jc ;
                for(int k=0;k<PC[jc].nbSet;k++) {
                    if((PC[jc].Set[k] & globalContact) == PC[jc].Set[k] ){
                        int l ;
                        for(l=nbSet-1;l>=0 && Set[l] > PC[jc].Set[k]   ;l--) {
                            Set[l+1] = Set[l] ;
                        }
                        Set[l+1] = PC[jc].Set[k] ;
                        nbSet++ ;
                    }
                }
            }
            int in ; // find the resulting PP (if missing BIB BUG !!)
            for(in=data701->firstPPByMask[jc];in<data701->firstPPByMask[jc+1];in++) {
                if(memcmp(Set,data701->PP[in].Set,sizeof(Set))==0) break ;
            }
            assert(in != data701->firstPPByMask[jc+1]) ;
            data701->PT[ip*data701->nbCol+jc].numPP = in ;
            // compute the transformation
            // delta[] is the size of set (number ob black case to add)
            for(int id=0;id < data701->PP[in].nbSet ; id++) data701->PT[ip*data701->nbCol+jc].delta[id] =  nbBit[data701->PP[in].Set[id]] ;
            for(int is=0;is<data701->PP[ip].nbSet;is++) {
                if(data701->PP[ip].Set[is] & jc) { // numCv[] is the correspondance of index from old PP to new PP.
                    for(int id=0;id < data701->PP[in].nbSet ; id++){
                        if( data701->PP[ip].Set[is] & data701->PP[in].Set[id] ) {
                            data701->PT[ip*data701->nbCol+jc].numCv[is] = id+1 ;
                            break ;
                        }
                        assert(id < data701->PP[in].nbSet) ;
                    }
                }
            }
        }
    }
    free(PC) ;
    free(nbBit) ;
    free(interSec);
}
int PB701(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    PB701_DATA data701 ;
    int sizeColumn = PB701_L ;
    int maxNbColumn =  PB701_H ;
    PB701_Init(&data701,sizeColumn) ; // precompute all transformation to add a column
    int nbCol = data701.nbCol ; // recover
    PB701_PP * PP = data701.PP ;
    PB701_TRANSF * PT = data701.PT;
    
    PB701_STATE *State= calloc(300000000,sizeof(State[0])) ;
    int nbState = 0;
    NsumInd * PPindState = calloc(data701.nbPP,sizeof(PPindState[0])) ;
    int firstNbState , lastNbState ;
    firstNbState = 0 ;
    State[0].prob = 1 ; State[0].numPP = 0 ;
    nbState = lastNbState = 1;
    double Sres ;
    for(int nbColumn=1;nbColumn< maxNbColumn ;nbColumn++) { // loops on column
        int N=sizeColumn*nbColumn ;
        Nsum * NS= NsumAlloc(N,PB701_MAXCONNECT+1) ;
        Prob701 S=0 , ST = 0 ;
        int maxHash = 0 ;
        int maxState = 0 ;
        for(int jc=0;jc<nbCol;jc++) { // loop on differents column the current added column
            int isDouble = 0 ;
            if(nbColumn==maxNbColumn-1) { // if last one, check symetry
                if(data701.isSymetric[jc]< 0) continue ;
                isDouble = data701.isSymetric[jc] ? 1 : 0 ;
            }
            // compute size of exact hash for current added column
            int SizeHash = 0 ;
            int firstStateForColumn = nbState ;
            for(int ip=data701.firstPPByMask[jc];ip<data701.firstPPByMask[jc+1];ip++) {
                PPindState[ip] = SizeHash ;
                SizeHash += NsumGetSize(NS,PP[ip].nbSet+1)  ;
            }
            // alloc hash states merge
            NsumInd * hashState = calloc(SizeHash,sizeof(hashState[0])) ;
            if(SizeHash > maxHash) maxHash = SizeHash ;
            for(int istate=firstNbState;istate < lastNbState;istate++) { // loop on state from preceedind column
                PB701_STATE *antState = State+istate ;
                int ip = antState->numPP ;
                PB701_TRANSF *trf = PT + ip*nbCol + jc ;
                int id = trf->numPP ;
                PB701_STATE *newState = State+nbState ;
                NsumVal * surf = newState->Surf;
                newState->numPP = id ;
                surf[0] = antState->Surf[0] ;
                for(int is=0;is<PP[id].nbSet;is++) {
                    surf[is+1] = trf->delta[is] ;
                }
                for(int is=0;is<PP[ip].nbSet;is++) {
                    if(trf->numCv[is]) {
                        surf[trf->numCv[is]] += antState->Surf[is+1];
                    } else {
                        if(antState->Surf[is+1] >  surf[0] ) surf[0] = antState->Surf[is+1] ;
                    }
                }
                NsumInd numHash = PPindState[id] + NsumGetIndex(NS, PP[id].nbSet+1, surf) ;
                int newIstate = hashState[numHash] ;
                if(newIstate == 0) { // not used state
                    hashState[numHash]= nbState++ ;
                    newState->prob = isDouble ? 2*antState->prob  : antState->prob ;
                } else { // used state , cumulate only prob
                    State[newIstate].prob += isDouble ? 2*antState->prob  : antState->prob  ;
                }
            }
            free(hashState);
            for(int istate = firstStateForColumn ; istate<nbState; istate++) {
                PB701_STATE * newState = State + istate  ;
                int maxS = newState->Surf[0];
                for(int is=0;is<PP[newState->numPP].nbSet;is++) {
                    if(newState->Surf[is+1] > maxS) maxS = newState->Surf[is+1]  ;
                }
                S += maxS * newState->prob ;
            }
            if(nbColumn==maxNbColumn-1) { // last column, special case not, necessary to save new states)
                for(int istate=firstStateForColumn;istate < nbState;istate++) {
                    PB701_STATE *antState = State + istate  ;
                    int ip = antState->numPP ;
                    for(int jc=0;jc<nbCol;jc++) { // [0,1[
                        PB701_TRANSF *trf = PT + ip*nbCol + jc ;
                        int id = trf->numPP ;
                        NsumVal surf[PB701_MAXCONNECT] ;
                        int maxS = antState->Surf[0] ;
                        for(int is=0;is<PP[id].nbSet;is++) {
                            surf[is] = trf->delta[is] ;
                        }
                        for(int is=0;is<PP[ip].nbSet;is++) {
                            if(trf->numCv[is]) {
                                surf[trf->numCv[is]-1] += antState->Surf[is+1];
                            } else {
                                if(antState->Surf[is+1] >  maxS ) maxS = antState->Surf[is+1] ;
                            }
                        }
                        for(int is=0;is<PP[id].nbSet;is++) {
                            if(surf[is] > maxS) maxS = surf[is]  ;
                        }
                        ST += maxS * antState->prob;
                    }
                }
                if(nbState > maxState) maxState = nbState ;
                nbState = firstStateForColumn ; // remins, the used states for last pass are not necessary
            }
        }
        
        int exp2 = (nbColumn)*sizeColumn ;
        Sres = (double)S / (double) (1LL << (exp2/2)) / (double)(1LL << ((exp2+1)/2)) ;
#if defined PB701_BIG
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s \t%.3fs %dx%d States [%d -> %d[ HashSize=%d S=%.8f S=%llu%llu\n",
                                   pbR->ident,(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC ,sizeColumn,nbColumn,
                                   firstNbState,lastNbState,maxHash
                                   ,Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL));
#else
        if(pbR->isVerbose) fprintf(stdout,"\tPB%s \t%.3fs %dx%d States [%d -> %d[ HashSize=%d S=%.8f S=%llu\n",
                                   pbR->ident,(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC ,sizeColumn,nbColumn,
                                   firstNbState,lastNbState,maxHash
                                   ,Sres,S);
#endif
        firstNbState = lastNbState ;
        lastNbState = nbState ;
        NS=NsumFree(NS);
        if(nbColumn == maxNbColumn-1) {
            exp2 = (nbColumn+1)*sizeColumn ;
            Sres = (double)ST / (double) (1LL << (exp2/2 )) / (double)(1LL << ((exp2+1)/2)) ;
#if defined PB701_BIG
            if(pbR->isVerbose) fprintf(stdout,"\tPB%s %.3fs LxH=%dx%d MaxStates=%d S=%.8f S=%llu%llu\n",pbR->ident
                                      ,(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC
                   ,sizeColumn,maxNbColumn,maxState
                   ,Sres,(int64_t)(ST / 1000000000000000LL), (int64_t)(ST % 1000000000000000LL));
#else
            if(pbR->isVerbose) fprintf(stdout,"\tPB%s %.3fs LxH=%dx%d MaxStates=%d S=%.8f S=%%llu\n",pbR->ident
                                       ,(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC
                                       ,sizeColumn,maxNbColumn,maxState
                                       ,Sres,ST);
#endif
        }
    }
    free(PPindState);
    free(State);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.8f",Sres);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB702_MAX   5
#define PB702_NBT   5

typedef struct PB702_POINT {
    uint32_t ix ;
    uint32_t iy ;
} PB702_POINT ;

static int CmpPoint(const void * el1, const void *el2) {
    const PB702_POINT * p1 =(PB702_POINT *)el1 ;
    const PB702_POINT * p2 =(PB702_POINT *)el2 ;
    int diff = p1->ix - p2->ix ;
    if(diff) return diff ;
    return p1->iy - p2->iy ;
}

int PB702b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nbPoints = 1 ;
    int nl= PB702_NBT ;
    for(int i=0;i<=PB702_MAX;i++) nbPoints *=6 ;
    nbPoints = (nbPoints-1)/5 ;
    PB702_POINT * points = malloc(nbPoints*sizeof(points[0])) ;
    points[0] = (PB702_POINT) {0,0} ;
    int nbPointsAnt = 0 ;
    int nbPointCur = 1;
    nbPoints = 1;
    int ixMax = 1 << (PB702_MAX+1) ;
    int ixMax2 = ixMax >> 1 ;
    int iyMax = 1 << PB702_MAX ;
    int8_t *pointsDesc = calloc((ixMax+1)*(iyMax+1),sizeof(pointsDesc[0])) ;
    int16_t *pointsTriangle = calloc((ixMax+1)*(iyMax+1),sizeof(pointsTriangle[0])) ;
    
    for(int iy=iyMax;iy>=0;iy--) {
        for(int ix=0;ix<=ixMax;ix++) {
            int i1,i2,i3,im = 1 ;
            i1 = (nl*iy)/ iyMax;
            if(nl*iy == i1*iyMax) im = 0 ;
            i2 = ((2*iyMax+iy - ix)*nl)/ixMax ;
            if((2*iyMax+iy - ix)*nl == i2 * ixMax) im = 0 ;
            i3 = nl*(ix+iy)/ixMax;
            if(nl*(ix+iy) == i3 * ixMax) im =0 ;
            if(im) im = (i1 ^ i2 ^ i3 ^ nl) & 1 ;
            if(im) {
                int il = (iy * nl)/iyMax ;
                im = 1 + il * nl ;
                int ir = (il & 1) ? ((ix*nl+ixMax2)/ ixMax) : ((ix  * nl) / ixMax) ;
                im += ir ;
                if((il & 1) && ir == 0) im = -im ;
            }
            pointsTriangle[iy*(ixMax+1)+ix] = im ;
        }
    }
    int *isOktriangle = calloc(nl*nl+1,sizeof(isOktriangle[0]));
    int64_t sumB = 0 ;
    pointsDesc[0] = 1 ;
    for(int ipass=0;ipass<PB702_MAX;ipass++) {
        for(int ip=nbPointsAnt;ip<nbPointCur;ip++) {
            PB702_POINT oldPoint = points[ip] ;
            PB702_POINT newPoint ;
            newPoint.ix = (ixMax - oldPoint.ix) >> 1 ;
            newPoint.iy = (oldPoint.iy) >> 1 ;
            points[nbPoints++] = newPoint ;
            newPoint.ix = (ixMax2 > oldPoint.ix ) ?  ((ixMax2 - oldPoint.ix) >> 1)  : ((oldPoint.ix - ixMax2) >> 1);
            newPoint.iy = (oldPoint.iy+iyMax) >> 1 ;
            points[nbPoints++] = newPoint ;
            newPoint.ix = (ixMax2 + oldPoint.ix) >> 1 ;
            newPoint.iy = (oldPoint.iy+iyMax) >> 1 ;
            points[nbPoints++] = newPoint ;
            newPoint.ix = (ixMax + oldPoint.ix) >> 1 ;
            newPoint.iy = (oldPoint.iy) >> 1 ;
            points[nbPoints++] = newPoint ;
            newPoint.ix = (ixMax2 + oldPoint.ix) >> 1 ;
            newPoint.iy = (-oldPoint.iy+iyMax) >> 1 ;
            points[nbPoints++] = newPoint ;
            newPoint.ix = (ixMax2 > oldPoint.ix ) ?  ((ixMax2 - oldPoint.ix) >> 1)  : ((oldPoint.ix - ixMax2) >> 1);
            newPoint.iy = (-oldPoint.iy+iyMax) >> 1 ;
            points[nbPoints++] = newPoint ;
        }
        nbPointsAnt = nbPointCur ;
        nbPointCur = nbPoints ;
        qsort(points+nbPointsAnt,nbPointCur-nbPointsAnt,sizeof(points[0]),CmpPoint) ;
        nbPoints = nbPointsAnt ;
        for(int i=nbPointsAnt+1;i<nbPointCur;i++) {
            if(CmpPoint(points+nbPoints,points+i)==0 ) continue ;
            nbPoints++ ;
            if(nbPoints == i) continue ;
            points[nbPoints] = points[i] ;
        }
        nbPointCur = ++nbPoints ;
        printf("%d =>\n",nbPointCur-nbPointsAnt);
        for(int i=nbPointsAnt;i<nbPointCur;i++) {
            if(pointsDesc[points[i].iy*(ixMax+1)+points[i].ix]==0) pointsDesc[points[i].iy*(ixMax+1)+points[i].ix]= ipass+1 ;
        }
        for(int iy=iyMax;iy>=0;iy--) {
            for(int ix=0;ix<=ixMax;ix++) printf("%c",(pointsDesc[iy*(ixMax+1)+ix] && pointsTriangle[iy*(ixMax+1)+ix]) ?
                                                ((pointsDesc[iy*(ixMax+1)+ix]==ipass+1) ? '0'+ipass+1 : '+'  ): ((pointsDesc[iy*(ixMax+1)+ix]==ipass+1)? '*': '.'));
            printf("\n") ;
        }
 
        for(int iy=iyMax;iy>=0;iy--) {
            for(int ix=0;ix<=ixMax;ix++) {
                if(pointsTriangle[iy*(ixMax+1)+ix] && (pointsDesc[iy*(ixMax+1)+ix]==ipass+1)) {
                    if(pointsTriangle[iy*(ixMax+1)+ix] > 0) {
                        if(isOktriangle[pointsTriangle[iy*(ixMax+1)+ix]] == 0) {
                            isOktriangle[pointsTriangle[iy*(ixMax+1)+ix]] = ipass + 1;
                            sumB += (ipass+1) * 2 ;
 //                           printf("T(%d,%d),n%d+(%d),sum=%lld ",ix,iy,pointsTriangle[iy*(ixMax+1)+ix],ipass+1,sumB) ;
                        } else {
 //                           printf("Ident(%d,%d)n%d ",ix,iy,pointsTriangle[iy*(ixMax+1)+ix] ) ;
                        }
                    } else {
                        if(isOktriangle[-pointsTriangle[iy*(ixMax+1)+ix]] == 0) {
                            isOktriangle[-pointsTriangle[iy*(ixMax+1)+ix]] = ipass + 1;
                            sumB += (ipass+1)  ;
  //                          printf("H(%d,%d),n%d+(%d),sum=%lld ",ix,iy,pointsTriangle[iy*(ixMax+1)+ix],ipass+1,sumB) ;
                        }
                        
                    }
                }
            }
        }
        
        
    }
    printf("\n");
    for(int iy=iyMax;iy>=0;iy--) {
        for(int ix=0;ix<=ixMax;ix++) printf("%c",( pointsDesc[iy*(ixMax+1)+ix] && pointsTriangle[iy*(ixMax+1)+ix])
                                            ? '0'+pointsDesc[iy*(ixMax+1)+ix] : (pointsTriangle[iy*(ixMax+1)+ix]? '+' :  '.'));
        printf("\n") ;
    }
 
    printf("SumB = %lld\n",sumB);
    
    uint8_t *isTok = malloc((nl+1)*sizeof(isTok[0]));
    int64_t sum = 0 ;
    for(int il=0;il<nl;il++) {
        int nbt = nl-(il+1)/2 ;
        int halfNbt = il & 1 ;
        int isHalfAlive = halfNbt ;
        int nbTalive = nbt  ;
//        printf("il=%d nbt=%d alive=%d ",il,nbt,nbTalive+isHalfAlive);
        memset(isTok,0,nbTalive *sizeof(isTok[0])) ;
        for(int ik = 2;nbTalive+isHalfAlive;ik++) {
            int kpow2 = 1 << ik ;
//            printf("P=%d ",kpow2);
            // facteur echell 2**in
            int ilklow = (il * kpow2)/ nl +1 ;
            //         if(ilklow * nl == il * kpow2) ilklow++ ;
            int ilkhigh = ((il+1)*kpow2) / nl +1 ;
//            printf("ilk=[%d,%d[ ",ilklow,ilkhigh);
            for(int ilk = ilklow;ilk<ilkhigh;ilk++) {
                int  idebT,deltaK, ifinT ; //  kpow2 * nl
                idebT = (nl * ilk  - il * kpow2) ;
                ifinT = (2 * kpow2 - ( nl * ilk -  il * kpow2)) ;
                if(ilk & 1) {
                    idebT += nl;
                    ifinT += nl ;
                    deltaK =  2* nl ;
                } else {
                    if((ilk & 2) == 0) {
                        idebT += 2* nl ;
                        ifinT += 2*nl ;
                    }
                    deltaK = 4 * nl ;
                }
                // ifinT - idebT = 2 * ( kpow2 - ( nl * ilk -  il * kpow2))
                if(halfNbt) {
                    // on traite le cas particulier du demi triangle
                    if(isHalfAlive) {
                        int ifinTH = ((ilk & 3) ==2) ? ifinT : ifinT - deltaK ;
                        if (ifinTH - kpow2 >=0) {
                            int ikfin = (ifinTH- kpow2)/ deltaK ;
                            if((ifinTH- kpow2)  == deltaK *ikfin) ikfin-- ;
                            if(ikfin >=0) {
                                if(ik==4) printf("H[%d]",il) ;
                                isHalfAlive = 0 ;
                                sum += ik ;
//                                printf("(ilk=%d,ifinT=%d,ikfin=%d)H[0],sum=%lld ",ilk,ifinTH-kpow2,ikfin,sum) ;
                                
                            }
                        }
                    }
                    
                    idebT += kpow2 ;
                    ifinT += kpow2 ;
                    if(ifinT-idebT > deltaK) {
                        sum += nbTalive * 2 * ik ;
 //                       printf("{%dx%d}=0,sum=%lld ",nbTalive,ik,sum) ;
                        nbTalive = 0 ;
                        break ;
                    } else {
                        for(int nT = 0; nT < nbt ; nT++) {
                            int ikdeb = idebT/deltaK;
                            if(idebT == ikdeb*deltaK) ikdeb++ ;
                            int ikfin = ifinT/ deltaK ;
                            if(ifinT == ikfin*deltaK) ikfin-- ;
 //                           printf("[%d->%d,%d->%d] ",idebT,ikdeb,ifinT,ikfin);
                            
                            if(isTok[nT] == 0 && ikdeb < ikfin) {
                                isTok[nT] = 1 ;
                                if(ik==4) printf("(%d,%d)",il,nT);
                                sum += 2*ik;
  //                              printf("{%d}=%d,sum=%lld ",nT,ik,sum) ;
                                nbTalive-- ;
                            }
                            idebT += 2 * kpow2 ;
                            ifinT += 2 * kpow2 ;
                        }
                    }
                } else {
                    if(ifinT-idebT > deltaK) {
 //                       printf("<%d->%d,%d->%d>",idebT,idebT/deltaK,ifinT,ifinT/deltaK);
                        sum += nbTalive * 2 * ik ;
  //                      printf("[%dx%d]=0,sum=%lld ",nbTalive,ik,sum) ;
                        nbTalive = 0 ;
                        break ;
                    } else
                    {
                        for(int nT = 0; nT < nbt ; nT++) {
                            int ikdeb = idebT/deltaK;
                            if(idebT == ikdeb*deltaK) ikdeb++ ;
                            int ikfin = ifinT/ deltaK ;
                            if(ifinT == ikfin*deltaK) ikfin-- ;
  //                          printf("<%d->%d,%d->%d> ",idebT,ikdeb,ifinT,ikfin);
                            if(isTok[nT] == 0 && ikdeb < ikfin) {
                                isTok[nT] = 1 ;
                                if(ik==4) printf("{%d,%d}",il,nT);
                               sum += 2*ik;
   //                             printf("[%d]=%d,sum=%lld ",nT,ik,sum) ;
                                nbTalive-- ;
                            }
                            idebT += 2 * kpow2 ;
                            ifinT += 2 * kpow2 ;
                        }
                    }
                }
 //               printf(" *** ");
            }
 //           printf("\n");
        }
    }
    printf("Sum=%lld \n",sum);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB705_MAX 100000000
#define PB705_MOD  1000000007LL


int PB705a(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int64_t SumD1=0, SumD2=0,SumD3=0,SumD4=0,SumD5=0,SumD6=0,SumD7=0,SumD8=0,SumD9=0 ;
    int64_t cost = 0 ;
    int64_t totalWeight = 1 ;
    int numDig = 0;
    pbR->nbClock = clock() ;
    for(int i=0;i<nbPrime;i++) {
        int32_t p = tbPrime[i] ;
        uint8_t dig[20] ;
        int nbDig = 0 ;
        while(p > 0) {
            int ip = p /10 ;
            int d = p - 10 * ip ;
            if(d)dig[nbDig++] = d ;
            p = ip ;
        }
        for(int nd=nbDig-1;nd>=0;nd--) {
            int d = dig[nd] ;
            //                        printf("%d :",d);
            switch(d) {
                case 1 :
                    cost = cost + SumD1    ;
                    break ;
                case 2 :
                    cost = 2*cost+SumD1 +SumD2   ;
                    SumD1 = 2*SumD1 + totalWeight ;
                    SumD2 *= 2 ;  SumD3 *= 2;  SumD4 *= 2 ;  SumD5 *= 2 ;
                    SumD6 *= 2 ;  SumD7 *= 2;  SumD8 *= 2 ;  SumD9 *= 2 ;
                    totalWeight = ( totalWeight * 2 );
                    break ;
                case 3 :
                    cost = 2*cost +SumD1 +SumD3    ;
                    SumD1 = 2*SumD1 + totalWeight ; SumD2 = 2*SumD2 + totalWeight ;
                    SumD3 *= 2 ;  SumD4 *= 2 ;  SumD5 *= 2 ; SumD6 *= 2;
                    SumD7 *= 2 ;  SumD8 *= 2 ;    SumD9 *= 2 ;
                    totalWeight = ( totalWeight * 2 ) ;
                    break ;
                case 4 :
                    cost = 3*cost+ +SumD1 +SumD2 +SumD4  ;
                    SumD1 = 3*SumD1 + 2*totalWeight ;
                    SumD2 = 3*SumD2 + totalWeight ; SumD3 = 3*SumD3 + totalWeight ;
                    SumD4 *= 3 ; SumD5 *= 3 ;   SumD6 *= 3 ;
                    SumD7 *= 3 ; SumD8 *= 3 ;   SumD9 *= 3 ;
                    totalWeight *= 3 ;
                    break ;
                case 5 :
                    cost = 2*cost +SumD1 +SumD5  ;
                    SumD1 = 2*SumD1 + totalWeight ;  SumD2 = 2*SumD2 + totalWeight ;
                    SumD3 = 2*SumD3 + totalWeight ;  SumD4 = 2*SumD4 + totalWeight ;
                    SumD5 *= 2 ; SumD6 *= 2 ;  SumD7 *= 2 ; SumD8 *= 2 ;SumD9 *= 2 ;
                    totalWeight *=  2  ;
                    break ;
                case 6 :
                    cost = 4*cost +SumD1 +SumD2 +SumD3 +SumD6 ;
                    SumD1 = 4*SumD1 + 3*totalWeight ;
                    SumD2 = 4*SumD2 + 2*totalWeight ;
                    SumD3 = 4*SumD3 + totalWeight ; SumD4 = 4*SumD4 + totalWeight ;
                    SumD5 = 4*SumD5 + totalWeight ;
                    SumD6 *= 4 ; SumD7 *= 4 ; SumD8 *= 4 ; SumD9 *= 4 ;
                    totalWeight *= 4 ;
                    break ;
                case 7 :
                    cost = 2*cost +SumD1 +SumD7   ;
                    SumD1 = 2*SumD1 + totalWeight ;  SumD2 = 2*SumD2 + totalWeight ;
                    SumD3 = 2*SumD3 + totalWeight ;  SumD4 = 2*SumD4 + totalWeight ;
                    SumD5 = 2*SumD5 + totalWeight ;  SumD6 = 2*SumD6 + totalWeight ;
                    SumD7 *= 2  ;  SumD8 *= 2  ; SumD9 *= 2  ;
                    totalWeight *=  2   ;
                    break ;
                case 8 :
                    cost = 4*cost +SumD1 +SumD2 +SumD4 +SumD8   ;
                    SumD1 = 4*SumD1 + 3*totalWeight ;
                    SumD2 = 4*SumD2 + 2*totalWeight ; SumD3 = 4*SumD3 + 2*totalWeight ;
                    SumD4 = 4*SumD4 + totalWeight ;  SumD5 = 4*SumD5 + totalWeight ;
                    SumD6 = 4*SumD6 + totalWeight ;  SumD7 = 4*SumD7 + totalWeight ;
                    SumD8 *= 4 ; SumD9 *= 4 ;
                    totalWeight = ( totalWeight * 4 )  ;
                    break ;
                case 9 :
                    cost = 3*cost +SumD1 +SumD3 +SumD9    ;
                    SumD1 = 3*SumD1 + 2*totalWeight ;
                    SumD2 = 3*SumD2 + 2*totalWeight ;
                    SumD3 = 3*SumD3 + totalWeight ; SumD4 = 3*SumD4 + totalWeight;
                    SumD5 = 3*SumD5 + totalWeight ; SumD6 = 3*SumD6 + totalWeight;
                    SumD7 = 3*SumD7 + totalWeight ; SumD8 = 3*SumD8 + totalWeight;
                    SumD9 *= 3 ;
                    totalWeight *=  3   ;
                    break ;
            }
            numDig++ ;
            if((numDig & 0x7) ==0  )
            {
                cost = cost %  PB705_MOD ;
                totalWeight = totalWeight % PB705_MOD ;
                SumD1 = SumD1  % PB705_MOD ;  SumD2 = SumD2  % PB705_MOD ;  SumD3 = SumD3  % PB705_MOD ;
                SumD4 = SumD4  % PB705_MOD ;  SumD5 = SumD5  % PB705_MOD ;  SumD6 = SumD6  % PB705_MOD ;
                SumD7 = SumD7  % PB705_MOD ;  SumD8 = SumD8  % PB705_MOD ;  SumD9 = SumD9  % PB705_MOD ;
            }
            /*
             printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld TW=%lld Cost=%lld\n"
             ,SumD1,SumD2,SumD3,SumD4,SumD5
             ,SumD6,SumD7,SumD8,SumD9,totalWeight,cost % PB705_MOD);
             */
        }
    }
    cost = cost %  PB705_MOD ;
    sprintf(pbR->strRes,"%lld",cost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB705(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    pbR->nbClock = clock() ;
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
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\tPB%s for Max=%d Nbinv=%lld (%%%lld) = %lld x 2**%lld x 3**%lld\n",
                              pbR->ident,PB705_MAX,costMod,PB705_MOD,cost,nbPow2-3,nbPow3-2);
    sprintf(pbR->strRes,"%lld",costMod) ;
    return 1 ;
}

#define PB708_MAX   100000000000000
// version brute n'aboutit pas pour la value demandee
int PB708(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t N2 = Sqrt64(N)+1 ;
    uint8_t * Nb2 = calloc(N2+1,sizeof(Nb2[0]));
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    int nbPrim = 0 ;
    int32_t *tbPrim = malloc(N2*sizeof(tbPrim[0])/10) ;
    for(int p=2;p<=N2;p++) {
        if(Nb2[p]) {
            histPow2[Nb2[p]]++ ;
            continue ;
        }
        histPow2[1]++ ; tbPrim[nbPrim++] = p ;
        for(int64_t pp=p;pp<=N2;pp *=p) {
            for(int64_t npp=pp;npp<=N2;npp+=pp) {
                Nb2[npp]++ ;
            }
        }
    }
    for(int64_t n=N2+1;n<=N;n++) {
        int nbP = 0 ;
        int64_t q = n ;
        for(int np=0;np<nbPrim && q > 1 ;np++) {
            int p = tbPrim[np] ;
            if(p*(int64_t)p > n) break ;
            while ( (q % p) == 0) {
                q /= p ;
                nbP++ ;
            }
        }
        if(q > 1) nbP++ ;
        histPow2[nbP]++ ;
        
    }
    int64_t S = 1 ;
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB708_MAXLEVEL  10
typedef struct LEVEL708 {
    int32_t ip ;        // prime index for the level
    int32_t ipow ;      // exponant for prime p minus 2. (as p*p is minima)
    int64_t p ;         // prime for the level
    int64_t inv_xp ;    // N / Product(pi ** expi)
} LEVEL708 ;

int64_t Alpha1(int64_t n) {
    int64_t alpha = 0 ;
    if(n < 0x7fffffff) {
        int32_t n32 = (int32_t) n ;
        for(int32_t d=1;d*(int64_t)d <= n32;d++ ) {
            alpha += n32/d ;
        }
        int32_t sqn = Sqrt32(n32) ;
        alpha = 2*alpha - sqn*sqn ;
    } else {
        for(int32_t d=1;d*(int64_t)d <= n;d++ ) {
            alpha += n/d ;
        }
        int64_t sqn = Sqrt64(n) ;
        alpha = 2*alpha - sqn*sqn ;
    }
    return alpha;
}

// version non recursive
// on developpe tous les produits avec au moins un carre pour chaque nombre premier
int PB708a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int32_t N2 = (int32_t) Sqrt64(N) ;
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int32_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    int32_t NL = (int32_t)pow(N/4,1/3.0) ;
    int64_t *alphaNL = malloc(NL*sizeof(alphaNL[0]));
    for(int in=1;in<NL;in++) { alphaNL[in] = Alpha1(in) ; }
    int64_t S = 0 ;
    LEVEL708 Lv[PB708_MAXLEVEL] ;
    Lv[0].inv_xp = N ;
    Lv[0].ipow = 0 ;
    Lv[0].ip = 0 ;
    Lv[0].p = 1 ;
    int is=0 ;
    while(1) {
        int64_t inv_xp = Lv[is].inv_xp ;
        int isUp = 0 ;
        if(inv_xp>=1) { // prime product don't exeed N
            if(inv_xp >= NL) {
                S += Alpha1(inv_xp) << Lv[is].ipow ;
            } else {
                S += alphaNL[inv_xp] <<  Lv[is].ipow ;
            }
            if(Lv[is].ip < nbPrime) { // not last prime < Sqrt(N)
                isUp = 1 ; // Up level
                Lv[is+1].ip = Lv[is].ip ; // next prime for up level
                is++ ;
            }
        }
        if (isUp || (Lv[is].ip < nbPrime && Lv[is].ipow > Lv[is-1].ipow) ) {
            // Up level or ( prime product exeed N and  ( Not last prime and p*p OK ))
            Lv[is].p = tbPrime[Lv[is].ip++] ; // next prime
            int64_t p2 = (Lv[is].p*Lv[is].p) ;
            if(p2 <= Lv[is-1].inv_xp) { // prime product don't exeed N
                // use next prime p => p*p
                Lv[is].inv_xp = Lv[is-1].inv_xp / p2 ;
                Lv[is].ipow = Lv[is-1].ipow ;
                continue ;
            }
        }
        // level is finished : last valid prime, and max exponant
        if(--is==0) break  ; // down level
        Lv[is].ipow++ ; // next exponant
        if(Lv[is].inv_xp >= Lv[is].p) Lv[is].inv_xp /= Lv[is].p ; // x p
        else Lv[is].inv_xp = 0 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int64_t Alpha_b(int64_t n) {
    static int32_t NL = 0 ;
    static int64_t *alphaL = NULL ;
    if(n < 0) {
        NL = (int32_t) (-n) ;
        alphaL = malloc(NL*sizeof(alphaL[0]));
        for(int i=1;i<=NL;i++) alphaL[i] = Alpha1(i);
        return 0 ;
    }
    if(n <= NL) {
        return alphaL[n] ;
    } else {
        return Alpha1(n) ;
    }
    
}

// version avec boucles imbriquees
int PB708b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int32_t N2 = (int32_t) Sqrt64(N) ;
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int32_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    int32_t NL = (int32_t)pow(N/4,1/3.0) ;
    Alpha_b(-NL);
    int64_t S = Alpha_b(N) ;
    for(int ip1=0;ip1<nbPrime;ip1++) {
        T_prime64 p1 = tbPrime[ip1] ;
        int64_t pn1 = p1*p1 ;
        int64_t IN1 = N/pn1 ;
        for(int ipowP1=0;IN1>=1;IN1 /=p1, ipowP1++) {
            S += Alpha_b(IN1) << ipowP1 ;
            for(int ip2=ip1+1;ip2<nbPrime;ip2++) {
                T_prime64 p2 = tbPrime[ip2] ;
                int64_t pn2 = p2*p2 ;
                if(pn2 > IN1) break ;
                int64_t IN2 = IN1/pn2 ;
                for(int ipowP2=ipowP1;IN2>=1;IN2 /=p2, ipowP2++) {
                    S += Alpha_b(IN2) << ipowP2 ;
                    for(int ip3=ip2+1;ip3<nbPrime;ip3++) {
                        T_prime64 p3 = tbPrime[ip3] ;
                        int64_t pn3 = p3*p3 ;
                        if(pn3 > IN2) break ;
                        int64_t IN3 = IN2/pn3 ;
                        for(int ipowP3=ipowP2;IN3>=1;IN3 /=p3, ipowP3++) {
                            S += Alpha_b(IN3) << ipowP3 ;
                            for(int ip4=ip3+1;ip4<nbPrime;ip4++) {
                                T_prime64 p4 = tbPrime[ip4] ;
                                int64_t pn4 = p4*p4 ;
                                if(pn4 > IN3) break ;
                                int64_t IN4 = IN3/pn4 ;
                                for(int ipowP4=ipowP3;IN4>=1;IN4 /=p4, ipowP4++) {
                                    S += Alpha_b(IN4) << ipowP4 ;
                                    for(int ip5=ip4+1;ip5<nbPrime;ip5++) {
                                        T_prime64 p5 = tbPrime[ip5] ;
                                        int64_t pn5 = p5*p5 ;
                                        if(pn5 > IN4) break ;
                                        int64_t IN5 = IN4/pn5 ;
                                        for(int ipowP5=ipowP4;IN5>=1;IN5 /=p5, ipowP5++) {
                                            S += Alpha_b(IN5) << ipowP5 ;
                                            for(int ip6=ip5+1;ip6<nbPrime;ip6++) {
                                                T_prime64 p6 = tbPrime[ip6] ;
                                                int64_t pn6 = p6*p6 ;
                                                if(pn6 > IN5) break ;
                                                int64_t IN6 = IN5/pn6 ;
                                                for(int ipowP6=ipowP5;IN6>=1;IN6 /=p6, ipowP6++) {
                                                    S += Alpha_b(IN6) << ipowP6 ;
                                                    for(int ip7=ip6+1;ip7<nbPrime;ip7++) {
                                                        T_prime64 p7 = tbPrime[ip7] ;
                                                        int64_t pn7 = p7*p7 ;
                                                        if(pn7 > IN6) break ;
                                                        int64_t IN7 = IN6/pn7 ;
                                                        for(int ipowP7=ipowP6;IN7>=1;IN7 /=p7, ipowP7++) {
                                                            S += Alpha_b(IN7) << ipowP7 ;
                                                            for(int ip8=ip7+1;ip8<nbPrime;ip8++) {
                                                                T_prime64 p8 = tbPrime[ip8] ;
                                                                int64_t pn8 = p8*p8 ;
                                                                if(pn8 > IN7) break ;
                                                                int64_t IN8 = IN7/pn8 ;
                                                                for(int ipowP8=ipowP7;IN8>=1;IN8 /=p8, ipowP8++) {
                                                                    S += Alpha_b(IN8) << ipowP8 ;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


// version avec fonction recursive
int64_t F708(const T_prime64 *tbPrime,int nbPrime,int64_t N,int32_t NP,int64_t divN) {
    int64_t IN0 = N/divN ;
    int64_t res = Alpha_b(IN0) ;
    int ip ;
    for(ip=0;ip<nbPrime;ip++) {
        int64_t p = tbPrime[ip] ;
        int64_t p2 = p*p ;
        int64_t IN = IN0 ;
        if(p2 > IN) break ;
        IN = IN/p2 ;
        for(int ipowP=0;IN>=1;p2 *=p,IN /=p, ipowP++) {
            res += F708(tbPrime+ip+1,nbPrime-ip-1,N,NP,p2*divN) << ipowP ;
        }
    }
    return res ;
}
int PB708c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int32_t N2 = (int32_t) Sqrt64(N) ;
    printf("%.3fs alphaH computed\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC);
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int32_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    printf("Nb primes=%d N2=%d\n",nbPrime,N2);
    int32_t NL = (int32_t)pow(N/4,1/3.0) ;
    Alpha_b(-NL);
    int64_t S = F708(tbPrime,nbPrime,N,N2,1);
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



