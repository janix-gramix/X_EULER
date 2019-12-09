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

#define M_PI 3.14159265358979323846264338327950288

#include "euler_utils.h"
#include "PB_other.h"

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

/*
int PB681(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t maxA2 = PB681_MAXA *  PB681_MAXA ;
    int64_t  sumS = 0 ;
    int64_t a,b,c,d,db,d0,Sb,Sc,S,A ;
    int32_t maxL = Sqrt64(PB681_MAXA) ;
    for(a=1;a<=maxL ;a++) {
      for(b=a;b<=PB681_MAXA ;b++) {
            Sb= a + b ;
            db = a^b ;
            for(c=b;c<=PB681_MAXA;c++) {
                Sc = Sb + c ;
                d0 = c  + ((db ) & 1) ;
                for(d=d0;d<=PB681_MAXA;d+=2) {
                    if((Sc+d)&1) printf(".");
                    S = (Sc + d)/ 2 ;
                    A = (S-a)*(S-b)*(S-c)*(S-d);
                    if(A<=0) continue ;
                    if(A > maxA2) { continue ;}
                    int32_t Asq = Sqrt64(A);
                    if(A == Asq*Asq) {
                        sumS += 2*S ;
                        printf("%c%04lld %lld+%lld+%lld+%lld S=%lld A=%lld  sumS=%lld\n"
                               ,(a!=b && b!=c && c!=d)?'D':' ',Asq,a,b,c,d,S,A,sumS) ;
                    }
                }
            }
        }
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",sumS);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
*/
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
//#define PB691_LG    5000000
#define PB691_LG    5000000


typedef int32_t DC3_TYPE ;
typedef int     DC3_INDEX ;

void RadixSort(const DC3_INDEX * inInd, DC3_INDEX * outInd,const DC3_TYPE *T, int n,DC3_TYPE nbType )
{ // count occurrences
    DC3_INDEX * count = calloc(nbType+1,sizeof(count[0])) ;
    int i;
    for (i = 0; i < n ; i++)  count[T[inInd[i]]]++;
    DC3_INDEX sumc ;
    for (i = 0, sumc = 0; i <= nbType; i++)  {// exclusive prefix sums
        int32_t tmp = count [i]; count[i] = sumc ; sumc += tmp;
    }
    for ( i = 0; i < n; i++) outInd[count[T[inInd[i]]]++ ] = inInd[i]; // sort
    free(count) ;

}

DC3_INDEX * CountRadixSort_uint16(const DC3_INDEX * inInd, DC3_INDEX * outInd,const uint16_t *T, int n,DC3_TYPE nbType )
{ // count occurrences
    DC3_INDEX * count = calloc(nbType+1,sizeof(count[0])) ;
    int i;
    for (i = 0; i < n ; i++)  count[T[inInd[i]]]++;
    DC3_INDEX sumc ;
    for (i = 0, sumc = 0; i <= nbType; i++)  {// exclusive prefix sums
        int32_t tmp = count [i]; count[i] = sumc ; sumc += tmp;
    }
    for ( i = 0; i < n; i++) outInd[count[T[inInd[i]]]++ ] = inInd[i]; // sort
    return count ;
    
}


typedef  int32_t T3[3] ;
typedef struct TT3 {
    int32_t V[3] ;
} TT3 ;

void SuffixSort(const DC3_TYPE * T, DC3_INDEX * SA, int n, DC3_TYPE nbType) {
    DC3_INDEX n2=n/3, n0=(n+2)/3, d1 = 0 , n02=n0+n2;
    DC3_INDEX * R = malloc((n02+3)*sizeof(R[0])) ; R[n02]=R[n02+1]=R[n02+2]=0;
    DC3_INDEX * SA12 = malloc((n02+3)*sizeof(SA12[0])); SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
    DC3_INDEX * R0 = malloc((n0+1)*sizeof(R0[0])) ;
    DC3_INDEX * SA0 = malloc((n0+1)*sizeof(SA0[0]));
    //******* Step 0: Construct sample ********
    // generate positions of mod 1 and mod 2 suffixes
    int i ; // R contains n02 index
    for(i=0;i< n2;i++) {   R[2*i] = 3*i+1;    R[2*i+1] = 3*i+2; }
    if (n == 3*n2+1) { R[2*n2] = n ;  d1=1 ; } // +1 (size of T)
    else if (n == 3*n2+2) R[2*n2] = n-1 ;
//---------- Sort sample suffixes
    // lsb radix sort the mod 1 and mod 2 triples in C = B1 u B2 (samples suffixes)
    RadixSort(R , SA12, T+2, n02, nbType); // +2 (size ofT) (with +1)=> +3
    RadixSort(SA12, R , T+1, n02, nbType);
    RadixSort(R , SA12, T , n02, nbType);
 
// find different triples and write them in R
    DC3_INDEX nb3uple = 0;
    DC3_TYPE T3_0[3] = {-1,-1,-1};
    const DC3_TYPE *antT3 =T3_0  ;
    for (int i = 0; i < n02; i++) {
        int ind = SA12[i] ;
        const DC3_TYPE  * T3 = T+ind ;
        if(T3[0] != antT3[0] || T3[1] != antT3[1] || T3[2] != antT3[2]) {
            nb3uple++;   antT3 = T3 ; // new triple
        }
        int ind_3 = ind/3 ; // recall, in SA12 ind = 1 or 2 mod[3]
        if (ind - 3*ind_3 == 1) { R[ind_3] = nb3uple; } // ind=1 mod[3]
        else { R[ind_3+n0] = nb3uple; } // ind=2 mod[3]
    }

    if (nb3uple < n02) { // triple not unique => recursion
        SuffixSort(R, SA12, n02, nb3uple);
        // compute reverse numerotation in R (add 1 as num in SA12 as +1)
        for (int i = 0; i < n02; i++) R[SA12[i]] = i + 1;
    } else {// all suffix in 12 are distinct, remove 1 to begin num T3 at 0
        for (int i = 0; i < n02; i++) SA12[R[i] - 1] = i;
    }
// ---------- Sort nonsample suffixes
    // radix sort the mod[3]= 0 suffixes from SA12
    for (int i=0, j=0; i < n02; i++) if (SA12[i] < n0) R0[j++] = 3*SA12[i];
    RadixSort(R0, SA0, T, n0, nbType);
//---------- Merge
    // merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p=0, t=d1, k=0; k < n; k++) {
// #define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
 //       int i = GetI(); // pos of current offset 12 suffix
        DC3_INDEX i_t, i_p = SA0[p]; // pos of current offset 0 suffix
        if((i_t = SA12[t] ) < n0) {
            DC3_INDEX i_3t = i_t * 3 + 1  ;
            //   Comparison (3t+1,R[t + n0]) and (p,R[p/3]
             if(T[i_3t] < T[i_p] || ( T[i_3t] == T[i_p] && R[i_t + n0] <= R[i_p/3] )) {
                SA[k] = i_3t ; t++;
                // add remaining  SA0 suffixes
                if (t == n02)  while(p < n0) SA[++k] = SA0[p++];
            } else {
                SA[k] = i_p; p++;
                // add remaining SA12 suffixes
                if (p == n0)  while(t < n02){ SA[++k] = ((i_t= SA12[t++]) < n0) ? (3*i_t+1) : (3*(i_t - n0)+2) ; }
            }
        } else {
           DC3_INDEX i_t3 =(i_t - n0) * 3 + 2 ;
            //  Comparison (3t+2,3t+3,R[t-n0+1]) and (p,p+1,R[p/3+n0]
            if(  T[i_t3] < T[i_p] || ( T[i_t3]==T[i_p]
                 &&( T[i_t3+1] < T[i_p+1] || (T[i_t3+1]==T[i_p+1] && R[i_t-n0+1] <=  R[i_p/3+n0]) )))
            {
                SA[k] = i_t3; t++;
                // add remaining  SA0 suffixes
                if (t == n02)  while(p < n0) SA[++k] = SA0[p++];
            } else {
                SA[k] = i_p; p++;
                // add remaining SA12 suffixes
                if (p == n0)  while(t < n02){ SA[++k] = ((i_t= SA12[t++]) < n0) ? (3*i_t+1) : (3*(i_t - n0)+2) ; }
            }
        }
    }
    free(R); free(SA12); free(SA0); free(R0);
}

static int AtomLCP(const uint8_t * text,DC3_INDEX len,DC3_INDEX ia,DC3_INDEX ib) {
    if(ia<=ib) { int32_t tmp = ia ; ia= ib ; ib =tmp ;    }
    int l = 0 ;
    while(ia < len && text[ia++] == text[ib++]) l++ ;
    return l ;
}
DC3_INDEX * GetLCP(const uint8_t * text,const DC3_INDEX * Sindex, int nbIndex,int len) {
    DC3_INDEX * lcp = malloc(nbIndex*sizeof(lcp[0])) ;
    lcp[0] = 0;
    int32_t *Sorder=malloc(nbIndex * sizeof(Sorder[0])) ; // inverse permutation of index
    for (int i = 0; i < nbIndex; i++) Sorder[Sindex[i]] = i; //inverse permutation
    if (Sorder[0]) { // chaine initiale, partie commune avec la chaine totale
        lcp[Sorder[0]] = AtomLCP(text,len,0,Sindex[Sorder[0] - 1]);
    }
    for (int i = 1; i < nbIndex; i++) {
        if (!Sorder[i]) continue; // rank=0
        if (lcp[Sorder[i - 1]] <= 1) { // si la chaine precedente n'a rien de commun calcul direct
            lcp[Sorder[i]] = AtomLCP(text ,len, i, Sindex[Sorder[i] - 1]);
        } else { // si la chaine precedente a un prefixe commun avec celle d'avant on regarde si l'on peut prolonger
            int L = lcp[Sorder[i - 1]] - 1;
           lcp[Sorder[i]] = L + AtomLCP(text,len, i + L,  Sindex[Sorder[i] - 1] + L);
        }
    }
    free(Sorder);
    return lcp ;
}



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
    DC3_TYPE *suit_dc3 = malloc((PB691_LG+3)*sizeof(suit_dc3[0])) ;
    for(i=0;i<n;i++) suit_dc3[i] = suit[i] + 1 ;  suit_dc3[n] = suit_dc3[n+1] = suit_dc3[n+2] = 0 ; //
    DC3_INDEX * sa= malloc((n+1)*sizeof(sa[0]));
    SuffixSort(suit_dc3, sa, n, 2);
    free(suit_dc3 ); // not necessary after
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


#define MAX_CHAR 3

uint8_t * suitUK ;
int suitUKLen ;

struct SuffixTreeNode {
    struct SuffixTreeNode *children[MAX_CHAR];
    
    //pointer to other node via suffix link
    struct SuffixTreeNode *suffixLink;
    
    /*(start, end) interval specifies the edge, by which the
     node is connected to its parent node. Each edge will
     connect two nodes,  one parent and one child, and
     (start, end) interval of a given edge  will be stored
     in the child node. Lets say there are two nods A and B
     connected by an edge with indices (5, 8) then this
     indices (5, 8) will be stored in node B. */
    int start;
    int *end;
    
    /*for leaf nodes, it stores the index of suffix for
     the path  from root to leaf*/
    int suffixIndex;
};

typedef struct SuffixTreeNode Node;

Node *root = NULL; //Pointer to root node

/*lastNewNode will point to newly created internal node,
 waiting for it's suffix link to be set, which might get
 a new suffix link (other than root) in next extension of
 same phase. lastNewNode will be set to NULL when last
 newly created internal node (if there is any) got it's
 suffix link reset to new internal node created in next
 extension of same phase. */
Node *lastNewNode = NULL;
Node *activeNode = NULL;

/*activeEdge is represeted as input string character
 index (not the character itself)*/
int activeEdge = -1;
int activeLength = 0;

// remainingSuffixCount tells how many suffixes yet to
// be added in tree
int remainingSuffixCount = 0;
int leafEnd = -1;
int *rootEnd = NULL;
int *splitEnd = NULL;
int size = -1; //Length of input string

Node *newNode(int start, int *end)
{
    Node *node =(Node*) malloc(sizeof(Node));
    int i;
    for (i = 0; i < MAX_CHAR; i++)
        node->children[i] = NULL;
    
    /*For root node, suffixLink will be set to NULL
     For internal nodes, suffixLink will be set to root
     by default in  current extension and may change in
     next extension*/
    node->suffixLink = root;
    node->start = start;
    node->end = end;
    
    /*suffixIndex will be set to -1 by default and
     actual suffix index will be set later for leaves
     at the end of all phases*/
    node->suffixIndex = -1;
    return node;
}
int edgeLength(Node *n) {
    return *(n->end) - (n->start) + 1;
}

int walkDown(Node *currNode)
{
    /*activePoint change for walk down (APCFWD) using
     Skip/Count Trick  (Trick 1). If activeLength is greater
     than current edge length, set next  internal node as
     activeNode and adjust activeEdge and activeLength
     accordingly to represent same activePoint*/
    if (activeLength >= edgeLength(currNode))
    {
        activeEdge += edgeLength(currNode);
        activeLength -= edgeLength(currNode);
        activeNode = currNode;
        return 1;
    }
    return 0;
}
void extendSuffixTree(int pos)
{
    /*Extension Rule 1, this takes care of extending all
     leaves created so far in tree*/
    leafEnd = pos;
    
    /*Increment remainingSuffixCount indicating that a
     new suffix added to the list of suffixes yet to be
     added in tree*/
    remainingSuffixCount++;
    
    /*set lastNewNode to NULL while starting a new phase,
     indicating there is no internal node waiting for
     it's suffix link reset in current phase*/
    lastNewNode = NULL;
    
    //Add all suffixes (yet to be added) one by one in tree
    while(remainingSuffixCount > 0) {
        
        if (activeLength == 0)
            activeEdge = pos; //APCFALZ
        
        // There is no outgoing edge starting with
        // activeEdge from activeNode
        if (activeNode->children[suitUK[activeEdge]] == NULL)
        {
            //Extension Rule 2 (A new leaf edge gets created)
            activeNode->children[suitUK[activeEdge]] = newNode(pos, &leafEnd);
            
            /*A new leaf edge is created in above line starting
             from  an existng node (the current activeNode), and
             if there is any internal node waiting for it's suffix
             link get reset, point the suffix link from that last
             internal node to current activeNode. Then set lastNewNode
             to NULL indicating no more node waiting for suffix link
             reset.*/
            if (lastNewNode != NULL)
            {
                lastNewNode->suffixLink = activeNode;
                lastNewNode = NULL;
            }
        }
        // There is an outgoing edge starting with activeEdge
        // from activeNode
        else
        {
            // Get the next node at the end of edge starting
            // with activeEdge
            Node *next = activeNode->children[suitUK[activeEdge]];
            if (walkDown(next))//Do walkdown
            {
                //Start from next node (the new activeNode)
                continue;
            }
            /*Extension Rule 3 (current character being processed
             is already on the edge)*/
            if (suitUK[next->start + activeLength] == suitUK[pos])
            {
                //If a newly created node waiting for it's
                //suffix link to be set, then set suffix link
                //of that waiting node to current active node
                if(lastNewNode != NULL && activeNode != root)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = NULL;
                }
                
                //APCFER3
                activeLength++;
                /*STOP all further processing in this phase
                 and move on to next phase*/
                break;
            }
            
            /*We will be here when activePoint is in middle of
             the edge being traversed and current character
             being processed is not  on the edge (we fall off
             the tree). In this case, we add a new internal node
             and a new leaf edge going out of that new node. This
             is Extension Rule 2, where a new leaf edge and a new
             internal node get created*/
            splitEnd = (int*) malloc(sizeof(int));
            *splitEnd = next->start + activeLength - 1;
            
            //New internal node
            Node *split = newNode(next->start, splitEnd);
            activeNode->children[suitUK[activeEdge]] = split;
            
            //New leaf coming out of new internal node
            split->children[suitUK[pos]] = newNode(pos, &leafEnd);
            next->start += activeLength;
            split->children[suitUK[next->start]] = next;
            
            /*We got a new internal node here. If there is any
             internal node created in last extensions of same
             phase which is still waiting for it's suffix link
             reset, do it now.*/
            if (lastNewNode != NULL)
            {
                /*suffixLink of lastNewNode points to current newly
                 created internal node*/
                lastNewNode->suffixLink = split;
            }
            
            /*Make the current newly created internal node waiting
             for it's suffix link reset (which is pointing to root
             at present). If we come across any other internal node
             (existing or newly created) in next extension of same
             phase, when a new leaf edge gets added (i.e. when
             Extension Rule 2 applies is any of the next extension
             of same phase) at that point, suffixLink of this node
             will point to that internal node.*/
            lastNewNode = split;
        }
        
        /* One suffix got added in tree, decrement the count of
         suffixes yet to be added.*/
        remainingSuffixCount--;
        if (activeNode == root && activeLength > 0) //APCFER2C1
        {
            activeLength--;
            activeEdge = pos - remainingSuffixCount + 1;
        }
        else if (activeNode != root) //APCFER2C2
        {
            activeNode = activeNode->suffixLink;
        }
    }
}

void print(int i, int j)
{
    int k;
    for (k=i; k<=j; k++)
        printf("%d", suitUK[k]-1);
}

//Print the suffix tree as well along with setting suffix index
//So tree will be printed in DFS manner
//Each edge along with it's suffix index will be printed
void setSuffixIndexByDFS(Node *n, int labelHeight)
{
    if (n == NULL)  return;
    
    if (n->start != -1) //A non-root node
    {
        //Print the label on edge from parent to current node
//        print(n->start, *(n->end));
    }
    int leaf = 1;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
 //           if (leaf == 1 && n->start != -1)printf(" [%d]\n", n->suffixIndex);
            
            //Current node is not a leaf as it has outgoing
            //edges from it.
            leaf = 0;
            setSuffixIndexByDFS(n->children[i], labelHeight +
                                edgeLength(n->children[i]));
        }
    }
    if (leaf == 1)
    {
        n->suffixIndex = size - labelHeight;
 //       printf(" [%d]\n", n->suffixIndex);
    }
}
void freeSuffixTreeByPostOrder(Node *n)
{
    if (n == NULL)
        return;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
            freeSuffixTreeByPostOrder(n->children[i]);
        }
    }
    if (n->suffixIndex == -1)
        free(n->end);
    free(n);
}

/*Build the suffix tree and print the edge labels along with
 suffixIndex. suffixIndex for leaf edges will be >= 0 and
 for non-leaf edges will be -1*/
void buildSuffixTree()
{
    size = suitUKLen ;
    int i;
    rootEnd = (int*) malloc(sizeof(int));
    *rootEnd = - 1;
    
    /*Root is a special node with start and end indices as -1,
     as it has no parent from where an edge comes to root*/
    root = newNode(-1, rootEnd);
    
    activeNode = root; //First activeNode will be root
    for (i=0; i<size; i++)
        extendSuffixTree(i);
    int labelHeight = 0;
    setSuffixIndexByDFS(root, labelHeight);
    
    //Free the dynamically allocated memory
//    freeSuffixTreeByPostOrder(root);
}

void doTraversal(Node *n, int suffixArray[], int *idx)
{
    if(n == NULL)
    {
        return;
    }
    int i=0;
    if(n->suffixIndex == -1) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != NULL)
            {
                doTraversal(n->children[i], suffixArray, idx);
            }
        }
    }
    //If it is Leaf node other than "$" label
    else if(n->suffixIndex > -1 && n->suffixIndex < size)
    {
        suffixArray[(*idx)++] = n->suffixIndex;
    }
}

void buildSuffixArray(int suffixArray[])
{
    int i = 0;
    for(i=0; i< size; i++)suffixArray[i] = -1;
    int idx = 0;
    doTraversal(root, suffixArray, &idx);
/*    printf("Suffix Array for String ");
    for(i=0; i<size; i++)printf("%d", suitUK[i]);
    printf(" is: \n");
    for(i=1; i<size; i++) {
        int j ;
        for(j=suffixArray[i];j<suitUKLen-1;j++) printf("%d",suitUK[j]-1) ;
       
        printf("\n");
    }
*/
}

int PB691c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i ;
    uint8_t * suit = InitSuit(PB691_LG+1) ;
    suitUK = suit ;
    int n = PB691_LG ;
    for(i=0;i<n;i++) suitUK[i] += 1 ;
    suit[n] = 0 ;
    
 //   memcpy(suitUK,"aabbaabbabb",11);
 //   n = 11 ;
    suitUKLen = n+1 ;
    buildSuffixTree();
    int *sa =(int*) malloc(n*sizeof(sa[0]));
    buildSuffixArray(sa);
  
    DC3_INDEX * lcp = GetLCP(suit, sa+1,n,n);
    free(sa); free(suit) ; // not necessary after
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
