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
    int i  ;
    uint64_t cnp = 1 ;
    for(i=0;i<p;i++) {
        cnp *= n-i ;
        cnp /= i+1 ;
    }
    return cnp ;
}
#define INDT2(n0,n1)  (((n0)+(n1))*((n0)+(n1)+1)/2+(n1))

typedef struct C3cut {
    int nbInsert ;
    int nbUseP ;
    int fact ;
    int deltaIndex[2] ;
} C3cut ;

typedef struct C3x {
    int nbPieces ;
    int nbCut ;
    const C3cut   *tbCut ;
} C3x ;



int Index3(int S,int i0,int i1) {
    static int Cn1[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91} ;
    int ind = Cn1[S-i0]+S-i0-i1 ;
    return ind ;
}


int Index5(int S,int i0,int i1,int i2,int i3) {
    static int Cn3[] = {0,1,5,15,35,70,126,210,330,495,715,1001,1365,1820} ;
    static int Cn2[] = {0,1,4,10,20,35,56,84,120,165,220,286,364,455} ;
    static int Cn1[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91} ;
    int ind = Cn3[S-i0]+Cn2[S-i0-i1]+Cn1[S-i0-i1-i2]+S-i0-i1-i2-i3 ;
    return ind ;
}


typedef struct C3i {
    int nc0 ;       // nombre de cesure diff de niveau 0
    int isNc1 ;     // y-a-t-l des cesure de niveau 1
    int deltaIndex[2] ;
} C3i ;

// #define SetIndex5(NR,N0,N1,N2,N3,val) printf("(%d,%d,%d,%d,%d)+=%lld ",N0,N1,N2,N3,NR-(N0)-(N1)-(N2)-(N3),(val)); cur5[Index5(NR,N0,N1,N2,N3)] += (val)

#define SetIndex5(NR,N0,N1,N2,N3,val)  cur5[Index5(NR,N0,N1,N2,N3)] += (val)

#define FLOAT687
#if defined(FLOAT687)
    typedef  double PROB687 ;
#   define MAX_NR   13
#else
    typedef  uint64_t PROB687 ;
#   define MAX_NR   7
#endif
// insertion ****
#define INS4_N0(val) SetIndex5(nr,n0,n1,n2,n3,(val)*nb0)
#define INS4_N1(val) if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3,(val)*n1)
#define INS4_N2(val) if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3,(val)*2*n2)
#define INS4_N3(val) if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3-1,(val)*2*n3)
#define INS4_N4_2(val) if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3,(val)*n4)
#define INS4_N4_3(val) if(n4>0)SetIndex5(nr,n0,n1,n2,n3+1,(val)*2*n4)

// insertion *** *
// *** * 0 x X different ranks
#define INS3_N0x0(val) SetIndex5(nr,n0,n1,n2,n3+1,(val)*nb0*(nb0-1))
#define INS3_N0x1(val) if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3+1,(val)*2*nb0*n1)
#define INS3_N0x2(val) if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3+1,(val)*2*nb0*2*n2)
#define INS3_N0x3(val) if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3,(val)*2*nb0*2*n3)
#define INS3_N0x4_2(val) if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3+1,(val)*2*nb0*n4)
#define INS3_N0x4_3(val) if(n4>0)SetIndex5(nr,n0,n1,n2,n3+2,(val)*2*nb0*2*n4)
// *** * 1 x X different ranks
#define INS3_N1x1(val) if(n1>1)SetIndex5(nr,n0+2,n1-2,n2,n3+1,(val)*n1*(n1-1))
#define INS3_N1x2(val) if(n1>0 && n2>0)SetIndex5(nr,n0+1,n1,n2-1,n3+1,(val)*2*n1*2*n2)
#define INS3_N1x3(val) if(n1>0 && n3>0)SetIndex5(nr,n0+1,n1,n2,n3,(val)*2*n1*2*n3)
#define INS3_N1x4_2(val) if(n1>0 && n4>0)SetIndex5(nr,n0+1,n1-1,n2+1,n3+1,(val)*2*n1*n4)
#define INS3_N1x4_3(val) if(n1>0 && n4>0)SetIndex5(nr,n0+1,n1-1,n2,n3+2,(val)*2*n1*2*n4)
// *** * 2 x X different ranks
#define INS3_N2x2(val) if(n2>1)SetIndex5(nr,n0,n1+2,n2-2,n3+1,(val)*2*n2*2*(n2-1))
#define INS3_N2x3(val) if(n2>0 && n3>0)SetIndex5(nr,n0,n1+2,n2-1,n3,(val)*2*2*n2*2*n3)
#define INS3_N2x4_2(val) if(n2>0 && n4>0)SetIndex5(nr,n0,n1+1,n2,n3+1,(val)*2*2*n2*n4)
#define INS3_N2x4_3(val) if(n2>0 && n4>0)SetIndex5(nr,n0,n1+1,n2-1,n3+2,(val)*2*2*n2*2*n4)
// *** * 3 x X different ranks
#define INS3_N3x3(val) if(n3>1)SetIndex5(nr,n0,n1+2,n2,n3-1,(val)*2*n3*2*(n3-1))
#define INS3_N3x4_2(val) if(n3>0 && n4>0)SetIndex5(nr,n0,n1+1,n2+1,n3,(val)*2*2*n3*n4)
#define INS3_N3x4_3(val) if(n3>0 && n4>0)SetIndex5(nr,n0,n1+1,n2,n3+1,(val)*2*2*n3*2*n4)
// *** * 4 x 4 different ranks
#define INS3_N4_2x4_2(val) if(n4>1)SetIndex5(nr,n0,n1,n2+2,n3+1,(val)*n4*(n4-1))
#define INS3_N4_2x4_3(val) if(n4>1)SetIndex5(nr,n0,n1,n2+1,n3+2,(val)*2*n4*2*(n4-1))
#define INS3_N4_3x4_3(val) if(n4>1)SetIndex5(nr,n0,n1,n2,n3+3,(val)*n4*2*2*(n4-1))
// *** * same rank
#define IND3_N2(val) if(n2>0)SetIndex5(nr,n0+1,n1,n2-1,n3+1,(val)*2*n2)
#define IND3_N3(val) if(n3>0)SetIndex5(nr,n0+1,n1,n2,n3,(val)*2*n3)
#define IND3_N4(val) if(n4>0)SetIndex5(nr,n0,n1+1,n2,n3+1,(val)*3*2*n4)

// insertion ** **
// ** ** 0 x X different ranks
#define INS2_N0x0(val) SetIndex5(nr,n0,n1,n2+1,n3,(val)*nb0*(nb0-1)/2)
#define INS2_N0x1(val) if(n1>0)SetIndex5(nr,n0+1,n1-1,n2+1,n3,(val)*nb0*n1)
#define INS2_N0x2(val) if(n2>0)SetIndex5(nr,n0,n1+1,n2,n3,(val)*nb0*2*n2)
#define INS2_N0x3(val) if(n3>0)SetIndex5(nr,n0,n1+1,n2+1,n3-1,(val)*nb0*2*n3)
#define INS2_N0x4_2(val) if(n4>0)SetIndex5(nr,n0,n1,n2+2,n3,(val)*nb0*n4)
#define INS2_N0x4_3(val) if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3+1,(val)*nb0*2*n4)
// ** ** 1 x X different ranks
#define INS2_N1x1(val) if(n1>1)SetIndex5(nr,n0+2,n1-2,n2+1,n3,(val)*n1*(n1-1)/2)
#define INS2_N1x2(val) if(n1>0 && n2>0)SetIndex5(nr,n0+1,n1,n2,n3,(val)*n1*2*n2)
#define INS2_N1x3(val) if(n1>0 && n3>0)SetIndex5(nr,n0+1,n1,n2+1,n3-1,(val)*n1*2*n3)
#define INS2_N1x4_2(val) if(n1>0 && n4>0)SetIndex5(nr,n0+1,n1-1,n2+2,n3,(val)*n1*n4)
#define INS2_N1x4_3(val) if(n1>0 && n4>0)SetIndex5(nr,n0+1,n1-1,n2+1,n3+1,(val)*n1*2*n4)
// ** ** 2 x X different ranks
#define INS2_N2x2(val) if(n2>1)SetIndex5(nr,n0,n1+2,n2-1,n3,(val)*n2*2*(n2-1))
#define INS2_N2x3(val) if(n2>0 && n3>0)SetIndex5(nr,n0,n1+2,n2,n3-1,(val)*2*n2*2*n3)
#define INS2_N2x4_2(val) if(n2>0 && n4>0)SetIndex5(nr,n0,n1+1,n2+1,n3,(val)*2*n2*n4)
#define INS2_N2x4_3(val) if(n2>0 && n4>0)SetIndex5(nr,n0,n1+1,n2,n3+1,(val)*2*n2*2*n4)
// ** ** 3 x X different ranks
#define INS2_N3x3(val) if(n3>1)SetIndex5(nr,n0,n1+2,n2+1,n3-2,(val)*n3*2*(n3-1))
#define INS2_N3x4_2(val) if(n3>0 && n4>0)SetIndex5(nr,n0,n1+1,n2+2,n3-1,(val)*2*n3*n4)
#define INS2_N3x4_3(val) if(n3>0 && n4>0)SetIndex5(nr,n0,n1+1,n2+1,n3,(val)*2*n3*2*n4)
// ** ** 4 x 4 different ranks
#define INS2_N4_2x4_2(val) if(n4>1)SetIndex5(nr,n0,n1,n2+3,n3,(val)*n4*(n4-1)/2)
#define INS2_N4_2x4_3(val) if(n4>1)SetIndex5(nr,n0,n1,n2+2,n3+1,(val)*n4*2*(n4-1))
#define INS2_N4_3x4_3(val) if(n4>1)SetIndex5(nr,n0,n1,n2+1,n3+2,(val)*n4*2*(n4-1))
// ** ** same ranks
#define IND2_N2(val) if(n2>0)SetIndex5(nr,n0+1,n1,n2,n3,(val)*n2)
#define IND2_N3(val) if(n3>0)SetIndex5(nr,n0+1,n1,n2+1,n3-1,(val)*n3)
#define IND2_N4(val) if(n4>0)SetIndex5(nr,n0,n1+1,n2+1,n3,(val)*3*n4)


int PB687a(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nr,nc ;

    PROB687 ant5[4000] ;
    PROB687 cur5[4000] ;
    
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
                        //                          int nb0 = 4*n0+3*n1+2*n2+2*n3+n4+1 ;//
                        int nb0 = 4*n0+3*n1+2*n2+2*n3+n4 ;
                        PROB687 ant = ant5[Index5(nr-1,n0,n1,n2,n3)] ;
                        if(ant <= 0) continue ;
                           // insert ****
                        INS4_N0(ant) ; INS4_N1(ant) ; INS4_N2(ant) ; INS4_N3(ant) ; INS4_N4_2(ant) ; INS4_N4_3(ant) ;
   /*                     // insert *** *
                        INS3_N0x0(ant) ;  INS3_N0x1(ant) ;  INS3_N0x2(ant) ; INS3_N0x3(ant) ; INS3_N0x4_2(ant) ; INS3_N0x4_3(ant);
                        INS3_N1x1(ant) ;  INS3_N1x2(ant) ; INS3_N1x3(ant) ; INS3_N1x4_2(ant) ; INS3_N1x4_3(ant) ;
                        INS3_N2x2(ant) ; INS3_N2x3(ant) ; INS3_N2x4_2(ant) ; INS3_N2x4_3(ant) ;
                        INS3_N3x3(ant) ; INS3_N3x4_2(ant) ; INS3_N3x4_3(ant) ;
                        INS3_N4_2x4_2(ant) ; INS3_N4_2x4_3(ant) ; INS3_N4_3x4_3(ant) ;
                        IND3_N2(ant) ;  IND3_N3(ant) ; IND3_N4(ant) ;
                        // insert ** **
                        INS2_N0x0(ant) ;  INS2_N0x1(ant) ;  INS2_N0x2(ant) ; INS2_N0x3(ant) ; INS2_N0x4_2(ant) ; INS2_N0x4_3(ant);
                        INS2_N1x1(ant) ;  INS2_N1x2(ant) ; INS2_N1x3(ant) ; INS2_N1x4_2(ant) ; INS2_N1x4_3(ant) ;
                        INS2_N2x2(ant) ; INS2_N2x3(ant) ; INS2_N2x4_2(ant) ; INS2_N2x4_3(ant) ;
                        INS2_N3x3(ant) ; INS2_N3x4_2(ant) ; INS2_N3x4_3(ant) ;
                        INS2_N4_2x4_2(ant) ; INS2_N4_2x4_3(ant) ; INS2_N4_3x4_3(ant) ;
                        IND2_N2(ant) ;  IND2_N3(ant) ; IND2_N4(ant) ;
     */
                        
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
//                        for(nbp=3;nbp<=4;nbp++) {
 //                           int d0 = (nbp == 4) ? 1 : 0 ;
                            for(i0=0;i0<=nbp && i0 <= nb0 ;i0++){
                                int64_t delta0 =  Cnp(nb0,i0) ;
                                if(d1 ) delta0 *= 3 ;
                                else if(d3) delta0 *= 2 ;
                                for(i1=0; i1<=nbp-i0 && i1<=n1;i1++){
                                    int64_t delta1 = delta0 * Cnp(n1,i1) ;
                                    int i2max = nbp-i0-i1 ;
                                    if(i2max==0) { SetIndex5( nr,n0+i1+d0,n1-i1+d1,n2+d2,n3+d3,ant*delta1) ; break ;}
                                    for(i22=0; i22 <= i2max/2 && i22<=n2;i22++){
                                        for(i21=0;i21<=(i2max-2*i22) && i21+i22 <= n2;i21++){
                                            int64_t delta2 = delta1 * Cnp(n2,i22) * Cnp(n2-i22,i21)*(1<<i21);
                                            int i3max = i2max-2*i22 - i21 ;
                                            if(i3max==0) {
                                                SetIndex5( nr,n0+d0+i1+i22,n1+d1-i1+i21,n2+d2-i22-i21,n3+d3,ant * delta2) ;
                                                break ;
                                            }
                                            for(i32=0; i32<=i3max/2 && i32 <= n3;i32++){
                                                for(i31=0; i31 <= i3max-2*i32 && i31+i32 <=n3;i31++){
                                                    int64_t delta3 = delta2 * Cnp(n3,i32) * Cnp(n3-i32,i31)*(1<<i31) ;
                                                    int i4max = i3max -2 * i32 - i31 ;
                                                    if(i4max == 0 ) { SetIndex5( nr,n0+d0+i1+i22+i32,n1+d1-i1+i21+i31,n2+d2-i22-i21,n3+d3-i32-i31,ant * delta3) ;break ; }
                                                    for(i43=0;i43<= (i4max)/3 && i43 <=n4 ;i43++){
                                                        for(i42=0;i42<=(i4max-3*i43)/2 && i43+i42<=n4;i42++){
                                                            int i41 = i4max-3*i43 -2*i42 ;
                                                            if(i41+i42+i43 >n4) continue ;
                                                            int64_t delta = delta3 ;
                                                            if(i43) delta *= n4 ;
                                                            if(i42) delta *= Cnp(n4-i43,i42)*pow3[i42]  ;
                                                            int j3 ;
                                                            for(j3=0;j3<=i41;j3++) {
                                                                SetIndex5( nr,n0+d0+i1+i22+i32+i43 ,n1+d1-i1+i21+i31+i42  ,n2+d2-i22-i21+i41-j3 ,n3+d3-i32-i31+j3,
                                                                          ant * delta * Cnp(n4-i43-i42,i41) * (1<<j3) * Cnp(i41,j3));
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
#if defined(FLOAT687)
        PROB687 g = 1 ;
#else
        PROB687 g =cur5[0];
        for(i=1;i<=Index5(nr,0,0,0,0);i++) { g = PGCD64(g,cur5[i]) ;} ;
#endif
        PROB687 den = 0;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) { /*cur5[i] /= g;*/ den += cur5[i]/g ; }
        PROB687 prob[20];
        for(i=0;i<=nr;i++)prob[i]=0 ;
        for(n0=0;n0<=nr;n0++) {
            for(n1=0;n1<=nr-n0;n1++) {
                for(n2=0;n2<=nr-n0-n1;n2++) {
                    for(n3=0;n3<=nr-n0-n1-n2;n3++) {
                        //                                 printf("(%d,%d,%d,%d,%d)=%lld ",n0,n1,n2,n3,nr-n0-n1-n2-n3,cur5[Index5(nr,n0,n1,n2,n3)]);
                        cur5[Index5(nr,n0,n1,n2,n3)] /= g ;
                        prob[n0]+= cur5[Index5(nr,n0,n1,n2,n3)] ;
                    }
                }
            }
        }
        
        printf("\n4 color, %d rank  => nb index=%d",nr,Index5(nr,0,0,0,0));
        printf("\n");
#if defined(FLOAT687)
        for(i=0;i<=nr;i++) { printf("%.11f ",prob[i]/(double)den) ; } ;
#else
        for(i=0;i<=nr;i++) { printf("%lld ",prob[i]) ; } ;
#endif
        if(nr==MAX_NR) {
            printf("\n answer=%.11f \n", (double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13])/den) ;
            snprintf(pbR->strRes, sizeof(pbR->strRes),"%.10f",(double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13])/den);
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
#if defined(FLOAT687)
        for(i=0;i<=Index5(nr,0,0,0,0);i++) ant5[i] = cur5[i]/den ;
#else
        for(i=0;i<=Index5(nr,0,0,0,0);i++) ant5[i] = cur5[i] ;
#endif
    } // end loop nr
    printf("\n*********\n") ;
    

    return 0;
}

int PB687(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nr,nc ;
    {
        PROB687 ant5[4000] ;
        PROB687 cur5[4000] ;
        
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
  //                          int nb0 = 4*n0+3*n1+2*n2+2*n3+n4+1 ;//
                            int nb0 = 4*n0+3*n1+2*n2+2*n3+n4 ;
                            PROB687 ant = ant5[Index5(nr-1,n0,n1,n2,n3)] ;
                            if(ant <= 0) continue ;
  //                          printf("\n(%d,%d,%d,%d,%d)->%lld ",n0,n1,n2,n3,nr-1-n0-n1-n2-n3,ant);
                            // ****
  //                          printf("\n Insert **** ") ;
                            int64_t S = 0;
                            SetIndex5(nr,n0,n1,n2,n3, ant * nb0) ; S+=  nb0 ; // (0,0,+1) ]****[
                            if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3, ant * n1) ; S+= n1 ;// (+1,-1,0,0,+1) [****]
                            if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3, ant * 2 * n2) ; S+=  2 * n2 ;// (0,+1,-1,0,+1) [[****]]
                            if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3-1, ant * 2 * n3) ; S+= 2*n3 ; // (0,+1,0,-1,+1) [[[****]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2,n3+1, ant * 2 * n4) ;S+= 2 * n4 ; // (0,0,0,+1,-1) [[[[****]]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3, ant  * n4); S+= n4 ; // (0,0,+1,0,-1) [[[[****]]]]
   //                         printf("S4=%lld ",S);
                           // *** *
   //                         printf("\n Insert *** * ") ;

                            S = 0 ;
                                    // (0,0,0,+1,0) ]*[ ]***[
                            SetIndex5(nr,n0,n1,n2,n3+1, ant * nb0 * (nb0-1)); S+= nb0 * (nb0-1);
                                    // (+1,-1,0,+1,0) ]*[ [***] || ]***[ [*]
                            if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3+1, ant * 2 * nb0 * n1) ; S+= 2 * nb0 * n1;
                                    // (0,+1,-1,+1,0) ]*[ [[***]] || ]***[ [[*]]
                            if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3+1, ant * 2 * nb0 * 2 * n2) ;S+= 2 * nb0 * 2 * n2 ;
                                    // (0,+1,0,0,0) ]*[ [[[***]]] || ]***[ [[[*]]]
                            if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3, ant * 2 * nb0 * 2 * n3 );S+= 2 * nb0 * 2 * n3 ;
                                    // (0,0,0,+2,0) ]*[ [[[[***]]]] || ]***[ [[[[*]]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2,n3+2, ant * 2 * nb0 * 2 * n4 );S+= 2 * nb0 * 2 * n4;
                                    // (0,+1,0,+1,0) ]*[ [[[[***]]]] || ]***[ [[[[*]]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3+1, ant * 2 * nb0  * n4 );S+= 2 * nb0  * n4;
                                // (+2,-2,0,+1,0) [*] [***] || [***] [*]
                            if(n1 > 1 ) SetIndex5(nr,n0+2,n1-2,n2,n3+1, ant * n1 * (n1-1)) ;S+= n1 * (n1-1);
                               // (+1,+1-1,-1,+1,0) [*] [[***]] || [***] [[*]]
                            if(n1 > 0 && n2 > 0) SetIndex5(nr,n0+1,n1,n2-1,n3+1, ant * 2 * n1 * 2 * n2) ;S+= 2 * n1 * 2 * n2;
                                    // (+1,+1-1,0,1-1,0) [*] [[[***]]] || [***] [[[*]]]
                            if(n1 > 0 && n3 > 0) SetIndex5(nr,n0+1,n1,n2,n3, ant * 2 * n1 * 2 * n3) ;S+= 2 * n1 * 2 * n3;
                                    // (+1,-1,0,+2,-1) [*] [[[[***]]]] || [***] [[[[*]]]]
                            if(n1 > 0 && n4 > 0) SetIndex5(nr,n0+1,n1-1,n2,n3+2, ant * 2 * n1 * 2 * n4) ;S+= 2 * n1 * 2 * n4;
                                    // (+1,-1,+1,+1,-1) [*] [[[[***]]]] || [***] [[[[*]]]]
                            if(n1 > 0 && n4 > 0) SetIndex5(nr,n0+1,n1-1,n2+1,n3+1, ant * 2 * n1 * n4) ;S+= 2 * n1 * n4;
                                    // (0,+2,-2,+1,0) [[*]] [[***]] || [[***]] [[*]] ; double insert *|* *|*
                            if(n2 > 1 ) SetIndex5(nr,n0,n1+2,n2-2,n3+1, ant * 4 * n2 * (n2-1) );S+= 4 * n2 * (n2-1);
                                    // (+1,0,-1,+1,0) [[*]][[***]] || [[***]][[*]]
                            if(n2 > 0 ) SetIndex5(nr,n0+1,n1,n2-1,n3+1, ant * 2 * n2 ) ;S+= 2 * n2;
                                    //!! (0,+2,-1,+1-1,0) [[*]] [[[***]]] || [[***]] [[[*]]]
                            if(n2 > 0 && n3 > 0) SetIndex5(nr,n0,n1+2,n2-1,n3, ant *2 * 2 * n2 * 2 * n3) ;S+= 2 * 2 * n2 * 2 * n3;
                                    // (0,+1,-1,+2,-1) [[*]] [[[[***]]]] || [[***]] [[[[*]]]]
                            if(n2 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2-1,n3+2, ant * 2 * 2 * n2 * 2 * n4) ;S+= 2 * 2 * n2 * 2 * n4;
                                    // (0,+1,1-1,+1,-1) [[*]] [[[[***]]]] || [[***]] [[[[*]]]]
                            if(n2 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2,n3+1, ant * 2 * 2 * n2  * n4 );S+= 2 * 2 * n2  * n4;
                                    // (0,+2,0,+1-2,0) [[[*]]] [[[***]]] || [[[***]]] [[[*]]] ; double insert *|** *|**
                            if(n3 > 1 ) SetIndex5(nr,n0,n1+2,n2,n3-1, ant * 4 * n3 * (n3-1) );S+= 4 * n3 * (n3-1);
                                    // (+1,0,0,+1-1,0) [[[*]]][[[***]]] || [[[***]]][[[*]]]
                            if(n3 > 0 ) SetIndex5(nr,n0+1,n1,n2,n3, ant * 2 * n3 ) ;S+= 2 * n3;
                                    // (0,+1,0,+2-1,-1) [[[*]]] [[[[***]]]] || [[[***]]] [[[[*]]]]
                            if(n3 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2,n3+1, ant * 2 * 2 * n3 * 2 * n4 );S+= 2 * 2 * n3 * 2 * n4;
                                    // (0,+1,+1,+1-1,-1) [[[*]]] [[[[***]]]] || [[[***]]] [[[[*]]]]
                            if(n3 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2+1,n3, ant * 2 * 2 * n3  * n4) ;S+= 2 * 2 * n3  * n4;
                                // (0,0,0,+3,-2) [[[[*]]]] |[[[***]]]] || [[[[***]]]] [[[[*]]]] ; double insert *|*** *|***
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2,n3+3, ant * 4 * n4 * (n4-1)) ;S+= 4 * n4 * (n4-1);
                                // !!! (0,0,+1,+2,-2) [[[[*]]]] |[[[***]]]] || [[[[***]]]] [[[[*]]]] ; insert *|*** **|**
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2+1,n3+2, ant * 4 * n4 * (n4-1)) ;S+= 4 * n4 * (n4-1);
                                // (0,0,+2,+1,-2) [[[[*]]]] [[[[***]]]] || [[[[***]]]] [[[[*]]]] ; double insert **|** **|**
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2+2,n3+1, ant * n4 * (n4-1) );S+= n4 * (n4-1);
                                    // (0,+1,0,+1,-1) [[[[*]]]][[[[***]]]] || [[[[***]]]][[[[*]]]] ; double insert *|***-*|***
                            if(n4 > 0 ) SetIndex5(nr,n0,n1+1,n2,n3+1, ant * 2 * n4 ) ;S+= 2 * n4 ;
                                    // (0,+1,0,+1,-1) [[[[*]]]][[[[***]]]] || [[[[***]]]][[[[*]]]] ; double insert *|***-**|**
                            if(n4 > 0 ) SetIndex5(nr,n0,n1+1,n2,n3+1, ant * 4 * n4 ) ;S+= 4 * n4;
   //                         printf("S3=%lld ",S);

                            // ** **
  //                          printf("\n Insert ** ** ") ;

                            S=0 ;
                            // (0,0,+1,0,0) ]**[ ]**[
                            SetIndex5(nr,n0,n1,n2+1,n3, ant * Cnp(nb0,2)); S+= Cnp(nb0,2);
                            // (+1,-1,+1,0,0) ]**[ [**]
                            if(n1>0)SetIndex5(nr,n0+1,n1-1,n2+1,n3, ant * nb0 * n1 );S+= nb0 * n1;
                            // (0,+1,+1-1,0,0) ]**[ [[**]]
                            if(n2>0)SetIndex5(nr,n0,n1+1,n2,n3, ant * nb0 * 2 * n2 );S+= nb0 * 2 * n2;
                            // (0,+1,+1,-1,0) ]**[ [[[**]]]
                            if(n3>0)SetIndex5(nr,n0,n1+1,n2+1,n3-1, ant * nb0 * 2 * n3 );S+= nb0 * 2 * n3;
                            // (0,0,+1,+1,-1) ]**[ [[[[**]]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3+1, ant * nb0 * 2 * n4 );S+= nb0 * 2 * n4;
                            // (0,+0,0,+2,-1) ]**[ [[[[**]]]]
                            if(n4>0)SetIndex5(nr,n0,n1,n2+2,n3, ant * nb0  * n4) ;S+= nb0  * n4;
                            // (+2,-2,+1,0,0) [**] [**] || [**] [**]
                            if(n1 > 1 ) SetIndex5(nr,n0+2,n1-2,n2+1,n3, ant * Cnp(n1,2)) ;S+= Cnp(n1,2);
                            // (+1,+1-1,+1-1,0,0) [**] [[**]]
                            if(n1 > 0 && n2 > 0) SetIndex5(nr,n0+1,n1,n2,n3, ant * n1 * 2 * n2 );S+= n1 * 2 * n2;
                            // (+1,+1-1,+1,-1,0) [**] [[[**]]]
                            if(n1 > 0 && n3 > 0) SetIndex5(nr,n0+1,n1,n2+1,n3-1, ant * n1 * 2 * n3) ;S+= n1 * 2 * n3;
                            // (+1,-1,+1,+1,-1) [**] [[[[**]]]]
                            if(n1 > 0 && n4 > 0) SetIndex5(nr,n0+1,n1-1,n2+1,n3+1, ant * n1 * 2 * n4) ;S+= n1 * 2 * n4;
                            // (+1,-1,+2,0,-1) [**] [[[[**]]]]
                            if(n1 > 0 && n4 > 0) SetIndex5(nr,n0+1,n1-1,n2+2,n3, ant  * n1 * n4) ;S+= n1 * n4;
                            // (0,+2,-2+1,0,0) [[**]] [[**]] ; double insert *|* *|*
                            if(n2 > 1 ) SetIndex5(nr,n0,n1+2,n2-1,n3, ant * 4 * Cnp(n2,2)) ;S+= 4 * Cnp(n2,2);
                            // (+1,0,1-1,0,0) [[**]][[**]]
                            if(n2 > 0 ) SetIndex5(nr,n0+1,n1,n2,n3, ant * n2 ) ;S+= n2 ;
                            // (0,+2,+1-1,-1,0) [[**]] [[[**]]]
                            if(n2 > 0 && n3 > 0) SetIndex5(nr,n0,n1+2,n2,n3-1, ant * 2 * n2 * 2 * n3 );S+= 2 * n2 * 2 * n3;
                            // (0,+1,+1-1,+1,-1) [[**]] [[[[**]]]]
                            if(n2 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2,n3+1, ant * 2* n2 * 2 * n4 );S+= 2 * n2 * 2 * n4;
                            // (0,+1,2-1,0,-1) [[**]] [[[[**]]]]
                            if(n2 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2+1,n3, ant  *2 * n2  * n4 );S+= n2 * 2 * n4;
                            // (0,+2,+1,-2,0) [[[**]]] [[[**]]]  ; double insert *|** *|**
                            if(n3 > 1 ) SetIndex5(nr,n0,n1+2,n2+1,n3-2, ant * 4 * Cnp(n3,2) ) ;S+= 4 * Cnp(n3,2);
                            // (+1,0,+1,-1,0) [[[**]]][[[**]]]
                            if(n3 > 0 ) SetIndex5(nr,n0+1,n1,n2+1,n3-1, ant  * n3 ) ;S+= n3 ;
                            // (0,+1,+1,+1-1,-1) [[[**]]] [[[[**]]]]
                            if(n3 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2+1,n3, ant * 2 * n3 * 2 * n4 );S+= 2 * n3 * 2 * n4 ;
                            // (0,+1,+2,-1,-1) [[[**]]] [[[[**]]]]
                            if(n3 > 0 && n4 > 0) SetIndex5(nr,n0,n1+1,n2+2,n3-1, ant  * 2 * n3  * n4 );S+= 2 * n3  * n4;
                            // (0,0,+1,+2,-2) [[[[**]]]] |[[[**]]]] ; double insert *|*** *|***
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2+1,n3+2, ant * 4 * Cnp(n4,2) );S+= 4 * Cnp(n4,2);
                            // (0,0,+2,+1,-2) [[[[**]]]] [[[[**]]]]  ; insert *|*** **|**
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2+2,n3+1, ant * 4 *  Cnp(n4,2));S+= 4 *  Cnp(n4,2);
                            // (0,0,+3,0,-2) [[[[**]]]] [[[[**]]]] || [[[[**]]]] [[[[**]]]] ; double insert **|** **|**
                            if(n4 > 1 ) SetIndex5(nr,n0,n1,n2+3,n3, ant * Cnp(n4,2)) ;S+=  Cnp(n4,2);
                            // (0,+1,+1,0,-1) [[[[**]]]][[[[**]]]] ; double insert *|***-*|***
                            if(n4 > 0 ) SetIndex5(nr,n0,n1+1,n2+1,n3, ant  * n4 ) ;S+= n4;
                            // (0,+1,+1,0,-1) [[[[**]]]][[[[**]]]] || [[[[**]]]][[[[**]]]] ; double insert *|***-**|**
                            if(n4 > 0 ) SetIndex5(nr,n0,n1+1,n2+1,n3, ant * 2 * n4 ) ;S+= 2 * n4 ;
   //                         printf("S2=%lld ",S);
                            // ** * *
    //                        printf("\n Insert ** * * ") ;

                            int pow3[4] = { 1,3,9,27} ;

                            S=0 ;
                            int i0,i1,i21,i21d,i22,i31,i32,i31d,i41_2,i41_3,i42,i42d     ,i41d,d0,d1,d2,d3,d4,d;
                            for(d=0;d<=4;d++) {
                                d0=(d==0) ?1 : 0 ; d1=(d==1) ?1 : 0 ;d2=(d==2) ?1 : 0 ;d3=(d==3) ?1 : 0 ;
                                d4=(d==4) ?1 : 0 ;
                            for(i0=2;i0>=0;i0--){
                                if(i0+d0 > nb0) continue ;
                                for(i1=2-i0;i1>=0;i1--){
                                    if(i1+d1 > n1) continue ;
                                    for(i22=(2-i0-i1)/2;i22>=0;i22--){
                                        for(i21d=0;i21d<=d2;i21d++){
                                            for(i21=2-i0-i1-2*i22-i21d;i21>=0;i21--){
                                                if(i21+i22+d2> n2) continue ;
                                                for(i32=(2-i0-i1-i21-i21d-2*i22)/2;i32>=0;i32--){
                                                    for(i31d=0;i31d<=d3;i31d++) {
                                                        for(i31=2-i0-i1-i21-i21d-2*i22-2*i32-i31d;i31>=0;i31--){
                                                            if(i31+i32+d3> n3) continue ;
                                                            for(i41d=0;i41d<=d4;i41d++) {
                                                                  for(i42d=0;i42d<=d4-i41d;i42d++) { for(i42=(2-i0-i1-i21-i21d-2*i22-i31-i31d-2*i32-i41d-2*i42d)/2;i42>=0;i42--){
                                                                for(i41_2=2-i0-i1-i21-i21d-2*i22-i31-i31d-2*i32-i41d-2*i42d-2*i42;i41_2>=0;i41_2--) {
                                                                                i41_3 = 2-i0-i1-i21-i21d-2*i22-i31-i31d-2*i32
                                                                                            -i41d-2*i42d-2*i42-i41_2 ;
                                                                    
                                                                                int64_t delta = 1 ;
                                                                                if((i41_2+i41_3+i42+d4) >n4) continue ;
                                                                                if(d0) delta *= nb0 ;
                                                                                if(i0) delta *= Cnp(nb0-d0,i0) ;
                                                                                if(d1) delta *= n1 ;
                                                                                if(i1) delta *= Cnp(n1-d1,i1) ;
                                                                                if(d2) delta *= n2*2 ;
                                                                                if(i22) delta *= Cnp(n2-d2,i22) ;
                                                                                if(i21){
                                                                                    /*if(i21d) {
                                                                                        delta *= Cnp(n2-i22-1,i21)*(1<<(i21)) ;
                                                                                    }else if(d2) {
                                                                                        delta *= Cnp(n2-i22-d2,i21)*(1<<i21) ;
                                                                                    } else {
                                                                                        delta *= Cnp(n2-i22,i21)*(1<<i21)  ;
                                                                                    }*/
                                                                                    delta *= Cnp(n2-i22-d2,i21)*(1<<i21) ;
                                                                                }
                                                                                if(d3) delta *= n3*2 ;
                                                                                if(i32) delta *= Cnp(n3-d3,i32) ;
                                                                                if(i31) {
                                                                                    delta *= Cnp(n3-i32-d3,i31)*(1<<i31) ;
                                                                                  }
                                                                                if(d4) {
                                                                                    //delta *= n4*3 ;
                                                                                    if(n4==0) continue ;
                                                                                    delta *= n4 ;
                                                                                    if(i42+i42d) {
                                                                                            SetIndex5( nr
                                                                                                      ,n0+i1+i22+i32+i42d
                                                                                                      ,n1-i1+i21+i31+i42+1
                                                                                                      ,n2-i22-i21+d4-i42d
                                                                                                      ,n3-i32-i31, ant * delta * Cnp(n4-d4,i42)*pow3[i42]);
                                                                                            SetIndex5( nr
                                                                                                      ,n0+i1+i22+i32+i42d
                                                                                                      ,n1-i1+i21+i31+i42+1
                                                                                                      ,n2-i22-i21
                                                                                                      ,n3-i32-i31+d4-i42d, ant * delta* 2 * Cnp(n4-d4,i42)*pow3[i42]);
                                                                                            delta *= 3 *  (Cnp(n4-d4,i42)*pow3[i42]) ;
 
                           
                                                                                    } else if(i41_2+i41_3+i41d){
                                                                                        int i41 = i41_2+i41_3 ;
                                                                                        SetIndex5( nr
                                                                                                ,n0+i1+i22+i32
                                                                                                ,n1-i1+i21+i31+1+i41d
                                                                                                ,n2-i22-i21+i41_2+d4-i41d
                                                                                                ,n3-i32-i31+i41_3,
                                                                                            ant * delta * Cnp(n4-d4,i41) * (1<<i41_3) * Cnp(i41,i41_3)* (1+i41d) );
                                                                                        SetIndex5( nr
                                                                                                ,n0+i1+i22+i32
                                                                                                ,n1-i1+i21+i31+1+i41d
                                                                                                ,n2-i22-i21+i41_2
                                                                                                ,n3-i32-i31+i41_3+d4-i41d,
                                                                                            ant * delta  * Cnp(n4-d4,i41) * (1<<i41_3) * Cnp(i41,i41_3)* 2 * (1+i41d));
                                                                                        delta *= Cnp(n4-d4,i41)* (1<<i41_3) * Cnp(i41,i41_3) * 3 *(1+i41d) ;
                                                                                    } else {
                                                                                        SetIndex5( nr
                                                                                                    ,n0+i1+i22+i32
                                                                                                    ,n1-i1+i21+i31+1
                                                                                                    ,n2-i22-i21+d4
                                                                                                    ,n3-i32-i31, ant * delta );
                                                                                        SetIndex5( nr
                                                                                                    ,n0+i1+i22+i32
                                                                                                    ,n1-i1+i21+i31+1
                                                                                                    ,n2-i22-i21
                                                                                                    ,n3-i32-i31+d4, ant * delta * 2 );
                                                                                        delta *= 3 ;
                                                                                    }
                                                                                } else {
                                                                                    if(i42) delta *= Cnp(n4,i42)*pow3[i42]  ;
                                                                                    int i41= i41_2 + i41_3 ;
                                                                                    if(i41) {
//                                                                                     int j3 ;
  //                                                                                      for(j3=0;j3<=i41;j3++)
                                                                                        {
                                                                                            SetIndex5( nr
                                                                                                        ,n0+i1+d1+i22+i21d+i32+i31d
                                                                                                        ,n1-i1-d1+i21+d2-i21d+i31+d3-i31d+i42+1
                                                                                                        ,n2-i22-i21-d2+i41_2
                                                                                                        ,n3-i32-i31-d3+i41_3, ant * delta * Cnp(n4-i42,i41) * (1<<i41_3) * Cnp(i41,i41_3));
                                                                                        }

                                                                                        delta *= Cnp(n4-i42,i41) * (1<<i41_3) * Cnp(i41,i41_3) ;
                                                                                    } else {
                                                                                        SetIndex5( nr
                                                                                                    ,n0+i1+d1+i22+i21d+i32+i31d
                                                                                                    ,n1-i1-d1+i21+d2-i21d+i31+d3-i31d+i42+1
                                                                                                    ,n2-i22-i21-d2
                                                                                                    ,n3-i32-i31-d3, ant * delta) ;
                                                                                    }
                                                                                }
                                                                                S += delta ;
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
    //                        printf("S1=%lld ",S) ;

                            
       //                     printf("\n Insert * * * * ") ;

                            // * * * *
                            int i43 ;
                            S=0 ;
   //                          int i0,i1,i21,i22,i31,i32,i43,i42;
                            for(i0=4;i0>=0;i0--){
                                for(i1=4-i0;i1>=0;i1--){
                                    for(i21=4-i0-i1;i21>=0;i21--){
                                        for(i22=(4-i0-i1-i21)/2;i22>=0;i22--){
                                            for(i31=4-i0-i1-i21-2*i22;i31>=0;i31--){
                                                for(i32=(4-i0-i1-i21-2*i22-i31)/2;i32>=0;i32--){
                                                    for(i43=(4-i0-i1-i21-2*i22-i31-2*i32)/3;i43>=0;i43--){
                                                        for(i42=(4-i0-i1-i21-2*i22-i31-2*i32-3*i43)/2;i42>=0;i42--){
                                                            int i41 = 4-i0-i1-i21-2*i22-i31-2*i32-3*i43 -2*i42 ;
  
                                                            int64_t delta = 1 ;
                                                            if(i0>nb0 || i1>n1 || (i21+i22) >n2
                                                               || (i31+i32) >n3 || (i41+i42+i43) >n4) continue ;
                                                            if(i0) delta *= Cnp(nb0,i0) ;
                                                            if(i1) delta *= Cnp(n1,i1) ;
                                                            if(i22) delta *= Cnp(n2,i22) ;
                                                            if(i21) delta *= Cnp(n2-i22,i21)*(1<<i21)  ;
                                                            if(i32) delta *= Cnp(n3,i32) ;
                                                            if(i31) delta *= Cnp(n3-i32,i31)*(1<<i31)  ;
                                                            if(i43) delta *= Cnp(n4,i43) ;
                                                            if(i42) delta *= Cnp(n4-i43,i42)*pow3[i42]  ;
                                                            if(i41) {
                                                                int j3 ;
                                                                for(j3=0;j3<=i41;j3++) {
                                                                    SetIndex5( nr
                                                                                ,n0+i1+i22+i32+i43+1
                                                                                ,n1-i1+i21+i31+i42
                                                                                ,n2-i22-i21+i41-j3
                                                                                ,n3-i32-i31+j3, ant * delta * Cnp(n4-i43-i42,i41) * (1<<j3) * Cnp(i41,j3));
                                                                }
                                                                delta *= Cnp(n4-i43-i42,i41)*pow3[i41]  ;
                                                                
                                                            } else {
                                                                SetIndex5( nr
                                                                        ,n0+i1+i22+i32+i43+1
                                                                        ,n1-i1+i21+i31+i42
                                                                        ,n2-i22-i21
                                                                        ,n3-i32-i31, ant * delta );
                                                            }
                                                            S += delta ;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                   }
                                }
                            }
 //                           printf("S0=%lld ",S) ;
                        }
                    }
                }
            }
#if defined(FLOAT687)
            PROB687 g = 1 ;
#else
            PROB687 g =cur5[0];
            for(i=1;i<=Index5(nr,0,0,0,0);i++) { g = PGCD64(g,cur5[i]) ;} ;
#endif
            PROB687 den = 0;
            for(i=0;i<=Index5(nr,0,0,0,0);i++) { /*cur5[i] /= g;*/ den += cur5[i]/g ; }
            PROB687 prob[20];
            for(i=0;i<=nr;i++)prob[i]=0 ;
            for(n0=0;n0<=nr;n0++) {
                for(n1=0;n1<=nr-n0;n1++) {
                    for(n2=0;n2<=nr-n0-n1;n2++) {
                        for(n3=0;n3<=nr-n0-n1-n2;n3++) {
                            //                                 printf("(%d,%d,%d,%d,%d)=%lld ",n0,n1,n2,n3,nr-n0-n1-n2-n3,cur5[Index5(nr,n0,n1,n2,n3)]);
                            cur5[Index5(nr,n0,n1,n2,n3)] /= g ;
                            prob[n0]+= cur5[Index5(nr,n0,n1,n2,n3)] ;
                        }
                    }
                }
            }
            
            printf("\n4 color, %d rank  => nb index=%d",nr,Index5(nr,0,0,0,0));
            printf("\n");
#if defined(FLOAT687)
            for(i=0;i<=nr;i++) { printf("%.11f ",prob[i]/(double)den) ; } ;
#else
            for(i=0;i<=nr;i++) { printf("%lld ",prob[i]) ; } ;
#endif
            if(nr==MAX_NR) {
                printf("\n answer=%.11f \n", (double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13])/den) ;
                snprintf(pbR->strRes, sizeof(pbR->strRes),"%.10f",(double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13])/den);
                pbR->nbClock = clock() - pbR->nbClock ;
                return 1 ;
            }
#if defined(FLOAT687)
            for(i=0;i<=Index5(nr,0,0,0,0);i++) ant5[i] = cur5[i]/den ;
#else
            for(i=0;i<=Index5(nr,0,0,0,0);i++) ant5[i] = cur5[i] ;
#endif
        }
        printf("\n*********\n") ;

    }
    uint64_t ant3[200] ;
    uint64_t cur3[200] ;
  
    ant3[Index3(1,0,0)] = 1 ;
    ant3[Index3(1,1,0)] = 0 ;
    ant3[Index3(1,0,1)] = 0 ;
    int n0,n1;
    n1 = 0 ;
    
    for(nr=2;nr<=4;nr++) {
        int i ;
        for(i=0;i<=Index3(nr,0,0);i++) cur3[i] = 0 ;
        // pour rajouter un rang : cela peut être : *,*,* || **,* || ***
        for(n0=0;n0<nr;n0++) {
            for(n1=0;n1<nr-n0;n1++) {
                int n2=nr-1-n0-n1  ;
//                int nb0 = 3*n0+2*n1+n2 ;//
               int nb0 = 3*n0+2*n1+n2 ;//
                 uint64_t ant = ant3[Index3(nr-1,n0,n1)] ;
                if(ant <= 0) continue ;
                // ***
                cur3[Index3(nr,n0,n1)] += ant * nb0 ; // (0,0,+1) ]***[
                if(n1>0)cur3[Index3(nr,n0+1,n1-1)] += ant * n1 ; // (+1,-1,+1) [***]
                if(n2>0)cur3[Index3(nr,n0,n1+1)] += ant * 2 * n2 ; // (0,+1,+1-1) [[***]]
                // ** *
                cur3[Index3(nr,n0,n1+1)] += ant * nb0 * (nb0-1); // (0,+1,0) ]*[ ]**[
                cur3[Index3(nr,n0+1,n1)] += ant * 2 * nb0 * n1 ; // (+1,+-0,0) ]*[ [**] || ]**[ [*]
                cur3[Index3(nr,n0,n1+2)] += ant * 2 * nb0 * 2 * n2 ; // (0,+2,-1) ]*[ [[**]] || ]**[ [[*]]
                if(n1 > 0 && n2 > 0) cur3[Index3(nr,n0+1,n1+1)] += ant * 2 * n1 * 2 * n2 ; // (+1,+2-1,-1) [*] [[**]] || [**] [[*]]
                if(n1 > 1 ) cur3[Index3(nr,n0+2,n1-1)] += ant * n1 * (n1-1) ; // (+2,+1-2,0) [*] [**]
                if(n2 > 0 ) cur3[Index3(nr,n0+1,n1+1)] += ant * 2 * n2; // (+1,+1,-1) [[*]][[**]]
                if(n2 > 1 ) cur3[Index3(nr,n0,n1+3)] += ant * 2* n2 * 2 * (n2-1) ; // (0,+3,-2) [[*]] [[**]]
                // * * *
                cur3[Index3(nr,n0+1,n1)] += ant * Cnp(nb0,3) ; // +1,0,0  ]*[,]*[,]*[
                if(n1 > 0) cur3[Index3(nr,n0+2,n1-1)] += ant * Cnp(nb0,2) * n1 ; // +2,-1,0  ]*[,]*[,[*]
                if(n2 > 0) cur3[Index3(nr,n0+1,n1+1)] += ant * Cnp(nb0,2) * 2 * n2 ; // +1,+1,-1  ]*[,]*[,[[*]]
                if(n1 > 0 && n2 > 0) cur3[Index3(nr,n0+2,n1)] += ant * nb0 * n1 * 2 * n2 ; // +2,+-0,-1  ]*[,[*],[[*]]
                if(n1 > 1)cur3[Index3(nr,n0+3,n1-2)] += ant * nb0 * Cnp(n1,2)  ; // +3,-2,0  ]*[,[*],[*]
                if(n2 > 0 ) cur3[Index3(nr,n0+2,n1)] += ant  * nb0 * n2 ; // (+2,0,-1) ]*[ [[*]][[*]]
                if(n2 > 1 ) cur3[Index3(nr,n0+1,n1+2)] += ant * nb0 * 4* Cnp(n2,2) ; // (+1,+2,-2) ]*[ [[*]] [[*]]
                if(n1 > 2)cur3[Index3(nr,n0+4,n1-3)] += ant *  Cnp(n1,3)  ; // +4,-3,0  [*],[*],[*]
                if(n1 > 1 && n2 > 0)cur3[Index3(nr,n0+3,n1-1)] += ant *  Cnp(n1,2) * 2 * n2 ; // +3,-1,-1  [*],[*],[[*]]
                if(n1 > 0 && n2 > 0)cur3[Index3(nr,n0+3,n1-1)] += ant *  n1 * n2 ; // +3,-1,-1  [*],[[*]][[*]]
                if(n1 > 0 && n2 > 1)cur3[Index3(nr,n0+2,n1+1)] += ant *  n1 * 4 * Cnp(n2,2) ; // +2,+1,-2  [*],[[*]] [[*]]
                if(n2 > 1) cur3[Index3(nr,n0+2,n1+1)] += ant *  n2 * 2 * (n2-1) ; // +2,+1,-2  [[*]][[*]] [[*]]
                if(n2 > 2) cur3[Index3(nr,n0+1,n1+3)] += ant *  8 * Cnp(n2,3) ; // +1,+3,-3  [[*]] [[*]] [[*]]
            }
        }
        int64_t g =cur3[0];
        for(i=1;i<=Index3(nr,0,0);i++) { g = PGCD64(g,cur3[i]) ;} ;
        int64_t den = 0;
        for(i=0;i<=Index3(nr,0,0);i++) { cur3[i] /= g; den += cur3[i] ; }
        int64_t prob[20];
        for(i=0;i<=nr;i++)prob[i]=0 ;
        for(n0=0;n0<=nr;n0++) {
            for(n1=0;n1<=nr-n0;n1++) {
                //               printf("(%d,%d,%d)=%lld ",n0,n1,nr-n0-n1,cur3[INDT2(n0,n1)]);
                prob[nr-n0]+= cur3[Index3(nr,n0,n1)] ;
            }
        }
        
        printf("\n3 color, %d rank  => %lld=%lld/%lld=",nr,den*g,den,g) ;
        //      g=prob[0] ;
        //       for(i=1;i<=nr;i++) g = PGCD64(g,prob[i]) ;
        for(i=0;i<=nr;i++) { printf("%lld,",prob[i]) ; } ;
        for(i=0;i<=Index3(nr,0,0);i++) ant3[i] = cur3[i] ;
    }
    return 0 ;
}
