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


#include "PB_other.h"

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
    u_int64_t nbVertice ;
    u_int64_t nbCube ;
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
    CC.nbVertice =(((u_int64_t)(N))*(N+1)*(3 * N*N*N*N*N + 39 * N*N*N*N + 213 * N*N*N + 627 * N*N + 640 * N + 158 ))/ 420 ;
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
                    Diff -= 2*n2+1 , n2++  ;
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
    sprintf(pbR->strRes,"%lld",CC.nbVertice);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB620_MAX   500 //
int Count620(int s, int p,int q) {
    int k ;
    int nbSol =0 ;
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
    sprintf (pbR->strRes,"%d",nbSol);
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
    sprintf (pbR->strRes,"%d",nbSol);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



typedef struct P_EXP {
    u_int32_t   p;
    u_int32_t   exp ;
    
} P_EXP ;
#define PB622_NBP   11 // number of prime factors for 2**60-1
int PB622(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int i, jp,k;
    u_int64_t S=0;
    P_EXP fact2exp60_m1[PB622_NBP] = {{3, 2}, {5, 2}, {7, 1}, {11, 1}, {13, 1}, {31, 1}, {41, 1}, {61,1}, {151, 1}, {331, 1}, {1321, 1}};
    P_EXP fact2div60byP[3][PB622_NBP] ;
    u_int32_t quotPrime60[3] = { (1<<(60/2))-1 , (1<<(60/3)) -1 , (1<<(60/5)) -1 } ;
    {
        for(i=0;i<PB622_NBP;i++) {
            u_int32_t p = fact2exp60_m1[i].p ;
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
                    u_int64_t N ;
                    for(k=0,N=1;k<PB622_NBP;k++) {
                        int ie ;
                        for(ie=exp[k];ie>0;ie--) N *= fact2exp60_m1[k].p ;
                    }
                    S += N+1 ;
                }
            } else {  ip-- ; }
        }
    }
    sprintf (pbR->strRes,"%lld",S);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB625_MAX   100000000000LL
// #define PB625_MAX   10000000LL //
// #define PB625_MAX   10LL //
// #define PB625_MAX   1000000L
#define PB625_MOD 998244353
// #define PB625_MOD 100000




int PB625(PB_RESULT *pbR) {
    int plim = Sqrt64(PB625_MAX) ;
    pbR->nbClock = clock() ;
    int i ,p ,j , N = PB625_MAX ;
#if 1
    if(PB625_MAX > 1000000000L) {
        printf("Too memory for %lld\n",PB625_MAX) ;
        return 0 ;
    }
    u_int64_t *nbDprime=calloc(PB625_MAX+1,sizeof(nbDprime[0])) ;
    u_int64_t nbP = 1 ;
    for(p=2;p<=PB625_MAX;p++)  {
        if(nbDprime[p]) {
            nbP += nbDprime[p] ;
            nbP =  nbP % PB625_MOD ;
            nbDprime[p] = nbP ;
//            printf("(%d,%d)",p,nbDprime[p]) ;
            continue ;
        } else {
            int pm = p-1 ;
            nbP += pm ;
            nbP =  nbP % PB625_MOD ;
            nbDprime[p] = nbP ;
            int k ;
//            printf("(%d,%d)",p,nbDprime[p]) ;
            for(j=2*p,k=2;j<=PB625_MAX;j+=p,k++) {
                if(nbDprime[j]) {
                    nbDprime[j] = (nbDprime[j]/p)*(p-1)  ;
                } else {
                    nbDprime[j] = k*(p-1)   ;
                }
                
            }
        }
        
    }
    
    u_int64_t isum =  (N+(N/2)+1)  ;
    u_int64_t idif =  (N-(N/2))  ;
    if(isum & 1) idif /=2 ;
    else isum /=2 ;
    idif = idif % PB625_MOD ;
    isum = isum % PB625_MOD ;

    u_int64_t S = (isum*idif) % PB625_MOD ;
    printf("SN=%lld ",S );
    
    //    printf("SN=%lld ",S );
    for(i=1;2*i<=N;i++) {
//    for(i=N/2;i>=1;i--) {
        if(i==plim) printf("(i=%d S=%lld)",i,S);
        S += ((u_int64_t)(i % PB625_MOD ) *  nbDprime[N/i] ) % PB625_MOD;
        S = S % PB625_MOD ;
//        printf("(%d,%lld,%lld)",i,nbDprime[N/i],S);
    }
    printf("\nMAX=%lld[%d] S=%lld Expected 317257140\n",PB625_MAX,PB625_MOD,S );
    free(nbDprime) ;
#else
    u_int64_t S = 0 ;
    for(i=1;i<=N;i++) {
        for(j=1;j<=i;j++) {
            S += PGCD(i,j);
        }
    }
    printf("\nMAX=%lld[%d] S=%lld Expected 317257140\n",PB625_MAX,PB625_MOD,S );
#endif
     pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
typedef struct P625  {
    u_int64_t p ;
    u_int64_t q ;
    u_int64_t ioffset ;
} P625 ;

int PB625a(PB_RESULT *pbR) {
    int pLim = Sqrt64(PB625_MAX)+1 ;
    int sizeD = pLim;
    pbR->nbClock = clock() ;
    int32_t i,j ;
    u_int64_t p , N = PB625_MAX ;
    u_int64_t *nbDprime=calloc(sizeD,sizeof(nbDprime[0])) ;
//    u_int32_t *lastPrime=calloc(sizeD,sizeof(lastPrime[0])) ;
    int64_t nbPrimeAlloc = (PB625_MAX / log((double)PB625_MAX) * (1.0 +4.0 / log((double)PB625_MAX)))+10 ;
    P625 * tbPrime = malloc(nbPrimeAlloc * sizeof(tbPrime[0])) ;
    int nbPrime = 0 ;
    //    u_int32_t *nbDprimev=calloc(PB625_MAX+1,sizeof(nbDprime[0])) ;
    
    u_int64_t nbP = 1 ;
    
    for(p=2;p<sizeD;p++)  {
        if(nbDprime[p-1]) {
            nbP += nbDprime[p-1] ;
            nbP = nbP % PB625_MOD ;
            nbDprime[p-1] = nbP ;
//            printf("[%lld,%lld]",p,nbDprime[p-1]) ;
            //            printf("(%d,%d)",p,nbDprime[p]-nbDprime[p-1]) ;
            continue ;
        } else {
            tbPrime[nbPrime].p = p ;
            u_int64_t pm = p-1 ;
            nbP += pm ;
            nbDprime[p-1] = nbP ;
            u_int64_t q ;
 //           printf("[%lld,%lld]",p,nbDprime[p-1]) ;
            for(j=2*p,q=2;j<sizeD;j+=p,q++) {
                if(nbDprime[j-1]) {
                    nbDprime[j-1] = (nbDprime[j-1]/p)*(p-1)  ;
                } else {
                    nbDprime[j-1] = q*(p-1)   ;
                }
                
            }
            tbPrime[nbPrime].q = q ;
            tbPrime[nbPrime++].ioffset = j - sizeD ;
        }
    }
    
    
    u_int64_t isum =  (N+(N/2)+1)  ;
    u_int64_t idif =  (N-(N/2))  ;
    if(isum & 1) idif /=2 ;
    else isum /=2 ;
    idif = idif % PB625_MOD ;
    isum = isum % PB625_MOD ;

    u_int64_t S = (isum*idif) % PB625_MOD ;
     printf("SN=%lld ",S );

    
    u_int64_t pl ;
   u_int64_t S1 = 0 ;
    {
        u_int64_t i ;
    for(i=sizeD;2*i<=N;i++) {
        //        if(i<4) printf("(%d->(%d,%d,%d ))",i,N/i,nbDprime[N/i-1],(N/i*(N/i+1))/2 - nbDprime[N/i-1]);
        S1 += ((u_int64_t)i * nbDprime[N/i-1]) % PB625_MOD ;
//        printf("[%lld,%lld=>%lld]",i,nbDprime[N/i-1],S1) ;
        S1 = S1 % PB625_MOD ;
    }
    }
    u_int64_t pgNext = sizeD-1 ;
//    int pgNext = N / sizeD + 2 ;
    u_int64_t invPgNext = N /pgNext ;
    while(invPgNext < sizeD) {
        S1 += (u_int64_t)pgNext * ( nbDprime[invPgNext-1] % PB625_MOD);
//        printf("{%lld,%lld,%lld}",pgNext,nbDprime[invPgNext-1],S1) ;
        S1 = S1 % PB625_MOD ;
        pgNext-- ;
        invPgNext = N /pgNext ;
    }
    for(pl=sizeD;pl<=PB625_MAX;pl+=sizeD ){
        memset(nbDprime,0,sizeD*sizeof(nbDprime[0])) ;
//        memset(lastPrime,0,sizeD*sizeof(lastPrime[0])) ;
        
        int32_t ip ;
        for(ip=0;ip<nbPrime;ip++) {
            if(tbPrime[ip].ioffset < sizeD) {
                p = tbPrime[ip].p ;
                u_int64_t q = tbPrime[ip].q ;
                u_int64_t j ;
                for(j=tbPrime[ip].ioffset;j<sizeD;j+=p,q++) {

                    if(nbDprime[j]) {
                        nbDprime[j] = ((nbDprime[j])/p)*(p-1)  ;
                    } else {
                        nbDprime[j] = q*(p-1)   ;
                    }
                }
                tbPrime[ip].q = q ;
                tbPrime[ip].ioffset = j - sizeD ;
            } else {
                tbPrime[ip].ioffset -= sizeD ;
            }
        }
        for(j=0;j<sizeD;j++) {
            if(nbDprime[j]) {
                nbP += nbDprime[j] ;
//                printf("(%lld,%lld)",j+pl,nbP);
            } else {
                p = pl+j ;
                nbP +=   p-1 ;
                tbPrime[nbPrime].p =p ;
//                printf("{%lld,%lld}",j+pl,nbP);

                int q ;
                //            printf("(%d,%d)",p,nbDprime[p]-nbDprime[p-1]) ;
                if(p < sizeD-j) {
                    u_int64_t ij ;
                    for(ij=j+p,q=2;ij<sizeD;ij+=p,q++) {
                        if(nbDprime[ij]) {
                            nbDprime[ij] = (nbDprime[ij]/p)*(p-1)  ;
                        } else {
                            nbDprime[ij] = q*(p-1)   ;
                        }
                    
                    }
                    tbPrime[nbPrime].q = q ;
                    tbPrime[nbPrime++].ioffset = ij - sizeD ;
                    if(nbPrime >= nbPrimeAlloc) {
                        printf("!!");
                    }
                } else {
                   tbPrime[nbPrime].q = 2 ;
                    tbPrime[nbPrime++].ioffset = j+p - sizeD ;
                    if(nbPrime >= nbPrimeAlloc) {
                        printf("!!");
                    }
                }
            }
            nbP = nbP % PB625_MOD ;
            if(pl+j == invPgNext){
                S += ( nbP * pgNext ) % PB625_MOD ;
                S = S % PB625_MOD ;
//               printf("(%d,%lld,%lld)",pgNext,nbP,S) ;
                pgNext-- ;
                if(pgNext== 0) break ;
                invPgNext = N / pgNext ;
            }
        }
        if(pgNext <= 0) break ;
    }
    printf("S=%lld",S) ;
    S =(S+S1) % PB625_MOD ;
    printf("\nMAX=%lld[%d] S=%lld Expected 317257140\n",PB625_MAX,PB625_MOD,S );
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

typedef struct N625 {
    u_int64_t   gcd ;
    u_int64_t   index ;
} N625 ;
int PB625b(PB_RESULT *pbR) {
//   int pLim = Sqrt64(PB625_MAX)+1 ;
    int pLim = (Sqrt64(PB625_MAX)+1) ;
    int sizeD = pLim ;
    pbR->nbClock = clock() ;
    int32_t j ;
    u_int64_t p , N = PB625_MAX ;
    u_int64_t *nbDprime=calloc(sizeD,sizeof(nbDprime[0])) ;
    u_int64_t *nbDprime2=malloc(sizeD*sizeof(nbDprime2[0])) ;
    N625 *nxt = malloc(sizeD*sizeof(nxt[0])) ;
    int64_t nbPrimeAlloc = sizeD ;
    P625 * tbPrime = malloc(nbPrimeAlloc * sizeof(tbPrime[0])) ;
    int nbPrime = 0 ;
    
    u_int64_t nbP = 1 ;
    nbDprime[1] = 1 ;
    {
        int j,p; // limited by sizeD, int3é is sufficient.
        for(p=2;p<sizeD;p++)  {
            if(nbDprime[p]) {
                nbP += nbDprime[p] ;
                nbP = nbP % PB625_MOD ;
                nbDprime[p] = nbP ;
    //           printf("[%lld,%lld]",p,nbDprime[p-1]) ;
                continue ;
            } else {
                tbPrime[nbPrime].p = p ;
                u_int64_t pm = p-1 ;
                nbP += pm ;
                nbP = nbP % PB625_MOD ;
                nbDprime[p] = nbP ;
    //            printf("[%lld,%lld]",p,nbDprime[p-1]) ;
                int q ;
                //            printf("(%d,%d)",p,nbDprime[p]-nbDprime[p-1]) ;
                for(j=2*p,q=2;j<sizeD;j+=p,q++) {
                    if(nbDprime[j]) {
                        nbDprime[j] = (nbDprime[j]/p)*(p-1)  ;
                    } else {
                        nbDprime[j] = q*(p-1)   ;
                    }
                    
                }
                tbPrime[nbPrime].q = q ;
                tbPrime[nbPrime++].ioffset = j - sizeD ;
            }
        }
    }
    u_int64_t isum =  (N+(N/2)+1) ;
    u_int64_t idif =  (N-(N/2))  ;
    if(isum & 1) idif /=2 ;
    else isum /=2 ;
    idif = idif % PB625_MOD ;
    isum = isum % PB625_MOD ;
    u_int64_t S = (isum*idif) % PB625_MOD ;
    printf("SN=%lld ",S );
    
    
    u_int64_t pl ;
    u_int64_t S1 = 0 , S2 = 0 ;
    u_int64_t i0,i1 ;

    {
        int32_t i ;
        for(i=2;i<=sizeD;i++) {
            i0 = N/i ;
            i1 = N/(i+1) ;
            u_int64_t isum = ((i0+i1+1) & 1 ) ? ((i0+i1+1) % PB625_MOD) : (((i0+i1+1)/2) % PB625_MOD) ;
            u_int64_t idif = ((i0+i1+1) & 1 ) ? (((i0-i1)/2) % PB625_MOD) : ((i0-i1) % PB625_MOD) ;
            u_int64_t nb= (isum*idif ) % PB625_MOD ;
            
            S1 += (nb * nbDprime[i]) % PB625_MOD ;
 //           printf("[%lld,%lld=>%lld]",i,nb,S1);
            S1 = S1 % PB625_MOD ;

        }
    }
    printf(" S1=%lld ",S1);
 //   u_int64_t pgNext = sizeD-1 ;
    u_int64_t pgNext = i0 ;

    u_int64_t invPgNext = N /pgNext ;
    while(invPgNext < sizeD) {
        S1 += ((u_int64_t)pgNext * nbDprime[invPgNext]) % PB625_MOD;
//        printf("{%lld,%lld}",pgNext,nbDprime[invPgNext-1]) ;
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
        memset(nbDprime2,0,sizeD*sizeof(nbDprime2[0])) ;
        int32_t ip ;
        for(ip=0;ip<nbPrime;ip++) {
            if(tbPrime[ip].ioffset < sizeD) {
                p = tbPrime[ip].p ;
                u_int64_t q = tbPrime[ip].q ;
                u_int64_t j ;
                for(j=tbPrime[ip].ioffset;j<sizeD;j+=p,q++) {
                    
                    if(nbDprime2[j]) {
                        nbDprime2[j] = ((nbDprime2[j])/p)*(p-1)  ;
                    } else {
                        nbDprime2[j] = q*(p-1)   ;
                    }
                }
                tbPrime[ip].q = q ;
                tbPrime[ip].ioffset = j - sizeD ;
            } else {
                tbPrime[ip].ioffset -= sizeD ;
            }
        }
        for(j=0;j<sizeD;j++) {
            if(nbDprime2[j]) {
                nbP += nbDprime2[j] % PB625_MOD ;
//                printf("(%lld,%lld)",pl+j,nbP);
            } else {
                p = pl+j ;
                nbP +=  ( p -1 ) % PB625_MOD;
//                printf("{%lld,%lld}",pl+j,nbP);
                int k ;
                for(k=indNxt;k<nbNxt;k++) {
                    int xp = nxt[k].index /p ;
                    S2 += ((nbDprime[xp]-1)*nxt[k].gcd ) % PB625_MOD;
                    S2 = S2 % PB625_MOD ;
//                    printf("[%lld,-%lld]",p,nbDprime[xp-1]-1);
                }
            }
            nbP = nbP % PB625_MOD ;
            if(pl+j == nxt[indNxt].index){
                S += ( nbP * nxt[indNxt].gcd ) % PB625_MOD ;
                S = S % PB625_MOD ;
//                printf("(%lld,%lld,%lld)",nxt[indNxt].gcd,nbP,S) ;
                indNxt++ ;
                if(indNxt >= nbNxt) break ;
            }
        }
        if(indNxt >= nbNxt)  break ;
    }
    printf("S=%lld",S) ;
    S =(S+S1-S2) % PB625_MOD ;
    printf("\nMAX=%lld[%d] S=%lld Expected 317257140\n",PB625_MAX,PB625_MOD,S );
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
    sprintf(pbR->strRes,"%u",tbPrime[PB1000_NUM-1]);
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

