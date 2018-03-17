//
//  PB_other.c
//  X_euler
//
//  Created by Jeannot on 13/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "PB_other.h"

#define PB579_SIZE  50
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

typedef struct cubeMin {
    int16_t Ay ;
    int16_t Az ;
    int16_t Bx ;
    int16_t Bz ;
    int16_t Cx ;
    int16_t Cy ;
} cubeMin ;


void MinimizeCube(cubeO *cube, cubeMin *cubeM) {

    cubeM->Ay = cubeM->Az = 0 ;
    cubeM->Bx = cubeM->Bz = 0 ;
    cubeM->Cx = cubeM->Cy = 0 ;
    if(cube->A.x< 0) {
        cubeM->Ay += cube->A.y ;
        cubeM->Az += cube->A.z ;
        cubeM->Bx -= cube->A.x ;
        cubeM->Cx -= cube->A.x ;

    }
    if(cube->B.x< 0) {
        cubeM->Ay += cube->B.y ;
        cubeM->Az += cube->B.z ;
        cubeM->Bx -= cube->B.x ;
        cubeM->Cx -= cube->B.x ;
    }
    if(cube->C.x< 0) {
        cubeM->Ay += cube->C.y ;
        cubeM->Az += cube->C.z ;
        cubeM->Bx -= cube->C.x ;
        cubeM->Cx -= cube->C.x ;
    }

    if(cube->A.y< 0) {
        cubeM->Bx += cube->A.x ;
        cubeM->Bz += cube->A.z ;
        cubeM->Ay -= cube->A.y ;
        cubeM->Cy -= cube->A.y ;
    }
    if(cube->B.y< 0) {
        cubeM->Bx += cube->B.x ;
        cubeM->Bz += cube->B.z ;
        cubeM->Ay -= cube->B.y ;
        cubeM->Cy -= cube->B.y ;
    }
    if(cube->C.y< 0) {
        cubeM->Bx += cube->C.x ;
        cubeM->Bz += cube->C.z ;
        cubeM->Ay -= cube->C.y ;
        cubeM->Cy -= cube->C.y ;
    }

    if(cube->A.z< 0) {
        cubeM->Cx += cube->A.x ;
        cubeM->Cy += cube->A.y ;
        cubeM->Az -= cube->A.z ;
        cubeM->Bz -= cube->A.z ;
    }
    if(cube->B.z< 0) {
        cubeM->Cx += cube->B.x ;
        cubeM->Cy += cube->B.y ;
        cubeM->Az -= cube->B.z ;
        cubeM->Bz -= cube->B.z ;
    }
    if(cube->C.z< 0) {
        cubeM->Cx += cube->C.x ;
        cubeM->Cy += cube->C.y ;
        cubeM->Az -= cube->C.z ;
        cubeM->Bz -= cube->C.z ;
    }
}

#define PB579_MAXT  20000

typedef struct CTX_Cube {
    int nbC ;
    int nbCand ;
    u_int64_t nbVertice ;
    u_int32_t nbCube ;
    u_int64_t nb8 ;
    u_int64_t nb12 ;
    u_int64_t nb24 ;
} CTX_Cube ;

static inline V3 PVect( V3 T1,V3 T2,int N) {
    V3 V ;
    V.x = (T1.y*T2.z - T1.z*T2.y)/N ;
    V.y = (T1.z*T2.x - T1.x*T2.z)/N ;
    V.z = (T1.x*T2.y - T1.y*T2.x)/N ;
    return V ;
}
static inline V3 PVectAbs( V3 T1,V3 T2,int N) {
    V3 V ;
    V.x = (T1.y*T2.z - T1.z*T2.y)/N ; if(V.x < 0) V.x = - V.x ;
    V.y = (T1.z*T2.x - T1.x*T2.z)/N ; if(V.y < 0) V.y = - V.y ;
    V.z = (T1.x*T2.y - T1.y*T2.x)/N ; if(V.z < 0) V.z = - V.z ;
    if(V.x > V.y) {
        int tmp = V.y ;
        V.y = V.x ;
        V.x = tmp ;
    }
    if(V.y > V.z) {
        int tmp = V.z ;
        V.z = V.y ;
        V.y = tmp ;
        if(V.x > V.y) {
            int tmp = V.y ;
            V.y = V.x ;
            V.x = tmp ;
        }
    }

    return V ;
}

#define IABS(x) ((x)>=0 ? (x) : -(x))
void AddCube(CTX_Cube * CC,V3 T1,V3 T2,int N,int fact) {
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
        int k ;
        int d1 = PGCD(PGCD(IABS(cube.A.x),IABS(cube.A.y)),IABS(cube.A.z)) ;
        int d2 = PGCD(PGCD(IABS(cube.B.x),IABS(cube.B.y)),IABS(cube.B.z)) ;
        int d3 = PGCD(PGCD(IABS(cube.C.x),IABS(cube.C.y)),IABS(cube.C.z)) ;
 //       printf("%dx%dx%d X%d (%d,%d,%d)+(%d,%d,%d)+(%d,%d,%d)\n"
 //              ,sizex,sizey,sizez,fact
 //              ,cube.A.x,cube.A.y,cube.A.z,cube.B.x,cube.B.y,cube.B.z,cube.C.x,cube.C.y,cube.C.z) ;
        CC->nbC++ ;
        for(k=1;k*sizexyz <= PB579_SIZE;k++) {
 //           int64_t nbCubes = (PB579_SIZE-k*sizex+1)*(PB579_SIZE-k*sizey+1)*(PB579_SIZE-k*sizez+1) * fact ;
            int64_t nbCubes =  (((PB579_SIZE-k*sizex+1)*(PB579_SIZE-k*sizey+1)) % PB579_MOD)
            * (((PB579_SIZE-k*sizez+1) * fact) % PB579_MOD) ;
            CC->nbCube += nbCubes ;
            // (k*N+1) (k*k*N*N+ (d1+d2+d3-N)*k +1 ) pour k=1
//            int64_t nbVertices = (k*N+1) * (k*k*N*N+ (d1+d2+d3-N)*k +1 ) ;
            int64_t nbVertices = ((k*N+1) * ( (int64_t)((k*k)% PB579_MOD ) * ((N*N)% PB579_MOD ) + (( (d1+d2+d3-N)*k )  % PB579_MOD) + 1)) % PB579_MOD ;
 //           printf("k=%d %dc x %lldv=%lld\n,",k,nbCubes,nbVertices,nbCubes * nbVertices) ;
            CC->nbVertice = ( CC->nbVertice+((nbCubes * nbVertices) % PB579_MOD)) % PB579_MOD ;
       }
    }

    
}


static inline V3 R1(V3 T) {
    V3 Tn   ; Tn.x = T.y    ; Tn.y = T.z    ; Tn.z = T.x    ; return Tn ;
}
static inline V3 R2(V3 T) {
    V3 Tn   ; Tn.x = T.z    ; Tn.y = T.x    ; Tn.z = T.y    ; return Tn ;
}
static inline V3 Pxy(V3 T) {
    V3 Tn   ; Tn.x = T.y   ; Tn.y = T.x    ; Tn.z = T.z    ; return Tn ;
}
static inline V3 Pyz(V3 T) {
    V3 Tn   ; Tn.x = T.x   ; Tn.y = T.z    ; Tn.z = T.y    ; return Tn ;
}
static inline V3 Pxz(V3 T) {
    V3 Tn   ; Tn.x = T.z   ; Tn.y = T.y    ; Tn.z = T.x    ; return Tn ;
}



static inline V3 Sx(V3 T) {
    V3 Tn   ; Tn.x = -T.x   ; Tn.y = T.y    ; Tn.z = T.z    ; return Tn ;
}
static inline V3 Sy(V3 T) {
    V3 Tn   ; Tn.x = T.x    ; Tn.y = -T.y   ; Tn.z = T.z    ; return Tn ;
}
static inline V3 Sz(V3 T) {
    V3 Tn   ; Tn.x = T.x    ; Tn.y = T.y    ; Tn.z = -T.z   ; return Tn ;
}



static inline int32_t SC (V3 S, V3 T) {
    return S.x*T.x + S.y*T.y + S.z*T.z ;
}

void MATCH(CTX_Cube *CC,V3GCD *tbTrip,int nbT, int i1, int i2,int N,V3 T1,V3 T2s, V3 T3) {
    int i3 ;
    for(i3=0;i3<nbT;i3++) {
        if(memcmp(&T3,&(tbTrip[i3].T),sizeof(T3)) == 0 ) {
            if(i3>= i2) {
                int fact ;
                printf("[%d,%d,%d]",i1,i2,i3) ;
                if(i1==i2) {
                    if(i3==i2) {
                        fact = (N==3) ? 4 : 8 ;
                    } else {//
                        fact = 12 ;
                    }
                } else if(i2==i3) {
                    fact = 12 ;
                } else {
                    fact = 24 ;
                }
                AddCube(CC,T1, T2s, N, fact);
            } else {
                printf("[*%d,%d,%d]",i1,i2,i3) ; AddCube(CC,T1, T2s, N, 11);
            }
            break ;
        }
    }
    if(i3==nbT) {
       printf("ERROR");
    }
}

void VerifCube(CTX_Cube *CC,V3GCD *tbTrip,int nbT, int i1, int i2,int N) {
    if(PGCD(tbTrip[i1].gcd,tbTrip[i2].gcd) > 1) return ;
    V3 T1 = tbTrip[i1].T ;
    V3 T2 = tbTrip[i2].T ;
    V3 T2r,T2s ;
    V3 T3 ;
    int i3 ;
/*    if(tbTrip[i1].cl >= V0ab) {
        if(tbTrip[i2].cl < V0ab ) return  ;
        if(tbTrip[i1].cl == tbTrip[i2].cl ) return  ;
        if(tbTrip[i1].cl == V0ab) {
            AddCube(CC,R1(T1),T2,N,6); // pas de test a faire cela est OK
        } else {
            AddCube(CC,R1(T2),T1,N,6);
        }
    } else if(tbTrip[i2].cl >= V0ab) {
        return ; // T1 n'est pas V0ab ou V00c
    }  else
*/
    if(tbTrip[i1].cl == V00c) {
        if(tbTrip[i2].cl == V0ab ) {
            AddCube(CC,R1(T2),T1,N,6);
        }
        return ;
    } else if(tbTrip[i2].cl == V00c) {
        if(tbTrip[i1].cl == V0ab ) {
            AddCube(CC,R1(T1),T2,N,6); // pas de test a faire cela est OK
        }
        return ;
    }
    {
        T2s = T2 ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        
        T2r = R1(T2);
        T2s = T2r ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}

        T2r = R2(T2);
        T2s = T2r ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}

        T2r = Pxy(T2);
        T2s = T2r ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        
        T2r = Pyz(T2);
        T2s = T2r ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}

        T2r = Pxz(T2);
        T2s = T2r ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sx(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sy(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}
        T2s =Sz(T2r) ; if( SC(T2s,T1) == 0) { T3 = PVectAbs(T1,T2s,N) ; goto MATCH ;}


    }
    return ;
    
MATCH:
    for(i3=0;i3<nbT;i3++) {
        if(memcmp(&T3,&(tbTrip[i3].T),sizeof(T3)) == 0 ) {
            if(i3>= i2) {
                int fact ;
//                printf("[%d,%d,%d]",i1,i2,i3) ;
                if(i1==i2) {
                    if(i3==i2) {
                        fact = (N==3) ? 4 : 8 ;
                    } else {//
                        fact = 12 ;
                    }
                } else if(i2==i3) {
                    fact = 12 ;
                } else {
                    fact = 24 ;
                }
                AddCube(CC,T1, T2s, N, fact);
            }
            break ;
        }
    }
    if(i3==nbT) {
        printf("ERROR");
    }
}

int PB579(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_Cube CC ;
    V3GCD tbTrip[PB579_MAXT] ;
    u_int32_t n1,n2,n3;
    u_int32_t N ;
    CC.nbC = 0 ;
    // on ajoute les cubes alignes sur le reseau
    N = PB579_SIZE ;
//    CC.nbVertice =(((u_int64_t)(N))*(N+1)*(3 * N*N*N*N*N + 39 * N*N*N*N + 213 * N*N*N + 627 * N*N + 640 * N + 158 ))/ 420 ;
    int k ;
    CC.nbVertice = 0 ;
    for(k=1;k<=N;k++) {
        CC.nbVertice = (CC.nbVertice + ( (((N+2-k)*(N+2-k)*(N+2-k)) % PB579_MOD) * ((k*k*k) % PB579_MOD) )) % PB579_MOD ;
    }
    CC.nbCube = (((int64_t)N)*N * (N+1) * (N+1)) / 4 ;

    CC.nb8 = CC.nb12 = CC.nb24 = 0L ;
    for(N=3;N<=PB579_SIZE;N+=2) {
        u_int32_t S1 = N*N ;
        int nbT = 0 ;
        for(n1=0;n1<N;n1++) {
            u_int32_t S2 = S1 - n1*n1 ;
            int n2Max =Sqrt32(S2>>1) ;
            for(n2=n1; n2<=n2Max;n2++) {
                n3 =Sqrt32(S2-n2*n2) ;
                if(n3*n3 == S2-n2*n2) {
                    int gcd  = PGCD(n1,n2);
                    tbTrip[nbT].T.x = n1 ;
                    tbTrip[nbT].T.y = n2 ;
                    tbTrip[nbT].T.z = n3 ;
                    tbTrip[nbT].gcd = PGCD(gcd,n3); ;
/*                    if(n2==n3) { // dans ce cas n1 != 0
                        tbTrip[nbT].T.x = n3 ;
                        tbTrip[nbT].T.z = n1 ;
                    }
 */                   ClassV cl = Vabc ;
                    if(tbTrip[nbT].T.x==0 ) {
                        if(tbTrip[nbT].T.y==0) cl = V00c ;
                        else cl = V0ab ;
                    }
                    tbTrip[nbT].cl = cl ;
                    nbT++;
                }
            }
        }
        if(nbT > 0) {
            int i1,i2 ;
            printf("\n%d ",N);
            for(i1=0;i1<nbT;i1++){
                for(i2=i1;i2<nbT;i2++) {
                    VerifCube(&CC,tbTrip,nbT,i1,i2,N) ;
                }
            }
        }
    }
    printf("NBcand=%d NBC=%d\n",CC.nbCand, CC.nbC) ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d NbCubes=%d NBVertices=%lld V8=%lld V12=%lld V24=%lld\n"
                              ,pbR->pbNum,CC.nbCube,CC.nbVertice,CC.nb8,CC.nb12,CC.nb24) ;
    sprintf(pbR->strRes,"%lld",CC.nbVertice);
    
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



#define PB1000_NUM  10000000
int PB1000(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrimeNb(PB1000_NUM)) == NULL) {
        fprintf(stdout,"\t PB%d Fail to alloc prime table\n",pbR->pbNum);
        return 0 ;
    }
    const T_prime * tbPrime = GetTbPrime(ctxP) ;
    
    fprintf(stdout,"\t PB%0.4d prime[%d] = %d\n",pbR->pbNum, PB1000_NUM,tbPrime[PB1000_NUM-1]) ;
    sprintf(pbR->strRes,"%d",tbPrime[PB1000_NUM-1]);
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

