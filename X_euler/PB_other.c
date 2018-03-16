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

typedef enum ClassV {
    Vabc = 0 , Vaab = 1 , V0ab = 2 , V00c = 3
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
} CTX_Cube ;

static inline V3 PVect( V3 T1,V3 T2,int N) {
    V3 V ;
    V.x = (T1.y*T2.z - T1.z*T2.y)/N ;
    V.y = (T1.z*T2.x - T1.x*T2.z)/N ;
    V.z = (T1.x*T2.y - T1.y*T2.x)/N ;
    return V ;
}
void AddCube(CTX_Cube * CC,V3 T1,V3 T2,int N,int fact) {
    cubeO cube ;
    cube.A = T1 ;
    cube.B = T2 ;
    cube.C = PVect(T1,T2,N) ;
    CC->nbCand++ ;
    printf("%dx%d (%d,%d,%d)+(%d,%d,%d)+(%d,%d,%d)\n"
           ,N,fact
           ,cube.A.x,cube.A.y,cube.A.z,cube.B.x,cube.B.y,cube.B.z,cube.C.x,cube.C.y,cube.C.z
    ) ;
    CC->nbC++ ;

    
}

void AddCube1(CTX_Cube * CC,int32_t a1,int32_t b1,int32_t c1, int32_t a2,int32_t b2,int32_t c2, int N) {
    cubeO cube ;
    int offx, offy,offz ;
    int sizex, sizey,sizez ;
    cube.A.x = a1 ; cube.A.y = b1 ; cube.A.z = c1 ;
    cube.B.x = a2 ; cube.B.y = b2 ; cube.B.z = c2 ;
    cube.C.x = -(b1*c2 - c1*b2) / N ;
    cube.C.y = -(c1*a2 - a1*c2) / N ;
    cube.C.z = -(a1*b2 - b1*a2) / N ;
    CC->nbCand++ ;
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
    if(sizex <= PB579_SIZE && sizey <= PB579_SIZE && sizez <= PB579_SIZE) {
        cubeMin cubeM ;
        MinimizeCube(&cube,&cubeM);
        printf("%dx%dx%d(%d,%d,%d)+(%d,%d,%d)+(%d,%d,%d)=>[0,%d,%d][%d,0,%d][%d,%d,0]\n"
               ,sizex,sizey,sizez
               ,cube.A.x,cube.A.y,cube.A.z,cube.B.x,cube.B.y,cube.B.z,cube.C.x,cube.C.y,cube.C.z
               ,cubeM.Ay,cubeM.Az,cubeM.Bx,cubeM.Bz,cubeM.Cx,cubeM.Cy) ;
        CC->nbC++ ;
    }
}

#define TST_ORTHO1(ptC,a1,b1,c1,a2,b2,c2,N) if(-a1*a2+b1*b2+c1*c2 == 0) AddCube(ptC,a1,b1,c1,-a2,b2,c2,N)
#define TST_ORTHO2(ptC,a1,b1,c1,a2,b2,c2,N) if(a1*a2-b1*b2+c1*c2 == 0)  AddCube(ptC,b1,a1,c1,-b2,a2,c2,N)
#define TST_ORTHO3(ptC,a1,b1,c1,a2,b2,c2,N) if(a1*a2+b1*b2-c1*c2 == 0)  AddCube(ptC,c1,a1,b1,-c2,a2,b2,N)

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

int PB579(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    CTX_Cube CC ;
    V3GCD tbTrip[PB579_MAXT] ;
        int nbO=0, nbp=0,nbq = 0 ;
    u_int32_t n1,n2,n3 ,N ;
    CC.nbC = 0 ;
    for(N=3;N<=PB579_SIZE;N+=2) {
        u_int32_t S1 = N*N ;
        int nbT = 0 ;
        for(n1=0;n1<N;n1++) {
            u_int32_t S2 = S1 - n1*n1 ;
            int n2Max =Sqrt32(S2>>1) ;
            for(n2=n1;n2<=n2Max;n2++) {
                n3 =Sqrt32(S2-n2*n2) ;
                if(n3*n3 == S2-n2*n2) {
                    int gcd  = PGCD(n1,n2);
                    tbTrip[nbT].T.x = n1 ;
                    tbTrip[nbT].T.y = n2 ;
                    tbTrip[nbT].T.z = n3 ;
                    tbTrip[nbT].gcd = PGCD(gcd,n3); ;
                    if(n2==n3) { // dans ce cas n1 != 0
                        tbTrip[nbT].T.x = n3 ;
                        tbTrip[nbT].T.z = n1 ;
                    }
                    ClassV cl = Vabc ;
                    if(tbTrip[nbT].T.x==tbTrip[nbT].T.y) {
                        if(tbTrip[nbT].T.x==0) cl = V00c ;
                        else cl = Vaab ;
                    } else if (tbTrip[nbT].T.x==0) cl=V0ab ;
                    tbTrip[nbT].cl = cl ;
                    nbT++;
                }
            }
        }
        if(nbT > 0) {
            int i1,i2 ;
            printf("\n%d ",N);
            for(i1=0;i1<nbT;i1++){
                int gcd1 = tbTrip[i1].gcd ;
                V3 T1=tbTrip[i1].T ;
                for(i2=i1;i2<nbT;i2++) {
                    int gcd2 = tbTrip[i2].gcd ;
                    V3 T2 = tbTrip[i2].T ;
                    if(PGCD(gcd1,gcd2) > 1) continue ;
                    if(tbTrip[i1].cl >= V0ab) {
                        if(tbTrip[i2].cl < V0ab ) continue ;
                        if(tbTrip[i1].cl == tbTrip[i2].cl ) continue ;
                        if(tbTrip[i1].cl == V0ab) {
                            AddCube(&CC,R1(T1),T2,N,6); // pas de test a faire cela est OK
                        } else {
                            AddCube(&CC,R1(T2),T1,N,6);
                        }
                    } else if(tbTrip[i2].cl >= V0ab) {
                       continue ; // T1 n'est pas V0ab ou V00c
                    } else if(tbTrip[i1].cl == Vaab) {
                        if(tbTrip[i2].cl == Vaab) {
                            if( SC(T1,T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,T1,T2,N,10) ;
                            if( SC(Sx(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sx(T1),T2,N,2) ;
                            if( SC(Sy(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sy(T1),T2,N,2) ;
                            if( SC(Sz(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sz(T1),T2,N,2) ;
                            
                            if( SC(R1(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,R1(T1),T2,N,4) ;
                            if( SC(Sx(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sx(R1(T1)),T2,N,4) ;
                            if( SC(Sy(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sy(R1(T1)),T2,N,4) ;
                            if( SC(Sz(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sz(R1(T1)),T2,N,4) ;

                        } else {
                            if( SC(T1,T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,T1,T2,N,2) ;
                            if( SC(Sx(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sx(T1),T2,N,2) ;
                            if( SC(Sy(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sy(T1),T2,N,2) ;
                            if( SC(Sz(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sz(T1),T2,N,2) ;
                            
                            if( SC(R1(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,R1(T1),T2,N,2) ;
                            if( SC(Sx(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sx(R1(T1)),T2,N,2) ;
                            if( SC(Sy(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sy(R1(T1)),T2,N,2) ;
                            if( SC(Sz(R1(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sz(R1(T1)),T2,N,2) ;
                            
                            if( SC(R2(T1),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,R2(T1),T2,N,2) ;
                            if( SC(Sx(R2(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sx(R2(T1)),T2,N,2) ;
                            if( SC(Sy(R2(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sy(R2(T1)),T2,N,2) ;
                            if( SC(Sz(R2(T1)),T2) == 0) printf("(%d,%d)",i1,i2),AddCube(&CC,Sz(R2(T1)),T2,N,2) ;
                            
                        }
                        
                    } else if(tbTrip[i2].cl == Vaab) {
                        if( SC(T2,T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,T2,T1,N,2) ;
                        if( SC(Sx(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(T2),T1,N,2) ;
                        if( SC(Sy(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(T2),T1,N,2) ;
                        if( SC(Sz(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(T2),T1,N,2) ;

                        if( SC(R1(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,R1(T2),T1,N,2) ;
                        if( SC(Sx(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(R1(T2)),T1,N,2) ;
                        if( SC(Sy(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(R1(T2)),T1,N,2) ;
                        if( SC(Sz(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(R1(T2)),T1,N,2) ;

                        if( SC(R2(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,R2(T2),T1,N,1) ;
                        if( SC(Sx(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(R2(T2)),T1,N,2) ;
                        if( SC(Sy(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(R2(T2)),T1,N,2) ;
                        if( SC(Sz(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(R2(T2)),T1,N,2) ;
                    }  else {
                        if( SC(T2,T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,T2,T1,N,2) ;
                        if( SC(Sx(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(T2),T1,N,2) ;
                        if( SC(Sy(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(T2),T1,N,2) ;
                        if( SC(Sz(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(T2),T1,N,2) ;
                        
                        if( SC(R1(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,R1(T2),T1,N,2) ;
                        if( SC(Sx(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(R1(T2)),T1,N,2) ;
                        if( SC(Sy(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(R1(T2)),T1,N,2) ;
                        if( SC(Sz(R1(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(R1(T2)),T1,N,2) ;
                        
                        if( SC(R2(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,R2(T2),T1,N,1) ;
                        if( SC(Sx(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(R2(T2)),T1,N,2) ;
                        if( SC(Sy(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(R2(T2)),T1,N,2) ;
                        if( SC(Sz(R2(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(R2(T2)),T1,N,2) ;

                        if( SC(Pxy(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Pxy(T2),T1,N,2) ;
                        if( SC(Sx(Pxy(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(Pxy(T2)),T1,N,2) ;
                        if( SC(Sy(Pxy(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(Pxy(T2)),T1,N,2) ;
                        if( SC(Sz(Pxy(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(Pxy(T2)),T1,N,2) ;

                        if( SC(Pyz(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Pyz(T2),T1,N,2) ;
                        if( SC(Sx(Pyz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(Pyz(T2)),T1,N,2) ;
                        if( SC(Sy(Pyz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(Pyz(T2)),T1,N,2) ;
                        if( SC(Sz(Pyz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(Pyz(T2)),T1,N,2) ;

                        if( SC(Pxz(T2),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Pxz(T2),T1,N,2) ;
                        if( SC(Sx(Pxz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sx(Pxz(T2)),T1,N,2) ;
                        if( SC(Sy(Pxz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sy(Pxz(T2)),T1,N,2) ;
                        if( SC(Sz(Pxz(T2)),T1) == 0) printf("(%d,%d)",i2,i1),AddCube(&CC,Sz(Pxz(T2)),T1,N,2) ;

                        
                    }

                
//                    if( SC(T1,Sx(R1(T2))) == 0) AddCube(&CC,T1,Sx(R1(T2)),N) ;
/*
                    TST_ORTHO1(&CC,T1.x,T1.y,T1.z,T2.x,T2.y,T2.z,N) ;
                    if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.x,T1.y,T1.z,T2.x,T2.y,T2.z,N) ;
                    TST_ORTHO3(&CC,T1.x,T1.y,T1.z,T2.x,T2.y,T2.z,N) ;
                    if(T1.y != T1.z)
                    {
                        TST_ORTHO1(&CC,T1.x,T1.z,T1.y,T2.x,T2.y,T2.z,N) ;
                        if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.x,T1.z,T1.y,T2.x,T2.y,T2.z,N) ;
                        TST_ORTHO3(&CC,T1.x,T1.z,T1.y,T2.x,T2.y,T2.z,N) ;
                    }
                    if(T1.x != T1.y) {
                        TST_ORTHO1(&CC,T1.y,T1.x,T1.z,T2.x,T2.y,T2.z,N) ;
                        if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.y,T1.x,T1.z,T2.x,T2.y,T2.z,N) ;
                        TST_ORTHO3(&CC,T1.y,T1.x,T1.z,T2.x,T2.y,T2.z,N) ;
                    }
                    if(T1.x != T1.y) {
                        TST_ORTHO1(&CC,T1.y,T1.z,T1.x,T2.x,T2.y,T2.z,N) ;
                        if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.y,T1.z,T1.x,T2.x,T2.y,T2.z,N) ;
                        TST_ORTHO3(&CC,T1.y,T1.z,T1.x,T2.x,T2.y,T2.z,N) ;
                    }
                    if(T1.y != T1.z)
                    {
                        TST_ORTHO1(&CC,T1.z,T1.x,T1.y,T2.x,T2.y,T2.z,N) ;
                        if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.z,T1.x,T1.y,T2.x,T2.y,T2.z,N) ;
                        TST_ORTHO3(&CC,T1.z,T1.x,T1.y,T2.x,T2.y,T2.z,N) ;
                    }
//                    if(T1.x != T1.y && T1.y != T1.z)
                    if(T1.y != T1.z )
                    {
                        TST_ORTHO1(&CC,T1.z,T1.y,T1.x,T2.x,T2.y,T2.z,N) ;
                        if(T2.x != T2.y && T2.y != T2.z) TST_ORTHO2(&CC,T1.z,T1.y,T1.x,T2.x,T2.y,T2.z,N) ;
                        TST_ORTHO3(&CC,T1.z,T1.y,T1.x,T2.x,T2.y,T2.z,N) ;
                     }
*/
                    
                }
            }
        }
    }
    printf("NBcand=%d NBC=%d\n",CC.nbCand, CC.nbC) ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%0.3d Nb=%d (O:%d,P:%d,Q:%d)\n",pbR->pbNum,nbO+nbp+nbq, nbO,nbp,nbq) ;
    sprintf(pbR->strRes,"%d",nbO+nbp+nbq);
    
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

