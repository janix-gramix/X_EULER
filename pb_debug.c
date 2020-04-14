
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


#define PB702_MAX  15
#define PB702_NBT   12345

/*
#define PB702_MAX  5
#define PB702_NBT   17
*/

typedef struct PB702_POINT {
    uint32_t ix ;
    uint32_t iy ;
} PB702_POINT ;

int CmpPoint(const void * el1, const void *el2) {
    const PB702_POINT * p1 =(PB702_POINT *)el1 ;
    const PB702_POINT * p2 =(PB702_POINT *)el2 ;
    int diff = p1->ix - p2->ix ;
    if(diff) return diff ;
    return p1->iy - p2->iy ;
}

int PB702(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nl= PB702_NBT ;
    int64_t ixMax = 1 << (PB702_MAX+1) ;
    int64_t ixMax2 = ixMax >> 1 ;
    int64_t iyMax = 1 << PB702_MAX ;
    
    uint8_t *isTok = malloc((nl+1)*sizeof(isTok[0]));
    int64_t sum = 0 ;
//    for(int il=0;il<nl;il++) {
    int ikMax = 0 ;
    int64_t *histoB = calloc(32,sizeof(histoB[0]));
    int64_t *histoBT = calloc(32,sizeof(histoB[0]));
//************** standard *******************
    for(int il=0;il<nl;il++) {
        int64_t nbt = nl-(il+1)/2 ;
        int64_t halfNbt = il & 1 ;
        int64_t isHalfAlive = halfNbt ;
        int64_t nbTalive = nbt  ;
        memset(isTok,0,nbTalive *sizeof(isTok[0])) ;
        int ik ;
        for(ik = 2;nbTalive+isHalfAlive;ik++) {
            int64_t kpow2 = 1 << ik ;
            int64_t ilklow = (il * kpow2)/ nl +1 ;
             int64_t ilkhigh = ((il+1)*kpow2) / nl +1 ;
            for(int64_t ilk = ilklow;ilk<ilkhigh;ilk++) {
                int64_t  idebT,deltaK, deltaT ; //  kpow2 * nl
                idebT = (nl * ilk  - il * kpow2) ;
                deltaT = 2 * (kpow2 - idebT) ;
                if(ilk & 1) {
                    idebT += nl;
                    deltaK =  2* nl ;
                } else {
                    if((ilk & 2) == 0) {
                        idebT += 2*nl ;
                    }
                    deltaK = 4 * nl ;
                }
                if(halfNbt) {
                   // on traite le cas particulier du demi triangle
                     if(isHalfAlive)
                     {
                       int64_t ifinTH = ((ilk & 3) ==2) ? idebT+deltaT : idebT+deltaT - deltaK ;
                         if (ifinTH - kpow2 >=0) {
                             int64_t ikfin = (ifinTH- kpow2)/ deltaK ;
                             if((ifinTH- kpow2)  == deltaK *ikfin) ikfin-- ;
                             if(ikfin >=0) {
                                 histoB[ik]++ ;
 //                                if(isHalfAlive)
                                 {
                                     isHalfAlive = 0 ;
                                     sum += ik ;
                                 }
                             }
                         }
                     }

                     idebT += kpow2 ;
                }
              if(deltaT > deltaK) {
                 sum += nbTalive * 2 * ik ;
                   histoB[ik] += 2*nbTalive  ;
                  if((ilk & 1) == 0 ) histoBT[ik] += 2*nbTalive  ;
                   nbTalive = 0 ;
                   break ;
               } else { // deltaT <= deltaK
                   int64_t modDebT = idebT % deltaK ;
                   for(int nT = 0; nT < nbt ; /*nT++ */) {
                       if(isTok[nT] == 0)
                       {
                           int64_t diffModT = modDebT + deltaT ;
                           if(modDebT == 0) diffModT -= deltaK ;
                           if(diffModT > deltaK)   {
 //                              if(isTok[nT] == 0)
                               {
                                   isTok[nT] = 1 ;
                                   sum += 2*ik;
                                   --nbTalive ;
                               }
                               histoB[ik] += 2 ;
                               if((ilk & 1) == 0 ) histoBT[ik] += 2 ;
                               
                           }
                       }
                       int64_t dnT = (deltaK-deltaT-modDebT) >> (ik+1) ;
                       if(dnT > 0) {
                           nT += dnT ;
                           if(nT>=nbt) break ;
                           modDebT += dnT * 2 * kpow2 ;
                       } else {
                           modDebT += 2 * kpow2 ;
                           nT++ ;
                       }
                
                       
                       while(modDebT >= deltaK) modDebT -= deltaK ;
                       //                       ifinT += 2 * kpow2 ;
                   }
               }
            }
        }
        if(ik > ikMax) ikMax = ik ;
        
    }
    printf("Sum=%lld ikMax=%d\n",sum,ikMax);
    //
#if 0
//********************** Half **********************
    int64_t sumH = 0 ;
    int64_t *histoH = calloc(32,sizeof(histoH[0]));
    for(int il=0;il<nl;il++) {
        int64_t nbt = nl-(il+1)/2 ;
        int64_t halfNbt = il & 1 ;
        if(halfNbt == 0) continue ;
        int64_t isHalfAlive = halfNbt ;
         int ik ;
        for(ik = 2;isHalfAlive;ik++) {
            int64_t kpow2 = 1 << ik ;
            int64_t ilklow = (il * kpow2)/ nl +1 ;
            int64_t ilkhigh = ((il+1)*kpow2) / nl +1 ;
            for(int64_t ilk = ilklow;ilk<ilkhigh;ilk++) {
                int64_t  idebT,deltaK,deltaT ; //  kpow2 * nl
                idebT = (nl * ilk  - il * kpow2) ;
                //               histoDelta[kpow2+idebT]++ ;
                deltaT = 2 * (kpow2 - idebT) ;
                if(ilk & 1) {
                    idebT += nl;
                    deltaK =  2* nl ;
                } else {
                    if((ilk & 2) == 0) {
                        idebT += 2*nl ;
                    }
                    deltaK = 4 * nl ;
                }
                // on traite le cas particulier du demi triangle
                if(isHalfAlive) {
                    int64_t ifinTH = ((ilk & 3) ==2) ? idebT+deltaT : idebT+deltaT - deltaK ;
                    if (ifinTH - kpow2 >=0) {
                        int64_t ikfin = (ifinTH- kpow2)/ deltaK ;
                        if((ifinTH- kpow2)  == deltaK *ikfin) ikfin-- ;
                        if(ikfin >=0) {
                            histoH[ik]++ ;
                              isHalfAlive = 0 ;
                            sumH += ik ;
                        }
                    }
                }
            }
        }
        if(ik > ikMax) ikMax = ik ;
    }
    
//***********************  Triangle half (sans le Half) + avec somme total => blanc)
    int64_t sumTB = 0 ;
    int64_t *histoTB = calloc(32,sizeof(histoTB[0]));
    for(int il=0;il<nl;il++) {
        int64_t nbt = (il)/2 ;
        int64_t halfNbt = il & 1 ;
        int64_t isHalfAlive = halfNbt ;
        int64_t nbTalive = nbt  ;
        memset(isTok,0,nbTalive *sizeof(isTok[0])) ;
        int ik ;
        for(ik = 2;nbTalive;ik++) {
            int64_t kpow2 = 1 << ik ;
            int64_t ilklow = (il * kpow2)/ nl +1 ;
            int64_t ilkhigh = ((il+1)*kpow2) / nl +1 ;
            for(int64_t ilk = ilklow;ilk<ilkhigh;ilk++) {
                int64_t  idebT,deltaK,deltaT ; //  kpow2 * nl
                idebT = (nl * ilk  - il * kpow2) ;
                 deltaT = 2 * (kpow2 - idebT) ;
                if(ilk & 1) {
                    idebT += nl;
                    deltaK =  2* nl ;
                } else {
                    if((ilk & 2) == 0) {
                        idebT += 2*nl ;
                    }
                    deltaK = 4 * nl ;
                }
                if(halfNbt) {
                    idebT += kpow2 ;
                }
                if(deltaT > deltaK) {
                    sumTB += nbTalive * ik ;
                    histoTB[ik] += nbTalive  ;
 //                   if(ik==4) printf("(il=%d,ilk=%lld,ali=%lld)",il,ilk,nbTalive);
                    nbTalive = 0 ;
                    break ;
                } else { // deltaT <= deltaK
                    int64_t modDebT = idebT % deltaK ;
                    for(int nT = 0; nT < nbt ; nT++ ) {
                        if(isTok[nT] == 0) {
                            int64_t diffModT = modDebT + deltaT ;
                            if(modDebT == 0) diffModT -= deltaK ;
                            if(diffModT > deltaK)   {
                                isTok[nT] = 1 ;
                                sumTB += ik;
                                --nbTalive ;
                                histoTB[ik] ++  ;
   //                             if(ik==4) printf("(il=%d,ilk=%lld,nT=%d)",il,ilk,nT);
                            }
                        }
                        modDebT += 2 * kpow2 ;
                        while(modDebT >= deltaK) modDebT -= deltaK ;
                        //                       ifinT += 2 * kpow2 ;
                    }
                }
            }
        }
        if(ik > ikMax) ikMax = ik ;
        
    }

 //***********************  Triangle half (sans le Half) + avec somme total => blanc calcul par les ik
    

    
    int64_t sumTBik = 0 ;
    int64_t *histoTBik = calloc(32,sizeof(histoTBik[0]));
 
    for(int ik = 2;ik<16;ik++) {
        int64_t kpow2 = 1 << ik ;
        int64_t mod2nl = (2 * kpow2) % (2 * nl) ;
        int64_t mod4nl = (2 * kpow2) % (4 * nl) ;

        int lastIl = 0 ;
        for(int ikl=1;ikl<kpow2;ikl++) {
            int64_t il =  (ikl * nl) / kpow2 ;
            if(il==lastIl) {
                continue;
            }
            lastIl = il ;
            int64_t nbt = (il)/2 ;
            int64_t halfNbt = il & 1 ;
            int64_t isHalfAlive = halfNbt ;
            int64_t nbTalive = nbt  ;
            int64_t  idebT,deltaK,deltaT ; //  kpow2 * nl
    //        idebT = (nl * ikl  - il * kpow2) ;
            idebT = il * kpow2 ;
            int64_t offSetDeb = (nl * ikl  - il * kpow2) ;
            if(offSetDeb ==0) continue ;
            deltaT = 2 * (kpow2 - offSetDeb) ;
            int64_t modnl = mod2nl ;
            if(ikl & 1) {
  //              idebT += nl;
                deltaK =  2* nl ;
            } else {
                if((ikl & 2) == 0) {
//                    idebT += 2*nl ;
                    idebT -= 2*nl ;
                }
                deltaK = 4 * nl ;
                modnl = mod4nl ;
            }
            if(halfNbt) {
 //               idebT += kpow2 ;
            }
            if(deltaT > deltaK) {
                sumTBik += nbTalive * ik ;
                histoTBik[ik] += nbTalive  ;
                nbTalive = 0 ;
//                break ;
            } else { // deltaT <= deltaK
                int64_t modDebT = idebT % deltaK ;
                for(int nT = 0; nT < nbt ; /*nT++ */) {
                    int64_t diffModT = modDebT - deltaT ;
                    if(modDebT == 0) diffModT += deltaK ;
                    if(diffModT > deltaK)   {
                        sumTBik += ik;
                        --nbTalive ;
                        histoTBik[ik] ++  ;
                    }
                    int64_t dnT = (deltaK-deltaT-modDebT) >> (ik+1) ;
                    if(dnT > 0) {
                        nT += dnT ;
                        if(nT>=nbt) break ;
                       modDebT += dnT * modnl ;
                    } else {
                       modDebT += modnl  ;
                        nT++ ;
                    }
                    
                    
                    
                    while(modDebT >= deltaK) modDebT -= deltaK ;
                    //                       ifinT += 2 * kpow2 ;
                }
            }
        }
    }
    
#endif
  
    printf("%.3fs\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC);
    //***********************  Triangle half inferieurx (sans le Half) + avec somme total => blanc calcul par les il
    
    
    
    int64_t sumTBil = 0 ;
    int64_t *histoTBil = calloc(32,sizeof(histoTBil[0]));
    int64_t *histoBTil = calloc(32,sizeof(histoBTil[0]));
    int ikm ;
    for(ikm=0;(1<<ikm) <nl;ikm++) ;
    ikm++ ;
    int64_t kmpow2  = 1 << ikm ;
    int32_t maskOdd[32] ,maskEven[32];
    int32_t im = 0xffffffff ;
    for(int ik=0; ik < ikm ;ik++ ) {
        maskOdd[ik] = (1 << (ik+1))-1 ;
        maskEven[ik] = (1 << (ik+2))-1 ;
    }
    
    
    for(int ik = 2;ik<ikm;ik++) {
        int kpow2 = 1 << ik ;
        for(int ilk=1; ilk<kpow2;ilk++) {
            int il = nl*ilk / kpow2 ;
            if(nl*ilk == il * kpow2) continue ;
            int64_t nbt = nl-(il+1)/2 ;
            int64_t halfNbt = il & 1 ;
 //           int64_t isHalfAlive = halfNbt ;
 //           int64_t nbTalive = nbt  ;
//            memset(isTok,0,nbTalive *sizeof(isTok[0])) ;
             int64_t  idebT,deltaK, deltaT ; //  kpow2 * nl
            idebT = (nl * ilk  - il * kpow2) ;
            deltaT = 2 * (kpow2 - idebT) ;
            if(ilk & 1) {
                idebT += nl;
                deltaK =  2* nl ;
            } else {
                if((ilk & 2) == 0) {
                    idebT += 2*nl ;
                }
                deltaK = 4 * nl ;
            }
            if(halfNbt) {
                // on traite le cas particulier du demi triangle
 //               if(isHalfAlive)
                {
                    int64_t ifinTH = ((ilk & 3) ==2) ? idebT+deltaT : idebT+deltaT - deltaK ;
                    if (ifinTH - kpow2 >=0) {
                        int64_t ikfin = (ifinTH- kpow2)/ deltaK ;
                        if((ifinTH- kpow2)  == deltaK *ikfin) ikfin-- ;
                        if(ikfin >=0) {
                            histoTBil[ik]++ ;
                            //                                if(isHalfAlive)
                            {
  //                              isHalfAlive = 0 ;
                                sum += ik ;
                            }
                        }
                    }
                }
                
                idebT += kpow2 ;
            }
            if(deltaT > deltaK) {
                //               sum += nbTalive * 2 * ik ; ;
 //               histoB[ik] += 2*nbTalive  ;
 //               nbTalive = 0 ;
                histoTBil[ik] += 2 * nbt ;
                if((ilk & 1) == 0)histoBTil[ik] += 2 * nbt ;
  //              printf("T>K!! ") ;
                continue ;
            } else { // deltaT <= deltaK
                int64_t modDebT = idebT % deltaK ;
                for(int nT = 0; nT < nbt ; /*nT++ */) {
 //                   if(isTok[nT] == 0)
                    {
                        int64_t diffModT = modDebT + deltaT ;
                        if(modDebT == 0) diffModT -= deltaK ;
                        if(diffModT > deltaK)   {
                            //                              if(isTok[nT] == 0)
                            {
  //                              isTok[nT] = 1 ;
                                sum += 2*ik;
                              }
                            histoTBil[ik] += 2 ;
                            if((ilk & 1) == 0) histoBTil[ik] += 2 ;
                            
                        }
                    }
                    int64_t dnT = (deltaK-deltaT-modDebT) >> (ik+1) ;
                    if(dnT > 0) {
                        nT += dnT ;
                        if(nT>=nbt) break ;
                        modDebT += dnT * 2 * kpow2 ;
                    } else {
                        modDebT += 2 * kpow2 ;
                        nT++ ;
                    }
                    
                    
                    while(modDebT >= deltaK) modDebT -= deltaK ;
                    //                       ifinT += 2 * kpow2 ;
                }
            }
        }
        if(ik > ikMax) ikMax = ik ;
        
    }
    printf("Sum=%lld ikMax=%d\n",sum,ikMax);

    
    
    int64_t SumByHistoB = 0 ;
    int64_t SumHand2TB = 0 ;
    int64_t SumBandTB = 0 ;
    int64_t SumTil = 0 ;
    int64_t nbT = 0 ;
    for(int il=0;il<nl;il++) {
        int64_t nbt = nl-(il+1)/2 ;
        int64_t halfNbt = il & 1 ;
        nbT += halfNbt + 2 * nbt ;

    }
    for(int ik=2;ik<ikMax;ik++){
        SumByHistoB += (histoB[ik])*ik ;
        if(histoTBil[ik]) {
            SumTil += histoTBil[ik]*ik ;
            nbT -= histoTBil[ik] ;
        } else {
            SumTil += nbT * ik ;
        }
 //       SumHand2TB +=  2*histoTB[ik]+histoH[ik];
//        int64_t exp = ((1LL<<ik)-1)* ((1LL<<ik)-2) - SumBandTB - 2*histoTB[ik]-histoH[ik] ;
//        SumBandTB += exp + 2*histoTB[ik]+histoH[ik];
        printf("%d=>(B=%lld/Bil=%lld, Triangl(%lld/%lld)  SumB=%lld SumTil=%lld) ",
               ik/*,histoH[ik]*/,histoB[ik],histoTBil[ik] , histoBT[ik],histoBTil[ik],
            SumByHistoB,SumTil) ;
    }
    printf("SumByB=%lld SumBandTB=%lld\n",SumByHistoB,SumBandTB);
    int64_t  antSum = 0 , SumC = 0 ;
 /*
    for(int ik=2;ik<=ikMax;ik++) {
        printf("%d=>(H=%lld,TB=%lld,B=%lld, B+2*T+H=%lld Exp=%lld) ",
               ik,histoH[ik],histoTB[ik],histoB[ik],histoB[ik]+2*histoTB[ik]+histoH[ik],
               (1LL<<(ik-1))* ((1LL<<ik)-1)) ;

    }
  */
    printf("SumC=%lld \n",SumC);
    pbR->nbClock = clock() - pbR->nbClock ;
 //   exit(0);
    return 1 ;
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
#define PB705_W 5
#define PB707_MAXF  7

int PB707(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t S = 0 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    int F0 =0 , F1 = 1 ;
    int w = PB705_W ;
    for(int jf=0;jf<PB707_MAXF;jf++) {
        if((w==F1) || (w & 3)==0 || (F1 & 3) == 0) {
            S += modPow(2,w*F1,1000000007) ;
        } else if (w<F1){
            S += modPow(2,w*(F1-1),1000000007) ;
        } else {
            S += modPow(2,(w-1)*F1,1000000007) ;

        }
        int tmp = F0 ;
        F0 = F1 ; F1 += tmp ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    return 1 ;
}

#define PB708_MAX   100000000
//#define PB708_MAX   17897130
#define PB708_MAX   100000000
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
        printf("h[%d]=%d ",i2,histPow2[i2]);

        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}




typedef uint64_t T_prime64 ;
typedef int(*TY_CPL_nxtPrime64)(void *ctx,T_prime64 nxtPrime);
uint32_t FindPrime64(T_prime64 maxValue,void *ctx,TY_CPL_nxtPrime64 nxtPrime) ;
typedef struct CTX_PRIMETABLE64 CTX_PRIMETABLE64 ;

struct  CTX_PRIMETABLE64 {
    uint32_t       nbPrime ;
    T_prime64     maxValue ;
    uint32_t    maxNbPrime ;
    T_prime64     *tbPrime ;
}  ;

int CPL_tablePrime64(void *ptCtx,T_prime64 nxtPrime) {
    CTX_PRIMETABLE64 * ctx = (CTX_PRIMETABLE64 *) ptCtx ;
    if(nxtPrime > ctx->maxValue) return 0 ;
    ctx->tbPrime[ctx->nbPrime] = nxtPrime ;
    ctx->nbPrime++ ;
    if(ctx->nbPrime >= ctx->maxNbPrime) {
        ctx->maxValue = nxtPrime * nxtPrime ;
        return 0 ;
    }
    return 1;
}

const T_prime64 * GetTbPrime64(CTX_PRIMETABLE64 * ctx) {
    return ctx->tbPrime ;
}

uint32_t GetNbPrime64(CTX_PRIMETABLE64 * ctx) {
    return ctx->nbPrime ;
}

CTX_PRIMETABLE64 * Free_tablePrime64(CTX_PRIMETABLE64 * ctx) {
    if(ctx != NULL) free(ctx->tbPrime) ;
    free(ctx);
    return NULL ;
}



//
// On replie la table.
// la taille de la table peut être quelconque
// Plus rapide pour taille nSqrt ou 32368 si trop grande valeurs
//
uint32_t FindPrime64_(T_prime64 nbMax,void *ctx,TY_CPL_nxtPrime64 nxtPrime) {
    uint32_t nSqrt = 1+ (int32_t)Sqrt64(nbMax ) ;
    int isEnd = 0;
    T_prime64 *tbPrime = malloc(nSqrt * sizeof(tbPrime[0])) ;
    int64_t *offSet = malloc(nSqrt * sizeof(offSet[0])) ;
    uint32_t sizeTable =  (nSqrt < 65536) ? nSqrt : 65536 ;
    int8_t *isComposed = calloc( sizeTable , sizeof(isComposed[0])) ;
    uint32_t nbPrime = 0 ;
    T_prime64 lastPrime = 0 ;
    T_prime64 offSetTable = 0 ;
    T_prime64 nbPrimeSqrt = 0 ;
    
    nbPrime ++ ;
    lastPrime = 2 ;
    if(nxtPrime(ctx,2) == 0) {
        return nbPrime ;
    }
    // remarque le nb premier 2 n'est pas stocke dans tbPrime car pas utilise dans le crible erasto.
    // pour commencer a 3 donc indPrime = (3>>1)
    T_prime64 indPrime = 1 ;
    
    while ( 1) {
        T_prime64 icp ;
        for(icp=indPrime - offSetTable ;icp < sizeTable; icp++ ) {
            if(!isComposed[icp] ) {
                lastPrime = 2*(icp+offSetTable)+1 ;
                
                if(lastPrime < nSqrt) {
                    T_prime64 icp2 ;
                    tbPrime[nbPrimeSqrt] = lastPrime ;
                    for(icp2 = icp + lastPrime; icp2 < sizeTable ; icp2 += lastPrime ) {
                        isComposed[icp2] = 1 ;
                    }
                    offSet[nbPrimeSqrt++] = icp2 - sizeTable ;
                }
                nbPrime++ ;
                if(nxtPrime(ctx,lastPrime) == 0) {
                    isEnd = 1;
                    break ; // demande d'arret
                }
            }
        }
        offSetTable += sizeTable ;
        if(isEnd || offSetTable >= nbMax) break ;
        indPrime = offSetTable ;
        if ( offSetTable + sizeTable > nbMax) { sizeTable = (int32_t)(nbMax - offSetTable) ; }
        memset(isComposed,0,sizeTable) ;
        {
            int np ;
            for(np=0;np<nbPrimeSqrt;np++) {
                T_prime64 p = tbPrime[np] ;
                int64_t indPrime = offSet[np] ;
                while ( indPrime < sizeTable) {
                    isComposed[indPrime] = 1 ;
                    indPrime += p ;
                }
                offSet[np] = indPrime - sizeTable ;
            }
        }
    }
    free(tbPrime);
    free(offSet) ;
    free(isComposed);
    return nbPrime ;
}

uint32_t FindPrime64(T_prime64 maxValue,void *ctx,TY_CPL_nxtPrime64 nxtPrime) {
    return FindPrime64_(maxValue,ctx,nxtPrime) ;
}

CTX_PRIMETABLE64 * Gen_tablePrime64(T_prime64 maxValue) {
    CTX_PRIMETABLE64 *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    ctx->maxValue = maxValue ;
    if(maxValue > 100) {
        ctx->maxNbPrime = (uint32_t) (1+ maxValue / (log((double)maxValue) - 4)) ;
    } else {
        ctx->maxNbPrime = 30 ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime64(ctx) ; }
    FindPrime64(maxValue,ctx,CPL_tablePrime64) ;
    return ctx ;
}

CTX_PRIMETABLE64 * Gen_tablePrimeNb64(T_prime64 maxNb) {
    CTX_PRIMETABLE64 *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    if(ctx == NULL) return ctx ;
    if(maxNb < 30)  {
        ctx->maxNbPrime = 30 ;
    } else {
        ctx->maxNbPrime = (uint32_t) maxNb ;
    }
    ctx->tbPrime = malloc(ctx->maxNbPrime * sizeof(ctx->tbPrime[0]));
    if(ctx->tbPrime == NULL) { return Free_tablePrime64(ctx) ; }
    //   ctx->maxValue = 0x7fffffff ;
    ctx->maxValue = ctx->maxNbPrime * (4 + log (ctx->maxNbPrime))  ;
    FindPrime64(ctx->maxValue,ctx,CPL_tablePrime64) ;
    return ctx ;
}

static int32_t CmpPrime(const void *e1,const void *e2) {
    T_prime64 * p1 = (T_prime64 *)e1 ;
    T_prime64 * p2 = (T_prime64 *)e2 ;
    if(*p1 > *p2) return 1;
    else if( *p1 < *p2) return -1 ;
    return 0 ;
}

int Search_TablePrime64(CTX_PRIMETABLE64 *ctxP, T_prime64 n) {
    return (bsearch(&n,ctxP->tbPrime,ctxP->nbPrime,sizeof(n),CmpPrime) != NULL) ;
}

int SearchRg_TablePrime64(CTX_PRIMETABLE64 *ctxP, T_prime64 n) {
    T_prime64 * pt= bsearch(&n,ctxP->tbPrime,ctxP->nbPrime,sizeof(n),CmpPrime) ;
    if(pt != NULL) return (int)(pt-ctxP->tbPrime) ;
    else return -1 ;
}
typedef struct PB708_N {
    int64_t val ;
    int8_t nbp2 ;
    
} PB708_N ;
int cmpPow2537(const void *el1,const void *el2) {
    PB708_N *n1 = (PB708_N *)el1 ;
    PB708_N *n2 = (PB708_N *)el2 ;
    if(n1->val > n2->val) return 1;
    else if(n1->val < n2->val) return -1 ;
    return 0 ;
}

int PB708b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
   
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;

//    int32_t N2 = (int32_t) Sqrt64(N+1) ;
    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
    int32_t invN2 = (int32_t)(N / N2) ;

    int32_t NP = N2 ;
//    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
//    int32_t N2 = (int32_t) N+1 ;
    
    uint8_t * Nb2_l = calloc(10*NP+1,sizeof(Nb2_l[0]));
    uint8_t * Nb2_h = calloc(N/NP+1,sizeof(Nb2_h[0]));

    uint64_t * S2_l = calloc(10*NP+1,sizeof(S2_l[0]));
    uint64_t * S2_h = calloc(N/NP+1,sizeof(S2_h[0]));

    uint64_t * NBNoP_l = calloc(10*NP+1,sizeof(NBNoP_l[0]));
    uint64_t * NBNoP_h = calloc(N/NP+1,sizeof(NBNoP_h[0]));
 
//    int64_t * countSQ = calloc(N2+1,sizeof(countSQ[0]));

   
    NBNoP_l[1] = 1 ;
    for(int p=2;p< /* 10* */ NP;p++) {
        if(Nb2_l[p]) {
            histPow2[Nb2_l[p]]++ ;
            NBNoP_l[p] = NBNoP_l[p-1]+1 ;
            continue ;
        } else {
            NBNoP_l[p] = (p <= N2) ? p : NBNoP_l[p-1] ;
            histPow2[1]++ ;
        }
        for(int64_t pp=p;pp</* 10* */ NP;pp *=p) {
            for(int64_t npp=pp;npp< /*10* */ NP;npp+=pp) {
                Nb2_l[npp]++ ;
            }
        }
    }
    S2_l[1] = 1 ;
    for(int i=2;i</* 10* */ NP;i++) {
        S2_l[i] = S2_l[i-1]+ (1LL << Nb2_l[i]) ;
    }
    
#define S2(i)  ( ( (i) < NP) ? S2_l[i] : S2_h[N/(i)])
#define NoP(i)  ( ( (i) < NP) ? NBNoP_l[i] : NBNoP_h[N/(i)])
/*
    printf("\nNBNoP_l=");
    for(int i=1;i<=NP;i++) {
        printf("[%d]=%lld ",i,NBNoP_l[i]);
    }
*/
    printf("\nN2=%d S2_l[%d]=%lld(%lld) \n",N2,NP-1,S2_l[NP-1],NBNoP_l[NP-1]);
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int64_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    for(int64_t in = N/(NP-1);in>0;in--) {
        int64_t Ni0 = N / in ;
        int64_t S2in = 1 ;
        int64_t NBnoPin = 1 ;
        int64_t p1 ;
        for(int ip1=0;ip1 < nbPrime && (p1=tbPrime[ip1]) <= Ni0  ; ip1++) {
            int64_t Ni1 = Ni0/p1 ;
            S2in += S2(Ni1) * 2 ;
            NBnoPin += NoP(Ni1) ;
 //           if(Ni1<N2)    countSQ[Ni1]++ ;
 //           else countSQ[N2]++ ;

            int64_t p2 ;
            for(int ip2=ip1+1; ip2 < nbPrime && (p2=tbPrime[ip2]) <=  Ni1; ip2++ ) {
                int64_t Ni2 = Ni1/p2 ;
                S2in -= S2(Ni2) * 4 ;
                NBnoPin -= NoP(Ni2) ;
 //               if(Ni2<N2)    countSQ[Ni2]++ ;
 //               else countSQ[N2]++ ;
                int64_t p3 ;
                for(int ip3=ip2+1; ip3 < nbPrime && (p3=tbPrime[ip3]) <=  Ni2; ip3++ ) {
                    int64_t Ni3 = Ni2/p3 ;
                    S2in += S2(Ni3) * 8 ;
                    NBnoPin += NoP(Ni3) ;
  //                  if(Ni3<N2)    countSQ[Ni3]++ ;
  //                  else countSQ[N2]++ ;
                    int64_t p4 ;
                    for(int ip4=ip3+1; ip4 < nbPrime && (p4=tbPrime[ip4]) <=  Ni3; ip4++ ) {
                        int64_t Ni4 = Ni3/p4 ;
                        S2in -= S2(Ni4) * 16 ;
                        NBnoPin -= NoP(Ni4) ;
  //                      if(Ni4<N2)    countSQ[Ni4]++ ;
  //                      else countSQ[N2]++ ;
                        int64_t p5 ;
                        for(int ip5=ip4+1; ip5 < nbPrime && (p5=tbPrime[ip5]) <=  Ni4; ip5++ ) {
                            int64_t Ni5 = Ni4/p5 ;
                            S2in += S2(Ni5) * 32 ;
                            NBnoPin += NoP(Ni5) ;
  //                          if(Ni5<N2)    countSQ[Ni5]++ ;
  //                          else countSQ[N2]++ ;
                           int64_t p6 ;
                            for(int ip6=ip5+1; ip6 < nbPrime && (p6=tbPrime[ip6]) <=  Ni5; ip6++ ) {
                                int64_t Ni6 = Ni5/p6 ;
                                S2in -= S2(Ni6) * 64 ;
                                NBnoPin -= NoP(Ni6) ;
  //                              if(Ni6<N2)    countSQ[Ni6]++ ;
  //                              else countSQ[N2]++ ;
                                int64_t p7 ;
                                for(int ip7=ip6+1; ip7 < nbPrime && (p7=tbPrime[ip7]) <=  Ni6; ip7++ ) {
                                    int64_t Ni7 = Ni6/p7 ;
                                    S2in += S2(Ni7) * 128 ;
                                    NBnoPin += NoP(Ni7) ;
  //                                  if(Ni7<N2)    countSQ[Ni7]++ ;
  //                                  else countSQ[N2]++ ;
                                    int64_t p8 ;
                                    for(int ip8=ip7+1; ip8 < nbPrime && (p8=tbPrime[ip8]) <=  Ni7; ip8++ ) {
                                        int64_t Ni8 = Ni7/p8 ;
                                        S2in -= S2(Ni8) * 256 ;
                                        NBnoPin -= NoP(Ni8) ;
  //                                      if(Ni8<N2)    countSQ[Ni8]++ ;
  //                                      else countSQ[N2]++ ;
                                   }
                                }
                           }
                        }
                    }
                }
            }
        }
        NBNoP_h[in] = NBnoPin ;
 //       S2_h[in] = S2in + 2*(Ni0 - NBnoPin);
        S2_h[in] = S2in  ;

//        if(Ni0 < N2*10) printf("(in=%lld)[%lld]=%lld/%lld(NoP=%lld)) ",in,Ni0,S2_h[in],S2_l[Ni0],NBnoPin);
    }
/*
    printf("\nS2_h=");
    for(int i=1;i<=NP;i++) {
        printf("[%d]=%lld ",i,S2_h[i]);
    }
 */
    printf("\n");
    int64_t nbN=NBNoP_h[1] ;
    for(int i=1;i<N2;i++) {
        NBNoP_h[i]= N/i - NBNoP_h[i]; // donc countSQ[n] contient les nb premiers p >N2 et  n *p < N
//        printf("(%d)%lld->%lld ",i,N/i - NBNoP_h[i],NBNoP_h[i]);
    }
    memset(histPow2,0,maxPow2*sizeof(histPow2[0]));
    for(int n=N2-1;n>1;n--) {
        //       printf("H[%d]+=(%lld/%d)%lld ",Nb2[n]+1,N,n,countSQ[n]);
        histPow2[ Nb2_l[n]+1] += NBNoP_h[n] ; // on rajoute les p * n
        nbN += NBNoP_h[n] ;
        int d ;
        int64_t pnc = NBNoP_h[n] ;
        for(d=2;d*d<n && d < N2 ;d++) {
            // on retranche si d diviseur strict de n , les p * n a countSQ[d]
            if((n % d) ==0) {
                NBNoP_h[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 < N2) NBNoP_h[d2] -=pnc ;
            }
        }
        if(d*d == n && d <N2 ) NBNoP_h[d] -=pnc ;
    }
    histPow2[1] += N-nbN ;

    
 
    
    int64_t S = S2_h[1] ;
    for(int i2=1;i2<maxPow2;i2++) {
 //       printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
  
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


int PB708b1(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    
    PB708_N *pow2357 = malloc(4000000*sizeof(pow2357[0])) ;
    int nbPow2357 = 0 ;
    pow2357[nbPow2357].val = 1 ;
    pow2357[nbPow2357++].nbp2 = 0 ;
    int n2=1 ;
    for(int64_t pow2=2;pow2<=N;pow2 *= 2 , n2++) {
        pow2357[nbPow2357].val = pow2 ;
        pow2357[nbPow2357++].nbp2 = n2 ;
    }
    
 
    
#define PB708_NBP   11
    int prim[PB708_NBP]= {3,5,7,11,13,17,19,23,29,31,37} ;
//#define PB708_NBP   3
//    int prim[PB708_NBP]= {3,5,7} ;
    for(int ip=0;ip<PB708_NBP;ip++) {
        int p = prim[ip];
        int nbp = 1 ;
        for(int64_t powp=p;powp<=N;powp *= p , nbp++) {
            int64_t pow2_p ;
            for(int ip2=0;(pow2_p=pow2357[ip2].val*powp)<=N;ip2++) {
                pow2357[nbPow2357].val = pow2_p ;
                pow2357[nbPow2357++].nbp2 = pow2357[ip2].nbp2 + nbp ;
            }
        }
        qsort(pow2357,nbPow2357,sizeof(pow2357[0]),cmpPow2537);
    }
    printf("Maxp2=%d nbp=%d\n",maxPow2,nbPow2357);
//    for(int i=0;i<nbPow2357;i++) printf("%lld(%d) ",pow2357[i].val,pow2357[i].nbp2);
    int64_t nbN = nbPow2357 ;
//    int32_t N2 = (int32_t) Sqrt64(N+1) ;
   int32_t N2 = (int32_t) Sqrt64(N)+1 ;
    int32_t invN2 = (int32_t)(N / N2) ;
    
    uint8_t * Nb2 = calloc(N2+1,sizeof(Nb2[0]));
    for(int p=2;p<N2;p++) {
        if(Nb2[p]) continue ;
        for(int64_t pp=p;pp<N2;pp *=p) {
            for(int64_t npp=pp;npp<N2;npp+=pp) {
                Nb2[npp]++ ;
            }
        }
    }
    
    int64_t * countSQ = calloc(N2+1,sizeof(countSQ[0]));
    int8_t * nb2SQ = calloc(N2,sizeof(nb2SQ[0])) ;
    for(int i=1;i<nbPow2357;i++)  {
        int64_t pn0 = pow2357[i].val ;
 /*       if(pn0 <N2 ) {
            nb2SQ[pn0]=pow2357[i].nbp2 ;
            countSQ[pn0] ++ ;
        }
       int d ;
        for(d=2;d*d<pn0 && d < N2 ;d++) {
//            printf("(%lld,%d)",pn0,d);
            if((pn0 % d) ==0) {
                countSQ[d] ++ ;
                int64_t d2 = pn0/d ;
                if(d2 < N2) countSQ[d2] ++ ;
            }
 
        }
        if(d*d == pn0 && d <N2 ) countSQ[d] ++ ;
  */
        histPow2[pow2357[i].nbp2]++;
    }
/*
    printf("\nP0 ");
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
    }
    printf("\n");
*/
//    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N+1) ;
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int64_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    int64_t nbPrimeSq ;
    printf("%.3fs N=%lld Nb primes=%lld N2=%d\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC,N, nbPrime,N2);
//    for(int n=2;n<N2;n++) {
     {
        int64_t count_n = 0 ;
        for(int i0=0;i0<nbPow2357;i0++) {
            int64_t p0 = pow2357[i0].val ;
              if(p0 > N) break ;
//            printf("%d+%lld ",n,p0) ;
//            count_n++ ;
            int nb0 = pow2357[i0].nbp2 ;
            if(2*p0<=N) {
                if(p0>=invN2) {
//                    printf("%lld->%lld ",p0,N/p0) ;
                    countSQ[N/p0]++ ;
                }
                else countSQ[N2]++ ;
            }
            int ip1 ;
            for(ip1=PB708_NBP+1;ip1<nbPrime;ip1++) {
                T_prime64 p1 = tbPrime[ip1] ;
                int64_t pn1 = p1*p0 ;
                if(pn1 > N) break ;
   //               printf("%d+%lld ",n,p1) ;
 //              count_n++ ;
                histPow2[nb0+1]++ ; nbN++ ;
 //                   printf("{%lld+%d}",pn1,nb0+1) ;
                if(2*pn1<=N) {
                    if(pn1 >=invN2)     {
 //                       printf("%lld->%lld ",pn1,N/pn1) ;
                        countSQ[N/pn1]++ ;
                        
                    }
                    else countSQ[N2]++ ;
                }
               for(int ip2=ip1;ip2<nbPrime;ip2++) {
                    T_prime64 p2 = tbPrime[ip2] ;
                    int64_t pn2 = pn1*p2 ;
                     if(pn2 > N) {  break ; }
 //                   count_n++ ;
                    histPow2[nb0+2]++ ; nbN++ ;
  //                      printf("{%lld+%d}",pn2,nb0+2) ;
                    if(2*pn2<=N) {
                       if(pn2 >=invN2 ) countSQ[N/pn2]++ ;
                       else countSQ[N2]++ ;
                   }

                
 //                    for(int ip3=PB708_NBP+1;ip3<=ip2;ip3++) {
                    for(int ip3=ip2;ip3<nbPrime;ip3++) {
                       T_prime64 p3 = tbPrime[ip3] ;
                        int64_t pn3 = pn2*p3 ;
                        if(pn3 > N) break ;
 //                       count_n++ ;
                        histPow2[nb0+3]++ ; nbN++ ;
                        if(2*pn3<=N) {
                            if(pn3 >=invN2) countSQ[N/pn3]++ ;
                            else countSQ[N2]++ ;
                        }

 //                       for(int ip4=PB708_NBP+1;ip4<=ip3;ip4++) {
                        for(int ip4=ip3;ip4<nbPrime;ip4++) {
                            T_prime64 p4 = tbPrime[ip4] ;
                            int64_t pn4 = pn3*p4 ;
                            if(pn4 > N) break ;
  //                          count_n++ ;
                              histPow2[nb0+4]++ ; nbN++ ;
                             if(2*pn4<=N) {
                                 if(pn4 >=invN2) countSQ[N/pn4]++ ;
                                 else countSQ[N2]++ ;
                             }

                            

 //                           for(int ip5=PB708_NBP+1;ip5<=ip4;ip5++) {
                            for(int ip5=ip4;ip5<nbPrime;ip5++) {
                               T_prime64 p5 = tbPrime[ip5] ;
                                int64_t pn5 = pn4*p5 ;
                                if(pn5 > N) break ;
  //                              count_n++ ;
                                histPow2[nb0+5]++ ; nbN++ ;
                                if(pn5 >=invN2) countSQ[N/pn5]++ ;
                                else countSQ[N2]++ ;
//                               for(int ip6=PB708_NBP+1;ip6<=ip5;ip6++) {
                                for(int ip6=ip5;ip6<nbPrime;ip6++) {
                                   T_prime64 p6 = tbPrime[ip6] ;
                                    int64_t pn6 = pn5*p6 ;
                                    if(pn6 > N) break ;
   //                                 count_n++ ;
                                    histPow2[nb0+6]++ ; nbN++ ;
                                    if(pn6 >=invN2) countSQ[N/pn6]++ ;
                                    else countSQ[N2]++ ;
                               }
                            }
                        }
                    }
               }
            }
       }
    }
    for(int n=N2-1;n>=2;n--) {
        // on cumule puisque si un p1*...*pk * n <N c'est vrai pour n plus petit.
         countSQ[n] +=  countSQ[n+1];
//        printf("(%d)=%lld ",n,countSQ[n]);
    }
    printf("N=%lld nbB=%lld \n",N,nbN);
//    int64_t S = 1 + (N-nbN)*2;
    int64_t S = 1 ;
    printf("\n Count=");
    for(int i=1;i<N2;i++) {
        countSQ[i]= N/i - countSQ[i]; // donc countSQ[n] contient les nb premiers p >N2 et  n *p < N
 //       printf("(%d)%lld->%lld ",i,N/i - countSQ[i],countSQ[i]);
    }
 /*
    printf("\nP1 ");
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
    }
    printf("\n");
*/
    for(int n=N2-1;n>1;n--) {
 //       printf("H[%d]+=(%lld/%d)%lld ",Nb2[n]+1,N,n,countSQ[n]);
        histPow2[ Nb2[n]+1] += countSQ[n] ; // on rajoute les p * n
        nbN += countSQ[n] ;
        int d ;
        int64_t pnc = countSQ[n] ;
        for(d=2;d*d<n && d < N2 ;d++) {
            // on retranche si d diviseur strict de n , les p * n a countSQ[d]
            if((n % d) ==0) {
                countSQ[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 < N2) countSQ[d2] -=pnc ;
            }
        }
        if(d*d == n && d <N2 ) countSQ[d] -=pnc ;
   }
        histPow2[1] += N-nbN ;
//    histPow2[2] += (N/2-nbN2) ;
//    histPow2[1] += (N-nbN) - (N/2-nbN2) ;
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

int PB708c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    
    PB708_N *pow2357 = malloc(4000000*sizeof(pow2357[0])) ;
    int nbPow2357 = 0 ;
    pow2357[nbPow2357].val = 1 ;
    pow2357[nbPow2357++].nbp2 = 0 ;
    int n2=1 ;
    for(int64_t pow2=2;pow2<=N;pow2 *= 2 , n2++) {
        pow2357[nbPow2357].val = pow2 ;
        pow2357[nbPow2357++].nbp2 = n2 ;
    }
    
    
    
#define PB708_NBP   11
    int prim[PB708_NBP]= {3,5,7,11,13,17,19,23,29,31,37} ;
    //#define PB708_NBP   3
    //    int prim[PB708_NBP]= {3,5,7} ;
    for(int ip=0;ip<PB708_NBP;ip++) {
        int p = prim[ip];
        int nbp = 1 ;
        for(int64_t powp=p;powp<=N;powp *= p , nbp++) {
            int64_t pow2_p ;
            for(int ip2=0;(pow2_p=pow2357[ip2].val*powp)<=N;ip2++) {
                pow2357[nbPow2357].val = pow2_p ;
                pow2357[nbPow2357++].nbp2 = pow2357[ip2].nbp2 + nbp ;
            }
        }
        qsort(pow2357,nbPow2357,sizeof(pow2357[0]),cmpPow2537);
    }
    printf("Maxp2=%d nbp=%d\n",maxPow2,nbPow2357);
    //    for(int i=0;i<nbPow2357;i++) printf("%lld(%d) ",pow2357[i].val,pow2357[i].nbp2);
    int64_t nbN = nbPow2357 ;
    //    int32_t N2 = (int32_t) Sqrt64(N+1) ;
    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
    
    
    uint8_t * Nb2 = calloc(N2+1,sizeof(Nb2[0]));
    for(int p=2;p<N2;p++) {
        if(Nb2[p]) continue ;
        for(int64_t pp=p;pp<N2;pp *=p) {
            for(int64_t npp=pp;npp<N2;npp+=pp) {
                Nb2[npp]++ ;
            }
        }
    }
    
    int64_t * countSQ = calloc(N2,sizeof(countSQ[0]));
    int8_t * nb2SQ = calloc(N2,sizeof(nb2SQ[0])) ;
    for(int i=1;i<nbPow2357;i++)  {
        int64_t pn0 = pow2357[i].val ;
        /*       if(pn0 <N2 ) {
         nb2SQ[pn0]=pow2357[i].nbp2 ;
         countSQ[pn0] ++ ;
         }
         int d ;
         for(d=2;d*d<pn0 && d < N2 ;d++) {
         //            printf("(%lld,%d)",pn0,d);
         if((pn0 % d) ==0) {
         countSQ[d] ++ ;
         int64_t d2 = pn0/d ;
         if(d2 < N2) countSQ[d2] ++ ;
         }
         
         }
         if(d*d == pn0 && d <N2 ) countSQ[d] ++ ;
         */
        histPow2[pow2357[i].nbp2]++;
    }
    /*
     printf("\nP0 ");
     for(int i2=1;i2<maxPow2;i2++) {
     printf("h[%d]=%lld ",i2,histPow2[i2]);
     }
     printf("\n");
     */
    //    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N+1) ;
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int64_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    int64_t nbPrimeSq ;
    printf("Nb primes=%lld N2=%lld\n",nbPrime,N2);
    //    for(int n=2;n<N2;n++) {
    for(int n=1;n<N2;n++) {
        for(int i0=0;i0<nbPow2357;i0++) {
            int64_t p0 = pow2357[i0].val*n ;
            if(p0 > N) break ;
            //            printf("%d+%lld ",n,p0) ;
            countSQ[n]++ ;
            int nb0 = pow2357[i0].nbp2+Nb2[n] ;
            for(int ip1=PB708_NBP+1;ip1<nbPrime;ip1++) {
                T_prime64 p1 = tbPrime[ip1] ;
                int64_t pn1 = p1*p0 ;
                if(pn1 > N) break ;
                //               printf("%d+%lld ",n,p1) ;
                countSQ[n]++ ;
                if(n==1)  {histPow2[nb0+1]++ ; nbN++ ;
                    //                   printf("{%lld+%d}",pn1,nb0+1) ;
                }
                for(int ip2=PB708_NBP+1;ip2<=ip1;ip2++) {
                    T_prime64 p2 = tbPrime[ip2] ;
                    int64_t pn2 = pn1*p2 ;
                    if(pn2 > N) break ;
                    countSQ[n]++ ;
                    if(n==1)  { histPow2[nb0+2]++ ; nbN++ ;
                        //                      printf("{%lld+%d}",pn2,nb0+2) ;
                    }
                    for(int ip3=PB708_NBP+1;ip3<=ip2;ip3++) {
                        T_prime64 p3 = tbPrime[ip3] ;
                        int64_t pn3 = pn2*p3 ;
                        if(pn3 > N) break ;
                        countSQ[n]++ ;
                        if(n==1)  {histPow2[nb0+3]++ ; nbN++ ;}
                        for(int ip4=PB708_NBP+1;ip4<=ip3;ip4++) {
                            T_prime64 p4 = tbPrime[ip4] ;
                            int64_t pn4 = pn3*p4 ;
                            if(pn4 > N) break ;
                            countSQ[n]++ ;
                            if(n==1)  {histPow2[nb0+4]++ ; nbN++ ;}
                            for(int ip5=PB708_NBP+1;ip5<=ip4;ip5++) {
                                T_prime64 p5 = tbPrime[ip5] ;
                                int64_t pn5 = pn4*p5 ;
                                if(pn5 > N) break ;
                                countSQ[n]++ ;
                                if(n==1)  {histPow2[nb0+5]++ ; nbN++ ;}
                                for(int ip6=PB708_NBP+1;ip6<=ip5;ip6++) {
                                    T_prime64 p6 = tbPrime[ip6] ;
                                    int64_t pn6 = pn5*p6 ;
                                    if(pn6 > N) break ;
                                    countSQ[n]++ ;
                                    if(n==1)  {histPow2[nb0+6]++ ; nbN++ ;}
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("N=%lld nbB=%lld \n",N,nbN);
    //    int64_t S = 1 + (N-nbN)*2;
    int64_t S = 1 ;
    printf("\n Count=");
    for(int i=1;i<N2;i++) {
        countSQ[i]= N/i - countSQ[i];
 //       printf("(%d)%lld->%lld ",i,N/i - countSQ[i],countSQ[i]);
    }
    /*
     printf("\nP1 ");
     for(int i2=1;i2<maxPow2;i2++) {
     printf("h[%d]=%lld ",i2,histPow2[i2]);
     }
     printf("\n");
     */
    for(int n=N2-1;n>1;n--) {
        //       printf("H[%d]+=(%lld/%d)%lld ",Nb2[n]+1,N,n,countSQ[n]);
        histPow2[ Nb2[n]+1] += countSQ[n] ;
        nbN += countSQ[n] ;
        int d ;
        int64_t pnc = countSQ[n] ;
        for(d=2;d*d<n && d < N2 ;d++) {
            if((n % d) ==0) {
                countSQ[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 < N2) countSQ[d2] -=pnc ;
            }
        }
        if(d*d == n && d <N2 ) countSQ[d] -=pnc ;
    }
    histPow2[1] += N-nbN ;
    //    histPow2[2] += (N/2-nbN2) ;
    //    histPow2[1] += (N-nbN) - (N/2-nbN2) ;
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


int PB708c1(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    
    PB708_N *pow2357 = malloc(4000000*sizeof(pow2357[0])) ;
    int nbPow2357 = 0 ;
    pow2357[nbPow2357].val = 1 ;
    pow2357[nbPow2357++].nbp2 = 0 ;
    int n2=1 ;
    for(int64_t pow2=2;pow2<=N;pow2 *= 2 , n2++) {
        pow2357[nbPow2357].val = pow2 ;
        pow2357[nbPow2357++].nbp2 = n2 ;
    }
//#define PB708_NBP   11
//    int prim[PB708_NBP]= {3,5,7,11,13,17,19,23,29,31,37} ;
    #define PB708_NBP   3
        int prim[PB708_NBP]= {3,5,7} ;
    for(int ip=0;ip<PB708_NBP;ip++) {
        int p = prim[ip];
        int nbp = 1 ;
        for(int64_t powp=p;powp<=N;powp *= p , nbp++) {
            int64_t pow2_p ;
            for(int ip2=0;(pow2_p=pow2357[ip2].val*powp)<=N;ip2++) {
                pow2357[nbPow2357].val = pow2_p ;
                pow2357[nbPow2357++].nbp2 = pow2357[ip2].nbp2 + nbp ;
            }
        }
        qsort(pow2357,nbPow2357,sizeof(pow2357[0]),cmpPow2537);
    }
    printf("Maxp2=%d nbp=%d\n",maxPow2,nbPow2357);
    //    for(int i=0;i<nbPow2357;i++) printf("%lld(%d) ",pow2357[i].val,pow2357[i].nbp2);
    int64_t nbN = nbPow2357 ;
    int32_t N2 = (int32_t) Sqrt64(N+1) ;
    int64_t * countSQ = calloc(N2,sizeof(countSQ[0]));
    int8_t * nb2SQ = calloc(N2,sizeof(nb2SQ[0])) ;
    for(int i=1;i<nbPow2357;i++)  {
        int64_t pn0 = pow2357[i].val ;
        if(pn0 <N2 ) {
            nb2SQ[pn0]=pow2357[i].nbp2 ;
            countSQ[pn0] ++ ;
        }
        int d ;
        for(d=2;d*d<pn0 && d < N2 ;d++) {
            //            printf("(%lld,%d)",pn0,d);
            if((pn0 % d) ==0) {
                countSQ[d] ++ ;
                int64_t d2 = pn0/d ;
                if(d2 < N2) countSQ[d2] ++ ;
            }
            
        }
        if(d*d == pn0 && d <N2 ) countSQ[d] ++ ;
        histPow2[pow2357[i].nbp2]++;
    }
/*    printf("\nP0 ");
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
    }
    printf("\n");
*/
    //    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N+1) ;
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int64_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    int64_t nbPrimeSq ;
    int64_t p0 ;
    printf("Nb primes=%lld \n",nbPrime);
    for(int ip1=PB708_NBP+1;ip1<nbPrime;ip1++) {
        T_prime64 p1 = tbPrime[ip1] ;
        //        if(2*p1>N) break ;
        int64_t n1 = p1 ;
        for(int i=0;i<nbPow2357;i++) {
            int64_t pn1 = pow2357[i].val*n1 ;
            if(pn1 > N) break ;
            if(pn1 <N2 ) {
                nb2SQ[pn1]=pow2357[i].nbp2+1 ;
                countSQ[pn1] ++ ;
            }
            int d ;
            for(d=2;d*d<pn1 && d < N2 ;d++) {
                if((pn1 % d) ==0) {
                    countSQ[d] ++ ;
                    int64_t d2 = pn1/d ;
                    if(d2 < N2) countSQ[d2] ++ ;
                }
            }
            if(d*d == pn1 && d <N2 ) countSQ[d] ++ ;
            histPow2[pow2357[i].nbp2+1]++ ; nbN++ ;
  //          printf("{%lld+%d}",pn1,pow2357[i].nbp2+1) ;
        }
        for(int ip2=PB708_NBP+1;ip2<=ip1;ip2++) {
            T_prime64 p2 = tbPrime[ip2] ;
            int64_t n2 = n1*p2 ;
            if(n2 > N) break ;
            for(int i=0;i<nbPow2357;i++) {
                int64_t pn2 = pow2357[i].val*n2 ;
                if(pn2 > N) break ;
                if(pn2 <N2 ) {
                    nb2SQ[pn2]=pow2357[i].nbp2+2 ;
                    countSQ[pn2] ++ ;
                }
                int d ;
                for(d=2;d*d<pn2 && d < N2 ;d++) {
                    if((pn2 % d) ==0) {
                        countSQ[d] ++ ;
                        int64_t d2 = pn2/d ;
                        if(d2 < N2) countSQ[d2] ++ ;
                    }
                }
                if(d*d == pn2 && d <N2 ) countSQ[d] ++ ;
                histPow2[pow2357[i].nbp2+2]++ ; nbN++ ;
 //               printf("{%lld+%d}",pn2,pow2357[i].nbp2+2) ;

            }
            for(int ip3=PB708_NBP+1;ip3<=ip2;ip3++) {
                T_prime64 p3 = tbPrime[ip3] ;
                int64_t n3 = n2*p3 ;
                if(n3 > N) break ;
                for(int i=0;i<nbPow2357;i++) {
                    int64_t pn3 = pow2357[i].val*n3 ;
                    if(pn3 > N) break ;
                    if(pn3 <N2 ) {
                        nb2SQ[pn3]=pow2357[i].nbp2+3 ;
                        countSQ[pn3] ++ ;
                        
                    }
                    int d ;
                    for(d=2;d*d<pn3 && d < N2 ;d++) {
                        if((pn3 % d) ==0) {
                            countSQ[d] ++ ;
                            int64_t d2 = pn3/d ;
                            if(d2 < N2) countSQ[d2] ++ ;
                        }
                    }
                    if(d*d == pn3 && d <N2 ) countSQ[d] ++ ;
                    histPow2[pow2357[i].nbp2+3]++ ; nbN++ ;
                }
                for(int ip4=PB708_NBP+1;ip4<=ip3;ip4++) {
                    T_prime64 p4 = tbPrime[ip4] ;
                    int64_t n4 = n3*p4 ;
                    if(n4 > N) break ;
                    for(int i=0;i<nbPow2357;i++) {
                        int64_t pn4 = pow2357[i].val*n4 ;
                        if(pn4 > N) break ;
                        if(pn4 <N2 ) {
                            nb2SQ[pn4]=pow2357[i].nbp2+4 ;
                            countSQ[pn4] ++ ;
                            
                        }
                        int d ;
                        for(d=2;d*d<pn4 && d < N2 ;d++) {
                            if((pn4 % d) ==0) {
                                countSQ[d] ++ ;
                                int64_t d2 = pn4/d ;
                                if(d2 < N2) countSQ[d2] ++ ;
                            }
                        }
                        if(d*d == pn4 && d <N2 ) countSQ[d] ++ ;
                        histPow2[pow2357[i].nbp2+4]++ ; nbN++ ;
                    }
                    for(int ip5=PB708_NBP+1;ip5<=ip4;ip5++) {
                        T_prime64 p5 = tbPrime[ip5] ;
                        int64_t n5 = n4*p5 ;
                        if(n5 > N) break ;
                        for(int i=0;i<nbPow2357;i++) {
                            int64_t pn5 = pow2357[i].val*n5 ;
                            if(pn5 > N) break ;
                            if(pn5 <N2 ) {
                                nb2SQ[pn5]=pow2357[i].nbp2+5 ;
                                countSQ[pn5] ++ ;
                                
                            }
                            int d ;
                            for(d=2;d*d<pn5 && d < N2 ;d++) {
                                if((pn5 % d) ==0) {
                                    countSQ[d] ++ ;
                                    int64_t d2 = pn5/d ;
                                    if(d2 < N2) countSQ[d2] ++ ;
                                }
                            }
                            if(d*d == pn5 && d <N2 ) countSQ[d] ++ ;
                            histPow2[pow2357[i].nbp2+5]++ ; nbN++ ;
                        }
                        for(int ip6=PB708_NBP+1;ip6<=ip5;ip6++) {
                            T_prime64 p6 = tbPrime[ip6] ;
                            int64_t n6 = n5*p6 ;
                            if(n6 > N) break ;
                            for(int i=0;i<nbPow2357;i++) {
                                int64_t pn6 = pow2357[i].val*n6 ;
                                if(pn6 > N) break ;
                                if(pn6 <N2 ) {
                                    nb2SQ[pn6]=pow2357[i].nbp2+6 ;
                                    countSQ[pn6] ++ ;
                                    
                                }
                                int d ;
                                for(d=2;d*d<pn6 && d < N2 ;d++) {
                                    if((pn6 % d) ==0) {
                                        countSQ[d] ++ ;
                                        int64_t d2 = pn6/d ;
                                        if(d2 < N2) countSQ[d2] ++ ;
                                    }
                                }
                                if(d*d == pn6 && d <N2 ) countSQ[d] ++ ;
                                histPow2[pow2357[i].nbp2+6]++ ; nbN++ ;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("N=%lld nbB=%lld \n",N,nbN);
    //    int64_t S = 1 + (N-nbN)*2;
    int64_t S = 1 ;
    printf("\n Count=");
    for(int i=1;i<N2;i++) {
        countSQ[i]= N/i - countSQ[i];
        printf("(%d)%lld->%lld ",i,N/i - countSQ[i],countSQ[i]);
    }
    
  /*
    printf("\nP1 ");
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
    }
    printf("\n");
*/
    for(int n=N2-1;n>1;n--) {
//        printf("H[%d]+=(%lld/%d)%lld ",nb2SQ[n]+1,N,n,countSQ[n]);
        histPow2[ nb2SQ[n]+1] += countSQ[n] ;
        nbN += countSQ[n] ;
        int d ;
        int64_t pnc = countSQ[n] ;
        for(d=2;d*d<n && d < N2 ;d++) {
            if((n % d) ==0) {
                countSQ[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 < N2) countSQ[d2] -=pnc ;
            }
        }
        if(d*d == n && d <N2 ) countSQ[d] -=pnc ;
    }
    histPow2[1] += N-nbN ;
    //    histPow2[2] += (N/2-nbN2) ;
    //    histPow2[1] += (N-nbN) - (N/2-nbN2) ;
 
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


int PB708a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    
    PB708_N *pow2357 = malloc(400000*sizeof(pow2357[0])) ;
    int nbPow2357 = 0 ;
    pow2357[nbPow2357].val = 1 ;
    pow2357[nbPow2357++].nbp2 = 0 ;
    int n2=1 ;
    for(int64_t pow2=2;pow2<=N;pow2 *= 2 , n2++) {
        pow2357[nbPow2357].val = pow2 ;
        pow2357[nbPow2357++].nbp2 = n2 ;
    }
#define PB708_NBP   11
    int prim[PB708_NBP]= {3,5,7,11,13,17,19,23,29,31,37} ;
    for(int ip=0;ip<PB708_NBP;ip++) {
        int p = prim[ip];
        int nbp = 1 ;
        for(int64_t powp=p;powp<=N;powp *= p , nbp++) {
            int64_t pow2_p ;
            for(int ip2=0;(pow2_p=pow2357[ip2].val*powp)<=N;ip2++) {
                pow2357[nbPow2357].val = pow2_p ;
                pow2357[nbPow2357++].nbp2 = pow2357[ip2].nbp2 + nbp ;
            }
        }
        qsort(pow2357,nbPow2357,sizeof(pow2357[0]),cmpPow2537);
    }
    printf("Maxp2=%d nbp=%d\n",maxPow2,nbPow2357);
//    for(int i=0;i<nbPow2357;i++) printf("%lld(%d) ",pow2357[i].val,pow2357[i].nbp2);
    for(int i=1;i<nbPow2357;i++) histPow2[pow2357[i].nbp2]++;
    int64_t nbN = nbPow2357 ;
    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
    
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N+1) ;
    int64_t nbPrime = GetNbPrime64(ctxP) ;
    printf("Nb primes=%lld\n",nbPrime);
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    for(int ip1=PB708_NBP+1;ip1<nbPrime;ip1++) {
        T_prime64 p1 = tbPrime[ip1] ;
        int64_t n1 = p1 ;
        for(int i=0;i<nbPow2357;i++) {
            if(pow2357[i].val*n1 > N) break ;
 //           printf("%lld<%d> ",pow2357[i].val*n1,pow2357[i].nbp2+1) ;
            histPow2[pow2357[i].nbp2+1]++ ; nbN++ ;
        }
        for(int ip2=ip1;ip2<nbPrime;ip2++) {
            T_prime64 p2 = tbPrime[ip2] ;
            int64_t n2 = n1*p2 ;
            if(n2 > N) break ;
            for(int i=0;i<nbPow2357;i++) {
                if(pow2357[i].val*n2 > N) break ;
 //               printf("%lld[%d] ",pow2357[i].val*n2,pow2357[i].nbp2+2) ;
                histPow2[pow2357[i].nbp2+2]++ ; nbN++ ;
            }
            for(int ip3=ip2;ip3<nbPrime;ip3++) {
                T_prime64 p3 = tbPrime[ip3] ;
                int64_t n3 = n2*p3 ;
                if(n3 > N) break ;
                for(int i=0;i<nbPow2357;i++) {
                    if(pow2357[i].val*n3 > N) break ;
 //                   printf("%lld{%d} ",pow2357[i].val*n3,pow2357[i].nbp2+3) ;
                    histPow2[pow2357[i].nbp2+3]++ ; nbN++ ;
                }
                for(int ip4=ip3;ip4<nbPrime;ip4++) {
                    T_prime64 p4 = tbPrime[ip4] ;
                    int64_t n4 = n3*p4 ;
                    if(n4 > N) break ;
                    for(int i=0;i<nbPow2357;i++) {
                        if(pow2357[i].val*n4 > N) break ;
                        //                   printf("%lld{%d} ",pow2357[i].val*n3,pow2357[i].nbp2+3) ;
                        histPow2[pow2357[i].nbp2+4]++ ; nbN++ ;
                    }
                    for(int ip5=ip4;ip5<nbPrime;ip5++) {
                        T_prime64 p5 = tbPrime[ip5] ;
                        int64_t n5 = n4*p5 ;
                        if(n5 > N) break ;
                        for(int i=0;i<nbPow2357;i++) {
                            if(pow2357[i].val*n5 > N) break ;
                            //                   printf("%lld{%d} ",pow2357[i].val*n3,pow2357[i].nbp2+3) ;
                            histPow2[pow2357[i].nbp2+5]++ ; nbN++ ;
                        }
                        for(int ip6=ip5;ip6<nbPrime;ip6++) {
                            T_prime64 p6 = tbPrime[ip6] ;
                            int64_t n6 = n5*p6 ;
                            if(n5 > N) break ;
                            for(int i=0;i<nbPow2357;i++) {
                                if(pow2357[i].val*n6 > N) break ;
                                //                   printf("%lld{%d} ",pow2357[i].val*n3,pow2357[i].nbp2+3) ;
                                histPow2[pow2357[i].nbp2+6]++ ; nbN++ ;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("N=%lld nbB=%lld \n",N,nbN);
    int64_t S = 1 + (N-nbN)*2;
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
