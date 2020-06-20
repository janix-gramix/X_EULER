
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
        printf("h[%d]=%lld ",i2,histPow2[i2]);

        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}




typedef uint64_t T_prime64 ;
typedef int(*TY_CPL_nxtPrime64)(void *ctx,T_prime64 nxtPrime);
uint32_t FindPrime64(T_prime64 maxValue,void *ctx,TY_CPL_nxtPrime64 nxtPrime) ;
// typedef int(*TY_CPL_nxtValExp)(void *ctx,T_prime64 nxtVal,uint8_t exp);
typedef struct CTX_PRIMETABLE64 CTX_PRIMETABLE64 ;

struct  CTX_PRIMETABLE64 {
    uint32_t       nbPrime ;
    T_prime64     maxValue ;
    uint32_t    maxNbPrime ;
    T_prime64     *tbPrime ;
}  ;

typedef int(*TY_CPL_nxtExp)(void *ctx,uint8_t exp);

typedef struct CTX_PRIMETABLE64exp CTX_PRIMETABLE64exp ;

struct  CTX_PRIMETABLE64exp {
    uint32_t       nbPrime ;
    T_prime64     maxValue ;
    uint32_t    maxNbPrime ;
    int32_t     sizeH ;
    int64_t     *histoP2 ;
//    T_prime64     *tbPrime ;
    
//    void *       contexte ;
//    TY_CPL_nxtExp nxtExp ;
}  ;
/*
int CPL_tablePrime64exp(void *ptCtx,T_prime64 nxtPrime,uint8_t exp) {
    CTX_PRIMETABLE64exp * ctx = (CTX_PRIMETABLE64exp *) ptCtx ;
    if(nxtPrime > ctx->maxValue) return 0 ;
    if(exp ==0) {
        ctx->tbPrime[ctx->nbPrime] = nxtPrime ;
        ctx->nbPrime++ ;
    }
//    printf("%lld(%d) ",nxtPrime,exp);
    if(ctx->nbPrime >= ctx->maxNbPrime) {
        ctx->maxValue = nxtPrime * nxtPrime ;
        return 0 ;
    }
    if(ctx->nxtExp) ctx->nxtExp(ctx->contexte,exp) ;
    return 1;
}
*/

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

CTX_PRIMETABLE64exp * Free_tablePrime64exp(CTX_PRIMETABLE64exp * ctx) {
    if(ctx != NULL) free(ctx->histoP2) ;
    free(ctx);
    return NULL ;
}


//
// On replie la table.
// la taille de la table peut être quelconque
// Plus rapide pour taille nSqrt ou 32368 si trop grande valeurs
//
uint32_t FindPrime64(T_prime64 nbMax,void *ctx,TY_CPL_nxtPrime64 nxtPrime) {
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

uint32_t FindPrime64exp(CTX_PRIMETABLE64exp *ctx) {
    uint32_t nSqrt = 1+ (int32_t)Sqrt64(ctx->maxValue ) ;
    int isEnd = 0;
    T_prime64 *tbPrime = malloc(nSqrt * sizeof(tbPrime[0])) ;
    int64_t *offSet = malloc(nSqrt * sizeof(offSet[0])) ;
    uint32_t sizeTable =  (nSqrt < 65536) ? nSqrt : 65536 ;
    int8_t *exp = calloc( sizeTable , sizeof(exp[0])) ;
    int64_t *decomp = malloc(sizeTable * sizeof(decomp[0]));
    uint32_t nbPrime = 0 ;
    T_prime64 lastPrime = 0 ;
    T_prime64 offSetTable = 0 ;
    T_prime64 nbPrimeSqrt = 0 ;
    
    tbPrime[nbPrime++] = lastPrime = 2 ;
     T_prime64 indPrime = 3 ;
//    for(int i =0;i<sizeTable;i++) decomp[i] = 1 ;
    
    T_prime64 icp ;
    for(int i =2 ;i<sizeTable;i++ ) {
        icp = offSetTable + i ;
        int8_t exp2 = 0 ;
        int64_t ipow2 ;
        for(ipow2=1;(ipow2 & icp)==0;ipow2 <<=1 ) exp2++;
        decomp[i] = ipow2  ;
        exp[i] = exp2 ;
    }

    
    while ( 1) {
        T_prime64 icp ;
        for(icp=indPrime - offSetTable ;icp < sizeTable; icp++ ) {
//            printf("(%lld#%d ",icp+offSetTable,exp[icp]);
            if(!exp[icp] ) {
                lastPrime = icp+offSetTable ;
                if(lastPrime < nSqrt) {
                    T_prime64 icp2 ;
                    tbPrime[nbPrimeSqrt] = lastPrime ;
 //                   printf("<%lld> ",lastPrime) ;
                    exp[icp] = 1 ;
                    for(icp2 = icp+lastPrime ; icp2 < sizeTable ; icp2 += lastPrime ) {
                        exp[icp2]++ ;
                        decomp[icp2] *= lastPrime ;
//                        printf("{%lld,%lld} ",icp2+offSetTable,decomp[icp2]);
                    }
                    offSet[nbPrimeSqrt++] = icp2 - sizeTable ;
                } else {
                    exp[icp]++ ;
                }
                nbPrime++ ;
            } else if(decomp[icp] != icp+offSetTable) {
//                printf("%lld->%lld ",icp+offSetTable,decomp[icp]);
                 exp[icp]++ ;
            }
 //           if(exp[icp]==1)printf("{%lld}",ctx->histoP2[exp[icp]]);
 //           printf("%lld(%d) ",icp+offSetTable,exp[icp]);
            ctx->histoP2[exp[icp]]++ ;
            
            if(icp+offSetTable>= ctx->maxValue ) {
                isEnd = 1 ;
                break ;
            }
        }
        offSetTable += sizeTable ;
        if(isEnd || offSetTable >= ctx->maxValue) break ;
        indPrime = offSetTable ;
        if ( offSetTable + sizeTable > ctx->maxValue) { sizeTable = (int32_t)(ctx->maxValue - offSetTable+1) ; }
        for(int i =0 ;i<sizeTable;i++ ) {
            icp = offSetTable + i ;
            int8_t exp2 = 0 ;
            int64_t ipow2 ;
            for(ipow2=1;(ipow2 & icp)==0;ipow2 <<=1 ) exp2++;
            decomp[i] = ipow2  ;
            exp[i] = exp2 ;
        }

        {
            int np ;
            for(np=0;np<nbPrimeSqrt;np++) {
                T_prime64 p = tbPrime[np] ;
                int64_t indPrime = offSet[np] ;
                while ( indPrime < sizeTable) {
                    exp[indPrime]++ ;
                    decomp[indPrime] *= p ;
                    indPrime += p ;
                }
                offSet[np] = indPrime - sizeTable ;
            }
        }
    }
    free(decomp);
    free(tbPrime);
    free(offSet) ;
    free(exp);
    return nbPrime ;
}
CTX_PRIMETABLE64exp * Gen_tablePrime64exp(T_prime64 maxValue) {
    CTX_PRIMETABLE64exp *ctx ;
    ctx = calloc(1,sizeof(ctx[0]));
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= maxValue  ) maxPow2++ ;
    ctx->histoP2 = calloc(maxPow2,sizeof(ctx->histoP2[0])) ;
    ctx->sizeH = maxPow2 ;
    if(ctx == NULL) return ctx ;
    ctx->maxValue = maxValue ;
    if(maxValue > 100) {
        ctx->maxNbPrime = (uint32_t) (1+ maxValue / (log((double)maxValue) - 4)) ;
    } else {
        ctx->maxNbPrime = 30 ;
    }
     FindPrime64exp(ctx) ;
    return ctx ;
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

    int32_t N2 = (int32_t) Sqrt64(N) ;
//    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
    int32_t invN2 = (int32_t)(N / N2) ;

    int32_t NP = N2 ;
//    int32_t N2 = (int32_t) Sqrt64(N)+1 ;
//    int32_t N2 = (int32_t) N+1 ;
    
    uint8_t * Nb2_l = calloc(NP+1,sizeof(Nb2_l[0]));
    uint8_t * Nb2_h = calloc(N/NP+1,sizeof(Nb2_h[0]));

    uint64_t * S2_l = calloc(NP+1,sizeof(S2_l[0]));
    uint64_t * S2_h = calloc(N/NP+1,sizeof(S2_h[0]));

    uint64_t * NBNoP_l = calloc(NP+1,sizeof(NBNoP_l[0]));
    uint64_t * NBNoP_h = calloc(N/NP+1,sizeof(NBNoP_h[0]));
 
//    int64_t * countSQ = calloc(N2+1,sizeof(countSQ[0]));


    
// NBNoP_l[i] compte les k <i  
   
    NBNoP_l[1] = 1 ;
    for(int p=2;p<= NP;p++) {
        if(Nb2_l[p]) {
            histPow2[Nb2_l[p]]++ ;
            NBNoP_l[p] = NBNoP_l[p-1]+1 ; //
            continue ;
        } else {
            // p premier
            NBNoP_l[p] = (p <= N2) ? p : NBNoP_l[p-1] ; //
            histPow2[1]++ ;
        }
        for(int64_t pp=p;pp<= NP;pp *=p) { // p, p**2, p**3
            for(int64_t npp=pp;npp<= NP;npp+=pp) {
                Nb2_l[npp]++ ; // on incremente
            }
        }
    }
    S2_l[1] = 1 ;
    for(int i=2;i<= NP;i++) { // calcul g(i) pour i<=NP
        S2_l[i] = S2_l[i-1]+ (1LL << Nb2_l[i]) ;
    }
    
    for(int i=1;i<=NP;i++) NBNoP_l[i] = i ;
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
    

     uint64_t * NB_lh = calloc(NP+1+N/NP+1,sizeof(NB_lh[0]));
    uint64_t * NB_l = NB_lh+NP+1 ;

/*

    
#define NB(i)  ( ( (i) <= NP) ? (i) : (NP+1+N/(i)) )

    NB_lh[1] = 1 ;
    int nbContinue = 0 ;
//    for(int64_t in = 1;in*NP<N;in++) {
    for(int64_t in = 1;in<=NP;in++) {
        if(NB_lh[in]==0) {
            nbContinue++ ;
            continue ;
        }
 //       printf("\n%lld->",in);
 //       for(int i=0;i<2*NP+1;i++) printf("%lld ",NB_lh[i]);
        int64_t nbIn = NB_lh[in] ;
        NB_lh[NP+2] += nbIn ;
        NB_lh[in] = 0 ;
        int64_t p1 ;
        for(int ip1=0;ip1 < nbPrime && (p1=tbPrime[ip1])*in <= N  ; ip1++) {
            int64_t Ni1 = in*p1 ;
            NB_lh[NB(Ni1)] += 2*nbIn ; // printf("p1=%lld,NB(%lld)+=%lld ",p1,NB(Ni1),2*nbIn);
            int64_t p2 ;
            for(int ip2=ip1+1; ip2 < nbPrime && (p2=tbPrime[ip2])*Ni1 <=  N; ip2++ ) {
                int64_t Ni2 = Ni1 * p2 ;
                NB_lh[NB(Ni2)] -= 4*nbIn ; // printf("p2=%lld,NB(%lld)-=%lld ",p2,NB(Ni2),4*nbIn);
                int64_t p3 ;
                for(int ip3=ip2+1; ip3 < nbPrime && (p3=tbPrime[ip3])*Ni2 <=  N; ip3++ ) {
                    int64_t Ni3 = Ni2*p3 ;
                    NB_lh[NB(Ni3)] += 8*nbIn ;
                    int64_t p4 ;
                    for(int ip4=ip3+1; ip4 < nbPrime && (p4=tbPrime[ip4])*Ni3 <=  N; ip4++ ) {
                        int64_t Ni4 = Ni3 * p4 ;
                        NB_lh[NB(Ni4)] -= 16*nbIn ;
                        int64_t p5 ;
                        for(int ip5=ip4+1; ip5 < nbPrime && (p5=tbPrime[ip5])*Ni4 <=  N; ip5++ ) {
                            int64_t Ni5 = Ni4 * p5 ;
                             int64_t p6 ;
                           NB_lh[NB(Ni5)] += 32*nbIn ;
                           for(int ip6=ip5+1; ip6 < nbPrime && (p6=tbPrime[ip6])*Ni5 <=  N; ip6++ ) {
                               int64_t Ni6 = Ni5 * p6 ;
                               NB_lh[NB(Ni6)] -= 64*nbIn ;
                               int64_t p7 ;
                                for(int ip7=ip6+1; ip7 < nbPrime && (p7=tbPrime[ip7])*Ni6 <=  N; ip7++ ) {
                                    int64_t Ni7 = Ni6*p7 ;
                                    NB_lh[NB(Ni7)] += 128*nbIn ;
                                    int64_t p8 ;
                                    for(int ip8=ip7+1; ip8 < nbPrime && (p8=tbPrime[ip8])*Ni7 <=  N; ip8++ ) {
                                        int64_t Ni8 = Ni7*p8 ;
                                        NB_lh[NB(Ni8)] -= 256*nbIn ;
                                        int64_t p9 ;
                                        for(int ip9=ip8+1; ip9 < nbPrime && (p9=tbPrime[ip9])*Ni8 <=  N; ip9++ ) {
                                            int64_t Ni9 = Ni8*p9 ;
                                            NB_lh[NB(Ni9)] += 512*nbIn ;
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
    printf("\n   ");
//    for(int i=0;i<2*NP+1;i++) printf("%lld ",NB_lh[i]);
    {
        int64_t S =0 ;

        for(int i=1;i<=NP;i++) {
             S += NB_lh[NP+1+i] * S2_l[i] ;
//            printf("(%d)=%lld=>%lld ",i,NB_lh[NP+1+i],S);
       }
        printf("\n%.3fs S=%lld nbCont=%d\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC, S,nbContinue);
    }
*/
    
#define S2(i)  ( ( (i) <= NP) ? S2_l[i] : S2_h[N/(i)])
    // #define NoP(i)  ( ( (i) <= NP) ? NBNoP_l[i] : NBNoP_h[N/(i)])
#define NoP(i)  ( ( (i) <= NP) ? (i) : NBNoP_h[N/(i)])

    for(int64_t in = N/NP;in>0;in--) {
        int64_t S2in = 1 ;
        int64_t NBnoPin = 1 ;
        int64_t Pi1 ;
        for(int ip1=0;ip1 < nbPrime && (Pi1=tbPrime[ip1]*in) <= N  ; ip1++) {
            if(N <= NP*Pi1) {
                if(Pi1*2 > N) {
                    NBnoPin++ ;
                    S2in += 2 ;
                } else {
                    int64_t Ni1 = N /Pi1 ;
                    NBnoPin +=  Ni1  ;
                    S2in += S2_l[Ni1] * 2 ;
                }
            } else {
                S2in +=  S2_h[Pi1] *2  ;
                NBnoPin +=  NBNoP_h[Pi1] ;
            }
            int64_t Pi2 ;
            for(int ip2=ip1+1; ip2 < nbPrime && (Pi2=tbPrime[ip2]*Pi1) <=  N; ip2++ ) {
                if(N <= NP*Pi2) {
                    if(Pi2*2 > N) {
                        NBnoPin-- ;
                        S2in -= 4 ;
                    } else {
                        int64_t Ni2 = N /Pi2 ;
                        NBnoPin -=  Ni2  ;
                        S2in -= S2_l[Ni2] * 4 ;
                    }
                } else {
                   S2in -=  S2_h[Pi2] *4  ;
                   NBnoPin -=  NBNoP_h[Pi2] ;
                }
                int64_t Pi3 ;
//                for(int ip3=ip2+1; ip3 < nbPrime && (p3=tbPrime[ip3]) <=  Ni2; ip3++ ) {
                for(int ip3=ip2+1; ip3 < nbPrime && (Pi3=tbPrime[ip3]*Pi2) <=  N ; ip3++ ) {
                    if(N <= NP*Pi3) {
                        if(Pi3*2 > N) {
                            NBnoPin++ ;
                            S2in += 8 ;
                        } else {
                            int64_t Ni3 = N /Pi3 ;
                            NBnoPin +=  Ni3  ;
                            S2in += S2_l[Ni3] * 8 ;
                        }
                    } else {
                        S2in +=  S2_h[Pi3] *8  ;
                        NBnoPin +=  NBNoP_h[Pi3] ;
                    }
                    int64_t Pi4 ;
                    for(int ip4=ip3+1; ip4 < nbPrime && (Pi4=tbPrime[ip4]*Pi3) <=  N; ip4++ ) {
                        if(N <= NP*Pi4) {
                            int64_t Ni4 = N /Pi4 ;
                            NBnoPin -=  Ni4  ;
                            S2in -= S2_l[Ni4] * 16 ;
                        } else {
                            S2in -=  S2_h[Pi4] * 16  ;
                            NBnoPin -=  NBNoP_h[Pi4] ;
                        }
                          int64_t Pi5 ;
                        for(int ip5=ip4+1; ip5 < nbPrime && (Pi5=tbPrime[ip5]*Pi4) <=  N; ip5++ ) {
                            if(N <= NP*Pi5) {
                                int64_t Ni5 = N /Pi5 ;
                                NBnoPin +=  Ni5  ;
                                S2in += S2_l[Ni5] * 32 ;
                            } else {
                                S2in +=  S2_h[Pi5] * 32  ;
                                NBnoPin +=  NBNoP_h[Pi5] ;
                            }
                            int64_t Pi6 ;
                            for(int ip6=ip5+1; ip6 < nbPrime && (Pi6=tbPrime[ip6]*Pi5) <=  N; ip6++ ) {
                                if(N <= NP*Pi6) {
                                    int64_t Ni6 = N /Pi6 ;
                                    NBnoPin -=  Ni6  ;
                                    S2in -= S2_l[Ni6] * 64 ;
                                } else {
                                    S2in -=  S2_h[Pi6] * 64  ;
                                    NBnoPin -=  NBNoP_h[Pi6] ;
                                }
                                int64_t Pi7 ;
                                for(int ip7=ip6+1; ip7 < nbPrime && (Pi7=tbPrime[ip7]*Pi6) <=  N; ip7++ ) {
                                    if(N <= NP*Pi7) {
                                        int64_t Ni7 = N /Pi7 ;
                                        NBnoPin +=  Ni7  ;
                                        S2in += S2_l[Ni7] * 128 ;
                                    } else {
                                        S2in +=  S2_h[Pi7] * 128  ;
                                        NBnoPin +=  NBNoP_h[Pi7] ;
                                    }
                                    int64_t Pi8 ;
                                    for(int ip8=ip7+1; ip8 < nbPrime && (Pi8=tbPrime[ip8]*Pi7) <=  N; ip8++ ) {
                                        if(N <= NP*Pi8) {
                                            int64_t Ni8 = N /Pi8 ;
                                            NBnoPin -=  Ni8  ;
                                            S2in -= S2_l[Ni8] * 256 ;
                                        } else {
                                            S2in -=  S2_h[Pi8] * 256  ;
                                            NBnoPin -=  NBNoP_h[Pi8] ;
                                        }

                                        int64_t Pi9 ;
                                        for(int ip9=ip8+1; ip9 < nbPrime && (Pi9=tbPrime[ip9]*Pi8) <=  N; ip9++ ) {
                                            if(N <= NP*Pi9) {
                                                int64_t Ni9 = N /Pi9 ;
                                                NBnoPin +=  Ni9  ;
                                                S2in += S2_l[Ni9] * 512 ;
                                            } else {
                                                S2in +=  S2_h[Pi9] * 512  ;
                                                NBnoPin +=  NBNoP_h[Pi9] ;
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
        NBNoP_h[in] = NBnoPin ;
 //       S2_h[in] = S2in + 2*(Ni0 - NBnoPin);
        S2_h[in] = S2in  ;

//        if(Ni0 < N2*10) printf("(in=%lld)[%lld]=%lld/%lld(NoP=%lld)) ",in,Ni0,S2_h[in],S2_l[Ni0],NBnoPin);
    }
    printf("Sh[1]=%lld\n",S2_h[1]);
/*
    printf("\nS2_h=");
    for(int i=1;i<=NP;i++) {
        printf("[%d]=%lld ",i,S2_h[i]);
    }
 */
    printf("\n");
    int64_t nbN=NBNoP_h[1] ;
    for(int i=1;i<=N2;i++) {
        NBNoP_h[i]= N/i - NBNoP_h[i]; // donc countSQ[n] contient les nb premiers p >N2 et  n *p < N
//        printf("(%d)%lld->%lld ",i,N/i - NBNoP_h[i],NBNoP_h[i]);
    }
    memset(histPow2,0,maxPow2*sizeof(histPow2[0]));
    for(int n=N2;n>1;n--) {
        //       printf("H[%d]+=(%lld/%d)%lld ",Nb2[n]+1,N,n,countSQ[n]);
        histPow2[ Nb2_l[n]+1] += NBNoP_h[n] ; // on rajoute les p * n
        nbN += NBNoP_h[n] ;
        int d ;
        int64_t pnc = NBNoP_h[n] ;
        for(d=2;d*d<n && d <= N2 ;d++) {
            // on retranche si d diviseur strict de n , les p * n a countSQ[d]
            if((n % d) ==0) {
                NBNoP_h[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 < N2) NBNoP_h[d2] -=pnc ;
            }
        }
        if(d*d == n && d <=N2 ) NBNoP_h[d] -=pnc ;
    }
    histPow2[1] += N-nbN ;

    
 
    
    int64_t S = S2_h[1] ;
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
        S += histPow2[i2] << i2 ;
    }
    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

#define PB708_MAX   100000000000000LL
//#define PB708_MAX   (1LL << 42)
int PB708b1(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t N = PB708_MAX ;
    int maxPow2 = 0 ;
    
    while( (1LL << maxPow2) <= N  ) maxPow2++ ;
    int64_t *histPow2 = calloc(maxPow2*2,sizeof(histPow2[0])) ;
    
    PB708_N *pow2357 = malloc(50000000*sizeof(pow2357[0])) ;
    int nbPow2357 = 0 ;
    pow2357[nbPow2357].val = 1 ;
    pow2357[nbPow2357++].nbp2 = 0 ;
    int n2=1 ;
    for(int64_t pow2=2;pow2<=N;pow2 *= 2 , n2++) {
        pow2357[nbPow2357].val = pow2 ;
        pow2357[nbPow2357++].nbp2 = n2 ;
    }
/*
    CTX_PRIMETABLE64exp * ctxExp= Gen_tablePrime64exp (100000000) ;
    {
        int64_t S = 1 ;
        for(int i2=1;i2<ctxExp->sizeH;i2++) {
            printf("H[%d]=%lld ",i2,ctxExp->histoP2[i2]);
            S += ctxExp->histoP2[i2] << i2 ;
        }
        
        printf("\n%.3fs S=%lld\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC,S) ;//
    }
*/
    
#define PB708_NBP   11
    int prim[PB708_NBP]= {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59} ;
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
    printf("%.3fs Maxp2=%d nbp=%d\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC, maxPow2,nbPow2357);
//    for(int i=0;i<nbPow2357;i++) printf("%lld(%d) ",pow2357[i].val,pow2357[i].nbp2);
//       int64_t nbN = nbPow2357 ;
//    int32_t N2 = (int32_t) Sqrt64(N+1) ;
   int32_t N2 = (int32_t) Sqrt64(N) ;
    int32_t invN2 = (int32_t)(N / N2) ;
    
    uint8_t * Nb2 = calloc(N2+1,sizeof(Nb2[0]));
    for(int p=2;p<=N2;p++) {
        if(Nb2[p]) continue ;
        for(int64_t pp=p;pp<=N2;pp *=p) {
            for(int64_t npp=pp;npp<=N2;npp+=pp) {
                Nb2[npp]++ ;
            }
        }
    }
/*    {
        int64_t S =1 ;
        for(int i =2;i<=400;i++) {
            S += 1 <<  Nb2[i] ;
            printf("%d->%lld%c",i,S,((i % 10) == 0) ? '\n' : ' ') ;
        }
    }
*/
    int64_t * countSQ = calloc(N2+1,sizeof(countSQ[0]));
    int8_t * nb2SQ = calloc(N2+1,sizeof(nb2SQ[0])) ;
    for(int i=1;i<nbPow2357;i++)  {
        int64_t pn0 = pow2357[i].val ;
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
    printf("%.3fs N=%lld Nb primes=%lld N2=%d nbPow2357=%d\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC,N, nbPrime,N2,nbPow2357);
//    for(int n=2;n<N2;n++) {
     {
        for(int i0=0;i0<nbPow2357;i0++) {
 //           printf("%d\n",i0) ;
            int64_t p0 = N/pow2357[i0].val ;
              if(p0 <2 ) break ;
//            printf("%d+%lld ",n,p0) ;
            int nb0 = pow2357[i0].nbp2 ;
            if(p0<N2) {
                  countSQ[p0]++ ;
            }
            else countSQ[N2]++ ;
                
            int ip1 ;
            for(ip1=PB708_NBP+1;ip1<nbPrime;ip1++) {
                T_prime64 p1 = tbPrime[ip1] ;
                if(p1 > p0) break ;
 //               histPow2[nb0+1]++ ; // nbN++ ;
                if(2*p1 > p0) continue ;
                int64_t pn1 = p0/p1 ;
                if(pn1 < N2)     {
                     countSQ[pn1]++ ;
 
                } else countSQ[N2]++ ;
               int ip2 ;
               for(ip2=ip1;ip2<nbPrime;ip2++) {
                    T_prime64 p2 = tbPrime[ip2] ;
                   if(p2 > pn1) break ;
//                   histPow2[nb0+2]++ ; //nbN++ ;
                   if(2*p2 > pn1) continue ;
                   int64_t pn2 = pn1/p2 ;
                   if(pn2 < N2 ) {
                        countSQ[pn2]++ ;
                   } else countSQ[N2]++ ;
                   int ip3 ;
                    for(ip3=ip2;ip3<nbPrime;ip3++) {
                       T_prime64 p3 = tbPrime[ip3] ;
                        if(p3 > pn2) break ;
 //                       histPow2[nb0+3]++ ; //nbN++ ;
                        if(2*p3 > pn2) continue ;
                        int64_t pn3 = pn2/p3 ;
                        if(pn3 < N2) {
                            countSQ[pn3]++ ;
                        } else countSQ[N2]++ ;
                        for(int ip4=ip3;ip4<nbPrime;ip4++) {
                            T_prime64 p4 = tbPrime[ip4] ;
                            if(p4 > pn3) break ;
                            histPow2[nb0+4]++ ; //nbN++ ;
                            if(2*p4 > pn3) continue ;
                            int64_t pn4 = pn3/p4 ;
                            if(pn4 < N2) {
                                countSQ[pn4]++ ;
                            } else countSQ[N2]++ ;
                            for(int ip5=ip4;ip5<nbPrime;ip5++) {
                               T_prime64 p5 = tbPrime[ip5] ;
                                if(p5 > pn4) break ;
                                histPow2[nb0+5]++ ; //nbN++ ;
                                if(2*p5 > pn4) continue ;
                                int64_t pn5 = pn4/p5 ;
                                if(pn5 < N2) {
                                    countSQ[pn5]++ ;
                                } else countSQ[N2]++ ;
                                for(int ip6=ip5;ip6<nbPrime;ip6++) {
                                   T_prime64 p6 = tbPrime[ip6] ;
                                    if(p6 > pn5) break ;
                                    histPow2[nb0+6]++ ; //nbN++ ;
                                    if(2*p6 > pn5) continue ;
                                   int64_t pn6 = pn5/p6 ;
                                    if(pn6 < N2) {
                                       countSQ[pn6]++ ;
                                    } else countSQ[N2]++ ;
                                    for(int ip7=ip6;ip7<nbPrime;ip7++) {
                                        T_prime64 p7 = tbPrime[ip7] ;
                                        if(p7 > pn6) break ;
                                        histPow2[nb0+7]++ ; //nbN++ ;
                                        if(2*p7 > pn6) continue ;
                                        int64_t pn7 = pn6/p7 ;
                                        if(pn7 < N2) {
                                            countSQ[pn7]++ ;
                                        } else countSQ[N2]++ ;
                                        for(int ip8=ip7;ip8<nbPrime;ip8++) {
                                            T_prime64 p8 = tbPrime[ip8] ;
                                            if(p8 > pn7) break ;
                                            histPow2[nb0+8]++ ; //nbN++ ;
                                            if(2*p8 > pn7) continue ;
                                            int64_t pn8 = pn7/p8 ;
                                            if(pn8 < N2) {
                                                countSQ[pn8]++ ;
                                            } else countSQ[N2]++ ;
                                            for(int ip9=ip8;ip9<nbPrime;ip9++) {
                                                T_prime64 p9 = tbPrime[ip9] ;
                                                if(p9 > pn8) break ;
                                                histPow2[nb0+9]++ ; //nbN++ ;
                                                if(2*p9 > pn8) continue ;
                                                int64_t pn9 = pn8/p9 ;
                                                if(pn9 < N2) {
                                                    countSQ[pn9]++ ;
                                                } else countSQ[N2]++ ;
                                            }
                                        }
                                    }
                               }
                            }
                        }
                    }
                    histPow2[nb0+3] += ip3-ip2 ;
               }
                histPow2[nb0+2] += ip2 - ip1 ;
            }
            histPow2[nb0+1] += ip1-PB708_NBP-1 ; // nbN++ ;

            
       }
    }
    for(int n=N2-1;n>=2;n--) {
        // on cumule puisque si un p1*...*pk * n <N c'est vrai pour n plus petit.
         countSQ[n] +=  countSQ[n+1];
//        printf("(%d)=%lld ",n,countSQ[n]);
    }
    printf("N=%lld  \n",N);
//    int64_t S = 1 + (N-nbN)*2;
    int64_t S = 1 ;
    printf("\n Count=");
    for(int i=1;i<=N2;i++) { // counSQ[n] contient le nombre de k tel que k * n < N et k n'a que des facteurs premiers <= N2
        countSQ[i]= N/i - countSQ[i]; // donc countSQ[n] contient les nombres de k tel k *n < N et k a un facteur premier p > N2
 //       printf("(%d)%lld->%lld ",i,N/i - countSQ[i],countSQ[i]);
    }
 /*
    printf("\nP1 ");
    for(int i2=1;i2<maxPow2;i2++) {
        printf("h[%d]=%lld ",i2,histPow2[i2]);
    }
    printf("\n");
*/
 //   int64_t nbN1 = 1 ;
//    for(int i=1;i<maxPow2;i++) nbN1 += histPow2[i] ;

    for(int n=N2;n>1;n--) {
 //       printf("H[%d]+=(%lld/%d)%lld ",Nb2[n]+1,N,n,countSQ[n]);
        histPow2[ Nb2[n]+1] += countSQ[n] ; // on rajoute les p * n
 //       nbN += countSQ[n] ;
        int d ;
        int64_t pnc = countSQ[n] ;
        for(d=2;d*d<n && d <= N2 ;d++) {
            // on retranche si d diviseur strict de n , car si k * n < N ,
            // (k * n/d) * d est une autre decomposition de k * n
            // on ne va garder que p * (k/p *  n) où p est le facteur premier > N2
            if((n % d) ==0) {
                countSQ[d] -= pnc ;
                int64_t d2 = n/d ;
                if(d2 <= N2) countSQ[d2] -=pnc ;
            }
        }
        if(d*d == n && d <=N2 ) countSQ[d] -=pnc ;
   }
    int64_t nbN1 = 1 ;
    for(int i=1;i<maxPow2;i++) nbN1 += histPow2[i] ;

    printf("H[1]=%lld nbN1=%lld \n",histPow2[1],nbN1);
        histPow2[1] += N-nbN1 ;
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

int64_t AlphaL(int32_t n) {
    static int64_t *alphaL = NULL ;
    if(n<0) {
        alphaL = calloc(-n,sizeof(alphaL[0]));
        return 0 ;
    }
    if(alphaL[n]) return alphaL[n] ;
    int64_t alpha = 0 ;
    for(int d=1;d*d <= n;d++ ) {
        alpha += n/d ;
    }
    int64_t sqn = Sqrt32(n) ;
    alpha = 2*alpha - sqn*sqn ;
    alphaL[n] = alpha ;
    return alpha;
}

int64_t AlphaH(int64_t in,int64_t N) {
    static int64_t *alphaH = NULL ;
    if(in<0) {
        alphaH = calloc(-in,sizeof(alphaH[0]));
        return 0 ;
    }
    if(alphaH[in]) return alphaH[in] ;
    int64_t n = N/in ;
    int64_t alpha = 0 ;
    for(int64_t d=1;d*d <= n;d++ ) {
        alpha += n/d ;
    }
    int64_t sqn = Sqrt64(n) ;
    alpha = 2*alpha - sqn*sqn ;
    alphaH[in] = alpha ;
    return alpha;
}

// #define Alpha(i)  ( ( (i) <= NP) ? alphaH[i] : alphaL[N/(i)])
#define Alpha(i)  ( ( (i) < NP) ? AlphaH(i,N) : AlphaL( (int32_t) (N/(i)) ))

int64_t F708(const T_prime64 *tbPrime,int nbPrime,int64_t N,int32_t NP,int64_t divN) {
    int64_t res = Alpha(divN) ;
    for(int ip=0;ip<nbPrime;ip++) {
        int64_t p = tbPrime[ip] ;
//        int64_t p2 = divN*p*p ;
       int64_t p2 = p*p ;
        int64_t IN = N/divN ;
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
    int64_t *histPow2 = calloc(maxPow2,sizeof(histPow2[0])) ;
    
    
    
    int32_t N2 = (int32_t) Sqrt64(N) ;
    int32_t NP = (int32_t) Sqrt64(N+1) ;
 /*
    int64_t * alphaL = malloc((N2+1)*sizeof(alphaL[0])) ;
    int64_t * alphaH = malloc((N2+1)*sizeof(alphaH[0])) ;
    for(int n=1;n<=N2;n++) {
        int64_t alpha = 0 ;
        for(int d=1;d*d <= n;d++ ) {
            alpha += n/d ;
        }
        int64_t sqn = Sqrt32(n) ;
        alpha = 2*alpha - sqn*sqn ;
        alphaL[n] = alpha ;
    }
    printf("%.3fs alphaL computed\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC);
    for(int in=1;in<=N2;in++) {
        int64_t alpha = 0 ;
        int64_t n = N/in ;
        for(int d=1;d*d <= n;d++ ) {
            alpha += n/d ;
        }
        int64_t sqn = Sqrt64(n) ;
        alpha = 2*alpha - sqn*sqn ;
        alphaH[in] = alpha ;
    }
  */
 printf("%.3fs alphaH computed\n",(float)(clock()-pbR->nbClock)/CLOCKS_PER_SEC);
    CTX_PRIMETABLE64 * ctxP = Gen_tablePrime64(N2) ;
    int32_t nbPrime = GetNbPrime64(ctxP) ;
    const T_prime64 *tbPrime = GetTbPrime64(ctxP) ;
    printf("Nb primes=%d N2=%d\n",nbPrime,N2);
    AlphaL(-N2-1);
    AlphaH(-N2-1,N);
    int64_t S = F708(tbPrime,nbPrime,N,NP,1);
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
