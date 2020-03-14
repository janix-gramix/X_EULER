
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

#define PB705_MAX 100000000
//#define PB705_MAX 50
#define PB705_MOD  1000000007

typedef struct PB705_DEC {
    uint8_t nbD[10] ;
    int32_t indHash ;
    int64_t weight ;
    int64_t cost ;
} PB705_DEC ;

typedef struct PB705_WC {
    int64_t weight ;
 //   int64_t cost ;
} PB705_WC ;
#if 0
int PB705(PB_RESULT *pbR) {
    int div1[] = { 1,1 } , div2[]= {2,2,1}, div3[]={2,3,1}, div4[]={3,4,2,1}
    ,div5[]={2,5,1},div6[]={4,6,3,2,1},div7[]={2,7,1},div8[]={4,8,4,2,1},div9[]={3,9,3,1};
    int *Div[10] = {
        NULL,div1,div2,div3,div4,div5,div6,div7,div8,div9
    };
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int64_t SumD[10] ;
    for(int i=0;i<10;i++) SumD[i] = 0 ;
    int64_t cost = 0 ;
    int64_t totalWeight = 1 ;
    pbR->nbClock = clock() ;
    int numDig = 0;
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
            int64_t oldCost = cost ;
            cost = 0 ;
            int nW = Div[d][0] ;
            cost = (cost + nW * oldCost) % PB705_MOD ;
            for(int id=0;id<nW;id++) {
                cost = ( cost + totalWeight * numDig -SumD[Div[d][id+1]]) % PB705_MOD  ;
            }
            SumD[1] = (nW*SumD[1] + totalWeight) %  PB705_MOD;
            switch(nW) {
                case 1 : {
                    for(int idw=2;idw<=9;idw++) {
                        SumD[idw] = (SumD[idw] + totalWeight) %  PB705_MOD;
                    }
                      break ;
                }
                case 2 : {
                    int32_t dd1 = Div[d][1];
                    for(int idw=2;idw<dd1;idw++) {
                        SumD[idw] = (2*SumD[idw] + totalWeight) % PB705_MOD ;
                    }
                    for(int idw=dd1;idw<=9;idw++) {
                        SumD[idw] = (2*SumD[idw] + 2 * totalWeight) % PB705_MOD ;
                    }
                    break ;
                }
                case 3 : {
                    int32_t dd1 = Div[d][1], dd2 = Div[d][2];
                    for(int idw=2;idw<dd2;idw++) {
                        SumD[idw] = (3*SumD[idw] + totalWeight) % PB705_MOD ;
                    }
                    for(int idw=dd2;idw<dd1;idw++) {
                        SumD[idw] = (3*SumD[idw] + 2*totalWeight) % PB705_MOD ;
                    }
                    for(int idw=dd1;idw<=9;idw++) {
                        SumD[idw] = (3*SumD[idw] + 3 * totalWeight) % PB705_MOD ;
                    }
                    break ;
                }
                case 4 : {
                    int32_t dd1 = Div[d][1], dd2 = Div[d][2], dd3 = Div[d][3]  ;
                  for(int idw=2;idw<dd3;idw++) {
                        SumD[idw] = (4*SumD[idw] + totalWeight) % PB705_MOD ;
                    }
                  for(int idw=dd3;idw<dd2;idw++) {
                        SumD[idw] = (4*SumD[idw] + 2*totalWeight) % PB705_MOD ;
                    }
                    for(int idw=dd2;idw<dd1;idw++) {
                        SumD[idw] = (4*SumD[idw] + 3*totalWeight) % PB705_MOD ;
                    }
                    for(int idw=dd1;idw<=9;idw++) {
                        SumD[idw] = (4*SumD[idw] + 4 * totalWeight) % PB705_MOD ;
                    }
                      break ;
                }

            }
            numDig++ ;
            if(nW > 1) totalWeight = ( totalWeight * nW ) % PB705_MOD ;
//            for(int i=1;i<=9;i++) printf("%lld ",SumD[i]);
//            printf("Cost=%lld\n",cost);

        }
    }
    sprintf(pbR->strRes,"%lld",cost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}
#endif
int PB705a(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int64_t SumD1=0, SumD2=0,SumD3=0,SumD4=0,SumD5=0,SumD6=0,SumD7=0,SumD8=0,SumD9=0 ;
    int64_t cost = 0 ;
    int64_t totalWeight = 1 ;
    pbR->nbClock = clock() ;
    int numDig = 0;
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
 //           printf("%d :",d);
            switch(d) {
                case 1 :
                    cost = (  cost + SumD1 ) % PB705_MOD  ;
 /*                   SumD1 = (SumD1)  %  PB705_MOD;
                    SumD2 = (SumD2 ) %  PB705_MOD;
                    SumD3 = (SumD3 ) %  PB705_MOD;
                    SumD4 = (SumD4 ) %  PB705_MOD;
                    SumD5 = (SumD5 ) %  PB705_MOD;
                    SumD6 = (SumD6 ) %  PB705_MOD;
                    SumD7 = (SumD7 ) %  PB705_MOD;
                    SumD8 = (SumD8 ) %  PB705_MOD;
                    SumD9 = (SumD9 ) %  PB705_MOD;
 */                  break ;
                case 2 :
                    cost = (  2*cost+SumD1 +SumD2 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 ) %  PB705_MOD;
                    SumD3 = (2*SumD3 ) %  PB705_MOD;
                    SumD4 = (2*SumD4 ) %  PB705_MOD;
                    SumD5 = (2*SumD5 ) %  PB705_MOD;
                    SumD6 = (2*SumD6 ) %  PB705_MOD;
                    SumD7 = (2*SumD7 ) %  PB705_MOD;
                    SumD8 = (2*SumD8 ) %  PB705_MOD;
                    SumD9 = (2*SumD9) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 3 :
                    cost = (  2*cost +SumD1 +SumD3 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*SumD3 ) %  PB705_MOD;
                    SumD4 = (2*SumD4 ) %  PB705_MOD;
                    SumD5 = (2*SumD5 ) %  PB705_MOD;
                    SumD6 = (2*SumD6 ) %  PB705_MOD;
                    SumD7 = (2*SumD7 ) %  PB705_MOD;
                    SumD8 = (2*SumD8 ) %  PB705_MOD;
                    SumD9 = (2*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 4 :
                    cost = (  3*cost+ +SumD1 +SumD2 +SumD4 ) % PB705_MOD  ;
                    SumD1 = (3*SumD1 + 2*totalWeight) %  PB705_MOD;
                    SumD2 = (3*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (3*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (3*SumD4 ) %  PB705_MOD;
                    SumD5 = (3*SumD5 ) %  PB705_MOD;
                    SumD6 = (3*SumD6 ) %  PB705_MOD;
                    SumD7 = (3*SumD7 ) %  PB705_MOD;
                    SumD8 = (3*SumD8 ) %  PB705_MOD;
                    SumD9 = (3*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 3 ) % PB705_MOD ;
                    break ;
                case 5 :
                    cost = (  2*cost +SumD1 +SumD5 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (2*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (2*SumD5 ) %  PB705_MOD;
                    SumD6 = (2*SumD6 ) %  PB705_MOD;
                    SumD7 = (2*SumD7 ) %  PB705_MOD;
                    SumD8 = (2*SumD8 ) %  PB705_MOD;
                    SumD9 = (2*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 6 :
                    cost = ( 4*cost +SumD1 +SumD2 +SumD3 +SumD6 ) % PB705_MOD  ;
                    SumD1 = (4*SumD1 + 3*totalWeight) %  PB705_MOD;
                    SumD2 = (4*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (4*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (4*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (4*SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (4*SumD6 ) %  PB705_MOD;
                    SumD7 = (4*SumD7 ) %  PB705_MOD;
                    SumD8 = (4*SumD8 ) %  PB705_MOD;
                    SumD9 = (4*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 4 ) % PB705_MOD ;
                    break ;
               case 7 :
                    cost = (  2*cost +SumD1 +SumD7 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (2*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (2*SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (2*SumD6 + totalWeight) %  PB705_MOD;
                    SumD7 = (2*SumD7 ) %  PB705_MOD;
                    SumD8 = (2*SumD8 ) %  PB705_MOD;
                    SumD9 = (2*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 8 :
                    cost = (  4*cost +SumD1 +SumD2 +SumD4 +SumD8 ) % PB705_MOD  ;
                    SumD1 = (4*SumD1 + 3*totalWeight) %  PB705_MOD;
                    SumD2 = (4*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (4*SumD3 + 2*totalWeight) %  PB705_MOD;
                    SumD4 = (4*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (4*SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (4*SumD6 + totalWeight) %  PB705_MOD;
                    SumD7 = (4*SumD7 + totalWeight) %  PB705_MOD;
                    SumD8 = (4*SumD8 ) %  PB705_MOD;
                    SumD9 = (4*SumD9 ) %  PB705_MOD;
                    totalWeight = ( totalWeight * 4 ) % PB705_MOD ;
                    break ;
                case 9 :
                    cost = (  3*cost +SumD1 +SumD3 +SumD9 ) % PB705_MOD  ;
                    SumD1 = (3*SumD1 + 2*totalWeight) %  PB705_MOD;
                    SumD2 = (3*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (3*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (3*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (3*SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (3*SumD6 + totalWeight) %  PB705_MOD;
                    SumD7 = (3*SumD7 + totalWeight) %  PB705_MOD;
                    SumD8 = (3*SumD8 + totalWeight) %  PB705_MOD;
                    SumD9 = (3*SumD9) %  PB705_MOD;
                    totalWeight = ( totalWeight * 3 ) % PB705_MOD ;
                    break ;
            }
            numDig++ ;
 //                       printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld TW=%lld Cost=%lld\n"
 //                             ,SumD1,SumD2,SumD3,SumD4,SumD5
 //                              ,SumD6,SumD7,SumD8,SumD9,totalWeight,cost);
            
        }
    }
    sprintf(pbR->strRes,"%lld",cost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


int PB705(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int64_t SumD1=0, SumD2=0,SumD3=0,SumD4=0,SumD5=0,SumD6=0,SumD7=0,SumD8=0,SumD9=0 ;
    int64_t cost = 0 ;
    int64_t totalWeight = 1 ;
    pbR->nbClock = clock() ;
    int numDig = 0;
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
 //           printf("%d :",d);
            switch(d) {
                case 1 :
                    cost = (  cost + totalWeight * numDig -SumD1 ) % PB705_MOD  ;
                    SumD1 = (SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (SumD6 + totalWeight) %  PB705_MOD;
                    SumD7 = (SumD7 + totalWeight) %  PB705_MOD;
                    SumD8 = (SumD8 + totalWeight) %  PB705_MOD;
                    SumD9 = (SumD9 + totalWeight) %  PB705_MOD;
                    break ;
                case 2 :
                    cost = (  2*cost+2*totalWeight * numDig -SumD1 -SumD2 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*(SumD2 + totalWeight)) %  PB705_MOD;
                    SumD3 = (2*(SumD3 + totalWeight)) %  PB705_MOD;
                    SumD4 = (2*(SumD4 + totalWeight)) %  PB705_MOD;
                    SumD5 = (2*(SumD5 + totalWeight)) %  PB705_MOD;
                    SumD6 = (2*(SumD6 + totalWeight)) %  PB705_MOD;
                    SumD7 = (2*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (2*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (2*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 3 :
                    cost = (  2*cost+2*totalWeight * numDig -SumD1 -SumD3 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*(SumD3 + totalWeight)) %  PB705_MOD;
                    SumD4 = (2*(SumD4 + totalWeight)) %  PB705_MOD;
                    SumD5 = (2*(SumD5 + totalWeight)) %  PB705_MOD;
                    SumD6 = (2*(SumD6 + totalWeight)) %  PB705_MOD;
                    SumD7 = (2*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (2*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (2*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 4 :
                    cost = (  3*cost+3 * totalWeight * numDig -SumD1 -SumD2 -SumD4 ) % PB705_MOD  ;
                    SumD1 = (3*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (3*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (3*SumD3 + 2*totalWeight) %  PB705_MOD;
                    SumD4 = (3*(SumD4 + totalWeight)) %  PB705_MOD;
                    SumD5 = (3*(SumD5 + totalWeight)) %  PB705_MOD;
                    SumD6 = (3*(SumD6 + totalWeight)) %  PB705_MOD;
                    SumD7 = (3*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (3*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (3*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 3 ) % PB705_MOD ;
                    break ;
                case 5 :
                    cost = (  2*cost+2*totalWeight * numDig -SumD1 -SumD5 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (2*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (2*(SumD5 + totalWeight)) %  PB705_MOD;
                    SumD6 = (2*(SumD6 + totalWeight)) %  PB705_MOD;
                    SumD7 = (2*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (2*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (2*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 6 :
                    cost = ( 4*cost+4*totalWeight * numDig -SumD1 -SumD2 -SumD3 -SumD6 ) % PB705_MOD  ;
                    SumD1 = (4*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (4*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (4*SumD3 + 3*totalWeight) %  PB705_MOD;
                    SumD4 = (4*SumD4 + 3*totalWeight) %  PB705_MOD;
                    SumD5 = (4*SumD5 + 3*totalWeight) %  PB705_MOD;
                    SumD6 = (4*(SumD6 + totalWeight)) %  PB705_MOD;
                    SumD7 = (4*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (4*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (4*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 4 ) % PB705_MOD ;
                    break ;
                case 7 :
                    cost = (  2*cost+2*totalWeight * numDig -SumD1 -SumD7 ) % PB705_MOD  ;
                    SumD1 = (2*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (2*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (2*SumD3 + totalWeight) %  PB705_MOD;
                    SumD4 = (2*SumD4 + totalWeight) %  PB705_MOD;
                    SumD5 = (2*SumD5 + totalWeight) %  PB705_MOD;
                    SumD6 = (2*SumD6 + totalWeight) %  PB705_MOD;
                    SumD7 = (2*(SumD7 + totalWeight)) %  PB705_MOD;
                    SumD8 = (2*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (2*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 2 ) % PB705_MOD ;
                    break ;
                case 8 :
                    cost = (  4*cost+4*totalWeight * numDig -SumD1 -SumD2 -SumD4 -SumD8 ) % PB705_MOD  ;
                    SumD1 = (4*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (4*SumD2 + 2*totalWeight) %  PB705_MOD;
                    SumD3 = (4*SumD3 + 2*totalWeight) %  PB705_MOD;
                    SumD4 = (4*SumD4 + 3*totalWeight) %  PB705_MOD;
                    SumD5 = (4*SumD5 + 3*totalWeight) %  PB705_MOD;
                    SumD6 = (4*SumD6 + 3*totalWeight) %  PB705_MOD;
                    SumD7 = (4*SumD7 + 3*totalWeight) %  PB705_MOD;
                    SumD8 = (4*(SumD8 + totalWeight)) %  PB705_MOD;
                    SumD9 = (4*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 4 ) % PB705_MOD ;
                    break ;
                case 9 :
                    cost = (  3*cost+3*totalWeight * numDig -SumD1 -SumD3 -SumD9 ) % PB705_MOD  ;
                    SumD1 = (3*SumD1 + totalWeight) %  PB705_MOD;
                    SumD2 = (3*SumD2 + totalWeight) %  PB705_MOD;
                    SumD3 = (3*SumD3 + 2*totalWeight) %  PB705_MOD;
                    SumD4 = (3*SumD4 + 2*totalWeight) %  PB705_MOD;
                    SumD5 = (3*SumD5 + 2*totalWeight) %  PB705_MOD;
                    SumD6 = (3*SumD6 + 2*totalWeight) %  PB705_MOD;
                    SumD7 = (3*SumD7 + 2*totalWeight) %  PB705_MOD;
                    SumD8 = (3*SumD8 + 2*totalWeight) %  PB705_MOD;
                    SumD9 = (3*(SumD9 + totalWeight)) %  PB705_MOD;
                    totalWeight = ( totalWeight * 3 ) % PB705_MOD ;
                    break ;
            }
            numDig++ ;
//            printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld TW=%lld Cost=%lld\n"
//                   ,SumD1,SumD2,SumD3,SumD4,SumD5
//                   ,SumD6,SumD7,SumD8,SumD9,totalWeight,cost);
            
        }
    }
    sprintf(pbR->strRes,"%lld",cost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#if 0

int PB705d(PB_RESULT *pbR) {
    int div1[] = { 1,1 } , div2[]= {2,2,1}, div3[]={2,3,1}, div4[]={3,4,2,1}
    ,div5[]={2,5,1},div6[]={4,6,3,2,1},div7[]={2,7,1},div8[]={4,8,4,2,1},div9[]={3,9,3,1};
    int *Div[10] = {
        NULL,div1,div2,div3,div4,div5,div6,div7,div8,div9
    };
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int N = 0 ;
    for(int i=0;i<nbPrime;i++) {
        int32_t p = tbPrime[i] ;
        while(p > 0) {
            int ip = p /10 ;
            int d = p - 10 * ip ;
            if(d) N++ ;
            p = ip ;
        }
    }
    printf("N= %d\n",N);
    PB705_WC * ptWC1[10] ;
    PB705_WC * ptWC2[10] ;

    for(int id=1;id<=9;id++) {
        ptWC1[id] =calloc(N+1,sizeof(ptWC1[id][0]));
        ptWC2[id] =calloc(N+1,sizeof(ptWC2[id][0]));
    }
    
    PB705_WC ** curWC = ptWC1 ;
    PB705_WC ** antWC = ptWC2 ;

    
    for(int id=1;id<=9;id++) {
        curWC[id][0].weight = 1 ;
    }
    int64_t cost = 0 ;
    int64_t totalWeight = 1 ;
    pbR->nbClock = clock() ;
    int numDig = 0;
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
            int64_t oldCost = cost ;
            cost = 0 ;
            numDig++ ;
            PB705_WC ** tmp = antWC ;
            antWC = curWC ;
            curWC = tmp ;
            for(int idw=1;idw<=9;idw++) {
                memset(curWC[idw],0,(N+1)*sizeof(curWC[idw][0])) ;
            }
           int nW = Div[d][0] ;
            for(int id=0;id<nW;id++) {
                int32_t dd = Div[d][id+1] ;
 //               printf("%d ",dd);
                for(int idw=1;idw<=9;idw++) {
                    for(int ih=0;ih<numDig;ih++) {
                        if(antWC[idw][ih].weight) {
                            if(idw < dd) {
                                 curWC[idw][ih].weight = (curWC[idw][ih].weight + antWC[idw][ih].weight) % PB705_MOD  ;
                            } else {
                                curWC[idw][ih+1].weight = (curWC[idw][ih+1].weight + antWC[idw][ih].weight) % PB705_MOD  ;
                            }
                        }
                   }
                }
                cost = (cost + oldCost) % PB705_MOD ;
                int64_t deltaCost = 0 ;
                for(int ih=0;ih<numDig;ih++) {
                    if(antWC[dd][ih].weight) {
                       deltaCost = (deltaCost + antWC[dd][ih].weight * ih) % PB705_MOD  ;
                    }
                }
                cost = (cost + totalWeight * (numDig -1)- deltaCost) % PB705_MOD  ;
            }
            if(nW > 1) totalWeight = ( totalWeight * nW ) % PB705_MOD ;
//           printf("Cost=%lld\n",cost);
/*            for(int k = 1;k<=9;k++) {
                for(int ih=0;ih<=N;ih++) printf("%lld ",curWC[k][ih].weight);
                printf("\n");
            }
 */
        }
    }
     sprintf(pbR->strRes,"%lld",cost) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int PB705b(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB705_MAX) ;
    int nbPrime = GetNbPrime(ctxP) ;
    const uint32_t *tbPrime = GetTbPrime(ctxP) ;
    int N = 0 ;
    for(int i=0;i<nbPrime;i++) {
        int32_t p = tbPrime[i] ;
        while(p > 0) {
            int ip = p /10 ;
            int d = p - 10 * ip ;
            if(d) N++ ;
            p = ip ;
        }
    }
    int8_t strDig[13] = { 2 , 3, 5, 7, 1, 1, 1, 3, 1, 7, 1, 9, 0} ;
  
    Nsum * NS = NsumAlloc(N, 9) ; // allocationof structure
    NsumInd size = NsumGetSize(NS,9) ; // get the number of decomposition
    // perfect hash, ks is the number of ni non null
//    NsumInd NsumGetIndex(Nsum *NS,int ks,NsumVal *sum) ;
    printf("N= %d Size=%d\n",N,size);
    
    PB705_DEC * P1dec = malloc(size*sizeof(P1dec[0])) ;
    PB705_DEC * P2dec = malloc(size*sizeof(P2dec[0])) ;

    PB705_DEC * antPdec = P1dec ;
    PB705_DEC * curPdec = P2dec ;

    int32_t *hashPdec = calloc(size+1,sizeof(hashPdec[0])) ;
    
    int nbPdec = 0 ;
    curPdec[nbPdec].weight = 1 ;
    curPdec[nbPdec].cost = 0 ;
    curPdec[nbPdec].indHash = 0 ;
    memset(curPdec[nbPdec].nbD,0,sizeof(curPdec[nbPdec].nbD));
    nbPdec ++ ;
    
    int div1[] = { 1,1 } , div2[]= {2,2,1}, div3[]={2,3,1}, div4[]={3,4,2,1}
    ,div5[]={2,5,1},div6[]={4,6,3,2,1},div7[]={2,7,1},div8[]={4,8,4,2,1},div9[]={3,9,3,1};
    int *Div[10] = {
        NULL,div1,div2,div3,div4,div5,div6,div7,div8,div9
    };
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
                PB705_DEC * tmp = antPdec ;
                antPdec = curPdec ;
                curPdec = tmp ;
                // on reset le hash
                int antNbpdec = nbPdec ;
                nbPdec = 0 ;
                for(int ia=0; ia<antNbpdec;ia++) {
                    hashPdec[antPdec[ia].indHash] = 0 ;
                }
                int nW = Div[d][0] ;
                for(int ia=0;ia<antNbpdec;ia++) {
                    uint8_t nbD[10] ;
                    memcpy(nbD,antPdec[ia].nbD,sizeof(antPdec[ia].nbD)) ;
                    for(int id=0;id<nW;id++) {
                        int32_t dd = Div[d][id+1] ;
                        nbD[dd]++ ;
                        int indH = NsumGetIndex(NS,9,nbD)  ;
                        int indNewP ;
                        if(hashPdec[indH]) {
                            indNewP = hashPdec[indH]-1 ;
                            curPdec[indNewP].weight = (curPdec[indNewP].weight + antPdec[ia].weight)  % PB705_MOD;
                            curPdec[indNewP].cost = (curPdec[indNewP].cost + antPdec[ia].cost) % PB705_MOD ;
                        } else {
                            indNewP = nbPdec ;
                            hashPdec[indH] = ++nbPdec ;
                            memcpy(curPdec[indNewP].nbD,nbD,sizeof(curPdec[indNewP].nbD));
                            curPdec[indNewP].indHash = indH ;
                            curPdec[indNewP].weight = antPdec[ia].weight ;
                            curPdec[indNewP].cost = antPdec[ia].cost ;
                        }
                        nbD[dd]-- ;
                        int dSum = 0 ;
                        for(int idd=9;idd>dd;idd--) {
                            dSum += nbD[idd] ;
                        }
                        curPdec[indNewP].cost = (curPdec[indNewP].cost + dSum * antPdec[ia].weight ) % PB705_MOD  ;
                    }
                }
            
        }
    }
    int64_t S = 0 ;
    for(int ia = 0;ia<nbPdec;ia++) {
        S = (S+curPdec[ia].cost) % PB705_MOD ;
    }

    sprintf(pbR->strRes,"%lld",S) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}



int PB705c(PB_RESULT *pbR) {
    int8_t strDig[] = { 2 , 3, 5, 7, 1, 1, 1, 3, 1, 7, 1, 9, 0} ;
    int nbDiv[9] = { 1 , 2 ,2 , 3 , 2 , 4 , 2 , 4, 3} ;
    int div1[] = { 1,1 } , div2[]= {2,2,1}, div3[]={2,3,1}, div4[]={3,4,2,1}
    ,div5[]={2,5,1},div6[]={4,6,3,2,1},div7[]={2,7,1},div8[]={4,8,4,2,1},div9[]={3,9,3,1};
    int *Div[10] = {
        NULL,div1,div2,div3,div4,div5,div6,div7,div8,div9
    };
    pbR->nbClock = clock() ;
    int64_t nbDig[10] ;
    int64_t oldNbDig[10] ;
    int64_t oldLgxDig[10] ;
    int64_t sumLgxDig[10] ;
    int64_t weight = 1 ;
    int64_t sum = 0 ;
    printf("%d %d %d \n",Div[1][0],Div[1][1],Div[2][0]);
    memset(nbDig,0,sizeof(nbDig));
    memset(sumLgxDig,0,sizeof(sumLgxDig));
    int8_t d ;
    for(int id=0;(d=strDig[id]) != 0;id++ ) {
 /*       memcpy(oldNbDig,nbDig,sizeof(nbDig));
        memcpy(oldLgxDig,sumLgxDig,sizeof(sumLgxDig));
        memset(nbDig,0,sizeof(nbDig));
        memset(sumLgxDig,0,sizeof(sumLgxDig));
  */
        int nW = Div[d][0] ;
        sum *= nW ;
        for(int n=0;n<Div[d][0];n++) {
            int dd = Div[d][n+1] ;
            for(int idd=dd+1;idd<=9;idd++) {
                sum += sumLgxDig[idd] ;
            }
        }

        for(int i=1;i<=9;i++) {
            sumLgxDig[i] = nW * (sumLgxDig[i]+nbDig[i]) ;
            nbDig[i] = nW * nbDig[i] ;
        }
         for(int n=0;n<Div[d][0];n++) {
            int dd = Div[d][n+1] ;
              nbDig[dd] += weight ;
             sumLgxDig[dd] += weight ;
        }
        weight *= nW ;
        printf("%d =>W=%lld,S=%lld",d,weight,sum) ;
        for(int i = 1;i<=9;i++) printf("[%lld x %lld]",sumLgxDig[i],nbDig[i]);
        printf("\n");
    }

    return 1 ;
}
#endif
