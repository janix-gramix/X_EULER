
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
#define PB702_NBT   123457

/*
#define PB702_MAX  15
#define PB702_NBT   123
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

int PB702a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nl= PB702_NBT ;
    int32_t ikMax  ;
    for(ikMax=0;(1<<ikMax) < nl; ikMax++);
    
    uint8_t *isTok = malloc((nl+1)*sizeof(isTok[0]));
    int64_t sumB = 0 ;
    int64_t sumW = 0 ;

    //    for(int il=0;il<nl;il++) {
    int64_t *histoB = calloc(32,sizeof(histoB[0]));
    int64_t *histoW = calloc(32,sizeof(histoW[0]));
    printf("nT==%x=%d =>",nl,nl);
    //************** standard *******************
    for(int il=0;il<nl;il++) { // boucle sur les lignes de triangle
        int32_t nbt = nl-il ;
        int32_t nbTblack = nbt  ;
        int32_t nbTwhite = nbt-1 ;
        memset(isTok,0,nbTblack *sizeof(isTok[0])) ;
//        printf("ligne %d:",il);
        int ik ;
        for(ik = 2;ik<=ikMax;ik++) { // nombre de sauts
            int64_t kpow2 = 1 << ik ;
            int64_t kpow2m1 = kpow2 - 1 ;
            int64_t ilklow = (il * kpow2)/ nl +1 ; // ligne basse
            int64_t ilkhigh = ((il+1)*kpow2) / nl +1 ; //ligne haute
            int nbBlack = 0 ;
            int nbWhite = 0 ;
            for(int64_t ilk = ilklow;ilk<ilkhigh;ilk++) {
 //               printf("[%d]",ilk);
                int64_t yi = ilk*nl ;
                int64_t xiEnd = (nl << ik) - yi;
                int64_t ry = yi & kpow2m1 ;
                int64_t  xi,deltaX ; //  kpow2 * nl
                xi = nl ;
                deltaX = (ilk & 1) ? nl : 2*nl  ;
 //               while(yi+xi < (nl << ik)) {
                while(xi < xiEnd) {
                   int64_t ti = xi >> ik ;
                   int isBlack = (ry+(xi & kpow2m1) <= kpow2) ? 1 : 0 ;
                    if(isBlack) {
                         if(! (isTok[ti] & 1)) {
                             nbBlack++ ;
                             isTok[ti] |= 1 ;
       //                      printf("B");
                         }
                    } else {
                        if(!(isTok[ti] & 2)) {
                            nbWhite++ ;
                            isTok[ti] |= 2 ;
    //                        printf("W");
                        }

                    }
                    xi += deltaX ;
                    
                }
  //              printf("|");
            }
  //          printf("%d=+%dB+%dW ",ik,nbBlack,nbWhite);
            sumB += nbBlack * ik;
            sumW += nbWhite * ik;
            nbTblack -= nbBlack ;
            nbTwhite -= nbWhite ;
 //           histoB[ik] += 2*nbBlack + nbWhite ;
            histoB[ik] += nbBlack  ;
            histoW[ik] +=  nbWhite ;
            if(nbTblack+nbTwhite==0) break ;
//           if(nbTblack==0) break ;
        }
        sumB += nbTblack * (ikMax+1);
        sumW += nbTwhite * (ikMax+1);
 //       printf("%d=+%dB+%dW\n",ikMax+1,nbTblack,nbTwhite);
        histoB[ikMax+1] += nbTblack  ;
        histoW[ikMax+1] +=  nbTwhite ;

    }
    for(int ik=2;ik<=ikMax+1;ik++) printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,2*histoB[ik]+histoW[ik], histoB[ik],histoW[ik]);
 //   for(int ik=2;ik<=ikMax+1;ik++) printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,3*(1<< (ik-2))* ((1<<(ik-1))-1)+ histoB[ik], histoB[ik]
//                                          ,3*(1<< (ik-2))* ((1<<(ik-1))-1) - histoB[ik]);
    printf("Sum=%lld\n",2*sumB+sumW);
    pbR->nbClock = clock() - pbR->nbClock ;

    return 1 ;
}

int PB702c(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nl= PB702_NBT ;
    int32_t ikMax  ;
    for(ikMax=0;(1<<ikMax) < nl; ikMax++);
    
    uint8_t *isTok = malloc((nl+1)*sizeof(isTok[0]));
    int64_t sumB = 0 ;
    int64_t sumW = 0 ;
    
    //    for(int il=0;il<nl;il++) {
    int64_t *histoB = calloc(32,sizeof(histoB[0]));
    int64_t *histoW = calloc(32,sizeof(histoW[0]));
    printf("nT==%x=%d =>",nl,nl);
    int ik = ikMax  ;
    int32_t kpow2 = 1 << ik ;
    int32_t kpow2m1 = kpow2 - 1 ;
    int32_t rl = nl & kpow2m1 ;
    int32_t rl2 = (2*nl) & kpow2m1 ;
    for(int iy=1;iy<kpow2;iy++) {
        int32_t y = (int32_t)(rl * (int64_t)iy) & kpow2m1 ;
        int32_t il = (int32_t) ((iy * nl) >> ik)  ;
        int isDouble = (y + nl < kpow2) ? 1: 0 ;
        if(isDouble) {
            if(iy&1) {
                y += nl ;
            }
            iy++ ; // two lines
            y = kpow2-y ;
            int nbB = 0 ;
            int nbW = 0 ;
           int32_t ix,ky,k,x ;
            for(ky=0;((1<<ky) & y)==0;ky++) ;
            for(ix=1,x=rl;ix+iy <kpow2;ix++ , x = (x+rl) & kpow2m1) {
//            for(ix=2,x= rl2 ;ix+iy <kpow2;ix+=2 , x = (x+rl2) & kpow2m1) {
//                if(ix & 1) k = 0 ;
                if(ix & 1) continue ;
//                else
                for(k=1;((ix >> k) & 1)==0 ;k++) {
                    if(k==ky) break ;
                }
                if(x < y ) { // Black
                    histoB[ik-k]++  ;
                    nbB++ ;
                } else { // white
                    histoW[ik-k]++  ;
                    nbW++ ;
                }
            }
            histoB[ik] += nl  - il - nbB  ;
            histoW[ik] += nl - il - 1 - nbW ;

        } else {
            int32_t ix,x ;
            int32_t ky  , k;
            int nbB = 0 ;
            int nbW = 0 ;
            y = kpow2-y ;
            for(ky=0;((1<<ky) & y)==0;ky++) ;
            if(ky==0) {
    //            for(ix=kpow2-iy-1,x=rl;ix > 0 ;ix-- , x = (x+rl) & kpow2m1) {
                for(x=rl; x != y ; x = (x+rl) & kpow2m1) {
                    if(x < y ) { // Black
                        histoB[ik]++  ;
                        nbB++ ;
                    } else { // white
                        histoW[ik]++  ;
                        nbW++ ;
                    }
                }
     
            } else {
                for(ix=1,x=rl;ix+iy <kpow2;ix++ , x = (x+rl) & kpow2m1) {
                    if(ix & 1) k = 0 ;
                    else for(k=1;((ix >> k) & 1)==0 ;k++) {
                        if(k==ky) break ;
                    }
                    if(x < y ) { // Black
                        histoB[ik-k]++  ;
                        nbB++ ;
                    } else { // white
                        histoW[ik-k]++  ;
                        nbW++ ;
                    }
                }
            }
            histoB[ik+1] += nl  - il - nbB  ;
            histoW[ik+1] += nl - il - 1 - nbW ;
        }
        
    }
    for(int ik=2;ik<=ikMax+1;ik++) {
        sumB += histoB[ik]*ik ;
        sumW += histoW[ik]*ik ;

        printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,2*histoB[ik]+histoW[ik], histoB[ik],histoW[ik]);
    }
    //   for(int ik=2;ik<=ikMax+1;ik++) printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,3*(1<< (ik-2))* ((1<<(ik-1))-1)+ histoB[ik], histoB[ik]
    //                                          ,3*(1<< (ik-2))* ((1<<(ik-1))-1) - histoB[ik]);
    printf("Sum=%lld\n",2*sumB+sumW);
    pbR->nbClock = clock() - pbR->nbClock ;
    
    return 1 ;
}

int PB702c0 (PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int64_t nl= PB702_NBT ;
    int32_t ikMax  ;
    for(ikMax=0;(1<<ikMax) < nl; ikMax++);
    
    uint8_t *isTok = malloc((nl+1)*sizeof(isTok[0]));
    int64_t sumB = 0 ;
    int64_t sumW = 0 ;
    
    //    for(int il=0;il<nl;il++) {
    int64_t *histoB = calloc(32,sizeof(histoB[0]));
    int64_t *histoW = calloc(32,sizeof(histoW[0]));
    printf("nT==%x=%d =>",nl,nl);
    int ik = ikMax - 1 ;
    int32_t kpow2 = 1 << ik ;
    int32_t kpow2m1 = kpow2 - 1 ;
    int32_t rl = nl & kpow2m1 ;
    for(int iy=1;iy<kpow2;iy++) {
        int32_t y = (int32_t)(rl * (int64_t)iy) & kpow2m1 ;
        int32_t ix,x ;
        int32_t ky  , k;
        y = kpow2-y ;
        for(ky=0;((1<<ky) & y)==0;ky++) ;
        if(ky==0) {
            //            for(ix=kpow2-iy-1,x=rl;ix > 0 ;ix-- , x = (x+rl) & kpow2m1) {
            for(x=rl; x != y ; x = (x+rl) & kpow2m1) {
                if(x < y ) { // Black
                    histoB[ik]++  ;
                } else { // white
                    histoW[ik]++  ;
                }
            }
            
        } else {
            for(ix=1,x=rl;ix+iy <kpow2;ix++ , x = (x+rl) & kpow2m1) {
                if(ix & 1) k = 0 ;
                else for(k=1;((ix >> k) & 1)==0 ;k++) {
                    if(k==ky) break ;
                }
                if(x < y ) { // Black
                    histoB[ik-k]++  ;
                } else { // white
                    histoW[ik-k]++  ;
                }
            }
            
        }
        
    }
    for(int ik=2;ik<=ikMax+1;ik++) {
        sumB += histoB[ik]*ik ;
        sumW += histoW[ik]*ik ;
        
        printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,2*histoB[ik]+histoW[ik], histoB[ik],histoW[ik]);
    }
    //   for(int ik=2;ik<=ikMax+1;ik++) printf("\t%d->%lld = 2 x %lldB + %lldW\n",ik,3*(1<< (ik-2))* ((1<<(ik-1))-1)+ histoB[ik], histoB[ik]
    //                                          ,3*(1<< (ik-2))* ((1<<(ik-1))-1) - histoB[ik]);
    printf("Sum=%lld\n",2*sumB+sumW);
    pbR->nbClock = clock() - pbR->nbClock ;
    
    return 1 ;
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
    for(int il=0;il<nl;il++) { // boucle sur les lignes de triangle
        int64_t nbt = nl-(il+1)/2 ;
        int64_t halfNbt = il & 1 ;
        int64_t isHalfAlive = halfNbt ;
        int64_t nbTalive = nbt  ;
        memset(isTok,0,nbTalive *sizeof(isTok[0])) ;
        int ik ;
        for(ik = 2;nbTalive+isHalfAlive;ik++) { // nombre de sauts
            int64_t kpow2 = 1 << ik ;
            int64_t ilklow = (il * kpow2)/ nl +1 ; // ligne basse
             int64_t ilkhigh = ((il+1)*kpow2) / nl +1 ; //ligne haute
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
static uint64_t modPow(int64_t a,int64_t exp,uint64_t mod) {
    uint64_t aPow2 = a % mod ;
    uint64_t aPowExp = (exp & 1) ? aPow2 : 1 ;
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

