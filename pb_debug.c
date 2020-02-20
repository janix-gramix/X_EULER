
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define M_PI 3.14159265358979323846264338327950288

#include "euler_utils.h"
#include "DC3.h"
#include "PB_other.h"
#include "ukkonen.h"

#define PB701_L 7
#define PB701_H 7
#define PB701_MAXCONNECT    5

typedef uint32_t NsumInd  ;
typedef uint8_t NsumVal ;
typedef struct Nsum {
    NsumInd N ;
    NsumInd k ;
    NsumInd *index ;
} Nsum ;

#define NS_IND(k,n) ((k) * NS->N +(n))
Nsum * NsumAlloc(NsumInd N, NsumInd k) {
    if(k<2) return NULL ;
    Nsum * NS = calloc(1,sizeof(NS[0])) ;
    NS->N = N+1 ;
    NS->k = k ;
    NS->index = malloc((k-1)*(N+1)*sizeof(NS->index[0]));
    // init k-1 serie
    NsumInd *ids = NS->index+ NS_IND(k-2,0) ;
    ids[0] = 0 ;
    for(int in=1;in<=N;in++){
        ids[in] = ids[in-1] + N + 2 - in ;
    }
    for(int ik=k-3;ik>=0;ik--) {
        NsumInd *idsAnt = ids ;
        ids = NS->index+ NS_IND(ik,0) ;
        ids[0] = 0 ;
        NsumInd indAnt = idsAnt[N]+1 ;
        for(int in=1;in<=N;in++){
            ids[in] = ids[in-1] + indAnt - idsAnt[in-1] ;
        }
    }
    return NS ;
}

Nsum * NsumFree(Nsum * NS) {
    if(NS != NULL) {
        free(NS->index);
        free(NS);
    }
    return NULL ;
}
NsumInd NsumGetSize(Nsum *NS,int ks) {
    if(ks<=1) {
        return ks ? NS->N : 1 ;
    }
    return NS->index[NS_IND(NS->k-ks,NS->N-1)]+1 ;
}

NsumInd NsumGetIndex(Nsum *NS,int ks,NsumVal *sum) {
   
    int k = NS->k ;
    int N = NS->N ;
    NsumInd * index = NS->index + (k-ks)*N ;
    switch(ks) {
        case 0: return 0;
        case 1: return sum[0] ;
        case 2: return index[sum[0]]+sum[1] ;
        case 3: return index[sum[0]] + index[N+sum[0]+sum[1]]-index[N+sum[0]]+sum[2];
        case 4: return index[sum[0]] + index[N+sum[0]+sum[1]]-index[N+sum[0]]
            + index[2*N+sum[0]+sum[1]+sum[2]]-index[2*N+sum[0]+sum[1]]  +sum[3];
        default :
            {
                NsumInd ksum = *sum++ ;
                NsumInd ind = index[ksum] ;
                for(int ik=1;ik< ks-1;ik++) {
                    NsumInd ksumNext = ksum + *sum++ ;
                    ind += index[N*ik+ksumNext]-index[N*ik+ksum];
                    ksum= ksumNext ;
                }
                return ind+ *sum ;
            }

    }
}

typedef union Contact {
 //   uint64_t   Glob ;
    uint16_t   Set[PB701_MAXCONNECT] ;
} Contact ;

typedef struct PB701_PP {
    int num ;
    int globalContact ;
    int nbSet ;
    Contact Cntc ;
 //   uint8_t sizeC[PB701_MAXCONNECT] ;
} PB701_PP ;

typedef struct PB701_CONN {
    int globalContact ;
    int nbSet ;
    Contact Cntc ;
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

typedef struct PB701_COH {
    uint16_t numPP ;
    NsumVal Surf[PB701_MAXCONNECT+1] ;
#if defined PB701_BIG
    __int128 prob ;
#else
    int64_t prob ;
#endif
} PB701_COH ;





int PB701a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int dim = 1 << PB701_L ;
    int FirstPPByFrontier[dim+1] ;
    PB701_PP *PP = malloc(4000*sizeof(PP[0])) ;
    PB701_CONN PC[dim] ;
    int nbPP = 0 ;
    int nbBit[512] ;
    PB701_INTER interSec[dim*dim] ;
    for(int i=0;i<dim;i++) {
        int nb = 0 ;
        for(int j=0;j<PB701_L;j++) {
            if(i & (1 << j)) nb++ ;
        }
        nbBit[i] = nb ;
    }
    
    for(int ic=0;ic<dim;ic++) {
        for(int j=0;j<dim;j++) {
            int k = ic & j ;
            if(k==0) {
                interSec[ic*dim+j].newContact = 0 ;
                interSec[ic*dim+j].surfIncrement = 0 ;
                continue ;
            }
            int n = 0 ;
            for(int m=0;m<PB701_L;m++) {
                if(m>0  && ( (1<<(m-1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if(m<PB701_L-1  && ( (1<<(m+1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
            }
            for(int m=PB701_L-1;m>=0;m--) {
                if(m>0  && ( (1<<(m-1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if(m<PB701_L-1  && ( (1<<(m+1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if((1 <<m) & k)  n++  ;
            }
            
            interSec[ic*dim+j].surfIncrement = n ;
            interSec[ic*dim+j].newContact = (uint16_t) k ;
            //           printf("I(%x,%x)=%d(%x) ",ic,j,n,k) ;
        }
    }
    
    FirstPPByFrontier[0] = 0 ;
    PP[nbPP].num  = nbPP ;
    PP[nbPP].globalContact = 0;
    PP[nbPP].nbSet = 0 ;
    memset(PP[nbPP].Cntc.Set,0,sizeof(PP[nbPP].Cntc.Set));
    nbPP++ ;
    int IsSymetric[dim] ;
    IsSymetric[0] = 0 ;
    for(int i=1;i<dim;i++) {
        {
            int iSym =0 ;
            for(int ib=0;ib<PB701_L;ib++) {
                iSym |= ((i >> ib) & 1) << (PB701_L -ib -1) ;
            }
            if(iSym == i) IsSymetric[i] = 0 ;
            else if(i > iSym) IsSymetric[i] = 1;
            else IsSymetric[i] = -1;
        }
        FirstPPByFrontier[i] = nbPP ;
        PB701_PP *curPP = PP + nbPP ;
        Contact  contact ;
        memset(contact.Set,0,sizeof(contact.Set));
        int nbConnect = 0 ;
        int n = i & 1 ;
        for(int j=1;j<=PB701_L;j++) {
            if(i & (1 << j)) {
                n++ ;
            } else {
                if(n) {
                    contact.Set[nbConnect] = (uint16_t)((1 << j) - (1<<(j-n))) ;
                    nbConnect++ ;
                }
                n = 0 ;
            }
        }
        if(nbConnect > PB701_MAXCONNECT ){
            nbConnect = PB701_MAXCONNECT ;
            printf("PB_maxconnect\n");
        }
        PC[i].globalContact = i ;
        PC[i].nbSet = nbConnect ;
        PC[i].Cntc = contact ;
#define INIT_PART(n)   curPP = PP + nbPP++ ; curPP->num = nbPP-1 ;  curPP->globalContact = i ; memset(curPP->Cntc.Set,0,sizeof(curPP->Cntc.Set)) ; curPP->nbSet = (n)
        
#define DST(i)  curPP->Cntc.Set[i]
#define SRC(i)  contact.Set[i]
        
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
        FirstPPByFrontier[i+1] = nbPP ;
        for(int ip = FirstPPByFrontier[i];ip<FirstPPByFrontier[i+1];ip++ ) {
            int nbs = PP[ip].nbSet ;
            for(int is = 1;is<nbs;is++) {
                for(int js=is;js<nbs;js++) {
                    if(PP[ip].Cntc.Set[is-1] > PP[ip].Cntc.Set[js]) {
                        uint16_t tmp = PP[ip].Cntc.Set[is-1] ;
                        PP[ip].Cntc.Set[is-1] = PP[ip].Cntc.Set[js] ;
                        PP[ip].Cntc.Set[js] = tmp ;
                    }
                }
            }
        }
    }
    FirstPPByFrontier[dim] = nbPP ;
    PB701_TRANSF *PT = calloc(nbPP*dim,sizeof(PT[0]));
    for(int ip=0;ip<nbPP;ip++) {
        for(int j=0;j<dim;j++) {
            int nbSet = 0 ;
            Contact  contact ;
            memset(contact.Set,0,sizeof(contact.Set));
            int globalContact = 0 ;
            //            printf("%x X n%2.2d %x=",j,PP[ip].num, PP[ip].globalContact) ;
            //            for(int iss=0;iss<PP[ip].nbSet;iss++) printf("%x+",PP[ip].Cntc.Set[iss]);
            //            printf("=> ");
            for(int is=0;is<PP[ip].nbSet;is++) {
                uint16_t newContact = interSec[PP[ip].Cntc.Set[is]*dim+j].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) {
                        if((newContact & contact.Set[jj])) {
                            newContact = interSec[(newContact | contact.Set[jj])*dim+j].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                contact.Set[l-1] = contact.Set[l] ;
                            }
                            nbSet--; jj-- ;
                            contact.Set[nbSet] = 0 ;
                        }
                    }
                    int l ;
                    for(l=nbSet-1;l>=0 && contact.Set[l] > newContact  ;l--) {
                        contact.Set[l+1] = contact.Set[l] ;
                    }
                    contact.Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }
            if(globalContact != j) {
                globalContact ^= j ;
                for(int k=0;k<PC[j].nbSet;k++) {
                    if((PC[j].Cntc.Set[k] & globalContact) == PC[j].Cntc.Set[k] ){
                        
                        int l ;
                        for(l=nbSet-1;l>=0 && contact.Set[l] > PC[j].Cntc.Set[k]   ;l--) {
                            contact.Set[l+1] = contact.Set[l] ;
                        }
                        contact.Set[l+1] = PC[j].Cntc.Set[k] ;
                        nbSet++ ;
                    }
                }
                
            }
            int in ;
            for(in=FirstPPByFrontier[j];in<FirstPPByFrontier[j+1];in++) {
                if(memcmp(contact.Set,PP[in].Cntc.Set,sizeof(contact.Set))==0) break ;
                //             if(contact.Glob==PP[in].Cntc.Glob) break ;
            }
            if(in == FirstPPByFrontier[j+1]) {
                printf("!!!!!!!! PB !!!!!!Ident\n");
                in = FirstPPByFrontier[j] ;
            }
            PT[ip*dim+j].numPP = in ;
            for(int id=0;id < PP[in].nbSet ; id++) PT[ip*dim+j].delta[id] = nbBit [PP[in].Cntc.Set[id]] ;
            for(int is=0;is<PP[ip].nbSet;is++) {
                if(PP[ip].Cntc.Set[is] & j) {
                    for(int id=0;id < PP[in].nbSet ; id++){
                        if( PP[ip].Cntc.Set[is] & PP[in].Cntc.Set[id] ) {
                            PT[ip*dim+j].numCv[is] = id+1 ;
                            break ;
                        }
                        if(id == PP[in].nbSet) printf("!") ;
                    }
                }
                
            }
        }
    }
    int firstCohByLevl[10] ;
    PB701_COH *Coh= calloc(500000000,sizeof(Coh[0])) ;
    int nbCoh = 0;
    
    //   Nsum * NS ;
    NsumInd * p1IndCoh = calloc(nbPP,sizeof(p1IndCoh[0])) ;
    NsumInd * p2IndCoh = calloc(nbPP,sizeof(p2IndCoh[0])) ;
    NsumInd *curPPindCoh = p1IndCoh ;
    NsumInd *antPPindCoh = p2IndCoh ;
    //    NsumInd *curUsedIch = calloc(1,sizeof(curUsedIch[0]));
    firstCohByLevl[0]= 0 ;
    //   curUsedIch[nbCoh++] = 0 ;
    nbCoh = 1;
    firstCohByLevl[1]= nbCoh ;
    Coh[0].prob = 1 ;
    Coh[0].numPP = 0 ;
    double S1 = 0 ;
    int level ;
    //    NS= NsumAlloc(0, 5) ;
    //    printf("SizeNS=%d\n",NsumGetSize(NS,1));
    //    for(level=0;level< PB701_L-1 ;level++) {
    for(level=0;level< PB701_H-1 ;level++) {
        //    for(level=0;level< PB701_L-2 ;level++) {
        S1 = 0 ;
        printf("**** %dx%d ",PB701_L,level+1) ;
        
        
        NsumInd *tmp = antPPindCoh ;
        antPPindCoh = curPPindCoh ;
        curPPindCoh = tmp ;
        int N=PB701_L*(level+1) ;
        
        Nsum * NS= NsumAlloc(N,PB701_MAXCONNECT+1) ;
        /*
         int SizeProb = 0 ;
         int curSize ;
         free(curUsedIch);
         for(int ip=0;ip<nbPP;ip++) {
         curPPindCoh[ip] = SizeProb ;
         curSize = NsumGetSize(NS, PP[ip].nbSet+1) ;
         SizeProb += curSize ;
         }
         curUsedIch = calloc(SizeProb,sizeof(curUsedIch[0])) ;
         */
        printf("Ich %d -> %d ",firstCohByLevl[level],firstCohByLevl[level+1]);
        
        int64_t Sprob = 0 ;
        int maxIch = 0 ;
        //       int maxj = (level) ? dim : dim/2 ;
        for(int j=0;j<dim;j++) {
            int isDouble = 0 ;
            if(level==PB701_H-2) {
                if(IsSymetric[j]< 0) continue ;
                isDouble = IsSymetric[j] ? 1 : 0 ;
            }
            int SizeIch = 0 ;
            int curSize ;
            //            for(int ip=0;ip<nbPP;ip++) {
            for(int ip=FirstPPByFrontier[j];ip<FirstPPByFrontier[j+1];ip++) {
                curPPindCoh[ip] = SizeIch ;
                curSize = NsumGetSize(NS, PP[ip].nbSet+1) ;
                SizeIch += curSize ;
            }
            NsumInd * curUsedIch = calloc(SizeIch,sizeof(curUsedIch[0])) ;
            if(SizeIch > maxIch) maxIch = SizeIch ;
            for(int ich=firstCohByLevl[level];ich < firstCohByLevl[level+1];ich++) {
                PB701_COH *antCoh = Coh+ich ;
                int ip = antCoh->numPP ;
                PB701_TRANSF *trf = PT + ip*dim + j ;
                int id = trf->numPP ;
                PB701_COH *newCoh = Coh+nbCoh ;
                NsumVal * surf = newCoh->Surf;
                newCoh->numPP = id ;
                surf[0] = antCoh->Surf[0] ;
                for(int is=0;is<PP[id].nbSet;is++) {
                    surf[is+1] = trf->delta[is] ;
                }
                for(int is=0;is<PP[ip].nbSet;is++) {
                    if(trf->numCv[is]) {
                        surf[trf->numCv[is]] += antCoh->Surf[is+1];
                    } else {
                        if(antCoh->Surf[is+1] >  surf[0] ) surf[0] = antCoh->Surf[is+1] ;
                    }
                }
                NsumInd numCoh = curPPindCoh[id] + NsumGetIndex(NS, PP[id].nbSet+1, surf) ;
                int newIch = curUsedIch[numCoh] ;
#if defined PB701_BIG
                __int128  prob = (isDouble) ? 2*antCoh->prob  : antCoh->prob ;
#else
                int64_t prob = (isDouble) ? 2*antCoh->prob  : antCoh->prob ;
#endif
                if(newIch == 0) {
                    curUsedIch[numCoh]= nbCoh++ ;
                    //                    newCoh->prob = antCoh->prob ;
                    newCoh->prob = prob ;
                } else {
                    //                     Coh[newIch].prob += antCoh->prob ;
                    Coh[newIch].prob += prob ;
                }
                //                Sprob += antCoh->prob ;
                Sprob += prob ;
            }
            free(curUsedIch);
            
        }
        
        firstCohByLevl[level+2] = nbCoh ;
        printf(" ->MaxIch=%d ",maxIch);
#if defined PB701_BIG
        __int128 S = 0 ;
#else
        uint64_t S = 0 ;
#endif
        
        for(int ich=firstCohByLevl[level+1];ich < firstCohByLevl[level+2];ich++) {
            PB701_COH * newCoh = Coh + ich  ;
            int maxS = newCoh->Surf[0];
            for(int is=0;is<PP[newCoh->numPP].nbSet;is++) {
                if(newCoh->Surf[is+1] > maxS) maxS = newCoh->Surf[is+1]  ;
            }
            S += maxS * newCoh->prob ;
        }
        NS=NsumFree(NS);
        double Sres ;
#if defined PB701_BIG
        Sres = (double)S / (double) (1LL << (level*PB701_L )) / (double)(1LL << PB701_L) ;
        printf("S=%.8f S=%llu%llu Sprob2=%llx\n",Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL),Sprob);
#else
        Sres = SHigh / (double) (1LL << ((level+1)*PB701_L-16 )) + SLow / (double) (1LL << (level*PB701_L )) / (double)(1LL << PB701_L) ;
        printf("S=%.8f S=%lld  Sprob2=%llx\n",Sres,SLow+ (SHigh << 16),Sprob);
#endif
    }
    
    double Sres ;
    {
        //       S = 0 ;
        
#if defined PB701_BIG
        __int128 S = 0 ;
#else
        uint64_t SLow = 0;
        uint64_t SHigh = 0 ;
#endif
        printf("**** %dx%d\n",PB701_L,level+1) ;
        //      free(antUsedIch) ;
        //        free(curUsedIch) ;
        int64_t Sprob2 = 0 ;
        printf("Ich %d -> %d \n",firstCohByLevl[level],firstCohByLevl[level+1]);
        for(int ich=firstCohByLevl[level];ich < firstCohByLevl[level+1];ich++) {
            PB701_COH *antCoh = Coh + ich  ;
            int ip = antCoh->numPP ;
            //           int isSymetric = IsSymetric[PP[ip].globalContact] ;
            for(int j=0;j<dim;j++) { // [0,1[
                //               if(IsSymetric[j]<=0) continue ;
                PB701_TRANSF *trf = PT + ip*dim + j ;
                int id = trf->numPP ;
                Sprob2 += antCoh->prob ;
                NsumVal surf[PB701_MAXCONNECT+1] ;
                surf[0] = antCoh->Surf[0] ;
                for(int is=0;is<PP[id].nbSet;is++) {
                    surf[is+1] = trf->delta[is] ;
                }
                for(int is=0;is<PP[ip].nbSet;is++) {
                    if(trf->numCv[is]) {
                        surf[trf->numCv[is]] += antCoh->Surf[is+1];
                    } else {
                        if(antCoh->Surf[is+1] >  surf[0] ) surf[0] = antCoh->Surf[is+1] ;
                    }
                }
                int maxS = surf[0];
                for(int is=0;is<PP[id].nbSet;is++) {
                    if(surf[is+1] > maxS) maxS = surf[is+1]  ;
                }
                
                
                /*               if(isSymetric) {
                 SLow += 2* maxS * (antCoh->prob & 0xffff) ;
                 SHigh += 2* maxS * (antCoh->prob >> 16) ;
                 } else
                 */
                
#if defined PB701_BIG
                S += maxS * antCoh->prob;
#else
                {
                    SLow += maxS * (antCoh->prob & 0xffff) ;
                    SHigh += maxS * (antCoh->prob >> 16) ;
                }
#endif
                //               S += maxS * antCoh->prob ;
            }
        }
        free(Coh);
#if defined PB701_BIG
        Sres = (double)S / (double) (1LL << (level*PB701_L )) / (double)(1LL << PB701_L) ;
        printf("LxH=%dx%d S=%.8f S=%llu%llu S1=%.8f  Sprob2=%llx\n",PB701_L,PB701_H
               ,Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL),S1,Sprob2);
#else
        Sres = SHigh / (double) (1LL << ((level+1)*PB701_L-16 )) + SLow / (double) (1LL << (level*PB701_L )) / (double)(1LL << PB701_L) ;
        printf("S=%.8f S=%lld S1=%.8f  Sprob2=%llx\n",Sres,SLow+ (SHigh << 16),S1,Sprob2);
#endif
    }
    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.8f",Sres);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

typedef struct PB701_COHa {
    uint16_t numPP ;
    NsumVal Surf[PB701_MAXCONNECT+1] ;
    int64_t prob ;
} PB701_COHa ;


typedef struct PB701_PPxPP {
//    int     numPP ;
    int nbSet ;
    int nbSet0 ;
    int nbSet1 ;
    int8_t delta[PB701_MAXCONNECT] ;
    uint8_t numCv0[PB701_MAXCONNECT] ;
    uint8_t numCv1[PB701_MAXCONNECT] ;
} PB701_PPxPP ;


int PB701b(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int dim = 1 << PB701_L ;
    int FirstPPByFrontier[dim+1] ;
    PB701_PP *PP = malloc(4000*sizeof(PP[0])) ;
    PB701_CONN PC[dim] ;
    int nbPP = 0 ;
    int nbBit[512] ;
    PB701_INTER interSec[dim*dim] ;
    for(int i=0;i<dim;i++) {
        int nb = 0 ;
        for(int j=0;j<PB701_L;j++) {
            if(i & (1 << j)) nb++ ;
        }
        nbBit[i] = nb ;
    }
    
    for(int ic=0;ic<dim;ic++) {
        for(int j=0;j<dim;j++) {
            int k = ic & j ;
            if(k==0) {
                interSec[ic*dim+j].newContact = 0 ;
                interSec[ic*dim+j].surfIncrement = 0 ;
                continue ;
            }
            int n = 0 ;
            for(int m=0;m<PB701_L;m++) {
                if(m>0  && ( (1<<(m-1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if(m<PB701_L-1  && ( (1<<(m+1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
            }
            for(int m=PB701_L-1;m>=0;m--) {
                if(m>0  && ( (1<<(m-1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if(m<PB701_L-1  && ( (1<<(m+1)) & k ) && ( (1<<m)& j) ) k |= 1 << m ;
                if((1 <<m) & k)  n++  ;
            }
            
            interSec[ic*dim+j].surfIncrement = n ;
            interSec[ic*dim+j].newContact = (uint16_t) k ;
            //           printf("I(%x,%x)=%d(%x) ",ic,j,n,k) ;
        }
    }
    
    FirstPPByFrontier[0] = 0 ;
    PP[nbPP].num  = nbPP ;
    PP[nbPP].globalContact = 0;
    PP[nbPP].nbSet = 0 ;
    memset(PP[nbPP].Cntc.Set,0,sizeof(PP[nbPP].Cntc.Set));
    nbPP++ ;
    int IsSymetric[dim] ;
    IsSymetric[0] = 0 ;
    PC[0].globalContact = 0 ;
    PC[0].nbSet = 0 ;
    memset(PC[0].Cntc.Set,0,sizeof(PC[0].Cntc.Set)) ;

    for(int i=1;i<dim;i++) {
        {
            int iSym =0 ;
            for(int ib=0;ib<PB701_L;ib++) {
                iSym |= ((i >> ib) & 1) << (PB701_L -ib -1) ;
            }
            if(iSym == i) IsSymetric[i] = 0 ;
            else if(i > iSym) IsSymetric[i] = 1;
            else IsSymetric[i] = -1;
        }
        FirstPPByFrontier[i] = nbPP ;
        PB701_PP *curPP = PP + nbPP ;
        Contact  contact ;
        memset(contact.Set,0,sizeof(contact.Set));
         int nbConnect = 0 ;
        int n = i & 1 ;
        for(int j=1;j<=PB701_L;j++) {
            if(i & (1 << j)) {
                n++ ;
            } else {
                if(n) {
                    contact.Set[nbConnect] = (uint16_t)((1 << j) - (1<<(j-n))) ;
                    nbConnect++ ;
                }
                n = 0 ;
            }
        }
        if(nbConnect > PB701_MAXCONNECT ){
            nbConnect = PB701_MAXCONNECT ;
            printf("PB_maxconnect\n");
        }
        PC[i].globalContact = i ;
        PC[i].nbSet = nbConnect ;
        PC[i].Cntc = contact ;
#define INIT_PART(n)   curPP = PP + nbPP++ ; curPP->num = nbPP-1 ;  curPP->globalContact = i ; memset(curPP->Cntc.Set,0,sizeof(curPP->Cntc.Set)) ; curPP->nbSet = (n)

#define DST(i)  curPP->Cntc.Set[i]
#define SRC(i)  contact.Set[i]

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
        FirstPPByFrontier[i+1] = nbPP ;
         for(int ip = FirstPPByFrontier[i];ip<FirstPPByFrontier[i+1];ip++ ) {
            int nbs = PP[ip].nbSet ;
            for(int is = 1;is<nbs;is++) {
                for(int js=is;js<nbs;js++) {
                    if(PP[ip].Cntc.Set[is-1] > PP[ip].Cntc.Set[js]) {
                        uint16_t tmp = PP[ip].Cntc.Set[is-1] ;
                        PP[ip].Cntc.Set[is-1] = PP[ip].Cntc.Set[js] ;
                        PP[ip].Cntc.Set[js] = tmp ;
                    }
                }
            }
        }
    }
    FirstPPByFrontier[dim] = nbPP ;
    PB701_TRANSF *PT = calloc(nbPP*dim,sizeof(PT[0]));
    for(int ip=0;ip<nbPP;ip++) {
        for(int j=0;j<dim;j++) {
            int nbSet = 0 ;
            Contact  contact ;
            memset(contact.Set,0,sizeof(contact.Set));
            int globalContact = 0 ;
            //            printf("%x X n%2.2d %x=",j,PP[ip].num, PP[ip].globalContact) ;
            //            for(int iss=0;iss<PP[ip].nbSet;iss++) printf("%x+",PP[ip].Cntc.Set[iss]);
            //            printf("=> ");
            for(int is=0;is<PP[ip].nbSet;is++) {
                uint16_t newContact = interSec[PP[ip].Cntc.Set[is]*dim+j].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) {
                        if((newContact & contact.Set[jj])) {
                            newContact = interSec[(newContact | contact.Set[jj])*dim+j].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                contact.Set[l-1] = contact.Set[l] ;
                            }
                            nbSet--; jj-- ;
                            contact.Set[nbSet] = 0 ;
                        }
                    }
                    int l ;
                    for(l=nbSet-1;l>=0 && contact.Set[l] > newContact  ;l--) {
                        contact.Set[l+1] = contact.Set[l] ;
                    }
                    contact.Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }
            if(globalContact != j) {
                globalContact ^= j ;
                for(int k=0;k<PC[j].nbSet;k++) {
                    if((PC[j].Cntc.Set[k] & globalContact) == PC[j].Cntc.Set[k] ){
                        
                        int l ;
                        for(l=nbSet-1;l>=0 && contact.Set[l] > PC[j].Cntc.Set[k]   ;l--) {
                            contact.Set[l+1] = contact.Set[l] ;
                        }
                        contact.Set[l+1] = PC[j].Cntc.Set[k] ;
                        nbSet++ ;
                    }
                }
                
            }
            int in ;
            for(in=FirstPPByFrontier[j];in<FirstPPByFrontier[j+1];in++) {
                if(memcmp(contact.Set,PP[in].Cntc.Set,sizeof(contact.Set))==0) break ;
   //             if(contact.Glob==PP[in].Cntc.Glob) break ;
            }
            if(in == FirstPPByFrontier[j+1]) {
                printf("!!!!!!!! PB !!!!!!Ident\n");
                in = FirstPPByFrontier[j] ;
            }
            PT[ip*dim+j].numPP = in ;
           for(int id=0;id < PP[in].nbSet ; id++) PT[ip*dim+j].delta[id] = nbBit [PP[in].Cntc.Set[id]] ;
            for(int is=0;is<PP[ip].nbSet;is++) {
                if(PP[ip].Cntc.Set[is] & j) {
                    for(int id=0;id < PP[in].nbSet ; id++){
                        if( PP[ip].Cntc.Set[is] & PP[in].Cntc.Set[id] ) {
                            PT[ip*dim+j].numCv[is] = id+1 ;
                            break ;
                        }
                        if(id == PP[in].nbSet) printf("!") ;
                    }
                }
                
            }
        }
    }
    PB701_PPxPP *PPxPP = calloc(nbPP*nbPP,sizeof(PPxPP[0]));
    for(int ip0=0;ip0<nbPP;ip0++) {
        int global0 = PP[ip0].globalContact ;
        for(int ip1=0;ip1<nbPP;ip1++) {
            int global1 = PP[ip1].globalContact ;
            if(global0 != global1) continue ;
            int nbSet = 0 ;
            Contact  contact ;
            memset(contact.Set,0,sizeof(contact.Set));
            int globalContact = 0 ;
            for(int is0=0;is0<PP[ip0].nbSet;is0++) {
                uint16_t newContact = interSec[PP[ip0].Cntc.Set[is0]*dim+global1].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) {
                        if((newContact & contact.Set[jj])) {
                            newContact = interSec[(newContact | contact.Set[jj])*dim+global1].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                contact.Set[l-1] = contact.Set[l] ;
                            }
                            nbSet--; jj-- ;
                            contact.Set[nbSet] = 0 ;
                        }
                    }
                    int l ;
                    for(l=nbSet-1;l>=0 && contact.Set[l] > newContact  ;l--) {
                        contact.Set[l+1] = contact.Set[l] ;
                    }
                    contact.Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }
   
            for(int is1=0;is1<PP[ip1].nbSet;is1++) {
                uint16_t newContact = interSec[PP[ip1].Cntc.Set[is1]*dim+global1].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) {
                        if((newContact & contact.Set[jj])) {
                            newContact = interSec[(newContact | contact.Set[jj])*dim+global1].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                contact.Set[l-1] = contact.Set[l] ;
                            }
                            nbSet--; jj-- ;
                            contact.Set[nbSet] = 0 ;
                        }
                    }
                    int l ;
                    for(l=nbSet-1;l>=0 && contact.Set[l] > newContact  ;l--) {
                        contact.Set[l+1] = contact.Set[l] ;
                    }
                    contact.Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }

            
            
            
            for(int k=0;k<PC[global1].nbSet;k++) {
                if((PC[global1].Cntc.Set[k] & globalContact) == 0 ){
                    
                    int l ;
                    for(l=nbSet-1;l>=0 && contact.Set[l] > PC[global1].Cntc.Set[k]   ;l--) {
                        contact.Set[l+1] = contact.Set[l] ;
                    }
                    contact.Set[l+1] = PC[global1].Cntc.Set[k] ;
                    nbSet++ ;
                    globalContact |= PC[global1].Cntc.Set[k] ;
                }
            }
 /*
            int in ;
            for(in=FirstPPByFrontier[globalContact];in<FirstPPByFrontier[globalContact+1];in++) {
                if(memcmp(contact.Set,PP[in].Cntc.Set,sizeof(contact.Set))==0) break ;
                //             if(contact.Glob==PP[in].Cntc.Glob) break ;
            }
            if(in == FirstPPByFrontier[globalContact+1]) {
                printf("!!! PB !!!!PPxPP \n");
                in = FirstPPByFrontier[globalContact] ;
            }
  */
            PPxPP[ip0*dim+ip1].nbSet = nbSet ;
            PPxPP[ip0*dim+ip1].nbSet0 = PP[ip0].nbSet ;
            PPxPP[ip0*dim+ip1].nbSet1 = PP[ip1].nbSet ;
           int is01[PB701_MAXCONNECT] ;
            memset(is01,0,sizeof(is01));
            for(int is0=0;is0<PP[ip0].nbSet;is0++) {
                if(PP[ip0].Cntc.Set[is0] & globalContact) {
                    for(int id=0;id < nbSet ; id++){
                        if( PP[ip0].Cntc.Set[is0] & contact.Set[id] ) {
                            PPxPP[ip0*dim+ip1].numCv0[is0] = id+1 ;
                            is01[id] |= 1 ;
                            break ;
                        }
                        if(id == nbSet) printf("!") ;
                    }
                }
            }
            for(int is1=0;is1<PP[ip1].nbSet;is1++) {
                if(PP[ip1].Cntc.Set[is1] & globalContact) {
                    for(int id=0;id < nbSet ; id++){
                        if( PP[ip1].Cntc.Set[is1] & contact.Set[id] ) {
                            PPxPP[ip0*dim+ip1].numCv1[is1] = id+1 ;
                            is01[id] |=2 ;
                            break ;
                        }
                        if(id == nbSet) printf("!") ;
                    }
                }
            }
            for(int id=0;id < nbSet ; id++){
                if(is01[id]==3) {
                    PPxPP[ip0*dim+ip1].delta[id] = -nbBit[contact.Set[id]] ;
                } else {
                    PPxPP[ip0*dim+ip1].delta[id] = 0 ;
                }
            }
/*            printf("%x=%x: PP0=",global0,global1);
            for(int is=0;is<PP[ip0].nbSet;is++) {printf("%x,",PP[ip0].Cntc.Set[is]);}
            printf(" PP1=");
            for(int is=0;is<PP[ip1].nbSet;is++) {printf("%x,",PP[ip1].Cntc.Set[is]);}
            printf(" PPfinal(%d)=",in);
            for(int is=0;is<PP[in].nbSet;is++) {printf("%x,",PP[in].Cntc.Set[is]);}
            printf(" cv0=") ;
            for(int is=0;is<PP[ip0].nbSet;is++) {printf("%d ",PPxPP[ip0*dim+ip1].numCv0[is]);}
            printf(" cv1=") ;
            for(int is=0;is<PP[ip1].nbSet;is++) {printf("%d ",PPxPP[ip0*dim+ip1].numCv1[is]);}
            printf(" delta=") ;
            for(int is=0;is<PP[ip1].nbSet;is++) {printf("%d ",PPxPP[ip0*dim+ip1].delta[is]);}
            printf("\n");
*/        }
    }

    

    int firstCohByLevl[10] ;
    PB701_COHa *Coh= calloc(100000000,sizeof(Coh[0])) ;
    int nbCoh = 0;
     NsumInd * p1IndCoh = calloc(nbPP,sizeof(p1IndCoh[0])) ;
    NsumInd * p2IndCoh = calloc(nbPP,sizeof(p2IndCoh[0])) ;
    NsumInd *curPPindCoh = p1IndCoh ;
    NsumInd *antPPindCoh = p2IndCoh ;
    firstCohByLevl[0]= 0 ;
    nbCoh = 1;
    firstCohByLevl[1]= nbCoh ;
    Coh[0].prob = 1 ;
    Coh[0].numPP = 0 ;
    double S1 = 0 ;
    int level ;
//    for(level=0;level< PB701_H-1 ;level++) {
   for(level=0;level< (PB701_H+1)/2 ;level++) {
        S1 = 0 ;
        printf("**** %dx%d ",PB701_L,level+1) ;
        
 
        NsumInd *tmp = antPPindCoh ;
        antPPindCoh = curPPindCoh ;
        curPPindCoh = tmp ;
         int N=PB701_L*(level+1) ;
   
        Nsum * NS= NsumAlloc(N,PB701_MAXCONNECT+1) ;

        printf("Ich %d -> %d ",firstCohByLevl[level],firstCohByLevl[level+1]);

        int64_t Sprob = 0 ;
         int maxIch = 0 ;
         for(int j=0;j<dim;j++) {
            int isDouble = 0 ;
            if(level==(PB701_H+1)/2-1) {
                if(IsSymetric[j]< 0) continue ;
                isDouble = IsSymetric[j] ? 1 : 0 ;
            }

            int SizeIch = 0 ;
            int curSize ;
            for(int ip=FirstPPByFrontier[j];ip<FirstPPByFrontier[j+1];ip++) {
                curPPindCoh[ip] = SizeIch ;
                curSize = NsumGetSize(NS, PP[ip].nbSet+1) ;
                SizeIch += curSize ;
            }
            NsumInd * curUsedIch = calloc(SizeIch,sizeof(curUsedIch[0])) ;
            if(SizeIch > maxIch) maxIch = SizeIch ;

            
            
            

            for(int ich=firstCohByLevl[level];ich < firstCohByLevl[level+1];ich++) {
                PB701_COHa *antCoh = Coh+ich ;
                int ip = antCoh->numPP ;
                PB701_TRANSF *trf = PT + ip*dim + j ;
                int id = trf->numPP ;
                PB701_COHa *newCoh = Coh+nbCoh ;
                NsumVal * surf = newCoh->Surf;
                newCoh->numPP = id ;
                surf[0] = antCoh->Surf[0] ;
                for(int is=0;is<PP[id].nbSet;is++) {
                    surf[is+1] = trf->delta[is] ;
                }
               for(int is=0;is<PP[ip].nbSet;is++) {
                    if(trf->numCv[is]) {
                        surf[trf->numCv[is]] += antCoh->Surf[is+1];
                    } else {
                        if(antCoh->Surf[is+1] >  surf[0] ) surf[0] = antCoh->Surf[is+1] ;
                    }
                }
                NsumInd numCoh = curPPindCoh[id] + NsumGetIndex(NS, PP[id].nbSet+1, surf) ;
                 int newIch = curUsedIch[numCoh] ;
//                int64_t prob = (isDouble) ? 2*antCoh->prob  : antCoh->prob ;
               int64_t prob =  antCoh->prob ;
                if(newIch == 0) {
                     curUsedIch[numCoh]= nbCoh++ ;
                    newCoh->prob = prob ;
                 } else {
                    Coh[newIch].prob += prob ;
                }
                Sprob += prob ;
            }
            free(curUsedIch);

        }
        
        firstCohByLevl[level+2] = nbCoh ;
        printf(" ->MaxIch=%d ",maxIch);
        uint64_t S = 0 ;
      
        for(int ich=firstCohByLevl[level+1];ich < firstCohByLevl[level+2];ich++) {
            PB701_COHa * newCoh = Coh + ich  ;
            int maxS = newCoh->Surf[0];
            for(int is=0;is<PP[newCoh->numPP].nbSet;is++) {
                if(newCoh->Surf[is+1] > maxS) maxS = newCoh->Surf[is+1]  ;
            }
            S += maxS * newCoh->prob ;
        }
        NS=NsumFree(NS);
       double Sres ;
       Sres = (double)S / (double) (1LL << (level*PB701_L )) / (double)(1LL << PB701_L) ;
       printf("S=%.8f S=%llu%llu Sprob2=%llx\n",Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL),Sprob);
    }
    
    double Sres ;
    {
        __int128 S = 0 ;
        int64_t Sprob2 = 0 ;
        int levelAnt = PB701_H - level ;
        printf("Ich %dx%d X %dx%d-> %d X %d\n",PB701_L,levelAnt,PB701_L,level,
               firstCohByLevl[levelAnt+1]-firstCohByLevl[levelAnt],
               firstCohByLevl[level+1]-firstCohByLevl[level]);
        int antGlobal  = 0 ;
        int lastich = firstCohByLevl[level] ;
        for(int ich=firstCohByLevl[level];ich < firstCohByLevl[level+1];ich++) {
            PB701_COHa *coh = Coh + ich  ;
            PB701_PP *PPi = PP + coh->numPP ;
            int ip = coh->numPP ;
             int global = PPi->globalContact ;
            if(global == antGlobal && ich+1 != firstCohByLevl[level+1]) continue;
            else {
                if(ich+1 == firstCohByLevl[level+1]) ich++ ;
                int isDouble = (IsSymetric[antGlobal]==0) ? 0 : 1 ;
 //               printf("%x[%d,%d[ ",antGlobal,lastich,ich);fflush(stdout) ;
                for(int ich0 = lastich;ich0<ich;ich0++) {
                    PB701_COHa *coh0 = Coh + ich0  ;
                    int ip0 = coh0->numPP ;
                    int64_t prob0 = coh0->prob ;
                    for(int ich1=ich0;ich1 < ich;ich1++) {
                        PB701_COHa *coh1 = Coh + ich1  ;
                        int ip1 = coh1->numPP ;
                        int64_t prob1 = coh1->prob ;
                        PB701_PPxPP * ppxpp = PPxPP+ip0*dim + ip1 ;
                        uint16_t surf[PB701_MAXCONNECT+1] ;
                          surf[0] = coh0->Surf[0] ;
                        for(int is=0;is<ppxpp->nbSet;is++) {
                            surf[is+1] = ppxpp->delta[is] ;
                        }
                        if(coh1->Surf[0] > surf[0] ) surf[0] = coh1->Surf[0] ;
                        for(int is0=0;is0<ppxpp->nbSet0;is0++) {
                            surf[ppxpp->numCv0[is0]] += coh0->Surf[is0+1];
                        }
                        for(int is1=0;is1<ppxpp->nbSet1;is1++) {
                            surf[ppxpp->numCv1[is1]] += coh1->Surf[is1+1];
                        }
                        int maxS = surf[0];
                          for(int is=0;is<ppxpp->nbSet;is++) {
//                              if(surf[is+1]-nbBit[PP[id].Cntc.Set[is]] > maxS) maxS = surf[is+1]-nbBit[PP[id].Cntc.Set[is]]  ;
                          if(surf[is+1] > maxS) maxS = surf[is+1]  ;
                        }
                        if(ich0==ich1) {
                            if(isDouble) {
                                S += 2 * maxS * prob0 * prob1;
                                Sprob2 += 2* prob0 * prob1 ;
                            } else {
                                S += maxS * prob0 * prob1;
                                Sprob2 += prob0 * prob1 ;

                            }
                        } else {
                            if(isDouble) {
                                S += 4 * maxS * prob0 * prob1;
                                Sprob2 += 4 * prob0 * prob1 ;
                           } else {
                                S += 2 * maxS *  prob0 * prob1;
                                Sprob2 += 2 * prob0 * prob1 ;
                            }

                        }
                  }
                    
                }
                antGlobal = global ;
                lastich = ich ;
            }
        }
        Sres = (double)S / (double) (1LL << (level*PB701_L )) / (double)(1LL <<(levelAnt* PB701_L)) ;
        printf("LxH=%dx%d S=%.8f S=%llu%llu S1=%.8f  Sprob2=%llx\n",PB701_L,PB701_H
               ,Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL),S1,Sprob2);
    }

    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.8f",Sres);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


