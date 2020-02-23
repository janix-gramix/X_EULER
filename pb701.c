//
//  pb701
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#define PB701_L 7
#define PB701_H 7
#define PB701_MAXCONNECT    5

#define PB701_BIG
#if defined(PB701_BIG)
    typedef __int128 Prob701 ;
    #define PB701_MAX_LENGTH_COL    9
#else
    typedef uint64_t Prob701 ;
    #define PB701_MAX_LENGTH_COL    7
#endif


// exact hash for decompostion N = n1+ n2 + n3 ... + nk.
typedef uint32_t NsumInd  ;
typedef uint8_t NsumVal ;
typedef struct Nsum {
    NsumInd N ;
    NsumInd k ;
    NsumInd *index ;
} Nsum ;

Nsum * NsumAlloc(NsumInd N, NsumInd k) ;
Nsum * NsumFree(Nsum * NS) ;
NsumInd NsumGetSize(Nsum *NS,int ks) ;
NsumInd NsumGetIndex(Nsum *NS,int ks,NsumVal *sum) ;



typedef struct PB701_PP {
    int num ;
    int globalContact ;
    int nbSet ;
    uint16_t   Set[PB701_MAXCONNECT]  ;
} PB701_PP ;


typedef struct PB701_TRANSF {
    int     numPP ;
    uint8_t numCv[PB701_MAXCONNECT] ;
    uint8_t delta[PB701_MAXCONNECT] ;
} PB701_TRANSF ;




typedef struct PB701_STATE {
    uint16_t numPP ;    // PP corresponding
    NsumVal  Surf[PB701_MAXCONNECT+1] ; // surface for each set (including hidden set)
    Prob701  prob ;     // weight of the state
} PB701_STATE ;

typedef struct PB701_DATA {
    int nbCol ;   // number of differents column ( 1<<sizeColumn)
    int nbPP ;  //  number of column partitions in connected parts
    int *isSymetric  ; // 0 <=> mask symetric, +1 im > sym(im) , -1 im < sym(im)
    PB701_PP *PP  ; // description for each column partition in connected masks
    int *firstPPByMask ; // index of first PP by column.
    PB701_TRANSF *PT ; // description of transformation to add a column to a partition
} PB701_DATA ;

void PB701_Init(PB701_DATA *data701, int sizeColumn) ;

int main(int argc, const char * argv[]) {
    clock_t nbClock = clock() ;
    PB701_DATA data701 ;
    int sizeColumn = PB701_L ;
    int maxNbColumn =  PB701_H ;
    if(argc > 1){
        sizeColumn =atoi(argv[1]) ;
        if(sizeColumn > PB701_MAX_LENGTH_COL) sizeColumn=PB701_MAX_LENGTH_COL ;
        if(argc > 2) {
            maxNbColumn = atoi(argv[2]);
            if(maxNbColumn > PB701_MAX_LENGTH_COL) maxNbColumn = PB701_MAX_LENGTH_COL ;
            if(maxNbColumn < sizeColumn) {
                int tmp = sizeColumn ;
                sizeColumn = maxNbColumn ;
                maxNbColumn = tmp ;
            }
        }
    } else {
        printf("Args : [sizeColumn (%d)] [nbColumn (%d)]\n",sizeColumn,maxNbColumn);
    }

    PB701_Init(&data701,sizeColumn) ; // precompute all transformation to add a column
    int nbCol = data701.nbCol ; // recover
    PB701_PP * PP = data701.PP ;
    PB701_TRANSF * PT = data701.PT;
    
    PB701_STATE *State= malloc(300000000*sizeof(State[0])) ;
    int nbState = 0;
    NsumInd * PPindState = calloc(data701.nbPP,sizeof(PPindState[0])) ;
    int firstNbState , lastNbState ;
    firstNbState = 0 ;
    State[0].prob = 1 ; State[0].numPP = 0 ;
    nbState = lastNbState = 1;
    double Sres ;
    for(int nbColumn=1;nbColumn< maxNbColumn ;nbColumn++) { // loops on column
        printf("\t%.3fs %dx%d ",(float)(clock()-nbClock)/CLOCKS_PER_SEC ,sizeColumn,nbColumn) ;
        int N=sizeColumn*nbColumn ;
        Nsum * NS= NsumAlloc(N,PB701_MAXCONNECT+1) ;
        printf("States [%d -> %d[ ",firstNbState,lastNbState);
        Prob701 S=0 , ST = 0 ;
        int maxHash = 0 ;
        int maxState = 0 ;
        for(int jc=0;jc<nbCol;jc++) { // loop on differents column the current added column
            int isDouble = 0 ;
            if(nbColumn==maxNbColumn-1) { // if last one, check symetry
                if(data701.isSymetric[jc]< 0) continue ;
                isDouble = data701.isSymetric[jc] ? 1 : 0 ;
            }
            // compute size of exact hash for current added column
            int SizeHash = 0 ;
            int firstStateForColumn = nbState ;
            for(int ip=data701.firstPPByMask[jc];ip<data701.firstPPByMask[jc+1];ip++) {
                PPindState[ip] = SizeHash ;
                SizeHash += NsumGetSize(NS,PP[ip].nbSet+1)  ;
            }
            // alloc hash states merge
            NsumInd * hashState = calloc(SizeHash,sizeof(hashState[0])) ;
            if(SizeHash > maxHash) maxHash = SizeHash ;
            for(int istate=firstNbState;istate < lastNbState;istate++) { // loop on state from preceedind column
                PB701_STATE *antState = State+istate ;
                int ip = antState->numPP ;
                PB701_TRANSF *trf = PT + ip*nbCol + jc ;
                int id = trf->numPP ;
                PB701_STATE *newState = State+nbState ;
                NsumVal * surf = newState->Surf;
                newState->numPP = id ;
                surf[0] = antState->Surf[0] ;
                for(int is=0;is<PP[id].nbSet;is++) {
                    surf[is+1] = trf->delta[is] ;
                }
                for(int is=0;is<PP[ip].nbSet;is++) {
                    if(trf->numCv[is]) {
                        surf[trf->numCv[is]] += antState->Surf[is+1];
                    } else {
                        if(antState->Surf[is+1] >  surf[0] ) surf[0] = antState->Surf[is+1] ;
                    }
                }
                NsumInd numHash = PPindState[id] + NsumGetIndex(NS, PP[id].nbSet+1, surf) ;
                int newIstate = hashState[numHash] ;
                if(newIstate == 0) { // not used state
                    hashState[numHash]= nbState++ ;
                    newState->prob = isDouble ? 2*antState->prob  : antState->prob ;
                } else { // used state , cumulate only prob
                    State[newIstate].prob += isDouble ? 2*antState->prob  : antState->prob  ;
                }
            }
            free(hashState);
            for(int istate = firstStateForColumn ; istate<nbState; istate++) {
                PB701_STATE * newState = State + istate  ;
                int maxS = newState->Surf[0];
                for(int is=0;is<PP[newState->numPP].nbSet;is++) {
                    if(newState->Surf[is+1] > maxS) maxS = newState->Surf[is+1]  ;
                }
                S += maxS * newState->prob ;
            }
            if(nbColumn==maxNbColumn-1) { // last column, special case not, necessary to save new states)
                for(int istate=firstStateForColumn;istate < nbState;istate++) {
                    PB701_STATE *antState = State + istate  ;
                    int ip = antState->numPP ;
                    for(int jc=0;jc<nbCol;jc++) { // [0,1[
                        PB701_TRANSF *trf = PT + ip*nbCol + jc ;
                        int id = trf->numPP ;
                        NsumVal surf[PB701_MAXCONNECT] ;
                        int maxS = antState->Surf[0] ;
                        for(int is=0;is<PP[id].nbSet;is++) {
                            surf[is] = trf->delta[is] ;
                        }
                        for(int is=0;is<PP[ip].nbSet;is++) {
                            if(trf->numCv[is]) {
                                surf[trf->numCv[is]-1] += antState->Surf[is+1];
                            } else {
                                if(antState->Surf[is+1] >  maxS ) maxS = antState->Surf[is+1] ;
                            }
                        }
                        for(int is=0;is<PP[id].nbSet;is++) {
                            if(surf[is] > maxS) maxS = surf[is]  ;
                        }
                        ST += maxS * antState->prob;
                    }
                }
                if(nbState > maxState) maxState = nbState ;
                nbState = firstStateForColumn ; // remins, the used states for last pass are not necessary
            }
        }
        
        firstNbState = lastNbState ;
        lastNbState = nbState ;
        printf(" ->MaxHash=%d ",maxHash);
        NS=NsumFree(NS);
        int exp2 = (nbColumn)*sizeColumn ;
        Sres = (double)S / (double) (1LL << (exp2/2)) / (double)(1LL << ((exp2+1)/2)) ;
#if defined PB701_BIG
        printf("S=%.8f S=%llu%llu\n",Sres,(int64_t)(S / 1000000000000000LL), (int64_t)(S % 1000000000000000LL));
#else
        printf("S=%.8f S=%llu\n",Sres,S);
#endif
        if(nbColumn == maxNbColumn-1) {
            exp2 = (nbColumn+1)*sizeColumn ;
            Sres = (double)ST / (double) (1LL << (exp2/2 )) / (double)(1LL << ((exp2+1)/2)) ;
#if defined PB701_BIG
            printf("%.3fs LxH=%dx%d MaxStates=%d S=%.8f S=%llu%llu\n",(float)(clock()-nbClock)/CLOCKS_PER_SEC
                   ,sizeColumn,maxNbColumn,maxState
                   ,Sres,(int64_t)(ST / 1000000000000000LL), (int64_t)(ST % 1000000000000000LL));
#else
            printf("S=%.8f S=%llu\n",Sres,ST);
#endif
        }
    }
    free(PPindState);
    free(State);
    return 0 ;
}

//***************************************************************************
typedef struct PB701_INTER {
    uint8_t surfIncrement ;
    uint16_t newContact ;
} PB701_INTER ;

typedef struct PB701_CONN {
    int globalContact ;
    int nbSet ;
    uint16_t Set[PB701_MAXCONNECT]  ;
    uint8_t sizeC[PB701_MAXCONNECT] ;
}PB701_CONN ;

// precompute transtions to add a column.
void PB701_Init(PB701_DATA *data701, int sizeColumn) {
    data701->nbCol = 1 << sizeColumn ;
    int *nbBit ;
    nbBit = malloc(data701->nbCol*sizeof(nbBit[0]));
    // compute nb bits by mask
    for(int ic=0;ic<data701->nbCol;ic++) {
        int nb = 0 ;
        for(int j=0;j < sizeColumn ;j++) {
            if(ic & (1 << j)) nb++ ;
        }
        nbBit[ic] = nb ;
    }
    // compute symetry of mask
    data701->isSymetric = malloc(data701->nbCol*sizeof(data701->isSymetric[0]));
    for(int ic=0;ic<data701->nbCol;ic++) {
        int iSym =0 ;
        for(int ib=0;ib<sizeColumn;ib++) {
            iSym |= ((ic >> ib) & 1) << (sizeColumn -ib -1) ;
        }
        if(iSym == ic) data701->isSymetric[ic] = 0 ;
        else if(ic > iSym) data701->isSymetric[ic] = 1;
        else data701->isSymetric[ic] = -1;
    }
    // compute intersection of 2 mask ,extended to connected parts for each mask.
    PB701_INTER *interSec ; // intersection of 2 masks extended to connexion parts
    interSec=malloc(data701->nbCol*data701->nbCol*sizeof(interSec[0]));
    for(int ic=0;ic<data701->nbCol;ic++) {
        for(int jc=0;jc<data701->nbCol;jc++) {
            int km = ic & jc ;
            if(km==0) {
                interSec[ic*data701->nbCol+jc].newContact = 0 ;
                interSec[ic*data701->nbCol+jc].surfIncrement = 0 ;
                continue ;
            }
            // extension to connected parts
            for(int m=0;m<sizeColumn;m++) {
                if(m>0  && ( (1<<(m-1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
                if(m<sizeColumn-1  && ( (1<<(m+1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
            }
            for(int m=sizeColumn-1;m>=0;m--) {
                if(m>0  && ( (1<<(m-1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
                if(m<sizeColumn-1  && ( (1<<(m+1)) & km ) && ( (1<<m)& jc) ) km |= 1 << m ;
            }
            interSec[ic*data701->nbCol+jc].surfIncrement = nbBit[km] ;
            interSec[ic*data701->nbCol+jc].newContact = (uint16_t) km ;
        }
    }
    // compute connected partition for each column
    PB701_CONN *PC = malloc(data701->nbCol*sizeof(PC[0])) ;
    int nbPPexp = 0 ;
    for(int ic=0;ic<data701->nbCol;ic++) {
        uint16_t Set[PB701_MAXCONNECT]  ;
        memset(Set,0,sizeof(Set));
        int nbConnect = 0 ;
        int n = ic & 1 ;
        for(int j=1;j<=sizeColumn;j++) {
            if(ic & (1 << j)) {
                n++ ;
            } else {
                if(n) {
                    Set[nbConnect] = (uint16_t)((1 << j) - (1<<(j-n))) ;
                    nbConnect++ ;
                    n=0 ;
                }
            }
        }
        assert(nbConnect <= sizeColumn);
        int ppByConnected[PB701_MAXCONNECT+1] = { 1,1,2,5,14,42} ;
        nbPPexp += ppByConnected[nbConnect] ;
        PC[ic].globalContact = ic ;
        PC[ic].nbSet = nbConnect ;
        memcpy(PC[ic].Set,Set,sizeof(Set)) ;
    }
    data701->PP = malloc(nbPPexp*sizeof(data701->PP[0]));
    data701->firstPPByMask = malloc((data701->nbCol+1)*sizeof(data701->firstPPByMask[0])) ;
    // compute different partitions of column in connected parts
    // dont add impossible configurations as aaabbbaabbb (entrelacement not planar)
    data701->nbPP = 0 ;
    data701->firstPPByMask[0] = data701->nbPP ;
    for(int ic=0;ic<data701->nbCol;ic++) {
        PB701_PP *curPP = data701->PP + data701->nbPP ;
        uint16_t Set[PB701_MAXCONNECT]  ;
        memcpy(Set,PC[ic].Set,sizeof(Set)) ;
        int nbConnect = PC[ic].nbSet ;
#define INIT_PART(n)   curPP = data701->PP + data701->nbPP ; curPP->num = data701->nbPP++ ;  curPP->globalContact = ic ; memset(curPP->Set,0,sizeof(curPP->Set)) ; curPP->nbSet = (n)
        
#define DST(i)  curPP->Set[i]
#define SRC(i)  Set[i]
        
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
            case 0 :
                INIT_PART(0) ;
                break ;
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
        // sort each partition set for future comparison between partitions
        data701->firstPPByMask[ic+1] = data701->nbPP ;
        for(int ip = data701->firstPPByMask[ic];ip<data701->firstPPByMask[ic+1];ip++ ) {
            int nbs = data701->PP[ip].nbSet ;
            for(int is = 1;is<nbs;is++) {
                for(int js=is;js<nbs;js++) {
                    if(data701->PP[ip].Set[is-1] > data701->PP[ip].Set[js]) {
                        uint16_t tmp = data701->PP[ip].Set[is-1] ;
                        data701->PP[ip].Set[is-1] = data701->PP[ip].Set[js] ;
                        data701->PP[ip].Set[js] = tmp ;
                    }
                }
            }
        }
        data701->firstPPByMask[ic+1] = data701->nbPP ;
    }
    // compute for each partition every transformation to add a new column
    data701->PT = calloc(data701->nbPP*data701->nbCol,sizeof(data701->PT[0]));
    for(int ip=0;ip<data701->nbPP;ip++) { // loop on partition
        for(int jc=0;jc<data701->nbCol;jc++) { // lopp on supplementary partition
            int nbSet = 0 ;
            uint16_t Set[PB701_MAXCONNECT]  ;
            memset(Set,0,sizeof(Set));
            int globalContact = 0 ;
            for(int is=0;is<data701->PP[ip].nbSet;is++) {
                // for each connected part search the new set of connected parts in new column
                // remove all intersecting precedent sets and add the new set.
                // respect ascending order in sets for future comparison
                uint16_t newContact = interSec[data701->PP[ip].Set[is]*data701->nbCol+jc].newContact ;
                if(newContact){
                    int jj ;
                    for(jj=0;jj<nbSet;jj++) { // remove intersecting sets
                        if((newContact & Set[jj])) {
                            newContact = interSec[(newContact | Set[jj]) * data701->nbCol+jc].newContact ;
                            int l ;
                            for(l=jj+1;l<nbSet ;l++) {
                                Set[l-1] = Set[l] ;
                            }
                            nbSet--; jj-- ;
                            Set[nbSet] = 0 ;
                        }
                    }
                    int l ; // add new set in good place
                    for(l=nbSet-1;l>=0 && Set[l] > newContact  ;l--) {
                        Set[l+1] = Set[l] ;
                    }
                    Set[l+1] = newContact ;
                    nbSet++ ;
                }
                globalContact |= newContact ;
            }
            if(globalContact != jc) { // add missing connected parts of new column
                globalContact ^= jc ;
                for(int k=0;k<PC[jc].nbSet;k++) {
                    if((PC[jc].Set[k] & globalContact) == PC[jc].Set[k] ){
                        int l ;
                        for(l=nbSet-1;l>=0 && Set[l] > PC[jc].Set[k]   ;l--) {
                            Set[l+1] = Set[l] ;
                        }
                        Set[l+1] = PC[jc].Set[k] ;
                        nbSet++ ;
                    }
                }
            }
            int in ; // find the resulting PP (if missing BIB BUG !!)
            for(in=data701->firstPPByMask[jc];in<data701->firstPPByMask[jc+1];in++) {
                if(memcmp(Set,data701->PP[in].Set,sizeof(Set))==0) break ;
            }
            assert(in != data701->firstPPByMask[jc+1]) ;
            data701->PT[ip*data701->nbCol+jc].numPP = in ;
            // compute the transformation
            // delta[] is the size of set (number ob black case to add)
            for(int id=0;id < data701->PP[in].nbSet ; id++) data701->PT[ip*data701->nbCol+jc].delta[id] =  nbBit[data701->PP[in].Set[id]] ;
            for(int is=0;is<data701->PP[ip].nbSet;is++) {
                if(data701->PP[ip].Set[is] & jc) { // numCv[] is the correspondance of index from old PP to new PP.
                    for(int id=0;id < data701->PP[in].nbSet ; id++){
                        if( data701->PP[ip].Set[is] & data701->PP[in].Set[id] ) {
                            data701->PT[ip*data701->nbCol+jc].numCv[is] = id+1 ;
                            break ;
                        }
                        assert(id < data701->PP[in].nbSet) ;
                    }
                }
            }
        }
    }
    free(PC) ;
    free(nbBit) ;
    free(interSec);
}

//***************************************************************
//  exact hash for decomposition N = n1+ n2 + n3 ... + nk.
// ni can be zero.

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
        default : {
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

