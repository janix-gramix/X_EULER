

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PB_MAX_STRLEN   100

typedef struct PB_RESULT {
    const char    *ident ;
    int     isVerbose ;
    char    strRes[PB_MAX_STRLEN] ;
    clock_t nbClock ;
} PB_RESULT ;

int PB151a(PB_RESULT *pbR) ;

int main(int argc, char**argv) {
    PB_RESULT pbR ;
    pbR.ident = "PB151a" ;
    pbR.isVerbose = 1 ;
    PB151a(&pbR);    
    fprintf(stdout,"\tPB%s(%.06fs) Sol=%s \n",pbR.ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes) ;
    return 0 ;
    
}

uint64_t PGCD64(uint64_t n1,uint64_t n2 ) {
    uint64_t r ;
    if (n1 > n2) {
        r = n2 ;
        n2 = n1 ;
    } else {
        r = n1 ;
    }
    while ( r > 0) {
        n1 = n2 ;
        n2 = r ;
        r = n1 % n2 ;
    }
    return n2 ;
}


typedef struct PB151a_STATE {
    int na[4] ; // number of sheet A2,A3,A4,A5
    int nbSheet ;  // total number of sheets
    uint64_t prob ;   // proba of the state
} PB151a_STATE ;

// number range for sheet
#define n2R   2
#define n3R   3
#define n4R   5
#define n5R   9

#define allR   (n2R*n3R*n4R*n5R)


int PB151a(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    PB151a_STATE st[allR] ;
    int stateToRg[allR] ;
    int rgLevel[16] ;
    uint64_t sum[16] ;
    memset(stateToRg,0,allR*sizeof(stateToRg[0]));
    int ns = 0 ; //
    rgLevel[0] = ns ; // initialize first state A2=A3=A4=A5=1
    sum[0] = 1 ;
    st[ns].na[0] =  st[ns].na[1] =  st[ns].na[2] = st[ns].na[3] = 1 ;
    st[ns].nbSheet = st[ns].na[0]+st[ns].na[1]+st[ns].na[2]+st[ns].na[3] ;
    st[ns].prob = 1 ;
    ns++ ;
    rgLevel[1] = ns ;
    int nl,i ;
    uint64_t alphaN=0 ; // proba for one sheet alphaN/alphaD
    uint64_t alphaD=1 ;
    for(nl=0;nl<14;nl++) { // loop far a week
        int no,nn ;
        int denFactor = 0 ; // bit mask for number of sheets 1,2,3,4,5,6,7,8
        for(no=rgLevel[nl];no<rgLevel[nl+1];no++) {
            denFactor |= 1 << st[no].nbSheet ; // add bit
        }
        uint64_t probDen = 1 ; // product of active number of sheets
        for(i=8;i>1;i--) { // so  prob contribution to Ax  => probDen / nbSheet * nb(Ax)
            if(denFactor & ( 1<< i)) probDen *= i ;
        }
        for(no=rgLevel[nl];no<rgLevel[nl+1];no++) { // developpement for active states
            uint64_t den = st[no].prob * probDen / st[no].nbSheet ; // common factor
            int index = st[no].na[3]+n5R*(st[no].na[2]+n4R*(st[no].na[1]+n3R*st[no].na[0]));
            int ia ;
            for(ia=0;ia<4;ia++) { // loop on Ax
                if(st[no].na[ia]==0) continue ;
                int idNext ; // compute next hash(exact) index
                if(ia==3)  idNext = index - 1 ;
                else if(ia==2) idNext = index - n5R +1 ;
                else if(ia==1) idNext = index - n5R*n4R + n5R +1 ;
                else idNext = index - n5R*n4R*n3R + n5R*n4R + n5R +1 ;
                int isn =stateToRg[idNext] ;
                if(isn==0) { // newstate
                    int ib ;
                    isn = ns++ ; stateToRg[idNext] = isn ;
                    memcpy(st[isn].na,st[no].na,4*sizeof(st[0].na[0])) ;
                    st[isn].na[ia]-- ; for(ib=ia+1;ib<4;ib++) st[isn].na[ib]++ ;
                    st[isn].nbSheet =st[no].nbSheet+2-ia ;
                    st[isn].prob =0 ;
                }
                st[isn].prob  +=  den *  st[no].na[ia] ; // add contribution
            }
        }
        rgLevel[nl+2] = ns ;
        uint64_t p = 0 ; // search pgcd of active probas
        for(nn=rgLevel[nl+1];nn<rgLevel[nl+2];nn++) {
            if (p==0) p = st[nn].prob ;
            else p = PGCD64(p,st[nn].prob);
            if(p==1) break ;
        }
        sum[nl+1] = 0 ;
        uint64_t sum1 = 0 ;
        for(nn=rgLevel[nl+1];nn<rgLevel[nl+2];nn++) {
            if(p>1)st[nn].prob /= p ; // remove pgcd to keep smaller probas
            //         printf("P(%d.%d.%d.%d:%d)=%lld ",st[nn].na[0],st[nn].na[1],st[nn].na[2],st[nn].na[3],st[nn].nbSheet,st[nn].prob );
            if(st[nn].nbSheet==1) sum1 += st[nn].prob ; // one sheet ?
            sum[nl+1] += st[nn].prob ;
        }
        if(nl<13 && sum1) { // add proba for one sheet
            uint64_t pgcdDS = PGCD64(alphaD,sum[nl+1]) ;
            alphaN = alphaN * (sum[nl+1]/pgcdDS) + sum1 * (alphaD/pgcdDS) ;
            alphaD *= (sum[nl+1]/pgcdDS) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t%s P=%lld/%lld=%.6f\n",pbR->ident,alphaN,alphaD, (double)alphaN/alphaD);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%.6f",(double)alphaN/alphaD) ;
    return 1 ;
}

