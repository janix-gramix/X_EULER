//
//  DC3.c
//  X_euler
//
//  Created by Jeannot on 17/12/2019.
//

#include "DC3.h"

void RadixSort(const DC3_INDEX * inInd, DC3_INDEX * outInd,const DC3_TYPE *T, int n,DC3_TYPE nbType )
{ // count occurrences
    DC3_INDEX * count = calloc(nbType+1,sizeof(count[0])) ;
    int i;
    for (i = 0; i < n ; i++)  count[T[inInd[i]]]++;
    DC3_INDEX sumc ;
    for (i = 0, sumc = 0; i <= nbType; i++)  {// exclusive prefix sums
        int32_t tmp = count [i]; count[i] = sumc ; sumc += tmp;
    }
    for ( i = 0; i < n; i++) outInd[count[T[inInd[i]]]++ ] = inInd[i]; // sort
    free(count) ;
    
}
typedef  int32_t T3[3] ;
typedef struct TT3 {
    int32_t V[3] ;
} TT3 ;

DC3_INDEX * CountRadixSort_uint16(const DC3_INDEX * inInd, DC3_INDEX * outInd,const uint16_t *T, int n,DC3_TYPE nbType )
{ // count occurrences
    DC3_INDEX * count = calloc(nbType+1,sizeof(count[0])) ;
    int i;
    for (i = 0; i < n ; i++)  count[T[inInd[i]]]++;
    DC3_INDEX sumc ;
    for (i = 0, sumc = 0; i <= nbType; i++)  {// exclusive prefix sums
        int32_t tmp = count [i]; count[i] = sumc ; sumc += tmp;
    }
    for ( i = 0; i < n; i++) outInd[count[T[inInd[i]]]++ ] = inInd[i]; // sort
    return count ;
    
}


void SuffixSort(const DC3_TYPE * T, DC3_INDEX * SA, int n, DC3_TYPE nbType) {
    DC3_INDEX n2=n/3, n0=(n+2)/3, d1 = 0 , n02=n0+n2;
    DC3_INDEX * R = malloc((n02+3)*sizeof(R[0])) ; R[n02]=R[n02+1]=R[n02+2]=0;
    DC3_INDEX * SA12 = malloc((n02+3)*sizeof(SA12[0])); SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
    DC3_INDEX * R0 = malloc((n0+1)*sizeof(R0[0])) ;
    DC3_INDEX * SA0 = malloc((n0+1)*sizeof(SA0[0]));
    //******* Step 0: Construct sample ********
    // generate positions of mod 1 and mod 2 suffixes
    int i ; // R contains n02 index
    for(i=0;i< n2;i++) {   R[2*i] = 3*i+1;    R[2*i+1] = 3*i+2; }
    if (n == 3*n2+1) { R[2*n2] = n ;  d1=1 ; } // +1 (size of T)
    else if (n == 3*n2+2) R[2*n2] = n-1 ;
    //---------- Sort sample suffixes
    // lsb radix sort the mod 1 and mod 2 triples in C = B1 u B2 (samples suffixes)
    switch(n-3*n2) {
        case 0 :
            RadixSort(R , SA12+2, T+2, n02-2, nbType);  SA12[1]=n-2 ;
            RadixSort(SA12+1, R+1 , T+1, n02-1, nbType); R[0] = n-1 ;
            RadixSort(R , SA12, T , n02, nbType);
            break ;
        case 1 :
            RadixSort(R , SA12+2, T+2, n02-2, nbType); SA12[1]=n-2 ;
            RadixSort(SA12+1, R+1 , T+1, n02-1, nbType);
            RadixSort(R+1 , SA12+1, T , n02-1, nbType);  SA12[0] = n ;
            break ;
        case 2:
            RadixSort(R , SA12+1, T+2, n02-1, nbType);
            RadixSort(SA12+1, R+1 , T+1, n02-1, nbType); R[0] = n-1 ;
            RadixSort(R , SA12, T , n02, nbType);
            break ;
    }
    //        RadixSort(R , SA12, T+2, n02, nbType); // +2 (size ofT) (with +1)=> +3
    //        RadixSort(SA12, R , T+1, n02, nbType);
    //        RadixSort(R , SA12, T , n02, nbType);
    
    // find different triples and write them in R
    DC3_INDEX nb3uple = 0;
    DC3_TYPE T3_0[3] = {-1,-1,-1};
    const DC3_TYPE *antT3 = T3_0 ;
    for (i=0 ; i < n02; i++) {
        int ind = SA12[i] ;
        if(ind < n-2) {
            const DC3_TYPE  * T3 = T+ind ;
            if(T3[0] != antT3[0] || T3[1] != antT3[1] || T3[2] != antT3[2]) {
                nb3uple++;   antT3 = T3 ; // new triple
            }
        } else {
            nb3uple++; antT3 = T3_0 ;
        }
        int ind_3 = ind/3 ; // recall, in SA12 ind = 1 or 2 mod[3]
        if (ind - 3*ind_3 == 1) { R[ind_3] = nb3uple; } // ind=1 mod[3]
        else { R[ind_3+n0] = nb3uple; } // ind=2 mod[3]
    }
    
    if (nb3uple < n02) { // triple not unique => recursion
        SuffixSort(R, SA12, n02, nb3uple);
        // compute reverse numerotation in R (add 1 as num in SA12 as +1)
        for (int i = 0; i < n02; i++) R[SA12[i]] = i + 1;
    } else {// all suffix in 12 are distinct, remove 1 to begin num T3 at 0
        for (int i = 0; i < n02; i++) SA12[R[i] - 1] = i;
    }
    // ---------- Sort nonsample suffixes
    // radix sort the mod[3]= 0 suffixes from SA12
    for (int i=0, j=0; i < n02; i++) if (SA12[i] < n0) R0[j++] = 3*SA12[i];
    RadixSort(R0, SA0, T, n0, nbType);
    //---------- Merge
    // merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p=0, t=d1, k=0; k < n; k++) {
        DC3_INDEX i_t, i_p = SA0[p]; // pos of current offset 0 suffix
        if((i_t = SA12[t] ) < n0) {
            DC3_INDEX i_3t = i_t * 3 + 1  ;
            //   Comparison (3t+1,R[t + n0]) and (p,R[p/3]
            if(T[i_3t] < T[i_p] || ( T[i_3t] == T[i_p] && R[i_t + n0] <= R[i_p/3] )) {
                SA[k] = i_3t ; t++;
                // add remaining  SA0 suffixes
                if (t == n02)  while(p < n0) SA[++k] = SA0[p++];
            } else {
                SA[k] = i_p; p++;
                // add remaining SA12 suffixes
                if (p == n0)  while(t < n02){ SA[++k] = ((i_t= SA12[t++]) < n0) ? (3*i_t+1) : (3*(i_t - n0)+2) ; }
            }
        } else {
            DC3_INDEX i_t3 =(i_t - n0) * 3 + 2 ;
            //  Comparison (3t+2,3t+3,R[t-n0+1]) and (p,p+1,R[p/3+n0]
            if(T[i_t3] < T[i_p]) {
                SA[k] = i_t3; t++;
                // add remaining  SA0 suffixes
                if (t == n02)  while(p < n0) SA[++k] = SA0[p++];
            } else if(T[i_t3] > T[i_p]) {
                SA[k] = i_p; p++;
                // add remaining SA12 suffixes
                if (p == n0)  while(t < n02){ SA[++k] = ((i_t= SA12[t++]) < n0) ? (3*i_t+1) : (3*(i_t - n0)+2) ; }
                
            } else {
                DC3_INDEX T31 = (i_t3 < n-1) ? T[i_t3+1] : -1 ;
                DC3_INDEX TP1 = (i_p < n-1) ? T[i_p+1] : - 1 ;
                if( T31 < TP1 || (T31==TP1 && R[i_t-n0+1] <=  R[i_p/3+n0]) ) {
                    SA[k] = i_t3; t++;
                    // add remaining  SA0 suffixes
                    if (t == n02)  while(p < n0) SA[++k] = SA0[p++];
                } else {
                    SA[k] = i_p; p++;
                    if (p == n0)  while(t < n02){ SA[++k] = ((i_t= SA12[t++]) < n0) ? (3*i_t+1) : (3*(i_t - n0)+2) ; }
                }
            }
        }
    }
    free(R); free(SA12); free(SA0); free(R0);
}

static int AtomLCP(const uint8_t * text,DC3_INDEX len,DC3_INDEX ia,DC3_INDEX ib) {
    if(ia<=ib) { int32_t tmp = ia ; ia= ib ; ib =tmp ;    }
    int l = 0 ;
    while(ia < len && text[ia++] == text[ib++]) l++ ;
    return l ;
}
DC3_INDEX * GetLCP(const uint8_t * text,const DC3_INDEX * Sindex, int nbIndex,int len) {
    DC3_INDEX * lcp = malloc(nbIndex*sizeof(lcp[0])) ;
    lcp[0] = 0;
    int32_t *Sorder=malloc(nbIndex * sizeof(Sorder[0])) ; // inverse permutation of index
    for (int i = 0; i < nbIndex; i++) Sorder[Sindex[i]] = i; //inverse permutation
    if (Sorder[0]) { // chaine initiale, partie commune avec la chaine totale
        lcp[Sorder[0]] = AtomLCP(text,len,0,Sindex[Sorder[0] - 1]);
    }
    for (int i = 1; i < nbIndex; i++) {
        if (!Sorder[i]) continue; // rank=0
        if (lcp[Sorder[i - 1]] <= 1) { // si la chaine precedente n'a rien de commun calcul direct
            lcp[Sorder[i]] = AtomLCP(text ,len, i, Sindex[Sorder[i] - 1]);
        } else { // si la chaine precedente a un prefixe commun avec celle d'avant on regarde si l'on peut prolonger
            int L = lcp[Sorder[i - 1]] - 1;
            lcp[Sorder[i]] = L + AtomLCP(text,len, i + L,  Sindex[Sorder[i] - 1] + L);
        }
    }
    free(Sorder);
    return lcp ;
}



