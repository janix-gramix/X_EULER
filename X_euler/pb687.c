//
//  pb687.c
//  
//
//  Created by Jeannot on 01/12/2019.
//  Copyright © 2019 Jeannot. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


static  uint64_t Cnp(int n,int p) {
    if(p==0) return 1 ;
    if(n<p) return 0 ;
    if(p==1) return n ;
    uint64_t cnp = n*(n-1)/2 ;
    if(p>2) {
        int i  ;
        for(i=2;i<p;i++) {
            cnp *= n-i ;
            cnp /= i+1 ;
        }
    }
    return cnp ;
}

static int Index5(int S,int i0,int i1,int i2,int i3) {
    static int Cn3[] = {0,1,5,15,35,70,126,210,330,495,715,1001,1365,1820} ;
    static int Cn2[] = {0,1,4,10,20,35,56,84,120,165,220,286,364,455} ;
    static int Cn1[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91} ;
    int ind = Cn3[S-i0]+Cn2[S-i0-i1]+Cn1[S-i0-i1-i2]+S-i0-i1-i2-i3 ;
    return ind ;
}



#define FLOAT687
#if defined(FLOAT687)
typedef  double PROB687 ;
#define SetIndex5(NR,N0,N1,N2,N3,ant,val)  cur5[Index5(NR,N0,N1,N2,N3)] += ant*(val)
#endif
#   define MAX_NR   13


// insertion ****
#define INS4_N0(ant) SetIndex5(nr,n0,n1,n2,n3,ant,nb0)
#define INS4_N1(ant) if(n1>0)SetIndex5(nr,n0+1,n1-1,n2,n3,ant,n1)
#define INS4_N2(ant) if(n2>0)SetIndex5(nr,n0,n1+1,n2-1,n3,ant,2*n2)
#define INS4_N3(ant) if(n3>0)SetIndex5(nr,n0,n1+1,n2,n3-1,ant,2*n3)
#define INS4_N4_2(ant) if(n4>0)SetIndex5(nr,n0,n1,n2+1,n3,ant,n4)
#define INS4_N4_3(ant) if(n4>0)SetIndex5(nr,n0,n1,n2,n3+1,ant,2*n4)

int main() {
    clock_t nbClock = clock()  ;
    int nr ;
    
    PROB687 buf0[2380],buf1[2380] ;
    PROB687 ant, *ant5 = buf0, *cur5=buf1 ;
    
    ant5[Index5(1,0,0,0,0)] = 1 ;
    ant5[Index5(1,1,0,0,0)] = 0 ;
    ant5[Index5(1,0,1,0,0)] = 0 ;
    ant5[Index5(1,0,0,1,0)] = 0 ;
    ant5[Index5(1,0,0,0,1)] = 0 ;
    
    int n0,n1,n2,n3;
    // n0 * * * *
    // n1 ** * *
    // n2 ** **
    // n3 *** *
    // n4 ****
    n1 = 0 ;
    
    for(nr=2;nr<=MAX_NR;nr++) {
        int i ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) cur5[i] = 0 ;
        // pour rajouter un rang : cela peut être : *,*,* || **,* || ***
        for(n0=0;n0<nr;n0++) {
            for(n1=0;n1<nr-n0;n1++) {
                for(n2=0;n2<nr-n0-n1;n2++) {
                    for(n3=0;n3<nr-n0-n1-n2;n3++) {
                        int n4=nr-1-n0-n1-n2-n3  ;
                        int nb0 = 4*n0+3*n1+2*n2+2*n3+n4 ; // not +1, by symetry on ranks we can impose not to insert at the end
                        ant = ant5[Index5(nr-1,n0,n1,n2,n3)] ;
                        if(ant <= 0) continue ;
                        // insert ****
                        INS4_N0(ant) ; INS4_N1(ant) ; INS4_N2(ant) ; INS4_N3(ant) ; INS4_N4_2(ant) ; INS4_N4_3(ant) ;
                        
                        int pow3[4] = { 1,3,9,27} ;
                        int i0,i1,i21,i22,i31,i32,i42,i43,nbp ;
                        int d,d0,d1,d2,d3 ;
                        // ** * *  || * * * *
                        for(d=0;d<4;d++) {
                            d0=d1=d2=d3=0 ;
                            switch(d) {
                                case 0 : d0=1; nbp = 4; break ; // * * * *
                                case 1 : d1=1; nbp = 3; break ; // ** * *
                                case 2 : d2=1; nbp = 2 ;break ; // ** **
                                case 3 : d3=1; nbp = 2 ;break ; // *** *
                            }
                            for(i0=0;i0<=nbp && i0 <= nb0 ;i0++){
                                int64_t delta0 =  Cnp(nb0,i0) ;
                                if(d1 ) delta0 *= 3 ;
                                else if(d3) delta0 *= 2 ;
                                for(i1=0; i1<=nbp-i0 && i1<=n1;i1++){
                                    int64_t delta1 = delta0 * Cnp(n1,i1) ;
                                    int i2max = nbp-i0-i1 ;
                                    if(i2max==0) { SetIndex5( nr,n0+i1+d0,n1-i1+d1,n2+d2,n3+d3,ant,delta1) ; break ;}
                                    for(i22=0; i22*2 <= i2max && i22<=n2;i22++){
                                        for(i21=0;i21<=(i2max-2*i22) && i21+i22 <= n2;i21++){
                                            int64_t delta2 = delta1 * Cnp(n2,i22) * Cnp(n2-i22,i21)*(1<<i21);
                                            int i3max = i2max-2*i22 - i21 ;
                                            if(i3max==0) {
                                                SetIndex5( nr,n0+d0+i1+i22,n1+d1-i1+i21,n2+d2-i22-i21,n3+d3,ant,delta2) ;
                                                break ;
                                            }
                                            for(i32=0; i32*2<=i3max && i32 <= n3;i32++){
                                                for(i31=0; i31 <= i3max-2*i32 && i31+i32 <=n3;i31++){
                                                    int64_t delta3 = delta2 * Cnp(n3,i32) * Cnp(n3-i32,i31)*(1<<i31) ;
                                                    int i4max = i3max -2 * i32 - i31 ;
                                                    if(i4max == 0 ) { SetIndex5( nr,n0+d0+i1+i22+i32,n1+d1-i1+i21+i31,n2+d2-i22-i21,n3+d3-i32-i31,ant, delta3) ;break ; }
                                                    for(i43=0;i43*3<= i4max && i43 <=n4 ;i43++){
                                                        for(i42=0;i42*2<=(i4max-3*i43) && i43+i42<=n4;i42++){
                                                            int i41 = i4max-3*i43 -2*i42 ;
                                                            if(i41+i42+i43 >n4) continue ;
                                                            int64_t delta = delta3 ;
                                                            if(i43) delta *= n4 ;
                                                            if(i42) delta *= Cnp(n4-i43,i42)*pow3[i42]  ;
                                                            int j3 ;
                                                            for(j3=0;j3<=i41;j3++) {
                                                                SetIndex5( nr,n0+d0+i1+i22+i32+i43 ,n1+d1-i1+i21+i31+i42  ,n2+d2-i22-i21+i41-j3 ,n3+d3-i32-i31+j3,
                                                                          ant , delta * Cnp(n4-i43-i42,i41) * (1<<j3) * Cnp(i41,j3));
                                                            }
                                                            delta *= Cnp(n4-i43-i42,i41)*pow3[i41]  ;
                                                        } // end i42
                                                    } // end i43
                                                } // end i31
                                            } // end i32
                                        } // end i21
                                    } // end i22
                                } // end i1
                            } // end i0
                        } // end nbp ** * *  or * * * *
                        
                    } // end n3
                } // end n2
            } // end n1
        } // end n0
        
        PROB687 den,prob[20];
        
        for(i=0;i<=nr;i++)prob[i] = 0 ;
        den = 0 ;
        for(i=0;i<=Index5(nr,0,0,0,0);i++) den += cur5[i] ;
        for(n0=0;n0<=nr;n0++) {
            for(n1=0;n1<=nr-n0;n1++) {
                for(n2=0;n2<=nr-n0-n1;n2++) {
                    for(n3=0;n3<=nr-n0-n1-n2;n3++) {
                        prob[n0]+= cur5[Index5(nr,n0,n1,n2,n3)] ;
                    }
                }
            }
        }
        
        fprintf(stdout,"\t 4 color, %d rank : ",nr);
        for(i=0;i<=nr;i++) { printf("%.11f ",prob[i]/(double)den) ; } ;
        printf("\n");
        if(nr==MAX_NR) {
            nbClock = clock() - nbClock ;
            printf("\tTime=%.6fs Prob(Prime)=%.10f\n",nbClock/(float)(CLOCKS_PER_SEC) ,((double)(prob[2]+prob[3]+prob[5]+prob[7]+prob[11]+prob[13]))/den);
            return 1 ;
        }
        PROB687 * tmp = ant5 ;
        ant5 = cur5 ;
        cur5 = tmp ;
    } // end loop nr
    printf("\n*********\n") ;
    return 0;
}
