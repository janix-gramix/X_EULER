//
//  PB001_050.c
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


#include "PB001_050.h"

int PB001(PB_RESULT *pbR) {
    int32_t i ;
    int64_t S = 0 ;
    int32_t i3 = 0;
    int32_t i5 = 0 ;
    pbR->nbClock = clock()  ;
    for(i=0;i<1000;i++) {
        if(i3 == 0  || i5 ==0) S += i ;
        if(++i3==3) i3 =0 ;
        if(++i5==5) i5 =0 ;
    };
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes), "%lld",S);
    return 1 ;
}

int PB002(PB_RESULT *pbR) {
    int64_t S = 0 ;
    int32_t u0 = 1 ;
    int32_t u1 = 2 ;
    pbR->nbClock = clock()  ;
    
    while ( u1 < 4000000) {
        if((u1 & 1) == 0) S += u1 ;
        u1 += u0 ;
        u0 = u1 -u0 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",S);
    return 1 ;
}

int PB003(PB_RESULT *pbR) {
    uint64_t n = 600851475143 ;
    uint64_t d = 2 ;
    uint64_t d2 = 2*2 ;
    uint64_t bp = 1 ;
    pbR->nbClock = clock() ;
    while ( d2 < n) {
        while ( (n % d) == 0) {
            bp = d ;
            n /= d ;
        }
        d++ ;
        d2 += d * 2 - 1 ;
    }
    if (n > bp) bp = n ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",bp);
    return 1 ;
    
}

int ChkPalyn(int n) {
    // on suppose que n > 100000
    if (   (((n % 100001) % 10010) % 1100) ==0  ){
        return 1 ;
    }
    return 0 ;
}

int PB004(PB_RESULT *pbR) {
    // on va assurer n1 >= n2
    // on balaye a somme constante (en partant du max S/2 * S/2)
    // on decremente la somme des que le produit est <  bestP
    // on termine quand le max S/2 * S/2 est < bestP
    int32_t n1 = 999 ;
    int32_t n2 = 999 ;
    int32_t P = n1*n2 ;
    int32_t S = n1+n2 ;
    int32_t bestP = 100000 ;
    int32_t bestn1 = 0 ;
    pbR->nbClock = clock() ;
    while ( 1 ) {
        //       printf("%d=%dx%d ",P,n1,n2);
        if ( (P > bestP) && ChkPalyn(P) ) {
            bestP = P ;
            bestn1 =n1 ;
        }
        if((P > bestP) && (n1 < 999) ) {
            P += n2 -n1 -1 ;
            n1++ ;
            n2-- ;
        } else {
            S-- ;
            n1 =(S+1)/2 ;
            n2 = S-n1 ;
            P = n1*n2 ;
            if(P <= bestP ) {
                break ;
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s P=%d = %d x %d \n",pbR->ident,bestP,bestn1,bestP/ bestn1);
    if ( bestn1 ) {
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",bestP);
        return 1 ;
    } else {
        return 0 ;
    }
}



int PB005(PB_RESULT *pbR) {
    int PPCM = 1;
    int i ;
    for(i=2;i<20;i++) {
        PPCM *= i/ PGCD(PPCM,i);
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",PPCM);
    return 1 ;
}

int PB006(PB_RESULT *pbR) {
    int64_t n = 100 ;
    int64_t Sn = (n *(n+1)) / 2 ;
    int64_t Sn2 = (n*(n+1)*(2*n+1)) / 6 ;
    int64_t diff = Sn * Sn - Sn2 ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",diff);
    return 1;
}

// #define PB007_NB_PRIME_ASK 100000000
// #define PB007_NB_MAX 2000000000

// #define PB007_NB_PRIME_ASK 100
// #define PB007_NB_MAX 542

// #define PB007_NB_PRIME_ASK 400001
#define PB007_NB_PRIME_ASK 10001
#define PB007_NB_MAX 2000000

#define PB007_GENTABLE  0
#define PB007_PRINTABLE 0




typedef struct CTX_PB007 {
    uint32_t nbPrime ;
    T_prime lastPrime ;
    uint32_t nbAsk ;
#if PB007_GENTABLE
    T_prime *tbPrime ;
#endif
    
} CTX_PB007 ;

int PB0007_nxtPrime(void *ptCtx,T_prime nxtPrime) {
    CTX_PB007 * ctx = (CTX_PB007 *) ptCtx ;
    ctx->lastPrime = nxtPrime ;
#if PB007_GENTABLE
    ctx->tbPrime[ctx->nbPrime] = nxtPrime ;
#endif
    ctx->nbPrime++ ;
    if(ctx->nbPrime >= ctx->nbAsk) {
        return 0 ;
    } else {
        return 1;
    }
}





int PB007(PB_RESULT *pbR) {
    CTX_PB007 ctx ;
    T_prime nbMax = PB007_NB_MAX ;
    ctx.nbPrime = 0 ;
    ctx.lastPrime = 0 ;
    ctx.nbAsk = PB007_NB_PRIME_ASK ;
    pbR->nbClock = clock() ;
#if PB007_GENTABLE
    ctx.tbPrime = malloc(ctx.nbAsk * sizeof(ctx.tbPrime[0])) ;
#endif
    FindPrime(nbMax,&ctx,PB0007_nxtPrime) ;
#if  PB007_PRINTABLE
    {
        int i;
        for(i=0;i<ctx.nbPrime ;i++) {
            printf("%d ",ctx.tbPrime[i]) ;
        }
    }
    printf("\n");
#endif
    if(pbR->isVerbose) {
        pbR->nbClock = clock() -  pbR->nbClock ;
        printf("\t PB%s(%.6fs)  prime n°%d  = %u\n",pbR->ident,(float)pbR->nbClock / CLOCKS_PER_SEC ,ctx.nbPrime,ctx.lastPrime) ;
    }
    ctx.nbPrime = 0 ;
    ctx.lastPrime = 0 ;
    pbR->nbClock = clock();
    FindPrime(nbMax,&ctx,PB0007_nxtPrime) ;
#if  PB007_PRINTABLE
    {
        int i;
        for(i=0;i<ctx.nbPrime ;i++) {
            printf("%d ",ctx.tbPrime[i]) ;
        }
    }
    printf("\n");
#endif
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%u",ctx.lastPrime);
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s  prime n°%d  = %u\n",pbR->ident,ctx.nbPrime,ctx.lastPrime) ;
    return 1 ;
    
    
}

#define PB008_NB    13
int PB008(PB_RESULT *pbR) {
    char strDig[] =
    "73167176531330624919225119674426574742355349194934"
    "96983520312774506326239578318016984801869478851843"
    "85861560789112949495459501737958331952853208805511"
    "12540698747158523863050715693290963295227443043557"
    "66896648950445244523161731856403098711121722383113"
    "62229893423380308135336276614282806444486645238749"
    "30358907296290491560440772390713810515859307960866"
    "70172427121883998797908792274921901699720888093776"
    "65727333001053367881220235421809751254540594752243"
    "52584907711670556013604839586446706324415722155397"
    "53697817977846174064955149290862569321978468622482"
    "83972241375657056057490261407972968652414535100474"
    "82166370484403199890008895243450658541227588666881"
    "16427171479924442928230863465674813919123162824586"
    "17866458359124566529476545682848912883142607690042"
    "24219022671055626321111109370544217506941658960408"
    "07198403850962455444362981230987879927244284909188"
    "84580156166097919133875499200524063689912560717606"
    "05886116467109405077541002256983155200055935729725"
    "71636269561882670428252483600823257530420752963450" ;
    int lg = (int) strlen(strDig) ;
    uint64_t bestProd = 0 ;
    int32_t bestIndex = 0  ;
    int32_t indNxt ;
    pbR->nbClock = clock() ;
    
    for(indNxt=PB008_NB;indNxt<lg;indNxt++) {
        int i ;
        uint64_t prod = 1 ;
        for(i=0;i<PB008_NB;i++) {
            prod *= ( strDig[indNxt -PB008_NB+i] - '0' ) ;
            if(prod == 0) {
                indNxt += i ;
                break ;
            }
        }
        if(prod > bestProd) {
            bestProd = prod ;
            bestIndex = indNxt - PB008_NB ;
        }
    }
    strDig[bestIndex+PB008_NB] = 0 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s Prod(%s)  =%lld\n",pbR->ident,strDig+bestIndex,bestProd) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",bestProd) ;
    return 1 ;
}

#define PB009_SUM   1000
int PB009(PB_RESULT *pbR) {
    // pour m > n > 0
    // a = m**2 - n**2 ; b = 2*m*n ; c = m**2 + n**2 ;
    // donc S = 2 * (m**2 + m*n) = 2 * m * (m+n)
    // 1000 = 2 * 25 * 20
    // m = 20 n=5
    int m = 20 ;
    int n = 5 ;
    int a = m*m - n*n ;
    int b = 2*m*n ;
    int c = m*m + n*n ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s %dx%dx%d=%d %d+%d+%d=%d %d**2+%d**2=%d %d**2=%d\n",pbR->ident, a,b,c,a*b*c,a,b,c,a+b+c,a,b,a*a+b*b,c,c*c) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",a*b*c) ;
    return 1;
}

typedef struct CTX_PB010 {
    uint64_t sum ;
} CTX_PB010 ;

#define PB010_NB_MAX    2000000
int PB010_nxtPrime(void *ptCtx,T_prime nxtPrime) {
    CTX_PB010 * ctx = (CTX_PB010 *) ptCtx ;
    if(nxtPrime > PB010_NB_MAX) return 0 ;
    ctx->sum += nxtPrime ;
    return 1;
}

int PB010(PB_RESULT *pbR) {    CTX_PB010 ctx ;
    T_prime nbMax = PB010_NB_MAX ;
    ctx.sum = 0 ;
    pbR->nbClock = clock() ;
    FindPrime(nbMax,&ctx,PB010_nxtPrime) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s SUM(prime<%u)  = %lld\n",pbR->ident,nbMax,ctx.sum) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",ctx.sum) ;
    return 1 ;
}

#define PB11_SIZE   20
int PB011(PB_RESULT *pbR) {
    static int tbIn[PB11_SIZE*PB11_SIZE] = {
        8,02,22,97,38,15,00,40,00,75,04,05,07,78,52,12,50,77,91, 8,
        49,49,99,40,17,81,18,57,60,87,17,40,98,43,69,48,04,56,62,00,
        81,49,31,73,55,79,14,29,93,71,40,67,53,88,30,03,49,13,36,65,
        52,70,95,23,04,60,11,42,69,24,68,56,01,32,56,71,37,02,36,91,
        22,31,16,71,51,67,63,89,41,92,36,54,22,40,40,28,66,33,13,80,
        24,47,32,60,99,03,45,02,44,75,33,53,78,36,84,20,35,17,12,50,
        32,98,81,28,64,23,67,10,26,38,40,67,59,54,70,66,18,38,64,70,
        67,26,20,68,02,62,12,20,95,63,94,39,63, 8,40,91,66,49,94,21,
        24,55,58,05,66,73,99,26,97,17,78,78,96,83,14,88,34,89,63,72,
        21,36,23, 9,75,00,76,44,20,45,35,14,00,61,33,97,34,31,33,95,
        78,17,53,28,22,75,31,67,15,94, 3,80,04,62,16,14, 9,53,56,92,
        16,39,05,42,96,35,31,47,55,58,88,24,00,17,54,24,36,29,85,57,
        86,56,00,48,35,71,89,07,05,44,44,37,44,60,21,58,51,54,17,58,
        19,80,81,68,05,94,47,69,28,73,92,13,86,52,17,77,04,89,55,40,
        04,52, 8,83,97,35,99,16,07,97,57,32,16,26,26,79,33,27,98,66,
        88,36,68,87,57,62,20,72,03,46,33,67,46,55,12,32,63,93,53,69,
        04,42,16,73,38,25,39,11,24,94,72,18, 8,46,29,32,40,62,76,36,
        20,69,36,41,72,30,23,88,34,62,99,69,82,67,59,85,74,04,36,16,
        20,73,35,29,78,31,90,01,74,31,49,71,48,86,81,16,23,57,05,54,
        01,70,54,71,83,51,54,69,16,92,33,48,61,43,52,01,89,19,67,48 } ;
    int ic,ir ;
    int bestProd = 0 ;
    int bestInd = 0 ;
    int bestDelta = 0;
    int prod ;
    pbR->nbClock = clock()  ;
    for(ir=0;ir<PB11_SIZE;ir++) {
        for(ic=0;ic<PB11_SIZE;ic++) {
            int ind = ir * PB11_SIZE + ic ;
            if(ic <= PB11_SIZE - 4) {
                prod = tbIn[ind] * tbIn[ind+1] * tbIn[ind+2] * tbIn[ind+3] ; // RIGHT
                if(prod > bestProd) {
                    bestProd = prod ;
                    bestInd = ind ;
                    bestDelta = 1;
                }
                if(ir <=  PB11_SIZE - 4) {
                    prod = tbIn[ind] * tbIn[ind+1+PB11_SIZE] * tbIn[ind+2+2*PB11_SIZE] * tbIn[ind+3+3*PB11_SIZE] ; // DIAG down-right
                    if(prod > bestProd) {
                        bestProd = prod ;
                        bestInd = ind ;
                        bestDelta = 1+PB11_SIZE;
                    }
                }
            }
            if(ir <=  PB11_SIZE - 4) {
                prod = tbIn[ind] * tbIn[ind+PB11_SIZE] * tbIn[ind+2*PB11_SIZE] * tbIn[ind+3*PB11_SIZE] ; // down
                if(prod > bestProd) {
                    bestProd = prod ;
                    bestInd = ind ;
                    bestDelta = PB11_SIZE;
                }
                if(ir >=  3) {
                    prod = tbIn[ind] * tbIn[ind-1+PB11_SIZE] * tbIn[ind-2+2*PB11_SIZE] * tbIn[ind-3+3*PB11_SIZE] ; // DIAG down-left
                    if(prod > bestProd) {
                        bestProd = prod ;
                        bestInd = ind ;
                        bestDelta = PB11_SIZE-1;
                    }
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s  [%d,%d]x...x[%d,%d]=%dx%dx%dx%d  = %d\n"
                              ,pbR->ident
                              , bestInd / PB11_SIZE ,bestInd % PB11_SIZE
                              , (bestInd+3*bestDelta) / PB11_SIZE ,(bestInd+3*bestDelta) % PB11_SIZE
                              , tbIn[bestInd],tbIn[bestInd+bestDelta],tbIn[bestInd+2*bestDelta],tbIn[bestInd+3*bestDelta]
                              ,bestProd ) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",bestProd);
    return 1 ;
}



#define PB012_MAX_SQPRIME   1000000 // permet de factoriser jusqu'a 10**12
int PB012(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB012_MAX_SQPRIME)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    } else {
        int64_t n = 0 ;
        int64_t ntr;
        int64_t nbDiv ;
        const T_prime * tbPrime = GetTbPrime(ctxP) ;
        do {
            n++ ;
            ntr =(n*(n+1))/2 ;
            nbDiv=FindNbDiv( ntr,tbPrime) ;
            
        } while (nbDiv <= 500) ;
        pbR->nbClock = clock() - pbR->nbClock ;
        if(pbR->isVerbose)fprintf(stdout,"\t PB%s triangle(%lld)=%lld a %lld diviseurs\n",pbR->ident,n,ntr,nbDiv);
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",ntr);
    }
    Free_tablePrime(ctxP);
    return 1;
}

#define PB013_NBNUM     100
#define PB013_LENNUM    50
int PB013(PB_RESULT *pbR) {
    static char strNumber[PB013_NBNUM*PB013_LENNUM] =
    "37107287533902102798797998220837590246510135740250"
    "46376937677490009712648124896970078050417018260538"
    "74324986199524741059474233309513058123726617309629"
    "91942213363574161572522430563301811072406154908250"
    "23067588207539346171171980310421047513778063246676"
    "89261670696623633820136378418383684178734361726757"
    "28112879812849979408065481931592621691275889832738"
    "44274228917432520321923589422876796487670272189318"
    "47451445736001306439091167216856844588711603153276"
    "70386486105843025439939619828917593665686757934951"
    "62176457141856560629502157223196586755079324193331"
    "64906352462741904929101432445813822663347944758178"
    "92575867718337217661963751590579239728245598838407"
    "58203565325359399008402633568948830189458628227828"
    "80181199384826282014278194139940567587151170094390"
    "35398664372827112653829987240784473053190104293586"
    "86515506006295864861532075273371959191420517255829"
    "71693888707715466499115593487603532921714970056938"
    "54370070576826684624621495650076471787294438377604"
    "53282654108756828443191190634694037855217779295145"
    "36123272525000296071075082563815656710885258350721"
    "45876576172410976447339110607218265236877223636045"
    "17423706905851860660448207621209813287860733969412"
    "81142660418086830619328460811191061556940512689692"
    "51934325451728388641918047049293215058642563049483"
    "62467221648435076201727918039944693004732956340691"
    "15732444386908125794514089057706229429197107928209"
    "55037687525678773091862540744969844508330393682126"
    "18336384825330154686196124348767681297534375946515"
    "80386287592878490201521685554828717201219257766954"
    "78182833757993103614740356856449095527097864797581"
    "16726320100436897842553539920931837441497806860984"
    "48403098129077791799088218795327364475675590848030"
    "87086987551392711854517078544161852424320693150332"
    "59959406895756536782107074926966537676326235447210"
    "69793950679652694742597709739166693763042633987085"
    "41052684708299085211399427365734116182760315001271"
    "65378607361501080857009149939512557028198746004375"
    "35829035317434717326932123578154982629742552737307"
    "94953759765105305946966067683156574377167401875275"
    "88902802571733229619176668713819931811048770190271"
    "25267680276078003013678680992525463401061632866526"
    "36270218540497705585629946580636237993140746255962"
    "24074486908231174977792365466257246923322810917141"
    "91430288197103288597806669760892938638285025333403"
    "34413065578016127815921815005561868836468420090470"
    "23053081172816430487623791969842487255036638784583"
    "11487696932154902810424020138335124462181441773470"
    "63783299490636259666498587618221225225512486764533"
    "67720186971698544312419572409913959008952310058822"
    "95548255300263520781532296796249481641953868218774"
    "76085327132285723110424803456124867697064507995236"
    "37774242535411291684276865538926205024910326572967"
    "23701913275725675285653248258265463092207058596522"
    "29798860272258331913126375147341994889534765745501"
    "18495701454879288984856827726077713721403798879715"
    "38298203783031473527721580348144513491373226651381"
    "34829543829199918180278916522431027392251122869539"
    "40957953066405232632538044100059654939159879593635"
    "29746152185502371307642255121183693803580388584903"
    "41698116222072977186158236678424689157993532961922"
    "62467957194401269043877107275048102390895523597457"
    "23189706772547915061505504953922979530901129967519"
    "86188088225875314529584099251203829009407770775672"
    "11306739708304724483816533873502340845647058077308"
    "82959174767140363198008187129011875491310547126581"
    "97623331044818386269515456334926366572897563400500"
    "42846280183517070527831839425882145521227251250327"
    "55121603546981200581762165212827652751691296897789"
    "32238195734329339946437501907836945765883352399886"
    "75506164965184775180738168837861091527357929701337"
    "62177842752192623401942399639168044983993173312731"
    "32924185707147349566916674687634660915035914677504"
    "99518671430235219628894890102423325116913619626622"
    "73267460800591547471830798392868535206946944540724"
    "76841822524674417161514036427982273348055556214818"
    "97142617910342598647204516893989422179826088076852"
    "87783646182799346313767754307809363333018982642090"
    "10848802521674670883215120185883543223812876952786"
    "71329612474782464538636993009049310363619763878039"
    "62184073572399794223406235393808339651327408011116"
    "66627891981488087797941876876144230030984490851411"
    "60661826293682836764744779239180335110989069790714"
    "85786944089552990653640447425576083659976645795096"
    "66024396409905389607120198219976047599490197230297"
    "64913982680032973156037120041377903785566085089252"
    "16730939319872750275468906903707539413042652315011"
    "94809377245048795150954100921645863754710598436791"
    "78639167021187492431995700641917969777599028300699"
    "15368713711936614952811305876380278410754449733078"
    "40789923115535562561142322423255033685442488917353"
    "44889911501440648020369068063960672322193204149535"
    "41503128880339536053299340368006977710650566631954"
    "81234880673210146739058568557934581403627822703280"
    "82616570773948327592232845941706525094512325230608"
    "22918802058777319719839450180888072429661980811197"
    "77158542502016545090413245809786882778948721859617"
    "72107838435069186155435662884062257473692284509516"
    "20849603980134001723930671666823555245252804609722"
    "53503534226472524250874054075591789781264330331690";
    uint8_t digitNumber[PB013_NBNUM*PB013_LENNUM] ;
    uint32_t sum[PB013_LENNUM] ;
    int i ;
    pbR->nbClock = clock()  ;
    for(i=0;i<PB013_NBNUM*PB013_LENNUM;i++) {
        digitNumber[i] = strNumber[i] - '0' ;
    }
    for(i=0;i<PB013_LENNUM;i++){
        int j ;
        sum[i] = 0 ;
        for(j=0;j<PB013_NBNUM;j++) {
            sum[i] += digitNumber[i+j*PB013_LENNUM] ;
        }
    }
    for(i=PB013_LENNUM-1;i>0;i--) {
        sum[i-1] += sum[i]/10 ;
        sum[i] = sum[i] % 10 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose){
        fprintf(stdout,"\t bP%s  Sum=",pbR->ident);
        for(i=0;i<PB013_LENNUM;i++) fprintf(stdout,"%d",sum[i]);
        fprintf(stdout,"\n");
    }
    {
        int lg  = 0 ;
        for(i=0;lg < 10;i++) {
            lg+= snprintf(pbR->strRes+lg,sizeof(pbR->strRes)-lg,"%d",sum[i]);
        }
        pbR->strRes[10] = '\0' ;
    }
    return 1 ;
}

typedef uint32_t TY_SYR   ;


void SyracuseVal(TY_SYR *valCur) {
    TY_SYR val = valCur[0] ;
    if(val == 1) {
        return ;
    } else {
        if(val & 1) {
            valCur[1] = (3 * val + 1) ;
        }else {
            valCur[1] = val / 2 ;
        }
        SyracuseVal(valCur+1);
    }
}

// ne retourne que la longueur, pas les valeurs
#define SYRACUSE_MEM    1
#if SYRACUSE_MEM
#define SYRACUSE_MEM_SIZE (1<<19)
static uint16_t lgS[SYRACUSE_MEM_SIZE] ;
static uint16_t privSyracuseLgMem(uint16_t lg,TY_SYR valInit,TY_SYR val) {
    if(val < SYRACUSE_MEM_SIZE ) {
        if(lgS[val]) {
            if(valInit < SYRACUSE_MEM_SIZE) lgS[valInit] = lg+lgS[val] ;
            return lg+lgS[val] ;
        } else if (val == 1) {
            lgS[1] = 1 ;
            if(valInit < SYRACUSE_MEM_SIZE) lgS[valInit] = lg+lgS[val] ;
            return lg+lgS[val] ;
        }
    }
    if(val == 1) {
        return lg;
    } else
    {
        if(val & 1) {
            val = (3 * val + 1) ;
        }else {
            val = val / 2 ;
        }
        return privSyracuseLgMem(lg+1,valInit,val);
    }
}
#else
static uint16_t privSyracuseLg(uint16_t lg,TY_SYR val) {
    if(val == 1) {
        return lg+1;
    } else {
        if(val & 1) {
            val = (3 * val + 1) ;
        }else {
            val = val / 2 ;
        }
        return privSyracuseLg(lg+1,val);
    }
}
#endif

static uint16_t SyracuseLg(TY_SYR val) {
#if SYRACUSE_MEM
    return privSyracuseLgMem(0,val,val) ;
#else
    return privSyracuseLg(0,val) ;
#endif
}

#define PB014_MAX_VALUE 1000000
#define PB014_PRINT 0
int PB014(PB_RESULT *pbR) {
    TY_SYR k ;
    TY_SYR kBest = 1 ;
    int bestLg = 0 ;
    pbR->nbClock = clock() ;
    for(k=1;k<PB014_MAX_VALUE;k++) {
        uint16_t lg = SyracuseLg(k) ;
        if(lg > bestLg ) {
            bestLg = lg ;
            kBest = k ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s  lg=%d for %d\n"
                              ,pbR->ident
                              ,bestLg,kBest);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",kBest);
#if PB014_PRINT
    {
        TY_SYR *chainVal = malloc((bestLg+1)*sizeof(chainVal[0])) ;
        chainVal[0] = kBest ;
        SyracuseVal(chainVal);
        
        {
            int i,lg ;
            for(i=0,lg=0 ; chainVal[i] != 1 ; i++) {
                printf("%d->",chainVal[i] );
                lg++ ;
            }
            printf("%d\n",chainVal[i]);
        }
        free(chainVal);
    }
#endif
    return 1 ;
}

int PB015(PB_RESULT *pbR) {
    // answer = C(20,40) = 40! / ( 20! * 20!) = 21x22x...x40 / 1x2x3...x20
    // on apparie les indices pairs en haut et les 10 derniers en bas (haut =2xbas)
    //    = 21x23x25x27x...x39x 2**10 / 1x2x3x4...x10
    // On fait les simplifications 25x2 = 5x10  ; 27=3x9 ; 21x2 = 6x7 ; 2**6=2x4x8 ;
    // Il reste en haut 23x29x31x33x35x37x39x2**2
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld", ((uint64_t) 23) * 29 * 31 * 33 * 35 *37 * 39 * 4 );
    return 1 ;
}



#define PB016_MAXL  1000/3
#define PB016_EXP   1000
int PB016(PB_RESULT *pbR) {
    uint8_t digLarge[PB016_MAXL] ;
    int i,ie ,len = 0 ;
    uint32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB016_MAXL;i++) {digLarge[i] = 0 ; }
    digLarge[0] = 1 ; len = 1 ;
    for(ie=0;ie<PB016_EXP;ie++) {
        for(i=0;i<len;i++) {
            digLarge[i] *= 2 ;
        }
        for(i=0;i<len;i++) {
            if(digLarge[i] >= 10) {
                digLarge[i] -= 10 ;
                digLarge[i+1]++ ;
            }
        }
        if(digLarge[len]) len++ ;
    }
    for(i=0;i<len;i++){
        S += digLarge[i] ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s  Sumdig(2**%d)=%d\n"
                               ,pbR->ident
                               ,PB016_EXP,S) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S);
    return 1 ;
}

int PB017(PB_RESULT *pbR) {
    char * digits[9] = { "one" , "two", "three" , "four", "five", "six", "seven" , "eight", "nine" } ;
    char * d11_19[9] = { "eleven", "twelve" , "thirteen" , "fourteen" , "fifteen" , "sixteen" , "seventeen" , "eighteen" , "nineteen" } ;
    char * d10[8] = { "twenty", "thirty", "forty", "fifty", "sixty", "seventy" , "eighty" , "ninety" } ;
    pbR->nbClock = clock()  ;
    int len,lenDigits,lenD11_19, lenD10 ;
    int i ;
    for(i=0,lenDigits=0;i<9;i++) lenDigits += strlen(digits[i]) ;
    for(i=0,lenD11_19=0;i<9;i++) lenD11_19 += strlen(d11_19[i]) ;
    for(i=0,lenD10=0;i<8;i++) lenD10 += strlen(d10[i]) ;
    
    len =
        lenDigits // 1..9
        + strlen("ten") // 10
        + lenD11_19 // 11..19
        + lenD10 * 10 + 8 * lenDigits // 20 21 ..29 30 31 .. 99
        + 100 * lenDigits + 900 * strlen("hundred") + 891 * strlen("and")
        + 9 * ( lenDigits + strlen("ten") + lenD11_19 + lenD10 * 10 + 8 * lenDigits  )
        + strlen("one" ) + strlen("thousand")
    ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s Len=%d\n",pbR->ident,len);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",len);
     return 1;
}


#define PB018_SIZE  15
#define PB018_TSIZE ((PB018_SIZE*(PB018_SIZE+1))/2)
int PB018(PB_RESULT *pbR) {
    int vals[PB018_TSIZE] = {
        75
        ,95 ,64
        ,17 ,47 ,82
        ,18 ,35 ,87 ,10
        ,20 ,04 ,82 ,47 ,65
        ,19 ,01 ,23 ,75 ,03 ,34
        ,88 ,02 ,77 ,73 ,07 ,63 ,67
        ,99 ,65 ,04 ,28 ,06 ,16 ,70 ,92
        ,41 ,41 ,26 ,56 ,83 ,40 ,80 ,70 ,33
        ,41 ,48 ,72 ,33 ,47 ,32 ,37 ,16 ,94 ,29
        ,53 ,71 ,44 ,65 ,25 ,43 ,91 ,52 ,97 ,51 ,14
        ,70 ,11 ,33 ,28 ,77 ,73 ,17 ,78 ,39 ,68 ,17 ,57
        ,91 ,71 ,52 ,38 ,17 ,14 ,91 ,43 ,58 ,50 ,27 ,29 ,48
        ,63 ,66 ,04 ,68 ,89 ,53 ,67 ,30 ,73 ,16 ,69 ,87 ,40 ,31
        ,04 ,62 ,98 ,27 ,23 ,9 ,70 ,98 ,73 ,93 ,38 ,53 ,60 ,04 ,23 } ;
    int ic,ir ;
    pbR->nbClock = clock() ;
    // on commernce a l'avant derniere ligne
    for(ir=PB018_SIZE-2;ir>=0;ir--) {
        int ic0 = (ir*(ir+1))/2 ;
        for(ic=0;ic<=ir;ic++) {
            int icnr = ic0+ir+1+ic ;
            vals[ic0+ic] += (vals[icnr] > vals[icnr+1]) ? vals[icnr] : vals[icnr+1] ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",vals[0]) ;
    return 1 ;
}

int PB019(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    int nbSunday = 0 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    int iy ;
    int d0 = (1 + 365) % 7 ; // monday + 365 car 1900 non bissextile
    for(iy=1901;iy<=2000;iy++) {
        if(d0 == 0) nbSunday++ ; // 1 janvier
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 fevrier
        d0 = (d0+28 +  (  ((iy & 3) == 0) ? 1 : 0 ) ) % 7  ;  if(d0 == 0) nbSunday++ ; // 1 mars
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 avril
        d0 = (d0+30) % 7 ;  if(d0 == 0) nbSunday++ ;
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 juin
        d0 = (d0+30) % 7 ;  if(d0 == 0) nbSunday++ ;
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 aout
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ;
        d0 = (d0+30) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 octobre
        d0 = (d0+31) % 7 ;  if(d0 == 0) nbSunday++ ;
        d0 = (d0+30) % 7 ;  if(d0 == 0) nbSunday++ ; // 1 decembre
        d0 = (d0+31) % 7 ; 
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s Number of sundays=%d\n",pbR->ident,nbSunday);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbSunday);
    return 1;
}


#define PB020_MAXL  (100*2)
#define PB020_FACT   100
int PB020(PB_RESULT *pbR) {
    uint32_t digLarge[PB020_MAXL] ;
    int i,ie ,len = 0 ;
    uint32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB020_MAXL;i++) {digLarge[i] = 0 ; }
    digLarge[0] = 1 ; len = 1 ;
    for(ie=1;ie<=PB020_FACT;ie++) {
        for(i=0;i<len;i++) {
            digLarge[i] *= ie ;
        }
        for(i=0;i<len;i++) {
            if(digLarge[i] >= 10) {
                digLarge[i+1] += digLarge[i] /10 ;
                digLarge[i] = digLarge[i] % 10 ;
            }
        }
        while(digLarge[len]) {
            if(digLarge[len] >= 10) {
                digLarge[len+1] = digLarge[len] /10 ;
                digLarge[len] = digLarge[len] % 10 ;
            }
            len++ ;
        }
        
    }
    for(i=0;i<len;i++){
        S += digLarge[i] ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s  Sumdig(%d!)=%d\n",pbR->ident, PB020_FACT,S) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}


#define PB021_MAXM  10000
int PB021(PB_RESULT *pbR) {
    int SumDiv[PB021_MAXM+1] ;
    int i ;
    uint32_t S = 0 ;
    pbR->nbClock = clock() ;
    for(i=1;i<=PB021_MAXM;i++) {
        SumDiv[i] = 0 ;
    }
    for(i=1;2*i<PB021_MAXM;i++) {
        int mxi ;
        for(mxi=2*i;mxi <= PB021_MAXM ; mxi += i) {
            SumDiv[mxi] += i ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    for(i=1;i<PB021_MAXM;i++){
        if((SumDiv[i] > i) && (SumDiv[i] <= PB021_MAXM)) {
            if(SumDiv[SumDiv[i]] == i ) {
                fprintf(stdout,"%d<->%d ",i,SumDiv[i]) ;
                S += i + SumDiv[i] ;
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",S) ;
    return 1 ;
}

#include "p022_data.h"

int CmpName022(const void *e1, const void *e2) {
    static const char ** ptNames =  NULL;
    if(ptNames== NULL) ptNames =   P022_GetData()  ;
    int i1 = ((int *)e1)[0] ;
    int i2 = ((int *)e2)[0] ;
    return strcmp(ptNames[i1], ptNames[i2] );
}

int PB022(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    int nameOrder[p022_size] ;
    int i;
    for(i=0;i<p022_size;i++) {
        nameOrder[i] = i ;
    }
    qsort(nameOrder,p022_size,sizeof(nameOrder[0]),CmpName022);
    int sum = 0 ;
    char * * ptNames =   P022_GetData()  ;
    for(i=0;i<p022_size;i++) {
        char * c = ptNames[nameOrder[i]] ;
        int cost = 0 ;
        while(*c) {
            cost += *c++ - '@' ;
        }
        if(i==937) printf("%s(%d) ",ptNames[nameOrder[i]],cost) ;
        sum += (i+1) * cost ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s Sum=%d\n",pbR->ident,sum) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",sum) ;
    return 1 ;
}


#define PB023_MAXM  28123
int PB023(PB_RESULT *pbR) {
    int SumDiv[PB023_MAXM+1] ;
    int *isSumAbun ;
    int abundant[PB023_MAXM] ;
    int nbAbund =0 ;
    int i ;
    uint32_t SumNo2a = 0 ;
    int nbNoSum2a = 0 ;
    pbR->nbClock = clock()  ;
    for(i=1;i<=PB023_MAXM;i++) {
        SumDiv[i] = 0 ;
    }
    for(i=1;2*i<=PB023_MAXM;i++) {
        int mxi ;
        for(mxi=2*i;mxi <= PB023_MAXM ; mxi += i) {
            SumDiv[mxi] += i ;
        }
    }
    isSumAbun = SumDiv ;
    for(i=1;i<=PB023_MAXM;i++){
        if(SumDiv[i] > i) {
            abundant[nbAbund++] = i ;
        }
        isSumAbun[i] = 0 ;
    }
    for(i=0;i<nbAbund;i++) {
        int j ;
        for(j=i;j<nbAbund;j++) {
            int S = abundant[i]+abundant[j] ;
            if(S<= PB023_MAXM) {
                isSumAbun[S] = 1;
            } else {
                break ;
            }
        }
    }
    for(i=1;i<=PB023_MAXM;i++){
        if(isSumAbun[i] == 0) {
            SumNo2a += i ;
            nbNoSum2a++ ;
        }
    }
    
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s  nbNo2ab=%d ,Sum=%d\n",pbR->ident,nbNoSum2a,SumNo2a) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",SumNo2a);
    return 1 ;
}


#define PB024_SIZE  10
#define PB024_RANG   1000000
int PB024(PB_RESULT *pbR) {
    int i  ;
    uint64_t factoriels[PB024_SIZE] ;
    uint8_t digits[PB024_SIZE] ;
    int rang = PB024_RANG - 1 ;
    pbR->nbClock = clock() ;
    factoriels[0] = 1 ;
    for(i=1;i<PB024_SIZE;i++) {
        factoriels[i] = factoriels[i-1] *i ;
    }
    for(i=0;i<PB024_SIZE;i++) { digits[i] = i ; }
    for(i=0;i<PB024_SIZE;i++) {
        int rgUnused = rang / factoriels[PB024_SIZE-1-i] ;
        rang -= rgUnused * factoriels[PB024_SIZE-1-i] ;
        if(rgUnused>0){
            int j ;
            uint8_t saveDig = digits[i+rgUnused] ;
            for(j=rgUnused+i; j > i ; j--) {
                digits[j] = digits[j-1] ;
            }
            digits[i] = saveDig ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    for(i=0;i<PB024_SIZE;i++)pbR->strRes[i] = '0' + digits[i];
    pbR->strRes[10] = 0 ;
    return 1 ;
    
}

#define PB025_MAXDIGIT 1000
int PB025(PB_RESULT *pbR) {
    uint8_t FA[PB025_MAXDIGIT+1] ;
    uint8_t FB[PB025_MAXDIGIT+1] ;
    int i,n,len ;
    pbR->nbClock = clock() ;
    for(i=0;i<PB025_MAXDIGIT+1;i++) {
        FA[i]=FB[i]= 0 ;
    }
    FA[0] = FB[0] = 1 ;
    len = 1 ;
    n=2 ;
    do {
        uint8_t *FD = FA ;
        if(++n & 1) FD =FB ;
        for(i=0;i<len;i++) {
            if(FA[i]+FB[i]< 10) {
                FD[i] = FA[i]+FB[i] ;
            } else {
                FD[i] = FA[i]+FB[i]-10 ;
                FD[i+1]++ ;
            }
        }
        if(FD[len]) len++ ;
    } while(len < PB025_MAXDIGIT) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s len>=%d pour Fibonnaci n°%d \n",pbR->ident , PB025_MAXDIGIT,n) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",n);
    return 1 ;
}

#define PB026_MAX 1000
int PB026(PB_RESULT *pbR) {
    int maxCycle = 0 ;
    int dmax = 0 ;
    int d ;
    uint16_t isReste[PB026_MAX] ;
    pbR->nbClock = clock() ;
    for(d=2;d<PB026_MAX;d++) {
        int r = 1;
        int lc = 1 ;
        memset(isReste,0,d*sizeof(isReste[0]));
        isReste[1] = 1 ;
        while(1) {
            r *= 10 ; lc++ ;
            if(r >= d) {
                r = r % d ;
                if(r==0) break ;
                if( isReste[r]) {
                    if(lc - isReste[r] > maxCycle) {
                        maxCycle = lc - isReste[r] ;
                        dmax = d ;
                    }
                    break ;
                } else {
                    isReste[r] = lc ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s Pour 1/%d cycle de %d \n",pbR->ident, dmax,maxCycle) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",dmax) ;
    return 1 ;
}

#define PB027_MAX_A 1000 // doit etre pair
#define PB027_MAX_B 1000


int PB027(PB_RESULT *pbR) {
    // b est positif et est un nbre premier (P(0)=B)
    // a est impair (P(1)-P(0) =  1 + a
    // P(b) n'est pas premier donc n<b
    // si P(b-a)=(b-a)(b-a)+a)+b n'est pas premier. Donc si a<b, n<b-a
    // n<b => P(n) < (b**2 + b*a+b < Max_b**2 + Max_b*Max_a
    pbR->nbClock = clock()  ;
    CTX_PRIMETABLE * ctxP = Gen_tablePrime(PB027_MAX_B*(PB027_MAX_A+PB027_MAX_B)) ;
    if(ctxP == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    {
        int32_t ib = 0;
        int32_t maxNbprime = 0 ;
        int32_t  a , b ;
        int32_t aMax = 0 , bMax = 0 ;
        const T_prime * tbPrime = GetTbPrime(ctxP) ;
        while ( (b=(int32_t) tbPrime[ib++]) < PB027_MAX_B) {
            for ( a= (- PB027_MAX_A+1 ); a < PB027_MAX_A ;  a+=2) {
                int nbPrime = 0 ;
                int32_t n = 0 ;
                while(1) {
                    int32_t pn = n*(n+a)+b ;
                    if(pn<=0) { break ;}
                    if(! Search_TablePrime(ctxP, pn) ) break ;
                    nbPrime ++ ;
                    n++ ;
                }
                if(nbPrime > maxNbprime) {
                    maxNbprime = nbPrime ;
                    aMax = a;
                    bMax = b ;
                }
            }
        }
        pbR->nbClock = clock() - pbR->nbClock ;
        
        if(pbR->isVerbose)fprintf(stdout,"\t PB%s %d nxn%c%dxn+%d produce %d premiers\n"
                                  ,pbR->ident
                                  ,aMax*bMax
                                  ,(aMax > 0) ? '+' : '-'
                                  ,(aMax > 0) ? aMax : -aMax
                                  , bMax, maxNbprime) ;
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",aMax*bMax) ;
        
    }
    Free_tablePrime(ctxP);
    return 1;
}

int PB028(PB_RESULT *pbR) {
    /*
     43 44 45 46 47 48 49
     42 21 22 23 24 25 26
     41 20  7  8  9 10 27
     40 19  6  1  2 11 28
     39 18  5  4  3 12 29
     38 17 16 15 14 13 30
     37 36 35 34 33 32 31
     
     Quite easy to demonstrate that the answer is a polynom of degree 3 of the half size
     (as the diagonal is odd squares (1 9 25 49 ...) and the other deduce by addition of a * p + b.
     So if size is 2k+1
     P(0) = 1
     P(1) = P(0)+3+5+7+9 = 25
     P(2) = P(1)+13+17+21+25 = 101
     P(3) = P(2)+31+37+43+49 = 261
     
     P(k) = A k**3+ B k**2+  C *p + D
     D = 1 ;
     (1) A+B+C = 24
     (2) 8A+4B+2C= 100 <=> (2') 4A + 2B + C = 50
     (3) 27A+9B+3C = 260
     
     (3)-4*(2)+(1) <=> 6A=260 - 3 * 100 + 3 * 24 = 32
     A = 16/3
     (2')-(1) 3A + B = 26 => B =10
     C =26/3
     
     donc P(k) = (16 k**3 + 30 k**2 + 26 k + 3) / 3
     */
    int k=500 ;
    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",(16*k*k*k + 30*k*k + 26*k + 3) / 3) ;
    return 1 ;
    
}



int PB029(PB_RESULT *pbR) {
    // les seules valeurs qui peuvent avoir des puissances communes sont 2, 3, 5, 6, 7 , 10
    // pour 2 il y a les puissances de 2, 4, 8,16,32,64 <100 soit exp multiple de 1,2,3,4,5,6
    // pour 3 il y a les puissances de 3 9 27 81 soit exp multiple de 1,2,3,4
    // pour 5 => 5 et 25 exp mult de 1,2
    // pour 6 => 6 et 36 exp mult de 1,2
    // pour 7 => 7 et 49 exp mult de 1,2
    // pour 10 => 10 et 100 exp muult de 1,2
    // 4,8,9 sont deja traites des 2 ou 3
    // on va compter les collisions et le total sera 99x99 - nbColl
    int nbCol = 0 ;
    int p ;
    pbR->nbClock = clock();
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    for(p=2;p<=3;p++) {
        int expM ; int n ; // p**expM <= 100
        for(n=p,expM=0; n < 100; expM++) {
            n *= p ;
        }
        {
            int i, exp ;
            uint8_t *isOccupied = calloc(expM*101,sizeof(isOccupied[0])) ;
            for(i=2;i<=100;i++) { isOccupied[i] = 1 ; }
            for(exp=2;exp <= expM; exp++) {
                for(i=2;i<=100;i++) {
                    if(isOccupied[i*exp]) nbCol++ ;
                    else isOccupied[i*exp] = 1 ;
                }
            }
            free(isOccupied) ;
            if(pbR->isVerbose)fprintf(stdout,"(%d,%d,%d)",p,expM,nbCol) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n");
    // pour p = 5 ,6 ,7 , 10 nbCol = 49 ;
    nbCol += 49 * 4 ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",99*99 - nbCol) ;
    return 1 ;
}

int PB030(PB_RESULT *pbR) {
    int32_t pow5[10] ;
    int i, n ;
    int SumRes = 0;
    int Maxn = 0 ;
    pbR->nbClock = clock() ;
    for(i=0;i<10;i++) {
        pow5[i] = i*i*i*i*i ;
    }
    Maxn = 6 * pow5[9] +1 ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s sum[ ",pbR->ident) ;
    for(n=2;n<Maxn;n++) {
        if(n == pow5[n % 10] + pow5[ (n/10) % 10] + pow5[(n/100) % 10] + pow5[(n/1000) % 10] + pow5[(n/10000) % 10] + + pow5[(n/100000)]) {
            SumRes += n ;
            if(pbR->isVerbose) fprintf(stdout,"%d ",n);
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose) fprintf(stdout,"] = %d\n",SumRes) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",SumRes);
    return 1 ;
}

#define PB031_NBV   8
int PB031(PB_RESULT *pbR) {
    int value[PB031_NBV] = { 200 , 100 , 50, 20 , 10 ,5 , 2 ,1};
    int nbPieces[PB031_NBV] ;
    int nbDecomp = 0 ;
    int np , sum ;
    pbR->nbClock = clock() ;
    np = -1 ;
    sum = 200 ;
    // on parcourt la decomposition dans l'ordre lexicographique
    while(1) {
        if(sum == 0) {
            nbDecomp++ ;
            //          { int i ;    for(i=0;i<PB031_NBV;i++) { printf("%c%dx%d",(i==0) ? '\n' : '+' ,nbPieces[i],value[i]) ; }   }
            while(sum < value[np]) {
                sum += nbPieces[np] * value[np] ;
                nbPieces[np] = 0 ;
                if(--np < 0 ) break ;
            }
            if(np < 0) {
                break ;
            }
            nbPieces[np]++ ;
            sum -= value[np] ;
        }
        while(++np < PB031_NBV-1) {
            nbPieces[np] = 0 ;
        }
        // derniere piece de 1
        nbPieces[np] = sum ;
        sum = 0 ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbDecomp) ;
    return 1 ;
}


int PB032(PB_RESULT *pbR) {
    uint8_t dg[9] ;
    int Sum = 0 ;
    pbR->nbClock = clock() ;
    { int i ; for(i=0;i<9;i++) {  dg[i] = i+1 ; } }
    int rg = 0 ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    do {
        int nb4 = 1000*dg[0]+100*dg[1]+10*dg[2]+dg[3]  ;
        do {
            // on verifie si decomposition OK ; c'est forcement les 4 premiers = dddd = d x dddd ou dddd = dd x ddd
            if(nb4 == dg[4] * (1000*dg[5]+100*dg[6]+10*dg[7]+dg[8]) ) {
                Sum += nb4 ;
                if(pbR->isVerbose)fprintf(stdout," %d=%dx%d%d%d%d ",nb4,dg[4],dg[5],dg[6],dg[7],dg[8]) ;
                while((rg = NextPermut(dg,9) ) >= 4) ;
                
            } else if(nb4 == (10*dg[4]+dg[5]) * (100*dg[6]+10*dg[7]+dg[8]) ) {
                Sum += nb4 ;
                if(pbR->isVerbose)fprintf(stdout," %d=%d%dx%d%d%d ",nb4,dg[4],dg[5],dg[6],dg[7],dg[8]) ;
                while((rg = NextPermut(dg,9) ) >= 4) ;
            } else {
                rg=NextPermut(dg,9) ;
            }
        } while ( rg >= 4) ;
    } while(rg >= 0 ) ;
    if(pbR->isVerbose)fprintf(stdout,"\n");
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Sum) ;
    return 1 ;
}

int PB033(PB_RESULT *pbR) {
    int a,n,d ;
    int Num = 1;
    int Den = 1;
    pbR->nbClock = clock();
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    for(n=1;n<9;n++) {
        for(d=n+1;d<10;d++) {
            // 2 cas (an/da = n/d ou na/ad = n/d)
            for(a=1;a<10;a++) {
                if ( (10*a+n)*d == (10*d+a)*n ) {
                    Num *= n ; Den *= d ;
                    if(pbR->isVerbose)printf(" %d%d/%d%d=%d/%d",a,n,d,a,n,d) ;
                } else if ( (10*n+a)*d == (10*a+d)*n ) {
                    Num *= n ; Den *= d ;
                    if(pbR->isVerbose)printf(" %d%d/%d%d=%d/%d",n,a,a,d,n,d) ;
                }
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    {
        int pgcd = PGCD(Num,Den);
        Num /= pgcd ; Den /=pgcd ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%s Prod=%d/%d\n",pbR->ident,Num,Den) ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Den);
    return 1 ;
}

// max value = 9! x 7 ;
#define PB034_MAX 2540160
#define PB034_MAXDIG    7
int PB034(PB_RESULT *pbR) {
    int fact[10] , diffFact[10] ;
    int dig[PB034_MAXDIG] ;
    int i, n , k , nbDig ;
    int  Diff,Sum = 0 ;
    pbR->nbClock = clock()  ;
    fact[0] = 1 ;
    for(i=1;i<=9;i++) {
        fact[i] =fact[i-1] * i ;
    }
    for(i=1;i<=9;i++) {
        diffFact[i] = fact[i] - fact[i-1] ;
    }
    diffFact[0] = fact[0] - fact[9] ;
    for(i=0;i<PB034_MAXDIG;i++) { dig[i] = 0 ; }
    n = 3 ;
    dig[0] = n ;
    Diff = fact[dig[0]] - n ;
    nbDig = 1 ;
    k = 0 ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    while (n < PB034_MAX){
        if(Diff == 0) {
            Sum += n ;
            if(pbR->isVerbose) {
                fprintf(stdout," %d",n); for(i=nbDig-1;i>=0;i--) { fprintf(stdout,"%c%d!",(i==nbDig-1) ? '=' : '+',dig[i]) ; }
            }
        }
        while(dig[k] == 9) {
            dig[k++] = 0;
            Diff += diffFact[0] ;
        }
        if(k >= nbDig) {
            dig[k] = 1 ;
            nbDig++ ;
        } else {
            dig[k]++ ; Diff += diffFact[dig[k]]  - 1;
        }
        k = 0 ;
        n++ ;
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Sum) ;
    return 1 ;
    
}

#define PB035_MAXN  1000000
int PB035(PB_RESULT *pbR) {
    int i , puis10 ;
    int nbFind = 0 ;
    uint8_t    *isPrimeUsed = NULL ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB035_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    uint32_t nbPrime = GetNbPrime(ctxP) ;
    const T_prime * tbPrime = GetTbPrime(ctxP);
    isPrimeUsed = calloc(nbPrime,sizeof(isPrimeUsed[0])) ;
    
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s",pbR->ident) ;
    for(i=0,puis10=1;i<nbPrime;i++) {
        T_prime n = tbPrime[i] ;
        if(isPrimeUsed[i] == 0) { // on met a zero ceux deja trouve
            int lgCyclePrime = 1 ;
            isPrimeUsed[i] = 1 ;
            T_prime nr = n ;
            while(10*puis10 <= n ) {
                puis10 *= 10 ;
            }
            // si dans le cycle on aboutit a une valeur inferieure c'est deja crame
            while ((nr = (nr % 10) * puis10  + nr /10) > n) {
                int rg  ;
                // il faur nr premier
                if((rg =SearchRg_TablePrime(ctxP,nr)) >= 0) {
                    lgCyclePrime++ ;
                    isPrimeUsed[rg] = 1 ;
                } else {
                    lgCyclePrime = 0 ;
                    break ;
                }
            }
            
            if(nr == n) {
                if(pbR->isVerbose)fprintf(stdout," %u",n) ;
                nbFind += lgCyclePrime ;
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbFind);
    free(isPrimeUsed) ;
    return 1 ;
}

static int IsPal2(uint32_t n) {
    uint32_t B=1, H = 1 << 31 ;
    if(n<=1) return 1 ;
    while(H>n) { H >>= 1 ; }
    while(H>B) {
        uint32_t cmp = (H+B) & n ;
        if(cmp == H || cmp == B) return 0  ;
        H >>= 1;
        B <<= 1 ;
    }
    return 1 ;
}

int PB036(PB_RESULT *pbR) {
    int i1,i2,i3 ;
    int Sum = 0 ;
    int nbFind = 0 ;
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s ",pbR->ident) ;
    for(i1=1;i1<10;i1++) {
        int n=i1;
        if(IsPal2(n)) {
            if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
            nbFind++ ;
            Sum += n;
        }
        n=i1*11 ;
        if(IsPal2(n)) {
            if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
            nbFind++ ;
            Sum += n;
        }
    }
    for(i1=1;i1<10;i1++) {
        for(i2=0;i2<10;i2++) {
            int n=i1*101 + i2*10 ;
            if(IsPal2(n)) {
                nbFind++ ;
                if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                Sum += n;
            }
            n= i1*1001 + i2 * 110 ;
            if(IsPal2(n)) {
                nbFind++ ;
                if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                Sum += n;
            }
        }
    }
    for(i1=1;i1<10;i1++) {
        for(i2=0;i2<10;i2++) {
            for(i3=0;i3<10;i3++) {
                int n=i1*10001 + i2*1010 + i3 * 100 ;
                if(IsPal2(n)) {
                    if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                    nbFind++ ;
                    Sum += n;
                }
                n= i1*100001 + i2 * 10010 + i3 * 1100 ;
                if(IsPal2(n)) {
                    if(pbR->isVerbose)fprintf(stdout,"%d ",n) ;
                    nbFind++ ;
                    Sum += n;
                }
            }
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%s Nbfind=%d\n",pbR->ident,nbFind) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Sum);
    return 1 ;
}

#define PB037_MAXP  1000000
#define PB037_MAX_NB    11
static int32_t CmpPrime037(const void *p1,const void *p2) {
    return *((int32_t *)p1) - *((int32_t *)p2) ;
}

int PB037(PB_RESULT *pbR) {
    int  truncPrime[PB037_MAX_NB] ;
    int nbAlloc = 1024 ;
    int *truncRight = malloc(nbAlloc*sizeof(truncRight[0])) ;
    int *truncLeft = malloc(nbAlloc*sizeof(truncLeft[0])) ;
    int nbRight = 0 , nbLeft = 0;
    int i,k,indm,indm1 = 0 , indrm, indrm1 = 0 , indlm , indlm1 = 0;
    int nbFind = 0 ;
    int puis10 = 10 ;
    int Sum = 0 ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB037_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    // on tabule ceux a 2 digits (composition de 2,3,5,7)
    truncPrime[nbFind] = 23 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 37 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 53 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 73 ; Sum += truncPrime[nbFind++] ;
    // on tbule les right a 2 digits
    truncRight[nbRight++] = 23 ;
    truncRight[nbRight++] = 29 ;
    truncRight[nbRight++] = 31 ;
    truncRight[nbRight++] = 37 ;
    truncRight[nbRight++] = 53 ;
    truncRight[nbRight++] = 59 ;
    truncRight[nbRight++] = 71 ;
    truncRight[nbRight++] = 73 ;
    truncRight[nbRight++] = 79 ;
    // on tabule les left a 2 digits
    truncLeft[nbLeft++] = 13 ;
    truncLeft[nbLeft++] = 17 ;
    truncLeft[nbLeft++] = 23 ;
    truncLeft[nbLeft++] = 37 ;
    truncLeft[nbLeft++] = 43 ;
    truncLeft[nbLeft++] = 47 ;
    truncLeft[nbLeft++] = 53 ;
    truncLeft[nbLeft++] = 73 ;
    truncLeft[nbLeft++] = 83 ;
    truncLeft[nbLeft++] = 97 ;
    do {
        indm = nbFind ;
        indrm = nbRight ;
        indlm = nbLeft ;
        puis10 *= 10 ;
        if(puis10*10 > PB037_MAXP) {
            free(truncRight); free(truncLeft) ;
            if(pbR->isVerbose) {
                fprintf(stdout,"\n\t PB%s ",pbR->ident) ;
                for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
            }
            return 0 ;
        }
        if(4*(indrm-indrm1)+ nbRight >= nbAlloc || 4*(indlm-indlm1)+ nbLeft >= nbAlloc) {
            nbAlloc *= 2 ;
            truncRight = realloc(truncRight,nbAlloc*sizeof(truncRight[0])) ;
            truncLeft = realloc(truncLeft,nbAlloc*sizeof(truncLeft[0])) ;
        }
        for(i=indrm1;i<indrm;i++) {
            // on va le completer avec des 1, 3, 7, 9
            uint32_t p = 10*truncRight[i]+1 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p += 3-1 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p+= 7-3 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p += 9 - 7 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
        }
        for(k=1;k<=9;k++) {
            for(i=indlm1;i<indlm;i++) {
                // on va le completer avec de 1 a 9
                uint32_t p = truncLeft[i] + k * puis10  ;
                if(Search_TablePrime(ctxP,p)) { truncLeft[nbLeft++] = p ; }
            }
        }
        for(i=indrm;i<nbRight;i++) {
            uint32_t p = truncRight[i] ;
            if(bsearch(&p,truncLeft+indlm,nbLeft-indlm,sizeof(p),CmpPrime037)) {
                truncPrime[nbFind++] = truncRight[i] ;
                Sum += truncRight[i]  ;
            }
        }
        indm1 = indm ;
        indrm1 = indrm ;
        indlm1 = indlm ;
    } while(nbFind < PB037_MAX_NB ) ;
    if(pbR->isVerbose) {
        fprintf(stdout,"\t PB%s right=%d left=%d RL=%d\n",pbR->ident,nbRight,nbLeft,nbFind) ;
        fprintf(stdout,"\t PB%s ",pbR->ident) ;
        for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Sum);
    free(truncRight);
    return 1 ;
}

int PB037r(PB_RESULT *pbR) {
    int  truncPrime[PB037_MAX_NB] ;
    int nbAlloc = 1024 ;
    int *truncRight = malloc(nbAlloc*sizeof(truncRight[0])) ;
    int nbRight = 0 ;
    int i,indm,indm1 = 0 , indrm, indrm1 = 0 ;
    int nbFind = 0 ;
    int puis10 = 10 ;
    int Sum = 0 ;
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB037_MAXP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    // on tabule ceux a 2 digits (composition de 2,3,5,7)
    truncPrime[nbFind] = 23 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 37 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 53 ; Sum += truncPrime[nbFind++] ;
    truncPrime[nbFind] = 73 ; Sum += truncPrime[nbFind++] ;
    // on tbule les right a 2 digits
    truncRight[nbRight++] = 23 ;
    truncRight[nbRight++] = 29 ;
    truncRight[nbRight++] = 31 ;
    truncRight[nbRight++] = 37 ;
    truncRight[nbRight++] = 53 ;
    truncRight[nbRight++] = 59 ;
    truncRight[nbRight++] = 71 ;
    truncRight[nbRight++] = 73 ;
    truncRight[nbRight++] = 79 ;
    do {
        indm = nbFind ;
        indrm = nbRight ;
        puis10 *= 10 ;
        if(puis10*10 > PB037_MAXP) {
            free(truncRight);
            if(pbR->isVerbose) {
                fprintf(stdout,"\n\t PB%s ",pbR->ident) ;
                for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
            }
            return 0 ;
        }
        if(4*(indrm-indrm1)+ nbRight >= nbAlloc ) {
            nbAlloc *= 2 ;
            truncRight = realloc(truncRight,nbAlloc*sizeof(truncRight[0])) ;
        }
        for(i=indrm1;i<indrm;i++) {
            // on va le completer avec des 1, 3, 7, 9
            uint32_t p = 10*truncRight[i]+1 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p += 3-1 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p+= 7-3 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
            p += 9 - 7 ;
            if(Search_TablePrime(ctxP,p)) { truncRight[nbRight++] = p ; }
        }
        for(i=indrm;i<nbRight;i++) {
            uint32_t k = puis10 ;
            uint32_t p = truncRight[i] ;
            while( k >= 10) {
                p = p % k ;
                if(!Search_TablePrime(ctxP,p))  break ;
                k /= 10 ;
            }
            if ( k == 1) {
                truncPrime[nbFind++] = truncRight[i] ;
                Sum += truncRight[i]  ;
            }
        }
        indm1 = indm ;
        indrm1 = indrm ;
    } while(nbFind < PB037_MAX_NB ) ;
    if(pbR->isVerbose) {
        fprintf(stdout,"\t PB%s right=%d RL=%d\n",pbR->ident,nbRight,nbFind) ;
        fprintf(stdout,"\t PB%s ",pbR->ident) ;
        for(i=0;i<nbFind;i++) fprintf(stdout,"%d%c",truncPrime[i], (i == nbFind-1) ? '\n' : ' ') ;
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Sum);
    free(truncRight);
    return 1 ;
}


// max value = 9! x 7 ;

int PB038(PB_RESULT *pbR) {
    uint8_t dig[10] ;
    int i , k ;
    int maxP = 0 ;
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose)fprintf(stdout,"\t PB%s (1,2)",pbR->ident) ;
    // XXXX par 1,2 XXXX =>4  XXXXx2 =>5
    for(i=5000;i<=10000;i++) {
        int n1 = i * 100002 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout," %d->%d",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%s (1,2,3)",pbR->ident) ;
    // XXX par 1 , 2 , 3  XXX => 3c XXXx2=>3c XXXx3=>3c
    for(i=100;i<=333;i++) {
        int n1 = i * 1002003 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%s (1,2,3,4)",pbR->ident) ;
    // XX par 1,2,3,4  XX=>2c XXx2=>2c XXx3=>2c XXx4=>3c
    for(i=25;i<=33;i++) {
        int n1 = i * 10203004 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }
    if(pbR->isVerbose)fprintf(stdout,"\n\t PB%s (1,2,3,4,5)",pbR->ident) ;
    // X par 1,2,3,4,5  X=>1c  Xx2=>2c Xx3=>2c Xx4=>2c Xx4=>2c
    for(i=1;i<=9;i++) {
        int n1 = i * 102030405 ;
        int n = n1 ;
        memset(dig,0,10) ;
        for(k=0;k<9;k++) {
            int r = n % 10 ; n /= 10 ;
            if(r==0 || dig[r] != 0) break ;
            dig[r] = 1 ;
        }
        if(k==9 ) {
            if(n1 > maxP) maxP = n1 ;
            if(pbR->isVerbose) fprintf(stdout,"%d->%d ",i,n1) ;
        }
    }
    
    
    if(pbR->isVerbose)fprintf(stdout,"\n") ;
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",maxP) ;
    return 1 ;
    
}

#define PB039_MAX   1000

int PB039(PB_RESULT *pbR) {
    // pour m > n > 0
    // a = (m**2 - n**2)*k ; b = 2*m*n*k ; c = (m**2 + n**2)*k ;
    // donc S = 2 * (m**2 + m*n)*k = 2 * m * (m+n) * k
    int32_t m , n , P , nMax =0 , Pbest=1 ;
    int32_t *nbSol = calloc(PB039_MAX/2+1,sizeof(nbSol[0])) ;
    pbR->nbClock = clock()  ;
    
    for(n=1;n*(2*n+1)<=PB039_MAX/2;n++) {
        for(m=n+1;(P=m*(m+n))<=PB039_MAX/2;m++) {
            if( PGCD(m,n)== 1 && (((m & n) & 1) == 0 )) {
                do {
                    nbSol[P]++ ;
                    P += m*(m+n) ;
                } while(P <= PB039_MAX/2 ) ;
            }
        }
    }
    for(n=1;n<PB039_MAX/2;n++) {
        if(nbSol[n]> nMax) {
            nMax = nbSol[n] ;
            Pbest = n ;
        }
    }
    if(pbR->isVerbose) {
        int mmn ;
        fprintf(stdout,"\t PB%s P=%d",pbR->ident,2*Pbest);
        for(mmn=2; mmn<Pbest;mmn++) {
            if(nbSol[mmn] == 1 && (Pbest % mmn) == 0) {
                int k = Pbest / mmn ;
                for(m=2;m*m< mmn; m++) {
                    if((mmn % m ) == 0) {
                        n = mmn / m  - m ;
                        if(m>n && PGCD(m,n)== 1 && (((m & n) & 1) == 0 )) {
                            fprintf(stdout," (%d,%d,%d)",k*(m*m - n*n) ,k*2*m*n,k*(m*m + n*n));
                        }
                    }
                }
            }
        }
        fprintf(stdout,"\n");
    }
    pbR->nbClock = clock() -  pbR->nbClock ;
    
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",2*Pbest) ;
    free(nbSol);
    return 1;
}


#define PB040_MAXN  7
#define PB040_IND_MAX   1000000
int PB040(PB_RESULT *pbR) {
    int i,nb,P=1,ind;
    uint32_t maxByDig[PB040_MAXN+1] ;
    uint32_t offsByDig[PB040_MAXN+1] ;
    pbR->nbClock = clock()  ;
    maxByDig[0] = 1 ; nb = 9 ; offsByDig[0] = -1 ;
    for(i=1;i<=PB040_MAXN;i++) {
        maxByDig[i] = maxByDig[i-1] + nb * i ;
        nb *= 10 ;
        offsByDig[i] = (offsByDig[i-1]+1)*10 ;
    }
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s P",pbR->ident);
    
    for(ind=1;ind<=PB040_IND_MAX;ind *=10) {
        int n, r ;
        for(i=1;maxByDig[i]<= ind;i++) ;
        n = (ind+offsByDig[i]) / i ;
        r = (ind+offsByDig[i]) % i ;
        while(++r < i) {
            n /= 10 ;
        }
        P *= n % 10 ;
        if(pbR->isVerbose) fprintf(stdout,"%c%d",(ind==1) ? '=' : 'x' , n % 10) ;
        
    }
    if(pbR->isVerbose) fprintf(stdout,"\n");
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",P);
    return 1 ;
}


// on commence a 7
// car permut de 8 et 9 sont multiples de 3
#define PB041_MAXN  7654321
int PB041(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    uint8_t permut[7] = {7,6,5,4,3,2,1} ;
    T_prime pBest = 0 ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB041_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    do {
        if((permut[6] & 1) && (permut[6] != 5)) { // pas pair ni multiple de 5
            T_prime p = permut[0]*1000000+permut[1]*100000+permut[2]*10000+permut[3]*1000+permut[4]*100+permut[5]*10+permut[6] ;
            if( Search_TablePrime(ctxP,p)) {
                pBest = p ;
                break ;
            }
        }
    } while(NextPermutRev(permut,7) >= 0 ) ;
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pBest) {
        snprintf(pbR->strRes, sizeof(pbR->strRes),"%u",pBest);
        return 1 ;
    } else {
        return 0 ;
        
    }
}

#include "p042_words.h"
int PB042(PB_RESULT *pbR) {
    const char ** P042_W = P042_GetData() ;
    const char *word ;
    int nbTriangle = 0 ;
    pbR->nbClock = clock()  ;
    while((word = *P042_W++) != NULL ) {
        int32_t val=0 ;
        int32_t n ;
        uint8_t c ;
        while((c=*word++)) {
            val += c - 'A' + 1 ;
        }
        val *= 2 ;
        n = (int32_t)Sqrt64(val) ;
        if(val == n*(n+1)) {
            nbTriangle++ ;
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",nbTriangle);
    return 1 ;
}






int PB043(PB_RESULT *pbR) {
    // attention 0 est le digit de poids faible, 9 celui de poids fort
    uint8_t dig[10] = {0,1,2,3,4,5,6,7,8,9} ;
    int rg  ;
    int nbL = 0 ;
    uint64_t Sum = 0 ;
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s ",pbR->ident);
    do  {
        if( ( (dig[0]+10*dig[1]+100*dig[2]) % 17 ) != 0 ) {
            // passer au suivant pour les 3 digits de poids faible
            rg = 2 ;
        } else if ( ( (dig[1]+10*dig[2]+100*dig[3]) % 13 ) != 0 ) {
            rg = 3 ;
        } else if ( ( (dig[2]+10*dig[3]+100*dig[4]) % 11 ) != 0 ) {
            rg = 4 ;
        } else if ( ( dig[4] % 5) != 0 ) {
            // test juste sur le dig4 car 10 et 100 divisible par 5
            rg = 4 ;
        } else if ( ( (dig[3]+10*dig[4]+100*dig[5]) % 7 ) != 0 ) {
            rg = 5 ;
        } else if ( ( dig[6] % 2 ) != 0 ) {
            // doit etre pair
            rg = 6 ;
        } else if ( ( (dig[5]+dig[6]+dig[7]) % 3 ) != 0 ) {
            rg = 7 ;
        } else {
            uint64_t N = dig[0]+10*(dig[1]+10*(dig[2]+10*(dig[3]+10*(dig[4]+10*(dig[5]+10*(dig[6]+10*(dig[7]+10*(dig[8]+10LL*dig[9])))))))) ;
            if(pbR->isVerbose) fprintf(stdout,"%lld ",N);
            rg = 8 ;
            Sum += N ;
        }
        nbL++ ;
        
    } while(NextPermutRg(dig, 10, rg) >= 0) ;
    if(pbR->isVerbose) fprintf(stdout,"nbLoop=%d\n",nbL);
    pbR->nbClock = clock() - pbR->nbClock ;
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    return 1 ;
}

uint32_t IsPentagonal ( uint64_t P) {
    uint32_t n = (uint32_t) Sqrt64((2*P)/3) + 1 ;
    if(P*2 == ((uint64_t)n )*(3*n -1)) {
        return n ;
    } else {
        return 0 ;
    }
}

typedef  int64_t PENT_T  ;
int PB044(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s ",pbR->ident);
    {
        PENT_T Pi,Pj ;
        int32_t  i,j,k,r ;
        for(k=2;;k++) {
            int32_t d ;
            int32_t Pk2 = k*(3*k-1) ;
            // on utilise l'idee que si Pk = Pi - Pj on a
            // k(3k-1) = (i-j) * (3i + 3j - 1)
            // donc on peut chercher les factorisation de k(3k-1) avec :
            // i-j = d ; 3i+3j-1 = Pk2/d <=> j = (Pk2/d-3*d+1)/6 et j > 0
            for(d=1; 0 < Pk2 - 3*d*d + d ;d++)  { // boucle sur les diviseurs
                if((Pk2 % d)==0 ) {
                    j = (Pk2/d-3*d+1)/6 ;
                    i = j + d ;
                    Pi = (i*(PENT_T)(3*i-1))/2 ;
                    Pj = (j*(PENT_T)(3*j-1))/2 ;
                    if((r=IsPentagonal(Pi+Pj)) && IsPentagonal(Pi-Pj) ) {
                        pbR->nbClock = clock() - pbR->nbClock ;
                        if(pbR->isVerbose)fprintf(stdout,"P(%d)=%d P(%d)=%lld P(%d)=%lld=%lld+%lld P(%d)=%lld=%lld+%lld\n",k,Pk2/2,j,Pj,i,Pi,Pi-Pj,Pj,r,Pi+Pj,Pi,Pj) ;
                        snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",Pk2/2);
                        return 1 ;
                        
                    }
                }
            }
        }
    }
}

// Tn = n(n+1)/2
// Hn = n(2n-1) = 2n(2n-1)/2
// Pn = n(3n-1)/2
// les hexagonaux sont les triangles de rang impair
// donc on a va chercher un pentagonal qui est aussi un trinagle impair
#define PB045_NBS   2
int PB045(PB_RESULT *pbR) {
    pbR->nbClock = clock()  ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s ",pbR->ident);
    {
        PENT_T Pk2 ;
        int32_t  k, nbSol=PB045_NBS;
        for(k=2;nbSol>0;k++) {
            Pk2 = k*(PENT_T)(3*k-1) ;
            int32_t d = (int32_t) Sqrt64(Pk2) ;
            if((d&1) && (d*(PENT_T)(d+1) == Pk2) ) {
                if(pbR->isVerbose)fprintf(stdout,"P(%d)=%lld=T(%d)=H(%d) ",k,Pk2/2,d,(d+1)/2) ;
                nbSol-- ;
                if(nbSol == PB045_NBS - 2) snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Pk2/2);
            }
        }
    }
    pbR->nbClock = clock() - pbR->nbClock ;
    if(pbR->isVerbose)fprintf(stdout,"\n");
    return 1 ;
}

#define PB046_MAXN  100000
int PB046(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    uint8_t * oddDec = calloc(PB046_MAXN/2, sizeof(oddDec[0])) ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB046_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    {
        int i ;
        int curOdd = 1;
        uint32_t nbPrime = GetNbPrime(ctxP);
        const T_prime * tbPrime= GetTbPrime(ctxP);
        for(i=1;i<nbPrime;i++) {
            int n ;
            int32_t p = (int32_t) tbPrime[i] ;
            int indP = p/2 ;
            while(curOdd < indP) {
                if(!oddDec[curOdd]) {
                    Free_tablePrime(ctxP) ;
                    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",2*curOdd+1) ;
                    pbR->nbClock = clock() - pbR->nbClock ;
                    return 1 ;
                }
                curOdd++ ;
            }
            curOdd++ ;
            indP += 1 ; // (2* 1*1 /2)
            for(n=1;indP<PB046_MAXN/2;n++) {
                oddDec[indP] = 1 ;
                indP += 2*n+1 ;
            }
        }
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}

// donc test jusqu'a 10000*10000
#define PB047_MAXN  10000
int PB047(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB047_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    {
        int n=2 ;
        const T_prime * tbPrime = GetTbPrime(ctxP);
        for(n=2;n<PB047_MAXN*PB047_MAXN;) {
            if(FindNbDivPrime(n,tbPrime) < 4 ) {
                n++ ; continue ;
            }
            if(FindNbDivPrime(n+1,tbPrime) < 4 ) {
                n += 2 ; continue ;
            }
            if(FindNbDivPrime(n+2,tbPrime) < 4 ) {
                n += 3 ; continue ;
            }
            if(FindNbDivPrime(n+3,tbPrime) < 4 ) {
                n += 4 ; continue ;
            }
            snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",n);
            Free_tablePrime(ctxP) ;
            pbR->nbClock = clock() - pbR->nbClock ;
            return 1 ;
        }
    }
    Free_tablePrime(ctxP) ;
    pbR->nbClock = clock() - pbR->nbClock ;
    return 0 ;
}


#define PB048_MAXN  1000
#define PB048_MOD   10000000000LL // dix derniet digits
#define PB048_MODHALF   100000
uint64_t mult10d(uint64_t a, uint64_t b) {
    uint64_t al = a % PB048_MODHALF ;
    uint64_t ah = ( a / PB048_MODHALF ) % PB048_MODHALF ;
    uint64_t bl = b % PB048_MODHALF ;
    uint64_t bh = ( b / PB048_MODHALF ) % PB048_MODHALF ;
    return (al * bl + (al * bh + ah * bl) * PB048_MODHALF) % PB048_MOD ;
}
int PB048(PB_RESULT *pbR) {
    uint64_t Sum = 0 ;
    int n;
    pbR->nbClock = clock()  ;
    for(n=1;n<=PB048_MAXN;n++) {
        int exp ;
        uint64_t pow = n ;
        uint64_t Npow= 1 ;
        for(exp=n ;exp != 0; exp >>=1) {
            if(exp & 1) Npow = mult10d(Npow,pow)  ;
            pow = mult10d(pow,pow)  ;
        }
        Sum = (Sum + Npow) % PB048_MOD ;
    }
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}

typedef struct perm4d {
    uint16_t nb ;
    uint16_t  perm[24] ;
} perm4d;

#define PB049_MAXN  10000
#define PB049_SOL1  1487
int PB049(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    perm4d *Permut ;
    int i;
    pbR->nbClock = clock()  ;
    if((ctxP = Gen_tablePrime(PB047_MAXN)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    Permut = calloc(PB049_MAXN,sizeof(Permut[0]));
    uint32_t nbPrime = GetNbPrime(ctxP);
    const T_prime * tbPrime= GetTbPrime(ctxP);
    
    for(i=0;tbPrime[i]<1000;i++) ;
    for(;i<nbPrime;i++) {
        int num ;
        uint8_t dig[4];
        int p = (int) tbPrime[i] ;
        dig[0] = p / 1000 ;
        dig[1] = (p/100) % 10 ;
        dig[2] = (p/10) % 10 ;
        dig[3] = p % 10 ;
        HeapSortUint8(dig,4);
        num = dig[0] * 1000 + dig[1] * 100 + dig[2]*10 + dig[3] ;
        Permut[num].perm[Permut[num].nb++] = p;
    }
    for(i=1000;i<PB049_MAXN;i++) {
        int nb = Permut[i].nb ;
        int i1,i2,i3 ;
        if(nb < 3 || Permut[i].perm[0] == PB049_SOL1) continue ;
        for(i1=0;i1<nb-2;i1++) {
            for(i2=i1+1;i2<nb-1; i2++) {
                int e3 = 2*Permut[i].perm[i2] - Permut[i].perm[i1] ;
                for(i3=i2+1;i3<nb; i3++) {
                    if( Permut[i].perm[i3] == e3) {
                        snprintf(pbR->strRes, sizeof(pbR->strRes),"%d%d%d",Permut[i].perm[i1],Permut[i].perm[i2],Permut[i].perm[i3]) ;
                        i = PB049_MAXN ; i3= i2 = i1 = nb ;
                    }
                }
            }
        }
        
    }
    //    snprintf(pbR->strRes, sizeof(pbR->strRes),"%lld",Sum);
    free(Permut);
    Free_tablePrime(ctxP);
    pbR->nbClock = clock() - pbR->nbClock ;
    return 1 ;
}


#define PB050_MAXN  1000000
int PB050(PB_RESULT *pbR) {
    CTX_PRIMETABLE * ctxP  ;
    int32_t indP0,indP1 ;
    int32_t S ;
    uint32_t bestIndP0 = 0 ;
    uint32_t bestLg = 0 ;
    uint32_t bestP = 0 ;
    pbR->nbClock = clock()  ;
    // we precompute N = sqrt(Max) * log2(Max)  primes
    uint32_t sizeP = (uint32_t)Sqrt64(PB050_MAXN) ;
    {   int i=1, nb = 1;
        while((i *= 2) < PB050_MAXN ) { nb++ ; }
        sizeP *= nb ;
    }
    if((ctxP = Gen_tablePrime(sizeP)) == NULL) {
        fprintf(stdout,"\t PB%s Fail to alloc prime table\n",pbR->ident);
        return 0 ;
    }
    const T_prime * tbPrime= GetTbPrime(ctxP);
    
    indP0 = 0; //first prime of the sum is index=0 (2)
    indP1 = indP0 ; S=0;
    do { // loop on the first prime
        int32_t indP2 ;
        int32_t S1 ;
        // on cherche lindex indP1 S > AMX
        while(S<PB050_MAXN) {
            S += tbPrime[indP1++] ;
        }
        // decrease the chain by the end until we reach a prime
        indP2 = indP1;
        S1 = 0 ;
        while(indP2>=indP0+bestLg) { //check not poor result
            S1 += tbPrime[--indP2] ;
            if(Is_Prime(S-S1,tbPrime)){
                bestLg =indP2 - indP0 ;
                bestIndP0 = indP0 ;
                bestP = S-S1 ;
                break ;
            }
        }
        // increment the index of the first prime in the sum
        S -= tbPrime[indP0++] ;
        //exit if the chain beginning at indP0 and length=bestLg
        // has a SUM exeeding MAX. It will be worse for superior indP0
    } while (indP0+bestLg < indP1) ;
    if(pbR->isVerbose) fprintf(stdout,"\t PB%s lg=%d %d=sum(%u..%u)\n",pbR->ident,bestLg,bestP,
                               tbPrime[bestIndP0],tbPrime[bestIndP0+bestLg-1]);
    
    pbR->nbClock = clock() - pbR->nbClock ;
    Free_tablePrime(ctxP);
    snprintf(pbR->strRes, sizeof(pbR->strRes),"%d",bestP);
    return 1 ;
}
