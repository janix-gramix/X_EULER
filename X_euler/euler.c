//
//  euler.c
//  EulerProject
//  Created by Jeannot on 26/01/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "euler_utils.h"
#include "PB001_050.h"
#include "PB051_100.h"
#include "PB101_150.h"
#include "PB_gmp.h"
#include "PB_other.h"





typedef struct TotalRun {
    clock_t TotalClock ;
    int isVerbose ;
    int nbPBOK ;
    int nbPBerror ;
} TotalRun ;


void Execute(TotalRun *ttr, PB_CALL *pbCall) {
    PB_RESULT pbR ;
    memset(&pbR,0,sizeof(pbR)) ;
    pbR.isVerbose = ttr->isVerbose ;
    if( pbCall->pbSolve != NULL) {
        pbR.pbNum = pbCall->pbNum ;
        pbR.nbClock = 0 ;
        if(pbCall->pbSolve(&pbR)) {
            ttr->TotalClock += pbR.nbClock ;
            if(strcmp(pbCall->Solution,pbR.strRes)==0) {
                ttr->nbPBOK++ ;
                fprintf(stdout,"OK\tPB%03d(%.06fs) Sol=%s '%s'\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,(pbCall->name) ? pbCall->name : "") ;
            } else {
                ttr->nbPBerror++ ;
                fprintf(stdout,"ERROR\tPB%03d(%.06fs) Find=%s != Exp=%s\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,pbCall->Solution) ;
            }
        } else {
            ttr->TotalClock += pbR.nbClock ;
            ttr->nbPBerror++ ;
            fprintf(stdout," BAD EXECUTION PB%03d(%.06fs)\n",pbCall->pbNum,(float)pbR.nbClock / CLOCKS_PER_SEC) ;
        }
        
    } else {
        fprintf(stdout,"SKIPPED\tPB%3d\n",pbCall->pbNum) ;
    }

    
}

static PB_CALL ALL_calls[] = {
    {  1,  PB001,  "233168"}
    ,{  2,  PB002,  "4613732"}
    ,{  3,  PB003,  "6857"}
    ,{  4,  PB004,  "906609"}
    ,{  5,  PB005,  "232792560"}
    ,{  6,  PB006,  "25164150"}
    ,{  7,  PB007,  "104743"}
    ,{  8,  PB008,  "23514624000"}
    ,{  9,  PB009,  "31875000"}
    ,{ 10,  PB010,  "142913828922","Summation of primes"}
    ,{ 11,  PB011,  "70600674"}
    ,{ 12,  PB012,  "76576500"}
    ,{ 13,  PB013,  "5537376230"}
    ,{ 14,  PB014,  "837799"}
    ,{ 15,  PB015,  "137846528820"}
    ,{ 16,  PB016,  "1366"}
    ,{ 16,  PB016_gmp,  "1366"}
    //        ,{ 17,  PB017,  ""}
    ,{ 18,  PB018,  "1074"}
    //        ,{ 19,  PB019,  ""}
    ,{ 20,  PB020,  "648"}
    ,{ 21,  PB021,  "31626"}
    //        ,{ 22,  PB022,  ""}
    ,{ 23,  PB023,  "4179871"}
    ,{ 24,  PB024,  "2783915460"}
    ,{ 25,  PB025,  "4782"}
    ,{ 26,  PB026,  "983"}
    ,{ 27,  PB027,  "-59231"}
    ,{ 28,  PB028,  "669171001"}
    ,{ 29,  PB029,  "9183"}
    ,{ 30,  PB030,  "443839"}
    ,{ 31,  PB031,  "73682"}
    ,{ 32,  PB032,  "45228"}
    ,{ 33,  PB033,  "100"}
    ,{ 34,  PB034,  "40730"}
    ,{ 35,  PB035,  "55"}
    ,{ 36,  PB036,  "872187"}
    ,{ 37,  PB037,  "748317"}
    ,{ 37,  PB037r, "748317"}
    ,{ 38,  PB038,  "932718654"}
    ,{ 39,  PB039,  "840","Integer right triangles"}
    ,{ 40,  PB040,  "210"}
    ,{ 41,  PB041,  "7652413"}
    ,{ 42,  PB042,  "162"}
    ,{ 43,  PB043,  "16695334890" , "Sub-string divisibility"}
    ,{ 44,  PB044,  "5482660" ,"Pentagon numbers"}
    ,{ 45,  PB045,  "1533776805" ,"Triangular, pentagonal, and hexagonal"}
    ,{ 46,  PB046,  "5777" ,    "Goldbach's other conjecture" }
    ,{ 47,  PB047,  "134043" ,  "Distinct primes factors" }
    ,{ 48,  PB048,  "9110846700" ,  "Self powers" }
    ,{ 49,  PB049,  "296962999629" ,"Prime permutations" }
    ,{ 50,  PB050,  "997651" ,  "Consecutive prime sum" }
    ,{ 51,  PB051,  "121313" ,   "Prime digit replacements" }
    ,{ 52,  PB052,  "142857", "Permuted multiples"}
    ,{ 52,  PB052a,  "142857", "Permuted multiples"}
    ,{ 53,  PB053,  "4075", "Combinatoric selections" }
    ,{ 55,  PB055,  "249" , "Lychrel numbers" }
    ,{ 56,  PB056_gmp,  "972",  "Powerful digit sum" }
    ,{ 57,  PB057,  "153",  "Square root convergents" }
    ,{ 58,  PB058,  "26241" ,   "Spiral primes" }
    ,{ 59,  PB059,  "107359" ,  "XOR decryption" }
    ,{ 60,  PB060,  "26033" ,  "Prime pair sets" }
    ,{ 61,  PB061,  "28684" , "Cyclical figurate numbers"}
    ,{ 62,  PB062,  "127035954683" , "Cubic permutations"}
    ,{ 63,  PB063,  "49" , "Powerful digit counts"}
    ,{ 64,  PB064,  "1322" , "Odd period square roots"}
    ,{ 65,  PB065,  "272" , "Convergents of e"}
    ,{ 66,  PB066,  "661", "Diophantine equation"}
    ,{ 67,  PB067,  "7273", "Maximum path sum II"}
    ,{ 68,  PB068,  "6531031914842725", "Magic 5-gon ring"}
//    ,{ 69,  PB069a,  "510510", "Totient maximum"}
    ,{ 69,  PB069b,  "510510", "Totient maximum"}
//    ,{ 70,  PB070,  "8319823", "Totient permutation"}
    ,{ 70,  PB070a,  "8319823", "Totient permutation"}
    ,{ 71,  PB071,  "428570", "Ordered fractions"}
    ,{ 72,  PB072,  "303963552391", "Counting fractions"}
    ,{ 73,  PB073,  "7295372", "Counting fractions in a range"}
//    ,{ 74,  PB074,  "402", "Digit factorial chains"}
    ,{ 74,  PB074a,  "402", "Digit factorial chains"}
    ,{ 75,  PB075,  "161667", "Singular integer right triangles"}
    ,{ 76,  PB076,  "190569291", "Counting summations"}
    ,{ 76,  PB076a,  "190569291", "Counting summations"}
    ,{ 77,  PB077,  "71", "Prime summations"}
    ,{ 78,  PB078,  "55374", "Counting summations"}
    ,{ 80,  PB080_gmp,  "40886", "Square root digital expansion"}
    ,{ 81,  PB081,  "427337", "Path sum: two ways"}
    ,{ 82,  PB082,  "260324", "Path sum: three ways"}
    ,{ 83,  PB083,  "425185", "Path sum: four ways"}
    ,{ 85,  PB085,  "2772", "Counting rectangles"}
    ,{ 86,  PB086,  "1818", "Cuboid route"}
    ,{ 87,  PB087,  "1097343", "Prime power triples"}
    ,{ 88,  PB088,  "7587457", "Product-sum numbers"}
    ,{ 90,  PB090,  "1217", "Cube digit pairs"}
    ,{ 91,  PB091,  "14234", "Right triangles with integer coordinates"}
//    ,{ 92,  PB092,  "8581146", "Square digit chains"}
    ,{ 92,  PB092a,  "8581146", "Square digit chains"}
    ,{ 93,  PB093,  "1258", "Arithmetic expressions"}
// version recursive parametrable en nombre de digits
//    ,{ 93,  PB093a,  "1258", "Arithmetic expressions"},
    ,{ 94,  PB094,  "518408346", "Almost equilateral triangles"}
    ,{ 95,  PB095,  "14316", "Amicable chains"}

    ,{ 579,  PB579,  "3805524", "Lattice points in lattice cubes"}

// a revoir beaucoup trop lent
//    ,{100,  PB100,  "756872327473"}
//    ,{ 1000,  PB1000,  "179424673", "Test 10 000 000 prime numbers"}
// version moins rapide pour grande valeurs
//    ,{ 597,  PB597_gmpa,  "50018178282", "Torpids"}
// version rapide mais calculs complets (pour debug)
//    ,{ 597,  PB597_gmp,  "50018178282", "Torpids"}
// version la plus rapide avec calcul FACT en u_int64_t et optimisation last loop
//    ,{ 597,  PB597_gmpx,  "50018178282", "Torpids"}
// version beaucoup plus rapide en O(n**2)
    ,{ 597,  PB597_gmpy,  "50018178282", "Torpids"}
    ,{ 622,  PB622,  "3010983666182123972", "Riffle Shuffles"}
    
    ,{  0,  NULL,   ""}
};


static PB_CALL CUR_calls[] = {
 //    { 51,  PB051,  "121313" ,   "Prime digit replacements" },
 //   {100,  PB100,  "756872327473"},
    { 101,  PB101,  "37076114526", "Optimum polynomial"},

    {  0,  NULL,   ""}
} ;


int main(int argc, const char * argv[]) {
    PB_CALL *ptCall ;
    TotalRun ttr ;
    

    ttr.isVerbose = 1 ;
    ttr.nbPBerror  = ttr.nbPBOK = 0 ;
    ttr.TotalClock = 0 ;
//    int isALL = 100 ;
    int isALL = 0 ;
    int pbMax = (isALL == 0) ? 1000 : isALL ;
    for(ptCall = (isALL==0) ? CUR_calls : ALL_calls ; ptCall->pbNum != 0 && ptCall->pbNum < pbMax ; ptCall++) {
        Execute(&ttr,ptCall);
    }
    fprintf(stdout,"### Execution of %d PB en %.06fs (%d en erreur)\n",ttr.nbPBOK+ttr.nbPBerror,(float) ttr.TotalClock / CLOCKS_PER_SEC,ttr.nbPBerror) ;
    return 0;
}
