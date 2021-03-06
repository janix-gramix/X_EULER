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
#include "PB151_200.h"
#include "PB201_250.h"
#include "PB_gmp.h"
#include "PB_other.h"
#include "pb_debug.h"





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
        pbR.ident = pbCall->ident ;
        pbR.nbClock = 0 ;
        if(pbCall->pbSolve(&pbR)) {
            ttr->TotalClock += pbR.nbClock ;
            if(strcmp(pbCall->Solution,pbR.strRes)==0) {
                ttr->nbPBOK++ ;
                fprintf(stdout,"OK\tPB%s(%.06fs) Sol=%s %s\n",pbCall->ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,(pbCall->name) ? pbCall->name : "") ;
            } else {
                ttr->nbPBerror++ ;
                fprintf(stdout,"ERROR\tPB%s(%.06fs) Find=%s != Exp=%s\n",pbCall->ident,(float)pbR.nbClock / CLOCKS_PER_SEC,pbR.strRes,pbCall->Solution) ;
            }
        } else {
            ttr->TotalClock += pbR.nbClock ;
            ttr->nbPBerror++ ;
            fprintf(stdout," BAD EXECUTION PB%s(%.06fs)\n",pbCall->ident,(float)pbR.nbClock / CLOCKS_PER_SEC) ;
        }
        
    } else {
        fprintf(stdout,"SKIPPED\tPB%s\n",pbCall->ident) ;
    }

    
}


#define PROTO_CALL(NUM,SOL,NAME)   {#NUM, PB##NUM , #SOL, #NAME}

/*
 Liens
 173->174->179
 76->181 (decompostion en partition et nombre pentagonaux)
 64 -> 94 -> 138 (fraction continue (Sqrt(n) n=3 pour 94, 138)
 
*/
static PB_CALL ALL_calls[] = {
    

    PROTO_CALL(001,233168,'Multiples of 3 and 5')
    ,PROTO_CALL(002,4613732,'Even Fibonacci numbers')
    ,PROTO_CALL(003,6857,'Largest prime factor')
    ,PROTO_CALL(004,906609,'Largest palindrome product')
    ,PROTO_CALL(005,232792560,'Smallest multipl')
    ,PROTO_CALL(006,25164150,'Smallest multiple')
    ,PROTO_CALL(007,104743,'Smallest multiple')
    ,PROTO_CALL(008,23514624000,'Largest product in a series')
    ,PROTO_CALL(009,31875000,'Special Pythagorean triplet')
    ,PROTO_CALL(010,142913828922,'Summation of primes')
    ,PROTO_CALL(011,70600674,'Largest product in a grid')
    ,PROTO_CALL(012,76576500,'Highly divisible triangular number')
    ,PROTO_CALL(013,5537376230,'Large sum')
    ,PROTO_CALL(014,837799,'Longest Collatz sequence')
    ,PROTO_CALL(015,137846528820,'Lattice paths')
    ,PROTO_CALL(016,1366,'Power digit sum')
    ,PROTO_CALL(016_gmp,1366,'Power digit sum')
    ,PROTO_CALL(017,21124,'Number letter counts')
    ,PROTO_CALL(018,1074,'Maximum path sum I')
    ,PROTO_CALL(019,171,'Counting Sundays')

    ,PROTO_CALL(020,648,'Factorial digit sum')
    ,PROTO_CALL(021,31626,'Amicable numbers')
    ,PROTO_CALL(022,871198282,'Names scores')

    ,PROTO_CALL(023,4179871,'Non-abundant sums')
    ,PROTO_CALL(024,2783915460,'Lexicographic permutations')
    ,PROTO_CALL(025,4782,'1000-digit Fibonacci number')
    ,PROTO_CALL(026,983,'Reciprocal cycles')
    ,PROTO_CALL(027,-59231,'Quadratic primes')
    ,PROTO_CALL(028,669171001,'Number spiral diagonals')
    ,PROTO_CALL(029,9183,'Distinct powers')
    ,PROTO_CALL(030,443839,'Digit fifth powers')
    ,PROTO_CALL(031,73682,'Coin sums')
    ,PROTO_CALL(032,45228,'Pandigital products')
    ,PROTO_CALL(033,100,'Digit cancelling fractions')
    ,PROTO_CALL(034,40730,'Digit factorials')
    ,PROTO_CALL(035,55,'Circular primes')
    ,PROTO_CALL(036,872187,'Double-base palindromes')
    ,PROTO_CALL(037,748317,'Truncatable primes')
    ,PROTO_CALL(037r,748317,'Truncatable primes')
    ,PROTO_CALL(038,932718654,'Truncatable primes')
    ,PROTO_CALL(039,840,'Integer right triangles')
    ,PROTO_CALL(040,210,'Champernowne\'s constant')
    ,PROTO_CALL(041,7652413,'Pandigital prime')
    ,PROTO_CALL(042,162,'Coded triangle numbers')
    ,PROTO_CALL(043,16695334890,'Sub-string divisibility')
    ,PROTO_CALL(044,5482660,'Pentagon numbers')
    ,PROTO_CALL(045,1533776805,'Triangular, pentagonal, and hexagonal')
    ,PROTO_CALL(046,5777,'Goldbach\'s other conjecture')
    ,PROTO_CALL(047,134043,'Distinct primes factors')
    ,PROTO_CALL(048,9110846700,'Self powers')
    ,PROTO_CALL(049,296962999629,'Prime permutations')
    ,PROTO_CALL(050,997651,'Consecutive prime sum')
    ,PROTO_CALL(051,121313,'Prime digit replacements')
    ,PROTO_CALL(052,142857,'Permuted multiples')
    ,PROTO_CALL(052a,142857,'Permuted multiples')
    ,PROTO_CALL(053,4075,'Combinatoric selections')
    ,PROTO_CALL(054,376,'Poker hands')
    ,PROTO_CALL(055,249,'Lychrel numbers')
    ,PROTO_CALL(056_gmp,972,'Powerful digit sum')
    ,PROTO_CALL(057,153,'Square root convergents')
    ,PROTO_CALL(058,26241,'Spiral primes')
    ,PROTO_CALL(059,107359,'XOR decryption')
    ,PROTO_CALL(060,26033,'Prime pair sets')
    ,PROTO_CALL(061,28684,'Cyclical figurate numbers')
    ,PROTO_CALL(062,127035954683,'Cubic permutations')
    ,PROTO_CALL(063,49,'Powerful digit counts')
    ,PROTO_CALL(064,1322,'Odd period square roots')
    ,PROTO_CALL(065,272,'Convergents of e')
    ,PROTO_CALL(066_gmp,661,'Diophantine equation')
    ,PROTO_CALL(067,7273,'Maximum path sum II')
    ,PROTO_CALL(068,6531031914842725,'Magic 5-gon ring')
//    ,PROTO_CALL(069a,510510,'Totient maximum')
    ,PROTO_CALL(069b,510510,'Totient maximum')
//    ,PROTO_CALL(070,8319823,'Totient permutation')
    ,PROTO_CALL(070a,8319823,'Totient permutation')
    ,PROTO_CALL(071,428570,'Ordered fractions')
    ,PROTO_CALL(072,303963552391,'Counting fractions')
    ,PROTO_CALL(073,7295372,'Counting fractions in a range')
//    ,PROTO_CALL(073a,7295372,'Counting fractions in a range')
//    ,PROTO_CALL(073b,7295372,'Counting fractions in a range')
//    ,PROTO_CALL(073c,7295372,'Counting fractions in a range')

//    ,PROTO_CALL(074,402,'Digit factorial chains')
    ,PROTO_CALL(074a,402,'Digit factorial chains')
    ,PROTO_CALL(075,161667,'Singular integer right triangles')
    ,PROTO_CALL(076,190569291,'Counting summations')
    ,PROTO_CALL(076a,190569291,'Counting summations')
    ,PROTO_CALL(077,71,'Prime summations')
    ,PROTO_CALL(078,55374,'Counting summations')
    ,PROTO_CALL(080_gmp,40886,'Square root  <digital expansion')
    ,PROTO_CALL(081,427337,'Path sum: two ways')
    ,PROTO_CALL(082,260324,'Path sum: three ways')
    ,PROTO_CALL(083,425185,'Path sum: four ways')
    //,PROTO_CALL(084,101524,'Monopoly odds')
    //,PROTO_CALL(084a,101524,'Monopoly odds')
    ,PROTO_CALL(084b,101524,'Monopoly odds')
    ,PROTO_CALL(085,2772,'Counting rectangles')
    ,PROTO_CALL(086,1818,'Cuboid route')
    ,PROTO_CALL(087,1097343,'Prime power triples')
    ,PROTO_CALL(088,7587457,'Product-sum numbers')
    ,PROTO_CALL(089,743,'Roman numerals')
    ,PROTO_CALL(090,1217,'Cube digit pairs')
    ,PROTO_CALL(091,14234,'Right triangles with integer coordinates')
//    ,PROTO_CALL(092,8581146,'Square digit chains')
    ,PROTO_CALL(092a,8581146,'Square digit chains')
    ,PROTO_CALL(093,1258,'Arithmetic expressions')
// version recursive parametrable en nombre de digits
//    ,PROTO_CALL(093a,1258,'Arithmetic expressions'),
    ,PROTO_CALL(094,518408346,'Almost equilateral triangles')
    ,PROTO_CALL(095,14316,'Amicable chains')
    ,PROTO_CALL(097,8739992577,'Large non-Mersenne prime')

    ,PROTO_CALL(098,18769,'Anagramic squares')
    ,PROTO_CALL(099,709,'Largest exponential')

    ,PROTO_CALL(101,37076114526,'Optimum polynomial')
    ,PROTO_CALL(102,228,'Triangle containment')

    ,PROTO_CALL(103g,20313839404245,'Special subset sums: optimum')
// diverse variante de moins en moins efficaces surtout pour des ordres superieurs (8,9,10 PB NP-complet)
//    ,PROTO_CALL(103f,20313839404245,'Special subset sums: optimum')
//    ,PROTO_CALL(103e,20313839404245,'Special subset sums: optimum')
//    ,PROTO_CALL(103c,20313839404245,'Special subset sums: optimum')
//    ,PROTO_CALL(103b,20313839404245,'Special subset sums: optimum')
//    ,PROTO_CALL(103a,20313839404245,'Special subset sums: optimum')
//    ,PROTO_CALL(103,20313839404245,'Special subset sums: optimum')

    
//    PROTO_CALL(104_gmp,329468,'Pandigital Fibonacci ends')
    ,PROTO_CALL(104,329468,'Pandigital Fibonacci ends')
    ,PROTO_CALL(105,73702,'Special subset sums: testing')
    ,PROTO_CALL(106,21384,'Special subset sums: meta-testing')
    ,PROTO_CALL(107,259679,'Minimal network')
     ,{"108a",PB110,"180180","Diophantine reciprocals I"}                // cas /    ,PROTO_CALL(108,180180,'Diophantine reciprocals I')
    ,PROTO_CALL(109,38182,'Darts')
    ,PROTO_CALL(110,9350130049860600,'Diophantine reciprocals II')
    ,PROTO_CALL(111,612407567715,'Primes with runs')
    ,PROTO_CALL(112,1587000,'Bouncy numbers')
    ,PROTO_CALL(113,51161058134250,'Non-bouncy numbers')
    ,PROTO_CALL(114,16475640049,'Counting block combinations I')
//    ,PROTO_CALL(114a,16475640049,'Counting block combinations I')
    ,PROTO_CALL(115,168,'Counting block combinations II')
    ,PROTO_CALL(116,20492570929,'Red, green or blue tiles')
    ,PROTO_CALL(117,100808458960497,'Red, green, and blue tiles')
    ,PROTO_CALL(118,44680,'Pandigital prime sets')
    ,PROTO_CALL(119,248155780267521,'Digit power sum')
    ,PROTO_CALL(120,333082500,'Square remainders')
    ,PROTO_CALL(121,2269,'Disc game prize fund')
    ,PROTO_CALL(123,21035,'Prime square remainders')
    ,PROTO_CALL(124,21417,'Ordered radicals')
    ,PROTO_CALL(125,2906969179,'Palindromic sums')
    ,PROTO_CALL(126,18522,'Cuboid layers')
//    ,PROTO_CALL(127,18407904,'abc-hits')
    ,PROTO_CALL(127a,18407904,'abc-hits')
//   ,PROTO_CALL(128,14516824220,'Hexagonal tile differences')
    ,PROTO_CALL(128a,14516824220,'Hexagonal tile differences')

    ,PROTO_CALL(129,1000023,'Repunit divisibility')
    ,PROTO_CALL(130,149253,'Composites with prime repunit property')
    ,PROTO_CALL(131,173,'Prime cube partnership')
    ,PROTO_CALL(132,843296,'Large repunit factors')
//    ,PROTO_CALL(133,453647705,'Repunit nonfactors')
    ,PROTO_CALL(133a,453647705,'Repunit nonfactors')
//    ,PROTO_CALL(134,18613426663617118,'Prime pair connection')
    ,PROTO_CALL(134a,18613426663617118,'Prime pair connection')
//    ,PROTO_CALL(135,4989,'Same differences')
    ,PROTO_CALL(135a,4989,'Same differences')
    ,PROTO_CALL(136,2544559,'Singleton difference')
    ,PROTO_CALL(137,1120149658760,'Fibonacci golden nuggets')

    ,PROTO_CALL(138,1118049290473932,'Special isosceles triangles')

//  ,PROTO_CALL(139,10057761,'Pythagorean tiles')
    ,PROTO_CALL(139a,10057761,'Pythagorean tiles')
    ,PROTO_CALL(140,5673835352990,'Modified Fibonacci golden nuggets')

//    ,PROTO_CALL(141c,878454337159,'Investigating progressive numbers, n, which are also square')
   ,PROTO_CALL(141b,878454337159,'Investigating progressive numbers, n, which are also square')
//    ,PROTO_CALL(141,878454337159,'Investigating progressive numbers, n, which are also square')
    ,PROTO_CALL(142,1006193,'Perfect Square Collection')

//    ,PROTO_CALL(143,30758397,'Investigating the Torricelli point of a triangle')
    ,PROTO_CALL(143a,30758397,'Investigating the Torricelli point of a triangle')
//    ,PROTO_CALL(144,354,'Investigating multiple reflections of a laser beam')
    ,PROTO_CALL(144a,354,'Investigating multiple reflections of a laser beam')

    ,PROTO_CALL(145a,608720,'How many reversible numbers are there below one-billion?')
 // ,PROTO_CALL(145,608720,'How many reversible numbers are there below one-billion?')
    ,PROTO_CALL(146a,676333270,'Investigating a Prime Pattern')
 //   ,PROTO_CALL(146,676333270,'Investigating a Prime Pattern')

    ,PROTO_CALL(147,846910284,'Rectangles in cross-hatched grids')
    ,PROTO_CALL(148,2129970655314432,'Exploring Pascal''s triangle')
    ,PROTO_CALL(149,52852124,'Searching for a maximum-sum subsequence')
    ,PROTO_CALL(150,-271248680,'Searching a triangular array for a sub-triangle having minimum-sum')

//    ,PROTO_CALL(151,0.464399,'Paper sheets of standard sizes: an expected-value problem')
    ,PROTO_CALL(151a,0.464399,'Paper sheets of standard sizes: an expected-value problem')

    ,PROTO_CALL(152c,301,'Writing 1/2 as a sum of inverse squares')
//    ,PROTO_CALL(152a,301,'Writing 1/2 as a sum of inverse squares')
//    ,PROTO_CALL(152,301,'Writing 1/2 as a sum of inverse squares')
    ,PROTO_CALL(153c,17971254122360635,'Investigating Gaussian Integers')
//    ,PROTO_CALL(153b,17971254122360635,'Investigating Gaussian Integers')
//    ,PROTO_CALL(153a,17971254122360635,'Investigating Gaussian Integers')
//    ,PROTO_CALL(153,17971254122360635,'Investigating Gaussian Integers')

     ,PROTO_CALL(154a,479742450,'Exploring Pascal''s pyramid')
//    ,PROTO_CALL(154,479742450,'Exploring Pascal''s pyramid')

    ,PROTO_CALL(155,3857447,'Counting Capacitor Circuits')
    ,PROTO_CALL(156,21295121502550,'Counting Digits')
    ,PROTO_CALL(157,53490,'Solving the diophantine equation 1/a+1/b= p/10n')
    ,PROTO_CALL(158,409511334375,'Exploring strings for which only one character comes lexicographically after its neighbour to the left')
//    ,PROTO_CALL(159,14489159,'Digital root sums of factorisations')
//    ,PROTO_CALL(159a,14489159,'Digital root sums of factorisations')
    ,PROTO_CALL(159b,14489159,'Digital root sums of factorisations')
    ,PROTO_CALL(160,16576,'Factorial trailing digits')
    ,PROTO_CALL(161,20574308184277971,'Triominoes')
    ,PROTO_CALL(162,3D58725572C62302,'Hexadecimal numbers')

    ,PROTO_CALL(164,378158756814587,'Numbers for which no three consecutive digits have a sum greater than a given value')
    ,PROTO_CALL(165,2868868,'Intersections')
    ,PROTO_CALL(166,7130034,'Criss Cross')
//    ,PROTO_CALL(167,3916160068885,'Investigating Ulam sequences')
    ,PROTO_CALL(167a,3916160068885,'Investigating Ulam sequences')
//    ,PROTO_CALL(168,59206,'Number Rotations')
    ,PROTO_CALL(168a,59206,'Number Rotations')

    ,PROTO_CALL(169b,178653872807,'Exploring the number of different ways a number can be expressed as a sum of powers of 2')
//   ,PROTO_CALL(169a,178653872807,'Exploring the number of different ways a number can be expressed as a sum of powers of 2')
//    ,PROTO_CALL(169,178653872807,'Exploring the number of different ways a number can be expressed as a sum of powers of 2')

    ,PROTO_CALL(173,1572729,'Using up to one million tiles how many different "hollow" square laminae can be formed?')
    ,PROTO_CALL(174,209566,'Counting the number of "hollow" square laminae that can form one, two, three, ... distinct arrangements')
    ,PROTO_CALL(175,'1,13717420,8','Fractions involving the number of different ways a number can be expressed as a sum of powers of 2')

    ,PROTO_CALL(179c,986262,'Consecutive positive divisors')
//    ,PROTO_CALL(179b,986262,'Consecutive positive divisors')
//    ,PROTO_CALL(179a,986262,'Consecutive positive divisors')
//    ,PROTO_CALL(179,986262,'Consecutive positive divisors')

    ,PROTO_CALL(181,83735848679360680,'Investigating in how many ways objects of two different colours can be grouped')
//    ,PROTO_CALL(181a,83735848679360680,'Investigating in how many ways objects of two different colours can be grouped')
    ,PROTO_CALL(182,399788195976,'RSA encryption')
    ,PROTO_CALL(183,48861552,'Maximum product of parts')

//    ,PROTO_CALL(184a,1725323624056,'Triangles containing the origin')
    ,PROTO_CALL(184b,1725323624056,'Triangles containing the origin')
//    ,PROTO_CALL(184,1725323624056,'Triangles containing the origin')
    ,PROTO_CALL(185,4640261571849533,'Number Mind')

//    ,PROTO_CALL(187,17427258,'Semiprimes')
    ,PROTO_CALL(187a,17427258,'Semiprimes')
//    ,PROTO_CALL(187b,17427258,'Semiprimes')
    
    ,PROTO_CALL(188,95962097,'The hyperexponentiation of a number')
//    ,PROTO_CALL(188_gmp,95962097,'The hyperexponentiation of a number')

//    ,PROTO_CALL(191,1918080160,'Prize Strings')
    ,PROTO_CALL(191a,1918080160,'Prize Strings')

    
    ,PROTO_CALL(192,57060635927998347,'Best Approximations')
//    ,PROTO_CALL(192_gmp,57060635927998347,'Best Approximations')
    ,PROTO_CALL(193,684465067343069,'Squarefree Numbers')
//    ,PROTO_CALL(193a,684465067343069,'Squarefree Numbers')
    ,PROTO_CALL(195a,75085391,'Inscribed circles of triangles with one angle of 60 degrees')
//   ,PROTO_CALL(195,75085391,'Inscribed circles of triangles with one angle of 60 degrees')


//    ,PROTO_CALL(198,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198a,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198b,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198c,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198d,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198e,52374425,'Ambiguous Numbers')
//    ,PROTO_CALL(198f,52374425,'Ambiguous Numbers')
    ,PROTO_CALL(198g,52374425,'Ambiguous Numbers')
    ,PROTO_CALL(199,0.00396087,'Iterative Circle Packing')

    ,PROTO_CALL(204b,2944730,'Generalised Hamming Numbers')
//    ,    PROTO_CALL(204a,2944730,'Generalised Hamming Numbers')
//    ,    PROTO_CALL(204,2944730,'Generalised Hamming Numbers')


    ,PROTO_CALL(206,1389019170,'Concealed Square')
    ,PROTO_CALL(206a,1389019170,'Concealed Square')
    
    ,PROTO_CALL(234,1259187438574927161,'Semidivisible numbers')
    ,PROTO_CALL(235,1.002322108633,'An Arithmetic Geometric sequence')
    ,PROTO_CALL(345a,13938,'Matrix Sum')

    // a revoir beaucoup trop lent
    //    ,{100,  PB100,756872327473,'XX')
    //    ,PROTO_CALL(1000,179424673,'Test 10 000 000 prime numbers')
    // version moins rapide pour grande valeurs
    //    ,PROTO_CALL(597_gmpa,50018178282,'Torpids')
    // version rapide mais calculs complets (pour debug)
    //    ,PROTO_CALL(597_gmp,50018178282',Torpids')
    // version la plus rapide avec calcul FACT en uint64_t et optimisation last loop
    //    ,PROTO_CALL(597_gmpx,50018178282,'Torpids')
    // version beaucoup plus rapide en O(n**2)
    
    ,PROTO_CALL(357c,1739023853137,'Prime generating integers')
//    ,PROTO_CALL(357b,1739023853137,'Prime generating integers')
//    ,PROTO_CALL(357a,1739023853137,'Prime generating integers')
//    ,PROTO_CALL(357,1739023853137,'Prime generating integers')

    ,PROTO_CALL(500,35407281,'Problem 500!!!')
    ,PROTO_CALL(540a,500000000002845,'Counting primitive Pythagorean triples')
    ,PROTO_CALL(540,500000000002845,'Counting primitive Pythagorean triples')

    ,PROTO_CALL(579,3805524,'Lattice points in lattice cubes')
    ,PROTO_CALL(597_gmpy,50018178282,'Torpids')
    ,PROTO_CALL(620a,1470337306,'Planetary Gears')
//    ,PROTO_CALL(620,1470337306,'Planetary Gears') //slow version
    ,PROTO_CALL(622,3010983666182123972,'Riffle Shuffles')
    ,PROTO_CALL(625a,551614306,'Gcd sum')
//    ,PROTO_CALL(625b,551614306,'Gcd sum') // too slow
//    ,PROTO_CALL(625,551614306,'Gcd sum') // verification for small values
//    ,PROTO_CALL(687a,0.3285320869,'Shuffling Cards')
//    ,PROTO_CALL(687,0.3285320869,'Shuffling Cards')

    //    ,PROTO_CALL(691c,11570761,'Long substring with many repetitions')
    ,PROTO_CALL(691d,11570761,'Long substring with many repetitions')
 //   ,PROTO_CALL(691a,11570761,'Long substring with many repetitions')
 //   ,PROTO_CALL(691b,11570761,'Long substring with many repetitions')
 //   ,PROTO_CALL(691,11570761,'Long substring with many repetitions')
     ,PROTO_CALL(693,699161,'Finite Sequence Generator')
//   ,PROTO_CALL(693b,699161,'Finite Sequence Generator')
//   ,PROTO_CALL(693a,699161,'Finite Sequence Generator')
    ,PROTO_CALL(696,436944244,'Mahjong')
    ,PROTO_CALL(697,4343871.06,'Randomly Decaying Sequence')
//    ,PROTO_CALL(697a,4343871.06,'Randomly Decaying Sequence')
    ,PROTO_CALL(701,13.51099836,'Random connected area')
    ,PROTO_CALL(705,480440153,'Total Inversion Count of Divided Sequences')
//    ,PROTO_CALL(705a,480440153,'Total Inversion Count of Divided Sequences')
    ,PROTO_CALL(708a,28874142998632109,'Twos are all you need')
//    ,PROTO_CALL(708b,28874142998632109,'Twos are all you need')
//    ,PROTO_CALL(708c,28874142998632109,'Twos are all you need')
//    ,PROTO_CALL(708,28874142998632109,'Twos are all you need')


    
    ,{ NULL,NULL,"","" }
};


static PB_CALL CUR_calls[] = {
 //    { 051,121313,'Prime digit replacements'),
 //   {100,756872327473,'XX'),
//    PROTO_CALL(103d,20313839404245,'Special subset sums: optimum')

//        PROTO_CALL(695,???,'Random Rectangles')

         PROTO_CALL(702c,???,'Jumping Flea')
        ,PROTO_CALL(702a,???,'Jumping Flea')
//        PROTO_CALL(707,???,'Lights Out')
  
    //        ,PROTO_CALL(690,???,'Tom and Jerry')
  //    ,PROTO_CALL(681a,2611227421428,'Maximal Area')
//   ,PROTO_CALL(681,2611227421428,'Maximal Area')
 //    PROTO_CALL(495,???,'Writing n as the product of k distinct positive integers')



    
//    PROTO_CALL(626,???,'Counting Binary Matrices')

//     ,PROTO_CALL(625,???,'Gcd sum')
//    ,PROTO_CALL(107,259679,'Minimal network')
    

    ,{NULL,NULL,"",""}
} ;


int main(int argc, const char * argv[]) {
    PB_CALL *ptCall ;
    TotalRun ttr ;
    

    ttr.isVerbose = 1 ;
    ttr.nbPBerror  = ttr.nbPBOK = 0 ;
    ttr.TotalClock = 0 ;
    clock_t debut = clock() ;
//    char * isALL = "151" ;
    char * isALL = NULL ;
    char  * pbMax = (isALL == NULL) ? "ZZZZ" : isALL ;
    for(ptCall = (isALL==NULL) ? CUR_calls : ALL_calls ; ptCall->ident != NULL && strcmp(ptCall->ident,pbMax) < 0 ; ptCall++) {
        Execute(&ttr,ptCall);
    }
    debut = clock() - debut ;
    fprintf(stdout,"### Execution of %d PB en %.03fs [%.03fs] (%d en erreur)\n",ttr.nbPBOK+ttr.nbPBerror,(float) ttr.TotalClock / CLOCKS_PER_SEC,(float) debut / CLOCKS_PER_SEC,ttr.nbPBerror) ;
    return 0;
}
