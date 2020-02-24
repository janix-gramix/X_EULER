//
//  ukkonen.h
//  X_euler
//
//  Created by Jeannot on 17/12/2019.
//

#ifndef ukkonen_h
#define ukkonen_h

#include <stdlib.h>

typedef struct STree STree ;

// build suffixTree
// character from [0..length[ must be in ]0,maxChar[
// suit must be terminated by 0 (suit[length]=0 )
// for a suit of 1,2 terminated by a 0, maxChar=3
STree * X_CompSuffixTree(const uint8_t *suit,int maxChar,int length) ;

// build suffixTree
// character from [0..length[ must be in [0,maxChar[
// For a suit of 0,1 maxChar=2
STree *  ComputeSuffixTree(const uint8_t *suit,int maxChar,int length) ;

STree *  FreeSuffixTree(STree * ST) ;

void CompSALCP(STree *ST,int32_t *sa,int32_t *lcp) ;
int32_t NbNodeUsed(STree * ST) ;
#endif /* ukkonen_h */
