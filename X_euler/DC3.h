//
//  DC3.h
//  X_euler
//
//  Created by Jeannot on 17/12/2019.
//  Copyright © 2019 Jeannot. All rights reserved.
//

#ifndef DC3_h
#define DC3_h
#include <stdlib.h>

typedef int32_t DC3_TYPE ;
typedef int     DC3_INDEX ;

void RadixSort(const DC3_INDEX * inInd, DC3_INDEX * outInd,const DC3_TYPE *T, int n,DC3_TYPE nbType ) ;
DC3_INDEX * CountRadixSort_uint16(const DC3_INDEX * inInd, DC3_INDEX * outInd,const uint16_t *T, int n,DC3_TYPE nbType );
void SuffixSort(const DC3_TYPE * T, DC3_INDEX * SA, int n, DC3_TYPE nbType) ;
DC3_INDEX * GetLCP(const uint8_t * text,const DC3_INDEX * Sindex, int nbIndex,int len) ;

#endif /* DC3_h */
