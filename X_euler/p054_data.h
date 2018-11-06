//
//  p054_data.h
//  X_euler
//
//  Created by Jeannot on 18/04/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef p054_data_h
#define p054_data_h

#include <stdio.h>
#include <stdint.h>
typedef struct Card {
    uint8_t h ;    // high
    uint8_t c ;    // color
    
} Card ;

typedef struct Play {
    Card Hand0[5] ;
    Card Hand1[5] ;
    
} Play ;

const Play *  P054_GetData(void) ;
#endif /* p054_data_h */
