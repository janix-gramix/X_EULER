//
//  p102_data.h
//  X_euler
//
//  Created by Jeannot on 05/04/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef p102_data_h
#define p102_data_h

#include <stdio.h>
#include <stdint.h>
#define PB102_NBT   1000
typedef struct Triangle {
    int32_t xa;
    int32_t ya;
    int32_t xb;
    int32_t yb;
    int32_t xc;
    int32_t yc;
    
} Triangle ;

const Triangle *  P102_GetData(void) ;

#endif /* p102_data_h */
