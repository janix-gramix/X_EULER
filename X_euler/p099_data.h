//
//  p099_data.h
//  X_euler
//
//  Created by Jeannot on 18/09/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef p099_data_h
#define p099_data_h
#define p099_size   1000

typedef struct P099_couple {
    int    val ;
    int    exp ;
} P099_couple ;

const P099_couple *  P099_GetData(void) ;

#endif /* p099_data_h */
