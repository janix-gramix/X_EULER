//
//  PB151_200.c
//  X_euler
//
//  Created by Jeannot on 01/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "PB151_200.h"

int PB198(PB_RESULT *pbR) {
    pbR->nbClock = clock() ;
    pbR->nbClock = clock() - pbR->nbClock ;
//    sprintf(pbR->strRes,"%d",num) ;
    return 1 ;
}
