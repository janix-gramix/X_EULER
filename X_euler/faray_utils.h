//
//  faray_utils.h
//  X_euler
//
//  Created by Jeannot on 03/10/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#ifndef faray_utils_h
#define faray_utils_h
#include <stdint.h>
typedef struct FRACTRED {
    int n ;
    int d ;
} FRACTRED ;

FRACTRED Besout(FRACTRED fr0) ; // return fr1 tel que fr0<fr1 et fr1.n*fr0.d-fr1.d*fr0.n=1
FRACTRED IBesout(FRACTRED fr1) ; // return fr0 tel que fr0<fr1 et fr1.n*fr0.d-fr1.d*fr0.n=1
// Implementation of Stern-Brocot Tree with stack
//

typedef struct SBTree {
    FRACTRED * stack ;
    int sizeStack ;
    int indS ;
    FRACTRED fr0 ;
    FRACTRED fr1 ;
} SBTree ;

SBTree * SBT_alloc(void) ;
void SBT_init(SBTree * sbt,FRACTRED fr0, FRACTRED fr1) ;

// obtient (dans sbt->fr1) la fraction suivante entre fr0 et fr1 d'init.
// retourne 0 quand le parcours est termine
int SBT_getNext(SBTree * sbt,int maxDen) ;

// is isOK=1 returne ds fr1 (et ds stack) la nouvelle fraction entre fr0 et fr1 . Si isOK=0 fr1->fr0, on depile fr1
int SBT_ValidNxt(SBTree * sbt, int isOK) ;

SBTree * SBT_free(SBTree * sbt) ;

// meme serie de fonction en gerant seulement les denominateurs
typedef struct SBdTree {
    int32_t * stack ;
    int sizeStack ;
    int indS ;
    int32_t d0 ;
    int32_t d1 ;
} SBdTree ;

typedef int (* STB_CB)(FRACTRED fr0,FRACTRED fr1) ;
typedef int (* STBDen_CB)(int d0,int d1) ;


void SBdT_init(SBdTree * sbdt,int32_t d0, int32_t d1) ;
SBdTree * SBdT_alloc(void) ;
int SBdT_ValidNxt(SBdTree * sbdt, int isOK) ;

SBdTree * SBdT_free(SBdTree * sbdt) ;
    
// implementation of Stern-Brocot by recursive fonction with CB
int STBrcv(FRACTRED fr0, FRACTRED fr1,STB_CB stbRcvCB) ;

// idem with denominator only
int STBrcvDen(int d0, int d1,STBDen_CB stbRcvDenCB ) ;

typedef struct FRACTrec {
    FRACTRED fr0 ;
    FRACTRED fr1 ;
    FRACTRED frEnd ;
    int maxDen ;
} FRACTrec ;
void FRC_init(FRACTrec *FRrec,int maxDen , FRACTRED frStart, FRACTRED frEnd) ;
// obtient (dans FRrec->fr1) la fraction suivante entre frStart et frEnd
// retourne 0 quand le parcours est termine
int FRC_getNext(FRACTrec *FRrec) ;


#endif /* faray_utils_h */
