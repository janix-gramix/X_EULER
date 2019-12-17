//
//  ukkonen.c
//  X_euler
//
//  Created by Jeannot on 17/12/2019.
//  Copyright © 2019 Jeannot. All rights reserved.
//

#include "ukkonen.h"

typedef int32_t indNode ;

typedef int32_t indNode ;
// suffix tree node
typedef struct STNode STNode ;
struct STNode {
    // [start,end[ edge label for edge connection with parent
    int32_t start;
    int32_t end ;
    //pointer to other node via suffix link
    indNode suffixLink;
    indNode children[1]; // pseudo size
} ;

#define unknownLeaf  (-1)
#define pseudoLeaf   (-2)

#define lastNewNULL -1
#define indNodeNULL  0

#define bitPseudoLeaf   0x80000000
#define MaskSuffixLink   0x7fffffff


#define OFFSTART        0
#define OFFEND          1
#define OFFSUFFLINK     2
#define OFFCHILD        3

#define NODE(indN)              ((STNode *) (ST->tbNode+(indN)))
#define START(indN)             NODE(indN)->start
#define END(indN)               NODE(indN)->end
#define SUFFLINK(indN)          NODE(indN)->suffixLink
// next macros only used in last pass where the field suffixLink is also use to mark pseudoLeaf
#define GETSUFFLINK(indN)       (NODE(indN)->suffixLink & MaskSuffixLink)
#define SETSUFFLINK(indN,newSuffLink)      NODE(indN)->suffixLink = (NODE(indN)->suffixLink & bitPseudoLeaf) | (newSuffLink)
#define SETPSDLEAF(indN)        NODE(indN)->suffixLink |= bitPseudoLeaf
#define CLEARPSDLEAF(indN)        NODE(indN)->suffixLink &= MaskSuffixLink

#define ISPSDLEAF(indN)         (NODE(indN)->suffixLink & bitPseudoLeaf)
#define CHILD(indN)             NODE(indN)->children

// suffix tree
struct STree {
    const uint8_t *suit ;
    indNode activeNode ;
    // activeEdge by its input string character
    int32_t activeEdge ;
    int32_t activeLength;
    // remainingSuffixCount tells how many suffixes yet to be added in tree
    int32_t remainingSuffixCount ;
    int32_t leafEnd ;
    int32_t nbNodeAlloc ;
    int32_t nxtFreeNode ;
    int32_t *tbNode  ;
    int32_t sizeNode ;
    int32_t maxChild ;
    int32_t minLCP ;
    indNode root ;
    int32_t *sa ;
    int32_t *lcp ;
    int32_t indSA ;
    int32_t minLabelHeight ;
}  ;

static indNode newNode(STree *ST,int start, int end) {
    indNode indN = ST->sizeNode * ST->nxtFreeNode++ ;
    STNode * newN =  (STNode *) (ST->tbNode + indN) ;
    memset(newN->children,0,ST->sizeNode-OFFCHILD);
    // For root node, suffixLink will be set to NULL
    // For internal nodes, suffixLink set to root by default
    // and can be may changed in next extension
    newN->suffixLink = ST->root;
    newN->start = start;
    newN->end = end;
    ST->nbNodeAlloc++ ;
    return indN;
}
static int edgeLength(STree *ST,indNode n) { return (END(n) < 0) ?(ST->leafEnd - START(n) ) : (END(n) - START(n)) ; }

static int walkDown(STree *ST,indNode currNode) {
    // (Trick 1) activePoint change for walk down (APCFWD) using
    //  Skip/Count Trick  . If activeLength is greater
    // than current edge length, set next  internal node as
    //  activeNode and adjust activeEdge and activeLength
    int len ;
    if (ST->activeLength >= (len=edgeLength(ST,currNode)) ) {
        ST->activeEdge += len;
        ST->activeLength -= len;
        ST->activeNode = currNode;
        return 1;
    }
    return 0;
}

static void extendSuffixTree(STree *ST,int pos) {
    // Rule 1, extend all leaves
    ST->leafEnd = pos+1;
    // Increment remainingSuffixCount to treat
    ST->remainingSuffixCount++;
    //set lastNewNode to NULL while starting a new phase,
    indNode lastNewNode = lastNewNULL ;
    //Add all suffixes (yet to be added) one by one in tree
    while(ST->remainingSuffixCount > 0) {
        if (ST->activeLength == 0) ST->activeEdge = pos; //APCFALZ
        // no outgoing edge starting with activeEdge from activeNode
        if (CHILD(ST->activeNode)[ST->suit[ST->activeEdge]] == indNodeNULL) {
            //Extension Rule 2 (A new leaf edge gets created)
            CHILD(ST->activeNode)[ST->suit[ST->activeEdge]] = newNode(ST,pos, -1);
            // new leaf edge is created in above line starting
            // no more node waiting for suffix link to reset
            if (lastNewNode != lastNewNULL) {
                SUFFLINK(lastNewNode) = ST->activeNode;
                lastNewNode = lastNewNULL ;
            }
        } else {
            // There is an outgoing edge starting with activeEdge from activeNode
            // Get the next node at the end of edge starting with activeEdge
            indNode next = CHILD(ST->activeNode)[ST->suit[ST->activeEdge]];
            if (walkDown(ST,next)) { //Do walkdown
                //Start from next node (the new activeNode)
                continue;
            }
            /*Extension Rule 3 (current character being processed
             is already on the edge)*/
            if (ST->suit[START(next) + ST->activeLength] == ST->suit[pos]) {
                //If a newly created node waiting for it's
                //suffix link to be set, then set suffix link
                //of that waiting node to current active node
                if(lastNewNode != lastNewNULL && ST->activeNode != ST->root) {
                    SUFFLINK(lastNewNode) = ST->activeNode;
                    lastNewNode = lastNewNULL ;
                }
                //APCFER3
                ST->activeLength++;
                // STOP all further processing in this phase and move on to next phase*/
                break;
            }
            //activePoint is in middle of the edge being  In this case, we add a new internal node
            // Extension Rule 2, new leaf edge and a newinternal node
            indNode split = newNode(ST,START(next), START(next) + ST->activeLength );
            CHILD(ST->activeNode)[ST->suit[ST->activeEdge]] = split;
            //New leaf coming out of new internal node
            CHILD(split)[ST->suit[pos]] = newNode(ST,pos, -1);
            START(next) += ST->activeLength;
            CHILD(split)[ST->suit[START(next)]] = next;
            
            // lasNewNode waiting ?
            if (lastNewNode != lastNewNULL) {
                // suffixLink of lastNewNode points to current newly created internal node
                SUFFLINK(lastNewNode) = split;
            }
            // Make the  newly cre internal node waiting for it's suffix link reset
            lastNewNode = split;
        }
        // One suffix got added in tree
        ST->remainingSuffixCount--;
        if (ST->activeNode == ST->root && ST->activeLength > 0) {//APCFER2C1
            ST->activeLength--;
            ST->activeEdge = pos - ST->remainingSuffixCount + 1;
        } else if (ST->activeNode != ST->root)  {//APCFER2C2
            ST->activeNode = SUFFLINK(ST->activeNode) ;
        }
    }
}

static void extendSuffixTreeLast(STree *ST) {
    // dont increment leafEnd as virtual add for "$"
    int pos  = ST->leafEnd ;
    ST->remainingSuffixCount++;
    indNode lastNewNode = lastNewNULL;
    while(ST->remainingSuffixCount > 0) {
        if (ST->activeLength == 0)
            ST->activeEdge = pos; //APCFALZ
        if (ST->activeEdge >= ST->leafEnd) { // vrtual "$" is current
            if(ST->activeNode != ST->root){ // mark activeNode as pseudoLeaf
                SETPSDLEAF(ST->activeNode);
            }
        } else if (CHILD(ST->activeNode)[ST->suit[ST->activeEdge]] != indNodeNULL) {
            // There is an outgoing edge starting with activeEdge
            indNode next = CHILD(ST->activeNode)[ST->suit[ST->activeEdge]];
            if (walkDown(ST,next)) { //Do walkdown
                continue;
            }
            //New internal node
            indNode split = newNode(ST,START(next), START(next) + ST->activeLength );
            CHILD(ST->activeNode)[ST->suit[ST->activeEdge]] = split;
            //New leaf coming out of new internal node
            SETPSDLEAF(split);
            START(next) += ST->activeLength;
            CHILD(split)[ST->suit[START(next)]] = next;
            if (lastNewNode != lastNewNULL) {
                SETSUFFLINK(lastNewNode,split) ;
            }
            lastNewNode = split;
        }
        ST->remainingSuffixCount--;
        if (ST->activeNode == ST->root && ST->activeLength > 0)  {//APCFER2C1
            ST->activeLength--;
            ST->activeEdge = pos - ST->remainingSuffixCount + 1;
        } else if (ST->activeNode != ST->root) { //APCFER2C2
            ST->activeNode = GETSUFFLINK(ST->activeNode) ;
        }
    }
}

static void GetSALCP_DP(STree *ST,indNode n, int labelHeight){
    if (n == indNodeNULL) return  ;
    int isLeaf = 1;
    int i;
    //    if(LEAFSTAT(n) != pseudoLeaf) {
    if(!ISPSDLEAF(n)) {
        for (i = 0; i < ST->maxChild; i++) {
            if (CHILD(n)[i] != indNodeNULL ) { isLeaf=0; break ; }
        }
    }
    if(isLeaf  && labelHeight >= ST->minLabelHeight ) {
        if(ST->lcp) { ST->lcp[ST->indSA]  = ST->minLCP ;  }
        if(ST->sa) {  ST->sa[ST->indSA]  = ST->leafEnd - labelHeight ; }
        ST->indSA++ ;
        ST->minLCP = labelHeight ;
    }
    for (i = 0; i < ST->maxChild; i++) {
        if (CHILD(n)[i] != indNodeNULL ) { //Current node not leaf
            if(ST->minLCP > labelHeight) ST->minLCP = labelHeight ;
            GetSALCP_DP(ST,CHILD(n)[i],labelHeight + edgeLength(ST,CHILD(n)[i])) ;
        }
    }
    return  ;
}


void CompSALCP(STree *ST,int32_t *sa,int32_t *lcp) {
    ST->sa = sa ;
    ST->lcp = lcp ;
    ST->indSA = 0 ;
    if(sa != NULL || lcp != NULL) {
        int labelHeight = 0;
        ST->minLCP = 0 ;
        GetSALCP_DP(ST,ST->root, labelHeight);
    }
    return ;
}

// build suffixTree
// character from [0..length[ must be in ]0,maxChar[
// suit must be terminated by 0 (suit[length]=0 )
// for a suit of 1,2 terminated by a 0, maxChar=3
STree * X_CompSuffixTree(const uint8_t *suit,int maxChar,int length) {
    STree * ST = calloc(1,sizeof(ST[0])) ;
    ST->activeNode = indNodeNULL;
    /*activeEdge is represeted as input string character
     index (not the character itself)*/
    ST->activeEdge = -1;
    ST->activeLength = 0;
    ST->suit = suit ;
    // remainingSuffixCount tells how many suffixes yet to
    // be added in tree
    ST->remainingSuffixCount = 0;
    ST->leafEnd = -1;
    ST->maxChild =maxChar+1 ;
    ST->sizeNode = OFFCHILD+ST->maxChild ;
    int i;
    ST->tbNode = malloc(2*length*ST->sizeNode*sizeof(ST->tbNode[0])) ;
    ST->nbNodeAlloc = 2*length ;
    ST->nxtFreeNode = 1 ;
    /*Root is a special node with start and end indices as -1,
     as it has no parent from where an edge comes to root*/
    ST->root = newNode(ST,-1, -1);
    
    ST->activeNode = ST->root; //First activeNode will be root
    for (i=0; i<length+1; i++) extendSuffixTree(ST,i);
    ST->minLabelHeight = 2 ;
    //    int labelHeight = 0;
    //    setSuffixIndexByDFS(ST,ST->root, labelHeight);
    return ST ;
}

// build suffixTree
// character from [0..length[ must be in [0,maxChar[
// For a suit of 0,1 maxChar=2
STree *  ComputeSuffixTree(const uint8_t *suit,int maxChar,int length) {
    STree * ST = calloc(1,sizeof(ST[0])) ;
    
    ST->activeNode = indNodeNULL;
    /*activeEdge is represeted as input string character
     index (not the character itself)*/
    ST->activeEdge = -1;
    ST->activeLength = 0;
    
    // remainingSuffixCount tells how many suffixes yet to
    // be added in tree
    ST->remainingSuffixCount = 0;
    ST->leafEnd = -1;
    
    ST->sa = NULL ;
    ST->lcp = NULL ;
    
    ST->suit = suit ;
    ST->maxChild =maxChar ;
    ST->sizeNode = OFFCHILD+ST->maxChild ;
    int i;
    ST->tbNode = malloc(2*length*ST->sizeNode*sizeof(ST->tbNode[0])) ;
    ST->nbNodeAlloc = 2*(length+1) ;
    ST->nxtFreeNode = 1 ;
    /*Root is a special node with start and end indices as -1,
     as it has no parent from where an edge comes to root*/
    ST->root = newNode(ST,-1, -1);
    ST->activeNode = ST->root; //First activeNode will be root
    for (i=0; i<length; i++) extendSuffixTree(ST,i);
    extendSuffixTreeLast(ST);
    ST->minLabelHeight = 1 ;
    return ST ;
}
STree *  FreeSuffixTree(STree * ST) {
    if(ST) {
        free(ST->tbNode) ;
        free(ST) ;
    }
    return NULL ;
}

int32_t NbNodeUsed(STree * ST) {
    return ST->nxtFreeNode ;
}


