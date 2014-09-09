/* -----  RSTJoin.h  ----- */
#ifndef __RSTJoin_h
#define __RSTJoin_h


#include "RStarTree.h"


/* declarations */

void JnCntNv(RSTREE R1, RSTREE R2,
             int depth,
             typrect R1rectangle1,
             typrect R1rectangle2,
             typrect R2rectangle1,
             typrect R2rectangle2,
             DirQueryProc Dir1Query,
             DataQueryProc Data1Query,
             DirQueryProc Dir2Query,
             DataQueryProc Data2Query,
             DirQueryProc DirJoin,
             DataQueryProc DataJoin,
             int *keysqualifying,
             int *mark);


void JnNv(RSTREE R1, RSTREE R2,
          int depth,
          typrect R1rectangle1,
          typrect R1rectangle2,
          typrect R2rectangle1,
          typrect R2rectangle2,
          DirQueryProc Dir1Query,
          DataQueryProc Data1Query,
          DirQueryProc Dir2Query,
          DataQueryProc Data2Query,
          DirQueryProc DirJoin,
          DataQueryProc DataJoin,
          JoinManageProc Manage,
          void *buf1,
          void *buf2,
          boolean *finish);


#endif /* !__RSTJoin_h */

