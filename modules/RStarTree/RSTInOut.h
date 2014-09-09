/* -----  RSTInOut.h  ----- */
#ifndef __RSTInOut_h
#define __RSTInOut_h


#include "RStarTree.h"


/* declarations */

void PutNode(RSTREE R,
             refnode nodeptr,
             int pagenr,
             int depth);

void GetNode(RSTREE R,
             refnode nodeptr,
             int pagenr,
             int depth);

void NewNode(RSTREE R,
             int depth);

void ReadPage(RSTREE R,
              typfiledesc fd,
              int pagenr,
              void *block);

void WritePage(RSTREE R,
               typfiledesc fd,
               int pagenr,
               void *block);

void GetPageNr(RSTREE R,
               int *pagenr,
               int depth);

void PutPageNr(RSTREE R,
               int pagenr,
               int depth);


#endif /* !__RSTInOut_h */
