/* -----  RSTUtil.h  ----- */
#ifndef __RSTUtil_h
#define __RSTUtil_h


#include "RStarTree.h"


/* declarations */

void     FalseArray(int *ptr,
                    int wordsqty);

void     CopyRect(RSTREE,
                  typrect from,
                  typrect to);

void     EvalCenter(RSTREE R,
                    typrect rectangle,
                    typpoint center);

double   RSTDistance(RSTREE R,
                     typpoint point1,
                     typpoint point2);

void     EvalDirEnclRect(RSTREE R,
                         refDIRnode node,
                         typrect rect);

void     EvalDataEnclRect(RSTREE R,
                          refDATAnode node,
                          typrect rect);

boolean  Overlaps(RSTREE R,
                  typrect rect1,
                  typrect rect2);

boolean  RSTEqual(RSTREE R,
                  typrect rect1,
                  typrect rect2);

boolean  Covers(RSTREE R,
                typrect crect,
                typrect rect);

void     GetOverlap(RSTREE R,
                    typrect r1,
                    typrect r2,
                    double *spc);

void     QuickSortValArr(int begin,
                         int end,
                         ValueArray value,
                         IndexArray I);

void     QuickSortDirEnt(int begin,
                         int end,
                         int dim,
                         Side side,
                         typDIRentries Ntosort,
                         IndexArray I);
                     
void     QuickSortDataEnt(int begin,
                          int end,
                          int dim,
                          Side side,
                          typDATAentries Ntosort,
                          IndexArray I);

void     ChooseSubtree(RSTREE R,
                       typrect newrect,
                       int depth,
                       refDIRnode node,
                       int *found);

void     AdjustChain(RSTREE R,
                     int depth,
                     typrect newrect);

void     AdjustChainAfterDeletion(RSTREE R,
                                  int depth);


#endif /* !__RSTUtil_h */
