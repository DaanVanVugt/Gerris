/* -----  RSTInterUtil.h  ----- */
#ifndef __RSTInterUtil_h
#define __RSTInterUtil_h


#include "RStarTree.h"


/* constants */

#define minimumFillPercentage     40
#define minFillPercSmallFanout    50
#define ReInsertPercentage        30
#define maximumNeighborsToExamine 32


/* declarations */

void CreateRSFiles(RSTREE R);
void OpenRSFiles(RSTREE R, int O_MODE);
void FastCloseRSFiles(RSTREE R);
void CloseRSFiles(RSTREE R);
void SetBase(RSTREE R, int pagelen, boolean unique);
void SetCheckDir(RSTREE R, boolean creation);
void SetCheckData(RSTREE R, boolean creation);
void InitChainFlags(RSTREE R);
void AllocBuffers(RSTREE R);
void DeallocBuffers(RSTREE R);
void InitCount(RSTREE R);


#endif /* !__RSTInterUtil_h */
