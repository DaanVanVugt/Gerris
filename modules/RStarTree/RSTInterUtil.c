/* ----- RSTInterUtil.c ----- */


#include "RStarTree.h"
#include "RSTInterUtil.h"


/* ----- declarations ----- */

/***********************************************************************/

void CreateRSFiles(RSTREE R)

{
  RSTName SufName;
  int O_MODE= O_RDWR | O_CREAT | O_EXCL;
  /* int O_MODE= O_RDWR | O_CREAT | O_TRUNC;     "forced Creation" */
    
  (*R).dir.f= open((*R).dirname,O_MODE,STDMODE);
  if ((*R).dir.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,datasuffix);
  (*R).data.f= open(SufName,O_MODE,STDMODE);
  if ((*R).data.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,dirPDsuffix);
  (*R).dirPD.f= open(SufName,O_MODE,STDMODE);
  if ((*R).dirPD.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,dataPDsuffix);
  (*R).dataPD.f= open(SufName,O_MODE,STDMODE);
  if ((*R).dataPD.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
}

/***********************************************************************/

void OpenRSFiles(RSTREE R, int O_MODE)

{
  RSTName SufName;
    
  (*R).dir.f= open((*R).dirname,O_MODE,STDMODE);
  if ((*R).dir.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,datasuffix);
  (*R).data.f= open(SufName,O_MODE,STDMODE);
  if ((*R).data.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,dirPDsuffix);
  (*R).dirPD.f= open(SufName,O_MODE,STDMODE);
  if ((*R).dirPD.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
  strcpy(SufName,(*R).dirname);
  strcat(SufName,dataPDsuffix);
  (*R).dataPD.f= open(SufName,O_MODE,STDMODE);
  if ((*R).dataPD.f == -1) {
    (*R).RSTDone= FALSE;
    return;
  }
}

/***********************************************************************/

void FastCloseRSFiles(RSTREE R)

{
  close((*R).dir.f);
  close((*R).data.f);
  close((*R).dirPD.f);
  close((*R).dataPD.f);
}

/***********************************************************************/

void CloseRSFiles(RSTREE R)

{
  if (close((*R).dir.f) == -1) {
    (*R).RSTDone= FALSE;
  }
  if (close((*R).data.f) == -1) {
    (*R).RSTDone= FALSE;
  }
  if (close((*R).dirPD.f) == -1) {
    (*R).RSTDone= FALSE;
  }
  if (close((*R).dataPD.f) == -1) {
    (*R).RSTDone= FALSE;
  }
}

/***********************************************************************/

void SetBase(RSTREE R, int pagelen, boolean unique)

{
  refparameters par;
  refpagedir dpd;
  int dirminfill, dataminfill;
  int i;
  
  par= &(*R).parameters._;
  
  (*par).minfillpercent= minimumFillPercentage;
  (*par).reinstpercent= ReInsertPercentage;
  (*par).nbrsexam= maximumNeighborsToExamine;
  (*par).maxdim= NumbOfDim-1;
  (*par).SIZEinfo= sizeof(typinfo);
  (*par).unique= unique;
  (*par).pagelen= pagelen;
  
  /* ----- set directory level parameters ----- */
  SetCheckDir(R,TRUE);
  (*par).dirM= ((*par).pagelen-(*par).SIZE_DIRnofentries) / (*par).direntrylen;
  if ((*par).dirM > maxentries) {
    (*R).RSTDone= FALSE;
    return;
  }
  if ((*par).dirM < 3) {
    (*R).RSTDone= FALSE;
    return;
  }
  if ((*par).dirM < 5) {
    dirminfill= minFillPercSmallFanout;
  }
  else {
    dirminfill= minimumFillPercentage;
  }
  (*par).dirMwords= (*par).dirM / sizeof(int) + 1;
  (*par).dirm= (dirminfill * (*par).dirM + 50) / 100;
  (*par).dirreinsertqty= ((*par).reinstpercent * (*par).dirM + 50) / 100;
  
  /* ----- set data level parameters ----- */
  SetCheckData(R,TRUE);
  (*par).dataM= ((*par).pagelen-(*par).SIZE_DATAnofentries) /
                                                           (*par).dataentrylen;
  if ((*par).dataM > maxentries) {
    (*R).RSTDone= FALSE;
    return;
  }
  if ((*par).dataM < 1) {
    (*R).RSTDone= FALSE;
    return;
  }
  if ((*par).dataM < 5) {
    dataminfill= minFillPercSmallFanout;
  }
  else {
    dataminfill= minimumFillPercentage;
  }
  (*par).dataMwords= (*par).dataM / sizeof(int) + 1;
  (*par).datam= (dataminfill * (*par).dataM + 50) / 100;
  (*par).datareinsertqty=
                          ((*par).reinstpercent * (*par).dataM + 50) / 100;
  
  /* ----- set defaults ----- */
  (*par).height= 1;
  (*par).dirpagecount= 1;
  (*par).datapagecount= 0;
  for (i= 0; i < chainlen; i++) {
    (*par).pagecountarr[i]= 0;
  }
  (*par).recordcount= 0;
  
  dpd= &(*R).dirpagedir._;
  
  (*dpd).childnr= firstPDblocknumb;
  (*dpd).nofnumbers= 0;
  (*dpd).number[0]= rootblocknumb;
  
  dpd= &(*R).datapagedir._;
  
  (*dpd).childnr= firstPDblocknumb;
  (*dpd).nofnumbers= 0;
  (*dpd).number[0]= rootblocknumb;
}

/***********************************************************************/

void SetCheckDir(RSTREE R, boolean creation)

{
          /* Alignment and compatibility checks */
             
  refparameters par;
  int SIZE_DIRnodeOf3, SIZE_DIRnodeOf2,
      PACKEDdirentrylen, REALdirentrylen,
      PACKEDSIZE_DIRnofentries, REALSIZE_DIRnofentries;
  
  par= &(*R).parameters._;
  
  if (creation) {
    (*par).direntrylen= sizeof(typDIRent);			/* set */
    PACKEDdirentrylen= sizeof(typrect) + sizeof(int);
#if 0
    if (PACKEDdirentrylen != (*par).direntrylen) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Directory entries are not packed!");
      fprintf(stderr,"%s %d\n","Packed directory entry length:",PACKEDdirentrylen);
      fprintf(stderr,"%s %d\n","            Real space needed:",(*par).direntrylen);
      fprintf(stderr,"%s\n",/* implicitly */"Applying the latter!");
    }
#endif
  }
  
  SIZE_DIRnodeOf3= sizeof(typDIRnodeOf3);
  SIZE_DIRnodeOf2= sizeof(typDIRnodeOf2);
  REALdirentrylen= SIZE_DIRnodeOf3 - SIZE_DIRnodeOf2;
  if (creation) {
    if (REALdirentrylen != (*par).direntrylen) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Directory nodes are not packed!");
      fprintf(stderr,"%s %d\n","Directory entry length:",(*par).direntrylen);
      fprintf(stderr,"%s %d\n","Space needed in a node:",REALdirentrylen);
      fprintf(stderr,"%s\n",/* explicitly */"Applying the latter!");
      (*par).direntrylen= REALdirentrylen;			/* reset */
    }
  }
  else {
    if ((*par).direntrylen != REALdirentrylen) {		/* check */
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"Incompatible R*-tree file '%s'\n", R->dirname);
      fprintf(stderr,"%s %d\n","Size of a directory entry:",(*par).direntrylen);
      fprintf(stderr,"%s %d\n","                Expecting:",REALdirentrylen);
    }
  }
  
  PACKEDSIZE_DIRnofentries= sizeof(int);
  REALSIZE_DIRnofentries= SIZE_DIRnodeOf3 - 3 * (*par).direntrylen;
  (*par).SIZE_DIRnofentries= REALSIZE_DIRnofentries;		/* set */
  if (creation) {
#if 0
    if (PACKEDSIZE_DIRnofentries != (*par).SIZE_DIRnofentries) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Gap before directory entries!");
    }
#endif
  }
  else {
    if ((*par).SIZE_DIRnofentries != REALSIZE_DIRnofentries) {
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"Incompatible R*-tree file '%s'\n", R->dirname);
      fprintf(stderr,"%s %d\n","Offset for entries:",(*par).SIZE_DIRnofentries);
      fprintf(stderr,"%s %d\n","         Expecting:",REALSIZE_DIRnofentries);
    }
  }
  
  if (! creation) {
    if ((*par).maxdim+1 != NumbOfDim) {
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"Incompatible R*-tree file '%s'\n", R->dirname);
      fprintf(stderr,"%s %d\n","Number of dimensions:",(*par).maxdim+1);
      fprintf(stderr,"%s %d\n","           Expecting:",NumbOfDim);
    }
  }
  
  (*R).dirnodelen= (*par).pagelen + (*par).direntrylen;
}

/***********************************************************************/

void SetCheckData(RSTREE R, boolean creation)

{
          /* Alignment and compatibility checks */
             
  refparameters par;
  int SIZE_DATAnodeOf3, SIZE_DATAnodeOf2,
      PACKEDdataentrylen, REALdataentrylen,
      PACKEDSIZE_DATAnofentries, REALSIZE_DATAnofentries;
  
  par= &(*R).parameters._;
  
  if (creation) {
    (*par).dataentrylen= sizeof(typDATAent);			/* set */
    PACKEDdataentrylen= sizeof(typrect) + sizeof(typinfo);
    if (PACKEDdataentrylen != (*par).dataentrylen) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Data entries are not packed!");
      fprintf(stderr,"%s %d\n","Packed data entry length:",PACKEDdataentrylen);
      fprintf(stderr,"%s %d\n","       Real space needed:",(*par).dataentrylen);
      fprintf(stderr,"%s\n",/* implicitly */"Applying the latter!");
    }
  }
  
  SIZE_DATAnodeOf3= sizeof(typDATAnodeOf3);
  SIZE_DATAnodeOf2= sizeof(typDATAnodeOf2);
  REALdataentrylen= SIZE_DATAnodeOf3 - SIZE_DATAnodeOf2;
  if (creation) {
    if (REALdataentrylen != (*par).dataentrylen) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Data nodes are not packed!");
      fprintf(stderr,"%s %d\n","     Data entry length:",(*par).dataentrylen);
      fprintf(stderr,"%s %d\n","Space needed in a node:",REALdataentrylen);
      fprintf(stderr,"%s\n",/* explicitly */"Applying the latter!");
      (*par).dataentrylen= REALdataentrylen;			/* reset */
    }
  }
  else {
    if ((*par).dataentrylen != REALdataentrylen) {		/* check */
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"Incompatible R*-tree file '%s'\n", R->dirname);
      fprintf(stderr,"%s %d\n","Size of a data entry:",(*par).dataentrylen);
      fprintf(stderr,"%s %d\n","           Expecting:",REALdataentrylen);
    }
  }
  
  PACKEDSIZE_DATAnofentries= sizeof(int);
  REALSIZE_DATAnofentries= SIZE_DATAnodeOf3 - 3 * (*par).dataentrylen;
  (*par).SIZE_DATAnofentries= REALSIZE_DATAnofentries;		/* set */
  if (creation) {
    if (PACKEDSIZE_DATAnofentries != (*par).SIZE_DATAnofentries) {
      fprintf(stderr,"\n%s\n","     -----  WARNING  -----");
      fprintf(stderr,"%s\n","Gap before data entries!");
    }
  }
  else {
    if ((*par).SIZE_DATAnofentries != REALSIZE_DATAnofentries) {
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"Incompatible R*-tree file '%s'\n", R->dirname);
      fprintf(stderr,"%s %d\n","Offset for entries:",(*par).SIZE_DATAnofentries);
      fprintf(stderr,"%s %d\n","         Expecting:",REALSIZE_DATAnofentries);
    }
  }
  
  if (! creation) {
    if ((*par).SIZEinfo != sizeof(typinfo)) {
      fprintf(stderr,"\n%s\n","FATAL ERROR:");
      fprintf(stderr,"%s %d\n","Size of an info part:",(*par).SIZEinfo);
      fprintf(stderr,"%s %d\n","           Expecting:",sizeof(typinfo));
    }
  }
  
  (*R).datanodelen= (*par).pagelen + (*par).dataentrylen;
}

/***********************************************************************/

void InitChainFlags(RSTREE R)

{
  int i;
  
  for (i= 1; i <= chainlen; i++) {
    (*R).N[i]= NULL; (*R).NInst[i]= NULL; (*R).NDel[i]= NULL;
    (*R).E[i]= -1; (*R).EInst[i]= -1;
    (*R).P[i]= 0;
    (*R).Nmodified[i]= FALSE;
    (*R).ReInsert[i]= FALSE;
  }
}

/***********************************************************************/

void AllocBuffers(RSTREE R)

{
  refparameters par;
  int i;
  
  par= &(*R).parameters._;
  
  for (i= 1; i < (*par).height; i++) {
    (*R).N[i]= (refnode)malloc((*R).dirnodelen);
  }
  (*R).N[(*par).height]= (refnode)malloc((*R).datanodelen);
  if ((*R).dirnodelen > (*R).datanodelen) {
    (*R).Nsibling= (refnode)malloc((*R).dirnodelen);
  }
  else {
    (*R).Nsibling= (refnode)malloc((*R).datanodelen);
  }
  (*R).helpdirnode= (refnode)malloc((*R).dirnodelen);
  (*R).helpdatanode= (refnode)malloc((*R).datanodelen);
}

/***********************************************************************/

void DeallocBuffers(RSTREE R)

{
  refparameters par;
  int i;
  
  par= &(*R).parameters._;
  
  for (i= 1; i < (*par).height; i++) {
    free((*R).N[i]); (*R).N[i]= NULL;
  }
  free((*R).N[(*par).height]); (*R).N[(*par).height]= NULL;
  free((*R).Nsibling);
  free((*R).helpdirnode);
  free((*R).helpdatanode);
}

/***********************************************************************/

void InitCount(RSTREE R)

{
  refcount c;
  
  c= &(*R).count;
  (*c).countflag= FALSE;
  (*c).dirvisitcount= 0; (*c).datavisitcount= 0;
  (*c).dirreadcount= 0; (*c).datareadcount= 0;
  (*c).dirmodifycount= 0; (*c).datamodifycount= 0;
  (*c).dirwritecount= 0; (*c).datawritecount= 0;
}

/***********************************************************************/
