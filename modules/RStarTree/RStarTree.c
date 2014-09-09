/* ----- RStarTree.c ----- */

#include <stdlib.h>

#include "RStarTree.h"
#include "RSTInterUtil.h"
#include "RSTInOut.h"
#include "RSTInstDel.h"
#include "RSTQuery.h"
#include "RSTJoin.h"


/* constants */


/* types */


/* declarations */

static void BasicCheck(void);
/*** --- Begin --- unused ***
static boolean InternEqual(RSTREE rst,
                           typrect RSTrect,
                           typrect qrect1,
                           typrect qrect2);
static boolean InternEncloses(RSTREE rst,
                              typrect RSTrect,
                              typrect qrect1,
                              typrect qrect2);
 *** ---   End --- unused ***/

/************************************************************************/

void NoRSTree(RSTREE *r)

{
  *r= NULL;
}

/************************************************************************/

boolean CreateRST(const 
char *name,
                  int pagelen,
                  boolean unique)

{
  RSTREE R;
  refparameters par;
  
  BasicCheck();
  
  R= (RSTREE)malloc(sizeof(rstree));
  (*R).readonly = FALSE;
  strcpy((*R).dirname,name);
  (*R).RSTDone= TRUE;
  CreateRSFiles(R);
  if (! (*R).RSTDone) {
    free(R); R= NULL;
    return FALSE;
  }
  SetBase(R,pagelen,unique);
  if (! (*R).RSTDone) {
    free(R); R= NULL;
    return FALSE;
  }
  (*R).dirPD.bl= SIZEfixblock;
  
  par= &(*R).parameters._;
  
  WritePage(R,(*R).dirPD,paramblocknumb,par);
  WritePage(R,(*R).dirPD,firstPDblocknumb,&(*R).dirpagedir);
  (*R).dataPD.bl= SIZEfixblock;
  WritePage(R,(*R).dataPD,paramblocknumb,par);     /* -- unused -- */
  WritePage(R,(*R).dataPD,firstPDblocknumb,&(*R).datapagedir);
  (*R).data.bl= (*par).pagelen;
  (*R).N[1]= (refnode)malloc((*R).datanodelen);
  (*(*R).N[1]).DATA.nofentries= 0;
  PutNode(R,(*R).N[1],rootblocknumb,1);
  free((*R).N[1]);
  if (! (*R).RSTDone) {
    free(R); R= NULL;
    return FALSE;
  }
  CloseRSFiles(R);
  if (! (*R).RSTDone) {
    free(R); R= NULL;
    return FALSE;
  }
  free(R); R= NULL;
  return TRUE;
}

/************************************************************************/

boolean RemoveRST(const char *name)

{
  RSTName SufName;
  boolean success= TRUE;

  if (unlink(name) != 0) {
    success= FALSE;
  }
  strcpy(SufName,name);
  strcat(SufName,datasuffix);
  if (unlink(SufName) != 0) {
    success= FALSE;
  }
  strcpy(SufName,name);
  strcat(SufName,dirPDsuffix);
  if (unlink(SufName) != 0) {
    success= FALSE;
  }
  strcpy(SufName,name);
  strcat(SufName,dataPDsuffix);
  if (unlink(SufName) != 0) {
    success= FALSE;
  }
  return success;
}

/************************************************************************/

boolean OpenRST(RSTREE *r,
                const char *name,
		const char *mode)

{
  RSTREE R;
  refparameters par;
  
  if (*r != NULL) {
    return FALSE;
  }
  *r= (RSTREE)malloc(sizeof(rstree));
  R= *r;
  strcpy((*R).dirname,name);
  (*R).RSTDone= TRUE;
  (*R).readonly = strcmp (mode, "rw");
  OpenRSFiles(R, (*R).readonly ? O_RDONLY : O_RDWR);
  if (! (*R).RSTDone) {
    free(R); *r= NULL;
    return FALSE;
  }
  InitChainFlags(R);
  (*R).dirPD.bl= SIZEfixblock;
  
  par= &(*R).parameters._;
  
  ReadPage(R,(*R).dirPD,paramblocknumb,par);
  ReadPage(R,(*R).dirPD,firstPDblocknumb,&(*R).dirpagedir);
  (*R).dataPD.bl= SIZEfixblock;
  ReadPage(R,(*R).dataPD,firstPDblocknumb,&(*R).datapagedir);
  if (! (*R).RSTDone) {
    FastCloseRSFiles(R);
    free(R); *r= NULL;
    return FALSE;
  }
  SetCheckDir(R,FALSE);
  SetCheckData(R,FALSE);
  AllocBuffers(R);
  (*R).dir.bl= (*par).pagelen;
  (*R).data.bl= (*par).pagelen;
  GetNode(R,(*R).N[1],rootblocknumb,1); (*R).P[1]= rootblocknumb;
  InitCount(R);
  if (! (*R).RSTDone) {
    FastCloseRSFiles(R);
    DeallocBuffers(R);
    free(R); *r= NULL;
    return FALSE;
  }
  return TRUE;
}

/************************************************************************/

boolean CloseRST(RSTREE *r)

{
  RSTREE R;
  refparameters par;
  boolean success;
  int i;
  
  if (*r == NULL) {
    return FALSE;
  }
  R= *r;
  (*R).RSTDone= TRUE;

  if (!(*R).readonly) {
    par= &(*R).parameters._;
    WritePage(R,(*R).dirPD,paramblocknumb,par);
    WritePage(R,(*R).dirPD,firstPDblocknumb,&(*R).dirpagedir);
    WritePage(R,(*R).dataPD,paramblocknumb,par);     /* -- unused -- */
    WritePage(R,(*R).dataPD,firstPDblocknumb,&(*R).datapagedir);
    for (i= 1; i <= (*par).height; i++) {
      if ((*R).Nmodified[i]) {
	PutNode(R,(*R).N[i],(*R).P[i],i);
      }
    }
  }
  if (! (*R).RSTDone) {
    return (*R).RSTDone;
  }
  CloseRSFiles(R);
  if (! (*R).RSTDone) {
    return (*R).RSTDone;
  }
  DeallocBuffers(R);
  success= (*R).RSTDone;
  free(R); *r= NULL;
  return success;
}

/************************************************************************/

boolean SetUnique(RSTREE R,
                  boolean mode)

{
  if (R == NULL) {
    return FALSE;
  }
  (*R).parameters._.unique= mode;
  return TRUE;
}

/************************************************************************/

boolean InsertRecord(RSTREE R,
                     typrect rectangle,
                     typinfo *info,
                     boolean *inserted)

{
  refparameters par;
  typentry entry;
  refinfo infoadr;
  int d;
  
  if (R == NULL) {
    *inserted= FALSE;
    return FALSE;
  }
  
  (*R).RSTDone= TRUE;
  
  par= &(*R).parameters._;
  
  if ((*par).unique) {
    *inserted= ! FoundRect(R,1,rectangle,TRUE,&infoadr);
  }
  else {
    *inserted= TRUE;
  }
  if (*inserted) {
    for (d= 0; d <= (*par).maxdim; d++) {
      entry.DATA.rect[d]= rectangle[d];      
    }
    entry.DATA.info= *info;
    (*R).ReInsert[(*par).height]= TRUE; /* general switch for Forced ReInsert */
    Insert(R,&entry,(*par).height);
    (*R).ReInsert[(*par).height]= FALSE;
    *inserted= (*R).RSTDone;
    if (*inserted) {
      (*par).recordcount++;
    }
  }
  return (*R).RSTDone;
}

/************************************************************************/

boolean DeleteRecord(RSTREE R,
                     typrect rectangle,
                     boolean *deleted)

{
  refinfo infoadr;
    
  if (R == NULL) {
    *deleted= FALSE;
    return FALSE;
  }
  
  (*R).RSTDone= TRUE;
  
  *deleted= FoundRect(R,1,rectangle,FALSE,&infoadr);
  if (*deleted) {
    DeleteOneRec(R);
    *deleted= (*R).RSTDone;
    if (*deleted) {
      (*R).parameters._.recordcount--;
    }
  }
  return (*R).RSTDone;
}

/************************************************************************/

boolean ExistsRegion(RSTREE R,
                     typrect rectangle1,
                     typrect rectangle2,
                     DirQueryProc DirQuery,
                     DataQueryProc DataQuery,
                     boolean *regionfound)

{
  if (R == NULL) {
    *regionfound= FALSE;
    return FALSE;
  }
  /*
  int i;
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  (*R).RSTDone= TRUE;
  *regionfound= FALSE;
  XstsRgn(R,1,rectangle1,rectangle2,DirQuery,DataQuery,regionfound);
  return (*R).RSTDone;
}

/************************************************************************/

boolean RegionCount(RSTREE R,
                    typrect rectangle1,
                    typrect rectangle2,
                    DirQueryProc DirQuery,
                    DataQueryProc DataQuery,
                    int *recordcount)

{
  if (R == NULL) {
    *recordcount= 0;
    return FALSE;
  }
  /*
  int i;  
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  (*R).RSTDone= TRUE;
  *recordcount= 0;
  RgnCnt(R,1,rectangle1,rectangle2,DirQuery,DataQuery,recordcount);
  return (*R).RSTDone;
}

/************************************************************************/

boolean RegionQuery(RSTREE R,
                    typrect rectangle1,
                    typrect rectangle2,
                    DirQueryProc DirQuery,
                    DataQueryProc DataQuery,
                    QueryManageProc ManageProc,
                    void *buf)

{
  boolean finish;
  
  if (R == NULL) {
    return FALSE;
  }
  /*
  int i;
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  (*R).RSTDone= TRUE;
  finish= FALSE;
  RgnQuery(R,1,rectangle1,rectangle2,DirQuery,DataQuery,
           ManageProc,buf,&finish);
  return (*R).RSTDone;
}

/************************************************************************/

boolean RegionQueryInfo(RSTREE R,
			Check includes,
			Check intersects,
			void * data,
			typrect rect,
			typdirinfo * info)

{
  if (R == NULL) {
    return FALSE;
  }
  /*
  int i;
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  (*R).RSTDone= TRUE;
  RgnQueryInfo(R,1,includes,intersects,data,rect,info);
  return (*R).RSTDone;
}

/************************************************************************/

boolean AllQuery(RSTREE R,
                 QueryManageProc ManageProc,
                 void *buf)

{
  boolean finish;
  
  if (R == NULL) {
    return FALSE;
  }
  /*
  int i;
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  (*R).RSTDone= TRUE;
  finish= FALSE;
  All(R,1,ManageProc,buf,&finish);
  return (*R).RSTDone;
}

/************************************************************************/

boolean Update(RSTREE R)

{
  typdirinfo info;

  if (R == NULL) {
    return FALSE;
  }
  /*
  int i;
  for (i= 2; i <= (*R).parameters._.height; i++) {
    if ((*R).Nmodified[i]) {
      PutNode(R,(*R).N[i],(*R).P[i],i);
      (*R).Nmodified[i]= FALSE;
    }
    (*R).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */

  typrect rect;
  rect[0].l = -0.5; rect[0].h = 0.5;
  rect[1].l = -0.5; rect[1].h = 0.5;
  (*R).RSTDone= TRUE;
  UpdateAll(R,1,rect,&info);
  return (*R).RSTDone;
}

/************************************************************************/

boolean JoinCountNv(RSTREE R1, RSTREE R2,
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
                    int *paircount)

{
  int mark;
  boolean twiceopen, success;
  int i;
  
  if (R1 == NULL || R2 == NULL) {
    *paircount= 0;
    return FALSE;
  }
  twiceopen= R1 == R2;
  if (twiceopen) {
    for (i= 1; i <= (*R1).parameters._.height; i++) {
      if ((*R1).Nmodified[i]) {
        PutNode(R1,(*R1).N[i],(*R1).P[i],i);
        (*R1).Nmodified[i]= FALSE;
      }
    } /* syncronize R1 */
    R2= NULL;
    success= OpenRST(&R2,(*R1).dirname,"rw"); /* NEW(R2) */
    if (! success) {
      fprintf(stderr,"%s\n","FATAL INTERNAL ERROR");
      fprintf(stderr,"%s\n","JoinCountNv 1");
      abort();
    }
  }
  /*
  for (i= 2; i <= (*R1).parameters._.height; i++) {
    if ((*R1).Nmodified[i]) {
      PutNode(R1,(*R1).N[i],(*R1).P[i],i);
      (*R1).Nmodified[i]= FALSE;
    }
    (*R1).P[i]= 0;
  }
  for (i= 2; i <= (*R2).parameters._.height; i++) {
    if ((*R2).Nmodified[i]) {
      PutNode(R2,(*R2).N[i],(*R2).P[i],i);
      (*R2).Nmodified[i]= FALSE;
    }
    (*R2).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  if (! ((*R1).RSTDone && (*R2).RSTDone) ) {
    *paircount= 0;
    return FALSE;
  }
  *paircount= 0;
  mark= 0;
  JnCntNv(R1,R2,
          1,
          R1rectangle1,R1rectangle2,R2rectangle1,R2rectangle2,
          Dir1Query,Data1Query,Dir2Query,Data2Query,DirJoin,DataJoin,
          paircount,
          &mark);
  success= (*R1).RSTDone && (*R2).RSTDone;
  if (twiceopen) {
    success= success && CloseRST(&R2);
  }
  return success;
}

/************************************************************************/

boolean JoinNv(RSTREE R1, RSTREE R2,
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
               void *buf2)

{
  boolean twiceopen, success;
  boolean finish;
  int i;
  
  if (R1 == NULL || R2 == NULL) {
    return FALSE;
  }
  twiceopen= R1 == R2;
  if (twiceopen) {
    for (i= 1; i <= (*R1).parameters._.height; i++) {
      if ((*R1).Nmodified[i]) {
        PutNode(R1,(*R1).N[i],(*R1).P[i],i);
        (*R1).Nmodified[i]= FALSE;
      }
    } /* syncronize R1 */
    R2= NULL;
    success= OpenRST(&R2,(*R1).dirname,"rw"); /* NEW(R2) */
    if (! success) {
      fprintf(stderr,"%s\n","FATAL INTERNAL ERROR");
      fprintf(stderr,"%s\n","JoinCountNv 1");
      abort();
    }
  }
  /*
  for (i= 2; i <= (*R1).parameters._.height; i++) {
    if ((*R1).Nmodified[i]) {
      PutNode(R1,(*R1).N[i],(*R1).P[i],i);
      (*R1).Nmodified[i]= FALSE;
    }
    (*R1).P[i]= 0;
  }
  for (i= 2; i <= (*R2).parameters._.height; i++) {
    if ((*R2).Nmodified[i]) {
      PutNode(R2,(*R2).N[i],(*R2).P[i],i);
      (*R2).Nmodified[i]= FALSE;
    }
    (*R2).P[i]= 0;
  }
  *//* to be inserted, if main memory path shall be initialized
         for test purpose */
  
  if (! ((*R1).RSTDone && (*R2).RSTDone) ) {
    return FALSE;
  }
  finish= FALSE;
  JnNv(R1,R2,
       1,
       R1rectangle1,R1rectangle2,R2rectangle1,R2rectangle2,
       Dir1Query,Data1Query,Dir2Query,Data2Query,DirJoin,DataJoin,
       Manage,buf1,buf2,&finish);
  success= (*R1).RSTDone && (*R2).RSTDone;
  if (twiceopen) {
    success= success && CloseRST(&R2);
  }
  return success;
}

/************************************************************************/

boolean InquireRSTDesc(RSTREE R,
                       char *name,
                       int *numbofdim,
                       int *pagesize,
                       int *sizedirentry,
                       int *sizedataentry,
                       int *sizeinfo,
                       int *maxdirfanout,
                       int *maxdatafanout,
                       int *numbofdirpages,
                       int *numbofdatapages,
                       int pagesperlevel[],
                       int *numbofrecords,
                       int *height,
                       boolean *unique)

{
  refparameters par;
  int i;
  
  if (R == NULL) {
    return FALSE;
  }
  
  strcpy(name,(*R).dirname);
  
  par= &(*R).parameters._;
  
  *numbofdim= (*par).maxdim+1;
  *pagesize= (*par).pagelen;
  *sizedirentry= (*par).direntrylen;
  *sizedataentry= (*par).dataentrylen;
  *sizeinfo= (*par).SIZEinfo;
  *maxdirfanout= (*par).dirM;
  *maxdatafanout= (*par).dataM;
  *numbofdirpages= (*par).dirpagecount;
  *numbofdatapages= (*par).datapagecount;
  pagesperlevel[0]= 1;
  for (i= 1; i < (*par).height; i++) {
    pagesperlevel[i]= (*par).pagecountarr[i+1];
  }
  *numbofrecords= (*par).recordcount;
  *height= (*par).height;
  *unique= (*par).unique;
  return TRUE;
}

/************************************************************************/

boolean CountsOn0(RSTREE R)

{
  refcount c;
  
  if (R == NULL) {
    return FALSE;
  }
  c= &(*R).count;
  (*c).countflag= TRUE;
  (*c).dirvisitcount= 0; (*c).datavisitcount= 0;
  (*c).dirreadcount= 0; (*c).datareadcount= 0;
  (*c).dirmodifycount= 0; (*c).datamodifycount= 0;
  (*c).dirwritecount= 0; (*c).datawritecount= 0;
  return TRUE;
}

/************************************************************************/

boolean CountsOn(RSTREE R)

{
  if (R == NULL) {
    return FALSE;
  }
  (*R).count.countflag= TRUE;
  return TRUE;
}

/************************************************************************/

boolean CountsOff(RSTREE R)

{
  if (R == NULL) {
    return FALSE;
  }
  (*R).count.countflag= FALSE;
  return TRUE;
}

/************************************************************************/

boolean GetCountRead(RSTREE R,
                     int *dirvis, int *datavis,
                     int *dirread, int *dataread)

{
  refcount c;
  
  if (R == NULL) {
    *dirvis= 0;
    *datavis= 0;
    *dirread= 0;
    *dataread= 0;
    return FALSE;
  }
  c= &(*R).count;
  *dirvis= (*c).dirvisitcount;
  *datavis= (*c).datavisitcount;
  *dirread= (*c).dirreadcount;
  *dataread= (*c).datareadcount;
  return TRUE;
}

/************************************************************************/

boolean GetCountWrite(RSTREE R,
                      int *dirmod, int *datamod,
                      int *dirwrite, int *datawrite)

{
  refcount c;
  
  if (R == NULL) {
    *dirmod= 0;
    *datamod= 0;
    *dirwrite= 0;
    *datawrite= 0;
    return FALSE;
  }
  c= &(*R).count;
  *dirmod= (*c).dirmodifycount;
  *datamod= (*c).datamodifycount;
  *dirwrite= (*c).dirwritecount;
  *datawrite= (*c).datawritecount;
  return TRUE;
}

/************************************************************************/

boolean GetMemory(RSTREE R,
                  int *dirpages, int *datapages)

{
  refparameters par;
  
  if (R == NULL) {
    *dirpages= 0;
    *datapages= 0;
    return FALSE;
  }
  par= &(*R).parameters._;
  *dirpages= (*par).dirpagecount;
  *datapages= (*par).datapagecount;
  return TRUE;
}

/************************************************************************/

boolean GetHeight(RSTREE R,
                  int *height)

{
  if (R == NULL) {
    *height= 0;
    return FALSE;
  }
  *height= (*R).parameters._.height;
  return TRUE;
}

/************************************************************************/

static void BasicCheck()
{
  if (sizeof(byte) != 1) {
    fprintf(stderr,"%s\n","FATAL ERROR:");
    fprintf(stderr,"%s\n","BasicCheck 1");
    fprintf(stderr,"%s\n","sizeof(byte) != 1");
    fprintf(stderr,"%s %ld\n","sizeof(byte):",sizeof(byte));
    abort();
    /* concerning application of type ByteArray */
  }
  if (sizeof(int) < 4) {
    fprintf(stderr,"%s\n","BasicCheck 2");
    fprintf(stderr,"%s\n","sizeof(int) < 4");
    fprintf(stderr,"%s %ld\n","sizeof(int):",sizeof(int));
    fprintf(stderr,"%s\n","WARNING: bigger int range assumed.");
  }
  if (sizeof(typinfo) < sizeof(int)) {
    fprintf(stderr,"%s\n","FATAL ERROR:");
    fprintf(stderr,"%s\n","BasicCheck 3");
    fprintf(stderr,"%s\n","sizeof(typinfo) < sizeof(int)");
    fprintf(stderr,"%s %ld\n","sizeof(typinfo):",sizeof(typinfo));
    fprintf(stderr,"%s %ld\n","    sizeof(int):",sizeof(int));
    abort();
  }
  if (sizeof(typpagedir) > sizeof(typfixblock)) {
    fprintf(stderr,"%s\n","FATAL ERROR:");
    fprintf(stderr,"%s\n","BasicCheck 4");
    fprintf(stderr,"%s\n","sizeof(typpagedir) > sizeof(typfixblock)");
    fprintf(stderr,"%s %ld\n"," sizeof(typpagedir):",sizeof(typpagedir));
    fprintf(stderr,"%s %ld\n","sizeof(typfixblock):",sizeof(typfixblock));
    abort();
  }
  if (sizeof(typparameters) > sizeof(typfixblock)) {
    fprintf(stderr,"%s\n","FATAL ERROR:");
    fprintf(stderr,"%s\n","BasicCheck 5");
    fprintf(stderr,"%s\n","sizeof(typparameters) > sizeof(typfixblock)");
    fprintf(stderr,"%s %ld\n","sizeof(typparameters):",sizeof(typparameters));
    fprintf(stderr,"%s %ld\n","  sizeof(typfixblock):",sizeof(typfixblock));
    abort();
  }
}

/************************************************************************/
/*** --- Begin --- unused ***

static boolean InternEqual(RSTREE R,
                           typrect RSTrect,
                           typrect queryrect,
                           typrect unused)

{
  boolean eql;
  int d;
  
  d= -1;
  do {
    d++;
    eql= RSTrect[d].l == queryrect[d].l &&
         RSTrect[d].h == queryrect[d].h;
  } while (eql && d != (*R).parameters._.maxdim);
  return eql;
}

static boolean InternEncloses(RSTREE R,
                              typrect RSTrect,
                              typrect queryrect,
                              typrect unused)

{
  int maxdim;
  boolean encl;
  int d;
  
  maxdim= (*R).parameters._.maxdim;
  d= -1;
  do {
    d++;
    encl= RSTrect[d].l <= queryrect[d].l &&
          RSTrect[d].h >= queryrect[d].h;
  } while (encl && d != maxdim);
  return encl;
}

 *** ---   End --- unused ***/
/***********************************************************************/

boolean Find(RSTREE R,
             typrect rectangle,
             boolean *found,
             void *buf,
             int nbytes)

{
  refinfo infoadr;
  
  if (R == NULL) {
    *found= FALSE;
    return FALSE;
  }
  
  *found= FoundRect(R,1,rectangle,FALSE,&infoadr);
  if (*found) {
    memcpy(buf,infoadr,nbytes);
  }
  return (*R).RSTDone;
}

/************************************************************************/
