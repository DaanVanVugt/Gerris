/* ----- RSTInstDel.c ----- */


#include "RStarTree.h"
#include "RSTInstDel.h"
#include "RSTUtil.h"
#include "RSTInOut.h"


/* declarations */

static void ExcludeDataEntries(RSTREE R,
                               refDATAent newentry,
                               typrect newrect,
                               int depth,
                               int M,
                               int reinsertqty);
static void ExcludeDirEntries(RSTREE R,
                              refDIRent newentry,
                              typrect newrect,
                              int depth,
                              int M,
                              int reinsertqty);
static void Split(RSTREE R,
                  refentry newentry,
                  int *depth,
                  int M,
                  int m);
static void SplitAndDistributDir(RSTREE R,
                                 int depth,
                                 refDIRent newentry,
                                 int M,
                                 int m);
static void SplitAndDistributData(RSTREE R,
                                  int depth,
                                  refDATAent newentry,
                                  int M,
                                  int m);
static void UnDistributData(RSTREE R,
                            refDATAent newentry);
static void GetInstChain(RSTREE R,
                         typrect newrect,
                         int depth);
static void NtoNDel(RSTREE R,
                    int depth);
static void ShrinkTree(RSTREE R);


/***********************************************************************/

static void GetInstChain(RSTREE R,
                         typrect newrect,
                         int depth)

{
  int i;
  refcount c;
  
  i= 1;
  while (i < depth) {
    if ((*R).NInst[i+1] != NULL) {
      /* already in path */
      (*R).E[i]= (*R).EInst[i]; (*R).EInst[i]= -1;
      i++;
      if ((*R).N[i] != (*R).NInst[i]) {
        (*R).P[i]= (*(*R).N[i-1]).DIR.entries[(*R).E[i-1]].ptrtosub;
        free((*R).N[i]); (*R).N[i]= NULL;
        (*R).N[i]= (*R).NInst[i];
      }
      (*R).NInst[i]= NULL;
    }
    else if ((*R).EInst[i] != -1) {
      /* known ... */
      (*R).E[i]= (*R).EInst[i]; (*R).EInst[i]= -1;
      i++;
      if ((*(*R).N[i-1]).DIR.entries[(*R).E[i-1]].ptrtosub != (*R).P[i]) {
        /* but not in path */
        NewNode(R,i);
      }
    }
    else {
      /* not known */
      ChooseSubtree(R,newrect,i,&(*(*R).N[i]).DIR,&(*R).E[i]);
      i++;
      if ((*(*R).N[i-1]).DIR.entries[(*R).E[i-1]].ptrtosub != (*R).P[i]) {
        /* and not in path */
        NewNode(R,i);
      }
    }
  }
  
  c= &(*R).count;
  if ((*c).countflag) {
    if (depth == (*R).parameters._.height) {
      (*c).dirvisitcount+= depth - 1;
      (*c).datavisitcount++;
    }
    else {
      (*c).dirvisitcount+= depth;
    }
  }
}

/***********************************************************************/

void Insert(RSTREE R, typentry *newentry, int depth)

{
  typrect newrect;
  int M, m;
  int reinsertqty, heightbefore;
  refparameters par;
  refDIRnode DIN, DINdel;
  refDATAnode DAN, DANdel;
  refcount c;
  int j;
  
  CopyRect(R,(*newentry).DIR.rect,newrect);     /* see RSTBase.h "hope ..." */
  GetInstChain(R,newrect,depth);
  
  par= &(*R).parameters._;
  
  for (;;) {
    if (depth == (*par).height) {
      M= (*par).dataM; m= (*par).datam;
      
      DAN= &(*(*R).N[depth]).DATA;
      
      if ((*DAN).nofentries < M) {
        (*DAN).entries[(*DAN).nofentries]= (*newentry).DATA; /* 0.. */
        (*DAN).nofentries++;
        (*R).Nmodified[depth]= TRUE;
        
        c= &(*R).count;
        if ((*c).countflag) {
          (*c).datamodifycount++;
        }
        
        depth--;
        AdjustChain(R,depth,newrect);
        break;
      }
    }
    else {
      M= (*par).dirM; m= (*par).dirm;
      
      DIN= &(*(*R).N[depth]).DIR;
      
      if ((*DIN).nofentries < M) {
        (*DIN).entries[(*DIN).nofentries]= (*newentry).DIR; /* 0.. */
        (*DIN).nofentries++;
        (*R).Nmodified[depth]= TRUE;
        
        c= &(*R).count;
        if ((*c).countflag) {
          (*c).dirmodifycount++;
        }
        
        depth--;
        AdjustChain(R,depth,newrect);
        break;
      }
    }
    
    /*** --- LOOP not EXITed by direct insert --- ***/
    
    if ((*R).ReInsert[depth] && depth != 1 && M > 1) {
      if (depth == (*par).height) {
        reinsertqty= (*par).datareinsertqty;
        (*R).NDel[depth]= (refnode)malloc((*R).datanodelen);
        ExcludeDataEntries(R,&(*newentry).DATA,newrect,depth,M,reinsertqty);
      }
      else {
        reinsertqty= (*par).dirreinsertqty;
        (*R).NDel[depth]= (refnode)malloc((*R).dirnodelen);
        ExcludeDirEntries(R,&(*newentry).DIR,newrect,depth,M,reinsertqty);
      }
      heightbefore= (*par).height;
      (*R).ReInsert[depth]= FALSE;
      (*R).ReInsert[depth-1]= TRUE;
      (*R).Nmodified[depth]= TRUE;
      AdjustChainAfterDeletion(R,depth);
      
      if (depth == (*par).height) {
      
        DANdel= &(*(*R).NDel[depth]).DATA;
        
        for (j= reinsertqty; j >= 0; j--) {
          Insert(R,(refentry)&(*DANdel).entries[j],depth);
          if ((*par).height > heightbefore) {
            depth+= (*par).height - heightbefore;
            heightbefore= (*par).height;
          }
        }
      }
      else {
      
        DINdel= &(*(*R).NDel[depth]).DIR;
        
        for (j= reinsertqty; j >= 0; j--) {
          Insert(R,(refentry)&(*DINdel).entries[j],depth);
          if ((*par).height > heightbefore) {
            depth+= (*par).height - heightbefore;
            heightbefore= (*par).height;
          }
        }
      }
      free((*R).NDel[depth]); (*R).NDel[depth]= NULL;
      if (depth != 1) {
        (*R).ReInsert[depth-1]= FALSE;
        /* ReInitialisation for the case, that the LOOP
           did not call ReInsert in the next higher level */
      }
      break;
    }
    else {
      Split(R,newentry,&depth,M,m);
    }
  }
}

/***********************************************************************/

static void ExcludeDataEntries(RSTREE R,
                               typDATAent *newentry,
                               typrect newrect,
                               int depth,
                               int M,
                               int reinsertqty)

{
  refparameters par;
  typpoint allcenter, newrectcenter, center;
  refDATAnode n, nd;
  ValueArray distarr;
  IndexArray I;
  ByteArray chosen;
  int j, k, sortind;
  
  par= &(*R).parameters._;
  
  EvalCenter(R,(*(*R).N[depth-1]).DIR.entries[(*R).E[depth-1]].rect,allcenter);
  EvalCenter(R,newrect,newrectcenter);
  
  n= &(*(*R).N[depth]).DATA;
  
  for (j= 0; j < M; j++) {
    EvalCenter(R,(*n).entries[j].rect,center);
    distarr[j]= RSTDistance(R,allcenter,center);
    I[j]= j;
  }
  distarr[M]= RSTDistance(R,allcenter,newrectcenter);
  I[M]= M;
  QuickSortValArr(0,M,distarr,I);
  
  nd= &(*(*R).NDel[depth]).DATA;
  
  FalseArray((int *)chosen,(*par).dataMwords);
  j= 0; k= M;
  while (j < reinsertqty) {
    sortind= I[k];
    chosen[sortind]= 1;
    if (sortind != M) {
      (*nd).entries[j]= (*n).entries[sortind];
    }
    else {
      (*nd).entries[j]= *newentry;
    }
    j++;
    k--;
  }
  /* entry referenced by distarr[I[M]], [I[M-1]], ...
     becomes NDel[depth]->entries[0], [1], ... */
  
  if (chosen[M] == 1) {
    (*nd).entries[reinsertqty]= (*n).entries[I[M-reinsertqty]];
    chosen[I[M-reinsertqty]]= 1;
  }
  else {
    (*nd).entries[reinsertqty]= *newentry;
  }
  /* entry referenced by I[M-reinsertqty] or newentry (referenced
     by I[?] == 0) becomes NDel[depth]->entries[reinsertqty] */
  /* Entries to reinsert now in NDel[depth] sorted by */
  /* decreasing distances to the midpoint.            */
  /*    !!! newentry is reinserted in any case !!!    */
  
  (*n).nofentries= M - reinsertqty;
  j= 0; k= M - 1;
  do {
    if (chosen[j] == 1) {
      while (chosen[k] == 1) {
        k--;
      }
      (*n).entries[j]= (*n).entries[k];
      chosen[k]= 1;
    }
    j++;
  } while (j < (*n).nofentries);
}

/***********************************************************************/

static void ExcludeDirEntries(RSTREE R,
                              typDIRent *newentry,
                              typrect newrect,
                              int depth,
                              int M,
                              int reinsertqty)

{
  refparameters par;
  typpoint allcenter, newrectcenter, center;
  refDIRnode n, nd;
  ValueArray distarr;
  IndexArray I;
  ByteArray chosen;
  int j, k, sortind;
  
  par= &(*R).parameters._;
  
  EvalCenter(R,(*(*R).N[depth-1]).DIR.entries[(*R).E[depth-1]].rect,allcenter);
  EvalCenter(R,newrect,newrectcenter);
  
  n= &(*(*R).N[depth]).DIR;
  
  for (j= 0; j < M; j++) {
    EvalCenter(R,(*n).entries[j].rect,center);
    distarr[j]= RSTDistance(R,allcenter,center);
    I[j]= j;
  }
  distarr[M]= RSTDistance(R,allcenter,newrectcenter);
  I[M]= M;
  QuickSortValArr(0,M,distarr,I);
  
  nd= &(*(*R).NDel[depth]).DIR;
  
  FalseArray((int *)chosen,(*par).dirMwords);
  j= 0; k= M;
  while (j < reinsertqty) {
    sortind= I[k];
    chosen[sortind]= 1;
    if (sortind != M) {
      (*nd).entries[j]= (*n).entries[sortind];
    }
    else {
      (*nd).entries[j]= *newentry;
    }
    j++;
    k--;
  }
  /* entry referenced by distarr[I[M]], [I[M-1]], ...
     becomes NDel[depth]->entries[0], [1], ... */
  
  if (chosen[M] == 1) {
    (*nd).entries[reinsertqty]= (*n).entries[I[M-reinsertqty]];
    chosen[I[M-reinsertqty]]= 1;
  }
  else {
    (*nd).entries[reinsertqty]= *newentry;
  }
  /* entry referenced by I[M-reinsertqty] or newentry (referenced
     by I[?] == 0) becomes NDel[depth]->entries[reinsertqty] */
  /* Entries to reinsert now in NDel[depth] sorted by */
  /* decreasing distances to the midpoint.            */
  /*    !!! newentry is reinserted in any case !!!    */
  
  (*n).nofentries= M - reinsertqty;
  j= 0; k= M - 1;
  do {
    if (chosen[j] == 1) {
      while (chosen[k] == 1) {
        k--;
      }
      (*n).entries[j]= (*n).entries[k];
      chosen[k]= 1;
    }
    j++;
  } while (j < (*n).nofentries);
}

/***********************************************************************/

static void Split(RSTREE R,
                  typentry *newentry,
                  int *depth,
                  int M,
                  int m)

{
  refparameters par;
  refcount c;
  int pagenr;
  boolean isdata;
  int i;
  
  par= &(*R).parameters._;
  
  if (*depth == (*par).height) {
    if (M == 1) {
      UnDistributData(R,&(*newentry).DATA);
    }
    else
    {
      SplitAndDistributData(R,*depth,&(*newentry).DATA,M,m);
    }
  }
  else {
    SplitAndDistributDir(R,*depth,&(*newentry).DIR,M,m);
  }
  if (*depth == 1) {
    (*par).height++;
    *depth= 2;
    for (i= (*par).height; i >= 2; i--) {
      (*R).N[i]= (*R).N[i-1];
      (*R).NDel[i]= (*R).NDel[i-1];         /* for ReInsert */
      (*R).ReInsert[i]= (*R).ReInsert[i-1]; /* for ReInsert */
      /*    following assignments not necessary for i=2: */
      (*R).P[i]= (*R).P[i-1];                      /* set(i=2) */
      (*R).Nmodified[i]= (*R).Nmodified[i-1];      /* set(i=2) */
      (*par).pagecountarr[i]= (*par).pagecountarr[i-1];
    }
    (*R).N[1]= (refnode)malloc((*R).dirnodelen);
    /* R->P[1] is RSTBase.rootblocknumb forever */
    /* R->Nmodified[1] set inserting the sibling */
    (*R).NDel[1]= NULL;                /* for ReInsert */
    (*R).ReInsert[1]= FALSE;           /* for ReInsert */
    /* initiate first root entry: */
    (*(*R).N[1]).DIR.nofentries= 1;
    (*R).E[1]= 0;
    GetPageNr(R,&pagenr,*depth);
    (*(*R).N[1]).DIR.entries[0].ptrtosub= pagenr;
    (*R).P[2]= pagenr;
  }
  isdata= *depth == (*par).height;
  (*R).Nmodified[*depth]= TRUE;
  GetPageNr(R,&pagenr,*depth);
  (*newentry).DIR.ptrtosub= pagenr;

  if (isdata) {
    EvalDataEnclRect(R,&(*(*R).Nsibling).DATA,(*newentry).DIR.rect);
    EvalDataEnclRect(R,&(*(*R).N[*depth]).DATA,
                       (*(*R).N[*depth-1]).DIR.entries[(*R).E[*depth-1]].rect);
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datamodifycount+= 2;
    }
  }
  else {
    EvalDirEnclRect(R,&(*(*R).Nsibling).DIR,(*newentry).DIR.rect);
    EvalDirEnclRect(R,&(*(*R).N[*depth]).DIR,
                       (*(*R).N[*depth-1]).DIR.entries[(*R).E[*depth-1]].rect);
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirmodifycount+= 2;
    }
  }
  PutNode(R,(*R).Nsibling,pagenr,*depth);
  (*depth)--;
}

/***********************************************************************/

static void SplitAndDistributDir(RSTREE R,
                                 int depth,
                                 typDIRent *newentry,
                                 int M,
                                 int m)

{
  ValueArray edgearray, spacearray;
  RectArray rectarray;
  refparameters par;
  refDIRnode n;
  refinterval re;
  
  IndexArray I;
  Side sortside, currsortside, axsortside;
  int axis;
  int mlow, mhigh,
      index, backind,
      axisindex;
  typrect leftrect, rightrect;
  double leftspace, rightspace,
         leftedge, rightedge, edge, validedge, backvaledge, axisedge,
         ovlp,
         value, validvalue, backvalvalue;
  int i, j, ii;
  
  par= &(*R).parameters._;
  
  (*R).Ntosplit= (*R).N[depth];
  (*R).N[depth]= (*R).helpdirnode;
  (*R).helpdirnode= (*R).Ntosplit;
  for (j= 0; j <= M; j++) {
    I[j]= j;
  }
  mhigh= M-m+1;
  mlow= m-1;
  
  n= &(*(*R).Ntosplit).DIR;
  
  (*n).entries[M]= *newentry;
  for (i= 0; i <= (*par).maxdim; i++) {
    sortside= low;
    currsortside= low;
    QuickSortDirEnt(0,M,i,low,(*n).entries,I);
    CopyRect(R,(*n).entries[I[M]].rect,rightrect);
    rightrect[i].l= (*n).entries[I[mhigh]].rect[i].l;
    for (j= M-1; j >= m; j--) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < rightrect[ii].l) {
          rightrect[ii].l= (*re).l;
        }
        if ((*re).h > rightrect[ii].h) {
          rightrect[ii].h= (*re).h;
        }
      }
      if (j <= mhigh) {
        rightedge= 0.0; rightspace= 1.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          rightedge+= rightrect[ii].h - rightrect[ii].l;
          rightspace*= rightrect[ii].h - rightrect[ii].l;
        }
        CopyRect(R,rightrect,rectarray[j]);
        edgearray[j]= rightedge;
        spacearray[j]= rightspace;
      }
    }
    /* fill array from M-m+1 down to m
       with rightrect, rightedge, rightspace */
    CopyRect(R,(*n).entries[I[0]].rect,leftrect);
    for (j= 1; j < mhigh; j++) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < leftrect[ii].l) {
          leftrect[ii].l= (*re).l;
        }
        if ((*re).h > leftrect[ii].h) {
          leftrect[ii].h= (*re).h;
        }
      }
      if (j >= mlow) {
        leftedge= 0.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          leftedge+= leftrect[ii].h - leftrect[ii].l;
        }
        if (Overlaps(R,leftrect,rectarray[j+1])) {
          GetOverlap(R,leftrect,rectarray[j+1],&ovlp);
          if (j == mlow) {
            validedge= leftedge+edgearray[j+1];
            validvalue= ovlp;
            index= m;
          }
          else {
            edge= leftedge+edgearray[j+1];
            validedge= validedge+edge;
            value= ovlp;
            if (value < validvalue) {
              validvalue= value;
              index= j+1;
            }
          }
        }
        else {
          leftspace= 1.0;
          for (ii= 0; ii <= (*par).maxdim; ii++) {
            leftspace*= leftrect[ii].h - leftrect[ii].l;
          }
          if (j == mlow) {
            validedge= leftedge+edgearray[j+1];
            validvalue= -1.0/(leftspace+spacearray[j+1]);
            index= m;
          }
          else {
            edge= leftedge+edgearray[j+1];
            validedge= validedge+edge;
            value= -1.0/(leftspace+spacearray[j+1]);
            if (value < validvalue) {
              validvalue= value;
              index= j+1;
            }
          }
        }
      }
    }
    /* Compute sum of the edges
       from leftedge[m-1]+rightedge[m]
       up to leftedge[M-m]+rightedge[M-m+1], -> axis;
       minimal area resp. overlap, -> index. */
    currsortside= high;
    QuickSortDirEnt(0,M,i,high,(*n).entries,I);
    CopyRect(R,(*n).entries[I[0]].rect,leftrect);
    leftrect[i].h= (*n).entries[I[mlow]].rect[i].h;
    for (j= 1; j < mhigh; j++) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < leftrect[ii].l) {
          leftrect[ii].l= (*re).l;
        }
        if ((*re).h > leftrect[ii].h) {
            leftrect[ii].h= (*re).h;
        }
      }
      if (j >= mlow) {
        leftedge= 0.0; leftspace= 1.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          leftedge+= leftrect[ii].h - leftrect[ii].l;
          leftspace*= leftrect[ii].h - leftrect[ii].l;
        }
        CopyRect(R,leftrect,rectarray[j]);
        edgearray[j]= leftedge;
        spacearray[j]= leftspace;
      }
    }
    /* fill array from m-1 up to M-m
       with leftrect, leftedge, leftspace */
    CopyRect(R,(*n).entries[I[M]].rect,rightrect);
    for (j= M-1; j >= m; j--) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re = &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < rightrect[ii].l) {
          rightrect[ii].l= (*re).l;
        }
        if ((*re).h > rightrect[ii].h) {
          rightrect[ii].h= (*re).h;
        }
      }
      if (j <= mhigh) {
        rightedge= 0.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          rightedge+= rightrect[ii].h- rightrect[ii].l;
        }
        if (Overlaps(R,rightrect,rectarray[j-1])) {
          GetOverlap(R,rightrect,rectarray[j-1],&ovlp);
          if (j == mhigh) {
            backvaledge= rightedge+edgearray[j-1];
            backvalvalue= ovlp;
            backind= j;
          }
          else {
            edge= rightedge+edgearray[j-1];
            backvaledge= backvaledge+edge;
            value= ovlp;
            if (value < backvalvalue) {
              backvalvalue= value;
              backind= j;
            }
          }
        }
        else {
          rightspace= 1.0;
          for (ii= 0; ii <= (*par).maxdim; ii++) {
            rightspace*= rightrect[ii].h - rightrect[ii].l;
          }
          if (j == mhigh) {
            backvaledge= rightedge+edgearray[j-1];
            backvalvalue= -1.0/(rightspace+spacearray[j-1]);
            backind= j;
          }
          else {
            edge= rightedge+edgearray[j-1];
            backvaledge= backvaledge+edge;
            value= -1.0/(rightspace+spacearray[j-1]);
            if (value < backvalvalue) {
              backvalvalue= value;
              backind= j;
            }
          }
        }
      }
    }
    /* Compute sum of the edges
       from rightedge[M-m+1]+leftedge[M-m]
       down to rightedge[m] +leftedge[m-1], -> axis;
       minimal area resp. overlap, -> index. */
    if (backvalvalue <= validvalue) {
      sortside= high;
      index= backind;
    }
    validedge= validedge+backvaledge;
    if (i == 0) {
      axsortside= sortside;
      axis= i;
      axisedge= validedge;
      axisindex= index;
    }
    else if (validedge < axisedge) {
      axsortside= sortside;
      axis= i;
      axisedge= validedge;
      axisindex= index;
    }
    
  }
  if (axis != (*par).maxdim || axsortside != currsortside) {
    QuickSortDirEnt(0,M,axis,axsortside,(*n).entries,I);
  }
  
  n= &(*(*R).N[depth]).DIR;
  
  (*n).nofentries= axisindex;
  for (j= 0; j < (*n).nofentries; j++) {
    (*n).entries[j]= (*(*R).Ntosplit).DIR.entries[I[j]];
  }
  /* printf("%3d",axisindex); */
  /* N[depth] gets entries 0 up to axisindex-1 */
  
  n= &(*(*R).Nsibling).DIR;
  
  (*n).nofentries= M-axisindex+1;
  for (j= 0; j < (*n).nofentries; j++) {
    (*n).entries[j]= (*(*R).Ntosplit).DIR.entries[I[axisindex+j]];
  }
  /* printf("%3d",M-axisindex+1); */
  /* printf("%3d\n",M+1); */
}

/***********************************************************************/

static void SplitAndDistributData(RSTREE R,
                                  int depth,
                                  typDATAent *newentry,
                                  int M,
                                  int m)

{
  ValueArray edgearray, spacearray;
  RectArray rectarray;
  refparameters par;
  refDATAnode n;
  refinterval re;
  
  IndexArray I;
  Side sortside, currsortside, axsortside;
  int axis;
  int mlow, mhigh,
      index, backind,
      axisindex;
  typrect leftrect, rightrect;
  double leftspace, rightspace,
         leftedge, rightedge, edge, validedge, backvaledge, axisedge,
         ovlp,
         value, validvalue, backvalvalue;
  int i, j, ii;
  
  par= &(*R).parameters._;
  
  (*R).Ntosplit= (*R).N[depth];
  (*R).N[depth]= (*R).helpdatanode;
  (*R).helpdatanode= (*R).Ntosplit;
  for (j= 0; j <= M; j++) {
    I[j]= j;
  }
  mhigh= M-m+1;
  mlow= m-1;
  
  n= &(*(*R).Ntosplit).DATA;
  
  (*n).entries[M]= *newentry;
  for (i= 0; i <= (*par).maxdim; i++) {
    sortside= low;
    currsortside= low;
    QuickSortDataEnt(0,M,i,low,(*n).entries,I);
    CopyRect(R,(*n).entries[I[M]].rect,rightrect);
    rightrect[i].l= (*n).entries[I[mhigh]].rect[i].l;
    for (j= M-1; j >= m; j--) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < rightrect[ii].l) {
          rightrect[ii].l= (*re).l;
        }
        if ((*re).h > rightrect[ii].h) {
          rightrect[ii].h= (*re).h;
        }
      }
      if (j <= mhigh) {
        rightedge= 0.0; rightspace= 1.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          rightedge+= rightrect[ii].h - rightrect[ii].l;
          rightspace*= rightrect[ii].h - rightrect[ii].l;
        }
        CopyRect(R,rightrect,rectarray[j]);
        edgearray[j]= rightedge;
        spacearray[j]= rightspace;
      }
    }
    /* fill array from M-m+1 down to m
       with rightrect, rightedge, rightspace */
    CopyRect(R,(*n).entries[I[0]].rect,leftrect);
    for (j= 1; j < mhigh; j++) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < leftrect[ii].l) {
          leftrect[ii].l= (*re).l;
        }
        if ((*re).h > leftrect[ii].h) {
          leftrect[ii].h= (*re).h;
        }
      }
      if (j >= mlow) {
        leftedge= 0.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          leftedge+= leftrect[ii].h - leftrect[ii].l;
        }
        if (Overlaps(R,leftrect,rectarray[j+1])) {
          GetOverlap(R,leftrect,rectarray[j+1],&ovlp);
          if (j == mlow) {
            validedge= leftedge+edgearray[j+1];
            validvalue= ovlp;
            index= m;
          }
          else {
            edge= leftedge+edgearray[j+1];
            validedge= validedge+edge;
            value= ovlp;
            if (value < validvalue) {
              validvalue= value;
              index= j+1;
            }
          }
        }
        else {
          leftspace= 1.0;
          for (ii= 0; ii <= (*par).maxdim; ii++) {
            leftspace*= leftrect[ii].h - leftrect[ii].l;
          }
          if (j == mlow) {
            validedge= leftedge+edgearray[j+1];
            validvalue= -1.0/(leftspace+spacearray[j+1]);
            index= m;
          }
          else {
            edge= leftedge+edgearray[j+1];
            validedge= validedge+edge;
            value= -1.0/(leftspace+spacearray[j+1]);
            if (value < validvalue) {
              validvalue= value;
              index= j+1;
            }
          }
        }
      }
    }
    /* Compute sum of the edges
       from leftedge[m-1]+rightedge[m]
       up to leftedge[M-m]+rightedge[M-m+1], -> axis;
       minimal area resp. overlap, -> index. */
    currsortside= high;
    QuickSortDataEnt(0,M,i,high,(*n).entries,I);
    CopyRect(R,(*n).entries[I[0]].rect,leftrect);
    leftrect[i].h= (*n).entries[I[mlow]].rect[i].h;
    for (j= 1; j < mhigh; j++) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re= &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < leftrect[ii].l) {
          leftrect[ii].l= (*re).l;
        }
        if ((*re).h > leftrect[ii].h) {
            leftrect[ii].h= (*re).h;
        }
      }
      if (j >= mlow) {
        leftedge= 0.0; leftspace= 1.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          leftedge+= leftrect[ii].h - leftrect[ii].l;
          leftspace*= leftrect[ii].h - leftrect[ii].l;
        }
        CopyRect(R,leftrect,rectarray[j]);
        edgearray[j]= leftedge;
        spacearray[j]= leftspace;
      }
    }
    /* fill array from m-1 up to M-m
       with leftrect, leftedge, leftspace */
    CopyRect(R,(*n).entries[I[M]].rect,rightrect);
    for (j= M-1; j >= m; j--) {
      for (ii= 0; ii <= (*par).maxdim; ii++) {
        
        re = &(*n).entries[I[j]].rect[ii];
        
        if ((*re).l < rightrect[ii].l) {
          rightrect[ii].l= (*re).l;
        }
        if ((*re).h > rightrect[ii].h) {
          rightrect[ii].h= (*re).h;
        }
      }
      if (j <= mhigh) {
        rightedge= 0.0;
        for (ii= 0; ii <= (*par).maxdim; ii++) {
          rightedge+= rightrect[ii].h- rightrect[ii].l;
        }
        if (Overlaps(R,rightrect,rectarray[j-1])) {
          GetOverlap(R,rightrect,rectarray[j-1],&ovlp);
          if (j == mhigh) {
            backvaledge= rightedge+edgearray[j-1];
            backvalvalue= ovlp;
            backind= j;
          }
          else {
            edge= rightedge+edgearray[j-1];
            backvaledge= backvaledge+edge;
            value= ovlp;
            if (value < backvalvalue) {
              backvalvalue= value;
              backind= j;
            }
          }
        }
        else {
          rightspace= 1.0;
          for (ii= 0; ii <= (*par).maxdim; ii++) {
            rightspace*= rightrect[ii].h - rightrect[ii].l;
          }
          if (j == mhigh) {
            backvaledge= rightedge+edgearray[j-1];
            backvalvalue= -1.0/(rightspace+spacearray[j-1]);
            backind= j;
          }
          else {
            edge= rightedge+edgearray[j-1];
            backvaledge= backvaledge+edge;
            value= -1.0/(rightspace+spacearray[j-1]);
            if (value < backvalvalue) {
              backvalvalue= value;
              backind= j;
            }
          }
        }
      }
    }
    /* Compute sum of the edges
       from rightedge[M-m+1]+leftedge[M-m]
       down to rightedge[m] +leftedge[m-1], -> axis;
       minimal area resp. overlap, -> index. */
    if (backvalvalue <= validvalue) {
      sortside= high;
      index= backind;
    }
    validedge= validedge+backvaledge;
    if (i == 0) {
      axsortside= sortside;
      axis= i;
      axisedge= validedge;
      axisindex= index;
    }
    else if (validedge < axisedge) {
      axsortside= sortside;
      axis= i;
      axisedge= validedge;
      axisindex= index;
    }
    
  }
  if (axis != (*par).maxdim || axsortside != currsortside) {
    QuickSortDataEnt(0,M,axis,axsortside,(*n).entries,I);
  }
  
  n= &(*(*R).N[depth]).DATA;
  
  (*n).nofentries= axisindex;
  for (j= 0; j < (*n).nofentries; j++) {
    (*n).entries[j]= (*(*R).Ntosplit).DATA.entries[I[j]];
  }
  /* printf("%3d",axisindex); */
  /* N[depth] gets entries 0 up to axisindex-1 */
  
  n= &(*(*R).Nsibling).DATA;
  
  (*n).nofentries= M-axisindex+1;
  for (j= 0; j < (*n).nofentries; j++) {
    (*n).entries[j]= (*(*R).Ntosplit).DATA.entries[I[axisindex+j]];
  }
  /* printf("%3d",M-axisindex+1); */
  /* printf("%3d\n",M+1); */
}

/***********************************************************************/

static void UnDistributData(RSTREE R,
                            typDATAent *newentry)

{
  refDATAnode DAN;
  
  DAN= &(*(*R).Nsibling).DATA;
  
  (*DAN).nofentries= 1;
  (*DAN).entries[0]= *newentry;
  /* printf(" Un\n"); */
}

/***********************************************************************/

void DeleteOneRec(RSTREE R)

{
  refparameters par;
  refDIRnode DIN;
  refDATAnode DAN;
  refcount c;
  int depth, heightbefore;
  int i, j;
  
  par = &(*R).parameters._;
      
  depth= (*par).height;
  for (;;) {
    if (depth == (*par).height) {
    
      DAN= &(*(*R).N[depth]).DATA;
      
      (*DAN).nofentries--;
      (*DAN).entries[(*R).E[depth]]= (*DAN).entries[(*DAN).nofentries];
      if ((*DAN).nofentries >= (*par).datam || depth == 1) {
        (*R).Nmodified[depth]= TRUE;
        
        c= &(*R).count;
        if ((*c).countflag) {
          (*c).datamodifycount++;
        }
          
        AdjustChainAfterDeletion(R,depth);
        break;
      }
    }
    else {
    
      DIN= &(*(*R).N[depth]).DIR;
      
      (*DIN).nofentries--;
      (*DIN).entries[(*R).E[depth]]= (*DIN).entries[(*DIN).nofentries];
      if ((*DIN).nofentries >= (*par).dirm || depth == 1) {
        (*R).Nmodified[depth]= TRUE;
        
        c= &(*R).count;
        if ((*c).countflag) {
          (*c).dirmodifycount++;
        }
        
        AdjustChainAfterDeletion(R,depth);
        break;
      }
    }
    /*** --- LOOP not EXITed by deletion without underflow --- ***/
    NtoNDel(R,depth);
    depth--;
  }
  /*** --- ReInsertion --- ***/
  i= 2; heightbefore= (*par).height;
  while (i <= (*par).height) {
    if ((*R).NDel[i] != NULL) {
    
      c= &(*R).count;
      if ((*c).countflag) {
        if (i == (*par).height) {
          (*c).datavisitcount++;
        }
        else {
          (*c).dirvisitcount++;
        }
      }
      if (i == (*par).height) {
      
        DAN= &(*(*R).NDel[i]).DATA;
        
        for (j= 0; j < (*DAN).nofentries; j++) {
          Insert(R,(refentry)&(*DAN).entries[j],i);
          if ((*par).height > heightbefore) {
            i++; heightbefore++;
          }
        }
      }
      else {
      
        DIN= &(*(*R).NDel[i]).DIR;
        
        for (j= 0; j < (*DIN).nofentries; j++) {
          Insert(R,(refentry)&(*DIN).entries[j],i);
          if ((*par).height > heightbefore) {
            i++; heightbefore++;
          }
        }
      }
      free((*R).NDel[i]); (*R).NDel[i]= NULL;
    }
    i++;
  }
  /*** --- Shrink the tree --- ***/
  if ((*par).height != 1 && (*(*R).N[1]).DIR.nofentries == 1) {
    ShrinkTree(R);
  }
}

/***********************************************************************/

static void NtoNDel(RSTREE R,
                    int depth)

{
  refparameters par;
  
  par= &(*R).parameters._;
  
  if ((*(*R).N[depth]).DIR.nofentries != 0) {     /* see RSTBase.h "hope ..." */
    (*R).NDel[depth]= (*R).N[depth];
    if (depth == (*par).height) {
      (*R).N[depth]= (refnode)malloc((*R).datanodelen);
    }
    else {
      (*R).N[depth]= (refnode)malloc((*R).dirnodelen);
    }
  }
  else {
    /* this is only possible for data nodes. NDel[height] stays NULL. */
    /* ==> no ReInsertion(height). P[height] set to 0 in any case.    */
    /* if height == 2 and shrink is done NewNode(R,2) is called.      */
  }
  PutPageNr(R,(*R).P[depth],depth);
  (*R).P[depth]= 0;
  (*R).Nmodified[depth]= FALSE;
}

/***********************************************************************/

static void ShrinkTree(RSTREE R)

{
  refparameters par;
  refcount c;
  int i;
  
  par= &(*R).parameters._;
  
  if ((*R).P[2] == 0) {
    (*R).E[1]= 0;
    NewNode(R,2);
  }
  free((*R).N[1]);
  for (i= 1; i <= (*par).height; i++) {
    (*R).N[i]= (*R).N[i+1];
  }
  (*R).Nmodified[1]= TRUE;
  
  c= &(*R).count;
  if ((*c).countflag) {
    (*c).dirmodifycount++;
  }
  
  PutPageNr(R,(*R).P[2],2);
  for (i= 2; i <= (*par).height; i++) {
    (*R).P[i]= (*R).P[i+1];
    (*R).Nmodified[i]= (*R).Nmodified[i+1];
    (*par).pagecountarr[i]= (*par).pagecountarr[i+1];
  }
  (*R).E[(*par).height]= -1;
  (*par).height--;
}

/***********************************************************************/
