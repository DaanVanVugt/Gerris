/* ----- RSTInOut.c ----- */


#include "RStarTree.h"
#include "RSTInOut.h"


/* declarations */

/************************************************************************/

void PutNode(RSTREE R,
             refnode nodeptr,
             int pagenr,
             int depth)

{
  refparameters par;
  refcount c;
  
  par= &(*R).parameters._;
  
  if (depth == (*par).height) {
    WritePage(R,(*R).data,pagenr,nodeptr);
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datawritecount++;
    }
  }
  else {
    WritePage(R,(*R).dir,pagenr,nodeptr);
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirwritecount++;
    }
  }
}
    
/************************************************************************/

void GetNode(RSTREE R,
             refnode nodeptr,
             int pagenr,
             int depth)

{
  refparameters par;
  refcount c;
  
  par= &(*R).parameters._;
  
  if (depth == (*par).height) {
    ReadPage(R,(*R).data,pagenr,nodeptr);
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datareadcount++;
    }
  }
  else {
    ReadPage(R,(*R).dir,pagenr,nodeptr);
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirreadcount++;
    }
  }
}

/************************************************************************/

void NewNode(RSTREE R,
             int depth)

{
  if ((*R).Nmodified[depth]) {
    PutNode(R,(*R).N[depth],(*R).P[depth],depth);
    (*R).Nmodified[depth]= FALSE;
  }
  (*R).P[depth]= (*(*R).N[depth-1]).DIR.entries[(*R).E[depth-1]].ptrtosub;
  GetNode(R,(*R).N[depth],(*R).P[depth],depth);
}

/************************************************************************/

void ReadPage(RSTREE R,
              typfiledesc fd,
              int pagenr,
              void *block)

{
  int pnb; /* position or nbytes */
  
  pnb= lseek(fd.f,(off_t)pagenr*(off_t)fd.bl,SEEK_SET);
  if (pnb != -1) {
    pnb= read(fd.f,block,fd.bl);
  }
  if (pnb <= 0) {
    (*R).RSTDone= FALSE;
  }
}

/************************************************************************/

void WritePage(RSTREE R,
               typfiledesc fd,
               int pagenr,
               void *block)

{
  int pnb; /* position or nbytes */
  
  pnb= lseek(fd.f,(off_t)pagenr*(off_t)fd.bl,SEEK_SET);
  if (pnb != -1) {
    pnb= write(fd.f,block,fd.bl);
  }
  if (pnb <= 0) {
    (*R).RSTDone= FALSE;
  }
}

/************************************************************************/

void GetPageNr(RSTREE R, int *pagenr, int depth)

{ 
  refparameters par;
  refpagedir dpd;
  
  par= &(*R).parameters._;
  
  if (depth == (*par).height) {
  
    dpd= &(*R).datapagedir._;
    
    if ((*dpd).nofnumbers == 0) {
      if ((*dpd).childnr == firstPDblocknumb) {
        (*dpd).number[0]++;
        *pagenr= (*dpd).number[0];
      }
      else {
        ReadPage(R,(*R).dataPD,(*dpd).childnr,&(*R).datapagedir._);
        (*dpd).childnr--;
        *pagenr= (*dpd).number[maxnum];
        (*dpd).nofnumbers= maxnum - 1;
      }
    }
    else {
      *pagenr= (*dpd).number[(*dpd).nofnumbers];
      (*dpd).nofnumbers--;
    }
    (*par).datapagecount++;
  }
  else {
  
    dpd= &(*R).dirpagedir._;
    
    if ((*dpd).nofnumbers == 0) {
      if ((*dpd).childnr == firstPDblocknumb) {
        (*dpd).number[0]++;
        *pagenr= (*dpd).number[0];
      }
      else {
        ReadPage(R,(*R).dirPD,(*dpd).childnr,&(*R).dirpagedir._);
        (*dpd).childnr--;
        *pagenr= (*dpd).number[maxnum];
        (*dpd).nofnumbers= maxnum - 1;
      }
    }
    else {
      *pagenr= (*dpd).number[(*dpd).nofnumbers];
      (*dpd).nofnumbers--;
    }
    (*par).dirpagecount++;
  }
  (*par).pagecountarr[depth]++;
}

/************************************************************************/

void PutPageNr(RSTREE R, int pagenr, int depth)

{
  refparameters par;
  refpagedir dpd;
  
  par= &(*R).parameters._;
  
  if (depth == (*par).height) {
  
    dpd= &(*R).datapagedir._;
    
    if ((*dpd).nofnumbers == maxnum) {
      (*dpd).childnr++;
      WritePage(R,(*R).dataPD,(*dpd).childnr,&(*R).datapagedir._);
      (*dpd).nofnumbers= 1;
      (*dpd).number[1]= pagenr;
    }
    else {
      (*dpd).nofnumbers++;
      (*dpd).number[(*dpd).nofnumbers]= pagenr;
    }
    (*par).datapagecount--;
  }
  else {
  
    dpd= &(*R).dirpagedir._;
    
    if ((*dpd).nofnumbers == maxnum) {
      (*dpd).childnr++;
      WritePage(R,(*R).dirPD,(*dpd).childnr,&(*R).dirpagedir._);
      (*dpd).nofnumbers= 1;
      (*dpd).number[1]= pagenr;
    }
    else {
      (*dpd).nofnumbers++;
      (*dpd).number[(*dpd).nofnumbers]= pagenr;
    }
    (*par).dirpagecount--;
  }
  (*par).pagecountarr[depth]--;
}

/************************************************************************/
