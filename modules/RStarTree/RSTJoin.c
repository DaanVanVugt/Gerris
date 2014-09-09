/* ----- RSTJoin.c ----- */


#include "RStarTree.h"
#include "RSTJoin.h"
#include "RSTUtil.h"
#include "RSTInOut.h"


/* declarations */

static void JnRgnCntNv(RSTREE R, RSTREE Rx,
                       boolean order,
                       int depth,
                       typrect dataqueryrect1,
                       typrect dataqueryrect2,
                       typrect joinrect,
                       DataQueryProc DataQuery,
                       DirQueryProc DirJoin,
                       DataQueryProc DataJoin,
                       int *keysqualifying);

static void JnRgnQueryNv(RSTREE R, RSTREE Rx,
                         boolean order,
                         int depth,
                         typrect dataqueryrect1,
                         typrect dataqueryrect2,
                         typrect joinrect,
                         refinfo ptr_to_info_of_other_tree,
                         DataQueryProc DataQuery,
                         DirQueryProc DirJoin,
                         DataQueryProc DataJoin,
                         JoinManageProc Manage,
                         void *buf1,
                         void *buf2,
                         boolean *finish);


/************************************************************************/

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
             int *mark)

{
  refDATAnode n;
  refcount c;
  boolean istoread1, istoread2;
  int downqualifying;
  typrect unused;
  int i, j;
  
  boolean verbosejoincount= TRUE;
  
  if (depth != (*R1).parameters._.height && depth !=
                                                   (*R2).parameters._.height) {
    for (i= 0; i < (*(*R1).N[depth]).DIR.nofentries; i++) {
      for (j= 0; j < (*(*R2).N[depth]).DIR.nofentries; j++) {
        if ( Dir1Query(R1,(*(*R1).N[depth]).DIR.entries[i].rect,
                                                  R1rectangle1,R1rectangle2) &&
             Dir2Query(R2,(*(*R2).N[depth]).DIR.entries[j].rect,
                                                 R2rectangle1,R2rectangle2) ) {
          if (DirJoin(R1,(*(*R1).N[depth]).DIR.entries[i].rect,
                               (*(*R2).N[depth]).DIR.entries[j].rect,unused)) {
          /* R1 sets the number of dimensions and type */
            (*R1).E[depth]= i; (*R2).E[depth]= j;
            istoread1= (*(*R1).N[depth]).DIR.entries[i].ptrtosub !=
                                                              (*R1).P[depth+1];
            if ( istoread1 ) {
              NewNode(R1,depth+1);
            }
            istoread2= (*(*R2).N[depth]).DIR.entries[j].ptrtosub !=
                                                              (*R2).P[depth+1];
            if ( istoread2 ) {
              NewNode(R2,depth+1);
            }
            JnCntNv(R1,R2,
                    depth+1,
                    R1rectangle1,R1rectangle2,R2rectangle1,R2rectangle2,
                    Dir1Query,Data1Query,Dir2Query,Data2Query,DirJoin,DataJoin,
                    keysqualifying,
                    mark);
          }
        }
      }
    }
    
    c= &(*R1).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
    
    c= &(*R2).count;
    if ( (*c).countflag ) {
      (*c).dirvisitcount++;
    }
  }
  else {
    if (depth == (*R1).parameters._.height) {
    
      n= &(*(*R1).N[depth]).DATA;
      
      for (i= 0; i < (*n).nofentries; i++) {
        downqualifying= 0;
        if (Data1Query(R1,(*n).entries[i].rect,R1rectangle1,R1rectangle2)) {
          JnRgnCntNv(R2,R1,FALSE,depth,
                     R2rectangle1,R2rectangle2,
                     (*n).entries[i].rect,
                     Data2Query,
                     DirJoin,DataJoin,
                     &downqualifying);
          *keysqualifying+= downqualifying;
          if (verbosejoincount) {
            if (*keysqualifying > *mark) {
              printf("%s%10d%s\n","More than",*mark," record pairs.");
              if (*mark < 1000) {
                *mark+= 100;
              }
              else if (*mark < 10000) {
                *mark+= 1000;
              }
              else {
                *mark+= 10000;
              }
            }
          }
        }
      }
        
      c= &(*R1).count;
      if ((*c).countflag) {
        (*c).datavisitcount++;
      }
    }
    else {
      
      n= &(*(*R2).N[depth]).DATA;
        
      for (i= 0; i < (*n).nofentries; i++) {
        downqualifying= 0;
        if (Data2Query(R2,(*n).entries[i].rect,R2rectangle1,R2rectangle2)) {
          JnRgnCntNv(R1,R2,TRUE,depth,
                     R1rectangle1,R1rectangle2,
                     (*n).entries[i].rect,
                     Data1Query,
                     DirJoin,DataJoin,
                     &downqualifying);
          *keysqualifying+= downqualifying;
          if (verbosejoincount) {
            if (*keysqualifying > *mark) {
              printf("%s%10d%s\n","More than",*mark," record pairs.");
              if (*mark < 1000) {
                *mark+= 100;
              }
              else if (*mark < 10000) {
                *mark+= 1000;
              }
              else {
                *mark+= 10000;
              }
            }
          }
        }
      }
        
      c= &(*R2).count;
      if ((*c).countflag) {
        (*c).datavisitcount++;
      }
        
    }
  }
}

/************************************************************************/

static void JnRgnCntNv(RSTREE R, RSTREE Rx,
                       boolean order,
                       int depth,
                       typrect dataqueryrect1,
                       typrect dataqueryrect2,
                       typrect joinrect,
                       DataQueryProc DataQuery,
                       DirQueryProc DirJoin,
                       DataQueryProc DataJoin,
                       int *keysqualifying)

{
  refDIRnode DIN;
  refDATAnode DAN;
  refcount c;
  boolean istoread;
  typrect unused;
  int i;
  
  
  if (depth != (*R).parameters._.height) {
  
    DIN= &(*(*R).N[depth]).DIR;

    if (order) {
      for (i= 0; i < (*DIN).nofentries; i++) {
        if (DirJoin(R,(*DIN).entries[i].rect,joinrect,unused)) {
          (*R).E[depth]= i;
          istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
          if ( istoread ) {
            NewNode(R,depth+1);
          }
          JnRgnCntNv(R,Rx,order,depth+1,
                     dataqueryrect1,dataqueryrect2,joinrect,
                     DataQuery,DirJoin,DataJoin,
                     keysqualifying);
        }
      }
    }
    else {
      for (i= 0; i < (*DIN).nofentries; i++) {
        if (DirJoin(Rx,joinrect,(*DIN).entries[i].rect,unused)) {
          (*R).E[depth]= i;
          istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
          if (istoread) {
            NewNode(R,depth+1);
          }
          JnRgnCntNv(R,Rx,order,depth+1,
                     dataqueryrect1,dataqueryrect2,joinrect,
                     DataQuery,DirJoin,DataJoin,
                     keysqualifying);
        }
      }
    }
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
  }
  else {
  
    DAN= &(*(*R).N[depth]).DATA;
    
    if (order) {
      for (i= 0; i < (*DAN).nofentries; i++) {
        if (DataQuery(R,(*DAN).entries[i].rect,dataqueryrect1,dataqueryrect2)) {
          if (DataJoin(R,(*DAN).entries[i].rect,joinrect,unused)) {
            (*R).E[depth]= i;
            (*keysqualifying)++;
          }
        }
      }
    }
    else {
      for (i= 0; i < (*DAN).nofentries; i++) {
        if (DataQuery(R,(*DAN).entries[i].rect,dataqueryrect1,dataqueryrect2)) {
          if (DataJoin(Rx,joinrect,(*DAN).entries[i].rect,unused)) {
            (*R).E[depth]= i;
            (*keysqualifying)++;
          }
        }
      }
    }
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
    
  }
}

/************************************************************************/

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
          boolean *finish)

{
  refDATAnode n;
  refcount c;
  boolean istoread1, istoread2;
  typrect unused;
  int i, j;
    
  if (depth != (*R1).parameters._.height && depth !=
                                                   (*R2).parameters._.height) {
    c= &(*R1).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
    c= &(*R2).count;
    if ( (*c).countflag ) {
      (*c).dirvisitcount++;
    }
    
    for (i= 0; i < (*(*R1).N[depth]).DIR.nofentries; i++) {
      if (*finish) {return;}
      for (j= 0; j < (*(*R2).N[depth]).DIR.nofentries; j++) {
        if (*finish) {return;}
        if ( Dir1Query(R1,(*(*R1).N[depth]).DIR.entries[i].rect,
                                                  R1rectangle1,R1rectangle2) &&
             Dir2Query(R2,(*(*R2).N[depth]).DIR.entries[j].rect,
                                                 R2rectangle1,R2rectangle2) ) {
          if (DirJoin(R1,(*(*R1).N[depth]).DIR.entries[i].rect,
                               (*(*R2).N[depth]).DIR.entries[j].rect,unused)) {
            /* R1 sets the number of dimensions and type */
            (*R1).E[depth]= i; (*R2).E[depth]= j;
            istoread1= (*(*R1).N[depth]).DIR.entries[i].ptrtosub !=
                                                              (*R1).P[depth+1];
            if ( istoread1 ) {
              NewNode(R1,depth+1);
            }
            istoread2= (*(*R2).N[depth]).DIR.entries[j].ptrtosub !=
                                                              (*R2).P[depth+1];
            if ( istoread2 ) {
              NewNode(R2,depth+1);
            }
            JnNv(R1,R2,
                 depth+1,
                 R1rectangle1,R1rectangle2,R2rectangle1,R2rectangle2,
                 Dir1Query,Data1Query,Dir2Query,Data2Query,DirJoin,DataJoin,
                 Manage,buf1,buf2,finish);
          }
        }
      }
    }
  }
  else {
    if (depth == (*R1).parameters._.height) {
    
      c= &(*R1).count;
      if ((*c).countflag) {
        (*c).datavisitcount++;
      }
      
      n= &(*(*R1).N[depth]).DATA;
      
      for (i= 0; i < (*n).nofentries; i++) {
        if (*finish) {return;}
        if (Data1Query(R1,(*n).entries[i].rect,R1rectangle1,R1rectangle2)) {
          JnRgnQueryNv(R2,R1,FALSE,depth,
                       R2rectangle1,R2rectangle2,
                       (*n).entries[i].rect,
                       &(*n).entries[i].info,
                       Data2Query,
                       DirJoin,DataJoin,
                       Manage,buf1,buf2,finish);
        }
      }
    }
    else {
      
      c= &(*R2).count;
      if ((*c).countflag) {
        (*c).datavisitcount++;
      }
      
      n= &(*(*R2).N[depth]).DATA;
        
      for (i= 0; i < (*n).nofentries; i++) {
        if (*finish) {return;}
        if (Data2Query(R2,(*n).entries[i].rect,R2rectangle1,R2rectangle2)) {
          JnRgnQueryNv(R1,R2,TRUE,depth,
                       R1rectangle1,R1rectangle2,
                       (*n).entries[i].rect,
                       &(*n).entries[i].info,
                       Data1Query,
                       DirJoin,DataJoin,
                       Manage,buf1,buf2,finish);
        }
      }
    }
  }
}

/************************************************************************/

static void JnRgnQueryNv(RSTREE R, RSTREE Rx,
                         boolean order,
                         int depth,
                         typrect dataqueryrect1,
                         typrect dataqueryrect2,
                         typrect joinrect,
                         typinfo *joininfo,
                         DataQueryProc DataQuery,
                         DirQueryProc DirJoin,
                         DataQueryProc DataJoin,
                         JoinManageProc Manage,
                         void *buf1,
                         void *buf2,
                         boolean *finish)
{
  refDIRnode DIN;
  refDATAnode DAN;
  refcount c;
  boolean istoread;
  typrect rectR1, rectR2, unused;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
    
    DIN= &(*(*R).N[depth]).DIR;

    if (order) {
      for (i= 0; i < (*DIN).nofentries; i++) {
        if (*finish) {return;}
        if (DirJoin(R,(*DIN).entries[i].rect,joinrect,unused)) {
          (*R).E[depth]= i;
          istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
          if ( istoread ) {
            NewNode(R,depth+1);
          }
          JnRgnQueryNv(R,Rx,order,depth+1,
                       dataqueryrect1,dataqueryrect2,
                       joinrect,
                       joininfo,
                       DataQuery,DirJoin,DataJoin,
                       Manage,buf1,buf2,finish);
        }
      }
    }
    else {
      for (i= 0; i < (*DIN).nofentries; i++) {
        if (*finish) {return;}
        if (DirJoin(Rx,joinrect,(*DIN).entries[i].rect,unused)) {
          (*R).E[depth]= i;
          istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
          if (istoread) {
            NewNode(R,depth+1);
          }
          JnRgnQueryNv(R,Rx,order,depth+1,
                       dataqueryrect1,dataqueryrect2,
                       joinrect,
                       joininfo,
                       DataQuery,DirJoin,DataJoin,
                       Manage,buf1,buf2,finish);
        }
      }
    }
  }
  else {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
    
    DAN= &(*(*R).N[depth]).DATA;
    
    if (order) {
      for (i= 0; i < (*DAN).nofentries; i++) {
        if (*finish) {return;}
        if (DataQuery(R,(*DAN).entries[i].rect,dataqueryrect1,dataqueryrect2)) {
          if (DataJoin(R,(*DAN).entries[i].rect,joinrect,unused)) {
            (*R).E[depth]= i;
            CopyRect(R,(*DAN).entries[i].rect,rectR1);
            CopyRect(R,joinrect,rectR2);
            Manage(R,Rx,
                   rectR1,rectR2,
                   &(*DAN).entries[i].info,joininfo,
                   buf1,buf2,
                   finish);
          }
        }
      }
    }
    else {
      for (i= 0; i < (*DAN).nofentries; i++) {
        if (*finish) {return;}
        if (DataQuery(R,(*DAN).entries[i].rect,dataqueryrect1,dataqueryrect2)) {
          if (DataJoin(Rx,joinrect,(*DAN).entries[i].rect,unused)) {
            (*R).E[depth]= i;
            CopyRect(Rx,joinrect,rectR1);
            CopyRect(Rx,(*DAN).entries[i].rect,rectR2);
            Manage(Rx,R,
                   rectR1,rectR2,
                   joininfo,&(*DAN).entries[i].info,
                   buf1,buf2,
                   finish);
          }
        }
      }
    }
  }
}

/************************************************************************/
