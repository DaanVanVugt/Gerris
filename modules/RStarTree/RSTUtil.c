/* ----- RSTUtil.c ----- */


#include "RStarTree.h"
#include "RSTUtil.h"


/* declarations */

static void ExChange(int *x, int *y);
static void CSWorkAround(RSTREE R,
                         typrect newrect,
                         refDIRnode node,
                         int index,
                         double space,
                         double newspace);


/***********************************************************************/

void FalseArray(int *ptr, int wordsqty)

{
  int i;
  
  i= 1;
  while (i <= wordsqty) {
    *ptr= FALSE;
    ptr++;
    i++;
  }
}

/***********************************************************************/

void CopyRect(RSTREE R,
              typrect from,
              typrect to)

{
  int d;
  for (d= 0; d <= (*R).parameters._.maxdim; d++) {
    to[d]= from[d];
  }
}

/***********************************************************************/

void EvalCenter(RSTREE R, typrect rectangle, typpoint center)

{
  int i;
  
  for (i= 0; i <= (*R).parameters._.maxdim; i++) {
    center[i]= (rectangle[i].l + rectangle[i].h) / 2.0;
  }
}

/***********************************************************************/

double RSTDistance(RSTREE R, typpoint point1, typpoint point2)

{
  double sum, factor;
  int i;

  sum= 0.0;
  for (i= 0; i <= (*R).parameters._.maxdim; i++) {
    factor= point1[i] - point2[i];
    sum= sum+factor*factor;
  }
  return sum; /* relativ distance (sqrt avoided) */
}

/***********************************************************************/

void EvalDirEnclRect(RSTREE R, typDIRnode *node, typrect rectangle)

{
  int maxdim;
  refinterval re;
  int i, j;
  
  maxdim= (*R).parameters._.maxdim;
  CopyRect(R,(*node).entries[0].rect,rectangle);
  for (i= 0; i < (*node).nofentries; i++) {
    for (j= 0; j <= maxdim; j++) {
    
      re= &(*node).entries[i].rect[j];
      
      if (rectangle[j].l > (*re).l) {
        rectangle[j].l= (*re).l;
      }
      if (rectangle[j].h < (*re).h) {
        rectangle[j].h= (*re).h;
      }
    }
  }
}

/***********************************************************************/

void EvalDataEnclRect(RSTREE R, typDATAnode *node, typrect rectangle)

{
  int maxdim;
  refinterval re;
  int i, j;
  
  maxdim= (*R).parameters._.maxdim;
  CopyRect(R,(*node).entries[0].rect,rectangle);
  for (i= 0; i < (*node).nofentries; i++) {
    for (j= 0; j <= maxdim; j++) {
    
      re= &(*node).entries[i].rect[j];
      
      if (rectangle[j].l > (*re).l) {
        rectangle[j].l= (*re).l;
      }
      if (rectangle[j].h < (*re).h) {
        rectangle[j].h= (*re).h;
      }
    }
  }
}

/***********************************************************************/

boolean Overlaps(RSTREE R,
                 typrect rect1,
                 typrect rect2)

{
  boolean ovlp;
  int maxdim;
  int d;
  
  maxdim= (*R).parameters._.maxdim;
  d= -1;
  do {
    d++;
    ovlp= rect1[d].l <= rect2[d].h &&
          rect1[d].h >= rect2[d].l;
  } while (ovlp && d != maxdim);
  return ovlp;
}

/***********************************************************************/

void GetOverlap(RSTREE R,
                typrect r1,
                typrect r2,
                double *spc)

{
  double low, high;
  int i;

  *spc= 1.0;
  for (i= 0; i <= (*R).parameters._.maxdim; i++) {
    if (r1[i].l < r2[i].l) {
      low= r2[i].l;
    }
    else {
      low= r1[i].l;
    }
    if (r1[i].h < r2[i].h) {
      high= r1[i].h;
    }
    else {
      high= r2[i].h;
    }
    *spc= *spc * (high-low);
  }
}

/***********************************************************************/

boolean RSTEqual(RSTREE R,
                 typrect rect1,
                 typrect rect2)

{
  boolean eql;
  int d;
  
  d= -1;
  do {
    d++;
    eql= rect1[d].l == rect2[d].l &&
         rect1[d].h == rect2[d].h;
  } while (eql && d != (*R).parameters._.maxdim);
  return eql;
}

/***********************************************************************/

boolean  Covers(RSTREE R,
                typrect crect,
                typrect rect)

{
  boolean cov;
  int maxdim;
  int d;
  
  maxdim= (*R).parameters._.maxdim;
  d= -1;
  do {
    d++;
    cov= crect[d].l <= rect[d].l &&
         crect[d].h >= rect[d].h;
  } while (cov && d != (*R).parameters._.maxdim);
  return cov;
}

/***********************************************************************/

void QuickSortValArr(int begin,
                     int end,
                     ValueArray value,
                     IndexArray I)
/* Sorts I
   by value[I[i]]. */

{
  double midelem;
  int i, j;
  
  i= begin; j= end;
  midelem= value[I[(i+j) / 2]];
  do {
    while (value[I[i]] < midelem) {
      i++;
    }
    while (value[I[j]]>midelem) {
      j--;
    }
    if (i < j) {
      ExChange(&I[i],&I[j]);
      i++; j--;
    }
    else if (i == j) {
      i++; j--;
    }
  } while (i <= j);
  if (begin < j) {
    if (j - begin > 1) {
      QuickSortValArr(begin,j,value,I);
    }
    else {
      if (value[I[begin]] > value[I[j]]) {
        ExChange(&I[begin],&I[j]);
      }
    }
  }
  if (i < end) {
    if (end - i > 1) {
      QuickSortValArr(i,end,value,I);
    }
    else {
      if (value[I[i]] > value[I[end]]) {
        ExChange(&I[i],&I[end]);
      }
    }
  }
}

/***********************************************************************/

void QuickSortDirEnt(int begin,
                     int end,
                     int dim,
                     Side side,
                     typDIRentries Ntosort,
                     IndexArray I)
/* Sorts I
   primarily   by Ntosort[I[i]].rect[dim].side,
   secondarily by Ntosort[I[i]].rect[dim]."otherside". */

{
  refinterval re;
  typatomkey midlow, midhigh;
  typinterval int1, int2;
  int i, j;
  
  i= begin; j= end;
  if (side == low) {
  
    re = &Ntosort[I[(i+j) / 2]].rect[dim];
    
    midlow= (*re).l;
    midhigh= (*re).h;
    do {
      while ( Ntosort[I[i]].rect[dim].l < midlow ||
              ( Ntosort[I[i]].rect[dim].l == midlow &&
                Ntosort[I[i]].rect[dim].h < midhigh ) ) {
        i++;
      }
      while ( Ntosort[I[j]].rect[dim].l > midlow ||
              ( Ntosort[I[j]].rect[dim].l == midlow &&
                Ntosort[I[j]].rect[dim].h > midhigh ) ) {
        j--;
      }
      if (i < j) {
        ExChange(&I[i],&I[j]);
        i++; j--;
      }
      else if (i == j) {
        i++; j--;
      }
    } while (i <= j);
    if (begin < j) {
      if (j - begin > 1) {
        QuickSortDirEnt(begin,j,dim,low,Ntosort,I);
      }
      else {
        int1= Ntosort[I[begin]].rect[dim];
        int2= Ntosort[I[j]].rect[dim];
        if ( int1.l > int2.l ||
             ( int1.l == int2.l &&
               int1.h > int2.h ) ) {
          ExChange(&I[begin],&I[j]);
        }
      }
    }
    if (i < end) {
      if (end - i > 1) {
        QuickSortDirEnt(i,end,dim,low,Ntosort,I);
      }
      else {
        int1= Ntosort[I[i]].rect[dim];
        int2= Ntosort[I[end]].rect[dim];
        if ( int1.l > int2.l ||
             ( int1.l == int2.l &&
               int1.h > int2.h ) ) {
          ExChange(&I[i],&I[end]);
        }
      }
    }
  }
  else {
  
    re = &Ntosort[I[(i+j) / 2]].rect[dim];
    
    midhigh= (*re).h;
    midlow= (*re).l;
    do {
      while ( Ntosort[I[i]].rect[dim].h < midhigh ||
              ( Ntosort[I[i]].rect[dim].h == midhigh &&
                Ntosort[I[i]].rect[dim].l < midlow ) ) {
        i++;
      };
      while ( Ntosort[I[j]].rect[dim].h > midhigh ||
              ( Ntosort[I[j]].rect[dim].h == midhigh &&
                Ntosort[I[j]].rect[dim].l > midlow ) ) {
        j--;
      };
      if (i < j) {
        ExChange(&I[i],&I[j]);
        i++; j--;
      }
      else if (i == j) {
        i++; j--;
      }
    } while (i <= j);
    if (begin < j) {
      if (j - begin > 1) {
        QuickSortDirEnt(begin,j,dim,high,Ntosort,I);
      }
      else {
        int1= Ntosort[I[begin]].rect[dim];
        int2= Ntosort[I[j]].rect[dim];
        if ( int1.h > int2.h ||
             ( int1.h == int2.h &&
               int1.l > int2.l ) ) {
          ExChange(&I[begin],&I[j]);
        }
      }
    }
    if (i < end) {
      if (end - i > 1) {
        QuickSortDirEnt(i,end,dim,high,Ntosort,I);
      }
      else {
        int1= Ntosort[I[i]].rect[dim];
        int2= Ntosort[I[end]].rect[dim];
        if ( int1.h > int2.h ||
             ( int1.h == int2.h &&
               int1.l > int2.l ) ) {
          ExChange(&I[i],&I[end]);
        }
      }     
    }
  }
}

/***********************************************************************/

void QuickSortDataEnt(int begin,
                      int end,
                      int dim,
                      Side side,
                      typDATAentries Ntosort,
                      IndexArray I)
/* Sorts I
   primarily   by Ntosort[I[i]].rect[dim].side,
   secondarily by Ntosort[I[i]].rect[dim]."otherside". */

{
  refinterval re;
  typatomkey midlow, midhigh;
  typinterval int1, int2;
  int i, j;
  
  i= begin; j= end;
  if (side == low) {
  
    re = &Ntosort[I[(i+j) / 2]].rect[dim];
    
    midlow= (*re).l;
    midhigh= (*re).h;
    do {
      while ( Ntosort[I[i]].rect[dim].l < midlow ||
              ( Ntosort[I[i]].rect[dim].l == midlow &&
                Ntosort[I[i]].rect[dim].h < midhigh ) ) {
        i++;
      }
      while ( Ntosort[I[j]].rect[dim].l > midlow ||
              ( Ntosort[I[j]].rect[dim].l == midlow &&
                Ntosort[I[j]].rect[dim].h > midhigh ) ) {
        j--;
      }
      if (i < j) {
        ExChange(&I[i],&I[j]);
        i++; j--;
      }
      else if (i == j) {
        i++; j--;
      }
    } while (i <= j);
    if (begin < j) {
      if (j - begin > 1) {
        QuickSortDataEnt(begin,j,dim,low,Ntosort,I);
      }
      else {
        int1= Ntosort[I[begin]].rect[dim];
        int2= Ntosort[I[j]].rect[dim];
        if ( int1.l > int2.l ||
             ( int1.l == int2.l &&
               int1.h > int2.h ) ) {
          ExChange(&I[begin],&I[j]);
        }
      }
    }
    if (i < end) {
      if (end - i > 1) {
        QuickSortDataEnt(i,end,dim,low,Ntosort,I);
      }
      else {
        int1= Ntosort[I[i]].rect[dim];
        int2= Ntosort[I[end]].rect[dim];
        if ( int1.l > int2.l ||
             ( int1.l == int2.l &&
               int1.h > int2.h ) ) {
          ExChange(&I[i],&I[end]);
        }
      }
    }
  }
  else {
  
    re = &Ntosort[I[(i+j) / 2]].rect[dim];
    
    midhigh= (*re).h;
    midlow= (*re).l;
    do {
      while ( Ntosort[I[i]].rect[dim].h < midhigh ||
              ( Ntosort[I[i]].rect[dim].h == midhigh &&
                Ntosort[I[i]].rect[dim].l < midlow ) ) {
        i++;
      };
      while ( Ntosort[I[j]].rect[dim].h > midhigh ||
              ( Ntosort[I[j]].rect[dim].h == midhigh &&
                Ntosort[I[j]].rect[dim].l > midlow ) ) {
        j--;
      };
      if (i < j) {
        ExChange(&I[i],&I[j]);
        i++; j--;
      }
      else if (i == j) {
        i++; j--;
      }
    } while (i <= j);
    if (begin < j) {
      if (j - begin > 1) {
        QuickSortDataEnt(begin,j,dim,high,Ntosort,I);
      }
      else {
        int1= Ntosort[I[begin]].rect[dim];
        int2= Ntosort[I[j]].rect[dim];
        if ( int1.h > int2.h ||
             ( int1.h == int2.h &&
               int1.l > int2.l ) ) {
          ExChange(&I[begin],&I[j]);
        }
      }
    }
    if (i < end) {
      if (end - i > 1) {
        QuickSortDataEnt(i,end,dim,high,Ntosort,I);
      }
      else {
        int1= Ntosort[I[i]].rect[dim];
        int2= Ntosort[I[end]].rect[dim];
        if ( int1.h > int2.h ||
             ( int1.h == int2.h &&
               int1.l > int2.l ) ) {
          ExChange(&I[i],&I[end]);
        }
      }     
    }
  }
}

/***********************************************************************/

static void ExChange(int *x, int *y)

{
  int z;
  
  z= *x; *x= *y; *y= z;
}

/***********************************************************************/

void ChooseSubtree(RSTREE R,
                   typrect newrect,
                   int depth,
                   typDIRnode *node,
                   int *found)

{
  int maxexam, inv;
  typatomkey low, high;
  double space, newspace,
  enlarge,
  validspace, ovlpspc, overlap, leastovlp;
  boolean didfit, first, ok;
  ValueArray enlarr;
  IndexArray I;
  typrect enlargedrect;
  int maxdim;
  refparameters par;
  refinterval re;
  refDIRent e;
  int i, j, k;
  
  didfit= FALSE; ok= TRUE;
  maxdim= (*R).parameters._.maxdim;
  for (i= 0; i < (*node).nofentries; i++) {
    space= 1.0; newspace= 1.0;
    for (j= 0; j <= maxdim; j++) {
    
      re= &(*node).entries[i].rect[j];
      
      low= (*re).l;
      high= (*re).h;
      space= space * (high - low);
      
      re= &newrect[j];
      
      if (low > (*re).l) {
        low= (*re).l;
      }
      if (high < (*re).h) {
        high= (*re).h;
      }
      newspace= newspace * (high - low);
    }
    if (space == 0.0) {
      CSWorkAround(R,newrect,node,i,space,newspace);
      ok= FALSE;
    } /* trap degenerated rectangles */
    if (didfit) {
      if (newspace == space) { /* it fits */
        if (space < validspace) {
          validspace= space;
          *found= i;
        }
      }
    }
    else {
      enlarge= newspace - space;
      if (enlarge == 0.0) {
        validspace= space;
        *found= i;
        didfit= TRUE;
      }
      else {
        enlarr[i]= enlarge;
      }
    }
  }
  if (! didfit) {
    for (i= 0; i < (*node).nofentries; i++) {
      I[i]= i;
    }
    QuickSortValArr(0,(*node).nofentries - 1,enlarr,I);
    if (ok && depth == (*R).parameters._.height - 1) { /* Subtrees are leafs */
    
      par = &(*R).parameters._;
      
      if ((*node).nofentries < (*par).nbrsexam) {
        maxexam= (*node).nofentries;
      }
      else {
        maxexam= (*par).nbrsexam;
      }

/*    maxexam= (*node).nofentries; *//* <- all  *) */

      first= TRUE;
      for (i= 0; i < maxexam; i++) {
        for (j= 0; j <= maxdim; j++) {
        
          re= &(*node).entries[I[i]].rect[j];
          
          low= (*re).l;
          high= (*re).h;
          
          re= &enlargedrect[j];
          
          (*re).l= low;
          (*re).h= high;
          if (low > newrect[j].l) {
            (*re).l= newrect[j].l;
          }
          if (high < newrect[j].h) {
            (*re).h= newrect[j].h;
          }
        } /* Create enlargedrect */
        overlap= 0.0;
        for (k= 0; k < (*node).nofentries; k++) {
          if (k != i) {
          
            e= &(*node).entries[I[k]];
            
            if (Overlaps(R,enlargedrect,(*e).rect)) {
              GetOverlap(R,enlargedrect,(*e).rect,&ovlpspc);
              overlap= overlap+ovlpspc;
              if (Overlaps(R,(*node).entries[I[i]].rect,(*e).rect)) {
                GetOverlap(R,(*node).entries[I[i]].rect,(*e).rect,&ovlpspc);
                overlap= overlap-ovlpspc;
              }
            }
          }
        }
        if (first) {
          leastovlp= overlap;
          *found= I[0];
          inv= 1;
          first= FALSE;
        }
        else {
          if (overlap < leastovlp) {
            leastovlp= overlap;
            *found= I[i];
            inv= i;
          }
        }
      }
    }
    else { /* Subtrees are not leafs */
      *found= I[0];
    }
  }
}

/***********************************************************************/

static void CSWorkAround(RSTREE R,
                         typrect newrect,
                         typDIRnode *node,
                         int index,
                         double space,
                         double newspace)

{
#define epsilon (1.0e-100)
  
  typatomkey low, high;
  double distance;
  refinterval re;
  int j;
  
  space= 1.0; newspace= 1.0;
  for (j= 0; j <= (*R).parameters._.maxdim; j++) {
  
    re= &(*node).entries[index].rect[j];
    
    low= (*re).l;
    high= (*re).h;
    distance= high - low;
    if (distance == 0.0) {
      distance= epsilon;
    }
    space= space * distance;
    
    re= &newrect[j];
    
    if (low > (*re).l) {
      low= (*re).l;
    }
    if (high < (*re).h) {
      high= (*re).h;
    }
    distance= high - low;
    if (distance == 0.0) {
      distance= epsilon;
    }
    newspace= newspace * distance;
  }

#undef epsilon
}

/***********************************************************************/

void AdjustChain(RSTREE R,
                 int depth,
                 typrect newrect)

{
  refinterval re;
  refcount c;
  int maxdim;
  boolean adjusting;
  int j;
  
  c= &(*R).count;
  maxdim= (*R).parameters._.maxdim;

  adjusting= TRUE;
  
  if (depth == (*R).parameters._.height) {
    adjusting= FALSE;
    for (j= 0; j <= maxdim; j++) {
      
      re = &(*(*R).N[depth]).DATA.entries[(*R).E[depth]].rect[j];
        
      if ((*re).l > newrect[j].l) {
        (*re).l= newrect[j].l;
        adjusting= TRUE;
      }
      if ((*re).h < newrect[j].h) {
        (*re).h= newrect[j].h;
        adjusting= TRUE;
      }
    }
    if (adjusting) {
      (*R).Nmodified[depth]= TRUE;
      if ((*c).countflag) {
        (*c).dirmodifycount++;
      }
    }
    depth--;
  }
  while (depth != 0 && adjusting) {
    adjusting= FALSE;
    for (j= 0; j <= maxdim; j++) {
      
      re = &(*(*R).N[depth]).DIR.entries[(*R).E[depth]].rect[j];
        
      if ((*re).l > newrect[j].l) {
        (*re).l= newrect[j].l;
        adjusting= TRUE;
      }
      if ((*re).h < newrect[j].h) {
        (*re).h= newrect[j].h;
        adjusting= TRUE;
      }
    }
    if (adjusting) {
      (*R).Nmodified[depth]= TRUE;
      if ((*c).countflag) {
        (*c).dirmodifycount++;
      }
    }
    depth--;
  }
}

/***********************************************************************/

void AdjustChainAfterDeletion(RSTREE R,
                              int depth)

{
  refDIRent e;
  refcount c;
  boolean changed;
  typrect helprect;
  
  c= &(*R).count;
  
  changed= TRUE;
  
  if (depth == (*R).parameters._.height && depth != 1) {
    EvalDataEnclRect(R,&(*(*R).N[depth]).DATA,helprect);
    depth--;
      
    e = &(*(*R).N[depth]).DIR.entries[(*R).E[depth]];
      
    changed= ! RSTEqual(R,helprect,(*e).rect);
    if (changed) {
      CopyRect(R,helprect,(*e).rect);
      (*R).Nmodified[depth]= TRUE;
      if ((*c).countflag) {
        (*c).dirmodifycount++;
      }
    }
  }
  while (depth != 1 && changed) {
    EvalDirEnclRect(R,&(*(*R).N[depth]).DIR,helprect);
    depth--;
      
    e = &(*(*R).N[depth]).DIR.entries[(*R).E[depth]];
      
    changed= ! RSTEqual(R,helprect,(*e).rect);
    if (changed) {
      CopyRect(R,helprect,(*e).rect);
      (*R).Nmodified[depth]= TRUE;
      if ((*c).countflag) {
        (*c).dirmodifycount++;
      }
    }
  }
}

/***********************************************************************/
