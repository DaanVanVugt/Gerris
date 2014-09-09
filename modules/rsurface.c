#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "RStarTree/RStarTree.h"
#include "rsurface.h"

/* RSurface */

struct _RSurface {
  RSTREE t;
  char * name;
};

RSurface * r_surface_open (const char * fname, const char * mode, int size)
{
  RSurface * rt = malloc (sizeof (RSurface));
  rt->t = NULL;
  if (!OpenRST (&rt->t, fname, "r")) {
    free (rt);
    return NULL;
  }
  rt->name = malloc (sizeof (char)*(strlen (fname) + 1));
  strcpy (rt->name, fname);
  return rt;
}

int r_surface_close (RSurface * rt)
{
  int status = 1;
  status = CloseRST (&rt->t);
  free (rt->name);
  free (rt);
  return status;
}

static boolean Intersects (RSTREE R,
			   typrect RSTrect,
			   typrect queryrect,
			   typrect unused)
{
  int maxdim= NumbOfDim -1;
  boolean inter;
  int d;
  
  d= -1;
  do {
    d++;
    inter= RSTrect[d].l <= queryrect[d].h &&
           RSTrect[d].h >= queryrect[d].l;
  } while (inter && d != maxdim);
  return inter;
}

static void ManageQuery (RSTREE R,
			 typrect rectangle,
			 refinfo infoptr,
			 void ** data,
			 boolean *modify,
			 boolean *finish)
{
  RSurfaceQuery q = data[0];
  double p[3];
  p[0] = rectangle[0].l; p[1] = rectangle[1].l;
  p[2] = infoptr->height;
  (*q) (p, data[1]);
  *modify = FALSE;
  *finish = FALSE;
}

void r_surface_query_region (RSurface * rt, 
			     double min[2], double max[2],
			     RSurfaceQuery q, void * user_data)
{
  typrect rect, unused;
  rect[0].l = min[0]; rect[0].h = max[0];
  rect[1].l = min[1]; rect[1].h = max[1];
  void * data[2];
  data[0] = q;
  data[1] = user_data;
  RegionQuery (rt->t, rect, unused, Intersects, Intersects, (QueryManageProc) ManageQuery, data);
}
