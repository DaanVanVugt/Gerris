/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "init.h"
#include "simulation.h"
#include "graphic.h"
#include "solid.h"
#include "adaptive.h"

#define DEBUG 0

static void merged_draw (GSList * merged, FILE * fp)
{
  if (merged->next != NULL) {
    GSList * i = merged;

    while (i) {
      FttCell * cell = i->data;
      FttCellNeighbors n;
      FttCellFace f;  

      f.cell = cell;
      ftt_cell_neighbors (cell, &n);
      for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++)
	if (!n.c[f.d] || !g_slist_find (merged, n.c[f.d]))
	  ftt_face_draw (&f, fp);
      i = i->next;
    }
  }
}

static gdouble local_size_ratio (GtsSegment * s, GfsDomain * domain)
{
  GtsPoint * p1 = GTS_POINT (s->v1);
  GtsPoint * p2 = GTS_POINT (s->v2);
  gdouble l = gts_point_distance (p1, p2);
  FttVector p;
  FttCell * cell;
  gdouble size = G_MAXDOUBLE;

  p.x = p1->x;
  p.y = p1->y;
  p.z = p1->z;
  cell = gfs_domain_locate (domain, p, -1, NULL);
  if (cell)
    size = ftt_cell_size (cell);

  p.x = p2->x;  
  p.y = p2->y;
  p.z = p2->z;
  cell = gfs_domain_locate (domain, p, -1, NULL);
  if (cell) {
    gdouble s = ftt_cell_size (cell);
    
    if (size == G_MAXDOUBLE || s > size)
      size = s;
  }
  
  return size/l;
}

static gboolean stop (gdouble cost, guint nedge)
{
  if (cost >= 1. || nedge > 50000)
    return TRUE;
  return FALSE;
}

static void draw_vector (FttCell * cell, gpointer * data)
{
  gdouble * scale = data[0];
  GfsVariable ** u = data[1];
  FILE * fp = stdout;
  FttVector pos, f;

  gfs_cell_cm (cell, &pos);
  
  f.x = GFS_VALUE (cell, u[0])*(*scale);
  f.y = GFS_VALUE (cell, u[1])*(*scale);
#if FTT_2D
  f.z = 0.;
#else
  f.z = GFS_VALUE (cell, u[2])*(*scale);
#endif
  fprintf (fp, "VECT 1 3 0 3 0 %g %g %g %g %g %g %g %g %g\n",
	   pos.x + f.x - (f.x - f.y/2.)/5.,
	   pos.y + f.y - (f.x/2. + f.y)/5.,
	   pos.z + f.z,
	   pos.x + f.x,
	   pos.y + f.y,
	   pos.z + f.z,
	   pos.x + f.x - (f.x + f.y/2.)/5.,
	   pos.y + f.y + (f.x/2. - f.y)/5.,
	   pos.z + f.z);
  fprintf (fp, "VECT 1 2 0 2 0 %g %g %g %g %g %g\n",
	   pos.x, pos.y, pos.z,
	   pos.x + f.x,
	   pos.y + f.y,
	   pos.z + f.z);
}

static void compute_mixed_vorticity (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * u = data[1];
  FttVector g;

  g_assert (((cell)->flags & GFS_FLAG_DIRICHLET) != 0);
  gfs_cell_dirichlet_gradient (cell, u->i, -1, GFS_STATE (cell)->solid->fv, &g);
  if (u->component == FTT_X)
    GFS_VALUE (cell, v) -= g.y;
  else
    GFS_VALUE (cell, v) += g.x;
}

static void output_mixed_vorticity (FttCell * cell, GfsVariable * v)
{
  gdouble size = ftt_cell_size (cell);
  GfsSolidVector * s = GFS_STATE (cell)->solid;

  printf ("%g %g %g %g\n", s->ca.x, s->ca.y, s->ca.z, 
	  GFS_VALUE (cell, v)/size);
}

static void output_mixed_pressure (FttCell * cell, GfsVariable * p)
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;

  printf ("%g %g %g %g\n", s->ca.x, s->ca.y, s->ca.z, 
	  gfs_dimensional_value (p, gfs_interpolate (cell, s->ca, p)));
}

static void output_mixed_variable (FttCell * cell, GfsVariable * v)
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;

  printf ("%g %g %g %g\n", s->ca.x, s->ca.y, s->ca.z,
	  gfs_dimensional_value (v, GFS_VALUE (cell, v)));
}

/* SVertex: Header */

typedef struct _SVertex         SVertex;

struct _SVertex {
  /*< private >*/
  GtsVertex parent;

  /*< public >*/
  gdouble s;
};

#define S_VERTEX(obj)            GTS_OBJECT_CAST (obj,\
					         SVertex,\
					         s_vertex_class ())
#define IS_S_VERTEX(obj)         (gts_object_is_from_class (obj,\
						 s_vertex_class ()))

/* SVertex: Object */

static GtsVertexClass * s_vertex_class (void)
{
  static GtsVertexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo s_vertex_info = {
      "SVertex",
      sizeof (SVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_vertex_class ()),
				  &s_vertex_info);
  }

  return klass;
}

typedef struct {
  GSList *** p;
  guint nx, ny;
  gdouble h;
  FttVector min, max;
} ClosestGrid;

static ClosestGrid * closest_box_grid_new (FttVector min, FttVector max,
					   gdouble h)
{
  ClosestGrid * g = g_malloc (sizeof (ClosestGrid));
  guint i;

  g->max = max;
  g->min = min;
  g->nx = (g->max.x - g->min.x)/h + 3;
  g->h = h;
  g->ny = (g->max.y - g->min.y)/h + 3;
  g->min.x -= h;
  g->min.y -= h;
  g->max.x = g->min.x + g->nx*g->h;
  g->max.y = g->min.y + g->ny*g->h;
  g->p = g_malloc (g->nx*sizeof (GSList **));
  for (i = 0; i < g->nx; i++)
    g->p[i] = g_malloc0 (g->ny*sizeof (GSList *));

  return g;
}

static void min_max_extent (FttCell * cell, ClosestGrid * g)
{
  FttVector pos;
  
  ftt_cell_pos (cell, &pos);
  if (pos.x > g->max.x) g->max.x = pos.x;
  if (pos.y > g->max.y) g->max.y = pos.y;
  if (pos.x < g->min.x) g->min.x = pos.x;
  if (pos.y < g->min.y) g->min.y = pos.y;
}

static ClosestGrid * closest_grid_new (GfsDomain * domain, gdouble h)
{
  ClosestGrid * g = g_malloc (sizeof (ClosestGrid));
  guint i;

  g->max.x = - G_MAXDOUBLE;
  g->max.y = - G_MAXDOUBLE;
  g->max.z = - G_MAXDOUBLE;
  g->min.x = G_MAXDOUBLE;
  g->min.y = G_MAXDOUBLE;
  g->min.z = G_MAXDOUBLE;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) min_max_extent, g);
  g->nx = (g->max.x - g->min.x)/h + 3;
  g->h = h;
  g->ny = (g->max.y - g->min.y)/h + 3;
  g->min.x -= h;
  g->min.y -= h;
  g->max.x = g->min.x + g->nx*g->h;
  g->max.y = g->min.y + g->ny*g->h;
  g->p = g_malloc (g->nx*sizeof (GSList **));
  for (i = 0; i < g->nx; i++)
    g->p[i] = g_malloc0 (g->ny*sizeof (GSList *));

  return g;
}

static void closest_grid_destroy (ClosestGrid * g)
{
  guint i, j;

  for (i = 0; i < g->nx; i++) {
    for (j = 0; j < g->ny; j++)
      g_slist_free (g->p[i][j]);
    g_free (g->p[i]);
  }
  g_free (g->p);
  g_free (g);
}

typedef enum {
  INSERTED, ALREADY_THERE, CLOSED, SELFCLOSED, FINISHED, OFFSIDE
} InsertStatus;

static gdouble point_distance2 (GtsPoint * p1, GtsPoint * p2)
{
  return ((p1->x - p2->x)*(p1->x - p2->x) + (p1->y - p2->y)*(p1->y - p2->y));
}

static InsertStatus closest_grid_is_insertable (ClosestGrid * g, GtsPoint * p,
						gdouble dmin, 
						gdouble ds)
{
  gint i = (p->x - g->min.x)/(g->max.x - g->min.x)*g->nx, i1;
  gint j = (p->y - g->min.y)/(g->max.y - g->min.y)*g->ny, j1;

  if (i < 0 || i >= g->nx || j < 0 || j >= g->ny)
    return OFFSIDE;

  for (i1 = i - 1; i1 <= i + 1; i1++)
    for (j1 = j - 1; j1 <= j + 1; j1++) 
      if (i1 >= 0 && i1 < g->nx && j1 >= 0 && j1 < g->ny) {
	GSList * k = g->p[i1][j1];
	while (k) {
	  GtsPoint * p1 = k->data;
	  
	  if (IS_S_VERTEX (p) && IS_S_VERTEX (p1) && S_VERTEX (p1)->s >= 0.) {
	    gdouble d = point_distance2 (p, p1);
	    
	    if (d == 0.) return ALREADY_THERE;
	    if (S_VERTEX (p)->s - S_VERTEX (p1)->s >= 3.*ds && d < 0.9*ds*ds)
	      return S_VERTEX (p1)->s == 0. ? CLOSED : SELFCLOSED;
	  }
	  else if ((!IS_S_VERTEX (p1) || S_VERTEX (p1)->s < 0.) &&
		   point_distance2 (p, p1) < dmin*dmin)
	    return FINISHED;
	  k = k->next;
	}
      }
  return INSERTED;
}

static InsertStatus closest_grid_add (ClosestGrid * g, GtsPoint * p,
				      gdouble dmin, gdouble ds)
{
  InsertStatus status;
  gint i = (p->x - g->min.x)/(g->max.x - g->min.x)*g->nx;
  gint j = (p->y - g->min.y)/(g->max.y - g->min.y)*g->ny;

  if ((status = closest_grid_is_insertable (g, p, dmin, ds)) != INSERTED)
    return status;
  g->p[i][j] = g_slist_prepend (g->p[i][j], p);
  return INSERTED;
}

static void closest_grid_remove (GtsPoint * p, ClosestGrid * g)
{
  gint i = (p->x - g->min.x)/(g->max.x - g->min.x)*g->nx;
  gint j = (p->y - g->min.y)/(g->max.y - g->min.y)*g->ny;

  if (i < 0 || i >= g->nx || j < 0 || j >= g->ny)
    return;
  g->p[i][j] = g_slist_remove (g->p[i][j], p);
}

static InsertStatus insert (FttVector p, ClosestGrid * g,
			    gdouble dmin, gdouble s, gdouble ds,
			    GtsVertex ** vertex)
{
  GtsVertex * v = gts_vertex_new (s_vertex_class (), p.x, p.y, p.z);
  InsertStatus status;

  S_VERTEX (v)->s = s;
  if ((status = closest_grid_add (g, GTS_POINT (v), dmin, ds)) != INSERTED) {
    gts_object_destroy (GTS_OBJECT (v));
    *vertex = NULL;
    return status;
  }
  *vertex = v;
  return INSERTED;
}

static gboolean advect (GfsDomain * domain,
			FttCell * cell,
			FttVector * p,
			gdouble ds,
			gint direction)
{
  FttComponent c;
  FttVector u, ph;
  gdouble nu = 0.;
  guint n = 10;
  gdouble h = ds/n;
  gboolean ad = TRUE;
  GfsVariable ** U = gfs_domain_velocity (domain);

  while (n-- > 0 && ad) {
    for (c = 0; c < 2/*FTT_DIMENSION*/; c++) {
      ((gdouble *) &u)[c] = direction*gfs_interpolate (cell, *p, U[c]);
      nu += ((gdouble *) &u)[c]*((gdouble *) &u)[c];
    }
    if (nu > 0.) {
      nu = sqrt (nu);
      ph = *p;
      for (c = 0; c < 2/*FTT_DIMENSION*/; c++)
	((gdouble *) &ph)[c] += h*((gdouble *) &u)[c]/(2.*nu);
      cell = gfs_domain_locate (domain, ph, -1, NULL);
      if (cell != NULL) {
	nu = 0.;
	for (c = 0; c < 2/*FTT_DIMENSION*/; c++) {
	  ((gdouble *) &u)[c] = direction*gfs_interpolate (cell, ph, U[c]);
	  nu += ((gdouble *) &u)[c]*((gdouble *) &u)[c];
	}
	if (nu > 0.) {
	  nu = sqrt (nu);
	  for (c = 0; c < 2/*FTT_DIMENSION*/; c++)
	    ((gdouble *) p)[c] += h*((gdouble *) &u)[c]/nu;
	}
	else
	  ad = FALSE;
      }
      else
	ad = FALSE;
    }
    else
      ad = FALSE;
  }
  return ad;
}

static InsertStatus grow_streamline (GfsDomain * domain,
				     ClosestGrid * grid,
				     FttVector p,
				     gdouble dmin,
				     gdouble rds,
				     gint direction,
				     GSList ** stream)
{
  FttCell * cell = gfs_domain_locate (domain, p, -1, NULL);
  GtsVertex * v, * vstart = NULL;
  gdouble s = 0.;
  InsertStatus status = cell ? INSERTED : OFFSIDE;

  while (status == INSERTED || status == ALREADY_THERE) {
    gdouble ds = rds*ftt_cell_size (cell);

    ds = MIN (dmin, ds);
    switch ((status = insert (p, grid, dmin, s, ds, &v))) {
    case INSERTED:
      *stream = g_slist_prepend (*stream, v);
#if DEBUG
      fprintf (stderr, "%g %g\n", GTS_POINT (v)->x, GTS_POINT (v)->y);
      fflush (stderr);
#endif
      if (!vstart) vstart = v;
    case ALREADY_THERE:
      if (advect (domain, cell, &p, ds, direction)) {
	s += ds;
	cell = gfs_domain_locate (domain, p, -1, NULL);
	if (cell == NULL)
	  status = OFFSIDE;
      }
      else
	status = OFFSIDE;
      break;
    default:
      ;
    }
  }
  if (direction > 0) {
    if (status == CLOSED)
      *stream = g_slist_prepend (*stream, vstart);
    *stream = g_slist_reverse (*stream);
  }
  return status;
}

static void set_not_current (SVertex * s)
{
  g_assert (IS_S_VERTEX (s));
  s->s = -1.;
}

static void streamline_destroy (GSList * s, ClosestGrid * grid)
{
  g_slist_foreach (s, (GFunc) closest_grid_remove, grid);
  g_slist_foreach (s, (GFunc) gts_object_destroy, NULL);
  g_slist_free (s);
}

static GSList * streamline (GfsDomain * domain,
			    ClosestGrid * grid,
			    FttVector p,
			    gdouble dmin,
			    gdouble rds,
			    gboolean closed)
{
  GSList * stream = NULL;
  InsertStatus status = grow_streamline (domain, grid, p, dmin, rds, 1, 
					 &stream);

  if (!closed) {
    if (status != CLOSED && status != SELFCLOSED)
      grow_streamline (domain, grid, p, dmin, rds, -1, &stream);
  }
  else {
    if (status != CLOSED && status != SELFCLOSED && status != OFFSIDE) {
      streamline_destroy (stream, grid);
      return NULL;
    }
    if (status == OFFSIDE) {
      status = grow_streamline (domain, grid, p, dmin, rds, -1, &stream);
      if (status != CLOSED && status != SELFCLOSED && status != OFFSIDE) {
	streamline_destroy (stream, grid);
	return NULL;
      }
    }
  }
  g_slist_foreach (stream, (GFunc) set_not_current, NULL);
#if DEBUG
  fprintf (stderr, "\n"); fflush (stderr);
#endif
  return stream;
}

static gboolean seed (GSList * i, 
		      GfsDomain * domain,
		      ClosestGrid * grid,
		      gdouble dsep,
		      gdouble dmin,
		      GList ** streams,
		      gboolean closed)
{
  GtsPoint * v = gts_point_new (gts_point_class (), 0., 0., 0.);
  FttVector p;

  p.z = 0.;/*-0.49;*/
  while (i) {
    GtsPoint * p1 = i->data;
    i = i->next;
    if (i) {
      GtsPoint * p2 = i->data;
      gdouble d = sqrt (point_distance2 (p1, p2));

      if (d > 1e-6) {
	v->x = p.x = (p1->x + p2->x)/2. - (p2->y - p1->y)*dsep/d;
	v->y = p.y = (p1->y + p2->y)/2. + (p2->x - p1->x)*dsep/d;
	if (gfs_domain_locate (domain, p, -1, NULL) &&
	    closest_grid_is_insertable (grid, v, dsep, 0.) == INSERTED) {
	  GSList * s = streamline (domain, grid, p, dmin, 0.25, closed);

	  if (s) {
	    *streams = g_list_prepend (*streams, s);
	    gts_object_destroy (GTS_OBJECT (v));
	    return TRUE;
	  }
	}
	v->x = p.x = (p1->x + p2->x)/2. + (p2->y - p1->y)*dsep/d;
	v->y = p.y = (p1->y + p2->y)/2. - (p2->x - p1->x)*dsep/d;
	if (gfs_domain_locate (domain, p, -1, NULL) &&
	    closest_grid_is_insertable (grid, v, dsep, 0.) == INSERTED) {
	  GSList * s = streamline (domain, grid, p, dmin, 0.25, closed);

	  if (s) {
	    *streams = g_list_prepend (*streams, s);
	    gts_object_destroy (GTS_OBJECT (v));
	    return TRUE;
	  }
	}
      }
      i = i->next;
    }
  }
  gts_object_destroy (GTS_OBJECT (v));
  return FALSE;
}

static void cell_center (FttCell * cell, gpointer * data)
{
  FttVector * p = data[0], pos;

  if (p->x == G_MAXDOUBLE) {
    GfsDomain * domain = data[1];
    ClosestGrid * grid = data[2];
    GtsPoint * v;

    ftt_cell_pos (cell, &pos);
    pos.z = 0.;
    v = gts_point_new (gts_point_class (), pos.x, pos.y, pos.z);
    if (gfs_domain_locate (domain, pos, -1, NULL) && 
	closest_grid_is_insertable (grid, v, 0., 0.) == INSERTED)
      *p = pos;
    gts_object_destroy (GTS_OBJECT (v));
  }
}

static GList * even_streamlines (GfsDomain * domain,
				 ClosestGrid * grid,
				 gdouble dsep, 
				 gdouble dmin,
				 gboolean closed)
{
  GList * streams = NULL, * current;
  FttVector p = {G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE};
  gboolean finished = FALSE;
  gpointer data[3];

  data[0] = &p;
  data[1] = domain;
  data[2] = grid;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) cell_center, data);
  current = streams = g_list_prepend (streams, 
			   streamline (domain, grid, p, dmin, 0.25, closed));
  do {
    if (!seed (current->data, domain, grid, dsep, dmin, &streams, closed)) {
      if (current == streams)
	finished = TRUE;
      else
	current = current->prev;
    }
  } while (!finished);
  return streams;
}

static gdouble curve_cost (GtsPoint *  p1, GtsPoint *  p2, GtsPoint *  p3)
{
  GtsVector v1, v2, a;

  v1[0] = p2->x - p1->x; v1[1] = p2->y - p1->y; v1[2] = p2->z - p1->z;
  v2[0] = p3->x - p2->x; v2[1] = p3->y - p2->y; v2[2] = p3->z - p2->z;
  gts_vector_cross (a, v1, v2);
  return gts_vector_norm (a)/2.;
}

static GSList * simplify_stream (GSList * stream,
				 gdouble maxcost)
{
  GSList * i = stream;
  GSList * s = NULL;
  GtsPoint * p1 = NULL, * p2 = NULL;
  gdouble cost = 0.;

  while (i) {
    GtsPoint * p = i->data;

    if (p1 == NULL) { 
      p1 = p;
      s = g_slist_prepend (s, p);
    }
    else if (p2 == NULL)
      p2 = p;
    else {
      cost += curve_cost (p1, p2, p);
      p1 = p2;
      p2 = p;
      if (cost > maxcost || !i->next) {
	s = g_slist_prepend (s, p);
	cost = 0.;
      }
      else
	GTS_OBJECT (p)->reserved = p;
    }
    i = i->next;
  }
  i = stream;
  while (i) {
    if (GTS_OBJECT (i->data)->reserved == i->data)
      gts_object_destroy (i->data);
    i = i->next;
  }
  g_slist_free (stream);
  return s;
}

static void write_stream (GSList * i, FILE * fp)
{
  guint n = g_slist_length (i);

  fprintf (fp, "VECT 1 %u 0 %u 0\n", n, n);
  while (i) {
    GtsPoint * p = i->data;
    fprintf (fp, "%g %g 0\n", p->x, p->y);
    i = i->next;
  }
}

static void update_var (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsFunction * f = data[1];

  GFS_VALUE (cell, v) = gfs_function_value (f, cell);
}

static void velocity_norm (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable ** u = data[1];
  GFS_VALUE (cell, v) = gfs_vector_norm (cell, u);
}

int main (int argc, char * argv[])
{
  int c = 0;
  GfsVariable * var = NULL;
  GtsFile * fp;
  GtsSurface * surface = NULL;
  gboolean draw_surface = FALSE;

  gboolean verbose = FALSE;
  gboolean refine = FALSE;

  gdouble vector = -1.;

  FILE * stream = NULL;
  gchar * streamname = NULL, * color = NULL;
  gboolean ribbon = FALSE, lines = FALSE, squares = FALSE;

  GtsBBox * box = NULL;

  gdouble min = 0., max = 0.;
  gboolean gnuplot = FALSE;

  gboolean merged = FALSE;
  gboolean reinit = FALSE;
  gboolean mixed = FALSE;
  gdouble even_stream = 0., rdmin = 0.5, maxcost = 2e-7;
  FttVector bmin, bmax = { -G_MAXDOUBLE, -G_MAXDOUBLE, -G_MAXDOUBLE };
  gboolean closed = FALSE;
  gint level = -1;
  gdouble iso = G_MAXDOUBLE;

  FILE * profile = NULL;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"iso", required_argument, NULL, 'u'},
      {"level", required_argument, NULL, 'L'},
      {"closed", no_argument, NULL, 'j'},
      {"stream", required_argument, NULL, 'f'},
      {"dmin", required_argument, NULL, 'I'},
      {"maxcost", required_argument, NULL, 'O'},
      {"box", required_argument, NULL, 'b'},
      {"profile", required_argument, NULL, 'p'},
      {"mixed", no_argument, NULL, 'o'},
      {"reinit", no_argument, NULL, 'i'},
      {"merged", no_argument, NULL, 'e'},
      {"min", required_argument, NULL, 'm'},
      {"max", required_argument, NULL, 'M'},
      {"squares", no_argument, NULL, 'S'},
      {"gnuplot", no_argument, NULL, 'g'},
      {"sx", required_argument, NULL, 'x'},
      {"sy", required_argument, NULL, 'y'},
      {"sz", required_argument, NULL, 'z'},
      {"color", required_argument, NULL, 'c'},
      {"streamlines", required_argument, NULL, 'l'},
      {"cylinder", required_argument, NULL, 'C'},
      {"ribbon", required_argument, NULL, 'R'},
      {"refine", no_argument, NULL, 'r'},
      {"surface", required_argument, NULL, 's'},
      {"vector", required_argument, NULL, 'V'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, 
			      "hvs:rV:C:R:c:x:y:z:Sm:M:eiop:f:I:O:b:jl:L:u:g",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, 
			 "hvs:rV:C:R:c:x:y:z:Sm:M:geiop:f:I:O:b:jl:L:u:g"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'u': /* isosurface */
      iso = atof (optarg);
      break;
    case 'L': /* level */
      level = atoi (optarg);
      break;
    case 'j': /* closed */
      closed = TRUE;
      rdmin = 0.9;
      break;
    case 'b': { /* box */
      gchar * s = strtok (optarg, ",");
      guint i = 0;

      while (i < 3 && s != NULL) {
	(&bmin.x)[i++] = atof (s);
	s = strtok (NULL, ",");
      }
      if (i != 3) {
	fprintf (stderr, "gfs2oogl: expecting a number for option `--box'\n");
	fprintf (stderr, "Try `gfs2oogl --help' for more information.\n");
	return 1;
      }
      i = 0;
      while (i < 3 && s != NULL) {
	(&bmax.x)[i++] = atof (s);
	s = strtok (NULL, ",");
      }
      if (i != 3) {
	fprintf (stderr, "gfs2oogl: expecting a number for option `--box'\n");
	fprintf (stderr, "Try `gfs2oogl --help' for more information.\n");
	return 1;
      }
      break;
    }
    case 'f': /* stream */
      even_stream = atof (optarg);
      break;
    case 'I': /* dmin */
      rdmin = atof (optarg);
      break;
    case 'O': /* maxcost */
      maxcost = atof (optarg);
      break;
    case 'p': /* profile */
      if ((profile = fopen (optarg, "rt")) == NULL) {
	fprintf (stderr, "gfs2oogl: cannot open file `%s'\n"
		 "Try `gfs2oogl --help' for more information.\n", optarg);
	return 1; /* failure */
      }
      break;
   case 'o': /* mixed */
      mixed = TRUE;
      break;
    case 'i': /* reinit */
      reinit = TRUE;
      break;
    case 'e': /* merged */
      merged = TRUE;
      break;
    case 'M': /* max */
      max = atof (optarg);
      break;
    case 'm': /* min */
      min = atof (optarg);
      break;
    case 'g': /* gnuplot */
      gnuplot = TRUE;
      break;
    case 'S': /* squares */
      squares = TRUE;
      break;
    case 'x': /* sx */
      box = gts_bbox_new (gts_bbox_class (), NULL,
			  atof (optarg), -G_MAXDOUBLE/2., -G_MAXDOUBLE/2.,
			  atof (optarg), G_MAXDOUBLE/2., G_MAXDOUBLE/2.);
      break;
    case 'y': /* sy */
      box = gts_bbox_new (gts_bbox_class (), NULL,
			  -G_MAXDOUBLE/2., atof (optarg), -G_MAXDOUBLE/2.,
			  G_MAXDOUBLE/2., atof (optarg), G_MAXDOUBLE/2.);
      break;
    case 'z': /* sz */
      box = gts_bbox_new (gts_bbox_class (), NULL,
			  -G_MAXDOUBLE/2., -G_MAXDOUBLE/2., atof (optarg),
			  G_MAXDOUBLE/2., G_MAXDOUBLE/2., atof (optarg));
      break;
    case 's': /* surface */
      draw_surface = TRUE;
      if (strcmp (optarg, "solid")) {
	FILE * fp = fopen (optarg, "rt");
	GtsFile * f;
	
	if (fp == NULL) {
	  fprintf (stderr, 
		   "gfs2oogl: cannot open file `%s'\n"
		   "Try `gfs2oogl --help' for more information.\n",
		   optarg);
	  return 1; /* failure */
	}
	f = gts_file_new (fp);
	surface = gts_surface_new (gts_surface_class (),
				   gts_face_class (),
				   gts_edge_class (),
				   gts_vertex_class ());
	if (gts_surface_read (surface, f)) {
	  fprintf (stderr, "gfs2oogl: file `%s' is not a valid GTS file\n", 
		   optarg);
	  fprintf (stderr, "%s:%d:%d: %s\n",
		   optarg, f->line, f->pos, f->error);
	  return 1; /* failure */
	}
	gts_file_destroy (f);
	fclose (fp);
      }
      break;
    case 'l': /* lines */
      lines = TRUE;
      /* fall through */
    case 'R': /* ribbon */
      ribbon = TRUE;
      /* fall through */
    case 'C': /* cylinder */
      stream = fopen (optarg, "rt");
      streamname = g_strdup (g_basename (optarg));
      if (stream == NULL) {
	fprintf (stderr, 
		 "gfs2oogl: cannot open file `%s'\n"
		 "Try `gfs2oogl --help' for more information.\n",
		 optarg);
	return 1; /* failure */
      }
      break;
    case 'V': /* vector */
      vector = atof (optarg);
      break;
    case 'r': /* refine */
      refine = TRUE;
      break;
    case 'c': /* color */
      color = g_strdup (optarg);
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: gfs2oogl [OPTION] < GFS_FILE\n"
     "Converts a Gerris simulation file to other (graphical) formats.\n"
     "\n"
     "  -u V    --iso=V       outputs a GTS file isosurface for value V\n"
     "                        the variable needs to be specified using -c\n"
     "  -f D    --stream=D    draw evenly-spaced streamlines (D is the spacing)\n"
     "  -I M    --dmin=M      controls length of evenly-spaced streamlines\n"
     "                        default is 0.5\n"
     "  -O M    --maxcost=M   controls compression of streamlines (default is 2e-7)\n"
     "  -b x,.. --box=x,y,..  specify bounding box for streamline calculation\n"
     "  -j      --closed      outputs only closed streamlines\n"
     "  -p F    --profile=F   output list of values for coordinates defined in F\n"
     "  -o      --mixed       output text values in mixed cells only\n"
     "  -L L    --level=L     use cells at level L only\n"
     "  -i      --reinit      reinitializes refinement and solid fractions\n"
     "  -e      --merged      draw boundaries of merged cells\n"
     "  -S      --squares     draw (colored) squares\n"
     "  -g      --gnuplot     output gnuplot data\n"
     "  -x VAL  --sx=VAL      outputs a GTS surface, cross section for x = VAL\n"
     "                        of the scalar variable\n"
     "  -y VAL  --sy=VAL      outputs a GTS surface, cross section for y = VAL\n"
     "                        of the scalar variable\n"
     "  -z VAL  --sz=VAL      outputs a GTS surface, cross section for z = VAL\n"
     "                        of the scalar variable\n"
     "  -s S    --surface=S   outputs the surface defined by file S (or the solid\n"
     "                        surface is S is equal to `solid')\n"
     "  -V S    --vector=S    output an OOGL representation of the velocity vector\n"
     "                        field in the mixed cells\n"
     "  -l F    --streamlines=F  draw streamlines starting from each point defined\n"
     "                        in file F\n"
     "  -C F    --cylinder=F  draw stream cylinders starting from each point defined\n"
     "                        in file F\n"
     "  -R F    --ribbon=F    draw stream ribbons starting from each point defined\n"
     "                        in file F\n"
     "  -r                    refines the solid surface according to the local\n"
     "                        resolution\n"
     "  -c V    --color=V     color surfaces, streamlines etc... according to the\n"
     "  -m V    --min=V       set minimum scalar value to V\n"
     "  -M V    --max=V       set maximum scalar value to V\n"
     "  -v      --verbose     display statistics and other info\n"
     "  -h      --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfs2oogl --help' for more information.\n");
      return 1; /* failure */
    }
  }

  fp = gts_file_new (stdin);

  while (fp->type == GTS_INT) {
    GfsSimulation * simulation;
    GfsDomain * domain;
    GtsRange stats;
      
    if (!(simulation = gfs_simulation_read (fp))) {
      fprintf (stderr, 
	       "gfs2oogl: file on standard input is not a valid simulation file\n"
	       "<stdin>:%d:%d: %s\n",
	       fp->line, fp->pos, fp->error);
      return 1;
    }
    gfs_simulation_init (simulation);

    domain = GFS_DOMAIN (simulation);

    if (color) {
      GtsFile * fp = gts_file_new_from_string (color);
      GfsFunction * f = gfs_function_new (gfs_function_class (), 0.);
 
      gfs_function_read (f, domain, fp);
      gfs_pending_functions_compilation (fp);
      if (fp->type == GTS_ERROR) {
	fprintf (stderr, 
		 "gfs2oogl: incorrect `color' argument\n"
		 "%d: %s\n",
		 fp->pos, fp->error);
	return 1;
      }
      gts_file_destroy (fp);
      g_free (color);
      
      if (!(var = gfs_function_get_variable (f))) {
	gpointer data[2];

	data[0] = var = gfs_temporary_variable (domain);
	data[1] = f;
	gfs_domain_cell_traverse (domain,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) update_var, data);
      }
      gts_object_destroy (GTS_OBJECT (f));
    }

    if (verbose)
      fprintf (stderr, "gfs2oogl: processing t = %10e\n", simulation->time.t);

    if (reinit) {
      gfs_clock_start (domain->timer);
      gfs_simulation_refine (simulation);
      gfs_simulation_init (simulation);
      gfs_clock_stop (domain->timer);
    }

    if (var != NULL) {
      if (min == max) {
	stats = gfs_domain_stats_variable (domain, var, FTT_TRAVERSE_ALL, -1,
					   NULL, NULL);
	if (verbose)
	  fprintf (stderr, 
		   "min: %g avg: %g| %g max: %g n: %7d\n",
		   stats.min, stats.mean, stats.stddev, stats.max, stats.n);
      }
      else {
	stats.min = min;
	stats.max = max;
      }	
    }
    else
      stats.min = stats.max = 0.;

    if (var != NULL && gnuplot) {
      if (mixed)
	gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				   (FttCellTraverseFunc) output_mixed_variable, var);
      else {
	if (level < 0)
	  gfs_write_gnuplot (domain, var, 
			     FTT_TRAVERSE_LEAFS, -1, box, stdout);
	else
	  gfs_write_gnuplot (domain, var, 
			     FTT_TRAVERSE_LEVEL, level, box, stdout);
      }
    }
    else if (box && var != NULL) {
      if (squares) {
	if (level < 0)
	  gfs_write_squares (domain, var, stats.min, stats.max, 
			     FTT_TRAVERSE_LEAFS, -1, box, stdout);
	else
	  gfs_write_squares (domain, var, stats.min, stats.max, 
			     FTT_TRAVERSE_LEVEL, level, box, stdout);
      }
      else
	gfs_write_gts (domain, var, FTT_TRAVERSE_LEAFS, -1, box, stdout);
    }
    else if (stream) {
      FttVector p;

      rewind (stream);
      printf ("(geometry \"%s-%g\" = LIST {\n",
	      streamname,
	      simulation->time.t);
      while (fscanf (stream, "%lf %lf %lf", &p.x, &p.y, &p.z) == 3) {
#if FTT_2D
	gfs_draw_streamline (domain, p, stdout);
#else /* 3D */
	if (lines)
	  gfs_draw_streamline (domain, p, stdout);
	else if (ribbon)
	  gfs_draw_stream_ribbon (domain, p, 2e-3,
				  var, stats.min, stats.max, stdout);
	else
	  gfs_draw_stream_cylinder (domain, p, 5e-4, 
				    var, stats.min, stats.max, stdout);
#endif /* 3D */
      }
      printf ("})\n");
    }
    else if (profile) {
      FttVector p;

      if (var)
	while (fscanf (profile, "%lf %lf %lf", &p.x, &p.y, &p.z) == 3) {
	  FttVector pm = p;
	  gfs_simulation_map (simulation, &pm);
	  FttCell * cell = gfs_domain_locate (domain, pm, -1, NULL);
	  if (cell)
	    printf ("%g %g %g %g\n", p.x, p.y, p.z, gfs_dimensional_value (var, gfs_interpolate (cell, pm, var)));
	}
      else {
	GSList * j;
	guint i = 4;

	printf ("# 1:X 2:Y 3:Z ");
	j = domain->variables;
	while (j) {
	  GfsVariable * v = j->data;
	  printf ("%d:%s ", i++, v->name);
	  j = j->next;
	}
	printf ("\n");
	while (fscanf (profile, "%lf %lf %lf", &p.x, &p.y, &p.z) == 3) {
	   FttVector pm = p;
	  gfs_simulation_map (simulation, &pm);
	  FttCell * cell = gfs_domain_locate (domain, pm, -1, NULL);
	  if (cell) {
	    printf ("%g %g %g ", p.x, p.y, p.z);
	    j = domain->variables;
	    while (j) {
	      GfsVariable * v = j->data;
	      printf ("%g ", gfs_dimensional_value (v, gfs_interpolate (cell, pm, v)));
	      j = j->next;
	    }
	    printf ("\n");
	  }
	}
      }
    }
    else if (vector > 0.) {
      GtsRange stats;
      gdouble scale = 1.;
      GfsVariable * norm = gfs_temporary_variable (domain);
      gpointer data[2];

      data[0] = norm;
      data[1] = gfs_domain_velocity (domain);
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) velocity_norm, data);
      stats = gfs_domain_stats_variable (domain, norm, FTT_TRAVERSE_LEAFS, -1,
					 NULL, NULL);
      gts_object_destroy (GTS_OBJECT (norm));
      if (verbose)
	fprintf (stderr, 
		 "min: %g avg: %g| %g max: %g n: %7d\n",
		 stats.min, stats.mean, stats.stddev, stats.max, stats.n);
      if (stats.max > 0.)
	scale = vector*ftt_level_size (gfs_domain_depth (domain))/stats.max;
      printf ("(geometry \"vector-%g\" = LIST {\n", simulation->time.t);
      data[0] = &scale;
#if FTT_2D
      if (box == NULL)
	box = gts_bbox_new (gts_bbox_class (), NULL,
			    0., -G_MAXDOUBLE/2., -G_MAXDOUBLE/2.,
			    0., G_MAXDOUBLE/2., G_MAXDOUBLE/2.);
#else /* 3D */
      if (box == NULL)
	gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				  (FttCellTraverseFunc) draw_vector, data);
      else
#endif /* 3D */
      gfs_domain_cell_traverse_box (domain, box,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttCellTraverseFunc) draw_vector, data);
      printf ("})\n");
    }
    else if (draw_surface) {
      GSList * l = gfs_simulation_get_solids (simulation), * i = l;
	
      while (i) {
	GtsSurface * s = GFS_IS_SURFACE (GFS_SOLID (i->data)->s) ?  
	  GFS_SURFACE (GFS_SOLID (i->data)->s)->s : NULL;
  
	if (s) {
	  if (refine)
	    gts_surface_refine (s, 
				(GtsKeyFunc) local_size_ratio, domain,
				NULL, NULL,
				(GtsStopFunc) stop, NULL);
	  gfs_draw_surface (domain, s, 
			    var, stats.min, stats.max,
			    stdout);
	}
	i = i->next;
      }
      g_slist_free (l);
    }
    else if (merged) {
      gfs_set_merged (domain);
      puts ("LIST {\n");
      gfs_domain_traverse_merged (domain,
				  (GfsMergedTraverseFunc) merged_draw, 
				  stdout);
      puts ("}\n");
    }
    else if (mixed && var->name && !strcmp (var->name, "Vorticity")) {
      FttComponent c;
      GfsVariable ** u, * vort = gfs_temporary_variable (domain);
      gpointer data[2];

      u = gfs_domain_velocity (domain);
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) gfs_cell_reset, vort);
      data[0] = vort;
      for (c = 0; c < FTT_DIMENSION; c++) {
	gfs_domain_surface_bc (domain, u[c]);
	data[1] = u[c];
	gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				   (FttCellTraverseFunc) compute_mixed_vorticity, data);
      }
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) output_mixed_vorticity, vort);
      gts_object_destroy (GTS_OBJECT (vort));
    }
    else if (mixed && var->name && !strcmp (var->name, "P"))
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) output_mixed_pressure, var);
    else if (even_stream > 0.) {
      GList * s, * i;
      ClosestGrid * grid;

      gfs_domain_cell_traverse (domain,
				FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) gfs_cell_coarse_init, domain);

      if (bmax.x == -G_MAXDOUBLE)
	grid = closest_grid_new (domain, even_stream);
      else
	grid = closest_box_grid_new (bmin, bmax, even_stream);
      i = s = even_streamlines (domain, grid, even_stream, 
				even_stream*rdmin, closed);
      closest_grid_destroy (grid);
      printf ("LIST {\n");
      while (i) {
	i->data = simplify_stream (i->data, maxcost);
	write_stream (i->data, stdout);
	i = i->next;
      }
      printf ("}\n");
    }
    else if (iso < G_MAXDOUBLE && var != NULL) {
      GtsSurface * s = gfs_isosurface (domain, var, iso, level);

      gts_surface_write (s, stdout);
      gts_object_destroy (GTS_OBJECT (s));
    }
    else {
      gfs_draw_refined_boundaries (domain, stdout);
      gfs_draw_solid_boundaries (domain, stdout);
      gfs_draw_boundary_conditions (domain, stdout);
    }

    gts_object_destroy (GTS_OBJECT (simulation));
  }

  gts_file_destroy (fp);

  return 0;
}
