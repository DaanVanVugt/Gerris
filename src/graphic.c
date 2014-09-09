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
/*! \file
 * \brief Graphical utility functions.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#  ifdef HAVE_MPI
#    include <mpi.h>
#  endif
#endif

#include <stdlib.h>
#include <math.h>
#include <gts.h>

#include "config.h"
#include "graphic.h"
#include "solid.h"
#include "variable.h"
#include "version.h"
#include "init.h"

#if !FTT_2D
#  include "isocube.h"
#endif /* 3D */

typedef struct {
  GPtrArray * colors;
  gboolean reversed;
} Colormap;

static GtsColor * color_new (gdouble r, gdouble g, gdouble b)
{
  GtsColor * c = g_malloc (sizeof (GtsColor));
  c->r = r; c->g = g; c->b = b;
  return c;
}

static void color_destroy (GtsColor * color)
{
  g_return_if_fail (color != NULL);

  g_free (color);
}

static Colormap * colormap_jet (void)
{
  Colormap * cmap = g_malloc (sizeof (Colormap));
  gint i;

  cmap->reversed = FALSE;
  cmap->colors = g_ptr_array_new ();
  for (i = 0; i < 127; i++) {
    gdouble r = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    gdouble g = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    gdouble b =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;

    g_ptr_array_add (cmap->colors, color_new (r, g, b));
  }
  return cmap;
}

static void colormap_destroy (Colormap * colormap)
{
  guint i;

  g_return_if_fail (colormap != NULL);

  for (i = 0; i < colormap->colors->len; i++)
    color_destroy (colormap->colors->pdata[i]);
  g_ptr_array_free (colormap->colors, TRUE);
  g_free (colormap);
}

static GtsColor colormap_color (Colormap * cmap, gdouble val)
{
  GtsColor c = {1., 1., 1.}, * c1, * c2;
  guint i, n;
  gdouble coef;

  g_return_val_if_fail (cmap != NULL, c);

  if (val > 1.0) val = 1.0;
  else if (val < 0.0) val = 0.0;
  if (cmap->reversed)
    val = 1.0 - val;

  n = cmap->colors->len;
  if (n == 0)
    return c;
  if (n == 1)
    return *((GtsColor *)cmap->colors->pdata[0]);

  i = floor ((gdouble)val*(gdouble)(n - 1));
  if (i == n - 1)
    return *((GtsColor *)cmap->colors->pdata[cmap->colors->len - 1]);
  coef = val*(gdouble)(n - 1) - (gdouble)i;
  c1 = cmap->colors->pdata[i];
  c2 = cmap->colors->pdata[i+1];
  c.r = c1->r + coef*(c2->r - c1->r);
  c.g = c1->g + coef*(c2->g - c1->g);
  c.b = c1->b + coef*(c2->b - c1->b);
  return c;
}

/* VertexCellFace: Header */

typedef struct _VertexCellFace         VertexCellFace;

struct _VertexCellFace {
  /*< private >*/
  GtsVertex parent;

  /*< public >*/
  guint index;
  FttCell * cell;
  FttCellFace face;
};

#define VERTEX_CELL_FACE(obj)            GTS_OBJECT_CAST (obj,\
					         VertexCellFace,\
					         vertex_cell_face_class ())

/* VertexCellFace: Object */

static GtsVertexClass * vertex_cell_face_class (void)
{
  static GtsVertexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo vertex_cell_face_info = {
      "VertexCellFace",
      sizeof (VertexCellFace),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_vertex_class ()),
				  &vertex_cell_face_info);
  }

  return klass;
}

static void triangulate (FttCell * cell, gpointer * data)
{
  GtsSurface * s = data[0];
  GfsVariable * var = data[1];
  GtsVertex * v;
  FttVector pos;

  if (var && var->centered)
    ftt_cell_pos (cell, &pos);
  else 
    gfs_cell_cm (cell, &pos);

  v = gts_vertex_new (s->vertex_class, pos.x, pos.y, pos.z);
  if (var) {
    GtsMatrix * transform = data[2];
    gdouble * z = data[3];
    GfsNorm * norm = data[4];

    gts_point_transform (GTS_POINT (v), transform);
    GTS_POINT (v)->z = *z + GFS_VALUE (cell, var)/(norm->infty*1000.);
  }
  g_assert (gts_delaunay_add_vertex (s, v, NULL) == NULL);
  VERTEX_CELL_FACE (v)->cell = cell;
}

static void triangulate_face (FttCell * cell, gpointer * data)
{
  GtsSurface * s = data[0];
  GfsVariable * var = data[1];
  FttDirection * d = data[5];
  GtsVertex * v;
  FttVector pos;
  FttCellFace face;

  face.cell = cell;
  face.d = *d;
  face.neighbor = ftt_cell_neighbor (cell, face.d);
  ftt_face_pos (&face, &pos);
  v = gts_vertex_new (s->vertex_class, pos.x, pos.y, pos.z);
  if (var) {
    GtsMatrix * transform = data[2];
    gdouble * z = data[3];
    GfsNorm * norm = data[4];

    gts_point_transform (GTS_POINT (v), transform);
    if (face.neighbor)
      GTS_POINT (v)->z = *z + 
	gfs_face_interpolated_value_generic  (&face, var)/(norm->infty*1000.);
    else
      GTS_POINT (v)->z = *z + GFS_VALUE (cell, var)/(norm->infty*1000.);
  }
  g_assert (gts_delaunay_add_vertex (s, v, NULL) == NULL);
  VERTEX_CELL_FACE (v)->face = face;
}

static void add_long_segment (GtsSegment * s, GSList ** list)
{
  FttCell * c1 = VERTEX_CELL_FACE (s->v1)->cell ? VERTEX_CELL_FACE (s->v1)->cell :
    VERTEX_CELL_FACE (s->v1)->face.cell;
  FttCell * c2 = VERTEX_CELL_FACE (s->v2)->cell ? VERTEX_CELL_FACE (s->v2)->cell :
    VERTEX_CELL_FACE (s->v2)->face.cell;
  gdouble s1 = ftt_cell_size (c1);
  gdouble s2 = ftt_cell_size (c2);
  gdouble size = MIN (s1, s2);

  if (gts_point_distance2 (GTS_POINT (s->v1), GTS_POINT (s->v2)) > 
      16.*size*size)
    *list = g_slist_prepend (*list, s);
}

void gfs_write_gts (GfsDomain * domain, 
		    GfsVariable * v, 
		    FttTraverseFlags flags,
		    gint level,
		    GtsBBox * box,
		    FILE * fp)
{
  GtsSurface * s;
  GtsVertex * v1, * v2, * v3;
  GtsEdge * e1, * e2, * e3;
  gpointer data[6];
  GtsMatrix * transform, * inv;
  gdouble z = 0.;
  GfsNorm norm;
  GSList * long_segments = NULL;

  g_return_if_fail (domain != NULL);
#if (!FTT_2D)
  g_return_if_fail (box != NULL);
#endif /* 3D */
  g_return_if_fail (fp != NULL);

  v1 = gts_vertex_new (gts_vertex_class (), -100., -100., 0.);
  v2 = gts_vertex_new (gts_vertex_class (), 100., -100., 0.);
  v3 = gts_vertex_new (gts_vertex_class (), 0., 100., 0.);
  e1 = gts_edge_new (gts_edge_class (), v1, v2);
  e2 = gts_edge_new (gts_edge_class (), v2, v3);
  e3 = gts_edge_new (gts_edge_class (), v3, v1);
  s = gts_surface_new (gts_surface_class (), 
		       gts_face_class (), 
		       gts_edge_class (), 
		       vertex_cell_face_class ());
  gts_surface_add_face (s, gts_face_new (gts_face_class (), e1, e2, e3));

  norm = gfs_domain_norm_variable (domain, v, NULL, flags, level, NULL, NULL);
  if (norm.infty == 0.)
    norm.infty = 1.;
#if FTT_2D
  transform = gts_matrix_identity (NULL);
#else /* 3D */
  if (box->x2 - box->x1 < box->z2 - box->z1 &&
      box->x2 - box->x1 < box->y2 - box->y1) {
    z = box->x2 = box->x1 = (box->x2 + box->x1)/2. + 1e-30;
    transform = gts_matrix_new (0., 1., 0., 0.,
				0., 0., 1., 0.,
				1., 0., 0., 0.,
				0., 0., 0., 0.);
  }
  else if (box->y2 - box->y1 < box->z2 - box->z1 &&
	   box->y2 - box->y1 < box->x2 - box->x1) {
    z = box->y2 = box->y1 = (box->y2 + box->y1)/2. + 1e-30;
    transform = gts_matrix_new (0., 0., 1., 0.,
				1., 0., 0., 0.,
				0., 1., 0., 0.,
				0., 0., 0., 0.);
  }
  else {
    z = box->z2 = box->z1 = (box->z2 + box->z1)/2. + 1e-30;
    transform = gts_matrix_new (1., 0., 0., 0.,
				0., 1., 0., 0.,
				0., 0., 1., 0.,
				0., 0., 0., 0.);
  }
#endif /* 3D */
  
  data[0] = s;
  data[1] = v;
  data[2] = transform;
  data[3] = &z;
  data[4] = &norm;
  if (box == NULL) {
    FttDirection d;

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level, 
			     (FttCellTraverseFunc) triangulate, data);
    data[5] = &d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      gfs_domain_cell_traverse_boundary (domain, d, 
					 FTT_PRE_ORDER, flags, level, 
		   (FttCellTraverseFunc) triangulate_face, data);
  }
  else
    gfs_domain_cell_traverse_box (domain, box, FTT_PRE_ORDER, flags, level, 
				 (FttCellTraverseFunc) triangulate, data);

  gts_allow_floating_vertices = TRUE;
  gts_object_destroy (GTS_OBJECT (v1));
  gts_object_destroy (GTS_OBJECT (v2));
  gts_object_destroy (GTS_OBJECT (v3));
  gts_allow_floating_vertices = FALSE;

  gts_surface_foreach_edge (s, (GtsFunc) add_long_segment, &long_segments);
  gts_allow_floating_edges = TRUE;
  g_slist_foreach (long_segments, (GFunc) gts_object_destroy, NULL);
  gts_allow_floating_edges = FALSE;
  g_slist_free (long_segments);

  inv = gts_matrix3_inverse (transform);
  gts_matrix_destroy (transform);
  transform = inv;
  gts_surface_foreach_vertex (s, (GtsFunc) gts_point_transform, transform);
  gts_surface_write (s, fp);

  gts_object_destroy (GTS_OBJECT (s));
  gts_matrix_destroy (transform);
}

static void extent (FttCell * cell, gpointer * data)
{
  FttVector * min = data[0];
  FttVector * max = data[1];
  FttVector pos;
  
  ftt_cell_pos (cell, &pos);
  if (pos.x > max->x) max->x = pos.x;
  if (pos.y > max->y) max->y = pos.y;
  if (pos.z > max->z) max->z = pos.z;
  if (pos.x < min->x) min->x = pos.x;
  if (pos.y < min->y) min->y = pos.y;
  if (pos.z < min->z) min->z = pos.z;
}

static void iso_func (gdouble ** a, GtsCartesianGrid g, guint k,
		      gpointer * data)
{
  GfsDomain * domain = data[0];
  guint * level = data[1], i, j;
  GfsVariable * v = data[2];
  FttVector p;
  fprintf (stderr, "\rslice %4d/%d", k + 1, g.nz);
  p.z = g.z;
  for (i = 0, p.x = g.x; i < g.nx; i++, p.x += g.dx)
    for (j = 0, p.y = g.y; j < g.ny; j++, p.y += g.dy) {
      FttCell * cell = gfs_domain_locate (domain, p, *level, NULL);

      if (cell == NULL)
	a[i][j] = 0.;
      else
	a[i][j] = gfs_interpolate (cell, p, v);
    }
}

GtsSurface * gfs_isosurface (GfsDomain * domain, 
			     GfsVariable * v, gdouble val,
			     gint level)
{
  FttVector cmax = { - G_MAXDOUBLE, - G_MAXDOUBLE, - G_MAXDOUBLE };
  FttVector cmin = { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE };
  guint depth;
  GtsCartesianGrid g;
  gpointer data[3];
  GtsSurface * s;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (v != NULL, NULL);

  if (level < 0)
    depth = gfs_domain_depth (domain);
  else
    depth = level;

  data[0] = &cmin;
  data[1] = &cmax;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, depth,
			    (FttCellTraverseFunc) extent, data);
  if (cmin.x == G_MAXDOUBLE)
    return NULL;

  g.dx = g.dy = g.dz = ftt_level_size (depth);
  g.x = cmin.x; g.y = cmin.y; g.z = cmin.z;
  g.nx = (cmax.x - cmin.x)/g.dx + 1;
  g.ny = (cmax.y - cmin.y)/g.dy + 1;
  g.nz = (cmax.z - cmin.z)/g.dz + 1;
  
  s = gts_surface_new (gts_surface_class (), 
		       gts_face_class (), 
		       gts_edge_class (), 
		       gts_vertex_class ());
  data[0] = domain;
  data[1] = &depth;
  data[2] = v;
  gts_isosurface_cartesian (s, g, (GtsIsoCartesianFunc) iso_func, data, val);

  return s;
}

static void write_gnuplot (FttCell * cell, gpointer * data)
{
  FILE * fp = data[0];
  GfsVariable * v = data[1];
  GtsBBox * bbox = data[2];
  FttVector pos;
  
  if (v->centered)
    ftt_cell_pos (cell, &pos);
  else
    gfs_cell_cm (cell, &pos);

  if (bbox == NULL || (pos.x >= bbox->x1 && pos.x <= bbox->x2 &&
		       pos.y >= bbox->y1 && pos.y <= bbox->y2 &&
		       pos.z >= bbox->z1 && pos.z <= bbox->z2)) {
    gfs_simulation_map_inverse (GFS_SIMULATION (v->domain), &pos);
    fprintf (fp, "%g %g %g %g\n", 
	     pos.x, pos.y, pos.z, GFS_VALUE (cell, v));
  }
}

void gfs_write_gnuplot (GfsDomain * domain, 
			GfsVariable * v, 
			FttTraverseFlags flags,
			gint level,
			GtsBBox * bbox,
			FILE * fp)
{  
  gpointer data[3];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  data[0] = fp;
  data[1] = v;
  data[2] = bbox;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level, 
			   (FttCellTraverseFunc) write_gnuplot, data);
}

typedef struct {
  guchar r, g, b;
} Color;

typedef struct {
  FttVector min;
  guint width, height, size;
  guchar * buf, *** data;
} Image;

static Image * image_new (FttVector min, FttVector max, guint size)
{
  Image * im = g_malloc0 (sizeof (Image));
  guint i;

  im->min = min;
  im->size = size;
  im->width = (max.x - min.x)*size;
  im->height = (max.y - min.y)*size;
  im->buf = g_malloc0 (sizeof (guchar)*3*im->width*im->height);
  im->data = g_malloc (sizeof (guchar **)*im->height);
  for (i = 0; i < im->height; i++) {
    guint j;

    im->data[i] = g_malloc (sizeof (guchar *)*im->width);
    for (j = 0; j < im->width; j++)
      im->data[i][j] = &im->buf[3*(i*im->width + j)];
  }
  return im;
}

static void image_write (Image * im, FILE * fp)
{
  fprintf (fp, 
	   "P6\n"
	   "# File generated by gerris using 2D libgfs version %s (%s)\n"
	   "# Origin: %d %d\n"
	   "%u %u 255\n",
	   GFS_VERSION,
	   GFS_BUILD_VERSION,
	   (gint) (im->min.x*im->size), 
	   (gint) (im->min.y*im->size),
	   im->width, im->height);
  fwrite (im->buf, sizeof (guchar), 3*im->width*im->height, fp);
}

static void image_destroy (Image * im)
{
  guint i;

  for (i = 0; i < im->height; i++)
    g_free (im->data[i]);
  g_free (im->data);
  g_free (im->buf);
  g_free (im);
}

static void image_draw_square (Image * im,
			       FttVector * p1, FttVector * p2,
			       Color c)
{
  gint i1, j1, i2, j2, i, j;

  i1 = (p1->x - im->min.x)*im->size;
  i2 = (p2->x - im->min.x)*im->size;
  j1 = (p1->y - im->min.y)*im->size;
  j2 = (p2->y - im->min.y)*im->size;

  j1 = im->height - 1 - j1;
  j2 = im->height - 1 - j2;
  for (i = i1; i <= i2; i++)
    for (j = j2; j <= j1; j++) 
      if (i >= 0 && i < im->width && j >= 0 && j < im->height) {
	im->data[j][i][0] = c.r;
	im->data[j][i][1] = c.g;
	im->data[j][i][2] = c.b;
      }
}

static void write_image_square (FttCell * cell, gpointer * data)
{
  Colormap * colormap = data[0];
  gdouble * min = data[1];
  gdouble * max = data[2];
  GfsVariable * v = data[3];
  Image * image = data[4];
  FttVector * lambda = data[5];
  FttVector p;
  GtsColor fc = { 0., 0., 0. }; /* nodata = black */
  if (GFS_HAS_DATA (cell, v))
    fc = colormap_color (colormap, (GFS_VALUE (cell, v) - *min)/(*max - *min));
  Color c;
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector p1, p2;

  ftt_cell_pos (cell, &p);
  c.r = fc.r*255;
  c.g = fc.g*255;
  c.b = fc.b*255;
  p1.x = (p.x - size)/lambda->x + 1e-9;
  p1.y = (p.y - size)/lambda->y + 1e-9;
  p2.x = (p.x + size)/lambda->x - 1e-9;
  p2.y = (p.y + size)/lambda->y - 1e-9;
  image_draw_square (image, &p1, &p2, c);
}

static void max_extent (FttCell * cell, FttVector extent[2])
{
  gdouble h = ftt_cell_size (cell)/2.;
  FttVector pos;
  
  ftt_cell_pos (cell, &pos);
  if (pos.x - h < extent[0].x) extent[0].x = pos.x - h;
  if (pos.y - h < extent[0].y) extent[0].y = pos.y - h;
  if (pos.z - h < extent[0].z) extent[0].z = pos.z - h;

  if (pos.x + h > extent[1].x) extent[1].x = pos.x + h;
  if (pos.y + h > extent[1].y) extent[1].y = pos.y + h;
  if (pos.z + h > extent[1].z) extent[1].z = pos.z + h;
}

static gboolean cell_condition (FttCell * cell, gpointer condition)
{
  return gfs_function_value (condition, cell);
}

void gfs_write_ppm (GfsDomain * domain, 
		    GfsFunction * condition,
		    GfsVariable * v, gdouble min, gdouble max,
		    FttTraverseFlags flags,
		    gint level,
		    FILE * fp,
		    gboolean parallel)
{
  Colormap * colormap;
  guint depth, size = 1;
  Image * image;
  FttVector extent[2] = {{ G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE },
			 { - G_MAXDOUBLE, - G_MAXDOUBLE, - G_MAXDOUBLE }};
  gpointer data[6];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  if (min == max)
    max = min + 1.;
  if (level < 0)
    depth = gfs_domain_depth (domain);
  else
    depth = level;
  while (depth-- > 0)
    size *= 2;

  if (condition) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, level,
					(FttCellTraverseFunc) max_extent, extent,
					cell_condition, condition);
    gfs_restore_fpe_for_function (condition);
  }
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			      (FttCellTraverseFunc) max_extent, extent);

  gfs_all_reduce (domain, extent[0].x, MPI_DOUBLE, MPI_MIN);
  gfs_all_reduce (domain, extent[0].y, MPI_DOUBLE, MPI_MIN);
  gfs_all_reduce (domain, extent[1].x, MPI_DOUBLE, MPI_MAX);
  gfs_all_reduce (domain, extent[1].y, MPI_DOUBLE, MPI_MAX);
    
  if (extent[0].x == G_MAXDOUBLE)
    return;

  extent[0].x /= domain->lambda.x; 
  extent[0].y /= domain->lambda.y;
  extent[1].x /= domain->lambda.x; 
  extent[1].y /= domain->lambda.y;

  colormap = colormap_jet ();
  image = image_new (extent[0], extent[1], size);

  data[0] = colormap;
  data[1] = &min;
  data[2] = &max;
  data[3] = v;
  data[4] = image;
  data[5] = &domain->lambda;
  if (condition) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, level,
					(FttCellTraverseFunc) write_image_square, data,
					cell_condition, condition);
    gfs_restore_fpe_for_function (condition);
  }
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			      (FttCellTraverseFunc) write_image_square, data);

#ifdef HAVE_MPI
  if (!parallel && domain->pid >= 0) {
    if (domain->pid == 0) {
      Image * im = image_new (extent[0], extent[1], size);
      int n, np;
      MPI_Comm_size (MPI_COMM_WORLD, &np);
      for (n = 1; n < np; n++) {
	MPI_Status status;
	MPI_Recv (im->buf, 3*image->width*image->height, MPI_BYTE, n, 0, MPI_COMM_WORLD, &status);
	int i, j;
	for (i = 0; i < im->height; i++)
	  for (j = 0; j < im->width; j++)
	    if (im->data[i][j][0] || im->data[i][j][1] || im->data[i][j][2]) {
	      image->data[i][j][0] = im->data[i][j][0];
	      image->data[i][j][1] = im->data[i][j][1];
	      image->data[i][j][2] = im->data[i][j][2];
	    }
      }
      image_destroy (im);
      image_write (image, fp);      
    }
    else if (domain->pid > 0)
      MPI_Send (image->buf, 3*image->width*image->height, MPI_BYTE, 0, 0, MPI_COMM_WORLD);    
  }
  else
#endif
    image_write (image, fp);
  image_destroy (image);
  colormap_destroy (colormap);
}

#define NODATA -9999

typedef struct {
  FttVector min;
  gdouble cellsize;
  guint width, height;
  gfloat * buf, ** data;
} Grid;

static Grid * grid_new (FttVector min, FttVector max, FttVector step) 
{
  Grid * im = g_malloc0 (sizeof (Grid));
  guint i;

  im->min = min;
  im->cellsize = MIN (step.x, step.y);
  im->width = (max.x - min.x)/im->cellsize;
  im->height = (max.y - min.y)/im->cellsize;
  im->buf = g_malloc (sizeof (gfloat)*im->width*im->height);
  for (i = 0; i < im->height*im->width; i++)
    im->buf[i] = NODATA;
  im->data = g_malloc (sizeof (gfloat *)*im->height);
  for (i = 0; i < im->height; i++)
    im->data[i] = &im->buf[i*im->width];
  return im;
}

static void grid_write (Grid * im, FILE * fp)
{
  fprintf (fp, 
	   "ncols\t\t%d\n"
	   "nrows\t\t%d\n"
	   "xllcorner\t%f\n"
	   "yllcorner\t%f\n"
	   "cellsize\t%.10f\n"
	   "nodata_value\t%d\n",
	   im->width, im->height, 
	   im->min.x, im->min.y, im->cellsize,
	   NODATA);
  guint i, j;
  for (i = 0; i < im->height; i++)
    for (j = 0; j < im->width; j++)
      fprintf (fp, "%g ", im->data[i][j]);
}

static void grid_destroy (Grid * im)
{
  g_free (im->data);
  g_free (im->buf);
  g_free (im);
}

static void max_physical_extent (FttCell * cell, gpointer * data)
{
  FttVector * extent = data[0];
  GfsSimulation * sim = data[1];
  gdouble h = ftt_cell_size (cell)/2.;
  FttVector p, 
    max = { -G_MAXDOUBLE, -G_MAXDOUBLE, -G_MAXDOUBLE },
    min = { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE };
  ftt_cell_pos (cell, &p);
  double x, y;
  for (x = -1; x <= 1.; x += 2.)
    for (y = -1; y <= 1.; y += 2.) {
      FttVector c = p;
      c.x += x*h; c.y += y*h;
      gfs_simulation_map_inverse (sim, &c);
      int i;
      for (i = 0; i < FTT_DIMENSION; i++) {
	if ((&c.x)[i] > (&max.x)[i]) (&max.x)[i] = (&c.x)[i];
	if ((&c.x)[i] < (&min.x)[i]) (&min.x)[i] = (&c.x)[i];
      }
    }

  int i;
  for (i = 0; i < 2; i++) {
    if ((&max.x)[i] > (&extent[1].x)[i]) (&extent[1].x)[i] = (&max.x)[i];
    if ((&min.x)[i] < (&extent[0].x)[i]) (&extent[0].x)[i] = (&min.x)[i];
    if ((&max.x)[i] - (&min.x)[i] < (&extent[2].x)[i]) 
      (&extent[2].x)[i] = (&max.x)[i] - (&min.x)[i];
  }
}

void gfs_write_grd (GfsSimulation * sim, 
		    GfsFunction * condition,
		    GfsVariable * v,
		    FttTraverseFlags flags,
		    gint level,
		    FILE * fp,
		    gboolean parallel,
		    gboolean interpolate)
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (fp != NULL);

  GfsDomain * domain = GFS_DOMAIN (sim);
  FttVector extent[3] = {{ G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE },
			 { - G_MAXDOUBLE, - G_MAXDOUBLE, - G_MAXDOUBLE },
			 { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE }};
  gpointer data[2] = { extent, sim };
  if (condition) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, level,
					(FttCellTraverseFunc) max_physical_extent, data,
					cell_condition, condition);
    gfs_restore_fpe_for_function (condition);
  }
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			      (FttCellTraverseFunc) max_physical_extent, data);
    
  gfs_all_reduce (domain, extent[0].x, MPI_DOUBLE, MPI_MIN);
  gfs_all_reduce (domain, extent[0].y, MPI_DOUBLE, MPI_MIN);
  gfs_all_reduce (domain, extent[1].x, MPI_DOUBLE, MPI_MAX);
  gfs_all_reduce (domain, extent[1].y, MPI_DOUBLE, MPI_MAX);
  gfs_all_reduce (domain, extent[2].x, MPI_DOUBLE, MPI_MIN);
  gfs_all_reduce (domain, extent[2].y, MPI_DOUBLE, MPI_MIN);
    
  if (extent[0].x == G_MAXDOUBLE)
    return;

  Grid * grid = grid_new (extent[0], extent[1], extent[2]);

  int i, j;
  for (i = 0; i < grid->width; i++)
    for (j = 0; j < grid->height; j++) {
      FttVector p = { grid->min.x + (0.5 + i)*grid->cellsize, 
		      grid->min.y + grid->height*grid->cellsize - (0.5 + j)*grid->cellsize, 
		      0. };
      gfs_simulation_map (sim, &p);
      FttCell * cell = gfs_domain_locate (domain, p, level, NULL);
      if (cell && GFS_HAS_DATA (cell, v))
	grid->data[j][i] = interpolate ? gfs_interpolate (cell, p, v) : GFS_VALUE (cell, v);
    }

#ifdef HAVE_MPI
  if (!parallel && domain->pid >= 0) {
    if (domain->pid == 0) {
      Grid * im = grid_new (extent[0], extent[1], extent[2]);
      int n, np;
      MPI_Comm_size (MPI_COMM_WORLD, &np);
      for (n = 1; n < np; n++) {
	MPI_Status status;
	MPI_Recv (im->buf, grid->width*grid->height, MPI_FLOAT, n, 0, MPI_COMM_WORLD, &status);
	int i, j;
	for (i = 0; i < im->height; i++)
	  for (j = 0; j < im->width; j++)
	    if (im->data[i][j] != NODATA)
	      grid->data[i][j] = im->data[i][j];
      }
      grid_destroy (im);
      grid_write (grid, fp);      
    }
    else if (domain->pid > 0)
      MPI_Send (grid->buf, grid->width*grid->height, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
  }
  else
#endif
    grid_write (grid, fp);
  grid_destroy (grid);
}

static gint gfs_combine_close (FILE ** f, Image ** im, gint n, gint ret)
{
  guint i;

  for (i = 0; i < n; i++) {
    if (f[i])
      fclose (f[i]);
    if (im[i])
      image_destroy (im[i]);
  }
  g_free (f);
  g_free (im);
  return ret;
}

static gint get_newline (FILE * fp)
{
  gint c;

  c = fgetc (fp);
  while (c != EOF && c != '\n')
    c = fgetc (fp);
  return c;
}

gint gfs_combine_ppm (gchar ** fname, guint nname, FILE * fp)
{
  FILE ** f;
  guint i;
  Image ** image;

  g_return_val_if_fail (fname != NULL, 0);
  g_return_val_if_fail (fp != NULL, 0);

  f = g_malloc0 (nname*sizeof (FILE *));
  image = g_malloc0 (nname*sizeof (Image *));
  for (i = 0; i < nname; i++) {
    f[i] = fopen (fname[i], "r");
    if (f[i] == NULL)
      return gfs_combine_close (f, image, nname, i);
  }
  
  while (TRUE) {
    gint x0 = G_MAXINT, y0 = G_MAXINT, x1 = - G_MAXINT, y1 = - G_MAXINT;
    FttVector min, max;
    Image * combo;

    for (i = 0; i < nname; i++) {
      gchar s[80];
      gint x, y;
      gint h, w;
      guint status;

      status = fscanf (f[i], "%79s", s);
      if (status != 1 && feof (f[i]))
	return gfs_combine_close (f, image, nname, -1);
      if (status != 1 ||
	  strcmp (s, "P6") ||
	  get_newline (f[i]) == EOF ||
	  get_newline (f[i]) == EOF ||
	  fscanf (f[i], "%*s %79s %d %d", s, &x, &y) != 3 ||
	  strcmp (s, "Origin:") ||
	  fscanf (f[i], "%d %d", &w, &h) != 2)
	return gfs_combine_close (f, image, nname, i);
      if (x < x0) x0 = x;
      if (y < y0) y0 = y;
      if (x + w > x1) x1 = x + w;
      if (y + h > y1) y1 = y + h;
      min.x = x;
      min.y = y;
      max.x = x + w;
      max.y = y + h;
      if (image[i] != NULL)
	image_destroy (image[i]);
      image[i] = image_new (min, max, 1);
      if (get_newline (f[i]) == EOF ||
	  fread (image[i]->buf, sizeof (guchar), 
		 3*image[i]->width*image[i]->height, f[i]) !=
	  3*image[i]->width*image[i]->height)
	return gfs_combine_close (f, image, nname, i);
    }

    min.x = x0;
    min.y = y0;
    max.x = x1;
    max.y = y1;
    combo = image_new (min, max, 1);
    for (i = 0; i < nname; i++) {
      guint x, y;

      for (y = 0; y < image[i]->height; y++)
	for (x = 0; x < image[i]->width; x++) {
	  gint x1 = x + image[i]->min.x - combo->min.x;
	  gint y1 = y + combo->min.y + combo->height - image[i]->min.y - image[i]->height;
	  guchar r = image[i]->data[y][x][0];
	  guchar g = image[i]->data[y][x][1];
	  guchar b = image[i]->data[y][x][2];

	  if (r || g || b) {
	    combo->data[y1][x1][0] = r;
	    combo->data[y1][x1][1] = g;
	    combo->data[y1][x1][2] = b;
	  }
	}
    }
    image_write (combo, fp);
    image_destroy (combo);
  }
  gfs_combine_close (f, image, nname, 0);
}

static void write_square (FttCell * cell, gpointer * data)
{
  Colormap * colormap = data[0];
  gdouble * min = data[1];
  gdouble * max = data[2];
  GfsVariable * v = data[3];
  FILE * fp = data[4];
  FttVector p;
  GtsColor c;
  gdouble size = ftt_cell_size (cell)/2.;

  ftt_cell_pos (cell, &p);
  c = colormap_color (colormap, 
		      (GFS_VALUE (cell, v) - *min)/(*max - *min));
#if FTT_2D    
  fprintf (fp, 
	   "OFF 4 1 4\n"
	   "%g %g 0\n%g %g 0\n%g %g 0\n%g %g 0\n"
	   "5 0 1 2 3 0 %g %g %g\n",
	   p.x - size, p.y - size,
	   p.x + size, p.y - size,
	   p.x + size, p.y + size,
	   p.x - size, p.y + size,
	   c.r, c.g, c.b);
#else  /* FTT_3D */
  fprintf (fp, 
	   "OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "4 3 2 1 0 %g %g %g\n"
	   "4 4 5 6 7 %g %g %g\n"
	   "4 2 3 7 6 %g %g %g\n"
	   "4 0 1 5 4 %g %g %g\n"
	   "4 0 4 7 3 %g %g %g\n"
	   "4 1 2 6 5 %g %g %g\n",
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size,
	   c.r, c.g, c.b,
	   c.r, c.g, c.b,
	   c.r, c.g, c.b,
	   c.r, c.g, c.b,
	   c.r, c.g, c.b,
	   c.r, c.g, c.b);
#endif /* FTT_3D */
}

void gfs_write_squares (GfsDomain * domain, 
			GfsVariable * v, gdouble min, gdouble max,
			FttTraverseFlags flags,
			gint level,
			GtsBBox * bbox,
			FILE * fp)
{
  Colormap * colormap;
  gpointer data[5];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  if (min == max)
    max = min + 1.;
  fputs ("LIST{\n", fp);
  colormap = colormap_jet ();
  data[0] = colormap;
  data[1] = &min;
  data[2] = &max;
  data[3] = v;
  data[4] = fp;  
  if (bbox == NULL)
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			     (FttCellTraverseFunc) write_square, data);
  else
    gfs_domain_cell_traverse_box (domain, bbox, FTT_PRE_ORDER, flags, level,
			     (FttCellTraverseFunc) write_square, data);
  fputs ("}\n", fp);
  colormap_destroy (colormap);
}

static void write_mac (FttCellFace * face, gpointer * data)
{
  gdouble * scale = data[0];
  FILE * fp = data[1];
  GtsBBox * bbox = data[2];
  FttVector p;

  ftt_face_pos (face, &p);
  if (bbox == NULL || (p.x >= bbox->x1 && p.x <= bbox->x2 &&
		       p.y >= bbox->y1 && p.y <= bbox->y2 &&
		       p.z >= bbox->z1 && p.z <= bbox->z2)) {
    FttVector f = {0., 0., 0.};
    gdouble un = GFS_FACE_NORMAL_VELOCITY (face)*(*scale);
    FttComponent c = face->d/2;

    switch (c) {
    case FTT_X: f.x = un; break;
    case FTT_Y: f.y = un; break;
#if (!FTT_2D)
    case FTT_Z: f.z = un; break;
#endif /* not FTT_2D */
    default: g_assert_not_reached ();
    }
    fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
	     p.x + f.x - (f.x - f.y/2.)/5.,
	     p.y + f.y - (f.x/2. + f.y)/5.,
	     p.z + f.z,
	     p.x + f.x,
	     p.y + f.y,
	     p.z + f.z,
	     p.x + f.x - (f.x + f.y/2.)/5.,
	     p.y + f.y + (f.x/2. - f.y)/5.,
	     p.z + f.z);
    fprintf (fp, "%g %g %g\n%g %g %g\n\n",
	     p.x, p.y, p.z,
	     p.x + f.x,
	     p.y + f.y,
	     p.z + f.z);
  }
}

void gfs_write_mac_velocity (GfsDomain * domain,
			     gdouble scale,
			     FttTraverseFlags flags,
			     gint level,
			     GtsBBox * bbox,
			     FILE * fp)
{
  gpointer data[3];
  GfsNorm norm;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  norm = gfs_domain_norm_velocity (domain, flags, level);
  scale = norm.infty > 0. ? 
    ftt_level_size (level < 0 ? gfs_domain_depth (domain) : level)*
    scale/norm.infty : scale;
  data[0] = &scale;
  data[1] = fp;
  data[2] = bbox;
  gfs_domain_face_traverse (domain, FTT_XYZ, FTT_PRE_ORDER, flags, level,
			   (FttFaceTraverseFunc) write_mac, data);
}

void gfs_draw_cells (FttCell * cell, 
		     FttTraverseFlags flags,
		     gint level,
		     FILE * fp)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp, "LIST {\n");
  ftt_cell_traverse (cell, FTT_PRE_ORDER, flags, level,
		     (FttCellTraverseFunc) ftt_cell_draw, fp);
  fprintf (fp, "}\n");
}

void gfs_draw_levels (FttCell * cell, FILE * fp)
{
  guint l, depth;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  depth = ftt_cell_depth (cell);
  for (l = 0; l <= depth; l++) {
    fprintf (fp, "(geometry \"level %d\" { = ", l);
    gfs_draw_cells (cell, FTT_TRAVERSE_LEVEL, l, fp);
    fputs ("})\n", fp);
  }
}

static void draw_box_boundaries (GfsBox * box, FILE * fp)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    GtsObject * o =  box->neighbor[d];

    if (GFS_IS_BOUNDARY (o)) {
      if (o->klass->color) {
	GtsColor c = (* o->klass->color) (o);

#if FTT_2D
	fprintf (fp, "appearance { material { edgecolor %g %g %g } }\n", 
		 c.r, c.g, c.b);
#else /* 3D */
	fprintf (fp, 
	  "appearance { material { ambient %g %g %g diffuse %g %g %g } }\n",
		 c.r, c.g, c.b, c.r, c.g, c.b);
#endif /* 3D */	
      }
      fputs ("LIST {\n", fp);
      ftt_face_traverse_boundary (box->root, d, 
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
				  (FttFaceTraverseFunc) ftt_face_draw, fp);
      fputs ("}\n", fp);
    }
    else if (GFS_IS_BOX (o) && box->pid != GFS_BOX (o)->pid) {
#if FTT_2D
      fputs ("appearance { material { edgecolor 1 0 0 } }\n", fp);
#else /* 3D */
      fputs ("appearance { material { ambient 1 0 0 diffuse 1 0 0 } }\n", fp);
#endif /* 3D */
      fputs ("LIST {\n", fp);
      ftt_face_traverse_boundary (box->root, d, 
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
				  (FttFaceTraverseFunc) ftt_face_draw, fp);
      fputs ("}\n", fp);
    }
  }
}

/**
 * gfs_draw_boundary_conditions:
 * @domain: a fluid domain.
 * @fp: a file pointer.
 *
 * Outputs in @fp an OOGL (geomview) representation of the boundary
 * conditions of the domain.  
 */
void gfs_draw_boundary_conditions (GfsDomain * domain, FILE * fp)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  fputs ("(geometry \"conditions\" = \n"
	 "LIST {\n", fp);
#if FTT_2D
  fputs ("appearance { linewidth 2 }\n", fp);
#endif /* 2D */
  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) draw_box_boundaries, fp);
  fputs ("})\n", fp);
}

static void draw_boundary_face (FttCell * cell, FILE * fp)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellFace face = ftt_cell_face (cell, d);

    if (ftt_face_type (&face) == FTT_BOUNDARY)
      ftt_face_draw (&face, fp);
  }
}

/**
 * gfs_draw_solid_boundaries:
 * @domain: a fluid domain.
 * @fp: a file pointer.
 *
 * Outputs in @fp an OOGL (geomview) representation of the solid boundaries
 * of the domain.
 */
void gfs_draw_solid_boundaries (GfsDomain * domain, FILE * fp)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  fputs ("(geometry \"solid\" = {\n", fp);
#if FTT_2D
  fputs ("appearance { linewidth 2 }\n", fp);
#endif /* 2D */
  fputs ("LIST{\n", fp);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) draw_boundary_face, fp);
  fputs ("}})\n", fp);
}

static void count_face (FttCell * cell, guint * count)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      FttCellFace face = ftt_cell_face (cell, d);
      
      if (ftt_face_type (&face) == FTT_FINE_COARSE && face.cell == cell)
	(*count)++;
    }
  }
}

static void draw_face (FttCell * cell, FILE * fp)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      FttCellFace face = ftt_cell_face (cell, d);
      
      if (ftt_face_type (&face) == FTT_FINE_COARSE && face.cell == cell)
	ftt_face_draw (&face, fp);
    }
  }
}

/**
 * gfs_draw_refined_boundaries:
 * @domain: a fluid domain.
 * @fp: a file pointer.
 *
 * Outputs in @fp an OOGL (geomview) representation of the boundaries
 * of the refined domains.
 */
void gfs_draw_refined_boundaries (GfsDomain * domain, FILE * fp)
{
  guint depth, level;
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  depth = gfs_domain_depth (domain);
  for (level = 1; level <= depth; level++) {
    guint count = 0;

    gfs_domain_cell_traverse (domain, 
			     FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, level,
			     (FttCellTraverseFunc) count_face, &count);
    if (count > 0) {
      fprintf (fp, "(geometry \"refine_%u_%u\" = \n", level - 1, level);
      fputs ("LIST{\n", fp);
      gfs_domain_cell_traverse (domain, 
			       FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, level,
			       (FttCellTraverseFunc) draw_face, fp);
      fputs ("}\n)\n", fp);
    }
  }
}

static gpointer color_data[5];

static GtsColor variable_color (GtsObject * o)
{
  GfsDomain * domain = color_data[0];
  GfsVariable * v = color_data[1];
  Colormap * colormap = color_data[2];
  gdouble * min = color_data[3];
  gdouble * max = color_data[4];
  FttCell * cell;

  GtsPoint * p = GTS_POINT (o);
  FttVector pos;
  gdouble val;
  GtsColor c;

  pos.x = p->x;
  pos.y = p->y;
  pos.z = p->z;

  cell = gfs_domain_locate (domain, pos, -1, NULL);
  if (cell) {
    val = gfs_interpolate (cell, pos, v);
    c = colormap_color (colormap, (val - *min)/(*max - *min));
  }
  else
    c.r = c.g = c.b = 1.;
  return c;
}

void gfs_draw_surface (GfsDomain * domain,
		       GtsSurface * s,
		       GfsVariable * v, 
		       gdouble min, gdouble max,
		       FILE * fp)
{
  GtsColor (*old_color) (GtsObject *);
  Colormap * colormap;

  g_return_if_fail (domain != NULL);  
  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);
  if (min == max)
    max = min + 1.;
  colormap = colormap_jet ();
  old_color = GTS_OBJECT_CLASS (s->vertex_class)->color;
  GTS_OBJECT_CLASS (s->vertex_class)->color = variable_color;
  color_data[0] = domain;
  color_data[1] = v;
  color_data[2] = colormap;
  color_data[3] = &min;
  color_data[4] = &max;
  gts_surface_write_oogl (s, fp);
  GTS_OBJECT_CLASS (s->vertex_class)->color = old_color;

  colormap_destroy (colormap);
}

/* GtsColoredVertex: Header */

typedef struct _GtsColoredVertex         GtsColoredVertex;

struct _GtsColoredVertex {
  /*< private >*/
  GtsVertex parent;

  /*< public >*/
  GtsColor c;
};

#define GTS_COLORED_VERTEX(obj)            GTS_OBJECT_CAST (obj,\
					         GtsColoredVertex,\
					         gts_colored_vertex_class ())
#define GTS_IS_COLORED_VERTEX(obj)         (gts_object_is_from_class (obj,\
						 gts_colored_vertex_class ()))

GtsVertexClass * gts_colored_vertex_class  (void);

/* GtsColoredVertex: Object */

static void gts_colored_vertex_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (c.r)");
    return;
  }
  GTS_COLORED_VERTEX (*o)->c.r = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (c.g)");
    return;
  }
  GTS_COLORED_VERTEX (*o)->c.g = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (c.b)");
    return;
  }
  GTS_COLORED_VERTEX (*o)->c.b = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gts_colored_vertex_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g %g %g", 
	   GTS_COLORED_VERTEX (o)->c.r, 
	   GTS_COLORED_VERTEX (o)->c.g, 
	   GTS_COLORED_VERTEX (o)->c.b);
}

static GtsColor gts_colored_vertex_color (GtsObject * o)
{
  return GTS_COLORED_VERTEX (o)->c;
}

static void gts_colored_vertex_class_init (GtsVertexClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gts_colored_vertex_read;
  GTS_OBJECT_CLASS (klass)->write = gts_colored_vertex_write;
  GTS_OBJECT_CLASS (klass)->color = gts_colored_vertex_color;
}

static void gts_colored_vertex_init (GtsColoredVertex * object)
{
  object->c.r = object->c.g = object->c.b = 1.;
}

GtsVertexClass * gts_colored_vertex_class (void)
{
  static GtsVertexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gts_colored_vertex_info = {
      "GtsColoredVertex",
      sizeof (GtsColoredVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) gts_colored_vertex_class_init,
      (GtsObjectInitFunc) gts_colored_vertex_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_vertex_class ()),
				  &gts_colored_vertex_info);
  }

  return klass;
}

/* GfsVertex: Header */

typedef struct _GfsVertex         GfsVertex;

struct _GfsVertex {
  /*< private >*/
  GtsColoredVertex parent;

  /*< public >*/
  gdouble v;
};

#define GFS_VERTEX(obj)            GTS_OBJECT_CAST (obj,\
					         GfsVertex,\
					         gfs_vertex_class ())
#define GFS_IS_VERTEX(obj)         (gts_object_is_from_class (obj,\
						 gfs_vertex_class ()))

GtsPointClass * gfs_vertex_class  (void);

/* GfsVertex: Object */

static void gfs_vertex_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v)");
    return;
  }
  GFS_VERTEX (*o)->v = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_vertex_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gts_colored_vertex_class ())->parent_class->write)
      (o, fp);
  fprintf (fp, " %g", GFS_VERTEX (o)->v);
}

static void gfs_vertex_class_init (GtsPointClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_vertex_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_vertex_write;
}

GtsPointClass * gfs_vertex_class (void)
{
  static GtsPointClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_vertex_info = {
      "GfsVertex",
      sizeof (GfsVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) gfs_vertex_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gts_colored_vertex_class ()),
			    &gfs_vertex_info);
  }

  return klass;
}

/* GfsTwistedVertex: Header */

typedef struct _GfsTwistedVertex         GfsTwistedVertex;

struct _GfsTwistedVertex {
  /*< private >*/
  GfsVertex parent;

  /*< public >*/
  gdouble theta;
};

#define GFS_TWISTED_VERTEX(obj)            GTS_OBJECT_CAST (obj,\
					         GfsTwistedVertex,\
					         gfs_twisted_vertex_class ())
#define GFS_IS_TWISTED_VERTEX(obj)         (gts_object_is_from_class (obj,\
						 gfs_twisted_vertex_class ()))

GtsPointClass * gfs_twisted_vertex_class  (void);

/* GfsTwistedVertex: Object */

static void gfs_twisted_vertex_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_twisted_vertex_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_twisted_vertex_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (theta)");
    return;
  }
  GFS_TWISTED_VERTEX (*o)->theta = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_twisted_vertex_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_twisted_vertex_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_twisted_vertex_class ())->parent_class->write)
      (o, fp);
  fprintf (fp, " %g", GFS_TWISTED_VERTEX (o)->theta);
}

static void gfs_twisted_vertex_class_init (GtsPointClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_twisted_vertex_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_twisted_vertex_write;
}

GtsPointClass * gfs_twisted_vertex_class (void)
{
  static GtsPointClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_twisted_vertex_info = {
      "GfsTwistedVertex",
      sizeof (GfsTwistedVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) gfs_twisted_vertex_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_vertex_class ()),
			    &gfs_twisted_vertex_info);
  }

  return klass;
}

static void matrix_transpose (GtsMatrix * m)
{
  guint i, j;

  for (i = 1; i < 3; i++)
    for (j = 0; j < i; j++) {
      gdouble t = m[i][j];

      m[i][j] = m[j][i];
      m[j][i] = t;
    }
}

static void base (GtsMatrix * b, GtsPoint * p1, GtsPoint * p2)
{
  GtsVector x, y;

  x[0] = b[0][0];
  x[1] = b[1][0];
  x[2] = b[2][0];
  gts_vector_init (b[2], p2, p1);
  gts_vector_normalize (b[2]);
  gts_vector_cross (y, b[2], x);
  if (gts_vector_norm (y) > 1e-2) {
    b[1][0] = y[0];
    b[1][1] = y[1];
    b[1][2] = y[2];
    gts_vector_normalize (b[1]);
  }
  gts_vector_cross (b[0], b[1], b[2]);
  gts_vector_normalize (b[0]);
  matrix_transpose (b);
}

static void edge_list (GtsMatrix * b, GtsPoint * o,
		       GSList * profile,
		       GtsSurface * s,
		       GtsEdge ** e, guint ne)
{
  guint i;
  GtsVertex * vold = NULL;
  GtsVertex * vfirst = NULL;
  gboolean colored = FALSE;
  GtsMatrix * c;

  if (GTS_IS_COLORED_VERTEX (o) && 
      gts_object_class_is_from_class (GTS_OBJECT_CLASS (s->vertex_class),
				      GTS_OBJECT_CLASS (gts_colored_vertex_class ())))
    colored = TRUE;
  if (GFS_IS_TWISTED_VERTEX (o)) {
    gdouble t = GFS_TWISTED_VERTEX (o)->theta;
    gdouble sint = sin (t), cost = cos (t);
    GtsMatrix * r = gts_matrix_new (cost, -sint, 0., 0.,
				    sint,  cost, 0., 0.,
				    0.,      0., 1., 0.,
				    0.,      0., 0., 0.);
    
    c = gts_matrix_product (b, r);
    gts_matrix_destroy (r);
  }
  else
    c = gts_matrix_new (b[0][0], b[0][1], b[0][2], 0.,
			b[1][0], b[1][1], b[1][2], 0.,
			b[2][0], b[2][1], b[2][2], 0.,
			0., 0., 0., 0.);

  for (i = 0; i <= ne && profile; i++, profile = profile->next) {
    GtsPoint * p = profile->data;
    GtsVertex * v = gts_vertex_new (s->vertex_class, p->x, p->y, 0.);

    if (colored)
      GTS_COLORED_VERTEX (v)->c = GTS_COLORED_VERTEX (o)->c;
    
    gts_point_transform (GTS_POINT (v), c);
    GTS_POINT (v)->x += o->x;
    GTS_POINT (v)->y += o->y;
    GTS_POINT (v)->z += o->z;

    if (vold != NULL)
      e[i-1] = gts_edge_new (s->edge_class, vold, v);
    vold = v;
    if (!vfirst) vfirst = v;
  }
  if (i <= ne && vold)
    e[i-1] = gts_edge_new (s->edge_class, vold, vfirst);
  gts_matrix_destroy (c);
}

static void add_face (GtsSurface * s, GtsEdge ** e1, GtsEdge ** e2,
		      guint ne)
{
  guint i;

  for (i = 0; i < ne; i++) {
    GtsVertex * v1 = GTS_SEGMENT (e1[i])->v1;
    GtsVertex * v2 = GTS_SEGMENT (e2[i])->v1;
    GtsVertex * v3 = GTS_SEGMENT (e2[i])->v2;
    GtsVertex * v4 = GTS_SEGMENT (e1[i])->v2;
    GtsEdge * e3 = gts_edge_new (s->edge_class, v1, v3);
    GtsEdge * e4 = GTS_EDGE (gts_vertices_are_connected (v1, v2));
    GtsEdge * e5 = GTS_EDGE (gts_vertices_are_connected (v3, v4));

    if (e4 == NULL)
      e4 = gts_edge_new (s->edge_class, v1, v2);
    if (e5 == NULL)
      e5 = gts_edge_new (s->edge_class, v3, v4);

    gts_surface_add_face (s, gts_face_new (s->face_class, e4, e2[i], e3));
    gts_surface_add_face (s, gts_face_new (s->face_class, e3, e5, e1[i]));
  }
}

static GList * next_far_enough (GList * p, gdouble size)
{
  GtsPoint * ps;
  GList * pf = NULL;

  if (p == NULL)
    return NULL;
  ps = p->data;
  p = p->next;
  size *= size;
  while (p && !pf) {
    if (gts_point_distance2 (ps, p->data) > size)
      pf = p;
    p = p->next;
  }
  return pf;
}

void gfs_extrude_profile (GtsSurface * s,
			  GSList * profile,
			  gboolean closed,
			  GList * path)
{
  GtsMatrix * r;
  GtsPoint * p0, * p1, * p2;
  GtsEdge ** e1, ** e2, ** tmp;
  GtsBBox * bb;
  gdouble size;
  guint ne;

  g_return_if_fail (s != NULL);
  g_return_if_fail (profile != NULL);
  g_return_if_fail (path != NULL);

  bb = gts_bbox_points (gts_bbox_class (), profile);
  size = bb->x2 - bb->x1;
  if (bb->y2 - bb->y1 > size)
    size = bb->y2 - bb->y1;
  gts_object_destroy (GTS_OBJECT (bb));

  size /= 4.;

  p0 = path->data;
  path = next_far_enough (path, size);
  if (path == NULL)
    return;
  p1 = path->data;
  r = gts_matrix_identity (NULL);
  ne = closed ? g_slist_length (profile) : g_slist_length (profile) - 1;
  e1 = g_malloc (sizeof (GtsEdge *)*ne);
  e2 = g_malloc (sizeof (GtsEdge *)*ne);

  base (r, p0, p1);
  edge_list (r, p0, profile, s, e1, ne);
  do {
    path = next_far_enough (path, size);
    p2 = path ? path->data : NULL;
    if (p2)
      base (r, p0, p2);
    else
      base (r, p0, p1);
    edge_list (r, p1, profile, s, e2, ne);
    add_face (s, e1, e2, ne);
    tmp = e1;
    e1 = e2;
    e2 = tmp;
    p0 = p1;
    p1 = p2;
  } while (p1);

  g_free (e1);
  g_free (e2);
  gts_matrix_destroy (r);
}

static gdouble triangle_area (FttVector p1, FttVector p2, FttVector p3)
{
  GtsVector v1, v2, a;

  v1[0] = p2.x - p1.x; v1[1] = p2.y - p1.y; v1[2] = p2.z - p1.z;
  v2[0] = p3.x - p2.x; v2[1] = p3.y - p2.y; v2[2] = p3.z - p2.z;
  gts_vector_cross (a, v1, v2);
  return gts_vector_norm (a)/2.;
}

static gdouble circumcircle_radius (FttVector p1, FttVector p2, FttVector p3)
{
  gdouble area = triangle_area (p1, p2, p3);

  if (area == 0.)
    return G_MAXDOUBLE;
  else {
    GtsVector a, b, c;
    gts_vector_init (a, &p1, &p2);
    gts_vector_init (b, &p2, &p3);
    gts_vector_init (c, &p3, &p1);
    return gts_vector_norm (a)*gts_vector_norm (b)*gts_vector_norm (c)/area;
  }
}

static GSList * circle_profile (GtsPointClass * klass, 
				gdouble radius, guint np)
{
  GSList * lp = NULL;
  guint i;

  for (i = 0; i < np; i++) {
    gdouble a = 2.*M_PI*i/(gdouble) np;

    lp = g_slist_prepend (lp, gts_point_new (klass, radius*cos (a), radius*sin(a), 0.));
  }
  return lp;
}

static GSList * ribbon_profile (GtsPointClass * klass, 
				gdouble half_width)
{
  GSList * lp = NULL;

  lp = g_slist_prepend (lp, gts_point_new (klass, 0., -half_width, 0.));
  lp = g_slist_prepend (lp, gts_point_new (klass, 0., half_width, 0.));
  return lp;
}

#if (!FTT_2D)
static void vorticity_vector (FttCell * cell, gpointer * data)
{
  gdouble size = ftt_cell_size (cell);
  GfsVariable ** g = data[0];
  GfsVariable ** v = data[1];

  GFS_VALUE (cell, g[0]) = (gfs_center_gradient (cell, FTT_Y, v[2]->i) -
			    gfs_center_gradient (cell, FTT_Z, v[1]->i))/size;
  GFS_VALUE (cell, g[1]) = (gfs_center_gradient (cell, FTT_Z, v[0]->i) -
			    gfs_center_gradient (cell, FTT_X, v[2]->i))/size;
  GFS_VALUE (cell, g[2]) = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
			    gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
}
#endif /* 3D */

static gdouble interpolated_velocity (FttCell * cell, FttVector p, GfsVariable ** U,
				      gdouble direction,
				      FttVector * u)
{
  FttComponent c;
  gdouble nu = 0.;
  gdouble (* interpolate) (FttCell *, FttVector, GfsVariable * v) = GFS_IS_MIXED (cell) ?
    gfs_mixed_cell_interpolate : gfs_interpolate;
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&u->x)[c] = direction* (*interpolate) (cell, p, U[c]);
    nu += (&u->x)[c]*(&u->x)[c];
  }
  return nu;
}

static GList * grow_curve (GfsDomain * domain,
			   GfsVariable ** U,
			   FttVector p,
			   GfsVariable * var,
			   gdouble min, 
			   gdouble max,
			   gboolean twist,
			   GList * path,
			   gdouble direction,
			   gboolean (* stop) (FttCell *, GList *, gpointer),
			   gpointer data)
{
  FttCell * cell;
  gdouble delta = 0.2;
  GtsPoint * oldp = NULL;
  FttVector p1, p2;
  gdouble cost = 0., theta = 0.;
  gdouble maxcost = 4e-9;
  guint nstep = 0, nmax = 10000;
  GtsPointClass * path_class = gfs_vertex_class ();
  Colormap * colormap = NULL;

  if (min < max)
    colormap = colormap_jet ();

#if (!FTT_2D)
  GfsVariable * vort[FTT_DIMENSION];
  if (twist) {
    FttComponent c;
    gpointer data[2];

    path_class = GTS_POINT_CLASS (gfs_twisted_vertex_class ());
    for (c = 0; c < FTT_DIMENSION; c++)
      vort[c] = gfs_temporary_variable (domain);
    gfs_variable_set_vector (vort, FTT_DIMENSION);
    data[0] = vort;
    data[1] = U;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) vorticity_vector, data);
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_domain_cell_traverse (domain,
				FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) vort[c]->fine_coarse, vort[c]);
  }
#else /* 2D */
  twist = FALSE;
#endif /* 2D */  

  p1 = p2 = p;
  while ((cell = gfs_domain_locate (domain, p, -1, NULL)) != NULL &&
	 circumcircle_radius (p1, p2, p) > ftt_cell_size (cell) &&
	 nmax--) {
    gdouble h = delta*ftt_cell_size (cell);
    FttVector u;
    FttComponent c;
    gdouble nu = 0.;

    cost += triangle_area (p1, p2, p);
    p1 = p2;
    p2 = p;
    if (oldp == NULL || cost > maxcost) {
      oldp = gts_point_new (path_class, p.x, p.y, p.z);
      if (var)
	GFS_VERTEX (oldp)->v = gfs_interpolate (cell, p, var);
      if (colormap)
	GTS_COLORED_VERTEX (oldp)->c = 
	  colormap_color (colormap, (GFS_VERTEX (oldp)->v - min)/(max - min));
      if (twist)
	GFS_TWISTED_VERTEX (oldp)->theta = theta;
      path = g_list_prepend (path, oldp);
      if (stop != NULL && (* stop) (cell, path, data))
	break;
      cost = 0.;
      nstep = 0;
    }

    nu = interpolated_velocity (cell, p, U, direction, &u);
    if (nu > 0) {
      FttVector p1 = p;
      FttCell * cell1;

      nu = 2.*sqrt (nu);
      for (c = 0; c < FTT_DIMENSION; c++)
	(&p1.x)[c] += h*(&u.x)[c]/nu;
      cell1 = gfs_domain_locate (domain, p1, -1, NULL);
      if (!cell1)
	break;
      nu = interpolated_velocity (cell1, p1, U, direction, &u);
    }
    else
      break;
    if (nu > 0. && nstep++ < nmax) {
      FttVector p1;

      p1 = p;
      nu = sqrt (nu);
      for (c = 0; c < FTT_DIMENSION; c++)
	((gdouble *) &p)[c] += h*((gdouble *) &u)[c]/nu;
#if (!FTT_2D)
      if (twist) {
	GtsVector rot;
	GtsVector dx;

	dx[0] = p1.x - p.x; dx[1] = p1.y - p.y; dx[2] = p1.z - p.z;
	for (c = 0; c < FTT_DIMENSION; c++)
	  rot[c] = gfs_interpolate (cell, p1, vort[c]);
	theta += gts_vector_scalar (rot, dx)/nu;
      }
#endif /* 3D */
    }
    else
      break;
  }
  if (oldp && (p2.x != oldp->x || p2.y != oldp->y || p2.z != oldp->z)) {
    cell = gfs_domain_locate (domain, p2, -1, NULL);
    if (cell) {
      oldp = gts_point_new (path_class, p2.x, p2.y, p2.z);
      if (var)
	GFS_VERTEX (oldp)->v = gfs_interpolate (cell, p2, var);
      if (twist)
	GFS_TWISTED_VERTEX (oldp)->theta = theta;
      path = g_list_prepend (path, oldp);
    }
  }

  if (colormap)
    colormap_destroy (colormap);

#if (!FTT_2D)
  if (twist) {
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      gts_object_destroy (GTS_OBJECT (vort[c]));
  }
#endif /* 3D */
  return direction > 0. ? g_list_reverse (path) : path;
}

GList * gfs_streamline_new (GfsDomain * domain,
			    GfsVariable ** U,
			    FttVector p,
			    GfsVariable * var,
			    gdouble min,
			    gdouble max,
			    gboolean twist,
			    gboolean (* stop) (FttCell *, 
					       GList *,
					       gpointer),
			    gpointer data)
{
  GList * i, * path;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (U != NULL, NULL);

  i = grow_curve (domain, U, p, var, min, max, twist, NULL, 1., stop, data);
  path = g_list_remove_link (i, i);
  if (i != NULL)
    gts_object_destroy (i->data);
  g_list_free_1 (i);
  path = grow_curve (domain, U, p, var, min, max, twist, path, -1., stop, data);
  return path;
}

void gfs_streamline_write (GList * stream, FILE * fp)
{
  g_return_if_fail (fp != NULL);

  fprintf (fp, "GfsStreamline %u\n", g_list_length (stream));
  while (stream) {
    (* GTS_OBJECT (stream->data)->klass->write) (stream->data, fp);
    fputc ('\n', fp);
    stream = stream->next;
  }
}

GList * gfs_streamline_read (GtsFile * fp)
{
  GList * stream = NULL;
  guint n = 0, nv;

  g_return_val_if_fail (fp != NULL, NULL);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsStreamline)");
    return NULL;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (number of vertices)");
    return NULL;
  }
  nv = atoi (fp->token->str);
  gts_file_first_token_after (fp, '\n');
  while (fp->type != GTS_ERROR && n < nv) {
    GtsObject * o = 
      gts_object_new (GTS_OBJECT_CLASS (gfs_twisted_vertex_class ()));

    (*o->klass->read) (&o, fp);
    gts_file_first_token_after (fp, '\n');
    stream = g_list_prepend (stream, o);
    n++;
  }

  if (fp->type == GTS_ERROR) {
    g_list_free (stream);
    return NULL;
  }

  return stream;
}

void gfs_streamline_draw (GList * stream, FILE * fp)
{
  guint n = g_list_length (stream);

  g_return_if_fail (fp != NULL);

  fprintf (fp, "VECT 1 %u 0 %u 0\n", n, n);
  while (stream) {
    fprintf (fp, "%g %g %g\n",
	     GTS_POINT (stream->data)->x,
	     GTS_POINT (stream->data)->y,
	     GTS_POINT (stream->data)->z);
    stream = stream->next;
  }
}

void gfs_streamline_destroy (GList * stream)
{
  g_list_foreach (stream, (GFunc) gts_object_destroy, NULL);
  g_list_free (stream);
}

void gfs_draw_stream_cylinder (GfsDomain * domain,
			       FttVector p,
			       gdouble radius,
			       GfsVariable * var,
			       gdouble min, gdouble max,
			       FILE * fp)
{
  GSList * profile;
  GList * path;
  GtsSurface * s;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  s = gts_surface_new (gts_surface_class (),
		       gts_face_class (),
		       gts_edge_class (),
		       min < max ? gts_colored_vertex_class () :
		       gts_vertex_class ());
  path = gfs_streamline_new (domain, gfs_domain_velocity (domain), p, var, min, max, FALSE, 
			     NULL, NULL);
  profile = circle_profile (gts_point_class (), radius, 10);
  gfs_extrude_profile (s, profile, TRUE, path);
  gts_surface_write_oogl (s, fp);
  gts_object_destroy (GTS_OBJECT (s));
  gfs_streamline_destroy (path);
  g_slist_foreach (profile, (GFunc) gts_object_destroy, NULL);
  g_slist_free (profile);
}

void gfs_draw_stream_ribbon (GfsDomain * domain,
			     FttVector p,
			     gdouble half_width,
			     GfsVariable * var,
			     gdouble min, gdouble max,
			     FILE * fp)
{
  GList * path;
  GSList * profile;
  GtsSurface * s;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  s = gts_surface_new (gts_surface_class (),
		       gts_face_class (),
		       gts_edge_class (),
		       min < max ? gts_colored_vertex_class () :
		       gts_vertex_class ());
  path = gfs_streamline_new (domain, gfs_domain_velocity (domain), p, var, min, max, TRUE, 
			     NULL, NULL);
  profile = ribbon_profile (gts_point_class (), half_width);
  gfs_extrude_profile (s, profile, FALSE, path);
  gts_surface_write_oogl (s, fp);
  gts_object_destroy (GTS_OBJECT (s));
  gfs_streamline_destroy (path);
  g_slist_foreach (profile, (GFunc) gts_object_destroy, NULL);
  g_slist_free (profile);
}

void gfs_draw_streamline (GfsDomain * domain,
			  FttVector p,
			  FILE * fp)
{
  GList * path;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  path = gfs_streamline_new (domain, gfs_domain_velocity (domain), p, NULL, 0., 0., FALSE, 
			     NULL, NULL);
  gfs_streamline_draw (path, fp);
  gfs_streamline_destroy (path);
}

#if !FTT_2D

static gdouble point_orientation (FttVector p[3], FttVector * c)
{
  gdouble adx, bdx, cdx;
  gdouble ady, bdy, cdy;
  gdouble adz, bdz, cdz;
  
  adx = p[0].x - c->x;
  bdx = p[1].x - c->x;
  cdx = p[2].x - c->x;
  ady = p[0].y - c->y;
  bdy = p[1].y - c->y;
  cdy = p[2].y - c->y;
  adz = p[0].z - c->z;
  bdz = p[1].z - c->z;
  cdz = p[2].z - c->z;
  
  return (adx * (bdy * cdz - bdz * cdy) +
	  bdx * (cdy * adz - cdz * ady) +
	  cdx * (ady * bdz - adz * bdy));
}

/**
 * gfs_plane_cuts_cell:
 * @plane: three points belonging to the plane.
 * @cell: a #FttCell.
 *
 * Returns: %TRUE if @plane cuts @cell, %FALSE otherwise.
 */
gboolean gfs_plane_cuts_cell (FttVector plane[3], FttCell * cell)
{
  FttVector o;
  gdouble h = ftt_cell_size (cell)*SLIGHTLY_LARGER;
  guint i;

  g_return_val_if_fail (cell != NULL, FALSE);

  ftt_cell_pos (cell, &o);
  o.x -= h/2.; o.y -= h/2.; o.z -= h/2.;
  for (i = 0; i < 12; i++) {
    FttVector e, d;
    gdouble a, b;
    d.x = o.x + h*edge[i][0].x; d.y = o.y + h*edge[i][0].y; d.z = o.z + h*edge[i][0].z;
    e.x = o.x + h*edge[i][1].x; e.y = o.y + h*edge[i][1].y; e.z = o.z + h*edge[i][1].z;
    a = point_orientation (plane, &e);
    b = point_orientation (plane, &d);
    if ((a <= 0. && b > 0.) || (a >= 0. && b < 0.))
      return TRUE;
  }
  return FALSE;
}

static void cube_plane_intersection (FttCell * cell,
				     FttVector * O,
				     FttVector * n,
				     FttVector p[12],
				     gint orient[12],
				     GfsVariable * var,
				     gdouble v[12],
				     gint max_level)
{
  FttVector o;
  gdouble h = ftt_cell_size (cell)*SLIGHTLY_LARGER, vc[8];
  guint i;

  if (var)
    for (i = 0; i < 8; i++)
      vc[i] = G_MAXDOUBLE;

  ftt_cell_pos (cell, &o);
  o.x -= h/2.; o.y -= h/2.; o.z -= h/2.;
  for (i = 0; i < 12; i++) {
    FttVector e, d;
    d.x = o.x + h*edge[i][0].x; d.y = o.y + h*edge[i][0].y; d.z = o.z + h*edge[i][0].z;
    e.x = o.x + h*edge[i][1].x; e.y = o.y + h*edge[i][1].y; e.z = o.z + h*edge[i][1].z;
    gdouble den = n->x*(e.x - d.x) + n->y*(e.y - d.y) + n->z*(e.z - d.z);
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      gdouble t = (n->x*(O->x - d.x) + n->y*(O->y - d.y) + n->z*(O->z - d.z))/den;
      if (t >= 0. && t < 1.) {
	p[i].x = d.x + t*(e.x - d.x); p[i].y = d.y + t*(e.y - d.y); p[i].z = d.z + t*(e.z - d.z);
	orient[i] = (n->x*(e.x - O->x) + n->y*(e.y - O->y) + n->z*(e.z - O->z) > 0.);
	if (var) {
	  guint j = edge1[i][0];
	  if (vc[j] == G_MAXDOUBLE)
	    vc[j] = gfs_cell_corner_value (cell, corner[j], var, max_level);
	  d.z = vc[j];
	  j = edge1[i][1];
	  if (vc[j] == G_MAXDOUBLE)
	    vc[j] = gfs_cell_corner_value (cell, corner[j], var, max_level);
	  e.z = vc[j];
	  v[i] = d.z + t*(e.z - d.z);
	}
      }
    }
  }
}

/**
 * gfs_cut_cube_vertices:
 * @cell: a #FttCell.
 * @maxlevel: the maximum level to consider (or -1).
 * @p: a point on the plane.
 * @n: the normal to the plane.
 * @v: where to return the vertices coordinates.
 * @d: where to return the direction.
 * @var: a #GfsVariable or %NULL.
 * @val: where to return the values of @var or %NULL.
 *
 * Fills @v, @d and @val with the coordinates/values of the vertices,
 * intersections of @cell with the plane defined by @p and @n.
 *
 * The vertices are ordered consistently to define a consistent,
 * oriented polygon.
 *
 * Returns: the number of vertices (0 if the plane does not cut the cell).
 */
guint gfs_cut_cube_vertices (FttCell * cell, gint maxlevel,
			     FttVector * p, FttVector * n,
			     FttVector v[12], FttDirection d[12],
			     GfsVariable * var,
			     gdouble val[12])
{
  FttVector a[12];
  gdouble vv[12];
  gint orient[12];
  guint i;

  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (p != NULL, 0);
  g_return_val_if_fail (n != NULL, 0);
  g_return_val_if_fail ((var == NULL && val == NULL) || (var != NULL && val != NULL), 0);

  cube_plane_intersection (cell, p, n, a, orient, var, vv, maxlevel);
  for (i = 0; i < 12; i++) {
    guint nv = 0, e = i;
    while (orient[e] >= 0) {
      guint m = 0, * ne = connect[e][orient[e]];
      d[nv] = ne[3];
      if (var)
	val[nv] = vv[e];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
	e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}

#endif /* 3D */
