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
 * \brief Surfaces.
 */

#include "simulation.h"
#include "surface.h"

/**
 * Abstract class for surfaces (an oriented surface (in 3D) or an
 * oriented curve (in 2D)).
 * \beginobject{GfsGenericSurface}
 */

GfsGenericSurfaceClass * gfs_generic_surface_class (void)
{
  static GfsGenericSurfaceClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_generic_surface_info = {
      "GfsGenericSurface",
      sizeof (GtsObject),
      sizeof (GfsGenericSurfaceClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_generic_surface_info);
  }

  return klass;
}

/**
 * gfs_cell_is_cut:
 * @cell: a #FttCell.
 * @s: a #GfsGenericSurface.
 * @flatten: if set to %TRUE, @cell is flattened in the z direction.
 * @maxlevel: the maximum (virtual) cell level to consider.
 *
 * Returns: a (possibly new) #GfsGenericSurface containing a subset of @s which may
 * intersect @cell or %NULL if @s does not intersect @cell.
 */
GfsGenericSurface * gfs_cell_is_cut (FttCell * cell, GfsGenericSurface * s, 
				     gboolean flatten, gint maxlevel)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (s != NULL, NULL);
  
  g_assert (GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->cell_is_cut);
  return (* GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->cell_is_cut) 
    (cell, s, flatten, maxlevel);
}

static void cell_traverse_cut (FttCell * cell,
			       GfsGenericSurface * s,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       FttCellTraverseCutFunc func,
			       gpointer data,
			       gboolean flatten)
{
  GfsGenericSurface * s1 = gfs_cell_is_cut (cell, s, flatten && FTT_CELL_IS_LEAF (cell), -1);

  if (s1 == NULL)
    return;
  if (order == FTT_PRE_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, s1, data);
  if (!FTT_CELL_IS_LEAF (cell)) {
    struct _FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if ((flags & FTT_TRAVERSE_DESTROYED) != 0 || !FTT_CELL_IS_DESTROYED (c))
	cell_traverse_cut (c, s1, order, flags, func, data, flatten);
    }
  }
  if (order == FTT_POST_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, s1, data);
  if (s1 != s) {
    if (GFS_IS_SURFACE (s1))
      GFS_SURFACE (s1)->bbtree = NULL;
    gts_object_destroy (GTS_OBJECT (s1));
  }
}

/**
 * gfs_cell_traverse_cut:
 * @root: the root #FttCell of the tree to traverse.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell cut by @s.
 */
void gfs_cell_traverse_cut (FttCell * root,
			    GfsGenericSurface * s,
			    FttTraverseType order,
			    FttTraverseFlags flags,
			    FttCellTraverseCutFunc func,
			    gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_cut (root, s, order, flags, func, data, FALSE);
}

/**
 * gfs_cell_traverse_cut_2D:
 * @root: the root #FttCell of the tree to traverse.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell cut by @s.
 *
 * The cells are "flattened" in the z-direction.
 */
void gfs_cell_traverse_cut_2D (FttCell * root,
			       GfsGenericSurface * s,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       FttCellTraverseCutFunc func,
			       gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_cut (root, s, order, flags, func, data, TRUE);
}

/**
 * gfs_generic_surface_read:
 * @s: a #GfsGenericSurface.
 * @sim: a #GfsSimulation.
 * @fp: a #GtsFile.
 * 
 * Calls the read() method of @s.
 */
void gfs_generic_surface_read (GfsGenericSurface * s, gpointer sim, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) s;

  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

  o->reserved = sim;
  (* GTS_OBJECT (s)->klass->read) (&o, fp);
}

/**
 * gfs_generic_surface_write:
 * @s: a #GfsGenericSurface.
 * @sim: a #GfsSimulation.
 * @fp: a file pointer.
 * 
 * Calls the write() method of @s.
 */
void gfs_generic_surface_write (GfsGenericSurface * s, gpointer sim, FILE * fp)
{
  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

  GTS_OBJECT (s)->reserved = sim;
  (* GTS_OBJECT (s)->klass->write) (GTS_OBJECT (s), fp);
}

/**
 * gfs_surface_segment_intersection:
 * @s: a #GfsGenericSurface.
 * @cell: a #FttCell containing @I.
 * @I: a GfsSegment.
 *
 * Fills @I with the intersection of @s and @I.
 *
 * Returns: the number of times @s intersects @I.
 */
guint gfs_surface_segment_intersection (GfsGenericSurface * s,
					FttCell * cell,
					GfsSegment * I)
{
  g_return_val_if_fail (s != NULL, 0);
  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (I != NULL, 0);

  g_assert (GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->segment_intersection);
  return (* GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->segment_intersection) (s, cell, I);
}

/**
 * gfs_surface_segment_normal:
 * @s: a #GfsGenericSurface.
 * @cell: a #FttCell containing @I.
 * @I: a GfsSegment.
 * @n: a #GtsVector.
 *
 * Fills @n with the normal to @s at the intersection of @s and @I.
 */
void gfs_surface_segment_normal (GfsGenericSurface * s,
				 FttCell * cell,
				 GfsSegment * I,
				 GtsVector n)
{
  g_return_if_fail (s != NULL);
  g_return_if_fail (cell != NULL);
  g_return_if_fail (I != NULL);
  g_return_if_fail (I->n > 0);
  g_return_if_fail (n != NULL);

  g_assert (GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->segment_normal);
  (* GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->segment_normal) (s, cell, I, n);
}

/**
 * gfs_surface_point_is_inside:
 * @s: a #GfsGenericSurface.
 * @p: a #FttVector.
 *
 * Returns: 1 if @p is inside @s, 0 if @p lies on the boundary of @s,
 * -1 otherwise.
 */
gint gfs_surface_point_is_inside (GfsGenericSurface * s,
				  FttVector * p)
{
  g_return_val_if_fail (s != NULL, 0);
  g_return_val_if_fail (p != NULL, 0);

  g_assert (GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->point_is_inside);
  return (* GFS_GENERIC_SURFACE_CLASS (GTS_OBJECT (s)->klass)->point_is_inside) (s, p);
}

/** \endobject{GfsGenericSurface} */

/**
 * A GfsSurface is an oriented surface (in 3D) or an oriented curve (in 2D).
 * \beginobject{GfsSurface}
 */

/**
 * gfs_surface_transformation:
 * @surface: a #GtsSurface
 * @rotate: a #GtsVector
 * @translate: a #GtsVector 
 * @scale: a #GtsVector
 * @flip: a #gboolean
 * @matrix: a pointer on a #GtsMatrix
 * 
 * Translates, rotates, scales and flips @surface.
 */
void gfs_surface_transformation (GtsSurface * surface, 
				 GtsVector rotate,
				 GtsVector translate, 
				 GtsVector scale,
				 gboolean flip,
				 GtsMatrix ** matrix)
{
  GtsMatrix * m;
  gint c;

  g_return_if_fail (matrix != NULL);
  
  m = gts_matrix_translate (NULL, translate);    
  for (c = 2; c >= 0; c--)
    if (rotate[c] != 0.) {
      GtsVector r = {0.,0.,0.};
      r[c] = 1.;
      GtsMatrix * mr = gts_matrix_rotate (NULL, r, rotate[c]*M_PI/180.);
      GtsMatrix * m1 = gts_matrix_product (m, mr);
      gts_matrix_destroy (m);
      gts_matrix_destroy (mr);
      m = m1;
    }

  GtsMatrix * ms = gts_matrix_scale (NULL, scale);
  if (*matrix)
    gts_matrix_destroy (*matrix);
  *matrix = gts_matrix_product (m, ms);
  gts_matrix_destroy (m);
  gts_matrix_destroy (ms);
  
  if (surface) {
    gts_surface_foreach_vertex (surface, (GtsFunc) gts_point_transform, *matrix);
    gts_matrix_destroy (*matrix);
    *matrix = NULL;
    if (flip)
      gts_surface_foreach_face (surface, (GtsFunc) gts_triangle_revert, NULL);
  }
  else {
    GtsMatrix * i = gts_matrix_inverse (*matrix);
    gts_matrix_destroy (*matrix);
    *matrix = i;
  }
}

static void point_map (GtsPoint * p, GfsSimulation * sim)
{
  gfs_simulation_map (sim, (FttVector *) &p->x);
}

static void surface_read (GtsObject ** o, GtsFile * fp)
{
  GfsSurface * surface = GFS_SURFACE (*o);
  gboolean dimensional = FALSE, implicit = FALSE;

  if (fp->type == '(') { /* implicit surface */
    gts_file_next_token (fp);
    if (surface->f)
      gts_object_destroy (GTS_OBJECT (surface->f));
    surface->f = gfs_function_new (gfs_function_spatial_class (), 0.);
    gfs_function_read (surface->f, gfs_object_simulation (*o), fp);
    if (fp->type == GTS_ERROR)
      return;
    if (fp->type != ')') {
      gts_file_error (fp, "expecting a closing bracket");
      return;
    }
  }
  else if (fp->type == '{') { /* embedded surface */
    fp->scope_max++;
    gts_file_next_token (fp);
    if (surface->s)
      gts_object_destroy (GTS_OBJECT (surface->s));
    surface->s = gts_surface_new (gts_surface_class (), 
				  surface->face_class, 
				  surface->edge_class, 
				  surface->vertex_class);
    if (gts_surface_read (surface->s, fp))
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
  }
  else { /* surface file name */
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (filename)");
      return;
    }
    FILE * fptr = fopen (fp->token->str, "rt");
    if (fptr == NULL) {
      gts_file_error (fp, "cannot open file `%s'", fp->token->str);
      return;
    }
    GtsFile * fp1 = gts_file_new (fptr);
    surface->s = gts_surface_new (gts_surface_class (), 
				  surface->face_class, 
				  surface->edge_class, 
				  surface->vertex_class);
    if (gts_surface_read (surface->s, fp1)) {
      gts_file_error (fp, 
		      "file `%s' is not a valid GTS file\n"
		      "%s:%d:%d: %s",
		      fp->token->str, fp->token->str,
		      fp1->line, fp1->pos, fp1->error);
      gts_file_destroy (fp1);
      fclose (fptr);
      return;
    }
    gts_file_destroy (fp1);
    fclose (fptr);
    dimensional = TRUE;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    gdouble scale = 1.;
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "rx",       TRUE, &surface->rotate[0]},
      {GTS_DOUBLE, "ry",       TRUE, &surface->rotate[1]},
      {GTS_DOUBLE, "rz",       TRUE, &surface->rotate[2]},
      {GTS_DOUBLE, "sx",       TRUE, &surface->scale[0]},
      {GTS_DOUBLE, "sy",       TRUE, &surface->scale[1]},
      {GTS_DOUBLE, "sz",       TRUE, &surface->scale[2]},
      {GTS_DOUBLE, "tx",       TRUE, &surface->translate[0]},
      {GTS_DOUBLE, "ty",       TRUE, &surface->translate[1]},
      {GTS_DOUBLE, "tz",       TRUE, &surface->translate[2]},
      {GTS_DOUBLE, "scale",    TRUE, &scale},
      {GTS_INT,    "flip",     TRUE, &surface->flip},
      {GTS_INT,    "twod",     TRUE, &surface->twod},
      {GTS_INT,    "implicit", TRUE, &implicit},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;

    gboolean set = FALSE;
    guint i;
    for (i = 0; i < 11 && !set; i++)
      if (var[i].set)
	set = TRUE;

    if (set) {
      if (var[9].set) {
	surface->scale[0] *= scale;
	surface->scale[1] *= scale;
	surface->scale[2] *= scale;
      }
      GtsVector translate;
      gdouble L = gfs_object_simulation (surface)->physical_params.L;
      for (i = 0; i < FTT_DIMENSION; i++)
	translate[i] = surface->translate[i]/L;
      gfs_surface_transformation (surface->s, 
				  surface->rotate, translate, surface->scale,
				  surface->flip, 
				  &surface->m);
    }
  }

  if (dimensional) {
    g_assert (surface->s);
    gts_surface_foreach_vertex (surface->s, (GtsFunc) point_map, gfs_object_simulation (*o));
  }

  if (surface->s) {
    if (implicit)
      surface->bbtree = gts_bb_tree_surface (surface->s);
    else if (!gts_surface_is_orientable (surface->s))
	gts_file_error (fp, "surface is not orientable");
  }
}

static void surface_write (GtsObject * o, FILE * fp)
{
  GfsSurface * surface = GFS_SURFACE (o);
  if (surface->s) {
    fputs (" { ", fp);
    GtsSurface * s = surface->s;
    if (GFS_DOMAIN (gfs_object_simulation (o))->binary) {
      gboolean binary = GTS_POINT_CLASS (s->vertex_class)->binary;
      GTS_POINT_CLASS (s->vertex_class)->binary = TRUE;
      gts_surface_write (s, fp);
      GTS_POINT_CLASS (s->vertex_class)->binary = binary;
    }
    else
      gts_surface_write (s, fp);
    fputc ('}', fp);
  }
  else if (surface->f) {
    fputs (" (", fp);
    gfs_function_write (surface->f, fp);
    fputs (" )", fp);
  }
  if (surface->m) {
    fputs (" {\n", fp);
    if (gts_vector_norm (surface->translate) > 0.)
      fprintf (fp, "  tx = %g ty = %g tz = %g\n",
	       surface->translate[0], surface->translate[1], surface->translate[2]);
    if (surface->scale[0] != 1. || surface->scale[1] != 1. || surface->scale[2] != 1.)
      fprintf (fp, "  sx = %g sy = %g sz = %g\n",
	       surface->scale[0], surface->scale[1], surface->scale[2]);
    if (surface->rotate[0] != 0.)
      fprintf (fp, "  rx = %g\n", surface->rotate[0]);
    if (surface->rotate[1] != 0.)
      fprintf (fp, "  ry = %g\n", surface->rotate[1]);
    if (surface->rotate[2] != 0.)
      fprintf (fp, "  rz = %g\n", surface->rotate[2]);
    if (surface->flip)
      fputs ("  flip = 1\n", fp);
    if (surface->twod)
      fputs ("  twod = 1\n", fp);
    if (surface->bbtree)
      fputs ("  implicit = 1\n", fp);
    fputc ('}', fp);
  }
  else if (surface->bbtree)
    fputs (" { implicit = 1 }", fp);
  else
    fputs (" {}", fp);
}

static void surface_destroy (GtsObject * object)
{
  GfsSurface * s = GFS_SURFACE (object);
  if (s->s)
    gts_object_destroy (GTS_OBJECT (s->s));
  if (s->f)
    gts_object_destroy (GTS_OBJECT (s->f));
  if (s->m)
    gts_matrix_destroy (s->m);
  if (s->bbtree)
    gts_bb_tree_destroy (s->bbtree, TRUE);

  (* GTS_OBJECT_CLASS (gfs_surface_class ())->parent_class->destroy) (object);
}


static void face_overlaps_box (GtsTriangle * t, gpointer * data)
{
  GtsBBox * bb = data[0];
  GtsSurface ** s1 = data[1];

  if (gts_bbox_overlaps_triangle (bb, t)) {
    if (*s1 == NULL)
      *s1 = gts_surface_new (gts_surface_class (),
			     gts_face_class (),
			     gts_edge_class (),
			     gts_vertex_class ());
    gts_surface_add_face (*s1, GTS_FACE (t));
  }
}

#define SIGN(v) ((v) > 0. ? 1 : (v) < 0. ? -1 : 0)

static GfsGenericSurface * cell_is_cut (FttCell * cell, GfsGenericSurface * s1, 
					gboolean flatten, gint maxlevel)
{
  GfsSurface * s = GFS_SURFACE (s1);
  if (s->f) {
    if (!FTT_CELL_IS_LEAF (cell))
      return GFS_GENERIC_SURFACE (s);
    FttVector p;
    gdouble h = ftt_cell_size (cell)/2.;
    ftt_cell_pos (cell, &p);
    gint i, j, k = 0, sign = 0, n = 1;
    i = maxlevel - ftt_cell_level (cell);
    while (i-- > 0)
      n *= 2;
#if !FTT_2D
    for (k = - n; k <= n; k += 2)
#endif
      for (i = - n; i <= n; i += 2)
	for (j = - n; j <= n; j += 2) {
	  GtsPoint o;
	  o.x = p.x + i*h/n; o.y = p.y + j*h/n; o.z = p.z + k*h/n;
	  gdouble v = gfs_surface_implicit_value (s, o);
	  if (sign && sign*SIGN(v) <= 0)
	    return GFS_GENERIC_SURFACE (s);
	  sign = SIGN(v);
	}
    return NULL;
  }
  else if (s->s) {
    GtsSurface * s1 = NULL;
    gpointer data[2];
    GtsBBox bb;
    ftt_cell_bbox (cell, &bb);
    if (flatten)
      bb.z1 = bb.z2 = 0.;
    data[0] = &bb;
    data[1] = &s1;
    gts_surface_foreach_face (s->s, (GtsFunc) face_overlaps_box, data);
    if (s1 == NULL)
      return NULL;
    GfsSurface * s2 = GFS_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
    s2->s = s1;
    s2->bbtree = s->bbtree;
    return GFS_GENERIC_SURFACE (s2);
  }
  g_assert_not_reached ();
  return NULL;
}

static gdouble segment_triangle_intersection (GtsPoint * E, GtsPoint * D,
					      GtsTriangle * t,
					      gboolean * inside)
{
  GtsVertex * vA, * vB, * vC;
  GtsPoint * A, * B, * C;
  gint ABCE, ABCD, ADCE, ABDE, BCDE;
  GtsEdge * AB, * BC, * CA;
  gdouble a, b;
  gboolean reversed = FALSE;

  gts_triangle_vertices_edges (t, NULL, &vA, &vB, &vC, &AB, &BC, &CA);
  A = GTS_POINT (vA);
  B = GTS_POINT (vB);
  C = GTS_POINT (vC);
  ABCE = gts_point_orientation_3d_sos (A, B, C, E);
  ABCD = gts_point_orientation_3d_sos (A, B, C, D);
  if (ABCE < 0 || ABCD > 0) {
    GtsPoint * tmpp;
    gint tmp;

    tmpp = E; E = D; D = tmpp;
    tmp = ABCE; ABCE = ABCD; ABCD = tmp;
    reversed = TRUE;
  }
  if (ABCE < 0 || ABCD > 0)
    return -1.;
  ADCE = gts_point_orientation_3d_sos (A, D, C, E);
  if (ADCE < 0)
    return -1.;
  ABDE = gts_point_orientation_3d_sos (A, B, D, E);
  if (ABDE < 0)
    return -1.;
  BCDE = gts_point_orientation_3d_sos (B, C, D, E);
  if (BCDE < 0)
    return -1.;
  *inside = reversed ? (ABCD < 0) : (ABCE < 0);
  a = gts_point_orientation_3d (A, B, C, E);
  b = gts_point_orientation_3d (A, B, C, D);
  if (a != b)
    return reversed ? 1. - a/(a - b) : a/(a - b);
  /* D and E are contained within ABC */
  g_assert (a == 0.);
  return 0.5;
}

static void triangle_face_intersection (GtsTriangle * t, GfsSegment * I)
{
  gboolean ins;
  gdouble x = segment_triangle_intersection (I->E, I->D, t, &ins);
  
  if (x != -1.) {
    I->x += x; I->n++;
    I->inside += ins ? 1 : -1;
  }
}

static GtsPoint segment_intersection (GfsSegment * I)
{
  GtsPoint p;
  p.x = I->E->x + I->x*(I->D->x - I->E->x);
  p.y = I->E->y + I->x*(I->D->y - I->E->y);
  p.z = I->E->z + I->x*(I->D->z - I->E->z);
  /* lines below just to prevent compiler warnings about uninitialised fields */
  p.object.flags = 0;
  p.object.reserved = NULL;
  p.object.klass = NULL;
  return p;
}

static guint surface_segment_intersection (GfsGenericSurface * s1,
					   FttCell * cell,
					   GfsSegment * I)
{
  GfsSurface * s = GFS_SURFACE (s1);

  I->n = 0;
  I->x = 0.;
  I->inside = 0;

  if (s->f || s->bbtree) {
    gdouble vE = gfs_surface_implicit_value (s, *I->E);
    gdouble vD = gfs_surface_implicit_value (s, *I->D);
    
    if ((vE > 0. && vD <= 0.) || (vE <= 0. && vD > 0.)) {
      I->n = 1;
      I->inside = vE > 0. ? -1 : 1;

      /* secant-bisection root-finding */
      gdouble t, t1, t2, v1, v2;
      if (vE > vD) {
	v1 = vD; t1 = 1.;
	v2 = vE; t2 = 0.;
      }
      else {
	v1 = vE; t1 = 0.;
	v2 = vD; t2 = 1.;
      }
      I->x = (v1*t2 - v2*t1)/(v1 - v2);
      guint n = 0;
      gint side = 0;
      do {
	t = I->x;
	gdouble v = gfs_surface_implicit_value (s, segment_intersection (I));
	if (v < 0.) {
	  v1 = v; t1 = I->x;
	  if (side == -1) v2 /= 2.;
	  side = -1;
	}
	else {
	  v2 = v; t2 = I->x;
	  if (side == +1) v1 /= 2.;
	  side = +1;
	}
	if (v2 > v1)
	  I->x = (v1*t2 - v2*t1)/(v1 - v2);
	n++;
      }
      while (fabs (t - I->x) > 1e-5 && n < 100);
      if (fabs (t - I->x) > 1e-5)
	g_warning ("gfs_surface_segment_intersection(): convergence could not be reached\n"
		   "after %d iterations, error is %g", n, fabs (t - I->x));
    }
  }
  else
    gts_surface_foreach_face (s->s, (GtsFunc) triangle_face_intersection, I);
  return I->n;
}

static void surface_normal (GtsTriangle * t, GtsVector n)
{
  GtsVector m;
  gts_triangle_normal (t, &m[0], &m[1], &m[2]);
  n[0] -= m[0];
  n[1] -= m[1];
  n[2] -= m[2];
}

static void surface_segment_normal (GfsGenericSurface * s1,
				    FttCell * cell,
				    GfsSegment * I,
				    GtsVector n)
{
  GfsSurface * s = GFS_SURFACE (s1);
  if (s->f) {
    FttComponent c;
    GtsPoint p = segment_intersection (I);
    for (c = 0; c < FTT_DIMENSION; c++) {
      GtsPoint p1 = p;
      (&p1.x)[c] -= 1e-4;
      gdouble v1 = gfs_surface_implicit_value (s, p1);
      (&p1.x)[c] += 2e-4;
      gdouble v2 = gfs_surface_implicit_value (s, p1);
      n[c] = v2 - v1;
    }
  }
  else {
    n[0] = n[1] = n[2] = 0.;
    gts_surface_foreach_face (s->s, (GtsFunc) surface_normal, n);
  }    
}

static void add_tetrahedron_volume (GtsTriangle * t, gpointer * data)
{
  GtsPoint * p = data[0];
  GtsVertex * v1, * v2, * v3;
  gdouble * vol = data[1];
  gts_triangle_vertices (t, &v1, &v2, &v3);
  *vol += gts_point_orientation_3d (GTS_POINT (v1), GTS_POINT (v2), GTS_POINT (v3), p);
}

static gint point_is_inside_surface (GfsGenericSurface * s1, FttVector * p)
{
  GfsSurface * s = GFS_SURFACE (s1);
  GtsPoint q;
  q.x = p->x; q.y = p->y; q.z = p->z;

  if (s->f || s->bbtree) {
    gdouble v = gfs_surface_implicit_value (s, q);
    return v == 0. ? 0 : v < 0. ? -1 : 1;
  }
  else {
    gdouble vol = 0.;
    gpointer data[2];
    data[0] = &q;
    data[1] = &vol;
    gts_surface_foreach_face (s->s, (GtsFunc) add_tetrahedron_volume, data);
    fprintf (stderr, "vol: %g\n", vol);
    return vol > 0. ? -1 : 1;
  }
  return 0;
}

static void gfs_surface_class_init (GtsObjectClass * klass)
{
  klass->read = surface_read;
  klass->write = surface_write;
  klass->destroy = surface_destroy;

  GFS_GENERIC_SURFACE_CLASS (klass)->cell_is_cut = cell_is_cut;
  GFS_GENERIC_SURFACE_CLASS (klass)->segment_intersection = surface_segment_intersection;
  GFS_GENERIC_SURFACE_CLASS (klass)->segment_normal = surface_segment_normal;
  GFS_GENERIC_SURFACE_CLASS (klass)->point_is_inside = point_is_inside_surface;
}

static void gfs_surface_init (GfsSurface * s)
{
  s->scale[0] = 1.; s->scale[1] = 1.; s->scale[2] = 1.;
  s->flip = FALSE;
  s->vertex_class = gts_vertex_class ();
  s->edge_class = gts_edge_class ();
  s->face_class = gts_face_class ();
}

GfsGenericSurfaceClass * gfs_surface_class (void)
{
  static GfsGenericSurfaceClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_info = {
      "GfsSurface",
      sizeof (GfsSurface),
      sizeof (GfsGenericSurfaceClass),
      (GtsObjectClassInitFunc) gfs_surface_class_init,
      (GtsObjectInitFunc) gfs_surface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_surface_class ()), 
				  &gfs_surface_info);
  }

  return klass;
}

/**
 * gfs_surface_implicit_value:
 * @s: an (implicit) #GfsSurface.
 * @p: a #GtsPoint.
 *
 * Returns: the value of the implicit surface a location @p.
 */
gdouble gfs_surface_implicit_value (GfsSurface * s, GtsPoint p)
{
  g_return_val_if_fail (s != NULL, 0.);
  g_return_val_if_fail (s->f != NULL || s->bbtree != NULL, 0.);

  if (s->f) { /* spatial function */
    if (s->m)
      gts_point_transform (&p, s->m);
    return (s->flip ? -1. : 1.)*gfs_function_spatial_value (s->f, (FttVector *)&p.x) + 1e-20;
  }
  else { /* bounding-box tree */
    GtsBBox * bbox;
    gdouble d2 = gts_bb_tree_point_distance (s->bbtree, &p, 
					     (GtsBBoxDistFunc) gts_point_triangle_distance2, 
					     &bbox);
    return gts_point_is_inside_surface (&p, s->bbtree, TRUE) ? d2 : - d2;
  }
}

/** \endobject{GfsSurface} */
