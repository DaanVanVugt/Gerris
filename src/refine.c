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
 * \brief Static mesh refinement.
 */

#include <stdlib.h>
#include "refine.h"
#include "solid.h"
#include "adaptive.h"
#include "init.h"

/**
 * Simple definition of the refinement levels.
 * \beginobject{GfsRefine}
 */

static gboolean refine_maxlevel (FttCell * cell, GfsFunction * maxlevel)
{
  return (ftt_cell_level (cell) < gfs_function_value (maxlevel, cell));
}

static void refine_box (GfsBox * box, GfsFunction * maxlevel)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_maxlevel, maxlevel,
		   (FttCellInitFunc) gfs_cell_fine_init, gfs_box_domain (box));
}

static void gfs_refine_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gfs_catch_floating_point_exceptions ();
  gts_container_foreach (GTS_CONTAINER (sim),
			 (GtsFunc) refine_box, refine->maxlevel);
  gfs_restore_fpe_for_function (refine->maxlevel);
}

static void gfs_refine_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_REFINE (o)->maxlevel));
  (* GTS_OBJECT_CLASS (gfs_refine_class ())->parent_class->destroy) (o);
}

static void gfs_refine_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, "%s", object->klass->info.name);
  gfs_function_write (GFS_REFINE (object)->maxlevel, fp);
}

static void gfs_refine_read (GtsObject ** o, GtsFile * fp)
{
  GfsRefine * refine = GFS_REFINE (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsRefineClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_refine_class ())) {
    gts_file_error (fp, "`%s' is not a GfsRefine", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (refine));
    refine = GFS_REFINE (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  gfs_function_read (refine->maxlevel, gfs_object_simulation (refine), fp);
  if (fp->type == GTS_ERROR)
    return;

  if (class_changed && fp->type != '\n' && klass->read)
    (* klass->read) (o, fp);
}

static void gfs_refine_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_refine;

  GTS_OBJECT_CLASS (klass)->destroy = gfs_refine_destroy;
  GTS_OBJECT_CLASS (klass)->write = gfs_refine_write;
  GTS_OBJECT_CLASS (klass)->read =  gfs_refine_read;
}

static void gfs_refine_init (GfsRefine * object)
{
  object->maxlevel = gfs_function_new (gfs_function_class (), 1.);
}

GfsRefineClass * gfs_refine_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_info = {
      "GfsRefine",
      sizeof (GfsRefine),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_class_init,
      (GtsObjectInitFunc) gfs_refine_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
			    &gfs_refine_info);
  }

  return klass;
}

GfsRefine * gfs_refine_new (GfsRefineClass * klass)
{
  GfsRefine * object;

  object = GFS_REFINE (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/** \endobject{GfsRefine} */

/**
 * Refine embedded solid surfaces.
 * \beginobject{GfsRefineSolid}
 */

typedef struct _GfsRefineSolid           GfsRefineSolid;

struct _GfsRefineSolid {
  GfsRefine parent;

  GfsDerivedVariable * v;
};
  
#define GFS_REFINE_SOLID(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineSolid,\
					           gfs_refine_solid_class ())

static void refine_solid_destroy (GtsObject * object)
{
  gfs_domain_remove_derived_variable (GFS_DOMAIN (gfs_object_simulation (object)), 
				      "SolidCurvature");

  (* GTS_OBJECT_CLASS (gfs_refine_solid_class ())->parent_class->destroy) (object);
}

typedef struct {
  GtsSurface * s;
  gdouble kappa;
} KappaData;

static void max_kappa (GtsVertex * v, KappaData * d)
{
  GtsVector Kh;

  if (gts_vertex_mean_curvature_normal (v, d->s, Kh)) {
    gdouble kappa = gts_vector_norm (Kh)/(FTT_DIMENSION - 1);
    if (kappa > d->kappa)
      d->kappa = kappa;
  }
}

static gdouble solid_curvature (FttCell * cell, FttCellFace * face, 
				GfsDomain * domain, GfsGenericSurface * s)
{
  KappaData d;
  d.kappa = gfs_solid_is_thin (cell, s) ? 1./ftt_cell_size (cell) : 0.;
  d.s = GFS_SURFACE (s)->s;
  gts_surface_foreach_vertex (d.s, (GtsFunc) max_kappa, &d);
  return d.kappa;
}

static void refine_solid_read (GtsObject ** o, GtsFile * fp)
{
  GfsRefineSolid * refine = GFS_REFINE_SOLID (*o);
  GfsDerivedVariableInfo v = { "SolidCurvature", "curvature of the solid boundary",
			       solid_curvature };
  refine->v = gfs_domain_add_derived_variable (GFS_DOMAIN (gfs_object_simulation (*o)), v);
  if (!refine->v) {
    gts_file_error (fp, "derived variable `SolidCurvature' already defined");
    return;
  }

  (* GTS_OBJECT_CLASS (gfs_refine_solid_class ())->parent_class->read) (o, fp);
}

typedef struct {
  GfsRefine * refine;
  GfsDomain * domain;
  GfsGenericSurface * surface;
} RefineCut;

static void refine_cut_cell (FttCell * cell, GfsGenericSurface * s, RefineCut * p)
{
  GTS_OBJECT (s)->reserved = p->surface;
  GFS_REFINE_SOLID (p->refine)->v->data = s;
  if (ftt_cell_level (cell) < gfs_function_value (p->refine->maxlevel, cell))
    ftt_cell_refine_single (cell, p->domain->cell_init, p->domain->cell_init_data);
  GFS_REFINE_SOLID (p->refine)->v->data = NULL;
}

static void refine_implicit_cell (FttCell * cell, RefineCut * p)
{
  guint maxlevel = gfs_function_value (p->refine->maxlevel, cell);
  if (ftt_cell_level (cell) < maxlevel && gfs_cell_is_cut (cell, p->surface, FALSE, maxlevel))
    ftt_cell_refine_single (cell, p->domain->cell_init, p->domain->cell_init_data);
}

static void gfs_refine_solid_refine (GfsRefine * refine, GfsSimulation * sim)
{
  if (sim->solids) {
    RefineCut p;
    p.refine = refine;
    p.domain = GFS_DOMAIN (sim);
    GSList * i = sim->solids->items;
    while (i) {
      p.surface = GFS_SOLID (i->data)->s;
      gfs_catch_floating_point_exceptions ();
      if (GFS_SURFACE (p.surface)->s)
	gfs_domain_traverse_cut (GFS_DOMAIN (sim), p.surface,
				 FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseCutFunc) refine_cut_cell, &p);
      else
	gfs_domain_cell_traverse (GFS_DOMAIN (sim),
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) refine_implicit_cell, &p);
      gfs_restore_fpe_for_function (refine->maxlevel);
      i = i->next;
    }
  }
}

static void gfs_refine_solid_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_solid_refine;
  GTS_OBJECT_CLASS (klass)->destroy = refine_solid_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_solid_read;
}

GfsRefineClass * gfs_refine_solid_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_solid_info = {
      "GfsRefineSolid",
      sizeof (GfsRefineSolid),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_solid_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_solid_info);
  }

  return klass;
}

/** \endobject{GfsRefineSolid} */

/**
 * Refine cells cut by a surface.
 * \beginobject{GfsRefineSurface}
 */

static void refine_surface_destroy (GtsObject * object)
{
  GfsRefineSurface * d = GFS_REFINE_SURFACE (object);

  gts_object_destroy (GTS_OBJECT (d->surface));

  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->destroy) (object);
}

static void refine_surface_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->write) (o, fp);
  gfs_generic_surface_write (GFS_REFINE_SURFACE (o)->surface, gfs_object_simulation (o), fp);
}

static void refine_surface_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_generic_surface_read (GFS_REFINE_SURFACE (*o)->surface, gfs_object_simulation (*o), fp);
}

static void gfs_refine_surface_refine (GfsRefine * refine, GfsSimulation * sim)
{
  RefineCut p;

  p.refine = refine;
  p.domain = GFS_DOMAIN (sim);
  p.surface = GFS_REFINE_SURFACE (refine)->surface;
  if (GFS_SURFACE (p.surface)->twod) {
    if (GFS_SURFACE (p.surface)->s)
      gfs_domain_traverse_cut_2D (GFS_DOMAIN (sim), GFS_REFINE_SURFACE (refine)->surface,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				  (FttCellTraverseCutFunc) refine_cut_cell, &p);
    else
      g_assert_not_implemented ();
  }
  else {
    if (GFS_SURFACE (p.surface)->s)
      gfs_domain_traverse_cut (GFS_DOMAIN (sim), GFS_REFINE_SURFACE (refine)->surface,
			       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseCutFunc) refine_cut_cell, &p);
    else
      gfs_domain_cell_traverse (GFS_DOMAIN (sim),
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) refine_implicit_cell, &p);
  }
}

static void gfs_refine_surface_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_surface_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_surface_destroy;
  GTS_OBJECT_CLASS (klass)->write = refine_surface_write;
  GTS_OBJECT_CLASS (klass)->read = refine_surface_read;
}

static void refine_surface_init (GfsRefineSurface * r)
{
  r->surface = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
}

GfsRefineClass * gfs_refine_surface_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_surface_info = {
      "GfsRefineSurface",
      sizeof (GfsRefineSurface),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_surface_class_init,
      (GtsObjectInitFunc) refine_surface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_surface_info);
  }

  return klass;
}

/** \endobject{GfsRefineSurface} */

/**
 * Refine cells as a function of the distance from a surface.
 * \beginobject{GfsRefineDistance}
 */

static void refine_distance_destroy (GtsObject * object)
{
  GfsRefineDistance * d = GFS_REFINE_DISTANCE (object);

  if (d->stree)
    gts_bb_tree_destroy (d->stree, TRUE);
  gfs_domain_remove_derived_variable (GFS_DOMAIN (gfs_object_simulation (object)), "Distance");

  (* GTS_OBJECT_CLASS (gfs_refine_distance_class ())->parent_class->destroy) (object);
}

static gdouble cell_distance (FttCell * cell, 
			      FttCellFace * face, 
			      GfsSimulation * sim,
			      GfsRefineDistance * refine)
{
  FttVector pos;
  gdouble h = GFS_DIAGONAL*ftt_cell_size (cell), d;
  GtsPoint p;

  ftt_cell_pos (cell, &pos);
  p.x = pos.x; p.y = pos.y; p.z = pos.z;
  d = gts_bb_tree_point_distance (refine->stree, &p,
				  (GtsBBoxDistFunc) gts_point_triangle_distance, NULL);
  return d > h ? d - h : 0.;
}

static void refine_distance_read (GtsObject ** o, GtsFile * fp)
{
  GfsDerivedVariableInfo v = { "Distance", "distance to the surface", 
			       cell_distance };

  v.data = *o;
  if (!gfs_domain_add_derived_variable (GFS_DOMAIN (gfs_object_simulation (*o)), v)) {
    gts_file_error (fp, "derived variable `Distance' already defined");
    return;
  }

  (* GTS_OBJECT_CLASS (gfs_refine_distance_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GtsSurface * s = GFS_SURFACE (GFS_REFINE_SURFACE (*o)->surface)->s;
  if (!s) {
    gts_file_error (fp, "RefineDistance only works with GTS surfaces");
    return;
  }

  GFS_REFINE_DISTANCE (*o)->stree = gts_bb_tree_surface (s);
}

static void gfs_refine_distance_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_distance_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_distance_read;
}

GfsRefineClass * gfs_refine_distance_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_distance_info = {
      "GfsRefineDistance",
      sizeof (GfsRefineDistance),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_distance_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_surface_class ()),
				  &gfs_refine_distance_info);
  }

  return klass;
}

/** \endobject{GfsRefineDistance} */

/**
 *
 * \beginobject{GfsRefineHeight}
 */

static void refine_height_destroy (GtsObject * object)
{
  gfs_domain_remove_derived_variable (GFS_DOMAIN (gfs_object_simulation (object)), "Height");

  (* GTS_OBJECT_CLASS (gfs_refine_height_class ())->parent_class->destroy) (object);
}

static gdouble interpolated_value (GtsSurface * s, FttVector * p)
{
  GtsPoint q;
  GtsFace * t;

  q.x = p->x; q.y = p->y;
  t = gts_point_locate (&q, s, NULL);
  if (t == NULL) {
    g_warning ("cannot locate point (%g,%g)", p->x, p->y);
    return 0.;
  }
  gts_triangle_interpolate_height (GTS_TRIANGLE (t), &q);
  return q.z;
}

static gdouble cell_height (FttCell * cell, 
			    FttCellFace * face, 
			    GfsSimulation * sim,
			    GfsRefineSurface * refine)
{
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  return interpolated_value (GFS_SURFACE (refine->surface)->s, &pos);
}

static void refine_height_read (GtsObject ** o, GtsFile * fp)
{
  GfsDerivedVariableInfo v = { "Height", "vertical distance to the surface", 
			       cell_height };

  v.data = *o;
  if (!gfs_domain_add_derived_variable (GFS_DOMAIN (gfs_object_simulation (*o)), v)) {
    gts_file_error (fp, "derived variable `Height' already defined");
    return;
  }

  (* GTS_OBJECT_CLASS (gfs_refine_height_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!GFS_SURFACE (GFS_REFINE_SURFACE (*o)->surface)->s) {
    gts_file_error (fp, "RefineHeight only works with GTS surfaces");
    return;
  }
}

static void gfs_refine_height_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_height_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_height_read;
}

GfsRefineClass * gfs_refine_height_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_height_info = {
      "GfsRefineHeight",
      sizeof (GfsRefineSurface),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_height_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_surface_class ()),
				  &gfs_refine_height_info);
  }

  return klass;
}

/** \endobject{GfsRefineHeight} */
