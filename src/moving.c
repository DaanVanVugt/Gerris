/* Gerris - The GNU Flow Solver
 * Copyright (C) 2005-2009 National Institute of Water and Atmospheric Research
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
 * \brief Moving solid boundaries.
 */

#include <stdlib.h>
#include <math.h>
#include <gts.h>
#include "moving.h"
#include "simulation.h"
#include "domain.h"
#include "utils.h"
#include "ftt.h"
#include "refine.h"
#include "adaptive.h"
#include "solid.h"
#include "vof.h"
#include "surface.h"
#include "advection.h"
#include "source.h"

#define OLD_SOLID(c) (*((GfsSolidVector **) &(GFS_VALUE (c, old_solid_v))))

typedef struct {
  GfsDomain * domain;
  gdouble dt;
  FttComponent c;
  GfsVariable * div;
  GfsVariable * v;
} DivergenceData;

#include "moving2.c"

/* GfsNumberedVertex: Object */

static void numbered_vertex_read (GtsObject ** o, GtsFile * f)
{
  static glong count = 0;

  (* GTS_OBJECT_CLASS (gfs_numbered_vertex_class ())->parent_class->read) (o, f);

  GfsNumberedVertex * v = GFS_NUMBERED_VERTEX (*o);
  v->num = count++;
}

static void numbered_vertex_class_init (GtsObjectClass * klass)
{
  klass->read = numbered_vertex_read;
}

GtsVertexClass * gfs_numbered_vertex_class (void)
{
  static GtsVertexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo numbered_vertex_info = {
      "GfsNumberedVertex",
      sizeof (GfsNumberedVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) numbered_vertex_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_vertex_class ()),
				  &numbered_vertex_info);
  }

  return klass;
}

/**
 * Moving solid boundaries.
 * \beginobject{GfsSolidMoving}
 */

typedef struct {
  GfsSimulation * sim;
  GfsSolidMoving * s;
  GfsVariable * old_solid_v, ** sold2, ** v;
  GArray * stmp, * sall;
} SolidInfo;

static double surface_value (FttCell * cell, GfsVariable * v, FttVector * ca)
{
  gdouble val = 0.;
  if (!v->surface_bc)
    /* default surface BC for velocity is zero */
    return 0.;
  else if (GFS_STATE (cell)->solid) {
    FttVector oldca;
    if (ca) {
      oldca = GFS_STATE (cell)->solid->ca;
      GFS_STATE (cell)->solid->ca = *ca;
    }
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    if (ca)
      GFS_STATE (cell)->solid->ca = oldca;
    val = GFS_STATE (cell)->solid->fv;
  }
  else {
    GfsSolidVector solid;
    if (ca)
      solid.ca = *ca;
    else
      ftt_cell_pos (cell, &solid.ca);
    solid.cm = solid.ca;
    GFS_STATE (cell)->solid = &solid;
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    GFS_STATE (cell)->solid = NULL;
    val = solid.fv;
  }
  if (!(cell->flags & GFS_FLAG_DIRICHLET))
    g_assert_not_implemented ();
  return val;
}

static void init_new_cell_velocity_from_solid (FttCell * cell, SolidInfo * p)
{
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, p->v[c]) = surface_value (cell, p->v[c], NULL);
}

static void update_neighbors (FttCell * cell)
{
  gint i;
  FttCellNeighbors neighbor;
  g_assert (cell);

  ftt_cell_neighbors (cell, &neighbor);
  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i] && neighbor.c[i]->children) {
      ftt_cell_neighbors_not_cached (neighbor.c[i], &(neighbor.c[i]->children->neighbors));}
}

static gboolean refine_maxlevel (FttCell * cell, gint * maxlevel)
{
  return (ftt_cell_level (cell) < *maxlevel);
}

static void moving_cell_coarse_fine (FttCell * cell, GfsVariable * v)
{
  FttCell * parent = ftt_cell_parent (cell);

  GFS_VALUE (cell, v) = GFS_VALUE (parent, v);
  if (!GFS_CELL_IS_BOUNDARY (parent)) {
    FttVector p;
    FttComponent c;

    ftt_cell_relative_pos (cell, &p);    
    for (c = 0; c < FTT_DIMENSION; c++)
      GFS_VALUE (cell, v) += (&p.x)[c]*gfs_center_van_leer_gradient (parent, c, v->i);
  }
}

static void moving_vof_cell_coarse_fine (FttCell * cell, GfsVariable * v)
{
  FttCell * parent = ftt_cell_parent (cell);
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gdouble f = GFS_VALUE (parent, v);
  FttComponent c;
  guint i;

  if (GFS_IS_FULL (f)) {
    GFS_VALUE (cell, v) = f;
    for (c = 1; c < FTT_DIMENSION; c++)
	GFS_VALUE (cell, t->m[c]) = 0.;
      GFS_VALUE (cell, t->m[0]) = 1.;
      GFS_VALUE (cell, t->alpha) = f;
  }
  else {
    gdouble alpha = GFS_VALUE (parent, t->alpha);
    FttVector m;

    for (i = 0; i < FTT_DIMENSION; i++)
      (&m.x)[i] = GFS_VALUE (parent, t->m[i]);

    gdouble alpha1 = alpha;

    FttVector p;

    ftt_cell_relative_pos (cell, &p);
    for (c = 0; c < FTT_DIMENSION; c++) {
      alpha1 -= (&m.x)[c]*(0.25 + (&p.x)[c]);
      GFS_VALUE (cell, t->m[c]) = (&m.x)[c];
    }
    GFS_VALUE (cell, v) = gfs_plane_volume (&m, 2.*alpha1);
    GFS_VALUE (cell, t->alpha) = 2.*alpha1;
  }
}

static void moving_cell_init (FttCell * cell, SolidInfo * solid_info)
{
  GSList * i;
  gint k;
  GfsDomain * domain = GFS_DOMAIN (solid_info->sim);
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (domain)->old_solid;

  gfs_cell_init (cell, domain);

  i = domain->variables;
  while (i) {
    GfsVariable * v = i->data;
    
    if (v->coarse_fine == (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine)
      moving_cell_coarse_fine (cell, v);

    if (GFS_IS_VARIABLE_TRACER_VOF (v))
      moving_vof_cell_coarse_fine (cell, v);

    i = i->next;
  }

  g_assert (OLD_SOLID (cell) == NULL);
  OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));
  OLD_SOLID (cell)->a = 0.;
  GfsVariable ** sold2 = solid_info->sold2;
  if (sold2)
    for (k = 0; k < FTT_NEIGHBORS; k++)
      SOLD2 (cell, k) = OLD_SOLID (cell)->s[k] = 0.;

  init_new_cell_velocity_from_solid (cell, solid_info);
}

static void moving_cell_fine_init (FttCell * cell, SolidInfo * solid_info)
{
  GfsDomain * domain = GFS_DOMAIN(solid_info->sim);
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (domain)->old_solid;
  GfsVariable ** sold2 = solid_info->sold2;
  FttCellChildren child;
  guint n;
 
  gfs_cell_fine_init (cell, domain);

  /* need to update the neighbors of the "undestroyed" parent cell */
  update_neighbors (cell);
 
  ftt_cell_children (cell, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    GfsSolidVector * solid = OLD_SOLID (child.c[n]);
    gint k;
    g_assert (!solid);
    solid = OLD_SOLID (child.c[n]) = g_malloc0 (sizeof (GfsSolidVector));
    solid->a = 0.;
    if (sold2)
      for (k = 0; k < FTT_NEIGHBORS; k++)
	SOLD2 (child.c[n], k) = solid->s[k] = 0.;
  }
}

static void create_new_cells (FttCell * cell, GfsSurface * s, SolidInfo * solid_info)
{
  GfsSolidMoving * solid = solid_info->s;
  gint maxlevel = gfs_function_value (solid->level, cell);

  if (FTT_CELL_IS_DESTROYED (cell) && ftt_cell_level (cell) <= maxlevel) {
    cell->flags &= ~FTT_FLAG_DESTROYED;
    moving_cell_init (cell, solid_info);
    if (ftt_cell_level (cell) < maxlevel)
      ftt_cell_refine (cell,
		       (FttCellRefineFunc) refine_maxlevel, &maxlevel,
		       (FttCellInitFunc) moving_cell_fine_init, solid_info);
  }
  else if (ftt_cell_level (cell) < maxlevel)
    ftt_cell_refine (cell,
		     (FttCellRefineFunc) refine_maxlevel, &maxlevel,
		     (FttCellInitFunc) gfs_cell_fine_init, solid_info->sim);
}

static FttVector rpos[FTT_NEIGHBORS] = {
#if FTT_2D
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}
#else  /* FTT_3D */
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}, {0.,0.,1.}, {0.,0.,-1.}
#endif /* FTT_3D */
};

static void match (FttCell * cell, GfsBoundary * boundary)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, boundary->d);
  FttCell * parent = ftt_cell_parent (cell);
  guint level = ftt_cell_level (cell);

  cell->flags |= GFS_FLAG_BOUNDARY;
  if (parent && GFS_CELL_IS_GRADIENT_BOUNDARY (parent))
    cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
  if (neighbor == NULL || ftt_cell_level (neighbor) < level) {
    if (FTT_CELL_IS_ROOT (cell))
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "root cell is entirely outside of the fluid domain\n"
	     "the solid surface orientation may be incorrect");
    ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup, gfs_box_domain (boundary->box));
    boundary->changed = TRUE;
    return;
  }
  if (ftt_cell_level (neighbor) == level) {
    GfsSolidVector * s = GFS_STATE (neighbor)->solid;

    if (s && s->s[FTT_OPPOSITE_DIRECTION (boundary->d)] == 0.) {
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	       "root cell is entirely outside of the fluid domain\n"
	       "the solid surface orientation may be incorrect");
      ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup,
			gfs_box_domain (boundary->box));
      boundary->changed = TRUE;
      return;
    }
    if (s) {
      FttDirection d;
      FttComponent c;
      GfsSolidVector * t;

      if (GFS_STATE (cell)->solid == NULL)
	GFS_STATE (cell)->solid = g_malloc0 (sizeof (GfsSolidVector));
      t = GFS_STATE (cell)->solid;
      t->a = s->a;
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (d/2 == boundary->d/2)
	  t->s[d] = s->s[FTT_OPPOSITE_DIRECTION (d)];
	else
	  t->s[d] = s->s[d];
      for (c = 0; c < FTT_DIMENSION; c++)
	if (c == boundary->d/2) {
	  FttVector p1, p2;
	  ftt_cell_pos (cell, &p1);
	  ftt_cell_pos (neighbor, &p2);
	  (&t->cm.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->ca.x)[c];
	}
	else {
	  (&t->cm.x)[c] = (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&s->ca.x)[c];
	}
    }
    else if (GFS_STATE (cell)->solid != NULL) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }      
    if (FTT_CELL_IS_LEAF (cell) && !FTT_CELL_IS_LEAF (neighbor)) {
      GfsDomain * domain = gfs_box_domain (boundary->box);
      ftt_cell_refine_single (cell, domain->cell_init, domain->cell_init_data);
      boundary->changed = TRUE;
    }
  }
  else
    g_assert_not_reached ();
  if (!FTT_CELL_IS_LEAF (cell))
    level++;
  if (level > boundary->depth)
    boundary->depth = level;
}

static void boundary_match (GfsBoundary * boundary)
{
  guint l = ftt_cell_level (boundary->root);

  boundary->changed = FALSE;
  boundary->depth = l;
  while (l <= boundary->depth) {
    ftt_cell_traverse_boundary (boundary->root, boundary->d,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
				(FttCellTraverseFunc) match, boundary);
    l++;
  }
  if (boundary->changed)
    ftt_cell_flatten (boundary->root, boundary->d, (FttCellCleanupFunc) gfs_cell_cleanup, 
		      gfs_box_domain (boundary->box));
}

static void renew_boundary (GfsBox * box, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN(sim);
  FttDirection d;
  FttVector pos;
  gdouble size;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (box->neighbor[d]);

      /* Expensive but case independent fix */

      ftt_cell_destroy (boundary->root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
      
      domain = gfs_box_domain (box);
      boundary->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);
      FTT_ROOT_CELL (boundary->root)->parent = box;
      ftt_cell_set_level (boundary->root, ftt_cell_level (box->root));
      ftt_cell_set_neighbor_match (boundary->root, box->root, boundary->d, (FttCellInitFunc) gfs_cell_init, domain);
      ftt_cell_pos (box->root, &pos);
      size = ftt_cell_size (box->root);
      pos.x += rpos[d].x*size;
      pos.y += rpos[d].y*size;
      pos.z += rpos[d].z*size;
      ftt_cell_set_pos (boundary->root, &pos);
      boundary_match (boundary);
    }
}

static void remesh_surface_moving (GfsSimulation * sim, GfsSolidMoving * s)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  SolidInfo solid_info;

  solid_info.sim = sim;
  solid_info.s = s;
  solid_info.sold2 = GFS_SIMULATION_MOVING (sim)->sold2;
  solid_info.v = gfs_domain_velocity (domain);
  gfs_domain_traverse_cut (domain, GFS_SOLID (s)->s,
			   FTT_POST_ORDER, FTT_TRAVERSE_LEAFS | FTT_TRAVERSE_DESTROYED,
			   (FttCellTraverseCutFunc) create_new_cells, &solid_info);

  if (domain->pid >= 0)
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) renew_boundary, sim);
}

static void solid_moving_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOLID_MOVING (object)->level));
  (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->destroy) (object);
}

static void solid_moving_read (GtsObject ** o, GtsFile * fp)
{
  GfsSolidMoving * solid = GFS_SOLID_MOVING (*o);

  GFS_SURFACE (GFS_SOLID (solid)->s)->vertex_class = gfs_numbered_vertex_class ();

  if (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  solid->nvertex = gts_surface_vertex_number (GFS_SURFACE (GFS_SOLID (solid)->s)->s);

  if (!GFS_IS_SURFACE (GFS_SOLID (solid)->s) || !GFS_SURFACE (GFS_SOLID (solid)->s)->s) {
    gts_file_error (fp, "moving implicit surfaces are not implemented yet");
    return;
  }

  if (!GFS_IS_SIMULATION_MOVING (gfs_object_simulation (*o))) {
    gts_file_error (fp, "GfsSolidMoving only makes sense with GfsSimulationMoving");
    return;
  }

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }
    if (!strcmp (fp->token->str, "level")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (solid->level, gfs_object_simulation (*o), fp);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void solid_moving_write (GtsObject * object, FILE * fp)
{
  GfsSolidMoving * solid = GFS_SOLID_MOVING (object);
  
  if (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->write) 
      (object, fp);
  fputs (" { level =", fp);
  gfs_function_write (solid->level, fp);
  fputs (" }", fp);
}

static void set_old_solid (FttCell * cell, GfsVariable * old_solid_v)
{
  g_free (OLD_SOLID (cell));
  OLD_SOLID (cell) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
  cell->flags &= ~GFS_FLAG_PERMANENT;
}

static void check_face (FttCellFace * f, guint * nf)
{
  GfsSolidVector * s = GFS_STATE (f->cell)->solid;

  if (s && !f->neighbor && s->s[f->d] > 0. && s->s[f->d] < 1.)
    (*nf)++;
}

static void check_solid_fractions (GfsBox * box, guint * nf)
{
  FttDirection d;

  gfs_cell_check_solid_fractions (box->root);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    ftt_face_traverse_boundary (box->root, d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) check_face, nf);
}

static void is_diffusion (GfsSource * s, gboolean * diffusion)
{
  *diffusion = (GFS_IS_SOURCE_DIFFUSION (s) != NULL);
}

static void set_permanent (FttCell * cell)
{
  cell->flags |= GFS_FLAG_PERMANENT;
}

typedef struct {
  GfsDomain * domain;
  GfsVariable * status;
  GfsVariable ** v;
} ReInitParams;

static void redistribute_destroyed_cells_content (FttCell * cell, ReInitParams * p)
{
  if (GFS_VALUE (cell, p->status) != 1.)
    return;

  GfsDomain * domain = p->domain;
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (domain)->old_solid;
  GSList * i;
  FttCell * merged, * next;
  gdouble s1, s2;
  gint c;

  if (!OLD_SOLID (cell) || !(merged = OLD_SOLID (cell)->merged))
    return;
  while (OLD_SOLID (merged) && (next = OLD_SOLID (merged)->merged))
    merged = next;

  s1 = ftt_cell_volume (cell);
  s2 = ftt_cell_volume (merged);

  /* redistribution of the velocity */
  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble a = OLD_SOLID (merged) ? OLD_SOLID (merged)->a : 1.;
    GfsVariable * var = p->v[c];
    GFS_VALUE (merged, var) = (s1*OLD_SOLID (cell)->a*GFS_VALUE (cell, var) +
			       s2*a*GFS_VALUE (merged, var))
      /(s1*OLD_SOLID (cell)->a + s2*a);
  }
  
  /* redistribution of tracers */
  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER (i->data)) {
      gdouble a = OLD_SOLID (merged) ? OLD_SOLID (merged)->a : 1.;
      GfsVariableTracer * t = GFS_VARIABLE_TRACER(i->data);
      GfsVariable * var = t->advection.v;
      GFS_VALUE (merged, var) = (s1*OLD_SOLID (cell)->a*GFS_VALUE (cell, var) +
				 s2*a*GFS_VALUE (merged, var))
	/(s1*OLD_SOLID (cell)->a + s2*a);
    }
    i = i->next;
  }
    
  if (!OLD_SOLID (merged)) {
    OLD_SOLID (merged) = g_malloc0 (sizeof (GfsSolidVector));
    OLD_SOLID (merged)->a = 1.; 
  }
  OLD_SOLID (merged)->a += s1/s2*OLD_SOLID (cell)->a;
  if (GFS_SIMULATION (domain)->advection_params.moving_order == 2)
    redistribute_old_face (cell, merged, GFS_SIMULATION_MOVING (domain)->old_solid);
}

/**
 * gfs_domain_reinit_solid_fractions:
 * @domain: a #GfsDomain.
 * @i: a list of #GfsSolids.
 *
 * Reinitializes the solid fractions of all the cells of @domain.
 *
 * If @destroy_solid is set to %TRUE, the cells entirely contained in
 * the solid are destroyed using @cleanup as cleanup function.  
 *
 * Destroy the cells that are not fluid cells anymore when the solid 
 * have moved.
 *
 * The fluid fractions of the destroyed is redistributed.
 *
 * Returns: the number of thin cells.
 */
static guint domain_reinit_solid_fractions (GfsSimulation * sim,
					    GSList * i)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * status;

  g_return_val_if_fail (sim != NULL, 0);

  status = gfs_temporary_variable (domain);
  guint thin = gfs_init_solid_fractions_leaves (domain, i, status);

  if (sim->time.t != 0.) {
    ReInitParams rp;
    rp.domain = domain;
    rp.status = status;
    rp.v = gfs_domain_velocity (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) redistribute_destroyed_cells_content, &rp);
  }

  gfs_init_solid_fractions_from_children (domain, TRUE,
					  (FttCellCleanupFunc) gfs_cell_cleanup, domain, 
					  status);
  gts_object_destroy (GTS_OBJECT (status));
  return thin;
}

/**
 * reinit_solid_fractions:
 * @sim: a #GfsSimulation.
 *
 * Calls the domain_reinit_solid_fractions(). Matches the
 * boundaries by calling gfs_domain_match().
 */
static void reinit_solid_fractions (GfsSimulation * sim)
{
  guint nf = 0;
  GfsDomain * domain = GFS_DOMAIN (sim);;
  GSList * solids = gfs_simulation_get_solids (sim);
  if (solids) {
    sim->thin = domain_reinit_solid_fractions (sim, solids);
    g_slist_free (solids);
    gfs_domain_match (domain);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) set_permanent, NULL);
  }
  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) check_solid_fractions, &nf);
  if (nf > 0) {
    GSList * i = domain->variables;
    gboolean diffusion = FALSE;
    
    while (i && !diffusion) {
      GfsVariable * v = i->data;

      if (v->sources)
	gts_container_foreach (v->sources, (GtsFunc) is_diffusion, &diffusion);
      i = i->next;
    }
    if (diffusion)
      g_warning ("the solid surface cuts %d boundary cells,\n"
		 "this may cause errors for diffusion terms\n", nf);
  }
}

/* see gfs_advection_update() for a description of what this function does */
static void moving_advection_update (GSList * merged, const GfsAdvectionParams * par)
{
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (par->v->domain)->old_solid;

  if (merged->next == NULL) { /* cell is not merged */
    FttCell * cell = merged->data;
    gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
    gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;

    if (GFS_IS_MIXED (cell))
      g_assert (!gfs_cell_is_small (cell));

    GFS_VALUE (cell, par->v) = (olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv))/a;
  }
  else if (1 /* par->average */) {
    /* average value */
    GSList * i = merged;
    gdouble w = 0., total_vol = 0.;

    while (i) {
      FttCell * cell = i->data;
      gdouble vol = ftt_cell_volume (cell);
      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
      gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
      
      total_vol += vol*a;
      w += vol*(olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv));
      i = i->next;
    }
    w /= total_vol;

    i = merged;
    while (i) {
      FttCell * cell = i->data;
      GFS_VALUE (cell, par->v) = w;
      i = i->next;
    }
  }
  else {
    GSList * i = merged;
    gdouble w = 0., total_vol = 0.;

    while (i) {
      FttCell * cell = i->data;
      gdouble vol = ftt_cell_volume (cell);
      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
      gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;

      total_vol += vol*a;
      if (a < GFS_SMALL) {
	GFS_VALUE (cell, par->v) = olda*GFS_VALUE (cell, par->v)/a + 
	  GFS_VALUE (cell, par->fv)/GFS_SMALL;
	w += vol*GFS_VALUE (cell, par->fv)*(1. - a/GFS_SMALL);   
      }
      else
	GFS_VALUE (cell, par->v) = (olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv))/a;
      i = i->next;
    }
    w /= total_vol;
    
    i = merged;
    while (i) {
      FttCell * cell = i->data;
      /* fixme: small cells should be excluded here?? 
	 (with corresponding modification in total_vol) */
      GFS_VALUE (cell, par->v) += w;
      i = i->next;
    }
  }
}

static void moving_init (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN(sim);
  GSList * i = domain->variables;

  if (sim->advection_params.moving_order == 2)
    sim->advection_params.flux = moving_face_velocity_advection_flux;
  else
    sim->advection_params.flux = gfs_face_velocity_advection_flux;
  sim->advection_params.update = (GfsMergedTraverseFunc) moving_advection_update;

  while (i) {
    if (GFS_IS_VARIABLE_TRACER (i->data)) {
      GfsAdvectionParams * par = &GFS_VARIABLE_TRACER (i->data)->advection;
      if (sim->advection_params.moving_order == 2)
	par->flux = moving_face_advection_flux;
      else
	par->flux = gfs_face_advection_flux;
      par->update = sim->advection_params.update;
      par->moving_order = sim->advection_params.moving_order;
    }
    i = i->next;
  }
}

static gboolean solid_moving_event (GfsEvent * event, GfsSimulation * sim)
{
  return (GFS_SOLID_MOVING (event)->active = 
	  (* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class)->event) 
	  (event, sim));
}

static void solid_moving_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = solid_moving_destroy;
  GTS_OBJECT_CLASS (klass)->read = solid_moving_read;
  GTS_OBJECT_CLASS (klass)->write = solid_moving_write;
  klass->event = solid_moving_event;
}

static void solid_moving_init (GfsSolidMoving * solid)
{
  gfs_event_set (GFS_EVENT (solid),
		 0., G_MAXDOUBLE/2., -1.,
		 0, G_MAXINT/2, 1);
  solid->level = gfs_function_new (gfs_function_class (), 0.);
}

GfsEventClass * gfs_solid_moving_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo solid_moving_info = {
      "GfsSolidMoving",
      sizeof (GfsSolidMoving),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) solid_moving_class_init,
      (GtsObjectInitFunc) solid_moving_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_class ()), &solid_moving_info);
  }

  return klass;
}

/** \endobject{GfsSolidMoving} */

#define MOVING_CFL 0.45

/**
 * Euler solver with moving solid boundaries.
 * \beginobject{GfsSimulationMoving}
 */

static void set_dtmax (FttCell * cell, SolidInfo * p)
{
  gdouble size = ftt_cell_size (cell);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble v = fabs (surface_value (cell, p->v[c], NULL));
    if (v != 0.) {
      gdouble dt = size*MOVING_CFL/v;
      if (dt < p->sim->time.dtmax)
	p->sim->time.dtmax = dt;
    }
  }
}

static void simulation_moving_set_timestep (GfsSimulation * sim)
{
  gdouble dtmax = sim->time.dtmax;
  SolidInfo p;
  p.sim = sim;
  p.v = gfs_domain_velocity (GFS_DOMAIN (sim));
  gfs_domain_traverse_mixed (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) set_dtmax, &p);
  gfs_simulation_set_timestep (sim);
  sim->time.dtmax = dtmax;
}

static void move_vertex (GfsNumberedVertex * v, SolidInfo * par)
{ 
  GtsPoint * p = &GTS_VERTEX(v)->p;
  FttVector pos = *((FttVector *) &p->x);
  FttCell * cell = gfs_domain_locate (GFS_DOMAIN (par->sim), pos, -2, NULL);
  FttComponent c;

  if (cell) {
    gdouble dt = par->sim->advection_params.dt;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&p->x)[c] += surface_value (cell, par->v[c], &pos)*dt;
  }
}

#ifdef HAVE_MPI

static void move_vertex_mpi (GfsNumberedVertex * v, SolidInfo * par)
{ 
  GtsPoint * p = &GTS_VERTEX(v)->p;
  FttVector pos = *((FttVector *) &p->x);
  FttCell * cell = gfs_domain_locate (GFS_DOMAIN (par->sim), pos, -2, NULL);
  FttComponent c;

  if (cell) {
    gdouble dt = par->sim->advection_params.dt;
    for (c = 0; c < FTT_DIMENSION; c++) { 
      (&p->x)[c] += surface_value (cell, par->v[c], &pos)*dt;
      g_array_index(par->stmp, double, FTT_DIMENSION*v->num+c) = (&p->x)[c];
    }  
  }
  else {
    for (c = 0; c < FTT_DIMENSION; c++)
      g_array_index(par->stmp,double,FTT_DIMENSION*v->num+c) = -G_MAXDOUBLE;
  }
}

static void synchronize_vertex (GfsNumberedVertex * v, SolidInfo * par)
{
  GtsPoint * p = &GTS_VERTEX(v)->p;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    (&p->x)[c] = g_array_index(par->sall,double,FTT_DIMENSION*v->num+c);
}

#endif /* HAVE_MPI */

static void solid_move_remesh (GfsSolidMoving * solid, GfsSimulation * sim)
{
  GfsSurface * surface = GFS_SURFACE (GFS_SOLID (solid)->s);
  if (surface->s) {
    SolidInfo p;
    p.sim = sim;
    p.s = solid;
    p.v = gfs_domain_velocity (GFS_DOMAIN (sim));

#ifdef HAVE_MPI
    if (GFS_DOMAIN (sim)->pid >= 0) { /* Parallel simulation */
      p.stmp = g_array_set_size ( g_array_new (FALSE, FALSE, sizeof (double)) , FTT_DIMENSION*solid->nvertex);
      p.sall = g_array_set_size ( g_array_new (FALSE, FALSE, sizeof (double)) , FTT_DIMENSION*solid->nvertex);

      gts_surface_foreach_vertex (surface->s, (GtsFunc) move_vertex_mpi, &p);

      MPI_Allreduce (p.stmp->data, p.sall->data, p.stmp->len, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      gts_surface_foreach_vertex (surface->s, (GtsFunc) synchronize_vertex, &p);
      
      g_array_free (p.stmp, FALSE);
      g_array_free (p.sall, FALSE);
    }
    else
#endif /* HAVE_MPI */
      gts_surface_foreach_vertex (surface->s, (GtsFunc) move_vertex, &p);
  }
  else
    /* implicit surface */
    g_assert_not_implemented ();
  remesh_surface_moving (sim, solid);
}

static void move_solids (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * old_solid = GFS_SIMULATION_MOVING (sim)->old_solid;
  GfsVariable * sold2[FTT_NEIGHBORS];

  gfs_domain_timer_start (domain, "move_solids");

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) set_old_solid, old_solid);

  if (sim->advection_params.moving_order == 2) {
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      sold2[d] = gfs_domain_add_variable (domain, NULL, NULL);
      sold2[d]->coarse_fine = sold2_fine_init;
    }
    GFS_SIMULATION_MOVING (sim)->sold2 = sold2;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) set_sold2, sim);
  }

  GSList * solids = gfs_simulation_get_solids (sim), * s = solids;
  while (s) {
    if (GFS_IS_SOLID_MOVING (s->data) && GFS_SOLID_MOVING (s->data)->active)
      solid_move_remesh (s->data, sim);
    s = s->next;
  }
  g_slist_free (solids);
  reinit_solid_fractions (sim);
  gfs_domain_reshape (domain, gfs_domain_depth (domain));

  if (sim->advection_params.moving_order == 2) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) second_order_face_fractions, sim);
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      gts_object_destroy (GTS_OBJECT (sold2[d]));    
    GFS_SIMULATION_MOVING (sim)->sold2 = NULL;
  }

  gfs_domain_timer_stop (domain, "move_solids");
}

static void moving_divergence_approx (FttCell * cell, DivergenceData * p)
{
  GFS_VALUE (cell, p->div) += 
    GFS_STATE (cell)->solid->fv*(GFS_STATE (cell)->solid->s[2*p->c + 1] -
				 GFS_STATE (cell)->solid->s[2*p->c])*ftt_cell_size (cell);
}

static void moving_divergence_distribution (GSList * merged, DivergenceData * p)
{
  if (merged->next != NULL && merged->next->data != merged->data) {
    gdouble total_volume = 0., total_div = 0.;
    GSList * i = merged;

    while (i) {
      FttCell * cell = i->data;
      g_assert (FTT_CELL_IS_LEAF (cell));
      gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
      total_volume += a*ftt_cell_volume (cell);
      total_div += GFS_VALUE (cell, p->div);
      i = i->next;
    }
    
    total_div /= total_volume;
    
    i = merged;
    while (i) {
      FttCell * cell = i->data;
      gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
      GFS_VALUE (cell, p->div) = total_div*a*ftt_cell_volume (cell);
      i = i->next;
    }
  }
}

static void divergence_approx_hook (GfsDomain * domain, 
				    gdouble dt,
				    GfsVariable * div)
{
  DivergenceData q;
  GfsVariable ** v = gfs_domain_velocity (domain);
  
  q.div = div;
  for (q.c = 0; q.c < FTT_DIMENSION; q.c++) {
    gfs_domain_surface_bc (domain, v[q.c]);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) moving_divergence_approx, &q);
  }
  gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) moving_divergence_distribution, &q);
}

static void moving_divergence_mac (FttCell * cell, DivergenceData * p)
{
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (p->domain)->old_solid;
  gdouble size = ftt_cell_size (cell);
  gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
  gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
  
  GFS_VALUE (cell, p->div) += (olda - a)*size*size/p->dt;
}

static void divergence_mac_hook_order_1 (GfsDomain * domain,
					 gdouble dt,
					 GfsVariable * div)
{
  DivergenceData q;

  q.dt = - 2.*dt;
  q.div = div;
  q.domain = domain;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) moving_divergence_mac, &q);
  gfs_domain_traverse_merged (domain,
			      (GfsMergedTraverseFunc) 
			      moving_divergence_distribution,
			      &q);
}

static void divergence_mac_hook_order_2 (GfsDomain * domain, 
					 gdouble dt,
					 GfsVariable * div)
{
  DivergenceData q;

  q.dt = 2.*dt;
  q.div = div;
  q.domain = domain;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) moving_divergence_mac, &q);
  gfs_domain_traverse_merged (domain,
			      (GfsMergedTraverseFunc) 
			       moving_divergence_distribution_second_order,
			      &q);
}

static void moving_mac_projection (GfsSimulation * sim,
				   GfsMultilevelParams * par,
				   GfsAdvectionParams * apar,
				   GfsVariable * p,
				   GfsFunction * alpha,
				   GfsVariable ** g)
{
  if (apar->moving_order == 2)
    swap_face_fractions (sim);
  gfs_mac_projection (GFS_DOMAIN (sim), par, apar->dt/2., p, alpha, g, 
		      (apar->moving_order == 2 ? 
		       divergence_mac_hook_order_2 : divergence_mac_hook_order_1));
  if (apar->moving_order == 2)
    swap_face_fractions_back (sim);
}

static void simulation_moving_run (GfsSimulation * sim)
{
  GfsVariable * p, * pmac, * res = NULL, * g[FTT_DIMENSION], * gmac[FTT_DIMENSION];
  GfsVariable ** gc = sim->advection_params.gc ? g : NULL;
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gmac[c] = gfs_temporary_variable (domain);
    if (sim->advection_params.gc)
      g[c] = gfs_temporary_variable (domain);
    else
      g[c] = gmac[c];
  }
  gfs_variable_set_vector (gmac, FTT_DIMENSION);
  gfs_variable_set_vector (g, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  moving_init (sim);

  simulation_moving_set_timestep (sim);
  if (sim->time.i == 0)
    gfs_approximate_projection (domain,
				&sim->approx_projection_params,
				sim->advection_params.dt,
				p, sim->physical_params.alpha, res, g,
				divergence_approx_hook);
  else if (sim->advection_params.gc)
    gfs_update_gradients (domain, p, sim->physical_params.alpha, g);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    move_solids (sim);
   
    gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);
    
    gfs_variables_swap (p, pmac);

    moving_mac_projection (sim,
			   &sim->projection_params, 
			   &sim->advection_params,
			   p, sim->physical_params.alpha, gmac);

    gfs_variables_swap (p, pmac);
    
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);
   
    gfs_centered_velocity_advection_diffusion (domain,
					       FTT_DIMENSION,
					       &sim->advection_params,
					       gmac,
					       sim->time.i > 0 || !gc ? gc : gmac,
					       sim->physical_params.alpha);
        
    gfs_advance_tracers (sim, sim->advection_params.dt);

    if (gc) {
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, sim->time.i > 0 ? gc : gmac, 
				       -sim->advection_params.dt);
    }
    else if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, sim->advection_params.dt);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, -sim->advection_params.dt);
    }

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    gfs_approximate_projection (domain,
				&sim->approx_projection_params, 
				sim->advection_params.dt, p, sim->physical_params.alpha, res, g,
				divergence_approx_hook);

    sim->time.t = sim->tnext;
    sim->time.i++;

    simulation_moving_set_timestep (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) {
    gts_object_destroy (GTS_OBJECT (gmac[c]));
    if (sim->advection_params.gc)
      gts_object_destroy (GTS_OBJECT (g[c]));
  }
}

static void simulation_moving_class_init (GfsSimulationClass * klass)
{
  klass->run = simulation_moving_run;
}

static void old_solid_cleanup (FttCell * cell, GfsVariable * old_solid_v)
{
  g_free (OLD_SOLID (cell));
  OLD_SOLID (cell) = NULL;
}

static void none (void) {}

static void simulation_moving_init (GfsDomain * domain)
{
  gfs_domain_add_variable (domain, "div", "Divergence")->centered = TRUE;

  /* old_solid will hold a pointer to a GfsSolidVector */
  GfsVariable * old_solid = gfs_domain_add_variable (domain, NULL, NULL); 
  GFS_SIMULATION_MOVING (domain)->old_solid = old_solid;
  /* pointers need to be "interpolated" correctly (i.e. not at all) */
  old_solid->coarse_fine = (GfsVariableFineCoarseFunc) none;
  old_solid->fine_coarse = (GfsVariableFineCoarseFunc) none;
  /* the memory needs to be freed when the cell is cleaned up */
  old_solid->cleanup = (FttCellCleanupFunc) old_solid_cleanup;
  /* switch off boundary conditions */
  GfsBc * bc = gfs_bc_new (gfs_bc_class (), old_solid, FALSE);
  bc->bc = bc->homogeneous_bc = bc->face_bc = (FttFaceTraverseFunc) none;
  gfs_variable_set_default_bc (old_solid, bc);
}

GfsSimulationClass * gfs_simulation_moving_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_simulation_moving_info = {
      "GfsSimulationMoving",
      sizeof (GfsSimulationMoving),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) simulation_moving_class_init,
      (GtsObjectInitFunc) simulation_moving_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), 
				  &gfs_simulation_moving_info);
  }

  return klass;
}

/** \endobject{GfsSimulationMoving} */
