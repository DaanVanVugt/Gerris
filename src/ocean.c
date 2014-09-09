/* Gerris - The GNU Flow Solver
 * Copyright (C) 2004-2012 Stéphane Popinet
 * National Institute of Water and Atmospheric Research
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
 * \brief The linearised shallow-water solver.
 */

#include <stdlib.h>

#include "ocean.h"
#include "timestep.h"
#include "adaptive.h"
#include "source.h"
#include "vof.h"
#include "graphic.h"

#include "solid.h"

/**
 * The linearised shallow-water solver.
 * \beginobject{GfsOcean}
 */

static void correct_normal_velocity (FttCellFace * face,
				     gpointer * data)
{
  GfsGradient g;
  gdouble dp;
  FttFaceType type;
  GfsVariable * p = data[0];
  GfsVariable ** gv = data[1];
  gdouble * dt = data[2];
  FttComponent c;

  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  type = ftt_face_type (face);
  c = face->d/2;

  gfs_face_gradient (face, &g, p->i, -1);
  dp = (g.b - g.a*GFS_VALUE (face->cell, p))/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    dp = - dp;

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= dp*(*dt);
  GFS_VALUE (face->cell, gv[c]) += dp;

  if (ftt_face_type (face) == FTT_FINE_COARSE)
    dp *= GFS_FACE_FRACTION_LEFT (face)/(GFS_FACE_FRACTION_RIGHT (face)*FTT_CELLS/2);
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp*(*dt);
  GFS_VALUE (face->neighbor, gv[c]) += dp;
}

static void scale_gradients (FttCell * cell, gpointer * data)
{
  GfsVariable ** g = data[0];
  guint * dimension = data[1];
  FttCellNeighbors n;
  FttComponent c;

  ftt_cell_neighbors (cell, &n);
  for (c = 0; c < *dimension; c++) {
    FttCell * c1 = n.c[2*c], * c2 = n.c[2*c + 1];
    
    if (c1 && c2 && !GFS_CELL_IS_GRADIENT_BOUNDARY (c1) && !GFS_CELL_IS_GRADIENT_BOUNDARY (c2))
      GFS_VALUE (cell, g[c]) /= 2.;
  }
}

/**
 * gfs_correct_normal_velocities_weighted:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @p: the pressure field.
 * @g: where to store the pressure gradient.
 * @dt: the timestep.
 * @weighted: whether to use fraction-weighting or not.
 *
 * Corrects the normal velocity field of @domain using @p and and @dt.
 *
 * Also allocates the @g variables and fills them with the centered gradient of @p.
 */
static void gfs_correct_normal_velocities_weighted (GfsDomain * domain,
						    guint dimension,
						    GfsVariable * p,
						    GfsVariable ** g,
						    gdouble dt,
						    gboolean weighted)
{
  FttComponent c;
    
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);
    
  for (c = 0; c < dimension; c++)
    g[c] = gfs_temporary_variable (domain);
  gfs_variable_set_vector (g, dimension);
  gfs_reset_gradients (domain, dimension, g);
  if (weighted) {
    gfs_correct_normal_velocities (domain, dimension, p, g, dt);
    gfs_scale_gradients (domain, dimension, g);
  }
  else {
    gpointer data[3];
    data[0] = p;
    data[1] = g;
    data[2] = &dt;
    gfs_domain_face_traverse (domain, dimension == 2 ? FTT_XY : FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) correct_normal_velocity, data);
    data[0] = g;
    data[1] = &dimension;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) scale_gradients, data);
    for (c = 0; c < dimension; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, g[c]);
  }
}

#define THETA 0.5

typedef struct {
  GfsVariable * pn, * div, * divn, * dia;
  gdouble dt, G;
} FreeSurfaceParams;

static void normal_divergence (FttCell * cell, FreeSurfaceParams * p)
{
  GFS_VALUE (cell, p->div) += (1. - THETA)*GFS_VALUE (cell, p->divn)/THETA;
}

static void scale_divergence_helmoltz (FttCell * cell, FreeSurfaceParams * p)
{
  gdouble h = ftt_cell_size (cell);
  gdouble c = 2.*h*h/(THETA*p->G*p->dt*p->dt);

#if FTT_2D
  c *= gfs_domain_cell_fraction (p->dia->domain, cell);
#else /* 3D */
  if (GFS_IS_MIXED (cell)) /* fixme: no metric yet */
    c *= GFS_STATE (cell)->solid->s[FTT_FRONT];
#endif /* 3D */

  GFS_VALUE (cell, p->dia) = c;
  GFS_VALUE (cell, p->div) = 2.*GFS_VALUE (cell, p->div)/p->dt -
    c*GFS_VALUE (cell, p->pn);
}

#if !FTT_2D
static void merge_pressures (GSList * merged, GfsVariable * v)
{
  if (merged->next != NULL) {
    /* average value */
    GSList * i = merged;
    gdouble w = 0., total_area = 0.;

    while (i) {
      FttCell * cell = i->data;
      GfsSolidVector * solid = GFS_STATE (cell)->solid;
      gdouble h = ftt_cell_size (cell);
      gdouble area = h*h*solid->s[FTT_FRONT];
      total_area += area;
      w += area*GFS_VALUE (cell, v);
      i = i->next;
    }
    w /= total_area;

    i = merged;
    while (i) {
      FttCell * cell = i->data;
      GFS_VALUE (cell, v) = w;
      i = i->next;
    }
  }
}
#endif /* 3D */

/**
 * gfs_free_surface_pressure:
 * @toplayer: a #GfsDomain.
 * @par: the multigrid paramaters.
 * @apar: the advection parameters.
 *
 */
static void gfs_free_surface_pressure (GfsDomain * toplayer,
				       GfsMultilevelParams * par,
				       GfsAdvectionParams * apar,
				       GfsVariable * p,
				       GfsVariable * div,
				       GfsVariable * divn,
				       GfsVariable * res,
				       gdouble G)
{
  FreeSurfaceParams fp;
  GfsVariable * res1;

  g_return_if_fail (toplayer != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (apar != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (div != NULL);
  g_return_if_fail (divn != NULL);
  g_return_if_fail (G > 0.);

  fp.pn = p;
  fp.div = div;
  fp.dia = gfs_temporary_variable (toplayer);
  res1 = res ? res : gfs_temporary_variable (toplayer);
  fp.divn = divn;
  fp.dt = apar->dt;
  fp.G = G;

  /* compute MAC divergence */
  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) normal_divergence, &fp);
  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
  			    (FttCellTraverseFunc) scale_divergence_helmoltz, &fp);
  
  /* solve for pressure */
  par->dimension = 2;
  par->poisson_solve (toplayer, par, p, fp.div, res1, fp.dia, apar->dt);
#if !FTT_2D
  gfs_domain_traverse_merged (toplayer, (GfsMergedTraverseFunc) merge_pressures, p);
#endif

  if (!res)
    gts_object_destroy (GTS_OBJECT (res1));
  gts_object_destroy (GTS_OBJECT (fp.dia));
}

#if FTT_2D

static void normal_velocities (GfsDomain * domain, GfsVariable ** u)
{
  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, u);
}

static void ocean_run (GfsSimulation * sim)
{
  GfsVariable * p, * div, * res = NULL;
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
#if 1
  GfsVariable * H = gfs_variable_from_name (domain->variables, "H");
  g_assert (H);
  GfsFunction * fH = gfs_function_new_from_variable (gfs_function_class (), H);
#else
  /* non-linear free surface */
  GtsFile * fp = gts_file_new_from_string ("(H + P)"); /* fixme: should be H + P/g */
  GfsFunction * fH = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (fH, domain, fp);
  g_assert (fp->type != GTS_ERROR);
  gts_file_destroy (fp);
#endif

  div = gfs_temporary_variable (domain);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    GfsVariable * g[2];
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    gfs_simulation_set_timestep (sim);

    normal_velocities (domain, gfs_domain_velocity (domain));
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) gfs_normal_divergence_2D, div);

    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_poisson_coefficients (domain, fH, TRUE, TRUE, TRUE);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., 
					    sim->approx_projection_params.weighted);
    gfs_centered_velocity_advection_diffusion (domain, 2,
					       &sim->advection_params,
					       g, g,
					       sim->physical_params.alpha);
    gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
    gfs_correct_centered_velocities (domain, 2, g, -sim->advection_params.dt/2.);
    gts_object_destroy (GTS_OBJECT (g[0]));
    gts_object_destroy (GTS_OBJECT (g[1]));

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_domain_timer_start (domain, "free_surface_pressure");
    GfsVariable * divn = gfs_temporary_variable (domain);
    normal_velocities (domain, gfs_domain_velocity (domain));
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) gfs_normal_divergence_2D, divn);
    gfs_poisson_coefficients (domain, fH, TRUE, TRUE, TRUE);
    gfs_free_surface_pressure (domain, &sim->approx_projection_params, &sim->advection_params,
			       p, divn, div, res, sim->physical_params.g);
    gts_object_destroy (GTS_OBJECT (divn));
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, sim->advection_params.dt/2., 
					    sim->approx_projection_params.weighted);
    gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt/2.);
    gts_object_destroy (GTS_OBJECT (g[0]));
    gts_object_destroy (GTS_OBJECT (g[1]));
    gfs_domain_timer_stop (domain, "free_surface_pressure");

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events),
			 (GtsFunc) gts_object_destroy, NULL);

  gts_object_destroy (GTS_OBJECT (div));
  gts_object_destroy (GTS_OBJECT (fH));
}

static void gfs_ocean_class_init (GfsSimulationClass * klass)
{
  klass->run = ocean_run;
}

static void gfs_ocean_init (GfsOcean * object)
{
  GFS_SIMULATION (object)->approx_projection_params.weighted = 1;
}

GfsSimulationClass * gfs_ocean_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_ocean_info = {
      "GfsOcean",
      sizeof (GfsOcean),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_ocean_class_init,
      (GtsObjectInitFunc) gfs_ocean_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_ocean_info);
  }

  return klass;
}

#else /* 3D */

#define MAC 0

static void ocean_destroy (GtsObject * object)
{
  guint i;
  GPtrArray * layer = GFS_OCEAN (object)->layer;

  for (i = 0; i < layer->len; i++) {
    GfsDomain * d = g_ptr_array_index (layer, i);
    d->allocated = g_array_new (FALSE, TRUE, sizeof (gboolean));
    gts_object_destroy (GTS_OBJECT (d));
  }
  g_ptr_array_free (layer, TRUE);

  (* GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class->destroy) (object);  
}

#define MAXLEVEL 16

static void ocean_read (GtsObject ** object, GtsFile * fp)
{
  /* fixme: lambda.z cannot be changed */
  GfsSimulation * sim = GFS_SIMULATION (*object);
  GFS_DOMAIN (sim)->lambda.z = 1./(1 << MAXLEVEL);

  (* GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class->read) (object, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_DOMAIN (*object)->refpos.z = -0.5;
  g_assert (GFS_DOMAIN (sim)->lambda.z == 1./(1 << MAXLEVEL));
  sim->physical_params.g /= sim->physical_params.L/* *GFS_DOMAIN (sim)->lambda.z*/;
  GfsVariable * H = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "H");
  g_assert (H);
  H->units = 1.;
  GFS_DOMAIN (sim)->lambda.z = 1./(1 << MAXLEVEL);
}

static void ocean_write (GtsObject * object, FILE * fp)
{
  FttVector * lambda = &GFS_DOMAIN (object)->lambda;
  GfsPhysicalParams * p = &GFS_SIMULATION (object)->physical_params;
  gdouble g = p->g;

  lambda->z *= 1 << MAXLEVEL;
  p->g *= p->L*lambda->z;
  (* GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class->write) (object, fp);
  lambda->z /= 1 << MAXLEVEL;
  p->g = g;
}

static void new_layer (GfsOcean * ocean)
{
  GfsDomain * domain = GFS_DOMAIN (ocean);
  GfsDomain * d = GFS_DOMAIN (gts_object_new (GTS_OBJECT_CLASS (gfs_domain_class ())));
  
  d->rootlevel = domain->rootlevel;
  d->refpos = domain->refpos;
  d->lambda = domain->lambda;
  g_array_free (d->allocated, TRUE);
  d->allocated = domain->allocated;
  g_ptr_array_add (ocean->layer, d);
}

static void add_layer (GfsBox * box, GfsDomain * domain)
{
  if (box->neighbor[FTT_FRONT] == NULL || GFS_IS_BOUNDARY (box->neighbor[FTT_FRONT])) {
    GPtrArray * layer = GFS_OCEAN (domain)->layer;
    GtsObject * n;
    guint l = 0;

    gts_container_add (GTS_CONTAINER (g_ptr_array_index (layer, l++)), GTS_CONTAINEE (box));
    n = box->neighbor[FTT_BACK];
    while (GFS_IS_BOX (n)) {
      if (l == layer->len)
	new_layer (GFS_OCEAN (domain));
      gts_container_add (GTS_CONTAINER (g_ptr_array_index (layer, l++)), GTS_CONTAINEE (n));
      n = GFS_BOX (n)->neighbor[FTT_BACK];
    }
  }
}

static void ocean_post_read (GfsDomain * domain, GtsFile * fp)
{
  (* GFS_DOMAIN_CLASS (GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class)->post_read) 
    (domain, fp);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) add_layer, domain);
  g_assert (GFS_OCEAN (domain)->layer->len > 0);
  GFS_OCEAN (domain)->toplayer = g_ptr_array_index (GFS_OCEAN (domain)->layer, 0);
}

static void compute_w (FttCell * c, GfsVariable * W)
{
  FttCell * n;
  guint level = ftt_cell_level (c);
  gdouble wf = 0., w = 0.;

  while ((n = ftt_cell_neighbor (c, FTT_BACK)))
    c = n;
  while (c) {
    GfsStateVector * s = GFS_STATE (c);

    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    s->f[FTT_BACK].un = w;
    wf += (s->f[FTT_LEFT].v*s->f[FTT_LEFT].un - s->f[FTT_RIGHT].v*s->f[FTT_RIGHT].un +
    	   s->f[FTT_BOTTOM].v*s->f[FTT_BOTTOM].un - s->f[FTT_TOP].v*s->f[FTT_TOP].un);
    if (GFS_IS_MIXED (c))
      s->f[FTT_FRONT].un = w = GFS_STATE (c)->solid->s[FTT_FRONT] > 0. ? 
	wf/GFS_STATE (c)->solid->s[FTT_FRONT] : 0.;
    else
      s->f[FTT_FRONT].un = w = wf;
    GFS_VALUE (c, W) = (s->f[FTT_BACK].un + s->f[FTT_FRONT].un)/2.;
    c = ftt_cell_neighbor (c, FTT_FRONT);
  }
}

static void compute_div (FttCell * c, GfsVariable * W)
{
  guint level = ftt_cell_level (c);
  gdouble wf = 0., size = ftt_cell_size (c);

  g_assert (level <= MAXLEVEL);
  size *= 1 << (MAXLEVEL - level);

  while (c) {
    GfsStateVector * s = GFS_STATE (c);
    GfsSolidVector * solid = s->solid;

    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    if (solid)
      wf += (solid->s[FTT_RIGHT]*s->f[FTT_RIGHT].un - solid->s[FTT_LEFT]*s->f[FTT_LEFT].un +
	     solid->s[FTT_TOP]*s->f[FTT_TOP].un - solid->s[FTT_BOTTOM]*s->f[FTT_BOTTOM].un);
    else
      wf += (s->f[FTT_RIGHT].un - s->f[FTT_LEFT].un +
	     s->f[FTT_TOP].un - s->f[FTT_BOTTOM].un);
    GFS_VALUE (c, W) = wf*size;
    c = ftt_cell_neighbor (c, FTT_BACK);
  }
}

/* fixme: this is ok for one layer but what about several? */
static gdouble height (FttCell * cell)
{
  if (!GFS_IS_MIXED (cell))
    return 1.;
  gdouble f = GFS_STATE (cell)->solid->s[FTT_FRONT];
  if (f == 0.)
    return 0.;
  guint level = ftt_cell_level (cell);
  g_assert (level <= MAXLEVEL);
  return GFS_STATE (cell)->solid->a/f*(1 << (MAXLEVEL - level));
}

static void compute_H (FttCell * cell, GfsVariable * H)
{
  GFS_VALUE (cell, H) = height (cell);
}

static void face_interpolated_normal_velocity (const FttCellFace * face, GfsVariable ** v)
{
  gdouble u;

  g_return_if_fail (face != NULL);
  g_return_if_fail (v != NULL);

  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  guint i = v[face->d/2]->i;
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    u = (GFS_VALUEI (face->cell, i) + GFS_VALUEI (face->neighbor, i))/2.; 
    break;
  case FTT_FINE_COARSE: {
    gdouble w1 = height (face->cell), w2 = height (face->neighbor);
    g_assert (w1 + w2);
    w1 = 2.*w1/(w1 + w2);
    u = w1*gfs_face_interpolated_value (face, i) + (1. - w1)*GFS_VALUEI (face->neighbor, i);
    break;
  }
  default:
     g_assert_not_reached ();
  }

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) = u;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) = u;
    break;
  case FTT_FINE_COARSE:
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) += 
      u*GFS_FACE_FRACTION_LEFT (face)/(GFS_FACE_FRACTION_RIGHT (face)*
				       FTT_CELLS_DIRECTION (face->d));
    break;
  default:
    g_assert_not_reached ();
  }
}

static void depth_integrated_divergence (GfsDomain * domain, GfsVariable * div)
{
  /* compute MAC velocities from centered velocities */
#if !MAC
  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) face_interpolated_normal_velocity,
			    gfs_domain_velocity (domain));
#endif
  /* barotropic divergence */
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				     (FttCellTraverseFunc) compute_div, div);
}

static void compute_coeff (FttCell * c)
{
  guint level = ftt_cell_level (c);
  gdouble wf[FTT_NEIGHBORS_2D] = {0.,0.,0.,0.}, size = 1.;

  g_assert (level <= MAXLEVEL);
  size = 1 << (MAXLEVEL - level);

  while (c) {
    GfsStateVector * s = GFS_STATE (c);
    FttDirection d;

    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    for (d = 0; d < FTT_NEIGHBORS_2D; d++) {
      wf[d] += s->f[d].v*size;
      s->f[d].v = wf[d];
    }
    c = ftt_cell_neighbor (c, FTT_BACK);
  }
}

static void face_coeff_from_below (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  guint neighbors = 0;

  for (d = 0; d < FTT_NEIGHBORS_2D; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
	f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= 2;

    FttCell * neighbor;
    if (f[d].v > 0. && (neighbor = ftt_cell_neighbor (cell, d)) && !GFS_CELL_IS_BOUNDARY (neighbor))
      neighbors++;
  }

  if (neighbors == 1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      f[d].v = 0.;
}

static void depth_integrated_coefficients (GfsDomain * domain)
{
  gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				     (FttCellTraverseFunc) compute_coeff, NULL);
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				     FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				     (FttCellTraverseFunc) face_coeff_from_below, NULL);
}

static void ocean_run (GfsSimulation * sim)
{
  GfsVariable * p, * div, * H, * res = NULL;
  GfsDomain * domain, * toplayer;
  GSList * i;

  domain = GFS_DOMAIN (sim);
  toplayer = GFS_OCEAN (sim)->toplayer;
  gfs_clock_start (toplayer->timer);
  
  gfs_simulation_refine (sim);

  H = gfs_variable_from_name (domain->variables, "H");
  g_assert (H);

  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) compute_H, H);

  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_init, sim);

  gfs_set_merged (domain);
  i = domain->variables;
  while (i) {
    gfs_event_init (i->data, sim);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, i->data);
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);

  div = gfs_temporary_variable (domain);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    GfsVariable * g[2];
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    gfs_simulation_set_timestep (sim);

    depth_integrated_divergence (domain, div);

    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p);

    /* baroclinic terms */
#if !MAC
    gfs_predicted_face_velocities (domain, 2, &sim->advection_params);

    gfs_domain_timer_start (domain, "correct_normal_velocities");
    gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, sim->advection_params.dt/2.,
					    sim->approx_projection_params.weighted);
    gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				       (FttCellTraverseFunc) compute_w, 
				       gfs_variable_from_name (domain->variables, "W"));
    gfs_domain_timer_stop (domain, "correct_normal_velocities");

    i = domain->variables;
    while (i) {
      if (GFS_IS_VARIABLE_TRACER_VOF (i->data)) {
	GfsVariableTracer * t = i->data;

	t->advection.dt = sim->advection_params.dt;
	gfs_tracer_vof_advection (domain, &t->advection);
	gfs_domain_variable_centered_sources (domain, i->data, i->data, t->advection.dt);
      }
      else if (GFS_IS_VARIABLE_TRACER (i->data)) {
	GfsVariableTracer * t = i->data;
	
	t->advection.dt = sim->advection_params.dt;
	gfs_tracer_advection_diffusion (domain, &t->advection, NULL);
      }
      i = i->next;
    }

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_centered_velocity_advection_diffusion (domain, 2,
					       &sim->advection_params,
					       g, g,
					       sim->physical_params.alpha);
    gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
    gfs_correct_centered_velocities (domain, 2, g, -sim->advection_params.dt/2.);
#else
    gfs_poisson_coefficients (domain, NULL);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, sim->advection_params.dt/2.,
					    sim->approx_projection_params.weighted);
    gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt/2.);
#endif
    gts_object_destroy (GTS_OBJECT (g[0]));
    gts_object_destroy (GTS_OBJECT (g[1]));

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_domain_timer_start (domain, "free_surface_pressure");
    GfsVariable * divn = gfs_temporary_variable (domain);
    depth_integrated_divergence (domain, divn);
    depth_integrated_coefficients (domain);
    gfs_free_surface_pressure (toplayer, &sim->approx_projection_params, &sim->advection_params,
			       p, divn, div, res, 
			       sim->physical_params.g/GFS_OCEAN (domain)->layer->len);
    gts_object_destroy (GTS_OBJECT (divn));

    gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, sim->advection_params.dt/2.,
					    sim->approx_projection_params.weighted);
    gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt/2.);
    gts_object_destroy (GTS_OBJECT (g[0]));
    gts_object_destroy (GTS_OBJECT (g[1]));
    
    gfs_domain_timer_stop (domain, "free_surface_pressure");

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events),
			 (GtsFunc) gts_object_destroy, NULL);

  gts_object_destroy (GTS_OBJECT (div));

  gfs_clock_stop (toplayer->timer);
}

static void gfs_ocean_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = ocean_destroy;
  GTS_OBJECT_CLASS (klass)->read = ocean_read;
  GTS_OBJECT_CLASS (klass)->write = ocean_write;
  GFS_DOMAIN_CLASS (klass)->post_read = ocean_post_read;
  klass->run = ocean_run;
}

static void depth_coarse_fine (FttCell * parent, GfsVariable * H)
{
  FttCellChildren child;
  guint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      GFS_VALUE (child.c[i], H) = height (child.c[i]);
      if (GFS_VALUE (child.c[i], H) <= 0.)
	ftt_cell_destroy (child.c[i], (FttCellCleanupFunc) gfs_cell_cleanup, H->domain);
    }
}

static void depth_fine_coarse (FttCell * parent, GfsVariable * H)
{
  GFS_VALUE (parent, H) = height (parent);
}

static void gfs_ocean_init (GfsOcean * object)
{
  GfsVariable * H = gfs_domain_add_variable (GFS_DOMAIN (object), "H", "Depth");
  H->coarse_fine = depth_coarse_fine;
  H->fine_coarse = depth_fine_coarse;
  GFS_SIMULATION (object)->approx_projection_params.weighted = 1;
  object->layer = g_ptr_array_new ();
  new_layer (object);
}

GfsSimulationClass * gfs_ocean_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_ocean_info = {
      "GfsOcean",
      sizeof (GfsOcean),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_ocean_class_init,
      (GtsObjectInitFunc) gfs_ocean_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_ocean_info);
  }

  return klass;
}

static void hydrostatic_pressure (FttCell * cell, gpointer * data)
{
  GfsVariable * vp = data[0];
  GfsVariable * rho = data[1];
  gdouble * g = data[2];
  gdouble r = GFS_VALUE (cell, rho), p = (*g)*r/2., r1;
  FttCellFace f;
  
  GFS_VALUE (cell, vp) = p;
  f.cell = cell;
  f.d = FTT_BACK;
  f.neighbor = ftt_cell_neighbor (f.cell, f.d);
  while (f.neighbor) {
    g_assert (ftt_face_type (&f) == FTT_FINE_FINE);
    r1 = gfs_face_interpolated_value_generic (&f, rho);
    /* g_assert (r1 >= r); */
    r = r1;
    GFS_VALUE (f.neighbor, vp) = p = p + (*g)*r;
    f.cell = f.neighbor;
    f.neighbor = ftt_cell_neighbor (f.cell, f.d);
  }
}

/**
 * gfs_hydrostatic_pressure:
 * @domain: a #GfsDomain.
 * @p: the hydrostatic pressure.
 * @rho: the density.
 * @g: the acceleration.
 *
 * Computes the hydrostatic pressure @p in @domain using the density
 * @rho.
 */
void gfs_hydrostatic_pressure (GfsDomain * domain,
			       GfsVariable * p,
			       GfsVariable * rho,
			       gdouble g)
{
  gpointer data[3];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (rho != NULL);
  g_return_if_fail (g >= 0.);

  g /= GFS_OCEAN (domain)->layer->len;
  data[0] = p;
  data[1] = rho;
  data[2] = &g;
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				     (FttCellTraverseFunc) hydrostatic_pressure, data);
}

/** \endobject{GfsOcean} */

/* GfsSourceHydrostatic: Object */

static void gfs_source_hydrostatic_destroy (GtsObject * o)
{
  if (GFS_SOURCE_HYDROSTATIC (o)->ph1)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_HYDROSTATIC (o)->ph1));

  (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->destroy) (o);
}


static void gfs_source_hydrostatic_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSourceHydrostatic * sh;

  if (GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  sh = GFS_SOURCE_HYDROSTATIC (*o);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (rho)");
    return;
  }
  sh->rho = gfs_variable_from_name (domain->variables, fp->token->str);
  if (sh->rho == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (ph)");
    return;
  }
  if (!(sh->ph = gfs_domain_get_or_add_variable (domain, fp->token->str, "Hydrostatic pressure"))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  sh->ph1 = gfs_temporary_variable (domain);
}

static void gfs_source_hydrostatic_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %s",
	   GFS_SOURCE_HYDROSTATIC (o)->rho->name, 
	   GFS_SOURCE_HYDROSTATIC (o)->ph->name);
}

static gdouble gfs_source_hydrostatic_mac_value (GfsSourceGeneric * s,
						 FttCell * cell,
						 GfsVariable * v)
{
  return - gfs_center_gradient (cell, v->component,
				GFS_SOURCE_HYDROSTATIC (s)->ph1->i)/ftt_cell_size (cell);
}

static gdouble gfs_source_hydrostatic_centered_value (GfsSourceGeneric * s,
						      FttCell * cell,
						      GfsVariable * v)
{
  GfsSourceHydrostatic * b = GFS_SOURCE_HYDROSTATIC (s);

  return - (gfs_center_gradient (cell, v->component, b->ph->i) + 
	    gfs_center_gradient (cell, v->component, b->ph1->i))/(2.*ftt_cell_size (cell));
}

static void copy_ph (FttCell * cell, GfsSourceHydrostatic * s)
{
  GFS_VALUE (cell, s->ph1) = GFS_VALUE (cell, s->ph);
}

static gboolean gfs_source_hydrostatic_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) 
      (event, sim)) {
    GfsSourceHydrostatic * s = GFS_SOURCE_HYDROSTATIC (event);

    if (s->not_first) {
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) copy_ph, s);
      gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph1);
    }
    else {
      gfs_hydrostatic_pressure (GFS_DOMAIN (sim), s->ph1, s->rho, sim->physical_params.g);
      gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph1);
      s->not_first = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_source_hydrostatic_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsSourceHydrostatic * s = GFS_SOURCE_HYDROSTATIC (event);

  gfs_hydrostatic_pressure (GFS_DOMAIN (sim), s->ph, s->rho, sim->physical_params.g);
  gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph);
}

static void gfs_source_hydrostatic_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_hydrostatic_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_hydrostatic_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_hydrostatic_write;

  GFS_EVENT_CLASS (klass)->event = gfs_source_hydrostatic_event;
  GFS_EVENT_CLASS (klass)->event_half = gfs_source_hydrostatic_event_half;
}

static void gfs_source_hydrostatic_init (GfsSourceGeneric * s)
{
  s->mac_value = gfs_source_hydrostatic_mac_value;
  s->centered_value = gfs_source_hydrostatic_centered_value;
}

GfsSourceGenericClass * gfs_source_hydrostatic_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_hydrostatic_info = {
      "GfsSourceHydrostatic",
      sizeof (GfsSourceHydrostatic),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_hydrostatic_class_init,
      (GtsObjectInitFunc) gfs_source_hydrostatic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_hydrostatic_info);
  }

  return klass;
}

#endif /* 3D */

/**
 *
 * \beginobject{GfsSourceFriction}
 */

static void gfs_source_friction_destroy (GtsObject * o)
{
  FttComponent c;

  for (c = 0; c <  FTT_DIMENSION; c++)
    if (GFS_SOURCE_FRICTION (o)->u[c])
      gts_object_destroy (GTS_OBJECT (GFS_SOURCE_FRICTION (o)->u[c]));

  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->destroy) (o);
}

static void gfs_source_friction_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsVariable h)");
    return;
  }
  GFS_SOURCE_FRICTION (*o)->h = gfs_variable_from_name (domain->variables, fp->token->str);
  if (GFS_SOURCE_FRICTION (*o)->h == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  /* fixme: units? */
  GFS_SOURCE_FRICTION (*o)->f = gfs_read_constant (fp, domain);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c <  FTT_DIMENSION; c++)
    GFS_SOURCE_FRICTION (*o)->u[c] = gfs_temporary_variable (domain);
}

static void gfs_source_friction_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", GFS_SOURCE_FRICTION (o)->h->name, GFS_SOURCE_FRICTION (o)->f);
}

static gdouble gfs_source_friction_saved_value (GfsSourceGeneric * s, 
						FttCell * cell, 
						GfsVariable * v)
{
  gdouble H = GFS_VALUE (cell, GFS_SOURCE_FRICTION (s)->h);

  g_assert (H > 0.);
  return - GFS_SOURCE_FRICTION (s)->f*
    GFS_VALUE (cell, GFS_SOURCE_FRICTION (s)->u[v->component])/H;
}

static void save_velocity (FttCell * cell, GfsSourceFriction * s)
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, s->u[c]) = GFS_VALUE (cell, GFS_SOURCE_VELOCITY (s)->v[c]);
}

static gboolean gfs_source_friction_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event)
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_velocity, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_source_friction_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_friction_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_friction_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_friction_write;
  GFS_EVENT_CLASS (klass)->event = gfs_source_friction_event;
}

static void gfs_source_friction_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = gfs_source_friction_saved_value;
}

GfsSourceGenericClass * gfs_source_friction_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_friction_info = {
      "GfsSourceFriction",
      sizeof (GfsSourceFriction),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_friction_class_init,
      (GtsObjectInitFunc) gfs_source_friction_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_friction_info);
  }

  return klass;
}

/** \endobject{GfsSourceFriction} */

/**
 * Flather boundary conditions.
 * \beginobject{GfsBcFlather}
 */

/* Also check whether modules/tide.mod needs upgrading when modyifing this class */

static void bc_flather_write (GtsObject * o, FILE * fp)
{
  GfsBcFlather * bc = GFS_BC_FLATHER (o);

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %s", bc->h->name, bc->p->name);
  if (bc->val)
    gfs_function_write (bc->val, fp);
}

static void set_gradient_boundary (FttCell * cell)
{
  cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
}

static void bc_flather_read (GtsObject ** o, GtsFile * fp)
{
  GfsBcFlather * bc = GFS_BC_FLATHER (*o);
  GfsDomain * domain = gfs_box_domain (GFS_BC (bc)->b->box);

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->read) (o, fp);

  gfs_function_set_units (GFS_BC_VALUE (bc)->val, 1.);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (h)");
    return;
  }
  bc->h = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->h == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (p)");
    return;
  }
  bc->p = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->p == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (bc->val == NULL)
    bc->val = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (bc->val, gfs_box_domain (GFS_BC (bc)->b->box), fp);
  gfs_function_set_units (bc->val, 1.);

  ftt_cell_traverse (GFS_BC (bc)->b->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) set_gradient_boundary, NULL);
}

static void bc_flather_destroy (GtsObject * o)
{
  if (GFS_BC_FLATHER (o)->val)
    gts_object_destroy (GTS_OBJECT (GFS_BC_FLATHER (o)->val));

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->destroy) (o);
}

static gdouble flather_value (FttCellFace * f, GfsBc * b)
{
  /* fixme: this will not work for multilayer domains */
  guint d, nb = 0;
  FttCellNeighbors n;
  gdouble H;

  ftt_cell_neighbors (f->neighbor, &n);
  for (d = 0; d < FTT_NEIGHBORS_2D; d++)
    if (n.c[d] != NULL && GFS_CELL_IS_BOUNDARY(n.c[d]) && nb++ > 0)
      /* if the boundary cell is bounded by more than one boundary -> no flux */
      return 0.;

  GfsSimulation * sim = GFS_SIMULATION (gfs_box_domain (b->b->box));
  H = gfs_face_interpolated_value (f, GFS_BC_FLATHER (b)->h->i);
  if (H > 10./sim->physical_params.L) { /* fixme: this a bit crappy */
    gdouble cg = sqrt (sim->physical_params.g*H);
    /* non-dimensional pressure at the boundary */
    gdouble lz = GFS_DOMAIN (sim)->lambda.z;
#if !FTT_2D
    lz *= 1 << MAXLEVEL;
#endif
    gdouble pb = gfs_function_face_value (GFS_BC_FLATHER (b)->val, f)*sim->physical_params.g*lz;
    
    return gfs_function_face_value (GFS_BC_VALUE (b)->val, f) +
      (FTT_FACE_DIRECT (f) ? -1. : 1.)*
      (GFS_VALUE (f->neighbor, GFS_BC_FLATHER (b)->p) - pb)*
      cg/sim->physical_params.g
#if !FTT_2D
      /H
#endif
      ;
  }
  else
    return 0.;
}

static void flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_VALUE (f->cell, b->v) = 2.*flather_value (f, b) - GFS_VALUE (f->neighbor, b->v);
}

static void homogeneous_flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_VALUE (f->cell, b->v) = - GFS_VALUE (f->neighbor, b->v);
}

static void face_flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_STATE (f->cell)->f[f->d].v = flather_value (f, b);
}

static void gfs_bc_flather_class_init (GtsObjectClass * klass)
{
  klass->write   = bc_flather_write;
  klass->read    = bc_flather_read;
  klass->destroy = bc_flather_destroy;
}

static void gfs_bc_flather_init (GfsBc * object)
{
  object->bc =             (FttFaceTraverseFunc) flather;
  object->homogeneous_bc = (FttFaceTraverseFunc) homogeneous_flather;
  object->face_bc =        (FttFaceTraverseFunc) face_flather;
}

GfsBcClass * gfs_bc_flather_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_flather_info = {
      "GfsBcFlather",
      sizeof (GfsBcFlather),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_flather_class_init,
      (GtsObjectInitFunc) gfs_bc_flather_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_flather_info);
  }

  return klass;
}

/** \endobject{GfsBcFlather} */
