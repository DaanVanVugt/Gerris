/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
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
 * \brief Timestepping.
 */

#include <math.h>
#include <stdlib.h>

#include "timestep.h"
#include "source.h"
#include "solid.h"
#include "tension.h"

typedef struct {
  GfsVariable ** g;
  guint dimension;
} GradientData;

static void reset_cell_gradients (FttCell * cell, GradientData * p)
{
  FttComponent c;
  for (c = 0; c < p->dimension; c++)
    GFS_VALUE (cell, p->g[c]) = 0.;
}

/**
 * gfs_reset_gradients:
 * @domain: a #GfsDomain.
 * @g: an array of @dimension #GfsVariables.
 * @dimension: the number of dimension.
 *
 * Resets the gradient vector @g.
 */
void gfs_reset_gradients (GfsDomain * domain, guint dimension, GfsVariable ** g)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (g != NULL);

  GradientData p = { g, dimension };
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) reset_cell_gradients, &p);
}

static void scale_cell_gradients (FttCell * cell, gpointer * data)
{
  GfsVariable ** g = data[0];
  guint * dimension = data[1];
  FttComponent c;

  /* fixme: mapping??? */
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;

    for (c = 0; c < *dimension; c++)
      if (s->s[2*c] + s->s[2*c + 1] > 0.)
	GFS_VALUE (cell, g[c]) /= s->s[2*c] + s->s[2*c + 1];
      else
	g_assert (GFS_VALUE (cell, g[c]) == 0.);
  }
  else {
    FttCellNeighbors n;
    
    ftt_cell_neighbors (cell, &n);
    for (c = 0; c < *dimension; c++) {
      FttCell * c1 = n.c[2*c], * c2 = n.c[2*c + 1];
      
      if (c1 && c2 && !GFS_CELL_IS_GRADIENT_BOUNDARY (c1) && !GFS_CELL_IS_GRADIENT_BOUNDARY (c2))
	GFS_VALUE (cell, g[c]) /= 2.;
    }
  }
}

/**
 * gfs_scale_gradients:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions.
 * @g: the components of the gradient.
 *
 * Scales the gradient accumulated in @g (typically using
 * gfs_correct_normal_velocities()).
 */
void gfs_scale_gradients (GfsDomain * domain, guint dimension, GfsVariable ** g)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (g != NULL);

  gpointer data[2];
  data[0] = g;
  data[1] = &dimension;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) scale_cell_gradients, data);
  FttComponent c;
  for (c = 0; c < dimension; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, g[c]);
}

typedef struct {
  GfsVariable * p, ** gv;
  gdouble dt;
} CorrectPar;

static void correct_normal_velocity (FttCellFace * face,
				     CorrectPar * par)
{
  GfsGradient g;
  gdouble dp, f;

  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  gfs_face_weighted_gradient (face, &g, par->p->i, -1);
  dp = (g.b - g.a*GFS_VALUE (face->cell, par->p))/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    dp = - dp;
  f = gfs_domain_face_fraction (par->p->domain, face);
  if (f > 0.)
    dp /= f;

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= dp*par->dt;
  if (par->gv)
    GFS_VALUE (face->cell, par->gv[face->d/2]) += dp*GFS_FACE_FRACTION_LEFT (face);

  if (ftt_face_type (face) == FTT_FINE_COARSE)
    dp *= GFS_FACE_FRACTION_LEFT (face)/(GFS_FACE_FRACTION_RIGHT (face)*FTT_CELLS/2);
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp*par->dt;
  if (par->gv)
    GFS_VALUE (face->neighbor, par->gv[face->d/2]) += dp*GFS_FACE_FRACTION_RIGHT (face);
}

/**
 * gfs_correct_normal_velocities:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @p: the pressure field.
 * @g: where to store the pressure gradient or %NULL.
 * @dt: the timestep.
 *
 * Corrects the normal velocity field of @domain using @p and and @dt.
 *
 * Assumes that the Poisson weighting coefficients have already been
 * computed using gfs_poisson_coefficients().
 *
 * Also fills the @g variables (if not %NULL) with the centered
 * gradient of @p.
 */
void gfs_correct_normal_velocities (GfsDomain * domain,
				    guint dimension,
				    GfsVariable * p,
				    GfsVariable ** g,
				    gdouble dt)
{
  CorrectPar par;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);

  par.p = p;
  par.gv = g;
  par.dt = dt;
  gfs_domain_face_traverse (domain, dimension == 2 ? FTT_XY : FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) correct_normal_velocity, &par);
}

static void scale_divergence (FttCell * cell, gpointer * data)
{
  GfsVariable * div = data[0];
  gdouble * dt = data[1];

  GFS_VALUE (cell, div) /= *dt;
}

typedef struct {
  GfsSourceGeneric * s;
  GfsVariable * v, ** g;
  FttComponent c;
  gdouble dt;
} FaceSource;

#define DIRICHLET_BOUNDARY (1 << GFS_FLAG_USER)

static void tag_dirichlet (FttCellFace * f)
{
  f->cell->flags |= DIRICHLET_BOUNDARY;
}

static void ignore_dirichlet_boundaries (GfsBox * box, FaceSource * f)
{
  FttDirection d;
  for (d = 2*f->c; d <= 2*f->c + 1; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, f->v);
      if (GFS_IS_BC_DIRICHLET (bc))
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) tag_dirichlet, NULL);
    }
}

static void add_face_source (FttCellFace * face,
			     FaceSource * f)
{
  gdouble dp;
  FttComponent c;

  if (face->neighbor->flags & DIRICHLET_BOUNDARY) {
    face->neighbor->flags &= ~DIRICHLET_BOUNDARY;
    return;
  }

  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  c = face->d/2;
  dp = (* f->s->face_value) (f->s, face, f->v);
  GFS_FACE_NORMAL_VELOCITY_LEFT (face) += dp*f->dt;
  if (f->g)
    GFS_VALUE (face->cell, f->g[c]) -= dp*GFS_FACE_FRACTION_LEFT (face);

  if (ftt_face_type (face) == FTT_FINE_COARSE)
    dp *= GFS_FACE_FRACTION_LEFT (face)/(GFS_FACE_FRACTION_RIGHT (face)*FTT_CELLS/2);
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) += dp*f->dt;
  if (f->g)
    GFS_VALUE (face->neighbor, f->g[c]) -= dp*GFS_FACE_FRACTION_RIGHT (face);
}

/**
 * gfs_velocity_face_sources:
 * @domain: a #GfsDomain.
 * @u: the velocity vector.
 * @dt: the timestep.
 * @alpha: the specific volume.
 * @g: the gradient vector.
 *
 * Add source terms on the velocity component of each cell faces.
 */
void gfs_velocity_face_sources (GfsDomain * domain,
				GfsVariable ** u,
				gdouble dt,
				GfsFunction * alpha,
				GfsVariable ** g)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);

  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (u[c]->sources) {
      GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;
      
      while (i) {
	GfsSourceGeneric * s = i->data;
	if (s->face_value) {
	  FaceSource f = { s, u[c], g, c, dt };
	  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) ignore_dirichlet_boundaries, &f);
	  gfs_domain_face_traverse (domain, c,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) add_face_source, &f);
	}	  
	i = i->next;
      }
    }
  if (u[0]->sources) {
    GSList * i = GTS_SLIST_CONTAINER (u[0]->sources)->items;
    
    while (i) {
      if (GFS_IS_SOURCE_TENSION (i->data)) {
	GfsSourceTension * s = i->data;
	gfs_source_tension_coefficients (s, domain, alpha);
	gfs_correct_normal_velocities (domain, FTT_DIMENSION,
				       GFS_SOURCE_TENSION_GENERIC (s)->c,
				       g, dt);
      }
      i = i->next;
    }
  }
}

/**
 * gfs_update_gradients:
 * @domain: a #GfsDomain.
 * @p: the pressure.
 * @alpha: the Poisson equation gradient weight.
 * @g: where to store the pressure gradient.
 *
 * Updates the gradients in @g using @p and @alpha.
 */
void gfs_update_gradients (GfsDomain * domain, 
			   GfsVariable * p,  
			   GfsFunction * alpha, 
			   GfsVariable ** g)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);

  /* Add face sources */
  gfs_reset_gradients (domain, FTT_DIMENSION, g);
  gfs_velocity_face_sources (domain, gfs_domain_velocity (domain), 0., alpha, g);
  /* Initialize face coefficients */
  gfs_poisson_coefficients (domain, alpha, TRUE, TRUE, TRUE);
  /* Add pressure gradient */
  gfs_correct_normal_velocities (domain, FTT_DIMENSION, p, g, 0.);
  gfs_scale_gradients (domain, FTT_DIMENSION, g);
}

typedef struct {
  GfsVariable * v, * div;
  gdouble L2;
} SourceVolumePar;   

static void cell_volume_source (FttCell * cell, SourceVolumePar * p)
{
  gdouble sum = 0.;
  GSList * i = GTS_SLIST_CONTAINER (p->v->sources)->items;

  while (i) {
    GfsSourceGeneric * s = i->data;
    if (s->centered_value)
      sum += (* s->centered_value) (s, cell, p->v);
    i = i->next;
  }

  /* div*h^2 (units source m^3/s/vol) */
  GFS_VALUE (cell, p->div) -= sum*p->L2*ftt_cell_volume (cell);
}

static void volume_sources (GfsDomain * domain, GfsVariable * p, GfsVariable * div) 
{
  SourceVolumePar par;
  par.div = div;
  par.v = p;
  par.L2 = GFS_SIMULATION (domain)->physical_params.L;
  par.L2 *= par.L2;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) cell_volume_source, &par);
}

static void mac_projection (GfsDomain * domain,
			    GfsMultilevelParams * par,
			    gdouble dt,
			    GfsVariable * p,
			    GfsFunction * alpha,
			    GfsVariable * res,
			    GfsVariable ** g,
			    void (* divergence_hook) (GfsDomain * domain, 
						      gdouble dt,
						      GfsVariable * div)
			    )
{
  /* Add face sources */
  gfs_reset_gradients (domain, FTT_DIMENSION, g);
  gfs_velocity_face_sources (domain, gfs_domain_velocity (domain), dt, alpha, g);

  GfsVariable * dia = gfs_temporary_variable (domain);
  GfsVariable * div = gfs_temporary_variable (domain);
  GfsVariable * res1 = res ? res : gfs_temporary_variable (domain);

  /* Initialize face coefficients */
  gfs_poisson_coefficients (domain, alpha, TRUE, TRUE, TRUE);

  /* hydrostatic pressure */
  GSList * i = domain->variables;
  while (i) {
    if (GFS_IS_HYDROSTATIC_PRESSURE (i->data))
      gfs_correct_normal_velocities (domain, FTT_DIMENSION, i->data, g, dt);
    i = i->next;
  }

  /* Initialize diagonal coefficient */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, dia);

  /* compute MAC divergence */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_normal_divergence, div);

  /* Divergence hook */
  if (divergence_hook)
    (* divergence_hook) (domain, dt, div);

  /* add volume sources (if any) */
  if (p->sources)
    volume_sources (domain, p, div);

  /* Scale divergence */
  gpointer data[2];
  data[0] = div;
  data[1] = &dt;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
  			    (FttCellTraverseFunc) scale_divergence, data);

#if 0
  {
    FILE * fp = fopen ("/tmp/mac", "wt");
    GfsNorm norm;

    gfs_write_mac_velocity (domain, 0.9, FTT_TRAVERSE_LEAFS, -1, NULL, fp);
    fclose (fp);
    norm = gfs_domain_norm_variable (domain, div, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "mac div before: %g %g %g\n",
	     norm.first, norm.second, norm.infty);
  }
#endif
  
  par->poisson_solve (domain, par, p, div, res1, dia, dt);

  gts_object_destroy (GTS_OBJECT (dia));
  gts_object_destroy (GTS_OBJECT (div));
  if (!res)
    gts_object_destroy (GTS_OBJECT (res1));

  gfs_correct_normal_velocities (domain, FTT_DIMENSION, p, g, dt);
  gfs_scale_gradients (domain, FTT_DIMENSION, g);
}

/**
 * gfs_mac_projection:
 * @domain: a #GfsDomain.
 * @par: the projection control parameters.
 * @dt: the timestep.
 * @p: the pressure.
 * @alpha: the Poisson equation gradient weight.
 * @g: where to store the pressure gradient.
 * @divergence_hook: a hook function or %NULL.
 *
 * Corrects the face-centered velocity field (MAC field) on the leaf
 * level of @domain using an exact (MAC) projection. The resulting
 * face-centered velocity field is (almost) exactly divergence
 * free. The (potential) pressure field is also obtained as a
 * by-product as well as its gradient at the center of the leaf cells
 * of the domain. The gradient is stored in newly-allocated @g[]
 * variables and is obtained by simple averaging from the face values
 * to the center. The newly-allocated @g[] variables should be freed
 * when not needed anymore.
 *
 * The @residual field of the @par projection parameters is set to the
 * norm of the residual after the projection. The @niter field of the
 * @par projection parameters is set to the number of iterations
 * performed to solve the Poisson equation. The other projection
 * parameters are not modified.
 */
void gfs_mac_projection (GfsDomain * domain,
			 GfsMultilevelParams * par,
			 gdouble dt,
			 GfsVariable * p,
			 GfsFunction * alpha,
			 GfsVariable ** g,
			 void (* divergence_hook) (GfsDomain * domain, 
						   gdouble dt,
						   GfsVariable * div)
			 )
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);

  gfs_domain_timer_start (domain, "mac_projection");

  mac_projection (domain, par, dt, p, alpha, NULL, g, divergence_hook);

  gfs_domain_timer_stop (domain, "mac_projection");

  if (par->residual.infty > par->tolerance)
    g_warning ("MAC projection: max residual %g > %g", par->residual.infty, par->tolerance);
}

static void correct (FttCell * cell, gpointer * data)
{
  FttComponent c;
  GfsVariable ** v = data[0];
  GfsVariable ** g = data[1];
  gdouble * dt = data[2];
  guint * dimension = data[3];

  for (c = 0; c < *dimension; c++)
    GFS_VALUE (cell, v[c]) -= GFS_VALUE (cell, g[c])*(*dt);
}

/**
 * gfs_correct_centered_velocities:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @g: the pressure gradient.
 * @dt: the timestep.
 *
 * Corrects the velocity field of @domain using the pressure gradient
 * stored in g[].
 *
 * The @g[] variables are freed by this function.
 */
void gfs_correct_centered_velocities (GfsDomain * domain,
				      guint dimension,
				      GfsVariable ** g,
				      gdouble dt)
{
  GfsVariable ** v;
  FttComponent c;
  gpointer data[4];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (g != NULL);

  data[0] = v = gfs_domain_velocity (domain);
  data[1] = g;
  data[2] = &dt;
  data[3] = &dimension;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) correct, data);
  for (c = 0; c < dimension; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v[c]);
}

/**
 * gfs_approximate_projection:
 * @domain: a #GfsDomain.
 * @par: the projection control parameters.
 * @dt: the timestep.
 * @p: the pressure.
 * @alpha: the Poisson equation gradient weight.
 * @res: the residual or %NULL.
 * @g: where to store the pressure gradient.
 * @divergence_hook: a hook function or %NULL.
 *
 * Corrects the centered velocity field on the leaf level of @domain
 * using an approximate projection. The resulting centered velocity
 * field is approximately divergence free. The (potential) pressure
 * field is also obtained as a by-product.
 *
 * The @residual field of the @par projection parameters is set to the
 * norm of the residual (on the MAC grid) after the projection. The
 * @niter field of the @par projection parameters is set to the number
 * of iterations performed to solve the Poisson equation. The other
 * projection parameters are not modified.
 *
 * The Poisson equation for the pressure is first solved on a MAC grid
 * where the MAC velocities are obtained from the centered velocities
 * by simple averaging. The resulting pressure gradients (defined on
 * the faces) are then averaged down on the center of the cells to
 * correct the centered velocity.  
 */
void gfs_approximate_projection (GfsDomain * domain,
				 GfsMultilevelParams * par,
				 gdouble dt,
				 GfsVariable * p,
				 GfsFunction * alpha,
				 GfsVariable * res,
				 GfsVariable ** g,
				 void (* divergence_hook) (GfsDomain * domain, 
							   gdouble dt,
							   GfsVariable * div)
				 )
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);

  gfs_domain_timer_start (domain, "approximate_projection");
  
  /* compute MAC velocities from centered velocities */
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, 
			    gfs_domain_velocity (domain));
  
  mac_projection (domain, par, dt, p, alpha, res, g, divergence_hook);

  gfs_correct_centered_velocities (domain, FTT_DIMENSION, g, dt);

  gfs_domain_timer_stop (domain, "approximate_projection");

  if (par->residual.infty > par->tolerance)
    g_warning ("approx projection: max residual %g > %g", par->residual.infty, par->tolerance);
}

static void save_face_values (FttCell * cell, GfsVariable * f[4])
{
  FttDirection d;
  for (d = 0; d < 4; d++)
    GFS_VALUE (cell, f[d]) = GFS_STATE (cell)->f[d].v;
}

static void face_values_init (FttCellTraverseFunc face_values, GfsAdvectionParams * par)
{
  if (par->scheme == GFS_GODUNOV &&
      par->v->component < FTT_DIMENSION && 
      par->v->domain->has_rotated_bc) {
    /* store normal and tangential face values for vector
       quantities. This is required to apply boundary conditions for
       rotated topologies (e.g. cubed sphere) */
    GfsVariable * v = par->v;
    GfsDomain * domain = v->domain;
    FttComponent c;

    for (c = 0; c < 2; c++) {
      FttDirection d;
      for (d = 0; d < 4; d++) {
	g_assert (v->vector[0]->face[c][d] == NULL);
	g_assert (v->vector[1]->face[c][d] == NULL);
	v->vector[0]->face[c][d] = v->vector[1]->face[c][d] = 
	  gfs_temporary_variable (domain);
      }
    }

    for (c = 0; c < 2; c++) {
      par->v = v->vector[c];
      gfs_domain_traverse_leaves (domain, face_values, par);
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) save_face_values, par->v->face[c]);
    }

    par->v = v;
  }
}

static void restore_face_values (FttCell * cell, GfsVariable * f[4])
{
  FttDirection d;
  for (d = 0; d < 4; d++)
    GFS_STATE (cell)->f[d].v = GFS_VALUE (cell, f[d]);
}

static void face_values_set (FttCellTraverseFunc face_values, GfsAdvectionParams * par)
{
  if (par->v->component < 2 && par->v->domain->has_rotated_bc)
    /* vector: restore stored face values */
    gfs_domain_traverse_leaves (par->v->domain, (FttCellTraverseFunc) restore_face_values, 
				par->v->face[par->v->component]);
  else
    /* scalar or z-component: compute face values */
    gfs_domain_traverse_leaves (par->v->domain, face_values, par);
  gfs_domain_face_bc (par->v->domain, FTT_XYZ, par->v);
}

static void face_values_free (GfsVariable * v)
{
  if (v->component < FTT_DIMENSION && v->domain->has_rotated_bc) {
    FttComponent c;
    for (c = 0; c < 2; c++) {
      FttDirection d;
      for (d = 0; d < 4; d++) {
	if (v->vector[0]->face[c][d])
	  gts_object_destroy (GTS_OBJECT (v->vector[0]->face[c][d]));
	v->vector[0]->face[c][d] = v->vector[1]->face[c][d] = NULL;
      }
    }
  }
}

/**
 * gfs_predicted_face_velocities:
 * @domain: a #GfsDomain.
 * @d: the number of dimensions (2 or 3).
 * @par: the advection parameters.
 *
 * Fills the face (MAC) normal velocities of each leaf cell of @domain
 * with the predicted values at time t + dt/2 using a godunov type
 * advection scheme.  
 */
void gfs_predicted_face_velocities (GfsDomain * domain,
				    guint d,
				    GfsAdvectionParams * par)
{
  FttComponent c;
  FttCellTraverseFunc face_values;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  gfs_domain_timer_start (domain, "predicted_face_velocities");

  gfs_domain_face_traverse (domain, d == 2 ? FTT_XY : FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  par->u = gfs_domain_velocity (domain);
  par->use_centered_velocity = TRUE;
  if (par->scheme == GFS_NONE) {
    face_values = (FttCellTraverseFunc) gfs_cell_non_advected_face_values;
    par->upwinding = GFS_NO_UPWINDING;
  }
  else {
    face_values = (FttCellTraverseFunc) gfs_cell_advected_face_values;
    par->upwinding = GFS_CENTERED_UPWINDING;
  }
  par->v = par->u[0];
  face_values_init (face_values, par);
  for (c = 0; c < d; c++) {
    par->v = par->u[c];
    face_values_set (face_values, par);
    gfs_domain_face_traverse (domain, c,
    			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_advected_normal_velocity, par);
  }
  face_values_free (par->u[0]);
  gfs_domain_timer_stop (domain, "predicted_face_velocities");
}

/**
 * gfs_diffusion:
 * @domain: a #GfsDomain.
 * @par: the multilevel parameters.
 * @v: a #GfsVariable.
 * @rhs: the right-hand side.
 * @rhoc: the mass.
 * @metric: the metric term.
 *
 * Solves a diffusion equation for variable @v using a Crank-Nicholson
 * scheme with multilevel relaxations.
 *
 * Diffusion coefficients must have been set using
 * gfs_diffusion_coefficients() and a right-hand side defined using
 * gfs_diffusion_rhs().
 */
void gfs_diffusion (GfsDomain * domain,
		    GfsMultilevelParams * par,
		    GfsVariable * v,
		    GfsVariable * rhs, 
		    GfsVariable * rhoc,
		    GfsVariable * metric)
{
  guint minlevel, maxlevel;
  GfsVariable * res;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);

  res = gfs_temporary_variable (domain);

  minlevel = domain->rootlevel;
  if (par->minlevel > minlevel)
    minlevel = par->minlevel;
  maxlevel = gfs_domain_depth (domain);
  gfs_diffusion_residual (domain, v, rhs, rhoc, metric, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_variable (domain, res, NULL, FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
  gdouble res_max_before = par->residual.infty;
  par->niter = 0;
  while (par->niter < par->nitermin ||
	 (par->residual.infty > par->tolerance && par->niter < par->nitermax)) {
    gfs_diffusion_cycle (domain, minlevel, maxlevel, par->nrelax, v, rhs, rhoc, metric, res);
    par->residual = gfs_domain_norm_variable (domain, res, NULL, 
					      FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
    if (par->residual.infty == res_max_before) /* convergence has stopped!! */
      break;
    if (par->residual.infty > res_max_before/1.1 && minlevel < maxlevel)
      minlevel++;
    res_max_before = par->residual.infty;
#if 0
    fprintf (stderr, "%d bias: %g first: %g second: %g infty: %g minlevel: %d\n",
	     par->niter, 
	     par->residual.bias, 
	     par->residual.first, 
	     par->residual.second, 
	     par->residual.infty,
	     minlevel);
#endif
    par->niter++;
  }

  gts_object_destroy (GTS_OBJECT (res));

  if (par->residual.infty > par->tolerance)
    g_warning ("gfs_diffusion(): max residual %g > %g", par->residual.infty, par->tolerance);
}

static GfsSourceDiffusion * source_diffusion (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (GFS_IS_SOURCE_DIFFUSION (o) && 
	  !GFS_IS_SOURCE_DIFFUSION_EXPLICIT (o) &&
	  !GFS_IS_SOURCE_VISCOSITY_EXPLICIT (o) &&
	  GFS_SOURCE_DIFFUSION (o)->phi == GFS_SOURCE_SCALAR (o)->v)
        return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static void add_pressure_gradient (FttCell * cell, GfsAdvectionParams * par)
{
  GFS_VALUE (cell, par->fv) -= GFS_VALUE (cell, par->g[par->v->component])*par->dt;
}

static void add_sinking (FttCellFace * face, GfsAdvectionParams * par)
{
  GFS_FACE_NORMAL_VELOCITY_LEFT (face) += gfs_function_face_value ( par->sink[face->d/2], face);
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) +=  gfs_function_face_value (par->sink[face->d/2], face);
}

/**
 * gfs_add_sinking_velocity:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 *
 * Adds the sinking velocity to the MAC velocity field of @domain.
 */
void gfs_add_sinking_velocity (GfsDomain * domain, GfsAdvectionParams * par)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  if (par->sink[0])
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) add_sinking, par);
}

static void remove_sinking (FttCellFace * face, GfsAdvectionParams * par)
{
  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= gfs_function_face_value (par->sink[face->d/2], face);
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= gfs_function_face_value (par->sink[face->d/2], face);
}

/**
 * gfs_remove_sinking_velocity:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 *
 * Removes the sinking velocity from the MAC velocity field of @domain.
 */
void gfs_remove_sinking_velocity (GfsDomain * domain, GfsAdvectionParams * par)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  if (par->sink[0])
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) remove_sinking, par);
}

static void variable_sources (GfsDomain * domain,
			      GfsAdvectionParams * par,
			      GfsVariable * sv,
			      GfsVariable ** gmac,
			      GfsVariable ** g)
{
  if (par->scheme == GFS_GODUNOV) {
    GfsVariable * v = par->v;

    par->u = gfs_domain_velocity (domain);
    par->g = gmac;
    par->fv = gfs_temporary_variable (domain);
    par->upwinding = GFS_FACE_UPWINDING;
    gfs_domain_face_traverse (domain, FTT_XYZ, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_reset, par->fv);
    gfs_add_sinking_velocity (domain, par);
    face_values_set ((FttCellTraverseFunc) gfs_cell_advected_face_values, par);
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) par->flux, par);
    gfs_remove_sinking_velocity (domain, par);
    par->v = sv;
    gfs_domain_traverse_merged (domain, par->update, par);
    par->v = v;
    par->u = par->g = NULL;
    gts_object_destroy (GTS_OBJECT (par->fv));
    par->fv = NULL;
  }
  if (g) {
    par->fv = sv;
    par->g = g;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_pressure_gradient, par);
    par->g = NULL;
    par->fv = NULL;
  }
  /* fixme: time should be set to t + dt/2 here for evaluation of
     source terms in the call below */
  par->fv = gfs_domain_variable_fluxes (domain, par->v, par->dt);
  if (par->fv) {
    GfsVariable * v = par->v;
    par->v = sv;
    /* fixme: for axi and moving should this be par->update? */
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, par);
    par->v = v;
    gts_object_destroy (GTS_OBJECT (par->fv));
    par->fv = NULL;
  }
  gfs_domain_variable_centered_sources (domain, par->v, sv, par->dt);
}

static void variable_diffusion (GfsDomain * domain,
				GfsSourceDiffusion * d,
				GfsAdvectionParams * par,
				GfsVariable * rhs,
				GfsFunction * alpha)
{
  GfsVariable * rhoc, * metric = NULL;

  rhoc = gfs_temporary_variable (domain);
  if (par->v->component < FTT_DIMENSION && domain->viscous_metric_implicit) {
    /* v is a component of a vector quantity and the domain has a metric */
    metric = gfs_temporary_variable (domain);
    metric->component = par->v->component;
  }

  gfs_domain_surface_bc (domain, par->v);
  gfs_diffusion_coefficients (domain, d, par->dt, rhoc, metric, alpha, d->D->par.beta);
  gfs_diffusion_rhs (domain, par->v, rhs, rhoc, metric, d->D->par.beta);
  /* fixme: time shoud be set to t + dt here in case boundary values are
     time-dependent in the call below */
  gfs_domain_surface_bc (domain, par->v);
  par->diffusion_solve (domain, &d->D->par, par->v, rhs, rhoc, metric);

  if (metric)
    gts_object_destroy (GTS_OBJECT (metric));
  gts_object_destroy (GTS_OBJECT (rhoc));
}

static void copy_v_rhs (FttCell * cell, GfsAdvectionParams * apar)
{
  GFS_VALUE (cell, apar->fv) = GFS_VALUE (cell, apar->v);
}

/**
 * gfs_centered_velocity_advection_diffusion:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @par: the advection parameters.
 * @gmac: the MAC pressure gradient.
 * @g: the pressure gradient.
 * @alpha: the inverse of density or %NULL.
 *
 * Advects the (centered) velocity field using the current
 * face-centered (MAC) velocity field and @par->flux to compute the
 * velocity flux through the faces of each cell.
 *
 * For each component of the velocity, before calling the @par->flux
 * function the face values are first defined (at time t + dt/2) and
 * can then be used within the @par->flux function.
 *
 * "Small" cut cells are treated using a cell-merging approach to
 * avoid any restrictive CFL stability condition.
 */
void gfs_centered_velocity_advection_diffusion (GfsDomain * domain,
						guint dimension,
						GfsAdvectionParams * par,
						GfsVariable ** gmac,
						GfsVariable ** g,
						GfsFunction * alpha)
{
  FttComponent c;
  GfsVariable ** v;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (gmac != NULL);

  gfs_domain_timer_start (domain, "centered_velocity_advection_diffusion");

  par->use_centered_velocity = FALSE;
  v = gfs_domain_velocity (domain);
  par->v = v[0];
  face_values_init ((FttCellTraverseFunc) gfs_cell_advected_face_values, par);
  for (c = 0; c < dimension; c++) {
    GfsSourceDiffusion * d = source_diffusion (v[c]);

    par->v = v[c];
    if (d) {
      GfsVariable * rhs;

      par->fv = rhs = gfs_temporary_variable (domain);
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) copy_v_rhs, par);
      variable_sources (domain, par, rhs, gmac, g);
      variable_diffusion (domain, d, par, rhs, alpha);
      gts_object_destroy (GTS_OBJECT (rhs));
    }
    else
      variable_sources (domain, par, par->v, gmac, g);
  }
  for (c = 0; c < dimension; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v[c]);
  face_values_free (par->v);

  gfs_domain_timer_stop (domain, "centered_velocity_advection_diffusion");
}

/**
 * gfs_tracer_advection_diffusion:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 * @alpha: the specific volume or %NULL.
 *
 * Advects the @v field of @par using the current face-centered (MAC)
 * velocity field.
 */
void gfs_tracer_advection_diffusion (GfsDomain * domain,
				     GfsAdvectionParams * par,
				     GfsFunction * alpha)
{
  GfsSourceDiffusion * d;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  gfs_domain_timer_start (domain, "tracer_advection_diffusion");

  if ((d = source_diffusion (par->v))) {
    GfsVariable * rhs;

    par->fv = rhs = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) copy_v_rhs, par);
    variable_sources (domain, par, rhs, NULL, NULL);
    variable_diffusion (domain, d, par, rhs, par->v->component < FTT_DIMENSION ? alpha : NULL);
    gts_object_destroy (GTS_OBJECT (rhs));
  }
  else {
    variable_sources (domain, par, par->v, NULL, NULL);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, par->v);
  }

  gfs_domain_timer_stop (domain, "tracer_advection_diffusion");
}

/**
 * Generic Surface boundary condition.
 * \beginobject{GfsSurfaceGenericBc}
 */

static void gfs_surface_generic_bc_destroy (GtsObject * o)
{
  if (GFS_SURFACE_GENERIC_BC (o)->v)
    GFS_SURFACE_GENERIC_BC (o)->v->surface_bc = NULL;

  (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->destroy) (o);
}

static void gfs_surface_generic_bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSurfaceGenericBc * bc = GFS_SURFACE_GENERIC_BC (*o);
  GtsObjectClass * klass;

  if (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a class name");
    return;
  }
  if (!(klass = gfs_object_class_from_name (fp->token->str))) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_surface_generic_bc_class ())) {
    gts_file_error (fp, "class `%s' is not a GfsSurfaceGenericClass", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  bc->v = gfs_variable_from_name (domain->variables, fp->token->str);
  if (!bc->v) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (bc->v->surface_bc) {
    gts_file_error (fp, "variable `%s' already has a surface boundary condition", 
		    fp->token->str);
    return;
  }
  bc->v->surface_bc = bc;
  gts_file_next_token (fp);
}

static void gfs_surface_generic_bc_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->write) (o, fp);
  fprintf (fp, "%s %s", o->klass->info.name, GFS_SURFACE_GENERIC_BC (o)->v->name);
}

static void gfs_surface_generic_bc_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_generic_bc_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_generic_bc_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_generic_bc_write;
}

GfsSurfaceGenericBcClass * gfs_surface_generic_bc_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_generic_bc_info = {
      "GfsSurfaceGenericBc",
      sizeof (GfsSurfaceGenericBc),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) gfs_surface_generic_bc_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_surface_generic_bc_info);
  }

  return klass;
}

/** \endobject{GfsSurfaceGenericBc} */

/**
 * Boundary conditions for diffusion and Poisson equations.
 * \beginobject{GfsSurfaceBc}
 */

static void gfs_surface_bc_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_SURFACE_BC (object)->type));
  gts_object_destroy (GTS_OBJECT (GFS_SURFACE_BC (object)->val));

  (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->destroy) (object);
}

static void gfs_surface_bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (!strcmp (fp->token->str, "Neumann")) {
    gfs_function_set_constant_value (bc->type, 0.);
    gts_file_next_token (fp);
  }
  else if (!strcmp (fp->token->str, "Dirichlet")) {
    gfs_function_set_constant_value (bc->type, 1.);
    gts_file_next_token (fp);
  }
  else {
    gfs_function_read (bc->type, gfs_object_simulation (bc), fp);
    if (fp->type == GTS_ERROR)
      return;
  }
  gfs_function_read (bc->val, gfs_object_simulation (bc), fp);
}

static void gfs_surface_bc_write (GtsObject * o, FILE * fp)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (o);
  gdouble val;

  if (GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->write) (o, fp);
  if ((val = gfs_function_get_constant_value (bc->type)) < G_MAXDOUBLE)
    fprintf (fp, " %s", val ? "Dirichlet" : "Neumann");
  else
    gfs_function_write (bc->type, fp);
  gfs_function_write (bc->val, fp);
}

static void gfs_surface_bc_bc (FttCell * cell, GfsSurfaceGenericBc * b)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (b);

  if (gfs_function_value (bc->type, cell) > 0.) {
    /* Dirichlet */
    cell->flags |= GFS_FLAG_DIRICHLET;
    gfs_function_set_units (bc->val, GFS_SURFACE_GENERIC_BC (bc)->v->units);
    GFS_STATE (cell)->solid->fv = gfs_function_value (bc->val, cell);
  }
  else {
    /* Neumann */
    cell->flags &= ~GFS_FLAG_DIRICHLET;
    gfs_function_set_units (bc->val, GFS_SURFACE_GENERIC_BC (bc)->v->units - 1.);
    FttVector n;
    gfs_solid_normal (cell, &n);
    GFS_STATE (cell)->solid->fv = gfs_function_value (bc->val, cell)*ftt_vector_norm (&n)
      *pow (ftt_cell_size (cell), FTT_DIMENSION - 1);;
  }
}

static void gfs_surface_bc_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_bc_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_bc_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_bc_destroy;
  klass->bc = gfs_surface_bc_bc;
}

static void gfs_surface_bc_init (GfsSurfaceBc * object)
{
  object->type = gfs_function_new (gfs_function_class (), 0.);
  object->val  = gfs_function_new (gfs_function_class (), 0.);
}

GfsSurfaceGenericBcClass * gfs_surface_bc_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_bc_info = {
      "GfsSurfaceBc",
      sizeof (GfsSurfaceBc),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) gfs_surface_bc_class_init,
      (GtsObjectInitFunc) gfs_surface_bc_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ()),
				  &gfs_surface_bc_info);
  }

  return klass;
}

/** \endobject{GfsSurfaceBc} */
