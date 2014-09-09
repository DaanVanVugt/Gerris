/* Gerris - The GNU Flow Solver
 * Copyright (C) 2012 S. Popinet, NIWA
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"

/* LayeredVariable: Header */

typedef struct {
  GfsVariable * v;
  GfsVariable ** vl;
} LayeredVariable;

/* GfsLayered: Header */

typedef struct _GfsLayered              GfsLayered;

struct _GfsLayered {
  /*< private >*/
  GfsSimulation parent;
  gdouble * ab;
  guint l;

  /*< public >*/
  LayeredVariable * u, * v, * lgmac[2], * lg[2];
  GSList * tracers, * variables;
  GfsVariable ** w, ** pr, ** un[FTT_NEIGHBORS], * gmac[FTT_DIMENSION], * g[FTT_DIMENSION];
  gdouble * dz, H;
  guint nl;        /**< number of layers */
  GfsFunction * b; /**< buoyancy */
};

#define GFS_LAYERED(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsLayered,	\
							 gfs_layered_class ())
#define GFS_IS_LAYERED(obj)         (gts_object_is_from_class (obj,	\
								   gfs_layered_class ()))

GfsSimulationClass * gfs_layered_class  (void);

/* LayeredVariable: Object */

static LayeredVariable * layered_variable_new (GfsVariable * v)
{
  LayeredVariable * lv = g_malloc (sizeof (LayeredVariable));  
  int l, nl = GFS_LAYERED (v->domain)->nl;
  lv->v = v;
  lv->vl = g_malloc (nl*sizeof (LayeredVariable));
  for (l = 0; l < nl; l++) {
    if (v->name) {
      gchar * name = g_strdup_printf ("%s%d", v->name, l);
      lv->vl[l] = gfs_variable_clone (v, name);
      g_free (name);
    }
    else
      lv->vl[l] = gfs_temporary_variable (v->domain);
  }
  lv->v = v;
  return lv;
}

static void layered_variable_average (FttCell * cell, LayeredVariable * lv)
{
  GfsLayered * layered = GFS_LAYERED (lv->v->domain);
  int l, nl = layered->nl;
  gdouble v = 0., * dz = layered->dz;
  for (l = 0; l < nl; l++)
    v += GFS_VALUE (cell, lv->vl[l])*dz[l];
  GFS_VALUE (cell, lv->v) = v;
}

static void layered_variable_swap (LayeredVariable * lv)
{
  GfsLayered * layered = GFS_LAYERED (lv->v->domain);
  gfs_variables_swap (lv->v, lv->vl[layered->l]);
}

static void layered_variable_destroy (LayeredVariable * lv)
{
  if (lv) {
    g_free (lv->vl);
    g_free (lv);
  }
}

/* GfsLayered: Object */

static void layered_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_layered_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }

  GfsLayered * layered = GFS_LAYERED (*o);
  GtsFileVariable var[] = {
    {GTS_UINT,   "nl", TRUE, &layered->nl},
    {GTS_DOUBLE, "H",  TRUE, &layered->H},
    {GTS_OBJ,    "b",  TRUE, &layered->b},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (layered);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    layered->un[d] = g_malloc (layered->nl*sizeof (GfsVariable *));
    int l;
    for (l = 0; l < layered->nl; l++)
      layered->un[d][l] = gfs_domain_add_variable (domain, NULL, NULL);
  }

  GfsVariable ** u = gfs_domain_velocity (domain);
  layered->u = layered_variable_new (u[0]);
  layered->v = layered_variable_new (u[1]);
  int l;
  for (l = 0; l < layered->nl; l++) {
    GfsVariable * u[2] = { layered->u->vl[l], layered->v->vl[l] };
    gfs_variable_set_vector (u, 2);
  }

  layered->gmac[0] = gfs_domain_add_variable (domain, NULL, NULL);
  layered->gmac[1] = gfs_domain_add_variable (domain, NULL, NULL);
  layered->g[0] = gfs_domain_add_variable (domain, NULL, NULL);
  layered->g[1] = gfs_domain_add_variable (domain, NULL, NULL);
  gfs_variable_set_vector (layered->gmac, FTT_DIMENSION);
  gfs_variable_set_vector (layered->g, FTT_DIMENSION);

  layered->lgmac[0] = layered_variable_new (layered->gmac[0]);
  layered->lgmac[1] = layered_variable_new (layered->gmac[1]);
  layered->lg[0] = layered_variable_new (layered->g[0]);
  layered->lg[1] = layered_variable_new (layered->g[1]);
  for (l = 0; l < layered->nl; l++) {
    GfsVariable * u[2] = { layered->lg[0]->vl[l], layered->lg[1]->vl[l] };
    gfs_variable_set_vector (u, 2);
    u[0] = layered->lgmac[0]->vl[l];
    u[1] = layered->lgmac[1]->vl[l];
    gfs_variable_set_vector (u, 2);
  }

  layered->w = g_malloc (layered->nl*sizeof (GfsVariable *));
  for (l = 0; l < layered->nl; l++) {
    gchar * name = g_strdup_printf ("W%d", l);
    layered->w[l] = gfs_domain_get_or_add_variable (domain, name, "z-component of the velocity");
    g_free (name);
  }
  
  layered->pr = g_malloc (layered->nl*sizeof (GfsVariable *));
  for (l = 0; l < layered->nl; l++) {
    gchar * name = g_strdup_printf ("Ph%d", l);
    layered->pr[l] = gfs_domain_get_or_add_variable (domain, name, "Hydrostatic potential");
    layered->pr[l]->units = 2.;
    g_free (name);
  }
  
  layered->dz = g_malloc (layered->nl*sizeof (gdouble));
  for (l = 0; l < layered->nl; l++)
    layered->dz[l] = 1./layered->nl;
  
  GSList * i = GFS_SIMULATION (layered)->events->items;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER (i->data))
      layered->tracers = g_slist_prepend (layered->tracers, layered_variable_new (i->data));
    else if (GFS_IS_SOURCE_CORIOLIS (i->data)) {
      GfsSourceCoriolis * s = i->data;
      if (s->u[0]) {
	layered->variables = g_slist_prepend (layered->variables, layered_variable_new (s->u[0]));
	layered->variables = g_slist_prepend (layered->variables, layered_variable_new (s->u[1]));
      }
    }
    else if (GFS_IS_EVENT_SUM (i->data)) {
      GfsEventSum * s = i->data;
      layered->variables = g_slist_prepend (layered->variables, layered_variable_new (s->sv));
    }
    i = i->next;
  }
}

static void layered_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_layered_class ())->parent_class->write) (o, fp);

  GfsLayered * layered = GFS_LAYERED (o);
  fprintf (fp, " { nl = %d H = %g b =", layered->nl, layered->H);
  gfs_function_write (layered->b, fp);
  fputs (" }", fp);
}

static void layered_destroy (GtsObject * object)
{
  GfsLayered * layered = GFS_LAYERED (object);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    g_free (layered->un[d]);
  
  layered_variable_destroy (layered->u);
  layered_variable_destroy (layered->v);
  layered_variable_destroy (layered->lgmac[0]);
  layered_variable_destroy (layered->lgmac[1]);
  layered_variable_destroy (layered->lg[0]);
  layered_variable_destroy (layered->lg[1]);
  g_free (layered->w);
  g_free (layered->pr);

  g_free (layered->dz);

  g_slist_foreach (layered->tracers, (GFunc) layered_variable_destroy, NULL);
  g_slist_free (layered->tracers);
  g_slist_foreach (layered->variables, (GFunc) layered_variable_destroy, NULL);
  g_slist_free (layered->variables);

  (* GTS_OBJECT_CLASS (gfs_layered_class ())->parent_class->destroy) (object);
}

static void swap_face_velocities (FttCell * cell, GfsLayered * layered)
{
  GfsStateVector * s = GFS_STATE (cell);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gdouble un = s->f[d].un;
    s->f[d].un = GFS_VALUE (cell, layered->un[d][layered->l]);
    GFS_VALUE (cell, layered->un[d][layered->l]) = un;
  }
}

static void swap_velocities (GfsLayered * layered)
{
  gfs_domain_traverse_leaves (GFS_DOMAIN (layered),
			      (FttCellTraverseFunc) swap_face_velocities, layered);
  layered_variable_swap (layered->u);
  layered_variable_swap (layered->v);
  g_slist_foreach (layered->variables, (GFunc) layered_variable_swap, NULL);
}

static void swap_gradients (GfsLayered * layered)
{
  layered_variable_swap (layered->lg[0]);
  layered_variable_swap (layered->lg[1]);
  layered_variable_swap (layered->lgmac[0]);
  layered_variable_swap (layered->lgmac[1]);
}

static void swap_layer (GfsLayered * layered)
{
  layered_variable_swap (layered->u);
  layered_variable_swap (layered->v);
  g_slist_foreach (layered->tracers, (GFunc) layered_variable_swap, NULL);
  g_slist_foreach (layered->variables, (GFunc) layered_variable_swap, NULL);
}

static void traverse_layers (GfsDomain * domain, FttCellTraverseFunc func, gpointer data)
{
  GfsLayered * layered = GFS_LAYERED (domain);
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_layer (layered);
    gfs_domain_traverse_leaves (domain, func, data);
    swap_layer (layered);
  }
}

static gdouble cell_z (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  GfsLayered * layered = GFS_LAYERED (sim);
  g_assert (layered->l < layered->nl);
  double z = layered->dz[layered->l]/2.;
  int i;
  for (i = 0; i < layered->l; i++)
    z += layered->dz[i];
  return z*layered->H;
}

static void layered_init (GfsLayered * object)
{
  GfsLayered * layered = GFS_LAYERED (object);
  layered->nl = 1;
  layered->H = 1.;
  GfsDomain * domain = GFS_DOMAIN (object);
  domain->traverse_layers = traverse_layers;
  GfsDerivedVariable * z = gfs_derived_variable_from_name (domain->derived_variables, "z");
  z->func = cell_z;
  layered->b = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (layered->b, 2.);
  gfs_object_simulation_set (layered->b, layered);
}

static void sum_face_velocities (FttCell * cell, GfsLayered * layered)
{
  GfsStateVector * s = GFS_STATE (cell);
  double dz = layered->dz[layered->l];
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    s->f[d].un += dz*GFS_VALUE (cell, layered->un[d][layered->l]);
}

static void compute_vertical_velocity (FttCell * cell, GfsLayered * layered)
{
  gdouble w = 0.;
  gdouble a = ftt_cell_size (cell)*gfs_domain_cell_fraction (GFS_DOMAIN (layered), cell)/layered->H;
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    FttCellFace face;
    gdouble div = 0.;
    face.cell = cell;
    for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++)
      div += (FTT_FACE_DIRECT (&face) ? 1. : -1.)*
	GFS_VALUE (cell, layered->un[face.d][layered->l])*
	gfs_domain_face_fraction (GFS_DOMAIN (layered), &face);
    w -= div*layered->dz[layered->l]/a;
    GFS_VALUE (cell, layered->w[layered->l]) = w;
  }
}

static void compute_hydrostatic_potential (FttCell * cell, GfsLayered * layered)
{
  double * ab = layered->ab, * dz = layered->dz, H = layered->H;
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    g_slist_foreach (layered->tracers, (GFunc) layered_variable_swap, NULL);
    ab[layered->l] = gfs_function_value (layered->b, cell);
    g_slist_foreach (layered->tracers, (GFunc) layered_variable_swap, NULL);
  }
  double * p = &GFS_VALUE (cell, layered->pr[0]), pr;
  int l, top = layered->nl - 1;
  pr = p[top] = 0.; //ab[top]*dz[top]*H/2.;
  for (l = top; l > 0; l--) {
    pr += (ab[l]*dz[l - 1] + ab[l - 1]*dz[l])*H/2.;
    p[l - 1] = pr;
  }
}

static void mac_projection (GfsLayered * layered, 
			    GfsMultilevelParams * par,
			    gdouble dt,
			    GfsVariable * p,
			    GfsVariable ** g)
{
  GfsSimulation * sim = GFS_SIMULATION (layered);
  GfsDomain * domain = GFS_DOMAIN (sim);

  /* project the summed MAC velocity field */
  gfs_mac_projection (domain, par, dt, p, sim->physical_params.alpha, g, NULL);
  
  /* apply barotropic pressure gradient to the MAC velocity field on each layer */
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_velocities (layered);
    gfs_correct_normal_velocities (domain, FTT_DIMENSION, p, NULL, dt);
    swap_velocities (layered);
  }

  /* compute vertical velocity */
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) compute_vertical_velocity, layered);
}

static void add_barotropic_gradient (FttCell * cell, GfsLayered * layered)
{
  int l;
  for (l = 0; l < layered->nl; l++) {
    GFS_VALUE (cell, layered->lg[0]->vl[l]) += GFS_VALUE (cell, layered->lg[0]->v);
    GFS_VALUE (cell, layered->lg[1]->vl[l]) += GFS_VALUE (cell, layered->lg[1]->v);
  }
}

static void add_barotropic_gmac (FttCell * cell, GfsLayered * layered)
{
  int l;
  for (l = 0; l < layered->nl; l++) {
    GFS_VALUE (cell, layered->lgmac[0]->vl[l]) += GFS_VALUE (cell, layered->lgmac[0]->v);
    GFS_VALUE (cell, layered->lgmac[1]->vl[l]) += GFS_VALUE (cell, layered->lgmac[1]->v);
  }
}

static void approximate_projection (GfsLayered * layered, GfsVariable * p)
{
  GfsSimulation * sim = GFS_SIMULATION (layered);
  GfsDomain * domain = GFS_DOMAIN (sim);

  /* compute face velocities for each layer */
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_velocities (layered);
    swap_gradients (layered);
    /* compute MAC velocities from centered velocities */
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, 
			      gfs_domain_velocity (domain));
    /* gradient of hydrostatic potential */
    gfs_reset_gradients (domain, FTT_DIMENSION, layered->g);
    gfs_correct_normal_velocities (domain, FTT_DIMENSION, layered->pr[layered->l],
				   layered->g, sim->advection_params.dt);
    gfs_scale_gradients (domain, FTT_DIMENSION, layered->g);
    gfs_correct_centered_velocities (domain, FTT_DIMENSION, layered->g, sim->advection_params.dt);
    swap_gradients (layered);
    swap_velocities (layered);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) sum_face_velocities, layered);
  }
  
  mac_projection (layered, &sim->approx_projection_params, sim->advection_params.dt, p, layered->g);

  /* apply barotropic pressure gradient to the centered velocity field on each layer */
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_velocities (layered);
    gfs_correct_centered_velocities (domain, FTT_DIMENSION, layered->g, sim->advection_params.dt);
    swap_velocities (layered);
  }

  /* add barotropic pressure gradient to hydrostatic potential gradient on each level */
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_barotropic_gradient, layered);

  /* store depth-averaged velocity in (U,V) */
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) layered_variable_average, layered->u);
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) layered_variable_average, layered->v);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, layered->u->v);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, layered->v->v);
}

typedef struct {
  GfsVariable * v;
  double dt;
  double * al, * ar;
} AdvectionParams;

static void cell_vertical_advection (FttCell * cell, AdvectionParams * p)
{
  double * al = p->al, * ar = p->ar, dt = p->dt;
  GfsLayered * layered = GFS_LAYERED (p->v->domain);
  double * a = &GFS_VALUE (cell, p->v);
  double * u = &GFS_VALUE (cell, layered->w[0]);
  double * dz = layered->dz, H = layered->H;

  int n = layered->nl, i;
  for (i = 0; i < n; i++) {
    double unorm = dt*((i > 0 ? u[i - 1] : 0.) + u[i])/(2.*dz[i]*H);
    if (fabs (unorm) > 1.)
      g_warning ("W CFL: %g", unorm);
    /* fixme: this gradient is correct only for dz[i] constant */
    double g = i == 0 ? a[i + 1] - a[i] : i == n - 1 ? a[i] - a[i-1] : (a[i + 1] - a[i - 1])/2.;
    al[i] = a[i] + MIN ((1. - unorm)/2., 0.5)*g;
    ar[i] = a[i] + MAX ((- 1. - unorm)/2., -0.5)*g;
  }
  for (i = 0; i < n - 1; i++) {
    double flux = (u[i] > 0. ? dt*u[i]*al[i] : 
		   u[i] < 0. ? dt*u[i]*ar[i + 1] :
		   dt*u[i]*(al[i] + ar[i + 1])/2.)/H;
    a[i] -= flux/dz[i];
    a[i + 1] += flux/dz[i + 1];
  }
}

static void vertical_advection (LayeredVariable * v, gdouble dt)
{
  GfsLayered * layered = GFS_LAYERED (v->v->domain);
  AdvectionParams p;
  p.v = v->vl[0]; p.dt = dt;
  p.al = g_malloc (layered->nl*sizeof (double));
  p.ar = g_malloc (layered->nl*sizeof (double));
  gfs_domain_traverse_leaves (v->v->domain, (FttCellTraverseFunc) cell_vertical_advection, &p);
  g_free (p.al);
  g_free (p.ar);
  int l;
  for (l = 0; l < layered->nl; l++)
    gfs_domain_bc (v->v->domain, FTT_TRAVERSE_LEAFS, -1, v->vl[l]);
}

static void advance_tracers (GfsLayered * layered, gdouble dt)
{
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_velocities (layered);
    g_slist_foreach (layered->tracers, (GFunc) layered_variable_swap, NULL);
    gfs_advance_tracers (GFS_SIMULATION (layered), dt);
    g_slist_foreach (layered->tracers, (GFunc) layered_variable_swap, NULL);
    swap_velocities (layered);
  }

  GfsDomain * domain = GFS_DOMAIN (layered);
  GSList * i = layered->tracers;
  while (i) {
    GfsVariable * v = ((LayeredVariable *) i->data)->v;
    if (GFS_VARIABLE_TRACER (v)->advection.scheme != GFS_NONE)
      vertical_advection (i->data, dt);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) layered_variable_average, i->data);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
    i = i->next;
  }

  layered->ab = g_malloc (layered->nl*sizeof (double));
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) compute_hydrostatic_potential, layered);
  g_free (layered->ab);
  int l;
  for (l = 0; l < layered->nl; l++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, layered->pr[l]);
}

static void layered_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsLayered * layered = GFS_LAYERED (sim);

  GfsVariable * p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  GfsVariable * pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  gfs_simulation_set_timestep (sim);
  if (sim->time.i == 0) {
    approximate_projection (layered, p);
    gfs_simulation_set_timestep (sim);
    advance_tracers (layered, sim->advection_params.dt/2.);
  }

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* compute face velocities at t + dt/2 for each layer */
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
    for (layered->l = 0; layered->l < layered->nl; layered->l++) {
      swap_velocities (layered);
      swap_gradients (layered);

      if (sim->advection_params.linear) {
	/* linearised advection */
	gfs_domain_face_traverse (domain, FTT_XYZ,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
	gfs_domain_face_traverse (domain, FTT_XYZ,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity,
				  sim->u0);
      }
      else
	gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);

      /* gradient of hydrostatic potential */
      gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
      gfs_reset_gradients (domain, FTT_DIMENSION, layered->gmac);
      gfs_correct_normal_velocities (domain, FTT_DIMENSION, layered->pr[layered->l], 
				     layered->gmac, sim->advection_params.dt/2.);
      gfs_scale_gradients (domain, FTT_DIMENSION, layered->gmac);
      swap_gradients (layered);
      swap_velocities (layered);
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) sum_face_velocities, layered);
    }

    gfs_variables_swap (p, pmac);
    mac_projection (layered, &sim->projection_params, sim->advection_params.dt/2., p, 
		    layered->gmac);
    /* add barotropic pressure gradient to hydrostatic potential gradient on each level */
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_barotropic_gmac, layered);
    int l;
    for (l = 0; l < layered->nl; l++) {
      /* we need to apply BC because
	 gfs_face_velocity_advection_flux() interpolates gmac on faces */
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, layered->lgmac[0]->vl[l]);
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, layered->lgmac[1]->vl[l]);
    }
    gfs_variables_swap (p, pmac);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    /* compute horizontal velocity at t + dt for each layer */
    for (layered->l = 0; layered->l < layered->nl; layered->l++) {
      swap_velocities (layered);
      swap_gradients (layered);
      gfs_centered_velocity_advection_diffusion (domain,
						 FTT_DIMENSION,
						 &sim->advection_params,
						 layered->gmac, 
						 sim->time.i > 0 ? layered->g : layered->gmac,
						 sim->physical_params.alpha);
      swap_gradients (layered);
      swap_velocities (layered);
    }

    if (sim->advection_params.scheme == GFS_GODUNOV) {
      vertical_advection (layered->u, sim->advection_params.dt);
      vertical_advection (layered->v, sim->advection_params.dt);
    }

    /* Coriolis */
    for (layered->l = 0; layered->l < layered->nl; layered->l++) {
      swap_velocities (layered);
      swap_gradients (layered);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, 
				       sim->time.i > 0 ? layered->g : layered->gmac, 
				       -sim->advection_params.dt);
      swap_gradients (layered);
      swap_velocities (layered);
    }

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    approximate_projection (layered, p);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    advance_tracers (layered, sim->advection_params.dt);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
}

typedef struct {
  GfsLayered * layered;
  gdouble cfl;
} CflData;

static void w_cfl (FttCell * cell, CflData * p)
{
  int l, nl = p->layered->nl;
  double * dz = p->layered->dz, H = p->layered->H;
  double * w = &GFS_VALUE (cell, p->layered->w[0]);
  for (l = 0; l < nl - 1; l++) 
    if (w[l] != 0.) {
      double wa = fabs(w[l])/H;
      double cfl = dz[l]/wa;
      if (cfl < p->cfl)
	p->cfl = cfl;
      cfl = dz[l + 1]/wa;
      if (cfl < p->cfl)
	p->cfl = cfl;
    }
}

static gdouble layered_cfl (GfsSimulation * sim)
{
  GfsLayered * layered = GFS_LAYERED (sim);
  CflData p = { layered, G_MAXDOUBLE };
  for (layered->l = 0; layered->l < layered->nl; layered->l++) {
    swap_velocities (layered);
    gdouble cfl = 
      (* GFS_SIMULATION_CLASS (GTS_OBJECT_CLASS (gfs_layered_class ())->parent_class)->cfl) (sim);
    if (cfl < p.cfl)
      p.cfl = cfl;
    swap_velocities (layered);
  }
  gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) w_cfl, &p);
  return p.cfl;
}

static void layered_class_init (GfsSimulationClass * klass) 
{
  GTS_OBJECT_CLASS (klass)->destroy = layered_destroy;
  GTS_OBJECT_CLASS (klass)->read =    layered_read;
  GTS_OBJECT_CLASS (klass)->write =   layered_write;
  klass->run =                        layered_run;
  klass->cfl =                        layered_cfl;
}

GfsSimulationClass * gfs_layered_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsLayered",
      sizeof (GfsLayered),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) layered_class_init,
      (GtsObjectInitFunc) layered_init
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "layered";
const gchar * g_module_check_init (void);
 
const gchar * g_module_check_init (void)
{
  gfs_layered_class ();
  return NULL;
} 
