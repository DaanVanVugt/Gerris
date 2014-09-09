/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010-2012 Jose M. López-Herrera Sánchez
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
#include "mpi_boundary.h"

/* GfsElectroHydro: Header */

typedef struct _GfsElectroHydro              GfsElectroHydro;

struct _GfsElectroHydro {
  /*< private >*/
  GfsSimulation parent;

  /*< public >*/
  GfsVariable * phi ;                             /* Electric potential */
  GfsVariable * E[FTT_DIMENSION] ;                /* Electric field; E=-Nabla Phi */
  GfsMultilevelParams electric_projection_params; /* Params for the electric potential */
  GfsFunction * perm ;                            /* electric permittivity */
  GfsFunction * charge ;                          /* charge density is user defined */
};

#define GFS_ELECTRO_HYDRO(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsElectroHydro,	\
							 gfs_electro_hydro_class ())
#define GFS_IS_ELECTRO_HYDRO(obj)         (gts_object_is_from_class (obj,	\
								   gfs_electro_hydro_class ()))

GfsSimulationClass * gfs_electro_hydro_class  (void);

/* GfsElectroHydro: Object */

static gboolean strmatch (const gchar * s, const gchar * s1)
{
  gboolean m = !strcmp (s, s1);

  if (!m) {
    gchar * s2 = g_strconcat ("Gfs", s, NULL);
    m = !strcmp (s2, s1);
    g_free (s2);
  }
  return m;
}

static GfsSourceDiffusion * source_implicit_ohmic (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    while (i) {
      GtsObject * o = i->data;
      if (GFS_IS_SOURCE_DIFFUSION (o) &&
	  !GFS_IS_SOURCE_DIFFUSION_EXPLICIT (o) &&
	  GFS_SOURCE_DIFFUSION (o)->phi == GFS_ELECTRO_HYDRO (v->domain)->phi)
	return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static GfsVariable * has_source_implicit_ohmic (GfsDomain * domain)
{
  GSList * i = domain->variables;
  while (i) {
    if (source_implicit_ohmic (i->data))
      return i->data;
    i = i->next;
  }
  return NULL;
}

static void gfs_electro_hydro_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_electro_hydro_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsElectroHydro * elec = GFS_ELECTRO_HYDRO (*o);
  GfsSimulation * sim = GFS_SIMULATION (elec);

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

    if (!strcmp (fp->token->str, "perm")) {
      gts_file_next_token (fp);
      if (fp->type != '=')
	gts_file_error (fp, "expecting `='");
      else {
	gts_file_next_token (fp);
	gfs_function_read (elec->perm, sim, fp);
      }
    }

    /* ------------ charge density defined by the function charge  ------------ */
    else if (!strcmp (fp->token->str, "charge")) {
      gts_file_next_token (fp);
      if (fp->type != '=')
	gts_file_error (fp, "expecting `='");
      else {
	gts_file_next_token (fp);
	gfs_function_read (elec->charge, sim, fp);
	GfsVariable * rhoe ;
	if (!gfs_function_get_variable (elec->charge) &&
	    (rhoe = has_source_implicit_ohmic (GFS_DOMAIN (sim))))
	  gts_file_error (fp, "for implicit charge diffusion, 'charge' must be equal to %s", 
			  rhoe->name);
      }
    }

    /* ------------ GfsElectricProjectionParams ------------ */
    else if (strmatch (fp->token->str, "GfsElectricProjectionParams")) {
      gts_file_next_token (fp);
      gfs_multilevel_params_read (&elec->electric_projection_params, fp);
    }

    else
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
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

static void gfs_electro_hydro_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_electro_hydro_class ())->parent_class->write) (o, fp);

  GfsElectroHydro * elect = GFS_ELECTRO_HYDRO (o);
  fputs (" {\n"
	 "  perm =", fp);
  gfs_function_write (elect->perm, fp);

  fputs (" \n"
	 "  charge =", fp);
  gfs_function_write (elect->charge, fp);

  fputs ("\n"
	 "  GfsElectricProjectionParams ", fp);
  gfs_multilevel_params_write (&elect->electric_projection_params, fp);
  fputs ("\n"
	 "}", fp);
}

static void gfs_electro_hydro_destroy (GtsObject * object)
{
  GfsElectroHydro * elec = GFS_ELECTRO_HYDRO (object);
  gts_object_destroy (GTS_OBJECT (elec->perm));
  gts_object_destroy (GTS_OBJECT (elec->charge));
  
  (* GTS_OBJECT_CLASS (gfs_electro_hydro_class ())->parent_class->destroy) (object);
}

/**
 * Electric field boundary condition.
 * \beginobject{GfsBcE}
 */

/* GfsBcE: Header */

#define GFS_IS_BC_E(obj)         (gts_object_is_from_class (obj,	\
					 gfs_bc_E_class ()))

GfsBcClass * gfs_bc_E_class  (void);

static void setting_E_from_phi (FttCellFace * f, GfsBc * b)
{
  if (b->v->component == f->d/2) {
    GfsVariable * phi = GFS_ELECTRO_HYDRO (gfs_object_simulation(b))->phi;
    GfsGradient g; 
    gfs_face_gradient (f, &g, phi->i, -1);
    double slope = (- g.b + g.a*GFS_VALUE (f->cell, phi))/ftt_cell_size (f->cell)
      *(FTT_FACE_DIRECT(f) ? 1 : -1);
    GFS_VALUE (f->cell, b->v) = - GFS_VALUE (f->neighbor, b->v) + 2.*slope;
  }
  else
    GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v);
}

static void face_setting_E_from_phi(FttCellFace *f, GfsBc * b)
{
  if (b->v->component == f->d/2) {
    GfsVariable * phi = GFS_ELECTRO_HYDRO (gfs_object_simulation(b))->phi;
    GfsGradient g; 
    gfs_face_gradient (f, &g, phi->i, -1);
    double slope = (- g.b + g.a*GFS_VALUE (f->cell, phi))/ftt_cell_size (f->cell)
      *(FTT_FACE_DIRECT(f) ? 1 : -1);
    GFS_STATE (f->cell)->f[f->d].v = 
      GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = slope;
  }
  else
    GFS_STATE (f->cell)->f[f->d].v = GFS_VALUE (f->neighbor, b->v);
}

static void gfs_bc_E_init (GfsBc * object)
{
  object->bc =                     (FttFaceTraverseFunc) setting_E_from_phi;
  object->homogeneous_bc =         (FttFaceTraverseFunc) setting_E_from_phi;
  object->face_bc =                (FttFaceTraverseFunc) face_setting_E_from_phi;
}

GfsBcClass * gfs_bc_E_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_E_info = {
      "GfsBcE",
      sizeof (GfsBc),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_bc_E_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_class ()),
				  &gfs_bc_E_info);
  }

  return klass;
}

/** \endobject{GfsBcE} */

static void box_set_efield_boundary (GfsBox * box, GfsElectroHydro * elec)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (!GFS_IS_BOUNDARY_MPI (box->neighbor[d]) &&
	GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++) {
	GfsBc * bc = gfs_bc_new (gfs_bc_E_class(), elec->E[c], FALSE);
	gfs_boundary_add_bc (b, bc);
      }
    }
}

static void gfs_electro_hydro_init (GfsElectroHydro * object)
{
  GfsDomain * domain = GFS_DOMAIN (object);
  static gchar name[][3] = {"Ex", "Ey", "Ez"};
  static gchar desc[][34] = {"x component of the electric field",
			     "y component of the electric field",
			     "z component of the electric field"};
  FttComponent c;  

  object->phi = gfs_domain_add_variable (domain, "Phi", "Electric potential");
  object->phi->centered = TRUE;

  for (c = 0; c < FTT_DIMENSION; c++) {
    object->E[c] = gfs_domain_add_variable (domain , name[c], desc[c]);
    object->E[c]->units = -1.;
  }
  gfs_variable_set_vector (object->E, FTT_DIMENSION);

  gfs_multilevel_params_init (&object->electric_projection_params);
  object->perm = gfs_function_new (gfs_function_class (), 1.);
  gfs_function_set_units (object->perm, -1.);

  object->charge = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (object->charge, -3.);
  gfs_object_simulation_set (object->charge, object);
}

static void gfs_electro_hydro_run (GfsSimulation * sim);

static void gfs_electro_hydro_class_init (GfsSimulationClass * klass) 
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_electro_hydro_destroy;
  GTS_OBJECT_CLASS (klass)->read =    gfs_electro_hydro_read;
  GTS_OBJECT_CLASS (klass)->write =   gfs_electro_hydro_write;
  klass->run =                        gfs_electro_hydro_run;
}

GfsSimulationClass * gfs_electro_hydro_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_electro_hydro_info = {
      "GfsElectroHydro",
      sizeof (GfsElectroHydro),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_electro_hydro_class_init,
      (GtsObjectInitFunc) gfs_electro_hydro_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()),
				  &gfs_electro_hydro_info);
  }

  return klass;
}

/* Setting div as - \int of rhoe on the cell volume */
static void rescale_div (FttCell * cell, GfsVariable * div)
{
  gdouble size = ftt_cell_size (cell);
  GFS_VALUE (cell, div) *= - size*size*gfs_domain_cell_fraction (div->domain, cell);
}

/* Calculates -gradient of @v and write it in vector @g  */
static void minus_gradient (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable ** g = data[1];
  gdouble size = ftt_cell_size (cell);
  FttVector gv;
  FttComponent c;

  gfs_cm_gradient (cell, v, &gv);
  for (c = 0; c < FTT_DIMENSION; c++) 
    GFS_VALUE (cell, g[c]) = - (&gv.x)[c]/size;
}

static void has_dirichlet (FttCell * cell, GfsVariable * p)
{
  if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    p->centered = FALSE;
}

typedef struct {
  GfsFunction * charge;
  GfsVariable * rhs;
} OhmicParams;

static void rhoe_update (FttCell * cell, gpointer * data)
{
  gdouble f, h, val;
  FttCellNeighbors neighbor;
  FttCellFace face;
  GfsVariable  * phi = data[0];
  GfsVariable * rhoe = data[1];
  
  if (GFS_IS_MIXED (cell)) {
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      f = gfs_cell_dirichlet_gradient_flux (cell, phi->i, -1, GFS_STATE (cell)->solid->fv);
    else
      f = GFS_STATE (cell)->solid->fv;
  }
  else
    f = 0.; /* Neumann condition by default */
  h = ftt_cell_size (cell);
  val = GFS_VALUE (cell, phi);
  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient g;

    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient (&face, &g, phi->i, -1);
    f += g.b - g.a*val;
  }
  GFS_VALUE (cell, rhoe) = -f/(h*h*gfs_domain_cell_fraction (rhoe->domain, cell));
}

static void charge_density_update (GfsDomain * domain, GfsVariable * phi, GfsVariable * rhoe)
{
  gpointer data[2];
  data[0] = phi;
  data[1] = rhoe;
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) rhoe_update, data);
}

static void set_dive (FttCell * cell, OhmicParams * p)
{
  GFS_VALUE (cell, p->rhs) = gfs_function_value (p->charge, cell);
}

static void poisson_electric (GfsElectroHydro * elec, gdouble dt)
{
  GfsMultilevelParams * par = &elec->electric_projection_params;
  GfsDomain * domain = GFS_DOMAIN (elec); 
  GfsVariable * diae, * dive, * res1e, * rhoe ;
  GfsVariable * phi = elec->phi; 
  GfsVariable ** e = elec->E;
  GfsSourceDiffusion * d;

  dive = gfs_temporary_variable (domain);

  OhmicParams p = { elec->charge, dive };
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) set_dive, &p);

  if ((rhoe = gfs_function_get_variable (elec->charge)))
    d = source_implicit_ohmic (rhoe);
  else
    d = NULL;

  if (d) {
    GfsVariable * rhoc = gfs_temporary_variable (domain);
    gfs_domain_surface_bc (domain, phi);
    gfs_diffusion_coefficients (domain, d, dt, rhoc, NULL, NULL, d->D->par.beta);
    gfs_diffusion_rhs (domain, phi, dive, rhoc, NULL, d->D->par.beta);
    gfs_poisson_coefficients (domain, elec->perm, FALSE, phi->centered, FALSE);
    gts_object_destroy (GTS_OBJECT (rhoc));
    par = &d->D->par;
  }
  else {
    gfs_domain_surface_bc (domain, phi);
    gfs_poisson_coefficients (domain, elec->perm, FALSE, phi->centered, TRUE);
  }
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) rescale_div, dive);

  res1e = gfs_temporary_variable (domain);
  diae = gfs_temporary_variable (domain);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, diae);
  par->poisson_solve (domain, par, phi, dive, res1e, diae, 1.);
  if (par->residual.infty > par->tolerance)
    g_warning ("poisson_electric: max residual %g > %g", par->residual.infty, par->tolerance);

  /* Set the electric field (-gradient of the potential) */
  gpointer data[2];
  data[0] = phi;
  data[1] = e;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) minus_gradient, data);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, e[c]);

  /* Compute the charge density from the electric potential */
  if (d) {
    gfs_poisson_coefficients (domain, elec->perm, FALSE, phi->centered, TRUE);
    charge_density_update (domain, phi, rhoe);
    /* fixme: update elec->rhoe bc ? */
  }

  gts_object_destroy (GTS_OBJECT (diae));
  gts_object_destroy (GTS_OBJECT (dive));
  gts_object_destroy (GTS_OBJECT (res1e));
}

static void gfs_electro_hydro_run (GfsSimulation * sim)
{
  GfsVariable * p, * pmac, * res = NULL, * g[FTT_DIMENSION], * gmac[FTT_DIMENSION];
  GfsVariable ** gc = sim->advection_params.gc ? g : NULL;
  GfsDomain * domain;
  GSList * i;
  GfsElectroHydro * elec;

  domain = GFS_DOMAIN (sim);
  elec = GFS_ELECTRO_HYDRO (sim) ;

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);

  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gmac[c] = gfs_temporary_variable (domain);
    if (sim->advection_params.gc) {
      g[c] = gfs_temporary_variable (domain);
    }
    else
      g[c] = gmac[c];
  }

  gfs_variable_set_vector (gmac, FTT_DIMENSION);
  gfs_variable_set_vector (g, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_domain_surface_bc (domain, elec->phi);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) has_dirichlet, elec->phi);
  if (!elec->phi->centered) {
    guint nf = gfs_check_solid_fractions (domain);
    if (nf > 0)
      g_warning ("the solid surface cuts %d boundary cells,\n"
		 "this may cause errors for the potential solution\n", nf);
  }
  gfs_simulation_init (sim);

  /* default BC for the electric field */
  gts_container_foreach (GTS_CONTAINER (domain),
  			 (GtsFunc) box_set_efield_boundary, elec);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  gfs_simulation_set_timestep (sim);
  if (sim->time.i == 0) {
    gfs_approximate_projection (domain,
				&sim->approx_projection_params,
				sim->advection_params.dt,
				p, sim->physical_params.alpha, res, g, NULL);
    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt/2.);
    poisson_electric (elec, sim->advection_params.dt/2.);
  }
  else if (sim->advection_params.gc)
    gfs_update_gradients (domain, p, sim->physical_params.alpha, g);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);
    
    gfs_variables_swap (p, pmac);
    gfs_mac_projection (domain,
    			&sim->projection_params, 
    			sim->advection_params.dt/2.,
			p, sim->physical_params.alpha, gmac, NULL);
    gfs_variables_swap (p, pmac);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_centered_velocity_advection_diffusion (domain,
					       FTT_DIMENSION,
					       &sim->advection_params,
					       gmac,
					       sim->time.i > 0 || !gc ? gc : gmac,
					       sim->physical_params.alpha);

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
    				sim->advection_params.dt, 
				p, sim->physical_params.alpha, res, g, NULL);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt);
    poisson_electric (elec, sim->advection_params.dt);

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

/* GfsSourceElectric: Header */

typedef struct _GfsSourceElectric         GfsSourceElectric;

struct _GfsSourceElectric {
  /*< private >*/
  GfsSourceVelocity parent;

  /*< public >*/
  GfsVariable * fe[FTT_DIMENSION];
};

#define GFS_SOURCE_ELECTRIC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceElectric,\
					         gfs_source_electric_class ())
#define GFS_IS_SOURCE_ELECTRIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_electric_class ()))

GfsSourceGenericClass * gfs_source_electric_class  (void);

/* GfsSourceElectric: Object */

static void gfs_source_electric_destroy (GtsObject * o)
{
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (GFS_SOURCE_ELECTRIC (o)->fe[c])
      gts_object_destroy (GTS_OBJECT (GFS_SOURCE_ELECTRIC (o)->fe[c]));
  
  (* GTS_OBJECT_CLASS (gfs_source_electric_class ())->parent_class->destroy) (o) ;
}

static void gfs_source_electric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_electric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  FttComponent c;
  for (c = 0 ; c < FTT_DIMENSION ; c++) {
    GfsVariable * v = GFS_SOURCE_VELOCITY (*o)->v[c];
    if (v->sources) {
      GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

      while (i) {
	if (i->data != *o && GFS_IS_SOURCE_ELECTRIC (i->data)) {
	  gts_file_error (fp, "variable '%s' cannot have multiple electric source terms", v->name);
	  return;
	}
	i = i->next;
      }
    }
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_SOURCE_ELECTRIC (*o)->fe[c] = gfs_temporary_variable (domain);
}

static void save_fe (FttCell * cell, GfsSourceElectric * s)
{
  GfsElectroHydro * elec = GFS_ELECTRO_HYDRO (gfs_object_simulation (s));

  GfsFunction * perm = elec->perm;
  GfsVariable ** e = elec->E;
  GfsVariable * phi = elec->phi;

  FttComponent c;
  gdouble h = ftt_cell_size (cell);

  FttCellFace f;
  FttCellNeighbors n;
  gdouble fe[FTT_DIMENSION];

  for (c = 0; c < FTT_DIMENSION; c++)
    fe[c] = 0.;

  f.cell = cell;
  ftt_cell_neighbors (cell, &n);

  gdouble radc = gfs_domain_cell_fraction (GFS_DOMAIN (elec), cell);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = n.c[f.d]; 
    gdouble permf = gfs_function_face_value (perm, &f);
    gdouble emod = 0.;
    GfsGradient g;
    /* fixme: should we use gfs_face_cm_weighted_gradient? */
    gfs_face_cm_gradient (&f, &g, phi->i, -1);
    gdouble en = (- g.b + g.a*GFS_VALUE (cell, phi))/h;
    gdouble sign = (FTT_FACE_DIRECT (&f) ? 1 : -1);

    gdouble radf = gfs_domain_face_fraction (GFS_DOMAIN (elec), &f);
    for (c = 0; c < FTT_DIMENSION; c++) {
      gdouble es = (c == f.d/2 ? sign*en : gfs_face_interpolated_value_generic (&f, e[c]));
      emod += es*es; 
      fe[c] += permf*es*en*radf;
    }
    fe[f.d/2] -= sign*emod*permf*radc/2.;
  }

  if (GFS_IS_MIXED (cell)) {
    if (((cell)->flags & GFS_FLAG_DIRICHLET) == 0)
      /* Neumann conditions for Phi */
      g_assert_not_implemented ();
    
    FttVector m = {1.,1.,1.};
    gfs_domain_solid_metric (GFS_DOMAIN (elec), cell, &m);
    gdouble permc = gfs_function_value (perm, cell);
    gdouble emod = 0., en = 0., a;
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector g, n;
    gfs_cell_dirichlet_gradient (cell, phi->i, -1, s->fv, &g);
    gfs_solid_normal (cell, &n);
    a = ftt_vector_norm (&n);
    for (c = 0; c < FTT_DIMENSION; c++) {
      (&n.x)[c] /= a;
      (&g.x)[c] /= h;
      emod += (&g.x)[c]*(&g.x)[c];
      en   += (&g.x)[c]*(&n.x)[c];
    }
    for (c = 0; c < FTT_DIMENSION; c++) 
      fe[c] += a*((&g.x)[c]*en*(&m.x)[c] - emod/2.*(&n.x)[c]*radc)*permc;
  }

  /* fixme: we need to rescale, not entirely clear why... */
  gdouble scale = pow (GFS_SIMULATION (elec)->physical_params.L, -5.);
  if (GFS_SIMULATION (elec)->physical_params.alpha)
    scale *= gfs_function_value (GFS_SIMULATION (elec)->physical_params.alpha, cell);
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, s->fe[c]) = scale*fe[c]/h/radc;
}

static gboolean gfs_source_electric_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_electric_class ())->parent_class)->event)
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_fe, event);
    return TRUE;
  }
  return FALSE;
}

static gdouble gfs_source_electric_centered_value (GfsSourceGeneric * s,
						   FttCell * cell,
						   GfsVariable * v)
{
  return GFS_VALUE (cell, GFS_SOURCE_ELECTRIC (s)->fe[v->component]);
}

static void gfs_source_electric_class_init (GfsSourceGenericClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_source_electric_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_electric_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_electric_destroy;
}

static void gfs_source_electric_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = gfs_source_electric_centered_value;
}

GfsSourceGenericClass * gfs_source_electric_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_electric_info = {
      "GfsSourceElectric",
      sizeof (GfsSourceElectric),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_electric_class_init,
      (GtsObjectInitFunc) gfs_source_electric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_electric_info);
  }

  return klass;
}

/* GfsElectroHydroAxi: Header */

GfsSimulationClass * gfs_electro_hydro_axi_class  (void);

/* GfsElectroHydroAxi: Object */

static void gfs_electro_hydro_axi_read (GtsObject ** o, GtsFile * fp)
{
  gfs_electro_hydro_read (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GFS_DOMAIN (*o)->refpos.y = 0.5;
}

static void gfs_electro_hydro_axi_class_init (GfsSimulationClass * klass) 
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_electro_hydro_destroy;
  GTS_OBJECT_CLASS (klass)->read =    gfs_electro_hydro_axi_read;
  GTS_OBJECT_CLASS (klass)->write =   gfs_electro_hydro_write;
  klass->run =                        gfs_electro_hydro_run;
}

GfsSimulationClass * gfs_electro_hydro_axi_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_electro_hydro_axi_info = {
      "GfsElectroHydroAxi",
      sizeof (GfsElectroHydro),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_electro_hydro_axi_class_init,
      (GtsObjectInitFunc) gfs_electro_hydro_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_axi_class ()),
				  &gfs_electro_hydro_axi_info);
  }

  return klass;
}

/* GfsOutputPotentialStats: Object */

static gboolean potential_stats_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsElectroHydro * elec = GFS_ELECTRO_HYDRO (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    if (elec->electric_projection_params.niter > 0) {
      fprintf (fp, "Electric potential    before     after       rate\n");
      gfs_multilevel_params_stats_write (&elec->electric_projection_params, fp);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_potential_stats_class_init (GfsEventClass * klass)
{
  klass->event = potential_stats_event;
}

GfsOutputClass * gfs_output_potential_stats_class (void);

GfsOutputClass * gfs_output_potential_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_potential_stats_info = {
      "GfsOutputPotentialStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_potential_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_potential_stats_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "electrohydro";
const gchar * g_module_check_init (void);
 
const gchar * g_module_check_init (void)
{
  gfs_electro_hydro_class ();
  gfs_electro_hydro_axi_class ();
  gfs_source_electric_class ();
  gfs_bc_E_class ();
  gfs_output_potential_stats_class ();
  return NULL;
} 
