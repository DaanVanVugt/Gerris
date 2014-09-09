/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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
 * \brief GfsWave spectral wave model.
 */

#include <stdlib.h>
#include "wave.h"
#include "adaptive.h"
#include "solid.h"
#include "init.h"

/**
 * Spectral wave model.
 * \beginobject{GfsWave}
 */

static double frequency (int ik)
{
  double gamma = GFS_WAVE_GAMMA;
  double f0 = GFS_WAVE_F0;
  return f0*pow(gamma, ik);
}
      
static double theta (guint ith, guint ntheta)
{
  return 2.*M_PI*ith/ntheta;
}

static void group_velocity (int ik, int ith, FttVector * u, guint ntheta, gdouble g)
{
  double cg = g/(4.*M_PI*frequency (ik));
  u->x = cg*cos (theta (ith, ntheta));
  u->y = cg*sin (theta (ith, ntheta));
  u->z = 0.;
}

static gdouble cell_E (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  GfsWave * wave = GFS_WAVE (domain);
  GfsVariable *** F = wave->F;
  guint ik, ith;
  gdouble E = 0., sigma = 2.*M_PI*GFS_WAVE_F0, sgamma = (GFS_WAVE_GAMMA - 1./GFS_WAVE_GAMMA)/2.;
  for (ik = 0; ik < wave->nk; ik++) {
    gdouble df = sigma*sgamma;
    gdouble dE = 0.;
    for (ith = 0; ith < wave->ntheta; ith++)
      dE += GFS_VALUE (cell, F[ik][ith]);
    E += dE*df;
    sigma *= GFS_WAVE_GAMMA;
  }
  return E*2.*M_PI/wave->ntheta;
}

static void set_group_velocity (const FttCellFace * face, FttVector * u)
{
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) = 
    GFS_FACE_NORMAL_VELOCITY_LEFT (face) = (&u->x)[face->d/2];
}

typedef struct {
  GfsAdvectionParams * p;
  GfsVariable * div, * fv;
} SolidFluxParams;

static void solid_flux (FttCell * cell, SolidFluxParams * par)
{
  gfs_normal_divergence (cell, par->div);
  if (GFS_VALUE (cell, par->div) < 0.) {
    gdouble h = ftt_cell_size (cell);
    GFS_VALUE (cell, par->fv) = GFS_VALUE (cell, par->div)*par->p->dt*
      GFS_VALUE (cell, par->p->v)/(h*h);
  }
  else
    GFS_VALUE (cell, par->fv) = 0.;
}

typedef struct {
  GfsVariable * F, * Fn, * dF;
  gdouble D[2][2];
} GSEData;

static void compute_gradient (FttCell * cell, GSEData * p)
{
  GFS_VALUE (cell, p->dF) = gfs_center_regular_gradient (cell, p->dF->component, p->Fn);
}

static void diffusion (FttCell * cell, GSEData * p)
{
  gdouble h2 = ftt_cell_size (cell);
  h2 *= h2;
  /* off-diagonal */
  FttComponent j;
  for (j = 0; j < 2; j++)
    if (j != p->dF->component)
      GFS_VALUE (cell, p->F) += 
	p->D[j][p->dF->component]*gfs_center_regular_gradient (cell, j, p->dF)/h2;
  /* diagonal */
  GFS_VALUE (cell, p->F) += 
    p->D[p->dF->component][p->dF->component]*
    gfs_center_regular_2nd_derivative (cell, p->dF->component, p->Fn)/h2;
}

static void copy_F (FttCell * cell, GSEData * p)
{
  GFS_VALUE (cell, p->Fn) = GFS_VALUE (cell, p->F);
}

static void gse_alleviation_diffusion (GfsDomain * domain, GfsVariable * F,
				       FttVector * cg, gdouble dt)
{
  gfs_domain_timer_start (domain, "gse_alleviation");

  gdouble ncg = sqrt (cg->x*cg->x + cg->y*cg->y);
  gdouble dcg = (GFS_WAVE_GAMMA - 1./GFS_WAVE_GAMMA)*ncg/2.;
  gdouble dtheta = 2.*M_PI/GFS_WAVE (domain)->ntheta;
#if 0
  gdouble Ts = 4.*GFS_WAVE (domain)->alpha_s*GFS_WAVE (domain)->alpha_s*dt;
  gdouble dtDss = dt*dcg*dcg*Ts/12.;
  gdouble dtDnn = dt*(ncg*dtheta)*(ncg*dtheta)*Ts/12.;
#else
  gdouble alpha = GFS_WAVE (domain)->alpha_s*dcg*dt;
  gdouble beta = GFS_WAVE (domain)->alpha_s*ncg*dtheta*dt;
  gdouble dtDss = alpha*alpha/3.;
  gdouble dtDnn = beta*beta/3.;
#endif
  GSEData p;
  gdouble cost = cg->x/ncg, sint = cg->y/ncg;
  p.D[0][0] = dtDss*cost*cost + dtDnn*sint*sint;
  p.D[1][1] = dtDss*sint*sint + dtDnn*cost*cost;
  p.D[0][1] = p.D[1][0] = (dtDss - dtDnn)*cost*sint;
  p.F = F;
  p.Fn = gfs_temporary_variable (domain);
  p.dF = gfs_temporary_variable (domain);
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) copy_F, &p);
  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) p.Fn->fine_coarse, p.Fn);
  for (p.dF->component = 0; p.dF->component < 2; p.dF->component++) {
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) compute_gradient, &p);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, p.dF);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) diffusion, &p);
  }
  gts_object_destroy (GTS_OBJECT (p.Fn));
  gts_object_destroy (GTS_OBJECT (p.dF));
  gfs_domain_timer_stop (domain, "gse_alleviation");
}

static void redo_some_events (GfsEvent * event, GfsSimulation * sim)
{
  if (GFS_IS_ADAPT (event) || GFS_IS_INIT (event))
    gfs_event_redo (event, sim);
}

static void wave_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsWave * wave = GFS_WAVE (sim);

  SolidFluxParams par;
  par.div = gfs_variable_from_name (domain->variables, "P");
  g_assert (par.div);
  par.p = &sim->advection_params;
  par.fv = gfs_temporary_variable (domain);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* get global timestep */
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
    gfs_simulation_set_timestep (sim);
    gdouble dt = sim->advection_params.dt;
    gdouble g = sim->physical_params.g/sim->physical_params.L;
    gdouble tnext = sim->tnext;
    
    /* spatial advection */
    guint ik, ith;
    for (ik = 0; ik < wave->nk; ik++) {
      FttVector cg;
      group_velocity (ik, 0, &cg, wave->ntheta, g);
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) set_group_velocity, &cg);
      if (wave->alpha_s > 0.) {
	/* stability criterion for GSE diffusion */
	gdouble cfl = sim->advection_params.cfl;
	sim->advection_params.cfl = MIN (cfl, 2./(4.*wave->alpha_s*M_PI/wave->ntheta));
	/* fixme: this should be:
	   sim->advection_params.cfl = MIN (cfl, sqrt(3.)/(wave->alpha_s*2.*M_PI/wave->ntheta));
	*/
	gfs_simulation_set_timestep (sim);
	sim->advection_params.cfl = cfl;
      }
      else
	gfs_simulation_set_timestep (sim);
      /* subcycling */
      guint n = rint (dt/sim->advection_params.dt);
      g_assert (fabs (sim->time.t + sim->advection_params.dt*n - tnext) < 1e-12);
      while (n--) {
	for (ith = 0; ith < wave->ntheta; ith++) {
	  FttVector cg;
	  group_velocity (ik, ith, &cg, wave->ntheta, g);
	  gfs_domain_face_traverse (domain, FTT_XYZ,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) set_group_velocity, &cg);
	  GfsVariable * t = GFS_WAVE (sim)->F[ik][ith];
	  sim->advection_params.v = t;
	  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) solid_flux, &par);
	  gfs_tracer_advection_diffusion (domain, &sim->advection_params, NULL);
	  sim->advection_params.fv = par.fv;
	  gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, 
	  			      &sim->advection_params);
	  if (wave->alpha_s > 0.)
	    gse_alleviation_diffusion (domain, t, &cg, sim->advection_params.dt);
	  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, t);
	  gfs_domain_cell_traverse (domain,
				    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				    (FttCellTraverseFunc) t->fine_coarse, t);
	}
	gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) redo_some_events, sim);
	gfs_simulation_adapt (sim);
      }
    }

    sim->advection_params.dt = dt;

    /* source terms */
    if (wave->source)
      (* wave->source) (wave);

    sim->time.t = sim->tnext = tnext;
    sim->time.i++;

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (par.fv));
}

static void wave_destroy (GtsObject * object)
{
  if (GFS_WAVE (object)->F)
    gfs_matrix_free (GFS_WAVE (object)->F);
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->destroy) (object);
}

static void wave_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsWave * wave = GFS_WAVE (*o);
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_UINT,   "nk",      TRUE, &wave->nk},
      {GTS_UINT,   "ntheta",  TRUE, &wave->ntheta},
      {GTS_DOUBLE, "alpha_s", TRUE, &wave->alpha_s},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
  }

  GfsDomain * domain = GFS_DOMAIN (wave);
  guint ik, ith;
  wave->F = gfs_matrix_new (wave->nk, wave->ntheta, sizeof (GfsVariable *));
  for (ik = 0; ik < wave->nk; ik++)
    for (ith = 0; ith < wave->ntheta; ith++) {
      gchar * name = g_strdup_printf ("F%d_%d", ik, ith);
      gchar * description = g_strdup_printf ("Action density for f = %g Hz and theta = %g degrees",
					     frequency (ik), theta (ith, wave->ntheta)*180./M_PI);
      wave->F[ik][ith] = gfs_domain_get_or_add_variable (domain, name, description);
      g_assert (wave->F[ik][ith]);
      g_free (name);
      g_free (description);
    }
}

static void wave_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->write) (o, fp);

  GfsWave * wave = GFS_WAVE (o);
  fprintf (fp, " {\n"
	   "  nk = %d\n"
	   "  ntheta = %d\n"
	   "  alpha_s = %g\n"
	   "}",
	   wave->nk, wave->ntheta, wave->alpha_s);
}

static void gfs_wave_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = wave_destroy;
  GTS_OBJECT_CLASS (klass)->read = wave_read;
  GTS_OBJECT_CLASS (klass)->write = wave_write;
  klass->run = wave_run;
}

static gdouble cell_hs (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble E = cell_E (cell, face, domain);
  return E > 0. ? 4.*sqrt (E) : 0.;
}

static gdouble cell_frequency (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return frequency (GFS_WAVE (domain)->ik);
}

static gdouble cell_direction (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return theta (GFS_WAVE (domain)->ith, GFS_WAVE (domain)->ntheta);
}

static void wave_init (GfsWave * wave)
{
  wave->nk = 25;
  wave->ntheta = 24;
  wave->alpha_s = 0.;
  /* default for g is acceleration of gravity on Earth with kilometres as
     spatial units, hours as time units and Hz as frequency units */
  GFS_SIMULATION (wave)->physical_params.g = 9.81/1000.*3600.;

  GfsAdvectionParams * par = &GFS_SIMULATION (wave)->advection_params;
  par->gradient = gfs_center_van_leer_gradient;
  par->flux = gfs_face_advection_flux;
  par->use_centered_velocity = FALSE;  

  static GfsDerivedVariableInfo derived_variable[] = {
    { "Hs", "Significant wave height", cell_hs },
    { "Energy", "Wave energy", cell_E },
    { "Frequency", "Wave frequency", cell_frequency },
    { "Direction", "Wave direction (angle)", cell_direction },
    { NULL, NULL, NULL}
  };
  GfsDerivedVariableInfo * v = derived_variable;
  while (v->name) {
    g_assert (gfs_domain_add_derived_variable (GFS_DOMAIN (wave), *v));
    v++;
  }
}

GfsSimulationClass * gfs_wave_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_wave_info = {
      "GfsWave",
      sizeof (GfsWave),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_wave_class_init,
      (GtsObjectInitFunc) wave_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_wave_info);
  }

  return klass;
}

/** \endobject{GfsWave} */

/**
 * Initial wave spectrum for wave model.
 * \beginobject{GfsInitWave}
 */

static void gfs_init_wave_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!GFS_IS_WAVE (domain)) {
    gts_file_error (fp, "GfsInitWave can only be used within a GfsWave simulation");
    return;
  }
  
  gfs_function_read (GFS_INIT_WAVE (*o)->d, domain, fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_read (GFS_INIT_WAVE (*o)->hs, domain, fp);
}

static void gfs_init_wave_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->write) (o, fp);

  gfs_function_write (GFS_INIT_WAVE (o)->d, fp);
  gfs_function_write (GFS_INIT_WAVE (o)->hs, fp);
}

static void gfs_init_wave_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_INIT_WAVE (object)->d));
  gts_object_destroy (GTS_OBJECT (GFS_INIT_WAVE (object)->hs));

  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->destroy) (object);
}

static void init_energy (FttCell * cell, GfsInitWave * event)
{
  GfsWave * wave = GFS_WAVE (gfs_object_simulation (event));
  for (wave->ik = 0; wave->ik < wave->nk; wave->ik++)
    for (wave->ith = 0; wave->ith < wave->ntheta; wave->ith++)
      GFS_VALUE (cell, wave->F[wave->ik][wave->ith]) = gfs_function_value (event->d, cell);
}

static void scale_energy (FttCell * cell, GfsInitWave * event)
{
  GfsWave * wave = GFS_WAVE (gfs_object_simulation (event));
  gdouble E = cell_E (cell, NULL, GFS_DOMAIN (wave));
  if (E > 0.) {
    gdouble Hs = gfs_function_value (event->hs, cell);
    gdouble scaling = Hs*Hs/(16.*E);
    guint ik, ith;
    for (ik = 0; ik < wave->nk; ik++)
      for (ith = 0; ith < wave->ntheta; ith++)
	GFS_VALUE (cell, wave->F[ik][ith]) *= scaling;
  }
}

static gboolean gfs_init_wave_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_energy, event);
    gfs_restore_fpe_for_function (GFS_INIT_WAVE (event)->d);
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) scale_energy, event);
    gfs_restore_fpe_for_function (GFS_INIT_WAVE (event)->hs);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_wave_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_wave_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_wave_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_wave_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_wave_destroy;
}

static void gfs_init_wave_init (GfsInitWave * object)
{
  object->d = gfs_function_new (gfs_function_class (), 0.);
  object->hs = gfs_function_new (gfs_function_class (), 0.);
}

GfsGenericInitClass * gfs_init_wave_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_wave_info = {
      "GfsInitWave",
      sizeof (GfsInitWave),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_wave_class_init,
      (GtsObjectInitFunc) gfs_init_wave_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_wave_info);
  }

  return klass;
}

/** \endobject{GfsInitWave} */
