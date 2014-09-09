/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2011 National Institute of Water and Atmospheric Research
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
 * \brief GfsSimulation, GfsPoisson, GfsAdvection models.
 */

#include <stdlib.h>
#include <math.h>
#include <gmodule.h>
#include "config.h"
#include "simulation.h"
#include "output.h"
#include "refine.h"
#include "solid.h"
#include "adaptive.h"
#include "source.h"
#include "vof.h"
#include "tension.h"
#include "map.h"
#include "river.h"
#include "version.h"

/**
 * The incompressible Euler solver.
 * \beginobject{GfsSimulation}
 */

static void module_close (GModule *module)
{
  if (!g_module_close (module))
    g_warning ("%s: %s", g_module_name (module), g_module_error ());
}

static void simulation_destroy (GtsObject * object)
{
  GfsSimulation * sim = GFS_SIMULATION (object);

  gts_container_foreach (GTS_CONTAINER (sim->refines), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (sim->refines));

  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (sim->events));

  gts_container_foreach (GTS_CONTAINER (sim->maps), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (sim->maps));

  gts_object_destroy (GTS_OBJECT (sim->adapts));
  gts_object_destroy (GTS_OBJECT (sim->solids));

  g_slist_foreach (sim->modules, (GFunc) module_close, NULL);
  g_slist_free (sim->modules);
  g_slist_foreach (sim->globals, (GFunc) gts_object_destroy, NULL);
  g_slist_free (sim->globals);
  g_slist_foreach (sim->preloaded_modules, (GFunc) module_close, NULL);
  g_slist_free (sim->preloaded_modules);

  (* GTS_OBJECT_CLASS (gfs_simulation_class ())->parent_class->destroy) (object);
}

static void simulation_write (GtsObject * object, FILE * fp)
{
  GfsSimulation * sim = GFS_SIMULATION (object);
  GfsDomain * domain = GFS_DOMAIN (sim);
  GSList * i;
  GfsVariable * v;

  (* GTS_OBJECT_CLASS (gfs_simulation_class ())->parent_class->write)
    (object, fp);

  fputs (" {\n", fp);

  i = sim->modules;
  while (i) {
    void (* module_write) (FILE *);
    const gchar * name = NULL;
    fprintf (fp, "  GModule %s", 
	     g_module_symbol (i->data, "gfs_module_name", (gpointer) &name) ? name : 
	     g_module_name (i->data));
    if (g_module_symbol (i->data, "gfs_module_write", (gpointer) &module_write))
      (* module_write) (fp);
    fputc ('\n', fp);
    i = i->next;
  }

  i = sim->globals;
  while (i) {
    fputs ("  ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    i = i->next;
  }

  fputs ("  GfsTime ", fp);
  gfs_time_write (&sim->time, fp);
  fputc ('\n', fp);

  i = sim->refines->items;
  while (i) {
    GtsObject * object = i->data;
    if (GFS_IS_LAYERS (object) || domain->max_depth_write < -1) {
      fputs ("  ", fp);
      g_assert (object->klass->write);
      (* object->klass->write) (object, fp);
      fputc ('\n', fp);
    }
    i = i->next;
  }

  i = sim->events->items;
  while (i) {
    GtsObject * object = i->data;
    GfsEvent * event = i->data;
    
    if (event->t < event->end && event->i < event->iend) {
      fputs ("  ", fp);
      g_assert (object->klass->write);
      (* object->klass->write) (object, fp);
      fputc ('\n', fp);
    }
    i = i->next;
  }

  i = domain->variables;
  while (i) {
    v = i->data;
    if (v->surface_bc) {
      fputs ("  ", fp);
      (* GTS_OBJECT (v->surface_bc)->klass->write) (GTS_OBJECT (v->surface_bc), fp);
      fputc ('\n', fp);
    }
    i = i->next;
  }

  fputs ("  GfsPhysicalParams ", fp);
  gfs_physical_params_write (&sim->physical_params, fp);
  fputs ("\n  GfsAdvectionParams ", fp);
  gfs_advection_params_write (&sim->advection_params, fp);
  fputs ("\n  GfsApproxProjectionParams ", fp);
  gfs_multilevel_params_write (&sim->approx_projection_params, fp);
  fputs ("\n  GfsProjectionParams ", fp);
  gfs_multilevel_params_write (&sim->projection_params, fp);
  fputc ('\n', fp);

  i = sim->maps->items;
  while (i) {
    GtsObject * object = i->data;
    g_assert (object->klass->write);
    (* object->klass->write) (object, fp);
    fputc ('\n', fp);
    i = i->next;
  }

  fputc ('}', fp);
}

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

static GModule * load_module (GtsFile * fp, GfsSimulation * sim)
{
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (module name)");
    return NULL;
  }
  if (!g_module_supported ()) {
    g_warning ("modules are not supported on this system");
    gts_file_next_token (fp);
    return NULL;
  }
  else {
    GModule * module;

    module = g_module_open (fp->token->str, 0);
    if (module == NULL) {
      gchar * name = g_strconcat (fp->token->str, 
#if FTT_2D
				  "2D"
#else
				  "3D"
#endif
				  , NULL);
#ifdef MODULES_SUFFIX
      gchar * path;
      if (name[0] != G_DIR_SEPARATOR)
	path = g_strconcat (GFS_MODULES_DIR, G_DIR_SEPARATOR_S, "lib", name, MODULES_SUFFIX, NULL);
      else
	path = g_strdup (name);
#else
      gchar * path = g_module_build_path (GFS_MODULES_DIR, name);
#endif
      g_free (name);
      module = g_module_open (path, 0);
      g_free (path);
    }
    if (module == NULL) {
      gts_file_error (fp, "cannot load module: %s", g_module_error ());
      return NULL;
    }
    g_module_make_resident (module);
    gts_file_next_token (fp);

    void (* module_read) (GtsFile *, GfsSimulation * sim);
    if (g_module_symbol (module, "gfs_module_read", (gpointer) &module_read)) {
      (* module_read) (fp, sim);
      if (fp->type == GTS_ERROR)
	return NULL;
    }
    return module;
  }
}

static void simulation_read (GtsObject ** object, GtsFile * fp)
{
  GfsSimulation * sim = GFS_SIMULATION (*object);

  (* GTS_OBJECT_CLASS (gfs_simulation_class ())->parent_class->read) (object, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);

  GfsDomain * domain = GFS_DOMAIN (sim);
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }

    /* ------ fixme: obsolete, just here for backward compatibility ------*/
    if (!strcmp (fp->token->str, "GfsDeferredCompilation"))
      gts_file_next_token (fp);

    /* ------------ GModule ------------ */
    else if (!strcmp (fp->token->str, "GModule")) {
      GModule * module = load_module (fp, sim);
      if (module == NULL)
	return;
      sim->modules = g_slist_prepend (sim->modules, module);
    }

    /* ------------ GfsTime ------------ */
    else if (strmatch (fp->token->str, "GfsTime")) {
      gts_file_next_token (fp);
      gfs_time_read (&sim->time, fp);
      if (fp->type == GTS_ERROR)
	return;
    }

    /* ------------ GfsPhysicalParams ------------ */
    else if (strmatch (fp->token->str, "GfsPhysicalParams")) {
      gts_file_next_token (fp);
      gfs_physical_params_read (&sim->physical_params, domain, fp);
      if (fp->type == GTS_ERROR)
	return;
    }

    /* ------------ GfsProjectionParams ------------ */
    else if (strmatch (fp->token->str, "GfsProjectionParams")) {
      gts_file_next_token (fp);
      gfs_multilevel_params_read (&sim->projection_params, fp);
      if (fp->type == GTS_ERROR)
	return;
    }

    /* ------------ GfsApproxProjectionParams ------------ */
    else if (strmatch (fp->token->str, "GfsApproxProjectionParams")) {
      gts_file_next_token (fp);
      gfs_multilevel_params_read (&sim->approx_projection_params, fp);
      if (fp->type == GTS_ERROR)
	return;
    }

    /* ------------ GfsAdvectionParams ------------ */
    else if (strmatch (fp->token->str, "GfsAdvectionParams")) {
      gts_file_next_token (fp);
      gfs_advection_params_read (&sim->advection_params, fp);
      if (fp->type == GTS_ERROR)
	return;
      if (sim->advection_params.linear) {
	sim->u0[0] = gfs_domain_get_or_add_variable (domain, 
						     "U0", "x-component of the base velocity");
	sim->u0[0]->units = 1.;
	sim->u0[1] = gfs_domain_get_or_add_variable (domain, 
						     "V0", "y-component of the base velocity");
	sim->u0[1]->units = 1.;
#if (!FTT_2D)
	sim->u0[2] = gfs_domain_get_or_add_variable (domain, 
						     "W0", "z-component of the base velocity");
	sim->u0[2]->units = 1.;
#endif /* FTT_3D */
	gfs_variable_set_vector (sim->u0, FTT_DIMENSION);
      }
    }

    /* ------------ GtsObject ------------ */
    else {
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      GtsObject * object;

      if (klass == NULL) {
	gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
	return;
      }
      if (gts_object_class_is_from_class (klass, gfs_box_class ())) {
	gts_file_error (fp, "parse error (unclosed statement?)");
	return;
      }

      object = gts_object_new (klass);
      gfs_object_simulation_set (object, sim);

      g_assert (klass->read);
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	return;
      }

      if (GFS_IS_GLOBAL (object))
	sim->globals = g_slist_append (sim->globals, object);
      else if (GFS_IS_REFINE (object))
	gts_container_add (GTS_CONTAINER (sim->refines), GTS_CONTAINEE (object));
      else if (GFS_IS_ADAPT (object)) {
	gts_container_add (GTS_CONTAINER (sim->adapts), GTS_CONTAINEE (object));
	gts_container_add (GTS_CONTAINER (sim->events), GTS_CONTAINEE (object));
      }
      else if (GFS_IS_SOLID (object)) {
	gts_container_add (GTS_CONTAINER (sim->solids), GTS_CONTAINEE (object));
	gts_container_add (GTS_CONTAINER (sim->events), GTS_CONTAINEE (object));
      }
      else if (GFS_IS_EVENT (object))
	gts_container_add (GTS_CONTAINER (sim->events), GTS_CONTAINEE (object));
      else if (GFS_IS_MAP (object))
	gts_container_add (GTS_CONTAINER (sim->maps), GTS_CONTAINEE (object));
      else if (GFS_IS_SURFACE_GENERIC_BC (object))
	;
      else
	g_assert_not_reached ();
    }
  }
  
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  sim->refines->items = g_slist_reverse (sim->refines->items);
  sim->adapts->items = g_slist_reverse (sim->adapts->items);
  sim->events->items = g_slist_reverse (sim->events->items);
  sim->solids->items = g_slist_reverse (sim->solids->items);
  sim->modules = g_slist_reverse (sim->modules);

  if (sim->advection_params.scheme == GFS_NONE)
    domain->advection_metric = NULL;
  else if (domain->advection_metric && !gfs_has_source_coriolis (domain)) {
    /* we need to add a pseudo-Coriolis term to include the metric advection terms */
    
  }
}

/**
 * gfs_advance_tracers:
 * @sim: a #GfsSimulation.
 * @dt: the timestep.
 *
 * Performs advection/difussion of tracers associated with @sim.
 */
void gfs_advance_tracers (GfsSimulation * sim, gdouble dt)
{
  g_return_if_fail (sim != NULL);

  GfsDomain * domain = GFS_DOMAIN (sim);
  GSList * i = sim->events->items;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER_VOF (i->data)) {
      GfsVariableTracer * t = i->data;
      
      t->advection.dt = dt;
      gfs_tracer_vof_advection (domain, &t->advection);
      gfs_domain_variable_centered_sources (domain, i->data, i->data, t->advection.dt);
    }
    else if (GFS_IS_VARIABLE_TRACER (i->data)) {
      GfsVariableTracer * t = i->data;
      
      t->advection.dt = dt;
      gfs_tracer_advection_diffusion (domain, &t->advection, sim->physical_params.alpha);
      gfs_domain_cell_traverse (domain,
				FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) GFS_VARIABLE (t)->fine_coarse, t);
    }
    i = i->next;
  }  
}

static void simulation_run (GfsSimulation * sim)
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

  gfs_simulation_set_timestep (sim);
  if (sim->time.i == 0) {
    gfs_approximate_projection (domain,
				&sim->approx_projection_params,
				sim->advection_params.dt,
				p, sim->physical_params.alpha, res, g, NULL);
    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt/2.);
  }
  else if (sim->advection_params.gc)
    gfs_update_gradients (domain, p, sim->physical_params.alpha, g);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

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

static gdouble simulation_cfl (GfsSimulation * sim)
{
  GSList * i = GFS_DOMAIN (sim)->variables;
  gdouble cflmin = G_MAXDOUBLE;
  
  while (i) {
    GfsVariable * v = i->data;

    if (GFS_IS_VARIABLE_TRACER (v) && 
 	GFS_VARIABLE_TRACER (v)->advection.scheme != GFS_NONE &&
	GFS_VARIABLE_TRACER (v)->advection.sink[0]) {
      gfs_add_sinking_velocity (GFS_DOMAIN (sim), &GFS_VARIABLE_TRACER (v)->advection);
      gdouble cfl = gfs_domain_cfl (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1);
      gfs_remove_sinking_velocity (GFS_DOMAIN (sim), &GFS_VARIABLE_TRACER (v)->advection);
      if (cfl < cflmin)
	cflmin = cfl;
    }
    i = i->next;
  }
  return cflmin < G_MAXDOUBLE ? cflmin : gfs_domain_cfl (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1);
}

static void gfs_simulation_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->write =   simulation_write;
  GTS_OBJECT_CLASS (klass)->read =    simulation_read;
  GTS_OBJECT_CLASS (klass)->destroy = simulation_destroy;

  klass->run = simulation_run;
  klass->cfl = simulation_cfl;
}

/* Derived variables */

static gdouble cell_x (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    gfs_face_ca (face, &p);
  else
    gfs_cell_cm (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.x;
}

static gdouble cell_y (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    gfs_face_ca (face, &p);
  else
    gfs_cell_cm (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.y;
}

static gdouble cell_z (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    gfs_face_ca (face, &p);
  else
    gfs_cell_cm (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.z;
}

static gdouble cell_ax (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  g_return_val_if_fail (cell != NULL, 0.);
  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    FttVector p = GFS_STATE (cell)->solid->ca;
    gfs_simulation_map_inverse (sim, &p);
    return p.x;
  }
}

static gdouble cell_ay (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  g_return_val_if_fail (cell != NULL, 0.);
  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    FttVector p = GFS_STATE (cell)->solid->ca;
    gfs_simulation_map_inverse (sim, &p);
    return p.y;
  }
}

static gdouble cell_az (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  g_return_val_if_fail (cell != NULL, 0.);
  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    FttVector p = GFS_STATE (cell)->solid->ca;
    gfs_simulation_map_inverse (sim, &p);
    return p.z;
  }
}

static gdouble cell_cx (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.x;
}

static gdouble cell_cy (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.y;
}

static gdouble cell_cz (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  gfs_simulation_map_inverse (sim, &p);
  return p.z;
}

static gdouble cell_rx (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  return p.x;
}

static gdouble cell_ry (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  return p.y;
}

static gdouble cell_rz (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  FttVector p;

  g_return_val_if_fail (cell != NULL || face != NULL, 0.);

  if (face)
    ftt_face_pos (face, &p);
  else
    ftt_cell_pos (cell, &p);
  return p.z;
}

static gdouble cell_dV (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  gdouble dV = gfs_cell_volume (cell, GFS_DOMAIN (sim));
  gdouble L = sim->physical_params.L;
#if FTT_2D
  dV *= L*L;
#else
  dV *= L*L*L;
#endif
  return dV;
}

static gdouble cell_dL (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  gdouble dL = ftt_cell_size (cell);
  gdouble L = sim->physical_params.L;
  return L*dL;
}

static gdouble cell_t (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  return sim->time.t;
}

static gdouble cell_dt (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  return sim->advection_params.dt;
}

static gdouble cell_vorticity (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return gfs_vorticity (cell, gfs_domain_velocity (domain));
}

static gdouble cell_divergence (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return gfs_divergence (cell, gfs_domain_velocity (domain));
}

static gdouble cell_velocity_norm (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble L = GFS_SIMULATION (domain)->physical_params.L;
  return L*gfs_vector_norm (cell, gfs_domain_velocity (domain));
}

static gdouble cell_velocity_norm2 (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble L = GFS_SIMULATION (domain)->physical_params.L;
  return L*L*gfs_vector_norm2 (cell, gfs_domain_velocity (domain));
}

static gdouble cell_level (FttCell * cell)
{
  return ftt_cell_level (cell);
}

static gdouble cell_fraction (FttCell * cell)
{
  g_return_val_if_fail (cell != NULL, 0.);
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
}

static gdouble cell_solid_area (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  FttVector n;
  gfs_solid_normal (cell, &n);
  return ftt_vector_norm (&n)*pow (ftt_cell_size (cell)*GFS_SIMULATION (domain)->physical_params.L,
				   FTT_DIMENSION - 1);
}

static gdouble cell_solid_sr (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_RIGHT] : 1.;
}

static gdouble cell_solid_sl (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_LEFT] : 1.;
}

static gdouble cell_solid_st (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_TOP] : 1.;
}

static gdouble cell_solid_sb (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_BOTTOM] : 1.;
}

#if !FTT_2D
static gdouble cell_solid_sf (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_FRONT] : 1.;
}

static gdouble cell_solid_sk (FttCell * cell)
{
  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_BACK] : 1.;
}
#endif /* 3D */

static gdouble cell_velocity_lambda2 (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return gfs_vector_lambda2 (cell, gfs_domain_velocity (domain));
}

static gdouble cell_streamline_curvature (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble L = GFS_SIMULATION (domain)->physical_params.L;
  return gfs_streamline_curvature (cell, gfs_domain_velocity (domain))/L;
}

static gdouble cell_2nd_principal_invariant (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return gfs_2nd_principal_invariant (cell, gfs_domain_velocity (domain))/ftt_cell_size (cell);
}

static gdouble cell_pid (FttCell * cell)
{
  while (!FTT_CELL_IS_ROOT (cell))
    cell = ftt_cell_parent (cell);
  return GFS_BOX (FTT_ROOT_CELL (cell)->parent)->pid;
}

static gdouble cell_id (FttCell * cell)
{
  while (!FTT_CELL_IS_ROOT (cell))
    cell = ftt_cell_parent (cell);
  return GFS_BOX (FTT_ROOT_CELL (cell)->parent)->id;
}

static gdouble cell_orthogonality (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble h;
  FttVector p;
  if (cell) {
    ftt_cell_pos (cell, &p);
    h = ftt_cell_size (cell)/2.;
  }
  else {
    ftt_face_pos (face, &p);
    h = ftt_cell_size (face->cell)/2.;
  }
  FttVector p1 = p, p2 = p, p3 = p, p4 = p;
  p1.x += h; p2.x -= h;
  p3.y += h; p4.y -= h;
  gfs_simulation_map_inverse (GFS_SIMULATION (domain), &p1);
  gfs_simulation_map_inverse (GFS_SIMULATION (domain), &p2);
  gfs_simulation_map_inverse (GFS_SIMULATION (domain), &p3);
  gfs_simulation_map_inverse (GFS_SIMULATION (domain), &p4);
  GtsVector p1p2, p3p4;
  gts_vector_init (p1p2, &p1, &p2);
  gts_vector_init (p3p4, &p3, &p4);
  return fabs (acos (gts_vector_scalar (p1p2, p3p4)/(gts_vector_norm (p1p2)*gts_vector_norm (p3p4)))
	       - M_PI/2.)*180./M_PI;
}

static void simulation_init (GfsSimulation * object)
{
  GfsDomain * domain = GFS_DOMAIN (object);
  /* Please update http://gfs.sourceforge.net/wiki/index.php/Domain_variables 
     when changing this list */
  static GfsDerivedVariableInfo derived_variable[] = {
    { "x", "x-coordinate of the center of mass of the cell", cell_x },
    { "y", "y-coordinate of the center of mass of the cell", cell_y },
    { "z", "z-coordinate of the center of mass of the cell", cell_z },
    { "ax", "x-coordinate of the center of area of the solid surface", cell_ax },
    { "ay", "y-coordinate of the center of area of the solid surface", cell_ay },
    { "az", "z-coordinate of the center of area of the solid surface", cell_az },
    { "cx", "x-coordinate of the center of the cell", cell_cx },
    { "cy", "y-coordinate of the center of the cell", cell_cy },
    { "cz", "z-coordinate of the center of the cell", cell_cz },
    { "rx", "x-coordinate of the center of the cell (internal)", cell_rx },
    { "ry", "y-coordinate of the center of the cell (internal)", cell_ry },
    { "rz", "z-coordinate of the center of the cell (internal)", cell_rz },
    { "dV", "volume of the cell", cell_dV},
    { "dL", "length of the cell", cell_dL},
    { "t",  "Physical time", cell_t },
    { "dt", "Timestep", cell_dt },
    { "Vorticity", "Norm of the vorticity vector of the velocity field", cell_vorticity },
    { "Divergence", "Divergence of the velocity field", cell_divergence },
    { "Velocity", "Norm of the velocity vector", cell_velocity_norm },
    { "Velocity2", "Squared norm of the velocity vector", cell_velocity_norm2 },
    { "Level", "Quad/octree level of the cell", cell_level },
    { "A", "Fluid fraction of the cell", cell_fraction },
    { "S", "Area of the solid contained in the cell", cell_solid_area },
    { "Sr", "Fluid fraction of the right cell face", cell_solid_sr },
    { "Sl", "Fluid fraction of the left cell face", cell_solid_sl },
    { "St", "Fluid fraction of the top cell face", cell_solid_st },
    { "Sb", "Fluid fraction of the bottom cell face", cell_solid_sb },
#if !FTT_2D
    { "Sf", "Fluid fraction of the front cell face", cell_solid_sf },
    { "Sk", "Fluid fraction of the back cell face", cell_solid_sk },    
#endif /* 3D */
    { "Lambda2", "Vortex-detection criterion of Jeong & Hussein", cell_velocity_lambda2 },
    { "Curvature",  "Curvature of the local streamline", cell_streamline_curvature },
    { "D2", "Second principal invariant of the deformation tensor", cell_2nd_principal_invariant },
    { "Pid", "Parent box process ID", cell_pid },
    { "Id", "Parent box ID", cell_id },
    { "Orthogonality", "Deviation from orthogonality", cell_orthogonality },
    { NULL, NULL, NULL}
  };

  GfsVariable * v;
  v = gfs_domain_add_variable (domain, "P",    "Approximate projection pressure");
  v->centered = TRUE;
  v->units = 2.;
  v = gfs_domain_add_variable (domain, "Pmac", "MAC projection pressure");
  v->centered = TRUE;  
  v->units = 2.;
  GfsVariable * u[FTT_DIMENSION];
  u[0] = gfs_domain_add_variable (domain, "U",    "x-component of the velocity");
  u[0]->units = 1.;
  u[0]->face_source = TRUE;
  u[1] = gfs_domain_add_variable (domain, "V",    "y-component of the velocity");
  u[1]->units = 1.;
  u[1]->face_source = TRUE;
#if (!FTT_2D)
  u[2] = gfs_domain_add_variable (domain, "W",    "z-component of the velocity");
  u[2]->units = 1.;
  u[2]->face_source = TRUE;
#endif /* FTT_3D */
  gfs_variable_set_vector (u, FTT_DIMENSION);

  GfsDerivedVariableInfo * dv = derived_variable;
  while (dv->name) {
    g_assert (gfs_domain_add_derived_variable (domain, *dv));
    dv++;
  }
  domain->derived_variables = g_slist_reverse (domain->derived_variables);

  gfs_time_init (&object->time);
  gfs_physical_params_init (&object->physical_params);

  gfs_advection_params_init (&object->advection_params);
  object->advection_params.flux = gfs_face_velocity_advection_flux;
  object->advection_params.average = TRUE;

  gfs_multilevel_params_init (&object->projection_params);
  gfs_multilevel_params_init (&object->approx_projection_params);

  object->solids = GTS_SLIST_CONTAINER (gts_container_new
					(GTS_CONTAINER_CLASS
					 (gts_slist_container_class ())));
  object->output_solid = TRUE;

  object->refines = GTS_SLIST_CONTAINER (gts_container_new
					 (GTS_CONTAINER_CLASS
					  (gts_slist_container_class ())));
  object->maps = GTS_SLIST_CONTAINER (gts_container_new
				      (GTS_CONTAINER_CLASS
				       (gts_slist_container_class ())));
  object->adapts = GTS_SLIST_CONTAINER (gts_container_new
					(GTS_CONTAINER_CLASS
					 (gts_slist_container_class ())));
  gfs_adapt_stats_init (&object->adapts_stats);
  object->events = GTS_SLIST_CONTAINER (gts_container_new
					(GTS_CONTAINER_CLASS
					 (gts_slist_container_class ())));
  object->modules = object->preloaded_modules = NULL;
  
  object->tnext = 0.;
}

GfsSimulationClass * gfs_simulation_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_simulation_info = {
      "GfsSimulation",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_simulation_class_init,
      (GtsObjectInitFunc) simulation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_domain_class ()),
				  &gfs_simulation_info);
  }

  return klass;
}

GfsSimulation * gfs_simulation_new (GfsSimulationClass * klass)
{
  GfsSimulation * object;

  object = GFS_SIMULATION (gts_graph_new (GTS_GRAPH_CLASS (klass),
					  GTS_GNODE_CLASS (gfs_box_class ()),
					  GTS_GEDGE_CLASS (gfs_gedge_class ())));

  return object;
}

static void init_non_variable (GfsEvent * event, GfsSimulation * sim)
{
  if (!GFS_IS_VARIABLE (event))
    gfs_event_init (event, sim);
}

static void redo_non_variable (GfsEvent * event, GfsSimulation * sim)
{
  if (!GFS_IS_VARIABLE (event))
    gfs_event_redo (event, sim);
}

/**
 * gfs_simulation_init:
 * @sim: a #GfsSimulation.
 *
 * Initialises @sim: matches boundary conditions, applies boundary
 * conditions and initialises all variables, etc...
 */
void gfs_simulation_init (GfsSimulation * sim)
{
  g_return_if_fail (sim != NULL);

  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) init_non_variable, sim);

  GfsDomain * domain = GFS_DOMAIN (sim);
  gfs_domain_match (domain);
  gfs_set_merged (domain);
  GSList * i = domain->variables;
  while (i) {
    gfs_event_init (GFS_EVENT (i->data), sim);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, i->data);
    i = i->next;
  }
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_coarse_init, domain);

  /* Iterate mesh adaptation and initialisation */
  if (domain->timer->start >= 0.) /* only if timing is on (i.e. within gfs_simulation_run()) */
    while (gfs_simulation_adapt (sim)) {
      gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) redo_non_variable, sim);
      i = domain->variables;
      while (i) {
	gfs_event_redo (GFS_EVENT (i->data), sim);
	gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, i->data);
	i = i->next;
      }
      gfs_domain_cell_traverse (domain,
				FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) gfs_cell_coarse_init, domain);
      if (sim->adapts_stats.depth_increase <= 0)
	break;
    }
}

static void refine_cell_corner (FttCell * cell, GfsDomain * domain)
{
  if (FTT_CELL_IS_LEAF (cell) && ftt_refine_corner (cell))
    ftt_cell_refine_single (cell, domain->cell_init, domain->cell_init_data);
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

/**
 * gfs_check_solid_fractions:
 * @domain: a #GfsDomain.
 *
 * Checks consistency of solid fractions.
 *
 * Returns: the number of partial fluid faces which do not have two
 * neighbors.
 */
guint gfs_check_solid_fractions (GfsDomain * domain)
{
  g_return_val_if_fail (domain != NULL, 0);

  guint nf = 0;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) check_solid_fractions, &nf);
  return nf;
}

static void is_diffusion (GfsSource * s, gboolean * diffusion)
{
  *diffusion = (GFS_IS_SOURCE_DIFFUSION (s) != NULL);
}

static void set_permanent (FttCell * cell)
{
  cell->flags |= GFS_FLAG_PERMANENT;
}

/**
 * gfs_simulation_get_solids:
 * @sim: a #GfsSimulation.
 *
 * Returns: a new list of #GfsSolid defining the solid boundaries
 * contained in @sim.
 */
GSList * gfs_simulation_get_solids (GfsSimulation * sim)
{
  g_return_val_if_fail (sim != NULL, NULL);

  GSList * solids = NULL, * i = sim->solids->items;
  while (i) {
    solids = g_slist_prepend (solids, i->data);
    i = i->next;
  }
  return solids;
}

static gboolean coarsen_cell (void)
{
  return TRUE;
}

static void refine_leaf_boxes (GfsBox * box, GfsDomain * domain)
{
  /* refine then coarsen boxes which are also leaf cells: this is the
   * simplest way to initialise the values of "automatic" GfsVariables
   * e.g. VariableTerrain, VariableMetric etc... */
  if (FTT_CELL_IS_LEAF (box->root)) {
    ftt_cell_refine_single (box->root, domain->cell_init, domain->cell_init_data);
    gfs_cell_coarse_init (box->root, domain);
    ftt_cell_coarsen (box->root, 
		      (FttCellCoarsenFunc) coarsen_cell, NULL,
		      (FttCellCleanupFunc) gfs_cell_cleanup, domain);
  }
}

/**
 * gfs_simulation_refine:
 * @sim: a #GfsSimulation.
 *
 * Calls the @refine() methods of the #GfsRefine of @sim. Initialize the
 * solid fractions by calling gfs_init_solid_fractions(). Matches the
 * boundaries by calling gfs_domain_match().
 */
void gfs_simulation_refine (GfsSimulation * sim)
{
  GSList * i;
  guint depth;
  gint l;
  GfsDomain * domain;

  g_return_if_fail (sim != NULL);

  domain = GFS_DOMAIN (sim);

  gfs_domain_timer_start (domain, "simulation_refine");
  i = sim->refines->items;
  while (i) {
    GfsRefine * refine = i->data;
    GSList * next = i->next;    
    if (GFS_REFINE_CLASS (GTS_OBJECT (refine)->klass)->refine)
      (* GFS_REFINE_CLASS (GTS_OBJECT (refine)->klass)->refine) (refine, sim);
    i = next;
  }

  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) refine_leaf_boxes, sim);

  depth = gfs_domain_depth (domain);
  for (l = depth - 2; l >= 0; l--)
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
			      (FttCellTraverseFunc) refine_cell_corner, 
			      domain);

  gfs_domain_match (domain);
  gfs_domain_timer_stop (domain, "simulation_refine");

  GSList * solids = gfs_simulation_get_solids (sim);
  if (solids) {
    gfs_domain_timer_start (domain, "solid_fractions");
    sim->thin = gfs_domain_init_solid_fractions (domain, solids, TRUE,
						 (FttCellCleanupFunc) gfs_cell_cleanup, domain,  
						 NULL);
    g_slist_free (solids);
    gfs_domain_match (domain);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) set_permanent, NULL);
    gfs_domain_timer_stop (domain, "solid_fractions");
  }

  guint nf = gfs_check_solid_fractions (domain);
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

/**
 * gfs_simulation_read:
 * @fp: a #GtsFile.
 *
 * Reads a simulation file from @fp.
 *
 * Returns: the #GfsSimulation or %NULL if an error occured, in which
 * case the @pos and @error fields of @fp are set.
 */
GfsSimulation * gfs_simulation_read (GtsFile * fp)
{
  GfsDomain * d;
  GSList * ml = NULL; /* list of preloaded modules */

  g_return_val_if_fail (fp != NULL, NULL);

  while (fp->type == '\n')
     gts_file_next_token (fp);

  while (fp->type == GTS_STRING && !strcmp (fp->token->str, "GModule")) { /* preloaded module */
    GModule * module = load_module (fp, NULL);
    if (module == NULL)
      return NULL;
    ml = g_slist_prepend (ml, module);
    while (fp->type == '\n')
      gts_file_next_token (fp);
  }
      
  d = gfs_domain_read (fp);
  if (d != NULL && !GFS_IS_SIMULATION (d)) {
    gts_file_error (fp, "parent graph is not a GfsSimulation");
    gts_object_destroy (GTS_OBJECT (d));
    g_slist_free (ml);
    return NULL;
  }
  if (d) {
    gfs_pending_functions_compilation (fp);
    if (fp->type == GTS_ERROR) {
      gts_object_destroy (GTS_OBJECT (d));
      g_slist_free (ml);
      return NULL;
    }
    else
      GFS_SIMULATION (d)->preloaded_modules = g_slist_reverse (ml);
  }
  else
    g_slist_free (ml);
  return GFS_SIMULATION (d);
}

static void write_preloaded_modules (GfsSimulation * sim, FILE * fp)
{
  GSList * i = sim->preloaded_modules;
  while (i) {
    void (* module_write) (FILE *);
    const gchar * name = NULL;
    fprintf (fp, "GModule %s", 
	     g_module_symbol (i->data, "gfs_module_name", (gpointer) &name) ? name : 
	     g_module_name (i->data));
    if (g_module_symbol (i->data, "gfs_module_write", (gpointer) &module_write))
      (* module_write) (fp);
    fputc ('\n', fp);
    i = i->next;
  }
}

/**
 * gfs_simulation_write:
 * @sim: a #GfsSimulation.
 * @max_depth: the maximum depth at which to stop writing cell tree
 * data (-1 means no limit).
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of @sim. If @max_depth is
 * smaller or equal to -2, no cell tree data is written.  
 */
void gfs_simulation_write (GfsSimulation * sim,
			   gint max_depth,		  
			   FILE * fp)
{
  gint depth;
  GfsDomain * domain;

  g_return_if_fail (sim != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp, "# Gerris Flow Solver %dD version %s (%s)\n",
	   FTT_DIMENSION, GFS_VERSION, GFS_BUILD_VERSION);
  write_preloaded_modules (sim, fp);
  domain = GFS_DOMAIN (sim);
  depth = domain->max_depth_write;
  domain->max_depth_write = max_depth;
  gts_graph_write (GTS_GRAPH (sim), fp);
  domain->max_depth_write = depth;
}

#ifdef HAVE_MPI
static void count_edges (GtsGEdge * e, guint * nedge)
{
  (*nedge)++;
}

static void write_node (GtsObject * node, gpointer * data)
{
  FILE * fp = data[0];
  guint * nnode = data[1];

  node->reserved = GUINT_TO_POINTER ((*nnode)++);
  if (node->klass->write)
    (* node->klass->write) (node, fp);
  fputc ('\n', fp);
}

static void write_edge (GtsGEdge * edge, FILE * fp)
{
  fprintf (fp, "%u %u", 
	   GPOINTER_TO_UINT (GTS_OBJECT (edge->n1)->reserved),
	   GPOINTER_TO_UINT (GTS_OBJECT (edge->n2)->reserved));
  if (GTS_OBJECT (edge)->klass->write)
    (* GTS_OBJECT (edge)->klass->write) (GTS_OBJECT (edge), fp);
  fputc ('\n', fp);
}
#endif /* HAVE_MPI */

/**
 * gfs_simulation_union_write:
 * @sim: a #GfsSimulation.
 * @max_depth: the maximum depth at which to stop writing cell tree
 * data (-1 means no limit).
 * @fp: a file pointer.
 *
 * Identical to gfs_simulation_write() for serial simulations. For
 * parallel simulations writes the union of the simulations on all
 * processes to @fp.
 */
void gfs_simulation_union_write (GfsSimulation * sim,
				 gint max_depth,		  
				 FILE * fp)
{
  GfsDomain * domain = GFS_DOMAIN (sim);

  g_return_if_fail (sim != NULL);
  g_return_if_fail (fp != NULL);

  if (domain->pid < 0)
    gfs_simulation_write (sim, max_depth, fp);
  else {
#ifdef HAVE_MPI
    int gsize;
    guint * nbox;

    MPI_Comm_size (MPI_COMM_WORLD, &gsize);
    nbox = g_malloc (sizeof (guint)*gsize);
    guint n = gts_container_size (GTS_CONTAINER (sim));
    MPI_Allgather (&n, 1, MPI_UNSIGNED, nbox, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    /* nbox[] now contains the number of boxes on each PE */

    /* see gts/src/graph.c:gts_graph_write() for the original (serial) implementation */
    GtsGraph * g = GTS_GRAPH (sim);
    guint nedge = 0;
    gts_graph_foreach_edge (g, (GtsFunc) count_edges, &nedge);
    gfs_all_reduce (domain, nedge, MPI_UNSIGNED, MPI_SUM);

    if (domain->pid == 0) {
      fprintf (fp, "# Gerris Flow Solver %dD version %s (%s)\n",
	       FTT_DIMENSION, GFS_VERSION, GFS_BUILD_VERSION);
      write_preloaded_modules (sim, fp);
      guint i, nboxes = 0;
      for (i = 0; i < gsize; i++)
	nboxes += nbox[i];
      fprintf (fp, "%u %u", nboxes, nedge);
      if (GTS_OBJECT (g)->klass->write)
	(* GTS_OBJECT (g)->klass->write) (GTS_OBJECT (g), fp);
      fputc ('\n', fp);
    }

    gint depth = domain->max_depth_write;
    guint i, nnode = 1;
    gpointer data[2];

    for (i = 0; i < domain->pid; i++)
      nnode += nbox[i];
    g_free (nbox);

    GfsUnionFile uf;
    FILE * fpp = gfs_union_open (fp, domain->pid, &uf);
    data[0] = fpp;
    data[1] = &nnode;
    domain->max_depth_write = max_depth;
    gts_container_foreach (GTS_CONTAINER (g), (GtsFunc) write_node, data);
    domain->max_depth_write = depth;
    gfs_union_close (fp, domain->pid, &uf);

    fpp = gfs_union_open (fp, domain->pid, &uf);
    gts_graph_foreach_edge (g, (GtsFunc) write_edge, fpp);
    gfs_union_close (fp, domain->pid, &uf);

    gts_container_foreach (GTS_CONTAINER (g), (GtsFunc) gts_object_reset_reserved, NULL);
#endif /* HAVE_MPI */
  }
}

static gdouble min_cfl (GfsSimulation * sim)
{
  gdouble cfl = (sim->advection_params.scheme == GFS_NONE ?
		 G_MAXDOUBLE :
		 sim->advection_params.cfl);
  GSList * i = GFS_DOMAIN (sim)->variables;

  while (i) {
    GfsVariable * v = i->data;

    if (GFS_IS_VARIABLE_TRACER (v) && 
	GFS_VARIABLE_TRACER (v)->advection.scheme != GFS_NONE &&
	GFS_VARIABLE_TRACER (v)->advection.cfl < cfl)
      cfl = GFS_VARIABLE_TRACER (v)->advection.cfl;
    i = i->next;
  }
  return cfl;
}

/**
 * gfs_simulation_set_timestep:
 * @sim: a #GfsSimulation.
 *
 * Sets the time step for the next iteration of @sim using the CFL
 * (computed using gfs_domain_cfl()), the stability conditions for
 * source terms and taking into account the timings of the various
 * #GfsEvent associated to @sim.
 *
 * More precisely, the time step is adjusted (if necessary) so that
 * the time of the closest event is exactly reached after the
 * iteration.  
 */
void gfs_simulation_set_timestep (GfsSimulation * sim)
{
  gdouble t, cfl;

  g_return_if_fail (sim != NULL);

  t = sim->time.t;
  if ((cfl = min_cfl (sim)) < G_MAXDOUBLE)
    sim->advection_params.dt = cfl*(* GFS_SIMULATION_CLASS (GTS_OBJECT (sim)->klass)->cfl) (sim);
  else
    sim->advection_params.dt = G_MAXINT;
  if (sim->advection_params.dt > sim->time.dtmax)
    sim->advection_params.dt = sim->time.dtmax;

  GSList *  i = GFS_DOMAIN (sim)->variables;
  while (i) {
    GfsVariable * v = i->data;
    if (v->sources) {
      GSList * j = GTS_SLIST_CONTAINER (v->sources)->items;
      while (j) {
	GfsSourceGeneric * s = j->data;
	if (GFS_SOURCE_GENERIC_CLASS (GTS_OBJECT (s)->klass)->stability) {
	  gdouble dt = (* GFS_SOURCE_GENERIC_CLASS (GTS_OBJECT (s)->klass)->stability) (s, sim);
	  if (dt < sim->advection_params.dt)
	    sim->advection_params.dt = dt;
	}
	j = j->next;
      }
    }
    i = i->next;
  }

  gfs_all_reduce (GFS_DOMAIN (sim), sim->advection_params.dt, MPI_DOUBLE, MPI_MIN);

  gdouble tnext = G_MAXINT;
  i = sim->events->items;
  while (i) {
    gdouble next = gfs_event_next (i->data, sim);
    if (t < next && next < tnext)
      tnext = next + 1e-9;
    i = i->next;
  }
  if (sim->time.end < tnext)
    tnext = sim->time.end;

  gdouble n = ceil ((tnext - t)/sim->advection_params.dt);
  if (n > 0. && n < G_MAXINT) {
    sim->advection_params.dt = (tnext - t)/n;
    if (n == 1.)
      sim->tnext = tnext;
    else
      sim->tnext = t + sim->advection_params.dt;
  }
  else
    sim->tnext = t + sim->advection_params.dt;

  if (sim->advection_params.dt < 1e-9)
    sim->advection_params.dt = 1e-9;

  if (sim->time.t < sim->time.end && sim->advection_params.dt == G_MAXINT)
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_WARNING,
	   "could not find a suitable timescale to set the timestep.\n"
	   "Please set a timescale explicitly (e.g. through the 'step' parameter\n"
	   "of a GfsEvent).");
}

/**
 * gfs_time_write:
 * @t: the time structure.
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of the time structure @t.
 */
void gfs_time_write (GfsTime * t, FILE * fp)
{
  g_return_if_fail (t != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp, "{ i = %u t = %g ", t->i, t->t);
  if (t->start != 0.)
    fprintf (fp, "start = %g ", t->start);
  if (t->istart != 0)
    fprintf (fp, "start = %u ", t->istart);
  if (t->end < G_MAXDOUBLE)
    fprintf (fp, "end = %g ", t->end);
  if (t->iend < G_MAXINT)
    fprintf (fp, "iend = %u ", t->iend);
  if (t->dtmax < G_MAXDOUBLE)
    fprintf (fp, "dtmax = %g ", t->dtmax);
  fputc ('}', fp);
}

/**
 * gfs_time_init:
 * @t: the #GfsTime.
 *
 * Initializes the time structure @t with default values.
 */
void gfs_time_init (GfsTime * t)
{
  g_return_if_fail (t != NULL);
  
  t->t = t->start = 0.;
  t->end = G_MAXDOUBLE;

  t->i = t->istart = 0;
  t->iend = G_MAXINT;

  t->dtmax = G_MAXDOUBLE;
}

/**
 * gfs_time_read:
 * @t: the #GfsTime.
 * @fp: the #GtsFile.
 *
 * Reads a time structure from @fp and puts it in @t.
 */
void gfs_time_read (GfsTime * t, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "t",      TRUE},
    {GTS_DOUBLE, "start",  TRUE},
    {GTS_DOUBLE, "end",    TRUE},
    {GTS_UINT,   "i",      TRUE},
    {GTS_UINT,   "istart", TRUE},
    {GTS_UINT,   "iend",   TRUE},
    {GTS_DOUBLE, "dtmax",  TRUE},
    {GTS_NONE}
  };

  g_return_if_fail (t != NULL);
  g_return_if_fail (fp != NULL);

  var[0].data = &t->t;
  var[1].data = &t->start;
  var[2].data = &t->end;
  var[3].data = &t->i;
  var[4].data = &t->istart;
  var[5].data = &t->iend;
  var[6].data = &t->dtmax;

  gts_file_assign_variables (fp, var);

  if (t->t < t->start)
    t->t = t->start;
  if (t->i < t->istart)
    t->i = t->istart;
}

/**
 * gfs_physical_params_write:
 * @p: the physical parameters structure.
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of the physical parameters
 * structure @p.  
 */
void gfs_physical_params_write (GfsPhysicalParams * p, FILE * fp)
{
  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp, "{ g = %g L = %g", p->g, p->L);
  if (p->alpha) {
    fputs (" alpha =", fp);
    gfs_function_write (p->alpha, fp);
  }
  fputs (" }", fp);
}

/**
 * gfs_physical_params_init:
 * @p: the #GfsPhysicalParams.
 *
 * Initializes the physical parameters structure @p with default values.
 */
void gfs_physical_params_init (GfsPhysicalParams * p)
{
  g_return_if_fail (p != NULL);
  
  p->g = p->L = 1.;
  p->alpha = NULL;
}

/**
 * gfs_physical_params_read:
 * @p: the #GfsPhysicalParams.
 * @domain: a #GfsDomain.
 * @fp: the #GtsFile.
 *
 * Reads a physical parameters structure from @fp and puts it in @p.
 */
void gfs_physical_params_read (GfsPhysicalParams * p, GfsDomain * domain, GtsFile * fp)
{
  g_return_if_fail (p != NULL);
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

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
    else {
      gchar * id = g_strdup (fp->token->str);

      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);

      if (!strcmp (id, "g")) {
	/* fixme: units? */
	p->g = gfs_read_constant (fp, domain);
	if (fp->type == GTS_ERROR) {
	  g_free (id);
	  return;
	}
      }
      else if (!strcmp (id, "L")) {
	p->L = gfs_read_constant (fp, domain);
	if (fp->type == GTS_ERROR) {
	  g_free (id);
	  return;
	}
	if (p->L == 0.) {
	  g_free (id);
	  gts_file_error (fp, "L must be different from zero");
	  return;
	}
      }
      else if (!strcmp (id, "alpha")) {
	p->alpha = gfs_function_new (gfs_function_class (), 0.);
	gfs_function_read (p->alpha, domain, fp);
	if (fp->type == GTS_ERROR) {
	  g_free (id);
	  gts_object_destroy (GTS_OBJECT (p->alpha));
	  return;
	}
      }
      else {
	g_free (id);
	gts_file_error (fp, "unknown keyword `%s'", id);
	return;
      }
      g_free (id);
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void error_handler (const gchar *log_domain,
			   GLogLevelFlags log_level,
			   const gchar *message,
			   gpointer user_data)
{
  GfsDomain * domain = user_data;
  g_slist_free (domain->variables_io);
  domain->variables_io = NULL;
  GSList * i = domain->variables;
  while (i) {
    if (GFS_VARIABLE (i->data)->name)
      domain->variables_io = g_slist_append (domain->variables_io, i->data);
    i = i->next;
  }
  gchar fname[20] = "error.gfs";
  if (domain->pid >= 0)
    snprintf (fname, 20, "error-%d.gfs", domain->pid);
  FILE * fp = fopen (fname, "w");
  if (fp) {
    gfs_simulation_write (GFS_SIMULATION (domain), -1, fp);
    fclose (fp);
  }

  g_log_default_handler (log_domain, log_level, message, NULL);
}

/**
 * gfs_simulation_run:
 * @sim: a #GfsSimulation.
 *
 * Runs @sim.
 */
void gfs_simulation_run (GfsSimulation * sim)
{
  g_return_if_fail (sim != NULL);

  guint id = g_log_set_handler ("Gfs", G_LOG_LEVEL_ERROR | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION,
				error_handler, sim);
  GfsDomain * domain = GFS_DOMAIN (sim);
  g_timer_start (domain->clock);
  gfs_clock_start (domain->timer);
  gts_range_init (&domain->mpi_wait);
  (* GFS_SIMULATION_CLASS (GTS_OBJECT (sim)->klass)->run) (sim);
  gfs_clock_stop (domain->timer);
  g_timer_stop (domain->clock);
  g_log_remove_handler ("Gfs", id);
}

/**
 * gfs_simulation_map:
 * @sim: a #GfsSimulation.
 * @p: a #FttVector.
 *
 * Maps the physical space coordinates @p into computational space.
 */
void gfs_simulation_map (GfsSimulation * sim, FttVector * p)
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (p != NULL);
  
  FttComponent c;
  for (c = 0; c < 3; c++)
    (&p->x)[c] *= (&GFS_DOMAIN (sim)->lambda.x)[c]/sim->physical_params.L;

  GSList * i = sim->maps->items;
  while (i) {
    (* GFS_MAP (i->data)->transform) (i->data, p, p);
    i = i->next;
  }
}

/**
 * gfs_simulation_map_inverse:
 * @sim: a #GfsSimulation.
 * @p: a #FttVector.
 *
 * Maps the computational coordinates @p into physical space.
 */
void gfs_simulation_map_inverse (GfsSimulation * sim, FttVector * p)
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (p != NULL);
  
  GSList * reverse = g_slist_reverse (sim->maps->items);
  GSList * i = reverse;
  while (i) {
    (* GFS_MAP (i->data)->inverse) (i->data, p, p);
    i = i->next;
  }
  sim->maps->items = g_slist_reverse (reverse);

  FttComponent c;
  for (c = 0; c < 3; c++)
    (&p->x)[c] *= sim->physical_params.L/(&GFS_DOMAIN (sim)->lambda.x)[c];
}

/**
 * gfs_simulation_map_vector:
 * @sim: a #GfsSimulation.
 * @p: the position.
 * @v: a #FttVector.
 *
 * Applies the mapping transformations associated with @sim at
 * location @p, to vector @v.
 */
void gfs_simulation_map_vector (GfsSimulation * sim, const FttVector * p, FttVector * v)
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (v != NULL);

  GSList * i = sim->maps->items;
  while (i) {
    (* GFS_MAP (i->data)->transform_vector) (i->data, p, v, v);
    i = i->next;
  }
}

/**
 * gfs_simulation_map_inverse_vector:
 * @sim: a #GfsSimulation.
 * @p: the position.
 * @v: a #FttVector.
 *
 * Applies the inverse mapping transformations associated with @sim at
 * location @p, to vector @v.
 */
void gfs_simulation_map_inverse_vector (GfsSimulation * sim, const FttVector * p, FttVector * v)
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (v != NULL);
  
  GSList * reverse = g_slist_reverse (sim->maps->items);
  GSList * i = reverse;
  while (i) {
    (* GFS_MAP (i->data)->inverse_vector) (i->data, p, v, v);
    i = i->next;
  }
  sim->maps->items = g_slist_reverse (reverse);
}

/**
 * gfs_simulation_map_inverse_cell:
 * @sim: a #GfsSimulation.
 * @p: an array of #FttVector.
 *
 * Array @p should contain the (x,y) coordinates of the four corners
 * of a #FttCell. This function will inverse transform these
 * coordinates taking into account any specificity of the mapping
 * (e.g. periodicity).
 */
void gfs_simulation_map_inverse_cell (GfsSimulation * sim, FttVector p[4])
{
  g_return_if_fail (sim != NULL);
  g_return_if_fail (p != NULL);
  
  GSList * reverse = g_slist_reverse (sim->maps->items);
  GSList * i = reverse;
  while (i) {
    (* GFS_MAP (i->data)->inverse_cell) (i->data, p, p);
    i = i->next;
  }
  sim->maps->items = g_slist_reverse (reverse);

  FttComponent j, c;
  for (j = 0; j < 4; j++)
    for (c = 0; c < 3; c++)
      (&p[j].x)[c] *= sim->physical_params.L/(&GFS_DOMAIN (sim)->lambda.x)[c];
}

/**
 * gfs_dimensional_value:
 * @v: a #GfsVariable.
 * @val: a non-dimensional value of @v.
 *
 * Returns: the dimensional value of @val according to the units of @v.
 */
gdouble gfs_dimensional_value (GfsVariable * v, gdouble val)
{
  g_return_val_if_fail (v != NULL, 0.);

  gdouble L;
  if (val == GFS_NODATA || v->units == 0. || 
      (L = GFS_SIMULATION (v->domain)->physical_params.L) == 1.)
    return val;
  return val*pow (L, v->units);
}

/**
 * gfs_variable_is_dimensional:
 * @v: a #GfsVariable.
 *
 * Returns: %TRUE if @v has dimensions, %FALSE otherwise.
 */
gboolean gfs_variable_is_dimensional (GfsVariable * v)
{
  g_return_val_if_fail (v != NULL, FALSE);

  if (v->units == 0. || GFS_SIMULATION (v->domain)->physical_params.L == 1.)
    return FALSE;
  return TRUE;
}

/** \endobject{GfsSimulation} */

/**
 * The advection solver.
 * \beginobject{GfsAdvection}
 */

static void event_do_adapt (GfsEvent * event, GfsSimulation * sim)
{
  if (GFS_IS_ADAPT (event))
    gfs_event_do (event, sim);
}

static void event_do_not_adapt (GfsEvent * event, GfsSimulation * sim)
{
  if (!GFS_IS_ADAPT (event))
    gfs_event_do (event, sim);
}

static void advection_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  gboolean streamfunction = FALSE;

#if FTT_2D
  GSList * i = domain->variables;
  while (i && !streamfunction) {
    streamfunction = (GFS_IS_VARIABLE_STREAM_FUNCTION (i->data) != NULL);
    i = i->next;
  }
#endif /* 2D */

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  /* ignore default advection scheme (use per tracer parameters only) */
  sim->advection_params.scheme = GFS_NONE;
  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) event_do_adapt, sim);

    gfs_domain_cell_traverse (domain,
    			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
    			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) event_do_not_adapt, sim);

    if (!streamfunction) {
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity,
				gfs_domain_velocity (domain));
    }

    gfs_simulation_set_timestep (sim);

    gfs_advance_tracers (sim, sim->advection_params.dt);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
}

static void gfs_advection_class_init (GfsSimulationClass * klass)
{
  klass->run = advection_run;
}

GfsSimulationClass * gfs_advection_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_advection_info = {
      "GfsAdvection",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_advection_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_advection_info);
  }

  return klass;
}

/** \endobject{GfsAdvection} */

/**
 * The Poisson solver.
 * \beginobject{GfsPoisson}
 */

typedef struct {
  GfsVariable * divu, * div;
  gdouble ddiv;
  GtsRange vol;
} DivData;

static void rescale_div (FttCell * cell, DivData * p)
{
  gdouble size = ftt_cell_size (cell);
  gdouble a = size*size*gfs_domain_cell_fraction (p->div->domain, cell);
  GFS_VALUE (cell, p->div) = GFS_VALUE (cell, p->divu)*a;
  gts_range_add_value (&p->vol, a);
}

static void add_ddiv (FttCell * cell, DivData * p)
{
  gdouble size = ftt_cell_size (cell);
  GFS_VALUE (cell, p->div) += size*size*p->ddiv*gfs_domain_cell_fraction (p->div->domain, cell);
}

static void correct_div (GfsDomain * domain, GfsVariable * divu, GfsVariable * div,
			 gboolean dirichlet)
{
  DivData p = { divu, div };
  gts_range_init (&p.vol);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) rescale_div, &p);
  if (!dirichlet) {
    gts_range_update (&p.vol);
    
    p.ddiv = - gfs_domain_stats_variable (domain, div, 
					  FTT_TRAVERSE_LEAFS, -1, 
					  NULL, NULL).mean/p.vol.mean;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_ddiv, &p);
  }
}

static void has_dirichlet (FttCell * cell, GfsVariable * p)
{
  if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    p->centered = FALSE;
}

typedef struct {
  GfsVariable * lhs;
  gboolean dirichlet;
} CompatPar;

static void check_box_dirichlet (GfsBox * box, CompatPar * p)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->lhs);
      if (GFS_IS_BC_DIRICHLET (bc)) {
	p->dirichlet = TRUE;
	return;
      }	
    }
}

static void poisson_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * dia, * div, * res = NULL, * res1, * p;
  GfsMultilevelParams * par = &sim->approx_projection_params;
  GSList * i;

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);
  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  p = gfs_variable_from_name (domain->variables, "P");
  gfs_domain_surface_bc (domain, p);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) has_dirichlet, p);
  if (!p->centered) {
    guint nf = gfs_check_solid_fractions (domain);
    if (nf > 0)
      g_warning ("the solid surface cuts %d boundary cells,\n"
		 "this may cause errors for the Poisson solution\n", nf);
  }
  gboolean dirichlet = !p->centered;
  gfs_all_reduce (domain, dirichlet, MPI_INT, MPI_MAX);
  if (!dirichlet) {
    /* check whether any boundary has a Dirichlet condition on @p */
    CompatPar q = { p, FALSE };
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) check_box_dirichlet, &q);
    gfs_all_reduce (domain, q.dirichlet, MPI_INT, MPI_MAX);
    dirichlet = q.dirichlet;
  }

  div = gfs_temporary_variable (domain);
  res1 = res ? res : gfs_temporary_variable (domain);
  dia = gfs_temporary_variable (domain);

  while (sim->time.i < sim->time.iend && sim->time.t < sim->time.end) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);
    
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    gfs_domain_surface_bc (domain, p);
    correct_div (domain, gfs_variable_from_name (domain->variables, "Div"), div, dirichlet);
    gfs_poisson_coefficients (domain, sim->physical_params.alpha, FALSE, p->centered, TRUE);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dia);

    par->poisson_solve (domain, par, p, div, res1, dia, 1.);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  gts_object_destroy (GTS_OBJECT (dia));
  gts_object_destroy (GTS_OBJECT (div));
  if (!res)
    gts_object_destroy (GTS_OBJECT (res1));
}

static void poisson_class_init (GfsSimulationClass * klass)
{
  klass->run = poisson_run;
}

static void poisson_init (GfsDomain * domain)
{
  gfs_domain_add_variable (domain, "Div", "Right-hand-side of the Poisson equation");
  GFS_SIMULATION (domain)->time.iend = 1;
}

GfsSimulationClass * gfs_poisson_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_poisson_info = {
      "GfsPoisson",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) poisson_class_init,
      (GtsObjectInitFunc) poisson_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_poisson_info);
  }

  return klass;
}

/** \endobject{GfsPoisson} */

/**
 * The axisymmetric Euler solver.
 * \beginobject{GfsAxi}
 */

static void axi_read (GtsObject ** object, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_axi_class ())->parent_class->read) (object, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (*object);
  domain->refpos.y = 0.5;
  GfsVariable * w = gfs_variable_from_name (domain->variables, "W");
  if (w) {
    if (!GFS_IS_VARIABLE_TRACER (w)) {
      gts_file_error (fp, "W (the azimuthal velocity) must be a tracer");
      return;
    }
    w->component = FTT_Y; /* to include the appropriate diffusion metric terms */
  }
}

static gdouble axi_face_metric (const GfsDomain * domain, const FttCellFace * face)
{
  FttVector p;
  ftt_face_pos (face, &p);
  return p.y;
}

static gdouble axi_cell_metric (const GfsDomain * domain, const FttCell * cell)
{
  FttVector p;
  gfs_cell_cm (cell, &p);
  return p.y;
}

static void axi_solid_metric (const GfsDomain * domain, const FttCell * cell, FttVector * m)
{
  g_assert (GFS_IS_MIXED (cell));
  m->x = m->y = GFS_STATE (cell)->solid->ca.y;
}

static gdouble axi_scale_metric (const GfsDomain * domain, const FttCell * cell, FttComponent c)
{
  if (c < 2)
    return 1.;
  FttVector p;
  gfs_cell_cm (cell, &p);
  return p.y;
}

static gdouble axi_viscous_metric_implicit (const GfsDomain * domain, 
					    FttCell * cell, FttComponent c)
{
  if (c != FTT_Y)
    return 0.;
  FttVector p;
  gfs_cell_cm (cell, &p);
  return 1./(p.y*p.y);
}

static void axi_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = axi_read;
}

static void axi_init (GfsDomain * domain)
{
  domain->face_metric  = axi_face_metric;
  domain->cell_metric  = axi_cell_metric;
  domain->solid_metric = axi_solid_metric;
  domain->scale_metric = axi_scale_metric;
  domain->viscous_metric_implicit = axi_viscous_metric_implicit;
}

GfsSimulationClass * gfs_axi_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_axi_info = {
      "GfsAxi",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) axi_class_init,
      (GtsObjectInitFunc) axi_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_axi_info);
  }

  return klass;
}

/** \endobject{GfsAxi} */

/**
 * The axisymmetric advection solver.
 * \beginobject{GfsAdvectionAxi}
 */

GfsSimulationClass * gfs_advection_axi_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_advection_axi_info = {
      "GfsAdvectionAxi",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_advection_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_axi_class ()),
				  &gfs_advection_axi_info);
  }

  return klass;
}

/** \endobject{GfsAdvectionAxi} */
