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
/** \file
 * \brief GfsVariable and its descendants.
 */

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "variable.h"
#include "vof.h"
#include "source.h"
#include "init.h"

/**
 * Base class for scalar fields.
 * \beginobject{GfsVariable}
 */

static void variable_init_domain (GfsVariable * v, GfsDomain * domain)
{
  v->i = gfs_domain_alloc (domain);
  v->domain = domain;
  GTS_OBJECT (v)->reserved = domain;
}

static void gfs_variable_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;
  GfsVariable * v, * old;

  if (GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  domain = (*o)->reserved;
  if (gfs_derived_variable_from_name (domain->derived_variables, fp->token->str)) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  v = GFS_VARIABLE (*o);
  v->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  if ((old = gfs_variable_from_name (domain->variables, v->name))) {
    GSList * i;
    if ((i = g_slist_find (domain->variables_io, old)))
      i->data = v;
    domain->variables = g_slist_remove (domain->variables, old);
    gts_object_destroy (GTS_OBJECT (old));
  }
  variable_init_domain (v, domain);
  domain->variables = g_slist_append (domain->variables, v);
}

static void gfs_variable_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", GFS_VARIABLE (o)->name);
}

static void gfs_variable_destroy (GtsObject * object)
{
  GfsVariable * v = GFS_VARIABLE (object);

  g_free (v->name);
  g_free (v->description);
  if (v->sources)
    gts_object_destroy (GTS_OBJECT (v->sources));
  if (v->surface_bc)
    gts_object_destroy (GTS_OBJECT (v->surface_bc));
  if (v->default_bc)
    gts_object_destroy (GTS_OBJECT (v->default_bc));
  if (v->domain) {
    gfs_domain_free (v->domain, v->i);
    v->domain->variables = g_slist_remove (v->domain->variables, v);
  }

  if (GFS_IS_VARIABLE_TRACER (v)) {
    FttComponent c;
    GfsVariableTracer * t = GFS_IS_VARIABLE_TRACER (v);
    for (c = 0; c < FTT_DIMENSION; c++)
      if (t->advection.sink[c])
	gts_object_destroy (GTS_OBJECT (t->advection.sink[c]));
  }

  (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->destroy) (object);
}

static void gfs_variable_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_variable_read;
  klass->write = gfs_variable_write;
  klass->destroy = gfs_variable_destroy;
}

static void gfs_variable_init (GfsVariable * v)
{
  GFS_EVENT (v)->istep = 1;
  v->centered = FALSE;
  v->component = FTT_DIMENSION;
  v->fine_coarse = (GfsVariableFineCoarseFunc) gfs_get_from_below_intensive;
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
}

GfsVariableClass * gfs_variable_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_info = {
      "GfsVariable",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) gfs_variable_class_init,
      (GtsObjectInitFunc) gfs_variable_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()), &gfs_variable_info);
  }

  return klass;
}

/**
 * gfs_variable_new:
 * @klass: a #GfsVariableClass.
 * @domain: a #GfsDomain.
 * @name: the name of the variable or %NULL.
 * @description: the variable description or %NULL.
 *
 * Returns: a newly allocated #GfsVariable or %NULL if a variable
 * named @name already exists in @domain.
 */
GfsVariable * gfs_variable_new (GfsVariableClass * klass,
				GfsDomain * domain,
				const gchar * name,
				const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (klass != NULL, NULL);
  g_return_val_if_fail (domain != NULL, NULL);

  if (name &&
      (gfs_variable_from_name (domain->variables, name) ||
       gfs_derived_variable_from_name (domain->derived_variables, name)))
    return NULL;

  v = GFS_VARIABLE (gts_object_new (GTS_OBJECT_CLASS (klass)));
  if (name)
    v->name = g_strdup (name);
  if (description)
    v->description = g_strdup (description);
  variable_init_domain (v, domain);

  return v;
}

/**
 * gfs_variable_from_name:
 * @i: the list of available #GfsVariable.
 * @name: the name of the variable to find.
 *
 * Returns: the #GfsVariable @name or %NULL if this variable name does
 * not exist.  
 */
GfsVariable * gfs_variable_from_name (GSList * i,
				      const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (i && (!GFS_VARIABLE (i->data)->name || strcmp (name, GFS_VARIABLE (i->data)->name)))
    i = i->next;
  return i ? GFS_VARIABLE (i->data) : NULL;
}

/**
 * gfs_variables_from_list:
 * @i: the list of available #GfsVariable.
 * @list: a malloc'ed string containing comma separated variable names.
 * @error: where to return the variable name in case of error.
 *
 * Returns: a list of variables or %NULL in case of error, in which
 * case *@error points to the name of the unknown variable.  
 */
GSList * gfs_variables_from_list (GSList * i,
				  gchar * list,
				  gchar ** error)
{
  gchar * s;
  GSList * var = NULL;

  g_return_val_if_fail (i != NULL, NULL);
  g_return_val_if_fail (error != NULL, NULL);

  s = strtok (list, ",");
  while (s) {
    GfsVariable * v = gfs_variable_from_name (i, s);

    if (v == NULL) {
      *error = s;
      g_slist_free (var);
      return NULL;
    }
    var = g_slist_append (var, v);
    s = strtok (NULL, ",");
  }
  return var;
}

/**
 * gfs_variables_swap:
 * @v1: a #GfsVariable.
 * @v2: a #GfsVariable.
 *
 * Swaps the values of @v1 and @v2, belonging to the same #GfsDomain.
 */
void gfs_variables_swap (GfsVariable * v1, GfsVariable * v2)
{
  guint i;

  g_return_if_fail (v1 != NULL);
  g_return_if_fail (v2 != NULL);
  g_return_if_fail (v1->domain == v2->domain);

  i = v1->i; v1->i = v2->i; v2->i = i;
}

/**
 * gfs_variable_set_vector:
 * @v: the components of the vector.
 * @n: the vector dimension.
 *
 * Sets @v[0],...,@v[n-1] as components of a vector quantity.
 */
void gfs_variable_set_vector (GfsVariable ** v, guint n)
{
  guint i, j;

  g_return_if_fail (v != NULL);
  g_return_if_fail (n > 1 && n <= FTT_DIMENSION);

  for (i = 0; i < n; i++) {
    g_return_if_fail (v[i] != NULL);
    v[i]->component = i;
    for (j = 0; j < n; j++)
      v[i]->vector[j] = v[j];
  }
  
  /* for rotated boundary conditions i.e. gfs_boundary_periodic_rotate() */
  v[0]->orientation =  1.;
  v[1]->orientation = -1.;
}

/**
 * gfs_variable_set_tensor:
 * @t: the components of the 2x2 tensor.
 *
 * Sets @t[0][0],...,@t[1][1] as components of a tensor.
 */
void gfs_variable_set_tensor (GfsVariable * t[2][2])
{
  g_return_if_fail (t != NULL);

  /* for rotated boundary conditions i.e. gfs_boundary_periodic_rotate() */
  t[0][0]->component = 0; t[0][0]->vector[0] = t[0][0]; t[0][0]->vector[1] = t[1][1];
  t[0][0]->orientation = + 1.; t[0][0]->even = TRUE;
  t[1][1]->component = 1; t[1][1]->vector[0] = t[0][0]; t[1][1]->vector[1] = t[1][1];
  t[1][1]->orientation = + 1.; t[1][1]->even = TRUE;
  t[0][1]->component = 0; t[0][1]->vector[0] = t[0][1]; t[0][1]->vector[1] = t[1][0];
  t[0][1]->orientation = - 1.; t[0][1]->even = TRUE;
  t[1][0]->component = 1; t[1][0]->vector[0] = t[0][1]; t[1][0]->vector[1] = t[1][0];
  t[1][0]->orientation = - 1.; t[1][0]->even = TRUE;
}

/**
 * gfs_variable_clone:
 * @v: a #GfsVariable.
 * @name: a name.
 *
 * Returns: a new #GfsVariable called @name and clone of @v.
 */
GfsVariable * gfs_variable_clone (GfsVariable * v, gchar * name)
{
  g_return_val_if_fail (v != NULL, NULL);
  g_return_val_if_fail (name != NULL, NULL);

  char * buf;
  size_t len;
  FILE * f = open_memstream (&buf, &len);
  if (f == NULL)
    g_error ("gfs_variable_clone(): could not open_memstream:\n%s", strerror (errno));
  gchar * s = v->name;
  v->name = name;
  GtsObject * o = GTS_OBJECT (v);
  (* o->klass->write) (o, f);
  fclose (f);
  v->name = s;
  GtsFile * fp = gts_file_new_from_buffer (buf, len);
  GtsObject * clone = gts_object_new (o->klass);
  gfs_object_simulation_set (clone, gfs_object_simulation (o));
  (* o->klass->read) (&clone, fp);
  if (fp->type == GTS_ERROR)
    g_error ("gfs_variable_clone:\n%d:%d:%s", fp->line, fp->pos, fp->error);
  gts_file_destroy (fp);
  free (buf);
  GFS_VARIABLE (clone)->units = v->units;
  GFS_VARIABLE (clone)->fine_coarse = v->fine_coarse;
  GFS_VARIABLE (clone)->coarse_fine = v->coarse_fine;
  return GFS_VARIABLE (clone);
}

/** \endobject{GfsVariable} */

/**
 * Boolean value consistent with adaptive refinement.
 * \beginobject{GfsVariableBoolean}
 */

static void boolean_fine_coarse (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  gint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && GFS_VALUE (child.c[i], v) < 0.) {
      GFS_VALUE (parent, v) = -1.;
      return;
    }

  gfs_get_from_below_intensive (parent, v);
}

static void boolean_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  gdouble val = GFS_VALUE (parent, v);
  gint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      GFS_VALUE (child.c[i], v) = val;
}

static void variable_boolean_init (GfsVariable * v)
{
  v->fine_coarse = boolean_fine_coarse;
  v->coarse_fine = boolean_coarse_fine;
  v->description = g_strdup ("Boolean");
}

GfsVariableClass * gfs_variable_boolean_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_boolean_info = {
      "GfsVariableBoolean",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) variable_boolean_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_boolean_info);
  }

  return klass;
}

/** \endobject{GfsVariableBoolean} */

/**
 * Advected scalar fields.
 * \beginobject{GfsVariableTracer}
 */

static void variable_tracer_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == '{')
    gfs_advection_params_read (&GFS_VARIABLE_TRACER (*o)->advection, fp);

  if (fp->type != GTS_ERROR && fp->type == '{')
    g_warning ("%d:%d: specifying diffusion parameters is not done here anymore!",
	       fp->line, fp->pos);
}

static void variable_tracer_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_class ())->parent_class->write) (o, fp);

  fputc (' ', fp);
  gfs_advection_params_write (&GFS_VARIABLE_TRACER (o)->advection, fp);
}

static void variable_tracer_class_init (GtsObjectClass * klass)
{
  klass->read = variable_tracer_read;
  klass->write = variable_tracer_write;
}

static void variable_tracer_init (GfsVariableTracer * v)
{
  gfs_advection_params_init (&v->advection);
  v->advection.gradient = gfs_center_van_leer_gradient;
  v->advection.flux = gfs_face_advection_flux;
  v->advection.v = GFS_VARIABLE (v);
  v->advection.fv = NULL;
  GFS_VARIABLE (v)->description = g_strdup ("Tracer");
}

GfsVariableClass * gfs_variable_tracer_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_tracer_info = {
      "GfsVariableTracer",
      sizeof (GfsVariableTracer),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_tracer_class_init,
      (GtsObjectInitFunc) variable_tracer_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_tracer_info);
  }

  return klass;
}

/** \endobject{GfsVariableTracer} */

/**
 * \beginobject{GfsVariableResidual}
 */

static void scale_residual (FttCell * cell, GfsVariable * res)
{
  gdouble size = ftt_cell_size (cell);
  gdouble dt = GFS_SIMULATION (res->domain)->advection_params.dt;
  GFS_VALUE (cell, res) *= dt*dt/(size*size);
}

static gboolean variable_residual_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_class ())->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) scale_residual, event);
    return TRUE;
  }
  return FALSE;
}

static void variable_residual_class_init (GfsEventClass * klass)
{
  klass->event = variable_residual_event;
}

static void variable_residual_init (GfsVariable * v)
{
  v->description = g_strdup ("Residual of the Poisson equation");
}

GfsVariableClass * gfs_variable_residual_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_residual_info = {
      "GfsVariableResidual",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_residual_class_init,
      (GtsObjectInitFunc) variable_residual_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_residual_info);
  }

  return klass;
}

/** \endobject{GfsVariableResidual} */

/**
 * Spatial filtering.
 * \beginobject{GfsVariableFiltered}
 */

static void variable_filtered_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableFiltered * v = GFS_VARIABLE_FILTERED (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting a number (niter)");
    return;
  }
  v->niter = atoi (fp->token->str);
  if (v->niter == 0) {
    gts_file_error (fp, "niter must be strictly positive");
    return;
  }
  gts_file_next_token (fp);  

  if (GFS_VARIABLE (v)->description)
    g_free (GFS_VARIABLE (v)->description);
  GFS_VARIABLE (v)->description = g_strjoin (" ", "Variable", v->v->name, "filtered", NULL);

  GFS_VARIABLE (v)->units = v->v->units;
}

static void variable_filtered_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %d", GFS_VARIABLE_FILTERED (o)->v->name, GFS_VARIABLE_FILTERED (o)->niter);
}

static void variable_filtered_event_half (GfsEvent * event, GfsSimulation * sim)
{
  guint n = GFS_VARIABLE_FILTERED (event)->niter;
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * v = GFS_VARIABLE (event);

  gfs_domain_filter (domain, GFS_VARIABLE_FILTERED (event)->v, v);
  while (--n)
    gfs_domain_filter (domain, v, NULL);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);
  gfs_domain_bc (domain, FTT_TRAVERSE_NON_LEAFS, -1, v);
}

static gboolean variable_filtered_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class)->event)
      (event, sim)) {
    variable_filtered_event_half (event, sim);
    return TRUE;
  }
  return FALSE;
}

static void variable_filtered_class_init (GtsObjectClass * klass)
{
  klass->read = variable_filtered_read;
  klass->write = variable_filtered_write;
  GFS_EVENT_CLASS (klass)->event = variable_filtered_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_filtered_event_half;
}

static void variable_filtered_init (GfsEvent * v)
{
  /* the variable/event may need to be initialised at the start */
  v->start = -1;
}

GfsVariableClass * gfs_variable_filtered_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_filtered_info = {
      "GfsVariableFiltered",
      sizeof (GfsVariableFiltered),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_filtered_class_init,
      (GtsObjectInitFunc) variable_filtered_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_filtered_info);
  }

  return klass;
}

/** \endobject{GfsVariableFiltered} */

/**
 * \beginobject{GfsVariableDiagonal}
 */

static void unity (FttCell * cell, GfsVariable * v)
{
  GFS_VALUE (cell, v) = 1.;
}

static void variable_diagonal (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * tmp = data[1];
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  GFS_VALUE (cell, tmp) = G_MAXDOUBLE;
  g.a = g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, tmp->i, -1);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a != 0.)
    GFS_VALUE (cell, v) = g.b/g.a;
  else
    GFS_VALUE (cell, v) = G_MAXDOUBLE;
  GFS_VALUE (cell, tmp) = 1.;
}

static gboolean variable_diagonal_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_diagonal_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * tmp = gfs_temporary_variable (domain);
    gpointer data[2];
    data[0] = event;
    data[1] = tmp;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) unity, tmp);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, tmp);
    gfs_poisson_coefficients (domain, sim->physical_params.alpha, TRUE, TRUE, TRUE);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) variable_diagonal, data);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE (event));
    gts_object_destroy (GTS_OBJECT (tmp));
    return TRUE;
  }
  return FALSE;
}

static void variable_diagonal_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_diagonal_event;
}

GfsVariableClass * gfs_variable_diagonal_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_diagonal_info = {
      "GfsVariableDiagonal",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_diagonal_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_diagonal_info);
  }

  return klass;
}

/** \endobject{GfsVariableDiagonal} */

/**
 * Optimising function evaluations .
 * \beginobject{GfsVariableFunction}
 */

static void variable_function_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_FUNCTION (o)->f));

  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->destroy) (o);
}

static void variable_function_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariableFunction * v = GFS_VARIABLE_FUNCTION (*o);
  gfs_function_read (v->f, GFS_DOMAIN (gfs_object_simulation (*o)), fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_set_units (v->f, GFS_VARIABLE (*o)->units);
}

static void variable_function_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->write) (o, fp);

  gfs_function_write (GFS_VARIABLE_FUNCTION (o)->f, fp);
}

static void variable_function_class_init (GtsObjectClass * klass)
{
  klass->destroy = variable_function_destroy;
  klass->read = variable_function_read;
  klass->write = variable_function_write;
}

static void variable_function_coarse_fine (FttCell * parent, GfsVariable * v)
{
  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], v) = gfs_function_value (f, child.c[n]);
}

static void variable_function_init (GfsVariableFunction * v)
{
  GFS_VARIABLE (v)->coarse_fine = variable_function_coarse_fine;
  v->f = gfs_function_new (gfs_function_class (), 0.);
}

GfsVariableClass * gfs_variable_function_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_function_info = {
      "GfsVariableFunction",
      sizeof (GfsVariableFunction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_function_class_init,
      (GtsObjectInitFunc) variable_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_function_info);
  }

  return klass;
}

/** \endobject{GfsVariableFunction} */

#if FTT_2D

/**
 * Maintaining a velocity field defined by a stream function.
 * \beginobject{GfsVariableStreamFunction}
 */

static gdouble face_metric (FttCell * cell, FttDirection d, GfsDomain * domain)
{
  if (domain->face_metric) {
    FttCellFace f;
    f.cell = cell;
    f.d = d;
    return (* domain->face_metric) (domain, &f);
  }
  else
    return 1.;
}

static gdouble face_metric_inverse (FttCell * cell, FttDirection d, GfsDomain * domain)
{
  gdouble fm = face_metric (cell, d, domain);
  /* for degenerate metric e.g. lon-lat */
  return fm > 1e-6 ? 1./fm : 0.;
}

static void init_mac_from_stream_function (FttCell * cell,
					   gdouble psi0, gdouble psi1, gdouble psi2, gdouble psi3,
					   gdouble h,
					   GfsDomain * domain,
					   GfsVariable ** u)
{
  GFS_STATE (cell)->f[0].un = (psi2 - psi1)*face_metric_inverse (cell, 0, domain)/h;
  GFS_STATE (cell)->f[1].un = (psi3 - psi0)*face_metric_inverse (cell, 1, domain)/h;
  GFS_STATE (cell)->f[2].un = (psi3 - psi2)*face_metric_inverse (cell, 2, domain)/h;
  GFS_STATE (cell)->f[3].un = (psi0 - psi1)*face_metric_inverse (cell, 3, domain)/h;

  GFS_VALUE (cell, u[0]) = (GFS_STATE (cell)->f[0].un + GFS_STATE (cell)->f[1].un)/2.;
  GFS_VALUE (cell, u[1]) = (GFS_STATE (cell)->f[2].un + GFS_STATE (cell)->f[3].un)/2.;
}

static void variable_stream_function_coarse_fine (FttCell * parent, GfsVariable * v)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  FttVector o, p;
  ftt_cell_pos (parent, &o);
  gdouble h = ftt_cell_size (parent)/2.;
  p.x = o.x - h; p.y = o.y - h;
  gdouble psi0 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y - h;
  gdouble psi1 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y + h;
  gdouble psi2 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y + h;
  gdouble psi3 = gfs_function_spatial_value (f, &p);
  p.x = o.x; p.y = o.y - h;
  gdouble psi4 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y;
  gdouble psi5 = gfs_function_spatial_value (f, &p);
  p.x = o.x; p.y = o.y + h;
  gdouble psi6 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y;
  gdouble psi7 = gfs_function_spatial_value (f, &p);
  gdouble psi8 = gfs_function_spatial_value (f, &o);
  GfsVariable ** u = gfs_domain_velocity (v->domain);
  init_mac_from_stream_function (child.c[0], psi7, psi8, psi6, psi3, h, v->domain, u);
  init_mac_from_stream_function (child.c[1], psi8, psi5, psi2, psi6, h, v->domain, u);
  init_mac_from_stream_function (child.c[2], psi0, psi4, psi8, psi7, h, v->domain, u);
  init_mac_from_stream_function (child.c[3], psi4, psi1, psi5, psi8, h, v->domain, u);
  guint n;
  for (n = 0; n < FTT_CELLS; n++) {
    ftt_cell_pos (child.c[n], &o);
    GFS_VALUE (child.c[n], v) = gfs_function_spatial_value (f, &o);
  }
}

static void variable_stream_function_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariableFunction * v = GFS_VARIABLE_FUNCTION (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (fp->type == GTS_STRING && gfs_variable_from_name (domain->variables, fp->token->str)) {
    /* variable */
    gfs_function_read (v->f, domain, fp);
    g_assert (gfs_function_get_variable (v->f));
    if (fp->type == GTS_ERROR)
      return;
    /* fixme: what about the equivalent coarse_fine() method? */
  }
  else { /* spatial function */
    gts_object_destroy (GTS_OBJECT (v->f));
    v->f = gfs_function_new (gfs_function_spatial_class (), 0.);
    gfs_function_read (v->f, domain, fp);
    if (fp->type == GTS_ERROR)
      return;
    gfs_function_set_units (v->f, GFS_VARIABLE (*o)->units);
    GFS_VARIABLE (v)->coarse_fine = variable_stream_function_coarse_fine;
  }
}

static void variable_stream_function_fine_coarse (FttCell * cell, GfsVariable * v)
{
  FttCellChildren child;
  ftt_cell_children (cell, &child);
  GFS_STATE (cell)->f[0].un = 
    face_metric_inverse (cell, 0, v->domain)/2.*
    (face_metric (child.c[1], 0, v->domain)*GFS_STATE (child.c[1])->f[0].un +
     face_metric (child.c[3], 0, v->domain)*GFS_STATE (child.c[3])->f[0].un);
  GFS_STATE (cell)->f[1].un = 
    face_metric_inverse (cell, 1, v->domain)/2.*
    (face_metric (child.c[0], 1, v->domain)*GFS_STATE (child.c[0])->f[1].un +
     face_metric (child.c[2], 1, v->domain)*GFS_STATE (child.c[2])->f[1].un);
  GFS_STATE (cell)->f[2].un = 
    face_metric_inverse (cell, 2, v->domain)/2.*
    (face_metric (child.c[0], 2, v->domain)*GFS_STATE (child.c[0])->f[2].un +
     face_metric (child.c[1], 2, v->domain)*GFS_STATE (child.c[1])->f[2].un);
  GFS_STATE (cell)->f[3].un = 
    face_metric_inverse (cell, 3, v->domain)/2.*
    (face_metric (child.c[3], 3, v->domain)*GFS_STATE (child.c[3])->f[3].un +
     face_metric (child.c[2], 3, v->domain)*GFS_STATE (child.c[2])->f[3].un);
  GFS_VALUE (cell, v) = (GFS_VALUE (child.c[0], v) + GFS_VALUE (child.c[1], v) + 
			 GFS_VALUE (child.c[2], v) + GFS_VALUE (child.c[3], v))/4.;
}

static void init_streamfunction (FttCell * cell, GfsVariable * v)
{
  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttVector o, p;
  ftt_cell_pos (cell, &o);
  gdouble h = ftt_cell_size (cell)/2.;
  p.x = o.x - h; p.y = o.y - h;
  gdouble psi0 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y - h;
  gdouble psi1 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y + h;
  gdouble psi2 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y + h;
  gdouble psi3 = gfs_function_spatial_value (f, &p);
  init_mac_from_stream_function (cell, psi0, psi1, psi2, psi3, 2.*h, 
				 v->domain, gfs_domain_velocity (v->domain));
}

static void init_streamfunction_from_variable (FttCell * cell, GfsVariable * v)
{
  FttDirection d[FTT_DIMENSION];
  d[0] = FTT_LEFT; d[1] = FTT_BOTTOM;
  gdouble psi0 = gfs_cell_corner_value (cell, d, v, -1);
  d[0] = FTT_RIGHT; d[1] = FTT_BOTTOM;
  gdouble psi1 = gfs_cell_corner_value (cell, d, v, -1);
  d[0] = FTT_RIGHT; d[1] = FTT_TOP;
  gdouble psi2 = gfs_cell_corner_value (cell, d, v, -1);
  d[0] = FTT_LEFT; d[1] = FTT_TOP;
  gdouble psi3 = gfs_cell_corner_value (cell, d, v, -1);
  init_mac_from_stream_function (cell, psi0, psi1, psi2, psi3, ftt_cell_size (cell),
				 v->domain, gfs_domain_velocity (v->domain));
}

static gboolean variable_stream_function_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_function_class ())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * v = gfs_function_get_variable (GFS_VARIABLE_FUNCTION (event)->f);
    if (v)
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) init_streamfunction_from_variable, 
				  v);
    else
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) init_streamfunction, event);
    GfsVariable ** u = gfs_domain_velocity (domain);
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, u[c]);
    return TRUE;
  }
  return FALSE;
}

static void variable_stream_function_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = variable_stream_function_read;
  klass->event = variable_stream_function_event;
}

static void variable_stream_function_init (GfsVariable * v)
{
  v->units = 2.;
  v->fine_coarse = variable_stream_function_fine_coarse;
  GFS_EVENT (v)->istep = G_MAXINT/2;
}

GfsVariableClass * gfs_variable_stream_function_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_stream_function_info = {
      "GfsVariableStreamFunction",
      sizeof (GfsVariableFunction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_stream_function_class_init,
      (GtsObjectInitFunc) variable_stream_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_function_class ()), 
				  &gfs_variable_stream_function_info);
  }

  return klass;
}

/** \endobject{GfsVariableStreamFunction} */

#endif /* FTT_2D */

/**
 * A variable, solution of a Poisson equation.
 * \beginobject{GfsVariablePoisson}
 */

static void variable_poisson_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_poisson_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == '{')
    gfs_multilevel_params_read (&GFS_VARIABLE_POISSON (*o)->par, fp);
}

static void variable_poisson_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_poisson_class ())->parent_class->write) (o, fp);

  fputc (' ', fp);
  gfs_multilevel_params_write (&GFS_VARIABLE_POISSON (o)->par, fp);
}

static void has_dirichlet (FttCell * cell, GfsVariable * p)
{
  if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    p->centered = FALSE;
}

typedef struct {
  GfsFunction * f;
  GfsVariable * div;
} DivData;

static void rescale_div (FttCell * cell, DivData * p)
{
  gdouble size = ftt_cell_size (cell);
  gdouble a = size*size*gfs_domain_cell_fraction (p->div->domain, cell);
  GFS_VALUE (cell, p->div) = gfs_function_value (p->f, cell)*a;
}

static gboolean variable_poisson_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_class ())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * dia = gfs_temporary_variable (domain);
    GfsVariable * div = gfs_temporary_variable (domain);
    GfsVariable * res = gfs_temporary_variable (domain);
    GfsVariable * v = GFS_VARIABLE (event);
    GfsVariablePoisson * pv = GFS_VARIABLE_POISSON (v);

    gfs_domain_surface_bc (domain, v);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) has_dirichlet, v);
    DivData p = { GFS_VARIABLE_FUNCTION (event)->f, div };
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) rescale_div, &p);
    gfs_poisson_coefficients (domain, NULL, FALSE, TRUE, TRUE);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dia);
    pv->par.poisson_solve (domain, &pv->par, v, div, res, dia, 1.);
    if (pv->par.residual.infty > pv->par.tolerance) {
      g_warning ("VariablePoisson %s: max residual %g > %g",
		 v->name, pv->par.residual.infty, pv->par.tolerance);
      gfs_multilevel_params_stats_write (&pv->par, stderr);
    }

    gts_object_destroy (GTS_OBJECT (dia));
    gts_object_destroy (GTS_OBJECT (div));
    gts_object_destroy (GTS_OBJECT (res));
    return TRUE;
  }
  return FALSE;  
}

static void variable_poisson_class_init (GtsObjectClass * klass)
{
  klass->read = variable_poisson_read;
  klass->write = variable_poisson_write;
  GFS_EVENT_CLASS (klass)->event = variable_poisson_event;
}

static void variable_poisson_init (GfsVariable * v)
{
  v->centered = TRUE;
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
  gfs_multilevel_params_init (&GFS_VARIABLE_POISSON (v)->par);
}

GfsVariableClass * gfs_variable_poisson_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_poisson_info = {
      "GfsVariablePoisson",
      sizeof (GfsVariablePoisson),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_poisson_class_init,
      (GtsObjectInitFunc) variable_poisson_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_function_class ()), 
				  &gfs_variable_poisson_info);
  }

  return klass;
}

/** \endobject{GfsVariablePoisson} */

/**
 * The Laplacian of a function.
 * \beginobject{GfsVariableLaplacian}
 */

typedef struct {
  GfsVariable * lap;
  GfsFunction * f;
  GfsVariable * u;
} LapData;

static void evaluate_function (FttCell * cell, LapData * p)
{
  GFS_VALUE (cell, p->u) = gfs_function_value (p->f, cell);
}

static void laplacian (FttCell * cell, LapData * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, p->u->i, -1);
      g.a += ng.a;
      g.b += ng.b;
    }
  }

  gdouble h = ftt_cell_size (cell);
  GFS_VALUE (cell, p->lap) = (g.b - g.a*GFS_VALUE (cell, p->u))/
    (h*h*gfs_domain_cell_fraction (p->lap->domain, cell));
}

static gboolean variable_laplacian_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_class ())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    LapData p = { GFS_VARIABLE (event), GFS_VARIABLE_FUNCTION (event)->f };
    p.u = gfs_function_get_variable (p.f);
    if (p.u == NULL) {
      p.u = gfs_temporary_variable (domain);
      gfs_function_set_units (p.f, 0.);
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) evaluate_function, &p);
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.u);
    }
    gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) laplacian, &p);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.lap);
    if (!gfs_function_get_variable (p.f))
      gts_object_destroy (GTS_OBJECT (p.u));
    return TRUE;
  }
  return FALSE;  
}

static void variable_laplacian_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_laplacian_event;
}

static void variable_laplacian_init (GfsVariable * v)
{
  v->units = -2.;
  v->centered = TRUE;
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
}

GfsVariableClass * gfs_variable_laplacian_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableLaplacian",
      sizeof (GfsVariableFunction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_laplacian_class_init,
      (GtsObjectInitFunc) variable_laplacian_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_function_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableLaplacian} */

/**
 * How old is this cell?.
 * \beginobject{GfsVariableAge}
 */

static void increment_age (FttCell * cell, GfsVariable * v)
{
  GFS_VALUE (cell, v) += 1.;
}

static gboolean variable_age_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_age_class ())->parent_class)->event)
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) increment_age, event);
    return TRUE;
  }
  return FALSE;
}

static void variable_age_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_age_event;
}

static void none (FttCell * parent, GfsVariable * v)
{
}

static void variable_age_init (GfsVariable * v)
{
  v->fine_coarse = none;
  v->coarse_fine = none;
}

GfsVariableClass * gfs_variable_age_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_age_info = {
      "GfsVariableAge",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_age_class_init,
      (GtsObjectInitFunc) variable_age_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_age_info);
  }

  return klass;
}

/** \endobject{GfsVariableAge} */

/**
 * Averaging along a coordinate direction.
 * \beginobject{GfsVariableAverage}
 */

static void variable_average_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_average_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (component)");
    return;
  }
  gchar s[][2] = {"x", "y", "z"};
  GfsVariableAverage * v = GFS_VARIABLE_AVERAGE (*o);
  for (v->c = 0; v->c < FTT_DIMENSION; v->c++)
    if (!strcmp (fp->token->str, s[v->c]))
      break;
  if (v->c == FTT_DIMENSION) {
    gts_file_error (fp, "unknown component '%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void variable_average_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_average_class ())->parent_class->write) (o, fp);

  gchar s[][2] = {"x", "y", "z"};
  fprintf (fp, " %s", s[GFS_VARIABLE_AVERAGE (o)->c]);
}

/* see also hydrostatic_pressure() */
static void average (FttCell * cell, GfsVariable * v)
{
  FttDirection d = 2*GFS_VARIABLE_AVERAGE (v)->c + 1;
  gdouble avg = 0., vol = 0.;
  GtsFifo * fifo = gts_fifo_new (), * column = gts_fifo_new ();

  gts_fifo_push (fifo, cell);
  while ((cell = gts_fifo_pop (fifo))) {
    double w = gfs_cell_volume (cell, v->domain);
    vol += w;
    avg += w*gfs_function_value (GFS_VARIABLE_FUNCTION (v)->f, cell);
    gts_fifo_push (column, cell);

    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    if (neighbor) {
      if (FTT_CELL_IS_LEAF (neighbor)) {
	if (ftt_cell_level (neighbor) == ftt_cell_level (cell))
	  /* neighbor at same level */
	  gts_fifo_push (fifo, neighbor);
	else {
	  /* coarser neighbour */
#if 1
	  g_assert_not_implemented ();
#else
	  if (gts_fifo_top (fifo) == NULL) { /* only consider the last fine cell */
	    FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	    double dp = GFS_STATE (neighbor)->f[od].un*ftt_cell_size (neighbor)/
	      GFS_STATE (neighbor)->f[od].v;
	    double p = 0.;
	    FttCellChildren child;
	    int i, n = ftt_cell_children_direction (ftt_cell_parent (cell), d, &child);
	    for (i = 0; i < n; i++)
	      p += GFS_VALUE (child.c[i], v);
	    GFS_VALUE (neighbor, v) = p/n - 3.*dp/4.;
	    gts_fifo_push (fifo, neighbor);
	  }
#endif
	}
      }
      else {
	/* finer neighbor */
#if 1
	  g_assert_not_implemented ();
#else
	FttCellChildren child;
	int i, n = ftt_cell_children_direction (neighbor, FTT_OPPOSITE_DIRECTION (d), &child);
	double dp = GFS_STATE (cell)->f[d].un*ftt_cell_size (cell)/GFS_STATE (cell)->f[d].v;
	double p = GFS_VALUE (cell, v) - 3.*dp/4.;
	for (i = 0; i < n; i++) {
	  GFS_VALUE (child.c[i], v) = p;
	  gts_fifo_push (fifo, child.c[i]);
	}
#endif
      }
    }
  }
  gts_fifo_destroy (fifo);

  if (vol > 0.)
    avg /= vol;

  while ((cell = gts_fifo_pop (column)))
    GFS_VALUE (cell, v) = avg;
  gts_fifo_destroy (column);
}

static gboolean variable_average_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class)->event)
      (event, sim)) {
    GfsVariable * v = GFS_VARIABLE (event);
    GfsVariableAverage * av = GFS_VARIABLE_AVERAGE (v);
    gfs_domain_cell_traverse_boundary (v->domain, 2*av->c,
				       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				       (FttCellTraverseFunc) average, av);
    gfs_domain_bc (v->domain, FTT_TRAVERSE_LEAFS, -1, v);
    return TRUE;
  }
  return FALSE;
}

static void variable_average_class_init (GtsObjectClass * klass)
{
  klass->read = variable_average_read;
  klass->write = variable_average_write;
  GFS_EVENT_CLASS (klass)->event = variable_average_event;
}

static void variable_average_init (GfsVariable * v)
{
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
}

GfsVariableClass * gfs_variable_average_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableAverage",
      sizeof (GfsVariableAverage),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_average_class_init,
      (GtsObjectInitFunc) variable_average_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_function_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableAverage} */

/**
 * The hydrostatic pressure
 * \beginobject{GfsHydrostaticPressure}
 */

static FttComponent hydrostatic_component (GfsDomain * domain)
{
  GfsVariable ** u = gfs_domain_velocity (domain);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (u[c]->sources) {
      GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;      
      while (i) {
	GfsSourceGeneric * s = i->data;
	if (s->face_value)
	  return c;
	i = i->next;
      }
    }
  return FTT_DIMENSION;
}

static void hydrostatic_pressure_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_hydrostatic_pressure_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsHydrostaticPressure * ph = GFS_HYDROSTATIC_PRESSURE (*o);
  ph->c = hydrostatic_component (GFS_VARIABLE (*o)->domain);
  if (ph->c == FTT_DIMENSION) {
    gts_file_error (fp, "could not find any velocity sources");
    return;
  }

  GFS_VARIABLE (*o)->units = 2.;
}

static gboolean hydrostatic_pressure_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_class ())->event) (event, sim)) {
    gfs_hydrostatic_pressure_update (GFS_HYDROSTATIC_PRESSURE (event), sim->physical_params.alpha);
    return TRUE;
  }
  return FALSE;
}

static void hydrostatic_pressure_class_init (GtsObjectClass * klass)
{
  klass->read = hydrostatic_pressure_read;
  GFS_EVENT_CLASS (klass)->event = hydrostatic_pressure_event;
}

static void hydrostatic_pressure_init (GfsVariable * v)
{
  GFS_EVENT (v)->start = -1;
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
}

GfsVariableClass * gfs_hydrostatic_pressure_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsHydrostaticPressure",
      sizeof (GfsHydrostaticPressure),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) hydrostatic_pressure_class_init,
      (GtsObjectInitFunc) hydrostatic_pressure_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), &info);
  }

  return klass;
}

static void hydrostatic_pressure (FttCell * cell, GfsVariable * v)
{
  FttDirection d = 2*GFS_HYDROSTATIC_PRESSURE (v)->c + 1;
  GtsFifo * fifo = gts_fifo_new ();

  GFS_VALUE (cell, v) = 0.;
  gts_fifo_push (fifo, cell);

  while ((cell = gts_fifo_pop (fifo))) {
    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    if (neighbor) {
      if (FTT_CELL_IS_LEAF (neighbor)) {
	if (ftt_cell_level (neighbor) == ftt_cell_level (cell)) {
	  /* neighbor at same level */
	  double dp = GFS_STATE (cell)->f[d].un*ftt_cell_size (cell)/GFS_STATE (cell)->f[d].v;
	  GFS_VALUE (neighbor, v) = GFS_VALUE (cell, v) - dp;
	  gts_fifo_push (fifo, neighbor);
	}
	else {
	  /* coarser neighbour */
	  if (gts_fifo_top (fifo) == NULL) { /* only consider the last fine cell */
	    FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	    double dp = GFS_STATE (neighbor)->f[od].un*ftt_cell_size (neighbor)/
	      GFS_STATE (neighbor)->f[od].v;
	    double p = 0.;
	    FttCellChildren child;
	    int i, n = ftt_cell_children_direction (ftt_cell_parent (cell), d, &child);
	    for (i = 0; i < n; i++)
	      p += GFS_VALUE (child.c[i], v);
	    GFS_VALUE (neighbor, v) = p/n - 3.*dp/4.;
	    gts_fifo_push (fifo, neighbor);
	  }
	}
      }
      else {
	/* finer neighbor */
	FttCellChildren child;
	int i, n = ftt_cell_children_direction (neighbor, FTT_OPPOSITE_DIRECTION (d), &child);
	double dp = GFS_STATE (cell)->f[d].un*ftt_cell_size (cell)/GFS_STATE (cell)->f[d].v;
	double p = GFS_VALUE (cell, v) - 3.*dp/4.;
	for (i = 0; i < n; i++) {
	  GFS_VALUE (child.c[i], v) = p;
	  gts_fifo_push (fifo, child.c[i]);
	}
      }
    }
  }

  gts_fifo_destroy (fifo);
}

/**
 * gfs_hydrostatic_pressure_update:
 * @p: a #GfsHydrostaticPressure.
 * @alpha: the reciprocal of density or %NULL.
 *
 * Updates the hydrostatic pressure field @p. Note that face
 * velocities are also reset.
 */
void gfs_hydrostatic_pressure_update (GfsHydrostaticPressure * p, GfsFunction * alpha)
{
  g_return_if_fail (p != NULL);

  GfsVariable * v = GFS_VARIABLE (p);
  gfs_domain_face_traverse (v->domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_velocity_face_sources (v->domain, gfs_domain_velocity (v->domain), 1., NULL, NULL);
  gfs_poisson_coefficients (v->domain, alpha, TRUE, TRUE, TRUE);
  gfs_domain_cell_traverse_boundary (v->domain, 2*p->c,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				     (FttCellTraverseFunc) hydrostatic_pressure, p);
  gfs_domain_bc (v->domain, FTT_TRAVERSE_LEAFS, -1, v);
}

/** \endobject{GfsHydrostaticPressure} */

/**
 * Spatially-constant but time-dependent variables.
 * \beginobject{GfsConstant}
 */

static void gfs_constant_destroy (GtsObject * object)
{
  if (GFS_CONSTANT (object)->derived)
    gfs_domain_remove_derived_variable (GFS_DOMAIN (gfs_object_simulation (object)), 
					GFS_CONSTANT (object)->derived->name);

  (* GTS_OBJECT_CLASS (gfs_constant_class ())->parent_class->destroy) (object);
}

static gdouble constant_func (FttCell * cell, FttCellFace * face, GfsDomain * domain,
			      GfsConstant * c)
{
  return c->val;
}

static void gfs_constant_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_constant_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GfsDerivedVariableInfo v = {
    fp->token->str, "Constant",
    constant_func, *o
  };
  GFS_CONSTANT (*o)->derived = 
    gfs_domain_add_derived_variable (GFS_DOMAIN (gfs_object_simulation (*o)), v);
  if (!GFS_CONSTANT (*o)->derived) {
    gts_file_error (fp, "'%s' keyword already used", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_constant_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_constant_class ())->parent_class->write) (o, fp);  
  fprintf (fp, " %s", GFS_CONSTANT (o)->derived->name);
}

static void gfs_constant_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_constant_destroy;
  klass->read =    gfs_constant_read;
  klass->write =   gfs_constant_write;
}

static void gfs_constant_init (GfsEvent * event)
{
  event->start = -1;
  event->istep = 1;
}

GfsEventClass * gfs_constant_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_constant_info = {
      "GfsConstant",
      sizeof (GfsConstant),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_constant_class_init,
      (GtsObjectInitFunc) gfs_constant_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()), &gfs_constant_info);
  }

  return klass;
}

/** \endobject{GfsConstant} */

/**
 * Compute the spatial sum of a GfsFunction.
 * \beginobject{GfsSpatialSum}
 */

static void gfs_spatial_sum_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SPATIAL_SUM (o)->v));
  (* GTS_OBJECT_CLASS (gfs_spatial_sum_class ())->parent_class->destroy) (o);
}

static void gfs_spatial_sum_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_spatial_sum_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_SPATIAL_SUM (o)->v, fp);
}

static void gfs_spatial_sum_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_spatial_sum_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR) 
    return;
  gfs_function_read (GFS_SPATIAL_SUM (*o)->v, GFS_DOMAIN (gfs_object_simulation (*o)), fp);
}

static void add (FttCell * cell, GfsSpatialSum * s)
{
  GFS_CONSTANT (s)->val += gfs_cell_volume (cell, GFS_DOMAIN (gfs_object_simulation (s)))*
    gfs_function_value (s->v, cell);
}

static gboolean gfs_spatial_sum_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_spatial_sum_class ())->parent_class)->event) 
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);

    GFS_CONSTANT (event)->val = 0.;
    gfs_catch_floating_point_exceptions ();
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add, event);
    gfs_restore_fpe_for_function (GFS_SPATIAL_SUM (event)->v);
    gfs_all_reduce (domain, GFS_CONSTANT (event)->val, MPI_DOUBLE, MPI_SUM);
    GFS_CONSTANT (event)->val *= pow (sim->physical_params.L, FTT_DIMENSION);
    return TRUE;
  }
  return FALSE;
}

static void gfs_spatial_sum_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_spatial_sum_read;
  klass->write = gfs_spatial_sum_write;
  klass->destroy = gfs_spatial_sum_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_spatial_sum_event;
}

static void gfs_spatial_sum_init (GfsSpatialSum * object)
{
  object->v = gfs_function_new (gfs_function_class (), 0.);
}

GfsEventClass * gfs_spatial_sum_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_spatial_sum_info = {
      "GfsSpatialSum",
      sizeof (GfsSpatialSum),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_spatial_sum_class_init,
      (GtsObjectInitFunc) gfs_spatial_sum_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_constant_class ()),
				  &gfs_spatial_sum_info);
  }

  return klass;
}

/** \endobject{GfsSpatialSum} */

/**
 * Derived variables.
 * \beginobject{GfsDerivedVariable}
 */

static void gfs_derived_variable_destroy (GtsObject * object)
{
  g_free (GFS_DERIVED_VARIABLE (object)->name);
  g_free (GFS_DERIVED_VARIABLE (object)->description);

  (* GTS_OBJECT_CLASS (gfs_derived_variable_class ())->parent_class->destroy) (object);
}

static void gfs_derived_variable_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_derived_variable_destroy;
}

GtsObjectClass * gfs_derived_variable_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_derived_variable_info = {
      "GfsDerivedVariable",
      sizeof (GfsDerivedVariable),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_derived_variable_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()), 
				  &gfs_derived_variable_info);
  }

  return klass;
}

/**
 * gfs_derived_variable_from_name:
 * @i: a list of #GfsDerivedVariable.
 * @name: a name.
 *
 * Returns: the #GfsDerivedVariable @name of @list or %NULL.
 */
GfsDerivedVariable * gfs_derived_variable_from_name (GSList * i, const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (i) {
    GfsDerivedVariable * v = i->data;
    if (!strcmp (v->name, name))
      return v;
    i = i->next;
  }
  return NULL;
}

/** \endobject{GfsDerivedVariable} */
