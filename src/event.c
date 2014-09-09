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
 * \brief GfsEvents.
 */

#include <stdlib.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "event.h"
#include "solid.h"
#include "output.h"
#include "init.h"

/**
 * Any action to be performed at a given time.
 * \beginobject{GfsEvent}
 */

/**
 * gfs_event_next:
 * @event: a #GfsEvent.
 * @sim: a #GfsSimulation.
 *
 * Returns: the next physical time at which @event will be realised.
 */
gdouble gfs_event_next (GfsEvent * event, GfsSimulation * sim)
{
  g_return_val_if_fail (event != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (sim != NULL, G_MAXDOUBLE);

  if (sim->time.t < event->t)
    return event->t;
  if (event->t >= event->end ||
      event->i >= event->iend ||
      sim->time.t > event->end || 
      sim->time.i > event->iend)
    return G_MAXDOUBLE;
  if (event->end_event)
    return G_MAXDOUBLE;
  if (sim->time.t >= event->t) {
    if (event->istep < G_MAXINT) {
      if (event->n == 0)
	return G_MAXDOUBLE;
    }
    else
      return event->start + (event->n + 1)*event->step;
  }
  if (sim->time.i >= event->i && event->step < G_MAXDOUBLE && event->n == 0)
    return sim->time.t + event->step;
  return G_MAXDOUBLE;
}

static gboolean gfs_event_event (GfsEvent * event, GfsSimulation * sim)
{
  if (event->redo) {
    event->redo = FALSE;
    return event->realised;
  }
  if (event->t >= event->end ||
      event->i >= event->iend ||
      sim->time.t > event->end || 
      sim->time.i > event->iend) {
    gts_object_destroy (GTS_OBJECT (event));
    return FALSE;
  }
  if (event->end_event) {
    if (event->n == 0 &&
	(sim->time.t >= sim->time.end ||
	 sim->time.i >= sim->time.iend)) {
      event->n = 1;
      return (event->realised = TRUE);
    }
    else
      return (event->realised = FALSE);
  }
  if (sim->time.t >= event->t) {
    if (event->istep < G_MAXINT) {
      if (event->n == 0) {
	event->i = sim->time.i + event->istep;
	event->n++;
	return (event->realised = TRUE);
      }
    }
    else {
      event->n++;
      event->t = event->start + event->n*event->step;
      return (event->realised = TRUE);
    }
  }
  if (sim->time.i >= event->i) {
    if (event->step < G_MAXDOUBLE) {
      if (event->n == 0) {
	event->start = sim->time.t;
	event->t = event->start + event->step;
	event->n = 1;
	return (event->realised = TRUE);
      }
    }
    else {
      event->n++;
      event->i += event->istep;
      return (event->realised = TRUE);
    }
  }
  return (event->realised = FALSE);
}

static void gfs_event_write (GtsObject * object, FILE * fp)
{
  GfsEvent * event = GFS_EVENT (object);

  fprintf (fp, "%s", object->klass->info.name);
  if (event->name)
    fprintf (fp, " *%s", event->name);
  fputs (" { ", fp);
  if (event->end_event)
    fputs ("start = end ", fp);
  else {
    if (event->start > 0. && event->start < G_MAXDOUBLE/2.)
      fprintf (fp, "start = %g ", event->start);
    if (event->step < G_MAXDOUBLE)
      fprintf (fp, "step = %g ", event->step);
    if (event->end < G_MAXDOUBLE)
      fprintf (fp, "end = %g ", event->end);
    if (event->istart > 0 && event->istart < G_MAXINT/2)
      fprintf (fp, "istart = %u ", event->istart);
    if (event->istep < G_MAXINT)
      fprintf (fp, "istep = %u ", event->istep);
    if (event->iend < G_MAXINT)
      fprintf (fp, "iend = %u ", event->iend);
  }
  fputc ('}', fp);
}

static void event_init (GfsEvent * object)
{
  object->t      = 0.;
  object->start  = 0.;
  object->end    = G_MAXDOUBLE;
  object->step   = G_MAXDOUBLE;

  object->i      = 0;
  object->istart = 0;
  object->iend   = G_MAXINT;
  object->istep  = G_MAXINT;

  object->n         = 0;
  object->end_event = FALSE;
}

static void gfs_event_read (GtsObject ** o, GtsFile * fp)
{
  GfsEvent * event = GFS_EVENT (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsEventClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_event_class ())) {
    gts_file_error (fp, "`%s' is not a GfsEvent", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (event));
    event = GFS_EVENT (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  /* optional name */
  if (fp->type == GTS_STRING && fp->token->str[0] == '*') {
    event->name = g_strdup (&(fp->token->str[1]));
    g_hash_table_insert (GFS_DOMAIN (gfs_object_simulation (*o))->objects, event->name, *o);
    gts_file_next_token (fp);  
  }
  
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "start",  TRUE},
      {GTS_DOUBLE, "end",    TRUE},
      {GTS_DOUBLE, "step",   TRUE},
      {GTS_UINT,   "istart", TRUE},
      {GTS_UINT,   "iend",   TRUE},
      {GTS_UINT,   "istep",  TRUE},
      {GTS_NONE}
    };
    gchar * start = NULL;

    var[0].data = &start;
    var[1].data = &event->end;
    var[2].data = &event->step;

    var[3].data = &event->istart;
    var[4].data = &event->iend;
    var[5].data = &event->istep;
 
    gts_file_assign_variables (fp, var);

    if (fp->type == GTS_ERROR)
      return;

    if (start) {
      if (!strcmp (start, "end")) {
	event->end_event = TRUE;
	if (var[1].set)
	  gts_file_variable_error (fp, var, "end", 
				   "end cannot be set for an `end' event");
	else if (var[2].set)
	  gts_file_variable_error (fp, var, "step", 
				   "step cannot be set for an `end' event");
	else if (var[3].set)
	  gts_file_variable_error (fp, var, "istart", 
				   "istart cannot be set for an `end' event");
	else if (var[4].set)
	  gts_file_variable_error (fp, var, "iend", 
				   "iend cannot be set for an `end' event");
	else if (var[5].set)
	  gts_file_variable_error (fp, var, "istep", 
				   "istep cannot be set for an `end' event");
      }
      else
	event->start = atof (start);
      g_free (start);
    }

    if (fp->type == GTS_ERROR)
      return;

    if (var[2].set && var[5].set) {
      gts_file_variable_error (fp, var, "istep", 
			       "step and istep cannot be set simultaneously");
      return;
    }
    if (var[2].set)
      event->istep = G_MAXINT;

    if (var[2].set && event->step <= 0.) {
      gts_file_variable_error (fp, var, "step",
			       "step `%g' must be strictly positive", 
			       event->step);
      return;
    }
    if (!var[2].set && !var[5].set && var[1].set) {
      gts_file_error (fp, "expecting a number (step or istep)");
      return;
    }
    if (var[1].set && event->end <= event->start) {
      gts_file_variable_error (fp, var, "end",
			       "end `%g' must be larger than start `%g'", 
			       event->end, event->start);
      return;
    }
    if (event->start < 0. && var[2].set)
      event->start = 0.;
    if (var[0].set || !var[3].set)
      event->t = event->start;
    else
      event->t = event->start = G_MAXDOUBLE/2.;

    if (!var[5].set && !var[2].set && var[4].set) {
      gts_file_error (fp, "expecting a number (istep or step)");
      return;
    }
    if (var[3].set && event->iend <= event->istart) {
      gts_file_variable_error (fp, var, "iend",
			       "iend `%u' must be larger than istart `%u'", 
			       event->iend, event->istart);
      return;
    }
    if (var[3].set || !var[0].set)
      event->i = event->istart;
    else
      event->i = event->istart = G_MAXINT/2;
  }

  if (class_changed && fp->type != '\n' && klass->read)
    (* klass->read) (o, fp);
}

static void gfs_event_destroy (GtsObject * o)
{
  GfsEvent * event = GFS_EVENT (o);
  
  if (event->name) {
    g_hash_table_remove (GFS_DOMAIN (gfs_object_simulation (o))->objects, event->name);
    g_free (event->name);
  }

  (* GTS_OBJECT_CLASS (gfs_event_class ())->parent_class->destroy) (o);
}

static void gfs_event_class_init (GfsEventClass * klass)
{
  klass->event = gfs_event_event;

  GTS_OBJECT_CLASS (klass)->write   = gfs_event_write;
  GTS_OBJECT_CLASS (klass)->read    = gfs_event_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_destroy;
}

GfsEventClass * gfs_event_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_info = {
      "GfsEvent",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_class_init,
      (GtsObjectInitFunc) event_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
			    &gfs_event_info);
  }

  return klass;
}

GfsEvent * gfs_event_new (GfsEventClass * klass)
{
  GfsEvent * object;

  g_return_val_if_fail (klass != NULL, NULL);

  object = GFS_EVENT (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/**
 * gfs_event_set:
 * @e: a #GfsEvent.
 * @start: start time.
 * @end: end time.
 * @step: time step.
 * @istart: start iteration.
 * @iend: end iteration.
 * @istep: iteration step.
 *
 * Sets the properties of event @e.
 *
 * If any of the arguments is negative, the corresponding value in @e
 * is unchanged.
 */
void gfs_event_set (GfsEvent * e,
		    gdouble start, gdouble end, gdouble step,
		    gint istart, gint iend, gint istep)
{
  g_return_if_fail (e != NULL);
  g_return_if_fail (end < 0. || start < 0. || start <= end);
  g_return_if_fail (istep >= 0 || step >= 0. || iend < 0);
  g_return_if_fail (istart < 0 || iend < 0 || istart <= iend);
  
  if (start >= 0.) e->start = start;
  if (end >= 0.)   e->end = end;
  if (step >= 0.)  e->step = step;
  if (istart >= 0) e->istart = istart;
  if (iend >= 0)   e->iend = iend;
  if (istep >= 0)  e->istep = istep;
  if (start >= 0. || istart < 0)
    e->t = e->start;
  else
    e->t = e->start = G_MAXDOUBLE/2.;
  if (istart >= 0 || start < 0.)
    e->i = e->istart;
  else
    e->i = e->istart = G_MAXINT/2;
}

/**
 * gfs_event_init:
 * @event: a #GfsEvent.
 * @sim: a #GfsSimulation.
 *
 * Initalizes @event associated with @sim. In particular, if @event is
 * an "init" event it is activated by this function.
 */
void gfs_event_init (GfsEvent * event,
		     GfsSimulation * sim)
{
  g_return_if_fail (event != NULL);
  g_return_if_fail (sim != NULL);

  if (event->start < 0.) { /* "init" event */
    g_assert (GFS_EVENT_CLASS (GTS_OBJECT (event)->klass)->event);
    (* GFS_EVENT_CLASS (GTS_OBJECT (event)->klass)->event) (event, sim);
  }
  else if (event->end_event)
    event->t = event->start = G_MAXDOUBLE/2.;
  else {
    if (event->istep < G_MAXINT)
      while (event->i < sim->time.i) {
	event->n++;
	event->i += event->istep;
      }
    else
      while (event->t < sim->time.t) {
	event->n++;
	event->t = event->start + event->n*event->step;
      }
  }
}

/**
 * gfs_event_do:
 * @event: a #GfsEvent:
 * @sim: a #GfsSimulation.
 * 
 * Realises the event if active.
 */
void gfs_event_do (GfsEvent * event, GfsSimulation * sim)
{
  GfsEventClass * klass;

  g_return_if_fail (event != NULL);
  g_return_if_fail (sim != NULL);
  
  gchar * name = GTS_OBJECT (event)->klass->info.name;
  gfs_domain_timer_start (GFS_DOMAIN (sim), name);

  klass = GFS_EVENT_CLASS (GTS_OBJECT (event)->klass);
  g_assert (klass->event);
  if ((* klass->event) (event, sim) && klass->post_event)
    (* klass->post_event) (event, sim);

  gfs_domain_timer_stop (GFS_DOMAIN (sim), name);
}

/**
 * gfs_event_half_do:
 * @event: a #GfsEvent:
 * @sim: a #GfsSimulation.
 * 
 * Realises the half-event if active.
 */
void gfs_event_half_do (GfsEvent * event, GfsSimulation * sim)
{
  g_return_if_fail (event != NULL);
  g_return_if_fail (sim != NULL);

  if (event->realised && GFS_EVENT_CLASS (GTS_OBJECT (event)->klass)->event_half)
    (* GFS_EVENT_CLASS (GTS_OBJECT (event)->klass)->event_half) (event, sim);
}

/**
 * gfs_event_redo:
 * @event: a #GfsEvent:
 * @sim: a #GfsSimulation.
 * 
 * Realises the event if it has just been realised.
 */
void gfs_event_redo (GfsEvent * event, GfsSimulation * sim)
{
  GfsEventClass * klass;

  g_return_if_fail (event != NULL);
  g_return_if_fail (sim != NULL);

  klass = GFS_EVENT_CLASS (GTS_OBJECT (event)->klass);
  g_assert (klass->event);
  event->redo = TRUE;
  if ((* klass->event) (event, sim) && klass->post_event)
    (* klass->post_event) (event, sim);
}

/** \endobject{GfsEvent} */

/**
 * Abstract class for initialisation events.
 * \beginobject{GfsGenericInit}
 */

static void gfs_generic_init_init (GfsEvent * event)
{
  event->start = -1;
}

GfsEventClass * gfs_generic_init_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_generic_init_info = {
      "GfsGenericInit",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_generic_init_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_generic_init_info);
  }

  return klass;
}

/** \endobject{GfsGenericInit} */

/**
 * Initialising variables.
 * \beginobject{GfsInit}
 */

typedef struct {
  GfsVariable * v[FTT_DIMENSION];
  GfsFunction * f[FTT_DIMENSION];
  guint n;
} VarFunc;

static VarFunc * var_func_new (GfsVariable ** v, GfsFunction ** f, int n)
{
  VarFunc * vf = g_malloc (sizeof (VarFunc));
  gint i;
  for (i = 0; i < n; i++) {
    vf->v[i] = v[i];
    vf->f[i] = f[i];
  }
  vf->n = n;
  return vf;
}

static void var_func_destroy (VarFunc * v)
{
  gint i;
  for (i = 0; i < v->n; i++)
    gts_object_destroy (GTS_OBJECT (v->f[i]));
  g_free (v);
}

static void gfs_init_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
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
    GfsInit * init = GFS_INIT (*o);
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
      
    if (fp->type == '(') {
      /* vector */
      GfsVariable * v[FTT_DIMENSION];
      if (!gfs_read_variable_vector (fp, v, domain))
	return;

      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);

      GfsFunction * f[FTT_DIMENSION];
      if (!gfs_read_function_vector (fp, v, f, domain))
	return;

      init->f = g_slist_append (init->f, var_func_new (v, f, FTT_DIMENSION));
    }
    else if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a variable name");
      return;
    }
    else {
      /* scalar */
      GfsVariable * v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL);
      GfsFunction * f;

      if (!v) {
	gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
	return;
      }
      gts_file_next_token (fp);

      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);

      f = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (f, gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (GTS_OBJECT (f));
	return;
      }
      gfs_function_set_units (f, v->units);
      init->f = g_slist_append (init->f, var_func_new (&v, &f, 1));
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  GfsEvent * event = GFS_EVENT (*o);
  if (event->start < 0. && (event->istep < G_MAXINT || event->step < G_MAXDOUBLE))
    event->start = 0.;
}

static void gfs_init_write (GtsObject * o, FILE * fp)
{
  GSList * i;
  
  if (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write) 
      (o, fp);

  fputs (" {\n", fp);
  i = GFS_INIT (o)->f;
  while (i) {
    VarFunc * v = i->data;
    if (v->n == 1) {
      /* scalar */
      fprintf (fp, "  %s =", v->v[0]->name);
      gfs_function_write (v->f[0], fp);
      fputc ('\n', fp);
    }
    else {
      /* vector */
      fprintf (fp, "  (%s", v->v[0]->name);
      gint i;
      for (i = 1; i < v->n; i++)
	fprintf (fp, ",%s", v->v[i]->name);
      fputs (") = (", fp);
      gfs_function_write (v->f[0], fp);
      for (i = 1; i < v->n; i++) {
	fputc (',', fp);
	gfs_function_write (v->f[i], fp);
      }
      fputs (")\n", fp);
    }
    i = i->next;
  }
  fputc ('}', fp);
}

static void gfs_init_destroy (GtsObject * object)
{
  GfsInit * i = GFS_INIT (object);

  g_slist_foreach (i->f, (GFunc) var_func_destroy, NULL);
  g_slist_free (i->f);

  (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->destroy) 
    (object);
}

static void init_scalar (FttCell * cell, VarFunc * vf)
{
  GFS_VALUE (cell, vf->v[0]) = gfs_function_value (vf->f[0], cell);
}

static void init_vector (FttCell * cell, VarFunc * vf)
{
  FttVector p, u;
  int i;
  for (i = 0; i < FTT_DIMENSION; i++)
    (&u.x)[i] = gfs_function_value (vf->f[i], cell);
  ftt_cell_pos (cell, &p);
  gfs_simulation_map_vector (GFS_SIMULATION (vf->v[0]->domain), &p, &u);
  for (i = 0; i < FTT_DIMENSION; i++)
    GFS_VALUE (cell, vf->v[i]) = (&u.x)[i];
}

static gboolean gfs_init_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class)->event) 
      (event, sim)) {
    GSList * i = GFS_INIT (event)->f;

    while (i) {
      VarFunc * vf = i->data;
      gfs_catch_floating_point_exceptions ();
      /* fixme: the check for "layered" variables is messy */
      if (!gfs_char_in_string (vf->v[0]->name[strlen (vf->v[0]->name) - 1], "0123456789"))
	gfs_domain_traverse_layers (GFS_DOMAIN (sim), (FttCellTraverseFunc) 
				    (vf->n == 1 ? init_scalar : init_vector), vf);
      else
	gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) 
				    (vf->n == 1 ? init_scalar : init_vector), vf);
      gfs_restore_fpe_for_function (vf->f[0]);
      if (vf->v[0]->component == FTT_DIMENSION)
	gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, vf->v[0]);
      i = i->next;
    }
    /* boundary conditions for vector quantities need to be called in
       a separate loop so that they are correctly applied for weird
       topologies (e.g. cubed sphere) */
    i = GFS_INIT (event)->f;
    while (i) {
      VarFunc * vf = i->data;
      gint j;
      for (j = 0; j < vf->n; j++)
	if (vf->v[j]->component < FTT_DIMENSION)
	  gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, vf->v[j]);
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_destroy;
}

GfsGenericInitClass * gfs_init_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_info = {
      "GfsInit",
      sizeof (GfsInit),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_info);
  }

  return klass;
}

/** \endobject{GfsInit} */

/**
 * "Masking out" part of the mesh.
 * \beginobject{GfsInitMask}
 */

static void gfs_init_mask_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_mask_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (GFS_INIT_MASK (*o)->mask, GFS_DOMAIN (gfs_object_simulation (*o)), fp);
}

static void gfs_init_mask_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_mask_class ())->parent_class->write) (o, fp);

  gfs_function_write (GFS_INIT_MASK (o)->mask, fp);
}

static void gfs_init_mask_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_INIT_MASK (object)->mask));

  (* GTS_OBJECT_CLASS (gfs_init_mask_class ())->parent_class->destroy) (object);
}

#define MASKED (1 << GFS_FLAG_USER)

static void mask_cells (FttCell * cell, GfsInitMask * m)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	mask_cells (child.c[i], m);
    if (FTT_CELL_IS_LEAF (cell))
      /* all the children have been destroyed i.e. the cell is masked */
      cell->flags |= MASKED;
  }
  else if (gfs_function_value (m->mask, cell) > 0.)
    cell->flags |= MASKED;
  
  if ((cell->flags & MASKED) && !FTT_CELL_IS_ROOT (cell))
    ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup, gfs_object_simulation (m));
}

static void foreach_box (GfsBox * box, GfsInitMask * m)
{
  mask_cells (box->root, m);
  if (box->root->flags & MASKED)
    m->masked_boxes = g_slist_prepend (m->masked_boxes, box);
}

static gboolean gfs_init_mask_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_mask_class ())->parent_class)->event) 
      (event, sim)) {
    GfsInitMask * m = GFS_INIT_MASK (event);
    m->masked_boxes = NULL;
    gfs_catch_floating_point_exceptions ();
    gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) foreach_box, m);
    gfs_restore_fpe_for_function (m->mask);
    g_slist_foreach (m->masked_boxes, (GFunc) gts_object_destroy, NULL);
    g_slist_free (m->masked_boxes);
    gfs_domain_match (GFS_DOMAIN (sim));
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_mask_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_mask_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_mask_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_mask_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_mask_destroy;
}

static void gfs_init_mask (GfsInitMask * m)
{
  m->mask = gfs_function_new (gfs_function_class (), 0.);
}

GfsGenericInitClass * gfs_init_mask_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_mask_info = {
      "GfsInitMask",
      sizeof (GfsInitMask),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_mask_class_init,
      (GtsObjectInitFunc) gfs_init_mask,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_mask_info);
  }

  return klass;
}

/** \endobject{GfsInitMask} */

/* GfsInitFlowConstant: Object: fixme: deprecated */

static void gfs_init_flow_constant_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_flow_constant_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_flow_constant_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_warning ("GfsInitFlowConstant is deprecated you should use GfsInit instead");
}

static void gfs_init_flow_constant_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_flow_constant_read;
}

GfsEventClass * gfs_init_flow_constant_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_flow_constant_info = {
      "GfsInitFlowConstant",
      sizeof (GfsInit),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_flow_constant_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_init_class ()),
				  &gfs_init_flow_constant_info);
  }

  return klass;
}

/**
 * Initialising a vorticity field.
 * \beginobject{GfsInitVorticity}
 */

static void gfs_init_vorticity_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
#if FTT_2D
  gfs_function_read (GFS_INIT_VORTICITY (*o)->f, gfs_object_simulation (*o), fp);
#else /* 3D */
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    gfs_function_read (GFS_INIT_VORTICITY (*o)->fv[c], gfs_object_simulation (*o), fp);
#endif /* 3D */
}

static void gfs_init_vorticity_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class->write) (o, fp);
#if FTT_2D
  gfs_function_write (GFS_INIT_VORTICITY (o)->f, fp);
#else /* 3D */
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    gfs_function_write (GFS_INIT_VORTICITY (o)->fv[c], fp);
#endif /* 3D */
}

static void gfs_init_vorticity_destroy (GtsObject * object)
{
#if FTT_2D
  gts_object_destroy (GTS_OBJECT (GFS_INIT_VORTICITY (object)->f));
#else /* 3D */
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    gts_object_destroy (GTS_OBJECT (GFS_INIT_VORTICITY (object)->fv[c]));
#endif /* 3D */

  (* GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class->destroy) (object);
}

static void sum_volume (FttCell * cell, GtsRange * vol)
{
  gdouble size = ftt_cell_size (cell);
  
  if (GFS_IS_MIXED (cell))
    gts_range_add_value (vol, size*size*GFS_STATE (cell)->solid->a);
  else
    gts_range_add_value (vol, size*size);
}

static void add_ddiv (FttCell * cell, gpointer * data)
{
  gdouble * ddiv = data[0];
  GfsVariable * div = data[1];
  gdouble size = ftt_cell_size (cell);
  
  if (GFS_IS_MIXED (cell))
    GFS_VALUE (cell, div) += size*size*GFS_STATE (cell)->solid->a*(*ddiv);
  else
    GFS_VALUE (cell, div) += size*size*(*ddiv);
}

static void correct_div (GfsDomain * domain, GfsVariable * v)
{
  GtsRange div, vol;
  gdouble ddiv;
  gpointer data[2];

  div = gfs_domain_stats_variable (domain, v, FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
  gts_range_init (&vol);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum_volume, &vol);
  gts_range_update (&vol);
  ddiv = - div.mean/vol.mean;

  data[0] = &ddiv;
  data[1] = v;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) add_ddiv, data);
}

static void stream_from_vorticity (GfsDomain * domain,
				   GfsVariable * stream,
				   GfsVariable * vorticity,
				   gdouble tolerance)
{
  GfsNorm norm;
  guint maxit = 100;
  GfsVariable * res, * dia;
  GfsMultilevelParams par;

  g_return_if_fail (domain != NULL);

  dia = gfs_temporary_variable (domain);
  gfs_poisson_coefficients (domain, NULL, TRUE, TRUE, TRUE);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, dia);
  correct_div (domain, vorticity); /* enforce solvability condition */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, stream);
  res = gfs_temporary_variable (domain);
  gfs_residual (domain, FTT_DIMENSION, FTT_TRAVERSE_LEAFS, -1, stream, vorticity, dia, res);
  norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1., res);
  gfs_multilevel_params_init (&par);
  par.depth = gfs_domain_depth (domain);
  while (norm.infty > tolerance && maxit) {
    gfs_poisson_cycle (domain, &par, stream, vorticity, dia, res);
    norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1., res);
    // fprintf (stderr, "maxit: %d infty: %g\n", maxit, norm.infty);
    maxit--;
  }
  if (maxit == 0)
    g_warning ("GfsInitVorticity: cannot solve streamfunction from vorticity\n"
	       "  (residual: %g)", norm.infty);
  gts_object_destroy (GTS_OBJECT (res));
  gts_object_destroy (GTS_OBJECT (dia));
}

static void init_from_streamfunction (FttCell * cell, GfsInitVorticity * init)
{
  gdouble size = ftt_cell_size (cell);
#if FTT_2D
  GFS_VALUE (cell, init->u[0]) = - gfs_center_gradient (cell, FTT_Y, init->stream->i)/size;
  GFS_VALUE (cell, init->u[1]) = gfs_center_gradient (cell, FTT_X, init->stream->i)/size;
#else /* 3D */
  GFS_VALUE (cell, init->u[0]) = (gfs_center_gradient (cell, FTT_Z, init->stream[1]->i) -
				  gfs_center_gradient (cell, FTT_Y, init->stream[2]->i))/size;
  GFS_VALUE (cell, init->u[1]) = (gfs_center_gradient (cell, FTT_X, init->stream[2]->i) -
				  gfs_center_gradient (cell, FTT_Z, init->stream[0]->i))/size;
  GFS_VALUE (cell, init->u[2]) = (gfs_center_gradient (cell, FTT_Y, init->stream[0]->i) -
				  gfs_center_gradient (cell, FTT_X, init->stream[1]->i))/size;
#endif /* 3D */
}

static void compute_vorticity (FttCell * cell, GfsInitVorticity * init)
{
  gdouble size = ftt_cell_size (cell);

  GFS_VALUE (cell, init->vort) = gfs_function_value (init->f, cell)*size*size;
}

static gboolean gfs_init_vorticity_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class)->event) 
      (event, sim)) {
    GfsInitVorticity * init = GFS_INIT_VORTICITY (event);
    GfsDomain * domain = GFS_DOMAIN (sim);

#if FTT_2D
    init->vort = gfs_temporary_variable (domain);
    init->stream = gfs_temporary_variable (domain);
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_vorticity, event);
    gfs_restore_fpe_for_function (init->f);
    stream_from_vorticity (domain, init->stream, init->vort, 1e-9);
    gts_object_destroy (GTS_OBJECT (init->vort));
    init->u = gfs_domain_velocity (domain);
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_from_streamfunction, init);
    gts_object_destroy (GTS_OBJECT (init->stream));
#else /* 3D */
    FttComponent c;
    init->vort = gfs_temporary_variable (domain);
    for (c = 0; c < FTT_DIMENSION; c++) {
      init->stream[c] = gfs_temporary_variable (domain);
      gfs_catch_floating_point_exceptions ();
      init->f = init->fv[c];
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) compute_vorticity, event);
      gfs_restore_fpe_for_function (init->f);
      stream_from_vorticity (domain, init->stream[c], init->vort, 1e-9);
    }
    gts_object_destroy (GTS_OBJECT (init->vort));
    init->u = gfs_domain_velocity (domain);
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_from_streamfunction, init);
    for (c = 0; c < FTT_DIMENSION; c++)
      gts_object_destroy (GTS_OBJECT (init->stream[c]));
#endif /* 3D */
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_vorticity_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_init_vorticity_read;
  klass->write = gfs_init_vorticity_write;
  klass->destroy = gfs_init_vorticity_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_init_vorticity_event;
}

static void gfs_init_vorticity_init (GfsInitVorticity * init)
{
#if FTT_2D
  init->f = gfs_function_new (gfs_function_class (), 0.);
#else /* 3D */
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    init->fv[c] = gfs_function_new (gfs_function_class (), 0.);
#endif /* 3D */
}

GfsGenericInitClass * gfs_init_vorticity_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_vorticity_info = {
      "GfsInitVorticity",
      sizeof (GfsInitVorticity),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_vorticity_class_init,
      (GtsObjectInitFunc) gfs_init_vorticity_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_vorticity_info);
  }

  return klass;
}

/** \endobject{GfsInitVorticity} */

/**
 * Time-integral of a variable or function.
 * \beginobject{GfsEventSum}
 */

static void gfs_event_sum_destroy (GtsObject * o)
{
  GfsEventSum * s = GFS_EVENT_SUM (o);

  gts_object_destroy (GTS_OBJECT (s->v));

  (* GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->destroy) (o);
}

static void gfs_event_sum_write (GtsObject * o, FILE * fp)
{
  GfsEventSum * s = GFS_EVENT_SUM (o);

  (* GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->write) (o, fp);

  gfs_function_write (s->v, fp);
  fprintf (fp, " %s", s->sv->name);
}

static void gfs_event_sum_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventSum * s = GFS_EVENT_SUM (*o);
  GfsDomain * domain =  GFS_DOMAIN (gfs_object_simulation (s));

  (* GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (s->v, domain, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (sv)");
    return;
  }
  if (!(s->sv = gfs_domain_get_or_add_variable (domain, fp->token->str, "Sum"))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static gboolean gfs_event_sum_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventSum * s = GFS_EVENT_SUM (event);

    if (s->last < 0.)
      gfs_domain_traverse_layers (GFS_DOMAIN (sim), (FttCellTraverseFunc) gfs_cell_reset, s->sv);
    else {
      s->dt = sim->time.t - s->last;
      gfs_catch_floating_point_exceptions ();
      gfs_domain_traverse_layers (GFS_DOMAIN (sim), s->sum, s);
      gfs_restore_fpe_for_function (s->v);
    }
    s->last = sim->time.t;
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_sum_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_event_sum_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_sum_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_sum_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_event_sum_event;
}

static void sum (FttCell * cell, GfsEventSum * s)
{
  GFS_VALUE (cell, s->sv) += s->dt*gfs_function_value (s->v, cell);
}

static void gfs_event_sum_init (GfsEventSum * object)
{
  object->last = -1.;
  object->v = gfs_function_new (gfs_function_class (), 0.);
  object->sum = (FttCellTraverseFunc) sum;
}

GfsEventClass * gfs_event_sum_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_sum_info = {
      "GfsEventSum",
      sizeof (GfsEventSum),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_sum_class_init,
      (GtsObjectInitFunc) gfs_event_sum_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_sum_info);
  }

  return klass;
}

/** \endobject{GfsEventSum} */

/**
 *
 * \beginobject{GfsEventSumDirection}
 */

static void gfs_event_sum_direction_write (GtsObject * o, FILE * fp)
{
  GfsEventSumDirection * s = GFS_EVENT_SUM_DIRECTION (o);

  (* GTS_OBJECT_CLASS (gfs_event_sum_direction_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", ftt_direction_name [s->d]);
}

static void gfs_event_sum_direction_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventSumDirection * s = GFS_EVENT_SUM_DIRECTION (*o);

  (* GTS_OBJECT_CLASS (gfs_event_sum_direction_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (direction)");
    return;
  }
  s->d = ftt_direction_from_name (fp->token->str);
  if (s->d >= FTT_NEIGHBORS) {
    gts_file_error (fp, "unknown direction `%s'", fp->token->str);
    s->d = 0;
    return;
  }
  gts_file_next_token (fp);
}

static gboolean gfs_event_sum_direction_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventSumDirection * s = GFS_EVENT_SUM_DIRECTION (event);

    gfs_domain_sum (GFS_DOMAIN (sim), s->d, GFS_EVENT_SUM (event)->v, GFS_EVENT_SUM (event)->sv);
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_sum_direction_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_event_sum_direction_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_sum_direction_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_sum_direction_event;
}

GfsEventClass * gfs_event_sum_direction_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_sum_direction_info = {
      "GfsEventSumDirection",
      sizeof (GfsEventSumDirection),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_sum_direction_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_sum_class ()),
				  &gfs_event_sum_direction_info);
  }

  return klass;
}

/** \endobject{GfsEventSumDirection} */

/**
 * Harmonic decomposition of a variable.
 * \beginobject{GfsEventHarmonic}
 */

static void gfs_event_harmonic_destroy (GtsObject * o)
{
  GfsEventHarmonic * s = GFS_EVENT_HARMONIC (o);

  if (s->Mn)
    gfs_matrix_free (s->Mn);
  if (s->M)
    gfs_matrix_free (s->M);
  if (s->iM)
    gfs_matrix_free (s->iM);

  g_free (s->A);
  g_free (s->B);
  g_free (s->vsin);
  g_free (s->vcos);
  g_free (s->x);
  g_free (s->a);

  g_free (s->Aname);
  g_free (s->Bname);

  g_array_free (s->omega, TRUE);

  (* GTS_OBJECT_CLASS (gfs_event_harmonic_class ())->parent_class->destroy) (o);
}

static void gfs_event_harmonic_write (GtsObject * o, FILE * fp)
{
  GfsEventHarmonic * s = GFS_EVENT_HARMONIC (o);
  guint i, j;

  (* GTS_OBJECT_CLASS (gfs_event_harmonic_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %s %s %s", s->v->name, s->Aname, s->Bname, s->z->name);
  if (s->e)
    fprintf (fp, " %s", s->e->name);
  for (i = 0; i < s->omega->len; i++)
    fprintf (fp, " %.12lf", g_array_index (s->omega, gdouble, i));
  fprintf (fp, " { %d", s->invertible);
  for (i = 0; i < 2*s->omega->len + 1; i++)
    for (j = 0; j < 2*s->omega->len + 1; j++)
      fprintf (fp, " %.12lf", s->M[i][j]);
  fputs (" }", fp);
}

static void gfs_event_harmonic_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventHarmonic * s = GFS_EVENT_HARMONIC (*o);
  GfsDomain * domain =  GFS_DOMAIN (gfs_object_simulation (s));
  guint i;

  (* GTS_OBJECT_CLASS (gfs_event_harmonic_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  if (!(s->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (A)");
    return;
  }
  s->Aname = g_strdup (fp->token->str);  
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (B)");
    return;
  }
  s->Bname = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (Z)");
    return;
  }
  if (!(s->z = gfs_domain_get_or_add_variable (domain, fp->token->str, 
					       "Offset of harmonic decomposition"))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  s->z->units = s->v->units;
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (E)");
      return;
    }
    if (!(s->e = gfs_domain_get_or_add_variable (domain, fp->token->str, 
						 "Remainder of harmonic decomposition"))) {
      gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
      return;
    }
    s->e->units = 2.*s->v->units;
    gts_file_next_token (fp);
  }

  do {
    gdouble omega = gfs_read_constant (fp, domain);
    if (fp->type == GTS_ERROR)      
      return;
    g_array_append_val (s->omega, omega);
  } while (fp->type != '\n' && fp->type != '{');

  s->Mn = gfs_matrix_new (2*s->omega->len + 1, 2*s->omega->len + 1, sizeof (gdouble));
  for (i = 0; i < 2*s->omega->len + 1; i++)
    s->Mn[i][i] = 1.;

  s->M  = gfs_matrix_new (2*s->omega->len + 1, 2*s->omega->len + 1, sizeof (gdouble));
  s->iM = gfs_matrix_new (2*s->omega->len + 1, 2*s->omega->len + 1, sizeof (gdouble));

  s->A =    g_malloc (sizeof (GfsVariable *)*s->omega->len);
  s->B =    g_malloc (sizeof (GfsVariable *)*s->omega->len);
  s->vsin = g_malloc (sizeof (gdouble)*s->omega->len);
  s->vcos = g_malloc (sizeof (gdouble)*s->omega->len);
  s->x    = g_malloc (sizeof (gdouble)*(2*s->omega->len + 1));
  s->a    = g_malloc (sizeof (gdouble)*(2*s->omega->len + 1));

  for (i = 0; i < s->omega->len; i++) {
    gchar * u;
    
    u = g_strdup_printf ("%s%d", s->Aname, i);
    if (!(s->A[i] = gfs_domain_get_or_add_variable (domain, u, 
					"In-phase component of the harmonic decomposition"))) {
      gts_file_error (fp, "`%s' is a reserved keyword", u);
      return;
    }
    s->A[i]->units = s->v->units;
    g_free (u);
    u = g_strdup_printf ("%s%d", s->Bname, i);
    if (!(s->B[i] = gfs_domain_get_or_add_variable (domain, u,
				       "Out-of-phase component of the harmonic decomposition"))) {
      gts_file_error (fp, "`%s' is a reserved keyword", u);
      return;
    }
    s->B[i]->units = s->v->units;
    g_free (u);
  }

  if (fp->type == '{') {
    guint n = 2*s->omega->len + 1;
    guint j;

    fp->scope_max++;
    gts_file_next_token (fp);
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting a number (invertible)");
      return;
    }
    s->invertible = atoi (fp->token->str);
    gts_file_next_token (fp);
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	  gts_file_error (fp, "expecting a number (M[%d][%d])", i, j);
	  return;
	}
	else {
	  s->M[i][j] = atof (fp->token->str);
	  gts_file_next_token (fp);
	}
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    gts_file_next_token (fp);
    fp->scope_max--;

    if (s->invertible)
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  s->Mn[i][j] = s->M[i][j];
  }
}

static void add_xsin_xcos (FttCell * cell, GfsEventHarmonic * h)
{
  gdouble x = GFS_VALUE (cell, h->v);
  guint i;

  for (i = 0; i < h->omega->len; i++) {
    GFS_VALUE (cell, h->A[i]) += x*h->vcos[i];
    GFS_VALUE (cell, h->B[i]) += x*h->vsin[i];
  }
  GFS_VALUE (cell, h->z) += x;
  if (h->e)
    GFS_VALUE (cell, h->e) += x*x;
}

static gdouble de (GfsEventHarmonic * h, gdouble ** M)
{
  guint n = h->omega->len;
  gdouble xm = h->a[2*n];
  gdouble e = xm*(M[2*n][2*n]*xm - 2.*h->x[2*n]);
  guint i, j;

  for (i = 0; i < n; i++) {
    e += 2.*(h->a[i]*(xm*M[i][2*n] - h->x[i]) +
	     h->a[n + i]*(xm*M[n + i][2*n] - h->x[n + i]));
    for (j = 0; j < n; j++)
      e += (h->a[i]*h->a[j]*M[j][i] + 
	    h->a[n + i]*h->a[n + j]*M[n + j][n + i] +
	    2.*h->a[i]*h->a[n + j]*M[n + j][i]);
  }
  return e;
}

static void update_A_B_Z (FttCell * cell, GfsEventHarmonic * h)
{
  gdouble x = GFS_VALUE (cell, h->v), sx2 = 0.;
  guint n = h->omega->len;
  guint i, j;

  /* A^n */
  for (i = 0; i < n; i++) {
    h->a[i] =     GFS_VALUE (cell, h->A[i]);
    h->a[i + n] = GFS_VALUE (cell, h->B[i]);
  }
  h->a[2*n] = GFS_VALUE (cell, h->z);

  /* X^n = M^n.A^n */
  for (i = 0; i < 2*n + 1; i++) {
    h->x[i] = 0.;
    for (j = 0; j < 2*n + 1; j++)
      h->x[i] += h->Mn[i][j]*h->a[j];
  }

  if (h->e) {
    if (h->invertible)
      sx2 = x*x + h->Mn[2*n][2*n]*GFS_VALUE (cell, h->e) - de (h, h->Mn);
    else
      sx2 = x*x + GFS_VALUE (cell, h->e);
  }
  
  /* X^n+1 = X^n + Delta^n */
  for (i = 0; i < n; i++) {
    h->x[i]     += x*h->vcos[i];
    h->x[i + n] += x*h->vsin[i];
  }
  h->x[2*n] += x;

  /* A^n+1 = (M^n+1)^-1.X^n+1 */
  for (i = 0; i < 2*n + 1; i++) {
    h->a[i] = 0.;
    for (j = 0; j < 2*n + 1; j++)
      h->a[i] += h->iM[i][j]*h->x[j];
  }

  for (i = 0; i < n; i++) {
    GFS_VALUE (cell, h->A[i]) = h->a[i];
    GFS_VALUE (cell, h->B[i]) = h->a[i + n];
  }
  GFS_VALUE (cell, h->z) = h->a[2*n];

  if (h->e)
    GFS_VALUE (cell, h->e) = (sx2 + de (h, h->M))/h->M[2*n][2*n];
}

static gboolean gfs_event_harmonic_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_harmonic_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventHarmonic * h = GFS_EVENT_HARMONIC (event);
    gdouble ** M = h->M, * vsin = h->vsin, * vcos = h->vcos;
    gdouble ** iM = h->iM, ** Mn = h->Mn;
    guint i, j, n = h->omega->len;

    for (i = 0; i < n; i++) {
      vsin[i] = sin (g_array_index (h->omega, gdouble, i)*sim->time.t);
      vcos[i] = cos (g_array_index (h->omega, gdouble, i)*sim->time.t);
    }

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        M[i][j]         += vcos[j]*vcos[i];
	M[i][n + j]     += vsin[j]*vcos[i];
        M[n + i][j]     += vcos[j]*vsin[i];
	M[n + i][n + j] += vsin[j]*vsin[i];
      }
      M[i][2*n]     += vcos[i];
      M[n + i][2*n] += vsin[i];
    }
    for (j = 0; j < n; j++) {
      M[2*n][j]     += vcos[j];
      M[2*n][n + j] += vsin[j];
    }
    M[2*n][2*n] += 1.;

    for (i = 0; i < 2*n + 1; i++)
      for (j = 0; j < 2*n + 1; j++)
	iM[i][j] = M[i][j];

    if (!gfs_matrix_inverse (iM, 2*n + 1, 1e-6)) {
      g_assert (!h->invertible);
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) add_xsin_xcos, h);
    }
    else {
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) update_A_B_Z, h);
      h->invertible = TRUE;
      for (i = 0; i < 2*n + 1; i++)
	for (j = 0; j < 2*n + 1; j++)
	  Mn[i][j] = M[i][j];
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_harmonic_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_event_harmonic_destroy;
  klass->read = gfs_event_harmonic_read;
  klass->write = gfs_event_harmonic_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_harmonic_event;
}

static void gfs_event_harmonic_init (GfsEventHarmonic * object)
{
  object->omega = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

GfsEventClass * gfs_event_harmonic_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_harmonic_info = {
      "GfsEventHarmonic",
      sizeof (GfsEventHarmonic),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_harmonic_class_init,
      (GtsObjectInitFunc) gfs_event_harmonic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_harmonic_info);
  }

  return klass;
}

/** \endobject{GfsEventHarmonic} */

/**
 * Stopping when steady-state is reached.
 * \beginobject{GfsEventStop}
 */

static void gfs_event_stop_write (GtsObject * o, FILE * fp)
{
  GfsEventStop * s = GFS_EVENT_STOP (o);

  if (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->write)
      (o, fp);

  fprintf (fp, " %s %g", s->v->name, s->max);
  if (s->diff)
    fprintf (fp, " %s", s->diff->name);
  if (s->relative)
    fputs (" { relative = 1 }", fp);
}

static void gfs_event_stop_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventStop * s = GFS_EVENT_STOP (*o);
  GfsDomain * domain =  GFS_DOMAIN (gfs_object_simulation (s));

  if (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  if (!(s->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  s->max = gfs_read_constant (fp, domain);
  if (fp->type == GTS_ERROR)      
    return;
  s->oldv = gfs_domain_add_variable (domain, NULL, NULL);
  /* fixme: the lines below are necessary in the general case (e.g. when dealing with a VOF tracer)
   * but will crash if s->oldv is not of the same class as s->v.
   * s->oldv->fine_coarse = s->v->fine_coarse;
   * s->oldv->coarse_fine = s->v->coarse_fine;
   */

  if (fp->type == GTS_STRING) {
    if (!(s->diff = gfs_domain_get_or_add_variable (domain, fp->token->str, 
						    "Stopping field difference"))) {
      gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
      return;
    }
    s->diff->units = s->v->units;
    gts_file_next_token (fp);
  }

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "relative", TRUE, &s->relative},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);    
  }
}

static void gfs_event_stop_destroy (GtsObject * o)
{
  if (GFS_EVENT_STOP (o)->oldv)
    gts_object_destroy (GTS_OBJECT (GFS_EVENT_STOP (o)->oldv));

  (* GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->destroy) (o);
}

static void diff (FttCell * cell, GfsEventStop * s)
{
  GFS_VALUE (cell, s->oldv) -= GFS_VALUE (cell, s->v);
}

static void copy (FttCell * cell, GfsEventStop * s)
{
  GFS_VALUE (cell, s->oldv) = GFS_VALUE (cell, s->v);
}

static gboolean gfs_event_stop_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class)->event) 
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsEventStop * s = GFS_EVENT_STOP (event);

    if (s->last >= 0.) {
      GfsNorm n;
      gdouble max = 1.;

      if (s->relative) {
	n = gfs_domain_norm_variable (domain, s->oldv, NULL, FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
	max = n.infty;
      }
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) diff, s);
      n = gfs_domain_norm_variable (domain, s->oldv, NULL, FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
      if (s->relative) {
	if (n.infty <= s->max*max)
	  sim->time.end = sim->time.t;
      }
      else if (gfs_dimensional_value (s->v, n.infty) <= s->max)
	sim->time.end = sim->time.t;
      if (s->diff) {
	gfs_variables_swap (s->diff, s->oldv);
	gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, s->diff);
      }
    }
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) copy, s);
    gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, s->v, s->oldv);
    s->last = sim->time.t;
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_stop_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_event_stop_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_stop_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_stop_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_event_stop_event;
}

static void gfs_event_stop_init (GfsEventStop * object)
{
  object->last = -1.;
}

GfsEventClass * gfs_event_stop_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_stop_info = {
      "GfsEventStop",
      sizeof (GfsEventStop),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_stop_class_init,
      (GtsObjectInitFunc) gfs_event_stop_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_stop_info);
  }

  return klass;
}

/** \endobject{GfsEventStop} */

/**
 * Executing a shell-script at given times.
 * \beginobject{GfsEventScript}
 */

static void gfs_event_script_destroy (GtsObject * o)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (o);

  g_free (s->script);

  (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->destroy) (o);
} 

static void gfs_event_script_write (GtsObject * o, FILE * fp)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (o);

  if (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->write)
      (o, fp);

  fputs (" {", fp);
  if (s->script)
    fputs (s->script, fp);
  fputc ('}', fp);
}

static void gfs_event_script_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (*o);

  if (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_free (s->script);
  if ((s->script = gfs_file_statement (fp)))
    gts_file_next_token (fp);
}

/**
 * gfs_popen:
 * @sim: a #GfsSimulation.
 * @command: a shell command.
 * @type: "r" for reading or "w" for writing.
 *
 * Returns: a file descriptor pointing to a pipe opened using the
 * standard popen() command but which also defines standard Gerris
 * shell variables. In case of error returns %NULL.
 */
FILE * gfs_popen (GfsSimulation * sim, const char * command, const char * type)
{
  g_return_val_if_fail (command != NULL, NULL);
  g_return_val_if_fail (type != NULL, NULL);

  gchar * sname = gfs_template ();
  if (!gfs_mkftemp (sname)) {
    g_warning ("gfs_popen() cannot create FIFO: %s", strerror (errno));
    g_free (sname);
    return NULL;
  }
  /* When adding pre-defined shell variables please update this page:
     http://gfs.sourceforge.net/wiki/index.php/GfsEventScript */
  gchar * scommand = g_strdup_printf ("GfsTime=%g GfsIter=%d GfsPid=%d "
				      "GFS_STOP=%d sh %s",
				      sim->time.t, sim->time.i, 
				      GFS_DOMAIN (sim)->pid,
				      GFS_EVENT_SCRIPT_STOP,
				      sname);
  fflush (stdout);
  fflush (stderr);
  FILE * fp = popen (scommand, type);
  g_free (scommand);
  if (fp != NULL) {
    FILE * f = fopen (sname, "w");
    fputs (command, f);
    fclose (f);
  }
  else
    g_warning ("popen() command failed: %s", strerror (errno));
  remove (sname);
  g_free (sname);
  return fp;
}

static gboolean gfs_event_script_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventScript * s = GFS_EVENT_SCRIPT (event);
    if (s->script && GFS_DOMAIN (sim)->pid <= 0) {
      FILE * fp = gfs_popen (sim, s->script, "w");
      if (fp == NULL) {
	g_warning ("GfsEventScript cannot start script");
	return TRUE;
      }
      int status = pclose (fp);
      if (status != -1)
	status = WEXITSTATUS (status);
      if (status == GFS_EVENT_SCRIPT_STOP)
	exit (1);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_script_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_script_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_event_script_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_script_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_script_event;
}

GfsEventClass * gfs_event_script_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_script_info = {
      "GfsEventScript",
      sizeof (GfsEventScript),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_script_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_script_info);
  }

  return klass;
}

/** \endobject{GfsEventScript} */

/**
 * Initialising volume fractions from interface shapes.
 * \beginobject{GfsInitFraction}
 */

static void gfs_init_fraction_destroy (GtsObject * object)
{
  GfsInitFraction * init = GFS_INIT_FRACTION (object);

  gts_object_destroy (GTS_OBJECT (init->surface));

  (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->destroy) 
    (object);
}

static void gfs_init_fraction_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitFraction * init;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  init = GFS_INIT_FRACTION (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (init));
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }
  if ((init->c = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  gfs_generic_surface_read (init->surface, gfs_object_simulation (*o), fp);
}

static void gfs_init_fraction_write (GtsObject * o, FILE * fp)
{
  GfsInitFraction * init = GFS_INIT_FRACTION (o);


  (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", init->c->name);
  gfs_generic_surface_write (init->surface, gfs_object_simulation (o), fp);
}

static gboolean gfs_init_fraction_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_init_fraction (GFS_DOMAIN (sim), 
			      GFS_INIT_FRACTION (event)->surface,
			      GFS_INIT_FRACTION (event)->c);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_fraction_class_init (GfsInitFractionClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_fraction_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_fraction_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_fraction_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_fraction_destroy;
}

static void gfs_init_fraction_init (GfsInitFraction * object)
{
  object->surface = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
}

GfsInitFractionClass * gfs_init_fraction_class (void)
{
  static GfsInitFractionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_fraction_info = {
      "GfsInitFraction",
      sizeof (GfsInitFraction),
      sizeof (GfsInitFractionClass),
      (GtsObjectClassInitFunc) gfs_init_fraction_class_init,
      (GtsObjectInitFunc) gfs_init_fraction_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_fraction_info);
  }

  return klass;
}

/** \endobject{GfsInitFraction} */

/**
 * Controlling the topology of interfaces.
 * \beginobject{GfsRemoveDroplets}
 */

static void gfs_remove_droplets_destroy (GtsObject * object)
{
  GfsRemoveDroplets * d = GFS_REMOVE_DROPLETS (object);
  if (d->fc)
    gts_object_destroy (GTS_OBJECT (d->fc));
  (* GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->destroy) (object);
}

static void compute_v (FttCell * cell, GfsRemoveDroplets * d)
{
  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
}

static gboolean gfs_remove_droplets_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class)->event) 
      (event, sim)) {
    GfsRemoveDroplets * d = GFS_REMOVE_DROPLETS (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    d->v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
    if (d->v)
      gfs_domain_remove_droplets (domain, d->v, d->c, d->min, d->val);
    else {
      d->v = gfs_temporary_variable (domain);
      gfs_catch_floating_point_exceptions ();
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_v, d);
      gfs_restore_fpe_for_function (d->fc);
      gfs_domain_remove_droplets (domain, d->v, d->c, d->min, d->val);
      gts_object_destroy (GTS_OBJECT (d->v));
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_remove_droplets_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsRemoveDroplets * r = GFS_REMOVE_DROPLETS (*o);
  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (min)");
    return;
  }
  r->min = atoi (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != '\n') {
    r->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
    if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
      r->val = atof (fp->token->str);
      gts_file_next_token (fp);
    }
  }
}

static void gfs_remove_droplets_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->write) 
      (o, fp);
  GfsRemoveDroplets * r = GFS_REMOVE_DROPLETS (o);
  fprintf (fp, " %s %d", r->c->name, r->min);
  if (r->fc) {
    gfs_function_write (r->fc, fp);
    if (r->val != 0.)
      fprintf (fp, " %g", r->val);
  }
}

static void gfs_remove_droplets_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_remove_droplets_event;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_remove_droplets_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_remove_droplets_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_remove_droplets_write;
}

GfsEventClass * gfs_remove_droplets_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_remove_droplets_info = {
      "GfsRemoveDroplets",
      sizeof (GfsRemoveDroplets),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_remove_droplets_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_remove_droplets_info);
  }

  return klass;
}

/** \endobject{GfsRemoveDroplets} */

/**
 * Controlling the domain topology.
 * \beginobject{GfsRemovePonds}
 */

static gboolean gfs_remove_ponds_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_remove_ponds_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_remove_ponds (GFS_DOMAIN (sim), GFS_REMOVE_PONDS (event)->min,
			     (FttCellCleanupFunc) gfs_cell_cleanup, sim);
    return TRUE;
  }
  return FALSE;
}

static void gfs_remove_ponds_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_remove_ponds_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_remove_ponds_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (min)");
    return;
  }
  GFS_REMOVE_PONDS (*o)->min = atoi (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_remove_ponds_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_remove_ponds_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_remove_ponds_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %d", GFS_REMOVE_PONDS (o)->min);
}

static void gfs_remove_ponds_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_remove_ponds_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_remove_ponds_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_remove_ponds_write;
}

GfsEventClass * gfs_remove_ponds_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_remove_ponds_info = {
      "GfsRemovePonds",
      sizeof (GfsRemovePonds),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_remove_ponds_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_remove_ponds_info);
  }

  return klass;
}

/** \endobject{GfsRemovePonds} */

/**
 * Exponential filtering of variables.
 * \beginobject{GfsEventFilter}
 */

static void filter (FttCell * cell, GfsEventFilter * f)
{
  FttDirection d[4*(FTT_DIMENSION - 1)][FTT_DIMENSION] = {
#if FTT_2D
    {FTT_RIGHT, FTT_TOP}, {FTT_RIGHT, FTT_BOTTOM}, {FTT_LEFT, FTT_TOP}, {FTT_LEFT, FTT_BOTTOM}
#else
    {FTT_RIGHT, FTT_TOP, FTT_FRONT}, {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT}, 
    {FTT_LEFT, FTT_TOP, FTT_FRONT}, {FTT_LEFT, FTT_BOTTOM, FTT_FRONT},
    {FTT_RIGHT, FTT_TOP, FTT_BACK}, {FTT_RIGHT, FTT_BOTTOM, FTT_BACK}, 
    {FTT_LEFT, FTT_TOP, FTT_BACK}, {FTT_LEFT, FTT_BOTTOM, FTT_BACK}
#endif
  };
  guint i;
  gdouble val = 0.;

  for (i = 0; i < 4*(FTT_DIMENSION - 1); i++)
    val += gfs_cell_corner_value (cell, d[i], f->v, -1);
  GFS_VALUE (cell, f->tmp) = val/(4*(FTT_DIMENSION - 1));
}

static void filtered (FttCell * cell, GfsEventFilter * f)
{
  gdouble dt = gfs_object_simulation (f)->advection_params.dt/f->scale;
  GFS_VALUE (cell, f->v) = ((1. - dt)*GFS_VALUE (cell, f->v) +
			    dt*GFS_VALUE (cell, f->tmp));
}

static gboolean gfs_event_filter_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_filter_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventFilter * f = GFS_EVENT_FILTER (event);

    f->tmp = gfs_temporary_variable (GFS_DOMAIN (sim));
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) filter, f);
    gfs_traverse_and_bc (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) filtered, f,
			 f->v, f->v);
    gts_object_destroy (GTS_OBJECT (f->tmp));
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_filter_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_event_filter_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_filter_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if ((GFS_EVENT_FILTER (*o)->v = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  GFS_EVENT_FILTER (*o)->scale = gfs_read_constant (fp, domain);
  if (fp->type == GTS_ERROR)
    return;
  if (GFS_EVENT_FILTER (*o)->scale <= 0.) {
    gts_file_error (fp, "time scale must be strictly positive");
    return;
  }
}

static void gfs_event_filter_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_event_filter_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_filter_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %s %g", GFS_EVENT_FILTER (o)->v->name, GFS_EVENT_FILTER (o)->scale);
}

static void gfs_event_filter_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_event_filter_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_event_filter_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_filter_write;
}

GfsEventClass * gfs_event_filter_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_filter_info = {
      "GfsEventFilter",
      sizeof (GfsEventFilter),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_filter_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_filter_info);
  }

  return klass;
}

/** \endobject{GfsEventFilter} */

/**
 * A list of events.
 * \beginobject{GfsEventList}
 */

static gboolean gfs_event_list_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_list_class ())->parent_class)->event) 
      (event, sim)) {
    gts_container_foreach (GTS_CONTAINER (GFS_EVENT_LIST (event)->list), 
			   (GtsFunc) gfs_event_do, sim);
    return TRUE;
  }
  return FALSE;
}

static void copy_event (GfsEvent * event, GfsEvent * list)
{
  gfs_event_set (event, list->start, list->end, list->step, list->istart, list->iend, list->istep);
}

static void gfs_event_list_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_event_list_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsEventList * l = GFS_EVENT_LIST (*o);
  if (fp->type == GTS_STRING) {
    l->klass = gfs_object_class_from_name (fp->token->str);
    if (l->klass == NULL) {
      gts_file_error (fp, "unknown class `%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  while (fp->type == '\n') gts_file_next_token (fp);

  GtsObjectClass * default_class = l->klass;
  GfsSimulation * sim = gfs_object_simulation (*o);
  while (fp->type != '}') {
    GtsObjectClass * klass = default_class;
    if (klass == NULL) {
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a keyword");
	return;
      }
      klass = gfs_object_class_from_name (fp->token->str);
      if (klass == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	return;
      }
      if (!gts_object_class_is_from_class (klass, gfs_event_class ())) {
	gts_file_error (fp, "'%s' is not a GfsEvent", fp->token->str);
	return;
      }
    }
    
    GtsObject * object = gts_object_new (klass);
    gfs_object_simulation_set (object, sim);
    (* klass->read) (&object, fp);
    if (fp->type == GTS_ERROR) {
      gts_object_destroy (object);
      return;
    }
    while (fp->type == '\n') gts_file_next_token (fp);

    gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
  }
  fp->scope_max--;
  gts_file_next_token (fp);
  
  l->list->items = g_slist_reverse (l->list->items);
  gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc) copy_event, l);
}

static void gfs_event_list_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_event_list_class ())->parent_class->write) (o, fp);
  GfsEventList * l = GFS_EVENT_LIST (o);
  if (l->klass)
    fprintf (fp, " %s", l->klass->info.name);
  fputs (" {\n", fp);
  GSList * i = l->list->items;
  while (i) {
    fputs ("    ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    fputc ('\n', fp);
    i = i->next;
  }
  fputc ('}', fp);
}

static void gfs_event_list_destroy (GtsObject * o)
{
  GfsEventList * l = GFS_EVENT_LIST (o);
  gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (l->list));

  (* GTS_OBJECT_CLASS (gfs_event_list_class ())->parent_class->destroy) (o);
}

static void gfs_event_list_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_event_list_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_event_list_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_list_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_list_destroy;
}

static void gfs_event_list_init (GfsEventList * l)
{
  l->list = 
    GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ())));
}

GfsEventClass * gfs_event_list_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_list_info = {
      "GfsEventList",
      sizeof (GfsEventList),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_list_class_init,
      (GtsObjectInitFunc) gfs_event_list_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_list_info);
  }

  return klass;
}

/** \endobject{GfsEventList} */
