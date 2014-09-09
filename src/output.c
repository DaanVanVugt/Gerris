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
 * \brief Various outputs.
 */

#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "output.h"
#include "graphic.h"
#include "adaptive.h"
#include "solid.h"
#include "ocean.h"
#include "unstructured.h"
#include "init.h"

/**
 * Writing simulation data.
 * \beginobject{GfsOutput}
 */

static void output_free (GfsOutput * output)
{
  if (output->format)
    g_free (output->format);
  output->format = NULL;
  gfs_format_destroy (output->formats);
  output->formats = NULL;
}

static void gfs_output_destroy (GtsObject * object)
{
  GfsOutput * output = GFS_OUTPUT (object);

  if (output->file)
    gfs_output_file_close (output->file);
  output_free (output);

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->destroy) 
    (object);
}

static gboolean gfs_output_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsOutput * output = GFS_OUTPUT (event);
    gchar * fname;

    if (!output->parallel && GFS_DOMAIN (sim)->pid > 0) {
      if (output->file)
	output->first_call = FALSE;
      else
	gfs_output_mute (output);
      return (output->file != NULL);
    }

    if (!output->dynamic) {
      if (output->file)
	output->first_call = FALSE;
      else {
	if (output->format[0] == '{') { /* script */
	  guint len = strlen (output->format);
	  g_assert (output->format[len - 1] == '}');
	  output->format[len - 1] = '\0';
	  FILE * fp = gfs_popen (sim, &output->format[1], "w");
	  if (fp == NULL) {
	    g_warning ("GfsOutput cannot start script");
	    return TRUE;
	  }
	  output->file = gfs_output_file_new (fp);
	  output->file->is_pipe = TRUE;
	  output->format[len - 1] = '}';
	}
	else { /* standard file */
	  fname = gfs_format_string (output->formats,
				 GFS_DOMAIN (sim)->pid,
				 sim->time.i,
				 sim->time.t);
	  output->file = gfs_output_file_open (fname, 
					       sim->time.i > 0 && gfs_event_is_repetitive (event) ? 
					       "a" : "w");
	  if (output->file == NULL)
	    g_warning ("could not open file `%s'", fname);
	  g_free (fname);
	}
      }
      return (output->file != NULL);
    }

    if (output->file)
      gfs_output_file_close (output->file);
    fname = gfs_format_string (output->formats, 
			   GFS_DOMAIN (sim)->pid,
			   sim->time.i,
			   sim->time.t);
    output->file = gfs_output_file_open (fname, "w");
    if (output->file == NULL)
      g_warning ("could not open file `%s'", fname);
    g_free (fname);
    return (output->file != NULL);
  }
  return FALSE;
}

static void gfs_output_post_event (GfsEvent * event, GfsSimulation * sim)
{
  GfsOutput * output = GFS_OUTPUT (event);
  if (output->file)
    fflush (output->file->fp);
}

static void gfs_output_write (GtsObject * o, FILE * fp)
{
  GfsOutput * output = GFS_OUTPUT (o);

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->write) (o, fp);

  if (output->format)
    fprintf (fp, " %s", output->format);
}

static void gfs_output_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutput * output;

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT (*o);
  if (output->file)
    gfs_output_file_close (output->file);
  output->file = NULL;
  if (output->format)
    g_free (output->format);
  output->format = NULL;
  output->dynamic = FALSE;
  output->first_call = TRUE;

  if (fp->type == '{') {
    gchar * script = gfs_file_statement (fp);
    if (script == NULL)
      return;
    output->format = g_strconcat ("{", script, "}", NULL);
    g_free (script);
    gts_file_next_token (fp);
  }
  else if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (format)");
    return;
  }
  else {
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (output));

    output->format = g_strdup (fp->token->str);
    gts_file_next_token (fp);
    
    if (!strcmp (output->format, "stderr") || !strcmp (output->format, "stdout")) {
      if (domain->pid > 0)
	gfs_output_mute (output);
      else {
	g_assert (!output->file);
	output->file = gfs_output_file_open (output->format, "w");
      }
      return;
    }
    
    output->formats = gfs_format_new (output->format, fp, &output->dynamic, &output->parallel);
    if (fp->type == GTS_ERROR) {
      output_free (output);
      return;
    }

    if (output->parallel || domain->pid <= 0) {
      gchar * fname = gfs_format_string (output->formats, domain->pid, 0, 0.);
      gchar * fnamebak = g_strconcat (fname, "~", NULL);
      g_free (fname);
      FILE * fptr = fopen (fnamebak, "w");
      if (fptr == NULL) {
	gts_file_error (fp, "cannot open file specified by format `%s'\n  %s",
			output->format, strerror (errno));
	g_free (fnamebak);
	output_free (output);
	return;
      }
      fclose (fptr);
      remove (fnamebak);
      g_free (fnamebak);
    }
  }
}

static void gfs_output_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_event;
  GFS_EVENT_CLASS (klass)->post_event = gfs_output_post_event;

  GTS_OBJECT_CLASS (klass)->write = gfs_output_write;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_destroy;
}

static void gfs_output_init (GfsOutput * object)
{
  object->file = NULL;
  object->format = NULL;
  object->formats = NULL;
  object->dynamic = FALSE;
  object->parallel = FALSE;
  object->first_call = TRUE;
}

GfsOutputClass * gfs_output_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_info = {
      "GfsOutput",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_class_init,
      (GtsObjectInitFunc) gfs_output_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_output_info);
  }

  return klass;
}

/**
 * gfs_output_mute:
 * @output: a #GfsOutput.
 *
 * "Mutes" the output defined by @output, the event associated with
 * @output still takes place but the output itself is redirected to
 * /dev/null.
 */
void gfs_output_mute (GfsOutput * output)
{
  g_return_if_fail (output != NULL);

  output->dynamic = FALSE;
  if (output->file)
    gfs_output_file_close (output->file);
  output->file = gfs_output_file_open ("/dev/null", "w");
}

static GHashTable * gfs_output_files = NULL;

/**
 * gfs_output_file_new:
 * @fp: a file pointer.
 *
 * Returns: a new #GfsOutputFile for @fp.
 */
GfsOutputFile * gfs_output_file_new (FILE * fp)
{
  GfsOutputFile * file = g_malloc (sizeof (GfsOutputFile));
  file->refcount = 1;
  file->name = NULL;
  file->fp = fp;
  file->is_pipe = FALSE;
  return file;
}

/**
 * gfs_output_file_open:
 * @name: the name of the file to open.
 * @mode: the fopen mode.
 *
 * Checks whether @name has already been opened. If it has, its
 * reference count is incremented and the corresponding #GfsOutputFile
 * is returned. If it has not, it is created and opened for writing.
 *
 * Returns: the #GfsOutputFile of file @name.  
 */
GfsOutputFile * gfs_output_file_open (const gchar * name, const gchar * mode)
{
  GfsOutputFile * file;
  FILE * fp;

  g_return_val_if_fail (name != NULL, NULL);

  if (gfs_output_files == NULL) {
    gfs_output_files = g_hash_table_new (g_str_hash, g_str_equal);
    file = g_malloc (sizeof (GfsOutputFile));
    file->refcount = 2;
    file->name = g_strdup ("stderr");
    file->fp = stderr;
    g_hash_table_insert (gfs_output_files, file->name, file);
    file = g_malloc (sizeof (GfsOutputFile));
    file->refcount = 2;
    file->name = g_strdup ("stdout");
    file->fp = stdout;
    g_hash_table_insert (gfs_output_files, file->name, file);
  }

  if ((file = g_hash_table_lookup (gfs_output_files, name))) {
    file->refcount++;
    return file;
  }

  fp = fopen (name, mode);
  if (fp == NULL)
    return NULL;

  file = gfs_output_file_new (fp);
  file->name = g_strdup (name);
  g_hash_table_insert (gfs_output_files, file->name, file);

  return file;  
}

/**
 * gfs_output_file_close:
 * @file: a #GfsOutputFile.
 * 
 * Decreases the reference count of @file. If it reaches zero the file
 * corresponding to @file is closed and @file is freed.
 */
void gfs_output_file_close (GfsOutputFile * file)
{
  g_return_if_fail (file);

  file->refcount--;
  if (file->refcount == 0) {
    if (file->name)
      g_hash_table_remove (gfs_output_files, file->name);
    if (file->is_pipe)
      pclose (file->fp);
    else
      fclose (file->fp);
    g_free (file->name);
    g_free (file);
  }
}

/** \endobject{GfsOutput} */

/**
 * Model time, timestep and computing times.
 * \beginobject{GfsOutputTime}
 */

static gboolean time_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    gdouble cpu = gfs_clock_elapsed (domain->timer);
#ifdef HAVE_MPI
    if (domain->pid >= 0) {
      gfs_all_reduce (domain, cpu, MPI_DOUBLE, MPI_SUM);
      int size;
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      cpu /= size;
    }
#endif
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "step: %7u t: %15.8f dt: %13.6e cpu: %15.8f real: %15.8f\n",
	     sim->time.i, sim->time.t, 
	     sim->advection_params.dt,
	     cpu,
	     g_timer_elapsed (domain->clock, NULL));
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_time_class_init (GfsEventClass * klass)
{
  klass->event = time_event;
}

GfsOutputClass * gfs_output_time_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_time_info = {
      "GfsOutputTime",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_time_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_time_info);
  }

  return klass;
}

/** \endobject{GfsOutputTime} */

/**
 *
 * \beginobject{GfsOutputProgress}
 */

static gboolean progress_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    gdouble idone = sim->time.i/(gdouble) sim->time.iend;
    gdouble tdone = sim->time.t/sim->time.end;

    if (idone > tdone) tdone = idone;
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "\r%3.0f%% done",
	     100.*tdone);
    if (tdone > 0.) {
      gdouble remaining = GFS_DOMAIN (sim)->timestep.sum*(1. - tdone)/tdone;
      gdouble hours = floor (remaining/3600.);
      gdouble mins = floor ((remaining - 3600.*hours)/60.);
      gdouble secs = floor (remaining - 3600.*hours - 60.*mins);
      fprintf (GFS_OUTPUT (event)->file->fp,
	       ", %02.0f:%02.0f:%02.0f remaining ",
	       hours, mins, secs);
    }
    if (tdone == 1.)
      fputc ('\n', GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_progress_class_init (GfsEventClass * klass)
{
  klass->event = progress_event;
}

GfsOutputClass * gfs_output_progress_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_progress_info = {
      "GfsOutputProgress",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_progress_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_progress_info);
  }

  return klass;
}

/** \endobject{GfsOutputProgress} */

/**
 * 
 * \beginobject{GfsOutputProjectionStats}
 */

static gboolean projection_stats_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    if (sim->projection_params.niter > 0) {
      fprintf (fp, "MAC projection        before     after       rate\n");
      gfs_multilevel_params_stats_write (&sim->projection_params, fp);
    }
    fprintf (fp, "Approximate projection\n");
    gfs_multilevel_params_stats_write (&sim->approx_projection_params, fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_projection_stats_class_init (GfsEventClass * klass)
{
  klass->event = projection_stats_event;
}

GfsOutputClass * gfs_output_projection_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_projection_stats_info = {
      "GfsOutputProjectionStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_projection_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_projection_stats_info);
  }

  return klass;
}

/** \endobject{GfsOutputProjectionStats} */

/**
 *
 * \beginobject{GfsOutputDiffusionStats}
 */

static gboolean diffusion_stats_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    GSList * l = NULL, * i;
    
    i = GFS_DOMAIN (sim)->variables;
    while (i) {
      GfsVariable * v = i->data;

      if (v->sources) {
	GSList * j = GTS_SLIST_CONTAINER (v->sources)->items;
    
	while (j) {
	  GtsObject * o = j->data;
      
	  if (GFS_IS_SOURCE_DIFFUSION (o) && GFS_SOURCE_DIFFUSION (o)->D->par.niter > 0 &&
	      !g_slist_find (l, o)) {
	    l = g_slist_prepend (l, o);
	    fprintf (fp, "%s diffusion\n", v->name);
	    gfs_multilevel_params_stats_write (&GFS_SOURCE_DIFFUSION (o)->D->par, fp);
	  }
	  j = j->next;
	}
      }
      i = i->next;
    }
    g_slist_free (l);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_diffusion_stats_class_init (GfsEventClass * klass)
{
  klass->event = diffusion_stats_event;
}

GfsOutputClass * gfs_output_diffusion_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_diffusion_stats_info = {
      "GfsOutputDiffusionStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_diffusion_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_diffusion_stats_info);
  }

  return klass;
}

/** \endobject{GfsOutputDiffusionStats} */

/**
 * 
 * \beginobject{GfsOutputSolidStats}
 */

static gboolean gfs_output_solid_stats_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_solid_stats_class ())->parent_class)->event)
      (event, sim)) {
    GtsRange stats = gfs_domain_stats_solid (GFS_DOMAIN (sim));
    GtsRange ma, mn;

    gfs_domain_stats_merged (GFS_DOMAIN (sim), &ma, &mn);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "Solid volume fraction\n"
	     "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n"
	     "Total merged solid volume fraction\n"
	     "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n"
	     "Number of cells merged per merged cell\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n"
	     "Number of \"thin\" cells removed: %10d\n",
	     stats.min, stats.mean, stats.stddev, stats.max, stats.n,
	     ma.min, ma.mean, ma.stddev, ma.max, ma.n,
	     mn.min, mn.mean, mn.stddev, mn.max, mn.n,
	     sim->thin);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_solid_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_solid_stats_event;
}

GfsOutputClass * gfs_output_solid_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_solid_stats_info = {
      "GfsOutputSolidStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_solid_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_solid_stats_info);
  }

  return klass;
}

/** \endobject{GfsOutputSolidStats} */

/**
 * Information about the mesh adaptation.
 * \beginobject{GfsOutputAdaptStats}
 */

static gboolean gfs_output_adapt_stats_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_adapt_stats_class ())->parent_class)->event)
      (event, sim)) {
    gfs_adapt_stats_update (&sim->adapts_stats);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "Adaptive mesh refinement statistics\n"
	     "  Cells removed: %10d\n"
	     "  Cells created: %10d\n"
	     "  Number of cells\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n",
	     sim->adapts_stats.removed,
	     sim->adapts_stats.created,
	     sim->adapts_stats.ncells.min,
	     sim->adapts_stats.ncells.mean,
	     sim->adapts_stats.ncells.stddev,
	     sim->adapts_stats.ncells.max,
	     sim->adapts_stats.ncells.n);
    if (sim->adapts_stats.cmax.n > 0)
      fprintf (GFS_OUTPUT (event)->file->fp,
	       "  Maximum cost\n"
	       "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n",
	       sim->adapts_stats.cmax.min,
	       sim->adapts_stats.cmax.mean,
	       sim->adapts_stats.cmax.stddev,
	       sim->adapts_stats.cmax.max,
	       sim->adapts_stats.cmax.n);
    gfs_adapt_stats_init (&sim->adapts_stats);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_adapt_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_adapt_stats_event;
}

GfsOutputClass * gfs_output_adapt_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_adapt_stats_info = {
      "GfsOutputAdaptStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_adapt_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_adapt_stats_info);
  }

  return klass;
}

/** \endobject{GfsOutputAdaptStats} */

/**
 *
 * \beginobject{GfsOutputTiming}
 */

typedef struct {
  GfsTimer * t;
  gchar * name;
} Timer;

static int compare_timer (const void * a, const void * b)
{
  Timer * t1 = (Timer *) a;
  Timer * t2 = (Timer *) b;
  return (t1->t->r.sum < t2->t->r.sum) ? -1 : 1 ;
}

static void get_timer (gchar * name, GfsTimer * t, gpointer * data)
{
  Timer * timing = data[0];
  gint * count = data[1];
  timing[*count].t = t;
  timing[*count].name = name;
  (*count)++;
}

static void print_timing (GHashTable * timers, GfsDomain * domain, FILE * fp)
{
  Timer * timing = g_malloc (sizeof (Timer)*g_hash_table_size (timers));
  gint count = 0;
  gpointer data[2];

  data[0] = timing;
  data[1] = &count;  
  g_hash_table_foreach (domain->timers, (GHFunc) get_timer, data);
  qsort (timing, count, sizeof (Timer), compare_timer);
  while (--count >= 0)
    if (timing[count].t->r.sum > 0.)
      fprintf (fp, 
	       "  %s:\n"
	       "      min: %9.3f avg: %9.3f (%4.1f%%) | %7.3f max: %9.3f\n",
	       timing[count].name,
	       timing[count].t->r.min,
	       timing[count].t->r.mean, 
	       domain->timestep.sum > 0. ? 100.*timing[count].t->r.sum/domain->timestep.sum : 0.,
	       timing[count].t->r.stddev, 
	       timing[count].t->r.max);
  g_free (timing);
}

static gboolean timing_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    
    if (domain->timestep.mean > 0.) {
      fprintf (fp,
	       "Timing summary: %u timesteps %.0f node.timestep/s\n"
	       "  timestep:\n"
	       "      min: %9.3f avg: %9.3f         | %7.3f max: %9.3f\n"
               "  domain size:\n"
	       "      min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n"
	       "  maximum number of variables: %d\n",
	       domain->timestep.n,
	       domain->size.mean/domain->timestep.mean,
	       domain->timestep.min,
	       domain->timestep.mean,
	       domain->timestep.stddev, 
	       domain->timestep.max,
	       domain->size.min,
	       domain->size.mean,
	       domain->size.stddev, 
	       domain->size.max,
	       gfs_domain_variables_number (domain));
      print_timing (domain->timers, domain, fp);
      if (domain->mpi_messages.n > 0)
	fprintf (fp,
		 "Message passing summary\n"
		 "  n: %10d size: %10.0f bytes\n",
		 domain->mpi_messages.n,
		 domain->mpi_messages.sum);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_timing_class_init (GfsEventClass * klass)
{
  klass->event = timing_event;
}

GfsOutputClass * gfs_output_timing_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_timing_info = {
      "GfsOutputTiming",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_timing_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_timing_info);
  }

  return klass;
}

/** \endobject{GfsOutputTiming} */

/**
 * Writing simulation size statistics.
 * \beginobject{GfsOutputBalance}
 */

static gboolean gfs_output_balance_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_balance_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    GtsRange size, boundary, mpiwait;
    
    gfs_domain_stats_balance (domain, &size, &boundary, &mpiwait);
    fprintf (fp, 
	     "Balance summary: %u PE\n"
	     "  domain   min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n",
	     size.n,
	     size.min, size.mean, size.stddev, size.max);
    if (boundary.max > 0.)
      fprintf (fp, 
	       "  boundary min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n",
	       boundary.min, boundary.mean, boundary.stddev, boundary.max);
    if (mpiwait.max > 0.)
      fprintf (fp,
	       "  average timestep MPI wait time:\n"
	       "      min: %9.3f avg: %9.3f         | %7.3f max: %9.3f\n",
	       mpiwait.min, mpiwait.mean, mpiwait.stddev, mpiwait.max);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_balance_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_balance_event;
}

GfsOutputClass * gfs_output_balance_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_balance_info = {
      "GfsOutputBalance",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_balance_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_balance_info);
  }

  return klass;
}

/** \endobject{GfsOutputBalance} */

/**
 * orces and moments on the embedded solid boundaries.
 * \beginobject{GfsOutputSolidForce}
 */

static void gfs_output_solid_force_destroy (GtsObject * object)
{
  if (GFS_OUTPUT_SOLID_FORCE (object)->weight)
    gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_SOLID_FORCE (object)->weight));

  (* GTS_OBJECT_CLASS (gfs_output_solid_force_class ())->parent_class->destroy) (object);
}

static void gfs_output_solid_force_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputSolidForce * l = GFS_OUTPUT_SOLID_FORCE (*o);

  (* GTS_OBJECT_CLASS (gfs_output_solid_force_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '\n') {
    if (!l->weight)
      l->weight = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (l->weight, gfs_object_simulation (l), fp);
  }
}

static void gfs_output_solid_force_write (GtsObject * o, FILE * fp)
{
  GfsOutputSolidForce * l = GFS_OUTPUT_SOLID_FORCE (o);
  (* GTS_OBJECT_CLASS (gfs_output_solid_force_class ())->parent_class->write) (o, fp);
  if (l->weight)
    gfs_function_write (l->weight, fp);
}

static gboolean gfs_output_solid_force_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_solid_force_class ())->parent_class)->event)
      (event, sim) &&
      sim->advection_params.dt > 0.) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    FttVector pf, vf, pm, vm;
    gdouble L = sim->physical_params.L, Ln = pow (L, 3. + FTT_DIMENSION - 2.);

    if (GFS_OUTPUT (event)->first_call)
      fputs ("# 1: T (2,3,4): Pressure force (5,6,7): Viscous force "
	     "(8,9,10): Pressure moment (11,12,13): Viscous moment\n", fp);
    
    gfs_domain_solid_force (domain, &pf, &vf, &pm, &vm, GFS_OUTPUT_SOLID_FORCE (event)->weight);
    fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	     sim->time.t,
	     pf.x*Ln, pf.y*Ln, pf.z*Ln,
	     vf.x*Ln, vf.y*Ln, vf.z*Ln,
	     pm.x*Ln*L, pm.y*Ln*L, pm.z*Ln*L,
	     vm.x*Ln*L, vm.y*Ln*L, vm.z*Ln*L);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_solid_force_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_output_solid_force_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_solid_force_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_solid_force_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_output_solid_force_event;
}

GfsOutputClass * gfs_output_solid_force_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_solid_force_info = {
      "GfsOutputSolidForce",
      sizeof (GfsOutputSolidForce),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_solid_force_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_solid_force_info);
  }

  return klass;
}

/** \endobject{GfsOutputSolidForce} */

/**
 * Writing the values of variables at specified locations.
 * \beginobject{GfsOutputLocation}
 */

static gchar default_precision[] = "%g";

static void gfs_output_location_destroy (GtsObject * object)
{
  GfsOutputLocation * l = GFS_OUTPUT_LOCATION (object);
  g_array_free (l->p, TRUE);
  g_free (l->label);
  if (l->precision != default_precision)
    g_free (l->precision);

  (* GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->destroy) (object);
}

static gboolean vector_read (GtsFile * fp, FttVector * p)
{
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return FALSE;
  }
  p->x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return FALSE;
  }
  p->y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return FALSE;
  }
  p->z = atof (fp->token->str);
  gts_file_next_token (fp);
  return TRUE;
}

static void gfs_output_location_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputLocation * l = GFS_OUTPUT_LOCATION (*o);

  if (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_STRING) {
    FILE * fptr = fopen (fp->token->str, "r");
    GtsFile * fp1;

    if (fptr == NULL) {
      gts_file_error (fp, "cannot open file `%s'", fp->token->str);
      return;
    }
    fp1 = gts_file_new (fptr);
    while (fp1->type != GTS_NONE) {
      FttVector p;
      if (!vector_read (fp1, &p)) {
	gts_file_error (fp, "%s:%d:%d: %s", fp->token->str, fp1->line, fp1->pos, fp1->error);
	return;
      }
      g_array_append_val (l->p, p);
      while (fp1->type == '\n')
	gts_file_next_token (fp1);
    }
    gts_file_destroy (fp1);
    fclose (fptr);
    gts_file_next_token (fp);
  }
  else if (fp->type == '{') {
    fp->scope_max++;
    do
      gts_file_next_token (fp);
    while (fp->type == '\n');
    while (fp->type != GTS_NONE && fp->type != '}') {
      FttVector p;
      if (!vector_read (fp, &p))
	return;
      g_array_append_val (l->p, p);
      while (fp->type == '\n')
	gts_file_next_token (fp);
    }
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
    gts_file_next_token (fp);
  }
  else {
    FttVector p;
    if (!vector_read (fp, &p))
      return;
    g_array_append_val (l->p, p);
  }

  if (fp->type == '{') {
    gchar * label = NULL, * precision = NULL;
    GtsFileVariable var[] = {
      {GTS_STRING, "label", TRUE, &label},
      {GTS_STRING, "precision", TRUE, &precision},
      {GTS_INT,    "interpolate", TRUE, &l->interpolate},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR) {
      g_free (label);
      g_free (precision);
      return;
    }

    if (precision != NULL) {
      if (l->precision != default_precision)
	g_free (l->precision);
      l->precision = precision;
    }

    if (label != NULL) {
      g_free (l->label);
      l->label = label;
    }
  }
}

static void gfs_output_location_write (GtsObject * o, FILE * fp)
{
  GfsOutputLocation * l = GFS_OUTPUT_LOCATION (o);
  guint i;

  (* GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->write) (o, fp);

  fputs (" {\n", fp);
  gchar * format = g_strdup_printf ("%s %s %s\n", l->precision, l->precision, l->precision);
  for (i = 0; i < l->p->len; i++) {
    FttVector p = g_array_index (l->p, FttVector, i);
    fprintf (fp, format, p.x, p.y, p.z);
  }
  g_free (format);
  fputc ('}', fp);

  if (l->precision != default_precision || l->label) {
    fputs (" {\n", fp);
    if (l->precision != default_precision)
      fprintf (fp, "  precision = %s\n", l->precision);
    if (l->label)
      fprintf (fp, "  label = \"%s\"\n", l->label);
    if (!l->interpolate)
      fputs ("  interpolate = 0\n", fp);
    fputc ('}', fp);
  }
}

static gboolean gfs_output_location_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputLocation * location = GFS_OUTPUT_LOCATION (event);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    GfsUnionFile uf;
    FILE * fpp = ((domain->pid < 0 || GFS_OUTPUT (event)->parallel) ? fp:
		  gfs_union_open (GFS_OUTPUT (event)->file->fp, domain->pid, &uf));
    guint i;

    if (GFS_OUTPUT (event)->first_call) {
      GSList * i = domain->variables;
      guint nv = 5;

      fputs ("# 1:t 2:x 3:y 4:z", fp);
      while (i) {
	if (GFS_VARIABLE (i->data)->name)
	  fprintf (fp, " %d:%s", nv++, GFS_VARIABLE (i->data)->name);
	i = i->next;
      }
      fputc ('\n', fp);
    }
    gchar * pformat = g_strdup_printf ("%s %s %s %s", 
				       location->precision, location->precision, 
				       location->precision, location->precision);
    gchar * vformat = g_strdup_printf (" %s", location->precision);
    for (i = 0; i < location->p->len; i++) {
      FttVector p = g_array_index (location->p, FttVector, i), pm = p;
      gfs_simulation_map (sim, &pm);
      FttCell * cell = gfs_domain_locate (domain, pm, -1, NULL);
      
      if (cell != NULL) {
	GSList * i = domain->variables;
	
	fprintf (fpp, pformat, sim->time.t, p.x, p.y, p.z);
	while (i) {
	  GfsVariable * v = i->data;
	  if (v->name)
	    fprintf (fpp, vformat, gfs_dimensional_value (v, 
							  location->interpolate ? 
							  gfs_interpolate (cell, pm, v) :
							  GFS_VALUE (cell, v)));
	  i = i->next;
	}
	fputc ('\n', fpp);
      }
    }

    g_free (pformat);
    g_free (vformat);
    fflush (fp);
    if (!(domain->pid < 0 || GFS_OUTPUT (event)->parallel))
      gfs_union_close (GFS_OUTPUT (event)->file->fp, domain->pid, &uf);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_location_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_location_event;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_location_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_location_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_location_write;
}

static void gfs_output_location_init (GfsOutputLocation * object)
{
  object->p = g_array_new (FALSE, FALSE, sizeof (FttVector));
  object->precision = default_precision;
  object->interpolate = TRUE;
}

GfsOutputClass * gfs_output_location_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_location_info = {
      "GfsOutputLocation",
      sizeof (GfsOutputLocation),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_location_class_init,
      (GtsObjectInitFunc) gfs_output_location_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_location_info);
  }

  return klass;
}

/** \endobject{GfsOutputLocation} */

/**
 * Tracking Lagrangian particles.
 * \beginobject{GfsOutputParticle}
 */

static gboolean gfs_output_particle_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  GfsOutputLocation * location = GFS_OUTPUT_LOCATION (event);
  gboolean ret = FALSE;
  guint i;

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class)->event)
      (event,sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    for (i = 0; i < location->p->len; i++) {
      FttVector p = g_array_index (location->p, FttVector, i);
      fprintf (fp, "%d %g %g %g %g\n", i, sim->time.t, p.x, p.y, p.z);
    }
    ret = TRUE;
  }
  
  for (i = 0; i < location->p->len; i++) {
    FttVector p = g_array_index (location->p, FttVector, i);
    gfs_simulation_map (sim, &p);
    gfs_domain_advect_point (GFS_DOMAIN (sim), &p, sim->advection_params.dt);
    gfs_simulation_map_inverse (sim, &p);
    g_array_index (location->p, FttVector, i) = p;
  }

  return ret;
}

static void gfs_output_particle_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_particle_event;
}

GfsOutputClass * gfs_output_particle_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_particle_info = {
      "GfsOutputParticle",
      sizeof (GfsOutputLocation),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_particle_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_location_class ()),
				  &gfs_output_particle_info);
  }

  return klass;
}

/** \endobject{GfsOutputParticle} */

/**
 * Writing the whole simulation.
 * \beginobject{GfsOutputSimulation}
 */

static void output_simulation_destroy (GtsObject * object)
{
  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (object);

  g_slist_free (output->var);
  if (output->precision != default_precision)
    g_free (output->precision);

  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->destroy) (object);
}

static void write_text (FttCell * cell, gpointer * data)
{
  GfsOutputSimulation * output = data[0];  
  FILE * fp = data[1];
  GSList * i = GFS_DOMAIN (gfs_object_simulation (output))->variables_io;
  FttVector p;

  gfs_cell_cm (cell, &p);
  gfs_simulation_map_inverse (gfs_object_simulation (output), &p);
  gchar * format = g_strdup_printf ("%s %s %s", 
				    output->precision, output->precision, output->precision);
  fprintf (fp, format, p.x, p.y, p.z);
  g_free (format);
  format = g_strdup_printf (" %s", output->precision);
  while (i) {
    if (GFS_VARIABLE (i->data)->name)
      fprintf (fp, format, gfs_dimensional_value (i->data, 
						  GFS_VALUE (cell, GFS_VARIABLE (i->data))));
    i = i->next;
  }
  g_free (format);
  fputc ('\n', fp);
}

static gboolean output_simulation_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (event);

    g_slist_free (domain->variables_io);
    if (output->var)
      domain->variables_io = output->var;
    else {
      GSList * i = domain->variables;
      domain->variables_io = NULL;
      while (i) {
	if (GFS_VARIABLE (i->data)->name)
	  domain->variables_io = g_slist_append (domain->variables_io, i->data);
	i = i->next;
      }
    }

    domain->binary =       output->binary;
    sim->output_solid   =  output->solid;
    switch (output->format) {

    case GFS:
      if (GFS_OUTPUT (output)->parallel)
	gfs_simulation_write (sim,
			      output->max_depth,
			      GFS_OUTPUT (event)->file->fp);
      else
	gfs_simulation_union_write (sim,
				    output->max_depth,
				    GFS_OUTPUT (event)->file->fp);
      break;

    case GFS_TEXT: {
      if (GFS_OUTPUT (output)->parallel || domain->pid <= 0) {
	FILE * fp = GFS_OUTPUT (event)->file->fp;
	GSList * i = domain->variables_io;
	guint nv = 4;

	fputs ("# 1:x 2:y 3:z", fp);
	while (i) {
	  g_assert (GFS_VARIABLE (i->data)->name);
	  fprintf (fp, " %d:%s", nv++, GFS_VARIABLE (i->data)->name);
	  i = i->next;
	}
	fprintf (fp, " %g\n", sim->time.t);
      }
      gpointer data[2];
      data[0] = output;
      if (GFS_OUTPUT (output)->parallel) {
	data[1] = GFS_OUTPUT (event)->file->fp;
	gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL|FTT_TRAVERSE_LEAFS,
				  output->max_depth,
				  (FttCellTraverseFunc) write_text, data);
      }
      else {
	GfsUnionFile uf;
	FILE * fpp = gfs_union_open (GFS_OUTPUT (event)->file->fp, domain->pid, &uf);
	data[1] = fpp;
	gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL|FTT_TRAVERSE_LEAFS,
				  output->max_depth,
				  (FttCellTraverseFunc) write_text, data);
	gfs_union_close (GFS_OUTPUT (event)->file->fp, domain->pid, &uf);
      }
      break;
    }

    case GFS_VTK: {
      gfs_domain_write_vtk (domain, output->max_depth, domain->variables_io, output->precision,
			    GFS_OUTPUT (event)->file->fp);
      break;
    }

    case GFS_TECPLOT: {
      gfs_domain_write_tecplot (domain, output->max_depth, domain->variables_io, output->precision,
				GFS_OUTPUT (event)->file->fp);
#if !FTT_2D 
      gfs_domain_write_tecplot_surface (domain, output->max_depth, domain->variables_io, 
					output->precision,
					GFS_OUTPUT (event)->file->fp);
#endif
      break;
    }

    default:
      g_assert_not_reached ();
    }
    if (!output->var)
      g_slist_free (domain->variables_io);
    domain->variables_io = NULL;
    domain->binary =       TRUE;
    sim->output_solid   =  TRUE;
    return TRUE;
  }
  return FALSE;
}

static void output_simulation_write (GtsObject * o, FILE * fp)
{
  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (o);
  GSList * i = output->var;

  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->write) (o, fp);

  fputs (" {", fp);
  if (output->max_depth != -1)
    fprintf (fp, " depth = %d", output->max_depth);
  if (i != NULL) {
    fprintf (fp, " variables = %s", GFS_VARIABLE (i->data)->name);
    i = i->next;
    while (i) {
      fprintf (fp, ",%s", GFS_VARIABLE (i->data)->name);
      i = i->next;
    }
  }
  if (!output->binary)
    fputs (" binary = 0", fp);
  if (!output->solid)
    fputs (" solid = 0", fp);
  switch (output->format) {
  case GFS_TEXT:    fputs (" format = text", fp);    break;
  case GFS_VTK:     fputs (" format = VTK", fp);     break;
  case GFS_TECPLOT: fputs (" format = Tecplot", fp); break;
  default: break;
  }
  if (output->precision != default_precision)
    fprintf (fp, " precision = %s", output->precision);
  fputs (" }", fp);
}

static void output_simulation_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (*o);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT,    "depth",     TRUE},
      {GTS_STRING, "variables", TRUE},
      {GTS_INT,    "binary",    TRUE},
      {GTS_INT,    "solid",     TRUE},
      {GTS_STRING, "format",    TRUE},
      {GTS_STRING, "precision", TRUE},
      {GTS_NONE}
    };
    gchar * variables = NULL, * format = NULL, * precision = NULL;

    var[0].data = &output->max_depth;
    var[1].data = &variables;
    var[2].data = &output->binary;
    var[3].data = &output->solid;
    var[4].data = &format;
    var[5].data = &precision;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR) {
      g_free (variables);
      g_free (format);
      g_free (precision);
      return;
    }

    if (variables != NULL) {
      gchar * error = NULL;
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (output));
      GSList * vars = gfs_variables_from_list (domain->variables, variables, &error);

      if (vars == NULL) {
	gts_file_variable_error (fp, var, "variables",
				 "unknown variable `%s'", error);
	g_free (variables);
	return;
      }
      g_slist_free (output->var);
      output->var = vars;
      g_free (variables);
    }

    if (format != NULL) {
      if (!strcmp (format, "gfs"))
	output->format = GFS;
      else if (!strcmp (format, "text"))
	output->format = GFS_TEXT;
      else if (!strcmp (format, "VTK"))
	output->format = GFS_VTK;
      else if (!strcmp (format, "Tecplot"))
	output->format = GFS_TECPLOT;
      else {
	gts_file_variable_error (fp, var, "format",
				 "unknown format `%s'", format);
	g_free (format);
	return;
      }
      g_free (format);
    }

    if (precision != NULL) {
      if (output->precision != default_precision)
	g_free (output->precision);
      output->precision = precision;
    }
  }
}

static void gfs_output_simulation_class_init (GfsEventClass * klass)
{
  klass->event = output_simulation_event;
  GTS_OBJECT_CLASS (klass)->destroy = output_simulation_destroy;
  GTS_OBJECT_CLASS (klass)->read = output_simulation_read;
  GTS_OBJECT_CLASS (klass)->write = output_simulation_write;
}

static void gfs_output_simulation_init (GfsOutputSimulation * object)
{
  object->max_depth = -1;
  object->var = NULL;
  object->binary = 1;
  object->solid = 1;
  object->format = GFS;
  object->precision = default_precision;
}

GfsOutputClass * gfs_output_simulation_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_simulation_info = {
      "GfsOutputSimulation",
      sizeof (GfsOutputSimulation),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_simulation_class_init,
      (GtsObjectInitFunc) gfs_output_simulation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_simulation_info);
  }

  return klass;
}

/** \endobject{GfsOutputSimulation} */

/**
 *
 * \beginobject{GfsOutputBoundaries}
 */

static gboolean output_boundaries_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    
    gfs_draw_refined_boundaries (domain, fp);
    gfs_draw_solid_boundaries (domain, fp);
    gfs_draw_boundary_conditions (domain, fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_boundaries_class_init (GfsEventClass * klass)
{
  klass->event = output_boundaries_event;
}

GfsOutputClass * gfs_output_boundaries_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_boundaries_info = {
      "GfsOutputBoundaries",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_boundaries_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_boundaries_info);
  }

  return klass;
}

/** \endobject{GfsOutputBoundaries} */

/**
 * Generic output of scalar fields.
 * \beginobject{GfsOutputScalar}
 */

static void gfs_output_scalar_destroy (GtsObject * o)
{
  GfsOutputScalar * output = GFS_OUTPUT_SCALAR (o);

  gts_object_destroy (GTS_OBJECT (output->f));
  g_free (output->name);
  if (output->condition)
    gts_object_destroy (GTS_OBJECT (output->condition));
  if (output->w)
    gts_object_destroy (GTS_OBJECT (output->w));
  if (output->format)
    g_free (output->format);
  
  (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->destroy) (o);
}

static void gfs_output_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputScalar * output;

  if (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT_SCALAR (*o);
  output->autoscale = TRUE;

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
    else if (!strcmp (fp->token->str, "v")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (output->f, gfs_object_simulation (*o), fp);
      output->name = gfs_function_description (output->f, TRUE);
    }
    else if (!strcmp (fp->token->str, "min")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      output->min = gfs_read_constant (fp, gfs_object_simulation (*o));
      if (fp->type == GTS_ERROR)
	return;
      if (output->min > output->max) {
	gts_file_error (fp, "min `%g' must be smaller than or equal to max `%g'", 
			output->min, output->max);
	return;
      }
      output->autoscale = FALSE;
    }
    else if (!strcmp (fp->token->str, "max")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      output->max = gfs_read_constant (fp, gfs_object_simulation (*o));
      if (fp->type == GTS_ERROR)
	return;
      if (output->max < output->min) {
	gts_file_error (fp, "max `%g' must be larger than or equal to min `%g'", 
			output->max, output->min);
	return;
      }
      output->autoscale = FALSE;
    }
    else if (!strcmp (fp->token->str, "maxlevel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_INT) {
	gts_file_error (fp, "expecting an integer (maxlevel)");
	return;
      }
      output->maxlevel = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "condition")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      if (!output->condition)
	output->condition = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (output->condition, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "w")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      if (!output->w)
	output->w = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (output->w, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "format")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a string");
	return;
      }
      if (fp->token->str[0] == '"') {
	output->format = g_strdup (&(fp->token->str[1]));
	output->format[strlen(output->format) - 1] = '\0';
      }
      else
	output->format = g_strdup (fp->token->str);
      gts_file_next_token (fp);
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

static void gfs_output_scalar_write (GtsObject * o, FILE * fp)
{
  GfsOutputScalar * output = GFS_OUTPUT_SCALAR (o);

  if (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->write) 
      (o, fp);

  fputs (" { v = ", fp);
  gfs_function_write (output->f, fp);
  if (output->maxlevel >= 0)
    fprintf (fp, " maxlevel = %d", output->maxlevel);
  if (output->condition) {
    fputs (" condition = ", fp);
    gfs_function_write (output->condition, fp);
  }
  if (output->w) {
    fputs (" w = ", fp);
    gfs_function_write (output->w, fp);
  }
  if (output->format)
    fprintf (fp, " format = %s", output->format);
  if (!output->autoscale)
    fprintf (fp, " min = %g max = %g }", output->min, output->max);
  else
    fputs (" }", fp);
}

static void update_v (FttCell * cell, GfsOutputScalar * output)
{
  GFS_VALUE (cell, output->v) = gfs_function_value (output->f, cell);
}

static gboolean cell_condition (FttCell * cell, gpointer condition)
{
  return gfs_function_value (condition, cell);
}

static void output_scalar_traverse (GfsOutputScalar * output, 
				    FttTraverseType order,
				    FttTraverseFlags flags,
				    gint max_depth,
				    FttCellTraverseFunc func,
				    gpointer data)
{
  if (output->condition)
    gfs_domain_cell_traverse_condition (GFS_DOMAIN (gfs_object_simulation (output)),
					order, flags, max_depth, 
					func, data,
					cell_condition, output->condition);
  else
    gfs_domain_cell_traverse (GFS_DOMAIN (gfs_object_simulation (output)),
			      order, flags, max_depth, 
			      func, data);
}

static gboolean gfs_output_scalar_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsDomain * domain = GFS_DOMAIN (sim);

    if (!(output->v = gfs_function_get_variable (output->f)) ||
	gfs_variable_is_dimensional (output->v)) {
      output->v = gfs_temporary_variable (domain);
      gfs_catch_floating_point_exceptions ();
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) update_v, output);
      gfs_restore_fpe_for_function (output->f);
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, output->v);
    }
    if (output->maxlevel >= 0) {
      gfs_domain_cell_traverse (domain,
				FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) output->v->fine_coarse,
				output->v);
      gfs_domain_bc (domain, FTT_TRAVERSE_NON_LEAFS, -1, output->v);
    }
    if (output->autoscale) {
      GtsRange stats = gfs_domain_stats_variable (domain, output->v, 
						  FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
						  output->maxlevel,
						  output->condition ? cell_condition : NULL,
						  output->condition);
      output->min = stats.min;
      output->max = stats.max;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_post_event (GfsEvent * event,
					  GfsSimulation * sim)
{
  GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
  if (output->v != gfs_function_get_variable (output->f)) {
    gts_object_destroy (GTS_OBJECT (output->v));
    output->v = NULL;
  }

  (* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class)->post_event) 
    (event, sim);
}

static void gfs_output_scalar_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_event;
  GFS_EVENT_CLASS (klass)->post_event = gfs_output_scalar_post_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_scalar_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_scalar_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_scalar_destroy;
}

static void gfs_output_scalar_init (GfsOutputScalar * object)
{
  object->f = gfs_function_new (gfs_function_class (), 0.);
  object->min = -G_MAXDOUBLE;
  object->max =  G_MAXDOUBLE;
  object->autoscale = TRUE;
  object->maxlevel = -1;
  object->condition = NULL;
}

GfsOutputClass * gfs_output_scalar_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_info = {
      "GfsOutputScalar",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_class_init,
      (GtsObjectInitFunc) gfs_output_scalar_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_scalar_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalar} */

/**
 * Computing the norms of a scalar field.
 * \beginobject{GfsOutputScalarNorm}
 */

static gboolean gfs_output_scalar_norm_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_norm_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsNorm norm = gfs_domain_norm_variable (GFS_DOMAIN (sim), 
					     output->v, NULL,
					     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
					     output->maxlevel,
					     output->condition ? cell_condition : NULL,
					     output->condition);
    fprintf (GFS_OUTPUT (event)->file->fp, 
	     "%s time: %g first: % 10.3e second: % 10.3e infty: % 10.3e\n",
	     output->name,
	     sim->time.t,
	     norm.first, norm.second, norm.infty);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_norm_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_norm_event;
}

GfsOutputClass * gfs_output_scalar_norm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_norm_info = {
      "GfsOutputScalarNorm",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_norm_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_norm_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalarNorm} */

/**
 * Computing simple statistics for a scalar field.
 * \beginobject{GfsOutputScalarStats}
 */

static gboolean gfs_output_scalar_stats_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_stats_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GtsRange stats = gfs_domain_stats_variable (GFS_DOMAIN (sim), 
						output->v,
						FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
						output->maxlevel,
						output->condition ? cell_condition : NULL,
						output->condition);
    gchar * format;
    if (output->format) {
      gchar * f = output->format;
      format = g_strdup_printf ("%%s time: %s min: %s avg: %s | %s max: %s\n",
				f, f, f, f, f);
    }
    else
      format = g_strdup ("%s time: %g min: %10.3e avg: %10.3e | %10.3e max: %10.3e\n");
    fprintf (GFS_OUTPUT (event)->file->fp, format,
	     output->name, sim->time.t,
	     stats.min, stats.mean, stats.stddev, stats.max);
    g_free (format);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_stats_event;
}

GfsOutputClass * gfs_output_scalar_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_stats_info = {
      "GfsOutputScalarStats",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_stats_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalarStats} */

/**
 * Computing the sum of a scalar field.
 * \beginobject{GfsOutputScalarSum}
 */

typedef struct {
  GfsVariable * v;
  GfsFunction * w;
  gdouble sum;
} SumData;

static void add (FttCell * cell, SumData * s)
{
  gdouble vol = s->w ? gfs_function_value (s->w, cell) : gfs_cell_volume (cell, s->v->domain);
  s->sum += vol*GFS_VALUE (cell, s->v);
}

static gboolean gfs_output_scalar_sum_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_sum_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    SumData s = { output->v, output->w, 0. };
    output_scalar_traverse (output,
			    FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			    output->maxlevel,
			    (FttCellTraverseFunc) add, &s);
    gfs_all_reduce (GFS_DOMAIN (sim), s.sum, MPI_DOUBLE, MPI_SUM);
    if (!output->w)
      s.sum *= pow (sim->physical_params.L, FTT_DIMENSION);
    gchar * format;
    if (output->format) {
      gchar * f = output->format;
      format = g_strdup_printf ("%%s time: %s sum: %s\n", f, f);
    }
    else
      format = g_strdup ("%s time: %g sum: % 15.6e\n");
    fprintf (GFS_OUTPUT (event)->file->fp, format,
	     output->name, sim->time.t, s.sum);
    g_free (format);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_sum_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_sum_event;
}

GfsOutputClass * gfs_output_scalar_sum_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_sum_info = {
      "GfsOutputScalarSum",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_sum_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_sum_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalarSum} */

/**
 *
 * \beginobject{GfsOutputScalarMaxima}
 */

static void gfs_output_scalar_maxima_destroy (GtsObject * o)
{
  guint i;

  for (i = 0; i < 4; i++)
    g_free (GFS_OUTPUT_SCALAR_MAXIMA (o)->m[i]);

  (* GTS_OBJECT_CLASS (gfs_output_scalar_maxima_class ())->parent_class->destroy) (o);
}

static void gfs_output_scalar_maxima_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputScalarMaxima * m;
  guint i;

  (* GTS_OBJECT_CLASS (gfs_output_scalar_maxima_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (N)");
    return;
  }
  m = GFS_OUTPUT_SCALAR_MAXIMA (*o);
  m->N = atoi (fp->token->str);
  gts_file_next_token (fp);

  for (i = 0; i < 4; i++)
    m->m[i] = g_malloc (sizeof (gdouble)*m->N);
}

static void gfs_output_scalar_maxima_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_scalar_maxima_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %d", GFS_OUTPUT_SCALAR_MAXIMA (o)->N);
}

static void maxima (FttCell * cell, GfsOutputScalarMaxima * m)
{
  guint i;

  for (i = 0; i < m->N; i++) {
    gdouble v = GFS_VALUE (cell, GFS_OUTPUT_SCALAR (m)->v);

    if (v > m->m[3][i]) {
      FttVector p;

      gfs_cell_cm (cell, &p);
      gfs_simulation_map_inverse (gfs_object_simulation (m), &p);
      m->m[0][i] = p.x; m->m[1][i] = p.y; m->m[2][i] = p.z;
      m->m[3][i] = v;
      return;
    }
  }
}

static gboolean gfs_output_scalar_maxima_event (GfsEvent * event, 
						GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_maxima_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsOutputScalarMaxima * m = GFS_OUTPUT_SCALAR_MAXIMA (event);
    guint i;

    for (i = 0; i < m->N; i++)
      m->m[3][i] = -G_MAXDOUBLE;

    output_scalar_traverse (output,
			    FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			    output->maxlevel,
			    (FttCellTraverseFunc) maxima, m);
    for (i = 0; i < m->N; i++)
      fprintf (GFS_OUTPUT (event)->file->fp, 
	       "%s time: %g #: %d x: %g y: %g z: %g value: %g\n", 
	       output->name, sim->time.t, i,
	       m->m[0][i], m->m[1][i], m->m[2][i],
	       m->m[3][i]);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_maxima_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_scalar_maxima_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_scalar_maxima_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_scalar_maxima_write;
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_maxima_event;
}

GfsOutputClass * gfs_output_scalar_maxima_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_maxima_info = {
      "GfsOutputScalarMaxima",
      sizeof (GfsOutputScalarMaxima),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_maxima_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_maxima_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalarMaxima} */

/**
 * Computing histograms.
 * \beginobject{GfsOutputScalarHistogram}
 */

static void gfs_output_scalar_histogram_destroy (GtsObject * o)
{
  GfsOutputScalarHistogram * output = GFS_OUTPUT_SCALAR_HISTOGRAM (o);

  g_free (output->x);
  g_free (output->w);
  if (output->wf)
    gts_object_destroy (GTS_OBJECT (output->wf));
  if (output->yf) {
    gts_object_destroy (GTS_OBJECT (output->yf));
    g_free (output->y);
  }

  (* GTS_OBJECT_CLASS (gfs_output_scalar_histogram_class ())->parent_class->destroy) (o);
}

static void gfs_output_scalar_histogram_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputScalarHistogram * output;

  (* GTS_OBJECT_CLASS (gfs_output_scalar_histogram_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT_SCALAR_HISTOGRAM (*o);
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
    else if (!strcmp (fp->token->str, "n")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_INT) {
	gts_file_error (fp, "expecting a number (n)");
	return;
      }
      output->n = atoi (fp->token->str);
      if (output->n <= 0) {
	gts_file_error (fp, "n `%d' must be strictly positive", output->n);
	return;
      }
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "w")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      output->wf = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (output->wf, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "y")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      output->yf = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (output->yf, gfs_object_simulation (*o), fp);
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

  output->x = g_malloc0 (output->n*sizeof (gdouble));
  output->w = g_malloc0 (output->n*sizeof (gdouble));
  if (output->yf)
    output->y = g_malloc0 (output->n*sizeof (gdouble));
}

static void gfs_output_scalar_histogram_write (GtsObject * o, FILE * fp)
{
  GfsOutputScalarHistogram * output = GFS_OUTPUT_SCALAR_HISTOGRAM (o);

  (* GTS_OBJECT_CLASS (gfs_output_scalar_histogram_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { n = %d", output->n);
  if (output->wf) {
    fputs (" w = ", fp);
    gfs_function_write (output->wf, fp);
  }
  if (output->yf) {
    fputs (" y = ", fp);
    gfs_function_write (output->yf, fp);
  }
  fputs (" }", fp);
}

static void update_histogram (FttCell * cell, GfsOutputScalar * h)
{
  GfsOutputScalarHistogram * hi = GFS_OUTPUT_SCALAR_HISTOGRAM (h);
  gdouble v = GFS_VALUE (cell, h->v);
  gint i = (v - h->min)/(h->max - h->min)*hi->n;

  if (i >= 0 && i < hi->n) {
    gdouble w = hi->dt;

    if (hi->wf)
      w *= gfs_function_value (hi->wf, cell);
    else
      w *= gfs_cell_volume (cell, h->v->domain);

    hi->W += w;
    hi->w[i] += w;
    hi->x[i] += v*w;
    if (hi->yf)
      hi->y[i] += w*gfs_function_value (hi->yf, cell);
  }
}

static gboolean gfs_output_scalar_histogram_event (GfsEvent * event,
						   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_histogram_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalarHistogram * h = GFS_OUTPUT_SCALAR_HISTOGRAM (event);

    if (gfs_event_is_repetitive (event))
      h->dt = h->last >= 0. ? sim->time.t - h->last : 0.;
    else
      h->dt = 1.;

    if (h->dt > 0.) {
      GfsOutput * output = GFS_OUTPUT (event);
      guint i;

      gfs_catch_floating_point_exceptions ();
      output_scalar_traverse (GFS_OUTPUT_SCALAR (output), FTT_PRE_ORDER, 
				FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				GFS_OUTPUT_SCALAR (output)->maxlevel,
				(FttCellTraverseFunc) update_histogram, output);
      if (gfs_restore_floating_point_exceptions ()) {
	gchar * s = g_strdup ("\n");
	if (h->wf)
	  s = g_strconcat (s, gfs_function_description (h->wf, FALSE), NULL);
	if (h->yf)
	  s = g_strconcat (s, "\n", gfs_function_description (h->yf, FALSE), NULL);
	/* fixme: memory leaks */
	g_message ("floating-point exception in user-defined function(s):%s", s);
	exit (1);
      }

      if (output->file && !output->dynamic)
	output->file->fp = freopen (output->format, "w", output->file->fp);
      for (i = 0; i < h->n; i++)
	if (h->w[i] > 0.) {
	  fprintf (output->file->fp, "%g %g", h->x[i]/h->w[i], h->w[i]/h->W);
	  if (h->yf)
	    fprintf (output->file->fp, " %g", h->y[i]/h->w[i]);
	  fputc ('\n', output->file->fp);
	}
    }
    h->last = sim->time.t;
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_histogram_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_histogram_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_scalar_histogram_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_scalar_histogram_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_scalar_histogram_destroy;
}

static void gfs_output_scalar_histogram_init (GfsOutputScalarHistogram * object)
{
  GFS_OUTPUT_SCALAR (object)->min = -1.;
  GFS_OUTPUT_SCALAR (object)->max =  1.;
  GFS_OUTPUT_SCALAR (object)->autoscale = FALSE;
  object->n = 100;
  object->W = 0.;
  object->last = -1.;
}

GfsOutputClass * gfs_output_scalar_histogram_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_histogram_info = {
      "GfsOutputScalarHistogram",
      sizeof (GfsOutputScalarHistogram),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_histogram_class_init,
      (GtsObjectInitFunc) gfs_output_scalar_histogram_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_histogram_info);
  }

  return klass;
}

/** \endobject{GfsOutputScalarHistogram} */

/**
 * Computing sums for each droplet.
 * \beginobject{GfsOutputDropletSums}
 */

static void gfs_output_droplet_sums_destroy (GtsObject * object)
{
  GfsOutputDropletSums * d = GFS_OUTPUT_DROPLET_SUMS (object);
  gts_object_destroy (GTS_OBJECT (d->c));
  if (d->tag)
    gts_object_destroy (GTS_OBJECT (d->tag));

  (* GTS_OBJECT_CLASS (gfs_output_droplet_sums_class ())->parent_class->destroy) (object);
}

static void gfs_output_droplet_sums_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_droplet_sums_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsOutputDropletSums * d = GFS_OUTPUT_DROPLET_SUMS (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  gfs_function_read (d->c, domain, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_STRING) {
    if (!(d->tag = gfs_domain_get_or_add_variable (domain, fp->token->str, "Droplet index"))) {
      gts_file_error (fp, "`%s' is a reserved variable name", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }
}

static void gfs_output_droplet_sums_write (GtsObject * o, FILE * fp)
{
  GfsOutputDropletSums * d = GFS_OUTPUT_DROPLET_SUMS (o);

  (* GTS_OBJECT_CLASS (gfs_output_droplet_sums_class ())->parent_class->write) (o, fp);

  gfs_function_write (d->c, fp);
  if (d->tag)
    fprintf (fp, " %s", d->tag->name);
}

typedef struct {
  double v, f;
} VolumePair;

typedef struct {
  GfsVariable * s, * c, * tag;
  VolumePair * v;
  guint n;
  GfsFunction * fc;
} DropSumsPar;

static void droplet_sums (FttCell * cell, DropSumsPar * p)
{
  guint i = GFS_VALUE (cell, p->tag);
  if (i > 0) {
    p->v[i - 1].v += GFS_VALUE (cell, p->c)*ftt_cell_volume (cell);
    p->v[i - 1].f += GFS_VALUE (cell, p->s);
  }
}

static void compute_c (FttCell * cell, DropSumsPar * p)
{
  GFS_VALUE (cell, p->c) = gfs_function_value (p->fc, cell);
}

static int volume_sort (const void * p1, const void * p2)
{
  VolumePair * a = (VolumePair *) p1;
  VolumePair * b = (VolumePair *) p2;
  return a->v < b->v ? 1 : -1;
}

static gboolean gfs_output_droplet_sums_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_droplet_sums_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputDropletSums * d = GFS_OUTPUT_DROPLET_SUMS (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    DropSumsPar p;
    p.s = GFS_OUTPUT_SCALAR (event)->v;
    p.c = gfs_function_get_variable (d->c);
    if (!p.c) {
      p.c = gfs_temporary_variable (domain);
      p.fc = d->c;
      gfs_catch_floating_point_exceptions ();
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_c, &p);
      gfs_restore_fpe_for_function (p.fc);
    }
    p.tag = d->tag ? d->tag : gfs_temporary_variable (domain);
    p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
    if (p.n > 0) {
      p.v = g_malloc0 (p.n*sizeof (VolumePair));
      output_scalar_traverse (GFS_OUTPUT_SCALAR (event), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) droplet_sums, &p);
#ifdef HAVE_MPI
      if (domain->pid >= 0) {
	VolumePair * gv = g_malloc0 (p.n*sizeof (VolumePair));
	MPI_Allreduce (p.v, gv, p.n*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	g_free (p.v);
	p.v = gv;
      }
#endif /* HAVE_MPI */
      qsort (p.v, p.n, sizeof (VolumePair), volume_sort);
      gchar * f = GFS_OUTPUT_SCALAR (event)->format;
      gchar * format;
      if (f)
	format = g_strdup_printf ("%s %%d %s\n", f, f);
      else
	format = g_strdup ("%g %d %.12g\n");
      guint i;
      for (i = 0; i < p.n; i++)
	fprintf (GFS_OUTPUT (event)->file->fp, format, sim->time.t, i + 1, p.v[i].f);
      g_free (p.v);
      g_free (format);
    }
    if (p.tag != d->tag)
      gts_object_destroy (GTS_OBJECT (p.tag));
    if (!gfs_function_get_variable (d->c))
      gts_object_destroy (GTS_OBJECT (p.c));
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_droplet_sums_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_droplet_sums_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_droplet_sums_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_droplet_sums_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_droplet_sums_destroy;
}

static void gfs_output_droplet_sums_init (GfsOutputDropletSums * d)
{
  d->c = gfs_function_new (gfs_function_class (), 0.);
}

GfsOutputClass * gfs_output_droplet_sums_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_droplet_sums_info = {
      "GfsOutputDropletSums",
      sizeof (GfsOutputDropletSums),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_droplet_sums_class_init,
      (GtsObjectInitFunc) gfs_output_droplet_sums_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_droplet_sums_info);
  }

  return klass;
}

/** \endobject{GfsOutputDropletSums} */

/**
 * Computing differences to a reference solution.
 * \beginobject{GfsOutputErrorNorm}
 */

static void output_error_norm_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_ERROR_NORM (o)->s));
  gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_ERROR_NORM (o)->w));

  (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->destroy) (o);
}

static void output_error_norm_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputErrorNorm * n;

  if (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  n = GFS_OUTPUT_ERROR_NORM (*o);
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
      gts_file_error (fp, "expecting a parameter");
      return;
    }
    else if (!strcmp (fp->token->str, "unbiased")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_INT) {
	gts_file_error (fp, "expecting an integer");
	return;
      }
      n->unbiased = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "relative")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_INT) {
	gts_file_error (fp, "expecting an integer");
	return;
      }
      n->relative = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "s")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (n->s, gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;
    }
    else if (!strcmp (fp->token->str, "w")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (n->w, gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;
    }
    else if (!strcmp (fp->token->str, "v")) {
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a variable name");
	return;
      }
      if (!(n->v = gfs_domain_get_or_add_variable (domain, fp->token->str, "Error field"))) {
	gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
	return;
      }
      gts_file_next_token (fp);
    }
    else {
      gts_file_error (fp, "unknown identifier `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void output_error_norm_write (GtsObject * o, FILE * fp)
{
  GfsOutputErrorNorm * n = GFS_OUTPUT_ERROR_NORM (o);

  if (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->write) 
      (o, fp);
  fputs (" { s = ", fp);
  gfs_function_write (n->s, fp);
  fputs (" w = ", fp);
  gfs_function_write (n->w, fp);
  fprintf (fp, " unbiased = %d relative = %d", n->unbiased, n->relative);
  if (n->v)
    fprintf (fp, " v = %s }", n->v->name);
  else
    fputs (" }", fp);
}

static void reference_solution (FttCell * cell, GfsOutputScalar * o)
{
  GFS_VALUE (cell, GFS_OUTPUT_ERROR_NORM (o)->v) = 
    gfs_function_value (GFS_OUTPUT_ERROR_NORM (o)->s, cell);
}

static void substract (FttCell * cell, GfsOutputScalar * o)
{
  GFS_VALUE (cell, GFS_OUTPUT_ERROR_NORM (o)->v) = GFS_VALUE (cell, o->v) -
    GFS_VALUE (cell, GFS_OUTPUT_ERROR_NORM (o)->v);
}

static void compute_error (FttCell * cell, GfsOutputScalar * o)
{
  GFS_VALUE (cell, GFS_OUTPUT_ERROR_NORM (o)->v) = GFS_VALUE (cell, o->v) -
    gfs_function_value (GFS_OUTPUT_ERROR_NORM (o)->s, cell);
}

static void remove_bias (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsNorm * norm = data[1];
  GFS_VALUE (cell, v) -= norm->bias;
}

static gboolean gfs_output_error_norm_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsOutputErrorNorm * enorm = GFS_OUTPUT_ERROR_NORM (event);
    GfsVariable * v = enorm->v;
    GfsNorm norm, snorm = { 0., 0., 0., 0. };

    if (v == NULL)
      enorm->v = gfs_temporary_variable (domain);
    if (enorm->relative) {
      gfs_catch_floating_point_exceptions ();
      output_scalar_traverse (output, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
			      output->maxlevel,
			      (FttCellTraverseFunc) reference_solution, output);
      gfs_restore_fpe_for_function (enorm->s);
      snorm = gfs_domain_norm_variable (domain, enorm->v, enorm->w,
					FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
					output->maxlevel,
					output->condition ? cell_condition : NULL,
					output->condition);
      output_scalar_traverse (output, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
			      output->maxlevel,
			      (FttCellTraverseFunc) substract, output);
    }
    else {
      gfs_catch_floating_point_exceptions ();
      output_scalar_traverse (output, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
			      output->maxlevel,
			      (FttCellTraverseFunc) compute_error, output);
      gfs_restore_fpe_for_function (enorm->s);
    }
    norm = gfs_domain_norm_variable (domain, enorm->v, enorm->w,
				     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				     output->maxlevel,
				     output->condition ? cell_condition : NULL,
				     output->condition);
    if (GFS_OUTPUT_ERROR_NORM (event)->unbiased) {
      gpointer data[2];

      data[0] = enorm->v;
      data[1] = &norm;
      output_scalar_traverse (output, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
			      output->maxlevel,
			      (FttCellTraverseFunc) remove_bias, data);
      norm = gfs_domain_norm_variable (domain, enorm->v, enorm->w,
				       FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				       output->maxlevel,
				       output->condition ? cell_condition : NULL,
				       output->condition);
    }
    if (v == NULL) {
      gts_object_destroy (GTS_OBJECT (enorm->v));
      enorm->v = NULL;
    }
    if (enorm->relative) {
      if (snorm.first > 0.)  norm.first  /= snorm.first;
      if (snorm.second > 0.) norm.second /= snorm.second;
      if (snorm.infty > 0.)  norm.infty  /= snorm.infty;
    }
    gchar * format;
    if (output->format) {
      gchar * f = output->format;
      format = g_strdup_printf ("%%s time: %s first: %s second: %s infty: %s bias: %s\n",
				f, f, f, f, f);
    }
    else
      format = g_strdup ("%s time: %g first: %10.3e second: %10.3e infty: %10.3e bias: %10.3e\n");
    fprintf (GFS_OUTPUT (event)->file->fp, format,
	     output->name, sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias);
    g_free (format);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_error_norm_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = output_error_norm_destroy;
  GTS_OBJECT_CLASS (klass)->read = output_error_norm_read;
  GTS_OBJECT_CLASS (klass)->write = output_error_norm_write;
  GFS_EVENT_CLASS (klass)->event = gfs_output_error_norm_event;
}

static void output_error_norm_init (GfsOutputErrorNorm * e)
{
  e->s = gfs_function_new (gfs_function_class (), 0.);
  e->w = gfs_function_new (gfs_function_class (), 1.);
}

GfsOutputClass * gfs_output_error_norm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_error_norm_info = {
      "GfsOutputErrorNorm",
      sizeof (GfsOutputErrorNorm),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_error_norm_class_init,
      (GtsObjectInitFunc) output_error_norm_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_error_norm_info);
  }

  return klass;
}

/** \endobject{GfsOutputErrorNorm} */

/**
 * 
 * \beginobject{GfsOutputCorrelation}
 */

static void compute_correlation (FttCell * cell, gpointer * data)
{
  GfsOutputScalar * o = data[0];
  gdouble * bias = data[1];
  gdouble * sum = data[2];
  gdouble * sumref = data[3];
  gdouble v, ref, w;

  ref = gfs_function_value (GFS_OUTPUT_ERROR_NORM (o)->s, cell);
  v = GFS_VALUE (cell, o->v) - *bias;
  w = gfs_cell_volume (cell, o->v->domain);
  *sumref += ref*ref*w;
  *sum += v*ref*w;
}

static gboolean gfs_output_correlation_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsOutputErrorNorm * enorm = GFS_OUTPUT_ERROR_NORM (event);
    GfsVariable * v = enorm->v;
    gdouble bias = 0., sum = 0., sumref = 0.;
    gpointer data[4];

    if (GFS_DOMAIN (sim)->pid != -1)
      g_assert_not_implemented ();

    if (v == NULL)
      enorm->v = gfs_temporary_variable (GFS_DOMAIN (sim));
    if (enorm->unbiased) {
      output_scalar_traverse (output, FTT_PRE_ORDER,
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			      output->maxlevel,
			      (FttCellTraverseFunc) compute_error, output);
      bias = gfs_domain_norm_variable (GFS_DOMAIN (sim), enorm->v, NULL,
				       FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				       output->maxlevel,
				       output->condition ? cell_condition : NULL,
				       output->condition).bias;
    }
    data[0] = output;
    data[1] = &bias;
    data[2] = &sum;
    data[3] = &sumref;
    gfs_catch_floating_point_exceptions ();
    output_scalar_traverse (output, FTT_PRE_ORDER,
			    FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			    output->maxlevel,
			    (FttCellTraverseFunc) compute_correlation, data);
    gfs_restore_fpe_for_function (enorm->s);
    if (v == NULL) {
      gts_object_destroy (GTS_OBJECT (enorm->v));
      enorm->v = NULL;
    }
    gchar * format;
    if (output->format) {
      gchar * f = output->format;
      format = g_strdup_printf ("%%s time: %s %s\n", f, f);
    }
    else
      format = g_strdup ("%s time: %g %10.3e\n");
    fprintf (GFS_OUTPUT (event)->file->fp, format,
	     output->name, sim->time.t, sumref > 0. ? sum/sumref : 0.);
    g_free (format);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_correlation_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_correlation_event;
}

GfsOutputClass * gfs_output_correlation_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_correlation_info = {
      "GfsOutputCorrelation",
      sizeof (GfsOutputErrorNorm),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_correlation_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_error_norm_class ()),
				  &gfs_output_correlation_info);
  }

  return klass;
}

/** \endobject{GfsOutputCorrelation} */

/**
 *
 * \beginobject{GfsOutputSquares}
 */

static gboolean gfs_output_squares_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_squares_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    
    gfs_write_squares (GFS_DOMAIN (sim), 
		       output->v, output->min, output->max,
		       FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
		       output->maxlevel, NULL, 
		       GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_squares_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_squares_event;
}

GfsOutputClass * gfs_output_squares_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_squares_info = {
      "GfsOutputSquares",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_squares_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_squares_info);
  }

  return klass;
}

/** \endobject{GfsOutputSquares} */

/**
 *
 * \beginobject{GfsOutputStreamline}
 */

static void gfs_output_streamline_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputStreamline * l = GFS_OUTPUT_STREAMLINE (*o);

  if (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return;
  }
  l->p.x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return;
  }
  l->p.y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return;
  }
  l->p.z = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_output_streamline_write (GtsObject * o, FILE * fp)
{
  GfsOutputStreamline * l = GFS_OUTPUT_STREAMLINE (o);

  if (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g %g %g", l->p.x, l->p.y, l->p.z);
}

static gboolean gfs_output_streamline_event (GfsEvent * event, 
					    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class)->event)
      (event,sim)) {
    FttVector p = GFS_OUTPUT_STREAMLINE (event)->p;
    gfs_simulation_map (sim, &p);
    GList * stream = gfs_streamline_new (GFS_DOMAIN (sim),
					 gfs_domain_velocity (GFS_DOMAIN (sim)),
					 p,
					 GFS_OUTPUT_SCALAR (event)->v,
					 0., 0.,
					 TRUE,
					 NULL, NULL);
    /* fixme: mapping is not taken into account */
    gfs_streamline_write (stream, GFS_OUTPUT (event)->file->fp);
    gfs_streamline_destroy (stream);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_streamline_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_streamline_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_streamline_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_streamline_write;
}

GfsOutputClass * gfs_output_streamline_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_streamline_info = {
      "GfsOutputStreamline",
      sizeof (GfsOutputStreamline),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_streamline_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_streamline_info);
  }

  return klass;
}

/** \endobject{GfsOutputStreamline} */

/**
 * Writing 2D images.
 * \beginobject{GfsOutputPPM}
 */

static void gfs_output_ppm_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
#if (!FTT_2D)
  if (!GFS_IS_OCEAN (gfs_object_simulation (*o))) {
    gts_file_error (fp, 
		    "In more than two dimensions PPM output is possible\n"
		    "only for GfsOcean simulations");
    return;
  }
#endif /* 3D */
}

static gboolean gfs_output_ppm_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class)->event) 
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
#if FTT_2D
    GfsDomain * domain = GFS_DOMAIN (sim);
#else /* 3D */
    GfsDomain * domain = GFS_IS_OCEAN (sim) ? GFS_OCEAN (sim)->toplayer : GFS_DOMAIN (sim);
#endif /* 3D */

    gfs_write_ppm (domain,
		   output->condition,
		   output->v, output->min, output->max,
		   FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, output->maxlevel,
		   GFS_OUTPUT (event)->file->fp,
		   GFS_OUTPUT (event)->parallel);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_ppm_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_output_ppm_read;
  GFS_EVENT_CLASS (klass)->event = gfs_output_ppm_event;
}

GfsOutputClass * gfs_output_ppm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_ppm_info = {
      "GfsOutputPPM",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_ppm_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_ppm_info);
  }

  return klass;
}

/** \endobject{GfsOutputPPM} */

/**
 * Writing 2D ESRI raster files.
 * \beginobject{GfsOutputGRD}
 */

static gboolean gfs_output_grd_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_grd_class ())->parent_class)->event) 
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    gfs_write_grd (sim,
		   output->condition,
		   output->v,
		   FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, output->maxlevel,
		   GFS_OUTPUT (event)->file->fp,
		   GFS_OUTPUT (event)->parallel,
		   TRUE);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_grd_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_grd_event;
}

GfsOutputClass * gfs_output_grd_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_grd_info = {
      "GfsOutputGRD",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_grd_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_grd_info);
  }

  return klass;
}

/** \endobject{GfsOutputGRD} */

/**
 * Calling the write method of a given object.
 * \beginobject{GfsOutputObject}
 */

static void output_object_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_class())->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }
  GfsOutputObject * output = GFS_OUTPUT_OBJECT (*o);
  output->object = gfs_object_from_name (GFS_DOMAIN (gfs_object_simulation (*o)), fp->token->str);
  if (output->object == NULL)
    gts_file_error (fp, "unknown object '%s'", fp->token->str);
  else
    gts_file_next_token (fp);
}

static void output_object_write (GtsObject * o, FILE * fp)
{
  GfsOutputObject * output = GFS_OUTPUT_OBJECT (o);

  (* GTS_OBJECT_CLASS (gfs_output_object_class ())->parent_class->write) (o, fp);
  
  fprintf (fp, " %s\n", GFS_EVENT (output->object)->name );
}

static gboolean output_object_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GtsObject * object = GFS_OUTPUT_OBJECT (event)->object;
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    object->klass->write (object, fp);
    fprintf (fp, "\n");

    return TRUE;
  }
  return FALSE;
}

static void gfs_output_object_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = output_object_read;
  GTS_OBJECT_CLASS (klass)->write = output_object_write;
  GFS_EVENT_CLASS(klass)->event = output_object_event;
}

GfsOutputClass * gfs_output_object_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_object_info = {
      "GfsOutputObject",
      sizeof (GfsOutputObject),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_object_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_object_info);
  }

  return klass;
}

/** \endobject{GfsOutputObject} */
