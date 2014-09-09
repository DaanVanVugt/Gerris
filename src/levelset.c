/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2006 National Institute of Water and Atmospheric Research
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
 * \brief GfsVariableDistance.
 */

#include <stdlib.h>
#include "levelset.h"
#include "vof.h"

/**
 * Signed distance to a VOF interface.
 * \beginobject{GfsVariableDistance}
 */

static void variable_distance_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (c)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(GFS_VARIABLE_DISTANCE (*o)->v = 
	gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (GFS_VARIABLE (*o)->description)
    g_free (GFS_VARIABLE (*o)->description);
  GFS_VARIABLE (*o)->description = g_strjoin (" ", 
					       "Distance to the interface defined by tracer",
					       fp->token->str, NULL);
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "stencil", TRUE},
      {GTS_NONE}
    };
    var[0].data = &GFS_VARIABLE_DISTANCE (*o)->stencil; 
    gts_file_assign_variables (fp, var);
  }
}

static void variable_distance_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", GFS_VARIABLE_DISTANCE (o)->v->name);
  if (GFS_VARIABLE_DISTANCE (o)->stencil)
    fputs (" { stencil = 1 }", fp);
}

static gdouble vof_distance2 (FttCell * cell, GtsPoint * t, gpointer v)
{
  gdouble f = GFS_VALUE (cell, GFS_VARIABLE (v));
  
  if (GFS_IS_FULL (f))
    return GFS_NODATA;
  if (!FTT_CELL_IS_LEAF (cell))
    return ftt_cell_point_distance2_min (cell, t);
  else
    return gfs_vof_facet_distance2 (cell, v, t);
}

static void distance_for_stencil (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * s2 = data[2];

  if (GFS_VALUE (cell, s2)) {
    GfsVariableDistance * l = GFS_VARIABLE_DISTANCE (v);
    GtsPoint p;
    gdouble d2;
    
    ftt_cell_pos (cell, (FttVector *) &p.x);
    d2 = gfs_domain_cell_point_distance2 (v->domain, &p, vof_distance2, l->v, NULL);
    GFS_VALUE (cell, v) = GFS_VALUE (cell, l->v) > 0.5 ? sqrt (d2) : -sqrt (d2);
  }
  else
    GFS_VALUE (cell, v) = GFS_NODATA;
}

static void distance (FttCell * cell, GfsVariable * v)
{
  GfsVariableDistance * l = GFS_VARIABLE_DISTANCE (v);
  GtsPoint p;
  gdouble d2;
    
  ftt_cell_pos (cell, (FttVector *) &p.x);
  d2 = gfs_domain_cell_point_distance2 (v->domain, &p, vof_distance2, l->v, NULL);
  GFS_VALUE (cell, v) = GFS_VALUE (cell, l->v) > 0.5 ? sqrt (d2) : -sqrt (d2);
}

static void stencil_interpolate (FttCell * cell, gpointer * data)
{
  GfsVariableDistance * v = data[0];
  gdouble f = GFS_VALUE (cell, v->v);
  
  if (!GFS_IS_FULL (f))
    gfs_interpolate_stencil (cell, data[1]);
}

static void stencil_gradient (FttCell * cell, gpointer * data)
{
  GfsVariable * s1 = data[1];
  GfsVariable * s2 = data[2];
  FttComponent c;

  if (GFS_VALUE (cell, s1))
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_center_gradient_stencil (cell, c, s2->i);
}

static void variable_distance_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariableDistance * v = GFS_VARIABLE_DISTANCE (event);

  gfs_domain_timer_start (domain, "distance");

  if (v->stencil) { /* fixme: this "acceleration technique"
		       i.e. computing distance only in a band around
		       the interface seems to be slower than computing
		       the distance function everywhere! */
    gpointer data[3], tmp;
    data[0] = v;
    data[1] = gfs_temporary_variable (domain);
    data[2] = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, data[1]);
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) stencil_interpolate, data);
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, data[2]);
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) stencil_gradient, data);
    tmp = data[1]; data[1] = data[2]; data[2] = tmp;
    gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) stencil_gradient, data);

    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) v->v->fine_coarse, v->v);
    gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) distance_for_stencil, data,
			 GFS_VARIABLE (event),GFS_VARIABLE (event));
    gts_object_destroy (data[1]);
    gts_object_destroy (data[2]);
  }
  else
    gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) distance, v,
			 GFS_VARIABLE (event),GFS_VARIABLE (event));

  gfs_domain_timer_stop (domain, "distance");
}

static gboolean variable_distance_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class)->event)
      (event, sim)) {
    if (!GFS_VARIABLE_DISTANCE (event)->first_done) {
      variable_distance_event_half (event, sim);
      GFS_VARIABLE_DISTANCE (event)->first_done = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void variable_distance_class_init (GtsObjectClass * klass)
{
  klass->read = variable_distance_read;
  klass->write = variable_distance_write;
  GFS_EVENT_CLASS (klass)->event = variable_distance_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_distance_event_half;
}

static void variable_distance_init (GfsVariable * v)
{
  v->units = 1.;
}

GfsVariableClass * gfs_variable_distance_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_distance_info = {
      "GfsVariableDistance",
      sizeof (GfsVariableDistance),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_distance_class_init,
      (GtsObjectInitFunc) variable_distance_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_distance_info);
  }

  return klass;
}

/** \endobject{GfsVariableDistance} */
