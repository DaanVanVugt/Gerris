/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009 National Institute of Water and Atmospheric Research
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
 * \brief Lagrangian particles.
 */

#include <stdlib.h>
#include "particle.h"

/**
 * Lagrangian particules.
 * \beginobject{GfsParticle}
 */

static gboolean gfs_particle_event (GfsEvent * event, 
				    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particle_class ())->parent_class)->event)
      (event, sim)) {
    GfsParticle * p = GFS_PARTICLE (event);
    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);
    gfs_domain_advect_point (GFS_DOMAIN (sim), &pos, sim->advection_params.dt);
    gfs_simulation_map_inverse (sim, &pos);
    p->pos = pos;
    return TRUE;
  }
  return FALSE;
}

static void gfs_particle_read (GtsObject ** o, GtsFile * fp)
{
  GfsParticle * p = GFS_PARTICLE(*o);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (Id)");
    return;
  }
  p->id = atoi (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return;
  }
  p->pos.x = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return;
  }
  p->pos.y = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return;
  }
  p->pos.z = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_particle_write (GtsObject * o, FILE * fp)
{
  GfsParticle * p = GFS_PARTICLE(o);
  fprintf (fp, " %d %g %g %g", p->id, p->pos.x, p->pos.y, p->pos.z);
}

static void gfs_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event =  gfs_particle_event;
  GTS_OBJECT_CLASS (klass)->read =  gfs_particle_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_particle_write;
}

GfsEventClass * gfs_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_info = {
      "GfsParticle",
      sizeof (GfsParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particle_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_particle_info);
  }
  return klass;
}

/** \endobject{GfsParticle} */
