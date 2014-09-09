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
#include <stdlib.h>
#include <math.h>
#include "event.h"
#include "solid.h"
#include "config.h"

#define CW260_F77 F77_FUNC (cw260, CW260)
#define KMTS_F77 F77_FUNC (kmts, KMTS)

void CW260_F77 (const float * zd,
		const float * zt, 
		const float * zh, 
		const float * zu, 
		const int * nverb, 
		int * nfun, 
		float * zel);

void KMTS_F77  (const int * nfun,
		const float * xx,
		const float * yy,
		const float * tt,
		float * uu, 
		float * vv,
		float * ut,
		float * vt,
		float * du,
		float * dv,
		float * etah);

/* GfsInitStokesWave: Header */

typedef struct _GfsInitStokesWave         GfsInitStokesWave;

struct _GfsInitStokesWave {
  /*< private >*/
  GfsGenericInit parent;
  
  /*< public >*/
  gdouble steepness, depth;
};

#define GFS_INIT_STOKES_WAVE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitStokesWave,\
					         gfs_init_stokes_wave_class ())

GfsGenericInitClass * gfs_init_stokes_wave_class  (void);

/* GfsInitStokesWave: Object */

#define WAVELENGTH 100. /* metres */

static int order = 0;

static void gfs_init_stokes_wave_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_stokes_wave_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_stokes_wave_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsInitStokesWave * w = GFS_INIT_STOKES_WAVE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "depth",     TRUE},
    {GTS_DOUBLE, "steepness", TRUE},
    {GTS_NONE}
  };
  var[0].data = &w->depth;
  var[1].data = &w->steepness;
  gts_file_assign_variables (fp, var);  
  if (fp->type == GTS_ERROR)
    return;

  float zd = WAVELENGTH*w->depth, zh = w->steepness*WAVELENGTH/M_PI, zu = 0., l;
  int nverb = 0;
  
  float min = 1., max = 100., zt = (min + max)/2.;
  do {
    CW260_F77 (&zd, &zt, &zh, &zu, &nverb, &order, &l);
    fprintf (stderr, "# order: %d wavelength: %g period: %g\n", 
	     order, l, zt);
    if (l > WAVELENGTH)
      max = zt;
    if (l < WAVELENGTH)
      min = zt;
    zt = (min + max)/2.;
  } while (fabs (l - WAVELENGTH) > 1e-4);
}

static void gfs_init_stokes_wave_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_stokes_wave_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_stokes_wave_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " { depth = %g steepness = %g }\n",
	   GFS_INIT_STOKES_WAVE (o)->depth, 
	   GFS_INIT_STOKES_WAVE (o)->steepness);
}

static void init_velocity (FttCell * cell, GfsVariable ** velocity)
{
  FttVector p;
  gfs_cell_cm (cell, &p);
  float x = (p.x + 0.5)*WAVELENGTH, y = p.y*WAVELENGTH, t = 0., u, v, ut, vt, du, dv, etah;
  KMTS_F77 (&order, &x, &y, &t, &u, &v, &ut, &vt, &du, &dv, &etah);
  GFS_VALUE (cell, velocity[0]) = u/sqrt(WAVELENGTH*9.81);
  GFS_VALUE (cell, velocity[1]) = v/sqrt(WAVELENGTH*9.81);
}

static gdouble stokes_height (gdouble x, gdouble y, gdouble z, gdouble t)
{
  float xx = (x + 0.5)*WAVELENGTH, yy = WAVELENGTH, u, v, ut, vt, du, dv, etah, t1 = 0.;
  KMTS_F77 (&order, &xx, &yy, &t1, &u, &v, &ut, &vt, &du, &dv, &etah);
  return etah/WAVELENGTH - y;
}

static gboolean gfs_init_stokes_wave_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_stokes_wave_class ())->parent_class)->event) 
      (event, sim)) {
    GfsVariable ** velocity = gfs_domain_velocity (GFS_DOMAIN (sim));
    GfsVariable * t = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "T");
    g_assert (velocity);
    g_assert (t);
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_velocity, velocity);
    GfsSurface * surface = GFS_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
    surface->f = gfs_function_spatial_new (gfs_function_spatial_class (), stokes_height);
    gfs_object_simulation_set (surface->f, sim);
    gfs_domain_init_fraction (GFS_DOMAIN (sim), GFS_GENERIC_SURFACE (surface), t);
    gts_object_destroy (GTS_OBJECT (surface));
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_stokes_wave_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_stokes_wave_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_stokes_wave_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_stokes_wave_write;
}

static void gfs_init_stokes_wave (GfsInitStokesWave * w)
{
  w->depth = 0.5;
  w->steepness = 0.3;
}

GfsGenericInitClass * gfs_init_stokes_wave_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_stokes_wave_info = {
      "GfsInitStokesWave",
      sizeof (GfsInitStokesWave),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_stokes_wave_class_init,
      (GtsObjectInitFunc) gfs_init_stokes_wave,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_stokes_wave_info);
  }

  return klass;
}

const gchar gfs_module_name[] = "stokes";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_init_stokes_wave_class ();
  return NULL;
}
