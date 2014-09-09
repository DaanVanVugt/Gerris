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

#include <ode/ode.h>
#include "moving.h"

static dWorldID world = 0;

/* GfsSurfaceBcODE: Header */

typedef struct _GfsSurfaceBcODE         GfsSurfaceBcODE;

struct _GfsSurfaceBcODE {
  /*< private >*/
  GfsSurfaceGenericBc parent;

  /*< public >*/
  dBodyID body;
  FttComponent c;
};

#define GFS_SURFACE_BC_ODE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceBcODE,\
					         gfs_surface_bc_ode_class ())
#define GFS_IS_SURFACE_BC_ODE(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_bc_ode_class ()))

static GfsSurfaceGenericBcClass * gfs_surface_bc_ode_class (void);

/* GfsSurfaceBcODE: Object */

static void surface_bc_ode_read (GtsObject ** o, GtsFile * fp)
{
  /* defined only through SolidMovingODE */
}

static void surface_bc_ode_write (GtsObject * o, FILE * fp)
{
  /* defined only through SolidMovingODE */
}

static void surface_bc_ode_bc (FttCell * cell, GfsSurfaceGenericBc * b)
{
  GfsSurfaceBcODE * bc = GFS_SURFACE_BC_ODE (b);
  FttVector * p = &GFS_STATE (cell)->solid->ca;
  dVector3 v;
  cell->flags |= GFS_FLAG_DIRICHLET;
  dBodyGetPointVel (bc->body, p->x, p->y, p->z, v);
  GFS_STATE (cell)->solid->fv = v[bc->c];
}

static void surface_bc_ode_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = surface_bc_ode_read;
  GTS_OBJECT_CLASS (klass)->write = surface_bc_ode_write;
  klass->bc = surface_bc_ode_bc;
}

static GfsSurfaceGenericBcClass * gfs_surface_bc_ode_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_bc_info = {
      "GfsSurfaceBcODE",
      sizeof (GfsSurfaceBcODE),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) surface_bc_ode_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ()),
				  &gfs_surface_bc_info);
  }

  return klass;
}

/* GfsSolidMovingODE: Header */

typedef struct _GfsSolidMovingODE         GfsSolidMovingODE;

struct _GfsSolidMovingODE {
  /*< private >*/
  GfsSolidMoving parent;
  dBodyID body;

  /*< public >*/  
};

static GfsEventClass * gfs_solid_moving_ode_class (void);

#define GFS_SOLID_MOVING_ODE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSolidMovingODE,\
					         gfs_solid_moving_ode_class ())
#define GFS_IS_SOLID_MOVING_ODE(obj)         (gts_object_is_from_class (obj,\
						 gfs_solid_moving_ode_class ()))

/* GfsSolidMovingODE: Object */

static void solid_moving_ode_destroy (GtsObject * object)
{
  dBodyDestroy (GFS_SOLID_MOVING_ODE (object)->body);

  (* GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class->destroy) (object);
}

static void solid_moving_ode_read (GtsObject ** o, GtsFile * fp)
{
  GfsSolidMovingODE * solid = GFS_SOLID_MOVING_ODE (*o);
  
  (* GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariable ** v = gfs_domain_velocity (GFS_DOMAIN (gfs_object_simulation (solid)));
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (v[c]->surface_bc) {
      gts_file_error (fp, "variable `%s' already has a surface boundary condition", v[c]->name);
      return;
    }
    else {
      GfsSurfaceGenericBc * bc = 
	GFS_SURFACE_GENERIC_BC (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_bc_ode_class ())));
      bc->v = v[c];
      bc->v->surface_bc = bc;
      GFS_SURFACE_BC_ODE (bc)->body = solid->body;
      GFS_SURFACE_BC_ODE (bc)->c = c;
    }

  if (fp->type == '{') {
    gdouble vx = 0., vy = 0., vz = 0.;
    gdouble gx = 0., gy = 0., gz = 0.;
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "vx", TRUE, &vx},
      {GTS_DOUBLE, "vy", TRUE, &vy},
      {GTS_DOUBLE, "vz", TRUE, &vz},
      {GTS_DOUBLE, "gx", TRUE, &gx},
      {GTS_DOUBLE, "gy", TRUE, &gy},
      {GTS_DOUBLE, "gz", TRUE, &gz},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    dBodySetLinearVel (solid->body, vx, vy, vz);
    dWorldSetGravity (world, gx, gy, gz);
  }
}

static gboolean solid_moving_ode_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class)->event) 
      (event, sim)) {
    FttVector pf, vf, pm, vm;
    dBodyID body = GFS_SOLID_MOVING_ODE (event)->body;
    gfs_domain_solid_force (GFS_DOMAIN (sim), &pf, &vf, &pm, &vm, NULL);
    dBodyAddForce (body, pf.x + vf.x, pf.y + vf.y, pf.z + vf.z);
    dBodyAddTorque (body, pm.x + vm.x, pm.y + vm.y, pm.z + vm.z);
    dWorldStep (world, sim->advection_params.dt);
    return TRUE;
  }
  return FALSE;
}

static void solid_moving_ode_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = solid_moving_ode_destroy;
  GTS_OBJECT_CLASS (klass)->read = solid_moving_ode_read;
  klass->event = solid_moving_ode_event;
}

static void solid_moving_ode_init (GfsSolidMovingODE * solid)
{
  solid->body = dBodyCreate (world);
  dMass mass;
  dMassSetSphereTotal (&mass, 1., 0.25*0.3);
  dBodySetMass (solid->body, &mass);
}

static GfsEventClass * gfs_solid_moving_ode_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo solid_moving_ode_info = {
      "GfsSolidMovingODE",
      sizeof (GfsSolidMovingODE),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) solid_moving_ode_class_init,
      (GtsObjectInitFunc) solid_moving_ode_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_moving_class ()), 
				  &solid_moving_ode_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "ode";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  dInitODE ();
  world = dWorldCreate ();
  gfs_solid_moving_ode_class ();
  return NULL;
}
