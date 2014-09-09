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

#include <proj_api.h>
#include "simulation.h"
#include "map.h"

/* GfsMapProjection: Header */

typedef struct _GfsMapProjection         GfsMapProjection;

struct _GfsMapProjection {
  /*< private >*/
  GfsMap parent;
  projPJ pj;
  gdouble cosa, sina;

  /*< public >*/
  gdouble lon, lat, angle;
};

#define GFS_MAP_PROJECTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsMapProjection,\
					         gfs_map_projection_class ())
#define GFS_IS_MAP_PROJECTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_projection_class ()))

GfsMapClass * gfs_map_projection_class  (void);

/* GfsMapProjection: Object */

static void gfs_map_projection_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_projection_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "lon",    TRUE},
    {GTS_DOUBLE, "lat",    TRUE},
    {GTS_DOUBLE, "angle",  TRUE},
    {GTS_NONE}
  };
  GfsMapProjection * map = GFS_MAP_PROJECTION (*o);
  var[0].data = &map->lon;
  var[1].data = &map->lat;
  var[2].data = &map->angle;

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  map->cosa = cos (map->angle*DEG_TO_RAD);
  map->sina = sin (map->angle*DEG_TO_RAD);

  char * parms[] = {
    "proj=lcc", /* Lambert Conformal Conic */
    NULL, NULL, NULL, NULL
  };
  parms[1] = g_strdup_printf ("lon_0=%lf", map->lon);
  parms[2] = g_strdup_printf ("lat_0=%lf", map->lat);
  parms[3] = g_strdup_printf ("lat_1=%lf", map->lat);
  parms[4] = g_strdup_printf ("lat_2=%lf", map->lat);
  map->pj = pj_init (sizeof(parms)/sizeof(char *), parms);
  if (!map->pj)
    gts_file_error (fp, "cannot initialise projection");
  g_free (parms[1]);
  g_free (parms[2]);
  g_free (parms[3]);
  g_free (parms[4]);
}

static void gfs_map_projection_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_projection_class ())->parent_class->write) (o, fp);
  GfsMapProjection * map = GFS_MAP_PROJECTION (o);
  fprintf (fp, " { lon = %.8g lat = %.8g angle = %g }",
	   map->lon, map->lat, map->angle);
}

static void gfs_map_projection_destroy (GtsObject * object)
{
  if (GFS_MAP_PROJECTION (object)->pj)
    pj_free (GFS_MAP_PROJECTION (object)->pj);
  (* GTS_OBJECT_CLASS (gfs_map_projection_class ())->parent_class->destroy) (object);
}

static void projection_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  projLP idata;
  projXY odata;
  GfsMapProjection * m = GFS_MAP_PROJECTION (map);
  gdouble L = gfs_object_simulation (map)->physical_params.L;
  idata.u = src->x*L*DEG_TO_RAD;
  idata.v = src->y*L*DEG_TO_RAD;
  odata = pj_fwd (idata, m->pj);
  dest->x = (odata.u*m->cosa - odata.v*m->sina)/L;
  dest->y = (odata.v*m->cosa + odata.u*m->sina)/L;
  dest->z = src->z;
}

static void projection_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  projLP odata;
  projXY idata;
  GfsMapProjection * m = GFS_MAP_PROJECTION (map);
  gdouble L = gfs_object_simulation (map)->physical_params.L;
  idata.u = (src->x*m->cosa + src->y*m->sina)*L;
  idata.v = (src->y*m->cosa - src->x*m->sina)*L;
  odata = pj_inv (idata, GFS_MAP_PROJECTION (map)->pj);
  dest->x = odata.u*RAD_TO_DEG/L;
  dest->y = odata.v*RAD_TO_DEG/L;
  dest->z = src->z;
}

static void projection_inverse_cell (GfsMap * map, const FttVector * src, FttVector * dest)
{
  gint i;
  FttVector o = { 0., 0., 0. };
  for (i = 0; i < 4; i++) {
    o.x += src[i].x;
    o.y += src[i].y;
    o.z += src[i].z;
    projection_inverse (map, &(src[i]), &(dest[i]));
  }
  o.x /= 4.; o.y /= 4.; o.z /= 4.;
  projection_inverse (map, &o, &o);
  /* make sure we do not cross periodic longitude boundary */
  gdouble L = gfs_object_simulation (map)->physical_params.L;
  for (i = 0; i < 4; i++)
    if (dest[i].x > o.x + 180./L)
      dest[i].x -= 360./L;
    else if (dest[i].x < o.x - 180./L)
      dest[i].x += 360./L;
}

static void gfs_map_projection_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_projection_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_projection_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_map_projection_destroy;
}

static void gfs_map_projection_init (GfsMapProjection * object)
{
  /* Wellington */
  object->lon = 174.777222;
  object->lat = -41.288889;
  object->angle = 0.; object->cosa = 1.; object->sina = 0.;
  object->pj = NULL;
  GFS_MAP (object)->transform = projection_transform;
  GFS_MAP (object)->inverse = projection_inverse;
  GFS_MAP (object)->inverse_cell = projection_inverse_cell;
}

GfsMapClass * gfs_map_projection_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_projection_info = {
      "GfsMapProjection",
      sizeof (GfsMapProjection),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_projection_class_init,
      (GtsObjectInitFunc) gfs_map_projection_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()),
				  &gfs_map_projection_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "map";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_map_projection_class ();
  return NULL;
}
