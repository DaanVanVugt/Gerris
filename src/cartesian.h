/* Gerris - The GNU Flow Solver
 * Copyright (C) 2007 National Institute of Water and Atmospheric Research
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

#ifndef __CARTESIAN_H__
#define __CARTESIAN_H__

#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* GfsCartesianGrid: Header */

typedef struct _GfsCartesianGrid      GfsCartesianGrid;

struct _GfsCartesianGrid {
  /*< private >*/
  GtsObject parent;
  guint N;       /* Number of dimensions */
  guint * n;     /* Size of each dimension */
  gdouble ** x;  /* Position of each point in the grid */
  gdouble * v;   /* Data */
  gchar ** name; /* Name of each dimension */


  /*< public >*/
  /* add extra data here (if public) */
};

#define GFS_CARTESIAN_GRID(obj)            GTS_OBJECT_CAST (obj,\
							    GfsCartesianGrid, \
							    gfs_cartesian_grid_class ())
#define GFS_IS_CARTESIAN_GRID(obj)         (gts_object_is_from_class (obj,\
								      gfs_cartesian_grid_class ()))

GtsObjectClass *    gfs_cartesian_grid_class         (void);
GfsCartesianGrid *  gfs_cartesian_grid_new           (GtsObjectClass * klass);
GfsCartesianGrid *  gfs_cartesian_grid_read          (const gchar * name, 
						      GtsFile * fp);
gboolean            gfs_cartesian_grid_interpolate   (GfsCartesianGrid * g, 
						      gdouble * p, 
						      gdouble * val);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CARTESIAN_H__ */
