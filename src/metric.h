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

#ifndef __METRIC_H__
#define __METRIC_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "map.h"

/* GfsGenericMetric: Header */

#define GFS_IS_GENERIC_METRIC(obj)         (gts_object_is_from_class (obj,\
						   gfs_generic_metric_class ()))

GfsEventClass * gfs_generic_metric_class  (void);

/* GfsVariableMetric: Header */

#define GFS_IS_VARIABLE_METRIC(obj)         (gts_object_is_from_class (obj,\
						   gfs_variable_metric_class ()))

GfsVariableClass * gfs_variable_metric_class  (void);

/* GfsStoredMetric: Header */

typedef struct _GfsStoredMetric GfsStoredMetric;

struct _GfsStoredMetric {
  /*< private >*/
  GfsVariable parent;
  GfsMap * map;

  /*< public >*/
  GfsVariable * h[4], * e;
  GfsMapClass * map_class;
};

#define GFS_STORED_METRIC(obj)            GTS_OBJECT_CAST (obj,\
					           GfsStoredMetric,\
					           gfs_stored_metric_class ())
#define GFS_IS_STORED_METRIC(obj)         (gts_object_is_from_class (obj,\
						   gfs_stored_metric_class ()))

GfsVariableClass * gfs_stored_metric_class  (void);

/* GfsMetric: Header */

typedef struct _GfsMetric GfsMetric;

struct _GfsMetric {
  /*< private >*/
  GfsStoredMetric parent;

  /*< public >*/
  GfsFunction * x, * y, * z;
};

#define GFS_METRIC(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetric,\
					           gfs_metric_class ())
#define GFS_IS_METRIC(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_class ()))

GfsVariableClass * gfs_metric_class  (void);

/* GfsMetricCubed: Header */

typedef struct _GfsMetricCubed GfsMetricCubed;

struct _GfsMetricCubed {
  /*< private >*/
  GfsStoredMetric parent;

  /*< public >*/
  gint level;
};

#define GFS_METRIC_CUBED(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetricCubed,\
					           gfs_metric_cubed_class ())
#define GFS_IS_METRIC_CUBED(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_cubed_class ()))

GfsVariableClass * gfs_metric_cubed_class  (void);

/* GfsMetricLonLat: Header */

typedef struct _GfsMetricLonLat GfsMetricLonLat;

struct _GfsMetricLonLat {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * h2, * h3;
  gdouble r;
};

#define GFS_METRIC_LON_LAT(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetricLonLat,\
					           gfs_metric_lon_lat_class ())
#define GFS_IS_METRIC_LON_LAT(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_lon_lat_class ()))

GfsVariableClass * gfs_metric_lon_lat_class  (void);

/* GfsMetricStretch: Header */

typedef struct _GfsMetricStretch GfsMetricStretch;

struct _GfsMetricStretch {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  gdouble sx, sy, sz;
};

#define GFS_METRIC_STRETCH(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetricStretch,\
					           gfs_metric_stretch_class ())
#define GFS_IS_METRIC_STRETCH(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_stretch_class ()))

GfsEventClass * gfs_metric_stretch_class  (void);

/* GfsMetricCubed1: Header */

GfsVariableClass * gfs_metric_cubed1_class  (void);

/* GfsMetricVariable: Header */

typedef struct _GfsMetricVariable GfsMetricVariable;

struct _GfsMetricVariable {
  /*< private >*/
  GfsStoredMetric parent;

  /*< public >*/
  GfsVariable * x[3];
};

#define GFS_METRIC_VARIABLE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetricVariable,\
					           gfs_metric_variable_class ())
#define GFS_IS_METRIC_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_variable_class ()))

GfsVariableClass * gfs_metric_variable_class  (void);

/* GfsMetricLaplace: Header */

typedef struct _GfsMetricLaplace GfsMetricLaplace;

struct _GfsMetricLaplace {
  /*< private >*/
  GfsMetricVariable parent;

  /*< public >*/
  gboolean conformal, spherical;
};

#define GFS_METRIC_LAPLACE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsMetricLaplace,\
					           gfs_metric_laplace_class ())
#define GFS_IS_METRIC_LAPLACE(obj)         (gts_object_is_from_class (obj,\
						   gfs_metric_laplace_class ()))

GfsVariableClass * gfs_metric_laplace_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __METRIC_H__ */
