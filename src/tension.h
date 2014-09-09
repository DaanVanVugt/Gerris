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

#ifndef __TENSION_H__
#define __TENSION_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "source.h"

/* GfsSourceTensionGeneric: Header */

typedef struct _GfsSourceTensionGeneric         GfsSourceTensionGeneric;

struct _GfsSourceTensionGeneric {
  /*< private >*/
  GfsSourceVelocity parent;
  
  /*< public >*/
  GfsVariable * c;
  GfsFunction * sigma;
};

#define GFS_SOURCE_TENSION_GENERIC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceTensionGeneric,\
					         gfs_source_tension_generic_class ())
#define GFS_IS_SOURCE_TENSION_GENERIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_tension_generic_class ()))

GfsSourceGenericClass * gfs_source_tension_generic_class (void);

/* GfsSourceTensionCSS: Header */

typedef struct _GfsSourceTensionCSS         GfsSourceTensionCSS;

struct _GfsSourceTensionCSS {
  /*< private >*/
  GfsSourceTensionGeneric parent;
  GfsVariable * g[3];
  
  /*< public >*/
  GfsVariable * t[FTT_DIMENSION];
};

#define GFS_SOURCE_TENSION_CSS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceTensionCSS,\
					         gfs_source_tension_css_class ())
#define GFS_IS_SOURCE_TENSION_CSS(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_tension_css_class ()))

GfsSourceGenericClass * gfs_source_tension_css_class (void);

/* GfsSourceTension: Header */

typedef struct _GfsSourceTension         GfsSourceTension;

struct _GfsSourceTension {
  /*< private >*/
  GfsSourceTensionGeneric parent;
  
  /*< public >*/
  GfsVariable * k;
};

#define GFS_SOURCE_TENSION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceTension,\
					         gfs_source_tension_class ())
#define GFS_IS_SOURCE_TENSION(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_tension_class ()))

GfsSourceGenericClass * gfs_source_tension_class        (void);
void                    gfs_source_tension_coefficients (GfsSourceTension * s,
							 GfsDomain * domain,
							 GfsFunction * alpha);

/* GfsVariableCurvature: header */

typedef struct _GfsVariableCurvature                GfsVariableCurvature;

struct _GfsVariableCurvature {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * f, * kmax;
};

#define GFS_VARIABLE_CURVATURE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableCurvature,\
					           gfs_variable_curvature_class ())
#define GFS_IS_VARIABLE_CURVATURE(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_curvature_class ()))

GfsVariableClass * gfs_variable_curvature_class  (void);

/* GfsVariablePosition: header */

typedef struct _GfsVariablePosition                GfsVariablePosition;

struct _GfsVariablePosition {
  /*< private >*/
  GfsVariableCurvature parent;

  /*< public >*/
  FttComponent c;
  gdouble ref;
};

#define GFS_VARIABLE_POSITION(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariablePosition,\
					           gfs_variable_position_class ())
#define GFS_IS_VARIABLE_POSITION(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_position_class ()))

GfsVariableClass * gfs_variable_position_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __TENSION_H__ */
