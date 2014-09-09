/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
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

#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsSurfaceGenericBc        GfsSurfaceGenericBc;

#include "timestep.h"
#include "event.h"

/* GfsVariable: Header */

typedef void (* GfsVariableFineCoarseFunc) (FttCell * cell, GfsVariable * v);

struct _GfsVariable {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  guint i;
  FttComponent component;
  GfsVariable * vector[FTT_DIMENSION];
  gchar * name, * description;
  gboolean centered;
  GfsVariableFineCoarseFunc fine_coarse, coarse_fine;
  GtsContainer * sources;
  GfsSurfaceGenericBc * surface_bc;
  GfsBc * default_bc;
  GfsDomain * domain;
  FttCellCleanupFunc cleanup;
  gdouble units;
  GfsVariable * face[2][4];
  gboolean face_source; /* whether source terms should be evaluated at cell faces only */
  gdouble orientation;   /* orientation of rotated variable */
  gboolean even;         /* "eveness" of rotated variable */
};

typedef struct _GfsVariableClass    GfsVariableClass;

struct _GfsVariableClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
};

#define GFS_VARIABLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsVariable,\
					         gfs_variable_class ())
#define GFS_VARIABLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsVariableClass,\
						 gfs_variable_class())
#define GFS_IS_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_variable_class ()))
#define GFS_VALUE(cell,v)            ((&GFS_STATE (cell)->place_holder)[(v)->i])

GfsVariableClass *    gfs_variable_class            (void);
GfsVariable *         gfs_variable_new              (GfsVariableClass * klass,
						     GfsDomain * domain,
						     const gchar * name,
						     const gchar * description);
#define               gfs_temporary_variable(d)     (gfs_variable_new (gfs_variable_class (),\
                                                                      (d), NULL, NULL))
GfsVariable *         gfs_variable_from_name        (GSList * i,
						     const gchar * name);
GSList *              gfs_variables_from_list       (GSList * i,
						     gchar * list,
						     gchar ** error);
void                  gfs_variables_swap            (GfsVariable * v1, 
						     GfsVariable * v2);
void                  gfs_variable_set_vector       (GfsVariable ** v,
						     guint n);
void                  gfs_variable_set_tensor       (GfsVariable * t[2][2]);
GfsVariable *         gfs_variable_clone            (GfsVariable * v, 
						     gchar * name);

/* GfsVariableBoolean: header */

#define GFS_IS_VARIABLE_BOOLEAN(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_boolean_class ()))

GfsVariableClass * gfs_variable_boolean_class  (void);

/* GfsVariableTracer: header */

typedef struct _GfsVariableTracer                GfsVariableTracer;

struct _GfsVariableTracer {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsAdvectionParams advection;
};

#define GFS_VARIABLE_TRACER(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableTracer,\
					           gfs_variable_tracer_class ())
#define GFS_IS_VARIABLE_TRACER(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_tracer_class ()))

GfsVariableClass * gfs_variable_tracer_class  (void);

/* GfsVariableResidual: header */

#define GFS_IS_VARIABLE_RESIDUAL(obj)         (gts_object_is_from_class (obj,\
					       gfs_variable_residual_class ()))

GfsVariableClass * gfs_variable_residual_class  (void);

/* GfsVariableFiltered: header */

typedef struct _GfsVariableFiltered                GfsVariableFiltered;

struct _GfsVariableFiltered {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * v;
  guint niter;
};

#define GFS_VARIABLE_FILTERED(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableFiltered,\
					           gfs_variable_filtered_class ())
#define GFS_IS_VARIABLE_FILTERED(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_filtered_class ()))

GfsVariableClass * gfs_variable_filtered_class  (void);

/* GfsVariableDiagonal: Header */

GfsVariableClass * gfs_variable_diagonal_class  (void);

/* GfsVariableFunction: header */

typedef struct _GfsVariableFunction                GfsVariableFunction;

struct _GfsVariableFunction {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsFunction * f;
};

#define GFS_VARIABLE_FUNCTION(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableFunction,\
					           gfs_variable_function_class ())
#define GFS_IS_VARIABLE_FUNCTION(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_function_class ()))

GfsVariableClass * gfs_variable_function_class  (void);

#if FTT_2D

/* GfsVariableStreamFunction: header */

#define GFS_IS_VARIABLE_STREAM_FUNCTION(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_stream_function_class ()))

GfsVariableClass * gfs_variable_stream_function_class  (void);

#endif /* FTT_2D */

/* GfsVariableAge: header */

#define GFS_IS_VARIABLE_AGE(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_age_class ()))

GfsVariableClass * gfs_variable_age_class  (void);

/* GfsVariableAverage: header */

typedef struct _GfsVariableAverage                GfsVariableAverage;

struct _GfsVariableAverage {
  /*< private >*/
  GfsVariableFunction parent;

  /*< public >*/
  FttComponent c;
};

#define GFS_VARIABLE_AVERAGE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableAverage,\
					           gfs_variable_average_class ())
#define GFS_IS_VARIABLE_AVERAGE(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_average_class ()))

GfsVariableClass * gfs_variable_average_class  (void);

/* GfsDerivedVariable: Header */

struct _GfsDerivedVariable {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  gchar * name, * description;
  gpointer func, data;
};

#define GFS_DERIVED_VARIABLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDerivedVariable,\
					         gfs_derived_variable_class ())
#define GFS_IS_DERIVED_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_derived_variable_class ()))

GtsObjectClass *     gfs_derived_variable_class            (void);
GfsDerivedVariable * gfs_derived_variable_from_name        (GSList * i, 
							    const gchar * name);

/* GfsConstant: Header */

typedef struct _GfsConstant GfsConstant;

/** \instance{GfsConstant} */
struct _GfsConstant {
  /*< private >*/
  GfsEvent parent;
  GfsDerivedVariable * derived;

  /*< public >*/
  gdouble val; /**< the value of the constant */
};

#define GFS_CONSTANT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsConstant,\
					         gfs_constant_class ())
#define GFS_IS_CONSTANT(obj)         (gts_object_is_from_class (obj,\
						 gfs_constant_class ()))

GfsEventClass *     gfs_constant_class            (void);

/* GfsSpatialSum: Header */

typedef struct _GfsSpatialSum         GfsSpatialSum;

/** \instance{GfsSpatialSum} */
struct _GfsSpatialSum {
  /*< private >*/
  GfsConstant parent;

  /*< public >*/
  GfsFunction * v; /**< the function to sum */
};

#define GFS_SPATIAL_SUM(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSpatialSum,\
					         gfs_spatial_sum_class ())
#define GFS_IS_SPATIAL_SUM(obj)         (gts_object_is_from_class (obj,\
						 gfs_spatial_sum_class ()))

GfsEventClass * gfs_spatial_sum_class  (void);

/* GfsVariablePoisson: Header */

typedef struct _GfsVariablePoisson                GfsVariablePoisson;

struct _GfsVariablePoisson {
  /*< private >*/
  GfsVariableFunction parent;

  /*< public >*/
  GfsMultilevelParams par;
};

#define GFS_VARIABLE_POISSON(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariablePoisson,\
					           gfs_variable_poisson_class ())
#define GFS_IS_VARIABLE_POISSON(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_poisson_class ()))

GfsVariableClass * gfs_variable_poisson_class  (void);

/* GfsVariableLaplacian: Header */

#define GFS_IS_VARIABLE_LAPLACIAN(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_laplacian_class ()))

GfsVariableClass * gfs_variable_laplacian_class  (void);

/* GfsHydrostaticPressure: Header */

typedef struct _GfsHydrostaticPressure GfsHydrostaticPressure;

/** \instance{GfsHydrostaticPressure} */
struct _GfsHydrostaticPressure {
  /*< private >*/
  GfsVariable parent;
  FttComponent c;
};

#define GFS_HYDROSTATIC_PRESSURE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsHydrostaticPressure,\
					         gfs_hydrostatic_pressure_class ())
#define GFS_IS_HYDROSTATIC_PRESSURE(obj)         (gts_object_is_from_class (obj,\
						 gfs_hydrostatic_pressure_class ()))

GfsVariableClass * gfs_hydrostatic_pressure_class  (void);
void               gfs_hydrostatic_pressure_update (GfsHydrostaticPressure * p, 
						    GfsFunction * alpha);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __VARIABLE_H__ */
