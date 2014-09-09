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

#ifndef __REFINE_H__
#define __REFINE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "simulation.h"

#if FTT_2D
# define GFS_DIAGONAL 0.707106781187
#else /* 3D */
# define GFS_DIAGONAL 0.866025403785
#endif /* 3D */

/* GfsRefine: Header */

typedef struct _GfsRefine             GfsRefine;
typedef struct _GfsRefineClass        GfsRefineClass;

struct _GfsRefine {
  GtsSListContainee parent;

  GfsFunction * maxlevel;
};

struct _GfsRefineClass {
  GtsSListContaineeClass parent_class;

  void (* refine) (GfsRefine * refine, GfsSimulation * simulation);
};

#define GFS_REFINE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefine,\
					           gfs_refine_class ())
#define GFS_REFINE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsRefineClass,\
						   gfs_refine_class())
#define GFS_IS_REFINE(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_class ()))
     
GfsRefineClass * gfs_refine_class  (void);
GfsRefine *      gfs_refine_new    (GfsRefineClass * klass);

/* GfsRefineSolid: Header */

#define GFS_IS_REFINE_SOLID(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_solid_class ()))
     
GfsRefineClass * gfs_refine_solid_class  (void);

/* GfsRefineSurface: Header */

typedef struct _GfsRefineSurface         GfsRefineSurface;

struct _GfsRefineSurface {
  GfsRefine parent;

  GfsGenericSurface * surface;
};

#define GFS_REFINE_SURFACE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineSurface,\
					           gfs_refine_surface_class ())
#define GFS_IS_REFINE_SURFACE(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_surface_class ()))
     
GfsRefineClass * gfs_refine_surface_class  (void);

/* GfsRefineDistance: Header */

typedef struct _GfsRefineDistance         GfsRefineDistance;

struct _GfsRefineDistance {
  GfsRefineSurface parent;

  GNode * stree;
};

#define GFS_REFINE_DISTANCE(obj)            GTS_OBJECT_CAST (obj,\
					          GfsRefineDistance,\
					          gfs_refine_distance_class ())
#define GFS_IS_REFINE_DISTANCE(obj)         (gts_object_is_from_class (obj,\
						  gfs_refine_distance_class ()))
     
GfsRefineClass * gfs_refine_distance_class  (void);

/* GfsRefineHeight: Header */

#define GFS_IS_REFINE_HEIGHT(obj)         (gts_object_is_from_class (obj,\
						  gfs_refine_height_class ()))
     
GfsRefineClass * gfs_refine_height_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __REFINE_H__ */
