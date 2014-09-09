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

#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include "fluid.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsBox                    GfsBox;
typedef struct _GfsBoxClass               GfsBoxClass;
typedef struct _GfsBoundary               GfsBoundary;
typedef struct _GfsBoundaryClass          GfsBoundaryClass;

/* GfsBc: Header */

typedef struct _GfsBc         GfsBc;

struct _GfsBc {
  /*< private >*/
  GtsObject parent;
  GfsLinearProblem * lp;

  /*< public >*/
  GfsBoundary * b;
  GfsVariable * v;
  gboolean extra;
  FttFaceTraverseFunc bc, homogeneous_bc;
  FttFaceTraverseFunc homogeneous_bc_stencil;
  FttFaceTraverseFunc face_bc;
};

typedef struct _GfsBcClass    GfsBcClass;

struct _GfsBcClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/  
};

#define GFS_BC(obj)                    GTS_OBJECT_CAST (obj,\
				                        GfsBc,\
				                        gfs_bc_class ())
#define GFS_BC_CLASS(klass)            GTS_OBJECT_CLASS_CAST (klass,\
						        GfsBcClass,\
						        gfs_bc_class())
#define GFS_IS_BC(obj)                 (gts_object_is_from_class (obj,\
					gfs_bc_class ()))

GfsBcClass * gfs_bc_class  (void);
GfsBc *      gfs_bc_new    (GfsBcClass * k, 
			    GfsVariable * v, 
			    gboolean extra);

/* GfsBcValue: Header */

typedef struct _GfsBcValue         GfsBcValue;

struct _GfsBcValue {
  /*< private >*/
  GfsBc parent;

  /*< public >*/
  GfsFunction * val;
};

#define GFS_BC_VALUE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsBcValue,\
					         gfs_bc_value_class ())
#define GFS_IS_BC_VALUE(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_value_class ()))

GfsBcClass * gfs_bc_value_class  (void);

/* GfsBcDirichlet: Header */

#define GFS_IS_BC_DIRICHLET(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_dirichlet_class ()))

GfsBcClass * gfs_bc_dirichlet_class  (void);

/* GfsBcNeumann: Header */

#define GFS_IS_BC_NEUMANN(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_neumann_class ()))

GfsBcClass * gfs_bc_neumann_class  (void);

/* GfsBcAngle: Header */

#define GFS_IS_BC_ANGLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_angle_class ()))

GfsBcClass * gfs_bc_angle_class  (void);

/* GfsBcNavier: Header */

typedef struct _GfsBcNavier         GfsBcNavier;

struct _GfsBcNavier {
  /*< private >*/
  GfsBcValue parent;

  /*< public >*/
  GfsFunction * lambda;
};

#define GFS_BC_NAVIER(obj)            GTS_OBJECT_CAST (obj,\
					         GfsBcNavier,\
					         gfs_bc_navier_class ())
#define GFS_IS_BC_NAVIER(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_navier_class ()))

GfsBcClass * gfs_bc_navier_class  (void);

/* GfsBoundary: Header */

typedef enum {
  GFS_BOUNDARY_CENTER_VARIABLE,
  GFS_BOUNDARY_FACE_VARIABLE,
  GFS_BOUNDARY_MATCH_VARIABLE,
  GFS_BOUNDARY_VARIABLE_NUMBER
} GfsBoundaryVariableType;

struct _GfsBoundary {
  /*< private >*/
  GtsObject parent;

  FttCell * root;
  GfsBox * box;
  FttDirection d;
  guint depth;
  GfsBc * default_bc;
  gboolean changed;

  /*< public >*/
  GfsVariable * v;
  GfsBoundaryVariableType type;
  GHashTable * bc;
};

struct _GfsBoundaryClass {
  GtsObjectClass parent_class;

  void (* match)             (GfsBoundary * boundary);
  void (* send)              (GfsBoundary * boundary);
  void (* receive)           (GfsBoundary * boundary,
			      FttTraverseFlags flags,
			      gint max_depth);
  void (* synchronize)       (GfsBoundary * boundary);
  void (* update)            (GfsBoundary * boundary);
};

#define GFS_BOUNDARY(obj)            GTS_OBJECT_CAST (obj,\
					           GfsBoundary,\
					           gfs_boundary_class ())
#define GFS_BOUNDARY_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsBoundaryClass,\
						   gfs_boundary_class())
#define GFS_IS_BOUNDARY(obj)         (gts_object_is_from_class (obj,\
						   gfs_boundary_class ()))
     
GfsBoundaryClass * gfs_boundary_class                (void);
GfsBoundary *      gfs_boundary_new                  (GfsBoundaryClass * klass,
						      GfsBox * box,
						      FttDirection d);
void               gfs_boundary_send                 (GfsBoundary * boundary);
void               gfs_boundary_update               (GfsBoundary * boundary);
void               gfs_boundary_receive              (GfsBoundary * boundary,
						      FttTraverseFlags flags,
						      gint max_depth);
void               gfs_boundary_synchronize          (GfsBoundary * boundary);
GfsBc *            gfs_boundary_lookup_bc            (GfsBoundary * b, 
						      GfsVariable * v);
void               gfs_boundary_set_default_bc       (GfsBoundary * b, 
						      GfsBc * bc);
void               gfs_variable_set_default_bc       (GfsVariable * v, 
						      GfsBc * bc);
void               gfs_boundary_add_bc               (GfsBoundary * b, 
						      GfsBc * bc);

/* GfsBoundaryInflowConstant: Header */

typedef struct _GfsBoundaryInflowConstant         GfsBoundaryInflowConstant;
typedef struct _GfsBoundaryInflowConstantClass    GfsBoundaryInflowConstantClass;

struct _GfsBoundaryInflowConstant {
  GfsBoundary parent;

  GfsFunction * un;
};

struct _GfsBoundaryInflowConstantClass {
  GfsBoundaryClass parent_class;
};

#define GFS_BOUNDARY_INFLOW_CONSTANT(obj)  GTS_OBJECT_CAST (obj,\
					 GfsBoundaryInflowConstant,\
					 gfs_boundary_inflow_constant_class ())
#define GFS_BOUNDARY_INFLOW_CONSTANT_CLASS(klass) GTS_OBJECT_CLASS_CAST (klass,\
					 GfsBoundaryInflowConstantClass,\
					 gfs_boundary_inflow_constant_class())
#define GFS_IS_BOUNDARY_INFLOW_CONSTANT(obj) (gts_object_is_from_class (obj,\
					 gfs_boundary_inflow_constant_class ()))
     
GfsBoundaryInflowConstantClass * gfs_boundary_inflow_constant_class  (void);

/* GfsBoundaryOutflow: Header */

typedef struct _GfsBoundaryOutflowClass    GfsBoundaryOutflowClass;

struct _GfsBoundaryOutflowClass {
  GfsBoundaryClass parent_class;
};

#define GFS_BOUNDARY_OUTFLOW_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsBoundaryOutflowClass,\
						   gfs_boundary_outflow_class())
#define GFS_IS_BOUNDARY_OUTFLOW(obj)         (gts_object_is_from_class (obj,\
					    gfs_boundary_outflow_class ()))
     
GfsBoundaryOutflowClass * gfs_boundary_outflow_class    (void);

/* GfsBoundaryGradient: Header */

#define GFS_IS_BOUNDARY_GRADIENT(obj)         (gts_object_is_from_class (obj,\
					    gfs_boundary_gradient_class ()))
     
GfsBoundaryClass * gfs_boundary_gradient_class    (void);

/* GfsBoundaryPeriodic: Header */

typedef struct _GfsBoundaryPeriodic         GfsBoundaryPeriodic;

struct _GfsBoundaryPeriodic {
  /*< private >*/
  GfsBoundary parent;

  GfsBox * matching;
  FttDirection d;
  GArray * sndbuf, * rcvbuf;
  guint sndcount, rcvcount;

  gdouble rotate;
};

#define GFS_BOUNDARY_PERIODIC(obj)            GTS_OBJECT_CAST (obj,\
					           GfsBoundaryPeriodic,\
					           gfs_boundary_periodic_class ())
#define GFS_IS_BOUNDARY_PERIODIC(obj)         (gts_object_is_from_class (obj,\
						   gfs_boundary_periodic_class ()))
     
GfsBoundaryClass *    gfs_boundary_periodic_class    (void);
GfsBoundaryPeriodic * gfs_boundary_periodic_new      (GfsBoundaryClass * klass,
						      GfsBox * box,
						      FttDirection d,
						      GfsBox * matching);
void                  gfs_boundary_periodic_rotate   (GfsBoundaryPeriodic * boundary,
						      FttDirection rotate,
						      gdouble orientation);
GfsBoundaryPeriodic * gfs_boundary_periodic_rotate_new (GfsBoundaryClass * klass,
							GfsBox * box,
							FttDirection d,
							GfsBox * matching,
							FttDirection rotate,
							gdouble orientation);

/* GfsGEdge: Header */
  
typedef struct _GfsGEdge         GfsGEdge;
typedef struct _GfsGEdgeClass    GfsGEdgeClass;

struct _GfsGEdge {
  GtsGEdge parent;

  FttDirection d, rotate;
};

struct _GfsGEdgeClass {
  GtsGEdgeClass parent_class;
};

#define GFS_GEDGE(obj)            GTS_OBJECT_CAST (obj,\
					          GfsGEdge,\
					          gfs_gedge_class ())
#define GFS_GEDGE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						  GfsGEdgeClass,\
						  gfs_gedge_class())
#define GFS_IS_GEDGE(obj)         (gts_object_is_from_class (obj,\
						  gfs_gedge_class ()))
     
GfsGEdgeClass * gfs_gedge_class          (void);
GfsGEdge *      gfs_gedge_new            (GfsGEdgeClass * klass,
					  GfsBox * b1, GfsBox * b2,
					  FttDirection d);
void            gfs_gedge_link_boxes     (GfsGEdge * edge);

struct _GfsBox {
  GtsGNode parent;

  FttCell * root;
  GtsObject * neighbor[FTT_NEIGHBORS];
  guint id;
  int pid;
  gint size;
};

struct _GfsBoxClass {
  GtsGNodeClass parent_class;
};

#define GFS_BOX(obj)            GTS_OBJECT_CAST (obj,\
					        GfsBox,\
					        gfs_box_class ())
#define GFS_BOX_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						      GfsBoxClass,\
						      gfs_box_class())
#define GFS_IS_BOX(obj)         (gts_object_is_from_class (obj,\
						          gfs_box_class ()))
     
GfsBoxClass *    gfs_box_class                (void);
GfsBox *         gfs_box_new                  (GfsBoxClass * klass);

static inline
GfsDomain * gfs_box_domain (GfsBox * box)
{
  GfsDomain * d;

  g_return_val_if_fail (box != NULL, NULL);
  
  d = GTS_OBJECT (box)->reserved;
  if (GTS_SLIST_CONTAINEE (box)->containers) {
    GSList * i = GTS_SLIST_CONTAINEE (box)->containers;

    while (i->next)
      i = i->next;
    d = i->data;
  }
  return d;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BOUNDARY_H__ */
