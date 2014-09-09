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

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include "boundary.h"
#include "surface.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* GfsLocateArray: Header */

typedef struct {
  GSList ** root;
  gdouble h, min[FTT_DIMENSION], max[FTT_DIMENSION];
  gint n[FTT_DIMENSION], size;
} GfsLocateArray;

GfsLocateArray * gfs_locate_array_new     (GfsDomain * domain);
GSList *         gfs_locate_array_locate  (GfsLocateArray * a, 
					   FttVector * p);
void             gfs_locate_array_destroy (GfsLocateArray * a);

/* GfsDomain: Header */

typedef struct _GfsDomainClass     GfsDomainClass;
typedef struct _GfsDiffusion       GfsDiffusion;
typedef struct _GfsSourceDiffusion GfsSourceDiffusion;
typedef struct _GfsTimer           GfsTimer;

struct _GfsTimer {
  GtsRange r;
  gdouble start;
};

struct _GfsDomain {
  GtsWGraph parent;

  int pid;
  GfsClock * timer;
  GHashTable * timers;

  GtsRange timestep;
  GtsRange size;

  gboolean profile_bc;

  GtsRange mpi_messages;
  GtsRange mpi_wait;

  guint rootlevel;
  FttVector refpos;
  FttVector lambda;

  GArray * allocated;
  GSList * variables;
  GSList * derived_variables;

  GfsVariable * velocity[FTT_DIMENSION];

  GSList * variables_io;
  gboolean binary;
  gint max_depth_write;

  FttCellInitFunc cell_init;
  gpointer cell_init_data;

  gint version;

  gpointer array;

  gboolean overlap; /* whether to overlap MPI communications with computation */

  /* coordinate metrics */
  gpointer metric_data;
  gdouble (* face_metric)       (const GfsDomain *, const FttCellFace *);
  gdouble (* cell_metric)       (const GfsDomain *, const FttCell *);
  void    (* solid_metric)      (const GfsDomain *, const FttCell *, FttVector *);
  gdouble (* scale_metric)      (const GfsDomain *, const FttCell *, FttComponent);
  gdouble (* face_scale_metric) (const GfsDomain *, const FttCellFace *, FttComponent);
  gdouble (* viscous_metric_implicit)    (const GfsDomain * domain, 
					  FttCell * cell,
					  FttComponent component);
  gdouble (* viscous_metric_explicit)    (const GfsDomain * domain, 
					  FttCell * cell,
					  GfsVariable * v,
					  GfsDiffusion * d);
  void    (* advection_metric)           (const GfsDomain * domain, 
					  FttCell * cell,
					  FttComponent c1,
					  gdouble m[2]);

  /* Object hash table for (optional) object IDs */
  GHashTable * objects;

  /* total number of parallel processes */
  int np;

  /* real time */
  GTimer * clock;

  GPtrArray * sorted; /**< array of sorted boxes */
  gboolean dirty;     /**< whether the sorted array needs updating */

  GSList * projections; /**< list of GfsDomainProjection associated with this domain */

  gboolean has_rotated_bc; /**< whether the domain uses "rotated" edges */

  void (* traverse_layers) (GfsDomain *, FttCellTraverseFunc, gpointer);
};

struct _GfsDomainClass {
  GtsWGraphClass parent_class;

  void    (* post_read) (GfsDomain *, GtsFile * fp);
};

#define GFS_DOMAIN(obj)            GTS_OBJECT_CAST (obj,\
					           GfsDomain,\
					           gfs_domain_class ())
#define GFS_DOMAIN_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsDomainClass,\
						   gfs_domain_class())
#define GFS_IS_DOMAIN(obj)         (gts_object_is_from_class (obj,\
						   gfs_domain_class ()))

#define gfs_domain_variables_number(d) ((d)->allocated->len - 1)
#define gfs_domain_variables_size(d)   (sizeof (GfsStateVector) +\
					sizeof (gdouble)*(MAX ((d)->allocated->len, 1) - 1))
     
GfsDomainClass * gfs_domain_class          (void);
void         gfs_domain_cell_traverse         (GfsDomain * domain,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data);
#define gfs_domain_traverse_leaves(d,f,data)  (gfs_domain_cell_traverse(d, \
					    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, f,data))
FttCellTraverse * gfs_domain_cell_traverse_new (GfsDomain * domain,
						FttTraverseType order,
						FttTraverseFlags flags,
						gint max_depth);
void         gfs_domain_traverse_layers       (GfsDomain * domain,
					       FttCellTraverseFunc func,
					       gpointer data);
void         gfs_domain_cell_traverse_condition (GfsDomain * domain,
						 FttTraverseType order,
						 FttTraverseFlags flags,
						 gint max_depth,
						 FttCellTraverseFunc func,
						 gpointer data,
						 gboolean (* condition) (FttCell *, gpointer),
						 gpointer cdata);
void         gfs_domain_cell_traverse_box     (GfsDomain * domain,
					       GtsBBox * box,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data);
void         gfs_domain_cell_traverse_boundary (GfsDomain * domain,
					       FttDirection d,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data);
void         gfs_domain_traverse_mixed        (GfsDomain * domain,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       FttCellTraverseFunc func,
					       gpointer data);
void         gfs_domain_traverse_cut          (GfsDomain * domain,
					       GfsGenericSurface * s,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       FttCellTraverseCutFunc func,
					       gpointer data);
void         gfs_domain_traverse_cut_2D       (GfsDomain * domain,
					       GfsGenericSurface * s,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       FttCellTraverseCutFunc func,
					       gpointer data);
void         gfs_domain_face_traverse         (GfsDomain * domain,
					       FttComponent c,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttFaceTraverseFunc func,
					       gpointer data);
void         gfs_domain_bc                    (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth,
					       GfsVariable * v);
void         gfs_domain_copy_bc               (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth,
					       GfsVariable * v,
					       GfsVariable * v1);
void         gfs_domain_homogeneous_bc        (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth,
					       GfsVariable * ov,
					       GfsVariable * v);
void         gfs_domain_homogeneous_bc_stencil (GfsDomain * domain,
						FttTraverseFlags flags,
						gint max_depth,
						GfsVariable * ov,
						GfsVariable * v,
						GfsLinearProblem * lp);
void         gfs_traverse_and_homogeneous_bc  (GfsDomain * domain,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data,
					       GfsVariable * ov,
					       GfsVariable * v);
void         gfs_traverse_and_bc              (GfsDomain * domain,
					       FttTraverseType order,
					       FttTraverseFlags flags,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data,
					       GfsVariable * v,
					       GfsVariable * v1);
void         gfs_domain_face_bc               (GfsDomain * domain,
					       FttComponent c,
					       GfsVariable * v);
void         gfs_domain_match                 (GfsDomain * domain);
void         gfs_domain_forget_boundary       (GfsDomain * domain, 
					       GfsBoundary * boundary);
void         gfs_domain_surface_bc            (GfsDomain * domain,
					       GfsVariable * v);
guint        gfs_domain_depth                 (GfsDomain * domain);
GtsRange     gfs_domain_stats_variable        (GfsDomain * domain,
					       GfsVariable * v,
					       FttTraverseFlags flags,
					       gint max_depth,
					       gboolean (* condition) (FttCell *, gpointer),
					       gpointer cdata);
GtsRange     gfs_domain_stats_solid           (GfsDomain * domain);
void         gfs_domain_stats_merged          (GfsDomain * domain,
					       GtsRange * solid,
					       GtsRange * number);
void         gfs_domain_stats_balance         (GfsDomain * domain,
					       GtsRange * size,
					       GtsRange * boundary,
					       GtsRange * mpiwait);
GfsNorm      gfs_domain_norm_variable         (GfsDomain * domain,
					       GfsVariable * v,
					       GfsFunction * w,
					       FttTraverseFlags flags,
					       gint max_depth,
					       gboolean (* condition) (FttCell *, gpointer),
					       gpointer cdata);
GfsNorm      gfs_domain_norm_residual         (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth,
					       gdouble dt,
					       GfsVariable * res);
GfsVariable ** gfs_domain_velocity            (GfsDomain * domain);
GfsNorm      gfs_domain_norm_velocity         (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth);
GfsDomain *  gfs_domain_read                  (GtsFile * fp);
void         gfs_domain_split                 (GfsDomain * domain,
					       gboolean one_box_per_pe);
FttCell *    gfs_domain_locate                (GfsDomain * domain,
					       FttVector target,
					       gint max_depth,
					       GfsBox ** where);
FttCell *    gfs_domain_boundary_locate       (GfsDomain * domain,
					       FttVector target,
					       gint max_depth,
					       GtsObject ** where);
gdouble      gfs_domain_cell_point_distance2  (GfsDomain * domain,
					       GtsPoint * p,
					       gdouble (* distance2) (FttCell *, 
								      GtsPoint *, 
								      gpointer),
					       gpointer data,
					       FttCell ** closest);
void         gfs_domain_advect_point          (GfsDomain * domain, 
					       FttVector * p,
					       gdouble dt);
guint        gfs_domain_size                  (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth);
gdouble      gfs_domain_cfl                   (GfsDomain * domain,
					       FttTraverseFlags flags,
					       gint max_depth);
void         gfs_cell_init                    (FttCell * cell,
					       GfsDomain * domain);
void         gfs_cell_reinit                  (FttCell * cell, 
					       GfsDomain * domain);
void         gfs_cell_fine_init               (FttCell * cell,
					       GfsDomain * domain);
void         gfs_cell_copy                    (const FttCell * from, 
					       FttCell * to,
					       GfsDomain * domain);
void         gfs_cell_read                    (FttCell * cell, 
					       GtsFile * fp,
					       GfsDomain * domain);
void         gfs_cell_write                   (const FttCell * cell, 
					       FILE * fp,
					       GSList * variables);
void         gfs_cell_read_binary             (FttCell * cell, 
					       GtsFile * fp,
					       GfsDomain * domain);
void         gfs_cell_write_binary            (const FttCell * cell, 
					       FILE * fp,
					       GSList * variables);
guint        gfs_domain_alloc                 (GfsDomain * domain);
void         gfs_domain_free                  (GfsDomain * domain, 
					       guint i);
GfsVariable * gfs_domain_add_variable         (GfsDomain * domain, 
					       const gchar * name,
					       const gchar * description);
GfsVariable * gfs_domain_get_or_add_variable  (GfsDomain * domain,
					       const gchar * name,
					       const gchar * description);
void         gfs_domain_solid_force           (GfsDomain * domain, 
					       FttVector * pf,
					       FttVector * vf,
					       FttVector * pm,
					       FttVector * vm,
					       GfsFunction * weight);
guint        gfs_domain_tag_droplets          (GfsDomain * domain,
					       GfsVariable * c,
					       GfsVariable * tag);
void         gfs_domain_remove_droplets       (GfsDomain * domain,
					       GfsVariable * c,
					       GfsVariable * v,
					       gint min,
					       gdouble val);
void         gfs_domain_remove_ponds          (GfsDomain * domain, 
					       gint min,
					       FttCellCleanupFunc cleanup,
					       gpointer data);
void         gfs_domain_remove_specks         (GfsDomain * domain);
void         gfs_domain_timer_start           (GfsDomain * domain, 
					       const gchar * name);
void         gfs_domain_timer_stop            (GfsDomain * domain, 
					       const gchar * name);
typedef
void      (* FttCellCombineTraverseFunc)      (FttCell * cell1, 
					       FttCell * cell2, 
					       gpointer data);
void         gfs_domain_combine_traverse      (GfsDomain * domain1,
					       GfsDomain * domain2,
					       FttCellCombineTraverseFunc inside,
					       gpointer idata,
					       FttCellTraverseFunc outside,
					       gpointer odata);

typedef struct {
  gchar * name, * description;
  gpointer func, data;
} GfsDerivedVariableInfo;

GfsDerivedVariable * gfs_domain_add_derived_variable  (GfsDomain * domain, 
						       GfsDerivedVariableInfo info);
gboolean     gfs_domain_remove_derived_variable (GfsDomain * domain, 
						 const gchar * name);
void         gfs_domain_sum                     (GfsDomain * domain, 
						 FttDirection d, 
						 GfsFunction * f, 
						 GfsVariable * v);
void         gfs_domain_filter                  (GfsDomain * domain, 
						 GfsVariable * v,
						 GfsVariable * fv);
typedef struct _GfsRequest GfsRequest;
GfsRequest * gfs_send_objects                   (GSList * list,
						 int dest);
void         gfs_wait                           (GfsRequest * r);
GSList *     gfs_receive_objects                (GfsDomain * domain, 
						 int src);
GfsRequest * gfs_send_boxes                     (GfsDomain * domain, 
						 GSList * boxes, 
						 int dest);
GSList *     gfs_receive_boxes                  (GfsDomain * domain, 
						 int src);

/**
 * gfs_domain_face_fraction:
 * @domain; a #GfsDomain.
 * @face: a #FttCellFace.
 *
 * Returns: the surface fraction of @face taking into account any
 * orthogonal metric of @domain.
 */
static inline
gdouble gfs_domain_face_fraction (const GfsDomain * domain, const FttCellFace * face)
{
  gdouble f = GFS_FACE_FRACTION (face);
  if (domain->face_metric)
    f *= (* domain->face_metric) (domain, face);
  return f;
}

/**
 * gfs_domain_face_fraction_right:
 * @domain; a #GfsDomain.
 * @face: a #FttCellFace.
 *
 * Returns: the surface fraction "to the right" of @face taking into account any
 * orthogonal metric of @domain.
 */
static inline
gdouble gfs_domain_face_fraction_right (const GfsDomain * domain, const FttCellFace * face)
{
  gdouble f = GFS_FACE_FRACTION_RIGHT (face);
  if (domain->face_metric) {
    FttCellFace face1;
    face1.cell = face->neighbor;
    face1.d = FTT_OPPOSITE_DIRECTION (face->d);
    f *= (* domain->face_metric) (domain, &face1);
  }
  return f;
}

/**
 * gfs_domain_cell_fraction:
 * @domain; a #GfsDomain.
 * @cell: a #FttCell.
 *
 * Returns: the volume fraction of @cell taking into account any
 * orthogonal metric of @domain.
 */
static inline
gdouble gfs_domain_cell_fraction (const GfsDomain * domain, const FttCell * cell)
{
  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
  if (domain->cell_metric)
    a *= (* domain->cell_metric) (domain, cell);
  return a;
}

/**
 * gfs_domain_solid_metric:
 * @domain; a #GfsDomain.
 * @cell: a mixed #FttCell.
 * @m: a #FttVector.
 *
 * Fills @m with the components of the coordinate metric at the center
 * of area of the solid surface contained within @cell.
 */
static inline
void gfs_domain_solid_metric (const GfsDomain * domain, const FttCell * cell, FttVector * m)
{
  if (domain->solid_metric)
    (* domain->solid_metric) (domain, cell, m);
}

/**
 * gfs_domain_face_scale_metric:
 * @domain; a #GfsDomain.
 * @face: a #FttCellFace.
 * @c: a #FttComponent.
 *
 * Returns: the @c component of the scale factor of the metric at the
 * center of the face.
 */
static inline
gdouble gfs_domain_face_scale_metric (const GfsDomain * domain, 
				      const FttCellFace * face, 
				      FttComponent c)
{
  if (domain->face_scale_metric)
    return (* domain->face_scale_metric) (domain, face, c);
  return 1.;
}

/**
 * gfs_cell_volume:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain.
 *
 * Returns: the volume of @cell.
 */
static inline
gdouble gfs_cell_volume (const FttCell * cell, const GfsDomain * domain)
{
  double v = ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.);
  return domain->cell_metric ? (* domain->cell_metric) (domain, cell)*v : v;
}

GtsObject * gfs_object_from_name        (GfsDomain * domain, 
					 const gchar * name);

/* GfsDomainProjection: Header */

typedef struct _GfsDomainProjection GfsDomainProjection;

struct _GfsDomainProjection {
  /*< private >*/
  GfsDomain parent;
  
  /*< public >*/
  GfsDomain * domain;
  FttComponent c;
};

#define GFS_DOMAIN_PROJECTION(obj)            GTS_OBJECT_CAST (obj,\
					           GfsDomainProjection,\
					           gfs_domain_projection_class ())

GfsDomainClass *      gfs_domain_projection_class    (void);
GfsDomainProjection * gfs_domain_projection_new      (GfsDomain * domain,
						      FttComponent c);
void                  gfs_domain_projection_reshape  (GfsDomainProjection * proj);
typedef void       (* GfsProjectionTraverseFunc)     (FttCell * cell, 
						      FttCell * proj, 
						      gpointer data);
void                  gfs_domain_projection_traverse (GfsDomainProjection * domain,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      GfsProjectionTraverseFunc func,
						      gpointer data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __DOMAIN_H__ */
