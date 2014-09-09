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

#ifndef __EVENT_H__
#define __EVENT_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <gts.h>

typedef struct _GfsEvent                GfsEvent;
typedef struct _GfsEventClass           GfsEventClass;

struct _GfsEvent {
  GtsSListContainee parent;

  gdouble t, start, end, step;
  guint i, istart, iend, istep;
  
  guint n;
  gboolean end_event, realised, redo;
  gchar * name;
};

typedef struct _GfsSimulation           GfsSimulation;

struct _GfsEventClass {
  GtsSListContaineeClass parent_class;

  gboolean (* event)      (GfsEvent * event, GfsSimulation * sim);
  void     (* post_event) (GfsEvent * event, GfsSimulation * sim);
  void     (* event_half) (GfsEvent * event, GfsSimulation * sim);
};

#include "simulation.h"

#define GFS_EVENT(obj)            GTS_OBJECT_CAST (obj,\
					           GfsEvent,\
					           gfs_event_class ())
#define GFS_EVENT_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsEventClass,\
						   gfs_event_class())
#define GFS_IS_EVENT(obj)         (gts_object_is_from_class (obj,\
						   gfs_event_class ()))

GfsEventClass * gfs_event_class       (void);
GfsEvent *      gfs_event_new         (GfsEventClass * klass);
void            gfs_event_set         (GfsEvent * e, 
				       gdouble start, 
				       gdouble end, 
				       gdouble step,
				       gint istart, 
				       gint iend, 
				       gint istep);
void            gfs_event_init        (GfsEvent * event,
				       GfsSimulation * sim);
void            gfs_event_do          (GfsEvent * event, 
				       GfsSimulation * sim);
void            gfs_event_redo        (GfsEvent * event, 
				       GfsSimulation * sim);
gdouble         gfs_event_next        (GfsEvent * event, 
				       GfsSimulation * sim);
void            gfs_event_half_do     (GfsEvent * event, 
				       GfsSimulation * sim);
#define         gfs_event_is_repetitive(e) ((e)->step < G_MAXDOUBLE || (e)->istep < G_MAXINT)

/* GfsGenericInit: Header */

typedef struct _GfsEvent      GfsGenericInit;
typedef struct _GfsEventClass GfsGenericInitClass;

#define GFS_IS_GENERIC_INIT(obj)         (gts_object_is_from_class (obj,\
						 gfs_generic_init_class ()))

GfsEventClass * gfs_generic_init_class         (void);

/* GfsInit: Header */

typedef struct _GfsInit         GfsInit;

struct _GfsInit {
  /*< private >*/
  GfsGenericInit parent;
  GSList * f;
};

#define GFS_INIT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInit,\
					         gfs_init_class ())
#define GFS_IS_INIT(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_class ()))

GfsGenericInitClass * gfs_init_class  (void);

/* GfsInitMask: Header */

typedef struct _GfsInitMask         GfsInitMask;

struct _GfsInitMask {
  /*< private >*/
  GfsGenericInit parent;
  GSList * masked_boxes;

  /*< public >*/
  GfsFunction * mask;
};

#define GFS_INIT_MASK(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitMask,\
					         gfs_init_mask_class ())
#define GFS_IS_INIT_MASK(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_mask_class ()))

GfsGenericInitClass * gfs_init_mask_class  (void);

/* GfsInitFlowConstant: Header */

GfsEventClass * gfs_init_flow_constant_class  (void);

/* GfsInitVorticity: Header */

typedef struct _GfsInitVorticity         GfsInitVorticity;

struct _GfsInitVorticity {
  /*< private >*/
  GfsGenericInit parent;
  GfsVariable * vort, ** u;
#if FTT_2D
  GfsVariable * stream;
#else
  GfsVariable * stream[3];
#endif

  /*< public >*/
  GfsFunction * f;
#if !FTT_2D
  GfsFunction * fv[3];
#endif
};

#define GFS_INIT_VORTICITY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitVorticity,\
					         gfs_init_vorticity_class ())
#define GFS_IS_INIT_VORTICITY(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_vorticity_class ()))

GfsGenericInitClass * gfs_init_vorticity_class  (void);

/* GfsEventSum: Header */

typedef struct _GfsEventSum         GfsEventSum;

struct _GfsEventSum {
  GfsEvent parent;

  GfsFunction * v;
  GfsVariable * sv;
  FttCellTraverseFunc sum;
  gdouble last, dt;
};

#define GFS_EVENT_SUM(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventSum,\
					         gfs_event_sum_class ())
#define GFS_IS_EVENT_SUM(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_sum_class ()))

GfsEventClass * gfs_event_sum_class  (void);

/* GfsEventSumDirection: Header */

typedef struct _GfsEventSumDirection         GfsEventSumDirection;

struct _GfsEventSumDirection {
  GfsEventSum parent;

  FttDirection d;
};

#define GFS_EVENT_SUM_DIRECTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventSumDirection,\
					         gfs_event_sum_direction_class ())
#define GFS_IS_EVENT_SUM_DIRECTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_sum_direction_class ()))

GfsEventClass * gfs_event_sum_direction_class  (void);

/* GfsEventHarmonic: Header */

typedef struct _GfsEventHarmonic         GfsEventHarmonic;

struct _GfsEventHarmonic {
  GfsEvent parent;

  GArray * omega;
  GfsVariable * v, * z, * e, ** A, ** B;
  gdouble * vsin, * vcos, ** M, ** iM, ** Mn, * x, * a;
  gchar * Aname, * Bname;
  gboolean invertible;
};

#define GFS_EVENT_HARMONIC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventHarmonic,\
					         gfs_event_harmonic_class ())
#define GFS_IS_EVENT_HARMONIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_harmonic_class ()))

GfsEventClass * gfs_event_harmonic_class  (void);

/* GfsEventStop: Header */

typedef struct _GfsEventStop         GfsEventStop;

struct _GfsEventStop {
  GfsEvent parent;

  GfsVariable * v, * oldv, * diff;
  gdouble last, max;
  gboolean relative;
};

#define GFS_EVENT_STOP(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventStop,\
					         gfs_event_stop_class ())
#define GFS_IS_EVENT_STOP(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_stop_class ()))

GfsEventClass * gfs_event_stop_class  (void);

/* GfsEventScript: Header */

FILE * gfs_popen (GfsSimulation * sim, 
		  const char * command, 
		  const char * type);

typedef struct _GfsEventScript         GfsEventScript;

struct _GfsEventScript {
  GfsEvent parent;

  gchar * script;
};

#define GFS_EVENT_SCRIPT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventScript,\
					         gfs_event_script_class ())
#define GFS_IS_EVENT_SCRIPT(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_script_class ()))
#define GFS_EVENT_SCRIPT_STOP            64

GfsEventClass * gfs_event_script_class  (void);

/* GfsInitFraction: Header */

typedef struct _GfsInitFraction         GfsInitFraction;

struct _GfsInitFraction {
  /*< private >*/
  GfsGenericInit parent;

  GfsVariable * c;
  GfsGenericSurface * surface;
};

typedef struct _GfsInitFractionClass    GfsInitFractionClass;

struct _GfsInitFractionClass {
  /*< private >*/
  GfsGenericInitClass parent_class;

  /*< public >*/
};

#define GFS_INIT_FRACTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitFraction,\
					         gfs_init_fraction_class ())
#define GFS_INIT_FRACTION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitFractionClass,\
						 gfs_init_fraction_class())
#define GFS_IS_INIT_FRACTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_fraction_class ()))

GfsInitFractionClass * gfs_init_fraction_class  (void);

/* GfsRemoveDroplets: Header */

typedef struct _GfsRemoveDroplets         GfsRemoveDroplets;

struct _GfsRemoveDroplets {
  /*< private >*/
  GfsEvent parent;
  GfsVariable * v;

  /*< public >*/
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble val;
};

#define GFS_REMOVE_DROPLETS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsRemoveDroplets,\
					         gfs_remove_droplets_class ())
#define GFS_IS_REMOVE_DROPLETS(obj)         (gts_object_is_from_class (obj,\
						 gfs_remove_droplets_class ()))

GfsEventClass * gfs_remove_droplets_class  (void);

/* GfsRemovePonds: Header */

typedef struct _GfsRemovePonds         GfsRemovePonds;

struct _GfsRemovePonds {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  gint min;
};

#define GFS_REMOVE_PONDS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsRemovePonds,\
					         gfs_remove_ponds_class ())
#define GFS_IS_REMOVE_PONDS(obj)         (gts_object_is_from_class (obj,\
						 gfs_remove_ponds_class ()))

GfsEventClass * gfs_remove_ponds_class  (void);

/* GfsEventFilter: Header */

typedef struct _GfsEventFilter         GfsEventFilter;

struct _GfsEventFilter {
  /*< private >*/
  GfsEvent parent;
  GfsVariable * tmp;

  /*< public >*/
  GfsVariable * v;
  gdouble scale;
};

#define GFS_EVENT_FILTER(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventFilter,\
					         gfs_event_filter_class ())
#define GFS_IS_EVENT_FILTER(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_filter_class ()))

GfsEventClass * gfs_event_filter_class  (void);

/* GfsEventList: Header */

typedef struct _GfsEventList         GfsEventList;

struct _GfsEventList {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GtsObjectClass * klass;
  GtsSListContainer * list;
};

#define GFS_EVENT_LIST(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventList,\
					         gfs_event_list_class ())
#define GFS_IS_EVENT_LIST(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_list_class ()))

GfsEventClass * gfs_event_list_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __EVENT_H__ */
