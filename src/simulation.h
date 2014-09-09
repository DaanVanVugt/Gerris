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

#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "domain.h"
#include "timestep.h"

#ifndef __EVENT_H__
  typedef struct _GfsSimulation         GfsSimulation;
#endif
typedef struct _GfsSimulationClass    GfsSimulationClass;
typedef struct _GfsTime               GfsTime;
typedef struct _GfsPhysicalParams     GfsPhysicalParams;
typedef struct _GfsAdaptStats         GfsAdaptStats;

struct _GfsTime {
  gdouble t, start, end;
  guint i, istart, iend;
  gdouble dtmax;
};

struct _GfsPhysicalParams {
  gdouble L, g;
  GfsFunction * alpha;
};

struct _GfsAdaptStats {
  guint removed, created;
  GtsRange cmax;
  GtsRange ncells;
  gint depth_increase;
};

struct _GfsSimulation {
  GfsDomain parent;

  GfsTime time;
  GfsPhysicalParams physical_params;

  GfsMultilevelParams projection_params;
  GfsMultilevelParams approx_projection_params;

  GfsAdvectionParams advection_params;

  GtsSListContainer * refines;

  GtsSListContainer * adapts;
  GfsAdaptStats adapts_stats;

  GtsSListContainer * events, * maps;
  GSList * modules, * globals, * preloaded_modules;

  GtsSListContainer * solids;
  guint thin;
  gboolean output_solid;

  gdouble tnext;

  GfsVariable * u0[FTT_DIMENSION];
};

struct _GfsSimulationClass {
  GfsDomainClass parent_class;

  void    (* run) (GfsSimulation *);
  gdouble (* cfl) (GfsSimulation *);
};

#define GFS_SIMULATION(obj)            GTS_OBJECT_CAST (obj,\
					           GfsSimulation,\
					           gfs_simulation_class ())
#define GFS_SIMULATION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsSimulationClass,\
						   gfs_simulation_class())
#define GFS_IS_SIMULATION(obj)         (gts_object_is_from_class (obj,\
						   gfs_simulation_class ()))

void                 gfs_advance_tracers         (GfsSimulation * sim, 
						  gdouble dt);
GfsSimulationClass * gfs_simulation_class        (void);
GfsSimulation *      gfs_simulation_new          (GfsSimulationClass * klass);
void                 gfs_simulation_init         (GfsSimulation * sim);
void                 gfs_simulation_write        (GfsSimulation * sim,
						  gint max_depth,  
						  FILE * fp);
void                 gfs_simulation_union_write  (GfsSimulation * sim,
						  gint max_depth,  
						  FILE * fp);
GfsSimulation *      gfs_simulation_read         (GtsFile * fp);
GSList *             gfs_simulation_get_solids   (GfsSimulation * sim);
guint                gfs_check_solid_fractions   (GfsDomain * domain);
void                 gfs_simulation_refine       (GfsSimulation * sim);
void                 gfs_simulation_set_timestep (GfsSimulation * sim);
void                 gfs_simulation_map          (GfsSimulation * sim, 
						  FttVector * p);
void                 gfs_simulation_map_inverse  (GfsSimulation * sim, 
						  FttVector * p);
void                 gfs_simulation_map_vector   (GfsSimulation * sim, 
						  const FttVector * p,
						  FttVector * v);
void                 gfs_simulation_map_inverse_vector  (GfsSimulation * sim, 
							 const FttVector * p,
							 FttVector * v);
void                 gfs_simulation_map_inverse_cell (GfsSimulation * sim, 
						      FttVector p[4]);
gdouble              gfs_dimensional_value       (GfsVariable * v, 
						  gdouble val);
gboolean             gfs_variable_is_dimensional (GfsVariable * v);
void                 gfs_time_init               (GfsTime * t);
void                 gfs_time_write              (GfsTime * t, 
						  FILE * fp);
void                 gfs_time_read               (GfsTime * t, 
						  GtsFile * fp);
void                 gfs_physical_params_init    (GfsPhysicalParams * p);
void                 gfs_physical_params_write   (GfsPhysicalParams * p, 
						  FILE * fp);
void                 gfs_physical_params_read    (GfsPhysicalParams * p,
						  GfsDomain * domain,
						  GtsFile * fp);
void                 gfs_simulation_run          (GfsSimulation * sim);
#define              gfs_object_simulation(o)     GFS_SIMULATION(GTS_OBJECT (o)->reserved)
#define              gfs_object_simulation_set(o,s) (GTS_OBJECT (o)->reserved = (s))

/* GfsAdvection: Header */

GfsSimulationClass * gfs_advection_class          (void);

/* GfsPoisson: Header */

#define GFS_IS_POISSON(obj)      (gts_object_is_from_class (obj, gfs_poisson_class ()))

GfsSimulationClass * gfs_poisson_class            (void);

/* GfsAxi: Header */

#define GFS_IS_AXI(obj)          (gts_object_is_from_class (obj, gfs_axi_class ()))

GfsSimulationClass * gfs_axi_class                (void);

/* GfsAdvectionAxi: Header */

#define GFS_IS_ADVECTION_AXI(obj)      (gts_object_is_from_class (obj,\
					           gfs_advection_axi_class ()))

GfsSimulationClass * gfs_advection_axi_class                (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SIMULATION_H__ */
