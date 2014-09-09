/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009-2012 National Institute of Water and Atmospheric Research
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

#include "particle.h"

/* GfsParticulate: Header */

typedef struct _GfsParticulate GfsParticulate;

struct _GfsParticulate {
  GfsParticle parent;
  FttVector vel;
  gdouble mass, volume;
  FttVector force;
  GtsSListContainer * forces;
};

#define GFS_PARTICULATE(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsParticulate, gfs_particulate_class ())
#define GFS_IS_PARTICULATE(obj)         (gts_object_is_from_class (obj, gfs_particulate_class ()))

GfsEventClass * gfs_particulate_class  (void);

/* GfsParticleList: Header */

typedef struct _GfsParticleList GfsParticleList;

struct _GfsParticleList {
  GfsEventList parent;
  gint idlast;
  GtsSListContainer * forces;
};

#define GFS_PARTICLE_LIST(obj)            GTS_OBJECT_CAST (obj,		\
							   GfsParticleList, \
							   gfs_particle_list_class ())

#define GFS_IS_PARTICLE_LIST(obj)         (gts_object_is_from_class (obj, \
								     gfs_particle_list_class ()))

GfsEventClass * gfs_particle_list_class  (void);

/* GfsParticleForce: header */

typedef struct _GfsParticleForce GfsParticleForce;

struct _GfsParticleForce{
  GtsSListContainee parent;
  FttVector (* force) (GfsParticle *p, GfsParticleForce *force);
};

#define GFS_PARTICLE_FORCE(obj)            GTS_OBJECT_CAST (obj,		\
							GfsParticleForce, \
							gfs_particle_force_class ())
#define GFS_IS_PARTICLE_FORCE(obj)         (gts_object_is_from_class (obj, \
								      gfs_particle_force_class ()))

GtsSListContaineeClass * gfs_particle_force_class  (void);

/* GfsForceCoeff: header */

typedef struct _GfsForceCoeff GfsForceCoeff;

struct _GfsForceCoeff{
  GfsParticleForce parent;
  GfsFunction * coefficient;
  GfsVariable *re_p, *u_rel, *v_rel, *w_rel, *pdia;
  GfsParticulate *p;
};

#define FORCE_COEFF(obj)            GTS_OBJECT_CAST (obj,		\
						    GfsForceCoeff,		\
						    gfs_force_coeff_class ())
#define GFS_IS_FORCE_COEFF(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_coeff_class ()))
GtsSListContaineeClass * gfs_force_coeff_class  (void);

/* GfsForceLift: header */

#define GFS_IS_FORCE_LIFT(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_lift_class ()))
GtsSListContaineeClass * gfs_force_lift_class  (void);

/* GfsForceDrag: header */

#define GFS_IS_FORCE_DRAG(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_drag_class ()))
GtsSListContaineeClass * gfs_force_drag_class  (void);

/* GfsForceBuoy: header */

#define GFS_IS_FORCE_BUOY(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_buoy_class ()))
GtsSListContaineeClass * gfs_force_buoy_class  (void);

/* GfsDropletToParticle: header */

typedef struct _GfsDropletToParticle GfsDropletToParticle;

struct _GfsDropletToParticle{
  /*< private >*/
  GfsParticleList parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble resetwith;
  gdouble density;
};

#define DROPLET_TO_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDropletToParticle,\
					         gfs_droplet_to_particle_class ())

#define IS_DROPLET_TO_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_droplet_to_particle_class ()))

GfsEventClass * gfs_droplet_to_particle_class  (void);

/* GfsParticulateField: header */

typedef struct _GfsParticulateField                GfsParticulateField;

struct _GfsParticulateField {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsParticleList * plist;
  void (* voidfraction_func) (FttCell *, 
			      GfsVariable *, 
			      GfsParticulate *);
};
 

#define GFS_PARTICULATE_FIELD(obj)            GTS_OBJECT_CAST (obj,\
                                                   GfsParticulateField,\
                                                   gfs_particulate_field_class ())
#define GFS_IS_PARTICULATE_FIELD(obj)         (gts_object_is_from_class (obj,\
                                                   gfs_particulate_field_class ()))

GfsVariableClass * gfs_particulate_field_class  (void);

/* GfsFeedParticle: header */

typedef struct _GfsFeedParticle GfsFeedParticle;

struct _GfsFeedParticle{
  /*< private >*/
  GfsParticleList parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * posx, * posy, * posz;
  GfsFunction * velx, * vely, * velz;
  GfsFunction * np,   * mass, * vol;
};

#define GFS_FEED_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsFeedParticle,\
					         gfs_feed_particle_class ())

#define GFS_IS_FEED_PARTICLE(obj)         (gts_object_is_from_class (obj,\
					         gfs_feed_particle_class ()))

GfsEventClass * gfs_feed_particle_class  (void);


