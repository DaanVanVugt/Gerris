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

#ifndef __WAVE_H__
#define __WAVE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "simulation.h"

/* GfsWave: Header */

#define GFS_WAVE_GAMMA 1.1
#define GFS_WAVE_F0 0.04

typedef struct _GfsWave    GfsWave;

struct _GfsWave {
  /*< private >*/
  GfsSimulation parent;
  guint ik, ith;
  void (* source) (GfsWave * wave);

  /*< public >*/
  guint nk, ntheta;
  gdouble alpha_s;
  GfsVariable *** F;
};

#define GFS_WAVE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsWave,\
					           gfs_wave_class ())
#define GFS_IS_WAVE(obj)         (gts_object_is_from_class (obj,\
							    gfs_wave_class ()))

GfsSimulationClass * gfs_wave_class        (void);

/* GfsInitWave: Header */

typedef struct _GfsInitWave         GfsInitWave;

struct _GfsInitWave {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  GfsFunction * d, * hs;
};

#define GFS_INIT_WAVE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitWave,\
					         gfs_init_wave_class ())
#define GFS_IS_INIT_WAVE(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_wave_class ()))

GfsGenericInitClass * gfs_init_wave_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __WAVE_H__ */
