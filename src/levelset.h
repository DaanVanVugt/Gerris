/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2006 National Institute of Water and Atmospheric Research
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

#ifndef __LEVELSET_H__
#define __LEVELSET_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "variable.h"

/* GfsVariableDistance: header */

typedef struct _GfsVariableDistance                GfsVariableDistance;

struct _GfsVariableDistance {
  /*< private >*/
  GfsVariable parent;
  gboolean first_done;

  /*< public >*/
  GfsVariable * v;
  gboolean stencil;
};

#define GFS_VARIABLE_DISTANCE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableDistance,\
					           gfs_variable_distance_class ())
#define GFS_IS_VARIABLE_DISTANCE(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_distance_class ()))

GfsVariableClass * gfs_variable_distance_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __LEVELSET_H__ */
