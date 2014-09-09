/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2009 National Institute of Water and Atmospheric Research
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

#ifndef __BALANCE_H__
#define __BALANCE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "event.h"

/* GfsEventBalance: Header */

typedef struct _GfsEventBalance         GfsEventBalance;

struct _GfsEventBalance {
  GfsEvent parent;

  gdouble max;
};

#define GFS_EVENT_BALANCE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsEventBalance,\
					         gfs_event_balance_class ())
#define GFS_IS_EVENT_BALANCE(obj)         (gts_object_is_from_class (obj,\
						 gfs_event_balance_class ()))

GfsEventClass * gfs_event_balance_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BALANCE_H__ */
