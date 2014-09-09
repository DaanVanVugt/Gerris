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

#include "particulatecommon.h"

/* Initialize module */

const gchar gfs_module_name[] = "particulates";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_particulate_class ();
  gfs_particle_list_class ();
  gfs_force_lift_class ();
  gfs_force_drag_class ();
  gfs_force_buoy_class ();
  gfs_particle_force_class ();

  gfs_droplet_to_particle_class ();
  gfs_feed_particle_class ();

  gfs_particulate_field_class ();

  return NULL; 
}
