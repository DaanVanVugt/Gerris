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

#ifndef __MPI_BOUNDARY_H__
#define __MPI_BOUNDARY_H__

#include "boundary.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsBoundaryMpi         GfsBoundaryMpi;

#ifdef HAVE_CONFIG_H
#  include "config.h"
#  ifdef HAVE_MPI
#    include <mpi.h>
#  endif

struct _GfsBoundaryMpi {
  /*< private >*/
  GfsBoundaryPeriodic parent;
  gint process, id;

#ifdef HAVE_MPI
  MPI_Comm comm;
  MPI_Request request[2];
  guint nrequest;
#endif /* HAVE_MPI */
};

#endif /* HAVE_CONFIG_H */

#define GFS_BOUNDARY_MPI(obj)            GTS_OBJECT_CAST (obj,\
					           GfsBoundaryMpi,\
					           gfs_boundary_mpi_class ())
#define GFS_IS_BOUNDARY_MPI(obj)         (gts_object_is_from_class (obj,\
						   gfs_boundary_mpi_class ()))
     
GfsBoundaryClass *    gfs_boundary_mpi_class    (void);
GfsBoundaryMpi *      gfs_boundary_mpi_new      (GfsBoundaryClass * klass,
						 GfsBox * box,
						 FttDirection d,
						 gint process,
						 gint id);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MPI_BOUNDARY_H__ */
