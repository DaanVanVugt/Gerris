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

#ifndef __GRAPHIC_H__
#define __GRAPHIC_H__

#include "simulation.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void               gfs_write_gts               (GfsDomain * domain, 
						GfsVariable * v, 
						FttTraverseFlags flags,
						gint level,
						GtsBBox * box,
						FILE * fp);
GtsSurface *       gfs_isosurface              (GfsDomain * domain, 
						GfsVariable * v, 
						gdouble val,
						gint level);
void               gfs_write_gnuplot           (GfsDomain * domain, 
						GfsVariable * v, 
						FttTraverseFlags flags,
						gint level,
						GtsBBox * bbox,
						FILE * fp);
void               gfs_write_ppm               (GfsDomain * domain, 
						GfsFunction * condition,
						GfsVariable * v, 
						gdouble min, 
						gdouble max,
						FttTraverseFlags flags,
						gint level,
						FILE * fp,
						gboolean parallel);
void               gfs_write_grd               (GfsSimulation * sim, 
						GfsFunction * condition,
						GfsVariable * v,
						FttTraverseFlags flags,
						gint level,
						FILE * fp,
						gboolean parallel,
						gboolean interpolate);
gint               gfs_combine_ppm             (gchar ** fname, 
						guint nname, 
						FILE * fp);
void               gfs_write_squares           (GfsDomain * domain, 
						GfsVariable * v, 
						gdouble min, 
						gdouble max,
						FttTraverseFlags flags,
						gint level,
						GtsBBox * bbox,
						FILE * fp);
void               gfs_write_mac_velocity      (GfsDomain * domain,
						gdouble scale,
						FttTraverseFlags flags,
						gint level,
						GtsBBox * bbox,
						FILE * fp);
void               gfs_draw_cells              (FttCell * cell,
						FttTraverseFlags flags,
						gint level,
					       FILE * fp);
void               gfs_draw_boundary_conditions (GfsDomain * domain, 
						 FILE * fp);
void               gfs_draw_solid_boundaries   (GfsDomain * domain, 
						FILE * fp);
void               gfs_draw_refined_boundaries (GfsDomain * domain, 
						FILE * fp);
void               gfs_draw_levels             (FttCell * cell, 
						FILE * fp);
void               gfs_draw_surface            (GfsDomain * domain,
						GtsSurface * s,
						GfsVariable * v, 
						gdouble min, 
						gdouble max,
						FILE * fp);
void               gfs_extrude_profile         (GtsSurface * s,
						GSList * profile,
						gboolean closed,
						GList * path);
GList *            gfs_streamline_new          (GfsDomain * domain,
						GfsVariable ** U,
						FttVector p,
						GfsVariable * var,
						gdouble min,
						gdouble max,
						gboolean twist,
						gboolean (* stop) (FttCell *, 
								   GList *, 
								   gpointer),
						gpointer data);
void               gfs_streamline_write        (GList * stream, 
						FILE * fp);
GList *            gfs_streamline_read         (GtsFile * fp);
void               gfs_streamline_draw         (GList * stream, 
						FILE * fp);
void               gfs_streamline_destroy      (GList * stream);
void               gfs_draw_stream_cylinder    (GfsDomain * domain,
						FttVector p, 
						gdouble radius,
						GfsVariable * var, 
						gdouble min, 
						gdouble max,
						FILE * fp);
void               gfs_draw_stream_ribbon      (GfsDomain * domain,
						FttVector p, 
						gdouble radius,
						GfsVariable * var, 
						gdouble min, 
						gdouble max,
						FILE * fp);
void               gfs_draw_streamline         (GfsDomain * domain,
						FttVector p,
						FILE * fp);
gboolean           gfs_plane_cuts_cell         (FttVector plane[3], 
						FttCell * cell);
guint              gfs_cut_cube_vertices       (FttCell * cell, 
						gint maxlevel,
						FttVector * p, FttVector * n,
						FttVector v[12], FttDirection d[12],
						GfsVariable * var,
						gdouble val[12]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __GRAPHIC_H__ */
