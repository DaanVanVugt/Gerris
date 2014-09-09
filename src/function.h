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

#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#define NODATA GFS_NODATA

static double Dirichlet = 1.;
static double Neumann = 0.;
static GfsSimulation * _sim = NULL;
static FttCell * _cell = NULL;

static double dd (const gchar * name, FttComponent c) {
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  return gfs_dimensional_value (v, gfs_center_gradient (_cell, c, v->i)/
				(_sim->physical_params.L*ftt_cell_size (_cell)));
}

static double dd2 (const gchar * name, FttComponent c) {
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  FttCellFace f1 = ftt_cell_face (_cell, 2*c);
  FttCellFace f2 = ftt_cell_face (_cell, 2*c + 1);
  if (f1.neighbor && f2.neighbor) {
    GfsGradient g1, g2;
    gfs_face_gradient (&f1, &g1, v->i, -1);
    gfs_face_gradient (&f2, &g2, v->i, -1);
    return gfs_dimensional_value (v, (g1.b + g2.b - (g1.a + g2.a)*GFS_VALUE (_cell, v))
				  /pow (_sim->physical_params.L*ftt_cell_size (_cell), 2.));
  }
  return 0.;
}

static double dx (const gchar * name) { return dd (name, FTT_X); }
static double dy (const gchar * name) { return dd (name, FTT_Y); }
static double dx2 (const gchar * name) { return dd2 (name, FTT_X); }
static double dy2 (const gchar * name) { return dd2 (name, FTT_Y); }
#if !FTT_2D
static double dz (const gchar * name) { return dd (name, FTT_Z); }
static double dz2 (const gchar * name) { return dd2 (name, FTT_Z); }
#endif /* 3D */

static double area (const gchar * name)
{
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL || !GFS_IS_VARIABLE_TRACER_VOF (v))
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  FttVector m, p;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (_cell, t->m[c]);
  return gfs_plane_area_center (&m, GFS_VALUE (_cell, t->alpha), &p)/
    (_sim->physical_params.L*ftt_cell_size (_cell));
}

static double correctness (const gchar * name)
{
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL || !GFS_IS_VARIABLE_TRACER_VOF (v))
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  return gfs_vof_correctness (_cell, GFS_VARIABLE_TRACER_VOF (v));
}

static double distance (double xo, double yo, double zo)
{
  /* fixme: this doesn't take mapping into account properly */
  GtsPoint o;
  o.x = xo; o.y = yo; o.z = zo;
  gfs_simulation_map (_sim, (FttVector *) &o.x);
  GtsBBox bb;
  ftt_cell_bbox (_cell, &bb);
  gdouble min, max;
  gts_bbox_point_distance2 (&bb, &o, &min, &max);
  return sqrt (min)*_sim->physical_params.L;
}

static gboolean is_velocity (GfsVariable * v, GfsDomain * domain)
{
  FttComponent c;
  GfsVariable ** u = gfs_domain_velocity (domain);

  for (c = 0; c < FTT_DIMENSION; c++)
    if (v == u[c])
      return TRUE;
  return FALSE;
}

static void dirichlet_bc (FttCell * cell)
{
  cell->flags |= GFS_FLAG_DIRICHLET;
  GFS_STATE (cell)->solid->fv = 0.;
}

static double dsd (const gchar * name, FttComponent c)
{
  g_return_val_if_fail (_cell != NULL, NODATA);
  if (!GFS_IS_MIXED (_cell))
    return NODATA;
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return NODATA;
  if (v->surface_bc)
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (_cell, v->surface_bc);
  else if (is_velocity (v, GFS_DOMAIN (_sim)))
    dirichlet_bc (_cell);
  else /* Neumann */
    return 0.;
  if ((_cell->flags & GFS_FLAG_DIRICHLET) == 0)
    return NODATA;
  FttVector g;
  gfs_cell_dirichlet_gradient (_cell, v->i, -1, GFS_STATE (_cell)->solid->fv, &g);
  return gfs_dimensional_value (v, (&g.x)[c]/(_sim->physical_params.L*ftt_cell_size (_cell)));
}

static double dsx (const gchar * name) { return dsd (name, FTT_X); }
static double dsy (const gchar * name) { return dsd (name, FTT_Y); }
#if !FTT_2D
static double dsz (const gchar * name) { return dsd (name, FTT_Z); }
#endif /* 3D */

static double flux (const gchar * name)
{
  g_return_val_if_fail (_cell != NULL, NODATA);
  if (!GFS_IS_MIXED (_cell))
    return 0.;
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return 0.;
  if (v->surface_bc)
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (_cell, v->surface_bc);
  else if (is_velocity (v, GFS_DOMAIN (_sim)))
    dirichlet_bc (_cell);
  else /* Neumann */
    return 0.;
  gdouble flux;
  if ((_cell->flags & GFS_FLAG_DIRICHLET) == 0)
    flux = GFS_STATE (_cell)->solid->fv;
  else {
    GfsSolidVector * s = GFS_STATE (_cell)->solid;
    FttVector m = {1.,1.,1.};
    gfs_domain_solid_metric (GFS_DOMAIN (_sim), _cell, &m);
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&s->v.x)[c] = (&m.x)[c]*(s->s[2*c + 1] - s->s[2*c]);
    flux = gfs_cell_dirichlet_gradient_flux (_cell, v->i, -1, GFS_STATE (_cell)->solid->fv);
  }
  return gfs_dimensional_value (v, flux*pow (_sim->physical_params.L*ftt_cell_size (_cell), 
					     FTT_DIMENSION - 2.));
}

static gboolean overlaps (double x1, double y1, double x2, double y2)
{
  double h = ftt_cell_size (_cell)/2.;
  FttVector p, min = { G_MAXDOUBLE, G_MAXDOUBLE }, max = { -G_MAXDOUBLE, -G_MAXDOUBLE };
  ftt_cell_pos (_cell, &p);
  FttVector q[4];
  q[0] = p; q[0].x += h; q[0].y += h; 
  q[1] = p; q[1].x -= h; q[1].y += h; 
  q[2] = p; q[2].x -= h; q[2].y -= h; 
  q[3] = p; q[3].x += h; q[3].y -= h; 
  gfs_simulation_map_inverse_cell (_sim, q);
  int i;
  for (i = 0; i < 4; i++) {
    if (q[i].x < min.x) min.x = q[i].x;
    if (q[i].y < min.y) min.y = q[i].y;
    if (q[i].x > max.x) max.x = q[i].x;
    if (q[i].y > max.y) max.y = q[i].y;
  }
  return (min.x <= x2 && min.y <= y2 && max.x >= x1 && max.y >= y1);
}

static double mapv (double u, double v, FttComponent c) {
  FttVector p, q = {u, v, 0.};
  g_return_val_if_fail (_cell != NULL, 0.);
  ftt_cell_pos (_cell, &p);
  gfs_simulation_map_inverse_vector (_sim, &p, &q);
  return (&q.x)[c];
}

static double mapvx (double u, double v) { return mapv (u, v, FTT_X); }
static double mapvy (double u, double v) { return mapv (u, v, FTT_Y); }

#endif /* __FUNCTION_H__ */
