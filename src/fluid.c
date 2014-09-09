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
/*! \file
 * \brief Low-level discrete operators.
 */

#include <math.h>

#include <stdlib.h>

#include "fluid.h"
#include "variable.h"
#include "domain.h"
#include "solid.h"

/**
 * gfs_cell_face:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * This function is different from ftt_cell_face() because it takes
 * into account the solid fractions.
 *
 * Returns: the face of @cell in direction @d.
 */
FttCellFace gfs_cell_face (FttCell * cell,
			   FttDirection d)
{
  FttCellFace f = {cell, NULL, d};

  g_return_val_if_fail (cell != NULL, f);

  if (!GFS_IS_MIXED (cell) || GFS_STATE (cell)->solid->s[d] > 0.)
    f.neighbor = ftt_cell_neighbor (cell, d);
  return f;
}

typedef struct _Gradient Gradient;

/* grad(p) = -a*p(cell) + b*p(neighbor) + c */
struct _Gradient {
  gdouble a, b, c;
};

/* get_average_neighbor_value_stencil() needs to be updated whenever 
 * this function is modified
 */
static gdouble average_neighbor_value (const FttCellFace * face,
				       guint v,
				       gdouble * x)
{
  /* check for corner refinement violation (topology.fig) */
  g_assert (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell));
  
  if (FTT_CELL_IS_LEAF (face->neighbor))
    return GFS_VALUEI (face->neighbor, v);
  else {
    FttCellChildren children;
    gdouble av = 0., a = 0.;
    FttDirection od = FTT_OPPOSITE_DIRECTION (face->d);
    guint i, n;
    
    n = ftt_cell_children_direction (face->neighbor, od, &children);
    for (i = 0; i < n; i++)
      if (children.c[i] && GFS_VALUEI (children.c[i], v) != GFS_NODATA) {
	gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
	a += w;
	av += w*GFS_VALUEI (children.c[i], v);
      }
    if (a > 0.) {
      *x = 3./4.;
      return av/a;
    }
    else
      return GFS_VALUEI (face->cell, v);
  }
}

static void average_neighbor_value_stencil (const FttCellFace * face, guint v)
{
  /* check for corner refinement violation (topology.fig) */
  g_assert (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell));
  
  if (FTT_CELL_IS_LEAF (face->neighbor))
    GFS_VALUEI (face->neighbor, v) = 1.;
  else {
    FttCellChildren children;
    gdouble a = 0.;
    guint i, n;
    
    n = ftt_cell_children_direction (face->neighbor,
				     FTT_OPPOSITE_DIRECTION (face->d),
				     &children);
    for (i = 0; i < n; i++)
      if (children.c[i]) {
	a += 1.;
	GFS_VALUEI (children.c[i], v) = 1.;
      }
    if (a == 0.)
      GFS_VALUEI (face->cell, v) = 1.;
  }
}

#if FTT_2D

/* v = a*v(cell) + b 
 * 
 * Second order 1D interpolation.
 */
static GfsGradient interpolate_1D2 (FttCell * cell,
				    FttDirection d,
				    gdouble x,
				    guint v)
{
  GfsGradient p;
  FttCellFace f1, f2;
  gdouble p1 = 0., p2 = 0.;
  gdouble x1 = 1., x2 = 1.;
  gdouble a1, a2;

  g_return_val_if_fail (cell != NULL, p);
  g_return_val_if_fail (!GFS_IS_MIXED (cell), p);
  
  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor)
    p1 = average_neighbor_value (&f1, v, &x1);
  f2 = gfs_cell_face (cell, d);
  if (f2.neighbor)
    p2 = average_neighbor_value (&f2, v, &x2);

  a1 = x*(x - x2)/(x1*(x1 + x2));
  a2 = x*(x + x1)/(x2*(x1 + x2));

  p.a = 1. - a1 - a2;
  p.b = 0.;
  if (f1.neighbor)
    p.b += a1*p1;
  else
    p.a += a1;
  if (f2.neighbor)
    p.b += a2*p2;
  else
    p.a += a2;

  return p;
}

/* v = a*v(cell) + b 
 * 
 * First order 1D interpolation.
 *
 * get_interpolate_1D1_stencil() needs to be updated whenever 
 * this function is modified
 */
static GfsGradient interpolate_1D1 (FttCell * cell,
				    FttDirection d,
				    gdouble x,
				    guint v)
{
  GfsGradient p = { 1., 0. };
  FttCellFace f;

  f = gfs_cell_face (cell, d);
  if (f.neighbor) {
    gdouble x2 = 1.;
    gdouble p2 = average_neighbor_value (&f, v, &x2);
    if (p2 != GFS_NODATA) {
      gdouble a2 = x/x2;
      p.b += a2*p2;
      p.a -= a2;
    }
  }

  return p;
}

static void interpolate_1D1_stencil (FttCell * cell,
				     FttDirection d,
				     guint v)
{
  FttCellFace f;

  GFS_VALUEI (cell, v) = 1.;
  f = ftt_cell_face (cell, d);
  if (f.neighbor)
    average_neighbor_value_stencil (&f, v);
}

#else /* not FTT_2D */

/* v = a*v(cell) + b 
 * 
 * First order 2D interpolation.
 *
 * get_interpolate_2D1_stencil() needs to be updated whenever 
 * this function is modified                                 
 */
static GfsGradient interpolate_2D1 (FttCell * cell,
				    FttDirection d1, FttDirection d2,
				    gdouble x, gdouble y,
				    guint v)
{
  GfsGradient p = { 1., 0. };
  FttCellFace f1, f2;

  f1 = gfs_cell_face (cell, d1);
  if (f1.neighbor) {
    gdouble y1 = 1.;
    gdouble p1 = average_neighbor_value (&f1, v, &y1);
    if (p1 != GFS_NODATA) {
      gdouble a1 = y/y1;
      p.b += a1*p1;
      p.a -= a1;
    }
  }

  f2 = gfs_cell_face (cell, d2);
  if (f2.neighbor) {
    gdouble x2 = 1.;
    gdouble p2 = average_neighbor_value (&f2, v, &x2);
    if (p2 != GFS_NODATA) {
      gdouble a2 = x/x2;
      p.b += a2*p2;
      p.a -= a2;
    }
  }
  
  return p;
}

static void interpolate_2D1_stencil (FttCell * cell,
				     FttDirection d1, FttDirection d2,
				     guint v)
{
  FttCellFace f1, f2;

  GFS_VALUEI (cell, v) = 1.;
  f1 = ftt_cell_face (cell, d1);
  if (f1.neighbor)
    average_neighbor_value_stencil (&f1, v);
  f2 = ftt_cell_face (cell, d2);
  if (f2.neighbor)
    average_neighbor_value_stencil (&f2, v);
}

#endif /* not FTT_2D */

#if FTT_2D
static gint perpendicular[FTT_NEIGHBORS][FTT_CELLS] = 
  {{-1,  2, -1,  3},
   { 2, -1,  3, -1},
   { 1,  0, -1, -1},
   {-1, -1,  1,  0}};
#else  /* FTT_3D */
static gint perpendicular[FTT_NEIGHBORS][FTT_CELLS][2] = 
  {{{-1,-1},{2,4},{-1,-1},{3,4},{-1,-1},{2,5},{-1,-1},{3,5}},
   {{2,4},{-1,-1},{3,4},{-1,-1},{2,5},{-1,-1},{3,5},{-1,-1}},
   {{4,1},{4,0},{-1,-1},{-1,-1},{5,1},{5,0},{-1,-1},{-1,-1}},
   {{-1,-1},{-1,-1},{4,1},{4,0},{-1,-1},{-1,-1},{5,1},{5,0}},
   {{1,2},{0,2},{1,3},{0,3},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},
   {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{1,2},{0,2},{1,3},{0,3}}};
#endif /* FTT_3D */

/* gradient_fine_coarse_stencil() needs to be updated whenever 
   this function is modified  */
static Gradient gradient_fine_coarse (const FttCellFace * face, guint v)
{
  Gradient g;
  GfsGradient p;
#if FTT_2D
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  g_assert (face != NULL);
  g_assert (ftt_face_type (face) == FTT_FINE_COARSE);

  dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if FTT_2D
  g_assert (dp >= 0);
  p = interpolate_1D1 (face->neighbor, dp, 1./4., v);
#else  /* FTT_3D */
  g_assert (dp[0] >= 0 && dp[1] >= 0);
  p = interpolate_2D1 (face->neighbor, dp[0], dp[1], 1./4., 1./4., v);
#endif /* FTT_3D */

  g.a = 2./3.;
  g.b = 2.*p.a/3.;
  g.c = 2.*p.b/3.;

  return g;
}

#define REFINE_CORNER(cell) {if (cell && FTT_CELL_IS_LEAF (cell) && \
                              ftt_cell_level (cell) < level - 1) \
	                    ftt_cell_refine_single (cell, init, data);}

void ftt_cell_refine_corners (FttCell * cell,
			      FttCellInitFunc init,
			      gpointer data)
{
  FttDirection d;
  FttCellNeighbors neighbor;
  guint level;

  g_return_if_fail (cell != NULL);

  level = ftt_cell_level (cell);
  ftt_cell_neighbors (cell, &neighbor);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (neighbor.c[d] && ftt_cell_level (neighbor.c[d]) < level) {
      if (GFS_CELL_IS_BOUNDARY (neighbor.c[d]))
	ftt_cell_refine_single (neighbor.c[d], init, data);
      else {
	FttCell * n;
#if FTT_2D
	gint dp;
#else  /* FTT_3D */
	gint * dp;
#endif /* FTT_3D */
	
	dp = perpendicular[d][FTT_CELL_ID (cell)];
#if FTT_2D
	g_assert (dp >= 0);
	n = ftt_cell_neighbor (neighbor.c[d], dp);
	REFINE_CORNER (n)
#else  /* FTT_3D */
	g_assert (dp[0] >= 0 && dp[1] >= 0);
	n = ftt_cell_neighbor (neighbor.c[d], dp[0]);
	REFINE_CORNER (n)
	n = ftt_cell_neighbor (neighbor.c[d], dp[1]);
	REFINE_CORNER (n)
#endif /* FTT_3D */
      }
    }
}

/**
 * gfs_neighbor_value:
 * @face: a #FttCellFace.
 * @v: a variable index.
 * @x: relative position of the neigboring value.
 *
 * Returns: the (averaged) neighboring value of variable @v at @face
 * and its position @x.
 */
gdouble gfs_neighbor_value (const FttCellFace * face,
			    guint v,
			    gdouble * x)
{
  g_return_val_if_fail (face != NULL, 0.);
  g_return_val_if_fail (x != NULL, 0.);

  GfsGradient vc;
#if FTT_2D
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  if (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell))
    /* neighbor at same level */
    return average_neighbor_value (face, v, x);
  else {
    if (GFS_VALUEI (face->neighbor, v) == GFS_NODATA)
      return GFS_NODATA;
    /* neighbor at coarser level */
    dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if FTT_2D
    g_assert (dp >= 0);
    vc = interpolate_1D1 (face->neighbor, dp, 1./4., v);
#else  /* FTT_3D */
    g_assert (dp[0] >= 0 && dp[1] >= 0);
    vc = interpolate_2D1 (face->neighbor, dp[0], dp[1], 1./4., 1./4., v);
#endif /* FTT_3D */
    *x = 3./2.;
    return vc.a*GFS_VALUEI (face->neighbor, v) + vc.b;
  }
}

static void neighbor_value_stencil (const FttCellFace * face, guint v)
{
#if FTT_2D
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  if (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell))
    /* neighbor at same level */
    average_neighbor_value_stencil (face, v);
  else {
    /* neighbor at coarser level */
    dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if FTT_2D
    g_assert (dp >= 0);
    interpolate_1D1_stencil (face->neighbor, dp, v);
#else  /* FTT_3D */
    g_assert (dp[0] >= 0 && dp[1] >= 0);
    interpolate_2D1_stencil (face->neighbor, dp[0], dp[1], v);
#endif /* FTT_3D */
    GFS_VALUEI (face->neighbor, v) = 1.;
  }
}

/**
 * gfs_center_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_gradient (FttCell * cell,
			     FttComponent c,
			     guint v)
{
  FttDirection d = 2*c;
  FttCellFace f1;
  gdouble v0;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  v0 = GFS_VALUEI (cell, v);
  if (f1.neighbor) {
    FttCellFace f2 = gfs_cell_face (cell, d);
    gdouble x1 = 1., v1;
    
    v1 = gfs_neighbor_value (&f1, v, &x1);
    if (f2.neighbor) {
      /* two neighbors: second-order differencing (parabola) */
      gdouble x2 = 1., v2;
      
      v2 = gfs_neighbor_value (&f2, v, &x2);
      return (x1*x1*(v2 - v0) + x2*x2*(v0 - v1))/(x1*x2*(x2 + x1));
    }
    else
      /* one neighbor: first-order differencing */
      return (v0 - v1)/x1;
  }
  else {
    FttCellFace f2 = gfs_cell_face (cell, d);

    if (f2.neighbor) {
      gdouble x2 = 1.;
      
      /* one neighbor: first-order differencing */
      return (gfs_neighbor_value (&f2, v, &x2) - v0)/x2;
    }
  }
  /* no neighbors */
  return 0.;
}

/**
 * gfs_center_gradient_stencil:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * Sets to 1. the @v variable of all the cells which would be used if
 * gfs_center_gradient() was called with identical arguments.
 */
void gfs_center_gradient_stencil (FttCell * cell,
				  FttComponent c,
				  guint v)
{
  FttDirection d = 2*c;
  FttCellFace f1, f2;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (c < FTT_DIMENSION);

  f1 = ftt_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor == cell) /* periodic */
    return;
  if (f1.neighbor) {
    GFS_VALUEI (cell, v) = 1.;
    neighbor_value_stencil (&f1, v);
  }
  f2 = ftt_cell_face (cell, d);
  if (f2.neighbor) {
    GFS_VALUEI (cell, v) = 1.;
    neighbor_value_stencil (&f2, v);
  }
}

/**
 * gfs_center_van_leer_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell and is limited
 * using van Leer's limiter.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_van_leer_gradient (FttCell * cell,
				      FttComponent c,
				      guint v)
{
  FttDirection d = 2*c;
  FttCellFace f1;
  
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor) {
    FttCellFace f2 = gfs_cell_face (cell, d);
    
    if (f2.neighbor) {
      /* two neighbors: second-order differencing (parabola)
	 + van Leer limiter. See http://en.wikipedia.org/wiki/Flux_limiter */
      gdouble x1 = 1., x2 = 1., v0, v1, v2;      
      v0 = GFS_VALUEI (cell, v);
      v1 = gfs_neighbor_value (&f1, v, &x1);
      v2 = gfs_neighbor_value (&f2, v, &x2);
      if (v1 == GFS_NODATA || v2 == GFS_NODATA)
	return 0.;

      gdouble s1 = 2.*(v0 - v1);
      gdouble s2 = 2.*(v2 - v0);

      if (s1*s2 <= 0.)
	return 0.;
      gdouble s0 = (x1*x1*(v2 - v0) + x2*x2*(v0 - v1))/(x1*x2*(x2 + x1));
      if (ABS (s2) < ABS (s1))
	s1 = s2;
      if (ABS (s0) < ABS (s1))
	return s0;
      return s1;
    }
  }
  /* only one or no neighbors */
  return 0.;
}

static gdouble generic_limiter (gdouble r, gdouble beta)
{
  gdouble v1 = MIN (r, beta), v2 = MIN (beta*r, 1.);
  v1 = MAX (0., v1);
  return MAX (v1, v2);
}

static gdouble minmod_limiter (gdouble r)
{
  return generic_limiter (r, 1.);
}

static gdouble superbee_limiter (gdouble r)
{
  return generic_limiter (r, 2.);
}

static gdouble sweby_limiter (gdouble r)
{
  return generic_limiter (r, 1.5);
}

static gdouble center_limited_gradient (FttCell * cell,
					FttComponent c,
					guint v,
					gdouble (* limiter) (gdouble))
{
  FttDirection d = 2*c;
  FttCellFace f1;
  gdouble v0;

  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  v0 = GFS_VALUEI (cell, v);
  if (f1.neighbor) {
    FttCellFace f2 = gfs_cell_face (cell, d);
    if (f2.neighbor) {
      /* two neighbors */
      gdouble x1 = 1., v1, x2 = 1., v2;
      v1 = gfs_neighbor_value (&f1, v, &x1);
      v2 = gfs_neighbor_value (&f2, v, &x2);

      gdouble g;
      if (v0 == v1)
	g = 0.;
      else
	g = (* limiter) ((v2 - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
      return g;
    }
  }
  /* only one or no neighbors */
  return 0.;
}

/**
 * gfs_center_minmod_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell and is limited
 * using a minmod limiter.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_minmod_gradient (FttCell * cell,
				    FttComponent c,
				    guint v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  return center_limited_gradient (cell, c, v, minmod_limiter);
}

/**
 * gfs_center_superbee_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell and is limited
 * using a superbee limiter.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_superbee_gradient (FttCell * cell,
				      FttComponent c,
				      guint v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  return center_limited_gradient (cell, c, v, superbee_limiter);
}

/**
 * gfs_center_sweby_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell and is limited
 * using a Sweby limiter (beta = 1.5).
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_sweby_gradient (FttCell * cell,
				   FttComponent c,
				   guint v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  return center_limited_gradient (cell, c, v, sweby_limiter);
}

/**
 * gfs_center_regular_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable.
 *
 * The gradient is normalized by the size of the cell. Only regular
 * Cartesian stencils are used to compute the gradient.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell. 
 */
gdouble gfs_center_regular_gradient (FttCell * cell,
				     FttComponent c,
				     GfsVariable * v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  FttCell * n1 = ftt_cell_neighbor (cell, 2*c);
  guint level = ftt_cell_level (cell);
  if (n1) {
    if (ftt_cell_level (n1) < level)
      return gfs_center_regular_gradient (ftt_cell_parent (cell), c, v)/2.;
    FttCell * n2 = ftt_cell_neighbor (cell, 2*c + 1);
    if (n2) {
      if (ftt_cell_level (n2) < level)
	return gfs_center_regular_gradient (ftt_cell_parent (cell), c, v)/2.;
      /* two neighbors: second-order differencing (parabola) */
      return (GFS_VALUE (n1, v) - GFS_VALUE (n2, v))/2.;
    }
    else
      /* one neighbor: first-order differencing */
      return GFS_VALUE (n1, v) - GFS_VALUE (cell, v);
  }
  else {
    FttCell * n2 = ftt_cell_neighbor (cell, 2*c + 1);
    if (n2) {
      if (ftt_cell_level (n2) < level)
	return gfs_center_regular_gradient (ftt_cell_parent (cell), c, v)/2.;
      /* one neighbor: first-order differencing */
      return GFS_VALUE (cell, v) - GFS_VALUE (n2, v);
    }
  }
  /* no neighbors */
  return 0.;
}

/**
 * gfs_center_regular_2nd_derivative:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable.
 *
 * The derivative is normalized by the size of the cell squared. Only
 * regular Cartesian stencils are used to compute the derivative.
 *
 * Returns: the value of the @c component of the 2nd derivative of
 * variable @v at the center of the cell.
 */
gdouble gfs_center_regular_2nd_derivative (FttCell * cell, 
					   FttComponent c, 
					   GfsVariable * v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  FttCell * n1 = ftt_cell_neighbor (cell, 2*c);
  FttCell * n2 = ftt_cell_neighbor (cell, 2*c + 1);
  if (n1 && n2) {
    guint level = ftt_cell_level (cell);
    if (ftt_cell_level (n1) < level || ftt_cell_level (n2) < level)
      return gfs_center_regular_2nd_derivative (ftt_cell_parent (cell), c, v)/4.;
    return GFS_VALUE (n1, v) - 2.*GFS_VALUE (cell, v) + GFS_VALUE (n2, v);
  }
  /* one or no neighbors */
  return 0.;
}

/**
 * gfs_face_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the
 * @face. The value returned is second order accurate in space and
 * conservative, in the sense that values at a coarse/fine cell
 * boundary are consistent.
 *
 * The value of the gradient (normalised by the size of @face->cell)
 * is given by: @g->b - @g->a*GFS_VALUE (@face->cell, @v).
 */
void gfs_face_gradient (const FttCellFace * face,
			GfsGradient * g,
			guint v,
			gint max_level)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL || GFS_FACE_FRACTION (face) == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    Gradient gcf;

    gcf = gradient_fine_coarse (face, v);
    g->a = gcf.a;
    g->b = gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c;
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      g->a = 1.;
      g->b = GFS_VALUEI (face->neighbor, v);
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      gdouble s;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++)
	if ((f.cell = children.c[i])) {
	  Gradient gcf;
	  gcf = gradient_fine_coarse (&f, v);
	  s = GFS_FACE_FRACTION (&f);
	  g->a += s*gcf.b;
	  g->b += s*(gcf.a*GFS_VALUEI (f.cell, v) - gcf.c);
	}
      s = GFS_FACE_FRACTION (face)*n/2.;
      g->a /= s;
      g->b /= s;
    }
  }
}

/* face_weighted_gradient_stencil() needs to be updated whenever
   this function is modified                                   */
static void face_weighted_gradient (const FttCellFace * face,
				    GfsGradient * g,
				    guint v,
				    gint max_level,
				    guint dimension)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    Gradient gcf;
    gdouble w = GFS_STATE (face->cell)->f[face->d].v;

    gcf = gradient_fine_coarse (face, v);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      gdouble w = GFS_STATE (face->cell)->f[face->d].v;

      g->a = w;
      g->b = w*GFS_VALUEI (face->neighbor, v);
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  Gradient gcf;
	  gdouble w = GFS_STATE (f.cell)->f[f.d].v;
	
	  gcf = gradient_fine_coarse (&f, v);
	  g->a += w*gcf.b;
	  g->b += w*(gcf.a*GFS_VALUEI (f.cell, v) - gcf.c);
	}
      if (dimension > 2) {
	g->a /= n/2.;
	g->b /= n/2.;
      }
    }
  }
}

/**
 * gfs_face_weighted_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the @face
 * weighted by the value of the @v field of the face state vector of the
 * corresponding cell. The value returned is second order accurate in
 * space and conservative, in the sense that values at a coarse/fine
 * cell boundary are consistent.  
 */

void gfs_face_weighted_gradient (const FttCellFace * face,
				 GfsGradient * g,
				 guint v,
				 gint max_level)
{
  face_weighted_gradient (face, g, v, max_level, FTT_DIMENSION);
}

void gfs_face_weighted_gradient_2D (const FttCellFace * face,
				    GfsGradient * g,
				    guint v,
				    gint max_level)
{
  face_weighted_gradient (face, g, v, max_level, 2);
}

static void fullest_directions (const FttCellFace * face,
				FttDirection d[FTT_DIMENSION])
{
  FttComponent c = face->d/2, i;
  FttCell * mixed = GFS_IS_MIXED (face->cell) ? face->cell : face->neighbor;
  GfsSolidVector * s = GFS_STATE (mixed)->solid;

  d[0] = face->d;
  for (i = 1; i < FTT_DIMENSION; i++) {
    guint cp = (c + i) % FTT_DIMENSION;
    d[i] = s->s[2*cp] > s->s[2*cp + 1] ? 2*cp : 2*cp + 1;
  }
}

static FttCell * cell_corner_neighbor1 (FttCell * cell,
					FttDirection * d,
					gint max_level)
{
  if (!cell)
    return NULL;
  FttCell * neighbor = ftt_cell_neighbor (cell, d[0]);
  if (!neighbor)
    return NULL;
  else {
    guint level = ftt_cell_level (cell);
    if (ftt_cell_level (neighbor) < level)
      /* neighbor is at a shallower level */
      return neighbor;
    else {
      if (level == max_level || FTT_CELL_IS_LEAF (neighbor))
	/* neighbor is at the same level */
	return neighbor;
      else {
	/* neighbor is at a deeper level */
	guint i;
	FttDirection d1[FTT_DIMENSION];
	d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
	for (i = 1; i < FTT_DIMENSION; i++)
	  d1[i] = d[i];
	return ftt_cell_child_corner (neighbor, d1);
      }
    }
  }
}

#if FTT_2D
# define N_CELLS 4
#else  /* 3D */
# define N_CELLS 8
#endif /* 3D */

static gboolean inverse (gdouble mi[N_CELLS - 1][N_CELLS - 1])
{
#if FTT_2D
  gdouble m[N_CELLS - 1][N_CELLS - 1], det;
  guint i, j;
  
  for (i = 0; i < N_CELLS - 1; i++)
    for (j = 0; j < N_CELLS - 1; j++)
      m[i][j] = mi[i][j];
  
  det = (m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) -
	 m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
	 m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]));
  if (det == 0.)
    return FALSE;
  
  mi[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det; 
  mi[0][1] = (m[2][1]*m[0][2] - m[0][1]*m[2][2])/det;
  mi[0][2] = (m[0][1]*m[1][2] - m[1][1]*m[0][2])/det; 
  mi[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det; 
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])/det; 
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])/det; 
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])/det; 
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])/det; 
  mi[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det; 
#else /* 3D */
  gint indxc[N_CELLS - 1], indxr[N_CELLS - 1], ipiv[N_CELLS - 1];
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp;
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

  for (j = 0; j < N_CELLS - 1; j++) 
    ipiv[j] = -1;

  for (i = 0; i < N_CELLS - 1; i++) {
    big = 0.0;
    for (j = 0; j < N_CELLS - 1; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < N_CELLS - 1; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (mi[j][k]) >= big) {
	      big = fabs (mi[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < N_CELLS - 1; l++) 
	SWAP (mi[irow][l], mi[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (mi[icol][icol] == 0.)
      return FALSE;
    pivinv = 1.0/mi[icol][icol];
    mi[icol][icol] = 1.0;
    for (l = 0; l < N_CELLS - 1; l++) mi[icol][l] *= pivinv;
    for (ll = 0; ll < N_CELLS - 1; ll++)
      if (ll != icol) {
	dum = mi[ll][icol];
	mi[ll][icol] = 0.0;
	for (l = 0; l < N_CELLS - 1; l++)
	  mi[ll][l] -= mi[icol][l]*dum;
      }
  }
  for (l = N_CELLS - 1 - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < N_CELLS - 1; k++)
	SWAP (mi[k][indxr[l]], mi[k][indxc[l]]);
  }
#endif /* 3D */
  return TRUE;
}

#if (!FTT_2D)
static void draw_cell (FttCell * cell, gdouble r, gdouble g, gdouble b,
		       const gchar * name)
{
  FttVector p;
  gdouble size = ftt_cell_size (cell)/2.;

  ftt_cell_pos (cell, &p);
  fprintf (stderr,
	   "(geometry \"%s\" = OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "4 3 2 1 0 %g %g %g\n"
	   "4 4 5 6 7 %g %g %g\n"
	   "4 2 3 7 6 %g %g %g\n"
	   "4 0 1 5 4 %g %g %g\n"
	   "4 0 4 7 3 %g %g %g\n"
	   "4 1 2 6 5 %g %g %g\n"
	   ")\n",
	   name,
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b);
  gfs_cell_cm (cell, &p);
  size /= 8.;
  fprintf (stderr,
	   "(geometry \"cm %s\" = OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "4 3 2 1 0 %g %g %g\n"
	   "4 4 5 6 7 %g %g %g\n"
	   "4 2 3 7 6 %g %g %g\n"
	   "4 0 1 5 4 %g %g %g\n"
	   "4 0 4 7 3 %g %g %g\n"
	   "4 1 2 6 5 %g %g %g\n"
	   ")\n",
	   name,
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b);
}

static void output_error_mesh (FttCell ** n)
{
  draw_cell (n[0], 0., 0., 0., "n0");
  draw_cell (n[1], 0.5, 0., 0., "n1");
  draw_cell (n[2], 0., 0.5, 0., "n2");
  draw_cell (n[3], 0., 0., 0.5, "n3");
  draw_cell (n[4], 1., 1., 1., "n4");
  draw_cell (n[5], 1., 0., 0., "n5");
  draw_cell (n[6], 0., 1., 0., "n6");
  draw_cell (n[7], 0., 0., 1., "n7");
  g_assert_not_reached ();
}
#endif /* 3D */

static gboolean face_bilinear (const FttCellFace * face,
			       FttCell ** n,
			       FttVector * o,
			       void (*cell_pos) (const FttCell *, FttVector *),
			       gint max_level,
			       gdouble m[N_CELLS - 1][N_CELLS - 1])
{
  FttDirection d[3], d1[3];
  guint i;
  gdouble size = ftt_cell_size (face->cell);

  fullest_directions (face, d);
  n[0] = face->cell; n[1] = face->neighbor;
  d1[0] = d[1]; d1[1] = d[0]; d1[2] = d[2];
  if ((n[2] = cell_corner_neighbor1 (n[0], d1, max_level)) == NULL)
    return FALSE;
  d1[1] = FTT_OPPOSITE_DIRECTION (d[0]);
  if ((n[3] = cell_corner_neighbor1 (n[1], d1, max_level)) == NULL)
    return FALSE;
  if (n[3] == n[2]) {
    d1[0] = d[0]; d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
    if ((n[3] = cell_corner_neighbor1 (n[2], d1, max_level)) == NULL)
      return FALSE;
  }

#if FTT_2D
  for (i = 0; i < 3; i++) {
    FttVector cm;
    guint j;

    (*cell_pos) (n[i + 1], &cm);

    for (j = 0; j < FTT_DIMENSION; j++) {
      (&cm.x)[j] -= (&o->x)[j];
      (&cm.x)[j] /= size;
      g_assert (fabs ((&cm.x)[j]) <= 4.);
    }

    m[i][0] = cm.x;
    m[i][1] = cm.y; 
    m[i][2] = cm.x*cm.y;
  }
  g_assert (inverse (m));
#else /* 3D */
  d1[0] = d[2]; d1[1] = d[0]; d1[2] = d[1];
  if ((n[4] = cell_corner_neighbor1 (n[0], d1, max_level)) == NULL)
    return FALSE;
  d1[1] = FTT_OPPOSITE_DIRECTION (d[0]);
  if ((n[5] = cell_corner_neighbor1 (n[1], d1, max_level)) == NULL)
    return FALSE;
  if (n[5] == n[4]) {
    d1[0] = d[0]; d1[1] = d[1]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[5] = cell_corner_neighbor1 (n[4], d1, max_level)) == NULL)
      return FALSE;
  }
  d1[0] = d[2]; d1[1] = d[0]; d1[2] = FTT_OPPOSITE_DIRECTION (d[1]);
  if ((n[6] = cell_corner_neighbor1 (n[2], d1, max_level)) == NULL)
    return FALSE;
  if (n[6] == n[4]) {
    d1[0] = d[1]; d1[1] = d[0]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[6] = cell_corner_neighbor1 (n[4], d1, max_level)) == NULL)
      return FALSE;
  }
  d1[0] = d[2]; d1[1] = FTT_OPPOSITE_DIRECTION (d[0]); d1[2] = FTT_OPPOSITE_DIRECTION (d[1]);
  if ((n[7] = cell_corner_neighbor1 (n[3], d1, max_level)) == NULL)
    return FALSE;
  if (n[7] == n[4] || n[7] == n[5]) {
    d1[0] = d[1]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[7] = cell_corner_neighbor1 (n[5], d1, max_level)) == NULL)
      return FALSE;
  }
  if (n[7] == n[6]) {
    d1[0] = d[0]; d1[1] = FTT_OPPOSITE_DIRECTION (d[1]); d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[7] = cell_corner_neighbor1 (n[6], d1, max_level)) == NULL)
      return FALSE;
  }

  for (i = 0; i < 7; i++) {
    FttVector cm;
    guint j;

    for (j = i + 1; j < 7; j++)
      if (n[i + 1] == n[j + 1])
	output_error_mesh (n);

    (*cell_pos) (n[i + 1], &cm);

    for (j = 0; j < FTT_DIMENSION; j++) {
      (&cm.x)[j] -= (&o->x)[j];
      (&cm.x)[j] /= size;
      if (fabs ((&cm.x)[j]) > 4.)
	output_error_mesh (n);
    }

    m[i][0] = cm.x;
    m[i][1] = cm.y; 
    m[i][2] = cm.z;
    m[i][3] = cm.x*cm.y;
    m[i][4] = cm.x*cm.z;
    m[i][5] = cm.y*cm.z;
    m[i][6] = cm.x*cm.y*cm.z;
  }
  if (!inverse (m))
    output_error_mesh (n);
#endif /* 3D */

  return TRUE;
}

/* mixed_face_gradient_stencil needs to be updated whenever 
 * this function is modified
 */
/* grad(v) = -a*v(cell) + b*v(neighbor) + c */
static gboolean mixed_face_gradient (const FttCellFace * face,
				     Gradient * g,
				     guint v,
				     gint max_level)
{
  FttCell * n[N_CELLS];
  gdouble m[N_CELLS - 1][N_CELLS - 1];
  FttVector o, cm;
  FttComponent c = face->d/2;
  gdouble h = ftt_cell_size (face->cell);

  gfs_cell_cm (face->cell, &o);
  if (!face_bilinear (face, n, &o, gfs_cell_cm, max_level, m))
    return FALSE;

  gfs_face_ca (face, &cm);

#if FTT_2D
  {
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
    gdouble vp;
    
    vp = ((&cm.x)[cp] - (&o.x)[cp])/h;
    g->a = ((m[c][0] + m[2][0]*vp) +
	    (m[c][1] + m[2][1]*vp) +
	    (m[c][2] + m[2][2]*vp));
    g->b = m[c][0] + m[2][0]*vp;
    g->c = ((m[c][1] + m[2][1]*vp)*GFS_VALUEI (n[2], v) +
	    (m[c][2] + m[2][2]*vp)*GFS_VALUEI (n[3], v));
  }
#else /* 3D */
  {
    guint j;

    cm.x = (cm.x - o.x)/h;
    cm.y = (cm.y - o.y)/h;
    cm.z = (cm.z - o.z)/h;
    g->c = 0.;
    
    switch (c) {
    case FTT_X:
      g->a = g->b = m[0][0] + cm.y*m[3][0] + cm.z*m[4][0] + cm.y*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[0][j] + cm.y*m[3][j] + cm.z*m[4][j] + cm.y*cm.z*m[6][j];
	g->a += a;
	g->c += a*GFS_VALUEI (n[j+1], v);
      }
      break;
    case FTT_Y:
      g->a = g->b = m[1][0] + cm.x*m[3][0] + cm.z*m[5][0] + cm.x*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[1][j] + cm.x*m[3][j] + cm.z*m[5][j] + cm.x*cm.z*m[6][j];
	g->a += a;
	g->c += a*GFS_VALUEI (n[j+1], v);
      }
      break;
    case FTT_Z:
      g->a = g->b = m[2][0] + cm.x*m[4][0] + cm.y*m[5][0] + cm.x*cm.y*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[2][j] + cm.x*m[4][j] + cm.y*m[5][j] + cm.x*cm.y*m[6][j];
	g->a += a;
	g->c += a*GFS_VALUEI (n[j+1], v);
      }
      break;
    default:
      g_assert_not_reached ();
    }
  }
#endif /* 3D */

  if (!FTT_FACE_DIRECT (face)) {
    g->a = - g->a;
    g->b = - g->b;
    g->c = - g->c;
  }

  return TRUE;
}

/* face_cm_gradient_stencil needs to be updated whenever 
 * this function is modified
 */
static void face_cm_gradient (const FttCellFace * face,
			      GfsGradient * g,
			      guint v,
			      gint max_level,
			      gboolean weighted)
{
  guint level;
  Gradient gcf;
  gdouble w = weighted ? GFS_STATE (face->cell)->f[face->d].v : 1.;

  g->a = g->b = 0.;
  
  if (face->neighbor == NULL || w == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    if (GFS_IS_MIXED (face->cell) || GFS_IS_MIXED (face->neighbor)) {
      if (!mixed_face_gradient (face, &gcf, v, max_level))
	gcf = gradient_fine_coarse (face, v);
    }
    else
      gcf = gradient_fine_coarse (face, v);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      if (!GFS_IS_MIXED (face->cell) && !GFS_IS_MIXED (face->neighbor)) {
	g->a = w;
	g->b = w*GFS_VALUEI (face->neighbor, v);
      }
      else if (mixed_face_gradient (face, &gcf, v, max_level)) {
	g->a = w*gcf.a;
	g->b = w*(gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c);
      }
      else {
	g->a = w;
	g->b = w*GFS_VALUEI (face->neighbor, v);
      }
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  w = weighted ? GFS_STATE (f.cell)->f[f.d].v : 1.;
	  if (GFS_IS_MIXED (f.cell) || GFS_IS_MIXED (f.neighbor)) {
	    if (!mixed_face_gradient (&f, &gcf, v, max_level))
	      gcf = gradient_fine_coarse (&f, v);
	  }
	  else
	    gcf = gradient_fine_coarse (&f, v);
	  g->a += w*gcf.b;
	  g->b += w*(gcf.a*GFS_VALUEI (f.cell, v) - gcf.c);
	}
#if !FTT_2D
      g->a /= n/2.;
      g->b /= n/2.;
#endif
    }
  }
}

/**
 * gfs_face_cm_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the
 * @face. Variable @v is defined at the center of mass of its
 * cell. Linear interpolation is used to evaluate the gradient in the
 * vicinity of cut cells.
 */
void gfs_face_cm_gradient (const FttCellFace * face,
			   GfsGradient * g,
			   guint v,
			   gint max_level)
{
  g_return_if_fail (face != NULL);
  g_return_if_fail (g != NULL);

  face_cm_gradient (face, g, v, max_level, FALSE);
}

/**
 * gfs_face_cm_weighted_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the @face
 * weighted by the value of the @v field of the face state vector of
 * the corresponding cell. Variable @v is defined at the center of
 * mass of its cell. Linear interpolation is used to evaluate the
 * gradient in the vicinity of cut cells.
 */
/* gfs_face_cm_weighted_gradient_stencil needs to be updated whenever 
 * this function is modified
 */
void gfs_face_cm_weighted_gradient (const FttCellFace * face,
				    GfsGradient * g,
				    guint v,
				    gint max_level)
{
  g_return_if_fail (face != NULL);
  g_return_if_fail (g != NULL);

  face_cm_gradient (face, g, v, max_level, TRUE);
}

static gboolean cell_bilinear (FttCell * cell,
			       FttCell ** n,
			       FttVector * o,
			       void (*cell_pos) (const FttCell *, FttVector *),
			       gint max_level,
			       gdouble m[N_CELLS - 1][N_CELLS - 1])
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;  
  FttCellFace f;
  FttDirection d[FTT_DIMENSION];
  FttComponent c;

  if ((s->s[FTT_RIGHT] == 0. && s->s[FTT_LEFT] == 0.) ||
      (s->s[FTT_RIGHT] == 1. && s->s[FTT_LEFT] == 1.))
    return FALSE;

  for (c = 0; c < FTT_DIMENSION; c++)
    d[c] = s->s[2*c] > s->s[2*c + 1] ? 2*c : 2*c + 1;
  f.cell = cell;
  f.d = d[0];
  f.neighbor = cell_corner_neighbor1 (cell, d, max_level);

  return face_bilinear (&f, n, o, cell_pos, max_level, m);
}

/**
 * gfs_cell_dirichlet_gradient:
 * @cell: a #FttCell.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @v0: the Dirichlet value on the boundary.
 * @grad: a #FttVector.
 *
 * Fills @grad with components of the gradient of variable @v
 * interpolated at the center of area of the solid boundary contained
 * in @cell. The gradient is scaled by the size of the cell.
 */
/* gfs_cell_dirichlet_gradient_stencil needs to be updated whenever 
   this function is modified  */
void gfs_cell_dirichlet_gradient (FttCell * cell,
				  guint v,
				  gint max_level,
				  gdouble v0,
				  FttVector * grad)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (grad != NULL);

  if (!GFS_IS_MIXED (cell))
    return;
  else {
    FttCell * n[N_CELLS];
    gdouble m[N_CELLS - 1][N_CELLS - 1];
    guint i, c;

    grad->x = grad->y = grad->z = 0.;
    if (!cell_bilinear (cell, n, &GFS_STATE (cell)->solid->ca, 
			gfs_cell_cm, max_level, m))
      return;

    for (i = 0; i < N_CELLS - 1; i++)
      for (c = 0; c < FTT_DIMENSION; c++)
	(&grad->x)[c] += m[c][i]*(GFS_VALUEI (n[i + 1], v) - v0);
  }
}

/**
 * gfs_mixed_cell_gradient:
 * @cell: a mixed #FttCell.
 * @v: a #GfsVariable.
 * @g: the gradient.
 *
 * Fills @g with the components of the gradient of @v at the center of
 * mass of @cell.
 *
 * The gradient is normalized by the size of the cell.
 */
void gfs_mixed_cell_gradient (FttCell * cell,
			      GfsVariable * v,
			      FttVector * g)
{
  FttCell * n[N_CELLS];
  gdouble m[N_CELLS - 1][N_CELLS - 1];
  gdouble v0, h;
  FttVector * o, cm;
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (GFS_IS_MIXED (cell));
  g_return_if_fail (v != NULL);
  g_return_if_fail (g != NULL);

  g->x = g->y = g->z = 0.;

  o = &GFS_STATE (cell)->solid->cm;
  v0 = GFS_VALUE (cell, v);
  cm = *o;

  if (v->surface_bc) {
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0) {
      o = &GFS_STATE (cell)->solid->ca;
      v0 = GFS_STATE (cell)->solid->fv;
    }
  }
  g_assert (cell_bilinear (cell, n, o, gfs_cell_cm, -1, m));
  
  h = ftt_cell_size (cell);
  cm.x = (cm.x - o->x)/h;
  cm.y = (cm.y - o->y)/h;
  cm.z = (cm.z - o->z)/h;
  for (i = 0; i < N_CELLS - 1; i++) {
    gdouble val = GFS_VALUE (n[i + 1], v) - v0;
#if FTT_2D
    g->x += (m[0][i] + m[2][i]*cm.y)*val;
    g->y += (m[1][i] + m[2][i]*cm.x)*val;
#else /* 3D */
    g->x += (m[0][i] + m[3][i]*cm.y + m[4][i]*cm.z + m[6][i]*cm.y*cm.z)*val;
    g->y += (m[1][i] + m[3][i]*cm.x + m[5][i]*cm.z + m[6][i]*cm.x*cm.z)*val;
    g->z += (m[2][i] + m[4][i]*cm.x + m[5][i]*cm.y + m[6][i]*cm.x*cm.y)*val;
#endif /* 3D */
  }
}

/**
 * gfs_mixed_cell_interpolate:
 * @cell: a mixed #FttCell.
 * @p: a #FttVector.
 * @v: a #GfsVariable.
 *
 * Returns: the value of variable @v interpolated at position @p within @cell.
 */
gdouble gfs_mixed_cell_interpolate (FttCell * cell,
				    FttVector p,
				    GfsVariable * v)
{
  FttCell * n[N_CELLS];
  gdouble m[N_CELLS - 1][N_CELLS - 1], a[N_CELLS - 1];
  gdouble v0, h;
  FttVector * o;
  guint i, j;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (GFS_IS_MIXED (cell), 0.);
  g_return_val_if_fail (v != NULL, 0.);

  o = &GFS_STATE (cell)->solid->cm;
  v0 = GFS_VALUE (cell, v);

  if (v->surface_bc) {
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0) {
      o = &GFS_STATE (cell)->solid->ca;
      v0 = GFS_STATE (cell)->solid->fv;
    }
  }
  g_assert (cell_bilinear (cell, n, o, gfs_cell_cm, -1, m));

  for (i = 0; i < N_CELLS - 1; i++) {
    a[i] = 0.;
    for (j = 0; j < N_CELLS - 1; j++)
      a[i] += m[i][j]*(GFS_VALUE (n[j + 1], v) - v0);
  }
  
  h = ftt_cell_size (cell);
  p.x = (p.x - o->x)/h;
  p.y = (p.y - o->y)/h;
#if FTT_2D
  return a[0]*p.x + a[1]*p.y + a[2]*p.x*p.y + v0;
#else /* 3D */
  p.z = (p.z - o->z)/h;
  return (a[0]*p.x + a[1]*p.y + a[2]*p.z + 
	  a[3]*p.x*p.y + a[4]*p.x*p.z + a[5]*p.y*p.z + 
	  a[6]*p.x*p.y*p.z + v0);
#endif /* 3D */
}

/**
 * gfs_cell_dirichlet_gradient_flux:
 * @cell: a #FttCell.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @v0: the Dirichlet value on the boundary.
 *
 * Returns: the flux of the gradient of variable @v through the solid
 * boundary contained in @cell.
 */
/* gfs_cell_dirichlet_gradient_flux_stencil needs to be updated whenever 
   this function is modified  */
gdouble gfs_cell_dirichlet_gradient_flux (FttCell * cell,
					  guint v,
					  gint max_level,
					  gdouble v0)
{
  g_return_val_if_fail (cell != NULL, 0.);

  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector g;    
    gfs_cell_dirichlet_gradient (cell, v, max_level, v0, &g);
    return g.x*s->v.x + g.y*s->v.y + g.z*s->v.z;
  }
}

static void gfs_cell_dirichlet_gradient_stencil (FttCell * cell,
						 gint max_level,
						 gdouble v0,
						 FttVector * grad,
						 GfsLinearProblem * lp,
						 GfsStencil * stencil,
						 FttVector * weight)
{
  if (GFS_IS_MIXED (cell)) {
    FttCell * n[N_CELLS];
    gdouble m[N_CELLS - 1][N_CELLS - 1];
    guint i, c;

    grad->x = grad->y = grad->z = 0.;
    if (!cell_bilinear (cell, n, &GFS_STATE (cell)->solid->ca, 
			gfs_cell_cm, max_level, m))
      return;

    for (i = 0; i < N_CELLS - 1; i++)
      for (c = 0; c < FTT_DIMENSION; c++) {
	gfs_stencil_add_element (stencil, n[i + 1], lp, m[c][i]*(&weight->x)[c]);
	(&grad->x)[c] += -m[c][i]*v0; /* Remaining Dirichlet contribution */
      }
  }
}

/**
 * gfs_cell_dirichlet_gradient_flux_stencil:
 * @cell: a #FttCell.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @v0: the Dirichlet value on the boundary.
 * @lp: the #GfsLinearProblem.
 * @stencil: a GfsStencil.
 *
 * Returns the stencil accounting for the flux of the gradient of
 * a given variable through the solid boundary contained in @cell.
 * Returns the stencil equivalent to the action of
 * gfs_cell_dirichlet_gradient_flux().
 *
 * Returns: the flux of the gradient of variable @v through the solid
 * boundary contained in @cell.
 */
gdouble gfs_cell_dirichlet_gradient_flux_stencil (FttCell * cell,
						  gint max_level,
						  gdouble v0,
						  GfsLinearProblem * lp,
						  GfsStencil * stencil)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (lp != NULL, 0.);
  g_return_val_if_fail (stencil != NULL, 0.);

  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector g;
    gfs_cell_dirichlet_gradient_stencil (cell, max_level, v0, &g, lp, stencil, &s->v);
    return g.x*s->v.x + g.y*s->v.y + g.z*s->v.z;
  }
}

/**
 * gfs_cell_dirichlet_value:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Returns: the value of variable @v interpolated at the center of
 * area of the solid boundary contained in @cell.
 */
gdouble gfs_cell_dirichlet_value (FttCell * cell,
				  GfsVariable * v,
				  gint max_level)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    FttCell * n[N_CELLS];
    FttVector p;
    gdouble m[N_CELLS - 1][N_CELLS - 1], a[N_CELLS - 1];
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    gdouble v0, size = ftt_cell_size (cell);
    void (* cell_pos) (const FttCell *, FttVector *) = v->centered ? ftt_cell_pos : gfs_cell_cm;
    guint i, j;

    (*cell_pos) (cell, &p);
    if (!cell_bilinear (cell, n, &p, cell_pos, max_level, m))
      return 0.;

    v0 = GFS_VALUE (cell, v);
    for (i = 0; i < N_CELLS - 1; i++) {
      a[i] = 0.;
      for (j = 0; j < N_CELLS - 1; j++)
	a[i] += m[i][j]*(GFS_VALUE (n[j + 1], v) - v0);
    }
    p.x = (s->ca.x - p.x)/size;
    p.y = (s->ca.y - p.y)/size;
#if FTT_2D
    return (a[0]*p.x + a[1]*p.y + a[2]*p.x*p.y + v0);
#else /* 3D */
    p.z = (s->ca.z - p.z)/size;
    return (a[0]*p.x + a[1]*p.y + a[2]*p.z + 
	    a[3]*p.x*p.y + a[4]*p.x*p.z + a[5]*p.y*p.z + 
	    a[6]*p.x*p.y*p.z + v0);
#endif /* 3D */
  }
}

/**
 * gfs_shear_strain_rate_tensor:
 * @cell: a #FttCell.
 * @u: the velocity.
 * @t: the shear strain rate tensor t[i][j] = (d_iu_j+d_ju_i)/2.
 *
 * Fills @t with the shear strain rate tensor at the center of mass of
 * @cell, normalised by the size of the cell.
 */
void gfs_shear_strain_rate_tensor (FttCell * cell, 
				   GfsVariable ** u,
				   gdouble t[FTT_DIMENSION][FTT_DIMENSION])
{
  guint i, j;
  FttVector g[FTT_DIMENSION];

  g_return_if_fail (cell != NULL);
  g_return_if_fail (u != NULL);

  for (i = 0; i < FTT_DIMENSION; i++)
    if (GFS_IS_MIXED (cell))
      gfs_mixed_cell_gradient (cell, u[i], &g[i]);
    else
      for (j = 0; j < FTT_DIMENSION; j++)
	(&g[i].x)[j] = gfs_center_gradient (cell, j, u[i]->i);			     

  for (i = 0; i < FTT_DIMENSION; i++) {
    t[i][i] = (&g[i].x)[i];
    for (j = i + 1; j < FTT_DIMENSION; j++)
      t[i][j] = ((&g[j].x)[i] + (&g[i].x)[j])/2.;
  }
  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < i; j++)
      t[i][j] = t[j][i];
}

/**
 * gfs_2nd_principal_invariant:
 * @cell: a #FttCell.
 * @u: the velocity.
 *
 * Returns: the second principal invariant of the shear strain rate
 * tensor of @cell: sqrt(D:D).
 */
gdouble gfs_2nd_principal_invariant (FttCell * cell, GfsVariable ** u)
{
  gdouble t[FTT_DIMENSION][FTT_DIMENSION];
  gdouble D = 0.;
  guint i, j;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (u != NULL, 0.);

  gfs_shear_strain_rate_tensor (cell, u, t);
  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < FTT_DIMENSION; j++)
      D += t[i][j]*t[i][j];
  return sqrt (D);
}

/**
 * gfs_get_from_below_intensive:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to "get from below".
 *
 * Sets the value of the "intensive" variable @v of @cell by taking
 * the volume weighted average of the values of its children cells.
 *
 * This functions fails if @cell is a leaf of the cell tree.
 */
void gfs_get_from_below_intensive (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0., sa = 0.;
  guint i;
  FttCellChildren child;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));
  g_return_if_fail (v != NULL);

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && GFS_HAS_DATA (child.c[i], v)) {
      gdouble a = gfs_domain_cell_fraction (v->domain, child.c[i]);
      val += GFS_VALUE (child.c[i], v)*a;
      sa += a;
    }
  if (sa > 0.)
    GFS_VALUE (cell, v) = val/sa;
  else
    GFS_VALUE (cell, v) = GFS_NODATA;
}

/**
 * gfs_cell_coarse_fine:
 * @parent: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Initializes @v on the children of @parent using interpolation.
 *
 * First-order interpolation (straight injection) is used for boundary
 * cells and second-order interpolation for the other cells.  
 */
void gfs_cell_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  g_return_if_fail (parent != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));
  g_return_if_fail (v != NULL);

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], v) = GFS_VALUE (parent, v);

  if (!GFS_CELL_IS_BOUNDARY (parent) && GFS_HAS_DATA (parent, v)) {
    FttVector g;
    FttComponent c;
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&g.x)[c] = gfs_center_van_leer_gradient (parent, c, v->i);

    if (v->domain->cell_metric) {
      gdouble a[FTT_CELLS], sa = 0.;
      for (n = 0; n < FTT_CELLS; n++) {
	a[n] = (* v->domain->cell_metric) (v->domain, child.c[n]);
	sa += a[n];
      }
      g_assert (sa > 0.);
      sa *= 2.;
#if FTT_2D
      double gx1 = g.x*(a[0] + a[2])/sa, gx2 = - g.x*(a[1] + a[3])/sa;
      double gy1 = g.y*(a[2] + a[3])/sa, gy2 = - g.y*(a[0] + a[1])/sa;
      GFS_VALUE (child.c[0], v) += gx2 + gy1;
      GFS_VALUE (child.c[1], v) += gx1 + gy1;
      GFS_VALUE (child.c[2], v) += gx2 + gy2;
      GFS_VALUE (child.c[3], v) += gx1 + gy2;
#else /* 3D */
      double gx1 = g.x*(a[0] + a[2])/sa, gx2 = - g.x*(a[1] + a[3])/sa;
      double gx3 = g.x*(a[4] + a[6])/sa, gx4 = - g.x*(a[5] + a[7])/sa;
      double gy1 = g.y*(a[2] + a[3])/sa, gy2 = - g.y*(a[0] + a[1])/sa;
      double gy3 = g.y*(a[6] + a[7])/sa, gy4 = - g.y*(a[4] + a[5])/sa;
      double gz1 = g.z*(a[2] + a[6])/sa, gz2 = - g.z*(a[0] + a[4])/sa;
      double gz3 = g.z*(a[3] + a[7])/sa, gz4 = - g.z*(a[1] + a[5])/sa;
      GFS_VALUE (child.c[0], v) += gx2 + gy1 + gz1;
      GFS_VALUE (child.c[1], v) += gx1 + gy1 + gz1;
      GFS_VALUE (child.c[2], v) += gx2 + gy2 + gz1;
      GFS_VALUE (child.c[3], v) += gx1 + gy2 + gz1;
      GFS_VALUE (child.c[4], v) += gx2 + gy1 + gz2;
      GFS_VALUE (child.c[5], v) += gx1 + gy1 + gz2;
      GFS_VALUE (child.c[6], v) += gx2 + gy2 + gz2;
      GFS_VALUE (child.c[7], v) += gx1 + gy2 + gz2;
#endif /* 3D */
    }
    else
      for (n = 0; n < FTT_CELLS; n++) 
	if (child.c[n]) {
	  FttVector p;
	  
	  ftt_cell_relative_pos (child.c[n], &p);
	  for (c = 0; c < FTT_DIMENSION; c++)
	    GFS_VALUE (child.c[n], v) += (&p.x)[c]*(&g.x)[c];
	}
  }
}

/**
 * gfs_cell_cleanup:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain.
 *
 * Frees the memory allocated for extra data associated with @cell.
 *
 * This function must be used as "cleanup function" when using
 * ftt_cell_destroy().
 */
void gfs_cell_cleanup (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (domain != NULL);

  if (cell->data) {
    GSList * i = domain->variables;
    while (i) {
      GfsVariable * v = i->data;
      if (v->cleanup)
	(* v->cleanup) (cell, v);
      i = i->next;
    }
    if (GFS_STATE (cell)->solid) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }    
  }
  g_free (cell->data);
  cell->data = NULL;
}

/**
 * gfs_cell_reset:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to reset.
 *
 * Sets the value of the variable @v of @cell to zero.
 */
void gfs_cell_reset (FttCell * cell, GfsVariable * v)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  GFS_VALUE (cell, v) = 0.;
}

/**
 * gfs_face_reset:
 * @face: a #FttCellFace.
 * @v: a #GfsVariable to reset.
 *
 * Sets the value of the variable @v of @face to zero (on both sides).
 */
void gfs_face_reset (FttCellFace * face, GfsVariable * v)
{
  g_return_if_fail (face != NULL);
  g_return_if_fail (v != NULL);

  GFS_VALUE (face->cell, v) = GFS_VALUE (face->neighbor, v) = 0.;
}

static void add_stats (const FttCell * cell, gpointer * data)
{
  GtsRange * s = data[0];
  gdouble v = GFS_VALUE (cell, GFS_VARIABLE (data[1]));

  if (v != GFS_NODATA)
    gts_range_add_value (s, v);
}

/**
 * gfs_stats_variable:
 * @root: the root #FttCell of the tree to obtain statistics from.
 * @v: the variable to consider for statistics.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the cell tree defined by @root using ftt_cell_traverse()
 * and gathers statistics about variable @v.
 *
 * Returns: a #GtsRange containing the statistics about @v.
 */
GtsRange gfs_stats_variable (FttCell * root,
			     GfsVariable * v,
			     FttTraverseFlags flags,
			     gint max_depth)
{
  GtsRange s;
  gpointer data[2];

  g_return_val_if_fail (root != NULL, s);
  g_return_val_if_fail (v != NULL, s);
  
  gts_range_init (&s);
  data[0] = &s;
  data[1] = v;
  ftt_cell_traverse (root, FTT_PRE_ORDER, flags, max_depth, 
		     (FttCellTraverseFunc) add_stats, data);
  gts_range_update (&s);

  return s;
}

static void add_norm (const FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];

  gfs_norm_add (n, GFS_VALUE (cell, v), gfs_cell_volume (cell, v->domain));
}

/**
 * gfs_norm_variable:
 * @root: the root #FttCell of the tree to obtain norm from.
 * @v: the variable to consider for norm statistics.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the cell tree defined by @root using ftt_cell_traverse()
 * and gathers norm statistics about variable @v.
 *
 * Returns: a #GfsNorm containing the norm statistics about @v.
 */
GfsNorm gfs_norm_variable (FttCell * root,
			   GfsVariable * v,
			   FttTraverseFlags flags,
			   gint max_depth)
{
  GfsNorm n;
  gpointer data[2];

  g_return_val_if_fail (root != NULL, n);
  g_return_val_if_fail (v != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = &n;
  data[1] = v;
  ftt_cell_traverse (root, FTT_PRE_ORDER, flags, max_depth, 
		     (FttCellTraverseFunc) add_norm, data);
  gfs_norm_update (&n);

  return n;
}

/**
 * gfs_norm_init:
 * @n: a #GfsNorm.
 *
 * Initializes a #GfsNorm.
 */
void gfs_norm_init (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  n->bias = n->first = n->second = 0.;
  n->infty = - G_MAXDOUBLE;
  n->w = 0.;
}

/**
 * gfs_norm_reset:
 * @n: a #GfsNorm.
 *
 * Sets all the fields of @n to 0.
 */
void gfs_norm_reset (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  n->bias = n->first = n->second = 0.;
  n->infty = 0.;
  n->w = 0.;
}

/**
 * gfs_norm_add:
 * @n: a #GfsNorm.
 * @val: a value to add to @n.
 * @weight: weight of @val.
 *
 * Adds @val to @n.
 */
void gfs_norm_add (GfsNorm * n, gdouble val, gdouble weight)
{
  g_return_if_fail (n != NULL);

  if (val != GFS_NODATA) {
    n->bias += weight*val;
    val = fabs (val);
    if (weight != 0. && val > n->infty)
      n->infty = val;
    n->first += weight*val;
    n->second += weight*val*val;
    n->w += weight;
  }
}

/**
 * gfs_norm_update:
 * @n: a #GfsNorm.
 * 
 * Updates the fields of @n.
 */
void gfs_norm_update (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  if (n->w > 0.0) {
    n->bias /= n->w;
    n->first /= n->w;
    n->second = sqrt (n->second/n->w);
  }
  else
    n->infty = 0.0;
}

/**
 * gfs_face_interpolated_value:
 * @face: a #FttFace.
 * @v: a #GfsVariable index.
 *
 * Computes the value of variable @v on the @face using second-order
 * interpolation from the cell-centered values.
 *
 * Note that this function cannot be called for coarse -> fine faces
 * (use gfs_face_interpolated_value_generic() instead).
 *
 * Returns: the value of variable @v on the face.  
 */
gdouble gfs_face_interpolated_value (const FttCellFace * face,
				     guint v)
{
  gdouble x1 = 1., v1;
#if 1
  g_return_val_if_fail (face != NULL, 0.);

  if (face->neighbor) {
    g_assert (FTT_CELL_IS_LEAF (face->neighbor) || 
	      ftt_cell_level (face->neighbor) < ftt_cell_level (face->cell));
    v1 = gfs_neighbor_value (face, v, &x1);
    return ((x1 - 0.5)*GFS_VALUEI (face->cell, v) + 0.5*v1)/x1;
  }
  else
    return GFS_VALUEI (face->cell, v);
#else
  gdouble v0;
  FttCellFace f2;

  g_return_val_if_fail (face != NULL, 0.);

  v0 = GFS_VALUEI (face->cell, v);
  v1 = gfs_neighbor_value (face, v, &x1);
  f2 = gfs_cell_face (face->cell, FTT_OPPOSITE_DIRECTION (face->d));
  if (f2.neighbor) {
    gdouble x2 = 1.;
    gdouble v2 = gfs_neighbor_value (&f2, v, &x2);

    return v0 + (x2*(v1 - v0)*(1. + 2.*x2) - x1*(v0 - v2)*(1. - 2.*x1))
      /(4.*x1*x2*(x1 + x2));
  }
  else
    return ((x1 - 0.5)*v0 + 0.5*v1)/x1;
#endif
}

/**
 * gfs_face_interpolated_value_generic:
 * @face: a #FttFace.
 * @v: a #GfsVariable.
 *
 * Same as gfs_face_interpolated_value() but also works for coarse ->
 * fine faces.
 *
 * Returns: the value of variable @v on the face.  
 */
gdouble gfs_face_interpolated_value_generic (const FttCellFace * face, const GfsVariable * v)
{  
  g_return_val_if_fail (face != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);
  
  if (!face->neighbor || FTT_CELL_IS_LEAF (face->neighbor) ||
      ftt_cell_level (face->neighbor) < ftt_cell_level (face->cell))
    return gfs_face_interpolated_value (face, v->i);
  else {
    /* finer neighbor */
    FttCellFace f = { NULL, face->cell, FTT_OPPOSITE_DIRECTION (face->d) };
    FttCellChildren child;
    int i, n = ftt_cell_children_direction (face->neighbor, f.d, &child);
    gdouble avg = 0.;
    for (i = 0; i < n; i++)
      if (child.c[i]) {
	f.cell = child.c[i];
	avg += gfs_face_interpolated_value (&f, v->i)*gfs_domain_face_fraction (v->domain, &f);
      }
    return avg == 0. ? 0. : avg/(gfs_domain_face_fraction (v->domain, face)*n);
  }
}

/**
 * gfs_face_weighted_interpolated_value:
 * @face: a #FttFace.
 * @v: a #GfsVariable index.
 *
 * Computes the value of variable @v on the @face weighted by the
 * value of the @v field of the face state vector using interpolation
 * from the cell-centered values. The value returned is second order
 * accurate in space and conservative, in the sense that values at a
 * coarse/fine cell boundary are consistent.
 *
 * Returns: the weighted value of variable @v on the face.  
 */
gdouble gfs_face_weighted_interpolated_value (const FttCellFace * face,
					      guint v)
{
  g_return_val_if_fail (face != NULL, 0.);

  if (face->neighbor) {
    if (FTT_CELL_IS_LEAF (face->neighbor)) {
      gdouble w = GFS_STATE (face->cell)->f[face->d].v, x1 = 1., v1;
      v1 = gfs_neighbor_value (face, v, &x1);
      return w*((x1 - 0.5)*GFS_VALUEI (face->cell, v) + 0.5*v1)/x1;
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      gdouble val = 0.;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++)
	if ((f.cell = children.c[i])) {
	  gdouble w = GFS_STATE (f.cell)->f[f.d].v, x1 = 1., v1;
	  v1 = gfs_neighbor_value (&f, v, &x1);
	  val += w*v1;
	}
      return val/n;
    }
  }
  else
    return GFS_STATE (face->cell)->f[face->d].v*GFS_VALUEI (face->cell, v);
}

/**
 * gfs_normal_divergence:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Adds to variable @v of @cell the integral of the divergence
 * of the (MAC) velocity field in this cell.  
 */
void gfs_normal_divergence (FttCell * cell,
			    GfsVariable * v)
{
  FttCellFace face;
  gdouble div = 0.;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  face.cell = cell;
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++)
    div += (FTT_FACE_DIRECT (&face) ? 1. : -1.)*
      GFS_STATE (cell)->f[face.d].un*gfs_domain_face_fraction (v->domain, &face);
  GFS_VALUE (cell, v) = div*ftt_cell_size (cell);
}

/**
 * gfs_normal_divergence_2D:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the integral of the 2D
 * divergence of the (MAC) velocity field in this cell.
 */
void gfs_normal_divergence_2D (FttCell * cell,
			       GfsVariable * v)
{
  FttCellFace face;
  gdouble div = 0.;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);
  
  face.cell = cell;
  for (face.d = 0; face.d < 4; face.d++)
    div += (FTT_FACE_DIRECT (&face) ? 1. : -1.)*
      GFS_STATE (cell)->f[face.d].un*gfs_domain_face_fraction (v->domain, &face);
  GFS_VALUE (cell, v) = div*ftt_cell_size (cell);
}

/**
 * gfs_divergence:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * Returns: the divergence of the (centered) vector field @v in @cell.
 */
gdouble gfs_divergence (FttCell * cell,
			GfsVariable ** v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  GfsDomain * domain = v[0]->domain;
  gdouble div = 0.;
  FttCellNeighbors n;
  ftt_cell_neighbors (cell, &n);
  FttCellFace f;
  f.cell = cell;
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++)
    if (n.c[f.d]) {
      f.neighbor = n.c[f.d];
      div += gfs_domain_face_fraction (domain, &f)*
	(FTT_FACE_DIRECT (&f) ? 1. : -1.)*gfs_face_interpolated_value_generic (&f, v[f.d/2]);
    }
  return div/(gfs_domain_cell_fraction (domain, cell)*ftt_cell_size (cell));
}

/**
 * gfs_vorticity:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * Returns: the vorticity (norm of the vorticity vector in 3D) of the
 * vector field @v in @cell.
 */
gdouble gfs_vorticity (FttCell * cell,
		       GfsVariable ** v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

#if FTT_2D
  if (GFS_IS_MIXED (cell) || GFS_IS_AXI (v[0]->domain))
    /* fixme: the formulation below should also work (better) for solid cells */
    return (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	    gfs_center_gradient (cell, FTT_Y, v[0]->i))/ftt_cell_size (cell);
  else {
    FttCellNeighbors n;
    ftt_cell_neighbors (cell, &n);
    FttCellFace f = { cell };
    double circ = 0.;
    f.d = FTT_RIGHT; f.neighbor = n.c[f.d];
    circ += gfs_face_interpolated_value_generic (&f, v[1])*
      gfs_domain_face_fraction (v[0]->domain, &f);
    f.d = FTT_LEFT; f.neighbor = n.c[f.d];
    circ -= gfs_face_interpolated_value_generic (&f, v[1])*
      gfs_domain_face_fraction (v[0]->domain, &f);
    f.d = FTT_TOP; f.neighbor = n.c[f.d];
    circ -= gfs_face_interpolated_value_generic (&f, v[0])*
      gfs_domain_face_fraction (v[0]->domain, &f);
    f.d = FTT_BOTTOM; f.neighbor = n.c[f.d];
    circ += gfs_face_interpolated_value_generic (&f, v[0])*
      gfs_domain_face_fraction (v[0]->domain, &f);
    /* vorticity is circulation over area */
    return circ*ftt_cell_size (cell)/gfs_cell_volume (cell, v[0]->domain);
  }
#else  /* FTT_3D */
  /* fixme: does not work with metric */
  FttVector vort;
  vort.x = (gfs_center_gradient (cell, FTT_Y, v[2]->i) -
	    gfs_center_gradient (cell, FTT_Z, v[1]->i));
  vort.y = (gfs_center_gradient (cell, FTT_Z, v[0]->i) -
	    gfs_center_gradient (cell, FTT_X, v[2]->i));
  vort.z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	    gfs_center_gradient (cell, FTT_Y, v[0]->i));
  return sqrt (vort.x*vort.x + vort.y*vort.y + vort.z*vort.z)/ftt_cell_size (cell);
#endif /* FTT_3D */
}

/**
 * gfs_vector_norm2:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * Returns: the squared norm of the vector field @v in @cell.
 */
gdouble gfs_vector_norm2 (FttCell * cell,
			  GfsVariable ** v)
{
  FttComponent c;
  gdouble n = 0.;
  
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  for (c = 0; c < FTT_DIMENSION; c++)
    n += GFS_VALUE (cell, v[c])*GFS_VALUE (cell, v[c]);
  return n;
}

/**
 * gfs_vector_norm:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * Returns: the norm of the vector field @v in @cell.
 */
gdouble gfs_vector_norm (FttCell * cell,
			 GfsVariable ** v)
{
  return sqrt (gfs_vector_norm2 (cell, v));
}

/**
 * gfs_vector_lambda2:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * Returns: The value of the lambda2 eigenvalue used by Jeong and
 * Hussain as vortex criterion (JFM 285, 69-94, 1995), normalized by
 * the square of the size of the cell.
 */
gdouble gfs_vector_lambda2 (FttCell * cell,
			    GfsVariable ** v)
{
  gdouble J[FTT_DIMENSION][FTT_DIMENSION];
  gdouble S2O2[FTT_DIMENSION][FTT_DIMENSION];
  gdouble lambda[FTT_DIMENSION], ev[FTT_DIMENSION][FTT_DIMENSION];
  guint i, j, k;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < FTT_DIMENSION; j++)
      J[i][j] = gfs_center_gradient (cell, j, v[i]->i);
  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < FTT_DIMENSION; j++) {
      S2O2[i][j] = 0.;
      for (k = 0; k < FTT_DIMENSION; k++)
	S2O2[i][j] += J[i][k]*J[k][j] + J[k][i]*J[j][k];
    }
  gfs_eigenvalues (S2O2, lambda, ev);
  return lambda[1]/2.;
}

/**
 * gfs_pressure_force:
 * @cell: a #FttCell.
 * @p: a #GfsVariable.
 * @f: a #FttVector.
 *
 * Fills @f with the pressure @p component of the force exerted by the
 * fluid on the fraction of embedded boundary contained in @cell.
 */
void gfs_pressure_force (FttCell * cell,
			 GfsVariable * p,
			 FttVector * f)
{
  GfsSolidVector * s;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (f != NULL);

  if ((s = GFS_STATE (cell)->solid)) {
    gdouble size = ftt_cell_size (cell);
    gdouble pv = gfs_cell_dirichlet_value (cell, p, -1)*size;
    FttComponent c;

#if (!FTT_2D)
    pv *= size;
#endif /* 3D */

    gfs_solid_normal (cell, f);
    for (c = 0; c < FTT_DIMENSION; c++)
      (&f->x)[c] *= pv;
  }
  else
    f->x = f->y = f->z = 0.;
}

static void cell_traverse_mixed (FttCell * cell,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  if (!GFS_IS_MIXED (cell))
    return;
  if (order == FTT_PRE_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, data);
  if (!FTT_CELL_IS_LEAF (cell)) {
    struct _FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_mixed (c, order, flags, func, data);
    }
  }
  if (order == FTT_POST_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, data);
}

/**
 * gfs_cell_traverse_mixed:
 * @root: the root #FttCell of the tree to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each mixed cell.
 */
void gfs_cell_traverse_mixed (FttCell * root,
			      FttTraverseType order,
			      FttTraverseFlags flags,
			      FttCellTraverseFunc func,
			      gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_mixed (root, order, flags, func, data);
}

#if FTT_2D
static FttDirection corner[4][FTT_DIMENSION] = {
  { FTT_LEFT,  FTT_BOTTOM },
  { FTT_RIGHT, FTT_BOTTOM },
  { FTT_RIGHT, FTT_TOP },
  { FTT_LEFT,  FTT_TOP }
};
#else  /* 3D */
static FttDirection corner[8][FTT_DIMENSION] = {
  { FTT_LEFT,  FTT_BOTTOM, FTT_FRONT },
  { FTT_RIGHT, FTT_BOTTOM, FTT_FRONT },
  { FTT_RIGHT, FTT_TOP,    FTT_FRONT },
  { FTT_LEFT,  FTT_TOP,    FTT_FRONT },
  { FTT_LEFT,  FTT_BOTTOM, FTT_BACK },
  { FTT_RIGHT, FTT_BOTTOM, FTT_BACK },
  { FTT_RIGHT, FTT_TOP,    FTT_BACK },
  { FTT_LEFT,  FTT_TOP,    FTT_BACK }
};
#endif /* 3D */

/**
 * gfs_cell_corner_values:
 * @cell: a #FttCell containing location @p.
 * @v: a #GfsVariable.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @f: an array with the correct size (4*(FTT_DIMENSION - 1) + 1).
 *
 * Fills @f with the values of @v interpolated at the corners of @cell.
 */
void gfs_cell_corner_values (FttCell * cell, 
			     GfsVariable * v, 
			     gint max_level,
			     gdouble f[4*(FTT_DIMENSION - 1) + 1])
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (f != NULL);

  int i;
  for (i = 0; i < 4*(FTT_DIMENSION - 1); i++)
    f[i] = gfs_cell_corner_value (cell, corner[i], v, max_level);
  f[4*(FTT_DIMENSION - 1)] = GFS_VALUE (cell, v);
}

/**
 * gfs_interpolate_from_corners:
 * @cell: a #FttCell containing location @p.
 * @p: the location at which to interpolate.
 * @f: values at the corners of @cell.
 *
 * Returns: the interpolated value at location @p.
 */
gdouble gfs_interpolate_from_corners (FttCell * cell,
				      FttVector p,
				      gdouble * f)
{
  FttVector o;
  gdouble size;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (f != NULL, 0.);

  ftt_cell_pos (cell, &o);
  size = ftt_cell_size (cell)/2.;
  p.x = (p.x - o.x)/size;
  p.y = (p.y - o.y)/size;
#if FTT_2D
  gdouble x = (p.x + p.y)/2., y = (p.y - p.x)/2., v = f[4];
  if (x > 0.)
    v += x*(f[2] - f[4]);
  else
    v -= x*(f[0] - f[4]);
  if (y > 0.)
    v += y*(f[3] - f[4]);
  else
    v -= y*(f[1] - f[4]);
  return v;
#else  /* 3D */
  gdouble c[8];
  
  p.z = (p.z - o.z)/size;
  c[0] = - f[0] + f[1] + f[2] - f[3] - f[4] + f[5] + f[6] - f[7];
  c[1] = - f[0] - f[1] + f[2] + f[3] - f[4] - f[5] + f[6] + f[7];
  c[2] =   f[0] + f[1] + f[2] + f[3] - f[4] - f[5] - f[6] - f[7];
  c[3] =   f[0] - f[1] + f[2] - f[3] + f[4] - f[5] + f[6] - f[7];
  c[4] = - f[0] + f[1] + f[2] - f[3] + f[4] - f[5] - f[6] + f[7];
  c[5] = - f[0] - f[1] + f[2] + f[3] + f[4] + f[5] - f[6] - f[7];
  c[6] =   f[0] - f[1] + f[2] - f[3] - f[4] + f[5] - f[6] + f[7];
  c[7] =   f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];

  return (c[0]*p.x + c[1]*p.y + c[2]*p.z + 
	  c[3]*p.x*p.y + c[4]*p.x*p.z + c[5]*p.y*p.z + 
	  c[6]*p.x*p.y*p.z + 
	  c[7])/8.;
#endif /* 3D */  
}

/**
 * gfs_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the location at which to interpolate.
 * @v: a #GfsVariable.
 *
 * Interpolates the @v variable of @cell, at location @p. Linear
 * interpolation is used and the boundaries of the domain are treated
 * as planes of symmetry for all variables.
 *
 * Returns: the interpolated value of variable @v at location @p.
 */
gdouble gfs_interpolate (FttCell * cell,
			 FttVector p,
			 GfsVariable * v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  if (GFS_VALUE (cell, v) == GFS_NODATA)
    return GFS_NODATA;

  gdouble f[4*(FTT_DIMENSION - 1) + 1];
  gfs_cell_corner_values (cell, v, -1, f);
  return gfs_interpolate_from_corners (cell, p, f);
}

/**
 * gfs_interpolate_stencil:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Sets to 1. the @v variable of all the cells which would be used by
 * a call to gfs_interpolate().
 */
void gfs_interpolate_stencil (FttCell * cell,
			      GfsVariable * v)
{
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  for (i = 0; i < (FTT_DIMENSION == 2 ? 4 : 8); i++) {
    GfsInterpolator inter;
    guint j;

    gfs_cell_corner_interpolator (cell, corner[i], -1, TRUE, &inter);
    for (j = 0; j < inter.n; j++)
      GFS_VALUE (inter.c[j], v) = 1.;
  }
}

/**
 * gfs_center_curvature:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The curvature is normalized by the square of the size of the cell.
 *
 * Returns: the value of the @c component of the curvature of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_curvature (FttCell * cell,
			      FttComponent c,
			      guint v)
{
  FttCellFace f;
  GfsGradient g = { 0., 0. };

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  if (GFS_IS_MIXED (cell))
    return 0.;

  f.cell = cell;
  for (f.d = 2*c; f.d <= 2*c + 1; f.d++)
    if ((f.neighbor = ftt_cell_neighbor (cell, f.d))) {
      GfsGradient e;

      gfs_face_gradient (&f, &e, v, -1);
      g.a += e.a;
      g.b += e.b;
    }

  return g.b - g.a*GFS_VALUEI (cell, v);
}

/**
 * gfs_streamline_curvature:
 * @cell: a #FttCell.
 * @v: the components of the vector.
 *
 * The curvature is normalized by the size of the cell.
 *
 * Returns: the value of the curvature of the streamline defined by @v passing
 * through the center of the cell.
 */
gdouble gfs_streamline_curvature (FttCell * cell,
				  GfsVariable ** v)
{
  gdouble u2;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  u2 = gfs_vector_norm2 (cell, v);

  if (u2 == 0.)
    return 0.;
  else {
    FttComponent i;
    gdouble ugu = 0.;

    for (i = 0; i < FTT_DIMENSION; i++) {
      FttComponent j;
      gdouble ugui = 0.;

      for (j = 0; j < FTT_DIMENSION; j++)
	ugui += GFS_VALUE (cell, v[j])*gfs_center_gradient (cell, j, v[i]->i);
      ugu += ugui*ugui;
    }
    return sqrt (ugu)/u2;
  }
}

static FttCell * cell_corner_neighbor (FttCell * cell,
				       FttDirection * d,
				       gint max_level,
				       gboolean * t_junction)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, d[0]);
  if (!neighbor)
    return NULL;
  else {
    guint level = ftt_cell_level (cell);
    if (ftt_cell_level (neighbor) < level) {
      /* neighbor is at a shallower level */
      if (ftt_cell_child_corner (ftt_cell_parent (cell), d) != cell)
	*t_junction = TRUE;
      return neighbor;
    }
    else {
      if (level == max_level || FTT_CELL_IS_LEAF (neighbor))
	/* neighbor is at the same level */
	return neighbor;
      else {
	/* neighbor is at a deeper level */
	FttCell * n;
	guint i;
	FttDirection d1[FTT_DIMENSION];
	d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
	for (i = 1; i < FTT_DIMENSION; i++)
	  d1[i] = d[i];
	n = ftt_cell_child_corner (neighbor, d1);
	return n ? n : neighbor;
      }
    }
  }
}

static void interpolator_merge (GfsInterpolator * a, GfsInterpolator * b)
{
  guint i;
  for (i = 0; i < b->n; i++) {
    FttCell * c = b->c[i];
    guint j;

    for (j = 0; j < a->n && c != a->c[j]; j++)
      ;
    if (j < a->n)
      a->w[j] += b->w[i];
    else {
#if FTT_2D
      g_assert (j < 7);
#else
      g_assert (j < 29);
#endif
      a->c[j] = c;
      a->w[j] = b->w[i];
      a->n++;
    }
  }
}

static void interpolator_scale (GfsInterpolator * a, gdouble b)
{
  guint i;
  for (i = 0; i < a->n; i++)
    a->w[i] *= b;
}

static void t_junction_interpolator (FttCell * cell,
				     FttDirection * d,
				     FttCell * n,
				     gint max_level,
				     gboolean centered,
				     GfsInterpolator * inter)
{
  FttDirection d1[FTT_DIMENSION];
  GfsInterpolator a;

  d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
#if FTT_2D
  d1[1] = d[1];
  gfs_cell_corner_interpolator (n, d1, max_level, centered, inter);
  d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
  gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
  interpolator_merge (inter, &a);
  interpolator_scale (inter, 0.5);
#else /* 3D */
  d1[1] = d[1]; d1[2] = d[2];
  gfs_cell_corner_interpolator (n, d1, max_level, centered, inter);
  if (ftt_cell_neighbor_is_brother (cell, d[1])) {
    d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
    gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
    interpolator_merge (inter, &a);
    if (ftt_cell_neighbor_is_brother (cell, d[2])) {
      d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
      gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
      interpolator_merge (inter, &a);
      d1[1] = d[1];
      gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
      interpolator_merge (inter, &a);
      interpolator_scale (inter, 0.25);
    }
    else
      interpolator_scale (inter, 0.5);
  }
  else {
    d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
    interpolator_merge (inter, &a);
    interpolator_scale (inter, 0.5);
  }
#endif /* 3D */
}

static gboolean do_path (FttCell * cell, gint i,
			 FttCell * n[N_CELLS],
			 FttDirection * d,
			 gint max_level,
			 gboolean centered,
			 GfsInterpolator * inter)
{
  /* paths from each cell to neighbors. Cell indices are as in
     doc/figures/indices.fig
     first index: cell index
     second index: path index
     third index: < FTT_DIMENSION: directions.
                  = FTT_DIMENSION: neighboring cell index */
  static gint path[N_CELLS][FTT_DIMENSION][FTT_DIMENSION + 1] = {
#if FTT_2D
    {{1,2,1},   {2,1,2}},
    {{2,-1,3},  {-1,2,0}},
    {{1,-2,3},  {-2,1,0}},
    {{-1,-2,2}, {-2,-1,1}}
#else /* 3D */
    {{1,2,3,1},    {2,1,3,2},    {3,1,2,4}},
    {{2,3,-1,3},   {3,2,-1,5},   {-1,2,3,0}},
    {{1,-2,3,3},   {3,-2,1,6},   {-2,3,1,0}},
    {{3,-1,-2,7},  {-2,-1,3,1},  {-1,-2,3,2}},
    {{2,1,-3,6},   {1,2,-3,5},   {-3,1,2,0}},
    {{2,-1,-3,7},  {-1,2,-3,4},  {-3,2,-1,1}},
    {{1,-2,-3,7},  {-2,1,-3,4},  {-3,1,-2,2}},
    {{-1,-2,-3,6}, {-2,-1,-3,5}, {-3,-1,-2,3}}
#endif /* 3D */
  };
  guint j;

  for (j = 0; j < FTT_DIMENSION; j++) {
    guint k = path[i][j][FTT_DIMENSION];

    if (n[k] == NULL) {
      gboolean t_junction = FALSE;
      FttDirection d1[FTT_DIMENSION];
      guint l;

      for (l = 0; l < FTT_DIMENSION; l++)
	d1[l] = path[i][j][l] < 0 ? FTT_OPPOSITE_DIRECTION (d[- path[i][j][l] - 1]) : 
	  d[path[i][j][l] - 1];
      n[k] = cell_corner_neighbor (cell, d1, max_level, &t_junction);
      if (t_junction) {
	t_junction_interpolator (cell, d1, n[k], max_level, centered, inter);
	return TRUE;
      }
      if (n[k]) {
	t_junction = do_path (n[k], k, n, d, max_level, centered, inter);
	if (t_junction)
	  return TRUE;
      }
    }
  }
  return FALSE;
}

static gdouble distance (FttVector * c, FttCell * cell, gboolean centered)
{
  if (centered || !GFS_IS_MIXED (cell))
    return ftt_cell_size (cell)*
#if FTT_2D
      0.707106781185
#else  /* 3D */
      0.866025403785
#endif /* 3D */
      ;
  else {
    FttVector cm;
    gfs_cell_cm (cell, &cm);
#if FTT_2D
    return sqrt ((cm.x - c->x)*(cm.x - c->x) + (cm.y - c->y)*(cm.y - c->y));
#else /* 3D */
    return sqrt ((cm.x - c->x)*(cm.x - c->x) + (cm.y - c->y)*(cm.y - c->y) +
		 (cm.z - c->z)*(cm.z - c->z));
#endif /* 3D */
  }
}

/**
 * gfs_cell_corner_interpolator:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @centered: %TRUE if the interpolator is cell-centered.
 * @inter: a #GfsInterpolator.
 *
 * Fills @inter with the interpolator for the corner of @cell defined by @d.
 */
void gfs_cell_corner_interpolator (FttCell * cell,
				   FttDirection d[FTT_DIMENSION],
				   gint max_level,
				   gboolean centered,
				   GfsInterpolator * inter)
{
  FttCell * n[N_CELLS];
  guint i;
  gboolean t_junction;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (inter != NULL);

  while (!FTT_CELL_IS_LEAF (cell) && 
	 ftt_cell_level (cell) != max_level &&
	 (n[0] = ftt_cell_child_corner (cell, d)))
    cell = n[0];
  n[0] = cell;
  for (i = 1; i < N_CELLS; i++)
    n[i] = NULL;
  t_junction = do_path (cell, 0, n, d, max_level, centered, inter);
  if (t_junction)
    return;

  FttVector c;
  gdouble w = 0.;
  int boundaries = 0;
  
  inter->n = 0;
  ftt_corner_pos (cell, d, &c);
  for (i = 0; i < N_CELLS; i++)
    if (n[i]) {
      gdouble a;
      a = 1./(distance (&c, n[i], centered) + 1e-12);
      inter->c[inter->n] = n[i];
      inter->w[inter->n++] = a;
      w += a;
      if (GFS_CELL_IS_BOUNDARY (n[i]))
	boundaries++;
    }

  /* corners of domain receive special treatment */
  if (inter->n == FTT_DIMENSION + 1 && boundaries == FTT_DIMENSION) {
    /* remove central cell from interpolator */
    w -= inter->w[0];
    for (i = 0; i < inter->n - 1; i++) {
      inter->c[i] = inter->c[i + 1];
      inter->w[i] = inter->w[i + 1];
    }
    inter->n--;
  }

  g_assert (w > 0.);
  interpolator_scale (inter, 1./w);
}

/**
 * gfs_cell_corner_value:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @v: a #GfsVariable.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Returns: the value of variable @v interpolated at the corner of
 * @cell defined by @d.
 */
gdouble gfs_cell_corner_value (FttCell * cell,
			       FttDirection d[FTT_DIMENSION],
			       GfsVariable * v,
			       gint max_level)
{
  GfsInterpolator inter;
  gdouble val = 0.;
  guint i;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  gfs_cell_corner_interpolator (cell, d, max_level, v->centered, &inter);
  for (i = 0; i < inter.n; i++) {
    gdouble v1 = GFS_VALUE (inter.c[i], v);
    if (v1 == GFS_NODATA)
      return GFS_VALUE (cell, v);
    val += inter.w[i]*v1;
  }
  return val;
}

/* GfsStencil: Object */

/* The convention is that the first element is the diagonal */

/**
 * gfs_stencil_new:
 * @cell: the #FttCell associated with this stencil.
 * @lp: the linear problem.
 * @coeff: the cell coefficient.
 * 
 * Returns: a new #GfsStencil.
 */
GfsStencil * gfs_stencil_new (FttCell * cell, GfsLinearProblem * lp, gdouble coeff)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (lp != NULL, NULL);

  GfsStencil * s = g_malloc (sizeof (GfsStencil));
  s->id = g_array_new (FALSE, FALSE, sizeof (int));
  s->coeff = g_array_new (FALSE, FALSE, sizeof (double));
  gfs_stencil_add_element (s, cell, lp, coeff);
  return s;
}

/**
 * gfs_stencil_add_element:
 * @stencil: the stencil
 * @cell: a #FttCell.
 * @lp: the linear problem
 * @coeff: the cell coefficient.
 *
 * Adds the contribution of a @cell to @stencil.
 */
void gfs_stencil_add_element (GfsStencil * stencil, 
			      FttCell * cell,
			      GfsLinearProblem * lp,
			      gdouble coeff)
{
  gint i;

  g_return_if_fail (stencil != NULL);
  g_return_if_fail (lp != NULL);

  int id = GFS_VALUE (cell, lp->id);
  if (id < 0) {
    g_assert (GFS_CELL_IS_BOUNDARY (cell));
    FttCell * neighbor = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, lp->neighbor));
    g_assert (neighbor);
    id = GFS_VALUE (neighbor, lp->id);
    g_assert (id >= 0);
    coeff *= GFS_VALUE (cell, lp->neighborw);
  }

  for (i = 0; i < stencil->id->len; i++)
    if (g_array_index (stencil->id, int, i) == id) {
      g_array_index (stencil->coeff, double, i) += coeff;
      return;
    }

  g_array_append_val (stencil->id, id);
  g_array_append_val (stencil->coeff, coeff);
}

/**
 * gfs_stencil_destroy:
 * @stencil: a stencil 
 *
 * Destroys a stencil.
 */
void gfs_stencil_destroy (GfsStencil * stencil)
{
  g_return_if_fail (stencil != NULL);
  g_array_free (stencil->id, TRUE);
  g_array_free (stencil->coeff, TRUE);
  g_free (stencil);
}

static void get_average_neighbor_value_stencil (const FttCellFace * face,
						gdouble * x,
						GfsLinearProblem * lp, GfsStencil * stencil,
						gdouble weight)
{
  /* check for corner refinement violation (topology.fig) */
  g_assert (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell));
  
  if (FTT_CELL_IS_LEAF (face->neighbor))
    gfs_stencil_add_element (stencil, face->neighbor, lp, weight/(*x));
  else {
    FttCellChildren children;
    gdouble a = 0.;
    FttDirection od = FTT_OPPOSITE_DIRECTION (face->d);
    guint i, n;
    
    n = ftt_cell_children_direction (face->neighbor, od, &children);
    for (i = 0; i < n; i++)
      if (children.c[i]) {
	gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
	a += w;
      }
    if (a > 0.) {
      *x = 3./4.;
      for (i = 0; i < n; i++)
	if (children.c[i]) {
	  gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;

	  gfs_stencil_add_element (stencil, children.c[i], lp, weight*w/a/(*x));
	}
    }
    else
      gfs_stencil_add_element (stencil, face->cell, lp, weight/(*x));
  }
}

#if (FTT_2D || FTT_2D3)
/* v = a*v(cell) + b 
 * 
 * First order 1D interpolation.
 */
static GfsGradient get_interpolate_1D1_stencil (FttCell * cell,
						FttDirection d,
						gdouble x,
						GfsLinearProblem * lp,
						GfsStencil * stencil,
						gdouble weight)
{
  GfsGradient p;
  FttCellFace f;

  f = gfs_cell_face (cell, d);
  if (f.neighbor) {
    gdouble x2 = 1.;
    get_average_neighbor_value_stencil (&f, &x2, lp, stencil, weight*x);
    p.a = 1. - x/x2;
  }
  else
    p.a = 1.;

  return p;
}

#else

/* v = a*v(cell) + b 
 * 
 * First order 2D interpolation.
 */
static GfsGradient get_interpolate_2D1_stencil (FttCell * cell,
						FttDirection d1, FttDirection d2,
						gdouble x, gdouble y,
						GfsLinearProblem * lp, GfsStencil * stencil,
						gdouble weight)
{
  GfsGradient p;
  gdouble y1 = 1.;
  gdouble x2 = 1.;
  gdouble a1 = y, a2 = x;
  FttCellFace f1, f2;

  p.a = 1.;

  f1 = gfs_cell_face (cell, d1);
  if (f1.neighbor) {
    get_average_neighbor_value_stencil (&f1, &y1, lp, stencil, a1*weight);
    a1 /= y1;
    p.a -= a1;
  }
    
  f2 = gfs_cell_face (cell, d2);
  if (f2.neighbor) {
    get_average_neighbor_value_stencil (&f2, &x2, lp, stencil, a2*weight);
    a2 /= x2;
    p.a -= a2;
  }
  
  return p;
}

#endif /* not FTT_2D */

static Gradient gradient_fine_coarse_stencil (const FttCellFace * face,
					      GfsLinearProblem * lp,
					      GfsStencil * stencil, gdouble weight)
{
  Gradient g;
  GfsGradient p;
#if (FTT_2D || FTT_2D3)
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  g_assert (face != NULL);
  g_assert (ftt_face_type (face) == FTT_FINE_COARSE);

  dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if (FTT_2D || FTT_2D3)
  g_assert (dp >= 0);
  p = get_interpolate_1D1_stencil (face->neighbor, dp, 1./4., lp, stencil, 2./3.*weight);
#else  /* FTT_3D */
  g_assert (dp[0] >= 0 && dp[1] >= 0);
  p = get_interpolate_2D1_stencil (face->neighbor, dp[0], dp[1], 1./4., 1./4., lp,
				   stencil, 2./3.*weight);
#endif /* FTT_3D */

  g.a = 2./3.;
  g.b = 2.*p.a/3.;
  g.c = 0.;

  return g;
}

static void face_weighted_gradient_stencil (const FttCellFace * face,
					    GfsGradient * g,
					    gint max_level,
					    guint dimension,
					    GfsLinearProblem * lp,
					    GfsStencil * stencil)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = 0.;
  
  if (face->neighbor == NULL)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    Gradient gcf;
    gdouble w = GFS_STATE (face->cell)->f[face->d].v;

    gcf = gradient_fine_coarse_stencil (face, lp, stencil, w);
    g->a = w*gcf.a;
    gfs_stencil_add_element (stencil, face->neighbor, lp,  w*gcf.b); 
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      gdouble w = GFS_STATE (face->cell)->f[face->d].v;

      g->a = w;
      gfs_stencil_add_element (stencil, face->neighbor, lp, w);
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++)
	if ((f.cell = children.c[i])) {
	  Gradient gcf;
	  gdouble w = GFS_STATE (f.cell)->f[f.d].v;
	  if (dimension > 2)
	    w /=  n/2; /* fixme??? */
	
	  gcf = gradient_fine_coarse_stencil (&f, lp, stencil, -w);
	  g->a += w*gcf.b;
	  
	  gfs_stencil_add_element (stencil, f.cell, lp, w*gcf.a);
	}
    }
  }
}

/**
 * gfs_face_weighted_gradient_stencil:
 * @face: a FttCellFace
 * @g: a GfsGradient
 * @max_level: an integer
 * @lp: the linear problem
 * @stencil: a stencil 
 *
 * Fills @stencil with the stencil corresponding to the
 * weighted gradient on @face.
 */
void gfs_face_weighted_gradient_stencil (const FttCellFace * face,
					 GfsGradient * g,
					 gint max_level,
					 GfsLinearProblem * lp,
					 GfsStencil * stencil)
{
  face_weighted_gradient_stencil (face, g, max_level, FTT_DIMENSION, lp, stencil);
}

static gboolean mixed_face_gradient_stencil (const FttCellFace * face,
					     Gradient * g,
					     gint max_level,
					     GfsLinearProblem * lp,
					     GfsStencil * stencil,
					     gdouble weight)
{
  FttCell * n[N_CELLS];
  gdouble m[N_CELLS - 1][N_CELLS - 1];
  FttVector o, cm;
  FttComponent c = face->d/2;
  gdouble h = ftt_cell_size (face->cell);
  gdouble w = weight;

  gfs_cell_cm (face->cell, &o);
  if (!face_bilinear (face, n, &o, gfs_cell_cm, max_level, m))
    return FALSE;

  gfs_face_ca (face, &cm);

  if (!FTT_FACE_DIRECT (face))
    w = -w;

#if FTT_2D
  {
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
    gdouble vp;
    
    vp = ((&cm.x)[cp] - (&o.x)[cp])/h;
    g->a = ((m[c][0] + m[2][0]*vp) +
	    (m[c][1] + m[2][1]*vp) +
	    (m[c][2] + m[2][2]*vp));
    g->b = m[c][0] + m[2][0]*vp;
    gfs_stencil_add_element (stencil, n[2], lp, w*(m[c][1] + m[2][1]*vp));
    gfs_stencil_add_element (stencil, n[3], lp, w*(m[c][2] + m[2][2]*vp));
  }
#else /* 3D */
  {
    guint j;

    cm.x = (cm.x - o.x)/h;
    cm.y = (cm.y - o.y)/h;
    cm.z = (cm.z - o.z)/h;
    g->c = 0.;
    
    switch (c) {
    case FTT_X:
      g->a = g->b = m[0][0] + cm.y*m[3][0] + cm.z*m[4][0] + cm.y*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
  	gdouble a = m[0][j] + cm.y*m[3][j] + cm.z*m[4][j] + cm.y*cm.z*m[6][j];
  	g->a += a;
	gfs_stencil_add_element (stencil, n[j+1], lp, w*a);
      }
      break;
    case FTT_Y:
      g->a = g->b = m[1][0] + cm.x*m[3][0] + cm.z*m[5][0] + cm.x*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
  	gdouble a = m[1][j] + cm.x*m[3][j] + cm.z*m[5][j] + cm.x*cm.z*m[6][j];
  	g->a += a;
	gfs_stencil_add_element (stencil, n[j+1], lp, w*a);
      }
      break;
    case FTT_Z:
      g->a = g->b = m[2][0] + cm.x*m[4][0] + cm.y*m[5][0] + cm.x*cm.y*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
  	gdouble a = m[2][j] + cm.x*m[4][j] + cm.y*m[5][j] + cm.x*cm.y*m[6][j];
  	g->a += a;
	gfs_stencil_add_element (stencil, n[j+1], lp, w*a);
      }
      break;
    default:
      g_assert_not_reached ();
    }
  }
#endif /* 3D */

  if (!FTT_FACE_DIRECT (face)) {
    g->a = - g->a;
    g->b = - g->b;
  }

  return TRUE;
}

static void face_cm_gradient_stencil (const FttCellFace * face,
				      GfsGradient * g,
				      gint max_level,
				      gboolean weighted,
				      GfsLinearProblem * lp,
				      GfsStencil * stencil)
{
  guint level;
  Gradient gcf;
  gdouble w = weighted ? GFS_STATE (face->cell)->f[face->d].v : 1.;

  g->a = g->b = 0.;
  
  if (face->neighbor == NULL || w == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    if (GFS_IS_MIXED (face->cell) || GFS_IS_MIXED (face->neighbor)) {
      if (!mixed_face_gradient_stencil (face, &gcf, max_level, lp, stencil, w))
	gcf = gradient_fine_coarse_stencil (face, lp, stencil, w);
    }
    else
      gcf = gradient_fine_coarse_stencil (face, lp, stencil, w);
    g->a = w*gcf.a;
    gfs_stencil_add_element (stencil, face->neighbor, lp, w*gcf.b);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      if (!GFS_IS_MIXED (face->cell) && !GFS_IS_MIXED (face->neighbor)) {
  	g->a = w;
	gfs_stencil_add_element (stencil, face->neighbor, lp,  w);
      }
      else if (mixed_face_gradient_stencil (face, &gcf, max_level, lp, stencil, w)) {
  	g->a = w*gcf.a;
	gfs_stencil_add_element (stencil, face->neighbor, lp,  w*gcf.b);
      }
      else {
  	g->a = w;
	gfs_stencil_add_element (stencil, face->neighbor, lp,  w);
      }
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++)
    	if ((f.cell = children.c[i])) {
    	  w = weighted ? GFS_STATE (f.cell)->f[f.d].v : 1.;
    	  if (GFS_IS_MIXED (f.cell) || GFS_IS_MIXED (f.neighbor)) {
    	    if (!mixed_face_gradient_stencil (&f, &gcf, max_level, lp, stencil, w))
    	      gcf = gradient_fine_coarse_stencil (&f, lp, stencil, -w);
    	  }
    	  else
    	    gcf = gradient_fine_coarse_stencil (&f, lp, stencil, -w);
    	  g->a += w*gcf.b;
	  gfs_stencil_add_element (stencil, f.cell, lp,  w*gcf.a);
    	}
    }
  }
}

/**
 * gfs_face_cm_weighted_gradient_stencil:
 * @face: a #FttCellFace.
 * @g: a GfsGradient.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @lp: the linear problem.
 * @stencil: a stencil.
 *
 * Fills @stencil with the stencil corresponding to the
 * weighted gradient on @face. The weights are defined at
 * the center of mass of the cell. Linear interpolation is used
 * to evaluate the gradient in the vicinity of cut cells.
 * This function returns the stencil equivalent to the action of
 * gfs_face_cm_weighted_gradient().
 */
void gfs_face_cm_weighted_gradient_stencil (const FttCellFace * face,
					    GfsGradient * g,
					    gint max_level,
					    GfsLinearProblem * lp,
					    GfsStencil * stencil)
{
  g_return_if_fail (face != NULL);
  g_return_if_fail (g != NULL);
  g_return_if_fail (lp != NULL);
  g_return_if_fail (stencil != NULL);

  face_cm_gradient_stencil (face, g, max_level, TRUE, lp, stencil);
}

/**
 * gfs_cm_gradient:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @g: a #FttVector.
 *
 * Fills @g with the components of the gradient of @v at the center of
 * mass of @cell.
 *
 * The gradient is normalized by the size of the cell.
 */
void gfs_cm_gradient (FttCell * cell,
		      GfsVariable * v,
		      FttVector * g)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (g != NULL);

  if (v->centered) {
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&g->x)[c] = gfs_center_gradient (cell, c, v->i);
  }
  else if (GFS_IS_MIXED (cell))
    gfs_mixed_cell_gradient (cell, v, g);
  else {
    FttComponent c;
    GfsGradient g1, g2;
    FttCellFace f;
    FttCellNeighbors n;
    gdouble val = GFS_VALUE (cell, v);
    ftt_cell_neighbors (cell, &n);
    f.cell = cell;
    for (c = 0; c < FTT_DIMENSION; c++) {
      f.d = 2*c; f.neighbor = n.c[f.d];
      gfs_face_cm_gradient (&f, &g1, v->i, -1);
      f.d = 2*c + 1; f.neighbor = n.c[f.d];
      gfs_face_cm_gradient (&f, &g2, v->i, -1);
      (&g->x)[c] = (g1.b - g2.b + (g2.a - g1.a)*val)/2.;
    }
  }
}
