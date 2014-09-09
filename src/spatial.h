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

#ifndef __SPATIAL_H__
#define __SPATIAL_H__

/* When modifying this file, please also update the documentation at:
 * http://gfs.sourceforge.net/wiki/index.php/GfsSurface
 */

#define intersection(a,b) MAX(a, b)
#define union(a,b)        MIN(a, b)
#define difference(a,b)   MAX(a, -(b))

#define X _x
#define Y _y
#define Z _z
#define rx _x
#define ry _y
#define rz _z

static double _x = 0., _y = 0., _z = 0.;

static double ellipse (double xc, double yc, double a, double b)
{
  g_return_val_if_fail (a != 0. && b != 0., 0.);
  return (_x - xc)*(_x - xc)/(a*a) + (_y - yc)*(_y - yc)/(b*b) - 1.;
}

static double sphere (double xc, double yc, double zc, double r)
{
  return (_x - xc)*(_x - xc) + (_y - yc)*(_y - yc) + (_z - zc)*(_z - zc) - r*r;
}

static double cube (double xc, double yc, double zc, double h)
{
  double vmax = (_x - xc)*(_x - xc) - h*h/4.;
  double v = (_y - yc)*(_y - yc) - h*h/4.;
  if (v > vmax) vmax = v;
  v = (_z - zc)*(_z - zc) - h*h/4.;
  if (v > vmax) vmax = v;
  return vmax;
}

#endif /* __SPATIAL_H__ */
