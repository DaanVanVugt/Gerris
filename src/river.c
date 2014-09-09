/* Gerris - The GNU Flow Solver
 * Copyright (C) 2008-2012 National Institute of Water and Atmospheric Research
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

/*
 * Relevant references:
 *
 * Saint-Venant:
 * 
 * [Audusse2005] E. Audusse and M.-O. Bristeau. A well-balanced,
 * positivity-preserving second-order scheme for shallow-water flows
 * on unstructured meshes, JCP, 2005, 311-333.
 *
 * [Popinet2011] S. Popinet. Quadtree-adaptive tsunami modelling. Ocean Dynamics
 * 61(9):1261-1285, 2011.
 *
 * [An2012] Hyunuk An, Soonyoung Yu. Well-balanced shallow water flow
 * simulation on quadtree cut cell grids. Advances in Water Resources
 * 39:60-70, 2012.
 *
 * Multi-layer Saint-Venant, constant density:
 *
 * [Audusse2011a] E. Audusse, M.-O. Bristeau, B. Perthame and J. Sainte-Marie. A
 * multilayer Saint-Venant system with mass exchanges for
 * shallow-water flows. Derivation and numerical
 * validation. Mathematical Modelling and Numerical analysis, 2011.
 * 
 * Multi-layer Saint-Venant, variable density:
 *
 * [Audusse2011b] E. Audusse, M.-O. Bristeau, M. Pelanti,
 * J. Sainte-Marie. Approximation of the hydrostatic Navier-Stokes
 * system for density stratified flows by a multilayer model. Kinetic
 * interpretation and numerical solution, JCP, 2011.
 */

/*! \file
 * \brief GfsRiver model.
 */

#include <stdlib.h>
#include "river.h"
#include "adaptive.h"
#include "source.h"
#include "solid.h"
#include "init.h"

/* generalisation of the limited gradients (in fluid.c) to mixed cells */

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

static gdouble center_limited_gradient_full (FttCell * cell,
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

static gdouble center_limited_gradient (FttCell * cell,
					FttComponent c,
					guint v,
					gdouble (* limiter) (gdouble))
{
  FttDirection d = 2*c;
  FttCellFace f1, f2;
  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  f2 = gfs_cell_face (cell, d);
  if (!GFS_IS_MIXED (cell) && 
      (!f1.neighbor || !GFS_IS_MIXED (f1.neighbor)) &&
      (!f2.neighbor || !GFS_IS_MIXED (f2.neighbor)))
    return center_limited_gradient_full (cell, c, v, limiter);
  gdouble h = ftt_cell_size (cell);
  FttVector cm;
  gfs_cell_cm (cell, &cm);  
  gdouble v0 = GFS_VALUEI (cell, v), g = 0.;
  
  if (f1.neighbor && f2.neighbor) {
    /* two neighbors */
    gdouble x1, x2;
    gdouble v1 = gfs_neighbor_value (&f1, v, &x1);
    gdouble v2 = gfs_neighbor_value (&f2, v, &x2);
    if (v0 != v1) {
      FttVector cm1, cm2;
      gfs_cell_cm (f1.neighbor, &cm1);
      gfs_cell_cm (f2.neighbor, &cm2);       
      /* fixme: this is not correct at coarse/fine boundaries */
      x1 = ((&cm.x)[c] - (&cm1.x)[c])/h;
      x2 = ((&cm2.x)[c] - (&cm.x)[c])/h;
      g = (* limiter) ((v2 - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
    }
  }
 
  /* mixed cells gradient following Causon et al. (2000) */
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector ca = s->ca;
    FttVector n;
    gdouble nn;

    gfs_solid_normal (cell, &n);
    nn = sqrt (n.x*n.x + n.y*n.y);
    n.x /= nn;
    n.y /= nn;
    
    /* solid is on the right side of the cell */
    if (s->s[2*c] < s->s[2*c + 1]) {
      if (f1.neighbor) {
	gdouble vr;
	/* fixme: this relies on specific indices for U and V. Not recommended. */
	if (v == 2) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.x;
	else if (v == 3) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.y;
	else return s->s[2*c]*g/s->s[2*c + 1];
	gdouble x1, v1 = gfs_neighbor_value (&f1, v, &x1);
	FttVector cm1;
	gfs_cell_cm (f1.neighbor, &cm1);
	/* fixme: this is not correct at coarse/fine boundaries */
	x1 = ((&cm.x)[c] - (&cm1.x)[c])/h;
	gdouble x2 = 2.*((&ca.x)[c] - (&cm.x)[c])/h;
	gdouble gs = ((v0 - v1)*x2 == 0. || x1 == 0.) ? 0. :
	  (* limiter) ((vr - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
	return (s->s[2*c]*g + (s->s[2*c + 1] - s->s[2*c])*gs)/s->s[2*c + 1];
      }
      else 
	return 0;
    }    
    /* solid is on the left side of the cell */
    else if (s->s[2*c] > s->s[2*c + 1]) {
      if (f2.neighbor) {
	gdouble vr;
	/* fixme: this relies on specific indices for U and V. Not recommended. */
	if (v == 2) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.x; 
	else if (v == 3) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.y; 
	else return s->s[2*c + 1]*g/s->s[2*c];
	gdouble x2, v2 = gfs_neighbor_value (&f2, v, &x2);
	FttVector cm2;
	gfs_cell_cm (f2.neighbor, &cm2);
	/* fixme: this is not correct at coarse/fine boundaries */
	gdouble x1 = 2.*((&cm.x)[c] - (&ca.x)[c])/h;
	x2 = ((&cm2.x)[c] - (&cm.x)[c])/h;
 	gdouble gs = ((v0 - vr)*x2 == 0. || x1 == 0.) ? 0. :
	  (* limiter) ((v2 - v0)*x1/((v0 - vr)*x2))*(v0 - vr)/x1;
	return (s->s[2*c + 1]*g + (s->s[2*c] - s->s[2*c + 1])*gs)/s->s[2*c];
      }
      else 
	return 0;
    }
  }
  return g;
}

static gdouble center_minmod_gradient (FttCell * cell,
				       FttComponent c,
				       guint v)
{
  return center_limited_gradient (cell, c, v, minmod_limiter);
}

static gdouble center_superbee_gradient (FttCell * cell,
					 FttComponent c,
					 guint v)
{
  return center_limited_gradient (cell, c, v, superbee_limiter);
}

static gdouble center_sweby_gradient (FttCell * cell,
				      FttComponent c,
				      guint v)
{
  return center_limited_gradient (cell, c, v, sweby_limiter);
}

/**
 * Solves the Saint-Venant equations.
 * \beginobject{GfsRiver}
 */

#define H  0
#define U  1
#define V  2
#define HL 3
#define T(i,l) ((2 + (i))*r->nlayers + 1 + (l))
#define RHO(l) ((2 + r->nt)*r->nlayers + 3*(l) + 1)
#define HP(l)  ((2 + r->nt)*r->nlayers + 3*(l) + 2)
#define HPT(l) ((2 + r->nt)*r->nlayers + 3*(l) + 3)
#define ZB (r->nvar)

static void flux (const gdouble * u, gdouble g, gdouble * f)
{
  f[H] = u[H]*u[U];                       /* h*u */
  f[U] = u[H]*(u[U]*u[U] + g*u[H]/2.);    /* h*(u*u + g*h/2) */
  f[V] = u[H]*u[U]*u[V];                  /* h*u*v */
}

static gdouble min (gdouble a, gdouble b)
{
  return a < b ? a : b;
}

static gdouble max (gdouble a, gdouble b)
{
  return a > b ? a : b;
}

/*
 * uL: left state vector [h,u,v,zb].
 * uR: right state vector.
 * g: acceleration of gravity.
 * f: flux vector.
 *
 * Fills @f by solving an approximate Riemann problem using the HLLC
 * scheme. See e.g. Liang, Borthwick, Stelling, IJNMF, 2004.
 */
static void riemann_hllc (const GfsRiver * r,
			  const gdouble * uL, const gdouble * uR,
			  gdouble * f)
{
  gdouble cL = sqrt (r->g*uL[H]), cR = sqrt (r->g*uR[H]);
  gdouble ustar = (uL[U] + uR[U])/2. + cL - cR;
  gdouble cstar = (cL + cR)/2. + (uL[U] - uR[U])/4.;
  gdouble SL = uL[H] == 0. ? uR[U] - 2.*cR : min (uL[U] - cL, ustar - cstar);
  gdouble SR = uR[H] == 0. ? uL[U] + 2.*cL : max (uR[U] + cR, ustar + cstar);

  if (0. <= SL)
    flux (uL, r->g, f);
  else if (0. >= SR)
    flux (uR, r->g, f);
  else {
    gdouble fL[3], fR[3];
    flux (uL, r->g, fL);
    flux (uR, r->g, fR);
    f[H] = (SR*fL[H] - SL*fR[H] + SL*SR*(uR[H] - uL[H]))/(SR - SL);
    f[U] = (SR*fL[U] - SL*fR[U] + SL*SR*(uR[H]*uR[U] - uL[H]*uL[U]))/(SR - SL);
    gdouble SM = ((SL*uR[H]*(uR[U] - SR) - SR*uL[H]*(uL[U] - SL))/
		  (uR[H]*(uR[U] - SR) - uL[H]*(uL[U] - SL)));
    if (SL <= 0. && 0. <= SM)
      f[V] = uL[V]*f[H];
    else if (SM <= 0. && 0. <= SR)
      f[V] = uR[V]*f[H];
    else {
      fprintf (stderr, "L: %g %g %g R: %g %g %g\n",
	       uL[H], uL[U], uL[V],
	       uR[H], uR[U], uR[V]);
      fprintf (stderr, "SL: %g SR: %g SM: %g\n", SL, SR, SM);
      g_assert_not_reached ();
    }
  }
}

/*
 * uL: left state vector [h,u,v,zb].
 * uR: right state vector.
 * g: acceleration of gravity.
 * f: flux vector.
 *
 * Fills @f by solving an approximate Riemann problem using the kinetic
 * scheme. See Audusse2011b and Audusse2005.
 */

#define SQRT3 1.73205080756888

#define PARENT_TRACER(v) ((v)->vector[0]) /* hack: use vector[0] to store parent tracer */

static double density (GfsRiver * r, int l, const gdouble * u, FttCell * cell)
{
  /* set parent tracer values to values for this level (stored in u) */
  GfsSimulation * sim = GFS_SIMULATION (r);
  for (int i = 0; i < r->nt; i++)
    GFS_VALUE (cell, PARENT_TRACER (r->v[T(i,l)])) = u[T(i,l)]/sim->physical_params.L;
  /* evaluate density for this level */
  r->l = l;
  return 1./gfs_function_value (sim->physical_params.alpha, cell);
}

static void hydrostatic_pressure (GfsRiver * r, double * u, FttCell * cell)
{
  double pa = 0.; /* pressure at the top */
  u[HPT(r->nlayers)] = pa;
  for (int l = r->nlayers - 1; l >= 0; l--) {
    u[RHO(l)] = density (r, l, u, cell);
    g_assert (u[RHO(l)] > 0.);
    double dp = r->g*u[RHO(l)]*u[H]*r->dz[l];
    u[HP(l)] = pa + dp/2.; /* midlayer pressure i.e. p_\alpha */
    pa += dp;
    u[HPT(l)] = pa; /* pressure at bottom of layer i.e. p_\alpha-1/2 */
    /* Boussinesq */
    u[RHO(l)] = 1.;
  }
}

static void riemann_kinetic (const GfsRiver * r,
			     const gdouble * uL, const gdouble * uR,
			     gdouble * f)
{
  gdouble * dz = r->dz;
  f[H] = 0.;
  int l;
  for (l = 0; l < r->nlayers; l++) {
    gdouble ci, Mp, Mm, cig, fHl;

    if (uL[H] > r->dry) {
      if (r->variable_density)
	ci = sqrt (uL[HP(l)]/uL[RHO(l)]);
      else
	ci = sqrt (r->g*uL[H]/2.);
      Mp = MAX (uL[U + 2*l] + ci*SQRT3, 0.);
      Mm = MAX (uL[U + 2*l] - ci*SQRT3, 0.);
      if (r->variable_density)
	cig = dz[l]*uL[H]/(12.*SQRT3*ci);
      else
	cig = dz[l]*ci/(6.*r->g*SQRT3);
      fHl = cig*3.*(Mp*Mp - Mm*Mm);
      f[U + 3*l] = cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
    }
    else
      fHl = f[U + 3*l] = 0.;

    if (uR[H] > r->dry) {
      if (r->variable_density)
	ci = sqrt (uR[HP(l)]/uR[RHO(l)]);
      else
	ci = sqrt (r->g*uR[H]/2.);
      Mp = MIN (uR[U + 2*l] + ci*SQRT3, 0.);
      Mm = MIN (uR[U + 2*l] - ci*SQRT3, 0.);
      if (r->variable_density)
	cig = dz[l]*uR[H]/(12.*SQRT3*ci);
      else
	cig = dz[l]*ci/(6.*r->g*SQRT3);
      fHl += cig*3.*(Mp*Mp - Mm*Mm);
      f[U + 3*l] += cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
    }

    f[V + 3*l] = (fHl > 0. ? uL[V + 2*l] : uR[V + 2*l])*fHl;

    f[HL + 3*l] = fHl;
    f[H] += fHl;
  }
}

typedef struct {
  FttComponent u;
  gdouble du;
  FttComponent v;
  gdouble dv;
} Sym;

#define CFL_CLAMP(u, umax) (fabs (u) <= (umax) ? (u) : (u) > 0. ? (umax) : - (umax))

static double left (const GfsRiver * r, const FttCellFace * f, int i, double a)
{
  return GFS_VALUE (f->cell, r->v1[i]) + a*GFS_VALUE (f->cell, r->dv[f->d/2][i]);
}

static double right (const GfsRiver * r, const FttCellFace * f, int i, double a)
{
  return GFS_VALUE (f->neighbor, r->v1[i]) - a*GFS_VALUE (f->neighbor, r->dv[f->d/2][i]);
}

static void face_fluxes (FttCellFace * face, GfsRiver * r)
{
  gdouble eta = GFS_VALUE (face->cell, r->v1[H]), etan = GFS_VALUE (face->neighbor, r->v1[H]);

  if (eta <= r->dry && etan <= r->dry)
    return;

  gdouble a = 1., b = 1.;
  if (GFS_IS_MIXED (face->cell)) {
    FttVector ca, cm;
    gfs_face_ca (face, &ca);
    gfs_cell_cm (face->cell, &cm);
    FttComponent c = face->d/2;
    a = fabs (2.*((&ca.x)[c] - (&cm.x)[c])/ftt_cell_size (face->cell));
  }
  if (GFS_IS_MIXED (face->neighbor)) {
    FttVector ca, cm;
    gfs_face_ca (face, &ca); /* fixme?: this is not symmetric with the above for face->cell */
    gfs_cell_cm (face->neighbor, &cm);
    FttComponent c = face->d/2;
    b = fabs (2.*((&ca.x)[c] - (&cm.x)[c])/ftt_cell_size (face->neighbor));
  } 

  static Sym sym[4] = {
    {U,  1., V,  1.},
    {U, -1., V, -1.},
    {V,  1., U, -1.},
    {V, -1., U,  1.}
  };
  Sym * s = &sym[face->d];

  gdouble etaL = (eta <= r->dry ? 0. : left (r, face, H, a*s->du));
  gdouble zbL = (GFS_VALUE (face->cell, r->zb)
		 + a*s->du*GFS_VALUE (face->cell, r->dv[face->d/2][ZB]));
  gdouble zbR = (GFS_VALUE (face->neighbor, r->zb)
		 - b*s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][ZB])); 
  gdouble zbLR = MAX (zbL, zbR);

  gdouble * uL = r->uL, * uR = r->uR;
  int l;
  if (etaL > r->dry)
    for (l = 0; l < r->nlayers; l++) {
      gdouble etal = etaL*r->dz[l];
      /* ul = uhl/hl */
      uL[U + 2*l] = s->du*left (r, face, s->u + 2*l, a*s->du)/etal;
      /* vl = vhl/hl */
      uL[V + 2*l] = s->dv*left (r, face, s->v + 2*l, a*s->du)/etal;
      /* tl = thl/hl */
      for (int i = 0; i < r->nt; i++)
	uL[T(i,l)] = left (r, face, T(i,l), a*s->du)/etal;
    }
  else
    for (l = 0; l < r->nlayers; l++) {
      uL[U + 2*l] = uL[V + 2*l] = 0.;
      for (int i = 0; i < r->nt; i++)
	uL[T(i,l)] = 0.; /* fixme! */
    }
  uL[H] = MAX (0., etaL + zbL - zbLR);

  gdouble etaR = (etan <= r->dry ? 0. : right (r, face, H, b*s->du));
  /* fixme: this is only first-order accurate for fine/coarse */
  if (etaR > r->dry)
    for (l = 0; l < r->nlayers; l++) {
      gdouble etal = etaR*r->dz[l];
      /* ul = uhl/hl */
      uR[U + 2*l] = s->du*right (r, face, s->u + 2*l, b*s->du)/etal;
      /* vl = vhl/hl */
      uR[V + 2*l] = s->dv*right (r, face, s->v + 2*l, b*s->du)/etal;
      /* tl = thl/hl */
      for (int i = 0; i < r->nt; i++)
	uR[T(i,l)] = right (r, face, T(i,l), a*s->du)/etal;
    }
  else
    for (l = 0; l < r->nlayers; l++) {
      uR[U + 2*l] = uR[V + 2*l] = 0.;
      for (int i = 0; i < r->nt; i++)
	uR[T(i,l)] = 0.; /* fixme! */
    }
  uR[H] = MAX (0., etaR + zbR - zbLR);

  gdouble h = ftt_cell_size (face->cell);
  gdouble umax = GFS_SIMULATION (r)->advection_params.cfl*h/r->dt;
  for (l = 0; l < r->nlayers; l++) {
    uL[U + 2*l] = CFL_CLAMP (uL[U + 2*l], umax);
    uR[U + 2*l] = CFL_CLAMP (uR[U + 2*l], umax);
    uL[V + 2*l] = CFL_CLAMP (uL[V + 2*l], umax);
    uR[V + 2*l] = CFL_CLAMP (uR[V + 2*l], umax);
  }

  gdouble * u, * un;
  if (r->variable_density) {
    hydrostatic_pressure (r, uL, face->cell);
    hydrostatic_pressure (r, uR, face->cell);

    u = g_malloc ((r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));
    un = g_malloc ((r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));
    u[H] = eta; un[H] = etan;
    for (l = 0; l < r->nlayers; l++) {
      gdouble etal = eta*r->dz[l];
      gdouble etanl = etan*r->dz[l];
      for (int i = 0; i < r->nt; i++) {
	u[T(i,l)] = etal > 0. ? GFS_VALUE (face->cell, r->v1[T(i,l)])/etal : 0.; /* fixme! */
	un[T(i,l)] = etanl > 0. ? GFS_VALUE (face->neighbor, r->v1[T(i,l)])/etanl : 0.; /* fixme! */
      }
    }
    hydrostatic_pressure (r, u, face->cell);
    hydrostatic_pressure (r, un, face->cell);
  }
  else
    u = un = NULL;

  /* Riemann solver */
  gdouble * f = r->f;
  (* r->scheme) (r, uL, uR, f);

  gdouble dt = gfs_domain_face_fraction (GFS_DOMAIN (r), face)*r->dt/h;
  GFS_VALUE (face->cell, r->flux[H]) -= dt*f[H];
  gdouble nn = (ftt_face_type (face) == FTT_FINE_COARSE ? FTT_CELLS : 1.);
  GFS_VALUE (face->neighbor, r->flux[H]) += dt*f[H]/nn;

  gdouble zb = GFS_VALUE (face->cell, r->zb);
  gdouble zbn = GFS_VALUE (face->neighbor, r->zb);
  gdouble SbL0, SbR0, SbL, SbR;
  if (r->variable_density) {
    /* eq. (89) of Audusse2011b (corrected) */
    SbL0 = (uL[HPT(0)] + u[HPT(0)])*(zbLR - zb)/2.;
    SbR0 = (uR[HPT(0)] + un[HPT(0)])*(zbLR - zbn)/2.;
    SbL = SbR = 0.;
  }
  else { /* constant density */
    if (eta <= r->dry) eta = 0.;
    if (etan <= r->dry) etan = 0.;
    /* see Audusse2005, equations 4.4 and 5.13 */
    /* Slope source term "S_{i,j,p}" and second-order correction for
       slope source term "Sc_{i,j,p}" of An2012, equations (11) and (12) */
    SbL = r->g/2.*(uL[H]*uL[H] - etaL*etaL - (etaL + eta)*(zbL - zb));
    SbR = r->g/2.*(uR[H]*uR[H] - etaR*etaR - (etaR + etan)*(zbR - zbn));
    SbL0 = SbR0 = 0.;
  }

  gdouble G = 0.;
  for (l = 0; l < r->nlayers; l++) {
    gdouble dz = r->dz[l];
    if (r->variable_density) {
      /* eq. (84-85) of Audusse2011b (corrected) */
      zbLR += (uL[H] + uR[H])*dz/2.;
      zb += u[H]*dz; zbn += un[H]*dz;
      /* eq. (89) of Audusse2011b (corrected) */
      gdouble SbL = (uL[HPT(l + 1)] + u[HPT(l + 1)])*(zbLR - zb)/2.;
      gdouble SbR = (uR[HPT(l + 1)] + un[HPT(l + 1)])*(zbLR - zbn)/2.;
      GFS_VALUE (face->cell, r->flux[s->u + 2*l]) -= s->du*dt*(
							       f[U + 3*l] 
							       - SbL + SbL0
							       );
      GFS_VALUE (face->neighbor, r->flux[s->u + 2*l]) += s->du*dt*(
								   f[U + 3*l] 
								   - SbR + SbR0
								   )/nn;
      SbL0 = SbL;
      SbR0 = SbR;
    }
    else { /* constant density */
      GFS_VALUE (face->cell, r->flux[s->u + 2*l]) -= s->du*dt*(f[U + 3*l] - dz*SbL);
      GFS_VALUE (face->neighbor, r->flux[s->u + 2*l]) += s->du*dt*(f[U + 3*l] - dz*SbR)/nn;
    }
    GFS_VALUE (face->cell, r->flux[s->v + 2*l]) -= s->dv*dt*f[V + 3*l];
    GFS_VALUE (face->neighbor, r->flux[s->v + 2*l]) += s->dv*dt*f[V + 3*l]/nn;
    /* mass flux between layers */
    G += dt*(f[HL + 3*l] - dz*f[H]); /* eq. (5.109) of Audusse2011a and (75) of Audusse2011b */
    if (l < r->nlayers - 1) {
      GFS_VALUE (face->cell, r->massflux[l]) += G;
      GFS_VALUE (face->neighbor, r->massflux[l]) -= G/nn;
    }
    /* horizontal tracer advection */
    for (int i = 0; i < r->nt; i++) {
      double flux = dt*f[HL + 3*l];
      flux *= flux > 0. ? uL[T(i,l)] : uR[T(i,l)];
      GFS_VALUE (face->cell, r->flux[T(i,l)]) -= flux;
      GFS_VALUE (face->neighbor, r->flux[T(i,l)]) += flux/nn;
    }
  }

  if (r->variable_density) {
    g_free (u);
    g_free (un);
  }
}

static gdouble limited_gradient (const FttCell * cell, const GfsRiver * r, 
				 int i0, int i1, int i2, 
				 int l,
				 gdouble (* limiter) (gdouble))
{
  if (l < 1 || l > r->nlayers - 2)
    return 0.;
  gdouble v0 = GFS_VALUE (cell, r->v1[i0]);
  gdouble v1 = GFS_VALUE (cell, r->v1[i1]);
  if (v0 == v1)
    return 0.;

  gdouble x1 = (r->dz[l] + r->dz[l - 1])/(2.*r->dz[l]);
  gdouble v2 = GFS_VALUE (cell, r->v1[i2]);
  gdouble x2 = (r->dz[l] + r->dz[l + 1])/(2.*r->dz[l]);
  return (* limiter) ((v2 - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
}

/* fixme: minmod limiter for vertical advection may be too diffusive */
#define limited_gradient_u(l) limited_gradient(cell, r, c+2*(l), c+2*((l)-1), c+2*((l)+1), l, \
					       minmod_limiter)
#define limited_gradient_t(l) limited_gradient(cell, r, T(i,l), T(i,(l)-1), T(i,(l)+1), l, \
					       minmod_limiter)

static void vertical_advection (FttCell * cell, const GfsRiver * r)
{
  double eta = GFS_VALUE (cell, r->v1[H]);
  if (eta > r->dry)
    for (int l = 0; l < r->nlayers - 1; l++) {
      double dz = eta*(r->dz[l] + r->dz[l + 1])/2.;
      double G = GFS_VALUE (cell, r->massflux[l])/dz;
      for (FttComponent c = U; c <= V; c++) {
	double flux = G < 0. ?
	  G*(GFS_VALUE (cell, r->v1[c + 2*l]) + limited_gradient_u (l)/2.) :
	  G*(GFS_VALUE (cell, r->v1[c + 2*(l + 1)]) - limited_gradient_u (l + 1)/2.);
	GFS_VALUE (cell, r->flux[c + 2*l]) += flux;
	GFS_VALUE (cell, r->flux[c + 2*(l + 1)]) -= flux;
      }
      for (int i = 0; i < r->nt; i++) {
	double flux = G < 0. ? 
	  G*(GFS_VALUE (cell, r->v1[T(i,l)]) +  limited_gradient_t (l)/2.) :
	  G*(GFS_VALUE (cell, r->v1[T(i,l + 1)]) - limited_gradient_t (l + 1)/2.);
	GFS_VALUE (cell, r->flux[T(i,l)]) += flux;
	GFS_VALUE (cell, r->flux[T(i,l + 1)]) -= flux;
      }
    }
}

static void reset_fluxes (FttCell * cell, const GfsRiver * r)
{
  int i;
  for (i = 0; i < r->nvar; i++)
    GFS_VALUE (cell, r->flux[i]) = 0.;
  for (i = 0; i < r->nlayers - 1; i++)
    GFS_VALUE (cell, r->massflux[i]) = 0.;
}

static void solid_boundary_fluxes (FttCell * cell, GfsRiver * r)
{
  gdouble h = ftt_cell_size (cell);
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttVector cm;

  gdouble hh = MAX (GFS_VALUE (cell, r->v1[H]), 0.);
  gfs_cell_cm (cell, &cm);
  gdouble hs = hh + 2.*((s->ca.x - cm.x)*GFS_VALUE (cell, r->dv[0][H]) + 
			(s->ca.y - cm.y)*GFS_VALUE (cell, r->dv[1][H]))/h;
  gdouble zbs = 2.*((s->ca.x - cm.x)*GFS_VALUE (cell, r->dv[0][ZB]) +
		    (s->ca.y - cm.y)*GFS_VALUE (cell, r->dv[1][ZB]))/h;
  gdouble hszbs = 
    r->dt/h*r->g/2.*(
		     /* the normal component of the velocity is zero at the solid boundary. */  
		     hs*hs +
		     /* Second-order correction for slope source term ("Sc_{i,j,s}" of 
			An2012, equation (27) */
		     (hs + hh)*zbs);
  FttVector n;
  gfs_solid_normal (cell, &n);
  GFS_VALUE (cell, r->flux[U]) -= n.x*hszbs;
  GFS_VALUE (cell, r->flux[V]) -= n.y*hszbs;
}

/* Metric source terms (see doc/figures/lonlat.tm) */
static void metric_sources (FttCell * cell, GfsRiver * r)
{
  if (GFS_VALUE (cell, r->v1[H]) > r->dry) {
    /* fixme: this will probably not work when combined with solids */
    GfsDomain * domain = GFS_DOMAIN (r);
    gdouble fm[FTT_NEIGHBORS], cm;
    FttCellFace face = { cell };
    for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++)
      fm[face.d] = (* domain->face_metric) (domain, &face);
    gdouble dh_dl = fm[FTT_RIGHT] - fm[FTT_LEFT];
    gdouble dh_dt = fm[FTT_TOP]   - fm[FTT_BOTTOM];
    cm = (* domain->cell_metric) (domain, cell)*ftt_cell_size (cell);
    gdouble dldh = cm*GFS_SIMULATION (r)->physical_params.L;
    gdouble 
      phiu = GFS_VALUE (cell, r->v1[U]), 
      phiv = GFS_VALUE (cell, r->v1[V]);
    gdouble fG = phiv*dh_dl - phiu*dh_dt;
    gdouble g = GFS_SIMULATION (r)->physical_params.g;

    gdouble eta = GFS_VALUE (cell, r->v1[H]);
    GFS_VALUE (cell, r->v[U]) += r->dt*(g*eta*eta/2.*dh_dl + fG*phiv)/dldh;
    GFS_VALUE (cell, r->v[V]) += r->dt*(g*eta*eta/2.*dh_dt - fG*phiu)/dldh;
  }
}

typedef struct {
  double * a; /* sub-diagonal indexed from 1..n-1 */
  double * b; /* diagonal (destroyed) */
  double * c; /* sup-diagonal indexed from 0..n-2 */
  double * v; /* rhs (destroyed) */
  int n;
} Tridiagonal;

static void tridiagonal_init (Tridiagonal * t, int n)
{
  t->a = g_malloc (sizeof (double)*n);
  t->b = g_malloc (sizeof (double)*n);
  t->c = g_malloc (sizeof (double)*(n - 1));
  t->v = g_malloc (sizeof (double)*n);
  t->n = n;
}

static void tridiagonal_solve (Tridiagonal * t, double * x)
{
  int n = t->n;
  double * a = t->a, * b = t->b, * c = t->c, * v = t->v;
  for (int i = 1; i < n; i++) {
    double m = a[i]/b[i-1];
    b[i] -= m*c[i-1];
    v[i] -= m*v[i-1];
  }
  x[n-1] = v[n-1]/b[n-1];  
  for (int i = n - 2; i >= 0; i--)
    x[i] = (v[i] - c[i]*x[i+1])/b[i];  
}

static void tridiagonal_free (Tridiagonal * t)
{
  g_free (t->a);
  g_free (t->b);
  g_free (t->c);
  g_free (t->v);
}

/* see doc/figures/diffusion.tm */
static void vertical_diffusion (double * u, 
				const double * mu, const double * dz,
				int N, double dt,
				double dut, 
				double lambdab, double ub, double k,
				Tridiagonal * t,
				/* working array of size N */
				double * a)
{
  if (k > 0.) { /* use Navier coefficient k rather than slip length */
    lambdab = mu[0]/k;
    ub = 0.;
  }
  for (int l = 0; l < N - 1; l++)
    a[l] = dt*(mu[l] + mu[l+1])/(dz[l]*(dz[l] + dz[l+1]));
  a[N-1] = dt*mu[N-1]/(dz[N-1]*dz[N-1]);
  double am = dt*mu[0]/(dz[0]*dz[0]);
  t->b[0] = 1. + a[0] + (1. - (2.*lambdab - dz[0])/(2.*lambdab + dz[0]))*am;
  t->c[0] = - a[0];
  t->v[0] = u[0] + 2.*dz[0]/(2.*lambdab + dz[0])*ub*am;
  for (int l = 1; l < N - 1; l++) {
    t->a[l] = - a[l-1];
    t->b[l] = 1. + a[l] + a[l-1];
    t->c[l] = - a[l];
    t->v[l] = u[l];
  }
  t->a[N-1] = - a[N-2];
  t->b[N-1] = 1. + a[N-2];
  t->v[N-1] = u[N-1] + dut*dz[N-1]*a[N-1];
  tridiagonal_solve (t, u);
}

/* bottom friction for a single layer. For more than one layer, bottom
   friction is a boundary condition for vertical_diffusion() above  */
static void bottom_friction (FttCell * cell, GfsRiver * r)
{
  double h = GFS_VALUE (cell, r->v[H]);
  if (h > r->dry) {
    double a = 1. + gfs_function_value (r->k, cell)/h*r->dt;
    GFS_VALUE (cell, r->v[U]) /= a;
    GFS_VALUE (cell, r->v[V]) /= a;
  }
  else
    GFS_VALUE (cell, r->v[U]) = GFS_VALUE (cell, r->v[V]) = 0.;
}

static void domain_vertical_diffusion (GfsRiver * r, double dt)
{
  Tridiagonal tri;
  int n = r->nlayers;
  tridiagonal_init (&tri, n);
  double * a = g_malloc (n*sizeof (double));
  double * u = g_malloc (n*sizeof (double));
  double * mu = g_malloc (n*sizeof (double));
  double * dz = g_malloc (n*sizeof (double));

  FttCellTraverse * t = gfs_domain_cell_traverse_new (GFS_DOMAIN (r), 
						      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1);
  FttCell * cell;
  while ((cell = ftt_cell_traverse_next (t))) {
    double h = GFS_VALUE (cell, r->v[H]);
    if (h > r->dry) {
      double nu = gfs_function_value (r->nu, cell);
      for (int l = 0; l < n; l++) {
	mu[l] = nu;
	dz[l] = r->dz[l]*h;
	u[l] = GFS_VALUE (cell, r->v[U + 2*l])/dz[l];
      }
      double lambdab = 0., ub = 0., dut = r->dut ? gfs_function_value (r->dut, cell) : 0.;
      double k = r->k ? gfs_function_value (r->k, cell) : 0.;
      vertical_diffusion (u, mu, dz, n, dt, dut, lambdab, ub, k, &tri, a);
      for (int l = 0; l < n; l++) {
	GFS_VALUE (cell, r->v[U + 2*l]) = u[l]*dz[l];
	u[l] = GFS_VALUE (cell, r->v[V + 2*l])/dz[l];
      }
      dut = 0.;
      vertical_diffusion (u, mu, dz, n, dt, dut, lambdab, ub, k, &tri, a);
      for (int l = 0; l < n; l++)
	GFS_VALUE (cell, r->v[V + 2*l]) = u[l]*dz[l];
    }
    else
      for (int l = 0; l < n; l++) {
	GFS_VALUE (cell, r->v[U + 2*l]) = 0.;
	GFS_VALUE (cell, r->v[V + 2*l]) = 0.;
      }
  }
  ftt_cell_traverse_destroy (t);

  g_free (a);
  g_free (u);
  g_free (mu);
  g_free (dz);
  tridiagonal_free (&tri);
}

static void advance (GfsRiver * r, gdouble dt)
{
  GfsDomain * domain = GFS_DOMAIN (r);
  guint i;

  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) reset_fluxes, r);
  r->dt = dt;
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) face_fluxes, r);
  if (r->nlayers > 1)
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) vertical_advection, r);
  gfs_domain_traverse_mixed (domain,
  			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
  			     (FttCellTraverseFunc) solid_boundary_fluxes, domain);
  if (domain->cell_metric)
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) metric_sources, r);
  for (i = 0; i < r->nvar; i++) {
    GfsAdvectionParams par;
    par.v = r->v[i];
    par.fv = r->flux[i];
    par.average = FALSE;
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, &par);
    gfs_domain_variable_centered_sources (domain, par.v, par.v, dt);
  }
  if (r->nlayers > 1) {
    /* also add "global" sources specified through U,V to momentum on each layer */
    GfsVariable ** u = gfs_domain_velocity (domain);
    for (int l = 0; l < r->nlayers; l++)
      for (FttComponent c = 0; c < 2; c++)
	gfs_domain_variable_centered_sources (domain, u[c], r->v[U + c + 2*l], dt);
    /* and "global" sources on tracers */
    for (int i = 0; i < r->nt; i++)
      for (r->l = 0; r->l < r->nlayers; r->l++)
	gfs_domain_variable_centered_sources (domain, 
					      PARENT_TRACER(r->v[T(i,r->l)]), r->v[T(i,r->l)], 
					      dt);
  }
  gfs_source_coriolis_implicit (domain, dt);
  if (r->nu)
    domain_vertical_diffusion (r, dt);
  else if (r->k) {
    g_assert (r->nlayers == 1);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) bottom_friction, r);
  }
  for (i = 0; i < r->nvar; i++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, r->v[i]);
}

static void copy (FttCell * cell, const GfsRiver * r)
{
  guint v;
  for (v = 0; v < r->nvar; v++)
    GFS_VALUE (cell, r->v1[v]) = GFS_VALUE (cell, r->v[v]);
}

static void cell_H (FttCell * cell, const GfsRiver * r)
{
  GFS_VALUE (cell, r->h) = GFS_VALUE (cell, r->zb) + GFS_VALUE (cell, r->v[H]);
  gdouble u = 0., v = 0.;
  int l;
  for (l = 0; l < r->nlayers; l++) {
    u += GFS_VALUE (cell, r->v[U + 2*l]);
    v += GFS_VALUE (cell, r->v[V + 2*l]);
  }
  GFS_VALUE (cell, r->qx) = u;
  GFS_VALUE (cell, r->qy) = v;
}

static void cell_gradients (FttCell * cell,
			    const GfsRiver * r)
{
  FttComponent c;
  guint v;

  if (GFS_VALUE (cell, r->v[H]) <= r->dry) {
    for (c = 0; c < FTT_DIMENSION; c++) {
      for (v = 0; v < r->nvar; v++)
	GFS_VALUE (cell, r->dv[c][v]) = 0.;
      GFS_VALUE (cell, r->dv[c][ZB]) = 0.;
    }
  }
  else { /* wet */
    for (c = 0; c < FTT_DIMENSION; c++) {
      for (v = 0; v < r->nvar; v++)
	GFS_VALUE (cell, r->dv[c][v]) = (* r->gradient) (cell, c, r->v[v]->i)/2.;
      /* recontruct Zb + eta rather than Zb: see Theorem 3.1 of Audusse et al, 2004 */
      GFS_VALUE (cell, r->dv[c][ZB]) =
	(* r->gradient) (cell, c, r->h->i)/2.
	- GFS_VALUE (cell, r->dv[c][H]);
    }
  }
}

typedef struct { 
  FttCellTraverseFunc func;
  FttDirection d;
  gpointer data;
} FaceTraverseData;

static void face_traverse (FttCell * cell, FaceTraverseData * p)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, p->d);
  if (neighbor)
    (* p->func) (neighbor, p->data);
}

static void domain_traverse_all_leaves (GfsDomain * domain,
					FttCellTraverseFunc func,
					gpointer data)
{
  FaceTraverseData p;

  gfs_domain_traverse_leaves (domain, func, data);
  /* now traverses boundary cells */
  p.func = func;
  p.data = data;
  for (p.d = 0; p.d < FTT_NEIGHBORS; p.d++)
    gfs_domain_cell_traverse_boundary (domain, p.d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				       (FttCellTraverseFunc) face_traverse, &p);
}

static void dirichlet_p (FttCellFace * f, GfsBc * b)
{
  GFS_VALUE (f->cell, b->v) = gfs_function_face_value (GFS_BC_VALUE (b)->val, f);
}

static void fix_box_bc (GfsBox * box, GfsRiver * r)
{
  for (FttDirection d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, r->v[0]);
      if (bc && GFS_IS_BC_DIRICHLET (bc)) {
	/* use first-order Dirichlet BC for P in GfsRiver */
	bc->bc = (FttFaceTraverseFunc) dirichlet_p;
      }
      if (r->nlayers > 1) {
	GfsBc * bc = g_hash_table_lookup (b->bc, r->qx->name);
	if (bc) {
	  /* applies explicit BC on global flux on all layer fluxes which have no explicit BCs */
	  for (int l = 0; l < r->nlayers; l++) {
	    GfsVariable * v = r->v[U + 2*l];
	    if (!g_hash_table_lookup (b->bc, v->name))
	      g_hash_table_insert (b->bc, v->name, bc);
	  }
	}
      }
    }
}

static void river_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRiver * r = GFS_RIVER (sim);

  r->v[ZB] = r->zb = gfs_variable_from_name (domain->variables, "Zb");

  r->g = sim->physical_params.g/sim->physical_params.L;
  r->variable_density = (r->nlayers > 1 && sim->physical_params.alpha != NULL);

  r->gradient = sim->advection_params.gradient;
  if (r->gradient == gfs_center_minmod_gradient)
    r->gradient = center_minmod_gradient;
  else if (r->gradient == gfs_center_superbee_gradient)
    r->gradient = center_superbee_gradient;
  else if (r->gradient == gfs_center_sweby_gradient)
    r->gradient = center_sweby_gradient;

  gts_container_foreach (GTS_CONTAINER (r), (GtsFunc) fix_box_bc, r);
  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  gfs_simulation_set_timestep (sim);

  domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    /* events */
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* update H */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);
    
    /* gradients */
    gfs_domain_timer_start (domain, "gradients");
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) cell_gradients, r);
    FttComponent c;
    guint v;
    for (c = 0; c < FTT_DIMENSION; c++)
      for (v = 0; v < r->nvar + 1; v++)
	gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, r->dv[c][v]);
    gfs_domain_timer_stop (domain, "gradients");

    /* predictor */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) copy, r);
    if (r->time_order == 2) {
      gfs_domain_timer_start (domain, "predictor");
      for (v = 0; v < r->nvar; v++)
	gfs_variables_swap (r->v[v], r->v1[v]);
      advance (r, sim->advection_params.dt/2.);
      for (v = 0; v < r->nvar; v++)
	gfs_variables_swap (r->v[v], r->v1[v]);
      gfs_domain_timer_stop (domain, "predictor");
    }
    /* corrector */
    gfs_domain_timer_start (domain, "corrector");
    advance (r, sim->advection_params.dt);
    gfs_domain_timer_stop (domain, "corrector");

    /* update H */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
}

static gdouble maximum_face_metric (FttCell * cell, GfsDomain * domain, FttComponent c)
{
  if (domain->face_metric) {
    FttCellFace f;
    f.cell = cell; f.d = 2*c;
    gdouble fm1 = (* domain->face_metric) (domain, &f);
    f.d = 2*c + 1;
    gdouble fm2 = (* domain->face_metric) (domain, &f);
    return MAX (fm1, fm2);
  }
  else
    return 1.;
}

static void minimum_cfl (FttCell * cell, GfsRiver * r)
{
  gdouble h = GFS_VALUE (cell, r->v[H]);
  if (h > r->dry) {
    GfsDomain * domain = GFS_DOMAIN (r);
    gdouble vol = ftt_cell_size (cell);
    if (domain->cell_metric)
      vol *= (* domain->cell_metric) (domain, cell);
    gdouble cg = sqrt (r->g*h);
    FttComponent c;
    for (c = FTT_X; c <= FTT_Y; c++) {
      gdouble fm = maximum_face_metric (cell, domain, c);
      int l;
      for (l = 0; l < r->nlayers; l++) {
	gdouble uh = fabs (GFS_VALUE (cell, r->v[c + 1 + 2*l]));
	gdouble cfl = vol/(fm*(uh/(r->dz[l]*h) + cg));
	if (cfl < r->cfl)
	  r->cfl = cfl;
      }
    }
  }
}

static gdouble river_cfl (GfsSimulation * sim)
{
  GfsRiver * r = GFS_RIVER (sim);
  r->cfl = G_MAXDOUBLE;
  gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) minimum_cfl, r);
  gfs_all_reduce (GFS_DOMAIN (sim), r->cfl, MPI_DOUBLE, MPI_MIN);
  return r->cfl;
}

static void (* default_tracer_read) (GtsObject ** o, GtsFile * fp) = NULL;

static void river_tracer_read (GtsObject ** o, GtsFile * fp)
{
  (* default_tracer_read) (o, fp);

  GfsVariable * v = GFS_VARIABLE (*o);
  GfsDomain * domain = v->domain;
  GfsRiver * r = GFS_RIVER (domain);
  if (fp->type == GTS_ERROR)
    return;

  int oldn = r->nvar;
  r->nvar += r->nlayers;
  r->v = g_realloc (r->v, (r->nvar + 1)*sizeof (GfsVariable *));
  r->v1 = g_realloc (r->v1, r->nvar*sizeof (GfsVariable *));
  for (FttComponent c = 0; c < FTT_DIMENSION; c++)
    r->dv[c] = g_realloc (r->dv[c], (r->nvar + 1)*sizeof (GfsVariable *));
  r->dv[0][ZB] = r->dv[0][oldn];
  r->dv[1][ZB] = r->dv[1][oldn];
  r->flux = g_realloc (r->flux, r->nvar*sizeof (GfsVariable *));
  r->uL = g_realloc (r->uL, (r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));
  r->uR = g_realloc (r->uR, (r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));

  v->units += 1.; /* depth-integrated concentration */
  if (r->nlayers > 1)
    for (int l = 0; l < r->nlayers; l++) {
      gchar * name = g_strdup_printf ("%s%d", v->name, l);
      gchar * description = g_strdup_printf ("%s for layer %d", v->description, l);
      r->v[T(r->nt,l)] = gfs_domain_get_or_add_variable (domain, name, description);
      r->v[T(r->nt,l)]->units = v->units;
      PARENT_TRACER (r->v[T(r->nt,l)]) = v;
      g_free (name);
      g_free (description);
    }
  else
    r->v[T(r->nt,0)] = v;

  for (int l = 0; l < r->nlayers; l++) {
    int i = T(r->nt,l);
    r->v1[i] = gfs_domain_add_variable (domain, NULL, NULL);
    GfsVariable * v[2];
    r->dv[0][i] = v[0] = gfs_domain_add_variable (domain, NULL, NULL);
    r->dv[1][i] = v[1] = gfs_domain_add_variable (domain, NULL, NULL);
    gfs_variable_set_vector (v, 2);
    r->flux[i] = gfs_domain_add_variable (domain, NULL, NULL);
  }
  r->nt++;
}

static void river_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsRiver * river = GFS_RIVER (*o);
  if (fp->type == '{') {
    double dry;
    gchar * scheme = NULL;
    if (!river->nu) {
      river->nu = gfs_function_new (gfs_function_class (), 1.);
      gfs_function_set_units (river->nu, 2.);
      gfs_object_simulation_set (river->nu, river);
    }
    if (!river->dut) {
      river->dut = gfs_function_new (gfs_function_class (), 0.);
      gfs_object_simulation_set (river->dut, river);
    }
    if (!river->k) {
      river->k = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_set_units (river->k, 1.);
      gfs_object_simulation_set (river->k, river);
    }
    GtsFileVariable var[] = {
      {GTS_UINT,   "time_order", TRUE, &river->time_order},
      {GTS_DOUBLE, "dry",        TRUE, &dry},
      {GTS_STRING, "scheme",     TRUE, &scheme},
      {GTS_OBJ,    "nu",         TRUE, &river->nu},
      {GTS_OBJ,    "dut",        TRUE, &river->dut},
      {GTS_OBJ,    "k",          TRUE, &river->k},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[1].set)
      river->dry = dry/GFS_SIMULATION (river)->physical_params.L;
    if (!var[3].set || river->nlayers < 2) {
      gts_object_destroy (GTS_OBJECT (river->nu));
      river->nu = NULL;
    }
    if (!var[4].set || river->nlayers < 2) {
      gts_object_destroy (GTS_OBJECT (river->dut));
      river->dut = NULL;
    }
    if (var[5].set) {
      if (river->nlayers > 1 && !river->nu) {
	gts_file_variable_error (fp, var, "k", "Navier condition requires viscosity to be set");
	return;
      }
    }
    else {
      gts_object_destroy (GTS_OBJECT (river->k));
      river->k = NULL;
    }
    if (scheme) {
      if (!strcmp (scheme, "hllc")) {
	if (river->nlayers > 1)
	  gts_file_error (fp, "HLLC solver can only be used for a single layer");
	else
	  river->scheme = riemann_hllc;
      }
      else if (!strcmp (scheme, "kinetic"))
	river->scheme = riemann_kinetic;
      else
	gts_file_error (fp, "unknown scheme '%s'", scheme);
      g_free (scheme);
    }
  }

  GfsSourceCoriolis * s = gfs_has_source_coriolis (GFS_DOMAIN (river));
  if (s)
    s->beta = 1.; /* backward Euler */
}

static void river_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->write) (o, fp);

  GfsRiver * river = GFS_RIVER (o);
  fprintf (fp, " {\n"
	   "  time_order = %d\n"
	   "  dry = %g\n"
	   "  scheme = %s\n",
	   river->time_order,
	   river->dry*GFS_SIMULATION (river)->physical_params.L,
	   river->scheme == riemann_hllc ? "hllc" : "kinetic");
  if (river->nu) {
    fputs ("  nu =", fp);
    gfs_function_write (river->nu, fp);
  }
  if (river->dut) {
    fputs ("  dut =", fp);
    gfs_function_write (river->dut, fp);
  }
  if (river->k) {
    fputs ("  k =", fp);
    gfs_function_write (river->k, fp);
  }
  fputs ("\n}", fp);
}

static void river_destroy (GtsObject * o)
{
  GfsRiver * r = GFS_RIVER (o);

  g_free (r->v);
  g_free (r->v1);
  g_free (r->flux);
  g_free (r->massflux);
  g_free (r->uL);
  g_free (r->uR);
  g_free (r->f);
  int i;
  for (i = 0; i < FTT_DIMENSION; i++)
    g_free (r->dv[i]);
  if (r->nu)
    gts_object_destroy (GTS_OBJECT (r->nu));
  if (r->dut)
    gts_object_destroy (GTS_OBJECT (r->dut));
  if (r->k)
    gts_object_destroy (GTS_OBJECT (r->k));

  /* restore the default read() method for tracers */
  GTS_OBJECT_CLASS (gfs_variable_tracer_class ())->read = default_tracer_read;

  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->destroy) (o);
}

static void river_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = river_destroy;
  GTS_OBJECT_CLASS (klass)->read = river_read;
  GTS_OBJECT_CLASS (klass)->write = river_write;
  klass->run = river_run;
  klass->cfl = river_cfl;
}

static gdouble cell_velocity (FttCell * cell, FttCellFace * face, GfsRiver * r)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble D = GFS_VALUE (cell, r->v[H]);
  gdouble L = GFS_SIMULATION (r)->physical_params.L;
  return D > r->dry ? L*gfs_vector_norm (cell, gfs_domain_velocity (GFS_DOMAIN (r)))/D : 0.;
}

static gdouble cell_velocity2 (FttCell * cell, FttCellFace * face, GfsRiver * r)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble D = GFS_VALUE (cell, r->v[H]);
  gdouble L = GFS_SIMULATION (r)->physical_params.L;
  return D > r->dry ? 
    L*L*gfs_vector_norm2 (cell, gfs_domain_velocity (GFS_DOMAIN (r)))/(D*D) : 0.;
}

static void momentum_coarse_fine (FttCell * parent, GfsVariable * v)
{
  /* Only initializes momentum when @parent is "deep enough". This
     assumes that shallow parents have just been submerged. For these
     shallow cells, childrens' initial momentum defaults to zero. This
     prevents creating spurious large velocities. */
  GfsRiver * r = GFS_RIVER (v->domain);
  if (GFS_VALUE (parent, r->v[H]) > 2.*r->dry)
    gfs_cell_coarse_fine (parent, v);
}

static GfsVariable * massflux (GfsDomain * domain, int l)
{
  gchar * name = g_strdup_printf ("G%d", l);
  gchar * description = g_strdup_printf ("Mass flux between layer %d and %d", l + 1, l);
  GfsVariable * v = gfs_domain_get_or_add_variable (domain, name, description);
  g_free (name);
  g_free (description);
  return v;
}

static void allocate_river (GfsRiver * r, int start, int nl)
{
  r->nlayers = nl;
  r->dz = g_realloc (r->dz, r->nlayers*sizeof (gdouble));
  int l;
  for (l = 0; l < r->nlayers; l++)
    r->dz[l] = 1./r->nlayers;
  r->nvar = 2*r->nlayers + 1;
  r->v = g_realloc (r->v, (r->nvar + 1)*sizeof (GfsVariable *));
  r->v1 = g_realloc (r->v1, r->nvar*sizeof (GfsVariable *));
  for (FttComponent c = 0; c < FTT_DIMENSION; c++)
    r->dv[c] = g_realloc (r->dv[c], (r->nvar + 1)*sizeof (GfsVariable *));
  r->dv[0][ZB] = r->dv[0][3];
  r->dv[1][ZB] = r->dv[1][3];
  r->flux = g_realloc (r->flux, r->nvar*sizeof (GfsVariable *));
  r->massflux = g_realloc (r->massflux, (r->nlayers - 1)*sizeof (GfsVariable *));

  r->uL = g_realloc (r->uL, (r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));
  r->uR = g_realloc (r->uR, (r->nvar + 3*(r->nlayers + 1))*sizeof (gdouble));
  r->f = g_realloc (r->f, (3*r->nlayers + 1)*sizeof (gdouble));

  GfsDomain * domain = GFS_DOMAIN (r);
  if (r->nlayers > 1)
    r->massflux[0] = massflux (domain, 0);

  for (l = start; l < r->nlayers; l++) {
    r->flux[U + 2*l] = gfs_domain_add_variable (domain, NULL, NULL);
    r->flux[V + 2*l] = gfs_domain_add_variable (domain, NULL, NULL);

    if (l < r->nlayers - 1)
      r->massflux[l] = massflux (domain, l);

    r->v1[U + 2*l] = gfs_domain_add_variable (domain, NULL, NULL);
    r->v1[V + 2*l] = gfs_domain_add_variable (domain, NULL, NULL);
    gfs_variable_set_vector (&r->v1[U + 2*l], 2);

    GfsVariable * tensor[2][2];
    r->dv[0][U + 2*l] = tensor[0][0] = gfs_domain_add_variable (domain, NULL, NULL);
    r->dv[1][U + 2*l] = tensor[0][1] = gfs_domain_add_variable (domain, NULL, NULL);
    r->dv[0][V + 2*l] = tensor[1][0] = gfs_domain_add_variable (domain, NULL, NULL);
    r->dv[1][V + 2*l] = tensor[1][1] = gfs_domain_add_variable (domain, NULL, NULL);
    gfs_variable_set_tensor (tensor);
  }
}

static void river_init (GfsRiver * r)
{
  GfsDomain * domain = GFS_DOMAIN (r);

  gts_object_destroy (GTS_OBJECT (gfs_variable_from_name (domain->variables, "Pmac")));

  allocate_river (r, 0, 1); /* one layer by default */

  r->v[H] = gfs_variable_from_name (domain->variables, "P");
  r->v[H]->units = 1.;
  g_free (r->v[H]->description);
  r->v[H]->description = g_strdup ("Fluid depth");
  r->v1[H] = gfs_domain_add_variable (domain, NULL, NULL);
  r->flux[H] = gfs_domain_add_variable (domain, NULL, NULL);

  r->zb = gfs_domain_add_variable (domain, "Zb", "Bed elevation above datum");
  r->zb->units = 1.;

  r->h = gfs_domain_add_variable (domain, "H", "Elevation above datum (Zb + P)");
  r->h->units = 1.;

  GfsVariable * u[2];
  r->dv[0][H] = u[0] = gfs_domain_add_variable (domain, "Px", "x-component of the depth gradient");
  r->dv[1][H] = u[1] = gfs_domain_add_variable (domain, "Py", "y-component of the depth gradient");
  gfs_variable_set_vector (u, 2);

  r->dv[0][ZB] = u[0] = gfs_domain_add_variable (domain, "Zbx", "x-component of the bed slope");
  r->dv[1][ZB] = u[1] = gfs_domain_add_variable (domain, "Zby", "y-component of the bed slope");
  gfs_variable_set_vector (u, 2);

  GfsVariable * v;
  r->v[U] = r->qx = v = gfs_variable_from_name (domain->variables, "U");
  v->units = 2.;
  v->face_source = FALSE;
  g_free (v->description);
  v->description = g_strdup ("x-component of the (depth-integrated) fluid flux");
  v->coarse_fine = momentum_coarse_fine;
  
  r->v[V] = r->qy = v = gfs_variable_from_name (domain->variables, "V");
  v->units = 2.;
  v->face_source = FALSE;
  g_free (v->description);
  v->description = g_strdup ("y-component of the (depth-integrated) fluid flux");
  v->coarse_fine = momentum_coarse_fine;

  GFS_SIMULATION (r)->advection_params.gradient = gfs_center_minmod_gradient;
  GFS_SIMULATION (r)->advection_params.cfl = 0.5;

  GfsDerivedVariable * dv = gfs_derived_variable_from_name (domain->derived_variables, "Velocity");
  dv->func = cell_velocity;
  dv = gfs_derived_variable_from_name (domain->derived_variables, "Velocity2");
  dv->func = cell_velocity2;

  gfs_domain_remove_derived_variable (domain, "Vorticity");
  gfs_domain_remove_derived_variable (domain, "Divergence");
  gfs_domain_remove_derived_variable (domain, "Lambda2");
  gfs_domain_remove_derived_variable (domain, "Curvature");
  gfs_domain_remove_derived_variable (domain, "D2");

  r->time_order = 2;
  r->dry = 1e-6;
  r->scheme = riemann_kinetic;

  /* overload the default read() method for tracers */
  if (!default_tracer_read) {
    GtsObjectClass * klass = GTS_OBJECT_CLASS (gfs_variable_tracer_class ());
    default_tracer_read = klass->read;
    klass->read = river_tracer_read;
  }
}

GfsSimulationClass * gfs_river_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsRiver",
      sizeof (GfsRiver),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) river_class_init,
      (GtsObjectInitFunc) river_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsRiver} */

/**
 * Add multiple layers.  
 * \beginobject{GfsLayers}
 */

static void traverse_layers (GfsDomain * domain, FttCellTraverseFunc func, gpointer data)
{
  GfsRiver * r = GFS_RIVER (domain);
  GfsVariable ** u = gfs_domain_velocity (domain);
  for (r->l = 0; r->l < r->nlayers; r->l++) {
    for (int i = 0; i < r->nt; i++)
      gfs_variables_swap (r->v[T(i,r->l)], PARENT_TRACER (r->v[T(i,r->l)]));
    for (FttComponent c = 0; c < FTT_DIMENSION; c++)
      gfs_variables_swap (r->v[U + c + 2*r->l], u[c]);
    gfs_domain_traverse_leaves (domain, func, data);
    for (FttComponent c = 0; c < FTT_DIMENSION; c++)
      gfs_variables_swap (r->v[U + c + 2*r->l], u[c]);
    for (int i = 0; i < r->nt; i++)
      gfs_variables_swap (r->v[T(i,r->l)], PARENT_TRACER (r->v[T(i,r->l)]));
  }
}

static gdouble cell_sigma (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  GfsRiver * r = GFS_RIVER (sim);
  g_assert (r->l < r->nlayers);
  double sigma = r->dz[r->l]/2.;
  for (int i = 0; i < r->l; i++)
    sigma += r->dz[i];
  return sigma;
}

static gdouble cell_z (FttCell * cell, FttCellFace * face, GfsSimulation * sim)
{
  GfsRiver * r = GFS_RIVER (sim);
  double zb = cell ? GFS_VALUE (cell, r->zb) : gfs_face_interpolated_value (face, r->zb->i);
  double h = cell ? GFS_VALUE (cell, r->v[H]) : gfs_face_interpolated_value (face, r->v[H]->i);
  return (zb + cell_sigma (cell, face, sim)*h)*sim->physical_params.L;
}

static void gfs_layers_read (GtsObject ** o, GtsFile * fp)
{
  gts_file_next_token (fp);
  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (number of layers)");
    return;
  }
  GfsLayers * layers = GFS_LAYERS (*o);
  layers->nl = atoi (fp->token->str);
  if (layers->nl < 1) {
    gts_file_error (fp, "number of layers must be > 0)");
    return;
  }
  gts_file_next_token (fp);

  /* this is specific to GfsRiver for the moment */
  GfsSimulation * sim = gfs_object_simulation (layers);
  if (!GFS_IS_RIVER (sim)) {
    gts_file_error (fp, "layering is only valid for GfsRiver");
    return;
  }
  if (layers->nl < 2)
    return;

  GfsRiver * r = GFS_RIVER (sim);
  GfsDomain * domain = GFS_DOMAIN (r);

  allocate_river (r, 1, layers->nl);

  /* allocate (Ul,Vl) in separate loops so that values are contiguous in memory */
  for (FttComponent c = 0; c < 2; c++)
    for (int l = 0; l < r->nlayers; l++) {
      gchar * name = g_strdup_printf ("%s%d", c ? "V" : "U", l);
      gchar * description = g_strdup_printf ("%s-component of the fluid flux for layer %d",
					     c ? "y" : "x", l);
      GfsVariable * v = gfs_domain_get_or_add_variable (domain, name, description);
      g_free (name);
      g_free (description);
      r->v[U + c + 2*l] = v;
      v->units = 2.;
      v->coarse_fine = momentum_coarse_fine;
    }
  for (int l = 0; l < r->nlayers; l++) {
    GfsVariable * u[2] = { r->v[U + 2*l], r->v[V + 2*l] };
    gfs_variable_set_vector (u, 2);
  }

  /* configure layers traversal */
  domain->traverse_layers = traverse_layers;
  GfsDerivedVariable * z = gfs_derived_variable_from_name (domain->derived_variables, "z");
  z->func = cell_z;
  GfsDerivedVariableInfo info = { "sigma", "vertical coordinate", cell_sigma };
  gfs_domain_add_derived_variable (domain, info);
}

static void gfs_layers_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s %d", o->klass->info.name, GFS_LAYERS (o)->nl);
}

static void gfs_layers_class_init (GtsObjectClass * klass)
{
  GFS_REFINE_CLASS (klass)->refine = NULL;
  klass->read = gfs_layers_read;
  klass->write = gfs_layers_write;
}

static void gfs_layers_init (GfsLayers * l)
{
  l->nl = 1;
}

GtsObjectClass * gfs_layers_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsLayers",
      sizeof (GfsLayers),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_layers_class_init,
      (GtsObjectInitFunc) gfs_layers_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()), &info);
  }

  return klass;
}

/**
 * \endobject{GfsLayers}
 */

/**
 * 
 * \beginobject{GfsBcSubcritical}
 */

static void subcritical (FttCellFace * f, GfsBc * b)
{
  gdouble hb = gfs_function_face_value (GFS_BC_VALUE (b)->val, f);
  GfsRiver * river = GFS_RIVER (b->v->domain);
  gdouble hi = GFS_VALUE (f->neighbor, river->v[H]);

  g_assert (hi >= 0.);
  GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v) + 
    (FTT_FACE_DIRECT (f) ? -1. : 1.)*2.*hi*(sqrt (river->g*hi) - sqrt (river->g*MAX (hb, 0.)));
}

static void bc_subcritical_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_subcritical_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_subcritical_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!GFS_IS_RIVER (bc->v->domain)) {
    gts_file_error (fp, "GfsBcSubcritical only makes sense for GfsRiver simulations");
    return;
  }

  gfs_function_set_units (GFS_BC_VALUE (bc)->val, 1.);
}

static void gfs_bc_subcritical_init (GfsBc * object)
{
  object->bc =  (FttFaceTraverseFunc) subcritical;
}

static void gfs_bc_subcritical_class_init (GtsObjectClass * klass)
{
  klass->read = bc_subcritical_read;
}

GfsBcClass * gfs_bc_subcritical_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_subcritical_info = {
      "GfsBcSubcritical",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_subcritical_class_init,
      (GtsObjectInitFunc) gfs_bc_subcritical_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_subcritical_info);
  }

  return klass;
}

/** \endobject{GfsBcSubcritical} */

/**
 *
 * \beginobject{GfsDischargeElevation}
 */

static void discharge_elevation_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_DISCHARGE_ELEVATION (o)->Q));
  gts_object_destroy (GTS_OBJECT (GFS_DISCHARGE_ELEVATION (o)->profile));

  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->destroy) (o);
}

static void discharge_elevation_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!GFS_IS_RIVER (domain)) {
    gts_file_error (fp, "GfsDischargeElevation only makes sense for GfsRiver simulations");
    return;
  }
  GfsDischargeElevation * bd = GFS_DISCHARGE_ELEVATION (*o);
  gfs_function_read (bd->Q, domain, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '\n')
    gfs_function_read (bd->profile, domain, fp);
  else
    gfs_object_simulation_set (bd->profile, domain);

  bd->P = GFS_RIVER (domain)->v[H];
  g_free (GFS_CONSTANT (bd)->derived->description);
  GFS_CONSTANT (bd)->derived->description = g_strdup ("Elevation for a given discharge");
}

static void discharge_elevation_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_DISCHARGE_ELEVATION (o)->Q, fp);
  if (gfs_function_get_constant_value (GFS_DISCHARGE_ELEVATION (o)->profile) != 0.)
    gfs_function_write (GFS_DISCHARGE_ELEVATION (o)->profile, fp);
}

static void boundary_flux (FttCellFace * f, GfsDischargeElevation * b)
{
  gdouble profile = gfs_function_face_value (b->profile, f);
  if (profile != GFS_NODATA) {
    GfsRiver * river = GFS_RIVER (gfs_object_simulation (b));
    GFS_VALUE (f->cell, river->flux[H]) = 0.;
    gdouble v1 = GFS_VALUE (f->cell, river->v1[H]);
    GFS_VALUE (f->cell, river->v1[H]) = MAX (0.,
					     profile + GFS_CONSTANT (b)->val - 
					     gfs_face_interpolated_value (f, river->zb->i));
    gdouble dt = river->dt;
    river->dt = 1.;
    face_fluxes (f, river);
    river->dt = dt;
    GFS_VALUE (f->cell, river->v1[H]) = v1;
    double h = ftt_cell_size (f->cell);
    b->flow -= GFS_VALUE (f->cell, river->flux[H])*h*h;
  }
}

static void traverse_dirichlet_boundaries (GfsBox * box, GfsDischargeElevation * bd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, bd->P);
      if (GFS_IS_BC_DIRICHLET (bc))
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) boundary_flux, bd);
    }
}

static gboolean discharge_elevation_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS 
			  (gfs_discharge_elevation_class ())->parent_class)->event)
      (event, sim)) {
    GfsConstant * c = GFS_CONSTANT (event);
    GfsDischargeElevation * bd = GFS_DISCHARGE_ELEVATION (event);
    GfsRiver * r = GFS_RIVER (sim);
    guint v;
      
    for (v = 0; v < r->nvar; v++)
      gfs_variables_swap (r->v[v], r->v1[v]);
      
    gfs_catch_floating_point_exceptions ();
    gdouble Q = gfs_function_value (bd->Q, NULL);
    gfs_restore_fpe_for_function (bd->Q);
    gdouble hmax, hmin = 0.;
    gdouble L = sim->physical_params.L;
      
    hmax = c->val*2./L;
    bd->flow = 0.;
    c->val = hmax;
    gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) traverse_dirichlet_boundaries, bd);
    gfs_all_reduce (GFS_DOMAIN (sim), bd->flow, MPI_DOUBLE, MPI_SUM);
    if (Q > bd->flow)
      hmax = 1.;
      
    guint n = 0, nitermin = 4, nitermax = 100;
    c->val = hmax/2.;
    do {
      bd->flow = 0.;
      gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) traverse_dirichlet_boundaries, bd);
      gfs_all_reduce (GFS_DOMAIN (sim), bd->flow, MPI_DOUBLE, MPI_SUM);
      if (bd->flow > Q) {
	hmax = c->val;
	c->val = (c->val + hmin)/2.;
      }
      else {
	hmin = c->val;
	c->val = (c->val + hmax)/2.;
      }
      n++;
    } while (n < nitermax && (n < nitermin || fabs (Q - bd->flow)/Q > bd->tolerance));
    if (n == nitermax)
      g_warning ("discharge_elevation_event() did not converge after %d iterations: %g", n, 
		 fabs (Q - bd->flow)/Q);
    c->val *= L;
    /* g_message ("### flow: %g H: %g nitermax: %d\n", bd->flow*pow(L,3), c->val, n); */

    for (v = 0; v < r->nvar; v++)
      gfs_variables_swap (r->v[v], r->v1[v]);

    return TRUE;
  }
  return FALSE;
}

static void gfs_discharge_elevation_class_init (GtsObjectClass * klass)
{
  klass->destroy  = discharge_elevation_destroy;
  klass->read     = discharge_elevation_read;
  klass->write    = discharge_elevation_write;
  GFS_EVENT_CLASS (klass)->event = discharge_elevation_event;
}

static void gfs_discharge_elevation_init (GfsDischargeElevation * b)
{
  GFS_EVENT (b)->start = 0.; /* this is not an "Init" event */
  b->tolerance = 1e-2;
  b->Q = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (b->Q, 3.);
  b->profile = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (b->profile, 1.);
}

GfsEventClass * gfs_discharge_elevation_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_discharge_elevation_info = {
      "GfsDischargeElevation",
      sizeof (GfsDischargeElevation),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_discharge_elevation_class_init,
      (GtsObjectInitFunc) gfs_discharge_elevation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_constant_class ()),
				  &gfs_discharge_elevation_info);
  }

  return klass;
}

/** \endobject{GfsDischargeElevation} */

/**
 * "Pipe" between two locations.
 * \beginobject{GfsSourcePipe}
 */

static gboolean read_position (GtsFile * fp, FttVector * p)
{
  gchar * start[FTT_DIMENSION];
  if (!gfs_read_vector (fp, start))
    return FALSE;
  FttComponent c;
  p->z = 0.;
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&p->x)[c] = atof (start[c]);
    g_free (start[c]);
  }
  return TRUE;
}

static void source_pipe_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_pipe_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSimulation * sim = gfs_object_simulation (*o);
  if (!GFS_IS_RIVER (sim)) {
    gts_file_error (fp, "%s only makes sense for GfsRiver simulations",
		    (*o)->klass->info.name);
    return;
  }
  GfsVariable * v = GFS_RIVER (sim)->v[H];
  if (v->sources == NULL)
    v->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (v->sources, GTS_CONTAINEE (*o));

  GfsSourcePipe * p = GFS_SOURCE_PIPE (*o);
  if (!read_position (fp, &p->start))
    return;
  if (!read_position (fp, &p->end))
    return;

  p->diameter = gfs_read_constant (fp, GFS_DOMAIN (gfs_object_simulation (p)));
  if (p->diameter == G_MAXDOUBLE)
    return;
}

static void source_pipe_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_pipe_class ())->parent_class->write) (o, fp);
  GfsSourcePipe * p = GFS_SOURCE_PIPE (o);
  fprintf (fp, " (%f,%f) (%f,%f) %f",
	   p->start.x, p->start.y,
	   p->end.x, p->end.y,
	   p->diameter);
}

#define DQ (1e-4/L3)

static double flow_rate_Q (double z1, double h1, double z2, double h2,
			   double l, double g, GfsSourcePipe * p,
			   double a1, double a2, double Q)
{
  double Q1 = (*p->flow_rate) (z1, h1 - Q/a1, z2, h2 + Q/a2, l, g, p);
  if (Q1 > 0.) Q1 = MIN (Q1, a1*h1);
  if (Q1 < 0.) Q1 = MAX (Q1, - a2*h2);
  return Q1;
}

static gboolean source_pipe_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourcePipe * p = GFS_SOURCE_PIPE (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    FttVector start = p->start, end = p->end;
    gfs_simulation_map (sim, &start);
    gfs_simulation_map (sim, &end);
    /* fixme: this won't work in parallel if the ends of the pipe are on different PEs */
    p->scell = gfs_domain_locate (domain, start, -1, NULL);
    p->ecell = gfs_domain_locate (domain, end, -1, NULL);
    p->Q = 0.;
    if (p->scell && p->ecell && p->scell != p->ecell) {
      gdouble L = sim->physical_params.L, g = sim->physical_params.g;
      GfsVariable * h = GFS_RIVER (sim)->v[H], * zb = GFS_RIVER (sim)->zb;
      gdouble h1 = MAX (L*GFS_VALUE (p->scell, h), 0.), z1 = L*GFS_VALUE (p->scell, zb);
      gdouble h2 = MAX (L*GFS_VALUE (p->ecell, h), 0.), z2 = L*GFS_VALUE (p->ecell, zb);      
      /* fixme: the length below does not take into account metric
	 properly (e.g. won't work for MetricLonLat) */
      gdouble l = L*sqrt ((start.x - end.x)*(start.x - end.x) +
			  (start.y - end.y)*(start.y - end.y));
      gdouble L2 = L*L, L3 = L*L*L;
      gdouble a1 = L2*gfs_cell_volume (p->scell, GFS_DOMAIN (sim))/sim->advection_params.dt;
      gdouble a2 = L2*gfs_cell_volume (p->ecell, GFS_DOMAIN (sim))/sim->advection_params.dt;

      /* secant-bisection root-finding: solves flow_rate(h, l, Q) - Q = 0 for the flow rate Q */
      p->Q = (*p->flow_rate) (z1, h1, z2, h2, l, g, p)/L3;
      gdouble Q1 = p->Q*2.;
      gdouble v1 = flow_rate_Q (z1, h1, z2, h2, l, g, p, a1, a2, Q1*L3)/L3 - Q1;
      gdouble Q2 = 0.;
      gdouble v2 = p->Q;
      if (fabs (v1) > DQ && fabs (v2) > DQ) {
	if (v1 > v2) {
	  gdouble v = v1;
	  v1 = v2; v2 = v;
	  v = Q1;
	  Q1 = Q2; Q2 = v;
	}
	if (v1*v2 >= 0.)
	  g_warning ("source_pipe_event: v1: %g v2: %g", v1*L3, v2*L3);
	else {
	  guint nitermax = 1000;
	  gdouble Qb;
	  p->Q = (v1*Q2 - v2*Q1)/(v1 - v2);
	  do {
	    Qb = p->Q;
	    gdouble v = flow_rate_Q (z1, h1, z2, h2, l, g, p, a1, a2, p->Q*L3)/L3 - p->Q;
	    if (v < 0.) {
	      v1 = v; Q1 = p->Q;
	    }
	    else {
	      v2 = v; Q2 = p->Q;
	    }
	    if (v2 > v1)
	      p->Q = (v1*Q2 - v2*Q1)/(v1 - v2);
	    nitermax--;
	  } while (fabs (p->Q - Qb) > DQ && nitermax);
	  if (nitermax == 0)
	    g_warning ("source_pipe_event: failed to converge! %g %g", 
		       p->Q*L3, fabs (p->Q - Qb)*L3);
	}
      }
    }
    return TRUE;
  }
  return FALSE;
}

static void source_pipe_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_pipe_read;
  GTS_OBJECT_CLASS (klass)->write = source_pipe_write;
  GFS_EVENT_CLASS (klass)->event = source_pipe_event;
}

static gdouble source_pipe_value (GfsSourceGeneric * s, 
				  FttCell * cell, 
				  GfsVariable * v)
{
  GfsSourcePipe * p = GFS_SOURCE_PIPE (s);
  if (cell == p->scell)
    return - p->Q/gfs_cell_volume (cell, v->domain);
  if (cell == p->ecell)
    return   p->Q/gfs_cell_volume (cell, v->domain);
  return 0.;
}

/* This is a simplistic flow rate model for a circular pipe. 
   The pipe is assumed to be always fully submerged. */
static double pipe_flow_rate (double z1, double h1, /* terrain elevation and flow depth at inlet */
			      double z2, double h2, /* terrain elevation and flow depth at outlet */
			      double l,             /* pipe length */
			      double g,             /* acceleration of gravity */
			      GfsSourcePipe * p)
{
  gdouble r = p->diameter/2.;
  gdouble A = M_PI*r*r; /* area */
  gdouble P = 2.*M_PI*r; /* perimeter */
  gdouble Rh = A/P; /* hydraulic radius */
  gdouble S = fabs (z1 + h1 - z2 - h2)/l; /* slope */
  gdouble n = 0.03; /* Gauckler-Manning coefficient */
  /* Gauckler-Manning-Strickler formula for the (signed) flow rate */
  return (z1 + h1 > z2 + h2 ? 1. : -1.)*A/n*pow (Rh, 2./3.)*sqrt (S);
}

static void source_pipe_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = source_pipe_value;
  GFS_SOURCE_PIPE (s)->flow_rate = pipe_flow_rate;
}

GfsSourceGenericClass * gfs_source_pipe_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsSourcePipe",
      sizeof (GfsSourcePipe),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_pipe_class_init,
      (GtsObjectInitFunc) source_pipe_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsSourcePipe} */
