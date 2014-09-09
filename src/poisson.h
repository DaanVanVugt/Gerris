/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2011 National Institute of Water and Atmospheric Research
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

#ifndef __POISSON_H__
#define __POISSON_H__

#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "domain.h"

typedef struct _GfsMultilevelParams GfsMultilevelParams;
typedef void (* GfsPoissonSolverFunc) (GfsDomain * domain,
				       GfsMultilevelParams * par,
				       GfsVariable * lhs,
				       GfsVariable * rhs,
				       GfsVariable * res,
				       GfsVariable * dia,
				       gdouble dt);
struct _GfsMultilevelParams {
  gdouble tolerance;
  guint nrelax, erelax;
  guint minlevel;
  guint nitermax, nitermin;

  guint dimension;
  guint niter;
  guint depth;
  gboolean weighted, function;
  gdouble beta, omega;
  GfsNorm residual_before, residual;
  GfsPoissonSolverFunc poisson_solve;
};

void                  gfs_multilevel_params_init     (GfsMultilevelParams * par);
void                  gfs_multilevel_params_write    (GfsMultilevelParams * par, 
						      FILE * fp);
void                  gfs_multilevel_params_read     (GfsMultilevelParams * par, 
						      GtsFile * fp);
void                  gfs_multilevel_params_stats_write (GfsMultilevelParams * par,
							 FILE * fp);
void                  gfs_relax                      (GfsDomain * domain,
						      guint d,
						      gint max_depth,
						      gdouble omega,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia);
void                  gfs_residual                   (GfsDomain * domain,
						      guint d,
						      FttTraverseFlags flags,
						      gint max_depth,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);
void                  gfs_poisson_coefficients       (GfsDomain * domain,
						      GfsFunction * alpha,
						      gboolean positive,
						      gboolean centered,
						      gboolean reset);
void                  gfs_poisson_cycle              (GfsDomain * domain,
						      GfsMultilevelParams * p,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);
gdouble               gfs_poisson_compatibility      (GfsDomain * domain, 
						      GfsVariable * lhs, 
						      GfsVariable * rhs,
						      gdouble dt);
void                  gfs_poisson_solve              (GfsDomain * domain,
						      GfsMultilevelParams * par,
						      GfsVariable * lhs,
						      GfsVariable * rhs,
						      GfsVariable * res,
						      GfsVariable * dia,
						      gdouble dt);
void                  gfs_diffusion_coefficients     (GfsDomain * domain,
						      GfsSourceDiffusion * d,
						      gdouble dt,
						      GfsVariable * rhoc,
						      GfsVariable * metric,
						      GfsFunction * alpha,
						      gdouble beta);
void                  gfs_diffusion_rhs              (GfsDomain * domain,
						      GfsVariable * v,
						      GfsVariable * rhs,
						      GfsVariable * rhoc,
						      GfsVariable * metric,
						      gdouble beta);
void                  gfs_diffusion_residual         (GfsDomain * domain,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * rhoc,
						      GfsVariable * metric,
						      GfsVariable * res);
void                  gfs_diffusion_cycle            (GfsDomain * domain,
						      guint levelmin,
						      guint depth,
						      guint nrelax,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * rhoc,
						      GfsVariable * metric,
						      GfsVariable * res);

/* GfsLinearProblem: Object */

struct _GfsLinearProblem {
  GPtrArray * LP;
  GArray * rhs, * lhs;
  GfsVariable * id, * neighbor, * neighborw;
  gint istart;
};

GfsLinearProblem * gfs_linear_problem_new            (GfsDomain * domain);
void               gfs_linear_problem_add_stencil    (GfsLinearProblem * lp, 
						      GfsStencil * stencil);
void               gfs_linear_problem_destroy        (GfsLinearProblem * lp);
GfsLinearProblem * gfs_get_poisson_problem           (GfsDomain * domain,
						      GfsVariable * rhs, 
						      GfsVariable * lhs,
						      GfsVariable * dia,
						      gint maxlevel,
						      GfsVariable * v);
GfsLinearProblem * gfs_get_diffusion_problem         (GfsDomain * domain,
						      GfsVariable * rhs, 
						      GfsVariable * lhs,
						      GfsVariable * rhoc,
						      GfsVariable * metric,
						      gint maxlevel,
						      GfsVariable * v);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __POISSON_H__ */
