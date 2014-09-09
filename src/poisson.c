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
/*! \file
 * \brief Poisson and diffusion solvers.
 */

#include <math.h>
#include <stdlib.h>
#include "poisson.h"
#include "solid.h"
#include "source.h"
#include "tension.h"
#include "init.h"
#include "mpi_boundary.h"

/**
 * gfs_multilevel_params_write:
 * @par: the multilevel parameters.
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of the multilevel parameters
 * @par.  
 */
void gfs_multilevel_params_write (GfsMultilevelParams * par, FILE * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp,
           "{\n"
	   "  tolerance = %g\n"
	   "  nrelax    = %u\n"
           "  erelax    = %u\n"
	   "  minlevel  = %u\n"
	   "  nitermax  = %u\n"
	   "  nitermin  = %u\n"
	   "  weighted  = %d\n"
	   "  beta      = %g\n",
	   par->tolerance,
	   par->nrelax,
	   par->erelax,
	   par->minlevel,
	   par->nitermax,
	   par->nitermin,
	   par->weighted,
	   par->beta);
  if (par->omega != 1.)
    fprintf (fp, "  omega     = %g\n", par->omega);
  if (par->function)
    fputs ("  function  = 1\n", fp);
  fputc ('}', fp);
}

void gfs_multilevel_params_init (GfsMultilevelParams * par)
{
  g_return_if_fail (par != NULL);

  par->tolerance = 1e-3;
  par->nrelax    = 4;
  par->erelax    = 1;
  par->minlevel  = 0;
  par->nitermax  = 100;
  par->nitermin  = 1;

  par->dimension = FTT_DIMENSION;
  par->weighted = FALSE;
  par->beta = 1.;
  par->omega = 1.;

  par->function = FALSE;

  par->poisson_solve = gfs_poisson_solve;
}

void gfs_multilevel_params_read (GfsMultilevelParams * par, GtsFile * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "tolerance", TRUE, &par->tolerance},
    {GTS_UINT,   "nrelax",    TRUE, &par->nrelax},
    {GTS_UINT,   "erelax",    TRUE, &par->erelax},
    {GTS_UINT,   "minlevel",  TRUE, &par->minlevel},
    {GTS_UINT,   "nitermax",  TRUE, &par->nitermax},
    {GTS_UINT,   "nitermin",  TRUE, &par->nitermin},
    {GTS_INT,    "weighted",  TRUE, &par->weighted},
    {GTS_DOUBLE, "beta",      TRUE, &par->beta},
    {GTS_DOUBLE, "omega",     TRUE, &par->omega},
    {GTS_INT,    "function",  TRUE, &par->function},
    {GTS_NONE}
  };

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (par->tolerance <= 0.) {
    gts_file_variable_error (fp, var, "tolerance",
			     "tolerance `%g' must be strictly positive",
			     par->tolerance);
    return;
  }
  if (par->nrelax == 0)
    gts_file_variable_error (fp, var, "nrelax", "nrelax must be non zero");
  if (par->erelax == 0)
    gts_file_variable_error (fp, var, "erelax", "erelax must be non zero");
  if (par->beta < 0.5 || par->beta > 1.)
    gts_file_variable_error (fp, var, "beta", "beta must be in [0.5,1]");
}

static gdouble rate (gdouble a, gdouble b, guint n)
{
  if (a > 0. && b > 0. && n > 0)
    return exp (log (b/a)/n);
  return 0.;
}

/**
 * gfs_multilevel_params_stats_write:
 * @par: the multilevel parameters.
 * @fp: a file pointer.
 *
 * Writes in @fp the statistics contained in @p.
 */
void gfs_multilevel_params_stats_write (GfsMultilevelParams * par,
					FILE * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp,
	   "    niter: %4d\n"
	   "    residual.bias:   % 10.3e % 10.3e\n"
	   "    residual.first:  % 10.3e % 10.3e %6.2g\n"
	   "    residual.second: % 10.3e % 10.3e %6.2g\n"
	   "    residual.infty:  % 10.3e % 10.3e %6.2g\n",
	   par->niter,
	   par->residual_before.bias,
	   par->residual.bias,
	   par->residual_before.first,
	   par->residual.first,
	   rate (par->residual.first,
		 par->residual_before.first,
		 par->niter),
	   par->residual_before.second,
	   par->residual.second,
	   rate (par->residual.second,
		 par->residual_before.second,
		 par->niter),
	   par->residual_before.infty,
	   par->residual.infty,
	   rate (par->residual.infty,
		 par->residual_before.infty,
		 par->niter));
}

/* GfsLinearProblem: Object */

/**
 * gfs_linear_problem_new:
 * @domain: a #GfsDomain.
 *
 * Returns: a new #GfsLinearProblem.
 */
GfsLinearProblem * gfs_linear_problem_new (GfsDomain * domain)
{
  g_return_val_if_fail (domain != NULL, NULL);

  GfsLinearProblem * lp = g_malloc (sizeof (GfsLinearProblem));

  lp->rhs = g_array_new (FALSE, FALSE, sizeof (gdouble));
  lp->lhs = g_array_new (FALSE, FALSE, sizeof (gdouble));
  lp->LP = g_ptr_array_new ();
  lp->id = gfs_temporary_variable (domain);
  lp->neighbor = gfs_temporary_variable (domain);
  lp->neighborw = gfs_temporary_variable (domain);
  lp->istart = 0;

  return lp;
}

/**
 * gfs_linear_problem_add_stencil:
 * @lp: a #GfsLinearProblem.
 * @stencil: a #GfsStencil. 
 *
 * Adds a stencil to the linear problem.
 */
void gfs_linear_problem_add_stencil (GfsLinearProblem * lp, 
				     GfsStencil * stencil)
{
  g_return_if_fail (lp != NULL);
  g_return_if_fail (stencil != NULL);

  g_ptr_array_add (lp->LP, stencil);
}

/**
 * gfs_linear_problem_destroy:
 * @lp: a #GfsLinearProblem.
 * 
 * Destroys a #GfsLinearProblem.
 */
void gfs_linear_problem_destroy (GfsLinearProblem * lp)
{
  g_return_if_fail (lp != NULL);

  gts_object_destroy (GTS_OBJECT (lp->id));
  gts_object_destroy (GTS_OBJECT (lp->neighbor));
  gts_object_destroy (GTS_OBJECT (lp->neighborw));

  g_array_free (lp->rhs, TRUE);  
  g_array_free (lp->lhs, TRUE);

  int i;
  for (i = 0; i < lp->LP->len; i++)
    gfs_stencil_destroy (g_ptr_array_index (lp->LP, i));

  g_ptr_array_free (lp->LP, TRUE);
  g_free (lp);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * dia;
  gint maxlevel;
  GfsVariable * metric, * rhs;
} RelaxStencilParams;

static void relax_stencil (FttCell * cell, RelaxStencilParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;
  
  GfsStencil * stencil = gfs_stencil_new (cell, p->lp, 0.);

  g.a = GFS_VALUE (cell, p->dia);
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_stencil (&f, &ng, p->maxlevel, p->lp, stencil);
      g.a += ng.a;
    }
  }

  if (g.a != 0.)
    gfs_stencil_add_element (stencil, cell, p->lp, -g.a);
  else {
    gfs_stencil_destroy (stencil);
    stencil = gfs_stencil_new (cell, p->lp, 1.);
    g_array_index (p->lp->rhs, gdouble, (gint) GFS_VALUE (cell, p->lp->id)) = 0.;
  }

  gfs_linear_problem_add_stencil (p->lp, stencil);
}

static void relax_dirichlet_stencil (FttCell * cell, RelaxStencilParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;
  GfsVariable * id = p->lp->id;
  GfsStencil * stencil = gfs_stencil_new (cell, p->lp, 0.);

  g.a = GFS_VALUE (cell, p->dia);
  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0) {
    g.b = gfs_cell_dirichlet_gradient_flux_stencil (cell, p->maxlevel, 0., p->lp, stencil);
    /* Dirichlet contribution */
    g_array_index (p->lp->lhs, gdouble, (gint) GFS_VALUE (cell, id)) -= g.b;
  }

  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    gfs_face_cm_weighted_gradient_stencil (&f, &ng, p->maxlevel, p->lp, stencil);
    g.a += ng.a;
  }
  if (g.a != 0.)
    gfs_stencil_add_element (stencil, cell, p->lp, -g.a);
  else {
    gfs_stencil_destroy (stencil);
    stencil = gfs_stencil_new (cell, p->lp, 1.);
    g_array_index (p->lp->rhs, gdouble, (gint) GFS_VALUE (cell, id)) = 0.;
  }

  gfs_linear_problem_add_stencil (p->lp, stencil);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * lhs, * rhs;
  guint nleafs;
  gint maxlevel;
  FttFaceTraverseFunc bc_number;
} NumberingParams;

static void leaves_numbering (FttCell * cell, NumberingParams * p)
{ 
  GFS_VALUE (cell, p->lp->id) = p->nleafs++;
  GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, p->lp->neighbor)) = NULL;
  g_array_append_val (p->lp->lhs, GFS_VALUE (cell, p->lhs));
  g_array_append_val (p->lp->rhs, GFS_VALUE (cell, p->rhs));
}

static void bc_number (FttCellFace * f, NumberingParams * p)
{
  /* rhs = 0 as all boundary conditions are homogeneous */
  GFS_VALUE (f->cell, p->rhs) = 0.;
  g_array_append_val (p->lp->lhs, GFS_VALUE (f->cell, p->lhs));
  g_array_append_val (p->lp->rhs, GFS_VALUE (f->cell, p->rhs));
  GFS_VALUE (f->cell, p->lp->id) = p->nleafs++;
}

static void reset_bc (FttCellFace * f, GfsLinearProblem * lp)
{
  GFS_VALUE (f->cell, lp->id) = -1;
  GFS_DOUBLE_TO_POINTER (GFS_VALUE (f->cell, lp->neighbor)) = NULL;
}

static void box_reset_bc (GfsBox * box, GfsLinearProblem * lp)
{ 
  FttDirection d;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]) && !GFS_IS_BOUNDARY_PERIODIC (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      ftt_face_traverse_boundary (b->root, b->d,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttFaceTraverseFunc) reset_bc, lp);
    }
}

#ifdef HAVE_MPI

static void set_mpi_domain_index (GfsDomain * domain, GfsLinearProblem * lp)
{
  int gsize, i;
  guint * mpi_domain_index;

  MPI_Comm_size (MPI_COMM_WORLD, &gsize);
  mpi_domain_index = g_malloc (sizeof (guint)*gsize);

  MPI_Allgather (&lp->rhs->len, 1, MPI_UNSIGNED, mpi_domain_index, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

  lp->istart = 0;
  for (i = 0; i < domain->pid; i++)
    lp->istart += mpi_domain_index[i];
  
  g_free (mpi_domain_index);
}

static void leaves_renumbering (FttCell * cell, NumberingParams * p)
{ 
  GFS_VALUE (cell, p->lp->id) = p->nleafs++;
}

#endif /* HAVE_MPI */

typedef struct {
  GfsVariable * lhs;
  gboolean dirichlet;
} CompatPar;

static void check_box_dirichlet (GfsBox * box, CompatPar * p)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->lhs);
      if (GFS_IS_BC_DIRICHLET (bc)) {
	p->dirichlet = TRUE;
	return;
      }	
    }
}

static void cell_numbering (GfsDomain * domain,
			    GfsLinearProblem * lp,
			    GfsVariable * rhs, GfsVariable * lhs,
			    gint maxlevel)
{
  /* fixme: should it be FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS */
  NumberingParams np = { lp, lhs, rhs, 0, maxlevel, (FttFaceTraverseFunc) bc_number };

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, maxlevel,
			    (FttCellTraverseFunc) leaves_numbering, &np);

#ifdef HAVE_MPI
  /* Renumbering of the different subdomains for parallel simulations */
  if (domain->pid >= 0) {
    set_mpi_domain_index (domain, lp);

    np.nleafs = lp->istart;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, maxlevel,
			      (FttCellTraverseFunc) leaves_renumbering, &np);
  }
#endif /* HAVE_MPI */  

  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, maxlevel, lp->id);
}

/**
 * gfs_get_poisson_problem:
 * @domain: the domain over which the poisson problem is defined
 * @rhs: the variable to use as right-hand side
 * @lhs: the variable to use as left-hand side
 * @dia: the diagonal weight
 * @maxlevel: the maximum level to consider (or -1).
 * @v: a #GfsVariable of which @lhs is an homogeneous version.
 *
 * Extracts the poisson problem associated with @lhs and @rhs.
 *
 * Returns: a #GfsLinearProblem.
 */
GfsLinearProblem * gfs_get_poisson_problem (GfsDomain * domain,
					    GfsVariable * rhs, GfsVariable * lhs,
					    GfsVariable * dia, gint maxlevel,
					    GfsVariable * v)
{
  gfs_domain_timer_start (domain, "get_poisson_problem");

  GfsLinearProblem * lp = gfs_linear_problem_new (domain);
 
  cell_numbering (domain, lp, rhs, lhs, maxlevel);
 
  /* Create stencils on the fly */
  RelaxStencilParams p = { lp, dia, maxlevel };

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_reset_bc, lp);
  gfs_domain_homogeneous_bc_stencil (domain, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
				     maxlevel, lhs, v, lp);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    maxlevel, (FttCellTraverseFunc) (v->centered ? relax_stencil:
							     relax_dirichlet_stencil), &p);

  /* If Neumann conditions are applied everywhere, the solution is
     defined to within a constant so we need to remove one degree of
     freedom to make the solution unique, assuming that the domain is
     simply connected... */
  if (v->centered &&
      /* fixme: this seems to reduce the convergence speed for pure Poisson problems */
      !GFS_IS_POISSON (domain)) {
    /* check whether any boundary has a Dirichlet condition on @v */
    CompatPar p = { v, FALSE };
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) check_box_dirichlet, &p);
    gfs_all_reduce (domain, p.dirichlet, MPI_INT, MPI_MAX);
    if (!p.dirichlet) {
      /* impose P = 0 for one unknown of the system */
      int i;
      for (i = 0; i < lp->rhs->len; i++) {
	GfsStencil * stencil = g_ptr_array_index (lp->LP, i);
	if (g_array_index (stencil->id, int, 0) == 0) {
	  stencil->id->len = 1;
	  stencil->coeff->len = 1;
	  g_array_index (stencil->coeff, double, 0) = -1.;
	  g_array_index (lp->rhs, double, i) = 0.;
	}
	else {
	  int j, ncols = stencil->id->len;
	  for (j = 0; j < ncols; j++)
	    if (g_array_index (stencil->id, int, j) == 0)
	      g_array_index (stencil->coeff, double, j) = 0.;
	}
      }
    }
  }
  
  gfs_domain_timer_stop (domain, "get_poisson_problem");

  return lp;
}

typedef struct {
  guint u, rhs, dia, res;
  gint maxlevel;
  gdouble beta, omega;
  guint metric;
} RelaxParams;

/* relax_stencil() needs to be updated whenever this
 * function is modified
 */
static void relax (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a != 0.)
    GFS_VALUEI (cell, p->u) = (g.b - GFS_VALUEI (cell, p->rhs))/g.a;
  else
    GFS_VALUEI (cell, p->u) = 0.;
}

static void relax2D (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS_2D; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_2D (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a != 0.)
    GFS_VALUEI (cell, p->u) = 
      (1. - p->omega)*GFS_VALUEI (cell, p->u) 
      + p->omega*(g.b - GFS_VALUEI (cell, p->rhs))/g.a;
  else
    GFS_VALUEI (cell, p->u) = 0.;
}

/* relax_dirichlet_stencil() needs to be updated whenever this
   function is modified */
static void relax_dirichlet (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, 0.);
  else
    g.b = 0.;

  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    gfs_face_cm_weighted_gradient (&f, &ng, p->u, p->maxlevel);
    g.a += ng.a;
    g.b += ng.b;
  }
  if (g.a != 0.)
    GFS_VALUEI (cell, p->u) = (g.b - GFS_VALUEI (cell, p->rhs))/g.a;
  else
    GFS_VALUEI (cell, p->u) = 0.;
}

/**
 * gfs_relax:
 * @domain: the domain to relax.
 * @d: number of dimensions (2 or 3).
 * @max_depth: the maximum depth of the domain to relax.
 * @omega: the over-relaxation parameter.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 *
 * Apply one pass of a Jacobi relaxation to all the leaf cells of
 * @domain with a level inferior or equal to @max_depth and to all the
 * cells at level @max_depth. The relaxation should converge (if the
 * right-hand-side @rhs verifies the solvability conditions) toward
 * the solution of a Poisson equation for @u at the maximum depth.
 */
void gfs_relax (GfsDomain * domain,
		guint d,
		gint max_depth,
		gdouble omega,
		GfsVariable * u,
		GfsVariable * rhs,
		GfsVariable * dia)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d > 1 && d <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.maxlevel = max_depth;
  p.omega = omega;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    max_depth,
			    (FttCellTraverseFunc) (u->centered ?
						   (d == 2 ? relax2D : relax) : 
						   relax_dirichlet),
			    &p);
}

static void residual_set (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  GFS_VALUEI (cell, p->res) = GFS_VALUEI (cell, p->rhs) - 
    (g.b - GFS_VALUEI (cell, p->u)*g.a);
}

static void residual_set2D (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS_2D; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_2D (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  GFS_VALUEI (cell, p->res) = GFS_VALUEI (cell, p->rhs) - 
    (g.b - GFS_VALUEI (cell, p->u)*g.a);
}

static void residual_set_dirichlet (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VALUEI (cell, p->dia);
  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, GFS_STATE (cell)->solid->fv);
  else
    g.b = 0.;

  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    gfs_face_cm_weighted_gradient (&f, &ng, p->u, p->maxlevel);
    g.a += ng.a;
    g.b += ng.b;
  }
  GFS_VALUEI (cell, p->res) = GFS_VALUEI (cell, p->rhs) - 
    (g.b - GFS_VALUEI (cell, p->u)*g.a);
}

/**
 * gfs_residual:
 * @domain: a domain.
 * @d: number of dimensions (2 or 3).
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 * @res: the variable to use to store the residual.
 *
 * For each cell of @domain, computes the sum of the residual over
 * the volume of the cell for a Poisson equation with @u as
 * left-hand-side and @rhs as right-hand-side. Stores the result in
 * @res.  
 */
void gfs_residual (GfsDomain * domain,
		   guint d,
		   FttTraverseFlags flags,
		   gint max_depth,
		   GfsVariable * u, GfsVariable * rhs, GfsVariable * dia,
		   GfsVariable * res)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d > 1 && d <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.res = res->i;
  p.maxlevel = max_depth;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth,
			    (FttCellTraverseFunc) (u->centered ? 
						   (d == 2 ? residual_set2D : residual_set) :
						   residual_set_dirichlet),
			    &p);
}

typedef struct {
  gdouble lambda2[FTT_DIMENSION];
  GfsFunction * alpha;
  GfsDomain * domain;
  gboolean positive;
} PoissonCoeff;

static void reset_coeff (FttCell * cell, PoissonCoeff * p)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  if (GFS_IS_MIXED (cell)) {
    FttVector v = {0.,0.,0.};
    GFS_STATE (cell)->solid->v = v;
  }
  for (d = 0; d < FTT_NEIGHBORS; d++)
    f[d].v = 0.;
}

static void poisson_coeff (FttCellFace * face,
			   PoissonCoeff * p)
{
  gdouble alpha = p->alpha ? gfs_function_face_value (p->alpha, face) : 1.;
  gdouble v = p->lambda2[face->d/2]*alpha*gfs_domain_face_fraction (p->domain, face)/
    gfs_domain_face_scale_metric (p->domain, face, face->d/2);

  if (alpha <= 0. && p->positive) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "alpha is negative (%g) at face (%g,%g,%g).\n"
	   "Please check your definition.",
	   alpha, p.x, p.y, p.z);
  }
  GFS_STATE (face->cell)->f[face->d].v += v;
  
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v += v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v +=
      v/FTT_CELLS_DIRECTION (face->d);
    break;
  default:
    g_assert_not_reached ();
  }
}

static void poisson_mixed_coeff (FttCell * cell, PoissonCoeff * p)
{
  if (GFS_IS_MIXED (cell)) {
    gdouble alpha = p->alpha ? gfs_function_value (p->alpha, cell) : 1.;
    if (((cell)->flags & GFS_FLAG_DIRICHLET) == 0) 
      /* Neumann condition (prescribed flux) */
      GFS_STATE (cell)->solid->v.x += alpha;
    else {
      /* Dirichlet */
      GfsSolidVector * s = GFS_STATE (cell)->solid;
      FttVector m = {1.,1.,1.};
      gfs_domain_solid_metric (p->domain, cell, &m);
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	(&s->v.x)[c] += alpha*(&m.x)[c]*(s->s[2*c + 1] - s->s[2*c]);
    }

    if (alpha <= 0. && p->positive) {
      FttVector p;
      ftt_cell_pos (cell, &p);
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "alpha is negative (%g) at cell (%g,%g,%g).\n"
	     "Please check your definition.",
	     alpha, p.x, p.y, p.z);
    }
  }
}

static void face_coeff_from_below (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  guint neighbors = 0;

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
	f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= n;
    /* fixme: this stuff may not be necessary anymore? The 'dumbell'
       test case seems to work fine without this */
    FttCell * neighbor;
    if (f[d].v != 0. && 
	(neighbor = ftt_cell_neighbor (cell, d)) && !GFS_CELL_IS_BOUNDARY (neighbor))
      neighbors++;
  }

  if (neighbors == 1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      f[d].v = 0.;
}

/**
 * gfs_poisson_coefficients:
 * @domain: a #GfsDomain.
 * @alpha: the inverse of density or %NULL.
 * @positive: if %TRUE, @alpha must be strictly positive.
 * @centered: %TRUE if solving for a centered variable.
 * @reset: %TRUE if resetting previous coefficients.
 *
 * Initializes the face coefficients for the Poisson equation
 * \f$\nabla\cdot\alpha\nabla p=\dots\f$.
 *
 * If @alpha is %NULL, it is taken to be unity.
 */
void gfs_poisson_coefficients (GfsDomain * domain,
			       GfsFunction * alpha,
			       gboolean positive,
			       gboolean centered,
			       gboolean reset)
{
  PoissonCoeff p;
  FttComponent i;

  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    p.lambda2[i] = lambda*lambda;
  }
  p.alpha = alpha;
  p.domain = domain;
  p.positive = positive;
  if (reset)
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) reset_coeff, &p);
  if (!centered)
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) poisson_mixed_coeff, &p);
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) poisson_coeff, &p);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, NULL);
}

static void tension_coeff (FttCellFace * face, gpointer * data)
{
  gdouble * lambda2 = data[0];
  GfsSourceTensionGeneric * t = data[1];
  GfsVariable * kappa = GFS_SOURCE_TENSION (data[1])->k;
  gdouble alpha = data[2] ? gfs_function_face_value (data[2], face) : 1.;
  gdouble v = lambda2[face->d/2]*alpha*gfs_domain_face_fraction (kappa->domain, face)*
    gfs_function_face_value (t->sigma, face);
  gdouble k1 = GFS_VALUE (face->cell, kappa);
  gdouble k2 = GFS_VALUE (face->neighbor, kappa);
#if 0
  gdouble c1 = GFS_VALUE (face->cell, t->c);
  gdouble c2 = GFS_VALUE (face->neighbor, t->c);
  gdouble w1 = c1*(1. - c1);
  gdouble w2 = c2*(1. - c2);

  if (w1 + w2 > 0.)
    v *= (w1*k1 + w2*k2)/(w1 + w2);
  else
#endif
  {
    if (k1 < G_MAXDOUBLE) {
      if (k2 < G_MAXDOUBLE)
	v *= (k1 + k2)/2.;
      else
	v *= k1;
    }
    else if (k2 < G_MAXDOUBLE)
      v *= k2;
    else /* the curvature is undefined: we assume that this is because
	    we are far enough from the interface and thus that the
	    surface tension force is zero */
      v = 0.;
  }

  if (alpha <= 0.) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "alpha is negative (%g) at face (%g,%g,%g).\n"
	   "Please check your definition.",
	   alpha, p.x, p.y, p.z);
  }
  GFS_STATE (face->cell)->f[face->d].v = v;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = G_MAXDOUBLE;
    break;
  default:
    g_assert_not_reached ();
  }
}

/**
 * gfs_source_tension_coefficients:
 * @s: a #GfsSourceTension.
 * @domain: a #GfsDomain.
 * @alpha: the inverse of density or %NULL.
 *
 * Initializes the face coefficients with the surface tension term
 * (interface curvature times surface tension coefficient).
 *
 * If @alpha is %NULL, it is taken to be unity.
 */
void gfs_source_tension_coefficients (GfsSourceTension * s,
				      GfsDomain * domain,
				      GfsFunction * alpha)
{
  gdouble lambda2[FTT_DIMENSION];
  gpointer data[3];
  FttComponent i;

  g_return_if_fail (s != NULL);
  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    lambda2[i] = lambda*lambda;
  }
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_coeff, NULL);
  data[0] = lambda2;
  data[1] = s;
  data[2] = alpha;
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) tension_coeff, data);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void get_from_above (FttCell * parent, GfsVariable * v)
{
  guint level = ftt_cell_level (parent);
  FttCellNeighbors n;
  FttCellChildren child;
  FttComponent c;
  FttVector h;
  guint i;

  ftt_cell_neighbors (parent, &n);
  for (c = 0; c < FTT_DIMENSION; c++) {
    FttCellFace f;
    GfsGradient g;
    gdouble g1, g2;
    
    f.cell = parent;
    f.d = 2*c;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v->i, level);
    g1 = g.b - g.a*GFS_VALUE (parent, v);
    f.d = 2*c + 1;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v->i, level);
    g2 = g.b - g.a*GFS_VALUE (parent, v);
    (&h.x)[c] = (g1 - g2)/2.;
  }

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++) 
    if (child.c[i]) {
      FttVector p;
      
      GFS_VALUE (child.c[i], v) = GFS_VALUE (parent, v);
      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VALUE (child.c[i], v) += (&p.x)[c]*(&h.x)[c];
    }
}

static void get_from_below_3D (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0.;
  guint i;
  FttCellChildren child;

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      val += GFS_VALUE (child.c[i], v);
  GFS_VALUE (cell, v) = val/2.;
}

static void get_from_below_2D (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0.;
  guint i;
  FttCellChildren child;

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      val += GFS_VALUE (child.c[i], v);
  GFS_VALUE (cell, v) = val;
}

static void relax_loop (GfsDomain * domain, 
			GfsVariable * dp, GfsVariable * u, 
			RelaxParams * q, guint nrelax,
			FttCellTraverseFunc relaxfunc)
{
  guint n;

  gfs_domain_homogeneous_bc (domain,
			     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel, 
			     dp, u);
  for (n = 0; n < nrelax - 1; n++)
    gfs_traverse_and_homogeneous_bc (domain, FTT_PRE_ORDER, 
				     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel,
				     relaxfunc, q, dp, u);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel,
			    relaxfunc, q);
}

/**
 * gfs_poisson_cycle:
 * @domain: the domain on which to solve the Poisson equation.
 * @p: the #GfsMultilevelParams.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 * @res: the residual.
 *
 * Apply one multigrid iteration to the Poisson equation defined by @u
 * and @rhs.
 *
 * The initial value of @res on the leaves of @root must be set to
 * the residual of the Poisson equation (using gfs_residual()).
 *
 * The face coefficients must be set using gfs_poisson_coefficients().
 *
 * The values of @u on the leaf cells are updated as well as the values
 * of @res (i.e. the cell tree is ready for another iteration).
 */
void gfs_poisson_cycle (GfsDomain * domain,
			GfsMultilevelParams * p,
			GfsVariable * u,
			GfsVariable * rhs,
			GfsVariable * dia,
			GfsVariable * res)
{
  guint l, nrelax, minlevel;
  GfsVariable * dp;
  gpointer data[2];
  
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (p->dimension > 1 && p->dimension <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  dp = gfs_temporary_variable (domain);
  minlevel = MAX (domain->rootlevel, p->minlevel);

  /* compute residual on non-leafs cells */
  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) (p->dimension == 2 ? 
						   get_from_below_2D : 
						   get_from_below_3D),
			    res);

  /* relax top level */
  nrelax = p->nrelax;
  for (l = minlevel; l < p->depth; l++)
    nrelax *= p->erelax;

  RelaxParams q;
  q.u = dp->i;
  q.rhs = res->i;
  q.dia = dia->i;
  q.maxlevel = minlevel;
  q.omega = p->omega;
  
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q.maxlevel,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  FttCellTraverseFunc relaxfunc = (FttCellTraverseFunc)
    (u->centered ? (p->dimension == 2 ? relax2D : relax) : relax_dirichlet);
  relax_loop (domain, dp, u, &q, nrelax, relaxfunc);
  nrelax /= p->erelax;

  /* relax from top to bottom */
  for (q.maxlevel = minlevel + 1; q.maxlevel <= p->depth; q.maxlevel++, nrelax /= p->erelax) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS, 
			      q.maxlevel - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    relax_loop (domain, dp, u, &q, nrelax, relaxfunc);
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) correct, data,
		       u, u);
  /* compute new residual on leaf cells */
  gfs_residual (domain, p->dimension, FTT_TRAVERSE_LEAFS, -1, u, rhs, dia, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

/**
 * gfs_poisson_compatibility:
 * @domain: the domain over which the poisson problem is solved.
 * @lhs: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dt: the timestep.
 *
 * Returns: the absolute value of the compatibility condition for a
 * Poisson equation using @lhs as left-hand-side and @rhs as
 * right-hand-side.
 */
gdouble gfs_poisson_compatibility (GfsDomain * domain,
				   GfsVariable * lhs,
				   GfsVariable * rhs,
				   gdouble dt)
{
  g_return_val_if_fail (domain != NULL, 0.);
  g_return_val_if_fail (lhs != NULL, 0.);
  g_return_val_if_fail (rhs != NULL, 0.);

  /* check whether any boundary has a Dirichlet condition on @lhs */
  CompatPar p = { lhs, FALSE };
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) check_box_dirichlet, &p);
  gfs_all_reduce (domain, p.dirichlet, MPI_INT, MPI_MAX);
  if (p.dirichlet)
    /* Assumes that the domain is simply connected in which case compatibility is guaranteed */
    return 0.;
  /* compute volume integral of right-hand-side (again we assume that
     the domain is simply connected) */
  return fabs (gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, rhs).bias);
}

/**
 * gfs_poisson_solve:
 * @domain: the domain over which the poisson problem is solved.
 * @par: the parameters of the poisson problem.
 * @lhs: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @res: the variable in which to store the residual
 * @dia: the diagonal weight.
 * @dt:  the length of the time-step.
 *
 * Solves the poisson problem over domain using Gerris' native
 * multigrid poisson solver.
 */
void gfs_poisson_solve (GfsDomain * domain, 
			GfsMultilevelParams * par,
			GfsVariable * lhs, GfsVariable * rhs, GfsVariable * res,
			GfsVariable * dia, gdouble dt)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (lhs != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (res != NULL);
  g_return_if_fail (dia != NULL);

  gfs_domain_timer_start (domain, "poisson_solve");

  guint minlevel = par->minlevel;
  par->depth = gfs_domain_depth (domain);
  par->niter = 0;

  /* calculates the initial residual and its norm */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

  gdouble res_max_before = par->residual.infty;

  while (par->niter < par->nitermin ||
	 (par->residual.infty > par->tolerance && par->niter < par->nitermax)) {

    /* Does one iteration */
    gfs_poisson_cycle (domain, par, lhs, rhs, dia, res);
    
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

    if (par->residual.infty == res_max_before) /* convergence has stopped!! */
      break;
    if (par->residual.infty > res_max_before/1.1 && par->minlevel < par->depth)
      par->minlevel++;
    res_max_before = par->residual.infty;
    par->niter++;
  }

  par->minlevel = minlevel;

  gfs_domain_timer_stop (domain, "poisson_solve");
}

typedef struct {
  GfsSourceDiffusion * d;
  gdouble lambda2[FTT_DIMENSION];
  gdouble dt;
  GfsVariable * rhoc, * metric;
  GfsFunction * alpha;
  GfsDomain * domain;
} DiffusionCoeff;

static void diffusion_coef (FttCellFace * face, DiffusionCoeff * c)
{
  gdouble v = 
    c->lambda2[face->d/2]*c->dt*
    gfs_source_diffusion_face (c->d, face)*
    gfs_domain_face_fraction (c->domain, face)/
    gfs_domain_face_scale_metric (c->domain, face, face->d/2);

  GFS_STATE (face->cell)->f[face->d].v = v;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v +=
      v/FTT_CELLS_DIRECTION (face->d);
    break;
  default:
    g_assert_not_reached ();
  }
}

static void diffusion_mixed_coeff (FttCell * cell, DiffusionCoeff * c)
{
  reset_coeff (cell, NULL);
  if (GFS_IS_MIXED (cell)) {
    gdouble diffusion = c->dt*gfs_source_diffusion_cell (c->d, cell);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) == 0)
      /* Neumann condition (prescribed flux) */
      GFS_STATE (cell)->solid->v.x = diffusion;
    else {
      /* Dirichlet */
      GfsSolidVector * s = GFS_STATE (cell)->solid;
      FttVector m = {1.,1.,1.};
      gfs_domain_solid_metric (c->domain, cell, &m);
      FttComponent i;
      for (i = 0; i < FTT_DIMENSION; i++)
	(&s->v.x)[i] = diffusion*(&m.x)[i]*(s->s[2*i + 1] - s->s[2*i]);
    }
  }
  if (c->rhoc) {
    gdouble rho = c->alpha ? 1./gfs_function_value (c->alpha, cell) : 1.;
    if (rho <= 0.) {
      FttVector p;
      ftt_cell_pos (cell, &p);
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "density is negative (%g) at cell (%g,%g,%g).\n"
	     "Please check your definition of alpha.",
	     rho, p.x, p.y, p.z);
    }
    GFS_VALUE (cell, c->rhoc) = rho*gfs_domain_cell_fraction (c->domain, cell);
  }
}

static void viscous_metric_coeff (FttCell * cell, DiffusionCoeff * c)
{
  GFS_VALUE (cell, c->metric) =
    c->dt*
    gfs_source_diffusion_cell (c->d, cell)*
    (* c->domain->viscous_metric_implicit) (c->domain, cell, c->metric->component)
    *gfs_domain_cell_fraction (c->domain, cell)
    /GFS_VALUE (cell, c->rhoc);
}

/**
 * gfs_diffusion_coefficients:
 * @domain: a #GfsDomain.
 * @d: a #GfsSourceDiffusion.
 * @dt: the time-step.
 * @rhoc: where to store the mass.
 * @metric: where to store the implicit metric term (or %NULL).
 * @alpha: the inverse of density or %NULL.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Initializes the face coefficients for the diffusion equation.
 */
void gfs_diffusion_coefficients (GfsDomain * domain,
				 GfsSourceDiffusion * d,
				 gdouble dt,
				 GfsVariable * rhoc,
				 GfsVariable * metric,
				 GfsFunction * alpha,
				 gdouble beta)
{
  DiffusionCoeff coef;
  FttComponent i;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    coef.lambda2[i] = lambda*lambda;
  }
  coef.d = d;
  coef.dt = beta*dt;
  coef.rhoc = rhoc;
  coef.alpha = alpha;
  coef.domain = domain;
  coef.metric = metric;
  gfs_catch_floating_point_exceptions ();
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) diffusion_mixed_coeff, &coef);
  if (coef.metric && coef.rhoc)
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) viscous_metric_coeff, &coef);
  gfs_restore_fpe_for_function (alpha);
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) diffusion_coef, &coef);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, 
			    NULL);
}

static void diffusion_rhs (FttCell * cell, RelaxParams * p)
{
  gdouble f, h, val;
  FttCellNeighbors neighbor;
  FttCellFace face;
  
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      f = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, solid->fv);
    else
      f = solid->fv*solid->v.x;
  }
  else
    f = 0.; /* Neumann condition by default */
  h = ftt_cell_size (cell);
  val = GFS_VALUEI (cell, p->u);
  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient g;

    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient (&face, &g, p->u, -1);
    f += g.b - g.a*val;
  }
  GFS_VALUEI (cell, p->rhs) += p->beta*f/(h*h*GFS_VALUEI (cell, p->dia));
  if (p->metric)
    GFS_VALUEI (cell, p->rhs) -= val*p->beta*GFS_VALUEI (cell, p->metric);
}

/**
 * gfs_diffusion_rhs:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @rhs: a #GfsVariable.
 * @rhoc: the mass.
 * @metric: the metric term.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Adds to the @rhs variable of @cell the right-hand side of the
 * diffusion equation for variable @v.
 *
 * The diffusion coefficients must have been already set using
 * gfs_diffusion_coefficients().
 */
void gfs_diffusion_rhs (GfsDomain * domain, 
			GfsVariable * v, GfsVariable * rhs, 
			GfsVariable * rhoc, GfsVariable * metric,
			gdouble beta)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  p.u = v->i;
  p.rhs = rhs->i;
  p.dia = rhoc->i;
  p.beta = (1. - beta)/beta;
  p.metric = metric ? metric->i : FALSE;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) diffusion_rhs, &p);
}

/* diffusion_relax_stencil() needs to be updated whenever this
   function is modified */
static void diffusion_relax (FttCell * cell, RelaxParams * p)
{
  GfsGradient g = { 0., 0. };
  gdouble h = ftt_cell_size (cell);
  FttCellNeighbors neighbor;
  FttCellFace face;

  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, 0.);

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient (&face, &ng, p->u, p->maxlevel);
    g.a += ng.a;
    g.b += ng.b;
  }
  gdouble a = GFS_VALUEI (cell, p->dia)*h*h;
  g_assert (a != 0.);
  g.a = 1. + g.a/a;
  if (p->metric)
    g.a += GFS_VALUEI (cell, p->metric);
  g_assert (g.a != 0.);
  GFS_VALUEI (cell, p->u) = (g.b/a + GFS_VALUEI (cell, p->res))/g.a;
}

static void diffusion_relax_stencil (FttCell * cell, RelaxStencilParams * p)
{
  GfsGradient g = { 0., 0. };
  gdouble h = ftt_cell_size (cell);
  FttCellNeighbors neighbor;
  FttCellFace face;
  GfsVariable * id = p->lp->id;
  GfsStencil * stencil = gfs_stencil_new (cell, p->lp, 0.);

  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0) {
    g.b = gfs_cell_dirichlet_gradient_flux_stencil (cell, p->maxlevel, 0., p->lp, stencil);
    /* fixme: not sure about this */
    g_array_index (p->lp->lhs, gdouble, (gint) GFS_VALUE (cell, id)) -= g.b;
  }

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;
    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient_stencil (&face, &ng, p->maxlevel, p->lp, stencil);
    g.a += ng.a;
  }
  gdouble a = GFS_VALUE (cell, p->dia)*h*h;
  g_assert (a != 0.);
  g.a = 1. + g.a/a;
  if (p->metric)
    g.a += GFS_VALUE (cell, p->metric);
  g_assert (g.a != 0.);
  gfs_stencil_add_element (stencil, cell, p->lp, - a*g.a);

  gfs_linear_problem_add_stencil (p->lp, stencil);
}

static void diffusion_residual (FttCell * cell, RelaxParams * p)
{
  gdouble a;
  GfsGradient g = { 0., 0. };
  gdouble h;
  FttCellNeighbors neighbor;
  FttCellFace face;

  h = ftt_cell_size (cell);
  a = GFS_VALUEI (cell, p->dia);
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, solid->fv);
    else
      g.b = solid->fv*solid->v.x;
  }

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient (&face, &ng, p->u, -1);
    g.a += ng.a;
    g.b += ng.b;
  }
  a *= h*h;
  g_assert (a != 0.);
  g.a = 1. + g.a/a;
  if (p->metric)
    g.a += GFS_VALUEI (cell, p->metric);
  g.b = GFS_VALUEI (cell, p->rhs) + g.b/a;
  GFS_VALUEI (cell, p->res) = g.b - g.a*GFS_VALUEI (cell, p->u);
}

/**
 * gfs_diffusion_residual:
 * @domain: a #GfsDomain.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @rhoc: the mass.
 * @metric: the metric term.
 * @res: the residual.
 *
 * Sets the @res variable of each leaf cell of @domain to the residual
 * of the diffusion equation for @v.
 *
 * The diffusion coefficients must have been set using
 * gfs_diffusion_coefficients() and the right-hand side using
 * gfs_diffusion_rhs().
 */
void gfs_diffusion_residual (GfsDomain * domain,
			     GfsVariable * u,
			     GfsVariable * rhs,
			     GfsVariable * rhoc,
			     GfsVariable * metric,
			     GfsVariable * res)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (res != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = rhoc->i;
  p.res = res->i;
  p.metric = metric ? metric->i : FALSE;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) diffusion_residual, &p);
}

/**
 * gfs_diffusion_cycle:
 * @domain: the domain on which to solve the diffusion equation.
 * @levelmin: the top level of the multigrid hierarchy.
 * @depth: the total depth of the domain.
 * @nrelax: the number of relaxations to apply at each level.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @rhoc: the mass.
 * @metric: the metric term.
 * @res: the residual.
 *
 * Apply one multigrid iteration to the diffusion equation for @u.
 *
 * The initial value of @res on the leaves of @root must be set to
 * the residual of the diffusion equation using gfs_diffusion_residual().
 *
 * The diffusion coefficients must be set using gfs_diffusion_coefficients().
 *
 * The values of @u on the leaf cells are updated as well as the values
 * of @res (i.e. the cell tree is ready for another iteration).
 */
void gfs_diffusion_cycle (GfsDomain * domain,
			  guint levelmin,
			  guint depth,
			  guint nrelax,
			  GfsVariable * u,
			  GfsVariable * rhs,
			  GfsVariable * rhoc,
			  GfsVariable * metric,
			  GfsVariable * res)
{
  GfsVariable * dp;
  RelaxParams p;
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (res != NULL);

  dp = gfs_temporary_variable (domain);

  /* compute residual on non-leafs cells */
  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_get_from_below_intensive, res);

  /* relax top level */
  p.maxlevel = levelmin;
  p.u = dp->i;
  p.res = res->i;
  p.dia = rhoc->i;
  p.metric = metric ? metric->i : FALSE;

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, levelmin,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  relax_loop (domain, dp, u, &p, 10*nrelax, (FttCellTraverseFunc) diffusion_relax);
  /* relax from top to bottom */
  for (p.maxlevel = levelmin + 1; p.maxlevel <= depth; p.maxlevel++) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS,
			      p.maxlevel - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    relax_loop (domain, dp, u, &p, nrelax, (FttCellTraverseFunc) diffusion_relax);
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) correct, data,
		       u, u);
  /* compute new residual on leaf cells */
  gfs_diffusion_residual (domain, u, rhs, rhoc, metric, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

static void scale_rhs (FttCell * cell, RelaxStencilParams * p)
{
  gdouble h = ftt_cell_size (cell);
  GFS_VALUE (cell, p->rhs) *= - GFS_VALUE (cell, p->dia)*h*h;
}

/**
 * gfs_get_diffusion_problem:
 * @domain: the domain over which the poisson problem is defined
 * @rhs: the variable to use as right-hand side
 * @lhs: the variable to use as left-hand side
 * @rhoc: the mass.
 * @metric: the metric term (or %NULL).
 * @maxlevel: the maximum level to consider (or -1).
 * @v: a #GfsVariable of which @lhs is an homogeneous version.
 *
 * Extracts the diffusion problem associated with @lhs and @rhs.
 *
 * Returns: a #GfsLinearProblem.
 */
GfsLinearProblem * gfs_get_diffusion_problem (GfsDomain * domain,
					      GfsVariable * rhs, 
					      GfsVariable * lhs,
					      GfsVariable * rhoc,
					      GfsVariable * metric,
					      gint maxlevel,
					      GfsVariable * v)
{
  gfs_domain_timer_start (domain, "get_diffusion_problem");

  GfsLinearProblem * lp = gfs_linear_problem_new (domain);
 
  RelaxStencilParams p = {
    lp, rhoc, maxlevel, 
    metric,
    rhs
  };

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    maxlevel, (FttCellTraverseFunc) scale_rhs, &p);
  cell_numbering (domain, lp, rhs, lhs, maxlevel);
 
  /* Creates stencils on the fly */
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_reset_bc, lp);
  gfs_domain_homogeneous_bc_stencil (domain, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
				     maxlevel, lhs, v, lp);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    maxlevel, (FttCellTraverseFunc) diffusion_relax_stencil, &p);
  gfs_domain_timer_stop (domain, "get_diffusion_problem");

  return lp;
}
