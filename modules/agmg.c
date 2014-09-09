/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2010 Stephane Popinet, CNRS
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

#include "variable.h"
#include "poisson.h"

extern void dagmg_ (const int * n, 
		    double * a, int * ja, int * ia, 
		    double * f, double * x, 
		    const int * ijob, const int * iprint, 
		    const int * nrest, int * iter, const double * tol);

static int ijob = 0, verbose = 0, nrest = -1;

static void solve_poisson_problem_using_agmg (GfsDomain * domain,
					      GfsLinearProblem * lp,
					      GfsMultilevelParams * par)
{
  /* fixme: doesn't work in parallel yet */
  g_assert (lp->istart == 0);

  /* create "compressed sparse row" (CSR) format arrays (with Fortran indexing) */
  int i, j, k, n = lp->rhs->len;
  GArray * a = g_array_new (FALSE, FALSE, sizeof (double));
  GArray * ia = g_array_new (FALSE, FALSE, sizeof (int));
  GArray * ja = g_array_new (FALSE, FALSE, sizeof (int));
  for (i = 0; i < n; i++) {
    GfsStencil * stencil = g_ptr_array_index (lp->LP, i);
    k = a->len + 1;
    g_array_append_val (ia, k);
    /* check that the system is properly ordered */
    g_assert (g_array_index (stencil->id, int, 0) == i);
    /* AGMG needs strictly positive diagonal components */
    g_assert (- g_array_index (stencil->coeff, double, 0) > 0.);
    for (j = 0; j < stencil->id->len; j++) {
      double coeff = - g_array_index (stencil->coeff, double, j);
      g_array_append_val (a, coeff);
      k = g_array_index (stencil->id, int, j) + 1;
      g_array_append_val (ja, k);
    }
  }
  k = a->len + 1;
  g_array_append_val (ia, k);

  /* solve */
  int iter = par->nitermax;
  double tol = MIN (par->tolerance/par->residual.infty, 0.99);
  int iprint = verbose ? 6 : -1;
  dagmg_ (&n, (double *) a->data, (int *) ja->data, (int *) ia->data, 
	  (double *) lp->rhs->data, (double *) lp->lhs->data, 
	  &ijob, &iprint, &nrest, &iter, &tol);
  par->niter = iter;

  /* cleanup */
  g_array_free (a, TRUE);
  g_array_free (ia, TRUE);
  g_array_free (ja, TRUE);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * lhs;
} CopyParams;

static void copy_poisson_solution (FttCell * cell, CopyParams * p)
{
   GFS_VALUE (cell, p->lhs) =
     /* we need to change sign because we changed the sign of the
	matrix to make the diagonal terms positive */
     - g_array_index (p->lp->lhs, gdouble, (int) GFS_VALUE (cell, p->lp->id) - p->lp->istart);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void poisson_solve (GfsDomain * domain,
			   GfsMultilevelParams * par,
			   GfsVariable * lhs,
			   GfsVariable * rhs,
			   GfsVariable * res,
			   GfsVariable * dia,
			   gdouble dt)
{
  /* calculates the initial residual and its norm */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

  if (par->nitermax > 0 && par->residual.infty > 0.) {
    GfsVariable * dp = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dp);
    GfsLinearProblem * lp = gfs_get_poisson_problem (domain, res, dp, dia, -1, lhs);
 
    solve_poisson_problem_using_agmg (domain, lp, par);

    CopyParams p = { lp, dp };
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      (FttCellTraverseFunc) copy_poisson_solution, &p);
    gfs_linear_problem_destroy (lp);

    /* correct on leaf cells */
    gpointer data[2];
    data[0] = lhs;
    data[1] = dp;
    gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) correct, data,
			 lhs, lhs);
    gts_object_destroy (GTS_OBJECT (dp));

    /* compute new residual on leaf cells */
    gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);
  }
}

/* Initialize module */
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);
void          gfs_module_write    (FILE * fp);

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "agmg";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);

  /* initialise the poisson cycle hook */
  sim->approx_projection_params.poisson_solve = poisson_solve;
  sim->projection_params.poisson_solve = poisson_solve;

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "ijob",    TRUE, &ijob},
      {GTS_INT, "verbose", TRUE, &verbose},
      {GTS_INT, "nrest",   TRUE, &nrest},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
  }
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);

  fprintf (fp, " { ijob = %d verbose = %d nrest = %d }", ijob, verbose, nrest);
}
