/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2010 CNRS
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

#include <lis.h>

#include "variable.h"
#include "poisson.h"

static GString * options = NULL;

static void solve_poisson_problem_using_lis (GfsDomain * domain,
					     GfsLinearProblem * lp,
					     GfsMultilevelParams * par)
{
  /* fixme: doesn't work in parallel yet */
  g_assert (lp->istart == 0);

  LIS_MATRIX A;
  lis_matrix_create (LIS_COMM_WORLD, &A);
  lis_matrix_set_size (A, 0, lp->rhs->len);
  int i;
  for (i = 0; i < lp->rhs->len; i++) {
    GfsStencil * stencil = g_ptr_array_index (lp->LP, i);
    int row = g_array_index (stencil->id, int, 0), j;
    for (j = 0; j < stencil->id->len; j++)
      lis_matrix_set_value (LIS_INS_VALUE, 
			    row, 
			    g_array_index (stencil->id, int, j),
			    g_array_index (stencil->coeff, double, j), 
			    A);
  }
  lis_matrix_set_type (A, LIS_MATRIX_CRS);
  lis_matrix_assemble (A);

  LIS_VECTOR b, x;
  lis_vector_duplicate (A, &b);
  lis_vector_duplicate (A, &x);
  for (i = 0; i < lp->rhs->len; i++) {
    lis_vector_set_value (LIS_INS_VALUE, i, g_array_index (lp->rhs, double, i), b);
    lis_vector_set_value (LIS_INS_VALUE, i, g_array_index (lp->lhs, double, i), x);
  }

  LIS_SOLVER solver;
  lis_solver_create (&solver);

  gchar * opt = g_strdup_printf ("%s-maxiter %d -tol %g", 
				 options->str, par->nitermax, 
				 MIN (par->tolerance/par->residual.infty, 0.99));
  lis_solver_set_option (opt, solver);
  g_free (opt);

  lis_solve (A, b, x, solver);
  int iter;
  lis_solver_get_iters (solver, &iter);
  par->niter = iter;

  lis_vector_get_values (x, 0, lp->lhs->len, (double *) lp->lhs->data);

  lis_solver_destroy (solver);
  lis_matrix_destroy (A);
  lis_vector_destroy (b);
  lis_vector_destroy (x);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * lhs;
} CopyParams;

static void copy_poisson_solution (FttCell * cell, CopyParams * p)
{
   GFS_VALUE (cell, p->lhs) =
     g_array_index (p->lp->lhs, gdouble, (int) GFS_VALUE (cell, p->lp->id) - p->lp->istart);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void lis_poisson_solve (GfsDomain * domain,
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

  if (par->nitermax > 0) {
    GfsVariable * dp = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dp);
    GfsLinearProblem * lp = gfs_get_poisson_problem (domain, res, dp, dia, -1, lhs);
 
    solve_poisson_problem_using_lis (domain, lp, par);

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
const gchar gfs_module_name[] = "lis";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);

  /* initialise lis */
  int argc = 0;
  char ** argv = NULL;
  lis_initialize (&argc, &argv);
  
  /* initialise the poisson cycle hook */
  sim->approx_projection_params.poisson_solve = lis_poisson_solve;
  sim->projection_params.poisson_solve = lis_poisson_solve;
  
  if (fp->type != GTS_STRING)
    options = g_string_new ("-i bicgstab");
  else {
    options = g_string_new ("");
    while (fp->type == GTS_STRING) {
      g_string_append (options, fp->token->str);
      g_string_append_c (options, ' ');
      gts_file_next_token (fp);
    }
  }
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);

  fprintf (fp, " %s", options->str);
}
