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

#include <math.h>
#include <stdlib.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include "variable.h"
#include "poisson.h"
#include "mpi_boundary.h"

/*#define DEBUG*/

typedef enum {
  HYPRE_BOOMER_AMG,
  HYPRE_PCG,
  HYPRE_HYBRID,
  HYPRE_LGMRES,
  HYPRE_BICGSTAB,
  HYPRE_GMRES,
  HYPRE_AMS,
  HYPRE_FLEXGMRES
} HypreSolverType;

typedef enum {
  HYPRE_AMG_PRECOND,
  HYPRE_PARASAILS_PRECOND,
  HYPRE_EUCLID_PRECOND,
  HYPRE_PILUT_PRECOND,
  HYPRE_AMS_PRECOND,
  NO_PRECOND
} HyprePrecondType;

typedef struct _HypreProblem HypreProblem;
typedef struct _HypreSolverParams HypreSolverParams;

struct _HypreSolverParams {
  HypreSolverType solver_type;
  HyprePrecondType precond_type;
  gint relax_type;
  gint coarsening_type;
  gint cycle_type;
  gint nlevel;
  gboolean verbose;
};

/* Parameters for the projection schemes are stored in proj_hp */
HypreSolverParams proj_hp;

struct _HypreProblem {
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
};

static HYPRE_PtrToParSolverFcn HYPRE_precond_solver (void)
{

  if (proj_hp.precond_type == HYPRE_AMG_PRECOND)
    return HYPRE_BoomerAMGSolve;
  else if (proj_hp.precond_type == HYPRE_PARASAILS_PRECOND)
    return HYPRE_ParaSailsSolve;
  else if (proj_hp.precond_type == HYPRE_EUCLID_PRECOND)
    return HYPRE_EuclidSolve;
  else if (proj_hp.precond_type == HYPRE_PILUT_PRECOND)
    return HYPRE_ParCSRPilutSolve;
  else if (proj_hp.precond_type == HYPRE_AMS_PRECOND)
    return HYPRE_AMSSolve;
  else
    g_assert_not_reached ();
  return 0;
}

static HYPRE_PtrToParSolverFcn HYPRE_precond_setup (void)
{
  if (proj_hp.precond_type == HYPRE_AMG_PRECOND)
    return HYPRE_BoomerAMGSetup;
  else if (proj_hp.precond_type == HYPRE_PARASAILS_PRECOND)
    return HYPRE_ParaSailsSetup;
  else if (proj_hp.precond_type == HYPRE_EUCLID_PRECOND)
    return HYPRE_EuclidSetup;
  else if (proj_hp.precond_type == HYPRE_PILUT_PRECOND)
    return HYPRE_ParCSRPilutSetup;
  else if (proj_hp.precond_type == HYPRE_AMS_PRECOND)
    return HYPRE_AMSSetup;
  else
    g_assert_not_reached ();
  return 0;
}

static void ParaSails_precond (HYPRE_Solver * precond)
{
  HYPRE_Solver pc;
  int      sai_max_levels = 1;
  double   sai_threshold = 0.1;
  double   sai_filter = 0.05;
  int      sai_sym = 1;
  
/* Set some parameters (See Reference Manual for more parameters) */
  HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &pc);
  HYPRE_ParaSailsSetParams(pc, sai_threshold, sai_max_levels);
  HYPRE_ParaSailsSetFilter(pc, sai_filter);
  HYPRE_ParaSailsSetSym(pc, sai_sym);
  HYPRE_ParaSailsSetLogging(pc, 3);
  *precond = pc;
}

static void AMG_precond (HYPRE_Solver * precond)
{
  HYPRE_Solver pc;
  
  HYPRE_BoomerAMGCreate(&pc);
  HYPRE_BoomerAMGSetPrintLevel(pc, 1); /* print amg solution info */
  HYPRE_BoomerAMGSetCoarsenType(pc, 6);
  HYPRE_BoomerAMGSetRelaxType(pc, 6); /* Sym G.S./Jacobi hybrid */ 
  HYPRE_BoomerAMGSetNumSweeps(pc, 1);
  HYPRE_BoomerAMGSetTol(pc, 0.0); /* conv. tolerance zero */
  HYPRE_BoomerAMGSetMaxIter(pc, 1); /* do only one iteration! */
  *precond = pc;
}

static void Euclid_precond (HYPRE_Solver * precond)
{
  HYPRE_Solver pc;
    
  HYPRE_EuclidCreate (MPI_COMM_WORLD, &pc);
  *precond = pc;
}

static void Pilut_precond (HYPRE_Solver * precond)
{
  HYPRE_Solver pc;

  HYPRE_ParCSRPilutCreate (MPI_COMM_WORLD, &pc);
  HYPRE_ParCSRPilutSetMaxIter (pc, 3);
  *precond = pc;
}

static void AMS_precond (HYPRE_Solver * precond)
{
  HYPRE_Solver pc;

  HYPRE_AMSCreate (&pc);
  *precond = pc;
}

static void set_precond (HYPRE_Solver * precond)
{
  if (proj_hp.precond_type == HYPRE_PARASAILS_PRECOND)
    ParaSails_precond (precond);
  else if (proj_hp.precond_type == HYPRE_AMG_PRECOND)
    AMG_precond (precond);
  else if (proj_hp.precond_type == HYPRE_EUCLID_PRECOND)
    Euclid_precond (precond);
  else if (proj_hp.precond_type == HYPRE_PILUT_PRECOND)
    Pilut_precond (precond);
  else if (proj_hp.precond_type == HYPRE_AMS_PRECOND)
    AMS_precond (precond);
}

static void destroy_precond (HYPRE_Solver precond)
{
  if (proj_hp.precond_type == HYPRE_PARASAILS_PRECOND)
    HYPRE_ParaSailsDestroy(precond);
  else if (proj_hp.precond_type == HYPRE_AMG_PRECOND)
    HYPRE_BoomerAMGDestroy(precond);
  else if (proj_hp.precond_type == HYPRE_EUCLID_PRECOND)
    HYPRE_EuclidDestroy(precond);
  else if (proj_hp.precond_type == HYPRE_PILUT_PRECOND)
    HYPRE_ParCSRPilutDestroy(precond);
  else if (proj_hp.precond_type == HYPRE_AMS_PRECOND)
    HYPRE_AMSDestroy(precond);
}

/***********************************************/
/*     Boomer Algebraic Multigrid Solver       */
/***********************************************/
static void call_AMG_Boomer_solver (GfsDomain * domain, GfsMultilevelParams * par,
				    HypreProblem * hp)
{
  HYPRE_Solver solver;
  int num_iterations;
  double final_res_norm;

  gfs_domain_timer_start (domain, "Hypre: AMG_Boomer_solver");

  if (proj_hp.nlevel == 0)
    proj_hp.nlevel = gfs_domain_depth (domain);
  
  /* Create solver */
  HYPRE_BoomerAMGCreate(&solver);

  if (proj_hp.verbose)
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_BoomerAMGSetCoarsenType(solver, proj_hp.coarsening_type);
  HYPRE_BoomerAMGSetRelaxType(solver, proj_hp.relax_type);
  HYPRE_BoomerAMGSetCycleType(solver, proj_hp.cycle_type);
  HYPRE_BoomerAMGSetNumSweeps(solver, par->nrelax);     /* Sweeps on each level */
  HYPRE_BoomerAMGSetMaxLevels(solver, proj_hp.nlevel);  /* maximum number of levels */
  HYPRE_BoomerAMGSetTol(solver, par->tolerance);        /* conv. tolerance */
  HYPRE_BoomerAMGSetMaxIter(solver, par->nitermax); /* maximum number of iterations */
  HYPRE_BoomerAMGSetMinIter(solver, par->nitermin); /* minimum number of iterations */

  /* Now setup and solve! */
  HYPRE_BoomerAMGSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_BoomerAMGSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;

  /* Prints informations on the residual */
  if (proj_hp.verbose && domain->pid <= 0) {      
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    printf("\n");
    printf("Iterations = %d\n", num_iterations);
    printf("Final Relative Residual Norm = %e\n", final_res_norm);
    printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_BoomerAMGDestroy(solver);
  gfs_domain_timer_stop (domain, "Hypre: AMG_Boomer_solver");
}

/******************************************/
/*       PreConjugateGradient Solver      */
/******************************************/
static void call_PCG_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: PCG_Solver");

  /* Create solver */
  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRPCGSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRPCGSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRPCGSetTol(solver, par->tolerance); /* conv. tolerance */
  HYPRE_ParCSRPCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  HYPRE_ParCSRPCGSetLogging(solver, 1); /* needed to get run info later */
  

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRPCGSetPrecond(solver,
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
			      precond);
  }

  HYPRE_ParCSRPCGSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRPCGSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRPCGGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRPCGDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: PCG_Solver");
}

/**********************************************/
/*          Hybrid DSCG/PCG Solver            */
/**********************************************/
static void call_Hybrid_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: Hybrid_Solver");

  /* Create solver */
  HYPRE_ParCSRHybridCreate(&solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRHybridSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRHybridSetDSCGMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRHybridSetPCGMaxIter(solver, par->nitermax);
  HYPRE_ParCSRHybridSetTol(solver, par->tolerance); /* conv. tolerance */
  HYPRE_ParCSRHybridSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  HYPRE_ParCSRHybridSetLogging(solver, 1); /* needed to get run info later */

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRHybridSetPrecond(solver,
				 (HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
				 (HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
				 precond);
  }

  HYPRE_ParCSRHybridSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRHybridSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRHybridGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_ParCSRHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRHybridDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: Hybrid_Solver");
}

/******************************************/
/*              BICGSTAB Solver           */
/******************************************/
static void call_BICGSTAB_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: BICGSTAB_Solver");

  /* Create solver */
  HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRBiCGSTABSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRBiCGSTABSetTol(solver, par->tolerance); /* conv. tolerance */
  HYPRE_ParCSRBiCGSTABSetLogging(solver, 1); /* needed to get run info later */

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRBiCGSTABSetPrecond(solver,
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
			      precond);
  }

  HYPRE_ParCSRBiCGSTABSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRBiCGSTABSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRBiCGSTABDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: BICGSTAB_Solver");
}

/******************************************/
/*              LGMRES Solver             */
/******************************************/
static void call_LGMRES_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: LGMRES_Solver");

  /* Create solver */
  HYPRE_ParCSRLGMRESCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRLGMRESSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRLGMRESSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRLGMRESSetTol(solver, par->tolerance); /* conv. tolerance */
  /* HYPRE_ParCSRLGMRESSetTwoNorm(solver, 1);  *//* use the two norm as the stopping criteria */
  HYPRE_ParCSRLGMRESSetLogging(solver, 1); /* needed to get run info later */

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRLGMRESSetPrecond(solver,
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
			      (HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
			      precond);
  }

  HYPRE_ParCSRLGMRESSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRLGMRESSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRLGMRESGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    
    HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRLGMRESDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: LGMRES_Solver");
}

/******************************************/
/*               GMRES Solver             */
/******************************************/
static void call_GMRES_solver (GfsDomain * domain, GfsMultilevelParams * par,
			       HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: GMRES_Solver");

  /* Create solver */
  HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRGMRESSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRGMRESSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRGMRESSetTol(solver, par->tolerance); /* conv. tolerance */
  /* HYPRE_ParCSRGMRESSetTwoNorm(solver, 1);  *//* use the two norm as the stopping criteria */
  HYPRE_ParCSRGMRESSetLogging(solver, 1); /* needed to get run info later */

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRGMRESSetPrecond(solver,
				(HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
				(HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
				precond);
  }

  HYPRE_ParCSRGMRESSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRGMRESSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRGMRESDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: GMRES_Solver");
}

/******************************************/
/*               AMS Solver               */
/******************************************/
static void call_AMS_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: AMS_Solver");

  /* Create solver */
  HYPRE_AMSCreate(&solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_AMSSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_AMSSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_AMSSetTol(solver, par->tolerance); /* conv. tolerance */

  HYPRE_AMSSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_AMSSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_AMSGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_AMSGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_AMSDestroy(solver);
  gfs_domain_timer_stop (domain, "Hypre: AMS_Solver");
}

/******************************************/
/*          FlexGMRES Solver              */
/******************************************/
static void call_FlexGMRES_solver (GfsDomain * domain, GfsMultilevelParams * par,
				   HypreProblem * hp)
{
  HYPRE_Solver solver, precond;
  int num_iterations;
  gfs_domain_timer_start (domain, "Hypre: FlexGMRES_Solver");

  /* Create solver */
  HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_ParCSRFlexGMRESSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_ParCSRFlexGMRESSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_ParCSRFlexGMRESSetTol(solver, par->tolerance); /* conv. tolerance */

  if (proj_hp.precond_type != NO_PRECOND) {
    set_precond (&precond);

    HYPRE_ParCSRFlexGMRESSetPrecond(solver,
				    (HYPRE_PtrToParSolverFcn) HYPRE_precond_solver (),
				    (HYPRE_PtrToParSolverFcn) HYPRE_precond_setup (),
				    precond);
  }

  HYPRE_ParCSRFlexGMRESSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRFlexGMRESSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_ParCSRFlexGMRESGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;
  /*  Run info - needed logging turned on */
  if (proj_hp.verbose && domain->pid <= 0) {
    double final_res_norm;
    HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRFlexGMRESDestroy(solver);
  if (proj_hp.precond_type != NO_PRECOND)
    destroy_precond (precond);
  gfs_domain_timer_stop (domain, "Hypre: FlexGMRES_Solver");
}

static void hypre_problem_new (HypreProblem * hp, GfsDomain * domain,
			       gint size, gint istart)
{
  gfs_domain_timer_start (domain, "HYPRE: Solver setup");
  
  /* Create the matrix.*/
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 
		       istart, istart + size - 1, istart, istart + size - 1, 
		       &hp->A);

  /* Create the vectors rhs and solution.*/
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, istart, istart + size-1, &hp->b);
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, istart, istart + size-1, &hp->x);
  
  /* Choose a parallel csr format storage */
  HYPRE_IJMatrixSetObjectType(hp->A, HYPRE_PARCSR);
  HYPRE_IJVectorSetObjectType(hp->b, HYPRE_PARCSR);
  HYPRE_IJVectorSetObjectType(hp->x, HYPRE_PARCSR);

  /* Initialize before setting coefficients */
  HYPRE_IJMatrixInitialize(hp->A);
  HYPRE_IJVectorInitialize(hp->b);
  HYPRE_IJVectorInitialize(hp->x);
}

static void hypre_problem_destroy (HypreProblem * hp)
{
  g_assert (hp->A);

  HYPRE_IJMatrixDestroy(hp->A);
  HYPRE_IJVectorDestroy(hp->b);
  HYPRE_IJVectorDestroy(hp->x);
}

static void extract_stencil (GfsStencil * stencil, HypreProblem * hp)
{
  int ncols = stencil->id->len;
  int rows = g_array_index (stencil->id, int, 0);
  HYPRE_IJMatrixSetValues (hp->A, 1, &ncols, &rows, 
			   (int *) stencil->id->data, (double *) stencil->coeff->data);
}

static void hypre_problem_init (HypreProblem * hp, GfsLinearProblem * lp,
				GfsDomain * domain)
{
  double *rhs_values, *x_values;
  int    *rows;
  gint i;

  /* Now go through my local rows and set the matrix entries.*/
  rhs_values = (double *) lp->rhs->data;
  x_values = (double *) lp->lhs->data;
  rows = g_malloc (lp->lhs->len*sizeof(int));

  for (i = 0; i < lp->rhs->len; i++) {
    extract_stencil (g_ptr_array_index (lp->LP, i), hp);
    rows[i] = lp->istart + i;
  }

  HYPRE_IJVectorSetValues(hp->b, lp->rhs->len, rows, rhs_values );
  HYPRE_IJVectorSetValues(hp->x, lp->lhs->len, rows, x_values);

  /* Assemble after setting the coefficients */
  HYPRE_IJMatrixAssemble(hp->A);
  HYPRE_IJVectorAssemble(hp->b);
  HYPRE_IJVectorAssemble(hp->x);

  /* Get the parcsr matrix object to use */
  HYPRE_IJMatrixGetObject(hp->A, (void **) &hp->parcsr_A);
  HYPRE_IJVectorGetObject(hp->b, (void **) &hp->par_b);
  HYPRE_IJVectorGetObject(hp->x, (void **) &hp->par_x);

#ifdef DEBUG
  HYPRE_IJMatrixPrint(hp->A, "Aij.dat");
  HYPRE_IJVectorPrint(hp->x, "xi.dat");
  HYPRE_IJVectorPrint(hp->b, "bi.dat");
#endif

  free(rows);
  gfs_domain_timer_stop (domain, "HYPRE: Solver setup");
}

static void hypre_problem_copy (HypreProblem * hp, GfsLinearProblem * lp)
{
  double *x_values;
  int    *rows;
  gint i;

  /* Copy the solution to the GfsLinearProblem structure */
  x_values = g_malloc (lp->lhs->len*sizeof (double));
  rows = g_malloc (lp->lhs->len*sizeof (int));
    
  for (i = 0; i < lp->lhs->len; i++) {
    x_values[i] = 0.;
    rows[i] = i + lp->istart;
  }
    
  HYPRE_IJVectorGetValues(hp->x, lp->lhs->len, rows, x_values);
    
  for (i = 0; i < lp->lhs->len; i++)
    g_array_index (lp->lhs, gdouble, i) = x_values[i];

  free(x_values);
  free(rows);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * lhs;
} CopyParams;

static void copy_solution (FttCell * cell, CopyParams * p)
{
  GFS_VALUE (cell, p->lhs) =
    g_array_index (p->lp->lhs, gdouble, (int) GFS_VALUE (cell, p->lp->id) - p->lp->istart);
}

static void solve_linear_problem (GfsDomain * domain,
				  GfsLinearProblem * lp,
				  GfsMultilevelParams * par,
				  gdouble tolerance)
{
  HypreProblem hp;
  gdouble old = par->tolerance;
  par->tolerance = tolerance;
 
  hypre_problem_new (&hp, domain, lp->rhs->len, lp->istart);
  hypre_problem_init (&hp, lp, domain);
  
  /* Choose a solver and solve the system */
  if (proj_hp.solver_type == HYPRE_BOOMER_AMG)
    call_AMG_Boomer_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_PCG)
    call_PCG_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_HYBRID)
    call_Hybrid_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_LGMRES)
    call_LGMRES_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_GMRES)
    call_GMRES_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_FLEXGMRES)
    call_FlexGMRES_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_AMS)
    call_AMS_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_BICGSTAB)
    call_BICGSTAB_solver (domain, par, &hp);
  else
    g_assert_not_reached();
  
  hypre_problem_copy (&hp, lp);
  hypre_problem_destroy (&hp);
  par->tolerance = old;
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void hypre_poisson_solve (GfsDomain * domain,
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
 
    solve_linear_problem (domain, lp, par, MIN (0.1*par->tolerance/par->residual.infty, 0.99));

    CopyParams p = { lp, dp };
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      (FttCellTraverseFunc) copy_solution, &p);
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

static void hypre_diffusion_solve (GfsDomain * domain,
				   GfsMultilevelParams * par,
				   GfsVariable * lhs,
				   GfsVariable * rhs, 
				   GfsVariable * rhoc,
				   GfsVariable * axi)
{
  /* calculates the initial residual and its norm */
  GfsVariable * res = gfs_temporary_variable (domain);
  gfs_diffusion_residual (domain, lhs, rhs, rhoc, axi, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_variable (domain, res, NULL, FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
  
  if (par->nitermax > 0 && par->residual.infty > 0.) {
    GfsVariable * dp = gfs_temporary_variable (domain);
    dp->component = lhs->component;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dp);
    GfsLinearProblem * lp = gfs_get_diffusion_problem (domain, res, dp, rhoc, axi, -1, lhs);
 
    solve_linear_problem (domain, lp, par, MIN (0.1*par->tolerance/par->residual.infty, 0.99));

    CopyParams p = { lp, dp };
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      (FttCellTraverseFunc) copy_solution, &p);
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
    gfs_diffusion_residual (domain, lhs, rhs, rhoc, axi, res);
    par->residual = gfs_domain_norm_variable (domain, res, NULL, 
					      FTT_TRAVERSE_LEAFS, -1, NULL, NULL);
  }

  gts_object_destroy (GTS_OBJECT (res));
}

static void hypre_solver_write (HypreSolverParams * par,FILE * fp)
{
  fprintf (fp," {\n");
  switch (par->solver_type) {
  case HYPRE_BOOMER_AMG: fputs ("  solver_type      = boomer_amg\n", fp); break;
  case HYPRE_PCG:        fputs ("  solver_type      = pcg\n", fp); break;
  case HYPRE_HYBRID:     fputs ("  solver_type      = hybrid\n", fp); break;
  case HYPRE_LGMRES:     fputs ("  solver_type      = lgmres\n", fp); break;
  case HYPRE_GMRES:      fputs ("  solver_type      = gmres\n", fp); break;
  case HYPRE_FLEXGMRES:  fputs ("  solver_type      = flexgmres\n", fp); break;
  case HYPRE_AMS:        fputs ("  solver_type      = ams\n", fp); break;
  case HYPRE_BICGSTAB:   fputs ("  solver_type      = bicgstab\n", fp); break;
  }

  switch (par->relax_type) {
  case 0: fputs ("  relax_type       = jacobi\n", fp); break;
  case 1: fputs ("  relax_type       = gauss_seidel\n", fp); break;
  case 3: fputs ("  relax_type       = sor-j-forward\n", fp); break;
  case 4: fputs ("  relax_type       = sor-j-backward\n", fp); break;
  case 5: fputs ("  relax_type       = gs-j\n", fp); break;
  case 6: fputs ("  relax_type       = ssor-j\n", fp); break;
  case 7: fputs ("  relax_type       = matvec-jacobi\n", fp); break;
  case 9: fputs ("  relax_type       = direct\n", fp); break;
  }

  switch (par->precond_type) {
  case HYPRE_AMG_PRECOND:      fputs ("  precond_type     = amg\n", fp); break;
  case HYPRE_PARASAILS_PRECOND:fputs ("  precond_type     = parasails\n", fp); break;
  case HYPRE_EUCLID_PRECOND:   fputs ("  precond_type     = euclid\n", fp); break;
  case HYPRE_PILUT_PRECOND:    fputs ("  precond_type     = pilut\n", fp); break;
  case HYPRE_AMS_PRECOND:      fputs ("  precond_type     = ams\n", fp); break;
  case NO_PRECOND:             fputs ("  precond_type     = none\n", fp); break;
  }
    
  switch (par->coarsening_type) {
  case 0: fputs  ("  coarsening_type  = cljp\n", fp); break;
  case 3: fputs  ("  coarsening_type  = ruge_stueben\n", fp); break;
  case 6: fputs  ("  coarsening_type  = falgout\n", fp); break;
  case 8: fputs  ("  coarsening_type  = pmis\n", fp); break;
  case 10: fputs ("  coarsening_type  = hmis\n", fp); break;
  case 21: fputs ("  coarsening_type  = cgc\n", fp); break;
  case 22: fputs ("  coarsening_type  = cgc_e\n", fp); break;
  }

  fprintf (fp,"  cycle_type       = %i\n", par->cycle_type);
    
  fprintf (fp,"  nlevel           = %i\n", par->nlevel);
    
  fprintf (fp,"  verbose          = %i\n", par->verbose);

  fputc ('}', fp);
}

static void hypre_solver_read (HypreSolverParams * par, GtsFile * fp)
{
  gchar * solver_type = NULL, * relax_type = NULL, * coarsening_type = NULL;
  gchar * precond_type = NULL;
  GtsFileVariable var[] = {
    {GTS_STRING,  "relax_type",      TRUE, &relax_type},
    {GTS_STRING,  "solver_type",     TRUE, &solver_type},
    {GTS_STRING,  "precond_type",    TRUE, &precond_type},
    {GTS_STRING,  "coarsening_type", TRUE, &coarsening_type},
    {GTS_INT,     "cycle_type",      TRUE, &par->cycle_type},
    {GTS_INT,     "nlevel",          TRUE, &par->nlevel},
    {GTS_INT,     "verbose",         TRUE, &par->verbose},
    {GTS_NONE}
  };

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (solver_type) {
    if (!strcmp (solver_type, "boomer_amg"))
      par->solver_type = HYPRE_BOOMER_AMG;
    else if (!strcmp (solver_type, "pcg"))
      par->solver_type = HYPRE_PCG;
    else if (!strcmp (solver_type, "hybrid"))
      par->solver_type = HYPRE_HYBRID;
    else if (!strcmp (solver_type, "lgmres"))
      par->solver_type = HYPRE_LGMRES;
    else if (!strcmp (solver_type, "gmres"))
      par->solver_type = HYPRE_GMRES;
    else if (!strcmp (solver_type, "ams"))
      par->solver_type = HYPRE_AMS;
    else if (!strcmp (solver_type, "flexgmres"))
      par->solver_type = HYPRE_FLEXGMRES;
    else if (!strcmp (solver_type, "bicgstab"))
      par->solver_type = HYPRE_BICGSTAB;
    else
      gts_file_variable_error (fp, var, "solver_type", "unknown solver type `%s'", solver_type);
    g_free (solver_type);
  }

  /* The default is no preconditioner */
  if (precond_type) {
    if (!strcmp (precond_type, "amg"))
      par->precond_type = HYPRE_AMG_PRECOND;
    else if (!strcmp (precond_type, "parasails"))
      par->precond_type = HYPRE_PARASAILS_PRECOND;
    else if (!strcmp (precond_type, "euclid"))
      par->precond_type = HYPRE_EUCLID_PRECOND;
    else if (!strcmp (precond_type, "pilut"))
      par->precond_type = HYPRE_PILUT_PRECOND;
    else if (!strcmp (precond_type, "ams"))
      par->precond_type = HYPRE_AMS_PRECOND;
    else if (!strcmp (precond_type, "none"))
      par->precond_type = NO_PRECOND;
    else
      gts_file_variable_error (fp, var, "precond_type", 
			       "unknown preconditioner `%s'", precond_type);
    g_free (precond_type);
  }
  if (fp->type == GTS_ERROR)
    return;

  if (relax_type) {
    if (!strcmp (relax_type, "jacobi"))
      par->relax_type = 0;
    else if (!strcmp (relax_type, "gauss_seidel"))
      par->relax_type = 1;
    else if (!strcmp (relax_type, "sor-j-forward"))
      par->relax_type = 3;
    else if (!strcmp (relax_type, "sor-j-backward"))
      par->relax_type = 4;
    else if (!strcmp (relax_type, "gs-j"))
      par->relax_type = 5;
    else if (!strcmp (relax_type, "ssor-j"))
      par->relax_type = 6;
    else if (!strcmp (relax_type, "matvec-jacobi"))
      par->relax_type = 7;
    else if (!strcmp (relax_type, "direct"))
      par->relax_type = 9;
    else
      gts_file_variable_error (fp, var, "relax_type", "unknown relax type `%s'", relax_type);
    g_free (relax_type);
  }
  if (fp->type == GTS_ERROR)
    return;

  if (coarsening_type) {
    if (par->solver_type != HYPRE_BOOMER_AMG)
      g_warning ("coarsening algorithms are only for the BoomerAMG Solver !!\n"
		 "none will be used with the selected solver");

    if (!strcmp (coarsening_type, "cljp"))
      par->coarsening_type = 0;
    else if (!strcmp (coarsening_type, "ruge_stueben"))
      par->coarsening_type = 3;
    else if (!strcmp (coarsening_type, "falgout"))
      par->coarsening_type = 6;
    else if (!strcmp (coarsening_type, "pmis"))
      par->coarsening_type = 8;
    else if (!strcmp (coarsening_type, "hmis"))
      par->coarsening_type = 10;
    else if (!strcmp (coarsening_type, "cgc"))
      par->coarsening_type = 21;
    else if (!strcmp (coarsening_type, "cgc_e"))
      par->coarsening_type = 22;
    else
      gts_file_variable_error (fp, var, "coarsening_type",
			       "unknown coarsening type `%s'", coarsening_type);
    g_free (coarsening_type);
  }
  if (fp->type == GTS_ERROR)
    return;

  if (par->cycle_type < 1 || (par->cycle_type > 8 && par->cycle_type < 11) ||
      par->cycle_type > 14)
    gts_file_variable_error (fp, var, "cycle_type",
			     "unknown cycle type `%i'", par->cycle_type);
  else if (par->nlevel < 0)
    gts_file_variable_error (fp, var, "nlevel", "nlevel cannot be < 0");
}

/* Initialize module */
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);
void          gfs_module_write    (FILE * fp);

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "hypre";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);

  /* Default hypre solver parameters */
  /* Boomer AMG is the default solver */
  proj_hp.solver_type = HYPRE_BOOMER_AMG;
  proj_hp.precond_type = NO_PRECOND;
  proj_hp.relax_type = 5;
  proj_hp.coarsening_type = 22;
  proj_hp.cycle_type = 1;
  proj_hp.nlevel = 0;
  proj_hp.verbose = FALSE;

  if (fp->type == '{')
    hypre_solver_read (&proj_hp, fp);
  
  if (fp->type != GTS_ERROR) {
    if (sim == NULL)
      gts_file_error (fp, "hypre module must be called within the GfsSimulation parameter block");
    else {
      /* initialise the poisson solver hook */
      sim->approx_projection_params.poisson_solve = hypre_poisson_solve;
      sim->projection_params.poisson_solve = hypre_poisson_solve;
      /* initialise the other diffusion and Poisson solver hook(s) */
      sim->advection_params.diffusion_solve = hypre_diffusion_solve;
      GSList * i = GFS_DOMAIN (sim)->variables;
      while (i) {
	if (GFS_IS_VARIABLE_TRACER (i->data))
	  GFS_VARIABLE_TRACER (i->data)->advection.diffusion_solve = hypre_diffusion_solve;
	else if (GFS_IS_VARIABLE_POISSON (i->data))
	  GFS_VARIABLE_POISSON (i->data)->par.poisson_solve = hypre_poisson_solve;
	i = i->next;
      }
    }
  }
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);

  hypre_solver_write (&proj_hp, fp);
}
