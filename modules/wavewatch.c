/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "init.h"
#include "wave.h"
#include "config.h"
#if WW3_VERSION == 222
#  include "wavewatch/wavewatch_222.h"
#elif WW3_VERSION == 312
#  include "wavewatch/wavewatch_312.h"
#else /* 3.14 */
#  include "wavewatch/wavewatch_314.h"
#endif

#define DEBUG 0

static double frequency (int ik)
{
  return GFS_WAVE_F0*pow(GFS_WAVE_GAMMA, ik);
}

#if DEBUG      
static double theta (guint ith, guint ntheta)
{
  return 2.*M_PI*ith/ntheta;
}
#endif

typedef struct {
  GfsWave * wave;
  GfsVariable * ustar, * fpi, * u10, * v10, * as;
  REAL * A;     /* Actions (NK*NTHETA) */
  REAL * CG;    /* Group velocities (NK) */
  REAL * WN;    /* Wavenumbers (NK) */
  REAL * ALPHA; /* Nondimensional 1-D spectrum (NK) */  
} SourceParams;

static void energy_to_action (FttCell * cell, SourceParams * p)
{
  guint i, j;
  REAL * A = p->A;
  for (i = 0; i < p->wave->nk; i++)
    for (j = 0; j < p->wave->ntheta; j++) {
      *A = GFS_VALUE (cell, p->wave->F[i][j])*p->CG[i]/(2.*M_PI*frequency (i));
      if (*A < 0.)
        *A = 0.;
      A++;
    }
}

static void action_to_energy (FttCell * cell, SourceParams * p)
{
  guint i, j;
  REAL * A = p->A;
  for (i = 0; i < p->wave->nk; i++)
    for (j = 0; j < p->wave->ntheta; j++) {
      GFS_VALUE (cell, p->wave->F[i][j]) = *A*(2.*M_PI*frequency (i))/p->CG[i];
      A++;
    }
}

static void stability_correction (FttCell * cell, SourceParams * p, REAL * U10ABS, REAL * U10DIR)
{
  double u10 = p->u10 ? GFS_VALUE (cell, p->u10) : 0.;
  double v10 = p->v10 ? GFS_VALUE (cell, p->v10) : 0.;
  *U10ABS = sqrt (u10*u10 + v10*v10);
  *U10DIR = atan2 (v10, u10);
  if (p->as) {
    /* see p.24 of wavewatch manual version 3.12 */
    double shstab = 1.4, ofstab = -0.01, ccng = -0.1, ccps = 0.1, ffng = -150., ffps = 150.;
    double zwind = 10.;  /* height at which the wind is defined */
    double grav = 9.81; /* acceleration of gravity */
    double stab0 = zwind*grav/273.;
    double max = MAX (5., *U10ABS);
    double stab = stab0*GFS_VALUE (cell, p->as)/(max*max);
    stab = MAX (-1., MIN (1., stab));
    double tharg1 = MAX (0., ffng*(stab - ofstab));
    double tharg2 = MAX (0., ffps*(stab - ofstab));
    double cor1 = ccng*tanh (tharg1);
    double cor2 = ccps*tanh (tharg2);
    double cor = sqrt ((1. + cor1 + cor2)/shstab);
    *U10ABS /= cor;
  }
}

static void source (FttCell * cell, SourceParams * p)
{
  energy_to_action (cell, p);

  REAL DEPTH = 1000.; /* fixme: depth is fixed at 1000 m for now */
  REAL U10ABS, U10DIR;

  stability_correction (cell, p, &U10ABS, &U10DIR);

  REAL USTAR = GFS_VALUE (cell, p->ustar);
  REAL FPI = GFS_VALUE (cell, p->fpi);
  REAL DTG = GFS_SIMULATION (p->wave)->advection_params.dt*3600.;

  REAL EMEAN, FMEAN, WMEAN, AMAX;
  REAL CD, Z0;
  REAL DTDYN, FCUT;
  
#if WW3_VERSION == 222
  REAL DTMIN = DTG/10., DTMAX = DTG;
  W3SRCE (p->A, p->A, p->ALPHA, p->WN, p->CG, &DEPTH, 
	  &U10ABS, &U10DIR, &USTAR,
	  &EMEAN, &FMEAN, &WMEAN, &AMAX, 
	  &FPI, &CD, &Z0, 
	  &DTDYN, &FCUT, &DTG, &DTMIN, &DTMAX);
#elif WW3_VERSION == 312
  INTEGER IX, IY, IMOD = 1;
  REAL USTDIR;
  W3SRCE (&IX, &IY, &IMOD, p->A, p->ALPHA, p->WN, p->CG, &DEPTH, 
	  &U10ABS, &U10DIR, &USTAR, &USTDIR,
	  &EMEAN, &FMEAN, &WMEAN, &AMAX, 
	  &FPI, &CD, &Z0, 
	  &DTDYN, &FCUT, &DTG);
#else /* 3.14 */
  INTEGER IX, IY, IMOD = 1;
  REAL AS = p->as ? GFS_VALUE (cell, p->as) : 0.;
  REAL USTDIR;
  REAL CX = 0., CY = 0.;
  W3SRCE (&IX, &IY, &IMOD, p->A, p->ALPHA, p->WN, p->CG, &DEPTH, 
	  &U10ABS, &U10DIR, &AS, &USTAR, &USTDIR,
	  &CX, &CY,
	  &EMEAN, &FMEAN, &WMEAN, &AMAX, 
	  &FPI, &CD, &Z0, 
	  &DTDYN, &FCUT, &DTG);
#endif /* 3.14 */

#if DEBUG
  guint i, j;
  for (i = 0; i < p->wave->nk; i++) {
    for (j = 0; j < p->wave->ntheta; j++)      
      fprintf (stderr, "%g %g %g\n", frequency (i), theta (j, p->wave->ntheta), 
	       p->A[j + i*p->wave->ntheta]);
    fprintf (stderr, "\n");
  }
#endif

  action_to_energy (cell, p);
  GFS_VALUE (cell, p->ustar) = USTAR;
  GFS_VALUE (cell, p->fpi) = FPI;
}

static void deletedir (const char * name)
{
  gchar * command = g_strconcat ("rm -r -f ", name, NULL);
  int status = system (command);
  if (status)
    g_warning ("could not cleanup wavewatch setup");
  g_free (command);
}

static void initialize (GfsWave * wave)
{
  static gboolean initialized = FALSE;
  if (!initialized) {

    /* Creates temporary directory */
    gchar * template = gfs_template ();
    if (!g_mkdtemp (template)) {
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	     "wavewatch module: could not create temporary directory\n%s",
	     strerror (errno));
      return;
    }

    /* Creates wavewatch ww3_grid.inp input file */
    gchar * sinput = g_strconcat (template, "/ww3_grid.inp", NULL);
    FILE * input = fopen (sinput, "w");
    g_free (sinput);
    if (input == NULL) {
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	     "wavewatch module: could not create input file\n"
	     "%s\n", strerror (errno));
      deletedir (template);
      return;      
    }
    fprintf (input,
	     "$\n"
	     "      'Gerris wavewatch module'\n"
	     "   %g  %g  %d  %d  0.\n"
	     "   F F F F F T\n"
	     "    900. 950. 900. 300.\n"
	     "END OF NAMELISTS\n",
	     GFS_WAVE_GAMMA, GFS_WAVE_F0, wave->nk, wave->ntheta);

    /* Dummy wavewatch parameters */
    static gchar constant_parameters[] =
      "      3      3\n"
      "      1000.     1000.     1.\n"
      "     -1000.    -1000.     1.\n"
      "     -0.1 2.50  10  -1000. 3 1 '(....)' 'NAME' 'bottom.inp'\n"
      "  1 1 1\n"
      "  1 1 1\n"
      "  1 1 1\n"
#if WW3_VERSION == 222
      "      0   0   F\n"
#else /* version 3.12 and 3.14 */
      "   10 3 1 '(....)' 'PART' 'mapsta.inp'\n"
      "      0   0   F\n"
      "      0   0   F\n"
      "      0   0\n"
#endif /* version 3.12 and 3.14 */
      "     0.    0.    0.    0.       0\n";
    fputs (constant_parameters, input);
    fclose (input);

    /* Calls 'ww3_grid' of wavewatch to generate 'mod_def.w3'
     * required to initialize wavewatch. */
    char * wdir = getcwd (NULL, 0);
    char * command = g_strconcat ("cd ",  template, " && "
				  "test -f $HOME/.wwatch3.env && "
				  "`grep WWATCH3_DIR $HOME/.wwatch3.env | "
				  "awk '{print $2}'`/exe/ww3_grid > ", wdir, "/log_grid.ww3 && "
				  "mv mod_def.ww3 ", wdir, " && "
				  "rm -r -f ", template,
				  NULL);
    int status = system (command);
    deletedir (template);
    g_free (template);
    free (wdir);
    g_free (command);
    if (status == -1 || WEXITSTATUS (status) != 0) {
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	     "wavewatch module: error when running ww3_grid\nsee log_grid.ww3 for details");
      return;
    }

    /* Initialize wavewatch */
    GFSW3INIT ();

    /* cleanup */
    remove ("mod_def.ww3");
    initialized = TRUE;
  }
}

static void wavewatch_source (GfsWave * wave)
{
  GfsDomain * domain = GFS_DOMAIN (wave);
  SourceParams p;

  gfs_domain_timer_start (domain, "wavewatch_source");

  initialize (wave);

  p.wave = wave;
  p.A = g_malloc (wave->nk*wave->ntheta*sizeof (REAL));
  p.CG = g_malloc (wave->nk*sizeof (REAL));
  p.WN = g_malloc (wave->nk*sizeof (REAL));
  p.ALPHA = g_malloc (wave->nk*sizeof (REAL));
  p.ustar = gfs_variable_from_name (domain->variables, "Ustar");
  p.fpi = gfs_variable_from_name (domain->variables, "Fpi");
  p.u10 = gfs_variable_from_name (domain->variables, "U10");
  p.v10 = gfs_variable_from_name (domain->variables, "V10");
  p.as = gfs_variable_from_name (domain->variables, "AS");

  guint i;
  for (i = 0; i < wave->nk; i++) {
    REAL omega = 2.*M_PI*frequency (i);
    p.WN[i] = omega*omega/9.81;
    p.CG[i] = 9.81/omega/2.;
  }

  gfs_catch_floating_point_exceptions ();
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) source, &p);
  if (gfs_restore_floating_point_exceptions ())
    g_warning ("floating-point exceptions caught in wavewatch module");

  g_free (p.A);
  g_free (p.CG);
  g_free (p.WN);
  g_free (p.ALPHA);

  gfs_domain_timer_stop (domain, "wavewatch_source");
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "wavewatch";
const gchar * g_module_check_init (void);
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);
  g_return_if_fail (sim != NULL);

  if (!GFS_IS_WAVE (sim)) {
    gts_file_error (fp, "wavewatch module can only be used with GfsWave");
    return;
  }

  GFS_WAVE (sim)->source = wavewatch_source;
  gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), "Ustar", "Friction velocity");
  gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), "Fpi",   "Peak-input frequency");
}
