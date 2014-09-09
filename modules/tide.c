/* Gerris - The GNU Flow Solver
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

#include <string.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>
#include "fes2004/fes2004_lib.h"
#include "variable.h"

/* Heavily based on GfsBcFlather */

static gchar * reference = "1950/01/01-00:00:00-UTC";
static gdouble deltat = 0.;

/* GfsBcTide: Header */

typedef struct _GfsBcTide         GfsBcTide;

struct _GfsBcTide {
  /*< private >*/
  GfsBcValue parent;
  gdouble ** amplitude, ** phase, x, size;

  /*< public >*/
  GfsVariable * h, * p;
};

#define GFS_BC_TIDE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsBcTide,\
					         gfs_bc_tide_class ())
#define GFS_IS_BC_TIDE(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_tide_class ()))

GfsBcClass * gfs_bc_tide_class  (void);

/* GfsBcTide: Object */

#define N  64 /* number of discretisation points */

#define NM 14 /* number of tidal modes (must match those of FES2004) */

static void bc_tide_write (GtsObject * o, FILE * fp)
{
  GfsBcTide * bc = GFS_BC_TIDE (o);
  guint i, j;

  (* GTS_OBJECT_CLASS (gfs_bc_tide_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %s {\n", bc->h->name, bc->p->name);
  for (i = 0; i < N; i++) {
    for (j = 0; j < NM; j++)
      fprintf (fp, "  %g %g\n", bc->amplitude[j][i], bc->phase[j][i]);
  }
  fputc ('}', fp);
}

static void set_gradient_boundary (FttCell * cell)
{
  cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
}

static void bc_tide_read (GtsObject ** o, GtsFile * fp)
{
  GfsBcTide * bc = GFS_BC_TIDE (*o);
  GfsDomain * domain = gfs_box_domain (GFS_BC (bc)->b->box);

  (* GTS_OBJECT_CLASS (gfs_bc_tide_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsBoundary * b = GFS_BC (bc)->b;
  if (b->d > FTT_BOTTOM) {
    gts_file_error (fp, "GfsBcTide cannot be used for 3D boundaries");
    return;
  }

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (h)");
    return;
  }
  bc->h = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->h == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (p)");
    return;
  }
  bc->p = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->p == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  ftt_cell_traverse (b->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) set_gradient_boundary, NULL);

  guint i;  
  gpointer tmp = g_malloc0 (N*NM*sizeof (double));
  bc->amplitude = g_malloc (N*sizeof (double *));
  for (i = 0; i < N; i++)
    bc->amplitude[i] = tmp + i*NM*sizeof (double);
  tmp = g_malloc0 (N*NM*sizeof (double));
  bc->phase = g_malloc (N*sizeof (double *));
  for (i = 0; i < N; i++)
    bc->phase[i] = tmp + i*NM*sizeof (double);

  FttCellFace face;
  face.cell = b->root;
  face.d = b->d;
  FttVector p;
  ftt_face_pos (&face, &p);
  FttComponent c = face.d < FTT_TOP ? FTT_Y : FTT_X;
  bc->size = ftt_cell_size (b->root);
  (&p.x)[c] -= bc->size/2.;
  bc->x = (&p.x)[c];

  if (fp->type == '{') {
    /* read embedded coefficients */
    guint j;
    fp->scope_max++;
    gts_file_next_token (fp);
    for (i = 0; i < N; i++) {
      for (j = 0; j < NM; j++) {
	while (fp->type == '\n') gts_file_next_token (fp);
	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	  gts_file_error (fp, "expecting a number (amplitude[%d][%d])", j, i);
	  return;
	}
	bc->amplitude[j][i] = atof (fp->token->str);
	gts_file_next_token (fp);
	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	  gts_file_error (fp, "expecting a number (phase[%d][%d])", j, i);
	  return;
	}
	bc->phase[j][i] = atof (fp->token->str);
	gts_file_next_token (fp);
      }
    }
    while (fp->type == '\n') gts_file_next_token (fp);
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
    gts_file_next_token (fp);
  }
  else {
    /* extract FES2004 tidal coefficients */
    gchar * fname = getenv ("GFS_FES2004") ?
      g_strdup (getenv ("GFS_FES2004"))
      :
      g_strjoin ("/", GFS_MODULES_DIR, "tide.fes2004.nc", NULL);
    FILE * f = fopen (fname, "r");
    if (f == NULL) {
      gts_file_error (fp, "cannot open file `%s': %s", fname, strerror (errno));
      g_free (fname);
      return;
    }
    fclose (f);

    double * lon = g_malloc (N*sizeof (double));
    double * lat = g_malloc (N*sizeof (double));
    gdouble dh = bc->size/(N - 1);
    for (i = 0; i < N; i++, (&p.x)[c] += dh) {
      FttVector mp = p;
      gfs_simulation_map_inverse (GFS_SIMULATION (gfs_box_domain (b->box)), &mp);
      lon[i] = mp.x; lat[i] = mp.y;
    }

    fes2004_extraction (fname, N, lat, lon, bc->amplitude, bc->phase, 1);

    g_free (fname);
    g_free (lon);
    g_free (lat);
  }
}

static void bc_tide_destroy (GtsObject * o)
{
  GfsBcTide * bc = GFS_BC_TIDE (o);

  if (bc->amplitude) {
    g_free (bc->amplitude[0]);
    g_free (bc->amplitude);
  }
  if (bc->phase) {
    g_free (bc->phase[0]);
    g_free (bc->phase);
  }
    
  (* GTS_OBJECT_CLASS (gfs_bc_tide_class ())->parent_class->destroy) (o);
}

static tidal_wave wave[NM];

static gdouble amplitude_value (FttCellFace * face, GfsBcTide * bc, gdouble t)
{
  FttComponent c = face->d < FTT_TOP ? FTT_Y : FTT_X;
  FttVector p;
  ftt_face_pos (face, &p);
  guint i = floor (((&p.x)[c] - bc->x)/bc->size*(N - 1));
  g_assert (i < N - 1);
  if (bc->amplitude[i][2] < 0. && bc->amplitude[i+1][2] < 0.)
    return G_MAXDOUBLE;

  gdouble dh = bc->size/(N - 1);
  gdouble a = ((&p.x)[c] - bc->x - dh*i)/dh;
  if (bc->amplitude[i][2] < 0.)
    a = 1.;
  if (bc->amplitude[i+1][2] < 0.)
    a = 0.;

  gdouble prediction = 0.;
  guint j;
  for (j = 0; j < NM; j++) {
    fcomplex Z1, Z2, Z;
    Z1.reel = bc->amplitude[i][j]*cos (- bc->phase[i][j]*M_PI/180.);
    Z1.imag = bc->amplitude[i][j]*sin (- bc->phase[i][j]*M_PI/180.);
    Z2.reel = bc->amplitude[i+1][j]*cos (- bc->phase[i+1][j]*M_PI/180.);
    Z2.imag = bc->amplitude[i+1][j]*sin (- bc->phase[i+1][j]*M_PI/180.);
    Z.reel = (1. - a)*Z1.reel + a*Z2.reel;
    Z.imag = (1. - a)*Z1.imag + a*Z2.imag;
    prediction += Tide_prediction (t, wave[j], Z, 0, 0);
  }
  return prediction;
}

static gdouble tide_value (FttCellFace * f, GfsBc * b)
{
  /* fixme: this will not work for multilayer domains */
  guint d, nb = 0;
  FttCellNeighbors n;
  gdouble H;

  ftt_cell_neighbors (f->neighbor, &n);
  for (d = 0; d < FTT_NEIGHBORS_2D; d++)
    if (n.c[d] != NULL && GFS_CELL_IS_BOUNDARY(n.c[d]) && nb++ > 0)
      /* if the boundary cell is bounded by more than one boundary -> no flux */
      return 0.;

  H = gfs_face_interpolated_value (f, GFS_BC_TIDE (b)->h->i);
  if (H > 2e-3) { /* fixme: 2e-3 is an arbitrary constant which should be a parameter or sthg*/
    GfsSimulation * sim = GFS_SIMULATION (gfs_box_domain (b->b->box));
    gdouble a = amplitude_value (f, GFS_BC_TIDE (b), sim->time.t + deltat);
    if (a < G_MAXDOUBLE) {
      gdouble cg = sqrt (sim->physical_params.g*H);
      
      a *= sim->physical_params.g/5000.; /* fixme: reference depth is fixed at 5000 meters */
      return gfs_function_face_value (GFS_BC_VALUE (b)->val, f) +
	(FTT_FACE_DIRECT (f) ? -1. : 1.)*
	(GFS_VALUE (f->neighbor, GFS_BC_TIDE (b)->p) - a)*
	cg/sim->physical_params.g
#if !FTT_2D
	/H
#endif
	;
    }
  }
  return 0.;
}

static void tide (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  g_assert (ftt_cell_neighbor (f->cell, f->d) == f->neighbor);
  GFS_VALUE (f->cell, b->v) = 2.*tide_value (f, b) - GFS_VALUE (f->neighbor, b->v);
}

static void homogeneous_tide (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_VALUE (f->cell, b->v) = - GFS_VALUE (f->neighbor, b->v);
}

static void face_tide (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_STATE (f->cell)->f[f->d].v = tide_value (f, b);
}

static void gfs_bc_tide_class_init (GtsObjectClass * klass)
{
  klass->write   = bc_tide_write;
  klass->read    = bc_tide_read;
  klass->destroy = bc_tide_destroy;
}

static void gfs_bc_tide_init (GfsBc * bc)
{
  bc->bc =             (FttFaceTraverseFunc) tide;
  bc->homogeneous_bc = (FttFaceTraverseFunc) homogeneous_tide;
  bc->face_bc =        (FttFaceTraverseFunc) face_tide;
}

GfsBcClass * gfs_bc_tide_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_tide_info = {
      "GfsBcTide",
      sizeof (GfsBcTide),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_tide_class_init,
      (GtsObjectInitFunc) gfs_bc_tide_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_tide_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "tide";
const gchar * g_module_check_init (void);
void          gfs_module_read     (GtsFile * fp);
void          gfs_module_write    (FILE * fp);

const gchar * g_module_check_init (void)
{
  wave[0] = w2N2;
  wave[1] = wK1;
  wave[2] = wK2;
  wave[3] = wM2;
  wave[4] = wM4;
  wave[5] = wMf;
  wave[6] = wMm;
  wave[7] = wMSqm;
  wave[8] = wMtm;
  wave[9] = wN2;
  wave[10] = wO1;
  wave[11] = wP1;
  wave[12] = wQ1;
  wave[13] = wS2;
  gfs_bc_tide_class ();
  return NULL;
}

#define TIME_FORMAT "%Y/%m/%d-%T"

void gfs_module_read (GtsFile * fp)
{
  g_return_if_fail (fp != NULL);
  
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "reference",  TRUE},
      {GTS_NONE}
    };
    var[0].data = &reference;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;

    if (var[0].set) {
      struct tm timeref;
      time_t tref, t0;
      memset (&timeref, 0, sizeof(struct tm));
      strptime ("1950/01/01-00:00:00-UTC", TIME_FORMAT, &timeref);
      tref = mktime (&timeref);
      memset (&timeref, 0, sizeof(struct tm));
      if (!strptime (reference, TIME_FORMAT, &timeref)) {
	gts_file_variable_error (fp, var, "reference", "error parsing date format");
	return;
      }
      t0 = mktime (&timeref);
      deltat = difftime (t0, tref)/3600.;
    }
  }
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);
 
  fprintf (fp, " { reference = %s }", reference);
}
