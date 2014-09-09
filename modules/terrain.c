/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
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
#include <glob.h>
#if GSL
# include <gsl/gsl_linalg.h>
#endif
#include "refine.h"
#include "solid.h"
#include "kdt/kdt.h"
#include "river.h"

static gchar * default_path = ".";

/* Kdtrees */

typedef struct {
  Kdt ** rs;
  gdouble * weight;
  guint nrs;
  gchar * path, * basename;
} Kdtrees;

static void kdtrees_destroy (Kdtrees * rs)
{
  g_free (rs->path);
  g_free (rs->basename);
  if (rs->rs) {
    guint i;
    for (i = 0; i < rs->nrs; i++)
      kdt_destroy (rs->rs[i]);
    g_free (rs->rs);
  }
  g_free (rs->weight);
}

static Kdt * open_kdt (const gchar * fname)
{
  Kdt * kdt = kdt_new ();
  if (kdt_open (kdt, fname)) {
    kdt_destroy (kdt);
    gchar * name = g_strconcat (fname, ".DataPD", NULL);
    FILE * fp = fopen (name, "r");
    g_free (name);
    if (fp != NULL) {
      fclose (fp);
      g_warning ("\nFound obsolete R*-tree terrain database. Use:\n"
		 "%% rsurface2kdt -v %s\n"
		 "to convert to the new KDT format.\n", fname);
    }
    return NULL;
  }
  return kdt;
}

static void kdtrees_read (Kdtrees * rs, GtsFile * fp)
{
  gchar * path = NULL;
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "basename", TRUE, &rs->basename},
      {GTS_STRING, "path",     TRUE, &rs->path},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    path = g_strconcat (default_path, ":", rs->path, NULL);
  }
  else
    path = g_strdup (default_path);

  if (!strcmp (rs->basename, "*")) { /* file globbing */
    gchar * pattern = g_strconcat (path, "/*.kdt", NULL);
    glob_t pglob;
    if (glob (pattern, GLOB_ERR, NULL, &pglob)) {
      gts_file_error (fp, "cannot find/open terrain databases in path:\n%s", pattern);
      g_free (pattern);
      g_free (path);
      return;
    }
    g_free (pattern);
    guint i;
    for (i = 0; i < pglob.gl_pathc; i++) {
      pglob.gl_pathv[i][strlen (pglob.gl_pathv[i]) - 5] = '\0';
      rs->rs = g_realloc (rs->rs, (rs->nrs + 1)*sizeof (Kdt *));
      rs->weight = g_realloc (rs->weight, (rs->nrs + 1)*sizeof (double));
      rs->rs[rs->nrs] = open_kdt (pglob.gl_pathv[i]);
      if (!rs->rs[rs->nrs]) {
	gts_file_error (fp, "cannot open terrain database `%s'", pglob.gl_pathv[i]);
	globfree (&pglob);
	g_free (path);
	return;
      }
      rs->weight[rs->nrs] = 1.;
      rs->nrs++;
    }
    globfree (&pglob);
  }
  else { /* basename is of the form: set1:w1,set2:w2,set3:w3... */
    gchar ** names = g_strsplit (rs->basename, ",", 0);
    gchar ** s = names;
    while (*s) {
      /* look for weight */
      gchar ** wname = g_strsplit (*s, ":", 2);
      rs->rs = g_realloc (rs->rs, (rs->nrs + 1)*sizeof (Kdt *));
      rs->weight = g_realloc (rs->weight, (rs->nrs + 1)*sizeof (double));
      if (path) {
	/* search path */
	gchar ** pathes = g_strsplit (path, ":", 0);
	gchar ** spath = pathes, * fname;
	g_assert (*spath);
	do {
	  fname = (*wname)[0] == '/' ? g_strdup (*wname) : g_strconcat (*spath, "/", *wname, NULL);
	  rs->rs[rs->nrs] = open_kdt (fname);
	} while (rs->rs[rs->nrs] == NULL && *(++spath));
	g_strfreev (pathes);
      }
      else
	rs->rs[rs->nrs] = open_kdt (*wname);
      if (!rs->rs[rs->nrs]) {
	if (path)
	  gts_file_error (fp, "cannot find/open terrain database `%s' in path:\n%s", *wname, path);
	else
	  gts_file_error (fp, "cannot open terrain database `%s'", *wname);
	g_strfreev (wname);
	g_strfreev (names);
	g_free (path);
	return;
      }
      if (wname[1])
	rs->weight[rs->nrs] = strtod (wname[1], NULL);
      else
	rs->weight[rs->nrs] = 1.;
      g_strfreev (wname);
      rs->nrs++;
      s++;
    }
    g_strfreev (names);
  }
  g_free (path);
}

static void kdtrees_write (Kdtrees * rs, FILE * fp)
{
  if (rs->path || rs->basename) {
    fputs (" {\n", fp);
    if (rs->path)
      fprintf (fp, "  path = %s\n", rs->path);
    if (rs->basename)
      fprintf (fp, "  basename = %s\n", rs->basename);
    fputc ('}', fp);
  }
}

/* GfsRefineTerrain: Header */

typedef struct _GfsRefineTerrain         GfsRefineTerrain;

#define NM 4

struct _GfsRefineTerrain {
  /*< private >*/
  GfsRefine parent;
  guint level;
  gboolean refined;
  GfsVariable * type;

#if !FTT_2D
  GfsVariable * min, * max;
  gdouble front, scale;
#endif

  Kdtrees rs;

  /*< public >*/
  gchar * name;
  GfsVariable * h[NM], * he, * hn, * hdmin, * hdmax;
  GfsFunction * criterion;
};

#define GFS_REFINE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineTerrain,\
					           gfs_refine_terrain_class ())
#define GFS_IS_REFINE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_terrain_class ()))
     
GfsRefineClass * gfs_refine_terrain_class  (void);

/* GfsRefineTerrain: Object */

typedef struct {
  FttVector c;
  FttVector p[4];
  gdouble min[2], max[2], h;
  Kdtrees * rs;
  FttCell * cell;
} Polygon;

static void polygon_init (GfsSimulation * sim, Polygon * poly, FttCell * cell, Kdtrees * rs)
{
  FttVector q;
  ftt_cell_pos (cell, &q);
  poly->cell = cell;
  poly->rs = rs;
  poly->h = ftt_cell_size (cell)/2.;
  poly->p[0].x = q.x + poly->h; poly->p[0].y = q.y + poly->h; poly->p[0].z = 0.;
  poly->p[1].x = q.x - poly->h; poly->p[1].y = q.y + poly->h; poly->p[1].z = 0.;
  poly->p[2].x = q.x - poly->h; poly->p[2].y = q.y - poly->h; poly->p[2].z = 0.;
  poly->p[3].x = q.x + poly->h; poly->p[3].y = q.y - poly->h; poly->p[3].z = 0.;
  gfs_simulation_map_inverse_cell (sim, poly->p);

  poly->c.x = poly->c.y = 0.;
  poly->min[0] = poly->min[1] = G_MAXDOUBLE;
  poly->max[0] = poly->max[1] = - G_MAXDOUBLE;
  gint i;
  FttVector * p = poly->p;
  for (i = 0; i < 4; i++, p++) {
    if (p->x < poly->min[0]) poly->min[0] = p->x;
    if (p->x > poly->max[0]) poly->max[0] = p->x;
    if (p->y < poly->min[1]) poly->min[1] = p->y;
    if (p->y > poly->max[1]) poly->max[1] = p->y;
    poly->c.x += p->x; poly->c.y += p->y;    
  }
  poly->c.x /= 4; poly->c.y /= 4;
  poly->h = MAX (poly->max[0] - poly->min[0], poly->max[1] - poly->min[1])/2.;
}

static gboolean right (const double a[2], const double b[2], const double c[2])
{
  return (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]) < 0.;
}

static gboolean polygon_contains (Polygon * p, gdouble q[2])
{
  if (right (&p->p[0].x, &p->p[1].x, q))
    return FALSE;
  if (right (&p->p[1].x, &p->p[2].x, q))
    return FALSE;
  if (right (&p->p[2].x, &p->p[3].x, q))
    return FALSE;
  if (right (&p->p[3].x, &p->p[0].x, q))
    return FALSE;
  return TRUE;
}

static gboolean polygon_includes (KdtRect rect, Polygon * p)
{
  gdouble q[2];
  q[0] = rect[0].l; q[1] = rect[1].l;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].l; q[1] = rect[1].h;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].h; q[1] = rect[1].l;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].h; q[1] = rect[1].h;
  if (!polygon_contains (p, q))
    return FALSE;
  return TRUE;
}

static gboolean polygon_intersects (KdtRect rect, Polygon * p)
{
  /* fixme: this could be improved? */
  return (rect[0].l <= p->max[0] && rect[0].h >= p->min[0] &&
	  rect[1].l <= p->max[1] && rect[1].h >= p->min[1]);
}

typedef struct {
  gdouble H[NM+1], m[NM][NM];
  gdouble h[NM], he, cond, min, max;
  GfsRefineTerrain * t;
  FttCell * cell;
  gboolean relative;
  int n;
} RMS;

static void rms_init (GfsRefineTerrain * t, RMS * rms, Polygon * p, gboolean relative)
{
  guint i, j;
  for (i = 0; i < NM + 1; i++)
    rms->H[i] = 0.;
  for (i = 0; i < NM; i++)
    for (j = 0; j < NM; j++)
      rms->m[i][j] = 0.;
  rms->t = t;
  rms->cell = p->cell;
  rms->relative = relative;
  rms->min = G_MAXDOUBLE;
  rms->max = - G_MAXDOUBLE;
}

static void function_from_corners (gdouble h[4], gdouble H[4])
{
  h[0] = (H[0] + H[1] + H[2] + H[3])/4.;
  h[1] = (H[0] - H[1] - H[2] + H[3])/4.;
  h[2] = (H[0] + H[1] - H[2] - H[3])/4.;
  h[3] = (H[0] - H[1] + H[2] - H[3])/4.;  
}

static gdouble rms_minimum (RMS * rms)
{
  if (rms->m[0][0] == 0.)
    return 0.;
  return sqrt (fabs (rms->h[0]*(rms->h[0]*rms->m[0][0] + 
				2.*(rms->h[1]*rms->m[0][1] + 
				    rms->h[2]*rms->m[0][2] +
				    rms->h[3]*rms->m[0][3] - rms->H[0])) +
		     rms->h[1]*(rms->h[1]*rms->m[1][1] + 
				2.*(rms->h[2]*rms->m[1][2] +
				    rms->h[3]*rms->m[1][3] - rms->H[1])) +
		     rms->h[2]*(rms->h[2]*rms->m[2][2] +
				2.*(rms->h[3]*rms->m[2][3] - rms->H[2])) +
		     rms->h[3]*(rms->h[3]*rms->m[3][3] - 2.*rms->H[3]) +
		     rms->H[4])/rms->m[0][0]);
}

static gdouble cell_value (FttCell * cell, GfsVariable * h[NM], FttVector p)
{
  if (GFS_VALUE (cell, h[0]) == GFS_NODATA)
    return GFS_NODATA;
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector q;
  ftt_cell_pos (cell, &q);
  p.x = (p.x - q.x)/size;
  p.y = (p.y - q.y)/size;
  return (GFS_VALUE (cell, h[0]) + 
	  GFS_VALUE (cell, h[1])*p.x + 
	  GFS_VALUE (cell, h[2])*p.y + 
	  GFS_VALUE (cell, h[3])*p.x*p.y);
}

static void corners_from_parent (FttCell * cell, GfsRefineTerrain * t, gdouble H[4])
{
  gdouble size = ftt_cell_size (cell);
  FttCell * parent = ftt_cell_parent (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += size/2.; p.y += size/2.;
  H[0] = cell_value (parent, t->h, p);
  p.x -= size;
  H[1] = cell_value (parent, t->h, p);
  p.y -= size;
  H[2] = cell_value (parent, t->h, p);
  p.x += size;
  H[3] = cell_value (parent, t->h, p);
}

static void variance_check (RMS * rms)
{
  g_assert (rms->n >= NM);
  gdouble H[4], h[4];
  guint i;
  h[0] = rms->h[0] + rms->h[1] + rms->h[2] + rms->h[3];
  h[1] = rms->h[0] - rms->h[1] + rms->h[2] - rms->h[3];
  h[2] = rms->h[0] - rms->h[1] - rms->h[2] + rms->h[3];
  h[3] = rms->h[0] + rms->h[1] - rms->h[2] - rms->h[3];
  if (rms->relative) {
    gdouble H0[4];
    corners_from_parent (rms->cell, rms->t, H0);
    for (i = 0; i < 4; i++)
      H[i] = CLAMP (h[i], rms->min - H0[i], rms->max - H0[i]);
  }
  else
    for (i = 0; i < 4; i++)
      H[i] = CLAMP (h[i], rms->min, rms->max);
  function_from_corners (rms->h, H);
}

static int rms_update (RMS * rms)
{
  guint i;
  if (rms->m[0][0] == 0.) {
    for (i = 1; i < NM; i++)
      rms->h[i] = 0.;
    rms->h[0] = GFS_NODATA;
    rms->he = 0.;
    rms->cond = GFS_NODATA;
    rms->n = 0;
    return 0;
  }
  else if (rms->n >= NM) {
    guint j;
    for (i = 1; i < NM; i++)
      for (j = 0; j < i; j++)
	rms->m[i][j] = rms->m[j][i];
#if GSL  
    double m[NM*NM], v[NM*NM], s[NM];
    gsl_matrix_view gv = gsl_matrix_view_array (v, NM, NM);
    gsl_vector_view gs = gsl_vector_view_array (s, NM);
    for (i = 0; i < NM; i++)
      for (j = 0; j < NM; j++)
	m[i+NM*j] = rms->m[i][j];
    gsl_matrix_view gm = gsl_matrix_view_array (m, NM, NM);
    gsl_linalg_SV_decomp_jacobi (&gm.matrix, &gv.matrix, &gs.vector);
    rms->cond = s[NM - 1] > 0. ? s[0]/s[NM - 1] : GFS_NODATA;
    if (rms->cond < 10000.) {
      gsl_vector_view gH = gsl_vector_view_array (rms->H, NM);
      gsl_vector_view gh = gsl_vector_view_array (rms->h, NM);
      gsl_linalg_SV_solve (&gm.matrix, &gv.matrix, &gs.vector, &gH.vector, &gh.vector);
      variance_check (rms);
      rms->he = rms_minimum (rms);
      return 1;
    }
#else
    gdouble ** m = gfs_matrix_new (NM, NM, sizeof (gdouble));
    for (i = 0; i < NM; i++)
      for (j = 0; j < NM; j++)
	m[i][j] = rms->m[i][j];
    if (gfs_matrix_inverse (m, NM, 1e-5)) {
      for (i = 0; i < NM; i++) {
	rms->h[i] = 0.;
	for (j = 0; j < NM; j++)
	  rms->h[i] += m[i][j]*rms->H[j];
      }
      gfs_matrix_free (m);
      variance_check (rms);
      rms->he = rms_minimum (rms);
      return 1;
    }
    gfs_matrix_free (m);
#endif
  }
  rms->h[0] = rms->H[0]/rms->m[0][0];
  for (i = 1; i < NM; i++)
    rms->h[i] = 0.;
  rms->he = rms_minimum (rms);
  return 0;
}

#if DEBUG
static gdouble rms_value (RMS * rms, FttVector * p)
{
  return rms->h[0] + rms->h[1]*p->x + rms->h[2]*p->y + rms->h[3]*p->x*p->y;
}

static void rms_write (RMS * rms, Polygon * p)
{
  FttVector q, r;
  q.x = p->c.x + p->h; q.y = p->c.y + p->h;
  r.x = 1.; r.y = 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x + p->h; q.y = p->c.y - p->h;
  r.x = 1.; r.y = - 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x - p->h; q.y = p->c.y - p->h;
  r.x = - 1.; r.y = - 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x - p->h; q.y = p->c.y + p->h;
  r.x = - 1.; r.y = 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
}

static int write_points (double p[3], Polygon * poly)
{
  if (polygon_contains (poly, p))
    fprintf (stderr, "aa %g %g %g\n", p[0], p[1], p[2]);
}
#endif

#define RAW        0. /* fitted but not continuous (C0) */
#define FAIR       1. /* fitted and C0 */
#define REFINED    2. /* non-fitted refined cell */
#define NEW_CHILD  3. /* non-fitted child extrapolated from its parent */
#define CONTAINS_SURFACE 4. /* 3D-only */
#define BOUNDARY 5.         /* 3D-only */

static void parent_cell_coefficients (FttCell * cell, GfsVariable ** v, gdouble hP[NM])
{
  FttCell * parent = ftt_cell_parent (cell);
  gdouble h[4];
  guint i;
  for (i = 0; i < 4; i++)
    h[i] = GFS_VALUE (parent, v[i]);

  FttVector p;
  ftt_cell_relative_pos (cell, &p);
  p.x *= 2.; p.y *= 2.;

  hP[0] = h[0] +  h[1]*p.x + (h[2] + h[3]*p.x)*p.y;
  hP[1] = (h[1] + h[3]*p.y)/2.;
  hP[2] = (h[2] + h[3]*p.x)/2.;
  hP[3] = h[3]/4.;
}

static void projection_matrix (Polygon * poly, double m[2][2])
{
  double x[4], y[4];
  guint i;
  for (i = 0; i < 4; i++) {
    x[i] = (poly->p[i].x - poly->c.x)/poly->h;
    y[i] = (poly->p[i].y - poly->c.y)/poly->h;
  }
  m[0][0] = (x[0] - x[1] + x[3] - x[2])/4.;
  m[0][1] = (x[0] + x[1] - x[3] - x[2])/4.;
  m[1][0] = (y[0] - y[1] + y[3] - y[2])/4.;
  m[1][1] = (y[0] + y[1] - y[3] - y[2])/4.;
}

static void add_weighted_kdt_sum  (KdtSum * s, const KdtSum * stmp, double w)
{
  s->m01 += w*stmp->m01;
  s->m02 += w*stmp->m02;
  s->m03 += w*stmp->m03;

  s->m11 += w*stmp->m11;
  s->m13 += w*stmp->m13;

  s->m22 += w*stmp->m22;
  s->m23 += w*stmp->m23;
  s->m33 += w*stmp->m33;

  s->m44 += w*stmp->m44;
  s->m55 += w*stmp->m55;
  s->m66 += w*stmp->m66;
  s->m77 += w*stmp->m77;

  s->m67 += w*stmp->m67;
  s->m76 += w*stmp->m76;

  s->H0 += w*stmp->H0;
  s->H1 += w*stmp->H1;
  s->H2 += w*stmp->H2;
  s->H3 += w*stmp->H3;
  s->H4 += w*stmp->H4;
  s->H5 += w*stmp->H5;
  s->H6 += w*stmp->H6;

  if (stmp->Hmax > s->Hmax)
    s->Hmax = stmp->Hmax;
  
  if (stmp->Hmin < s->Hmin)
    s->Hmin = stmp->Hmin;
  s->n += stmp->n;
  s->w += w*stmp->w;
  s->coverage += stmp->coverage;
}

static void update_terrain_rms (GfsRefineTerrain * t, Polygon * poly, gboolean relative, RMS * rms)
{
  rms_init (t, rms, poly, relative);
  KdtSum s;
  kdt_sum_init (&s);
  guint i;
  KdtRect rect;
  rect[0].l = poly->min[0]; rect[0].h = poly->max[0];
  rect[1].l = poly->min[1]; rect[1].h = poly->max[1];
  for (i = 0; i < poly->rs->nrs; i++) {
    KdtSum stmp;
    kdt_sum_init (&stmp);
    kdt_query_sum (poly->rs->rs[i],
		   (KdtCheck) polygon_includes,
		   (KdtCheck) polygon_intersects, poly, 
		   rect, &stmp);
    add_weighted_kdt_sum (&s, &stmp, poly->rs->weight[i]);
  }
  
  rms->m[0][0] = s.w;
  rms->n = s.n;
  if (s.w > 0.) {
    KdtSum sp;

    sp.H0 = s.H0;
    sp.H4 = s.H4;
    sp.Hmin = s.Hmin;
    sp.Hmax = s.Hmax;

    /* The sums returned by kdt_query_region_sum are defined in
       a (lon,lat) coordinate system, we need to project these into
       the local Cartesian coordinate system. The corresponding
       transform is given by matrix p below. */
    double p[2][2];
    projection_matrix (poly, p);

    /* This is the transformation of the sums */
    double det = p[0][0]*p[1][1] - p[0][1]*p[1][0], det1 = det;
    g_assert (det > 0.1);

    sp.m01 = (s.m01*p[1][1] - s.m02*p[0][1])/det;
    sp.m02 = (s.m02*p[0][0] - s.m01*p[1][0])/det;
    sp.H1 =  (s.H1*p[1][1] - s.H2*p[0][1])/det;
    sp.H2 =  (s.H2*p[0][0] - s.H1*p[1][0])/det;

    det *= det1;
    sp.m11 = (p[1][1]*p[1][1]*s.m11 + p[0][1]*p[0][1]*s.m22 - 2.*p[0][1]*p[1][1]*s.m03)/det;
    sp.m22 = (p[1][0]*p[1][0]*s.m11 + p[0][0]*p[0][0]*s.m22 - 2.*p[0][0]*p[1][0]*s.m03)/det;
    sp.m03 = - (p[1][0]*p[1][1]*s.m11 + p[0][0]*p[0][1]*s.m22 
		- (p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.m03)/det;
    sp.H3  = - (p[1][0]*p[1][1]*s.H5 + p[0][0]*p[0][1]*s.H6 
		- (p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.H3)/det;

    det *= det1;
    sp.m13 = (- p[1][0]*p[1][1]*p[1][1]*s.m44 
	      + p[0][0]*p[0][1]*p[0][1]*s.m55
	      + p[1][1]*(p[0][0]*p[1][1] + 2.*p[0][1]*p[1][0])*s.m13
	      - p[0][1]*(2.*p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.m23)/det;
    sp.m23 = (+ p[1][0]*p[1][0]*p[1][1]*s.m44 
	      - p[0][0]*p[0][0]*p[0][1]*s.m55
	      - p[1][0]*(2.*p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.m13
	      + p[0][0]*(p[0][0]*p[1][1] + 2.*p[0][1]*p[1][0])*s.m23)/det;

    det *= det1;
    sp.m33 = (+ p[1][0]*p[1][0]*p[1][1]*p[1][1]*s.m66
	      + p[0][0]*p[0][0]*p[0][1]*p[0][1]*s.m77
	      - 2.*p[1][0]*p[1][1]*(p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.m67 
	      - 2.*p[0][0]*p[0][1]*(p[0][0]*p[1][1] + p[0][1]*p[1][0])*s.m76
	      + (p[0][0]*p[0][0]*p[1][1]*p[1][1]
		 + 4.*p[0][0]*p[0][1]*p[1][0]*p[1][1]
		 + p[0][1]*p[0][1]*p[1][0]*p[1][0])*s.m33)/det;

    rms->m[0][1] = sp.m01;
    rms->m[0][2] = sp.m02;
    rms->m[1][1] = sp.m11;
    rms->m[2][2] = sp.m22;
    rms->m[0][3] = sp.m03;
    rms->m[1][3] = sp.m13;
    rms->m[2][3] = sp.m23;
    rms->m[3][3] = sp.m33;
    rms->m[1][2] = rms->m[0][3];

    if (rms->relative) {
      double hp[NM];

      parent_cell_coefficients (rms->cell, rms->t->h, hp);
      rms->H[0] = sp.H0 - s.w*hp[0] - sp.m01*hp[1] - sp.m02*hp[2] - sp.m03*hp[3];

      /* See terrain.mac for a "maxima" derivation of the terms below */
      rms->H[1] = sp.H1 - hp[0]*sp.m01 - hp[1]*sp.m11 - hp[2]*sp.m03 - hp[3]*sp.m13;
      rms->H[2] = sp.H2 - hp[0]*sp.m02 - hp[1]*sp.m03 - hp[2]*sp.m22 - hp[3]*sp.m23;
      rms->H[3] = sp.H3 - hp[0]*sp.m03 - hp[1]*sp.m13 - hp[2]*sp.m23 - hp[3]*sp.m33;
      rms->H[4] = (sp.H4 - 2.*hp[3]*sp.H3 - 2.*hp[2]*sp.H2 - 2.*hp[1]*sp.H1 - 2.*hp[0]*sp.H0
		   + hp[3]*hp[3]*sp.m33
		   + 2.*hp[2]*hp[3]*sp.m23
		   + hp[2]*hp[2]*sp.m22
		   + 2.*hp[1]*hp[3]*sp.m13
		   + hp[1]*hp[1]*sp.m11
		   + 2.*(hp[0]*hp[3] + hp[1]*hp[2])*sp.m03
		   + 2.*hp[0]*hp[2]*sp.m02
		   + 2.*hp[0]*hp[1]*sp.m01
		   + hp[0]*hp[0]*s.w);
    }
    else {
      rms->H[0] = sp.H0;
      rms->H[1] = sp.H1;
      rms->H[2] = sp.H2;
      rms->H[3] = sp.H3;
      rms->H[4] = sp.H4;
    }
    rms->max = sp.Hmax;
    rms->min = sp.Hmin;
  }
}

static void update_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  RMS rms;
  guint i;
  g_assert (GFS_VALUE (cell, t->type) == REFINED);
  Polygon poly;
  polygon_init (gfs_object_simulation (t), &poly, cell, &t->rs);
  update_terrain_rms (t, &poly, ftt_cell_parent (cell) != NULL, &rms);
  rms_update (&rms);

  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = rms.h[i];
  GFS_VALUE (cell, t->he) = rms.he;
  GFS_VALUE (cell, t->hn) = rms.n;
  GFS_VALUE (cell, t->hdmin) = rms.min <   G_MAXDOUBLE ? rms.min : GFS_NODATA;
  GFS_VALUE (cell, t->hdmax) = rms.max > - G_MAXDOUBLE ? rms.max : GFS_NODATA;
  GFS_VALUE (cell, t->type) = RAW;
}

static void function_from_parent (FttCell * cell, GfsRefineTerrain * t, gdouble h[4])
{
  gdouble H[4];
  corners_from_parent (cell, t, H);
  function_from_corners (h, H);
}

static void cell_fine_init (FttCell * parent, GfsRefineTerrain * t)
{
  gfs_cell_fine_init (parent, GFS_DOMAIN (gfs_object_simulation (t)));
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  guint i;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      gdouble h[NM];
      function_from_parent (child.c[i], t, h);
      guint j;
      for (j = 0; j < NM; j++)
	GFS_VALUE (child.c[i], t->h[j]) = h[j];
      GFS_VALUE (child.c[i], t->he) = GFS_VALUE (parent, t->he);
      GFS_VALUE (child.c[i], t->hn) = GFS_VALUE (parent, t->hn)/FTT_CELLS;
      GFS_VALUE (child.c[i], t->hdmin) = GFS_VALUE (parent, t->hdmin);
      GFS_VALUE (child.c[i], t->hdmax) = GFS_VALUE (parent, t->hdmax);
      GFS_VALUE (child.c[i], t->type) = NEW_CHILD;
    }
}

static gdouble corner_value (GfsRefineTerrain * t, FttVector * p, gdouble eps, guint level)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
  gdouble v = 0., w = 0.;
  gint i, j;
  for (i = -1; i <= 1; i += 2)
    for (j = -1; j <= 1; j += 2) {
      FttVector q;
      q.x = p->x + eps*i; q.y = p->y + eps*j; q.z = p->z;
      FttCell * cell = gfs_domain_locate (domain, q, level, NULL);
      if (cell) {
	if (ftt_cell_level (cell) < level)
	    return 0.;
	else if (GFS_VALUE (cell, t->type) == FAIR)
	  return cell_value (cell, t->h, *p);
	gdouble n = GFS_VALUE (cell, t->hn);
	if (n > 0.) {
	  g_assert (GFS_VALUE (cell, t->type) == RAW);
	  v += cell_value (cell, t->h, *p);
	  w += 1.;
	}
      }
    }
  return w > 0 ? v/w : 0.;
}

static void update_error_estimate (FttCell * cell, GfsRefineTerrain * t, gboolean relative)
{
  if (GFS_VALUE (cell, t->hn) > 0.) {
    RMS rms;
    guint i;
    Polygon poly;
    polygon_init (gfs_object_simulation (t), &poly, cell, &t->rs);
    update_terrain_rms (t, &poly, relative, &rms);
    for (i = 0; i < NM; i++)
      rms.h[i] = GFS_VALUE (cell, t->h[i]);
    GFS_VALUE (cell, t->he) = rms_minimum (&rms);
  }
  else
    GFS_VALUE (cell, t->he) = 0.;
}

static void remove_knots (FttCell * cell, GfsRefineTerrain * t)
{
  gdouble size = ftt_cell_size (cell), eps = size/1000.;
  guint level = ftt_cell_level (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h[4], H[4];
  p.x += size/2.; p.y += size/2.;
  H[0] = corner_value (t, &p, eps, level);
  p.x -= size;
  H[1] = corner_value (t, &p, eps, level);
  p.y -= size;
  H[2] = corner_value (t, &p, eps, level);
  p.x += size;
  H[3] = corner_value (t, &p, eps, level);
  function_from_corners (h, H);
  guint i;
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = h[i];
  GFS_VALUE (cell, t->type) = FAIR;

  update_error_estimate (cell, t, ftt_cell_parent (cell) != NULL);
}

static void update_height_and_check_for_refinement (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) == FAIR) {
    if (ftt_cell_parent (cell)) {
      gdouble h[4];
      function_from_parent (cell, t, h);
      guint i;
      for (i = 0; i < NM; i++)
	GFS_VALUE (cell, t->h[i]) += h[i];
    }

    if (ftt_cell_level (cell) < gfs_function_value (GFS_REFINE (t)->maxlevel, cell) &&
	gfs_function_value (t->criterion, cell)) {
      g_assert (FTT_CELL_IS_LEAF (cell));
      ftt_cell_refine_single (cell, (FttCellInitFunc) cell_fine_init, t);
      FttCellChildren child;
      guint i;
      ftt_cell_children (cell, &child);
      for (i = 0; i < FTT_CELLS; i++)
	GFS_VALUE (child.c[i], t->type) = REFINED;
    }

    if (!FTT_CELL_IS_LEAF (cell))
      t->refined = TRUE;
  }
  else
    g_assert (GFS_VALUE (cell, t->type) == NEW_CHILD);
}

static void reset_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  guint i;
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = 0.;
  GFS_VALUE (cell, t->type) = REFINED;
  if (FTT_CELL_IS_LEAF (cell) && ftt_cell_level (cell) < t->level)
    t->level = ftt_cell_level (cell);
}

#if FTT_2D
# define traverse_boundary(domain,order,flags,depth,func,data) \
         gfs_domain_cell_traverse(domain,order,flags,depth,func,data)
#else /* 3D */
# define traverse_boundary(domain,order,flags,depth,func,data) \
         gfs_domain_cell_traverse_boundary(domain,FTT_FRONT,order,flags,depth,func,data)

static void terrain_min_max (gdouble H[NM], gdouble minmax[2], gdouble scale)
{
  gdouble dx, dy;
  minmax[0] = G_MAXDOUBLE; minmax[1] = - G_MAXDOUBLE;
  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      gdouble v = H[0] + dx*H[1] + dy*H[2] + dx*dy*H[3];
      if (v < minmax[0]) minmax[0] = v;
      if (v > minmax[1]) minmax[1] = v;
    }
  minmax[0] *= scale;
  minmax[1] *= scale;
}

static void min_max (FttCell * cell, GfsRefineTerrain * t)
{
  gdouble minmax[2] = { G_MAXDOUBLE, - G_MAXDOUBLE };
  if (FTT_CELL_IS_LEAF (cell)) {
    gdouble h[4];
    h[0] = GFS_VALUE (cell, t->h[0]);
    h[1] = GFS_VALUE (cell, t->h[1]);
    h[2] = GFS_VALUE (cell, t->h[2]);
    h[3] = GFS_VALUE (cell, t->h[3]);
    terrain_min_max (h, minmax, t->scale);

    FttVector p;
    ftt_cell_pos (cell, &p);
    if (p.z > t->front)
      t->front = p.z;
  }
  else {
    FttCellChildren child;
    guint i, n = ftt_cell_children_direction (cell, FTT_FRONT, &child);
    for (i = 0; i < n; i++)
      if (child.c[i]) {
	if (GFS_VALUE (child.c[i], t->max) > minmax[1])
	  minmax[1] = GFS_VALUE (child.c[i], t->max);
	if (GFS_VALUE (child.c[i], t->min) < minmax[0])
	  minmax[0] = GFS_VALUE (child.c[i], t->min);	
      }
  }
  GFS_VALUE (cell, t->min) = minmax[0];
  GFS_VALUE (cell, t->max) = minmax[1];
  GFS_VALUE (cell, t->type) = BOUNDARY;
}

static gboolean refine_terrain_from_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h = ftt_cell_size (cell)/2., zmin = p.z - h, zmax = p.z + h;
  p.z = t->front;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
  FttCell * boundary = gfs_domain_locate (domain, p, ftt_cell_level (cell), NULL);
  g_assert (boundary);
  if (GFS_VALUE (boundary, t->min) > zmax || GFS_VALUE (boundary, t->max) < zmin)
    return FALSE;
  GFS_VALUE (cell, t->type) = CONTAINS_SURFACE;
  return !FTT_CELL_IS_LEAF (boundary);
}

static void refine_box (GfsBox * box, GfsRefineTerrain * t)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_terrain_from_boundary, t,
		   (FttCellInitFunc) gfs_cell_fine_init, gfs_box_domain (box));
}

static void init_terrain_from_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) == CONTAINS_SURFACE) {
    FttVector p;
    ftt_cell_pos (cell, &p);
    p.z = t->front;
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
    FttCell * boundary = gfs_domain_locate (domain, p, -1, NULL);
    g_assert (boundary);
    g_assert (ftt_cell_level (cell) == ftt_cell_level (boundary));
    guint i;
    for (i = 0; i < NM; i++)
      GFS_VALUE (cell, t->h[i]) = GFS_VALUE (boundary, t->h[i]);
    GFS_VALUE (cell, t->he) = GFS_VALUE (boundary, t->he);
    GFS_VALUE (cell, t->hn) = GFS_VALUE (boundary, t->hn);
  }
}

static gboolean coarsen_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  return (GFS_VALUE (cell, t->type) != CONTAINS_SURFACE);
}

static void coarsen_box (GfsBox * box, GfsRefineTerrain * t)
{
  ftt_cell_coarsen (box->root,
		    (FttCellCoarsenFunc) coarsen_boundary, t,
		    (FttCellCleanupFunc) gfs_cell_cleanup, gfs_box_domain (box));
}

static void reset_empty_cell (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) != CONTAINS_SURFACE) {
    guint i;
    for (i = 0; i < NM; i++)
      GFS_VALUE (cell, t->h[i]) = GFS_NODATA;
    GFS_VALUE (cell, t->he) = 0.;
    GFS_VALUE (cell, t->hn) = 0.;
  }
}
#endif /* 3D */

#if DEBUG
static void draw_terrain (FttCell * cell, gpointer * data)
{
  GfsRefineTerrain * t = data[0];
  FILE * fp = data[1];
  gdouble h = ftt_cell_size (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += h/2.; p.y += h/2.;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.x -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.y -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.x += h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
}

static void draw_level (GfsDomain * domain, GfsRefine * refine, guint level, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  fprintf (data[1], "QUAD\n");
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, level,
		     (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}

static void draw_all (GfsDomain * domain, GfsRefine * refine, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  //  fprintf (data[1], "QUAD\n");
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}
#endif

#define ASCII_ZERO 48 /* ASCII value for character "0" */

static void terrain_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint len = strlen (v->name) - 1;
  gint c = v->name[len] - ASCII_ZERO;
  guint n;
  gdouble h[NM];

  g_assert (c >= 0 && c < NM);
  for (n = 0; n < NM; n++) {
    GSList * i = v->domain->variables;
    while (i && (!GFS_VARIABLE (i->data)->name || 
		 strncmp (v->name, GFS_VARIABLE (i->data)->name, len) ||
		 GFS_VARIABLE (i->data)->name[len] != ASCII_ZERO + n))
      i = i->next;
    g_assert (i);
    h[n] = GFS_VALUE (parent, GFS_VARIABLE (i->data));
  }

  ftt_cell_children (parent, &child);
  if (h[0] == GFS_NODATA) {
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n])
	GFS_VALUE (child.c[n], v) = GFS_NODATA;
  }
  else {
#if !FTT_2D
    gdouble size = ftt_cell_size (parent)/4.;
#endif
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n]) {
	gdouble hc[NM];
	FttVector p;
	ftt_cell_relative_pos (child.c[n], &p);
	p.x *= 2.; p.y *= 2.;
	hc[0] = h[0] + h[1]*p.x + h[2]*p.y + h[3]*p.x*p.y;
	hc[1] = (h[1] + h[3]*p.y)/2.;
	hc[2] = (h[2] + h[3]*p.x)/2.;
	hc[3] = h[3]/4.;
#if !FTT_2D
	ftt_cell_pos (child.c[n], &p);
	gdouble zmin = p.z - size, zmax = p.z + size, minmax[2];
	p.z = 1.;
	gfs_simulation_map (GFS_SIMULATION (v->domain), &p);
	terrain_min_max (hc, minmax, p.z);
	if (minmax[0] > zmax || minmax[1] < zmin)
	  GFS_VALUE (child.c[n], v) = GFS_NODATA;
	else
#endif
	  GFS_VALUE (child.c[n], v) = hc[c];
      }
  }
}

static void hn_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  guint i, n = 0;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      n++;
  if (n > 0) {
    guint hn = GFS_VALUE (parent, v)/n;
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	GFS_VALUE (child.c[i], v) = hn;
  }
}

static void terrain_refine (GfsRefine * refine, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (refine);
  t->type = gfs_temporary_variable (domain);
  t->level = G_MAXINT/2;
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) reset_terrain, refine);
  do {
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) update_terrain, refine);
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) remove_knots, refine);
    t->refined = FALSE;
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) update_height_and_check_for_refinement,
		       refine);
#if DEBUG
    GfsNorm norm = gfs_domain_norm_variable (domain, t->he, NULL, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "level: %d bias: %g 1: %g 2: %g inf: %g\n", 
	     t->level, norm.bias, norm.first, norm.second, norm.infty);
    fprintf (stderr, "level: %d depth: %d\n", t->level, gfs_domain_depth (domain));
    gchar name[] = "/tmp/level-x";
    name[11] = ASCII_ZERO + t->level;
    draw_level (domain, refine, t->level, name);
#endif
    t->level++;
  } while (t->refined);
#if DEBUG
  draw_all (domain, refine, "/tmp/all");
#endif
#if !FTT_2D
  /* The height field is only defined on the front boundary, we need
     to define it volumetrically */
  t->min = gfs_temporary_variable (domain);
  t->max = gfs_temporary_variable (domain);
  t->front = - G_MAXDOUBLE;
  FttVector p = {0.,0.,1.};
  gfs_simulation_map (sim, &p);
  t->scale = p.z;
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT, FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
				     (FttCellTraverseFunc) min_max, t);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) refine_box, t);
  gts_object_destroy (GTS_OBJECT (t->min));
  gts_object_destroy (GTS_OBJECT (t->max));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) init_terrain_from_boundary, t);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) coarsen_box, t);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_empty_cell, t);
#endif /* 3D */  
  gts_object_destroy (GTS_OBJECT (t->type));
  guint i;
  for (i = 0; i < NM; i++)
    t->h[i]->coarse_fine = terrain_coarse_fine;
}

static void refine_terrain_destroy (GtsObject * object)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (object);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (object));

  if (t->name) {
    gchar * dname = g_strconcat (t->name, "min", NULL);
    gfs_domain_remove_derived_variable (domain, dname);
    g_free (dname);
    
    dname = g_strconcat (t->name, "max", NULL);
    gfs_domain_remove_derived_variable (domain, dname);
    g_free (dname);
  }
  g_free (t->name);

  kdtrees_destroy (&t->rs);
  
  gts_object_destroy (GTS_OBJECT (t->criterion));  
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->destroy) (object);
}

static gdouble terrain_hmin (FttCell * cell, FttCellFace * face, 
			     GfsDomain * domain, GfsRefineTerrain * t)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble dx, dy, min = G_MAXDOUBLE;
  gdouble H0 = GFS_VALUE (cell, t->h[0]), H1 = GFS_VALUE (cell, t->h[1]);
  gdouble H2 = GFS_VALUE (cell, t->h[2]), H3 = GFS_VALUE (cell, t->h[3]);

  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      double v = H0 + dx*H1 + dy*H2 + dx*dy*H3;
      if (v < min) min = v;
    }
  return min;
}

static gdouble terrain_hmax (FttCell * cell, FttCellFace * face, 
			     GfsDomain * domain, GfsRefineTerrain * t)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble dx, dy, max = - G_MAXDOUBLE;
  gdouble H0 = GFS_VALUE (cell, t->h[0]), H1 = GFS_VALUE (cell, t->h[1]);
  gdouble H2 = GFS_VALUE (cell, t->h[2]), H3 = GFS_VALUE (cell, t->h[3]);

  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      double v = H0 + dx*H1 + dy*H2 + dx*dy*H3;
      if (v > max) max = v;
    }
  return max;
}

static void none (FttCell * parent, GfsVariable * v)
{
}

static void refine_terrain_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (*o);
  t->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  kdtrees_read (&t->rs, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  guint i;
  for (i = 0; i < NM; i++) {
    gchar * name = g_strdup_printf ("%s%d", t->name, i);
    t->h[i] = gfs_domain_get_or_add_variable (domain, name, "Terrain height");
    t->h[i]->coarse_fine = none;
    g_free (name);
  }
  gchar * name = g_strjoin (NULL, t->name, "e", NULL);
  t->he = gfs_domain_get_or_add_variable (domain, name, "Terrain RMS error");
  t->he->coarse_fine = none;
  g_free (name);
  name = g_strjoin (NULL, t->name, "n", NULL);
  t->hn = gfs_domain_get_or_add_variable (domain, name, "Terrain samples #");
  t->hn->coarse_fine = hn_coarse_fine;
  g_free (name);
  name = g_strjoin (NULL, t->name, "dmin", NULL);
  t->hdmin = gfs_domain_get_or_add_variable (domain, name, "Minimum data height");
  t->hdmin->coarse_fine = none;
  g_free (name);
  name = g_strjoin (NULL, t->name, "dmax", NULL);
  t->hdmax = gfs_domain_get_or_add_variable (domain, name, "Maximum data height");
  t->hdmax->coarse_fine = none;
  g_free (name);

  GfsDerivedVariableInfo v;

  v.name = g_strjoin (NULL, t->name, "min", NULL);
  v.description = "Minimum terrain height";
  v.func = terrain_hmin;
  v.data = t;
  if (!gfs_domain_add_derived_variable (domain, v)) {
    gts_file_error (fp, "derived variable `%s' already defined", v.name);
    g_free (v.name);
    return;
  }
  g_free (v.name);

  v.name = g_strjoin (NULL, t->name, "max", NULL);
  v.description = "Maximum terrain height";
  v.func = terrain_hmax;
  v.data = t;
  if (!gfs_domain_add_derived_variable (domain, v)) {
    gts_file_error (fp, "derived variable `%s' already defined", v.name);
    g_free (v.name);
    return;
  }
  g_free (v.name);

  gfs_function_read (t->criterion, domain, fp);
}

static void refine_terrain_write (GtsObject * o, FILE * fp)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (o);
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", t->name);
  kdtrees_write (&t->rs, fp);
  gfs_function_write (t->criterion, fp);
}

static void gfs_refine_terrain_class_init (GfsRefineClass * klass)
{
  klass->refine = terrain_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_terrain_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_terrain_read;
  GTS_OBJECT_CLASS (klass)->write = refine_terrain_write;
}

static void gfs_refine_terrain_init (GfsRefineTerrain * t)
{
  t->criterion = gfs_function_new (gfs_function_class (), 0.);
  t->rs.basename = g_strdup ("*");
}

GfsRefineClass * gfs_refine_terrain_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_terrain_info = {
      "GfsRefineTerrain",
      sizeof (GfsRefineTerrain),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_terrain_class_init,
      (GtsObjectInitFunc) gfs_refine_terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_terrain_info);
  }

  return klass;
}

/* GfsSurfaceTerrain: Header */

typedef struct _GfsSurfaceTerrain         GfsSurfaceTerrain;

struct _GfsSurfaceTerrain {
  /*< private >*/
  GfsGenericSurface parent;
  GfsVariable * h[NM];
  gdouble scale;

  /*< public >*/
  gchar * name;
};

#define GFS_SURFACE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceTerrain,\
					         gfs_surface_terrain_class ())
#define GFS_IS_SURFACE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_terrain_class ()))

GfsGenericSurfaceClass * gfs_surface_terrain_class  (void);

/* GfsSurfaceTerrain: Object */

static void gfs_surface_terrain_read (GtsObject ** o, GtsFile * fp)
{
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsSurfaceTerrain * t = GFS_SURFACE_TERRAIN (*o);
  t->name = g_strdup (fp->token->str);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  guint i;
  for (i = 0; i < NM; i++) {
    gchar * name = g_strdup_printf ("%s%d", t->name, i);
    t->h[i] = gfs_variable_from_name (domain->variables, name);
    if (!t->h[i]) {
      gts_file_error (fp, "%s is not a valid variable name", name);
      g_free (name);
      return;
    }
    t->h[i]->coarse_fine = terrain_coarse_fine;
    g_free (name);
  }
  gts_file_next_token (fp);
}

static void gfs_surface_terrain_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, " %s", GFS_SURFACE_TERRAIN (o)->name);
}

static void gfs_surface_terrain_destroy (GtsObject * object)
{
  g_free (GFS_SURFACE_TERRAIN (object)->name);
  (* GTS_OBJECT_CLASS (gfs_surface_terrain_class ())->parent_class->destroy)
    (object);
}

static GfsGenericSurface * cell_is_cut (FttCell * cell, GfsGenericSurface * s1,
					gboolean flatten, gint maxlevel)
{
  g_assert (!flatten); /* not implemented */
  if (!FTT_CELL_IS_LEAF (cell))
    return s1;
  return GFS_VALUE (cell, GFS_SURFACE_TERRAIN (s1)->h[0]) != GFS_NODATA ? s1 : NULL;
}

static gdouble zscale (GfsSurfaceTerrain * t)
{
  if (t->scale == 0.) {
    FttVector p = {0.,0.,1.};
    gfs_simulation_map (gfs_object_simulation (t), &p);
    t->scale = p.z;
  }
  return t->scale;
}

static guint surface_segment_intersection (GfsGenericSurface * s1,
					   FttCell * cell,
					   GfsSegment * I)
{
  I->n = 0;
  I->x = 0.;
  I->inside = 0;

  FttVector pE, pD;
  pE.x = I->E->x; pE.y = I->E->y;
  pD.x = I->D->x; pD.y = I->D->y;
  GfsSurfaceTerrain * t = GFS_SURFACE_TERRAIN (s1);
  gdouble vE = I->E->z - cell_value (cell, t->h, pE)*zscale (t);
  gdouble vD = I->D->z - cell_value (cell, t->h, pD)*zscale (t);
  
  if ((vE > 0. && vD <= 0.) || (vE <= 0. && vD > 0.)) {
    I->n = 1;
    I->inside = vE > 0. ? -1 : 1;
    I->x = vE/(vE - vD);
#if DEBUG
    gdouble size = ftt_cell_size (cell)/2.;
    FttVector q;
    ftt_cell_pos (cell, &q);
    pE.x = (pE.x - q.x)/size;
    pE.y = (pE.y - q.y)/size;
    pD.x = (pD.x - q.x)/size;
    pD.y = (pD.y - q.y)/size;
    fprintf (stderr, "p %g %g %g %g %g %g %g %d %g %g %g %g\n", 
	     I->D->x, I->D->y, I->D->z,
	     I->E->x, I->E->y, I->E->z,
	     I->x,
	     ftt_cell_level (cell),
	     pE.x, pE.y, pD.x, pD.y);
    fprintf (stderr, "q %g %g %g\nq %g %g %g\nq\nq\n",
	     I->D->x, I->D->y, I->D->z,
	     I->E->x, I->E->y, I->E->z);
    fprintf (stderr, "i %g %g %g\n",
	     I->E->x + I->x*(I->D->x - I->E->x),
	     I->E->y + I->x*(I->D->y - I->E->y),
	     I->E->z + I->x*(I->D->z - I->E->z));
#endif
  }
  return I->n;
}

static void surface_segment_normal (GfsGenericSurface * s1,
				    FttCell * cell,
				    GfsSegment * I,
				    GtsVector n)
{
  GfsVariable ** h = GFS_SURFACE_TERRAIN (s1)->h;
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector p, q;
  ftt_cell_pos (cell, &q);
  p.x = I->E->x + I->x*(I->D->x - I->E->x);
  p.y = I->E->y + I->x*(I->D->y - I->E->y);
  p.x = (p.x - q.x)/size;
  p.y = (p.y - q.y)/size;
  n[0] = - (GFS_VALUE (cell, h[1]) + GFS_VALUE (cell, h[3])*p.y)/size;
  n[1] = - (GFS_VALUE (cell, h[2]) + GFS_VALUE (cell, h[3])*p.x)/size;
  n[2] = 1./zscale (GFS_SURFACE_TERRAIN (s1));
}

static void gfs_surface_terrain_class_init (GfsGenericSurfaceClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_terrain_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_terrain_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_terrain_destroy;

  klass->cell_is_cut = cell_is_cut;
  klass->segment_intersection = surface_segment_intersection;
  klass->segment_normal = surface_segment_normal;
}

GfsGenericSurfaceClass * gfs_surface_terrain_class (void)
{
  static GfsGenericSurfaceClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_terrain_info = {
      "GfsSurfaceTerrain",
      sizeof (GfsSurfaceTerrain),
      sizeof (GfsGenericSurfaceClass),
      (GtsObjectClassInitFunc) gfs_surface_terrain_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_surface_class ()),
				  &gfs_surface_terrain_info);
  }

  return klass;
}

/* GfsTerrain: Header */

#define GFS_IS_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						 gfs_terrain_class ()))

GfsEventClass * gfs_terrain_class  (void);

/* GfsTerrain: Object */

static void terrain_init (GfsSolid * s)
{
  gts_object_destroy (GTS_OBJECT (s->s));
  s->s = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_terrain_class ())));
}

GfsEventClass * gfs_terrain_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_terrain_info = {
      "GfsTerrain",
      sizeof (GfsSolid),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_class ()),
				  &gfs_terrain_info);
  }

  return klass;
}

/* GfsVariableTerrain: header */

typedef struct _GfsVariableTerrain                GfsVariableTerrain;

struct _GfsVariableTerrain {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * p, * H, * n, * dmin, * dmax;
  Kdtrees rs;
};

#define GFS_VARIABLE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableTerrain,\
					           gfs_variable_terrain_class ())
#define GFS_IS_VARIABLE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_terrain_class ()))

GfsVariableClass * gfs_variable_terrain_class  (void);

/* GfsVariableTerrain: Object */

static void variable_terrain_destroy (GtsObject * o)
{
  kdtrees_destroy (&GFS_VARIABLE_TERRAIN (o)->rs);

  (* GTS_OBJECT_CLASS (gfs_variable_terrain_class ())->parent_class->destroy) (o);
}

static double reconstruct_terrain (FttCell * cell, GfsVariableTerrain * t)
{
  GfsVariable * v = GFS_VARIABLE (t);
  GfsSimulation * sim = GFS_SIMULATION (v->domain);
  Polygon poly;
  KdtRect rect;
  KdtSum s;
  guint i;
  polygon_init (sim, &poly, cell, &t->rs);
  kdt_sum_init (&s);
  rect[0].l = poly.min[0]; rect[0].h = poly.max[0];
  rect[1].l = poly.min[1]; rect[1].h = poly.max[1]; 
  for (i = 0; i < poly.rs->nrs; i++) {
    KdtSum stmp;
    kdt_sum_init (&stmp);
    kdt_query_sum (poly.rs->rs[i],
		   (KdtCheck) polygon_includes,
		   (KdtCheck) polygon_intersects, &poly, 
		   rect, &stmp);
    add_weighted_kdt_sum (&s, &stmp, poly.rs->weight[i]);
  }
  GFS_VALUE (cell, t->n) = s.n;
  if (s.w > 0.) {
    GfsVariable * v = GFS_VARIABLE (t);
    GFS_VALUE (cell, v) = s.H0/s.w/sim->physical_params.L;
    GFS_VALUE (cell, t->dmin) = s.Hmin;
    GFS_VALUE (cell, t->dmax) = s.Hmax;
    return s.w;
  }
  return 0.;
}

static void variable_terrain_coarse_fine (FttCell * parent, GfsVariable * v)
{
  GfsVariableTerrain * t = GFS_VARIABLE_TERRAIN (v);
  GfsSimulation * sim = GFS_SIMULATION (v->domain);
  FttCellChildren child;
  gdouble f[4*(FTT_DIMENSION - 1) + 1] = { GFS_NODATA };
  guint n;

  /* Reconstruct terrain */
  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n] && !reconstruct_terrain (child.c[n], t)) {
      GFS_VALUE (child.c[n], t->dmin) = GFS_NODATA;
      GFS_VALUE (child.c[n], t->dmax) = GFS_NODATA;
      if (GFS_CELL_IS_BOUNDARY (parent))
	GFS_VALUE (child.c[n], v) = GFS_VALUE (parent, v);
      else {
	FttVector p;
	ftt_cell_pos (child.c[n], &p);
	if (f[0] == GFS_NODATA)
	  gfs_cell_corner_values (parent, v, ftt_cell_level (parent), f);
	GFS_VALUE (child.c[n], v) = gfs_interpolate_from_corners (parent, p, f);
      }
    }

  /* If we are part of GfsRiver, reconstruct H and P */
  if (t->H) {
    /* Reconstruct H */
    double dry = GFS_RIVER (sim)->dry;
    if (GFS_VALUE (parent, t->p) < dry) {
      /* Dry cell */
      FttCellNeighbors neighbor;
      ftt_cell_neighbors (parent, &neighbor);
      FttDirection d;
      gdouble H = 0., s = 0.;

      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (neighbor.c[d]) {
	  if (FTT_CELL_IS_LEAF (neighbor.c[d])) {
	    if (GFS_VALUE (neighbor.c[d], t->p) >= dry) {
	       H += GFS_VALUE (neighbor.c[d], t->p)*GFS_VALUE (neighbor.c[d], t->H);
	       s += GFS_VALUE (neighbor.c[d], t->p);
	    }
	  }
	  else {
	    FttCellChildren child;
	    guint i, n = ftt_cell_children_direction (neighbor.c[d],
						      FTT_OPPOSITE_DIRECTION (d), &child);
	    for (i = 0; i < n; i++)
	      if (child.c[i] && GFS_VALUE (child.c[i], t->p) >= dry) {
		H += GFS_VALUE (child.c[i], t->p)*GFS_VALUE (child.c[i], t->H);
		s += GFS_VALUE (child.c[i], t->p);
	      }
	  }
	}
      
      if (s > 0.) {
	H /= s; /* average H of neighbouring wet cells */
	for (n = 0; n < FTT_CELLS; n++)
	  if (child.c[n])
	    GFS_VALUE (child.c[n], t->H) = H;
      }
      else { /* surrounded by dry cells */
	for (n = 0; n < FTT_CELLS; n++)
	  if (child.c[n])
	    GFS_VALUE (child.c[n], t->H) = 0.; /* default "sealevel" */
      }
    }
    else {
      /* wet cell */
      GfsVariable * v = t->H;
      for (n = 0; n < FTT_CELLS; n++)
	if (child.c[n])
	  GFS_VALUE (child.c[n], v) = GFS_VALUE (parent, v);
      
      if (!GFS_CELL_IS_BOUNDARY (parent)) {
	FttVector g;
	FttComponent c;
	FttCellNeighbors neighbor;
	ftt_cell_neighbors (parent, &neighbor);

	for (c = 0; c < FTT_DIMENSION; c++)
	  if (neighbor.c[2*c] && GFS_VALUE (neighbor.c[2*c], t->p) >= dry &&
	      neighbor.c[2*c + 1] && GFS_VALUE (neighbor.c[2*c + 1], t->p) >= dry)
	    (&g.x)[c] = gfs_center_minmod_gradient (parent, c, v->i);
	  else
	    (&g.x)[c] = 0.;
	
	for (n = 0; n < FTT_CELLS; n++) 
	  if (child.c[n]) {
	    FttVector p;
	    
	    ftt_cell_relative_pos (child.c[n], &p);
	    for (c = 0; c < FTT_DIMENSION; c++)
	      GFS_VALUE (child.c[n], v) += (&p.x)[c]*(&g.x)[c];
	  }
      }
    }
    /* Deduce P from the reconstruction of Zb and H */
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n]) {
	GFS_VALUE (child.c[n], t->p) = MAX (0., (GFS_VALUE (child.c[n], t->H) - 
						 GFS_VALUE (child.c[n], v)));
	GFS_VALUE (child.c[n], t->H) = GFS_VALUE (child.c[n], t->p) + GFS_VALUE (child.c[n], v);
      }
  }
}

static void variable_terrain_fine_coarse (FttCell * parent, GfsVariable * v)
{
  GfsVariableTerrain * t = GFS_VARIABLE_TERRAIN (v);
  FttCellChildren child;
  guint n;

  /* Reconstruct terrain (weighted average) */
  gdouble Zb = 0., sa = 0.;
  gdouble N = 0., dmin = G_MAXDOUBLE, dmax = - G_MAXDOUBLE;
  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++) 
    if (child.c[n]) {
      gdouble a = GFS_IS_MIXED (child.c[n]) ? GFS_STATE (child.c[n])->solid->a : 1.;
      Zb += GFS_VALUE (child.c[n], v)*a;
      sa += a;
      N += GFS_VALUE (child.c[n], t->n);
      if (GFS_VALUE (child.c[n], t->n) > 0) {
	dmax = MAX (dmax, GFS_VALUE (child.c[n], t->dmax));
	dmin = MIN (dmin, GFS_VALUE (child.c[n], t->dmin));
      }
    }
  if (sa > 0.)
    GFS_VALUE (parent, v) = Zb/sa;
  GFS_VALUE (parent, t->n) = N;
  GFS_VALUE (parent, t->dmax) = dmax > - G_MAXDOUBLE ? dmax : GFS_NODATA;
  GFS_VALUE (parent, t->dmin) = dmin <   G_MAXDOUBLE ? dmin : GFS_NODATA;

  /* If we are part of GfsRiver, reconstruct H and P */
  if (t->H) {
    /* Reconstruct H */
    gdouble H = 0., sa = 0.;
    for (n = 0; n < FTT_CELLS; n++) 
      if (child.c[n] && GFS_VALUE (child.c[n], t->p) >= GFS_RIVER (v->domain)->dry) {
	gdouble a = GFS_IS_MIXED (child.c[n]) ? GFS_STATE (child.c[n])->solid->a : 1.;
	H += GFS_VALUE (child.c[n], t->H)*a;
	sa += a;
      }
    if (sa > 0.) {
      GFS_VALUE (parent, t->H) = H/sa;
      GFS_VALUE (parent, t->p) = MAX (0., GFS_VALUE (parent, t->H) - GFS_VALUE (parent, v));
    }
    else { /* dry cell */
      GFS_VALUE (parent, t->p) = 0.;
      GFS_VALUE (parent, t->H) = GFS_VALUE (parent, v);
    }
  }
}

static void variable_terrain_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_terrain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariableTerrain * v = GFS_VARIABLE_TERRAIN (*o);
  kdtrees_read (&v->rs, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariable * v1 = GFS_VARIABLE (*o);
  v1->units = 1.;
  g_free (v1->description);
  v1->description = g_strdup ("Terrain");
  v1->coarse_fine = variable_terrain_coarse_fine;
  v1->fine_coarse = variable_terrain_fine_coarse;

  GfsSimulation * sim = gfs_object_simulation (*o);
  gchar * name = g_strjoin (NULL, v1->name, "n", NULL);
  v->n = gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), name, "Terrain samples # (weighted)");
  v->n->coarse_fine = none;
  v->n->fine_coarse = none;
  g_free (name);
  name = g_strjoin (NULL, v1->name, "dmin", NULL);
  v->dmin = gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), name, "Minimum data height");
  v->dmin->coarse_fine = none;
  v->dmin->fine_coarse = none;
  g_free (name);
  name = g_strjoin (NULL, v1->name, "dmax", NULL);
  v->dmax = gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), name, "Maximum data height");
  v->dmax->coarse_fine = none;
  v->dmax->fine_coarse = none;
  g_free (name);

  if (GFS_IS_RIVER (sim) && fp->type == '{') {
    gboolean reconstruct = FALSE;
    GtsFileVariable var[] = {
      {GTS_INT, "reconstruct", TRUE, &reconstruct},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (reconstruct) {
      v->p = GFS_RIVER (sim)->v[0];
      v->H = GFS_RIVER (sim)->h;
      /* the coarse -> fine and fine -> coarse interpolations of p and H
	 are taken over by variable_terrain_coarse_fine (below )*/
      v->p->coarse_fine = none;
      v->H->coarse_fine = none;
      v->p->fine_coarse = none;
      v->H->fine_coarse = none;
    }
  }
}

static void variable_terrain_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_terrain_class ())->parent_class->write) (o, fp);

  kdtrees_write (&GFS_VARIABLE_TERRAIN (o)->rs, fp);
  if (GFS_VARIABLE_TERRAIN (o)->H)
    fputs (" { reconstruct = 1 }", fp);
}

static void variable_terrain_class_init (GtsObjectClass * klass)
{
  klass->destroy = variable_terrain_destroy;
  klass->read = variable_terrain_read;
  klass->write = variable_terrain_write;
}

static void variable_terrain_init (GfsVariableTerrain * v)
{
  v->rs.basename = g_strdup ("*");
}

GfsVariableClass * gfs_variable_terrain_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_terrain_info = {
      "GfsVariableTerrain",
      sizeof (GfsVariableTerrain),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_terrain_class_init,
      (GtsObjectInitFunc) variable_terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_terrain_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "terrain";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gchar * path = getenv ("GFS_TERRAIN_PATH");
  if (path && path[0] != '\0')
    default_path = path;
  gfs_refine_terrain_class ();
  gfs_terrain_class ();
  gfs_variable_terrain_class ();
  return NULL;
}
