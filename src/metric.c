/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009-2012 National Institute of Water and Atmospheric Research
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
 * \brief Orthogonal curvilinear metric.
 */

#include <stdlib.h>
#include "metric.h"
#include <complex.h>
#include "map.h"
#include "solid.h"
#include "source.h"

#if USE_GSL
# include <gsl/gsl_integration.h>
#endif

/** \beginobject{GfsGenericMetric} */

static gdouble face_metric_direction (const GfsDomain * domain, FttCell * cell, FttDirection d)
{
  FttCellFace f;
  f.cell = cell;
  f.d = d;
  return (* domain->face_metric) (domain, &f);
}

/* see: doc/figures/viscous-metric.tm equation (4) */
static gdouble viscous_metric_implicit (const GfsDomain * domain,
					FttCell * cell,
					FttComponent component)
{
  g_assert (component < FTT_DIMENSION);
  /* fixme: 2D only */
  if (component > 1)
    return 0.;
  /* fixme: this does not include the "curvature" of the metric yet */
  FttComponent c1 = component;
  FttComponent c2 = (c1 + 1) % 2;
  double h1h2 = (* domain->cell_metric) (domain, cell);
  double size = ftt_cell_size (cell);
  double h1_2 = (face_metric_direction (domain, cell, 2*c2) - 
		 face_metric_direction (domain, cell, 2*c2 + 1));
  double h2_1 = (face_metric_direction (domain, cell, 2*c1) - 
		 face_metric_direction (domain, cell, 2*c1 + 1));
  return (h1_2*h1_2 + h2_1*h2_1)/(size*size*h1h2*h1h2);
}

/* see: doc/figures/viscous-metric.tm equation (4) */
static gdouble viscous_metric_explicit (const GfsDomain * domain, 
					FttCell * cell,
					GfsVariable * v,
					GfsDiffusion * d)
{
  g_assert (v->component < FTT_DIMENSION);
  /* fixme: 2D only */
  if (v->component > 1)
    return 0.;
  FttComponent c1 = v->component;
  FttComponent c2 = (c1 + 1) % 2;
  double h1h2 = (* domain->cell_metric) (domain, cell);
  double h1 = (* domain->scale_metric) (domain, cell, c1);
  double h2 = (* domain->scale_metric) (domain, cell, c2);
  double size = ftt_cell_size (cell);
  double h1_2 = (face_metric_direction (domain, cell, 2*c2) - 
		 face_metric_direction (domain, cell, 2*c2 + 1));
  double h2_1 = (face_metric_direction (domain, cell, 2*c1) - 
		 face_metric_direction (domain, cell, 2*c1 + 1));
  double u2_1 = gfs_center_gradient (cell, c1, v->vector[c2]->i);
  double u2_2 = gfs_center_gradient (cell, c2, v->vector[c2]->i);
  double eta = gfs_diffusion_cell (d, cell);
  /* fixme: this does not include the terms with derivatives of the viscosity yet */
  /* fixme: this does not include the "curvature" of the metric yet */
  /* fixme: this does not take into account density */
  return eta*(
	      + 2.*(u2_1*h1_2/h1 - u2_2*h2_1/h2)/(size*size)
	      )/h1h2;
}

static void advection_metric (const GfsDomain * domain, 
			      FttCell * cell,
			      FttComponent c1,
			      gdouble m[2])
{
  g_assert (c1 < FTT_DIMENSION);
  /* fixme: 2D only */
  g_assert (c1 <= 1);
  FttComponent c2 = (c1 + 1) % 2;
  double h1h2 = (* domain->cell_metric) (domain, cell);
  double size = ftt_cell_size (cell);
  double h1_2 = (face_metric_direction (domain, cell, 2*c2) - 
		 face_metric_direction (domain, cell, 2*c2 + 1));
  double h2_1 = (face_metric_direction (domain, cell, 2*c1) - 
		 face_metric_direction (domain, cell, 2*c1 + 1));

  m[0] = h1_2/(size*h1h2);
  m[1] = h2_1/(size*h1h2);
}

static void set_default_metric (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  if (domain->metric_data || domain->face_metric || domain->cell_metric) {
    gts_file_error (fp, "cannot use multiple metrics (yet)");
    return;
  }

  domain->viscous_metric_implicit = viscous_metric_implicit;
  domain->viscous_metric_explicit = viscous_metric_explicit;
  domain->advection_metric = advection_metric;
}

static void generic_metric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_generic_metric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  set_default_metric (o, fp);
}

static void generic_metric_class_init (GtsObjectClass * klass)
{
  klass->read = generic_metric_read;
}

GfsEventClass * gfs_generic_metric_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsGenericMetric",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) generic_metric_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsGenericMetric} */

/**
 * A generic class for metrics which require storage.
 * \beginobject{GfsVariableMetric}
 */

static void variable_metric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_metric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  set_default_metric (o, fp);
}

static void variable_metric_class_init (GtsObjectClass * klass)
{
  klass->read = variable_metric_read;
}

GfsVariableClass * gfs_variable_metric_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableMetric",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_metric_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableMetric} */

/** \beginobject{GfsStoredMetric} */

#define USE_GSL 0

/* Coefficients from http://sparse-grids.de/ */
/* KPU: Univariate nested quadrature rules as basis - delayed
   Kronrod-Patterson rules, see Knut Petras (2003): "Smolyak cubature
   of given polynomial degree with few nodes for increasing
   dimension." Numerische Mathematik 93, 729-753. */

static double kpu_1d_l3[3][2] = {
  {.1127017, .2777778},
  {.5,       .4444444},
  {.8872983, .2777778}
};

static double kpu_2d_l2[5][3] = {
  {.1127017, .5,       .2777778},
  {.5,       .1127017, .2777778},
  {.5,       .5,      -.1111112},
  {.5,       .8872983, .2777778},
  {.8872983, .5,       .2777778}
};

static double kpu_2d_l3[9][3] = {
  {.1127017, .1127017, .07716050617284001},
  {.1127017, .5,       .12345678765432},
  {.1127017, .8872983, .07716050617284001},
  {.5,       .1127017, .12345678765432},
  {.5,       .5,       .19753082469136002},
  {.5,       .8872983, .12345678765432},
  {.8872983, .1127017, .07716050617284001},
  {.8872983, .5,       .12345678765432},
  {.8872983, .8872983, .07716050617284001}
};

static double kpu_2d_l4[17][3] = {
  { .0197544, .5, .052328105232810528 },
  { .1127017, .1127017, .077160506172840024 },
  { .1127017, .5, -.020076998921278625 },
  { .1127017, .8872983, .077160506172840024 },
  { .2828781, .5, .20069872006987202 },
  { .5, .0197544, .052328105232810528 },
  { .5, .1127017, -.020076998921278667 },
  { .5, .2828781, .20069872006987202 },
  { .5, .5, -.24044133021697545 },
  { .5, .7171219, .20069872006987202 },
  { .5, .8872983, -.020076998921278667 },
  { .5, .9802456, .052328105232810528 },
  { .7171219, .5, .20069872006987202 },
  { .8872983, .1127017, .077160506172840024 },
  { .8872983, .5, -.020076998921278625 },
  { .8872983, .8872983, .077160506172840024 },
  { .9802456, .5, .052328105232810528 },
};

#define EPS 1e-6

static double ru_rv (FttVector r, GfsMap * map)
{
  FttVector ru = { r.x + EPS, r.y, 0. };
  FttVector rv = { r.x, r.y + EPS, 0. };
  (* map->inverse) (map, &r, &r);
  (* map->inverse) (map, &ru, &ru);
  (* map->inverse) (map, &rv, &rv);
  ru.x -= r.x; ru.y -= r.y; ru.z -= r.z;
  rv.x -= r.x; rv.y -= r.y; rv.z -= r.z;
  return sqrt ((ru.x*ru.x + ru.y*ru.y + ru.z*ru.z)*(rv.x*rv.x + rv.y*rv.y + rv.z*rv.z)
	       /* the cross-term should be zero for an orthogonal
		  metric but we keep them for clarity*/
	       - (ru.x*rv.x + ru.y*rv.y + ru.z*rv.z)*(ru.x*rv.x + ru.y*rv.y + ru.z*rv.z)
	       );
}

static double integration2d (GfsMap * map,
			     double u1, double v1,
			     double u2, double v2)
{
  int i;
  FttVector r;
  double du = u2 - u1;
  double dv = v2 - v1;
  double a = 0.;
  r.z = 0.;
  for (i = 0; i < 9; i++) {
    r.x = u1 + du*kpu_2d_l3[i][0];
    r.y = v1 + dv*kpu_2d_l3[i][1];
    a += kpu_2d_l3[i][2]*ru_rv (r, map);
  }
  return a*du*dv;
}

typedef struct {
  FttVector * r;
  GfsMap * map;
  gdouble v1, v2;
} RuRvData;

static double ru (FttVector r, GfsMap * map)
{
  FttVector dr = { r.x + EPS, r.y, 0. };
  (* map->inverse) (map, &r, &r);
  (* map->inverse) (map, &dr, &dr);
  dr.x -= r.x; dr.y -= r.y; dr.z -= r.z;
  return sqrt ((dr.x*dr.x + dr.y*dr.y + dr.z*dr.z));
}

static double rv (FttVector r, GfsMap * map)
{
  FttVector dr = { r.x, r.y + EPS, 0. };
  (* map->inverse) (map, &r, &r);
  (* map->inverse) (map, &dr, &dr);
  dr.x -= r.x; dr.y -= r.y; dr.z -= r.z;
  return sqrt ((dr.x*dr.x + dr.y*dr.y + dr.z*dr.z));
}

/* Returns: \sqrt{(r_u.r_u)} */

#if USE_GSL
static double ru_gsl (double u, void * data)
{
  RuRvData * p = data;
  p->r->x = u;
  return ru (*(p->r), p->map)/EPS;
}

/* Returns: \sqrt{(r_v.r_v)} */
static double rv_gsl (double v, void * data)
{
  RuRvData * p = data;
  p->r->y = v;
  return rv (*(p->r), p->map)/EPS;
}

static gdouble integration (const gsl_function * f, double a, double b)
{
  double result, abserr;
  size_t neval;
  /* we set the error to one but the QNG code will use at least
     21 points which is enough */
  gsl_integration_qng (f, a, b, 1., 0., &result, &abserr, &neval);
  //  fprintf (stderr, "neval: %d abserr: %g result: %g\n", neval, abserr, result);
  return result;
}
#endif /* USE_GSL */

/* Returns: \int \sqrt{(r_u.r_u)} du */
static double length_u (GfsMap * map, double u1, double u2, double v)
{
#if USE_GSL
  FttVector r;
  RuRvData p = { &r, map };
  r.y = v;
  r.z = 0.;
  gsl_function f;
  f.function = ru_gsl;
  f.params = &p;
  return integration (&f, u1, u2);
#else
  int i;
  FttVector r;
  double du = u2 - u1;
  double a = 0.;
  r.y = v;
  r.z = 0.;
  for (i = 0; i < 3; i++) {
    r.x = u1 + du*kpu_1d_l3[i][0];
    a += kpu_1d_l3[i][1]*ru (r, map);
  }
  return a*du/EPS;
#endif
}

/* Returns: \int \sqrt{(r_v.r_v)} dv */
static double length_v (GfsMap * map, double v1, double v2, double u)
{
#if USE_GSL
  FttVector r;
  gsl_function f;
  RuRvData p = { &r, map };
  f.function = rv_gsl;
  f.params = &p;
  r.x = u;
  r.z = 0.;
  return integration (&f, v1, v2);
#else
  int i;
  FttVector r;
  double dv = v2 - v1;
  double a = 0.;
  r.x = u;
  r.z = 0.;
  for (i = 0; i < 3; i++) {
    r.y = v1 + dv*kpu_1d_l3[i][0];
    a += kpu_1d_l3[i][1]*rv (r, map);
  }
  return a*dv/EPS;
#endif
}

/* Returns: \sqrt{(r_u.r_u)(r_v.r_v) - (r_u.r_v)^2} */
#if USE_GSL
static double ru_rv_gsl (double v, void * data)
{
  RuRvData * p = data;
  p->r->y = v;
  return ru_rv (*(p->r), p->map)/(EPS*EPS);
}

/* Returns: \int \sqrt{(r_u.r_u)(r_v.r_v) - (r_u.r_v)^2} dv */
static double ru_rv_dv (double u, void * data)
{
  RuRvData * p = data;
  gsl_function f;
  f.function = ru_rv_gsl;
  f.params = p;
  p->r->x = u;
  return integration (&f, p->v1, p->v2);
}
#endif /* USE_GSL */

/* Returns: \int\int \sqrt{(r_u.r_u)(r_v.r_v) - (r_u.r_v)^2} du dv */
static double area (GfsMap * map,
		    double u1, double v1,
		    double u2, double v2)
{
#if USE_GSL
  FttVector r;
  gsl_function f;
  RuRvData p = { &r, map, v1, v2 };
  f.function = ru_rv_dv;
  f.params = &p;
  r.z = 0.;
  return integration (&f, u1, u2);
#else
  return integration2d (map, u1, v1, u2, v2)/(EPS*EPS);
#endif
}

static void metric_coarse_fine (FttCell * parent, GfsVariable * a)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  FttVector p;
  FttCellChildren child;
  gdouble h = ftt_cell_size (parent)/2., sa = 0.;
  ftt_cell_children (parent, &child);
  int i;

  GfsStoredMetric * m = GFS_STORED_METRIC (a);
  GfsMap * map = m->map;
  for (i = 0; i < FTT_CELLS; i++) {
    ftt_cell_pos (child.c[i], &p);
    GFS_VALUE (child.c[i], a) = area (map, 
				      p.x - h/2., p.y - h/2., 
				      p.x + h/2., p.y + h/2.)/(h*h);
    sa += GFS_VALUE (child.c[i], a);
  }

  if (m->e) {
    double e = GFS_VALUE (parent, a) - sa/4.;
    for (i = 0; i < FTT_CELLS; i++)
      GFS_VALUE (child.c[i], m->e) = e;
  }

  ftt_cell_pos (parent, &p);

  GFS_VALUE (child.c[0], m->h[0]) = GFS_VALUE (child.c[1], m->h[1]) = 
    length_v (map, p.y, p.y + h, p.x)/h;
  GFS_VALUE (child.c[0], m->h[3]) = GFS_VALUE (child.c[2], m->h[2]) = 
    length_u (map, p.x - h, p.x, p.y)/h;
  GFS_VALUE (child.c[2], m->h[0]) = GFS_VALUE (child.c[3], m->h[1]) = 
    length_v (map, p.y - h, p.y, p.x)/h;
  GFS_VALUE (child.c[1], m->h[3]) = GFS_VALUE (child.c[3], m->h[2]) = 
    length_u (map, p.x, p.x + h, p.y)/h;

  GFS_VALUE (child.c[0], m->h[2]) = length_u (map, p.x - h, p.x, p.y + h)/h;
  GFS_VALUE (child.c[0], m->h[1]) = length_v (map, p.y, p.y + h, p.x - h)/h;
  GFS_VALUE (child.c[1], m->h[2]) = length_u (map, p.x, p.x + h, p.y + h)/h;
  GFS_VALUE (child.c[1], m->h[0]) = length_v (map, p.y, p.y + h, p.x + h)/h;
  GFS_VALUE (child.c[2], m->h[3]) = length_u (map, p.x - h, p.x, p.y - h)/h;
  GFS_VALUE (child.c[2], m->h[1]) = length_v (map, p.y - h, p.y, p.x - h)/h;
  GFS_VALUE (child.c[3], m->h[3]) = length_u (map, p.x, p.x + h, p.y - h)/h;
  GFS_VALUE (child.c[3], m->h[0]) = length_v (map, p.y - h, p.y, p.x + h)/h;
}

static void metric_fine_coarse (FttCell * parent, GfsVariable * a)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  gdouble va = 0.;
  for (n = 0; n < 4; n++)
    va += GFS_VALUE (child.c[n], a);
  GFS_VALUE (parent, a) = va/4.;

  GfsStoredMetric * m = GFS_STORED_METRIC (a);
  GFS_VALUE (parent, m->h[0]) = (GFS_VALUE (child.c[1], m->h[0]) +
				 GFS_VALUE (child.c[3], m->h[0]))/2.;
  GFS_VALUE (parent, m->h[1]) = (GFS_VALUE (child.c[0], m->h[1]) +
				 GFS_VALUE (child.c[2], m->h[1]))/2.;
  GFS_VALUE (parent, m->h[2]) = (GFS_VALUE (child.c[0], m->h[2]) +
				 GFS_VALUE (child.c[1], m->h[2]))/2.;
  GFS_VALUE (parent, m->h[3]) = (GFS_VALUE (child.c[2], m->h[3]) +
				 GFS_VALUE (child.c[3], m->h[3]))/2.;
}

static gdouble face_metric (const GfsDomain * domain, const FttCellFace * face)
{ 
  if (face->d/2 > FTT_Y)
    return 1.;
  return GFS_VALUE (face->cell, GFS_STORED_METRIC (domain->metric_data)->h[face->d]);
}

static gdouble cell_metric (const GfsDomain * domain, const FttCell * cell)
{
  return GFS_VALUE (cell, GFS_VARIABLE (domain->metric_data));
}

static void solid_metric (const GfsDomain * domain, const FttCell * cell, FttVector * m)
{
  g_assert (GFS_IS_MIXED (cell));
  g_assert_not_implemented ();
}

static gdouble scale_metric (const GfsDomain * domain, const FttCell * cell, FttComponent c)
{
  /* fixme: this does not allow for Z-metric */
  if (c > FTT_Y)
    return 1.;
  FttComponent d = FTT_ORTHOGONAL_COMPONENT (c);
  return (GFS_VALUE (cell, GFS_STORED_METRIC (domain->metric_data)->h[2*d]) +
	  GFS_VALUE (cell, GFS_STORED_METRIC (domain->metric_data)->h[2*d + 1]))/2.;
}

static gdouble face_scale_metric (const GfsDomain * domain, const FttCellFace * face,
				  FttComponent c)
{
  /* fixme: this does not allow for Z-metric */
  if (c > FTT_Y)
    return 1.;
  /* fixme: this is not second-order for fine/coarse faces */
  return (scale_metric (domain, face->cell, c) + scale_metric (domain, face->neighbor, c))/2.;
}

static void none (FttCell * parent, GfsVariable * v)
{
}

static void stored_metric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_stored_metric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsStoredMetric * m = GFS_STORED_METRIC (*o);
  if (fp->type == GTS_STRING) {
    if (!(m->e = gfs_domain_get_or_add_variable (domain, fp->token->str, "Metric error"))) {
      gts_file_error (fp, "`%s' is a reserved variable name", fp->token->str);
      return;
    }
    m->e->fine_coarse = m->e->coarse_fine = none;
    gts_file_next_token (fp);
  }

  GfsVariable * a = GFS_VARIABLE (*o);
  FttDirection d;
  for (d = 0; d < 4; d++) {
    gchar * name = g_strdup_printf ("%sh%d", a->name, d);
    m->h[d] = gfs_domain_get_or_add_variable (domain, name, "Face metric");
    m->h[d]->fine_coarse = m->h[d]->coarse_fine = none;
    g_free (name);
  }

  g_free (a->description);
  a->description = g_strdup ("Cell metric");
  a->coarse_fine = metric_coarse_fine;
  a->fine_coarse = metric_fine_coarse;

  m->map = GFS_MAP (gts_object_new (GTS_OBJECT_CLASS (m->map_class)));
  gfs_object_simulation_set (m->map, domain);
  gts_container_add (GTS_CONTAINER (GFS_SIMULATION (domain)->maps), GTS_CONTAINEE (m->map));

  domain->metric_data = *o;
  domain->face_metric  = face_metric;
  domain->cell_metric  = cell_metric;
  domain->solid_metric = solid_metric;
  domain->scale_metric = scale_metric;
  domain->face_scale_metric = face_scale_metric;
}

static void stored_metric_class_init (GtsObjectClass * klass)
{
  klass->read = stored_metric_read;
}

static void stored_metric_init (GfsStoredMetric * m)
{
  m->map_class = gfs_map_class ();
}

GfsVariableClass * gfs_stored_metric_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsStoredMetric",
      sizeof (GfsStoredMetric),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) stored_metric_class_init,
      (GtsObjectInitFunc) stored_metric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_metric_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsStoredMetric} */

/* GfsMapMetric: Header */

#define GFS_IS_MAP_METRIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_metric_class ()))

static GfsMapClass * gfs_map_metric_class      (void);

/* GfsMapMetric: Object */

static void gfs_map_metric_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetric */
}

static void gfs_map_metric_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetric */
}

static void gfs_map_metric_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_metric_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_metric_write;
}

static void map_metric_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMetric * m = GFS_DOMAIN (gfs_object_simulation (map))->metric_data;
  FttVector src1 = *src; /* just in case src == dest */
  FttComponent c;
  for (c = 0; c < 3; c++)
    if ((&m->x)[c])
      (&dest->x)[c] = gfs_function_spatial_value ((&m->x)[c], &src1);
    else
      (&dest->x)[c] = (&src1.x)[c];
}

static void gfs_map_metric_init (GfsMap * map)
{
  map->inverse =   map_metric_inverse;
}

static GfsMapClass * gfs_map_metric_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_metric_info = {
      "GfsMapMetric",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_metric_class_init,
      (GtsObjectInitFunc) gfs_map_metric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_metric_info);
  }

  return klass;
}

/** \beginobject{GfsMetric} */

static void metric_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_class ())->parent_class->write) (o, fp);
  
  GfsMetric * m = GFS_METRIC (o);
  fputs (" {", fp);
  FttComponent c;
  static gchar name[3][2] = {"x", "y", "z"};
  for (c = 0; c < 3; c++)
    if ((&m->x)[c]) {
      fprintf (fp, "\n    %s = ", name[c]);
      gfs_function_write ((&m->x)[c], fp);
    }
  fputs ("\n  }", fp);
}

static void metric_destroy (GtsObject * o)
{
  GfsMetric * m = GFS_METRIC (o);
  FttComponent c;
  for (c = 0; c < 3; c++)
    if ((&m->x)[c])
      gts_object_destroy (GTS_OBJECT ((&m->x)[c]));

  (* GTS_OBJECT_CLASS (gfs_metric_class ())->parent_class->destroy) (o);
}

static void metric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting a parameter block");
    return;
  }

  GfsMetric * m = GFS_METRIC (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (m));
  GtsFileVariable var[] = {
    {GTS_OBJ, "x", TRUE, &m->x},
    {GTS_OBJ, "y", TRUE, &m->y},
    {GTS_OBJ, "z", TRUE, &m->z},
    {GTS_NONE}
  };
  FttComponent c;
  for (c = 0; c < 3; c++)
    gfs_object_simulation_set ((&m->x)[c], domain);
  
  gts_file_assign_variables (fp, var);
  
  for (c = 0; c < 3; c++)
    if (!var[c].set) {
      gts_object_destroy (GTS_OBJECT ((&m->x)[c]));
      (&m->x)[c] = NULL;
    }
  
  if (fp->type == GTS_ERROR)
    return;
}

static void metric_class_init (GtsObjectClass * klass)
{
  klass->destroy = metric_destroy;
  klass->read = metric_read;
  klass->write = metric_write;
}

static void metric_init (GfsMetric * m)
{
  GFS_STORED_METRIC (m)->map_class = gfs_map_metric_class ();
  m->x = gfs_function_new (gfs_function_map_class (), 1.);
  m->y = gfs_function_new (gfs_function_map_class (), 1.);
  m->z = gfs_function_new (gfs_function_map_class (), 1.);
}

GfsVariableClass * gfs_metric_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_metric_info = {
      "GfsMetric",
      sizeof (GfsMetric),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) metric_class_init,
      (GtsObjectInitFunc) metric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_stored_metric_class ()),
				  &gfs_metric_info);
  }

  return klass;
}

/** \endobject{GfsMetric} */

/* "Expanded spherical cube" metric */

#define N 30

#if 1
/* Conformal mapping Taylor coefficients: from Rancic et al, 1996, Table B.1 */

static double A[N] = {
   1.47713062600964, -0.38183510510174, -0.05573058001191, -0.00895883606818, -0.00791315785221,
  -0.00486625437708, -0.00329251751279, -0.00235481488325, -0.00175870527475, -0.00135681133278,
  -0.00107459847699, -0.00086944475948, -0.00071607115121, -0.00059867100093, -0.00050699063239,
  -0.00043415191279, -0.00037541003286, -0.00032741060100, -0.00028773091482, -0.00025458777519,
  -0.00022664642371, -0.00020289261022, -0.00018254510830, -0.00016499474461, -0.00014976117168,
  -0.00013646173946, -0.00012478875823, -0.00011449267279, -0.00010536946150, -0.00009725109376
};

static double B[N] = {
  0.67698819751739, 0.11847293456554, 0.05317178134668, 0.02965810434052, 0.01912447304028,
  0.01342565621117, 0.00998873323180, 0.00774868996406, 0.00620346979888, 0.00509010874883,
  0.00425981184328, 0.00362308956077, 0.00312341468940, 0.00272360948942, 0.00239838086555,
  0.00213001905118, 0.00190581316131, 0.00171644156404, 0.00155493768255, 0.00141600715207,
  0.00129556597754, 0.00119042140226, 0.00109804711790, 0.00101642216628, 0.00094391366522,
  0.00087919021224, 0.00082115710311, 0.00076890728775, 0.00072168382969, 0.00067885087750
};

#else
/* Conformal mapping Taylor coefficients: from map_xy2xyz.m from mitgcm */

static double A[N] = {
  1.47713057321600, -0.38183513110512, -0.05573055466344, -0.00895884801823, -0.00791314396396,
  -0.00486626515498, -0.00329250387158, -0.00235482619663, -0.00175869000970, -0.00135682443774,
  -0.00107458043205, -0.00086946107050, -0.00071604933286, -0.00059869243613, -0.00050696402446,
  -0.00043418115349, -0.00037537743098, -0.00032745130951, -0.00028769063795, -0.00025464473946,
  -0.00022659577923, -0.00020297175587, -0.00018247947703, -0.00016510295548, -0.00014967258633,
  -0.00013660647356, -0.00012466390509, -0.00011468147908, -0.00010518717478, -0.00009749136078
};

static double B[N] = {
  0.67698822171341, 0.11847295533659, 0.05317179075349, 0.02965811274764, 0.01912447871071,
  0.01342566129383, 0.00998873721022, 0.00774869352561, 0.00620347278164, 0.00509011141874,
  0.00425981415542, 0.00362309163280, 0.00312341651697, 0.00272361113245, 0.00239838233411,
  0.00213002038153, 0.00190581436893, 0.00171644267546, 0.00155493871562, 0.00141600812949,
  0.00129556691848, 0.00119042232809, 0.00109804804853, 0.00101642312253, 0.00094391466713,
  0.00087919127990, 0.00082115825576, 0.00076890854394, 0.00072168520663, 0.00067885239089
};
#endif

static complex double WofZ (complex double Z)
{
  complex double W = 0.;
  int n = N;
  while (n-- > 0)
    W = (W + A[n])*Z;
  return W;
}

static complex double ZofW (complex double W)
{
  complex double Z = 0.;
  int n = N;
  while (n-- > 0)
    Z = (Z + B[n])*W;
  return Z;
}

/* I^(1/3) */
#define I3 (0.86602540378444 + I/2.)
/* sqrt (3.) - 1. */
#define RA 0.73205080756888

/* Conformal mapping of a cube face onto a sphere. Maps (x,y) on the
 * north-pole face of a cube to (X,Y,Z) coordinates in physical space.
 *
 * Based on f77 code from Jim Purser & Misha Rancic.
 *
 * Face is oriented normal to Z-axis with X and Y increasing with x
 * and y.
 */
static void fmap_xy2XYZ (double x, double y, double * X, double * Y, double * Z)
{
  int kx = x < 0., ky = y < 0.;
  x = fabs (x); y = fabs (y);
  int kxy = y > x;

  if (kxy) {
    double tmp = x;
    x = 1. - y;
    y = 1. - tmp;
  }
  else {
    x = 1. - x;
    y = 1. - y;
  }

  complex double z = (x + I*y)/2.;
  complex double W;
  if (cabs (z) > 0.) {
    W = WofZ (z*z*z*z);
    W = I3*cpow (W*I, 1./3.);
  }
  else
    W = 0.;
  complex double cb = I - 1.;
  complex double cc = RA*cb/2.;
  W = (W - RA)/(cb + cc*W);
  *X = creal (W);
  *Y = cimag (W);
  double H = 2./(1. + (*X)*(*X) + (*Y)*(*Y));
  *X *= H;
  *Y *= H;
  *Z = H - 1.;
  
  if (kxy) {
    double tmp = *X;
    *X = *Y;
    *Y = tmp;
  }
  if (kx)
    *X = - *X;
  if (ky)
    *Y = - *Y;
}

/* Conformal mapping of a cube onto a sphere. Maps (x,y) on the
 * 6 faces of the cube to (X,Y,Z) coordinates in physical space.
 *
 * Based on f77 code from Jim Purser & Misha Rancic.
 *
 * Face 1 is oriented normal to Z-axis with X and Y increasing with x
 * and y (see doc/figures/cubed.fig).
 *
 * returns: FALSE if the input coordinates are invalid, TRUE otherwise.
 */
static void cmap_xy2XYZ (double x, double y, double * X, double * Y, double * Z)
{
  x *= 2.; y *= 2.;

  /* fixme: causes crash in gfsview when saving in gnuplot format */
  //  g_assert (x >= -1. && x <= 7. && y >= -1. && y <= 5.);

  /* symmetries: see doc/figures/cubed.fig */
  double tmp;
  if (y <= 1. && x <= 3.) {
    if (x <= 1.) /* face 1 */
      fmap_xy2XYZ (x, y, X, Y, Z);
    else { /* face 2 */
      fmap_xy2XYZ (x - 2., y, X, Y, Z);
      tmp = *X;
      *X = *Z; *Z = - tmp;
    }
  }
  else if (y <= 3. && x <= 5.) {
    if (x <= 3.) { /* face 3 */
      fmap_xy2XYZ (x - 2., y - 2., X, Y, Z);
      tmp = *X;
      *X = -*Y; *Y = *Z; *Z = - tmp;
    }
    else { /* face 4 */
      fmap_xy2XYZ (x - 4., y - 2., X, Y, Z);
      tmp = *Y;
      *Z = - *Z; *Y = - *X; *X = - tmp;
    }
  }
  else {
    if (x <= 5.) { /* face 5 */
      fmap_xy2XYZ (x - 4., y - 4., X, Y, Z);
      tmp = *Z;
      *Z = *Y; *Y = - *X; *X = - tmp;
    }
    else { /* face 6 */
      fmap_xy2XYZ (x - 6., y - 4., X, Y, Z);
      tmp = *Y;
      *Y = - *Z; *Z = tmp;
    }
  }
}

/* Conformal mapping of a sphere onto a cube face. Maps (X,Y,Z) coordinates
 * in physical space to (x,y) on the north-pole face of a cube.
 *
 * This is the inverse transform of fmap_xy2XYZ().
 */
static void fmap_XYZ2xy (double X, double Y, double Z, double * x, double * y)
{
  int kx = X < 0., ky = Y < 0.;
  X = fabs (X); Y = fabs (Y);
  int kxy = Y > X;

  if (kxy) {
    double tmp = X;
    X = Y;
    Y = tmp;
  }

  double H = Z + 1.;
  X /= H; Y /= H;
  complex double W = X + Y*I;
  complex double cb = I - 1.;
  complex double cc = RA*cb/2.;
  W = (W*cb + RA)/(1. - W*cc);
  W = W/I3;
  W = W*W*W;
  W /= I;
  complex double z = ZofW (W);
  z = cpow (z, 1./4.)*2.;
  *x = fabs (creal (z));
  *y = fabs (cimag (z));

  if (kxy) {
    *x = 1. - *x;
    *y = 1. - *y;
  }
  else {
    double tmp = *x;
    *x = 1. - *y;
    *y = 1. - tmp;
  }
  if (kx)
    *x = - *x;
  if (ky)
    *y = - *y;
}

static double angle_between_vectors (const double v1[3], const double v2[3])
{
  return acos (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

static void plane_normal (const double p1[3], const double p2[3], double plane[3])
{
  plane[0] = p1[1]*p2[2] - p1[2]*p2[1];
  plane[1] = p1[2]*p2[0] - p1[0]*p2[2];
  plane[2] = p1[0]*p2[1] - p1[1]*p2[0];
  double mag = sqrt (plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
  plane[0] /= mag;
  plane[1] /= mag;
  plane[2] /= mag;
}

static double excess_of_quad (const double v1[3], const double v2[3],
			      const double v3[3], const double v4[3])
{
  double plane1[3], plane2[3], plane3[3], plane4[3];

  plane_normal (v1, v2, plane1);
  plane_normal (v2, v3, plane2);
  plane_normal (v3, v4, plane3);
  plane_normal (v4, v1, plane4);
  
  return 2.*M_PI -
    angle_between_vectors (plane2, plane1) -
    angle_between_vectors (plane3, plane2) -
    angle_between_vectors (plane4, plane3) -
    angle_between_vectors (plane1, plane4);
}

/* GfsMapCubed: Object */

static void gfs_map_cubed_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed */
}

static void gfs_map_cubed_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed */
}

static void gfs_map_cubed_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_cubed_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_cubed_write;
}

/* Returns the index of the face of the cube containing */
/* the point of coordinates X,Y,Z */
static int face_num (gdouble X, gdouble Y, gdouble Z)
{
  if (fabs(X) < Z && fabs(Y) < Z)
    return 1;
  else if (fabs(X) > fabs(Y) && X > fabs(Z))
    return 2;
  else if (fabs(X) < fabs(Y) && Y > fabs(Z))
    return 3;
  else if (-fabs(X) > Z && -fabs(Y) > Z)
    return 4;
  else if (fabs(X) > fabs(Y) && -X > fabs(Z))
    return 5;
  else
    return 6;
}

static void map_cubed_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsSimulation * sim = gfs_object_simulation (map);
  double lon = src->x*sim->physical_params.L*M_PI/180.;
  double lat = src->y*sim->physical_params.L*M_PI/180.;
  double coslat = cos (lat), sinlat = sin (lat), coslon = cos (lon);
  double X = coslat*sin (lon), Y = sinlat, Z = coslat*coslon;
  double x, y;

  /* Maybe not the most elegant but works */
  switch (face_num(X,Y,Z)) {
  case 1:
    fmap_XYZ2xy (X, Y, fabs(Z), &x, &y);
    dest->x = x/2.;
    dest->y = y/2.;
    dest->z = src->z;
    break;
  case 2:
    X = - coslat*coslon;
    Z = sqrt (1. - X*X - Y*Y);
    fmap_XYZ2xy (X, Y, Z, &x, &y);
    dest->x = (1. + x/2.);
    dest->y = y/2.;
    dest->z = src->z;
    break;
  case 3:
    if (M_PI/4. < fabs(lon) && fabs(lon) < 3.*M_PI/4.) {
      X = - coslat*coslon;
      Z = sqrt (1. - X*X - Y*Y);
      fmap_XYZ2xy (X, Y, Z, &x, &y);
      dest->x = (1. + x/2.);
      if (lon < 0.)
	dest->y = (1. + y/2.);
      else
	dest->y = (1. - y/2.);
      dest->z = src->z;
    }
    else {
      fmap_XYZ2xy (X, Y, fabs(Z), &x, &y);
      if (lon > -3.*M_PI/4. && lon < 3.*M_PI/4.)
	dest->x = (1. - y/2.);
      else
	dest->x = (1. + y/2.);
      dest->y = (1. - x/2.);
      dest->z = src->z;
    }
    break;
  case 4:
    fmap_XYZ2xy (X, Y, fabs(Z), &x, &y);
    dest->x = (2. - y/2.);
    dest->y = (1. - x/2.);
    dest->z = src->z;
    break;
  case 5:
    X = coslat*coslon;
    Z = sqrt (1. - X*X - Y*Y);
    fmap_XYZ2xy (X, Y, Z, &x, &y);
    dest->x = (2. - y/2.);
    dest->y = (2. + x/2.);
    dest->z = src->z;
    break;
  case 6:
    if (M_PI/4. < fabs(lon) && fabs(lon) < 3.*M_PI/4.) {
      X = - coslat*coslon;
      Z = sqrt (1. - X*X - Y*Y);
      fmap_XYZ2xy (X, Y, Z, &x, &y);
      dest->y = (2. - x/2.);
      if (lon < 0.)
	dest->x = (3. + y/2.);
      else
	dest->x = (3. - y/2.);
      dest->z = src->z;
    }
    else {
      fmap_XYZ2xy (X, Y, fabs(Z), &x, &y);
      if (lon > -3.*M_PI/4. && lon < 3.*M_PI/4.)
	dest->y = (2. - y/2.);
      else
	dest->y = (2. + y/2.);
      dest->x = (3. + x/2.);
      dest->z = src->z;
    }
    break;
  default:
    g_assert_not_reached ();
  }
}

static void map_cubed_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsSimulation * sim = gfs_object_simulation (map);
  double X, Y, Z;
  cmap_xy2XYZ (src->x, src->y, &X, &Y, &Z);
  dest->x = atan2 (X, Z)*180./M_PI/sim->physical_params.L;
  dest->y = asin (Y)*180./M_PI/sim->physical_params.L;
  dest->z = src->z;
}

static void map_cubed_inverse_cell (GfsMap * map, const FttVector * src, FttVector * dest)
{
  gint i;
  FttVector o = { 0., 0., 0. };
  for (i = 0; i < 4; i++) {
    o.x += src[i].x;
    o.y += src[i].y;
    o.z += src[i].z;
    map_cubed_inverse (map, &(src[i]), &(dest[i]));
  }
  o.x /= 4.; o.y /= 4.; o.z /= 4.;
  map_cubed_inverse (map, &o, &o);
  /* make sure we do not cross periodic longitude boundary */
  gdouble L = gfs_object_simulation (map)->physical_params.L;
  for (i = 0; i < 4; i++)
    if (dest[i].x > o.x + 180./L)
      dest[i].x -= 360./L;
    else if (dest[i].x < o.x - 180./L)
      dest[i].x += 360./L;
}

static double evaluate (const FttVector * x, const FttVector * rhs, FttVector * f)
{
  gdouble delta = 0.;
  cmap_xy2XYZ (x->x, x->y, &f->x, &f->y, &f->z);
  int i;
  for (i = 0; i < 3; i++) {
    (&f->x)[i] -= (&rhs->x)[i];
    delta += (&f->x)[i]*(&f->x)[i];
  }
  return delta;
}

#define DELTA 1e-6

static void jacobian (const FttVector * x, const FttVector * rhs, FttVector * f,
		      GtsMatrix * J)
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      FttVector df, dx = *x;
      (&dx.x)[j] += DELTA;
      cmap_xy2XYZ (dx.x, dx.y, &df.x, &df.y, &df.z);
      df.x *= 1. + dx.z;
      df.y *= 1. + dx.z;
      df.z *= 1. + dx.z;
      J[i][j] = ((&df.x)[i] - (&rhs->x)[i] - (&f->x)[i])/DELTA;
    }
}

static void normalized_jacobian (const FttVector * p, GtsMatrix * J)
{
  FttVector f, rhs = {0., 0., 0.};
  g_assert (p->z == 0.);
  evaluate (p, &rhs, &f);
  jacobian (p, &rhs, &f, J);
  /* normalize */
  int i, j;
  for (i = 0; i < 3; i++) {
    gdouble h = 0.;
    for (j = 0; j < 3; j++)
      h += J[j][i]*J[j][i];
    h = sqrt (h);
    for (j = 0; j < 3; j++)
      J[j][i] /= h;
  }
}

static void map_cubed_inverse_vector (GfsMap * map, const FttVector * p,
				      const FttVector * src, FttVector * dest)
{
  GtsMatrix J[4];
  normalized_jacobian (p, J);
  FttVector src1;

  int i, j;
  for (i = 0; i < 3; i++) {
    (&src1.x)[i] = 0.;
    for (j = 0; j < 3; j++)
      (&src1.x)[i] += (&src->x)[j]*J[i][j];
  }
  
  FttVector p1;
  map_cubed_inverse (map, p, &p1);
  GfsSimulation * sim = gfs_object_simulation (map);
  gdouble lon = p1.x*sim->physical_params.L*M_PI/180.;
  gdouble lat = p1.y*sim->physical_params.L*M_PI/180.;
  dest->x = (src1.x + tan (lat)*sin (lon)*src1.y)/cos (lon);
  dest->y = src1.y/cos (lat);
  dest->z = src1.z;
}

static void map_cubed_transform_vector (GfsMap * map, const FttVector * p,
					const FttVector * src, FttVector * dest)
{
  GtsMatrix J[4];
  normalized_jacobian (p, J);
  GtsMatrix * iJ = gts_matrix3_inverse (J);
  if (!iJ) {
    gts_matrix_print (J, stderr);
    g_assert_not_reached ();
  }

  FttVector src1 = *src;
  FttVector p1;
  map_cubed_inverse (map, p, &p1);
  GfsSimulation * sim = gfs_object_simulation (map);
  gdouble lon = p1.x*sim->physical_params.L*M_PI/180.;
  gdouble lat = p1.y*sim->physical_params.L*M_PI/180.;
  gdouble coslon = cos (lon), sinlon = sin (lon);
  gdouble sinlat = sin (lat);
  src1.x = coslon*src->x - sinlat*sinlon*src->y;
  src1.y = cos(lat)*src->y;
  src1.z = - sinlon*src->x - sinlat*coslon*src->y;

  int i, j;
  for (i = 0; i < 3; i++) {
    (&dest->x)[i] = 0.;
    for (j = 0; j < 3; j++)
      (&dest->x)[i] += (&src1.x)[j]*iJ[i][j];
  }
  gts_matrix_destroy (iJ);
}

static void gfs_map_cubed_init (GfsMap * map)
{
  map->transform = map_cubed_transform;
  map->inverse =   map_cubed_inverse;
  map->inverse_cell = map_cubed_inverse_cell;
  map->inverse_vector = map_cubed_inverse_vector;
  map->transform_vector = map_cubed_transform_vector;
}

static GfsMapClass * gfs_map_cubed_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_cubed_info = {
      "GfsMapCubed",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_cubed_class_init,
      (GtsObjectInitFunc) gfs_map_cubed_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_cubed_info);
  }

  return klass;
}

/**
 * The 'cubed sphere' metric.
 * \beginobject{GfsMetricCubed}
 */

static gdouble cubed_face_scale_metric (const GfsDomain * domain, const FttCellFace * face,
					FttComponent c)
{
  if (c > FTT_Y)
    return 1.;
  /* The transformation is conformal so the metric is isotropic */
  return GFS_VALUE (face->cell, GFS_STORED_METRIC (domain->metric_data)->h[face->d]);
}

typedef struct {
  double x, y, z;
  double x1, y1, z1;
  double a;
} Point;

static void point_new (Point * p, double x, double y)
{
  p->x = x; p->y = y;
  cmap_xy2XYZ (x, y, &p->x1, &p->y1, &p->z1);
}

static Point ** matrix_refine (Point ** m, int n)
{
  int n1 = 2*n - 1, i, j;
  Point ** r = gfs_matrix_new (n1, n1, sizeof (Point));
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      r[2*i][2*j] = m[i][j];
  for (i = 0; i < n - 1; i++)
    for (j = 0; j < n - 1; j++) {
      point_new (&r[2*i+1][2*j], (m[i][j].x + m[i+1][j].x)/2., m[i][j].y);
      point_new (&r[2*i][2*j+1], m[i][j].x, (m[i][j].y + m[i][j+1].y)/2.);
      point_new (&r[2*i+1][2*j+1], (m[i][j].x + m[i+1][j].x)/2., (m[i][j].y + m[i][j+1].y)/2.);
    }
  i = n - 1;
  for (j = 0; j < n - 1; j++)
    point_new (&r[2*i][2*j+1], m[i][j].x, (m[i][j].y + m[i][j+1].y)/2.);
  j = n - 1;
  for (i = 0; i < n - 1; i++)
    point_new (&r[2*i+1][2*j], (m[i][j].x + m[i+1][j].x)/2., m[i][j].y);
  gfs_matrix_free (m);
  return r;
}

static Point ** matrix_from_cell (FttCell * cell)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  double h = ftt_cell_size (cell)/2.;
  Point ** r = gfs_matrix_new (2, 2, sizeof (Point));
  point_new (&r[0][0], p.x - h, p.y - h);
  point_new (&r[1][0], p.x + h, p.y - h);
  point_new (&r[1][1], p.x + h, p.y + h);
  point_new (&r[0][1], p.x - h, p.y + h);
  return r;
}

static double matrix_a (Point ** r, int m, int i0, int j0)
{
  int i, j;
  double a = 0.;
  double h = r[m][0].x - r[0][0].x;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      a += excess_of_quad (&r[i0+i][j0+j].x1, &r[i0+i+1][j0+j].x1,
			   &r[i0+i+1][j0+j+1].x1, &r[i0+i][j0+j+1].x1);
  return 4.*a/(M_PI*M_PI*h*h);
}

static double matrix_hx (Point ** r, int m, int i0, int j0)
{
  int i;
  double hx = 0.;
  double h = r[m][0].x - r[0][0].x;
  for (i = 0; i < m; i++)
    hx += angle_between_vectors (&r[i0+i][j0].x1, &r[i0+i+1][j0].x1);
  return 2.*hx/(M_PI*h);
}

static double matrix_hy (Point ** r, int m, int i0, int j0)
{
  int j;
  double hy = 0.;
  double h = r[m][0].x - r[0][0].x;
  for (j = 0; j < m; j++)
    hy += angle_between_vectors (&r[i0][j0+j].x1, &r[i0][j0+j+1].x1);
  return 2.*hy/(M_PI*h);
}

static void cubed_coarse_fine (FttCell * parent, GfsVariable * a)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  GfsMetricCubed * cubed = GFS_METRIC_CUBED (a);
  Point ** r = matrix_from_cell (parent);
  r = matrix_refine (r, 2);
  int n = 3, level = cubed->level - (ftt_cell_level (parent) + 1);
  while (level-- > 0) {
    r = matrix_refine (r, n);
    n = 2*n - 1;
  }

  FttCellChildren child;
  ftt_cell_children (parent, &child);
  int m = n/2;

  GFS_VALUE (child.c[0], a) = matrix_a (r, m, 0, m);
  GFS_VALUE (child.c[1], a) = matrix_a (r, m, m, m);
  GFS_VALUE (child.c[2], a) = matrix_a (r, m, 0, 0);
  GFS_VALUE (child.c[3], a) = matrix_a (r, m, m, 0);

  GfsStoredMetric * metric = GFS_STORED_METRIC (a);
  GFS_VALUE (child.c[0], metric->h[0]) = GFS_VALUE (child.c[1], metric->h[1]) = 
    matrix_hy (r, m, m, m);
  GFS_VALUE (child.c[0], metric->h[3]) = GFS_VALUE (child.c[2], metric->h[2]) = 
    matrix_hx (r, m, 0, m);
  GFS_VALUE (child.c[2], metric->h[0]) = GFS_VALUE (child.c[3], metric->h[1]) = 
    matrix_hy (r, m, m, 0);
  GFS_VALUE (child.c[1], metric->h[3]) = GFS_VALUE (child.c[3], metric->h[2]) = 
    matrix_hx (r, m, m, m);

  GFS_VALUE (child.c[0], metric->h[2]) = matrix_hx (r, m, 0, n - 1);
  GFS_VALUE (child.c[0], metric->h[1]) = matrix_hy (r, m, 0, m);
  GFS_VALUE (child.c[1], metric->h[2]) = matrix_hx (r, m, m, n - 1);
  GFS_VALUE (child.c[1], metric->h[0]) = matrix_hy (r, m, n - 1, m);
  GFS_VALUE (child.c[2], metric->h[3]) = matrix_hx (r, m, 0, 0);
  GFS_VALUE (child.c[2], metric->h[1]) = matrix_hy (r, m, 0, 0);
  GFS_VALUE (child.c[3], metric->h[0]) = matrix_hy (r, m, n - 1, 0);
  GFS_VALUE (child.c[3], metric->h[3]) = matrix_hx (r, m, m, 0);

  gfs_matrix_free (r);
}

static void metric_cubed_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_cubed_class ())->parent_class->write) (o, fp);
  if (GFS_METRIC_CUBED (o)->level != 0)
    fprintf (fp, " %d", GFS_METRIC_CUBED (o)->level);
}

static void metric_cubed_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_cubed_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsMetricCubed * cubed = GFS_METRIC_CUBED (*o);
  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (level)");
    return;
  }
  cubed->level = atoi (fp->token->str);
  gts_file_next_token (fp);

  GfsVariable * a = GFS_VARIABLE (*o);
  a->coarse_fine = cubed_coarse_fine;

  a->domain->face_scale_metric = cubed_face_scale_metric;
}

static void metric_cubed_class_init (GtsObjectClass * klass)
{
  klass->read = metric_cubed_read;
  klass->write = metric_cubed_write;
}

static void metric_cubed_init (GfsStoredMetric * m)
{
  m->map_class = gfs_map_cubed_class ();
}

GfsVariableClass * gfs_metric_cubed_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMetricCubed",
      sizeof (GfsMetricCubed),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) metric_cubed_class_init,
      (GtsObjectInitFunc) metric_cubed_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_stored_metric_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsMetricCubed} */

/* GfsMapLonLat: Header */

typedef struct _GfsMapLonLat         GfsMapLonLat;

struct _GfsMapLonLat {
  /*< private >*/
  GfsMap parent;

  /*< public >*/
  gdouble r;
};

#define GFS_MAP_LONLAT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsMapLonLat,\
					         gfs_map_lonlat_class ())
#define GFS_IS_MAP_LONLAT(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_lonlat_class ()))

static GfsMapClass * gfs_map_lonlat_class      (void);

static void gfs_map_lonlat_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricLonLat */
}

static void gfs_map_lonlat_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricLonLat */
}

static void gfs_map_lonlat_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_lonlat_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_lonlat_write;
}

static void map_lonlat_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  dest->x = src->x*M_PI/180.*GFS_MAP_LONLAT (map)->r;
  dest->y = src->y*M_PI/180.*GFS_MAP_LONLAT (map)->r;
  dest->z = src->z;
}

static void map_lonlat_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  double x = src->x*180./(M_PI*GFS_MAP_LONLAT (map)->r);

  dest->x = x < -180. ? x + 360. : x > 180. ? x - 360. : x;
  dest->y = src->y*180./(M_PI*GFS_MAP_LONLAT (map)->r);
  dest->z = src->z;
}

static void map_lonlat_inverse_cell (GfsMap * map, const FttVector * src, FttVector * dest)
{
  gint i;
  FttVector q;

  q.x = (src[0].x + src[1].x)/2.;
  q.y = (src[0].y + src[2].y)/2.;
  q.z = 0.;

  for (i = 0; i < 4; i++)
    (* map->inverse) (map, &(src[i]), &(dest[i]));
  
  /* Fix for cells that contain the -180/180 degrees longitude line */
  if (dest[0].x < dest[1].x)  {
    (* map->inverse) (map, &q, &q);
    if (q.x > dest[0].x) {
      dest[0].x = 2.*q.x - dest[1].x;
      dest[3].x = dest[0].x;
    }
    else {
      dest[1].x = 2.*q.x - dest[0].x;
      dest[2].x = dest[1].x;
    }
  }
}

static void gfs_map_lonlat_init (GfsMap * map)
{
  map->transform = map_lonlat_transform;
  map->inverse =   map_lonlat_inverse;
  map->inverse_cell = map_lonlat_inverse_cell;
  GFS_MAP_LONLAT (map)->r = 1.;
}

static GfsMapClass * gfs_map_lonlat_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_lonlat_info = {
      "GfsMapLonLat",
      sizeof (GfsMapLonLat),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_lonlat_class_init,
      (GtsObjectInitFunc) gfs_map_lonlat_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_lonlat_info);
  }

  return klass;
}

/**
 * The longitude/latitude metric.
 * \beginobject{GfsMetricLonLat}
 */

static void metric_lon_lat_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_lon_lat_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %g", GFS_METRIC_LON_LAT (o)->r);
}

static gdouble lon_lat_face_metric (const GfsDomain * domain, const FttCellFace * face)
{
  if (face->d/2 != FTT_Y)
    return 1.;
  return face->d == 2 ? 
    GFS_VALUE (face->cell, GFS_METRIC_LON_LAT (domain->metric_data)->h2) :
    GFS_VALUE (face->cell, GFS_METRIC_LON_LAT (domain->metric_data)->h3);
}

static gdouble lon_lat_cell_metric (const GfsDomain * domain, const FttCell * cell)
{
  return GFS_VALUE (cell, GFS_VARIABLE (domain->metric_data));
}

static gdouble lon_lat_scale_metric (const GfsDomain * domain, const FttCell * cell, FttComponent c)
{
  if (c != FTT_X)
    return 1.;
  return GFS_VALUE (cell, GFS_VARIABLE (domain->metric_data));
}

static gdouble lon_lat_face_scale_metric (const GfsDomain * domain, const FttCellFace * face, 
					  FttComponent c)
{
  if (c != FTT_X)
    return 1.;
  return gfs_face_interpolated_value (face, GFS_VARIABLE (domain->metric_data)->i);
}

static void lonlat_coarse_fine (FttCell * parent, GfsVariable * a)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  GfsMetricLonLat * lonlat = GFS_METRIC_LON_LAT (a);
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  FttVector p;
  ftt_cell_pos (parent, &p);
  double theta0 = gfs_object_simulation (lonlat)->physical_params.L/lonlat->r;
  double theta = p.y*theta0;
  double h = ftt_cell_size (parent);
  double dtheta = h*theta0/2.;
  double theta1 = theta + dtheta;
  double theta2 = theta - dtheta;
  double sintheta = sin (theta);

  GFS_VALUE (child.c[0], a) = GFS_VALUE (child.c[1], a) = (sin (theta1) - sintheta)/dtheta;
  GFS_VALUE (child.c[2], a) = GFS_VALUE (child.c[3], a) = (sintheta - sin (theta2))/dtheta;

  GFS_VALUE (child.c[0], lonlat->h2) = GFS_VALUE (child.c[1], lonlat->h2) = cos (theta1);
  GFS_VALUE (child.c[0], lonlat->h3) = GFS_VALUE (child.c[1], lonlat->h3) = 
    GFS_VALUE (child.c[2], lonlat->h2) = GFS_VALUE (child.c[3], lonlat->h2) = 
    cos (theta);
  GFS_VALUE (child.c[2], lonlat->h3) = GFS_VALUE (child.c[3], lonlat->h3) = cos (theta2);
}

static void lonlat_fine_coarse (FttCell * parent, GfsVariable * a)
{
  GfsMetricLonLat * lonlat = GFS_METRIC_LON_LAT (a);
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  gdouble va = 0.;
  for (n = 0; n < 4; n++)
    /* fixme: won't work with solid boundaries */
    va += GFS_VALUE (child.c[n], a);
  GFS_VALUE (parent, a) = va/4.;

  GFS_VALUE (parent, lonlat->h2) = (GFS_VALUE (child.c[0], lonlat->h2) +
				    GFS_VALUE (child.c[1], lonlat->h2))/2.;
  GFS_VALUE (parent, lonlat->h3) = (GFS_VALUE (child.c[3], lonlat->h3) +
				    GFS_VALUE (child.c[2], lonlat->h3))/2.;
}

static void metric_lon_lat_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_lon_lat_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_METRIC_LON_LAT (*o)->r = gfs_read_constant (fp, gfs_object_simulation (*o));
  if (fp->type == GTS_ERROR)
    return;
  if (GFS_METRIC_LON_LAT (*o)->r <= 0.) {
    gts_file_error (fp, "radius must be strictly positive");
    return;
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsVariable * a = GFS_VARIABLE (*o);
  GfsMetricLonLat * lonlat = GFS_METRIC_LON_LAT (a);
  gchar * name = g_strdup_printf ("%sh2", a->name);
  lonlat->h2 = gfs_domain_get_or_add_variable (domain, name, "LonLat face metric");
  lonlat->h2->coarse_fine = lonlat->h2->fine_coarse = none;
  g_free (name);
  name = g_strdup_printf ("%sh3", a->name);
  lonlat->h3 = gfs_domain_get_or_add_variable (domain, name, "LonLat face metric");
  lonlat->h3->coarse_fine = lonlat->h3->fine_coarse = none;
  g_free (name);
  g_free (a->description);
  a->description = g_strdup ("LonLat cell metric");
  a->coarse_fine = lonlat_coarse_fine;
  a->fine_coarse = lonlat_fine_coarse;

  GtsObject * map = gts_object_new (GTS_OBJECT_CLASS (gfs_map_lonlat_class ()));
  gfs_object_simulation_set (map, domain);
  gts_container_add (GTS_CONTAINER (GFS_SIMULATION (domain)->maps), GTS_CONTAINEE (map));
  GFS_MAP_LONLAT (map)->r = GFS_METRIC_LON_LAT (*o)->r;

  domain->metric_data = *o;
  domain->face_metric  = lon_lat_face_metric;
  domain->cell_metric  = lon_lat_cell_metric;
  domain->solid_metric = solid_metric;
  domain->scale_metric = lon_lat_scale_metric;
  domain->face_scale_metric = lon_lat_face_scale_metric;
}

static void metric_lon_lat_class_init (GtsObjectClass * klass)
{
  klass->read = metric_lon_lat_read;
  klass->write = metric_lon_lat_write;
}

static void metric_lon_lat_init (GfsMetricLonLat * m)
{
  m->r = 1.;
}

GfsVariableClass * gfs_metric_lon_lat_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_metric_lon_lat_info = {
      "GfsMetricLonLat",
      sizeof (GfsMetricLonLat),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) metric_lon_lat_class_init,
      (GtsObjectInitFunc) metric_lon_lat_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_metric_class ()), 
				  &gfs_metric_lon_lat_info);
  }

  return klass;
}

/** \endobject{GfsMetricLonLat} */

/* GfsMapStretch: Header */

#define GFS_IS_MAP_STRETCH(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_stretch_class ()))

static GfsMapClass * gfs_map_stretch_class      (void);

static void gfs_map_stretch_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricStretch */
}

static void gfs_map_stretch_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricStretch */
}

static void gfs_map_stretch_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_stretch_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_stretch_write;
}

static void map_stretch_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMetricStretch * m = GFS_DOMAIN (gfs_object_simulation (map))->metric_data;
  dest->x = src->x/m->sx;
  dest->y = src->y/m->sy;
#if !FTT_2D
  dest->z = src->z/m->sz;
#endif
}

static void map_stretch_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMetricStretch * m = GFS_DOMAIN (gfs_object_simulation (map))->metric_data;
  dest->x = src->x*m->sx;
  dest->y = src->y*m->sy;
#if !FTT_2D
  dest->z = src->z*m->sz;
#endif
}

static void gfs_map_stretch_init (GfsMap * map)
{
  map->transform = map_stretch_transform;
  map->inverse =   map_stretch_inverse;
}

static GfsMapClass * gfs_map_stretch_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_stretch_info = {
      "GfsMapStretch",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_stretch_class_init,
      (GtsObjectInitFunc) gfs_map_stretch_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_stretch_info);
  }

  return klass;
}

/**
 * The "stretch" metric.
 * \beginobject{GfsMetricStretch}
 */

static void metric_stretch_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_stretch_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { sx = %g sy = %g sz = %g }",
	   GFS_METRIC_STRETCH (o)->sx, 
	   GFS_METRIC_STRETCH (o)->sy,
	   GFS_METRIC_STRETCH (o)->sz);
}

static gdouble stretch_face_metric (const GfsDomain * domain, const FttCellFace * face)
{ 
  GfsMetricStretch * s = GFS_METRIC_STRETCH (domain->metric_data);
  switch (face->d/2) {
#if FTT_2D
  case FTT_X: return s->sy;
  case FTT_Y: return s->sx;
#else
  case FTT_X: return s->sy*s->sz;
  case FTT_Y: return s->sx*s->sz;
  case FTT_Z: return s->sx*s->sy;
#endif
  default: g_assert_not_reached ();
  }
  return 1.;
}

static gdouble stretch_cell_metric (const GfsDomain * domain, const FttCell * cell)
{
  GfsMetricStretch * s = GFS_METRIC_STRETCH (domain->metric_data);
#if FTT_2D
  return s->sx*s->sy;
#else
  return s->sx*s->sy*s->sz;
#endif
}

static void stretch_solid_metric (const GfsDomain * domain, const FttCell * cell, FttVector * m)
{
  GfsMetricStretch * s = GFS_METRIC_STRETCH (domain->metric_data);
  g_assert (GFS_IS_MIXED (cell));
#if FTT_2D
  m->x = s->sy/s->sx;
  m->y = s->sx/s->sy;
#else
  m->x = s->sy*s->sz/s->sx;
  m->y = s->sx*s->sz/s->sy;
  m->z = s->sx*s->sy/s->sz;
#endif
}

static gdouble stretch_scale_metric (const GfsDomain * domain, const FttCell * cell, FttComponent c)
{
  return (&GFS_METRIC_STRETCH (domain->metric_data)->sx)[c];
}

static gdouble stretch_face_scale_metric (const GfsDomain * domain, const FttCellFace * face,
					  FttComponent c)
{
  return (&GFS_METRIC_STRETCH (domain->metric_data)->sx)[c];
}

static void metric_stretch_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_stretch_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsMetricStretch * s = GFS_METRIC_STRETCH (*o);
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "sx", TRUE, &s->sx},
      {GTS_DOUBLE, "sy", TRUE, &s->sy},
      {GTS_DOUBLE, "sz", TRUE, &s->sz},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (s->sx <= 0. || s->sy <= 0. || s->sz <= 0.) {
      gts_file_error (fp, "stretching factors must be strictly positive");
      return;
    }
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GtsObject * map = gts_object_new (GTS_OBJECT_CLASS (gfs_map_stretch_class ()));
  gfs_object_simulation_set (map, domain);
  gts_container_add (GTS_CONTAINER (GFS_SIMULATION (domain)->maps), GTS_CONTAINEE (map));

  domain->metric_data = *o;
  domain->face_metric  = stretch_face_metric;
  domain->cell_metric  = stretch_cell_metric;
  domain->solid_metric = stretch_solid_metric;
  domain->scale_metric = stretch_scale_metric;
  domain->face_scale_metric = stretch_face_scale_metric;
}

static void metric_stretch_class_init (GtsObjectClass * klass)
{
  klass->read = metric_stretch_read;
  klass->write = metric_stretch_write;
}

static void metric_stretch_init (GfsMetricStretch * m)
{
  GFS_EVENT (m)->istep = 1;
  m->sx = m->sy = m->sz = 1.;
}

GfsEventClass * gfs_metric_stretch_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMetricStretch",
      sizeof (GfsMetricStretch),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) metric_stretch_class_init,
      (GtsObjectInitFunc) metric_stretch_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_metric_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsMetricStretch} */

/* GfsMetricCubed1 is a reimplementation of GfsMetricCubed using
   GfsStoredMetric. This is left here as an example of how to
   implement a complex metric relatively simply. */

/* GfsMapCubed1: Header */

#define GFS_IS_MAP_CUBED1(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_cubed1_class ()))

static GfsMapClass * gfs_map_cubed1_class      (void);

/* GfsMapCubed1: Object */

static void gfs_map_cubed1_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed1 */
}

static void gfs_map_cubed1_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed1 */
}

static void gfs_map_cubed1_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_cubed1_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_cubed1_write;
}

static void map_cubed1_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  cmap_xy2XYZ (src->x, src->y, &dest->x, &dest->y, &dest->z);
  dest->x *= 2./M_PI;
  dest->y *= 2./M_PI; 
  dest->z *= 2./M_PI;
}

static void gfs_map_cubed1_init (GfsMap * map)
{
  map->inverse = map_cubed1_inverse;
}

static GfsMapClass * gfs_map_cubed1_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMapCubed1",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_cubed1_class_init,
      (GtsObjectInitFunc) gfs_map_cubed1_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &info);
  }

  return klass;
}

/* GfsMetricCubed1: Object */

static void metric_cubed1_init (GfsStoredMetric * m)
{
  m->map_class = gfs_map_cubed1_class ();
}

GfsVariableClass * gfs_metric_cubed1_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMetricCubed1",
      sizeof (GfsStoredMetric),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) metric_cubed1_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_stored_metric_class ()), &info);
  }

  return klass;
}

/* GfsMapVariable: Header */

#define GFS_IS_MAP_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_variable_class ()))

static GfsMapClass * gfs_map_variable_class      (void);

/* GfsMapVariable: Object */

static void gfs_map_variable_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricVariable */
}

static void gfs_map_variable_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricVariable */
}

static void gfs_map_variable_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_variable_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_variable_write;
}

static void map_variable_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (map));
  FttVector p = *src;
  FttCell * cell = gfs_domain_boundary_locate (domain, p, -1, NULL);
  if (cell == NULL)
    fprintf (stderr, "%g %g\n", p.x, p.y);
  g_return_if_fail (cell != NULL);
  GfsMetricVariable * m = domain->metric_data;
  FttComponent c;
  for (c = 0; c < 2; c++)
    (&dest->x)[c] = gfs_interpolate (cell, p, m->x[c]);
  dest->z = m->x[2] ? gfs_interpolate (cell, p, m->x[2]) : src->z;
}

static void gfs_map_variable_init (GfsMap * map)
{
  map->inverse = map_variable_inverse;
}

static GfsMapClass * gfs_map_variable_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMapVariable",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_variable_class_init,
      (GtsObjectInitFunc) gfs_map_variable_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &info);
  }

  return klass;
}

/* GfsMetricVariable: Object */

static void metric_variable_coarse_fine (FttCell * parent, GfsVariable * a)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  GfsStoredMetric * m = GFS_STORED_METRIC (a);
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  int i;
  for (i = 0; i < FTT_CELLS; i++) {
    FttComponent d;
    for (d = 0; d < 4; d++)
      GFS_VALUE (child.c[i], m->h[d]) = 1.;
    GFS_VALUE (child.c[i], a) = 1.;
  }
}

static void metric_variable_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_variable_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsMetricVariable * m = GFS_METRIC_VARIABLE (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (m));
  GfsVariable * v = GFS_VARIABLE (*o);
  FttComponent c;
  for (c = 0; c < 3; c++) {
    gchar cname[3][2] = {"x", "y", "z"};
    gchar * name = g_strdup_printf ("%s%s", v->name, cname[c]);
    m->x[c] = gfs_domain_get_or_add_variable (domain, name, "coordinate for variable metric");
    g_free (name);
  }

  v->coarse_fine = metric_variable_coarse_fine;
}

typedef struct {
  GfsMetricVariable * m;
  gdouble max;
  GfsVariable * div;
} UpdateData;

static void update_metric (FttCell * cell, UpdateData * p)
{
  gdouble h = ftt_cell_size (cell);
  gdouble x[4*(FTT_DIMENSION - 1) + 1], y[4*(FTT_DIMENSION - 1) + 1], dx, dy, dz;
  gfs_cell_corner_values (cell, p->m->x[0], -1, x);
  gfs_cell_corner_values (cell, p->m->x[1], -1, y);
  gdouble z[4*(FTT_DIMENSION - 1) + 1];
  if (p->m->x[2])
    gfs_cell_corner_values (cell, p->m->x[2], -1, z);
  else
    z[0] = z[1] = z[2] = z[3] = 0.;

  GfsStoredMetric * m = GFS_STORED_METRIC (p->m);
  dx = x[1] - x[2]; dy = y[1] - y[2]; dz = z[1] - z[2]; 
  GFS_VALUE (cell, m->h[0]) = sqrt(dx*dx + dy*dy + dz*dz)/h;
  dx = x[3] - x[0]; dy = y[3] - y[0]; dz = z[3] - z[0];
  GFS_VALUE (cell, m->h[1]) = sqrt(dx*dx + dy*dy + dz*dz)/h;
  dx = x[2] - x[3]; dy = y[2] - y[3]; dz = z[2] - z[3];
  GFS_VALUE (cell, m->h[2]) = sqrt(dx*dx + dy*dy + dz*dz)/h;
  dx = x[0] - x[1]; dy = y[0] - y[1]; dz = z[0] - z[1];
  GFS_VALUE (cell, m->h[3]) = sqrt(dx*dx + dy*dy + dz*dz)/h;

  GtsVector v02, v01, v03;
  v01[0] = x[1] - x[0]; v01[1] = y[1] - y[0]; v01[2] = z[1] - z[0];
  v02[0] = x[2] - x[0]; v02[1] = y[2] - y[0]; v02[2] = z[2] - z[0];
  v03[0] = x[3] - x[0]; v03[1] = y[3] - y[0]; v03[2] = z[3] - z[0];
  GtsVector c012, c023;
  gts_vector_cross (c012, v01, v02);
  gts_vector_cross (c023, v02, v03);
  gdouble am = (gts_vector_norm (c012) + gts_vector_norm (c023))/(2.*h*h);
  gdouble change = fabs (GFS_VALUE (cell, GFS_VARIABLE (m)) - am)/am;
  if (change > p->max)
    p->max = change;
  GFS_VALUE (cell, GFS_VARIABLE (m)) = am;
}

static gboolean metric_variable_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_stored_metric_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsMetricVariable * m = GFS_METRIC_VARIABLE (event);
    UpdateData p = { m, 1e6 };
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) update_metric, &p);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) GFS_VARIABLE (event)->fine_coarse, event);
    return TRUE;
  }
  return FALSE;
}

static void metric_variable_class_init (GtsObjectClass * klass)
{
  klass->read = metric_variable_read;
  GFS_EVENT_CLASS (klass)->event = metric_variable_event;
}

static void metric_variable_init (GfsStoredMetric * m)
{
  GFS_EVENT (m)->start = 0.;
  GFS_EVENT (m)->istep = G_MAXINT/2;
  m->map_class = gfs_map_variable_class ();
}

GfsVariableClass * gfs_metric_variable_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMetricVariable",
      sizeof (GfsMetricVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) metric_variable_class_init,
      (GtsObjectInitFunc) metric_variable_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_stored_metric_class ()), &info);
  }

  return klass;
}

/* GfsMetricLaplace: Object */

static gdouble conformal_face_scale_metric (const GfsDomain * domain, const FttCellFace * face,
					    FttComponent c)
{
  if (c > FTT_Y)
    return 1.;
  return GFS_VALUE (face->cell, GFS_STORED_METRIC (domain->metric_data)->h[face->d]);
}

static void metric_laplace_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_laplace_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsMetricLaplace * m = GFS_METRIC_LAPLACE (*o);
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "conformal", TRUE, &m->conformal},
      {GTS_INT, "spherical", TRUE, &m->spherical},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
  }

  GfsSimulation * sim = gfs_object_simulation (*o);
  if (!m->spherical) {
    gts_object_destroy (GTS_OBJECT (GFS_METRIC_VARIABLE (m)->x[2]));
    GFS_METRIC_VARIABLE (m)->x[2] = NULL;
  }
  if (m->conformal)
    GFS_DOMAIN (sim)->face_scale_metric = conformal_face_scale_metric;
}

static void metric_laplace_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_laplace_class ())->parent_class->write) (o, fp);
  GfsMetricLaplace * m = GFS_METRIC_LAPLACE (o);
  fprintf (fp, " { conformal = %d spherical = %d }", m->conformal, m->spherical);
}

static void cartesian_coordinates (FttCell * cell, UpdateData * d)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  GFS_VALUE (cell, d->m->x[0]) = p.x;
  GFS_VALUE (cell, d->m->x[1]) = p.y;
  if (d->m->x[2])
    GFS_VALUE (cell, d->m->x[2]) = p.z;
  GFS_VALUE (cell, d->div) = 0.;
}

static void spherical_laplacian (FttCell * cell, GfsVariable * dia)
{
  gdouble h = ftt_cell_size (cell);
  GFS_VALUE (cell, dia) = -2.*h*h*gfs_domain_cell_fraction (dia->domain, cell);
}

static gboolean metric_laplace_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_stored_metric_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * dia = gfs_temporary_variable (domain);
    GfsVariable * div = gfs_temporary_variable (domain);
    GfsVariable * res = gfs_temporary_variable (domain);
    GfsMetricLaplace * m = GFS_METRIC_LAPLACE (event);
    GfsMetricVariable * mv = GFS_METRIC_VARIABLE (event);

    GfsMultilevelParams par;
    gfs_multilevel_params_init (&par);
    /* convergence is dominated by non-linear terms, 1 iteration of the multigrid is enough */
    par.tolerance = 1.;
    par.nitermin = 1;

    UpdateData p = { mv, 1e6, div };
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) cartesian_coordinates, &p);
    if (!m->spherical)
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) gfs_cell_reset, dia);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, mv->x[0]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, mv->x[1]);
    if (mv->x[2])
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, mv->x[2]);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) update_metric, &p);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) GFS_VARIABLE (event)->fine_coarse, event);

    int niter = 1000;
    while (niter-- && p.max > 1e-3) {
      gfs_poisson_coefficients (domain, NULL, FALSE, TRUE, TRUE);
      if (m->spherical)
	gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				  (FttCellTraverseFunc) spherical_laplacian, dia);
      par.poisson_solve (domain, &par, mv->x[0], div, res, dia, 1.);
      par.poisson_solve (domain, &par, mv->x[1], div, res, dia, 1.);
      if (m->spherical)
	par.poisson_solve (domain, &par, mv->x[2], div, res, dia, 1.);
      p.max = 0.;
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) update_metric, &p);
      gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				(FttCellTraverseFunc) GFS_VARIABLE (event)->fine_coarse, event);
      //      gfs_multilevel_params_stats_write (&par, stderr);
      //      fprintf (stderr, "%d %g %g\n", niter, p.max, asym);
    }

    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) GFS_VARIABLE (event)->fine_coarse, event);

    gts_object_destroy (GTS_OBJECT (dia));
    gts_object_destroy (GTS_OBJECT (div));
    gts_object_destroy (GTS_OBJECT (res));
    return TRUE;
  }
  return FALSE;
}

static void metric_laplace_class_init (GtsObjectClass * klass)
{
  klass->read = metric_laplace_read;
  klass->write = metric_laplace_write;
  GFS_EVENT_CLASS (klass)->event = metric_laplace_event;
}

GfsVariableClass * gfs_metric_laplace_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsMetricLaplace",
      sizeof (GfsMetricLaplace),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) metric_laplace_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_metric_variable_class ()), &info);
  }

  return klass;
}
