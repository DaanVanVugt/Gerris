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
 * \brief Coordinates transformations.
 */

#include "map.h"

/**
 * Coordinates transformations.
 * \beginobject{GfsMap}
 */

static void gfs_map_read (GtsObject ** o, GtsFile * fp)
{
  GfsMap * map = GFS_MAP (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsMapClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_map_class ())) {
    gts_file_error (fp, "`%s' is not a GfsMap", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (map));
    map = GFS_MAP (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);
}

static void gfs_map_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "  %s", o->klass->info.name);
}

static void gfs_map_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_write;
}

static void not_implemented (GfsMap * map, const FttVector * src, FttVector * dest)
{
  g_assert_not_implemented ();
}

static double evaluate (GfsMap * map, const FttVector * x, const FttVector * rhs, FttVector * f)
{
  gdouble delta = 0.;
  (* map->inverse) (map, x, f);
  int i;
  for (i = 0; i < 3; i++) {
    (&f->x)[i] -= (&rhs->x)[i];
    delta += (&f->x)[i]*(&f->x)[i];
  }
  return delta;
}

#define DELTA 1e-6
#define NMAX 100

static void jacobian (GfsMap * map, const FttVector * x, const FttVector * rhs, FttVector * f,
		      GtsMatrix * J)
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      FttVector df, dx = *x;
      (&dx.x)[j] += DELTA;
      (* map->inverse) (map, &dx, &df);
      J[i][j] = ((&df.x)[i] - (&rhs->x)[i] - (&f->x)[i])/DELTA;
    }
}

static void map_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  FttVector f, rhs = *src;
  GtsMatrix J[4];
  int n = 0;
  /* use multidimensional Newton iterations to invert map(dest) = src */
  *dest = *src; /* use src as initial guess */
  gdouble delta = evaluate (map, dest, &rhs, &f);
  while (delta > 1e-12 && n < NMAX) {
    jacobian (map, dest, &rhs, &f, J);
    GtsMatrix * iJ = gts_matrix3_inverse (J);
    if (!iJ) {
      gts_matrix_print (J, stderr);
      g_assert_not_reached ();
    }
    int i, j;
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
	(&dest->x)[i] -= iJ[i][j]*(&f.x)[j];
    gts_matrix_destroy (iJ);
    delta = evaluate (map, dest, &rhs, &f);
    n++;
  }
  g_assert (n < NMAX);
}

static void normalized_jacobian (GfsMap * map, const FttVector * p, GtsMatrix * J)
{
  FttVector f, rhs = {0., 0., 0.};
  evaluate (map, p, &rhs, &f);
  jacobian (map, p, &rhs, &f, J);
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

static void map_inverse_vector (GfsMap * map, const FttVector * p,
				const FttVector * src, FttVector * dest)
{
  GtsMatrix J[4];
  normalized_jacobian (map, p, J);
  FttVector src1 = *src; /* in case src and dest are identical */
  int i, j;
  for (i = 0; i < 3; i++) {
    (&dest->x)[i] = 0.;
    for (j = 0; j < 3; j++)
      (&dest->x)[i] += (&src1.x)[j]*J[i][j];
  }
}

static void map_transform_vector (GfsMap * map, const FttVector * p,
				  const FttVector * src, FttVector * dest)
{
  GtsMatrix J[4];
  normalized_jacobian (map, p, J);
  GtsMatrix * iJ = gts_matrix3_inverse (J);
  if (!iJ) {
    gts_matrix_print (J, stderr);
    g_assert_not_reached ();
  }
  FttVector src1 = *src; /* in case src and dest are identical */
  int i, j;
  for (i = 0; i < 3; i++) {
    (&dest->x)[i] = 0.;
    for (j = 0; j < 3; j++)
      (&dest->x)[i] += (&src1.x)[j]*iJ[i][j];
  }
  gts_matrix_destroy (iJ);
}

static void inverse_cell (GfsMap * map, const FttVector * src, FttVector * dest)
{
  gint i;
  for (i = 0; i < 4; i++)
    (* map->inverse) (map, &(src[i]), &(dest[i]));
}

static void gfs_map_init (GfsMap * map)
{
  map->inverse = not_implemented;
  map->transform = map_transform;
  map->inverse_vector = map_inverse_vector;
  map->transform_vector = map_transform_vector;
  map->inverse_cell = inverse_cell;
}

GfsMapClass * gfs_map_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_info = {
      "GfsMap",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_class_init,
      (GtsObjectInitFunc) gfs_map_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_map_info);
  }

  return klass;
}

/** \endobject{GfsMap} */

/**
 *
 * \beginobject{GfsMapFunction}
 */

static void gfs_map_function_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting a parameter block");
    return;
  }

  GfsMapFunction * m = GFS_MAP_FUNCTION (*o);
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

static void gfs_map_function_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->write) (o, fp);
 
  GfsMapFunction * m = GFS_MAP_FUNCTION (o);
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

static void gfs_map_function_destroy (GtsObject * o)
{
  GfsMapFunction * m = GFS_MAP_FUNCTION (o);
  FttComponent c;
  for (c = 0; c < 3; c++)
    if ((&m->x)[c])
      gts_object_destroy (GTS_OBJECT ((&m->x)[c]));

  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->destroy) (o);
}

static void gfs_map_function_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_function_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_function_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_map_function_destroy;
}

static void map_function_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMapFunction * m = GFS_MAP_FUNCTION (map);
  gdouble L = gfs_object_simulation (m)->physical_params.L;
  FttVector src1;
  FttComponent c;
  for (c = 0; c < 3; c++)
    (&src1.x)[c] = (&src->x)[c]*L;
  for (c = 0; c < 3; c++)
    if ((&m->x)[c])
      (&dest->x)[c] = gfs_function_spatial_value ((&m->x)[c], &src1)/L;
    else
      (&dest->x)[c] = (&src->x)[c];
}

static void gfs_map_function_init (GfsMapFunction * m)
{
  GFS_MAP (m)->inverse = map_function_inverse;
  m->x = gfs_function_new (gfs_function_map_class (), 1.);
  m->y = gfs_function_new (gfs_function_map_class (), 1.);
  m->z = gfs_function_new (gfs_function_map_class (), 1.);
}

GfsMapClass * gfs_map_function_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_function_info = {
      "GfsMapFunction",
      sizeof (GfsMapFunction),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_function_class_init,
      (GtsObjectInitFunc) gfs_map_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_function_info);
  }

  return klass;
}

/** \endobject{GfsMapFunction} */

/**
 * Isometric coordinates transformations.
 * \beginobject{GfsMapTransform}
 */

static void gfs_map_transform_destroy (GtsObject * o)
{
  gts_matrix_destroy (GFS_MAP_TRANSFORM (o)->m);
  gts_matrix_destroy (GFS_MAP_TRANSFORM (o)->im);

  (* GTS_OBJECT_CLASS (gfs_map_transform_class ())->parent_class->destroy) (o);
}

static void gfs_map_transform_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_transform_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsMapTransform * map = GFS_MAP_TRANSFORM (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "tx", TRUE, &map->translate[0]},
    {GTS_DOUBLE, "ty", TRUE, &map->translate[1]},
    {GTS_DOUBLE, "tz", TRUE, &map->translate[2]},
    {GTS_DOUBLE, "rx", TRUE, &map->rotate[0]},
    {GTS_DOUBLE, "ry", TRUE, &map->rotate[1]},
    {GTS_DOUBLE, "rz", TRUE, &map->rotate[2]},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);

  if (map->rotate[0] != 0.) { /* rotate around x-axis */
    gdouble angle = map->rotate[0]*M_PI/180.;
    gdouble cosa = cos (angle), sina = sin (angle);
    GtsMatrix * rot = gts_matrix_identity (NULL);
    rot[1][1] = cosa; rot[1][2] = -sina;
    rot[2][1] = sina; rot[2][2] = cosa;
    GtsMatrix * p = gts_matrix_product (map->m, rot);
    gts_matrix_destroy (rot);
    gts_matrix_destroy (map->m);
    map->m = p;
  }
  if (map->rotate[1] != 0.) { /* rotate around y-axis */
    gdouble angle = map->rotate[1]*M_PI/180.;
    gdouble cosa = cos (angle), sina = sin (angle);
    GtsMatrix * rot = gts_matrix_identity (NULL);
    rot[0][0] = cosa; rot[0][2] = sina;
    rot[2][0] = -sina; rot[2][2] = cosa;
    GtsMatrix * p = gts_matrix_product (map->m, rot);
    gts_matrix_destroy (rot);
    gts_matrix_destroy (map->m);
    map->m = p;
  }
  if (map->rotate[2] != 0.) { /* rotate around z-axis */
    gdouble angle = map->rotate[2]*M_PI/180.;
    gdouble cosa = cos (angle), sina = sin (angle);
    GtsMatrix * rot = gts_matrix_identity (NULL);
    rot[0][0] = cosa; rot[0][1] = -sina;
    rot[1][0] = sina; rot[1][1] = cosa;
    GtsMatrix * p = gts_matrix_product (map->m, rot);
    gts_matrix_destroy (rot);
    gts_matrix_destroy (map->m);
    map->m = p;
  }

  gdouble L = gfs_object_simulation (map)->physical_params.L;
  map->m[0][3] += map->translate[0]/L;
  map->m[1][3] += map->translate[1]/L;
  map->m[2][3] += map->translate[2]/L;

  gts_matrix_destroy (map->im);
  map->im = gts_matrix_inverse (map->m);
}

static void gfs_map_transform_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_transform_class ())->parent_class->write) (o, fp);
  GfsMapTransform * map = GFS_MAP_TRANSFORM (o);
  fputs (" {\n", fp);
  if (gts_vector_norm (map->translate) > 0.)
    fprintf (fp, "  tx = %g ty = %g tz = %g\n",
	     map->translate[0], map->translate[1], map->translate[2]);
  if (gts_vector_norm (map->rotate) > 0.)
    fprintf (fp, "  rx = %g ry = %g rz = %g\n",
	     map->rotate[0], map->rotate[1], map->rotate[2]);
  fputc ('}', fp);
}

static void gfs_map_transform_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_map_transform_destroy;
  klass->read = gfs_map_transform_read;
  klass->write = gfs_map_transform_write;
}

static void map_transform_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GtsPoint p;
  p.x = src->x; p.y = src->y; p.z = src->z;
  gts_point_transform (&p, GFS_MAP_TRANSFORM (map)->im);
  dest->x = p.x; dest->y = p.y; dest->z = p.z;
}

static void map_transform_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GtsPoint p;
  p.x = src->x; p.y = src->y; p.z = src->z;
  gts_point_transform (&p, GFS_MAP_TRANSFORM (map)->m);
  dest->x = p.x; dest->y = p.y; dest->z = p.z;
}

static void gfs_map_transform_init (GfsMap * map)
{
  map->transform = map_transform_transform;
  map->inverse =   map_transform_inverse;
  GFS_MAP_TRANSFORM (map)->m =  gts_matrix_identity (NULL);
  GFS_MAP_TRANSFORM (map)->im = gts_matrix_identity (NULL);
}

GfsMapClass * gfs_map_transform_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_transform_info = {
      "GfsMapTransform",
      sizeof (GfsMapTransform),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_transform_class_init,
      (GtsObjectInitFunc) gfs_map_transform_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_transform_info);
  }

  return klass;
}

/** \endobject{GfsMapTransform} */

