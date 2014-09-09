/* Gerris - The GNU Flow Solver
 * Copyright (C) 2007-2012 National Institute of Water and Atmospheric Research
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
 * \brief Data defined on Cartesian grids.
 */

#include <stdlib.h>
#include "cartesian.h"

/**
 * Storage for Cartesian grid data.
 * \beginobject{GfsCartesianGrid}
 */

static void cartesian_grid_read (GtsObject ** o, GtsFile * fp)
{
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (*o);
  guint i, j, size = 1;

  if (GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  while (fp->type == '\n') 
    gts_file_next_token (fp);
  if (fp->type != GTS_INT) {
     gts_file_error (fp, "expecting an integer (N)");
     return;
  }
  cgd->N = atoi (fp->token->str);
  gts_file_next_token (fp);

  cgd->name = g_malloc0 (cgd->N*sizeof (char *));
  for (i = 0; i < cgd->N; i++) {
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (name[%d])", i);
      return;
    }
    cgd->name[i] = g_strdup (fp->token->str);
    gts_file_next_token (fp);
  }

  cgd->n = g_malloc (cgd->N*sizeof (guint));  
  for (i = 0; i < cgd->N; i++) {
    while (fp->type == '\n') 
      gts_file_next_token (fp);
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting an integer (n[%d])", i);
      return;
    }
    cgd->n[i] = atoi (fp->token->str);
    gts_file_next_token (fp);
    size *= cgd->n[i];
  }

  cgd->x = g_malloc0 (cgd->N*sizeof (gdouble *));
  for (i = 0; i < cgd->N; i++) {
    cgd->x[i] = g_malloc (cgd->n[i]*sizeof (gdouble));
    gdouble last = - G_MAXDOUBLE;
    for (j = 0; j < cgd->n[i]; j++) {
      if (fp->type == '\n')
	gts_file_next_token (fp);
      if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number (x[%d][%d])", i, j);
        return;
      }
      cgd->x[i][j] = atof (fp->token->str);
      if (cgd->x[i][j] < last) {
	gts_file_error (fp, "coordinates must be in increasing order", i, j);
        return;
      }
      last = cgd->x[i][j];
      gts_file_next_token (fp);
    }
  }

  cgd->v = g_malloc (size*sizeof (gdouble));  
  for (i = 0; i < size; i++) {
    if (fp->type == '\n')
      gts_file_next_token (fp);
    if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
      gts_file_error (fp, "expecting a number");
      return;
    }
    cgd->v[i] = atof (fp->token->str);
    gts_file_next_token (fp);
  }
}

static void cartesian_grid_write (GtsObject * o, FILE * fp)
{
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (o);
  guint i, j, size = 1;

  if (GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->write) 
      (o, fp);

  for (i = 0; i < cgd->N; i++)
    size *= cgd->n[i];

  fprintf (fp, "%d ", cgd->N);
  for (i = 0; i < cgd->N; i++)
    fprintf (fp, "%s ", cgd->name[i]);
  fputc ('\n', fp);
  for (i = 0; i < cgd->N; i++)
    fprintf (fp, "%d\n", cgd->n[i]);

  for (i = 0; i < cgd->N; i++)
    for (j = 0; j < cgd->n[i]; j++)
      fprintf (fp, "%f\n", cgd->x[i][j]);

  for (i = 0; i < size; i++)
    fprintf (fp, "%f\n", cgd->v[i]);  
}

static void cartesian_grid_destroy (GtsObject * object)
{
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (object);  

  guint i;
  if (cgd->name) {
    for (i = 0; i < cgd->N; i++)
      g_free (cgd->name[i]);
    g_free (cgd->name);
  }
  g_free (cgd->n);
  if (cgd->x) {
    for (i = 0; i < cgd->N; i++)
      g_free (cgd->x[i]);
    g_free (cgd->x);
  }
  g_free (cgd->v);
 
  (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->destroy) (object);
}

static void cartesian_grid_class_init (GtsObjectClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = cartesian_grid_read;
  GTS_OBJECT_CLASS (klass)->write = cartesian_grid_write;
  GTS_OBJECT_CLASS (klass)->destroy = cartesian_grid_destroy;
}

GtsObjectClass * gfs_cartesian_grid_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsCartesianGrid",
      sizeof (GfsCartesianGrid),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) cartesian_grid_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()), &info);
  }

  return klass;
}

GfsCartesianGrid * gfs_cartesian_grid_new (GtsObjectClass * klass)
{
  GfsCartesianGrid * object;

  object = GFS_CARTESIAN_GRID (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

static void slice (GfsCartesianGrid * g, guint p, GfsCartesianGrid * s)
{
  s->N = g->N - 1;
  s->n = &g->n[1];
  s->x = &g->x[1];
  guint i;
  gulong size = 1;
  for (i = 1; i < g->N; i++)
    size *= g->n[i];
  s->v = &g->v[size*p];
}

static gint lookup (GfsCartesianGrid * g, gdouble x)
{
  guint min = 0, max = g->n[0] - 1;
  if (x < g->x[0][min] || x > g->x[0][max])
    return -1;
  while (max > min + 1) {
    guint n = (min + max)/2;
    if (x > g->x[0][n])
      min = n;
    else
      max = n;
  }
  return min;
}

/**
 * gfs_cartesian_grid_interpolate:
 * @g: a Cartesian grid.
 * @p: a position vector of dimension @g->N.
 * @val: the interpolated value at position @p.
 *
 * Returns: %TRUE if @val has been computed, %FALSE if @p is not
 * contained within @g.
 */
gboolean gfs_cartesian_grid_interpolate (GfsCartesianGrid * g, gdouble * p, gdouble * val)
{
  g_return_val_if_fail (g != NULL, FALSE);
  g_return_val_if_fail (g->N > 0, FALSE);
  g_return_val_if_fail (p != NULL, FALSE);
  g_return_val_if_fail (val != NULL, FALSE);

  gint i = lookup (g, p[0]);
  if (i < 0)
    return FALSE;
  gdouble v1, v2;
  if (g->N > 1) {
    GfsCartesianGrid g1;
    slice (g, i, &g1);
    if (!gfs_cartesian_grid_interpolate (&g1, &p[1], &v1))
      return FALSE;
    slice (g, i + 1, &g1);
    if (!gfs_cartesian_grid_interpolate (&g1, &p[1], &v2))
      return FALSE;
  }
  else {
    v1 = g->v[i];
    v2 = g->v[i + 1];
  }

  g_assert (g->x[0][i + 1] -  g->x[0][i] != 0.);
  *val = v1 + (v2 - v1)*(p[0] - g->x[0][i])/(g->x[0][i + 1] -  g->x[0][i]);
  return TRUE;
}

/**
 * gfs_cartesian_grid_read:
 * @name: a .cgd filename.
 * @fp: a #GtsFile or %NULL.
 *
 * Returns: the #GfsCartesianGrid stored in @name, or %NULL if an
 * error occured, in which case an error message is set in @fp (if not
 * %NULL).
 */
GfsCartesianGrid * gfs_cartesian_grid_read (const gchar * name, GtsFile * fp)
{
  g_return_val_if_fail (name != NULL, NULL);

  FILE * fptr = fopen (name, "r");
  GtsFile * fp1;
  GfsCartesianGrid * grid;
  GtsObjectClass * klass;

  if (fptr == NULL) {
    if (fp)
      gts_file_error (fp, "cannot open file `%s'", name);
    return NULL;
  }

  fp1 = gts_file_new (fptr);

  klass = gfs_cartesian_grid_class ();

  grid = gfs_cartesian_grid_new (klass);
  GtsObject * o = GTS_OBJECT (grid);
  (* klass->read) (&o, fp1);

  if (fp1->type == GTS_ERROR) {
    if (fp)
      gts_file_error (fp, "%s:%d:%d: %s", name, fp1->line, fp1->pos, fp1->error);
    gts_object_destroy (GTS_OBJECT(grid));
    grid = NULL;
  }
  gts_file_destroy (fp1);
  fclose (fptr);
  return grid;
}

/** \endobject{GfsCartesianGrid} */
