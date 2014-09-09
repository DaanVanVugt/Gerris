/* Gerris - The GNU Flow Solver
 * Copyright (C) 2012 National Institute of Water and Atmospheric Research
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

#include "river.h"
#include "culvert/boyd87.h"

/* GfsSourceCulvert: Header */

typedef struct _GfsSourceCulvert         GfsSourceCulvert;

struct _GfsSourceCulvert {
  /*< private >*/
  GfsSourcePipe parent;

  /*< public >*/
  gint entrance;
  gdouble B, n, ke;
};

#define GFS_SOURCE_CULVERT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceCulvert,\
					         gfs_source_culvert_class ())
#define GFS_IS_SOURCE_CULVERT(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_culvert_class ()))

GfsSourceGenericClass * gfs_source_culvert_class  (void);

/* GfsSourceCulvert: Object */

static double box_flow_rate (double z1, double h1, /* terrain elevation and flow depth at inlet */
			     double z2, double h2, /* terrain elevation and flow depth at outlet */
			     double l,             /* pipe length */
			     double g,             /* acceleration of gravity */
			     GfsSourcePipe * p)
{
  GfsSourceCulvert * c = GFS_SOURCE_CULVERT (p);
  if (z1 + h1 > z2 + h2)
    return + Q_box (h1, h2, c->B, p->diameter, c->entrance, (z1 - z2)/l, l, c->n, c->ke, g);
  else
    /* reverse flow: assume the culvert is symmetrical */
    return - Q_box (h2, h1, c->B, p->diameter, c->entrance, (z2 - z1)/l, l, c->n, c->ke, g);
}

static double pipe_flow_rate (double z1, double h1,
			      double z2, double h2,
			      double l, double g,
			      GfsSourcePipe * p)
{
  GfsSourceCulvert * c = GFS_SOURCE_CULVERT (p);
  if (z1 + h1 > z2 + h2)
    return + Q_pipe (h1, h2, p->diameter, c->entrance, (z1 - z2)/l, l, c->n, c->ke, g);
  else
    /* reverse flow: assume the culvert is symmetrical */
    return - Q_pipe (h2, h1, p->diameter, c->entrance, (z2 - z1)/l, l, c->n, c->ke, g);
}

static void gfs_source_culvert_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_culvert_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSourceCulvert * c = GFS_SOURCE_CULVERT (*o);
  gchar * type = NULL;
  GtsFileVariable var[] = {
    {GTS_STRING, "type",     TRUE, &type},
    {GTS_INT,    "entrance", TRUE, &c->entrance},
    {GTS_DOUBLE, "B",        TRUE, &c->B},
    {GTS_DOUBLE, "n",        TRUE, &c->n},
    {GTS_DOUBLE, "ke",       TRUE, &c->ke},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (type) {
    if (!strcmp (type, "box"))
      GFS_SOURCE_PIPE (c)->flow_rate = box_flow_rate;
    else if (!strcmp (type, "pipe"))
      GFS_SOURCE_PIPE (c)->flow_rate = pipe_flow_rate;
    else {
      gts_file_variable_error (fp, var, "type", "unknown culvert type '%s'", type);
      g_free (type);
      return;
    }
    g_free (type);
  }

  if (c->entrance < 1 || c->entrance > 3)
    gts_file_variable_error (fp, var, "entrance", "entrance type must be 1,2 or 3");
  else if (GFS_SOURCE_PIPE (c)->flow_rate == pipe_flow_rate && var[2].set)
    gts_file_variable_error (fp, var, "B", "box width is irrelevant for a pipe culvert");
  else if (c->B <= 0.)
    gts_file_variable_error (fp, var, "B", "box width must be greater than zero");
  else if (c->n < 0.)
    gts_file_variable_error (fp, var, "n", "Manning coefficient must be greater than zero");
  else if (c->ke < 0.)
    gts_file_variable_error (fp, var, "ke", "entrance loss coefficient must be greater than zero");
}

static void gfs_source_culvert_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_culvert_class ())->parent_class->write) (o, fp);
  GfsSourceCulvert * c = GFS_SOURCE_CULVERT (o);
  if (GFS_SOURCE_PIPE (c)->flow_rate == box_flow_rate)
    fprintf (fp, " { type = box B = %g entrance = %d n = %g ke = %g }",
	     c->B, c->entrance, c->n, c->ke);
  else
    fprintf (fp, " { type = pipe entrance = %d n = %g ke = %g }",
	     c->entrance, c->n, c->ke);
}

static void gfs_source_culvert_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_source_culvert_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_culvert_write;
}

static void gfs_source_culvert_init (GfsSourceCulvert * c)
{
  c->entrance = 1;
  c->B = 1.;
  c->n = 0.011;
  c->ke = 0.;
  GFS_SOURCE_PIPE (c)->flow_rate = box_flow_rate;
}

GfsSourceGenericClass * gfs_source_culvert_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsSourceCulvert",
      sizeof (GfsSourceCulvert),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_culvert_class_init,
      (GtsObjectInitFunc) gfs_source_culvert_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_pipe_class ()), &info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "culvert";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_source_culvert_class ();
  return NULL;
}
