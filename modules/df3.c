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

#include <stdlib.h>
#include "simulation.h"
#include "map.h"
#include "output.h"
#include "init.h"

typedef struct Vector_ {
  gint x, y, z;
} Vector;

#define POW2(x) (1 << (int)(x))

/* GfsOutputPovrayDF3: Header */

GfsOutputClass * gfs_output_povray_DF3_class (void);

#define GFS_IS_OUTPUT_POVRAY_DF3(obj) (gts_object_is_from_class (obj, \
								 gfs_output_povray_DF3_class ()))

static void write_density_value (FttCell * cell, gpointer * data)
{
  unsigned char * buf = data[0];
  gdouble * min = data[1];
  gdouble * max = data[2];
  GfsVariable * v = data[3];
  guint * min_depth = data[4];
  guint * max_depth = data[5];
  Vector * dimensions = (Vector *) data[6];
  guint xc = dimensions->x;
  guint yc = dimensions->y;
  guint zc = dimensions->z;
  guint position_size = *((guint *) data[7]);
  GtsBBox * bb_density = data[8];

  FttVector p;
  guint pos = 0;
  gdouble value;
  guint value_i;
  guint level = ftt_cell_level(cell);

  if (level > *max_depth || level < *min_depth) {
    return;
  }

  ftt_cell_pos (cell, &p);

  gint xp = ((gdouble)(p.x - bb_density->x1))*POW2(*max_depth);
  gint yp = ((gdouble)(p.y - bb_density->y1))*POW2(*max_depth);
  gint zp = ((gdouble)(p.z - bb_density->z1))*POW2(*max_depth);

  pos = zp*(yc*xc) + yp*xc + xp;
  if (pos > xc*yc*zc) {
    return;
  }

  value = (GFS_VALUE (cell, v) - *min)/(*max - *min);

  if (level < *max_depth) {
    gdouble size = ftt_cell_size(cell)/2.0;
    gdouble x, y, z;
    int n = POW2(*max_depth - level);
    double d = ftt_cell_size(cell)/n;
    gdouble xmin = p.x - size;
    gdouble ymin = p.y - size;
    gdouble zmin = p.z - size;
    int i, j, k;
    gdouble t;
    FttVector p1;
    gint pos1;

    for (i = 0; i < n; i++) {
      x = xmin + (i + 0.5) * d;
      if (x < bb_density->x1 || x > bb_density->x2) {
	continue;
      }
      for (j = 0; j < n; j++) {
	y = ymin + (j + 0.5) * d;
	if (y < bb_density->y1 || y > bb_density->y2) {
	  continue;
	}
	for (k = 0; k < n; k++) {
	  z = zmin + (k + 0.5) * d;
	  if (z < bb_density->z1 || z > bb_density->z2) {
	    continue;
	  }

	  p1.x = x;
	  p1.y = y;
	  p1.z = z;

	  pos1 = pos + (k - n/2)*(xc*yc) + (j - n/2)*xc + (i - n/2);

	  if (pos1 < 0 || pos1 > xc*yc*zc) {
	     continue;
	  }

	  t = gfs_interpolate(cell, p1, v);
	  value = (t - *min)/(*max - *min);

	  if (position_size == 2) {
	    value_i = value * 0xffff;
	    buf[2*pos1] = (unsigned char) ((value_i >> 8) & 0xff);
	    buf[2*pos1 + 1] = (unsigned char) (value_i & 0xff);
	  }
	  else {
	    value_i = value * 0xff;
	    buf[pos1] = (unsigned char) (value_i & 0xff);
	  }
	}
      }
    }
  }
  else {
    if (position_size == 2) {
      value_i = value * 0xffff;
      buf[2*pos] = (unsigned char) ((value_i >> 8) & 0xff);
      buf[2*pos + 1] = (unsigned char) (value_i & 0xff);
    }
    else {
      value_i = value * 0xff;
      buf[pos] = (unsigned char) (value_i & 0xff);
    }
  }
}

static void write_density_file (FILE * file, guint xc, guint yc, guint zc,
				char * density_buf, size_t density_buf_s)
{
  unsigned char header[] = "\0\0\0\0\0\0";
  g_return_if_fail (file != NULL);
  g_return_if_fail (density_buf != NULL);
  /* header: x,y,z dimensions, two bytes each, big endian */
  header[0] = (unsigned char) ((xc & 0xff00) >> 8);
  header[1] = (unsigned char) (xc & 0xff);
  header[2] = (unsigned char) ((yc & 0xff00) >> 8);
  header[3] = (unsigned char) (yc & 0xff);
  header[4] = (unsigned char) ((zc & 0xff00) >> 8);
  header[5] = (unsigned char) (zc & 0xff);
  fwrite(header, sizeof(char), 6, file);
  /* write prepared data */
  fwrite(density_buf, sizeof(char), density_buf_s, file);
}

static void max_extent (FttCell * cell, FttVector extent[2])
{
  gdouble h = ftt_cell_size (cell)/2.;
  FttVector pos;
  
  ftt_cell_pos (cell, &pos);
  if (pos.x - h < extent[0].x) extent[0].x = pos.x - h;
  if (pos.y - h < extent[0].y) extent[0].y = pos.y - h;
  if (pos.z - h < extent[0].z) extent[0].z = pos.z - h;

  if (pos.x + h > extent[1].x) extent[1].x = pos.x + h;
  if (pos.y + h > extent[1].y) extent[1].y = pos.y + h;
  if (pos.z + h > extent[1].z) extent[1].z = pos.z + h;
}

static gboolean cell_condition (FttCell * cell, gpointer condition)
{
  return gfs_function_value (condition, cell);
}

static void gfs_write_povray_density(GfsDomain * domain, 
				     GfsFunction * condition,
				     GfsVariable * v, gdouble min, gdouble max,
				     FttTraverseFlags flags,
				     gint level,
				     FILE * fp)
{
  FttVector extent[2] = {{ G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE },
			 { - G_MAXDOUBLE, - G_MAXDOUBLE, - G_MAXDOUBLE }};
  Vector dimensions;
  GtsBBox bb_density;
  gpointer data[9];
  /* Povray would support 1, 2, 4, but here we use only 1 byte
   * per voxel. */
  guint position_size = 1;
  guint max_depth = 0;
  guint min_depth = 0;
  guint xc;
  guint yc;
  guint zc;
  size_t density_buf_s;
  char *density_buf = NULL;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  if (min == max) {
    max = min + 1.;
  }

  if (level < 0) {
    max_depth = gfs_domain_depth(domain);
  }
  else {
    max_depth = level;
  }

  if (condition) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, level,
					(FttCellTraverseFunc) max_extent, extent,
					cell_condition, condition);
    gfs_restore_fpe_for_function (condition);
  }
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			      (FttCellTraverseFunc) max_extent, extent);
    
  if (extent[0].x == G_MAXDOUBLE)
    return;

  bb_density.x1 = extent[0].x;
  bb_density.y1 = extent[0].y;
  bb_density.z1 = extent[0].z;
  bb_density.x2 = extent[1].x;
  bb_density.y2 = extent[1].y;
  bb_density.z2 = extent[1].z;

  /* total memory size */
  xc = (extent[1].x - extent[0].x)*POW2(max_depth);
  yc = (extent[1].y - extent[0].y)*POW2(max_depth);
  zc = (extent[1].z - extent[0].z)*POW2(max_depth);

  dimensions.x = xc;
  dimensions.y = yc;
  dimensions.z = zc;

  /*
  if (verbose == TRUE) {
    fprintf (stdout, "level: %d\n", level);
    fprintf (stdout, "max_depth: %d\n", max_depth);
    fprintf (stdout, "xc, yc, zc: <%d, %d, %d>\n", xc, yc, zc);
  }
  */

  density_buf_s = position_size * xc * yc * zc;
  density_buf = g_malloc (density_buf_s);

  if (density_buf == NULL) {
    g_warning ("GfsOutputPovrayDF3: Failed to allocate %ld bytes of memory",
	       (long) density_buf_s);
    return; /* failure */
  }

  memset(density_buf, 0, density_buf_s);

  data[0] = density_buf;
  data[1] = &min;
  data[2] = &max;
  data[3] = v;
  data[4] = &min_depth;
  data[5] = &max_depth;
  data[6] = &dimensions;
  data[7] = &position_size;
  data[8] = &bb_density;

  if (condition) {
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, level,
					(FttCellTraverseFunc) write_density_value, data,
					cell_condition, condition);
    gfs_restore_fpe_for_function (condition);
  }
  else {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, level,
			      (FttCellTraverseFunc) write_density_value, data);
  }

  /* write buffer to DF3 file */
  write_density_file (fp, xc, yc, zc, density_buf, density_buf_s);

  g_free(density_buf);
}

/* ------------------------------------------------------ */
/* GfsOutputPovrayDF3: Object */

static gboolean gfs_output_povray_density_event (GfsEvent * event,
						 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS
			  (gfs_output_povray_DF3_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar *output = GFS_OUTPUT_SCALAR (event);
    GfsDomain *domain = GFS_DOMAIN (sim);
    
    gfs_write_povray_density (domain,
			      output->condition,
			      output->v, output->min, output->max,
			      FTT_TRAVERSE_LEAFS, output->maxlevel,
			      GFS_OUTPUT (event)->file->fp);

    fflush (GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_povray_DF3_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_povray_density_event;
}

GfsOutputClass * gfs_output_povray_DF3_class (void)
{
  static GfsOutputClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_povray_density_info = {
      "GfsOutputPovrayDF3",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_povray_DF3_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS
				  (gfs_output_scalar_class ()),
				  &gfs_output_povray_density_info);
  }

  return klass;
}

/* Initialize module */

const gchar gfs_module_name[] = "df3";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_output_povray_DF3_class ();
  return NULL;
}
