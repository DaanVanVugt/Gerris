/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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
#include <string.h>
#include <math.h>
#include <gts.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#ifndef M_PI
# define M_PI 3.14159265359
#endif

static GtsSurface * triangulate (GSList * vertices,
				 GSList * constraints)
{
  GtsVertex * v1, * v2, * v3;
  GtsSurface * s;
  GSList * i;

  v1 = gts_vertex_new (gts_vertex_class (), -1e10, -1e10, 0.);
  v2 = gts_vertex_new (gts_vertex_class (), 1e10, -1e10, 0.);
  v3 = gts_vertex_new (gts_vertex_class (), 0., 1e10, 0.);
  s = gts_surface_new (gts_surface_class (), 
		       gts_face_class (),
		       gts_edge_class (),
		       gts_vertex_class ());
  gts_surface_add_face (s, gts_face_new (gts_face_class (), 
				gts_edge_new (gts_edge_class (), v1, v2),
				gts_edge_new (gts_edge_class (), v2, v3),
				gts_edge_new (gts_edge_class (), v3, v1)));
  i = vertices;
  while (i) {
    if (gts_delaunay_add_vertex (s, i->data, NULL) != NULL) {
      gts_object_destroy (GTS_OBJECT (s));
      return NULL;
    }
    i = i->next;
  }
  
  i = constraints;
  while (i) {
    if (gts_delaunay_add_constraint (s, i->data) != NULL) {
      gts_object_destroy (GTS_OBJECT (s));
      return NULL;
    }
    i = i->next;
  }
  
  gts_delaunay_remove_hull (s);

  return s;
}

static GSList * contour (GSList * i, GtsSurface * s, 
			 gdouble z, gboolean closed)
{
  GSList * edges = NULL;
  GtsVertex * vold = NULL, * vin = NULL;

  if (i == NULL || i->next == NULL)
    return NULL;

  while (i) {
    GtsPoint * p = i->data;
    GtsVertex * v;

    v = gts_vertex_new (s->vertex_class, p->x, p->y, z);
    if (vold)
      edges = g_slist_prepend (edges, 
	  gts_edge_new (GTS_EDGE_CLASS (gts_constraint_class ()), v, vold));
    else
      vin = v;
    vold = v;
    i = i->next;
  }
  if (closed)
    edges = g_slist_prepend (edges, 
	  gts_edge_new (GTS_EDGE_CLASS (gts_constraint_class ()), vin, vold));
  return edges;
}

static void surface_add_shape (GtsSurface * s, 
			       GSList * shape,
			       gdouble z1,
			       gdouble z2,
			       guint nz,
			       gboolean closed,
			       gboolean ends)
{
  gdouble z, dz = (z2 - z1)/nz;
  guint i;
  GSList * bottom = NULL, * e1 = NULL, * e2 = NULL;

  for (i = 0, z = z1; i <= nz; i++, z += dz) {
    if (e1 == NULL)
      bottom = e1 = contour (shape, s, z, closed);
    else {
      GSList * i, * j;
      GtsEdge * eold = NULL, * ein = NULL;

      e2 = contour (shape, s, z, closed);
      i = e1; j = e2;
      while (i && j) {
	GtsEdge * e;
	GtsEdge * ee = gts_edge_new (gts_edge_class (),
				     GTS_SEGMENT (i->data)->v2,
				     GTS_SEGMENT (j->data)->v1);

	if (!closed || i->next != NULL)
	  e = gts_edge_new (gts_edge_class (),
			    GTS_SEGMENT (i->data)->v2,
			    GTS_SEGMENT (j->data)->v2);
	else
	  e = ein;
	if (eold == NULL)
	  eold = ein = gts_edge_new (gts_edge_class (), 
				     GTS_SEGMENT (i->data)->v1,
				     GTS_SEGMENT (j->data)->v1);
	gts_surface_add_face (s, gts_face_new (s->face_class, 
					       i->data, eold, ee));
	gts_surface_add_face (s, gts_face_new (s->face_class, 
					       e, ee, j->data));
	eold = e;
	i = i->next; 
	j = j->next;
      }
      if (e1 != bottom) g_slist_free (e1);
      e1 = e2;
    }
  }

  if (ends) {
    GSList * vertices1 = gts_vertices_from_segments (bottom);
    GSList * vertices2 = gts_vertices_from_segments (e1);
    GtsSurface * s1 = triangulate (vertices1, bottom);
    GtsSurface * s2 = triangulate (vertices2, e1);

    if (s1 == NULL || s2 == NULL) {
      GSList * i = shape;

      while (i) {
	fprintf (stderr, "%g %g\n", 
		 GTS_POINT (i->data)->x, 
		 GTS_POINT (i->data)->y); 
	i = i->next;
      }
      fprintf (stderr, "\n");
      if (s1) gts_object_destroy (GTS_OBJECT (s1));
      if (s2) gts_object_destroy (GTS_OBJECT (s2));
    }
    else {
      gts_surface_foreach_face (s2, (GtsFunc) gts_triangle_revert, NULL);
      
      gts_surface_merge (s, s1);
      gts_surface_merge (s, s2);
    
      gts_object_destroy (GTS_OBJECT (s1));
      gts_object_destroy (GTS_OBJECT (s2));
    }
    g_slist_free (vertices1);
    g_slist_free (vertices2);
  }

  g_slist_free (bottom);
  g_slist_free (e1);
}

static void surface_add_ellipse_shape (GtsSurface * s,
				       gdouble x, gdouble y,
				       gdouble radius,
				       gdouble theta,
				       gdouble thetamax,
				       gdouble e,
				       gdouble z1, gdouble z2,
				       guint np,
				       gboolean closed)
{
  GSList * shape = NULL;
  guint i, npm = np;

  g_return_if_fail (s != NULL);
  g_return_if_fail (np >= 3);

  if (thetamax < 2.*M_PI)
    npm = np + 1;
  for (i = 0; i < npm; i++) {
    gdouble theta1 = theta + i*thetamax/(gdouble) np;
    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (),
					    x + radius*e*cos (theta1),
					    y + radius/e*sin (theta1),
					    0.));
  }
  surface_add_shape (s, shape, z1, z2, 10, TRUE, closed);
  g_slist_free (shape);
}

static void surface_add_star_shape (GtsSurface * s,
				    gdouble dr,
				    gdouble z1, gdouble z2,
				    guint np,
				    gboolean closed)
{
  GSList * shape = NULL;
  guint i;

  g_return_if_fail (s != NULL);
  g_return_if_fail (np >= 3);

  for (i = 0; i < np; i++) {
    gdouble theta = .001 + 2.*i*M_PI/np;
    gdouble radius = 0.45 - dr + dr*cos (6.*theta);

    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (),
					    radius*cos (theta),
					    radius*sin (theta),
					    0.));
  }
  surface_add_shape (s, shape, z1, z2, 10, TRUE, closed);
  g_slist_free (shape);
}

static gdouble shape_func_bottom (gdouble x)
{
  gdouble y1 = 0.2/4.;
  gdouble y2 = 1e-6/4.;

  if (x <= -0.25)
    return y1;
  if (x < 0.25)
    return y2 + 0.5*(y1 - y2)*(1. + cos (2.*M_PI*(x + 0.25)));
  return y2;
}

static void surface_add_channel_shape (GtsSurface * s,
				       gdouble z1, gdouble z2,
				       guint np,
				       gboolean closed)
{
  GSList * shape = NULL;
  gint i;

  g_return_if_fail (s != NULL);

  for (i = np - 1; i >= 0; i--) {
    gdouble x = - 0.501 + 1.002*(gdouble) i/(gdouble) (np - 1);
    gdouble y = shape_func_bottom (x) - 0.125;

    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (), x, y, 0.));
  }
  for (i = 0; i < np; i++) {
    gdouble x = - 0.501 + 1.002*(gdouble) i/(gdouble) (np - 1);
    gdouble y = 0.25 - shape_func_bottom (x) - 0.125;

    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (), x, y, 0.));
  }
  shape = g_slist_reverse (shape);
  surface_add_shape (s, shape, z1, z2, 10, TRUE, closed);
  g_slist_free (shape);
}

static void surface_add_rectangular_channel_shape (GtsSurface * s,
						   gdouble z1, gdouble z2,
						   gboolean closed)
{
  GSList * shape = NULL;

  g_return_if_fail (s != NULL);

  shape = g_slist_prepend (shape, 
		  gts_point_new (gts_point_class (), -0.5001, -0.125001, 0.));
  shape = g_slist_prepend (shape, 
		  gts_point_new (gts_point_class (), -0.5001,  0.125001, 0.));
  shape = g_slist_prepend (shape, 
		  gts_point_new (gts_point_class (),  0.5001, 0.125001, 0.));
  shape = g_slist_prepend (shape, 
		  gts_point_new (gts_point_class (),  0.5001, -0.125001, 0.));
  surface_add_shape (s, shape, z1, z2, 10, TRUE, closed);
  g_slist_free (shape);
}

static void surface_add_witch_shape (GtsSurface * s,
				     gdouble xo, gdouble h, gdouble lh,
				     gdouble z1, gdouble z2,
				     guint np)
{
  GSList * shape = NULL;
  guint i;

  g_return_if_fail (s != NULL);
  g_return_if_fail (np >= 2);
  g_return_if_fail (lh > 0.);

  shape = g_slist_prepend (shape, 
			   gts_point_new (gts_point_class (),
					  -0.5001, -0.5002, 0.));
  for (i = 0; i < np; i++) {
    gdouble x = -0.5001 + 1.0002*i/(gdouble) (np - 1);

    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (),
					    x,
				 h/(1. + (x - xo)*(x - xo)/(lh*lh)) - 0.5001,
					    0.));
  }
  shape = g_slist_prepend (shape, 
			   gts_point_new (gts_point_class (),
					  0.5001, -0.5002, 0.));
  shape = g_slist_reverse (shape);
  surface_add_shape (s, shape, z1, z2, 10, TRUE, TRUE);
  g_slist_free (shape);
}

static void surface_add_rayleigh_taylor_shape (GtsSurface * s,
					       gdouble yo, gdouble a,
					       gdouble z1, gdouble z2,
					       guint np)
{
  GSList * shape = NULL;
  guint i;

  g_return_if_fail (s != NULL);
  g_return_if_fail (np >= 2);

  shape = g_slist_prepend (shape, 
			   gts_point_new (gts_point_class (),
					  -0.5001, -10., 0.));
  for (i = 0; i < np; i++) {
    gdouble x = -0.5001 + 1.0002*i/(gdouble) (np - 1);

    shape = g_slist_prepend (shape, 
			     gts_point_new (gts_point_class (),
					    x,
			     yo + a*cos (i*2.*M_PI/(np - 1)),
					    0.));
  }
  shape = g_slist_prepend (shape, 
			   gts_point_new (gts_point_class (),
					  0.5001, -10., 0.));
  shape = g_slist_reverse (shape);
  surface_add_shape (s, shape, z1, z2, 10, TRUE, TRUE);
  g_slist_free (shape);
}

static gboolean surface_add_file_shape (GtsSurface * s,
					gdouble z1, gdouble z2,
					gboolean zextrude,
					FILE * fp)
{
  gdouble xs = G_MAXDOUBLE, ys = G_MAXDOUBLE, x, y;
  GSList * shape = NULL;
  guint length;
  GtsSurface * s1, * self = NULL;

  if (zextrude && fscanf (fp, "%lf %lf", &z1, &z2) != 2)
    return FALSE;

  while (fscanf (fp, "%lf %lf", &x, &y) == 2) {
    if (x != xs || y != ys)
      shape = g_slist_prepend (shape, 
			       gts_point_new (gts_point_class (), x, y, 0));
    if (xs == G_MAXDOUBLE) { xs = x; ys = y; }
  }
  fgetc (fp);
  length = g_slist_length (shape);
  if (length > 0 && length < 3) {
    g_slist_free (shape);
    return FALSE;
  }
  s1 = gts_surface_new (gts_surface_class (), 
			s->face_class, s->edge_class, s->vertex_class);
  surface_add_shape (s1, shape, z1, z2, zextrude ? 1 : 10, TRUE, TRUE);
  if (!gts_surface_is_orientable (s1) || 
      !gts_surface_is_closed (s1) ||
      (self = gts_surface_is_self_intersecting (s1))) {
    GSList * i = shape;
    
    while (i) {
      fprintf (stderr, "%g %g\n", 
	       GTS_POINT (i->data)->x, 
	       GTS_POINT (i->data)->y); 
      i = i->next;
    }
    fprintf (stderr, "\n");
    if (self)
      gts_object_destroy (GTS_OBJECT (self));
  }
  else
    gts_surface_merge (s, s1);
  gts_object_destroy (GTS_OBJECT (s1));
  g_slist_free (shape);
  return TRUE;
}

int main (int argc, char * argv[])
{
  GtsSurface * s;
  guint number = 100;
  int c = 0;
  gchar * shape = NULL;
  gboolean verbose = FALSE;
  gdouble dr = 0.15;
  gdouble ratio = 1.;
  gboolean closed = TRUE;
  gboolean zextrude = FALSE;
  FILE * sfp = NULL;

  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"zextrude", no_argument, NULL, 'z'},
      {"number", required_argument, NULL, 'n'},
      {"open", no_argument, NULL, 'o'},
      {"dr", required_argument, NULL, 'd'},
      {"ratio", required_argument, NULL, 'r'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hvn:d:r:oz",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hvn:d:r:oz"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'z': /* zextrude */
      zextrude = TRUE;
      break;
    case 'o': /* open */
      closed = FALSE;
      break;
    case 'n': /* number */
      number = atoi (optarg);
      break;
    case 'd': /* dr */
      dr = atof (optarg);
      break;
    case 'r' :/* ratio */
      ratio = atof (optarg);
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: shapes [OPTIONS] SHAPE\n"
     "Generates various shapes. SHAPE can be one of:\n"
     "  ellipse, star, 4ellipses, square, almgren, channel, half-cylinder,\n"
     "  rayleigh-taylor, FILE\n"
     "\n"
     "  -z    --zextrude    shape file contains z coordinate\n"
     "  -n N  --number=N    set number of points for polar surfaces (default is 100)\n"
     "  -o    --open        generate open surfaces\n"
     "  -d R  --dr=R        set inner radius for star to R (default is 0.15)\n"
     "  -r R  --ratio=R     ratio x/y of the ellipse (default is 1)\n"
     "  -v    --verbose     display surface statistics\n"
     "  -h    --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       GTS_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `shapes --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing SHAPE */  
    fprintf (stderr, 
	     "shapes: missing SHAPE\n"
	     "Try `shapes --help' for more information.\n");
    return 1; /* failure */
  }
  shape = argv[optind++];

  s = gts_surface_new (gts_surface_class (),
		       gts_face_class (),
		       gts_edge_class (),
		       gts_vertex_class ());
  if (!strcmp (shape, "ellipse"))
    surface_add_ellipse_shape (s, 0., 0., 0.25, 0.001, 2.*M_PI,
			       sqrt (1./ratio), - 1., 1., number, closed);
  else if (!strcmp (shape, "star"))
    surface_add_star_shape (s, dr, - 1., 1., number, closed);
  else if (!strcmp (shape, "4ellipses")) {
    surface_add_ellipse_shape (s, 0.25, 0.25, 5./32./sqrt (2.), 0.001, 2.*M_PI,
			       sqrt (2.), - 1., 1., number, closed);
    surface_add_ellipse_shape (s, -0.25, 0.25, 5./32./sqrt (2.), 0.001, 2.*M_PI,
			       sqrt (2.), - 1., 1., number, closed);
    surface_add_ellipse_shape (s, 0.25, -0.25, 5./32./sqrt (2.), 0.001, 2.*M_PI,
			       sqrt (2.), - 1., 1., number, closed);
    surface_add_ellipse_shape (s, -0.25, -0.25, 5./32./sqrt (2.), 0.001, 2.*M_PI,
			       sqrt (2.), - 1., 1., number, closed);
  }
  else if (!strcmp (shape, "square"))
    surface_add_ellipse_shape (s, 0., 0., 0.25*sqrt(2.), M_PI/4.,  2.*M_PI, 1.,
			       - 1., 1., 4, closed);
  else if (!strcmp (shape, "almgren")) {
    surface_add_ellipse_shape (s, 0.25, 0.25, 0.1, 0.001,  2.*M_PI, 1., -1., 1., 
			       number, closed);
    surface_add_ellipse_shape (s, -0.25, 0.125, sqrt (0.15*0.1),
			       0.001, 2.*M_PI,
			       0.15/sqrt (0.15*0.1),
			       -1., 1., number, closed);
    surface_add_ellipse_shape (s, 0., -0.25, sqrt (0.2*0.1),
			       0.001, 2.*M_PI,
			       0.2/sqrt (0.2*0.1),
			       -1., 1., number, closed);
  }
  else if (!strcmp (shape, "channel"))
    surface_add_channel_shape (s, -1., 1., number, closed);
  else if (!strcmp (shape, "half-cylinder")) {
    //    surface_add_rectangular_channel_shape (s, -1., 1., closed);
    surface_add_ellipse_shape (s, -0.375001, 0., 0.03125001, 
			       M_PI/2., M_PI, 1., -1., 1., 
			       number, closed);
  }
  else if (!strcmp (shape, "witch"))
    surface_add_witch_shape (s, -0.25, 0.05, 0.05, -1., 1., number);
  else if (!strcmp (shape, "rayleigh-taylor"))
    surface_add_rayleigh_taylor_shape (s, 0., 0.025, -1., 1., number);
  else if (!strcmp (shape, "annulus")) {
    surface_add_ellipse_shape (s, 0., 0., 0.5, 0.001, 2.*M_PI,
			       1., - 2., 2., number, closed);
    gts_surface_foreach_face (s, (GtsFunc) gts_triangle_revert, NULL);
    surface_add_ellipse_shape (s, 0., 0., 0.25, 0.001, 2.*M_PI,
			       1., - 1., 1., number, closed);
  }
  else if ((!strcmp (shape, "-") && (sfp = stdin)) || 
	   (sfp = fopen (shape, "rt")) != NULL) {
    while (!feof (sfp))
      if (!surface_add_file_shape (s, -1., 1., zextrude, sfp)) {
	fprintf (stderr, 
		 "shapes: file `%s' is not a valid shape\n"
		 "Try `shapes --help' for more information.\n",
		 shape);
	return 1; /* failure */
      }
    fclose (sfp);
  }
  else {
    fprintf (stderr, 
	     "shapes: unknown shape `%s'\n"
	     "Try `shapes --help' for more information.\n",
	     shape);
    return 1; /* failure */
  }
  if (verbose)
    gts_surface_print_stats (s, stderr);
  gts_surface_write (s, stdout);

  return 0;
}
				      
