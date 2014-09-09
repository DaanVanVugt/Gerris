#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "graphic.h"
#include "solid.h"
#include "init.h"
#include "simulation.h"

#if FTT_2D

/* GfsVertex: Header */

typedef struct _GfsVertex         GfsVertex;

struct _GfsVertex {
  /*< private >*/
  GtsVertex parent;

  /*< public >*/
  FttCell * cell;
};

#define GFS_VERTEX(obj)            GTS_OBJECT_CAST (obj,\
					         GfsVertex,\
					         gfs_vertex_class ())
#define IS_GFS_VERTEX(obj)         (gts_object_is_from_class (obj,\
						 gfs_vertex_class ()))

static GtsVertexClass * gfs_vertex_class  (void);

/* GfsVertex: Object */

GtsVertexClass * gfs_vertex_class (void)
{
  static GtsVertexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_vertex_info = {
      "GfsVertex",
      sizeof (GfsVertex),
      sizeof (GtsVertexClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_vertex_class ()),
				  &gfs_vertex_info);
  }

  return klass;
}

static void add_vertex (GSList * merged, GtsSurface * s)
{
  FttVector cm = {0., 0., 0.};
  gdouble ta = 0.;
  GtsVertex * v;
  GSList * i = merged;
  
  while (i) {
    FttVector p;
    FttCell * cell = i->data;
    gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;

    gfs_cell_cm (i->data, &p);
    cm.x += a*p.x; cm.y += a*p.y; cm.z += a*p.z;
    ta += a;
    i = i->next;
  }

  v = gts_vertex_new (s->vertex_class, cm.x/ta, cm.y/ta, cm.z/ta);
  g_assert (gts_delaunay_add_vertex (s, v, NULL) == NULL);
  GFS_VERTEX (v)->cell = merged->data;
}

static GtsSurface * surface_from_domain (GfsDomain * domain)
{
  GtsSurface * s = gts_surface_new (gts_surface_class (),
				    gts_face_class (),
				    gts_edge_class (),
				    GTS_VERTEX_CLASS (gfs_vertex_class ()));
  GtsVertex * v1 = gts_vertex_new (s->vertex_class, -100., -100., 0.);
  GtsVertex * v2 = gts_vertex_new (s->vertex_class,  100., -100., 0.);
  GtsVertex * v3 = gts_vertex_new (s->vertex_class,    0.,  100., 0.);
  GtsEdge * e1 = gts_edge_new (s->edge_class, v1, v2);
  GtsEdge * e2 = gts_edge_new (s->edge_class, v2, v3);
  GtsEdge * e3 = gts_edge_new (s->edge_class, v3, v1);
  
  gts_surface_add_face (s, gts_face_new (s->face_class, e1, e2, e3));
  gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) add_vertex, s);
  gts_allow_floating_vertices = TRUE;
  gts_object_destroy (GTS_OBJECT (v1));
  gts_object_destroy (GTS_OBJECT (v2));
  gts_object_destroy (GTS_OBJECT (v3));
  gts_allow_floating_vertices = FALSE;

  return s;
}

static void difference_triangulated (GfsVertex * v, gpointer * data)
{
  GtsSurface * s = data[0];
  GfsVariable * var1 = data[1];
  GfsVariable * var2 = data[2];
  GfsVariable * e = data[3];
  GtsFace * f = gts_point_locate (GTS_POINT (v), s, NULL);

  if (f != NULL && gts_triangle_quality (GTS_TRIANGLE (f)) > 0.8) {
    GtsVertex * v1, * v2, * v3;
    gdouble x, x1, x2, y, y1, y2, a, b, det;
    gdouble fv3, fv1, fv2;

    gts_triangle_vertices (GTS_TRIANGLE (f), &v1, &v2, &v3);
    x = GTS_POINT (v)->x - GTS_POINT (v1)->x;
    y = GTS_POINT (v)->y - GTS_POINT (v1)->y;
    x1 = GTS_POINT (v2)->x - GTS_POINT (v1)->x;
    y1 = GTS_POINT (v2)->y - GTS_POINT (v1)->y;
    x2 = GTS_POINT (v3)->x - GTS_POINT (v1)->x;
    y2 = GTS_POINT (v3)->y - GTS_POINT (v1)->y;
    det = x1*y2 - x2*y1;
    g_assert (det != 0.);
    a = (x*y2 - y*x2)/det;
    b = (y*x1 - x*y1)/det;
    fv1 = GFS_VALUE (GFS_VERTEX (v1)->cell, var2);
    fv2 = GFS_VALUE (GFS_VERTEX (v2)->cell, var2);
    fv3 = GFS_VALUE (GFS_VERTEX (v3)->cell, var2);
    GTS_POINT (v)->z = GFS_VALUE (v->cell, e) = 
      GFS_VALUE (v->cell, var1) -
      (fv1 + a*(fv2 - fv1) + b*(fv3 - fv1));
  }
}

#endif /* FTT_2D */

static gboolean is_mixed (FttCell * cell, guint level)
{
  if (GFS_IS_MIXED (cell))
    return TRUE;
  if (!FTT_CELL_IS_ROOT (cell) && ftt_cell_level (cell) > level)
    return is_mixed (ftt_cell_parent (cell), level);
  return FALSE;
}

static void inject (FttCell * cell, GfsVariable * e)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	GFS_VALUE (child.c[i], e) = GFS_VALUE (cell, e);
	inject (child.c[i], e);
      }
  }
}

static gboolean difference_tree (FttCell * cell,
				 GfsDomain * ref,
				 GfsVariable * v1,
				 GfsVariable * v2,
				 GfsVariable * e,
				 gdouble period)
{
  guint level = ftt_cell_level (cell);
  FttVector pos;
  FttCell * locate;
  gboolean added = FALSE;
  
  ftt_cell_pos (cell, &pos);
  pos.x += period;
  locate = gfs_domain_locate (ref, pos, level, NULL);
  if (locate == NULL) {
    pos.x -= 2.*period;
    locate = gfs_domain_locate (ref, pos, level, NULL);
  }
  if (locate == NULL) {
    fprintf (stderr, "gfscompare: the files are not comparable\n");
    exit (1);
  }
  if (ftt_cell_level (locate) != level)
    return FALSE;
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && difference_tree (child.c[i], ref, v1, v2, e, period))
	added = TRUE;
  }
  if (!added) {
    if (GFS_HAS_DATA (cell, v1) && GFS_HAS_DATA (locate, v2))
      GFS_VALUE (cell, e) = (GFS_VALUE (cell, v1) - GFS_VALUE (locate, v2));
    else if (!GFS_HAS_DATA (cell, v1) && GFS_HAS_DATA (locate, v2))
      GFS_VALUE (cell, e) = GFS_VALUE (locate, v2);
    else if (GFS_HAS_DATA (cell, v1) && !GFS_HAS_DATA (locate, v2))
      GFS_VALUE (cell, e) = GFS_VALUE (cell, v1);
    else
      GFS_VALUE (cell, e) = GFS_NODATA;
    inject (cell, e);
  }
  return TRUE;
}

static void difference_box (GfsBox * box, gpointer * data)
{
  gdouble * period = data[4];

  difference_tree (box->root, data[0], data[1], data[2], data[3], *period);
}

static void difference_constant (FttCell * cell, gpointer * data)
{
  gint full = *((gint *) data[0]);
  gdouble * sum = data[1];
  gboolean * centered = data[3];
  gboolean * weighted = data[4];
  gdouble * weight = data[5];
  GfsVariable * e = data[6];
  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;

  if ((full == -2 || 
       (full == -1 && !GFS_IS_MIXED (cell)) ||
       (full >= 0 && !is_mixed (cell, full))) &&
      (!(*centered) || a >= 0.5)) {
    gdouble w = *weighted ? ftt_cell_volume (cell)*a : 1.;

    *sum += w*GFS_VALUE (cell, e);
    *weight += w;
  }
}

static void difference (FttCell * cell, gpointer * data)
{
  gint full = *((gint *) data[0]);
  GfsNorm * norm = data[1];
  gboolean * histogram = data[2];
  gboolean * centered = data[3];
  gboolean * weighted = data[4];
  gdouble * constant = data[5];
  gboolean * mixed = data[6];
  GfsVariable * e = data[7];
  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;

  if ((!(*mixed) || a < 1.) &&
      (full == -2 || 
       (full == -1 && !GFS_IS_MIXED (cell)) ||
       (full >= 0 && !is_mixed (cell, full))) &&
      (!(*centered) || a >= 0.5)) {
    gfs_norm_add (norm, GFS_VALUE (cell, e) - *constant,
		  *weighted ? ftt_cell_volume (cell)*a : 1.);
    if (*histogram)
      printf ("%g %g\n", GFS_VALUE (cell, e), a);
  }
  else
    GFS_VALUE (cell, e) = 0.;
}

static void compute_gradient (FttCell * cell, gpointer * data) 
{
  GfsVariable * v = data[0];
  FttComponent * c = data[1];
  GfsVariable * g = data[2];

  GFS_VALUE (cell, g) = 
    gfs_center_gradient (cell, *c, v->i)/ftt_cell_size (cell);
}

static void compute_log (FttCell * cell, GfsVariable * e) 
{
  GFS_VALUE (cell, e) = log10 (fabs (GFS_VALUE (cell, e)) + 1e-10);
}

static void compute_absolute (FttCell * cell, GfsVariable * e)
{
  GFS_VALUE (cell, e) =  fabs (GFS_VALUE (cell, e));
}

static void difference_centered (FttCell * cell, gpointer * data)
{
  GfsDomain * ref = data[0];
  GfsVariable * v1 = data[1];
  GfsVariable * v2 = data[2];
  GfsVariable * e = data[3];
  FttVector p;
  FttCell * locate;

  gfs_cell_cm (cell, &p);
  locate = gfs_domain_locate (ref, p, -1, NULL);
  if (locate == NULL || ftt_cell_level (locate) < ftt_cell_level (cell)) {
    fprintf (stderr, "gfscompare: the files are not comparable\n");
    exit (1);
  }
  GFS_VALUE (cell, e) = GFS_VALUE (cell, v1) - gfs_interpolate (locate, p, v2);
}

int main (int argc, char * argv[])
{
  GtsFile * fp;
  FILE * f;
  int c = 0;
  gchar * name;
  GfsVariable * var1, * var2, * e;
  GfsSimulation * s1, * s2;
  
  gboolean verbose = FALSE;
  gint full = -2;
  gboolean no_check = FALSE;
  gboolean output = FALSE;
  gboolean squares = FALSE;
  gboolean take_log = FALSE;
  gchar * fname1, * fname2;
  gdouble period = 0.;

  FttComponent gradient = FTT_DIMENSION;

  GfsNorm norm;
  gpointer data[8];

  gboolean refined_error = FALSE;
  gboolean histogram = FALSE;
  gboolean centered = FALSE;
  gboolean weighted = TRUE;
  gdouble constant = 0.;
  gboolean absolute = FALSE;
#if FTT_2D
  gboolean gnuplot = FALSE;
  gboolean triangulate = FALSE;
#endif /* FTT_2D */
  gdouble min = G_MAXDOUBLE, max = - G_MAXDOUBLE;
  gboolean mixed = FALSE;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
#if FTT_2D
      {"gnuplot", no_argument, NULL, 'G'},
      {"triangulate", no_argument, NULL, 't'},
#endif /* FTT_2D */
      {"mixed", no_argument, NULL, 'x'},
      {"min", required_argument, NULL, 'm'},
      {"max", required_argument, NULL, 'M'},
      {"period", required_argument, NULL, 'p'},
      {"histogram", no_argument, NULL, 'H'},
      {"refined", no_argument, NULL, 'r'},
      {"log", no_argument, NULL, 'l'},
      {"abs", no_argument, NULL, 'a'},
      {"full", required_argument, NULL, 'f'},
      {"gradient", required_argument, NULL, 'g'},
      {"output", no_argument, NULL, 'o'},
      {"squares", no_argument, NULL, 'S'},
      {"nocheck", no_argument, NULL, 'n'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      {"centered", no_argument, NULL, 'c'},
      {"not-weighted", no_argument, NULL, 'w'},
      {"constant", no_argument, NULL, 'C'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hvnog:f:lSrHp:cwCeaGm:M:xt",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hvnog:f:lSrHp:cwCeaGm:M:xt"))) {
#endif /* not HAVE_GETOPT_LONG */
#if FTT_2D
    case 'G': /* gnuplot */
      gnuplot = TRUE;
      break;
    case 't': /* triangulate */
      triangulate = TRUE;
      break;
#endif /* FTT_2D */
    case 'x': /* mixed */
      mixed = TRUE;
      break;
    case 'm': /* min */
      min = atof (optarg);
      break;
    case 'M': /* max */
      max = atof (optarg);
      break;
    case 'a': /* abs */
      absolute = TRUE;
      break;
    case 'C': /* constant */
      constant = TRUE;
      break;
    case 'w': /* not-weighted */
      weighted = FALSE;
      break;
    case 'c': /* centered */
      centered = TRUE;
      break;
    case 'p': /* period */
      period = atof (optarg);
      break;
    case 'H': /* histogram */
      histogram = TRUE;
      break;
    case 'r': /* refined */
      refined_error = TRUE;
      break;
    case 'l': /* log */
      take_log = TRUE;
      break;
    case 'f': /* full */
      full = atoi (optarg);
      break;
    case 'g': /* gradient */
      gradient = atoi (optarg);
      if (gradient >= FTT_DIMENSION) {
	fprintf (stderr, 
		 "gfscompare: invalid argument for option `gradient'.\n"
		 "Try `gfscompare --help' for more information.\n");
	return 1; /* failure */
      }
      break;
    case 'S': /* squares */
      squares = TRUE;
      break;
    case 'o': /* output */
      output = TRUE;
      break;
    case 'n': /* nocheck */
      no_check = TRUE;
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: gfscompare [OPTION] FILE1 FILE2 VAR\n"
     "Computes the difference between the solutions in FILE1 and FILE2\n"
     "for variable VAR.\n"
     "\n"
     "  -x    --mixed       compute error only in mixed cells\n"
     "  -m V  --min=V       set minimum of color scale to V (used with -S)\n"
     "  -M V  --max=V       set maximum of color scale to V\n"
     "  -a    --abs         output the absolute value of the error field\n"
     "  -C    --constant    apply a constant shift to one of the field, minimizing\n"
     "                      the error between the two fields (useful for pressure)\n"
     "  -w    --not-weighted do not use area-weighted norm estimation\n"
     "  -c    --centered    use error estimation for cell-centered variables\n"
     "  -p P  --period=P    shifts FILE1 by P along the x axis\n"
     "  -H    --histogram   output (error,volume) pairs for each cell used\n"
     "                      to compute the error norms\n"
     "  -o    --output      output a GTS representation of the error field\n"
     "  -S    --squares     output an OOGL representation of the error field\n"
#if FTT_2D
     "  -G    --gnuplot     output a gnuplot representation of the error field\n"
     "  -t    --triangulate use center of mass triangulation\n"
#endif /* FTT_2D */
     "  -l    --log         output the log10 of the absolute value of the error field\n"
     "  -f L  --full=L      compare only leaf cells descendants of a cell full at level L\n"
     "                      or all full leaf cells if L = -1\n"
     "  -r    --refined     display error norm on the finest grid\n"
     "  -n    --nocheck     do not check solid fractions\n"
     "  -g C  --gradient=C  use the C component of the gradient of VAR\n"
     "  -v    --verbose     display difference statistics and other info\n"
     "  -h    --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfscompare --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing FILE1 */  
    fprintf (stderr, 
	     "gfscompare: missing FILE1\n"
	     "Try `gfscompare --help' for more information.\n");
    return 1; /* failure */
  }
  fname1 = argv[optind++];

  if (optind >= argc) { /* missing FILE2 */  
    fprintf (stderr, 
	     "gfscompare: missing FILE2\n"
	     "Try `gfscompare --help' for more information.\n");
    return 1; /* failure */
  }
  fname2 = argv[optind++];

  if (optind >= argc) { /* missing VAR */  
    fprintf (stderr, 
	     "gfscompare: missing VAR\n"
	     "Try `gfscompare --help' for more information.\n");
    return 1; /* failure */
  }
  name = argv[optind++];

  f = fopen (fname1, "rt");
  if (f == NULL) {
    fprintf (stderr, "gfscompare: cannot open file `%s'\n", fname1);
    return 1;
  }
  fp = gts_file_new (f);
  if (!(s1 = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gfscompare: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     fname1, fname1, fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);
  fclose (f);
  gfs_simulation_init (s1);

  f = fopen (fname2, "rt");
  if (f == NULL) {
    fprintf (stderr, "gfscompare: cannot open file `%s'\n", fname2);
    return 1;
  }
  fp = gts_file_new (f);
  if (!(s2 = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gfscompare: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     fname2, fname2, fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);
  fclose (f);
  gfs_simulation_init (s2);

  var1 = gfs_variable_from_name (GFS_DOMAIN (s1)->variables, name);
  if (var1 == NULL) {
    fprintf (stderr, 
	     "gfscompare: unknown variable `%s' for `%s'\n"
	     "Try `gfscompare --help' for more information.\n",
	     name, fname1);
    return 1; /* failure */
  }

  var2 = gfs_variable_from_name (GFS_DOMAIN (s2)->variables, name);
  if (var2 == NULL) {
    fprintf (stderr, 
	     "gfscompare: unknown variable `%s' for `%s'\n"
	     "Try `gfscompare --help' for more information.\n",
	     name, fname2);
    return 1; /* failure */
  }

  if (verbose) {
    GtsRange s;

    norm = gfs_domain_norm_variable (GFS_DOMAIN (s1),
				     var1, NULL, FTT_TRAVERSE_LEAFS, -1,
				     NULL, NULL);
    s = gfs_domain_stats_variable (GFS_DOMAIN (s1),
				   var1, FTT_TRAVERSE_LEAFS, -1,
				   NULL, NULL);
    gdouble f = pow (s1->physical_params.L, var1->units);
    fprintf (stderr, 
	     "%s:\n"
	     "  first: %g second: %g infty: %g w: %g\n"
	     "  min: %g avg: %g | %g max: %g\n",
	     fname1, 
	     norm.first*f, norm.second*f, norm.infty*f, norm.w,
	     s.min*f, s.mean*f, s.stddev*f, s.max*f);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (s2),
				     var2, NULL, FTT_TRAVERSE_LEAFS, -1,
				     NULL, NULL);
    s = gfs_domain_stats_variable (GFS_DOMAIN (s2),
				   var2, FTT_TRAVERSE_LEAFS, -1,
				   NULL, NULL);
    fprintf (stderr, 
	     "%s:\n"
	     "  first: %g second: %g infty: %g w: %g\n"
	     "  min: %g avg: %g | %g max: %g\n",
	     fname2, 
	     norm.first*f, norm.second*f, norm.infty*f, norm.w,
	     s.min*f, s.mean*f, s.stddev*f, s.max*f);
  }

  if (gradient < FTT_DIMENSION) {
    gpointer data[3];
    GfsVariable * g1 = gfs_temporary_variable (GFS_DOMAIN (s1));
    GfsVariable * g2 = gfs_temporary_variable (GFS_DOMAIN (s2));

    data[0] = var1;
    data[1] = &gradient;
    data[2] = g1;
    gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_gradient, data);
    data[0] = var2;
    data[2] = g2;
    gfs_domain_cell_traverse (GFS_DOMAIN (s2), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_gradient, data);
    var1 = g1;
    var2 = g2;
  }

  data[0] = s2;
  data[1] = var1;
  data[2] = var2;
  data[3] = e = gfs_temporary_variable (GFS_DOMAIN (s1));
  if (centered)
    gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) difference_centered, data);
#if FTT_2D
  else if (triangulate) {
    GtsSurface * ss1, * ss2;
    gpointer data[4];

    gfs_clock_start (GFS_DOMAIN (s1)->timer);
    gfs_clock_start (GFS_DOMAIN (s2)->timer);
    gfs_simulation_refine (s1);
    gfs_simulation_refine (s2);
    gfs_set_merged (GFS_DOMAIN (s1));
    gfs_set_merged (GFS_DOMAIN (s2));
    gfs_clock_stop (GFS_DOMAIN (s1)->timer);
    gfs_clock_stop (GFS_DOMAIN (s2)->timer);
    ss1 = surface_from_domain (GFS_DOMAIN (s1));
    ss2 = surface_from_domain (GFS_DOMAIN (s2));
    data[0] = ss2;
    data[1] = var1;
    data[2] = var2;
    data[3] = e;
    gts_surface_foreach_vertex (ss1, (GtsFunc) difference_triangulated, data);
  }
#endif /* FTT_2D */
  else {
    gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_get_from_below_intensive, var1);
    gfs_domain_cell_traverse (GFS_DOMAIN (s2), 
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_get_from_below_intensive, var2);
    data[4] = &period;
    gts_container_foreach (GTS_CONTAINER (s1), (GtsFunc) difference_box, data);
  }

  data[0] = &full;
  data[2] = &histogram;
  data[3] = &centered;
  data[4] = &weighted;
  if (constant) {
    gdouble sum = 0., weight = 0.;

    data[1] = &sum;
    data[5] = &weight;
    data[6] = e;
    gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) difference_constant, data);
    constant = weight > 0. ? sum/weight : 0.;
  }
  
  gfs_norm_init (&norm);
  data[1] = &norm;
  data[5] = &constant;
  data[6] = &mixed;
  data[7] = e;
  gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) difference, data);
  gfs_norm_update (&norm);
  if (verbose) {
    gdouble f = pow (s1->physical_params.L, var1->units);
    fprintf (stderr, 
	     "total err first: %10.3e second: %10.3e infty: %10.3e w: %g\n",
	     norm.first*f, norm.second*f, norm.infty*f, norm.w);
    if (refined_error) {
      norm = gfs_domain_norm_variable (GFS_DOMAIN (s1),
				       e, NULL, FTT_TRAVERSE_LEVEL,
				       gfs_domain_depth (GFS_DOMAIN (s1)),
				       NULL, NULL);
      fprintf (stderr, 
	       "refined err first: %10.3e second: %10.3e infty: %10.3e w: %g\n",
	       norm.first*f, norm.second*f, norm.infty*f, norm.w);
    }
  }

  if (output ||
#if FTT_2D
      gnuplot ||
#endif /* FTT_2D */
      squares) {
    if (take_log)
      gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			       (FttCellTraverseFunc) compute_log, e);
    else if (absolute)
      gfs_domain_cell_traverse (GFS_DOMAIN (s1), 
			       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			       (FttCellTraverseFunc) compute_absolute, e);
    if (squares) {
      GtsRange stats = gfs_domain_stats_variable (GFS_DOMAIN (s1), e, FTT_TRAVERSE_LEAFS, -1,
						  NULL, NULL);

      gfs_write_squares (GFS_DOMAIN (s1), e, 
			 min < G_MAXDOUBLE ? min : stats.min, 
			 max > - G_MAXDOUBLE ? max : stats.max,
			 FTT_TRAVERSE_LEAFS, -1, 
			 NULL, stdout);
    }
#if FTT_2D
    else if (gnuplot)
      gfs_write_gnuplot (GFS_DOMAIN (s1), e,
			 FTT_TRAVERSE_LEAFS, -1, 
			 NULL, stdout);
#endif /* FTT_2D */
    else
	gfs_write_gts (GFS_DOMAIN (s1), e, FTT_TRAVERSE_LEAFS, -1, NULL, stdout);
  }

  return 0;
}
