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

#include "init.h"
#include "simulation.h"

typedef struct {
  GfsSimulation * s1, * s2;
  GfsVariable * v1, * v2;
} CombineData;

static void relativepos (FttCell * cell, GfsSimulation * s1, GfsSimulation * s2, FttVector * pos)
{
  ftt_cell_pos (cell, pos);
  gfs_simulation_map_inverse (s1, pos);
  gfs_simulation_map (s2, pos);
}

static void refine (FttCell * cell, CombineData * p)
{
  FttVector pos;
  relativepos (cell, p->s1, p->s2, &pos);
  FttCell * matching = gfs_domain_locate (GFS_DOMAIN (p->s2), pos, -1, NULL);
  if (matching && ftt_cell_level (cell) < ftt_cell_level (matching))
    ftt_cell_refine_single (cell, GFS_DOMAIN (p->s1)->cell_init, 
			    GFS_DOMAIN (p->s1)->cell_init_data);
}

static void combine (FttCell * cell, CombineData * p)
{
  FttVector pos;
  relativepos (cell, p->s1, p->s2, &pos);
  FttCell * matching = gfs_domain_locate (GFS_DOMAIN (p->s2), pos, -1, NULL);
  if (matching) {
    g_assert (ftt_cell_level (cell) >= ftt_cell_level (matching));
    GFS_VALUE (cell, p->v1) = MAX (GFS_VALUE (cell, p->v1), GFS_VALUE (matching, p->v2));
  }
}

int main (int argc, char * argv[])
{
  GtsFile * fp;
  FILE * f;
  int c = 0;
  gchar * name;
  GfsVariable * var1, * var2;
  GfsSimulation * s1, * s2;
  
  gboolean verbose = FALSE;
  gchar * fname1, * fname2;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hv",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hv"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: gfscombine [OPTION] FILE1 FILE2 VAR\n"
     "Computes the maximum of VAR between the solutions in FILE1 and FILE2\n"
     "and outputs the corresponding simulation.\n"
     "\n"
     "  -v    --verbose     display statistics and other info\n"
     "  -h    --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfscombine --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing FILE1 */  
    fprintf (stderr, 
	     "gfscombine: missing FILE1\n"
	     "Try `gfscombine --help' for more information.\n");
    return 1; /* failure */
  }
  fname1 = argv[optind++];

  if (optind >= argc) { /* missing FILE2 */  
    fprintf (stderr, 
	     "gfscombine: missing FILE2\n"
	     "Try `gfscombine --help' for more information.\n");
    return 1; /* failure */
  }
  fname2 = argv[optind++];

  if (optind >= argc) { /* missing VAR */  
    fprintf (stderr, 
	     "gfscombine: missing VAR\n"
	     "Try `gfscombine --help' for more information.\n");
    return 1; /* failure */
  }
  name = argv[optind++];

  f = fopen (fname1, "rt");
  if (f == NULL) {
    fprintf (stderr, "gfscombine: cannot open file `%s'\n", fname1);
    return 1;
  }
  fp = gts_file_new (f);
  if (!(s1 = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gfscombine: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     fname1, fname1, fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);
  fclose (f);
  gfs_simulation_init (s1);

  f = fopen (fname2, "rt");
  if (f == NULL) {
    fprintf (stderr, "gfscombine: cannot open file `%s'\n", fname2);
    return 1;
  }
  fp = gts_file_new (f);
  if (!(s2 = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gfscombine: file `%s' is not a valid simulation file\n"
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
	     "gfscombine: unknown variable `%s' for `%s'\n"
	     "Try `gfscombine --help' for more information.\n",
	     name, fname1);
    return 1; /* failure */
  }

  var2 = gfs_variable_from_name (GFS_DOMAIN (s2)->variables, name);
  if (var2 == NULL) {
    fprintf (stderr, 
	     "gfscombine: unknown variable `%s' for `%s'\n"
	     "Try `gfscombine --help' for more information.\n",
	     name, fname2);
    return 1; /* failure */
  }

  if (verbose) {
    GfsNorm norm = gfs_domain_norm_variable (GFS_DOMAIN (s1),
					     var1, NULL, FTT_TRAVERSE_LEAFS, -1,
					     NULL, NULL);
    GtsRange s = gfs_domain_stats_variable (GFS_DOMAIN (s1),
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

  CombineData p = { s1, s2, var1, var2 };
  gfs_domain_traverse_leaves (GFS_DOMAIN (s1), (FttCellTraverseFunc) refine, &p);
  gfs_domain_traverse_leaves (GFS_DOMAIN (s1), (FttCellTraverseFunc) combine, &p);

  if (verbose) {
    GfsNorm norm = gfs_domain_norm_variable (GFS_DOMAIN (s1),
					     var1, NULL, FTT_TRAVERSE_LEAFS, -1,
					     NULL, NULL);
    GtsRange s = gfs_domain_stats_variable (GFS_DOMAIN (s1),
					    var1, FTT_TRAVERSE_LEAFS, -1,
					    NULL, NULL);
    gdouble f = pow (s1->physical_params.L, var1->units);
    fprintf (stderr, 
	     "max(%s,%s):\n"
	     "  first: %g second: %g infty: %g w: %g\n"
	     "  min: %g avg: %g | %g max: %g\n",
	     fname1, fname2,
	     norm.first*f, norm.second*f, norm.infty*f, norm.w,
	     s.min*f, s.mean*f, s.stddev*f, s.max*f);
  }

  gfs_simulation_write (s1, -1, stdout);

  gts_object_destroy (GTS_OBJECT (s1));
  gts_object_destroy (GTS_OBJECT (s2));

  return 0;
}
