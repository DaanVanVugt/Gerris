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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include "kdt/kdt.h"
#include "rsurface.h"

static int verbose = 0;

static int query (double p[3], void * data)
{
  static long n = 0;
  KdtHeap * h = data;
  KdtPoint q = { p[0], p[1], p[2] };
  kdt_heap_put (h, &q);
  if (verbose && n % 3571 == 0)
    fprintf (stderr, "rsurface2kdt: reading points... %ld\r", n);
  n++;
  return 0;
}

static int includes_true (KdtRect rect)
{
  return 1;
}

static void progress (float complete, void * data)
{
#if 0 /* this doesn't work yet */
  struct timeval * start = data, now;
  gettimeofday (&now, NULL);
  double remaining = ((double) (now.tv_usec - start->tv_usec) + 
		      1e6*(double) (now.tv_sec - start->tv_sec))*(1. - complete);
  int hours = remaining/1e6/3600.;
  int mins = remaining/1e6/60. - 60.*hours;
  int secs = remaining/1e6 - 3600.*hours - 60.*mins;
  fprintf (stderr, "\rrsurface2kdt: %3.0f%% complete %02d:%02d:%02d remaining", 
	   (complete > 1. ? 1. : complete)*100.,
	   hours, mins, secs);
#else
  fprintf (stderr, "\rrsurface2kdt: %3.0f%% complete    ", 
	   (complete > 1. ? 1. : complete)*100.);
#endif
}

int main (int argc, char * argv[])
{
  int c = 0, pagesize = 4096;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"pagesize", required_argument, NULL, 'p'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "p:hv",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "p:hv"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'v': /* verbose */
      verbose = 1;
      break;
    case 'p': /* pagesize */
      pagesize = atoi (optarg);
      break;
    case 'h': /* help */
      fprintf (stderr,
	       "Usage: rsurface2kdt [OPTION] BASENAME\n"
	       "\n"
	       "Converts an existing R-tree database to the KDT format.\n"
	       "\n"
	       "  -p N  --pagesize=N  sets the pagesize in bytes (default is 4096)\n"
	       "  -v    --verbose     display progress bar\n"
	       "  -h    --help        display this help and exit\n"
	       "\n"
	       "Report bugs to %s\n",
	       "popinet@users.sf.net");
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `rsurface2kdt -h' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing BASENAME */
    fprintf (stderr, 
	     "rsurface2kdt: missing BASENAME\n"
	     "Try `rsurface2kdt -h' for more information.\n");
    return 1; /* failure */
  }

  RSurface * rs = r_surface_open (argv[optind], "r", 0);
  if (rs == NULL) {
    fprintf (stderr, "rsurface2kdt: could not open `%s'\n", argv[optind]);
    return 1;
  }

  if (verbose)
    fprintf (stderr, "rsurface2kdt: reading points...\r");
  KdtHeap h;
  kdt_heap_create (&h, kdt_tmpfile (), 0, -1, 1000000);
  double min[2] = { -1e100, -1e100 }, max[2] = { 1e100, 1e100 };
  r_surface_query_region (rs, min, max, query, &h);
  kdt_heap_flush (&h);
  r_surface_close (rs);

  if (verbose)
    fprintf (stderr, "rsurface2kdt:   0%% complete              ");

  struct timeval start;
  gettimeofday (&start, NULL);
  Kdt * kdt = kdt_new ();
  kdt_create (kdt, argv[optind], pagesize, &h, verbose ? progress : NULL, &start);
  kdt_destroy (kdt);

  if (verbose) {
    Kdt * kdt = kdt_new ();
    KdtRect rect = {{-1e30,1e30},{-1e30,1e30}};
    KdtSum sum;

    kdt_sum_init (&sum);
    assert (!kdt_open (kdt, argv[optind]));
    long n = kdt_query_sum (kdt, 
			    (KdtCheck) includes_true, (KdtCheck) includes_true, NULL,
			    rect, &sum);
    fprintf (stderr,
	     "\r%ld points Height min: %g average: %g max: %g\n",
	     n, sum.Hmin, sum.H0/sum.w, sum.Hmax);
    kdt_destroy (kdt);
  }

  char * name = malloc (strlen (argv[optind]) + strlen (".DataPD") + 1);
  sprintf (name, "%s.DataPD", argv[optind]);
  remove (name);
  sprintf (name, "%s.Data", argv[optind]);
  remove (name);
  sprintf (name, "%s.DirPD", argv[optind]);
  remove (name);
  remove (argv[optind]);
  free (name);

  return 0;
}
