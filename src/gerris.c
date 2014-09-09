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

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "init.h"
#include "simulation.h"
#include "refine.h"
#include "output.h"
#include "adaptive.h"
#include "solid.h"
#include "version.h"

static void set_box_pid (GfsBox * box, gint * pid)
{
  box->pid = *pid;
}

static void setup_binary_IO (GfsDomain * domain)
{
  /* make sure that all the variables are sent */
  g_slist_free (domain->variables_io);
  domain->variables_io = NULL;
  GSList * i = domain->variables;
  while (i) {
    if (GFS_VARIABLE (i->data)->name)
      domain->variables_io = g_slist_append (domain->variables_io, i->data);
    i = i->next;
  }
  domain->binary = TRUE;	
}

static gboolean set_macros ()
{
#ifndef HAVE_M4
  gfs_error (0, "gerris: macros are not supported on this system\n");
  return 1;
#endif /* not HAVE_M4 */
  return 0;
}

int main (int argc, char * argv[])
{
  GfsSimulation * simulation;
  GfsDomain * domain;
  FILE * fptr;
  GtsFile * fp;
  int c = 0;
  guint split = 0;
  guint npart = 0;
  gboolean profile = FALSE, macros = FALSE, one_box_per_pe = TRUE, bubble = FALSE, verbose = FALSE;
  gchar * m4_options = g_strdup (M4_OPTIONS);
  GPtrArray * events = g_ptr_array_new ();
  gint maxlevel = -2;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"split", required_argument, NULL, 's'},
      {"pid", no_argument, NULL, 'i'},
      {"partition", required_argument, NULL, 'p'},
      {"profile", no_argument, NULL, 'P'},
      {"define", required_argument, NULL, 'D'},
      {"include", required_argument, NULL, 'I'},
      {"macros", no_argument, NULL, 'm'},
      {"data", no_argument, NULL, 'd'},
      {"event", required_argument, NULL, 'e'},
      {"bubble", required_argument, NULL, 'b'},
      {"debug", no_argument, NULL, 'B'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hVs:ip:PD:I:mde:b:vB",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hVs:ip:PD:I:mde:b:vB"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'P': /* profile */
      profile = TRUE;
      break;
    case 'p': /* partition */
      npart = atoi (optarg);
      break;
    case 'b': /* "bubble" partition */
      npart = atoi (optarg);
      bubble = TRUE;
      break;
    case 's': /* split */
      split = atoi (optarg);
      break;
    case 'i': /* pid */
      one_box_per_pe = FALSE;
      break;
    case 'I': { /* include */
      gchar * tmp = g_strjoin (" ", m4_options, "-I", optarg, NULL);
      g_free (m4_options);
      m4_options = tmp;
      if (set_macros ())
	return 1;
      macros = TRUE;
      break;
    }
    case 'D': { /* define */
      gchar * tmp = g_strjoin (" ", m4_options, "-D", optarg, NULL);
      g_free (m4_options);
      m4_options = tmp;
      if (set_macros ())
	return 1;
      macros = TRUE;
      break;
    }
    case 'm': /* macros */
      if (set_macros ())
	return 1;
      macros = TRUE;
      break;
    case 'd': /* data */
      maxlevel = -1;
      break;
    case 'e': /* event */
      g_ptr_array_add (events, g_strdup (optarg));
      break;
    case 'B': /* debug */
      gfs_debug_enabled (TRUE);
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': { /* help */
      gchar * usage = 
	"Usage: gerris [OPTION] FILE\n"
	"The Gerris flow solver simulation engine.\n"
	"\n"
	"  -s N   --split=N     splits the domain N times and returns\n"
	"                       the corresponding simulation\n"
	"  -i     --pid         keep box pids when splitting\n"
	"  -p N   --partition=N partition the domain in 2^N subdomains and returns\n" 
	"                       the corresponding simulation\n"
	"  -b N   --bubble=N    partition the domain in N subdomains and returns\n" 
	"                       the corresponding simulation\n"
	"  -d     --data        when splitting or partitioning, output all data\n"
	"  -P     --profile     profiles calls to boundary conditions\n"
#ifdef HAVE_M4
	"  -m     --macros      Turn macros support on\n"
	"  -DNAME               Defines NAME as a macro expanding to VALUE\n"
	"  -DNAME=VALUE         (macro support is implicitly turned on)\n"
	"         --define=NAME\n"
	"         --define=NAME=VALUE\n"
	"  -IDIR --include=DIR  Append DIR to macro include path\n"
#endif /* HAVE_M4 */
	"  -eEV   --event=EV    Evaluates GfsEvent EV and returns the simulation\n"
	"  -B     --debug       Enables debugging messages\n"
	"  -v     --verbose     Display more messages\n"
	"  -h     --help        display this help and exit\n"
	"  -V     --version     output version information and exit\n"
	"\n"
	"Reports bugs to %s\n";
      gfs_error (0, usage, FTT_MAINTAINER);
      return 0; /* success */
      break;
    }
    case 'V': { /* version */
      gchar * mpi = 
#ifdef HAVE_MPI
	"yes";
#else
        "no";
#endif
      gchar * pkgconfig =
#ifdef HAVE_PKG_CONFIG
	"yes";
#else
        "no";
#endif
      gchar * m4 =
#ifdef HAVE_M4
	"yes";
#else
        "no";
#endif
      gfs_error (0,
		 "gerris: using %dD libgfs version %s (%s)\n"
		 "  compiled with flags: %s\n"
		 "  MPI:          %s\n"
		 "  pkg-config:   %s\n"
		 "  m4:           %s\n"
		 "Copyright (C) 2001-2011 NIWA.\n"
		 "This is free software; see the source for copying conditions.  There is NO\n"
		 "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n",
		 FTT_DIMENSION,
		 GFS_VERSION,
		 GFS_BUILD_VERSION,
		 GFS_COMPILATION_FLAGS,
		 mpi,
		 pkgconfig,
		 m4);
      return 0; /* succes */
      break;
    }
    case '?': /* wrong options */
      gfs_error (0, "Try `gerris --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing FILE */
    gfs_error (0, 
	     "gerris: missing FILE\n"
	     "Try `gerris --help' for more information.\n");
    return 1; /* failure */
  }

  if (macros) {
    gchar * awk = g_strconcat ("awk -v prefix=",
			       strstr (m4_options, "-P") ? "m4_" : "",
			       " -f ", GFS_DATA_DIR, "/m4.awk ", 
			       NULL);
    gchar * command;

    if (!strcmp (argv[optind], "-"))
      command = g_strjoin (NULL, awk, "| m4 ", m4_options, NULL);
    else
      command = g_strjoin (NULL, awk, argv[optind], " | m4 ", m4_options, NULL);
    fptr = popen (command, "r");
    g_free (command);
    g_free (awk);
  }
  else { /* no macros */
    if (!strcmp (argv[optind], "-"))
      fptr = stdin;
    else {
#ifdef HAVE_MPI
      gint pid = -1;
      int size;
      
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      if (size > 1)
	MPI_Comm_rank (MPI_COMM_WORLD, &pid);
      
      gboolean dynamic = FALSE, parallel = FALSE;
      GSList * format = gfs_format_new (argv[optind], NULL, &dynamic, &parallel);
      if (dynamic) {
	gfs_error (-1, "gerris: simulation file cannot be time-dependent\n");
	return 1;
      }
      if (parallel) {
	gchar * fpname = gfs_format_string (format, pid, 0, 0.);
	fptr = fopen (fpname, "r");
	g_free (fpname);
      }
      else
	fptr = fopen (argv[optind], "r");
      gfs_format_destroy (format);
#else
      fptr = fopen (argv[optind], "r");
#endif
    }
  }
  g_free (m4_options);

  if (fptr == NULL) {
    gfs_error (-1, "gerris: unable to open file `%s'\n", argv[optind]);
    return 1;
  }

  fp = gts_file_new (fptr);
  if (!(simulation = gfs_simulation_read (fp))) {
    gfs_error (-1, 
	     "gerris: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     argv[optind], argv[optind],
	     fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);

  if (macros)
    pclose (fptr);
  else
    if (fptr != stdin)
      fclose (fptr);

  domain = GFS_DOMAIN (simulation);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    int size;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    if (size > 1 && domain->np < size) {
      if (npart != 0) {
	gfs_error (0, "gerris: manual partitioning is not a valid option in parallel\n");
	return 1;
      }
      /* automatic bubble partitioning if the simulation is not
	 partitioned correctly (or at all) */
      npart = size;
      bubble = TRUE;
    }
  }
#endif /* HAVE_MPI */

  if (split) {
    int pid = domain->pid;
    domain->pid = -1; /* force serial */
    gfs_clock_start (domain->timer);
    gfs_simulation_refine (simulation);
    gfs_clock_stop (domain->timer);
    while (split) {
      gfs_domain_split (domain, one_box_per_pe);
      split--;
    }
    if (npart == 0) {
      gfs_simulation_write (simulation, maxlevel, stdout);
      return 0;
    }
    domain->pid = pid;
  }

  if ((bubble && npart > 1) || (!bubble && npart > 0)) {
    guint nmin = 1000;
    guint mmax = 10000;
    guint ntry = 10000;
    guint np = bubble ? npart : pow (2., npart);
    gfloat imbalance = 0.0;
    GSList * partition, * i;

    if (verbose && domain->pid <= 0)
      gts_graph_print_stats (GTS_GRAPH (simulation), stderr);
    if (gts_container_size (GTS_CONTAINER (simulation)) < np) {
      gfs_error (0,
		 "gerris: the number of boxes in the domain to partition should be >= %d\n"
		 "Use option '-s' to split the domain first\n"
		 "Try `gerris --help' for more information.\n",
		 np);
      return 1;
    }
    if (bubble)
      partition = gts_graph_bubble_partition (GTS_GRAPH (simulation), npart, 100, 
					      verbose ? 
					      (GtsFunc) gts_graph_partition_print_stats : NULL, 
					      stderr);
    else
      partition = gts_graph_recursive_bisection (GTS_WGRAPH (simulation),
						 npart, 
						 ntry, mmax, nmin, imbalance);

    gint pid = 0;
    i = partition;
    while (i) {
      if (gts_container_size (GTS_CONTAINER (i->data)) == 0) {
	fprintf (stderr, "gerris: partitioning failed: empty partition\n");
	if (!bubble)
	  fprintf (stderr, 
		   "Try using the '-b' option\n"
		   "Try `gerris --help' for more information.\n");
	return 1;
      }
      gts_container_foreach (GTS_CONTAINER (i->data), (GtsFunc) set_box_pid, &pid);
      pid++;
      i = i->next;
    }

    if (pid != np)
      fprintf (stderr, "gerris: warning: only %d partitions were created\n", pid);

    if (verbose && domain->pid <= 0)
      gts_graph_partition_print_stats (partition, stderr);
    gts_graph_partition_destroy (partition);
      
    if (domain->pid >= 0) { /* we are running a parallel job */
      /* write partitioned simulation in a temporary file */
      gchar * partname = gfs_template ();
      gint fd = g_mkstemp (partname);
      remove (partname);
      g_free (partname);
      FILE * fptr = fdopen (fd, "w+");
      gfs_simulation_write (simulation, maxlevel, fptr);
      gts_object_destroy (GTS_OBJECT (simulation));
      /* replace the simulation with its partitioned version */
      rewind (fptr);
      fp = gts_file_new (fptr);
      simulation = gfs_simulation_read (fp);
      domain = GFS_DOMAIN (simulation);
      /* cleanup */
      gts_file_destroy (fp);
      fclose (fptr);
      close (fd);
      g_assert (simulation);
    }
    else { /* just a serial job */
      gfs_simulation_write (simulation, maxlevel, stdout);
      return 0;
    }
  }

  if (events->len > 0) {
    GSList * l = NULL;
    guint i;
    for (i = 0; i < events->len; i++) {
      GtsFile * fp = gts_file_new_from_string (g_ptr_array_index (events, i));
      if (fp->type != GTS_STRING) {
	gfs_error (-1, 
		   "gerris: invalid event: '%s'\n"
		   "expecting a GfsEvent name\n",
		   (char *) g_ptr_array_index (events, i));
	return 1;
      }
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      if (klass == NULL) {
	gfs_error (-1, "gerris: unknown event class `%s'\n", fp->token->str);
	return 1;
      }
      if (!gts_object_class_is_from_class (klass, gfs_event_class ())) {
	gfs_error (-1, "gerris: class `%s' is not a GfsEvent\n", fp->token->str);
	return 1;
      }
      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, simulation);
      g_assert (klass->read);
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR) {
	gfs_error (-1,
		   "gerris: invalid event: '%s'\n"
		   "%d:%d: %s\n",
		   (char *) g_ptr_array_index (events, i),
		   fp->line, fp->pos, fp->error);
	return 1;
      }
      if (GFS_IS_ADAPT (object))
	gts_container_add (GTS_CONTAINER (simulation->adapts), GTS_CONTAINEE (object));
      else if (GFS_IS_SOLID (object))
	gts_container_add (GTS_CONTAINER (simulation->solids), GTS_CONTAINEE (object));
      gts_container_add (GTS_CONTAINER (simulation->events), GTS_CONTAINEE (object));
      l = g_slist_append (l, object);
      gts_file_destroy (fp);
    }
    GtsFile * fp = gts_file_new_from_string ("");
    gfs_pending_functions_compilation (fp);
    if (fp->type == GTS_ERROR) {
      gfs_error (-1,
		 "gerris: invalid event\n"
		 "%d:%d: %s\n",
		 fp->line, fp->pos, fp->error);
      return 1;
    }
    gts_file_destroy (fp);
    gfs_clock_start (domain->timer);
    GSList * j = domain->variables;
    while (j) {
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, j->data);
      j = j->next;
    }
    g_slist_foreach (l, (GFunc) gfs_event_do, simulation);    
    gfs_clock_stop (domain->timer);
    setup_binary_IO (domain);
    gfs_simulation_write (simulation, -1, stdout);
    return 0;
  }

  domain->profile_bc = profile;

  gfs_simulation_run (simulation);

  gts_object_destroy (GTS_OBJECT (simulation));

  return 0;
}
