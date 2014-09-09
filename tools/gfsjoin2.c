/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2009 National Institute of Water and Atmospheric Research
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
#include <unistd.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include "init.h"
#include "simulation.h"

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi_boundary.h"
#endif /* HAVE_MPI */

static void add_box (GfsBox * box, GfsSimulation * sim)
{
  gts_container_add (GTS_CONTAINER (sim), GTS_CONTAINEE (box));
}

static void add_id (GfsBox * box, GPtrArray * ids)
{
  if (box->id > ids->len)
    g_ptr_array_set_size (ids, box->id);
  g_ptr_array_index (ids, box->id - 1) = box;
}

static void convert_boundary_mpi_into_edges (GfsBox * box, GPtrArray * ids)
{
#ifdef HAVE_MPI
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundaryMpi * b = GFS_BOUNDARY_MPI (box->neighbor[d]);
      if (b->id < 0) {
	fprintf (stderr, 
		 "gfsjoin: id < 0, you maybe trying to join old parallel simulation files\n");
	exit (1);
      }
      GfsBox * nbox = g_ptr_array_index (ids, b->id - 1);
      if (nbox) {
	if (nbox->pid != b->process) {
	  fprintf (stderr, "gfsjoin: inconsistent MPI boundary: pid = %d, nbox->pid = %d\n",
		   b->process, nbox->pid);
	  exit (1);
	}
	if (!GFS_IS_BOUNDARY_MPI (nbox->neighbor[FTT_OPPOSITE_DIRECTION (d)])) {
	  fprintf (stderr, "gfsjoin: inconsistent MPI boundary: nbox[%d] is not an MPI boundary\n",
		   FTT_OPPOSITE_DIRECTION (d));
	  exit (1);
	}
	GfsBoundaryMpi * nb = GFS_BOUNDARY_MPI (nbox->neighbor[FTT_OPPOSITE_DIRECTION (d)]);
	if (box->pid != nb->process || box->id != nb->id) {
	  fprintf (stderr, "gfsjoin: inconsistent MPI boundary\n"
		   "box->pid != nb->process || box->id != nb->id\n");
	  exit (1);
	}
	gts_object_destroy (GTS_OBJECT (b));
	gts_object_destroy (GTS_OBJECT (nb));
	gfs_gedge_new (gfs_gedge_class (), box, nbox, d);
      }
    }
#endif /* HAVE_MPI */
}

int main (int argc, char * argv[])
{
  int c = 0;
  gboolean verbose = FALSE, keep = FALSE;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"keep", no_argument, NULL, 'k'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, 
			      "hvk",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, 
			 "hvk"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'k': /* keep */
      keep = TRUE;
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
	       "Usage: gfsjoin [OPTION] NAME1 NAME2... > JOINED\n"
	       "Joins several parallel Gerris simulation files\n"
	       "\n"
	       "  -k      --keep        keep MPI boundaries\n"
	       "  -v      --verbose     display statistics and other info\n"
	       "  -h      --help        display this help and exit\n"
	       "\n"
	       "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfsjoin --help' for more information.\n");
      return 1; /* failure */
    }
  }

#ifndef HAVE_MPI
  fprintf (stderr, 
	   "gfsjoin: only works for MPI version of Gerris\n"
	   "Try `gfsjoin --help' for more information.\n");
  return 1; /* failure */
#endif

  if (optind >= argc) { /* missing NAME */
    fprintf (stderr, 
	     "gfsjoin: missing NAME\n"
	     "Try `gfsjoin --help' for more information.\n");
    return 1; /* failure */
  }

  GfsSimulation ** sim = g_malloc (sizeof (GfsSimulation *)*(argc - optind));
  for (c = optind; c < argc; c++) {
    FILE * fptr = fopen (argv[c], "r");
    if (fptr == NULL) {
      fprintf (stderr, "gfsjoin: cannot open file `%s'\n", argv[c]);
      return 1; /* failure */
    }
    GtsFile * fp = gts_file_new (fptr);    
    if (!(sim[c - optind] = gfs_simulation_read (fp))) {
      fprintf (stderr, 
	       "gfsjoin: file `%s' is not a valid simulation file\n"
	       "%s:%d:%d: %s\n",
	       argv[c], argv[c], fp->line, fp->pos, fp->error);
      return 1;
    }
    if (verbose)
      fprintf (stderr, "%s: %d box(es), %d cells\n", argv[c],
	       gts_container_size (GTS_CONTAINER (sim[c - optind])),
	       gfs_domain_size (GFS_DOMAIN (sim[c - optind]), FTT_TRAVERSE_LEAFS, -1));
  }

  /* Add all boxes to first simulation */
  for (c = 1; c < argc - optind; c++)
    gts_container_foreach (GTS_CONTAINER (sim[c]), (GtsFunc) add_box, sim[0]);

  if (!keep) {
    /* Create array for fast linking of ids to GfsBox pointers */
    GPtrArray * ids = g_ptr_array_new ();
    gts_container_foreach (GTS_CONTAINER (sim[0]), (GtsFunc) add_id, ids);
    
    /* Convert GfsBoundaryMpi into graph edges */
    gts_container_foreach (GTS_CONTAINER (sim[0]), (GtsFunc) convert_boundary_mpi_into_edges, ids);

    g_ptr_array_free (ids, TRUE);
  }

  gfs_simulation_write (sim[0], -1, stdout);

  return 0;
}
