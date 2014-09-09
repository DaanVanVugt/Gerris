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

int main (int argc, char * argv[])
{
  int c = 0;
  gboolean verbose = FALSE;

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
     "Usage: ppmcombine [OPTION] FILE1 FILE2...\n"
     "Combines several PPM files produced by a parallel run of Gerris.\n"
     "\n"
     "  -v      --verbose     display statistics and other info\n"
     "  -h      --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `ppmcombine --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) {
    fprintf (stderr, 
	     "ppmcombine: missing FILE1\n"
	     "Try `ppmcombine --help' for more information.\n");
    return 1; /* failure */
  }

#if FTT_2D
  c = gfs_combine_ppm (&argv[optind], argc - optind, stdout);
  if (c >= 0) {
    fprintf (stderr,
	     "ppmcombine: format error in file `%s'\n",
	     argv[optind + c]);
    return 1; /* failure */
  }
#endif /* 2D only */

  return 0;
}
