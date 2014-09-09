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

#include "init.h"
#include "simulation.h"
#include "graphic.h"

static void streamline_draw (GList * s, FILE * fp)
{
  guint np = g_list_length (s);

  fprintf (fp, "VECT 1 %u 0 %u 0\n", np, np);
  while (s) {
    GtsPoint * p = s->data;

    fprintf (fp, "%g %g %g\n", p->x, p->y, p->z);
    s = s->next;
  }
}

int main (int argc, char * argv[])
{
  int c = 0;
  GtsFile * fp;
  GfsTime time;
  guint ns = 0;

  gboolean verbose = FALSE;

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
     "Usage: streamanime [OPTION] < STREAMLINE_FILE\n"
     "Converts a Gerris streamline file to other (graphical) formats.\n"
     "\n"
     "  -v      --verbose     display statistics and other info\n"
     "  -h      --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `streamanime --help' for more information.\n");
      return 1; /* failure */
    }
  }

  fp = gts_file_new (stdin);
  gfs_time_init (&time);
  while (fp->type == GTS_STRING) {
    if (!strcmp (fp->token->str, "GfsTime")) {
      gts_file_next_token (fp);
      gfs_time_read (&time, fp);
      if (fp->type == GTS_ERROR) {
	fprintf (stderr, 
		 "streamanime: file on standard input is not a valid streamline file\n"
		 "<stdin>:%d:%d: %s\n",
		 fp->line, fp->pos, fp->error);
	return 1;
      }
      if (verbose)
	fprintf (stderr, "\rstreamanime: processing t: %7.3f n: %5u", 
		 time.t, ns);
      gts_file_first_token_after (fp, '\n');
      ns = 0;
      printf ("(redraw focus)\n(freeze focus)\n");
    }
    else if (!strcmp (fp->token->str, "GfsStreamline")) {
      GList * streamline = gfs_streamline_read (fp);

      if (fp->type == GTS_ERROR) {
	fprintf (stderr, 
		 "streamanime: file on standard input is not a valid streamline file\n"
		 "<stdin>:%d:%d: %s\n",
		 fp->line, fp->pos, fp->error);
	return 1;
      }
      printf ("(geometry \"stream-%u\" = {\n", ns++);
      streamline_draw (streamline, stdout);
      printf ("})\n");
      gfs_streamline_destroy (streamline);
    }
    else {
      gts_file_error (fp, "unknown identifier `%s'", fp->token->str);
      fprintf (stderr, 
	       "streamanime: file on standard input is not a valid streamline file\n"
	       "<stdin>:%d:%d: %s\n",
	       fp->line, fp->pos, fp->error);
      return 1;
    }
  }
  gts_file_destroy (fp);
  printf ("(redraw focus)\n");
  if (verbose)
    fputc ('\n', stderr);

  return 0;
}
