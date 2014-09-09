/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
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

// #include <glib.h>
#include <stdlib.h>
#include <stdio.h>

#include "kdt.h"

int main (int argc, char * argv[])
{
  if (argc != 2) {
    fprintf (stderr, "Usage: %s basename\n", argv[0]);
    return -1;
  }

  Kdt * kdt = kdt_new ();
  if (kdt_open (kdt, argv[1])) {
    fprintf (stderr, "%s: could not open `%s'\n", argv[0], argv[1]);
    return -1;
  }

  KdtRect query;
  int count = 0;
  //  GTimer * t = g_timer_new ();
  while (scanf ("%f %f %f %f", 
		&query[0].l, &query[1].l, 
		&query[0].h, &query[1].h) == 4) {
#if 1
    fprintf (stderr, "%ld\n", kdt_query (kdt, query));
#else
    KdtSum s;
    kdt_sum_init (&s);
    //    g_timer_start (t);
    long n = kdt_query_sum (kdt, (KdtCheck) kdt_includes, (KdtCheck) kdt_intersects, 
			    query, query, &s);
    //    g_timer_stop (t);
    //    fprintf (stderr, "%d %g %g %g %g\n", n, s.H0, s.Hmax, s.Hmin, g_timer_elapsed (t, NULL));
    printf ("%ld %g %g %g\n", n, s.H0, s.H1, s.H2);
#endif
    if (count > 0 && count % 1000 == 0)
      fprintf (stderr, "\r%d", count);
    count++;
  }
  if (count >= 1000)
    fputc ('\n', stderr);

  kdt_destroy (kdt);
  return 0;
}
