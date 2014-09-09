#include <stdio.h>
#include "boyd87.h"

int main (int argc, char * argv[])
{
  double TW, D;
  for (D = 0.6; D <= 3.; D += 0.6) {
    for (TW = 0.; TW <= 1.5*D; TW += 0.001)
      printf ("%g %g %g %g\n", D, TW, 
	      Q_outlet_box  (1.5*D, TW, 1., D, 0., 10., 0.01, 0., 9.81),
	      Q_outlet_pipe (1.5*D, TW,     D, 0., 10., 0.01, 0., 9.81));
    printf ("\n");
  }
  return 0;
}
