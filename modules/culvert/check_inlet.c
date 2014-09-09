#include <stdio.h>
#include "boyd87.h"

int main (int argc, char * argv[])
{
  double HW, D;
  for (D = 0.6; D <= 3.; D += 0.6) {
    for (HW = 0.1; HW <= 20.; HW += 0.1)
      printf ("%g %g %g %g\n", D, HW, Q_inlet_box (HW, 1., D, 1, 9.81), Q_inlet_pipe (HW, D, 1, 9.81));
    printf ("\n");
  }
  return 0;
}
