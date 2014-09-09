#include <stdio.h>
#include "boyd87.h"

int main (int argc, char * argv[])
{
  double HW, TW, B, D, S0, L, N, KE, G;
  int type;
  B = 1.0;
  D = 1.0;
  type = 1;
  S0 = 0.01;
  L = 10.0;
  N = 0.013;
  KE = 1.0;
  G = 9.81;

  printf ("B D S0 L N KE G\n");
  printf ("%g %g %g %g %g %g %g\n", B, D, S0, L, N, KE, G);
  printf ("type\n");
  printf ("%d\n", type);
  printf ("HW TW Q_Box Q_Pipe\n");
  for (TW = 0.; TW <= 1.5*D; TW += 0.001) {
    for (HW = 0.1; HW <= 20.; HW += 0.1)
      /*printf ("--->%g %g %g\n", HW, TW, Q_box (HW, TW, B, D, type, S0, L, N, KE, G));*/
      printf ("--->%g %g %g\n", HW, TW, Q_pipe (HW, TW, D, type, S0, L, N, KE, G));
    /*printf ("%g %g %g %g\n", HW, TW, Q_box (HW, TW, B, D, type, S0, L, N, KE, G), Q_pipe (HW, TW, D, type, S0, L, N, KE, G));*/
    printf ("\n");
  }
  return 0;
}
