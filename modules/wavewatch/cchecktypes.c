/* Simple checks for compatibility of C and Fortran types */

#include <stdio.h>

extern void check_int_ (void * a, void * b, void * c);
extern void check_float_ (void * a, void * b, void * c);

int main (int argc, char * argv[])
{
  int ia, ib, ic;
  long long lc, la, lb;
  float fa, fb, fc;
  double dc, da, db;

  printf ("/* Automatically generated using 'checktypes' */\n");
  ia = -123; ib = 875;
  check_int_ (&ia, &ib, &ic);
  if (ic == ia*ib)
    printf ("typedef int INTEGER;\n");
  else {
    la = -123; lb = 875;
    check_int_ (&la, &lb, &lc);
    if (lc != la*lb) {
      fprintf (stderr, "checktypes: could not find compatible C and Fortran integers\n");
      return 1;
    }
    printf ("typedef long long INTEGER;\n");
  }

  fa = -123.; fb = 875.;
  check_float_ (&fa, &fb, &fc);
  if (fc == fa*fb)
    printf ("typedef float REAL;\n");
  else {
    da = -123.; db = 875.;
    check_float_ (&da, &db, &dc);
    if (dc != da*db) {
      fprintf (stderr, "checktypes: could not find compatible C and Fortran floats\n");
      return 1;
    }
    printf ("typedef double REAL;\n");
  }

  return 0;
}
