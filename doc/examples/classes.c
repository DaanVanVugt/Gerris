#include <string.h>
#include "init.h"

static void key_value_pair (const char * key)
{
  printf ("%s\n", key);
  /* keywords must start with Gfs */
  g_assert (strstr (key, "Gfs") == key);
  printf ("%s\n", &(key[3]));
}

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass;

  klass = gfs_classes ();

  key_value_pair ("GfsDefine");
  key_value_pair ("GfsTime");
  key_value_pair ("GfsDeferredCompilation");
  key_value_pair ("GfsProjectionParams");
  key_value_pair ("GfsApproxProjectionParams");
  key_value_pair ("GfsPhysicalParams");
  key_value_pair ("GfsAdvectionParams");

  while (*klass) {
    key_value_pair ((*klass)->info.name);
    klass++;
  }
  
  return 0;
}
