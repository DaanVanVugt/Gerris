#include <string.h>
#include "init.h"

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass;
  int status = 0;

  klass = gfs_classes ();

  printf ("/** \\file\n"
	  " * \\brief Class hierarchy definition for doxygen */\n");

  char s[80];
  while (fgets (s, 80, stdin)) {
    char * file = strtok (s, " ");
    char * name = strtok (NULL, "\n");
    GtsObjectClass * klass = gts_object_class_from_name (name);
    if (klass == NULL) {
      fprintf (stderr, "autoclass:%s: unknown class '%s'\n", file, name);
      status = 1;
    }
    else
      printf ("/** \n"
	      " * \\defgroup %s %s\n"
	      " * \\ingroup %s\n"
	      " */\n", 
	      klass->info.name, klass->info.name, klass->parent_class->info.name);
  }

  return status;
}
