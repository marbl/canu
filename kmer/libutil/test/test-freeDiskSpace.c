#include "util.h"

int
main(int argc, char **argv) {
  int i;

  if (argc == 1) {
    fprintf(stderr, "usage: %s file [...]\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc; i++)
    fprintf(stderr, "%s: %d\n", argv[i], (int)freeDiskSpace(argv[i]));

  return(0);
}
