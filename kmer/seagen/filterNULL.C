#include "posix.H"
#include "aHit.H"
#include "hitReader.H"

//  A NULL filter.  What comes in, comes out.  Seems useless, but the
//  hitReader merges overlapping hits which would otherwise screw up
//  your mapping.

int
main(int argc, char **argv) {
  hitReader    HR(argc);

  if (argc < 2)
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n"), exit(1);

  int arg = 1;
  while (arg < argc)
    HR.addInputFile(argv[arg++]);

  while (HR.loadHits())
    for (u32bit i=0; i < HR.numHits(); i++)
      ahit_printASCII(&HR[i].a, stdout);

  return(0);
}
