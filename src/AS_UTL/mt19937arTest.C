#include "mt19937ar.H"

int
main(int argc, char **argv) {
  mtRandom  mt;

  if (argc != 4)
    fprintf(stderr, "usage: %s <iterations> <lambda> <rho>\n", argv[0]), exit(1);

  uint32  number  = atoi(argv[1]);
  double  mode    = atof(argv[2]);
  double  scale   = atof(argv[3]);

  for (uint32 ii=0; ii<number; ii++)
    fprintf(stdout, "%f\n", mt.mtRandomExponential(mode, scale));

  exit(0);
}

