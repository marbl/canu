#include "posix.H"
#include "searchGENOME.H"

int
main(int argc, char **argv) {

  if (argc == 0) {
  }

  u32bit  zero = 0;
  u32bit  totl = 0;

  FastABuffer    B;
  FastA         *F = new FastA(argv[1]);
  encodedQuery  *Q = 0L;

  for (F->first(B); !F->eof(); F->next(B)) {
    if ((totl & 0xfff) == 0xfff) {
      fprintf(stderr, "%9lu / %9lu\r", totl, zero);
      fflush(stderr);
    }

    Q = new encodedQuery(B.sequence(),
                         B.sequenceLength(),
                         20,
                         false);

    totl++;

    if (Q->numberOfMers() == 0) {
      zero++;
    }

    delete Q;
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "Total: %9lu\n", totl);
  fprintf(stdout, "Zero:  %9lu\n", zero);

  return(0);
}
