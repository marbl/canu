#include "bio++.H"
#include "merTable.H"


int
main(int argc, char **argv) {
  merTable X;

  chainedSequence *CS = new chainedSequence();
  CS->setSource(argv[1]);
  CS->finish();

  X.build(CS, 8);
}

