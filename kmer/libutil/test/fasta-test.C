#include "libbri.H"

//  Tests the fasta index by opening a file twice, once with an index
//  and once without, then reading all the stuff using both methods.
//

int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    exit(1);
  }

  fprintf(stderr, "Opening for sequential access.\n");
  FastABuffer  Sb;
  FastA        S(argv[1], false, true);

  fprintf(stderr, "Opening for random access.\n");
  FastABuffer  Rb;
  FastA        R(argv[1], true, true);

  u32bit i = 0;;
  u64bit p = 0;

  S.first(Sb);
  printf("%s\n", Sb.sequence());
  printf("%d\n", Sb.sequenceLength());

  exit(1);

  for (S.first(Sb); !S.eof(); S.next(Sb)) {
    if ((i & 0xfff) == 0x000) {
      fprintf(stderr, "testing %u (%lu)\r", i, p);
      fflush(stderr);
    }

    R.seek(Rb, i);

    if (Sb.headerLength() != Rb.headerLength()) {
      fprintf(stderr, "Header Length Error in sequence %u\n", i);
    }

    if (Sb.sequenceLength() != Rb.sequenceLength()) {
      fprintf(stderr, "Sequence Length Error in sequence %u\n", i);
    }

    if (strcmp((char *)Sb.header(), (char *)Rb.header()) != 0) {
      fprintf(stderr, "Header Error in sequence %u\n", i);
    }

    if (strcmp((char *)Sb.sequence(), (char *)Rb.sequence()) != 0) {
      fprintf(stderr, "Sequence Error in sequence %u\n", i);
    }

    p += Rb.sequenceLength();

    i++;
  }
}
