#include "bio++.H"
#include "encodedQuery.H"

int
main(int argc, char **argv) {

  if (argc == 1) {
    mt_s  *mt = mtInit(time(0L));

    fprintf(stderr, "Building random sequences for testing.\n");

    for (u32bit i=0; i<100000; i++) {
      char *seq = new char [10000];
      char *hdr = new char [128];

      for (u32bit j=0; j<10000; j++) {
        seq[j] = decompressSymbol[mtRandom32(mt) % 4];
        if (mtRandomRealOpen(mt) < 0.01)
          seq[j] = 'n';
      }
      seq[9999] = 0;

      sprintf(hdr, ">"u32bitFMT, i);

      seqInCore            *S = new seqInCore(i, hdr, strlen(hdr), seq, 9999);
      encodedQuery         *Q = new encodedQuery(S, 22);
      Q->test(S);
      delete Q;
      delete S;
    }

  } else {
    seqFile *F = openSeqFile(argv[1]);

    while (F->eof() == false) {
      seqInCore            *S = F->getSequenceInCore();
      encodedQuery         *Q = new encodedQuery(S, 22);
      Q->test(S);
      delete Q;
      delete S;
    }

    delete F;
  }


  exit(0);
}
