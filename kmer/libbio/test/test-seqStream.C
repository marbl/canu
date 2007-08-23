#include "bio++.H"
#include "seq.H"

#include "alphabet.c"

#include "fasta.C"

#include "seqInCore.C"
#include "seqOnDisk.C"
#include "seqStream.C"
//#include "seqStore.C"
#include "seqFactory.C"

int
main(int argc, char **argv) {
  seqStream    A;
  seqStream   *B = 0L;
  seqStream   *C = 0L;

  seqFile     *F = openSeqFile(argv[1]);

  B = new seqStream(argv[1], true);
  C = new seqStream(F, true);

  A.setFile(argv[1]);
  A.parse("0,2,1");
  A.finish();

  A.setSeparator('.', 10);

  fprintf(stderr, "A: numSeqs "u32bitFMT" len "u64bitFMT"\n", A.numberOfSequences(), A.lengthOfSequences());
  fprintf(stderr, "B: numSeqs "u32bitFMT" len "u64bitFMT"\n", B->numberOfSequences(), B->lengthOfSequences());
  fprintf(stderr, "C: numSeqs "u32bitFMT" len "u64bitFMT"\n", C->numberOfSequences(), C->lengthOfSequences());

  for (u32bit i=0; i<A.numberOfSequences(); i++) {
    printf("seq:"u32bitFMT" len:"u32bitFMT" start:"u32bitFMT" iid:"u32bitFMT"\n",
           i, A.lengthOf(i), A.startOf(i), A.IIDOf(i));
  }

  char   ch = 0;
  u32bit xx = 0;
  while (ch = A.get()) {
    u64bit sp = A.seqPos();
    u64bit si = A.seqIID();
    u64bit st = A.strPos();
    printf("%c -- seqPos:"u64bitFMT" seqIID:"u64bitFMT" strPos:"u64bitFMT" lookup:"u32bitFMT"\n",
           ch, sp, si, st,
           A.sequenceNumberOfPosition(xx));
    xx++;
  }

  A.rewind();

  while (!A.eof()) {
    printf("%c", A.get());
  }
  printf("\n");

  delete C;
  delete B;
}
