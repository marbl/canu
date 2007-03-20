#include <stdio.h>
#include <stdlib.h>

#include "fasta.H"


int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s <query-filename>\n", argv[0]);
    exit(1);
  }

  FastAFile  seqs(argv[1]);

#if 0
  while (!seqs.eof()) {
    FastASequenceInCore *aSeq = seqs.getSequence();

    fprintf(stdout, "%s\n%s\n\n", aSeq->header(), aSeq->sequence());

    delete aSeq;
  }
#endif


#if 1
  seqs.openIndex();
  for (int id=seqs.getNumberOfSequences() - 1; id>=0; id--) {
    seqs.find(id);

    FastASequenceOnDisk *aSeq = seqs.getSequenceOnDisk();

    fprintf(stdout, "%s\n", aSeq->header());

#if 0
    //  Test the char-by-char method
    //
    while (aSeq->get()) {
      fputc(aSeq->get(), stdout);
      aSeq->next();
    }
    fputc('\n', stdout);
#endif

#if 0
    //  Test the block-copy method
    //
    char *b = new char [aSeq->sequenceLength() + 1];
    aSeq->getChars(b, 0, aSeq->sequenceLength());
    fprintf(stdout, "%s\n", b);
    delete [] b;
#endif

#if 1
    //  Test the block-copy method on large files
    //
    char *b = new char [80 + 1];
    aSeq->getChars(b, 0, 80);
    fprintf(stdout, "%s\n", b);
    delete [] b;
#endif


    delete aSeq;
  }
#endif
}

