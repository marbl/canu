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

  while (!seqs.eof()) {
    FastASequenceInCore *aSeq = seqs.getSequence();

    fprintf(stderr, "%8d  '%s'\n", aSeq->getIID(), aSeq->header());

    delete aSeq;
  }

  fprintf(stderr, "\n");

  seqs.openIndex();
  for (int id=seqs.getNumberOfSequences() - 1; id>=0; id--) {
    seqs.find(id);

    FastASequenceOnDisk *aSeq = seqs.getSequenceOnDisk();

    fprintf(stderr, "%8d  '%s'\n", aSeq->getIID(), aSeq->header());

    delete aSeq;
  }
}

