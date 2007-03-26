#include <stdio.h>
#include <stdlib.h>

#include "fasta.H"


int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s <query-filename>\n", argv[0]);
    exit(1);
  }

  seqFile *seqs = openSeqFile(argv[1]);

  while (!seqs.eof()) {
    seqInCore *aSeq = seqs.getSequenceInCore();

    fprintf(stderr, "%8d  '%s'\n", aSeq->getIID(), aSeq->header());

    delete aSeq;
  }

  fprintf(stderr, "\n");

  seqs.openIndex();
  for (int id=seqs.getNumberOfSequences() - 1; id>=0; id--) {
    seqs.find(id);

    seqOnDisk *aSeq = seqs.getSequenceOnDisk();

    fprintf(stderr, "%8d  '%s'\n", aSeq->getIID(), aSeq->header());

    delete aSeq;
  }

  delete seqs;
}

