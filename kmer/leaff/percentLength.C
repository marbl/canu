#include "bio++.H"


struct info {
  u32bit  realLength;
  u32bit  acgtLength;
  u32bit  iid;
  char    name[32];
};


int
info_compare(const void *a, const void *b) {
  const info  *A = (const info *)a;
  const info  *B = (const info *)b;

  if (A->realLength < B->realLength)
    return(-1);
  if (A->realLength > B->realLength)
    return(1);
  return(0);
}


int
main(int argc, char **argv) {
  char *filename = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      filename = argv[++arg];
    } else if (strcmp(argv[arg], "-x") == 0) {

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (filename == 0L)) {
    fprintf(stderr, "usage: %s -f some.fasta\n", argv[0]);
    exit(1);
  }

  seqFile *F = openSeqFile(filename);
  F->openIndex();

  u32bit    N = F->getNumberOfSequences();

  info     *I = new info [N];

  u64bit    realLengthTotal = 0;
  u64bit    acgtLengthTotal = 0;

  seqOnDisk  *seq;

  while ((seq = F->getSequenceOnDisk()) != 0L) {
    seqIID      iid = seq->getIID();

    strncpy(I[iid].name, seq->header(), 32);
    for (u32bit x=0; x<32; x++)
      if (isspace(I[iid].name[x]))
        I[iid].name[x] = 0;

    I[iid].realLength = seq->sequenceLength();
    I[iid].acgtLength = validSymbol[seq->get()];

    while (seq->next())
      I[iid].acgtLength += validSymbol[seq->get()];

    realLengthTotal += I[iid].realLength;
    acgtLengthTotal += I[iid].acgtLength;

    delete seq;
  }

  delete F;

  qsort(I, N, sizeof(info), info_compare);

  u64bit  realCumulative = 0;
  u64bit  acgtCumulative = 0;

  for (u32bit i=0; i<N; i++) {
    realCumulative += I[i].realLength;
    acgtCumulative += I[i].acgtLength;

    fprintf(stdout, "%s\t"u32bitFMT"\t%6.2f\t"u32bitFMT"\t%6.2f\n",
            I[i].name,
            I[i].realLength, 100.0 * realCumulative / realLengthTotal,
            I[i].acgtLength, 100.0 * acgtCumulative / acgtLengthTotal);
  }
}
