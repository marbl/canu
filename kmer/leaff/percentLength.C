#include "bio++.H"
#include "seqCache.H"

bool  sortOnN  = false;

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

  if (sortOnN) {
    if (A->realLength > B->realLength)
      return(-1);
    if (A->realLength < B->realLength)
      return(1);
  } else {
    if (A->acgtLength > B->acgtLength)
      return(-1);
    if (A->acgtLength < B->acgtLength)
      return(1);
  }

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
    } else if (strcmp(argv[arg], "-n") == 0) {
      sortOnN = true;
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

  seqCache *F = new seqCache(filename);
  u32bit    N = F->getNumberOfSequences();
  info     *I = new info [N];

  u64bit    realLengthTotal = 0;
  u64bit    acgtLengthTotal = 0;

  u32bit      iid = 0;
  seqInCore  *seq = F->getSequenceInCore(iid);

#warning old style seq iteration

  while (seq != 0L) {
    strncpy(I[iid].name, seq->header() + (seq->header()[0] == '>'), 32);
    for (u32bit x=0; x<32; x++)
      if (isspace(I[iid].name[x]))
        I[iid].name[x] = 0;

    I[iid].realLength = seq->sequenceLength();

    for (u32bit xx=0; xx<seq->sequenceLength(); xx++)
      I[iid].acgtLength += (letterToBits[seq->sequence()[xx]] != 0xff);

    realLengthTotal += I[iid].realLength;
    acgtLengthTotal += I[iid].acgtLength;

    delete seq;
    seq = F->getSequenceInCore(++iid);
  }

  delete F;

  qsort(I, N, sizeof(info), info_compare);

  u64bit  realCumulative = 0;
  u64bit  acgtCumulative = 0;

  fprintf(stdout, "name\tlength\tcumulativeLength\t%% covered by seq >= given length\n");

  for (u32bit i=0; i<N; i++) {
    realCumulative += I[i].realLength;
    acgtCumulative += I[i].acgtLength;

    if (sortOnN)
      fprintf(stdout, "%s\t"u32bitFMT"\t"u64bitFMT"\t%f\n",
              I[i].name,
              I[i].realLength,
              realCumulative,
              100.0 * realCumulative / realLengthTotal);
    else
      fprintf(stdout, "%s\t"u32bitFMT"\t"u64bitFMT"\t%f\n",
              I[i].name,
              I[i].acgtLength,
              acgtCumulative,
              100.0 * acgtCumulative / acgtLengthTotal);
  }
}
