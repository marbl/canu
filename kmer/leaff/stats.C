#include "bio++.H"
#include "seqCache.H"


static
int
u32bit_compare(const void *a, const void *b) {
  const u32bit A = *((const u32bit *)a);
  const u32bit B = *((const u32bit *)b);
  if (A < B)  return(-1);
  if (A > B)  return(1);
  return(0);
}


void
stats(char *filename) {
  seqCache    *F = new seqCache(filename);

  bool                  V[256];
  for (u32bit i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

  u32bit  numSeq = F->getNumberOfSequences();

  //  sum of spans
  //  lengths of spans (for N50 computation)
  //  
  u64bit    Ss = 0;
  u32bit   *Ls = new u32bit [numSeq];

  //  sum of bases
  //  lengths of bases (for N50 computation)
  //  
  u64bit    Sb = 0;
  u32bit   *Lb = new u32bit [numSeq];

  for (u32bit i=0; i<numSeq; i++)
    Ls[i] = Lb[i] = 0;

  for (u32bit s=0; s<F->getNumberOfSequences(); s++) {
    seqInCore  *S      = F->getSequenceInCore(s);
    u32bit      len    = S->sequenceLength();
    u32bit      span   = len;
    u32bit      base   = len;

    for (u32bit pos=1; pos<len; pos++) {
      if (V[S->sequence()[pos]])
        base--;
    }

    Ss += span;
    Sb += base;

    Ls[S->getIID()] = span;
    Lb[S->getIID()] = base;

    delete S;
  }

  qsort(Ls, numSeq, sizeof(u32bit), u32bit_compare);
  qsort(Lb, numSeq, sizeof(u32bit), u32bit_compare);

  u32bit  n50s=0, n50b=0;

  for (u32bit i=0, sum=0; sum < Ss/2; i++) {
    n50s = Ls[i];
    sum += Ls[i];
  }

  for (u32bit i=0, sum=0; sum < Sb/2; i++) {
    n50b = Lb[i];
    sum += Lb[i];
  }

  fprintf(stderr, "%s\n", F->getSourceName());
  fprintf(stderr, "\n");
  fprintf(stderr, "numSeqs  "u32bitFMT"\n", F->getNumberOfSequences());
  fprintf(stderr, "\n");
  fprintf(stderr, "SPAN\n");
  fprintf(stderr, "n50      "u32bitFMTW(10)"\n", n50s);
  fprintf(stderr, "smallest "u32bitFMTW(10)"\n", Ls[0]);
  fprintf(stderr, "largest  "u32bitFMTW(10)"\n", Ls[numSeq-1]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Ss);
  fprintf(stderr, "\n");
  fprintf(stderr, "BASES\n");
  fprintf(stderr, "n50      "u32bitFMTW(10)"\n", n50b);
  fprintf(stderr, "smallest "u32bitFMTW(10)"\n", Lb[0]);
  fprintf(stderr, "largest  "u32bitFMTW(10)"\n", Lb[numSeq-1]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Sb);

  delete [] Ls;
  delete [] Lb;
}





#if 0

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

#endif
