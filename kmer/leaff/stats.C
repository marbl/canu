#include "bio++.H"
#include "seqCache.H"

#include <algorithm>

using namespace std;


void
stats(char *filename, uint64 refLen) {
  seqCache    *F = new seqCache(filename);

  bool                  V[256];
  for (uint32 i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

  uint32  numSeq = F->getNumberOfSequences();

  uint64    Ss = 0;  //  actual length of span
  uint64    Rs = 0;  //  reference length of span
  uint32   *Ls = new uint32 [numSeq];

  uint64    Sb = 0;
  uint64    Rb = 0;
  uint32   *Lb = new uint32 [numSeq];

  for (uint32 i=0; i<numSeq; i++)
    Ls[i] = Lb[i] = 0;

  for (uint32 s=0; s<numSeq; s++) {
    seqInCore  *S      = F->getSequenceInCore(s);
    uint32      len    = S->sequenceLength();
    uint32      span   = len;
    uint32      base   = len;

    for (uint32 pos=1; pos<len; pos++) {
      if (V[S->sequence()[pos]])
        base--;
    }

    Ss += span;
    Sb += base;

    Ls[S->getIID()] = span;
    Lb[S->getIID()] = base;

    delete S;
  }

  if (refLen > 0) {
    Rs = refLen;
    Rb = refLen;
  } else {
    Rs = Ss;
    Rb = Sb;
  }

  //qsort(Ls, numSeq, sizeof(uint32), uint32_compare);
  //qsort(Lb, numSeq, sizeof(uint32), uint32_compare);

  sort(Ls, Ls + numSeq);
  sort(Lb, Lb + numSeq);

  reverse(Ls, Ls + numSeq);
  reverse(Lb, Lb + numSeq);

  uint32  n50s[11] = {0};
  uint32  l50s[11] = {0};

  uint32  n50b[11] = {0};
  uint32  l50b[11] = {0};

  uint32  sizes[11] = {0};
  uint32  sizeb[11] = {0};

  for (uint32 i=0; i<11; i++) {
    sizes[i] = i * Rs / 10;
    sizeb[i] = i * Rb / 10;
    //fprintf(stderr, "SIZE %2d  s=%d b=%d\n", i, sizes[i], sizeb[i]);
  }

  for (uint32 i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizes[n]) && (sizes[n] <= sum + Ls[i])) {
      n50s[n]  = Ls[i];
      l50s[n]  = i;
      n++;
    }

    sum += Ls[i];
  }


  for (uint32 i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizeb[n]) && (sizeb[n] <= sum + Lb[i])) {
      n50b[n]  = Ls[i];
      l50b[n]  = i;
      n++;
    }

    sum += Lb[i];
  }

  //for (uint32 i=0, sum=0; sum < Rb/2; i++) {
  //}

  fprintf(stdout, "%s\n", F->getSourceName());
  fprintf(stdout, "\n");
  fprintf(stdout, "numSeqs  "uint32FMT"\n", numSeq);
  fprintf(stdout, "\n");
  fprintf(stdout, "SPAN (smallest "uint32FMT" largest "uint32FMT")\n", Ls[numSeq-1], Ls[0]);
  for (uint32 i=1; i<10; i++)
    fprintf(stdout, "n"uint32FMT"      "uint32FMT" at index "uint32FMT"\n", 10 * i, n50s[i], l50s[i]);
  fprintf(stdout, "totLen   "uint64FMTW(10)"\n", Ss);
  fprintf(stdout, "refLen   "uint64FMTW(10)"\n", Rs);
  fprintf(stdout, "\n");
  fprintf(stdout, "BASES (smallest "uint32FMT" largest "uint32FMT")\n", Lb[numSeq-1], Lb[0]);
  for (uint32 i=1; i<10; i++)
    fprintf(stdout, "n"uint32FMT"      "uint32FMT" at index "uint32FMT"\n", 10 * i, n50b[i], l50b[i]);
  fprintf(stdout, "totLen   "uint64FMTW(10)"\n", Sb);
  fprintf(stdout, "refLen   "uint64FMTW(10)"\n", Rb);

  delete [] Ls;
  delete [] Lb;
}
