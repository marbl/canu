#include "bio++.H"
#include "seqCache.H"

#include <algorithm>

using namespace std;


void
stats(char *filename, u64bit refLen) {
  seqCache    *F = new seqCache(filename);

  bool                  V[256];
  for (u32bit i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

  u32bit  numSeq = F->getNumberOfSequences();

  u64bit    Ss = 0;  //  actual length of span
  u64bit    Rs = 0;  //  reference length of span
  u32bit   *Ls = new u32bit [numSeq];

  u64bit    Sb = 0;
  u64bit    Rb = 0;
  u32bit   *Lb = new u32bit [numSeq];

  for (u32bit i=0; i<numSeq; i++)
    Ls[i] = Lb[i] = 0;

  for (u32bit s=0; s<numSeq; s++) {
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

  if (refLen > 0) {
    Rs = refLen;
    Rb = refLen;
  } else {
    Rs = Ss;
    Rb = Sb;
  }

  //qsort(Ls, numSeq, sizeof(u32bit), u32bit_compare);
  //qsort(Lb, numSeq, sizeof(u32bit), u32bit_compare);

  sort(Ls, Ls + numSeq);
  sort(Lb, Lb + numSeq);

  reverse(Ls, Ls + numSeq);
  reverse(Lb, Lb + numSeq);

  u32bit  n50s[11] = {0};
  u32bit  l50s[11] = {0};

  u32bit  n50b[11] = {0};
  u32bit  l50b[11] = {0};

  u32bit  sizes[11] = {0};
  u32bit  sizeb[11] = {0};

  for (u32bit i=0; i<11; i++) {
    sizes[i] = i * Rs / 10;
    sizeb[i] = i * Rb / 10;
    fprintf(stderr, "SIZE %2d  s=%d b=%d\n", i, sizes[i], sizeb[i]);
  }

  for (u32bit i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizes[n]) && (sizes[n] <= sum + Ls[i])) {
      n50s[n]  = Ls[i];
      l50s[n]  = i;
      n++;
    }

    sum += Ls[i];
  }


  for (u32bit i=0, sum=0, n=1; (i < numSeq) && (n < 11); i++) {
    if ((sum <  sizeb[n]) && (sizeb[n] <= sum + Lb[i])) {
      n50b[n]  = Ls[i];
      l50b[n]  = i;
      n++;
    }

    sum += Lb[i];
  }

  //for (u32bit i=0, sum=0; sum < Rb/2; i++) {
  //}

  fprintf(stderr, "%s\n", F->getSourceName());
  fprintf(stderr, "\n");
  fprintf(stderr, "numSeqs  "u32bitFMT"\n", numSeq);
  fprintf(stderr, "\n");
  fprintf(stderr, "SPAN (smallest "u32bitFMT" largest "u32bitFMT")\n", Ls[numSeq-1], Ls[0]);
  for (u32bit i=1; i<10; i++)
    fprintf(stderr, "n"u32bitFMT"      "u32bitFMT" at index "u32bitFMT"\n", 10 * i, n50s[i], l50s[i]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Ss);
  fprintf(stderr, "refLen   "u64bitFMTW(10)"\n", Rs);
  fprintf(stderr, "\n");
  fprintf(stderr, "BASES (smallest "u32bitFMT" largest "u32bitFMT")\n", Lb[numSeq-1], Lb[0]);
  for (u32bit i=1; i<10; i++)
    fprintf(stderr, "n"u32bitFMT"      "u32bitFMT" at index "u32bitFMT"\n", 10 * i, n50b[i], l50b[i]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Sb);
  fprintf(stderr, "refLen   "u64bitFMTW(10)"\n", Rb);

  delete [] Ls;
  delete [] Lb;
}
