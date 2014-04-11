#include "bio++.H"
#include "seqCache.H"

void
dumpBlocks(char *filename) {
seqCache      *F     = 0L;
  seqInCore   *S     = 0L;

  bool                  V[256] = {0};

  for (uint32 i=0; i<256; i++)
    V[i] = false;

  V['n'] = true;
  V['N'] = true;

  F = new seqCache(filename);

  for (uint32 s=0; s<F->getNumberOfSequences(); s++) {
    seqInCore *S = F->getSequenceInCore(s);

    uint32  len    = S->sequenceLength();
    char    begseq = S->sequence()[0];
    bool    nnn    = V[begseq];
    uint32  begpos = 0;
    uint32  pos    = 0;

    for (pos=0; pos<len; pos++) {
      char seq = S->sequence()[pos];

      if (nnn != V[seq]) {
        fprintf(stdout, "%c "uint32FMT" "uint32FMT" "uint32FMT" "uint32FMT"\n",
                begseq, s, begpos, pos, pos - begpos);
        nnn = V[seq];
        begpos = pos;
        begseq = seq;
      }
    }

    fprintf(stdout, "%c "uint32FMT" "uint32FMT" "uint32FMT" "uint32FMT"\n",
            begseq, s, begpos, pos, pos - begpos);
    fprintf(stdout, ". "uint32FMT" "uint32FMT" "uint32FMT"\n", s, pos, uint32ZERO);

    delete S;
  }

  delete F;
}


