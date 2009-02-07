#include "bio++.H"
#include "seqCache.H"

void
dumpBlocks(char *filename) {
seqCache      *F     = 0L;
  seqInCore   *S     = 0L;

  bool                  V[256] = {0};

  for (u32bit i=0; i<256; i++)
    V[i] = false;

  V['n'] = true;
  V['N'] = true;

  F = new seqCache(filename);

  for (u32bit s=0; s<F->getNumberOfSequences(); s++) {
    seqInCore *S = F->getSequenceInCore(s);

    u32bit  len    = S->sequenceLength();
    char    begseq = S->sequence()[0];
    bool    nnn    = V[begseq];
    u32bit  begpos = 0;
    u32bit  pos    = 0;

    for (pos=0; pos<len; pos++) {
      char seq = S->sequence()[pos];

      if (nnn != V[seq]) {
        fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                begseq, s, begpos, pos, pos - begpos);
        nnn = V[seq];
        begpos = pos;
        begseq = seq;
      }
    }

    fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            begseq, s, begpos, pos, pos - begpos);
    fprintf(stdout, ". "u32bitFMT" "u32bitFMT" "u32bitFMT"\n", s, pos, u32bitZERO);

    delete S;
  }

  delete F;
}


