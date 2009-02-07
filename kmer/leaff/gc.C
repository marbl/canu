#include "bio++.H"
#include "seqCache.H"


void
computeGCcontent(char *filename) {
  seqCache   *A = new seqCache(filename);

  for (u32bit idx=0; idx < A->getNumberOfSequences(); idx++) {
    seqInCore *S = A->getSequenceInCore(idx);
    char      *s = S->sequence();
    u32bit     genomeLength = S->sequenceLength();

    fprintf(stdout, ">%s\n", S->header());

    int gc[256] = {0};
    gc['c'] = 1;
    gc['C'] = 1;
    gc['g'] = 1;
    gc['G'] = 1;

    //  Replace the sequence with "g or c".  We can't do this inline,
    //  since output reports the sequence too.  The extra 1000 at the
    //  end is important, since we do not bother checking for the end
    //  of the valid data, just assume that it's zero.
    //
    char                *g = new char [S->sequenceLength() + 1000];
    for (u32bit i=0; i<genomeLength+1000; i++)
      g[i] = 0;
    for (u32bit i=0; i<genomeLength; i++)
      g[i] = gc[s[i]];

    //  This stolen from depthOfPolishes.C

    u32bit  ave3    = 0;
    u32bit  ave5    = 0;
    u32bit  ave11   = 0;
    u32bit  ave51   = 0;
    u32bit  ave101  = 0;
    u32bit  ave201  = 0;
    u32bit  ave501  = 0;
    u32bit  ave1001 = 0;
    u32bit  ave2001 = 0;

    //  Preload the averages
    ave3   += g[0];
    ave5   += g[0] + g[1];

    for (u32bit i=0; i<5; i++)     ave11   += g[i];
    for (u32bit i=0; i<25; i++)    ave51   += g[i];
    for (u32bit i=0; i<50; i++)    ave101  += g[i];
    for (u32bit i=0; i<100; i++)   ave201  += g[i];
    for (u32bit i=0; i<250; i++)   ave501  += g[i];
    for (u32bit i=0; i<500; i++)   ave1001 += g[i];
    for (u32bit i=0; i<1000; i++)  ave2001 += g[i];

    for (u32bit i=0; i<genomeLength; i++) {
      ave3    += g[i+1]    - ((i >    1) ? g[i-2]    : 0);
      ave5    += g[i+2]    - ((i >    2) ? g[i-3]    : 0);
      ave11   += g[i+5]    - ((i >    5) ? g[i-6]    : 0);
      ave51   += g[i+25]   - ((i >   25) ? g[i-25]   : 0);
      ave101  += g[i+50]   - ((i >   50) ? g[i-51]   : 0);
      ave201  += g[i+100]  - ((i >  100) ? g[i-101]  : 0);
      ave501  += g[i+250]  - ((i >  250) ? g[i-251]  : 0);
      ave1001 += g[i+500]  - ((i >  500) ? g[i-501]  : 0);
      ave2001 += g[i+1000] - ((i > 1000) ? g[i-1001] : 0);

      fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
              i,
              s[i],
              ave3    / (double)((i >=   1)  ? 3    - ((i < genomeLength -   1) ? 0 : i +    2 - genomeLength) : i+2),
              ave5    / (double)((i >=   2)  ? 5    - ((i < genomeLength -   2) ? 0 : i +    3 - genomeLength) : i+3),
              ave11   / (double)((i >=   5)  ? 11   - ((i < genomeLength -   4) ? 0 : i +    5 - genomeLength) : i+6),
              ave51   / (double)((i >=  25)  ? 51   - ((i < genomeLength -  24) ? 0 : i +   25 - genomeLength) : i+26),
              ave101  / (double)((i >=  50)  ? 101  - ((i < genomeLength -  49) ? 0 : i +   50 - genomeLength) : i+51),
              ave201  / (double)((i >= 100)  ? 201  - ((i < genomeLength -  99) ? 0 : i +  100 - genomeLength) : i+101),
              ave501  / (double)((i >= 250)  ? 501  - ((i < genomeLength - 249) ? 0 : i +  250 - genomeLength) : i+251),
              ave1001 / (double)((i >= 500)  ? 1001 - ((i < genomeLength - 499) ? 0 : i +  500 - genomeLength) : i+501),
              ave2001 / (double)((i >= 1000) ? 2001 - ((i < genomeLength - 999) ? 0 : i + 1000 - genomeLength) : i+1001));
    }

    delete [] g;
    delete    S;
  }
}
  
