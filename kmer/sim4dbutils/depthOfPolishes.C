#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

// ./genomics/sim4dbutils/depthOfPolishes -v < runA.1.ms12.filtered.sim4db > depth-out
// plot [112000:113000][] "depth-out" using 2 with lines

int
main(int argc, char **argv) {
  bool     beVerbose    = false;
  u32bit   genomeLength = 0;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-l", 2) == 0) {
      genomeLength = strtou32bit(argv[arg], 0L);
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  intervalList  IL;

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if (p != 0L) {
      u32bit  beg = p->genLo + p->exons[0].genFrom - 1;
      u32bit  end = p->genLo + p->exons[p->numExons-1].genTo;

      if (end > genomeLength)
        genomeLength = end;

      IL.add(beg, end-beg);

      s4p_destroyPolish(p);
    }
  }

  intervalDepth ID(IL);

  //  The extra 1000 here is so we can be lazy in the
  //  output section when computing averages.
  //
  u32bit       *DD = new u32bit [genomeLength + 1000];
  for (u32bit i=0; i<genomeLength + 1000; i++)
    DD[i] = 0;

  for (u32bit i=0; i<ID.numberOfIntervals(); i++) {
    u32bit l = ID.lo(i);
    u32bit h = ID.hi(i);
    u32bit d = ID.de(i);

    while (l < h) {
      DD[l] = d;
      l++;
    }
  }

  //  This stolen to leaff.C for %GC computation

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
  ave3   += DD[0];
  ave5   += DD[0] + DD[1];

  for (u32bit i=0; i<5; i++)     ave11   += DD[i];
  for (u32bit i=0; i<25; i++)    ave51   += DD[i];
  for (u32bit i=0; i<50; i++)    ave101  += DD[i];
  for (u32bit i=0; i<100; i++)   ave201  += DD[i];
  for (u32bit i=0; i<250; i++)   ave501  += DD[i];
  for (u32bit i=0; i<500; i++)   ave1001 += DD[i];
  for (u32bit i=0; i<1000; i++)  ave2001 += DD[i];

  for (u32bit i=0; i<genomeLength; i++) {
    ave3    += DD[i+1]    - ((i >    1) ? DD[i-2]    : 0);
    ave5    += DD[i+2]    - ((i >    2) ? DD[i-3]    : 0);
    ave11   += DD[i+5]    - ((i >    5) ? DD[i-6]    : 0);
    ave51   += DD[i+25]   - ((i >   25) ? DD[i-25]   : 0);
    ave101  += DD[i+50]   - ((i >   50) ? DD[i-51]   : 0);
    ave201  += DD[i+100]  - ((i >  100) ? DD[i-101]  : 0);
    ave501  += DD[i+250]  - ((i >  250) ? DD[i-251]  : 0);
    ave1001 += DD[i+500]  - ((i >  500) ? DD[i-501]  : 0);
    ave2001 += DD[i+1000] - ((i > 1000) ? DD[i-1001] : 0);

    fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
            i,
            DD[i],
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

  return(0);
}
