#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"


int
main(int argc, char **argv) {
  bool beVerbose = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
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

      IL.add(beg, end-beg);

      s4p_destroyPolish(p);
    }
  }

  intervalDepth ID(IL);

  for (u32bit i=0; i<ID.numberOfIntervals(); i++) {
    if (ID.de(i) > 0)
      fprintf(stdout, u64bitFMT"\t"u64bitFMT"\t"u32bitFMT"\n", ID.lo(i), ID.hi(i), ID.de(i));
  }

  return(0);
}

