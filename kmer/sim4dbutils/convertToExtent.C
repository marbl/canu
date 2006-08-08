#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes polishes from stdin as a one-line-per-match format, space-based!

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

  char          E[1024];
  char          G[1024];
  splitToWords  W;

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if (p != 0L) {

      W.split(p->estDefLine);
      strcpy(E, W[0] + ((W[0][0] == '>') ? 1 : 0));

      W.split(p->genDefLine);
      strcpy(G, W[0] + ((W[0][0] == '>') ? 1 : 0));

      u32bit  beg = p->exons[0].estFrom - 1;
      u32bit  end = p->exons[p->numExons-1].estTo;
      if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
        beg = p->estLen - beg;
        end = p->estLen - end;
      }

      if (p->exons[0].estAlignment)
        fprintf(stdout, "%s\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%s\t"u32bitFMT"\t"u32bitFMT"\t%6.3f\t%6.3f\n",
                E, p->estLen, beg, end,
                G, p->genLo + p->exons[0].genFrom - 1, p->genLo + p->exons[p->numExons-1].genTo,
                s4p_percentIdentityExact(p),
                s4p_percentCoverageExact(p));
      else
        fprintf(stdout, "%s\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%s\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                E, p->estLen, beg, end,
                G, p->genLo + p->exons[0].genFrom - 1, p->genLo + p->exons[p->numExons-1].genTo,
                p->percentIdentity,
                p->querySeqIdentity);


      s4p_destroyPolish(p);
    }
  }

  return(0);
}

