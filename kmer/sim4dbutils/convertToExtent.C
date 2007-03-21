#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes polishes from stdin as a one-line-per-match format, space-based!


void
output(sim4polish *p,
       char       *Ep,
       char       *Gp,
       u32bit      a,
       u32bit      b,
       bool        isExon) {
  u32bit  beg = p->exons[a].estFrom - 1;
  u32bit  end = p->exons[b].estTo;

  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
    beg = p->estLen - beg;
    end = p->estLen - end;
  }

  double  ident = p->exons[a].percentIdentity;
  double  cover = 0.0;

  //  If we're not a single exon, compute the real identity of the whole thing.
  //
  if (isExon == false) {
    if (p->exons[a].estAlignment) {
      s4p_percentIdentityExact(p);
      s4p_percentCoverageExact(p);
    } else {
      ident = p->percentIdentity;
      cover = p->querySeqIdentity;
    }
  }

  fprintf(stdout, "%s\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%s\t"u32bitFMT"\t"u32bitFMT"\t%6.3f\t%6.3f\n",
          Ep, p->estLen, a, beg, end,
          Gp, p->exons[a].genFrom - 1, p->exons[b].genTo,
          ident, cover);
}


int
main(int argc, char **argv) {
  bool beVerbose     = false;
  bool wholeEDefLine = false;
  bool wholeGDefLine = false;
  bool doExons       = false;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-q", 2) == 0) {
      wholeEDefLine = true;
    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      wholeGDefLine = true;
    } else if (strncmp(argv[arg], "-exons", 2) == 0) {
      doExons = true;
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (isatty(fileno(stdin)) || (err)) {
    fprintf(stderr, "usage: %s [-v] [-q] [-g] [-exons]\n", argv[0]);
    exit(1);
  }

  fprintf(stdout, "cDNAid\tcDNAlen\texonNum\tbegin\tend\tgenomicid\tbegin\tend\tidentity\tcoverage\n");

  char          E[1024], *Ep;
  char          G[1024], *Gp;
  splitToWords  W;

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if (p != 0L) {

      if (wholeEDefLine == true) {
        Ep = p->estDefLine;
      } else {
        W.split(p->estDefLine);
        strcpy(E, W[0] + ((W[0][0] == '>') ? 1 : 0));
        Ep = E;
      }

      if (wholeGDefLine == true) {
        Gp = p->genDefLine;
      } else {
        W.split(p->genDefLine);
        strcpy(G, W[0] + ((W[0][0] == '>') ? 1 : 0));
        Gp = G;
      }

      if (doExons == false) {
        output(p, Ep, Gp, 0, p->numExons-1, false);
      } else {
        for (u32bit i=0; i<p->numExons; i++)
          output(p, Ep, Gp, i, i, true);
      }

      s4p_destroyPolish(p);
    }
  }

  return(0);
}

