#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes polishes from stdin as a one-line-per-match format, space-based!

bool extendedFormat = false;

void
output(sim4polish *p,
       char       *Ep,
       char       *Gp,
       uint32      a,
       uint32      b,
       bool        isExon) {
  uint32  beg = p->_exons[a]._estFrom - 1;
  uint32  end = p->_exons[b]._estTo;

  if (p->_matchOrientation == SIM4_MATCH_COMPLEMENT) {
    beg = p->_estLen - beg;
    end = p->_estLen - end;
  }

  double  ident = p->_exons[a]._percentIdentity;
  double  cover = 0.0;

  //  If we're not a single exon, compute the real identity of the whole thing.
  //
  if (isExon == false) {
    if (p->_exons[a]._estAlignment) {
      ident = p->s4p_percentIdentityExact();
      cover = p->s4p_percentCoverageExact();
    } else {
      ident = p->_percentIdentity;
      cover = p->_querySeqIdentity;
    }
  }

  if (extendedFormat)
    fprintf(stdout, "%s\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%s\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%6.3f\t%6.3f\n",
            Ep, p->_estID,
            p->_estLen, a, beg, end,
            Gp, p->_genID,
            p->_exons[a]._genFrom - 1, p->_exons[b]._genTo,
            ident, cover);
  else
    fprintf(stdout, "%s\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%s\t"uint32FMT"\t"uint32FMT"\t%6.3f\t%6.3f\n",
            Ep, p->_estLen, a, beg, end,
            Gp, p->_exons[a]._genFrom - 1, p->_exons[b]._genTo,
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
    if        (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-fullquery") == 0) {
      wholeEDefLine = true;
    } else if (strcmp(argv[arg], "-fullgenomic") == 0) {
      wholeGDefLine = true;
    } else if (strcmp(argv[arg], "-exons") == 0) {
      doExons = true;
    } else if (strcmp(argv[arg], "-extended") == 0) {
      extendedFormat = true;
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (isatty(fileno(stdin)) || (err)) {
    fprintf(stderr, "usage: %s [options] < IN > OUT\n", argv[0]);
    fprintf(stderr, "  -v             be chatty\n");
    fprintf(stderr, "  -fullquery     output the whole query def line\n");
    fprintf(stderr, "  -fullgenomic   output the whole genomic def line\n");
    fprintf(stderr, "  -exons         include exons\n");
    fprintf(stderr, "  -extended      include the IDX of each sequence\n");
    exit(1);
  }

  if (extendedFormat)
    fprintf(stdout, "cDNAid\tcDNAidx\tcDNAlen\texonNum\tbegin\tend\tgenomicid\tgenomicidx\tbegin\tend\tidentity\tcoverage\n");
  else
    fprintf(stdout, "cDNAid\tcDNAlen\texonNum\tbegin\tend\tgenomicid\tbegin\tend\tidentity\tcoverage\n");

  char          E[1024], *Ep;
  char          G[1024], *Gp;
  splitToWords  W;

  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *p = 0L;

  while (R->nextAlignment(p)) {
    if (wholeEDefLine == true) {
      Ep = p->_estDefLine;
    } else {
      W.split(p->_estDefLine);
      strcpy(E, W[0] + ((W[0][0] == '>') ? 1 : 0));
      Ep = E;
    }

    if (wholeGDefLine == true) {
      Gp = p->_genDefLine;
    } else {
      W.split(p->_genDefLine);
      strcpy(G, W[0] + ((W[0][0] == '>') ? 1 : 0));
      Gp = G;
    }

    if (doExons == false) {
      output(p, Ep, Gp, 0, p->_numExons-1, false);
    } else {
      for (uint32 i=0; i<p->_numExons; i++)
        output(p, Ep, Gp, i, i, true);
    }
  }

  return(0);
}

