#include "sim4polish.h"

void
s4p_printPolish(FILE *O, sim4polish *p) {
  int          i;
  char const  *mOri;
  char const  *sOri;

  if (p) {
    fprintf(O, "sim4begin\n");

    switch (p->matchOrientation) {
      case SIM4_MATCH_FORWARD:
        mOri = "forward";
        break;
      case SIM4_MATCH_COMPLEMENT:
        mOri = "complement";
        break;
      default:
        fprintf(stderr, "sim4reader: Unknown matchOrientation '%d' in printPolish()\n", p->matchOrientation);
        mOri = "UNKNOWN";
        break;
    }

    switch (p->strandOrientation) {
      case SIM4_STRAND_POSITIVE:
        sOri = "forward";
        break;
      case SIM4_STRAND_NEGATIVE:
        sOri = "reverse";
        break;
      case SIM4_STRAND_UNKNOWN:
        sOri = "unknown";
        break;
      case SIM4_STRAND_INTRACTABLE:
        sOri = "intractable";
        break;
      case SIM4_STRAND_FAILED:
        sOri = "aborted";
        break;
      default:
        fprintf(stderr, "sim4reader: Unknown strandOrientation '%d' in printPolish()\n", p->matchOrientation);
        sOri = "UNKNOWN";
        break;
    }

    fprintf(O, "%d[%d-%d-%d] %d[%d-%d] <%d-%d-%d-%s-%s>\n",
            p->estID, p->estLen, p->estPolyA, p->estPolyT,
            p->genID, p->genLo, p->genHi,
            p->numMatches, p->numMatchesN, p->percentIdentity, mOri, sOri);

    if (p->comment)
      fprintf(O, "comment=%s\n", p->comment);
    if (p->estDefLine)
      fprintf(O, "edef=%s\n", p->estDefLine);
    if (p->genDefLine)
      fprintf(O, "ddef=%s\n", p->genDefLine);
    for (i=0; i<p->numExons; i++) {
      fprintf(O, "%d-%d (%d-%d) <%d-%d-%d>",
              p->exons[i].estFrom, p->exons[i].estTo,
              p->exons[i].genFrom, p->exons[i].genTo,
              p->exons[i].numMatches, p->exons[i].numMatchesN, p->exons[i].percentIdentity);

      switch (p->exons[i].intronOrientation) {
        case SIM4_INTRON_POSITIVE:
          fprintf(O, " ->\n");
          break;
        case SIM4_INTRON_NEGATIVE:
          fprintf(O, " <-\n");
          break;
        case SIM4_INTRON_AMBIGUOUS:
          fprintf(O, " --\n");
          break;
        case SIM4_INTRON_GAP:
          fprintf(O, " ==\n");
          break;
        case SIM4_INTRON_NONE:
          fprintf(O, "\n");
          break;
        default:
          fprintf(O, "\n");
          break;
      }
    }
    for (i=0; i<p->numExons; i++) {
      if (p->exons[i].estAlignment)
        fprintf(O, "%s\n", p->exons[i].estAlignment);
      if (p->exons[i].genAlignment)
        fprintf(O, "%s\n", p->exons[i].genAlignment);
    }
    fprintf(O, "sim4end\n");
  }
}
