#include "sim4polish.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char*
s4p_polishToString(sim4polish *p) {
  u32bit  spaceNeeded = 0;
  u32bit  i;
  char   *outc=NULL, *outs=NULL;
  char   *mOri="", *sOri="", *iOri="";

  if (p == 0L)
    return(0L);

  //  Make a decent estimate of how much space we'll need to
  //  store the string
  //
  spaceNeeded = 1024 + 128 * p->numExons;

  if (p->comment)     spaceNeeded += strlen(p->comment);
  if (p->estDefLine)  spaceNeeded += strlen(p->estDefLine);
  if (p->genDefLine)  spaceNeeded += strlen(p->genDefLine);

  for (i=0; i<p->numExons; i++) {
    if (p->exons[i].estAlignment)
      spaceNeeded += 2 * strlen(p->exons[i].estAlignment);
  }

  outc = outs = (char *)malloc(spaceNeeded * sizeof(char));

  switch (p->matchOrientation) {
    case SIM4_MATCH_FORWARD:
      mOri = "forward";
      break;
    case SIM4_MATCH_COMPLEMENT:
      mOri = "complement";
      break;
    case SIM4_MATCH_ERROR:
      mOri = "error";
      break;
    default:
      fprintf(stderr, "sim4reader: Unknown matchOrientation '"u32bitFMT"' in printPolish()\n", p->matchOrientation);
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
    case SIM4_STRAND_ERROR:
      sOri = "error";
      break;
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"u32bitFMT"' in printPolish()\n", p->matchOrientation);
      sOri = "UNKNOWN";
      break;
  }

  sprintf(outc, "sim4begin\n"u32bitFMT"["u32bitFMT"-"u32bitFMT"-"u32bitFMT"] "u32bitFMT"["u32bitFMT"-"u32bitFMT"] <"u32bitFMT"-"u32bitFMT"-"u32bitFMT"-%s-%s>\n",
          p->estID, p->estLen, p->estPolyA, p->estPolyT,
          p->genID, p->genLo, p->genHi,
          p->numMatches, p->numMatchesN, p->percentIdentity, mOri, sOri);
  while (*outc)  outc++;

  if (p->comment) {
    sprintf(outc, "comment=%s\n", p->comment);
    while (*outc)  outc++;
  }

  if (p->estDefLine) {
    sprintf(outc, "edef=%s\n", p->estDefLine);
    while (*outc)  outc++;
  }

  if (p->genDefLine) {
    sprintf(outc, "ddef=%s\n", p->genDefLine);
    while (*outc)  outc++;
  }

  for (i=0; i<p->numExons; i++) {
    switch (p->exons[i].intronOrientation) {
      case SIM4_INTRON_POSITIVE:
        iOri = " ->";
        break;
      case SIM4_INTRON_NEGATIVE:
        iOri = " <-";
        break;
      case SIM4_INTRON_AMBIGUOUS:
        iOri = " --";
        break;
      case SIM4_INTRON_GAP:
        iOri = " ==";
        break;
      case SIM4_INTRON_ERROR:
        iOri = " ??";
        break;
      default:
        iOri = "";
        break;
    }

    sprintf(outc, ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">%s\n",
            p->exons[i].estFrom, p->exons[i].estTo,
            p->exons[i].genFrom, p->exons[i].genTo,
            p->exons[i].numMatches, p->exons[i].numMatchesN, p->exons[i].percentIdentity, iOri);

    while (*outc)  outc++;
  }

  //fputs(outs, stderr);

  for (i=0; i<p->numExons; i++) {
    if (p->exons[i].estAlignment) {
      strcpy(outc, p->exons[i].estAlignment);
      while (*outc)  outc++;
      *outc = '\n';
      outc++;
    }
    if (p->exons[i].genAlignment) {
      strcpy(outc, p->exons[i].genAlignment);
      while (*outc)  outc++;
      *outc = '\n';
      outc++;
    }
  }

  strcpy(outc, "sim4end\n");

  return(outs);
}

