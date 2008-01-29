#include "sim4polish.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void
s4p_deleteExon(sim4polish *p, int a) {
  char  *ed, *gd;
  int    i;
  int    editDistance    = 0;
  int    alignmentLength = 0;

  //  Warn if we don't have alignments -- this is now done by the
  //  driver (e.g., cleanPolishes.C)
  //
#if 0
  if ((p->exons[0].estAlignment == 0L) || (p->exons[0].genAlignment == 0L))
    fprintf(stderr, "s4p_deleteExon()-- Need alignments to recompute scores correctly!\n");
#endif

  //  Set the intron orientation for the exon before the one we are
  //  deleting:
  //    If we are deleting the first exon, there is no previous exon
  //    If we are deleting the last exon, set the previous to SIM4_INTRON_NONE
  //    Otherwise, set the previous to SIM4_INTRON_GAP
  //
  if (p->numExons > 1) {
    if (a == p->numExons-1)
      p->exons[a-1].intronOrientation = SIM4_INTRON_NONE;
    else if (a > 0)
      p->exons[a-1].intronOrientation = SIM4_INTRON_GAP;
  }

  //  Update the match scores
  //
  p->numMatches  -= p->exons[a].numMatches;
  p->numMatchesN -= p->exons[a].numMatchesN;

  //  Delete any alignment in the soon to be deleted exon
  //
  free(p->exons[a].genAlignment);
  free(p->exons[a].estAlignment);

  //  Shift all the exons down by one, and decrement the number of
  //  exons present in the list.
  //
  for (i=a+1; i<p->numExons; i++)
    memcpy(p->exons+i-1, p->exons+i, sizeof(sim4polishExon));

  p->numExons--;

  //  To be safe, zero out the space for the (still allocated, but
  //  unused) last exon
  //
  memset(p->exons+p->numExons, 0, sizeof(sim4polishExon));

  //  The strand orientation becomes unknown if we delete internal
  //  exons, or we end up with only one exon.
  //
  if (((0 < a) && (a < p->numExons)) ||
      (p->numExons == 1))
    p->strandOrientation = SIM4_STRAND_UNKNOWN;


  //  Compute the alignment length and the number of edits.
  //
  alignmentLength = 0;
  editDistance    = 0;

  p->numCovered   = 0;

  for (i=0; i<p->numExons; i++) {
    ed = p->exons[i].estAlignment;
    gd = p->exons[i].genAlignment;

    if (ed && gd) {
      alignmentLength += 2 * strlen(ed);
      for (; *ed && *gd; ed++, gd++) {
        if (*ed != *gd)
          editDistance++;
      }
    } else {
      int len = p->exons[i].estTo - p->exons[i].estFrom + 1 + p->exons[i].estTo - p->exons[i].estFrom + 1;

      alignmentLength += len;
      editDistance    += len / 2 - p->exons[i].numMatches - p->exons[i].numMatchesN;
    }

    p->numCovered += p->exons[i].genTo - p->exons[i].genFrom + 1;
  }

#if 0
  fprintf(stdout, "Found (new)alignLen = %d\n", alignmentLength);
  fprintf(stdout, "Found (new)editDist = %d\n", editDistance);
#endif

  //  Fix the scores for the match.  Special case; if there is only
  //  one exon left, the score for the exon is the score for the
  //  match.
  //
  if (p->numExons == 1)
    p->percentIdentity = p->exons[0].percentIdentity;
  else
    p->percentIdentity = s4p_percentIdentityApprox(editDistance, alignmentLength);

  //  Update the query sequence identity
  //
  p->querySeqIdentity = s4p_percentCoverageApprox(p);
}

