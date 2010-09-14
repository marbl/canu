#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void
sim4polish::s4p_deleteExon(u32bit a) {
  char  *ed, *gd;
  int    editDistance    = 0;
  int    alignmentLength = 0;

  //  Warn if we don't have alignments -- this is now done by the
  //  driver (e.g., cleanPolishes.C)
  //
#if 0
  if ((p->exons[0]._estAlignment == 0L) || (p->exons[0]._genAlignment == 0L))
    fprintf(stderr, "s4p_deleteExon()-- Need alignments to recompute scores correctly!\n");
#endif

  //  Set the intron orientation for the exon before the one we are
  //  deleting:
  //    If we are deleting the first exon, there is no previous exon
  //    If we are deleting the last exon, set the previous to SIM4_INTRON_NONE
  //    Otherwise, set the previous to SIM4_INTRON_GAP
  //
  if (_numExons > 1) {
    if (a == _numExons - 1)
      _exons[a-1]._intronOrientation = SIM4_INTRON_NONE;
    else if (a > 0)
      _exons[a-1]._intronOrientation = SIM4_INTRON_GAP;
  }

  //  Update the match scores
  //
  _numMatches  -= _exons[a]._numMatches;
  _numMatchesN -= _exons[a]._numMatchesN;

  //  Erase the exon we're removing, but save a copy so we can stash it in the
  //  soon-to-be-emptied last location.
  //
  _exons[a].s4p_clearExon();

  sim4polishExon d = _exons[a];

  //  Shift all the exons down by one, and decrement the number of
  //  exons present in the list.
  //
  for (u32bit i=a+1; i<_numExons; i++)
    _exons[i-1] = _exons[i];

  _numExons--;

  //  Stash the now deleted exon in the last spot, just to clear out the old contents.
  //
  _exons[_numExons] = d;

  //  The strand orientation becomes unknown if we delete internal
  //  exons, or we end up with only one exon.
  //
  if (((0 < a) && (a < _numExons)) ||
      (_numExons == 1))
    _strandOrientation = SIM4_STRAND_UNKNOWN;


  //  Compute the alignment length and the number of edits.
  //
  alignmentLength = 0;
  editDistance    = 0;

  _numCovered   = 0;

  for (u32bit i=0; i<_numExons; i++) {
    ed = _exons[i]._estAlignment;
    gd = _exons[i]._genAlignment;

    if (ed && gd) {
      alignmentLength += 2 * strlen(ed);
      for (; *ed && *gd; ed++, gd++) {
        if (*ed != *gd)
          editDistance++;
      }
    } else {
      int len = _exons[i]._estTo - _exons[i]._estFrom + 1 + _exons[i]._estTo - _exons[i]._estFrom + 1;

      alignmentLength += len;
      editDistance    += len / 2 - _exons[i]._numMatches - _exons[i]._numMatchesN;
    }

    _numCovered += _exons[i]._genTo - _exons[i]._genFrom + 1;
  }

#if 0
  fprintf(stdout, "Found (new)alignLen = %d\n", alignmentLength);
  fprintf(stdout, "Found (new)editDist = %d\n", editDistance);
#endif

  //  Fix the scores for the match.  Special case; if there is only
  //  one exon left, the score for the exon is the score for the
  //  match.
  //
  if (_numExons == 1)
    _percentIdentity = _exons[0]._percentIdentity;
  else
    _percentIdentity = s4p_percentIdentityApprox(editDistance, alignmentLength);

  //  Update the query sequence identity
  //
  _querySeqIdentity = s4p_percentCoverageApprox();
}

