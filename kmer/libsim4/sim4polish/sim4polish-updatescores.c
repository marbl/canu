#include <math.h>
#include "sim4polish.h"

void
s4p_updateAlignmentScores(sim4polish *p) {

  int  exon;

  int  ni = 0, numInDel    = 0;
  int  ne = 0, numEdits    = 0;
  int  nn = 0, numMatchesN = 0;
  int  nm = 0, numMatches  = 0;
  int  al = 0, alignmentLength = 0;
  int  nc = 0, numCovered  = 0;

  int  estn = 0;
  int  genn = 0;

  for (exon=0; exon<p->numExons; exon++) {
    char *est = p->exons[exon].estAlignment;
    char *gen = p->exons[exon].genAlignment;

    al = 0;

    ni = 0;
    ne = 0;
    nn = 0;
    nm = 0;

    if (est && gen) {
      while (*est && *gen) {
        estn = (*est == 'N') || (*est == 'n');
        genn = (*gen == 'N') || (*gen == 'n');

        if        ((*est == '-') || (*gen == '-')) {
          ni++;
          ne++;
        } else if (estn && genn) {
          //  Both are N.  It isn't a match and it isn't an edit.
          //
          nn++;
        } else if (estn || genn) {
          //  One is an N.  Someone has low quality sequence, and we
          //  should penalize.  We need to special case this because
          //  IUPACidentity[][] claims N matches all.
          //
          ne++;
        } else if (IUPACidentity[(int)*est][(int)*gen]) {
          //  Got a match.
          nm++;
        } else {
          //  Got a substitution
          ne++;
        }

        est++;
        gen++;
      }
    }

    p->exons[exon].numMatches  = nm;
    p->exons[exon].numMatchesN = nn;

    al = (p->exons[exon].genTo - p->exons[exon].genFrom + 1 +
          p->exons[exon].estTo - p->exons[exon].estFrom + 1 +
          ne);
    nc = (p->exons[exon].genTo - p->exons[exon].genFrom + 1);

    p->exons[exon].percentIdentity = (int)floor(100*(1 - 2.0 * ne / (double)(al)));

    numInDel        += ni;
    numEdits        += ne;
    numMatchesN     += nn;
    numMatches      += nm;
    alignmentLength += al;
    numCovered      += nc;
  }

  p->numMatches  = numMatches;
  p->numMatchesN = numMatchesN;
  p->numCovered  = numCovered;

#if 0
  fprintf(stderr, "numInDel    = %d\n", numInDel);
  fprintf(stderr, "numEdits    = %d\n", numEdits);
  fprintf(stderr, "numMatchesN = %d\n", numMatchesN);
  fprintf(stderr, "numMatches  = %d\n", numMatches);
  fprintf(stderr, "alignLen    = %d\n", alignmentLength);
  fprintf(stderr, "numCovered  = %d\n", numCovered);
#endif

  p->percentIdentity  = 0;
  if (alignmentLength > 0)
    p->percentIdentity  = (int)floor(100*(1 - 2.0 * numEdits / (double)(alignmentLength)));

  p->querySeqIdentity = (int)floor(100 * (double)(p->numCovered) / (double)(p->estLen - p->estPolyA - p->estPolyT));
}



double
s4p_percentCoverage(sim4polish *p) {
  return( 100 * (double)(p->numCovered) / (double)(p->estLen - p->estPolyA - p->estPolyT) );
}



//  A very expensive and accurate calculation of the percent identity.
double
s4p_percentIdentity(sim4polish *p) {

  int  exon;

  int  ni = 0, numInDel    = 0;
  int  ne = 0, numEdits    = 0;
  int  nn = 0, numMatchesN = 0;
  int  nm = 0, numMatches  = 0;
  int  al = 0, alignmentLength = 0;
  int  nc = 0, numCovered  = 0;

  int  estn = 0;
  int  genn = 0;

  double ret = 0.0;

  for (exon=0; exon<p->numExons; exon++) {
    char *est = p->exons[exon].estAlignment;
    char *gen = p->exons[exon].genAlignment;

    al = 0;

    ni = 0;
    ne = 0;
    nn = 0;
    nm = 0;

    if (est && gen) {
      while (*est && *gen) {
        estn = (*est == 'N') || (*est == 'n');
        genn = (*gen == 'N') || (*gen == 'n');

        if        ((*est == '-') || (*gen == '-')) {
          ni++;
          ne++;
        } else if (estn && genn) {
          //  Both are N.  It isn't a match and it isn't an edit.
          //
          nn++;
        } else if (estn || genn) {
          //  One is an N.  Someone has low quality sequence, and we
          //  should penalize.  We need to special case this because
          //  IUPACidentity[][] claims N matches all.
          //
          ne++;
        } else if (IUPACidentity[(int)*est][(int)*gen]) {
          //  Got a match.
          nm++;
        } else {
          //  Got a substitution
          ne++;
        }

        est++;
        gen++;
      }
    }

#if 0
    p->exons[exon].numMatches  = nm;
    p->exons[exon].numMatchesN = nn;
#endif

    al = (p->exons[exon].genTo - p->exons[exon].genFrom + 1 +
          p->exons[exon].estTo - p->exons[exon].estFrom + 1 +
          ne);
    nc = (p->exons[exon].genTo - p->exons[exon].genFrom + 1);

#if 0
    p->exons[exon].percentIdentity = (int)floor(100*(1 - 2.0 * ne / (double)(al)));
#endif

    numInDel        += ni;
    numEdits        += ne;
    numMatchesN     += nn;
    numMatches      += nm;
    alignmentLength += al;
    numCovered      += nc;
  }

#if 0
  p->numMatches  = numMatches;
  p->numMatchesN = numMatchesN;
  p->numCovered  = numCovered;
#endif

#if 0
  fprintf(stderr, "numInDel    = %d\n", numInDel);
  fprintf(stderr, "numEdits    = %d\n", numEdits);
  fprintf(stderr, "numMatchesN = %d\n", numMatchesN);
  fprintf(stderr, "numMatches  = %d\n", numMatches);
  fprintf(stderr, "alignLen    = %d\n", alignmentLength);
  fprintf(stderr, "numCovered  = %d\n", numCovered);
#endif

  if (alignmentLength > 0)
    ret = 100.0 * (1 - 2.0 * numEdits / (double)(alignmentLength));

  return(ret);
}
