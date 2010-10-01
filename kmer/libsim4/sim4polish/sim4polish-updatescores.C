#include <math.h>
#include "sim4polish.H"


void
sim4polish::s4p_updateAlignmentScores(void) {
  u32bit  ni = 0, numInDel    = 0;
  u32bit  ne = 0, numEdits    = 0;
  u32bit  nn = 0, numMatchesN = 0;
  u32bit  nm = 0, numMatches  = 0;
  u32bit  al = 0, alignmentLength = 0;
  u32bit  nc = 0, numCovered  = 0;

  u32bit  estn = 0;
  u32bit  genn = 0;

  for (u32bit exon=0; exon<_numExons; exon++) {
    char *est = _exons[exon]._estAlignment;
    char *gen = _exons[exon]._genAlignment;

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
          *est = toupper(*est);
          *gen = toupper(*gen);
        } else if (estn && genn) {
          //  Both are N.  It isn't a match and it isn't an edit.
          //
          nn++;
          *est = toupper(*est);
          *gen = toupper(*gen);
        } else if (estn || genn) {
          //  One is an N.  Someone has low quality sequence, and we
          //  should penalize.  We need to special case this because
          //  IUPACidentity[][] claims N matches all.
          //
          ne++;
          *est = toupper(*est);
          *gen = toupper(*gen);
        } else if (IUPACidentity[(int)*est][(int)*gen]) {
          //  Got a match.
          nm++;
          *est = tolower(*est);
          *gen = tolower(*gen);
        } else {
          //  Got a substitution
          ne++;
          *est = toupper(*est);
          *gen = toupper(*gen);
        }

        est++;
        gen++;
      }
    }

    _exons[exon]._numMatches  = nm;
    _exons[exon]._numMatchesN = nn;

    al = (_exons[exon]._genTo - _exons[exon]._genFrom + 1 +
          _exons[exon]._estTo - _exons[exon]._estFrom + 1 +
          ne);
    nc = (_exons[exon]._estTo - _exons[exon]._estFrom + 1);

    _exons[exon]._percentIdentity = s4p_percentIdentityApprox(ne, al);

    numInDel        += ni;
    numEdits        += ne;
    numMatchesN     += nn;
    numMatches      += nm;
    alignmentLength += al;
    numCovered      += nc;
  }

  _numMatches  = numMatches;
  _numMatchesN = numMatchesN;
  _numCovered  = numCovered;

#if 0
  fprintf(stderr, "numInDel    = %d\n", numInDel);
  fprintf(stderr, "numEdits    = %d\n", numEdits);
  fprintf(stderr, "numMatchesN = %d\n", numMatchesN);
  fprintf(stderr, "numMatches  = %d\n", numMatches);
  fprintf(stderr, "alignLen    = %d\n", alignmentLength);
  fprintf(stderr, "numCovered  = %d\n", numCovered);
#endif

  _percentIdentity  = s4p_percentIdentityApprox(numEdits, alignmentLength);
  _querySeqIdentity = s4p_percentCoverageApprox();
}


int
sim4polish::s4p_percentCoverageApprox(void) {
  int ret;

  if (_numCovered == _estLen - _estPolyA - _estPolyT)
    return(100); 

  return(((ret=(int)round(100.0 * _numCovered / (double)(_estLen - _estPolyA - _estPolyT))) < 100) ? ret : 99);
}


int
sim4polish::s4p_percentIdentityApprox(int numEdits, int alignmentLength) {
  int ret;

  if (alignmentLength == 0)
    return(0);

  if (numEdits == 0)
    return(100);

  return(((ret=(int)round(100.0 * (1 - 2.0 * numEdits / alignmentLength))) < 100) ? ret : 99);
}


double
sim4polish::s4p_percentCoverageExact(void) {
  return( 100 * (double)(_numCovered) / (double)(_estLen - _estPolyA - _estPolyT) );
}


double
sim4polish::s4p_percentIdentityExact(void) {
  u32bit  ni = 0, numInDel    = 0;
  u32bit  ne = 0, numEdits    = 0;
  u32bit  nn = 0, numMatchesN = 0;
  u32bit  nm = 0, numMatches  = 0;
  u32bit  al = 0, alignmentLength = 0;
  u32bit  nc = 0, numCovered  = 0;

  u32bit  estn = 0;
  u32bit  genn = 0;

  double ret = 0.0;

  for (u32bit exon=0; exon<_numExons; exon++) {
    char *est = _exons[exon]._estAlignment;
    char *gen = _exons[exon]._genAlignment;

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
    _exons[exon]._numMatches  = nm;
    _exons[exon]._numMatchesN = nn;
#endif

    al = (_exons[exon]._genTo - _exons[exon]._genFrom + 1 +
          _exons[exon]._estTo - _exons[exon]._estFrom + 1 +
          ne);
    nc = (_exons[exon]._genTo - _exons[exon]._genFrom + 1);

#if 0
    _exons[exon]._percentIdentity = s4p_percentIdentityApprox(ne, al);
#endif

    numInDel        += ni;
    numEdits        += ne;
    numMatchesN     += nn;
    numMatches      += nm;
    alignmentLength += al;
    numCovered      += nc;
  }

#if 0
  _numMatches  = numMatches;
  _numMatchesN = numMatchesN;
  _numCovered  = numCovered;
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
