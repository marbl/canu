#include "sim4polish.H"
#include "memory.h"

#include <errno.h>
#include <string.h>


void
sim4polishExon::s4p_copyExon(sim4polishExon *orig) {

  if (orig == 0L)
    return;

  _estFrom            = orig->_estFrom;
  _estTo              = orig->_estTo;
  _genFrom            = orig->_genFrom;
  _genTo              = orig->_genTo;
  _numMatches         = orig->_numMatches;
  _numMatchesN        = orig->_numMatchesN;
  _percentIdentity    = orig->_percentIdentity;
  _intronOrientation  = orig->_intronOrientation;

  delete [] _estAlignment;
  delete [] _genAlignment;

  _estAlignment       = NULL;
  _genAlignment       = NULL;

  if (orig->_estAlignment) {
    u32bit len = strlen(orig->_estAlignment) + 1;
    _estAlignment = new char [len];
    memcpy(_estAlignment, orig->_estAlignment, sizeof(char) * len);
  }

  if (orig->_genAlignment) {
    u32bit len = strlen(orig->_genAlignment) + 1;
    _genAlignment = new char [len];
    memcpy(_genAlignment, orig->_genAlignment, sizeof(char) * len);
  }
}


void
sim4polish::s4p_copyPolish(sim4polish *orig, u32bit exonNum) {

  if (orig == 0L)
    return;

  _estID             = orig->_estID;
  _estLen            = orig->_estLen;
  _estPolyA          = orig->_estPolyA;
  _estPolyT          = orig->_estPolyT;

  _genID             = orig->_genID;
  _genRegionOffset   = orig->_genRegionOffset;
  _genRegionLength   = orig->_genRegionLength;

  _numMatches        = orig->_numMatches;
  _numMatchesN       = orig->_numMatchesN;
  _numCovered        = orig->_numCovered;
  _percentIdentity   = orig->_percentIdentity;
  _querySeqIdentity  = orig->_querySeqIdentity;
  _matchOrientation  = orig->_matchOrientation;
  _strandOrientation = orig->_strandOrientation;

  delete [] _comment;
  delete [] _estDefLine;
  delete [] _genDefLine;

  _comment           = NULL;
  _estDefLine        = NULL;
  _genDefLine        = NULL;

  delete [] _exons;
    
  _numExons          = 0;
  _exons             = NULL;

  //  Well, that was easy.  Onto the deep copy!

  if (orig->_comment) {
    u32bit len = strlen(orig->_comment) + 1;
    _comment = new char [len];
    memcpy(_comment, orig->_comment, sizeof(char) * len);
  }

  if (orig->_estDefLine) {
    u32bit len = strlen(orig->_estDefLine) + 1;
    _estDefLine = new char [len];
    memcpy(_estDefLine, orig->_estDefLine, sizeof(char) * len);
  }

  if (orig->_genDefLine) {
    u32bit len = strlen(orig->_genDefLine) + 1;
    _genDefLine = new char [len];
    memcpy(_genDefLine, orig->_genDefLine, sizeof(char) * len);
  }

  //  No exons?  We're done here.  Should never happen...

  if (orig->_numExons == 0)
    return;

  //  If told to copy one exon, just copy one exon....and then rebuild statistics.

  if (exonNum < orig->_numExons) {
    _numExons = 1;
    _exons    = new sim4polishExon [_numExons];

    _exons[0].s4p_copyExon(orig->_exons + exonNum);

    //  Rebuild stats
    _numMatches       = _exons[0]._numMatches;
    _numMatchesN      = _exons[0]._numMatchesN;
    _numCovered       = _exons[0]._estTo - _exons[0]._estFrom + 1;
    _percentIdentity  = _exons[0]._percentIdentity;
    _querySeqIdentity = s4p_percentCoverageApprox();

    return;
  }

  //  Otherwise, copy all exons into the new polish

  _numExons = orig->_numExons;
  _exons    = new sim4polishExon [_numExons];

  for (u32bit i=0; i<_numExons; i++)
    _exons[i].s4p_copyExon(orig->_exons + i);
}
