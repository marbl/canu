#ifndef ENCODEDQUERY_H
#define ENCODEDQUERY_H

#include <stdio.h>
#include <stdlib.h>
#include "bio++.H"
#include "seqCache.H"

class encodedQuery {
public:
  encodedQuery(seqInCore  *S,
               uint32      k);
  ~encodedQuery();

  void             test(seqInCore  *S);

  uint32           IID(void)                 { return(_iid); };

  uint32           numberOfMers(void)        { return(_mersTotal); };
  uint32           numberOfValidMers(void)   { return(_mersAvail); };

  //  Return either an approximation or the exact number of bp covered by valid mers.
  //
  uint32           bpCovered(bool doCompute);
  uint32           bpTotal(void);

  uint64           getMer(uint32 i, bool isReverse);
  bool             getSkip(uint32 i, bool isReverse);
  void             setSkip(uint32 i, bool isReverse);

  void             addOutput(void *output, uint32 size);

  char            *theOutput(void)       { return(_output); };
  uint32           theOutputLength(void) { return(_outputLen); };
  uint32           numberOfResults(void) { return(_numberOfResults); };
private:
  uint32          _iid;
  uint32          _sequenceLength;

  uint32          _merSize;

  uint32          _mersTotal;  //  Number of mers possible in the query
  uint32          _mersAvail;  //  Number of mers not masked out
  uint64         *_mers;       //  List of mers
  uint8          *_skip;

  uint32          _numberOfResults;
  char           *_output;
  uint32          _outputLen;
  uint32          _outputMax;
};


inline
uint32
encodedQuery::bpCovered(bool doCompute) {
  uint32  bp = numberOfValidMers();

  if (doCompute) {
    merCovering   *IL = new merCovering(_merSize);

    for (uint32 qi=0; qi<numberOfMers(); qi++) {
      if (getSkip(qi, false) == false)
        IL->addMer(qi);
    }

    bp = IL->sumOfLengths();
    delete IL;
  }

  return(bp);
}

inline
uint32
encodedQuery::bpTotal(void) {
  return(_sequenceLength);
}


//  XXX:  We need to extend get*() to also take a isReverse flag,
//  and to then return the reverse-complement mer, skip, etc.

inline
uint64
encodedQuery::getMer(uint32 i, bool isReverse) {

  if (isReverse) {
    i = _mersTotal - i - 1;
    return(reverseComplementMer(_merSize, _mers[i]));
  } else {
    return(_mers[i]);
  }
}

inline
bool
encodedQuery::getSkip(uint32 i, bool isReverse) {

  if (isReverse)
    i = _mersTotal - i - 1;

  return(_skip[i]);
}

inline
void
encodedQuery::setSkip(uint32 i, bool isReverse) {

  if (isReverse)
    i = _mersTotal - i - 1;

  //  If skip[i] is already set, don't adjust; otherwise, subtract
  //  one from the available.
  _mersAvail -= 1 - _skip[i];
  _skip[i]    = 1;
}


#endif  //  ENCODEDQUERY_H
