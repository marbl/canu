#ifndef ENCODEDQUERY_H
#define ENCODEDQUERY_H

#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include "bio++.H"


class encodedQuery {
public:
  encodedQuery(FastASequenceInCore  *S,
               u32bit                k);
  ~encodedQuery();

  void             test(FastASequenceInCore  *S);

  u32bit           IID(void)                 { return(_iid); };

  u32bit           numberOfMers(void)        { return(_mersTotal); };
  u32bit           numberOfValidMers(void)   { return(_mersAvail); };

  //  Return either an approximation or the exact number of bp covered by valid mers.
  //
  u32bit           bpCovered(bool doCompute);
  u32bit           bpTotal(void);

  u64bit           getMer(u32bit i, bool isReverse);
  bool             getSkip(u32bit i, bool isReverse);
  void             setSkip(u32bit i, bool isReverse);

  void             addOutput(void *output, u32bit size);

  char            *theOutput(void)       { return(_output); };
  u32bit           theOutputLength(void) { return(_outputLen); };
  u32bit           numberOfResults(void) { return(_numberOfResults); };
private:
  u32bit          _iid;
  u32bit          _sequenceLength;

  u32bit          _merSize;

  u32bit          _mersTotal;  //  Number of mers possible in the query
  u32bit          _mersAvail;  //  Number of mers not masked out
  u64bit         *_mers;       //  List of mers
  u8bit          *_skip;

  u32bit          _numberOfResults;
  char           *_output;
  u32bit          _outputLen;
  u32bit          _outputMax;
};


inline
u32bit
encodedQuery::bpCovered(bool doCompute) {
  u32bit  bp = numberOfValidMers();

  if (doCompute) {
    merCovering   *IL = new merCovering(_merSize);

    for (u32bit qi=0; qi<numberOfMers(); qi++) {
      if (getSkip(qi, false) == false)
        IL->addMer(qi);
    }

    bp = IL->sumLengths();
    delete IL;
  }

  return(bp);
}

inline
u32bit
encodedQuery::bpTotal(void) {
  return(_sequenceLength);
}


//  XXX:  We need to extend get*() to also take a isReverse flag,
//  and to then return the reverse-complement mer, skip, etc.

inline
u64bit
encodedQuery::getMer(u32bit i, bool isReverse) {

  if (isReverse) {
    i = _mersTotal - i - 1;
    return(reverseComplementMer(_merSize, _mers[i]));
  } else {
    return(_mers[i]);
  }
}

inline
bool
encodedQuery::getSkip(u32bit i, bool isReverse) {

  if (isReverse)
    i = _mersTotal - i - 1;

  return(_skip[i]);
}

inline
void
encodedQuery::setSkip(u32bit i, bool isReverse) {

  if (isReverse)
    i = _mersTotal - i - 1;

  //  If skip[i] is already set, don't adjust; otherwise, subtract
  //  one from the available.
  _mersAvail -= 1 - _skip[i];
  _skip[i]    = 1;
}


#endif  //  ENCODEDQUERY_H
