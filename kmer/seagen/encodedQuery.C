#include <new>

#include "encodedQuery.H"

encodedQuery::encodedQuery(seqInCore  *S,
                           uint32      k) {

  _iid             = S->getIID();
  _sequenceLength  = S->sequenceLength();
  _merSize         = k;
  _mersTotal       = 0;
  _mersAvail       = 0;
  _mers            = 0L;
  _skip            = 0L;
  _numberOfResults = 0;
  _output          = 0L;
  _outputLen       = 0;
  _outputMax       = 0;

  if (k > _sequenceLength)
    return;

  char const           *seq    = S->sequence();
  uint32                seqLen = S->sequenceLength();

  _mersTotal = seqLen - k + 1;
  _mersAvail = 0;
  _mers      = new uint64 [_mersTotal];
  _skip      = new uint8  [_mersTotal];

  uint64   substring      = uint64ZERO;
  uint64   mermask        = uint64MASK(2 * k);
  int32   timeUntilValid = k;

  for (uint32 i=0; i<seqLen; ) {
    substring <<= 2;
    substring  &= mermask;

    if (letterToBits[seq[i]] != 0xff) {
      substring |= letterToBits[ seq[i] ];
      timeUntilValid--;
    } else {
      timeUntilValid = k;
    }

    i++;

    if (i >= k) {
      _mers[i-k]  = substring;
      _skip[i-k]  = timeUntilValid > 0;
      _mersAvail += 1 - _skip[i-k];
    }
  }
}

encodedQuery::~encodedQuery() {
  delete [] _mers;
  delete [] _skip;
  delete [] _output;
}




void
encodedQuery::test(seqInCore *S) {

  //  We assume we've been initialized with the forward version!

  uint32                k = _merSize;

  char const           *seq    = S->sequence();
  uint32                seqLen = S->sequenceLength();

  uint64   substring      = uint64ZERO;
  uint64   mermask        = uint64MASK(2 * k);
  int32   timeUntilValid = k;

  //  Compute the complement version; we'll iterate through all data
  //  in us, comparing against what the original method would say.

  uint32    _r_mersAvail = 0;
  uint64   *_r_mers      = new uint64 [_mersTotal];
  uint8    *_r_skip      = new uint8  [_mersTotal];

  substring      = uint64ZERO;
  mermask        = uint64MASK(2 * k);
  timeUntilValid = k;

  for (uint32 i=0; i<seqLen; ) {
    substring <<= 2;
    substring  &= mermask;

    if (letterToBits[seq[seqLen - 1 - i]] != 0xff) {
      substring |= letterToBits[ complementSymbol[ seq[seqLen - 1 - i] ]];
      timeUntilValid--;
    } else {
      timeUntilValid = k;
    }

    i++;

    if (i >= k) {
      _r_mers[i-k]  = substring;
      _r_skip[i-k]  = timeUntilValid > 0;
      _r_mersAvail += 1 - _r_skip[i-k];
    }
  }

#if 0
  //  For comparison, this is the original code used to compute the
  //  reverse complement mers.

  for (uint32 i=0; i<seqLen; ) {
    substring <<= 2;
    substring  &= mermask;

    if (validSymbol[seq[seqLen - 1 - i]]) {
      substring |= compressSymbol[ complementSymbol[ seq[seqLen - 1 - i] ]];
      timeUntilValid--;
    } else {
      timeUntilValid = k;
    }

    i++;

    if (i >= k) {
      _mers[i-k]  = substring;
      _skip[i-k]  = timeUntilValid > 0;
      _mersAvail += 1 - _skip[i-k];
    }
  }
#endif





  //  CHECK!
  //
  if (_r_mersAvail != _mersAvail) {
    fprintf(stderr, "encodedQuery::test()-- mersAvail incorrect:  Recomputed:"uint32FMT"  Real:"uint32FMT"\n", _mersAvail, _r_mersAvail);
  }

  char  mer1[65];
  char  mer2[65];
  bool  fail = false;

  for (uint32 i=0; i<_mersTotal; i++) {

    if (getSkip(i, true) != _r_skip[i]) {
      fprintf(stderr, "encodedQuery::test()-- skip["uint32FMTW(4)"] incorrect:  Acc:%d  Real:%d\n", i, getSkip(i, true), _r_skip[i]);
      fail = true;
    }

    if (getSkip(i, true) == false) {
      if (getMer(i, true) != _r_mers[i]) {
        uint64ToMerString(_merSize, getMer(i, true), mer1);
        uint64ToMerString(_merSize, _r_mers[i], mer2);
        fprintf(stderr, "encodedQuery::test()-- mers["uint32FMTW(4)"] incorrect:  Acc:"uint64HEX" %s   Real:"uint64HEX" %s\n",
                i,
                getMer(i, true), mer1,
                _r_mers[i], mer2);
        fail = true;
      }
    }

    if (fail) {
      char   rev[2048];
      strcpy(rev, seq);
      fprintf(stderr, "seq='%s'\n", seq);
      fprintf(stderr, "rev='%s'\n", reverseComplementSequence(rev, seqLen));
      exit(1);
    }
  }

  //fprintf(stderr, "encodedQuery::test()-- %s\n", seq);
  //fprintf(stderr, "encodedQuery::test()-- tested avail:"uint32FMT" total:"uint32FMT"\n", _mersAvail, _mersTotal);

  delete [] _r_mers;
  delete [] _r_skip;
}


void
encodedQuery::addOutput(void *newout, uint32 size) {

  //  Allocate space for the output -- 1MB should be enough for about
  //  29000 signals.  Make it 32K -> 900 signals.
  //
  //  You probably do not want to move this into the query
  //  constructor, as that will just waste a lot of memory with
  //  thousands of these in the input queue.
  //
  if (_output == 0L) {
    _outputLen = 0;
    _outputMax = 32 * 1024;
    _output    = new char [_outputMax];
  }

  if (_outputLen + 128 >= _outputMax) {
    _outputMax <<= 1;
    char *o = 0L;

    try {
      o = new char [_outputMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "encodedQuery::addOutput()-- out of memory, tried to extend output string\n");
      fprintf(stderr, "encodedQuery::addOutput()-- from "uint32FMT" to "uint32FMT" bytes.\n",
              _outputLen, _outputMax);
      exit(1);
    }

    memcpy(o, _output, _outputLen);
    delete [] _output;
    _output = o;
  }

  if (size > 0) {
    memcpy(_output + _outputLen, newout, size);
    _outputLen += size;
  } else {
    char *n = (char *)newout;

    while (*n)
      _output[_outputLen++] = *n++;
    _output[_outputLen] = 0;
  }

  _numberOfResults++;
};
