#include "bri++.H"

merStream::merStream(u32bit merSize, const char *filename) {
  _theFile         = new FastAstream(filename);
  _theSeq          = 0L;
  _theLen          = 0;
  _thePos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = u64bitMASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, char const *seq, u32bit len) {
  _theFile         = 0L;
  _theSeq          = seq;
  _theLen          = len;
  _thePos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = u64bitMASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  loadMer(_merSize - 1);
}

merStream::~merStream() {
  delete _theFile;
}

char const *
merStream::theFMerString(void) {
  for (u32bit i=0; i<_merSize; i++)
    _theMerString[_merSize-i-1] = decompressSymbol[(_theFMer >> (2*i)) & 0x03];
  _theMerString[_merSize] = 0;
  return(_theMerString);
}

char const *
merStream::theRMerString(void) {
  for (u32bit i=0; i<_merSize; i++)
    _theMerString[_merSize-i-1] = decompressSymbol[(_theRMer >> (2*i)) & 0x03];
  _theMerString[_merSize] = 0;
  return(_theMerString);
}
