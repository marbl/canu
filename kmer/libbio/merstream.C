#include "bri++.H"

merStream::merStream(merStreamFileReader *msf) {
  _theMers         = msf;
  _theFile         = 0L;
  _theString       = 0L;
  _theStringPos    = 0;
  _theStringLen    = 0;
  _thePos          = 0;
  _theNum          = 0;

  _merSize         = msf->merSize();
  _timeUntilValid  = 0;
  _theMerMask      = u64bitMASK(_merSize << 1);

  _theFMer         = _theMers->theFMer();
  _theRMer         = _theMers->theRMer();

  _theRMerShift    = (_merSize << 1) - 2;
}

merStream::merStream(u32bit merSize, const char *filename) {
  _theMers         = 0L;
  _theFile         = new FastAstream(filename);
  _theString       = 0L;
  _theStringPos    = 0;
  _theStringLen    = 0;
  _thePos          = 0;
  _theNum          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = u64bitMASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, const char *seq, u32bit len) {
  _theMers         = 0L;
  _theFile         = 0L;
  _theString       = seq;
  _theStringPos    = 0;
  _theStringLen    = len;
  _thePos          = 0;
  _theNum          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = u64bitMASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  //  If the bloody user gave us no length, reset _theLen to
  //  be maximum.  nextSymbol() will then stop when it hits
  //  the end of string marker \0.
  //
  if (_theStringLen == 0)
    _theStringLen = ~u32bitZERO;

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
