#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bio++.H"

FastAWrapper::FastAWrapper(const char *filename, u32bit bufferSize) {
  _theSeqs               = 0L;
  _theNamesLen           = 0;
  _theNames              = 0L;
  _filebuffer            = 0L;
  _filename              = 0L;
  _indexname             = 0L;
  _currentSequenceNumber = 0;
  _isStreamInput         = false;
  _isRandomAccess        = false;
  _isRandomAccessOpt     = false;

  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    _filename = new char [6];
    strcpy(_filename, "stdin");

    if (bufferSize != ~u32bitZERO)
      _filebuffer    = new readBuffer(fileno(stdin), _filename, bufferSize);
    else
      _filebuffer    = new readBuffer(fileno(stdin), _filename);
    _isStreamInput = true;
  } else {
    _filename = new char [strlen(filename)+1];
    strcpy(_filename, filename);

    if (bufferSize != ~u32bitZERO)
      _filebuffer    = new readBuffer(_filename, bufferSize);
    else
#ifdef _AIX
      //  Open with mmap()
      _filebuffer    = new readBuffer(_filename, 0);
#else
      //  Open with the default buffersize
      _filebuffer    = new readBuffer(_filename);
#endif
    _isStreamInput = false;
  }
}

FastAWrapper::~FastAWrapper() {
  delete [] _filename;
  delete [] _indexname;
  delete    _filebuffer;
  delete [] _theSeqs;
  delete [] _theNames;
}


u32bit
FastAWrapper::getNumberOfSequences(void) {
  if (_isRandomAccess)
    return(_theGlobalDesc._numberOfSequences);
  else
    return(0);
}

const char *
FastAWrapper::getSourceName(void) {
  return(_filename);
}

u32bit
FastAWrapper::sequenceLength(IID_t iid) {
  if (_isRandomAccess)
    return(_theSeqs[iid]._seqLen);
  else
    return(0);
}

u32bit
FastAWrapper::headerLength(IID_t iid) {
  if (_isRandomAccess)
    return(_theSeqs[iid]._headerLen);
  else
    return(0);
}

u32bit
FastAWrapper::currentIID(void) {
  return(_currentSequenceNumber);
}
