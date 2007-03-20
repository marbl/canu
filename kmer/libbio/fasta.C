#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bio++.H"

FastAFile::FastAFile(const char *filename, u32bit bufferSize) {
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

    _indexname = 0L;

    if (bufferSize != ~u32bitZERO)
      _filebuffer    = new readBuffer(fileno(stdin), _filename, bufferSize);
    else
      _filebuffer    = new readBuffer(fileno(stdin), _filename);
    _isStreamInput = true;
  } else {
    u32bit l = strlen(filename);

    _filename = new char [l + 1];
    strcpy(_filename, filename);

    //  If the filename ends in '.fasta' then append a 'idx',
    //  otherwise, append '.fastaidx'.  Rather complicated for such a
    //  trivial thing....
    //
    if ((l > 5) && (strcmp(filename + l - 6, ".fasta") == 0)) {
      _indexname = new char [l + 4];
      strcpy(_indexname, _filename);
      strcat(_indexname, "idx");
    } else {
      _indexname = new char [l + 10];
      strcpy(_indexname, _filename);
      strcat(_indexname, ".fastaidx");
    }

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

FastAFile::~FastAFile() {
  delete [] _filename;
  delete [] _indexname;
  delete    _filebuffer;
  delete [] _theSeqs;
  delete [] _theNames;
}


u32bit
FastAFile::getNumberOfSequences(void) {
  if (_isRandomAccess)
    return(_theGlobalDesc._numberOfSequences);
  else
    return(0);
}

const char *
FastAFile::getSourceName(void) {
  return(_filename);
}

u32bit
FastAFile::sequenceLength(IID_t iid) {
  if (_isRandomAccess)
    return(_theSeqs[iid]._seqLen);
  else
    return(0);
}

u32bit
FastAFile::headerLength(IID_t iid) {
  if (_isRandomAccess)
    return(_theSeqs[iid]._headerLen);
  else
    return(0);
}

u32bit
FastAFile::currentIID(void) {
  return(_currentSequenceNumber);
}
