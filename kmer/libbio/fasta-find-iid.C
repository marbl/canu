#include "bio++.H"

#include <stdio.h>
#include <string.h>


bool
FastAWrapper::find(IID_t  iid) {

  if (!_isRandomAccess) {
    fprintf(stderr, "FastAWrapper::find(IID)-- ERROR: '%s' is not open for random access; find() failed.\n", _filename);
    return(false);
  }

  if (iid >= _theGlobalDesc._numberOfSequences) {
    fprintf(stderr, "FastAWrapper::find(IID)-- ERROR: index of "u32bitFMT" too large for '%s' (only "u32bitFMT" sequences).\n",
            iid, _filename, _theGlobalDesc._numberOfSequences);
    return(false);
  }

  if (iid != _currentSequenceNumber) {
    _filebuffer->seek(_theSeqs[iid]._headerStart);
    _currentSequenceNumber = iid;
  }

  return(true);
}
