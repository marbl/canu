#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bio++.H"


FastASequenceInCore::FastASequenceInCore(IID_t iid) {
  _idx       = iid;

  _headerLen = 0;
  _headerMax = 0;   //  1024;
  _header    = 0L;  //  new char [_headerMax];

  _seqLen = 0;
  _seqMax = 0;   //  4 * 1024;
  _seq    = 0L;  //  new char [_seqMax];
}


FastASequenceInCore::FastASequenceInCore(IID_t iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen) {
  _idx       = iid;

  _headerLen = hdrlen;
  _headerMax = hdrlen;
  _header    = hdr;

  _seqLen = seqlen;
  _seqMax = seqlen;
  _seq    = seq;
}


FastASequenceInCore::~FastASequenceInCore() {
  delete [] _header;
  delete [] _seq;
}


FastASequenceInCore*
FastAWrapper::getSequence(void) {
  FastASequenceInCore *c = 0L;

  //  Skip whitespace at the start of the sequence.
  //
  while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
    _filebuffer->next();

  if (_filebuffer->eof())
    return(0L);

  //  We should be at a '>' character now.  Fail if not.
  //
  if (_filebuffer->get() != '>') {
    fprintf(stderr, "getSequence()-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, _filebuffer->get());
    exit(1);
  }

  //
  //  If we have an index, we can be quicker about reading things.
  //    If squeezed, suck in the whole thing in two reads
  //    If fixed width, suck in whole lines
  //    Else do character by character
  //

  if        ((_isRandomAccess && _theGlobalDesc._squeezedSequences) ||
             (_isRandomAccess && _theGlobalDesc._fixedWidth)) {

    //  Yea!  We have an index!  Preallocate the space, and do some quick
    //  block reads!
    //
    char *h = new char [ _theSeqs[_currentSequenceNumber]._headerLen + 1 ];
    char *s = new char [ _theSeqs[_currentSequenceNumber]._seqLen + 1 ];
    
    _filebuffer->read(h, _theSeqs[_currentSequenceNumber]._headerLen);
  
    //  Skip any whitespace between the defline and the start of the sequence
    //  (Rather than a seek, we just skip the spaces)
    //
    while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
      _filebuffer->next();

    if (_theGlobalDesc._squeezedSequences) {
      _filebuffer->read(s, _theSeqs[_currentSequenceNumber]._seqLen);
    } else {
      u32bit  bytesRead = 0;
      while (bytesRead < _theSeqs[_currentSequenceNumber]._seqLen) {

        _idxfa_len  lengthToRead = _theGlobalDesc._seqlineLength;

        if (bytesRead + lengthToRead > _theSeqs[_currentSequenceNumber]._seqLen)
          lengthToRead = _theSeqs[_currentSequenceNumber]._seqLen - bytesRead;

        bytesRead += (u32bit)_filebuffer->read(s + bytesRead, lengthToRead);

        //  Skip any whitespace between the defline and the start of the sequence
        //  (Rather than a seek, we just skip the spaces)
        //
        while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
          _filebuffer->next();
      }
    }

    h[_theSeqs[_currentSequenceNumber]._headerLen] = 0;
    s[_theSeqs[_currentSequenceNumber]._seqLen] = 0;

    c = new FastASequenceInCore(_currentSequenceNumber,
                                h, _theSeqs[_currentSequenceNumber]._headerLen,
                                s, _theSeqs[_currentSequenceNumber]._seqLen);
  } else {
    //  Shucks!  We know nothing about the file.  Go character by
    //  character.
    //
    u32bit hLen = 0;
    u32bit hMax = 128;
    u32bit sLen = 0;
    u32bit sMax = 16 * 1024 * 1024;

    //  If we have an index, we know the sizes.  The file must not
    //  be fixed width (or squeezed).
    //
    if (_isRandomAccess) {
      hMax = _theSeqs[_currentSequenceNumber]._headerLen + 1;
      sMax = _theSeqs[_currentSequenceNumber]._seqLen + 1;
    }

    char   *h   = new char [ hMax + 1 ];
    char   *s   = new char [ sMax + 1 ];

    while ((!_filebuffer->eof()) && (_filebuffer->get() != '\r') && (_filebuffer->get() != '\n')) {
      if (hLen >= hMax) {
        hMax += 128;
        char *htmp = new char [hMax + 1];
        memcpy(htmp, h, sizeof(char) * hLen);
        delete [] h;
        h = htmp;
      }
      h[hLen++] = _filebuffer->get();

      _filebuffer->next();
    }
    h[hLen] = 0;

    //  Skip any whitespace between the defline and the start of the sequence
    //
    while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
      _filebuffer->next();

    while ((!_filebuffer->eof()) &&
           (_filebuffer->get() != '>')) {
      if (!whitespaceSymbol[_filebuffer->get()]) {
        if (sLen >= sMax) {
          sMax += 32 * 1024 * 1024;
          char *stmp = new char [sMax + 1];
          memcpy(stmp, s, sizeof(char) * sLen);
          delete [] s;
          s = stmp;
        }
        s[sLen++] = _filebuffer->get();
      }

      _filebuffer->next();
    }
    s[sLen] = 0;

    c = new FastASequenceInCore(_currentSequenceNumber,
                                h, hLen, s, sLen);
  }


  //  Skip whitespace at the end of the sequence.
  //
  while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
    _filebuffer->next();


  //  By reading the sequence, we have implicitly moved to the next
  //  sequence
  //
  _currentSequenceNumber++;

  return(c);
}

