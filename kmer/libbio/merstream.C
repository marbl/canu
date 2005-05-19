#include "bio++.H"

merStream::merStream(merStreamFileReader *msf) {
  _ms_mers         = msf;

  _fs_file         = 0L;

  _cs_chainedSeq   = 0L;

  _st_string       = 0L;
  _st_stringPos    = 0;
  _st_stringLen    = 0;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = _ms_mers->merSize();
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  _fMer            = _ms_mers->theFMer();
  _rMer            = _ms_mers->theRMer();
}

merStream::merStream(u32bit merSize, FastAstream *str) {
  _ms_mers         = 0L;

  _fs_file         = str;

  _cs_chainedSeq   = 0L;

  _st_string       = 0L;
  _st_stringPos    = 0;
  _st_stringLen    = 0;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, chainedSequence *cs) {
  _ms_mers         = 0L;

  _fs_file         = 0L;

  _cs_chainedSeq   = cs;

  _st_string       = 0L;
  _st_stringPos    = 0;
  _st_stringLen    = 0;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, const char *seq, u32bit len) {
  _ms_mers         = 0L;

  _fs_file         = 0L;

  _cs_chainedSeq   = 0L;

  _st_string       = seq;
  _st_stringPos    = 0;
  _st_stringLen    = len;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  //  If the bloody user gave us no length, reset _st_stringLen to
  //  be maximum.  nextSymbol() will then stop when it hits
  //  the end of string marker \0.
  //
  if (_st_stringLen == 0)
    _st_stringLen = ~u32bitZERO;

  //  We want to start numbering sequences and stream positions at
  //  zero.  nextSymbol() always increments the sequence number and
  //  the stream position on a defline.  So, we skip the first defline
  //  before we do anything else.
  //
  while (whitespaceSymbol[_st_string[_st_stringPos]] && (_st_string[_st_stringPos]) && (_st_stringPos < _st_stringLen))
    _st_stringPos++;
  if (_st_string[_st_stringPos] == '>')
    while ((_st_string[_st_stringPos] != '\n') && (_st_string[_st_stringPos]) && (_st_stringPos < _st_stringLen))
      _st_stringPos++;

  loadMer(_merSize - 1);
}

merStream::~merStream() {
}


bool
merStream::rewind(void) {
  bool  ret = false;

  if (_ms_mers) {
    ret = _ms_mers->seekToMer(0);
    _fMer = _ms_mers->theFMer();
    _rMer = _ms_mers->theRMer();
  } else if (_fs_file) {
    ret = _fs_file->rewind();
    loadMer(_merSize - 1);
  } else if (_cs_chainedSeq) {
    ret = _cs_chainedSeq->rewind();
    loadMer(_merSize - 1);
  } else {
    ret = true;
    _st_stringPos = 0;
    _st_posInSeq  = 0;
    _st_posInStr  = 0;
    _st_num       = 0;
    loadMer(_merSize - 1);
    _st_posInStr--;
    ret = true;
  }

  return(ret);
}
