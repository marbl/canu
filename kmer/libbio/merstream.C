#include "bio++.H"

merStream::merStream(merStreamFileReader *msf) {
  _ms_mers         = msf;

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

  _fMerP           = &_ms_mers->theFMer();
  _rMerP           = &_ms_mers->theRMer();
}

merStream::merStream(u32bit merSize, seqStream *cs) {
  _ms_mers         = 0L;

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

  _fMerP = &_fMer;
  _rMerP = &_rMer;

  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, seqInCore *seq, u32bit offset, u32bit length) {
  _ms_mers         = 0L;

  _cs_chainedSeq   = 0L;

  _st_string       = seq->sequence() + offset;
  _st_stringPos    = 0;
  _st_stringLen    = length;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  _fMerP = &_fMer;
  _rMerP = &_rMer;

  if (_st_stringLen == 0)
    _st_stringLen = seq->sequenceLength() - offset;

  loadMer(_merSize - 1);
}

//  A specialization of the last constructor, used in wgs-assembler's
//  overmerry to avoid a string allocation, copy and creation of a
//  seqInCore.
merStream::merStream(u32bit merSize, char *seq, u32bit offset, u32bit length) {
  _ms_mers         = 0L;

  _cs_chainedSeq   = 0L;

  _st_string       = seq + offset;
  _st_stringPos    = 0;
  _st_stringLen    = length;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);

  _fMer.clear();
  _rMer.clear();

  _fMerP = &_fMer;
  _rMerP = &_rMer;

  loadMer(_merSize - 1);
}

merStream::~merStream() {
}

bool
merStream::rewind(void) {
  bool  ret = false;

  if (_ms_mers) {
    ret = _ms_mers->rewind();
    _fMer = _ms_mers->theFMer();
    _rMer = _ms_mers->theRMer();
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
