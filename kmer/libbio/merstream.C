#include "bio++.H"

void
merStream::initialize(u32bit merSize) {
  _ss_file         = 0L;
  _cs_stream       = 0L;
  _st_string       = 0L;
  _st_stringPos    = 0;
  _st_stringLen    = 0;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;

  _fMer.clear();
  _rMer.clear();

  _fMer.setMerSize(_merSize);
  _rMer.setMerSize(_merSize);
}

merStream::merStream(u32bit merSize, seqStore *ss) {
  initialize(merSize);
  _ss_file         = ss;
  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, seqStream *cs) {
  initialize(merSize);
  _cs_stream       = cs;
  loadMer(_merSize - 1);
}

merStream::merStream(u32bit merSize, seqInCore *seq, u32bit offset, u32bit length) {
  initialize(merSize);
  _st_string       = seq->sequence() + offset;
  _st_stringLen    = length;
  if (_st_stringLen == 0)
    _st_stringLen = seq->sequenceLength() - offset;
  loadMer(_merSize - 1);
}

//  A specialization of the last constructor, used in wgs-assembler's
//  overmerry to avoid a string allocation, copy and creation of a
//  seqInCore.
merStream::merStream(u32bit merSize, char *seq, u32bit offset, u32bit length) {
  initialize(merSize);
  _st_string       = seq + offset;
  _st_stringLen    = length;
  loadMer(_merSize - 1);
}

merStream::~merStream() {
}

bool
merStream::rewind(void) {
  bool  ret = false;

  if (_ss_file) {
    ret = _ss_file->rewind();
    loadMer(_merSize - 1);
  } else if (_cs_stream) {
    ret = _cs_stream->rewind();
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
