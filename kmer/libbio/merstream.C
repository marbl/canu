#include "bio++.H"

void
merStream::initialize(kMerBuilder *kb) {
  _kb              = kb;

  _ss_file         = 0L;
  _ss_beg          =  u64bitZERO;
  _ss_end          = ~u64bitZERO;
  _ss_pos          =  u64bitZERO;
  _cs_stream       = 0L;
  _st_string       = 0L;
  _st_stringPos    = 0;
  _st_stringLen    = 0;
  _st_posInSeq     = 0;
  _st_posInStr     = 0;
  _st_num          = 0;

  _merSize         = kb->theFMer().getMerSize();

  _kb->clear();
}

merStream::merStream(kMerBuilder *kb, seqStore *ss) {
  initialize(kb);
  _ss_file         = ss;
}

merStream::merStream(kMerBuilder *kb, seqStream *cs) {
  initialize(kb);
  _cs_stream       = cs;
}

merStream::merStream(kMerBuilder *kb, seqInCore *seq, u32bit offset, u32bit length) {
  initialize(kb);
  _st_string       = seq->sequence() + offset;
  _st_stringLen    = length;
  if (_st_stringLen == 0)
    _st_stringLen = seq->sequenceLength() - offset;
}

//  A specialization of the last constructor, used in wgs-assembler's
//  overmerry to avoid a string allocation, copy and creation of a
//  seqInCore.
merStream::merStream(kMerBuilder *kb, char *seq, u32bit offset, u32bit length) {
  initialize(kb);
  _st_string       = seq + offset;
  _st_stringLen    = length;
}

merStream::~merStream() {
}

bool
merStream::rewind(void) {
  bool  ret = false;

  if (_ss_file) {
    ret = _ss_file->rewind();
  } else if (_cs_stream) {
    ret = _cs_stream->rewind();
  } else {
    _st_stringPos = 0;
    _st_posInSeq  = 0;
    _st_posInStr  = 0;
    _st_num       = 0;
    ret = true;
  }

  _kb->clear();

  return(ret);
}

bool
merStream::setRange(u64bit beg, u64bit end) {
  if (_ss_file) {
    //  We can't tell the seqStore when to stop; while we could
    //  compute the span of a spaced seed, we cannot compute it for a
    //  compressed seed.  We need to stop iterating when the beginning
    //  of the mer reaches the requested end.
    //
    _ss_file->setRange(beg, ~u64bitZERO);
    _ss_beg = beg;
    _ss_end = end;
    _ss_pos = beg;
    _kb->clear();
    return(true);
  }
  return(false);
}
