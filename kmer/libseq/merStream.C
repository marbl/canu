#include "merStream.H"


merStream::merStream(kMerBuilder *kb, seqStream *ss, bool kbown, bool ssown) {
  _kb       = kb;
  _ss       = ss;

  _kbdelete = kbown;
  _ssdelete = ssown;

  _beg      =  u64bitZERO;
  _end      = ~u64bitZERO;

  _kb->clear();

  _invalid = true;
}


merStream::~merStream() {
  if (_kbdelete)  delete _kb;
  if (_ssdelete)  delete _ss;
}


void
merStream::rewind(void) {
  _ss->rewind();
  _kb->clear();
  _invalid = true;
}


void
merStream::rebuild(void) {
  _ss->setPosition(_ss->strPos() - _kb->theFMer().getMerSpan());
  _kb->clear();
  _invalid = true;
}


void
merStream::setBaseRange(u64bit beg, u64bit end) {

  assert(beg < end);

  //fprintf(stderr, "merStream::setBaseRange()-- from "u64bitFMT" to "u64bitFMT".\n", beg, end);

  //  We can't tell the seqStore when to stop; while we could compute the span of a spaced seed, we
  //  cannot compute it for a compressed seed.  We need to stop iterating when the beginning of the
  //  mer reaches the requested end.

  _ss->setRange(beg, ~u64bitZERO);

  _beg = beg;
  _end = end;

  _kb->clear();

  _invalid = true;
}


u64bit
merStream::approximateNumberOfMers(void) {
  u64bit  approx = _end - _beg;
  u64bit  k      = _kb->merSize();

  //  If we don't know the range, sum all the sequence lengths, otherwise, it's just the length from
  //  begin to end.

  if (_end == ~u64bitZERO) {
    approx = u64bitZERO;

    for (u32bit s=0; s<_ss->numberOfSequences(); s++) {
      u32bit l = _ss->lengthOf(s);

      if (l > k)
        approx += l - k + 1;
    }
  }

  return(approx);
}
