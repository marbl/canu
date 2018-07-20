
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "merStream.H"


merStream::merStream(kMerBuilder *kb, seqStream *ss, bool kbown, bool ssown) {
  _kb       = kb;
  _ss       = ss;

  _kbdelete = kbown;
  _ssdelete = ssown;

  _beg      =  uint64ZERO;
  _end      = ~uint64ZERO;

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
merStream::setBaseRange(uint64 beg, uint64 end) {

  assert(beg < end);

  //fprintf(stderr, "merStream::setBaseRange()-- from "uint64FMT" to "uint64FMT".\n", beg, end);

  //  We can't tell the seqStore when to stop; while we could compute the span of a spaced seed, we
  //  cannot compute it for a compressed seed.  We need to stop iterating when the beginning of the
  //  mer reaches the requested end.

  _ss->setRange(beg, ~uint64ZERO);

  _beg = beg;
  _end = end;

  _kb->clear();

  _invalid = true;
}


uint64
merStream::approximateNumberOfMers(void) {
  uint64  approx = _end - _beg;
  uint64  k      = _kb->merSize();

  //  If we don't know the range, sum all the sequence lengths, otherwise, it's just the length from
  //  begin to end.

  if (_end == ~uint64ZERO) {
    approx = uint64ZERO;

    for (uint32 s=0; s<_ss->numberOfSequences(); s++) {
      uint32 l = _ss->lengthOf(s);

      if (l > k)
        approx += l - k + 1;
    }
  }

  return(approx);
}
