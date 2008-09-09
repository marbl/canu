#ifndef MERSTREAM_H
#define MERSTREAM_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"
#include "seqStream.H"

//
//  merStream needs exclusive use of a kMerBuilder and a seqStream.
//
//  The kMerBuilder can be used over and over.  I think snapper is the
//  only one that does this though.
//
//  The seqStream can be used elsewhere, but ONLY for looking up
//  positions.
//
//  The merStream does NOT assume ownership of either of these, unless
//  the own flags are set.
//

class merStream {
public:
  merStream(kMerBuilder *kb, seqStream *ss, bool kbown=false, bool ssown=false);
  ~merStream();

  kMer const &           theFMer(void)    { return(_kb->theFMer()); };
  kMer const &           theRMer(void)    { return(_kb->theRMer()); };
  kMer const &           theCMer(void)    { return(_kb->theCMer()); };

  bool                   nextMer(u32bit skip=0) {
    char  ch;

    do {
      ch = _ss->get();
      if (ch == 0)
        return(false);
    } while ((_kb->addBase(ch) == true) || (skip-- > 0));

    _kb->mask();

    //  _ss->acgtPosForMerStream(theFMer().getMerSpan())
    if (_ss->seqPos() - theFMer().getMerSpan() >= _end)
      return(false);

    return(true);
  };

  //  Everyone can rewind back to the start, and it's boring.
  //  Slightly more efficient than just opening a new stream.
  //
  //  If we're backed by a seqStore, we can reposition the merstream
  //  anywhere we like to.  The positions refer to the N'th ACGT in
  //  the file, not the N'th mer in the file!  Note that rewind()ing a
  //  file with a range set rewinds to the beginning of the _range_
  //  not the file.
  //
  void                   rewind(void);
  void                   setRange(u64bit beg, u64bit end);

  u64bit                 thePositionInSequence(void)   { return(_ss->seqPos() + 1 - theFMer().getMerSpan()); };
  u64bit                 thePositionInStream(void)     { return(_ss->strPos() + 1 - theFMer().getMerSpan()); };
  u64bit                 theSequenceNumber(void)       { return(_ss->seqIID()); };

  u64bit                 approximateNumberOfMers(void);

private:
  kMerBuilder          *_kb;
  seqStream            *_ss;

  bool                  _kbdelete;
  bool                  _ssdelete;

  //  We can't tell the seqStore when to stop; while we could compute
  //  the span of a spaced seed, we cannot compute it for a compressed
  //  seed.  We need to stop iterating when the beginning of the mer
  //  reaches the requested end.
  //
  u64bit                _beg;
  u64bit                _end;
  u64bit                _pos;
};


#endif  //  MERSTREAM_H
