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
//  The stream is not valid until nextMer is called; allowing loops of
//     while (MS->nextMer()) {
//         process(MS->theFMer());
//     }
//
//  setRange() positions refer to ACGT letters in the input, NOT mers.
//  rewind() repositions the file to the start of the range.
//

class merStream {
public:
  merStream(kMerBuilder *kb, seqStream *ss, bool kbown=false, bool ssown=false);
  ~merStream();

  kMer const &           theFMer(void)    { assert(_invalid == false); return(_kb->theFMer()); };
  kMer const &           theRMer(void)    { assert(_invalid == false); return(_kb->theRMer()); };
  kMer const &           theCMer(void)    { assert(_invalid == false); return(_kb->theCMer()); };

  bool                   nextMer(u32bit skip=0) {
    char  ch;

    do {
      ch = _ss->get();
      if (ch == 0)
        return(false);
    } while ((_kb->addBase(ch) == true) || (skip-- > 0));

    _kb->mask();
    _invalid = false;

#if 0
    char   merstring[256];

    fprintf(stderr, "merStream::nextMer()-- seqPos="u64bitFMT" merPos="u64bitFMT" span="u32bitFMT" base0span="u32bitFMT" end="u64bitFMT" %s %s\n",
            _ss->strPos(),
            _ss->strPos() - theFMer().getMerSpan(),
            theFMer().getMerSpan(),
            _kb->baseSpan(0),
            _end,
            _kb->theFMer().merToString(merstring),
            (_ss->strPos() - theFMer().getMerSpan() < _end) ? "" : "STOP");
#endif

    //  The mer is out of range if:
    //    o it begins at or past the _end
    //    o the span of the first base ends at or past the _end
    //
    //  If the mer isn't spaced, the base span is always 1.  If it is spaced, the span will be
    //  between 1 and ... who knows.

    return(_ss->strPos() - theFMer().getMerSpan() + _kb->baseSpan(0) - 1 < _end);
  };

  void                   rewind(void);
  void                   rebuild(void);
  void                   setBaseRange(u64bit beg, u64bit end);

  u64bit                 thePositionInSequence(void)   { assert(_invalid == false); return(_ss->seqPos() - theFMer().getMerSpan()); };
  u64bit                 thePositionInStream(void)     { assert(_invalid == false); return(_ss->strPos() - theFMer().getMerSpan()); };
  u64bit                 theSequenceNumber(void)       { assert(_invalid == false); return(_ss->seqIID()); };

  u64bit                 approximateNumberOfMers(void);

private:
  kMerBuilder          *_kb;
  seqStream            *_ss;

  bool                  _kbdelete;
  bool                  _ssdelete;

  bool                  _invalid;

  u64bit                _beg;
  u64bit                _end;
};


#endif  //  MERSTREAM_H
