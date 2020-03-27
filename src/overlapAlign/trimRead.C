
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"

#include "intervalList.H"
#include "system.H"
#include "sequence.H"

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "alignStats.H"
#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"



bool
testAlignment(char   *aRead,  int32  &abgn,  int32  &aend,  int32  UNUSED(alen),  uint32 Aid,
              char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
              double  maxAlignErate,
              double  maxAcceptErate,
              double &erate);

int32
extend5(char   *aRead,  int32  &abgn,  int32  &aend,  int32         alen,   uint32 Aid,
        char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
        double  maxAlignErate,
        double  maxAcceptErate,
        double &erate);

int32
extend3(char   *aRead,  int32  &abgn,  int32  &aend,  int32         alen,   uint32 Aid,
        char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
        double  maxAlignErate,
        double  maxAcceptErate,
        double &erate);


//  Align the overlap in small blocks to find the largest region that aligns.
//
void
maComputation::trimRead(uint32   minOverlapLength,
                        double   maxErate) {

  intervalList<int32>   clearRange;
  intervalList<int32>   failedRange;

  _readData[_aID].clrBgn = 0;  //  Default to fully trimming the read.
  _readData[_aID].clrEnd = 0;

  if (_overlapsLen == 0)
    return;

  fetchUntrimmedRead(_aID);

  //  PARAMETERS
  //
  //  maxAlignErate  -- allow edlib to generate alignments up to this error rate.
  //  maxAcceptErate -- accept edlib alignments if they are below this error rate.
  //
  //  extension      -- 'step size' for extending the validated overlap region; probably should be even.
  //  anchor         -- extend the next region this many bases into the previous (validated) region
  //  overhang       -- require overlaps to overlap by this much before merging regions.
  //
  //  bextra         -- for the initial alignment, extend the B region by this amount on each end
  //  slop           --

  double   maxAlignErate  = maxErate + 0.2;
  double   maxAcceptErate = maxErate;

  int32    extension      = 500;
  int32    anchor         = 500;
  int32    overhang       = 500;

  int32    bextra         = maxAcceptErate * extension + 1250;
  int32    slop           = maxAcceptErate * extension;


  if (_verboseTrim > 1) {
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "BEGIN trimming for read %u\n", _aID);
  }

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    assert(_aID == ov->a_iid);
    _bID = ov->b_iid;

    double  erate     = 0;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;
    int32   abgn      =                             a5;
    int32   aend      = _readData[_aID].rawLength - a3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;
    int32   bbgn      =                             b5;
    int32   bend      = _readData[_bID].rawLength - b3;

    //  Log.

    //
    //  Decide if this overlap will help us validate the A read.  If it's contained in
    //  any range, it's not useful.
    //

    bool  useful = true;

    for (uint32 ii=0; ii<clearRange.numberOfIntervals(); ii++) {
      int32  lo = clearRange.lo(ii) - overhang;
      int32  hi = clearRange.hi(ii) + overhang;

      if ((abgn < lo) ||   //  Useful, extends this range.
          (hi < aend))
        continue;

      useful = false;      //  Not useful.  Contained in lo .. hi.

      if (_verboseTrim > 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "  notU A %8u %7d-%-7d   B %8u %7d-%-7d (contained in range %d-%d)\n",
                _aID, abgn, aend, _bID, bbgn, bend, lo, hi);
      }

      break;
    }

    if (useful == false)
      continue;

    //
    //  Decide if this overlap is spanning a known to be bad region.  If many
    //  other overlaps have failed through here, just give up.
    //
    //  This isn't perfect, since this overlap could be spanning a bad region
    //  AND extending into good sequence.
    //

    bool  giveup = false;

    for (uint32 ii=0; ii<failedRange.numberOfIntervals(); ii++) {
      int32  lo = failedRange.lo(ii);
      int32  hi = failedRange.hi(ii);

      if (_verboseTrim > 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "  ALIGN FAIL          %7d-%-7d\n", lo, hi);
      }

      if ((abgn < lo) &&
          (hi   < aend))
        giveup = true;
    }

    if (giveup == true) {
      if (_verboseTrim > 1) {
      fprintf(stderr, "\n");
        fprintf(stderr, "  BAD  %8u A %7d-%-7d B %8u %7d-%-7d\n", _aID, abgn, aend, _bID, bbgn, bend);
      }
      continue;
    }

    //  It's useful!

    if (_verboseTrim > 1) {
      fprintf(stderr, "\n");
      fprintf(stderr, "  OLAP A %8u %7d-%-7d   B %8u %7d-%-7d\n", _aID, abgn, aend, _bID, bbgn, bend);
    }

    //  Grab the sequence and orient it for this overlap.

    fetchUntrimmedRead(ov->b_iid, false, ov->flipped());

    //  Find the expected middle of the alignment in both the A and B reads.

    int32  amid = abgn + (aend - abgn) / 2;
    int32  bmid = bbgn + (bend - bbgn) / 2;

    //  Extend the A region on each side, and the B region by a little bit more (since the end gaps are free).

    abgn = max(amid - extension / 2, 0);
    aend = min(amid + extension / 2, _readData[_aID].rawLength);

    bbgn = max(bmid - extension / 2 - bextra, 0);
    bend = min(bmid + extension / 2 + bextra, _readData[_bID].rawLength);

    //  Align to find the precise region we align to in the B read.  If the alignment passes,
    //  bbgn and bend are updated to those coordinates.

    if (_verboseTrim > 1) {
      fprintf(stderr, "  test align      %7d-%-7d   B %8u %7d-%-7d\n", abgn, aend, _bID, bbgn, bend);
    }

    //  XXX  if it fails, shift left/right until we find the seed.
    if (testAlignment(_aRead, abgn, aend, _readData[_aID].rawLength, _aID,
                      _bRead, bbgn, bend, _readData[_bID].rawLength, _bID,
                      maxAlignErate,
                      maxAcceptErate,
                      erate) == 0) {
      if (_verboseTrim > 1)
        fprintf(stderr, "  trim %8u fails.\n", _aID);
      continue;
    }

    if (_verboseTrim > 1) {
      fprintf(stderr, "  INIT align      %7d-%-7d   B %8u %7d-%-7d\n", abgn, aend, _bID, bbgn, bend);
    }

    //
    //
    //  Flag regions that fail alignment.  When enough overlaps fail there, just give up.
    //
    //




    //  Until we hit a block of low quality or we exhaust the supply of
    //  bases, extend out from that seed in both directions.  Save the
    //  'clear' region as the last (extended) region with a good alignment.
    //
    //  This fails at the ends of the B read; we're trying to align an A
    //  chunk of ~1000 bases into a B chunk of less than ~1000 bases.
    //  However, since this is generally internal to the A read, we don't
    //  care.


    //  Extend toward the 5' end.
    erate = 0.0;

    while ((erate < maxAcceptErate) &&
           (abgn > 0) &&
           (bbgn > 0)) {
      int32   ab = max(abgn - extension, 0);
      int32   ae = min(abgn + anchor,    _readData[_aID].rawLength);

      int32   bb = bbgn - extension - slop;
      int32   be = bbgn + anchor    + slop;

      //  If the B chunk is too small for the A chunk, stop.  We should align
      //  B into A, but that's a lot of work.

      if ((bb < 0) ||
          (be > _readData[_bID].rawLength)) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext5 test     %7d-%-7d   B %8u %7d-%-7d    too-small STOP\n", ab, ae, _bID, bb, be);
        break;
      }

      if (_verboseTrim > 1)
        fprintf(stderr, "    ext5 test     %7d-%-7d   B %8u %7d-%-7d\n", ab, ae, _bID, bb, be);

      //  Returns -1 if the whole alignment is crap
      //  Returns  0 if some was good
      //  Returns  1 if all was good

      int32  ext = extend5(_aRead, ab, ae, _readData[_aID].rawLength, _aID,
                           _bRead, bb, be, _readData[_bID].rawLength, _bID,
                           maxAlignErate,
                           maxAcceptErate,
                           erate);

      if (ext == -1) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext5          %7d-%-7d   B %8u %7d-%-7d error %6.2f FAIL\n", ab, ae, _bID, bb, be, erate);
        failedRange.add(abgn, aend - abgn);
      }

      if (ext == 0) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext5          %7d-%-7d   B %8u %7d-%-7d error %6.2f GOOD\n", ab, ae, _bID, bb, be, erate);
        abgn = ab;
        bbgn = bb;
      }

      if (ext == 1) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext5          %7d-%-7d   B %8u %7d-%-7d error %6.2f STOP\n", ab, ae, _bID, bb, be, erate);
        abgn = ab;
        bbgn = bb;
        break;
      }

      assert(0 <= ab);
      assert(0 <= bb);

      assert(ab <= _readData[_aID].rawLength);
      assert(bb <= _readData[_bID].rawLength);
    }

    //if (_verboseTrim > 1)
    //  fprintf(stderr, "    ext5          %7d-%-7d\n", abgn, aend);



    //  Extend toward the 3' end.
    erate = 0.0;

    while ((erate < maxAcceptErate) &&
           (aend < _readData[_aID].rawLength) &&
           (bend < _readData[_bID].rawLength)) {
      int32   ab = max(aend - anchor,    0);
      int32   ae = min(aend + extension, _readData[_aID].rawLength);

      int32   bb = bend - anchor    - slop;
      int32   be = bend + extension + slop;

      if ((bb < 0) ||
          (be > _readData[_bID].rawLength)) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext3 test     %7d-%-7d   B %8u %7d-%-7d    too-small STOP\n", ab, ae, _bID, bb, be);
        break;
      }

      if (_verboseTrim > 1)
        fprintf(stderr, "    ext3 test     %7d-%-7d   B %8u %7d-%-7d\n", ab, ae, _bID, bb, be);

      int32  ext = extend3(_aRead, ab, ae, _readData[_aID].rawLength, _aID,
                           _bRead, bb, be, _readData[_bID].rawLength, _bID,
                           maxAlignErate,
                           maxAcceptErate,
                           erate);

      if (ext == -1) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext3          %7d-%-7d   B %8u %7d-%-7d error %6.2f FAIL\n", ab, ae, _bID, bb, be, erate);
        failedRange.add(abgn, aend - abgn);
      }

      if (ext == 0) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext3          %7d-%-7d   B %8u %7d-%-7d error %6.2f GOOD\n", ab, ae, _bID, bb, be, erate);
        aend = ae;
        bend = be;
      }

      if (ext == 1) {
        if (_verboseTrim > 1)
          fprintf(stderr, "    ext3          %7d-%-7d   B %8u %7d-%-7d error %6.2f STOP\n", ab, ae, _bID, bb, be, erate);
        aend = ae;
        bend = be;
        break;
      }

      assert(0 <= ae);
      assert(0 <= be);

      assert(ae <= _readData[_aID].rawLength);
      assert(be <= _readData[_bID].rawLength);
    }

    //if (_verboseTrim > 1)
    //  fprintf(stderr, "    ext3          %7d-%-7d\n", abgn, aend);


    //  Extend the saved A clear range based on this overlap.  We can't use
    //  the B read alignment (for trimming the B read) without storing a set
    //  of disjoint clear ranges for each read.

    abgn += overhang;
    aend -= overhang;

    if (abgn < aend) {
      if (_verboseTrim > 1)
        fprintf(stderr, "  clr  %8u A %7d-%-7d (excluding overhang requirement)\n", _aID, abgn, aend);

      clearRange.add(abgn, aend - abgn);
      clearRange.merge();

      if (_verboseTrim > 1)
        for (uint32 ii=0; ii<clearRange.numberOfIntervals(); ii++)
          fprintf(stderr, "  CLR             %7d-%-7d\n", clearRange.lo(ii) - overhang, clearRange.hi(ii) + overhang);
    }

    else {
      if (_verboseTrim > 1)
        fprintf(stderr, "  clr  %8u A %7d-%-7d  TOO SMALL!\n", _aID, abgn, aend);
    }
  }

  //  Pick the largest clear region and save it in the globals.  We could let
  //  the 'writer' thread do this, but we only care about this in memory at
  //  the moment, it's small, and it's easier to write in one bug chunk.

  for (uint32 ii=0; ii<clearRange.numberOfIntervals(); ii++) {
    uint32 bgn = clearRange.lo(ii) - overhang;
    uint32 end = clearRange.hi(ii) + overhang;

    if ((_verboseTrim > 0) && (clearRange.numberOfIntervals() > 1))
      fprintf(stderr, "  CLR %8u   %7d-%-7d - %6.2f%% clear\n",
              _aID,
              bgn, end,
              100.0 * (end - bgn) / _readData[_aID].rawLength);

    if (end - bgn > _readData[_aID].clrEnd - _readData[_aID].clrBgn)
      _readData[_aID].setClear(bgn, end);
  }

  if (_verboseTrim > 0)
    fprintf(stderr, "  CLR %8u   %7d-%-7d - %6.2f%% clear (final)\n",
            _aID,
            _readData[_aID].clrBgn, _readData[_aID].clrEnd,
            100.0 * (_readData[_aID].clrEnd - _readData[_aID].clrBgn) / _readData[_aID].rawLength);
}
