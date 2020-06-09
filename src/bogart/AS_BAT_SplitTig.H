
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

#ifndef INCLUDE_AS_BAT_SPLIT_TIG
#define INCLUDE_AS_BAT_SPLIT_TIG

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_MarkRepeatReads.H"

#include "intervalList.H"
#include "stddev.H"

#include <vector>

using namespace std;


//  We can split a tig based on a list of repeat/unique intervals, or based
//  on a list of read ends.
//
//  A list of breakPointCoords will split a tig into multiple tigs (and label
//  those new tigs as repeat or unique) based on how reads intersect the supplied
//  intervals.  A read must be fully contained in a repeat interval for it to be
//  moved to the new tig, while it only needs to touch a unique interval for it
//  to end up there.
//
//  A list of breakPointEnd will slice the tig at that read end.
//    _bgn == true  -> reads that begin at after position are in the region
//    _end == false -> reads that end before position are in the region

class breakPointCoords {
public:
  breakPointCoords(int32 bgn, int32 end, bool rpt=false) {
    _bgn = bgn;
    _end = end;
    _rpt = rpt;
  };
  ~breakPointCoords() {
  };

  bool    operator<(breakPointCoords const &that) const {
    return(_bgn < that._bgn);
  };

  int32   _bgn;
  int32   _end;
  bool    _rpt;
};



class breakReadEnd {
public:
  breakReadEnd(uint32 rid,
               bool   highend,
               int32  point,
               int32  regionBgn,
               int32  regionEnd) {
    readID     = rid;
    breakHigh  = highend;
    splitCoord = point;
    repeatBgn  = regionBgn;
    repeatEnd  = regionEnd;
  };
  ~breakReadEnd() {
  };

  uint32  readID;
  bool    breakHigh;
  int32   splitCoord;
  int32   repeatBgn;
  int32   repeatEnd;
};



class breakPointEnd {
public:
  breakPointEnd(uint32 tigID, uint32 pos, bool bgn) {
    _tigID    = tigID;
    _pos      = pos;
    _bgn      = bgn;
  };
  ~breakPointEnd() {
  };

  bool     operator<(breakPointEnd const &that) const {
    uint64  a =      _tigID;  a <<= 32;  a |=      _pos;  a <<= 1;  a |=      _bgn;  //  Because _tigID is 32-bit
    uint64  b = that._tigID;  b <<= 32;  b |= that._pos;  b <<= 1;  b |= that._bgn;

    return(a < b);
  };

  bool     operator==(breakPointEnd const &that) const {
    uint64  a =      _tigID;  a <<= 32;  a |=      _pos;  a <<= 1;  a |=      _bgn;  //  Because _tigID is 32-bit
    uint64  b = that._tigID;  b <<= 32;  b |= that._pos;  b <<= 1;  b |= that._bgn;

    return(a == b);
  };

  uint32  _tigID;
  uint32  _pos;
  bool    _bgn;
};



uint32
splitTigByRegions(TigVector                &tigs,
                  Unitig                   *tig,
                  vector<breakPointCoords> &BP,
                  Unitig                  **newTigs,
                  int32                    *lowCoord,
                  uint32                   *nRepeat,
                  uint32                   *nUnique,
                  bool                      doMove);

uint32
splitTigAtReadEnds(TigVector                &tigs,
                   Unitig                   *tig,
                   vector<breakPointEnd>    &BP,
                   Unitig                  **newTigs,
                   int32                    *lowCoord,
                   bool                      doMove);

void
splitTigAtReadEnds(TigVector                &tigs,
                   Unitig                   *tig,
                   vector<breakReadEnd>     &BP,
                   intervalList<int32>      &tigMarksR);

#endif  //  INCLUDE_AS_BAT_SPLIT_TIG
