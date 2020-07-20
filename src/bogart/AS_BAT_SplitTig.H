
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


//  Specifies a read end to break a tig at.  No other reads are allowed
//  to pass this read end in the post-split tigs.
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



void
splitTigAtReadEnds(TigVector                &tigs,
                   Unitig                   *tig,
                   vector<breakReadEnd>     &BP,
                   intervalList<int32>      &tigMarksR);

#endif  //  INCLUDE_AS_BAT_SPLIT_TIG
