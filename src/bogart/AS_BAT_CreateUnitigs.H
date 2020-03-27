
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

#ifndef AS_BAT_CREATEUNITIGS_H
#define AS_BAT_CREATEUNITIGS_H

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_MarkRepeatReads.H"  //  confusedEdge

#include "AS_BAT_TigVector.H"


class tigLoc {
public:
  tigLoc() {
    cID  = UINT32_MAX;
    cBgn = 0;
    cEnd = 0;
    uID  = UINT32_MAX;
  };

  uint32   cID;
  uint32   cBgn;
  uint32   cEnd;
  uint32   uID;    //  Debugging.
};



void
createUnitigs(TigVector             &contigs,
              TigVector             &unitigs,
              uint32                 minIntersectLen,
              uint32                 maxPlacements,
              vector<confusedEdge>  &confusedEdges,
              vector<tigLoc>        &unitigSource);

#endif  //  AS_BAT_CREATEUNITIGS_H
