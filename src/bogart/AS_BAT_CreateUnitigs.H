
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
 *    Brian P. Walenz beginning on 2016-OCT-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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
