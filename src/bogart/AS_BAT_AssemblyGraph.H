
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
 *    Brian P. Walenz beginning on 2016-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_ASSEMBLYGRAPH
#define INCLUDE_AS_BAT_ASSEMBLYGRAPH

#include "AS_global.H"
#include "AS_BAT_OverlapCache.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"


class BestPlacement {
public:
  BestPlacement() {
    tigID       = UINT32_MAX;

    placedBgn   = INT32_MIN;
    placedEnd   = INT32_MAX;

    olapBgn     = INT32_MIN;
    olapEnd     = INT32_MAX;

    isRepeat    = false;
    isBubble    = false;
    isUnitig    = false;

  };
  ~BestPlacement() {
  };

  uint32            tigID;        //  Which tig this is placed in.

  int32             placedBgn;    //  Position in the tig.  Can extend negative.
  int32             placedEnd;    //

  uint32            olapBgn;      //  Position in the tig covered by overlaps.
  uint32            olapEnd;      //

  bool              isRepeat;     //  This placement is to a repeat
  bool              isBubble;     //  This placement is to a bubble
  bool              isUnitig;     //  This placement is in a unitig

  BAToverlapInt     bestC;
  BAToverlapInt     best5;        //  For space reasons, overlaps are stored in the
  BAToverlapInt     best3;        //  internal format.  It's time inefficient though.
};



class BestReverse {
public:
  BestReverse() {
    readID  = 0;
    placeID = 0;
  };
  BestReverse(uint32 id, uint32 pp) {
    readID  = id;
    placeID = pp;
  };
  ~BestReverse() {
  };

  uint32    readID;    //  readID we have an overlap from; Index into _pForward
  uint32    placeID;   //  index into the vector for _pReverse[readID]
};



class AssemblyGraph {
public:
  AssemblyGraph(const char   *prefix,
                double        deviationGraph,
                double        deviationBubble,
                double        deviationRepeat,
                UnitigVector &unitigs);
  
  ~AssemblyGraph() {
    delete [] _pForward;
    delete [] _pReverse;
  };

  void                      reportGraph(const char *prefix, const char *label);

  vector<BestPlacement>    &getForward(uint32 fi) {
    return(_pForward[fi]);
  };

  vector<BestReverse>      &getReverse(uint32 fi) {
    return(_pReverse[fi]);
  };




private:
  vector<BestPlacement>  *_pForward;   //  Where each read is placed in other tigs
  vector<BestReverse>    *_pReverse;   //  What reads overlap to me
};



extern AssemblyGraph *AG;

#endif  //  INCLUDE_AS_BAT_ASSEMBLYGRAPH