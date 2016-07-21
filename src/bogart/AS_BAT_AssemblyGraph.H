
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

    thickestC   = 0;
    thickest5   = 0;
    thickest3   = 0;
  };
  ~BestPlacement() {
  };

  uint32            tigID;        //  Which tig this is in

  int32             placedBgn;    //  Position in the tig.  Can extend negative.
  int32             placedEnd;    //

  uint32            olapBgn;      //  Position in the tig covered by overlaps.
  uint32            olapEnd;      //

  uint32            thickestC;    //  Read ID of best container
  uint32            thickest5;    //  Read ID of thickest overlap on the 5' end
  uint32            thickest3;    //  Read ID of thickest overlap on the 5' end

  BAToverlapInt     bestC;
  BAToverlapInt     best5;        //  For space reasons, overlaps are stored in the
  BAToverlapInt     best3;        //  internal format.  It's time inefficient though.
};


class AssemblyGraph {
public:
  AssemblyGraph(const char   *prefix,
                double        deviationGraph,
                double        deviationBubble,
                double        deviationRepeat,
                UnitigVector &unitigs);
  
  ~AssemblyGraph() {
    delete [] _placements;
  };

  void          reportGraph(const char *prefix, const char *label);

private:
  vector<BestPlacement>  *_placements;
};



extern AssemblyGraph *AG;

#endif  //  INCLUDE_AS_BAT_ASSEMBLYGRAPH
