
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

#ifndef INCLUDE_AS_BAT_ASSEMBLYGRAPH
#define INCLUDE_AS_BAT_ASSEMBLYGRAPH

#include "runtime.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  For ReadEnd
#include "AS_BAT_PlaceReadUsingOverlaps.H"
#include "AS_BAT_Unitig.H"            //  For SeqInterval
#include "AS_BAT_TigVector.H"


class BestPlacement {
public:
  BestPlacement(overlapPlacement  &pl,
                BAToverlap        *ovl,
                uint32             tc,
                uint32             t5,
                uint32             t3) : tigID(pl.tigID),
                                         olapMin(pl.verified.min()),
                                         olapMax(pl.verified.max()) {
    if (tc < UINT32_MAX)   bestC = ovl[tc];
    if (t5 < UINT32_MAX)   best5 = ovl[t5];
    if (t3 < UINT32_MAX)   best3 = ovl[t3];

    if (bestC.b_iid != 0) {                 //  Simple sanity check.  Ensure that contained
      assert(best5.b_iid == 0);             //  edges have no dovetail edges.  This screws up
      assert(best3.b_iid == 0);             //  the logic when outputting the graph.
    }

    assert((bestC.a_hang <= 0) && (bestC.b_hang >= 0));  //  ALL contained edges should be this.
    assert((best5.a_hang <= 0) && (best5.b_hang <= 0));  //  ALL 5' edges should be this.
    assert((best3.a_hang >= 0) && (best3.b_hang >= 0));  //  ALL 3' edges should be this.
  };

  uint32       tigID   = UINT32_MAX;   //  Which tig this is placed in.

  int32        olapMin = INT32_MIN;    //  Position in the tig covered by overlaps.
  int32        olapMax = INT32_MAX;    //

  BAToverlap   bestC;                  //  Used internally to generate back edges,
  BAToverlap   best5;                  //  and in AS_BAT_DropDeadEnds.C
  BAToverlap   best3;
};



class BestReverse {
public:
  BestReverse(uint32 id=0, uint32 pp=0) : readID(id),
                                          placeID(pp) {
  };

  uint32       readID;    //  readID we have an overlap from; Index into _pForward
  uint32       placeID;   //  index into the vector for _pForward[readID]
};



class AssemblyGraph {
public:
  AssemblyGraph(const char   *prefix,
                double        deviationRepeat,
                double        repeatLimit,
                TigVector    &tigs) {
    buildGraph(prefix, deviationRepeat, repeatLimit, tigs);
  }

  ~AssemblyGraph() {
    delete [] _pForward;
    delete [] _pReverse;
  };

  vector<BestPlacement>    &getForward(uint32 fi) const  { return(_pForward[fi]); };
  vector<BestReverse>      &getReverse(uint32 fi) const  { return(_pReverse[fi]); };

  void                      buildGraph(const char   *prefix,
                                       double        deviationRepeat,
                                       double        repeatLimit,
                                       TigVector    &tigs);

private:
  vector<BestPlacement>  *_pForward;   //  Where each read is placed in other tigs
  vector<BestReverse>    *_pReverse;   //  What reads overlap to me
};




#endif  //  INCLUDE_AS_BAT_ASSEMBLYGRAPH
