
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

#include "runtime.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  For ReadEnd
#include "AS_BAT_Unitig.H"            //  For SeqInterval
#include "AS_BAT_TigVector.H"


class BestPlacement {
public:
  BestPlacement() {
    tigID       = UINT32_MAX;

    olapMin     = INT32_MIN;
    olapMax     = INT32_MAX;
  };
  ~BestPlacement() {
  };

  uint32            tigID;        //  Which tig this is placed in.

  int32             olapMin;      //  Position in the tig covered by overlaps.
  int32             olapMax;      //

  BAToverlap        bestC;        //  Used internally to generate back edges,
  BAToverlap        best5;        //  and in AS_BAT_DropDeadEnds.C
  BAToverlap        best3;
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
  uint32    placeID;   //  index into the vector for _pForward[readID]
};



class AssemblyGraph {
public:
  AssemblyGraph(const char   *prefix,
                double        deviationRepeat,
                TigVector    &tigs) {
    buildGraph(prefix, deviationRepeat, tigs);
  }

  ~AssemblyGraph() {
    delete [] _pForward;
    delete [] _pReverse;
  };

public:
  vector<BestPlacement>    &getForward(uint32 fi)  { return(_pForward[fi]); };
  vector<BestReverse>      &getReverse(uint32 fi)  { return(_pReverse[fi]); };

public:
  void                      buildGraph(const char   *prefix,
                                       double        deviationRepeat,
                                       TigVector    &tigs);

private:
  vector<BestPlacement>  *_pForward;   //  Where each read is placed in other tigs
  vector<BestReverse>    *_pReverse;   //  What reads overlap to me
};




#endif  //  INCLUDE_AS_BAT_ASSEMBLYGRAPH
