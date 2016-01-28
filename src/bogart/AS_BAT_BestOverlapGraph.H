
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_BestOverlapGraph.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-29
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-21 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_BESTOVERLAPGRAPH
#define INCLUDE_AS_BAT_BESTOVERLAPGRAPH

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"

class BestOverlaps {
public:
  BestEdgeOverlap     _best5;
  BestEdgeOverlap     _best3;
  BestContainment     _bestC;
};


class BestScores {
public:
  BestScores() {
    _best5score = 0;
    _best3score = 0;
    _bestCscore = 0;
  };

  uint64  _best5score;
  uint64  _best3score;
  uint64  _bestCscore;
};



class BestOverlapGraph {
private:
  void   removeSuspicious(void);
  void   examineOnlyTopN(void);
  void   removeSpurs(void);
  void   removeFalseBest(void);
  void   removeWeak(double threshold);

public:
  BestOverlapGraph(double      erate,
                   const char *prefix,
                   double      doRemoveWeakThreshold,
                   bool        doRemoveSuspicious,
                   bool        doRemoveSpurs);

  BestOverlapGraph(double erate,
                   set<uint32> *restrict);

  ~BestOverlapGraph() {
    delete [] _bestA;
    delete [] _scorA;
  };

  //  Given a fragment UINT32 and which end, returns pointer to
  //  BestOverlap node.
  BestEdgeOverlap *getBestEdgeOverlap(uint32 fragid, bool threePrime) {
    if (_bestA)
      return((threePrime) ? (&_bestA[fragid]._best3) : (&_bestA[fragid]._best5));
    return((threePrime) ? (&_bestM[fragid]._best3) : (&_bestM[fragid]._best5));
  };

  // given a FragmentEnd sets it to the next FragmentEnd after following the
  // best edge
  FragmentEnd   followOverlap(FragmentEnd end) {
    if (end.fragId() == 0)
      return(FragmentEnd());

    BestEdgeOverlap *edge = getBestEdgeOverlap(end.fragId(), end.frag3p());

    return(FragmentEnd(edge->fragId(), !edge->frag3p()));
  };

  bool isContained(const uint32 fragid) {
    if (_bestA)
      return(_bestA[fragid]._bestC.isContained);
    return(_bestM[fragid]._bestC.isContained);
  };

  bool isSuspicious(const uint32 fragid) {
    return(_suspicious.count(fragid) > 0);
  };

  // Given a containee, returns pointer to BestContainment record
  BestContainment *getBestContainer(const uint32 fragid) {
    if (_bestA)
      return(&_bestA[fragid]._bestC);
    return(&_bestM[fragid]._bestC);
  };

  void      reportBestEdges(const char *prefix);

public:
  void      rebuildBestContainsWithoutSingletons(UnitigVector  &unitigs,
                                                 double         erate,
                                                 const char    *prefix);

private:
  bool     isOverlapBadQuality(const BAToverlap& olap);
  bool     isOverlapRestricted(const BAToverlap &olap);
  uint64   scoreOverlap(const BAToverlap& olap);

private:
  void     scoreContainment(const BAToverlap& olap);
  void     scoreEdge(const BAToverlap& olap);

private:
  uint64  &bestCscore(uint32 id) {
    if (_restrictEnabled == false)
      return(_scorA[id]._bestCscore);
    return(_scorM[id]._bestCscore);
  };

  uint64  &best5score(uint32 id) {
    if (_restrictEnabled == false)
      return(_scorA[id]._best5score);
    return(_scorM[id]._best5score);
  };

  uint64  &best3score(uint32 id) {
    if (_restrictEnabled == false)
      return(_scorA[id]._best3score);
    return(_scorM[id]._best3score);
  };

private:
  BestOverlaps              *_bestA;
  BestScores                *_scorA;

  set<uint32>                _suspicious;

  map<uint32, BestOverlaps>  _bestM;
  map<uint32, BestScores>    _scorM;

  set<uint32>               *_restrict;
  bool                       _restrictEnabled;

public:
  double                     _erate;
}; //BestOverlapGraph



#endif //INCLUDE_AS_BAT_BESTOVERLAPGRAPH
