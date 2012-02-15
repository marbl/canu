
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef INCLUDE_AS_BAT_BESTOVERLAPGRAPH
#define INCLUDE_AS_BAT_BESTOVERLAPGRAPH

static const char *rcsid_INCLUDE_AS_BAT_BESTOVERLAPGRAPH = "$Id: AS_BAT_BestOverlapGraph.H,v 1.6 2012-02-15 03:41:08 brianwalenz Exp $";

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
public:
  BestOverlapGraph(double erate,
                   double elimit,
                   const char *prefix);

  BestOverlapGraph(double erate,
                   double elimit,
                   set<AS_IID> *restrict);

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

  void      reportBestEdges(void);

private:
  bool     isOverlapBadQuality(const BAToverlap& olap);
  bool     isOverlapRestricted(const BAToverlap &olap);
  uint64   scoreOverlap(const BAToverlap& olap);

private:
  void     scoreContainment(const BAToverlap& olap);
  void     scoreEdge(const BAToverlap& olap);

private:
  uint64  &bestCscore(AS_IID id) {
    if (_scorA)
      return(_scorA[id]._bestCscore);
    return(_scorM[id]._bestCscore);
  };

  uint64  &best5score(AS_IID id) {
    if (_scorA)
      return(_scorA[id]._best5score);
    return(_scorM[id]._best5score);
  };

  uint64  &best3score(AS_IID id) {
    if (_scorA)
      return(_scorA[id]._best3score);
    return(_scorM[id]._best3score);
  };

private:
  BestOverlaps              *_bestA;
  BestScores                *_scorA;

  set<AS_IID>                _suspicious;

  map<AS_IID, BestOverlaps>  _bestM;
  map<AS_IID, BestScores>    _scorM;

  set<AS_IID>               *_restrict;

public:
  double                      mismatchCutoff;
  double                      mismatchLimit;
}; //BestOverlapGraph



#endif //INCLUDE_AS_BAT_BESTOVERLAPGRAPH
