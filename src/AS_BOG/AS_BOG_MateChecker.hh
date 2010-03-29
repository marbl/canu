
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

#ifndef INCLUDE_AS_BOG_MATECHEKER
#define INCLUDE_AS_BOG_MATECHEKER

static const char *rcsid_INCLUDE_AS_BOG_MATECHEKER = "$Id: AS_BOG_MateChecker.hh,v 1.35 2010-03-29 16:52:40 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"

#if 0
typedef std::map<uint32,uint32> IdMap;

typedef IdMap::iterator IdMapIter;
typedef IdMap::const_iterator IdMapConstIter;

typedef std::vector<int> DistanceList;
typedef DistanceList::const_iterator DistanceListCIter;

typedef std::map<uint32,DistanceList> LibraryDistances;
typedef LibraryDistances::const_iterator LibDistsConstIter;
#endif

static const SeqInterval NULL_SEQ_LOC = {0,0};



struct DistanceCompute {
  double  stddev;
  double  mean;
  uint32  samples;

  uint32  distancesLen;
  uint32  distancesMax;
  uint32 *distances;

  DistanceCompute() {
    stddev       = 0.0;
    mean         = 0.0;
    samples      = 0;

    distancesLen = 0;
    distancesMax = 1048576;
    distances    = new uint32 [distancesMax];
  };

  ~DistanceCompute() {
    delete [] distances;
  }
};




struct MateChecker{
  MateChecker(FragmentInfo *fi) {
    _globalStats = NULL;
    _fi          = fi;
  };
  ~MateChecker() {
    delete [] _globalStats;
  };

  void checkUnitigGraph(UnitigGraph &tigGraph, int badMateBreakThreshold);

private:
  void  accumulateLibraryStats(Unitig *utg);
  void  computeGlobalLibStats(UnitigGraph &tigGraph);

  UnitigBreakPoints* computeMateCoverage(Unitig *utg,
                                         BestOverlapGraph *bog,
                                         int badMateBreakThreshold);

  void evaluateMates(UnitigGraph &tigGraph);

  void moveContains(UnitigGraph &tigGraph);
  void splitDiscontinuousUnitigs(UnitigGraph &tigGraph);

private:
  DistanceCompute  *_globalStats;
  FragmentInfo     *_fi;
};



#if 1
inline
bool
operator==(SeqInterval a, SeqInterval b) {
  return((a.bgn == b.bgn) && (a.end == b.end) ||
         (a.bgn == b.end) && (a.end == b.bgn));
};

inline
bool
operator<(SeqInterval a, SeqInterval b) {
  if (isReverse(a)) {
    if (isReverse(b)) return a.end < b.end;
    else              return a.end < b.bgn;
  } else {
    if (isReverse(b)) return a.bgn < b.end;
    else              return a.bgn < b.bgn;
  }
};
#endif

class MateLocationEntry {
public:
  SeqInterval mlePos1;
  SeqInterval mlePos2;
  uint32      mleFrgID1;
  uint32      mleFrgID2;
  uint32      mleUtgID1;
  uint32      mleUtgID2; // in the future the table might be across unitigs
  bool        isGrumpy;

  //bool operator==(MateLocation &that) {
  //  return((mlePos1 == that.mlePos1) && (mlePos2 == that.mlePos2));
  //};

  bool operator<(MateLocationEntry const &that) const {

#if 0
    if (mlePos1  < that.mlePos1)                           return true;
    if (mlePos1 == that.mlePos1 && mlePos2 < that.mlePos2) return true;
    else                                                   return false;
#endif

    int32  m1 = MIN(     mlePos1.bgn,      mlePos1.end);
    int32  t1 = MIN(that.mlePos1.bgn, that.mlePos1.end);
    int32  m2 = MIN(     mlePos2.bgn,      mlePos2.end);
    int32  t2 = MIN(that.mlePos2.bgn, that.mlePos2.end);

    return((m1 < t1) || ((m1 == t1) && (m2 < t2)));
  };
};



//  The MateLocation table builds a table of positions of mated reads.
//    o  Unmated reads are NOT in the table.
//    o  Mates in other unitigs are not in the table.  The fragment
//       in this unitig is present, but the mate is NULL.
//    o  Mates in the same unitig are in the table.
//
class MateLocation {
public:
  MateLocation(FragmentInfo *fi, Unitig *utg, DistanceCompute *dc) {
    MateLocationEntry   mle;

    mle.mlePos1.bgn = mle.mlePos1.end = 0;
    mle.mlePos2.bgn = mle.mlePos2.end = 0;

    mle.mleFrgID1 = 0;
    mle.mleFrgID2 = 0;

    mle.mleUtgID1 = 0;
    mle.mleUtgID2 = 0;

    mle.isGrumpy = false;

    _table.clear();
    _table.push_back(mle);

    _tigLen = utg->getLength();

    goodGraph   = new int32 [_tigLen + 1];
    badFwdGraph = new int32 [_tigLen + 1];
    badRevGraph = new int32 [_tigLen + 1];

    memset(goodGraph,   0, sizeof(int32) * (_tigLen + 1));
    memset(badFwdGraph, 0, sizeof(int32) * (_tigLen + 1));
    memset(badRevGraph, 0, sizeof(int32) * (_tigLen + 1));

    _fi = fi;

    buildTable(utg);
    buildHappinessGraphs(utg, dc);
  };

  ~MateLocation() {
    delete [] goodGraph;
    delete [] badFwdGraph;
    delete [] badRevGraph;
  };
            
  MateLocationEntry getById(uint32 fragId) {
    map<uint32,uint32>::const_iterator  e = _iidToTableEntry.find(fragId);

    if (e == _iidToTableEntry.end())
      return _table[0];
    else
      return _table[e->second];
  };

  int32  *goodGraph;
  int32  *badFwdGraph;
  int32  *badRevGraph;

private:
  void buildTable(Unitig *utg);
  void buildHappinessGraphs(Unitig *utg, DistanceCompute *);

  void incrRange(int32 *graph, int32 val, int32 n, int32 m) {
    n = MAX(n, 0);
    m = MIN(m, _tigLen);

    //  Earlier versions asserted n<m (and even earlier versions used i<=m in the loop below, which
    //  made this far more complicated than necessary).  Now, we don't care.  We'll adjust n and m
    //  to the min/max possible, and ignore out of bounds cases.  Those happen when, for example,
    //  fragments are the same orientation.  If one of those is the last fragment in the unitig,
    //  we'll call incrRange with n=(the higher coord)=(_tigLen), and m=(the lower coord + max
    //  insert size).  We threshold m to _tigLen, and correctly do nothing in the loop.

    for(uint32 i=n; i<m; i++)
      graph[i] += val;
  };

  uint32                     _tigLen;

  vector<MateLocationEntry>  _table;
  map<uint32,uint32>         _iidToTableEntry;
  FragmentInfo              *_fi;
};


#endif
