
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

#ifndef INCLUDE_AS_BAT_DATATYPES
#define INCLUDE_AS_BAT_DATATYPES

static const char *rcsid_INCLUDE_AS_BAT_DATATYPES = "$Id: AS_BAT_Datatypes.H,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include <map>
#include <set>
#include <list>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

const uint32  noUnitig = 0xffffffff;

////////////////////////////////////////

#define BADMATE_INTRA_STDDEV 3  //  Mates more than this stddev away in the same unitig are bad
#define BADMATE_INTER_STDDEV 5  //  Mates more than this stddev away from the end of the unitig are bad

//  Historical.  Eli had problems with STL and max.
#undef max

////////////////////////////////////////

class FragmentInfo;
class BestOverlapGraph;
class ChunkGraph;
class InsertSizes;

extern FragmentInfo     *FI;
extern BestOverlapGraph *OG;
extern ChunkGraph       *CG;
extern InsertSizes      *IS;

////////////////////////////////////////

void  setLogFile(const char *prefix, const char *name);

#define logFileFlagSet(L) ((logFileFlags & L) == L)

extern FILE   *logFile;
extern uint32  logFileOrder;
extern uint64  logFileFlags;

extern uint64 LOG_OVERLAP_QUALITY;
extern uint64 LOG_OVERLAPS_USED;
extern uint64 LOG_CHUNK_GRAPH;
extern uint64 LOG_INTERSECTIONS;
extern uint64 LOG_POPULATE_UNITIG;
extern uint64 LOG_INTERSECTION_BREAKING;
extern uint64 LOG_INTERSECTION_BUBBLES;
extern uint64 LOG_INTERSECTION_BUBBLES_DEBUG;
extern uint64 LOG_INTERSECTION_JOINING;
extern uint64 LOG_INTERSECTION_JOINING_DEBUG;
extern uint64 LOG_INITIAL_CONTAINED_PLACEMENT;
extern uint64 LOG_HAPPINESS;
extern uint64 LOG_INTERMEDIATE_UNITIGS;
extern uint64 LOG_MATE_SPLIT_ANALYSIS;
extern uint64 LOG_MATE_SPLIT_DISCONTINUOUS;
extern uint64 LOG_MATE_SPLIT_UNHAPPY_CONTAINS;
extern uint64 LOG_MATE_SPLIT_COVERAGE_PLOT;
extern uint64 LOG_STDERR;

extern uint64 LOG_PLACE_FRAG;

////////////////////////////////////////

static const SeqInterval NULL_SEQ_LOC = {0,0};

inline
bool
isReverse(SeqInterval pos) {
  return(pos.bgn > pos.end);
}

inline
bool
operator==(SeqInterval a, SeqInterval b) {
  return((a.bgn == b.bgn) && (a.end == b.end) ||
         (a.bgn == b.end) && (a.end == b.bgn));
}

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
}

////////////////////////////////////////

class FragmentEnd {
public:
  FragmentEnd() {
    _id  = 0;
    _e3p = false;
  };
  FragmentEnd(uint32 id, bool e3p) {
    _id  = id;
    _e3p = e3p;
  };

  uint32  fragId(void)  const { return(_id); };
  bool    frag3p(void)  const { return(_e3p == true);  };
  bool    frag5p(void)  const { return(_e3p == false); };
  uint32  index(void)   const { return(_id * 2 + _e3p); };

  bool operator==(FragmentEnd const that) const {
    return((fragId() == that.fragId()) && (frag3p() == that.frag3p()));
  };

  bool operator!=(FragmentEnd const that) const {
    return((fragId() != that.fragId()) || (frag3p() != that.frag3p()));
  };

  bool operator<(FragmentEnd const that) const {
    if (fragId() != that.fragId())
      return fragId() < that.fragId();
    else
      return frag3p() < that.frag3p();
  };

private:
  uint32   _id:31;
  uint32   _e3p:1;
};


class BestEdgeOverlap {
public:
  BestEdgeOverlap() {
    _id    = 0;
    _e3p   = 0;
    _ahang = 0;
    _bhang = 0;
  };
  ~BestEdgeOverlap() {
  };

  void    set(OVSoverlap const &olap) {
    _id    = olap.b_iid;
    _e3p   = AS_OVS_overlapBEndIs3prime(olap);
    _ahang = olap.dat.ovl.a_hang;
    _bhang = olap.dat.ovl.b_hang;
  };

  void    set(uint32 id, bool e3p, int32 ahang, int32 bhang) {
    _id    = id;
    _e3p   = e3p;
    _ahang = ahang;
    _bhang = bhang;
  };


  uint32  fragId(void)  const { return(_id); };
  bool    frag3p(void)  const { return(_e3p == true);  };
  bool    frag5p(void)  const { return(_e3p == false); };

  int32   ahang(void)   const { return(_ahang); };
  int32   bhang(void)   const { return(_bhang); };

private:
  uint32            _id:31;
  uint32            _e3p:1;    //  Overlap with the 3' end of that fragment
  int32             _ahang;
  int32             _bhang;
};


// Contains what kind of containment relationship exists between
// fragment a and fragment b
//
// This needs 93 bits (plus a 64 bit pointer) and must be 8-byte
// aligned, so the best we can hope for is 24 bytes.
//
class BestContainment{
public:
  BestContainment() {
  };
  ~BestContainment() {
    delete [] olaps;
  };

  uint32  container:31;
  uint32  isContained:1;  //  See comments in AS_BAT_Unitig_PlaceFragmentsUsingEdges.cc

  int32   a_hang;
  int32   b_hang;

  uint32  sameOrientation:1;

  uint32  olapsSorted:1;
  uint32  olapsLen:29;
  uint32 *olaps;
};




class FragmentInfo {
public:
  FragmentInfo(gkStore *gkpStore, const char *prefix);
  ~FragmentInfo();

  uint32  numFragments(void) { return(_numFragments); };
  uint32  numLibraries(void) { return(_numLibraries); };

  uint32  fragmentLength(uint32 iid) { return(_fragLength[iid]); };
  uint32  mateIID(uint32 iid)        { return(_mateIID[iid]); };
  uint32  libraryIID(uint32 iid)     { return(_libIID[iid]);  };

  double  mean(uint32 iid)   { return(_mean[iid]); };
  double  stddev(uint32 iid) { return(_stddev[iid]); };

  uint32  numMatesInLib(uint32 iid) { return(_numMatesInLib[iid]); };

private:
  void      save(const char *prefix);
  bool      load(const char *prefix);

  uint32   _numFragments;
  uint32   _numLibraries;

  uint32  *_fragLength;
  uint32  *_mateIID;
  uint32  *_libIID;

  double  *_mean;
  double  *_stddev;

  uint32  *_numFragsInLib;
  uint32  *_numMatesInLib;
};


#endif



