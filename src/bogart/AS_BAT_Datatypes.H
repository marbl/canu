
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
 *    src/AS_BAT/AS_BAT_Datatypes.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-28
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-11 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_DATATYPES
#define INCLUDE_AS_BAT_DATATYPES

#include "AS_global.H"
#include "ovStore.H"
#include "gkStore.H"

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

#ifndef BROKEN_CLANG_OpenMP
#include <omp.h>
#endif

using namespace std;

#include "AS_BAT_Logging.H"
#include "AS_BAT_OverlapCache.H"

const uint32  noUnitig = 0xffffffff;

////////////////////////////////////////

class FragmentInfo;
class OverlapCache;
class BestOverlapGraph;
class ChunkGraph;

extern FragmentInfo     *FI;
extern OverlapCache     *OC;
extern BestOverlapGraph *OG;
extern ChunkGraph       *CG;

////////////////////////////////////////

////////////////////////////////////////

//  These MUST be signed integers.  Read placement -- Unitig::placeFrag() for sure -- depends on
//  being able to set coordinates to negative values to indicate they extend the unitig off the
//  beginning.
//
struct SeqInterval {
  int32           bgn;
  int32           end;
};

static const SeqInterval NULL_SEQ_LOC = {0,0};

inline
bool
isReverse(SeqInterval pos) {
  return(pos.bgn > pos.end);
}

inline
bool
operator==(SeqInterval a, SeqInterval b) {
  return(((a.bgn == b.bgn) && (a.end == b.end)) ||
         ((a.bgn == b.end) && (a.end == b.bgn)));
}

inline
bool
operator!=(SeqInterval a, SeqInterval b) {
  return(((a.bgn != b.bgn) || (a.end != b.end)) &&
         ((a.bgn != b.end) || (a.end != b.bgn)));
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


//  Swiped from AS_OVS_overlap.h, modified to take a BAToverlap instead of an OVSoverlap.

static
uint32
AS_BAT_overlapAEndIs5prime(const BAToverlap& olap) {
  return((olap.a_hang < 0) && (olap.b_hang < 0));
}

static
uint32
AS_BAT_overlapAEndIs3prime(const BAToverlap& olap) {
  return((olap.a_hang > 0) && (olap.b_hang > 0));
}

static
uint32
AS_BAT_overlapBEndIs3prime(const BAToverlap& olap) {
  return((AS_BAT_overlapAEndIs5prime(olap) && (olap.flipped == false)) ||
         (AS_BAT_overlapAEndIs3prime(olap) && (olap.flipped == true)));
}


class BestEdgeOverlap {
public:
  BestEdgeOverlap() {
    clear();
  };
  ~BestEdgeOverlap() {
  };

  void    clear(void) {
    _id    = 0;
    _e3p   = 0;
    _ahang = 0;
    _bhang = 0;
  };

  void    set(BAToverlap const &olap) {
    _id    = olap.b_iid;
    _e3p   = AS_BAT_overlapBEndIs3prime(olap);
    _ahang = olap.a_hang;
    _bhang = olap.b_hang;
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
  uint32            _id;
  uint64            _e3p   : 1;    //  Overlap with the 3' end of that fragment
  int64             _ahang : AS_MAX_READLEN_BITS+1;
  int64             _bhang : AS_MAX_READLEN_BITS+1;
};


// Contains what kind of containment relationship exists between fragment a and fragment b
//
class BestContainment{
public:
  BestContainment() {
    clear();
  };
  ~BestContainment() {
  };

  void    clear(void) {
    container       = 0;
    isContained     = false;
    sameOrientation = false;
    a_hang          = 0;
    b_hang          = 0;
  };

  uint32  container;
  uint64  isContained     : 1;
  uint64  sameOrientation : 1;
  int64   a_hang          : AS_MAX_READLEN_BITS+1;
  int64   b_hang          : AS_MAX_READLEN_BITS+1;
};




class FragmentInfo {
public:
  FragmentInfo(gkStore *gkp, const char *prefix, uint32 minReadLen);
  ~FragmentInfo();

  uint64  memoryUsage(void) {
    return((3 * sizeof(uint32) * _numFragments) +
           (2 * sizeof(double) * _numLibraries) +
           (2 * sizeof(uint32) * _numLibraries));
  };

  uint32  numFragments(void) { return(_numFragments); };
  uint32  numLibraries(void) { return(_numLibraries); };

  uint32  fragmentLength(uint32 iid) { return(_fragLength[iid]); };
  uint32  libraryIID(uint32 iid)     { return(_libIID[iid]);  };

  uint32  overlapLength(uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang) {
    int32  alen = fragmentLength(a_iid);
    int32  blen = fragmentLength(b_iid);
    int32  aovl = 0;
    int32  bovl = 0;

    assert(alen > 0);
    assert(blen > 0);

    if (a_hang < 0) {
      //  b_hang < 0      ?     ----------  :     ----
      //                  ?  ----------     :  ----------
      //
      aovl = (b_hang < 0) ? (alen + b_hang) : (alen);
      bovl = (b_hang < 0) ? (blen + a_hang) : (blen + a_hang - b_hang);
    } else {
      //  b_hang < 0      ?  ----------              :  ----------
      //                  ?     ----                 :     ----------
      //
      aovl = (b_hang < 0) ? (alen - a_hang + b_hang) : (alen - a_hang);
      bovl = (b_hang < 0) ? (blen)                   : (blen - b_hang);
    }

    if ((aovl <= 0) || (bovl <= 0) || (aovl > alen) || (bovl > blen)) {
      fprintf(stderr, "WARNING: bogus overlap found for A="F_U32" B="F_U32"\n", a_iid, b_iid);
      fprintf(stderr, "WARNING:                     A len="F_S32" hang="F_S32" ovl="F_S32"\n", alen, a_hang, aovl);
      fprintf(stderr, "WARNING:                     B len="F_S32" hang="F_S32" ovl="F_S32"\n", blen, b_hang, bovl);
    }

    if (aovl < 0)     aovl = 0;
    if (bovl < 0)     bovl = 0;

    if (aovl > alen)  aovl = alen;
    if (bovl > blen)  bovl = blen;

    assert(aovl > 0);
    assert(bovl > 0);
    assert(aovl <= alen);
    assert(bovl <= blen);

    //  AVE does not work.      return((uint32)((aovl, bovl)/2));
    //  MAX does not work.      return((uint32)MAX(aovl, bovl));

    return(aovl);
  };

  //  For DIAGNOSTIC ONLY.  Mark this fragment as 'deleted' so checkUnitigMembership() will pass
  //  when we ignore fragments.  The resuting assembly will fail in CGW.
  void    markAsIgnore(uint32 iid) {
    _fragLength[iid] = 0;
  };

private:
  void      save(const char *prefix);
  bool      load(const char *prefix);

  uint32   _numFragments;
  uint32   _numLibraries;

  uint32  *_fragLength;
  uint32  *_libIID;

  uint32  *_numFragsInLib;
};

#endif



