
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

#ifndef INCLUDE_AS_BAT_UNITIG
#define INCLUDE_AS_BAT_UNITIG

#include "runtime.H"
#include "AS_BAT_TigVector.H"
#include "AS_BAT_OverlapCache.H"

#include "stddev.H"

#include <vector>
#include <set>
#include <algorithm>

using namespace std;


class BestEdgeOverlap;
class optPos;


class SeqInterval {
public:
  SeqInterval() {
    bgn = 0;
    end = 0;
  };
  ~SeqInterval() {
  };


  int32  min(void) const         {  return(::min(bgn, end));  };
  int32  max(void) const         {  return(::max(bgn, end));  };

  bool   isForward(void) const   {  return(bgn < end);        };
  bool   isReverse(void) const   {  return(bgn > end);        };

  bool   operator==(SeqInterval const that) const {
    return(((bgn == that.bgn) && (end == that.end)) ||
           ((bgn == that.end) && (end == that.bgn)));
  };

  bool   operator!=(SeqInterval const that) const {
    return(((bgn != that.bgn) || (end != that.end)) &&
           ((bgn != that.end) || (end != that.bgn)));
  };

  bool   operator<(SeqInterval const that) const {
    return(min() < that.min());
  };


public:
  int32  bgn;  //  MUST be signed!  Read placement needs to set coordinates to negative
  int32  end;  //  coordinates to indicate the read extends off the start of the tig.
};



//  True if A is contained in B.
inline
bool
isContained(int32 Abgn, int32 Aend,
            int32 Bbgn, int32 Bend) {
  assert(Abgn < Aend);
  assert(Bbgn < Bend);
  return((Bbgn <= Abgn) &&
         (Aend <= Bend));
}

inline
bool
isContained(SeqInterval &A, SeqInterval &B) {
  return((B.min() <= A.min()) &&
         (A.max() <= B.max()));
}



//  True if the A and B intervals overlap
inline
bool
isOverlapping(int32 Abgn, int32 Aend,
              int32 Bbgn, int32 Bend) {
  assert(Abgn < Aend);
  assert(Bbgn < Bend);
  return((Abgn < Bend) &&
         (Bbgn < Aend));
}

inline
bool
isOverlapping(SeqInterval &A, SeqInterval &B) {
  return((A.min() < B.max()) &&
         (B.min() < A.max()));
}



//  Derived from IntMultiPos, but removes some of the data (48b in IntMultiPos, 32b in struct
//  ufNode).  The minimum size (bit fields, assuming maximum limits, not using the contained
//  field) seems to be 24b, and is more effort than it is worth (just removing 'contained' would be
//  a chore).
//
//  ufNode is, of course, 'unitig fragment node'.
//
class ufNode {
public:
  ufNode(uint32 fi, SeqInterval pos) {
    ident        = fi;
    contained    = 0;
    parent       = 0;

    ahang        = 0;
    bhang        = 0;

    position     = pos;
  };
  ufNode(uint32 fi=0, uint32 bgn=0, uint32 end=0) {
    ident        = fi;
    contained    = 0;
    parent       = 0;

    ahang        = 0;
    bhang        = 0;

    position.bgn = bgn;
    position.end = end;
  };
  ~ufNode() {
  };

public:
  bool   isForward(void) const   {  return(position.isForward());  };
  bool   isReverse(void) const   {  return(position.isReverse());  };

  //  These two functions return the coordinates of an overlap on this read.
  //  The hangs are as they are from the overlap; do not do anything fancy to
  //  them.
  //              MinCoord   MaxCoord                    min       max
  //                  v         v                         v         v
  //           -----------A----->    bhang          "end" <-----A------------ "bgn"
  //           ahang  ----B---------------       ---------------B----   ahang
  //                                             bhang
  //         "end"            "bgn"
  //           <----------A------   -ahang  <-- minCoord isReverse(), add -(-bhang)
  //           -bhang  ---B---------------      or zero (picture above)
  //
  int32  hangToMinCoord(int32 ahang, int32 bhang) {
    if (position.isForward())
      return(position.bgn + max(0, ahang));
    else
      return(position.end - min(0, bhang));
  };

  int32  hangToMaxCoord(int32 ahang, int32 bhang) {
    if (position.isForward())
      return(position.end + min(0, bhang));
    else
      return(position.bgn - max(0, ahang));
  };

  bool  operator<(ufNode const &that) const {
    int32 abgn = (position.bgn < position.end) ? position.bgn : position.end;
    int32 aend = (position.bgn < position.end) ? position.end : position.bgn;

    int32 bbgn = (that.position.bgn < that.position.end) ? that.position.bgn : that.position.end;
    int32 bend = (that.position.bgn < that.position.end) ? that.position.end : that.position.bgn;

    if (abgn < bbgn) return(true);    //  A starts before B!
    if (abgn > bbgn) return(false);   //  B starts before A!

    if (aend < bend) return(false);   //  A contained in B, not less than.
    if (aend > bend) return(true);    //  B contained in A, is  less than.

    return(false);                    //  Equality, not less than.
  };

public:
  uint32           ident;
  uint32           contained;
  uint32           parent;     //  IID of the read we align to

  int32            ahang;       //  If parent defined, these are relative
  int32            bhang;       //  that read

  SeqInterval      position;
};



class Unitig {
private:
  Unitig(TigVector *v) {
    _vector        = v;
    _length        = 0;
    _id            = 0;

    _isUnassembled = false;
    _isRepeat      = false;
    _isCircular    = false;
    _isBubble      = false;
  };

public:
  ~Unitig(void) {
  };

  friend class TigVector;

  void sort(void) {
    std::sort(ufpath.begin(), ufpath.end());

    for (uint32 fi=0; fi<ufpath.size(); fi++)
      _vector->registerRead(ufpath[fi].ident, _id, fi);
  };
  //void   bubbleSortLastRead(void);
  void reverseComplement(bool doSort=true);

  //  Ensure that the children are sorted by begin position,
  //  and that unitigs start at position zero.
  void cleanUp(void);

  //  Recompute bgn/end positions using all overlaps.
  bool optimize_isCompatible(uint32       ii,
                             uint32       jj,
                             BAToverlap  &olap,
                             bool         inInit,
                             bool         secondPass,
                             bool         beVerbose);

  void optimize_initPlace(uint32        pp,
                          optPos       *op,
                          optPos       *np,
                          bool          firstPass,
                          set<uint32>  &failed,
                          bool          beVerbose);
  void optimize_recompute(uint32        ii,
                          optPos       *op,
                          optPos       *np,
                          bool          beVerbose);
  void optimize_expand(optPos       *op, bool beVerbose);
  void optimize_setPositions(optPos       *op,
                             bool          beVerbose);
  void optimize(const char *prefix, const char *label);


  uint32 id(void)                 { return(_id); };

  int32  getLength(void)          { return(_length);       };
  uint32 getNumReads(void)        { return(ufpath.size()); };

  //  Place 'read' using an edge to some read in this tig.  The edge is from 'read3p' end.
  //
  bool   placeRead(ufNode          &read,     //  resulting placement
                   uint32           readId,   //  read we want to place
                   bool             read3p,   //  end that the edge is from
                   BestEdgeOverlap *edge);    //  edge to something in this tig

  void   addRead(ufNode node, int offset=0, bool report=false);


public:
  class epValue {
  public:
    epValue(uint32 b, uint32 e) {
      bgn    = b;
      end    = e;
      mean   = 0;
      stddev = 0;
    };

    epValue(uint32 b, uint32 e, float m, float s) {
      bgn    = b;
      end    = e;
      mean   = m;
      stddev = s;
    };

    double    max(double deviations) {
      double val = mean + deviations * stddev;
      return (val < 1e-10 ? 1e-10 : val);
    };

    bool operator<(const epValue &that) const { return(bgn < that.bgn); };
    bool operator<(const uint32 &that)  const { return(bgn < that);     };

    uint32          bgn;
    uint32          end;

    float           mean;
    float           stddev;
  };

  static size_t epValueSize(void) { return(sizeof(epValue)); };

  void   computeArrivalRate(const char *prefix,
                            const char *label,
                            vector<int32> *hist);

  void   computeErrorProfile(const char *prefix, const char *label);
  void   reportErrorProfile(const char *prefix, const char *label);
  void   clearErrorProfile(void)       { errorProfile.clear(); };

  double overlapConsistentWithTig(double deviations,
                                  uint32 bgn, uint32 end,
                                  double erate);


  //  Returns the read that is touching the start of the tig.
  ufNode *firstRead(void) {
    ufNode  *rd5 = &ufpath.front();

    for (uint32 fi=1; (fi < ufpath.size()) && (rd5->position.min() != 0); fi++)
      rd5 = &ufpath[fi];

    if (rd5->position.min() != 0)
      fprintf(stderr, "ERROR: firstRead() in tig %u doesn't start at the start\n", id());
    assert(rd5->position.min() == 0);

    return(rd5);
  };


  //  Returns the read that is touching the end of the tig.
  ufNode *lastRead(void) {
    ufNode  *rd3 = &ufpath.back();

    for (uint32 fi=ufpath.size()-1; (fi-- > 0) && (rd3->position.max() != getLength()); )
      rd3 = &ufpath[fi];

    if (rd3->position.max() != getLength())
      fprintf(stderr, "ERROR: lastRead() in tig %u doesn't end at the end\n", id());
    assert(rd3->position.max() == getLength());

    return(rd3);
  };

  // Public Member Variables
public:
  vector<ufNode>     ufpath;
  vector<epValue>    errorProfile;
  vector<uint32>     errorProfileIndex;

public:
  //  r > 0 guards against calling these from Idx's, while r < size guards
  //  against calling with Id's.
  //
  uint32   inUnitig(uint32 r)     { assert(r > 0);             return(_vector->inUnitig(r));   };
  uint32   ufpathIdx(uint32 r)    { assert(r > 0);             return(_vector->ufpathIdx(r));  };

  ufNode  *readFromId(uint32 r)   { assert(r > 0);             return(&ufpath[ ufpathIdx(r) ]);  };
  ufNode  *readFromIdx(uint32 r)  { assert(r < ufpath.size()); return(&ufpath[ r ]);             };

private:
  TigVector        *_vector;   //  For updating the read map.

private:
  int32             _length;
  uint32            _id;

public:
  //  Classification.

  bool              _isUnassembled;  //  Is a single read or a pseudo singleton.
  bool              _isRepeat;       //  Is from an identified repeat region.
  bool              _isCircular;     //  Is (probably) a circular tig.
  bool              _isBubble;       //  Is (probably) a bubble.
};


#endif  //  INCLUDE_AS_BAT_UNITIG
