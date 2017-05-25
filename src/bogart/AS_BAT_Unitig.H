
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
 *    src/AS_BAT/AS_BAT_Unitig.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_UNITIG
#define INCLUDE_AS_BAT_UNITIG

#include "AS_global.H"
#include "AS_BAT_TigVector.H"

#include "stddev.H"

#include <vector>
#include <algorithm>

class BestEdgeOverlap;



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







//  Derived from IntMultiPos, but removes some of the data (48b in IntMultiPos, 32b in struct
//  ufNode).  The minimum size (bit fields, assuming maximum limits, not using the contained
//  field) seems to be 24b, and is more effort than it is worth (just removing 'contained' would be
//  a chore).
//
//  ufNode is, of course, 'unitig fragment node'.
//
class ufNode {
public:
  uint32           ident;
  uint32           contained;
  uint32           parent;     //  IID of the read we align to

  int32            ahang;       //  If parent defined, these are relative
  int32            bhang;       //  that read

  SeqInterval      position;

  bool   isForward(void) const   {  return(position.isForward());  };
  bool   isReverse(void) const   {  return(position.isReverse());  };

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
};





class Unitig {
private:
  Unitig(TigVector *v) {
    _vector        = v;
    _length        = 0;
    _id            = 0;

    _isUnassembled = false;
    _isBubble      = false;
    _isRepeat      = false;
    _isCircular    = false;
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

  void cleanUp(void);

  // getNumRandomRead() is a placeholder, random reads should not
  // contain guides, or other reads that are not randomly sampled
  // across the whole genome.

  uint32 id(void)                 { return(_id); };      //  ID internal to bogart

  int32  getLength(void)          { return(_length);       };
  uint32 getNumReads(void)        { return(ufpath.size()); };
  uint32 getNumRandomReads(void)  { return(getNumReads()); };

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
    };

    double    max(double deviations) {
      return(dev.mean() + deviations * dev.stddev());
    };

    bool operator<(const epValue &that) const { return(bgn < that.bgn); };
    bool operator<(const uint32 &that)  const { return(bgn < that);     };



    uint32          bgn;
    uint32          end;

    stdDev<float>   dev;
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
  //  Classification.  The output is in three files: 'unassembled', 'bubbles', 'contigs' (defined as
  //  not unassembled and not bubble).

  uint32            _isUnassembled;  //  Is a single read or a pseudo singleton
  uint32            _isBubble;       //  Annotation: from a failed bubble pop
  uint32            _isRepeat;       //  Annotation: from an identified repeat region
  uint32            _isCircular;     //  Annotation: has overlap to self

  char   type(void) {
    if (_isUnassembled)  return('U');
    if (_isBubble)       return('B');
    if (_isRepeat)       return('R');
    if (_isCircular)     return('C');
    return('N');
  }
};


#endif  //  INCLUDE_AS_BAT_UNITIG
