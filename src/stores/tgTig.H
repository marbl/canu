
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
 *    Brian P. Walenz from 2014-DEC-22 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef TG_TIG_H
#define TG_TIG_H

#include "AS_global.H"

#include "gkStore.H"


typedef enum {
  tgTig_noclass      = 0x00,    //  Could use a tgTig_read for corrections
  tgTig_unassembled  = 0x01,
  tgTig_bubble       = 0x02,
  tgTig_contig       = 0x03
} tgTig_class;

static
const
char *
toString(tgTig_class c) {
  switch (c) {
    case tgTig_noclass:      return("unsetc");  break;
    case tgTig_unassembled:  return("unassm");  break;
    case tgTig_bubble:       return("bubble");  break;
    case tgTig_contig:       return("contig");  break;
  }

  return("undefined-class");
}

static
tgTig_class
toClass(char const *c) {
  if (strcmp(c, "unsetc") == 0)
    return(tgTig_noclass);
  if (strcmp(c, "unassm") == 0)
    return(tgTig_unassembled);
  if (strcmp(c, "bubble") == 0)
    return(tgTig_bubble);
  if (strcmp(c, "contig") == 0)
    return(tgTig_contig);

  fprintf(stderr, "WARNING:  tiTig_class toClass('%s') is not a valid class.\n", c);

  return(tgTig_noclass);
}




//  Info about the placement of an object in a tig.  For unitigs, this
//  will be just reads.  For contigs, this could be unitigs and reads.
//
//  These coordinates are ALWAYS gapped coordinates.
//
class tgPosition {
public:
  tgPosition();

  //  Accessors

  uint32               ident(void)     const { return(_objID); };

  bool                 isRead(void)    const { return(_isRead   == true); };
  bool                 isUnitig(void)  const { return(_isUnitig == true); };
  bool                 isContig(void)  const { return(_isContig == true); };

  uint32               anchor(void)          { return(_anchor); };
  int32                aHang(void)           { return(_ahang);  };
  int32                bHang(void)           { return(_bhang);  };

  bool                 isForward(void) const { return(_isReverse == false); };
  bool                 isReverse(void) const { return(_isReverse == true);  };

  //  Position in the parent, both oriented (bgn/end) and unoriented (min/max).

  int32                bgn(void) const       { return((_isReverse == false) ? _min : _max); };
  int32                end(void) const       { return((_isReverse == false) ? _max : _min); };

  int32                min(void) const       { return(_min); };
  int32                max(void) const       { return(_max); };

  //  Amount of this object to ignore; e.g., trim from the ends.

  int32                askip(void) const       { return(_askip); };
  int32                bskip(void) const       { return(_bskip); };

  //  Delta values for the alignment to the parent.

  uint32               deltaOffset(void)     { return(_deltaOffset); };
  uint32               deltaLength(void)     { return(_deltaLen);    };

  //  Set just the anchor and hangs, leave positions alone.
  void                 setAnchor(uint32 anchor, int32 ahang, int32 bhang) {
    _anchor = anchor;
    _ahang  = ahang;
    _bhang  = bhang;
  };

  //  Set everything.  This is to be used by unitigger.
  void                 set(uint32 id, uint32 anchor, int32 ahang, int32 bhang, int32  bgn, int32  end) {

    _objID     = id;

    _isRead    = true;
    _isUnitig  = false;
    _isContig  = false;

    _anchor    = anchor;
    _ahang     = ahang;
    _bhang     = bhang;

    _askip     = 0;
    _bskip     = 0;

    if (bgn < end) {
      _min       = bgn;
      _max       = end;
      _isReverse = false;
    } else {
      _min       = end;
      _max       = bgn;
      _isReverse = true;
    }

    _deltaOffset = 0;
    _deltaLen    = 0;
  };

  //  Set the coords, ignoring orientation.
  void                 setMinMax(int32 min, int32 max) {
    _min = min;
    _max = max;
  };

  //  Operators

  bool                 operator<(tgPosition const &that) const {
    int32 al = min();
    int32 bl = that.min();

    if (al != bl)
      return(al < bl);

    if (that._anchor == _objID)  //  I'm his anchor, I must be before it.
      return(true);

    if (_anchor == that._objID)  //  My anchor is the other tig; I must
      return(false);             //  be after it.

    int32 ah = max();
    int32 bh = that.max();

    return(ah > bh);
  }

private:
public:
  uint32              _objID;            //  ID of this object

  uint32              _isRead      : 1;  //  A full length read alignment
  uint32              _isUnitig    : 1;  //
  uint32              _isContig    : 1;  //

  uint32              _isReverse   : 1;  //  Child is oriented forward relative to parent, used during consensus.

  uint32              _spare       : 28;

  uint32              _anchor;           //  ID of the like-type object we align to
  int32               _ahang;            //  Placement relative to anchor object
  int32               _bhang;            //

  int32               _askip;            //  Amount of sequence to not align on each end
  int32               _bskip;            //

  //  Must be signed, utgcns can push reads negative.
  //int32               _bgn;              //  Coords in the parent object we are part of
  //int32               _end;              //  (for a read, this will be the position in the unitig)

  //  Must be signed, utgcns can push reads negative.
  int32               _min;
  int32               _max;

  uint32              _deltaOffset;      //  Our delta alignment to the parent object consensus
  uint32              _deltaLen;
};


class tgTig;  //  Early declaration, for use in tgTigRecord operator=

//  On-disk tig, same as tgTig without the pointers
class tgTigRecord {
public:
  tgTigRecord();
  tgTigRecord(tgTig &tg) { *this = tg; };   //  to support tgTigRecord tr = tgtig

  tgTigRecord         &operator=(tgTig & tg);

private:
public:
  uint32              _tigID;

  double              _coverageStat;
  double              _microhetProb;   //  placeholder for some future multialignment based uniqueness probability

  tgTig_class         _class           : 2;    //  Output classification: unassembled, bubble, contig
  uint32              _suggestRepeat   : 1;    //  Bogart made this from detected repeat.
  uint32              _suggestCircular : 1;    //  Bogart found overlaps making a circle.

  uint32              _spare           : 32 - 2 - 2;

  uint32              _layoutLen;       //  Max coord
  uint32              _gappedLen;       //  Gapped consensus
  uint32              _childrenLen;
  uint32              _childDeltasLen;
};


//  Former MultiAlignT
//  In core data
class tgTig {
public:
  tgTig();    //  All data unallocated, lengths set to zero
  ~tgTig();   //  Releases memory

  //  Accessors

  uint32               tigID(void)                         { return(_tigID); };

  char const          *coordinateType(bool useGapped=true) { return( (consensusExists() == false) ? ("layout")       : ((useGapped == true) ? "gapped"       : "ungapped")       ); };
  uint32               length(bool useGapped=true)         { return( (consensusExists() == false) ? (layoutLength()) : ((useGapped == true) ? gappedLength() : ungappedLength()) ); };

  bool                 consensusExists(void)               { return(_gappedLen > 0); };

  char                *bases(bool useGapped=true)          { return( (useGapped == true) ? gappedBases() : ungappedBases() ); };
  char                *quals(bool useGapped=true)          { return( (useGapped == true) ? gappedQuals() : ungappedQuals() ); };

private:
  uint32               layoutLength(void)                  { return(_layoutLen); };

  uint32               gappedLength(void)                  { return(_gappedLen);   };
  char                *gappedBases(void)                   { return(_gappedBases); };
  char                *gappedQuals(void)                   { return(_gappedQuals); };

  void                 buildUngapped(void);

  uint32               ungappedLength(void)                { buildUngapped();  return(_ungappedLen);   };
  char                *ungappedBases(void)                 { buildUngapped();  return(_ungappedBases); };
  char                *ungappedQuals(void)                 { buildUngapped();  return(_ungappedQuals); };

public:

  //  This function needs to be hidden.  Coordinates in reads need to be transparent.
  uint32               mapGappedToUngapped(uint32 p)       {
    if (consensusExists() == false)
      return(p);

    buildUngapped();

    return(_gappedToUngapped[p]);
  };

  uint32               numberOfChildren(void)              {                            return(_childrenLen); };
  tgPosition          *getChild(uint32 c)                  { assert(c < _childrenLen);  return(_children + c);  };
  tgPosition          *addChild(void)                      {                            return(_children + _childrenLen++); };

  //  Operators

  void                 clear(void);    //  Clears data but doesn't release memory.

  tgTig               &operator=(tgTigRecord & tr);
  tgTig               &operator=(tgTig &tg);

  bool                 loadFromStreamOrLayout(FILE *F);

  void                 saveToStream(FILE *F);
  bool                 loadFromStream(FILE *F);

  void                 dumpLayout(FILE *F);
  bool                 loadLayout(FILE *F);

  void                 dumpFASTA(FILE *F, bool useGapped);
  void                 dumpFASTQ(FILE *F, bool useGapped);

  //  There are two multiAlign displays; this one, and one in abMultiAlign.
  void                 display(FILE     *F,
                               gkStore  *gkp,
                               uint32    displayWidth    = 100,    //  Width of display
                               uint32    displaySpacing  = 3,      //  Space between reads on the same line
                               bool      withQV          = false,
                               bool      withDots        = false);



private:
public:
  uint32              _tigID;           //  ID in the store, or UINT32_MAX if not set

  double              coverageStat(void)  { return(_coverageStat); };
  //double              microhetProb(void)  { return(_microhetProb); };

  double              _coverageStat;
  double              _microhetProb;  //  Placeholder.

  //  A variety of flags to suggest what type of unitig this is

  tgTig_class         _class           : 2;    //  Output classification: unassembled, bubble, contig
  uint32              _suggestRepeat   : 1;    //  Bogart made this from detected repeat.
  uint32              _suggestCircular : 1;    //  Bogart found overlaps making a circle.

  uint32              _spare           : 32 - 2 - 2;

  uint32              _layoutLen;         //  The max coord in the layout.  Same as gappedLen if it exists.

  char               *_gappedBases;       //  Gapped consensus - used by the multialignment.  NUL terminated.
  char               *_gappedQuals;
  uint32              _gappedLen;         //  Doesn't include the NUL.
  uint32              _gappedMax;         //  Does    include the NUL.

  char               *_ungappedBases;     //  Ungapped consensus - not used by the assember, only output.
  char               *_ungappedQuals;
  uint32              _ungappedLen;
  uint32              _ungappedMax;

  uint32             *_gappedToUngapped;  //  Map a gapped position to an ungapped posision, only output.

  tgPosition         *_children;          //  positions of objects that make up this tig
  uint32              _childrenLen;
  uint32              _childrenMax;

  int32              *_childDeltas;       //  deltas for all objects in the _children list
  uint32              _childDeltasLen;
  uint32              _childDeltasMax;

  //  Flags for computing consensus/multialignments

  uint32              _utgcns_verboseLevel;
};

#endif
