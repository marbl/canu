

#ifndef TG_TIG_H
#define TG_TIG_H

static const char *rcsid_TG_TIG_H = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"


//  Former IntMultiPos
//
//  Info about the placement of an object in a tig.  For unitigs, this
//  will be just reads.  For contigs, this could be unitigs and reads.
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

  //int32               &bgnset(void)             { return(_bgn); };
  //int32               &endset(void)             { return(_end); };

  int32          const bgn(void) const       { return((_isReverse == false) ? _min : _max); };
  int32          const end(void) const       { return((_isReverse == false) ? _max : _min); };

  int32                min(void) const       { return(_min); };
  int32                max(void) const       { return(_max); };

  //uint32               min(void) const       { return((_bgn < _end) ? _bgn : _end); };
  //uint32               max(void) const       { return((_bgn < _end) ? _end : _bgn); };

  uint32               deltaOffset(void)     { return(_deltaOffset); };
  uint32               deltaLength(void)     { return(_deltaLen);    };

  //  Set just the anchor and hangs, leave positions alone.
  void                 set(uint32 anchor, int32 ahang, int32 bhang) {
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
  void                 set(int32 bgn, int32 end) {
    if (bgn < end) {
      //assert(_isReverse == false);
      _min = bgn;
      _max = end;
    } else {
      //assert(_isReverse == true);
      _min = bgn;
      _max = end;
    }
  };

  //  Operators

  bool                 operator<(tgPosition const &that) const {
    int32 al = min();
    int32 bl = that.min();

    if (al != bl)
      return(al < bl);

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
  double              _microhetProb;

  uint32              _suggestRepeat    : 1;
  uint32              _suggestUnique    : 1;
  uint32              _suggestCircular  : 1;
  uint32              _suggestHaploid   : 1;
  uint32              _spare            : 28;

  uint32              _layoutLen;
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

  uint32               tigID(void)             { return(_tigID); };

  uint32               layoutLength(void)      { return(_layoutLen); };

  uint32               gappedLength(void)      { return(_gappedLen);   };
  char                *gappedBases(void)       { return(_gappedBases); }; 
  char                *gappedQuals(void)       { return(_gappedBases); }; 

private:
  void                 buildUngapped(void);
public:
  uint32               ungappedLength(void)    { buildUngapped();  return(_ungappedLen);   };
  char                *ungappedBases(void)     { buildUngapped();  return(_ungappedBases); }; 
  char                *ungappedQuals(void)     { buildUngapped();  return(_ungappedBases); }; 

  uint32               numberOfChildren(void)  { return(_childrenLen); };
  tgPosition          *getChild(uint32 c)      {
    assert(c < _childrenLen);
    return(_children + c);
  };
  tgPosition          *addChild(void)          {
    return(_children + _childrenLen++);
  };

  //  Operators

  void                 clear(void);    //  Clears data but doesn't release memory.

  tgTig               &operator=(tgTigRecord & tr);
  tgTig               &operator=(tgTig &tg);

  void                 saveToStream(FILE *F);
  void                 loadFromStream(FILE *F);

  void                 dumpLayout(FILE *F);
  bool                 loadLayout(FILE *F);

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
  double              microhetProb(void)  { return(_microhetProb); };

  double              _coverageStat;
  double              _microhetProb;

  //  Flags for computing consensus/multialignments

  uint32              _utgcns_verboseLevel;
  uint32              _utgcns_smoothWindow;
  bool                _utgcns_splitAlleles;
  bool                _utgcns_doPhasing;

  //  A variety of flags to suggest what type of unitig this is

  uint32              _suggestRepeat    : 1;
  uint32              _suggestUnique    : 1;
  uint32              _suggestCircular  : 1;
  uint32              _suggestHaploid   : 1;
  uint32              _spare:28;

  uint32              _layoutLen;       //  The max coord in the layout.  Same as gappedLen if it exists.

  char               *_gappedBases;     //  Gapped consensus - used by the multialignment
  char               *_gappedQuals;
  uint32              _gappedLen;
  uint32              _gappedMax;

  char               *_ungappedBases;   //  Ungapped consensus - not used by the assember, only output.
  char               *_ungappedQuals;
  uint32              _ungappedLen;
  uint32              _ungappedMax;

  tgPosition         *_children;        //  positions of objects that make up this tig
  uint32              _childrenLen;
  uint32              _childrenMax;

  int32              *_childDeltas;     //  deltas for all objects in the _children list
  uint32              _childDeltasLen;
  uint32              _childDeltasMax;
};


bool
generateMultiAlignment(tgTig     *tig,
                       gkStore   *gkpStore,
                       uint32    *failed,
                       double     errorRate,     //  AS_CNS_ERROR_RATE
                       double     errorRateMax,  //  AS_MAX_ERROR_RATE
                       uint32     minOverlap);   //  AS_OVERLAP_MIN_LEN


#endif
