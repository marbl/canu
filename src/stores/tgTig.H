

#ifndef TG_TIG_H
#define TG_TIG_H

static const char *rcsid_TG_TIG_H = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"

//  Old, for consensus
//  Former IntVarAllele
class tgAllelePosition {
public:
  tgAllelePosition() {};
  ~tgAllelePosition() {};

private:
public:
  uint32           num_reads;
  uint32           weight;
  uint32           var_seq_offset;  //  offset into var_seq_memory for this allele
  uint32           read_id_offset;  //  offset into read_id_memory for this allele
};


//  Former IntMultiVar
class tgVariantPosition {
public:
  tgVariantPosition() {};
  ~tgVariantPosition() {};

  void                 saveToStream(FILE *F);    //  These do NOT save the object, just the arrays
  void                 loadFromStream(FILE *F);


private:
public:
  int32              var_id;                 // id of current VAR record
  int32              phased_id;              // id of the VAR record phased with the current one

  uint32             bgn;                    // position of the var region
  uint32             end;

  int32              num_reads;              // num reads total
  int32              num_alleles;            // num alleles total
  int32              num_alleles_confirmed;  // num confirmed alleles
  int32              min_anchor_size;        //
  int32              var_length;             // bases in the var region

  tgAllelePosition  *alleles;                // list of num_alleles alleles

  char              *var_seq_memory;         // single allocation for all memory
  int32             *read_id_memory;         // single allocation for all memory

  char              *enc_num_reads;          //  the nra field
  char              *enc_weights;            //  the wgt field
  char              *enc_var_seq;            //  the seq field
  char              *enc_read_ids;           //  the rid field
};


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

  uint32              &anchor(void)          { return(_anchor); };
  int32               &aHang(void)           { return(_ahang);  };
  int32               &bHang(void)           { return(_bhang);  };

  bool                 isForward(void) const { return(_bgn < _end); };
  bool                 isReverse(void) const { return(_bgn > _end); };

  uint32              &bgn(void)             { return(_bgn); };
  uint32              &end(void)             { return(_end); };

  uint32               bgn(void) const       { return(_bgn); };
  uint32               end(void) const       { return(_end); };

  uint32               min(void) const       { return((_bgn < _end) ? _bgn : _end); };
  uint32               max(void) const       { return((_bgn < _end) ? _end : _bgn); };

  uint32               deltaOffset(void)     { return(_deltaOffset); };
  uint32               deltaLength(void)     { return(_deltaLen);    };

  //  Operators

  void                 saveToStream(FILE *F);    //  These do NOT save the object, just the deltas
  void                 loadFromStream(FILE *F);


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

  uint32              _isRead      : 1;  //  Type of object this position represents, former FragType
  uint32              _isUnitig    : 1;  //  
  uint32              _isContig    : 1;  //

  uint32              _spare       : 29;

  uint32              _anchor;           //  ID of the like-type object we align to
  int32               _ahang;            //  Placement relative to anchor object
  int32               _bhang;            //

  uint32              _bgn;              //  Coords in the parent object we are part of
  uint32              _end;              //  (for a read, this will be the position in the unitig)

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
  uint32              _spare:28;

  uint32              _layoutLen;
  uint32              _gappedLen;       //  Gapped consensus
  uint32              _ungappedLen;     //  Ungapped consensus
  uint32              _childrenLen;
  uint32              _childDeltasLen;
  uint32              _variantsLen;
  uint32              _variantDataLen;
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

  uint32               ungappedLength(void)    { return(_ungappedLen);   };
  char                *ungappedBases(void)     { return(_ungappedBases); }; 
  char                *ungappedQuals(void)     { return(_ungappedBases); }; 

  uint32               numberOfChildren(void)  { return(_childrenLen); };
  tgPosition          *getChild(uint32 c)      {
    assert(c < _childrenLen);
    return(_children + c);
  };
  tgPosition          *addChild(uint32 c)      {
    assert(c == _childrenLen);
    _childrenLen++;
    return(_children + c);
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

  tgVariantPosition  *_variants;        //  variations
  uint32              _variantsLen;
  uint32              _variantsMax;

  uint8              *_variantData;
  uint64              _variantDataLen;
  uint64              _variantDataMax;
};


bool
generateMultiAlignment(tgTig                     *tig,
                       gkStore                   *gkpStore,
                       uint32                    *failed);


#endif
