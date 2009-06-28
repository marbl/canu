


////////////////////////////////////////
//
//  Mate orientations.  UNFORTUNATELY, these are used in CGW too.  The
//  ANTI orientation is not a valid mate orientation.
//
#define AS_READ_ORIENT_UNKNOWN    0x00
#define AS_READ_ORIENT_INNIE      0x01
#define AS_READ_ORIENT_OUTTIE     0x02
#define AS_READ_ORIENT_NORMAL     0x03
#define AS_READ_ORIENT_ANTINORMAL 0x04

static const char *AS_READ_ORIENT_NAMES[4] = {
  "U", "I", "O", "N",
};


////////////////////////////////////////
//
//  Clear ranges.  The symbols are used internally to the code, while
//  the names are used for files and interacting with users.
//
#define AS_READ_CLEAR_LATEST         0
#define AS_READ_CLEAR_CLR            1  //  Original FRG message clear range
#define AS_READ_CLEAR_VEC            2  //  Original FRG message vector clear
#define AS_READ_CLEAR_MAX            3  //  Original FRG message maximum clear
#define AS_READ_CLEAR_TNT            4  //  Original FRG message contaminated region (tainted)
#define AS_READ_CLEAR_OBTINITIAL     5
#define AS_READ_CLEAR_OBTMERGE       6
#define AS_READ_CLEAR_OBTCHIMERA     7
#define AS_READ_CLEAR_ECR_0          8  //  ECR starting clear range; either OBTCHIMERA or CLR
#define AS_READ_CLEAR_ECR_1          9
#define AS_READ_CLEAR_ECR_2         10
#define AS_READ_CLEAR_ECR_3         11
#define AS_READ_CLEAR_ECR_4         12
#define AS_READ_CLEAR_ECR_5         13
#define AS_READ_CLEAR_ECR_6         14
#define AS_READ_CLEAR_ECR_7         15
#define AS_READ_CLEAR_ECR_8         16

#define AS_READ_CLEAR_ERROR        255

#define AS_READ_CLEAR_NUM           17

static const char *AS_READ_CLEAR_NAMES[AS_READ_CLEAR_NUM] = {
  "LATEST",
  "CLR", "VEC", "MAX", "TAINT",
  "OBTINITIAL", "OBTMERGE", "OBTCHIMERA",
  "ECR_0", "ECR_1", "ECR_2", "ECR_3", "ECR_4", "ECR_5", "ECR_6", "ECR_7", "ECR_8"
};


////////////////////////////////////////
//
//  gkFragment types -- the first group tells which data to load, the
//  second is an internal flag of what type of record this is.
//
#define GKFRAGMENT_INF     0x00
#define GKFRAGMENT_SEQ     0x01
#define GKFRAGMENT_QLT     0x02

#define GKFRAGMENT_ERROR   0x00
#define GKFRAGMENT_SHORT   0x01
#define GKFRAGMENT_MEDIUM  0x02
#define GKFRAGMENT_LONG    0x03

#define AS_READ_TYPE_NUM   4

static const char *AS_READ_TYPE_NAMES[4] = { "ERROR", "SHORT", "MEDIUM", "LONG" };



////////////////////////////////////////
//
//  Fragments
//
//  There are three types of fragments:
//
//  Short  -- less than 256 bp.  These fragments are special in that they
//            do not have a clear range.
//
//  Medium -- up to 2047 bases.  These are allowed to be trimmed and extended.
//
//  Long   -- larger than 2047 bases.  These are allowed to be split.
//
//  The three types are encapsulated in the gkFragment, which is the
//  only supported access mechanism.
//
class gkShortFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad         : 4;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;
  uint32           seqLen      : AS_READ_MAX_SHORT_LEN_BITS;
  uint32           clearBeg    : AS_READ_MAX_SHORT_LEN_BITS;
  uint32           clearEnd    : AS_READ_MAX_SHORT_LEN_BITS;

#if 4 + 1 + 1 + 2 + 3 * AS_READ_MAX_SHORT_LEN_BITS != 32
#error gkShortFragment size wrong
#endif

  char             enc[AS_READ_MAX_SHORT_LEN];

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


class  gkMediumFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint64           pad2        : 27;
  uint64           deleted     : 1;
  uint64           nonrandom   : 1;
  uint64           orientation : 2;
  uint64           seqLen      : AS_READ_MAX_MEDIUM_LEN_BITS;
  uint64           clearBeg    : AS_READ_MAX_MEDIUM_LEN_BITS;
  uint64           clearEnd    : AS_READ_MAX_MEDIUM_LEN_BITS;

#if 27 + 1 + 1 + 2 + 3 * AS_READ_MAX_MEDIUM_LEN_BITS != 64
#error gkMediumFragment size wrong
#endif

  uint64           seqOffset;
  uint64           qltOffset;

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


class gkLongFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad2        : 28;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;

#if 28 + 1 + 1 + 2 != 32
#error gkLongFragment size wrong
#endif

  uint64           pad1        : 4;
  uint64           seqLen      : AS_READ_MAX_LONG_LEN_BITS;
  uint64           clearBeg    : AS_READ_MAX_LONG_LEN_BITS;
  uint64           clearEnd    : AS_READ_MAX_LONG_LEN_BITS;

#if 4 + 3 * AS_READ_MAX_LONG_LEN_BITS != 64
#error gkLongFragment size wrong
#endif

  uint64           seqOffset;
  uint64           qltOffset;

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};






////////////////////////////////////////
//
//  gkFragment -- the accessor to the store.
//

#define gkFragment_get(ITEM)                    \
  switch (type) {                           \
    case GKFRAGMENT_SHORT:                      \
      r = fr.sm.ITEM;                       \
      break;                                    \
    case GKFRAGMENT_MEDIUM:                     \
      r = fr.md.ITEM;                       \
      break;                                    \
    case GKFRAGMENT_LONG:                       \
      r = fr.lg.ITEM;                       \
      break;                                    \
  }


#define gkFragment_set(ITEM, VAL)          \
  switch (type) {                           \
    case GKFRAGMENT_SHORT:                      \
      fr.sm.ITEM = (VAL);                 \
      break;                                    \
    case GKFRAGMENT_MEDIUM:                     \
      fr.md.ITEM = (VAL);                 \
      break;                                    \
    case GKFRAGMENT_LONG:                       \
      fr.lg.ITEM = (VAL);                 \
      break;                                    \
  }

class gkFragment {
public:
  gkFragment() {
   type = 0;
   tiid = 0;

   hasSEQ = 0;
   hasQLT = 0;

   isGKP = 0;

   gkp = NULL;

   enc = NULL;
   seq = NULL;
   qlt = NULL;
  };
  ~gkFragment() {
    safe_free(enc);
    safe_free(seq);
    safe_free(qlt);
  };


  uint32 gkFragment_getClearRegionBegin(uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r, e;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearBeg);
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, r, e);
    }
    return(r);
  };
  uint32 gkFragment_getClearRegionEnd  (uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r, b;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearEnd);
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, b, r);
    }
    return(r);
  };

  uint32 gkFragment_getClearRegionLength(uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r, b, e;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearBeg);  b = r;
      gkFragment_get(clearEnd);  e = r;
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, b, e);
    }
    return(e - b);
  };

  void   gkFragment_getClearRegion(uint32 &begin, uint32 &end, uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearBeg);  begin = r;
      gkFragment_get(clearEnd);  end   = r;
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, begin, end);
    }
  };
  void   gkFragment_setClearRegion(uint32  begin, uint32  end, uint32 which) {
    assert(gkp->isReadOnly == 0);
    gkFragment_set(clearBeg, begin);
    gkFragment_set(clearEnd, end);
    if (which != AS_READ_CLEAR_LATEST)
      gkp->clearRange[which]->gkClearRange_setClearRegion(this, begin, end);
  }


  AS_UID      gkFragment_getReadUID(void) {
    AS_UID r = AS_UID_undefined();
    gkFragment_get(readUID);
    return(r);
  };

  AS_IID      gkFragment_getReadIID(void) {
    AS_IID r = -1;
    gkFragment_get(readIID);
    return(r);
  };

  AS_IID      gkFragment_getMateIID(void) {
    AS_IID r = -1;
    gkFragment_get(mateIID);
    return(r);
  };

  AS_IID      gkFragment_getLibraryIID(void) {
    AS_IID r = -1;
    gkFragment_get(libraryIID);
    return(r);
  };

  uint32      gkFragment_getIsDeleted(void) {
    uint32 r = 0;
    gkFragment_get(deleted);
    return(r);
  };

  uint32      gkFragment_getIsNonRandom(void) {
    uint32 r = 0;
    gkFragment_get(nonrandom);
    return(r);
  };

  uint32      gkFragment_getOrientation(void) {
    uint32 r = 0;
    gkFragment_get(orientation);
    return(r);
  };

  uint32      gkFragment_getSequenceLength(void) {
    uint32 r = 0;
    gkFragment_get(seqLen);
    return(r);
  }

  uint32      gkFragment_getQualityLength(void) {
    uint32 r = 0;
    gkFragment_get(seqLen);
    return(r);
  }

  uint64      gkFragment_getSequenceOffset(void) {
    uint64 r = 0;
    switch (type) {
      case GKFRAGMENT_SHORT:
        r = 0;
        break;
      case GKFRAGMENT_MEDIUM:
        r = fr.md.seqOffset;
        //assert(fr.md.qltOffset == 4 * tiid + 4 * fr.md.seqOffset);
        break;
      case GKFRAGMENT_LONG:
        r = fr.lg.seqOffset;
        break;
    }
    return(r);
  }

  uint64      gkFragment_getQualityOffset(void) {
    uint64 r = 0;
    switch (type) {
      case GKFRAGMENT_SHORT:
        r = 0;
        break;
      case GKFRAGMENT_MEDIUM:
        r = fr.md.qltOffset;
        //assert(fr.md.qltOffset == 4 * tiid + 4 * fr.md.seqOffset);
        break;
      case GKFRAGMENT_LONG:
        r = fr.lg.qltOffset;
        break;
    }
    return(r);
  }

  //  Both of these; isGKP bypasses the usual check of if we have
  //  called getFragmentData with the correct flags since all it wants
  //  is the memory.
  char       *gkFragment_getSequence(void) {
    assert(isGKP || hasSEQ);
    return(seq);
  }
  char       *gkFragment_getQuality(void) {
    assert(isGKP || hasQLT);
    return(qlt);
  }


  //  For use ONLY by AS_GKP

  void        gkFragment_enableGatekeeperMode(gkStore *g) {
    if (isGKP == 0) {
      isGKP = 1;

      hasSEQ = 1;
      hasQLT = 1;

      gkp = g;
      enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
    }
  };

  void        gkFragment_setType(uint32 t)        { assert(isGKP);  type = t; };
  void        gkFragment_setReadUID(AS_UID u)     { assert(isGKP);  gkFragment_set(readUID, u); };
  void        gkFragment_setLength(uint32 l)      { assert(isGKP);  gkFragment_set(seqLen, l); };
  void        gkFragment_setLibraryIID(AS_IID l)  { assert(isGKP);  gkFragment_set(libraryIID, l); };
  void        gkFragment_setMateIID(AS_IID i)     { assert(isGKP);  gkFragment_set(mateIID, i); };

  void        gkFragment_setOrientation(uint32 i) { assert(isGKP);  gkFragment_set(orientation, i); };
  void        gkFragment_setIsDeleted(uint32 i)   { assert(isGKP);  gkFragment_set(deleted, i); };
  void        gkFragment_setIsNonRandom(uint32 i) { assert(isGKP);  gkFragment_set(nonrandom, i); };

private:
  uint32   type;
  uint32   tiid;

  uint32   hasSEQ;
  uint32   hasQLT;

  uint32   isGKP;

public:
  uint32   clrBgn, clrEnd;  //  For use by gatekeeper and sffToCA ONLY.
  uint32   vecBgn, vecEnd;  //  These are the initial clear ranges to be
  uint32   maxBgn, maxEnd;  //  loaded into the store.
  uint32   tntBgn, tntEnd;  //  DO NOT USE!

private:
  gkStore *gkp;

  union {
    gkShortFragment    sm;
    gkMediumFragment   md;
    gkLongFragment     lg;
  } fr;

  char   *enc;
  char   *seq;
  char   *qlt;

  friend class gkStore;
  friend class gkStream;
  friend class gkClearRange;
};

