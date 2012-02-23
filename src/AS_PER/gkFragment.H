/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

#ifndef GKFRAGMENT_H
#define GKFRAGMENT_H

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

static const char *AS_READ_ORIENT_NAMES[5] = {
  "U", "I", "O", "N", "A"
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
#define GKFRAGMENT_PACKED  0x01
#define GKFRAGMENT_NORMAL  0x02
#define GKFRAGMENT_STROBE  0x03

#define AS_READ_TYPE_NUM   4

static const char *AS_READ_TYPE_NAMES[4] = { "ERROR", "PACKED", "NORMAL", "STROBE" };



////////////////////////////////////////
//
//  Fragments
//
//  There are three types of fragments:
//
//  Packed -- less than 256 bp.
//  Normal -- up to 1048575 bases, but set at compile time.
//  Strobe -- for future use.
//
//  The three types are encapsulated in the gkFragment, which is the only supported access
//  mechanism.
//

class gkPackedFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad         : 4;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;
  uint32           seqLen      : AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearBeg    : AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearEnd    : AS_READ_MAX_PACKED_LEN_BITS;

#if 4 + 1 + 1 + 2 + 3 * AS_READ_MAX_PACKED_LEN_BITS != 32
#error gkPackedFragment size wrong
#endif

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


class  gkNormalFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad2        : 32 - 1 - 1 - 2;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;

  uint64           pad1        : 64 - 3 * AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           seqLen      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           clearBeg    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           clearEnd    : AS_READ_MAX_NORMAL_LEN_BITS;

  uint64           seqOffset;
  uint64           qltOffset;

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


class gkStrobeFragment {
private:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad2        : 32 - 1 - 1 - 2;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;

  uint64           pad1        : 64 - 3 * AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           seqLen      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           clearBeg    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64           clearEnd    : AS_READ_MAX_NORMAL_LEN_BITS;

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

#define gkFragment_get(ITEM)                \
  switch (type) {                           \
    case GKFRAGMENT_PACKED:                 \
      r = fr.packed.ITEM;                   \
      break;                                \
    case GKFRAGMENT_NORMAL:                 \
      r = fr.normal.ITEM;                   \
      break;                                \
    case GKFRAGMENT_STROBE:                 \
      r = fr.strobe.ITEM;                   \
      break;                                \
  }


#define gkFragment_set(ITEM, VAL)           \
  switch (type) {                           \
    case GKFRAGMENT_PACKED:                 \
      fr.packed.ITEM = (VAL);               \
      break;                                \
    case GKFRAGMENT_NORMAL:                 \
      fr.normal.ITEM = (VAL);               \
      break;                                \
    case GKFRAGMENT_STROBE:                 \
      fr.strobe.ITEM = (VAL);               \
      break;                                \
  }

class gkFragment {
public:
  gkFragment() {
   type = 0;
   tiid = 0;

   hasSEQ = 0;
   hasQLT = 0;

   isGKP = 0;

   clrBgn = clrEnd = 0;
   vecBgn = vecEnd = 0;
   maxBgn = maxEnd = 0;
   tntBgn = tntEnd = 0;

   gkp = NULL;

   memset(&fr, 0, sizeof(gkFragmentData));

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
    uint32 r=0, e=0;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearBeg);
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, r, e);
    }
    return(r);
  };
  uint32 gkFragment_getClearRegionEnd  (uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r=0, b=0;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearEnd);
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, b, r);
    }
    return(r);
  };

  uint32 gkFragment_getClearRegionLength(uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r=0, b=0, e=0;
    if (which == AS_READ_CLEAR_LATEST) {
      gkFragment_get(clearBeg);  b = r;
      gkFragment_get(clearEnd);  e = r;
    } else {
      gkp->clearRange[which]->gkClearRange_getClearRegion(this, b, e);
    }
    return(e - b);
  };

  void   gkFragment_getClearRegion(uint32 &begin, uint32 &end, uint32 which=AS_READ_CLEAR_LATEST) {
    uint32 r=0;
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
    AS_IID r = 0;
    gkFragment_get(readIID);
    return(r);
  };

  AS_IID      gkFragment_getMateIID(void) {
    AS_IID r = 0;
    gkFragment_get(mateIID);
    return(r);
  };

  AS_IID      gkFragment_getLibraryIID(void) {
    AS_IID r = 0;
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
      case GKFRAGMENT_PACKED:
        r = 0;
        break;
      case GKFRAGMENT_NORMAL:
        r = fr.normal.seqOffset;
        //assert(fr.normal.qltOffset == 4 * tiid + 4 * fr.normal.seqOffset);
        break;
      case GKFRAGMENT_STROBE:
        r = fr.strobe.seqOffset;
        break;
    }
    return(r);
  }

  uint64      gkFragment_getQualityOffset(void) {
    uint64 r = 0;
    switch (type) {
      case GKFRAGMENT_PACKED:
        r = 0;
        break;
      case GKFRAGMENT_NORMAL:
        r = fr.normal.qltOffset;
        //assert(fr.normal.qltOffset == 4 * tiid + 4 * fr.normal.seqOffset);
        break;
      case GKFRAGMENT_STROBE:
        r = fr.strobe.qltOffset;
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
      enc = (char *)safe_malloc(sizeof(char) * (AS_READ_MAX_NORMAL_LEN + 1));
      seq = (char *)safe_malloc(sizeof(char) * (AS_READ_MAX_NORMAL_LEN + 1));
      qlt = (char *)safe_malloc(sizeof(char) * (AS_READ_MAX_NORMAL_LEN + 1));
    }
  };

  void        gkFragment_setType(uint32 t)        { assert(isGKP);  type = t; };
  void        gkFragment_setReadUID(AS_UID u)     {                 gkFragment_set(readUID, u); };
  void        gkFragment_setLength(uint32 l)      { assert(isGKP);  gkFragment_set(seqLen, l); };
  void        gkFragment_setLibraryIID(AS_IID l)  {                 gkFragment_set(libraryIID, l); };
  void        gkFragment_setMateIID(AS_IID i)     {                 gkFragment_set(mateIID, i); };

  void        gkFragment_setOrientation(uint32 i) {                 gkFragment_set(orientation, i); };
  void        gkFragment_setIsDeleted(uint32 i)   {                 gkFragment_set(deleted, i); };
  void        gkFragment_setIsNonRandom(uint32 i) {                 gkFragment_set(nonrandom, i); };

  void        gkFragment_clear(void)              { memset(&fr, 0, sizeof(gkFragmentData)); };

private:
public:
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

  union gkFragmentData {
    gkPackedFragment   packed;  //  was sm
    gkNormalFragment   normal;  //  was md
    gkStrobeFragment   strobe;  //  was ;g
  } fr;

  char   *enc;
  char   *seq;
  char   *qlt;

  friend class gkStore;
  friend class gkStream;
  friend class gkClearRange;
};

#endif
