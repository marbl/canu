
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

#ifndef GKSTORE_H
#define GKSTORE_H

#include "AS_UTL_Hash.H"

////////////////////////////////////////
//
//  The gkpStore
//
//  gkStoreInfo is saved on disk.
//  gkStore is the in memory structure used to access the data.
//
class gkStoreInfo {
private:
public:
  uint64    gkMagic;
  uint64    gkVersion;

  uint32    gkLibrarySize;
  uint32    gkPackedSequenceSize;
  uint32    gkPackedFragmentSize;
  uint32    gkNormalFragmentSize;
  uint32    gkStrobeFragmentSize;
  uint32    gkPlacementSize;

  uint32    readMaxLenBits;

  //  Statistics on our load

  uint32    libInput;
  uint32    libLoaded;
  uint32    libErrors;
  uint32    libWarnings;

  uint32    frgInput;
  uint32    frgLoaded;
  uint32    frgErrors;
  uint32    frgWarnings;

  uint32    lkgInput;
  uint32    lkgLoaded;
  uint32    lkgErrors;
  uint32    lkgWarnings;
   
  uint32    sffInput;
  uint32    sffLoaded;
  uint32    sffErrors;
  uint32    sffWarnings;

  uint32    sffLibCreated;

  uint32    plcInput;
  uint32    plcLoaded;
  uint32    plcErrors;
  uint32    plcWarnings;

  //  Counts of types of things we have loaded

  uint32    numRandom;

  //  The IID space is broken into three classes.  See comments in AS_PER_gkStore_IID.C

  uint32    numPacked;
  uint32    numNormal;
  uint32    numStrobe;

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


//gkStore *
//gkStoreConstruct(const char *path, uint32 packedLength);

class gkStore {
private:
  void         gkStore_open(int writable);
  //void         gkStore_create(const char *path, uint32 packedLength);

  void         gkStore_clear(void);
  
public:
  gkStore();
  gkStore(const char *path, uint32 packedLength);
  gkStore(const char *path, int creatable, int writable);
  ~gkStore();

  //  The 'name' of this gkpStore.
  const char  *gkStore_path(void) { return(storePath); };

  ////////////////////////////////////////
  //
  //  Operatons on the whole store.  Loading into memory, removing from disk.
  //
  ////////////////////////////////////////

  void         gkStore_load(AS_IID firstElem, AS_IID lastElem, int flags);

  void         gkStore_loadPartition(uint32 partition);
  void         gkStore_buildPartitions(short *partitionMap, uint32 maxPart);

  void         gkStore_delete(void);

  uint64       gkStore_metadataSize(void);
  void         gkStore_metadataCaching(bool enable=true);

  ////////////////////////////////////////
  //
  //  AS_PER_gkStore.c
  //
  ////////////////////////////////////////

  uint32      gkStore_getNumLibraries(void)       { return(inf.libLoaded); }
  uint32      gkStore_getNumFragments(void)       { return(inf.frgLoaded); }
  uint32      gkStore_getNumRandomFragments(void) { return(inf.numRandom); }

  gkLibrary *gkStore_getLibrary(int id)                { return((id == 0) ? libNull : (gkLibrary *)getIndexStorePtr(lib, id)); }
  void       gkStore_getLibrary(int id, gkLibrary *lb) { getIndexStore(lib, id, lb); };
  void       gkStore_setLibrary(int id, gkLibrary *lb) { setIndexStore(lib, id, lb); }
  void       gkStore_addLibrary(gkLibrary *lb) {
    appendIndexStore(lib, lb);
    inf.libLoaded++;
  };

  void      gkStore_decodeTypeFromIID(AS_IID iid, uint32& type, uint32& tiid);
  void      gkStore_addIIDtoTypeMap(AS_IID iid, uint32 type, uint32 tiid);
  void      gkStore_computeRanges(AS_IID  bgnIID, AS_IID  endIID,
                                  int64 &bgnPK,  int64 &endPK,  int64 &valPK,
                                  int64 &bgnNM,  int64 &endNM,  int64 &valNM,
                                  int64 &bgnSB,  int64 &endSB,  int64 &valSB);

  void      gkStore_getFragmentData(gkStream *gst, gkFragment *fr, uint32 flags);

  void      gkStore_getFragment(AS_IID iid, gkFragment *fr, int32 flags);
  void      gkStore_setFragment(gkFragment *fr);
  void      gkStore_delFragment(AS_IID iid, bool deleteMateFrag=false);
  void      gkStore_addFragment(gkFragment *fr);

  ////////////////////////////////////////
  //
  //  Clear range support.
  //
  ////////////////////////////////////////

  //  Add a 'new' clear range.  Initial values are set to the version
  //  in the fragment.  Without enabling a range, only
  //  getClearRegion() works, and returns [1,0].
  //
  void      gkStore_enableClearRange(uint32 which);
  void      gkStore_purgeClearRange(uint32 which);

  ////////////////////////////////////////
  //
  //
  //
  ////////////////////////////////////////

  void   gkStore_flagWarning(char t);
  void   gkStore_flagError(char t);

private:
  char                     storePath[FILENAME_MAX];

  uint32                   isReadOnly;
  uint32                   isCreating;

public:    //  Sigh, public needed for AS_GKP.
  gkStoreInfo              inf;

private:
  StoreStruct             *fpk;  //  Packed fragments
  StoreStruct             *qpk;

  StoreStruct             *fnm;  //  Normal fragments
  StoreStruct             *snm;
  StoreStruct             *qnm;

  StoreStruct             *fsb;  //  Strobe fragments
  StoreStruct             *ssb;
  StoreStruct             *qsb;

  StoreStruct             *lib;
  gkLibrary               *libNull;  //  Pointer, because we don't know gkLibrary here
 
  //  Clear ranges
  //
  gkClearRange           **clearRange;

  //  If this is defined (not NULL) then the store contains reads of mixed types, AND, it was
  //  created 'out of order' (in order: packed first, then normal, then strobe).  See comments in
  //  AS_PER_gkStore_IID.C
  //
  uint64                   IIDmax;
  uint8                   *IIDtoTYPE;
  uint32                  *IIDtoTIID;

  //  The rest are for a partitioned fragment store.
  //
  //  We load all frg and qlt in this partition into memory.  The map converts an iid (global
  //  fragment iid) into a pointer to the correct frg record, which we can then use to grab the
  //  encoded seq/qlt.
  //
  int32                    partnum;

  StoreStruct             *partfpk, *partqpk;
  StoreStruct             *partfnm, *partqnm;
  StoreStruct             *partfsb, *partqsb;

  HashTable_AS            *partmap;

  friend class gkStream;
  friend class gkFragment;  //  for clearRange
  friend class gkClearRange;
};


class gkStoreStats {
public:
  gkStoreStats(char const *gkStoreName);
  gkStoreStats(gkStore *gkp);
  ~gkStoreStats();

  void    init(gkStore *gkp);

  uint32  numActiveFrag;
  uint32  numDeletedFrag;
  uint32  numMatedFrag;
  uint64  readLength;
  uint64  clearLength;

  uint32  *lowestIID;
  uint32  *highestIID;

  uint32  *numActivePerLib;
  uint32  *numDeletedPerLib;
  uint32  *numMatedPerLib;
  uint64  *readLengthPerLib;
  uint64  *clearLengthPerLib;
};


#endif
