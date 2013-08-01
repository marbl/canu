
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

#ifndef MULTIALIGNSTORE_H
#define MULTIALIGNSTORE_H

static const char *rcsid_MULTIALIGNSTORE_H = "$Id$";

#include "AS_global.H"
#include "MultiAlign.H"

//
//  The MultiAlignStore is a disk-resident (with memory cache) database of MultiAlign (MA)
//  structures.
//
//  The directory structure looks like:
//    x.maStore/
//    x.maStore/v001.dat       x.maStore/v001.utg       x.maStore/v001.ctg
//    x.maStore/v002.p001.dat  x.maStore/v002.p001.dat  x.maStore/v002.p001.dat
//    x.maStore/v002.p002.dat  x.maStore/v002.p001.utg  x.maStore/v002.p002.ctg
//    x.maStore/v002.p003.dat  x.maStore/v002.p001.utg  x.maStore/v002.p003.ctg
//
//  Showing two "versions" of data (v001 and v002), with the second version being "partitioned" into
//  three sets (p001, p002, p003).
//
//  The MA structures are stored in the 'dat' files, in the order they are written.  Multiple copies
//  of the same MA can be present in each file, for example, if the same MA is changed twice.
//
//  The 'utg' and 'ctg' files store an array of metadata (the MAR struct below) for each MA.  The
//  primary information in the metadata is where the latest version of a MA structure is stored --
//  the version, partition and position in the file.
//
//  For partitioned data, each 'utg' and 'ctg' file contains metadata for ALL MAs, even those not in
//  the partition.  The metadata is only valid for the current partition.  The store explicitly
//  disallows access to an MA not in the current partition.  For example, v002.p003.utg contains
//  metadata for all unitigs, but only unitigs in partition 3 are guaranteed to be up-to-date.  When
//  the store is next opened 'unpartitioned' it will consolidate the metadata from all partitions.
//

class MultiAlignStore {
public:

  //  Create a MultiAlignStore (first constructor).
  //
  //    If the partitionMap is supplied, a partitioned store is created by default, placing MA i
  //    into partition partitionMap[i].  nextVersion() is NOT ALLOWED here.
  //
  //    If the partitionMap is not supplied, an unpartitioned store is created.  nextVersion() is
  //    allowed.
  //
  //  Open a MultiAlignStore (second constructor).
  //
  //    If 'partition' is non-zero, then only MAs in that partition are allowed to be accessed, and
  //    any writes will maintain the partitioning.  In particular, writes to partitions are
  //    independent.
  //
  //    If 'partition' is zero, any previous partitioning is merged to form a single partition.  If
  //    writable, the next version will be unpartitioned.  Note that data is still stored in
  //    partitioned files, it is not copied to an unpartitioned file.
  //
  MultiAlignStore(const char *path);
  MultiAlignStore(const char *path,
                  uint32      version,
                  uint32      unitigPartition,
                  uint32      contigPartition,
                  bool        writable=false,
                  bool        inplace=false,
                  bool        append=false);
  ~MultiAlignStore();

  //  Update to the next version.  Fails if the store is opened partitioned -- there is no decent
  //  way to ensure that all partitions will be at the same version.
  //
  void           nextVersion(void);

  //  Switch from writing non-partitioned data to writing partitioned data.  As usual, calling
  //  nextVersion() after this will fail.  Contigs that do not get placed into a partition will
  //  still exist in the (unpartitioned) store, but any clients opening a specific partition will
  //  not see them.
  //
  //  Suppose we have three contigs, A, B and C.  We place A and B in partition 1, but do not touch
  //  C.  Clients open partitions and process contigs.  Since C is not in a partition, it is never
  //  processed.  Later, the store is opened unpartitioned.  We now see all three contigs.
  //
  //
  void           writeToPartitioned(uint32 *unitigPartMap, uint32 unitigPartMapLen,
                                    uint32 *contigPartMap, uint32 contigPartMapLen);

  
  //  Add or update a MA in the store.  If keepInCache, we keep a pointer to the MultiAlignT.  THE
  //  STORE NOW OWNS THE OBJECT.
  //
  void           insertMultiAlign(MultiAlignT *ma, bool isUnitig, bool keepInCache);

  //  delete() removes the tig from the cache, and marks it as deleted in the store.
  //
  void           deleteMultiAlign(int32 maID, bool isUnitig);

  //  load() will load and cache the MA.  THE STORE OWNS THIS OBJECT.
  //  copy() will load and copy the MA.  It will not cache.  YOU OWN THIS OBJECT.
  //
  MultiAlignT   *loadMultiAlign(int32 maID, bool isUnitig);
  void           unloadMultiAlign(int32 maID, bool isUnitig);

  void           copyMultiAlign(int32 maID, bool isUnitig, MultiAlignT *ma);

  //  Flush to disk any cached MAs.  This is called by flushCache().
  //
  void           flushDisk(int32 maID, bool isUnitig);
  void           flushDisk(void);

  //  Flush the cache of loaded MAs.  Be aware that this is expensive in that the flushed things
  //  usually just get loaded back into core.
  //
  void           flushCache(void);

  uint32         numUnitigs(void) { return(utgLen); };
  uint32         numContigs(void) { return(ctgLen); };

  //  Accessors to MultiAlignD data; these do not load the multialign.

  bool           isDeleted(int32 maID, bool isUnitig);

  int32          getUnitigCoverageStat(int32 maID);
  double         getUnitigMicroHetProb(int32 maID);
  UnitigStatus   getUnitigStatus(int32 maID);
  UnitigFUR      getUnitigFUR(int32 maID);
  ContigStatus   getContigStatus(int32 maID);

  uint32         getNumFrags(int32 maID, bool isUnitig);
  uint32         getNumUnitigs(int32 maID, bool isUnitig);

  void           setUnitigCoverageStat(int32 maID, double cs);
  void           setUnitigMicroHetProb(int32 maID, double mp);
  void           setUnitigStatus(int32 maID, UnitigStatus status);
  void           setUnitigFUR(int32 maID, UnitigFUR fur);

  void           setContigStatus(int32 maID, ContigStatus status);

  uint32         getUnitigVersion(int32 maID);
  uint32         getContigVersion(int32 maID);

  void           dumpMultiAlignR(int32 maID, bool isUnitig);
  void           dumpMultiAlignRTable(bool isUnitig);

private:
  struct MultiAlignR {
    MultiAlignD  mad;
    uint64       unusedFlags : 1;   //  One whole bit for future use.
    uint64       flushNeeded : 1;   //  If true, this MAR and associated MultiAlign are NOT saved to disk.
    uint64       isPresent   : 1;   //  If true, this MAR is present in this partition.
    uint64       isDeleted   : 1;   //  If true, this MAR has been deleted from the assembly.
    uint64       ptID        : 10;  //  10 -> 1024 partitions
    uint64       svID        : 10;  //  10 -> 1024 versions
    uint64       fileOffset  : 40;  //  40 -> 1 TB file size; offset in file where MA is stored
  };

  void                    init(const char *path_, uint32 version_, bool writable_, bool inplace_, bool append_);

  void                    writeTigToDisk(MultiAlignT *ma, MultiAlignR *maRecord);

  void                    dumpMASRfile(char *name, MultiAlignR *R, uint32 L, uint32 M, uint32 part);
  bool                    loadMASRfile(char *name, MultiAlignR *R, uint32 L, uint32 M, uint32 part, bool onlyThisV);
  uint32                  numTigsInMASRfile(char *name);

  void                    dumpMASR(MultiAlignR* &R, uint32& L, uint32& M, uint32 V, bool isUnitig);
  void                    loadMASR(MultiAlignR* &R, uint32& L, uint32& M, uint32 V, bool isUnitig, bool onlyThisV);

  void                    purgeVersion(int version);
  void                    purgeCurrentVersion(void);

  friend void operationCompress(char *tigName, int tigVers);

  FILE                   *openDB(uint32 V, uint32 P);

  char                    path[FILENAME_MAX];
  char                    name[FILENAME_MAX];

  bool                    writable;               //  We are able to write
  bool                    inplace;                //  We read and write to the same version
  bool                    append;                 //  Do not nuke an existing partition

  bool                    newTigs;                //  internal flag, set if tigs were added

  uint32                  originalVersion;        //  Version we started from (see newTigs in code)
  uint32                  currentVersion;         //  Version we are writing to

  //  Creating or changing the partitioning.  These act independently, though it (currently) makes
  //  little sense to change the unitig partitioning when the contigs are partitioned.
  //
  uint32                 *unitigPartMap;          //  Unitig partitioning
  uint32                  unitigPartMapLen;       //  Unitig partitioning

  uint32                 *contigPartMap;          //  Contig partitioning
  uint32                  contigPartMapLen;       //  Contig partitioning

  //  Loading restrictions; if these are non-zero, we can only load tigs in this partition.
  //  Attempts to load other tigs result in NULL returns.  On write, these will set or override the
  //  setting of ptID in the MultiAlignR.
  //
  //  Only one may be set at a time.
  //
  //    If unitigPart is set, no contigs may be loaded.  Only unitigs from the same partition are
  //    loaded, others simply return NULL pointers.
  //
  //    If contigPart is set, all unitigs may be loaded, but any writes will repartition the unitig
  //    to this partition.  Like unitigs, loading a contig from a different partition results in
  //    a NULL pointer.
  //
  uint32                  unitigPart;             //  Partition we are restricted to read from
  uint32                  contigPart;             //  Partition we are restricted to read from


  uint32                  utgMax;
  uint32                  utgLen;
  MultiAlignR            *utgRecord;
  MultiAlignT           **utgCache;

  uint32                  ctgMax;
  uint32                  ctgLen;
  MultiAlignR            *ctgRecord;
  MultiAlignT           **ctgCache;

  struct dataFileT {
    FILE   *FP;
    bool    atEOF;
  };

  dataFileT             **dataFile;       //  dataFile[version][partition] = FP
};


inline
bool
MultiAlignStore::isDeleted(int32 maID, bool isUnitig) {
  assert(maID >= 0);
  return((isUnitig) ? utgRecord[maID].isDeleted : ctgRecord[maID].isDeleted);
}

inline
int
MultiAlignStore::getUnitigCoverageStat(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  return((int)utgRecord[maID].mad.unitig_coverage_stat);
}

inline
double
MultiAlignStore::getUnitigMicroHetProb(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  return(utgRecord[maID].mad.unitig_microhet_prob);
}

inline
UnitigStatus
MultiAlignStore::getUnitigStatus(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  return(utgRecord[maID].mad.unitig_status);
}

inline
UnitigFUR
MultiAlignStore::getUnitigFUR(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  return(utgRecord[maID].mad.unitig_unique_rept);
}

inline
ContigStatus
MultiAlignStore::getContigStatus(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)ctgLen);
  return(ctgRecord[maID].mad.contig_status);
}

inline
uint32
MultiAlignStore::getNumFrags(int32 maID, bool isUnitig) {
  assert(maID >= 0);
  return((isUnitig) ? utgRecord[maID].mad.num_frags : ctgRecord[maID].mad.num_frags);  
}

inline
uint32
MultiAlignStore::getNumUnitigs(int32 maID, bool isUnitig) {
  assert(maID >= 0);
  return((isUnitig) ? utgRecord[maID].mad.num_unitigs : ctgRecord[maID].mad.num_unitigs);
}

inline
void
MultiAlignStore::setUnitigCoverageStat(int32 maID, double cs) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  utgRecord[maID].mad.unitig_coverage_stat = cs;
  if (utgCache[maID])
    utgCache[maID]->data.unitig_coverage_stat = cs;
}

inline
void
MultiAlignStore::setUnitigMicroHetProb(int32 maID, double mp) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  utgRecord[maID].mad.unitig_microhet_prob = mp;
  if (utgCache[maID])
    utgCache[maID]->data.unitig_microhet_prob = mp;
}

inline
void
MultiAlignStore::setUnitigStatus(int32 maID, UnitigStatus status) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  utgRecord[maID].mad.unitig_status = status;
  if (utgCache[maID])
    utgCache[maID]->data.unitig_status = status;
}

inline
void
MultiAlignStore::setUnitigFUR(int32 maID, UnitigFUR fur) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  utgRecord[maID].mad.unitig_unique_rept = fur;
  if (utgCache[maID])
    utgCache[maID]->data.unitig_unique_rept = fur;
}

inline
void
MultiAlignStore::setContigStatus(int32 maID, ContigStatus status) {
  assert(maID >= 0);
  assert(maID < (int32)ctgLen);
  ctgRecord[maID].mad.contig_status = status;
  if (ctgCache[maID])
    ctgCache[maID]->data.contig_status = status;
}

inline
uint32
MultiAlignStore::getUnitigVersion(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)utgLen);
  return(utgRecord[maID].svID);
}

inline
uint32
MultiAlignStore::getContigVersion(int32 maID) {
  assert(maID >= 0);
  assert(maID < (int32)ctgLen);
  return(ctgRecord[maID].svID);
}

#endif
