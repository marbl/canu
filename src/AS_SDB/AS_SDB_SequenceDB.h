
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
/* 	$Id: AS_SDB_SequenceDB.h,v 1.4 2005-03-22 19:49:26 jason_miller Exp $	 */
/*
  This SequenceDB is a more sophisticated version of the MultiAlignStore idea.
  We have the following design requirements:
  1) Keep unitig and contig multi-alignments on disk, loading into memory as needed
  2) For contigs containing a single unitig, simply reference the unitig multi-alignment
  3) Incrementally generate multi-alignments as needed, storing the the new multi-alignments in
     a new file.  This reduces the I/O requirements of cgw, and makes it difficult to corrupt stores,
     since most are opened read only.
  4) Compress the sequence and quality of the multi-alignments to consume a single byte per seq/quality pair.
     This is accomplished by upgrading the streaming of MultiAlignT data structure.
  5) Maintain a limited cache of multi-alignments, and flush as required to meet memory budget

  The stores in this SequenceDB are names:
      SequenceDB.index.i    // The cumulative index of stores 0-i
      SequenceDB.data.i

      Where i is the index of the substore, starting from 0.
      The 0-th substore holds Unitigs only.
      The 1-st substore holds the contigs produced the first set of contigs.
*/


#ifndef AS_SEB_SEQUENCEDB_H
#define AS_SEB_SEQUENCEDB_H

#include "MultiAlignStore_CNS.h"

/* MARecord records the location of each Unitig/Contig in the Assembly */
typedef struct MARecord_tag {
  int32 storeID;   // in a sequenceDB this is the index of a file, in a sequenceDBPartition it is overloaded as the id of a multiAlignment
  
   union{
   struct{
     unsigned int deleted:1;
     unsigned int spare:31;
   }bits;
     int32 all;
   }flags;
  off_t offset;   // offset in store from which to load this multi-alignment
}tMARecord;


VA_DEF(tMARecord)


typedef struct SequenceDB_tag{
  char *path;      
#ifdef i386
  int32 ptrPad1;
#endif

  VA_TYPE(tMARecord) *Unitigs;
#ifdef i386
  int32 ptrPad2;
#endif

  VA_TYPE(tMARecord) *Contigs;
#ifdef i386
  int32 ptrPad3;
#endif

  // The following are caches of MultiAlignTs.
  MultiAlignStoreT *UnitigStore;
#ifdef i386
  int32 ptrPad4;
#endif

  MultiAlignStoreT *ContigStore;
#ifdef i386
  int32 ptrPad5;
#endif

  // All but the last store are opened read only
  // The last store is opened r/w
  VA_TYPE(PtrT) *SubStores;
#ifdef i386
  int32 ptrPad6;
#endif

  int32 currentRevision;
  int positionedAtEnd;

  size_t totalCacheSize;
#ifdef i386
  int32 sztPad1;
#endif

  size_t totalCacheLimit;
#ifdef i386
  int32 sztPad2;
#endif

  off_t offsetOfEOF;
}tSequenceDB;


/* CreateSequenceDB:
   Create a SequenceDB with a single data file, and the index in memory..
   If force = true, this will remove an existing store with the same path.
   Fails if cannot open both the data and index files.
*/
tSequenceDB *CreateSequenceDB(char *path, int initialSize, int force);

/* DeleteSequenceDB:
   Destructor.
*/
void DeleteSequenceDB(tSequenceDB *db);

static int32 GetRevisionSequenceDB(tSequenceDB *db){
  return db->currentRevision;
}


/* OpenSequenceDB:
   Opens an existing SequenceDB.  The index file that will be opened corresponds
   to revision.  If this does not exist, failure.  If successful, the files
   corresponding to revisions 0,revision are opened read-only, and the files
   for revision+1 is opened r/w.
*/
tSequenceDB *OpenSequenceDB(char *path, int readWrite, int revision);

/* SaveSequenceDB
   Save the current revision of the indices.  Used by the checkpointing code.
   The indicies are maintained in memory.
*/
void SaveSequenceDB(tSequenceDB *db);  // Save the current revision of the indices

/* InsertMultiAlignTInSequenceDB
   Inserts a new MultiAlignT into the appropriate cache and indices and
   appends it to the data file.
*/
void InsertMultiAlignTInSequenceDB(tSequenceDB *db, int index, int isUnitig, MultiAlignT *ma, int keepInCache);

/* UpdateMultiAlignTInSequenceDB
   Inserts a new MultiAlignT for a given chunk that will replace an old MultiAlignT.  The old data is actually
   left on disk, but the substore and offset indicators for the chunk index portion of the (latest version of)
   the SDB will be updated to point to the new ma.
*/
void UpdateMultiAlignTInSequenceDB(tSequenceDB *db, int index, int isUnitig, MultiAlignT *ma, int keepInCache);

/* DuplicateEntryInSequenceDB
*/
void DuplicateEntryInSequenceDB(tSequenceDB *db, int fromIndex, int fromIsUnitig, int toIndex, int toIsUnitig, int keepInCache);

/* DeleteMultiAlignTFromSequenceDB
   Mark the appropriate entry as deleted in an index.  This does not change the data file
*/
void DeleteMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig);

/* LoadMultiAlignTFromSequenceDB
   Loads a MultiAlignT into the cache (unless it is already present.
   On Failure, returns NULL.
*/
MultiAlignT *LoadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig);

/* ReLoadMultiAlignTFromSequenceDB
   Loads a MultiAlignT into an existing ma
   Does not load cache
   Returns 0 on success, returns 1 if index is deleted or not present.
*/
int32 ReLoadMultiAlignTFromSequenceDB(tSequenceDB *db, MultiAlignT *ma, int index, int isUnitig);

/* UnLoadMultiAlignTFromSequenceDB
   Unloads a MultiAlignT from a cache.
*/
void UnloadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig);

/* ClearCacheSequenceDB
   Clears all the entries from a cache.
*/
void ClearCacheSequenceDB(tSequenceDB *db, int isUnitig);

/* ClearCacheSequenceDBConditionally
   Clears all the entries from a cache if it is bigger than a certain size
*/
void ClearCacheSequenceDBConditionally(tSequenceDB *db, size_t maxSize);

/* Return the size of the current data file being written for the sequence DB */
size_t GetSizeOfCurrentSequenceDB(tSequenceDB *db);

static int32 NumUnitigsInSequenceDB(tSequenceDB *db){
  return GetNumtMARecords(db->Unitigs);
}

static int32 NumContigsInSequenceDB(tSequenceDB *db){
  return GetNumtMARecords(db->Contigs);
}


#endif
