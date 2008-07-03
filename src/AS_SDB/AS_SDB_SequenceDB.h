
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

#ifndef AS_SEB_SEQUENCEDB_H
#define AS_SEB_SEQUENCEDB_H

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"

#include "MultiAlignStore_CNS.h"

//  File structure is:
//    utg -> tMARecord Unitigs
//    ctg -> tMARecord Contigs
//    dat -> MultiAlignT data for both unitigs and contigs


typedef struct {
  uint64       storeID:12;
  uint64       multiAlignID:28;
  uint64       isDeleted:1;
  off_t        offset;   // offset in store from which to load this multi-alignment
} tMARecord;
VA_DEF(tMARecord)

typedef struct {
  char                *path;
  int32                currentRevision;

  VA_TYPE(tMARecord)  *Unitigs;
  VA_TYPE(tMARecord)  *Contigs;

  MultiAlignStoreT    *UnitigStore;  // cache of loaded MultiAligns
  MultiAlignStoreT    *ContigStore;

  //  VA's of pointers are awkward, and we have a reasonably good feel
  //  for the upper limit here, since we want to keep files open.
  //
  int32                dataFileLen;
  int32                dataFileMax;
  FILE               **dataFile;

  //  Data for partitioned seqDB
  VA_TYPE(tMARecord)  *multiAligns;
  HashTable_AS        *multiAlignLookup;
} tSequenceDB;


tSequenceDB *createSequenceDB(char *path);
tSequenceDB *openSequenceDB(char *path, int readWrite, int revision);
void         openSequenceDBPartition(tSequenceDB *db, int32 partition);
void         saveSequenceDB(tSequenceDB *db);  // Save the current revision of the indices
void         deleteSequenceDB(tSequenceDB *db);



//  update() requires that the thing already be there.
//
//  insert() requires that the thing isn't there, or, it was there,
//  but has been deleted.
//
//  If keepInCache, we keep a pointer to the MultiAlignT.  THE STORE
//  NOW OWNS THE MULTIALIGN.
//
void         updateMultiAlignTInSequenceDB(tSequenceDB *db, int index, int isUnitig, MultiAlignT *ma, int keepInCache);
void         insertMultiAlignTInSequenceDB(tSequenceDB *db, int index, int isUnitig, MultiAlignT *ma, int keepInCache);

//  Remove a multialign from the sequenceDB, and delete any cached
//  copies.  You should have no pointers to that multialign.
//
void         deleteMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig);

//  load() will load and cache the multialign.  THE STORE OWNS THIS MULTIALIGN.
//  DO NOT DeleteMultiAlignT() it.
//
//  copy() will load and copy the multialign into the provided ma.  It
//  will not cache, YOU OWN THE MULTIALIGN.
//
MultiAlignT *loadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig);
void         copyMultiAlignTFromSequenceDB(tSequenceDB *db, MultiAlignT *ma, int index, int isUnitig);


static
void
clearCacheSequenceDB(tSequenceDB *db) {
  fflush(db->dataFile[db->currentRevision]);
  fprintf(stderr, "clearCacheSequenceDB() disabled.\n");
  return;
  ClearMultiAlignStoreT(db->UnitigStore);
  ClearMultiAlignStoreT(db->ContigStore);
}

static
off_t
getSizeOfCurrentSequenceDB(tSequenceDB *db) {
  return(AS_UTL_ftell(db->dataFile[db->currentRevision]));
}

static
int32
numUnitigsInSequenceDB(tSequenceDB *db){
  return GetNumtMARecords(db->Unitigs);
}

static
int32
numContigsInSequenceDB(tSequenceDB *db){
  return GetNumtMARecords(db->Contigs);
}

#endif
