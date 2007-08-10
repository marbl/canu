
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

#ifndef AS_PER_GENERICSTORE_H
#define AS_PER_GENERICSTORE_H


#ifdef GENERIC_STORE_USE_LONG_STRINGS
#define VLSTRING_SIZE_T uint32
#define F_VLS F_U32
#define VLSTRING_MAX_SIZE (4 * 1024 * 1024 - 1)
#else
#define VLSTRING_SIZE_T uint16
#define F_VLS F_U16
#define VLSTRING_MAX_SIZE (64 * 1024 - 1)
#endif

#include <time.h>

typedef enum { UnAllocatedStore = 0, 
	       UnInitializedStore, 
	       ActiveStore} StoreStatus;

#define INVALID_STORE  0
#define INDEX_STORE    1
#define VLRECORD_STORE 2

#define STREAM_UNTILEND  -1 
#define STREAM_FROMSTART -1

typedef struct{
  uint32        isDeleted:1;
  uint32        type:3;
  uint32        p1:28;          // padding field
  uint32        p2:32;          // padding field
  char          storeType[8];
  int64         firstElem;      // Initially -1.  If >0, index of first allocated element
  int64         lastElem;       // Initially -1.  If >0, index of last allocated element
  int32         unused_version; // Was the 'version', now just to keep the struct the same size
  int32         elementSize;  
  int64         creationTime;
  int64         lastUpdateTime;
} StoreStat;

typedef struct{
  FILE        *fp;            //  For a file-based store
  char        *buffer;        //  For a memory-based store, also holds setbuffer() buffer for disk stores
  int64        allocatedSize;
  StoreStatus  status; 
  StoreStat    header;
  int64        lastCommittedElem;  //  Initially -1.  If >0, index of last committed element
  int          isMemoryStore;
  int          isDirty;
  //  The "lastWasRead" field allows us to flush the stream between read
  //  and write events, as per ANSI 4.9.5.3.  It's not clear if this is
  //  really needed though.
  //
  int          lastWasRead;
} StoreStruct;

typedef struct{
  StoreStruct  *store;
  void         *buffer;
  int32         bufferSize;
  int64         startIndex;
  int64         endIndex;
  StoreStatus   status; 
} StreamStruct;



//  Make all stores use a system buffer of size wb bytes.  Bigger
//  buffers work great for writing sequentially to lots of files
//  (e.g., partitioning a store) but are TERRIBLE for updating a store
//  (nearly everything else).
//
void          AS_PER_setBufferSize(int wb);

StoreStruct  *openStore(const char *StorePath, const char *rw);
int           statsStore(StoreStruct *store, StoreStat *stats);
int           closeStore(StoreStruct *sh);

int64         getLastElemStore(StoreStruct *store);
int64         getFirstElemStore(StoreStruct *store);

StoreStruct  *createIndexStore(const char *StorePath, const char *storeType, int32 elementSize, int64 firstID);
int           getIndexStore(StoreStruct *fs, int64 indx, void *buffer);
void         *getIndexStorePtr(StoreStruct *fs, int64 indx);
int           setIndexStore(StoreStruct *store, int64 indx, void *element);
int           appendIndexStore(StoreStruct *store, void *element);

StoreStruct  *createVLRecordStore(const char *StorePath, const char *storeType, int32 expectedRecordSize);
int           getVLRecordStore(StoreStruct *s, int64 offset, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength);
int           appendVLRecordStore(StoreStruct *store, void *element, VLSTRING_SIZE_T length);

StreamStruct *openStream(StoreStruct *sh, void *buffer, int32 bufferSize);
int           resetStream(StreamStruct *sh, int64 startIndex, int64 endIndex);
int           nextStream(StreamStruct *sh, void *buffer);
int           closeStream(StreamStruct *sh);

int           nextVLRecordStream(StreamStruct *sh, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength);


//  "Convert" the loadStore into a new memory store.  The loadStore is
//  closed.
//
StoreStruct *
convertStoreToPartialMemoryStore(StoreStruct *loadStore, int64 firstElem, int64 lastElem);


//  Open an existing file-based Store, and load a portion of its
//  contents into A newly created memory-based Store.
///
static
StoreStruct *
loadStorePartial(const char *StorePath, int64 firstElem, int64 lastElem) {
  return(convertStoreToPartialMemoryStore(openStore(StorePath, "r"), firstElem, lastElem));
}


#endif /* AS_PER_GENERICSTORE_H */
