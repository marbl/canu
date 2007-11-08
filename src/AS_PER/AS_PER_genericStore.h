
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

#define STREAM_FROMSTART -1
#define STREAM_UNTILEND  -1 

typedef struct{
  char          storeLabel[8];  //
  uint32        storeType;      //
  uint32        elementSize;    //
  int64         firstElem;      // Initially -1.  If >0, index of first allocated element
  int64         lastElem;       // Initially -1.  If >0, index of last allocated element

  //  Everything else is only valid when the store is loaded.

  int64         allocatedSize; //  size of that buffer
  char         *memoryBuffer;  //  Non-NULL if we have a copy of the store in memory

  FILE         *fp;            //  Non-NULL if we have a disk-backed store (can also be in memory)
  char         *diskBuffer;    //  system buffer for disk-based store

#ifdef TRUE32BIT
  char         *ptrs[3];
  int           pad;           //  64-bit wants to make this 80 bytes, not 76
#endif

  int           isDirty;       //  True if we need to flush on close
  int           readOnly;      //  True if we're a partial memory store
  int           lastWasWrite;  //  True if the last disk op was a write
} StoreStruct;

//  The "lastWasWrite" field allows us to flush the stream between
//  read and write events, as per ANSI 4.9.5.3.  It's not clear if
//  this is really needed though.  We always flush, even if there is a
//  previous seek, since our AS_UTL seek notices when the fp is
//  correctly placed, and skips the seek.
//
//  "A change of input/output direction on an update file is only
//  allowed following a fsetpos, fseek, rewind, or fflush operation,
//  since these are precisely the functions which assure that the I/O
//  buffer has been flushed."

typedef struct{
  StoreStruct  *store;
  int64         startIndex;
  int64         endIndex;
} StreamStruct;


//  Make all stores use a system buffer of size wb bytes.  Bigger
//  buffers work great for writing sequentially to lots of files
//  (e.g., partitioning a store) but are TERRIBLE for updating a store
//  (nearly everything else).
//
void          AS_PER_setBufferSize(int wb);

StoreStruct  *openStore(const char *StorePath, const char *rw);
void          closeStore(StoreStruct *sh);

static int64  getLastElemStore(StoreStruct *s)  { return(s->lastElem);  }
static int64  getFirstElemStore(StoreStruct *s) { return(s->firstElem); }

StoreStruct  *createIndexStore(const char *StorePath, const char *storeType, int32 elementSize, int64 firstID);
void          getIndexStore(StoreStruct *fs, int64 indx, void *buffer);
void         *getIndexStorePtr(StoreStruct *fs, int64 indx);
void          setIndexStore(StoreStruct *store, int64 indx, void *element);
void          appendIndexStore(StoreStruct *store, void *element);

StoreStruct  *createStringStore(const char *StorePath, const char *storeType);
void          getStringStore(StoreStruct *s, int64 offset, void *buffer, uint32 maxLength, uint32 *actualLength);
char         *getStringStorePtr(StoreStruct *s, int64 offset, uint32 *actualLength);
void          appendStringStore(StoreStruct *store, void *element, uint32 length);

StreamStruct *openStream(StoreStruct *sh);
void          resetStream(StreamStruct *sh, int64 startIndex, int64 endIndex);
int           nextStream(StreamStruct *sh, void *buffer, uint32 maxLength, uint32 *actualLength);
void          closeStream(StreamStruct *sh);


//  "Convert" the loadStore into a new memory store.  The loadStore is
//  closed.  A partial memory store does not allow writes.
//
StoreStruct *convertStoreToMemoryStore(StoreStruct *source);
StoreStruct *convertStoreToPartialMemoryStore(StoreStruct *loadStore, int64 firstElem, int64 lastElem);


//  Open an existing file-based Store, and load a portion of its
//  contents into a newly created memory-based Store.
///
static
StoreStruct *loadStore(const char *StorePath) {
  return(convertStoreToMemoryStore(openStore(StorePath, "r")));
}

static
StoreStruct *loadStorePartial(const char *StorePath, int64 firstElem, int64 lastElem) {
  return(convertStoreToPartialMemoryStore(openStore(StorePath, "r"), firstElem, lastElem));
}


#endif /* AS_PER_GENERICSTORE_H */
