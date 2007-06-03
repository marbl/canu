
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

static char CM_ID[] = "$Id: AS_PER_genericStore.c,v 1.19 2007-06-03 08:13:22 brianwalenz Exp $";

// Module:  AS_PER_genericStore
// Description:
//
//     This module defines the interface and implementation of the index
//     and string "stores" -- binary files used to store fixed and
//     variable size data for the assembler.  An additional type of
//     store needs to be implemented, for storing variable length
//     records.  Strings should probably be implemented as a variable
//     length record prefixed by its length.  The Fragment Store is
//     built using the building blocks of the index and string stores.
//     Both types of stores provide both a file-based and memory-based
//     implementations.  The intent is for the client to build up
//     buffers of records/strings in a memory-based store, and
//     concatenate them to the file-based store using the (as yet
//     unimplemented) concat operation.  A store consists of a fixed
//     length header, followed by data.  The only types of modifications
//     to the store that are currently permitted are to append a new
//     record/string to the store, or to mark a record in an indexStore
//     as deleted.  In principle, there is no reason why the index store
//     could not support a replace operation, that would replace a
//     record with new data, as long as the index of the record in
//     question was within range.
//
//     Client code relates to a store through an opaque handle.  This
//     will faciliate changes to the store structure as we learn more
//     about requirements and optimization.
//
//     Each type of store supports a "stream" operation, for read access
//     to successive elements of the store.
//
// Assumptions:
//
//      To support the delete operation on an index store, the Most
//      Significant Bit of the stored data record must be a deleted
//      flag.  Currently, only the client pays attention to the deleted
//      bit.
//
// Document:
//      GenericStore.rtf

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_fileIO.h"


//  Sigh.  As much as BPW likes the BSD, he is rather perturbed that
//  this is needed.  When adding LKG's to the gkpStore, which involves
//  lots of random access IO, FreeBSD 6.1 (and, IIRC, the 4.x series)
//  _occasionally_ lose a couple bytes here and there, during the
//  writes.  Shouldn't be too hard to figure out based on what this
//  define enables.
//
//  What is harder to figure out, is why the fix enabled by this
//  works.
//
//#ifdef __FreeBSD__
//#warning FreeBSD random access IO is BROKEN.  Performance degraded.
//#define CHECK_SET_INDEX_STORE
//#endif


//  This is the structure maintained in memory for each open Store
//
//  The "lastWasRead" field allows us to flush the stream between read
//  and write events, as per ANSI 4.9.5.3.  It's not clear if this is
//  really needed though.
//
typedef struct{
  FILE        *fp;            //  For a file-based store
  char        *buffer;        //  For a memory-based store, also holds setbuffer() buffer for disk stores
  int64        allocatedSize;
  StoreStatus  status; 
  StoreStat    header;
  int64        lastCommittedElem;  //  Initially -1.  If >0, index of last committed element
  int          isMemoryStore;
  int          isDirty;
  int          lastWasRead;
}StoreStruct;


//  Utility for writing a store Header
//
static void
writeHeader(StoreStruct *myStore){

  myStore->header.lastUpdateTime = time(0);

  if(myStore->isMemoryStore){
    memcpy(myStore->buffer,&myStore->header, sizeof(StoreStat));
  }else{
    rewind(myStore->fp);
    AS_UTL_safeWrite(myStore->fp,(void *)&myStore->header,"writeHeader",sizeof(StoreStat), 1);
  }
}



static StoreStruct *gStores = NULL;
static int gNumStores       = 0;
static int gMaxNumStores    = 0;

static StreamStruct *gStreams = NULL;
static int gNumStreams        = 0;
static int gMaxNumStreams     = 0;

static int INITIAL_ALLOCATION = 4096;
static int WRITING_BUFFER     = 8 * 1024;


#define  Store_isMemoryStore(S)    ((S)->isMemoryStore)
#define  Store_isDirty(S)          ((S)->isDirty)
#define  Store_setType(S, T)       ((S)->header.type = (T))
#define  Store_getType(S)          ((S)->header.type)
#define  Store_isIndexStore(S)     ((S)->header.type == INDEX_STORE)
#define  Store_isStringStore(S)    ((S)->header.type == STRING_STORE)
#define  Store_isVLRecordStore(S)  ((S)->header.type == VLRECORD_STORE)

#define  Store_myHandle(S)  ((S) - gStores)
#define  Stream_myHandle(S) ((S) - gStreams)

#define  Store_myStruct(S)  (((S) < 0 || (S) > gMaxNumStores)  ? NULL : (gStores + (S)))
#define  Stream_myStruct(S) (((S) < 0 || (S) > gMaxNumStreams) ? NULL : (gStreams + (S)))

#define  computeOffset(S, I)  (((I) - (S)->header.firstElem) * (int64)((S)->header.elementSize) + sizeof(StoreStat))




void AS_PER_setBufferSize(int wb) {
  WRITING_BUFFER = wb;
}


StoreStruct  *allocateStore(void){
  int i;
  StoreStruct *ss = NULL;

  if (gStores == NULL) {
    gMaxNumStores = 64;
    gStores       = (StoreStruct *)safe_calloc(gMaxNumStores, sizeof(StoreStruct));
    for  (i=gNumStores;  i<gMaxNumStores;  i++)
      gStores[i].status = UnAllocatedStore;
  }

  for (i=0; i<gMaxNumStores; i++) {
    if (gStores[i].status == UnAllocatedStore) {
      ss = gStores + i;
      break;
    }
  }

  if (ss == NULL) {
    assert(gNumStores >= gMaxNumStores);

    gMaxNumStores *= 2;
    gStores = (StoreStruct *) safe_realloc(gStores, gMaxNumStores * sizeof(StoreStruct));
    for  (i=gNumStores;  i<gMaxNumStores;  i++)
      gStores[i].status = UnAllocatedStore;

    ss = gStores + gNumStores;
  }

  ss->fp                    = NULL;
  ss->buffer                = NULL;
  ss->allocatedSize         = 0;
  ss->status                = UnInitializedStore;
  ss->header.firstElem      = -1;
  ss->header.lastElem       = -1;
  ss->lastCommittedElem     = -1;

  ss->header.isDeleted      = 0;
  ss->header.type           = 0;
  ss->header.p1             = 0;
  ss->header.p2             = 0;
  ss->header.storeType[0]   = 0;
  ss->header.storeType[1]   = 0;
  ss->header.storeType[2]   = 0;
  ss->header.storeType[3]   = 0;
  ss->header.storeType[4]   = 0;
  ss->header.storeType[5]   = 0;
  ss->header.storeType[6]   = 0;
  ss->header.storeType[7]   = 0;
  ss->header.firstElem      = -1;
  ss->header.lastElem       = -1;
  ss->header.version        = -1;
  ss->header.elementSize    = -1;
  ss->header.creationTime   = 0;
  ss->header.lastUpdateTime = 0;

  ss->lastCommittedElem     = -1;
  ss->isMemoryStore         = 0;
  ss->isDirty               = 0;
  ss->lastWasRead           = 1;

  gNumStores++;
  return(ss);
}



StreamStruct  *allocateStream(StoreHandle s, void *buffer, int32 bufferSize  ){
  int i;
  StreamStruct *ss = NULL;

  if (gStreams == NULL) {
    gMaxNumStreams = 64;
    gStreams       = (StreamStruct *)safe_calloc(gMaxNumStreams, sizeof(StreamStruct));
    for  (i=gNumStreams;  i<gMaxNumStreams;  i++)
      gStreams[i].status = UnAllocatedStore;
  }

  for (i=0; i<gMaxNumStreams; i++) {
    if (gStreams[i].status == UnAllocatedStore) {
      ss = gStreams + i;
      break;
    }
  }

  if (ss == NULL) {
    assert(gNumStreams >= gMaxNumStreams);

    gMaxNumStreams *= 2;
    gStreams = (StreamStruct *) safe_realloc(gStreams, gMaxNumStreams * sizeof(StreamStruct));
    for  (i=gNumStreams;  i<gMaxNumStreams;  i++)
      gStreams[i].status = UnAllocatedStore;

    ss = gStreams + gNumStreams;
  }

  ss->store       = s;
  ss->status      = ActiveStore;
  ss->buffer      = buffer;
  ss->bufferSize  = bufferSize;
  ss->startIndex  = gStores[s].header.firstElem;
  ss->endIndex    = gStores[s].header.lastElem;

  gNumStreams++;
  return(ss);
}


StreamHandle openStream(StoreHandle fs,/* handle to a fragment store */
                        void *buffer,  /* User supplied buffer for prefetching */
                        int32 bufferSize) {
  StreamStruct *stream = allocateStream(fs,buffer,bufferSize);
  return Stream_myHandle(stream);
}


int resetStream(StreamHandle sh, int64 startIndex, int64 endIndex){
  StreamStruct *ss = Stream_myStruct(sh);
  StoreHandle s = ss->store;

  if(startIndex == STREAM_FROMSTART)
    ss->startIndex = gStores[s].header.firstElem;
  else
    ss->startIndex = startIndex;
  if(endIndex == STREAM_UNTILEND)
    ss->endIndex = gStores[s].header.lastElem;
  else
    ss->endIndex = endIndex;

  return(0);
}


StreamHandle openStringStream(StoreHandle fs, /* handle to a fragment store */
                              int64 startOffset, /* Offset to start streaming */
                              void *buffer, /* User supplied buffer for prefetching */
                              int32 bufferSize) {
  StreamStruct *stream = allocateStream(fs,buffer,bufferSize);
  resetStream(Stream_myHandle(stream), startOffset, STREAM_UNTILEND);
  return Stream_myHandle(stream);
}


int closeStream(StreamHandle sh){
  StreamStruct *stream = Stream_myStruct(sh);

  stream->status = UnAllocatedStore;
  gNumStreams--;
  return 0;
}


StoreHandle createIndexStore(const char *path, const char *storeType, int32 elementSize, 
                             int32 version, int64 firstID) {
  StoreStruct *myStore;
  myStore = allocateStore();

  AssertPtr(myStore);

  Store_setType(myStore,INDEX_STORE);
  strncpy(myStore->header.storeType, storeType,7);
  myStore->header.elementSize = elementSize;
  myStore->header.version = version;
  myStore->header.firstElem = firstID;
  myStore->header.lastElem =  firstID -1; /* Sentinel */
  myStore->lastCommittedElem = -1;
  myStore->status = UnInitializedStore;
  myStore->header.creationTime = time(0);
  myStore->status = ActiveStore;

  if(path){ /* File based store */
    assert(strlen(path) < FILENAME_MAX);
    errno = 0;
    myStore->fp = fopen(path, "w+");
    if(errno){
      fprintf(stderr,"createIndexStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
      assert(errno == 0);
    }
    myStore->buffer = NULL;
    if (WRITING_BUFFER > 0) {
      //fprintf(stderr, "createIndexStore()--  Allocating fp buffer of "F_SIZE_T" bytes for '%s'.\n", WRITING_BUFFER, path);
      myStore->buffer = (char *)safe_malloc(WRITING_BUFFER);
      setbuffer(myStore->fp, myStore->buffer, WRITING_BUFFER);
    }
    myStore->isMemoryStore = 0;
  }else{
    myStore->allocatedSize = INITIAL_ALLOCATION * elementSize;
    myStore->buffer = (char *)safe_malloc(myStore->allocatedSize);
    myStore->isMemoryStore = 1;
    myStore->fp = NULL;
  }
  writeHeader(myStore); 
  myStore->isDirty = 0;
  return Store_myHandle(myStore);
}

StoreHandle resetIndexStore( StoreHandle sh, int64 firstID){
  StoreStruct *myStore = Store_myStruct(sh);

  AssertPtr(myStore);
  assert(myStore->status == ActiveStore);
  myStore->header.firstElem = firstID;
  myStore->header.lastElem =  firstID -1; /* Sentinel */
  myStore->lastCommittedElem = -1;
  myStore->status = UnInitializedStore;
  myStore->header.creationTime = time(0);
  myStore->status = ActiveStore;

  writeHeader(myStore); 
  myStore->isDirty = 0;
  return Store_myHandle(myStore);
}


StoreHandle createStringStore(const char *path, const char *storeType, 
                              int32 expectedStringSize, int32 version) {
  StoreStruct *myStore;
  myStore = allocateStore();

  AssertPtr(myStore);

  Store_setType(myStore,STRING_STORE);

  strncpy(myStore->header.storeType, storeType,7);
  myStore->header.elementSize = expectedStringSize;
  myStore->header.version = version;
  myStore->header.firstElem = 0;
  myStore->header.lastElem =  0; /* We use this to track total length */
  myStore->lastCommittedElem = -1;
  myStore->status = UnInitializedStore;
  myStore->header.creationTime = time(0);

  if(path){ /* File based store */
    myStore->isMemoryStore = 0;
    assert(strlen(path) < FILENAME_MAX);
    myStore->isMemoryStore = 0;
    errno = 0;
    myStore->fp = fopen(path,"w+");
    if(errno){
      fprintf(stderr,"createStringStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
      assert(errno == 0);
    }
    myStore->buffer = NULL;
    if (WRITING_BUFFER > 0) {
      //fprintf(stderr, "createStringStore()--  Allocating fp buffer of "F_SIZE_T" bytes for '%s'.\n", WRITING_BUFFER, path);
      myStore->buffer = (char *)safe_malloc(WRITING_BUFFER);
      setbuffer(myStore->fp, myStore->buffer, WRITING_BUFFER);
    }
  }else{
    myStore->allocatedSize = expectedStringSize;
    myStore->buffer = (char *)safe_malloc(myStore->allocatedSize);
    AssertPtr(myStore->buffer);
    myStore->isMemoryStore = 1;
    myStore->fp = NULL;
  }

  myStore->status = ActiveStore;

  writeHeader(myStore); 

  return Store_myHandle(myStore);

}


StoreHandle resetStringStore( StoreHandle sh){
  StoreStruct *myStore = Store_myStruct(sh);

  AssertPtr(myStore);
  assert(myStore->status == ActiveStore);
  myStore->header.firstElem = 0;
  myStore->header.lastElem =  0; /* We use this to track total length */
  myStore->lastCommittedElem = -1;
  myStore->status = UnInitializedStore;
  myStore->header.creationTime = time(0);
  myStore->status = ActiveStore;

  writeHeader(myStore); 
  myStore->isDirty = 0;
  return Store_myHandle(myStore);
}


StoreHandle createVLRecordStore(const char *path, const char *storeType, 
                                int32 expectedVLRecordSize, int32 version) {
  StoreHandle newStore = createStringStore(path, storeType, expectedVLRecordSize, VLRECORDSTORE_VERSION);
  StoreStruct *s = Store_myStruct(newStore);

  AssertPtr(s);
  Store_setType(s,VLRECORD_STORE);

  return newStore;

}


StoreHandle resetVLRecordStore( StoreHandle sh){
  return resetStringStore(sh);
}


int closeStore(StoreHandle s){
  StoreStruct *myStore = Store_myStruct(s);

  assert(myStore->status == ActiveStore);

  if(myStore->isMemoryStore){
    myStore->allocatedSize = 0;
    safe_free(myStore->buffer);
  }else{
    if(myStore->isDirty)
      writeHeader(myStore);
    if(fclose(myStore->fp) != 0)
      assert(0);
    safe_free(myStore->buffer);
  }
  myStore->status = UnAllocatedStore;
  return 0;
}


int setIndexStore(StoreHandle s, int64 index, void *element){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  //if (index != myStore->header.lastElem)
  //  fprintf(stderr, "setIndexStore()--  Write index %d of size %d\n", (int)index, myStore->header.elementSize);

  assert(myStore->status == ActiveStore);
  assert(myStore->header.firstElem <= index);
  assert(myStore->header.lastElem >= index);

  myStore->isDirty = 1;

  offset = computeOffset(myStore, index);

  if(myStore->isMemoryStore){
    memcpy(myStore->buffer + offset, element, myStore->header.elementSize);
  }else{

    if (myStore->lastWasRead)
      fflush(myStore->fp);
    myStore->lastWasRead = 0;

    AS_UTL_fseek(myStore->fp, (off_t) offset, SEEK_SET);
    AS_UTL_safeWrite(myStore->fp, element, "setIndexStore", myStore->header.elementSize, 1);

#ifdef CHECK_SET_INDEX_STORE
    {
      int    i;
      char  *checkbuffer = (char *)safe_malloc(sizeof(char) * myStore->header.elementSize);
      for (i=0; i<myStore->header.elementSize; i++)
        checkbuffer[i] = 0xff;

      fflush(myStore->fp);
      AS_UTL_fseek(myStore->fp, (off_t) offset, SEEK_SET);
      if (1 != AS_UTL_safeRead(myStore->fp, checkbuffer, "setIndexStore", myStore->header.elementSize, 1))
        assert(0);
      fflush(myStore->fp);
      if (memcmp(element, checkbuffer, myStore->header.elementSize) != 0) {
        fprintf(stderr, "failed on index %d\n", (int)index);
        fprintf(stderr, "element\n");
        for (i=0; i<myStore->header.elementSize; i++)
          fprintf(stderr, "%02x", ((char *)(element))[i]);
        fprintf(stderr, "\n");
        fprintf(stderr, "readback\n");
        for (i=0; i<myStore->header.elementSize; i++)
          fprintf(stderr, "%02x", checkbuffer[i]);
        fprintf(stderr, "\n");
        assert(0);
      }

      safe_free(checkbuffer);
    }
#endif

  }

  return(0);
}


int appendIndexStore(StoreHandle s, void *element){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  assert(myStore->status == ActiveStore);
  myStore->isDirty = 1;

  offset = computeOffset(myStore, myStore->header.lastElem + 1);

  if (myStore->isMemoryStore){
    int64 desiredSize = offset + (myStore->header.lastElem - myStore->header.firstElem + 1) * myStore->header.elementSize;

    if (myStore->allocatedSize <= desiredSize) {
      while(myStore->allocatedSize <= desiredSize)
	myStore->allocatedSize *= 2;
      myStore->buffer = (char *)safe_realloc(myStore->buffer, myStore->allocatedSize);
    }
  }

  myStore->header.lastElem++;

  setIndexStore(s, myStore->header.lastElem, element);

  return(0);
}


StoreHandle openStore(const char *path, const char *rw) {
  StoreStruct *myStore = allocateStore();
  errno = 0;
  myStore->fp = fopen(path,rw);
  if(errno){
    fprintf(stderr,"openStore()-- Failure opening Store %s for %s: %s\n", path, rw, strerror(errno));
    return(NULLSTOREHANDLE);
  }

  myStore->status = ActiveStore;

  if (1 != AS_UTL_safeRead(myStore->fp, (void *)&myStore->header, "openStore", sizeof(StoreStat), 1))
    assert(0);

  return myStore -gStores;
}


StoreHandle convertStoreToPartialMemoryStore(StoreHandle loadStoreHandle,
                                             int64       firstElem,
                                             int64       lastElem) {
  StoreStruct *source       = Store_myStruct(loadStoreHandle);
  StoreStruct *target         = NULL;
  int64        sourceOffset    = 0;
  int64        sourceMaxOffset = 0;

  if (loadStoreHandle == NULLSTOREHANDLE)
    return(NULLSTOREHANDLE);

  if (firstElem <= 0)
    firstElem = source->header.firstElem;

  if (lastElem <= 0)
    lastElem = source->header.lastElem;

  if (Store_isStringStore(source)) {
    sourceOffset    = sizeof(StoreStat) + firstElem;
    sourceMaxOffset = sizeof(StoreStat) + lastElem;
  }else{
    sourceOffset    = computeOffset(source, firstElem);
    sourceMaxOffset = computeOffset(source, lastElem + 1);
  }

  target  = allocateStore();
  *target = *source;

  target->allocatedSize    = sizeof(StoreStat) + sourceMaxOffset - sourceOffset;
  target->buffer           = (char *)safe_malloc(target->allocatedSize);
  target->isMemoryStore    = 1;
  target->fp               = NULL;
  target->header.firstElem = firstElem;
  target->header.lastElem  = lastElem;

  //  Copy the store into the memory store.  The first
  //  sizeof(StoreStat) bytes are for a copy of target->header
  //  (written by writeHeader()), which is never used.
  //
  AS_UTL_fseek(source->fp, (off_t)sourceOffset, SEEK_SET);

  if (source->lastWasRead == 0)
    fflush(source->fp);
  source->lastWasRead = 1;

  fprintf(stderr, "reading %d bytes\n", (int)(sourceMaxOffset - sourceOffset));

  if (sourceMaxOffset - sourceOffset != AS_UTL_safeRead(source->fp, target->buffer + sizeof(StoreStat),
                                                        "convertStoreToPartialMemoryStore",
                                                        sizeof(char), sourceMaxOffset - sourceOffset))
    assert(0);

  closeStore(loadStoreHandle);

  return(target - gStores);
}




int getIndexStore(StoreHandle s, int64 index, void *buffer){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  assert(myStore->status == ActiveStore);
  offset = computeOffset(myStore,index);

  if(myStore->isMemoryStore){
    memcpy(buffer, myStore->buffer + offset, myStore->header.elementSize);
    return(0);
  }else{
    AS_UTL_fseek(myStore->fp, (off_t) offset, SEEK_SET);

    if (myStore->lastWasRead == 0)
      fflush(myStore->fp);
    myStore->lastWasRead = 1;

    return(1 == AS_UTL_safeRead(myStore->fp,buffer,"getIndexStore",myStore->header.elementSize, 1));
  }
}



void *getIndexStorePtr(StoreHandle s, int64 index) {
  StoreStruct *myStore = Store_myStruct(s);

  assert(myStore->status == ActiveStore);
  assert(myStore->isMemoryStore);

  return((void *)(myStore->buffer + computeOffset(myStore,index)));
}


int getStringStore(StoreHandle s, int64 offset, char *buffer, int32 maxLength){
  StoreStruct *myStore = Store_myStruct(s);
  int32 length = maxLength;
  assert(Store_isStringStore(myStore));
  if(offset > myStore->header.lastElem){
    fprintf(stderr," ERROR: StringStoreGet at offset "F_S64"...lastElem = "F_S64"\n",
	    offset, myStore->header.lastElem);
  }
  if(offset + length > myStore->header.lastElem){
    length = myStore->header.lastElem - offset;
  }
  
  if(myStore->isMemoryStore){
    strncpy(buffer, myStore->buffer + offset + sizeof(StoreStat), length);
  }else{
    AS_UTL_fseek(myStore->fp, (off_t) (offset + sizeof(StoreStat)), SEEK_SET);
    
    if (myStore->lastWasRead == 0)
      fflush(myStore->fp);
    myStore->lastWasRead = 1;

    if (length != AS_UTL_safeRead(myStore->fp,buffer,"getStringStore", sizeof(char), length))
      assert(0);
  }
  return(0);
}


int getVLRecordStore(StoreHandle s, int64 offset, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength){
  StoreStruct *myStore = Store_myStruct(s);
  VLSTRING_SIZE_T length = 0;
  int64 actualOffset;

  assert(Store_isVLRecordStore(myStore));
  if(offset > myStore->header.lastElem){
    fprintf(stderr," ERROR: VLRecordStoreGet at offset "F_S64"...lastElem = "F_S64"\n",
	    offset, myStore->header.lastElem);
    return(1);
  }
  
  actualOffset = offset - myStore->header.firstElem;

  if(myStore->isMemoryStore){
    memcpy(&length,(myStore->buffer + actualOffset + sizeof(StoreStat)),
           sizeof(VLSTRING_SIZE_T));
    if(!(length + actualOffset <= myStore->header.lastElem &&
	 length <= maxLength)){
      fprintf(stderr,
              "getVLRecordStore(memory) FAILURE at offset "F_S64" actualOffset "F_S64"\n\t"
	      "length = " F_VLS " implied max offset = "F_S64"  "
              "lastElem = "F_S64" maxLength = " F_VLS "\n",
	      offset, actualOffset,
	      length, length + actualOffset,
              myStore->header.lastElem, maxLength);
      assert(0);
    }
    if(length > 0)
      memcpy(buffer, myStore->buffer + actualOffset +
             sizeof(StoreStat) + sizeof(VLSTRING_SIZE_T), length);

    
  }else{
    AS_UTL_fseek(myStore->fp, (off_t) (actualOffset + sizeof(StoreStat)), SEEK_SET);

    if (myStore->lastWasRead == 0)
      fflush(myStore->fp);
    myStore->lastWasRead = 1;

    if (1 != AS_UTL_safeRead(myStore->fp,&length,"getVLRecordStore",sizeof(VLSTRING_SIZE_T), 1))
      assert(0);

    if(!(length + offset + sizeof(VLSTRING_SIZE_T) <= myStore->header.lastElem &&
	 length <= maxLength)){
      fprintf(stderr,"* getVLRecordStore(file) FAILURE: Length = " F_VLS " offset = "F_S64" "
	      "maxLength = " F_VLS " lastElem = "F_S64"\n",
	      length,offset, maxLength, myStore->header.lastElem);
      assert(0);
    }

    if (length != AS_UTL_safeRead(myStore->fp,buffer,"getVLRecordStore",sizeof(char), length))
      assert(0);
  }
  *actualLength = length;
  return(0);
}


int appendStringStore(StoreHandle s, char *string){
  StoreStruct *myStore = Store_myStruct(s);
  int32 length;
  int64 offset;

  assert(myStore->status == ActiveStore);
  assert(Store_isStringStore(myStore));

  offset = sizeof(StoreStat) + myStore->header.lastElem;

  length = strlen(string);
  myStore->header.lastElem += length + 1;
  myStore->isDirty = 1;

  if(myStore->isMemoryStore){
    if(myStore->allocatedSize <= offset + length + 1 ){
      while(myStore->allocatedSize <= offset + length + 1){
	myStore->allocatedSize *= 2;
      }
      myStore->buffer = (char *)
	safe_realloc(myStore->buffer, myStore->allocatedSize);
    }
    memcpy(myStore->buffer + offset, string, length + 1);
  }else{
    AS_UTL_fseek(myStore->fp, (off_t) offset, SEEK_SET);

    if (myStore->lastWasRead)
      fflush(myStore->fp);
    myStore->lastWasRead = 0;

    AS_UTL_safeWrite(myStore->fp, string, "appendStringStore", sizeof(char), length + 1);
  }

  return(0);
}


int appendVLRecordStore(StoreHandle s, void *vlr, VLSTRING_SIZE_T length){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  assert(myStore->status == ActiveStore);
  assert(Store_isVLRecordStore(myStore));
  

  offset = sizeof(StoreStat) + myStore->header.lastElem;

  if(myStore->isMemoryStore){
    if(myStore->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T) ){
      while(myStore->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T)){
	myStore->allocatedSize *= 2;
      }

      myStore->buffer = (char *)
	safe_realloc(myStore->buffer, myStore->allocatedSize);
    }
    memcpy(myStore->buffer + offset, &length, sizeof(VLSTRING_SIZE_T));
    if(length > 0)
      memcpy(myStore->buffer + offset + sizeof(VLSTRING_SIZE_T), vlr, length);
  }else{
    AS_UTL_fseek(myStore->fp, (off_t) offset, SEEK_SET);

    if (myStore->lastWasRead)
      fflush(myStore->fp);
    myStore->lastWasRead = 0;

    AS_UTL_safeWrite(myStore->fp, &length, "appendVLRecordStore", sizeof(VLSTRING_SIZE_T), 1);
    AS_UTL_safeWrite(myStore->fp, vlr, "appendVLRecordStore", sizeof(char), length);
  }

  myStore->header.lastElem += length + sizeof(VLSTRING_SIZE_T);
  myStore->isDirty = 1;

  return(0);
}


int nextStream(StreamHandle sh, void *buffer){
  StreamStruct *ss= Stream_myStruct(sh);

  if(ss->status != ActiveStore){
    fprintf(stderr,"nextStream: status = %d start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->status,ss->startIndex, ss->endIndex);
    assert(0);
  }

  if(ss->startIndex >  ss->endIndex){
    return(0);
  }

  getIndexStore(ss->store, ss->startIndex++, buffer);
  return(1);

}


int kNextStream(StreamHandle sh, void *buffer, int skipNum){
  StreamStruct *ss= Stream_myStruct(sh);

  if(ss->status != ActiveStore){
    fprintf(stderr,"nextStream: status = %d start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->status,ss->startIndex, ss->endIndex);
    assert(0);
  }

  if(ss->startIndex >  ss->endIndex){
    return(0);
  }

  getIndexStore(ss->store, ss->startIndex, buffer);
  ss->startIndex+=skipNum;
  return(1);

}


int nextStringStream(StreamHandle sh, char *buffer, int32 maxLength){

  StreamStruct *ss= Stream_myStruct(sh);

  assert(ss->status == ActiveStore);

  if(ss->startIndex >  ss->endIndex){
    return(0);
  }

  getStringStore(ss->store, ss->startIndex, buffer, maxLength);
  ss->startIndex += strlen(buffer) + 1;
  return(1);

}


int nextVLRecordStream(StreamHandle sh, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength){

  StreamStruct *ss= Stream_myStruct(sh);
  int result;
  assert(ss->status == ActiveStore);

  if(ss->startIndex >  ss->endIndex){
    return(0);
  }

  result = getVLRecordStore(ss->store, ss->startIndex, buffer, maxLength, actualLength);
  assert(result == 0);
  ss->startIndex += (*actualLength + sizeof(VLSTRING_SIZE_T));

  return(1);

}


int statsStore(StoreHandle s, StoreStat *stats){
  StoreStruct *myStore = Store_myStruct(s);
  /*** HACK!!!! This should be done properly  */
  *((StoreStat *)stats) = myStore->header;
  return(0);
}


int64 getLastElemStore(StoreHandle store){
  StoreStruct *myStore = Store_myStruct(store);
  return myStore->header.lastElem;
}


int64 getFirstElemStore(StoreHandle store){
  StoreStruct *myStore = Store_myStruct(store);
  return myStore->header.firstElem;
}


int getStartIndexStream(StreamHandle stream){
  StreamStruct *ss= Stream_myStruct(stream);

  assert(ss->status == ActiveStore);

  return(ss->startIndex);
}
