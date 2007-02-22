
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
static char CM_ID[] = "$Id: AS_PER_genericStore.c,v 1.12 2007-02-22 14:44:40 brianwalenz Exp $";
/*************************************************************************
 Module:  AS_PER_genericStore
 Description:
     This module defines the interface and implementation of the index and
 string "stores" -- binary files used to store fixed and variable size data
 for the assembler.  An additional type of store needs to be implemented, for
 storing variable length records.  Strings should probably be implemented as
 a variable length record prefixed by its length.
     The Fragment Store is built using the building blocks of the index and string
 stores.
     Both types of stores provide both a file-based and memory-based implementations.
 The intent is for the client to build up buffers of records/strings in a memory-based
 store, and concatenate them to the file-based store using the (as yet unimplemented)
 concat operation.
     A store consists of a fixed length header, followed by data.  The only types of
 modifications to the store that are currently permitted are to append a new record/string
 to the store, or to mark a record in an indexStore as deleted.  In principle, there is no
 reason why the index store could not support a replace operation, that would replace a record
 with new data, as long as the index of the record in question was within range.

     Client code relates to a store through an opaque handle.  This will faciliate changes to the
 store structure as we learn more about requirements and optimization.

     Each type of store supports a "stream" operation, for read access to successive elements of the store.

 Assumptions:
      
 Document:
      GenericStore.rtf

 *************************************************************************/

/* RCS Info
 * $Id: AS_PER_genericStore.c,v 1.12 2007-02-22 14:44:40 brianwalenz Exp $
 * $Revision: 1.12 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_fileIO.h"

#define INITIAL_ALLOCATION (2048)
#define WRITING_BUFFER     (16 * 1024)


/* This is the structure maintained in memory for each open Store */
typedef struct{
  FILE *fp;            /* For a file-based store */
  char FileName[FILENAME_MAX];
  char *buffer;        /* For a memory-based store, also holds setbuffer() buffer for disk stores */
  int64 allocatedSize;
  StoreStatus status; 
  StoreStat header;
  int64 lastCommittedElem; /* Initially -1.  If >0, index of last committed element */
  int isMemoryStore;
  int isDirty;
}StoreStruct;


/* Utility routines for reading, writing, and dumping a store Header */
static int readHeader(StoreStruct *myStore);
static int writeHeader(StoreStruct *myStore);
static int dumpHeader(StoreStruct *s);

/* Utility routines for copying stores */
int  copyStore(StoreStruct *source, int64 sourceOffset, int64 sourceMaxOffset, 
	       StoreStruct *target, int64 targetOffset);
int readBufFromStore(StoreStruct *source, char *buffer, int64 offset, 
		     int64 numBytes, char **result);
int writeBufToStore(StoreStruct *target, char *srcPtr, int64 offset, 
		    int64 numBytes);

/* Functions to be inlined or macroified to hide StoreStruct */

static int Store_isMemoryStore(StoreStruct *s){
  return s->isMemoryStore;
}

static int Store_isDirty(StoreStruct *s){
  return s->isDirty;
}

static void Store_setType(StoreStruct *s, int32 type){
  assert(type == INDEX_STORE ||
         type == STRING_STORE ||
         type == VLRECORD_STORE);
  s->header.type = type;
}
static int32 Store_getType(StoreStruct *s){
  return s->header.type;
}
static int32 Store_isIndexStore(StoreStruct *s){
  return s->header.type == INDEX_STORE;
}
static int32 Store_isStringStore(StoreStruct *s){
  return s->header.type == STRING_STORE;
}
static int32 Store_isVLRecordStore(StoreStruct *s){
  //  fprintf(stderr,"* Header type is %d\n", s->header.type);
  return s->header.type == VLRECORD_STORE;
}



/** Globals **/

static StoreStruct *gStores = NULL;
static int gNumStores = 0;
static int gMaxNumStores = 0;

static StreamStruct *gStreams = NULL;
static int gNumStreams = 0;
static int gMaxNumStreams = 0;

static int Store_myHandle(StoreStruct *s){
  return s - gStores; /* Pointer Arithmetic */
}
static int Stream_myHandle(StreamStruct *s){
  return s - gStreams; /* Pointer Arithmetic */
}

static StoreStruct *Store_myStruct(StoreHandle s){
  if(s < 0 || s > gMaxNumStores)
    return NULL;
  return gStores + s; /* Pointer Arithmetic */
}
static StreamStruct * Stream_myStruct(StreamHandle s){
  if(s < 0 || s > gMaxNumStreams)
    return NULL;
  return gStreams + s; /* Pointer Arithmetic */
}




/* Function: computeOffset
 *  Description:  computes offset into store of element with index
 *
 */

static int64 computeOffset(StoreStruct *myStore, int64 index ){
  int64 offset;
  assert(Store_isIndexStore(myStore));

  assert( index >= myStore->header.firstElem);

  offset = sizeof(StoreStat) + 
	((int64)index - (int64)myStore->header.firstElem) * (int64)myStore->header.elementSize;

  assert(offset >= (int64)sizeof(StoreStat));

  return offset;
}


/************************************************************************/
/* Memory allocation for stores and streams                             */
/************************************************************************/


StoreStruct  *allocateStore(void){
  int i;
  StoreStruct *ss = NULL;
  for(i = 0; i < gMaxNumStores;i++){
    if(gStores[i].status == UnAllocatedStore){
      ss = gStores + i;
      break;
    }
  }
  if(ss == NULL){
    if(gNumStores  >= gMaxNumStores){
      if(gNumStores == 0)
	gMaxNumStores = 16;
      else
	gMaxNumStores *= 2;
      gStores = (StoreStruct *) safe_realloc(gStores, gMaxNumStores * sizeof(StoreStruct));
      for  (i = gNumStores;  i < gMaxNumStores;  i ++)
        gStores [i] . status = UnAllocatedStore;
    }
    ss = gStores + gNumStores;
  }

  ss->fp = (FILE *)NULL;
  ss->FileName[0] = 0;
  ss->buffer = NULL;
  ss->allocatedSize = 0;
  ss->status = UnInitializedStore;
  ss->header.firstElem = -1;
  ss->header.lastElem = -1;
  ss->lastCommittedElem = -1;

  ss->header.isDeleted = 0;
  ss->header.type = 0;
  ss->header.p1 = 0;
  ss->header.p2 = 0;
  ss->header.storeType[0] = 0;
  ss->header.firstElem = -1;
  ss->header.lastElem = -1;
  ss->header.version = -1;
  ss->header.elementSize = -1;
  ss->header.creationTime = 0;
  ss->header.lastUpdateTime = 0;

  ss->lastCommittedElem = -1;
  ss->isMemoryStore = 0;
  ss->isDirty = 0;

  memset(ss->FileName, 0, FILENAME_MAX);
  memset(ss->header.storeType, 0, 8);

  gNumStores++;
  return ss;
}





/************************************************************************/
StreamStruct  *allocateStream(StoreHandle s, void *buffer, int32 bufferSize  ){
  int i;
  StreamStruct *ss = NULL;

  for(i = 0; i < gMaxNumStreams;i++){
    if(gStreams[i].status == UnAllocatedStore){
      ss = gStreams + i;
      gNumStreams++;
      break;
    }
  }
  if(ss == NULL){
    if(gNumStreams  >= gMaxNumStreams){
      if(gNumStreams == 0)
	gMaxNumStreams = 16;
      else
	gMaxNumStreams *= 2;
      gStreams = (StreamStruct *)
	safe_realloc(gStreams, gMaxNumStreams * sizeof(StreamStruct));
      for  (i = gNumStreams;  i < gMaxNumStreams;  i ++)
        gStreams [i] . status = UnAllocatedStore;
    }
    ss = gStreams + gNumStreams++;
  }
  ss->store = s;
  ss->status = ActiveStore;
  ss->buffer = buffer;
  ss->bufferSize = bufferSize;
  /* Default is to stream from beginning to end */
  ss->startIndex = gStores[s].header.firstElem;
  ss->endIndex = gStores[s].header.lastElem;

#ifdef DEBUG
  fprintf(stderr," Allocating Stream %d with ("F_S64","F_S64")\n",
	  gNumStreams, ss->startIndex, ss->endIndex);
#endif  
  return ss;
}

/****************************************************************************/

StreamHandle openStream
( StoreHandle fs,/* handle to a fragment store */
  void *buffer,  /* User supplied buffer for prefetching */
  int32 bufferSize)
{
  StreamStruct *stream = allocateStream(fs,buffer,bufferSize);
  return Stream_myHandle(stream);
}

/****************************************************************************/

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
  
/************************************************************************/


StreamHandle openStringStream
( StoreHandle fs, /* handle to a fragment store */
  int64 startOffset, /* Offset to start streaming */
  void *buffer, /* User supplied buffer for prefetching */
  int32 bufferSize)
{

  StreamStruct *stream = allocateStream(fs,buffer,bufferSize);
  resetStream(Stream_myHandle(stream), startOffset, STREAM_UNTILEND);
  return Stream_myHandle(stream);
}


/****************************************************************************/
int closeStream(StreamHandle sh){
  StreamStruct *stream = Stream_myStruct(sh);

  stream->status = UnAllocatedStore;
  gNumStreams--;
  return 0;
}

/****************************************************************************/

StoreHandle createIndexStore
( const char *path, const char *storeType, int32 elementSize, 
  int32 version, int64 firstID)
{
  //  int myStoreIndex;
  StoreStruct *myStore;
  myStore = allocateStore();

  AssertPtr(myStore);

  Store_setType(myStore,INDEX_STORE);
  strncpy(myStore->header.storeType, storeType,7);
  myStore->header.elementSize = elementSize;
  myStore->header.version = version;
  myStore->header.firstElem = firstID;
  myStore->header.lastElem =  firstID -1; /* Sentinel */
#ifdef DEBUG
 fprintf(stderr,
         "createIndexStore with firstID = "F_S64" lastElem = "F_S64"\n",
	 firstID, myStore->header.lastElem);
#endif
  myStore->lastCommittedElem = -1;
  myStore->status = UnInitializedStore;
  myStore->header.creationTime = time(0);
  myStore->status = ActiveStore;

  if(path){ /* File based store */
#ifdef DEBUG
    fprintf(stderr," Creating file-based index store %s\n", path);
#endif
    assert(strlen(path) < FILENAME_MAX);
    myStore->fp = fopen(path,"w+b");
    AssertPtr(myStore->fp);
    myStore->buffer = (char *)safe_calloc(1, WRITING_BUFFER);
    setbuffer(myStore->fp, myStore->buffer, WRITING_BUFFER);
    myStore->isMemoryStore = 0;
    strcpy(myStore->FileName, path);
  }else{
#ifdef DEBUG
    fprintf(stderr," Creating memory-based index store \n");
#endif
    myStore->allocatedSize = INITIAL_ALLOCATION * elementSize;
    myStore->buffer = (char *)safe_calloc(1, myStore->allocatedSize);
    myStore->isMemoryStore = 1;
    myStore->fp = NULL;
    myStore->FileName[0]='\0';
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


/****************************************************************************/

/* createStringStore:
      Create a new String Store and open it for random read and append
*/
StoreHandle createStringStore
( const char *path, const char *storeType, 
  int32 expectedStringSize, int32 version)
{
  // int myStoreIndex;
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
#ifdef DEBUG
    fprintf(stderr," Creating file-based string store %s\n", path);
#endif
    myStore->isMemoryStore = 0;
    assert(strlen(path) < FILENAME_MAX);
    myStore->isMemoryStore = 0;
    strcpy(myStore->FileName, path);
    myStore->fp = fopen(path,"w+b");
    AssertPtr(myStore->fp);
    myStore->buffer = (char *)safe_calloc(1, WRITING_BUFFER);
    setbuffer(myStore->fp, myStore->buffer, WRITING_BUFFER);
  }else{
#ifdef DEBUG
    fprintf(stderr," Creating memory-based string store \n");
#endif
    myStore->allocatedSize = expectedStringSize;
    myStore->buffer = (char *)safe_calloc(1, myStore->allocatedSize);
    AssertPtr(myStore->buffer);
    myStore->isMemoryStore = 1;
    myStore->fp = NULL;
    myStore->FileName[0] = '\0';
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

/****************************************************************************/

/* createVLRecordStore:
      Create a new VLRecord Store and open it for random read and append
*/
StoreHandle createVLRecordStore
( const char *path, const char *storeType, 
  int32 expectedVLRecordSize, int32 version)
{
  StoreHandle newStore 
    = createStringStore(path, storeType, expectedVLRecordSize, VLRECORDSTORE_VERSION);
  StoreStruct *s = Store_myStruct(newStore);

  AssertPtr(s);
  Store_setType(s,VLRECORD_STORE);

  return newStore;

}


StoreHandle resetVLRecordStore( StoreHandle sh){
  return resetStringStore(sh);
}


/****************************************************************************/
int closeStore(StoreHandle s){
  StoreStruct *myStore = Store_myStruct(s);


#ifdef DEBUG
 fprintf(stderr,"*** closeStore %d  status = %d fileno = %d\n", 
	 s, myStore->status, (myStore->fp?fileno(myStore->fp):-1));
#endif
  assert(myStore->status == ActiveStore);

  if(myStore->isMemoryStore){
    myStore->allocatedSize = 0;
    safe_free(myStore->buffer);
  }else{
    if(myStore->isDirty)
    {
      writeHeader(myStore);
      fflush(myStore->fp);
    }
    if(fclose(myStore->fp) != 0)
      assert(0);
    safe_free(myStore->buffer);
  }
  myStore->status = UnAllocatedStore;
  return 0;
}
/****************************************************************************/
int commitStore(StoreHandle s){
  StoreStruct *myStore = Store_myStruct(s);

#ifdef DEBUG_COMMIT
  fprintf(stderr,"*** commitStore %d  dirty:%d\n", s, myStore->isDirty);
#endif

  assert(myStore->status == ActiveStore);
  if(myStore->isDirty){
    myStore->lastCommittedElem = myStore->header.lastElem;
    writeHeader(myStore);
    if(!myStore->isMemoryStore){
      if(fflush(myStore->fp) != 0)
        assert(0);
    }
    myStore->isDirty = 0;
  }
  return 0;
}



/****************************************************************************/
int appendIndexStore(StoreHandle s, void *element){
  StoreStruct *myStore = Store_myStruct(s);
  //  int elem;
  int64 offset;

  assert(myStore->status == ActiveStore);
  myStore->isDirty = 1;

  offset = computeOffset(myStore, myStore->header.lastElem + 1);

  #ifdef DEBUG
  fprintf(stderr,
          " *** IndexStore_append -- Seeking element "F_S64" at "F_S64"\n",
	  myStore->header.lastElem, offset);
  #endif

  if(myStore->isMemoryStore){
    if(myStore->allocatedSize <= offset + 
       (myStore->header.lastElem - myStore->header.firstElem + 1) * myStore->header.elementSize){
      while(myStore->allocatedSize <= offset + 
	    (myStore->header.lastElem - myStore->header.firstElem + 1) * myStore->header.elementSize){
	myStore->allocatedSize *= 2;
      }
      myStore->buffer = (char *)
	safe_realloc(myStore->buffer, myStore->allocatedSize);
    }
    memcpy(myStore->buffer + offset, 
	   element, myStore->header.elementSize);
  }else{
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    AS_UTL_safeWrite(myStore->fp, element, "appendIndexStore", myStore->header.elementSize);
  }
  myStore->header.lastElem ++;

#ifdef DEBUG1
  fprintf(stderr,"*** appendElements "F_S64"\n",
	  myStore->header.lastElem);
#endif

  return(0);
}

/****************************************************************************/
int setIndexStore(StoreHandle s, int64 index, void *element){
  StoreStruct *myStore = Store_myStruct(s);
  //  int elem;
  int64 offset;

  assert(myStore->status == ActiveStore);
  assert(myStore->header.firstElem <= index &&
         myStore->header.lastElem >= index);

  myStore->isDirty = 1;

  offset = computeOffset(myStore, index);

#ifdef DEBUG
  fprintf(stderr,
          " *** IndexStore_set -- Seeking element "F_S64" at "F_S64"\n",
	  index, offset);
#endif

  if(myStore->isMemoryStore){
    memcpy(myStore->buffer + offset, 
	   element, myStore->header.elementSize);
  }else{
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    AS_UTL_safeWrite(myStore->fp, element, "setIndexStore", myStore->header.elementSize);
  }

  return(0);
}

/****************************************************************************/
StoreHandle openStore
( const char *path, /* Path to file */
  const char *rw    /* "r" or "r+" */)
{
  StoreStruct *myStore;

  myStore = allocateStore();
  myStore->fp = fopen(path,rw);
  if(myStore->fp == NULL){
    fprintf(stderr,"* Failure opening Store %s for %s\n", path, rw);
    return(NULLSTOREHANDLE);
  }

  readHeader(myStore);
  myStore->status = ActiveStore;

  return myStore -gStores;
}

/****************************************************************************/


StoreHandle convertStoreToPartialMemoryStore(StoreHandle loadStoreHandle,
                                             int64       firstElem,
                                             int64       lastElem) {
  StoreStruct *loadStore       = Store_myStruct(loadStoreHandle);
  StoreStruct *myStore         = NULL;
  int64        sourceOffset    = 0;
  int64        sourceMaxOffset = 0;

  if (firstElem == STREAM_FROMSTART)
    firstElem = loadStore->header.firstElem;

  if (lastElem == STREAM_UNTILEND)
    lastElem = loadStore->header.lastElem;

  if(Store_isStringStore(loadStore)){
    sourceOffset    = sizeof(StoreStat) + firstElem;
    sourceMaxOffset = sizeof(StoreStat) + lastElem + 1;
  }else{
    sourceOffset    = computeOffset(loadStore, firstElem);
    sourceMaxOffset = computeOffset(loadStore, lastElem + 1);
  }

  myStore  = allocateStore();
  *myStore = *loadStore;

  myStore->allocatedSize    = sourceMaxOffset - sourceOffset + 1;
  myStore->buffer           = (char *)safe_calloc(1, myStore->allocatedSize);
  myStore->isMemoryStore    = 1;
  myStore->fp               = NULL;
  myStore->header.firstElem = firstElem;
  myStore->header.lastElem  = lastElem;

  fprintf(stderr, "* Copying "F_S64" bytes to memory store.\n", myStore->allocatedSize);

  copyStore(loadStore, sourceOffset, sourceMaxOffset,
            myStore, sizeof(StoreStat));

  closeStore(loadStoreHandle);

  return(myStore - gStores);
}

StoreHandle loadStorePartial(const char *StorePath,
                             int64 firstElem,
                             int64 lastElem) {
  return(convertStoreToPartialMemoryStore(openStore(StorePath, "r"), firstElem, lastElem));
}

StoreHandle loadStore(const char *StorePath) {
  return(loadStorePartial(StorePath, STREAM_FROMSTART, STREAM_UNTILEND));
}


/****************************************************************************/
int deleteIndexStore(StoreHandle s, int64 index){
  StoreStruct *myStore = Store_myStruct(s);
  char flags;
  int64 offset = computeOffset(myStore,index);
  myStore->isDirty = 1;
  if(myStore->isMemoryStore){
    *(myStore->buffer + offset) |= 0x80;
  }else{
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    AS_UTL_safeRead(myStore->fp, &flags, "deleteIndexStore", 1);
    flags |= 0x80;
    AS_UTL_safeWrite(myStore->fp, (void *)&flags, "deleteIndexStore", 1);
  }
  return 0;
}


/****************************************************************************/
int getIndexStore(StoreHandle s, int64 index, void *buffer){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  assert(myStore->status == ActiveStore);
#ifdef DEBUG
  fprintf(stderr," *** IndexStoreGet -- Seeking element "F_S64" ("F_S64","F_S64")\n",
	  index, myStore->header.firstElem, myStore->header.lastElem);
  fprintf(stderr," *** IndexStoreGet -- fileno = %d\n",
	  fileno(myStore->fp));
#endif
  offset = computeOffset(myStore,index);

  if(myStore->isMemoryStore){
    memcpy(buffer, myStore->buffer + offset, myStore->header.elementSize);
    return(0);
  }else{
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    
    return AS_UTL_safeRead(myStore->fp,buffer,"getIndexStore",myStore->header.elementSize);
  }
}



/****************************************************************************/
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
    if (CDS_FTELL(myStore->fp) != (off_t)(offset + sizeof(StoreStat)))
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) (offset + sizeof(StoreStat)), SEEK_SET))
        assert(0);
    
    if(AS_UTL_safeRead(myStore->fp,buffer,"getStringStore",length) != FALSE)
      assert(0);
    
  }
  return(0);


  }
/****************************************************************************/
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
    if (CDS_FTELL(myStore->fp) != (off_t)(actualOffset + sizeof(StoreStat)))
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) (actualOffset + sizeof(StoreStat)), SEEK_SET))
        assert(0);
    if(AS_UTL_safeRead(myStore->fp,&length,"getVLRecordStore",sizeof(VLSTRING_SIZE_T)) != FALSE)
      assert(0);

    if(!(length + offset + sizeof(VLSTRING_SIZE_T) <= myStore->header.lastElem &&
	 length <= maxLength)){
      fprintf(stderr,"* getVLRecordStore(file) FAILURE: Length = " F_VLS " offset = "F_S64" "
	      "maxLength = " F_VLS " lastElem = "F_S64"\n",
	      length,offset, maxLength, myStore->header.lastElem);
      assert(0);
    }
    if(length > 0)
      if(AS_UTL_safeRead(myStore->fp,buffer,"getVLRecordStore",length) != FALSE)
        assert(0);
  }
  *actualLength = length;
  return(0);


  }

/****************************************************************************/
int appendStringStore(StoreHandle s, char *string){
  StoreStruct *myStore = Store_myStruct(s);
  int32 length;
  int64 offset;

  assert(myStore->status == ActiveStore);
  assert(Store_isStringStore(myStore));

  offset = sizeof(StoreStat) + myStore->header.lastElem;
#ifdef DEBUG
  fprintf(stderr," *** StringStore_append -- Seeking offset "F_S64"\n",
	  offset);
#endif
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
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    AS_UTL_safeWrite(myStore->fp, string, "appendStringStore", length + 1);
#ifdef DEBUG
    fprintf(stderr,"appendString: wrote %s of length " F_SIZE_T "\n",
	    string, strlen(string));
#endif
  }

  return(0);
}

/****************************************************************************/
int appendVLRecordStore(StoreHandle s, void *vlr, VLSTRING_SIZE_T length){
  StoreStruct *myStore = Store_myStruct(s);
  int64 offset;

  assert(myStore->status == ActiveStore);
  assert(Store_isVLRecordStore(myStore));
  

  offset = sizeof(StoreStat) + myStore->header.lastElem;
#ifdef DEBUG
  fprintf(stderr," *** appendVLRecordStore -- length = " F_VLS " lastElem = "F_S64" offset = "F_S64"\n",
	  length, myStore->header.lastElem, offset);
#endif
  if(myStore->isMemoryStore){
    if(myStore->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T) ){
      while(myStore->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T)){
	myStore->allocatedSize *= 2;
      }
#ifdef DEBUG
      fprintf(stderr,"* Reallocating VLRecord Store %d to "F_S64"\n",
              s, myStore->allocatedSize);
#endif
      myStore->buffer = (char *)
	safe_realloc(myStore->buffer, myStore->allocatedSize);
    }
    memcpy(myStore->buffer + offset, &length, sizeof(VLSTRING_SIZE_T));
    if(length > 0)
      memcpy(myStore->buffer + offset + sizeof(VLSTRING_SIZE_T), vlr, length);
  }else{
    if (CDS_FTELL(myStore->fp) != offset)
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) offset, SEEK_SET))
        assert(0);
    AS_UTL_safeWrite(myStore->fp, &length, "appendVLRecordStore", sizeof(VLSTRING_SIZE_T));
    if (CDS_FTELL(myStore->fp) != (off_t)(offset + sizeof(VLSTRING_SIZE_T)))
      if(-1 == CDS_FSEEK(myStore->fp, (off_t) (offset + sizeof(VLSTRING_SIZE_T)), SEEK_SET))
        assert(0);
    if(length > 0)
      AS_UTL_safeWrite(myStore->fp, vlr, "appendVLRecordStore", length);
  }

  myStore->header.lastElem += length + sizeof(VLSTRING_SIZE_T);
  myStore->isDirty = 1;
#ifdef DEBUG
    fprintf(stderr,
            "appendVLRecord: wrote  length " F_VLS " lastElem = "F_S64"\n",
	    length, myStore->header.lastElem);
#endif

  return(0);
}


/****************************************************************************/
/* nextStream
   Streaming Operation.
*/
	
int nextStream(StreamHandle sh, void *buffer){
  StreamStruct *ss= Stream_myStruct(sh);
  //  int next;

  if(ss->status != ActiveStore){
    fprintf(stderr,"nextStream: status = %d start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->status,ss->startIndex, ss->endIndex);
    assert(0);
  }

  if(ss->startIndex >  ss->endIndex){
    #ifdef DEBUG
    fprintf(stderr,"nextStream: start = "F_S64" > end = "F_S64".... bailing out\n", 
	    ss->startIndex, ss->endIndex);
    #endif
    return(0);
  }

  getIndexStore(ss->store, ss->startIndex++, buffer);
  return(1);

}


int kNextStream(StreamHandle sh, void *buffer, int skipNum){
  StreamStruct *ss= Stream_myStruct(sh);
  //  int next;

  if(ss->status != ActiveStore){
    fprintf(stderr,"nextStream: status = %d start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->status,ss->startIndex, ss->endIndex);
    assert(0);
  }

  if(ss->startIndex >  ss->endIndex){
    #ifdef DEBUG
    fprintf(stderr,"nextStream: start = "F_S64" > end = "F_S64".... bailing out\n", 
	    ss->startIndex, ss->endIndex);
    #endif
    return(0);
  }

  getIndexStore(ss->store, ss->startIndex, buffer);
  ss->startIndex+=skipNum;
  return(1);

}



/****************************************************************************/
int nextStringStream(StreamHandle sh, char *buffer, int32 maxLength){

  StreamStruct *ss= Stream_myStruct(sh);

  assert(ss->status == ActiveStore);

  if(ss->startIndex >  ss->endIndex){
#ifdef DEBUG
    fprintf(stderr,"StringStoreNext: start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->startIndex, ss->endIndex);
#endif
    return(0);
  }

#ifdef DEBUG
  fprintf(stderr,"StringStoreNext: offset = "F_S64"\n", ss->startIndex);
#endif
  getStringStore(ss->store, ss->startIndex, buffer, maxLength);
  ss->startIndex += strlen(buffer) + 1;
  return(1);

}

/****************************************************************************/
int nextVLRecordStream(StreamHandle sh, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength){

  StreamStruct *ss= Stream_myStruct(sh);
  int result;
  assert(ss->status == ActiveStore);

  if(ss->startIndex >  ss->endIndex){
#ifdef DEBUG
    fprintf(stderr,"StringVLRecordNext: start = "F_S64" >= end = "F_S64" .... bailing out\n", 
	    ss->startIndex, ss->endIndex);
#endif
    return(0);
  }

#ifdef DEBUG
  fprintf(stderr,"nextVLRecordStream: offset = "F_S64"\n", ss->startIndex);
#endif
  result = getVLRecordStore(ss->store, ss->startIndex, buffer, maxLength, actualLength);
  assert(result == 0);
  ss->startIndex += (*actualLength + sizeof(VLSTRING_SIZE_T));

#ifdef DEBUG
  fprintf(stderr,
          "nextVLRecordStream: read:  " F_VLS " startIndex now:"F_S64"\n",
	  *actualLength, ss->startIndex);

#endif

  return(1);

}


/****************************************************************************/
int statsStore(StoreHandle s, StoreStat *stats){
  StoreStruct *myStore = Store_myStruct(s);

#ifdef DEBUG
  dumpHeader(myStore);
#endif
  /*** HACK!!!! This should be done properly  */
  *((StoreStat *)stats) = myStore->header;

#ifdef DEBUG
  fprintf(stderr,"header.firstElem = "F_S64" header.lastElem = "F_S64"\n",
	  myStore->header.firstElem, myStore->header.lastElem);
  fprintf(stderr,"stats.firstElem = "F_S64" stats.lastElem = "F_S64"\n",
	  stats->firstElem, stats->lastElem);
#endif
  return(0);

}
/***********************************************************************************/

int64 getLastElemStore(StoreHandle store){
  StoreStruct *myStore = Store_myStruct(store);
  return myStore->header.lastElem;
}

 /***********************************************************************************/

int64 getFirstElemStore(StoreHandle store){
  StoreStruct *myStore = Store_myStruct(store);
  return myStore->header.firstElem;
}

/***********************************************************************************/
int getStartIndexStream(StreamHandle stream){
  StreamStruct *ss= Stream_myStruct(stream);

  assert(ss->status == ActiveStore);

  return(ss->startIndex);

}



/****************************************************************************/

int writeBufToStore(StoreStruct *target, char *srcPtr, int64 offset, int64 numBytes){
#ifdef DEBUG_COPY
  fprintf(stderr,
          "*writeBufToStore (%s)write "F_S64" bytes at offset "F_S64"\n",
	  (Store_isMemoryStore(target)?"Memory":"File"), numBytes, offset);
#endif
   assert(target->status == ActiveStore &&
          srcPtr &&
          offset > 0 &&
          numBytes > 0);

   target->isDirty = 1;
   if(Store_isMemoryStore(target)){
      memcpy(target->buffer + offset, (void *)srcPtr, numBytes);
   }else{
     if (CDS_FTELL(target->fp) != offset)
       if(-1 == CDS_FSEEK(target->fp, (off_t) offset, SEEK_SET))
         assert(0);
      if(CDS_FTELL(target->fp) != offset)
        assert(0);
      AS_UTL_safeWrite(target->fp, srcPtr, "writeBufToStore", numBytes);
   }
   return(0);
}


/****************************************************************************/
int readBufFromStore(StoreStruct *source, char *buffer, int64 offset, int64 numBytes, char **result){
#ifdef DEBUG_COPY
  fprintf(stderr,
          "*readBufFromStore (%s) read "F_S64" bytes at offset "F_S64"\n",
	  (Store_isMemoryStore(source)?"Memory":"File"), numBytes, offset);
#endif


   if(Store_isMemoryStore(source)){
      *result = source->buffer + offset;
   }else{

     if (CDS_FTELL(source->fp) != offset)
       CDS_FSEEK(source->fp, (off_t) offset, SEEK_SET); /* set to offset */
     if(CDS_FTELL(source->fp) != offset)
       assert(0);
     AS_UTL_safeRead(source->fp, (void *)buffer, "readBufFromStore", numBytes);
     *result = buffer;
   }
   return(0);
}
   

/****************************************************************************/
#define COPY_SIZE BUFSIZ           /* BUFSIZE is in stdio.h */
  /* This does a low level copy of data */
int  copyStore(StoreStruct *source, int64 sourceOffset, int64 sourceMaxOffset, 
	       StoreStruct *target, int64 targetOffset){

   char buffer[COPY_SIZE];
   char *srcPtr;
   int64 offset = 0;
   int64 maxOffset = sourceMaxOffset - sourceOffset;
   int64 numReads = maxOffset/COPY_SIZE;
   int64 lastRead = maxOffset%COPY_SIZE;
   int64 read;
   int64 newSize;

   assert(source->status == ActiveStore &&
	  target->status == ActiveStore);

   // Unlike disk stores, memory stores cannot automatically resize...so we make sure there
   // is enough space by reallocing.
   if(target->isMemoryStore){
     newSize = targetOffset + (sourceMaxOffset - sourceOffset);
     if(target->allocatedSize <= newSize){
       while(target->allocatedSize <= newSize){
	 target->allocatedSize *= 2;
       }
       target->buffer = (char *) 
	 safe_realloc(target->buffer, target->allocatedSize);
     }
   }

#ifdef DEBUG_COPY
     fprintf(stderr,"CopyStore From %ld to %ld:  sourceOffset:"F_S64"  sourceMaxOffset:"F_S64" numReads:"F_S64" lastRead "F_S64"\n",
	     source - gStores, target - gStores,
             sourceOffset, sourceMaxOffset,
             numReads, lastRead);
#endif

     if(numReads){
#ifdef DEBUG_COPY
     fprintf(stderr," copyStore:  numReads = "F_S64"\n", numReads);
#endif
       for(read = 1; read <= numReads; read++){
	 readBufFromStore(source, buffer, sourceOffset + offset, COPY_SIZE, &srcPtr);
	 writeBufToStore(target, srcPtr, targetOffset + offset, COPY_SIZE);
	 offset += COPY_SIZE;
       }
     }
   if(lastRead){
#ifdef DEBUG_COPY
     fprintf(stderr," copyStore:  lastRead = "F_S64"\n", lastRead);
#endif
      readBufFromStore(source, buffer, sourceOffset + offset, lastRead, &srcPtr);
      writeBufToStore(target, srcPtr, targetOffset + offset, lastRead);
   }

   return(0);
}




/****************************************************************************/
static int writeHeader(StoreStruct *myStore){

  myStore->header.lastUpdateTime = time(0);
#ifdef DEBUG
  fprintf(stderr,"*** Begin WriteHeader:\n");
  dumpHeader(myStore);
#endif
  if(myStore->isMemoryStore){
    memcpy(myStore->buffer,&myStore->header, sizeof(StoreStat));
  }else{
#ifdef DEBUG
  fprintf(stderr,"\t*** WriteHeader: fileno %d status = %d err = %d\n",
	  fileno(myStore->fp), myStore->status, ferror(myStore->fp));
#endif    
    rewind(myStore->fp); /* Rewind to beginning of file */
    AS_UTL_safeWrite(myStore->fp,(void *)&myStore->header,"writeHeader",sizeof(StoreStat));
  }
#ifdef DEBUG
  fprintf(stderr,"*** End WriteHeader:\n");
#endif
  return(0);
}

/****************************************************************************/
static int readHeader(StoreStruct *myStore){

  if(myStore->isMemoryStore){
    memcpy(&myStore->header,&myStore->buffer, sizeof(StoreStat));
  }else{
    CDS_FSEEK(myStore->fp, (off_t) 0, SEEK_SET); /* Rewind to beginning of file */
    AS_UTL_safeRead(myStore->fp, (void *)&myStore->header, "readHeader", sizeof(StoreStat));
  }
#ifdef DEBUG
  fprintf(stderr,"*** ReadHeader: last = "F_S64"\n"
          , myStore->header.lastElem);
  dumpHeader(myStore);
  fprintf(stderr,"*** End ReadHeader:\n");
#endif
  return(0);
}


/****************************************************************************/

int dumpHeader(StoreStruct *s){
  char *t;
  fprintf(stderr,
          "*** Dumping StoreHandle %d's header (status = %d dirty = %d)***\n",
	  Store_myHandle(s), s->status, s->isDirty);
  fprintf(stderr,"\t*** isDeleted %d\n", s->header.isDeleted);
  fprintf(stderr,"\t*** isIndexStore %d\n", Store_isIndexStore(s));
  fprintf(stderr,"\t*** storeType %s\n", s->header.storeType);
  fprintf(stderr,"\t*** firstElem "F_S64"   lastElem "F_S64"\n", 
	  s->header.firstElem, s->header.lastElem);
  fprintf(stderr,"\t*** version %d elementSize %d\n", s->header.version,
	  s->header.elementSize);
  
  t = ctime(&(s->header.creationTime)); 
  fprintf(stderr,"\t*** created %s\n", t);
  t = ctime(&(s->header.lastUpdateTime));
  fprintf(stderr,"\t*** updated %s\n", t);
  fprintf(stderr,"*** ********************** ***\n");

  return 0;
}




#ifdef DEBUG_TEST

char *testString = "abcdefghijklmnopqrstuvwxyz"
                   "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                   "0123456789";


main(int argc, char **argv){
  StoreHandle  mine;
  char buffer[64], tmp[64];
  StreamHandle stream;

  int elements[100];
  int i,x;
  
  for(i = 0; i < 100; i ++)
    memcpy(elements + i, &i, sizeof(int));

  mine = createIndexStore("junk.db", "test", sizeof(int), 3, 22);
  sleep(1);
  IndexStore_appendElements(mine, elements, 100);

  closeStore(mine);

  sleep(1);
  mine = openStore("junk.db","r+");
  fprintf(stderr," \n\n***Opened store junk.db %d\n", mine);

#if 0
  for(i = 121; i >= 22; i--){
    IndexStoreGet(mine, i, &x);
    /*    fprintf(stderr,"Element %d is %d\n", i, x);*/
  }
#endif
  fprintf(stderr,"\n\n********* Streaming Store *************\n");
  stream = streamStore(mine, 22, NULL, 0);
  fprintf(stderr,"# Opened stream %d\n", stream);
  i = 0;
  fprintf(stderr,"\n\n********* Before While *************\n");
  while(StoreNext(stream, &x)){
    fprintf(stderr,"Element %d is %d\n", i++, x);
  }
  
  
  fprintf(stderr," ***Closing store %d\n", mine);
  closeStore(mine);
  fprintf(stderr,"\n\n");

  mine = createStringStore("sjunk.db", "stest", 64, 3);
  fprintf(stderr," Created String store %d\n", mine);
  sleep(1);
  for(i = 0; i < 100; i++){
    strncpy(buffer,testString + i%32, 30);
    /* fprintf(stderr,"\t\t%d %s\n", i, buffer);*/
    StringStore_appendString(mine, buffer);
  }
  fprintf(stderr,"*** Closing store %d\n\n", mine);
  closeStore(mine);

  sleep(1);
  mine = openStore("sjunk.db","r+");
  fprintf(stderr,"*** Opened String store %d\n", mine);
#if 0
  for(i = 0; i < 100; i++){
    char buf[2048];
    strncpy(buffer,testString + i%32, 30);
    StringStoreGet(mine, i * 31, buf, 2048);
    if(strcmp(buf,buffer)){
      fprintf(stderr,"Comparison Failure on string %d:  wrote %s read %s\n",
	      i, buffer, buf);
    }
  }
#else
  i = 0;
  stream = openStringStream(mine,0,NULL, 0);
  fprintf(stderr,"# Opened stream %d\n");
  while(nextStringStore(stream, tmp, 63)){
    strncpy(buffer,testString + i%32, 30);
    if(strcmp(tmp,buffer)){
      fprintf(stderr,"Comparison Failure on string %d:  wrote %s read %s\n",
	      i, buffer, tmp);
    }
    i++;
  }    
#endif
  closeStore(mine);

}
#endif
