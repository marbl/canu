
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
#define VLSTRING_MAX_SIZE (4 * 1024 * 1024  -1)
#define VLRECORDSTORE_VERSION 2
#else
#define VLSTRING_SIZE_T uint16
#define F_VLS F_U16
#define VLSTRING_MAX_SIZE (64 * 1024 - 1)
#define VLRECORDSTORE_VERSION 1
#endif

#include <time.h>


typedef int StoreHandle;   /* The handle returned by open/create operations */
typedef int StreamHandle ; /* The iterator used by the streaming operation */
 
#define NULLSTOREHANDLE (-1)

typedef enum { UnAllocatedStore = 0, 
	       UnInitializedStore, 
	       ActiveStore} StoreStatus;


#define INDEX_STORE 1
#define STRING_STORE 2
#define VLRECORD_STORE 2
#define INVALID_STORE 0

/*  Structure returned by statStore */

typedef struct{
  unsigned int  isDeleted:1;
  unsigned int  type:3;
  unsigned int  p1:28;  // padding field
  unsigned int  p2:32;  // padding field
  char          storeType[8];
  int64         firstElem; /* Initially -1.  If >0, index of first allocated element */
  int64         lastElem;  /* Initially -1.  If >0, index of last allocated element */
  int32         version;        /* For user information only */
  int32         elementSize;  
  int64         creationTime;
  int64         lastUpdateTime;
}StoreStat;


/* This is the structure maintained for each store Stream */
typedef struct{
  StoreHandle   store;
  void         *buffer;
  int32         bufferSize;
  int64         startIndex;
  int64         endIndex;
  StoreStatus   status; 
}StreamStruct;

/******************************************************************************
 * Function: createIndexStore:
 * Description:
 *     Allocates an index store, and returns its handle.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     storeType   Currently for user use only
 *     elementSize Size of the fixed length records
 *     version     Currently for user use, 
 *                 we may need this to do version conversions in the future.
 *     firstID     Index of the first element in the store.
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle createIndexStore
( const char *StorePath, const char *storeType, 
  int32 elementSize, int32 version, int64 firstID);


/******************************************************************************
 * Function: resetIndexStore:
 * Description:
 *     Recycles a string store, nuking its data.
 *
 * Inputs:
 *     sh           handle of open Index Store
 *     firstID      First ID for reset Store
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle resetIndexStore(StoreHandle sh, int64 firstID);


/******************************************************************************
 * Function: createStringStore:
 * Description:
 *     Allocates a string store, and returns its handle.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     storeType   Currently for user use only
 *     expectedStringSize   May be used in the future for optimization
 *     version     Currently for user use, 
 *                 we may need this to do version conversions in the future.
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle createStringStore
( const char *StorePath, const char *storeType, 
  int32 expectedStringSize, int32 version);




/******************************************************************************
 * Function: resetStringStore:
 * Description:
 *     Recycles a string store, nuking its data.
 *
 * Inputs:
 *     sh           handle of open string Store
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle resetStringStore(StoreHandle sh);


/******************************************************************************
 * Function: createVLRecordStore:
 * Description:
 *     Allocates a VLRecord store, and returns its handle.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     storeType   Currently for user use only
 *     expectedRecordSize   May be used in the future for optimization
 *     version     Currently for user use, 
 *                 we may need this to do version conversions in the future.
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle createVLRecordStore
( const char *StorePath, const char *storeType, 
  int32 expectedRecordSize, int32 version);




/******************************************************************************
 * Function: resetVLRecordStore:
 * Description:
 *     Recycles a VLRecord store, nuking its data.
 *
 * Inputs:
 *     sh           handle of open VLRecord  Store
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

StoreHandle resetVLRecordStore(StoreHandle sh);


/******************************************************************************
 * Function: openStore:
 * Description:
 *     Opens an existing, file-based store.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     rw          file access mode
 *
 * Return Value:
 *     StoreHandle != NULLSTOREHANDLE on success
 *****************************************************************************/

StoreHandle openStore
( const char *StorePath, /* Path to file */
  const char *rw         /* "r" or "rw" */
);



/******************************************************************************
 * Function: closeStore:
 * Description:
 *     Close an open store.  Updates lastElem and lastCommittedElem.
 *
 * Inputs:
 *     sh    Handle to store we want to close
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int closeStore(StoreHandle sh);


/******************************************************************************
 * Function: commitStore:
 * Description:
 *     Like closeStore, but store remains open.
 *
 * Inputs:
 *     sh    Handle to store we want to commit
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int commitStore(StoreHandle sh);



/******************************************************************************
 * Function: getIndexStore
 * Description:
 *     Random access to records in an index store
 *
 * Inputs:
 *     fs         Handle of index store
 *     index      index of record
 * Outputs:
 *     buffer     Buffer for element 
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
	
int getIndexStore(StoreHandle fs, int64 indx, void *buffer);

	
/******************************************************************************
 * Function: getStringStore
 * Description:
 *     Random access to strings in a string store
 *
 * Inputs:
 *     fs         Handle of string store
 *     offset     offset of string
 * Outputs:
 *     buffer     Buffer for element 
 *     maxLength  size of buffer
 *
 * Return Value:
 *     Zero if success.
 *     Hitting an EOF before an EOS is a no-no.
 *****************************************************************************/
	
int getStringStore(StoreHandle s, int64 offset, char *buffer, int32 maxLength);

/******************************************************************************
 * Function: getVLRecordStore
 * Description:
 *     Random access to VLRecords in a VLRecord store
 *
 * Inputs:
 *     fs         Handle of VLRecord store
 *     offset     offset of VLRecord
 * Outputs:
 *     buffer     Buffer for element 
 *     maxLength  size of buffer (if maxLength is less than record size --> error )
 *
 * Return Value:
 *     Zero if success.
 *     If maxLength < recordLength, return -1
 *****************************************************************************/
	
int getVLRecordStore(StoreHandle s, int64 offset, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength);


/******************************************************************************
 * Function: appendVLRecordStore
 * Description:
 *     Append an element to an index store
 *
 * Inputs:
 *     store      Handle of index store
 *     element    Pointer to record for append
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int appendVLRecordStore(StoreHandle store, void *element, VLSTRING_SIZE_T length);


/******************************************************************************
 * Function: setIndexStore
 * Description:
 *     Overwrite an existing  element of an index store
 *
 * Inputs:
 *     store      Handle of index store
 *     index      index of element to overwrite
 *     element    Pointer to record for overwrite
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int setIndexStore(StoreHandle store, int64 indx, void *element);


/******************************************************************************
 * Function: appendStringStore
 * Description:
 *     Append a string to a string store
 *
 * Inputs:
 *     store      Handle of index store
 *     string     String to be appended
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
int appendStringStore(StoreHandle s, char *string);



#if 0
/******************************************************************************
 * Function: deleteIndexStore
 * Description:
 *     Delete an element from an index store
 *
 * Inputs:
 *     store      Handle of index store
 *     index      Index of element to be deleted
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int deleteIndexStore(StoreHandle store, int64 indx);
#endif

/******************************************************************************
 * Function: statsStore
 * Description:
 *     Get statistics on a store
 *
 * Inputs:
 *     store      Handle of index store
 * Outputs:
 *     stats      Pointer to StoreStat structure 
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int statsStore(StoreHandle store, StoreStat *stats);

/******************************************************************************
 * Function: getLastElemStore
 * Description:
 *   Get the index of the last Element of the Store
 *
 * Inputs:
 *     store           FragStore Handle
 * Return Value
 *     lastElem           Pointer to int to hold result.
 *
 *****************************************************************************/

int64 getLastElemStore(StoreHandle store);


/******************************************************************************
 * Function: getFirstElemStore
 * Description:
 *   Get the index of the first Element of the Store
 *
 * Inputs:
 *     store           Store Handle
 * ReturnValue
 *     lastElem           Pointer to int to hold result.
 *
 *****************************************************************************/

int64 getFirstElemStore(StoreHandle store);


/******************************************************************************
 * Function: getStartIndexStream
 * Description:
 *   Get the startindex of the Stream
 *
 * Inputs:
 *     stream           Stream Handle
 * ReturnValue
 *     StartIndex           
 *
 *****************************************************************************/
int getStartIndexStream(StreamHandle stream);


//  "Convert" the loadStore into a new memory store.  The loadStore is
//  closed.
//
StoreHandle convertStoreToPartialMemoryStore(StoreHandle loadStore,
                                             int64       firstElem,
                                             int64       lastElem);

//  Open an existing file-based Store, and load its entire contents
//  into A newly created memory-based Store.  This is equivalent to
//  opening the file-based store, and concatting it to a newly created
//  memory based store.
//
StoreHandle loadStore(const char *StorePath);

//  Open an existing file-based Store, and load a portion of its
//  contents into A newly created memory-based Store.
///
StoreHandle loadStorePartial(const char *StorePath,
                             int64 firstElem,
                             int64 lastElem);


	
/******************************************************************************
 * Function: getStore
 * Description:
 *     Random access to records in an index store
 *
 * Inputs:
 *     fs         Handle of index store
 *     index      index of record
 * Outputs:
 *     buffer     Buffer for element 
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
	
int getIndexStore(StoreHandle fs, int64 indx, void *buffer);

	
/******************************************************************************
 * Function: getStringStore
 * Description:
 *     Random access to strings in a string store
 *
 * Inputs:
 *     fs         Handle of string store
 *     offset     offset of string
 * Outputs:
 *     buffer     Buffer for element 
 *     maxLength  size of buffer
 *
 * Return Value:
 *     Zero if success.
 *     Hitting an EOF before an EOS is a no-no.
 *****************************************************************************/
	
int getStringStore(StoreHandle fs, int64 offset, char *buffer, int32 maxLength);


/******************************************************************************
 * Function: appendIndexStore
 * Description:
 *     Append an element to an index store
 *
 * Inputs:
 *     store      Handle of index store
 *     element    Pointer to record for append
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int appendIndexStore(StoreHandle store, void *element);





/********************************* STREAMS ********************************/

#define STREAM_UNTILEND -1 
#define STREAM_FROMSTART -1

/******************************************************************************
 * Function: openStream:
 * Description:
 *     Open a stream on an open Index Store
 *
 * Inputs:
 *     sh         Handle of store we want to stream
 *     buffer     Memory buffer for prefetching or optimization (currently unused) 
 *     bufferSize Size of Memory buffer (currently unused)
 *
 * Return Value:
 *     Handle of opened Stream
 *****************************************************************************/
StreamHandle openStream(StoreHandle sh, /* handle to a fragment store */
			 void *buffer,  /* User supplied buffer for prefetching */
			 int32 bufferSize);

/******************************************************************************
 * Function: openStringStream:
 * Description:
 *     Opens a stream on an open String Store
 *
 * Inputs:
 *     fs         Handle of store we want to stream
 *     startOffset Offest of first string
 *     buffer     Memory buffer for prefetching or optimization (currently unused) 
 *     bufferSize Size of Memory buffer (currently unused)
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
StreamHandle openStringStream(StoreHandle fs,        /* handle to a fragment store */
			       int64 startOffset,     /* Offset to start streaming */
			       void *buffer,       /* User supplied buffer for prefetching */
				 int32 bufferSize);

#if 0
/******************************************************************************
 * Function: openVLRecordStream:
 * Description:
 *     Opens a stream on an open VLRecord Store
 *
 * Inputs:
 *     sh         Handle of store we want to stream
 *     startOffset Offest of first VLRecord
 *     buffer     Memory buffer for prefetching or optimization (currently unused) 
 *     bufferSize Size of Memory buffer (currently unused)
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
StreamHandle openVLRecordStream
(StoreHandle fs,    /* handle to a VlRecord store */
 int64 startoffset, /* Offset to start streaming */
 void *buffer,      /* User supplied buffer for prefetching */
 int32 bufferSize);
#endif
	


/******************************************************************************
 * Function: resetStream:
 * Description:
 *     Reset a stream with  new start/end indices
 *
 * Inputs:
 *     sh         Handle of store we want to stream
 *     startIndex First element of the store we want to see
 *     endIndex   Last element of the store we want to see
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
int resetStream(StreamHandle sh, /* handle to a fragment store */
			 int64 startIndex, /*  First record to stream */
			 int64 endIndex); /*  Last record to stream, -1 if stream to EOF*/


/******************************************************************************
 * Function: getStartIndexStream
 * Description:
 *   Get the startindex of the frag Stream
 *
 * Inputs:
 *     stream           FragStream Handle
 * ReturnValue
 *     StartIndex           
 *
 *****************************************************************************/
int getStartIndexStream(StreamHandle stream);



/******************************************************************************
 * Function: nextStream
 * Description:
 *     Read the next record from an open index stream
 *
 * Inputs:
 *     sh         Handle of stream
 * Outputs:
 *     buffer     Buffer big enough to hold an element.
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
	
int nextStream(StreamHandle sh, void *buffer);

int kNextStream(StreamHandle sh, void *buffer, int skipNum);

/******************************************************************************
 * Function: nextStringStream
 * Description:
 *     Read the next string from an open string stream
 *
 * Inputs:
 *     sh         Handle of string stream
 * Outputs:
 *     buffer     Buffer for string
 *     maxLength  Size of string that buffer can accomodate
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
int nextStringStream(StreamHandle sh, char *buffer, int32 maxLength);



/******************************************************************************
 * Function: nextVLRecordStream
 * Description:
 *     Read the next VLRecord from an open VLRecord stream
 *
 * Inputs:
 *     sh         Handle of VLRecord Stream
 * Outputs:
 *     buffer     Buffer for VLRecord
 *     maxLength  Size of VLRecord that buffer can accomodate
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/
int nextVLRecordStream(StreamHandle sh, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength);

	

/******************************************************************************
 * Function: closeStream
 * Description:
 *     Closes a stream
 *
 * Inputs:
 *     sh         Handle of stream we want to close
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int closeStream(StreamHandle sh);


#endif /* AS_PER_GENERICSTORE_H */
