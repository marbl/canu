
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
#ifndef AS_PER_FRAGSTORE_H
#define AS_PER_FRAGSTORE_H
/*************************************************************************
 Module:  AS_PER_fragStore

 Description:

     This module defines the interface and implementation of the
     Assembler Fragment Store.  The Fragment Store is built using the
     building blocks of the index and string stores in the
     AS_PER_genericStore module.  Currently, the fragment store
     manages a directory of files containing 3 files:

        - db.frg      Index Store with Fragment data
        - db.seq      String Store with Sequence/Quality data
        - db.src      String Store with source data.

 Assumptions:

      Implementation of Index/String stores.

 Document:

      FragmentStore.rtf

 *************************************************************************/

/* RCS Info
 * $Id: AS_PER_fragStore.h,v 1.3 2005-03-22 19:07:31 jason_miller Exp $
 * $Revision: 1.3 $
 *
 */

// History:
//         Version 1    Base
//         Version 2    added branch-point info
//         Version 3    nuked branch-point info, and added Screen Match info
//         Version 4    added indexed store features
//         Version 5    added multiple clear ranges (Oct 2001 - Jason)

#include "AS_PER_genericStore.h"
#include "AS_PER_ReadStruct.h"

// Compiler -D option sets fragStore version.
// The code supports only one version per compilation.
// Stores older than version 4 are not supported at all.
// All stores were version 4 until Oct 2001.
#ifndef  FRAGSTORE_VERSION
#define  FRAGSTORE_VERSION 5
#endif     // Default fragStore version is 5 (as of Oct 2001).
#if FRAGSTORE_VERSION<=4
#define FRAGSTORE_VERSION 4
#endif     // This is the lowest version supported by this code.
#if FRAGSTORE_VERSION>=5
#define FRAGSTORE_VERSION 5
#endif     // This is the highest version supported by this code.
#define  VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES 5


// In 1999, there was an effort to support re-sorting the frag store.
// Code to support that effort was embedded in tests like
// #if FRAGSTORE_VERSION_PRODUCTION < FRAGSTORE_VERSION.
// All that code was never used. It should be removed completely.
// For now, I will set the values such that the code never compiles.
// - Jason, Oct 2001.
#define  FRAGSTORE_VERSION_PRODUCTION   FRAGSTORE_VERSION


typedef int FragStoreHandle;    /* The handle returned by open/create operations */
typedef int FragStreamHandle;   /* The handle returned by stream operations */

#define NULLFRAGSTOREHANDLE -1


/******************************************************************************
 * Function: statsFragStore
 * Description:
 *     get stats on a fragstore fragments
 *
 * Inputs:
 *     handle    FragStoreHandle
 *
 * Return Value:
 *     None
 *****************************************************************************/
int statsFragStore(FragStoreHandle handle, StoreStat *stats);

/******************************************************************************
 * Function: DumpFragStoreStats
 * Description:
 *     Dump stats on a fragstore to a stream
 *
 * Inputs:
 *     handle    FragStoreHandle
 *
 * Return Value:
 *     None
 *****************************************************************************/
void DumpFragStoreStats(FILE *stream, FragStoreHandle handle);

/******************************************************************************
 * Function: existsFragStore:
 * Description:
 *     Test for fragment store existence
 *
 * Inputs:
 *     FragStorePath    Path to FragStore. 
 *
 * Return Value:
 *     1 if Fragment Store Files Exists
 *     0 if Fragment Store Files are absent
 *     -1 if some Fragment Store files exist
 *****************************************************************************/
int existsFragStore(const char *FragStorePath);

/******************************************************************************
 * Function: removeFragStore:
 * Description:
 *     Remove fragment store files
 *
 * Inputs:
 *     FragStorePath    Path to FragStore. 
 *
 * Return Value:
 *     0 if success
 *     1 if failure
 *****************************************************************************/
int removeFragStore(const char *FragStorePath);

/******************************************************************************
 * Function: createFragStore:
 * Description:
 *     Allocates a new FragStore and open it for random read and append
 *
 * Inputs:
 *     FragStorePath    Path to FragStore. If NULL, then this is a memory-based store
 *     name             Arbitrary name.  Currently stored, but not used
 *     firstID          ID of first element in store
 *
 * Return Value:
 *     Handle of created FragStore.
 *****************************************************************************/


FragStoreHandle createFragStore
(const char *FragStorePath, const char *name, int64 firstID);
FragStoreHandle createPartitionedFragStore
( const char *FragStorePath, const char *name, int64 firstID, int32 numPartitions);


#if FRAGSTORE_VERSION > FRAGSTORE_VERSION_PRODUCTION
/******************************************************************************
 * Function: createIndexedFragStore:
 * Description:
 *     Allocates a new indexed FragStore and open it for random read and append
 *
 * Inputs:
 *     FragStorePath    Path to FragStore. If NULL, then this is a memory-based store
 *     name             Arbitrary name.  Currently stored, but not used
 *     firstID          ID of first element in store
 *
 * Return Value:
 *     Handle of created FragStore.
 *****************************************************************************/


FragStoreHandle createIndexedFragStore
(const char *FragStorePath, const char *name, int64 firstID);
#endif

/******************************************************************************
 * Function: resetFragStore:
 * Description:
 *     Recycles an active FragStore and open it for random read and append
 *
 * Inputs:
 *     fs               FragStoreHandle
 *     firstID          ID of first element in store
 *
 * Return Value:
 *     Handle of created FragStore.
 *****************************************************************************/


FragStoreHandle resetFragStore(FragStoreHandle fs, int64 firstID);


/******************************************************************************
 * Function: closeFragStore:
 * Description:
 *     Commits all changes and closes the store.
 *
 * Inputs:
 *     FragStoreHandle  Handle of FragStore.
 *
 * Return Value:
 *     Zero if success.
 *****************************************************************************/

int closeFragStore(FragStoreHandle);


/******************************************************************************
 * Function: openFragStore:
 * Description:
 *     Open an existing FragStore.
 *
 * Inputs:
 *     FragStorePath  Path to fragStore
 *     rw             Mode for open (r or rw)
 *
 * Return Value:
 *     FragStoreHandle of opened FragStore
 *     NULLSTOREHANDLE on failure
 *****************************************************************************/

FragStoreHandle openFragStore
( const char *FragStorePath, /* Path to directory */
  const char *rw             /* "r" or "rw" */
);


/******************************************************************************
 * Function: copyFragStore:
 * Description:
 *     Create a copy of an existing FragStore.
 *
 * Inputs:
 *     SourceFragStorePath  Path to fragStore
 *     TargetFragStorePath  Path to fragStore
 *     move    (if true, mv don't cp)
 * Return Value:
 *     0 if success
 *     1 if failure
 *****************************************************************************/

FragStoreHandle copyFragStore
( const char *SourceFragStorePath, /* Path to directory */
  const char *TargetFragStorePath,  /* Path to directory */
  int move
);

/******************************************************************************
 * Function: testOpenFragStore:
 * Description:
 *     Open an existing FragStore.
 *
 * Inputs:
 *     FragStorePath  Path to fragStore
 *     rw             Mode for open (r or rw)
 *
 * Return Value:
 *     0 if success
 *     1 if failure
 *****************************************************************************/

int testOpenFragStore
( const char *FragStorePath, /* Path to directory */
  const char *rw             /* "r" or "rw" */
);



/******************************************************************************
 * Function: loadFragStore:
 * Description:
 *     Open an existing file-based fragStore, and load its entire contents into
 * a newly created memory-based fragStore.  This is equivalent to creating
 * a new memory-based fragStore based on the firstElem of the file-based store,
 * and concatting the file-based store to the memory based store.  This is also
 * how this should be implemented.
 *
 * Inputs:
 *     FragStorePath  Path to fragStore
 *
 * Return Value:
 *     FragStoreHandle of opened FragStore
 *****************************************************************************/

FragStoreHandle loadFragStore
( const char *FragStorePath /* Path to directory */);
/******************************************************************************
 * Function: loadFragStore:
 * Description:
 *     Open an existing file-based fragStore, and load elements (first,last) of its contents into
 * a newly created memory-based fragStore. 
 *
 * Inputs:
 *     FragStorePath  Path to fragStore
 *
 * Return Value:
 *     FragStoreHandle of opened FragStore
 *****************************************************************************/

FragStoreHandle loadFragStorePartial
( const char *FragStorePath /* Path to directory */,
  int64 firstElem,
  int64 lastElem);

/******************************************************************************
 * Function: storeFragStore:
 * Description:
 *     Store an open FragStore (file or memory based) to a file.
 * This is equivalent to creating a new file-based store based on the
 * firstElem of store, and concatting store to the new file-based store.
 * This is also how this should be implemented.
 *
 * Inputs:
 *     fs             FragStore Handle
 *     FragStorePath  Path to fragStore
 *
 * Return Value:
 *     0 if OK
 *****************************************************************************/

int storeFragStore
( FragStoreHandle fs,
  const char *FragStorePath /* Path to directory */);

/******************************************************************************
 * Function: getFragStore
 * Description:
 *     Random Access read.
 *
 * Inputs:
 *     fs              A FragStore Handle
 *     indx            Index of element to get
 *     getFlags        Flag word, as per streamFlags above
 * Outputs:
 *     rs              ReadStructp to hold the retreived data
 * Return Value:
 *     Zero if OK.
 *****************************************************************************/

int getFragStore(FragStoreHandle fs, int64 indx, int32 getFlags, ReadStructp rs);


/******************************************************************************
 * Function: getNDumpFragStore
 * Description:
 *     Random Access read and dump for Oracle load
 *
 * Inputs:
 *     fs              A FragStore Handle
 *     indx            Index of element to get
 *     outfp           FILE *
 * Outputs:
 *     None
 * Return Value:
 *     Zero if OK.
 *****************************************************************************/
int getNDumpFragStore(FragStoreHandle fs, int64 index, ReadStructp rs, FILE *outfp);

/******************************************************************************
 * Function: appendFragStore
 * Description:
 *     Append a ReadStruct to the Frag Store
 *     Index of element must follow the current last element in the store.
 * Inputs:
 *     store           FragStore Handle
 *     rs              ReadStruct with the data to append
 *
 * Return Value:
 *     Zero if OK.
 *
 *****************************************************************************/

int appendFragStore(FragStoreHandle store,ReadStructp rs);  

/******************************************************************************
 * Function: appendFragStorePartition
 * Description:
 *     Append a ReadStruct to the Frag Store
 *     Index of element must follow the current last element in the store.
 * Inputs:
 *     store           FragStore Handle
 *     rs              ReadStruct with the data to append
 *     partition       Index of partition
 * Return Value:
 *     Zero if OK.
 *
 *****************************************************************************/

int appendFragStorePartition(FragStoreHandle store,ReadStructp rs, int32 partition);  



/******************************************************************************
 * Function: setFragStore
 * Description:
 *     Write a ReadStruct to the Frag Store
 *     *NOTE*  This DOES NOT REPLACE SOURCE AND SEQUENCE/QUALITY DATA!!!!!
 *     *NOTE*  Intent is to use this for read/modify/write operations on FragStore
 * Inputs:
 *     store           FragStore Handle
 *     indx           index of element to replace
 *     rs              ReadStruct with the data to append
 *
 * Return Value:
 *     Zero if OK.
 *
 *****************************************************************************/

int setFragStore(FragStoreHandle store,int64 indx,  ReadStructp rs);  



/******************************************************************************
 * Function: commitFragStore
 * Description:
 *     Commit changes to file by updating lastCommittedElem and writing header.
 *
 * Inputs:
 *     store           FragStore Handle
 *
 * Return Value:
 *     Zero if OK.
 *
 *****************************************************************************/

int commitFragStore(FragStoreHandle store);


/******************************************************************************
 * Function: concatFragStore
 * Description:
 *     Append the source FragStore to the target.  Fails if frag IDs are not
 *     consecutive, or if target is not open for append.
 *
 * Inputs:
 *     target           FragStore Handle for target
 *     source           FragStore Handle for source
 *
 * Return Value:
 *     Zero if OK.
 *
 *****************************************************************************/

int concatFragStore(FragStoreHandle target, FragStoreHandle source);

/******************************************************************************
 * Function: getLastElemFragStore
 * Description:
 *   Get the index of the last Element of the Frag Store
 *
 * Inputs:
 *     store           FragStore Handle
 * Return Value
 *     lastElem           Pointer to int to hold result.
 *
 *****************************************************************************/

int64 getLastElemFragStore(FragStoreHandle store);

/******************************************************************************
 * Function: firstElemFragStore
 * Description:
 *   Get the index of the first Element of the Frag Store
 *
 * Inputs:
 *     store           FragStore Handle
 * ReturnValue
 *     lastElem           Pointer to int to hold result.
 *
 *****************************************************************************/

int64 getFirstElemFragStore(FragStoreHandle store);


/******************************************************************************
 * Function: GetStartIndexFragStream
 * Description:
 *   Get the startindex of the frag Stream
 *
 * Inputs:
 *     stream           FragStream Handle
 * ReturnValue
 *     StartIndex           
 *
 *****************************************************************************/
int64 getStartIndexFragStream(FragStreamHandle stream);



/***********************************************************************************
 * Function: deleteFragStore
 * Description:
 *     Delete an element (or at least mark for deletion) from a FragStore.
 *
 * Inputs:
 *     store           FragStore Handle
 *     index           index of element to be deleted
 *
 * Return Value:
 *     Zero if OK.
 *
 ***********************************************************************************/

int deleteFragStore(FragStoreHandle store, int64 indx);



/*************************************    FRAGSTREAM *********************************************/


#define FRAG_S_FIXED 0x1
#define FRAG_S_SOURCE 0x2
#define FRAG_S_SEQUENCE 0x4
#define FRAG_S_ALL (FRAG_S_FIXED | FRAG_S_SOURCE | FRAG_S_SEQUENCE)


/***********************************************************************************
 * Function: openFragStream
 * Description:
 *     Open a stream on an open FragStore.
 *
 * Inputs:
 *     fs              A FragStore Handle for the store to stream
 *     buffer          A memory buffer to be used for prefetching (NOT USED)
 *     bufferSize      Extent of the buffer (NOT USED)
 *
 * Return Value:
 *     FragStreamHandle of FragStore stream
 ***********************************************************************************/
	
FragStreamHandle openFragStream(FragStoreHandle fs, /* handle to a fragment store */
				void *buffer,       /* User supplied buffer for prefetching */
				int bufferSize);

/***********************************************************************************
 * Function: resetFragStream
 * Description:
 *     reset the open Stream to a new start/stop
 *
 * Inputs:
 *     fs              A FragStream 
 *     startIndex      First element for streaming
 *     endIndex        Last element for streaming (-1 if EOF is desired)
 *
 * Return Value:
 *     0 if OK
 *****************************************************************************/
	
int resetFragStream
(FragStreamHandle fs, /* handle to a fragment store */
 int64 startIndex,   /* First record to stream */
 int64 endIndex);     /* Last record to stream (or -1)   */

/******************************************************************************
 * Function: closeFragStream
 * Description:
 *     Close a FragStream
 *
 * Inputs:
 *     fs              An open FragStream Handle 
 *
 * Return Value:
 *     Zero if OK
 *****************************************************************************/
	
int closeFragStream(FragStreamHandle fs);


/******************************************************************************
 * Function: nextFragStream
 * Description:
 *     Get the next element from the specified FragStream
 *
 * Inputs:
 *     fs              A FragStream Handle
 *     streamFlags     Flags for specifying streaming a subset of available data
 *                     (see FRAG_S_* be
 * Outputs:
 *     rs              A ReadStruct in which the next element will be deposited
 * Return Value:
 *     FragStreamHandle of FragStore stream
 *****************************************************************************/
int nextFragStream(FragStreamHandle fs, ReadStructp rs, int streamFlags);

int kNextFragStream(FragStreamHandle fs, ReadStructp rs, int streamFlags, int skipNum);

/******************************************************************************
 * WORK IN PROGRESS : JASON, OCT 2001
 * Function: updateFragStore
 * Description:
 *     Update disk FragStore from an in-memory fragStore.
 *     Only the fragments' fixed-length record is updated;
 *     the fragments' sequence, quality, and source records are not affected.
 *     Assume fragStore was copied to RAM, as with loadFragStore().
 *     Assume fragStore was modified in RAM, as with setClearRegion_ReadStruct().
 *     
 * 
 * Inputs:
 *     fs              A FragStream Handle
 *     streamFlags     Flags for specifying streaming a subset of available data
 *                     (see FRAG_S_* be
 * Outputs:
 *     rs              A ReadStruct in which the next element will be deposited
 * Return Value:
 *     FragStreamHandle of FragStore stream
 *****************************************************************************/
// int nextFragStream(FragStreamHandle fs, ReadStructp rs, int streamFlags);


#endif  //  AS_PER_FRAGSTORE_H
