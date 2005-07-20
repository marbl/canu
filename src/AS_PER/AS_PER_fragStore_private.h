
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
#ifndef AS_PER_FRAGSTORE_PRIVATE_H
#define AS_PER_FRAGSTORE_PRIVATE_H
/*************************************************************************
 Module:  AS_PER_fragStore
 Description:
     This file contains private declarations for the AS_PER_FragStore module

 *************************************************************************/

/* RCS Info
 * $Date: 2005-07-20 19:55:39 $
 * $Id: AS_PER_fragStore_private.h,v 1.5 2005-07-20 19:55:39 eliv Exp $
 * $Revision: 1.5 $
 *
 */

#include <sys/types.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"


/*   FragStore
 *   The Fragment store is compsed of an index store and two string stores.
 *   This data structure ties them all together.
 *
 *
 *   if GENERIC_STORE_USE_LONG_STRINGS is #defined, the max fragment length is extended to
 *   4Gb.
 *
 */


static uint8 checkSumBlob(char *blob, int length){
  uint8 checksum = 0;
  int i;
  for(i = 0; i < length; i++){
    checksum ^= blob[i];
  }
  return checksum;

}

typedef struct{
  StoreHandle fragStore;
#if 0
  StoreHandle sequenceStore;

  StoreHandle sourceStore;
  int spare;

#else
  int numPartitions;     // length of the following two arrays

  StoreHandle *sequenceStore;
 #ifdef i386
  int32 ptrPad1;
 #endif

  StoreHandle *sourceStore;
 #ifdef i386
  int32 ptrPad2;
 #endif

  StoreHandle *partitionStore; // a fragStore structure per partition
 #ifdef i386
  int32 ptrPad3;
 #endif

#endif

  StoreStatus status;
  int32 stillMorePadding;

  //  char storePath[FILENAME_MAX];
  char storePath[1024];
}FragStore;



/*   FragStream
 *   The Fragment store is compsed of an index store and two string stores.
 *   This data structure ties them all together for streaming.
 */
typedef struct{
  StoreHandle fragStore;
  StreamHandle fragStream;
  StreamHandle sequenceStream;
  StreamHandle sourceStream;
  StoreStatus status;
}FragStream;


/*   ShortFragmentRecord
 *    This is the data structure that is archived in the index store.
 *    The store's fixed-length db.frg file holds an array of these structs.
 *    The store's load() functions copy the db.frg file to memory in bulk.
 *    All clients load at least one of these structs into memroy,
 *    but they access its fields indirectly, through the ReadStruct API.
 *
 *   2005/07/19
 *     Remove old FRAGSTORE_VERSION stuff, and make frag store binary
 *     compatible on i386, alpha and x86_64.
 *
 *   Modified Oct 2001 by Jason
 *    Added fields to support modified clear ranges.
 *    We anticipate 3 systems will modify ranges: OVL, CNS, CGW.
 *    We must retain each changed, as well as the original, value.
 *    Adding fields to the db.frg file was not the only option, but it 
 *    was the easy option, and this functionality was needed in a hurry.
 *    To support backward compatibility with old stores, 
 *    we observe the compiler directive FRAGSTORE_VERSION.
 */

typedef struct {
  // BLANK  LINES SHOW GROUPING INTO 8-byte BLOCKS

  uint   deleted:1;
  uint   readType:8;
  uint   hasQuality:1;
  uint   numScreenMatches:16; /* number of screen matches */
  uint   hasModifiedClearRegion:1;  // never used as of Oct 2001 - Jason
  uint   hasOVLClearRegion:1; 
  uint   hasCNSClearRegion:1; 
  uint   hasCGWClearRegion:1; 
  uint   spare1:2;
  VLSTRING_SIZE_T clearRegionStart;
  VLSTRING_SIZE_T clearRegionEnd;

  VLSTRING_SIZE_T ovlRegionStart; 
  VLSTRING_SIZE_T ovlRegionEnd; 
  VLSTRING_SIZE_T cnsRegionStart; 
  VLSTRING_SIZE_T cnsRegionEnd; 

  VLSTRING_SIZE_T cgwRegionStart; 
  VLSTRING_SIZE_T cgwRegionEnd; 
  CDS_IID_t readIndex;         /* Internal ID of this read */

  CDS_UID_t accID;             /* Accession ID of this read */

  uint64 sequenceOffset;    /* Offset of the sequence/quality data in the seq Store */

  uint64 sourceOffset;      /* Offset of the source in the source, localePos, and screen Matches */

}ShortFragRecord;


#define FILEOFFSET_SIZE_FRG (48)

//#define GET_FILEOFFSET(X) (X  & FILEOFFSET_MASK)
//#define GET_FILEID(X) 


static uint64 GET_FILEOFFSET(uint64 X){
  uint64 result = X & FILEOFFSET_MASK;
  //  fprintf(stderr,"* GET_FILEOFFSET 0x" F_X64 " & 0x" F_X64 " ==> 0x" F_X64 "\n", X, FILEOFFSET_MASK, result);
  return result;

}

static uint16 GET_FILEID(uint64 X){
  uint16 result = (uint16)((X  & FILEID_MASK)>>FILEOFFSET_SIZE_FRG);
  //  fprintf(stderr,"* Get_FILEID 0x" F_X64 " ==> %d\n",	  X, result);
  return result;

}

static StoreHandle GET_FILEHANDLE(StoreHandle *x, uint64 offset){
  uint8 fileID = GET_FILEID(offset);
  return x[fileID];

}



static uint64 MERGE_FILEID_OFFSET(uint16 fileID, uint64 offset){
  uint64 fileIDPart = ((uint64)fileID << FILEOFFSET_SIZE_FRG);
  uint64 fileOffsetPart = offset & FILEOFFSET_MASK;
  uint64 result = fileIDPart | fileOffsetPart;
  //  fprintf(stderr,"* MERGE_FILEID_OFFSET(file:0x%x offset:0x" F_X64 ")  fileIDPart = 0x" F_X64 "   fileOffsetPart = 0x" F_X64 " result = 0x" F_X64 "\n",
  //	  fileID, offset, fileIDPart, fileOffsetPart, result);

  return result;

}

static uint64 SET_FILEOFFSET(uint64 X, uint64 offset){
  uint64 fileIDPart = X & FILEID_MASK;
  uint64 fileOffsetPart = offset & FILEOFFSET_MASK;

  return fileIDPart | fileOffsetPart;
}

static uint64 SET_FILEID(uint64 X, uint16 fileID){
  uint64 fileIDPart = ((uint64)fileID << FILEOFFSET_SIZE_FRG);
  uint64 fileOffsetPart = X & FILEOFFSET_MASK;

  return fileIDPart | fileOffsetPart;
}


#if VLSTRING_MAX_SIZE < 64 * 2048
 #define MAX_SEQUENCE_LENGTH (64 * 2048 - 1)
#else
 #define MAX_SEQUENCE_LENGTH (2048 * 2048 - 1)
#endif

#define MAX_SOURCE_LENGTH 512
#define MAX_SCREEN_MATCH 2048
#define MAX_SOURCE_BUFFER_LENGTH ( MAX_SOURCE_LENGTH + MAX_SCREEN_MATCH*sizeof(IntScreenMatch)+sizeof(int32) + sizeof(int64))

/* This is the (opaque) data structure that the user's code manipulates */
/* On the user side, all that is visible is a pointer to this type, which is
   case to a (void *).
*/
typedef struct {
  ShortFragRecord frag;
  uint flags;
  // These are only used by UBAC and FBAC fragments and stored in the var length source field
  CDS_UID_t localeID;
  uint32 localePosStart;
  uint32 localePosEnd;
  //
  char source[MAX_SOURCE_BUFFER_LENGTH];
  char sequence[MAX_SEQUENCE_LENGTH];
  char quality[MAX_SEQUENCE_LENGTH];
  IntScreenMatch matches[MAX_SCREEN_MATCH];
}FragRecord;

int getSourceOffset_ReadStruct(ReadStructp rs, int64 *offset);
int getSequenceOffset_ReadStruct(ReadStructp rs, int64 *offset);

int setSourceOffset_ReadStruct(ReadStructp rs, int64 offset);
int setSequenceOffset_ReadStruct(ReadStructp rs, int64 offset);



/******************************************************************************
 * Function: loadDumpFragRecord
 * Description:
 *    Load a record from a binary dump
 *
 * Inputs:
 * Outputs:
 * Return Value:
 *****************************************************************************/
int loadDumpFragRecord(FILE *infp, ShortFragRecord *fr, VA_TYPE(char ) *sequence, VA_TYPE(char) *source);
int appendDumpToFragStore(FragStoreHandle store, ShortFragRecord *fr, VA_TYPE(char) *sequence,  VA_TYPE(char) *source);

void unloadFragRecordPartition(StoreHandle seqStore, StoreHandle srcStore, FragRecord *fr, int32 getFlags);

FragStore *FragStore_myStruct(FragStoreHandle s);

#endif
