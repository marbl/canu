
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
static char CM_ID[] = "$Id: AS_PER_fragStore.c,v 1.10 2006-10-08 08:47:40 brianwalenz Exp $";

/*************************************************************************
 Module:  AS_PER_fragStore
 Description:
     This module defines the interface and implementation of the Assembler
 Fragment Store.
     The Fragment Store is built using the building blocks of the index and string
 stores in the AS_PER_genericStore module.
     Currently, the fragment store manages a directory of files containing 3 files:
        - db.frg      Index Store with Fragment data
        - db.seq      VLRecord Store with Sequence/Quality data
        - db.src      VLRecord Store with source data.
        - db.idx      Persistent VA_int32 for fragment store index by iid (not always present) 
 All memory for FragStore and FragStream data structures is managed in global arrays allocated from the 
 heap.

 Assumptions:
      Implementation of Index/VLRecord stores.
 Document:
      FragmentStore.rtf

 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStore_private.h"
#include "AS_PER_SafeIO.h"
#include "AS_UTL_Var.h"

/************** Globals ********************/
static FragStore *gFragStores = NULL;
static int gNumFragStores = 0;
static int gMaxNumFragStores = 0;

static FragStream *gFragStreams = NULL;
static int gNumFragStreams = 0;
static int gMaxNumFragStreams = 0;

extern int  FragStore_Version;

/* Compute the handle for a FragStore using pointer arithmetic */
static int FragStore_myHandle(FragStore *s){
  return s - gFragStores; /* Pointer Arithmetic */
}

/* Compute the handle for a FragStore using pointer arithmetic */
static int FragStream_myHandle(FragStream *s){
  return s - gFragStreams; /* Pointer Arithmetic */
}

/* Compute the struct for a FragStore using pointer arithmetic */
FragStore *FragStore_myStruct(FragStoreHandle s){
  if(s < 0 || s > gMaxNumFragStores)
    return NULL;
  return gFragStores + s; /* Pointer Arithmetic */
}

/* Compute the struct for a FragStore using pointer arithmetic */
static FragStream *FragStream_myStruct(FragStreamHandle s){
  if(s < 0 || s > gMaxNumFragStreams)
    return NULL;
  return gFragStreams + s; /* Pointer Arithmetic */
}


/***************************************************************************
 Allocation functions
***************************************************************************/

/***********************************************************************************
 * Function allocateFragStore
 ***********************************************************************************/
static FragStore  *allocateFragStore(void){
  int i;
  FragStore *fs = NULL;
#ifdef DEBUG
  fprintf(stderr,"*** allocateFragStore\n");
#endif

  for(i = 0; i < gMaxNumFragStores; i++){
    if(gFragStores[i].status == UnAllocatedStore){
      gFragStores[i].status = UnInitializedStore;
      gNumFragStores++;
      fs = gFragStores + i;
      break;
    }
  }
  if(fs == NULL){
    /* We haven't found an available entry */
    if(gNumFragStores == 0)
      gMaxNumFragStores = 16;
    else
      gMaxNumFragStores *= 2;

    gFragStores = (FragStore *)
      realloc(gFragStores, gMaxNumFragStores * sizeof(FragStore));
    for(i = gNumFragStores; i < gMaxNumFragStores; i++)
      gFragStores[i].status = UnAllocatedStore;

    fs = gFragStores + gNumFragStores++;
  }
  fs->status = ActiveStore;

  return fs;
}


/***********************************************************************************
 * Function allocateFragStream
 ***********************************************************************************/
static FragStream  *allocateFragStream(FragStoreHandle sh, void *buffer, int32 bufferSize){
  int i;
  FragStore *s = FragStore_myStruct(sh);
  // int sourceOffset, sequenceOffset;
  FragStream *fs = NULL;

#ifdef DEBUG
  fprintf(stderr,"*** allocateFragStream\n");
#endif
  for(i = 0; i < gMaxNumFragStreams; i++){
    if(gFragStreams[i].status == UnAllocatedStore){
      gFragStreams[i].status = UnInitializedStore;
      gNumFragStreams++;
      fs = FragStream_myStruct(i);
      break;
    }
  }
  if(fs == NULL){
  /* We haven't found an available entry */
    if(gNumFragStreams == 0)
      gMaxNumFragStreams = 16;
    else
      gMaxNumFragStreams *= 2;

    gFragStreams = (FragStream *) 
      realloc(gFragStreams, gMaxNumFragStreams * sizeof(FragStore));
    for  (i = gNumFragStreams;  i < gMaxNumFragStreams;  i ++)
      gFragStreams [i] . status = UnAllocatedStore;

    fs = FragStream_myStruct(gNumFragStreams++);
  }
  fs->status = UnInitializedStore;

  /* Start the encapsulated stores a-streaming... */

  fs->fragStore = sh;
  fs->fragStream = openStream(s->fragStore,  NULL, 0);
 
  /* If we are not streaming from the beginning of the store, we have to
     get the startIndex-th element of the fragStore to get the start points
     for the other two streaming operations
  */
  fs->sequenceStream = openStream(s->sequenceStore[0], NULL, 0);
  fs->sourceStream   = openStream(s->sourceStore[0]  , NULL, 0);
  fs->status = ActiveStore;
  
  return fs ;
}

/*****************************************************************************/
#include "AS_PER_encodeSequenceQuality.h"


/***********************************************************************************
 * Function unloadFragRecord
 *    Utility routine for unloading a FragRecord.  used by getFragStore and nextFragStream.
 ***********************************************************************************/
void unloadFragRecord(FragStore *myStore, FragRecord *fr, int32 getFlags){
  static char encodeBuffer[MAX_SEQUENCE_LENGTH];
  uint16 localeLength=0, screenLength;
  VLSTRING_SIZE_T actualLength;

  /*  looks like checking AS_FA_READ and AS_FA_SHREDDED FT isolated VR */
  switch(fr->frag.readType){
  case AS_READ:
  case AS_B_READ:
  case AS_EXTR:
  case AS_TRNR:
    localeLength = 0;
    break;
  case AS_STS:
  case AS_EBAC:
  case AS_LBAC:
  case AS_FULLBAC:
  case AS_BACTIG:
    localeLength = 1 + sizeof(uint64);
    break;
  case AS_UBAC:
  case AS_FBAC:
    localeLength = 1 + sizeof(uint64) + 2 * sizeof(uint32);
    break;
  default:
    fprintf(stderr,"fr->frag.readType=%c\n",(char) fr->frag.readType);
    assert(0);
  }

 screenLength = fr->frag.numScreenMatches * sizeof(IntScreenMatch);

#ifdef DEBUG
 fprintf(stderr,"* screenLength = %u localeLength = %u\n",
	 screenLength, localeLength);
#endif
  if((getFlags & FRAG_S_SOURCE) ||
     fr->frag.numScreenMatches > 0 || 
     localeLength != 0){

     getVLRecordStore(GET_FILEHANDLE(myStore->sourceStore, fr->frag.sourceOffset), 
		     GET_FILEOFFSET(fr->frag.sourceOffset), fr->source, 
		     (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, 
		     &actualLength);
    
#ifdef DEBUG
    fprintf(stderr,"* Read " F_VLS " chars %c\n", actualLength, fr->source[actualLength - localeLength - screenLength]);
#endif
    //    fr->source[actualLength - localeLength - screenLength] = '\0';
    if(localeLength > 0){
      size_t offset = actualLength - localeLength + 1;
#ifdef DEBUG
      fprintf(stderr,"* Reading Locale at offset " F_SIZE_T "\n",
	      offset);
#endif
      memcpy(&fr->localeID, fr->source + offset, sizeof(uint64));
      offset += sizeof(uint64);
      if(localeLength > sizeof(uint64) + 1){
	memcpy(&fr->localePosStart, fr->source + offset, sizeof(uint32));
	offset += sizeof(int32);
	memcpy(&fr->localePosEnd, fr->source + offset, sizeof(uint32));
      }
    }
      if(screenLength > 0){
	size_t offset = actualLength - screenLength - localeLength + 1;
#ifdef DEBUG
      fprintf(stderr,"* Reading screenMatches at offset " F_SIZE_T "\n",
	      offset);
#endif
	memcpy(fr->matches, fr->source + offset, fr->frag.numScreenMatches * sizeof(IntScreenMatch));
	{ // What the hell, link the suckers together
	  int i;
	  for(i = 0; i < fr->frag.numScreenMatches - 1; i++){
	    //fr->matches[i].next = fr->matches + i + 1;
	    fr->matches[i].next = &(fr->matches[i + 1]);
	  }
	  fr->matches[i].next = NULL;
	}
      }
  }
  if(getFlags & FRAG_S_SEQUENCE){
    getVLRecordStore(GET_FILEHANDLE(myStore->sequenceStore, fr->frag.sequenceOffset), 
		     GET_FILEOFFSET(fr->frag.sequenceOffset), encodeBuffer, (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, &actualLength);
    
    encodeBuffer[actualLength] = '\0';
    decodeSequenceQuality(encodeBuffer, actualLength, fr->sequence, fr->quality, fr->frag.hasQuality);
  }
}

/***********************************************************************************
 * Function unloadFragRecordPartition
 *    Utility routine for unloading a FragRecord.  used by getFragStore and nextFragStream.
 ***********************************************************************************/
void unloadFragRecordPartition(StoreHandle seqStore, StoreHandle srcStore, FragRecord *fr, int32 getFlags){
  static char encodeBuffer[MAX_SEQUENCE_LENGTH];
  uint16 localeLength=0, screenLength;
  VLSTRING_SIZE_T actualLength;

  /* FT isolated VR */
  switch(fr->frag.readType){
  case AS_READ:
  case AS_B_READ:
  case AS_EXTR:
  case AS_TRNR:
    localeLength = 0;
    break;
  case AS_STS:
  case AS_EBAC:
  case AS_LBAC:
  case AS_FULLBAC:
  case AS_BACTIG:
    localeLength = 1 + sizeof(uint64);
    break;
  case AS_UBAC:
  case AS_FBAC:
    localeLength = 1 + sizeof(uint64) + 2 * sizeof(uint32);
    break;
  default:
    assert(0);
  }

 screenLength = fr->frag.numScreenMatches * sizeof(IntScreenMatch);

#ifdef DEBUG
 fprintf(stderr,"* screenLength = %u localeLength = %u\n",
	 screenLength, localeLength);
#endif
  if((getFlags & FRAG_S_SOURCE) ||
     fr->frag.numScreenMatches > 0 || 
     localeLength != 0){

     getVLRecordStore(srcStore,
		     GET_FILEOFFSET(fr->frag.sourceOffset), fr->source, 
		     (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, 
		     &actualLength);
    
#ifdef DEBUG
    fprintf(stderr,"* Read " F_VLS " chars %c\n", actualLength, fr->source[actualLength - localeLength - screenLength]);
#endif
    //    fr->source[actualLength - localeLength - screenLength] = '\0';
    if(localeLength > 0){
      size_t offset = actualLength - localeLength + 1;
#ifdef DEBUG
      fprintf(stderr,"* Reading Locale at offset " F_SIZE_T "\n",
	      offset);
#endif
      memcpy(&fr->localeID, fr->source + offset, sizeof(uint64));
      offset += sizeof(uint64);
      if(localeLength > sizeof(uint64) + 1){
	memcpy(&fr->localePosStart, fr->source + offset, sizeof(uint32));
	offset += sizeof(int32);
	memcpy(&fr->localePosEnd, fr->source + offset, sizeof(uint32));
      }
    }
      if(screenLength > 0){
	size_t offset = actualLength - screenLength - localeLength + 1;
#ifdef DEBUG
      fprintf(stderr,"* Reading screenMatches at offset " F_SIZE_T "\n",
	      offset);
#endif
	memcpy(fr->matches, fr->source + offset, fr->frag.numScreenMatches * sizeof(IntScreenMatch));
	{ // What the hell, link the suckers together
	  int i;
	  for(i = 0; i < fr->frag.numScreenMatches - 1; i++){
	    //fr->matches[i].next = fr->matches + i + 1;
	    fr->matches[i].next = &(fr->matches[i + 1]);
	  }
	  fr->matches[i].next = NULL;
	}
      }
  }
  if(getFlags & FRAG_S_SEQUENCE){
    getVLRecordStore(seqStore, 
		     GET_FILEOFFSET(fr->frag.sequenceOffset), encodeBuffer, (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, &actualLength);
    
    encodeBuffer[actualLength] = '\0';
    decodeSequenceQuality(encodeBuffer, actualLength, fr->sequence, fr->quality, fr->frag.hasQuality);
  }
}


/***********************************************************************************
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
 ***********************************************************************************/
int existsFragStore(const char *FragStorePath){
  char frgFile[FILENAME_MAX];
  char srcFile[FILENAME_MAX];
  char seqFile[FILENAME_MAX];
  char parFile[FILENAME_MAX];
  FILE *fp, *parfp;
  int fileCount = 0;

  AssertPtr(FragStorePath);

  sprintf(parFile,"%s/db.par", FragStorePath);
  if((parfp = fopen(parFile,"r")) == NULL){

    sprintf(frgFile,"%s/db.frg",FragStorePath);
    sprintf(srcFile,"%s/db.src",FragStorePath);
    sprintf(seqFile,"%s/db.seq",FragStorePath);



    if((fp = fopen(frgFile,"r"))){
      fileCount++;
      fclose(fp);
    }
    if((fp = fopen(srcFile,"r"))){
      fileCount++;
      fclose(fp);
    }
    if((fp = fopen(seqFile,"r"))){
      fileCount++;
      fclose(fp);
    }
    

    if(fileCount == 0)
      return 0;
    else if(fileCount == 3)
      return 1;
    else
      return -1;
  }else{
    int i;
    int numPartitions;
    fscanf(parfp,"%d",&numPartitions);
    assert(numPartitions > 1);
    fclose(parfp);

    sprintf(frgFile,"%s/db.frg", FragStorePath);
    if(((fp = fopen(frgFile, "r")) )){
      fileCount++;
      fclose(fp);
    }

    for(i = 0; i < numPartitions; i++){
      // Check if the required files are all there.
      sprintf(seqFile,"%s/db.seq.%d", FragStorePath,i);
      sprintf(srcFile,"%s/db.src.%d", FragStorePath,i);
      sprintf(frgFile,"%s/db.src.%d", FragStorePath,i);

      if((fp = fopen(srcFile, "r")) != NULL){
	fileCount++;
	fclose(fp);
      }
      if((fp = fopen(seqFile, "r")) != NULL){
	fileCount++;
	fclose(fp);
      }
      if((fp = fopen(frgFile, "r")) != NULL){
	fileCount++;
	fclose(fp);
      }
    }
    
    if(fileCount == 0)
      return 0;
    else if(fileCount == (1 + 3*numPartitions))
      return 1;
    else
      return -1;

  }

}

/***********************************************************************************
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
 ***********************************************************************************/
int removeFragStore(const char *FragStorePath){
  char cmd[FILENAME_MAX * 3];

  AssertPtr(FragStorePath);
  sprintf(cmd,"rm -rf %s/db.frg* %s/db.src* %s/db.seq* %s/db.par", FragStorePath, FragStorePath, FragStorePath, FragStorePath);
  if(system(cmd) != 0) assert(0);
  return 0;
}


//#define DEBUG_OPEN
/****************************************************************************
 * testOpenFragStore:
 * Description:
 *    See if all of fragstore constituents are present and openable.
 *
 * Inputs:
 *     FragStorePath    Path to FragStore. 
 *
 * Return Value:
 *     0 if success
 *     1 if failure
 ***********************************************************************************/
      
int testOpenFragStore(const char *FragStorePath, const char *rw){
  DIR *dbDir;
  char parbuffer[FILENAME_MAX];
  char frgbuffer[FILENAME_MAX];
  char srcbuffer[FILENAME_MAX];
  char seqbuffer[FILENAME_MAX];
  char *errorFile;
  FILE *fragFile, *sourceFile, *seqFile, *parFile;

  if(FragStorePath){
    dbDir = opendir(FragStorePath);
    if(dbDir == NULL){
      fprintf(stderr,"*** Couldn't find Fragment Store directory %s\n", FragStorePath);
      return 1; // failure
    }
    closedir(dbDir);
#ifdef DEBUG_OPEN
    fprintf(stderr," Verified existence of directory %s\n", FragStorePath);
#endif
  }else{
    fprintf(stderr,"Can't Open a Memory-based frag Store....\n");
    return 1; // failure
  }
  sourceFile = fragFile = seqFile = parFile = NULL;
  sprintf(parbuffer,"%s/db.par", FragStorePath);
  if((parFile = fopen(parbuffer,rw)) == NULL){

   fprintf(stderr,"* Unpartitioned frag store no file %s\n", parbuffer);
   fflush(stderr);
   sprintf(frgbuffer,"%s/db.frg", FragStorePath);
   sprintf(seqbuffer,"%s/db.seq", FragStorePath);
   sprintf(srcbuffer,"%s/db.src", FragStorePath);
   if(((fragFile = fopen(frgbuffer, rw)) != NULL) &&
      ((sourceFile = fopen(srcbuffer, rw)) != NULL) &&
      ((seqFile = fopen(seqbuffer, rw)) != NULL)
      ){
#ifdef DEBUG_OPEN
     fprintf(stderr,"** Successfully verified openability of db files\n");
#endif
     fclose(fragFile);
     fclose(sourceFile);
     fclose(seqFile);
     return 0; // success
   }else{
     if(fragFile){
       fclose(fragFile);
       if(sourceFile){
	 fclose(sourceFile);
	 if(seqFile){
	   assert(0);
	 }else{
	   errorFile = seqbuffer;
	 }
       }else {
	 errorFile = srcbuffer;
       }
     }else{
       errorFile = frgbuffer;
     }


     fprintf(stderr," * Couldn't open file %s\n", errorFile);
     return 1; // failure
   }
 }else{
   int i;
   int numPartitions;
   fscanf(parFile,"%d",&numPartitions);
   assert(numPartitions > 1);
   fclose(parFile);

   sprintf(frgbuffer,"%s/db.frg", FragStorePath);
   if(((fragFile = fopen(frgbuffer, rw)) == NULL)){
     fprintf(stderr," * Couldn't open file %s\n", frgbuffer);
     return 1;
   }
   fclose(fragFile);

   for(i = 0; i < numPartitions; i++){

   sprintf(seqbuffer,"%s/db.seq.%d", FragStorePath,i);
   sprintf(srcbuffer,"%s/db.src.%d", FragStorePath,i);

   if(((sourceFile = fopen(srcbuffer, rw)) != NULL) &&
       ((seqFile = fopen(seqbuffer, rw)) != NULL)    ){
#ifdef DEBUG_OPEN
     fprintf(stderr,"** Successfully verified openability of db files\n");
#endif
     fclose(sourceFile);
     fclose(seqFile);
   }else{
       if(sourceFile){
	 fclose(sourceFile);
	 errorFile = seqbuffer;
       } else {
	 errorFile = srcbuffer;
       }

     fprintf(stderr," * Couldn't open file %s\n", errorFile);
     return 1; // failure
   }
   }
     return 0;
 }
}
      
      
/* openFragStore:
      Open an existing fragStore.  
*/
FragStoreHandle openFragStoreCommon
( const char *FragStorePath, /* Path to directory */
  const char *rw,             /* "r" or "r+" */
  const int needIndex){

  char frgbuffer[FILENAME_MAX];
  char parbuffer[FILENAME_MAX];
  char srcbuffer[FILENAME_MAX];
  char seqbuffer[FILENAME_MAX];
  // char *nameBuffer = NULL;
  // int myStoreIndex;
  FragStore *myStore;
  StoreStat  status;

#ifdef DEBUG_OPEN
  fprintf(stderr,"*** openFragStore\n");
#endif


  myStore = allocateFragStore();

  AssertPtr(myStore);

  /* First we need to check that the directory and the constituent filed exist */
   
 if(testOpenFragStore(FragStorePath, rw)){
   return NULLSTOREHANDLE;
 }

 sprintf(parbuffer,"%s/db.par", FragStorePath);
 sprintf(frgbuffer,"%s/db.frg", FragStorePath);

#ifdef DEBUG_OPEN
 fprintf(stderr," *** Opening Stores ***\n");
#endif
 myStore->fragStore = openStore(frgbuffer,rw);

 {
   FILE *parfp = fopen(parbuffer,"r");
   if(parfp){
     fscanf(parfp,"%d",&myStore->numPartitions);
     fclose(parfp);
   }else{
     myStore->numPartitions = 1;
   }
 }

 fprintf(stderr,"* Number of partitions in %s is %d\n", parbuffer,myStore->numPartitions);
 myStore->sourceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
 myStore->sequenceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
 myStore->partitionStore = NULL;

 if(myStore->numPartitions == 1){
   sprintf(seqbuffer,"%s/db.seq", FragStorePath);
   sprintf(srcbuffer,"%s/db.src", FragStorePath);

   myStore->sourceStore[0] = openStore(srcbuffer, rw);
   myStore->sequenceStore[0] = openStore(seqbuffer, rw);

 }else{
   int i;
   for(i = 0; i < myStore->numPartitions; i++){

     sprintf(seqbuffer,"%s/db.seq.%d", FragStorePath,i);
     sprintf(srcbuffer,"%s/db.src.%d", FragStorePath,i);

     myStore->sourceStore[i] = openStore(srcbuffer, rw);
     myStore->sequenceStore[i] = openStore(seqbuffer, rw);

   }

 }
#ifdef DEBUG_OPEN
 fprintf(stderr," *** Done! Opening Stores ***\n");
#endif
 statsStore (myStore -> fragStore, & status);
 FragStore_Version = status . version;
 if  (status . version != FRAGSTORE_VERSION)
     {
      fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf (stderr, "ERROR!!!  The fragStore code version you are running (%2d) is\n"
	               "incompatible with the version of the code that was used to create\n"
	               "the store you are trying to open (%2d).  You must recreate the\n"
	               "fragStore with the current version of the code.  Sorry for this\n"
	               "inconvenience. Call Saul with questions/complaints/invective.\n",
                  FRAGSTORE_VERSION, status . version);
      fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      return NULLSTOREHANDLE;
     }

 return FragStore_myHandle(myStore);
}

FragStoreHandle openFragStore
( const char *FragStorePath, /* Path to directory */
  const char *rw){             /* "r" or "r+" */
  return openFragStoreCommon(FragStorePath, rw, FALSE);
}

/******************************************************************************
 * Function: copyFragStore:
 * Description:
 *     Create a copy of an existing FragStore.
 *
 * Inputs:
 *     SourceFragStorePath  Path to fragStore
 *     TargetFragStorePath  Path to fragStore
 *
 * Return Value:
 *     0 if success
 *     1 if failure
 *****************************************************************************/

int copyFragStore
( const char *SourceFragStorePath, /* Path to directory */
  const char *TargetFragStorePath,  /* Path to directory */
  int move  // if true, move, don't copy
  ){
  char cmd[FILENAME_MAX * 30];
  char buffer[2048];

  if(!SourceFragStorePath)
    return 1;

  if(testOpenFragStore(SourceFragStorePath,"r")){
    return 1;
  }
    /* Create a directory if needed */
    {
      DIR *dbDir = opendir(TargetFragStorePath);
      if (dbDir) {
        closedir(dbDir);
      } else {
        if(mkdir(TargetFragStorePath, S_IRWXU | S_IRWXG | S_IROTH)){
          sprintf(buffer,"%sFragStore: Failure to create directory %s", (move?"move":"copy"),TargetFragStorePath);
          perror(buffer);
          return 1;
        }
      }
    }

  sprintf(cmd,"%s %s/db.frg %s/db.src %s/db.seq %s", (move?"mv":"cp"),SourceFragStorePath, SourceFragStorePath, SourceFragStorePath, TargetFragStorePath);
  return(system(cmd));

}


/* loadFragStore:
      Load an existing fragStore into a memory based fragStore
*/
FragStoreHandle loadFragStorePartial
( const char *FragStorePath, /* Path to directory */
  int64 firstElem,
  int64 lastElem)
{
  DIR *dbDir;
  FILE *fragFile, *sourceFile, *seqFile;
  char frgbuffer[2048];
  char srcbuffer[2048];
  char seqbuffer[2048];
  char idxbuffer[2048];
  char *errorFile;
  // char *nameBuffer = NULL;
  // int myStoreIndex;
  FragStore *myStore;
  StoreStat  status;

#ifdef DEBUG_OPEN
  fprintf(stderr,"*** loadFragStore\n");
#endif


  /* First we need to check that the directory and the constituent filed exist */

  if(FragStorePath){
  dbDir = opendir(FragStorePath);
  if (dbDir == NULL) {
    fprintf(stderr,"*** Couldn't find Fragment Store directory %s\n", FragStorePath);
    return NULLSTOREHANDLE;
  }
  closedir(dbDir);
#ifdef DEBUG_OPEN
  fprintf(stderr," Verified existence of directory %s\n", FragStorePath);
#endif
  }else{
    fprintf(stderr,"Can't Open a Memory-based frag Store....\n");
    exit(1);
  }
 sourceFile = fragFile = seqFile = NULL;
 sprintf(frgbuffer,"%s/db.frg", FragStorePath);
 sprintf(seqbuffer,"%s/db.seq", FragStorePath);
 sprintf(srcbuffer,"%s/db.src", FragStorePath);
 if(((fragFile = fopen(frgbuffer, "r")) != NULL) &&
    ((sourceFile = fopen(srcbuffer, "r")) != NULL) &&
    ((seqFile = fopen(seqbuffer, "r")) != NULL) 
    ){
#ifdef DEBUG_OPEN
   fprintf(stderr,"** Successfully verified openability of db files\n");
#endif
   fclose(fragFile);
   fclose(sourceFile);
   fclose(seqFile);
 }else{
   if(fragFile){
     fclose(fragFile);
     if(sourceFile){
       fclose(sourceFile);
       if(seqFile){
	 fclose(seqFile);
	 errorFile = idxbuffer;
       }else{
	 errorFile = seqbuffer;
       }
     }else {
       errorFile = srcbuffer;
     }
   }else{
     errorFile = frgbuffer;
   }

   fprintf(stderr," * Couldn't open file %s\n", errorFile);
   return NULLSTOREHANDLE;
 }
   
if( firstElem == STREAM_FROMSTART &&
    lastElem  == STREAM_UNTILEND){
#ifdef DEBUG_OPEN
 fprintf(stderr," *** Opening Stores ***\n");
#endif
  myStore = allocateFragStore();
  AssertPtr(myStore);
 myStore->fragStore = loadStore(frgbuffer);
 myStore->numPartitions = 1;
 myStore->partitionStore = NULL;    // added by Art.
 myStore->sourceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
 myStore->sequenceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
 myStore->sourceStore[0] = loadStore(srcbuffer);
 myStore->sequenceStore[0] = loadStore(seqbuffer);
}else{ /* LoadPartial */
  ReadStructp frag = new_ReadStruct();
  int64 firstSourceElem, lastSourceElem;
  int64 firstSeqElem, lastSeqElem;
  ShortFragRecord sfr;
  StoreHandle sh = openStore(frgbuffer,"r");
  StoreStat stats;
  int64 lastStoreElem;

  statsStore(sh,&stats);
  lastStoreElem = stats.lastElem;

  myStore = allocateFragStore();
  myStore->numPartitions = 1;
  myStore->partitionStore = NULL;    // added by Art.
  myStore->sourceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
  myStore->sequenceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
  myStore->fragStore = loadStorePartial(frgbuffer,firstElem, lastElem);
 // Fetch firstElem and lastElem
  getIndexStore(sh, firstElem,&sfr);
  firstSourceElem = sfr.sourceOffset;
  firstSeqElem = sfr.sequenceOffset;
  //  fprintf(stderr,"* LoadFragStorePartial (" F_S64 "," F_S64 ") lastStoreElem " F_S64 "\n", firstElem, lastElem, lastStoreElem);
  if(lastElem == STREAM_UNTILEND ||
     lastElem == lastStoreElem){
    lastSourceElem = STREAM_UNTILEND;
    lastSeqElem = STREAM_UNTILEND;
  }else{
    fprintf(stderr,"* loading element " F_S64 " to get source/seq offsets\n",
            lastElem + 1);
    getIndexStore(sh, lastElem + 1,&sfr);
    lastSourceElem = sfr.sourceOffset - 1;
    lastSeqElem = sfr.sequenceOffset - 1;
  }
  delete_ReadStruct(frag);
  closeStore(sh);
  //  fprintf(stderr,"* LoadFragStorePartial Source (" F_S64 "," F_S64 ")\n", firstSourceElem, lastSourceElem);
  //  fprintf(stderr,"* LoadFragStorepartial Seq (" F_S64 "," F_S64 ")\n", firstSeqElem, lastSeqElem);

 // Extract first and last source/seq elem
 myStore->sourceStore[0] = loadStorePartial(srcbuffer, firstSourceElem, lastSourceElem);
 myStore->sequenceStore[0] = loadStorePartial(seqbuffer, firstSeqElem, lastSeqElem);
}
#ifdef DEBUG_OPEN
 fprintf(stderr," *** Done! Opening Stores ***\n");
#endif
 statsStore (myStore -> fragStore, & status);
 FragStore_Version = status . version;
 // if  (status . version < FRAGSTORE_VERSION)  // changed by Jason Oct 2001
  if  (status . version != FRAGSTORE_VERSION)
     {
      fprintf (stderr, "\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a");
      fprintf (stderr,
               "ERROR!!!  Code version = %d, save store version = %d\n",
               FRAGSTORE_VERSION, status . version);
      fprintf (stderr, "\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a\a");
      return NULLSTOREHANDLE;
     }

 return FragStore_myHandle(myStore);
}



FragStoreHandle loadFragStore
( const char *FragStorePath) /* Path to directory */
{
  return loadFragStorePartial(FragStorePath, STREAM_FROMSTART, STREAM_UNTILEND);
}

/***********************************************************************/
/*  Private create functions */
/* createFragStoreCommon:
   Description: 
   Creates a new fragment store in the directory specified by FragStorePath.
*/  
FragStoreHandle createFragStoreCommon
( const char *FragStorePath, const char *name, int64 firstID, int indexed, int numPartitions){
  char buffer[2048];
  char *nameBuffer = NULL;
  // int myStoreIndex;
  FragStore *myStore;


#ifdef DEBUG
  fprintf(stderr,"*** createFragStore\n");
#endif


  assert(FragStorePath || numPartitions == 1);

  myStore = allocateFragStore();
  myStore->numPartitions = numPartitions;
  myStore->sourceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
  myStore->sequenceStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));
  myStore->partitionStore = NULL;

  AssertPtr(myStore);

  /* First we need to create a directory for the database */


  if(FragStorePath){
#ifdef DEBUG
    fprintf(stderr," Creating file-based fragment store\n");
#endif
    DIR *dbDir = opendir(FragStorePath);

    if (dbDir) {
      closedir(dbDir);
    } else {
      if (mkdir(FragStorePath, S_IRWXU | S_IRWXG | S_IROTH)) {
	sprintf(buffer,"createFragStore: Failure to create directory %s",
                FragStorePath);
	perror(buffer);
	exit(1);
      }
    }

    nameBuffer = buffer;
  } else{
    nameBuffer = NULL;
#ifdef DEBUG
    fprintf(stderr," Creating memory-based fragment store\n");
#endif
  }

 
  if  (FragStorePath)
    sprintf(buffer,"%s/db.frg", FragStorePath);
   
  myStore->fragStore = createIndexStore(nameBuffer,name,
					sizeof(ShortFragRecord), FRAGSTORE_VERSION, firstID);

  if(myStore->numPartitions == 1){
  if  (FragStorePath)
    sprintf(buffer,"%s/db.seq", FragStorePath);
  myStore->sequenceStore[0] = createVLRecordStore(nameBuffer,name,MAX_SEQUENCE_LENGTH, 1);
  if  (FragStorePath)
    sprintf(buffer,"%s/db.src", FragStorePath);
  myStore->sourceStore[0] = createVLRecordStore(nameBuffer,name,MAX_SOURCE_LENGTH, 1);
  }else{
    int i;
    FILE *fpPar;

    myStore->partitionStore = (StoreHandle *)malloc(myStore->numPartitions * sizeof(StoreHandle));

    for(i = 0; i < myStore->numPartitions; i++){
      sprintf(buffer,"%s/db.frg.%d", FragStorePath,i);
      myStore->partitionStore[i] = createIndexStore(nameBuffer,name,
					sizeof(ShortFragRecord), FRAGSTORE_VERSION, 1);
      sprintf(buffer,"%s/db.seq.%d", FragStorePath,i);
      myStore->sequenceStore[i] = createVLRecordStore(nameBuffer,name,MAX_SEQUENCE_LENGTH, 1);

      sprintf(buffer,"%s/db.src.%d", FragStorePath,i);
      myStore->sourceStore[i] = createVLRecordStore(nameBuffer,name,MAX_SOURCE_LENGTH, 1);

    }
    sprintf(buffer,"%s/db.par", FragStorePath);
    fpPar = fopen(buffer,"w");
    fprintf(fpPar,"%d\n", myStore->numPartitions);
    fclose(fpPar);
  }
  return FragStore_myHandle(myStore);

}
/***********************************************************************/
/*  Public create functions */

FragStoreHandle createPartitionedFragStore
( const char *FragStorePath, const char *name, int64 firstID, int32 numPartitions){
  return createFragStoreCommon(FragStorePath, name, firstID, FALSE, numPartitions);
}

FragStoreHandle createFragStore
( const char *FragStorePath, const char *name, int64 firstID){
  return createFragStoreCommon(FragStorePath, name, firstID, FALSE, 1);
}

FragStoreHandle createIndexedFragStore
( const char *FragStorePath, const char *name, int64 firstID){
  return createFragStoreCommon(FragStorePath, name, firstID, TRUE,1);
}

/***********************************************************************/

FragStoreHandle resetFragStore(FragStoreHandle fs, int64 firstID){
  FragStore *myStore = FragStore_myStruct(fs);

#ifdef DEBUG
  fprintf(stderr,"*** restFragStore\n");
#endif
  assert(myStore->numPartitions == 1);
  myStore->fragStore = resetIndexStore(myStore->fragStore, firstID);

  myStore->sequenceStore[0] = resetVLRecordStore(myStore->sequenceStore[0]);

  myStore->sourceStore[0] = resetVLRecordStore(myStore->sourceStore[0]);

 return FragStore_myHandle(myStore);


}

/***********************************************************************/
// #define DEBUG_CONCAT 

int concatFragStore(FragStoreHandle targetH, FragStoreHandle sourceH){
 ReadStructp myRead;
  FragStreamHandle sourceStream;
  FragStore *target = FragStore_myStruct(targetH);
  FragStore *source = FragStore_myStruct(sourceH);
  StoreStat stats;
  int64 sequenceOffset, sourceOffset;

  assert(target->numPartitions ==1);
  assert(source->numPartitions ==1);
#ifdef DEBUG_CONCAT
  fprintf(stderr,
          "*** concatFragStore source:%d status:%d  target:%d status:%d\n",
	  sourceH, source->status, targetH, target->status);
#endif

  assert(source->status == ActiveStore &&
	 target->status == ActiveStore);

  sourceStream = openStream(source->fragStore, NULL, 0);
  /* Simply concatenate the two string stores.  We'll patch
     the offsets in the index store */

  statsStore(target->sequenceStore[0] , &stats);
  sequenceOffset = stats.lastElem;
  statsStore(target->sourceStore[0] , &stats);
  sourceOffset = stats.lastElem;

#ifdef DEBUG_CONCAT
  fprintf(stderr,"\tTarget sequenceOffset:" F_S64 " sourceoffset:" F_S64 "\n",
	  sequenceOffset, sourceOffset);
  fprintf(stderr,"Concat Sequence Store\n");
#endif

  concatStore(target->sequenceStore[0], source->sequenceStore[0]);
#ifdef DEBUG_CONCAT
  fprintf(stderr,"Concat Source Store\n");
#endif
  concatStore(target->sourceStore[0], source->sourceStore[0]);

  myRead = new_ReadStruct();
#ifdef DEBUG_CONCAT
  fprintf(stderr,"Appending Frag Records\n");
#endif
  while(nextStream(sourceStream, myRead)){
    int64 oldSourceOffset, oldSequenceOffset;
    getSourceOffset_ReadStruct(myRead, &oldSourceOffset);
    getSequenceOffset_ReadStruct(myRead, &oldSequenceOffset);
#ifdef DEBUG_CONCAT
  fprintf(stderr,"Appending Record source:" F_S64 "," F_S64 " sequence:" FS_64 "," F_S64 "\n",
	  oldSourceOffset, sourceOffset + oldSourceOffset,
	  oldSequenceOffset, sequenceOffset + oldSequenceOffset);
	  
#endif
    setSourceOffset_ReadStruct(myRead, sourceOffset + oldSourceOffset);
    setSequenceOffset_ReadStruct(myRead, sequenceOffset + oldSequenceOffset);
    appendIndexStore(target->fragStore, myRead);
  }

  closeStream(sourceStream);
  delete_ReadStruct(myRead);
  return(0);
}


/***********************************************************************/
/* closeFragStore
      Close Store: commit all changes and close the store.
*/
int closeFragStore(FragStoreHandle f){
  FragStore *myStore = FragStore_myStruct(f);
  int i;
#ifdef DEBUG
  fprintf(stderr,"*** closeFragStore\n");
#endif
  myStore->status = UnAllocatedStore;
  gNumFragStores--;
  closeStore(myStore->fragStore);

  for(i = 0; i < myStore->numPartitions; i++){
    if(myStore->partitionStore){
      closeStore(myStore->partitionStore[i]);
    }
    closeStore(myStore->sequenceStore[i]);
    closeStore(myStore->sourceStore[i]);
  }
  free(myStore->sequenceStore);
  free(myStore->sourceStore);
  return(0);

}

/***********************************************************************/
int commitFragStore(FragStoreHandle store){
  FragStore *myStore = FragStore_myStruct(store);

#ifdef DEBUG
  fprintf(stderr,"*** commitFragStore\n");
#endif
  assert(myStore->numPartitions == 1);

  commitStore(myStore->fragStore);
  commitStore(myStore->sequenceStore[0]);
  commitStore(myStore->sourceStore[0]);
  return(0);
}




/***********************************************************************/
/* Private interfaces */
/* getFragStoreCommon  
   Random Access Read.
   The data is read into the previously allocated ReadStruct rs.
   The types of data read can be restricted using the getFlags.
*/
	
int getFragStoreCommon(FragStoreHandle fs, int64 index, int32 getFlags, ReadStructp rs, int direct){
  FragStore *myStore = FragStore_myStruct(fs);
  FragRecord *fr = (FragRecord *)rs;


    if(getIndexStore(myStore->fragStore, index, &fr->frag))
      return 1;
#ifdef DEBUG
  fprintf(stderr," getFragStore:  seqOffset = " F_U64 "  srcOffset = " F_U64 "\n",
	  fr->frag.sequenceOffset, fr->frag.sourceOffset);
#endif  
  unloadFragRecord(myStore, fr, getFlags);

#ifdef DEBUG
 fprintf(stderr,"* GetFragStore " F_S64 " with id " F_IID " src: %s  seq: %s \n qu: %s\n",
	 index, fr->frag.readIndex, fr->source, fr->sequence, fr->quality);
#endif
 return(0);
}


/***********************************************************************/
/* Public interfaces */
int getFragStore(FragStoreHandle fs, int64 index, int32 getFlags, ReadStructp rs){
  return getFragStoreCommon(fs,index,getFlags, rs, FALSE);
}


int getFragStoreDirect(FragStoreHandle fs, int64 index, int32 getFlags, ReadStructp rs){
  return getFragStoreCommon(fs,index,getFlags, rs, TRUE);
}

/***********************************************************************/
/* appendFragStorePartition
Append a ReadStruct to the store */
int appendFragStorePartition(FragStoreHandle store, ReadStructp rs, int32 partition){
  FragStore *myStore = FragStore_myStruct(store);
  FragRecord *fr = (FragRecord *)rs;
  static char encodeBuffer[MAX_SEQUENCE_LENGTH];
  StoreStat stats;
  uint16 sourceLength, screenMatchLength;
  VLSTRING_SIZE_T length;

  //  fprintf(stderr,"* appendFragStorePartition frag " F_IID " partition:%d\n", fr->frag.readIndex,partition);
  assert(myStore->numPartitions > 1);
  assert(partition >= 0 && partition < myStore->numPartitions);
  statsStore(myStore->sequenceStore[partition] , &stats);

  fr->frag.sequenceOffset = MERGE_FILEID_OFFSET(partition,stats.lastElem);
#if 0
  fprintf(stderr,"* partition %d offset " F_S64 "  extracted partition " F_S64 " extracted Offset " F_S64 "\n",
	  partition, stats.lastElem, 
	  GET_FILEID(fr->frag.sequenceOffset),
	  GET_FILEOFFSET(fr->frag.sequenceOffset));
#endif	  
  statsStore(myStore->sourceStore[partition] , &stats);
  fr->frag.sourceOffset = MERGE_FILEID_OFFSET(partition,stats.lastElem);
#if 0
  fprintf(stderr,"* partition %d offset " F_S64 "  extracted partition " F_S64 " extracted Offset " F_S64 "\n",
	  partition, stats.lastElem, 
	  GET_FILEID(fr->frag.sourceOffset),
	  GET_FILEOFFSET(fr->frag.sourceOffset));
#endif
  statsStore(myStore->fragStore , &stats);

#ifdef DEBUG_APPEND
  fprintf(stderr,"*** appendFragStore lastElem = " F_S64 " readIndex = " F_IID "\n",
	  stats.lastElem, fr->frag.readIndex);
#endif

  if(!((stats.lastElem + 1 == fr->frag.readIndex) ||
     ((stats.lastElem == stats.firstElem-1) && 
      (stats.firstElem == fr->frag.readIndex)))){
    fprintf(stderr,"***Oops: firstElem = " F_S64 " lastElem = " F_S64 "   fragIndex = " F_IID "\n",
	    stats.firstElem, stats.lastElem, fr->frag.readIndex);
    exit(1);
  }
  assert(strlen(fr->sequence) < MAX_SEQUENCE_LENGTH);
  encodeSequenceQuality(encodeBuffer, fr->sequence, fr->quality, fr->frag.hasQuality);


  appendIndexStore(myStore->fragStore, (void *)(&fr->frag));  // global index
  appendIndexStore(myStore->partitionStore[partition], (void *)(&fr->frag)); // per partition collection of fixed records

  screenMatchLength = fr->frag.numScreenMatches * sizeof(IntScreenMatch);
  sourceLength = (uint16)strlen(fr->source) + 1;

  /* Tack the screen matches onto the end of the source */
  if(fr->frag.numScreenMatches > 0){
    int32  offset = sourceLength + 1;
    memcpy(fr->source + offset, fr->matches, screenMatchLength);
  }
  /* Tack the locID and localPos onto the end of the source */
  if(!AS_FA_READ(fr->frag.readType)){ 
    int32  offset = sourceLength +  screenMatchLength + 1;
    int type = fr->frag.readType;

    memcpy(fr->source + offset, &fr->localeID, sizeof(uint64));
    offset += sizeof(uint64);

    if(AS_FA_SHREDDED(type)){ 
      memcpy(fr->source + offset, &fr->localePosStart, sizeof(uint32));
      offset += sizeof(uint32);
      memcpy(fr->source + offset, &fr->localePosEnd, sizeof(uint32));
      offset += sizeof(uint32);
    }
    length = offset;
  }else{
    length = sourceLength + screenMatchLength;
  }

#if 0
  fprintf(stderr,"* Appending source field of length " F_VLS " screenLength = %u\n",
	  length, screenMatchLength);
#endif
  assert(length <= VLSTRING_MAX_SIZE);
  appendVLRecordStore(myStore->sourceStore[partition]   , fr->source, length);

  /*** NOTE: encodeBuffer is NOT a null terminated string.  Therefore
   *** we use the elngth of the sequence data as the length of the VL
   *** record.
   ***/
  assert(strlen(fr->sequence) <= VLSTRING_MAX_SIZE);
  length = (VLSTRING_SIZE_T)strlen(fr->sequence);
  appendVLRecordStore(myStore->sequenceStore[partition]   , encodeBuffer, length);


  statsStore(myStore->fragStore , &stats);

  return(0);
}
/***********************************************************************/
/* appendFragStore
Append a ReadStruct to the store */
int appendFragStore(FragStoreHandle store, ReadStructp rs){
  FragStore *myStore = FragStore_myStruct(store);
  FragRecord *fr = (FragRecord *)rs;
  static char encodeBuffer[MAX_SEQUENCE_LENGTH];
  StoreStat stats;
  uint16 sourceLength, screenMatchLength;
  VLSTRING_SIZE_T length;

  assert(myStore->numPartitions == 1);
  statsStore(myStore->sequenceStore[0] , &stats);

  fr->frag.sequenceOffset = stats.lastElem;

  statsStore(myStore->sourceStore[0] , &stats);
  fr->frag.sourceOffset = stats.lastElem;

  statsStore(myStore->fragStore , &stats);

#ifdef DEBUG_APPEND
  fprintf(stderr,"*** appendFragStore lastElem = " F_S64 " readIndex = " F_IID "\n",
	  stats.lastElem, fr->frag.readIndex);
#endif

  if(!((stats.lastElem + 1 == fr->frag.readIndex) ||
     ((stats.lastElem == stats.firstElem-1) && 
      (stats.firstElem == fr->frag.readIndex)))){
    fprintf(stderr,"***Oops: firstElem = " F_S64 " lastElem = " F_S64 "   fragIndex = " F_IID "\n",
	    stats.firstElem, stats.lastElem, fr->frag.readIndex);
    exit(1);
  }
  assert(strlen(fr->sequence) < MAX_SEQUENCE_LENGTH);
  encodeSequenceQuality(encodeBuffer, fr->sequence, fr->quality, fr->frag.hasQuality);


  appendIndexStore(myStore->fragStore, (void *)(&fr->frag));

  screenMatchLength = fr->frag.numScreenMatches * sizeof(IntScreenMatch);
  sourceLength = (uint16)strlen(fr->source) + 1;

  /* Tack the screen matches onto the end of the source */
  if(fr->frag.numScreenMatches > 0){
    int32  offset = sourceLength + 1;
    memcpy(fr->source + offset, fr->matches, screenMatchLength);
  }
  /* Tack the locID and localPos onto the end of the source */
  if (!AS_FA_READ(fr->frag.readType)){ 
    int32  offset = sourceLength +  screenMatchLength + 1;
    int type = fr->frag.readType;

    memcpy(fr->source + offset, &fr->localeID, sizeof(uint64));
    offset += sizeof(uint64);

    if(AS_FA_SHREDDED(type)){ 
      memcpy(fr->source + offset, &fr->localePosStart, sizeof(uint32));
      offset += sizeof(uint32);
      memcpy(fr->source + offset, &fr->localePosEnd, sizeof(uint32));
      offset += sizeof(uint32);
    }
    length = offset;
  }else{
    length = sourceLength + screenMatchLength;
  }
  assert(length <= VLSTRING_MAX_SIZE);

#if 0
  fprintf(stderr,"* Appending source field of length " F_VLS " screenLength = %u\n",
	  length, screenMatchLength);
#endif
  appendVLRecordStore(myStore->sourceStore[0]   , fr->source, length);

  /*** NOTE: encodeBuffer is NOT a null terminated string.  Therefore
   *** we use the elngth of the sequence data as the length of the VL
   *** record.
   ***/
  assert(strlen(fr->sequence) <= VLSTRING_MAX_SIZE);
  length = (VLSTRING_SIZE_T)strlen(fr->sequence);
  appendVLRecordStore(myStore->sequenceStore[0]   , encodeBuffer, length);


  statsStore(myStore->fragStore , &stats);

  return(0);
}




/***********************************************************************/
/* appendFragStore
Append a ReadStruct to the store */
int appendDumpToFragStore(FragStoreHandle store, ShortFragRecord *fr, VA_TYPE(char) *sequence,  VA_TYPE(char) *source){
  FragStore *myStore = FragStore_myStruct(store);
  StoreStat stats;

  assert(myStore->numPartitions == 1);
  statsStore(myStore->sequenceStore[0] , &stats);

  fr->sequenceOffset = stats.lastElem;

  statsStore(myStore->sourceStore[0] , &stats);
  fr->sourceOffset = stats.lastElem;

  fprintf(stderr,"* sourceOffset = " F_U64 " sequenceOffset = " F_U64 "\n",
	  fr->sourceOffset, fr->sequenceOffset);

  statsStore(myStore->fragStore , &stats);

#ifdef DEBUG_APPEND
  fprintf(stderr,"*** appendFragStore lastElem = " F_S64 " readIndex = " F_IID "\n",
	  stats.lastElem, fr->readIndex);
#endif

  if(!((stats.lastElem + 1 == fr->readIndex) ||
     ((stats.lastElem == stats.firstElem-1) && 
      (stats.firstElem == fr->readIndex)))){
    fprintf(stderr,"***Oops: firstElem = " F_S64 " lastElem = " F_S64 "   fragIndex = " F_IID "\n",
	    stats.firstElem, stats.lastElem, fr->readIndex);
    exit(1);
  }
  assert(GetNumchars(sequence)  < MAX_SEQUENCE_LENGTH);

  appendIndexStore(myStore->fragStore, (void *)(fr));

  assert(GetNumchars(source) <= VLSTRING_MAX_SIZE);

  appendVLRecordStore(myStore->sourceStore[0]   , Getchar(source,0), GetNumchars(source));

  assert(GetNumchars(sequence) <= VLSTRING_MAX_SIZE);

  appendVLRecordStore(myStore->sequenceStore[0]   , Getchar(sequence,0), GetNumchars(sequence));

  statsStore(myStore->fragStore , &stats);

  return(0);
}

/***********************************************************************/
/* setFragStore
   You have to be careful not to modify the source/sequence offsets when
   you set.  If you modify 'em, you're on your own!
*/
int setFragStore(FragStoreHandle store, int64 index, ReadStructp rs){
  FragStore *myStore = FragStore_myStruct(store);
  FragRecord *fr = (FragRecord *)rs;

  setIndexStore(myStore->fragStore, index, (void *)(&fr->frag));

  return(0);
}

/*********************************************************************************/
int deleteFragStore(FragStoreHandle store, int64 index){
  FragRecord fr;
  FragStore *myStore = FragStore_myStruct(store);
  int64 lastElem;

  assert(myStore->status == ActiveStore);
  lastElem = getLastElemFragStore(store);

  assert(index >= 0 && index <= lastElem);
  
  // Read, modify and write the
  getIndexStore(myStore->fragStore, index, &fr);
  fr.frag.deleted = TRUE;
  setIndexStore(myStore->fragStore, index, (void *)(&fr));


  return(0);
}

/*********************************************************************************/
FragStreamHandle  openFragStream( FragStoreHandle fs, /* handle to a fragment store */
  void *buffer,    /* User supplied buffer for prefetching */
  int32 bufferSize){

  FragStream *myStream;

  myStream = allocateFragStream(fs, buffer, bufferSize);
  
  return FragStream_myHandle(myStream);

}

/*********************************************************************************/

int resetFragStream(FragStreamHandle fsh , int64 startIndex, int64 endIndex){
   int64 sourceOffset, sequenceOffset;

  FragStream *fs = FragStream_myStruct(fsh);
  assert(fs->status == ActiveStore);

#ifdef DEBUG
 fprintf(stderr,"*** resetFragStream %d " F_S64 "," F_S64 "\n",
	 fsh, startIndex, endIndex);
#endif

  if(startIndex == STREAM_FROMSTART){
     sourceOffset = STREAM_FROMSTART;
     sequenceOffset = STREAM_FROMSTART;
  }else{
     /* We want to start streaming from the middle of the store, so...
	we need to find the corresponding offsets in the auxiliary string/sequence
	stores.  This involves reading a frag from the index store, and retrieveing the 
	offsets */
     ReadStructp myRead;
     myRead =  new_ReadStruct();
     getFragStore(fs->fragStore, startIndex, FRAG_S_FIXED, myRead);
     getSourceOffset_ReadStruct(myRead, &sourceOffset);
     getSequenceOffset_ReadStruct(myRead, &sequenceOffset);
     delete_ReadStruct(myRead);
  }
  /* If we are not streaming from the beginning of the store, we have to
     get the startIndex-th element of the fragStore to get the start points
     for the other two streaming operations
  */
  resetStream(fs->fragStream, startIndex, endIndex);
  resetStream(fs->sequenceStream,sequenceOffset, STREAM_UNTILEND);
  resetStream(fs->sourceStream,sourceOffset, STREAM_UNTILEND);

  return(0);
}


/*********************************************************************************/

int  closeFragStream( FragStreamHandle fs){
  FragStream *myStream = FragStream_myStruct(fs);

#ifdef DEBUG
  fprintf(stderr,"*** closeFragStream\n");
#endif
  closeStream(myStream->fragStream);
  closeStream(myStream->sequenceStream);
  closeStream(myStream->sourceStream);

  myStream->status = UnAllocatedStore;
  gNumFragStreams--;

  return 0;

}
/* #define DEBUGNEXT */
/*********************************************************************************/
/* nextFragStream
    Streaming Operation.
*/
	
int nextFragStream(FragStreamHandle fh, ReadStructp rs, int streamFlags){
  FragStream *fs = FragStream_myStruct(fh);
  FragRecord *fr = (FragRecord *)rs;
  FragStore *myStore = FragStore_myStruct(fs->fragStore);
  int result;

  assert(fs->status == ActiveStore);

#ifdef DEBUGNEXT
  fprintf(stderr,"*** nextFragStream\n");
#endif

  result = 0;
  fr->frag.sequenceOffset = 0;
  fr->frag.sourceOffset = 0;
  if(nextStream(fs->fragStream, &fr->frag) == 0)
    return 0;
#ifdef DEBUGNEXT
  fprintf(stderr,"\tRead Frag, source @ " F_U64 "  sequence@ " F_U64 "\n",
	  fr->frag.sourceOffset, fr->frag.sequenceOffset);
#endif  
  unloadFragRecord(myStore,fr, streamFlags);

#ifdef DEBUG
  fprintf(stderr,"* nextFragStream Got Frag with id " F_IID " src: %s  seq: %s \n qu: %s\n",
	 fr->frag.readIndex, fr->source, fr->sequence, fr->quality);
#endif  
  return (1);

}

int64 getFirstElemFragStore(FragStoreHandle store){
  //  StoreStat stats;
  FragStore *myStore = FragStore_myStruct(store);

  return getFirstElemStore(myStore->fragStore);
}


int64 getLastElemFragStore(FragStoreHandle store){
  //  StoreStat stats;
  FragStore *myStore = FragStore_myStruct(store);

  return getLastElemStore(myStore->fragStore);
}


/***********************************************************************************/
int64 getStartIndexFragStream(FragStreamHandle stream){
  FragStream *myStream = FragStream_myStruct(stream);

  return getStartIndexStream(myStream->fragStream);;

}


/***** Private Accessors for FragStore ****/

int getSourceOffset_ReadStructp(ReadStructp rs, int64 *offset){
  FragRecord *fr = (FragRecord *)rs;
  *offset = fr->frag.sourceOffset;
  return(0);
}  
int getSequenceOffset_ReadStructp(ReadStructp rs, int64 *offset){
  FragRecord *fr = (FragRecord *)rs;
  *offset = fr->frag.sequenceOffset;
  return(0);
}





/***** Dump Frag Store Stats *****/
void DumpFragStoreStats(FILE *stream, FragStoreHandle handle){
  FragStore *store = FragStore_myStruct(handle);
  StoreStat stats;
  int64 firstFrag, lastFrag;
  int64 seqSize=0, srcSize=0;
  int i;

  statsStore(store->fragStore, &stats);
  firstFrag = stats.firstElem;  
  lastFrag = stats.lastElem;

  for(i = 0; i < store->numPartitions; i++){
    statsStore(store->sequenceStore[i], &stats);
    seqSize += stats.lastElem;

    statsStore(store->sourceStore[i], &stats);
    srcSize += stats.lastElem;
  }
  fprintf(stream," ** FragStore has %d partitions**\tfragments (" F_S64 "," F_S64 ") for a total of " F_S64 " (" F_S64 " chars)\n*\tsequence  " F_S64 " chars\n*\tsource " F_S64 " chars\n",
          store->numPartitions,
          firstFrag, lastFrag, lastFrag - firstFrag,
          (int64) ((lastFrag - firstFrag) * sizeof(ShortFragRecord)),
          seqSize, srcSize);

}
  
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
int statsFragStore(FragStoreHandle handle, StoreStat *stats){
  FragStore *store = FragStore_myStruct(handle);

  statsStore(store->fragStore, stats);
  return 0;
}

/***********************************************************************************
 * dumpDelimiter
 *
 *
 ***********************************************************************************/

#define FIELD_DELIM "~~~,~~~"

static char FieldDelim[] = FIELD_DELIM;

void dumpFieldDelimiter(FILE *outfp){
    safeWrite(outfp,FieldDelim, strlen(FieldDelim));  
}

void loadFieldDelimiter(FILE *outfp){
  char buffer[25];
    safeRead(outfp,buffer, strlen(FieldDelim));  
}

#define RECORD_DELIM "~~~.~~~"
static char RecordDelim[] = RECORD_DELIM;

void dumpRecordDelimiter(FILE *outfp){
    safeWrite(outfp,RecordDelim, strlen(RecordDelim));  
}

void loadRecordDelimiter(FILE *outfp){
  char buffer[25];
    safeRead(outfp,buffer, strlen(RecordDelim));  
}



/***********************************************************************************
 * Function unloadNDumpFragRecord
 *    Utility routine for unloading a FragRecord.  used by getNDumpFragStore
 ***********************************************************************************/
void unloadNDumpFragRecord(FragStore *myStore, FragRecord *fr, FILE *outfp){
  char encodeBuffer[MAX_SEQUENCE_LENGTH];
#ifdef DEBUG
  uint16 localeLength=0, screenLength;
#endif
  VLSTRING_SIZE_T actualLength;
  uint8 checksum;

  fprintf(outfp,"%u", fr->frag.deleted);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,"%u",fr->frag.readType);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,"%u",fr->frag.hasQuality);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,"%u",fr->frag.numScreenMatches);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,F_VLS, fr->frag.clearRegionStart);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,F_VLS, fr->frag.clearRegionEnd);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,F_UID,fr->frag.accID);
  dumpFieldDelimiter(outfp);
  fprintf(outfp,F_IID,fr->frag.readIndex);
  dumpFieldDelimiter(outfp);

    getVLRecordStore(myStore->sourceStore[0], 
		     fr->frag.sourceOffset, fr->source, 
		     (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, 
		     &actualLength);
    
#ifdef DEBUG
    fprintf(stderr,"* Read " F_VLS " chars %c\n",
            actualLength,
            fr->source[actualLength - localeLength - screenLength]);
#endif

    safeWrite(outfp,&actualLength, sizeof(actualLength));
    checksum = checkSumBlob(fr->source, actualLength);
    safeWrite(outfp,&checksum, sizeof(checksum));
    //    fprintf(stderr,"* source length = " F_VLS " checksum = %u\n",
    //    actualLength, checksum);
    safeWrite(outfp,fr->source, actualLength);


    getVLRecordStore(myStore->sequenceStore[0], 
		     fr->frag.sequenceOffset, encodeBuffer, (VLSTRING_SIZE_T)VLSTRING_MAX_SIZE, &actualLength);
    

    safeWrite(outfp,&actualLength, sizeof(actualLength));
    checksum = checkSumBlob(encodeBuffer, actualLength);
    safeWrite(outfp,&checksum, sizeof(checksum));
    //    fprintf(stderr,"* sequence length = " F_VLS " checksum = %u\n",
    //	    actualLength, checksum);
    safeWrite(outfp,encodeBuffer, actualLength);

  dumpRecordDelimiter(outfp);


}

void getVLDumpField(FILE *infp, VA_TYPE(char) *data, int fieldDelim){
  int state = 0;
  int ch;
  char c;
  ResetVA_char(data);
  
  //  while(ch = fgetc(infp), c = ch, ch != EOF){
  while(1){
    ch = fgetc(infp);
    c = ch;
    switch(state){
    case 0:
    case 1:
    case 2:
    case 4:
    case 5:
      if(c == '~')
	state++;
      else
	state = 0;
      break;

    case 6:
      if(c == '~'){
	//	fprintf(stderr,"* Found delimiter %d...returning\n", fieldDelim);
	Appendchar(data, &c);
	//		fprintf(stderr,"* Field is %d chars long including the delimiter\n",
	//			(int) GetNumchars(data));
	return;
      }
      else
	state = 0;
      break;

    case 3:
      if((fieldDelim && c == ',') ||
	 (!fieldDelim && c == '.'))
	state++;
      else
	state = 0;
      break;
    default:
      assert(0);
    }

    Appendchar(data, &c);



  }
  //  if(ch == EOF){
  //    fprintf(stderr,"* READ EOF\n");
  //  }
}

/***********************************************************************************
 * Function loadDumpFragRecord
 *    Utility routine for unloading a FragRecord.  used by getNDumpFragStore
 ***********************************************************************************/
int loadDumpFragRecord(FILE *infp, ShortFragRecord *fr, VA_TYPE(char ) *sequence, VA_TYPE(char) *source){
  int scratch;
  int result;
  uint8 recordChecksum, computedChecksum;
  VLSTRING_SIZE_T length;

  // char c;

  if(EOF == (result = fscanf(infp,"%d~~~,~~~", &scratch)))
    return FALSE;
  assert(1 == result);

  fr->deleted = scratch;
  if(1 != fscanf(infp,"%d~~~,~~~",&scratch)) assert(0);
  fr->readType = scratch;
  fscanf(infp,"%d~~~,~~~",&scratch);
  fr->hasQuality = scratch;
  fscanf(infp,"%d~~~,~~~",&scratch);
  fr->numScreenMatches = scratch;
  fscanf(infp,"%d~~~,~~~", &scratch);
  fr->clearRegionStart = scratch;
  if(1 != fscanf(infp,"%d~~~,~~~", &scratch)) assert(0);
  fr->clearRegionEnd = scratch;
  if(1 != fscanf(infp,F_UID "~~~,~~~",&fr->accID)) assert(0);
  if(1 != fscanf(infp,F_IID "~~~,~~~",&fr->readIndex)) assert(0);

  fprintf(stderr,"* loading frag with readIndex " F_IID "\n", fr->readIndex);

  safeRead(infp, &length, sizeof(length));
  safeRead(infp, &recordChecksum, sizeof(recordChecksum));

    fprintf(stderr,"* source length = " F_VLS " checksum = %u\n",
	    length, recordChecksum);
    ResetVA_char(source);
  {
    int i;
    for(i = 0; i < length; i++){
      char c = fgetc(infp);
      Appendchar(source, &c);
    }
  }
  fprintf(stderr,"* source = %s\n", Getchar(source,0));

  computedChecksum = checkSumBlob(Getchar(source,0), length);
  fprintf(stderr,"* recordChecksum %u computedChecksum %u\n",
	  recordChecksum, computedChecksum);
  assert(recordChecksum == computedChecksum);

  safeRead(infp, &length, sizeof(length));
  safeRead(infp, &recordChecksum, sizeof(recordChecksum));
    //  c = fgetc(infp);
    //  fprintf(stderr,"* Next char is %c\n", c);
    //  ungetc(c,infp);

    fprintf(stderr,"* sequence length = " F_VLS " checksum = %u\n",
	    length, recordChecksum);
    ResetVA_char(sequence);
  {
    int i;
    for(i = 0; i < length; i++){
      char c = fgetc(infp);
      Appendchar(sequence, &c);
    }
  }
  computedChecksum = checkSumBlob(Getchar(sequence,0), length);
  fprintf(stderr,"* recordChecksum %u computedChecksum %u\n",
	  recordChecksum, computedChecksum);
  assert(recordChecksum == computedChecksum);

  fscanf(infp,"%d~~~.~~~",&scratch);
  loadFieldDelimiter(infp);

  fflush(NULL);
  return TRUE;

}



/***********************************************************************/
/* getNDumpFragStore  
   Random Access Read.
   The data is read into the previously allocated ReadStruct rs.
   The types of data read can be restricted using the getFlags.
*/
	
int getNDumpFragStore(FragStoreHandle fs, int64 index, ReadStructp rs, FILE *outfp){
  FragStore *myStore = FragStore_myStruct(fs);
  FragRecord *fr = (FragRecord *)rs;

  getIndexStore(myStore->fragStore, index, &fr->frag);
#ifdef DEBUG
  fprintf(stderr," getFragStore:  seqOffset = " F_U64 "  srcOffset = " F_U64 "\n",
	  fr->frag.sequenceOffset, fr->frag.sourceOffset);
#endif  
  assert(fr->frag.readIndex == index);
  unloadNDumpFragRecord(myStore, fr, outfp);

#ifdef DEBUG
 fprintf(stderr,"* GetFragStore " F_S64 " with id " F_IID " src: %s  seq: %s \n qu: %s\n",
	 index, fr->frag.readIndex, fr->source, fr->sequence, fr->quality);
#endif
 return(0);
}



int kNextFragStream(FragStreamHandle fh, ReadStructp rs, int streamFlags,
		    int skipNum){
  FragStream *fs = FragStream_myStruct(fh);
  FragRecord *fr = (FragRecord *)rs;
  FragStore *myStore = FragStore_myStruct(fs->fragStore);
  int result;

  assert(fs->status == ActiveStore);

#ifdef DEBUGNEXT
  fprintf(stderr,"*** kNextFragStream\n");
#endif

  result = 0;
  
  fr->frag.sequenceOffset = 0;
  fr->frag.sourceOffset = 0;
  if(kNextStream(fs->fragStream, &fr->frag, skipNum) == 0)
    return 0;
#ifdef DEBUGNEXT
  fprintf(stderr,"\tRead Frag, source @ " F_U64 "  sequence@ " F_U64 "\n",
	  fr->frag.sourceOffset, fr->frag.sequenceOffset);
#endif  
  unloadFragRecord(myStore,fr, streamFlags);

#ifdef DEBUG
  fprintf(stderr,"* kNextFragStream Got Frag with id " F_IID " src: %s  seq: %s \n qu: %s\n",
	 fr->frag.readIndex, fr->source, fr->sequence, fr->quality);
#endif  
  return (1);

}
