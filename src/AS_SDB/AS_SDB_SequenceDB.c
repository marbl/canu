
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
static char CM_ID[] = "$Id: AS_SDB_SequenceDB.c,v 1.2 2004-09-23 20:25:28 mcschatz Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_SafeIO.h"
#include "AS_SDB_SequenceDB.h"

// #define DEBUG_OUTPUT

/* CreateSequenceDB:
   Create a SequenceDB with a single data file, and the index in memory..
   If force = true, this will remove an existing store with the same path.
   Fails if cannot open both the data and index files.
*/
tSequenceDB *CreateSequenceDB(char *path, int initialSize, int force){
  int32 size;
  DIR *dbDir;
  FILE *indexfp, *datafp;
  char buffer[FILENAME_MAX + 30];
  tSequenceDB *sequenceDB = (tSequenceDB *)malloc(sizeof(tSequenceDB));
  assert(path && sequenceDB);

  dbDir = opendir(path);
  if(dbDir == NULL){
    mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
      if(mkdir(path, mode)){
	sprintf(buffer,"createSequenceDB: Failure to create directory %s", path);
	perror(buffer);
	exit(1);
      }
  }else{
    closedir (dbDir);
  }

  sequenceDB->path = strdup(path);
  fprintf(stderr,"* Creating db with path %s\n", sequenceDB->path);
  sprintf(buffer,"%s/seqDB.unitigs.0", path);
  indexfp = fopen(buffer,"w");
  if(!indexfp){
	sprintf(buffer,"createSequenceDB: couldn't open unitig index file %s", buffer);
	perror(buffer);
	exit(1);
  }
  fclose(indexfp);

  sprintf(buffer,"%s/seqDB.contigs.0", path);
  indexfp = fopen(buffer,"w");
  if(!indexfp){
	sprintf(buffer,"createSequenceDB: couldn't open contig index file %s", buffer);
	perror(buffer);
	exit(1);
  }
  fclose(indexfp);

  sprintf(buffer,"%s/seqDB.data.0", path); 
  datafp = fopen(buffer,"w");  // Initial open is "w", subsequent are w+
  if(!datafp){
	sprintf(buffer,"createSequenceDB: couldn't open data file %s", buffer);
	perror(buffer);
	exit(1);
  }
 
  sequenceDB->SubStores = CreateVA_PtrT(10);
  AppendPtrT(sequenceDB->SubStores, (const void *)&datafp);
  size = max(2048,initialSize);
  sequenceDB->Unitigs = CreateVA_tMARecord(size);
  sequenceDB->Contigs = CreateVA_tMARecord(size);
  sequenceDB->UnitigStore = CreateMultiAlignStoreT(size);
  sequenceDB->ContigStore = CreateMultiAlignStoreT(size);
  sequenceDB->currentRevision = 0;
  sequenceDB->totalCacheSize = 0;
  sequenceDB->totalCacheLimit = 0;
  fprintf(stderr,"* Reseting offsetOfEOF\n");
  sequenceDB->offsetOfEOF = 0;
  sequenceDB->positionedAtEnd = 1;
  return sequenceDB;
}

void DeleteSequenceDB(tSequenceDB *db){
  int i;

  free(db->path);
  DeleteMultiAlignStoreT(db->UnitigStore);
  DeleteMultiAlignStoreT(db->ContigStore);
  DeleteVA_tMARecord(db->Unitigs);
  DeleteVA_tMARecord(db->Contigs);


  for(i = 0; i < GetNumPtrTs(db->SubStores) ; i++){
    FILE *storefp = (FILE *) *GetPtrT(db->SubStores, i);
    fclose(storefp);
  }
  DeleteVA_PtrT(db->SubStores);

  free(db);
}

  

/* OpenSequenceDB:
   Opens an existing SequenceDB.  The index file that will be opened corresponds
   to revision.  If this does not exist, failure.  If successful, the files
   corresponding to revisions 0,revision are opened read-only, and the files
   for revision+1 is opened r/w.
*/
tSequenceDB *OpenSequenceDB(char *path, int readWrite, int revision){
  tSequenceDB *sequenceDB = (tSequenceDB *)malloc(sizeof(tSequenceDB));
  char buffer[FILENAME_MAX + 30];
  FILE *datafp;
  FILE *testfp;
  int i;

  sequenceDB->path = strdup(path);

  // Load the indicies


  sprintf(buffer,"%s/seqDB.unitigs.%d",path,revision);
  testfp = fopen(buffer,"r");
  if(!testfp){
	char newBuffer[FILENAME_MAX + 30];
	fprintf( stderr, "FILENAME_MAX = %d\n", FILENAME_MAX);
	sprintf(newBuffer,"createSequenceDB: couldn't open unitig index file %s", buffer);
	perror(newBuffer);
	exit(1);
  }    

  sequenceDB->Unitigs = CreateFromFileVA_tMARecord(testfp,10);
  fprintf(stderr,"* Loaded " F_SIZE_T " tMARecords from %s\n",
          GetNumtMARecords(sequenceDB->Unitigs), buffer);
  fclose(testfp);

  sprintf(buffer,"%s/seqDB.contigs.%d",path,revision);
  testfp = fopen(buffer,"r");
  if(!testfp){
	sprintf(buffer,"createSequenceDB: couldn't open contig index file %s", buffer);
	perror(buffer);
	exit(1);
  }    

  sequenceDB->Contigs = CreateFromFileVA_tMARecord(testfp,10);
  fprintf(stderr,"* Loaded " F_SIZE_T " tMARecords from %s\n",
          GetNumtMARecords(sequenceDB->Contigs), buffer);
  fclose(testfp);


  sequenceDB->UnitigStore = CreateMultiAlignStoreT(GetNumtMARecords(sequenceDB->Unitigs));
  sequenceDB->ContigStore = CreateMultiAlignStoreT(GetNumtMARecords(sequenceDB->Contigs));

  sequenceDB->SubStores = CreateVA_PtrT(revision+1);

  for(i = 0; i <= revision; i++){

    sprintf(buffer,"%s/seqDB.data.%d",path,i);
    datafp = fopen(buffer,"r");
    AssertPtr(datafp);
    AppendPtrT(sequenceDB->SubStores, (const void *)&datafp);
  }

  sequenceDB->currentRevision = revision;
  fprintf(stderr,"* Reseting offsetOfEOF\n");
  sequenceDB->offsetOfEOF = 0;
  sequenceDB->positionedAtEnd = 1;
  if(readWrite){
    sprintf(buffer,"%s/seqDB.data.%d",path,revision + 1);
    datafp = fopen(buffer,"w+");
    AssertPtr(datafp);
    AppendPtrT(sequenceDB->SubStores, (const void *)&datafp);
    sequenceDB->currentRevision = revision+1;
  }
    return sequenceDB;
}

/* SaveSequenceDB
   Save the current revision of the indices.  Used by the checkpointing code.
   The indicies are maintained in memory.
*/
void SaveSequenceDB(tSequenceDB *db){  // Save the current revision of the indices
  char buffer[FILENAME_MAX + 30];
  FILE *currentDatafp, *indexfp;
  int64 end = db->offsetOfEOF;

  fprintf(stderr,"* SaveSequenceDB  end = " F_S64 "\n", end);
  sprintf(buffer,"%s/seqDB.unitigs.%d",db->path, db->currentRevision);
  indexfp = fopen(buffer,"w");
  CopyToFileVA_tMARecord(db->Unitigs, indexfp);
  fclose(indexfp);
  fprintf(stderr,"* Saved " F_SIZE_T " Unitig records to %s\n",
          GetNumtMARecords(db->Unitigs),buffer);

  sprintf(buffer,"%s/seqDB.contigs.%d",db->path, db->currentRevision);
  indexfp = fopen(buffer,"w");
  CopyToFileVA_tMARecord(db->Contigs, indexfp);
  fclose(indexfp);
  fprintf(stderr,"* Saved " F_SIZE_T " Contig records to %s\n",
          GetNumtMARecords(db->Contigs),buffer);

  /* Close the current data file, and reopen it as read only */
  currentDatafp = (FILE *) *GetPtrT(db->SubStores, db->currentRevision);
  fflush(currentDatafp);
  fclose(currentDatafp);

  sprintf(buffer,"%s/seqDB.data.%d",db->path,db->currentRevision);
  // Now we have to make sure the data has been written to the file before proceeding
  // if(end > 0)
  {
    while(1){
      int result;
      currentDatafp = fopen(buffer,"r");
      result = CDS_FSEEK(currentDatafp,
                         (off_t) (((end == 0) ? 0 : (end - (off_t)1))),SEEK_SET);
      if(result  == 0){ // couldn't seek until the end...still flushing
	break;
      }
      fprintf(stderr,"* Sleeping 15 seconds while buffers flush. Seeking to offset " F_S64 " returned %d\n",
	      end, result);
      perror("Error: ");
      CDS_FSEEK(currentDatafp, (off_t) 0, SEEK_END);
      fprintf(stderr,
              "* EOF is at offset " F_OFF_T " offsetOfOEF = " F_SIZE_T "\n",
              CDS_FTELL(currentDatafp), end);
      fflush(stderr);
      fclose(currentDatafp);
      sleep(15);
    }
  }
  /// SetPtrT(db->SubStores, db->currentRevision, (void *)&currentDatafp);
  SetPtrT(db->SubStores, db->currentRevision, (const void *) &currentDatafp);
  fprintf(stderr,"* Reopened %s as read only\n", buffer);

  db->currentRevision++;
  fprintf(stderr,"* Reseting offsetOfEOF\n");
  db->offsetOfEOF = 0;
  db->positionedAtEnd = 1;

  sprintf(buffer,"%s/seqDB.data.%d",db->path,db->currentRevision);
  currentDatafp = fopen(buffer,"w+");

  /// SetPtrT(db->SubStores, db->currentRevision, (void *)&currentDatafp);
  SetPtrT(db->SubStores, db->currentRevision, (const void *)&currentDatafp);

  fprintf(stderr,"* Opened %s as write\n", buffer);
}

void DuplicateEntryInSequenceDB(tSequenceDB *db, int fromIndex, int fromIsUnitig, int toIndex, int toIsUnitig, int keepInCache){
 tMARecord *fromMaRecord = GettMARecord(fromIsUnitig?db->Unitigs:db->Contigs, fromIndex);
 MultiAlignT *ma;
 MultiAlignStoreT *fromStore = (fromIsUnitig?db->UnitigStore:db->ContigStore);
 MultiAlignStoreT *toStore = (toIsUnitig?db->UnitigStore:db->ContigStore);
 SettMARecord(toIsUnitig?db->Unitigs:db->Contigs, toIndex, fromMaRecord);
 
 if(keepInCache){
   /* If this entry is in the from store's cache, add it to the to store's cache */
   ma = GetMultiAlignInStore(fromStore, fromIndex);
 }else{
   ma = NULL;
 }
 SetMultiAlignInStore(toStore,toIndex,ma);
 
}


/* InsertMultiAlignTInSequenceDB
   Inserts a new MultiAlignT into the appropriate cache and indices and
   appends it to the data file.
*/

void InsertMultiAlignTInSequenceDB(tSequenceDB *db,
                                   int index,
                                   int isUnitig,
                                   MultiAlignT *ma,
                                   int keepInCache){
 tMARecord *maRecord = GettMARecord(isUnitig?db->Unitigs:db->Contigs, index);
 MultiAlignStoreT *maStore = (isUnitig?db->UnitigStore:db->ContigStore);
 int32 fileID = GetNumPtrTs(db->SubStores) - 1;
 FILE *file = (FILE *)*GetPtrT(db->SubStores, fileID);
 size_t totalSize = 0;
 // off_t pos;

 assert((maRecord  == NULL )||
	(!maRecord->flags.bits.deleted &&
	 NULL == GetMultiAlignInStore(maStore,index)) );

 if(keepInCache){
   if(GetReferenceCountMultiAlignT(ma) == 0){ // we are the original owner
     db->totalCacheSize+= GetMemorySize(ma);
   }
  SetMultiAlignInStore(maStore,index,ma);
 }
 //  fflush(file);
 // fprintf(stderr,"* Appending at " F_SIZE_T "\n", db->offsetOfEOF);
 if(!db->positionedAtEnd)
   CDS_FSEEK(file, db->offsetOfEOF, SEEK_SET);
  //  pos = CDS_FTELL(file);

  if(maRecord){
    maRecord->storeID = fileID;
    maRecord->offset = db->offsetOfEOF;
  }else{
    tMARecord mar;
    mar.flags.all = 0;
    mar.storeID = fileID;
    mar.offset = db->offsetOfEOF;
    SettMARecord(isUnitig?db->Unitigs:db->Contigs, index, &mar);
  }
  //  fprintf(stderr,"* ma %d saved at offset " F_SIZE_T "\n", index, db->offsetOfEOF);
  totalSize += SaveMultiAlignTToStream(ma,file);
  db->offsetOfEOF += totalSize;
  db->positionedAtEnd = 1;
  //  fprintf(stderr,"* Appending " F_SIZE_T " bytes offset now " F_SIZE_T "\n", totalSize, db->offsetOfEOF);
}

/* DeleteMultiAlignTFromSequenceDB
   Mark the appropriate entry as deleted in an index.  This does not change the data file
*/
void DeleteMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig){
 tMARecord *maRecord = GettMARecord(isUnitig?db->Unitigs:db->Contigs, index);

 if(maRecord->flags.bits.deleted)
     return;

     maRecord->flags.bits.deleted = TRUE;
     UnloadMultiAlignTFromSequenceDB(db,index,isUnitig);
}

/* LoadMultiAlignTFromSequenceDB
   Loads a MultiAlignT into the cache (unless it is already present.
   On Failure, returns NULL.
*/
MultiAlignT *LoadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig){
 MultiAlignStoreT *maStore = (isUnitig?db->UnitigStore:db->ContigStore);
 MultiAlignT *ma = GetMultiAlignInStore(maStore,index);
 tMARecord *maRecord;
 int reference;
 FILE *file;
 if(ma){
#ifdef DEBUG_OUTPUT
   fprintf(stderr,"* Loaded %s ma %d from cache\n", (isUnitig?"Unitig":"Contig"),index);
#endif
   //   AddReferenceMultiAlignT(ma);
   return ma;
 }

 maRecord = GettMARecord(isUnitig?db->Unitigs:db->Contigs, index);

 // fprintf(stderr,"* Loading from file %d offset " F_U64 "\n", maRecord->storeID, maRecord->offset);
 // fflush(stderr);

 if(maRecord->flags.bits.deleted)
   return NULL;

 file = (FILE*) *GetPtrT(db->SubStores, maRecord->storeID);

 db->positionedAtEnd = 0;
 CDS_FSEEK(file,maRecord->offset,SEEK_SET);
 ma = LoadMultiAlignTFromStream(file,&reference);
 // AddReferenceMultiAlignT(ma);

 if(ma == NULL){
   fprintf(stderr,"** LoadMultiAlignTFromSequenceDB FAILED for %s %d in file %d at offset " F_OFF_T " (reference:%d)\n",
	   (isUnitig?"Unitig":"Contig"), index, maRecord->storeID, maRecord->offset, reference);
   AssertPtr(ma);
 }

 db->totalCacheSize+= GetMemorySize(ma);
#ifdef DEBUG_OUTPUT
 fprintf(stderr,
         "* %s ma %d loaded from offset " F_OFF_T " cache = " F_SIZE_T "\n",
         (isUnitig?"Unitig":"Contig"), index,
         CDS_FTELL(file), db->totalCacheSize);
#endif
  SetMultiAlignInStore(maStore,index,ma);
 
 return ma;
}

/* ReLoadMultiAlignTFromSequenceDB
   If the multiAlign is in cache, copy it into the reusema
   Otherwise, load reusema from disk, and DO NOT insert in cache
   On Failure, returns Empty reusema.
*/
//#define DEBUG_OUTPUT

int32 ReLoadMultiAlignTFromSequenceDB(tSequenceDB *db, MultiAlignT *reusema, int index, int isUnitig){
 MultiAlignStoreT *maStore = (isUnitig?db->UnitigStore:db->ContigStore);
 MultiAlignT *ma = GetMultiAlignInStore(maStore,index);
 tMARecord *maRecord;
 int reference;
 FILE *file;
 if(ma){
#ifdef DEBUG_OUTPUT
   fprintf(stderr,"* Loaded %s ma %d from cache\n",
           (isUnitig?"Unitig":"Contig"),index);
#endif

   CopyMultiAlignT(reusema, ma);
   return 0;
 }

 maRecord = GettMARecord(isUnitig?db->Unitigs:db->Contigs, index);

 if(maRecord->flags.bits.deleted){
   //   fprintf(stderr,"* ReLoadMultiAlignTFromSequenceDB -- trying to load deleted %s %d\n",
   //	   (isUnitig?"Unitig":"Contig"), index);
   return 1;
 }
 file = (FILE*) *GetPtrT(db->SubStores, maRecord->storeID);

 db->positionedAtEnd = 0;
#ifdef DEBUG_OUTPUT
 fprintf(stderr,"* Loading from file %d at offset " F_U64 "\n",
	 maRecord->storeID, maRecord->offset);
#endif
 CDS_FSEEK(file,maRecord->offset,SEEK_SET);
 ReLoadMultiAlignTFromStream(file,reusema,&reference);
 return 0;
}


/* UnLoadMultiAlignTFromSequenceDB
   Unloads a MultiAlignT from a cache.
*/
void UnloadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig){
 MultiAlignStoreT *maStore = (isUnitig?db->UnitigStore:db->ContigStore);
 size_t redeemed;
 redeemed = RemoveMultiAlignFromStore(maStore,index);
 db->totalCacheSize -= redeemed;
#ifdef DEBUG_OUTPUT
 fprintf(stderr,
         "* Unloading %s %d redeemed " F_SIZE_T " bytes cache = " F_SIZE_T "\n",
         (isUnitig?"Unitig":"Contig"), index, redeemed, db->totalCacheSize);
#endif
}

/* ClearCacheSequenceDB
   Clears all the entries from a cache.
*/
void ClearCacheSequenceDB(tSequenceDB *db, int isUnitig){
  MultiAlignStoreT *maStore = (isUnitig?db->UnitigStore:db->ContigStore);
  size_t redeemed = ClearMultiAlignStoreT(maStore);
  db->totalCacheSize -=  redeemed;
  fprintf(stderr,"* Clear %s Cache redeemed = " F_SIZE_T "\n",
          (isUnitig?"Unitig":"Contig"), redeemed);
}

void ClearCacheSequenceDBConditionally(tSequenceDB *db, size_t maxSize){
  
  if(db->totalCacheSize > maxSize){
    ClearCacheSequenceDB(db,TRUE);
    ClearCacheSequenceDB(db,FALSE);
    db->totalCacheSize = 0;
  }

}

/* Return the size of the current data file being written for the sequence DB */
size_t GetSizeOfCurrentSequenceDB(tSequenceDB *db){
  return(db->offsetOfEOF);
}

