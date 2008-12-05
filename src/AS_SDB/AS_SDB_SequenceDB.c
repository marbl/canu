
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
static char *rcsid = "$Id: AS_SDB_SequenceDB.c,v 1.21 2008-12-05 19:06:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_SDB_SequenceDB.h"


tSequenceDB *
createSequenceDB(char *path) {
  FILE        *F;
  char         N[FILENAME_MAX];

  tSequenceDB *db = (tSequenceDB *)safe_calloc(1, sizeof(tSequenceDB));

  db->path              = strdup(path);
  db->currentRevision   = 0;
  db->Unitigs           = CreateVA_tMARecord(2048);
  db->Contigs           = CreateVA_tMARecord(2048);
  db->UnitigStore       = CreateMultiAlignStoreT();
  db->ContigStore       = CreateMultiAlignStoreT();
  db->dataFileLen       = 0;
  db->dataFileMax       = 1024;
  db->dataFile          = (FILE **)safe_calloc(db->dataFileMax, sizeof(FILE *));

  AS_UTL_mkdir(db->path);

  sprintf(N,"%s/seqDB.v%03d.utg", db->path, 0);
  errno = 0;
  F = fopen(N,"w");
  if (errno)
    fprintf(stderr, "CreateSequenceDB()-- Failed to create '%s': %s\n",
            N, strerror(errno)), exit(1);
  fclose(F);

  sprintf(N,"%s/seqDB.v%03d.ctg", db->path, 0);
  errno = 0;
  F = fopen(N,"w");
  if (errno)
    fprintf(stderr, "CreateSequenceDB()-- Failed to create '%s': %s\n",
            N, strerror(errno)), exit(1);
  fclose(F);

  sprintf(N,"%s/seqDB.v%03d.dat", db->path, 0);
  errno = 0;
  F = fopen(N,"w+");
  if (errno)
    fprintf(stderr, "CreateSequenceDB()-- Failed to create '%s': %s\n",
            N, strerror(errno)), exit(1);

  db->dataFile[db->dataFileLen++] = F;

  return(db);
}


tSequenceDB *
openSequenceDB(char *path, int readWrite, int revision){
  FILE *F;
  char N[FILENAME_MAX + 30];
  FILE *datafp;
  tSequenceDB *db = (tSequenceDB *)safe_calloc(1, sizeof(tSequenceDB));
  int i;

  db->path              = strdup(path);
  db->currentRevision   = revision;
  db->Unitigs           = NULL;
  db->Contigs           = NULL;
  db->UnitigStore       = CreateMultiAlignStoreT();
  db->ContigStore       = CreateMultiAlignStoreT();
  db->dataFileLen       = 0;
  db->dataFileMax       = 1024;
  db->dataFile          = (FILE **)safe_calloc(db->dataFileMax, sizeof(FILE *));

  sprintf(N,"%s/seqDB.v%03d.utg", db->path, revision);
  errno = 0;
  F = fopen(N,"r");
  if (errno)
    fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
            N, strerror(errno)), exit(1);
  db->Unitigs = CreateFromFileVA_tMARecord(F);
  fclose(F);

  sprintf(N,"%s/seqDB.v%03d.ctg", db->path, revision);
  errno = 0;
  F = fopen(N,"r");
  if (errno)
    fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
            N, strerror(errno)), exit(1);
  db->Contigs = CreateFromFileVA_tMARecord(F);
  fclose(F);


  for (i=0; i<=revision; i++) {
    sprintf(N, "%s/seqDB.v%03d.dat", db->path, i);
    errno = 0;
    F = fopen(N, "r");
    if (errno)
      fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
              N, strerror(errno)), exit(1);
    db->dataFile[db->dataFileLen++] = F;
  }

  if(readWrite){
    db->currentRevision = revision + 1;
    sprintf(N, "%s/seqDB.v%03d.dat", db->path, db->currentRevision);
    errno = 0;
    F = fopen(N, "w+");
    if (errno)
      fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for write+append (is it read only?): %s\n",
              N, strerror(errno)), exit(1);
    db->dataFile[db->dataFileLen++] = F;
  }

  return db;
}



void
openSequenceDBPartition(tSequenceDB *db, int32 partition){
  FILE *F = NULL;
  char  N[FILENAME_MAX];
  int   i;

  sprintf(N, "%s/seqDB.v%03d.dat.i%03d", db->path, db->currentRevision, partition);

  errno = 0;
  F = fopen(N, "r");
  if (errno) {
    fprintf(stderr, "openSequenceDBPartition()-- couldn't open '%s': %s\n", N, strerror(errno));
    exit(1);
  }

  db->multiAligns      = CreateFromFileVA_tMARecord(F);
  db->multiAlignLookup = CreateScalarHashTable_AS();

  fclose(F);

  //  Build a hash of multiAlignID to multiAlignT pointer.

  sprintf(N,"%s/seqDB.v%03d.dat.p%03d", db->path, db->currentRevision, partition);

  errno = 0;
  F = fopen(N,"r");
  if (errno) {
    fprintf(stderr, "openSequenceDBPartition()-- couldn't open '%s': %s\n", N, strerror(errno));
    exit(1);
  }

  for (i=0; i<GetNumtMARecords(db->multiAligns); i++){
    tMARecord     *maRecord = GettMARecord(db->multiAligns, i);
    MultiAlignT   *ma       = NULL;

    AS_UTL_fseek(F, maRecord->offset, SEEK_SET);
    ma = LoadMultiAlignTFromStream(F);

    if(InsertInHashTable_AS(db->multiAlignLookup,
                            (uint64)maRecord->multiAlignID, 0,
                            (uint64)(INTPTR)(ma), 0) != HASH_SUCCESS) {
      fprintf(stderr, "Failed to insert multiAlign %d into the lookup table.  Already there?\n", maRecord->multiAlignID);
      assert(0);
    }
  }

  fclose(F);
}



// Save the current revision of the indices
void
saveSequenceDB(tSequenceDB *db) {
  FILE   *F = NULL;
  char    N[FILENAME_MAX + 30];

  sprintf(N,"%s/seqDB.v%03d.utg",db->path, db->currentRevision);
  errno = 0;
  F = fopen(N,"w");
  if (errno)
    fprintf(stderr, "SaveSequenceDB()-- Failed to open '%s' for write: %s\n", N, strerror(errno)), exit(1);
  CopyToFileVA_tMARecord(db->Unitigs, F);
  fclose(F);

  sprintf(N,"%s/seqDB.v%03d.ctg",db->path, db->currentRevision);
  errno = 0;
  F = fopen(N,"w");
  if (errno)
    fprintf(stderr, "SaveSequenceDB()-- Failed to open '%s' for write: %s\n", N, strerror(errno)), exit(1);
  CopyToFileVA_tMARecord(db->Contigs, F);
  fclose(F);

  // Close the current data file, and reopen it as read only

  errno = 0;
  fclose(db->dataFile[db->currentRevision]);
  if (errno)
    fprintf(stderr, "SaveSequenceDB()-- Failed to close '%s': %s\n", N, strerror(errno)), exit(1);

  sprintf(N,"%s/seqDB.v%03d.dat",db->path,db->currentRevision);

  errno = 0;
  F = fopen(N, "r");
  if (errno)
    fprintf(stderr, "SaveSequenceDB()-- Failed to open '%s' for read: %s\n", N, strerror(errno)), exit(1);
  db->dataFile[db->currentRevision] = F;

  db->currentRevision++;

  sprintf(N,"%s/seqDB.v%03d.dat", db->path, db->currentRevision);
  errno = 0;
  F = fopen(N,"w+");
  if (errno)
    fprintf(stderr, "SaveSequenceDB()-- Failed to open '%s' for write: %s\n", N, strerror(errno)), exit(1);
  db->dataFile[db->currentRevision] = F;
}


void
deleteSequenceDB(tSequenceDB *db){
  int i;

  DeleteVA_tMARecord(db->Unitigs);
  DeleteVA_tMARecord(db->Contigs);
  DeleteMultiAlignStoreT(db->UnitigStore);
  DeleteMultiAlignStoreT(db->ContigStore);

  for (i=0; i<db->dataFileLen; i++)
    fclose(db->dataFile[i]);

  safe_free(db->path);
  safe_free(db->dataFile);
  safe_free(db);
}









#define MAStore(isUnitig) ((isUnitig) ? db->Unitigs : db->Contigs)




void
updateMultiAlignTInSequenceDB(tSequenceDB *db,
                              int index,
                              int isUnitig,
                              MultiAlignT *ma,
                              int keepInCache){
  MultiAlignStoreT *maStore  = (isUnitig) ? db->UnitigStore : db->ContigStore;
  tMARecord        *maRecord = GettMARecord(MAStore(isUnitig), index);

  //fprintf(stderr, "updateMultiAlignTFromSequenceDB()--  ma 0x%016p index=%d utg=%d\n", ma, index, isUnitig);

  assert(maRecord != NULL);

  if (keepInCache)
    SetMultiAlignInStore(maStore, index, ma);

  AS_UTL_fseek(db->dataFile[db->currentRevision], 0, SEEK_END);

  maRecord->storeID      = db->currentRevision;
  maRecord->multiAlignID = -1;
  maRecord->isDeleted    = 0;
  maRecord->offset       = AS_UTL_ftell(db->dataFile[db->currentRevision]);

  SaveMultiAlignTToStream(ma, db->dataFile[db->currentRevision]);
}


void
insertMultiAlignTInSequenceDB(tSequenceDB *db,
                              int index,
                              int isUnitig,
                              MultiAlignT *ma,
                              int keepInCache){

  MultiAlignStoreT *maStore  = (isUnitig) ? db->UnitigStore : db->ContigStore;
  tMARecord        *maRecord = GettMARecord(MAStore(isUnitig), index);

  //fprintf(stderr, "insertMultiAlignTFromSequenceDB()--  ma 0x%016p index=%d utg=%d\n", ma, index, isUnitig);

  //  We can either have:
  //
  //    no maRecord -- it is not in the store
  //
  //    with maRecord -- it is in the store -- and not deleted -- and with no
  //    existing multialign stored.
  assert((maRecord  == NULL )||
         (!maRecord->isDeleted && (NULL == GetMultiAlignInStore(maStore,index))));

  if (maRecord == NULL) {
    tMARecord mar = {0};
    SettMARecord(MAStore(isUnitig), index, &mar);
  }

  updateMultiAlignTInSequenceDB(db, index, isUnitig, ma, keepInCache);
}


void
deleteMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig){
  tMARecord        *maRecord = GettMARecord(MAStore(isUnitig), index);
  if ((maRecord == NULL) || (maRecord->isDeleted))
    return;
  //fprintf(stderr, "deleteMultiAlignTFromSequenceDB()--  index=%d utg=%d\n", index, isUnitig);
  maRecord->isDeleted = TRUE;
  RemoveMultiAlignFromStore(isUnitig ? db->UnitigStore : db->ContigStore, index);
}


MultiAlignT *
loadMultiAlignTFromSequenceDB(tSequenceDB *db, int index, int isUnitig){
  MultiAlignStoreT *maStore  = (isUnitig) ? db->UnitigStore : db->ContigStore;
  MultiAlignT      *ma       = NULL;
  tMARecord        *maRecord = NULL;

  //  If it's in the partition, return that -- only valid for contig?
  if (db->multiAlignLookup)
    ma = (MultiAlignT *)(INTPTR)LookupValueInHashTable_AS(db->multiAlignLookup, index, 0);
  if (ma)
    return(ma);

  //  Otherwise, grab it from the store
  ma = GetMultiAlignInStore(maStore, index);
  if (ma)
    return(ma);

  maRecord = GettMARecord(MAStore(isUnitig), index);

  if (maRecord == NULL)
    fprintf(stderr, "loadMultiAlignTFromSequenceDB()-- Unable to extract MA Record with iid #%d\n", index);
  assert(maRecord != NULL);

  if (maRecord->isDeleted)
    return(NULL);

  AS_UTL_fseek(db->dataFile[maRecord->storeID], maRecord->offset, SEEK_SET);
  ma = LoadMultiAlignTFromStream(db->dataFile[maRecord->storeID]);

  //fprintf(stderr, "loadMultiAlignTFromSequenceDB()--  ma 0x%016p index=%d utg=%d\n", ma, index, isUnitig);

  if (ma == NULL)
    fprintf(stderr,"loadMultiAlignTFromSequenceDB()-- FAILED for %s %d in file %d at offset "F_OFF_T"\n",
            (isUnitig?"Unitig":"Contig"), index, maRecord->storeID, maRecord->offset);
  assert(ma != NULL);

  SetMultiAlignInStore(maStore, index, ma);

  return(ma);
}



void
copyMultiAlignTFromSequenceDB(tSequenceDB *db, MultiAlignT *reusema, int index, int isUnitig) {
  MultiAlignStoreT *maStore  = (isUnitig?db->UnitigStore:db->ContigStore);
  MultiAlignT      *ma       = NULL;
  tMARecord        *maRecord = NULL;

  //  If it's in the partition, return that -- only valid for contig?
  //  Otherwise, grab it from the store
  if (db->multiAlignLookup)
    ma = (MultiAlignT *)(INTPTR)LookupValueInHashTable_AS(db->multiAlignLookup, index, 0);
  if (!ma)
    ma = GetMultiAlignInStore(maStore,index);
  if (ma) {
    CopyMultiAlignT(reusema, ma);
    return;
  }

  maRecord = GettMARecord(MAStore(isUnitig), index);

  if (maRecord->isDeleted) {
    ClearMultiAlignT(reusema);
    return;
  }

  AS_UTL_fseek(db->dataFile[maRecord->storeID], maRecord->offset, SEEK_SET);
  ReLoadMultiAlignTFromStream(db->dataFile[maRecord->storeID], reusema);
}
