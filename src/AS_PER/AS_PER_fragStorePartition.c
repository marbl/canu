
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
static char CM_ID[] = "$Id: AS_PER_fragStorePartition.c,v 1.4 2005-03-22 19:49:20 jason_miller Exp $";

/*************************************************************************
 Module:  AS_PER_fragStorePartition
 Description:
   Access a partition of a partitioned frag store

 Assumptions:

 Document:

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
#include "AS_PER_fragStorePartition.h"


int testOpenFragStorePartition(char *fragStorePath, int32 partition){
  char buffer[FILENAME_MAX];
  FILE *seqfp, *srcfp, *frgfp;
  int32 count = 0;


  /* Open the global files so we can handle misses */
  sprintf(buffer,"%s/db.frg", fragStorePath);
  if(!(frgfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openFragStorePartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(frgfp);
    count++;
  }

  /* Open the specific partition we want loaded in its entirety */
  sprintf(buffer,"%s/db.src.%d", fragStorePath,partition);
  if(!(srcfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openFragStorePartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(srcfp);
    count++;
  }
  sprintf(buffer,"%s/db.seq.%d", fragStorePath,partition);
  if(!(seqfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openFragStorePartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(seqfp);
    count++;
  }
  sprintf(buffer,"%s/db.frg.%d", fragStorePath,partition);
  if(!(frgfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openFragStorePartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(frgfp);
    count++;
  }

  switch(count){
  case 4:
    return 1;

  case 0:
    return 0;

  default:
    //return -1;
    break;
  }
  return -1;
}


static void buildFragStorePartitionHash(tFragStorePartition *partition){
  StoreStat stats;
  int i;

  statsStore(partition->partitionStore, &stats);
  partition->index = CreateHashTable_int32_AS(stats.lastElem + 1);

  for(i = stats.firstElem; i <= stats.lastElem; i++){
    void *elemPtr;
    getIndexStorePtr(partition->partitionStore, i, &elemPtr);
    if(InsertInHashTable_AS(partition->index,
                            (void *)(&((ShortFragRecord *)elemPtr)->readIndex),
                            sizeof(int32), elemPtr) != HASH_SUCCESS)
      assert(0);
  }
}


tFragStorePartition *openFragStorePartition(char *fragStorePath, int32 partition, int loadData){
  char buffer[FILENAME_MAX];
  tFragStorePartition *fragStorePartition= NULL;

  fragStorePartition = (tFragStorePartition *)malloc(sizeof(tFragStorePartition));
  fragStorePartition->fragStorePath = strdup(fragStorePath);

  if(1 != testOpenFragStorePartition(fragStorePath, partition)){
    return NULL;
  }
  
  sprintf(buffer,"%s/db.frg.%d", fragStorePath,partition);
  fragStorePartition->partitionStore = loadStore(buffer);

  sprintf(buffer,"%s/db.seq.%d", fragStorePath,partition);
  if(!loadData){
    fragStorePartition->sequenceStore = openStore(buffer,"r");
  }else{    
    fragStorePartition->sequenceStore = loadStore(buffer);
  }
  sprintf(buffer,"%s/db.src.%d", fragStorePath,partition);
  if(!loadData){
    fragStorePartition->sourceStore = openStore(buffer,"r");
  }else{    
    fragStorePartition->sourceStore = loadStore(buffer);
  }

  sprintf(buffer,"%s/db.frg", fragStorePath);
  fragStorePartition->globalFrgStore = openStore(buffer,"r");

  buildFragStorePartitionHash(fragStorePartition);

  return fragStorePartition;

}


void closeFragStorePartition(tFragStorePartition *partition){
  closeStore(partition->globalFrgStore);
  closeStore(partition->sourceStore);
  closeStore(partition->sequenceStore);
  closeStore(partition->partitionStore);
  free(partition->fragStorePath);
  free(partition);
}

int isMemberFragStorePartition(tFragStorePartition *partition, int32 indx){
  ShortFragRecord *sfr = (ShortFragRecord *)
    LookupInHashTable_AS(partition->index, (void *)&indx , sizeof(int32));

  return (sfr != NULL);

}

int getFragStorePartition(tFragStorePartition *partition, int32 indx, int32 getFlags, ReadStructp rs){
  ShortFragRecord *sfr = (ShortFragRecord *)
    LookupInHashTable_AS(partition->index, (void *)&indx , sizeof(int32));
  FragRecord *fr = (FragRecord * )rs;
  AssertPtr(rs);

  if(!sfr){ // not in our partition, open the relevant stores, get the data, and close them
    char buffer[FILENAME_MAX];
    StoreHandle fSourceStore;
    StoreHandle fSequenceStore;
    int partition_to_open;
    
    getIndexStore(partition->globalFrgStore, indx, &fr->frag);
    partition_to_open = GET_FILEID(fr->frag.sourceOffset);
    sprintf(buffer,"%s/db.src.%d",partition->fragStorePath, partition_to_open);
    fSourceStore = openStore(buffer,"r");
    sprintf(buffer,"%s/db.seq.%d",partition->fragStorePath, partition_to_open);
    fSequenceStore = openStore(buffer,"r");
    unloadFragRecordPartition(fSequenceStore, fSourceStore, fr, getFlags);
    closeStore(fSourceStore);
    closeStore(fSequenceStore);
  }else{
    fr->frag = *sfr;
    unloadFragRecordPartition(partition->sequenceStore, partition->sourceStore, fr, getFlags);
  }

  return (0);

}
