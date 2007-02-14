
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
static char CM_ID[] = "$Id: AS_SDB_SequenceDBPartition.c,v 1.5 2007-02-14 07:20:13 brianwalenz Exp $";

/*************************************************************************
 Module:  AS_SDB_SequenceDBPartition
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
#include "AS_SDB_SequenceDBPartition.h"


// The partition files for a given revision r of an SDB s are of the form:
// s/db.r.n

int testOpenSequenceDBPartition(char *sequenceDBPath, int32 revision, int32 partition){
  char buffer[FILENAME_MAX];
  FILE *seqfp, *srcfp;
  int32 count = 0;



  /* Open the specific partition we want loaded in its entirety */
  sprintf(buffer,"%s/seqDB.%d.%d", sequenceDBPath,revision, partition);
  if(!(srcfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openSequenceDBPartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(srcfp);
    count++;
  }
  sprintf(buffer,"%s/seqDB.data.%d.%d", sequenceDBPath,revision,partition);
  if(!(seqfp = fopen(buffer,"r"))){
    fprintf(stderr,"* openSequenceDBPartition: Couldn't open %s...\n",buffer);
  }else{
    fclose(seqfp);
    count++;
  }

  switch(count){
  case 2:
    return 1;

  case 0:
    return 0;

  default:
    //return -1;
    break;
  }
  return -1;
}


static void buildSequenceDBPartitionHash(tSequenceDBPartition *partition){
  int i;
  int numMultiAligns;

  numMultiAligns = GetNumtMARecords(partition->multiAligns);
  partition->index = CreateHashTable_int32_AS(numMultiAligns + 1);

  for(i = 0; i < numMultiAligns; i++){
    tMARecord *maRecord = GettMARecord(partition->multiAligns, i);
    MultiAlignT *ma;
    int reference;
    //fprintf(stderr,"* maRecord storeID %d  offset " F_U64 "\n",
    //	    maRecord->storeID, maRecord->offset);
    CDS_FSEEK(partition->datafp,maRecord->offset,SEEK_SET);
    ma = LoadMultiAlignTFromStream(partition->datafp,&reference);
    if(InsertInHashTable_AS(partition->index,
                            (void *)(&maRecord->storeID),
                            sizeof(int32), ma) != HASH_SUCCESS)
      assert(0);
  }
}


tSequenceDBPartition *openSequenceDBPartition(char *sequenceDBPath, int32 revision, int32 partition){
  FILE *testfp;
  char buffer[FILENAME_MAX];
  tSequenceDBPartition *SequenceDBPartition= NULL;

  SequenceDBPartition = (tSequenceDBPartition *)safe_malloc(sizeof(tSequenceDBPartition));
  SequenceDBPartition->sequenceDBPath = strdup(sequenceDBPath);

  if(1 != testOpenSequenceDBPartition(sequenceDBPath, revision, partition)){
    return NULL;
  }
  
  sprintf(buffer,"%s/seqDB.%d.%d", sequenceDBPath,revision,partition);
  testfp = fopen(buffer,"r");
  if(!testfp){
	sprintf(buffer,"openSequenceDBPartition: couldn't open index file %s", buffer);
	perror(buffer);
	exit(1);
  }    

  SequenceDBPartition->multiAligns = CreateFromFileVA_tMARecord(testfp,10);
  fclose(testfp);

  sprintf(buffer,"%s/seqDB.data.%d.%d", sequenceDBPath,revision,partition);
  SequenceDBPartition->datafp = fopen(buffer,"r");

  buildSequenceDBPartitionHash(SequenceDBPartition);

  return SequenceDBPartition;

}


void closeSequenceDBPartition(tSequenceDBPartition *partition){
  fclose(partition->datafp);
  Delete_VA(partition->multiAligns);
  safe_free(partition->sequenceDBPath);
  safe_free(partition);
}

int isMemberSequenceDBPartition(tSequenceDBPartition *partition, int32 indx){
  MultiAlignT *ma = (MultiAlignT *)LookupInHashTable_AS(partition->index, (void *)&indx , sizeof(int32));

  return (ma != NULL);

}

MultiAlignT *loadFromSequenceDBPartition(tSequenceDBPartition *partition, int32 indx){
  MultiAlignT *ma = (MultiAlignT *)LookupInHashTable_AS(partition->index, (void *)&indx , sizeof(int32));

  if(!ma){ // not in our partition, open the relevant stores, get the data, and close them
    fprintf(stderr,"* Multialign %d is NOT in this partition!   Exiting\n", indx);
    exit(1);
  }
  return ma;
}


VA_TYPE(int32) *GetContentSequenceDBPartition(tSequenceDBPartition *partition){
  HashTable_Iterator_AS iterator;
  VA_TYPE(int32) *members = CreateVA_int32(100);
  void *keyp, *valuep;
  InitializeHashTable_Iterator_AS(partition->index, &iterator);

  fprintf(stderr,"* SDBPartition has following members\n");
  while(NextHashTable_Iterator_AS(&iterator, &keyp, &valuep)){
    int32 key = *(int32 *)keyp;
    ///Appendint32(members, keyp);  /// C++PROJECT Next line makes this safer
    Appendint32(members, &key);
    fprintf(stderr,"* %d\n", key);
  }
  return members;
}
