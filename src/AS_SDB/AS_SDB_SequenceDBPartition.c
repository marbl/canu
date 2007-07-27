
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

static char CM_ID[] = "$Id: AS_SDB_SequenceDBPartition.c,v 1.10 2007-07-27 14:16:07 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_SDB_SequenceDBPartition.h"

tSequenceDBPartition *
openSequenceDBPartition(char *sequenceDBPath, int32 revision, int32 partition){
  FILE *testfp;
  char buffer[FILENAME_MAX];
  tSequenceDBPartition *sdbp= NULL;

  sdbp = (tSequenceDBPartition *)safe_malloc(sizeof(tSequenceDBPartition));
  sdbp->sequenceDBPath = strdup(sequenceDBPath);

  sprintf(buffer,"%s/seqDB.%d.%d", sequenceDBPath,revision,partition);
  testfp = fopen(buffer,"r");
  if(!testfp){
    sprintf(buffer,"openSequenceDBPartition: couldn't open index file %s", buffer);
    perror(buffer);
    exit(1);
  }    
  sdbp->multiAligns = CreateFromFileVA_tMARecord(testfp,10);
  fclose(testfp);

  sprintf(buffer,"%s/seqDB.data.%d.%d", sequenceDBPath,revision,partition);
  sdbp->datafp = fopen(buffer,"r");

  //  Build a hash of storeID to ma pointer.

  int i;
  int numMultiAligns;

  numMultiAligns = GetNumtMARecords(sdbp->multiAligns);
  sdbp->index = CreateScalarHashTable_AS(numMultiAligns + 1);

  for(i = 0; i < numMultiAligns; i++){
    tMARecord *maRecord = GettMARecord(sdbp->multiAligns, i);
    MultiAlignT *ma;
    int reference;
    AS_UTL_fseek(sdbp->datafp,maRecord->offset,SEEK_SET);
    ma = LoadMultiAlignTFromStream(sdbp->datafp,&reference);
    if(InsertInHashTable_AS(sdbp->index,
                            (uint64)maRecord->storeID, 0,
                            (uint64)ma, 0) != HASH_SUCCESS)
      assert(0);
  }

  return(sdbp);
}


void
closeSequenceDBPartition(tSequenceDBPartition *partition){
  fclose(partition->datafp);
  Delete_VA(partition->multiAligns);
  safe_free(partition->sequenceDBPath);
  safe_free(partition);
}


MultiAlignT *
loadFromSequenceDBPartition(tSequenceDBPartition *partition, int32 indx){
  MultiAlignT *ma = (MultiAlignT *)LookupValueInHashTable_AS(partition->index, indx, 0);

  if (ma == NULL)
    fprintf(stderr,"loadFromSequenceDBPartition()-- Multialign %d is NOT in this partition!\n", indx);
  assert(ma != NULL);

  return ma;
}

