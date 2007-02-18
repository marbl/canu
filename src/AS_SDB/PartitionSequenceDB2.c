
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
static char CM_ID[] = "$Id: PartitionSequenceDB2.c,v 1.8 2007-02-18 14:04:50 brianwalenz Exp $";

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
#include "AS_SDB_SequenceDB.h"
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_MSG_pmesg.h"

/*
  PartitionSequenceDB:
  
  Given a SequenceDB and a partition size, measured in number of fragments,
  produce the following:
  1) A set of SequenceDBs called partition0-partitionN where:
  - partition0 contains all live contigs containing a single unitig (can grow to arbitrary size)
  Paritition 0 is really not very interesting and should be omitted later
  - partition1-partitionN contain contigs with >1 unitig, together will all unitigs and
  unitig surrogates referenced by the contig. This may cause a unitig to appear in more than
  one partition
  2) A fragment partition file specifying for each fragment which partition (possibly partition 0) its
  contig belongs.  This file is of the form
  <partitionID>  <fragID>


  This tool will operate in two phases, first partitioning the Contigs and generating a fragment partition file,
  and an in-memory unitig assignment list.  The second phase will read the unitigs and partition them according the the
  assignment list.  Partition 0 contigs/unitigs can be dealt with on the first pass.

  This file implements phase 2.

*/

void usage(void){

  fprintf(stderr,"* usage:  PartitionSDB2 <sdbname> <sdbRevision> <partitionFileName>\n");
  exit(1);
}


VA_TYPE(PtrT) *outputFiles = NULL;
VA_TYPE(PtrT) *outputIndices = NULL;

FILE *GetOutputFile(char *path, int revision, int partition){
  FILE **outputFilep;
  VA_TYPE(tMARecord) *index;

  if(!outputFiles){
    outputFiles = CreateVA_PtrT(100);
    outputIndices = CreateVA_PtrT(100);
  }
  outputFilep = ((FILE **)GetPtrT(outputFiles, partition));
  if(!outputFilep || !*outputFilep){
    char buffer[8096];
    FILE *outputFile;
    sprintf(buffer,"%s/seqDB.data.%d.%d",path, revision, partition);
    fprintf(stderr,"* Opening file %s\n", buffer);
    outputFile = fopen(buffer,"w");
    SetPtrT(outputFiles, partition,(void *)&outputFile);
    index = CreateVA_tMARecord(100);
    SetPtrT(outputIndices, partition, (void *)&index);
    fprintf(stderr,"* Inserted index %d as %p\n", partition, index);
    return outputFile;
  }else{
    return *outputFilep;
  }
}


VA_TYPE(tMARecord) *GetIndex(int partition){
  VA_TYPE(tMARecord) *ret, **retp;
  retp = (VA_TYPE(tMARecord) **)GetPtrT(outputIndices, partition);
  ret = *retp;
  fprintf(stderr,"* Index %d is %p %p\n", partition,retp, ret);

  return ret;
}

int main(int argc, char **argv){
  char *storeName;
  char *partitionFileName;
  int storeVersion;
  FILE *partitionfp;
  tSequenceDB *sequenceDB;
  char buffer [2048];
  CDS_CID_t binID, unitigID;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  CDS_CID_t lastBin, lastUnitig;
  tMARecord mar;
  int i;
  int forced=0;

  if(argc < 4)
    usage();


  storeName = argv[1];
  storeVersion = atoi(argv[2]);
  partitionFileName = argv[3];
  if(argc>4){
    if(strcmp("-F",argv[4])==0){
      forced=1;
    }
  }

  partitionfp = fopen(partitionFileName,"r");
  AssertPtr(partitionfp);


  fprintf(stderr,"* PartitionSequenceDB2 %s version %d partition %s\n", storeName, storeVersion, partitionFileName);

  sequenceDB = OpenSequenceDB(storeName, FALSE, storeVersion);

  lastBin = lastUnitig = -1;
  while(fgets(buffer,1999,partitionfp)){
    if(2 == sscanf(buffer,F_CID " " F_CID, &unitigID, &binID)){
      FILE *fp = GetOutputFile(storeName, storeVersion, binID);
      VA_TYPE(tMARecord) *index = GetIndex(binID);
      fprintf(stderr,"* Unitig " F_CID " ==> bin " F_CID "\n",
              unitigID, binID);
      fflush(NULL);

      if(lastBin == binID && lastUnitig == unitigID){
	fprintf(stderr,
                "* Attempt to Duplicate unitig " F_CID " > 1 time in bin " F_CID "...skipping\n",
		unitigID, binID);
	continue;
      }
      lastUnitig = unitigID;
      lastBin = binID;
      mar.flags.all = 0;
      mar.storeID = unitigID;
      mar.offset = CDS_FTELL(fp);
      AppendtMARecord(index, &mar);
      ReLoadMultiAlignTFromSequenceDB(sequenceDB, ma, unitigID, TRUE);    // This will load the ma from the disk file
      SaveMultiAlignTToStream(ma,fp);
    }
  }

  //  If the input was empty, we should not crash and burn, but we still exit
  //  ungracefully.
  //
  if (outputIndices == NULL) {
    fprintf(stderr, "* Found no partitions in the input (empty input file?) -- I fail!\n");
    exit(forced!=1);
  } else {
    for(i = 1; i < GetNumPtrTs(outputIndices); i++){
      char buffer[2048];
      FILE *indexfp;
      VA_TYPE(tMARecord) *index = GetIndex(i);
      fclose(GetOutputFile(storeName,storeVersion, i));
      sprintf(buffer,"%s/seqDB.%d.%d",storeName, storeVersion, i);
      indexfp = fopen(buffer,"w");
      CopyToFileVA_tMARecord(index, indexfp);
      fclose(indexfp);
    }
  }
  return 0;
}
