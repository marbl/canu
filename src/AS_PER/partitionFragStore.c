
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
#include <stdlib.h> 
#include <stdio.h> 
#include <assert.h>
#include <unistd.h>

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_fragStore_private.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_UTL_Var.h"

int32 ProcessPartitionFile(FILE *partitionfp,
                           CDS_IID_t firstFragID, CDS_IID_t lastFragID,
                           int32 *fragBins);

int main(int argc, char *argv[]){

  FILE *partitionfp;
  FragStoreHandle source;
  FragStoreHandle target;
  FragStreamHandle jane;
  ReadStructp myRead;
  int numPartitions;
  CDS_IID_t firstFragID, lastFragID;
  StoreStat stats;
  char *sourceStorePath, *targetStorePath, *partitionFilePath;
  FragStore *sourceStore;
  int32 *fragBins;
  int32 handleBinZeroSpecial = FALSE;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
                                    "z")) != EOF)){
      switch(ch){
        case 'z':
          handleBinZeroSpecial = TRUE;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }
  }
  if(argc - optind != 3){
    fprintf(stderr,"Usage: %s <PartitionFilePath> <SourceStorePath> <PartitionedStorePath>\n",
            argv[0]);
    exit(1);
  }
  partitionFilePath = argv[optind++];   
  sourceStorePath = argv[optind++];
  targetStorePath = argv[optind++];

  partitionfp = fopen(partitionFilePath, "r");
  if(!partitionfp){
    fprintf(stderr,"* Couldn't open partition file %s\n", partitionFilePath);
    exit(1);
  }
  source = openFragStore(sourceStorePath,"r");

  sourceStore =  FragStore_myStruct(source);
  statsStore(sourceStore->fragStore, &stats);
  firstFragID = (CDS_IID_t) stats.firstElem;
  lastFragID = (CDS_IID_t) stats.lastElem;

  fprintf(stderr,"* Allocated fragBins to size " F_IID "\n", lastFragID + 1);
  fragBins = (int32 *)calloc(lastFragID + 1, sizeof(int32));
  numPartitions = ProcessPartitionFile(partitionfp, firstFragID, lastFragID, fragBins);

  target = createPartitionedFragStore(targetStorePath,"partitioned",stats.firstElem,numPartitions);

  fprintf(stderr,"* Opened all files\n");
  fprintf(stderr,"* Processing partition file\n");

  myRead =  new_ReadStruct();

  fprintf(stdout,"* Partitioning fragStore %s (" F_S64 "," F_S64 ") (each dot represents 64k fragments)\n",
          sourceStorePath,getFirstElemFragStore(source), getLastElemFragStore(source));

  jane = openFragStream(source,NULL,0);

  while(nextFragStream(jane, myRead, FRAG_S_ALL)){
    uint32 id;
    int32 bin;
    getReadIndex_ReadStruct(myRead, &id);
    bin = fragBins[id];
    if((id & 0x0000ffff) == 0){
      fprintf(stderr,".");
      if((id & 0x000fffff) == 0){
	fprintf(stderr,"\n");
      }
      fflush(stderr);
    }
    //     fprintf(stderr,"* Inserting frag %u in bin %d\n", id,bin);
    //     fflush(stderr);
    if(bin == 0 && handleBinZeroSpecial){
      FragRecord *fr = (FragRecord *)myRead;
      fr->frag.deleted = TRUE;
      setReadType_ReadStruct(myRead, (FragType) 0);
      setSequence_ReadStruct(myRead, NULL, NULL);
      setSource_ReadStruct(myRead, "");
      setEntryTime_ReadStruct(myRead,0);
      setClearRegion_ReadStruct(myRead,0,0,READSTRUCT_ORIGINAL);
      setLocID_ReadStruct(myRead, 0);
      setLocalePos_ReadStruct(myRead, 0,0);
    }
    appendFragStorePartition(target,myRead,bin);  
  }

  fprintf(stdout,"*\n Closing target, source\n");

  closeFragStore(target);
  closeFragStream(jane);
  closeFragStore(source);
  fprintf(stdout,"* Bye Bye\n");

  exit(0);
}


int32 ProcessPartitionFile(FILE *partitionfp,
                           CDS_IID_t firstID, CDS_IID_t lastID,
                           int32 *fragBins){
  int numFrags = 0;
  CDS_IID_t fragID;
  int32 binID;
  int maxBin = 0;
  VA_TYPE(int32) *binFrags = CreateVA_int32(100);
  char buffer[2000];

  fprintf(stderr,"* Processing Partition File (first:%d  last:%d) (each . represents 64k lines)\n", firstID, lastID);
  fflush(stderr);
  while(fgets(buffer,1999,partitionfp)){
    if(2 == sscanf(buffer,"%d " F_IID, &binID, &fragID)){
      int32 *cnt;
      numFrags++;
      if((numFrags & 0x0000ffff) == 0){
        fprintf(stderr,".");
        if((numFrags & 0x000fffff) == 0){
          fprintf(stderr,"\n%d",numFrags);
        }
        fflush(stderr);
      }
    
      cnt = Getint32(binFrags, binID);
      if(!cnt){
        int32 dummy = 1;
        Setint32(binFrags,binID, &dummy);
      }else{
        (*cnt)++;
      }
      if(fragID > lastID || fragID < firstID){
        fprintf(stderr,"*** FATAL ERROR: Bad line %d:   %s\n*** fragID " F_IID " is OUT OF RANGE [" F_IID "," F_IID "] ...exiting\n",
                numFrags, buffer, fragID, firstID, lastID);
        exit(1);
      }

      fragBins[fragID] = binID;
      if(binID>maxBin)
        maxBin = binID;
    }else{
      fprintf(stderr,"* Skipped line %s\n", buffer);
      fflush(stderr);
    }
  }
  fprintf(stderr,"\n* Read %d frags first:" F_IID " last " F_IID " in %d partitions\n",
	  numFrags, firstID, lastID, maxBin);
  {
    int i;
    for(i = 0; i < GetNumint32s(binFrags); i++){
      fprintf(stderr,"* Bin %d has %d fragments\n",
	      i, (int) *Getint32(binFrags, i));
    }
  }
  return maxBin+1;
}

