
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
static char CM_ID[] = "$Id: PartitionSequenceDB.c,v 1.4 2005-03-22 19:49:27 jason_miller Exp $";

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

*/

void usage(void){

  fprintf(stderr,"* usage:  PartitionSDB <sdbname> <sdbRevision> <fragsPerPartition> <inputFileName>\n");
  exit(1);
}



typedef struct {
  int elemID;
  int partitionID;
}tPartitionElement;

VA_DEF(tPartitionElement)

static int ComparePartitionElements(const void *p1, const void *p2){
  tPartitionElement *pe1 = (tPartitionElement *)p1;
  tPartitionElement *pe2 = (tPartitionElement *)p2;
  int diff;
  diff = pe1->elemID - pe2->elemID;
  if(diff)
    return diff;
  diff = pe1->partitionID - pe2->partitionID;

  return diff;

}


int main(int argc, char **argv){
  char *storeName;
  char *inputFileName;
  char outputFileName[4096];
  int storeVersion;
  GenericMesg *pmesg;
  FILE *inputFile, *outputFile;
  tSequenceDB *sequenceDB;
  char buffer [2048];
  int fragsPerPartition;
  int numUnitigs, numContigs;
  int i = 0;
  int deleted = 0;
  int totalFragments = 0;
  int totalContigs = 0;
  size_t totalContigSize = 0;
  int totalUnitigs = 0;
  int deletedContigs = 0;
  int partition0Contigs = 0;
  size_t partition0ContigSize = 0;
  int partition0Unitigs = 0;
  int partition0Fragments = 0;
  int partitionFragments = 0;
  int quanta;
  int currentPartition = 0;

  VA_TYPE(tPartitionElement) *unitigPartitionElems;
  VA_TYPE(tPartitionElement) *fragPartitionElems;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  MesgReader reader;

  if(argc < 4)
    usage();


  storeName = argv[1];
  storeVersion = atoi(argv[2]);
  fragsPerPartition = atoi(argv[3]);
  inputFileName = argv[4];

  sprintf(buffer,"%s.part", storeName);
  outputFile = fopen(buffer,"w");

  currentPartition = 1;

  inputFile = fopen(inputFileName,"r");
  assert(inputFile);

  reader = InputFileType_AS(inputFile);

  sprintf(outputFileName,"%s.%d", inputFileName,currentPartition);
  outputFile = fopen(outputFileName,"w");
  assert(outputFile);

  fprintf(stderr,"* OpenSequenceDB %s revision %d fragsPerPartition %d\n", storeName, storeVersion, fragsPerPartition);

  sequenceDB = OpenSequenceDB(storeName, FALSE, storeVersion);

  numUnitigs = NumUnitigsInSequenceDB(sequenceDB);
  numContigs = NumContigsInSequenceDB(sequenceDB);


  unitigPartitionElems = CreateVA_tPartitionElement(numUnitigs + 100);
  fragPartitionElems = CreateVA_tPartitionElement(numUnitigs * 5);

  quanta = numContigs/100;

 while (reader(inputFile,&pmesg) != EOF){
    int nfrags;
    int nunitigs;
    int failure;
    IntConConMesg *ma;

    if(pmesg->t != MESG_ICM)
      continue;

    ma = pmesg->m;
    i++;

    if(i%quanta == 0){
      fprintf(stderr,"* Loading contig %d\n", i);
      fflush(stderr);
    }
#if 0
    	failure = ReLoadMultiAlignTFromSequenceDB(sequenceDB, ma, i, FALSE);

	if(failure){
	  // ma doesn't exist
	  deletedContigs++;
	  continue;
	}
#endif	
	nunitigs = ma->num_unitigs;
	nfrags = ma->num_pieces;

	if(nunitigs > 1){
	  int u,f;
	  for(u = 0; u < nunitigs; u++){
	    tPartitionElement partElem;
	    IntUnitigPos *iup = ma->unitigs + u;

	    partElem.elemID = iup->ident;
	    partElem.partitionID = currentPartition;
	    AppendtPartitionElement(unitigPartitionElems, &partElem);
	  }	    

	  for(f = 0; f < nfrags; f++){
	    tPartitionElement partElem;
	    IntMultiPos *imp = ma->pieces + f;
	    partElem.elemID = imp->ident;
	    partElem.partitionID = currentPartition;
	    AppendtPartitionElement(fragPartitionElems, &partElem);
	  }


    WriteProtoMesg_AS(outputFile,pmesg);
    fflush(outputFile);

	  totalFragments += nfrags;
	  partitionFragments += nfrags;
	  totalUnitigs += nunitigs;
	  totalContigs++;

	  fprintf(stderr,"* totalFragments %d partitionFragments %d contigs:%d unitigs %d\n",
		  totalFragments, partitionFragments, totalContigs, totalUnitigs);

	  if(partitionFragments > fragsPerPartition){

	    currentPartition++;
	    fprintf(stderr,"Opening partition %d\n", currentPartition);
	    partitionFragments=0;
	    sprintf(outputFileName,"%s.%d", inputFileName,currentPartition);
	    outputFile = fopen(outputFileName,"w");
	    assert(outputFile);
	  }

	}else{
	  partition0Fragments += nfrags;
	  partition0Unitigs += nunitigs;
	  partition0Contigs++;
	}

  }
      qsort((void *)GettPartitionElement(unitigPartitionElems, 0), GetNumtPartitionElements(unitigPartitionElems), sizeof(tPartitionElement), ComparePartitionElements);
      qsort((void *)GettPartitionElement(fragPartitionElems, 0), GetNumtPartitionElements(fragPartitionElems), sizeof(tPartitionElement), ComparePartitionElements);

 {
   FILE *output;
   int i;
   output = fopen("UnitigPartition.txt", "w");
   assert(output);
   for(i = 0; i < GetNumtPartitionElements(unitigPartitionElems); i++){
     tPartitionElement *pElem = GettPartitionElement(unitigPartitionElems,i);
     fprintf(output,"%d %d\n", pElem->elemID, pElem->partitionID);
   }
 }


 {
   FILE *output;
   int i;
   output = fopen("FragPartition.txt", "w");
   assert(output);
   for(i = 0; i < GetNumtPartitionElements(fragPartitionElems); i++){
     tPartitionElement *pElem = GettPartitionElement(fragPartitionElems,i);
     fprintf(output,"%d %d\n", pElem->elemID, pElem->partitionID);
   }
 }

  fprintf(stderr,"* Processed %d contigs %d deleted %d to partition0 and %d to process\n",
	  numContigs, deletedContigs, partition0Contigs, totalContigs);
  fprintf(stderr,"* Contig Size   Partition0:" F_SIZE_T "     Process:" F_SIZE_T "\n",
	  partition0ContigSize, totalContigSize);
  fprintf(stderr,"* Processed %d unitigs    %d to partition 0 and %d to process\n",
	  totalUnitigs + partition0Unitigs, partition0Unitigs, totalUnitigs);
  fprintf(stderr,"* Processed %d fragments    %d to partition 0 and %d to process\n",
	  totalFragments + partition0Fragments, partition0Fragments, totalFragments);

}
