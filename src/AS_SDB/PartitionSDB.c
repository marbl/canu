
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"

//  Given a SequenceDB and a partition size, measured in number of
//  fragments, produce the following:
//
//  1) A set of SequenceDBs called partition1-partitionN where:
//
//    - partition1 - partitionN contain contigs with more than one
//    unitig (see -all below), together will all unitigs and unitig
//    surrogates referenced by the contig. This may cause a unitig to
//    appear in more than one partition
//
//    - If -all is present, then all live contigs containing a single
//    unitig are included.  This is only interesting if consensus is
//    run to output VAR records.
//
//    - The partition file, UnitigPartition.txt, is the input to the
//    second stage.
//
//  2) A fragment partition file specifying for each fragment which
//  partition (possibly partition 0) its contig belongs.  This file,
//  FragPartition.txt, is of the form <partitionID> <fragID>
//
//  This tool will operate in two phases, first partitioning the
//  Contigs and generating a fragment partition file, and an in-memory
//  unitig assignment list.  The second phase will read the unitigs
//  and partition them according the the assignment list.

typedef struct {
  int elemID;
  int partitionID;
} tPartitionElement;
VA_DEF(tPartitionElement)


FILE *
openInputFile(char *inputFileName) {
  FILE *I;

  errno = 0;
  I = fopen(inputFileName, "r");
  if (errno)
    fprintf(stderr, "Can't open inputFileName '%s': %s\n", inputFileName, strerror(errno)), exit(1);
  return(I);
}

FILE *
openOutputPartition(char *inputFileName, int currentPartition) {
  FILE *O;
  char  N[1024];

  sprintf(N, "%s.%d", inputFileName, currentPartition);
  errno = 0;
  O = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Can't open output partition '%s': %s\n", N, strerror(errno)), exit(1);
  return(O);
}


int
ComparePartitionElements(const void *p1, const void *p2) {
  const tPartitionElement *pe1 = (const tPartitionElement *)p1;
  const tPartitionElement *pe2 = (const tPartitionElement *)p2;

  int diff = pe1->elemID - pe2->elemID;
  if (diff)
    return(diff);
  return(pe1->partitionID - pe2->partitionID);
}



int
main(int argc, char **argv) {
  char  *storeName         = NULL;
  int    storeVersion      = 0;
  int    fragsPerPartition = 0;
  char  *inputFileName     = NULL;
  int    includeAllUnitigs = 0;

  int     skippedFragments = 0;
  int     skippedUnitigs = 0;
  int     skippedContigs = 0;

  int     totalFragments  = 0;
  int     totalUnitigs    = 0;
  int     totalContigs    = 0;

  int     partitionFragments = 0;
  int     currentPartition = 1;

  VA_TYPE(tPartitionElement) *utgPartitionElems;
  VA_TYPE(tPartitionElement) *frgPartitionElems;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-seqstore") == 0) {
      storeName = argv[++arg];
    } else if (strcmp(argv[arg], "-version") == 0) {
      storeVersion = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-fragsper") == 0) {
      fragsPerPartition = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-input") == 0) {
      inputFileName     = argv[++arg];
    } else if (strcmp(argv[arg], "-all") == 0) {
      includeAllUnitigs = 1;
    }
    arg++;
  }

  if ((storeName == NULL) ||
      (storeVersion == 0) ||
      (fragsPerPartition == 0) ||
      (inputFileName == NULL)) {
    fprintf(stderr, "usage: %s [opts] -seqstore s.SeqStore -version N -fragsper M -i asm.cgw_contigs\n", argv[0]);
    fprintf(stderr, "  -all      Include all unitigs, including unitigs in single-unitig-contigs.\n");
    exit(1);
  }

  FILE         *inputFile         = openInputFile(inputFileName);
  FILE         *outputPartition   = openOutputPartition(inputFileName, currentPartition);
  tSequenceDB  *sequenceDB        = openSequenceDB(storeName, FALSE, storeVersion);

  if (((numContigsInSequenceDB(sequenceDB) == 0) && (includeAllUnitigs == 0)) ||
      ((numUnitigsInSequenceDB(sequenceDB) == 0))) {
    fprintf(stderr, "Nothing to do: %d unitigs, %d contigs.\n",
            numUnitigsInSequenceDB(sequenceDB),
            numContigsInSequenceDB(sequenceDB));
    exit(0);
  }

  utgPartitionElems = CreateVA_tPartitionElement(numUnitigsInSequenceDB(sequenceDB) + 100);
  frgPartitionElems = CreateVA_tPartitionElement(numUnitigsInSequenceDB(sequenceDB) * 5);

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Step 1
  //
  GenericMesg  *pmesg;

  while (ReadProtoMesg_AS(inputFile, &pmesg) != EOF) {
    if (pmesg->t == MESG_ICM) {
      IntConConMesg *ma = pmesg->m;
      int nunitigs = ma->num_unitigs;
      int nfrags   = ma->num_pieces;

      if ((includeAllUnitigs) || (nunitigs > 1)) {
        int u, f;



        for (u=0; u<nunitigs; u++) {
          tPartitionElement partElem;
          partElem.elemID      = ma->unitigs[u].ident;
          partElem.partitionID = currentPartition;
          AppendtPartitionElement(utgPartitionElems, &partElem);
        }	    

        for (f=0; f<nfrags; f++) {
          tPartitionElement partElem;
          partElem.elemID      = ma->pieces[f].ident;
          partElem.partitionID = currentPartition;
          AppendtPartitionElement(frgPartitionElems, &partElem);
        }

        WriteProtoMesg_AS(outputPartition, pmesg);

        partitionFragments += nfrags;

        totalFragments     += nfrags;
        totalUnitigs       += nunitigs;
        totalContigs++;

        if (partitionFragments > fragsPerPartition) {
          currentPartition++;

          fclose(outputPartition);
          outputPartition = openOutputPartition(inputFileName, currentPartition);

          partitionFragments = 0;
        }
      } else {
        skippedFragments += nfrags;
        skippedUnitigs   += nunitigs;
        skippedContigs++;
      }
    }
  }

  fclose(outputPartition);

  qsort(GettPartitionElement(utgPartitionElems, 0),
        GetNumtPartitionElements(utgPartitionElems),
        sizeof(tPartitionElement),
        ComparePartitionElements);
  qsort(GettPartitionElement(frgPartitionElems, 0),
        GetNumtPartitionElements(frgPartitionElems),
        sizeof(tPartitionElement),
        ComparePartitionElements);

  {
    FILE *output = NULL;
    int i;

    errno = 0;
    output = fopen("UnitigPartition.txt", "w");
    if (errno)
      fprintf(stderr, "Can't open output 'UnitigPartition.txt': %s\n", strerror(errno)), exit(1);

    for (i=0; i<GetNumtPartitionElements(utgPartitionElems); i++) {
      tPartitionElement *pElem = GettPartitionElement(utgPartitionElems,i);
      fprintf(output, "%d %d\n", pElem->elemID, pElem->partitionID);
    }

    errno = 0;
    output = fopen("FragPartition.txt", "w");
    if (errno)
      fprintf(stderr, "Can't open output 'FragPartition.txt': %s\n", strerror(errno)), exit(1);

    for (i=0; i<GetNumtPartitionElements(frgPartitionElems); i++) {
      tPartitionElement *pElem = GettPartitionElement(frgPartitionElems,i);
      fprintf(output, "%d %d\n", pElem->partitionID, pElem->elemID);
    }
  }

  fprintf(stderr, "Fragments:  %8d skipped   %8d partitioned\n", skippedFragments, totalFragments);
  fprintf(stderr, "Unitigs:    %8d skipped   %8d partitioned   %8d input\n", skippedUnitigs, totalUnitigs, numUnitigsInSequenceDB(sequenceDB));
  fprintf(stderr, "Contigs:    %8d skipped   %8d partitioned   %8d input\n", skippedContigs, totalContigs, numContigsInSequenceDB(sequenceDB));


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Step 2
  //
  int i;

  //  We know there are 'currentPartition + 1' partitions, so we can
  //  allocate all our output files ahead of time.
  //
  currentPartition++;
  FILE               **pOutput = (FILE **)              safe_calloc(sizeof(FILE *),               currentPartition);
  VA_TYPE(tMARecord) **pIndex  = (VA_TYPE(tMARecord) **)safe_calloc(sizeof(VA_TYPE(tMARecord) *), currentPartition);

  for (i=1; i<currentPartition; i++) {
    char  N[1024];

    sprintf(N, "%s/seqDB.v%03d.dat.p%03d", storeName, storeVersion, i);

    errno = 0;
    pOutput[i] = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Can't open SDB partition '%s': %s\n", N, strerror(errno)), exit(1);

    pIndex[i] = CreateVA_tMARecord(100);
  }

  int lastBin = -1;
  int lastUtg = -1;

  MultiAlignT *ma = CreateEmptyMultiAlignT();

  for (i=0; i<GetNumtPartitionElements(utgPartitionElems); i++) {
    tPartitionElement *pElem = GettPartitionElement(utgPartitionElems, i);

    int binID = pElem->partitionID;
    int utgID = pElem->elemID;

    if ((lastBin == binID) && 
        (lastUtg == utgID)) {
      fprintf(stderr, "* Attempt to Duplicate unitig "F_CID" in bin "F_CID"...skipping\n", utgID, binID);
      continue;
    }

    assert(binID > 0);
    assert(binID < currentPartition);

    tMARecord mar = {0};

    mar.storeID      = 0;
    mar.multiAlignID = utgID;
    mar.isDeleted    = 0;
    mar.offset       = AS_UTL_ftell(pOutput[binID]);

    AppendtMARecord(pIndex[binID], &mar);

    copyMultiAlignTFromSequenceDB(sequenceDB, ma, utgID, TRUE);
    SaveMultiAlignTToStream(ma, pOutput[binID]);

    lastUtg = utgID;
    lastBin = binID;
  }


  for (i=1; i<currentPartition; i++) {
    char N[1024];

    sprintf(N, "%s/seqDB.v%03d.dat.i%03d", storeName, storeVersion, i);
    errno = 0;
    FILE *F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Can't open SDB partition '%s': %s\n", N, strerror(errno)), exit(1);
    CopyToFileVA_tMARecord(pIndex[i], F);
    fclose(F);

    fclose(pOutput[i]);
  }

  return(0);
}
