
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

//  scaffold2fasta - dumps scaffolds from a checkpoint as fasta
//
//  [-s scaffoldNum] -c ckpName -n ckpNum -f frgStore -g gkpStore

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "Globals_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

#include "CommonREZ.h"
#include "GapWalkerREZ.h"  // FindGapLength

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-s begin-scaffold] [-e end-scaffold] -c ckpName -n ckpNum -f frgStore -g gkpStore\n",
          name);
  exit(1);
}

int
main(int argc, char **argv) {
  int   startScaff   = 0;
  int   endScaff     = CDS_INT32_MAX;
  int   ckptNum      = 0;
  char *ckptFileName = 0L;
  char *frgStoreName = 0L;
  char *gkpStoreName = 0L;

  int arg=1;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-s") == 0) {
      startScaff = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      endScaff = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      ckptFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-f") == 0) {
      frgStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      usage(argv[0]);
    }
    arg++;
  }

  if ((ckptNum == 0) || (ckptFileName == 0L) || (frgStoreName == 0L) || (gkpStoreName == 0L))
    usage(argv[0]);

  GlobalData            = CreateGlobal_CGW();
  GlobalData->stderrc   = stderr;
  GlobalData->stderro   = stderr;
  GlobalData->stderrfp  = stderr;

  strcpy(GlobalData->File_Name_Prefix,      ckptFileName);
  strcpy(GlobalData->Frag_Store_Name,       frgStoreName);
  strcpy(GlobalData->Gatekeeper_Store_Name, gkpStoreName);

  //  Load the scaffold graph
  //
  fprintf(stderr, "Loading scaffold graph\n");
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix, ckptNum, FALSE);

  //  Get the scaffold we care about
  //
  NodeCGW_T  *scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, startScaff);

  while ((scaff) && (startScaff < endScaff)) {

    // not interested in dead, not real, or singleton scaffolds
    //
    int  skip = 0;

    if (isDeadCIScaffoldT(scaff)) {
      fprintf(stderr, "scaffold " F_CID " is dead\n", scaff->id);
      skip = 1;
    }
    if (scaff->type != REAL_SCAFFOLD) {
      fprintf(stderr, "scaffold " F_CID " is not real\n", scaff->id);
      skip = 1;
    }

    if (!skip) {
      CIScaffoldTIterator    CIsTemp;
      ChunkInstanceT        *nextContig      = NULL;
      static VA_TYPE(char)  *contigConsensus = NULL;
      static VA_TYPE(char)  *contigQuality   = NULL;
      char                  *sequence;

      contigConsensus = CreateVA_char(1024);
      contigQuality   = CreateVA_char(1024);
    
      fprintf(stdout, ">scaffold"F_CID" scaffold_length=%.f\n", 
              scaff->id,
              scaff->bpLength.mean);
    
      InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
      while (NextCIScaffoldTIterator(&CIsTemp)) {
        ChunkInstanceT  *contig = GetGraphNode( ScaffoldGraph->ContigGraph, CIsTemp.curr);      

        assert(contig != NULL);

        if (CIsTemp.next != NULLINDEX) {
          nextContig = GetGraphNode(ScaffoldGraph->ContigGraph, CIsTemp.next);
          assert(nextContig != NULL);
        }

        GetConsensus(ScaffoldGraph->ContigGraph, contig->id, contigConsensus, contigQuality);
        sequence = Getchar(contigConsensus, 0);

        fwrite(sequence, sizeof(char), strlen(sequence), stdout);
     
        if (CIsTemp.next != NULLINDEX) {
          char  *gapN = 0L;
          int    gapS = (int)FindGapLength(contig, nextContig, FALSE).mean;

          if (gapS < 20)
            gapS = 20;

          gapN = (char *)malloc(sizeof(char) * gapS);
          memset(gapN, 'N', gapS);
          fwrite(gapN, sizeof(char), gapS, stdout);
          free(gapN);
        }
      }
      fprintf(stdout, "\n");
    }

    startScaff++;
    scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, startScaff);
  }

  return(0);
}
