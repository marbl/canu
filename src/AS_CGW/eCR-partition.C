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

const char *mainid = "$Id: eCR-partition.C,v 1.3 2009-10-05 22:49:42 brianwalenz Exp $";

#include "eCR.h"
#include "ScaffoldGraph_CGW.h"


int
main(int argc, char **argv) {
  int   ckptNum            = -1;
  int   numPartRequested   = 0;
  char *partInfoName       = NULL;
  FILE *partInfoFile       = NULL;

  argc = AS_configure(argc, argv);

  GlobalData = new Globals_CGW();

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-N") == 0) {
      numPartRequested = atoi(argv[++arg]);
      
    } else if (strcmp(argv[arg], "-p") == 0) {
      partInfoName = argv[++arg];

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }

  if (numPartRequested == 0)
    err++;
  if (partInfoName == NULL)
    err++;

  if ((GlobalData->outputPrefix[0] == 0) ||
      (GlobalData->gkpStoreName[0] == 0) ||
      (GlobalData->tigStoreName[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s [opts] -g gkpStore -n ckpNumber -c ckpName -N numPart -M maxFrag\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g gkpStore  The gatekeeper store\n");
    fprintf(stderr, "  -n ckpNumber The checkpoint to use\n");
    fprintf(stderr, "  -c ckpName   Use ckpName as the checkpoint name\n");
    fprintf(stderr, "  -N numPart   Number of partitions to make\n");
    fprintf(stderr, "  -M maxFrag   Maximum fragments per partition\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p partOut   Partition information output file\n");
    fprintf(stderr, "\n");

  if (numPartRequested == 0)
    fprintf(stderr, "%s: ERROR!  No number of partitions (-N) supplied.\n", argv[0]);

  if (partInfoName == 0)
    fprintf(stderr, "%s: ERROR!  No partition information output file (-p) supplied.\n", argv[0]);

    exit(1);
  }

  errno = 0;
  partInfoFile = fopen(partInfoName, "w");
  if (errno)
    fprintf(stderr, "%s: Failed to open partition information output file '%s': %s\n",
            argv[0], partInfoName, strerror(errno)), exit(1);

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, ckptNum, TRUE);

  //
  //  Scan all the scaffolds, build the partition mapping.
  //

  gkStore  *gkp       = new gkStore(GlobalData->gkpStoreName, FALSE, FALSE);
  short    *partition = new short [gkp->gkStore_getNumFragments() + 1];

  for (uint32 i=0; i<gkp->gkStore_getNumFragments() + 1; i++)
    partition[i] = -1;

  uint32    totFrags  = 0;
  uint32   *frgPerScf = new uint32 [GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) + 1];

  for (int32 sid=1; sid<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldT    *scf      = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    ContigT        *ctg      = NULL;
    uint32          nf       = 0;

    frgPerScf[sid] = 0;

    if ((isDeadCIScaffoldT(scf)) ||
        (scf->type != REAL_SCAFFOLD) ||
        (scf->info.Scaffold.numElements < 2))
      continue;

    for (ctg = GetGraphNode(ScaffoldGraph->ContigGraph, scf->info.Scaffold.AEndCI);
         ctg;
         ctg = (ctg->BEndNext == -1) ? NULL : GetGraphNode(ScaffoldGraph->ContigGraph, ctg->BEndNext)) {
      MultiAlignT  *ma = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);

      nf += GetNumIntMultiPoss(ma->f_list);
    }

    totFrags       += nf;
    frgPerScf[sid]  = nf;
  }

  short     curPart   = 0;

  uint32    aveFrags  = totFrags / numPartRequested + 1;
  uint32    curFrags  = totFrags;
  int32     lastID    = 0;

  fprintf(stderr, "TOTAL %d fragments\n", gkp->gkStore_getNumFragments());
  fprintf(stderr, "PLACE %d fragments into %d partitions -> %d fragments per partition.\n",
          totFrags, numPartRequested, aveFrags);

  if (totFrags == 0)
    goto allDone;

  totFrags = 0;

  for (int32 sid=1; sid<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldT    *scf      = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    ContigT        *ctg      = NULL;

    if ((isDeadCIScaffoldT(scf)) ||
        (scf->type != REAL_SCAFFOLD) ||
        (scf->info.Scaffold.numElements < 2))
      continue;

    //  Too many for the current partition?  We could optimize this to give us optimality (some
    //  partitions smaller than average, some larger, but exactly N partitions)...but this is close
    //  enough.  We'll let the early partitions be big.
    //
    //  To make them small, and probably return N+1 partitions, add in "frgPerScf[sid] to curFrags below.

    if (curFrags >= aveFrags) {
      if (curPart > 0) {
        fprintf(stderr,       "PARTITION %d ends with scffold %d and contains %d fragments (tot=%d).\n", curPart, lastID, curFrags, totFrags);
        fprintf(partInfoFile, "\t%d\t%d\t%d\n", curPart, lastID, curFrags);
      }
      curPart++;
      curFrags = 0;

      fprintf(stderr, "PARTITION %d starts with scaffold %d (tot=%d).\n", curPart, scf->id, totFrags);
      fprintf(partInfoFile, "%d\t%d", curPart, scf->id);
    }

    //  Set the partitions.

    totFrags += frgPerScf[sid];

    for (ctg = GetGraphNode(ScaffoldGraph->ContigGraph, scf->info.Scaffold.AEndCI);
         ctg;
         ctg = (ctg->BEndNext == -1) ? NULL : GetGraphNode(ScaffoldGraph->ContigGraph, ctg->BEndNext)) {
      MultiAlignT  *ma = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);

      for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++)
        partition[GetIntMultiPos(ma->f_list, i)->ident] = curPart;
    }

    //  Update stats

    curFrags += frgPerScf[sid];
    lastID    = scf->id;
  }

  assert(curFrags > 0);

  fprintf(stderr,       "PARTITION %d ends with scffold %d and contains %d fragments (tot=%d).\n", curPart, lastID, curFrags, totFrags);
  fprintf(partInfoFile, "\t%d\t%d\t%d\n", curPart, lastID, curFrags);

  //  Build partition

  //  BROKEN UNTIL GKPSTORE ALLOWS UPDATES TO PARTITIONS
  //
  //fprintf(stderr, "Building %d partitions.\n", curPart);
  //gkp->gkStore_buildPartitions(partition, curPart);

 allDone:
  delete gkp;

  DestroyScaffoldGraph(ScaffoldGraph);
  delete GlobalData;

  fclose(partInfoFile);

  exit(0);
}


