
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

static char CM_ID[] = "$Id: AS_GKP_buildPartition.c,v 1.1 2007-03-05 05:57:16 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"



int
Build_Partition(char      *gatekeeperName,
                char      *partitionFile,
                int32      flags) {

  GateKeeperStore  *gkp    = NULL;
  int               iid    = 0;
  int               maxIID = 0;

  short            *partition;
  int               maxPart;

  FILE             *F   = NULL;

  GateKeeperStore **gkpart;

  fragRecord       *fr  = new_fragRecord();
  StoreStat         stats;
  char              encodedsequence[AS_FRAG_MAX_LEN];

  //  Open the gatekeeper, read-only
  //
  gkp    = openGateKeeperStore(gatekeeperName, FALSE);
  maxIID = getLastElemFragStore(gkp) + 1;

  //  Read the partition file, remembering the highest partition
  //
  errno = 0;
  F = fopen(partitionFile, "r");
  if (errno) {
    fprintf(stderr, "Build_Partition()-- failed to open '%s': %s\n", partitionFile, strerror(errno));
    exit(1);
  }

  partition = (short *)safe_malloc(sizeof(short) * maxIID);

  for (iid=0; iid<maxIID; iid++)
    partition[iid] = -1;

  while (!feof(F)) {
    int  i, p;
    if (2 == fscanf(F, " %d %d ", &i, &p)) {
      partition[i] = p;
      if (p > maxPart)
        maxPart = p;
    }
  }
  fclose(F);

  //  No, really, it's one more.
  maxPart++;

  //  Create the partitions by opening N copies of the gatekeeper store,
  //  and telling each one to make a partition.
  //
  gkpart = (GateKeeperStore **)safe_malloc(sizeof(GateKeeperStore *) * maxPart);

  for (iid=0; iid<maxPart; iid++) {
    gkpart[iid] = openGateKeeperStore(gatekeeperName, FALSE);
    createGateKeeperPartition(gkpart[iid], iid);
  }


  //  And, finally, add stuff to each partition.
  //
  for (iid=1; iid<maxIID; iid++) {
    int p;

    getFrag(gkp, iid, fr, flags);

    p = partition[iid];

    //  update pointers

    fr->gkfr.seqOffset = -1;

    statsStore(gkpart[p]->partqlt, &stats);
    fr->gkfr.qltOffset = stats.lastElem;

    statsStore(gkpart[p]->parthps, &stats);
    fr->gkfr.hpsOffset = stats.lastElem;

    statsStore(gkpart[p]->partsrc, &stats);
    fr->gkfr.srcOffset = stats.lastElem;

    //  append the elements

    appendGateKeeperFragmentStore(gkpart[p]->partfrg, &fr->gkfr);

    encodeSequenceQuality(encodedsequence, fr->seq, fr->qlt);
    appendVLRecordStore(gkpart[p]->partqlt, encodedsequence, fr->gkfr.seqLen);

    appendVLRecordStore(gkpart[p]->parthps, NULL,    0);

    appendVLRecordStore(gkpart[p]->partsrc, fr->src, fr->gkfr.srcLen);
  }


  //  cleanup -- close all the stores

  for (iid=0; iid<maxPart; iid++)
    closeGateKeeperStore(gkpart[iid]);

  closeGateKeeperStore(gkp);

  del_fragRecord(fr);
}


