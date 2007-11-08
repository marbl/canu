
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

static char const *rcsid = "$Id: AS_GKP_buildPartition.c,v 1.9 2007-11-08 12:38:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
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

  short            *partition = NULL;
  int               maxPart   = 0;

  FILE             *F   = NULL;

  GateKeeperStore **gkpart = NULL;

  fragRecord        fr;
  char              encodedsequence[AS_FRAG_MAX_LEN];

  uint64            seqLen = 0;
  uint64            hpsLen = 0;
  uint64            srcLen = 0;

  //  Open the gatekeeper, read-only
  //
  gkp    = openGateKeeperStore(gatekeeperName, FALSE);
  maxIID = getLastElemFragStore(gkp) + 1;

  //  Read the partition file, remembering the highest partition
  //
  errno = 0;
  F = fopen(partitionFile, "r");
  if (errno) {
    fprintf(stderr, "GKP Error: Build_Partition()-- failed to open '%s': %s\n", partitionFile, strerror(errno));
    exit(1);
  }

  partition = (short *)safe_malloc(sizeof(short) * maxIID);

  for (iid=0; iid<maxIID; iid++)
    partition[iid] = -1;

  while (!feof(F)) {
    int  i, p;
    if (2 == fscanf(F, " %d %d ", &p, &i)) {
      partition[i] = p;
      if (p > maxPart)
        maxPart = p;
    }
  }
  fclose(F);


  //  Create the partitions by opening N copies of the gatekeeper store,
  //  and telling each one to make a partition.
  //
  gkpart = (GateKeeperStore **)safe_calloc(sizeof(GateKeeperStore *), maxPart + 1);

  AS_PER_setBufferSize(512 * 1024);

  for (iid=1; iid<=maxPart; iid++)
    gkpart[iid] = createGateKeeperPartition(gatekeeperName, iid);

  //  And, finally, add stuff to each partition.
  //
  for (iid=1; iid<maxIID; iid++) {
    int p;

    getFrag(gkp, iid, &fr, flags);

    p = partition[iid];

    //  Check it's actually partitioned.  Deleted reads won't get
    //  assigned to a partition.
    if (p == -1)
      continue;

    //  update pointers

    fr.gkfr.seqOffset = -1;

    fr.gkfr.qltOffset = getLastElemStore(gkpart[p]->partqlt);
    fr.gkfr.hpsOffset = getLastElemStore(gkpart[p]->parthps);
    fr.gkfr.srcOffset = getLastElemStore(gkpart[p]->partsrc);

    //  append the elements

    appendIndexStore(gkpart[p]->partfrg, &fr.gkfr);

    encodeSequenceQuality(encodedsequence, fr.seq, fr.qlt);
    appendStringStore(gkpart[p]->partqlt, encodedsequence, fr.gkfr.seqLen);

    appendStringStore(gkpart[p]->parthps, NULL,    0);

    appendStringStore(gkpart[p]->partsrc, fr.src, fr.gkfr.srcLen);

    seqLen += fr.gkfr.seqLen;
    hpsLen += 0;
    srcLen += fr.gkfr.srcLen;

    if ((iid % 100000) == 0) {
      fprintf(stderr, "frags:%6d  seqLen:%9.3fMbp  hpsLen:%9.3fM  srcLen:%9.3fM\n",
              iid,
              (double)seqLen / 1000000.0,
              (double)hpsLen / 1000000.0,
              (double)srcLen / 1000000.0);
    }
  }


  //  cleanup -- close all the stores

  for (iid=0; iid<=maxPart; iid++)
    closeGateKeeperStore(gkpart[iid]);

  closeGateKeeperStore(gkp);
}


