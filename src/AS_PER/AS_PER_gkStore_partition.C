
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

static char *rcsid = "$Id: AS_PER_gkStore_partition.C,v 1.1 2009-10-28 17:27:29 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_fileIO.h"

void
gkStore::gkStore_loadPartition(uint32 partition) {
  char       name[FILENAME_MAX];
  int        i, f, e;

  assert(partmap    == NULL);
  assert(isCreating == 0);

  if (isReadOnly == 0)
    fprintf(stderr, "WARNING:  loading a partition from a writable gkpStore.\n");
  assert(isReadOnly == 1);

  partnum = partition;

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  if (AS_UTL_fileExists(name, FALSE, FALSE) == 0) {
    fprintf(stderr, "gkStore_loadPartition()--  Partition %d doesn't exist; normal store used instead.\n", partnum);
    return;
  }

  //  load all our data

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  partfpk = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/fnm.%03d", storePath, partnum);
  partfnm = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qnm.%03d", storePath, partnum);
  partqnm = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/fsb.%03d", storePath, partnum);
  partfsb = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qsb.%03d", storePath, partnum);
  partqsb = loadStorePartial(name, 0, 0);

  //  zip through the frags and build a map from iid to the frag record

  partmap = CreateScalarHashTable_AS();

  f = getFirstElemStore(partfpk);
  e = getLastElemStore(partfpk);
  for (i=f; i<=e; i++) {
    gkPackedFragment *p = (gkPackedFragment *)getIndexStorePtr(partfpk, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfnm);
  e = getLastElemStore(partfnm);
  for (i=f; i<=e; i++) {
    gkNormalFragment *p = (gkNormalFragment *)getIndexStorePtr(partfnm, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfsb);
  e = getLastElemStore(partfsb);
  for (i=f; i<=e; i++) {
    gkStrobeFragment *p = (gkStrobeFragment *)getIndexStorePtr(partfsb, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }
}



void
gkStore::gkStore_buildPartitions(short *partitionMap, uint32 maxPart) {
  gkFragment        fr;

  assert(partmap    == NULL);
  assert(isReadOnly == 1);
  assert(isCreating == 0);

  //  Create the partitions by opening N copies of the gatekeeper store,
  //  and telling each one to make a partition.
  //
  gkStore **gkpart = new gkStore * [maxPart + 1];

  AS_PER_setBufferSize(512 * 1024);

  for (uint32 i=1; i<=maxPart; i++)
    gkpart[i] = new gkStore (storePath, i);

  //  And, finally, add stuff to each partition.
  //
  for (uint32 iid=1; iid<=gkStore_getNumFragments(); iid++) {
    int p;

    gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    p = partitionMap[iid];

    //  Check it's actually partitioned.  Deleted reads won't get
    //  assigned to a partition.
    if (p < 1)
      continue;


    if (fr.type == GKFRAGMENT_PACKED) {
      appendIndexStore(gkpart[p]->partfpk, &fr.fr.packed);
    }


    if (fr.type == GKFRAGMENT_NORMAL) {
      fr.fr.normal.seqOffset = -1;
      fr.fr.normal.qltOffset = getLastElemStore(gkpart[p]->partqnm) + 1;

      appendIndexStore(gkpart[p]->partfnm, &fr.fr.normal);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqnm, fr.enc, fr.fr.normal.seqLen);
    }


    if (fr.type == GKFRAGMENT_STROBE) {
      fr.fr.strobe.seqOffset = -1;
      fr.fr.strobe.qltOffset = getLastElemStore(gkpart[p]->partqsb) + 1;

      appendIndexStore(gkpart[p]->partfsb, &fr.fr.strobe);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqsb, fr.enc, fr.fr.strobe.seqLen);
    }
  }

  //  cleanup -- close all the stores

  for (uint32 i=1; i<=maxPart; i++)
    delete gkpart[i];

  delete gkpart;
}
