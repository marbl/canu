
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

static char *rcsid = "$Id: AS_PER_gkStore_partition.C,v 1.5 2012-02-29 00:23:15 brianwalenz Exp $";

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
  int64      i, f, e;

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
  sprintf(name,"%s/qpk.%03d", storePath, partnum);
  partqpk = loadStorePartial(name, 0, 0);

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
    if (InsertInHashTable_AS(partmap, (uint64)p->readIID, 0, i, 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfnm);
  e = getLastElemStore(partfnm);
  for (i=f; i<=e; i++) {
    gkNormalFragment *p = (gkNormalFragment *)getIndexStorePtr(partfnm, i);
    if (InsertInHashTable_AS(partmap, (uint64)p->readIID, 0, i, 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfsb);
  e = getLastElemStore(partfsb);
  for (i=f; i<=e; i++) {
    gkStrobeFragment *p = (gkStrobeFragment *)getIndexStorePtr(partfsb, i);
    if (InsertInHashTable_AS(partmap, (uint64)p->readIID, 0, i, 0) != HASH_SUCCESS)
      assert(0);
  }
}



void
gkStore::gkStore_buildPartitions(short *partitionMap, uint32 maxPart) {
  char              name[FILENAME_MAX];
  gkFragment        fr;

  assert(partmap    == NULL);
  assert(isReadOnly == 1);
  assert(isCreating == 0);

  //  Create the partitions by opening N copies of the data stores,
  //  and writing data to each.

  StoreStruct  **partfpk = new StoreStruct * [maxPart + 1];
  StoreStruct  **partqpk = new StoreStruct * [maxPart + 1];

  StoreStruct  **partfnm = new StoreStruct * [maxPart + 1];
  StoreStruct  **partqnm = new StoreStruct * [maxPart + 1];

  StoreStruct  **partfsb = new StoreStruct * [maxPart + 1];
  StoreStruct  **partqsb = new StoreStruct * [maxPart + 1];

  AS_PER_setBufferSize(512 * 1024);

  for (uint32 i=0; i<=maxPart; i++) {
    sprintf(name,"%s/fpk.%03d", storePath, i);
    partfpk[i] = createIndexStore(name, "partfpk", sizeof(gkPackedFragment), 1);
    sprintf(name,"%s/qpk.%03d", storePath, i);
    partqpk[i] = createIndexStore(name, "partqpk", sizeof(char) * inf.gkPackedSequenceSize, 1);

    sprintf(name,"%s/fnm.%03d", storePath, i);
    partfnm[i] = createIndexStore(name, "partfnm", sizeof(gkNormalFragment), 1);
    sprintf(name,"%s/qnm.%03d", storePath, i);
    partqnm[i] = createStringStore(name, "partqnm");

    sprintf(name,"%s/fsb.%03d", storePath, i);
    partfsb[i] = createIndexStore(name, "partfsb", sizeof(gkStrobeFragment), 1);
    sprintf(name,"%s/qsb.%03d", storePath, i);
    partqsb[i] = createStringStore(name, "partqsb");
  }

  for (uint32 iid=1; iid<=gkStore_getNumFragments(); iid++) {
    gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    int32 p = partitionMap[iid];

    if (p < 1)
      //  Deleted reads are not assigned a partition; skip them
      continue;

    if (fr.type == GKFRAGMENT_PACKED) {
      appendIndexStore(partfpk[p], &fr.fr.packed);
      appendIndexStore(partqpk[p],  fr.enc);
    }


    if (fr.type == GKFRAGMENT_NORMAL) {
      fr.fr.normal.seqOffset = -1;
      fr.fr.normal.qltOffset = getLastElemStore(partqnm[p]) + 1;

      appendIndexStore(partfnm[p], &fr.fr.normal);
      appendStringStore(partqnm[p], fr.enc, fr.fr.normal.seqLen);
    }


    if (fr.type == GKFRAGMENT_STROBE) {
      fr.fr.strobe.seqOffset = -1;
      fr.fr.strobe.qltOffset = getLastElemStore(partqsb[p]) + 1;

      appendIndexStore(partfsb[p], &fr.fr.strobe);
      appendStringStore(partqsb[p], fr.enc, fr.fr.strobe.seqLen);
    }
  }

  //  cleanup -- close all the stores

  for (uint32 i=0; i<=maxPart; i++) {
    closeStore(partfpk[i]);
    closeStore(partqpk[i]);
    closeStore(partfnm[i]);
    closeStore(partqnm[i]);
    closeStore(partfsb[i]);
    closeStore(partqsb[i]);
  }

  delete [] partfpk;
  delete [] partqpk;
  delete [] partfnm;
  delete [] partqnm;
  delete [] partfsb;
  delete [] partqsb;
}
