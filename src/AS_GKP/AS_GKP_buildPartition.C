
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

static char const *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.H"
#include "AS_GKP_include.H"
#include "AS_PER_gkpStore.H"
#include "AS_PER_encodeSequenceQuality.H"



void
Build_Partition(char      *gatekeeperName,
                char      *partitionFile,
                int32      flags) {

  gkStore *gkp       = new gkStore(gatekeeperName, FALSE, FALSE);
  uint32   maxIID    = gkp->gkStore_getNumFragments() + 1;
  uint32   maxPart   = 0;
  short   *partition = new short [maxIID];

  for (uint32 i=0; i<maxIID; i++)
    partition[i] = -1;

  //  Read the partition file, remembering the highest partition
  //
  errno = 0;
  FILE *F = fopen(partitionFile, "r");
  if (errno) {
    fprintf(stderr, "GKP Error: Build_Partition()-- failed to open '%s': %s\n", partitionFile, strerror(errno));
    exit(1);
  }

  while (!feof(F)) {
    int  i, p;
    if (2 == fscanf(F, " %d %d ", &p, &i)) {
      partition[i] = p;
      if (p > maxPart)
        maxPart = p;
    }
  }
  fclose(F);

  gkp->gkStore_buildPartitions(partition, maxPart);

  delete [] partition;
  delete    gkp;
}


