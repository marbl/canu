
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

static char *rcsid = "$Id$";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"


gkStoreStats::gkStoreStats(const char *gkStoreName) {
  gkStore   *gkp = new gkStore(gkStoreName, FALSE, FALSE, TRUE);

  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkStoreName);
    exit(1);
  }

  init(gkp);

  delete gkp;
}


gkStoreStats::gkStoreStats(gkStore *gkp) {
  init(gkp);
}


void
gkStoreStats::init(gkStore *gkp) {
  gkFragment    fr;
  gkStream     *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  numActiveFrag     = 0;
  numDeletedFrag    = 0;
  numMatedFrag      = 0;
  readLength        = 0;
  clearLength       = 0;

  lowestIID         = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  highestIID        = new uint32 [gkp->gkStore_getNumLibraries() + 1];

  numActivePerLib   = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numDeletedPerLib  = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numMatedPerLib    = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  readLengthPerLib  = new uint64 [gkp->gkStore_getNumLibraries() + 1];
  clearLengthPerLib = new uint64 [gkp->gkStore_getNumLibraries() + 1];

  while (fs->next(&fr)) {
    AS_IID     lib = fr.gkFragment_getLibraryIID();
    AS_IID     iid = fr.gkFragment_getReadIID();

    if (lowestIID[lib] == 0) {
      lowestIID[lib]  = iid;
      highestIID[lib] = iid;
    }
    if (highestIID[lib] < iid) {
      highestIID[lib] = iid;
    }

    if (fr.gkFragment_getIsDeleted()) {
      numDeletedFrag++;
      numDeletedPerLib[lib]++;
    } else {
      numActiveFrag++;
      numActivePerLib[lib]++;

      if (fr.gkFragment_getMateIID() > 0) {
        numMatedFrag++;
        numMatedPerLib[lib]++;
      }

      readLength             += fr.gkFragment_getSequenceLength();
      readLengthPerLib[lib]  += fr.gkFragment_getSequenceLength();

      clearLength            += fr.gkFragment_getClearRegionLength();
      clearLengthPerLib[lib] += fr.gkFragment_getClearRegionLength();
    }
  }

  delete fs;
}

gkStoreStats::~gkStoreStats() {
  delete [] lowestIID;
  delete [] highestIID;
  delete [] numActivePerLib;
  delete [] numDeletedPerLib;
  delete [] numMatedPerLib;
  delete [] readLengthPerLib;
  delete [] clearLengthPerLib;
}


