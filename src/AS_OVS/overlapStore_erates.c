
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

static char CM_ID[] = "$Id: overlapStore_erates.c,v 1.1 2007-03-13 22:38:51 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"



void
updateErates(char *storeName, char *eratesName) {
  int             i;
  char            name[FILENAME_MAX];

  OverlapStore   *store = NULL;
  OverlapStore   *orig  = NULL;
  OVSoverlap      ovl;

  int32          eFirst;
  int32          eLast;
  int32          eNum;

  FILE           *eF   = NULL;
  int             eLen = 0;
  int             eMax = 0;
  int16          *e    = NULL;

  //  Open the erates, read in the header.
  //
  errno = 0;
  eF = fopen(eratesName, "r");
  if (errno) {
    fprintf(stderr, "failed to open erates file '%s': %s\n", eratesName, strerror(errno));
    exit(1);
  }
  AS_UTL_safeRead(eF, &eFirst, "updateErates read header 0", sizeof(int32), 1);
  AS_UTL_safeRead(eF, &eLast,  "updateErates read header 1", sizeof(int32), 1);
  AS_UTL_safeRead(eF, &eNum,   "updateErates read header 2", sizeof(int32), 1);


  //  Open the two stores so we can read overlaps from them.  The
  //  store we're merging into is opened as a "copy".  The store we
  //  merge from is opened first, so that if it fails, we haven't
  //  turned the original store into a backup.
  //
  orig = AS_OVS_openOverlapStorePrivate(storeName, TRUE,  TRUE);

  assert(eNum == orig->ovs.numOverlapsTotal);

  //  Recreate a store in the same place as the original store.
  //
  store = AS_OVS_createOverlapStore(storeName, FALSE);

  //  Grab some space for our cache of erates
  e = (int16 *)safe_malloc(sizeof(int16) * 1048576);

  while (AS_OVS_readOverlapFromStore(orig, &ovl)) {
    if (eLen >= eMax) {
      eLen = 0;
      eMax = AS_UTL_safeRead(eF, e, "updateErates read erates", sizeof(uint16), 1048576);
      if (eMax == 0) {
        fprintf(stderr, "failed to read more erates from '%s'.\n", eratesName);
        exit(1);
      }
    }

    ovl.dat.ovl.corr_erate = e[eLen++];
    eNum--;

    fprintf(stderr, "%d-vs-%d from %d(%6.3f) to %d(%6.2f)\n",
            ovl.a_iid, ovl.b_iid,
            ovl.dat.ovl.orig_erate, Expand_Quality(ovl.dat.ovl.orig_erate),
            ovl.dat.ovl.corr_erate, Expand_Quality(ovl.dat.ovl.corr_erate));

    AS_OVS_writeOverlapToStore(store, &ovl);
  }

  //  check that both files are empty now

  assert(eNum == 0);
  

  //  ALL DONE!  Close the stores, nuke the backups and get outta here.
  //
  AS_OVS_closeOverlapStore(orig);
  AS_OVS_closeOverlapStore(store);

  exit(0);
}
