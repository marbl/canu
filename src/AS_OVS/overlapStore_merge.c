
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

static char CM_ID[] = "$Id: overlapStore_merge.c,v 1.2 2007-03-09 22:00:02 brianwalenz Exp $";

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
mergeStore(char *storeName, char *mergeName) {
  int             i;
  char            name[FILENAME_MAX];

  OverlapStore   *store = NULL;
  OverlapStore   *a     = NULL;
  OverlapStore   *b     = NULL;

  OVSoverlap      aovl;
  OVSoverlap      bovl;

  int             aValid = 1;
  int             bValid = 1;


  //  Open the two stores so we can read overlaps from them.  The
  //  store we're merging into is opened as a "copy".  The store we
  //  merge from is opened first, so that if it fails, we haven't
  //  turned the original store into a backup.
  //
  a = AS_OVS_openOverlapStorePrivate(mergeName, FALSE, FALSE);
  b = AS_OVS_openOverlapStorePrivate(storeName, TRUE,  TRUE);

  //  Recreate a store in the same place as the original store.
  //
  store = AS_OVS_createOverlapStore(storeName, FALSE);
  
  //  Now just add stuff to the new store.
  //
  aValid = AS_OVS_readOverlapFromStore(a, &aovl);
  bValid = AS_OVS_readOverlapFromStore(b, &bovl);

  while (aValid && bValid) {
    if (OVSoverlap_sort(&aovl, &bovl) < 0) {
      AS_OVS_writeOverlapToStore(store, &aovl);
      aValid = AS_OVS_readOverlapFromStore(a, &aovl);
    } else {
      AS_OVS_writeOverlapToStore(store, &bovl);
      bValid = AS_OVS_readOverlapFromStore(b, &bovl);
    }
  }

  while (aValid) {
    AS_OVS_writeOverlapToStore(store, &aovl);
    aValid = AS_OVS_readOverlapFromStore(a, &aovl);
  }

  while (bValid) {
    AS_OVS_writeOverlapToStore(store, &bovl);
    bValid = AS_OVS_readOverlapFromStore(b, &bovl);
  }

  //  ALL DONE!  Close the stores, nuke the backups and get outta here.
  //
  AS_OVS_closeOverlapStore(a);
  AS_OVS_closeOverlapStore(b);
  AS_OVS_closeOverlapStore(store);

  exit(0);
}
