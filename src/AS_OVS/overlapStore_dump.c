
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

static char CM_ID[] = "$Id: overlapStore_dump.c,v 1.1 2007-03-09 07:29:23 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"

void
dumpStore(char *storeName, uint32 dumpBinary, uint32 bgnIID, uint32 endIID) {

  OverlapStore  *storeFile = AS_OVS_openOverlapStore(storeName);
  OVSoverlap     overlap;

  while (AS_OVS_readOverlapFromStore(storeFile, &overlap) == TRUE) {
    if (dumpBinary) {
      AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(OVSoverlap), 1);
    } else {
      switch (overlap.dat.ovl.type) {
        case AS_OVS_TYPE_OVL:
          fprintf(stdout, "    %8d %8d %c %5d %5d %4.1f %4.1f\n",
                  overlap.a_iid,
                  overlap.b_iid,
                  overlap.dat.ovl.flipped ? 'I' : 'N',
                  overlap.dat.ovl.a_hang,
                  overlap.dat.ovl.b_hang,
                  Expand_Quality(overlap.dat.ovl.orig_erate) * 100.0,
                  Expand_Quality(overlap.dat.ovl.corr_erate) * 100.0);
          break;
        case AS_OVS_TYPE_OBT:

          //  compatible with OBT convert
          //fprintf(stdout, "%8d %8d %c %4d %4d %4d %4d %4d %4d    %5.2f\n",
                  
          fprintf(stdout, "%7d %7d  %c %4d %4d %4d  %4d %4d %4d  %5.2f\n",
                  overlap.a_iid, overlap.b_iid,
                  overlap.dat.obt.fwd ? 'f' : 'r',
                  overlap.dat.obt.a_beg,
                  overlap.dat.obt.a_end,
                  666,
                  overlap.dat.obt.b_beg,
                  overlap.dat.obt.b_end,
                  666,
                  Expand_Quality(overlap.dat.obt.erate));
          break;
        case AS_OVS_TYPE_MER:
          break;
        default:
          assert(0);
          break;
      }
    }
  }

  AS_OVS_closeOverlapStore(storeFile);
}

