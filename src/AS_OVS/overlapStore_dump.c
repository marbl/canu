
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

static const char *rcsid = "$Id: overlapStore_dump.c,v 1.13 2008-11-13 09:14:33 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"

void
dumpStore(char *storeName, uint32 dumpBinary, double dumpERate, uint32 bgnIID, uint32 endIID, uint32 qryIID) {

  OverlapStore  *storeFile = AS_OVS_openOverlapStore(storeName);
  OVSoverlap     overlap;

  uint64         erate = AS_OVS_encodeQuality(dumpERate / 100.0);

  AS_OVS_setRangeOverlapStore(storeFile, bgnIID, endIID);

  while (AS_OVS_readOverlapFromStore(storeFile, &overlap, AS_OVS_TYPE_ANY) == TRUE) {
    switch (overlap.dat.ovl.type) {
      case AS_OVS_TYPE_OVL:
        if ((overlap.dat.ovl.corr_erate <= erate) &&
            ((qryIID == 0) || (qryIID == overlap.b_iid)))
          if (dumpBinary)
            AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(OVSoverlap), 1);
          else
            fprintf(stdout, "%8d %8d %c %5d %5d %5.2f %5.2f %5d\n",
                    overlap.a_iid,
                    overlap.b_iid,
                    overlap.dat.ovl.flipped ? 'I' : 'N',
                    overlap.dat.ovl.a_hang,
                    overlap.dat.ovl.b_hang,
                    AS_OVS_decodeQuality(overlap.dat.ovl.orig_erate) * 100.0,
                    AS_OVS_decodeQuality(overlap.dat.ovl.corr_erate) * 100.0,
                    overlap.dat.ovl.seed_value);
        break;
      case AS_OVS_TYPE_OBT:
        if ((overlap.dat.obt.erate <= erate) &&
            ((qryIID == 0) || (qryIID == overlap.b_iid)))
          if (dumpBinary)
            AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(OVSoverlap), 1);
          else
            fprintf(stdout, "%7d %7d %c %4d %4d %4d %4d %5.2f\n",
                    overlap.a_iid,
                    overlap.b_iid,
                    overlap.dat.obt.fwd ? 'f' : 'r',
                    overlap.dat.obt.a_beg,
                    overlap.dat.obt.a_end,
                    overlap.dat.obt.b_beg,
                    overlap.dat.obt.b_end,
                    AS_OVS_decodeQuality(overlap.dat.obt.erate) * 100.0);
        break;
      case AS_OVS_TYPE_MER:
        if ((qryIID == 0) || (qryIID == overlap.b_iid))
          if (dumpBinary)
            AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(OVSoverlap), 1);
          else
            fprintf(stdout, "%7d %7d %c %4d %4d %4d %4d\n",
                    overlap.a_iid,
                    overlap.b_iid,
                    overlap.dat.mer.palindrome ? 'p' : (overlap.dat.mer.fwd ? 'f' : 'r'),
                    overlap.dat.mer.a_pos,
                    overlap.dat.mer.b_pos,
                    overlap.dat.mer.k_count,
                    overlap.dat.mer.k_len);
        break;
      default:
        assert(0);
        break;
    }
  }

  AS_OVS_closeOverlapStore(storeFile);
}

