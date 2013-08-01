
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

static char *rcsid = "$Id: AS_PER_gkStore_load.C,v 1.4 2012-02-03 08:57:49 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.H"
#include "AS_PER_genericStore.H"
#include "AS_PER_gkpStore.H"
#include "AS_PER_encodeSequenceQuality.H"
#include "AS_UTL_fileIO.H"



void
gkStore::gkStore_load(AS_IID bgnIID, AS_IID endIID, int flags) {
  int64  bgnPK=0, endPK=0, valPK=0;
  int64  bgnNM=0, endNM=0, valNM=0;
  int64  bgnSB=0, endSB=0, valSB=0;

  //  Position the metadata stream -- this code is similar to gkStream::reset.

  assert(partmap    == NULL);
  assert(isReadOnly == 1);
  assert(isCreating == 0);

  assert(bgnIID <= gkStore_getNumFragments());
  assert(endIID <= gkStore_getNumFragments());

  if (bgnIID == 0)   bgnIID = 1;
  if (endIID == 0)   endIID = gkStore_getNumFragments();

  gkStore_computeRanges(bgnIID, endIID,
                        bgnPK, endPK, valPK,
                        bgnNM, endNM, valNM,
                        bgnSB, endSB, valSB);

  //  Load the stores.  If we're loading all the way till the end, the last+1 fragment doesn't
  //  exist.  In this case, we (ab)use the fact that convertStoreToPartialMemoryStore() treats 0 as
  //  meaning "from the start" or "till the end".

  if (valPK) {
    fpk = convertStoreToPartialMemoryStore(fpk, bgnPK, endPK);
    qpk = convertStoreToPartialMemoryStore(qpk, bgnPK, endPK);
  }

  if (valNM) {
    gkNormalFragment nmbeg;
    gkNormalFragment nmend;

    nmbeg.seqOffset = nmend.seqOffset = 0;  //  Abuse.
    nmbeg.qltOffset = nmend.qltOffset = 0;

    if (bgnNM != STREAM_FROMSTART)
      getIndexStore(fnm, bgnNM, &nmbeg);

    if ((endNM != STREAM_UNTILEND) && (endNM + 1 <= getLastElemStore(fnm)))
      getIndexStore(fnm, endNM+1, &nmend);

    fnm = convertStoreToPartialMemoryStore(fnm, bgnNM, endNM);

    if (flags == GKFRAGMENT_SEQ)
      snm = convertStoreToPartialMemoryStore(snm, nmbeg.seqOffset, nmend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qnm = convertStoreToPartialMemoryStore(qnm, nmbeg.qltOffset, nmend.qltOffset);
  }

  if (valSB) {
    gkStrobeFragment sbbeg;
    gkStrobeFragment sbend;

    sbbeg.seqOffset = sbend.seqOffset = 0;  //  Abuse.
    sbbeg.qltOffset = sbend.qltOffset = 0;

    if (bgnSB != STREAM_FROMSTART)
      getIndexStore(fsb, bgnSB, &sbbeg);

    if ((endSB != STREAM_UNTILEND) && (endSB + 1 <= getLastElemStore(fsb)))
      getIndexStore(fsb, endSB+1, &sbend);

    fsb = convertStoreToPartialMemoryStore(fsb, bgnSB, endSB);

    if (flags == GKFRAGMENT_SEQ)
      ssb = convertStoreToPartialMemoryStore(ssb, sbbeg.seqOffset, sbend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qsb = convertStoreToPartialMemoryStore(qsb, sbbeg.qltOffset, sbend.qltOffset);
  }
}
