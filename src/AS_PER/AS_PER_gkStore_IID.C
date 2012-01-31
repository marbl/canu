
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

static char *rcsid = "$Id: AS_PER_gkStore_IID.C,v 1.5 2012-01-31 07:09:21 brianwalenz Exp $";

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


//  If IIDtoTYPE and IIDtoTIID are defined, these will map from a fragment IID to a (sub)store and
//  element in that store.
//
//  In the usual case, with only one type of fragment (packed, normal or strobe) the mapping is the
//  identity, and we do not allocate these two arrays.  Even with multiple types, if they are added
//  in the correct order, we do not need these arrays.
//
//  Early versions of this new gkpStore assumed that all packed fragments were added first, then all
//  normal fragments, then all strobe fragments.  This let us use the rules in decodeTypeFromIID to
//  determine the type of a fragment given only its iid.  Hopefully, nothing still assumes this.
//  The gkStream and gkClearRange were making this assumption.
//

void
gkStore::gkStore_decodeTypeFromIID(AS_IID iid, uint32& type, uint32& tiid) {
  type = 0;
  tiid = 0;

  if (IIDtoTYPE) {
    type = IIDtoTYPE[iid];
    tiid = IIDtoTIID[iid];

  } else if (iid <= inf.numPacked) {
    type = GKFRAGMENT_PACKED;
    tiid = iid;

  } else if (iid <= inf.numPacked + inf.numNormal) {
    type = GKFRAGMENT_NORMAL;
    tiid = iid - inf.numPacked;

  } else if (iid <= inf.numPacked + inf.numNormal + inf.numStrobe) {
    type = GKFRAGMENT_STROBE;
    tiid = iid - inf.numPacked - inf.numNormal;

  } else {
    fprintf(stderr, "gkStore_decodeTypeFromIID()-- ERROR:  fragment iid %u is out of range.\n", iid);
    fprintf(stderr, "gkStore_decodeTypeFromIID()--         numPacked=%u numNormal=%u numStrobe=%u\n",
            inf.numPacked, inf.numNormal, inf.numStrobe);
    assert(0);
  }
}


void
gkStore::gkStore_addIIDtoTypeMap(AS_IID iid, uint32 type, uint32 tiid) {

  //  Bail if we've got nothing to do.

  if (((IIDtoTYPE == NULL) && (type == GKFRAGMENT_PACKED) && (inf.numNormal + inf.numStrobe == 0)) ||
      ((IIDtoTYPE == NULL) && (type == GKFRAGMENT_NORMAL) && (inf.numStrobe == 0)) ||
      ((IIDtoTYPE == NULL) && (type == GKFRAGMENT_STROBE)))
    return;

  //  Ugh.  Something to do.

  if (IIDtoTYPE == NULL) {
    gkFragment  fr;

    fprintf(stderr, "gkStore_addIIDtoTypeMap()--  Creating type map.  This is an inefficient store now.\n");

    //  In all cases, we need to scan ALL reads already in the store to create the initial map.
    //  NOTE that fr.type needs to be set for the gkFragment_ calls to work.

    IIDmax    = MIN(inf.numPacked + inf.numNormal + inf.numStrobe + (uint64)1048576, (uint64)UINT_MAX);
    IIDtoTYPE = (uint8  *)safe_malloc(sizeof(uint8)  * IIDmax);
    IIDtoTIID = (uint32 *)safe_malloc(sizeof(uint32) * IIDmax);

    memset(IIDtoTYPE, 0xff, sizeof(uint8)  * IIDmax);
    memset(IIDtoTIID, 0xff, sizeof(uint32) * IIDmax);

    fr.type = GKFRAGMENT_PACKED;
    for (uint32 i=1; i<=inf.numPacked; i++) {
      getIndexStore(fpk, i, &fr.fr.packed);
      IIDtoTYPE[fr.gkFragment_getReadIID()] = GKFRAGMENT_PACKED;
      IIDtoTIID[fr.gkFragment_getReadIID()] = i;
    }

    fr.type = GKFRAGMENT_NORMAL;
    for (uint32 i=1; i<=inf.numNormal; i++) {
      getIndexStore(fnm, i, &fr.fr.normal);
      IIDtoTYPE[fr.gkFragment_getReadIID()] = GKFRAGMENT_NORMAL;
      IIDtoTIID[fr.gkFragment_getReadIID()] = i;
    }

    fr.type = GKFRAGMENT_STROBE;
    for (uint32 i=1; i<=inf.numStrobe; i++) {
      getIndexStore(fsb, i, &fr.fr.strobe);
      IIDtoTYPE[fr.gkFragment_getReadIID()] = GKFRAGMENT_STROBE;
      IIDtoTIID[fr.gkFragment_getReadIID()] = i;
    }
  }

  //  Well, at least the task we came here to do is simple.

  if (IIDmax <= iid) {
    IIDmax = MIN(IIDmax * 2, (uint64)UINT_MAX);
    IIDtoTYPE = (uint8  *)safe_realloc(IIDtoTYPE, sizeof(uint8)  * IIDmax);
    IIDtoTIID = (uint32 *)safe_realloc(IIDtoTIID, sizeof(uint32) * IIDmax);
  }

  IIDtoTYPE[iid] = type;
  IIDtoTIID[iid] = tiid;
}



//  For a range of IIDs, this function will return the ranges of the internal stores for each type.
//  This function is straightforward if fragments are added in the correct order (packed then normal
//  then strobe); just figure out the type and id of the first and last reads.
//
//  If not conforming, we need to scan the entire range and remember the min/max ids.  Heh, which
//  turns out to be simpler than the conformaing case, but much slower.
//
void
gkStore::gkStore_computeRanges(AS_IID  bgnIID, AS_IID  endIID,
                               int64 &bgnPK,  int64 &endPK,  int64 &valPK,
                               int64 &bgnNM,  int64 &endNM,  int64 &valNM,
                               int64 &bgnSB,  int64 &endSB,  int64 &valSB) {
  uint32  stType    = 0;
  uint32  stTiid    = 0;
  uint32  edType    = 0;
  uint32  edTiid    = 0;

  bgnPK = INT64_MAX;  endPK = INT64_MAX;  valPK = 0;
  bgnNM = INT64_MAX;  endNM = INT64_MAX;  valNM = 0;
  bgnSB = INT64_MAX;  endSB = INT64_MAX;  valSB = 0;

  if (IIDtoTYPE == NULL) {
    gkStore_decodeTypeFromIID(bgnIID, stType, stTiid);
    gkStore_decodeTypeFromIID(endIID, edType, edTiid);

    if (stType == edType) {
      switch (stType) {
        case GKFRAGMENT_PACKED:
          bgnPK = stTiid;
          endPK = edTiid;
          valPK = 1;
          break;
        case GKFRAGMENT_NORMAL:
          bgnNM = stTiid;
          endNM = edTiid;
          valNM = 1;
          break;
        case GKFRAGMENT_STROBE:
          bgnSB = stTiid;
          endSB = edTiid;
          valSB = 1;
          break;
      }
    } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_NORMAL)) {
      bgnPK = stTiid;
      endPK = STREAM_UNTILEND;
      valPK = 1;

      bgnNM = STREAM_FROMSTART;
      endNM = edTiid;
      valNM = 1;
    } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_STROBE)) {
      bgnPK = stTiid;
      endPK = STREAM_UNTILEND;
      valPK = 1;
      bgnNM = STREAM_FROMSTART;
      endNM = STREAM_UNTILEND;
      valNM = 1;
      bgnSB = STREAM_FROMSTART;
      endSB  = edTiid;
      valSB = 1;
    } else if ((stType == GKFRAGMENT_NORMAL) && (edType == GKFRAGMENT_STROBE)) {
      bgnNM = stTiid;
      endNM = STREAM_UNTILEND;
      valNM = 1;
      bgnSB = STREAM_FROMSTART;
      endSB = edTiid;
      valSB = 1;
    } else {
      assert(0);
    }

    return;
  }

  //  Else, not nice.  Reads added in the wrong order.

  for (AS_IID iid=bgnIID; iid<=endIID; iid++) {
    uint32  type = IIDtoTYPE[iid];
    uint32  tiid = IIDtoTIID[iid];

    switch (type) {
      case GKFRAGMENT_PACKED:
        if (tiid < bgnPK)  bgnPK = tiid;
        if (endPK < tiid)  endPK = tiid;
        valPK = 1;
        break;
      case GKFRAGMENT_NORMAL:
        if (tiid < bgnNM)  bgnNM = tiid;
        if (endNM < tiid)  endNM = tiid;
        valNM = 1;
        break;
      case GKFRAGMENT_STROBE:
        if (tiid < bgnSB)  bgnSB = tiid;
        if (endSB < tiid)  endSB = tiid;
        valSB = 1;
        break;
      default:
        fprintf(stderr, "gkStore_computeRanges()-- invalid type "F_U32" for iid "F_U32"\n", type, iid);
        assert(0);
        break;
    }
  }
}
