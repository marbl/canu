
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

static char *rcsid = "$Id: AS_PER_gkStore_fragments.C,v 1.9 2013-01-08 02:39:15 brianwalenz Exp $";

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
gkStore::gkStore_getFragmentData(gkStream *gst, gkFragment *fr, uint32 flags) {

  assert((flags == GKFRAGMENT_INF) || (flags == GKFRAGMENT_SEQ) || (flags == GKFRAGMENT_QLT));

  fr->gkp = this;

  if (fr->enc == NULL) {
    fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
    fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
    fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
  }

  uint32  seqLen = fr->gkFragment_getSequenceLength();

  if ((fr->type == GKFRAGMENT_PACKED) &&
      ((flags == GKFRAGMENT_SEQ) ||
       (flags == GKFRAGMENT_QLT))) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (partmap)
      getIndexStore(partqpk, (int32)LookupValueInHashTable_AS(partmap, fr->fr.packed.readIID, 0), fr->enc);
    else if (gst == NULL)
      getIndexStore(qpk, fr->tiid, fr->enc);
    else
      nextStream(gst->qpk, fr->enc, 0, NULL);

    decodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
    assert(fr->seq[seqLen] == 0);
    assert(fr->qlt[seqLen] == 0);

    return;
  }

  uint32  actLen = 0;
  int64   nxtOff = 0;

  int64   seqOff = fr->gkFragment_getSequenceOffset();
  int64   qltOff = fr->gkFragment_getQualityOffset();

  if (flags == GKFRAGMENT_SEQ) {
    fr->hasSEQ = 1;

    assert(partmap == NULL);

    if (gst == NULL)
      getStringStore((fr->type == GKFRAGMENT_NORMAL) ? snm : ssb,
                     seqOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    else
      nextStream((fr->type == GKFRAGMENT_NORMAL) ? gst->snm : gst->ssb,
                 fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen);

    decodeSequence(fr->enc, fr->seq, seqLen);

    assert(fr->seq[seqLen] == 0);
  }


  if (flags == GKFRAGMENT_QLT) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (partmap)
      getStringStore((fr->type == GKFRAGMENT_NORMAL) ? partqnm : partqsb,
                     qltOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    else if (gst == NULL)
      getStringStore((fr->type == GKFRAGMENT_NORMAL) ? qnm : qsb,
                     qltOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    else
      nextStream((fr->type == GKFRAGMENT_NORMAL) ? gst->qnm : gst->qsb,
                 fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen);

    decodeSequenceQuality(fr->enc, fr->seq, fr->qlt);

    if ((fr->seq[seqLen] != 0) ||
        (fr->qlt[seqLen] != 0)) {
      fprintf(stderr, "gkStore_getFragmentData()- Potential gkpStore corruption.\n");
      fprintf(stderr, "gkStore_getFragmentData()- Fragment "F_IID" reports length %d, but seq/qlt report length "F_SIZE_T"/"F_SIZE_T".\n",
              fr->gkFragment_getReadIID(),
              seqLen,
              strlen(fr->seq),
              strlen(fr->qlt));
    }
    assert(fr->seq[seqLen] == 0);
    assert(fr->qlt[seqLen] == 0);
  }
}



void
gkStore::gkStore_getFragment(AS_IID iid, gkFragment *fr, int32 flags) {

  fr->hasSEQ = 0;
  fr->hasQLT = 0;

  fr->gkp = this;

  if (iid == 0)
    fprintf(stderr, "gkStore_getFragment()-- ERROR, attempt to retrieve non-existent fragment 0.\n");
  assert(iid > 0);

  gkStore_decodeTypeFromIID(iid, fr->type, fr->tiid);

  //fprintf(stderr, "gkStore_getFragment()--  Retrieving IID=%d from %d,%d\n", iid, fr->type, fr->tiid);

  if (partmap) {
    //  If partitioned, we have everything in memory.  This is keyed
    //  off of the global IID.

    if (ExistsInHashTable_AS(partmap, iid, 0) == FALSE)
      fprintf(stderr, "getFrag()-- ERROR!  IID "F_IID" not in partition!\n", iid);
    assert(ExistsInHashTable_AS(partmap, iid, 0) == TRUE);

    assert(fr->isGKP == 0);

    switch (fr->type) {
      case GKFRAGMENT_PACKED:
        memcpy(&fr->fr.packed, getIndexStorePtr(partfpk, (int32)LookupValueInHashTable_AS(partmap, iid, 0)), sizeof(gkPackedFragment));
        break;
      case GKFRAGMENT_NORMAL:
        memcpy(&fr->fr.normal, getIndexStorePtr(partfnm, (int32)LookupValueInHashTable_AS(partmap, iid, 0)), sizeof(gkNormalFragment));
        break;
      case GKFRAGMENT_STROBE:
        memcpy(&fr->fr.strobe, getIndexStorePtr(partfsb, (int32)LookupValueInHashTable_AS(partmap, iid, 0)), sizeof(gkStrobeFragment));
        break;
    }

  } else {
    //  Not paritioned, load from disk store.
    switch (fr->type) {
      case GKFRAGMENT_PACKED:
        getIndexStore(fpk, fr->tiid, &fr->fr.packed);
        break;
      case GKFRAGMENT_NORMAL:
        getIndexStore(fnm, fr->tiid, &fr->fr.normal);
        break;
      case GKFRAGMENT_STROBE:
        getIndexStore(fsb, fr->tiid, &fr->fr.strobe);
        break;
    }
  }

  //  GatekeeperMode assumes these are set.  Set them.  Currently only used by sffToCA, when it
  //  loads frags from the temporary store before detecting mate pairs.
  if (fr->isGKP) {
    fr->gkFragment_getClearRegion(fr->clrBgn, fr->clrEnd, AS_READ_CLEAR_CLR);
    fr->gkFragment_getClearRegion(fr->vecBgn, fr->vecEnd, AS_READ_CLEAR_VEC);
    fr->gkFragment_getClearRegion(fr->maxBgn, fr->maxEnd, AS_READ_CLEAR_MAX);
    fr->gkFragment_getClearRegion(fr->tntBgn, fr->tntEnd, AS_READ_CLEAR_TNT);
  }

  //fprintf(stderr, "gkStore_getFragment()--  Retrieved IID=%d (asked for %d) from %d,%d\n", fr->gkFragment_getReadIID(), iid, fr->type, fr->tiid);

  gkStore_getFragmentData(NULL, fr, flags);
}




void
gkStore::gkStore_setFragment(gkFragment *fr) {
  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  //  Sanity check that type and type IID agree.  These change while the store is being built.
  {
    uint32  type;
    uint32  tiid;

    gkStore_decodeTypeFromIID(fr->gkFragment_getReadIID(), type, tiid);

    assert(type == fr->type);
    assert(tiid == fr->tiid);
  }

  //fprintf(stderr, "gkStore_setFragment()--  Setting IID=%d to %d,%d\n", fr->gkFragment_getReadIID(), fr->type, fr->tiid);

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      setIndexStore(fpk, fr->tiid, &fr->fr.packed);
      break;
    case GKFRAGMENT_NORMAL:
      setIndexStore(fnm, fr->tiid, &fr->fr.normal);
      break;
    case GKFRAGMENT_STROBE:
      setIndexStore(fsb, fr->tiid, &fr->fr.strobe);
      break;
  }
}




//  Delete fragment with iid from the store.  If the fragment has a
//  mate, remove the mate relationship from both fragmentss.
//
//  If the mate is supplied, delete it too.
//
void
gkStore::gkStore_delFragment(AS_IID iid, bool deleteMateFrag) {
  gkFragment   fr;
  int32        mid = 0;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
  switch (fr.type) {
    case GKFRAGMENT_PACKED:
      mid = fr.fr.packed.mateIID;
      fr.fr.packed.deleted = 1;
      fr.fr.packed.mateIID = 0;
      break;
    case GKFRAGMENT_NORMAL:
      mid = fr.fr.normal.mateIID;
      fr.fr.normal.deleted = 1;
      fr.fr.normal.mateIID = 0;
      break;
    case GKFRAGMENT_STROBE:
      mid = fr.fr.strobe.mateIID;
      fr.fr.strobe.deleted = 1;
      fr.fr.strobe.mateIID = 0;
      break;
  }
  gkStore_setFragment(&fr);

  //  No mate, we're done.
  if (mid == 0)
    return;

  gkStore_getFragment(mid, &fr, GKFRAGMENT_INF);
  switch (fr.type) {
    case GKFRAGMENT_PACKED:
      fr.fr.packed.deleted = fr.fr.packed.deleted || deleteMateFrag;
      fr.fr.packed.mateIID = 0;
      break;
    case GKFRAGMENT_NORMAL:
      fr.fr.normal.deleted = fr.fr.normal.deleted || deleteMateFrag;
      fr.fr.normal.mateIID = 0;
      break;
    case GKFRAGMENT_STROBE:
      fr.fr.strobe.deleted = fr.fr.strobe.deleted || deleteMateFrag;
      fr.fr.strobe.mateIID = 0;
      break;
  }
  gkStore_setFragment(&fr);
}



void
gkStore::gkStore_addFragment(gkFragment *fr) {
  int encLen;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);
  assert(isCreating == 1);

  assert(fr->type != GKFRAGMENT_ERROR);

  int32 iid = gkStore_getNumFragments() + 1;

  if (fr->clrBgn <= fr->clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_enableCreate();
  if (fr->vecBgn <  fr->vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_enableCreate();
  if (fr->maxBgn <  fr->maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_enableCreate();
  if (fr->tntBgn <  fr->tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_enableCreate();

  assert(strlen(fr->gkFragment_getSequence()) == fr->gkFragment_getSequenceLength());
  assert(strlen(fr->gkFragment_getQuality())  == fr->gkFragment_getQualityLength());

  assert((fr->gkFragment_getIsDeleted() == 1) || (fr->clrBgn <= fr->clrEnd));

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      fr->tiid          = ++inf.numPacked;
      fr->fr.packed.readIID = iid;

      assert(fr->tiid == getLastElemStore(fpk) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpacePacked(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpacePacked(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpacePacked(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpacePacked(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.packed.clearBeg = fr->clrBgn;
      fr->fr.packed.clearEnd = fr->clrEnd;
      break;

    case GKFRAGMENT_NORMAL:
      fr->tiid          = ++inf.numNormal;
      fr->fr.normal.readIID = iid;

      assert(fr->tiid == getLastElemStore(fnm) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceNormal(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceNormal(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceNormal(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceNormal(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.normal.clearBeg = fr->clrBgn;
      fr->fr.normal.clearEnd = fr->clrEnd;
      break;

    case GKFRAGMENT_STROBE:
      fr->tiid           = ++inf.numStrobe;
      fr->fr.strobe.readIID = iid;

      assert(fr->tiid == getLastElemStore(fsb) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceStrobe(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceStrobe(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceStrobe(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceStrobe(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.strobe.clearBeg = fr->clrBgn;
      fr->fr.strobe.clearEnd = fr->clrEnd;
      break;
  }

  assert(fr->tiid != 0);

  //  Set clear ranges...if they're defined.  !!NOTE!!  The "CLR"
  //  range MUST be set last since it is "the" clear range, and
  //  setClearRange() always updates the latest clear range along with
  //  the named region.  (That is, setting VEC will set both VEC and
  //  LATEST.)
  //
  if (fr->vecBgn <  fr->vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_setClearRegion(fr, fr->vecBgn, fr->vecEnd);
  if (fr->maxBgn <  fr->maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_setClearRegion(fr, fr->maxBgn, fr->maxEnd);
  if (fr->tntBgn <  fr->tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_setClearRegion(fr, fr->tntBgn, fr->tntEnd);
  if (fr->clrBgn <= fr->clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_setClearRegion(fr, fr->clrBgn, fr->clrEnd);

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      assert(fr->seq[fr->fr.packed.seqLen] == 0);
      assert(fr->qlt[fr->fr.packed.seqLen] == 0);

      gkStore_setUIDtoIID(fr->fr.packed.readUID, fr->fr.packed.readIID, AS_IID_FRG);
      appendIndexStore(fpk, &fr->fr.packed);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendIndexStore(qpk, fr->enc);

      gkStore_addIIDtoTypeMap(iid, GKFRAGMENT_PACKED, fr->tiid);
      break;

    case GKFRAGMENT_NORMAL:
      assert(fr->seq[fr->fr.normal.seqLen] == 0);
      assert(fr->qlt[fr->fr.normal.seqLen] == 0);

      fr->fr.normal.seqOffset = getLastElemStore(snm) + 1;
      fr->fr.normal.qltOffset = getLastElemStore(qnm) + 1;

      gkStore_setUIDtoIID(fr->fr.normal.readUID, fr->fr.normal.readIID, AS_IID_FRG);
      appendIndexStore(fnm, &fr->fr.normal);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(snm, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qnm, fr->enc, fr->fr.normal.seqLen);

      gkStore_addIIDtoTypeMap(iid, GKFRAGMENT_NORMAL, fr->tiid);
      break;

    case GKFRAGMENT_STROBE:
      assert(fr->seq[fr->fr.strobe.seqLen] == 0);
      assert(fr->qlt[fr->fr.strobe.seqLen] == 0);

      fr->fr.strobe.seqOffset = getLastElemStore(ssb) + 1;
      fr->fr.strobe.qltOffset = getLastElemStore(qsb) + 1;

      gkStore_setUIDtoIID(fr->fr.strobe.readUID, fr->fr.strobe.readIID, AS_IID_FRG);
      appendIndexStore(fsb, &fr->fr.strobe);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(ssb, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qsb, fr->enc, fr->fr.strobe.seqLen);

      gkStore_addIIDtoTypeMap(iid, GKFRAGMENT_STROBE, fr->tiid);
      break;
  }

  //  We loaded a fragment regardless of its deleted status.  This is
  //  just the count of fragments in the store.
  inf.frgLoaded++;

  if (fr->gkFragment_getIsDeleted()) {
    //  Errors are defined as a fragment loaded but marked as deleted.
    inf.frgErrors++;
  } else {
    //  Only living fragments count towards random fragments.
    if (fr->gkFragment_getIsNonRandom() == 0)
      inf.numRandom++;
  }
}
