
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

static char *rcsid = "$Id: AS_PER_gkStream.C,v 1.4 2009-10-26 13:20:26 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"


gkStream::gkStream(gkStore *gkp_, AS_IID beginIID_, AS_IID endIID_, uint32 flags_) {

  gkp     = gkp_;
  flags   = flags_;
  bgnIID  = 0;
  curIID  = 0;
  endIID  = 0;

  fpk = NULL;

  fnm = NULL;
  snm = NULL;
  qnm = NULL;

  fsb = NULL;
  ssb = NULL;
  qsb = NULL;

  reset(beginIID_, endIID_);
}

gkStream::~gkStream() {
  closeStream(fpk);

  closeStream(fnm);
  closeStream(snm);
  closeStream(qnm);

  closeStream(fsb);
  closeStream(ssb);
  closeStream(qsb);
}


void
gkStream::reset(AS_IID beginIID_, AS_IID endIID_) {
  uint32  stType    = 0;
  uint32  stTiid    = 0;
  uint32  edType    = 0;
  uint32  edTiid    = 0;

  uint32  firstmd   = 0;
  uint32  firstlg   = 0;

  if (beginIID_ <= 0)
    beginIID_ = 1;
  if (endIID_ <= 0)
    endIID_ = gkp->gkStore_getNumFragments();

  //  Close any open streams, we'll open them again as needed.

  closeStream(fpk);  fpk = NULL;

  closeStream(fnm);  fnm = NULL;
  closeStream(snm);  snm = NULL;
  closeStream(qnm);  qnm = NULL;

  closeStream(fsb);  fsb = NULL;
  closeStream(ssb);  ssb = NULL;
  closeStream(qsb);  qsb = NULL;

  //  Position the metadata stream -- this code is similar to
  //  gkStore_loadPartition.

  bgnIID = beginIID_;
  curIID = beginIID_ - 1;
  endIID = endIID_;

  gkp->gkStore_decodeTypeFromIID(bgnIID, stType, stTiid);
  gkp->gkStore_decodeTypeFromIID(endIID, edType, edTiid);

  if (stType == edType) {
    switch (stType) {
      case GKFRAGMENT_PACKED:
        fpk = openStream(gkp->fpk);
        resetStream(fpk, stTiid, edTiid);
        break;
      case GKFRAGMENT_NORMAL:
        fnm = openStream(gkp->fnm);
        snm = openStream(gkp->snm);
        qnm = openStream(gkp->qnm);

        firstmd = stTiid;

        resetStream(fnm, stTiid, edTiid);
        break;
      case GKFRAGMENT_STROBE:
        fsb = openStream(gkp->fsb);
        ssb = openStream(gkp->ssb);
        qsb = openStream(gkp->qsb);

        firstlg = stTiid;

        resetStream(fsb, stTiid, edTiid);
        break;
    }
  } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_NORMAL)) {
    fpk = openStream(gkp->fpk);

    fnm = openStream(gkp->fnm);
    snm = openStream(gkp->snm);
    qnm = openStream(gkp->qnm);

    resetStream(fpk, stTiid, STREAM_UNTILEND);
    resetStream(fnm, STREAM_FROMSTART, edTiid);

  } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_STROBE)) {
    fpk = openStream(gkp->fpk);

    fnm = openStream(gkp->fnm);
    snm = openStream(gkp->snm);
    qnm = openStream(gkp->qnm);

    fsb = openStream(gkp->fsb);
    ssb = openStream(gkp->ssb);
    qsb = openStream(gkp->qsb);

    resetStream(fpk, stTiid, STREAM_UNTILEND);
    resetStream(fnm, STREAM_FROMSTART, STREAM_UNTILEND);
    resetStream(fsb, STREAM_FROMSTART, edTiid);

  } else if ((stType == GKFRAGMENT_NORMAL) && (edType == GKFRAGMENT_STROBE)) {
    fnm = openStream(gkp->fnm);
    snm = openStream(gkp->snm);
    qnm = openStream(gkp->qnm);

    fsb = openStream(gkp->fsb);
    ssb = openStream(gkp->ssb);
    qsb = openStream(gkp->qsb);

    firstmd = stTiid;

    resetStream(fnm, stTiid, STREAM_UNTILEND);
    resetStream(fsb, STREAM_FROMSTART, edTiid);
  } else {
    assert(0);
  }

  //  Position the streams.  This requires us to grab a couple fragments.

  if (firstmd) {
    gkNormalFragment md;

    getIndexStore(gkp->fnm, firstmd, &md);

    resetStream(snm, md.seqOffset, STREAM_UNTILEND);
    resetStream(qnm, md.qltOffset, STREAM_UNTILEND);
  }

  if (firstlg) {
    gkStrobeFragment lg;

    getIndexStore(gkp->fsb, firstlg, &lg);

    resetStream(ssb, lg.seqOffset, STREAM_UNTILEND);
    resetStream(qsb, lg.qltOffset, STREAM_UNTILEND);
  }
}



int
gkStream::next(gkFragment *fr) {
  int           loaded = 0;

  if ((bgnIID == 0) && (curIID == 0) && (endIID == 0))
    reset(0, 0);

  if (++curIID > endIID)
    return(0);

  fr->gkp = gkp;

  //  If nextStream() on whatever stream is open returns 1, we've got
  //  a valid read.  Otherwise, move to the next type.

  if ((!loaded) && (fpk)) {
    if (nextStream(fpk, &fr->fr.packed, 0, NULL)) {
      fr->type = GKFRAGMENT_PACKED;
      fr->tiid = curIID;
      loaded = 1;
    } else {
      closeStream(fpk);
      fpk = NULL;
    }
  }

  if ((!loaded) && (fnm)) {
    if (nextStream(fnm, &fr->fr.normal, 0, NULL)) {
      fr->type = GKFRAGMENT_NORMAL;
      fr->tiid = curIID - gkp->inf.numShort;
      loaded = 1;
    } else {
      closeStream(fnm);
      fnm = NULL;
    }
  }

  if ((!loaded) && (fsb)) {
    if (nextStream(fsb, &fr->fr.strobe, 0, NULL)) {
      fr->type = GKFRAGMENT_STROBE;
      fr->tiid = curIID - gkp->inf.numShort - gkp->inf.numMedium;
      loaded = 1;
    } else {
      closeStream(fsb);
      fsb = NULL;
    }
  }

  if (!loaded) {
    assert(0);
  }

  gkp->gkStore_getFragmentData(this, fr, flags);

  return(1);
}
