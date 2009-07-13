
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

static char *rcsid = "$Id: AS_PER_gkStream.C,v 1.3 2009-07-13 02:31:33 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"


gkStream::gkStream(gkStore *gkp_, AS_IID beginIID_, AS_IID endIID_, uint32 flags_) {

  gkp     = gkp_;
  flags   = flags_;
  bgnIID  = 0;
  curIID  = 0;
  endIID  = 0;

  fsm = NULL;

  fmd = NULL;
  smd = NULL;
  qmd = NULL;

  flg = NULL;
  slg = NULL;
  qlg = NULL;

  reset(beginIID_, endIID_);
}

gkStream::~gkStream() {
  closeStream(fsm);

  closeStream(fmd);
  closeStream(smd);
  closeStream(qmd);

  closeStream(flg);
  closeStream(slg);
  closeStream(qlg);
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

  closeStream(fsm);  fsm = NULL;

  closeStream(fmd);  fmd = NULL;
  closeStream(smd);  smd = NULL;
  closeStream(qmd);  qmd = NULL;

  closeStream(flg);  flg = NULL;
  closeStream(slg);  slg = NULL;
  closeStream(qlg);  qlg = NULL;

  //  Position the metadata stream -- this code is similar to
  //  gkStore_loadPartition.

  bgnIID = beginIID_;
  curIID = beginIID_ - 1;
  endIID = endIID_;

  gkp->gkStore_decodeTypeFromIID(bgnIID, stType, stTiid);
  gkp->gkStore_decodeTypeFromIID(endIID, edType, edTiid);

  if (stType == edType) {
    switch (stType) {
      case GKFRAGMENT_SHORT:
        fsm = openStream(gkp->fsm);
        resetStream(fsm, stTiid, edTiid);
        break;
      case GKFRAGMENT_MEDIUM:
        fmd = openStream(gkp->fmd);
        smd = openStream(gkp->smd);
        qmd = openStream(gkp->qmd);

        firstmd = stTiid;

        resetStream(fmd, stTiid, edTiid);
        break;
      case GKFRAGMENT_LONG:
        flg = openStream(gkp->flg);
        slg = openStream(gkp->slg);
        qlg = openStream(gkp->qlg);

        firstlg = stTiid;

        resetStream(flg, stTiid, edTiid);
        break;
    }
  } else if ((stType == GKFRAGMENT_SHORT) && (edType == GKFRAGMENT_MEDIUM)) {
    fsm = openStream(gkp->fsm);

    fmd = openStream(gkp->fmd);
    smd = openStream(gkp->smd);
    qmd = openStream(gkp->qmd);

    resetStream(fsm, stTiid, STREAM_UNTILEND);
    resetStream(fmd, STREAM_FROMSTART, edTiid);

  } else if ((stType == GKFRAGMENT_SHORT) && (edType == GKFRAGMENT_LONG)) {
    fsm = openStream(gkp->fsm);

    fmd = openStream(gkp->fmd);
    smd = openStream(gkp->smd);
    qmd = openStream(gkp->qmd);

    flg = openStream(gkp->flg);
    slg = openStream(gkp->slg);
    qlg = openStream(gkp->qlg);

    resetStream(fsm, stTiid, STREAM_UNTILEND);
    resetStream(fmd, STREAM_FROMSTART, STREAM_UNTILEND);
    resetStream(flg, STREAM_FROMSTART, edTiid);

  } else if ((stType == GKFRAGMENT_MEDIUM) && (edType == GKFRAGMENT_LONG)) {
    fmd = openStream(gkp->fmd);
    smd = openStream(gkp->smd);
    qmd = openStream(gkp->qmd);

    flg = openStream(gkp->flg);
    slg = openStream(gkp->slg);
    qlg = openStream(gkp->qlg);

    firstmd = stTiid;

    resetStream(fmd, stTiid, STREAM_UNTILEND);
    resetStream(flg, STREAM_FROMSTART, edTiid);
  } else {
    assert(0);
  }

  //  Position the streams.  This requires us to grab a couple fragments.

  if (firstmd) {
    gkMediumFragment md;

    getIndexStore(gkp->fmd, firstmd, &md);

    resetStream(smd, md.seqOffset, STREAM_UNTILEND);
    resetStream(qmd, md.qltOffset, STREAM_UNTILEND);
  }

  if (firstlg) {
    gkLongFragment lg;

    getIndexStore(gkp->flg, firstlg, &lg);

    resetStream(slg, lg.seqOffset, STREAM_UNTILEND);
    resetStream(qlg, lg.qltOffset, STREAM_UNTILEND);
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

  if ((!loaded) && (fsm)) {
    if (nextStream(fsm, &fr->fr.sm, 0, NULL)) {
      fr->type = GKFRAGMENT_SHORT;
      fr->tiid = curIID;
      loaded = 1;
    } else {
      closeStream(fsm);
      fsm = NULL;
    }
  }

  if ((!loaded) && (fmd)) {
    if (nextStream(fmd, &fr->fr.md, 0, NULL)) {
      fr->type = GKFRAGMENT_MEDIUM;
      fr->tiid = curIID - gkp->inf.numShort;
      loaded = 1;
    } else {
      closeStream(fmd);
      fmd = NULL;
    }
  }

  if ((!loaded) && (flg)) {
    if (nextStream(flg, &fr->fr.lg, 0, NULL)) {
      fr->type = GKFRAGMENT_LONG;
      fr->tiid = curIID - gkp->inf.numShort - gkp->inf.numMedium;
      loaded = 1;
    } else {
      closeStream(flg);
      flg = NULL;
    }
  }

  if (!loaded) {
    assert(0);
  }

  gkp->gkStore_getFragmentData(this, fr, flags);

  return(1);
}
