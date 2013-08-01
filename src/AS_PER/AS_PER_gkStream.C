
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

#include "AS_PER_gkpStore.H"


gkStream::gkStream(gkStore *gkp_, AS_IID beginIID_, AS_IID endIID_, uint32 flags_) {

  gkp     = gkp_;
  flags   = flags_;
  bgnIID  = 0;
  curIID  = 0;
  endIID  = 0;

  fpk = NULL;
  qpk = NULL;

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
  closeStream(qpk);

  closeStream(fnm);
  closeStream(snm);
  closeStream(qnm);

  closeStream(fsb);
  closeStream(ssb);
  closeStream(qsb);
}


void
gkStream::reset(AS_IID bgnIID_, AS_IID endIID_) {
  int64  bgnPK=0, endPK=0, valPK=0;
  int64  bgnNM=0, endNM=0, valNM=0;
  int64  bgnSB=0, endSB=0, valSB=0;

  if (bgnIID_ <= 0)
    bgnIID_ = 1;
  if (endIID_ <= 0)
    endIID_ = gkp->gkStore_getNumFragments();

  //  Close any open streams, we'll open them again as needed.

  closeStream(fpk);  fpk = NULL;
  closeStream(qpk);  qpk = NULL;

  closeStream(fnm);  fnm = NULL;
  closeStream(snm);  snm = NULL;
  closeStream(qnm);  qnm = NULL;

  closeStream(fsb);  fsb = NULL;
  closeStream(ssb);  ssb = NULL;
  closeStream(qsb);  qsb = NULL;

  //  Position the metadata stream -- this code is similar to gkStore_load.

  bgnIID = bgnIID_;
  curIID = bgnIID_ - 1;
  endIID = endIID_;

  gkp->gkStore_computeRanges(bgnIID, endIID,
                             bgnPK, endPK, valPK,
                             bgnNM, endNM, valNM,
                             bgnSB, endSB, valSB);

  if (valPK) {
    fpk = openStream(gkp->fpk);
    qpk = openStream(gkp->qpk);

    resetStream(fpk, bgnPK, endPK);
    resetStream(qpk, bgnPK, endPK);
  }

  if (valNM) {
    gkNormalFragment nm;

    fnm = openStream(gkp->fnm);
    snm = openStream(gkp->snm);
    qnm = openStream(gkp->qnm);

    resetStream(fnm, bgnNM, endNM);

    getIndexStore(gkp->fnm, bgnNM, &nm);

    resetStream(snm, nm.seqOffset, STREAM_UNTILEND);
    resetStream(qnm, nm.qltOffset, STREAM_UNTILEND);
  }

  if (valSB) {
    gkStrobeFragment sb;

    fsb = openStream(gkp->fsb);
    resetStream(fsb, bgnSB, endSB);

    getIndexStore(gkp->fsb, bgnSB, &sb);

    resetStream(ssb, sb.seqOffset, STREAM_UNTILEND);
    resetStream(qsb, sb.qltOffset, STREAM_UNTILEND);
  }
}



int
gkStream::next(gkFragment *fr) {
  int           loaded = 0;
  uint32        type   = 0;
  uint32        tiid   = 0;

  if ((bgnIID == 0) && (curIID == 0) && (endIID == 0))
    reset(0, 0);

  if (++curIID > endIID)
    return(0);

  fr->gkp = gkp;

  if (gkp->IIDtoTYPE == NULL) {

    //  A conforming store.  Reads were loaded in the correct order, and we can just return the
    //  first store that gives us a read.

    if ((!loaded) && (fpk)) {
      if (nextStream(fpk, &fr->fr.packed, 0, NULL)) {
        type   = GKFRAGMENT_PACKED;
        tiid   = curIID;
        loaded = 1;
      } else {
        closeStream(fpk);
        fpk = NULL;
      }
    }

    if ((!loaded) && (fnm)) {
      if (nextStream(fnm, &fr->fr.normal, 0, NULL)) {
        type   = GKFRAGMENT_NORMAL;
        tiid   = curIID - gkp->inf.numPacked;
        loaded = 1;
      } else {
        closeStream(fnm);
        fnm = NULL;
      }
    }

    if ((!loaded) && (fsb)) {
      if (nextStream(fsb, &fr->fr.strobe, 0, NULL)) {
        type   = GKFRAGMENT_STROBE;
        tiid   = curIID - gkp->inf.numPacked - gkp->inf.numNormal;
        loaded = 1;
      } else {
        closeStream(fsb);
        fsb = NULL;
      }
    }

    if (!loaded) {
      fprintf(stderr, "gkStream::next()-- Failed to load the next (conforming) fragment.\n");
      fprintf(stderr, "                   bgnIID=%d curIID=%d endIID=%d\n", bgnIID, curIID, endIID);
      fprintf(stderr, "                   numPacked=%d numNormal=%d numStrobe=%d\n", gkp->inf.numPacked, gkp->inf.numNormal, gkp->inf.numStrobe);
    }
    assert(loaded);

  } else {

    //  Otherwise, not a conforming store.  Use the type map to decide which store to read from.

    type = gkp->IIDtoTYPE[curIID];
    tiid = gkp->IIDtoTIID[curIID];

    switch (type) {
      case GKFRAGMENT_PACKED:
        loaded = nextStream(fpk, &fr->fr.packed, 0, NULL);
        break;
      case GKFRAGMENT_NORMAL:
        loaded = nextStream(fnm, &fr->fr.normal, 0, NULL);
        break;
      case GKFRAGMENT_STROBE:
        loaded = nextStream(fsb, &fr->fr.strobe, 0, NULL);
        break;
      default:
        fprintf(stderr, "gkStream::next()-- invalid type "F_U32" for iid "F_U32"\n", type, curIID);
        assert(0);
        break;
    }

    if (!loaded) {
      fprintf(stderr, "gkStream::next()-- Failed to load the next (nonconforming) fragment.\n");
      fprintf(stderr, "                   bgnIID=%d curIID=%d endIID=%d\n", bgnIID, curIID, endIID);
      fprintf(stderr, "                   numPacked=%d numNormal=%d numStrobe=%d\n", gkp->inf.numPacked, gkp->inf.numNormal, gkp->inf.numStrobe);
    }
    assert(loaded);
  }


  fr->type = type;
  fr->tiid = tiid;

  gkp->gkStore_getFragmentData(this, fr, flags);

  return(1);
}
