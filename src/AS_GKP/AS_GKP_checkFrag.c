
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

static char const *rcsid = "$Id: AS_GKP_checkFrag.c,v 1.36 2008-02-20 10:53:30 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"


static int    checkfraginitialized = 0;
static int    isspacearray[256] = {0};
static int    isvalidACGTN[256] = {0};
static char   encodedsequence[AS_FRAG_MAX_LEN+1] = {0};




int
checkSequenceAndQuality(FragMesg *frg_mesg) {
  char    *s = frg_mesg->sequence;
  char    *q = frg_mesg->quality;
  int      sl = 0;
  int      ql = 0;
  int      p  = 0;
  int      failed = 0;

  //  Check all letters are valid -- we change invalid ones to valid
  //  ones, just so we can load the sequence into the store.  We'll
  //  remove spaces later.
  //
  for (p = 0; s[p]; p++) {
    if ((isspacearray[s[p]]) || (isvalidACGTN[s[p]])) {
    } else {
      AS_GKP_reportError(AS_GKP_FRG_INVALID_CHAR_SEQ,
                         AS_UID_toString(frg_mesg->eaccession), s[p], p);
      failed = 1;
    }
  }

  for (p = 0; q[p]; p++) {
    if ((isspacearray[q[p]]) || ((q[p] >= '0') && (q[p] <= 'l'))) {
    } else {
      AS_GKP_reportError(AS_GKP_FRG_INVALID_CHAR_QLT,
                         AS_UID_toString(frg_mesg->eaccession), q[p], p);
      failed = 1;
    }
  }


  //  The reader should get rid of white space, but it doesn't hurt
  //  (too much) to do it again.
  //
  for (p = 0; s[p]; ) {
    if (isspacearray[s[p]])
      p++;
    else
      s[sl++] = isvalidACGTN[s[p++]];
  }
  s[sl] = 0;

  for (p = 0; q[p]; ) {
    if (isspacearray[q[p]])
      p++;
    else
      q[ql++] = q[p++];
  }
  q[ql] = 0;


  //  Check that the two sequence lengths are the same.  If not,
  //  adjust to the minimum, just so we can load in the store.
  //
  if (sl != ql) {
    AS_GKP_reportError(AS_GKP_FRG_INVALID_LENGTH,
                       AS_UID_toString(frg_mesg->eaccession), sl, ql);
    sl = MIN(sl, ql);
    ql = sl;
    s[sl] = 0;
    q[ql] = 0;
    failed = 1;
  }

  return(failed);
}



static
int
checkClearRanges(FragMesg *frg_mesg,
                 int       assembler) {
  int      failed = 0;
  int      sl     = strlen(frg_mesg->sequence);

  //  Check that lengths are legit -- not too long, not too short.  If
  //  they are too long, we trim them, again, just so we can load them
  //  into the store.
  //
  if (sl > AS_FRAG_MAX_LEN) {
    AS_GKP_reportError(AS_GKP_FRG_SEQ_TOO_LONG,
                       AS_UID_toString(frg_mesg->eaccession), sl, AS_FRAG_MAX_LEN);
    sl = AS_FRAG_MAX_LEN;
    frg_mesg->sequence[sl] = 0;
    frg_mesg->quality[sl]  = 0;
    failed = 1;
  }

  if (sl < AS_FRAG_MIN_LEN) {
    AS_GKP_reportError(AS_GKP_FRG_SEQ_TOO_SHORT,
                       AS_UID_toString(frg_mesg->eaccession), sl, AS_FRAG_MIN_LEN);
    failed = 1;
  }


  //  Check clear range bounds, adjust
  //
  if (frg_mesg->clear_rng.bgn < 0) {
    AS_GKP_reportError(AS_GKP_FRG_CLR_BGN,
                       AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end, sl);
    if (frg_mesg->action == AS_ADD)
      gkpStore->gkp.frgWarnings++;
    frg_mesg->clear_rng.bgn = 0;
  }
  if (sl < frg_mesg->clear_rng.end) {
    AS_GKP_reportError(AS_GKP_FRG_CLR_END,
                       AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end, sl);
    if (frg_mesg->action == AS_ADD)
      gkpStore->gkp.frgWarnings++;
    frg_mesg->clear_rng.end = sl;
  }


  //  Check clear range length -- only if we are not OBT
  //
  if (assembler != AS_ASSEMBLER_OBT) {
    if ((frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn) < AS_FRAG_MIN_LEN) {
      AS_GKP_reportError(AS_GKP_FRG_CLR_TOO_SHORT,
                         AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn, AS_FRAG_MIN_LEN);
      failed = 1;
    }
  }

  return(failed);
}




static
int
setLibrary(GateKeeperFragmentRecord *gkf, FragMesg *frg_mesg) {
  static  uint32    libOrientationMax = 0;
  static  uint64   *libOrientation    = NULL;

  //  Version 2 comes with library information.  Get it.
  //
  //  Version 1 sets orient and library when mates are added.
  //  Version 1 libraries NEVER have a valid orientation.

  gkf->libraryIID  = 0;
  gkf->orientation = AS_READ_ORIENT_UNKNOWN;

  if (AS_UID_isDefined(frg_mesg->library_uid) == FALSE)
    return(0);

  gkf->libraryIID = getGatekeeperUIDtoIID(gkpStore, frg_mesg->library_uid, NULL);

  if (AS_IID_isDefined(gkf->libraryIID) == FALSE) {
    AS_GKP_reportError(AS_GKP_FRG_UNKNOWN_LIB,
                       AS_UID_toString1(frg_mesg->eaccession),
                       AS_UID_toString2(frg_mesg->library_uid));
    return(1);
  }

  //
  //  Get the library orientation; cache values across invocations.
  //

  if (libOrientation == NULL) {
    libOrientationMax = 256;
    libOrientation    = (uint64 *)safe_calloc(libOrientationMax, sizeof(uint64));
  }

  if (libOrientationMax <= gkf->libraryIID) {
    libOrientation    = safe_realloc(libOrientation, 2 * gkf->libraryIID);
    while (libOrientationMax < 2 * gkf->libraryIID)
      libOrientation[libOrientationMax++] = 0;
  }

  if (libOrientation[gkf->libraryIID] == 0) {
    GateKeeperLibraryRecord  gkpl;
    getIndexStore(gkpStore->lib, gkf->libraryIID, &gkpl); 
    libOrientation[gkf->libraryIID] = gkpl.orientation;
  }

  gkf->orientation = libOrientation[gkf->libraryIID];

  return(0);
}


static
int
setClearRanges(GateKeeperFragmentRecord *gkf, FragMesg *frg_mesg) {
  int which;

  for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
    gkf->clearBeg[which] = frg_mesg->clear_rng.bgn;
    gkf->clearEnd[which] = frg_mesg->clear_rng.end;
  }

  if (frg_mesg->clear_qlt.bgn <= frg_mesg->clear_qlt.end) {
    gkf->hasQualityClear = 1;
    gkf->clearBeg[AS_READ_CLEAR_QLT] = frg_mesg->clear_qlt.bgn;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = frg_mesg->clear_qlt.end;
  } else {
    gkf->clearBeg[AS_READ_CLEAR_QLT] = 1;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = 0;
  }

  if (frg_mesg->clear_vec.bgn <= frg_mesg->clear_vec.end) {
    gkf->hasVectorClear = 1;
    gkf->clearBeg[AS_READ_CLEAR_VEC] = frg_mesg->clear_vec.bgn;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = frg_mesg->clear_vec.end;
  } else {
    gkf->clearBeg[AS_READ_CLEAR_VEC] = 1;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = 0;
  }

  return(0);
}



int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   assembler) {

  int               failed = 0;

  if (frg_mesg->action == AS_ADD)
    gkpStore->gkp.frgInput++;

  if (frg_mesg->action == AS_IGNORE)
    return 0;

  if (checkfraginitialized == 0) {
    int i;

    for (i=0; i<256; i++)
      isspacearray[i] = isspace(i);

    isvalidACGTN['a'] = 'A';
    isvalidACGTN['c'] = 'C';
    isvalidACGTN['g'] = 'G';
    isvalidACGTN['t'] = 'T';
    isvalidACGTN['n'] = 'N';
    isvalidACGTN['A'] = 'A';
    isvalidACGTN['C'] = 'C';
    isvalidACGTN['G'] = 'G';
    isvalidACGTN['T'] = 'T';
    isvalidACGTN['N'] = 'N';

    checkfraginitialized = 1;
  }


  if (frg_mesg->action == AS_ADD) {
    GateKeeperFragmentRecord gkf = {0};
    char                    *s = NULL;

    clearGateKeeperFragmentRecord(&gkf);

    //  Make sure we haven't seen this frag record before... if so
    //  it is a fatal error
    //
    if (AS_UID_isDefined(frg_mesg->eaccession) == FALSE) {
      AS_GKP_reportError(AS_GKP_FRG_ZERO_UID);
      gkpStore->gkp.frgErrors++;
      return(1);
    }
    if (getGatekeeperUIDtoIID(gkpStore, frg_mesg->eaccession, NULL)) {
      AS_GKP_reportError(AS_GKP_FRG_EXISTS,
                         AS_UID_toString(frg_mesg->eaccession));
      gkpStore->gkp.frgErrors++;
      return(1);
    }


    //  Check sequence and quality for invalid letters, lengths, etc.
    //
    failed |= checkSequenceAndQuality(frg_mesg);
    failed |= checkClearRanges(frg_mesg, assembler);

    //  Set the library and clear ranges
    //
    failed |= setLibrary(&gkf, frg_mesg);
    failed |= setClearRanges(&gkf, frg_mesg);

    //  Report if this fragment is dead
    //
    if (failed) {
      AS_GKP_reportError(AS_GKP_FRG_LOADED_DELETED,
                         AS_UID_toString(frg_mesg->eaccession));
      gkf.deleted = TRUE;
    }

    //  Now add the fragment to the store
    //
    gkf.readUID = frg_mesg->eaccession;
    gkf.readIID = getLastElemStore(gkpStore->frg) + 1;

    gkf.seqLen = strlen(frg_mesg->sequence);
    gkf.hpsLen = 0;
    gkf.srcLen = strlen(frg_mesg->source);

    gkf.seqOffset = getLastElemStore(gkpStore->seq) + 1;
    gkf.qltOffset = getLastElemStore(gkpStore->qlt) + 1;
    gkf.hpsOffset = getLastElemStore(gkpStore->hps) + 1;
    gkf.srcOffset = getLastElemStore(gkpStore->src) + 1;

    setGatekeeperUIDtoIID(gkpStore, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkpStore->frg, &gkf);

    appendStringStore(gkpStore->seq, frg_mesg->sequence, gkf.seqLen);

    encodeSequenceQuality(encodedsequence, frg_mesg->sequence, frg_mesg->quality);
    appendStringStore(gkpStore->qlt, encodedsequence,    gkf.seqLen);

    appendStringStore(gkpStore->hps, NULL,               0);
    appendStringStore(gkpStore->src, frg_mesg->source,   gkf.srcLen);

    if (failed)
      gkpStore->gkp.frgErrors++;
    else
      gkpStore->gkp.frgLoaded++;

  } else if (frg_mesg->action == AS_DELETE) {

    GateKeeperFragmentRecord gkf;
    AS_IID                   iid = getGatekeeperUIDtoIID(gkpStore, frg_mesg->eaccession, NULL);

    if (iid == 0) {
      AS_GKP_reportError(AS_GKP_FRG_DOESNT_EXIST,
                         AS_UID_toString(frg_mesg->eaccession));
      return(1);
    }

    getGateKeeperFragment(gkpStore, iid, &gkf);

    if (gkf.mateIID > 0) {
      AS_GKP_reportError(AS_GKP_FRG_HAS_MATE,
                         AS_UID_toString(frg_mesg->eaccession));
      return(1);
    }      

    GateKeeperFragmentRecord dr;
    getIndexStore(gkpStore->frg, iid, &dr);
    dr.deleted = TRUE;
    setIndexStore(gkpStore->frg, iid, &dr);

  } else {
    AS_GKP_reportError(AS_GKP_FRG_UNKNOWN_ACTION);
    return 1;
  }

  return(failed);
}
