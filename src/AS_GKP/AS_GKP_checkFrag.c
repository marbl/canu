
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

static char const *rcsid = "$Id: AS_GKP_checkFrag.c,v 1.47 2009-06-10 18:05:13 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fasta.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"


static int*   isspacearray;
static int*   isValidACGTN;


int
checkSequenceAndQuality(FragMesg *frg_mesg, int *seqLen) {
  char    *s = frg_mesg->sequence;
  char    *q = frg_mesg->quality;
  char    *S = gkFrag1->gkFragment_getSequence();
  char    *Q = gkFrag1->gkFragment_getQuality();
  int      sl = 0;
  int      ql = 0;
  int      p  = 0;
  int      failed = 0;

  //  Check all letters are valid -- we change invalid ones to valid
  //  ones, just so we can load the sequence into the store.  We'll
  //  remove spaces later.
  //
  for (p = 0; s[p]; p++) {
    if ((isspacearray[s[p]]) || (isValidACGTN[s[p]])) {
    } else {
      AS_GKP_reportError(AS_GKP_FRG_INVALID_CHAR_SEQ,
                         AS_UID_toString(frg_mesg->eaccession), s[p], p);
      s[p] = 'N';
      q[p] = '0';
      failed = 1;
    }
  }

  for (p = 0; q[p]; p++) {
    if ((isspacearray[q[p]]) || ((q[p] >= '0') && (q[p] <= 'l'))) {
    } else {
      AS_GKP_reportError(AS_GKP_FRG_INVALID_CHAR_QLT,
                         AS_UID_toString(frg_mesg->eaccession), q[p], p);

      q[p] = '0';
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
      S[sl++] = isValidACGTN[s[p++]];
  }
  S[sl] = 0;

  for (p = 0; q[p]; ) {
    if (isspacearray[q[p]])
      p++;
    else
      Q[ql++] = q[p++];
  }
  Q[ql] = 0;


  //  Check that the two sequence lengths are the same.  If not,
  //  adjust to the minimum, just so we can load in the store.
  //
  if (sl != ql) {
    AS_GKP_reportError(AS_GKP_FRG_INVALID_LENGTH,
                       AS_UID_toString(frg_mesg->eaccession), sl, ql);
    sl = MIN(sl, ql);
    ql = sl;
    S[sl] = 0;
    Q[ql] = 0;
    failed = 1;
  }

  if        (sl < AS_READ_MIN_LEN) {
    AS_GKP_reportError(AS_GKP_FRG_SEQ_TOO_SHORT,
                       AS_UID_toString(frg_mesg->eaccession), sl, AS_READ_MIN_LEN);
    failed = 1;

#if 0
  } else if (sl < AS_READ_MAX_SHORT_LEN) {
    gkFrag1->gkFragment_setType(GKFRAGMENT_SHORT);
#endif

  } else if (sl < AS_READ_MAX_MEDIUM_LEN) {
    gkFrag1->gkFragment_setType(GKFRAGMENT_MEDIUM);

#if 0
  } else if (sl < AS_READ_MAX_LONG_LEN) {
    gkFrag1->gkFragment_setType(GKFRAGMENT_LONG);
#endif

  } else {
    gkFrag1->gkFragment_setType(GKFRAGMENT_LONG);

    AS_GKP_reportError(AS_GKP_FRG_SEQ_TOO_LONG,
                       AS_UID_toString(frg_mesg->eaccession), sl, AS_READ_MAX_LONG_LEN);

    sl = AS_READ_MAX_LONG_LEN - 1;

    S[sl] = 0;
    Q[sl] = 0;

    failed = 1;
  }

  gkFrag1->gkFragment_setLength(sl);

  *seqLen = sl;

  return(failed);
}




static
int
checkClearRanges(FragMesg   *frg_mesg,
                 int         seqLen,
                 int         assembler) {
  int      failed = 0;

  //  Check clear range bounds, adjust
  //
  if (frg_mesg->clear_rng.bgn < 0) {
    AS_GKP_reportError(AS_GKP_FRG_CLR_BGN,
                       AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end, seqLen);
    if (frg_mesg->action == AS_ADD)
      gkpStore->inf.frgWarnings++;
    frg_mesg->clear_rng.bgn = 0;
  }
  if (seqLen < frg_mesg->clear_rng.end) {
    AS_GKP_reportError(AS_GKP_FRG_CLR_END,
                       AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end, seqLen);
    if (frg_mesg->action == AS_ADD)
      gkpStore->inf.frgWarnings++;
    frg_mesg->clear_rng.end = seqLen;
  }

  if (frg_mesg->clear_vec.bgn > frg_mesg->clear_vec.end) {
    frg_mesg->clear_vec.bgn = 1;
    frg_mesg->clear_vec.end = 0;
  }

  if (frg_mesg->clear_max.bgn > frg_mesg->clear_max.end) {
    frg_mesg->clear_max.bgn = 1;
    frg_mesg->clear_max.end = 0;
  }

  //  If not OBT, check the clear range.
  //
  if (assembler != AS_ASSEMBLER_OBT) {
    if (frg_mesg->clear_rng.bgn > frg_mesg->clear_rng.end) {
      AS_GKP_reportError(AS_GKP_FRG_CLR_INVALID,
                         AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end);
      failed = 1;
    }

    if ((frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn) < AS_READ_MIN_LEN) {
      AS_GKP_reportError(AS_GKP_FRG_CLR_TOO_SHORT,
                         AS_UID_toString(frg_mesg->eaccession), frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn, AS_READ_MIN_LEN);
      failed = 1;
    }
  }

  return(failed);
}




static
int
setLibrary(FragMesg *frg_mesg) {
  static  uint32    libOrientationMax = 0;
  static  uint64   *libOrientation    = NULL;

  //  Version 2 comes with library information.  Get it.
  //
  //  Version 1 sets orient and library when mates are added.
  //  Version 1 libraries NEVER have a valid orientation.

  gkFrag1->gkFragment_setLibraryIID(0);
  gkFrag1->gkFragment_setMateIID(0);
  gkFrag1->gkFragment_setOrientation(0);

  if (AS_UID_isDefined(frg_mesg->library_uid) == FALSE)
    return(0);

  AS_IID  libIID = gkpStore->gkStore_getUIDtoIID(frg_mesg->library_uid, NULL);
  gkFrag1->gkFragment_setLibraryIID(libIID);

  if (AS_IID_isDefined(libIID) == FALSE) {
    AS_GKP_reportError(AS_GKP_FRG_UNKNOWN_LIB,
                       AS_UID_toString(frg_mesg->eaccession),
                       AS_UID_toString(frg_mesg->library_uid));
    return(1);
  }

  //  Get the library orientation; cache values across invocations.
  //
  if (libOrientation == NULL) {
    libOrientationMax = 2;
    libOrientation    = (uint64 *)safe_calloc(libOrientationMax, sizeof(uint64));
  }
  if (libOrientationMax <= libIID) {
    libOrientation = (uint64 *)safe_realloc(libOrientation, sizeof(uint64) * 2 * libIID);
    while (libOrientationMax < 2 * libIID)
      libOrientation[libOrientationMax++] = 0;
  }
  if (libOrientation[libIID] == 0) {
    gkLibrary  gkpl;
    gkpStore->gkStore_getLibrary(libIID, &gkpl);
    libOrientation[libIID] = gkpl.orientation;
  }

  //  Set the orientation, finally.
  //
  gkFrag1->gkFragment_setOrientation(libOrientation[libIID]);

  return(0);
}



static
int
setUID(FragMesg *frg_mesg) {

  //  Make sure we haven't seen this frag record before... if so
  //  it is a fatal error
  //
  if (AS_UID_isDefined(frg_mesg->eaccession) == FALSE) {
    AS_GKP_reportError(AS_GKP_FRG_ZERO_UID);
    return(1);
  }

  if (gkpStore->gkStore_getUIDtoIID(frg_mesg->eaccession, NULL)) {
    AS_GKP_reportError(AS_GKP_FRG_EXISTS,
                       AS_UID_toString(frg_mesg->eaccession));
    return(1);
  }

  //  NOTE!  We cannot call this until gkFrag1->type is set, done in checkSequenceAndQuality()!
  gkFrag1->gkFragment_setReadUID(frg_mesg->eaccession);

  return(0);
}



int
Check_FragMesg(FragMesg   *frg_mesg,
               int         assembler) {
  int  failed = 0;
  int  seqLen = 0;

  if (frg_mesg->action == AS_ADD)
    gkpStore->inf.frgInput++;

  if (frg_mesg->action == AS_IGNORE)
    return 0;

  if (isspacearray == NULL) {
    isspacearray    = AS_UTL_getSpaceArray();
    isValidACGTN    = AS_UTL_getValidACGTN();
  }

  if (frg_mesg->action == AS_ADD) {
    failed  = checkSequenceAndQuality(frg_mesg, &seqLen);
    failed |= checkClearRanges(frg_mesg, seqLen, assembler);
    failed |= setLibrary(frg_mesg);
    failed |= setUID(frg_mesg);

    if (failed) {
      AS_GKP_reportError(AS_GKP_FRG_LOADED_DELETED,
                         AS_UID_toString(frg_mesg->eaccession));
      gkFrag1->gkFragment_setIsDeleted(1);
    }

    gkpStore->gkStore_addFragment(gkFrag1,
                                  frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end,
                                  frg_mesg->clear_vec.bgn, frg_mesg->clear_vec.end,
                                  frg_mesg->clear_max.bgn, frg_mesg->clear_max.end,
                                  1,                       0);

  } else if (frg_mesg->action == AS_DELETE) {
    AS_IID       iid = gkpStore->gkStore_getUIDtoIID(frg_mesg->eaccession, NULL);

    if (iid == 0) {
      AS_GKP_reportError(AS_GKP_FRG_DOESNT_EXIST,
                         AS_UID_toString(frg_mesg->eaccession));
      return(1);
    }

    gkpStore->gkStore_getFragment(iid, gkFrag1, GKFRAGMENT_INF);

    if (gkFrag1->gkFragment_getMateIID() > 0) {
      AS_GKP_reportError(AS_GKP_FRG_HAS_MATE,
                         AS_UID_toString(frg_mesg->eaccession));
      return(1);
    }

    gkpStore->gkStore_delFragment(iid);

  } else {
    AS_GKP_reportError(AS_GKP_FRG_UNKNOWN_ACTION);
    return 1;
  }

  return(failed);
}
