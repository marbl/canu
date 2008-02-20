
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

static char const *rcsid = "$Id: AS_GKP_checkLink.c,v 1.18 2008-02-20 10:53:30 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int
Check_LinkMesg(LinkMesg *lkg_mesg) {
  AS_IID                     frag1IID;
  AS_IID                     frag2IID;
  GateKeeperFragmentRecord   gkFrag1;
  GateKeeperFragmentRecord   gkFrag2;

  if (lkg_mesg->action == AS_ADD)
    gkpStore->gkp.lkgInput++;

  if (lkg_mesg->action == AS_IGNORE)
    return 0;

  //  Check that the fragments are different
  //
  if (AS_UID_compare(lkg_mesg->frag1, lkg_mesg->frag2) == 0) {
    AS_GKP_reportError(AS_GKP_LKG_SELF_LINK,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }

  //  Check that it's a currently supported link type
  //
  if (lkg_mesg->type != AS_MATE) {
    AS_GKP_reportError(AS_GKP_LKG_UNSUPPORTED_TYPE,
                       lkg_mesg->type, AS_UID_toString1(lkg_mesg->frag1), AS_UID_toString2(lkg_mesg->frag2), AS_UID_toString3(lkg_mesg->distance));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }


  //  Check that the two fragments are alive and well
  //
  frag1IID = getGatekeeperUIDtoIID(gkpStore, lkg_mesg->frag1, NULL);
  if (frag1IID == 0) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DOESNT_EXIST,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }

  frag2IID = getGatekeeperUIDtoIID(gkpStore, lkg_mesg->frag2, NULL);
  if (frag2IID == 0) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DOESNT_EXIST,
                       AS_UID_toString(lkg_mesg->frag2));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }


  //  Now grab the reads, we'll need them soon enough anyway.
  //
  getGateKeeperFragment(gkpStore, frag1IID, &gkFrag1);
  getGateKeeperFragment(gkpStore, frag2IID, &gkFrag2);


  //  Make sure they're not deleted
  //
  if (gkFrag1.deleted) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DELETED,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }

  if (gkFrag2.deleted) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DELETED,
                       AS_UID_toString(lkg_mesg->frag2));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }

  //  Make sure they're not already mated
  //
  if (lkg_mesg->action == AS_ADD) {
    int err = 0;

    if ((gkFrag1.mateIID > 0) && (gkFrag1.mateIID != frag2IID)) {
      AS_GKP_reportError(AS_GKP_LKG_ALREADY_MATED,
                         AS_UID_toString1(gkFrag1.readUID), gkFrag1.readIID,
                         gkFrag1.mateIID,
                         AS_UID_toString2(lkg_mesg->frag2), frag2IID);
      err++;
    }
    if ((gkFrag2.mateIID > 0) && (gkFrag2.mateIID != frag1IID)) {
      AS_GKP_reportError(AS_GKP_LKG_ALREADY_MATED,
                         AS_UID_toString1(gkFrag2.readUID), gkFrag2.readIID,
                         gkFrag2.mateIID,
                         AS_UID_toString2(lkg_mesg->frag1), frag1IID);
      err++;
    }
    if (err) {
      gkpStore->gkp.lkgErrors++;
      return(1);
    }

    //  If only one of the frags has the mate relationship, and the
    //  other has mateIID == 0, well, this should be difficult to do,
    //  and your store is hosed anyway.  Something else should check
    //  the store for sanity, not gatekeeper.


    //  And if they ARE mated, just quit because their mates are
    //  correct.
    //
    if ((gkFrag1.mateIID > 0) &&
        (gkFrag2.mateIID > 0))
      return(0);
  }


  //  Version 1 encodes the library in the mate, not the read.  We
  //  need to check that the library (from a distance record) is
  //  there, and get the library IID to set in the reads.
  //
  if (AS_UID_isDefined(lkg_mesg->distance) == TRUE) {
    gkFrag1.libraryIID = getGatekeeperUIDtoIID(gkpStore, lkg_mesg->distance, NULL);
    gkFrag2.libraryIID = gkFrag1.libraryIID;

    if (gkFrag1.libraryIID == 0) {
      AS_GKP_reportError(AS_GKP_LKG_LIB_DOESNT_EXIST,
                         AS_UID_toString(lkg_mesg->distance));
      if (lkg_mesg->action == AS_ADD)
        gkpStore->gkp.lkgErrors++;
      return(1);
    }

    gkFrag1.orientation = lkg_mesg->link_orient;
    gkFrag2.orientation = lkg_mesg->link_orient;
  }

  //  Now make absolutely sure the two reads are in the same library.
  //
  if (gkFrag1.libraryIID != gkFrag2.libraryIID) {
    AS_GKP_reportError(AS_GKP_LKG_DIFFERENT_LIB,
            frag1IID, gkFrag1.libraryIID,
            frag2IID, gkFrag2.libraryIID);
    if (lkg_mesg->action == AS_ADD)
      gkpStore->gkp.lkgErrors++;
    return(1);
  }

  //  And commit the change.

  if (lkg_mesg->action == AS_ADD) {
    gkFrag1.mateIID    = frag2IID;
    gkFrag2.mateIID    = frag1IID;

    setIndexStore(gkpStore->frg, frag1IID, &gkFrag1);
    setIndexStore(gkpStore->frg, frag2IID, &gkFrag2);

    gkpStore->gkp.lkgLoaded++;
  } else if (lkg_mesg->action == AS_DELETE) {
    gkFrag1.mateIID    = 0;
    gkFrag2.mateIID    = 0;

    setIndexStore(gkpStore->frg, frag1IID, &gkFrag1);
    setIndexStore(gkpStore->frg, frag2IID, &gkFrag2);
  } else {
    AS_GKP_reportError(AS_GKP_LKG_UNKNOWN_ACTION);
    return 1;
  }

  return 0;
}
