
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

static char const *rcsid = "$Id: AS_GKP_checkPlace.C,v 1.1 2009-08-14 13:37:06 skoren Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int
Check_PlacementMesg(PlacementMesg *plc_mesg) {
  gkPlacement                gkpl;

  AS_IID                     fragIID;
  AS_IID                     frag1IID;
  AS_IID                     frag2IID;

  if (plc_mesg->action == AS_ADD)
    gkpStore->inf.plcInput++;

  if (plc_mesg->action == AS_IGNORE)
    return 0;

  //  Check that the bounding fragments are different
  //
  if (AS_UID_compare(plc_mesg->bound1, plc_mesg->bound2) == 0) {
    AS_GKP_reportError(AS_GKP_PLC_SAME_CONSTRAINT, 0,
                       AS_UID_toString(plc_mesg->frag), AS_UID_toString(plc_mesg->bound1));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  
  // Check that the bounded fragment is different from the bounding fragments
  //
  if (AS_UID_compare(plc_mesg->frag, plc_mesg->bound1) == 0 || 
         AS_UID_compare(plc_mesg->frag, plc_mesg->bound2) == 0) {
    AS_GKP_reportError(AS_GKP_PLC_SELF_CONSTRAINT, 0,
                       AS_UID_toString(plc_mesg->frag));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }  

  //  Check that the three fragments are alive and well
  //
  fragIID  = gkpStore->gkStore_getUIDtoIID(plc_mesg->frag, NULL);
  if (fragIID == 0) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DOESNT_EXIST, 0,
                       AS_UID_toString(plc_mesg->frag));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  frag1IID = gkpStore->gkStore_getUIDtoIID(plc_mesg->bound1, NULL);
  if (frag1IID == 0) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DOESNT_EXIST, 0,
                       AS_UID_toString(plc_mesg->bound1));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  frag2IID = gkpStore->gkStore_getUIDtoIID(plc_mesg->bound2, NULL);
  if (frag2IID == 0) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DOESNT_EXIST, 0,
                       AS_UID_toString(plc_mesg->bound2));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }

  //  Now grab the reads, we'll need them soon enough anyway.
  //
  //  Make sure they're not deleted
  //
  gkpStore->gkStore_getFragment(fragIID, gkFrag1, GKFRAGMENT_INF);
  if (gkFrag1->gkFragment_getIsDeleted()) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DELETED, 0,
                       AS_UID_toString(plc_mesg->frag));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  gkpStore->gkStore_getFragment(frag1IID, gkFrag1, GKFRAGMENT_INF);
  gkpStore->gkStore_getFragment(frag2IID, gkFrag2, GKFRAGMENT_INF);
  if (gkFrag1->gkFragment_getIsDeleted()) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DELETED, 0,
                       AS_UID_toString(plc_mesg->bound1));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  if (gkFrag2->gkFragment_getIsDeleted()) {
    AS_GKP_reportError(AS_GKP_PLC_FRG_DELETED, 0,
                       AS_UID_toString(plc_mesg->bound2));
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  //  Make sure the fragment is not already bounded
  //
  // if already bounded, error
  gkPlacement *oldGkpl = gkpStore->gkStore_getReadPlacement(fragIID); 
  if (oldGkpl != NULL) {
    gkpStore->gkStore_getFragment(oldGkpl->bound1, gkFrag1, GKFRAGMENT_INF);
    gkpStore->gkStore_getFragment(oldGkpl->bound2, gkFrag2, GKFRAGMENT_INF);    
    AS_GKP_reportError(AS_GKP_PLC_ALREADY_CONSTRAINED, 0,
                       AS_UID_toString(plc_mesg->frag), fragIID,
                       AS_UID_toString(gkFrag1->gkFragment_getReadUID()), gkFrag1->gkFragment_getReadIID(),
                       AS_UID_toString(gkFrag2->gkFragment_getReadUID()), gkFrag1->gkFragment_getReadIID(),
                       AS_UID_toString(plc_mesg->bound1), frag1IID,
                       AS_UID_toString(plc_mesg->bound2), frag2IID);    
    if (plc_mesg->action == AS_ADD)
      gkpStore->inf.plcErrors++;
    return(1);
  }
  
  //  And commit the change.
  if (plc_mesg->action == AS_ADD) {
    // store in gkp store
    gkpl.frag = fragIID;
    gkpl.bound1 = frag1IID;
    gkpl.bound2 = frag2IID;
    gkpStore->gkStore_addPlacement(&gkpl);
  } else {
    AS_GKP_reportError(AS_GKP_PLC_UNKNOWN_ACTION, 0);
    return 1;
  }

  return 0;
}
