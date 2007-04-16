
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

static char CM_ID[] = "$Id: AS_GKP_checkLink.c,v 1.6 2007-04-16 22:26:39 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int
Check_LinkMesg(LinkMesg *lkg_mesg,
               CDS_CID_t batchID,
               time_t currentTime,
               int verbose) {

  PHashValue_AS              value;
  CDS_IID_t                  frag1IID;
  CDS_IID_t                  frag2IID;
  GateKeeperFragmentRecord   gkFrag1;
  GateKeeperFragmentRecord   gkFrag2;


  //  Check that the fragments are different
  //
  if (lkg_mesg->frag1 == lkg_mesg->frag2) {
    fprintf(stderr,"# Link Message:  Can't make a link from fragment " F_UID " to itself!!!" ,
	    lkg_mesg->frag1);
    return(GATEKEEPER_FAILURE);
  }


  //  Check that the two fragments are alive and well
  //
  if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                UID_NAMESPACE_AS,
                                                lkg_mesg->frag1, 
                                                AS_IID_FRG, 
                                                FALSE,
                                                stderr,
                                                &value)) {
    printGKPError(stderr, GKPError_MissingFRG);
    fprintf(stderr,"# Link Message:  Fragment " F_UID " DOES NOT exist!!!" 
	    "Can't reference it... bye\n",	    lkg_mesg->frag1);
    return(GATEKEEPER_FAILURE);
  }
  frag1IID = value.IID;

  if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                UID_NAMESPACE_AS,
                                                lkg_mesg->frag2, 
                                                AS_IID_FRG, 
                                                FALSE,
                                                stderr,
                                                &value)){
    printGKPError(stderr, GKPError_MissingFRG);
    fprintf(stderr,"# Link Message:  Fragment " F_UID " DOES NOT exist!!!" 
	    "Can't reference it... bye\n",	    lkg_mesg->frag2);
    return(GATEKEEPER_FAILURE);
  }
  frag2IID = value.IID;


  //  Now grab the reads, we'll need them soon enough anyway.
  //
  getGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
  getGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);


  if (lkg_mesg->distance != 0) {
    //  Version 1 encodes the library in the mate, not the read.  We
    //  need to check that the library (from a distance record) is
    //  there, and get the library IID to set in the reads.
    //
    if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private,
                                                  UID_NAMESPACE_AS,
                                                  lkg_mesg->distance,
                                                  AS_IID_DST,
                                                  FALSE,
                                                  stderr,
                                                  &value)) {
      printGKPError(stderr, GKPError_MissingLIB);
      fprintf(stderr,"# Link Message:  Library " F_UID " DOES NOT exist!!!" 
              "Can't reference it... bye\n",	    lkg_mesg->distance);
      return(GATEKEEPER_FAILURE);
    }

    //  One could also check that the two reads are version 1...but,
    //  heck, if they aren't, we're still good to go!

    gkFrag1.libraryIID = value.IID;
    gkFrag2.libraryIID = value.IID;
  }

  //  Now make absolutely sure the two reads are in the same library.
  //
  if (gkFrag1.libraryIID != gkFrag2.libraryIID) {
    printGKPError(stderr, GKPError_LNKFragLibMismatch);
    fprintf(stderr,"# Link Message:  Fragment "F_IID" in lib "F_IID" -- fragment "F_IID" in lib "F_IID".\n",
            frag1IID, gkFrag1.libraryIID,
            frag2IID, gkFrag2.libraryIID);
    return(GATEKEEPER_FAILURE);
  }

  //  And now set the library in the reads.

  if (lkg_mesg->action == AS_ADD) {
    gkFrag1.mateIID    = frag2IID;
    gkFrag2.mateIID    = frag1IID;

    setGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
    setGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);

    AddRefPHashTable_AS(gkpStore->phs_private, UID_NAMESPACE_AS, gkFrag2.readUID);
    AddRefPHashTable_AS(gkpStore->phs_private, UID_NAMESPACE_AS, gkFrag1.readUID);
  } else if (lkg_mesg->action == AS_DELETE) {
    gkFrag1.mateIID    = 0;
    gkFrag2.mateIID    = 0;

    setGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
    setGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);

    UnRefPHashTable_AS(gkpStore->phs_private, UID_NAMESPACE_AS, gkFrag2.readUID);
    UnRefPHashTable_AS(gkpStore->phs_private, UID_NAMESPACE_AS, gkFrag1.readUID);
  } else {
    fprintf(stderr,"# LinkMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }

  return GATEKEEPER_SUCCESS;
}

