
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

static char CM_ID[] = "$Id: AS_GKP_checkLink.c,v 1.2 2007-02-18 14:04:49 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <ctype.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

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
  if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs, 
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

  if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs, 
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


  //  Version 1 encodes the library in the mate, not the read.  We
  //  need to check that the library (from a distance record) is
  //  there, and get the library IID to set in the reads.
  //
  if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs, 
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

  //  And now set, not check, the library in the reads.
  //
  getGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
  getGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);

  gkFrag1.libraryIID = value.IID;
  gkFrag2.libraryIID = value.IID;

#if 0
  //  Keep this one around for version 2....
  if (gkFrag1.libraryIID != gkFrag2.libraryIID) {
    printGKPError(stderr, GKPError_LNKFragLibMismatch);
    fprintf(stderr,"# Link Message:  Fragments " F_UID " and "F_UID" are in different libraries!" 
	    "Can't reference it... bye\n",
	    frag1IID, frag2IID);
    return(GATEKEEPER_FAILURE);
  }
#endif

  if (lkg_mesg->action == AS_ADD) {
    gkFrag1.mateIID = frag2IID;
    gkFrag2.mateIID = frag1IID;

    setGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
    setGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);

    AddRefPHashTable_AS(gkpStore->phs, UID_NAMESPACE_AS, gkFrag2.readUID);
    AddRefPHashTable_AS(gkpStore->phs, UID_NAMESPACE_AS, gkFrag1.readUID);
  } else if (lkg_mesg->action == AS_DELETE) {
    gkFrag1.mateIID = 0;
    gkFrag2.mateIID = 0;

    setGateKeeperFragmentStore(gkpStore->frg, frag1IID, &gkFrag1);
    setGateKeeperFragmentStore(gkpStore->frg, frag2IID, &gkFrag2);

    UnRefPHashTable_AS(gkpStore->phs, UID_NAMESPACE_AS, gkFrag2.readUID);
    UnRefPHashTable_AS(gkpStore->phs, UID_NAMESPACE_AS, gkFrag1.readUID);
  } else {
    fprintf(stderr,"# LinkMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }

  return GATEKEEPER_SUCCESS;
}

