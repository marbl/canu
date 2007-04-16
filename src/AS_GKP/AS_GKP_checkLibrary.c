
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

static char CM_ID[] = "$Id: AS_GKP_checkLibrary.c,v 1.5 2007-04-16 22:26:39 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"


int
Check_DistanceMesg(DistanceMesg    *dst_mesg,
                   CDS_CID_t        currentBatchID,
                   int              verbose) {
  LibraryMesg  lmesg;

  //  Upconvert to a real LibraryMesg, then pass it on to the library
  //  check.

  lmesg.action       = dst_mesg->action;
  lmesg.eaccession   = dst_mesg->eaccession;
  lmesg.mean         = dst_mesg->mean;
  lmesg.stddev       = dst_mesg->stddev;
  lmesg.entry_time   = 0;
#ifdef AS_ENABLE_SOURCE
  lmesg.source       = NULL;
#endif
  lmesg.link_orient  = AS_READ_ORIENT_INNIE;
  lmesg.num_features = 0;
  lmesg.features     = NULL;
  lmesg.values       = NULL;

  return(Check_LibraryMesg(&lmesg, currentBatchID, verbose));
}


int
Check_LibraryMesg(LibraryMesg      *lib_mesg,
                  CDS_CID_t         currentBatchID,
                  int               verbose) {

  GateKeeperLibraryRecord  gkpl;
  PHashValue_AS            value;

  clearGateKeeperLibraryRecord(&gkpl);

  if (lib_mesg->action == AS_ADD) {

    if (HASH_FAILURE != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                  UID_NAMESPACE_AS,
                                                  lib_mesg->eaccession, 
                                                  AS_IID_DST, 
                                                  FALSE,
                                                  stderr,
                                                  &value)) {
      printGKPError(stderr, GKPError_BadUniqueLIB);
      fprintf(stderr,"# Check_DistanceMessage:  A message with UID " F_UID " exists %s!!! Can't reuse the UID... bye\n",
	      lib_mesg->eaccession, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }

    if ((lib_mesg->mean   <= 0.0) ||
        (lib_mesg->stddev <= 0.0) ||
        (lib_mesg->mean - 3.0 * lib_mesg->stddev < 0.0)) {
      printGKPError(stderr, GKPError_DSTValues);
      fprintf(stderr,"# Check_DistanceMessage:  Illegal Mean %g and/or Standard Deviation %g\n",
	      lib_mesg->mean, lib_mesg->stddev);
      return(GATEKEEPER_FAILURE);
    }

    value.type = AS_IID_DST;
    InsertInPHashTable_AS(&gkpStore->phs_private,
                          UID_NAMESPACE_AS,
                          lib_mesg->eaccession,
                          &value,
                          FALSE,
                          TRUE);

    gkpl.libraryUID   = lib_mesg->eaccession;
    gkpl.mean         = lib_mesg->mean;
    gkpl.stddev       = lib_mesg->stddev;
    gkpl.birthBatch   = currentBatchID;

#warning not using full library message

    appendIndexStore(gkpStore->lib, &gkpl);

  } else if (lib_mesg->action == AS_REDEFINE) {

    if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                  UID_NAMESPACE_AS,
                                                  lib_mesg->eaccession, 
                                                  AS_IID_DST, 
                                                  FALSE,
                                                  stderr,
                                                  &value)) {
      printGKPError(stderr, GKPError_MissingLIB);
      fprintf(stderr,"# Check_DistanceMessage:  Distance " F_UID " does NOT exist %s!!! Can't redefine it...\n",
	      lib_mesg->eaccession, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }

    if ((lib_mesg->mean   <= 0.0) ||
        (lib_mesg->stddev <= 0.0) ||
        (lib_mesg->mean - 3.0 * lib_mesg->stddev < 0.0)) {
      printGKPError(stderr, GKPError_DSTValues);
      fprintf(stderr,"# Check_DistanceMessage:  Illegal Mean %g and/or Standard Deviation %g\n",
	      lib_mesg->mean, lib_mesg->stddev);
      return(GATEKEEPER_FAILURE);
    }

    getIndexStore(gkpStore->lib, value.IID, &gkpl); 

    gkpl.redefined      = TRUE;
    gkpl.prevID         = value.IID;
    gkpl.prevInstanceID = 0;
    gkpl.deathBatch     = currentBatchID;

    appendGateKeeperLibraryStore(gkpStore->lis, &gkpl); 

    gkpl.prevID = 0;
    gkpl.prevInstanceID = getNumGateKeeperLibrarys(gkpStore->lis) - 1;

    // Update these two fields ONLY
    gkpl.mean   = lib_mesg->mean;
    gkpl.stddev = lib_mesg->stddev;

    setGateKeeperLibraryStore(gkpStore->lib, value.IID, &gkpl); // save the IID->UID mapping

  } else if (lib_mesg->action == AS_DELETE) {

    if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
						 UID_NAMESPACE_AS,
						 lib_mesg->eaccession, 
						 AS_IID_DST, 
						 FALSE,
						 stderr,
						 &value)) {
      printGKPError(stderr, GKPError_MissingLIB);
      fprintf(stderr,"# Check_DistanceMessage:  Distance " F_UID " DOES NOT exist!!!" 
	      "Can't delete it... bye\n", lib_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    if (value.refCount > 1) {
      printGKPError(stderr, GKPError_DeleteLIB);
      fprintf(stderr, "# Check_DistanceMessage: There are %d references outstanding to Distance "F_UID ".\n"
	      "Can't delete it...\n", value.refCount, lib_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    deleteAndMarkGateKeeperLibraryStore(gkpStore->lib, value.IID, currentBatchID);

    DeleteFromPHashTable_AS(gkpStore->phs_private,UID_NAMESPACE_AS, lib_mesg->eaccession);

  } else {
    fprintf(stderr,"# Check_DistanceMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }

  return GATEKEEPER_SUCCESS;
}
