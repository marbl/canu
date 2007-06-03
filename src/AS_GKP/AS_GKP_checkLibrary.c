
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

static char CM_ID[] = "$Id: AS_GKP_checkLibrary.c,v 1.11 2007-06-03 08:13:22 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

//  Define this to allow libraries to be deleted from the store.  Not
//  recommended.
#undef ALLOW_LIBRARY_DELETE

int
Check_DistanceMesg(DistanceMesg    *dst_mesg) {
  LibraryMesg  lmesg;

  //  Upconvert to a real LibraryMesg, then pass it on to the library
  //  check.

  lmesg.action       = dst_mesg->action;
  lmesg.eaccession   = dst_mesg->eaccession;
  lmesg.mean         = dst_mesg->mean;
  lmesg.stddev       = dst_mesg->stddev;
#ifdef AS_ENABLE_SOURCE
  lmesg.source       = NULL;
#endif
  lmesg.link_orient  = AS_READ_ORIENT_INNIE;
  lmesg.num_features = 0;
  lmesg.features     = NULL;
  lmesg.values       = NULL;

  return(Check_LibraryMesg(&lmesg));
}


int
Check_LibraryMesg(LibraryMesg      *lib_mesg) {

  GateKeeperLibraryRecord  gkpl;

  if (lib_mesg->action == AS_IGNORE)
    return GATEKEEPER_SUCCESS;

  clearGateKeeperLibraryRecord(&gkpl);

  if (lib_mesg->action == AS_ADD) {
    CDS_IID_t  iid = getGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, NULL);
    if (iid) {
      fprintf(errorFP, "# LIB Error: Library "F_UID","F_IID" already exists; can't add it again.\n",
              lib_mesg->eaccession, iid);
      return(GATEKEEPER_FAILURE);
    }

    gkpl.libraryUID   = lib_mesg->eaccession;
    gkpl.comment[0]   = 0;

    switch (lib_mesg->link_orient) {
      case 'U': gkpl.orientation = AS_READ_ORIENT_UNKNOWN;    break;
      case 'I': gkpl.orientation = AS_READ_ORIENT_INNIE;      break;
      case 'O': gkpl.orientation = AS_READ_ORIENT_OUTTIE;     break;
      case 'N': gkpl.orientation = AS_READ_ORIENT_NORMAL;     break;
      case 'A': gkpl.orientation = AS_READ_ORIENT_ANTINORMAL; break;
      default:  gkpl.orientation = AS_READ_ORIENT_UNKNOWN;    break;
    }

    gkpl.deleted      = FALSE;
    gkpl.mean         = lib_mesg->mean;
    gkpl.stddev       = lib_mesg->stddev;

    if (gkpl.orientation != AS_READ_ORIENT_UNKNOWN) {
      if ((gkpl.mean   <= 0.0) ||
          (gkpl.stddev <= 0.0) ||
          (gkpl.mean - 3.0 * gkpl.stddev < 0.0)) {
        fprintf(errorFP, "# LIB Error:  Library "F_UID" has lllegal mean (%g) and standard deviation (%g).\n",
                gkpl.libraryUID, gkpl.mean, gkpl.stddev);
        return(GATEKEEPER_FAILURE);
      }
    }

#ifdef AS_ENABLE_SOURCE
    if (lib_mesg->source)
      strncpy(gkpl.comment, lib_mesg->source, AS_PER_COMMENT_LEN);
#endif

    AS_PER_decodeLibraryFeatures(&gkpl, lib_mesg);

    appendIndexStore(gkpStore->lib, &gkpl);
    setGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, getLastElemStore(gkpStore->lib), AS_IID_LIB);

  } else if (lib_mesg->action == AS_UPDATE) {
    CDS_IID_t  iid = getGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, NULL);

    if (iid == 0) {
      fprintf(errorFP, "# LIB Error:  Library "F_UID" does not exist, can't update it.\n",
	      lib_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    if ((lib_mesg->mean   <= 0.0) ||
        (lib_mesg->stddev <= 0.0) ||
        (lib_mesg->mean - 3.0 * lib_mesg->stddev < 0.0)) {
      fprintf(errorFP, "# LIB Error:  Library "F_UID" has lllegal mean (%g) and standard deviation (%g).\n",
              gkpl.libraryUID, gkpl.mean, gkpl.stddev);
      return(GATEKEEPER_FAILURE);
    }

    getIndexStore(gkpStore->lib, iid, &gkpl); 
    gkpl.mean   = lib_mesg->mean;
    gkpl.stddev = lib_mesg->stddev;
    setIndexStore(gkpStore->lib, iid, &gkpl);

  } else if (lib_mesg->action == AS_DELETE) {
    CDS_IID_t  iid = getGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, NULL);

    if (iid == 0) {
      fprintf(errorFP, "# LIB Error:  Library "F_UID" does not exist, can't delete it.\n",
              lib_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

#ifdef ALLOW_LIBRARY_DELETE
    if (HASH_SUCCESS == delGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession)) {
      GateKeeperLibraryRecord dr;
      getIndexStore(gkpStore->lib, iid, &dr);
      dr.deleted = TRUE;
      setIndexStore(gkpStore->lib, iid, &dr);
    } else {
      assert(0);
    }
#else
    fprintf(errorFP, "# LIB Error:  Library "F_UID" exists, but we don't allow libraries to be deleted.\n",
            lib_mesg->eaccession);
    return(GATEKEEPER_FAILURE);
#endif

  } else {
    fprintf(errorFP, "# LIB Error: invalid action %c.\n", lib_mesg->action);
    return GATEKEEPER_FAILURE;
  }

  return GATEKEEPER_SUCCESS;
}
