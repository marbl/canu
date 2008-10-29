
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

static char const *rcsid = "$Id: AS_GKP_checkLibrary.c,v 1.28 2008-10-29 10:53:25 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int
Check_DistanceMesg(DistanceMesg    *dst_mesg,
                   int              fixInsertSizes) {
  LibraryMesg  lmesg;

  //  Upconvert to a real LibraryMesg, then pass it on to the library
  //  check.

  lmesg.action       = dst_mesg->action;
  lmesg.eaccession   = dst_mesg->eaccession;
  lmesg.mean         = dst_mesg->mean;
  lmesg.stddev       = dst_mesg->stddev;
  lmesg.source       = NULL;
  lmesg.link_orient  = AS_READ_ORIENT_INNIE;
  lmesg.num_features = 0;
  lmesg.features     = NULL;
  lmesg.values       = NULL;

  return(Check_LibraryMesg(&lmesg, fixInsertSizes));
}



void
checkLibraryDistances(LibraryMesg *lib_mesg,
                      int          fixInsertSizes) {

  if (lib_mesg->link_orient == 'U')
    return;


  if ((lib_mesg->mean   <= 0.0) &&
      (lib_mesg->stddev <= 0.0)) {
    AS_GKP_reportError(AS_GKP_LIB_ILLEGAL_MEAN_STDDEV,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, lib_mesg->stddev);
    if (lib_mesg->action == AS_ADD)
      gkpStore->gkp.libWarnings++;
    lib_mesg->mean   = 3000.0;
    lib_mesg->stddev = 300.0;
  }

  if (lib_mesg->mean   <= 0.0) {
    AS_GKP_reportError(AS_GKP_LIB_INVALID_MEAN,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, 10.0 * lib_mesg->stddev);
    if (lib_mesg->action == AS_ADD)
      gkpStore->gkp.libWarnings++;
    lib_mesg->mean = 10.0 * lib_mesg->stddev;
  }

  if (lib_mesg->stddev <= 0.0) {
    AS_GKP_reportError(AS_GKP_LIB_INVALID_STDDEV,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->stddev, 0.1 * lib_mesg->mean);
    if (lib_mesg->action == AS_ADD)
      gkpStore->gkp.libWarnings++;
    lib_mesg->stddev = 0.1 * lib_mesg->mean;
  }

  if (fixInsertSizes) {
    if (lib_mesg->mean < 3.0 * lib_mesg->stddev) {
      AS_GKP_reportError(AS_GKP_LIB_STDDEV_TOO_BIG,
                         AS_UID_toString(lib_mesg->eaccession), lib_mesg->stddev, lib_mesg->mean, 0.1 * lib_mesg->mean);
      if (lib_mesg->action == AS_ADD)
        gkpStore->gkp.libWarnings++;
      lib_mesg->stddev = 0.1 * lib_mesg->mean;
    }

    //  What's the 0.001 for?  If we reset the stddev in any of the
    //  blocks above, we can still fail the test below, because
    //  floating point math sucks.

    if (lib_mesg->stddev + 0.001 < 0.1 * lib_mesg->mean) {
      AS_GKP_reportError(AS_GKP_LIB_STDDEV_TOO_SMALL,
                         AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, lib_mesg->stddev, 0.1 * lib_mesg->mean);
      if (lib_mesg->action == AS_ADD)
        gkpStore->gkp.libWarnings++;
      lib_mesg->stddev = 0.1 * lib_mesg->mean;
    }
  }
}


int
Check_LibraryMesg(LibraryMesg      *lib_mesg,
                  int                fixInsertSizes) {

  GateKeeperLibraryRecord  gkpl;

  if (lib_mesg->action == AS_ADD)
    gkpStore->gkp.libInput++;

  if (lib_mesg->action == AS_IGNORE)
    return 0;

  clearGateKeeperLibraryRecord(&gkpl);

  checkLibraryDistances(lib_mesg, fixInsertSizes);

  if (lib_mesg->action == AS_ADD) {
    AS_IID     iid = getGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, NULL);
    if (iid) {
      AS_GKP_reportError(AS_GKP_LIB_EXISTS,
                         AS_UID_toString(lib_mesg->eaccession), iid);
      gkpStore->gkp.libErrors++;
      return(1);
    }
    if (AS_UID_isDefined(lib_mesg->eaccession) == FALSE) {
      AS_GKP_reportError(AS_GKP_LIB_ZERO_UID);
      gkpStore->gkp.libErrors++;
      return(1);
    }

    gkpl.libraryUID   = lib_mesg->eaccession;
    gkpl.comment[0]   = 0;

    //  This crud is documented in AS_PER/AS_PER_gkpStore.h
    //  Zero is the default, we set to make it explicit

    gkpl.spare2                     = 0;
    gkpl.spare1                     = 0;

    gkpl.forceBOGunitigger          = 0;

    gkpl.doNotQVTrim                = 0;
    gkpl.goodBadQVThreshold         = 12;

    gkpl.unused1                    = 0;
    gkpl.doNotTrustHomopolymerRuns  = 0;
    gkpl.doNotOverlapTrim           = 0;
    gkpl.isNotRandom                = 0;

    gkpl.hpsIsSomethingElse         = 0;
    gkpl.hpsIsFlowGram              = 0;
    gkpl.hpsIsPeakSpacing           = 0;

    gkpl.orientation                = AS_READ_ORIENT_UNKNOWN;

    switch (lib_mesg->link_orient) {
      case 'U': gkpl.orientation = AS_READ_ORIENT_UNKNOWN;    break;
      case 'I': gkpl.orientation = AS_READ_ORIENT_INNIE;      break;
      case 'O': gkpl.orientation = AS_READ_ORIENT_OUTTIE;     break;
      case 'N': gkpl.orientation = AS_READ_ORIENT_NORMAL;     break;
      case 'A': gkpl.orientation = AS_READ_ORIENT_ANTINORMAL; break;
      default:  gkpl.orientation = AS_READ_ORIENT_UNKNOWN;    break;
    }

    gkpl.mean                       = lib_mesg->mean;
    gkpl.stddev                     = lib_mesg->stddev;

    if (lib_mesg->source)
      strncpy(gkpl.comment, lib_mesg->source, AS_PER_COMMENT_LEN);

    AS_PER_decodeLibraryFeatures(&gkpl, lib_mesg);

    appendIndexStore(gkpStore->lib, &gkpl);
    setGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, getLastElemStore(gkpStore->lib), AS_IID_LIB);

    gkpStore->gkp.libLoaded++;

  } else if (lib_mesg->action == AS_UPDATE) {
    AS_IID     iid = getGatekeeperUIDtoIID(gkpStore, lib_mesg->eaccession, NULL);

    if (iid == 0) {
      AS_GKP_reportError(AS_GKP_LIB_DOESNT_EXIST_UPDATE,
                         AS_UID_toString(lib_mesg->eaccession));
      return(1);
    }

    getIndexStore(gkpStore->lib, iid, &gkpl);

    if ((gkpl.mean   != lib_mesg->mean) ||
        (gkpl.stddev != lib_mesg->stddev)) {
      gkpl.mean   = lib_mesg->mean;
      gkpl.stddev = lib_mesg->stddev;

      setIndexStore(gkpStore->lib, iid, &gkpl);
    }

  } else {
    AS_GKP_reportError(AS_GKP_LIB_UNKNOWN_ACTION);
    return 1;
  }

  return 0;
}
