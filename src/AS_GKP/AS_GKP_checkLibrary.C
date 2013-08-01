
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

static char const *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.H"
#include "AS_GKP_include.H"
#include "AS_PER_gkpStore.H"

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
  lmesg.link_orient.setIsInnie();
  lmesg.num_features = 0;
  lmesg.features     = NULL;
  lmesg.values       = NULL;

  return(Check_LibraryMesg(&lmesg, fixInsertSizes, 0));
}



void
checkLibraryDistances(LibraryMesg *lib_mesg,
                      int          fixInsertSizes) {

  if (lib_mesg->link_orient.isUnknown())
    return;


  if ((lib_mesg->mean   <= 0.0) &&
      (lib_mesg->stddev <= 0.0)) {
    AS_GKP_reportError(AS_GKP_LIB_ILLEGAL_MEAN_STDDEV, 0,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, lib_mesg->stddev);
    if (lib_mesg->action == AS_ADD)
      gkpStore->inf.libWarnings++;
    lib_mesg->mean   = 3000.0;
    lib_mesg->stddev = 300.0;
  }

  if (lib_mesg->mean   <= 0.0) {
    AS_GKP_reportError(AS_GKP_LIB_INVALID_MEAN, 0,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, 10.0 * lib_mesg->stddev);
    if (lib_mesg->action == AS_ADD)
      gkpStore->inf.libWarnings++;
    lib_mesg->mean = 10.0 * lib_mesg->stddev;
  }

  if (lib_mesg->stddev <= 0.0) {
    AS_GKP_reportError(AS_GKP_LIB_INVALID_STDDEV, 0,
                       AS_UID_toString(lib_mesg->eaccession), lib_mesg->stddev, 0.1 * lib_mesg->mean);
    if (lib_mesg->action == AS_ADD)
      gkpStore->inf.libWarnings++;
    lib_mesg->stddev = 0.1 * lib_mesg->mean;
  }

  if (fixInsertSizes) {
    if (lib_mesg->mean < 3.0 * lib_mesg->stddev) {
      AS_GKP_reportError(AS_GKP_LIB_STDDEV_TOO_BIG, 0,
                         AS_UID_toString(lib_mesg->eaccession), lib_mesg->stddev, lib_mesg->mean, 0.1 * lib_mesg->mean);
      if (lib_mesg->action == AS_ADD)
        gkpStore->inf.libWarnings++;
      lib_mesg->stddev = 0.1 * lib_mesg->mean;
    }

    //  What's the 0.001 for?  If we reset the stddev in any of the
    //  blocks above, we can still fail the test below, because
    //  floating point math sucks.

    if (lib_mesg->stddev + 0.001 < 0.1 * lib_mesg->mean) {
      AS_GKP_reportError(AS_GKP_LIB_STDDEV_TOO_SMALL, 0,
                         AS_UID_toString(lib_mesg->eaccession), lib_mesg->mean, lib_mesg->stddev, 0.1 * lib_mesg->mean);
      if (lib_mesg->action == AS_ADD)
        gkpStore->inf.libWarnings++;
      lib_mesg->stddev = 0.1 * lib_mesg->mean;
    }
  }
}


int
Check_LibraryMesg(LibraryMesg      *lib_mesg,
                  int                fixInsertSizes,
                  uint32             packedLength) {

  gkLibrary  gkpl;

  if (lib_mesg->action == AS_ADD)
    gkpStore->inf.libInput++;

  if (lib_mesg->action == AS_IGNORE)
    return 0;

  checkLibraryDistances(lib_mesg, fixInsertSizes);

  if (lib_mesg->action == AS_ADD) {
    AS_IID     iid = gkpStore->gkStore_getUIDtoIID(lib_mesg->eaccession, NULL);
    if (iid) {
      AS_GKP_reportError(AS_GKP_LIB_EXISTS, 0,
                         AS_UID_toString(lib_mesg->eaccession), iid);
      gkpStore->inf.libErrors++;
      checkLibraryForFastQPointers(lib_mesg, iid, packedLength);
      return(1);
    }
    if (AS_UID_isDefined(lib_mesg->eaccession) == FALSE) {
      AS_GKP_reportError(AS_GKP_LIB_ZERO_UID, 0);
      gkpStore->inf.libErrors++;
      return(1);
    }

    gkpl.libraryUID                 = lib_mesg->eaccession;

    strcpy(gkpl.libraryName, AS_UID_toString(gkpl.libraryUID));

    gkpl.mean                       = lib_mesg->mean;
    gkpl.stddev                     = lib_mesg->stddev;

    gkpl.orientation                = AS_READ_ORIENT_UNKNOWN;

    if (lib_mesg->link_orient.isInnie())
      gkpl.orientation = AS_READ_ORIENT_INNIE;

    if (lib_mesg->link_orient.isOuttie())
      gkpl.orientation = AS_READ_ORIENT_OUTTIE;

    if (lib_mesg->link_orient.isNormal())
      gkpl.orientation = AS_READ_ORIENT_NORMAL;

    if (lib_mesg->link_orient.isAnti())
      gkpl.orientation = AS_READ_ORIENT_ANTINORMAL;

    gkpl.gkLibrary_decodeFeatures(lib_mesg);

    gkpStore->gkStore_addLibrary(lib_mesg->eaccession, &gkpl);

    //  If this library specifies fastq reads, load them now.

    iid = gkpStore->gkStore_getUIDtoIID(lib_mesg->eaccession, NULL);
    assert(iid > 0);

    checkLibraryForFastQPointers(lib_mesg,
                                 iid,
                                 packedLength);

  } else if (lib_mesg->action == AS_UPDATE) {
    AS_IID     iid = gkpStore->gkStore_getUIDtoIID(lib_mesg->eaccession, NULL);

    if (iid == 0) {
      AS_GKP_reportError(AS_GKP_LIB_DOESNT_EXIST_UPDATE, 0,
                         AS_UID_toString(lib_mesg->eaccession));
      return(1);
    }

    gkpStore->gkStore_getLibrary(iid, &gkpl);

    if ((gkpl.mean   != lib_mesg->mean) ||
        (gkpl.stddev != lib_mesg->stddev)) {
      gkpl.mean   = lib_mesg->mean;
      gkpl.stddev = lib_mesg->stddev;

      gkpStore->gkStore_setLibrary(iid, &gkpl);
    }

  } else {
    AS_GKP_reportError(AS_GKP_LIB_UNKNOWN_ACTION, 0);
    return 1;
  }

  return 0;
}
