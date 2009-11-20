
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

#ifndef UID_CLIENT_H
#define UID_CLIENT_H

static const char *rcsid_UID_CLIENT_H = "$Id: SYS_UIDclient.h,v 1.14 2009-11-20 22:19:06 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"

#define UID_MAX_REQUEST_SIZE 1000000

typedef struct {
  int      useDummy;

  uint64   curUID;
  uint64   maxUID;
  uint64   incUID;
} UIDserver;


uint64  getGUIDBlock(int guidRequestSize);
void    SYS_UIDset_euid_server(const char * servers);
void    SYS_UIDset_euid_namespace(const char * namespaceName);


//  Call UIDserverInitialize() to initializes the UID client.
//
//    blockSize is how many UIDs to request at once from the server.
//
//    If firstUID is non-zero, the client will return numbers starting with
//    firstUID, otherwise, the UID server is contacted to get real UIDs.
//
//  getUID() returns one UID, either a fake or a real UID.


static
UIDserver *
UIDserverInitialize(uint32 blockSize, uint64 firstuid) {
  UIDserver  *u = (UIDserver *)safe_calloc(1, sizeof(UIDserver));

  if (firstuid) {
    u->useDummy = 1;
    u->curUID   = firstuid;
    u->incUID   = 0;
    u->maxUID   = UINT64_MAX;

  } else {
    u->incUID = MIN(blockSize, UID_MAX_REQUEST_SIZE);

    //  Grab an interval

    uint64  newUID = getGUIDBlock(u->incUID);

    if (newUID > 0) {
      u->curUID = newUID;
      u->maxUID = newUID + u->incUID;
    } else {
      fprintf(stderr, "UIDserverInitialize()-- failed to get a new UID block from the server.\n");
      exit(1);
    }
  }

  return(u);
}



static
uint64
getUID(UIDserver *u) {

  if (u->useDummy)
    return(u->curUID++);

  if (u->curUID < u->maxUID)
    return(u->curUID++);

  //  Grab another interval

  uint64  newUID = getGUIDBlock(u->incUID);

  if (newUID > 0) {
    u->curUID = newUID;
    u->maxUID = newUID + u->incUID;
  } else {
    fprintf(stderr, "getUID()-- failed to get another UID block from the server.\n");
    exit(1);
  }

  if (u->curUID < u->maxUID)
    return(u->curUID++);

  //  Should never occur.
  fprintf(stderr, "getUID()-- failed to get another UID block from the server.\n");
  exit(1);

  return(0);
}

#endif
