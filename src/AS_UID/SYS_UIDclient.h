
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "SYS_UIDcommon.h"

//
//  The complicated UID client interface -- probably not what you
//  want.
//

int32    SYS_UIDgetLastUIDInterval(uint64* interval);
int32    SYS_UIDgetNewUIDInterval(uint64* interval);
int32    SYS_UIDgetMaxUIDSize(uint64* size);
void     SYS_UIDsetUIDSize(uint64 block_size);
int32    SYS_UIDgetNextUID(uint64* uid);
int32    SYS_UIDgetLastUID(uint64* uid);
void     SYS_UIDset_euid_server(const char * servers);
void     SYS_UIDset_euid_namespace(const char * namespaceName);


//
//  The simple UID client interface -- probably what you want.
//
//  Call UIDserverInitialize() to initializes the UID client.  If
//  firstUID is non-zero, the client will return numbers starting with
//  firstUID, otherwise, the UID server is contacted to get real UIDs.
//  getUID() returns one UID, either a fake or a real UID.
//

typedef struct {
  int      useDummy;
  uint64   startUID;
  uint64   interval[4];
} UIDserver;


//  firstuid is an integer (and not an AS_UID) on purpose.  The UID
//  server currently returns only integer UIDs; it makes no sense
//  to pass in a string uid.

static
UIDserver   *UIDserverInitialize(uint32 blockSize, uint64 firstuid) {
  UIDserver  *u = (UIDserver *)safe_calloc(1, sizeof(UIDserver));
  int32       s;

  if (firstuid) {
    u->useDummy    = 1;
    u->startUID    = firstuid;
  } else {
    uint64      maxBlockSize;

    // First check whether the UID server can accomodate for our buffer
    s = SYS_UIDgetMaxUIDSize(&maxBlockSize);
    if (s != UID_CODE_OK) {
      fprintf(stderr, "UIDserverInitialize()-- UID blocksize query failed (%d).\n", s);
      exit(1);
    }
    if (maxBlockSize < blockSize)
      blockSize = maxBlockSize;

    // Now set the blocksize of th UID server appropriately
    SYS_UIDsetUIDSize(blockSize);

    // Finally, get the actual interval
    s = SYS_UIDgetNewUIDInterval(u->interval);
    if (s != UID_CODE_OK) {
      fprintf(stderr, "UIDserverInitialize()-- SYS_UIDgetNewUIDInterval failed (%d) blockSize=%d startUID="F_U64", interval="F_U64","F_U64","F_U64","F_U64".\n",
              s, blockSize, u->startUID, u->interval[0], u->interval[1], u->interval[2], u->interval[3]);
      exit(1);
    }
  }

  return(u);
}

static
uint64
getUID(UIDserver *u) {
  uint64   uid = 0;
  int32    s;

  if (u->useDummy)
    return(u->startUID++);

  s = SYS_UIDgetNextUID(&uid);
  if (s != UID_CODE_OK) {
    s = SYS_UIDgetNewUIDInterval(u->interval);
    if (s != UID_CODE_OK) {
      fprintf(stderr, "getUID()-- SYS_UIDgetNewUIDInterval failed (%d) startUID="F_U64", interval="F_U64","F_U64","F_U64","F_U64".\n",
              s, u->startUID, u->interval[0], u->interval[1], u->interval[2], u->interval[3]);
      exit(1);
    }
    s = SYS_UIDgetNextUID(&uid);
  }
  if (s != UID_CODE_OK) {
    fprintf(stderr, "getUID()-- SYS_UIDgetNextUID failed (%d).\n", s);
    exit(1);
  }

  return(uid);
}


#endif






