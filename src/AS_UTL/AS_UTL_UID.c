
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

// $Id: AS_UTL_UID.c,v 1.4 2008-02-20 10:51:09 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#define MAX_UID_LENGTH (128)


static  GateKeeperStore *AS_UID_gkp = NULL;

static  size_t           AS_UID_stringsLen = 0;
static  size_t           AS_UID_stringsMax = 0;
static  char            *AS_UID_strings    = NULL;

static  HashTable_AS    *AS_UID_STRtoUID   = NULL;

void
AS_UID_setGatekeeper(void *gkp) {
  AS_UID_gkp = gkp;
}

//  If we don't have a gkpStore we are likely processing a protoIO
//  file, like the asm file.  To do this, we need to convert external
//  UIDs (strings) to internal UIDs (offsets).  We duplicate most of
//  the functionality, with memory leaks.
//
//  CODE DUPLICATION:
//
//  AS_GKP_getUIDfromString() == AS_UID_getUIDfromString()
//  AS_GKP_getUID()           == AS_UID_getUID()
//  AS_GKP_addUID()           == AS_UID_addUID()
//  

AS_UID
AS_UID_getUIDfromString(AS_UID uid) {
  uint64  loc = 0;

  
  if ((AS_UID_STRtoUID != NULL) &&
      (LookupInHashTable_AS(AS_UID_STRtoUID,
                            (INTPTR)uid.UIDstring, strlen(uid.UIDstring),
                            &loc, 0))) {
    uid.isString  = 1;
    uid.UID       = loc;
    uid.UIDstring = AS_UID_strings + uid.UID;
  } else {
    uid = AS_UID_undefined();
  }
  return(uid);
}

AS_UID
AS_UID_getUID(AS_UID uid) {
  uid.UIDstring = NULL;
  if (uid.isString) {
    assert(AS_UID_strings != NULL);
    uid.UIDstring = AS_UID_strings + uid.UID;
  }
  return(uid);
}

AS_UID
AS_UID_addUID(AS_UID uid) {

  assert((uid.isString == 0) && (uid.UID == 0) && (uid.UIDstring != NULL));

  uint64     loc    = 0;
  uint64     len    = strlen(uid.UIDstring);

  if (AS_UID_STRtoUID == NULL) {
    AS_UID_STRtoUID = CreateStringHashTable_AS(32 * 1024);
    AS_UID_stringsLen = 0;
    AS_UID_stringsMax = 16 * 1024 * 1024;
    AS_UID_strings    = (char *)safe_malloc(sizeof(char) * AS_UID_stringsMax);
  }

  if (LookupInHashTable_AS(AS_UID_STRtoUID, (INTPTR)uid.UIDstring, len, &loc, 0) == FALSE) {
    loc = AS_UID_stringsLen;

    if (AS_UID_stringsLen + len + 1 >= AS_UID_stringsMax) {
      AS_UID_stringsMax *= 2;
      AS_UID_strings     = (char *)safe_realloc(AS_UID_strings, sizeof(char) * AS_UID_stringsMax);
    }

    memcpy(AS_UID_strings + AS_UID_stringsLen, uid.UIDstring, len + 1);
    AS_UID_stringsLen += len + 1;

    if (InsertInHashTable_AS(AS_UID_STRtoUID,
                             (INTPTR)AS_UID_strings + loc, len,
                             loc, 0) == HASH_FAILURE) {
      fprintf(stderr, "AS_UID_addUID()-- failed to insert uid '%s' into store; already there?!\n", uid.UIDstring);
      assert(0);
    }
  }

  uid.isString  = 1;
  uid.UID       = loc;
  uid.UIDstring = NULL;  //  its not valid until we get().

  return(uid);
}





static
inline
char *
AS_UID_toStringBuffer(AS_UID uid, char *buffer) {

  if (uid.UIDstring)
    return(uid.UIDstring);

  assert(buffer != NULL);

  if (uid.isString) {
    if ((uid.UIDstring == NULL) && (AS_UID_gkp))
      uid = AS_GKP_getUID(AS_UID_gkp, uid);

    if ((uid.UIDstring == NULL) && (AS_UID_STRtoUID))
      uid = AS_UID_getUID(uid);

    if (uid.UIDstring) {
      sprintf(buffer, "%s", uid.UIDstring);
    } else {
      sprintf(buffer, "CAx%llu", uid.UID);
    }
  } else {
    sprintf(buffer, "%llu", uid.UID);
  }

  return(buffer);
}

//  A very common operation is to print a bunch of UIDs at the same
//  time.  We allow printing of one UID (AS_UID_toString()), or up to
//  three.
//
//  This is HUGELY error prone.  What was I thinking?

char *
AS_UID_toString(AS_UID uid) {
  static  char  localbuffer[MAX_UID_LENGTH + 1];
  return(AS_UID_toStringBuffer(uid, localbuffer));
}

char *
AS_UID_toString1(AS_UID uid) {
  static  char  localbuffer[MAX_UID_LENGTH + 1];
  return(AS_UID_toStringBuffer(uid, localbuffer));
}

char *
AS_UID_toString2(AS_UID uid) {
  static  char  localbuffer[MAX_UID_LENGTH + 1];
  return(AS_UID_toStringBuffer(uid, localbuffer));
}

char *
AS_UID_toString3(AS_UID uid) {
  static  char  localbuffer[MAX_UID_LENGTH + 1];
  return(AS_UID_toStringBuffer(uid, localbuffer));
}




//  Converts the next word in a string to a uid, returns a pointer to
//  the string after the uid.  If the uid is known to the assembler, a
//  valid uid is returned.  Otherwise, the uid holds a pointer to the
//  string containing the uid, and the uid should be passed to
//  AS_UID_load() to save it in the assembler.
//
AS_UID
AS_UID_lookup(char *uidstr, char **nxtstr) {
  AS_UID  uid = AS_UID_undefined();

  assert(uidstr != NULL);
  //  Bump up to the next word, then bump past the word.
  //
  while (*uidstr && isspace(*uidstr))  uidstr++;

  //  If we are all numeric, and smaller than what will fit in
  //  63-bits, then we're numeric, otherwise, we're a string.

  int   len = 0;
  int   dig = 0;

  while (uidstr[len] && !isspace(uidstr[len])) {
    if (('0' <= uidstr[len]) && (uidstr[len] <= '9'))
      dig++;
    len++;
  }

  if (len > MAX_UID_LENGTH) {
    fprintf(stderr, "ERROR:  UID '%s' larger than maximum allowed (%d letters).\n",
            uidstr, MAX_UID_LENGTH);
    exit(1);
  }

  int  isInteger  = 0;
  int  isInternal = 0;

  //  If all digits, and waaay below 2^64 - 1, we're integer.
  if ((len == dig) &&
      (len < 19))
    isInteger = 1;

  //  If all digits, and barely below 2^64 - 1, we're integer.
  if ((len == dig) &&
      (len == 19) &&
      (strcmp(uidstr, "9223372036854775807") <= 0))
    isInteger = 1;

  //  Integer?  Convert.  Otherwise, look it up in gatekeeper.
  if (isInteger) {
    uid.isString  = 0;
    uid.UID       = strtoull(uidstr, NULL, 10);
    uid.UIDstring = NULL;
  } else if (AS_UID_gkp) {
    uid.isString  = 0;
    uid.UID       = 0;
    uid.UIDstring = uidstr;
    uid = AS_GKP_getUIDfromString(AS_UID_gkp, uid);
  } else {
    uid.isString  = 0;
    uid.UID       = 0;
    uid.UIDstring = uidstr;
    uid = AS_UID_getUIDfromString(uid);
  }

  //  Now bump past the uid, almost like strtoull does.
  if (nxtstr) {
    while (*uidstr && !isspace(*uidstr))  uidstr++;
    while (*uidstr &&  isspace(*uidstr))  uidstr++;
    *nxtstr = uidstr;
  }

  return(uid);
};



//  Adds a uid to the gatekeeper store.  If the uid is already known
//  (or doesn't need to be saved -- if it's just an integer), the
//  function does nothing.  Do NOT use this for general UID lookups!
//  Use AS_UID_lookup() instead.
//
AS_UID
AS_UID_load(char *uidstr) {
  AS_UID  uid = AS_UID_lookup(uidstr, NULL);

  if ((uid.UID > 0) || (uid.isString == 1))
    return(uid);

  //  Brand new uid.  Add it.

  uid.isString  = 0;
  uid.UID       = 0;
  uid.UIDstring = uidstr;

  if (AS_UID_gkp)
    uid = AS_GKP_addUID(AS_UID_gkp, uid);
  else
    uid = AS_UID_addUID(uid);

  return(uid);
}
