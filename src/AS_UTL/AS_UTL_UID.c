
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

// $Id: AS_UTL_UID.c,v 1.1 2007-11-08 12:38:16 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#define MAX_UID_LENGTH (128)


static
GateKeeperStore *AS_UID_gkp = NULL;

void
AS_UID_setGatekeeper(void *gkp) {
  AS_UID_gkp = gkp;
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

  //  Bump up to the next word, then bump past the word.
  //
  while (uidstr && isspace(*uidstr))  uidstr++;

  //  If we are all numeric, and smaller than what will fit in
  //  63-bits, then we're numeric, otherwise, we're a string.

  int   len = 0;
  int   dig = 0;

  while (uidstr[len]) {
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
  } else {
    if (AS_UID_gkp == NULL) {
      fprintf(stderr, "AS_UID_lookup()-- cannot lookup string UID from GateKeeper; AS_UID_gkp is not set.\n");
      exit(1);
    }
    uid.isString  = 0;
    uid.UID       = 0;
    uid.UIDstring = uidstr;
    uid = AS_GKP_getUIDfromString(AS_UID_gkp, uid);
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

  if (AS_UID_gkp == NULL) {
    fprintf(stderr, "AS_UID_load()-- cannot load string UID into GateKeeper; AS_UID_gkp is not set.\n");
    exit(1);
  }

  uid = AS_GKP_addUID(AS_UID_gkp, uid);

  return(uid);
}
