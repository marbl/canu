
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

static char *rcsid = "$Id: AS_PER_gkStore_UID.C,v 1.5 2011-08-01 22:33:01 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"


//  For numeric UIDs (those with both isString and UID valid), you can
//  do a UID to IID mapping.
//
//  String UIDs must be added to the store first.  The string itself
//  is stored in a StringStore just like sequence and source.  The
//  index into this store is saved in an AS_UID.  The string you pass
//  in can be a temporary.
//
//  A hash must be maintained that maps the strings to the location in
//  the store (which is also exactly the (isString, UID) pair.  This
//  will be used to lookup a string UID.
//
//  The gkp files are:
//    'u2i' -> the UID to IID mapping
//    'uid' -> the string store itself


////////////////////////////////////////////////////////////////////////////////
//
//  UID to IID lookups.
//
void
gkStore::gkStore_loadUIDtoIID(void) {
  if (UIDtoIID == NULL) {
     if (doNotLoadUIDs == FALSE) {
          char  name[FILENAME_MAX];
          sprintf(name,"%s/u2i", storePath);
          UIDtoIID = LoadUIDtoIIDHashTable_AS(name);
     } else {
       UIDtoIID = CreateScalarHashTable_AS();
     }
  }
  assert(UIDtoIID != NULL);
}

AS_IID
gkStore::gkStore_getUIDtoIID(AS_UID uid, uint32 *type) {
  uint64   iid = 0;
  gkStore_loadUIDtoIID();
  if (AS_UID_isDefined(uid))
    LookupInHashTable_AS(UIDtoIID, AS_UID_toInteger(uid), 0, &iid, type);
  return((AS_IID)iid);
}

int
gkStore::gkStore_setUIDtoIID(AS_UID uid, AS_IID iid, uint32 type) {
  gkStore_loadUIDtoIID();
  assert(AS_UID_isDefined(uid) == TRUE);
  assert(AS_IID_isDefined(iid) == TRUE);
  return(InsertInHashTable_AS(UIDtoIID, AS_UID_toInteger(uid), 0, (uint64)iid, type));
}



////////////////////////////////////////////////////////////////////////////////
//
//  IID to UID lookups
//
void
gkStore::gkStore_loadIIDtoUID(void) {

  if (frgUID)
    return;

  frgUID = (uint64 *)safe_calloc(inf.frgLoaded + 1, sizeof(uint64));

  HashTable_Iterator_AS   iterator  = {0};
  uint64                  key       = 0;
  uint64                  value     = 0;
  uint32                  valuetype = 0;
  uint32                  added     = 0;

  gkStore_loadUIDtoIID();
  InitializeHashTable_Iterator_AS(UIDtoIID, &iterator);

  while (NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)) {
    if (valuetype == AS_IID_FRG)
      frgUID[value] = key;
  }
}

AS_UID
gkStore::gkStore_getIIDtoUID(AS_IID iid, uint32 type) {
  AS_UID  uid = AS_UID_undefined();

  gkStore_loadIIDtoUID();

  switch (type) {
    case AS_IID_FRG:
      uid = AS_UID_fromInteger(frgUID[iid]);
      break;
    case AS_IID_LIB:
      uid = gkStore_getLibrary(iid)->libraryUID;
      break;
    default:
      break;
  }

  return(uid);
}



////////////////////////////////////////////////////////////////////////////////
//
//  UIDtoIID rebuild
//
void
gkStore::gkStore_rebuildUIDtoIID(void) {
  AS_IID      i, f, l;

  if (UIDtoIID) {
    ResetHashTable_AS(UIDtoIID);
  } else {
    char name[FILENAME_MAX];
    sprintf(name,"%s/u2i", storePath);
    UIDtoIID = CreateScalarHashTable_AS();
    SaveHashTable_AS(name, UIDtoIID);
  }

  //
  //  Insert library info
  //

  i = 0;
  f = getFirstElemStore(lib);
  l = getLastElemStore(lib);

  for (i=f; i<=l; i++) {
    gkLibrary *L = gkStore_getLibrary(i);

    if (InsertInHashTable_AS(UIDtoIID, AS_UID_toInteger(L->libraryUID), 0, (uint64)i, AS_IID_LIB) == HASH_FAILURE)
      fprintf(stderr, "Error inserting library %s,"F_IID" into hash table.\n",
              AS_UID_toString(L->libraryUID), i);
  }

  //
  //  Insert fragment info
  //

  gkFragment   fr;
  gkStream     gs(this, 0, 0, GKFRAGMENT_INF);

  while (gs.next(&fr)) {
    if (InsertInHashTable_AS(UIDtoIID, AS_UID_toInteger(fr.gkFragment_getReadUID()), 0, (uint64)fr.gkFragment_getReadIID(), AS_IID_FRG) == HASH_FAILURE)
      fprintf(stderr, "Error inserting UID %s and IID "F_IID"into hash table.\n",
              AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getReadIID());
  }
}



////////////////////////////////////////////////////////////////////////////////
//
//  AS_MSG utility functions
//
void
gkStore::gkStore_loadSTRtoUID(void) {  
  if (doNotLoadUIDs == FALSE) {
     if (STRtoUID == NULL) {
       uid = convertStoreToMemoryStore(uid);
   
       STRtoUID = CreateStringHashTable_AS();
   
       char          *uidptr = NULL;
       int64          uidoff = 1;
       int64          nxtoff = 1;
       uint32         actlen = 0;
   
       while ((uidptr = getStringStorePtr(uid, uidoff, &actlen, &nxtoff)) != NULL) {
         if (strlen(uidptr) != actlen) {
           int i;
           fprintf(stderr, "gkStore_loadSTRtoUID()-- string '%s' length %d != stored actlen = %d\n",
                   uidptr, strlen(uidptr), actlen);
           for (i=0; i<strlen(uidptr); i++)
             fprintf(stderr, "[%2d] %3d '%c'\n", i, uidptr[i], uidptr[i]);
         }
         assert(strlen(uidptr) == actlen);
   
         if (InsertInHashTable_AS(STRtoUID,
                                  (INTPTR)uidptr, actlen,
                                  uidoff, 0) == HASH_FAILURE) {
           fprintf(stderr, "gkStore_loadSTRtoUID()-- failed to insert uid '%s' into store; already there?!\n", uidptr);
           assert(0);
         }
         uidoff = nxtoff;
       }
     }
     assert(STRtoUID != NULL);
  }
}


//  Given a string, returns the AS_UID for it.
//
AS_UID
gkStore::gkStore_getUIDfromString(char *uidstr) {
  AS_UID  uid = AS_UID_undefined();
  
  if (doNotLoadUIDs == TRUE) {
    return uid;
  }
  
  uint64  loc = 0;
  uint64  len = 0;
  char    end = 0;

  gkStore_loadSTRtoUID();

  //  A common error (especially when reading from files) it to leave whitespace at the ends.  We
  //  temporarily trim it off.

  while (*uidstr && isspace(*uidstr))
    uidstr++;

  while (uidstr[len] && !isspace(uidstr[len]))
    len++;

  end         = uidstr[len];
  uidstr[len] = 0;

  if (LookupInHashTable_AS(STRtoUID, (INTPTR)uidstr, len, &loc, 0)) {
    uid.isString  = 1;
    uid.UID       = loc;
  }

  uidstr[len] = end;

  return(uid);
}


//  Given an AS_UID, returns a pointer to the string.
//
char *
gkStore::gkStore_getUIDstring(AS_UID u) {
  char *retval = NULL;

  if (!doNotLoadUIDs && u.isString) {
    uint32  actlen = 0;
    int64   uidoff = 0;

    gkStore_loadSTRtoUID();
    retval = getStringStorePtr(uid, u.UID, &actlen, &uidoff);
  }

  return(retval);
}


//  Given a string, creates a new AS_UID for it.  If the UID already
//  exists, a duplicate is NOT added.
//
AS_UID
gkStore::gkStore_addUID(char *uidstr) {
  uint64     loc = 0;
  uint64     len = 0;
  char       end = 0;

  if (doNotLoadUIDs == TRUE) {
    fprintf(stderr, "gkStore_addUID: UID string store is turned off but it is being added to.\n");
    assert(0); 
  }
  
  //  A common error (especially when reading from files) it to leave whitespace at the ends.  We
  //  temporarily trim it off.

  while (*uidstr && isspace(*uidstr))
    uidstr++;

  while (uidstr[len] && !isspace(uidstr[len]))
    len++;

  end         = uidstr[len];
  uidstr[len] = 0;

  gkStore_loadSTRtoUID();

  //  If the UID is already in the store, just return as if it was
  //  new.  Otherwise, add it to the store.

  if (LookupInHashTable_AS(STRtoUID, (INTPTR)uidstr, len, &loc, 0) == FALSE) {
    char    *str = NULL;
    uint32   act = 0;
    int64    off = 0;

    loc = getLastElemStore(uid) + 1;

    //  Stash the UID on disk.
    off = appendStringStore(uid, uidstr, len);

    //  If our string store changed, update all the pointers in our hash table.
    if (off)
      UpdatePointersInHashTable_AS(STRtoUID, off);

    str = getStringStorePtr(uid, loc, &act, &off);

    if (InsertInHashTable_AS(STRtoUID,
                             (INTPTR)str, len,
                             loc, 0) == HASH_FAILURE) {
      fprintf(stderr, "gkStore_addUID()-- failed to insert uid '%s' into store; already there?!\n", uidstr);
      assert(0);
    }
  }

  AS_UID     u;

  u.isString  = 1;
  u.UID       = loc;

  uidstr[len] = end;

  return(u);
}
