
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

static char *rcsid = "$Id$";

#include "AS_PER_gkpStore.H"


//  A hash must be maintained that maps the read IIDS to the location in
//  the PLC store.  This will be used to lookup a PLC record.
//
//  The gkp files are:
//    'f2p' -> the IID to PLC mapping


////////////////////////////////////////////////////////////////////////////////
//
//  FRG to PLC lookups.
//
void
gkStore::gkStore_loadFRGtoPLC(void) {
  if (FRGtoPLC == NULL) {
    char  name[FILENAME_MAX];
    sprintf(name,"%s/f2p", storePath);
    // even though this isn't a UID to IID hash table, we use the same (int-based) indexing functions
    FRGtoPLC = LoadUIDtoIIDHashTable_AS(name);
  }
  assert(FRGtoPLC != NULL);
}

int 
gkStore::gkStore_getFRGtoPLC(AS_IID iid) {
  uint64   plid = 0;
  
  gkStore_loadFRGtoPLC();
  if (AS_IID_isDefined(iid))
    LookupInHashTable_AS(FRGtoPLC, iid, 0, &plid, 0);
    
  return ((int)plid);
}

int
gkStore::gkStore_setFRGtoPLC(AS_IID iid, int plid) {
  gkStore_loadFRGtoPLC();
  assert(AS_IID_isDefined(iid) == TRUE);
  assert(plid != 0);
  return(InsertInHashTable_AS(FRGtoPLC, (uint64)iid, 0, (uint64)plid, 0));
}

void
gkStore::gkStore_addPlacement(gkPlacement *pl) {
  appendIndexStore(plc, pl);
  gkStore_setFRGtoPLC(pl->frag, getLastElemStore(plc));
  inf.plcLoaded++;
}
