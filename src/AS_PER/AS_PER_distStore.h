
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
#ifndef AS_PER_DISTSTORE_H
#define AS_PER_DISTSTORE_H
/*************************************************************************
 Module:  AS_PER_distStore
 Description:
    A thin layer on top of the IndexStore supporing the storage and
 retrieval of DST records.
    The idea is to provide easier to use shortcuts for the common
 operations, and let the other operations be accessed through the
 generic Index Store API.

 Assumptions:
    Nothing special beyond genericStore.rtf

 Document:
      GenericStore.rtf

 *************************************************************************/

/* RCS Info
 * $Id: AS_PER_distStore.h,v 1.1.1.1 2004-04-14 13:52:50 catmandew Exp $
 * $Revision: 1.1.1.1 $
 *
 */


#include <time.h>


#include "AS_global.h"
#include "AS_PER_genericStore.h"

typedef struct{
  unsigned int   deleted :1;
  unsigned int   spare   :31;
  uint32 IID;            /* Internal ID of this read */
  uint64 UID;             /* Accession ID of this read */
  float32 mean;
  float32 stddev;
}DistRecord;


typedef StoreHandle DistStore;

/***********************************************************************************
 * Function: createDistStore:
 * Description:
 *     Allocates an index store for the DST Store, and returns its handle.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     firstID     Index of the first element in the store.
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static DistStore createDistStore(char *StorePath, int firstID){
  StoreHandle distStore = createIndexStore(StorePath, "dst", sizeof(DistRecord), 1, firstID);
  return(distStore);
}


/***********************************************************************************
 * Function: openDistStore:
 * Description:
 *     Opens an existing, file-based store.
 *
 * Inputs:
 *     StorePath   path to the indexStore.  If NULL this is a memory Store.
 *     rw          file access mode
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static DistStore openDistStore(char *StorePath, /* Path to file */
		      char *rw             /* "r" or "rw" */
			  ){
  return openStore(StorePath, rw);

}


/***********************************************************************************
 * Function: commitDistStore:
 * Description:
 *     Like closeStore, but store remains open.
 *
 * Inputs:
 *     sh    Handle to store we want to commit
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static int commitDistStore(DistStore sh){
  return commitStore(sh);
}

/***********************************************************************************
 * Function: resetDistStore:
 * Description:
 *     Recycles a Dist store, nuking its data.
 *
 * Inputs:
 *     sh           handle of open Dist Store
 *     firstID      First ID for reset Store
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static DistStore resetDistStore(DistStore sh, int firstID){
  return resetIndexStore(sh, firstID);
}


/***********************************************************************************
 * Function: closeDistStore:
 * Description:
 *     Close an open store.  Updates lastElem and lastCommittedElem.
 *
 * Inputs:
 *     sh    Handle to store we want to close
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static int closeDistStore(DistStore sh){
  return closeStore(sh);
}


/***********************************************************************************
 * Function: getDistStore
 * Description:
 *     Random access to records in a Dist store
 *
 * Inputs:
 *     fs         Handle of Dist store
 *     index      index of record
 * Outputs:
 *     dr     Buffer for element 
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/
	
static int getDistStore(DistStore fs, int index, DistRecord *dr){
  return getIndexStore(fs,index,dr);
}

/***********************************************************************************
 * Function: setDistStore
 * Description:
 *     Overwrite an existing  element of a Dist store
 *
 * Inputs:
 *     store      Handle of Dist store
 *     index      index of element to overwrite
 *     element    Pointer to record for overwrite
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static int setDistStore(DistStore store, int index, DistRecord *element){

  return setIndexStore(store,index, element);
}


/***********************************************************************************
 * Function: deleteDistStore
 * Description:
 *     Delete an element from an index store
 *
 * Inputs:
 *     store      Handle of index store
 *     index      Index of element to be deleted
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static int deleteDistStore(DistStore store, int index){
  DistRecord dr;

  // Read, modify and write the
  getIndexStore(store, index, &dr);
  dr.deleted = TRUE;
  setIndexStore(store, index, (void *)(&dr));

  return 0;
}


/***********************************************************************************
 * Function: appendDistStore
 * Description:
 *     Append an element to an Dist store
 *
 * Inputs:
 *     store      Handle of Dist store
 *     element    Pointer to record for append
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/

static int appendDistStore(DistStore store, DistRecord *element){
  return appendIndexStore(store,element);
}






	

#endif
