
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
static char CM_ID[] = "$Id: AS_GKP_checkBac.c,v 1.1.1.1 2004-04-14 13:51:36 catmandew Exp $";

//#define DEBUG_GKP 1
#define DEBUG_GKP_VERBOSE 1
/*************************************************
* Module:  AS_GKP_check.c
* Description:
*    Gatekeeper check routines.
* 
*    Programmer:  S. Kravitz
*       Written:  Jan 1999
* 
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <ctype.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"
#include "AS_UTL_Hash.h"

/****************************************************************************/
/* These are the workhorse checker routines for the BAC messages.
   We limit the size and complexity of this routine by dealing with Adds, Deletes,
   and redefines in separate routines.
*/
int Check_AddBacMesg(BacMesg *bac_mesg,
		     InternalBacMesg *ibc_mesg,
		     VA_TYPE(int32) *addedBactigs,
		     CDS_CID_t batchID,
		     time_t currentTime,
		     int32 assembler,
		     int32 strict,
		     int verbose);

int Check_DeleteBacMesg(BacMesg *bac_mesg,
			InternalBacMesg *ibc_mesg,
			CDS_CID_t batchID,
			time_t currentTime,
			int32 assembler,
			int32 strict,
			int verbose);

int Check_RedefineBacMesg(BacMesg *bac_mesg,
			  InternalBacMesg *ibc_mesg,
			  VA_TYPE(int32) *addedBactigs,
			  CDS_CID_t batchID,
			  time_t currentTime,
			  int32 assembler,
			  int32 strict,
			  int verbose);

/* Utility routines for output of delete messages required for the overlay assembler */
void outputBactigFragDeletes(GateKeeperLocaleRecord bacRecord);
void outputBacDelete(CDS_IID_t bacIid, CDS_UID_t bacUid);


/****************************************************************************/

int Check_BacMesg(BacMesg *bac_mesg,
		  InternalBacMesg *ibc_mesg,
		  VA_TYPE(int32) *addedBactigs,
		  CDS_CID_t batchID,
		  time_t currentTime,
		  int assembler,
		  int strict,
		  int verbose){



  *ibc_mesg = *bac_mesg;
  ibc_mesg->ibac_id = 0;
  ibc_mesg->iseq_id = 0;
  ibc_mesg->ilength = 0;

  switch(bac_mesg->action){
  case AS_DELETE:
    /*************************** DELETE *******************************/
    return Check_DeleteBacMesg(bac_mesg, ibc_mesg,  batchID,currentTime, assembler, strict, verbose);
    

    break;
  case AS_ADD:
    /********************** ADD ***********************/
    return Check_AddBacMesg(bac_mesg, ibc_mesg,  addedBactigs, batchID,currentTime, assembler, strict, verbose);
    break;
  case AS_REDEFINE:
    /*************************** oREDEFINE *******************************/
    return Check_RedefineBacMesg(bac_mesg, ibc_mesg,  addedBactigs, batchID,currentTime, assembler, strict, verbose);
    break;
  default:
    fprintf(Msgfp,"# Check_BacMesg: Invalid action %c\n",
	    bac_mesg->action);
  }

  return GATEKEEPER_FAILURE;
}

/*****************************************************************************
 *  Special Handling for FBACs for OVERLAY
 *
 *  If an FBAC is input, we convert it into a UBAC as follows:
 *     - allocate a Bactig with UID 0 in the bactig store
 *     - set numbactigs to 1 and add a bactig to the list
 *     - change the type to UNFINISHED
 *  When we see the subsequent FULL_BAC fragment, we convert it
 *  to a BACTIG fragment, and number it with the bactig ID.
 ****************************************************************************/

/****************************************************************************/
int Check_AddBacMesg(BacMesg *bac_mesg,
		     InternalBacMesg *ibc_mesg,
		     VA_TYPE(int32) *addedBactigs,
		     CDS_CID_t batchID,
		     time_t currentTime,
		     int assembler,
		     int strict,
		     int verbose)
{
  PHashValue_AS value;

  switch(bac_mesg->type){
  case AS_ENDS:
  case AS_LIGHT_SHOTGUN:
  case AS_UNFINISHED:
  case AS_FINISHED:
    break;

  default:
    fprintf(Msgfp,"# Check_BacMessage:  Trying to add unknown BAC type for Bac " F_UID " type %c \n",
	    ibc_mesg->ebac_id,
	    ibc_mesg->type);
    return(GATEKEEPER_FAILURE);
  }

  // Make sure bac_id and seq_id are unique
  /* Make sure we haven't seen this frag record before... if so
     it is a fatal error */
  if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bac_mesg->ebac_id, 
					       AS_IID_LOC, 
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# Check_BacMessage:  A message with UID " F_UID " exists %s!!! Can't add it... bye\n",
	    bac_mesg->ebac_id, (value.deleted?"and was deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  ibc_mesg->ibac_id = value.IID;

  if(bac_mesg->type != AS_ENDS &&
     bac_mesg->type != AS_LIGHT_SHOTGUN){
    if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						 UID_NAMESPACE_AS,
						 bac_mesg->eseq_id, 
						 AS_IID_SEQ, 
						 FALSE,
						 Msgfp,
						 &value)){

      fprintf(Msgfp,"# Check_BacMessage:  A message with UID " F_UID " exists %s!!! Can't define it... bye\n",
	      bac_mesg->ebac_id, (value.deleted?"and was deleted":""));
      return(GATEKEEPER_FAILURE);
    }
    ibc_mesg->iseq_id = value.IID;
  }else{
    if(verbose)
      fprintf(stderr,"* seqID = 0\n");
    ibc_mesg->iseq_id = 0;
  }      
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bac_mesg->elength, 
					       AS_IID_DST, 
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# Check_BacMessage:  Distance " F_UID " DOES NOT exist %s!!! Can't reference it... bye\n",
	    bac_mesg->elength, (value.deleted?"and was deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  ibc_mesg->ilength = value.IID;

  switch(bac_mesg->type){
  
  case AS_LIGHT_SHOTGUN:
    if(strict &&
       assembler == AS_ASSEMBLER_OVERLAY){
      printGKPError(Msgfp,GKPError_BACWrongTypeForOverlay);
      return(GATEKEEPER_FAILURE);
    }
    /*** INTENTIONAL FALL THROUGH ***/
  case AS_ENDS:  // we now handle these for overlay
  case AS_FINISHED:   //We now handle these as a special case
    if(bac_mesg->num_bactigs != 0){
      printGKPError(Msgfp,GKPError_BACNumBactigsShouldBeZero);
      fprintf(Msgfp,"# Check_BacMessage:  num_bactigs is %d must be 0\n", bac_mesg->num_bactigs);
      return(GATEKEEPER_FAILURE);
    }
    break;
  case AS_UNFINISHED:
    {
      int i;
      HashTable_AS *btgHash = CreateHashTable_uint64_AS(32);
      if(bac_mesg->num_bactigs == 0){
	printGKPError(Msgfp,GKPError_BACNumBactigs);
	fprintf(Msgfp,"# Check_BacMessage:  num_bactigs is %d must > 0\n", bac_mesg->num_bactigs);
	DeleteHashTable_AS(btgHash);
	return(GATEKEEPER_FAILURE);
      }
      for(i = 0; i < bac_mesg->num_bactigs; i++){
	InternalBactigMesg *bactig = bac_mesg->bactig_list + i;
	
	if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						     UID_NAMESPACE_AS,
						     bactig->eaccession, 
						     AS_IID_BTG, 
						     FALSE,
						     Msgfp,
						     &value)){

	  printGKPError(Msgfp,GKPError_BadUniqueBAC);
	  fprintf(Msgfp,"# Check_BacMessage:  A message with UID " F_UID " ALREADY exists %s!!! Can't define it... bye\n",
		  bactig->eaccession, (value.deleted?"and was deleted":""));
	  DeleteHashTable_AS(btgHash);
	  return(GATEKEEPER_FAILURE);
	}
	if(HASH_SUCCESS != InsertInHashTable_AS(btgHash, &bactig->eaccession, sizeof(bactig->eaccession), &bactig)){
	  printGKPError(Msgfp,GKPError_BACMultiplyDefinedBactig);
	  fprintf(Msgfp,"# Check_BacMessage:  Bactig " F_UID " defined more than once in this BAC!!! Can't define it... bye\n",
		  bactig->eaccession);
	  DeleteHashTable_AS(btgHash);
	  
	  return(GATEKEEPER_FAILURE);
	}
      }
      DeleteHashTable_AS(btgHash);
    }
    break;
  default:
    assert(0);
  }

  // OK, we passed all of the tests.  Add symbols and references
  {
    PHashValue_AS locValue;
    GateKeeperLocaleRecord bacRecord;
    GateKeeperSequenceRecord seqRecord;

    locValue.type = AS_IID_LOC;
    if(HASH_SUCCESS != InsertInPHashTable_AS(&(GkpStore.hashTable),
                                             UID_NAMESPACE_AS,
                                             bac_mesg->ebac_id,
                                             &locValue, FALSE, TRUE))
      assert(0);

    ibc_mesg->ibac_id = locValue.IID;
    bacRecord.hasSequence = FALSE;
    bacRecord.deleted = FALSE;
    bacRecord.redefined = FALSE;
    bacRecord.birthBatch = batchID;
    bacRecord.deathBatch = 0;
    bacRecord.prevID = 0;
    bacRecord.prevInstanceID = 0;

    bacRecord.UID = bac_mesg->ebac_id;
    bacRecord.isBac = TRUE;
    bacRecord.type = bac_mesg->type;
    bacRecord.numBactigs = 0;
    bacRecord.firstBactig = 0;
    bacRecord.sequenceID = 0;
    if(bac_mesg->type == AS_FINISHED ||
       bac_mesg->type == AS_UNFINISHED){
      locValue.type = AS_IID_SEQ;
      if(HASH_SUCCESS != InsertInPHashTable_AS(&(GkpStore.hashTable),
                                               UID_NAMESPACE_AS,
                                               bac_mesg->eseq_id,
                                               &locValue, FALSE, TRUE))
        assert(0);
      ibc_mesg->iseq_id = locValue.IID;
      if(verbose)
	fprintf(stderr,"* Inserted seq " F_UID " in hashtable with IID " F_IID "\n",
		bac_mesg->eseq_id, ibc_mesg->iseq_id);
      seqRecord.deleted = FALSE;
      seqRecord.UID = bac_mesg->eseq_id;
      seqRecord.localeID = ibc_mesg->ibac_id;
      bacRecord.sequenceID = locValue.IID;
      if(verbose)
	fprintf(stderr,"* Sequence " F_IID " associated with bac " F_IID "\n",
		bacRecord.sequenceID, seqRecord.localeID);
    }
    bacRecord.lengthID = ibc_mesg->ilength;
    AddRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, bac_mesg->elength);
    {
      int i;
      for(i = 0; i < bac_mesg->num_bactigs; i++){
	InternalBactigMesg *bactig = bac_mesg->bactig_list + i;
	GateKeeperBactigRecord btgRecord;

	PHashValue_AS btgValue;

	btgValue.type = AS_IID_BTG;
	if(HASH_SUCCESS != InsertInPHashTable_AS(&(GkpStore.hashTable),
                                                 UID_NAMESPACE_AS,
                                                 bactig->eaccession,
                                                 &btgValue, FALSE, TRUE))
          assert(0);
	bactig->iaccession = btgValue.IID;
	btgRecord.hasSequence = FALSE;
	btgRecord.deleted = FALSE;
	btgRecord.bacID = ibc_mesg->ibac_id;
	btgRecord.seqID = ibc_mesg->iseq_id;
	btgRecord.UID = bactig->eaccession;
	btgRecord.length = bactig->length;
	if(verbose)
	  fprintf(stderr,"* Inserting bactig " F_UID " in hashtable bactigID = " F_IID "  localeID=" F_IID "\n", 
		  bactig->eaccession, bactig->iaccession, btgRecord.bacID);
	appendGateKeeperBactigStore(GkpStore.btgStore,&btgRecord);
	if(assembler == AS_ASSEMBLER_OVERLAY &&
	   bac_mesg->type == AS_UNFINISHED){
	  
	  Appendint32(addedBactigs, (int32 *)&(bactig->iaccession)); // Add to list of bactigs that were added
	}
      }

      bacRecord.firstBactig = (bac_mesg->bactig_list?bac_mesg->bactig_list[0].iaccession:0);
      bacRecord.numBactigs = bac_mesg->num_bactigs;
    }
    /*******vvvvvvv****** SPECIAL HANDLING *****************vvvvvvv**********/

    if(/* assembler = AS_ASSEMBLER_OVERLAY && */ bac_mesg->type == AS_FINISHED){
      static InternalBactigMesg dummy;
      GateKeeperBactigRecord btgRecord;

      // Create a dummy bactig 
      dummy.eaccession = 0;
      dummy.iaccession = AllocateCountPHashTable_AS(GkpStore.hashTable,AS_IID_BTG);
      dummy.length = 0;
      ibc_mesg->bactig_list = &dummy;
      ibc_mesg->num_bactigs = 1;
      bac_mesg->type = ibc_mesg->type = AS_UNFINISHED;
      bacRecord.numBactigs = 1;
      bacRecord.firstBactig = dummy.iaccession;
      btgRecord.hasSequence = FALSE;
      btgRecord.deleted = FALSE;
      btgRecord.bacID = ibc_mesg->ibac_id;
      btgRecord.seqID = ibc_mesg->iseq_id;
      btgRecord.UID = 0;
      btgRecord.length = 0;
      if(verbose)
	fprintf(stderr,"* Inserting *** DUMMY **** bactig " F_UID " in hashtable bactigID = " F_IID "  localeID=" F_IID "\n", 
		dummy.eaccession, dummy.iaccession, btgRecord.bacID);
      appendGateKeeperBactigStore(GkpStore.btgStore,&btgRecord);

      Appendint32(addedBactigs, (int32 *)&(dummy.iaccession)); // Add to list of bactigs that were added
    }
    /*******^^^^^^^^^******** SPECIAL HANDLING *************^^^^^^^^**********/

    appendGateKeeperLocaleStore(GkpStore.locStore, &bacRecord);
    if(bac_mesg->type == AS_FINISHED ||
       bac_mesg->type == AS_UNFINISHED){
      appendGateKeeperSequenceStore(GkpStore.seqStore, &seqRecord);
    }

  }

#ifdef DEBUG_LOC
  fprintf(stderr,"* Inserted loc UID " F_UID " got iid " F_IID "\n",
	  frg_mesg->locale, ifg_mesg->ilocale);
#endif
  return GATEKEEPER_SUCCESS;
   
}

/****************************************************************************/
int Check_DeleteBacMesg(BacMesg *bac_mesg,
			InternalBacMesg *ibc_mesg,
			CDS_CID_t batchID,
			time_t currentTime,
			int32 assembler,
			int32 strict,
			int verbose)
{
  PHashValue_AS value;
  GateKeeperLocaleRecord gkpbac;
  // Check that UID exists and is right type
  /* Make sure we haven't seen this frag record before... if so
     it is a fatal error */
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bac_mesg->ebac_id, 
					       AS_IID_LOC, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp,GKPError_MissingBAC);
    fprintf(Msgfp,"# Check_BacMessage:  Bac " F_UID " DOES NOT exist!!!" 
	    "# Can't delete it...\n",	    bac_mesg->ebac_id);
    return(GATEKEEPER_FAILURE);

  }else{
    // Check that reference count is zero
    if(value.refCount > 1){
      /*************** REPORT THE REFERENCES **********************/
      
      printGKPError(Msgfp,GKPError_DeleteBAC);
      fprintf(Msgfp,"# Check_BacMessage: There are %d references outstanding to Bac " F_UID ".\n"
	      "#  Can't delete it...\n",
	      value.refCount, bac_mesg->ebac_id);
      return(GATEKEEPER_FAILURE);
    }

    // Mark it as deleted
    getGateKeeperLocaleStore(GkpStore.locStore, value.IID, &gkpbac);

    if(assembler == AS_ASSEMBLER_OVERLAY){
      outputBactigFragDeletes(gkpbac);
    }

    UnRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, gkpbac.UID);
    deleteAndMarkGateKeeperLocaleStore(GkpStore.locStore, value.IID, batchID);
    DeleteFromPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, bac_mesg->ebac_id);


    // Do we need to invalidate all of the consitutent bactigs?

  }
  return GATEKEEPER_SUCCESS;

}


int Check_RedefineBacMesg(BacMesg *bac_mesg,
			  InternalBacMesg *ibc_mesg,
			  VA_TYPE(int32) *addedBactigs,
			  CDS_CID_t batchID,
			  time_t currentTime,
			  int assembler,
			  int strict,
			  int verbose)
{
  // Make sure bac_id the bacID has been previously defined
  // Make sure the seq_id is new
  PHashValue_AS value;

  GateKeeperLocaleRecord localeRecord;


  //  debug_hash(); // ####

  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bac_mesg->ebac_id, 
					       AS_IID_LOC, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp,GKPError_MissingBAC);
    fprintf(Msgfp,"# Check_BacMessage:  Bac " F_UID " does NOT exist %s!!! Can't redefine it..\n",
	    bac_mesg->ebac_id, (value.deleted?"and was deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  ibc_mesg->ibac_id = value.IID;
  getGateKeeperLocaleStore(GkpStore.locStore, value.IID, &localeRecord);
  if(verbose)
    fprintf(stderr,"* Redefining bac " F_UID ", " F_IID " to from  %c->%c\n",
	    bac_mesg->ebac_id, ibc_mesg->ibac_id, localeRecord.type, ibc_mesg->type);


  switch(localeRecord.type){
  case AS_ENDS:
    // All upgrade types are allowed
    break;
  case AS_LIGHT_SHOTGUN:
    if(ibc_mesg->type != AS_FINISHED &&
       ibc_mesg->type != AS_UNFINISHED){
      printGKPError(Msgfp,GKPError_BACRedefinition);
      fprintf(Msgfp,"Check_BacMessage:  Can't redefine a light shotgun BAC to type %c!\n", ibc_mesg->type);
      return GATEKEEPER_FAILURE;
    }
    break;
  case AS_UNFINISHED:
    if(ibc_mesg->type != AS_FINISHED &&
       ibc_mesg->type != AS_UNFINISHED){
      printGKPError(Msgfp,GKPError_BACRedefinition);
      fprintf(Msgfp,"Check_BacMessage:  Can't redefine an unfinished BAC to type %c!\n", ibc_mesg->type);
      return GATEKEEPER_FAILURE;
    }
    break;
  case AS_FINISHED:
    if(ibc_mesg->type != AS_FINISHED){
      printGKPError(Msgfp,GKPError_BACRedefinition);
      fprintf(Msgfp,"Check_BacMessage:  Can't redefine a finished BAC to type %c!\n", ibc_mesg->type);
      return GATEKEEPER_FAILURE;
    }

  default:
    printGKPError(Msgfp,GKPError_Scalar);
    fprintf(Msgfp,"#Check_BacMessage: Unknown BAC type for Bac " F_UID "," F_IID " is %c %d (trying to redefine to %c)\n",
	    ibc_mesg->ebac_id,
	    ibc_mesg->ibac_id,
	    localeRecord.type,
	    localeRecord.type,
	    ibc_mesg->type);
    return GATEKEEPER_FAILURE;
  }



  /* Check sequence IDs */
  if(bac_mesg->type != AS_ENDS &&
     bac_mesg->type != AS_LIGHT_SHOTGUN){
    // Seq id must be new
    if(verbose)
      fprintf(stderr,"* bac_mesg->type = %c UID:" F_UID "\n",
	      bac_mesg->type, ibc_mesg->ebac_id);
    if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						 UID_NAMESPACE_AS,
						 bac_mesg->eseq_id, 
						 AS_IID_SEQ, 
						 FALSE,
						 Msgfp,
						 &value)){

      printGKPError(Msgfp,GKPError_BadUniqueSEQ);
      fprintf(Msgfp,"# Check_BacMessage:  A message with UID " F_UID " exists %s!!! Can't define it... bye\n",
	      bac_mesg->eseq_id, (value.deleted?"and was deleted":""));
      return(GATEKEEPER_FAILURE);
    }
    ibc_mesg->iseq_id = value.IID;
  }else{
    ibc_mesg->iseq_id = 0;
  }      


  /* Check Length */
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bac_mesg->elength, 
					       AS_IID_DST, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp,GKPError_MissingDST);
    fprintf(Msgfp,"# Check_BacMessage:  Distance " F_UID " DOES NOT exist %s!!! Can't reference it... bye\n",
	    bac_mesg->elength, (value.deleted?"and was deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  ibc_mesg->ilength = value.IID;


      /* Check out Bactigs */

  switch(bac_mesg->type){
  case AS_ENDS:
  case AS_LIGHT_SHOTGUN:
  case AS_FINISHED:
    if(bac_mesg->num_bactigs != 0){
      printGKPError(Msgfp,GKPError_BACNumBactigsShouldBeZero);
      fprintf(Msgfp,"# Check_BacMessage:  num_bactigs is %d must be 0\n", bac_mesg->num_bactigs);
      return(GATEKEEPER_FAILURE);
    }
    break;
  case AS_UNFINISHED:
    {
      int i;
      if(bac_mesg->num_bactigs == 0){
	printGKPError(Msgfp,GKPError_BACNumBactigs);
	fprintf(Msgfp,"# Check_BacMessage:  num_bactigs is %d must > 0\n", bac_mesg->num_bactigs);
	return(GATEKEEPER_FAILURE);
      }
      for(i = 0; i < bac_mesg->num_bactigs; i++){
	InternalBactigMesg *bactig = bac_mesg->bactig_list + i;

	if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						     UID_NAMESPACE_AS,
						     bactig->eaccession, 
						     AS_IID_BTG, 
						     FALSE,
						     Msgfp,
						     &value)){

	  printGKPError(Msgfp,GKPError_BadUniqueBTG);
	  fprintf(Msgfp,"# Check_BacMessage:  A message with UID " F_UID " ALREADY exists %s!!! Can'tx define it... bye\n",
		  bactig->eaccession, (value.deleted?"and was deleted":""));
	  return(GATEKEEPER_FAILURE);
	}
	 

      }
    }
    break;
  default:
    assert(0);
  }

  // OK, we passed all of the tests.  Add symbols and references
  {
    PHashValue_AS locValue;
    GateKeeperLocaleRecord bacRecord;
    GateKeeperSequenceRecord seqRecord;

    locValue.type = AS_IID_LOC;

    /* Fetch existing bacRecord and update relevant fields */
    /* Append this record to the store and chain it to the currently
       active record */
    getGateKeeperLocaleStore(GkpStore.locStore, ibc_mesg->ibac_id, &bacRecord);
    if(assembler == AS_ASSEMBLER_OVERLAY){
      outputBacDelete(ibc_mesg->ibac_id,ibc_mesg->ebac_id);
      outputBactigFragDeletes(bacRecord);
      ibc_mesg->action = AS_ADD;
    }
    bacRecord.prevID = ibc_mesg->ibac_id;
    bacRecord.hasSequence = FALSE;
    bacRecord.prevInstanceID = 0;
    bacRecord.redefined = TRUE;
    bacRecord.deathBatch = batchID;  // Record when we were redefined
    appendGateKeeperLocaleStore(GkpStore.s_locStore,  &bacRecord);

    /* The newly revised bacRecord will reference the one we just stashed in
       the table */
    bacRecord.prevInstanceID = getNumGateKeeperLocales(GkpStore.s_locStore) - 1; // index in s_locStore
    if(verbose)
      fprintf(stderr,"* Appended record " F_IID " to s_locStore\n",
	    bacRecord.prevInstanceID);

    if(verbose)
      fprintf(stderr,"* Saved Locale " F_IID " in new entry " F_IID "\n",
	      ibc_mesg->ibac_id, bacRecord.prevInstanceID);
    bacRecord.prevID = 0;
    bacRecord.birthBatch = batchID;
    bacRecord.deathBatch = 0;  // Record when we were redefined
	
    ibc_mesg->ibac_id = ibc_mesg->ibac_id;
    bacRecord.type = bac_mesg->type;
    bacRecord.numBactigs = 0;
    bacRecord.firstBactig = 0;
    bacRecord.sequenceID = 0;
    if(bac_mesg->type == AS_FINISHED ||
       bac_mesg->type == AS_UNFINISHED){
      locValue.type = AS_IID_SEQ;
      if(HASH_SUCCESS != InsertInPHashTable_AS(&(GkpStore.hashTable),
                                               UID_NAMESPACE_AS,
                                               bac_mesg->eseq_id,
                                               &locValue, FALSE, TRUE))
         assert(0);

      ibc_mesg->iseq_id = locValue.IID;
      if(verbose)
	fprintf(stderr,"* Inserted seq " F_UID " in hashtable with IID " F_IID "\n",
		bac_mesg->eseq_id, ibc_mesg->iseq_id);
      seqRecord.deleted = FALSE;
      seqRecord.UID = bac_mesg->eseq_id;
      seqRecord.localeID = ibc_mesg->ibac_id;
      bacRecord.sequenceID = locValue.IID;
      if(verbose)
	fprintf(stderr,"* Sequence " F_IID " associated with bac " F_IID "\n",
		bacRecord.sequenceID, seqRecord.localeID);
    }
    AddRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, bac_mesg->elength);
    bacRecord.lengthID = ibc_mesg->ilength;
    {
      int i;
      for(i = 0; i < bac_mesg->num_bactigs; i++){
	InternalBactigMesg *bactig = bac_mesg->bactig_list + i;
	GateKeeperBactigRecord btgRecord;

	PHashValue_AS btgValue;

	btgValue.type = AS_IID_BTG;
	btgValue.IID = 0;

	if(HASH_SUCCESS != InsertInPHashTable_AS(&(GkpStore.hashTable),
                                                 UID_NAMESPACE_AS,
                                                 bactig->eaccession,
                                                 &btgValue, FALSE, TRUE))
          assert(0);

	bactig->iaccession = btgValue.IID;
	btgRecord.hasSequence = FALSE;
	btgRecord.deleted = FALSE;
	btgRecord.bacID = ibc_mesg->ibac_id;
	btgRecord.seqID = ibc_mesg->iseq_id;
	btgRecord.UID = bactig->eaccession;
	btgRecord.length = bactig->length;
	if(verbose)
	  fprintf(stderr,"* Inserting bactig " F_UID " in hashtable bactigID = " F_IID "  localeID=" F_IID "\n",
		  bactig->eaccession, bactig->iaccession, btgRecord.bacID);
	appendGateKeeperBactigStore(GkpStore.btgStore,&btgRecord);
	if(assembler == AS_ASSEMBLER_OVERLAY &&
	   bac_mesg->type == AS_UNFINISHED){
	  
	  Appendint32(addedBactigs, (int32 *)&(bactig->iaccession)); // Add to list of bactigs that were added
	}
      }
      bacRecord.firstBactig = (bac_mesg->bactig_list?bac_mesg->bactig_list[0].iaccession:0);
      bacRecord.numBactigs = bac_mesg->num_bactigs;
    }
    /*******vvvvvvv****** SPECIAL HANDLING *****************vvvvvvv**********/

    if(/* assembler = AS_ASSEMBLER_OVERLAY && */ bac_mesg->type == AS_FINISHED){
      static InternalBactigMesg dummy;
      GateKeeperBactigRecord btgRecord;

      // Create a dummy bactig 
      dummy.eaccession = 0;
      dummy.iaccession = AllocateCountPHashTable_AS(GkpStore.hashTable,AS_IID_BTG);
      dummy.length = 0;
      ibc_mesg->bactig_list = &dummy;
      ibc_mesg->num_bactigs = 1;
      bac_mesg->type = ibc_mesg->type = AS_UNFINISHED;
      bacRecord.numBactigs = 1;
      bacRecord.firstBactig = dummy.iaccession;
      btgRecord.hasSequence = FALSE;
      btgRecord.deleted = FALSE;
      btgRecord.bacID = ibc_mesg->ibac_id;
      btgRecord.seqID = ibc_mesg->iseq_id;
      btgRecord.UID = 0;
      btgRecord.length = 0;
      if(verbose)
	fprintf(stderr,"* Inserting *** DUMMY **** bactig " F_UID " in hashtable bactigID = " F_IID "  localeID=" F_IID "\n", 
		dummy.eaccession, dummy.iaccession, btgRecord.bacID);
      appendGateKeeperBactigStore(GkpStore.btgStore,&btgRecord);
      Appendint32(addedBactigs, (int32 *)&(dummy.iaccession)); // Add to list of bactigs that were added
    }
    /*******^^^^^^^^^******** SPECIAL HANDLING *************^^^^^^^^**********/

    setGateKeeperLocaleStore(GkpStore.locStore, ibc_mesg->ibac_id, &bacRecord);
    appendGateKeeperSequenceStore(GkpStore.seqStore, &seqRecord);
  }

#ifdef DEBUG_LOC
  fprintf(stderr,"* Inserted loc UID " F_UID " got iid " F_IID "\n",
	  frg_mesg->locale, ifg_mesg->ilocale);
#endif

   
  return GATEKEEPER_SUCCESS;

}



void outputBactigFragDeletes(GateKeeperLocaleRecord bacRecord){
  int i;
  InternalFragMesg ifg_mesg;
  GenericMesg pmesg;

  if(bacRecord.type != AS_UNFINISHED)
    return;

  ifg_mesg.action = AS_DELETE;
  pmesg.t = MESG_IFG;
  pmesg.m = &ifg_mesg;

  fprintf(stderr,"* Generating %d fragment delete messages to output for Bac " F_UID "'s bactigs\n",bacRecord.numBactigs,bacRecord.UID);
  for(i = 0; i < bacRecord.numBactigs; i++){
    GateKeeperBactigRecord btgRecord;
    getGateKeeperBactigStore(GkpStore.btgStore, i + bacRecord.firstBactig, &btgRecord);
    ifg_mesg.eaccession = btgRecord.UID;
    ifg_mesg.iaccession = i + bacRecord.firstBactig;
    Writer(Outfp, &pmesg);
  }
}

void outputBacDelete(CDS_IID_t bacIid, CDS_UID_t bacUid){
  InternalBacMesg ibc_mesg;
  GenericMesg pmesg;

  fprintf(stderr,"* Generating IBC delete messages to output for Bac " F_UID "\n",bacUid);
  pmesg.t = MESG_IBC;
  pmesg.m = &ibc_mesg;

  ibc_mesg.action = AS_DELETE;
  ibc_mesg.ebac_id = bacUid;
  ibc_mesg.ibac_id = bacIid;

  Writer(Ibcfp, &pmesg);
}
