
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
static char CM_ID[] = "$Id: AS_GKP_check.c,v 1.4 2005-03-22 19:48:52 jason_miller Exp $";

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



/***********************************************************************************/
int Check_DistanceMesg(DistanceMesg *dst_mesg,  
		       InternalDistMesg *idt_mesg,  
		       CDS_CID_t currentBatchID,
		       int verbose){
  /* Check validity of UID
     Check legal range of std/mean
     
     If OK, add/update a GateKeeperDistanceRecord
  */
	
  PHashValue_AS value;
  GateKeeperDistanceRecord gkpdst;

  Transfer_DST_to_IDT_AS(dst_mesg, idt_mesg);

  switch(dst_mesg->action){

  case AS_ADD:    /****************** ADD *************************/

    /* Make sure we haven't seen this distance record before... if so
       it is a fatal error */
    if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						 UID_NAMESPACE_AS,
						 dst_mesg->eaccession, 
						 AS_IID_DST, 
						 FALSE,
						 Msgfp,
						 &value)){

      printGKPError(Msgfp, GKPError_BadUniqueDST);
      fprintf(Msgfp,"# Check_DistanceMessage:  A message with UID " F_UID " exists %s!!! Can't reuse the UID... bye\n",
	      dst_mesg->eaccession, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }

    if(dst_mesg->mean <= 0 || dst_mesg->stddev <= 0 ||
       (dst_mesg->mean - 3.0 * dst_mesg->stddev < 0)){
      printGKPError(Msgfp, GKPError_DSTValues);
      fprintf(Msgfp,"# Check_DistanceMessage:  Illegal Mean %g and/or Standard Deviation %g\n",
	      dst_mesg->mean, dst_mesg->stddev);
      return(GATEKEEPER_FAILURE);
    }

    value.type = AS_IID_DST;
    InsertInPHashTable_AS(&(GkpStore.hashTable), UID_NAMESPACE_AS, dst_mesg->eaccession, &value, FALSE, TRUE);

    idt_mesg->iaccession = value.IID;
    gkpdst.UID = dst_mesg->eaccession;
    gkpdst.birthBatch = (uint16) currentBatchID;
    gkpdst.deathBatch = 0;
    gkpdst.deleted = FALSE;
    gkpdst.mean = idt_mesg->mean;
    gkpdst.stddev = idt_mesg->stddev;
    gkpdst.prevID = 0;
    gkpdst.prevInstanceID = 0;

    /*
    fprintf(stderr,"* Appended DST (" F_IID "," F_UID ") mean:%g std:%g\n",
            value.IID, dst_mesg->eaccession, gkpdst.mean, gkpdst.stddev);
    */

    appendIndexStore(GkpStore.dstStore, &gkpdst); // save the IID->UID mapping

    /* Also add a spot for the librayr to the libStore,
       to keep consistent IIDs/indices */
    {
      GateKeeperLibDonorRecord gkpldr;
      memset(&gkpldr, 0, sizeof(GateKeeperLibDonorRecord));
      appendGateKeeperLibDonorStore(GkpStore.libStore, &gkpldr);
    }
    break;

  case AS_REDEFINE:/****************** REDEFINE *************************/

    /* Make sure we HAVE seen this distance record before... if we haven't
       it is a fatal error */
    if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						 UID_NAMESPACE_AS,
						 dst_mesg->eaccession, 
						 AS_IID_DST, 
						 FALSE,
						 Msgfp,
						 &value)){

      printGKPError(Msgfp, GKPError_MissingDST);
      fprintf(Msgfp,"# Check_DistanceMessage:  Distance " F_UID " does NOT exist %s!!! Can't redefine it...\n",
	      dst_mesg->eaccession, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }

    if(dst_mesg->mean <= 0 || dst_mesg->stddev <= 0 ||
       (dst_mesg->mean - 3.0 * dst_mesg->stddev < 0)){
      printGKPError(Msgfp, GKPError_DSTValues);
      fprintf(Msgfp,"# Check_DistanceMessage:  Illegal Mean %g and/or Standard Deviation %g\n",
	      dst_mesg->mean, dst_mesg->stddev);
      return(GATEKEEPER_FAILURE);
    }

    idt_mesg->iaccession = value.IID;

    // Retrieve the record for this distance and update it

    getIndexStore(GkpStore.dstStore, value.IID, &gkpdst); 
    gkpdst.redefined = TRUE;
    gkpdst.prevID = value.IID;
    gkpdst.prevInstanceID = 0;
    gkpdst.deathBatch = (uint16) currentBatchID;

    appendGateKeeperDistanceStore(GkpStore.s_dstStore, &gkpdst); 

    gkpdst.prevID = 0;
    gkpdst.prevInstanceID = getNumGateKeeperDistances(GkpStore.s_dstStore) - 1; // index in s_dstStore

    // Update these two fields ONLY
    gkpdst.mean = idt_mesg->mean;
    gkpdst.stddev = idt_mesg->stddev;

    setGateKeeperDistanceStore(GkpStore.dstStore, value.IID, &gkpdst); // save the IID->UID mapping
    break;

  case AS_DELETE: /****************** Delete *************************/
    /* Make sure we HAVE seen this distance record before... if not
       it is a fatal error. */
    if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						 UID_NAMESPACE_AS,
						 dst_mesg->eaccession, 
						 AS_IID_DST, 
						 FALSE,
						 Msgfp,
						 &value)){

      printGKPError(Msgfp, GKPError_MissingDST);
      fprintf(Msgfp,"# Check_DistanceMessage:  Distance " F_UID " DOES NOT exist!!!" 
	      "Can't delete it... bye\n", dst_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }
    if(value.refCount > 1){
      printGKPError(Msgfp, GKPError_DeleteDST);
      fprintf(Msgfp,
	      "# Check_DistanceMessage: There are %d references outstanding to Distance "
	      F_UID ".\n"
	      "Can't delete it...\n",
	      value.refCount, dst_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    idt_mesg->iaccession = value.IID;
    deleteAndMarkGateKeeperDistanceStore(GkpStore.dstStore, value.IID, currentBatchID);
    DeleteFromPHashTable_AS(GkpStore.hashTable,UID_NAMESPACE_AS, dst_mesg->eaccession);
    break;

  default:
    fprintf(Msgfp,"# Check_DistanceMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }
  return GATEKEEPER_SUCCESS;
}
/***********************************************************************************************/
int Check_BatchMesg(BatchMesg *bat_mesg,  
		    InternalBatchMesg *iba_mesg,  
		    time_t currentTime,
		    int verbose){

  PHashValue_AS value;
  GateKeeperBatchRecord gkpbat;
  char *comment = (char *)malloc(strlen(bat_mesg->comment) + 2048);
  /* 
     Check validity of UID
     Check validity of Entrytime

     Add a GateKeeperBatchRecord
  */

  *iba_mesg = *bat_mesg;
  sprintf(comment,"Input File is %s\n%s",
	  Input_File_Name,
	  bat_mesg->comment);
  //  fprintf(stderr,"* in comment = %s\n", iba_mesg->comment);
  fprintf(stderr,"* comment = %s\n", comment);
  iba_mesg->comment = comment;


  /* Make sure we haven't seen this distance record before... if so
     it is a fatal error */
  if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       bat_mesg->eaccession, 
					       AS_IID_BAT, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp, GKPError_BadUniqueBAT);
    fprintf(Msgfp,"# Check_BatchMessage:  A message with UID " F_UID " exists %s!!! Can't reuse it... bye\n",
	    bat_mesg->eaccession, (value.deleted?"and has been deleted":""));
    return(GATEKEEPER_FAILURE);
  }

  /* Check entryTime */
  if(bat_mesg->created > currentTime){
    printGKPError(Msgfp, GKPError_Time);
    fprintf(Msgfp,"# Check_BatchMessage: invalid entry time " F_TIME_T " > current time (" F_TIME_T ")\n",
	    bat_mesg->created, currentTime);
    return GATEKEEPER_FAILURE;
  }


  value.type = AS_IID_BAT;
  InsertInPHashTable_AS(&(GkpStore.hashTable), UID_NAMESPACE_AS, bat_mesg->eaccession, &value, FALSE, TRUE);

  iba_mesg->iaccession = value.IID;

  gkpbat.UID = bat_mesg->eaccession;
  gkpbat.created = bat_mesg->created;
  strncpy(gkpbat.name, iba_mesg->name, 255);
  strncpy(gkpbat.comment, iba_mesg->comment, 255);
  gkpbat.numFragments = getNumGateKeeperFragments(GkpStore.frgStore);
  gkpbat.numLocales = getNumGateKeeperLocales(GkpStore.locStore);
  gkpbat.num_s_Locales = getNumGateKeeperLocales(GkpStore.s_locStore);
  gkpbat.numBactigs = getNumGateKeeperBactigs(GkpStore.btgStore);
  gkpbat.numDistances = getNumGateKeeperDistances(GkpStore.dstStore);
  gkpbat.num_s_Distances = getNumGateKeeperDistances(GkpStore.s_dstStore);
  gkpbat.numScreens = getNumGateKeeperScreens(GkpStore.scnStore);
  gkpbat.numRepeats = getNumGateKeeperRepeats(GkpStore.rptStore);
  gkpbat.numPlates = getNumGateKeeperSequencePlates(GkpStore.sqpStore);
  gkpbat.numWells = getNumGateKeeperWells(GkpStore.welStore);
  gkpbat.numLinks = getNumGateKeeperLinks(GkpStore.lnkStore);
  gkpbat.numSequences = getNumGateKeeperSequences(GkpStore.seqStore);

  if(verbose)
    fprintf(stderr,"* Batch " F_IID "  name:%s  comment:%s\n",
	  iba_mesg->iaccession, gkpbat.name, gkpbat.comment);

  appendGateKeeperBatchStore(GkpStore.batStore, &gkpbat); // save the IID->UID mapping

  return GATEKEEPER_SUCCESS;
}





/************************************************************************************************/

int Check_LinkMesg(LinkMesg *lkg_mesg,
		   InternalLinkMesg *ilk_mesg,
		   CDS_CID_t batchID,
  	           time_t currentTime,
		   int verbose,
                   int matchBAC_DstsInLkgMsgs){

  CDS_IID_t frag1IID, frag2IID;
  CDS_IID_t distIID;
  PHashValue_AS value;
  GateKeeperLinkRecord newLink;

  Transfer_LKG_to_ILK_AS(lkg_mesg, ilk_mesg);

#ifdef DEBUG_GKP
  fprintf(Msgfp,"* Link Message " F_UID " " F_UID "\n",
	  lkg_mesg->frag1, lkg_mesg->frag2);
#endif
  /* First check that the two referenced frags are distinct */
  if(lkg_mesg->frag1 == lkg_mesg->frag2){

    fprintf(Msgfp,"# Link Message:  Can't make a link from fragment " F_UID " to itself!!!" ,
	    lkg_mesg->frag1);
    return(GATEKEEPER_FAILURE);
  }

  /* Check that the two referenced fragments are alive and well */
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       lkg_mesg->frag1, 
					       AS_IID_FRAG, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp, GKPError_MissingFRG);
    fprintf(Msgfp,"# Link Message:  Fragment " F_UID " DOES NOT exist!!!" 
	    "Can't reference it... bye\n",	    lkg_mesg->frag1);
    return(GATEKEEPER_FAILURE);
  }
  ilk_mesg->ifrag1 = frag1IID = value.IID;
#if 0
  fprintf(stderr,"* Frag " F_UID " has %d references\n",
	  lkg_mesg->frag1, value.refCount);
#endif
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       lkg_mesg->frag2, 
					       AS_IID_FRAG, 
					       FALSE,
					       Msgfp,
					       &value)){

    printGKPError(Msgfp, GKPError_MissingFRG);
    fprintf(Msgfp,"# Link Message:  Fragment " F_UID " DOES NOT exist!!!" 
	    "Can't reference it... bye\n",	    lkg_mesg->frag2);
    return(GATEKEEPER_FAILURE);
  }
#if 0
  fprintf(stderr,"* Frag " F_UID " has %d references\n",
	  lkg_mesg->frag2, value.refCount);
#endif
  ilk_mesg->ifrag2 = frag2IID = value.IID;
  
    /* Make a canonical order for the link records (frag1 < frag2), to simplify
       some operations */
  if(frag1IID > frag2IID){
    CDS_IID_t tmp = frag1IID;
    frag1IID = frag2IID;
    frag2IID = tmp;
  }

  newLink.frag1 = frag1IID;
  newLink.frag2 = frag2IID;
  newLink.deleted = 0;
  newLink.frag1Next = 0;
  newLink.frag2Next = 0;

  switch(lkg_mesg->action){
  case AS_ADD:
    {
      GateKeeperFragmentRecord gkFrag1;
      GateKeeperFragmentRecord gkFrag2;

#ifdef DEBUG_GKP
      fprintf(Msgfp,"**** Checking a Link ADD\n");
#endif
      newLink.birthBatch = (uint16) batchID;
      newLink.deathBatch = 0;

      /* Check entryTime */
      if(lkg_mesg->entry_time > currentTime){
	printGKPError(Msgfp, GKPError_Time);
	fprintf(Msgfp,"# Check_LinkMessage: invalid entry time " F_TIME_T " > current time (" F_TIME_T ")\n",
		lkg_mesg->entry_time, currentTime);
	return GATEKEEPER_FAILURE;
      }
      /* Get the GateKeeperFragmentRecord for frag1 */
      getGateKeeperFragmentStore(GkpStore.frgStore, frag1IID, &gkFrag1);
      /* Get the GateKeeperFragmentRecord for frag2 */
      getGateKeeperFragmentStore(GkpStore.frgStore, frag2IID, &gkFrag2);

#ifdef DEBUG_GKP
      fprintf(stderr,"* frag1 = " F_IID " linkHead = " F_IID "   frag2 = " F_IID "  linkHead = " F_IID "\n",
	      frag1IID, gkFrag1.linkHead,
	      frag2IID, gkFrag2.linkHead);
#endif

      /* For Mate, Guide, and Join links, check that the distance record
	 is defined */
      distIID = 0;
      newLink.distance = 0;
      ilk_mesg->idistance = 0;

      if(lkg_mesg->type != AS_REREAD){
	if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						     UID_NAMESPACE_AS,
						     lkg_mesg->distance, 
						     AS_IID_DST, 
						     FALSE,
						     Msgfp,
						     &value)){

	  printGKPError(Msgfp, GKPError_MissingDST);
	  fprintf(Msgfp,"# Link Message:  Distance " F_UID " DOES NOT exist!!!" 
		  "Can't add link... bye\n",	    lkg_mesg->distance);
	  return(GATEKEEPER_FAILURE);
	}
	// newLink.distance = lkg_mesg->distance;
	// CHANGED by Knut Reinert 02/02/00
	newLink.distance    = value.IID;
	ilk_mesg->idistance = value.IID;

      }

      switch(lkg_mesg->type){
      case AS_MATE:
	  {
	/* Check if it has a mate -- if so squawk */
#ifdef DEBUG_GKP
	fprintf(Msgfp," Mate\n");
#endif
	newLink.type = lkg_mesg->type;

	if(gkFrag1.type == AS_READ && gkFrag2.type == AS_READ &&
           (char)lkg_mesg->link_orient != AS_INNIE){
	  fprintf(Msgfp,"# Read Mate Links must have orientation 'I' for AS_INNIE\n");
	  return GATEKEEPER_FAILURE;
	}

        if(gkFrag1.type == AS_TRNR && gkFrag2.type == AS_TRNR &&
           (char)lkg_mesg->link_orient != AS_OUTTIE){
	  fprintf(Msgfp,"# Transposon Read Mate Links must have orientation 'O' for AS_OUTTIE\n");
	  return GATEKEEPER_FAILURE;
	}

	newLink.orientation = lkg_mesg->link_orient;


	/* First check that both frags are READS or UBACs */
	if((gkFrag1.type != AS_READ || gkFrag2.type != AS_READ) &&
           (gkFrag1.type != AS_EXTR || gkFrag2.type != AS_EXTR) &&
           (gkFrag1.type != AS_TRNR || gkFrag2.type != AS_TRNR) &&
           (gkFrag1.type != AS_LBAC || gkFrag2.type != AS_LBAC) &&
           (gkFrag1.type != AS_UBAC || gkFrag2.type != AS_UBAC)){
	  printGKPError(Msgfp, GKPError_LNKFragtypeMismatch);
	  fprintf(Msgfp,
		  "# Fragments " F_UID " (type = %d) and " F_UID " (type = %d) are inconsistent with a Mate link\n",
		  gkFrag1.readUID, gkFrag1.type, gkFrag2.readUID, gkFrag2.type);
	  return GATEKEEPER_FAILURE;
	}

	break;
	}
	case AS_B_MATE:
	{
	/* Check if it has a mate -- if so squawk */
#ifdef DEBUG_GKP
	fprintf(Msgfp," Mate\n");
#endif
	newLink.type = lkg_mesg->type;

	if(gkFrag1.type == AS_B_READ && gkFrag2.type == AS_B_READ &&
           (char)lkg_mesg->link_orient != AS_OUTTIE){
	  fprintf(Msgfp,"# Read BGLII Links must have orientation 'O' for AS_OUTTIE\n");
	  return GATEKEEPER_FAILURE;
	}
	/* First check that both frags are BGLII */
	if((gkFrag1.type != AS_B_READ || gkFrag2.type != AS_B_READ))
	{
	  printGKPError(Msgfp, GKPError_LNKFragtypeMismatch);
	  fprintf(Msgfp,
		  "# Fragments " F_UID " (type = %d) and " F_UID " (type = %d) are inconsistent with a BGLII link\n",
		  gkFrag1.readUID, gkFrag1.type, gkFrag2.readUID, gkFrag2.type);
	  return GATEKEEPER_FAILURE;
	 }
	 break;
	}
      case AS_BAC_GUIDE:
#ifdef DEBUG_GKP
	fprintf(Msgfp," BAC_GUIDE\n");
#endif

	newLink.type = AS_BAC_GUIDE;
	newLink.orientation = AS_GKP_INNIE;
#if 0
	/* For now, we simply ignore the orientation on these */
	if((char)lkg_mesg->link_orient != AS_INNIE){
	  fprintf(Msgfp,"# BAC_GUIDE Links must have orientation 'I' for AS_INNIE\n");
	  return GATEKEEPER_FAILURE;
	}
#endif
	/* First check that both frags are READS */
	if(gkFrag1.type != AS_EBAC || 
	   gkFrag2.type != AS_EBAC){
	  printGKPError(Msgfp, GKPError_LNKFragtypeMismatch);
	  fprintf(Msgfp,
		  "# Fragments " F_UID " (type = %d) and " F_UID " (type = %d) are inconsistent with a Bac Guide link\n",
		  gkFrag1.readUID, gkFrag1.type, gkFrag2.readUID, gkFrag2.type);
	  return GATEKEEPER_FAILURE;
	}
	/* Second check that they are from the same locale*/
	if(gkFrag1.localeID != gkFrag2.localeID ){
	  printGKPError(Msgfp, GKPError_LNKLocaleMismatch);
	  fprintf(Msgfp,
		  "# Fragments " F_UID " (locale = " F_IID ") and " F_UID " (locale = " F_IID ") have locales inconsistent with a Bac Guide link\n",
		  gkFrag1.readUID, gkFrag1.localeID, gkFrag2.readUID, gkFrag2.localeID);
	  return GATEKEEPER_FAILURE;
	}
	/* Third check that the distance ID specified is the same as the lengthID specified for the BAC */
        if(matchBAC_DstsInLkgMsgs)
        {
	  GateKeeperLocaleRecord gkpl;
	  
	  getGateKeeperLocaleStore(GkpStore.locStore, gkFrag1.localeID, &gkpl);

	  if(gkpl.lengthID != newLink.distance ||
	     gkpl.lengthID != newLink.distance ){
	  printGKPError(Msgfp, GKPError_LNKLocaleDistanceMismatch);
	  fprintf(Msgfp,
		  "# Link specifies distance id " F_UID " that is inconsistent with its associated locale " F_UID "\n",
		 lkg_mesg->distance , gkFrag1.readUID);
	  return GATEKEEPER_FAILURE;
	  }
	     
	}
        // validity of distance UID is checked earlier in this function
	break;

      case AS_STS_GUIDE:
#ifdef DEBUG_GKP
	fprintf(Msgfp," STS Guide\n");
#endif
	if(gkFrag1.type != (unsigned int) AS_STS || 
	   gkFrag2.type != (unsigned int) AS_STS){
	  printGKPError(Msgfp, GKPError_LNKFragtypeMismatch);
	  fprintf(Msgfp,
		  "# Fragments " F_UID " (type = %d) and " F_UID " (type = %d) are inconsistent with an STS Guide link(%d)\n",
		  gkFrag1.readUID, gkFrag1.type, gkFrag2.readUID, gkFrag2.type, AS_STS);
	  return GATEKEEPER_FAILURE;
	}
	newLink.type = AS_STS_GUIDE;
	newLink.orientation = AS_GKP_UNKNOWN;
	newLink.distance = 0;

	break;

      case AS_REREAD:
#ifdef DEBUG_GKP
	fprintf(Msgfp," ReRead\n");
#endif
	newLink.type = AS_REREAD;
	newLink.orientation = AS_GKP_UNKNOWN;
	newLink.distance = 0;

	break;
        
      case AS_MAY_JOIN:
#ifdef DEBUG_GKP
	fprintf(Msgfp," MayJoin\n");
#endif
	newLink.type = AS_MAY_JOIN;
	newLink.orientation = AS_GKP_UNKNOWN;
	newLink.distance = 0;
        
        break;
        
      case AS_MUST_JOIN:
#ifdef DEBUG_GKP
	fprintf(Msgfp," MustJoin\n");
#endif
	newLink.type = AS_MUST_JOIN;
	newLink.orientation = AS_GKP_UNKNOWN;
	newLink.distance = 0;
        
        break;
        
      default:
	printGKPError(Msgfp, GKPError_Scalar);
	fprintf(Msgfp,"# Check_LinkMessage: invalid link type\n");
	return GATEKEEPER_FAILURE;
      }


      /* Find any links incompatible with the new link. */
      if(findLinksInCompatibleWith(GkpStore.lnkStore, frag1IID, gkFrag1.linkHead, &newLink, Msgfp, verbose)){
	return GATEKEEPER_FAILURE;
      }
      if(findLinksInCompatibleWith(GkpStore.lnkStore, frag2IID, gkFrag2.linkHead, &newLink, Msgfp, verbose)){
	return GATEKEEPER_FAILURE;
      }


      switch((char)lkg_mesg->link_orient){
      case AS_UNKNOWN:
	newLink.orientation = AS_GKP_UNKNOWN;
	break;
      case AS_INNIE:
	newLink.orientation = AS_GKP_INNIE;
	break;
      case AS_OUTTIE:
	newLink.orientation = AS_GKP_OUTTIE;
	break;
      case AS_NORMAL:
	newLink.orientation = AS_GKP_NORMAL;
	break;
      case AS_ANTI:
	newLink.orientation = AS_GKP_ANTINORMAL;
	break;
      default:
	printGKPError(Msgfp, GKPError_Scalar);
	fprintf(Msgfp,"# Check_LinkMessage: invalid link_orient %c\n",
		lkg_mesg->link_orient);
	return GATEKEEPER_FAILURE;
      }



      if(HASH_SUCCESS != AddRefPHashTable_AS(GkpStore.hashTable,
                                             UID_NAMESPACE_AS,
                                             lkg_mesg->distance) ||
         HASH_SUCCESS != AddRefPHashTable_AS(GkpStore.hashTable,
                                             UID_NAMESPACE_AS,
                                             lkg_mesg->frag1) ||
         HASH_SUCCESS != AddRefPHashTable_AS(GkpStore.hashTable,
                                             UID_NAMESPACE_AS,
                                             lkg_mesg->frag2))
        assert(0);


      linkLink_GKP(GkpStore.lnkStore, GkpStore.frgStore, &newLink, frag1IID, frag2IID, &gkFrag1, &gkFrag2);

    }
    break;
  case AS_DELETE:
    {
      GateKeeperFragmentRecord gkFrag1;
      GateKeeperFragmentRecord gkFrag2;
      int deleteLink1 = 0;

#ifdef DEBUG_GKP
      fprintf(Msgfp,"**** Checking a Link DELETE\n");
      fprintf(stderr,"*** Deleting link (" F_IID "," F_IID ")\n",
	      frag1IID, frag2IID);
#endif
      /* Get the GateKeeperFragmentRecord for frag1 */
      getGateKeeperFragmentStore(GkpStore.frgStore, frag1IID, &gkFrag1);
      /* Get the GateKeeperFragmentRecord for frag2 */
      getGateKeeperFragmentStore(GkpStore.frgStore, frag2IID, &gkFrag2);
      newLink.type = (char)lkg_mesg->type;
      switch(lkg_mesg->type){
      case AS_MATE:
	  case AS_B_MATE:
      case AS_BAC_GUIDE:
      case AS_STS_GUIDE:
      case AS_MAY_JOIN:
      case AS_MUST_JOIN:
      case AS_REREAD:
	break;
      default:
	fprintf(Msgfp,"# Link Message: invalid link type\n");
	return GATEKEEPER_FAILURE;
      }
      newLink.distance = 0;          /* Matches any */
      newLink.orientation = AS_GKP_UNKNOWN; /* Matches any */
      deleteLink1 = findLink(GkpStore.lnkStore, frag1IID, gkFrag1.linkHead, &newLink, &newLink);

      /*
      fprintf(stderr,"* Looking from frag " F_IID " head " F_IID "\n",
              frag1IID, gkFrag1.linkHead);
      */
      
      if(deleteLink1 == 0){
	fprintf(Msgfp,"# Link Mesg:  Can't find link of type %c (" F_UID ", " F_UID ").. can't delete it\n",
		lkg_mesg->type, lkg_mesg->frag1, lkg_mesg->frag2);
	return GATEKEEPER_FAILURE;
      }else{
	unlinkLink_GKP(GkpStore.lnkStore, GkpStore.frgStore, frag1IID, frag2IID, 
		       &gkFrag1, &gkFrag2, &newLink, deleteLink1);
	UnRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, newLink.distance);
	UnRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, gkFrag2.readUID);
	UnRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, gkFrag1.readUID);

      }


    }
    break;
  default:
    fprintf(Msgfp,"# LinkMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }    

  return GATEKEEPER_SUCCESS;

}

/***********************************************************************************/
int Check_RepeatItemMesg(RepeatItemMesg *rpt_mesg,
			 InternalRepeatItemMesg *irp_mesg,
			 CDS_CID_t currentBatchID, 
			 int verbose){

  PHashValue_AS value;
  GateKeeperRepeatRecord gkprpt;

  *irp_mesg = *rpt_mesg;

  if(rpt_mesg->length < 0){
    fprintf(Msgfp,"# RepeatItemMessage: invalid length " F_COORD "\n",
            rpt_mesg->length);
    return GATEKEEPER_FAILURE;
  }
  if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       rpt_mesg->erepeat_id, 
					       AS_IID_RPT, 
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# A message with UID " F_UID " already exists!!! Can't add it... \n",
	    rpt_mesg->erepeat_id);
    return(GATEKEEPER_FAILURE);
  }
  value.type = AS_IID_RPT;
  InsertInPHashTable_AS(&GkpStore.hashTable, UID_NAMESPACE_AS, rpt_mesg->erepeat_id, &value, FALSE, TRUE);
  irp_mesg->irepeat_id = value.IID;

  /*
    fprintf(stderr,"* RPT " F_UID " ==> " F_IID "\n",
            rpt_mesg->erepeat_id, irp_mesg->irepeat_id);
  */

  gkprpt.UID = irp_mesg->erepeat_id;
  gkprpt.deleted = FALSE;
  strncpy(gkprpt.which,irp_mesg->which,255);
  appendGateKeeperRepeatStore(GkpStore.rptStore, &gkprpt);
  return GATEKEEPER_SUCCESS;

}




/***********************************************************************************/

int Check_ScreenItemMesg(ScreenItemMesg *scn_mesg,  
			 InternalScreenItemMesg *isn_mesg,  
			 CDS_CID_t currentBatchID, 
			 int verbose){
  PHashValue_AS value;
  GateKeeperScreenRecord gkpscn;

  gkpscn.deleted = FALSE;
  gkpscn.spare = 0;
  gkpscn.UID = scn_mesg->eaccession;
  gkpscn.birthBatch = (uint16) currentBatchID;
  gkpscn.deathBatch = 0;

  Transfer_SCN_to_ISN_AS(scn_mesg, isn_mesg);

  switch(scn_mesg->action){
  case AS_ADD:

  if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       scn_mesg->eaccession, 
					       AS_IID_SCN, 
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# ScreenItem: A message with UID " F_UID " already exists!!! Can't add it... \n",
	    scn_mesg->eaccession);
    return(GATEKEEPER_FAILURE);
  }
  isn_mesg->iaccession = value.IID;
  break;

  case AS_DELETE:
  case AS_UPDATE:
  if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       scn_mesg->eaccession, 
					       AS_IID_SCN, 
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# ScreenItem " F_UID " does NOT exist!!! Can't %s it... \n",
	    scn_mesg->eaccession, (scn_mesg->action == AS_DELETE?"delete":"update"));
    return(GATEKEEPER_FAILURE);
  }
  break;

  default:
    fprintf(Msgfp,"# ScreenItem " F_UID " has an illegal action\n",
	    scn_mesg->eaccession);
    return GATEKEEPER_FAILURE;
  }


  if(scn_mesg->action == AS_DELETE){
    deleteAndMarkGateKeeperScreenStore(GkpStore.scnStore, value.IID, currentBatchID);
    DeleteFromPHashTable_AS(GkpStore.hashTable,UID_NAMESPACE_AS, scn_mesg->eaccession);
    return GATEKEEPER_SUCCESS;
  }


  if(scn_mesg->min_length < GATEKEEPER_SCREENER_MIN_LENGTH){
    fprintf(Msgfp,"# ScreenItemMessage: invalid min_length " F_COORD "\n", scn_mesg->min_length);
    return GATEKEEPER_FAILURE;
  }

  if(scn_mesg->variation > 1.0 ||
     scn_mesg->variation < 0.0){
    fprintf(Msgfp,"# Check_ScreenItemMessage: invalid variation %g\n", scn_mesg->variation);
    return GATEKEEPER_FAILURE;
  }
  switch(scn_mesg->type){
    
  case AS_UBIQREP:
    break;
  case AS_CONTAMINANT:
    break;
  default:
    fprintf(Msgfp,"# Check_ScreenItemMessage: invalid type %c\n", scn_mesg->type);
    return GATEKEEPER_FAILURE;
  }

  if(HASH_FAILURE == LookupTypeInPHashTable_AS(GkpStore.hashTable, 
					       UID_NAMESPACE_AS,
					       scn_mesg->erepeat_id, 
					       AS_IID_RPT,
					       FALSE,
					       Msgfp,
					       &value)){

    fprintf(Msgfp,"# ScreenItem references undefined repeat_id " F_UID "!!! Can't add it... \n",
	    scn_mesg->erepeat_id);
    return(GATEKEEPER_FAILURE);
  }
  isn_mesg->irepeat_id = value.IID;
  gkpscn.repeatID = value.IID;

  {
    char *s;
    int seqLength;
    if(GATEKEEPER_FAILURE == checkSequence(scn_mesg->sequence, &s, &seqLength)){
      fprintf(Msgfp,"# ScreenItem: Invalid char %c at position " F_SIZE_T " in sequence\n",
	      *s, s - scn_mesg->sequence);
      return GATEKEEPER_FAILURE;
    }
  }
	
  if(scn_mesg->action == AS_ADD){
    value.type = AS_IID_SCN;
    InsertInPHashTable_AS(&GkpStore.hashTable, UID_NAMESPACE_AS, scn_mesg->eaccession, &value, FALSE, TRUE);
    isn_mesg->iaccession = value.IID;
    AddRefPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, scn_mesg->erepeat_id);
    appendGateKeeperScreenStore(GkpStore.scnStore, &gkpscn);
    
    fprintf(stderr,
            "* Appended scn uid:" F_UID " iid:" F_IID " rptid:" F_IID "\n",
    	    gkpscn.UID, isn_mesg->iaccession, gkpscn.repeatID);
    
  }else{
    setGateKeeperScreenStore(GkpStore.scnStore,isn_mesg->iaccession,&gkpscn);
  }
  return GATEKEEPER_SUCCESS;

}


/******************************************************************************
 * Function: AddPlateWells
 * Inputs:
 *  PlateMesg * pla_mesg - pointer to plate message containing wells
 *  IntPlate_ID iplate   - internal ID of plate
 *  CDS_IID_t  * first_well - to be assigned within function
 * I/O
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
******************************************************************************/
static int AddPlateWells(PlateMesg * pla_mesg,
                         IntPlate_ID iplate,
                         CDS_IID_t * first_well){
  WellMesg * wel_mesg;
  GateKeeperWellRecord gkpwel;
  GateKeeperAuxFragRecord gkpafr;
  PHashValue_AS value;
  Well_ID well_id = 0;
  int i;
  
  gkpwel.deleted = FALSE;
  
  // check library & fragment for each well
  for(i = 0, wel_mesg = pla_mesg->well_list;
      i < pla_mesg->num_wells;
      i++, wel_mesg++){
    
    if(i > 0 && wel_mesg->ewell <= well_id){
      printGKPError(Msgfp, GKPError_BadUniqueWEL);
      fprintf(Msgfp,"# Well Message: A message with ID %u on plate " F_UID " is out of order. IDs must be in ascending order.\n",
              wel_mesg->ewell, pla_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }
    well_id = wel_mesg->ewell;
    
    // check that well number is okay
    if(wel_mesg->ewell > GATEKEEPER_MAX_WELL_NUMBER)
    {
      printGKPError(Msgfp, GKPError_WellNumberOutOfRange);
      fprintf(Msgfp,"# Well Message:  Well Number is out of range " F_U16 " (must be [0," F_U16 "]\n",
              wel_mesg->ewell, GATEKEEPER_MAX_WELL_NUMBER);
      return(GATEKEEPER_FAILURE);
    }
    
    // NOTE - this is a well number on the plate, not a UID
    gkpwel.ewell = wel_mesg->ewell;
    
    /* Check that library exists */
    if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                 UID_NAMESPACE_AS,
                                                 wel_mesg->elibrary, 
                                                 AS_IID_DST, 
                                                 FALSE,
                                                 Msgfp,
                                                 &value)){
      
      printGKPError(Msgfp, GKPError_MissingDST);
      fprintf(Msgfp,"# Plate well Message:  Distance " F_UID " DOES NOT exist %s!!! Can't add plate well!\n",
              wel_mesg->elibrary, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }
    gkpwel.ilib = value.IID;
    
    /* Check that the fragment exists */
    if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                 UID_NAMESPACE_AS,
                                                 wel_mesg->efrag, 
                                                 AS_IID_FRAG, 
                                                 FALSE,
                                                 Msgfp,
                                                 &value)){
      printGKPError(Msgfp, GKPError_MissingFRG);
      fprintf(Msgfp,"# Plate well Message:  Fragment " F_UID " DOES NOT exist %s!!! Can't add plate well!\n",
              wel_mesg->efrag, (value.deleted?"and has been deleted":""));
      return(GATEKEEPER_FAILURE);
    }
    gkpwel.ifrag = value.IID;

    /* check that the fragment type is acceptable */
    {
      GateKeeperFragmentRecord gkpfr;
      getGateKeeperFragmentStore(GkpStore.frgStore, gkpwel.ifrag, &gkpfr);
      if(gkpfr.type != AS_READ &&
         gkpfr.type != AS_EXTR &&
         gkpfr.type != AS_TRNR){
        fprintf(Msgfp,"# Plate well Message: Fragment " F_UID " type %c is not supported.\n",
                wel_mesg->efrag, gkpfr.type);
        return(GATEKEEPER_FAILURE);
      }
    }
    
    /* check that aux frag store doesn't already have a well for the frag */
    getGateKeeperAuxFragStore(GkpStore.auxStore, gkpwel.ifrag, &gkpafr);
    if(gkpafr.set != 0 || gkpafr.iwell != 0){
      printGKPError(Msgfp, GKPError_BadUniqueFRG);
      fprintf(Msgfp,"# Plate well Message:  Fragment " F_UID " is referenced by another well!!!" 
              "Can't add plate well!\n", wel_mesg->efrag);
      return(GATEKEEPER_FAILURE);
    }else{
      gkpafr.iplate = iplate;
      gkpafr.deleted = FALSE;
      gkpafr.set = TRUE;
      gkpafr.iwell = getNumGateKeeperWells(GkpStore.welStore) + 1;
      gkpafr.ilib = gkpwel.ilib;
      setGateKeeperAuxFragStore(GkpStore.auxStore, gkpwel.ifrag, &gkpafr);
    }

    appendGateKeeperWellStore(GkpStore.welStore, &gkpwel);
    if( i == 0 )
      *first_well = getNumGateKeeperWells(GkpStore.welStore);
  }
  return(GATEKEEPER_SUCCESS);
}

  
/******************************************************************************
 * Function: Check_PlateMesg
 * Inputs:
 * I/O
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
 *****************************************************************************/
int Check_PlateMesg(PlateMesg *pla_mesg,
                    CDS_CID_t batchID,
                    time_t currentTime,
                    int assembler,
                    int strict,
                    int verbose){
  GateKeeperSequencePlateRecord gkpsqp;
  PHashValue_AS value;

  switch(pla_mesg->action){
    case AS_ADD:
    {
      /* Check that plate doesn't exist */
      if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   pla_mesg->eaccession, 
                                                   AS_IID_PLA, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        
        printGKPError(Msgfp, GKPError_BadUniquePLA);
        fprintf(Msgfp,"# Plate Message:  A message with UID " F_UID " EXISTS already, can't define it!!!" ,
                pla_mesg->eaccession);
        return(GATEKEEPER_FAILURE);
      }

      memset(&gkpsqp, 0, sizeof(GateKeeperSequencePlateRecord));
      if(GATEKEEPER_SUCCESS !=
         AddPlateWells(pla_mesg,
                       getNumGateKeeperSequencePlates(GkpStore.sqpStore) + 1,
                       &(gkpsqp.firstWell)))
        return(GATEKEEPER_FAILURE);
      
      value.type = AS_IID_PLA;
      InsertInPHashTable_AS(&GkpStore.hashTable, UID_NAMESPACE_AS,
                            pla_mesg->eaccession, &value, FALSE, TRUE);
      
      /* Append record for the plate */
      gkpsqp.UID = pla_mesg->eaccession;
      gkpsqp.deleted = FALSE;
      gkpsqp.numWells = pla_mesg->num_wells;
      // gkpsqp.firstWell is set in loop over wells, above
      gkpsqp.mate = 0;
      gkpsqp.birthBatch = (uint16) batchID;
      appendGateKeeperSequencePlateStore(GkpStore.sqpStore, &gkpsqp);
    }
    break;
    case AS_DELETE:
    case AS_REDEFINE:
    {
      /* Checks & actions:
         1. The plate must exist.
         2. Get plate record
            if( AS_DELETE )
              delete it.
         3. Loop over (old) wells in well Store
            disassociate auxFrag from well
            delete well
         if( AS_REDEFINE )
           4. Loop over new wells
              look up fragments - must exist
              look up aux fragments - must not be set
              associate auxFrag with well
              add well
           5. Set new plate record firstWell & num_wells
      */
      GateKeeperSequencePlateRecord gkpsqp;
      IntWell_ID i;
      
      /* Check that plate exists */
      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   pla_mesg->eaccession, 
                                                   AS_IID_PLA, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        
        printGKPError(Msgfp, GKPError_MissingPLA);
        fprintf(Msgfp,"# Plate Message:  Plate UID " F_UID " DOES NOT exist %s! Can't delete/redefine it!!!" ,
                pla_mesg->eaccession,
                (value.deleted?"and has been deleted":""));
        return(GATEKEEPER_FAILURE);
      }
      getGateKeeperSequencePlateStore(GkpStore.sqpStore, value.IID, &gkpsqp);

      /* Disassociate aux frags from old wells
         and delete old wells
       */
      for(i = gkpsqp.firstWell; i < gkpsqp.firstWell + gkpsqp.numWells; i++){
        GateKeeperWellRecord gkpwel;

        getGateKeeperWellStore(GkpStore.welStore, i, &gkpwel);
        if(gkpwel.ifrag > 0){
          GateKeeperAuxFragRecord gkpafr;
          getGateKeeperAuxFragStore(GkpStore.auxStore, gkpwel.ifrag, &gkpafr);
          memset(&gkpafr, 0, sizeof(GateKeeperAuxFragRecord));
          setGateKeeperAuxFragStore(GkpStore.auxStore, gkpwel.ifrag, &gkpafr);
        }
        gkpwel.deleted = TRUE;
        setGateKeeperWellStore(GkpStore.welStore, i, &gkpwel);
      }
      
      if(pla_mesg->action == AS_DELETE){
        /* Link message must have been deleted first */
        if(gkpsqp.mate){
          fprintf(Msgfp,"# Plate Message:  Plate UID " F_UID " has mate. Link must be deleted before plate!\n",
                pla_mesg->eaccession);
          return(GATEKEEPER_FAILURE);
        }
      
        // delete from PHashTable & store
        deleteAndMarkGateKeeperSequencePlateStore(GkpStore.sqpStore, value.IID,batchID);
        DeleteFromPHashTable_AS(GkpStore.hashTable,UID_NAMESPACE_AS,
                                pla_mesg->eaccession);
      }else{
        // move the old plate message to the s_sqp store
        gkpsqp.prevID = value.IID;
        gkpsqp.prevInstanceID = 0;
        gkpsqp.redefined = TRUE;
        gkpsqp.deathBatch = (uint16) batchID;
        appendGateKeeperSequencePlateStore(GkpStore.s_sqpStore, &gkpsqp);

        gkpsqp.prevInstanceID = getNumGateKeeperSequencePlates(GkpStore.s_sqpStore) - 1;
        gkpsqp.prevID = 0;
        gkpsqp.birthBatch = (uint16) batchID;
        gkpsqp.deathBatch = 0;
        
        // add new wells & replace previous gatekeeper plate record
        gkpsqp.numWells = pla_mesg->num_wells;
        if(GATEKEEPER_SUCCESS !=
           AddPlateWells(pla_mesg, value.IID, &(gkpsqp.firstWell)))
          return(GATEKEEPER_FAILURE);

        // update the plate message
        setGateKeeperSequencePlateStore(GkpStore.sqpStore, value.IID, &gkpsqp);
      }
    }
    break;
    default:
    {
      printGKPError(Msgfp, GKPError_Action);
      fprintf(Msgfp,"# Plate Message: Plate " F_UID ". Action %c not supported\n", pla_mesg->eaccession, pla_mesg->action);
      return(GATEKEEPER_FAILURE);
    }
    break;
  }

  return(GATEKEEPER_SUCCESS);
}


/******************************************************************************
 * Function: Check_LinkPlateMesg
 * Description:
 * Checks on add:
 * Checks on delete:
 * Inputs:
 *
 * I/O
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
******************************************************************************/
int Check_LinkPlateMesg(LinkPlateMesg *lkp_mesg,
                        CDS_CID_t batchID,
                        time_t currentTime,
                        int assembler,
                        int strict,
                        int verbose){
  GateKeeperSequencePlateRecord gkpsqp_for;
  GateKeeperSequencePlateRecord gkpsqp_rev;
  PHashValue_AS value;
  IntPlate_ID iplate_for;
  IntPlate_ID iplate_rev;

  if(lkp_mesg->action != AS_ADD &&
     lkp_mesg->action != AS_DELETE){
    printGKPError(Msgfp, GKPError_Action);
    fprintf(Msgfp,"# Link Plate Mesg: action %c not supported.\n",
            lkp_mesg->action);
    return(GATEKEEPER_FAILURE);
  }

  /* Check that both plates exist */
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                               UID_NAMESPACE_AS,
                                               lkp_mesg->eplate_for, 
                                               AS_IID_PLA, 
                                               FALSE,
                                               Msgfp,
                                               &value)){
    
    printGKPError(Msgfp, GKPError_MissingPLA);
    fprintf(Msgfp,"# Link Plate Mesg:  Plate " F_UID " DOES NOT exist %s!!! Can't add link plate!\n",
            lkp_mesg->eplate_for,(value.deleted?"and has been deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  iplate_for = value.IID;
  
  if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                               UID_NAMESPACE_AS,
                                               lkp_mesg->eplate_rev, 
                                               AS_IID_PLA, 
                                               FALSE,
                                               Msgfp,
                                               &value)){
    
    printGKPError(Msgfp, GKPError_MissingPLA);
    fprintf(Msgfp,"# Link Plate Mesg:  Plate " F_UID " DOES NOT exist %s!!! Can't add link plate!\n",
            lkp_mesg->eplate_rev,(value.deleted?"and has been deleted":""));
    return(GATEKEEPER_FAILURE);
  }
  iplate_rev = value.IID;
      
  getGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_for, &gkpsqp_for);
  getGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_rev, &gkpsqp_rev);

  if(lkp_mesg->action == AS_ADD){
    /* check that neither plate has a mate yet */
    if(gkpsqp_for.mate != 0){
      printGKPError(Msgfp, GKPError_BadUniquePLA);
      fprintf(Msgfp,"# Link Plate Mesg: Forward plate " F_UID " ALREADY has reverse mate\n",
              lkp_mesg->eplate_for);
      return(GATEKEEPER_FAILURE);
    }else
    {
      gkpsqp_for.mate = iplate_rev;
      setGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_for, &gkpsqp_for);
    }
    
    if(gkpsqp_rev.mate != 0){
      printGKPError(Msgfp, GKPError_BadUniquePLA);
      fprintf(Msgfp,"# Link Plate Mesg: Reverse plate " F_UID " ALREADY has forward mate\n",
              lkp_mesg->eplate_rev);
      return(GATEKEEPER_FAILURE);
    }else
    {
      gkpsqp_rev.mate = iplate_for;
      setGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_rev, &gkpsqp_rev);
    }
  }else{
    /* check that both plates have a mate */
    if(gkpsqp_for.mate != iplate_rev){
      fprintf(Msgfp,"# Link Plate Mesg:  Can't find link between " F_UID " and " F_UID ". Can't delete it\n",
              lkp_mesg->eplate_for, lkp_mesg->eplate_rev);
      return(GATEKEEPER_FAILURE);
    }else{
      gkpsqp_for.mate = 0;
      setGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_for, &gkpsqp_for);
    }

    if(gkpsqp_rev.mate != iplate_for){
      fprintf(Msgfp,"# Link Plate Mesg:  Can't find link between " F_UID " and " F_UID ". Can't delete it\n",
              lkp_mesg->eplate_for, lkp_mesg->eplate_rev);
      return(GATEKEEPER_FAILURE);
    }else{
      gkpsqp_rev.mate = 0;
      setGateKeeperSequencePlateStore(GkpStore.sqpStore, iplate_rev, &gkpsqp_rev);
    }
  }
  
  return(GATEKEEPER_SUCCESS);
}


/******************************************************************************
 * Function: Check_LibDonorMesg
 * Description:
 * Checks on add:
 * Checks on delete:
 * Inputs:
 *
 * I/O
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
******************************************************************************/
int Check_LibDonorMesg(LibDonorMesg *lib_mesg,  
                       CDS_CID_t currentBatchID,
                       time_t currentTime,
                       int assembler,
                       int strict,
                       int verbose){
  PHashValue_AS value;
  IntLibrary_ID lib_iid;
  GateKeeperLibDonorRecord gkpldr;

  switch(lib_mesg->action){
    case AS_ADD:
    {
      /* check that library (distance) is present in distance store */
      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   lib_mesg->eaccession,
                                                   AS_IID_DST, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        printGKPError(Msgfp, GKPError_MissingDST);
        fprintf(Msgfp,"# LibDonor Mesg:  Distance " F_UID " does NOT exist %s!!! Can't associate it with a donor...\n",
                lib_mesg->eaccession,(value.deleted?"and has been deleted":""));
        return(GATEKEEPER_FAILURE);
      }
      lib_iid = value.IID;
      
      /* check that the library is not already set in the library store */
      getGateKeeperLibDonorStore(GkpStore.libStore, lib_iid, &gkpldr);
      if(gkpldr.set != 0 || gkpldr.idonor != 0){
        printGKPError(Msgfp, GKPError_BadUniqueDST);
        fprintf(Msgfp,"# LibDonor Mesg: A message with UID " F_UID " exists!!! Can't reuse the UID!\n",
                lib_mesg->eaccession);
        return(GATEKEEPER_FAILURE);
      }
      
      /* check whether the donor is already present */
      /* if not, add it */
      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   lib_mesg->donor,
                                                   AS_IID_DON, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        GateKeeperDonorRecord gkpdon;
        value.type = AS_IID_DON;
        InsertInPHashTable_AS(&(GkpStore.hashTable), UID_NAMESPACE_AS, lib_mesg->donor, &value, FALSE, TRUE);
        gkpdon.deleted = 0;
        gkpdon.UID = lib_mesg->donor;
        appendIndexStore(GkpStore.donStore, &gkpdon);
      }
      
      gkpldr.idonor = value.IID;
      gkpldr.set = 1;
      strcpy(gkpldr.source, lib_mesg->source);
      setGateKeeperLibDonorStore(GkpStore.libStore, lib_iid, &gkpldr);
    }
    break;
    case AS_DELETE:
    {
      // check that library (distance) exists
      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   lib_mesg->eaccession,
                                                   AS_IID_DST, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        printGKPError(Msgfp, GKPError_MissingDST);
        fprintf(Msgfp,"# LibDonor Mesg:  Distance " F_UID " does NOT exist %s!!! Can't associate it with a donor.\n",
                lib_mesg->eaccession,(value.deleted?"and has been deleted":""));
        return(GATEKEEPER_FAILURE);
      }
      lib_iid = value.IID;
      
      // check that the library is set in the library store
      getGateKeeperLibDonorStore(GkpStore.libStore, lib_iid, &gkpldr);
      if(gkpldr.set == 0 || gkpldr.idonor == 0){
        fprintf(Msgfp,"# LibDonor Mesg: Library UID " F_UID " is not associated with a donor! Can't delete.\n",
                lib_mesg->eaccession);
        return(GATEKEEPER_FAILURE);
      }

      gkpldr.idonor = 0;
      // leave set field as is
      setGateKeeperLibDonorStore(GkpStore.libStore, lib_iid, &gkpldr);
    }
    break;
    default:
    {
      printGKPError(Msgfp, GKPError_Action);
      fprintf(Msgfp,"# LibDonor Mesg: Library UID " F_UID " action %c is not supported.\n", lib_mesg->eaccession, lib_mesg->action);
      return(GATEKEEPER_FAILURE);
    }
    break;
  }
  
  return(GATEKEEPER_SUCCESS);
}


/***********************************************************************************/
FILE *  File_Open
(const char * Filename, const char * Mode, int exitOnFailure)

     /* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

{
  FILE  *  fp;

  fp = fopen (Filename, Mode);
  if  (fp == NULL && exitOnFailure)
    {
      fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
      exit (EXIT_FAILURE);
    }

  return  fp;
}


/***********************************************************************************/
/* Returns link index of incompatible link */
int findLinksInCompatibleWith(GateKeeperLinkStore store, 
			      CDS_IID_t followFragIID,
			      CDS_IID_t linkHead, 
			      GateKeeperLinkRecord *newlink,
			      FILE *Msgfp,
			      int verbose){
  CDS_IID_t frag1IID = newlink->frag1;
  CDS_IID_t frag2IID = newlink->frag2;
  int linktype = newlink->type;

  GateKeeperLinkRecordIterator iterator;
  GateKeeperLinkRecord link;

#ifdef DEBUG_GKP
  fprintf(stderr,"*FindIncompatible from link " F_IID " with link (" F_IID "," F_IID ",t%d)\n",
	  linkHead, frag1IID, frag2IID, linktype);
#endif
  if(linkHead == NULL_LINK)
    return NULL_LINK;

  CreateGateKeeperLinkRecordIterator(store, linkHead,
				     followFragIID, &iterator);

  while(NextGateKeeperLinkRecordIterator(&iterator, &link)){

    if(link.deleted){
      continue;
    }
#ifdef DEBUG_GKP
    fprintf(stderr,"* link found (" F_IID "," F_IID ") type = %d\n",
	    link.frag1, link.frag2, link.type);
#endif
    switch(linktype){  // Branch depending on the linktype we are looking for

      /************* Looking for a Join ***********/

    case AS_MAY_JOIN:
    case AS_MUST_JOIN:
      switch(link.type){
      case AS_MATE:
      case AS_B_MATE:
      case AS_BAC_GUIDE:
      case AS_STS_GUIDE:
	continue;

      case AS_REREAD:
	fprintf(Msgfp,"# Found a previous reread link between these two fragments\n");
	return(iterator.prevLinkRecord);
	break;
      case AS_MAY_JOIN:
      case AS_MUST_JOIN:
	if(link.frag1 != frag1IID || 
	   link.frag2 != frag2IID)
	  continue;
	fprintf(Msgfp,"# Found a previous join link between these two fragments\n");
	return(iterator.prevLinkRecord);
      default:
	assert(0);
      }


      /************* Looking for a Mate ***********/
	 
    case AS_MATE:
    case AS_B_MATE:
      switch(link.type){
      case AS_MAY_JOIN:
      case AS_MUST_JOIN:
	if(link.frag1 != frag1IID || 
	   link.frag2 != frag2IID)
	  continue;
	fprintf(Msgfp,"# Warning: a %s join link exists between (" F_IID "," F_IID ")\n",
		(link.type == AS_MAY_JOIN?"MAY":"MUST"),
		frag1IID, frag2IID);
	break;
      case AS_MATE:
      case AS_B_MATE:
	/* Any mate link is a no no */
	fprintf(Msgfp,"# Found a previous mate link between these two fragments\n");
	return(iterator.prevLinkRecord);
      case AS_REREAD:
	/* Re reads with the same pair of fragIDs */
	if(link.frag1 != frag1IID || 
	   link.frag2 != frag2IID)
	  continue;
	fprintf(Msgfp,"# Found a previous reread link between these two fragments\n");
	return(iterator.prevLinkRecord);
      case AS_BAC_GUIDE:
      case AS_STS_GUIDE:
	/* By definition, we should have caught this before */
	/* Fragments with mates should not have guides */
	/* ABORT!!! *** INTENTIONAL FALLTHROUGH */
      default:
	assert(0);
	break;
      }
      /************* Looking for a Guide ***********/
    case AS_STS_GUIDE:
    case AS_BAC_GUIDE:
      switch(link.type){
      case AS_MAY_JOIN:  /* Ignore */
      case AS_MUST_JOIN:  /* Ignore */
	fprintf(Msgfp,"# Warning: a %s join link exists between (" F_IID "," F_IID ")\n",
		(link.type == AS_MAY_JOIN?"MAY":"MUST"),
		frag1IID, frag2IID);
	break;
      case AS_STS_GUIDE:
      case AS_BAC_GUIDE:
      case AS_REREAD:
	/* Guides with same fragIDs, with different distance/orientation */
	if(link.frag1 != frag1IID && link.frag2 != frag2IID)
	  continue;
	fprintf(Msgfp,"# Found a previous link between these two fragments\n");
	return(iterator.prevLinkRecord);

      case AS_MATE:
      case AS_B_MATE:
	/* By definition, we should have caught this before */
	/* Fragments with mates should not have guides */
	/* ABORT!!! *** INTENTIONAL FALLTHROUGH */
      default:
	assert(0);
	break;
      }
      break;
      /************* Looking for a Reread ***********/
    case AS_REREAD:
      if(link.frag1 == frag1IID && link.frag2 == frag2IID){
	fprintf(Msgfp,"# Found a previous  link between these two fragments\n");
	return(iterator.prevLinkRecord);
      }
      continue;

    default:
      assert(0);
    }
    
  }
  return NULL_LINK;


}






/*****************/
int dumpFrag_GKP(GateKeeperLinkStore gkplStore, 
		 GateKeeperFragmentStore gkpStore, 
		 CDS_IID_t frag1IID){
  GateKeeperFragmentRecord gkFrag1;
  GateKeeperLinkRecord link;
  GateKeeperLinkRecordIterator iterator;

  fprintf(stderr,"***** DumpFrag " F_IID " ********\n", frag1IID);
  getGateKeeperFragmentStore(gkpStore, frag1IID, &gkFrag1);
  fprintf(stderr,"* IID " F_IID " UID " F_UID " linkHead " F_IID "\n",
	  frag1IID, gkFrag1.readUID, gkFrag1.linkHead);

  if(gkFrag1.linkHead != 0){
    fprintf(stderr,"* Has the following links:\n");
    CreateGateKeeperLinkRecordIterator(gkplStore, gkFrag1.linkHead,
				       frag1IID, &iterator);

    while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
      fprintf(stderr,"\t* link " F_IID " (" F_IID "," F_IID ") next = " F_IID " type = %d %s\n",
	      iterator.prevLinkRecord, link.frag1, link.frag2, iterator.linkRecord,
	      link.type, (link.deleted?"DELETED":""));
    }
  }
  return TRUE;
}


/* Verify that the link with id linkID is alive and well and
   on the lists of both of its fragments */
int verifyLink_GKP(GateKeeperLinkStore gkplStore, 
		   GateKeeperFragmentStore     gkpStore, 
		   CDS_IID_t linkID){
  GateKeeperFragmentRecord gkFrag1, gkFrag2;
  GateKeeperLinkRecord link;
  GateKeeperLinkRecord slink;
  GateKeeperLinkRecordIterator iterator;
  CDS_IID_t foundLink1 = 0, foundLink2 = 0;

  getGateKeeperLinkStore(gkplStore, linkID, &link);

#ifdef DEBUG_GKP_VERBOSE
  fprintf(stderr,"* Verify Link " F_IID " (" F_IID "," F_IID ") type %d next (" F_IID "," F_IID ")\n",
	  linkID, link.frag1, link.frag2, link.type, link.frag1Next, link.frag2Next);

  dumpFrag_GKP( gkplStore, gkpStore, link.frag1);
  dumpFrag_GKP( gkplStore, gkpStore, link.frag2);
#endif
  getGateKeeperFragmentStore(gkpStore, link.frag1, &gkFrag1);
  getGateKeeperFragmentStore(gkpStore, link.frag2, &gkFrag2);

  if(gkFrag1.linkHead == linkID){
    foundLink1 = linkID;
  }else{
    CreateGateKeeperLinkRecordIterator(gkplStore, gkFrag1.linkHead,
				       link.frag1, &iterator);

    while(NextGateKeeperLinkRecordIterator(&iterator, &slink)){
      if(iterator.linkRecord == linkID)
	foundLink1 = linkID;
    }

  }

  if(gkFrag2.linkHead == linkID){
    foundLink2 = linkID;
  }else{
    CreateGateKeeperLinkRecordIterator(gkplStore, gkFrag2.linkHead,
				       link.frag2, &iterator);

    while(NextGateKeeperLinkRecordIterator(&iterator, &slink)){
      if(iterator.linkRecord == linkID)
	foundLink2 = linkID;
    }

  }

  assert(foundLink2 == foundLink1 && foundLink1 == linkID);

  return TRUE;
}



