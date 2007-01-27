
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
static char CM_ID[] = "$Id: AS_GKP_check.c,v 1.8 2007-01-27 00:30:10 brianwalenz Exp $";

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

  memset(&gkpdst, 0, sizeof(gkpdst));

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

    appendIndexStore(GkpStore.dstStore, &gkpdst); // save the IID->UID mapping

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

  memset(&gkpbat, 0, sizeof(gkpbat));

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

  memset(&newLink, 0, sizeof(newLink));

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

      memset(&gkFrag1, 0, sizeof(gkFrag1));
      memset(&gkFrag2, 0, sizeof(gkFrag2));

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
        /* This check is that both fragment are of the same type and can have mate.  VR */
        if((gkFrag1.type != gkFrag2.type) || !AS_FA_CAN_HAVE_MATE(gkFrag1.type)) {
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
	  memset(&gkpl, 0, sizeof(gkpl));
	  
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

      memset(&gkFrag1, 0, sizeof(gkFrag1));
      memset(&gkFrag2, 0, sizeof(gkFrag2));

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

int Check_ScreenItemMesg(ScreenItemMesg *scn_mesg,  
			 InternalScreenItemMesg *isn_mesg,  
			 CDS_CID_t currentBatchID, 
			 int verbose){
  PHashValue_AS value;
  GateKeeperScreenRecord gkpscn;
  
  memset(&gkpscn, 0, sizeof(gkpscn));
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

  memset(&link, 0, sizeof(link));

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



