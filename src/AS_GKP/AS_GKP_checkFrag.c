
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
static char CM_ID[] = "$Id: AS_GKP_checkFrag.c,v 1.13 2007-02-09 21:17:40 brianwalenz Exp $";

//#define DEBUG_GKP 1
//#define DEBUG_GKP_VERBOSE 1
/*************************************************
* Module:  AS_GKP_check.c
* Description:
*    Gatekeeper check routines.
* 
*    Programmer:  S. Kravitz
*       Written:  Jan 1999
*       Revised   Mar 2000 
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

// Map of quality value to fraction error
double qualityToFractionError[61];

/************************************************************************************************/
int Check_FragMesg(FragMesg *frg_mesg,  
		   InternalFragMesg *ifg_mesg,
		   int check_nmers,
                   int check_qvs,
		   int32 batchID,
  	           time_t currentTime,
		   int assembler,
		   int verbose){
  

  PHashValue_AS value;

  Transfer_FRG_to_IFG_AS(frg_mesg, ifg_mesg);

  switch(frg_mesg->action){

  case AS_ADD:
    {
      GateKeeperFragmentRecord gkf;
      int seqLength, quaLength;
      char *s;
      
      memset(&gkf, 0, sizeof(gkf));
      gkf.deleted = FALSE;
      gkf.birthBatch = batchID;
      gkf.deathBatch = 0;

      /* Check read types */
      switch((char)frg_mesg->type){
      case AS_READ:
      case AS_EXTR:
      case AS_TRNR:
	gkf.type = (char)frg_mesg->type;
	break;
      default:
	printGKPError(Msgfp, GKPError_Scalar);
	fprintf(Msgfp,"# Check_FragMessage: invalid type  %c \n",
		frg_mesg->type);
	return GATEKEEPER_FAILURE;
      }

      /* Check entryTime */
      if(frg_mesg->entry_time > currentTime){
	printGKPError(Msgfp, GKPError_Time);
	fprintf(Msgfp,"# Check_FragMessage: invalid entry time " F_TIME_T " > current time (" F_TIME_T ")\n",
		frg_mesg->entry_time, currentTime);
	return GATEKEEPER_FAILURE;
      }
      
      /* Make sure we haven't seen this frag record before... if so
	 it is a fatal error */
      if(HASH_FAILURE != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   frg_mesg->eaccession, 
                                                   AS_IID_FRG, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){

        printGKPError(Msgfp, GKPError_BadUniqueFRG);
        fprintf(Msgfp,"# Check_FragMessage:  A message with UID " F_U64 " exists with %d refs %s!!! Can't add it... \n",
                frg_mesg->eaccession, value.refCount, (value.deleted?"and was deleted":""));
        return(GATEKEEPER_FAILURE);
      }

      seqLength = 0;
      if(GATEKEEPER_FAILURE == checkSequence(frg_mesg->sequence, &s, &seqLength)){
	printGKPError(Msgfp, GKPError_FRGSequence);
	fprintf(Msgfp,"# Check_FragMessage: invalid char %c at position " F_SIZE_T " in sequence\n",
		*s, s - frg_mesg->sequence);
	return GATEKEEPER_FAILURE;
      }
      
      /* Compute quality length and check for legal characters */
      quaLength = 0;
      if(GATEKEEPER_FAILURE == checkQuality(frg_mesg->quality, &s, &quaLength)){
	printGKPError(Msgfp, GKPError_FRGQuality);
	fprintf(Msgfp,"# Check_FragMessage: invalid char %c at position " F_SIZE_T " in quality\n",
		*s, s - frg_mesg->quality);
	return GATEKEEPER_FAILURE;
      }

      //  Check the length before checking the quality, to handle
      //  fragments with too short of a clear range.
      //
      /* Check lengths and sequence intervals */
      if(GATEKEEPER_FAILURE == CheckLengthsIntervals(frg_mesg, ifg_mesg,
                                                     seqLength, quaLength,
                                                     assembler,
                                                     verbose)){
        /*
        fprintf(stderr,"# Check_FragMessage: ID: " F_U64 " lengths and intervals are incompatible\n",
                frg_mesg->eaccession);
        */
	return GATEKEEPER_FAILURE;
      }


      if( check_qvs ){
	double fractionError = checkOverallQuality(frg_mesg->quality, frg_mesg->clear_rng);
	if(fractionError > (1.01 * GATEKEEPER_MAX_ERROR_RATE)){
	  printGKPError(Msgfp, GKPError_FRGQualityGlobal);
	  fprintf(Msgfp,"# Global quality failed: ID: " F_U64 " eprob: %g in (" F_S32 "," F_S32 ")\n",
		  frg_mesg->eaccession,fractionError, frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end);
	  return GATEKEEPER_FAILURE;
	}
	fractionError = checkWindowQuality(frg_mesg, Msgfp);
	if(fractionError > 0){
	  return GATEKEEPER_FAILURE;
	}
      }

      value.type = AS_IID_FRG;
      InsertInPHashTable_AS(&(GkpStore.hashTable), UID_NAMESPACE_AS, frg_mesg->eaccession, &value, FALSE, TRUE);
      ifg_mesg->iaccession = value.IID;

      gkf.readUID  = ifg_mesg->eaccession;
      gkf.linkHead = 0;
      gkf.numLinks = 0;
      if(verbose)
	fprintf(stderr,"* Appending frag " F_U64 " " F_U32 " of type %c to store\n",
		gkf.readUID, value.IID, gkf.type);

      switch(ifg_mesg->type){
      case AS_READ:
      case AS_EXTR:
      case AS_TRNR:
	if(check_nmers){
	  char *seq5p = ifg_mesg->sequence + ifg_mesg->clear_rng.bgn; // 5p
	  char *seq3p = ifg_mesg->sequence + ifg_mesg->clear_rng.end - 8; // 3p
	  char *seqSanity = ifg_mesg->sequence + ifg_mesg->clear_rng.end - 50; // Sanity
          int clearRangeLength = ifg_mesg->clear_rng.end - ifg_mesg->clear_rng.bgn;

	  IncrementSequenceBucketArrayPrefix(LinkerDetector_READ, seq5p);
	  IncrementSequenceBucketArrayPrefix(Linker3pDetector_READ, seq3p); // 3p
	  IncrementSequenceBucketArrayPrefix(SanityDetector_READ, seqSanity);
	  IncrementBucketTotal(SequenceProbabilities, ifg_mesg->sequence + ifg_mesg->clear_rng.bgn, 
				       ifg_mesg->clear_rng.end - ifg_mesg->clear_rng.bgn);
	  IncrementSequenceLengthHistogram(Linker5pHistogram, seq5p, clearRangeLength, frg_mesg->eaccession);
	  IncrementSequenceLengthHistogram(Linker3pHistogram, seq3p, clearRangeLength, frg_mesg->eaccession);
	  IncrementSequenceLengthHistogram(LinkerSanityHistogram,seqSanity, clearRangeLength, frg_mesg->eaccession );
	}else{
	  // fprintf(stderr,"*** No nmer checking\n");
	}
        appendGateKeeperFragmentStore(GkpStore.frgStore, &gkf);
        break;
      default:
	  appendGateKeeperFragmentStore(GkpStore.frgStore, &gkf);
	  break;

      }
    }
    break;
  case AS_DELETE:
    {
      GateKeeperFragmentRecord gkf;

      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
						   UID_NAMESPACE_AS,
						   frg_mesg->eaccession, 
						   AS_IID_FRG, 
						   FALSE,
						   Msgfp,
						   &value)){

	printGKPError(Msgfp, GKPError_MissingFRG);
	fprintf(Msgfp,"# Check_FragMessage:  Fragment " F_U64 " DOES NOT exist!!!" 
		"# Can't delete it... bye\n",	    frg_mesg->eaccession);
	return(GATEKEEPER_FAILURE);
      }
      ifg_mesg->iaccession = value.IID;
      if(value.refCount > 1){
	printGKPError(Msgfp, GKPError_DeleteFRG);
	fprintf(Msgfp,"# Check_FragMessage: There are %d references outstanding to Fragment "
		F_U64 ".\n"
		"#  Can't delete it...\n",
		value.refCount, frg_mesg->eaccession);
	/*************** REPORT THE REFERENCES **********************/
	return(GATEKEEPER_FAILURE);
      }
      getGateKeeperFragmentStore(GkpStore.frgStore, value.IID, &gkf);
      if(gkf.numLinks > 0){
	printGKPError(Msgfp, GKPError_DeleteFRG);
	fprintf(Msgfp,"# Check_FragMessage:  Fragment has %d remaining links %d references\n", gkf.numLinks, value.refCount);
	return(GATEKEEPER_FAILURE);
      }      
      getGateKeeperFragmentStore(GkpStore.frgStore, value.IID, &gkf);
      if(verbose)
	fprintf(stderr,"* Deleting fragment...refcount = %d\n", value.refCount);
      if(HASH_SUCCESS == DeleteFromPHashTable_AS(GkpStore.hashTable, UID_NAMESPACE_AS, frg_mesg->eaccession)){
	deleteAndMarkGateKeeperFragmentStore(GkpStore.frgStore, value.IID, batchID);
      }else{
	assert(0);
      }

    }
    break;
    
#if 0
  case AS_UPDATE:
    {
      int seqLength, quaLength;

      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(GkpStore.hashTable, 
                                                   UID_NAMESPACE_AS,
                                                   frg_mesg->eaccession, 
                                                   AS_IID_FRG, 
                                                   FALSE,
                                                   Msgfp,
                                                   &value)){
        fprintf(Msgfp,"# Check_FragMessage:  Fragment " F_U64 " DOES NOT exist!!!" 
                "# Can't update it... bye\n",	    frg_mesg->eaccession);
        return(GATEKEEPER_FAILURE);
      }

      seqLength = strlen( frg_mesg->sequence );
      quaLength = strlen( frg_mesg->quality );
      if(GATEKEEPER_FAILURE == CheckLengthsIntervals(frg_mesg, ifg_mesg,
                                                     seqLength, quaLength,
                                                     verbose)){
        fprintf(Msgfp,"# Check_FragMessage: ID: " F_U64 " lengths and intervals are incompatible\n",
                frg_mesg->eaccession);
        return GATEKEEPER_FAILURE;
      }
	    
    }
    break;
#endif
  default:
    fprintf(Msgfp,"# Check_FragMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }
  return GATEKEEPER_SUCCESS;
}



/***************************************************************************/

int CheckLengthsIntervals(FragMesg *frg_mesg,
                          InternalFragMesg *ifg_mesg,  
                          int seqLength,
                          int quaLength,
                          int assembler,
                          int verbose)
{
  /* Check that sequence and quality lengths
     are the same on all fragments */
  if(verbose){
    fprintf(stderr,"* VERBOSE\n");
  }
  if(quaLength != seqLength){
    printGKPError(Msgfp, GKPError_FRGLength);
    fprintf(Msgfp,"# Check_FragMessage: sequence length (%d) != (%d) quality length\n",
            seqLength, quaLength);
    return GATEKEEPER_FAILURE;
  }

  {
    int maxlen;
    int minlen;

    if (assembler == AS_ASSEMBLER_OBT) {
      minlen = 0;
      maxlen = AS_FRAG_MAX_LEN;
    } else {
      minlen = AS_FRAG_MIN_LEN;
      maxlen = AS_FRAG_MAX_LEN;
    }

    /* Check that lengths are legit -- not too long */
    if(seqLength > maxlen){
      printGKPError(Msgfp, GKPError_FRGLength);
      fprintf(Msgfp,"# Check_FragMessage: sequence length (%d) > (%d) Max allowed sequence length\n",
	      seqLength, maxlen);
      return GATEKEEPER_FAILURE;
    }
    /* Check that lengths are legit -- not too short */
    if(seqLength < minlen){
      printGKPError(Msgfp, GKPError_FRGLength);
      fprintf(Msgfp,"# Check_FragMessage: sequence length (%d) < (%d) Min allowed sequence length\n",
	      seqLength, minlen);
      return GATEKEEPER_FAILURE;
    }
  

    /* Check clear range */
    if(0 >  frg_mesg->clear_rng.bgn ||
       frg_mesg->clear_rng.end > seqLength){
      printGKPError(Msgfp, GKPError_FRGClrRange);
      fprintf(Msgfp,"# Check_FragMessage: Invalid clear range (" F_S32 "," F_S32 ") valid range is (0,%d)\n",
	      frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end,
	      seqLength);
      return GATEKEEPER_FAILURE;
    }
  
    /* Check clear range */
    if((frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn) < minlen){
      printGKPError(Msgfp, GKPError_FRGClrRange);
      fprintf(Msgfp,"# Check_FragMessage: clear range length (" F_S32 ") < (%d) Min allowed length\n",
	      frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn, minlen);
      return GATEKEEPER_FAILURE;
    }
  
  }

  return GATEKEEPER_SUCCESS;
}



void InitQualityToFractionError(void){
  int i;
  for(i = 0; i <= 60; i++){
    double qualDiv10 = (double)i/10.0;
    double fError = pow(10.0, -qualDiv10);

    qualityToFractionError[i] = fError;
  }

}

				

/***********************************************************************************/
int checkSequence(char *input, char **errorChar, int *length){
  char *s;
	
  *length = 0;
  *errorChar = NULL;

  /* Compute sequence length and check for legal characters */
  for(s = input; *s != '\0'; s++){
    if(isspace(*s))
      continue;
    switch(tolower(*s)){
    case 'a':
    case 't':
    case 'c':
    case 'g':
    case 'n':
      (*length)++;
      break;
    default:
      *errorChar = s;
      return GATEKEEPER_FAILURE;
    }
  }
  return GATEKEEPER_SUCCESS;
}

/***********************************************************************************/
int checkQuality(char *input, char **errorChar, int *length){
  char *s;
	
  *length = 0;
  *errorChar = NULL;

  /* Compute quality length and check for legal characters */
  for(s = input; *s != '\0'; s++){
    if(isspace(*s))
      continue;
    if(*s >= '0' && *s <= 'l'){
      (*length)++;
    }else{
      *errorChar = s;
      return GATEKEEPER_FAILURE;
    }
  }
  return GATEKEEPER_SUCCESS;
}

/***********************************************************************************/

double checkOverallQuality(char *input, SeqInterval clearRange){
  int i;
  CDS_COORD_t length = clearRange.end - clearRange.bgn;
  double cumError = 0.0;
  double normError = 0.0;

  assert(length > 0);

  for(i = clearRange.bgn; i < clearRange.end; i++){
    cumError += qualityToFractionError[(int)(input[i] - '0')];
  }

  normError = cumError / (double)length;

  //  fprintf(stderr,"* checkOverallQuality on length " F_COORD " returns %g\n",
  //	  length, normError);

  return normError;

}

/* Checks every window of width GATEKEEPER_QV_WINDOW_WIDTH 
   We no longer do something special on the tails */
double checkWindowQuality(FragMesg *frg_mesg, FILE *Msgfp) {
  char *input=frg_mesg->quality;
  SeqInterval clearRange=frg_mesg->clear_rng;
  int i;
  CDS_COORD_t length = clearRange.end - clearRange.bgn;
  double cumError = 0.0;
  double normError = 0.0;
  SeqInterval window;
  // clearRange
  //   bgn                                                         end
  //   <------------>---------------------------------|------------|
  //   |------------|<----->--------------------------|------------|
  //   |------------|-<----->-------------------------|------------|
  //   |------------|--<----->------------------------|------------|
  //   |------------|------------ ...-----------------|------------|
  //   |------------|------------------------<----->--|------------|
  //   |------------|-------------------------<----->-|------------|
  //   |------------|--------------------------<----->|------------|
  //   |------------|---------------------------------<------------>
  //    TAIL_WIDTH     <     >
  //                      \ WINDOW_WIDTH

  assert(length > 0);

  cumError=0.0;
  // prime the pump
  window.bgn =  clearRange.bgn;
  window.end =  window.bgn + GATEKEEPER_QV_WINDOW_WIDTH;
  for(i = window.bgn; i< window.end; i++){
    cumError += qualityToFractionError[(int)(input[i] - '0')];
  }
  normError = cumError / (double)GATEKEEPER_QV_WINDOW_WIDTH;
  if (normError > (1.01 * GATEKEEPER_QV_WINDOW_THRESH))  { 
    printGKPError(Msgfp, GKPError_FRGQualityWindow);
    fprintf(Msgfp,"# Window quality failed: ID: " F_U64 " eprob: %g > %g in (" F_S32 "," F_S32 ") clr_rng: (" F_S32 "," F_S32 ")\n",
	    frg_mesg->eaccession,normError,
	    GATEKEEPER_QV_WINDOW_THRESH,
	    window.bgn,window.end,clearRange.bgn,clearRange.end);
    return normError; 
  }
  for(;window.end <clearRange.end;) {
    cumError-= qualityToFractionError[(int)(input[window.bgn] - '0')];
    cumError+= qualityToFractionError[(int)(input[window.end] - '0')];
    window.bgn++;
    window.end++;
    normError = cumError / (double)GATEKEEPER_QV_WINDOW_WIDTH;
    if (normError > (1.01 * GATEKEEPER_QV_WINDOW_THRESH))  { 
      printGKPError(Msgfp, GKPError_FRGQualityWindow);
      fprintf(Msgfp,"# Window quality failed: ID: " F_U64 " eprob: %g > %g in (" F_S32 "," F_S32 ") clr_rng: (" F_S32 "," F_S32 ")\n",
	    frg_mesg->eaccession,normError,
	    GATEKEEPER_QV_WINDOW_THRESH,
	    window.bgn,window.end,clearRange.bgn,clearRange.end);
      return normError; 
    }
  }
  /*
  fprintf(stderr,"* checkOverallQuality on length " F_COORD " returns %g\n",
  	  length, normError);
  */

  return 0; // zero return indicates passing the QV window test.
}
