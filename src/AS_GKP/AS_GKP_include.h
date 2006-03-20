
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
/* $Id: AS_GKP_include.h,v 1.5 2006-03-20 15:39:18 mhayton Exp $ */

/*************************************************
* Module:  AS_GKP_include.h
* Description:
*    Gatekeeper Header
*
* Reference:
*    GateKeeper.rtf
* 
*    Programmer:  S. Kravitz
*       Written:  Jan 1999
* 
*************************************************/

#ifndef AS_GKP_INCLUDE_H
#define AS_GKP_INCLUDE_H

#define GATEKEEPER_SUCCESS 0
#define GATEKEEPER_WARNING 1
#define GATEKEEPER_FAILURE 2

#include <stdio.h>
#include <errno.h>
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"
#include "AS_UTL_PHash.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_SequenceBucket.h"

#define GATEKEEPER_SCREENER_MIN_LENGTH 40
#define GATEKEEPER_MAX_ERROR_RATE 0.025 
#define GATEKEEPER_QV_WINDOW_WIDTH 50
#define GATEKEEPER_QV_WINDOW_THRESH 0.03

#define GATEKEEPER_MAX_WELL_NUMBER  384

#define AS_ASSEMBLER_GRANDE ((int)'A')
#define AS_ASSEMBLER_OVERLAY ((int)'O')


/* The persistent symbol table has 3 name spaces, UIDs,RPDIDs,LOCALEIDs */
#define UID_NAMESPACE_AS 'U'
#define RPTID_NAMESPACE_AS 'R'
#define LOCALEID_NAMESPACE_AS 'L'

extern MesgReader Reader;
extern MesgWriter Writer, ErrorWriter;
extern FILE *Infp, *Outfp, *Msgfp, *Ignfp, *Ibcfp, *Msgfp;
extern GateKeeperStore GkpStore;
extern SequenceBucketArrayT *LinkerDetector_READ;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence
extern SequenceBucketArrayT *Linker3pDetector_READ;   // Used to collect stats for 1-mers - 8-mers on 3-prime ends of sequence
extern SequenceBucketArrayT *LinkerDetector_EBAC;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence
extern SequenceBucketArrayT *LinkerDetector_LBAC;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence
extern SequenceBucketArrayT *SanityDetector_READ;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence
extern SequenceBucketArrayT *SanityDetector_EBAC;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence
extern SequenceBucketArrayT *SanityDetector_LBAC;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence
extern SequenceBucketT *SequenceProbabilities; // Used to collect stats on all sequence within clear ranges
extern SequenceLengthHistogramT *Linker5pHistogram ;
extern SequenceLengthHistogramT *Linker3pHistogram ;
extern SequenceLengthHistogramT *LinkerSanityHistogram;


extern char  File_Name_Prefix [FILENAME_MAX];
extern   char  Bac_File_Name [FILENAME_MAX];
extern char  Ignore_File_Name [FILENAME_MAX];
extern   char  Error_File_Name [FILENAME_MAX];
extern   char  Input_File_Name [FILENAME_MAX];
extern   char  Output_File_Name [FILENAME_MAX];

/***********************************************************************************
 * Function: Check_RepeatItemMesg
 * Description:
 *     Check RepeatItem message for correctness, and if correct, add it to the hashTable
 *
 * Checks:
 *    1) Repeat item has not been previously assigned
 *    2) Repeat length > 0
 * Inputs:
 *     rpt_mesg      *RepeatItemMesg
 *
 * I/O
 *     hashtable     PHashTable_AS
 *     msgFile       FILE * for diagnostic output.
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_FAILURE if failue.
 ***********************************************************************************/
int Check_RepeatItemMesg(RepeatItemMesg *rpt_msg,
			 InternalRepeatItemMesg *irp_msg,
			 CDS_CID_t batchID,
			 int verbose);


/***********************************************************************************
 * Function: Check_DistanceMesg
 * Description:
 *     Check Distance message for correctness, and if correct, add it to the GateKeeperStore,
 *     and hashtable.
 *     Maps UIDs to IIDs.
 * Checks:
 *    1) UID of message hasn't been previously seen
 *    2) Action is either add or delete
 *    3) If ADD, median > 0 and delta >= 0
 * Inputs:
 *     dst_mesg      *JoinMesg
 *
 * I/O
 *     idt_mesg      *InternalDistanceMesg -- filled in by this call, valid if return is SUCCESS
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_FAILURE if failue.
 ***********************************************************************************/
int Check_DistanceMesg(DistanceMesg *dst_mesg,  
		       InternalDistMesg *idt_mesg,  
		       CDS_CID_t batchID,
		       int verbose);


/***********************************************************************************
 * Function: Check_FragMesg
 * Description:
 *     Check Fragment message for correctness, and if correct, add it to the GateKeeperFragmentStore,
 *     and hashtable.
 *     Maps UIDs to IIDs.
 * Checks on add:
 *    1) UID of message hasn't been previously seen
 *    2) Action is either add or delete
 *    3) If ADD, median > 0 and delta >= 0
 *    4) length of sequence GATEKEEPER_MAX_SEQUENCE_LENGTH
 *    5) sequence is a string on [actgnACTGN]
 *    6) quality is a string on ASCII 0x30-0x6c
 *    7) length(sequence) == length(quality)
 *    8) clear range is bounded by  [0,sequence length]
 *    9) fragment entry time is < currentTime
 * Checks on delete:
 *    1) UID of referenced fragment have been previously defined
 *    2) The fragment has no outstanding references
 * 
 * Inputs:
 *     frg_mesg      *FragMesg
 *     check_qvs     flag to check or not check quality values
 *     currentTime   current time
 *
 * I/O
 *     ifg_mesg      *InternalFragMesg -- filled in by this call, valid if return is SUCCESS
 *     hashtable     PHashTable_AS
 *     msgFile       FILE * for diagnostic output.
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_FAILURE if failue.
 ***********************************************************************************/
int Check_FragMesg(FragMesg *frg_mesg,  
                   InternalFragMesg *ifg_mesg,
		   int check_nmers,
                   int check_qvs,
		   CDS_CID_t batchID,
                   time_t currentTime,
		   int assembler,
		   int verbose);


/***********************************************************************************
 * Function: Check_ScreenItemMesg
 * Description:
 *     Check ScreenItem message for correctness.
 *     Maps UIDs to IIDs.
 * Checks:
 *    1) UID of message hasn't been previously seen
 *    2) If variation in range [0,1]  
 *    3) Type is either UR or Contaminent
 *    4) Repeat ID has been previously defined
 *    5) Sequence is a string on [atcgnACTGN]
 *    6) min_length is >= GATEKEEPER_SCREENER_MIN_LENGTH
 * Inputs:
 *     scn_mesg      *ScreenItemMesg
 *
 * I/O
 *     isn_mesg      *InternalScreenItemMesg -- filled in by this call, valid if return is SUCCESS
 *     hashtable     PHashTable_AS
 *     msgFile       FILE * for diagnostic output.
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_FAILURE if failue.
 ***********************************************************************************/
int Check_ScreenItemMesg(ScreenItemMesg *scn_mesg,  
		         InternalScreenItemMesg *isn_mesg,  
			 CDS_CID_t batchID,
			 int verbose);

/***********************************************************************************
 * Function: Check_LinkMesg
 * Description:
 *     Check Link message for correctness. and if correct, add it to the GateKeeperFragmentStore,
 *     and hashtable.  Add symbol table references to frags and dist.
 *     Maps UIDs to IIDs.
 * Checks on add:
 *    1) UID of referenced fragments have been previously defined
 *    2) referenced fragments are unique
 *    3) referenced fragments are are of appropriate type for type of link
 *    4) Entry time is in the past
 *    5) UID of distance record has been previously defined
 *    6) Mates:  Have no previous mate, no reread relation between the pari
 *    7) Guides: No previous guide link between the pair, no reread relation
 *    8) ReReads: No guide/mate/reread between the pair.
 * Checks on delete:
 *    1) UID of referenced fragments have been previously defined
 *    2) The link in question exists.
 * Inputs:
 *     lnk_mesg      *LinkMesg
 *     currentTime   current time
 *
 * I/O
 *     ilk_mesg      *InternalLinkMesg -- filled in by this call, valid if return is SUCCESS
 *     hashtable     PHashTable_AS
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
 ***********************************************************************************/
int Check_LinkMesg(LinkMesg *lnk_mesg,
 	           InternalLinkMesg *ilk_mesg,  
		   CDS_CID_t batchID,
  	           time_t currentTime,
		   int verbose,
                   int matchBAC_DstsInLkgMsgs);





/***********************************************************************************
 * Function: Check_BacMesg
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
 ***********************************************************************************/
int Check_BacMesg(BacMesg *bac_mesg,
		  InternalBacMesg *ibc_mesg,
		  VA_TYPE(int32) *addedBacs,
		  CDS_CID_t batchID,
		  time_t currentTime,
		  int assembler,
		  int strict,
		  int verbose);

/******************************************************************************
 * Function: Check_WellMesg
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
 *****************************************************************************/
/*
int Check_WellMesg(WellMesg *wel_mesg,
		  InternalWellMesg *iwe_mesg,
		  CDS_CID_t batchID,
		  time_t currentTime,
		  int assembler,
		  int strict,
		  int verbose);
*/

/******************************************************************************
 * Function: Check_SeqPlateMesg
 * Inputs:
 * I/O
 * Return Value:
 *     GATEKEEPER_SUCCESS if success. 
 *     GATEKEEPER_WARNING if success but warning
 *     GATEKEEPER_FAILURE if failure
******************************************************************************/
/*
int Check_SeqPlateMesg(SeqPlateMesg *sqp_mesg,
		       InternalSeqPlateMesg *isq_mesg,
		       CDS_CID_t batchID,
		       time_t currentTime,
		       int assembler,
		       int strict,
		       int verbose);
*/
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
                    int verbose);

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
                        int verbose);

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
                       int verbose);

/***********************************************************************************
 * Function: Check_BatchMesg
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
 ***********************************************************************************/
int Check_BatchMesg(BatchMesg *bat_mesg,
		    InternalBatchMesg *iba_mesg,
		    time_t currentTime,
		    int verbose);



/*******************************************************************************/
/*                       Utility Functions                                     */
/*******************************************************************************/


/***********************************************************************************
 * Function: findlinksInCompatibleWith
 * Description:
 *     Searches the GateKeeperLinkStore gkplStore from the record with index linkHead
 *     for links that conflict with link.  If found, returns the index of the offending
 *     link.
 ***********************************************************************************/
int findLinksInCompatibleWith(GateKeeperLinkStore gkplStore, 
			      CDS_IID_t followFrag, 
			      CDS_IID_t linkHead, 
			      GateKeeperLinkRecord *link,
			      FILE *msgFile,
			      int verbose);




/* Check the sequence and quality data and report errors */
int checkSequence(char *input, char **errorChar, int *length);
int checkQuality(char *input, char **errorChar, int *length);
double checkOverallQuality(char *input, SeqInterval clearRange);
double checkWindowQuality(FragMesg *frg_mesg, FILE *msgFile);

int CheckLengthsIntervalsLocales(FragMesg *frg_mesg,
                                 InternalFragMesg *ifg_mesg,  
                                 int seqLength,
                                 int quaLength,
				 int assembler,
				 int verbose);

FILE *  File_Open(const char * Filename, const char * Mode, int exitOnFailure);


#define OUTPUT_FILENAME_EXTENSION ".inp"
#define INITIAL_HASH_SIZE (256 * 2048)



/***********************************************************************************/
/* ReadFile:  Main driver routine for the gatekeeper                               */
/***********************************************************************************/
int ReadFile(int check_qvs,
	     int check_nmers,
	     VA_TYPE(int32) *addedBacs,
	     CDS_CID_t batchID,
	     int32 assembler,
	     int32 strict,
	     int argc, char **argv,
	     int32 verbose,
             int matchBAC_DstsInLkgMsgs);

/***********************************************************************************/
/* CheckOverlapInputComplete                                                       */
/***********************************************************************************/
int CheckOverlayInputComplete(VA_TYPE(int32) *addedBacs, int verbose);

// Initialization Function
void InitQualityToFractionError(void);


/**************************************************************************************/
/* Gatekeeper Errors                                                                  */

typedef enum GKPErrorType_tag{
  GKPError_Invalid = -1,
  GKPError_FirstMessageBAT = 0,
  GKPError_BadUniqueBAT=100,
  GKPError_BadUniqueFRG,
  GKPError_BadUniqueDST,
  GKPError_BadUniqueBAC,
  GKPError_BadUniqueBTG,
  GKPError_BadUniqueSEQ,
  GKPError_BadUniqueRPT,
  GKPError_BadUniqueSCN,
  GKPError_BadUniqueWEL,
  GKPError_BadUniquePLA,

  GKPError_MissingFRG=200,
  GKPError_MissingDST,
  GKPError_MissingBAC,
  GKPError_MissingBTG,
  GKPError_MissingSEQ,
  GKPError_MissingRPT,
  GKPError_MissingPLA,

  GKPError_DeleteFRG=300,
  GKPError_DeleteDST,
  GKPError_DeleteBAC,
  GKPError_DeleteLNK,

  GKPError_Time=400,
  GKPError_Action,
  GKPError_Scalar,
  GKPError_IncompleteOAInput,

  GKPError_FRGSequence=500,
  GKPError_FRGQuality,
  GKPError_FRGLength,
  GKPError_FRGClrRange,
  GKPError_FRGLocalPos,
  GKPError_FRGAccession,
  GKPError_FRGBacSeqMismatch,
  GKPError_FRGBacFragAlreadyDefined,
  GKPError_FRGBactigFragAlreadyDefined,
  GKPError_FRGBactigInWrongBac,
  GKPError_FRGWrongTypeForOverlay,
  GKPError_FRGWrongTypeForBAC,
  GKPError_FRGQualityWindow,
  GKPError_FRGQualityGlobal,
  GKPError_FRGQualityTail,

  GKPError_LNKLocaleMismatch=600,
  GKPError_LNKLocaleDistanceMismatch,
  GKPError_LNKFragtypeMismatch,
  GKPError_LNKOneLink,

  GKPError_BACNumBactigs = 700,
  GKPError_BACNumBactigsShouldBeZero,
  GKPError_BACRedefinition,
  GKPError_BACWrongTypeForOverlay,
  GKPError_BACMultiplyDefinedBactig,

  GKPError_DSTValues = 800,
  
  GKPError_SCNVariation=900,
  GKPError_SCNminLength,

  GKPError_RPTLength = 1000,

  GKPError_WellNumberOutOfRange = 1100,
  GKPError_MAX
  
}GKPErrorType;

#define MAX_GKPERROR ((int)(GKPError_MAX )-1)

void printGKPError(FILE *fout, GKPErrorType type);
void printAllGKPErrors(FILE *fout);

int32 CheckNmerProbabilities(FILE *, double threshhold);

#endif
