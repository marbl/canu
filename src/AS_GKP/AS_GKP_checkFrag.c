
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

static char CM_ID[] = "$Id: AS_GKP_checkFrag.c,v 1.22 2007-05-02 09:30:15 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"


static int    checkfraginitialized = 0;
static double qualityToFractionError[61] = {0.0};
static int    isspacearray[256] = {0};
static int    isvalidACGTN[256] = {0};
static char   encodedsequence[AS_FRAG_MAX_LEN+1] = {0};

int
CheckLengthsIntervals(FragMesg *frg_mesg,
                      int seqLength,
                      int quaLength,
                      int assembler,
                      int verbose) {
  int minlen = AS_FRAG_MIN_LEN;
  int maxlen = AS_FRAG_MAX_LEN;

  if (assembler == AS_ASSEMBLER_OBT)
    minlen = 0;

  //  Check sequence and quality are same length.
  if (quaLength != seqLength) {
    printGKPError(stderr, GKPError_FRGLength);
    fprintf(stderr,"# Check_FragMessage: sequence length (%d) != (%d) quality length\n",
            seqLength, quaLength);
    return GATEKEEPER_FAILURE;
  }

  //  Check that lengths are legit -- not too long */
  if (seqLength > maxlen){
    printGKPError(stderr, GKPError_FRGLength);
    fprintf(stderr,"# Check_FragMessage: sequence length (%d) > (%d) Max allowed sequence length\n",
            seqLength, maxlen);
    return GATEKEEPER_FAILURE;
  }

  //  Check that lengths are legit -- not too short */
  if (seqLength < minlen){
    printGKPError(stderr, GKPError_FRGLength);
    fprintf(stderr,"# Check_FragMessage: sequence length (%d) < (%d) Min allowed sequence length\n",
            seqLength, minlen);
    return GATEKEEPER_FAILURE;
  }

  //  Check clear range bounds
  if ((0 > frg_mesg->clear_rng.bgn) || (frg_mesg->clear_rng.end > seqLength)) {
    printGKPError(stderr, GKPError_FRGClrRange);
    fprintf(stderr,"# Check_FragMessage: Invalid clear range (" F_S32 "," F_S32 ") valid range is (0,%d)\n",
            frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end,
            seqLength);
    return GATEKEEPER_FAILURE;
  }
  
  //  Check clear range length
  if ((frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn) < minlen) {
    printGKPError(stderr, GKPError_FRGClrRange);
    fprintf(stderr,"# Check_FragMessage: clear range length (" F_S32 ") < (%d) Min allowed length\n",
            frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn, minlen);
    return GATEKEEPER_FAILURE;
  }

  return GATEKEEPER_SUCCESS;
}



int
checkSequence(char *input, char **errorChar, int *length) {
  char *s;

  *length = 0;
  *errorChar = NULL;

  /* Compute sequence length and check for legal characters */
  for(s = input; *s != '\0'; s++){
    if(isspacearray[*s])
      continue;

    if (isvalidACGTN[*s]) {
      *s = isvalidACGTN[*s];
      (*length)++;
    } else {
      *errorChar = s;
      return GATEKEEPER_FAILURE;
    }
  }
  return GATEKEEPER_SUCCESS;
}



int
checkQuality(char *input, char **errorChar, int *length) {
  char *s;

  *length = 0;
  *errorChar = NULL;

  /* Compute quality length and check for legal characters */
  for(s = input; *s != '\0'; s++){
    if(isspacearray[*s])
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


double
checkOverallQuality(char *input, SeqInterval clearRange) {
  int i;
  CDS_COORD_t length = clearRange.end - clearRange.bgn;
  double cumError = 0.0;
  double normError = 0.0;

  assert(length > 0);

  for(i = clearRange.bgn; i < clearRange.end; i++){
    cumError += qualityToFractionError[(int)(input[i] - '0')];
  }

  normError = cumError / (double)length;

  return normError;

}



// Checks every window of width GATEKEEPER_QV_WINDOW_WIDTH 
// We no longer do something special on the tails
//
double checkWindowQuality(FragMesg *frg_mesg, FILE *err) {
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
    printGKPError(err, GKPError_FRGQualityWindow);
    fprintf(err,"# Window quality failed: ID: " F_U64 " eprob: %g > %g in (" F_S32 "," F_S32 ") clr_rng: (" F_S32 "," F_S32 ")\n",
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
      printGKPError(err, GKPError_FRGQualityWindow);
      fprintf(err,"# Window quality failed: ID: " F_U64 " eprob: %g > %g in (" F_S32 "," F_S32 ") clr_rng: (" F_S32 "," F_S32 ")\n",
              frg_mesg->eaccession,normError,
              GATEKEEPER_QV_WINDOW_THRESH,
              window.bgn,window.end,clearRange.bgn,clearRange.end);
      return normError; 
    }
  }

  return 0; // zero return indicates passing the QV window test.
}







int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   check_qvs,
               int                   assembler,
               int                   verbose) {
  static  uint32    libOrientationMax = 0;
  static  uint64   *libOrientation    = NULL;

  PHashValue_AS value;


  if (frg_mesg->action == AS_IGNORE)
    return GATEKEEPER_SUCCESS;


  if (checkfraginitialized == 0) {
    int i;

    for (i=0; i<256; i++)
      isspacearray[i] = isspace(i);

    isvalidACGTN['a'] = 'A';
    isvalidACGTN['c'] = 'C';
    isvalidACGTN['g'] = 'G';
    isvalidACGTN['t'] = 'T';
    isvalidACGTN['n'] = 'N';
    isvalidACGTN['A'] = 'A';
    isvalidACGTN['C'] = 'C';
    isvalidACGTN['G'] = 'G';
    isvalidACGTN['T'] = 'T';
    isvalidACGTN['N'] = 'N';

    for (i=0; i<61; i++)
      qualityToFractionError[i] = pow(10.0, -(double)i/10.0);

    checkfraginitialized = 1;
  }

  if (frg_mesg->action == AS_ADD) {
    GateKeeperFragmentRecord gkf;
    int                      seqLength = 0;
    int                      quaLength = 0;
    char                    *s = NULL;

    clearGateKeeperFragmentRecord(&gkf);

    //  Make sure we haven't seen this frag record before... if so
    //  it is a fatal error
    //
    if (HASH_FAILURE != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                  UID_NAMESPACE_AS,
                                                  frg_mesg->eaccession, 
                                                  AS_IID_FRG, 
                                                  FALSE,
                                                  stderr,
                                                  &value)) {
      printGKPError(stderr, GKPError_BadUniqueFRG);
      fprintf(stderr, "# Check_FragMessage:  A message with UID " F_U64 " exists with %d refs %s!!! Can't add it... \n",
              frg_mesg->eaccession, value.refCount, (value.deleted?"and was deleted":""));
      return(GATEKEEPER_FAILURE);
    }

    //  Check sequence and quality
    //
    if (GATEKEEPER_FAILURE == checkSequence(frg_mesg->sequence, &s, &seqLength)) {
      printGKPError(stderr, GKPError_FRGSequence);
      fprintf(stderr, "# Check_FragMessage: invalid char %c at position " F_SIZE_T " in sequence\n",
              *s, s - frg_mesg->sequence);
      return GATEKEEPER_FAILURE;
    }
      
    if (GATEKEEPER_FAILURE == checkQuality(frg_mesg->quality, &s, &quaLength)) {
      printGKPError(stderr, GKPError_FRGQuality);
      fprintf(stderr, "# Check_FragMessage: invalid char %c at position " F_SIZE_T " in quality\n",
              *s, s - frg_mesg->quality);
      return GATEKEEPER_FAILURE;
    }

    //  Check the length before checking the quality, to handle
    //  fragments with too short of a clear range.
    //
    if (GATEKEEPER_FAILURE == CheckLengthsIntervals(frg_mesg,
                                                    seqLength, quaLength,
                                                    assembler,
                                                    verbose)) {
      fprintf(stderr, "# Check_FragMessage: ID " F_U64 " lengths and intervals are incompatible\n",
              frg_mesg->eaccession);
      return GATEKEEPER_FAILURE;
    }


    if (check_qvs) {
      double fractionError = checkOverallQuality(frg_mesg->quality, frg_mesg->clear_rng);

      if(fractionError > (1.01 * GATEKEEPER_MAX_ERROR_RATE)){
        printGKPError(stderr, GKPError_FRGQualityGlobal);
        fprintf(stderr,"# Global quality failed: ID: " F_U64 " eprob: %g in (" F_S32 "," F_S32 ")\n",
                frg_mesg->eaccession,fractionError, frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end);
        return GATEKEEPER_FAILURE;
      }

      fractionError = checkWindowQuality(frg_mesg, stderr);
      if (fractionError > 0)
        return GATEKEEPER_FAILURE;
    }

    //  Version 2 comes with library information.  Get it.
    //
    //  Version 1 sets orient and library when mates are added.
    //  Version 1 libraries NEVER have a valid orientation.

    gkf.libraryIID  = 0;
    gkf.orientation = AS_READ_ORIENT_UNKNOWN;

    if (frg_mesg->library_uid > 0) {
      value.type = AS_IID_DST;
      value.IID  = 0;

      if (HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private,
                                                    UID_NAMESPACE_AS,
                                                    frg_mesg->library_uid,
                                                    AS_IID_DST,
                                                    FALSE,
                                                    stderr,
                                                    &value)) {
        printGKPError(stderr, GKPError_MissingLIB);
        fprintf(stderr,"# Check_FragMEssage: ID "F_UID" references unknown library "F_UID"\n",
                frg_mesg->eaccession, frg_mesg->library_uid);
        return(GATEKEEPER_FAILURE);
      }

      gkf.libraryIID  = value.IID;

      //  Get the library orientation; cache values across invocations.
      //
      if (libOrientation == NULL) {
        libOrientationMax = 256;
        libOrientation    = (uint64 *)safe_calloc(libOrientationMax, sizeof(uint64));
      }
      if (libOrientationMax <= gkf.libraryIID) {
        libOrientation    = safe_realloc(libOrientation, 2 * gkf.libraryIID);
        while (libOrientationMax < 2 * gkf.libraryIID)
          libOrientation[libOrientationMax++] = 0;
      }
      if (libOrientation[gkf.libraryIID] == 0) {
        GateKeeperLibraryRecord  gkpl;
        getIndexStore(gkpStore->lib, gkf.libraryIID, &gkpl); 
        libOrientation[gkf.libraryIID] = gkpl.orientation;
      }
      gkf.orientation = libOrientation[gkf.libraryIID];
    }


    {
      int which;

      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        gkf.clearBeg[which] = frg_mesg->clear_rng.bgn;
        gkf.clearEnd[which] = frg_mesg->clear_rng.end;
      }

      if (frg_mesg->clear_qlt.bgn <= frg_mesg->clear_qlt.end) {
        gkf.hasQualityClear = 1;
        gkf.clearBeg[AS_READ_CLEAR_QLT] = frg_mesg->clear_qlt.bgn;
        gkf.clearEnd[AS_READ_CLEAR_QLT] = frg_mesg->clear_qlt.end;
      } else {
        gkf.clearBeg[AS_READ_CLEAR_QLT] = 1;
        gkf.clearEnd[AS_READ_CLEAR_QLT] = 0;
      }

      if (frg_mesg->clear_vec.bgn <= frg_mesg->clear_vec.end) {
        gkf.hasVectorClear = 1;
        gkf.clearBeg[AS_READ_CLEAR_VEC] = frg_mesg->clear_vec.bgn;
        gkf.clearEnd[AS_READ_CLEAR_VEC] = frg_mesg->clear_vec.end;
      } else {
        gkf.clearBeg[AS_READ_CLEAR_VEC] = 1;
        gkf.clearEnd[AS_READ_CLEAR_VEC] = 0;
      }
    }

    //  Now add the fragment to the store
    //
    value.type = AS_IID_FRG;
    value.IID  = 0;

    InsertInPHashTable_AS(&gkpStore->phs_private,
                          UID_NAMESPACE_AS,
                          frg_mesg->eaccession,
                          &value,
                          FALSE,
                          TRUE);

    gkf.readUID = frg_mesg->eaccession;
    gkf.readIID = value.IID;

    gkf.seqLen = strlen(frg_mesg->sequence);
    gkf.hpsLen = 0;
    gkf.srcLen = strlen(frg_mesg->source);

    {
      StoreStat   stats;

      statsStore(gkpStore->seq, &stats);
      gkf.seqOffset = stats.lastElem;

      statsStore(gkpStore->qlt, &stats);
      gkf.qltOffset = stats.lastElem;

      statsStore(gkpStore->hps, &stats);
      gkf.hpsOffset = stats.lastElem;

      statsStore(gkpStore->src, &stats);
      gkf.srcOffset = stats.lastElem;
    }

    appendIndexStore(gkpStore->frg, &gkf);

    appendVLRecordStore(gkpStore->seq, frg_mesg->sequence, gkf.seqLen);

    encodeSequenceQuality(encodedsequence, frg_mesg->sequence, frg_mesg->quality);
    appendVLRecordStore(gkpStore->qlt, encodedsequence,    gkf.seqLen);

    appendVLRecordStore(gkpStore->hps, NULL,               0);
    appendVLRecordStore(gkpStore->src, frg_mesg->source,   gkf.srcLen);

  } else if (frg_mesg->action == AS_DELETE) {

    GateKeeperFragmentRecord gkf;

    if(HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore->phs_private, 
                                                 UID_NAMESPACE_AS,
                                                 frg_mesg->eaccession, 
                                                 AS_IID_FRG, 
                                                 FALSE,
                                                 stderr,
                                                 &value)){
      printGKPError(stderr, GKPError_MissingFRG);
      fprintf(stderr,"# Check_FragMessage:  Fragment " F_U64 " DOES NOT exist!!!" 
              "# Can't delete it... bye\n",	    frg_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    if (value.refCount > 1){
      printGKPError(stderr, GKPError_DeleteFRG);
      fprintf(stderr,"# Check_FragMessage: There are %d references outstanding to Fragment "F_U64".\n"
              "#  Can't delete it...\n",
              value.refCount, frg_mesg->eaccession);

      /*************** REPORT THE REFERENCES **********************/

      return(GATEKEEPER_FAILURE);
    }

    getGateKeeperFragment(gkpStore, value.IID, &gkf);

    if (gkf.mateIID > 0) {
      printGKPError(stderr, GKPError_DeleteFRG);
      fprintf(stderr,"# Check_FragMessage:  Fragment has mate pair relationship; %d references\n", value.refCount);
      return(GATEKEEPER_FAILURE);
    }      

    getGateKeeperFragment(gkpStore, value.IID, &gkf);

    if(verbose)
      fprintf(stderr,"* Deleting fragment...refcount = %d\n", value.refCount);

    if(HASH_SUCCESS == DeleteFromPHashTable_AS(gkpStore->phs_private, UID_NAMESPACE_AS, frg_mesg->eaccession)){
      //deleteGateKeeperFragment(gkpStore, value.IID);
      GateKeeperFragmentRecord dr;
      getIndexStore(gkpStore->frg, value.IID, &dr);
      dr.deleted = TRUE;
      setIndexStore(gkpStore->frg, value.IID, &dr);
    }else{
      assert(0);
    }

  } else {
    fprintf(stderr,"# Check_FragMessage: invalid action\n");
    return GATEKEEPER_FAILURE;
  }
  return GATEKEEPER_SUCCESS;
}



