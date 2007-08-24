
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

static char CM_ID[] = "$Id: AS_GKP_checkFrag.c,v 1.29 2007-08-24 15:29:48 brianwalenz Exp $";

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
double checkWindowQuality(FragMesg *frg_mesg) {
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
    fprintf(errorFP, "# FRG Error: Fragment "F_UID" failed window quality; eprob: %g > %g in range (" F_S32 "," F_S32 ")\n",
	    frg_mesg->eaccession,normError,
	    GATEKEEPER_QV_WINDOW_THRESH,
	    window.bgn,window.end);
    return normError; 
  }
  for(;window.end <clearRange.end;) {
    cumError-= qualityToFractionError[(int)(input[window.bgn] - '0')];
    cumError+= qualityToFractionError[(int)(input[window.end] - '0')];
    window.bgn++;
    window.end++;
    normError = cumError / (double)GATEKEEPER_QV_WINDOW_WIDTH;
    if (normError > (1.01 * GATEKEEPER_QV_WINDOW_THRESH))  { 
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" failed window quality; eprob: %g > %g in range (" F_S32 "," F_S32 ")\n",
              frg_mesg->eaccession,
              normError,
              GATEKEEPER_QV_WINDOW_THRESH,
              window.bgn,window.end);
      return normError; 
    }
  }

  return 0; // zero return indicates passing the QV window test.
}







int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   check_qvs,
               int                   assembler) {
  static  uint32    libOrientationMax = 0;
  static  uint64   *libOrientation    = NULL;

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
    GateKeeperFragmentRecord gkf = {0};
    int                      seqLength = 0;
    int                      quaLength = 0;
    char                    *s = NULL;

    clearGateKeeperFragmentRecord(&gkf);

    //  Make sure we haven't seen this frag record before... if so
    //  it is a fatal error
    //
    if (frg_mesg->eaccession == 0) {
      fprintf(errorFP, "# FRG Error: Fragment has zero or no UID; can't add it.\n");
      return(GATEKEEPER_FAILURE);
    }
    if (getGatekeeperUIDtoIID(gkpStore, frg_mesg->eaccession, NULL)) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" exists, can't add it again.\n",
              frg_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }


    //  Check sequence and quality
    //
    if (GATEKEEPER_FAILURE == checkSequence(frg_mesg->sequence, &s, &seqLength)) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" invalid char %c at position "F_SIZE_T" in sequence.\n",
              frg_mesg->eaccession, *s, s - frg_mesg->sequence);
      return GATEKEEPER_FAILURE;
    }
      
    if (GATEKEEPER_FAILURE == checkQuality(frg_mesg->quality, &s, &quaLength)) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" invalid char %c at position " F_SIZE_T " in quality.\n",
              frg_mesg->eaccession, *s, s - frg_mesg->quality);
      return GATEKEEPER_FAILURE;
    }


    //  Check the length before checking the quality, to handle
    //  fragments with too short of a clear range.

    //  Check sequence and quality are same length.
    if (quaLength != seqLength) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" sequence length %d != %d quality length\n",
              frg_mesg->eaccession, seqLength, quaLength);
      return GATEKEEPER_FAILURE;
    }

    //  Check that lengths are legit -- not too long */
    if (seqLength > AS_FRAG_MAX_LEN){
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" sequence length %d > %d max allowed sequence length\n",
              frg_mesg->eaccession, seqLength, AS_FRAG_MAX_LEN);
      return GATEKEEPER_FAILURE;
    }

    //  Check that lengths are legit -- not too short */
    if (seqLength < AS_FRAG_MIN_LEN){
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" sequence length %d < %d min allowed sequence length\n",
              frg_mesg->eaccession, seqLength, AS_FRAG_MIN_LEN);
      return GATEKEEPER_FAILURE;
    }

    //  Check clear range bounds
    if ((frg_mesg->clear_rng.bgn < 0) || (seqLength < frg_mesg->clear_rng.end)) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" invalid clear range (%d,%d) valid range is (0,%d)\n",
              frg_mesg->eaccession, frg_mesg->clear_rng.bgn, frg_mesg->clear_rng.end, seqLength);
      return GATEKEEPER_FAILURE;
    }
  
    //  Check clear range length
    if ((assembler == AS_ASSEMBLER_GRANDE) &&
        ((frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn) < AS_FRAG_MIN_LEN)) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" clear range length %d < %d min allowed length\n",
              frg_mesg->eaccession, frg_mesg->clear_rng.end - frg_mesg->clear_rng.bgn, AS_FRAG_MIN_LEN);
      return GATEKEEPER_FAILURE;
    }



    if (check_qvs) {
      double fractionError = checkOverallQuality(frg_mesg->quality, frg_mesg->clear_rng);

      if(fractionError > (1.01 * GATEKEEPER_MAX_ERROR_RATE)){
        fprintf(errorFP, "# FRG Error: Fragment "F_UID" eprob: %g in (" F_S32 "," F_S32 ")\n",
                frg_mesg->eaccession,
                fractionError,
                frg_mesg->clear_rng.bgn,
                frg_mesg->clear_rng.end);
        return GATEKEEPER_FAILURE;
      }

      fractionError = checkWindowQuality(frg_mesg);
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
      gkf.libraryIID = getGatekeeperUIDtoIID(gkpStore, frg_mesg->library_uid, NULL);

      if (gkf.libraryIID == 0) {
        fprintf(errorFP, "# FRG Error: Fragment "F_UID" references unknown library "F_UID"\n",
                frg_mesg->eaccession, frg_mesg->library_uid);
        return(GATEKEEPER_FAILURE);
      }

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
    gkf.readUID = frg_mesg->eaccession;
    gkf.readIID = getLastElemStore(gkpStore->frg) + 1;

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

    setGatekeeperUIDtoIID(gkpStore, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkpStore->frg, &gkf);

    appendVLRecordStore(gkpStore->seq, frg_mesg->sequence, gkf.seqLen);

    encodeSequenceQuality(encodedsequence, frg_mesg->sequence, frg_mesg->quality);
    appendVLRecordStore(gkpStore->qlt, encodedsequence,    gkf.seqLen);

    appendVLRecordStore(gkpStore->hps, NULL,               0);
    appendVLRecordStore(gkpStore->src, frg_mesg->source,   gkf.srcLen);

  } else if (frg_mesg->action == AS_DELETE) {

    GateKeeperFragmentRecord gkf;
    CDS_IID_t                iid = getGatekeeperUIDtoIID(gkpStore, frg_mesg->eaccession, NULL);

    if (iid == 0) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" does not exist, can't delete it.\n",
              frg_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }

    getGateKeeperFragment(gkpStore, iid, &gkf);

    if (gkf.mateIID > 0) {
      fprintf(errorFP, "# FRG Error: Fragment "F_UID" has mate pair relationship, can't delete it.\n",
              frg_mesg->eaccession);
      return(GATEKEEPER_FAILURE);
    }      

    if (HASH_SUCCESS == delGatekeeperUIDtoIID(gkpStore, frg_mesg->eaccession)) {
      GateKeeperFragmentRecord dr;
      getIndexStore(gkpStore->frg, iid, &dr);
      dr.deleted = TRUE;
      setIndexStore(gkpStore->frg, iid, &dr);
    }else{
      assert(0);
    }

  } else {
    fprintf(errorFP, "# FRG Error: invalid action %c\n", frg_mesg->action);
    return GATEKEEPER_FAILURE;
  }
  return GATEKEEPER_SUCCESS;
}
