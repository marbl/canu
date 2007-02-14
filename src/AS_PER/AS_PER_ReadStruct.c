
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

static char CM_ID[] = "$Id: AS_PER_ReadStruct.c,v 1.11 2007-02-14 07:20:13 brianwalenz Exp $";

/*************************************************************************
 Module:  AS_PER_ReadStruct
 Description:
     This module defines the interface and implementation of the 
 opaque datatype used by the Fragment Store.

 Assumptions:
      
 Document:
      FragStore.rtf

*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"


//  Enable this to get dumpFragStore, and anything that calls
//  dump_ReadStruct(), to print quality in the non-internal two-digits
//  format.  E.g., 03 01 03 10, etc.
//#define DUMP_QUALITY_AS_NUMBERS


ReadStruct *new_ReadStruct(void){
  ReadStruct *newFR = (ReadStruct *)safe_calloc( 1, sizeof(ReadStruct));
  clear_ReadStruct(newFR);
  return(newFR);
}

void        delete_ReadStruct(ReadStruct *rs){
  safe_free(rs);
}


void        clear_ReadStruct(ReadStruct *rs){
  clearGateKeeperFragmentRecord(&rs->gkfr);
  setClearRegion_ReadStruct(rs, 0, 0, READSTRUCT_ORIGINAL); 
}




int
dump_ReadStruct(ReadStruct *rs,
                FILE *fout,
                int clearRangeOnly){

  //  XXX  This should just call a gatekeeper function to dump the gatekeeper fragment.

#if 0

  fprintf(fout,"Dumping ReadStruct at 0x%p\n", fr);

  fprintf(fout, "\tDeleted: %d ReadType:%c hasQLT:%d spare1:%d\n",
          rs->gkfr.deleted, rs->gkfr.readType, rs->gkfr.hasQLT, rs->gkfr.spare1);

  fprintf(fout, 
          "\taccID:" F_UID "  readIdx:" F_IID "  ClearRanges(start,stop,modified):\n", 
          rs->gkfr.accID, rs->gkfr.readIndex);
  fprintf(fout, "\tOrig(" F_VLS "," F_VLS ") Ovl(" F_VLS "," F_VLS ",%c) Cns(" F_VLS "," F_VLS ",%c) Cgw(" F_VLS "," F_VLS ",%c)\n",
          rs->gkfr.clrSta, rs->gkfr.clrEnd,
          rs->gkfr.ovlSta, rs->gkfr.ovlEnd, (rs->gkfr.hasOVLclr)+'0',
          rs->gkfr.cnsSta, rs->gkfr.cnsEnd, (rs->gkfr.hasCNSclr)+'0',
          rs->gkfr.cgwSta, rs->gkfr.cgwEnd, (rs->gkfr.hasCGWclr)+'0');

  fprintf(fout, "\tseqFile:%u seqOffset:" F_U64 " srcFile:%u srcOffset:" F_U64 "\n",
          GET_FILEID(rs->gkfr.sequenceOffset), 
          GET_FILEOFFSET(rs->gkfr.sequenceOffset), 
          GET_FILEID(rs->gkfr.sourceOffset), 
          GET_FILEOFFSET(rs->gkfr.sourceOffset)
          );
#endif

  if( strlen(rs->seq) > AS_READ_MAX_LEN )
    fprintf(fout,"LONG FRAGMENT !!!\n");

  fprintf(fout,"\tlength     %d\n", strlen(rs->src));
  fprintf(fout,"\tsource     %s\n", rs->src);

  if(clearRangeOnly){
    uint32  hold_start;
    uint32  hold_end;
    int     length;

    getClearRegion_ReadStruct(rs,&hold_start,&hold_end,READSTRUCT_LATEST);     

    length = hold_end - hold_start;

    fprintf(fout,"\tlength   %d\n", length);
    fprintf(fout,"\tsequence %*s\n", length, rs->seq + rs->gkfr.clrSta);
    fprintf(fout,"\tquality  %*s\n", length, rs->qlt  + rs->gkfr.clrSta);
  }else{
    fprintf(fout,"\tlength   %d\n", strlen(rs->seq));
    fprintf(fout,"\tsequence %s\n", rs->seq);
    fprintf(fout,"\tquality  %s\n", rs->qlt);

#ifdef DUMP_QUALITY_AS_NUMBERS
    fprintf(fout,"\tquality :");
    {
      char *q;
      for (q = rs->qlt; *q; q++)
        fprintf(fout, " %02d", *q - '0');
      fprintf(fout, "\n");
    }
#endif

  }

  return(0);
}


int setClearRegion_ReadStruct(ReadStruct *rs, 
                              uint32 start,
                              uint32 end,
                              uint32 flags){

  uint32 tempFlag = flags;

  // Not valid to set the latest -- which would that be?
  assert (flags!=READSTRUCT_LATEST);
  assert (flags==READSTRUCT_ORIGINAL || flags==READSTRUCT_OVL 
          || flags==READSTRUCT_CGW || flags==READSTRUCT_CNS); 

  // An ordering is encoded here as
  // ORIGINAL >> OVL >> CNS >> CGW.
  // An update to any one propagates up >>
  // such that CGW always has the LATEST.
  // See corresponding get() function.

  if (tempFlag==READSTRUCT_ORIGINAL) {
    rs->gkfr.clrSta = start;
    rs->gkfr.clrEnd = end;
    tempFlag = READSTRUCT_OVL; // fall through
  }
  if (tempFlag==READSTRUCT_OVL) {
    rs->gkfr.ovlSta = start;
    rs->gkfr.ovlEnd = end;
    tempFlag = READSTRUCT_CNS; // fall through
  }
  if (tempFlag==READSTRUCT_CNS) {
    rs->gkfr.cnsSta = start;
    rs->gkfr.cnsEnd = end;
    tempFlag = READSTRUCT_CGW; // fall through
  }
  if (tempFlag==READSTRUCT_CGW) {
    rs->gkfr.cgwSta = start;
    rs->gkfr.cgwEnd = end;
  }

  // Set flags to indicate the last program to modify values
  if (flags==READSTRUCT_OVL) {
    rs->gkfr.hasOVLclr = 1;
    rs->gkfr.hasCNSclr = 0;
    rs->gkfr.hasCGWclr = 0;
  } else if (flags==READSTRUCT_CNS) {
    rs->gkfr.hasCNSclr = 1;
    rs->gkfr.hasCGWclr = 0;
  } else if (flags==READSTRUCT_CGW) {
    rs->gkfr.hasCGWclr = 1;
  }
  return(0);
}





int getAccID_ReadStruct(ReadStruct *rs, CDS_UID_t *accID){
  *accID = rs->gkfr.UID;
  return(0);
}

int getReadIndex_ReadStruct(ReadStruct *rs, CDS_IID_t *readIndex){
  *readIndex = rs->gkfr.readIID;
  return(0);
}

int getClearRegion_ReadStruct(ReadStruct *rs, 
			      uint32 *start, uint32 *end, uint32 flags){

  assert (flags==READSTRUCT_LATEST || flags==READSTRUCT_ORIGINAL 
	  || flags==READSTRUCT_OVL || flags==READSTRUCT_CGW 
	  || flags==READSTRUCT_CNS);

  // An ordering is assumed here as
  // ORIGINAL >> OVL >> CNS >> CGW.
  // Changes to any one were propagated up >>
  // such that LATEST is same as CGW.
  // See corresponding set() function.

  if (flags==READSTRUCT_LATEST || flags==READSTRUCT_CGW) {
    *start = rs->gkfr.cgwSta;
    *end   = rs->gkfr.cgwEnd;
  } else if (flags==READSTRUCT_CNS) {
    *start = rs->gkfr.cnsSta;
    *end   = rs->gkfr.cnsEnd;
  } else if (flags==READSTRUCT_OVL) {
    *start = rs->gkfr.ovlSta;
    *end   = rs->gkfr.ovlEnd;
  } else { // if (flags==READSTRUCT_ORIGINAL)
    *start = rs->gkfr.clrSta;
    *end   = rs->gkfr.clrEnd;
  }
  return(0);
}

int getSource_ReadStruct(ReadStruct *rs, char *src, int length){

  if (strlen(rs->src) + 1 > length) {
    if (src)
      src[0] = 0;
    return(strlen(rs->src) + 1);
  }
  strcpy(src, rs->src);
  return(0);
}

int getSequence_ReadStruct(ReadStruct *rs,
                           char *seq,
                           char *qua,
                           int   len){

  if(strlen(rs->seq) + 1> len){
    if (seq)  seq[0] = 0;
    if (qua)  qua[0] = 0;
    return (strlen(rs->seq) + 1);
  }
  strcpy(seq, rs->seq);
  strcpy(qua, rs->qlt);
  return(0);
}

int getIsDeleted_ReadStruct(ReadStruct *rs, uint32 *isDeleted){
  *isDeleted = rs->gkfr.deleted;
  return(0);
}
