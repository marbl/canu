
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
static char CM_ID[] = "$Id: AS_PER_ReadStruct.c,v 1.3 2005-03-22 19:07:14 jason_miller Exp $";
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
#include "AS_PER_genericStore.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStore_private.h"

int  FragStore_Version = FRAGSTORE_VERSION;


/*****************************************************************************/

ReadStructp new_ReadStruct(void){
  FragRecord *newFR = (FragRecord *)calloc( 1, sizeof(FragRecord));
  AssertPtr(newFR);
  clear_ReadStruct(newFR);
  return ((ReadStructp) newFR);
}
/***************************************************************************/
void        delete_ReadStruct(ReadStructp r){
  FragRecord *FR = (FragRecord *)r;
  free(FR);
}

/***************************************************************************/
void        clear_ReadStruct(ReadStructp r){
#if 0
 Branch_Info_t  branch_info;
#endif

  FragRecord *FR = (FragRecord *)r;
 setAccID_ReadStruct(FR, 0);
 setReadIndex_ReadStruct(FR, 0);
 setReadType_ReadStruct(FR, (FragType)0);
 setSequence_ReadStruct(FR, NULL, NULL);
 setSource_ReadStruct(FR, "");
 setEntryTime_ReadStruct(FR,0);
 setLocID_ReadStruct(FR, 0);
 setLocalePos_ReadStruct(FR, 0,0);
 setScreenMatches_ReadStruct(FR, 0, NULL);
 // setHasModifiedClearRange(FR,0);  // not used
 FR->frag.deleted=0;
 FR->frag.hasQuality=0;
#if FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES 
 FR->frag.hasOVLClearRegion=0;
 FR->frag.hasCNSClearRegion=0;
 FR->frag.hasCGWClearRegion=0;
#endif
 // Reset ORIGINAL clear range.
 // Assume this resets the other clear ranges too.
 setClearRegion_ReadStruct(FR,0,0,READSTRUCT_ORIGINAL); 
}


/***************************************************************************/
int dumpShortFragRecord(ShortFragRecord *sfr, FILE *fout){
   fprintf(fout, "\tDeleted: %d ReadType:%c hasQuality:%d spare1:%d\n",
	   sfr->deleted, sfr->readType, sfr->hasQuality, sfr->spare1);

#if FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES 
   fprintf(fout, 
      "\taccID:" F_UID "  readIdx:" F_IID "  ClearRanges(start,stop,modified):\n", 
      sfr->accID, sfr->readIndex);
   fprintf(fout, "\tOrig(" F_VLS "," F_VLS ") Ovl(" F_VLS "," F_VLS ",%c) Cns(" F_VLS "," F_VLS ",%c) Cgw(" F_VLS "," F_VLS ",%c)\n",
      sfr->clearRegionStart, sfr->clearRegionEnd,
      sfr->ovlRegionStart, sfr->ovlRegionEnd, (sfr->hasOVLClearRegion)+'0',
      sfr->cnsRegionStart, sfr->cnsRegionEnd, (sfr->hasCNSClearRegion)+'0',
      sfr->cgwRegionStart, sfr->cgwRegionEnd, (sfr->hasCGWClearRegion)+'0');
#else
    fprintf(fout, "\taccID:" F_UID "  readIdx:" F_IID " clear(" F_VLS "," F_VLS ")\n",
      sfr->accID, sfr->readIndex, 
      sfr->clearRegionStart, sfr->clearRegionEnd);
#endif

   fprintf(fout, "\tseqFile:%u seqOffset:" F_U64 " srcFile:%u srcOffset:" F_U64 " entryTime:%s\n",
	   GET_FILEID(sfr->sequenceOffset), 
	   GET_FILEOFFSET(sfr->sequenceOffset), 
	   GET_FILEID(sfr->sourceOffset), 
	   GET_FILEOFFSET(sfr->sourceOffset), 
	   ctime(&(sfr->entryTime)));
   return(0);
}

/***************************************************************************/
int dump_ReadStruct(ReadStructp rs, FILE *fout, int clearRangeOnly){
  FragRecord *fr = (FragRecord *)rs;
  int type = fr->frag.readType;
  fprintf(fout,"Dumping FragRecord at 0x%p\n", fr);
  dumpShortFragRecord(&(fr->frag),fout);
  if(type != AS_READ &&
     type != AS_B_READ &&
     type != AS_EXTR &&
     type != AS_TRNR){
    fprintf(fout,"\tlocale : " F_UID " type %c %d\n",
            fr->localeID, (char)type, type);
    if(type == AS_UBAC ||
       type == AS_FBAC){
    fprintf(fout,"\tlocale_pos : [" F_VLS "," F_VLS "]\n", fr->localePosStart, fr->localePosEnd);

    }
  }
  if( strlen(fr->sequence) > AS_READ_MAX_LEN )
	fprintf(fout,"LONG FRAGMENT !!!\n");
  fprintf(fout,"\tsource (solength=" F_SIZE_T ")  : %s\n",
          strlen(fr->source),fr->source);
  if(clearRangeOnly){

    // As of Oct 2001, take advantage of the LATEST clear range modification.
    uint32 hold_start, hold_end;
    int length;
    getClearRegion_ReadStruct(rs,&hold_start,&hold_end,READSTRUCT_LATEST);     
    length = hold_end - hold_start;
    //int length = fr->frag.clearRegionEnd - fr->frag.clearRegionStart;

    fprintf(fout,"\tsequence (selength=%d) : %*s\n",
            length, length, fr->sequence + fr->frag.clearRegionStart);
    fprintf(fout,"\tquality                : %*s\n",
            length, fr->quality + fr->frag.clearRegionStart);
  }else{
    fprintf(fout,"\tsequence (selength=" F_SIZE_T ") : %s\n",
            strlen(fr->sequence), fr->sequence);
    fprintf(fout,"\tquality : %s\n", fr->quality);
  }
  if(fr->frag.numScreenMatches > 0){
    IntScreenMatch *m = fr->matches;
    int i;
    fprintf(fout,"\tScreen Matches (%d)\n",
	    fr->frag.numScreenMatches);
    for(i = 0, m = fr->matches; m != NULL; m = m->next, i++){
      fprintf(fout,"\t %d (" F_COORD "," F_COORD ") w:" F_IID " r:%d (" F_COORD "," F_COORD ") %c\n",
	      i, 
	      m->where.bgn,
	      m->where.end,
	      m->iwhat,
	      m->relevance,
	      m->portion_of.bgn,
	      m->portion_of.end,
	      m->direction);
    }
  }
   return(0);

}


/***************************************************************************/
/* Accessors */
int setAccID_ReadStruct(ReadStructp rs, CDS_UID_t accID){
  FragRecord *FR = (FragRecord *)rs;

  FR->frag.accID = accID;
   return(0);
}
/***************************************************************************/
int setReadIndex_ReadStruct(ReadStructp rs, CDS_IID_t readIndex){
  FragRecord *FR = (FragRecord *)rs;

  FR->frag.readIndex = readIndex;
   return(0);
}

/***************************************************************************/
int setReadType_ReadStruct(ReadStructp rs, FragType r){
  //  int rv=0;
  FragRecord *FR = (FragRecord *)rs;

  FR->frag.readType = r;

   return(0);
}


/***************************************************************************/
int setEntryTime_ReadStruct(ReadStructp rs, time_t entryTime){
  FragRecord *FR = (FragRecord *)rs;

  FR->frag.entryTime = entryTime;
   return(0);
}



/***************************************************************************/
int setClearRegion_ReadStruct(ReadStructp rs, 
                              uint32 start, uint32 end, uint32 flags){
  FragRecord *FR = (FragRecord *)rs;
  uint32 tempFlag = flags;

  // Not valid to set the latest -- which would that be?
  assert (flags!=READSTRUCT_LATEST);
  assert (flags==READSTRUCT_ORIGINAL || flags==READSTRUCT_OVL 
  || flags==READSTRUCT_CGW || flags==READSTRUCT_CNS); 

#if FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES 
  // An ordering is encoded here as
  // ORIGINAL >> OVL >> CNS >> CGW.
  // An update to any one propagates up >>
  // such that CGW always has the LATEST.
  // See corresponding get() function.

  if (tempFlag==READSTRUCT_ORIGINAL) {
    FR->frag.clearRegionStart = start;
    FR->frag.clearRegionEnd = end;
    tempFlag = READSTRUCT_OVL; // fall through
  }
  if (tempFlag==READSTRUCT_OVL) {
    FR->frag.ovlRegionStart = start;
    FR->frag.ovlRegionEnd = end;
    tempFlag = READSTRUCT_CNS; // fall through
  }
  if (tempFlag==READSTRUCT_CNS) {
    FR->frag.cnsRegionStart = start;
    FR->frag.cnsRegionEnd = end;
    tempFlag = READSTRUCT_CGW; // fall through
  }
  if (tempFlag==READSTRUCT_CGW) {
    FR->frag.cgwRegionStart = start;
    FR->frag.cgwRegionEnd = end;
  }

  // Set flags to indicate the last program to modify values
    if (flags==READSTRUCT_OVL) {
      FR->frag.hasOVLClearRegion = 1;
      FR->frag.hasCNSClearRegion = 0;
      FR->frag.hasCGWClearRegion = 0;
    } else if (flags==READSTRUCT_CNS) {
      FR->frag.hasCNSClearRegion = 1;
      FR->frag.hasCGWClearRegion = 0;
    } else if (flags==READSTRUCT_CGW) {
      FR->frag.hasCGWClearRegion = 1;
    }
#else
    // In the old days, there was only one clear region.
    FR->frag.clearRegionStart = start;
    FR->frag.clearRegionEnd = end;
#endif
   return(0);
}


/***************************************************************************/
int setSource_ReadStruct(ReadStructp rs, const char *src){
  FragRecord *FR = (FragRecord *)rs;

  FR->flags |= FRAG_S_SOURCE;
  strcpy(FR->source,src);
   return(0);
}

/***************************************************************************/
int setScreenMatches_ReadStruct(ReadStructp rs, int numScreenMatches, IntScreenMatch *matches){
 
  int i;
  FragRecord *FR = (FragRecord *)rs;
  IntScreenMatch *in_matches=matches;
  
  assert(numScreenMatches < MAX_SCREEN_MATCH);
  FR->frag.numScreenMatches = numScreenMatches;

  // Chain down the list, don't fix up the next pointer value.  We do this on a get
  for (i=0;i<numScreenMatches;i++) {
    AssertPtr(in_matches);
    FR->matches[i] = *in_matches;
    in_matches = in_matches->next;
  }

   return(0);
}


/***************************************************************************/
int setSequence_ReadStruct(ReadStructp rs, char *sequence, char *quality){
  FragRecord *FR = (FragRecord *)rs;
  
  FR->flags |= FRAG_S_SEQUENCE;

  if( !quality || strlen(quality) == 0)
    FR->frag.hasQuality = 0;
  else
    FR->frag.hasQuality = 1;

  strcpy(FR->quality,(quality?quality:""));
  strcpy(FR->sequence,(sequence?sequence:""));

  return(0);

}
/***************************************************************************/
/********** Get ***********/

int getAccID_ReadStruct(ReadStructp rs, CDS_UID_t *accID){
  FragRecord *FR = (FragRecord *)rs;

  *accID =   FR->frag.accID;
   return(0);

}
/**************************************************************************/
int getReadIndex_ReadStruct(ReadStructp rs, CDS_IID_t *readIndex){
  FragRecord *FR = (FragRecord *)rs;
  *readIndex = FR->frag.readIndex;
  return(0);
}

/**************************************************************************/
int getReadType_ReadStruct(ReadStructp rs, FragType *r){
  FragRecord *FR = (FragRecord *)rs;

  *r = (FragType) FR->frag.readType;
  return(0);

}


/***************************************************************************/
int getEntryTime_ReadStruct(ReadStructp rs, time_t *entryTime){
  FragRecord *FR = (FragRecord *)rs;

  *entryTime = FR->frag.entryTime;
  return(0);

}




/****************************************************************************/
int getClearRegion_ReadStruct(ReadStructp rs, 
			      uint32 *start, uint32 *end, uint32 flags){
  FragRecord *FR = (FragRecord *)rs;

  assert (flags==READSTRUCT_LATEST || flags==READSTRUCT_ORIGINAL 
	  || flags==READSTRUCT_OVL || flags==READSTRUCT_CGW 
	  || flags==READSTRUCT_CNS);

#if FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES 
  // An ordering is assumed here as
  // ORIGINAL >> OVL >> CNS >> CGW.
  // Changes to any one were propagated up >>
  // such that LATEST is same as CGW.
  // See corresponding set() function.

  if (flags==READSTRUCT_LATEST || flags==READSTRUCT_CGW) {
    *start = FR->frag.cgwRegionStart;
    *end =   FR->frag.cgwRegionEnd;
  } else if (flags==READSTRUCT_CNS) {
    *start = FR->frag.cnsRegionStart;
    *end =   FR->frag.cnsRegionEnd;
  } else if (flags==READSTRUCT_OVL) {
    *start = FR->frag.ovlRegionStart;
    *end =   FR->frag.ovlRegionEnd;
  } else { // if (flags==READSTRUCT_ORIGINAL)
    *start = FR->frag.clearRegionStart;
    *end =   FR->frag.clearRegionEnd;
  }
#else
  *start = FR->frag.clearRegionStart;
  *end =   FR->frag.clearRegionEnd;
#endif
  return(0);
}

/***************************************************************************/
int getSource_ReadStruct(ReadStructp rs, char *src, int length){
  FragRecord *FR = (FragRecord *)rs;

  if(strlen(FR->source) + 1> length){
    /* strcpy(src,"");   <- this had undesirable behavior when routine was
                            used with a src=NULL length=0 call to determine required
                            buffer space for subsequent call                 */
    if (src) strcpy(src,"");
    return (strlen(FR->source) + 1);
  }
  if((FR->flags & FRAG_S_SOURCE) == 0){
    *src = '\0';
    return 0;
  }
  strcpy(src,FR->source);
  return(0);

}

/***************************************************************************/

int getNumScreenMatches_ReadStruct(ReadStructp rs){
  FragRecord *FR = (FragRecord *)rs;
  
  return(FR->frag.numScreenMatches);
}

/***************************************************************************/

int getScreenMatches_ReadStruct
(ReadStructp rs, IntScreenMatch *matches, int length)
{
  int i;
  FragRecord *fr = (FragRecord *)rs;
  
  if(fr->frag.numScreenMatches == 0)
    return 0;

  if(length < fr->frag.numScreenMatches){
    return fr->frag.numScreenMatches;
  }

  // Copy the data
  memcpy((char *)matches, (char *)fr->matches,  fr->frag.numScreenMatches * sizeof(IntScreenMatch));

  // Adjust the next pointers to be meaningful
  for (i=0;i<fr->frag.numScreenMatches - 1;i++) {
    matches[i].next = matches + i + 1;
  }
  // Last one is NULL
  matches[fr->frag.numScreenMatches - 1].next = NULL;

  return(0);

}


/***************************************************************************/
int getSequence_ReadStruct(ReadStructp rs, char *sequence,
                           char *quality, int length){
  FragRecord *FR = (FragRecord *)rs;

  if(strlen(FR->sequence) + 1> length){
    if (sequence) strcpy(sequence,"");
    if (quality) strcpy(quality,"");
    return (strlen(FR->sequence) + 1);
  }
  if((FR->flags & FRAG_S_SEQUENCE) == 0){
    *sequence = '\0';
    *quality = '\0';
    return 0;
  }

  strcpy(quality, FR->quality);
  strcpy(sequence,FR->sequence);
  return(0);
}

/***************************************************************************/
int getSourceOffset_ReadStruct(ReadStructp rs, int64 *offset){
  FragRecord *FR = (FragRecord *)rs;
  *offset = FR->frag.sourceOffset;
  return(0);
}

/***************************************************************************/
int getSequenceOffset_ReadStruct(ReadStructp rs, int64 *offset){
  FragRecord *FR = (FragRecord *)rs;
  *offset = FR->frag.sequenceOffset;
  return(0);
}
/***************************************************************************/
int setSourceOffset_ReadStruct(ReadStructp rs, int64 offset){
  FragRecord *FR = (FragRecord *)rs;
  FR->frag.sourceOffset = offset;
  return(0);
}

/***************************************************************************/
int setSequenceOffset_ReadStruct(ReadStructp rs, int64 offset){
  FragRecord *FR = (FragRecord *)rs;
  FR->frag.sequenceOffset = offset;
  return(0);
}

/***************************************************************************/
int setLocID_ReadStruct(ReadStructp rs, CDS_UID_t locID){
  FragRecord *FR = (FragRecord *)rs;
  
  FR->localeID = locID;
  return(0);
}
/***************************************************************************/
  int setLocalePos_ReadStruct(ReadStructp rs, uint32 start, uint32 end){
  FragRecord *FR = (FragRecord *)rs;
  FR->localePosStart = start;
  FR->localePosEnd = end;
  return(0);
  }
/***************************************************************************/

int getLocID_ReadStruct(ReadStructp rs, CDS_UID_t *locID){
  FragRecord *FR = (FragRecord *)rs;

  *locID = FR->localeID;
  //  fprintf(stderr,"* Get locID " F_UID "\n", *locID);
  return(0);
  }
/***************************************************************************/
int getLocalePos_ReadStruct(ReadStructp rs, uint32 *start, uint32 *end){
  FragRecord *FR = (FragRecord *)rs;
  *start = FR->localePosStart;
  *end   = FR->localePosEnd;
  return(0);
}

/***************************************************************************/
int getIsDeleted_ReadStruct(ReadStructp rs, uint32 *isDeleted){
  FragRecord *FR = (FragRecord *)rs;
  *isDeleted = FR->frag.deleted;
  return(0);
}
