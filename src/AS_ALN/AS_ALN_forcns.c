
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"

/* 
   Wrappers for finding fragment overlaps with moderate sized indels
   ("bubbles" in CGB that break up unitigging in the presence of
   moderate polymorphisms).

   The wrappers are around the routines Local_Overlap_AS  ([originally]
   defined in AS_ALN_loverlapper.c) and AS_ALN_affine_overlap ([originally]
   defined in AS_ALN_qvaligner.c).  See those functions for more details.

*/


#define AFFINE_QUALITY   /* overlap diff and length reported in affine terms */


/* The following functions make copies of the a and b sequences, returning a
   pointer to the location of the first character of the sequence; Gene's
   local alignment code (like DP_Compare?) seems to assume that the character
   before the beginning of the sequence will be '\0' -- so we copy the
   sequence into a character array starting at the second position, and
   return the location of the second position. 

   Need two versions of the function because of the use of static memory. */

static char *safe_copy_Astring_with_preceding_null(char *in){
  static char* out=NULL;
  static int outsize=0;
  int length=strlen(in);
  if(outsize<length+2){
    if(outsize==0){
      outsize=length * 2 + 2;
      out=(char*)safe_malloc(outsize*sizeof(char));
    } else {
      outsize=length * 2 + 2;
      out=(char*)safe_realloc(out,outsize*sizeof(char));
    }
  }
  out[0]='\0';
  strcpy(out+1,in);
  return(out+1);
}

static char *safe_copy_Bstring_with_preceding_null(char *in){
  static char* out=NULL;
  static int outsize=0;
  int length=strlen(in);
  if(outsize<length+2){
    if(outsize==0){
      outsize=length * 2 + 2;
      out=(char*)safe_malloc(outsize*sizeof(char));
    } else {
      outsize=length * 2 + 2;
      out=(char*)safe_realloc(out,outsize*sizeof(char));
    }
  }
  out[0]='\0';
  strcpy(out+1,in);
  return(out+1);
}

  
Overlap *Local_Overlap_AS_forCNS(char *a, char *b,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
				 CompareOptions what){

  InternalFragMesg  A = {0}, B = {0};
  OverlapMesg  *O;
  int alen,blen,del,sub,ins,affdel,affins,blockdel,blockins;
  double errRate,errRateAffine;
  int AFFINEBLOCKSIZE=4;
  static Overlap o;
  int where=0;

  assert((0.0 <= erate) && (erate <= AS_MAX_ERROR_RATE));

  A.sequence = safe_copy_Astring_with_preceding_null(a);
  A.quality = NULL;
  A.iaccession = 1;
  A.eaccession = AS_UID_fromInteger(A.iaccession);

  B.sequence = safe_copy_Bstring_with_preceding_null(b);
  B.quality = NULL;
  B.iaccession = 2;
  B.eaccession = AS_UID_fromInteger(B.iaccession);

  O = Local_Overlap_AS(&A,&B,beg,end,
		       opposite,
		       erate,
		       thresh,
		       minlen,
		       AS_FIND_LOCAL_ALIGN,
		       &where);
  if(O==NULL){
    return(NULL);
  }
 
  Analyze_Affine_Overlap_AS(&A,&B,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
		    &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, NULL);
  
  errRate = (sub+ins+del)/(double)(alen+ins);
  
  errRateAffine = (sub+affins+affdel)/ (double)(alen-del+affins+affdel);

  o.begpos=O->ahg;
  o.endpos=O->bhg;
#ifdef AFFINE_QUALITY
#ifdef REVISED_QUALITY
  o.length=O->length;
  o.diffs=O->diffs;
#else
  o.length=(alen+blen+affins+affdel)/2;
  o.diffs=sub+affins+affdel;
#endif
#else
  o.length=(alen+blen+ins+del)/2;
  o.diffs=sub+ins+del;
#endif
  o.comp=opposite;
  o.trace=Unpack_Alignment_AS(O);
  if(O->aifrag==2){/*The OverlapMesg gives b first for nonnegative ahang*/
    int i=0;
    while(o.trace[i]!=0){
      o.trace[i++]*=-1;
    }
    o.begpos*=-1;
    o.endpos*=-1;
  }


  /* At this point, we may have gaps at the ends of the sequences (either
     before the first base or after the last).

     These may occur, e.g., due to disagreements in optimal non-affine 
     alignment (as determined by Boundary) and optimal affine alignment
     (as determined by OKNAffine).

     These gaps give consensus hiccups, so we want to eliminate them, and
     change the hangs appropriately.

  */

  { int i=0;
    int j=0;
    int changeahang=0;
    int changebhang=0;
    int fullLenA = strlen(a);
    int fullLenB = strlen(b);
    char c;
    //printf("Trace (lens %d %d):",fullLenA,fullLenB);
    while(o.trace[i]!=0){
      c='*';
      if(o.trace[i]<-fullLenA){
	changebhang++;
      } else if (o.trace[i]>fullLenB){
	changebhang--;
      } else if (o.trace[i]==-1){
	changeahang--;
      } else if (o.trace[i]==1){
	changeahang++;
      } else {
	c=' ';
	o.trace[j++]=o.trace[i];
      }
      //printf(" %c%d",c,o.trace[i]);
      i++;
    }
    //printf("\n");
    o.trace[j]=0;
    o.begpos+=changeahang;
    o.endpos+=changebhang;
  }


  return(&o);
}



Overlap *Affine_Overlap_AS_forCNS(char *a, char *b,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
				 CompareOptions what){
  InternalFragMesg  A, B;
  OverlapMesg  *O;
  int alen,blen,del,sub,ins,affdel,affins,blockdel,blockins;
  double errRate,errRateAffine;
  extern int AS_ALN_TEST_NUM_INDELS;
  int orig_TEST_NUM_INDELS;
  int AFFINEBLOCKSIZE=4;
  int where=0;
  static Overlap o;

  assert((0.0 <= erate) && (erate <= AS_MAX_ERROR_RATE));

  orig_TEST_NUM_INDELS = AS_ALN_TEST_NUM_INDELS;
  AS_ALN_TEST_NUM_INDELS = 0;

  A.sequence = safe_copy_Astring_with_preceding_null(a);
  A.quality = NULL;
  A.iaccession = 1;
  A.eaccession = AS_UID_fromInteger(A.iaccession);

  B.sequence = safe_copy_Bstring_with_preceding_null(b);
  B.quality = NULL;
  B.iaccession = 2;
  B.eaccession = AS_UID_fromInteger(B.iaccession);

  O =AS_ALN_affine_overlap(&A,&B,beg,end,
		       opposite,
		       erate,
		       thresh,
		       minlen,
		       AS_FIND_AFFINE_ALIGN,
		       &where);
  if(O==NULL){
    AS_ALN_TEST_NUM_INDELS = orig_TEST_NUM_INDELS;
    return(NULL);
  }
 
  Analyze_Affine_Overlap_AS(&A,&B,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
		    &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, NULL);
  
  errRate = (sub+ins+del)/(double)(alen+ins);
  
  errRateAffine = (sub+affins+affdel)/ (double)(alen-del+affins+affdel);

  o.begpos=O->ahg;
  o.endpos=O->bhg;
#ifdef AFFINE_QUALITY
  o.length=(alen+blen+affins+affdel)/2;
  o.diffs=sub+affins+affdel;
#else
  o.length=(alen+blen+ins+del)/2;
  o.diffs=sub+ins+del;
#endif
  o.comp=opposite;
  o.trace=Unpack_Alignment_AS(O);
  if(O->aifrag==2){/*The OverlapMesg gives b first for nonnegative ahang*/
    int i=0;
    while(o.trace[i]!=0){
      o.trace[i++]*=-1;
    }
    o.begpos*=-1;
    o.endpos*=-1;
  }

  /* At this point, we may have gaps at the ends of the sequences (either
     before the first base or after the last).

     These seem to occur due to disagreements in optimal non-affine alignment
     (as determined by Boundary) and optimal affine alignment
     (as determined by OKNAffine).

     These gaps give consensus hiccups, so we want to eliminate them, and
     change the hangs appropriately.

  */

  { int i=0;
    int j=0;
    int changeahang=0;
    int changebhang=0;
    int fullLenA = strlen(a);
    int fullLenB = strlen(b);
    while(o.trace[i]!=0){
      if(o.trace[i]<-fullLenA){
	changebhang++;
      } else if (o.trace[i]>fullLenB){
	changebhang--;
      } else if (o.trace[i]==-1){
	changeahang--;
      } else if (o.trace[i]==1){
	changeahang++;
      } else {
	o.trace[j++]=o.trace[i];
      }
      i++;
    }
    o.trace[j]=0;
    o.begpos+=changeahang;
    o.endpos+=changebhang;
  }

  AS_ALN_TEST_NUM_INDELS = orig_TEST_NUM_INDELS;

  return(&o);
}
