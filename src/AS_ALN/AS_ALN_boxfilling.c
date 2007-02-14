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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "AS_global.h"
#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"


#define ENDHACK
#ifdef ENDHACK
  static int ENDGAPHACK=3;
  static int MINCORESEG=30;
#endif


#define MINUS_ONE 0   /* whether "zero" is one before the string or the start */

#define AFFINE_QUALITY   /* overlap diff and length reported in affine terms */

static void safe_suffix(char **Dest,int *DestLen,char *src,int start){
  int len=strlen(src+start);
  if(*DestLen==0||*DestLen<len+1){
    *DestLen=len+1;
    if(*Dest==NULL){
      *Dest=(char*)safe_malloc(sizeof(char)*(*DestLen));
    } else {
      *Dest=(char*)safe_realloc(*Dest,sizeof(char)*(*DestLen));
    }
  }
  strcpy(*Dest,src+start);
  assert(strlen(*Dest)==len);
}


/* Given fragments a and b, find the best overlap between them using
   local alignment code.  The purpose is to permit "bubbles" in the fragment
   graph due to multi-base polymorphisms to be smoothed out, such that a 
   unitig can be constructed.

   The method relies on Myers' local overlap code.  MORE DETAILS????

   The local overlap code returns the best overlap available as a chain of
   matching segments.

   The current function is concerned with three things:
   
   1) Providing a wrapper around the local alignment code so that users
   can carry out operations on the level of the "message" (IFM, OVL, etc).

   2) Once a set of local alignments has been found, determine the best chaining
   using the local overlap code

   3) Resolve any overlaps within the best chain

   The following values for 'what' may or may not be relevant here:
   AS_FIND_LOCAL_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_LOCAL_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.
   AS_FIND_LOCAL_ALIGN_NO_TRACE
      A hack by SAK...don't compute delta encoding.

   As with DP_Compare_AS, the resulting overlap message occupies
   memory which is persistently owned (reused) by the function; it
   must be copied if it is to be retained beyond the given call.

*/


OverlapMesg *BoxFill_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where)
{ int   alen,  blen;
  int   atop,  btop;
  char *aseq, *bseq;
  int   pos1,  pos2;
  int   dif1,  dif2;
  int NumSegs=0;
  int seglen,sumseglen=0;
  int noCoreSeg=1;
  int i;
  int forcenext=0;
  int *trace;

  static char *Ausable=NULL, *Busable=NULL;
  static int AuseLen=0, BuseLen=0;

#ifdef ENDHACK
  int coreseglen=MIN(MINCORESEG,minlen);
#else
  int coreseglen=minlen;
#endif

  double avgerror=0.;

  static OverlapMesg QVBuffer;

  Local_Segment *local_results=NULL;
  Local_Overlap *O=NULL;
  
  assert(what == AS_FIND_LOCAL_ALIGN || what == AS_FIND_LOCAL_ALIGN_NO_TRACE
	 || what == AS_FIND_LOCAL_OVERLAP );

  aseq = a->sequence;  /* Setup sequence access */
  bseq = b->sequence;
  alen = strlen(aseq);
  blen = strlen(bseq);

  aseq -= MINUS_ONE;
  bseq -= MINUS_ONE;

  if (opposite){                 /* Compare in opposite orientation. */
    Complement_Fragment_AS(b);
  }


  /* now generate the substrings of the input sequences which contain the
     portion usable according to the user's begin/end parameters */
#define EXPAND_BAND 10 /* amount to pad ends of  user's begin..end range */
#define STRETCH 1.5  /* amount of "stretching" allowed in bounding 
			overlap end */
#define BIGPAD 500 /* amount of extra slop allowed in bounding overlap end */


  // not really using the "suffix" part of next calls; just a safe string copy
  safe_suffix(&Ausable,&AuseLen,aseq+MINUS_ONE,0);
  safe_suffix(&Busable,&BuseLen,bseq+MINUS_ONE,0);

  for(i=0;i<MAX(0,beg-EXPAND_BAND);i++){
    Ausable[i]='N';
  }
  assert(alen>=i);
  for(i=0;i<MAX(0,-end-EXPAND_BAND);i++){
    Busable[i]='N';
  }
  assert(blen>=i);

  for(i=MAX(0,end)+(blen+MIN(0,end))*STRETCH+BIGPAD;i<alen;i++){
    Ausable[i]='N';
  }

  i=-MIN(0,beg)+(alen-MAX(0,beg))*STRETCH+BIGPAD;
  assert(i>=0);

  for(;i<blen;i++){
    Busable[i]='N';
  }


  /* Notes on handling of reverse complement overlap searching:
     1) The implementation here uses Find_Local_Segments and 
        Find_Local_Overlaps only in the forward direction; the interaction
	of these two routines on reverse orientation searches is non-obvious
	[even Gene agrees], so we avoid it.
     2) Thus, it is the responsibility of this routine to reverse the B
        sequence prior to running these routines.
  */



  local_results=Find_Local_Segments(Ausable-MINUS_ONE,strlen(Ausable),
				    Busable-MINUS_ONE,strlen(Busable),
				    LOCAL_FORW, 16, erate, &NumSegs);

  if(NumSegs==0){
    goto nooverlap;    
  }

  O=Find_Local_Overlap(alen,blen,
		       0 /*comp==0 -> fwd orientation */,
		       0 /*nextbest==0 -> best overlap*/,
		       local_results,NumSegs,
		       minlen
#ifdef ENDHACK
		       -2*ENDGAPHACK
#endif
		        ,1.);

  if(O==NULL){
    goto nooverlap;
  }

  // thoughts for clark: test for O!=NULL?  minimum segment size?  kmer size?

#if MINUS_ONE == 0 
  // coordinates from Find_Local routines will be one off from
  // those expected by the trace routines, so adjust them!
  //	O->begpos+= (O->begpos>=0 ? 1 : -1);
  //	O->endpos+= (O->endpos>=0 ? 1 : -1);
  for(i=0;i<=O->num_pieces;i++){
    if(i<O->num_pieces){
      O->chain[i].piece.abpos++;
      O->chain[i].piece.bbpos++;
      O->chain[i].piece.aepos++;
      O->chain[i].piece.bepos++;
    }
  }
#endif



  //AS_Local_Trace assumes string pointer one before start of string!
  trace = AS_Local_Trace(O,Ausable-1,Busable-1);


  QVBuffer.quality=avgerror;

  /* Now that you have finish point, first compute complete alignment
     by calling AS_ALN_OKNAlign or extrapolate start point as dictated by
     `what'.  Then build an overlap record modelling the overlap.   */

  { int ahang, bhang;

    //printf("Begpos for OverlapMesg == %d\n",O->begpos);

    ahang=O->begpos; /* N.B.: can be changed in AS_Local_Trace*/
    //printf("ahang for OverlapMesg set to %d\n",ahang);

    bhang=O->endpos;

    { 

      int i=0;
      int j=0;
      int changeahang=0;
      int changebhang=0;
      int fullLenA = strlen(a->sequence);
      int fullLenB = strlen(b->sequence);
      while(trace[i]!=0){

	if(trace[i]<-fullLenA){
	  changebhang++;
	} else if (trace[i]>fullLenB){
	  changebhang--;
	} else if (trace[i]==-1){
	  changeahang--;
	} else if (trace[i]==1){
	  changeahang++;
	} else {
	  trace[j++]=trace[i];
	}
	i++;
      }
      trace[j]=0;
      ahang+=changeahang;
      bhang+=changebhang;
    }


    if (ahang < 0 || ahang == 0 && bhang > 0)
      { if (bhang >= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = -ahang;
        QVBuffer.bhg = -bhang;
        if (opposite)
          QVBuffer.orientation = AS_OUTTIE;
        else
          QVBuffer.orientation = AS_NORMAL;
        QVBuffer.aifrag = b->iaccession;
        QVBuffer.bifrag = a->iaccession;
        if (trace != NULL)
          { int i;
            for (i = 0; trace[i] != 0; i++)
              trace[i] = -trace[i];
          }
      }
    else
      { if (bhang <= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = ahang;
        QVBuffer.bhg = bhang;
        if (opposite)
          QVBuffer.orientation = AS_INNIE;
        else
          QVBuffer.orientation = AS_NORMAL;
        QVBuffer.aifrag = a->iaccession;
        QVBuffer.bifrag = b->iaccession;
      }
    assert(trace!=NULL);
    QVBuffer.delta = Pack_Alignment_AS(trace,QVBuffer.ahg);

    *where = ahang;

  }

  if (opposite){
    Complement_Fragment_AS(b);
  }
  safe_free(O);

  return (&QVBuffer);

nooverlap:
  if (opposite)
    Complement_Fragment_AS(b);
  safe_free(O);
  return ((OverlapMesg*)NULL);
  
}
